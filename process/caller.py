from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
from tabulate import tabulate

import process.constraints as constraints
from process import data_structure
from process.final import finalise
from process.io.mfile import MFile
from process.iteration_variables import set_scaled_iteration_variable
from process.objectives import objective_function
from process.process_output import OutputFileManager

if TYPE_CHECKING:
    from process.main import Models

logger = logging.getLogger(__name__)


class Caller:
    """Calls physics and engineering models."""

    def __init__(self, models: Models) -> None:
        """Initialise all physics and engineering models.

        To ensure that, at the start of a run, all physics/engineering
        variables are fully initialised with consistent values, the models are
        called with the initial optimisation paramters, x.

        :param models: physics and engineering model objects
        :type models: Models
        """
        self.models = models

    @staticmethod
    def check_agreement(
        previous: float | np.ndarray, current: float | np.ndarray
    ) -> bool:
        """Compare previous and current arrays for agreement within a tolerance.

        :param previous: value(s) from previous models evaluation
        :type previous: float | np.ndarray
        :param current: value(s) from current models evaluation
        :type current: float | np.ndarray
        :return: whether values agree or not
        :rtype: bool
        """
        # Check for same shape: mfile length can change between iterations
        if isinstance(previous, float) or previous.shape == current.shape:
            return np.allclose(previous, current, rtol=1.0e-6, equal_nan=True)
        return False

    def call_models(self, xc: np.ndarray, m: int) -> tuple[float, np.ndarray]:
        """Evalutate models until results are idempotent.

        Ensure objective function and constraints are idempotent before returning.

        :param xc: optimisation parameters
        :type xc: np.ndarray
        :param m: number of constraints
        :type m: int
        :raises RuntimeError: if values are non-idempotent after successive
        evaluations
        :return: objective function and constraints
        :rtype: Tuple[float, np.ndarray]
        """
        objf_prev = None
        conf_prev = None

        # Evaluate models up to 10 times; any more implies non-converging values
        for _ in range(10):
            self._call_models_once(xc)
            # Evaluate objective function and constraints
            objf = objective_function(data_structure.numerics.minmax)
            conf, _, _, _, _ = constraints.constraint_eqns(m, -1)

            if objf_prev is None and conf_prev is None:
                # First run: run again to check idempotence
                logger.debug("New optimisation parameter vector being evaluated")
                objf_prev = objf
                conf_prev = conf
                continue

            # Check for idempotence
            if self.check_agreement(objf_prev, objf) and self.check_agreement(
                conf_prev, conf
            ):
                # Idempotent: no longer changing, so return
                logger.debug(
                    "Model evaluations idempotent, returning objective "
                    "function and constraints"
                )
                return objf, conf

            # Not idempotent: still changing, so evaluate models again
            logger.debug("Model evaluations not idempotent: evaluating again")
            objf_prev = objf
            conf_prev = conf

        raise RuntimeError(
            "After 10 model evaluations at the current optimisation parameter "
            "vector, values for the objective function and constraints haven't "
            "converged (don't produce idempotent values)."
        )

    def call_models_and_write_output(self, xc: np.ndarray, ifail: int) -> None:
        """Evaluate models until results are idempotent, then write output files.

        Ensure all outputs in mfile are idempotent before returning, by
        evaluating models multiple times. Typically used at the end of an
        optimisation, or in a non-optimising evaluation. Writes OUT.DAT and
        MFILE.DAT with final results.

        :param xc: optimisation parameter
        :type xc: np.ndarray
        :param ifail: return code of solver
        :type ifail: int
        :raises RuntimeError: if values are non-idempotent after successive
        evaluations
        """
        # TODO The only way to ensure idempotence in all outputs is by comparing
        # mfiles at this stage
        previous_mfile_data = None

        try:
            # Evaluate models up to 10 times; any more implies non-converging values
            for _ in range(10):
                # Divert OUT.DAT and MFILE.DAT output to scratch files for
                # idempotence checking
                OutputFileManager.open_idempotence_files()
                self._call_models_once(xc)
                # Write mfile
                finalise(self.models, ifail)

                # Extract data from intermediate idempotence-checking mfile
                mfile_path = (
                    data_structure.global_variables.output_prefix
                ) + "IDEM_MFILE.DAT"
                mfile = MFile(mfile_path)
                # Create mfile dict of float values: only compare floats
                mfile_data = {
                    var: val
                    for var in mfile.data
                    if isinstance(val := mfile.data[var].get_scan(-1), float)
                }

                if previous_mfile_data is None:
                    # First run: need another run to compare with
                    logger.debug(
                        "New mfile created: evaluating models again to check idempotence"
                    )
                    previous_mfile_data = mfile_data.copy()
                    continue

                # Compare previous and current mfiles for agreement
                nonconverged_vars = {}
                for var in previous_mfile_data:
                    previous_value = previous_mfile_data[var]
                    current_value = mfile_data.get(var, np.nan)
                    if self.check_agreement(previous_value, current_value):
                        continue
                    # Value has changed between previous and current mfiles
                    nonconverged_vars[var] = [
                        previous_value,
                        current_value,
                    ]

                if len(nonconverged_vars) == 0:
                    # Previous and current mfiles agree (idempotent)
                    logger.debug("Mfiles idempotent, returning")
                    # Divert OUT.DAT and MFILE.DAT output back to original files
                    # now idempotence checking complete
                    OutputFileManager.close_idempotence_files()
                    # Write final output file and mfile
                    finalise(self.models, ifail)
                    return

                # Mfiles not yet idempotent: need to re-evaluate models
                logger.debug("Mfiles not idempotent, evaluating models again")
                previous_mfile_data = mfile_data.copy()

            # Values haven't all stabilised after 10 evaluations
            # Which variables are still changing?
            non_idempotent_warning = (
                "Model evaluations at the current optimisation parameter vector "
                "don't produce idempotent values in the final output."
            )
            non_idempotent_table = tabulate(
                [[k, v[0], v[1]] for k, v in nonconverged_vars.items()],
                headers=["Variable", "Previous value", "Current value"],
            )

            warnings.warn(
                f"\033[93m{non_idempotent_warning}\n{non_idempotent_table}\033[0m",
                stacklevel=2,
            )

            # Close idempotence files, write final output file and mfile
            OutputFileManager.close_idempotence_files()
            finalise(
                self.models,
                ifail,
                non_idempotent_msg=non_idempotent_warning + "\n" + non_idempotent_table,
            )
            return

        except Exception:
            # If exception in model evaluations delete intermediate idempotence
            # files to clean up
            OutputFileManager.close_idempotence_files()
            raise

    def _call_models_once(self, xc: np.ndarray) -> None:
        """Call the physics and engineering models.

        This method is the principal caller of all the physics and
        engineering models. Some are Fortran subroutines within modules, others
        will be methods on Python model objects.
        :param xc: Array of optimisation parameters
        :type xc: np.array
        """
        # Number of active iteration variables
        nvars = len(xc)

        # Increment the call counter
        data_structure.numerics.ncalls = data_structure.numerics.ncalls + 1

        # Convert variables
        set_scaled_iteration_variable(xc, nvars)

        # Perform the various function calls
        # Stellarator caller
        if data_structure.stellarator_variables.istell != 0:
            self.models.stellarator.run(output=False)
            # TODO Is this return safe?
            return

        # Inertial Fusion Energy calls
        if data_structure.ife_variables.ife != 0:
            self.models.ife.run(output=False)
            return

        # Tokamak calls
        # Plasma geometry model
        self.models.plasma_geom.plasma_geometry()

        # Machine Build Model
        # Radial build
        self.models.build.run()

        self.models.physics.physics()

        # Toroidal field coil model

        # Toroidal field coil resistive model
        if data_structure.tfcoil_variables.i_tf_sup == 0:
            self.models.copper_tf_coil.run(output=False)

        # Toroidal field coil superconductor model
        if data_structure.tfcoil_variables.i_tf_sup == 1:
            self.models.sctfcoil.run(output=False)

        if data_structure.tfcoil_variables.i_tf_sup == 2:
            self.models.aluminium_tf_coil.run(output=False)

        # Poloidal field and central solenoid model
        self.models.pfcoil.run()

        # Pulsed reactor model
        self.models.pulse.run(output=False)

        # First wall model
        self.models.fw.run()

        # Blanket model
        """Blanket switch values
        No.  |  model
        ---- | ------
        1    |  CCFE HCPB model
        2    |  KIT HCPB model
        3    |  CCFE HCPB model with Tritium Breeding Ratio calculation
        4    |  KIT HCLL model
        5    |  DCLL model
        """
        if data_structure.fwbs_variables.i_blanket_type == 1:
            # CCFE HCPB model
            self.models.ccfe_hcpb.run(output=False)
        # i_blanket_type = 2, KIT HCPB removed
        # i_blanket_type = 3, CCFE HCPB with TBR calculation removed
        # i_blanket_type = 4, KIT HCLL removed
        elif data_structure.fwbs_variables.i_blanket_type == 5:
            # DCLL model
            self.models.dcll.run(output=False)

        self.models.divertor.run(output=False)

        self.models.cryostat.run()

        # Structure Model
        self.models.structure.run(output=False)

        # Tight aspect ratio machine model
        if (
            data_structure.physics_variables.itart == 1
            and data_structure.tfcoil_variables.i_tf_sup != 1
        ):
            self.models.tfcoil.cntrpst()

        # Toroidal field coil power model
        self.models.power.tfpwr(output=False)

        # Poloidal field coil power model
        self.models.power.pfpwr(output=False)

        # Plant heat transport part 1
        self.models.power.component_thermal_powers()

        # Cryoplant loads
        self.models.power.calculate_cryo_loads()

        # Vacuum model
        self.models.vacuum.run(output=False)

        # Buildings model
        self.models.buildings.run(output=False)

        # Plant AC power requirements
        self.models.power.acpow(output=False)

        # Plant heat transport pt 2 & 3
        self.models.power.plant_electric_production()

        # Availability model
        self.models.availability.run(output=False)

        # Water usage in secondary cooling system
        self.models.water_use.run(output=False)

        # Costs model
        """Cost switch values
        No.  |  model
        ---- | ------
        0    |  1990 costs model
        1    |  2015 Kovari model
        2    |  Custom model
        """
        self.models.costs.run()

        # FISPACT and LOCA model (not used)- removed


def write_output_files(models: Models, ifail: int) -> None:
    """Evaluate models and write output files (OUT.DAT and MFILE.DAT).

    :param models: physics and engineering models
    :type models: Models
    :param ifail: solver return code
    :type ifail: int
    """
    n = data_structure.numerics.nvar
    x = data_structure.numerics.xcm[:n]
    # Call models, ensuring output mfiles are fully idempotent
    caller = Caller(models)
    caller.call_models_and_write_output(
        xc=x,
        ifail=ifail,
    )
