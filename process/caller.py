from __future__ import annotations
from process import fortran as ft
import numpy as np
import logging
from process.final import finalise
from process.io.mfile import MFile
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from typing import Union, Tuple, TYPE_CHECKING

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
        previous: Union[float, np.ndarray], current: Union[float, np.ndarray]
    ) -> bool:
        """Compare previous and current arrays for agreement within a tolerance.

        :param previous: value(s) from previous models evaluation
        :type previous: float | np.ndarray
        :param current: value(s) from current models evaluation
        :type current: float | np.ndarray
        :return: whether values agree or not
        :rtype: bool
        """
        return np.allclose(previous, current, rtol=1.0e-6)

    def call_models(self, xc: np.ndarray, m: int) -> Tuple[float, np.ndarray]:
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
            objf = ft.function_evaluator.funfom()
            conf, _, _, _, _ = ft.constraints.constraint_eqns(m, -1)

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
        previous_mfile_arr = None

        try:
            # Evaluate models up to 10 times; any more implies non-converging values
            for _ in range(10):
                # Divert OUT.DAT and MFILE.DAT output to scratch files for
                # idempotence checking
                ft.init_module.open_idempotence_files()
                self._call_models_once(xc)
                # Write mfile
                finalise(self.models, ifail)

                # Extract data from intermediate idempotence-checking mfile
                mfile_path = (
                    f2py_compatible_to_string(ft.global_variables.output_prefix)
                    + "IDEM_MFILE.DAT"
                )
                mfile = MFile(mfile_path)
                mfile_data = {}
                for var in mfile.data.keys():
                    mfile_data[var] = mfile.data[var].get_scan(-1)

                # Extract floats from mfile dict into array for straightforward
                # comparison: only compare floats
                current_mfile_arr = np.array(
                    [val for val in mfile_data.values() if isinstance(val, float)]
                )
                if previous_mfile_arr is None:
                    # First run: need another run to compare with
                    logger.debug(
                        "New mfile created: evaluating models again to check idempotence"
                    )
                    previous_mfile_arr = np.copy(current_mfile_arr)
                    continue

                if self.check_agreement(previous_mfile_arr, current_mfile_arr):
                    # Previous and current mfiles agree (idempotent)
                    logger.debug("Mfiles idempotent, returning")
                    # Divert OUT.DAT and MFILE.DAT output back to original files
                    # now idempotence checking complete
                    ft.init_module.close_idempotence_files()
                    # Write final output file and mfile
                    finalise(self.models, ifail)
                    return

                # Mfiles not yet idempotent: re-evaluate models
                logger.debug("Mfiles not idempotent, evaluating models again")
                previous_mfile_arr = np.copy(current_mfile_arr)

            raise RuntimeError(
                "Model evaluations at the current optimisation parameter vector "
                "don't produce idempotent values in the final output."
            )
        except Exception:
            # If exception in model evaluations or idempotence can't be
            # achieved, delete intermediate idempotence files to clean up
            ft.init_module.close_idempotence_files()
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
        nvars = xc.shape[0]

        # Increment the call counter
        ft.numerics.ncalls = ft.numerics.ncalls + 1

        # Convert variables
        ft.define_iteration_variables.convxc(xc, nvars)

        # Perform the various function calls
        # Stellarator caller
        if ft.stellarator_variables.istell != 0:
            self.models.stellarator.run(output=False)
            # TODO Is this return safe?
            return

        # Inertial Fusion Energy calls
        if ft.ife_variables.ife != 0:
            self.models.ife.run(output=False)
            return

        # Tokamak calls
        # Plasma geometry model
        self.models.plasma_geom.geomty()

        # Machine Build Model
        # Radial build
        if ft.build_variables.tf_in_cs == 1:
            self.models.build.tf_in_cs_bore_calc()
        self.models.build.calculate_radial_build(output=False)

        # Vertical build
        self.models.build.vbuild(output=False)

        self.models.physics.physics()

        # Toroidal field coil model

        self.models.tfcoil.run()

        # Toroidal field coil superconductor model
        if ft.tfcoil_variables.i_tf_sup == 1:
            self.models.sctfcoil.run(output=False)

        # Poloidal field and central solenoid model
        self.models.pfcoil.run()

        # Pulsed reactor model
        self.models.pulse.run(output=False)

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
        if ft.fwbs_variables.iblanket == 1:
            # CCFE HCPB model
            self.models.ccfe_hcpb.run(output=False)
        # iblanket = 2, KIT HCPB removed
        elif ft.fwbs_variables.iblanket == 3:
            # CCFE HCPB model with Tritium Breeding Ratio calculation
            self.models.ccfe_hcpb.run(output=False)
            ft.fwbs_variables.tbr = self.models.ccfe_hcpb.tbr_shimwell(
                ft.fwbs_variables.breeder_f,
                ft.fwbs_variables.li6enrich,
                ft.fwbs_variables.iblanket_thickness,
                output=False,
            )
        # iblanket = 4, KIT HCLL removed
        elif ft.fwbs_variables.iblanket == 5:
            # DCLL model
            self.models.dcll.run(output=False)

        self.models.divertor.run(output=False)

        # Structure Model
        self.models.structure.run(output=False)

        # Tight aspect ratio machine model
        if ft.physics_variables.itart == 1 and ft.tfcoil_variables.i_tf_sup != 1:
            self.models.tfcoil.cntrpst()

        # Toroidal field coil power model
        self.models.power.tfpwr(output=False)

        # Poloidal field coil power model
        self.models.power.pfpwr(output=False)

        # Plant heat transport part 1
        self.models.power.power1()

        # Vacuum model
        self.models.vacuum.run(output=False)

        # Buildings model
        self.models.buildings.run(output=False)

        # Plant AC power requirements
        self.models.power.acpow(output=False)

        # Plant heat transport pt 2 & 3
        self.models.power.power2(output=False)

        self.models.power.power3(output=False)

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
    n = ft.numerics.nvar
    x = ft.numerics.xcm[:n]
    # Call models, ensuring output mfiles are fully idempotent
    caller = Caller(models)
    caller.call_models_and_write_output(
        xc=x,
        ifail=ifail,
    )
