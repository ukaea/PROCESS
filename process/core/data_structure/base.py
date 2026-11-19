"""The PROCESS data structure dataclass."""

from __future__ import annotations

import copy
import datetime as dt
import getpass
import logging
import socket
import subprocess  # noqa: S404
from dataclasses import dataclass, fields
from pathlib import Path
from typing import TYPE_CHECKING, get_type_hints

import numpy as np

import process
from process.core import constants, process_output
from process.core.data_structure.check_input import check_process
from process.core.data_structure.input import (
    INPUT_VARIABLES,
    InputVariable,
    ValidInputTypes,
)
from process.core.exceptions import ProcessValidationError, ProcessValueError
from process.core.io.in_dat.base import InDat
from process.core.solver import iteration_variables
from process.core.solver.constraints import ConstraintManager
from process.data_structure.blanket_variables import BlanketData  # noqa: TC001
from process.data_structure.build_variables import BuildData  # noqa: TC001
from process.data_structure.buildings_variables import BuildingsData  # noqa: TC001
from process.data_structure.ccfe_hcpb_variables import CCFEHCPBData  # noqa: TC001
from process.data_structure.constraint_variables import ConstraintData  # noqa: TC001
from process.data_structure.cost_2015_variables import Cost2015Data  # noqa: TC001
from process.data_structure.cost_variables import CostData  # noqa: TC001
from process.data_structure.cs_fatigue_variables import CSFatigueData  # noqa: TC001
from process.data_structure.current_drive_variables import CurrentDriveData  # noqa: TC001
from process.data_structure.dcll_variables import DCLLData  # noqa: TC001
from process.data_structure.divertor_variables import DivertorData  # noqa: TC001
from process.data_structure.first_wall_variables import FirstWallData  # noqa: TC001
from process.data_structure.fwbs_variables import FWBSData  # noqa: TC001
from process.data_structure.global_variables import GlobalData  # noqa: TC001
from process.data_structure.heat_transport_variables import HeatTransportData  # noqa: TC001
from process.data_structure.ife_variables import IFEData  # noqa: TC001
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData  # noqa: TC001
from process.data_structure.neoclassics_variables import NeoclassicsData  # noqa: TC001
from process.data_structure.numerics import FiguresOfMerit, NumericsData, PROCESSRunMode
from process.data_structure.pf_power_variables import PFPowerData  # noqa: TC001
from process.data_structure.pfcoil_variables import PFCoilData  # noqa: TC001
from process.data_structure.physics_variables import PhysicsData  # noqa: TC001
from process.data_structure.power_variables import PowerData  # noqa: TC001
from process.data_structure.primary_pumping_variables import PrimaryPumpingData  # noqa: TC001
from process.data_structure.pulse_variables import PulseData  # noqa: TC001
from process.data_structure.rebco_variables import RebcoData  # noqa: TC001
from process.data_structure.reinke_variables import ReinkeData  # noqa: TC001
from process.data_structure.scan_variables import ScanData  # noqa: TC001
from process.data_structure.stellarator_configuration import StellaratorConfigData  # noqa: TC001
from process.data_structure.stellarator_variables import (
    StellaratorData,
    StellaratorModel,
)
from process.data_structure.structure_variables import StructureData  # noqa: TC001
from process.data_structure.superconducting_tf_coil_variables import (
    SuperconductingTFData,  # noqa: TC001
)
from process.data_structure.tfcoil_variables import TFData  # noqa: TC001
from process.data_structure.times_variables import TimesData  # noqa: TC001
from process.data_structure.vacuum_variables import VacuumData  # noqa: TC001
from process.data_structure.water_usage_variables import WaterUseData  # noqa: TC001
from process.models.stellarator.initialization import st_init

logger = logging.getLogger(__name__)
if TYPE_CHECKING:
    from process.core.data_structure.base import DataStructure


initialise_later = object()


@dataclass(kw_only=True)
class DataStructure:
    """Dataclass holding the data structure"""

    water_use: WaterUseData = initialise_later
    costs_2015: Cost2015Data = initialise_later
    cs_fatigue: CSFatigueData = initialise_later
    vacuum: VacuumData = initialise_later
    costs: CostData = initialise_later
    first_wall: FirstWallData = initialise_later
    fwbs: FWBSData = initialise_later
    blanket: BlanketData = initialise_later
    structure: StructureData = initialise_later
    times: TimesData = initialise_later
    reinke: ReinkeData = initialise_later
    ccfe_hcpb: CCFEHCPBData = initialise_later
    pulse: PulseData = initialise_later
    build: BuildData = initialise_later
    primary_pumping: PrimaryPumpingData = initialise_later
    buildings: BuildingsData = initialise_later
    constraints: ConstraintData = initialise_later
    dcll: DCLLData = initialise_later
    current_drive: CurrentDriveData = initialise_later
    heat_transport: HeatTransportData = initialise_later
    ife: IFEData = initialise_later
    divertor: DivertorData = initialise_later
    pf_coil: PFCoilData = initialise_later
    power: PowerData = initialise_later
    stellarator: StellaratorData = initialise_later
    stellarator_config: StellaratorConfigData = initialise_later
    pf_power: PFPowerData = initialise_later
    neoclassics: NeoclassicsData = initialise_later
    impurity_radiation: ImpurityRadiationData = initialise_later
    physics: PhysicsData = initialise_later
    rebco: RebcoData = initialise_later
    tfcoil: TFData = initialise_later
    superconducting_tfcoil: SuperconductingTFData = initialise_later
    globals: GlobalData = initialise_later
    scan: ScanData = initialise_later
    numerics: NumericsData = initialise_later

    def __post_init__(self):
        """Set up the data structure"""
        for i, f in enumerate(fields(self)):
            if getattr(self, f.name) is initialise_later:
                # need to resolve type otherwise it is read as a str
                r = get_type_hints(DataStructure)
                ff = fields(DataStructure)[i]
                ftype = r[ff.name]
                setattr(self, f.name, ftype())

    def init_process(self, update_obsolete: bool = False):
        """Routine that calls the initialisation routines

        This routine calls the main initialisation routines that set
        the default values for the global variables, reads in data from
        the input file, and checks the run parameters for consistency.
        """
        # Initialise the program variables
        iteration_variables.initialise_iteration_variables(self)

        # Creating and open the files MFile and OUTFile
        process_output.OutputFileManager.open_files(self.globals.output_prefix)

        filename = self.globals.fileprefix
        # InDat reads in and updates obsolete variables if requested
        in_dat = InDat(filename=filename, update_obsolete=update_obsolete)

        inputs = self.parse_input_file(in_dat)

        # Set active constraints
        self.set_active_constraints()

        # set the device type (icase)
        self.set_device_type()

        # Check input data for errors/ambiguities
        check_process(self, inputs)

        self.write_outfile_header()

    def set_active_constraints(self):
        """Set constraints provided in the input file as 'active'"""
        num_constraints = 0
        for i in range(ConstraintManager.num_constraints()):
            if self.numerics.icc[i] != 0:
                self.numerics.active_constraints[self.numerics.icc[i] - 1] = True
                num_constraints += 1

        if self.numerics.n_equality_constraints < 0:
            # The value of n_equality_constraints has not been set in the input file.
            # Default = 0.
            self.numerics.n_equality_constraints = (
                num_constraints - self.numerics.n_inequality_constraints
            )
        else:
            self.numerics.n_inequality_constraints = (
                num_constraints - self.numerics.n_equality_constraints
            )

    def set_device_type(self):
        """Set the fusion device type on the data structure"""
        if self.ife.ife == 1:
            self.globals.icase = "Inertial Fusion model"
        elif self.stellarator.istell != StellaratorModel.DISABLED:
            self.globals.icase = "Stellarator model"
            st_init(self)

    def write_outfile_header(self):
        """Write a summary of the PROCESS run to the output file and MFile"""
        # Outfile and terminal #
        for outfile in [constants.NOUT, constants.IOTTY]:
            # PROCESS code header
            process_output.oblnkl(outfile)
            process_output.ostars(outfile, 110)
            process_output.ocentr(outfile, "PROCESS", 110)
            process_output.ocentr(outfile, "Power Reactor Optimisation Code", 110)
            process_output.ostars(outfile, 110)
            process_output.oblnkl(outfile)

            # Run execution details
            version = process.__version__
            process_output.ocmmnt(outfile, f"Version : {version}")

            git_branch, git_tag = get_git_summary()

            process_output.ocmmnt(outfile, f"Git Tag : {git_tag}")
            process_output.ocmmnt(outfile, f"Git Branch : {git_branch}")

            date_string = dt.datetime.now(dt.UTC).strftime("%d/%m/%Y %Z")
            time_string = dt.datetime.now(dt.UTC).strftime("%H:%M")

            process_output.ocmmnt(outfile, f"Date : {date_string}")
            process_output.ocmmnt(outfile, f"Time : {time_string}")

            user = getpass.getuser()
            machine = socket.gethostname()

            process_output.ocmmnt(outfile, f"User : {user}")
            process_output.ocmmnt(outfile, f"Computer : {machine}")
            process_output.ocmmnt(outfile, f"Directory : {Path.cwd()}")

            fileprefix = self.globals.fileprefix
            process_output.ocmmnt(
                outfile,
                f"Input : {fileprefix}",
            )
            runtitle = self.globals.runtitle
            process_output.ocmmnt(
                outfile,
                f"Run title : {runtitle}",
            )

            process_output.ocmmnt(
                outfile,
                f"Run type : Reactor concept design: {self.globals.icase},"
                " (c) UK Atomic Energy Authority",
            )

            process_output.oblnkl(outfile)
            process_output.ostars(outfile, 110)
            process_output.oblnkl(outfile)

            process_output.ocmmnt(
                outfile, f"Equality constraints : {self.numerics.n_equality_constraints}"
            )
            process_output.ocmmnt(
                outfile,
                f"Inequality constraints : {self.numerics.n_inequality_constraints}",
            )
            process_output.ocmmnt(
                outfile,
                f"Total constraints : "
                f"{self.numerics.n_inequality_constraints + self.numerics.n_equality_constraints}",  # noqa: E501
            )
            process_output.ocmmnt(
                outfile, f"Iteration variables : {self.numerics.n_iteration_variables}"
            )
            # If optimising, write objective function and convergence parameter
            if self.numerics.i_process_run_mode == PROCESSRunMode.OPTIMISATION:
                process_output.ocmmnt(
                    outfile,
                    f"Max iterations : {self.globals.maxcal}",
                )

                if self.numerics.i_figure_merit > 0:
                    figure_merit_string = "  -- minimise "
                    figure_merit_sign = "+"
                else:
                    figure_merit_string = "  -- maximise "
                    figure_merit_sign = "-"

                fom_string = FiguresOfMerit(
                    abs(self.numerics.i_figure_merit)
                ).description
                process_output.ocmmnt(
                    outfile,
                    f"Figure of merit : "
                    f"{figure_merit_sign}{abs(self.numerics.i_figure_merit)}"
                    f"{figure_merit_string}{fom_string}",
                )
                process_output.ocmmnt(
                    outfile,
                    f"Convergence parameter : {self.numerics.epsvmc}",
                )

            process_output.oblnkl(outfile)
            process_output.ostars(outfile, 110)

        # MFile #
        mfile = constants.MFILE

        process_output.ovarre(mfile, "PROCESS version", "(procver)", f'"{version}"')
        process_output.ovarre(mfile, "Date of run", "(date)", f'"{date_string}"')
        process_output.ovarre(mfile, "Time of run", "(time)", f'"{time_string}"')
        process_output.ovarre(mfile, "User", "(username)", f'"{user}"')
        process_output.ovarre(mfile, "PROCESS run title", "(runtitle)", f'"{runtitle}"')
        process_output.ovarre(mfile, "PROCESS git tag", "(tagno)", f'"{git_tag}"')
        process_output.ovarre(
            mfile, "PROCESS git branch", "(branch_name)", f'"{git_branch}"'
        )
        process_output.ovarre(mfile, "Input filename", "(fileprefix)", f'"{fileprefix}"')

        process_output.ovarre(
            mfile,
            "Optimisation switch",
            "(i_process_run_mode)",
            self.numerics.i_process_run_mode,
        )
        # If optimising, write figure of merit switch
        if self.numerics.i_process_run_mode == PROCESSRunMode.OPTIMISATION:
            process_output.ovarre(
                mfile,
                "Figure of merit switch",
                "(i_figure_merit)",
                self.numerics.i_figure_merit,
            )

    def parse_input_file(self, in_dat: InDat):
        """Create a dictionary of variables from the input file
        and set the values on the data structure

        Parameters
        ----------
        in_dat: InDat
            The input file data from InDat

        Raises
        ------
        ProcessValidationError
            If bounds are not of the correct form of boundl(index) = value

        """
        # These need to be set to zero so the iteration variables are not counted
        # multiple times (this would otherwise occur when using
        # get_mfile_initial_ixc_values())
        self.numerics.n_iteration_variables = 0
        self.numerics.n_constraints = 0

        variables = {}

        for variable_name, info in in_dat.data.items():
            # the upper and lower bounds already exist as boundu and
            # boundl arrays, so no need to process the bounds dictionary
            if variable_name == "bounds":
                continue
            variables = self.set_on_datastructure(
                variable_name=variable_name,
                info=info,
                variables=variables,
            )

        return variables

    def set_on_datastructure(self, variable_name, info, variables):
        """Set the values from the IN.DAT on the data structure, and store in
        variables dictionary

        Parameters
        ----------
        variable_name: str
            The variable's name
        info: INVariable
            Information of a single variable from the IN.DAT file
        variables: dict
            Dictionary of variables to be updated
        """
        variable_config = copy.copy(INPUT_VARIABLES.get(variable_name))
        # string indicates it should be set on the new object data structure
        if isinstance(variable_config.module, str):
            module = self
            for name in variable_config.module.split("."):
                module = getattr(module, name)

            variable_config.module = module
        variable_value = info.value

        # Validate the variable and also clean it (cast to correct type)
        # If the variable value (after the = sign) contains a ',' or is a list (len > 0)
        # then it defines the whole array so needs to be split down into its elements
        # and the parsed like an array defined as 'my_array(<index>) = <value>'
        if "," in variable_value:
            variable_value = list(variable_value.split(","))
        if len(np.shape(variable_value)) > 0 or "," in variable_value:
            clean_variable_value = [
                validate_variable(
                    variable_name,
                    v,
                    i,
                    variable_config,
                )
                for i, v in enumerate(variable_value, start=1)
            ]
        else:
            clean_variable_value = validate_variable(
                variable_name,
                variable_value,
                None,
                variable_config,
            )

        # perform any additional actions before setting values on the data structure
        if variable_config.additional_actions is not None:
            # intentionally passing the variable name as in the input file,
            # not the module target name.
            variable_config.additional_actions(
                variable_name,
                clean_variable_value,
                variable_config,
                self,
            )

        # check if the target name in the module is different to the variable name
        # in the input file
        variable_name_in_module = variable_config.target_name or variable_name

        if variable_config.set_variable:
            if isinstance(clean_variable_value, list):
                for idx, value in enumerate(clean_variable_value, start=1):
                    set_array_variable(
                        variable_name_in_module,
                        value,
                        idx,
                        variable_config,
                    )
            else:
                set_scalar_variable(
                    variable_name_in_module, clean_variable_value, variable_config
                )

            # add the variable to a dictionary indexed by the variable name
            # (in the input file)
            variables[variable_name] = {
                "value": clean_variable_value,
                "config": variable_config,
            }

        return variables


def validate_variable(
    name: str,
    value: str,
    array_index: int | None,
    config: InputVariable,
) -> ValidInputTypes:
    """Validate an input.

    Parameters
    ----------
    name :
        the name of the input.
    value :
        the value of the input variable.
    array_index :
        the array index of the variable in the input file.
    config :
        the config of the variable that describes how to validate and process it.

    Returns
    -------
    :
        the input value with the correct type.

    Raises
    ------
    ProcessValidationError
        Variable validation failure
    """
    # check that if the variable should be an array, then an array index is provided

    if array_index is None and config.array:
        error_msg = f"Expected '{name}' to be an array."
        raise ProcessValidationError(error_msg)

    if array_index is not None and not config.array:
        error_msg = f"Not expecting '{name}' to be an array."
        raise ProcessValidationError(error_msg)

    if config.type in {float, int} and isinstance(value, str):
        value = value.lower().replace("d", "e")

    try:
        clean_value = config.type(value)
    except ValueError as e:
        error_msg = (
            f"Cannot cast variable name '{name}' to a {config.type} (value = {value})"
        )
        raise ProcessValidationError(error_msg) from e

    if (
        config.range is not None
        and not config.range[0] <= clean_value <= config.range[1]
    ):
        error_msg = (
            f"Variable '{name}' is not on the prescribed range"
            f" {config.range} (value = {value})"
        )
        raise ProcessValidationError(error_msg)

    if config.choices is not None and clean_value not in config.choices:
        error_msg = f"Variable '{name}' is not one of {config.choices} (value = {value})"
        raise ProcessValidationError(error_msg)

    if config.additional_validation is not None:
        clean_value = config.additional_validation(name, clean_value, config)

    return clean_value


def set_scalar_variable(name: str, value: ValidInputTypes, config: InputVariable):
    """Set a scalar (not part of an array) variable in the `config.module`.

    Parameters
    ----------
    name :
        the name of the input.
    value :
        the value of the input variable.
    config :
        the config of the variable that describes how to validate and process it.

    Raises
    ------
    ProcessValueError
        variable not found
    """
    current_value = getattr(config.module, name, ...)

    # use ... sentinel because None is probably a valid return from Fortran
    # and definitely will be when moving to a Python data structure
    if current_value is ...:
        error_msg = (
            f"Fortran module '{config.module}' does not have a variable '{name}'."
        )
        raise ProcessValueError(error_msg)

    setattr(config.module, name, value)


def set_array_variable(name: str, value: str, array_index: int, config: InputVariable):
    """Set an array variable in the `config.module`.

    The way PROCESS input files are structured, each element of the array is provided on
    one line
    so this function just needs to set the `value` at `array_index` (-1)
    because of Fortran-based indexing.

    Parameters
    ----------
    name :
        the name of the input.
    value :
        the value of the input variable.
    array_index :
        the array index of the variable in the input file.
    config :
        the config of the variable that describes how to validate and process it.

    Raises
    ------
    ProcessValueError
        variable not found
    """
    current_array = getattr(config.module, name, ...)
    shape = current_array.shape

    # use ... sentinel for same reason as above
    if current_array is ...:
        error_msg = f"Data structure '{config.module}' does not have an array '{name}'."
        raise ProcessValueError(error_msg)

    if config.type is str:
        raise NotImplementedError

    new_array = copy.deepcopy(current_array)

    if len(shape) > 1:
        new_array = new_array.T.ravel()

    new_array[array_index - 1] = value

    if len(shape) > 1:
        new_array = new_array.reshape(shape, order="F")

    setattr(config.module, name, new_array)


def get_git_summary() -> tuple[str, str]:
    """Get git summary"""
    directory = Path(process.__file__).parent
    try:
        git_branch = (
            subprocess  # noqa: S602
            .run(
                "git rev-parse --abbrev-ref HEAD",  # noqa: S607
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

        git_tag = (
            subprocess  # noqa: S602
            .run(
                "git describe --tags",  # noqa: S607
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

    except (subprocess.CalledProcessError, AttributeError):
        return "", ""
    else:
        return git_branch, git_tag
