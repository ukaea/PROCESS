"""Run Process by calling into the Fortran.

This uses a Python module called fortran.py, which uses an extension module
called "_fortran.cpython... .so", which are both generated from
process_module.f90. The process_module module contains the code to actually run
Process.

This file, process.py, is now analogous to process.f90, which contains the
Fortran "program" statement. This Python module effectively acts as the Fortran
"program".

Power Reactor Optimisation Code for Environmental and Safety Studies

This is a systems code that evaluates various physics and
engineering aspects of a fusion power plant subject to given
constraints, and can optimise these parameters by minimising
or maximising a function of them, such as the fusion power or
cost of electricity.

This program is derived from the TETRA and STORAC codes produced by
Oak Ridge National Laboratory, Tennessee, USA. The main authors in
the USA were J.D.Galambos and P.C.Shipe.

The code was transferred to Culham Laboratory, Oxfordshire, UK, in
April 1992, and the physics models were updated by P.J.Knight to
include the findings of the Culham reactor studies documented in
Culham Report AEA FUS 172 (1992). The standard of the Fortran has
been thoroughly upgraded since that time, and a number of additional
models have been added.

During 2012, PROCESS was upgraded from FORTRAN 77 to Fortran 95,
to facilitate the restructuring of the code into proper modules
(with all the benefits that modern software practices bring), and to
aid the inclusion of more advanced physics and engineering models under
development as part of a number of EFDA-sponsored collaborations.

Box file F/RS/CIRE5523/PWF (up to 15/01/96)
Box file F/MI/PJK/PROCESS and F/PL/PJK/PROCESS (15/01/96 to 24/01/12)
Box file T&amp;M/PKNIGHT/PROCESS (from 24/01/12)
"""

import argparse
import logging
import os
from pathlib import Path
from typing import Any, Protocol

import process
import process.data_structure as data_structure
import process.init as init
from process import constants
from process.io import (
    mfile,
    plot_plotly_sankey,
    plot_proc,
)
from process.io import obsolete_vars as ov

# For VaryRun
from process.io.process_config import RunProcessConfig
from process.io.process_funcs import (
    check_input_error,
    get_neqns_itervars,
    get_variable_range,
    no_unfeasible_mfile,
    process_stopped,
    process_warnings,
    vary_iteration_variables,
)
from process.log import logging_model_handler, show_errors
from process.models.availability import Availability
from process.models.blankets.blanket_library import BlanketLibrary
from process.models.blankets.dcll import DCLL
from process.models.blankets.hcpb import CCFE_HCPB
from process.models.build import Build
from process.models.buildings import Buildings
from process.models.costs.costs import Costs
from process.models.costs.costs_2015 import Costs2015
from process.models.cryostat import Cryostat
from process.models.cs_fatigue import CsFatigue
from process.models.divertor import Divertor
from process.models.fw import FirstWall
from process.models.ife import IFE
from process.models.pfcoil import PFCoil
from process.models.physics.current_drive import (
    CurrentDrive,
    ElectronBernstein,
    ElectronCyclotron,
    IonCyclotron,
    LowerHybrid,
    NeutralBeam,
)
from process.models.physics.impurity_radiation import initialise_imprad
from process.models.physics.physics import (
    DetailedPhysics,
    Physics,
    PlasmaBeta,
    PlasmaInductance,
)
from process.models.physics.plasma_geometry import PlasmaGeom
from process.models.physics.plasma_profiles import PlasmaProfile
from process.models.power import Power
from process.models.pulse import Pulse
from process.models.shield import Shield
from process.models.structure import Structure
from process.models.tfcoil.resistive_tf_coil import (
    AluminiumTFCoil,
    CopperTFCoil,
    ResistiveTFCoil,
)
from process.models.tfcoil.superconducting_tf_coil import SuperconductingTFCoil
from process.models.tfcoil.tf_coil import TFCoil
from process.process_output import OutputFileManager, oheadr
from process.scan import Scan
from process.stellarator.neoclassics import Neoclassics
from process.stellarator.stellarator import Stellarator
from process.vacuum import Vacuum, VacuumVessel
from process.water_use import WaterUse

os.environ["PYTHON_PROCESS_ROOT"] = os.path.join(os.path.dirname(__file__))

PACKAGE_LOGGING = True
"""Can be set False to disable package-level logging, e.g. in the test suite"""
logger = logging.getLogger("process")


class Process:
    """The main Process class."""

    def __init__(self, args: list[Any] | None = None):
        """Run Process.

        :param args: Arguments to parse, defaults to None
        """
        self.parse_args(args)
        self.run_mode()
        self.post_process()

    def parse_args(self, args: list[Any] | None):
        """Parse the command-line arguments, such as the input filename.

        Parameters
        ----------
        args :
            Arguments to parse
        """
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(
                "PROCESS\n"
                "Power Reactor Optimisation Code\n"
                "Copyright (c) [2023] [United Kingdom Atomic Energy Authority]\n"
                "\n"
                "Contact\n"
                "James Morris  : james.morris2@ukaea.uk\n"
                "Jonathan Maddock : jonathan.maddock@ukaea.uk\n"
                "\n"
                "GitHub        : https://github.com/ukaea/PROCESS\n"
            ),
        )

        # Optional args
        parser.add_argument(
            "-i",
            "--input",
            default="IN.DAT",
            metavar="input_file_path",
            type=str,
            help="The path to the input file that Process runs on",
        )
        parser.add_argument(
            "-s",
            "--solver",
            default="vmcon",
            metavar="solver_name",
            type=str,
            help="Specify which solver to use: only 'vmcon' at the moment",
        )
        parser.add_argument(
            "-v",
            "--varyiterparams",
            action="store_true",
            help="Vary iteration parameters",
        )
        parser.add_argument(
            "-c",
            "--varyiterparamsconfig",
            metavar="config_file",
            default="run_process.conf",
            help="configuration file for varying iteration parameters",
        )
        parser.add_argument(
            "-m",
            "--mfile",
            default="MFILE.DAT",
            help="mfile for post-processing/plotting",
        )
        parser.add_argument(
            "-mj",
            "--mfilejson",
            action="store_true",
            help="Produce a filled json from --mfile arg in working dir",
        )
        parser.add_argument(
            "--version",
            action="store_true",
            help="Print the version of PROCESS to the terminal",
        )
        parser.add_argument(
            "--update-obsolete",
            action="store_true",
            help="Automatically update obsolete variables in the IN.DAT file",
        )
        parser.add_argument(
            "--full-output",
            action="store_true",
            help="Run all summary plotting scripts for the output",
        )

        # If args is not None, then parse the supplied arguments. This is likely
        # to come from the test suite when testing command-line arguments; the
        # method is being run from the test suite.
        # If args is None, then use actual command-line arguments (e.g.
        # sys.argv), as the method is being run from the command-line.
        self.args = parser.parse_args(args)
        # Store namespace object of the args

    def run_mode(self):
        """Determine how to run Process."""
        if self.args.version:
            print(process.__version__)
            return
        # Store run object: useful for testing
        if self.args.varyiterparams:
            self.run = VaryRun(self.args.varyiterparamsconfig, self.args.solver)
        else:
            self.run = SingleRun(
                self.args.input,
                self.args.solver,
                update_obsolete=self.args.update_obsolete,
            )
        self.run.run()

    def post_process(self):
        """Perform post-run actions, like plotting the mfile."""
        # TODO Currently, Process will always run on an input file beforehand.
        # It would be better to not require this, so just plot_proc could be
        # run, for example.
        if self.args.mfilejson:
            # Produce a json file containing mfile output, useful for VVUQ work.
            mfile_path = Path(self.args.mfile)
            mfile_data = mfile.MFile(filename=mfile_path)
            mfile_data.open_mfile()
            mfile_data.write_to_json()
        if self.args.full_output:
            # Run all summary plotting scripts for the output
            mfile_path = Path(str(self.args.input).replace("IN.DAT", "MFILE.DAT"))
            mfile_str = str(mfile_path.resolve())
            print(f"Plotting mfile {mfile_str}")
            if mfile_path.exists():
                plot_proc.main(args=["-f", mfile_str])
                plot_plotly_sankey.main(args=["-m", mfile_str])

            else:
                logger.error("mfile to be used for plotting doesn't exist")


class VaryRun:
    """Vary iteration parameters until a solution is found.

    This is the old run_process.py utility.

    Code to run PROCESS with a variation of the iteration parameters
    until a feasible solution is found.
    If running in sweep mode, the allowed number of unfeasible solutions
    can be changed in the config file.

    Input files:
    run_process.conf (config file, in the same directory as this file)
    An IN.DAT file as specified in the config file

    Output files:
    All of them in the work directory specified in the config file
    OUT.DAT     -  PROCESS output
    MFILE.DAT   -  PROCESS output
    process.log - logfile of PROCESS output to stdout
    README.txt  - contains comments from config file
    """

    def __init__(self, config_file: str, solver: str = "vmcon"):
        """Initialise and perform a VaryRun.

        Parameters
        ----------
        config_file:
            config file for run parameters
        solver:
            which solver to use, as specified in solver.py
        """
        # Store the absolute path to the config file immediately: various
        # dir changes happen in old run_process code
        self.config_file = Path(config_file).resolve()
        self.solver = solver

    def run(self):
        """Perform a VaryRun by running multiple SingleRuns.

        Raises
        ------
        FileNotFoundError
            if input file doesn't exist
        """
        # The input path for the varied input file
        input_path = self.config_file.parent / "IN.DAT"

        # Taken without much modification from the original run_process.py
        # Something changes working dir in config lines below
        config = RunProcessConfig(self.config_file)
        config.setup()

        setup_loggers(Path(config.wdir) / "process.log")

        init.init_all_module_vars()
        init.init_process()

        _neqns, itervars = get_neqns_itervars()
        lbs, ubs = get_variable_range(itervars, config.factor)

        # If config file contains WDIR, use that. Otherwise, use the directory
        # containing the config file (used when running regression tests in
        # temp dirs)
        # TODO Not sure this is required any more
        wdir = config.wdir or Path(self.config_file).parent

        # Check IN.DAT exists
        if not input_path.exists():
            raise FileNotFoundError

        # TODO add diff ixc summary part
        for i in range(config.niter):
            print(i, end=" ")

            # Run single runs (SingleRun()) of process as subprocesses. This
            # is the only way to deal with Fortran "stop" statements when
            # running VaryRun(), which otherwise cause the Python
            # interpreter to exit, when we want to vary the parameters and
            # run again
            # TODO Don't do this; remove stop statements from Fortran and
            # handle error codes
            # Run process on an IN.DAT file
            config.run_process(input_path, self.solver)

            check_input_error(wdir=wdir)

            if not process_stopped():
                no_unfeasible = no_unfeasible_mfile()
                if no_unfeasible <= config.no_allowed_unfeasible:
                    if no_unfeasible > 0:
                        print(
                            "WARNING: Non feasible point(s) in sweep, "
                            f"But finished anyway! {no_unfeasible} "
                        )
                    if process_warnings():
                        print(
                            "\nThere were warnings in the final PROCESS run. "
                            "Please check the log file!\n"
                        )
                    # This means success: feasible solution found
                    break
                print(
                    f"WARNING: {no_unfeasible} non-feasible point(s) in sweep! Rerunning!"
                )
            else:
                print("PROCESS has stopped without finishing!")

            vary_iteration_variables(itervars, lbs, ubs, config.generator)

        config.error_status2readme()


class SingleRun:
    """Perform a single run of PROCESS."""

    def __init__(
        self, input_file: str, solver: str = "vmcon", *, update_obsolete: bool = False
    ):
        """Read input file and initialise variables.

        Parameters
        ----------
        input_file:
            input file named <optional_name>IN.DAT
        solver:
            which solver to use, as specified in solver.py
        """
        self.input_file = input_file

        self.validate_input(update_obsolete)
        self.init_module_vars()
        self.set_filenames()
        self.initialise()
        self.models = Models()
        self.solver = solver

    def run(self):
        """Run PROCESS

        This is separate from init to allow model instances to be modified before a run.
        """
        self.validate_user_model()
        self.run_scan()
        self.finish()
        self.append_input()

    @staticmethod
    def init_module_vars():
        """Initialise all module variables in the Fortran.

        This "resets" all module variables to their initialised values, so each
        new run doesn't have any side-effects from previous runs.
        """
        init.init_all_module_vars()

    def set_filenames(self):
        """Validate the input filename and create other filenames from it."""
        self.set_input()
        self.set_output()
        self.set_mfile()

    def set_input(self):
        """Validate and set the input file path."""
        # Check input file ends in "IN.DAT", then save prefix
        # (the part before the IN.DAT)
        if self.input_file[-6:] != "IN.DAT":
            raise ValueError("Input filename must end in IN.DAT.")

        self.filename_prefix = self.input_file[:-6]

        # Check input file exists (path specified as CLI argument)
        input_path = Path(self.input_file)
        if input_path.exists():
            self.input_path = input_path
            # Set input as Path object
        else:
            print("-- Info -- run `process --help` for usage")
            raise FileNotFoundError(
                "Input file not found on this path. There is no input file named",
                self.input_file,
                "in the analysis folder",
            )

        # Set the input file in the Fortran
        data_structure.global_variables.fileprefix = str(self.input_path.resolve())

    def set_output(self):
        """Set the output file name.

        Set Path object on the Process object, and set the prefix in the Fortran.
        """
        self.output_path = Path(self.filename_prefix + "OUT.DAT")
        data_structure.global_variables.output_prefix = self.filename_prefix

    def set_mfile(self):
        """Set the mfile filename."""
        self.mfile_path = Path(self.filename_prefix + "MFILE.DAT")

    def initialise(self):
        """Run the init module to call all initialisation routines."""
        setup_loggers(
            Path(self.output_path.as_posix().replace("OUT.DAT", "process.log"))
        )

        initialise_imprad()
        # Reads in input file
        init.init_process()

        # Order optimisation parameters (arbitrary order in input file)
        # Ensures consistency and makes output comparisons more straightforward
        n = int(data_structure.numerics.nvar)
        # [:n] as array always at max size: contains 0s
        data_structure.numerics.ixc[:n].sort()

    def run_scan(self):
        """Create scan object if required."""
        # TODO Move this solver logic up to init?
        # ioptimz == 1: optimisation
        if data_structure.numerics.ioptimz == 1:
            pass
        # ioptimz == -2: evaluation
        elif data_structure.numerics.ioptimz == -2:
            # No optimisation: solve equality (consistency) constraints only using fsolve (HYBRD)
            self.solver = "fsolve"
        else:
            raise ValueError(
                f"Invalid ioptimz value: {data_structure.numerics.ioptimz}. Please "
                "select either 1 (optimise) or -2 (no optimisation)."
            )
        self.scan = Scan(self.models, self.solver)

    def show_errors(self):
        """Report all informational/error messages encountered."""
        show_errors(constants.NOUT)

    def finish(self):
        """Run the finish subroutine to close files open in the Fortran.

        Files being handled by Fortran must be closed before attempting to
        write to them using Python, otherwise only parts are written.
        """
        oheadr(constants.NOUT, "End of PROCESS Output")
        oheadr(constants.IOTTY, "End of PROCESS Output")
        oheadr(constants.NOUT, "Copy of PROCESS Input Follows")
        OutputFileManager.finish()

    def append_input(self):
        """Append the input file to the output file and mfile."""
        # Read IN.DAT input file
        with open(self.input_path, encoding="utf-8") as input_file:
            input_lines = input_file.readlines()

        # Append the input file to the output file
        with open(self.output_path, "a", encoding="utf-8") as output_file:
            output_file.writelines(input_lines)

        # Append the input file to the mfile
        with open(self.mfile_path, "a", encoding="utf-8") as mfile_file:
            mfile_file.write("***********************************************")
            mfile_file.writelines(input_lines)

    def validate_input(self, replace_obsolete: bool = False):
        """Checks the input IN.DAT file for any obsolete variables in the OBS_VARS dict contained
        within obsolete_variables.py. If obsolete variables are found, and if `replace_obsolete`
        is set to True, they are either removed or replaced by their updated names as specified
        in the OBS_VARS dictionary.
        """

        obsolete_variables = ov.OBS_VARS
        obsolete_vars_help_message = ov.OBS_VARS_HELP

        filename = self.input_file
        variables_in_in_dat = []
        modified_lines = []
        changes_made = []  # To store details of the changes

        with open(filename) as file:
            for line in file:
                # Skip comment lines or lines without an assignment
                if line.startswith("*") or "=" not in line:
                    modified_lines.append(line)
                    continue

                # Extract the variable name before the separator
                raw_variable_name = line.split("=", 1)[0].strip()
                # handle cases where the variable name might have parentheses
                variable_name = (
                    raw_variable_name.split("(", 1)[0]
                    if "(" in raw_variable_name
                    else raw_variable_name
                )

                # Check if the variable is obsolete and needs replacing
                if variable_name in obsolete_variables:
                    replacement = obsolete_variables.get(variable_name)
                    if replace_obsolete:
                        # Prepare replacement or removal
                        if replacement is None:
                            # If no replacement is defined, comment out the line
                            modified_lines.append(f"* Obsolete: {line}")
                            changes_made.append(
                                f"Commented out obsolete variable: {variable_name}"
                            )
                        else:
                            if isinstance(replacement, list):
                                # Raise an error if replacement is a list
                                replacement_str = ", ".join(replacement)
                                raise ValueError(
                                    f"The variable '{variable_name}' is obsolete and should be replaced by the following variables: {replacement_str}. "
                                    "Please set their values accordingly."
                                )
                            # Replace obsolete variable
                            modified_line = line.replace(variable_name, replacement, 1)
                            modified_lines.append(
                                f"* Replaced '{variable_name}' with '{replacement}'\n{modified_line}"
                            )
                            changes_made.append(
                                f"Replaced '{variable_name}' with '{replacement}'"
                            )
                            variables_in_in_dat.append(variable_name)
                    else:
                        # If replacement is False, add the line as-is
                        modified_lines.append(line)
                else:
                    modified_lines.append(line)

        obs_vars_in_in_dat = [
            var for var in variables_in_in_dat if var in obsolete_variables
        ]
        if obs_vars_in_in_dat:
            if replace_obsolete:
                # If replace_obsolete is True, write the modified content to the file
                with open(filename, "w") as file:
                    file.writelines(modified_lines)
                print(
                    "The IN.DAT file has been updated to replace or comment out obsolete variables."
                )
                print("Summary of changes made:")
                for change in changes_made:
                    print(f" - {change}")
            else:
                # Only print the report if replace_obsolete is False
                message = (
                    "The IN.DAT file contains obsolete variables from the OBS_VARS dictionary. "
                    f"The obsolete variables in your IN.DAT file are: {obs_vars_in_in_dat}. "
                    "Either remove these or replace them with their updated variable names. "
                )
                for obs_var in obs_vars_in_in_dat:
                    replacement = obsolete_variables.get(obs_var)
                    if replacement is None:
                        message += f"\n\n{obs_var} is an obsolete variable and needs to be removed."
                    else:
                        message += f"\n\n{obs_var} is an obsolete variable and needs to be replaced by {replacement}."
                    message += f" {obsolete_vars_help_message.get(obs_var, '')}"
                raise ValueError(message)

        else:
            print("The IN.DAT file does not contain any obsolete variables.")

    def validate_user_model(self):
        """Checks that a user-created model has been injected correctly

        Ensures that the corresponding model variable in Models is defined
        and that any relevant switches are set correctly.
        """
        # try and get costs model
        try:
            _tmp = self.models.costs
        except ValueError as err:
            raise ValueError("User-created model not injected correctly") from err


class CostsProtocol(Protocol):
    """Protocol layout for costs models"""

    def run(self):
        """Run the model"""

    def output(self):
        """write model output"""


class Models:
    """Creates instances of physics and engineering model classes.

    Creates objects to interface with corresponding Fortran physics and
    engineering modules.
    """

    def __init__(self):
        """Create physics and engineering model objects.

        This also initialises module variables in the Fortran for that module.
        """
        self._costs_custom = None
        self._costs_1990 = Costs()
        self._costs_2015 = Costs2015()
        self.cs_fatigue = CsFatigue()
        self.pfcoil = PFCoil(cs_fatigue=self.cs_fatigue)
        self.power = Power()
        self.cryostat = Cryostat()
        self.build = Build()
        self.sctfcoil = SuperconductingTFCoil()
        self.tfcoil = TFCoil(build=self.build)
        self.resistive_tf_coil = ResistiveTFCoil()
        self.copper_tf_coil = CopperTFCoil()
        self.aluminium_tf_coil = AluminiumTFCoil()
        self.divertor = Divertor()
        self.structure = Structure()
        self.plasma_geom = PlasmaGeom()
        self.availability = Availability()
        self.buildings = Buildings()
        self.vacuum = Vacuum()
        self.vacuum_vessel = VacuumVessel()
        self.water_use = WaterUse()
        self.pulse = Pulse()
        self.shield = Shield()
        self.ife = IFE(availability=self.availability, costs=self.costs)
        self.plasma_profile = PlasmaProfile()
        self.fw = FirstWall()
        self.blanket_library = BlanketLibrary(fw=self.fw)
        self.ccfe_hcpb = CCFE_HCPB(fw=self.fw)
        self.current_drive = CurrentDrive(
            plasma_profile=self.plasma_profile,
            electron_cyclotron=ElectronCyclotron(plasma_profile=self.plasma_profile),
            ion_cyclotron=IonCyclotron(plasma_profile=self.plasma_profile),
            lower_hybrid=LowerHybrid(plasma_profile=self.plasma_profile),
            neutral_beam=NeutralBeam(plasma_profile=self.plasma_profile),
            electron_bernstein=ElectronBernstein(plasma_profile=self.plasma_profile),
        )
        self.plasma_beta = PlasmaBeta()
        self.plasma_inductance = PlasmaInductance()
        self.physics = Physics(
            plasma_profile=self.plasma_profile,
            current_drive=self.current_drive,
            plasma_beta=self.plasma_beta,
            plasma_inductance=self.plasma_inductance,
        )
        self.physics_detailed = DetailedPhysics(
            plasma_profile=self.plasma_profile,
        )
        self.neoclassics = Neoclassics()
        if data_structure.stellarator_variables.istell != 0:
            self.stellarator = Stellarator(
                availability=self.availability,
                buildings=self.buildings,
                vacuum=self.vacuum,
                costs=self.costs,
                power=self.power,
                plasma_profile=self.plasma_profile,
                hcpb=self.ccfe_hcpb,
                current_drive=self.current_drive,
                physics=self.physics,
                neoclassics=self.neoclassics,
                plasma_beta=self.plasma_beta,
            )

        self.dcll = DCLL(fw=self.fw)

    @property
    def costs(self) -> CostsProtocol:
        if data_structure.cost_variables.cost_model == 0:
            return self._costs_1990
        if data_structure.cost_variables.cost_model == 1:
            return self._costs_2015
        if data_structure.cost_variables.cost_model == 2:
            if self._costs_custom is not None:
                return self._costs_custom
            raise ValueError("Custom costs model not initialised")
        # Probably overkill but makes typing happy
        raise ValueError("Unknown costs model")

    @costs.setter
    def costs(self, value: CostsProtocol):
        self._costs_custom = value


# setup handlers for writing to terminal (on warnings+)
# or writing to the log file (on info+)
logging_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
logging_stream_handler = logging.StreamHandler()
logging_stream_handler.setLevel(logging.CRITICAL)
logging_stream_handler.setFormatter(logging_formatter)

logging_file_handler = logging.FileHandler("process.log", mode="a")
logging_file_handler.setLevel(logging.INFO)
logging_file_handler.setFormatter(logging_formatter)

logging_model_handler.setLevel(logging.WARNING)
logging_model_handler.setFormatter(logging_formatter)


def setup_loggers(working_directory_log_path: Path | None = None):
    """A function that adds our handlers to the appropriate logger object.

    Parameters
    ----------
    working_directory_log_path: Path | None :
         (Default value = None)
    """
    # Remove all of the existing handlers from the 'process' package logger

    logger.handlers.clear()

    # we always want to add this handler because otherwise PROCESS' error
    # handling system won't work properly
    logger.addHandler(logging_model_handler)

    if not PACKAGE_LOGGING:
        return

    # (Re)add the loggers to the 'process' package logger (and its children)
    logger.addHandler(logging_stream_handler)
    logger.addHandler(logging_file_handler)

    if working_directory_log_path is not None:
        logging_file_input_location_handler = logging.FileHandler(
            working_directory_log_path.as_posix(), mode="w"
        )
        logging_file_input_location_handler.setLevel(logging.INFO)
        logging_file_input_location_handler.setFormatter(logging_formatter)
        logger.addHandler(logging_file_input_location_handler)


def main(args: list[Any] | None = None):
    """Run Process.

    The args parameter is used to control command-line arguments when running
    tests. Optional args can be supplied by different tests, which are then
    used instead of command-line arguments by argparse. This allows testing of
    different command-line arguments from the test suite.

    Parameters
    ----------
    args :
        Arguments to parse, defaults to None
    """

    Process(args)


if __name__ == "__main__":
    main()
