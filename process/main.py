"""Run Process by calling into the Fortran.

This uses a Python module called fortran.py, which uses an extension module
called "_fortran.cpython... .so", which are both generated from
process_module.f90. The process_module module contains the code to actually run
Process.

This file, process.py, is now analogous to process.f90, which contains the
Fortran "program" statement. This Python module effectively acts as the Fortran
"program".

Power Reactor Optimisation Code for Environmental and Safety Studies
P J Knight, CCFE, Culham Science Centre
J Morris, CCFE, Culham Science Centre

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

from typing import Protocol
from process import fortran
from process.buildings import Buildings
from process.costs import Costs
from process.io import plot_proc
from process.io import mfile
from process.plasma_geometry import PlasmaGeom
from process.pulse import Pulse
from process.scan import Scan
from process.stellarator import Stellarator
from process.structure import Structure
from process.build import Build
from process.utilities.f2py_string_patch import string_to_f2py_compatible
import argparse
from process.pfcoil import PFCoil
from process.tfcoil import TFcoil
from process.divertor import Divertor
from process.availability import Availability
from process.ife import IFE
from process.costs_2015 import Costs2015
from process.power import Power
from process.cs_fatigue import CsFatigue
from process.physics import Physics
from process.io import obsolete_vars as ov
from process.plasma_profiles import PlasmaProfile
from process.hcpb import CCFE_HCPB
from process.dcll import DCLL
from process.blanket_library import BlanketLibrary
from process.fw import Fw
from process.current_drive import CurrentDrive
from process.impurity_radiation import initialise_imprad
from process.caller import write_output_files


from pathlib import Path
import os
import logging

# For VaryRun
from process.io.process_config import RunProcessConfig
from process.io.process_funcs import (
    get_neqns_itervars,
    get_variable_range,
    check_input_error,
    process_stopped,
    no_unfeasible_mfile,
    vary_iteration_variables,
    process_warnings,
)
from process.vacuum import Vacuum
from process.water_use import WaterUse
from process.sctfcoil import Sctfcoil


os.environ["PYTHON_PROCESS_ROOT"] = os.path.join(os.path.dirname(__file__))

# Define parent logger
logger = logging.getLogger("process")
# Ensure every log goes through to a handler
logger.setLevel(logging.DEBUG)
# Handler for logging to stderr (and hence the terminal by default)
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.WARNING)
# Handler for logging to file
f_handler = logging.FileHandler("process.log", mode="w")
f_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
s_handler.setFormatter(formatter)
f_handler.setFormatter(formatter)
logger.addHandler(s_handler)
logger.addHandler(f_handler)


class Process:
    """The main Process class."""

    def __init__(self, args=None):
        """Run Process.

        :param args: Arguments to parse, defaults to None
        :type args: list, optional
        """
        self.parse_args(args)
        self.run_mode()
        self.post_process()

    def parse_args(self, args):
        """Parse the command-line arguments, such as the input filename.

        :param args: Arguments to parse
        :type args: list
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
        parser.add_argument("-p", "--plot", action="store_true", help="plot an mfile")
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

        # If args is not None, then parse the supplied arguments. This is likely
        # to come from the test suite when testing command-line arguments; the
        # method is being run from the test suite.
        # If args is None, then use actual command-line arguments (e.g.
        # sys.argv), as the method is being run from the command-line.
        self.args = parser.parse_args(args)
        # Store namespace object of the args

    def run_mode(self):
        """Determine how to run Process."""
        # Store run object: useful for testing
        if self.args.varyiterparams:
            self.run = VaryRun(self.args.varyiterparamsconfig, self.args.solver)
        else:
            self.run = SingleRun(self.args.input, self.args.solver)
        self.run.run()

    def post_process(self):
        """Perform post-run actions, like plotting the mfile."""
        # TODO Currently, Process will always run on an input file beforehand.
        # It would be better to not require this, so just plot_proc could be
        # run, for example.
        if self.args.plot:
            # Check mfile exists, then plot
            mfile_path = Path(self.args.mfile)
            mfile_str = str(mfile_path.resolve())
            if mfile_path.exists():
                # TODO Get --show arg to work: actually show the plot, don't
                # just save it
                plot_proc.main(args=["-f", mfile_str])
            else:
                logger.error("mfile to be used for plotting doesn't exist")
        if self.args.mfilejson:
            # Produce a json file containing mfile output, useful for VVUQ work.
            mfile_path = Path(self.args.mfile)
            mfile_data = mfile.MFile(filename=mfile_path)
            mfile_data.open_mfile()
            mfile_data.write_to_json()


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

    def __init__(self, config_file, solver="vmcon"):
        """Initialise and perform a VaryRun.

        :param config_file: config file for run parameters
        :type config_file: str
        :param solver: which solver to use, as specified in solver.py
        :type solver: str, optional
        """
        # Store the absolute path to the config file immediately: various
        # dir changes happen in old run_process code
        self.config_file = Path(config_file).resolve()
        self.solver = solver

    def run(self):
        """Perform a VaryRun by running multiple SingleRuns.

        :raises FileNotFoundError: if input file doesn't exist
        """
        # The input path for the varied input file
        input_path = self.config_file.parent / "IN.DAT"

        # Taken without much modification from the original run_process.py
        # Something changes working dir in config lines below
        config = RunProcessConfig(self.config_file)
        config.setup()

        fortran.init_module.init_all_module_vars()
        fortran.init_module.init()

        neqns, itervars = get_neqns_itervars()
        lbs, ubs = get_variable_range(itervars, config.factor)

        # If config file contains WDIR, use that. Otherwise, use the directory
        # containing the config file (used when running regression tests in
        # temp dirs)
        # TODO Not sure this is required any more
        if config.wdir:
            wdir = config.wdir
        else:
            wdir = Path(self.config_file).parent

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
                            "But finished anyway! {} ".format(no_unfeasible)
                        )
                    if process_warnings():
                        print(
                            "\nThere were warnings in the final PROCESS run. "
                            "Please check the log file!\n"
                        )
                    # This means success: feasible solution found
                    break
                else:
                    print(
                        "WARNING: {} non-feasible point(s) in sweep! "
                        "Rerunning!".format(no_unfeasible)
                    )
            else:
                print("PROCESS has stopped without finishing!")

            vary_iteration_variables(itervars, lbs, ubs)

        config.error_status2readme()


class SingleRun:
    """Perform a single run of PROCESS."""

    def __init__(self, input_file, solver="vmcon"):
        """Read input file and initialise variables.

        :param input_file: input file named <optional_name>IN.DAT
        :type input_file: str
        :param solver: which solver to use, as specified in solver.py
        :type solver: str, optional
        """
        self.input_file = input_file

        self.validate_input()
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
        self.run_tests()
        self.call_solver()
        self.run_scan(self.solver)
        self.finish()
        self.append_input()

    @staticmethod
    def init_module_vars():
        """Initialise all module variables in the Fortran.

        This "resets" all module variables to their initialised values, so each
        new run doesn't have any side-effects from previous runs.
        """
        fortran.init_module.init_all_module_vars()

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
                "Input file not found on this path. There " "is no input file named",
                self.input_file,
                "in the analysis " "folder",
            )

        # Set the input file in the Fortran
        fortran.global_variables.fileprefix = string_to_f2py_compatible(
            fortran.global_variables.fileprefix,
            str(self.input_path.resolve()),
            except_length=True,
        )

    def set_output(self):
        """Set the output file name.

        Set Path object on the Process object, and set the prefix in the Fortran.
        """
        self.output_path = Path(self.filename_prefix + "OUT.DAT")
        fortran.global_variables.output_prefix = string_to_f2py_compatible(
            fortran.global_variables.output_prefix, self.filename_prefix
        )

    def set_mfile(self):
        """Set the mfile filename."""
        self.mfile_path = Path(self.filename_prefix + "MFILE.DAT")

    @staticmethod
    def initialise():
        """Run the init module to call all initialisation routines."""
        initialise_imprad()
        # Reads in input file
        fortran.init_module.init()

        # Order optimisation parameters (arbitrary order in input file)
        # Ensures consistency and makes output comparisons more straightforward
        n = int(fortran.numerics.nvar)
        # [:n] as array always at max size: contains 0s
        fortran.numerics.ixc[:n].sort()

    def run_tests(self):
        """Run tests if required to by input file."""
        # TODO This would do better in a separate input validation module.
        if fortran.global_variables.run_tests == 1:
            fortran.main_module.runtests()

    def call_solver(self):
        """Call the equation solver (HYBRD)."""
        # If no HYBRD (non-optimisation) runs are required, return
        if (fortran.numerics.ioptimz > 0) or (fortran.numerics.ioptimz == -2):
            return
        else:
            # eqslv() has been temporarily commented out. Please see the comment
            # in fortran.function_evaluator.fcnhyb() for an explanation.
            # Original call:
            # self.ifail = fortran.main_module.eqslv()
            raise NotImplementedError(
                "HYBRD non-optimisation solver is not " "implemented"
            )

    def run_scan(self, solver):
        """Create scan object if required.

        :param solver: which solver to use, as specified in solver.py
        :type solver: str
        """
        if fortran.numerics.ioptimz == 1:
            # VMCON optimisation
            self.scan = Scan(self.models, solver)
        elif fortran.numerics.ioptimz == -2:
            # No optimisation: compute the output variables now
            # Get optimisation parameters x, evaluate models
            fortran.define_iteration_variables.loadxc()
            self.ifail = 6
            write_output_files(models=self.models, ifail=self.ifail)
            self.show_errors()
        else:
            raise ValueError(
                f"Invalid ioptimz value: {fortran.numerics.ioptimz}. Please "
                "select either 1 (optimise) or -2 (no optimisation)."
            )

    def show_errors(self):
        """Report all informational/error messages encountered."""
        fortran.error_handling.show_errors()

    def finish(self):
        """Run the finish subroutine to close files open in the Fortran.

        Files being handled by Fortran must be closed before attempting to
        write to them using Python, otherwise only parts are written.
        """
        fortran.init_module.finish()

    def append_input(self):
        """Append the input file to the output file and mfile."""
        # Read IN.DAT input file
        with open(self.input_path, "r", encoding="utf-8") as input_file:
            input_lines = input_file.readlines()

        # Append the input file to the output file
        with open(self.output_path, "a", encoding="utf-8") as output_file:
            output_file.writelines(input_lines)

        # Append the input file to the mfile
        with open(self.mfile_path, "a", encoding="utf-8") as mfile_file:
            mfile_file.write("***********************************************")
            mfile_file.writelines(input_lines)

    def validate_input(self):
        """Checks the input IN.DAT file for any obsolete variables in the OBS_VARS dict contained
        within obsolete_variables.py.
        Then will print out what the used obsolete variables are (if any) before continuing the proces run.
        """

        obsolete_variables = ov.OBS_VARS
        obsolete_vars_help_message = ov.OBS_VARS_HELP

        filename = self.input_file

        variables_in_in_dat = []
        with open(filename, "r") as file:
            for line in file:
                if line[0] == "*" or "=" not in line:
                    continue

                else:
                    sep = " "
                    variables = line.strip().split(sep, 1)[0]
                    variables_in_in_dat.append(variables)

        obs_vars_in_in_dat = []
        replace_hints = {}
        for var in variables_in_in_dat:
            if var in obsolete_variables:
                obs_vars_in_in_dat.append(var)
                replace_hints[var] = obsolete_variables.get(var)

        if len(obs_vars_in_in_dat) > 0:
            message = (
                "The IN.DAT file contains obsolete variables from the OBS_VARS dictionary. The obsolete variables in your IN.DAT file are: "
                f"{obs_vars_in_in_dat}. "
                "Either remove these or replace them with their updated variable names. "
            )
            for obs_var in obs_vars_in_in_dat:
                if replace_hints[obs_var] is None:
                    message += f"\n\n {obs_var} is an obsolete variable and needs to be removed. "
                else:
                    message += f"\n \n {obs_var} is an obsolete variable and needs to be replaced by {str(replace_hints[obs_var])}. "
                message += f"{obsolete_vars_help_message.get(obs_var, '')}"

            raise ValueError(message)

        else:
            print("The IN.DAT file does not contain any obsolete variables.")

    def validate_user_model(self):
        """Checks that a user-created model has been injected correctly

        Ensures that the corresponding model variable in Models is defined
        and that any relevant switches are set correctly.
        """
        # try and get costs model
        self.models.costs


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
        self.build = Build()
        self.sctfcoil = Sctfcoil()
        self.tfcoil = TFcoil(build=self.build, sctfcoil=self.sctfcoil)
        self.divertor = Divertor()
        self.structure = Structure()
        self.plasma_geom = PlasmaGeom()
        self.availability = Availability()
        self.buildings = Buildings()
        self.vacuum = Vacuum()
        self.water_use = WaterUse()
        self.pulse = Pulse()
        self.ife = IFE(availability=self.availability, costs=self.costs)
        self.plasma_profile = PlasmaProfile()
        self.fw = Fw()
        self.blanket_library = BlanketLibrary(fw=self.fw)
        self.ccfe_hcpb = CCFE_HCPB(blanket_library=self.blanket_library)
        self.current_drive = CurrentDrive()
        self.physics = Physics(
            plasma_profile=self.plasma_profile, current_drive=self.current_drive
        )
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
        )
        self.dcll = DCLL(blanket_library=self.blanket_library)

    @property
    def costs(self) -> CostsProtocol:
        if fortran.cost_variables.cost_model == 0:
            return self._costs_1990
        if fortran.cost_variables.cost_model == 1:
            return self._costs_2015
        if fortran.cost_variables.cost_model == 2:
            if self._costs_custom is not None:
                return self._costs_custom
            raise ValueError("Custom costs model not initialised")
        # Probably overkill but makes typing happy
        raise ValueError("Unknown costs model")

    @costs.setter
    def costs(self, value: CostsProtocol):
        self._costs_custom = value


def main(args=None):
    """Run Process.

    The args parameter is used to control command-line arguments when running
    tests. Optional args can be supplied by different tests, which are then
    used instead of command-line arguments by argparse. This allows testing of
    different command-line arguments from the test suite.

    :param args: Arguments to parse, defaults to None
    :type args: list, optional
    """
    Process(args)


if __name__ == "__main__":
    main()
