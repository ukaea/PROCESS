"""
Author: Hanni Lux (Hanni.Lux@ccfe.ac.uk)

Interfaces for Configuration values for programs
- run_process.py
- test_process.py
- evaluate_uncertainties.py

Compatible with PROCESS version 382

24/11/2021: Global dictionary variables moved within the functions
            to avoid cyclic dependencies. This is because the dicts
            generation script imports, and inspects, process.
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from sys import stderr
from time import sleep
from typing import ClassVar

from numpy import argsort, argwhere, logical_or
from numpy.random import default_rng

from process.io.configuration import Config
from process.io.in_dat import InDat
from process.io.mfile import MFile
from process.io.process_funcs import (
    check_in_dat,
    get_from_indat_or_default,
    set_variable_in_indat,
)
from process.io.python_fortran_dicts import get_dicts

logger = logging.getLogger(__name__)


class ProcessConfig:
    """
    Configuration parameters for PROCESS runs

    filename - Configuration file name

    wdir     - Working directory

    or_in_dat  - Original IN.DAT file

    process  - PROCESS binary

    niter    - (Maximum) number of iterations

    u_seed   - User specified seed value for the random number generator

    factor   - Multiplication factor adjusting the range in which the original
               iteration variables should get varied

    comment  - additional comment to be written into README.txt

    """

    filename = None
    configfileexists = True
    wdir = "."
    or_in_dat = "IN.DAT"
    process = "process.exe"
    niter = 10
    u_seed = None
    factor = 1.5
    comment = " "

    def echo(self):
        """echos the attributes of the class"""

        if self.wdir != ".":
            print(f"Working directory:   {self.wdir}")
        print(f"Original IN.DAT:     {self.or_in_dat}")
        print(f"PROCESS binary:      {self.process}")
        print(f"Number of iterations {self.niter}")

        if self.u_seed is not None:
            print(f"random seed          {self.u_seed}")
        print(f"variable range factor {self.factor}")
        if self.filename is not None:
            print(f"Config file          {self.filename}")
        if self.comment != "":
            print(f"Comment  {self.comment}")

    def prepare_wdir(self):
        """prepares the working directory"""

        try:
            os.stat(self.wdir)
        except FileNotFoundError:
            os.mkdir(self.wdir)

        if self.or_in_dat != "IN.DAT" or self.wdir != ".":
            subprocess.call(["cp", self.or_in_dat, self.wdir + "/IN.DAT"])
        else:
            subprocess.call(["cp", self.or_in_dat, "Input_IN.DAT"])

        if self.configfileexists:
            subprocess.call(["cp", self.filename, self.wdir])
        os.chdir(self.wdir)
        subprocess.call(
            [
                "rm -f OUT.DAT MFILE.DAT README.txt\
        SolverTest.out process.log *.pdf  uncertainties.nc time.info"
            ],
            shell=True,
        )

    def create_readme(self, directory="."):
        """creates README.txt containing comment"""

        if self.comment != "":
            with open(directory + "/README.txt", "w") as readme:
                readme.write(self.comment)

    def error_status2readme(self, directory="."):
        """appends PROCESS outcome to README.txt"""

        if os.path.isfile("MFILE.DAT"):
            m_file = MFile(filename=directory + "/MFILE.DAT")

            error_status = (
                f"Error status: {m_file.data['error_status'].get_scan(-1)}  "
                f"Error ID: {m_file.data['error_id'].get_scan(-1)}\n"
            )

            if self.comment != "":
                with open(directory + "/README.txt", "a") as readme:
                    readme.write(error_status)
            else:
                with open(directory + "/README.txt", "w") as readme:
                    readme.write(error_status)

    def modify_in_dat(self):
        """modifies the original IN.DAT file"""

    def setup(self):
        """sets up the program for running"""

        self.echo()

        self.prepare_wdir()

        self.create_readme()

        self.modify_in_dat()

        check_in_dat()

        self.generator = default_rng(seed=self.u_seed)

    def get_comment(self):
        """gets the comment line from the configuration file"""

        if not self.configfileexists:
            return False

        try:
            with open(self.filename) as configfile:
                for line in configfile:
                    condense = line.replace(" ", "")
                    condense = condense.rstrip()
                    lcase = condense.lower()
                    if (
                        len(condense) > 0
                        and (condense[0] != "*")
                        and lcase[:7] == "comment"
                    ):
                        self.comment = line[line.find("=") + 1 :]
                        return True

                self.configfileexists = False
                return False
        except FileNotFoundError:
            print(f"Error: No config file named {self.filename}", file=stderr)

        return False

    def get_attribute(self, attributename):
        """gets class attribute from configuration file"""

        if not self.configfileexists:
            return None

        try:
            with open(self.filename) as configfile:
                for line in configfile:
                    condense = line.replace(" ", "")
                    if (condense[0] != "*") and len(condense) > 1 and "=" in line:
                        varname = line[: line.find("=")]
                        varname = varname.replace(" ", "")
                        varname = varname.upper()
                        auxvar = line[line.find("=") + 1 :]
                        auxvar = auxvar.replace(" ", "")
                        auxvar = auxvar.rstrip()
                        if varname == attributename.upper() and auxvar == "":
                            return None
                        if varname == attributename.upper() and auxvar != "":
                            return auxvar
        except FileNotFoundError:
            print(f"Error: No config file named {self.filename}", file=stderr)
            self.configfileexists = False
            return None

        return None

    def set_attributes(self):
        """sets the attributes of the class"""
        # TODO This Path/str filename type confusion is a mess and must be
        # sorted out once a run_process regression test is running. Using Paths
        # or strs to decide if we're running a test or not is a bad idea, but
        # works in the short term to get a test running.

        # If self.filename is an instance of Path, this is being run from the
        # test suite in a temp dir. Set the working dir and input filename
        # accordingly
        if isinstance(self.filename, Path):
            self.wdir = str(self.filename.parent)
            # TODO It appears that the input file for run_process shouldn't be
            # called IN.DAT. Calling it ref_IN.DAT for now
            original_in_dat = self.get_attribute("original_in_dat") or "ref_IN.DAT"
            self.or_in_dat = str(self.filename.parent / original_in_dat)
        else:
            # run_process is not being run by the test suite
            # Get working dir and input files from run_process.conf
            buf = self.get_attribute("wdir")
            if buf is not None:
                self.wdir = buf

            buf = self.get_attribute("ORIGINAL_IN_DAT")
            if buf is not None:
                self.or_in_dat = buf

        try:
            with open(self.or_in_dat) as indatfile:
                indatfile.close()
        except FileNotFoundError:
            print(
                f"Error: {self.or_in_dat} does not exist! Create file or modify config file!",
                file=stderr,
            )
            raise

        buf = self.get_attribute("process")
        if buf is not None:
            self.process = buf

        buf = self.get_attribute("niter")
        if buf is not None:
            self.niter = int(buf)

        buf = self.get_attribute("seed")
        if buf is not None:
            if buf == "None":
                self.u_seed = None
            else:
                self.u_seed = int(buf)

        buf = self.get_attribute("factor")
        if buf is not None:
            self.factor = float(buf)

        if not self.get_comment():
            print(f"No comment in config file {self.filename}")

    def run_process(self, input_path, solver="vmcon"):
        """Perform a single run of PROCESS, catching any errors.

        :param input_path: the input file to run on
        :type input_path: pathlib.Path
        :param solver: which solver to use, as specified in solver.py, defaults
        to "vmcon"
        :type solver: str, optional
        """
        with open("process.log", "w") as logfile:
            print("PROCESS run started ...", end="")

            try:
                # Run process on an IN.DAT file
                subprocess.check_call(
                    ["process", "-i", str(input_path), "--solver", solver],
                    stdout=logfile,
                    stderr=logfile,
                )
            except subprocess.CalledProcessError as err:
                # Process has exited with a non-zero exit code. This is likely to
                # be due to hitting a "stop 1" in the Fortran. Catch this exception
                # to allow execution to continue without exiting
                print(
                    "Error: There was a problem with the PROCESS execution!",
                    err,
                    file=sys.stderr,
                )
                print("Refer to the logfile for more information!", file=sys.stderr)

        print("finished.")


################################################################################
# class TestProcessConfig(ProcessConfig)
################################################################################


class TestProcessConfig(ProcessConfig):
    """
    Configuration parameter of the test_process.py program

    ioptimz - sets ioptimz (optimisation solver) in IN.DAT

    epsvmc  - sets epsvmc (VMCON error tolerance) in IN.DAT

    epsfcn  - sets epsfcn (finite diff. steplength) in IN.DAT

    minmax  - sets minmax (figure of merit switch) in IN.DAT

    """

    ioptimz = "None"
    epsvmc = "None"
    epsfcn = "None"
    minmax = "None"

    def __init__(self, filename="test_process.conf"):
        """create configuration instance"""

        self.filename = filename

        super().set_attributes()

        buf = self.get_attribute("ioptimz")
        if buf is not None:
            self.ioptimz = buf

        buf = self.get_attribute("epsvmc")
        if buf is not None:
            self.epsvmc = buf

        buf = self.get_attribute("epsfcn")
        if buf is not None:
            self.epsfcn = buf

        buf = self.get_attribute("minmax")
        if buf is not None:
            self.minmax = buf

    def echo(self):
        """echos the values of the current class"""

        print()
        super().echo()

        if self.ioptimz != "None":
            print(f"ioptimz              {self.ioptimz}")
        if self.epsvmc != "None":
            print(f"epsvmc               {self.epsvmc}")
        if self.epsfcn != "None":
            print(f"epsfcn               {self.epsfcn}")
        if self.minmax != "None":
            print(f"minmax               {self.minmax}")
        print()
        sleep(1)

    def modify_in_dat(self):
        """modifies IN.DAT using the configuration parameters"""

        in_dat = InDat()

        # by convention all variablenames are lower case
        if self.ioptimz != "None":
            in_dat.add_parameter("ioptimz", self.ioptimz)
        if self.epsvmc != "None":
            in_dat.add_parameter("epsvmc", self.epsvmc)

        if self.epsfcn != "None":
            in_dat.add_parameter("epsfcn", self.epsfcn)

        if self.minmax != "None":
            in_dat.add_parameter("minmax", self.minmax)

        in_dat.write_in_dat(output_filename="IN.DAT")


################################################################################
# class RunProcessConfig(ProcessConfig)
################################################################################


class RunProcessConfig(ProcessConfig):
    """
    Configuration parameters of the run_process.py program

    no_allowed_unfeasible - the number of allowed unfeasible points in a sweep

    create_itervar_diff - boolean to indicate the creation of a summary file
                          of the iteration variable values at each stage

    add_ixc - List of iteration variables to be added to IN.DAT

    del_ixc - List of iteration variables to be deleted from IN.DAT

    add_icc - List of constrained equations to be added to IN.DAT

    del_icc - List of constrained equations to be deleted from IN.DAT

    dictvar - Dictionary mapping variable name to new value (replaces old
              or gets appended)

    del_var - List of variables to be deleted from IN.DAT

    """

    no_allowed_unfeasible = 0
    create_itervar_diff = False
    add_ixc: ClassVar = []
    del_ixc: ClassVar = []
    add_icc: ClassVar = []
    del_icc: ClassVar = []
    dictvar: ClassVar = {}
    del_var: ClassVar = []

    def __init__(self, filename="run_process.conf"):
        """
        creates an instance of the RunProcessConfig class
        that stores all configuration parameters of the
        run_process.py
        """
        # TODO filename can be a Path object or a str here
        self.filename = filename

        super().set_attributes()

        buf = self.get_attribute("no_allowed_unfeasible")
        if buf is not None:
            self.no_allowed_unfeasible = int(buf)

        buf = self.get_attribute("create_itervar_diff")
        if buf is not None:
            if buf.lower() in ["true", "y", "yes"]:
                self.create_itervar_diff = True
            elif buf.lower() in ["false", "n", "no"]:
                self.create_itervar_diff = False
            else:
                print(
                    "WARNING: Value for create_itervar_diff\
 is not defined!",
                    file=stderr,
                )

        self.add_ixc = self.get_attribute_csv_list("add_ixc")

        self.del_ixc = self.get_attribute_csv_list("del_ixc")

        self.add_icc = self.get_attribute_csv_list("add_icc")

        self.del_icc = self.get_attribute_csv_list("del_icc")

        self.set_del_var()

        self.set_dictvar()

    def get_attribute_csv_list(self, attributename):
        """
        get class attribute list from configuration file
        expects comma separated values
        """

        if not self.configfileexists:
            return []

        try:
            with open(self.filename) as configfile:
                attribute_list = []

                for line in configfile:
                    condense = line.replace(" ", "")
                    condense = condense.rstrip()
                    lcase = condense.lower()
                    if (
                        (len(condense) > 0)
                        and (condense[0] != "*")
                        and (attributename == lcase[: len(attributename)])
                    ):
                        buf = condense[condense.find("=") + 1 :].split(",")
                        if buf[-1] == "":  # if last value has ended on comma
                            buf = buf[:-1]
                        attribute_list += buf
        except FileNotFoundError:
            print(f"Error: No config file named {self.filename}", file=stderr)
            self.configfileexists = False
            return []
        return attribute_list

    def set_del_var(self):
        """sets the del_var attribute from the config file"""

        if not self.configfileexists:
            return

        try:
            with open(self.filename) as configfile:
                for line in configfile:
                    condense = line.replace(" ", "")
                    condense = condense.rstrip()
                    lcase = condense.lower()
                    if (
                        (len(condense) > 0)
                        and (condense[0] != "*")
                        and (lcase[:8] == "del_var_")
                        and (len(condense) > 8)
                    ):
                        self.del_var += [condense[8:]]
        except FileNotFoundError:
            print(f"Error: No config file named {self.filename}", file=stderr)
            self.configfileexists = False
            return

    def set_dictvar(self):
        """sets the dictvar attribute from config file"""

        if not self.configfileexists:
            return

        try:
            with open(self.filename) as configfile:
                for line in configfile:
                    condense = line.replace(" ", "")
                    condense = condense.rstrip()
                    lcase = condense.lower()
                    if len(condense) > 0 and (condense[0] != "*") and "=" in lcase:
                        varname = lcase[: lcase.find("=")]
                        auxvar = condense[condense.find("=") + 1 :]
                        if varname[:4] == "var_" and auxvar != "":
                            self.dictvar[varname[4:]] = auxvar
        except FileNotFoundError:
            print(f"Error: No config file named {self.filename}", file=stderr)
            self.configfileexists = False
            return

    def echo(self):
        """echos the values of the current class"""

        print()
        super().echo()

        print(f"no. allowed UNFEASIBLE points {self.no_allowed_unfeasible:d}")
        if self.create_itervar_diff:
            print(
                "Set to create a summary file of the iteration variable\
 values!"
            )

        if self.add_ixc != []:
            print("add_ixc", self.add_ixc)
        if self.del_ixc != []:
            print("del_ixc", self.del_ixc)
        if self.add_icc != []:
            print("add_icc", self.add_icc)
        if self.del_icc != []:
            print("del_icc", self.del_icc)
        for key, value in self.dictvar.items():
            print(f"set {key}  to {value}")
        if self.del_var != []:
            print("del_var", self.del_var)

        print()
        sleep(1)

    def modify_in_dat(self):
        """modifies IN.DAT using the configuration parameters"""

        # Need to keep this order!
        # If bounds are modified in vars, but ixc is newly added,
        # bounds do not get put into IN.DAT. Hence, vars needs to be modified
        # first.
        self.modify_ixc()
        self.modify_icc()
        self.modify_vars()

    def modify_vars(self):
        """modifies IN.DAT by adding, deleting and modifiying variables"""

        in_dat = InDat()

        # add and modify variables
        for key in self.dictvar:
            set_variable_in_indat(in_dat, key, self.dictvar[key])

        # delete variables
        for key in self.del_var:
            key = key.lower()
            if "bound" in key:
                number = (key.split("("))[1].split(")")[0]
                if "boundu" in key:
                    in_dat.remove_bound(number, "u")
                else:
                    in_dat.remove_bound(number, "l")
            else:
                in_dat.remove_parameter(key)

        in_dat.write_in_dat(output_filename="IN.DAT")

    def modify_ixc(self):
        """modifies the array of iteration variables in IN.DAT"""

        # check that there is no variable in both lists
        if set(self.add_ixc).intersection(self.del_ixc) != set():
            print(
                "Error: You are trying to add and delete \
                   the same variable from ixc!",
                file=stderr,
            )
            exit()

        in_dat = InDat()

        for iter_var in self.add_ixc:
            in_dat.add_iteration_variable(int(iter_var))

        for iter_var in self.del_ixc:
            in_dat.remove_iteration_variable(int(iter_var))

        in_dat.write_in_dat(output_filename="IN.DAT")

    def modify_icc(self):
        """modifies the array of constraint equations in IN.DAT"""

        # check that there is no variable in both lists
        if set(self.add_icc).intersection(self.del_icc) != set():
            print(
                "Error: You are trying to add and delete the same \
                  variable from icc!",
                file=stderr,
            )
            exit()

        in_dat = InDat()

        for constr in self.add_icc:
            in_dat.add_constraint_equation(int(constr))

        for constr in self.del_icc:
            in_dat.remove_constraint_equation(int(constr))

        in_dat.write_in_dat(output_filename="IN.DAT")


################################################################################
# class UncertaintiesConfig(RunProcessConfig)
################################################################################


class UncertaintiesConfig(ProcessConfig, Config):
    """
    Configuration parameters for evaluate_uncertainties.py program
    """

    no_allowed_unfeasible = 2
    no_scans = 5
    no_samples = 1000
    output_mean = 8056.98
    uncertainties: ClassVar = []
    morris_uncertainties: ClassVar = []
    sobol_uncertainties: ClassVar = []
    output_vars: ClassVar = []
    dict_results: ClassVar = {}
    ncdf_writer = None
    figure_of_merit = "rmajor"
    latin_hypercube_level = 4
    vary_iteration_variables = False

    def __init__(self, configfilename="config_evaluate_uncertainties.json"):
        """
        creates and instance of the UncertaintiesConfig class
        """

        # TODO: Once ProcessConfig has been ported to only use json
        # use the ProcessConfig __init__ routine to set up these
        # parameters
        super().__init__(configfilename)
        self.filename = configfilename

        self.wdir = os.path.abspath(
            self.get("config", "working_directory", default=self.wdir)
        )
        self.or_in_dat = os.path.abspath(
            self.get("config", "IN.DAT_path", default=self.or_in_dat)
        )
        self.niter = self.get("config", "no_iter", default=self.niter)
        self.u_seed = self.get("config", "pseudorandom_seed", default=self.u_seed)
        self.factor = self.get("config", "factor", default=self.factor)
        self.comment = self.get("config", "runtitle", default=self.comment)

        # additional new parameters
        self.no_scans = self.get("no_scans", default=self.no_scans)
        self.no_samples = self.get("no_samples", default=self.no_samples)
        self.uncertainties = self.get("uncertainties", default=self.uncertainties)
        self.morris_uncertainties = self.get(
            "morris_uncertainties", default=self.morris_uncertainties
        )
        self.sobol_uncertainties = self.get(
            "sobol_uncertainties", default=self.sobol_uncertainties
        )
        self.output_vars = self.get("output_vars", default=self.output_vars)
        self.output_mean = self.get("output_mean", default=self.output_mean)
        self.figure_of_merit = self.get("figure_of_merit", default=self.figure_of_merit)
        self.latin_hypercube_level = self.get(
            "latin_hypercube_level", default=self.latin_hypercube_level
        )
        self.vary_iteration_variables = self.get(
            "vary_iteration_variables", default=self.vary_iteration_variables
        )
        # setup the output_vars
        for u_dict in self.uncertainties:
            if u_dict["varname"] not in self.output_vars:
                self.output_vars += [u_dict["varname"]]

        # add normalised constraints/iteration variables to output
        in_dat = InDat(self.or_in_dat)
        nvar = in_dat.number_of_itvars
        for i in range(1, nvar + 1):
            nitvar = f"nitvar{i:03}"
            if nitvar not in self.output_vars:
                self.output_vars += [nitvar]
        neqns = in_dat.number_of_constraints
        for i in range(1, neqns + 1):
            normres = f"normres{i:03}"
            if normres not in self.output_vars:
                self.output_vars += [normres]

        # treat special cases (bounds, fimp, zref)
        add_bounds = False
        add_zref = False
        del_list = []
        for i, varname in enumerate(self.output_vars):
            if "bound" in varname:
                del_list += [varname]
                add_bounds = True
            elif "fimp(" in varname:
                # has different format in MFILE!!
                fimpno = int(varname.split("(")[1].split(")")[0])
                self.output_vars[i] = f"fimp({fimpno:02}"
            elif "zref" in varname:
                del_list += [varname]
                add_zref = True

        for varname in del_list:
            self.output_vars.remove(varname)

        if add_bounds:
            self.output_vars += ["bounds"]
        if add_zref:
            self.output_vars += ["zref"]

        for varname in self.output_vars:
            self.dict_results[varname] = []

    def echo(self):
        """echos the values of the current class"""

        print()
        super().echo()

        print(f"No scans            {self.no_scans:d}")
        print(f"No samples          {self.no_samples:d}")
        if self.uncertainties != []:
            print("uncertainties:")
            for item in self.uncertainties:
                print("     ", item["varname"])
                for key in item:
                    if key not in ["varname"]:
                        print("     ", key, item[key])
                print(" -------")
        if self.output_vars != []:
            print("output vars        ", self.output_vars)
        print()
        sleep(1)

    def modify_in_dat(self):
        """modifies IN.DAT before running uncertainty evaluation"""

        # Load dicts from dicts JSON file
        dicts = get_dicts()

        in_dat = InDat()

        # Is IN.DAT already having a scan?
        isweep = get_from_indat_or_default(in_dat, "isweep")
        if isweep > 0:
            # check that sweep variable is not an iteration variable
            # and all values are the same value
            nsweep = get_from_indat_or_default(in_dat, "nsweep")
            ixc = get_from_indat_or_default(in_dat, "ixc")
            sweep = get_from_indat_or_default(in_dat, "sweep")
            if (
                (str(nsweep) in dicts["DICT_NSWEEP2IXC"])
                and (dicts["DICT_NSWEEP2IXC"][str(nsweep)] in ixc)
                and (not all(sweep[0] == item for item in sweep))
            ):
                if self.no_scans != isweep:
                    # Change no of sweep points to correct value!!
                    set_variable_in_indat(in_dat, "isweep", self.no_scans)
                    value = sweep[0]
                    set_variable_in_indat(in_dat, "sweep", [value] * self.no_scans)
                # Else: we can actually use this scan

            else:
                print(
                    "Error: Inbuild sweep is not compatible with uncertainty\
 evaluation! Edit IN.DAT file!",
                    file=stderr,
                )
                exit()
        else:
            # create a scan!
            nsweep = "3"
            if nsweep in dicts["DICT_NSWEEP2IXC"]:
                # TODO: if this ever happens, program this testing whether
                # a certain variable is used as iteration variable, if not
                # choose another
                print(
                    "Error: The developer should program this more wisely\
 using a sweep variable that is not an iteration variable!",
                    file=stderr,
                )
                exit()
            else:
                set_variable_in_indat(in_dat, "nsweep", nsweep)
                set_variable_in_indat(in_dat, "isweep", self.no_scans)
                value = get_from_indat_or_default(
                    in_dat, dicts["DICT_NSWEEP2VARNAME"][nsweep]
                )
                set_variable_in_indat(in_dat, "sweep", [value] * self.no_scans)

        # write comment in runtitle!
        runtitle = self.comment.replace(",", " ")
        runtitle = runtitle.replace("\n", " ")
        set_variable_in_indat(in_dat, "runtitle", runtitle)

        # set epsvmc to appropriate value!
        # recommendation from solver work!
        set_variable_in_indat(in_dat, "epsvmc", 1e-8)

        in_dat.write_in_dat(output_filename="IN.DAT")

    def checks_before_run(self):
        """run several checks before you start running"""

        # Load dicts from dicts JSON file
        dicts = get_dicts()

        if self.uncertainties == {}:
            print(
                "Error: No uncertain parameter specified in config file!", file=stderr
            )
            exit()

        # check uncertainties are not conflicting
        u_var_list = []
        for u_dict in self.uncertainties:
            u_var_list += [u_dict["varname"]]
        if len(u_var_list) != len(set(u_var_list)):
            print("Error: You have multiply defined uncertain variables!", file=stderr)
            exit()

        if self.output_vars == []:
            print("Error: No output variables specified in config file!", file=stderr)
            exit()

        in_dat = InDat()
        ixc_list = in_dat.data["ixc"].get_value
        assert isinstance(ixc_list, list)

        ixc_varname_list = [dicts["DICT_IXC_SIMPLE"][str(x)] for x in ixc_list]

        for u_dict in self.uncertainties:
            varname = u_dict["varname"].lower()
            if varname in ixc_varname_list:
                print(
                    "Error: an uncertain variable should never be an\
 iteration variable at the same time!",
                    varname,
                    file=stderr,
                )
                exit()

            if u_dict["errortype"].lower() == "uniform":
                lbound = u_dict["lowerbound"]
                ubound = u_dict["upperbound"]
                if lbound > ubound:
                    print(
                        "Error: the lower bound of the uncertain variable",
                        varname,
                        "is higher than its upper bound!",
                        file=stderr,
                    )
                    exit()

                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    # check bounds are inside input bounds
                    if lbound < dicts["DICT_INPUT_BOUNDS"][varname]["lb"]:
                        u_dict["lowerbound"] = dicts["DICT_INPUT_BOUNDS"][varname]["lb"]
                        print(
                            "Warning: The lower bound of the uncertain variable",
                            varname,
                            "is lower than the lowest allowed input",
                            "value! \n Corrected value to",
                            u_dict["lowerbound"],
                            file=stderr,
                        )
                    if ubound > dicts["DICT_INPUT_BOUNDS"][varname]["ub"]:
                        u_dict["upperbound"] = dicts["DICT_INPUT_BOUNDS"][varname]["ub"]
                        print(
                            "Warning: The upper bound of the uncertain variable",
                            varname,
                            "is higher than the highest allowed input",
                            "value! \n Corrected value to",
                            u_dict["upperbound"],
                            file=stderr,
                        )

            elif u_dict["errortype"].lower() == "relative":
                err = u_dict["percentage"] / 100.0
                lbound = u_dict["mean"] * (1.0 - err)
                ubound = u_dict["mean"] * (1.0 + err)

                if err < 0.0:
                    print(
                        "Error: The percentage of the uncertain variable",
                        varname,
                        "should never be negative!",
                        file=stderr,
                    )
                    exit()

                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    # check mean is inside input bounds
                    if u_dict["mean"] < dicts["DICT_INPUT_BOUNDS"][varname]["lb"]:
                        print(
                            "Error: The mean of the uncertain variable",
                            varname,
                            "is lower than the lowest allowed input!",
                            dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                            file=stderr,
                        )
                        exit()
                    elif u_dict["mean"] > dicts["DICT_INPUT_BOUNDS"][varname]["ub"]:
                        print(
                            "Error: The mean of the uncertain variable",
                            varname,
                            "is higher than the highest allowed input!",
                            dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                            file=stderr,
                        )
                        exit()

                    # check bounds are inside input bounds
                    if lbound < dicts["DICT_INPUT_BOUNDS"][varname]["lb"]:
                        u_dict["errortype"] = "uniform"
                        u_dict["lowerbound"] = dicts["DICT_INPUT_BOUNDS"][varname]["lb"]
                        u_dict["upperbound"] = ubound
                        print(
                            "Warning: The lower bound of the uncertain variable",
                            varname,
                            "is lower than the lowest allowed input",
                            "value! \n Corrected value to",
                            u_dict["lowerbound"],
                            file=stderr,
                        )
                    if ubound > dicts["DICT_INPUT_BOUNDS"][varname]["ub"]:
                        if u_dict["errortype"] != "uniform":
                            u_dict["errortype"] = "uniform"
                            u_dict["lowerbound"] = lbound
                        u_dict["upperbound"] = dicts["DICT_INPUT_BOUNDS"][varname]["ub"]
                        print(
                            "Warning: The upper bound of the uncertain variable",
                            varname,
                            "is higher than the highest allowed input",
                            "value! \n Corrected value to",
                            u_dict["upperbound"],
                            file=stderr,
                        )

            elif u_dict["errortype"].lower() in [
                "gaussian",
                "lowerhalfgaussian",
                "upperhalfgaussian",
            ]:
                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    # check mean is inside input bounds
                    if u_dict["mean"] < dicts["DICT_INPUT_BOUNDS"][varname]["lb"]:
                        print(
                            "Error: The mean of the uncertain variable",
                            varname,
                            "is lower than the lowest allowed input!",
                            u_dict["mean"],
                            "<",
                            dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                            file=stderr,
                        )
                        exit()
                    elif u_dict["mean"] > dicts["DICT_INPUT_BOUNDS"][varname]["ub"]:
                        print(
                            "Error: The mean of the uncertain variable",
                            varname,
                            "is higher than the highest allowed input!",
                            u_dict["mean"],
                            ">",
                            dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                            file=stderr,
                        )
                        exit()

    def set_sample_values(self):
        """determines the values of each sample point and orders them"""

        # Load dicts from dicts JSON file
        dicts = get_dicts()

        for u_dict in self.uncertainties:
            varname = u_dict["varname"].lower()
            if u_dict["errortype"].lower() == "gaussian":
                mean = u_dict["mean"]
                std = u_dict["std"]
                values = self.generator.normal(mean, std, self.no_samples)
                # assures values are inside input bounds!
                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    args = argwhere(
                        logical_or(
                            values < dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                            values > dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                        )
                    )
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(
                            logical_or(
                                values < dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                                values > dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                            )
                        )
                else:  # cutoff at 0 - typically negative values are meaningless
                    args = argwhere(values < 0.0)
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(values < 0)

            elif u_dict["errortype"].lower() == "uniform":
                lbound = u_dict["lowerbound"]
                ubound = u_dict["upperbound"]
                values = self.generator.uniform(lbound, ubound, self.no_samples)
            elif u_dict["errortype"].lower() == "relative":
                err = u_dict["percentage"] / 100.0
                lbound = u_dict["mean"] * (1.0 - err)
                ubound = u_dict["mean"] * (1.0 + err)
                values = self.generator.uniform(lbound, ubound, self.no_samples)
            elif u_dict["errortype"].lower() == "lowerhalfgaussian":
                mean = u_dict["mean"]
                std = u_dict["std"]
                values = self.generator.normal(mean, std, self.no_samples)
                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    args = argwhere(
                        logical_or(
                            values < dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                            values > mean,
                        )
                    )
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(
                            logical_or(
                                values < dicts["DICT_INPUT_BOUNDS"][varname]["lb"],
                                values > mean,
                            )
                        )
                else:
                    args = argwhere(logical_or(values < 0.0, values > mean))
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(logical_or(values < 0.0, values > mean))
            elif u_dict["errortype"].lower() == "upperhalfgaussian":
                mean = u_dict["mean"]
                std = u_dict["std"]
                values = self.generator.normal(mean, std, self.no_samples)
                if varname in dicts["DICT_INPUT_BOUNDS"]:
                    args = argwhere(
                        logical_or(
                            values < mean,
                            values > dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                        )
                    )
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(
                            logical_or(
                                values < mean,
                                values > dicts["DICT_INPUT_BOUNDS"][varname]["ub"],
                            )
                        )
                else:
                    args = argwhere(values < mean)
                    while len(args) > 0:
                        values[args] = self.generator.normal(mean, std, args.shape)
                        args = argwhere(values < mean)

            u_dict["samples"] = values

        # order by one parameter
        # TODO: Find a more cunning way to determine which is a good
        # variable to sort by! e.g. largest range?
        # always try to choose a uniform case?
        arr = self.uncertainties[0]["samples"]
        sorted_index = argsort(arr)
        for u_dict in self.uncertainties:
            u_dict["samples"] = u_dict["samples"][sorted_index]

    def go2newsamplepoint(self, sample_index):
        """create a new sample point from uncertainty distributions"""

        in_dat = InDat()

        for u_dict in self.uncertainties:
            value = u_dict["samples"][sample_index]
            varname = u_dict["varname"]
            set_variable_in_indat(in_dat, varname, value)

        in_dat.write_in_dat(output_filename="IN.DAT")

    def write_error_summary(self, sample_index):
        """reads current MFILE and IN.DAT and adds specified output variables
        to error summary file"""

        m_file = MFile(filename="MFILE.DAT")
        nvar = int(m_file.data["nvar"].get_scan(-1))

        if sample_index == 0:
            header = "#sample_index"

            # Uncertain input variables
            for u_dict in self.uncertainties:
                header += " {u_dict['varname']:10s}"

            # normalised iteration varialbes
            for i in range(1, nvar + 1):
                label = m_file.data[f"nitvar{i:03}"].var_description
                header += f" n_{label.replace('_(range_normalised)', ''):8s}"

            # error status, id and ifail
            header += " error_status error_id ifail\n"
            with open(self.wdir + "/UQ_error_summary.txt", "w") as err_summary:
                err_summary.write(header)

        # Uncertain input variables
        output = f"{sample_index:12d}"
        for u_dict in self.uncertainties:
            output += f" {u_dict['samples'][sample_index]:10f}"

        # normalised iteration variables
        for i in range(1, nvar + 1):
            output += f" {m_file.data[f'nitvar{i:03}'].get_scan(-1):10f}"

        # error status and id
        output += (
            f" {int(m_file.data['error_status'].get_scan(-1)):13d} "
            f"{int(m_file.data['error_id'].get_scan(-1)):8d}"
        )
        # ifail
        if m_file.data["error_status"].get_scan(-1) < 3:
            output += f" {int(m_file.data['ifail'].get_scan(-1)):5d}\n"
        else:
            output += "   -1\n"

        with open(self.wdir + "/UQ_error_summary.txt", "a+") as err_summary:
            err_summary.write(output)
