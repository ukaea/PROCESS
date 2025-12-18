"""
Author: Hanni Lux (Hanni.Lux@ccfe.ac.uk)

Interfaces for Configuration values for programs
- run_process.py
- test_process.py

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

from numpy.random import default_rng

from process.io.in_dat import InDat
from process.io.mfile import MFile
from process.io.process_funcs import (
    check_in_dat,
    set_variable_in_indat,
)

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
            Path(directory + "/README.txt").write_text(self.comment)

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
                Path(directory + "/README.txt").write_text(error_status)

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
