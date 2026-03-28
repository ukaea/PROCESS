"""
Interfaces for Configuration values for programs
- run_process.py
- test_process.py
"""

import logging
import os
from contextlib import suppress
from dataclasses import dataclass, field
from pathlib import Path
from shutil import SameFileError, copy
from sys import stderr
from time import sleep

import click
from numpy.random import default_rng

from process.core.io.in_dat import InDat
from process.core.io.mfile import MFile
from process.core.io.process_funcs import (
    check_in_dat,
    set_variable_in_indat,
)

logger = logging.getLogger(__name__)


@dataclass
class ProcessConfig:
    """Configuration parameters for PROCESS runs"""

    filename: str | Path = ""
    """Configuration file name"""
    wdir: Path = field(default_factory=Path.cwd)
    """Working directory"""
    or_in_dat: Path = Path("IN.DAT")
    """Original IN.DAT file"""
    niter: int = 10
    """(Maximum) number of iterations"""
    u_seed: int | None = None
    """User specified seed value for the random number generator"""
    factor: float = 1.5
    """Multiplication factor adjusting the range in which the original
iteration variables should get varied"""
    comment: str = " "
    """additional comment to be written into README.txt"""

    def __post_init__(self):
        if isinstance(self.filename, str) and not self.filename:
            # run_process is not being run by the test suite
            # Get working dir and input files from run_process.conf
            self.wdir = (
                Path(self.wdir)
                if (buf := self.get_attribute("wdir")) is None
                else Path(buf)
            )

            buf = self.get_attribute("ORIGINAL_IN_DAT")
            if buf is not None:
                self.or_in_dat = Path(buf)
        else:
            self.filename = Path(self.filename)
            self.wdir = self.filename.parent
            # TODO It appears that the input file for run_process shouldn't be
            # called IN.DAT. Calling it ref_IN.DAT for now
            original_in_dat = self.get_attribute("original_in_dat") or "ref_IN.DAT"
            self.or_in_dat = self.filename.parent / original_in_dat

        if not isinstance(self.or_in_dat, Path) or not self.or_in_dat.is_file():
            raise FileNotFoundError(
                f"Error: {self.or_in_dat} does not exist! Create file or modify config file!",
            )

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

    @property
    def configfileexists(self):
        if self.filename is None:
            return False
        return Path(self.filename).is_file()

    def echo(self):
        """echos the attributes of the class"""

        if self.wdir != Path.cwd():
            print(f"Working directory:   {self.wdir}")
        print(f"Original IN.DAT:     {self.or_in_dat}")
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

        if not self.wdir.is_dir():
            self.wdir.mkdir(exist_ok=True, parents=True)

        if self.or_in_dat != Path("IN.DAT") or self.wdir != Path.cwd():
            copy(self.or_in_dat, self.wdir / "IN.DAT")
        else:
            copy(self.or_in_dat, "Input_IN.DAT")

        if self.configfileexists:
            with suppress(SameFileError):
                copy(self.filename, self.wdir)

        for file in (
            "OUT.DAT",
            "MFILE.DAT",
            "README.txt",
            "SolverTest.out",
            "process.log",
            "uncertainties.nc",
            "time.info",
        ):
            Path(self.wdir, file).unlink(missing_ok=True)
        for f in self.wdir.glob("*.pdf"):
            f.unlink()

        os.chdir(self.wdir)

    def create_readme(self):
        """creates README.txt containing comment"""

        if self.comment != "":
            Path(self.wdir, "README.txt").write_text(self.comment)

    def error_status2readme(self):
        """appends PROCESS outcome to README.txt"""

        m_file = MFile(filename=self.wdir / "MFILE.DAT")

        error_status = (
            f"Error status: {m_file.data['error_status'].get_scan(-1)}  "
            f"Error ID: {m_file.data['error_id'].get_scan(-1)}\n"
        )

        if self.comment != "":
            with open(Path(self.wdir, "README.txt"), "a") as readme:
                readme.write(error_status)
        else:
            Path(self.wdir, "README.txt").write_text(error_status)

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

        with open(self.filename) as configfile:
            for line in configfile:
                condense = line.replace(" ", "")
                condense = condense.rstrip()
                lcase = condense.lower()
                if len(condense) > 0 and (condense[0] != "*") and lcase[:7] == "comment":
                    self.comment = line[line.find("=") + 1 :]
                    return True

            return False

    def get_attribute(self, attributename):
        """gets class attribute from configuration file"""

        if not self.configfileexists:
            return None

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

        return None

    def run_process(self, input_path: Path, solver: str = "vmcon"):
        """Perform a single run of PROCESS, catching any errors.

        Parameters
        ----------
        input_path :
            the input file to run on
        solver :
            which solver to use, as specified in solver.py, defaults to "vmcon"
        """
        # TODO should call SingleRun directly...at least this is not a subprocess!
        from process.main import process_cli

        print("PROCESS run started ...", end="")

        try:
            # Run process on an IN.DAT file
            process_cli.main(
                args=["-i", str(input_path), "--solver", solver], standalone_mode=False
            )
        except SystemExit as err:
            if err.code != 0:
                # Process has exited with a non-zero exit code.
                # Catch this exception to allow execution to continue without exiting
                logger.error(
                    f"There was a problem with the PROCESS execution! {err}",
                )
        except (KeyboardInterrupt, click.exceptions.Abort):
            raise KeyboardInterrupt from None
        except Exception as e:
            print(e)
            logger.error(
                f"There was a problem with the PROCESS execution! {e}",
            )

        print("finished.")


@dataclass
class RunProcessConfig(ProcessConfig):
    """Configuration parameters of the run_process.py program"""

    filename: str | Path | None = "run_process.conf"
    dictvar: dict[str, str] = field(default_factory=dict, init=False)
    """Dictionary mapping variable name to new value"""
    del_var: list[str] = field(default_factory=list, init=False)
    """List of variables to be deleted from IN.DAT"""
    no_allowed_unfeasible: int = 0
    """the number of allowed unfeasible points in a sweep"""
    create_itervar_diff: bool = False
    """boolean to indicate the creation of a summary file of the iteration variable values at each stage"""

    def __post_init__(self):
        super().__post_init__()

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
                    "WARNING: Value for create_itervar_diff is not defined!",
                    file=stderr,
                )

        self.add_ixc = self.get_attribute_csv_list("add_ixc")
        """ List of iteration variables to be added to IN.DAT"""

        self.del_ixc = self.get_attribute_csv_list("del_ixc")
        """List of iteration variables to be deleted from IN.DAT"""

        self.add_icc = self.get_attribute_csv_list("add_icc")
        """List of constrained equations to be added to IN.DAT"""

        self.del_icc = self.get_attribute_csv_list("del_icc")
        """List of constrained equations to be deleted from IN.DAT"""

        self.set_del_var()

        self.set_dictvar()

    def get_attribute_csv_list(self, attributename):
        """get class attribute list from configuration file
        expects comma separated values
        """
        attribute_list = []

        if not self.configfileexists:
            return attribute_list

        with open(self.filename) as configfile:
            for line in configfile:
                condense = line.replace(" ", "").rstrip()
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
        return attribute_list

    def set_del_var(self):
        """sets the del_var attribute from the config file"""

        if not self.configfileexists:
            return

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

    def set_dictvar(self):
        """sets the dictvar attribute from config file"""

        if not self.configfileexists:
            return

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

    def echo(self):
        """echos the values of the current class"""

        print()
        super().echo()

        print(f"no. allowed UNFEASIBLE points {self.no_allowed_unfeasible:d}")
        if self.create_itervar_diff:
            print("Set to create a summary file of the iteration variable values!")

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

        in_dat.write_in_dat(output_filename=self.wdir / "IN.DAT")

    def modify_ixc(self):
        """modifies the array of iteration variables in IN.DAT"""

        # check that there is no variable in both lists
        if set(self.add_ixc).intersection(self.del_ixc) != set():
            print(
                "Error: You are trying to add and delete the same variable from ixc!",
                file=stderr,
            )
            exit()

        in_dat = InDat()

        for iter_var in self.add_ixc:
            in_dat.add_iteration_variable(int(iter_var))

        for iter_var in self.del_ixc:
            in_dat.remove_iteration_variable(int(iter_var))

        in_dat.write_in_dat(output_filename=self.wdir / "IN.DAT")

    def modify_icc(self):
        """modifies the array of constraint equations in IN.DAT"""

        # check that there is no variable in both lists
        if set(self.add_icc).intersection(self.del_icc) != set():
            print(
                "Error: You are trying to add and delete the same variable from icc!",
                file=stderr,
            )
            exit()

        in_dat = InDat()

        for constr in self.add_icc:
            in_dat.add_constraint_equation(int(constr))

        for constr in self.del_icc:
            in_dat.remove_constraint_equation(int(constr))

        in_dat.write_in_dat(output_filename=self.wdir / "IN.DAT")
