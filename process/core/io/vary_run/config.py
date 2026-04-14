"""
Interfaces for Configuration values for programs
- run_process.py
- test_process.py
"""

import logging
import os
import sys
from contextlib import suppress
from dataclasses import dataclass, field
from pathlib import Path
from shutil import SameFileError, copy

import click
from numpy.random import default_rng

from process.core.io.in_dat import InDat
from process.core.io.mfile import MFile
from process.core.io.vary_run.tools import (
    check_in_dat,
    check_input_error,
    get_neqns_itervars,
    get_variable_range,
    no_unfeasible_mfile,
    process_stopped,
    process_warnings,
    set_variable_in_indat,
)

logger = logging.getLogger(__name__)


@dataclass
class ProcessConfig:
    """Configuration parameters for PROCESS runs"""

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
    solver: str = "vmcon"
    """Solver to use"""
    comment: str = " "
    """additional comment to be written into README.txt"""
    _filename: Path | None = None
    """filename if supplied"""

    @classmethod
    def from_file(cls, filename: str | Path, solver: str = "vmcon"):
        if isinstance(filename, str):
            filename = Path(filename)

        if not filename.is_file():
            raise FileNotFoundError(f"Config file '{filename}' not found")

        wdir = (
            None if (buf := cls.get_attribute(filename, "wdir")) is None else Path(buf)
        )

        if wdir in {None, Path()}:
            wdir = filename.parent

        or_in_dat = (
            None
            if (buf := cls.get_attribute(filename, "ORIGINAL_IN_DAT")) is None
            else Path(buf)
        )
        or_in_dat = (
            cls.get_attribute(filename, "original_in_dat")
            if or_in_dat is None
            else or_in_dat
        )
        or_in_dat = Path("ref_IN.DAT") if or_in_dat is None else Path(or_in_dat)
        if not or_in_dat.is_file():
            or_in_dat = wdir / or_in_dat

        niter = 10 if (buf := cls.get_attribute(filename, "niter")) is None else int(buf)
        u_seed = (
            None
            if (buf := cls.get_attribute(filename, "seed")) in {None, "None"}
            else int(buf)
        )
        factor = (
            1.5 if (buf := cls.get_attribute(filename, "factor")) is None else float(buf)
        )

        comment = cls.get_comment(filename)
        if not comment:
            print(f"No comment in config file {filename}")

        return cls(
            wdir, or_in_dat, niter, u_seed, factor, solver, comment or " ", filename
        )

    @staticmethod
    def get_attribute(filename: Path, attributename: str):
        """Gets class attribute from configuration file"""
        with open(filename) as configfile:
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

    @staticmethod
    def get_comment(filename):
        """Gets the comment line from the configuration file"""
        with open(filename) as configfile:
            for line in configfile:
                condense = line.replace(" ", "")
                condense = condense.rstrip()
                lcase = condense.lower()
                if len(condense) > 0 and (condense[0] != "*") and lcase[:7] == "comment":
                    return line[line.find("=") + 1 :]
        return ""

    def __post_init__(self):
        self._current_iteration = 0
        self._base_input = "{}_IN.DAT"
        self._base_output = "{}_MFILE.DAT"

        if not isinstance(self.or_in_dat, Path) or not self.or_in_dat.is_file():
            raise FileNotFoundError(
                f"Error: {self.or_in_dat} does not exist! Create file or modify config file!",
            )

    def __iter__(self):
        return self

    def __next__(self):
        _neqns, itervars = get_neqns_itervars(wdir=self.wdir)
        lbs, ubs = get_variable_range(itervars, self.factor, self.wdir)
        indat, mfile = self.infile, self.outfile
        self._current_iteration += 1
        if self._current_iteration >= self.niter:
            self.error_status2readme(mfile=mfile)
            raise StopIteration
        self.run_process(self.wdir / indat, self.solver)
        check_input_error(wdir=self.wdir, mfile=mfile)

        return indat, mfile, itervars, lbs, ubs

    @property
    def infile(self):
        return self._base_input.format(self._current_iteration)

    @property
    def outfile(self):
        return self._base_output.format(self._current_iteration)

    def echo(self):
        """Echos the attributes of the class"""
        if self.wdir != Path.cwd():
            print(f"Working directory:   {self.wdir}")
        print(f"Original IN.DAT:     {self.or_in_dat}")
        print(f"Number of iterations {self.niter}")

        if self.u_seed is not None:
            print(f"random seed          {self.u_seed}")
        print(f"variable range factor {self.factor}")
        if self._filename is not None:
            print(f"Config file          {self._filename}")
        if self.comment != "":
            print(f"Comment  {self.comment}")

    def prepare_wdir(self):
        """Prepares the working directory"""
        if not self.wdir.is_dir():
            self.wdir.mkdir(exist_ok=True, parents=True)

        copy(self.or_in_dat, self.wdir / "0_IN.DAT")

        if self._filename is not None:
            with suppress(SameFileError):
                copy(self._filename, self.wdir)

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
        """Creates README.txt containing comment"""
        if self.comment != "":
            Path(self.wdir, "README.txt").write_text(self.comment)

    def error_status2readme(self, mfile="MFILE.DAT"):
        """Appends PROCESS outcome to README.txt"""
        m_file = MFile(filename=self.wdir / mfile)

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
        """Modifies the original IN.DAT file"""

    def setup(self):
        """Sets up the program for running"""
        self.echo()

        self.prepare_wdir()

        self.create_readme()

        self.modify_in_dat()

        check_in_dat("0_IN.DAT")

        self.generator = default_rng(seed=self.u_seed)

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
        from process.main import process_cli  # noqa:PLC0415

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
                logger.exception(
                    "There was a problem with the PROCESS execution!",
                )
        except (KeyboardInterrupt, click.exceptions.Abort):
            raise KeyboardInterrupt from None
        except Exception:
            logger.exception(
                "There was a problem with the PROCESS execution!",
            )

        print("finished.")


@dataclass
class RunProcessConfig(ProcessConfig):
    """Configuration parameters of the run_process.py program"""

    no_allowed_unfeasible: int = 0
    """the number of allowed unfeasible points in a sweep"""
    create_itervar_diff: bool = False
    """boolean to indicate the creation of a summary file of the iteration variable values at each stage"""
    dictvar: dict[str, str] = field(default_factory=dict)
    """Dictionary mapping variable name to new value"""
    del_var: list[str] = field(default_factory=list)
    """List of variables to be deleted from IN.DAT"""
    add_ixc: list[str] = field(default_factory=list)
    """ List of iteration variables to be added to IN.DAT"""
    del_ixc: list[str] = field(default_factory=list)
    """List of iteration variables to be deleted from IN.DAT"""
    add_icc: list[str] = field(default_factory=list)
    """List of constrained equations to be added to IN.DAT"""
    del_icc: list[str] = field(default_factory=list)
    """List of constrained equations to be deleted from IN.DAT"""

    @classmethod
    def from_file(cls, filename: str | Path = "run_process.conf", solver: str = "vmcon"):
        self = super().from_file(filename, solver)

        no_allowed_unfeasible = (
            0
            if (buf := cls.get_attribute(self._filename, "no_allowed_unfeasible"))
            is None
            else int(buf)
        )

        buf = cls.get_attribute(self._filename, "create_itervar_diff")

        create_itervar_diff = False
        if buf is not None:
            if buf.lower() in {"true", "y", "yes"}:
                create_itervar_diff = True
            elif buf.lower() in {"false", "n", "no"}:
                create_itervar_diff = False
            else:
                logger.warning("Value for create_itervar_diff is not defined!")

        add_ixc = cls.get_attribute_csv_list(self._filename, "add_ixc")
        del_ixc = cls.get_attribute_csv_list(self._filename, "del_ixc")
        add_icc = cls.get_attribute_csv_list(self._filename, "add_icc")
        del_icc = cls.get_attribute_csv_list(self._filename, "del_icc")

        return cls(
            self.wdir,
            self.or_in_dat,
            self.niter,
            self.u_seed,
            self.factor,
            self.solver,
            self.comment,
            self._filename,
            no_allowed_unfeasible,
            create_itervar_diff,
            cls.get_dictvar(self._filename),
            cls.get_del_var(self._filename),
            add_ixc,
            del_ixc,
            add_icc,
            del_icc,
        )

    @staticmethod
    def get_attribute_csv_list(filename, attributename):
        """Get class attribute list from configuration file
        expects comma separated values
        """
        attribute_list = []

        with open(filename) as configfile:
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

    @staticmethod
    def get_del_var(filename):
        """Sets the del_var attribute from the config file"""
        del_var = []
        with open(filename) as configfile:
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
                    del_var += [condense[8:]]
        return del_var

    @staticmethod
    def get_dictvar(filename):
        """Sets the dictvar attribute from config file"""
        dictvar = {}
        with open(filename) as configfile:
            for line in configfile:
                condense = line.replace(" ", "")
                condense = condense.rstrip()
                lcase = condense.lower()
                if len(condense) > 0 and (condense[0] != "*") and "=" in lcase:
                    varname = lcase[: lcase.find("=")]
                    auxvar = condense[condense.find("=") + 1 :]
                    if varname[:4] == "var_" and auxvar != "":
                        dictvar[varname[4:]] = auxvar
        return dictvar

    def __next__(self):
        indat, mfile, itervars, lbs, ubs = super().__next__()

        if not process_stopped(wdir=self.wdir, mfile=mfile):
            no_unfeasible = no_unfeasible_mfile(self.wdir, mfile)
            if no_unfeasible <= self.no_allowed_unfeasible:
                if no_unfeasible > 0:
                    logger.warning(
                        "Non feasible point(s) in sweep, But finished anyway! %s ",
                        no_unfeasible,
                    )
                if process_warnings(self.wdir, mfile):
                    print(
                        "\nThere were warnings in the final PROCESS run. "
                        "Please check the log file!\n"
                    )
                    # This means success: feasible solution found
                    raise StopIteration
                logger.warning(
                    "%s non-feasible point(s) in sweep! Rerunning!", no_unfeasible
                )
        else:
            print("PROCESS has stopped without finishing!")

        return indat, mfile, itervars, lbs, ubs

    def echo(self):
        """Echos the values of the current class"""
        super().echo()

        print(f"no. allowed UNFEASIBLE points {self.no_allowed_unfeasible:d}")
        if self.create_itervar_diff:
            print("Set to create a summary file of the iteration variable values!")

        for n, v in (
            ("add_ixc", self.add_ixc),
            ("del_ixc", self.del_ixc),
            ("add_icc", self.add_icc),
            ("del_icc", self.del_icc),
        ):
            if v != []:
                print(n, v)

        for key, value in self.dictvar.items():
            print(f"set {key}  to {value}")
        if self.del_var != []:
            print("del_var", self.del_var)

    def modify_in_dat(self):
        """Modifies IN.DAT using the configuration parameters"""
        # Need to keep this order!
        # If bounds are modified in vars, but ixc is newly added,
        # bounds do not get put into IN.DAT. Hence, vars needs to be modified
        # first.
        self.modify_ixc()
        self.modify_icc()
        self.modify_vars()

    def modify_vars(self):
        """Modifies IN.DAT by adding, deleting and modifiying variables"""
        in_dat = InDat(filename="0_IN.DAT")

        # add and modify variables
        for key in self.dictvar:
            set_variable_in_indat(in_dat, key, self.dictvar[key])

        # delete variables
        for key_ in self.del_var:
            key = key_.lower()
            if "bound" in key:
                number = (key.split("("))[1].split(")")[0]
                if "boundu" in key:
                    in_dat.remove_bound(number, "u")
                else:
                    in_dat.remove_bound(number, "l")
            else:
                in_dat.remove_parameter(key)

        in_dat.write_in_dat(output_filename=self.wdir / "0_IN.DAT")

    def modify_ixc(self):
        """Modifies the array of iteration variables in IN.DAT"""
        # check that there is no variable in both lists
        if set(self.add_ixc).intersection(self.del_ixc) != set():
            logger.error(
                "You are trying to add and delete the same variable from ixc!",
                file=sys.stderr,
            )
            sys.exit()

        in_dat = InDat(filename="0_IN.DAT")

        for iter_var in self.add_ixc:
            in_dat.add_iteration_variable(int(iter_var))

        for iter_var in self.del_ixc:
            in_dat.remove_iteration_variable(int(iter_var))

        in_dat.write_in_dat(output_filename=self.wdir / "0_IN.DAT")

    def modify_icc(self):
        """Modifies the array of constraint equations in IN.DAT"""
        # check that there is no variable in both lists
        if set(self.add_icc).intersection(self.del_icc) != set():
            logger.error(
                "You are trying to add and delete the same variable from icc!",
                file=sys.stderr,
            )
            sys.exit()

        in_dat = InDat(filename="0_IN.DAT")

        for constr in self.add_icc:
            in_dat.add_constraint_equation(int(constr))

        for constr in self.del_icc:
            in_dat.remove_constraint_equation(int(constr))

        in_dat.write_in_dat(output_filename=self.wdir / "0_IN.DAT")
