import logging
from contextlib import contextmanager

from process.core import constants, process_output
from process.core.solver.evaluators import Evaluators
from process.core.solver.iteration_variables import (
    load_iteration_variables,
    load_scaled_bounds,
)
from process.core.solver.solver import get_solver
from process.data_structure.numerics import FiguresOfMerit, PROCESSRunMode

logger = logging.getLogger(__name__)


class SolverHandler:
    """Creates and runs a solver instance.

    This may be an optimiser (e.g. VMCON) or an equation solver (e.g. fsolve).

    Parameters
    ----------
    models : process.main.Models
        physics and engineering model objects
    solver_name : str
        which solver to use, as specified in solver.py
    data: DataStructure
        data structure object for providing objective/constraint data to
        the solver
    """

    def __init__(self, models, solver_name, data):
        self.models = models
        self.solver_name = solver_name
        self.data = data

    def run(self):
        """Run solver and retry if it fails in certain ways."""
        # Initialise iteration variables and bounds in Fortran
        load_iteration_variables(self.data)
        load_scaled_bounds(self.data)

        # Initialise iteration variables and bounds in Python: relies on Fortran
        # iteration variables being defined above
        # Trim maximum size arrays down to actually used size
        x = self.data.numerics.xcm[: self.data.numerics.nvar]
        bndl = self.data.numerics.itv_scaled_lower_bounds[: self.data.numerics.nvar]
        bndu = self.data.numerics.itv_scaled_upper_bounds[: self.data.numerics.nvar]

        # Evaluators() calculates the objective and constraint functions and
        # their gradients for a given vector x
        evaluators = Evaluators(self.models, self.data, x)

        # Configure solver for problem
        self.solver = get_solver(self.data, self.solver_name)
        self.solver.set_evaluators(evaluators)
        self.solver.set_bounds(bndl, bndu)
        self.solver.set_opt_params(x)
        # Define total number of constraints and equality constraints
        self.solver.set_constraints(
            m=self.data.numerics.neqns + self.data.numerics.nineqns,
            meq=self.data.numerics.neqns,
        )
        ifail = self.solver.solve()

        # If VMCON optimisation has failed then try altering value of epsfcn
        if self.solver_name == "vmcon":
            if ifail != 1:
                with epsfcn_context(self.data.numerics):
                    self.solver.solve()

            # If VMCON has exited with error code 5 try another run using a multiple
            # of the identity matrix as input for the Hessian b(n,n)
            # Only do this if VMCON has not iterated (nviter=1)
            if ifail == 5 and self.data.numerics.nviter < 2:
                print(
                    "VMCON error code = 5.  Rerunning VMCON with a new initial "
                    "estimate of the second derivative matrix."
                )
                self.solver.set_b(2.0)
                ifail = self.solver.solve()

        self.output()

        return ifail

    def output(self):
        """Store results back in Fortran self.data.numerics module.

        Objective function value, solution vector and constraints vector.
        """
        self.data.numerics.norm_objf = self.solver.objf
        # Slicing required due to Fortran arrays being maximum possible, rather
        # than required, size
        self.data.numerics.xcm[: self.solver.x.shape[0]] = self.solver.x
        self.data.numerics.rcm[: self.solver.conf.shape[0]] = self.solver.conf

        nums = self.data.numerics

        process_output.oheadr(constants.NOUT, "Numerics")
        process_output.ocmmnt(
            constants.NOUT,
            f"PROCESS has performed a {'fsolve' if self.solver == 'fsolve' else 'VMCON'} (optimisation) run.",
        )
        ifail = self.solver.info
        if ifail != 1:
            process_output.ovarin(constants.NOUT, "Error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.IOTTY, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            print()

            logger.critical("Solver returns with ifail /= 1. %s", ifail)

            if self.solver_name == "vmcon":
                self.solver.verror()

            process_output.oblnkl(constants.NOUT)
            print()
        else:
            # Solution found
            descr = "consistent" if self.solver == "fsolve" else "feasible"
            process_output.ocmmnt(
                constants.NOUT, f"and found a {descr} set of parameters."
            )
            process_output.oheadr(constants.IOTTY, f"PROCESS found a {descr} solution")
            process_output.oblnkl(constants.NOUT)
            process_output.ovarin(constants.NOUT, "Error flag", "(ifail)", ifail)

            if nums.sqsumsq >= 1.0e-2:
                string = (
                    "WARNING: Constraint residues are HIGH; consider re-running\n"
                    "   with lower values of EPSVMC to confirm convergence...\n"
                    "   (should be able to get down to about 1.0E-8 okay)\n"
                )
                process_output.ocmmnt(constants.NOUT, ("\n" + string))
                print(string)

                logger.warning(f"High final constraint residues. {nums.sqsumsq=}")

        for d, var, v in (
            ("Number of iteration variables", "(nvar)", nums.nvar),
            (
                "Number of constraints (total)",
                "(neqns+nineqns)",
                nums.neqns + nums.nineqns,
            ),
            ("Optimisation switch", "(ioptimz)", nums.ioptimz),
        ):
            process_output.ovarin(constants.NOUT, d, var, v)

        process_output.ocmmnt(
            constants.NOUT,
            f"     {PROCESSRunMode(nums.ioptimz).description}",
        )

        # Objective function output: none for fsolve
        if self.solver_name != "fsolve":
            process_output.ovarin(
                constants.NOUT,
                "Figure of merit switch",
                "(minmax)",
                nums.minmax,
            )

            nums.objf_name = f'"{FiguresOfMerit(abs(nums.minmax)).description}"'

            for d, var, v, o in (
                ("Objective function name", "(objf_name)", nums.objf_name, ""),
                ("Normalised objective function", "(norm_objf)", nums.norm_objf, "OP "),
                (
                    "VMCON convergence parameter",
                    "(convergence_parameter)",
                    self.data.globals.convergence_parameter,
                    "OP ",
                ),
                (
                    "Number of optimising solver iterations",
                    "(nviter)",
                    nums.nviter,
                    "OP ",
                ),
                (
                    "Square root of the sum of squares of the constraint residuals",
                    "(sqsumsq)",
                    nums.sqsumsq,
                    "OP ",
                ),
            ):
                process_output.ovarre(constants.NOUT, d, var, v, o)

        process_output.oblnkl(constants.NOUT)

        if self.solver_name == "fsolve":
            process_output.write(
                constants.NOUT,
                "PROCESS has solved using fsolve.\n"
                if ifail == 1
                else "PROCESS failed to solve using fsolve.\n",
            )
        else:
            process_output.write(
                constants.NOUT,
                (
                    (
                        "PROCESS has successfully optimised"
                        if ifail == 1
                        else "PROCESS has failed to optimise"
                    )
                    + " the optimisation parameters to"
                    + ("minimise" if nums.minmax > 0 else "maximise")
                    + f" the objective function: {nums.objf_name}\n"
                ),
            )


@contextmanager
def epsfcn_context(numerics):
    print("Trying again with new epsfcn")
    # epsfcn is only used in evaluators.Evaluators()
    # TODO epsfcn could be set in Evaluators instance now, don't need to
    # set/unset in numerics module
    numerics.epsfcn *= 10  # try new larger value
    print("new epsfcn = ", numerics.epsfcn)
    try:
        yield
    finally:
        numerics.epsfcn /= 10  # reset value
