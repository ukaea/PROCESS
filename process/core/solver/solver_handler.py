from contextlib import contextmanager

from process.core.solver.evaluators import Evaluators
from process.core.solver.iteration_variables import (
    load_iteration_variables,
    load_scaled_bounds,
)
from process.core.solver.solver import get_solver


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
            m=self.data.numerics.neqns + self.data.numerics.nineqns, meq=self.data.numerics.neqns
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
