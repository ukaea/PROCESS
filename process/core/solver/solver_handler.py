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
        load_iteration_variables(self.data)
        load_scaled_bounds(self.data)

        n = self.data.numerics.nvar
        x = self.data.numerics.xcm[:n]
        bndl = self.data.numerics.itv_scaled_lower_bounds[:n]
        bndu = self.data.numerics.itv_scaled_upper_bounds[:n]

        # Define total number of constraints and equality constraints
        m = self.data.numerics.neqns + self.data.numerics.nineqns
        meq = self.data.numerics.neqns

        # Evaluators() calculates the objective and constraint functions and
        # their gradients for a given vector x
        evaluators = Evaluators(self.models, self.data, x)

        # Configure solver for problem
        self.solver = get_solver(self.data, self.solver_name)
        self.solver.set_evaluators(evaluators)
        self.solver.set_bounds(bndl, bndu)
        self.solver.set_opt_params(x)
        self.solver.set_constraints(m, meq)
        ifail = self.solver.solve()

        # If VMCON optimisation has failed then try altering value of epsfcn
        if self.solver_name == "vmcon":
            if ifail != 1:
                print("Trying again with new epsfcn")
                # epsfcn is only used in evaluators.Evaluators()
                # TODO epsfcn could be set in Evaluators instance now, don't need to
                # set/unset in self.data.numerics module
                self.data.numerics.epsfcn *= 10  # try new larger value
                print("new epsfcn = ", self.data.numerics.epsfcn)

                ifail = self.solver.solve()
                # First solution attempt failed (ifail != 1): supply ifail value
                # to next attempt
                self.data.numerics.epsfcn /= 10  # reset value

            if ifail != 1:
                print("Trying again with new epsfcn")
                self.data.numerics.epsfcn /= 10  # try new smaller value
                print("new epsfcn = ", self.data.numerics.epsfcn)
                ifail = self.solver.solve()
                self.data.numerics.epsfcn *= 10  # reset value

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
        self.data.numerics.xcm[: self.solver.x.shape[0]] = self.solver.x
        self.data.numerics.rcm[: self.solver.conf.shape[0]] = self.solver.conf
