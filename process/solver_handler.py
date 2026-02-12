from process.data_structure import numerics
from process.evaluators import Evaluators
from process.iteration_variables import load_iteration_variables, load_scaled_bounds
from process.solver import get_solver


class SolverHandler:
    """Creates and runs a solver instance.

    This may be an optimiser (e.g. VMCON) or an equation solver (e.g. fsolve).

    Parameters
    ----------
    models : process.main.Models
        physics and engineering model objects
    solver_name : str
        which solver to use, as specified in solver.py
    """

    def __init__(self, models, solver_name):
        self.models = models
        self.solver_name = solver_name

    def run(self):
        """Run solver and retry if it fails in certain ways."""
        # Initialise iteration variables and bounds in Fortran
        load_iteration_variables()
        load_scaled_bounds()

        # Initialise iteration variables and bounds in Python: relies on Fortran
        # iteration variables being defined above
        # Trim maximum size arrays down to actually used size
        n = numerics.nvar
        x = numerics.xcm[:n]
        bndl = numerics.itv_scaled_lower_bounds[:n]
        bndu = numerics.itv_scaled_upper_bounds[:n]

        # Define total number of constraints and equality constraints
        m = numerics.neqns + numerics.nineqns
        meq = numerics.neqns

        # Evaluators() calculates the objective and constraint functions and
        # their gradients for a given vector x
        evaluators = Evaluators(self.models, x)

        # Configure solver for problem
        self.solver = get_solver(self.solver_name)
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
                # set/unset in numerics module
                numerics.epsfcn = numerics.epsfcn * 10  # try new larger value
                print("new epsfcn = ", numerics.epsfcn)

                ifail = self.solver.solve()
                # First solution attempt failed (ifail != 1): supply ifail value
                # to next attempt
                numerics.epsfcn = numerics.epsfcn / 10  # reset value

            if ifail != 1:
                print("Trying again with new epsfcn")
                numerics.epsfcn = numerics.epsfcn / 10  # try new smaller value
                print("new epsfcn = ", numerics.epsfcn)
                ifail = self.solver.solve()
                numerics.epsfcn = numerics.epsfcn * 10  # reset value

            # If VMCON has exited with error code 5 try another run using a multiple
            # of the identity matrix as input for the Hessian b(n,n)
            # Only do this if VMCON has not iterated (nviter=1)
            if ifail == 5 and numerics.nviter < 2:
                print(
                    "VMCON error code = 5.  Rerunning VMCON with a new initial "
                    "estimate of the second derivative matrix."
                )
                self.solver.set_b(2.0)
                ifail = self.solver.solve()

        self.output()

        return ifail

    def output(self):
        """Store results back in Fortran numerics module.

        Objective function value, solution vector and constraints vector.
        """
        numerics.norm_objf = self.solver.objf
        # Slicing required due to Fortran arrays being maximum possible, rather
        # than required, size
        numerics.xcm[: self.solver.x.shape[0]] = self.solver.x
        numerics.rcm[: self.solver.conf.shape[0]] = self.solver.conf
