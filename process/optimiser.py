from process.fortran import numerics
from process.solver import get_solver
from process.fortran import define_iteration_variables
from process.evaluators import Evaluators


class Optimiser:
    def __init__(self, models, solver_name):
        """Creates and runs a Vmcon instance.

        This routine calls the minimisation/maximisation routine VMCON,
        developed by Argonne National Laboratory.
        On exit, the (normalised) value of the variable being maximised
        or minimised (i.e. the figure of merit) is returned in argument f.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code.

        This represents the old optimiz subroutine in the numerics module.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param solver_name: which solver to use, as specified in solver.py
        :type solver_name: str
        """
        self.models = models
        self.solver_name = solver_name

    def run(self):
        """Run vmcon solver and retry if it fails in certain ways."""
        # Initialise iteration variables and bounds in Fortran
        define_iteration_variables.loadxc()
        define_iteration_variables.boundxc()

        # Initialise iteration variables and bounds in Python: relies on Fortran
        # iteration variables being defined above
        # Trim maximum size arrays down to actually used size
        n = numerics.nvar
        x = numerics.xcm[:n]
        bndl = numerics.bondl[:n]
        bndu = numerics.bondu[:n]

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

        # If fail then alter value of epsfcn - this can be improved
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
