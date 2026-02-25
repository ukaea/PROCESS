"""An adapter for different solvers."""

import importlib
import logging
from abc import ABC, abstractmethod

import cvxpy
import numpy as np
from pyvmcon import (
    AbstractProblem,
    LineSearchConvergenceException,
    QSPSolverException,
    Result,
    VMCONConvergenceException,
    solve,
)
from scipy.optimize import fsolve

from process.core.solver.evaluators import Evaluators
from process.data_structure import global_variables, numerics
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class _Solver(ABC):
    """Base class for different solver implementations."""

    def __init__(self):
        """Initialise a solver."""
        # Exit code for the solver
        self.ifail = 0
        self.tolerance = numerics.epsvmc
        self.b: float | None = None

    def set_evaluators(self, evaluators: Evaluators):
        """Set objective and constraint functions and their gradient evaluators.

        Parameters
        ----------
        evaluators : Evaluators
            objective and constraint evaluators

        """
        self.evaluators = evaluators

    def set_opt_params(self, x_0: np.ndarray):
        """Define the initial optimisation parameters.

        Parameters
        ----------
        x_0 : np.ndarray
            optimisation parameters vector

        """
        self.x_0 = x_0

    def set_bounds(
        self,
        bndl: np.ndarray,
        bndu: np.ndarray,
        ilower: np.ndarray | None = None,
        iupper: np.ndarray | None = None,
    ):
        """Set the bounds on the optimisation parameters.

        Parameters
        ----------
        bndl : np.ndarray
            lower bounds for the optimisation parameters
        bndu : np.ndarray
            upper bounds for the optimisation parameters
        ilower : np.ndarray, optional
            array of 0s and 1s to activate lower bounds on
            optimsation parameters in x
        iupper : np.ndarray, optional
            array of 0s and 1s to activate upper bounds on
            optimsation parameters in x
        """
        self.bndl = bndl
        self.bndu = bndu

        # TODO Remove ilower/iupper and use finite vs. infinite values in bndl/bndu
        # instead to determine defined vs. undefined bounds
        # If lower and upper bounds switch arrays aren't specified, set all bounds
        # to defined
        if ilower is None and iupper is None:
            ilower = np.ones(len(bndl), dtype=int)
            iupper = np.ones(len(bndu), dtype=int)

        self.ilower = ilower
        self.iupper = iupper

    def set_constraints(self, m: int, meq: int):
        """Set the total number of constraints and equality constraints.

        Parameters
        ----------
        m : int
            number of constraint equations
        meq : int
            of the constraint equations, how many are equalities
        """
        self.m = m
        self.meq = meq

    def set_tolerance(self, tolerance: float):
        """Set tolerance for solver termination.

        Parameters
        ----------
        tolerance : float
            tolerance for solver termination
        """
        self.tolerance = tolerance

    def set_b(self, b: float):
        """Set the multiplier for the Hessian approximation.

        Parameters
        ----------
        b : float
            multiplier for an identity matrix as input for the Hessian b(n,n)
        """
        self.b = b

    @abstractmethod
    def solve(self) -> int:
        """Run the optimisation.

        Returns
        -------
        int
            solver error code
        """


class VmconProblem(AbstractProblem):
    def __init__(self, evaluator, nequality, ninequality):
        self._evaluator = evaluator
        self._nequality = nequality
        self._ninequality = ninequality

    def __call__(self, x: np.ndarray) -> Result:
        n = x.shape[0]
        objf, conf = self._evaluator.fcnvmc1(n, self.total_constraints, x, 0)
        fgrd, cnorm = self._evaluator.fcnvmc2(n, self.total_constraints, x, n)

        return Result(
            objf,
            fgrd,
            conf[: self.num_equality],
            cnorm[:, : self.num_equality].T,
            conf[self.num_equality :],
            cnorm[:, self.num_equality :].T,
        )

    @property
    def num_equality(self) -> int:
        return self._nequality

    @property
    def num_inequality(self) -> int:
        return self._ninequality


class Vmcon(_Solver):
    """New VMCON implementation."""

    def solve(self) -> int:
        """Optimise using new VMCON.

        Returns
        -------
        int
            solver error code
        """
        problem = VmconProblem(self.evaluators, self.meq, self.m - self.meq)

        bb = None
        if self.b is not None:
            bb = np.identity(numerics.nvar) * self.b

        def _solver_callback(i: int, _result, _x, convergence_param: float):
            numerics.nviter = i + 1
            global_variables.convergence_parameter = convergence_param
            print(
                f"{i + 1} | Convergence Parameter: {convergence_param:.3E}",
                end="\r",
                flush=True,
            )

        def _ineq_cons_satisfied(
            result: Result,
            _x: np.ndarray,
            _delta: np.ndarray,
            _lambda_eq: np.ndarray,
            _lambda_in: np.ndarray,
        ) -> bool:
            """Check that inequality constraints are satisfied.

            This additional convergence criterion ensures that solutions have
            satisfied inequality constraints.

            Parameters
            ----------
            result : Result
                evaluation of current optimisation parameter vector
            _x : np.ndarray
                current optimisation parameter vector
            _delta : np.ndarray
                search direction for line search
            _lambda_eq : np.ndarray
                equality Lagrange multipliers
            _lambda_in : np.ndarray
                inequality Lagrange multipliers

            Returns
            -------
            bool
                True if inequality constraints satisfied

            """
            # negative constraint value = violated
            # Check all ineqs are satisfied to within the tolerance
            # E.g. the relative violations are no more than v=0-tolerance
            return bool(np.all(result.ie >= -numerics.force_vmcon_inequality_tolerance))

        try:
            x, _, _, res = solve(
                problem,
                np.array(self.x_0),
                np.array(self.bndl),
                np.array(self.bndu),
                max_iter=global_variables.maxcal,
                epsilon=self.tolerance,
                qsp_options={"solver": cvxpy.CLARABEL},
                initial_B=bb,
                callback=_solver_callback,
                additional_convergence=_ineq_cons_satisfied
                if numerics.force_vmcon_inequality_satisfication
                else None,
            )
        except VMCONConvergenceException as e:
            if isinstance(e, LineSearchConvergenceException):
                self.info = 3
            elif isinstance(e, QSPSolverException):
                self.info = 5
            else:
                self.info = 2

            logger.critical(str(e))

            x = e.x
            res = e.result

        except ValueError as e:
            itervar_name_list = ""
            for count, iter_var in enumerate(numerics.ixc[: numerics.nvar]):
                itervar_name = numerics.lablxc[iter_var - 1]
                itervar_name_list += f"{count}: {itervar_name} \n"

            logger.warning(f"Active iteration variables are : \n{itervar_name_list}")
            raise e

        else:
            self.info = 1

        # print a blank line because of the carridge return
        # in the callback
        print()

        self.x = x
        self.objf = res.f
        self.conf = np.hstack((res.eq, res.ie))

        return self.info


class VmconBounded(Vmcon):
    """A solver that uses VMCON but checks x is in bounds before running"""

    def set_opt_params(self, x_0: np.ndarray):
        lower_violated = np.less(x_0, self.bndl)
        upper_violated = np.greater(x_0, self.bndu)

        for index, entry in enumerate(lower_violated):
            if entry:
                x_0[index] = self.bndl[index]
        for index, entry in enumerate(upper_violated):
            if entry:
                x_0[index] = self.bndu[index]
        self.x_0 = x_0


class FSolve(_Solver):
    """Solve equality constraints to ensure model consistency."""

    def evaluate_eq_cons(self, x: np.ndarray) -> np.ndarray:
        """Evaluate equality constraints.

        Parameters
        ----------
        x : np.ndarray
            parameter vector

        Returns
        -------
        np.ndarray
            equality constraint vector
        """
        # Evaluate equality constraints only
        _, conf = self.evaluators.fcnvmc1(x.shape[0], self.meq, x, 0)

        return conf

    def solve(self) -> int:
        """Solve equality constraints.

        Returns
        -------
        int
            solver error code
        """
        print("Solving equality constraints using fsolve")
        self.x, _info, err, msg = fsolve(
            self.evaluate_eq_cons, self.x_0, full_output=True
        )

        # Evaluate equality and inequality constraints at equality-satisfying solution
        # (or at last iteration of x if solution not found)
        _, self.conf = self.evaluators.fcnvmc1(self.x.shape[0], self.m, self.x, 0)

        # err == 1 for successful solve
        if err != 1:
            print(f"fsolve error code {err}: {msg}")
        self.info = err
        # No objective function
        self.objf = None
        return self.info


def get_solver(solver_name: str = "vmcon") -> _Solver:
    """Return a solver instance.

    Parameters
    ----------
    solver_name : str, optional
        solver to create, defaults to "vmcon"

    Returns
    -------
    _Solver
        solver to use for optimisation
    """
    solver: _Solver

    if solver_name == "vmcon":
        solver = Vmcon()
    elif solver_name == "vmcon_bounded":
        solver = VmconBounded()
    elif solver_name == "fsolve":
        solver = FSolve()
    else:
        try:
            solver = load_external_solver(solver_name)
        except Exception as e:
            raise ProcessValueError(
                f'Solver name is not an inbuilt PROCESS solver or recognised package "{solver_name}"'
            ) from e

    return solver


def load_external_solver(package: str):
    """Attempts to load a package of name `package`.

    If a package of the name is available, return the `__process_solver__`
    attribute of that package or raise an `AttributeError`.

    Parameters
    ----------
    package: str :

    """
    module = importlib.import_module(package)

    solver = getattr(module, "__process_solver__", None)

    if solver is None:
        raise AttributeError(
            f"Module {module.__name__} does not have a '__process_solver__' attribute."
        )

    return solver()
