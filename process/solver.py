"""An adapter for different solvers."""

import importlib
import logging
from abc import ABC, abstractmethod

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

from process.data_structure import global_variables, numerics
from process.evaluators import Evaluators
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class _Solver(ABC):
    """Base class for different solver implementations.

    :param ABC: abstract base class
    :type ABC: ABC
    """

    def __init__(self) -> None:
        """Initialise a solver."""
        # Exit code for the solver
        self.ifail = 0
        self.tolerance = numerics.epsvmc
        self.b: float | None = None

    def set_evaluators(self, evaluators: Evaluators) -> None:
        """Set objective and constraint functions and their gradient evaluators.

        :param evaluators: objective and constraint evaluators
        :type evaluators: Evaluators
        """
        self.evaluators = evaluators

    def set_opt_params(self, x_0: np.ndarray) -> None:
        """Define the initial optimisation parameters.

        :param x_0: optimisation parameters vector
        :type x_0: np.ndarray
        """
        self.x_0 = x_0

    def set_bounds(
        self,
        bndl: np.ndarray,
        bndu: np.ndarray,
        ilower: np.ndarray | None = None,
        iupper: np.ndarray | None = None,
    ) -> None:
        """Set the bounds on the optimisation parameters.

        :param bndl: lower bounds for the optimisation parameters
        :type bndl: np.ndarray
        :param bndu: upper bounds for the optimisation parameters
        :type bndu: np.ndarray
        :param ilower: array of 0s and 1s to activate lower bounds on
        optimsation parameters in x
        :type ilower: np.ndarray, optional
        :param iupper: array of 0s and 1s to activate upper bounds on
        optimsation parameters in x
        :type iupper: np.ndarray, optional
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

    def set_constraints(self, m: int, meq: int) -> None:
        """Set the total number of constraints and equality constraints.

        :param m: number of constraint equations
        :type m: int
        :param meq: of the constraint equations, how many are equalities
        :type meq: int
        """
        self.m = m
        self.meq = meq

    def set_tolerance(self, tolerance: float) -> None:
        """Set tolerance for solver termination.

        :param tolerance: tolerance for solver termination
        :type tolerance: float
        """
        self.tolerance = tolerance

    def set_b(self, b: float) -> None:
        """Set the multiplier for the Hessian approximation.

        :param b: multiplier for an identity matrix as input for the Hessian b(n,n)
        :type b: float
        """
        self.b = b

    @abstractmethod
    def solve(self) -> int:
        """Run the optimisation.

        :return: solver error code
        :rtype: int
        """


class VmconProblem(AbstractProblem):
    def __init__(self, evaluator, nequality, ninequality) -> None:
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
    """New VMCON implementation.

    :param _Solver: Solver base class
    :type _Solver: _Solver
    """

    def solve(self) -> int:
        """Optimise using new VMCON.

        :return: solver error code
        :rtype: int
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

            :param result: evaluation of current optimisation parameter vector
            :type result: Result
            :param _x: current optimisation parameter vector
            :type _x: np.ndarray
            :param _delta: search direction for line search
            :type _delta: np.ndarray
            :param _lambda_eq: equality Lagrange multipliers
            :type _lambda_eq: np.ndarray
            :param _lambda_in: inequality Lagrange multipliers
            :type _lambda_in: np.ndarray
            :return: True if inequality constraints satisfied
            :rtype: bool
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
                qsp_options={"adaptive_rho_interval": 25},
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

    def set_opt_params(self, x_0: np.ndarray) -> None:
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
    """Solve equality constraints to ensure model consistency.

    :param _Solver: Solver base class
    :type _Solver: _Solver
    """

    def evaluate_eq_cons(self, x: np.ndarray) -> np.ndarray:
        """Evaluate equality constraints.

        :param x: parameter vector
        :type x: np.ndarray
        :return: equality constraint vector
        :rtype: np.ndarray
        """
        # Evaluate equality constraints only
        _, conf = self.evaluators.fcnvmc1(x.shape[0], self.meq, x, 0)

        return conf

    def solve(self) -> int:
        """Solve equality constraints.

        :return: solver error code
        :rtype: int
        """
        print("Solving equality constraints using fsolve")
        self.x, info, err, msg = fsolve(
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

    :param solver_name: solver to create, defaults to "vmcon"
    :type solver_name: str, optional
    :return: solver to use for optimisation
    :rtype: _Solver
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
    attribute of that package or raise an `AttributeError`."""
    module = importlib.import_module(package)

    solver = getattr(module, "__process_solver__", None)

    if solver is None:
        raise AttributeError(
            f"Module {module.__name__} does not have a '__process_solver__' attribute."
        )

    return solver()
