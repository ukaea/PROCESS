"""An adapter for different solvers."""

import logging
from process.fortran import numerics, global_variables
from process.utilities.f2py_string_patch import f2py_compatible_to_string
import numpy as np
from process.evaluators import Evaluators
from abc import ABC, abstractmethod
from typing import Optional, Union
import importlib
from pyvmcon import (
    AbstractProblem,
    Result,
    solve,
    QSPSolverException,
    VMCONConvergenceException,
    LineSearchConvergenceException,
)

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
        self.b: Union[float, None] = None

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
        ilower: Optional[np.ndarray] = None,
        iupper: Optional[np.ndarray] = None,
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
            ilower = np.ones(bndl.shape[0], dtype=int)
            iupper = np.ones(bndu.shape[0], dtype=int)

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
        pass


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

        :raises NotImplementedError: not currently implemented
        :return: solver error code
        :rtype: int
        """
        problem = VmconProblem(self.evaluators, self.meq, self.m - self.meq)

        B = None
        if self.b is not None:
            B = np.identity(numerics.nvar) * self.b

        def _solver_callback(i: int, _result, _x, convergence_param: float):
            numerics.nviter = i + 1
            global_variables.convergence_parameter = convergence_param
            print(
                f"{i+1} | Convergence Parameter: {convergence_param:.3E}",
                end="\r",
                flush=True,
            )

        try:
            x, _, _, res = solve(
                problem,
                self.x_0,
                self.bndl,
                self.bndu,
                max_iter=global_variables.maxcal,
                epsilon=self.tolerance,
                qsp_options={"eps_rel": 1e-1, "adaptive_rho_interval": 25},
                initial_B=B,
                callback=_solver_callback,
            )
        except VMCONConvergenceException as e:
            if isinstance(e, LineSearchConvergenceException):
                self.info = 3
            elif isinstance(e, QSPSolverException):
                self.info = 5
            else:
                self.info = 2

            logger.warning(str(e))

            x = e.x
            res = e.result

        except ValueError as e:
            itervar_name_list = ""
            for count, iter_var in enumerate(numerics.ixc[: numerics.nvar]):
                itervar_name = f2py_compatible_to_string(numerics.lablxc[iter_var - 1])
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
    else:
        try:
            solver = load_external_solver(solver_name)
        except Exception as e:
            raise ValueError(
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
