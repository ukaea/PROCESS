"""Integration tests for VMCON.

Cases 1 to 3 are recommended, others are included to probe the code's
behaviour with different initial guesses for the solution vector x
Expected answers for tests 1 to 3 are given in
VMCON documentation ANL-80-64
"""

import logging
from abc import ABC, abstractmethod

import numpy as np
import pytest

from process.evaluators import Evaluators
from process.init import init_all_module_vars
from process.solver import get_solver

# Debug-level terminal output logging
logger = logging.getLogger(__name__)


@pytest.fixture(autouse=True)
def reinit():
    """Re-initialise Fortran module variables before each test is run."""
    init_all_module_vars()


class Case:
    """A Vmcon test case."""

    def __init__(self, name, evaluator):
        """Initialise name, Vmcon and expected result objects.

        :param name: name of test case
        :type name: str
        :param vmcon: custom Vmcon object for this test
        :type vmcon: Vmcon
        """
        self.name = name
        self.solver_args = SolverArgs()
        self.evaluator = evaluator
        self.exp = ExpectedResult()


class SolverArgs:
    def __init__(self, n=2):
        """Initialise some common arguments to the solver adapter.

        :param n: number of input variables (default 2)
        :type n: int

        These arguments are shared between some of the test cases.
        """
        # No bounds on x values set
        self.x = np.zeros(n)
        self.ilower = np.zeros(n)
        self.iupper = np.zeros(n)
        self.bndl = np.zeros(n)
        self.bndu = np.full(n, 5.0)
        self.tolerance = 1.0e-8


class ExpectedResult:
    """Expected result class for comparing an observed solver result against."""

    def __init__(self):
        """Initialise expected attributes."""
        # TODO Some of these are Vmcon-specific
        self.x = []
        self.c = []
        self.vlam = []
        self.objf = None
        self.errlg = 0.0
        self.errlm = 0.0
        self.errcom = 0.0
        self.errcon = 0.0
        self.ifail = 1


class CustomFunctionEvaluator(ABC, Evaluators):
    """Abstract function evaluator for different solver test cases to override.

    :param ABC: abstract base class
    :type ABC: abc.ABC
    :param Evaluators: objective and constraint function and gradient evaluator
    :type Evaluators: process.evaluators.Evaluators
    """

    def __init__(self):
        """Override to prevent Caller() to physics and engineering models being
        initialised.
        """

    @abstractmethod
    def fcnvmc1(self):
        """Function evaluator."""

    @abstractmethod
    def fcnvmc2(self):
        """Gradient function evaluator."""


class Evaluator1(CustomFunctionEvaluator):
    """Override fcnvmc1 and 2 methods for test case 1.

    This allows a test to be run using custom function and gradient function
    evaluators specific to this test.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 - 2*x2 + 1 = 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0
    :param CustomFunctionEvaluator: abstract function evaluator
    :type CustomFunctionEvaluator: CustomFunctionEvaluator
    """

    def fcnvmc1(self, n, m, x, ifail):
        """Function evaluator.


        Calculates the objective and constraint functions at the
        n-dimensional point of interest x.
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple(float, np.ndarray)
        """
        objf = (x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2
        conf = np.array([
            x[0] - 2.0 * x[1] + 1.0,
            -0.25 * x[0] ** 2 - x[1] * x[1] + 1.0,
        ])

        return objf, conf

    def fcnvmc2(self, n, m, x, lcnorm):
        """Gradient function evaluator.

        Calculates the gradients of the objective and constraint functions at
        the n-dimensional point of interest x. The constraint gradients or
        normals are returned as the columns of cnorm.

        :param n: number of iteration variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrd (np.ndarray (n)) gradient of the objective function
        cnorm (np.ndarray (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fgrd = np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)])

        cnorm = np.array([[1.0, -0.5 * x[0]], [-2.0, -2.0 * x[1]]], order="F")
        return fgrd, cnorm


class Evaluator2(CustomFunctionEvaluator):
    """Override fcnvmc1 and 2 methods for test case 2.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 - 2*x2 + 1 >= 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0
    :param CustomFunctionEvaluator: abstract function evaluator
    :type CustomFunctionEvaluator: CustomFunctionEvaluator
    """

    def fcnvmc1(self, n, m, x, ifail):
        """Function evaluator.


        Calculates the objective and constraint functions at the
        n-dimensional point of interest x.
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple(float, np.ndarray)
        """
        objf = (x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2
        conf = np.array([
            x[0] - 2.0 * x[1] + 1.0,
            -0.25 * x[0] ** 2 - x[1] * x[1] + 1.0,
        ])

        return objf, conf

    def fcnvmc2(self, n, m, x, lcnorm):
        """Gradient function evaluator.

        Calculates the gradients of the objective and constraint functions at
        the n-dimensional point of interest x. The constraint gradients or
        normals are returned as the columns of cnorm.

        :param n: number of iteration variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrd (np.ndarray (n)) gradient of the objective function
        cnorm (np.ndarray (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fgrd = np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)])
        cnorm = np.array([[1.0, -0.5 * x[0]], [-2.0, -2.0 * x[1]]])

        return fgrd, cnorm


class Evaluator3(CustomFunctionEvaluator):
    """Override fcnvmc1 and 2 methods for test case 3.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 + x2 - 3 = 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0

    :param CustomFunctionEvaluator: abstract function evaluator
    :type CustomFunctionEvaluator: CustomFunctionEvaluator
    """

    def fcnvmc1(self, n, m, x, ifail):
        """Function evaluator.


        Calculates the objective and constraint functions at the
        n-dimensional point of interest x.
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple(float, np.ndarray)
        """
        objf = (x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2
        conf = np.array([x[0] + x[1] - 3.0, -0.25 * x[0] ** 2 - x[1] * x[1] + 1.0])

        return objf, conf

    def fcnvmc2(self, n, m, x, lcnorm):
        """Gradient function evaluator.

        Calculates the gradients of the objective and constraint functions at
        the n-dimensional point of interest x. The constraint gradients or
        normals are returned as the columns of cnorm.

        :param n: number of iteration variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrd (np.ndarray (n)) gradient of the objective function
        cnorm (np.ndarray (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fgrd = np.array([2.0, 2.0 * (x[1] - 1.0)])
        cnorm = np.array([[1.0, -0.5 * x[0]], [1.0, -2.0 * x[1]]])

        return fgrd, cnorm


class Evaluator4(CustomFunctionEvaluator):
    """Override fcnvmc1 and 2 methods for test case 4.

    From Wikipedia: Lagrange Multiplier article
    Maximise f(x1,x2) = x1 + x2
    subject to the following constraint:
    c1(x1,x2) = x1**2 + x2**2 - 1 = 0

    :param CustomFunctionEvaluator: abstract function evaluator
    :type CustomFunctionEvaluator: CustomFunctionEvaluator
    """

    def fcnvmc1(self, n, m, x, ifail):
        """Function evaluator.

        Calculates the objective and constraint functions at the
        n-dimensional point of interest x.
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple(float, np.ndarray)
        """
        objf = x[0] + x[1]
        conf = np.array([x[0] * x[0] + x[1] * x[1] - 1.0])

        return objf, conf

    def fcnvmc2(self, n, m, x, lcnorm):
        """Gradient function evaluator.

        Calculates the gradients of the objective and constraint functions at
        the n-dimensional point of interest x. The constraint gradients or
        normals are returned as the columns of cnorm.

        :param n: number of iteration variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrd (np.ndarray (n)) gradient of the objective function
        cnorm (np.ndarray (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fgrd = np.array([1.0, 1.0])
        cnorm = np.array([[2.0 * x[0]], [2.0 * x[1]]])

        return fgrd, cnorm


class Evaluator5(CustomFunctionEvaluator):
    """Override fcnvmc1 and 2 methods for test case 5.

    :param CustomFunctionEvaluator: abstract function evaluator
    :type CustomFunctionEvaluator: CustomFunctionEvaluator
    """

    def fcnvmc1(self, n, m, x, ifail):
        """Function evaluator.

        Calculates the objective and constraint functions at the
        n-dimensional point of interest x.
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple(float, np.ndarray)
        """
        objf = x[0] ** 2
        conf = np.array([x[0] ** 2 - 2.0 * x[0] - 3.0])

        return objf, conf

    def fcnvmc2(self, n, m, x, lcnorm):
        """Gradient function evaluator.

        Calculates the gradients of the objective and constraint functions at
        the n-dimensional point of interest x. The constraint gradients or
        normals are returned as the columns of cnorm.

        :param n: number of iteration variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param x: iteration variable array, length n
        :type x: np.ndarray
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrd (np.ndarray (n)) gradient of the objective function
        cnorm (np.ndarray (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fgrd = np.array([2.0 * x[0]])
        cnorm = np.array([[2.0 * x[0] - 2.0]])

        return fgrd, cnorm


def get_case1():
    """Create test case 1 for the solver.

    Set up the problem and define the expected result.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 - 2*x2 + 1 = 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0

    VMCON documentation ANL-80-64
    """
    # Create a case-specific Evaluator object with overridden fcnvmc1 and 2
    case = Case("1", Evaluator1())

    # Set up solver args for this case
    neqns = 1
    nineqns = 1
    case.solver_args.x[0:2] = 2.0e0
    case.solver_args.n = 2
    case.solver_args.m = neqns + nineqns
    case.solver_args.meq = neqns

    # Expected values
    case.exp.x = np.array([8.228756e-1, 9.114378e-1])
    case.exp.objf = 1.393464
    case.exp.c = np.array([1.387778e-17, -7.671641e-13])
    case.exp.vlam = np.array([-1.594491, 1.846591])
    case.exp.errlg = 3.345088e-12
    case.exp.errcom = 1.416660e-12
    case.exp.errcon = 7.671779e-13

    return case


def get_case2():
    """Create test case 2 for the solver.

    Set up the problem and define the expected result.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 - 2*x2 + 1 >= 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0

    VMCON documentation ANL-80-64
    """
    case = Case("2", Evaluator2())

    # Solver args for this case
    neqns = 0
    nineqns = 2
    case.solver_args.n = 2
    case.solver_args.m = neqns + nineqns
    case.solver_args.meq = neqns
    case.solver_args.x[0:2] = 2.0e0

    # Expected values
    case.exp.x = np.array([1.664968, 5.540486e-1])
    case.exp.objf = 3.111186e-1
    case.exp.c = np.array([1.556871, -1.021405e-14])
    case.exp.vlam = np.array([0.0, 8.048955e-1])
    case.exp.errlg = 2.343333e-11
    case.exp.errcom = 8.221245e-15
    case.exp.errcon = 1.021405e-14

    return case


def get_case3():
    """Create test case 3 for the solver.

    Set up the problem and define the expected result.

    Minimise f(x1,x2) = (x1 - 2)**2 + (x2 - 1)**2
    subject to the following constraints:
    c1(x1,x2) = x1 + x2 - 3 = 0
    c2(x1,x2) = -x1**2/4 - x2**2 + 1 >= 0

    Note that this test is supposed to fail with ifail=5
    as there is no feasible solution
    VMCON documentation ANL-80-64
    """
    # Create a case-specific Vmcon object with overridden fcnvmc1 and 2
    case = Case("3", Evaluator3())

    # Solver args for this case
    neqns = 1
    nineqns = 1
    case.solver_args.n = 2
    case.solver_args.m = neqns + nineqns
    case.solver_args.meq = neqns
    case.solver_args.x[0:2] = 2.0e0

    # Expected values
    case.exp.x = np.array([2.0, 2.0])
    case.exp.objf = 6.203295e-1
    case.exp.c = np.array([-6.6613381477509392e-16, -8.000000000004035e-01])
    case.exp.vlam = np.array([0.0, 0.0])
    case.exp.errlg = 1.599997724349894
    case.exp.errcon = 8.0000000000040417e-01
    case.exp.ifail = 5

    return case


def get_case4():
    """Create test case 4 for the solver.

    Set up the problem and define the expected result.

    Maximise f(x1,x2) = x1 + x2
    subject to the following constraint:
    c1(x1,x2) = x1**2 + x2**2 - 1 = 0

    http://en.wikipedia.org/wiki/Lagrange_multiplier
    """
    # Create a case-specific Vmcon object with overridden fcnvmc1 and 2
    case = Case("4", Evaluator4())

    # Set up vmcon values for this case
    neqns = 1
    nineqns = 0
    case.solver_args.n = 2
    case.solver_args.m = neqns + nineqns
    case.solver_args.meq = neqns
    case.solver_args.xtol = 2.0e-8
    case.solver_args.x[0:2] = 1.0e0
    # N.B. results can flip to minimum instead of maximum
    # if x(1), x(2) are initialised at different points...

    # Expected values
    case.exp.x = np.array([0.5 * 2.0 ** (1 / 2), 0.5 * 2.0 ** (1 / 2)])
    case.exp.objf = 2.0 ** (1 / 2)
    case.exp.c = np.array([0.0])
    case.exp.vlam = np.array([1.0 * 2.0 ** (1 / 2)])

    return case


def get_case5():
    """Create test case 5 for the solver.

    Set up the problem and define the expected result.

    Intersection of parabola x^2 with straight line 2x+3
    Unorthodox (and not recommended) method to find the root
    of an equation.

    Maximise f(x1) = x1**2
    subject to the following constraint:
    c1(x1) = x1**2 - 2.0D0*x1 - 3 = 0

    Solutions to c1(x1) are x1 = -1 and x1 = 3, and depending on
    the initial guess for x1 either (or neither...) solution might
    be found. Since there is one constraint equation with one unknown
    the code cannot optimise properly.
    """
    # Create a case-specific Vmcon object with overridden fcnvmc1 and 2
    case = Case("5", Evaluator5())
    case.solver_args = SolverArgs(n=1)

    # Set up vmcon values for this case
    neqns = 1
    nineqns = 0
    case.solver_args.n = 1
    case.solver_args.m = neqns + nineqns
    case.solver_args.meq = neqns
    case.solver_args.x = np.array([
        5.0
    ])  # Try different values, e.g. 5.0, 2.0, 1.0, 0.0...
    case.solver_args.bndl = np.zeros(1)
    case.solver_args.bndu = np.full(1, 5.0)
    case.solver_args.ilower = np.zeros(1)
    case.solver_args.iupper = np.zeros(1)

    # Expected values
    case.exp.x = np.array([3.0])
    case.exp.objf = 9.0
    case.exp.c = np.array([0.0])
    case.exp.vlam = np.array([1.5])

    return case


@pytest.fixture(params=[get_case1, get_case2, get_case3, get_case4, get_case5])
def case(request):
    """Parameterised fixture for providing Vmcon test cases to run.

    :param request: provides access to various test case functions
    :type request: SubRequest
    :return: Vmcon scenario to run and its expected result
    :rtype: test_vmcon.Case
    """
    case_fn = request.param
    return case_fn()


def test_vmcon(case, solver_name):
    """Integration test for Vmcon.

    :param case: a Vmcon scenario and its expected result
    :type case: test_vmcon.Case
    :param solver_name: solver to use (from pytest CLI option)
    :type solver_name: str
    """
    if case.name == "3":
        pytest.skip("Test case 3 currently fails with new VMCON.")

    logger.debug("Initial solution estimate:")
    for i in range(case.solver_args.n):
        logger.debug(f"x[{i}] = {case.solver_args.x[i]}")

    # Configure solver for problem
    solver = get_solver(solver_name)
    solver.set_evaluators(case.evaluator)
    solver.set_opt_params(case.solver_args.x)
    solver.set_bounds(
        case.solver_args.bndl,
        case.solver_args.bndu,
        ilower=case.solver_args.ilower,
        iupper=case.solver_args.iupper,
    )
    solver.set_constraints(case.solver_args.m, case.solver_args.meq)
    solver.set_tolerance(case.solver_args.tolerance)

    # Run solver for this case
    info = solver.solve()
    x = solver.x
    objf = solver.objf

    # Assert result
    try:
        # Assert ifail is expected value (converged: ifail == 1)
        assert info == case.exp.ifail

        # Assert final objective function value
        assert objf == pytest.approx(case.exp.objf)

        # Final solution estimate
        assert x == pytest.approx(case.exp.x)

    except AssertionError:
        logger.exception(f"Vmcon test {case.name} failed")

        # Detailed debugging only
        # log_failure(case)
        raise


def log_failure(case):
    """Write extra logs in the case of a test case failure.

    :param case: a failing test case
    :type case: Case
    """
    # TODO Reinstate detailed debugging now solver adapter is used rather than
    # Vmcon() directly
    logger.debug(f"ifail = {case.vmcon.ifail} (expected value = {case.exp.ifail}")
    logger.debug(f"Number of function evaluations = {case.vmcon.fcnvmc1_calls}")

    logger.debug("Final solution estimate: calculated vs expected")
    for i in range(case.vmcon.n):
        logger.debug(f"x[{i}] = {case.vmcon.x[i]}, {case.exp.x[i]}")

    logger.debug(
        f"Final objective function value: calculated vs expected: "
        f"{case.vmcon.objf}, {case.exp.objf}"
    )

    logger.debug("Constraints evaluated at x: calculated vs expected':")
    for i in range(case.vmcon.m):
        logger.debug(f"{case.vmcon.conf[i]}, {case.exp.c[i]}")

    logger.debug("Lagrange multiplier estimates: calculated vs expected")
    for i in range(case.vmcon.m):
        logger.debug(f"{case.vmcon.vlam[i]}, {case.exp.vlam[i]}")

    logger.debug("Lagrangian gradient error: calculated vs expected")
    errlg = 0.0
    for i in range(case.vmcon.n):
        summ = case.vmcon.fgrd[i]
        for j in range(case.vmcon.m):
            summ = summ - case.vmcon.vlam[j] * case.vmcon.cnorm[i, j]
        errlg = errlg + abs(summ)
    logger.debug(f"{errlg}, {case.exp.errlg}")

    logger.debug("Lagrange multiplier error: calculated vs expected")
    errlm = 0.0
    for i in range(case.vmcon.m):
        if (i <= case.vmcon.meq) or (case.vmcon.vlam[i] >= 0.0):
            continue
        errlm = errlm + abs(case.vmcon.vlam[i])
    logger.debug(f"{errlm}, {case.exp.errlm}")

    logger.debug("Complementarity error: calculated vs expected")
    errcom = 0.0
    for i in range(case.vmcon.m):
        errcom = errcom + abs(case.vmcon.vlam[i] * case.vmcon.conf[i])
    logger.debug(f"{errcom}, {case.exp.errcom}")

    logger.debug("Constraint error: calculated vs expected")
    errcon = 0.0
    for i in range(case.vmcon.m):
        if (i > case.vmcon.meq) and (case.vmcon.conf[i] >= 0.0):
            continue
        errcon = errcon + abs(case.vmcon.conf[i])
    logger.debug(f"{errcon}, {case.exp.errcon}")
