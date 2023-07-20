"""Python counterpart to the Fortran maths library.

It is possible that some functions here are duplicated in both Python and
Fortran; this is a temporary measure during the Python conversion process.
"""

import numpy as np


def secant_solve(f, x1, x2, opt_tol=None):
    """Solve an equation f(x)=0.

    Only requires one function evaluation per iteration (plus two to start with)
    Does not require any derivatives.
    https://en.wikipedia.org/wiki/Secant_method
    Requires two initial values, x0 and x1, which should ideally be chosen to lie close to the root.

    TODO This is a duplicate of maths_library.f90:secant_solve(). The reason for
    this is because Python functions can't be passed to Fortran subroutines via
    the f2py interface, hence a Python function can't be passed into a Fortran
    secant_solve(). Implementing a Python version (and temporarily leaving the
    Fortran one) allows Python and Fortran functions to both use secant_solve()
    as the Python conversion progresses.

    :param f: function to solve
    :type f: funnction, check
    :param x1: first initial value of x
    :type x1: float
    :param x2: second initial value of x
    :type x2: float
    :param opt_tol: optional tolerance, defaults to None
    :type opt_tol: float, optional
    :return: the solution, error boolean and the residual
    :rtype: tuple[float, bool, float]
    """
    error = False
    tol = 0.001
    if opt_tol:
        tol = opt_tol

    x = np.zeros(20)
    x[0] = x1
    x[1] = x2

    # Calculate the first two values before the loop
    fximinus1 = f(x[1])
    fximinus2 = f(x[0])
    for i in range(2, 20):
        # Kludge to avoid zero in denominator
        if fximinus1 == fximinus2:
            fximinus1 = fximinus1 + (fximinus1 * 1e-6)

        x[i] = x[i - 1] - fximinus1 * ((x[i - 1] - x[i - 2]) / (fximinus1 - fximinus2))
        # Store values for the *next iteration
        fximinus2 = fximinus1
        fximinus1 = f(x[i])
        residual = fximinus1

        # test for convergence
        if abs(residual) < tol:
            solution = x[i]
            return solution, error, residual

    # Convergence not achieved.  Return the best value and a flag.
    error = True
    solution = x[i]
    print(f"Secant solver not converged. {solution=} {residual=}")

    return solution, error, residual
