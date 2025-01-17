"""Unit tests for maths_library.f90."""

import pytest

from process import fortran


def binomial_param(**kwargs):
    """Make parameters for a single binomial() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"n": 1, "k": 1, "expected": 1.0}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def binomial_params():
    """Create a list of parameter dicts for the calc_u_planned fixture.

    Case 1: 1, 1
    Case 2: 3, 1
    Case 3: 3, 3

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        binomial_param(),
        binomial_param(n=3, expected=3.0),
        binomial_param(n=3, k=3, expected=1.0),
    ]


@pytest.fixture(params=binomial_params(), ids=["11", "31", "33"])
def binomial_fix(request):
    """Fixture for the binomial() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :return: Parameterised arguments and expected return value for binomial()
    :rtype: dict
    """
    return request.param


def test_binomial(binomial_fix):
    """Test binomial().

    :param binomial_fix: Parameterised arguments and expected return value for
    binomial()
    :type binomial_fix: dict
    """
    # Run binomial() with the current fixture,
    # then assert the result is the expected one
    n = binomial_fix["n"]
    k = binomial_fix["k"]
    expected = binomial_fix["expected"]

    result = fortran.maths_library.binomial(n, k)
    assert result == expected
