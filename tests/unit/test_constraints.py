import contextlib
from operator import eq as python_eq
from operator import ge as python_geq
from operator import le as python_leq

import pytest

from process.core.exceptions import ProcessValueError
from process.core.init import init_all_module_vars
from process.core.solver.constraints import (
    ConstraintManager,
    ConstraintRegistration,
    eq,
    geq,
    leq,
)


@pytest.mark.parametrize(
    ["value", "bound"],
    [
        [10, 5],
        [5, 10],
        [1, 1],
        [10, -5],
        [-5, 10],
        [-1, -1],
        [1.0, 1.001],
        [1.001, 1.0],
        [-1.1, -1.0],
        [-1.0, -1.1],
        # [10, 0],
        # [0, 10],
        # [0, 0],
    ],
)
@pytest.mark.parametrize(
    ["constraint_function", "python_function"],
    [[eq, python_eq], [leq, python_leq], [geq, python_geq]],
)
def test_constraints(value, bound, constraint_function, python_function):
    mock_registration = ConstraintRegistration("test", ..., ..., ...)

    result = constraint_function(value, bound, mock_registration)

    if constraint_function is eq:
        constraint_function_implies_true = result.normalised_residual == 0
    else:
        constraint_function_implies_true = result.normalised_residual <= 0
    python_implies_true = python_function(value, bound)

    assert bool(constraint_function_implies_true) is python_implies_true


@pytest.mark.parametrize(
    "constraint_registration",
    ConstraintManager._constraint_registry.values(),  # noqa: SLF001
    ids=[f"constraint_{i}" for i in ConstraintManager._constraint_registry],  # noqa: SLF001
)
def test_constraint_functions(constraint_registration):
    """A simple test that runs a constraint with PROCESS' default initialisation to check
    that no attribute or type errors (etc) occur.
    """
    init_all_module_vars()

    # Allow zero division errors because this means all of the attributes exist
    # and are initialised with a number.
    # Allow ProcessValueError because that happens when a constraint is incompatible with the
    # default flags (i_pulsed_plant=0 or itart=0).
    with contextlib.suppress(ZeroDivisionError, ProcessValueError):
        # call the constraint equation and check no error occurs
        constraint_registration.constraint_equation(constraint_registration)
