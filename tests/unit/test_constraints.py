import contextlib

import pytest

from process.core.exceptions import ProcessValueError
from process.core.solver.constraints import ConstraintManager
from process.init import init_all_module_vars


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
