"""It is well-known at this point that PROCESS produces ever so slight different results depending upon
the operating system, system architecture, and the software versions installed on the system. This is
not surprising nor massively concerning for numeric scientific software. However, it can be problematic
for tests where the results we assert on change from system to system.

We identify systems that will be susceptible to these tiny changes and skip the tests so as not to
create failures. We do this by checking the users systems results on an algorithm compared to some
known result.
This algorithm solves a system of equations that is very ill-conditioned where tiny numeric differences
can be identified and compared against a known system.
"""

import numpy as np
import pytest


def _hilbert_matrix(n):
    """A matrix that is known to be ill-conditioned"""
    return np.array(
        [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)], dtype=np.float64
    )


def system_compatible():
    """Check the system compatibility with numerically unstable tests."""
    # create an ill-conditioned system
    a = _hilbert_matrix(10)
    b = np.ones(10)

    # solve the ill-conditioned system
    x = np.linalg.solve(a, b)

    # create a similar ill-conditioned system
    a_perturbed = a.copy()
    a_perturbed[0, 0] += np.finfo(float).eps
    x_perturbed = np.linalg.solve(a_perturbed, b)

    # return the relative difference between the two solutions
    rel_diff = np.linalg.norm(x - x_perturbed) / np.linalg.norm(x)

    # check that relative difference against a known solution
    return rel_diff == pytest.approx(5.8813403984318865e-05)
