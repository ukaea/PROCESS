import pytest
import numpy as np
from numpy import typing as npt

from process.neutronics_data import scattering_weight_matrix, _convolved_scattering_fraction

np.random.seed(1)
@pytest.mark.parametrize("group_structure, atomic_mass",
    [
        (np.cumsum(100**abs(np.random.randn(5)))[::-1], np.random.rand()*200),
        (np.cumsum(100**abs(np.random.randn(10)))[::-1], np.random.rand()*239),
        (np.cumsum(100**abs(np.random.randn(20)))[::-1], np.random.rand()*150),
        (np.cumsum(100**abs(np.random.randn(30)))[::-1], 1),
        (np.cumsum(100**abs(np.random.randn(40)))[::-1], 2),
    ]
)
def test_scattering_matrix(group_structure: npt.NDArray[np.float64], atomic_mass: float):
    """
    Check that even randomly generated scattering matrix are normalized
    and non-negative.
    """
    scattering_matrix = scattering_weight_matrix(group_structure, atomic_mass)
    row_sum = scattering_matrix.sum(axis=1)
    assert np.logical_or(np.isclose(row_sum, 1, rtol=0), row_sum<1).all(), (
        "Every row must be unitary or less."
    )
    assert (scattering_matrix>=0).all(), "Must be all non-negative."
