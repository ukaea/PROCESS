import numpy as np
import pytest
from numpy import typing as npt

from process.neutronics_data import n2n_weight_matrix, scattering_weight_matrix

np.random.seed(1)


@pytest.mark.parametrize(
    "group_structure, atomic_mass",
    [
        (
            np.cumsum(100 ** abs(np.random.randn(5)))[::-1],
            np.random.rand() * 200,
        ),
        (
            np.cumsum(100 ** abs(np.random.randn(10)))[::-1],
            np.random.rand() * 239,
        ),
        (
            np.cumsum(100 ** abs(np.random.randn(20)))[::-1],
            np.random.rand() * 150,
        ),
        (np.cumsum(100 ** abs(np.random.randn(30)))[::-1], 1),
        (np.cumsum(100 ** abs(np.random.randn(40)))[::-1], 2),
    ],
)
def test_scattering_matrix(
    group_structure: npt.NDArray[np.float64], atomic_mass: float
):
    """
    Check that even randomly generated scattering matrix are normalized
    and non-negative.
    """
    scattering_matrix = scattering_weight_matrix(group_structure, atomic_mass)
    row_sum = scattering_matrix.sum(axis=1)
    assert np.logical_or(np.isclose(row_sum, 1, rtol=0), row_sum < 1).all(), (
        "Every row must be unitary or less."
    )
    assert (scattering_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg1():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = n2n_weight_matrix(fibo_gs, -1)
    row_sum = n2n_matrix.sum(axis=1)
    assert np.isclose(
        row_sum, [1, 1, 1, np.log(2 / 1.5) / np.log(2 / 1)], rtol=0
    ).all(), "Expected bin 4 to be partially out-of-bounds"
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg2():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = n2n_weight_matrix(fibo_gs, -2)
    row_sum = n2n_matrix.sum(axis=1)
    assert np.isclose(row_sum, [1, 1, 1, 0], rtol=0).all(), (
        "Expected bin 4 to be out-of-bounds"
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg3():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = n2n_weight_matrix(fibo_gs, -3)
    row_sum = n2n_matrix.sum(axis=1)
    assert np.isclose(
        row_sum, [1, 1, np.log(3 / 2.5) / np.log(3 / 2), 0], rtol=0
    ).all(), "Expected bin 3 to be partially out-of-bounds"
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_pos3():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = n2n_weight_matrix(fibo_gs, 3)
    row_sum = n2n_matrix.sum(axis=1)
    assert np.isclose(
        row_sum, [np.log(6.5 / 5) / np.log(8 / 5), 1, 1, 1], rtol=0
    ).all(), "Expected bin 1 to be partially out-of-bounds"
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_pos6():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = n2n_weight_matrix(fibo_gs, 6)
    row_sum = n2n_matrix.sum(axis=1)
    assert np.isclose(row_sum, [0, 1, 1, 1], rtol=0).all(), (
        "Expected bin 1 to be completely out-of-bounds"
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."
