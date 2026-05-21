import numpy as np
import pytest
from numpy import typing as npt

from process.models.neutronics.data import (
    nXn_weight_matrix,
    scattering_weight_matrix,
    calculate_mean_energy_and_incident_bin,
    NEUTRON_MASS_AS_ENERGY,
    relativistic_energy_to_speed,
    SPEED_OF_LIGHT,
)

rng = np.random.default_rng(1)

def fraction_of_speed_of_light(ratio_of_ke_to_re):
    """
    Calculate the speed of the object (as a fraction of speed of light) given its kinetic
    energy: rest energy ratio.
    """
    return relativistic_energy_to_speed(ratio_of_ke_to_re * NEUTRON_MASS_AS_ENERGY)/SPEED_OF_LIGHT

@pytest.mark.parametrize(
    "kinetic_energy_to_rest_energy_ratio, speed_of_light_fraction_expected",
    [
        (0.0, 0.0),
        (1.0, 0.86603),
        (2.0, 0.94281),
        (3.0, 0.96825),
    ],
)
def test_relativistic_energy_to_speed(
        kinetic_energy_to_rest_energy_ratio: float,
        speed_of_light_fraction_expected: float,
    ):
    assert np.isclose(
        fraction_of_speed_of_light(kinetic_energy_to_rest_energy_ratio),
        speed_of_light_fraction_expected,
        atol=0.0, rtol=0.0001,
    )


@pytest.mark.parametrize(
    "group_structure, atomic_mass",
    [
        (np.cumsum(100 ** abs(rng.normal(size=5)))[::-1], rng.random() * 200),
        (np.cumsum(100 ** abs(rng.normal(size=10)))[::-1], rng.random() * 239),
        (np.cumsum(100 ** abs(rng.normal(size=20)))[::-1], rng.random() * 150),
        (np.cumsum(100 ** abs(rng.normal(size=30)))[::-1], 1),
        (np.cumsum(100 ** abs(rng.normal(size=40)))[::-1], 2),
    ],
)
def test_scattering_matrix(group_structure: npt.NDArray[np.float64], atomic_mass: float):
    """
    Check that even randomly generated scattering matrix are normalisize=zed
    and non-negative.
    """
    incident_energy = np.mean(sorted(group_structure)[-2:])
    group_energy, _ = calculate_mean_energy_and_incident_bin(group_structure, incident_energy)
    scattering_matrix = scattering_weight_matrix(group_structure, group_energy, atomic_mass)
    row_sum = scattering_matrix.sum(axis=1)
    assert np.logical_or(np.isclose(row_sum, 1, rtol=0), row_sum < 1).all(), (
        "Every row must be unitary or less, given no up-scatter."
    )
    assert (scattering_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg1():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = nXn_weight_matrix(fibo_gs, -1, 2)
    row_sum = n2n_matrix.sum(axis=1)
    np.testing.assert_allclose(
        row_sum,
        [1, 1, 1, np.log(2 / 1.5) / np.log(2 / 1)],
        err_msg="Expected bin 4 to be partially out-of-bounds",
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg2():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = nXn_weight_matrix(fibo_gs, -2, 2)
    row_sum = n2n_matrix.sum(axis=1)
    np.testing.assert_allclose(
        row_sum, [1, 1, 1, 0], err_msg="Expected bin 4 to be out-of-bounds"
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_neg3():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = nXn_weight_matrix(fibo_gs, -3, 2)
    row_sum = n2n_matrix.sum(axis=1)
    np.testing.assert_allclose(
        row_sum,
        [1, 1, np.log(3 / 2.5) / np.log(3 / 2), 0],
        err_msg="Expected bin 3 to be partially out-of-bounds",
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_pos3():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = nXn_weight_matrix(fibo_gs, 3, 2)
    row_sum = n2n_matrix.sum(axis=1)
    np.testing.assert_allclose(
        row_sum,
        [np.log(6.5 / 5) / np.log(8 / 5), 1, 1, 1],
        err_msg="Expected bin 1 to be partially out-of-bounds",
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."


def test_n2n_matrix_pos6():
    fibo_gs = np.array([8, 5, 3, 2, 1])
    n2n_matrix = nXn_weight_matrix(fibo_gs, 6, 2)
    row_sum = n2n_matrix.sum(axis=1)
    np.testing.assert_allclose(
        row_sum,
        [0, 1, 1, 1],
        err_msg="Expected bin 1 to be completely out-of-bounds",
    )
    assert (n2n_matrix >= 0).all(), "Must be all non-negative."
