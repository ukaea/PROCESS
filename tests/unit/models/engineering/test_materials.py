import numpy as np
import pytest

from process.models.engineering.materials import (
    calculate_tresca_stress,
    calculate_von_mises_stress,
    eurofer97_thermal_conductivity,
)


def test_eurofer97_thermal_conductivity(monkeypatch, process_models):
    monkeypatch.setattr(process_models.data.fwbs, "fw_th_conductivity", 28.9)

    assert eurofer97_thermal_conductivity(
        1900.0, process_models.data.fwbs.fw_th_conductivity
    ) == pytest.approx(326.70406785462256)


def test_calculate_tresca_stress_floats():
    """Test Tresca stress calculation with float inputs."""
    result = calculate_tresca_stress(100.0, 50.0, 25.0)
    assert result == pytest.approx(75.0)


def test_calculate_tresca_stress_ndarrays():
    """Test Tresca stress calculation with ndarray inputs."""
    stress_x = np.array([100.0, 200.0])
    stress_y = np.array([50.0, 100.0])
    stress_z = np.array([25.0, 50.0])

    result = calculate_tresca_stress(stress_x, stress_y, stress_z)
    expected = np.array([75.0, 150.0])

    np.testing.assert_array_almost_equal(result, expected)


def test_calculate_tresca_stress_mixed_float_ndarray():
    """Test Tresca stress calculation with mixed float and ndarray inputs."""
    stress_x = 100.0
    stress_y = np.array([50.0, 75.0])
    stress_z = 25.0

    result = calculate_tresca_stress(stress_x, stress_y, stress_z)
    expected = np.array([75.0, 75.0])

    np.testing.assert_array_almost_equal(result, expected)


def test_calculate_tresca_stress_zeros():
    """Test Tresca stress calculation with zero values."""
    result = calculate_tresca_stress(0.0, 0.0, 0.0)
    assert result == pytest.approx(0.0)


def test_calculate_tresca_stress_zeros_ndarray():
    """Test Tresca stress calculation with zero ndarray values."""
    stress_x = np.array([0.0, 100.0])
    stress_y = np.array([0.0, 50.0])
    stress_z = np.array([0.0, 25.0])

    result = calculate_tresca_stress(stress_x, stress_y, stress_z)
    expected = np.array([0.0, 75.0])

    np.testing.assert_array_almost_equal(result, expected)


def test_calculate_von_mises_stress_floats():
    """Test von Mises stress calculation with float inputs."""
    result = calculate_von_mises_stress(100.0, 50.0, 25.0, 10.0, 5.0, 2.0)
    expected = np.sqrt(
        0.5
        * (
            (100.0 - 50.0) ** 2
            + (50.0 - 25.0) ** 2
            + (25.0 - 100.0) ** 2
            + 6 * (10.0**2 + 5.0**2 + 2.0**2)
        )
    )
    assert result == pytest.approx(expected)


def test_calculate_von_mises_stress_ndarrays():
    """Test von Mises stress calculation with ndarray inputs."""
    stress_x = np.array([100.0, 200.0])
    stress_y = np.array([50.0, 100.0])
    stress_z = np.array([25.0, 50.0])
    stress_shear_xy = np.array([10.0, 20.0])
    stress_shear_yz = np.array([5.0, 10.0])
    stress_shear_zx = np.array([2.0, 4.0])

    result = calculate_von_mises_stress(
        stress_x, stress_y, stress_z, stress_shear_xy, stress_shear_yz, stress_shear_zx
    )
    expected = np.sqrt(
        0.5
        * (
            (stress_x - stress_y) ** 2
            + (stress_y - stress_z) ** 2
            + (stress_z - stress_x) ** 2
            + 6 * (stress_shear_xy**2 + stress_shear_yz**2 + stress_shear_zx**2)
        )
    )

    np.testing.assert_array_almost_equal(result, expected)


def test_calculate_von_mises_stress_mixed_float_ndarray():
    """Test von Mises stress calculation with mixed float and ndarray inputs."""
    stress_x = 100.0
    stress_y = np.array([50.0, 75.0])
    stress_z = 25.0
    stress_shear_xy = 10.0
    stress_shear_yz = np.array([5.0, 7.5])
    stress_shear_zx = 2.0

    result = calculate_von_mises_stress(
        stress_x, stress_y, stress_z, stress_shear_xy, stress_shear_yz, stress_shear_zx
    )
    expected = np.sqrt(
        0.5
        * (
            (stress_x - stress_y) ** 2
            + (stress_y - stress_z) ** 2
            + (stress_z - stress_x) ** 2
            + 6 * (stress_shear_xy**2 + stress_shear_yz**2 + stress_shear_zx**2)
        )
    )

    np.testing.assert_array_almost_equal(result, expected)


def test_calculate_von_mises_stress_zeros():
    """Test von Mises stress calculation with zero values."""
    result = calculate_von_mises_stress(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert result == pytest.approx(0.0)


def test_calculate_von_mises_stress_zeros_ndarray():
    """Test von Mises stress calculation with zero ndarray values."""
    stress_x = np.array([0.0, 100.0])
    stress_y = np.array([0.0, 50.0])
    stress_z = np.array([0.0, 25.0])
    stress_shear_xy = np.array([0.0, 10.0])
    stress_shear_yz = np.array([0.0, 5.0])
    stress_shear_zx = np.array([0.0, 2.0])

    result = calculate_von_mises_stress(
        stress_x, stress_y, stress_z, stress_shear_xy, stress_shear_yz, stress_shear_zx
    )
    expected = np.sqrt(
        0.5
        * (
            (stress_x - stress_y) ** 2
            + (stress_y - stress_z) ** 2
            + (stress_z - stress_x) ** 2
            + 6 * (stress_shear_xy**2 + stress_shear_yz**2 + stress_shear_zx**2)
        )
    )

    np.testing.assert_array_almost_equal(result, expected)
