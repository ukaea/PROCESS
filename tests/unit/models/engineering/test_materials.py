import numpy as np
import pytest

from process.models.engineering.materials import (
    calculate_tresca_stress,
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
