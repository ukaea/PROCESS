import pytest

from process.core.coolprop_interface import FluidProperties


class WATER_PROPERTIES:
    temperature = 293.15
    pressure = 101325
    density = 998.20609246795
    enthalpy = 84013.058152596
    entropy = 296.48264655725
    specific_heat_const_p = 4184.7940947755
    specific_heat_const_v = 4157.4368513145
    viscosity = 0.0010016053256479
    thermal_conductivity = 0.5984608954958


def test_water_with_pressure_temp():
    fluid_properties = FluidProperties.of("Water", temperature=293.15, pressure=101325)

    # trivial
    assert pytest.approx(fluid_properties.temperature) == WATER_PROPERTIES.temperature
    assert pytest.approx(fluid_properties.pressure) == WATER_PROPERTIES.pressure

    # calculated
    assert pytest.approx(fluid_properties.density, rel=1e-3) == WATER_PROPERTIES.density
    assert (
        pytest.approx(fluid_properties.viscosity, rel=1e-3) == WATER_PROPERTIES.viscosity
    )
    assert (
        pytest.approx(fluid_properties.enthalpy, rel=1e-3) == WATER_PROPERTIES.enthalpy
    )
    assert pytest.approx(fluid_properties.entropy, rel=1e-3) == WATER_PROPERTIES.entropy
    assert (
        pytest.approx(fluid_properties.specific_heat_const_p, rel=1e-3)
        == WATER_PROPERTIES.specific_heat_const_p
    )
    assert (
        pytest.approx(fluid_properties.specific_heat_const_v, rel=1e-3)
        == WATER_PROPERTIES.specific_heat_const_v
    )
    assert (
        pytest.approx(fluid_properties.thermal_conductivity, rel=1e-3)
        == WATER_PROPERTIES.thermal_conductivity
    )


def test_water_enthalpy_from_entropy_pressure():
    fluid_properties = FluidProperties.of(
        "Water", pressure=WATER_PROPERTIES.pressure, entropy=WATER_PROPERTIES.entropy
    )

    # trivial
    assert pytest.approx(fluid_properties.pressure) == WATER_PROPERTIES.pressure
    assert pytest.approx(fluid_properties.entropy) == WATER_PROPERTIES.entropy

    # calculated
    assert pytest.approx(fluid_properties.enthalpy) == WATER_PROPERTIES.enthalpy


def test_find_water_saturation_temperature():
    fluid_properties = FluidProperties.of(
        "Water", pressure=WATER_PROPERTIES.pressure, vapor_quality=0
    )

    assert pytest.approx(fluid_properties.pressure) == WATER_PROPERTIES.pressure
    assert pytest.approx(fluid_properties.temperature, rel=1e-3) == 373.12
