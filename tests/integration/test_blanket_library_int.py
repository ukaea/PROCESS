import pytest

from process.blanket_library import BlanketLibrary
from process.fortran import (
    build_variables as bv,
)
from process.fortran import (
    fwbs_variables as fwbs,
)
from process.fortran import (
    physics_variables as pv,
)
from process.fw import Fw


@pytest.fixture
def blanket_library_fixture():
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return BlanketLibrary(Fw())


def test_hydraulic_diameter(monkeypatch, blanket_library_fixture):
    """
    Test for hydraulic_diameter function.
    """
    # Set var values
    monkeypatch.setattr(fwbs, "afw", 1.0)
    monkeypatch.setattr(fwbs, "a_bz_liq", 1.0)
    monkeypatch.setattr(fwbs, "b_bz_liq", 1.0)

    # hydraulic_diameter input = i_channel_shape: 1 = circle, 2 = rectangle
    assert blanket_library_fixture.hydraulic_diameter(1) == 2.0  # 2.0D0*afw
    assert (
        blanket_library_fixture.hydraulic_diameter(2) == 1.0
    )  # 2*a_bz_liq*b_bz_liq/(a_bz_liq+b_bz_liq)


def test_elbow_coeff(blanket_library_fixture):
    """
    Test for elbow_coeff function.
    """
    # input = r_elbow, ang_elbow, lambda, dh
    assert blanket_library_fixture.elbow_coeff(1, 0, 1, 1) == pytest.approx(
        0.0, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 1, 1) == pytest.approx(
        1.785, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 180, 1, 1) == pytest.approx(
        3.3, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 1, 0.1) == pytest.approx(
        15.816, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(0.1, 90, 1, 1) == pytest.approx(
        66.57, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 0.1, 1) == pytest.approx(
        0.3675, rel=1e-3
    )


def test_flow_velocity(monkeypatch, blanket_library_fixture):
    """
    Test for flow_velocity function.
    """
    # Set var values
    monkeypatch.setattr(fwbs, "afw", 1.0)
    monkeypatch.setattr(fwbs, "a_bz_liq", 1.0)
    monkeypatch.setattr(fwbs, "b_bz_liq", 1.0)

    # input = i_channel_shape, mass_flow_rate, flow_density
    assert blanket_library_fixture.flow_velocity(1, 1, 1) == pytest.approx(
        0.318, rel=1e-3
    )
    assert blanket_library_fixture.flow_velocity(2, 1, 1) == 1.0
    assert blanket_library_fixture.flow_velocity(1, 0, 1) == 0.0
    assert blanket_library_fixture.flow_velocity(2, 0, 1) == 0.0


def test_liquid_breeder_properties_part_1(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure.
    PbLi or Li, with inboard blanket, no inlet/outlet temp difference.
    """
    # Set var values
    monkeypatch.setattr(fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(pv, "bt", 6.0)
    monkeypatch.setattr(pv, "rmajor", 8.0)
    monkeypatch.setattr(pv, "aspect", 3.0)
    monkeypatch.setattr(bv, "blnkith", 0.1)
    monkeypatch.setattr(bv, "blnkoth", 0.2)
    monkeypatch.setattr(fwbs, "iblnkith", 1)
    monkeypatch.setattr(fwbs, "inlet_temp_liq", 1.0)
    monkeypatch.setattr(fwbs, "outlet_temp_liq", 1.0)

    # PbLi - see [Fer2020] for relavent equations
    monkeypatch.setattr(fwbs, "i_bb_liq", 0)

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(fwbs.den_liq, rel=1e-3) == 1.052e4
    assert pytest.approx(fwbs.specific_heat_liq, rel=1e-3) == 195.0
    assert pytest.approx(fwbs.thermal_conductivity_liq, rel=1e-3) == -3.384
    assert pytest.approx(fwbs.dynamic_viscosity_liq, rel=1e-3) == 0.0155
    assert pytest.approx(fwbs.electrical_conductivity_liq, rel=1e-3) == 9.71e5

    assert pytest.approx(fwbs.b_mag_blkt, rel=1e-3) == (9.085, 4.458)
    assert pytest.approx(fwbs.hartmann_liq, rel=1e-3) == (7.189e3, 3.528e3)

    # Li - see [Lyublinski et al., 2009] for relavent equations
    monkeypatch.setattr(fwbs, "i_bb_liq", 1)

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(fwbs.den_liq, rel=1e-3) == 504.0
    assert pytest.approx(fwbs.specific_heat_liq, rel=1e-3) == 2.833e6
    assert pytest.approx(fwbs.dynamic_viscosity_liq, rel=1e-3) == 1.051e112
    assert pytest.approx(fwbs.electrical_conductivity_liq, rel=1e-3) == 9.27e8
    assert pytest.approx(fwbs.b_mag_blkt, rel=1e-3) == (9.085, 4.458)
    assert pytest.approx(fwbs.hartmann_liq, rel=1e-3) == (2.7e-53, 1.3e-53)

    # con_vsc_rat = electrical_conductivity_liq/dynamic_viscosity_liq
    # hartmann_liq = b_mag_blkt * a_bz_liq/2.0D0 * sqrt(con_vsc_rat)


def test_liquid_breeder_properties_part_2(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure. No inboard blanket.
    """
    # Set var values
    monkeypatch.setattr(fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(pv, "bt", 6.0)
    monkeypatch.setattr(pv, "rmajor", 8.0)
    monkeypatch.setattr(pv, "aspect", 3.0)
    monkeypatch.setattr(bv, "blnkith", 0.0)
    monkeypatch.setattr(bv, "blnkoth", 0.2)
    monkeypatch.setattr(fwbs, "iblnkith", 0)
    monkeypatch.setattr(fwbs, "i_bb_liq", 0)
    monkeypatch.setattr(fwbs, "inlet_temp_liq", 0.0)
    monkeypatch.setattr(fwbs, "outlet_temp_liq", 0.0)

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(fwbs.b_mag_blkt, rel=1e-3) == (8.999, 4.458)


def test_liquid_breeder_properties_part_3(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure.
    With inlet/outlet temp difference.
    """
    # Set var values
    monkeypatch.setattr(fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(pv, "bt", 6.0)
    monkeypatch.setattr(pv, "rmajor", 8.0)
    monkeypatch.setattr(pv, "aspect", 3.0)
    monkeypatch.setattr(bv, "blnkith", 0.1)
    monkeypatch.setattr(bv, "blnkoth", 0.2)
    monkeypatch.setattr(fwbs, "iblnkith", 1)
    monkeypatch.setattr(fwbs, "inlet_temp_liq", 0.0)
    monkeypatch.setattr(fwbs, "outlet_temp_liq", 1.0)

    # PbLi - see [Fer2020] for relavent equations
    monkeypatch.setattr(fwbs, "i_bb_liq", 0)

    blanket_library_fixture.liquid_breeder_properties()
    assert pytest.approx(fwbs.den_liq, rel=1e-3) == 1.052e4
    # Li - see [Lyublinski et al., 2009] for relavent equations
    monkeypatch.setattr(fwbs, "i_bb_liq", 1)

    blanket_library_fixture.liquid_breeder_properties()
    assert pytest.approx(fwbs.den_liq, rel=1e-3) == 504.0


def test_pressure_drop(monkeypatch, blanket_library_fixture):
    """
    Test for pressure_drop function.
    """

    # Set var values
    monkeypatch.setattr(fwbs, "afw", 1.0)
    monkeypatch.setattr(fwbs, "a_bz_liq", 1.0)
    monkeypatch.setattr(fwbs, "b_bz_liq", 1.0)
    monkeypatch.setattr(fwbs, "roughness", 1.0e-6)

    # input = ip, ofile, i_ps, num_90, num_180, l_pipe, den, vsc, vv, label
    assert blanket_library_fixture.pressure_drop(
        2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, "label"
    ) == pytest.approx(1.438, rel=1e-3)


# Should add test_liquid_breeder_pressure_drop_mhd
