from typing import Any, NamedTuple

import pytest

from process.models.engineering.pumping import (
    CoolantFrictionLossParameters,
    CoolantType,
    calculate_required_mass_flow_rate,
    calculate_reynolds_number,
    coolant_friction_pressure_drop,
    coolant_pumping_power,
    darcy_friction_haaland,
    elbow_coeff,
    gnielinski_heat_transfer_coefficient,
    pipe_hydraulic_diameter,
)


@pytest.fixture
def blanket_library(process_models):
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return process_models.blanket_library


def test_darcy_friction_haaland():
    assert darcy_friction_haaland(
        reynolds=5500, roughness_channel=1e-6, radius_channel=0.1
    ) == pytest.approx(0.0366668931278784)


def test_gnielinski_heat_transfer_coefficient():
    assert gnielinski_heat_transfer_coefficient(
        mflux_coolant=112.19853108876258,
        den_coolant=8.8673250601290707,
        radius_channel=0.0060000000000000001,
        heatcap_coolant=5184.9330299967578,
        visc_coolant=4.0416219836935569e-05,
        thermcond_coolant=0.3211653052986152,
        roughness_channel=6e-8,
    ) == pytest.approx(1929.2042015869506)


def test_calculate_reynolds_number():

    assert calculate_reynolds_number(
        den_coolant=8.8673250601290707,
        vel_coolant=12.649110769896881,
        radius_channel=0.0060000000000000001,
        visc_coolant=4.0416219836935569e-05,
    ) == pytest.approx(33302.602975971815)


def test_elbow_coeff():
    """
    Test for elbow_coeff function.
    """
    # input = r_elbow, ang_elbow, lambda, dh
    assert elbow_coeff(1, 0, 1, 1) == pytest.approx(0.0, rel=1e-3)
    assert elbow_coeff(1, 90, 1, 1) == pytest.approx(1.7807963267948965, rel=1e-3)
    assert elbow_coeff(1, 180, 1, 1) == pytest.approx(3.291157766597427, rel=1e-3)
    assert elbow_coeff(1, 90, 1, 0.1) == pytest.approx(15.774371098812502, rel=1e-3)
    assert elbow_coeff(0.1, 90, 1, 1) == pytest.approx(66.57, rel=1e-3)
    assert elbow_coeff(1, 90, 0.1, 1) == pytest.approx(0.3670796326794896, rel=1e-3)


def test_calculate_required_mass_flow_rate():
    assert calculate_required_mass_flow_rate(
        p_heat_total=1000.0,
        heatcap_coolant=100.0,
        temp_in_coolant=300.0,
        temp_out_coolant=310.0,
    ) == pytest.approx(1.0)


def test_calculate_required_mass_flow_rate_with_realistic_values():
    assert calculate_required_mass_flow_rate(
        p_heat_total=50000.0,
        heatcap_coolant=4180.0,
        temp_in_coolant=293.15,
        temp_out_coolant=313.15,
    ) == pytest.approx(0.5980861244019139)


def test_calculate_required_mass_flow_rate_zero_temperature_rise():
    with pytest.raises(ZeroDivisionError):
        calculate_required_mass_flow_rate(
            p_heat_total=1000.0,
            heatcap_coolant=4180.0,
            temp_in_coolant=300.0,
            temp_out_coolant=300.0,
        )


def test_hydraulic_diameter():
    """
    Test for hydraulic_diameter function.
    """

    # hydraulic_diameter input = i_channel_shape: 1 = circle, 2 = rectangle
    # 2.0D0*radius_fw_channel
    assert pipe_hydraulic_diameter(
        1, radius_fw_channel=1.0, a_bz_liq=1.0, b_bz_liq=1.0
    ) == pytest.approx(2.0)
    # 2*a_bz_liq*b_bz_liq/(a_bz_liq+b_bz_liq)
    assert pipe_hydraulic_diameter(
        2, radius_fw_channel=1.0, a_bz_liq=1.0, b_bz_liq=1.0
    ) == pytest.approx(1.0)


class CoolantFrictionLossParam(NamedTuple):
    radius_channel: Any = None
    radius_pipe_90_deg_bend: Any = None
    radius_pipe_180_deg_bend: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    roughness_channel: Any = None
    i_ps: Any = None
    n_pipe_90_deg_bends: Any = None
    n_pipe_180_deg_bends: Any = None
    len_pipe: Any = None
    den_coolant: Any = None
    visc_coolant: Any = None
    vel_coolant: Any = None
    label: Any = None
    expected_pressure_drop_out: Any = None


@pytest.mark.parametrize(
    "coolantfrictionlossparam",
    [
        CoolantFrictionLossParam(
            radius_channel=0.0060000000000000001,
            radius_pipe_90_deg_bend=0.018,
            radius_pipe_180_deg_bend=0.09,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            roughness_channel=9.9999999999999995e-07,
            i_ps=1,
            n_pipe_90_deg_bends=2,
            n_pipe_180_deg_bends=0,
            len_pipe=4,
            den_coolant=10.405276820718059,
            visc_coolant=3.604452999475736e-05,
            vel_coolant=32.753134225223164,
            label="Inboard first wall",
            expected_pressure_drop_out=36213.58989742931,
        ),
        CoolantFrictionLossParam(
            radius_channel=1.0,
            radius_pipe_90_deg_bend=1.0,
            radius_pipe_180_deg_bend=1.0,
            a_bz_liq=1.0,
            b_bz_liq=1.0,
            roughness_channel=1e-6,
            i_ps=2,
            n_pipe_90_deg_bends=1.0,
            n_pipe_180_deg_bends=1.0,
            len_pipe=1.0,
            den_coolant=1.0,
            visc_coolant=1.0,
            vel_coolant=1.0,
            label="label",
            expected_pressure_drop_out=1.4325633520224854,
        ),
    ],
)
def test_coolant_friction_loss(coolantfrictionlossparam, monkeypatch, blanket_library):
    """
    Automatically generated Regression Unit Test for coolant_friction_loss.

    This test was generated using data from
    blanket_files/large_tokamak_primary_pumping2.IN.DAT.

    Parameters
    ----------
    coolantfrictionlossparam : CoolantFrictionLossParam
        the data used to mock and assert in this test.
    monkeypatch : _pytest.monkeypatch.monkeypatch
        pytest fixture used to mock module/class variables
    blanket_library : BlanketLibrary
        the blanket library instance used in this test.

    """
    monkeypatch.setattr(
        blanket_library.data.fwbs, "a_bz_liq", coolantfrictionlossparam.a_bz_liq
    )
    monkeypatch.setattr(
        blanket_library.data.fwbs, "b_bz_liq", coolantfrictionlossparam.b_bz_liq
    )

    pressure_params: CoolantFrictionLossParameters = coolant_friction_pressure_drop(
        i_ps=coolantfrictionlossparam.i_ps,
        radius_pipe_90_deg_bend=(coolantfrictionlossparam.radius_pipe_90_deg_bend),
        radius_pipe_180_deg_bend=(coolantfrictionlossparam.radius_pipe_180_deg_bend),
        n_pipe_90_deg_bends=coolantfrictionlossparam.n_pipe_90_deg_bends,
        n_pipe_180_deg_bends=coolantfrictionlossparam.n_pipe_180_deg_bends,
        len_pipe=coolantfrictionlossparam.len_pipe,
        den_coolant=coolantfrictionlossparam.den_coolant,
        visc_coolant=coolantfrictionlossparam.visc_coolant,
        vel_coolant=coolantfrictionlossparam.vel_coolant,
        roughness_channel=coolantfrictionlossparam.roughness_channel,
        radius_channel=coolantfrictionlossparam.radius_channel,
        a_bz_liq=coolantfrictionlossparam.a_bz_liq,
        b_bz_liq=coolantfrictionlossparam.b_bz_liq,
    )

    assert pressure_params.dpres_total == pytest.approx(
        coolantfrictionlossparam.expected_pressure_drop_out
    )


def test_pumppower_primary_helium():

    data = {
        "i_liquid_breeder": 2,
        "temp_coolant_pump_outlet": 570,
        "temp_coolant_pump_inlet": 720,
        "pres_coolant_pump_inlet": 1700000,
        "dpres_coolant": 303517.3,
        "mflow_coolant_total": 35677.7,
        "i_coolant_type": CoolantType.HELIUM,
        "den_coolant": 9753.25,
        "etaiso": 0.9,
        "etaiso_liq": 0.85,
    }

    assert pytest.approx(coolant_pumping_power(**data)) == 1.8251284651310427


def test_pumppower_secondary_pb_li():

    data = {
        "i_liquid_breeder": 1,
        "temp_coolant_pump_outlet": 573,
        "temp_coolant_pump_inlet": 773,
        "pres_coolant_pump_inlet": 8000000,
        "dpres_coolant": 20088.23,
        "mflow_coolant_total": 956.3,
        "i_coolant_type": CoolantType.HELIUM,
        "den_coolant": 5.64,
        "etaiso": 0.9,
        "etaiso_liq": 0.85,
    }

    assert pytest.approx(coolant_pumping_power(**data), rel=1e-4) == 3.2374845432302464
