from typing import Any, NamedTuple

import pytest

from process.data_structure import (
    build_variables,
    ccfe_hcpb_module,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    physics_variables,
    primary_pumping_variables,
    tfcoil_variables,
)
from process.fw import Fw
from process.hcpb import CCFE_HCPB


@pytest.fixture
def ccfe_hcpb():
    """Provides CCFE_HCPB object for testing.

    :returns: initialised CCFE_HCPB object
    :rtype: process.hcpb.CCFE_HCPB
    """
    return CCFE_HCPB(Fw())


class NuclearHeatingMagnetsParam(NamedTuple):
    dr_fw_inboard: Any = None

    dr_vv_inboard: Any = None

    dr_vv_outboard: Any = None

    dr_fw_outboard: Any = None

    dr_blkt_inboard: Any = None

    dr_blkt_outboard: Any = None

    dr_shld_inboard: Any = None

    dr_shld_outboard: Any = None

    radius_fw_channel: Any = None

    dx_fw_module: Any = None

    den_steel: Any = None

    m_blkt_total: Any = None

    vol_blkt_total: Any = None

    whtshld: Any = None

    vol_shld_total: Any = None

    m_vv: Any = None

    vol_vv: Any = None

    fw_armour_thickness: Any = None

    p_tf_nuclear_heat_mw: Any = None

    f_a_fw_coolant_inboard: Any = None

    f_a_fw_coolant_outboard: Any = None

    p_fusion_total_mw: Any = None

    itart: Any = None

    m_tf_coils_total: Any = None

    whttflgs: Any = None

    verbose: Any = None

    armour_density: Any = None

    fw_density: Any = None

    blanket_density: Any = None

    shield_density: Any = None

    vv_density: Any = None

    x_blanket: Any = None

    x_shield: Any = None

    tfc_nuc_heating: Any = None

    expected_p_tf_nuclear_heat_mw: Any = None

    expected_f_a_fw_coolant_inboard: Any = None

    expected_f_a_fw_coolant_outboard: Any = None

    expected_armour_density: Any = None

    expected_fw_density: Any = None

    expected_blanket_density: Any = None

    expected_shield_density: Any = None

    expected_vv_density: Any = None

    expected_x_blanket: Any = None

    expected_x_shield: Any = None

    expected_tfc_nuc_heating: Any = None


@pytest.mark.parametrize(
    "nuclearheatingmagnetsparam",
    (
        NuclearHeatingMagnetsParam(
            dr_fw_inboard=0.018000000000000002,
            dr_vv_inboard=0.30000000000000004,
            dr_vv_outboard=0.30000000000000004,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            dr_blkt_outboard=0.98199999999999998,
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            radius_fw_channel=0.0060000000000000001,
            dx_fw_module=0.02,
            den_steel=7800,
            m_blkt_total=3501027.3252278985,
            vol_blkt_total=1397.9003011502937,
            whtshld=2294873.8131476045,
            vol_shld_total=735.53647857295027,
            m_vv=9043937.8018644415,
            vol_vv=1159.4792053672361,
            fw_armour_thickness=0.0050000000000000001,
            p_tf_nuclear_heat_mw=0,
            f_a_fw_coolant_inboard=0,
            f_a_fw_coolant_outboard=0,
            p_fusion_total_mw=1986.0623241661431,
            itart=0,
            m_tf_coils_total=19649856.627845347,
            whttflgs=0,
            verbose=0,
            armour_density=0,
            fw_density=0,
            blanket_density=0,
            shield_density=0,
            vv_density=0,
            x_blanket=0,
            x_shield=0,
            tfc_nuc_heating=0,
            expected_p_tf_nuclear_heat_mw=0.044541749095475737,
            expected_f_a_fw_coolant_inboard=0.31415926535897931,
            expected_f_a_fw_coolant_outboard=0.31415926535897931,
            expected_armour_density=13202.434141839649,
            expected_fw_density=5349.557730199961,
            expected_blanket_density=2504.4899999999998,
            expected_shield_density=3119.9999999999995,
            expected_vv_density=7800,
            expected_x_blanket=2.3374537748527975,
            expected_x_shield=4.056,
            expected_tfc_nuc_heating=22427.165831352642,
        ),
        NuclearHeatingMagnetsParam(
            dr_fw_inboard=0.018000000000000002,
            dr_vv_inboard=0.30000000000000004,
            dr_vv_outboard=0.30000000000000004,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            dr_blkt_outboard=0.98199999999999998,
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            radius_fw_channel=0.0060000000000000001,
            dx_fw_module=0.02,
            den_steel=7800,
            m_blkt_total=3507503.3737008357,
            vol_blkt_total=1400.4860764869636,
            whtshld=2297808.3935174854,
            vol_shld_total=736.47704920432227,
            m_vv=9056931.558219457,
            vol_vv=1161.1450715665972,
            fw_armour_thickness=0.0050000000000000001,
            p_tf_nuclear_heat_mw=0.044184461825198453,
            f_a_fw_coolant_inboard=0.31415926535897931,
            f_a_fw_coolant_outboard=0.31415926535897931,
            p_fusion_total_mw=1985.4423932312809,
            itart=0,
            m_tf_coils_total=19662548.210142396,
            whttflgs=0,
            verbose=0,
            armour_density=13202.434141839649,
            fw_density=5349.557730199961,
            blanket_density=2504.4899999999998,
            shield_density=3119.9999999999995,
            vv_density=7800,
            x_blanket=2.3374537748527975,
            x_shield=4.056,
            tfc_nuc_heating=22427.165831352642,
            expected_p_tf_nuclear_heat_mw=0.044556605747797934,
            expected_f_a_fw_coolant_inboard=0.31415926535897931,
            expected_f_a_fw_coolant_outboard=0.31415926535897931,
            expected_armour_density=13202.434141839649,
            expected_fw_density=5349.557730199961,
            expected_blanket_density=2504.4900000000002,
            expected_shield_density=3120,
            expected_vv_density=7799.9999999999991,
            expected_x_blanket=2.3374537748527979,
            expected_x_shield=4.056,
            expected_tfc_nuc_heating=22441.651240901861,
        ),
    ),
)
def test_nuclear_heating_magnets(nuclearheatingmagnetsparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_magnets.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingmagnetsparam: the data used to mock and assert in this test.
    :type nuclearheatingmagnetsparam: nuclearheatingmagnetsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        build_variables, "dr_fw_inboard", nuclearheatingmagnetsparam.dr_fw_inboard
    )

    monkeypatch.setattr(
        build_variables, "dr_vv_inboard", nuclearheatingmagnetsparam.dr_vv_inboard
    )

    monkeypatch.setattr(
        build_variables, "dr_vv_outboard", nuclearheatingmagnetsparam.dr_vv_outboard
    )

    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", nuclearheatingmagnetsparam.dr_fw_outboard
    )

    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", nuclearheatingmagnetsparam.dr_blkt_inboard
    )

    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", nuclearheatingmagnetsparam.dr_blkt_outboard
    )

    monkeypatch.setattr(
        build_variables, "dr_shld_inboard", nuclearheatingmagnetsparam.dr_shld_inboard
    )

    monkeypatch.setattr(
        build_variables, "dr_shld_outboard", nuclearheatingmagnetsparam.dr_shld_outboard
    )

    monkeypatch.setattr(
        fwbs_variables,
        "radius_fw_channel",
        nuclearheatingmagnetsparam.radius_fw_channel,
    )

    monkeypatch.setattr(
        fwbs_variables, "dx_fw_module", nuclearheatingmagnetsparam.dx_fw_module
    )

    monkeypatch.setattr(
        fwbs_variables, "den_steel", nuclearheatingmagnetsparam.den_steel
    )

    monkeypatch.setattr(
        fwbs_variables, "m_blkt_total", nuclearheatingmagnetsparam.m_blkt_total
    )

    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_total", nuclearheatingmagnetsparam.vol_blkt_total
    )

    monkeypatch.setattr(fwbs_variables, "whtshld", nuclearheatingmagnetsparam.whtshld)

    monkeypatch.setattr(
        fwbs_variables, "vol_shld_total", nuclearheatingmagnetsparam.vol_shld_total
    )

    monkeypatch.setattr(fwbs_variables, "m_vv", nuclearheatingmagnetsparam.m_vv)

    monkeypatch.setattr(fwbs_variables, "vol_vv", nuclearheatingmagnetsparam.vol_vv)

    monkeypatch.setattr(
        fwbs_variables,
        "fw_armour_thickness",
        nuclearheatingmagnetsparam.fw_armour_thickness,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_tf_nuclear_heat_mw",
        nuclearheatingmagnetsparam.p_tf_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "f_a_fw_coolant_inboard",
        nuclearheatingmagnetsparam.f_a_fw_coolant_inboard,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "f_a_fw_coolant_outboard",
        nuclearheatingmagnetsparam.f_a_fw_coolant_outboard,
    )

    monkeypatch.setattr(
        physics_variables,
        "p_fusion_total_mw",
        nuclearheatingmagnetsparam.p_fusion_total_mw,
    )

    monkeypatch.setattr(physics_variables, "itart", nuclearheatingmagnetsparam.itart)

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coils_total",
        nuclearheatingmagnetsparam.m_tf_coils_total,
    )

    monkeypatch.setattr(
        tfcoil_variables, "whttflgs", nuclearheatingmagnetsparam.whttflgs
    )

    monkeypatch.setattr(global_variables, "verbose", nuclearheatingmagnetsparam.verbose)

    monkeypatch.setattr(
        ccfe_hcpb_module, "armour_density", nuclearheatingmagnetsparam.armour_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "fw_density", nuclearheatingmagnetsparam.fw_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "blanket_density", nuclearheatingmagnetsparam.blanket_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "shield_density", nuclearheatingmagnetsparam.shield_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "vv_density", nuclearheatingmagnetsparam.vv_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "x_blanket", nuclearheatingmagnetsparam.x_blanket
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "x_shield", nuclearheatingmagnetsparam.x_shield
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "tfc_nuc_heating", nuclearheatingmagnetsparam.tfc_nuc_heating
    )

    ccfe_hcpb.nuclear_heating_magnets(False)

    assert fwbs_variables.p_tf_nuclear_heat_mw == pytest.approx(
        nuclearheatingmagnetsparam.expected_p_tf_nuclear_heat_mw
    )

    assert fwbs_variables.f_a_fw_coolant_inboard == pytest.approx(
        nuclearheatingmagnetsparam.expected_f_a_fw_coolant_inboard
    )

    assert fwbs_variables.f_a_fw_coolant_outboard == pytest.approx(
        nuclearheatingmagnetsparam.expected_f_a_fw_coolant_outboard
    )

    assert ccfe_hcpb_module.armour_density == pytest.approx(
        nuclearheatingmagnetsparam.expected_armour_density
    )

    assert ccfe_hcpb_module.fw_density == pytest.approx(
        nuclearheatingmagnetsparam.expected_fw_density
    )

    assert ccfe_hcpb_module.blanket_density == pytest.approx(
        nuclearheatingmagnetsparam.expected_blanket_density
    )

    assert ccfe_hcpb_module.shield_density == pytest.approx(
        nuclearheatingmagnetsparam.expected_shield_density
    )

    assert ccfe_hcpb_module.vv_density == pytest.approx(
        nuclearheatingmagnetsparam.expected_vv_density
    )

    assert ccfe_hcpb_module.x_blanket == pytest.approx(
        nuclearheatingmagnetsparam.expected_x_blanket
    )

    assert ccfe_hcpb_module.x_shield == pytest.approx(
        nuclearheatingmagnetsparam.expected_x_shield
    )

    assert ccfe_hcpb_module.tfc_nuc_heating == pytest.approx(
        nuclearheatingmagnetsparam.expected_tfc_nuc_heating
    )


class NuclearHeatingFwParam(NamedTuple):
    p_fw_nuclear_heat_total_mw: Any = None

    m_fw_total: Any = None

    p_fusion_total_mw: Any = None

    fw_armour_u_nuc_heating: Any = None

    expected_p_fw_nuclear_heat_total_mw: Any = None

    expected_fw_armour_u_nuc_heating: Any = None


@pytest.mark.parametrize(
    "nuclearheatingfwparam",
    (
        NuclearHeatingFwParam(
            p_fw_nuclear_heat_total_mw=0,
            m_fw_total=224802.80270851994,
            p_fusion_total_mw=1986.0623241661431,
            fw_armour_u_nuc_heating=6.2500000000000005e-07,
            expected_p_fw_nuclear_heat_total_mw=279.04523551646628,
        ),
        NuclearHeatingFwParam(
            p_fw_nuclear_heat_total_mw=276.80690153753221,
            m_fw_total=182115.83467868491,
            p_fusion_total_mw=1985.4423932312809,
            fw_armour_u_nuc_heating=6.2500000000000005e-07,
            expected_p_fw_nuclear_heat_total_mw=225.98781165610032,
        ),
    ),
)
def test_nuclear_heating_fw(nuclearheatingfwparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_fw.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingfwparam: the data used to mock and assert in this test.
    :type nuclearheatingfwparam: nuclearheatingfwparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        fwbs_variables,
        "p_fw_nuclear_heat_total_mw",
        nuclearheatingfwparam.p_fw_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(fwbs_variables, "m_fw_total", nuclearheatingfwparam.m_fw_total)

    monkeypatch.setattr(
        physics_variables, "p_fusion_total_mw", nuclearheatingfwparam.p_fusion_total_mw
    )

    monkeypatch.setattr(
        ccfe_hcpb_module,
        "fw_armour_u_nuc_heating",
        nuclearheatingfwparam.fw_armour_u_nuc_heating,
    )

    p_fw_nuclear_heat_total_mw = ccfe_hcpb.nuclear_heating_fw(
        m_fw_total=nuclearheatingfwparam.m_fw_total,
        fw_armour_u_nuc_heating=nuclearheatingfwparam.fw_armour_u_nuc_heating,
        p_fusion_total_mw=nuclearheatingfwparam.p_fusion_total_mw,
    )

    assert p_fw_nuclear_heat_total_mw == pytest.approx(
        nuclearheatingfwparam.expected_p_fw_nuclear_heat_total_mw
    )


class NuclearHeatingBlanketParam(NamedTuple):
    m_blkt_total: Any = None

    p_blkt_nuclear_heat_total_mw: Any = None

    p_fusion_total_mw: Any = None

    exp_blanket: Any = None

    expected_p_blkt_nuclear_heat_total_mw: Any = None

    expected_exp_blanket: Any = None


@pytest.mark.parametrize(
    "nuclearheatingblanketparam",
    (
        NuclearHeatingBlanketParam(
            m_blkt_total=3501027.3252278985,
            p_blkt_nuclear_heat_total_mw=0,
            p_fusion_total_mw=1986.0623241661431,
            exp_blanket=0,
            expected_p_blkt_nuclear_heat_total_mw=1517.0907688379014,
            expected_exp_blanket=0.99982809071915879,
        ),
        NuclearHeatingBlanketParam(
            m_blkt_total=3507503.3737008357,
            p_blkt_nuclear_heat_total_mw=1504.9215740808861,
            p_fusion_total_mw=1985.4423932312809,
            exp_blanket=0.99982809071915879,
            expected_p_blkt_nuclear_heat_total_mw=1516.6213709741428,
            expected_exp_blanket=0.99983082524994527,
        ),
    ),
)
def test_nuclear_heating_blanket(nuclearheatingblanketparam, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_blanket.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingblanketparam: the data used to mock and assert in this test.
    :type nuclearheatingblanketparam: nuclearheatingblanketparam

    """

    p_blkt_nuclear_heat_total_mw, exp_blanket = ccfe_hcpb.nuclear_heating_blanket(
        m_blkt_total=nuclearheatingblanketparam.m_blkt_total,
        p_fusion_total_mw=nuclearheatingblanketparam.p_fusion_total_mw,
    )

    assert p_blkt_nuclear_heat_total_mw == pytest.approx(
        nuclearheatingblanketparam.expected_p_blkt_nuclear_heat_total_mw
    )

    assert exp_blanket == pytest.approx(nuclearheatingblanketparam.expected_exp_blanket)


class NuclearHeatingShieldParam(NamedTuple):
    dr_shld_inboard: Any = None

    dr_shld_outboard: Any = None

    whtshld: Any = None

    p_shld_nuclear_heat_mw: Any = None

    p_fusion_total_mw: Any = None

    itart: Any = None

    shield_density: Any = None

    x_blanket: Any = None

    shld_u_nuc_heating: Any = None

    exp_shield1: Any = None

    exp_shield2: Any = None

    expected_p_shld_nuclear_heat_mw: Any = None

    expected_shld_u_nuc_heating: Any = None

    expected_exp_shield1: Any = None

    expected_exp_shield2: Any = None


@pytest.mark.parametrize(
    "nuclearheatingshieldparam",
    (
        NuclearHeatingShieldParam(
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            whtshld=2294873.8131476045,
            p_shld_nuclear_heat_mw=0,
            p_fusion_total_mw=1986.0623241661431,
            itart=0,
            shield_density=3119.9999999999995,
            x_blanket=2.3374537748527975,
            shld_u_nuc_heating=0,
            exp_shield1=0,
            exp_shield2=0,
            expected_p_shld_nuclear_heat_mw=1.3721323841005297,
            expected_shld_u_nuc_heating=690880.82856444374,
            expected_exp_shield1=0.0017209365527675318,
            expected_exp_shield2=0.25426760591013942,
        ),
        NuclearHeatingShieldParam(
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            whtshld=2297808.3935174854,
            p_shld_nuclear_heat_mw=1.3611259588044891,
            p_fusion_total_mw=1985.4423932312809,
            itart=0,
            shield_density=3120,
            x_blanket=2.3374537748527979,
            shld_u_nuc_heating=690880.82856444374,
            exp_shield1=0.0017209365527675318,
            exp_shield2=0.25426760591013942,
            expected_p_shld_nuclear_heat_mw=1.3734581585671393,
            expected_shld_u_nuc_heating=691764.29557941214,
            expected_exp_shield1=0.0017209365527675303,
            expected_exp_shield2=0.25426760591013942,
        ),
    ),
)
def test_nuclear_heating_shield(nuclearheatingshieldparam, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_shield.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingshieldparam: the data used to mock and assert in this test.
    :type nuclearheatingshieldparam: nuclearheatingshieldparam

    """

    p_shld_nuclear_heat_mw, exp_shield1, exp_shield2, shld_u_nuc_heating = (
        ccfe_hcpb.nuclear_heating_shield(
            itart=nuclearheatingshieldparam.itart,
            dr_shld_inboard=nuclearheatingshieldparam.dr_shld_inboard,
            dr_shld_outboard=nuclearheatingshieldparam.dr_shld_outboard,
            shield_density=nuclearheatingshieldparam.shield_density,
            whtshld=nuclearheatingshieldparam.whtshld,
            x_blanket=nuclearheatingshieldparam.x_blanket,
            p_fusion_total_mw=nuclearheatingshieldparam.p_fusion_total_mw,
        )
    )

    assert p_shld_nuclear_heat_mw == pytest.approx(
        nuclearheatingshieldparam.expected_p_shld_nuclear_heat_mw
    )

    assert shld_u_nuc_heating == pytest.approx(
        nuclearheatingshieldparam.expected_shld_u_nuc_heating
    )

    assert exp_shield1 == pytest.approx(nuclearheatingshieldparam.expected_exp_shield1)

    assert exp_shield2 == pytest.approx(nuclearheatingshieldparam.expected_exp_shield2)


class PowerflowCalcParam(NamedTuple):
    a_fw_outboard: Any = None

    a_fw_total: Any = None

    p_beam_orbit_loss_mw: Any = None

    f_ster_div_single: Any = None

    p_div_rad_total_mw: Any = None

    p_fw_hcd_rad_total_mw: Any = None

    f_a_fw_outboard_hcd: Any = None

    p_fw_rad_total_mw: Any = None

    i_blkt_coolant_type: Any = None

    temp_blkt_coolant_out: Any = None

    pres_blkt_coolant: Any = None

    i_p_coolant_pumping: Any = None

    p_fw_nuclear_heat_total_mw: Any = None

    p_blkt_nuclear_heat_total_mw: Any = None

    p_div_nuclear_heat_total_mw: Any = None

    p_shld_nuclear_heat_mw: Any = None

    etaiso: Any = None

    p_cp_shield_nuclear_heat_mw: Any = None

    psurffwi: Any = None

    psurffwo: Any = None

    p_fw_coolant_pump_mw: Any = None

    f_p_fw_coolant_pump_total_heat: Any = None

    p_blkt_coolant_pump_mw: Any = None

    f_p_blkt_coolant_pump_total_heat: Any = None

    p_shld_coolant_pump_mw: Any = None

    f_p_shld_coolant_pump_total_heat: Any = None

    p_div_coolant_pump_mw: Any = None

    f_p_div_coolant_pump_total_heat: Any = None

    n_divertors: Any = None

    p_plasma_rad_mw: Any = None

    p_fw_alpha_mw: Any = None

    p_plasma_separatrix_mw: Any = None

    p_he: Any = None

    dp_he: Any = None

    gamma_he: Any = None

    t_in_bb: Any = None

    t_out_bb: Any = None

    p_fw_blkt_coolant_pump_mw: Any = None

    expected_p_div_rad_total_mw: Any = None

    expected_p_fw_rad_total_mw: Any = None

    expected_psurffwi: Any = None

    expected_psurffwo: Any = None

    expected_p_shld_coolant_pump_mw: Any = None

    expected_p_div_coolant_pump_mw: Any = None

    expected_p_fw_blkt_coolant_pump_mw: Any = None


@pytest.mark.parametrize(
    "powerflowcalcparam",
    (
        PowerflowCalcParam(
            a_fw_outboard=988.92586580655245,
            a_fw_total=1601.1595634509963,
            p_beam_orbit_loss_mw=0,
            f_ster_div_single=0.115,
            p_fw_hcd_rad_total_mw=0,
            f_a_fw_outboard_hcd=0,
            p_fw_rad_total_mw=0,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_out=823,
            pres_blkt_coolant=15500000,
            i_p_coolant_pumping=3,
            p_div_rad_total_mw=33.056596978820579,
            p_fw_nuclear_heat_total_mw=276.80690153753221,
            p_blkt_nuclear_heat_total_mw=1504.9215740808861,
            p_div_nuclear_heat_total_mw=182.71773382328519,
            p_shld_nuclear_heat_mw=1.3611259588044891,
            etaiso=0.90000000000000002,
            p_cp_shield_nuclear_heat_mw=0,
            psurffwi=0,
            psurffwo=0,
            p_fw_coolant_pump_mw=0,
            f_p_fw_coolant_pump_total_heat=0.0050000000000000001,
            p_blkt_coolant_pump_mw=0,
            f_p_blkt_coolant_pump_total_heat=0.0050000000000000001,
            p_shld_coolant_pump_mw=0,
            f_p_shld_coolant_pump_total_heat=0.0050000000000000001,
            p_div_coolant_pump_mw=0,
            f_p_div_coolant_pump_total_heat=0.0050000000000000001,
            n_divertors=1,
            p_plasma_rad_mw=287.44866938104849,
            p_fw_alpha_mw=19.835845058655043,
            p_plasma_separatrix_mw=143.6315222649435,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            p_fw_blkt_coolant_pump_mw=0,
            expected_p_fw_rad_total_mw=254.39207240222791,
            expected_psurffwi=97.271629070225231,
            expected_psurffwo=176.95628839065773,
            expected_p_shld_coolant_pump_mw=0.0068056297940224456,
            expected_p_div_coolant_pump_mw=1.7970292653352464,
            expected_p_fw_blkt_coolant_pump_mw=202.00455086503842,
        ),
        PowerflowCalcParam(
            a_fw_outboard=1168.1172772224481,
            a_fw_total=1891.2865102700493,
            p_beam_orbit_loss_mw=0,
            f_ster_div_single=0.115,
            p_fw_hcd_rad_total_mw=0,
            f_a_fw_outboard_hcd=0,
            p_fw_rad_total_mw=254.39207240222791,
            p_div_rad_total_mw=33.056596978820579,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_out=823,
            pres_blkt_coolant=15500000,
            i_p_coolant_pumping=3,
            p_fw_nuclear_heat_total_mw=230.98304919926957,
            p_blkt_nuclear_heat_total_mw=1550.1447895848396,
            p_div_nuclear_heat_total_mw=182.66070017727785,
            p_shld_nuclear_heat_mw=1.4038170956592293,
            etaiso=0.90000000000000002,
            p_cp_shield_nuclear_heat_mw=0,
            psurffwi=97.271629070225231,
            psurffwo=176.95628839065773,
            p_fw_coolant_pump_mw=0,
            f_p_fw_coolant_pump_total_heat=0.0050000000000000001,
            p_blkt_coolant_pump_mw=0,
            f_p_blkt_coolant_pump_total_heat=0.0050000000000000001,
            p_shld_coolant_pump_mw=0.0068056297940224456,
            f_p_shld_coolant_pump_total_heat=0.0050000000000000001,
            p_div_coolant_pump_mw=1.7970292653352464,
            f_p_div_coolant_pump_total_heat=0.0050000000000000001,
            n_divertors=1,
            p_plasma_rad_mw=287.44866938104849,
            p_fw_alpha_mw=19.829653483586444,
            p_plasma_separatrix_mw=143.51338080047339,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            p_fw_blkt_coolant_pump_mw=202.00455086503842,
            expected_p_fw_rad_total_mw=254.39207240222791,
            expected_psurffwi=97.271629070225259,
            expected_psurffwo=176.95009681558912,
            expected_p_shld_coolant_pump_mw=0.007019085478296147,
            expected_p_div_coolant_pump_mw=1.7961533897828594,
            expected_p_fw_blkt_coolant_pump_mw=201.94492795635171,
        ),
    ),
)
def test_powerflow_calc(powerflowcalcparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for powerflow_calc.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param powerflowcalcparam: the data used to mock and assert in this test.
    :type powerflowcalcparam: powerflowcalcparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        build_variables, "a_fw_outboard", powerflowcalcparam.a_fw_outboard
    )

    monkeypatch.setattr(build_variables, "a_fw_total", powerflowcalcparam.a_fw_total)

    monkeypatch.setattr(
        current_drive_variables,
        "p_beam_orbit_loss_mw",
        powerflowcalcparam.p_beam_orbit_loss_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "f_ster_div_single", powerflowcalcparam.f_ster_div_single
    )

    monkeypatch.setattr(
        fwbs_variables, "p_div_rad_total_mw", powerflowcalcparam.p_div_rad_total_mw
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_fw_hcd_rad_total_mw",
        powerflowcalcparam.p_fw_hcd_rad_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "f_a_fw_outboard_hcd", powerflowcalcparam.f_a_fw_outboard_hcd
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_rad_total_mw", powerflowcalcparam.p_fw_rad_total_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "i_blkt_coolant_type", powerflowcalcparam.i_blkt_coolant_type
    )

    monkeypatch.setattr(
        fwbs_variables,
        "temp_blkt_coolant_out",
        powerflowcalcparam.temp_blkt_coolant_out,
    )

    monkeypatch.setattr(
        fwbs_variables, "pres_blkt_coolant", powerflowcalcparam.pres_blkt_coolant
    )

    monkeypatch.setattr(
        fwbs_variables, "i_p_coolant_pumping", powerflowcalcparam.i_p_coolant_pumping
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_fw_nuclear_heat_total_mw",
        powerflowcalcparam.p_fw_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_blkt_nuclear_heat_total_mw",
        powerflowcalcparam.p_blkt_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_div_nuclear_heat_total_mw",
        powerflowcalcparam.p_div_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_shld_nuclear_heat_mw",
        powerflowcalcparam.p_shld_nuclear_heat_mw,
    )

    monkeypatch.setattr(fwbs_variables, "etaiso", powerflowcalcparam.etaiso)

    monkeypatch.setattr(
        fwbs_variables,
        "p_cp_shield_nuclear_heat_mw",
        powerflowcalcparam.p_cp_shield_nuclear_heat_mw,
    )

    monkeypatch.setattr(fwbs_variables, "psurffwi", powerflowcalcparam.psurffwi)

    monkeypatch.setattr(fwbs_variables, "psurffwo", powerflowcalcparam.psurffwo)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_fw_coolant_pump_mw",
        powerflowcalcparam.p_fw_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_fw_coolant_pump_total_heat",
        powerflowcalcparam.f_p_fw_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_blkt_coolant_pump_mw",
        powerflowcalcparam.p_blkt_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_blkt_coolant_pump_total_heat",
        powerflowcalcparam.f_p_blkt_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shld_coolant_pump_mw",
        powerflowcalcparam.p_shld_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_shld_coolant_pump_total_heat",
        powerflowcalcparam.f_p_shld_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_div_coolant_pump_mw",
        powerflowcalcparam.p_div_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_div_coolant_pump_total_heat",
        powerflowcalcparam.f_p_div_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        divertor_variables, "n_divertors", powerflowcalcparam.n_divertors
    )

    monkeypatch.setattr(
        physics_variables, "p_plasma_rad_mw", powerflowcalcparam.p_plasma_rad_mw
    )

    monkeypatch.setattr(
        physics_variables, "p_fw_alpha_mw", powerflowcalcparam.p_fw_alpha_mw
    )

    monkeypatch.setattr(
        physics_variables,
        "p_plasma_separatrix_mw",
        powerflowcalcparam.p_plasma_separatrix_mw,
    )

    monkeypatch.setattr(primary_pumping_variables, "p_he", powerflowcalcparam.p_he)

    monkeypatch.setattr(primary_pumping_variables, "dp_he", powerflowcalcparam.dp_he)

    monkeypatch.setattr(
        primary_pumping_variables, "gamma_he", powerflowcalcparam.gamma_he
    )

    monkeypatch.setattr(primary_pumping_variables, "t_in_bb", powerflowcalcparam.t_in_bb)

    monkeypatch.setattr(
        primary_pumping_variables, "t_out_bb", powerflowcalcparam.t_out_bb
    )

    monkeypatch.setattr(
        primary_pumping_variables,
        "p_fw_blkt_coolant_pump_mw",
        powerflowcalcparam.p_fw_blkt_coolant_pump_mw,
    )

    ccfe_hcpb.powerflow_calc(False)

    assert fwbs_variables.p_fw_rad_total_mw == pytest.approx(
        powerflowcalcparam.expected_p_fw_rad_total_mw
    )

    assert fwbs_variables.psurffwi == pytest.approx(powerflowcalcparam.expected_psurffwi)

    assert fwbs_variables.psurffwo == pytest.approx(powerflowcalcparam.expected_psurffwo)

    assert heat_transport_variables.p_shld_coolant_pump_mw == pytest.approx(
        powerflowcalcparam.expected_p_shld_coolant_pump_mw
    )

    assert heat_transport_variables.p_div_coolant_pump_mw == pytest.approx(
        powerflowcalcparam.expected_p_div_coolant_pump_mw
    )

    assert primary_pumping_variables.p_fw_blkt_coolant_pump_mw == pytest.approx(
        powerflowcalcparam.expected_p_fw_blkt_coolant_pump_mw
    )


class StCpAngleFractionParam(NamedTuple):
    z_cp_top: Any = None

    r_cp_top: Any = None

    r_cp_mid: Any = None

    rmajor: Any = None

    expected_f_geom_cp: Any = None


@pytest.mark.parametrize(
    "stcpanglefractionparam",
    (
        StCpAngleFractionParam(
            z_cp_top=2.6714285714285717,
            r_cp_top=0.92643571428571436,
            r_cp_mid=0.20483000000000001,
            rmajor=1.7000000000000002,
            expected_f_geom_cp=0.08375588625302606,
        ),
        StCpAngleFractionParam(
            z_cp_top=2.6714285714285717,
            r_cp_top=0.92643571428571436,
            r_cp_mid=0.20483000000000001,
            rmajor=1.7000000000000002,
            expected_f_geom_cp=0.08375588625302606,
        ),
    ),
)
def test_st_cp_angle_fraction(stcpanglefractionparam, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for st_cp_angle_fraction.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param stcpanglefractionparam: the data used to mock and assert in this test.
    :type stcpanglefractionparam: stcpanglefractionparam
    """

    f_geom_cp = ccfe_hcpb.st_cp_angle_fraction(
        z_cp_top=stcpanglefractionparam.z_cp_top,
        r_cp_top=stcpanglefractionparam.r_cp_top,
        r_cp_mid=stcpanglefractionparam.r_cp_mid,
        rmajor=stcpanglefractionparam.rmajor,
    )

    assert f_geom_cp == pytest.approx(stcpanglefractionparam.expected_f_geom_cp)


class StTfCentrepostFastNeutFluxParam(NamedTuple):
    i_tf_sup: Any = None

    p_neutron_total_mw: Any = None

    sh_width: Any = None

    rmajor: Any = None

    expected_neut_flux_cp: Any = None


@pytest.mark.parametrize(
    "sttfcentrepostfastneutfluxparam",
    (
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            p_neutron_total_mw=400.65875490746737,
            sh_width=0.60000000000000009,
            rmajor=3,
            expected_neut_flux_cp=144701998710998.5,
        ),
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            p_neutron_total_mw=409.82485143909827,
            sh_width=0.60000000000000009,
            rmajor=3,
            expected_neut_flux_cp=148012428028364.28,
        ),
    ),
)
def test_st_tf_centrepost_fast_neut_flux(
    sttfcentrepostfastneutfluxparam, monkeypatch, ccfe_hcpb
):
    """
    Automatically generated Regression Unit Test for st_tf_centrepost_fast_neut_flux.

    This test was generated using data from tests/regression/scenarios/Menard_HTS-PP/IN.DAT.

    :param sttfcentrepostfastneutfluxparam: the data used to mock and assert in this test.
    :type sttfcentrepostfastneutfluxparam: sttfcentrepostfastneutfluxparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_sup", sttfcentrepostfastneutfluxparam.i_tf_sup
    )

    neut_flux_cp = ccfe_hcpb.st_tf_centrepost_fast_neut_flux(
        p_neutron_total_mw=sttfcentrepostfastneutfluxparam.p_neutron_total_mw,
        sh_width=sttfcentrepostfastneutfluxparam.sh_width,
        rmajor=sttfcentrepostfastneutfluxparam.rmajor,
    )

    assert neut_flux_cp == pytest.approx(
        sttfcentrepostfastneutfluxparam.expected_neut_flux_cp
    )


class StCentrepostNuclearHeatingParam(NamedTuple):
    rmajor: Any = None

    i_tf_sup: Any = None

    pneut: Any = None

    sh_width: Any = None

    expected_pnuc_cp_tf: Any = None

    expected_p_cp_shield_nuclear_heat_mw: Any = None

    expected_pnuc_cp: Any = None


@pytest.mark.parametrize(
    "stcentrepostnuclearheatingparam",
    (
        StCentrepostNuclearHeatingParam(
            rmajor=3,
            i_tf_sup=1,
            pneut=400.65875490746737,
            sh_width=0.60000000000000009,
            expected_pnuc_cp_tf=0.0073082167825651127,
            expected_p_cp_shield_nuclear_heat_mw=111.82272602156291,
            expected_pnuc_cp=111.83003423834548,
        ),
        StCentrepostNuclearHeatingParam(
            rmajor=3,
            i_tf_sup=1,
            pneut=409.82485143909827,
            sh_width=0.60000000000000009,
            expected_pnuc_cp_tf=0.0074754109838213612,
            expected_p_cp_shield_nuclear_heat_mw=114.3809576553144,
            expected_pnuc_cp=114.38843306629822,
        ),
    ),
)
def test_st_centrepost_nuclear_heating(
    stcentrepostnuclearheatingparam, monkeypatch, ccfe_hcpb
):
    """
    Automatically generated Regression Unit Test for st_centrepost_nuclear_heating.

    This test was generated using data from tests/regression/scenarios/Menard_HTS-PP/IN.DAT.

    :param stcentrepostnuclearheatingparam: the data used to mock and assert in this test.
    :type stcentrepostnuclearheatingparam: stcentrepostnuclearheatingparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        physics_variables, "rmajor", stcentrepostnuclearheatingparam.rmajor
    )

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_sup", stcentrepostnuclearheatingparam.i_tf_sup
    )

    pnuc_cp_tf, p_cp_shield_nuclear_heat_mw, pnuc_cp = (
        ccfe_hcpb.st_centrepost_nuclear_heating(
            pneut=stcentrepostnuclearheatingparam.pneut,
            sh_width=stcentrepostnuclearheatingparam.sh_width,
        )
    )

    assert pnuc_cp_tf == pytest.approx(
        stcentrepostnuclearheatingparam.expected_pnuc_cp_tf
    )

    assert p_cp_shield_nuclear_heat_mw == pytest.approx(
        stcentrepostnuclearheatingparam.expected_p_cp_shield_nuclear_heat_mw
    )

    assert pnuc_cp == pytest.approx(stcentrepostnuclearheatingparam.expected_pnuc_cp)


class ComponentMassesParam(NamedTuple):
    a_div_surface_total: Any = None
    f_vol_div_coolant: Any = None
    dx_div_plate: Any = None
    fdiva: Any = None
    m_div_plate: Any = None
    den_div_structure: Any = None
    rminor: Any = None
    rmajor: Any = None
    n_divertors: Any = None
    a_plasma_surface: Any = None
    dr_blkt_inboard: Any = None
    blbuith: Any = None
    blbmith: Any = None
    blbpith: Any = None
    dr_blkt_outboard: Any = None
    blbuoth: Any = None
    blbmoth: Any = None
    blbpoth: Any = None
    a_fw_inboard: Any = None
    dr_fw_inboard: Any = None
    a_fw_outboard: Any = None
    dr_fw_outboard: Any = None
    a_fw_total: Any = None
    vol_blkt_total: Any = None
    f_a_blkt_cooling_channels: Any = None
    m_blkt_beryllium: Any = None
    m_blkt_steel_total: Any = None
    den_steel: Any = None
    m_blkt_total: Any = None
    vol_shld_total: Any = None
    vfshld: Any = None
    m_fw_blkt_div_coolant_total: Any = None
    fwclfr: Any = None
    breeder_f: Any = None
    breeder_multiplier: Any = None
    m_blkt_tibe12: Any = None
    m_blkt_li4sio4: Any = None
    m_blkt_li2o: Any = None
    vfcblkt: Any = None
    vfpblkt: Any = None
    whtshld: Any = None
    wpenshld: Any = None
    m_fw_total: Any = None
    fw_armour_vol: Any = None
    fw_armour_thickness: Any = None
    fw_armour_mass: Any = None
    armour_fw_bl_mass: Any = None
    vol_blkt_inboard: Any = None
    vol_blkt_outboard: Any = None
    i_blkt_inboard: Any = None
    fblhebmi: Any = None
    fblhebpi: Any = None
    fblhebmo: Any = None
    fblhebpo: Any = None
    fblss: Any = None
    fblbe: Any = None
    whtblbreed: Any = None
    densbreed: Any = None
    fblbreed: Any = None
    i_blanket_type: Any = None
    f_a_fw_coolant_inboard: Any = None
    f_a_fw_coolant_outboard: Any = None
    vol_fw_total: Any = None
    f_vol_blkt_steel: Any = None
    f_vol_blkt_li4sio4: Any = None
    f_vol_blkt_tibe12: Any = None
    expected_a_div_surface_total: Any = None
    expected_m_div_plate: Any = None
    expected_m_blkt_beryllium: Any = None
    expected_m_blkt_steel_total: Any = None
    expected_m_blkt_total: Any = None
    expected_m_fw_blkt_div_coolant_total: Any = None
    expected_fwclfr: Any = None
    expected_m_blkt_tibe12: Any = None
    expected_whtblli4sio4: Any = None
    expected_m_blkt_li2o: Any = None
    expected_whtshld: Any = None
    expected_wpenshld: Any = None
    expected_m_fw_total: Any = None
    expected_fw_armour_vol: Any = None
    expected_fw_armour_mass: Any = None
    expected_armour_fw_bl_mass: Any = None
    expected_fblss_ccfe: Any = None
    expected_f_vol_blkt_li4sio4: Any = None
    expected_f_vol_blkt_tibe12: Any = None


@pytest.mark.parametrize(
    "componentmassesparam",
    (
        ComponentMassesParam(
            a_div_surface_total=0,
            f_vol_div_coolant=0.29999999999999999,
            dx_div_plate=0.035000000000000003,
            fdiva=1.1100000000000001,
            m_div_plate=0,
            den_div_structure=10000,
            rminor=2.6666666666666665,
            rmajor=8,
            n_divertors=1,
            a_plasma_surface=1173.8427771245592,
            dr_blkt_inboard=0.70000000000000007,
            blbuith=0.36499999999999999,
            blbmith=0.17000000000000001,
            blbpith=0.29999999999999999,
            dr_blkt_outboard=1,
            blbuoth=0.46500000000000002,
            blbmoth=0.27000000000000002,
            blbpoth=0.34999999999999998,
            a_fw_inboard=505.96109565204046,
            dr_fw_inboard=0.018000000000000002,
            a_fw_outboard=838.00728058362097,
            dr_fw_outboard=0.018000000000000002,
            a_fw_total=1343.9683762356615,
            vol_blkt_total=1182.5433772195902,
            f_a_blkt_cooling_channels=0.25,
            m_blkt_beryllium=0,
            m_blkt_steel_total=0,
            den_steel=7800,
            m_blkt_total=0,
            vol_shld_total=783.69914576548854,
            vfshld=0.60000000000000009,
            m_fw_blkt_div_coolant_total=0,
            fwclfr=0.14999999999999999,
            breeder_f=0.5,
            breeder_multiplier=0.75,
            m_blkt_tibe12=0,
            m_blkt_li4sio4=0,
            m_blkt_li2o=0,
            vfcblkt=0.052949999999999997,
            vfpblkt=0.10000000000000001,
            whtshld=0,
            wpenshld=0,
            m_fw_total=0,
            fw_armour_vol=0,
            fw_armour_thickness=0.0050000000000000001,
            fw_armour_mass=0,
            armour_fw_bl_mass=0,
            vol_blkt_inboard=315.83946385183026,
            vol_blkt_outboard=866.70391336775992,
            i_blkt_inboard=1,
            fblhebmi=0.40000000000000002,
            fblhebpi=0.65949999999999998,
            fblhebmo=0.40000000000000002,
            fblhebpo=0.67130000000000001,
            fblss=0.097049999999999997,
            fblbe=0.59999999999999998,
            whtblbreed=0,
            densbreed=0,
            fblbreed=0.154,
            i_blanket_type=1,
            f_a_fw_coolant_inboard=0,
            f_a_fw_coolant_outboard=0,
            vol_fw_total=0,
            f_vol_blkt_steel=0,
            f_vol_blkt_li4sio4=0,
            f_vol_blkt_tibe12=0,
            expected_a_div_surface_total=148.78582807401261,
            expected_m_div_plate=36452.527878133093,
            expected_m_blkt_beryllium=1002205.5121936026,
            expected_m_blkt_steel_total=895173.51112145756,
            expected_m_blkt_total=2961668.0628126911,
            expected_m_fw_blkt_div_coolant_total=1161.8025382862772,
            expected_fwclfr=0,
            expected_m_blkt_tibe12=1002205.5121936026,
            expected_whtblli4sio4=1064289.0394976311,
            expected_m_blkt_li2o=1064289.0394976311,
            expected_whtshld=2445141.3347883238,
            expected_wpenshld=2445141.3347883238,
            expected_m_fw_total=188693.16002348688,
            expected_fw_armour_vol=5.8692138856227967,
            expected_fw_armour_mass=112982.36729823884,
            expected_armour_fw_bl_mass=3263343.5901344167,
            expected_fblss_ccfe=0.097049999999999997,
            expected_f_vol_blkt_li4sio4=0.375,
            expected_f_vol_blkt_tibe12=0.375,
        ),
    ),
)
def test_component_masses(componentmassesparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for component_masses.

    This test was generated using data from tests/regression/input_files/large_tokamak_eval.IN.DAT.

    :param componentmassesparam: the data used to mock and assert in this test.
    :type componentmassesparam: componentmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        divertor_variables,
        "a_div_surface_total",
        componentmassesparam.a_div_surface_total,
    )
    monkeypatch.setattr(
        divertor_variables, "f_vol_div_coolant", componentmassesparam.f_vol_div_coolant
    )
    monkeypatch.setattr(
        divertor_variables, "dx_div_plate", componentmassesparam.dx_div_plate
    )
    monkeypatch.setattr(divertor_variables, "fdiva", componentmassesparam.fdiva)
    monkeypatch.setattr(
        divertor_variables, "m_div_plate", componentmassesparam.m_div_plate
    )
    monkeypatch.setattr(
        divertor_variables, "den_div_structure", componentmassesparam.den_div_structure
    )
    monkeypatch.setattr(physics_variables, "rminor", componentmassesparam.rminor)
    monkeypatch.setattr(physics_variables, "rmajor", componentmassesparam.rmajor)
    monkeypatch.setattr(
        divertor_variables, "n_divertors", componentmassesparam.n_divertors
    )
    monkeypatch.setattr(
        physics_variables, "a_plasma_surface", componentmassesparam.a_plasma_surface
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", componentmassesparam.dr_blkt_inboard
    )
    monkeypatch.setattr(build_variables, "blbuith", componentmassesparam.blbuith)
    monkeypatch.setattr(build_variables, "blbmith", componentmassesparam.blbmith)
    monkeypatch.setattr(build_variables, "blbpith", componentmassesparam.blbpith)
    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", componentmassesparam.dr_blkt_outboard
    )
    monkeypatch.setattr(build_variables, "blbuoth", componentmassesparam.blbuoth)
    monkeypatch.setattr(build_variables, "blbmoth", componentmassesparam.blbmoth)
    monkeypatch.setattr(build_variables, "blbpoth", componentmassesparam.blbpoth)
    monkeypatch.setattr(
        build_variables, "a_fw_inboard", componentmassesparam.a_fw_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_inboard", componentmassesparam.dr_fw_inboard
    )
    monkeypatch.setattr(
        build_variables, "a_fw_outboard", componentmassesparam.a_fw_outboard
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", componentmassesparam.dr_fw_outboard
    )
    monkeypatch.setattr(build_variables, "a_fw_total", componentmassesparam.a_fw_total)
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_total", componentmassesparam.vol_blkt_total
    )
    monkeypatch.setattr(
        fwbs_variables,
        "f_a_blkt_cooling_channels",
        componentmassesparam.f_a_blkt_cooling_channels,
    )
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_beryllium", componentmassesparam.m_blkt_beryllium
    )
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_steel_total", componentmassesparam.m_blkt_steel_total
    )
    monkeypatch.setattr(fwbs_variables, "den_steel", componentmassesparam.den_steel)
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_total", componentmassesparam.m_blkt_total
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_shld_total", componentmassesparam.vol_shld_total
    )
    monkeypatch.setattr(fwbs_variables, "vfshld", componentmassesparam.vfshld)
    monkeypatch.setattr(
        fwbs_variables,
        "m_fw_blkt_div_coolant_total",
        componentmassesparam.m_fw_blkt_div_coolant_total,
    )
    monkeypatch.setattr(fwbs_variables, "fwclfr", componentmassesparam.fwclfr)
    monkeypatch.setattr(fwbs_variables, "breeder_f", componentmassesparam.breeder_f)
    monkeypatch.setattr(
        fwbs_variables, "breeder_multiplier", componentmassesparam.breeder_multiplier
    )
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_tibe12", componentmassesparam.m_blkt_tibe12
    )
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_li4sio4", componentmassesparam.m_blkt_li4sio4
    )
    monkeypatch.setattr(fwbs_variables, "m_blkt_li2o", componentmassesparam.m_blkt_li2o)
    monkeypatch.setattr(fwbs_variables, "vfcblkt", componentmassesparam.vfcblkt)
    monkeypatch.setattr(fwbs_variables, "vfpblkt", componentmassesparam.vfpblkt)
    monkeypatch.setattr(fwbs_variables, "whtshld", componentmassesparam.whtshld)
    monkeypatch.setattr(fwbs_variables, "wpenshld", componentmassesparam.wpenshld)
    monkeypatch.setattr(fwbs_variables, "m_fw_total", componentmassesparam.m_fw_total)
    monkeypatch.setattr(
        fwbs_variables, "fw_armour_vol", componentmassesparam.fw_armour_vol
    )
    monkeypatch.setattr(
        fwbs_variables, "fw_armour_thickness", componentmassesparam.fw_armour_thickness
    )
    monkeypatch.setattr(
        fwbs_variables, "fw_armour_mass", componentmassesparam.fw_armour_mass
    )
    monkeypatch.setattr(
        fwbs_variables, "armour_fw_bl_mass", componentmassesparam.armour_fw_bl_mass
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_inboard", componentmassesparam.vol_blkt_inboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_outboard", componentmassesparam.vol_blkt_outboard
    )
    monkeypatch.setattr(
        fwbs_variables, "i_blkt_inboard", componentmassesparam.i_blkt_inboard
    )
    monkeypatch.setattr(fwbs_variables, "fblhebmi", componentmassesparam.fblhebmi)
    monkeypatch.setattr(fwbs_variables, "fblhebpi", componentmassesparam.fblhebpi)
    monkeypatch.setattr(fwbs_variables, "fblhebmo", componentmassesparam.fblhebmo)
    monkeypatch.setattr(fwbs_variables, "fblhebpo", componentmassesparam.fblhebpo)
    monkeypatch.setattr(fwbs_variables, "fblss", componentmassesparam.fblss)
    monkeypatch.setattr(fwbs_variables, "fblbe", componentmassesparam.fblbe)
    monkeypatch.setattr(fwbs_variables, "whtblbreed", componentmassesparam.whtblbreed)
    monkeypatch.setattr(fwbs_variables, "densbreed", componentmassesparam.densbreed)
    monkeypatch.setattr(fwbs_variables, "fblbreed", componentmassesparam.fblbreed)
    monkeypatch.setattr(
        fwbs_variables, "i_blanket_type", componentmassesparam.i_blanket_type
    )
    monkeypatch.setattr(
        fwbs_variables,
        "f_a_fw_coolant_inboard",
        componentmassesparam.f_a_fw_coolant_inboard,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "f_a_fw_coolant_outboard",
        componentmassesparam.f_a_fw_coolant_outboard,
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_fw_total", componentmassesparam.vol_fw_total
    )
    monkeypatch.setattr(
        fwbs_variables, "f_vol_blkt_steel", componentmassesparam.f_vol_blkt_steel
    )
    monkeypatch.setattr(
        fwbs_variables, "f_vol_blkt_li4sio4", componentmassesparam.f_vol_blkt_li4sio4
    )
    monkeypatch.setattr(
        fwbs_variables, "f_vol_blkt_tibe12", componentmassesparam.f_vol_blkt_tibe12
    )

    ccfe_hcpb.component_masses()

    assert divertor_variables.a_div_surface_total == pytest.approx(
        componentmassesparam.expected_a_div_surface_total
    )
    assert divertor_variables.m_div_plate == pytest.approx(
        componentmassesparam.expected_m_div_plate
    )
    assert fwbs_variables.m_blkt_beryllium == pytest.approx(
        componentmassesparam.expected_m_blkt_beryllium
    )
    assert fwbs_variables.m_blkt_steel_total == pytest.approx(
        componentmassesparam.expected_m_blkt_steel_total
    )
    assert fwbs_variables.m_blkt_total == pytest.approx(
        componentmassesparam.expected_m_blkt_total
    )
    assert fwbs_variables.m_fw_blkt_div_coolant_total == pytest.approx(
        componentmassesparam.expected_m_fw_blkt_div_coolant_total
    )
    assert fwbs_variables.fwclfr == pytest.approx(componentmassesparam.expected_fwclfr)
    assert fwbs_variables.m_blkt_tibe12 == pytest.approx(
        componentmassesparam.expected_m_blkt_tibe12
    )
    assert fwbs_variables.m_blkt_li4sio4 == pytest.approx(
        componentmassesparam.expected_whtblli4sio4
    )
    assert fwbs_variables.m_blkt_li2o == pytest.approx(
        componentmassesparam.expected_m_blkt_li2o
    )
    assert fwbs_variables.whtshld == pytest.approx(componentmassesparam.expected_whtshld)
    assert fwbs_variables.wpenshld == pytest.approx(
        componentmassesparam.expected_wpenshld
    )
    assert fwbs_variables.m_fw_total == pytest.approx(
        componentmassesparam.expected_m_fw_total
    )
    assert fwbs_variables.fw_armour_vol == pytest.approx(
        componentmassesparam.expected_fw_armour_vol
    )
    assert fwbs_variables.fw_armour_mass == pytest.approx(
        componentmassesparam.expected_fw_armour_mass
    )
    assert fwbs_variables.armour_fw_bl_mass == pytest.approx(
        componentmassesparam.expected_armour_fw_bl_mass
    )
    assert fwbs_variables.f_vol_blkt_steel == pytest.approx(
        componentmassesparam.expected_fblss_ccfe
    )
    assert fwbs_variables.f_vol_blkt_li4sio4 == pytest.approx(
        componentmassesparam.expected_f_vol_blkt_li4sio4
    )
    assert fwbs_variables.f_vol_blkt_tibe12 == pytest.approx(
        componentmassesparam.expected_f_vol_blkt_tibe12
    )
