import pytest
from typing import NamedTuple, Any

from process.hcpb import CCFE_HCPB
from process.blanket_library import BlanketLibrary
from process.fw import Fw
from process.fortran import (
    fwbs_variables,
    build_variables,
    global_variables,
    tfcoil_variables,
    physics_variables,
    ccfe_hcpb_module,
    primary_pumping_variables,
    current_drive_variables,
    heat_transport_variables,
    constraint_variables,
    divertor_variables,
)


@pytest.fixture
def ccfe_hcpb():
    """Provides CCFE_HCPB object for testing.

    :returns: initialised CCFE_HCPB object
    :rtype: process.hcpb.CCFE_HCPB
    """
    return CCFE_HCPB(BlanketLibrary(Fw()))


class NuclearHeatingMagnetsParam(NamedTuple):
    fwith: Any = None

    d_vv_in: Any = None

    d_vv_out: Any = None

    fwoth: Any = None

    blnkith: Any = None

    blnkoth: Any = None

    shldith: Any = None

    shldoth: Any = None

    afw: Any = None

    pitch: Any = None

    denstl: Any = None

    whtblkt: Any = None

    volblkt: Any = None

    whtshld: Any = None

    volshld: Any = None

    vvmass: Any = None

    vdewin: Any = None

    fw_armour_thickness: Any = None

    ptfnuc: Any = None

    denw: Any = None

    vffwi: Any = None

    vffwo: Any = None

    fusion_power: Any = None

    itart: Any = None

    whttf: Any = None

    whttflgs: Any = None

    verbose: Any = None

    ip: Any = None

    ofile: Any = None

    armour_density: Any = None

    fw_density: Any = None

    blanket_density: Any = None

    shield_density: Any = None

    vv_density: Any = None

    x_blanket: Any = None

    x_shield: Any = None

    tfc_nuc_heating: Any = None

    expected_ptfnuc: Any = None

    expected_vffwi: Any = None

    expected_vffwo: Any = None

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
            fwith=0.018000000000000002,
            d_vv_in=0.30000000000000004,
            d_vv_out=0.30000000000000004,
            fwoth=0.018000000000000002,
            blnkith=0.75500000000000012,
            blnkoth=0.98199999999999998,
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            afw=0.0060000000000000001,
            pitch=0.02,
            denstl=7800,
            whtblkt=3501027.3252278985,
            volblkt=1397.9003011502937,
            whtshld=2294873.8131476045,
            volshld=735.53647857295027,
            vvmass=9043937.8018644415,
            vdewin=1159.4792053672361,
            fw_armour_thickness=0.0050000000000000001,
            ptfnuc=0,
            denw=19250,
            vffwi=0,
            vffwo=0,
            fusion_power=1986.0623241661431,
            itart=0,
            whttf=19649856.627845347,
            whttflgs=0,
            verbose=0,
            ip=0,
            ofile=11,
            armour_density=0,
            fw_density=0,
            blanket_density=0,
            shield_density=0,
            vv_density=0,
            x_blanket=0,
            x_shield=0,
            tfc_nuc_heating=0,
            expected_ptfnuc=0.044541749095475737,
            expected_vffwi=0.31415926535897931,
            expected_vffwo=0.31415926535897931,
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
            fwith=0.018000000000000002,
            d_vv_in=0.30000000000000004,
            d_vv_out=0.30000000000000004,
            fwoth=0.018000000000000002,
            blnkith=0.75500000000000012,
            blnkoth=0.98199999999999998,
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            afw=0.0060000000000000001,
            pitch=0.02,
            denstl=7800,
            whtblkt=3507503.3737008357,
            volblkt=1400.4860764869636,
            whtshld=2297808.3935174854,
            volshld=736.47704920432227,
            vvmass=9056931.558219457,
            vdewin=1161.1450715665972,
            fw_armour_thickness=0.0050000000000000001,
            ptfnuc=0.044184461825198453,
            denw=19250,
            vffwi=0.31415926535897931,
            vffwo=0.31415926535897931,
            fusion_power=1985.4423932312809,
            itart=0,
            whttf=19662548.210142396,
            whttflgs=0,
            verbose=0,
            ip=0,
            ofile=11,
            armour_density=13202.434141839649,
            fw_density=5349.557730199961,
            blanket_density=2504.4899999999998,
            shield_density=3119.9999999999995,
            vv_density=7800,
            x_blanket=2.3374537748527975,
            x_shield=4.056,
            tfc_nuc_heating=22427.165831352642,
            expected_ptfnuc=0.044556605747797934,
            expected_vffwi=0.31415926535897931,
            expected_vffwo=0.31415926535897931,
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

    monkeypatch.setattr(build_variables, "fwith", nuclearheatingmagnetsparam.fwith)

    monkeypatch.setattr(build_variables, "d_vv_in", nuclearheatingmagnetsparam.d_vv_in)

    monkeypatch.setattr(
        build_variables, "d_vv_out", nuclearheatingmagnetsparam.d_vv_out
    )

    monkeypatch.setattr(build_variables, "fwoth", nuclearheatingmagnetsparam.fwoth)

    monkeypatch.setattr(build_variables, "blnkith", nuclearheatingmagnetsparam.blnkith)

    monkeypatch.setattr(build_variables, "blnkoth", nuclearheatingmagnetsparam.blnkoth)

    monkeypatch.setattr(build_variables, "shldith", nuclearheatingmagnetsparam.shldith)

    monkeypatch.setattr(build_variables, "shldoth", nuclearheatingmagnetsparam.shldoth)

    monkeypatch.setattr(fwbs_variables, "afw", nuclearheatingmagnetsparam.afw)

    monkeypatch.setattr(fwbs_variables, "pitch", nuclearheatingmagnetsparam.pitch)

    monkeypatch.setattr(fwbs_variables, "denstl", nuclearheatingmagnetsparam.denstl)

    monkeypatch.setattr(fwbs_variables, "whtblkt", nuclearheatingmagnetsparam.whtblkt)

    monkeypatch.setattr(fwbs_variables, "volblkt", nuclearheatingmagnetsparam.volblkt)

    monkeypatch.setattr(fwbs_variables, "whtshld", nuclearheatingmagnetsparam.whtshld)

    monkeypatch.setattr(fwbs_variables, "volshld", nuclearheatingmagnetsparam.volshld)

    monkeypatch.setattr(fwbs_variables, "vvmass", nuclearheatingmagnetsparam.vvmass)

    monkeypatch.setattr(fwbs_variables, "vdewin", nuclearheatingmagnetsparam.vdewin)

    monkeypatch.setattr(
        fwbs_variables,
        "fw_armour_thickness",
        nuclearheatingmagnetsparam.fw_armour_thickness,
    )

    monkeypatch.setattr(fwbs_variables, "ptfnuc", nuclearheatingmagnetsparam.ptfnuc)

    monkeypatch.setattr(fwbs_variables, "denw", nuclearheatingmagnetsparam.denw)

    monkeypatch.setattr(fwbs_variables, "vffwi", nuclearheatingmagnetsparam.vffwi)

    monkeypatch.setattr(fwbs_variables, "vffwo", nuclearheatingmagnetsparam.vffwo)

    monkeypatch.setattr(
        physics_variables, "fusion_power", nuclearheatingmagnetsparam.fusion_power
    )

    monkeypatch.setattr(physics_variables, "itart", nuclearheatingmagnetsparam.itart)

    monkeypatch.setattr(tfcoil_variables, "whttf", nuclearheatingmagnetsparam.whttf)

    monkeypatch.setattr(
        tfcoil_variables, "whttflgs", nuclearheatingmagnetsparam.whttflgs
    )

    monkeypatch.setattr(global_variables, "verbose", nuclearheatingmagnetsparam.verbose)

    monkeypatch.setattr(ccfe_hcpb_module, "ip", nuclearheatingmagnetsparam.ip)

    monkeypatch.setattr(ccfe_hcpb_module, "ofile", nuclearheatingmagnetsparam.ofile)

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

    assert fwbs_variables.ptfnuc == pytest.approx(
        nuclearheatingmagnetsparam.expected_ptfnuc
    )

    assert fwbs_variables.vffwi == pytest.approx(
        nuclearheatingmagnetsparam.expected_vffwi
    )

    assert fwbs_variables.vffwo == pytest.approx(
        nuclearheatingmagnetsparam.expected_vffwo
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
    p_fw_nuclear_heat_mw: Any = None

    fwmass: Any = None

    fusion_power: Any = None

    fw_armour_u_nuc_heating: Any = None

    expected_p_fw_nuclear_heat_mw: Any = None

    expected_fw_armour_u_nuc_heating: Any = None


@pytest.mark.parametrize(
    "nuclearheatingfwparam",
    (
        NuclearHeatingFwParam(
            p_fw_nuclear_heat_mw=0,
            fwmass=224802.80270851994,
            fusion_power=1986.0623241661431,
            fw_armour_u_nuc_heating=0,
            expected_p_fw_nuclear_heat_mw=279.04523551646628,
            expected_fw_armour_u_nuc_heating=6.2500000000000005e-07,
        ),
        NuclearHeatingFwParam(
            p_fw_nuclear_heat_mw=276.80690153753221,
            fwmass=182115.83467868491,
            fusion_power=1985.4423932312809,
            fw_armour_u_nuc_heating=6.2500000000000005e-07,
            expected_p_fw_nuclear_heat_mw=225.98781165610032,
            expected_fw_armour_u_nuc_heating=6.2500000000000005e-07,
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
        "p_fw_nuclear_heat_mw",
        nuclearheatingfwparam.p_fw_nuclear_heat_mw,
    )

    monkeypatch.setattr(fwbs_variables, "fwmass", nuclearheatingfwparam.fwmass)

    monkeypatch.setattr(
        physics_variables, "fusion_power", nuclearheatingfwparam.fusion_power
    )

    monkeypatch.setattr(
        ccfe_hcpb_module,
        "fw_armour_u_nuc_heating",
        nuclearheatingfwparam.fw_armour_u_nuc_heating,
    )

    ccfe_hcpb.nuclear_heating_fw()

    assert fwbs_variables.p_fw_nuclear_heat_mw == pytest.approx(
        nuclearheatingfwparam.expected_p_fw_nuclear_heat_mw
    )

    assert ccfe_hcpb_module.fw_armour_u_nuc_heating == pytest.approx(
        nuclearheatingfwparam.expected_fw_armour_u_nuc_heating
    )


class NuclearHeatingBlanketParam(NamedTuple):
    whtblkt: Any = None

    p_blanket_nuclear_heat_mw: Any = None

    fusion_power: Any = None

    exp_blanket: Any = None

    expected_p_blanket_nuclear_heat_mw: Any = None

    expected_exp_blanket: Any = None


@pytest.mark.parametrize(
    "nuclearheatingblanketparam",
    (
        NuclearHeatingBlanketParam(
            whtblkt=3501027.3252278985,
            p_blanket_nuclear_heat_mw=0,
            fusion_power=1986.0623241661431,
            exp_blanket=0,
            expected_p_blanket_nuclear_heat_mw=1517.0907688379014,
            expected_exp_blanket=0.99982809071915879,
        ),
        NuclearHeatingBlanketParam(
            whtblkt=3507503.3737008357,
            p_blanket_nuclear_heat_mw=1504.9215740808861,
            fusion_power=1985.4423932312809,
            exp_blanket=0.99982809071915879,
            expected_p_blanket_nuclear_heat_mw=1516.6213709741428,
            expected_exp_blanket=0.99983082524994527,
        ),
    ),
)
def test_nuclear_heating_blanket(nuclearheatingblanketparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_blanket.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingblanketparam: the data used to mock and assert in this test.
    :type nuclearheatingblanketparam: nuclearheatingblanketparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(fwbs_variables, "whtblkt", nuclearheatingblanketparam.whtblkt)

    monkeypatch.setattr(
        fwbs_variables,
        "p_blanket_nuclear_heat_mw",
        nuclearheatingblanketparam.p_blanket_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        physics_variables, "fusion_power", nuclearheatingblanketparam.fusion_power
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "exp_blanket", nuclearheatingblanketparam.exp_blanket
    )

    ccfe_hcpb.nuclear_heating_blanket()

    assert fwbs_variables.p_blanket_nuclear_heat_mw == pytest.approx(
        nuclearheatingblanketparam.expected_p_blanket_nuclear_heat_mw
    )

    assert ccfe_hcpb_module.exp_blanket == pytest.approx(
        nuclearheatingblanketparam.expected_exp_blanket
    )


class NuclearHeatingShieldParam(NamedTuple):
    shldith: Any = None

    shldoth: Any = None

    whtshld: Any = None

    p_shield_nuclear_heat_mw: Any = None

    fusion_power: Any = None

    itart: Any = None

    shield_density: Any = None

    x_blanket: Any = None

    shld_u_nuc_heating: Any = None

    exp_shield1: Any = None

    exp_shield2: Any = None

    expected_p_shield_nuclear_heat_mw: Any = None

    expected_shld_u_nuc_heating: Any = None

    expected_exp_shield1: Any = None

    expected_exp_shield2: Any = None


@pytest.mark.parametrize(
    "nuclearheatingshieldparam",
    (
        NuclearHeatingShieldParam(
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            whtshld=2294873.8131476045,
            p_shield_nuclear_heat_mw=0,
            fusion_power=1986.0623241661431,
            itart=0,
            shield_density=3119.9999999999995,
            x_blanket=2.3374537748527975,
            shld_u_nuc_heating=0,
            exp_shield1=0,
            exp_shield2=0,
            expected_p_shield_nuclear_heat_mw=1.3721323841005297,
            expected_shld_u_nuc_heating=690880.82856444374,
            expected_exp_shield1=0.0017209365527675318,
            expected_exp_shield2=0.25426760591013942,
        ),
        NuclearHeatingShieldParam(
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            whtshld=2297808.3935174854,
            p_shield_nuclear_heat_mw=1.3611259588044891,
            fusion_power=1985.4423932312809,
            itart=0,
            shield_density=3120,
            x_blanket=2.3374537748527979,
            shld_u_nuc_heating=690880.82856444374,
            exp_shield1=0.0017209365527675318,
            exp_shield2=0.25426760591013942,
            expected_p_shield_nuclear_heat_mw=1.3734581585671393,
            expected_shld_u_nuc_heating=691764.29557941214,
            expected_exp_shield1=0.0017209365527675303,
            expected_exp_shield2=0.25426760591013942,
        ),
    ),
)
def test_nuclear_heating_shield(nuclearheatingshieldparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_shield.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingshieldparam: the data used to mock and assert in this test.
    :type nuclearheatingshieldparam: nuclearheatingshieldparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "shldith", nuclearheatingshieldparam.shldith)

    monkeypatch.setattr(build_variables, "shldoth", nuclearheatingshieldparam.shldoth)

    monkeypatch.setattr(fwbs_variables, "whtshld", nuclearheatingshieldparam.whtshld)

    monkeypatch.setattr(
        fwbs_variables,
        "p_shield_nuclear_heat_mw",
        nuclearheatingshieldparam.p_shield_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        physics_variables, "fusion_power", nuclearheatingshieldparam.fusion_power
    )

    monkeypatch.setattr(physics_variables, "itart", nuclearheatingshieldparam.itart)

    monkeypatch.setattr(
        ccfe_hcpb_module, "shield_density", nuclearheatingshieldparam.shield_density
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "x_blanket", nuclearheatingshieldparam.x_blanket
    )

    monkeypatch.setattr(
        ccfe_hcpb_module,
        "shld_u_nuc_heating",
        nuclearheatingshieldparam.shld_u_nuc_heating,
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "exp_shield1", nuclearheatingshieldparam.exp_shield1
    )

    monkeypatch.setattr(
        ccfe_hcpb_module, "exp_shield2", nuclearheatingshieldparam.exp_shield2
    )

    ccfe_hcpb.nuclear_heating_shield()

    assert fwbs_variables.p_shield_nuclear_heat_mw == pytest.approx(
        nuclearheatingshieldparam.expected_p_shield_nuclear_heat_mw
    )

    assert ccfe_hcpb_module.shld_u_nuc_heating == pytest.approx(
        nuclearheatingshieldparam.expected_shld_u_nuc_heating
    )

    assert ccfe_hcpb_module.exp_shield1 == pytest.approx(
        nuclearheatingshieldparam.expected_exp_shield1
    )

    assert ccfe_hcpb_module.exp_shield2 == pytest.approx(
        nuclearheatingshieldparam.expected_exp_shield2
    )


class NuclearHeatingDivertorParam(NamedTuple):
    fdiv: Any = None

    pnucdiv: Any = None

    pnuchcd: Any = None

    idivrt: Any = None

    fusion_power: Any = None

    ip: Any = None

    expected_pnucdiv: Any = None


@pytest.mark.parametrize(
    "nuclearheatingdivertorparam",
    (
        NuclearHeatingDivertorParam(
            fdiv=0.115,
            pnucdiv=0,
            pnuchcd=0,
            idivrt=1,
            fusion_power=1986.0623241661431,
            ip=0,
            expected_pnucdiv=182.71773382328519,
        ),
        NuclearHeatingDivertorParam(
            fdiv=0.115,
            pnucdiv=182.71773382328519,
            pnuchcd=0,
            idivrt=1,
            fusion_power=1985.4423932312809,
            ip=0,
            expected_pnucdiv=182.66070017727785,
        ),
    ),
)
def test_nuclear_heating_divertor(nuclearheatingdivertorparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for nuclear_heating_divertor.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param nuclearheatingdivertorparam: the data used to mock and assert in this test.
    :type nuclearheatingdivertorparam: nuclearheatingdivertorparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(fwbs_variables, "fdiv", nuclearheatingdivertorparam.fdiv)

    monkeypatch.setattr(fwbs_variables, "pnucdiv", nuclearheatingdivertorparam.pnucdiv)

    monkeypatch.setattr(fwbs_variables, "pnuchcd", nuclearheatingdivertorparam.pnuchcd)

    monkeypatch.setattr(physics_variables, "idivrt", nuclearheatingdivertorparam.idivrt)

    monkeypatch.setattr(
        physics_variables, "fusion_power", nuclearheatingdivertorparam.fusion_power
    )

    monkeypatch.setattr(ccfe_hcpb_module, "ip", nuclearheatingdivertorparam.ip)

    ccfe_hcpb.nuclear_heating_divertor()

    assert fwbs_variables.pnucdiv == pytest.approx(
        nuclearheatingdivertorparam.expected_pnucdiv
    )


class PowerflowCalcParam(NamedTuple):
    fwareaob: Any = None

    fwarea: Any = None

    p_nb_orbit_loss_mw: Any = None

    fdiv: Any = None

    praddiv: Any = None

    pradhcd: Any = None

    fhcd: Any = None

    p_fw_radiation_mw: Any = None

    coolwh: Any = None

    outlet_temp: Any = None

    blpressure: Any = None

    primary_pumping: Any = None

    p_fw_nuclear_heat_mw: Any = None

    p_blanket_nuclear_heat_mw: Any = None

    pnucdiv: Any = None

    p_shield_nuclear_heat_mw: Any = None

    etaiso: Any = None

    pnuc_cp_sh: Any = None

    psurffwi: Any = None

    psurffwo: Any = None

    p_fw_pumping_mw: Any = None

    fpumpfw: Any = None

    p_blanket_pumping_mw: Any = None

    fpumpblkt: Any = None

    p_shield_pumping_mw: Any = None

    fpumpshld: Any = None

    p_div_pump_cool_mw: Any = None

    fpumpdiv: Any = None

    idivrt: Any = None

    pradmw: Any = None

    p_fw_alpha_mw: Any = None

    pdivt: Any = None

    p_he: Any = None

    dp_he: Any = None

    gamma_he: Any = None

    t_in_bb: Any = None

    t_out_bb: Any = None

    p_fw_blanket_pumping_mw: Any = None

    ip: Any = None

    ofile: Any = None

    expected_praddiv: Any = None

    expected_p_fw_radiation_mw: Any = None

    expected_psurffwi: Any = None

    expected_psurffwo: Any = None

    expected_p_shield_pumping_mw: Any = None

    expected_p_div_pump_cool_mw: Any = None

    expected_p_fw_blanket_pumping_mw: Any = None


@pytest.mark.parametrize(
    "powerflowcalcparam",
    (
        PowerflowCalcParam(
            fwareaob=988.92586580655245,
            fwarea=1601.1595634509963,
            p_nb_orbit_loss_mw=0,
            fdiv=0.115,
            praddiv=0,
            pradhcd=0,
            fhcd=0,
            p_fw_radiation_mw=0,
            coolwh=1,
            outlet_temp=823,
            blpressure=15500000,
            primary_pumping=3,
            p_fw_nuclear_heat_mw=276.80690153753221,
            p_blanket_nuclear_heat_mw=1504.9215740808861,
            pnucdiv=182.71773382328519,
            p_shield_nuclear_heat_mw=1.3611259588044891,
            etaiso=0.90000000000000002,
            pnuc_cp_sh=0,
            psurffwi=0,
            psurffwo=0,
            p_fw_pumping_mw=0,
            fpumpfw=0.0050000000000000001,
            p_blanket_pumping_mw=0,
            fpumpblkt=0.0050000000000000001,
            p_shield_pumping_mw=0,
            fpumpshld=0.0050000000000000001,
            p_div_pump_cool_mw=0,
            fpumpdiv=0.0050000000000000001,
            idivrt=1,
            pradmw=287.44866938104849,
            p_fw_alpha_mw=19.835845058655043,
            pdivt=143.6315222649435,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            p_fw_blanket_pumping_mw=0,
            ip=0,
            ofile=11,
            expected_praddiv=33.056596978820579,
            expected_p_fw_radiation_mw=254.39207240222791,
            expected_psurffwi=97.271629070225231,
            expected_psurffwo=176.95628839065773,
            expected_p_shield_pumping_mw=0.0068056297940224456,
            expected_p_div_pump_cool_mw=1.7970292653352464,
            expected_p_fw_blanket_pumping_mw=202.00455086503842,
        ),
        PowerflowCalcParam(
            fwareaob=1168.1172772224481,
            fwarea=1891.2865102700493,
            p_nb_orbit_loss_mw=0,
            fdiv=0.115,
            praddiv=33.056596978820579,
            pradhcd=0,
            fhcd=0,
            p_fw_radiation_mw=254.39207240222791,
            coolwh=1,
            outlet_temp=823,
            blpressure=15500000,
            primary_pumping=3,
            p_fw_nuclear_heat_mw=230.98304919926957,
            p_blanket_nuclear_heat_mw=1550.1447895848396,
            pnucdiv=182.66070017727785,
            p_shield_nuclear_heat_mw=1.4038170956592293,
            etaiso=0.90000000000000002,
            pnuc_cp_sh=0,
            psurffwi=97.271629070225231,
            psurffwo=176.95628839065773,
            p_fw_pumping_mw=0,
            fpumpfw=0.0050000000000000001,
            p_blanket_pumping_mw=0,
            fpumpblkt=0.0050000000000000001,
            p_shield_pumping_mw=0.0068056297940224456,
            fpumpshld=0.0050000000000000001,
            p_div_pump_cool_mw=1.7970292653352464,
            fpumpdiv=0.0050000000000000001,
            idivrt=1,
            pradmw=287.44866938104849,
            p_fw_alpha_mw=19.829653483586444,
            pdivt=143.51338080047339,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            p_fw_blanket_pumping_mw=202.00455086503842,
            ip=0,
            ofile=11,
            expected_praddiv=33.056596978820579,
            expected_p_fw_radiation_mw=254.39207240222791,
            expected_psurffwi=97.271629070225259,
            expected_psurffwo=176.95009681558912,
            expected_p_shield_pumping_mw=0.007019085478296147,
            expected_p_div_pump_cool_mw=1.7961533897828594,
            expected_p_fw_blanket_pumping_mw=201.94492795635171,
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

    monkeypatch.setattr(build_variables, "fwareaob", powerflowcalcparam.fwareaob)

    monkeypatch.setattr(build_variables, "fwarea", powerflowcalcparam.fwarea)

    monkeypatch.setattr(
        current_drive_variables,
        "p_nb_orbit_loss_mw",
        powerflowcalcparam.p_nb_orbit_loss_mw,
    )

    monkeypatch.setattr(fwbs_variables, "fdiv", powerflowcalcparam.fdiv)

    monkeypatch.setattr(fwbs_variables, "praddiv", powerflowcalcparam.praddiv)

    monkeypatch.setattr(fwbs_variables, "pradhcd", powerflowcalcparam.pradhcd)

    monkeypatch.setattr(fwbs_variables, "fhcd", powerflowcalcparam.fhcd)

    monkeypatch.setattr(
        fwbs_variables, "p_fw_radiation_mw", powerflowcalcparam.p_fw_radiation_mw
    )

    monkeypatch.setattr(fwbs_variables, "coolwh", powerflowcalcparam.coolwh)

    monkeypatch.setattr(fwbs_variables, "outlet_temp", powerflowcalcparam.outlet_temp)

    monkeypatch.setattr(fwbs_variables, "blpressure", powerflowcalcparam.blpressure)

    monkeypatch.setattr(
        fwbs_variables, "primary_pumping", powerflowcalcparam.primary_pumping
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_nuclear_heat_mw", powerflowcalcparam.p_fw_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_blanket_nuclear_heat_mw",
        powerflowcalcparam.p_blanket_nuclear_heat_mw,
    )

    monkeypatch.setattr(fwbs_variables, "pnucdiv", powerflowcalcparam.pnucdiv)

    monkeypatch.setattr(
        fwbs_variables,
        "p_shield_nuclear_heat_mw",
        powerflowcalcparam.p_shield_nuclear_heat_mw,
    )

    monkeypatch.setattr(fwbs_variables, "etaiso", powerflowcalcparam.etaiso)

    monkeypatch.setattr(fwbs_variables, "pnuc_cp_sh", powerflowcalcparam.pnuc_cp_sh)

    monkeypatch.setattr(fwbs_variables, "psurffwi", powerflowcalcparam.psurffwi)

    monkeypatch.setattr(fwbs_variables, "psurffwo", powerflowcalcparam.psurffwo)

    monkeypatch.setattr(
        heat_transport_variables, "p_fw_pumping_mw", powerflowcalcparam.p_fw_pumping_mw
    )

    monkeypatch.setattr(heat_transport_variables, "fpumpfw", powerflowcalcparam.fpumpfw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_blanket_pumping_mw",
        powerflowcalcparam.p_blanket_pumping_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpblkt", powerflowcalcparam.fpumpblkt
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shield_pumping_mw",
        powerflowcalcparam.p_shield_pumping_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpshld", powerflowcalcparam.fpumpshld
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_div_pump_cool_mw",
        powerflowcalcparam.p_div_pump_cool_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpdiv", powerflowcalcparam.fpumpdiv
    )

    monkeypatch.setattr(physics_variables, "idivrt", powerflowcalcparam.idivrt)

    monkeypatch.setattr(physics_variables, "pradmw", powerflowcalcparam.pradmw)

    monkeypatch.setattr(
        physics_variables, "p_fw_alpha_mw", powerflowcalcparam.p_fw_alpha_mw
    )

    monkeypatch.setattr(physics_variables, "pdivt", powerflowcalcparam.pdivt)

    monkeypatch.setattr(primary_pumping_variables, "p_he", powerflowcalcparam.p_he)

    monkeypatch.setattr(primary_pumping_variables, "dp_he", powerflowcalcparam.dp_he)

    monkeypatch.setattr(
        primary_pumping_variables, "gamma_he", powerflowcalcparam.gamma_he
    )

    monkeypatch.setattr(
        primary_pumping_variables, "t_in_bb", powerflowcalcparam.t_in_bb
    )

    monkeypatch.setattr(
        primary_pumping_variables, "t_out_bb", powerflowcalcparam.t_out_bb
    )

    monkeypatch.setattr(
        primary_pumping_variables,
        "p_fw_blanket_pumping_mw",
        powerflowcalcparam.p_fw_blanket_pumping_mw,
    )

    monkeypatch.setattr(ccfe_hcpb_module, "ip", powerflowcalcparam.ip)

    monkeypatch.setattr(ccfe_hcpb_module, "ofile", powerflowcalcparam.ofile)

    ccfe_hcpb.powerflow_calc(False)

    assert fwbs_variables.praddiv == pytest.approx(powerflowcalcparam.expected_praddiv)

    assert fwbs_variables.p_fw_radiation_mw == pytest.approx(
        powerflowcalcparam.expected_p_fw_radiation_mw
    )

    assert fwbs_variables.psurffwi == pytest.approx(
        powerflowcalcparam.expected_psurffwi
    )

    assert fwbs_variables.psurffwo == pytest.approx(
        powerflowcalcparam.expected_psurffwo
    )

    assert heat_transport_variables.p_shield_pumping_mw == pytest.approx(
        powerflowcalcparam.expected_p_shield_pumping_mw
    )

    assert heat_transport_variables.p_div_pump_cool_mw == pytest.approx(
        powerflowcalcparam.expected_p_div_pump_cool_mw
    )

    assert primary_pumping_variables.p_fw_blanket_pumping_mw == pytest.approx(
        powerflowcalcparam.expected_p_fw_blanket_pumping_mw
    )


class StCpAngleFractionParam(NamedTuple):
    h_cp_top: Any = None

    r_cp_top: Any = None

    r_cp_mid: Any = None

    rmajor: Any = None

    expected_f_geom_cp: Any = None


@pytest.mark.parametrize(
    "stcpanglefractionparam",
    (
        StCpAngleFractionParam(
            h_cp_top=2.6714285714285717,
            r_cp_top=0.92643571428571436,
            r_cp_mid=0.20483000000000001,
            rmajor=1.7000000000000002,
            expected_f_geom_cp=0.08375588625302606,
        ),
        StCpAngleFractionParam(
            h_cp_top=2.6714285714285717,
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
        h_cp_top=stcpanglefractionparam.h_cp_top,
        r_cp_top=stcpanglefractionparam.r_cp_top,
        r_cp_mid=stcpanglefractionparam.r_cp_mid,
        rmajor=stcpanglefractionparam.rmajor,
    )

    assert f_geom_cp == pytest.approx(stcpanglefractionparam.expected_f_geom_cp)


class StTfCentrepostFastNeutFluxParam(NamedTuple):
    i_tf_sup: Any = None

    neutron_power_total: Any = None

    sh_width: Any = None

    rmajor: Any = None

    expected_neut_flux_cp: Any = None


@pytest.mark.parametrize(
    "sttfcentrepostfastneutfluxparam",
    (
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            neutron_power_total=400.65875490746737,
            sh_width=0.60000000000000009,
            rmajor=3,
            expected_neut_flux_cp=144701998710998.5,
        ),
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            neutron_power_total=409.82485143909827,
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
        neutron_power_total=sttfcentrepostfastneutfluxparam.neutron_power_total,
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

    expected_pnuc_cp_sh: Any = None

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
            expected_pnuc_cp_sh=111.82272602156291,
            expected_pnuc_cp=111.83003423834548,
        ),
        StCentrepostNuclearHeatingParam(
            rmajor=3,
            i_tf_sup=1,
            pneut=409.82485143909827,
            sh_width=0.60000000000000009,
            expected_pnuc_cp_tf=0.0074754109838213612,
            expected_pnuc_cp_sh=114.3809576553144,
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

    pnuc_cp_tf, pnuc_cp_sh, pnuc_cp = ccfe_hcpb.st_centrepost_nuclear_heating(
        pneut=stcentrepostnuclearheatingparam.pneut,
        sh_width=stcentrepostnuclearheatingparam.sh_width,
    )

    assert pnuc_cp_tf == pytest.approx(
        stcentrepostnuclearheatingparam.expected_pnuc_cp_tf
    )

    assert pnuc_cp_sh == pytest.approx(
        stcentrepostnuclearheatingparam.expected_pnuc_cp_sh
    )

    assert pnuc_cp == pytest.approx(stcentrepostnuclearheatingparam.expected_pnuc_cp)


class TbrShimwellParam(NamedTuple):
    fwith: Any = None

    fwoth: Any = None

    tbrmin: Any = None

    fw_armour_thickness: Any = None

    ip: Any = None

    iprint: Any = None

    outfile: Any = None

    iblanket_thickness: Any = None

    breeder_f: Any = None

    li6enrich: Any = None

    expected_tbr: Any = None


@pytest.mark.parametrize(
    "tbrshimwellparam",
    (
        TbrShimwellParam(
            fwith=0.018000000000000002,
            fwoth=0.018000000000000002,
            tbrmin=1.1499999999999999,
            fw_armour_thickness=0.0030000000000000001,
            ip=0,
            iprint=0,
            outfile=11,
            iblanket_thickness=1,
            breeder_f=0.56366688384345121,
            li6enrich=82.131743925121199,
            expected_tbr=1.1284864235692258,
        ),
    ),
)
def test_tbr_shimwell(tbrshimwellparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for tbr_shimwell.

    This test was generated using data from tests/regression/scenarios/vacuum_model/IN.DAT.

    :param tbrshimwellparam: the data used to mock and assert in this test.
    :type tbrshimwellparam: tbrshimwellparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "fwith", tbrshimwellparam.fwith)

    monkeypatch.setattr(build_variables, "fwoth", tbrshimwellparam.fwoth)

    monkeypatch.setattr(constraint_variables, "tbrmin", tbrshimwellparam.tbrmin)

    monkeypatch.setattr(
        fwbs_variables, "fw_armour_thickness", tbrshimwellparam.fw_armour_thickness
    )

    monkeypatch.setattr(ccfe_hcpb_module, "ip", tbrshimwellparam.ip)

    tbr = ccfe_hcpb.tbr_shimwell(
        iblanket_thickness=tbrshimwellparam.iblanket_thickness,
        breeder_f=tbrshimwellparam.breeder_f,
        li6enrich=tbrshimwellparam.li6enrich,
        output=False,
    )

    assert tbr == pytest.approx(tbrshimwellparam.expected_tbr)


class ComponentMassesParam(NamedTuple):
    divsur: Any = None
    divclfr: Any = None
    divplt: Any = None
    fdiva: Any = None
    divmas: Any = None
    divdens: Any = None
    rminor: Any = None
    rmajor: Any = None
    idivrt: Any = None
    sarea: Any = None
    blnkith: Any = None
    blbuith: Any = None
    blbmith: Any = None
    blbpith: Any = None
    blnkoth: Any = None
    blbuoth: Any = None
    blbmoth: Any = None
    blbpoth: Any = None
    fwareaib: Any = None
    fwith: Any = None
    fwareaob: Any = None
    fwoth: Any = None
    fwarea: Any = None
    volblkt: Any = None
    vfblkt: Any = None
    whtblbe: Any = None
    whtblss: Any = None
    denstl: Any = None
    whtblkt: Any = None
    volshld: Any = None
    vfshld: Any = None
    coolmass: Any = None
    fwclfr: Any = None
    breeder_f: Any = None
    breeder_multiplier: Any = None
    whtbltibe12: Any = None
    whtblli4sio4: Any = None
    wtblli2o: Any = None
    vfcblkt: Any = None
    vfpblkt: Any = None
    whtshld: Any = None
    wpenshld: Any = None
    fwmass: Any = None
    fw_armour_vol: Any = None
    fw_armour_thickness: Any = None
    fw_armour_mass: Any = None
    armour_fw_bl_mass: Any = None
    volblkti: Any = None
    volblkto: Any = None
    iblnkith: Any = None
    fblhebmi: Any = None
    fblhebpi: Any = None
    fblhebmo: Any = None
    fblhebpo: Any = None
    fblss: Any = None
    fblbe: Any = None
    whtblbreed: Any = None
    densbreed: Any = None
    fblbreed: Any = None
    iblanket: Any = None
    denw: Any = None
    vffwi: Any = None
    vffwo: Any = None
    volfw: Any = None
    fblss_ccfe: Any = None
    fblli2sio4: Any = None
    fbltibe12: Any = None
    expected_divsur: Any = None
    expected_divmas: Any = None
    expected_whtblbe: Any = None
    expected_whtblss: Any = None
    expected_whtblkt: Any = None
    expected_coolmass: Any = None
    expected_fwclfr: Any = None
    expected_whtbltibe12: Any = None
    expected_whtblli4sio4: Any = None
    expected_wtblli2o: Any = None
    expected_whtshld: Any = None
    expected_wpenshld: Any = None
    expected_fwmass: Any = None
    expected_fw_armour_vol: Any = None
    expected_fw_armour_mass: Any = None
    expected_armour_fw_bl_mass: Any = None
    expected_fblss_ccfe: Any = None
    expected_fblli2sio4: Any = None
    expected_fbltibe12: Any = None


@pytest.mark.parametrize(
    "componentmassesparam",
    (
        ComponentMassesParam(
            divsur=0,
            divclfr=0.29999999999999999,
            divplt=0.035000000000000003,
            fdiva=1.1100000000000001,
            divmas=0,
            divdens=10000,
            rminor=2.6666666666666665,
            rmajor=8,
            idivrt=1,
            sarea=1173.8427771245592,
            blnkith=0.70000000000000007,
            blbuith=0.36499999999999999,
            blbmith=0.17000000000000001,
            blbpith=0.29999999999999999,
            blnkoth=1,
            blbuoth=0.46500000000000002,
            blbmoth=0.27000000000000002,
            blbpoth=0.34999999999999998,
            fwareaib=505.96109565204046,
            fwith=0.018000000000000002,
            fwareaob=838.00728058362097,
            fwoth=0.018000000000000002,
            fwarea=1343.9683762356615,
            volblkt=1182.5433772195902,
            vfblkt=0.25,
            whtblbe=0,
            whtblss=0,
            denstl=7800,
            whtblkt=0,
            volshld=783.69914576548854,
            vfshld=0.60000000000000009,
            coolmass=0,
            fwclfr=0.14999999999999999,
            breeder_f=0.5,
            breeder_multiplier=0.75,
            whtbltibe12=0,
            whtblli4sio4=0,
            wtblli2o=0,
            vfcblkt=0.052949999999999997,
            vfpblkt=0.10000000000000001,
            whtshld=0,
            wpenshld=0,
            fwmass=0,
            fw_armour_vol=0,
            fw_armour_thickness=0.0050000000000000001,
            fw_armour_mass=0,
            armour_fw_bl_mass=0,
            volblkti=315.83946385183026,
            volblkto=866.70391336775992,
            iblnkith=1,
            fblhebmi=0.40000000000000002,
            fblhebpi=0.65949999999999998,
            fblhebmo=0.40000000000000002,
            fblhebpo=0.67130000000000001,
            fblss=0.097049999999999997,
            fblbe=0.59999999999999998,
            whtblbreed=0,
            densbreed=0,
            fblbreed=0.154,
            iblanket=1,
            denw=19250,
            vffwi=0,
            vffwo=0,
            volfw=0,
            fblss_ccfe=0,
            fblli2sio4=0,
            fbltibe12=0,
            expected_divsur=148.78582807401261,
            expected_divmas=36452.527878133093,
            expected_whtblbe=1002205.5121936026,
            expected_whtblss=895173.51112145756,
            expected_whtblkt=2961668.0628126911,
            expected_coolmass=1161.8025382862772,
            expected_fwclfr=0,
            expected_whtbltibe12=1002205.5121936026,
            expected_whtblli4sio4=1064289.0394976311,
            expected_wtblli2o=1064289.0394976311,
            expected_whtshld=2445141.3347883238,
            expected_wpenshld=2445141.3347883238,
            expected_fwmass=188693.16002348688,
            expected_fw_armour_vol=5.8692138856227967,
            expected_fw_armour_mass=112982.36729823884,
            expected_armour_fw_bl_mass=3263343.5901344167,
            expected_fblss_ccfe=0.097049999999999997,
            expected_fblli2sio4=0.375,
            expected_fbltibe12=0.375,
        ),
    ),
)
def test_component_masses(componentmassesparam, monkeypatch, ccfe_hcpb):
    """
    Automatically generated Regression Unit Test for component_masses.

    This test was generated using data from tests/regression/input_files/large_tokamak_once_through.IN.DAT.

    :param componentmassesparam: the data used to mock and assert in this test.
    :type componentmassesparam: componentmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(divertor_variables, "divsur", componentmassesparam.divsur)
    monkeypatch.setattr(divertor_variables, "divclfr", componentmassesparam.divclfr)
    monkeypatch.setattr(divertor_variables, "divplt", componentmassesparam.divplt)
    monkeypatch.setattr(divertor_variables, "fdiva", componentmassesparam.fdiva)
    monkeypatch.setattr(divertor_variables, "divmas", componentmassesparam.divmas)
    monkeypatch.setattr(divertor_variables, "divdens", componentmassesparam.divdens)
    monkeypatch.setattr(physics_variables, "rminor", componentmassesparam.rminor)
    monkeypatch.setattr(physics_variables, "rmajor", componentmassesparam.rmajor)
    monkeypatch.setattr(physics_variables, "idivrt", componentmassesparam.idivrt)
    monkeypatch.setattr(physics_variables, "sarea", componentmassesparam.sarea)
    monkeypatch.setattr(build_variables, "blnkith", componentmassesparam.blnkith)
    monkeypatch.setattr(build_variables, "blbuith", componentmassesparam.blbuith)
    monkeypatch.setattr(build_variables, "blbmith", componentmassesparam.blbmith)
    monkeypatch.setattr(build_variables, "blbpith", componentmassesparam.blbpith)
    monkeypatch.setattr(build_variables, "blnkoth", componentmassesparam.blnkoth)
    monkeypatch.setattr(build_variables, "blbuoth", componentmassesparam.blbuoth)
    monkeypatch.setattr(build_variables, "blbmoth", componentmassesparam.blbmoth)
    monkeypatch.setattr(build_variables, "blbpoth", componentmassesparam.blbpoth)
    monkeypatch.setattr(build_variables, "fwareaib", componentmassesparam.fwareaib)
    monkeypatch.setattr(build_variables, "fwith", componentmassesparam.fwith)
    monkeypatch.setattr(build_variables, "fwareaob", componentmassesparam.fwareaob)
    monkeypatch.setattr(build_variables, "fwoth", componentmassesparam.fwoth)
    monkeypatch.setattr(build_variables, "fwarea", componentmassesparam.fwarea)
    monkeypatch.setattr(fwbs_variables, "volblkt", componentmassesparam.volblkt)
    monkeypatch.setattr(fwbs_variables, "vfblkt", componentmassesparam.vfblkt)
    monkeypatch.setattr(fwbs_variables, "whtblbe", componentmassesparam.whtblbe)
    monkeypatch.setattr(fwbs_variables, "whtblss", componentmassesparam.whtblss)
    monkeypatch.setattr(fwbs_variables, "denstl", componentmassesparam.denstl)
    monkeypatch.setattr(fwbs_variables, "whtblkt", componentmassesparam.whtblkt)
    monkeypatch.setattr(fwbs_variables, "volshld", componentmassesparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vfshld", componentmassesparam.vfshld)
    monkeypatch.setattr(fwbs_variables, "coolmass", componentmassesparam.coolmass)
    monkeypatch.setattr(fwbs_variables, "fwclfr", componentmassesparam.fwclfr)
    monkeypatch.setattr(fwbs_variables, "breeder_f", componentmassesparam.breeder_f)
    monkeypatch.setattr(
        fwbs_variables, "breeder_multiplier", componentmassesparam.breeder_multiplier
    )
    monkeypatch.setattr(fwbs_variables, "whtbltibe12", componentmassesparam.whtbltibe12)
    monkeypatch.setattr(
        fwbs_variables, "whtblli4sio4", componentmassesparam.whtblli4sio4
    )
    monkeypatch.setattr(fwbs_variables, "wtblli2o", componentmassesparam.wtblli2o)
    monkeypatch.setattr(fwbs_variables, "vfcblkt", componentmassesparam.vfcblkt)
    monkeypatch.setattr(fwbs_variables, "vfpblkt", componentmassesparam.vfpblkt)
    monkeypatch.setattr(fwbs_variables, "whtshld", componentmassesparam.whtshld)
    monkeypatch.setattr(fwbs_variables, "wpenshld", componentmassesparam.wpenshld)
    monkeypatch.setattr(fwbs_variables, "fwmass", componentmassesparam.fwmass)
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
    monkeypatch.setattr(fwbs_variables, "volblkti", componentmassesparam.volblkti)
    monkeypatch.setattr(fwbs_variables, "volblkto", componentmassesparam.volblkto)
    monkeypatch.setattr(fwbs_variables, "iblnkith", componentmassesparam.iblnkith)
    monkeypatch.setattr(fwbs_variables, "fblhebmi", componentmassesparam.fblhebmi)
    monkeypatch.setattr(fwbs_variables, "fblhebpi", componentmassesparam.fblhebpi)
    monkeypatch.setattr(fwbs_variables, "fblhebmo", componentmassesparam.fblhebmo)
    monkeypatch.setattr(fwbs_variables, "fblhebpo", componentmassesparam.fblhebpo)
    monkeypatch.setattr(fwbs_variables, "fblss", componentmassesparam.fblss)
    monkeypatch.setattr(fwbs_variables, "fblbe", componentmassesparam.fblbe)
    monkeypatch.setattr(fwbs_variables, "whtblbreed", componentmassesparam.whtblbreed)
    monkeypatch.setattr(fwbs_variables, "densbreed", componentmassesparam.densbreed)
    monkeypatch.setattr(fwbs_variables, "fblbreed", componentmassesparam.fblbreed)
    monkeypatch.setattr(fwbs_variables, "iblanket", componentmassesparam.iblanket)
    monkeypatch.setattr(fwbs_variables, "denw", componentmassesparam.denw)
    monkeypatch.setattr(fwbs_variables, "vffwi", componentmassesparam.vffwi)
    monkeypatch.setattr(fwbs_variables, "vffwo", componentmassesparam.vffwo)
    monkeypatch.setattr(fwbs_variables, "volfw", componentmassesparam.volfw)
    monkeypatch.setattr(fwbs_variables, "fblss_ccfe", componentmassesparam.fblss_ccfe)
    monkeypatch.setattr(fwbs_variables, "fblli2sio4", componentmassesparam.fblli2sio4)
    monkeypatch.setattr(fwbs_variables, "fbltibe12", componentmassesparam.fbltibe12)

    ccfe_hcpb.component_masses()

    assert divertor_variables.divsur == pytest.approx(
        componentmassesparam.expected_divsur
    )
    assert divertor_variables.divmas == pytest.approx(
        componentmassesparam.expected_divmas
    )
    assert fwbs_variables.whtblbe == pytest.approx(
        componentmassesparam.expected_whtblbe
    )
    assert fwbs_variables.whtblss == pytest.approx(
        componentmassesparam.expected_whtblss
    )
    assert fwbs_variables.whtblkt == pytest.approx(
        componentmassesparam.expected_whtblkt
    )
    assert fwbs_variables.coolmass == pytest.approx(
        componentmassesparam.expected_coolmass
    )
    assert fwbs_variables.fwclfr == pytest.approx(componentmassesparam.expected_fwclfr)
    assert fwbs_variables.whtbltibe12 == pytest.approx(
        componentmassesparam.expected_whtbltibe12
    )
    assert fwbs_variables.whtblli4sio4 == pytest.approx(
        componentmassesparam.expected_whtblli4sio4
    )
    assert fwbs_variables.wtblli2o == pytest.approx(
        componentmassesparam.expected_wtblli2o
    )
    assert fwbs_variables.whtshld == pytest.approx(
        componentmassesparam.expected_whtshld
    )
    assert fwbs_variables.wpenshld == pytest.approx(
        componentmassesparam.expected_wpenshld
    )
    assert fwbs_variables.fwmass == pytest.approx(componentmassesparam.expected_fwmass)
    assert fwbs_variables.fw_armour_vol == pytest.approx(
        componentmassesparam.expected_fw_armour_vol
    )
    assert fwbs_variables.fw_armour_mass == pytest.approx(
        componentmassesparam.expected_fw_armour_mass
    )
    assert fwbs_variables.armour_fw_bl_mass == pytest.approx(
        componentmassesparam.expected_armour_fw_bl_mass
    )
    assert fwbs_variables.fblss_ccfe == pytest.approx(
        componentmassesparam.expected_fblss_ccfe
    )
    assert fwbs_variables.fblli2sio4 == pytest.approx(
        componentmassesparam.expected_fblli2sio4
    )
    assert fwbs_variables.fbltibe12 == pytest.approx(
        componentmassesparam.expected_fbltibe12
    )
