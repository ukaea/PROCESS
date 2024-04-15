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

    fw_armour_thicknessi: Any = None

    fw_armour_thicknesso: Any = None

    ptfnuc: Any = None

    denw: Any = None

    vffwi: Any = None

    vffwo: Any = None

    powfmw: Any = None

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
            blnkith=0.70000000000000000,
            blnkoth=1.00000000000000000,
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            afw=0.0060000000000000001,
            pitch=0.02,
            denstl=7800,
            whtblkt=2967600.0000000000,
            volblkt=1184.900112000000,
            whtshld=2449800.000000000,
            volshld=785.2000000000000,
            vvmass=7938800.0000000000,
            vdewin=1017.8000000000000,
            fw_armour_thicknessi=0.0030000000000000000,
            fw_armour_thicknesso=0.0030000000000000000,
            ptfnuc=0,
            denw=19250,
            vffwi=0,
            vffwo=0,
            powfmw=1663.1000000000000,
            itart=0,
            whttf=10065000.000000000,
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
            expected_ptfnuc=0.02347076094771759,
            expected_vffwi=0.31415926535897931,
            expected_vffwo=0.31415926535897931,
            expected_armour_density=13202.434141839649,
            expected_fw_density=5349.557730199961,
            expected_blanket_density=2504.5149122240946,
            expected_shield_density=3119.9694345389707,
            expected_vv_density=7799.960699548045,
            expected_x_blanket=2.264737016959599,
            expected_x_shield=4.055971398860848,
            expected_tfc_nuc_heating=14112.657656014428,
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
            fw_armour_thicknessi=0.0050000000000000001,
            fw_armour_thicknesso=0.0050000000000000001,
            ptfnuc=0.044184461825198453,
            denw=19250,
            vffwi=0.31415926535897931,
            vffwo=0.31415926535897931,
            powfmw=1985.4423932312809,
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
        "fw_armour_thicknessi",
        nuclearheatingmagnetsparam.fw_armour_thicknessi,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "fw_armour_thicknesso",
        nuclearheatingmagnetsparam.fw_armour_thicknesso,
    )

    monkeypatch.setattr(fwbs_variables, "ptfnuc", nuclearheatingmagnetsparam.ptfnuc)

    monkeypatch.setattr(fwbs_variables, "denw", nuclearheatingmagnetsparam.denw)

    monkeypatch.setattr(fwbs_variables, "vffwi", nuclearheatingmagnetsparam.vffwi)

    monkeypatch.setattr(fwbs_variables, "vffwo", nuclearheatingmagnetsparam.vffwo)

    monkeypatch.setattr(physics_variables, "powfmw", nuclearheatingmagnetsparam.powfmw)

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
    pnucfw: Any = None

    fwmass: Any = None

    powfmw: Any = None

    fw_armour_u_nuc_heating: Any = None

    expected_pnucfw: Any = None

    expected_fw_armour_u_nuc_heating: Any = None


@pytest.mark.parametrize(
    "nuclearheatingfwparam",
    (
        NuclearHeatingFwParam(
            pnucfw=0,
            fwmass=224802.80270851994,
            powfmw=1986.0623241661431,
            fw_armour_u_nuc_heating=0,
            expected_pnucfw=279.04523551646628,
            expected_fw_armour_u_nuc_heating=6.2500000000000005e-07,
        ),
        NuclearHeatingFwParam(
            pnucfw=276.80690153753221,
            fwmass=182115.83467868491,
            powfmw=1985.4423932312809,
            fw_armour_u_nuc_heating=6.2500000000000005e-07,
            expected_pnucfw=225.98781165610032,
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

    monkeypatch.setattr(fwbs_variables, "pnucfw", nuclearheatingfwparam.pnucfw)

    monkeypatch.setattr(fwbs_variables, "fwmass", nuclearheatingfwparam.fwmass)

    monkeypatch.setattr(physics_variables, "powfmw", nuclearheatingfwparam.powfmw)

    monkeypatch.setattr(
        ccfe_hcpb_module,
        "fw_armour_u_nuc_heating",
        nuclearheatingfwparam.fw_armour_u_nuc_heating,
    )

    ccfe_hcpb.nuclear_heating_fw()

    assert fwbs_variables.pnucfw == pytest.approx(nuclearheatingfwparam.expected_pnucfw)

    assert ccfe_hcpb_module.fw_armour_u_nuc_heating == pytest.approx(
        nuclearheatingfwparam.expected_fw_armour_u_nuc_heating
    )


class NuclearHeatingBlanketParam(NamedTuple):
    whtblkt: Any = None

    pnucblkt: Any = None

    powfmw: Any = None

    exp_blanket: Any = None

    expected_pnucblkt: Any = None

    expected_exp_blanket: Any = None


@pytest.mark.parametrize(
    "nuclearheatingblanketparam",
    (
        NuclearHeatingBlanketParam(
            whtblkt=3501027.3252278985,
            pnucblkt=0,
            powfmw=1986.0623241661431,
            exp_blanket=0,
            expected_pnucblkt=1517.0907688379014,
            expected_exp_blanket=0.99982809071915879,
        ),
        NuclearHeatingBlanketParam(
            whtblkt=3507503.3737008357,
            pnucblkt=1504.9215740808861,
            powfmw=1985.4423932312809,
            exp_blanket=0.99982809071915879,
            expected_pnucblkt=1516.6213709741428,
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

    monkeypatch.setattr(fwbs_variables, "pnucblkt", nuclearheatingblanketparam.pnucblkt)

    monkeypatch.setattr(physics_variables, "powfmw", nuclearheatingblanketparam.powfmw)

    monkeypatch.setattr(
        ccfe_hcpb_module, "exp_blanket", nuclearheatingblanketparam.exp_blanket
    )

    ccfe_hcpb.nuclear_heating_blanket()

    assert fwbs_variables.pnucblkt == pytest.approx(
        nuclearheatingblanketparam.expected_pnucblkt
    )

    assert ccfe_hcpb_module.exp_blanket == pytest.approx(
        nuclearheatingblanketparam.expected_exp_blanket
    )


class NuclearHeatingShieldParam(NamedTuple):
    shldith: Any = None

    shldoth: Any = None

    whtshld: Any = None

    pnucshld: Any = None

    powfmw: Any = None

    itart: Any = None

    shield_density: Any = None

    x_blanket: Any = None

    shld_u_nuc_heating: Any = None

    exp_shield1: Any = None

    exp_shield2: Any = None

    expected_pnucshld: Any = None

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
            pnucshld=0,
            powfmw=1986.0623241661431,
            itart=0,
            shield_density=3119.9999999999995,
            x_blanket=2.3374537748527975,
            shld_u_nuc_heating=0,
            exp_shield1=0,
            exp_shield2=0,
            expected_pnucshld=1.3721323841005297,
            expected_shld_u_nuc_heating=690880.82856444374,
            expected_exp_shield1=0.0017209365527675318,
            expected_exp_shield2=0.25426760591013942,
        ),
        NuclearHeatingShieldParam(
            shldith=0.30000000000000004,
            shldoth=0.80000000000000004,
            whtshld=2297808.3935174854,
            pnucshld=1.3611259588044891,
            powfmw=1985.4423932312809,
            itart=0,
            shield_density=3120,
            x_blanket=2.3374537748527979,
            shld_u_nuc_heating=690880.82856444374,
            exp_shield1=0.0017209365527675318,
            exp_shield2=0.25426760591013942,
            expected_pnucshld=1.3734581585671393,
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

    monkeypatch.setattr(fwbs_variables, "pnucshld", nuclearheatingshieldparam.pnucshld)

    monkeypatch.setattr(physics_variables, "powfmw", nuclearheatingshieldparam.powfmw)

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

    assert fwbs_variables.pnucshld == pytest.approx(
        nuclearheatingshieldparam.expected_pnucshld
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

    powfmw: Any = None

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
            powfmw=1986.0623241661431,
            ip=0,
            expected_pnucdiv=182.71773382328519,
        ),
        NuclearHeatingDivertorParam(
            fdiv=0.115,
            pnucdiv=182.71773382328519,
            pnuchcd=0,
            idivrt=1,
            powfmw=1985.4423932312809,
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

    monkeypatch.setattr(physics_variables, "powfmw", nuclearheatingdivertorparam.powfmw)

    monkeypatch.setattr(ccfe_hcpb_module, "ip", nuclearheatingdivertorparam.ip)

    ccfe_hcpb.nuclear_heating_divertor()

    assert fwbs_variables.pnucdiv == pytest.approx(
        nuclearheatingdivertorparam.expected_pnucdiv
    )


class PowerflowCalcParam(NamedTuple):
    fwareaob: Any = None

    fwarea: Any = None

    porbitlossmw: Any = None

    fdiv: Any = None

    praddiv: Any = None

    pradhcd: Any = None

    fhcd: Any = None

    pradfw: Any = None

    coolwh: Any = None

    outlet_temp: Any = None

    blpressure: Any = None

    primary_pumping: Any = None

    pnucfw: Any = None

    pnucblkt: Any = None

    pnucdiv: Any = None

    pnucshld: Any = None

    etaiso: Any = None

    pnuc_cp_sh: Any = None

    psurffwi: Any = None

    psurffwo: Any = None

    htpmw_fw: Any = None

    fpumpfw: Any = None

    htpmw_blkt: Any = None

    fpumpblkt: Any = None

    htpmw_shld: Any = None

    fpumpshld: Any = None

    htpmw_div: Any = None

    fpumpdiv: Any = None

    idivrt: Any = None

    pradmw: Any = None

    palpfwmw: Any = None

    pdivt: Any = None

    p_he: Any = None

    dp_he: Any = None

    gamma_he: Any = None

    t_in_bb: Any = None

    t_out_bb: Any = None

    htpmw_fw_blkt: Any = None

    ip: Any = None

    ofile: Any = None

    expected_praddiv: Any = None

    expected_pradfw: Any = None

    expected_psurffwi: Any = None

    expected_psurffwo: Any = None

    expected_htpmw_shld: Any = None

    expected_htpmw_div: Any = None

    expected_htpmw_fw_blkt: Any = None


@pytest.mark.parametrize(
    "powerflowcalcparam",
    (
        PowerflowCalcParam(
            fwareaob=988.92586580655245,
            fwarea=1601.1595634509963,
            porbitlossmw=0,
            fdiv=0.115,
            praddiv=0,
            pradhcd=0,
            fhcd=0,
            pradfw=0,
            coolwh=1,
            outlet_temp=823,
            blpressure=15500000,
            primary_pumping=3,
            pnucfw=276.80690153753221,
            pnucblkt=1504.9215740808861,
            pnucdiv=182.71773382328519,
            pnucshld=1.3611259588044891,
            etaiso=0.90000000000000002,
            pnuc_cp_sh=0,
            psurffwi=0,
            psurffwo=0,
            htpmw_fw=0,
            fpumpfw=0.0050000000000000001,
            htpmw_blkt=0,
            fpumpblkt=0.0050000000000000001,
            htpmw_shld=0,
            fpumpshld=0.0050000000000000001,
            htpmw_div=0,
            fpumpdiv=0.0050000000000000001,
            idivrt=1,
            pradmw=287.44866938104849,
            palpfwmw=19.835845058655043,
            pdivt=143.6315222649435,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            htpmw_fw_blkt=0,
            ip=0,
            ofile=11,
            expected_praddiv=33.056596978820579,
            expected_pradfw=254.39207240222791,
            expected_psurffwi=97.271629070225231,
            expected_psurffwo=176.95628839065773,
            expected_htpmw_shld=0.0068056297940224456,
            expected_htpmw_div=1.7970292653352464,
            expected_htpmw_fw_blkt=202.00455086503842,
        ),
        PowerflowCalcParam(
            fwareaob=1168.1172772224481,
            fwarea=1891.2865102700493,
            porbitlossmw=0,
            fdiv=0.115,
            praddiv=33.056596978820579,
            pradhcd=0,
            fhcd=0,
            pradfw=254.39207240222791,
            coolwh=1,
            outlet_temp=823,
            blpressure=15500000,
            primary_pumping=3,
            pnucfw=230.98304919926957,
            pnucblkt=1550.1447895848396,
            pnucdiv=182.66070017727785,
            pnucshld=1.4038170956592293,
            etaiso=0.90000000000000002,
            pnuc_cp_sh=0,
            psurffwi=97.271629070225231,
            psurffwo=176.95628839065773,
            htpmw_fw=0,
            fpumpfw=0.0050000000000000001,
            htpmw_blkt=0,
            fpumpblkt=0.0050000000000000001,
            htpmw_shld=0.0068056297940224456,
            fpumpshld=0.0050000000000000001,
            htpmw_div=1.7970292653352464,
            fpumpdiv=0.0050000000000000001,
            idivrt=1,
            pradmw=287.44866938104849,
            palpfwmw=19.829653483586444,
            pdivt=143.51338080047339,
            p_he=8000000,
            dp_he=550000,
            gamma_he=1.667,
            t_in_bb=573.13,
            t_out_bb=773.13,
            htpmw_fw_blkt=202.00455086503842,
            ip=0,
            ofile=11,
            expected_praddiv=33.056596978820579,
            expected_pradfw=254.39207240222791,
            expected_psurffwi=97.271629070225259,
            expected_psurffwo=176.95009681558912,
            expected_htpmw_shld=0.007019085478296147,
            expected_htpmw_div=1.7961533897828594,
            expected_htpmw_fw_blkt=201.94492795635171,
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
        current_drive_variables, "porbitlossmw", powerflowcalcparam.porbitlossmw
    )

    monkeypatch.setattr(fwbs_variables, "fdiv", powerflowcalcparam.fdiv)

    monkeypatch.setattr(fwbs_variables, "praddiv", powerflowcalcparam.praddiv)

    monkeypatch.setattr(fwbs_variables, "pradhcd", powerflowcalcparam.pradhcd)

    monkeypatch.setattr(fwbs_variables, "fhcd", powerflowcalcparam.fhcd)

    monkeypatch.setattr(fwbs_variables, "pradfw", powerflowcalcparam.pradfw)

    monkeypatch.setattr(fwbs_variables, "coolwh", powerflowcalcparam.coolwh)

    monkeypatch.setattr(fwbs_variables, "outlet_temp", powerflowcalcparam.outlet_temp)

    monkeypatch.setattr(fwbs_variables, "blpressure", powerflowcalcparam.blpressure)

    monkeypatch.setattr(
        fwbs_variables, "primary_pumping", powerflowcalcparam.primary_pumping
    )

    monkeypatch.setattr(fwbs_variables, "pnucfw", powerflowcalcparam.pnucfw)

    monkeypatch.setattr(fwbs_variables, "pnucblkt", powerflowcalcparam.pnucblkt)

    monkeypatch.setattr(fwbs_variables, "pnucdiv", powerflowcalcparam.pnucdiv)

    monkeypatch.setattr(fwbs_variables, "pnucshld", powerflowcalcparam.pnucshld)

    monkeypatch.setattr(fwbs_variables, "etaiso", powerflowcalcparam.etaiso)

    monkeypatch.setattr(fwbs_variables, "pnuc_cp_sh", powerflowcalcparam.pnuc_cp_sh)

    monkeypatch.setattr(fwbs_variables, "psurffwi", powerflowcalcparam.psurffwi)

    monkeypatch.setattr(fwbs_variables, "psurffwo", powerflowcalcparam.psurffwo)

    monkeypatch.setattr(
        heat_transport_variables, "htpmw_fw", powerflowcalcparam.htpmw_fw
    )

    monkeypatch.setattr(heat_transport_variables, "fpumpfw", powerflowcalcparam.fpumpfw)

    monkeypatch.setattr(
        heat_transport_variables, "htpmw_blkt", powerflowcalcparam.htpmw_blkt
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpblkt", powerflowcalcparam.fpumpblkt
    )

    monkeypatch.setattr(
        heat_transport_variables, "htpmw_shld", powerflowcalcparam.htpmw_shld
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpshld", powerflowcalcparam.fpumpshld
    )

    monkeypatch.setattr(
        heat_transport_variables, "htpmw_div", powerflowcalcparam.htpmw_div
    )

    monkeypatch.setattr(
        heat_transport_variables, "fpumpdiv", powerflowcalcparam.fpumpdiv
    )

    monkeypatch.setattr(physics_variables, "idivrt", powerflowcalcparam.idivrt)

    monkeypatch.setattr(physics_variables, "pradmw", powerflowcalcparam.pradmw)

    monkeypatch.setattr(physics_variables, "palpfwmw", powerflowcalcparam.palpfwmw)

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
        primary_pumping_variables, "htpmw_fw_blkt", powerflowcalcparam.htpmw_fw_blkt
    )

    monkeypatch.setattr(ccfe_hcpb_module, "ip", powerflowcalcparam.ip)

    monkeypatch.setattr(ccfe_hcpb_module, "ofile", powerflowcalcparam.ofile)

    ccfe_hcpb.powerflow_calc(False)

    assert fwbs_variables.praddiv == pytest.approx(powerflowcalcparam.expected_praddiv)

    assert fwbs_variables.pradfw == pytest.approx(powerflowcalcparam.expected_pradfw)

    assert fwbs_variables.psurffwi == pytest.approx(
        powerflowcalcparam.expected_psurffwi
    )

    assert fwbs_variables.psurffwo == pytest.approx(
        powerflowcalcparam.expected_psurffwo
    )

    assert heat_transport_variables.htpmw_shld == pytest.approx(
        powerflowcalcparam.expected_htpmw_shld
    )

    assert heat_transport_variables.htpmw_div == pytest.approx(
        powerflowcalcparam.expected_htpmw_div
    )

    assert primary_pumping_variables.htpmw_fw_blkt == pytest.approx(
        powerflowcalcparam.expected_htpmw_fw_blkt
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

    pneutmw: Any = None

    sh_width: Any = None

    rmajor: Any = None

    expected_neut_flux_cp: Any = None


@pytest.mark.parametrize(
    "sttfcentrepostfastneutfluxparam",
    (
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            pneutmw=400.65875490746737,
            sh_width=0.60000000000000009,
            rmajor=3,
            expected_neut_flux_cp=144701998710998.5,
        ),
        StTfCentrepostFastNeutFluxParam(
            i_tf_sup=1,
            pneutmw=409.82485143909827,
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
        pneutmw=sttfcentrepostfastneutfluxparam.pneutmw,
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

    fw_armour_thicknessi: Any = None

    fw_armour_thicknesso: Any = None

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
            fw_armour_thicknessi=0.0030000000000000001,
            fw_armour_thicknesso=0.0030000000000000001,
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
        fwbs_variables, "fw_armour_thicknessi", tbrshimwellparam.fw_armour_thicknessi
    )

    monkeypatch.setattr(
        fwbs_variables, "fw_armour_thicknesso", tbrshimwellparam.fw_armour_thicknesso
    )

    monkeypatch.setattr(ccfe_hcpb_module, "ip", tbrshimwellparam.ip)

    tbr = ccfe_hcpb.tbr_shimwell(
        iblanket_thickness=tbrshimwellparam.iblanket_thickness,
        breeder_f=tbrshimwellparam.breeder_f,
        li6enrich=tbrshimwellparam.li6enrich,
        output=False,
    )

    assert tbr == pytest.approx(tbrshimwellparam.expected_tbr)
