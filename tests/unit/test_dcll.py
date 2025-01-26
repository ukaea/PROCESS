from typing import Any, NamedTuple

import pytest

from process.blanket_library import BlanketLibrary
from process.dcll import DCLL
from process.fortran import (
    build_variables,
    current_drive_variables,
    dcll_module,
    fwbs_variables,
    physics_variables,
)
from process.fw import Fw


@pytest.fixture
def dcll():
    """Provides DCLL object for testing.

    :returns: initialised DCLL object
    :rtype: process.dcll.DCLL
    """
    return DCLL(BlanketLibrary(Fw()))


class DcllNeutronicsAndPowerParam(NamedTuple):
    fwareaob: Any = None

    fwarea: Any = None

    porbitlossmw: Any = None

    fdiv: Any = None

    praddiv: Any = None

    pnucdiv: Any = None

    fhcd: Any = None

    pradhcd: Any = None

    pnuchcd: Any = None

    pnucshld: Any = None

    pradfw: Any = None

    pnucfw: Any = None

    psurffwi: Any = None

    psurffwo: Any = None

    pnucblkt: Any = None

    pnuc_fw_ratio_dcll: Any = None

    pnuc_blkt_ratio_dcll: Any = None

    emult: Any = None

    emultmw: Any = None

    ptfnuc: Any = None

    idivrt: Any = None

    neutron_power_total: Any = None

    p_plasma_rad_mw: Any = None

    palpfwmw: Any = None

    expected_praddiv: Any = None

    expected_pnucdiv: Any = None

    expected_pradfw: Any = None

    expected_pnucfw: Any = None

    expected_pnucblkt: Any = None

    expected_emultmw: Any = None


@pytest.mark.parametrize(
    "dcllneutronicsandpowerparam",
    (
        DcllNeutronicsAndPowerParam(
            fwareaob=988.92586580655245,
            fwarea=1601.1595634509963,
            porbitlossmw=0,
            fdiv=0.115,
            praddiv=0,
            pnucdiv=0,
            fhcd=0,
            pradhcd=0,
            pnuchcd=0,
            pnucshld=0,
            pradfw=0,
            pnucfw=0,
            psurffwi=0,
            psurffwo=0,
            pnucblkt=0,
            pnuc_fw_ratio_dcll=0.14000000000000001,
            pnuc_blkt_ratio_dcll=0.85999999999999999,
            emult=1.2689999999999999,
            emultmw=0,
            ptfnuc=0,
            idivrt=1,
            neutron_power_total=1587.7386535917431,
            p_plasma_rad_mw=287.44866938104849,
            palpfwmw=19.835845058655043,
            expected_praddiv=33.056596978820579,
            expected_pnucdiv=182.58994516305046,
            expected_pradfw=254.39207240222791,
            expected_pnucfw=196.72081918001697,
            expected_pnucblkt=1533.4949914565693,
            expected_emultmw=325.06710220789364,
        ),
        DcllNeutronicsAndPowerParam(
            fwareaob=1168.1172772224481,
            fwarea=1891.2865102700493,
            porbitlossmw=0,
            fdiv=0.115,
            praddiv=33.056596978820579,
            pnucdiv=182.58994516305046,
            fhcd=0,
            pradhcd=0,
            pnuchcd=0,
            pnucshld=0,
            pradfw=254.39207240222791,
            pnucfw=196.72081918001697,
            psurffwi=97.271629070225231,
            psurffwo=176.95628839065773,
            pnucblkt=1533.4949914565693,
            pnuc_fw_ratio_dcll=0.14000000000000001,
            pnuc_blkt_ratio_dcll=0.85999999999999999,
            emult=1.2689999999999999,
            emultmw=325.06710220789364,
            ptfnuc=0,
            idivrt=1,
            neutron_power_total=1587.2430556964196,
            p_plasma_rad_mw=287.44866938104849,
            palpfwmw=19.829653483586444,
            expected_praddiv=33.056596978820579,
            expected_pnucdiv=182.53295140508826,
            expected_pradfw=254.39207240222791,
            expected_pnucfw=196.65941460078642,
            expected_pnucblkt=1533.0163252173013,
            expected_emultmw=324.96563552675644,
        ),
    ),
)
def test_dcll_neutronics_and_power(dcllneutronicsandpowerparam, monkeypatch, dcll):
    """
    Automatically generated Regression Unit Test for dcll_neutronics_and_power.

    This test was generated using data from ./dcll/dcll_mms_lt_IN.DAT.

    :param dcllneutronicsandpowerparam: the data used to mock and assert in this test.
    :type dcllneutronicsandpowerparam: dcllneutronicsandpowerparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        build_variables, "fwareaob", dcllneutronicsandpowerparam.fwareaob
    )

    monkeypatch.setattr(build_variables, "fwarea", dcllneutronicsandpowerparam.fwarea)

    monkeypatch.setattr(
        current_drive_variables,
        "porbitlossmw",
        dcllneutronicsandpowerparam.porbitlossmw,
    )

    monkeypatch.setattr(fwbs_variables, "fdiv", dcllneutronicsandpowerparam.fdiv)

    monkeypatch.setattr(fwbs_variables, "praddiv", dcllneutronicsandpowerparam.praddiv)

    monkeypatch.setattr(fwbs_variables, "pnucdiv", dcllneutronicsandpowerparam.pnucdiv)

    monkeypatch.setattr(fwbs_variables, "fhcd", dcllneutronicsandpowerparam.fhcd)

    monkeypatch.setattr(fwbs_variables, "pradhcd", dcllneutronicsandpowerparam.pradhcd)

    monkeypatch.setattr(fwbs_variables, "pnuchcd", dcllneutronicsandpowerparam.pnuchcd)

    monkeypatch.setattr(
        fwbs_variables, "pnucshld", dcllneutronicsandpowerparam.pnucshld
    )

    monkeypatch.setattr(fwbs_variables, "pradfw", dcllneutronicsandpowerparam.pradfw)

    monkeypatch.setattr(fwbs_variables, "pnucfw", dcllneutronicsandpowerparam.pnucfw)

    monkeypatch.setattr(
        fwbs_variables, "psurffwi", dcllneutronicsandpowerparam.psurffwi
    )

    monkeypatch.setattr(
        fwbs_variables, "psurffwo", dcllneutronicsandpowerparam.psurffwo
    )

    monkeypatch.setattr(
        fwbs_variables, "pnucblkt", dcllneutronicsandpowerparam.pnucblkt
    )

    monkeypatch.setattr(
        fwbs_variables,
        "pnuc_fw_ratio_dcll",
        dcllneutronicsandpowerparam.pnuc_fw_ratio_dcll,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "pnuc_blkt_ratio_dcll",
        dcllneutronicsandpowerparam.pnuc_blkt_ratio_dcll,
    )

    monkeypatch.setattr(fwbs_variables, "emult", dcllneutronicsandpowerparam.emult)

    monkeypatch.setattr(fwbs_variables, "emultmw", dcllneutronicsandpowerparam.emultmw)

    monkeypatch.setattr(fwbs_variables, "ptfnuc", dcllneutronicsandpowerparam.ptfnuc)

    monkeypatch.setattr(physics_variables, "idivrt", dcllneutronicsandpowerparam.idivrt)

    monkeypatch.setattr(
        physics_variables,
        "neutron_power_total",
        dcllneutronicsandpowerparam.neutron_power_total,
    )

    monkeypatch.setattr(
        physics_variables,
        "p_plasma_rad_mw",
        dcllneutronicsandpowerparam.p_plasma_rad_mw,
    )

    monkeypatch.setattr(
        physics_variables, "palpfwmw", dcllneutronicsandpowerparam.palpfwmw
    )

    dcll.dcll_neutronics_and_power(False)

    assert fwbs_variables.praddiv == pytest.approx(
        dcllneutronicsandpowerparam.expected_praddiv
    )

    assert fwbs_variables.pnucdiv == pytest.approx(
        dcllneutronicsandpowerparam.expected_pnucdiv
    )

    assert fwbs_variables.pradfw == pytest.approx(
        dcllneutronicsandpowerparam.expected_pradfw
    )

    assert fwbs_variables.pnucfw == pytest.approx(
        dcllneutronicsandpowerparam.expected_pnucfw
    )

    assert fwbs_variables.pnucblkt == pytest.approx(
        dcllneutronicsandpowerparam.expected_pnucblkt
    )

    assert fwbs_variables.emultmw == pytest.approx(
        dcllneutronicsandpowerparam.expected_emultmw
    )


class DcllMassesParam(NamedTuple):
    fwareaib: Any = None

    dr_fw_inboard: Any = None

    fwareaob: Any = None

    dr_fw_outboard: Any = None

    dr_blkt_inboard: Any = None

    blbuith: Any = None

    blbmith: Any = None

    dr_blkt_outboard: Any = None

    blbuoth: Any = None

    blbmoth: Any = None

    a_plasma_surface: Any = None

    a_plasma_surface_outboard: Any = None

    i_blkt_inboard: Any = None

    volblkt: Any = None

    volblkti: Any = None

    volblkto: Any = None

    whtblkt: Any = None

    fwmass: Any = None

    fw_armour_vol: Any = None

    fw_armour_thickness: Any = None

    fw_armour_mass: Any = None

    volfw: Any = None

    armour_fw_bl_mass: Any = None

    denstl: Any = None

    den_liq: Any = None

    ifci: Any = None

    den_ceramic: Any = None

    th_wall_secondary: Any = None

    nopol: Any = None

    r_f_liq_ib: Any = None

    r_f_liq_ob: Any = None

    w_f_liq_ib: Any = None

    w_f_liq_ob: Any = None

    vfblkt: Any = None

    icooldual: Any = None

    rhof_fw: Any = None

    rhof_bl: Any = None

    nblktmodti: Any = None

    nblktmodto: Any = None

    r_fci: Any = None

    r_backwall: Any = None

    bz_r_ib: Any = None

    bz_r_ob: Any = None

    f_vol_stff_plates: Any = None

    f_vol_stl_bz_struct: Any = None

    f_vol_stl_back_wall: Any = None

    f_vol_stl_fw: Any = None

    f_vol_mfbss_stl: Any = None

    f_vol_mfbss_he: Any = None

    f_vol_mfbss_pbli: Any = None

    vol_fci: Any = None

    vol_bz_struct: Any = None

    vol_bz_liq: Any = None

    vol_bz_liq_ib: Any = None

    vol_bz_liq_ob: Any = None

    vol_bw: Any = None

    vol_bss: Any = None

    wht_cer: Any = None

    wht_stl_struct: Any = None

    wht_cool_struct: Any = None

    wht_bw_stl: Any = None

    wht_bw_cool: Any = None

    wht_mfbss_stl: Any = None

    wht_mfbss_cool: Any = None

    wht_mfbss_pbli: Any = None

    fwmass_stl: Any = None

    fwmass_cool: Any = None

    mass_cool_blanket: Any = None

    mass_liq_blanket: Any = None

    mass_stl_blanket: Any = None

    mass_segm_ib: Any = None

    mass_segm_ob: Any = None

    expected_blbmith: Any = None

    expected_blbmoth: Any = None

    expected_whtblkt: Any = None

    expected_fwmass: Any = None

    expected_fw_armour_vol: Any = None

    expected_fw_armour_mass: Any = None

    expected_volfw: Any = None

    expected_armour_fw_bl_mass: Any = None

    expected_r_f_liq_ib: Any = None

    expected_w_f_liq_ib: Any = None

    expected_w_f_liq_ob: Any = None

    expected_vfblkt: Any = None

    expected_r_fci: Any = None

    expected_r_backwall: Any = None

    expected_bz_r_ib: Any = None

    expected_bz_r_ob: Any = None

    expected_f_vol_stff_plates: Any = None

    expected_f_vol_stl_bz_struct: Any = None

    expected_f_vol_stl_back_wall: Any = None

    expected_f_vol_stl_fw: Any = None

    expected_f_vol_mfbss_stl: Any = None

    expected_f_vol_mfbss_he: Any = None

    expected_f_vol_mfbss_pbli: Any = None

    expected_vol_fci: Any = None

    expected_vol_bz_struct: Any = None

    expected_vol_bz_liq: Any = None

    expected_vol_bz_liq_ib: Any = None

    expected_vol_bz_liq_ob: Any = None

    expected_vol_bw: Any = None

    expected_vol_bss: Any = None

    expected_wht_cer: Any = None

    expected_wht_stl_struct: Any = None

    expected_wht_cool_struct: Any = None

    expected_wht_bw_stl: Any = None

    expected_wht_bw_cool: Any = None

    expected_wht_mfbss_stl: Any = None

    expected_wht_mfbss_cool: Any = None

    expected_wht_mfbss_pbli: Any = None

    expected_fwmass_stl: Any = None

    expected_fwmass_cool: Any = None

    expected_mass_cool_blanket: Any = None

    expected_mass_liq_blanket: Any = None

    expected_mass_stl_blanket: Any = None

    expected_mass_segm_ib: Any = None

    expected_mass_segm_ob: Any = None


@pytest.mark.parametrize(
    "dcllmassesparam",
    (
        DcllMassesParam(
            fwareaib=612.23369764444396,
            dr_fw_inboard=0.018000000000000002,
            fwareaob=988.92586580655245,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            blbuith=0.36499999999999999,
            blbmith=0.17000000000000001,
            dr_blkt_outboard=0.98199999999999998,
            blbuoth=0.46500000000000002,
            blbmoth=0.27000000000000002,
            a_plasma_surface=1403.2719775669307,
            a_plasma_surface_outboard=949.22962703393853,
            i_blkt_inboard=1,
            volblkt=1397.9003011502937,
            volblkti=401.90579863726225,
            volblkto=995.99450251303142,
            whtblkt=0,
            fwmass=0,
            fw_armour_vol=0,
            fw_armour_thickness=0.0050000000000000001,
            fw_armour_mass=0,
            volfw=0,
            armour_fw_bl_mass=0,
            denstl=7800,
            den_liq=9753.2497999999996,
            ifci=1,
            den_ceramic=3210,
            th_wall_secondary=0.012500000000000001,
            nopol=2,
            r_f_liq_ib=0.5,
            r_f_liq_ob=0.5,
            w_f_liq_ib=0.5,
            w_f_liq_ob=0.5,
            vfblkt=0.25,
            icooldual=2,
            rhof_fw=5.6389735407435868,
            rhof_bl=5.6389735407435868,
            nblktmodti=32,
            nblktmodto=48,
            r_fci=0,
            r_backwall=0,
            bz_r_ib=0,
            bz_r_ob=0,
            f_vol_stff_plates=0,
            f_vol_stl_bz_struct=0,
            f_vol_stl_back_wall=0,
            f_vol_stl_fw=0,
            f_vol_mfbss_stl=0,
            f_vol_mfbss_he=0,
            f_vol_mfbss_pbli=0,
            vol_fci=0,
            vol_bz_struct=0,
            vol_bz_liq=0,
            vol_bz_liq_ib=0,
            vol_bz_liq_ob=0,
            vol_bw=0,
            vol_bss=0,
            wht_cer=0,
            wht_stl_struct=0,
            wht_cool_struct=0,
            wht_bw_stl=0,
            wht_bw_cool=0,
            wht_mfbss_stl=0,
            wht_mfbss_cool=0,
            wht_mfbss_pbli=0,
            fwmass_stl=0,
            fwmass_cool=0,
            mass_cool_blanket=0,
            mass_liq_blanket=0,
            mass_stl_blanket=0,
            mass_segm_ib=0,
            mass_segm_ob=0,
            expected_blbmith=0.37000000000000011,
            expected_blbmoth=0.49699999999999994,
            expected_whtblkt=10654509.24412049,
            expected_fwmass=193353.16636179245,
            expected_fw_armour_vol=7.0163598878346534,
            expected_fw_armour_mass=135064.92784081708,
            expected_volfw=28.820872142117942,
            expected_armour_fw_bl_mass=10982927.3383231,
            expected_r_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ob=0.79000002145767212,
            expected_vfblkt=0.082598954955828252,
            expected_r_fci=0.050000000000000003,
            expected_r_backwall=0.02,
            expected_bz_r_ib=0.315,
            expected_bz_r_ob=0.41500000000000004,
            expected_f_vol_stff_plates=0.9100000262260437,
            expected_f_vol_stl_bz_struct=0.52999997138977051,
            expected_f_vol_stl_back_wall=0.86000001430511475,
            expected_f_vol_stl_fw=0.86000001430511475,
            expected_f_vol_mfbss_stl=0.51289999485015869,
            expected_f_vol_mfbss_he=0.04349999874830246,
            expected_f_vol_mfbss_pbli=0.44359999895095825,
            expected_vol_fci=77.328829099899536,
            expected_vol_bz_struct=245.67041910375485,
            expected_vol_bz_liq=342.92630631451561,
            expected_vol_bz_liq_ib=132.46921948004106,
            expected_vol_bz_liq_ob=210.45708683447458,
            expected_vol_bw=30.931531639959815,
            expected_vol_bss=701.04321499216383,
            expected_wht_cer=248225.54141067751,
            expected_wht_stl_struct=1015601.4577511634,
            expected_wht_cool_struct=651.10466637722732,
            expected_wht_bw_stl=207488.71769218749,
            expected_wht_bw_cool=24.419089893808916,
            expected_wht_mfbss_stl=2804607.478601912,
            expected_wht_mfbss_cool=171.96263515308456,
            expected_wht_mfbss_pbli=3033092.6337963375,
            expected_fwmass_stl=193330.41354515703,
            expected_fwmass_cool=22.752816635408795,
            expected_mass_cool_blanket=870.23920805952946,
            expected_mass_liq_blanket=6377738.5622731261,
            expected_mass_stl_blanket=4221028.0675904201,
            expected_mass_segm_ib=99402.417142669714,
            expected_mass_segm_ob=162542.70811995145,
        ),
        DcllMassesParam(
            fwareaib=723.16923304760132,
            dr_fw_inboard=0.018000000000000002,
            fwareaob=1168.1172772224481,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            blbuith=0.36499999999999999,
            blbmith=0.37000000000000011,
            dr_blkt_outboard=0.98199999999999998,
            blbuoth=0.46500000000000002,
            blbmoth=0.49699999999999994,
            a_plasma_surface=1403.2719775669307,
            a_plasma_surface_outboard=949.22962703393853,
            i_blkt_inboard=1,
            volblkt=1400.4860764869636,
            volblkti=402.02180553751157,
            volblkto=998.46427094945204,
            whtblkt=10654509.24412049,
            fwmass=193353.16636179245,
            fw_armour_vol=7.0163598878346534,
            fw_armour_thickness=0.0050000000000000001,
            fw_armour_mass=135064.92784081708,
            volfw=28.820872142117942,
            armour_fw_bl_mass=10982927.3383231,
            denstl=7800,
            den_liq=9753.2497999999996,
            ifci=1,
            den_ceramic=3210,
            th_wall_secondary=0.012500000000000001,
            nopol=2,
            r_f_liq_ib=0.79000002145767212,
            r_f_liq_ob=0.5,
            w_f_liq_ib=0.79000002145767212,
            w_f_liq_ob=0.79000002145767212,
            vfblkt=0.082598954955828252,
            icooldual=2,
            rhof_fw=5.6389735407435868,
            rhof_bl=5.6389735407435868,
            nblktmodti=32,
            nblktmodto=48,
            r_fci=0.050000000000000003,
            r_backwall=0.02,
            bz_r_ib=0.315,
            bz_r_ob=0.41500000000000004,
            f_vol_stff_plates=0.9100000262260437,
            f_vol_stl_bz_struct=0.52999997138977051,
            f_vol_stl_back_wall=0.86000001430511475,
            f_vol_stl_fw=0.86000001430511475,
            f_vol_mfbss_stl=0.51289999485015869,
            f_vol_mfbss_he=0.04349999874830246,
            f_vol_mfbss_pbli=0.44359999895095825,
            vol_fci=77.328829099899536,
            vol_bz_struct=245.67041910375485,
            vol_bz_liq=342.92630631451561,
            vol_bz_liq_ib=132.46921948004106,
            vol_bz_liq_ob=210.45708683447458,
            vol_bw=30.931531639959815,
            vol_bss=701.04321499216383,
            wht_cer=248225.54141067751,
            wht_stl_struct=1015601.4577511634,
            wht_cool_struct=651.10466637722732,
            wht_bw_stl=207488.71769218749,
            wht_bw_cool=24.419089893808916,
            wht_mfbss_stl=2804607.478601912,
            wht_mfbss_cool=171.96263515308456,
            wht_mfbss_pbli=3033092.6337963375,
            fwmass_stl=193330.41354515703,
            fwmass_cool=22.752816635408795,
            mass_cool_blanket=870.23920805952946,
            mass_liq_blanket=6377738.5622731261,
            mass_stl_blanket=4221028.0675904201,
            mass_segm_ib=99402.417142669714,
            mass_segm_ob=162542.70811995145,
            expected_blbmith=0.37000000000000011,
            expected_blbmoth=0.49699999999999994,
            expected_whtblkt=10673841.813263938,
            expected_fwmass=228388.37777659783,
            expected_fw_armour_vol=7.0163598878346534,
            expected_fw_armour_mass=135064.92784081708,
            expected_volfw=34.043157184860888,
            expected_armour_fw_bl_mass=11037295.118881352,
            expected_r_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ob=0.79000002145767212,
            expected_vfblkt=0.082624998748551323,
            expected_r_fci=0.050000000000000003,
            expected_r_backwall=0.02,
            expected_bz_r_ib=0.315,
            expected_bz_r_ob=0.41500000000000004,
            expected_f_vol_stff_plates=0.9100000262260437,
            expected_f_vol_stl_bz_struct=0.52999997138977051,
            expected_f_vol_stl_back_wall=0.86000001430511475,
            expected_f_vol_stl_fw=0.86000001430511475,
            expected_f_vol_mfbss_stl=0.51289999485015869,
            expected_f_vol_mfbss_he=0.04349999874830246,
            expected_f_vol_mfbss_pbli=0.44359999895095825,
            expected_vol_fci=77.462263633122888,
            expected_vol_bz_struct=246.20245377274514,
            expected_vol_bz_liq=343.48641311892817,
            expected_vol_bz_liq_ib=132.50745566270484,
            expected_vol_bz_liq_ob=210.97895745622333,
            expected_vol_bw=30.984905453249151,
            expected_vol_bss=702.35004050891826,
            expected_wht_cer=248653.86626232447,
            expected_wht_stl_struct=1017800.8889540406,
            expected_wht_cool_struct=652.51472729102318,
            expected_wht_bw_stl=207846.74923768779,
            expected_wht_bw_cool=24.461226182430458,
            expected_wht_mfbss_stl=2809835.5908482568,
            expected_wht_mfbss_cool=172.28319336510418,
            expected_wht_mfbss_pbli=3038746.6687598871,
            expected_fwmass_stl=228361.50219457873,
            expected_fwmass_cool=26.875582019101881,
            expected_mass_cool_blanket=876.13472885765964,
            expected_mass_liq_blanket=6388855.4588147905,
            expected_mass_stl_blanket=4263844.7312345635,
            expected_mass_segm_ib=99845.314502560635,
            expected_mass_segm_ob=163380.10530832107,
        ),
    ),
)
def test_dcll_masses(dcllmassesparam, monkeypatch, dcll):
    """
    Automatically generated Regression Unit Test for dcll_masses.

    This test was generated using data from ./dcll/dcll_mms_lt_IN.DAT.

    :param dcllmassesparam: the data used to mock and assert in this test.
    :type dcllmassesparam: dcllmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "fwareaib", dcllmassesparam.fwareaib)

    monkeypatch.setattr(build_variables, "dr_fw_inboard", dcllmassesparam.dr_fw_inboard)

    monkeypatch.setattr(build_variables, "fwareaob", dcllmassesparam.fwareaob)

    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", dcllmassesparam.dr_fw_outboard
    )

    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", dcllmassesparam.dr_blkt_inboard
    )

    monkeypatch.setattr(build_variables, "blbuith", dcllmassesparam.blbuith)

    monkeypatch.setattr(build_variables, "blbmith", dcllmassesparam.blbmith)

    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", dcllmassesparam.dr_blkt_outboard
    )

    monkeypatch.setattr(build_variables, "blbuoth", dcllmassesparam.blbuoth)

    monkeypatch.setattr(build_variables, "blbmoth", dcllmassesparam.blbmoth)

    monkeypatch.setattr(
        physics_variables, "a_plasma_surface", dcllmassesparam.a_plasma_surface
    )

    monkeypatch.setattr(
        physics_variables,
        "a_plasma_surface_outboard",
        dcllmassesparam.a_plasma_surface_outboard,
    )

    monkeypatch.setattr(
        fwbs_variables, "i_blkt_inboard", dcllmassesparam.i_blkt_inboard
    )

    monkeypatch.setattr(fwbs_variables, "volblkt", dcllmassesparam.volblkt)

    monkeypatch.setattr(fwbs_variables, "volblkti", dcllmassesparam.volblkti)

    monkeypatch.setattr(fwbs_variables, "volblkto", dcllmassesparam.volblkto)

    monkeypatch.setattr(fwbs_variables, "whtblkt", dcllmassesparam.whtblkt)

    monkeypatch.setattr(fwbs_variables, "fwmass", dcllmassesparam.fwmass)

    monkeypatch.setattr(fwbs_variables, "fw_armour_vol", dcllmassesparam.fw_armour_vol)

    monkeypatch.setattr(
        fwbs_variables, "fw_armour_thickness", dcllmassesparam.fw_armour_thickness
    )

    monkeypatch.setattr(
        fwbs_variables, "fw_armour_mass", dcllmassesparam.fw_armour_mass
    )

    monkeypatch.setattr(fwbs_variables, "volfw", dcllmassesparam.volfw)

    monkeypatch.setattr(
        fwbs_variables, "armour_fw_bl_mass", dcllmassesparam.armour_fw_bl_mass
    )

    monkeypatch.setattr(fwbs_variables, "denstl", dcllmassesparam.denstl)

    monkeypatch.setattr(fwbs_variables, "den_liq", dcllmassesparam.den_liq)

    monkeypatch.setattr(fwbs_variables, "ifci", dcllmassesparam.ifci)

    monkeypatch.setattr(fwbs_variables, "den_ceramic", dcllmassesparam.den_ceramic)

    monkeypatch.setattr(
        fwbs_variables, "th_wall_secondary", dcllmassesparam.th_wall_secondary
    )

    monkeypatch.setattr(fwbs_variables, "nopol", dcllmassesparam.nopol)

    monkeypatch.setattr(fwbs_variables, "r_f_liq_ib", dcllmassesparam.r_f_liq_ib)

    monkeypatch.setattr(fwbs_variables, "r_f_liq_ob", dcllmassesparam.r_f_liq_ob)

    monkeypatch.setattr(fwbs_variables, "w_f_liq_ib", dcllmassesparam.w_f_liq_ib)

    monkeypatch.setattr(fwbs_variables, "w_f_liq_ob", dcllmassesparam.w_f_liq_ob)

    monkeypatch.setattr(fwbs_variables, "vfblkt", dcllmassesparam.vfblkt)

    monkeypatch.setattr(fwbs_variables, "icooldual", dcllmassesparam.icooldual)

    monkeypatch.setattr(fwbs_variables, "rhof_fw", dcllmassesparam.rhof_fw)

    monkeypatch.setattr(fwbs_variables, "rhof_bl", dcllmassesparam.rhof_bl)

    monkeypatch.setattr(fwbs_variables, "nblktmodti", dcllmassesparam.nblktmodti)

    monkeypatch.setattr(fwbs_variables, "nblktmodto", dcllmassesparam.nblktmodto)

    monkeypatch.setattr(dcll_module, "r_fci", dcllmassesparam.r_fci)

    monkeypatch.setattr(dcll_module, "r_backwall", dcllmassesparam.r_backwall)

    monkeypatch.setattr(dcll_module, "bz_r_ib", dcllmassesparam.bz_r_ib)

    monkeypatch.setattr(dcll_module, "bz_r_ob", dcllmassesparam.bz_r_ob)

    monkeypatch.setattr(
        dcll_module, "f_vol_stff_plates", dcllmassesparam.f_vol_stff_plates
    )

    monkeypatch.setattr(
        dcll_module, "f_vol_stl_bz_struct", dcllmassesparam.f_vol_stl_bz_struct
    )

    monkeypatch.setattr(
        dcll_module, "f_vol_stl_back_wall", dcllmassesparam.f_vol_stl_back_wall
    )

    monkeypatch.setattr(dcll_module, "f_vol_stl_fw", dcllmassesparam.f_vol_stl_fw)

    monkeypatch.setattr(dcll_module, "f_vol_mfbss_stl", dcllmassesparam.f_vol_mfbss_stl)

    monkeypatch.setattr(dcll_module, "f_vol_mfbss_he", dcllmassesparam.f_vol_mfbss_he)

    monkeypatch.setattr(
        dcll_module, "f_vol_mfbss_pbli", dcllmassesparam.f_vol_mfbss_pbli
    )

    monkeypatch.setattr(dcll_module, "vol_fci", dcllmassesparam.vol_fci)

    monkeypatch.setattr(dcll_module, "vol_bz_struct", dcllmassesparam.vol_bz_struct)

    monkeypatch.setattr(dcll_module, "vol_bz_liq", dcllmassesparam.vol_bz_liq)

    monkeypatch.setattr(dcll_module, "vol_bz_liq_ib", dcllmassesparam.vol_bz_liq_ib)

    monkeypatch.setattr(dcll_module, "vol_bz_liq_ob", dcllmassesparam.vol_bz_liq_ob)

    monkeypatch.setattr(dcll_module, "vol_bw", dcllmassesparam.vol_bw)

    monkeypatch.setattr(dcll_module, "vol_bss", dcllmassesparam.vol_bss)

    monkeypatch.setattr(dcll_module, "wht_cer", dcllmassesparam.wht_cer)

    monkeypatch.setattr(dcll_module, "wht_stl_struct", dcllmassesparam.wht_stl_struct)

    monkeypatch.setattr(dcll_module, "wht_cool_struct", dcllmassesparam.wht_cool_struct)

    monkeypatch.setattr(dcll_module, "wht_bw_stl", dcllmassesparam.wht_bw_stl)

    monkeypatch.setattr(dcll_module, "wht_bw_cool", dcllmassesparam.wht_bw_cool)

    monkeypatch.setattr(dcll_module, "wht_mfbss_stl", dcllmassesparam.wht_mfbss_stl)

    monkeypatch.setattr(dcll_module, "wht_mfbss_cool", dcllmassesparam.wht_mfbss_cool)

    monkeypatch.setattr(dcll_module, "wht_mfbss_pbli", dcllmassesparam.wht_mfbss_pbli)

    monkeypatch.setattr(dcll_module, "fwmass_stl", dcllmassesparam.fwmass_stl)

    monkeypatch.setattr(dcll_module, "fwmass_cool", dcllmassesparam.fwmass_cool)

    monkeypatch.setattr(
        dcll_module, "mass_cool_blanket", dcllmassesparam.mass_cool_blanket
    )

    monkeypatch.setattr(
        dcll_module, "mass_liq_blanket", dcllmassesparam.mass_liq_blanket
    )

    monkeypatch.setattr(
        dcll_module, "mass_stl_blanket", dcllmassesparam.mass_stl_blanket
    )

    monkeypatch.setattr(dcll_module, "mass_segm_ib", dcllmassesparam.mass_segm_ib)

    monkeypatch.setattr(dcll_module, "mass_segm_ob", dcllmassesparam.mass_segm_ob)

    dcll.dcll_masses(False)

    assert build_variables.blbmith == pytest.approx(dcllmassesparam.expected_blbmith)

    assert build_variables.blbmoth == pytest.approx(dcllmassesparam.expected_blbmoth)

    assert fwbs_variables.whtblkt == pytest.approx(dcllmassesparam.expected_whtblkt)

    assert fwbs_variables.fwmass == pytest.approx(dcllmassesparam.expected_fwmass)

    assert fwbs_variables.fw_armour_vol == pytest.approx(
        dcllmassesparam.expected_fw_armour_vol
    )

    assert fwbs_variables.fw_armour_mass == pytest.approx(
        dcllmassesparam.expected_fw_armour_mass
    )

    assert fwbs_variables.volfw == pytest.approx(dcllmassesparam.expected_volfw)

    assert fwbs_variables.armour_fw_bl_mass == pytest.approx(
        dcllmassesparam.expected_armour_fw_bl_mass
    )

    assert fwbs_variables.r_f_liq_ib == pytest.approx(
        dcllmassesparam.expected_r_f_liq_ib
    )

    assert fwbs_variables.w_f_liq_ib == pytest.approx(
        dcllmassesparam.expected_w_f_liq_ib
    )

    assert fwbs_variables.w_f_liq_ob == pytest.approx(
        dcllmassesparam.expected_w_f_liq_ob
    )

    assert fwbs_variables.vfblkt == pytest.approx(dcllmassesparam.expected_vfblkt)

    assert dcll_module.r_fci == pytest.approx(dcllmassesparam.expected_r_fci)

    assert dcll_module.r_backwall == pytest.approx(dcllmassesparam.expected_r_backwall)

    assert dcll_module.bz_r_ib == pytest.approx(dcllmassesparam.expected_bz_r_ib)

    assert dcll_module.bz_r_ob == pytest.approx(dcllmassesparam.expected_bz_r_ob)

    assert dcll_module.f_vol_stff_plates == pytest.approx(
        dcllmassesparam.expected_f_vol_stff_plates
    )

    assert dcll_module.f_vol_stl_bz_struct == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_bz_struct
    )

    assert dcll_module.f_vol_stl_back_wall == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_back_wall
    )

    assert dcll_module.f_vol_stl_fw == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_fw
    )

    assert dcll_module.f_vol_mfbss_stl == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_stl
    )

    assert dcll_module.f_vol_mfbss_he == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_he
    )

    assert dcll_module.f_vol_mfbss_pbli == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_pbli
    )

    assert dcll_module.vol_fci == pytest.approx(dcllmassesparam.expected_vol_fci)

    assert dcll_module.vol_bz_struct == pytest.approx(
        dcllmassesparam.expected_vol_bz_struct
    )

    assert dcll_module.vol_bz_liq == pytest.approx(dcllmassesparam.expected_vol_bz_liq)

    assert dcll_module.vol_bz_liq_ib == pytest.approx(
        dcllmassesparam.expected_vol_bz_liq_ib
    )

    assert dcll_module.vol_bz_liq_ob == pytest.approx(
        dcllmassesparam.expected_vol_bz_liq_ob
    )

    assert dcll_module.vol_bw == pytest.approx(dcllmassesparam.expected_vol_bw)

    assert dcll_module.vol_bss == pytest.approx(dcllmassesparam.expected_vol_bss)

    assert dcll_module.wht_cer == pytest.approx(dcllmassesparam.expected_wht_cer)

    assert dcll_module.wht_stl_struct == pytest.approx(
        dcllmassesparam.expected_wht_stl_struct
    )

    assert dcll_module.wht_cool_struct == pytest.approx(
        dcllmassesparam.expected_wht_cool_struct
    )

    assert dcll_module.wht_bw_stl == pytest.approx(dcllmassesparam.expected_wht_bw_stl)

    assert dcll_module.wht_bw_cool == pytest.approx(
        dcllmassesparam.expected_wht_bw_cool
    )

    assert dcll_module.wht_mfbss_stl == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_stl
    )

    assert dcll_module.wht_mfbss_cool == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_cool
    )

    assert dcll_module.wht_mfbss_pbli == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_pbli
    )

    assert dcll_module.fwmass_stl == pytest.approx(dcllmassesparam.expected_fwmass_stl)

    assert dcll_module.fwmass_cool == pytest.approx(
        dcllmassesparam.expected_fwmass_cool
    )

    assert dcll_module.mass_cool_blanket == pytest.approx(
        dcllmassesparam.expected_mass_cool_blanket
    )

    assert dcll_module.mass_liq_blanket == pytest.approx(
        dcllmassesparam.expected_mass_liq_blanket
    )

    assert dcll_module.mass_stl_blanket == pytest.approx(
        dcllmassesparam.expected_mass_stl_blanket
    )

    assert dcll_module.mass_segm_ib == pytest.approx(
        dcllmassesparam.expected_mass_segm_ib
    )

    assert dcll_module.mass_segm_ob == pytest.approx(
        dcllmassesparam.expected_mass_segm_ob
    )
