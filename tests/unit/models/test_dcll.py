from typing import Any, NamedTuple

import pytest

from process.data_structure.build_variables import InboardBlanketConfiguration


@pytest.fixture
def dcll(process_models):
    """Fixture to get the DCLL instance from process_models.

    :returns: initialised DCLL object
    :rtype: process.dcll.DCLL
    """
    return process_models.dcll


class DcllNeutronicsAndPowerParam(NamedTuple):
    a_fw_outboard: Any = None

    a_fw_total: Any = None

    p_beam_orbit_loss_mw: Any = None

    f_ster_div_single: Any = None

    p_div_rad_total_mw: Any = None

    f_a_fw_outboard_hcd: Any = None

    p_fw_hcd_rad_total_mw: Any = None

    p_fw_hcd_nuclear_heat_mw: Any = None

    p_shld_nuclear_heat_mw: Any = None

    p_fw_rad_total_mw: Any = None

    p_fw_nuclear_heat_total_mw: Any = None

    p_fw_inboard_surface_heat_mw: Any = None

    p_fw_outboard_surface_heat_mw: Any = None

    p_blkt_nuclear_heat_total_mw: Any = None

    pnuc_fw_ratio_dcll: Any = None

    f_p_blkt_multiplication: Any = None

    n_divertors: Any = None

    p_neutron_total_mw: Any = None

    p_plasma_rad_mw: Any = None

    p_fw_alpha_surface_total_mw: Any = None

    expected_p_fw_nuclear_heat_total_mw: Any = None

    expected_p_blkt_nuclear_heat_total_mw: Any = None

    expected_p_blkt_multiplication_mw: Any = None


@pytest.mark.parametrize(
    "dcllneutronicsandpowerparam",
    [
        DcllNeutronicsAndPowerParam(
            a_fw_outboard=988.92586580655245,
            a_fw_total=1601.1595634509963,
            p_beam_orbit_loss_mw=0,
            f_ster_div_single=0.115,
            p_div_rad_total_mw=33.056596978820579,
            f_a_fw_outboard_hcd=0,
            p_fw_hcd_rad_total_mw=0,
            p_fw_hcd_nuclear_heat_mw=0,
            p_shld_nuclear_heat_mw=0,
            p_fw_rad_total_mw=0,
            p_fw_nuclear_heat_total_mw=0,
            p_fw_inboard_surface_heat_mw=0,
            p_fw_outboard_surface_heat_mw=0,
            p_blkt_nuclear_heat_total_mw=0,
            pnuc_fw_ratio_dcll=0.14000000000000001,
            f_p_blkt_multiplication=1.2689999999999999,
            n_divertors=1,
            p_neutron_total_mw=1587.7386535917431,
            p_plasma_rad_mw=287.44866938104849,
            p_fw_alpha_surface_total_mw=19.835845058655043,
            expected_p_fw_nuclear_heat_total_mw=196.72081918001697,
            expected_p_blkt_nuclear_heat_total_mw=1533.4949914565693,
            expected_p_blkt_multiplication_mw=325.06710220789364,
        ),
        DcllNeutronicsAndPowerParam(
            a_fw_outboard=1168.1172772224481,
            a_fw_total=1891.2865102700493,
            p_beam_orbit_loss_mw=0,
            f_ster_div_single=0.115,
            p_div_rad_total_mw=33.056596978820579,
            f_a_fw_outboard_hcd=0,
            p_fw_hcd_rad_total_mw=0,
            p_fw_hcd_nuclear_heat_mw=0,
            p_shld_nuclear_heat_mw=0,
            p_fw_rad_total_mw=254.39207240222791,
            p_fw_nuclear_heat_total_mw=196.72081918001697,
            p_fw_inboard_surface_heat_mw=97.271629070225231,
            p_fw_outboard_surface_heat_mw=176.95628839065773,
            p_blkt_nuclear_heat_total_mw=1533.4949914565693,
            pnuc_fw_ratio_dcll=0.14000000000000001,
            f_p_blkt_multiplication=1.2689999999999999,
            n_divertors=1,
            p_neutron_total_mw=1587.2430556964196,
            p_plasma_rad_mw=287.44866938104849,
            p_fw_alpha_surface_total_mw=19.829653483586444,
            expected_p_fw_nuclear_heat_total_mw=196.65941460078642,
            expected_p_blkt_nuclear_heat_total_mw=1533.0163252173013,
            expected_p_blkt_multiplication_mw=324.96563552675644,
        ),
    ],
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
        dcll.data.first_wall, "a_fw_outboard", dcllneutronicsandpowerparam.a_fw_outboard
    )

    monkeypatch.setattr(
        dcll.data.first_wall, "a_fw_total", dcllneutronicsandpowerparam.a_fw_total
    )

    monkeypatch.setattr(
        dcll.data.current_drive,
        "p_beam_orbit_loss_mw",
        dcllneutronicsandpowerparam.p_beam_orbit_loss_mw,
    )

    for field in [
        "f_ster_div_single",
        "p_div_rad_total_mw",
        "f_a_fw_outboard_hcd",
        dcllneutronicsandpowerparam.f_a_fw_outboard_hcd,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_hcd_rad_total_mw",
        dcllneutronicsandpowerparam.p_fw_hcd_rad_total_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_hcd_nuclear_heat_mw",
        dcllneutronicsandpowerparam.p_fw_hcd_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_shld_nuclear_heat_mw",
        dcllneutronicsandpowerparam.p_shld_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_rad_total_mw",
        dcllneutronicsandpowerparam.p_fw_rad_total_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_nuclear_heat_total_mw",
        dcllneutronicsandpowerparam.p_fw_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_inboard_surface_heat_mw",
        dcllneutronicsandpowerparam.p_fw_inboard_surface_heat_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_fw_outboard_surface_heat_mw",
        dcllneutronicsandpowerparam.p_fw_outboard_surface_heat_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "p_blkt_nuclear_heat_total_mw",
        dcllneutronicsandpowerparam.p_blkt_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        dcll.data.fwbs,
        "pnuc_fw_ratio_dcll",
        "f_p_blkt_multiplication",
    ]:
        monkeypatch.setattr(
            dcll.data.fwbs,
            field,
            getattr(dcllneutronicsandpowerparam, field),
        )

    monkeypatch.setattr(
        dcll.data.divertor, "n_divertors", dcllneutronicsandpowerparam.n_divertors
    )

    for field in [
        "p_neutron_total_mw",
        dcllneutronicsandpowerparam.p_neutron_total_mw,
    )

    monkeypatch.setattr(
        dcll.data.physics,
        "p_plasma_rad_mw",
        dcllneutronicsandpowerparam.p_plasma_rad_mw,
    )

    monkeypatch.setattr(
        dcll.data.physics,
        "p_fw_alpha_surface_total_mw",
        dcllneutronicsandpowerparam.p_fw_alpha_surface_total_mw,
    )

    dcll.dcll_neutronics_and_power(False)

    assert dcll.data.fwbs.p_fw_nuclear_heat_total_mw == pytest.approx(
        dcllneutronicsandpowerparam.expected_p_fw_nuclear_heat_total_mw
    )

    assert dcll.data.fwbs.p_blkt_nuclear_heat_total_mw == pytest.approx(
        dcllneutronicsandpowerparam.expected_p_blkt_nuclear_heat_total_mw
    )

    assert dcll.data.fwbs.p_blkt_multiplication_mw == pytest.approx(
        dcllneutronicsandpowerparam.expected_p_blkt_multiplication_mw
    )


class DcllMassesParam(NamedTuple):
    a_fw_inboard: Any = None

    dr_fw_inboard: Any = None

    a_fw_outboard: Any = None

    dr_fw_outboard: Any = None

    dr_blkt_inboard: Any = None

    blbuith: Any = None

    dr_blkt_outboard: Any = None

    blbuoth: Any = None

    a_plasma_surface: Any = None

    a_plasma_surface_outboard: Any = None

    i_blkt_inboard: Any = None

    vol_blkt_total: Any = None

    vol_blkt_inboard: Any = None

    vol_blkt_outboard: Any = None

    fw_armour_thickness: Any = None

    den_steel: Any = None

    den_liq: Any = None

    i_blkt_liquid_breeder_channel_type: Any = None

    den_ceramic: Any = None

    th_wall_secondary: Any = None

    nopol: Any = None

    r_f_liq_ob: Any = None

    i_blkt_dual_coolant: Any = None

    den_fw_coolant: Any = None

    den_blkt_coolant: Any = None

    n_blkt_inboard_modules_toroidal: Any = None

    n_blkt_outboard_modules_toroidal: Any = None

    expected_blbmith: Any = None

    expected_blbmoth: Any = None

    expected_m_blkt_total: Any = None

    expected_m_fw_total: Any = None

    expected_fw_armour_vol: Any = None

    expected_fw_armour_mass: Any = None

    expected_vol_fw_total: Any = None

    expected_armour_fw_bl_mass: Any = None

    expected_r_f_liq_ib: Any = None

    expected_w_f_liq_ib: Any = None

    expected_w_f_liq_ob: Any = None

    expected_f_a_blkt_cooling_channels: Any = None

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
    [
        DcllMassesParam(
            a_fw_inboard=612.23369764444396,
            dr_fw_inboard=0.018000000000000002,
            a_fw_outboard=988.92586580655245,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            blbuith=0.36499999999999999,
            dr_blkt_outboard=0.98199999999999998,
            blbuoth=0.46500000000000002,
            a_plasma_surface=1403.2719775669307,
            a_plasma_surface_outboard=949.22962703393853,
            i_blkt_inboard=InboardBlanketConfiguration.INBOARD_BLANKET_PRESENT,
            vol_blkt_total=1397.9003011502937,
            vol_blkt_inboard=401.90579863726225,
            vol_blkt_outboard=995.99450251303142,
            fw_armour_thickness=0.0050000000000000001,
            den_steel=7800,
            den_liq=9753.2497999999996,
            i_blkt_liquid_breeder_channel_type=1,
            den_ceramic=3210,
            th_wall_secondary=0.012500000000000001,
            nopol=2,
            r_f_liq_ob=0.5,
            i_blkt_dual_coolant=2,
            den_fw_coolant=5.6389735407435868,
            den_blkt_coolant=5.6389735407435868,
            n_blkt_inboard_modules_toroidal=32,
            n_blkt_outboard_modules_toroidal=48,
            expected_blbmith=0.37000000000000011,
            expected_blbmoth=0.49699999999999994,
            expected_m_blkt_total=10654509.24412049,
            expected_m_fw_total=193353.16636179245,
            expected_fw_armour_vol=7.0163598878346534,
            expected_fw_armour_mass=135064.92784081708,
            expected_vol_fw_total=28.820872142117942,
            expected_armour_fw_bl_mass=10982927.3383231,
            expected_r_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ob=0.79000002145767212,
            expected_f_a_blkt_cooling_channels=0.082598954955828252,
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
            a_fw_inboard=723.16923304760132,
            dr_fw_inboard=0.018000000000000002,
            a_fw_outboard=1168.1172772224481,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            blbuith=0.36499999999999999,
            dr_blkt_outboard=0.98199999999999998,
            blbuoth=0.46500000000000002,
            a_plasma_surface=1403.2719775669307,
            a_plasma_surface_outboard=949.22962703393853,
            i_blkt_inboard=InboardBlanketConfiguration.INBOARD_BLANKET_PRESENT,
            vol_blkt_total=1400.4860764869636,
            vol_blkt_inboard=402.02180553751157,
            vol_blkt_outboard=998.46427094945204,
            fw_armour_thickness=0.0050000000000000001,
            den_steel=7800,
            den_liq=9753.2497999999996,
            i_blkt_liquid_breeder_channel_type=1,
            den_ceramic=3210,
            th_wall_secondary=0.012500000000000001,
            nopol=2,
            r_f_liq_ob=0.5,
            i_blkt_dual_coolant=2,
            den_fw_coolant=5.6389735407435868,
            den_blkt_coolant=5.6389735407435868,
            n_blkt_inboard_modules_toroidal=32,
            n_blkt_outboard_modules_toroidal=48,
            expected_blbmith=0.37000000000000011,
            expected_blbmoth=0.49699999999999994,
            expected_m_blkt_total=10673841.813263938,
            expected_m_fw_total=228388.37777659783,
            expected_fw_armour_vol=7.0163598878346534,
            expected_fw_armour_mass=135064.92784081708,
            expected_vol_fw_total=34.043157184860888,
            expected_armour_fw_bl_mass=11037295.118881352,
            expected_r_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ib=0.79000002145767212,
            expected_w_f_liq_ob=0.79000002145767212,
            expected_f_a_blkt_cooling_channels=0.082624998748551323,
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
    ],
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
    for field in [
        "a_fw_inboard",
        "a_fw_outboard",
    ]:
        monkeypatch.setattr(dcll.data.first_wall, field, getattr(dcllmassesparam, field))

    for field in [
        "dr_fw_inboard",
        "dr_fw_outboard",
        "dr_blkt_inboard",
        "blbuith",
        "dr_blkt_outboard",
        "blbuoth",
        "i_blkt_inboard",
    ]:
        monkeypatch.setattr(dcll.data.build, field, getattr(dcllmassesparam, field))

    for field in [
        "a_plasma_surface",
        "a_plasma_surface_outboard",
    ]:
        monkeypatch.setattr(dcll.data.physics, field, getattr(dcllmassesparam, field))

    for field in [
        "vol_blkt_total",
        "vol_blkt_inboard",
        "vol_blkt_outboard",
        "fw_armour_thickness",
        "den_steel",
        "den_liq",
        "i_blkt_liquid_breeder_channel_type",
        "den_ceramic",
        "th_wall_secondary",
        "nopol",
        "r_f_liq_ob",
        "i_blkt_dual_coolant",
        "den_fw_coolant",
        "den_blkt_coolant",
        "n_blkt_inboard_modules_toroidal",
        "n_blkt_outboard_modules_toroidal",
    ]:
        monkeypatch.setattr(dcll.data.fwbs, field, getattr(dcllmassesparam, field))

    dcll.dcll_masses(output=False)

    for field in ["blbmith", "blbmoth"]:
        assert getattr(dcll.data.build, field) == pytest.approx(
            getattr(dcllmassesparam, f"expected_{field}")
        )

    assert dcll.data.fwbs.m_blkt_total == pytest.approx(
        dcllmassesparam.expected_m_blkt_total
    )

    assert dcll.data.fwbs.m_fw_total == pytest.approx(
        dcllmassesparam.expected_m_fw_total
    )

    assert dcll.data.fwbs.fw_armour_vol == pytest.approx(
        dcllmassesparam.expected_fw_armour_vol
    )

    assert dcll.data.fwbs.fw_armour_mass == pytest.approx(
        dcllmassesparam.expected_fw_armour_mass
    )

    assert dcll.data.fwbs.vol_fw_total == pytest.approx(
        dcllmassesparam.expected_vol_fw_total
    )

    assert dcll.data.fwbs.armour_fw_bl_mass == pytest.approx(
        dcllmassesparam.expected_armour_fw_bl_mass
    )

    assert dcll.data.fwbs.r_f_liq_ib == pytest.approx(
        dcllmassesparam.expected_r_f_liq_ib
    )

    assert dcll.data.fwbs.w_f_liq_ib == pytest.approx(
        dcllmassesparam.expected_w_f_liq_ib
    )

    assert dcll.data.fwbs.w_f_liq_ob == pytest.approx(
        dcllmassesparam.expected_w_f_liq_ob
    )

    assert dcll.data.fwbs.f_a_blkt_cooling_channels == pytest.approx(
        dcllmassesparam.expected_f_a_blkt_cooling_channels
    )

    assert dcll.data.dcll.r_fci == pytest.approx(dcllmassesparam.expected_r_fci)

    assert dcll.data.dcll.r_backwall == pytest.approx(
        dcllmassesparam.expected_r_backwall
    )

    assert dcll.data.dcll.bz_r_ib == pytest.approx(dcllmassesparam.expected_bz_r_ib)

    assert dcll.data.dcll.bz_r_ob == pytest.approx(dcllmassesparam.expected_bz_r_ob)

    assert dcll.data.dcll.f_vol_stff_plates == pytest.approx(
        dcllmassesparam.expected_f_vol_stff_plates
    )

    assert dcll.data.dcll.f_vol_stl_bz_struct == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_bz_struct
    )

    assert dcll.data.dcll.f_vol_stl_back_wall == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_back_wall
    )

    assert dcll.data.dcll.f_vol_stl_fw == pytest.approx(
        dcllmassesparam.expected_f_vol_stl_fw
    )

    assert dcll.data.dcll.f_vol_mfbss_stl == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_stl
    )

    assert dcll.data.dcll.f_vol_mfbss_he == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_he
    )

    assert dcll.data.dcll.f_vol_mfbss_pbli == pytest.approx(
        dcllmassesparam.expected_f_vol_mfbss_pbli
    )

    assert dcll.data.dcll.vol_fci == pytest.approx(dcllmassesparam.expected_vol_fci)

    assert dcll.data.dcll.vol_bz_struct == pytest.approx(
        dcllmassesparam.expected_vol_bz_struct
    )

    assert dcll.data.dcll.vol_bz_liq == pytest.approx(
        dcllmassesparam.expected_vol_bz_liq
    )

    assert dcll.data.dcll.vol_bz_liq_ib == pytest.approx(
        dcllmassesparam.expected_vol_bz_liq_ib
    )

    assert dcll.data.dcll.vol_bz_liq_ob == pytest.approx(
        dcllmassesparam.expected_vol_bz_liq_ob
    )

    assert dcll.data.dcll.vol_bw == pytest.approx(dcllmassesparam.expected_vol_bw)

    assert dcll.data.dcll.vol_bss == pytest.approx(dcllmassesparam.expected_vol_bss)

    assert dcll.data.dcll.wht_cer == pytest.approx(dcllmassesparam.expected_wht_cer)

    assert dcll.data.dcll.wht_stl_struct == pytest.approx(
        dcllmassesparam.expected_wht_stl_struct
    )

    assert dcll.data.dcll.wht_cool_struct == pytest.approx(
        dcllmassesparam.expected_wht_cool_struct
    )

    assert dcll.data.dcll.wht_bw_stl == pytest.approx(
        dcllmassesparam.expected_wht_bw_stl
    )

    assert dcll.data.dcll.wht_bw_cool == pytest.approx(
        dcllmassesparam.expected_wht_bw_cool
    )

    assert dcll.data.dcll.wht_mfbss_stl == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_stl
    )

    assert dcll.data.dcll.wht_mfbss_cool == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_cool
    )

    assert dcll.data.dcll.wht_mfbss_pbli == pytest.approx(
        dcllmassesparam.expected_wht_mfbss_pbli
    )

    assert dcll.data.dcll.fwmass_stl == pytest.approx(
        dcllmassesparam.expected_fwmass_stl
    )

    assert dcll.data.dcll.fwmass_cool == pytest.approx(
        dcllmassesparam.expected_fwmass_cool
    )

    assert dcll.data.dcll.mass_cool_blanket == pytest.approx(
        dcllmassesparam.expected_mass_cool_blanket
    )

    assert dcll.data.dcll.mass_liq_blanket == pytest.approx(
        dcllmassesparam.expected_mass_liq_blanket
    )

    assert dcll.data.dcll.mass_stl_blanket == pytest.approx(
        dcllmassesparam.expected_mass_stl_blanket
    )

    assert dcll.data.dcll.mass_segm_ib == pytest.approx(
        dcllmassesparam.expected_mass_segm_ib
    )

    assert dcll.data.dcll.mass_segm_ob == pytest.approx(
        dcllmassesparam.expected_mass_segm_ob
    )
