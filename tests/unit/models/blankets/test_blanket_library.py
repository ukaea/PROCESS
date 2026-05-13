from typing import Any, NamedTuple

import numpy as np
import pytest

from process.data_structure import (
    divertor_variables,
    physics_variables,
)


@pytest.fixture
def blanket_library_fixture(process_models):
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return process_models.blanket_library


class PrimaryCoolantPropertiesParam(NamedTuple):
    i_fw_coolant_type: Any = None

    temp_fw_coolant_in: Any = None

    temp_fw_coolant_out: Any = None

    pres_fw_coolant: Any = None

    den_fw_coolant: Any = None

    cp_fw: Any = None

    cv_fw: Any = None

    i_blkt_coolant_type: Any = None

    temp_blkt_coolant_in: Any = None

    temp_blkt_coolant_out: Any = None

    pres_blkt_coolant: Any = None

    den_blkt_coolant: Any = None

    i_blkt_dual_coolant: Any = None

    visc_blkt_coolant: Any = None

    cp_bl: Any = None

    cv_bl: Any = None

    visc_fw_coolant: Any = None

    i_fw_blkt_shared_coolant: Any = None

    expected_den_fw_coolant: Any = None

    expected_cp_fw: Any = None

    expected_cv_fw: Any = None

    expected_den_blkt_coolant: Any = None

    expected_visc_blkt_coolant: Any = None

    expected_cp_bl: Any = None

    expected_cv_bl: Any = None

    expected_visc_fw_coolant: Any = None


@pytest.mark.parametrize(
    "primarycoolantpropertiesparam",
    [
        PrimaryCoolantPropertiesParam(
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            temp_fw_coolant_out=773,
            pres_fw_coolant=8000000,
            den_fw_coolant=0,
            cp_fw=0,
            cv_fw=0,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_in=573,
            temp_blkt_coolant_out=773,
            pres_blkt_coolant=8000000,
            den_blkt_coolant=0,
            i_blkt_dual_coolant=2,
            visc_blkt_coolant=0,
            cp_bl=0,
            cv_bl=0,
            visc_fw_coolant=0,
            i_fw_blkt_shared_coolant=0,
            expected_den_fw_coolant=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_den_blkt_coolant=5.6389735407435868,
            expected_visc_blkt_coolant=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw_coolant=3.5036293160410249e-05,
        ),
        PrimaryCoolantPropertiesParam(
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            temp_fw_coolant_out=773,
            pres_fw_coolant=8000000,
            den_fw_coolant=5.6389735407435868,
            cp_fw=5188.5588430173211,
            cv_fw=3123.5687263525392,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_in=573,
            temp_blkt_coolant_out=773,
            pres_blkt_coolant=8000000,
            den_blkt_coolant=5.6389735407435868,
            i_blkt_dual_coolant=2,
            visc_blkt_coolant=3.5036293160410249e-05,
            cp_bl=5188.5588430173211,
            cv_bl=3123.5687263525392,
            visc_fw_coolant=3.5036293160410249e-05,
            i_fw_blkt_shared_coolant=0,
            expected_den_fw_coolant=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_den_blkt_coolant=5.6389735407435868,
            expected_visc_blkt_coolant=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw_coolant=3.5036293160410249e-05,
        ),
    ],
)
def test_primary_coolant_properties(
    primarycoolantpropertiesparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for primary_coolant_properties.

    This test was generated using data from dcll/dcll_mms_lt_IN.DAT.

    :param primarycoolantpropertiesparam: the data used to mock and assert in this test.
    :type primarycoolantpropertiesparam: primarycoolantpropertiesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    # monkeypatch doesnt work for strings
    # but helium is the default
    # monkeypatch.setattr(
    #     blanket_library_fixture.data.fwbs, "i_fw_coolant_type", primarycoolantpropertiesparam.i_fw_coolant_type
    # )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "temp_fw_coolant_in",
        primarycoolantpropertiesparam.temp_fw_coolant_in,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "temp_fw_coolant_out",
        primarycoolantpropertiesparam.temp_fw_coolant_out,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "pres_fw_coolant",
        primarycoolantpropertiesparam.pres_fw_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "den_fw_coolant",
        primarycoolantpropertiesparam.den_fw_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "cp_fw", primarycoolantpropertiesparam.cp_fw
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "cv_fw", primarycoolantpropertiesparam.cv_fw
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_coolant_type",
        primarycoolantpropertiesparam.i_blkt_coolant_type,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "temp_blkt_coolant_in",
        primarycoolantpropertiesparam.temp_blkt_coolant_in,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "temp_blkt_coolant_out",
        primarycoolantpropertiesparam.temp_blkt_coolant_out,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "pres_blkt_coolant",
        primarycoolantpropertiesparam.pres_blkt_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "den_blkt_coolant",
        primarycoolantpropertiesparam.den_blkt_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_dual_coolant",
        primarycoolantpropertiesparam.i_blkt_dual_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "visc_blkt_coolant",
        primarycoolantpropertiesparam.visc_blkt_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "cp_bl", primarycoolantpropertiesparam.cp_bl
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "cv_bl", primarycoolantpropertiesparam.cv_bl
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "visc_fw_coolant",
        primarycoolantpropertiesparam.visc_fw_coolant,
    )

    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_fw_blkt_shared_coolant",
        primarycoolantpropertiesparam.i_fw_blkt_shared_coolant,
    )

    blanket_library_fixture.primary_coolant_properties(output=False)

    assert blanket_library_fixture.data.fwbs.den_fw_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_den_fw_coolant, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.cp_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_fw, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.cv_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_fw, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.den_blkt_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_den_blkt_coolant, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.visc_blkt_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_blkt_coolant, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.cp_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_bl, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.cv_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_bl, rel=1e-4
    )

    assert blanket_library_fixture.data.fwbs.visc_fw_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_fw_coolant, rel=1e-4
    )


def test_deltap_tot_inboard_first_wall(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "radius_fw_channel", 0.006)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 0.22481)

    data = {
        "icoolpump": 1,
        "vel_coolant": 15.9,
        "len_pipe": 4,
        "n_pipe_90_deg_bends": 2,
        "n_pipe_180_deg_bends": 0,
        "den_coolant": 5.6,
        "visc_coolant_dynamic": 3.5e-5,
        "coolant_electrical_conductivity": 0.0,
        "pol_channel_length": 1.89,
        "nopolchan": 0,
        "label": "Inboard first wall",
    }

    assert (
        pytest.approx(blanket_library_fixture.total_pressure_drop(False, **data))
        == 5884.982168510442
    )


def test_deltap_tot_outboard_blanket_breeder_liquid(
    monkeypatch, blanket_library_fixture
):
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "radius_fw_channel", 0.006)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 0.22481)
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_channel_type", 1
    )
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "b_bz_liq", 0.11625)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "b_mag_blkt", [8.393, 3.868])
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "bz_channel_conduct_liq", 833000
    )
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "th_wall_secondary", 0.0125)

    data = {
        "icoolpump": 2,
        "vel_coolant": 0.06,
        "len_pipe": 4.7,
        "n_pipe_90_deg_bends": 2,
        "n_pipe_180_deg_bends": 1,
        "den_coolant": 9753.2,
        "visc_coolant_dynamic": 0.0017,
        "coolant_electrical_conductivity": 861800.8,
        "pol_channel_length": 1.89,
        "nopolchan": 0,
        "label": "Outboard blanket breeder liquid",
    }

    assert (
        pytest.approx(blanket_library_fixture.total_pressure_drop(False, **data))
        == 56.95922064419226
    )


def test_pumppower_primary_helium(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "etaiso", 0.9)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "etaiso_liq", 0.85)

    data = {
        "i_liquid_breeder": 2,
        "temp_coolant_pump_outlet": 570,
        "temp_coolant_pump_inlet": 720,
        "pres_coolant_pump_inlet": 1700000,
        "dpres_coolant": 303517.3,
        "mflow_coolant_total": 35677.7,
        "primary_coolant_switch": 1,
        "den_coolant": 9753.25,
        "label": "Liquid Metal Breeder/Coolant",
    }

    assert (
        pytest.approx(blanket_library_fixture.coolant_pumping_power(False, **data))
        == 1.8251284651310427
    )


def test_pumppower_secondary_pb_li(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "etaiso", 0.9)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "etaiso_liq", 0.85)

    data = {
        "i_liquid_breeder": 1,
        "temp_coolant_pump_outlet": 573,
        "temp_coolant_pump_inlet": 773,
        "pres_coolant_pump_inlet": 8000000,
        "dpres_coolant": 20088.23,
        "mflow_coolant_total": 956.3,
        "primary_coolant_switch": "Helium",
        "den_coolant": 5.64,
        "label": "First Wall and Blanket",
    }

    assert (
        pytest.approx(
            blanket_library_fixture.coolant_pumping_power(False, **data), rel=1e-4
        )
        == 3.2374845432302464
    )


class ComponentHalfHeightParam(NamedTuple):
    z_tf_inside_half: Any = None
    dz_xpoint_divertor: Any = None
    dz_shld_vv_gap: Any = None
    dz_blkt_upper: Any = None
    dz_shld_upper: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_outboard: Any = None
    dz_vv_lower: Any = None
    dz_vv_upper: Any = None
    z_plasma_xpoint_lower: Any = None
    z_plasma_xpoint_upper: Any = None
    n_divertors: Any = None
    dz_divertor: Any = None
    expected_icomponent: Any = None
    expected_half_height: Any = None


@pytest.mark.parametrize(
    "componenthalfheightparam",
    [
        ComponentHalfHeightParam(
            z_tf_inside_half=8.8182171641274945,
            dz_xpoint_divertor=2.0018838307941582,
            dz_shld_vv_gap=0.16300000000000001,
            dz_blkt_upper=0.85000000000000009,
            dz_shld_upper=0.59999999999999998,
            dr_fw_plasma_gap_inboard=0.25,
            dr_fw_plasma_gap_outboard=0.25,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_outboard=0.018000000000000002,
            dz_vv_lower=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            z_plasma_xpoint_lower=4.93333333333333333,
            z_plasma_xpoint_upper=4.93333333333333333,
            n_divertors=1,
            dz_divertor=0.62000000000000011,
            expected_half_height=5.9532752487304119,
        ),
    ],
)
def test_calculate_blkt_half_height(componenthalfheightparam, blanket_library_fixture):
    """
    Regression Unit Test for component_half_height.

    This test was generated using data from large_tokamak.IN.DAT.

    :param componenthalfheightparam: the data used in this test.
    :type componenthalfheightparam: componenthalfheightparam
    """
    half_height = blanket_library_fixture.calculate_blkt_half_height(
        z_plasma_xpoint_lower=componenthalfheightparam.z_plasma_xpoint_lower,
        dz_xpoint_divertor=componenthalfheightparam.dz_xpoint_divertor,
        dz_divertor=componenthalfheightparam.dz_divertor,
        z_plasma_xpoint_upper=componenthalfheightparam.z_plasma_xpoint_upper,
        dr_fw_plasma_gap_inboard=componenthalfheightparam.dr_fw_plasma_gap_inboard,
        dr_fw_plasma_gap_outboard=componenthalfheightparam.dr_fw_plasma_gap_outboard,
        dr_fw_inboard=componenthalfheightparam.dr_fw_inboard,
        dr_fw_outboard=componenthalfheightparam.dr_fw_outboard,
        dz_blkt_upper=componenthalfheightparam.dz_blkt_upper,
        n_divertors=componenthalfheightparam.n_divertors,
    )

    assert half_height == pytest.approx(componenthalfheightparam.expected_half_height)


class ApplyCoverageFactorsParam(NamedTuple):
    a_blkt_outboard_surface: Any = None
    a_blkt_outboard_surface_full_coverage: Any = None
    a_blkt_total_surface: Any = None
    a_blkt_total_surface_full_coverage: Any = None
    a_blkt_inboard_surface: Any = None
    a_blkt_inboard_surface_full_coverage: Any = None
    f_ster_div_single: Any = None
    f_a_fw_outboard_hcd: Any = None
    vol_blkt_outboard: Any = None
    vol_blkt_inboard: Any = None
    vol_blkt_inboard_full_coverage: Any = None
    vol_blkt_total: Any = None
    vol_blkt_total_full_coverage: Any = None
    fvolsi: Any = None
    fvolso: Any = None
    vol_shld_total: Any = None
    n_divertors: Any = None
    vol_shld_inboard: Any = None
    vol_shld_outboard: Any = None
    expected_a_blkt_outboard_surface: Any = None
    expected_a_blkt_total_surface: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_volblkt: Any = None
    expected_vol_vv: Any = None


@pytest.mark.parametrize(
    "applycoveragefactorsparam",
    [
        ApplyCoverageFactorsParam(
            a_blkt_outboard_surface=1101.3666396424403,
            a_blkt_outboard_surface_full_coverage=1101.3666396424403,
            a_blkt_total_surface=1766.3354109399943,
            a_blkt_total_surface_full_coverage=1766.3354109399943,
            a_blkt_inboard_surface=664.9687712975541,
            a_blkt_inboard_surface_full_coverage=664.9687712975541,
            f_ster_div_single=0.115,
            f_a_fw_outboard_hcd=0,
            vol_blkt_outboard=1020.3677420460117,
            vol_blkt_inboard=315.83946385183026,
            vol_blkt_inboard_full_coverage=315.83946385183026,
            vol_blkt_total=1336.207205897842,
            vol_blkt_total_full_coverage=1336.207205897842,
            fvolsi=1,
            fvolso=0.64000000000000001,
            n_divertors=1,
            expected_a_blkt_outboard_surface=898.23806738434075,
            expected_a_blkt_total_surface=1563.2068386818949,
            expected_vol_blkt_outboard=866.70391336775992,
            expected_volblkt=1182.5433772195902,
            expected_vol_vv=1016.2876250857248,
        ),
    ],
)
def test_apply_coverage_factors(
    applycoveragefactorsparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for apply_coverage_factors.

    This test was generated using data from tests/regression/input_files/large_tokamak_eval.IN.DAT.

    :param applycoveragefactorsparam: the data used to mock and assert in this test.
    :type applycoveragefactorsparam: applycoveragefactorsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_outboard_surface",
        applycoveragefactorsparam.a_blkt_outboard_surface,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_outboard_surface_full_coverage",
        applycoveragefactorsparam.a_blkt_outboard_surface_full_coverage,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_total_surface",
        applycoveragefactorsparam.a_blkt_total_surface,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_total_surface_full_coverage",
        applycoveragefactorsparam.a_blkt_total_surface_full_coverage,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_inboard_surface",
        applycoveragefactorsparam.a_blkt_inboard_surface,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "a_blkt_inboard_surface_full_coverage",
        applycoveragefactorsparam.a_blkt_inboard_surface_full_coverage,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "f_ster_div_single",
        applycoveragefactorsparam.f_ster_div_single,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "f_a_fw_outboard_hcd",
        applycoveragefactorsparam.f_a_fw_outboard_hcd,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "vol_blkt_outboard",
        applycoveragefactorsparam.vol_blkt_outboard,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "vol_blkt_inboard",
        applycoveragefactorsparam.vol_blkt_inboard,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "vol_blkt_inboard_full_coverage",
        applycoveragefactorsparam.vol_blkt_inboard_full_coverage,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "vol_blkt_total",
        applycoveragefactorsparam.vol_blkt_total,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "vol_blkt_total_full_coverage",
        applycoveragefactorsparam.vol_blkt_total_full_coverage,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "fvolsi", applycoveragefactorsparam.fvolsi
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "fvolso", applycoveragefactorsparam.fvolso
    )
    monkeypatch.setattr(
        divertor_variables, "n_divertors", applycoveragefactorsparam.n_divertors
    )

    blanket_library_fixture.apply_coverage_factors()

    assert blanket_library_fixture.data.build.a_blkt_outboard_surface == pytest.approx(
        applycoveragefactorsparam.expected_a_blkt_outboard_surface
    )
    assert blanket_library_fixture.data.build.a_blkt_total_surface == pytest.approx(
        applycoveragefactorsparam.expected_a_blkt_total_surface
    )
    assert blanket_library_fixture.data.fwbs.vol_blkt_outboard == pytest.approx(
        applycoveragefactorsparam.expected_vol_blkt_outboard
    )
    assert blanket_library_fixture.data.fwbs.vol_blkt_total == pytest.approx(
        applycoveragefactorsparam.expected_volblkt
    )


class DshapedInboardBlktSegmentParam(NamedTuple):
    dz_blkt_half: Any = None
    n_blkt_inboard_modules_poloidal: Any = None
    expected_len_blkt_inboard_segment_poloidal: Any = None


@pytest.mark.parametrize(
    "dshaped_inboard_param",
    [
        DshapedInboardBlktSegmentParam(
            dz_blkt_half=8.25,
            n_blkt_inboard_modules_poloidal=7,
            expected_len_blkt_inboard_segment_poloidal=2.3571428571428572,
        ),
    ],
)
def test_calculate_dshaped_inboard_blkt_segment_poloidal(
    dshaped_inboard_param, blanket_library_fixture
):
    """Test for calculate_dshaped_inboard_blkt_segment_poloidal."""
    result = blanket_library_fixture.calculate_dshaped_inboard_blkt_segment_poloidal(
        dz_blkt_half=dshaped_inboard_param.dz_blkt_half,
        n_blkt_inboard_modules_poloidal=dshaped_inboard_param.n_blkt_inboard_modules_poloidal,
    )
    assert result == pytest.approx(
        dshaped_inboard_param.expected_len_blkt_inboard_segment_poloidal
    )


class DshapedOutboardBlktSegmentParam(NamedTuple):
    n_blkt_outboard_modules_poloidal: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    rminor: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dz_blkt_half: Any = None
    n_divertors: Any = None
    f_ster_div_single: Any = None
    expected_len_blkt_outboard_segment_poloidal: Any = None


@pytest.mark.parametrize(
    "dshaped_outboard_param",
    [
        DshapedOutboardBlktSegmentParam(
            n_blkt_outboard_modules_poloidal=8,
            dr_fw_plasma_gap_inboard=0.10000000000000001,
            rminor=2.5,
            dr_fw_plasma_gap_outboard=0.10000000000000001,
            dz_blkt_half=8.25,
            n_divertors=2,
            f_ster_div_single=0.115,
            expected_len_blkt_outboard_segment_poloidal=2.0597205347177807,
        ),
    ],
)
def test_calculate_dshaped_outboard_blkt_segment_poloidal(
    dshaped_outboard_param, blanket_library_fixture
):
    """Test for calculate_dshaped_outboard_blkt_segment_poloidal."""
    result = blanket_library_fixture.calculate_dshaped_outboard_blkt_segment_poloidal(
        n_blkt_outboard_modules_poloidal=dshaped_outboard_param.n_blkt_outboard_modules_poloidal,
        dr_fw_plasma_gap_inboard=dshaped_outboard_param.dr_fw_plasma_gap_inboard,
        rminor=dshaped_outboard_param.rminor,
        dr_fw_plasma_gap_outboard=dshaped_outboard_param.dr_fw_plasma_gap_outboard,
        dz_blkt_half=dshaped_outboard_param.dz_blkt_half,
        n_divertors=dshaped_outboard_param.n_divertors,
        f_ster_div_single=dshaped_outboard_param.f_ster_div_single,
    )
    assert result == pytest.approx(
        dshaped_outboard_param.expected_len_blkt_outboard_segment_poloidal
    )


class EllipticalInboardBlktSegmentParam(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    dz_blkt_half: Any = None
    n_blkt_inboard_modules_poloidal: Any = None
    n_divertors: Any = None
    f_ster_div_single: Any = None
    expected_len_blkt_inboard_segment_poloidal: Any = None


@pytest.mark.parametrize(
    "elliptical_inboard_param",
    [
        EllipticalInboardBlktSegmentParam(
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            dr_fw_plasma_gap_inboard=0.25,
            dz_blkt_half=5.9532752487304119,
            n_blkt_inboard_modules_poloidal=7,
            n_divertors=1,
            f_ster_div_single=0.115,
            expected_len_blkt_inboard_segment_poloidal=1.6252823720672551,
        ),
    ],
)
def test_calculate_elliptical_inboard_blkt_segment_poloidal(
    elliptical_inboard_param, blanket_library_fixture
):
    """Test for calculate_elliptical_inboard_blkt_segment_poloidal."""
    result = blanket_library_fixture.calculate_elliptical_inboard_blkt_segment_poloidal(
        rmajor=elliptical_inboard_param.rmajor,
        rminor=elliptical_inboard_param.rminor,
        triang=elliptical_inboard_param.triang,
        dr_fw_plasma_gap_inboard=elliptical_inboard_param.dr_fw_plasma_gap_inboard,
        dz_blkt_half=elliptical_inboard_param.dz_blkt_half,
        n_blkt_inboard_modules_poloidal=elliptical_inboard_param.n_blkt_inboard_modules_poloidal,
        n_divertors=elliptical_inboard_param.n_divertors,
        f_ster_div_single=elliptical_inboard_param.f_ster_div_single,
    )
    assert result == pytest.approx(
        elliptical_inboard_param.expected_len_blkt_inboard_segment_poloidal
    )


class EllipticalOutboardBlktSegmentParam(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    dz_blkt_half: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    n_blkt_outboard_modules_poloidal: Any = None
    n_divertors: Any = None
    f_ster_div_single: Any = None
    expected_len_blkt_outboard_segment_poloidal: Any = None


@pytest.mark.parametrize(
    "elliptical_outboard_param",
    [
        EllipticalOutboardBlktSegmentParam(
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            dz_blkt_half=5.9532752487304119,
            dr_fw_plasma_gap_outboard=0.25,
            n_blkt_outboard_modules_poloidal=8,
            n_divertors=1,
            f_ster_div_single=0.115,
            expected_len_blkt_outboard_segment_poloidal=1.7853902013340495,
        ),
    ],
)
def test_calculate_elliptical_outboard_blkt_segment_poloidal(
    elliptical_outboard_param, blanket_library_fixture
):
    """Test for calculate_elliptical_outboard_blkt_segment_poloidal."""
    result = blanket_library_fixture.calculate_elliptical_outboard_blkt_segment_poloidal(
        rmajor=elliptical_outboard_param.rmajor,
        rminor=elliptical_outboard_param.rminor,
        triang=elliptical_outboard_param.triang,
        dz_blkt_half=elliptical_outboard_param.dz_blkt_half,
        dr_fw_plasma_gap_outboard=elliptical_outboard_param.dr_fw_plasma_gap_outboard,
        n_blkt_outboard_modules_poloidal=elliptical_outboard_param.n_blkt_outboard_modules_poloidal,
        n_divertors=elliptical_outboard_param.n_divertors,
        f_ster_div_single=elliptical_outboard_param.f_ster_div_single,
    )
    assert result == pytest.approx(
        elliptical_outboard_param.expected_len_blkt_outboard_segment_poloidal
    )


class LiquidBreederPropertiesParam(NamedTuple):
    inlet_temp_liq: Any = None
    outlet_temp_liq: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    den_liq: Any = None
    specific_heat_liq: Any = None
    thermal_conductivity_liq: Any = None
    dynamic_viscosity_liq: Any = None
    electrical_conductivity_liq: Any = None
    i_blkt_liquid_breeder_type: Any = None
    hartmann_liq: Any = None
    b_mag_blkt: Any = None
    i_blkt_inboard: Any = None
    i_blkt_dual_coolant: Any = None
    b_plasma_toroidal_on_axis: Any = None
    aspect: Any = None
    rmajor: Any = None
    dr_blkt_inboard: Any = None
    dr_blkt_outboard: Any = None
    ip: Any = None
    expected_den_liq: Any = None
    expected_specific_heat_liq: Any = None
    expected_thermal_conductivity_liq: Any = None
    expected_dynamic_viscosity_liq: Any = None
    expected_electrical_conductivity_liq: Any = None
    expected_hartmann_liq: Any = None
    expected_b_mag_blkt: Any = None


@pytest.mark.parametrize(
    "liquidbreederpropertiesparam",
    [
        LiquidBreederPropertiesParam(
            inlet_temp_liq=570,
            outlet_temp_liq=720,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            den_liq=9500,
            specific_heat_liq=190,
            thermal_conductivity_liq=30,
            dynamic_viscosity_liq=0,
            electrical_conductivity_liq=0,
            i_blkt_liquid_breeder_type=0,
            hartmann_liq=np.array(
                np.array((0.0, 0.0), order="F"), order="F"
            ).transpose(),
            b_mag_blkt=np.array(np.array((5.0, 5.0), order="F"), order="F").transpose(),
            i_blkt_inboard=1,
            i_blkt_dual_coolant=0,
            b_plasma_toroidal_on_axis=5.7000000000000002,
            aspect=3,
            rmajor=8,
            dr_blkt_inboard=0.70000000000000007,
            dr_blkt_outboard=1,
            ip=0,
            expected_den_liq=9753.2497999999996,
            expected_specific_heat_liq=189.12018,
            expected_thermal_conductivity_liq=9.238260167312621,
            expected_dynamic_viscosity_liq=0.0017477589255667965,
            expected_electrical_conductivity_liq=861800.80431007256,
            expected_hartmann_liq=np.array(
                np.array((20319.245984309102, 9067.8426109080938), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=np.array(
                np.array((9.1505016722408019, 4.0835820895522392), order="F"),
                order="F",
            ).transpose(),
        ),
        LiquidBreederPropertiesParam(
            inlet_temp_liq=570,
            outlet_temp_liq=720,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            den_liq=9500,
            specific_heat_liq=190,
            thermal_conductivity_liq=30,
            expected_thermal_conductivity_liq=30,  # doesn't change when i_blkt_liquid_breeder_type=1
            dynamic_viscosity_liq=0,
            electrical_conductivity_liq=0,
            i_blkt_liquid_breeder_type=1,
            hartmann_liq=np.array(
                np.array((0.0, 0.0), order="F"), order="F"
            ).transpose(),
            b_mag_blkt=np.array(np.array((5.0, 5.0), order="F"), order="F").transpose(),
            i_blkt_inboard=1,
            i_blkt_dual_coolant=0,
            b_plasma_toroidal_on_axis=5.7000000000000002,
            aspect=3,
            rmajor=8,
            dr_blkt_inboard=0.70000000000000007,
            dr_blkt_outboard=1,
            ip=0,
            expected_den_liq=305.30702851374997,
            expected_specific_heat_liq=34.640761200690406,
            expected_dynamic_viscosity_liq=0.00037298826343426359,
            expected_electrical_conductivity_liq=596562356750.5,
            expected_hartmann_liq=np.array(
                np.array((36595294.326740541, 16331332.841336451), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=np.array(
                np.array((9.1505016722408019, 4.0835820895522392), order="F"),
                order="F",
            ).transpose(),
        ),
    ],
)
def test_liquid_breeder_properties(
    liquidbreederpropertiesparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for liquid_breeder_properties.

    This test was generated using data from blanket_files/large_tokamak_dcll.IN.DAT.

    :param liquidbreederpropertiesparam: the data used to mock and assert in this test.
    :type liquidbreederpropertiesparam: liquidbreederpropertiesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "inlet_temp_liq",
        liquidbreederpropertiesparam.inlet_temp_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "outlet_temp_liq",
        liquidbreederpropertiesparam.outlet_temp_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "a_bz_liq",
        liquidbreederpropertiesparam.a_bz_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "b_bz_liq",
        liquidbreederpropertiesparam.b_bz_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "den_liq",
        liquidbreederpropertiesparam.den_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "specific_heat_liq",
        liquidbreederpropertiesparam.specific_heat_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "thermal_conductivity_liq",
        liquidbreederpropertiesparam.thermal_conductivity_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "dynamic_viscosity_liq",
        liquidbreederpropertiesparam.dynamic_viscosity_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "electrical_conductivity_liq",
        liquidbreederpropertiesparam.electrical_conductivity_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_liquid_breeder_type",
        liquidbreederpropertiesparam.i_blkt_liquid_breeder_type,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "hartmann_liq",
        liquidbreederpropertiesparam.hartmann_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "b_mag_blkt",
        liquidbreederpropertiesparam.b_mag_blkt,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_inboard",
        liquidbreederpropertiesparam.i_blkt_inboard,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_dual_coolant",
        liquidbreederpropertiesparam.i_blkt_dual_coolant,
    )
    monkeypatch.setattr(
        physics_variables,
        "b_plasma_toroidal_on_axis",
        liquidbreederpropertiesparam.b_plasma_toroidal_on_axis,
    )
    monkeypatch.setattr(physics_variables, "aspect", liquidbreederpropertiesparam.aspect)
    monkeypatch.setattr(physics_variables, "rmajor", liquidbreederpropertiesparam.rmajor)
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "dr_blkt_inboard",
        liquidbreederpropertiesparam.dr_blkt_inboard,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.build,
        "dr_blkt_outboard",
        liquidbreederpropertiesparam.dr_blkt_outboard,
    )

    blanket_library_fixture.liquid_breeder_properties()

    assert blanket_library_fixture.data.fwbs.den_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_den_liq
    )
    assert blanket_library_fixture.data.fwbs.specific_heat_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_specific_heat_liq
    )
    assert blanket_library_fixture.data.fwbs.thermal_conductivity_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_thermal_conductivity_liq
    )
    assert blanket_library_fixture.data.fwbs.dynamic_viscosity_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_dynamic_viscosity_liq
    )
    assert (
        blanket_library_fixture.data.fwbs.electrical_conductivity_liq
        == pytest.approx(
            liquidbreederpropertiesparam.expected_electrical_conductivity_liq
        )
    )
    assert blanket_library_fixture.data.fwbs.hartmann_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_hartmann_liq
    )
    assert blanket_library_fixture.data.fwbs.b_mag_blkt == pytest.approx(
        liquidbreederpropertiesparam.expected_b_mag_blkt
    )


class PressureDropParam(NamedTuple):
    radius_fw_channel: Any = None
    radius_pipe_90_deg_bend: Any = None
    radius_pipe_180_deg_bend: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    roughness_fw_channel: Any = None
    ip: Any = None
    i_ps: Any = None
    num_90: Any = None
    num_180: Any = None
    l_pipe: Any = None
    den: Any = None
    vsc: Any = None
    vv: Any = None
    label: Any = None
    expected_pressure_drop_out: Any = None


@pytest.mark.parametrize(
    "pressuredropparam",
    [
        PressureDropParam(
            radius_fw_channel=0.0060000000000000001,
            radius_pipe_90_deg_bend=0.018,
            radius_pipe_180_deg_bend=0.09,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            roughness_fw_channel=9.9999999999999995e-07,
            ip=0,
            i_ps=1,
            num_90=2,
            num_180=0,
            l_pipe=4,
            den=10.405276820718059,
            vsc=3.604452999475736e-05,
            vv=32.753134225223164,
            label="Inboard first wall",
            expected_pressure_drop_out=36213.58989742931,
        ),
        PressureDropParam(
            radius_fw_channel=1.0,
            radius_pipe_90_deg_bend=1.0,
            radius_pipe_180_deg_bend=1.0,
            a_bz_liq=1.0,
            b_bz_liq=1.0,
            roughness_fw_channel=1e-6,
            ip=0,
            i_ps=2,
            num_90=1.0,
            num_180=1.0,
            l_pipe=1.0,
            den=1.0,
            vsc=1.0,
            vv=1.0,
            label="label",
            expected_pressure_drop_out=1.4325633520224854,
        ),
    ],
)
def test_pressure_drop(pressuredropparam, monkeypatch, blanket_library_fixture):
    """
    Automatically generated Regression Unit Test for pressure_drop.

    This test was generated using data from blanket_files/large_tokamak_primary_pumping2.IN.DAT.

    :param pressuredropparam: the data used to mock and assert in this test.
    :type pressuredropparam: pressuredropparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "radius_fw_channel",
        pressuredropparam.radius_fw_channel,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "a_bz_liq", pressuredropparam.a_bz_liq
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "b_bz_liq", pressuredropparam.b_bz_liq
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "roughness_fw_channel",
        pressuredropparam.roughness_fw_channel,
    )

    pressure_drop_out = blanket_library_fixture.coolant_friction_pressure_drop(
        i_ps=pressuredropparam.i_ps,
        radius_pipe_90_deg_bend=pressuredropparam.radius_pipe_90_deg_bend,
        radius_pipe_180_deg_bend=pressuredropparam.radius_pipe_180_deg_bend,
        n_pipe_90_deg_bends=pressuredropparam.num_90,
        n_pipe_180_deg_bends=pressuredropparam.num_180,
        len_pipe=pressuredropparam.l_pipe,
        den_coolant=pressuredropparam.den,
        visc_coolant=pressuredropparam.vsc,
        vel_coolant=pressuredropparam.vv,
        label=pressuredropparam.label,
    )

    assert pressure_drop_out == pytest.approx(
        pressuredropparam.expected_pressure_drop_out
    )


class LiquidBreederPressureDropMhdParam(NamedTuple):
    i_blkt_liquid_breeder_channel_type: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    b_mag_blkt: Any = None
    bz_channel_conduct_liq: Any = None
    th_wall_secondary: Any = None
    ip: Any = None
    vel: Any = None
    vsc: Any = None
    conduct_liq: Any = None
    l_channel: Any = None
    num_pol: Any = None
    label: Any = None
    expected_liquid_breeder_pressure_drop_mhd_out: Any = None


@pytest.mark.parametrize(
    "liquidbreederpressuredropmhdparam",
    [
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=0,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=282697824.60502106,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=1,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=160029.28473931071,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=282697824.60502106,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            vel=0.042500391943592931,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.7176027768600395,
            num_pol=1792,
            label="Inboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=487177576.30010164,
        ),
    ],
)
def test_liquid_breeder_pressure_drop_mhd(
    liquidbreederpressuredropmhdparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for liquid_breeder_pressure_drop_mhd.

    This test was generated using data from blanket_files/dcll_mms_lt_IN.DAT.

    :param liquidbreederpressuredropmhdparam: the data used to mock and assert in this test.
    :type liquidbreederpressuredropmhdparam: liquidbreederpressuredropmhdparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "i_blkt_liquid_breeder_channel_type",
        liquidbreederpressuredropmhdparam.i_blkt_liquid_breeder_channel_type,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "a_bz_liq",
        liquidbreederpressuredropmhdparam.a_bz_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "b_bz_liq",
        liquidbreederpressuredropmhdparam.b_bz_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "b_mag_blkt",
        liquidbreederpressuredropmhdparam.b_mag_blkt,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "bz_channel_conduct_liq",
        liquidbreederpressuredropmhdparam.bz_channel_conduct_liq,
    )
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs,
        "th_wall_secondary",
        liquidbreederpressuredropmhdparam.th_wall_secondary,
    )

    liquid_breeder_pressure_drop_mhd_out = (
        blanket_library_fixture.liquid_breeder_mhd_pressure_drop(
            vel=liquidbreederpressuredropmhdparam.vel,
            vsc=liquidbreederpressuredropmhdparam.vsc,
            conduct_liq=liquidbreederpressuredropmhdparam.conduct_liq,
            l_channel=liquidbreederpressuredropmhdparam.l_channel,
            num_pol=liquidbreederpressuredropmhdparam.num_pol,
            label=liquidbreederpressuredropmhdparam.label,
        )
    )

    assert liquid_breeder_pressure_drop_mhd_out == pytest.approx(
        liquidbreederpressuredropmhdparam.expected_liquid_breeder_pressure_drop_mhd_out
    )


class CalculateDshapedBlktAreasParam(NamedTuple):
    r_shld_inboard_inner: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    rminor: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_outboard: Any = None
    dz_blkt_half: Any = None
    expected_a_blkt_inboard_surface: Any = None
    expected_a_blkt_outboard_surface: Any = None
    expected_a_blkt_total_surface: Any = None


@pytest.mark.parametrize(
    "calculatedshapedblktareasparam",
    [
        CalculateDshapedBlktAreasParam(
            r_shld_inboard_inner=1.5,
            dr_shld_inboard=0.4,
            dr_blkt_inboard=0.0,
            dr_fw_inboard=0.018,
            dr_fw_plasma_gap_inboard=0.1,
            rminor=2.5,
            dr_fw_plasma_gap_outboard=0.1,
            dr_fw_outboard=0.018,
            dz_blkt_half=8.25,
            expected_a_blkt_inboard_surface=196.97785938008002,
            expected_a_blkt_outboard_surface=852.24160940262459,
            expected_a_blkt_total_surface=1049.2194687827046,
        ),
    ],
)
def test_calculate_dshaped_blkt_areas(
    calculatedshapedblktareasparam, blanket_library_fixture
):
    """
    Regression Unit Test for calculate_dshaped_blkt_areas.

    This test was generated using data from tests/regression/input_files/st_regression.IN.DAT.

    :param calculatedshapedblktareasparam: the data used in this test.
    :type calculatedshapedblktareasparam: CalculateDshapedBlktAreasParam
    """
    (
        a_blkt_inboard_surface,
        a_blkt_outboard_surface,
        a_blkt_total_surface,
    ) = blanket_library_fixture.calculate_dshaped_blkt_areas(
        r_shld_inboard_inner=calculatedshapedblktareasparam.r_shld_inboard_inner,
        dr_shld_inboard=calculatedshapedblktareasparam.dr_shld_inboard,
        dr_blkt_inboard=calculatedshapedblktareasparam.dr_blkt_inboard,
        dr_fw_inboard=calculatedshapedblktareasparam.dr_fw_inboard,
        dr_fw_plasma_gap_inboard=calculatedshapedblktareasparam.dr_fw_plasma_gap_inboard,
        rminor=calculatedshapedblktareasparam.rminor,
        dr_fw_plasma_gap_outboard=calculatedshapedblktareasparam.dr_fw_plasma_gap_outboard,
        dr_fw_outboard=calculatedshapedblktareasparam.dr_fw_outboard,
        dz_blkt_half=calculatedshapedblktareasparam.dz_blkt_half,
    )

    assert a_blkt_inboard_surface == pytest.approx(
        calculatedshapedblktareasparam.expected_a_blkt_inboard_surface
    )
    assert a_blkt_outboard_surface == pytest.approx(
        calculatedshapedblktareasparam.expected_a_blkt_outboard_surface
    )
    assert a_blkt_total_surface == pytest.approx(
        calculatedshapedblktareasparam.expected_a_blkt_total_surface
    )


class CalculateDshapedBlktVolumesParam(NamedTuple):
    r_shld_inboard_inner: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    rminor: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_outboard: Any = None
    dz_blkt_half: Any = None
    dr_blkt_outboard: Any = None
    dz_blkt_upper: Any = None
    expected_vol_blkt_inboard: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_vol_blkt_total: Any = None


@pytest.mark.parametrize(
    "calculatedshapedblktvolumesparam",
    [
        CalculateDshapedBlktVolumesParam(
            r_shld_inboard_inner=1.5,
            dr_shld_inboard=0.4,
            dr_blkt_inboard=0.6,
            dr_fw_inboard=0.018,
            dr_fw_plasma_gap_inboard=0.1,
            rminor=2.5,
            dr_fw_plasma_gap_outboard=0.1,
            dr_fw_outboard=0.018,
            dz_blkt_half=8.25,
            dr_blkt_outboard=1.0,
            dz_blkt_upper=0.85,
            expected_vol_blkt_inboard=150.94724381968237,
            expected_vol_blkt_outboard=869.2500537130913,
            expected_vol_blkt_total=1020.1972975327737,
        ),
    ],
)
def test_calculate_dshaped_blkt_volumes(
    calculatedshapedblktvolumesparam, blanket_library_fixture
):
    """
    Regression Unit Test for calculate_dshaped_blkt_volumes.

    This test was generated using data from tests/regression/input_files/st_regression.IN.DAT.

    :param calculatedshapedblktvolumesparam: the data used in this test.
    :type calculatedshapedblktvolumesparam: CalculateDshapedBlktVolumesParam
    """
    (
        vol_blkt_inboard,
        vol_blkt_outboard,
        vol_blkt_total,
    ) = blanket_library_fixture.calculate_dshaped_blkt_volumes(
        r_shld_inboard_inner=calculatedshapedblktvolumesparam.r_shld_inboard_inner,
        dr_shld_inboard=calculatedshapedblktvolumesparam.dr_shld_inboard,
        dr_blkt_inboard=calculatedshapedblktvolumesparam.dr_blkt_inboard,
        dr_fw_inboard=calculatedshapedblktvolumesparam.dr_fw_inboard,
        dr_fw_plasma_gap_inboard=calculatedshapedblktvolumesparam.dr_fw_plasma_gap_inboard,
        rminor=calculatedshapedblktvolumesparam.rminor,
        dr_fw_plasma_gap_outboard=calculatedshapedblktvolumesparam.dr_fw_plasma_gap_outboard,
        dr_fw_outboard=calculatedshapedblktvolumesparam.dr_fw_outboard,
        dz_blkt_half=calculatedshapedblktvolumesparam.dz_blkt_half,
        dr_blkt_outboard=calculatedshapedblktvolumesparam.dr_blkt_outboard,
        dz_blkt_upper=calculatedshapedblktvolumesparam.dz_blkt_upper,
    )

    assert vol_blkt_inboard == pytest.approx(
        calculatedshapedblktvolumesparam.expected_vol_blkt_inboard
    )
    assert vol_blkt_outboard == pytest.approx(
        calculatedshapedblktvolumesparam.expected_vol_blkt_outboard
    )
    assert vol_blkt_total == pytest.approx(
        calculatedshapedblktvolumesparam.expected_vol_blkt_total
    )


class CalculateEllipticalBlktAreasParam(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    r_shld_inboard_inner: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    r_shld_outboard_outer: Any = None
    dr_shld_outboard: Any = None
    dr_blkt_outboard: Any = None
    dz_blkt_half: Any = None
    expected_a_blkt_inboard_surface: Any = None
    expected_a_blkt_outboard_surface: Any = None
    expected_a_blkt_total_surface: Any = None


@pytest.mark.parametrize(
    "calculateellipticalblktareasparam",
    [
        CalculateEllipticalBlktAreasParam(
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            r_shld_inboard_inner=4.0833333333333339,
            dr_shld_inboard=0.30000000000000004,
            dr_blkt_inboard=0.70000000000000007,
            r_shld_outboard_outer=12.716666666666667,
            dr_shld_outboard=0.80000000000000004,
            dr_blkt_outboard=1,
            dz_blkt_half=5.9532752487304119,
            expected_a_blkt_inboard_surface=664.9687712975541,
            expected_a_blkt_outboard_surface=1101.3666396424403,
            expected_a_blkt_total_surface=1766.3354109399943,
        ),
    ],
)
def test_calculate_elliptical_blkt_areas(
    calculateellipticalblktareasparam, blanket_library_fixture
):
    """
    Regression Unit Test for calculate_elliptical_blkt_areas.

    This test was generated using data from tests/regression/input_files/large_tokamak_eval.IN.DAT.

    :param calculateellipticalblktareasparam: the data used in this test.
    :type calculateellipticalblktareasparam: CalculateEllipticalBlktAreasParam
    """
    (
        a_blkt_inboard_surface,
        a_blkt_outboard_surface,
        a_blkt_total_surface,
    ) = blanket_library_fixture.calculate_elliptical_blkt_areas(
        rmajor=calculateellipticalblktareasparam.rmajor,
        rminor=calculateellipticalblktareasparam.rminor,
        triang=calculateellipticalblktareasparam.triang,
        r_shld_inboard_inner=calculateellipticalblktareasparam.r_shld_inboard_inner,
        dr_shld_inboard=calculateellipticalblktareasparam.dr_shld_inboard,
        dr_blkt_inboard=calculateellipticalblktareasparam.dr_blkt_inboard,
        r_shld_outboard_outer=calculateellipticalblktareasparam.r_shld_outboard_outer,
        dr_shld_outboard=calculateellipticalblktareasparam.dr_shld_outboard,
        dr_blkt_outboard=calculateellipticalblktareasparam.dr_blkt_outboard,
        dz_blkt_half=calculateellipticalblktareasparam.dz_blkt_half,
    )

    assert a_blkt_inboard_surface == pytest.approx(
        calculateellipticalblktareasparam.expected_a_blkt_inboard_surface
    )
    assert a_blkt_outboard_surface == pytest.approx(
        calculateellipticalblktareasparam.expected_a_blkt_outboard_surface
    )
    assert a_blkt_total_surface == pytest.approx(
        calculateellipticalblktareasparam.expected_a_blkt_total_surface
    )


class CalculateEllipticalBlktVolumesParam(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    r_shld_inboard_inner: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    r_shld_outboard_outer: Any = None
    dr_shld_outboard: Any = None
    dr_blkt_outboard: Any = None
    dz_blkt_half: Any = None
    dz_blkt_upper: Any = None
    expected_vol_blkt_inboard: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_vol_blkt_total: Any = None


@pytest.mark.parametrize(
    "calculateellipticalblktvolumesparam",
    [
        CalculateEllipticalBlktVolumesParam(
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            r_shld_inboard_inner=4.0833333333333339,
            dr_shld_inboard=0.30000000000000004,
            dr_blkt_inboard=0.70000000000000007,
            r_shld_outboard_outer=12.716666666666667,
            dr_shld_outboard=0.80000000000000004,
            dr_blkt_outboard=1,
            dz_blkt_half=5.9532752487304119,
            dz_blkt_upper=0.85000000000000009,
            expected_vol_blkt_inboard=315.83946385183026,
            expected_vol_blkt_outboard=1020.3677420460117,
            expected_vol_blkt_total=1336.207205897842,
        ),
    ],
)
def test_calculate_elliptical_blkt_volumes(
    calculateellipticalblktvolumesparam, blanket_library_fixture
):
    """
    Regression Unit Test for calculate_elliptical_blkt_volumes.

    This test was generated using data from tests/regression/input_files/large_tokamak_eval.IN.DAT.

    :param calculateellipticalblktvolumesparam: the data used in this test.
    :type calculateellipticalblktvolumesparam: CalculateEllipticalBlktVolumesParam
    """
    (
        vol_blkt_inboard,
        vol_blkt_outboard,
        vol_blkt_total,
    ) = blanket_library_fixture.calculate_elliptical_blkt_volumes(
        rmajor=calculateellipticalblktvolumesparam.rmajor,
        rminor=calculateellipticalblktvolumesparam.rminor,
        triang=calculateellipticalblktvolumesparam.triang,
        r_shld_inboard_inner=calculateellipticalblktvolumesparam.r_shld_inboard_inner,
        dr_shld_inboard=calculateellipticalblktvolumesparam.dr_shld_inboard,
        dr_blkt_inboard=calculateellipticalblktvolumesparam.dr_blkt_inboard,
        r_shld_outboard_outer=calculateellipticalblktvolumesparam.r_shld_outboard_outer,
        dr_shld_outboard=calculateellipticalblktvolumesparam.dr_shld_outboard,
        dr_blkt_outboard=calculateellipticalblktvolumesparam.dr_blkt_outboard,
        dz_blkt_half=calculateellipticalblktvolumesparam.dz_blkt_half,
        dz_blkt_upper=calculateellipticalblktvolumesparam.dz_blkt_upper,
    )

    assert vol_blkt_inboard == pytest.approx(
        calculateellipticalblktvolumesparam.expected_vol_blkt_inboard
    )
    assert vol_blkt_outboard == pytest.approx(
        calculateellipticalblktvolumesparam.expected_vol_blkt_outboard
    )
    assert vol_blkt_total == pytest.approx(
        calculateellipticalblktvolumesparam.expected_vol_blkt_total
    )


def test_hydraulic_diameter(monkeypatch, blanket_library_fixture):
    """
    Test for hydraulic_diameter function.
    """
    # Set var values
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "radius_fw_channel", 1.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 1.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "b_bz_liq", 1.0)

    # hydraulic_diameter input = i_channel_shape: 1 = circle, 2 = rectangle
    # 2.0D0*radius_fw_channel
    assert blanket_library_fixture.pipe_hydraulic_diameter(1) == pytest.approx(2.0)
    # 2*a_bz_liq*b_bz_liq/(a_bz_liq+b_bz_liq)
    assert blanket_library_fixture.pipe_hydraulic_diameter(2) == pytest.approx(1.0)


def test_elbow_coeff(blanket_library_fixture):
    """
    Test for elbow_coeff function.
    """
    # input = r_elbow, ang_elbow, lambda, dh
    assert blanket_library_fixture.elbow_coeff(1, 0, 1, 1) == pytest.approx(
        0.0, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 1, 1) == pytest.approx(
        1.7807963267948965, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 180, 1, 1) == pytest.approx(
        3.291157766597427, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 1, 0.1) == pytest.approx(
        15.774371098812502, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(0.1, 90, 1, 1) == pytest.approx(
        66.57, rel=1e-3
    )
    assert blanket_library_fixture.elbow_coeff(1, 90, 0.1, 1) == pytest.approx(
        0.3670796326794896, rel=1e-3
    )


def test_flow_velocity(monkeypatch, blanket_library_fixture):
    """
    Test for flow_velocity function.
    """
    # Set var values
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "radius_fw_channel", 1.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 1.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "b_bz_liq", 1.0)

    # input = i_channel_shape, mass_flow_rate, flow_density
    assert blanket_library_fixture.flow_velocity(1, 1, 1) == pytest.approx(
        0.318, rel=1e-3
    )
    assert blanket_library_fixture.flow_velocity(2, 1, 1) == pytest.approx(1.0)
    assert blanket_library_fixture.flow_velocity(1, 0, 1) == pytest.approx(0.0)
    assert blanket_library_fixture.flow_velocity(2, 0, 1) == pytest.approx(0.0)


def test_liquid_breeder_properties_part_1(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure.
    PbLi or Li, with inboard blanket, no inlet/outlet temp difference.
    """
    # Set var values
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(physics_variables, "b_plasma_toroidal_on_axis", 6.0)
    monkeypatch.setattr(physics_variables, "rmajor", 8.0)
    monkeypatch.setattr(physics_variables, "aspect", 3.0)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_inboard", 0.1)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_outboard", 0.2)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "i_blkt_inboard", 1)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "inlet_temp_liq", 1.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "outlet_temp_liq", 1.0)

    # PbLi - see [Fer2020] for relavent equations
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_type", 0
    )

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(blanket_library_fixture.data.fwbs.den_liq, rel=1e-3) == 1.052e4
    assert (
        pytest.approx(blanket_library_fixture.data.fwbs.specific_heat_liq, rel=1e-3)
        == 195.0
    )
    assert (
        pytest.approx(
            blanket_library_fixture.data.fwbs.thermal_conductivity_liq, rel=1e-3
        )
        == -3.384
    )
    assert (
        pytest.approx(blanket_library_fixture.data.fwbs.dynamic_viscosity_liq, rel=1e-3)
        == 0.0155
    )
    assert (
        pytest.approx(
            blanket_library_fixture.data.fwbs.electrical_conductivity_liq, rel=1e-3
        )
        == 9.71e5
    )

    assert pytest.approx(blanket_library_fixture.data.fwbs.b_mag_blkt, rel=1e-3) == (
        9.085,
        4.458,
    )
    assert pytest.approx(blanket_library_fixture.data.fwbs.hartmann_liq, rel=1e-3) == (
        7.189e3,
        3.528e3,
    )

    # Li - see [Lyublinski et al., 2009] for relavent equations
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_type", 1
    )

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(blanket_library_fixture.data.fwbs.den_liq, rel=1e-3) == 504.0
    assert (
        pytest.approx(blanket_library_fixture.data.fwbs.specific_heat_liq, rel=1e-3)
        == 2.833e6
    )
    assert (
        pytest.approx(blanket_library_fixture.data.fwbs.dynamic_viscosity_liq, rel=1e-3)
        == 1.051e112
    )
    assert (
        pytest.approx(
            blanket_library_fixture.data.fwbs.electrical_conductivity_liq, rel=1e-3
        )
        == 9.27e8
    )
    assert pytest.approx(blanket_library_fixture.data.fwbs.b_mag_blkt, rel=1e-3) == (
        9.085,
        4.458,
    )
    assert pytest.approx(blanket_library_fixture.data.fwbs.hartmann_liq, rel=1e-3) == (
        2.7e-53,
        1.3e-53,
    )

    # con_vsc_rat = electrical_conductivity_liq/dynamic_viscosity_liq
    # hartmann_liq = b_mag_blkt * a_bz_liq/2.0D0 * sqrt(con_vsc_rat)


def test_liquid_breeder_properties_part_2(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure. No inboard blanket.
    """
    # Set var values
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(physics_variables, "b_plasma_toroidal_on_axis", 6.0)
    monkeypatch.setattr(physics_variables, "rmajor", 8.0)
    monkeypatch.setattr(physics_variables, "aspect", 3.0)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_inboard", 0.0)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_outboard", 0.2)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "i_blkt_inboard", 0)
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_type", 0
    )
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "inlet_temp_liq", 0.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "outlet_temp_liq", 0.0)

    blanket_library_fixture.liquid_breeder_properties()

    assert pytest.approx(blanket_library_fixture.data.fwbs.b_mag_blkt, rel=1e-3) == (
        8.999,
        4.458,
    )


def test_liquid_breeder_properties_part_3(monkeypatch, blanket_library_fixture):
    """
    Test for liquid_breeder_properties procedure.
    With inlet/outlet temp difference.
    """
    # Set var values
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "a_bz_liq", 0.2)
    monkeypatch.setattr(physics_variables, "b_plasma_toroidal_on_axis", 6.0)
    monkeypatch.setattr(physics_variables, "rmajor", 8.0)
    monkeypatch.setattr(physics_variables, "aspect", 3.0)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_inboard", 0.1)
    monkeypatch.setattr(blanket_library_fixture.data.build, "dr_blkt_outboard", 0.2)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "i_blkt_inboard", 1)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "inlet_temp_liq", 0.0)
    monkeypatch.setattr(blanket_library_fixture.data.fwbs, "outlet_temp_liq", 1.0)

    # PbLi - see [Fer2020] for relavent equations
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_type", 0
    )

    blanket_library_fixture.liquid_breeder_properties()
    assert pytest.approx(blanket_library_fixture.data.fwbs.den_liq, rel=1e-3) == 1.052e4
    # Li - see [Lyublinski et al., 2009] for relavent equations
    monkeypatch.setattr(
        blanket_library_fixture.data.fwbs, "i_blkt_liquid_breeder_type", 1
    )

    blanket_library_fixture.liquid_breeder_properties()
    assert pytest.approx(blanket_library_fixture.data.fwbs.den_liq, rel=1e-3) == 504.0
