from typing import Any, NamedTuple

import numpy as np
import pytest

from process.data_structure import (
    build_variables,
    constraint_variables,
    divertor_variables,
    fwbs_variables,
    global_variables,
    physics_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
)
from process.models.tfcoil import superconducting_tf_coil as sctf
from process.models.tfcoil.superconducting_tf_coil import SuperconductingTFCoil


@pytest.fixture
def sctfcoil():
    """Provides SuperconductingTFCoil object for testing.

    :returns: initialised SuperconductingTFCoil object
    :rtype: process.sctfcoil.SuperconductingTFCoil
    """
    return SuperconductingTFCoil()


class ProtectParam(NamedTuple):
    aio: Any = None

    tfes: Any = None

    acs: Any = None

    aturn: Any = None

    tdump: Any = None

    fcond: Any = None

    fcu: Any = None

    tba: Any = None

    tmax: Any = None

    peak_field: Any = None

    cu_rrr: Any = None

    detection_time: Any = None

    fluence: Any = None

    expected_ajwpro: Any = None

    expected_vd: Any = None


@pytest.mark.parametrize(
    "protectparam",
    (
        ProtectParam(
            aio=74026.751437500003,
            tfes=9561415368.8360519,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            tdump=25.829000000000001,
            fcond=0.63927285511442711,
            fcu=0.80884,
            tba=4.75,
            tmax=150,
            # These are picked to more or less correspond to the previous model
            peak_field=0.0,
            cu_rrr=33.0,
            fluence=0.0,
            detection_time=0.0,
            expected_ajwpro=17785745.149250004,
            expected_vd=10001.287165953383,
        ),
        # Test new model features
        ProtectParam(
            aio=74026.751437500003,
            tfes=9561415368.8360519,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            tdump=25.829000000000001,
            fcond=0.63927285511442711,
            fcu=0.80884,
            tba=4.75,
            tmax=150,
            peak_field=11.0,
            cu_rrr=200.0,
            fluence=3.2e21,
            detection_time=3.0,
            expected_ajwpro=15248071.694109693,
            expected_vd=10001.287165953383,
        ),
    ),
)
def test_protect(protectparam, sctfcoil):
    """
    Automatically generated Regression Unit Test for protect.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param protectparam: the data used to mock and assert in this test.
    :type protectparam: protectparam

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    ajwpro, vd = sctfcoil.quench_heat_protection_current_density(
        c_tf_turn=protectparam.aio,
        e_tf_coil_magnetic_stored=protectparam.tfes,
        a_tf_turn_cable_space=protectparam.acs,
        a_tf_turn=protectparam.aturn,
        t_tf_quench_dump=protectparam.tdump,
        f_a_tf_turn_cable_space_conductor=protectparam.fcond,
        f_a_tf_turn_cable_copper=protectparam.fcu,
        temp_tf_coolant_peak_field=protectparam.tba,
        temp_tf_conductor_quench_max=protectparam.tmax,
        b_tf_inboard_peak=protectparam.peak_field,
        cu_rrr=protectparam.cu_rrr,
        t_tf_quench_detection=protectparam.detection_time,
        nflutfmax=protectparam.fluence,
    )

    assert ajwpro == pytest.approx(protectparam.expected_ajwpro)

    assert vd == pytest.approx(protectparam.expected_vd)


class SuperconParam(NamedTuple):
    temp_tf_superconductor_margin_min: Any = None

    n_tf_coils: Any = None

    temp_margin: Any = None

    dia_tf_turn_coolant_channel: Any = None

    c_tf_turn: Any = None

    b_tf_inboard_peak_with_ripple: Any = None

    str_tf_con_res: Any = None

    a_tf_turn_cable_space_effective: Any = None

    b_crit_upper_nbti: Any = None

    i_str_wp: Any = None

    str_wp: Any = None

    t_crit_nbti: Any = None

    tf_fit_t: Any = None

    tf_fit_z: Any = None

    f_b_tf_inboard_peak_ripple_symmetric: Any = None

    run_tests: Any = None

    i_tf_superconductor: Any = None

    iprint: Any = None

    outfile: Any = None

    a_tf_turn_cable_space: Any = None

    a_tf_turn: Any = None

    b_tf_inboard_peak_symmetric: Any = None

    f_a_tf_turn_cable_copper: Any = None

    f_a_tf_turn_cooling_extra: Any = None

    f_strain_scale: Any = None

    j_tf_wp: Any = None

    t_tf_quench_dump: Any = None

    f_a_tf_turn_cable_space_cooling: Any = None

    e_tf_coil_magnetic_stored: Any = None

    temp_tf_coolant_peak_field: Any = None

    temp_tf_conductor_peak_quench: Any = None

    cu_rrr: Any = None

    detection_time: Any = None

    fluence: Any = None

    bcritsc: Any = None

    tcritsc: Any = None

    expected_temp_margin: Any = None

    expected_j_tf_wp_quench_heat_max: Any = None

    expected_j_tf_wp_critical: Any = None

    expected_vd: Any = None

    expected_j_superconductor_critical: Any = None

    expected_f_c_tf_turn_operating_critical: Any = None

    expected_j_tf_coil_turn: Any = None

    expected_bc20m: Any = None

    expected_tc0m: Any = None

    expected_c_turn_cables_critical: Any = None

    expected_j_superconductor: Any = None


@pytest.mark.parametrize(
    "superconparam",
    (
        SuperconParam(
            temp_tf_superconductor_margin_min=1.5,
            n_tf_coils=16,
            temp_margin=0,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            c_tf_turn=74026.751437500003,
            b_tf_inboard_peak_with_ripple=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            f_a_tf_turn_cable_space_cooling=0.3,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            f_b_tf_inboard_peak_ripple_symmetric=1.0658869305062604,
            run_tests=0,
            i_tf_superconductor=5,
            iprint=0,
            outfile=11,
            a_tf_turn_cable_space=0.001293323051622732,
            a_tf_turn=0.0032012300777680192,
            a_tf_turn_cable_space_effective=0.001,
            b_tf_inboard_peak_symmetric=12.48976756562082,
            f_a_tf_turn_cable_copper=0.80884,
            f_a_tf_turn_cooling_extra=0.30000000000000004,
            f_strain_scale=0.5,
            j_tf_wp=23124470.793774806,
            t_tf_quench_dump=25.829000000000001,
            e_tf_coil_magnetic_stored=9548964780.4287167,
            temp_tf_coolant_peak_field=4.75,
            temp_tf_conductor_peak_quench=150,
            # These are picked to more or less correspond to the previous model
            # (but now we are bringing in magnetoresistivity)
            cu_rrr=33.0,
            fluence=0.0,
            detection_time=0.0,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.864553846654988,
            expected_j_tf_wp_quench_heat_max=17213147.288375787,
            expected_j_tf_wp_critical=49719296.722920775,
            expected_vd=9988.2637896807955,
            expected_j_superconductor_critical=832616175.5329928,
            expected_f_c_tf_turn_operating_critical=0.46510052068203006,
            expected_j_tf_coil_turn=23124470.793774802,
            expected_bc20m=32.97,
            expected_tc0m=16.06,
            expected_c_turn_cables_critical=159162.9081148869,
            expected_j_superconductor=387250216.7686755,
        ),
        SuperconParam(
            temp_tf_superconductor_margin_min=1.5,
            n_tf_coils=16,
            temp_margin=2.3431632224075836,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            c_tf_turn=74026.751437500003,
            b_tf_inboard_peak_with_ripple=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            f_a_tf_turn_cable_space_cooling=0.3,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            f_b_tf_inboard_peak_ripple_symmetric=1.0658869305062604,
            run_tests=0,
            i_tf_superconductor=5,
            iprint=0,
            outfile=11,
            a_tf_turn_cable_space=0.001293323051622732,
            a_tf_turn=0.0032012300777680192,
            a_tf_turn_cable_space_effective=0.001,
            b_tf_inboard_peak_symmetric=12.48976756562082,
            f_a_tf_turn_cable_copper=0.80884,
            f_a_tf_turn_cooling_extra=0.30000000000000004,
            f_strain_scale=0.5,
            j_tf_wp=23124470.793774806,
            t_tf_quench_dump=25.829000000000001,
            e_tf_coil_magnetic_stored=9561415368.8360519,
            temp_tf_coolant_peak_field=4.75,
            temp_tf_conductor_peak_quench=150,
            cu_rrr=33.0,
            fluence=0.0,
            detection_time=0.0,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.864553846654988,
            expected_j_tf_wp_quench_heat_max=17213147.288375787,
            expected_j_tf_wp_critical=49719296.722920775,
            expected_vd=10001.287165953383,
            expected_j_superconductor_critical=832616175.5329928,
            expected_f_c_tf_turn_operating_critical=0.46510052068203006,
            expected_j_tf_coil_turn=23124470.793774802,
            expected_bc20m=32.97,
            expected_tc0m=16.06,
            expected_c_turn_cables_critical=159162.9081148869,
            expected_j_superconductor=387250216.7686755,
        ),
    ),
)
def test_supercon(superconparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for supercon.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param superconparam: the data used to mock and assert in this test.
    :type superconparam: superconparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    monkeypatch.setattr(
        tfcoil_variables,
        "temp_tf_superconductor_margin_min",
        superconparam.temp_tf_superconductor_margin_min,
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", superconparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "temp_margin", superconparam.temp_margin)

    monkeypatch.setattr(tfcoil_variables, "rrr_tf_cu", superconparam.cu_rrr)

    monkeypatch.setattr(
        tfcoil_variables, "t_tf_quench_detection", superconparam.detection_time
    )

    monkeypatch.setattr(constraint_variables, "nflutfmax", superconparam.fluence)

    monkeypatch.setattr(
        tfcoil_variables,
        "dia_tf_turn_coolant_channel",
        superconparam.dia_tf_turn_coolant_channel,
    )

    monkeypatch.setattr(tfcoil_variables, "c_tf_turn", superconparam.c_tf_turn)

    monkeypatch.setattr(
        tfcoil_variables,
        "b_tf_inboard_peak_with_ripple",
        superconparam.b_tf_inboard_peak_with_ripple,
    )

    monkeypatch.setattr(tfcoil_variables, "str_tf_con_res", superconparam.str_tf_con_res)

    monkeypatch.setattr(
        tfcoil_variables, "b_crit_upper_nbti", superconparam.b_crit_upper_nbti
    )

    monkeypatch.setattr(tfcoil_variables, "i_str_wp", superconparam.i_str_wp)

    monkeypatch.setattr(tfcoil_variables, "str_wp", superconparam.str_wp)

    monkeypatch.setattr(tfcoil_variables, "t_crit_nbti", superconparam.t_crit_nbti)

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "tf_fit_t", superconparam.tf_fit_t
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "tf_fit_z", superconparam.tf_fit_z
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "f_b_tf_inboard_peak_ripple_symmetric",
        superconparam.f_b_tf_inboard_peak_ripple_symmetric,
    )

    monkeypatch.setattr(global_variables, "run_tests", superconparam.run_tests)

    (
        j_tf_wp_critical,
        j_superconductor_critical,
        f_c_tf_turn_operating_critical,
        j_superconductor,
        j_tf_coil_turn,
        bc20m,
        tc0m,
        c_turn_cables_critical,
    ) = sctfcoil.tf_cable_in_conduit_superconductor_properties(
        i_tf_superconductor=superconparam.i_tf_superconductor,
        a_tf_turn_cable_space=superconparam.a_tf_turn_cable_space,
        a_tf_turn=superconparam.a_tf_turn,
        a_tf_turn_cable_space_effective=superconparam.a_tf_turn_cable_space_effective,
        f_a_tf_turn_cable_space_cooling=superconparam.f_a_tf_turn_cable_space_cooling,
        b_tf_inboard_peak=superconparam.b_tf_inboard_peak_symmetric,
        f_a_tf_turn_cable_copper=superconparam.f_a_tf_turn_cable_copper,
        f_strain_scale=superconparam.f_strain_scale,
        c_tf_turn=superconparam.c_tf_turn,
        j_tf_wp=superconparam.j_tf_wp,
        temp_tf_coolant_peak_field=superconparam.temp_tf_coolant_peak_field,
        bcritsc=superconparam.bcritsc,
        tcritsc=superconparam.tcritsc,
    )

    assert j_superconductor == pytest.approx(superconparam.expected_j_superconductor)

    assert j_tf_wp_critical == pytest.approx(superconparam.expected_j_tf_wp_critical)

    assert j_superconductor_critical == pytest.approx(
        superconparam.expected_j_superconductor_critical
    )

    assert f_c_tf_turn_operating_critical == pytest.approx(
        superconparam.expected_f_c_tf_turn_operating_critical
    )

    assert j_tf_coil_turn == pytest.approx(superconparam.expected_j_tf_coil_turn)

    assert bc20m == pytest.approx(superconparam.expected_bc20m)

    assert tc0m == pytest.approx(superconparam.expected_tc0m)

    assert c_turn_cables_critical == pytest.approx(
        superconparam.expected_c_turn_cables_critical
    )


class PeakTfWithRippleParam(NamedTuple):
    tf_fit_t: Any = None

    tf_fit_z: Any = None

    f_b_tf_inboard_peak_ripple_symmetric: Any = None

    n_tf_coils: Any = None

    dx_tf_wp_primary_toroidal: Any = None

    dr_tf_wp_with_insulation: Any = None

    tfin: Any = None

    b_tf_inboard_peak_symmetric: Any = None

    expected_tf_fit_t: Any = None

    expected_tf_fit_z: Any = None

    expected_f_b_tf_inboard_peak_ripple_symmetric: Any = None

    expected_b_tf_inboard_peak_with_ripple: Any = None


@pytest.mark.parametrize(
    "peaktfwithrippleparam",
    (
        PeakTfWithRippleParam(
            tf_fit_t=0,
            tf_fit_z=0,
            f_b_tf_inboard_peak_ripple_symmetric=0,
            n_tf_coils=16,
            dx_tf_wp_primary_toroidal=1.299782604942499,
            dr_tf_wp_with_insulation=0.50661087836601015,
            tfin=3.789896624292115,
            b_tf_inboard_peak_symmetric=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_f_b_tf_inboard_peak_ripple_symmetric=1.0658869305062604,
            expected_b_tf_inboard_peak_with_ripple=12.48976756562082,
        ),
        PeakTfWithRippleParam(
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            f_b_tf_inboard_peak_ripple_symmetric=1.0658869305062604,
            n_tf_coils=16,
            dx_tf_wp_primary_toroidal=1.299782604942499,
            dr_tf_wp_with_insulation=0.50661087836601015,
            tfin=3.789896624292115,
            b_tf_inboard_peak_symmetric=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_f_b_tf_inboard_peak_ripple_symmetric=1.0658869305062604,
            expected_b_tf_inboard_peak_with_ripple=12.48976756562082,
        ),
    ),
)
def test_peak_tf_with_ripple(peaktfwithrippleparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for peak_tf_with_ripple.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param peaktfwithrippleparam: the data used to mock and assert in this test.
    :type peaktfwithrippleparam: peaktfwithrippleparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "tf_fit_t", peaktfwithrippleparam.tf_fit_t
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "tf_fit_z", peaktfwithrippleparam.tf_fit_z
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "f_b_tf_inboard_peak_ripple_symmetric",
        peaktfwithrippleparam.f_b_tf_inboard_peak_ripple_symmetric,
    )

    b_tf_inboard_peak_with_ripple = sctfcoil.peak_b_tf_inboard_with_ripple(
        n_tf_coils=peaktfwithrippleparam.n_tf_coils,
        dx_tf_wp_primary_toroidal=peaktfwithrippleparam.dx_tf_wp_primary_toroidal,
        dr_tf_wp_no_insulation=peaktfwithrippleparam.dr_tf_wp_with_insulation,
        r_tf_wp_inboard_centre=peaktfwithrippleparam.tfin,
        b_tf_inboard_peak_symmetric=peaktfwithrippleparam.b_tf_inboard_peak_symmetric,
    )

    assert superconducting_tf_coil_variables.tf_fit_t == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_t
    )

    assert superconducting_tf_coil_variables.tf_fit_z == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_z
    )

    assert (
        superconducting_tf_coil_variables.f_b_tf_inboard_peak_ripple_symmetric
        == pytest.approx(
            peaktfwithrippleparam.expected_f_b_tf_inboard_peak_ripple_symmetric
        )
    )

    assert b_tf_inboard_peak_with_ripple == pytest.approx(
        peaktfwithrippleparam.expected_b_tf_inboard_peak_with_ripple
    )


class TfWpGeomParam(NamedTuple):
    dr_tf_inboard: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    dr_tf_wp_with_insulation: Any = None

    dr_tf_wp_no_insulation: Any = None

    dr_tf_plasma_case: Any = None

    dr_tf_nose_case: Any = None

    dx_tf_side_case_min: Any = None

    dx_tf_wp_primary_toroidal: Any = None

    dx_tf_wp_secondary_toroidal: Any = None

    dx_tf_wp_insulation: Any = None

    dx_tf_wp_insertion_gap: Any = None

    a_tf_wp_with_insulation: Any = None

    a_tf_wp_no_insulation: Any = None

    r_tf_wp_inboard_inner: Any = None

    r_tf_wp_inboard_outer: Any = None

    r_tf_wp_inboard_centre: Any = None

    dx_tf_wp_toroidal_min: Any = None

    dx_tf_wp_toroidal_average: Any = None

    a_tf_wp_ground_insulation: Any = None

    rad_tf_coil_inboard_toroidal_half: Any = None

    tan_theta_coil: Any = None

    i_tf_wp_geom: Any = None

    expected_dx_tf_wp_primary_toroidal: Any = None

    expected_dx_tf_wp_secondary_toroidal: Any = None

    expected_a_tf_wp_with_insulation: Any = None

    expected_a_tf_wp_no_insulation: Any = None

    expected_dr_tf_wp_no_insulation: Any = None

    expected_r_tf_wp_inboard_inner: Any = None

    expected_r_tf_wp_inboard_outer: Any = None

    expected_r_tf_wp_inboard_centre: Any = None

    expected_t_wp_toroidal: Any = None

    expected_dx_tf_wp_toroidal_average: Any = None

    expected_a_tf_wp_ground_insulation: Any = None


@pytest.mark.parametrize(
    "tfwpgeomparam",
    (
        TfWpGeomParam(
            dr_tf_inboard=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp_with_insulation=0.54261087836601019,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_nose_case=0.52465000000000006,
            dx_tf_side_case_min=0.05000000000000001,
            dx_tf_wp_primary_toroidal=0,
            dx_tf_wp_secondary_toroidal=0,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            a_tf_wp_with_insulation=0,
            a_tf_wp_no_insulation=0,
            r_tf_wp_inboard_inner=0,
            r_tf_wp_inboard_outer=0,
            r_tf_wp_inboard_centre=0,
            dx_tf_wp_toroidal_min=0,
            dx_tf_wp_toroidal_average=0,
            a_tf_wp_ground_insulation=0,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            expected_dx_tf_wp_primary_toroidal=1.299782604942499,
            expected_dx_tf_wp_secondary_toroidal=1.299782604942499,
            expected_a_tf_wp_with_insulation=0.70527618095271016,
            expected_a_tf_wp_no_insulation=0.64024601555360383,
            expected_r_tf_wp_inboard_inner=3.5185911851091101,
            expected_dr_tf_wp_no_insulation=0.5066108783660102,
            expected_r_tf_wp_inboard_outer=4.06120206347512,
            expected_r_tf_wp_inboard_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_dx_tf_wp_toroidal_average=1.299782604942499,
            expected_a_tf_wp_ground_insulation=0.028582295732936136,
        ),
        TfWpGeomParam(
            dr_tf_inboard=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp_with_insulation=0.54261087836601019,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_nose_case=0.52465000000000006,
            dx_tf_side_case_min=0.05000000000000001,
            dx_tf_wp_primary_toroidal=1.299782604942499,
            dx_tf_wp_secondary_toroidal=0,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            a_tf_wp_with_insulation=0.70527618095271016,
            a_tf_wp_no_insulation=0.64024601555360383,
            r_tf_wp_inboard_inner=3.5185911851091101,
            r_tf_wp_inboard_outer=4.06120206347512,
            r_tf_wp_inboard_centre=3.789896624292115,
            dx_tf_wp_toroidal_min=1.299782604942499,
            dx_tf_wp_toroidal_average=1.299782604942499,
            a_tf_wp_ground_insulation=0.028582295732936136,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            expected_dx_tf_wp_primary_toroidal=1.299782604942499,
            expected_dx_tf_wp_secondary_toroidal=1.299782604942499,
            expected_a_tf_wp_with_insulation=0.70527618095271016,
            expected_a_tf_wp_no_insulation=0.64024601555360383,
            expected_dr_tf_wp_no_insulation=0.5066108783660102,
            expected_r_tf_wp_inboard_inner=3.5185911851091101,
            expected_r_tf_wp_inboard_outer=4.06120206347512,
            expected_r_tf_wp_inboard_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_dx_tf_wp_toroidal_average=1.299782604942499,
            expected_a_tf_wp_ground_insulation=0.028582295732936136,
        ),
        # New test case for i_tf_wp_geom=1
        TfWpGeomParam(
            dr_tf_inboard=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp_with_insulation=0.54261087836601019,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_nose_case=0.52465000000000006,
            dx_tf_side_case_min=0.05000000000000001,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=1,
            expected_dx_tf_wp_primary_toroidal=1.4077146193242378,
            expected_dx_tf_wp_secondary_toroidal=1.299782604942499,
            expected_a_tf_wp_with_insulation=0.7345587235164542,
            expected_a_tf_wp_no_insulation=0.6675857818584766,
            expected_dr_tf_wp_no_insulation=0.5066108783660102,
            expected_r_tf_wp_inboard_inner=3.5185911851091101,
            expected_r_tf_wp_inboard_outer=4.06120206347512,
            expected_r_tf_wp_inboard_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_dx_tf_wp_toroidal_average=1.3537486121333684,
            expected_a_tf_wp_ground_insulation=0.02944575184799003,
        ),
        # New test case for i_tf_wp_geom=2
        TfWpGeomParam(
            dr_tf_inboard=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp_with_insulation=0.54261087836601019,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_nose_case=0.52465000000000006,
            dx_tf_side_case_min=0.05000000000000001,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=2,
            expected_dx_tf_wp_primary_toroidal=1.5156466337059764,
            expected_dx_tf_wp_secondary_toroidal=1.299782604942499,
            expected_a_tf_wp_with_insulation=0.7638412660801982,
            expected_a_tf_wp_no_insulation=0.6949255481633493,
            expected_dr_tf_wp_no_insulation=0.5066108783660102,
            expected_r_tf_wp_inboard_inner=3.5185911851091101,
            expected_r_tf_wp_inboard_outer=4.06120206347512,
            expected_r_tf_wp_inboard_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_dx_tf_wp_toroidal_average=1.4077146193242376,
            expected_a_tf_wp_ground_insulation=0.030309207963043927,
        ),
    ),
)
def test_superconducting_tf_wp_geometry(tfwpgeomparam, sctfcoil):
    """
    Automatically generated Regression Unit Test for superconducting_tf_wp_geometry.

    This test was generated using data from
    tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfwpgeomparam: the data used to mock and assert in this test.
    :type tfwpgeomparam: tfwpgeomparam
    """

    (
        r_tf_wp_inboard_inner,
        r_tf_wp_inboard_outer,
        r_tf_wp_inboard_centre,
        dx_tf_wp_toroidal_min,
        dr_tf_wp_no_insulation,
        dx_tf_wp_primary_toroidal,
        dx_tf_wp_secondary_toroidal,
        dx_tf_wp_toroidal_average,
        a_tf_wp_with_insulation,
        a_tf_wp_no_insulation,
        a_tf_wp_ground_insulation,
    ) = sctfcoil.superconducting_tf_wp_geometry(
        i_tf_wp_geom=tfwpgeomparam.i_tf_wp_geom,
        r_tf_inboard_in=tfwpgeomparam.r_tf_inboard_in,
        dr_tf_nose_case=tfwpgeomparam.dr_tf_nose_case,
        dr_tf_wp_with_insulation=tfwpgeomparam.dr_tf_wp_with_insulation,
        tan_theta_coil=tfwpgeomparam.tan_theta_coil,
        dx_tf_side_case_min=tfwpgeomparam.dx_tf_side_case_min,
        dx_tf_wp_insulation=tfwpgeomparam.dx_tf_wp_insulation,
        dx_tf_wp_insertion_gap=tfwpgeomparam.dx_tf_wp_insertion_gap,
    )

    assert dx_tf_wp_primary_toroidal == pytest.approx(
        tfwpgeomparam.expected_dx_tf_wp_primary_toroidal
    )

    assert dx_tf_wp_secondary_toroidal == pytest.approx(
        tfwpgeomparam.expected_dx_tf_wp_secondary_toroidal
    )

    assert a_tf_wp_with_insulation == pytest.approx(
        tfwpgeomparam.expected_a_tf_wp_with_insulation
    )

    assert a_tf_wp_no_insulation == pytest.approx(
        tfwpgeomparam.expected_a_tf_wp_no_insulation
    )

    assert dr_tf_wp_no_insulation == pytest.approx(
        tfwpgeomparam.expected_dr_tf_wp_no_insulation
    )

    assert r_tf_wp_inboard_inner == pytest.approx(
        tfwpgeomparam.expected_r_tf_wp_inboard_inner
    )

    assert r_tf_wp_inboard_outer == pytest.approx(
        tfwpgeomparam.expected_r_tf_wp_inboard_outer
    )

    assert r_tf_wp_inboard_centre == pytest.approx(
        tfwpgeomparam.expected_r_tf_wp_inboard_centre
    )

    assert dx_tf_wp_toroidal_min == pytest.approx(tfwpgeomparam.expected_t_wp_toroidal)

    assert dx_tf_wp_toroidal_average == pytest.approx(
        tfwpgeomparam.expected_dx_tf_wp_toroidal_average
    )

    assert a_tf_wp_ground_insulation == pytest.approx(
        tfwpgeomparam.expected_a_tf_wp_ground_insulation
    )


class TfCaseGeomParam(NamedTuple):
    a_tf_leg_outboard: Any = None

    a_tf_inboard_total: Any = None

    n_tf_coils: Any = None

    dx_tf_side_case_min: Any = None

    dr_tf_plasma_case: Any = None

    dr_tf_wp_with_insulation: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    a_tf_wp_with_insulation: Any = None

    r_tf_wp_inboard_inner: Any = None

    r_tf_wp_inboard_outer: Any = None

    rad_tf_coil_inboard_toroidal_half: Any = None

    tan_theta_coil: Any = None

    i_tf_wp_geom: Any = None

    i_tf_case_geom: Any = None

    expected_a_tf_coil_inboard_case: Any = None

    expected_a_tf_coil_outboard_case: Any = None

    expected_dx_tf_side_case_average: Any = None

    expected_dx_tf_side_case_peak: Any = None

    expected_a_tf_plasma_case: Any = None

    expected_a_tf_coil_nose_case: Any = None


@pytest.mark.parametrize(
    "tfcasegeomparam",
    (
        TfCaseGeomParam(
            a_tf_leg_outboard=1.9805354702921749,
            a_tf_inboard_total=27.308689677971632,
            n_tf_coils=16,
            dx_tf_side_case_min=0.05000000000000001,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_wp_with_insulation=0.54261087836601019,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            a_tf_wp_with_insulation=0.70527618095271016,
            r_tf_wp_inboard_inner=3.5185911851091101,
            r_tf_wp_inboard_outer=4.06120206347512,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            i_tf_case_geom=0,
            expected_a_tf_coil_inboard_case=1.0015169239205168,
            expected_a_tf_coil_outboard_case=1.2752592893394648,
            expected_dx_tf_side_case_average=0.10396600719086938,
            expected_dx_tf_side_case_peak=0.15793201438173876,
            expected_a_tf_plasma_case=0.18607458590131154,
            expected_a_tf_coil_nose_case=0.70261616505511615,
        ),
        TfCaseGeomParam(
            a_tf_leg_outboard=1.9805354702921749,
            a_tf_inboard_total=27.308689677971632,
            n_tf_coils=16,
            dx_tf_side_case_min=0.05000000000000001,
            dr_tf_plasma_case=0.060000000000000012,
            dr_tf_wp_with_insulation=0.54261087836601019,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            a_tf_wp_with_insulation=0.70527618095271016,
            r_tf_wp_inboard_inner=3.5185911851091101,
            r_tf_wp_inboard_outer=4.06120206347512,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            i_tf_case_geom=0,
            expected_a_tf_coil_inboard_case=1.0015169239205168,
            expected_a_tf_coil_outboard_case=1.2752592893394648,
            expected_dx_tf_side_case_average=0.10396600719086938,
            expected_dx_tf_side_case_peak=0.15793201438173876,
            expected_a_tf_plasma_case=0.18607458590131154,
            expected_a_tf_coil_nose_case=0.70261616505511615,
        ),
    ),
)
def test_superconducting_tf_case_geometry(tfcasegeomparam, sctfcoil):
    """
    Automatically generated Regression Unit Test for superconducting_tf_case_geometry.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcasegeomparam: the data used to mock and assert in this test.
    :type tfcasegeomparam: tfcasegeomparam
    """
    (
        a_tf_coil_inboard_case,
        a_tf_coil_outboard_case,
        a_tf_plasma_case,
        a_tf_coil_nose_case,
        dx_tf_side_case_average,
        dx_tf_side_case_peak,
    ) = sctfcoil.superconducting_tf_case_geometry(
        i_tf_wp_geom=tfcasegeomparam.i_tf_wp_geom,
        i_tf_case_geom=tfcasegeomparam.i_tf_case_geom,
        a_tf_inboard_total=tfcasegeomparam.a_tf_inboard_total,
        n_tf_coils=tfcasegeomparam.n_tf_coils,
        a_tf_wp_with_insulation=tfcasegeomparam.a_tf_wp_with_insulation,
        a_tf_leg_outboard=tfcasegeomparam.a_tf_leg_outboard,
        rad_tf_coil_inboard_toroidal_half=tfcasegeomparam.rad_tf_coil_inboard_toroidal_half,
        r_tf_inboard_out=tfcasegeomparam.r_tf_inboard_out,
        tan_theta_coil=tfcasegeomparam.tan_theta_coil,
        r_tf_wp_inboard_outer=tfcasegeomparam.r_tf_wp_inboard_outer,
        dr_tf_plasma_case=tfcasegeomparam.dr_tf_plasma_case,
        r_tf_wp_inboard_inner=tfcasegeomparam.r_tf_wp_inboard_inner,
        r_tf_inboard_in=tfcasegeomparam.r_tf_inboard_in,
        dx_tf_side_case_min=tfcasegeomparam.dx_tf_side_case_min,
        dr_tf_wp_with_insulation=tfcasegeomparam.dr_tf_wp_with_insulation,
    )

    assert a_tf_coil_inboard_case == pytest.approx(
        tfcasegeomparam.expected_a_tf_coil_inboard_case
    )

    assert a_tf_coil_outboard_case == pytest.approx(
        tfcasegeomparam.expected_a_tf_coil_outboard_case
    )

    assert dx_tf_side_case_average == pytest.approx(
        tfcasegeomparam.expected_dx_tf_side_case_average
    )

    assert dx_tf_side_case_peak == pytest.approx(
        tfcasegeomparam.expected_dx_tf_side_case_peak
    )

    assert a_tf_plasma_case == pytest.approx(tfcasegeomparam.expected_a_tf_plasma_case)

    assert a_tf_coil_nose_case == pytest.approx(
        tfcasegeomparam.expected_a_tf_coil_nose_case
    )


class TfIntegerTurnGeomParam(NamedTuple):
    dr_tf_wp_with_insulation: Any = None

    dx_tf_wp_insulation: Any = None

    dx_tf_wp_insertion_gap: Any = None

    t_conductor: Any = None

    dx_tf_turn_general: Any = None

    c_tf_coil: Any = None

    dx_tf_wp_toroidal_min: Any = None

    t_conductor_radial: Any = None

    t_conductor_toroidal: Any = None

    dr_tf_turn_cable_space: Any = None

    dx_tf_turn_cable_space: Any = None

    dr_tf_turn: Any = None

    dx_tf_turn: Any = None

    dx_tf_turn_cable_space_average: Any = None

    n_tf_wp_layers: Any = None

    n_tf_wp_pancakes: Any = None

    dx_tf_turn_steel: Any = None

    dx_tf_turn_insulation: Any = None

    expected_t_conductor: Any = None

    expected_dx_tf_turn_general: Any = None

    expected_t_conductor_radial: Any = None

    expected_t_conductor_toroidal: Any = None

    expected_dr_tf_turn_cable_space: Any = None

    expected_dx_tf_turn_cable_space: Any = None

    expected_t_turn_radial: Any = None

    expected_dx_tf_turn: Any = None

    expected_t_cable: Any = None

    expected_a_tf_turn_cable_space: Any = None

    expected_a_tf_turn_steel: Any = None

    expected_a_tf_turn_insulation: Any = None

    expected_cpttf: Any = None

    expected_n_tf_coil_turns: Any = None


@pytest.mark.parametrize(
    "tfintegerturngeomparam",
    (
        TfIntegerTurnGeomParam(
            dr_tf_wp_with_insulation=0.54261087836601019,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            t_conductor=0,
            dx_tf_turn_general=0,
            c_tf_coil=14805350.287500001,
            dx_tf_wp_toroidal_min=1.299782604942499,
            t_conductor_radial=0,
            t_conductor_toroidal=0,
            dr_tf_turn_cable_space=0,
            dx_tf_turn_cable_space=0,
            dr_tf_turn=0,
            dx_tf_turn=0,
            dx_tf_turn_cable_space_average=0,
            n_tf_wp_layers=10,
            n_tf_wp_pancakes=20,
            dx_tf_turn_steel=0.0080000000000000002,
            dx_tf_turn_insulation=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_dx_tf_turn_general=0.056579413904423038,
            expected_t_conductor_radial=0.046661087836601015,
            expected_t_conductor_toroidal=0.059189130247124938,
            expected_dr_tf_turn_cable_space=0.030661087836601014,
            expected_dx_tf_turn_cable_space=0.043189130247124938,
            expected_t_turn_radial=0.050661087836601018,
            expected_dx_tf_turn=0.063189130247124942,
            expected_t_cable=0.036389912284773368,
            expected_a_tf_turn_cable_space=0.001293323051622732,
            expected_a_tf_turn_steel=0.0014685061538103825,
            expected_a_tf_turn_insulation=0.00043940087233490435,
            expected_cpttf=74026.751437500003,
            expected_n_tf_coil_turns=200,
        ),
        TfIntegerTurnGeomParam(
            dr_tf_wp_with_insulation=0.54261087836601019,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_wp_insertion_gap=0.01,
            t_conductor=0.052553108427885735,
            dx_tf_turn_general=0.056579413904423038,
            c_tf_coil=14805350.287500001,
            dx_tf_wp_toroidal_min=1.299782604942499,
            t_conductor_radial=0.046661087836601015,
            t_conductor_toroidal=0.059189130247124938,
            dr_tf_turn_cable_space=0.030661087836601014,
            dx_tf_turn_cable_space=0.043189130247124938,
            dr_tf_turn=0.050661087836601018,
            dx_tf_turn=0.063189130247124942,
            dx_tf_turn_cable_space_average=0.036389912284773368,
            n_tf_wp_layers=10,
            n_tf_wp_pancakes=20,
            dx_tf_turn_steel=0.0080000000000000002,
            dx_tf_turn_insulation=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_dx_tf_turn_general=0.056579413904423038,
            expected_t_conductor_radial=0.046661087836601015,
            expected_t_conductor_toroidal=0.059189130247124938,
            expected_dr_tf_turn_cable_space=0.030661087836601014,
            expected_dx_tf_turn_cable_space=0.043189130247124938,
            expected_t_turn_radial=0.050661087836601018,
            expected_dx_tf_turn=0.063189130247124942,
            expected_t_cable=0.036389912284773368,
            expected_a_tf_turn_cable_space=0.001293323051622732,
            expected_a_tf_turn_steel=0.0014685061538103825,
            expected_a_tf_turn_insulation=0.00043940087233490435,
            expected_cpttf=74026.751437500003,
            expected_n_tf_coil_turns=200,
        ),
    ),
)
def test_tf_cable_in_conduit_integer_turn_geometry(tfintegerturngeomparam, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_cable_in_conduit_integer_turn_geometry.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfintegerturngeomparam: the data used to mock and assert in this test.
    :type tfintegerturngeomparam: tfintegerturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        radius_tf_turn_cable_space_corners,
        dr_tf_turn,
        dx_tf_turn,
        a_tf_turn_cable_space_no_void,
        a_tf_turn_steel,
        a_tf_turn_insulation,
        c_tf_turn,
        n_tf_coil_turns,
        t_conductor_radial,
        t_conductor_toroidal,
        t_conductor,
        dr_tf_turn_cable_space,
        dx_tf_turn_cable_space,
        dx_tf_turn_cable_space_average,
    ) = sctfcoil.tf_cable_in_conduit_integer_turn_geometry(
        dr_tf_wp_with_insulation=tfintegerturngeomparam.dr_tf_wp_with_insulation,
        dx_tf_wp_insulation=tfintegerturngeomparam.dx_tf_wp_insulation,
        dx_tf_wp_insertion_gap=tfintegerturngeomparam.dx_tf_wp_insertion_gap,
        n_tf_wp_layers=tfintegerturngeomparam.n_tf_wp_layers,
        dx_tf_wp_toroidal_min=tfintegerturngeomparam.dx_tf_wp_toroidal_min,
        n_tf_wp_pancakes=tfintegerturngeomparam.n_tf_wp_pancakes,
        c_tf_coil=tfintegerturngeomparam.c_tf_coil,
        dx_tf_turn_steel=tfintegerturngeomparam.dx_tf_turn_steel,
        dx_tf_turn_insulation=tfintegerturngeomparam.dx_tf_turn_insulation,
    )

    assert radius_tf_turn_cable_space_corners == pytest.approx(
        0.75 * tfintegerturngeomparam.dx_tf_turn_steel
    )

    assert t_conductor == pytest.approx(tfintegerturngeomparam.expected_t_conductor)

    assert tfcoil_variables.dx_tf_turn_general == pytest.approx(
        tfintegerturngeomparam.expected_dx_tf_turn_general
    )

    assert t_conductor_radial == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_radial
    )

    assert t_conductor_toroidal == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_toroidal
    )

    assert dr_tf_turn_cable_space == pytest.approx(
        tfintegerturngeomparam.expected_dr_tf_turn_cable_space
    )

    assert dx_tf_turn_cable_space == pytest.approx(
        tfintegerturngeomparam.expected_dx_tf_turn_cable_space
    )

    assert dr_tf_turn == pytest.approx(tfintegerturngeomparam.expected_t_turn_radial)

    assert dx_tf_turn == pytest.approx(tfintegerturngeomparam.expected_dx_tf_turn)

    assert dx_tf_turn_cable_space_average == pytest.approx(
        tfintegerturngeomparam.expected_t_cable
    )

    assert a_tf_turn_cable_space_no_void == pytest.approx(
        tfintegerturngeomparam.expected_a_tf_turn_cable_space
    )

    assert a_tf_turn_steel == pytest.approx(
        tfintegerturngeomparam.expected_a_tf_turn_steel
    )

    assert a_tf_turn_insulation == pytest.approx(
        tfintegerturngeomparam.expected_a_tf_turn_insulation
    )

    assert c_tf_turn == pytest.approx(tfintegerturngeomparam.expected_cpttf)

    assert n_tf_coil_turns == pytest.approx(
        tfintegerturngeomparam.expected_n_tf_coil_turns
    )


class TfAveragedTurnGeomParam(NamedTuple):
    layer_ins: Any = None

    t_conductor: Any = None

    dx_tf_turn_general: Any = None

    i_dx_tf_turn_general_input: Any = None

    c_tf_turn: Any = None

    dx_tf_turn_cable_space_general: Any = None

    i_dx_tf_turn_cable_space_general_input: Any = None

    a_tf_wp_no_insulation: Any = None

    dr_tf_turn: Any = None

    dx_tf_turn: Any = None

    dx_tf_turn_cable_space_average: Any = None

    i_tf_sc_mat: Any = None

    j_tf_wp: Any = None

    dx_tf_turn_steel: Any = None

    dx_tf_turn_insulation: Any = None

    expected_t_conductor: Any = None

    expected_dx_tf_turn_general: Any = None

    expected_t_turn_radial: Any = None

    expected_dx_tf_turn: Any = None

    expected_t_cable: Any = None

    expected_a_tf_turn_cable_space: Any = None

    expected_a_tf_turn_steel: Any = None

    expected_a_tf_turn_insulation: Any = None

    expected_n_tf_coil_turns: Any = None


@pytest.mark.parametrize(
    "tfaveragedturngeomparam",
    (
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0,
            dx_tf_turn_general=0,
            i_dx_tf_turn_general_input=False,
            c_tf_turn=65000,
            dx_tf_turn_cable_space_general=0,
            i_dx_tf_turn_cable_space_general_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0,
            dx_tf_turn=0,
            dx_tf_turn_cable_space_average=0,
            i_tf_sc_mat=5,
            j_tf_wp=26493137.688284047,
            dx_tf_turn_steel=0.0080000000000000019,
            dx_tf_turn_insulation=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_dx_tf_turn_general=0.049532469413859428,
            expected_t_turn_radial=0.049532469413859428,
            expected_dx_tf_turn=0.049532469413859428,
            expected_t_cable=0.031932469413859424,
            expected_a_tf_turn_cable_space=0.00098877993839630008,
            expected_a_tf_turn_steel=0.0013087416857142699,
            expected_a_tf_turn_insulation=0.00015594390212434958,
            expected_n_tf_coil_turns=246.63461538461544,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0.047932469413859431,
            dx_tf_turn_general=0.049532469413859428,
            i_dx_tf_turn_general_input=False,
            c_tf_turn=65000,
            dx_tf_turn_cable_space_general=0,
            i_dx_tf_turn_cable_space_general_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.049532469413859428,
            dx_tf_turn=0.049532469413859428,
            dx_tf_turn_cable_space_average=0.031932469413859424,
            i_tf_sc_mat=5,
            j_tf_wp=26493137.688284047,
            dx_tf_turn_steel=0.0080000000000000019,
            dx_tf_turn_insulation=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_dx_tf_turn_general=0.049532469413859428,
            expected_t_turn_radial=0.049532469413859428,
            expected_dx_tf_turn=0.049532469413859428,
            expected_t_cable=0.031932469413859424,
            expected_a_tf_turn_cable_space=0.00098877993839630008,
            expected_a_tf_turn_steel=0.0013087416857142699,
            expected_a_tf_turn_insulation=0.00015594390212434958,
            expected_n_tf_coil_turns=246.63461538461544,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=5.712e-02,
            dx_tf_turn_general=0.05872,
            i_dx_tf_turn_general_input=True,
            c_tf_turn=0,
            dx_tf_turn_cable_space_general=0,
            i_dx_tf_turn_cable_space_general_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.05872,
            dx_tf_turn=0.05872,
            dx_tf_turn_cable_space_average=0.04109,
            i_tf_sc_mat=1,
            j_tf_wp=2.301e07,
            dx_tf_turn_steel=8.015e-03,
            dx_tf_turn_insulation=8.0e-4,
            expected_t_conductor=5.712e-02,
            expected_dx_tf_turn_general=0.05872,
            expected_t_turn_radial=0.05872,
            expected_dx_tf_turn=0.05872,
            expected_t_cable=0.04109,
            expected_a_tf_turn_cable_space=0.001657369442,
            expected_a_tf_turn_steel=0.001605324958,
            expected_a_tf_turn_insulation=0.000185344,
            expected_n_tf_coil_turns=175.49384787,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0.058296,
            dx_tf_turn_general=0,
            i_dx_tf_turn_general_input=False,
            c_tf_turn=0,
            dx_tf_turn_cable_space_general=0.042,
            i_dx_tf_turn_cable_space_general_input=True,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.05872,
            dx_tf_turn=0.05872,
            dx_tf_turn_cable_space_average=0.04109,
            i_tf_sc_mat=1,
            j_tf_wp=2.673e07,
            dx_tf_turn_steel=8.148e-03,
            dx_tf_turn_insulation=8.0e-4,
            expected_t_conductor=0.058296,
            expected_dx_tf_turn_general=0.059896,
            expected_t_turn_radial=0.059896,
            expected_dx_tf_turn=0.059896,
            expected_t_cable=0.042,
            expected_a_tf_turn_cable_space=0.001731943361,
            expected_a_tf_turn_steel=0.001666480255,
            expected_a_tf_turn_insulation=0.00018910719999999962,
            expected_n_tf_coil_turns=168.6701961481806,
        ),
    ),
)
def test_tf_cable_in_conduit_averaged_turn_geometry(
    tfaveragedturngeomparam, monkeypatch, sctfcoil
):
    """
    Automatically generated Regression Unit Test for tf_cable_in_conduit_averaged_turn_geometry.

    This test was generated using data from tests/regression/scenarios/i_mode/IN.DAT.

    :param tfaveragedturngeomparam: the data used to mock and assert in this test.
    :type tfaveragedturngeomparam: tfaveragedturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        a_tf_turn_cable_space_no_void,
        a_tf_turn_steel,
        a_tf_turn_insulation,
        n_tf_coil_turns,
        dx_tf_turn_general,
        c_tf_turn,
        dx_tf_turn_general2,
        dr_tf_turn,
        dx_tf_turn,
        t_conductor,
        radius_tf_turn_cable_space_corners,
        dx_tf_turn_cable_space_average,
        _a_tf_turn_cable_space_effective,
        f_a_tf_turn_cable_space_cooling,
    ) = sctfcoil.tf_cable_in_conduit_averaged_turn_geometry(
        j_tf_wp=tfaveragedturngeomparam.j_tf_wp,
        dx_tf_turn_steel=tfaveragedturngeomparam.dx_tf_turn_steel,
        dx_tf_turn_insulation=tfaveragedturngeomparam.dx_tf_turn_insulation,
        i_tf_sc_mat=tfaveragedturngeomparam.i_tf_sc_mat,
        dx_tf_turn_general=tfaveragedturngeomparam.dx_tf_turn_general,
        c_tf_turn=tfaveragedturngeomparam.c_tf_turn,
        i_dx_tf_turn_general_input=tfaveragedturngeomparam.i_dx_tf_turn_general_input,
        i_dx_tf_turn_cable_space_general_input=tfaveragedturngeomparam.i_dx_tf_turn_cable_space_general_input,
        dx_tf_turn_cable_space_general=tfaveragedturngeomparam.dx_tf_turn_cable_space_general,
        layer_ins=tfaveragedturngeomparam.layer_ins,
        a_tf_wp_no_insulation=tfaveragedturngeomparam.a_tf_wp_no_insulation,
        dia_tf_turn_coolant_channel=0.004,
        f_a_tf_turn_cable_space_extra_void=0.3,
    )

    # Existing checks
    assert t_conductor == pytest.approx(tfaveragedturngeomparam.expected_t_conductor)
    assert dx_tf_turn_general == pytest.approx(
        tfaveragedturngeomparam.expected_dx_tf_turn_general
    )
    assert dr_tf_turn == pytest.approx(tfaveragedturngeomparam.expected_t_turn_radial)
    assert dx_tf_turn == pytest.approx(tfaveragedturngeomparam.expected_dx_tf_turn)
    assert dx_tf_turn_cable_space_average == pytest.approx(
        tfaveragedturngeomparam.expected_t_cable
    )
    assert a_tf_turn_cable_space_no_void == pytest.approx(
        tfaveragedturngeomparam.expected_a_tf_turn_cable_space
    )
    assert a_tf_turn_steel == pytest.approx(
        tfaveragedturngeomparam.expected_a_tf_turn_steel
    )
    assert a_tf_turn_insulation == pytest.approx(
        tfaveragedturngeomparam.expected_a_tf_turn_insulation
    )
    assert n_tf_coil_turns == pytest.approx(
        tfaveragedturngeomparam.expected_n_tf_coil_turns
    )

    # Expanded checks for unchecked variables
    assert radius_tf_turn_cable_space_corners == pytest.approx(
        0.75 * tfaveragedturngeomparam.dx_tf_turn_steel
    )
    # c_tf_turn is an input, so just check it matches input if input is used
    if (
        not tfaveragedturngeomparam.i_dx_tf_turn_general_input
        and tfaveragedturngeomparam.c_tf_turn != 0
    ):
        assert c_tf_turn == pytest.approx(tfaveragedturngeomparam.c_tf_turn)
    # dx_tf_turn_general2 should match dx_tf_turn_general
    assert dx_tf_turn_general2 == pytest.approx(dx_tf_turn_general)
    # f_a_tf_turn_cable_space_cooling should be a float between 0 and 1
    assert 0.0 <= f_a_tf_turn_cable_space_cooling <= 1.0


class TfWpCurrentsParam(NamedTuple):
    c_tf_total: Any = None

    n_tf_coils: Any = None

    j_tf_wp: Any = None

    a_tf_wp_no_insulation: Any = None

    expected_j_tf_wp: Any = None


@pytest.mark.parametrize(
    "tfwpcurrentsparam",
    (
        TfWpCurrentsParam(
            c_tf_total=256500000.00000003,
            n_tf_coils=16,
            j_tf_wp=0,
            a_tf_wp_no_insulation=0.60510952642236249,
            expected_j_tf_wp=26493137.688284047,
        ),
        TfWpCurrentsParam(
            c_tf_total=256500000.00000003,
            n_tf_coils=16,
            j_tf_wp=26493137.688284047,
            a_tf_wp_no_insulation=0.60510952642236249,
            expected_j_tf_wp=26493137.688284047,
        ),
    ),
)
def test_tf_wp_currents(tfwpcurrentsparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_wp_currents.

    This test was generated using data from tests/regression/scenarios/i_mode/IN.DAT.

    :param tfwpcurrentsparam: the data used to mock and assert in this test.
    :type tfwpcurrentsparam: tfwpcurrentsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "c_tf_total", tfwpcurrentsparam.c_tf_total)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", tfwpcurrentsparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "j_tf_wp", tfwpcurrentsparam.j_tf_wp)

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "a_tf_wp_no_insulation",
        tfwpcurrentsparam.a_tf_wp_no_insulation,
    )

    sctfcoil.tf_wp_currents()

    assert tfcoil_variables.j_tf_wp == pytest.approx(tfwpcurrentsparam.expected_j_tf_wp)


def test_vv_stress_on_quench():
    """Tests the VV stress on TF quench model presented in Itoh et al using the
    values they use to test the model for JA DEMO concept.
    """
    assert (
        pytest.approx(
            sctf.vv_stress_on_quench(
                # TF shape
                H_coil=9.5,
                ri_coil=3.55,
                ro_coil=15.62,
                rm_coil=7.66,
                ccl_length_coil=51.1,
                theta1_coil=48,
                # VV shape
                H_vv=7.9,
                ri_vv=4.45,
                ro_vv=13.09,
                rm_vv=7.88,
                theta1_vv=1,
                # TF properties
                n_tf_coils=18,
                n_tf_coil_turns=192,
                s_rp=0.55,
                s_cc=0.94,
                taud=30,
                i_op=83200,
                # VV properties
                d_vv=0.12,  # for 6.6 restistance -> lambda2 = 2.1
            )
        )
        == 57045874.69917925
    )


def test_vv_stress_on_quench_integration(sctfcoil, monkeypatch):
    """Tests the VV stress on TF quench model presented in Itoh et al using the
    values they use to test the model for JA DEMO concept. Includes the assumptions
    and approximations in the models integration with PROCESS.
    """
    monkeypatch.setattr(build_variables, "dr_tf_inboard", 1.4)  # Baseline 2018 value
    monkeypatch.setattr(build_variables, "z_tf_inside_half", 8.8)  # Table 2
    monkeypatch.setattr(build_variables, "r_tf_inboard_mid", 3.55)  # Table 2
    monkeypatch.setattr(build_variables, "r_tf_outboard_mid", 15.62)  # Table 2
    monkeypatch.setattr(tfcoil_variables, "theta1_coil", 48)  # Table 2
    monkeypatch.setattr(tfcoil_variables, "theta1_vv", 1)  # Table 2
    monkeypatch.setattr(
        build_variables,
        "r_tf_inboard_out",
        build_variables.r_tf_inboard_mid + (build_variables.dr_tf_inboard / 2),
    )
    monkeypatch.setattr(build_variables, "dr_tf_outboard", 0)  # simplifies

    monkeypatch.setattr(
        build_variables, "z_plasma_xpoint_upper", 5.47008
    )  # Baseline 2018

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "a_tf_coil_inboard_steel", 0.55
    )  # Section 3

    # Sum from Section 3
    monkeypatch.setattr(superconducting_tf_coil_variables, "a_tf_plasma_case", 0.42)
    monkeypatch.setattr(superconducting_tf_coil_variables, "a_tf_coil_nose_case", 0.42)
    monkeypatch.setattr(
        superconducting_tf_coil_variables, "dx_tf_side_case_average", 0.05
    )

    monkeypatch.setattr(build_variables, "dz_xpoint_divertor", 0.05)  # Baseline 2018
    monkeypatch.setattr(build_variables, "dz_shld_upper", 0.3)  # Baseline 2018
    monkeypatch.setattr(
        divertor_variables, "dz_divertor", 2.05
    )  # chosen to achieve H_vv in Table 2

    monkeypatch.setattr(build_variables, "dr_tf_shld_gap", 0.05)  # Baseline 2018
    monkeypatch.setattr(
        build_variables, "dr_shld_thermal_outboard", 0.05
    )  # Baseline 2018
    monkeypatch.setattr(build_variables, "dr_tf_outboard", 1.4)  # Baseline 2018
    monkeypatch.setattr(
        build_variables, "dr_shld_vv_gap_outboard", 1.7
    )  # chosen to achieve Ro_vv in Table 2

    monkeypatch.setattr(build_variables, "dr_vv_outboard", 0.06)  # Section 3
    monkeypatch.setattr(build_variables, "dr_vv_inboard", 0.06)  # Section 3
    monkeypatch.setattr(build_variables, "dz_vv_upper", 0.06)  # Section 3

    monkeypatch.setattr(tfcoil_variables, "len_tf_coil", 51.1)  # Table 2
    monkeypatch.setattr(
        tfcoil_variables, "tfa", [3.41, 7.77, 7.77, 3.41]
    )  # chosen to achieve Rm_coil in Table 2
    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", 18)  # Section 3
    monkeypatch.setattr(tfcoil_variables, "n_tf_coil_turns", 192)  # Section 3
    monkeypatch.setattr(tfcoil_variables, "t_tf_superconductor_quench", 30)  # Figure 6
    monkeypatch.setattr(
        superconducting_tf_coil_variables, "c_tf_coil", 83200 * 192
    )  # Section 3

    monkeypatch.setattr(
        build_variables, "r_vv_inboard_out", 4.45 + (build_variables.dr_vv_inboard / 2)
    )  # Table 2

    sctfcoil.vv_stress_on_quench()

    assert (
        pytest.approx(superconducting_tf_coil_variables.vv_stress_quench)
        == 56893800.120420754
    )


class TfCoilAreaAndMassesParam(NamedTuple):
    hr1: Any = None

    r_tf_outboard_mid: Any = None

    dr_tf_inboard: Any = None

    r_tf_inboard_mid: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    z_tf_inside_half: Any = None

    den_steel: Any = None

    m_tf_wp_steel_conduit: Any = None

    m_tf_coils_total: Any = None

    m_tf_coil_case: Any = None

    tficrn: Any = None

    tfcryoarea: Any = None

    m_tf_coil_wp_insulation: Any = None

    tfocrn: Any = None

    m_tf_coil_superconductor: Any = None

    m_tf_coil_copper: Any = None

    m_tf_coil_conductor: Any = None

    m_tf_coil_wp_turn_insulation: Any = None

    f_a_tf_turn_cable_space_extra_void: Any = None

    dcond: Any = None

    den_tf_wp_turn_insulation: Any = None

    len_tf_coil: Any = None

    den_tf_coil_case: Any = None

    a_tf_turn_steel: Any = None

    n_tf_coil_turns: Any = None

    n_tf_coils: Any = None

    a_tf_coil_wp_turn_insulation: Any = None

    a_tf_coil_outboard_case: Any = None

    a_tf_coil_inboard_case: Any = None

    f_a_tf_turn_cable_copper: Any = None

    a_tf_wp_coolant_channels: Any = None

    a_tf_turn_cable_space_no_void: Any = None

    whttflgs: Any = None

    whtcp: Any = None

    whtconal: Any = None

    vol_cond_cp: Any = None

    i_tf_sup: Any = None

    i_tf_sc_mat: Any = None

    a_tf_leg_outboard: Any = None

    dr_tf_nose_case: Any = None

    voltfleg: Any = None

    cplen: Any = None

    itart: Any = None

    a_tf_wp_with_insulation: Any = None

    a_tf_wp_no_insulation: Any = None

    vol_ins_cp: Any = None

    vol_gr_ins_cp: Any = None

    vol_case_cp: Any = None

    a_leg_ins: Any = None

    a_leg_gr_ins: Any = None

    a_leg_cond: Any = None

    rad_tf_coil_inboard_toroidal_half: Any = None

    tan_theta_coil: Any = None

    expected_m_tf_wp_steel_conduit: Any = None

    expected_m_tf_coil_case: Any = None

    expected_tficrn: Any = None

    expected_tfcryoarea: Any = None

    expected_m_tf_coil_wp_insulation: Any = None

    expected_tfocrn: Any = None

    expected_m_tf_coil_superconductor: Any = None

    expected_m_tf_coil_copper: Any = None

    expected_whtcon: Any = None

    expected_m_tf_coil_wp_turn_insulation: Any = None

    expected_cplen: Any = None


@pytest.mark.parametrize(
    "tfcoilareaandmassesparam",
    (
        TfCoilAreaAndMassesParam(
            hr1=0.0,
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_inboard=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            z_tf_inside_half=9.0730900215620327,
            den_steel=7800.0,
            m_tf_wp_steel_conduit=0.0,
            m_tf_coils_total=0.0,
            m_tf_coil_case=0.0,
            tficrn=0.0,
            tfcryoarea=0.0,
            m_tf_coil_wp_insulation=0.0,
            tfocrn=0.0,
            m_tf_coil_superconductor=0.0,
            m_tf_coil_copper=0.0,
            m_tf_coil_conductor=0.0,
            m_tf_coil_wp_turn_insulation=0.0,
            f_a_tf_turn_cable_space_extra_void=0.30000000000000004,
            dcond=np.array(
                np.array(
                    (
                        6080.0,
                        6080.0,
                        6070.0,
                        6080.0,
                        6080.0,
                        8500.0,
                        6070.0,
                        8500.0,
                        8500.0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            den_tf_wp_turn_insulation=1800.0,
            len_tf_coil=50.483843027201402,
            den_tf_coil_case=8000.0,
            a_tf_turn_steel=0.0014685061538103825,
            n_tf_coil_turns=200.0,
            n_tf_coils=16.0,
            a_tf_coil_wp_turn_insulation=0.087880174466980876,
            a_tf_coil_outboard_case=1.2752592893394648,
            a_tf_coil_inboard_case=1.0015169239205168,
            f_a_tf_turn_cable_copper=0.80884,
            a_tf_wp_coolant_channels=0.015707963267948974,
            a_tf_turn_cable_space_no_void=0.001293323051622732,
            whttflgs=0.0,
            whtcp=0.0,
            whtconal=0.0,
            vol_cond_cp=0.0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            a_tf_leg_outboard=1.9805354702921749,
            dr_tf_nose_case=0.52465000000000006,
            voltfleg=0.0,
            cplen=0.0,
            itart=0,
            a_tf_wp_with_insulation=0.70527618095271016,
            a_tf_wp_no_insulation=0.64024601555360383,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_m_tf_wp_steel_conduit=115651.90127937049,
            expected_m_tf_coil_case=1034021.9996272125,
            expected_m_tf_coil_wp_insulation=5909.3507916745702,
            expected_m_tf_coil_superconductor=5802.5700395134345,
            expected_m_tf_coil_copper=58744.465423173802,
            expected_whtcon=188184.68882144717,
            expected_m_tf_coil_wp_turn_insulation=7985.7520793894437,
            expected_cplen=20.562180043124066,
        ),
        TfCoilAreaAndMassesParam(
            hr1=0,
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_inboard=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            z_tf_inside_half=9.0730900215620327,
            den_steel=7800.0,
            m_tf_wp_steel_conduit=115651.90127937049,
            m_tf_coils_total=19649856.627845347,
            m_tf_coil_case=1034021.9996272125,
            tficrn=0.8197580588957678,
            tfcryoarea=6381.2092203414386,
            m_tf_coil_wp_insulation=5909.3507916745702,
            tfocrn=0.59553192892551199,
            m_tf_coil_superconductor=5802.5700395134345,
            m_tf_coil_copper=58744.465423173802,
            m_tf_coil_conductor=0.0,
            m_tf_coil_wp_turn_insulation=0.0,
            f_a_tf_turn_cable_space_extra_void=0.30000000000000004,
            dcond=np.array(
                np.array(
                    (
                        6080.0,
                        6080.0,
                        6070.0,
                        6080.0,
                        6080.0,
                        8500.0,
                        6070.0,
                        8500.0,
                        8500.0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            den_tf_wp_turn_insulation=1800.0,
            len_tf_coil=50.514015976170839,
            den_tf_coil_case=8000.0,
            a_tf_turn_steel=0.0014685061538103825,
            n_tf_coil_turns=200,
            n_tf_coils=16,
            a_tf_coil_wp_turn_insulation=0.087880174466980876,
            a_tf_coil_outboard_case=1.2752592893394648,
            a_tf_coil_inboard_case=1.0015169239205168,
            f_a_tf_turn_cable_copper=0.80884,
            a_tf_wp_coolant_channels=0.015707963267948974,
            a_tf_turn_cable_space_no_void=0.001293323051622732,
            whttflgs=0.0,
            whtcp=0.0,
            whtconal=0.0,
            vol_cond_cp=0.0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            a_tf_leg_outboard=1.9805354702921749,
            dr_tf_nose_case=0.52465000000000006,
            voltfleg=0.0,
            cplen=20.562180043124066,
            itart=0,
            a_tf_wp_with_insulation=0.70527618095271016,
            a_tf_wp_no_insulation=0.64024601555360383,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_m_tf_wp_steel_conduit=115721.02357090525,
            expected_m_tf_coil_case=1034699.2182961091,
            expected_m_tf_coil_wp_insulation=5912.8826650262808,
            expected_m_tf_coil_superconductor=5806.038092640837,
            expected_m_tf_coil_copper=58779.575542593491,
            expected_whtcon=188297.16217276,
            expected_m_tf_coil_wp_turn_insulation=7990.5249666247555,
            expected_cplen=20.562180043124066,
        ),
    ),
)
def test_superconducting_tf_coil_area_and_masses(
    tfcoilareaandmassesparam, monkeypatch, sctfcoil
):
    """
    Automatically generated Regression Unit Test for tf_coil_area_and_masses.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcoilareaandmassesparam: the data used to mock and assert in this test.
    :type tfcoilareaandmassesparam: tfcoilareaandmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "hr1", tfcoilareaandmassesparam.hr1)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfcoilareaandmassesparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_inboard", tfcoilareaandmassesparam.dr_tf_inboard
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_mid", tfcoilareaandmassesparam.r_tf_inboard_mid
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcoilareaandmassesparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcoilareaandmassesparam.r_tf_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "z_tf_inside_half", tfcoilareaandmassesparam.z_tf_inside_half
    )

    monkeypatch.setattr(fwbs_variables, "den_steel", tfcoilareaandmassesparam.den_steel)

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_wp_steel_conduit",
        tfcoilareaandmassesparam.m_tf_wp_steel_conduit,
    )

    monkeypatch.setattr(
        tfcoil_variables, "m_tf_coils_total", tfcoilareaandmassesparam.m_tf_coils_total
    )

    monkeypatch.setattr(
        tfcoil_variables, "m_tf_coil_case", tfcoilareaandmassesparam.m_tf_coil_case
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coil_wp_insulation",
        tfcoilareaandmassesparam.m_tf_coil_wp_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coil_superconductor",
        tfcoilareaandmassesparam.m_tf_coil_superconductor,
    )

    monkeypatch.setattr(
        tfcoil_variables, "m_tf_coil_copper", tfcoilareaandmassesparam.m_tf_coil_copper
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coil_conductor",
        tfcoilareaandmassesparam.m_tf_coil_conductor,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coil_wp_turn_insulation",
        tfcoilareaandmassesparam.m_tf_coil_wp_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables, "tfcryoarea", tfcoilareaandmassesparam.tfcryoarea
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "f_a_tf_turn_cable_space_extra_void",
        tfcoilareaandmassesparam.f_a_tf_turn_cable_space_extra_void,
    )

    monkeypatch.setattr(tfcoil_variables, "dcond", tfcoilareaandmassesparam.dcond)

    monkeypatch.setattr(
        tfcoil_variables,
        "den_tf_wp_turn_insulation",
        tfcoilareaandmassesparam.den_tf_wp_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables, "len_tf_coil", tfcoilareaandmassesparam.len_tf_coil
    )

    monkeypatch.setattr(
        tfcoil_variables, "den_tf_coil_case", tfcoilareaandmassesparam.den_tf_coil_case
    )

    monkeypatch.setattr(
        tfcoil_variables, "a_tf_turn_steel", tfcoilareaandmassesparam.a_tf_turn_steel
    )

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coil_turns", tfcoilareaandmassesparam.n_tf_coil_turns
    )

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coils", tfcoilareaandmassesparam.n_tf_coils
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_coil_wp_turn_insulation",
        tfcoilareaandmassesparam.a_tf_coil_wp_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_coil_outboard_case",
        tfcoilareaandmassesparam.a_tf_coil_outboard_case,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_coil_inboard_case",
        tfcoilareaandmassesparam.a_tf_coil_inboard_case,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "f_a_tf_turn_cable_copper",
        tfcoilareaandmassesparam.f_a_tf_turn_cable_copper,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_wp_coolant_channels",
        tfcoilareaandmassesparam.a_tf_wp_coolant_channels,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_turn_cable_space_no_void",
        tfcoilareaandmassesparam.a_tf_turn_cable_space_no_void,
    )

    monkeypatch.setattr(tfcoil_variables, "whttflgs", tfcoilareaandmassesparam.whttflgs)

    monkeypatch.setattr(tfcoil_variables, "whtcp", tfcoilareaandmassesparam.whtcp)

    monkeypatch.setattr(tfcoil_variables, "whtconal", tfcoilareaandmassesparam.whtconal)

    monkeypatch.setattr(
        tfcoil_variables, "vol_cond_cp", tfcoilareaandmassesparam.vol_cond_cp
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfcoilareaandmassesparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_sc_mat", tfcoilareaandmassesparam.i_tf_sc_mat
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_leg_outboard",
        tfcoilareaandmassesparam.a_tf_leg_outboard,
    )

    monkeypatch.setattr(
        tfcoil_variables, "dr_tf_nose_case", tfcoilareaandmassesparam.dr_tf_nose_case
    )

    monkeypatch.setattr(tfcoil_variables, "voltfleg", tfcoilareaandmassesparam.voltfleg)

    monkeypatch.setattr(tfcoil_variables, "cplen", tfcoilareaandmassesparam.cplen)

    monkeypatch.setattr(physics_variables, "itart", tfcoilareaandmassesparam.itart)

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "a_tf_wp_with_insulation",
        tfcoilareaandmassesparam.a_tf_wp_with_insulation,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "a_tf_wp_no_insulation",
        tfcoilareaandmassesparam.a_tf_wp_no_insulation,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "vol_ins_cp",
        tfcoilareaandmassesparam.vol_ins_cp,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "vol_gr_ins_cp",
        tfcoilareaandmassesparam.vol_gr_ins_cp,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "vol_case_cp",
        tfcoilareaandmassesparam.vol_case_cp,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "rad_tf_coil_inboard_toroidal_half",
        tfcoilareaandmassesparam.rad_tf_coil_inboard_toroidal_half,
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "tan_theta_coil",
        tfcoilareaandmassesparam.tan_theta_coil,
    )

    sctfcoil.superconducting_tf_coil_areas_and_masses()

    assert tfcoil_variables.m_tf_wp_steel_conduit == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_wp_steel_conduit
    )

    assert tfcoil_variables.m_tf_coil_case == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_coil_case
    )

    assert tfcoil_variables.m_tf_coil_wp_insulation == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_coil_wp_insulation
    )

    assert tfcoil_variables.m_tf_coil_superconductor == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_coil_superconductor
    )

    assert tfcoil_variables.m_tf_coil_copper == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_coil_copper
    )

    assert tfcoil_variables.m_tf_coil_conductor == pytest.approx(
        tfcoilareaandmassesparam.expected_whtcon
    )

    assert tfcoil_variables.m_tf_coil_wp_turn_insulation == pytest.approx(
        tfcoilareaandmassesparam.expected_m_tf_coil_wp_turn_insulation
    )

    assert tfcoil_variables.cplen == pytest.approx(
        tfcoilareaandmassesparam.expected_cplen
    )


@pytest.mark.parametrize(
    "i_tf_superconductor, j_superconductor, b_tf_inboard_peak, strain, bc20m, tc0m, c0, temp_tf_coolant_peak_field, expected_margin",
    [
        # ITER Nb3Sn, standard parameters
        (1, 1e8, 12.0, 0.0, 32.97, 16.06, 1e10, 4.5, 5.679499736095401),
        # NbTi
        (3, 1e8, 8.0, 0.0, 15.0, 9.3, 1e10, 4.5, 1.3048296694055175),
        # User-defined Nb3Sn
        (4, 1e8, 10.0, 0.0, 30.0, 15.0, 1e10, 4.5, 5.539631803535094),
        # WST Nb3Sn
        (5, 1e8, 13.0, 0.0, 32.97, 16.06, 1e10, 4.5, 5.221287311831414),
        # Durham Ginzburg-Landau Nb-Ti
        (7, 1e8, 7.0, 0.0, 14.85, 9.04, 1e10, 4.5, 1.263064155425198),
        # Durham Ginzburg-Landau REBCO
        (8, 1e8, 10.0, 0.0, 430, 185, 1e10, 20.0, 31.82616792800119),
        # Hazelton-Zhai REBCO
        (9, 1e8, 10.0, 0.0, 138, 92, 1e10, 20.0, 48.363687012510425),
    ],
)
def test_calculate_superconductor_temperature_margin(
    i_tf_superconductor,
    j_superconductor,
    b_tf_inboard_peak,
    strain,
    bc20m,
    tc0m,
    c0,
    temp_tf_coolant_peak_field,
    expected_margin,
):
    sctfcoil = SuperconductingTFCoil()
    margin = sctfcoil.calculate_superconductor_temperature_margin(
        i_tf_superconductor=i_tf_superconductor,
        j_superconductor=j_superconductor,
        b_tf_inboard_peak=b_tf_inboard_peak,
        strain=strain,
        bc20m=bc20m,
        tc0m=tc0m,
        c0=c0,
        temp_tf_coolant_peak_field=temp_tf_coolant_peak_field,
    )
    # The expected_margin values are illustrative; in real tests, use values from reference calculations.
    assert margin == pytest.approx(expected_margin)
