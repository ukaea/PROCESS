from typing import Any, NamedTuple

import pytest

from process.data_structure import (
    build_variables,
    physics_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
)
from process.models.tfcoil.resistive_tf_coil import ResistiveTFCoil


@pytest.fixture
def resistive_tf_coil():
    """Provides SuperconductingTFCoil object for testing.

    :returns: initialised SuperconductingTFCoil object
    :rtype: process.sctfcoil.SuperconductingTFCoil
    """
    return ResistiveTFCoil()


class ResTfInternalGeomParam(NamedTuple):
    n_tf_coil_turns: Any = None

    dx_tf_turn_insulation: Any = None

    dr_tf_nose_case: Any = None

    dr_tf_wp_with_insulation: Any = None

    dx_tf_inboard_out_toroidal: Any = None

    a_tf_inboard_total: Any = None

    c_tf_total: Any = None

    fcoolcp: Any = None

    c_tf_turn: Any = None

    cdtfleg: Any = None

    dr_tf_plasma_case: Any = None

    a_tf_coil_wp_turn_insulation: Any = None

    a_tf_coil_inboard_case: Any = None

    dx_tf_wp_insulation: Any = None

    n_tf_coils: Any = None

    dr_tf_outboard: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    r_cp_top: Any = None

    itart: Any = None

    expected_n_tf_coil_turns: Any = None

    expected_cpttf: Any = None

    expected_cdtfleg: Any = None

    expected_a_tf_coil_wp_turn_insulation: Any = None

    expected_a_tf_coil_inboard_case: Any = None


@pytest.mark.parametrize(
    "restfinternalgeomparam",
    (
        ResTfInternalGeomParam(
            n_tf_coil_turns=0,
            dx_tf_turn_insulation=0.00080000000000000004,
            dr_tf_nose_case=0,
            dr_tf_wp_with_insulation=0.15483000000000002,
            dx_tf_inboard_out_toroidal=0.45367650933034859,
            a_tf_inboard_total=0.0753112923616783,
            c_tf_total=25500000,
            fcoolcp=0.12725,
            c_tf_turn=70000,
            cdtfleg=0,
            dr_tf_plasma_case=0.0077415000000000019,
            a_tf_coil_wp_turn_insulation=0,
            a_tf_coil_inboard_case=0,
            dx_tf_wp_insulation=0,
            n_tf_coils=12,
            dr_tf_outboard=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            itart=1,
            expected_n_tf_coil_turns=1,
            expected_cpttf=2125000,
            expected_cdtfleg=421788350.27812088,
            expected_a_tf_coil_wp_turn_insulation=0.00030678028680367151,
            expected_a_tf_coil_inboard_case=0.00061190425043863676,
        ),
        ResTfInternalGeomParam(
            n_tf_coil_turns=1,
            dx_tf_turn_insulation=0.00080000000000000004,
            dr_tf_nose_case=0,
            dr_tf_wp_with_insulation=0.14708850000000001,
            dx_tf_inboard_out_toroidal=0.44435902370665786,
            a_tf_inboard_total=0.0753112923616783,
            c_tf_total=25500000,
            fcoolcp=0.12725,
            c_tf_turn=2125000,
            cdtfleg=421788350.27812088,
            dr_tf_plasma_case=0.0077415000000000019,
            a_tf_coil_wp_turn_insulation=0.00030678028680367151,
            a_tf_coil_inboard_case=0.00061190425043863676,
            dx_tf_wp_insulation=0,
            n_tf_coils=12,
            dr_tf_outboard=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            itart=1,
            expected_n_tf_coil_turns=1,
            expected_cpttf=2125000,
            expected_cdtfleg=430664525.98439038,
            expected_a_tf_coil_wp_turn_insulation=0.00029439388680367086,
            expected_a_tf_coil_inboard_case=0.00061190425043863676,
        ),
    ),
)
def test_res_tf_internal_geom(restfinternalgeomparam, monkeypatch, resistive_tf_coil):
    """
    Automatically generated Regression Unit Test for res_tf_internal_geom.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param restfinternalgeomparam: the data used to mock and assert in this test.
    :type restfinternalgeomparam: restfinternalgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coil_turns", restfinternalgeomparam.n_tf_coil_turns
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_turn_insulation",
        restfinternalgeomparam.dx_tf_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables, "dr_tf_nose_case", restfinternalgeomparam.dr_tf_nose_case
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dr_tf_wp_with_insulation",
        restfinternalgeomparam.dr_tf_wp_with_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_inboard_out_toroidal",
        restfinternalgeomparam.dx_tf_inboard_out_toroidal,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_inboard_total",
        restfinternalgeomparam.a_tf_inboard_total,
    )

    monkeypatch.setattr(
        tfcoil_variables, "c_tf_total", restfinternalgeomparam.c_tf_total
    )

    monkeypatch.setattr(tfcoil_variables, "fcoolcp", restfinternalgeomparam.fcoolcp)

    monkeypatch.setattr(tfcoil_variables, "c_tf_turn", restfinternalgeomparam.c_tf_turn)

    monkeypatch.setattr(tfcoil_variables, "cdtfleg", restfinternalgeomparam.cdtfleg)

    monkeypatch.setattr(
        tfcoil_variables, "dr_tf_plasma_case", restfinternalgeomparam.dr_tf_plasma_case
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_coil_wp_turn_insulation",
        restfinternalgeomparam.a_tf_coil_wp_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_coil_inboard_case",
        restfinternalgeomparam.a_tf_coil_inboard_case,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_wp_insulation",
        restfinternalgeomparam.dx_tf_wp_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coils", restfinternalgeomparam.n_tf_coils
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_outboard", restfinternalgeomparam.dr_tf_outboard
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", restfinternalgeomparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", restfinternalgeomparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "r_cp_top", restfinternalgeomparam.r_cp_top)

    monkeypatch.setattr(physics_variables, "itart", restfinternalgeomparam.itart)

    resistive_tf_coil.res_tf_internal_geom()

    assert tfcoil_variables.n_tf_coil_turns == pytest.approx(
        restfinternalgeomparam.expected_n_tf_coil_turns
    )

    assert tfcoil_variables.c_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.a_tf_coil_wp_turn_insulation == pytest.approx(
        restfinternalgeomparam.expected_a_tf_coil_wp_turn_insulation
    )

    assert tfcoil_variables.a_tf_coil_inboard_case == pytest.approx(
        restfinternalgeomparam.expected_a_tf_coil_inboard_case
    )

    assert tfcoil_variables.n_tf_coil_turns == pytest.approx(
        restfinternalgeomparam.expected_n_tf_coil_turns
    )

    assert tfcoil_variables.c_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.a_tf_coil_wp_turn_insulation == pytest.approx(
        restfinternalgeomparam.expected_a_tf_coil_wp_turn_insulation
    )

    assert tfcoil_variables.a_tf_coil_inboard_case == pytest.approx(
        restfinternalgeomparam.expected_a_tf_coil_inboard_case
    )


class TfResHeatingParam(NamedTuple):
    rho_cp: Any = None
    temp_tf_legs_outboard: Any = None
    dx_tf_turn_insulation: Any = None
    th_joint_contact: Any = None
    rho_tf_leg: Any = None
    vol_cond_cp: Any = None
    n_tf_coil_turns: Any = None
    dr_tf_nose_case: Any = None
    dx_tf_inboard_out_toroidal: Any = None
    len_tf_coil: Any = None
    res_tf_leg: Any = None
    temp_cp_average: Any = None
    a_tf_leg_outboard: Any = None
    c_tf_total: Any = None
    rho_tf_joints: Any = None
    p_tf_leg_resistive: Any = None
    p_cp_resistive: Any = None
    p_tf_joints_resistive: Any = None
    n_tf_joints_contact: Any = None
    n_tf_joints: Any = None
    n_tf_coils: Any = None
    i_tf_sup: Any = None
    frholeg: Any = None
    frhocp: Any = None
    fcoolcp: Any = None
    dr_tf_plasma_case: Any = None
    a_cp_cool: Any = None
    f_a_tf_cool_outboard: Any = None
    i_cp_joints: Any = None
    dx_tf_wp_insulation: Any = None
    dr_tf_outboard: Any = None
    dr_tf_inboard: Any = None
    r_cp_top: Any = None
    z_tf_inside_half: Any = None
    r_tf_inboard_in: Any = None
    r_tf_inboard_out: Any = None
    itart: Any = None
    z_cp_top: Any = None
    is_leg_cp_temp_same: Any = None
    expected_rho_cp: Any = None
    expected_rho_tf_leg: Any = None
    expected_vol_cond_cp: Any = None
    expected_res_tf_leg: Any = None
    expected_p_tf_leg_resistive: Any = None
    expected_p_cp_resistive: Any = None
    expected_pres_joints: Any = None
    expected_a_cp_cool: Any = None
    expected_is_leg_cp_temp_same: Any = None


@pytest.mark.parametrize(
    "tfresheatingparam",
    (
        TfResHeatingParam(
            rho_cp=0,
            temp_tf_legs_outboard=-1,
            dx_tf_turn_insulation=0.00080000000000000004,
            th_joint_contact=0.029999999999999999,
            rho_tf_leg=0,
            vol_cond_cp=0,
            n_tf_coil_turns=1,
            dr_tf_nose_case=0,
            dx_tf_inboard_out_toroidal=0.45367650933034859,
            len_tf_coil=15.582502857142856,
            res_tf_leg=0,
            temp_cp_average=347.13,
            a_tf_leg_outboard=0.070242733939617885,
            c_tf_total=25500000,
            rho_tf_joints=2.5000000000000002e-10,
            p_tf_leg_resistive=0,
            p_cp_resistive=0,
            p_tf_joints_resistive=0,
            n_tf_joints_contact=6,
            n_tf_joints=4,
            n_tf_coils=12,
            i_tf_sup=0,
            frholeg=1,
            frhocp=1,
            fcoolcp=0.12725,
            dr_tf_plasma_case=0.0077415000000000019,
            a_cp_cool=0,
            f_a_tf_cool_outboard=0.20000000000000001,
            i_cp_joints=1,
            dx_tf_wp_insulation=0,
            dr_tf_outboard=0.15483000000000002,
            dr_tf_inboard=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            z_tf_inside_half=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            itart=1,
            z_cp_top=2.6714285714285717,
            is_leg_cp_temp_same=0,
            expected_rho_cp=2.0721414e-08,
            expected_rho_tf_leg=2.0721414e-08,
            expected_vol_cond_cp=12.020160732580297,
            expected_res_tf_leg=5.826541893267926e-06,
            expected_p_tf_leg_resistive=315725738.84145576,
            expected_p_cp_resistive=446175692.00121117,
            expected_pres_joints=1944336.7995005273,
            expected_a_cp_cool=0.00068328705812121333,
            expected_is_leg_cp_temp_same=1,
        ),
        TfResHeatingParam(
            rho_cp=2.0721414e-08,
            temp_tf_legs_outboard=-1,
            dx_tf_turn_insulation=0.00080000000000000004,
            th_joint_contact=0.029999999999999999,
            rho_tf_leg=2.0721414e-08,
            vol_cond_cp=12.020160732580297,
            n_tf_coil_turns=1,
            dr_tf_nose_case=0,
            dx_tf_inboard_out_toroidal=0.44435902370665786,
            len_tf_coil=15.654502857142857,
            res_tf_leg=5.647653956699231e-06,
            temp_cp_average=347.13,
            a_tf_leg_outboard=0.068800107640501845,
            c_tf_total=25500000,
            rho_tf_joints=2.5000000000000002e-10,
            p_tf_leg_resistive=332643748.67243439,
            p_cp_resistive=432477095.0716282,
            p_tf_joints_resistive=1944336.7995005273,
            n_tf_joints_contact=6,
            n_tf_joints=4,
            n_tf_coils=12,
            i_tf_sup=0,
            frholeg=1,
            frhocp=1,
            fcoolcp=0.12725,
            dr_tf_plasma_case=0.0077415000000000019,
            a_cp_cool=0.00068328705812121333,
            f_a_tf_cool_outboard=0.20000000000000001,
            i_cp_joints=1,
            dx_tf_wp_insulation=0,
            dr_tf_outboard=0.15483000000000002,
            dr_tf_inboard=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            z_tf_inside_half=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            itart=1,
            z_cp_top=2.6714285714285717,
            is_leg_cp_temp_same=1,
            expected_rho_cp=2.0721414e-08,
            expected_rho_tf_leg=2.0721414e-08,
            expected_vol_cond_cp=11.545770024935592,
            expected_res_tf_leg=5.9766449694251195e-06,
            expected_p_tf_leg_resistive=323859449.2807237,
            expected_p_cp_resistive=451516213.33864105,
            expected_pres_joints=1944336.7995005273,
            expected_a_cp_cool=0.00068328705812121333,
            expected_is_leg_cp_temp_same=1,
        ),
    ),
)
def test_tf_res_heating(tfresheatingparam, monkeypatch, resistive_tf_coil):
    """
    Automatically generated Regression Unit Test for tf_res_heating.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param tfresheatingparam: the data used to mock and assert in this test.
    :type tfresheatingparam: tfresheatingparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "rho_cp", tfresheatingparam.rho_cp)

    monkeypatch.setattr(
        tfcoil_variables,
        "temp_tf_legs_outboard",
        tfresheatingparam.temp_tf_legs_outboard,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_turn_insulation",
        tfresheatingparam.dx_tf_turn_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables, "th_joint_contact", tfresheatingparam.th_joint_contact
    )

    monkeypatch.setattr(tfcoil_variables, "rho_tf_leg", tfresheatingparam.rho_tf_leg)

    monkeypatch.setattr(tfcoil_variables, "vol_cond_cp", tfresheatingparam.vol_cond_cp)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coil_turns", tfresheatingparam.n_tf_coil_turns
    )

    monkeypatch.setattr(
        tfcoil_variables, "dr_tf_nose_case", tfresheatingparam.dr_tf_nose_case
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_inboard_out_toroidal",
        tfresheatingparam.dx_tf_inboard_out_toroidal,
    )

    monkeypatch.setattr(tfcoil_variables, "len_tf_coil", tfresheatingparam.len_tf_coil)

    monkeypatch.setattr(tfcoil_variables, "res_tf_leg", tfresheatingparam.res_tf_leg)

    monkeypatch.setattr(
        tfcoil_variables, "temp_cp_average", tfresheatingparam.temp_cp_average
    )

    monkeypatch.setattr(
        tfcoil_variables, "a_tf_leg_outboard", tfresheatingparam.a_tf_leg_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "c_tf_total", tfresheatingparam.c_tf_total)

    monkeypatch.setattr(
        tfcoil_variables, "rho_tf_joints", tfresheatingparam.rho_tf_joints
    )

    monkeypatch.setattr(
        tfcoil_variables, "p_tf_leg_resistive", tfresheatingparam.p_tf_leg_resistive
    )

    monkeypatch.setattr(
        tfcoil_variables, "p_cp_resistive", tfresheatingparam.p_cp_resistive
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "p_tf_joints_resistive",
        tfresheatingparam.p_tf_joints_resistive,
    )

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_joints_contact", tfresheatingparam.n_tf_joints_contact
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf_joints", tfresheatingparam.n_tf_joints)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", tfresheatingparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfresheatingparam.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "frholeg", tfresheatingparam.frholeg)

    monkeypatch.setattr(tfcoil_variables, "frhocp", tfresheatingparam.frhocp)

    monkeypatch.setattr(tfcoil_variables, "fcoolcp", tfresheatingparam.fcoolcp)

    monkeypatch.setattr(
        tfcoil_variables, "dr_tf_plasma_case", tfresheatingparam.dr_tf_plasma_case
    )

    monkeypatch.setattr(tfcoil_variables, "a_cp_cool", tfresheatingparam.a_cp_cool)

    monkeypatch.setattr(
        tfcoil_variables, "f_a_tf_cool_outboard", tfresheatingparam.f_a_tf_cool_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "i_cp_joints", tfresheatingparam.i_cp_joints)

    monkeypatch.setattr(
        tfcoil_variables, "dx_tf_wp_insulation", tfresheatingparam.dx_tf_wp_insulation
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_outboard", tfresheatingparam.dr_tf_outboard
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_inboard", tfresheatingparam.dr_tf_inboard
    )

    monkeypatch.setattr(build_variables, "r_cp_top", tfresheatingparam.r_cp_top)

    monkeypatch.setattr(
        build_variables, "z_tf_inside_half", tfresheatingparam.z_tf_inside_half
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfresheatingparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfresheatingparam.r_tf_inboard_out
    )

    monkeypatch.setattr(physics_variables, "itart", tfresheatingparam.itart)

    monkeypatch.setattr(
        superconducting_tf_coil_variables, "z_cp_top", tfresheatingparam.z_cp_top
    )

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "is_leg_cp_temp_same",
        tfresheatingparam.is_leg_cp_temp_same,
    )

    resistive_tf_coil.tf_res_heating()

    assert tfcoil_variables.rho_cp == pytest.approx(tfresheatingparam.expected_rho_cp)

    assert tfcoil_variables.rho_tf_leg == pytest.approx(
        tfresheatingparam.expected_rho_tf_leg
    )

    assert tfcoil_variables.vol_cond_cp == pytest.approx(
        tfresheatingparam.expected_vol_cond_cp
    )

    assert tfcoil_variables.res_tf_leg == pytest.approx(
        tfresheatingparam.expected_res_tf_leg
    )

    assert tfcoil_variables.p_tf_leg_resistive == pytest.approx(
        tfresheatingparam.expected_p_tf_leg_resistive
    )

    assert tfcoil_variables.p_cp_resistive == pytest.approx(
        tfresheatingparam.expected_p_cp_resistive
    )

    assert tfcoil_variables.p_tf_joints_resistive == pytest.approx(
        tfresheatingparam.expected_pres_joints
    )

    assert tfcoil_variables.a_cp_cool == pytest.approx(
        tfresheatingparam.expected_a_cp_cool
    )

    assert superconducting_tf_coil_variables.is_leg_cp_temp_same == pytest.approx(
        tfresheatingparam.expected_is_leg_cp_temp_same
    )


class CpostParam(NamedTuple):
    n_tf_coils: Any = None

    z_tf_inside_half: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    r_cp_top: Any = None

    ztop: Any = None

    hmaxi: Any = None

    cas_in_th: Any = None

    cas_out_th: Any = None

    gr_ins_th: Any = None

    ins_th: Any = None

    n_tf_coil_turns: Any = None

    curr: Any = None

    rho: Any = None

    fcool: Any = None

    expected_vol_ins_cp: Any = None

    expected_vol_gr_ins_cp: Any = None

    expected_vol_case_cp: Any = None

    expected_respow: Any = None

    expected_vol_cond_cp: Any = None

    expected_a_cp_cool: Any = None


@pytest.mark.parametrize(
    "cpostparam",
    (
        CpostParam(
            n_tf_coils=12,
            z_tf_inside_half=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            ztop=2.6714285714285717,
            hmaxi=4.5762585714285713,
            cas_in_th=0,
            cas_out_th=0.0077415000000000019,
            gr_ins_th=0,
            ins_th=0.00080000000000000004,
            n_tf_coil_turns=1,
            curr=25500000,
            rho=2.1831760869565221e-08,
            fcool=0.12725,
            expected_vol_ins_cp=0.12917075053120922,
            expected_vol_gr_ins_cp=0,
            expected_vol_case_cp=0.12791418544773489,
            expected_respow=470083798.99090022,
            expected_vol_cond_cp=12.020160732580297,
            expected_a_cp_cool=0.00068328705812121333,
        ),
        CpostParam(
            n_tf_coils=12,
            z_tf_inside_half=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            ztop=2.6714285714285717,
            hmaxi=4.5762585714285713,
            cas_in_th=0,
            cas_out_th=0.0077415000000000019,
            gr_ins_th=0,
            ins_th=0.00080000000000000004,
            n_tf_coil_turns=1,
            curr=25500000,
            rho=2.1831760869565221e-08,
            fcool=0.12725,
            expected_vol_ins_cp=0.12679799009998483,
            expected_vol_gr_ins_cp=0,
            expected_vol_case_cp=0.12648575512245444,
            expected_respow=475710489.56122422,
            expected_vol_cond_cp=11.545770024935592,
            expected_a_cp_cool=0.00068328705812121333,
        ),
    ),
)
def test_cpost(cpostparam, monkeypatch, resistive_tf_coil):
    """
    Automatically generated Regression Unit Test for cpost.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param cpostparam: the data used to mock and assert in this test.
    :type cpostparam: cpostparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "z_tf_inside_half", cpostparam.z_tf_inside_half)

    (
        a_cp_cool,
        vol_cond_cp,
        respow,
        vol_ins_cp,
        vol_case_cp,
        vol_gr_ins_cp,
    ) = resistive_tf_coil.cpost(
        r_tf_inboard_in=cpostparam.r_tf_inboard_in,
        r_tf_inboard_out=cpostparam.r_tf_inboard_out,
        r_cp_top=cpostparam.r_cp_top,
        ztop=cpostparam.ztop,
        hmaxi=cpostparam.hmaxi,
        cas_in_th=cpostparam.cas_in_th,
        cas_out_th=cpostparam.cas_out_th,
        gr_ins_th=cpostparam.gr_ins_th,
        ins_th=cpostparam.ins_th,
        n_tf_coil_turns=cpostparam.n_tf_coil_turns,
        curr=cpostparam.curr,
        rho=cpostparam.rho,
        fcool=cpostparam.fcool,
        n_tf_coils=cpostparam.n_tf_coils,
    )

    assert vol_ins_cp == pytest.approx(cpostparam.expected_vol_ins_cp)

    assert vol_gr_ins_cp == pytest.approx(cpostparam.expected_vol_gr_ins_cp)

    assert vol_case_cp == pytest.approx(cpostparam.expected_vol_case_cp)

    assert respow == pytest.approx(cpostparam.expected_respow)

    assert vol_cond_cp == pytest.approx(cpostparam.expected_vol_cond_cp)

    assert a_cp_cool == pytest.approx(cpostparam.expected_a_cp_cool)
