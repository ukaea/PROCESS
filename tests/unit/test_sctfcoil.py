from typing import Any, NamedTuple

import pytest

from process import superconducting_tf_coil as sctf
from process.data_structure import divertor_variables
from process.fortran import (
    build_variables,
    global_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.superconducting_tf_coil import SuperconductingTFCoil


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
            expected_ajwpro=17475706.393616617,
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

    ajwpro, vd = sctfcoil.protect(
        aio=protectparam.aio,
        tfes=protectparam.tfes,
        acs=protectparam.acs,
        aturn=protectparam.aturn,
        tdump=protectparam.tdump,
        fcond=protectparam.fcond,
        fcu=protectparam.fcu,
        tba=protectparam.tba,
        tmax=protectparam.tmax,
    )

    assert ajwpro == pytest.approx(protectparam.expected_ajwpro)

    assert vd == pytest.approx(protectparam.expected_vd)


class SuperconParam(NamedTuple):
    tmargmin_tf: Any = None

    n_tf_coils: Any = None

    temp_margin: Any = None

    jwdgpro: Any = None

    dia_tf_turn_coolant_channel: Any = None

    c_tf_turn: Any = None

    bmaxtfrp: Any = None

    str_tf_con_res: Any = None

    b_crit_upper_nbti: Any = None

    i_str_wp: Any = None

    str_wp: Any = None

    t_crit_nbti: Any = None

    tf_fit_t: Any = None

    tf_fit_z: Any = None

    tf_fit_y: Any = None

    run_tests: Any = None

    i_tf_superconductor: Any = None

    iprint: Any = None

    outfile: Any = None

    a_tf_turn_cable_space: Any = None

    a_tf_turn: Any = None

    b_tf_inboard_peak: Any = None

    f_a_tf_turn_cable_copper: Any = None

    f_a_tf_turn_cooling_extra: Any = None

    f_strain_scale: Any = None

    j_tf_wp: Any = None

    t_tf_quench_dump: Any = None

    e_tf_coil_magnetic_stored: Any = None

    temp_tf_coolant_peak_field: Any = None

    temp_tf_conductor_peak_quench: Any = None

    bcritsc: Any = None

    tcritsc: Any = None

    expected_temp_margin: Any = None

    expected_jwdgpro: Any = None

    expected_j_tf_wp_critical: Any = None

    expected_vd: Any = None

    expected_tmarg: Any = None


@pytest.mark.parametrize(
    "superconparam",
    (
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf_coils=16,
            temp_margin=0,
            jwdgpro=0,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            c_tf_turn=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            i_tf_superconductor=5,
            iprint=0,
            outfile=11,
            a_tf_turn_cable_space=0.001293323051622732,
            a_tf_turn=0.0032012300777680192,
            b_tf_inboard_peak=12.48976756562082,
            f_a_tf_turn_cable_copper=0.80884,
            f_a_tf_turn_cooling_extra=0.30000000000000004,
            f_strain_scale=0.5,
            j_tf_wp=23124470.793774806,
            t_tf_quench_dump=25.829000000000001,
            e_tf_coil_magnetic_stored=9548964780.4287167,
            temp_tf_coolant_peak_field=4.75,
            temp_tf_conductor_peak_quench=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.34312129,
            expected_jwdgpro=17475706.393616617,
            expected_j_tf_wp_critical=41107234.360397324,
            expected_vd=9988.2637896807955,
            expected_tmarg=2.34312129,
        ),
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf_coils=16,
            temp_margin=2.3431632224075836,
            jwdgpro=17475706.393616617,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            c_tf_turn=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            i_tf_superconductor=5,
            iprint=0,
            outfile=11,
            a_tf_turn_cable_space=0.001293323051622732,
            a_tf_turn=0.0032012300777680192,
            b_tf_inboard_peak=12.48976756562082,
            f_a_tf_turn_cable_copper=0.80884,
            f_a_tf_turn_cooling_extra=0.30000000000000004,
            f_strain_scale=0.5,
            j_tf_wp=23124470.793774806,
            t_tf_quench_dump=25.829000000000001,
            e_tf_coil_magnetic_stored=9561415368.8360519,
            temp_tf_coolant_peak_field=4.75,
            temp_tf_conductor_peak_quench=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.34312129,
            expected_jwdgpro=17475706.393616617,
            expected_j_tf_wp_critical=41107234.360397324,
            expected_vd=10001.287165953383,
            expected_tmarg=2.34312129,
        ),
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf_coils=16,
            temp_margin=2.3431632224075836,
            jwdgpro=17475706.393616617,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            c_tf_turn=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            i_tf_superconductor=5,
            iprint=0,
            outfile=11,
            a_tf_turn_cable_space=0.001293323051622732,
            a_tf_turn=0.0032012300777680192,
            b_tf_inboard_peak=12.48976756562082,
            f_a_tf_turn_cable_copper=0.80884,
            f_a_tf_turn_cooling_extra=0.30000000000000004,
            f_strain_scale=0.5,
            j_tf_wp=23124470.793774806,
            t_tf_quench_dump=25.829000000000001,
            e_tf_coil_magnetic_stored=9561415368.8360519,
            temp_tf_coolant_peak_field=4.75,
            temp_tf_conductor_peak_quench=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.34312129,
            expected_jwdgpro=17475706.393616617,
            expected_j_tf_wp_critical=41107234.360397324,
            expected_vd=10001.287165953383,
            expected_tmarg=2.34312129,
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

    monkeypatch.setattr(tfcoil_variables, "tmargmin_tf", superconparam.tmargmin_tf)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", superconparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "temp_margin", superconparam.temp_margin)

    monkeypatch.setattr(tfcoil_variables, "jwdgpro", superconparam.jwdgpro)

    monkeypatch.setattr(
        tfcoil_variables,
        "dia_tf_turn_coolant_channel",
        superconparam.dia_tf_turn_coolant_channel,
    )

    monkeypatch.setattr(tfcoil_variables, "c_tf_turn", superconparam.c_tf_turn)

    monkeypatch.setattr(tfcoil_variables, "bmaxtfrp", superconparam.bmaxtfrp)

    monkeypatch.setattr(
        tfcoil_variables, "str_tf_con_res", superconparam.str_tf_con_res
    )

    monkeypatch.setattr(
        tfcoil_variables, "b_crit_upper_nbti", superconparam.b_crit_upper_nbti
    )

    monkeypatch.setattr(tfcoil_variables, "i_str_wp", superconparam.i_str_wp)

    monkeypatch.setattr(tfcoil_variables, "str_wp", superconparam.str_wp)

    monkeypatch.setattr(tfcoil_variables, "t_crit_nbti", superconparam.t_crit_nbti)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_t", superconparam.tf_fit_t)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_z", superconparam.tf_fit_z)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_y", superconparam.tf_fit_y)

    monkeypatch.setattr(global_variables, "run_tests", superconparam.run_tests)

    j_tf_wp_critical, vd, tmarg = sctfcoil.supercon(
        i_tf_superconductor=superconparam.i_tf_superconductor,
        a_tf_turn_cable_space=superconparam.a_tf_turn_cable_space,
        a_tf_turn=superconparam.a_tf_turn,
        b_tf_inboard_peak=superconparam.b_tf_inboard_peak,
        f_a_tf_turn_cable_copper=superconparam.f_a_tf_turn_cable_copper,
        f_a_tf_turn_cooling_extra=superconparam.f_a_tf_turn_cooling_extra,
        f_strain_scale=superconparam.f_strain_scale,
        c_tf_turn=superconparam.c_tf_turn,
        j_tf_wp=superconparam.j_tf_wp,
        t_tf_quench_dump=superconparam.t_tf_quench_dump,
        e_tf_coil_magnetic_stored=superconparam.e_tf_coil_magnetic_stored,
        temp_tf_coolant_peak_field=superconparam.temp_tf_coolant_peak_field,
        temp_tf_conductor_peak_quench=superconparam.temp_tf_conductor_peak_quench,
        bcritsc=superconparam.bcritsc,
        tcritsc=superconparam.tcritsc,
        output=False,
    )

    assert tfcoil_variables.temp_margin == pytest.approx(
        superconparam.expected_temp_margin
    )

    assert tfcoil_variables.jwdgpro == pytest.approx(superconparam.expected_jwdgpro)

    assert j_tf_wp_critical == pytest.approx(superconparam.expected_j_tf_wp_critical)

    assert vd == pytest.approx(superconparam.expected_vd)

    assert tmarg == pytest.approx(superconparam.expected_tmarg)


class PeakTfWithRippleParam(NamedTuple):
    tf_fit_t: Any = None

    tf_fit_z: Any = None

    tf_fit_y: Any = None

    n_tf_coils: Any = None

    dx_tf_wp_primary_toroidal: Any = None

    dr_tf_wp_with_insulation: Any = None

    tfin: Any = None

    b_tf_inboard_peak: Any = None

    expected_tf_fit_t: Any = None

    expected_tf_fit_z: Any = None

    expected_tf_fit_y: Any = None

    expected_bmaxtfrp: Any = None

    expected_flag: Any = None


@pytest.mark.parametrize(
    "peaktfwithrippleparam",
    (
        PeakTfWithRippleParam(
            tf_fit_t=0,
            tf_fit_z=0,
            tf_fit_y=0,
            n_tf_coils=16,
            dx_tf_wp_primary_toroidal=1.299782604942499,
            dr_tf_wp_with_insulation=0.50661087836601015,
            tfin=3.789896624292115,
            b_tf_inboard_peak=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_tf_fit_y=1.0658869305062604,
            expected_bmaxtfrp=12.48976756562082,
            expected_flag=0,
        ),
        PeakTfWithRippleParam(
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            n_tf_coils=16,
            dx_tf_wp_primary_toroidal=1.299782604942499,
            dr_tf_wp_with_insulation=0.50661087836601015,
            tfin=3.789896624292115,
            b_tf_inboard_peak=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_tf_fit_y=1.0658869305062604,
            expected_bmaxtfrp=12.48976756562082,
            expected_flag=0,
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

    monkeypatch.setattr(sctfcoil_module, "tf_fit_t", peaktfwithrippleparam.tf_fit_t)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_z", peaktfwithrippleparam.tf_fit_z)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_y", peaktfwithrippleparam.tf_fit_y)

    bmaxtfrp, flag = sctfcoil.peak_tf_with_ripple(
        n_tf_coils=peaktfwithrippleparam.n_tf_coils,
        dx_tf_wp_primary_toroidal=peaktfwithrippleparam.dx_tf_wp_primary_toroidal,
        dr_tf_wp_with_insulation=peaktfwithrippleparam.dr_tf_wp_with_insulation,
        tfin=peaktfwithrippleparam.tfin,
        b_tf_inboard_peak=peaktfwithrippleparam.b_tf_inboard_peak,
    )

    assert sctfcoil_module.tf_fit_t == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_t
    )

    assert sctfcoil_module.tf_fit_z == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_z
    )

    assert sctfcoil_module.tf_fit_y == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_y
    )

    assert bmaxtfrp == pytest.approx(peaktfwithrippleparam.expected_bmaxtfrp)

    assert flag == pytest.approx(peaktfwithrippleparam.expected_flag)


class TfWpGeomParam(NamedTuple):
    dr_tf_inboard: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    dr_tf_wp_with_insulation: Any = None

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

    assert a_tf_plasma_case == pytest.approx(tfcasegeomparam.expected_a_tf_plasma_case)

    assert a_tf_coil_nose_case == pytest.approx(
        tfcasegeomparam.expected_a_tf_coil_nose_case
    )


class TfIntegerTurnGeomParam(NamedTuple):
    dr_tf_wp_with_insulation: Any = None

    dx_tf_wp_insulation: Any = None

    dx_tf_wp_insertion_gap: Any = None

    t_conductor: Any = None

    t_turn_tf: Any = None

    c_tf_coil: Any = None

    dx_tf_wp_toroidal_min: Any = None

    t_conductor_radial: Any = None

    t_conductor_toroidal: Any = None

    dr_tf_turn_cable_space: Any = None

    dx_tf_turn_cable_space: Any = None

    dr_tf_turn: Any = None

    dx_tf_turn: Any = None

    dx_tf_turn_cable_space_average: Any = None

    n_layer: Any = None

    n_pancake: Any = None

    dx_tf_turn_steel: Any = None

    dx_tf_turn_insulation: Any = None

    expected_t_conductor: Any = None

    expected_t_turn_tf: Any = None

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
            t_turn_tf=0,
            c_tf_coil=14805350.287500001,
            dx_tf_wp_toroidal_min=1.299782604942499,
            t_conductor_radial=0,
            t_conductor_toroidal=0,
            dr_tf_turn_cable_space=0,
            dx_tf_turn_cable_space=0,
            dr_tf_turn=0,
            dx_tf_turn=0,
            dx_tf_turn_cable_space_average=0,
            n_layer=10,
            n_pancake=20,
            dx_tf_turn_steel=0.0080000000000000002,
            dx_tf_turn_insulation=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_t_turn_tf=0.056579413904423038,
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
            t_turn_tf=0.056579413904423038,
            c_tf_coil=14805350.287500001,
            dx_tf_wp_toroidal_min=1.299782604942499,
            t_conductor_radial=0.046661087836601015,
            t_conductor_toroidal=0.059189130247124938,
            dr_tf_turn_cable_space=0.030661087836601014,
            dx_tf_turn_cable_space=0.043189130247124938,
            dr_tf_turn=0.050661087836601018,
            dx_tf_turn=0.063189130247124942,
            dx_tf_turn_cable_space_average=0.036389912284773368,
            n_layer=10,
            n_pancake=20,
            dx_tf_turn_steel=0.0080000000000000002,
            dx_tf_turn_insulation=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_t_turn_tf=0.056579413904423038,
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
def test_tf_integer_turn_geom(tfintegerturngeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_integer_turn_geom.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfintegerturngeomparam: the data used to mock and assert in this test.
    :type tfintegerturngeomparam: tfintegerturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        tfcoil_variables,
        "dr_tf_wp_with_insulation",
        tfintegerturngeomparam.dr_tf_wp_with_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_wp_insulation",
        tfintegerturngeomparam.dx_tf_wp_insulation,
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_wp_insertion_gap",
        tfintegerturngeomparam.dx_tf_wp_insertion_gap,
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_conductor", tfintegerturngeomparam.t_conductor
    )

    monkeypatch.setattr(tfcoil_variables, "t_turn_tf", tfintegerturngeomparam.t_turn_tf)

    monkeypatch.setattr(sctfcoil_module, "c_tf_coil", tfintegerturngeomparam.c_tf_coil)

    monkeypatch.setattr(
        sctfcoil_module,
        "dx_tf_wp_toroidal_min",
        tfintegerturngeomparam.dx_tf_wp_toroidal_min,
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_conductor_radial", tfintegerturngeomparam.t_conductor_radial
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "t_conductor_toroidal",
        tfintegerturngeomparam.t_conductor_toroidal,
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "dr_tf_turn_cable_space",
        tfintegerturngeomparam.dr_tf_turn_cable_space,
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "dx_tf_turn_cable_space",
        tfintegerturngeomparam.dx_tf_turn_cable_space,
    )

    monkeypatch.setattr(
        sctfcoil_module, "dr_tf_turn", tfintegerturngeomparam.dr_tf_turn
    )

    monkeypatch.setattr(
        sctfcoil_module, "dx_tf_turn", tfintegerturngeomparam.dx_tf_turn
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "dx_tf_turn_cable_space_average",
        tfintegerturngeomparam.dx_tf_turn_cable_space_average,
    )

    (
        a_tf_turn_cable_space_no_void,
        a_tf_turn_steel,
        a_tf_turn_insulation,
        c_tf_turn,
        n_tf_coil_turns,
    ) = sctfcoil.tf_integer_turn_geom(
        n_layer=tfintegerturngeomparam.n_layer,
        n_pancake=tfintegerturngeomparam.n_pancake,
        dx_tf_turn_steel=tfintegerturngeomparam.dx_tf_turn_steel,
        dx_tf_turn_insulation=tfintegerturngeomparam.dx_tf_turn_insulation,
    )

    assert tfcoil_variables.t_conductor == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor
    )

    assert tfcoil_variables.t_turn_tf == pytest.approx(
        tfintegerturngeomparam.expected_t_turn_tf
    )

    assert sctfcoil_module.t_conductor_radial == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_radial
    )

    assert sctfcoil_module.t_conductor_toroidal == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_toroidal
    )

    assert sctfcoil_module.dr_tf_turn_cable_space == pytest.approx(
        tfintegerturngeomparam.expected_dr_tf_turn_cable_space
    )

    assert sctfcoil_module.dx_tf_turn_cable_space == pytest.approx(
        tfintegerturngeomparam.expected_dx_tf_turn_cable_space
    )

    assert sctfcoil_module.dr_tf_turn == pytest.approx(
        tfintegerturngeomparam.expected_t_turn_radial
    )

    assert sctfcoil_module.dx_tf_turn == pytest.approx(
        tfintegerturngeomparam.expected_dx_tf_turn
    )

    assert sctfcoil_module.dx_tf_turn_cable_space_average == pytest.approx(
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

    t_turn_tf: Any = None

    t_turn_tf_is_input: Any = None

    c_tf_turn: Any = None

    t_cable_tf: Any = None

    t_cable_tf_is_input: Any = None

    a_tf_wp_no_insulation: Any = None

    dr_tf_turn: Any = None

    dx_tf_turn: Any = None

    dx_tf_turn_cable_space_average: Any = None

    i_tf_sc_mat: Any = None

    j_tf_wp: Any = None

    dx_tf_turn_steel: Any = None

    dx_tf_turn_insulation: Any = None

    expected_t_conductor: Any = None

    expected_t_turn_tf: Any = None

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
            t_turn_tf=0,
            t_turn_tf_is_input=False,
            c_tf_turn=65000,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0,
            dx_tf_turn=0,
            dx_tf_turn_cable_space_average=0,
            i_tf_sc_mat=5,
            j_tf_wp=26493137.688284047,
            dx_tf_turn_steel=0.0080000000000000019,
            dx_tf_turn_insulation=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_t_turn_tf=0.049532469413859428,
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
            t_turn_tf=0.049532469413859428,
            t_turn_tf_is_input=False,
            c_tf_turn=65000,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.049532469413859428,
            dx_tf_turn=0.049532469413859428,
            dx_tf_turn_cable_space_average=0.031932469413859424,
            i_tf_sc_mat=5,
            j_tf_wp=26493137.688284047,
            dx_tf_turn_steel=0.0080000000000000019,
            dx_tf_turn_insulation=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_t_turn_tf=0.049532469413859428,
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
            t_turn_tf=0.05872,
            t_turn_tf_is_input=True,
            c_tf_turn=0,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.05872,
            dx_tf_turn=0.05872,
            dx_tf_turn_cable_space_average=0.04109,
            i_tf_sc_mat=1,
            j_tf_wp=2.301e07,
            dx_tf_turn_steel=8.015e-03,
            dx_tf_turn_insulation=8.0e-4,
            expected_t_conductor=5.712e-02,
            expected_t_turn_tf=0.05872,
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
            t_turn_tf=0,
            t_turn_tf_is_input=False,
            c_tf_turn=0,
            t_cable_tf=0.042,
            t_cable_tf_is_input=True,
            a_tf_wp_no_insulation=0.60510952642236249,
            dr_tf_turn=0.05872,
            dx_tf_turn=0.05872,
            dx_tf_turn_cable_space_average=0.04109,
            i_tf_sc_mat=1,
            j_tf_wp=2.673e07,
            dx_tf_turn_steel=8.148e-03,
            dx_tf_turn_insulation=8.0e-4,
            expected_t_conductor=0.058296,
            expected_t_turn_tf=0.059896,
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
def test_tf_averaged_turn_geom(tfaveragedturngeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_averaged_turn_geom.

    This test was generated using data from tests/regression/scenarios/i_mode/IN.DAT.

    :param tfaveragedturngeomparam: the data used to mock and assert in this test.
    :type tfaveragedturngeomparam: tfaveragedturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        tfcoil_variables, "layer_ins", tfaveragedturngeomparam.layer_ins
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_conductor", tfaveragedturngeomparam.t_conductor
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_turn_tf", tfaveragedturngeomparam.t_turn_tf
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "t_turn_tf_is_input",
        tfaveragedturngeomparam.t_turn_tf_is_input,
    )

    monkeypatch.setattr(
        tfcoil_variables, "c_tf_turn", tfaveragedturngeomparam.c_tf_turn
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_cable_tf", tfaveragedturngeomparam.t_cable_tf
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "t_cable_tf_is_input",
        tfaveragedturngeomparam.t_cable_tf_is_input,
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "a_tf_wp_no_insulation",
        tfaveragedturngeomparam.a_tf_wp_no_insulation,
    )

    monkeypatch.setattr(
        sctfcoil_module, "dr_tf_turn", tfaveragedturngeomparam.dr_tf_turn
    )

    monkeypatch.setattr(
        sctfcoil_module, "dx_tf_turn", tfaveragedturngeomparam.dx_tf_turn
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "dx_tf_turn_cable_space_average",
        tfaveragedturngeomparam.dx_tf_turn_cable_space_average,
    )

    (
        a_tf_turn_cable_space_no_void,
        a_tf_turn_steel,
        a_tf_turn_insulation,
        n_tf_coil_turns,
    ) = sctfcoil.tf_averaged_turn_geom(
        i_tf_sc_mat=tfaveragedturngeomparam.i_tf_sc_mat,
        j_tf_wp=tfaveragedturngeomparam.j_tf_wp,
        dx_tf_turn_steel=tfaveragedturngeomparam.dx_tf_turn_steel,
        dx_tf_turn_insulation=tfaveragedturngeomparam.dx_tf_turn_insulation,
    )

    assert tfcoil_variables.t_conductor == pytest.approx(
        tfaveragedturngeomparam.expected_t_conductor
    )

    assert tfcoil_variables.t_turn_tf == pytest.approx(
        tfaveragedturngeomparam.expected_t_turn_tf
    )

    assert sctfcoil_module.dr_tf_turn == pytest.approx(
        tfaveragedturngeomparam.expected_t_turn_radial
    )

    assert sctfcoil_module.dx_tf_turn == pytest.approx(
        tfaveragedturngeomparam.expected_dx_tf_turn
    )

    assert sctfcoil_module.dx_tf_turn_cable_space_average == pytest.approx(
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
        sctfcoil_module,
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

    monkeypatch.setattr(sctfcoil_module, "a_tf_coil_inboard_steel", 0.55)  # Section 3

    # Sum from Section 3
    monkeypatch.setattr(sctfcoil_module, "a_tf_plasma_case", 0.42)
    monkeypatch.setattr(sctfcoil_module, "a_tf_coil_nose_case", 0.42)
    monkeypatch.setattr(sctfcoil_module, "dx_tf_side_case_average", 0.05)

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
    monkeypatch.setattr(tfcoil_variables, "tdmptf", 30)  # Figure 6
    monkeypatch.setattr(sctfcoil_module, "c_tf_coil", 83200 * 192)  # Section 3

    monkeypatch.setattr(
        build_variables, "r_vv_inboard_out", 4.45 + (build_variables.dr_vv_inboard / 2)
    )  # Table 2

    sctfcoil.vv_stress_on_quench()

    assert pytest.approx(sctfcoil_module.vv_stress_quench) == 56893800.120420754
