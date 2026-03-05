import pytest

from process.models.physics.confinement_time import PlasmaConfinementTime


@pytest.mark.parametrize(
    "func, args, expected",
    [
        (
            PlasmaConfinementTime().neo_alcator_confinement_time,
            (1.0, 1.0, 1.0, 1.0),
            0.07,
        ),
        (PlasmaConfinementTime().mirnov_confinement_time, (1.0, 1.0, 1.0), 0.2),
        (
            PlasmaConfinementTime().merezhkin_muhkovatov_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.011067971810589328,
        ),
        (
            PlasmaConfinementTime().shimomura_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0),
            0.045,
        ),
        (
            PlasmaConfinementTime().kaye_goldston_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.04490731195102493,
        ),
        (
            PlasmaConfinementTime().iter_89p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.048,
        ),
        (
            PlasmaConfinementTime().iter_89_0_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.104,
        ),
        (
            PlasmaConfinementTime().rebut_lallia_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.18434273785533292,
        ),
        (
            PlasmaConfinementTime().goldston_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.03021037349432586,
        ),
        (
            PlasmaConfinementTime().t10_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.07307692307692307,
        ),
        (
            PlasmaConfinementTime().jaeri_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.13187536611732242,
        ),
        (
            PlasmaConfinementTime().kaye_big_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.105,
        ),
        (
            PlasmaConfinementTime().iter_h90_p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.064,
        ),
        (
            PlasmaConfinementTime().riedel_l_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.044,
        ),
        (
            PlasmaConfinementTime().christiansen_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.24,
        ),
        (
            PlasmaConfinementTime().lackner_gottardi_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.13120344887319335,
        ),
        (
            PlasmaConfinementTime().neo_kaye_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.063,
        ),
        (
            PlasmaConfinementTime().riedel_h_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.1,
        ),
        (
            PlasmaConfinementTime().iter_h90_p_amended_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.082,
        ),
        (
            PlasmaConfinementTime().sudo_et_al_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0),
            0.17,
        ),
        (
            PlasmaConfinementTime().gyro_reduced_bohm_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0),
            0.25,
        ),
        (
            PlasmaConfinementTime().lackner_gottardi_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.17,
        ),
        (
            PlasmaConfinementTime().iter_93h_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.036,
        ),
        (
            PlasmaConfinementTime().iter_h97p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.031,
        ),
        (
            PlasmaConfinementTime().iter_h97p_elmy_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.029,
        ),
        (
            PlasmaConfinementTime().iter_96p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.023,
        ),
        (
            PlasmaConfinementTime().valovic_elmy_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.067,
        ),
        (
            PlasmaConfinementTime().kaye_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.021,
        ),
        (
            PlasmaConfinementTime().iter_pb98py_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0615,
        ),
        (
            PlasmaConfinementTime().iter_ipb98y_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0365,
        ),
        (
            PlasmaConfinementTime().iter_ipb98y1_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0503,
        ),
        (
            PlasmaConfinementTime().iter_ipb98y2_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0562,
        ),
        (
            PlasmaConfinementTime().iter_ipb98y3_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0564,
        ),
        (
            PlasmaConfinementTime().iter_ipb98y4_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0587,
        ),
        (
            PlasmaConfinementTime().iss95_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.079,
        ),
        (
            PlasmaConfinementTime().iss04_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.134,
        ),
        (
            PlasmaConfinementTime().ds03_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.028,
        ),
        (
            PlasmaConfinementTime().murari_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.03669697337069055,
        ),
        (
            PlasmaConfinementTime().petty08_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.052,
        ),
        (
            PlasmaConfinementTime().lang_high_density_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            1.0970434976417428e-103,
        ),
        (
            PlasmaConfinementTime().hubbard_nominal_confinement_time,
            (1.0, 1.0, 1.0, 1.0),
            0.014,
        ),
        (
            PlasmaConfinementTime().hubbard_lower_confinement_time,
            (1.0, 1.0, 1.0, 1.0),
            0.014,
        ),
        (
            PlasmaConfinementTime().hubbard_upper_confinement_time,
            (1.0, 1.0, 1.0, 1.0),
            0.014,
        ),
        (
            PlasmaConfinementTime().menard_nstx_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.095,
        ),
        (
            PlasmaConfinementTime().menard_nstx_petty08_hybrid_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.095,
        ),
        (
            PlasmaConfinementTime().nstx_gyro_bohm_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0),
            0.21,
        ),
        (
            PlasmaConfinementTime().itpa20_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.06802157257083392,
        ),
        (
            PlasmaConfinementTime().itpa20_il_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.09877603755850378,
        ),
    ],
)
def test_confinement_time(func, args, expected):
    result = func(*args)
    assert result == pytest.approx(expected)
