import pytest

from process.models.physics import confinement_time as conf


@pytest.mark.parametrize(
    "func, args, expected",
    [
        (conf.neo_alcator_confinement_time, (1.0, 1.0, 1.0, 1.0), 0.07),
        (conf.mirnov_confinement_time, (1.0, 1.0, 1.0), 0.2),
        (
            conf.merezhkin_muhkovatov_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.011067971810589328,
        ),
        (conf.shimomura_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0), 0.045),
        (
            conf.kaye_goldston_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.04490731195102493,
        ),
        (
            conf.iter_89p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.048,
        ),
        (
            conf.iter_89_0_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.104,
        ),
        (
            conf.rebut_lallia_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.18434273785533292,
        ),
        (
            conf.goldston_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.03021037349432586,
        ),
        (
            conf.t10_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.07307692307692307,
        ),
        (
            conf.jaeri_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.13187536611732242,
        ),
        (
            conf.kaye_big_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.105,
        ),
        (
            conf.iter_h90_p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.064,
        ),
        (conf.riedel_l_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.044),
        (
            conf.christiansen_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.24,
        ),
        (
            conf.lackner_gottardi_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.13120344887319335,
        ),
        (conf.neo_kaye_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.063),
        (conf.riedel_h_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.1),
        (
            conf.iter_h90_p_amended_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.082,
        ),
        (conf.sudo_et_al_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0), 0.17),
        (conf.gyro_reduced_bohm_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0), 0.25),
        (
            conf.lackner_gottardi_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.17,
        ),
        (
            conf.iter_93h_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.036,
        ),
        (
            conf.iter_h97p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.031,
        ),
        (
            conf.iter_h97p_elmy_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.029,
        ),
        (
            conf.iter_96p_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.023,
        ),
        (
            conf.valovic_elmy_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.067,
        ),
        (conf.kaye_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.021),
        (
            conf.iter_pb98py_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0615,
        ),
        (
            conf.iter_ipb98y_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0365,
        ),
        (
            conf.iter_ipb98y1_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0503,
        ),
        (
            conf.iter_ipb98y2_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0562,
        ),
        (
            conf.iter_ipb98y3_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0564,
        ),
        (
            conf.iter_ipb98y4_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.0587,
        ),
        (
            conf.iss95_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.079,
        ),
        (
            conf.iss04_stellarator_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.134,
        ),
        (conf.ds03_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.028),
        (
            conf.murari_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.03669697337069055,
        ),
        (conf.petty08_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 0.052),
        (
            conf.lang_high_density_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            1.0970434976417428e-103,
        ),
        (conf.hubbard_nominal_confinement_time, (1.0, 1.0, 1.0, 1.0), 0.014),
        (conf.hubbard_lower_confinement_time, (1.0, 1.0, 1.0, 1.0), 0.014),
        (conf.hubbard_upper_confinement_time, (1.0, 1.0, 1.0, 1.0), 0.014),
        (
            conf.menard_nstx_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.095,
        ),
        (
            conf.menard_nstx_petty08_hybrid_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.095,
        ),
        (conf.nstx_gyro_bohm_confinement_time, (1.0, 1.0, 1.0, 1.0, 1.0), 0.21),
        (
            conf.itpa20_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.06802157257083392,
        ),
        (
            conf.itpa20_il_confinement_time,
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            0.09877603755850378,
        ),
    ],
)
def test_confinement_time(func, args, expected):
    result = func(*args)
    assert result == pytest.approx(expected)
