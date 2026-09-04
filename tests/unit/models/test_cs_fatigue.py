from typing import NamedTuple

import pytest


@pytest.fixture
def cs_fatigue_python(process_models):
    """Fixture to create a CsFatigue object.

    :return: an instance of CsFatigue
    :rtype: process.cs_fatigue.CsFatigue
    """
    return process_models.cs_fatigue


class NcycleParam(NamedTuple):
    max_hoop_stress: float
    residual_stress: float
    t_crack_vertical: float
    dz_cs_turn_conduit: float
    dr_cs_turn_conduit: float
    t_crack_radial: float
    expected_n_cycle: float
    expected_t_crack_radial: float


@pytest.mark.parametrize(
    "ncycleparam",
    [
        NcycleParam(
            max_hoop_stress=659999225.25370133,
            residual_stress=240000000,
            t_crack_vertical=0.00088999999999999995,
            dz_cs_turn_conduit=0.0063104538380405924,
            dr_cs_turn_conduit=0.0063104538380405924,
            t_crack_radial=0.0026699999999999996,
            expected_n_cycle=1113.5875631615095,
            expected_t_crack_radial=0.0026699999999999996,
        ),
    ],
)
def test_ncycle(ncycleparam, monkeypatch, cs_fatigue_python):
    """
    Automatically generated Unit Test for ncycle.

    This test was generated using data from baseline_2018_IN.DAT
    (no longer exists in the PROCESS repository).

    :param ncycleparam: the data used to mock and assert in this test.
    :type ncycleparam: ncycleparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    n_cycle, t_crack_radial = cs_fatigue_python.ncycle(
        max_hoop_stress=ncycleparam.max_hoop_stress,
        residual_stress=ncycleparam.residual_stress,
        t_crack_vertical=ncycleparam.t_crack_vertical,
        dz_cs_turn_conduit=ncycleparam.dz_cs_turn_conduit,
        dr_cs_turn_conduit=ncycleparam.dr_cs_turn_conduit,
    )

    assert n_cycle == pytest.approx(ncycleparam.expected_n_cycle)

    assert t_crack_radial == pytest.approx(ncycleparam.expected_t_crack_radial)


@pytest.mark.parametrize(
    ("hoop_stress", "t", "w", "a", "c", "phi", "expected_k"),
    [
        (
            659.99351867335338,
            0.0063104538380405924,
            0.0063104538380405924,
            0.00088999999999999995,
            0.0026699999999999996,
            1.5707963267948966,
            31.96412802853516,
        )
    ],
)
def test_embedded_stress_intensity_factor(
    hoop_stress, t, w, a, c, phi, expected_k, cs_fatigue_python
):
    """Tests `embedded_stress_intensity_factor` function.

    :param hoop_stress: change in hoop stress over cycle.
    :type hoop_stress: float

    :param t: plate thickness.
    :type t: float

    :param w: plate width.
    :type w: float

    :param a: crack depth (t -direction).
    :type a: float

    :param c: crack length (w - direction).
    :type c: float
    """
    k = cs_fatigue_python.embedded_stress_intensity_factor(hoop_stress, t, w, a, c, phi)

    assert pytest.approx(k) == expected_k


@pytest.mark.parametrize(
    ("hoop_stress", "t", "w", "a", "c", "phi", "expected_k"),
    [
        (
            659.99351867335338,
            0.0063104538380405924,
            0.0063104538380405924,
            0.00088999999999999995,
            0.0026699999999999996,
            1.5707963267948966,
            35.744426954844926,
        ),
        # a > c branch, phi = pi/2 and phi = 0. Expected values computed
        # independently from the Newman-Raju surface-crack equations
        # (NASA TM, "Stress-Intensity Factor Equations for Cracks in
        # Three-Dimensional Finite Bodies Subjected to Tension and Bending
        # Loads") for a/c > 1 with zero bending:
        # Q = 1 + 1.464(c/a)^1.65, M1 = sqrt(c/a)(1 + 0.04 c/a),
        # M2 = 0.2(c/a)^4, M3 = -0.11(c/a)^4,
        # g = 1 + [0.1 + 0.35 (c/a)(a/t)^2](1 - sin(phi))^2,
        # f_phi = [(c/a)^2 sin^2(phi) + cos^2(phi)]^0.25,
        # f_w = sec[pi c/(2w) sqrt(a/t)]^0.5,
        # K = sigma (M1 + M2 (a/t)^2 + M3 (a/t)^4) g f_phi f_w sqrt(pi a / Q)
        (
            659.99351867335338,
            0.0063104538380405924,
            0.0063104538380405924,
            0.0026699999999999996,
            0.00088999999999999995,
            1.5707963267948966,
            18.451605120658158,
        ),
        (
            659.99351867335338,
            0.0063104538380405924,
            0.0063104538380405924,
            0.0026699999999999996,
            0.00088999999999999995,
            0.0,
            35.82251644569079,
        ),
    ],
)
def test_surface_stress_intensity_factor(
    hoop_stress, t, w, a, c, phi, expected_k, cs_fatigue_python
):
    """Tests `surface_stress_intensity_factor` function.

    :param hoop_stress: change in hoop stress over cycle.
    :type hoop_stress: float

    :param t: plate thickness.
    :type t: float

    :param w: plate width.
    :type w: float

    :param a: crack depth (t -direction).
    :type a: float

    :param c: crack length (w - direction).
    :type c: float
    """
    k = cs_fatigue_python.surface_stress_intensity_factor(hoop_stress, t, w, a, c, phi)

    assert pytest.approx(k) == expected_k
