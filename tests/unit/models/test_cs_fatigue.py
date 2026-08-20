from typing import NamedTuple

import pytest

from process.models.engineering.materials import (
    PARIS_COEFFICIENT_SS_316LN,
)


@pytest.fixture
def cs_fatigue_python(process_models):
    """Fixture to create a CSFatigue object.

    :return: an instance of CSFatigue
    :rtype: process.cs_fatigue.CSFatigue
    """
    return process_models.cs_fatigue


class NcycleParam(NamedTuple):
    max_hoop_stress: float
    residual_stress: float
    dz_cs_turn_crack_initial: float
    dz_cs_turn_conduit: float
    dr_cs_turn_conduit: float
    dr_cs_turn_crack_initial: float
    paris_coefficient_cs_turn: float
    expected_n_cycle: float
    expected_t_crack_radial: float


@pytest.mark.parametrize(
    "ncycleparam",
    [
        NcycleParam(
            max_hoop_stress=659999225.25370133,
            residual_stress=240000000,
            dz_cs_turn_crack_initial=0.00088999999999999995,
            dz_cs_turn_conduit=0.0063104538380405924,
            dr_cs_turn_conduit=0.0063104538380405924,
            dr_cs_turn_crack_initial=0.0026699999999999996,
            paris_coefficient_cs_turn=PARIS_COEFFICIENT_SS_316LN,
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

    n_cycle, dr_cs_turn_crack_initial = cs_fatigue_python.ncycle(
        max_hoop_stress=ncycleparam.max_hoop_stress,
        residual_stress=ncycleparam.residual_stress,
        dz_cs_turn_crack_initial=ncycleparam.dz_cs_turn_crack_initial,
        dz_cs_turn_conduit=ncycleparam.dz_cs_turn_conduit,
        dr_cs_turn_conduit=ncycleparam.dr_cs_turn_conduit,
        paris_coefficient_cs_turn=ncycleparam.paris_coefficient_cs_turn,
    )

    assert n_cycle == pytest.approx(ncycleparam.expected_n_cycle)

    assert dr_cs_turn_crack_initial == pytest.approx(ncycleparam.expected_t_crack_radial)


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
        )
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
