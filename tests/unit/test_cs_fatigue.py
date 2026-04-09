from typing import Any, NamedTuple

import pytest

from process.models.cs_fatigue import CsFatigue


@pytest.fixture
def cs_fatigue_python():
    """Fixture to create a CsFatigue object.

    :return: an instance of CsFatigue
    :rtype: process.cs_fatigue.CsFatigue
    """
    return CsFatigue()


class NcycleParam(NamedTuple):
    max_hoop_stress: Any = None

    residual_stress: Any = None

    t_crack_vertical: Any = None

    dz_cs_turn_conduit: Any = None

    dr_cs_turn_conduit: Any = None

    n_cycle: Any = None

    t_crack_radial: Any = None

    expected_n_cycle: Any = None

    expected_t_crack_radial: Any = None


@pytest.mark.parametrize(
    "ncycleparam",
    (
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
    ),
)
def test_ncycle(ncycleparam, monkeypatch, cs_fatigue_python):
    """
    Automatically generated Regression Unit Test for ncycle.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

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
    "hoop_stress, t, w, a, c, phi, expected_k",
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
    "hoop_stress, t, w, a, c, phi, expected_k",
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
