"""Unit tests for physics_functions.f90."""


from typing import Any, NamedTuple
from process.fortran import physics_variables as pv
from process import physics_functions
import numpy as np
from process import physics_functions
import numpy as np
import pytest
from pytest import approx


class SetFusionPowersParam(NamedTuple):
    f_alpha_plasma: Any = None

    fdeut: Any = None

    ifalphap: Any = None

    bp: Any = None

    bt: Any = None

    dene: Any = None

    deni: Any = None

    dnitot: Any = None

    falpe: Any = None

    falpi: Any = None

    alpha_power_beams: Any = None

    charged_power_density: Any = None

    ten: Any = None

    tin: Any = None

    plasma_volume: Any = None

    alpha_power_density_plasma: Any = None

    neutron_power_density_plasma: Any = None

    expected_alpha_power_density: Any = None

    expected_neutron_power_density: Any = None

    expected_alpha_power_total: Any = None

    expected_neutron_power_total: Any = None

    expected_non_alpha_charged_power: Any = None

    expected_betaft: Any = None

    expected_alpha_power_electron_density: Any = None

    expected_alpha_power_ion_density: Any = None

    expected_charged_particle_power: Any = None

    expected_fusion_power: Any = None


@pytest.mark.parametrize(
    "setfusionpowersparam",
    (
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            fdeut=0.5,
            ifalphap=1,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            alpha_power_beams=0,
            charged_power_density=0.00066,
            ten=13.84,
            tin=13.84,
            plasma_volume=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.163,
            expected_neutron_power_density=0.654,
            expected_alpha_power_total=395.47875,
            expected_neutron_power_total=1586.7675,
            expected_non_alpha_charged_power=1.601325,
            expected_betaft=0.00423788,
            expected_alpha_power_ion_density=0.049552,
            expected_alpha_power_electron_density=0.105298,
            expected_charged_particle_power=397.080075,
            expected_fusion_power=1983.847575,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            fdeut=0.5,
            ifalphap=1,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            ten=13.84,
            tin=13.84,
            plasma_volume=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.81968779,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1988.7675,
            expected_non_alpha_charged_power=1.601325,
            expected_betaft=0.00531482,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2486.347575,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            fdeut=0.5,
            ifalphap=0,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            ten=13.84,
            tin=13.84,
            plasma_volume=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.81968779,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1988.7675,
            expected_non_alpha_charged_power=1.601325,
            expected_betaft=0.00701622,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2486.347575,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            fdeut=2.5,
            ifalphap=0,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            ten=13.84,
            tin=13.84,
            plasma_volume=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.81968779,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1988.7675,
            expected_non_alpha_charged_power=1.601325,
            expected_betaft=0.0,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2486.347575,
        ),
    ),
)
def test_set_fusion_powers(setfusionpowersparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for set_fusion_powers().

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param setfusionpowersparam: the data used to mock and assert in this test.
    :type setfusionpowersparam: setfusionpowersparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(pv, "f_alpha_plasma", setfusionpowersparam.f_alpha_plasma)
    monkeypatch.setattr(pv, "fdeut", setfusionpowersparam.fdeut)

    (
        neutron_power_density,
        alpha_power_total,
        neutron_power_total,
        non_alpha_charged_power,
        betaft,
        alpha_power_density,
        alpha_power_electron_density,
        alpha_power_ions_density,
        charged_particle_power,
        fusion_power,
    ) = physics_functions.set_fusion_powers(
        ifalphap=setfusionpowersparam.ifalphap,
        bp=setfusionpowersparam.bp,
        bt=setfusionpowersparam.bt,
        dene=setfusionpowersparam.dene,
        deni=setfusionpowersparam.deni,
        dnitot=setfusionpowersparam.dnitot,
        falpe=setfusionpowersparam.falpe,
        falpi=setfusionpowersparam.falpi,
        alpha_power_beams=setfusionpowersparam.alpha_power_beams,
        charged_power_density=setfusionpowersparam.charged_power_density,
        ten=setfusionpowersparam.ten,
        tin=setfusionpowersparam.tin,
        plasma_volume=setfusionpowersparam.plasma_volume,
        alpha_power_density_plasma=setfusionpowersparam.alpha_power_density_plasma,
        neutron_power_density_plasma=setfusionpowersparam.neutron_power_density_plasma,
    )

    assert alpha_power_density == pytest.approx(setfusionpowersparam.expected_alpha_power_density)
    assert neutron_power_density == pytest.approx(setfusionpowersparam.expected_neutron_power_density)
    assert alpha_power_total == pytest.approx(setfusionpowersparam.expected_alpha_power_total)
    assert neutron_power_total == pytest.approx(setfusionpowersparam.expected_neutron_power_total)
    assert non_alpha_charged_power == pytest.approx(setfusionpowersparam.expected_non_alpha_charged_power)
    assert betaft == pytest.approx(setfusionpowersparam.expected_betaft)
    assert alpha_power_electron_density == pytest.approx(setfusionpowersparam.expected_alpha_power_electron_density)
    assert alpha_power_ions_density == pytest.approx(setfusionpowersparam.expected_alpha_power_ion_density)
    assert charged_particle_power== pytest.approx(setfusionpowersparam.expected_charged_particle_power)
    assert fusion_power == pytest.approx(setfusionpowersparam.expected_fusion_power)


@pytest.mark.parametrize(
    "t, reaction, expected_bosch_hale",
    (
        (55.73, physics_functions.REACTION_CONSTANTS_DT, 8.832857074192583e-22),
        (55.73, physics_functions.REACTION_CONSTANTS_DHE3, 7.067916724597656e-23),
        (55.73, physics_functions.REACTION_CONSTANTS_DD1, 1.3127277533210717e-23),
        (55.73, physics_functions.REACTION_CONSTANTS_DD2, 1.1329338540436287e-23),
    ),
    ids=["DT", "DHE3", "DD1", "DD2"],
)
def test_bosch_hale(t, reaction, expected_bosch_hale):
    """
    Unit test for the bosch_hale function.

    :param t: input Maxwellian density-weighted ion temperature
    :type t: float
    :param reaction: input flag for fusion reaction to use
    :type reaction: int
    :param expected_bosch_hale: expected return value from the bosch_hale function
    :type expected_bosch_hale: float
    """
    bosch_hale = physics_functions.bosch_hale_reactivity(
        np.array([t]), physics_functions.BoschHaleConstants(**reaction)
    )

    assert bosch_hale == approx(expected_bosch_hale, abs=1e-23)


def test_beam_fusion():
    betanb, dnbeam2, alpha_power_beams = physics_functions.beam_fusion(
        1.0,
        1.5,
        0.85,
        5.3,
        130,
        7.8e19,
        6.6e19,
        17.8,
        3520.0,
        1000.0,
        0.5,
        0.5,
        1e-06,
        2.8e-22,
        13.5,
        13.5,
        1888.0,
        0.425,
    )

    assert betanb == pytest.approx(0.002616169278788316)
    assert dnbeam2 == pytest.approx(4.2028390908892986e17)
    assert alpha_power_beams == pytest.approx(11.506114015489336)


def test_beamcalc():
    palfdb, palftb, nhot, ehot = physics_functions.beamcalc(
        3.3e19,
        3.3e19,
        3520.0,
        1000.0,
        276.7,
        415.0,
        1.42,
        1e-06,
        130,
        13.5,
        1888.0,
        2.8e-22,
    )

    assert palfdb == pytest.approx(11.489365278680932)
    assert palftb == pytest.approx(1.0379265294979434e-05)
    assert nhot == pytest.approx(4.1968331737565126e17)
    assert ehot == pytest.approx(445.05787301616635)


def test__fast_ion_pressure_integral():
    xbrak = physics_functions._fast_ion_pressure_integral(1000.0, 276.7)

    assert xbrak == pytest.approx(1.1061397270783706)


def test_alpha_power_beam():
    palphabm = physics_functions.alpha_power_beam(
        3520.0, 316000000000, 3.3e19, 7.5e-22, 1888.0, 13.5, 2.8e-22
    )

    assert palphabm == pytest.approx(1.0413228502045627e-05)


def test_beam_reaction_rate():
    sgvhot = physics_functions.beam_reaction_rate(3, 5140000.0, 1000.0)

    assert sgvhot == pytest.approx(7.465047902975452e-18)
    assert palfdb == pytest.approx(11.489365278680932)
    assert palftb == pytest.approx(1.0379265294979434e-05)
    assert nhot == pytest.approx(4.1968331737565126e17)
    assert ehot == pytest.approx(445.05787301616635)


def test__fast_ion_pressure_integral():
    xbrak = physics_functions._fast_ion_pressure_integral(1000.0, 276.7)

    assert xbrak == pytest.approx(1.1061397270783706)


def test_alpha_power_beam():
    palphabm = physics_functions.alpha_power_beam(
        3520.0, 316000000000, 3.3e19, 7.5e-22, 1888.0, 13.5, 2.8e-22
    )

    assert palphabm == pytest.approx(1.0413228502045627e-05)


def test_beam_reaction_rate():
    sgvhot = physics_functions.beam_reaction_rate(3, 5140000.0, 1000.0)

    assert sgvhot == pytest.approx(7.465047902975452e-18)
