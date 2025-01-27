"""Unit tests for physics_functions.f90."""

from typing import Any, NamedTuple

import numpy as np
import pytest
from pytest import approx

from process import physics_functions
from process.fortran import physics_variables as pv


class SetFusionPowersParam(NamedTuple):
    f_alpha_plasma: Any = None

    f_deuterium: Any = None

    i_beta_fast_alpha: Any = None

    bp: Any = None

    bt: Any = None

    dene: Any = None

    deni: Any = None

    dnitot: Any = None

    f_alpha_electron: Any = None

    f_alpha_ion: Any = None

    alpha_power_beams: Any = None

    charged_power_density: Any = None

    ten: Any = None

    tin: Any = None

    vol_plasma: Any = None

    alpha_power_density_plasma: Any = None

    neutron_power_density_plasma: Any = None

    expected_alpha_power_density: Any = None

    expected_neutron_power_density: Any = None

    expected_alpha_power_total: Any = None

    expected_neutron_power_total: Any = None

    expected_non_alpha_charged_power: Any = None

    expected_beta_fast_alpha: Any = None

    expected_alpha_power_electron_density: Any = None

    expected_alpha_power_ion_density: Any = None

    expected_charged_particle_power: Any = None

    expected_fusion_power: Any = None


@pytest.mark.parametrize(
    "setfusionpowersparam",
    (
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            f_deuterium=0.5,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            alpha_power_beams=0,
            charged_power_density=0.00066,
            vol_plasma=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.163,
            expected_neutron_power_density=0.654,
            expected_alpha_power_total=395.47875,
            expected_neutron_power_total=1586.7675,
            expected_non_alpha_charged_power=1.601325,
            expected_alpha_power_ion_density=0.049552,
            expected_alpha_power_electron_density=0.105298,
            expected_charged_particle_power=397.080075,
            expected_fusion_power=1983.847575,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            f_deuterium=0.5,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            vol_plasma=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.8183263050336705,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1985.464197587943,
            expected_non_alpha_charged_power=1.601325,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2483.04427258794345,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            f_deuterium=0.5,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            vol_plasma=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.8183263050336705,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1985.464197587943,
            expected_non_alpha_charged_power=1.601325,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2483.0442725879434,
        ),
        SetFusionPowersParam(
            f_alpha_plasma=0.95,
            f_deuterium=2.5,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            alpha_power_beams=100.5,
            charged_power_density=0.00066,
            vol_plasma=2426.25,
            alpha_power_density_plasma=0.163,
            neutron_power_density_plasma=0.654,
            expected_alpha_power_density=0.20442195,
            expected_neutron_power_density=0.8183263050336705,
            expected_alpha_power_total=495.97875,
            expected_neutron_power_total=1985.464197587943,
            expected_non_alpha_charged_power=1.601325,
            expected_alpha_power_ion_density=0.062144272,
            expected_alpha_power_electron_density=0.132056578,
            expected_charged_particle_power=497.580075,
            expected_fusion_power=2483.0442725879434,
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
    monkeypatch.setattr(pv, "f_deuterium", setfusionpowersparam.f_deuterium)

    (
        neutron_power_density_total,
        alpha_power_plasma,
        alpha_power_total,
        neutron_power_plasma,
        neutron_power_total,
        non_alpha_charged_power,
        alpha_power_density_total,
        alpha_power_electron_density,
        alpha_power_ions_density,
        charged_particle_power,
        fusion_power,
    ) = physics_functions.set_fusion_powers(
        f_alpha_electron=setfusionpowersparam.f_alpha_electron,
        f_alpha_ion=setfusionpowersparam.f_alpha_ion,
        alpha_power_beams=setfusionpowersparam.alpha_power_beams,
        charged_power_density=setfusionpowersparam.charged_power_density,
        neutron_power_density_plasma=setfusionpowersparam.neutron_power_density_plasma,
        vol_plasma=setfusionpowersparam.vol_plasma,
        alpha_power_density_plasma=setfusionpowersparam.alpha_power_density_plasma,
    )

    assert alpha_power_density_total == pytest.approx(
        setfusionpowersparam.expected_alpha_power_density
    )
    assert neutron_power_density_total == pytest.approx(
        setfusionpowersparam.expected_neutron_power_density
    )
    assert alpha_power_total == pytest.approx(
        setfusionpowersparam.expected_alpha_power_total
    )
    assert neutron_power_total == pytest.approx(
        setfusionpowersparam.expected_neutron_power_total
    )
    assert non_alpha_charged_power == pytest.approx(
        setfusionpowersparam.expected_non_alpha_charged_power
    )
    assert alpha_power_electron_density == pytest.approx(
        setfusionpowersparam.expected_alpha_power_electron_density
    )
    assert alpha_power_ions_density == pytest.approx(
        setfusionpowersparam.expected_alpha_power_ion_density
    )
    assert charged_particle_power == pytest.approx(
        setfusionpowersparam.expected_charged_particle_power
    )
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
    beta_beam, beam_density_out, alpha_power_beams = physics_functions.beam_fusion(
        1.0,
        1.5,
        0.85,
        5.3,
        130,
        7.8e19,
        6.6e19,
        17.8,
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

    assert beta_beam == pytest.approx(0.0026264022466211366)
    assert beam_density_out == pytest.approx(4.2133504058678246e17)
    assert alpha_power_beams == pytest.approx(11.593221085189192)


def test_beamcalc():
    (
        deuterium_beam_alpha_power,
        tritium_beam_alpha_power,
        hot_beam_density,
        beam_deposited_energy,
    ) = physics_functions.beamcalc(
        3.3e19,
        3.3e19,
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

    assert deuterium_beam_alpha_power == pytest.approx(11.561197668655383)
    assert tritium_beam_alpha_power == pytest.approx(1.0434445041616093e-05)
    assert hot_beam_density == pytest.approx(4.1968331737565126e17)
    assert beam_deposited_energy == pytest.approx(445.05787301616635)


def test__fast_ion_pressure_integral():
    pressure_integral = physics_functions.fast_ion_pressure_integral(1000.0, 276.7)

    assert pressure_integral == pytest.approx(1.1061397270783706)


def test_alpha_power_beam():
    alpha_power_beam = physics_functions.alpha_power_beam(
        316000000000, 3.3e19, 7.5e-22, 1888.0, 13.5, 2.8e-22
    )

    assert alpha_power_beam == pytest.approx(1.047572705194316e-05)


def test_beam_reaction_rate():
    beam_reaction_rate = physics_functions.beam_reaction_rate(
        3.01550071597, 5140000.0, 1000.0
    )

    assert beam_reaction_rate == pytest.approx(7.465047902975452e-18)
