"""Unit tests for fusion_reactions.py."""

from typing import Any, NamedTuple

import numpy as np
import pytest

from process.core import constants
from process.core.data_structure.base import DataStructure
from process.data_structure.physics_variables import PlasmaIgnitionModel
from process.models.physics import fusion_reactions as reactions


class SetFusionPowersParam(NamedTuple):
    f_p_alpha_plasma_deposited: Any = None

    f_alpha_electron: Any = None

    f_alpha_ion: Any = None

    p_beam_alpha_mw: Any = None

    pden_non_alpha_charged_mw: Any = None

    vol_plasma: Any = None

    pden_plasma_alpha_vol_avg_mw: Any = None

    pden_plasma_neutron_vol_avg_mw: Any = None

    expected_alpha_power_density: Any = None

    expected_neutron_power_density: Any = None

    expected_alpha_power_total: Any = None

    expected_neutron_power_total: Any = None

    expected_non_alpha_charged_power: Any = None

    expected_alpha_power_electron_density: Any = None

    expected_alpha_power_ion_density: Any = None

    expected_charged_particle_power: Any = None

    expected_fusion_power: Any = None


@pytest.mark.parametrize(
    "setfusionpowersparam",
    [
        SetFusionPowersParam(
            f_p_alpha_plasma_deposited=0.95,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            p_beam_alpha_mw=0,
            pden_non_alpha_charged_mw=0.00066,
            vol_plasma=2426.25,
            pden_plasma_alpha_vol_avg_mw=0.163,
            pden_plasma_neutron_vol_avg_mw=0.654,
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
            f_p_alpha_plasma_deposited=0.95,
            f_alpha_electron=0.68,
            f_alpha_ion=0.32,
            p_beam_alpha_mw=100.5,
            pden_non_alpha_charged_mw=0.00066,
            vol_plasma=2426.25,
            pden_plasma_alpha_vol_avg_mw=0.163,
            pden_plasma_neutron_vol_avg_mw=0.654,
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
    ],
)
def test_set_fusion_powers(setfusionpowersparam):
    """
    Automatically generated Unit Test for set_fusion_powers().

    This test was generated using data from baseline_2018_IN.DAT
    (no longer exists in the PROCESS repository).

    :param setfusionpowersparam: the data used to mock and assert in this test.
    :type setfusionpowersparam: setfusionpowersparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        pden_neutron_total_vol_avg_mw,
        _p_plasma_alpha_mw,
        p_alpha_total_mw,
        _p_plasma_neutron_mw,
        p_neutron_total_mw,
        p_non_alpha_charged_mw,
        pden_alpha_total_vol_avg_mw,
        f_pden_alpha_electron_mw,
        f_pden_alpha_ions_mw,
        p_charged_particle_mw,
        p_fusion_total_mw,
    ) = reactions.set_fusion_powers(
        f_alpha_electron=setfusionpowersparam.f_alpha_electron,
        f_alpha_ion=setfusionpowersparam.f_alpha_ion,
        p_beam_alpha_mw=setfusionpowersparam.p_beam_alpha_mw,
        pden_non_alpha_charged_mw=setfusionpowersparam.pden_non_alpha_charged_mw,
        pden_plasma_neutron_vol_avg_mw=setfusionpowersparam.pden_plasma_neutron_vol_avg_mw,
        vol_plasma=setfusionpowersparam.vol_plasma,
        pden_plasma_alpha_vol_avg_mw=setfusionpowersparam.pden_plasma_alpha_vol_avg_mw,
        f_p_alpha_plasma_deposited=setfusionpowersparam.f_p_alpha_plasma_deposited,
    )

    assert pden_alpha_total_vol_avg_mw == pytest.approx(
        setfusionpowersparam.expected_alpha_power_density
    )
    assert pden_neutron_total_vol_avg_mw == pytest.approx(
        setfusionpowersparam.expected_neutron_power_density
    )
    assert p_alpha_total_mw == pytest.approx(
        setfusionpowersparam.expected_alpha_power_total
    )
    assert p_neutron_total_mw == pytest.approx(
        setfusionpowersparam.expected_neutron_power_total
    )
    assert p_non_alpha_charged_mw == pytest.approx(
        setfusionpowersparam.expected_non_alpha_charged_power
    )
    assert f_pden_alpha_electron_mw == pytest.approx(
        setfusionpowersparam.expected_alpha_power_electron_density
    )
    assert f_pden_alpha_ions_mw == pytest.approx(
        setfusionpowersparam.expected_alpha_power_ion_density
    )
    assert p_charged_particle_mw == pytest.approx(
        setfusionpowersparam.expected_charged_particle_power
    )
    assert p_fusion_total_mw == pytest.approx(setfusionpowersparam.expected_fusion_power)


@pytest.mark.parametrize(
    ("t", "reaction", "expected_bosch_hale"),
    [
        (55.73, reactions.REACTION_CONSTANTS_DT, 8.832857074192583e-22),
        (55.73, reactions.REACTION_CONSTANTS_DHE3, 7.067916724597656e-23),
        (55.73, reactions.REACTION_CONSTANTS_DD1, 1.3127277533210717e-23),
        (55.73, reactions.REACTION_CONSTANTS_DD2, 1.1329338540436287e-23),
    ],
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
    bosch_hale = reactions.bosch_hale_reactivity(
        np.array([t]), reactions.BoschHaleConstants(**reaction)
    )

    assert bosch_hale == pytest.approx(expected_bosch_hale, abs=1e-23)


def test_beam_fusion():
    beta_beam, nd_beam_ions_out, p_beam_alpha_mw = reactions.beam_fusion(
        1.0,
        1.5,
        5.367727,
        130,
        7.8e19,
        6.6e19,
        17.8,
        1000.0,
        0.5,
        0.5,
        1e-06,
        13.5,
        1888.0,
        0.425,
    )

    assert beta_beam == pytest.approx(0.0026264022466211366)
    assert nd_beam_ions_out == pytest.approx(4.2133504058678246e17)
    assert p_beam_alpha_mw == pytest.approx(9.271206216041564)


def test_beam_slowing_down_state():
    beam_state = reactions.beam_slowing_down_state(
        1000.0,
        276.7,
        415.0,
        1.42,
        1e-06,
        130,
        1888.0,
    )

    assert beam_state.deuterium_beam_density == pytest.approx(4.1968331737565126e17)
    assert beam_state.tritium_beam_density == pytest.approx(316553077182.4059, rel=1e-6)
    assert beam_state.deuterium_critical_energy_speed == pytest.approx(
        5.1495e6, rel=1e-4
    )
    assert beam_state.tritium_critical_energy_speed == pytest.approx(5.1534e6, rel=1e-4)
    assert beam_state.nd_beam_hot == pytest.approx(4.1968331737565126e17)
    assert beam_state.e_beam_deposited_kev == pytest.approx(445.05787301616635)


def test__fast_ion_pressure_integral():
    pressure_integral = reactions.fast_ion_pressure_integral(1000.0, 276.7)

    assert pressure_integral == pytest.approx(1.1061397270783706)


def test_beam_target_reaction_rate():
    reaction_rate = reactions.beam_target_reaction_rate(
        nd_beam_ion=3.16e11,
        nd_target_ion=3.3e19,
        sigv_beam=7.5e-22,
        vol_plasma=1888.0,
    )

    assert reaction_rate == pytest.approx(1.4766048e13)


def test_alpha_power_beam():
    beam_target_rate = 1.0e13  # s^-1
    result = reactions.alpha_power_beam(beam_target_rate)
    expected = beam_target_rate * constants.DT_ALPHA_ENERGY / 1.0e6

    assert result == pytest.approx(expected)


def test_beam_reaction_rate_coefficient():
    beam_reaction_rate = reactions.beam_reaction_rate_coefficient(
        3.01550071597, 5140000.0, 1000.0
    )

    assert beam_reaction_rate == pytest.approx(7.465047902975452e-22)


def _make_beam_fusion_reactions():
    return reactions.BeamReactions(data=DataStructure())


def test_calculate_beam_fusion():
    beam_reactions = _make_beam_fusion_reactions()
    data = beam_reactions.data

    data.physics.beamfus0 = 1.0
    data.physics.betbm0 = 1.5
    data.physics.b_plasma_total = 5.367727
    data.current_drive.c_beam_total = 130
    data.physics.nd_plasma_electrons_vol_avg = 7.8e19
    data.physics.nd_plasma_fuel_ions_vol_avg = 6.6e19
    data.physics.dlamie = 17.8
    data.current_drive.e_beam_kev = 1000.0
    data.physics.f_plasma_fuel_deuterium = 0.5
    data.physics.f_plasma_fuel_tritium = 0.5
    data.current_drive.f_beam_tritium = 1e-06
    data.physics.temp_plasma_electron_density_weighted_kev = 13.5
    data.physics.vol_plasma = 1888.0
    data.physics.n_charge_plasma_effective_mass_weighted_vol_avg = 0.425
    data.physics.i_plasma_ignited = PlasmaIgnitionModel.NON_IGNITED

    beam_reactions.calculate_beam_fusion()

    assert beam_reactions.beta_beam == pytest.approx(0.0026264022466211366)
    assert beam_reactions.nd_beam_ions_out == pytest.approx(4.2133504058678246e17)
    assert beam_reactions.p_beam_alpha_mw == pytest.approx(9.271206216041564)
    assert beam_reactions.p_beam_neutron_mw == pytest.approx(
        beam_reactions.p_beam_alpha_mw
        * (
            constants.DT_NEUTRON_ENERGY_FRACTION
            / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION)
        )
    )
    assert beam_reactions.p_beam_dt_mw == pytest.approx(
        beam_reactions.p_beam_alpha_mw / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION)
    )


@pytest.mark.parametrize(
    ("c_beam_total", "i_plasma_ignited"),
    [
        (0.0, PlasmaIgnitionModel.NON_IGNITED),
        (130.0, PlasmaIgnitionModel.IGNITED),
    ],
    ids=["no_beam_current", "ignited_plasma"],
)
def test_calculate_beam_fusion_neglected(c_beam_total, i_plasma_ignited):
    beam_reactions = _make_beam_fusion_reactions()
    beam_reactions.data.current_drive.c_beam_total = c_beam_total
    beam_reactions.data.physics.i_plasma_ignited = i_plasma_ignited

    beam_reactions.calculate_beam_fusion()

    assert beam_reactions.beta_beam == pytest.approx(0.0)
    assert beam_reactions.nd_beam_ions_out == pytest.approx(0.0)
    assert beam_reactions.p_beam_alpha_mw == pytest.approx(0.0)
    assert beam_reactions.p_beam_neutron_mw == pytest.approx(0.0)
    assert beam_reactions.p_beam_dt_mw == pytest.approx(0.0)


def test_calculate_beam_fusion_resets_stale_state():
    """A BeamReactions instance is composed once and reused across many
    calls, so beam results from a previous call must not persist once the beam
    is switched off.
    """
    beam_reactions = _make_beam_fusion_reactions()
    data = beam_reactions.data

    data.physics.beamfus0 = 1.0
    data.physics.betbm0 = 1.5
    data.physics.b_plasma_total = 5.367727
    data.current_drive.c_beam_total = 130
    data.physics.nd_plasma_electrons_vol_avg = 7.8e19
    data.physics.nd_plasma_fuel_ions_vol_avg = 6.6e19
    data.physics.dlamie = 17.8
    data.current_drive.e_beam_kev = 1000.0
    data.physics.f_plasma_fuel_deuterium = 0.5
    data.physics.f_plasma_fuel_tritium = 0.5
    data.current_drive.f_beam_tritium = 1e-06
    data.physics.temp_plasma_electron_density_weighted_kev = 13.5
    data.physics.vol_plasma = 1888.0
    data.physics.n_charge_plasma_effective_mass_weighted_vol_avg = 0.425
    data.physics.i_plasma_ignited = PlasmaIgnitionModel.NON_IGNITED

    beam_reactions.calculate_beam_fusion()
    assert beam_reactions.p_beam_alpha_mw != pytest.approx(0.0)

    data.current_drive.c_beam_total = 0.0
    beam_reactions.calculate_beam_fusion()

    assert beam_reactions.beta_beam == pytest.approx(0.0)
    assert beam_reactions.nd_beam_ions_out == pytest.approx(0.0)
    assert beam_reactions.p_beam_alpha_mw == pytest.approx(0.0)
    assert beam_reactions.p_beam_neutron_mw == pytest.approx(0.0)
    assert beam_reactions.p_beam_dt_mw == pytest.approx(0.0)


def test_beam_fusion_reactions_set_physics_variables():
    beam_reactions = _make_beam_fusion_reactions()
    beam_reactions.beta_beam = 1.1
    beam_reactions.nd_beam_ions_out = 2.2
    beam_reactions.p_beam_alpha_mw = 3.3
    beam_reactions.p_beam_neutron_mw = 4.4
    beam_reactions.p_beam_dt_mw = 5.5

    beam_reactions.set_physics_variables()

    data = beam_reactions.data
    assert data.physics.beta_beam == pytest.approx(1.1)
    assert data.physics.nd_beam_ions_out == pytest.approx(2.2)
    assert data.physics.p_beam_alpha_mw == pytest.approx(3.3)
    assert data.physics.p_beam_neutron_mw == pytest.approx(4.4)
    assert data.physics.p_beam_dt_mw == pytest.approx(5.5)
