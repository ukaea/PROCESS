"""Unit tests for the p-B11 ion-electron equilibration power."""

import numpy as np
import pytest

from process.models.physics.physics import rether, rether_pb

# Uniform case: ne = 1e20 m^-3, Te = 10 keV, Ti = 11 keV, <Z^2/A> = 0.7.
# lnΛ = 31.3 - ln(ne)/2 + ln(Te in eV) = 17.484489442035727.
# P = 2.423340921690355e-41 * lnΛ * ne**2 * <Z^2/A> * (Ti - Te) / Te**1.5.
_UNIFORM_POWER_MW = 0.09379193834055811


def test_rether_pb_uniform_profile():
    """A uniform profile returns the local Spitzer power density."""
    rho = np.linspace(0.0, 1.0, 21)
    power = rether_pb(
        np.full_like(rho, 1.0e20),
        np.full_like(rho, 10.0),
        np.full_like(rho, 11.0),
        0.7,
        0,
        None,
        rho,
    )

    assert power == pytest.approx(_UNIFORM_POWER_MW)


def test_rether_pb_matches_rether_coefficient_ratio():
    """Flat p-B11 rate differs from rether only by the stored prefactor.

    rether uses 2.42165e-41. The prefactor from the physical constants is
    2.423340921690355e-41.
    """
    rho = np.linspace(0.0, 1.0, 21)
    ne = 1.0e20
    te = 10.0
    ti = 11.0
    charge_mass = 0.7
    dlamie = 31.3 - np.log(ne) / 2.0 + np.log(te * 1000.0)

    power = rether_pb(
        np.full_like(rho, ne),
        np.full_like(rho, te),
        np.full_like(rho, ti),
        charge_mass,
        0,
        None,
        rho,
    )
    reference = rether(0.0, 0.0, ne, dlamie, te, ti, charge_mass)

    assert power / reference == pytest.approx(2.423340921690355e-41 / 2.42165e-41)


def test_rether_pb_sign_follows_temperature_difference():
    """Ti < Te transfers power from electrons to ions."""
    rho = np.linspace(0.0, 1.0, 21)
    power = rether_pb(
        np.full_like(rho, 1.0e20),
        np.full_like(rho, 10.0),
        np.full_like(rho, 9.0),
        0.7,
        0,
        None,
        rho,
    )

    assert power == pytest.approx(-_UNIFORM_POWER_MW)


def test_rether_pb_nonpositive_profiles_are_zero():
    """No exchange is counted where density and electron temperature are zero."""
    rho = np.linspace(0.0, 1.0, 21)
    power = rether_pb(
        np.zeros_like(rho),
        np.zeros_like(rho),
        np.full_like(rho, 10.0),
        0.7,
        0,
        None,
        rho,
    )

    assert power == 0.0


def test_rether_pb_omits_nonpositive_separatrix():
    """A zero separatrix sample is left out of the 2 rho d(rho) average.

    Five radial points, the last with ne = Te = 0. Simpson quadrature of
    the remaining uniform power gives the expected average.
    """
    rho = np.linspace(0.0, 1.0, 5)
    ne = np.full_like(rho, 1.0e20)
    te = np.full_like(rho, 10.0)
    ti = np.full_like(rho, 11.0)
    ne[-1] = 0.0
    te[-1] = 0.0

    power = rether_pb(ne, te, ti, 0.7, 0, None, rho)

    assert power == pytest.approx(0.07815994861713176)


def test_rether_pb_floors_coulomb_logarithm():
    """A negative local Coulomb logarithm is replaced by 1.

    ne = 1e30 m^-3, Te = 0.01 keV, Ti = 0.02 keV, <Z^2/A> = 0.7.
    The unfloored logarithm is -0.9361913019166375.
    """
    rho = np.linspace(0.0, 1.0, 21)
    power = rether_pb(
        np.full_like(rho, 1.0e30),
        np.full_like(rho, 0.01),
        np.full_like(rho, 0.02),
        0.7,
        0,
        None,
        rho,
    )

    assert power == pytest.approx(1.6963386451832485e20)


def test_rether_pb_parabolic_profile():
    """Parabolic profiles use the local density squared and local Coulomb log.

    n_e = 2e20 (1-rho^2) m^-3, T_e = 40 (1-rho^2) keV,
    T_i = 80 (1-rho^2) keV, <Z^2/A> = 0.7, on 21 points from 0 to 1.
    The separatrix values are zero and are omitted.
    """
    rho = np.linspace(0.0, 1.0, 21)
    shape = 1.0 - rho**2
    power = rether_pb(
        2.0e20 * shape,
        40.0 * shape,
        80.0 * shape,
        0.7,
        0,
        None,
        rho,
    )

    assert power == pytest.approx(0.7864653163555486)


def test_rether_pb_equilibrium_volume_average():
    """A solved equilibrium averages a uniform profile to the local value."""

    class _Equilibrium:
        def __init__(self):
            self.rho = np.linspace(0.0, 1.0, 30)
            weight = np.linspace(1.0, 2.0, 30)
            self.R = weight[:, None] * np.ones((1, 4))
            self.J = np.ones((30, 4))
            self.grid = self

        def integrate(self, fun):
            radial = np.asarray(fun, dtype=float).mean(axis=1)
            return float(np.trapezoid(radial * 2.0 * self.rho, self.rho))

    rho = np.linspace(0.0, 1.0, 21)
    power = rether_pb(
        np.full_like(rho, 1.0e20),
        np.full_like(rho, 10.0),
        np.full_like(rho, 11.0),
        0.7,
        1,
        _Equilibrium(),
        rho,
    )

    assert power == pytest.approx(_UNIFORM_POWER_MW)
