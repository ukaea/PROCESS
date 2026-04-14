import numpy as np
import pytest
from scipy.optimize import approx_fprime

from process.models.tfcoil.quench import (
    _copper_irradiation_resistivity,
    _copper_magneto_resistivity,
    _copper_rrr_resistivity,
    _copper_specific_heat_capacity,
    calculate_quench_protection_current_density,
)


def test_copper_irradiation():
    """
    Data as read off from:

    M. Nakagawa et al., "High-dose neutron-irradiation
    effects in fcc metals at 4.6 K", *Phys. Rev. B*, 16, 5285 (1977).
    https://doi.org/10.1103/PhysRevB.16.5285  Figure 6

    and listed in:

    M. Kovari, 09/11/2012, internal notes (Excel / Mathcad), Technology Program,
    WP12, PEX, Super-X Divertor for DEMO
    """
    fluence = 1e20 * np.array([
        0.24,
        0.45,
        0.60,
        0.79,
        1.06,
        1.26,
        1.56,
        1.89,
        2.23,
        2.68,
        2.92,
        3.13,
        3.30,
        3.59,
        3.90,
        4.20,
        4.71,
        5.04,
        5.34,
        5.56,
        5.82,
        6.03,
        6.23,
        6.58,
        6.86,
        7.38,
        7.75,
        8.11,
        8.40,
        8.51,
    ])

    nu_irr_data = 1e-11 * np.array([
        16.32,
        33.80,
        44.12,
        58.44,
        76.54,
        89.38,
        106.13,
        124.96,
        142.84,
        163.72,
        173.59,
        181.75,
        188.98,
        200.51,
        211.92,
        221.67,
        238.77,
        249.15,
        257.36,
        262.83,
        268.51,
        273.56,
        278.82,
        286.39,
        293.10,
        303.27,
        309.51,
        316.05,
        320.58,
        322.19,
    ])

    nu_irr = _copper_irradiation_resistivity(fluence)
    assert np.allclose(nu_irr, nu_irr_data)


def test_copper_specific_heat_capacity():
    """
    Data from NIST:

    J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of Copper and Copper Alloys
    at Cryogenic Temperatures", U.S. Government Printing Office, February 1992.
    """
    temperature_k = np.array([
        4,
        6,
        8,
        10,
        14,
        20,
        30,
        50,
        80,
        100,
        120,
        160,
        200,
        260,
    ])
    heat_capacity_j_per_kg_k = np.array([
        0.09,
        0.23,
        0.48,
        0.92,
        2.35,
        7.49,
        26.7,
        96.9,
        203,
        252,
        287,
        331,
        356,
        377,
    ])

    cp = _copper_specific_heat_capacity(temperature_k)
    # NIST quotes errors of 5 - 20 % below 100 K... (comparable to data scatter)
    assert np.allclose(cp, heat_capacity_j_per_kg_k, rtol=0.24)


@pytest.mark.parametrize("rrr", [100, 300, 1000])
def test_copper_rrr_resistivity(rrr):
    t = 4.0
    rrr1 = _copper_rrr_resistivity(t, rrr=1)
    rrr100 = _copper_rrr_resistivity(t, rrr=rrr)
    ratio = rrr1 / rrr100
    assert np.isclose(ratio, rrr, rtol=1e-3)


@pytest.mark.parametrize("rrr", [1, 100, 300, 1000])
def test_copper_magneto_resistivity(rrr):
    t = 4.0
    rrr = _copper_rrr_resistivity(t, rrr=rrr)
    nu = _copper_magneto_resistivity(rrr, 0.0)
    assert rrr == nu
    nu = _copper_magneto_resistivity(rrr, 10.0)
    assert rrr < nu
    nu2 = _copper_magneto_resistivity(rrr, 20.0)
    assert nu < nu2


def test_calculate_quench_protection_intuitive_gradient():
    x = [23.0, 11.0, 0.7, 0.2, 4.75, 150, 100, 1.0, 3e21]
    grad = approx_fprime(
        x, lambda x: calculate_quench_protection_current_density(*x), epsilon=1e-6
    )
    # Increasing discharge time will decrease the protection current
    assert grad[0] < 0.0
    # Increasing field will decrease the protection current
    assert grad[1] < 0.0
    # Increasing copper fraction will increase the required protection current
    assert grad[2] > 0.0
    # Increasing helium fraction will decrease the protection current
    assert grad[3] < 0.0
    # Increasing start temperature will decrease the protection current
    # This is because this decreases the range over which the integral is performed
    assert grad[4] < 0.0
    # Increasing max temperature will increase the protection current
    # This is because this increases the range over which the integral is performed
    assert grad[5] > 0.0
    # Increasing RRR will increase protection current
    assert grad[6] > 0.0
    # Increasing detection time will decrease the protection current
    assert grad[7] < 0.0
    assert np.isclose(grad[7], 2.0 * grad[0])
    # Increasing fluence will decrease the required protection current
    assert grad[8] < 0.0
