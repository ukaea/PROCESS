import logging

import numpy as np

logger = logging.getLogger(__name__)


def darcy_friction_haaland(
    reynolds: float, roughness_channel: float, radius_channel: float
) -> float:
    """Calculate Darcy friction factor using the Haaland equation.

    Parameters
    ----------
    reynolds:
        Reynolds number.
    roughness_channel:
        Roughness of the first wall coolant channel (m).
    radius_channel:
        Radius of the first wall coolant channel (m).

    Returns
    -------
    :
        Darcy friction factor.

    Notes
    -----
        The Haaland equation is an approximation to the implicit Colebrook-White equation.
        It is used to calculate the Darcy friction factor for turbulent flow in pipes.

    References
    ----------
        - https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
    """
    # Bracketed term in Haaland equation
    bracket = (roughness_channel / radius_channel / 3.7) ** 1.11 + 6.9 / reynolds

    # Calculate Darcy friction factor
    return (1.8 * np.log10(bracket)) ** (-2)


def gnielinski_heat_transfer_coefficient(
    mflux_coolant: float,
    den_coolant: float,
    radius_channel: float,
    heatcap_coolant: float,
    visc_coolant: float,
    thermcond_coolant: float,
    roughness_channel: float,
) -> float:
    """Calculate heat transfer coefficient using Gnielinski correlation.

    Parameters
    ----------
    mflux_coolant:
        Coolant mass flux in a single channel (kg/m²/s).
    den_coolant:
        Coolant density (average of inlet and outlet) (kg/m³).
    radius_channel:
        Coolant pipe radius (m).
    heatcap_coolant:
        Coolant specific heat capacity (average of inlet and outlet) (J/kg/K).
    visc_coolant:
        Coolant viscosity (average of inlet and outlet) (Pa.s).
    thermcond_coolant:
        Thermal conductivity of coolant (average of inlet and outlet) (W/m.K).
    roughness_channel:
        Roughness of the coolant channel (m).

    Returns
    -------
    :
        Heat transfer coefficient (W/m²K).

    Notes
    -----
    Gnielinski correlation. Ignore the distinction between wall and
    bulk temperatures. Valid for: 3000 < Re < 5e6, 0.5 < Pr < 2000

    References
    ----------
        - https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation

    """
    # Calculate pipe diameter (m)
    diameter = 2 * radius_channel

    # Calculate flow velocity (m/s)
    vel_coolant = mflux_coolant / den_coolant

    # Calculate Reynolds number
    reynolds = calculate_reynolds_number(
        den_coolant=den_coolant,
        vel_coolant=vel_coolant,
        radius_channel=diameter / 2,
        visc_coolant=visc_coolant,
    )

    # Calculate Prandtl number
    pr = heatcap_coolant * visc_coolant / thermcond_coolant

    # Calculate Darcy friction factor, using Haaland equation
    f = darcy_friction_haaland(
        reynolds=reynolds,
        roughness_channel=roughness_channel,
        radius_channel=radius_channel,
    )

    # Calculate the Nusselt number
    nusselt = (
        (f / 8.0)
        * (reynolds - 1000.0)
        * pr
        / (1 + 12.7 * np.sqrt(f / 8.0) * (pr ** (2 / 3) - 1.0))
    )

    # Calculate the heat transfer coefficient (W/m^2K)
    heat_transfer_coefficient = nusselt * thermcond_coolant / (2.0 * radius_channel)

    # Check that Reynolds number is in valid range for the Gnielinski correlation
    if (reynolds <= 3000.0) or (reynolds > 5.0e6):
        logger.error("Reynolds number out of range : [3e3-5000e3]. %s", reynolds)

    # Check that Prandtl number is in valid range for the Gnielinski correlation
    if (pr < 0.5) or (pr > 2000.0):
        logger.error("Prandtl number out of range : [0.5-2000]. %s", pr)

    # Check that the Darcy friction factor is in valid range for the Gnielinski correlation
    if f <= 0.0:
        logger.error("Negative Darcy friction factor (f). %s", f)

    return heat_transfer_coefficient


def calculate_reynolds_number(
    den_coolant: float, vel_coolant: float, radius_channel: float, visc_coolant: float
) -> float:
    """Calculate Reynolds number for flow in a pipe.

    Parameters
    ----------
    den_coolant:
        Coolant density (average of inlet and outlet) (kg/m³).
    vel_coolant:
        Coolant velocity in a single channel (m/s).
    radius_channel:
        Coolant pipe radius (m).
    visc_coolant:
        Coolant viscosity (average of inlet and outlet) (Pa.s).

    Returns
    -------
    :
        Reynolds number.

    """
    # Calculate pipe diameter (m)
    diameter = 2 * radius_channel

    # Calculate Reynolds number
    return den_coolant * vel_coolant * diameter / visc_coolant
