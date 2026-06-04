"""Engineering models for pumping system analysis."""

import logging
from enum import IntEnum
from types import DynamicClassAttribute

import numpy as np


class CoolantType(IntEnum):
    """Enum for coolant types."""

    HELIUM = (1, "Helium")
    WATER = (2, "Water")

    def __new__(cls, value: int, full_name: str):
        """Create a new CoolantType enum member with value and full_name.

        Parameters
        ----------
        value : int
            The numeric value of the enum member.
        full_name : str
            The full name description of the coolant type.
            This should match the CoolProp name for the coolant type.

        Returns
        -------
        CoolantType
            A new enum member with the specified value and full_name.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        """The full name of the coolant type."""
        return self._full_name_


from process.core.exceptions import ProcessValueError

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
    [1] https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
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
    [1] https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation

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

    # Check that the Darcy friction factor is in valid range for the Gnielinski
    # correlation
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


def elbow_coeff(
    radius_pipe_elbow: float,
    deg_pipe_elbow: float,
    darcy_friction: float,
    dia_pipe: float,
) -> float:
    """Calculates elbow bend coefficients for pressure drop calculations.

    Parameters
    ----------
    radius_pipe_elbow : float
        Pipe elbow radius (m)
    deg_pipe_elbow : float
        Pipe elbow angle (degrees)
    darcy_friction : float
        Darcy friction factor
    dia_pipe : float
        Pipe diameter (m)

    Returns
    -------
    float
        Elbow coefficient for pressure drop calculation

    References
    ----------
    [1] Idel'Cik, I. E. (1969), Memento des pertes de charge,
    Collection de la Direction des Etudes et Recherches d'Electricité de France.
    """
    if deg_pipe_elbow == 90:
        a = 1.0
    elif deg_pipe_elbow < 70:
        a = 0.9 * np.sin(deg_pipe_elbow * np.pi / 180.0)
    elif deg_pipe_elbow > 100:
        a = 0.7 + (0.35 * np.sin((deg_pipe_elbow / 90.0) * (np.pi / 180.0)))
    else:
        raise ProcessValueError(
            "No formula for 70 <= elbow angle(deg) <= 100, only 90 deg option available in this range."
        )

    r_ratio = radius_pipe_elbow / dia_pipe

    if r_ratio > 1:
        b = 0.21 / r_ratio**0.5
    elif r_ratio < 1:
        b = 0.21 / r_ratio**2.5
    else:
        b = 0.21

    # Singularity
    ximt = a * b

    # Friction
    xift = (
        (np.pi / 180.0)
        * darcy_friction
        * (radius_pipe_elbow / dia_pipe)
        * deg_pipe_elbow
    )

    # Elbow Coefficient
    return ximt + xift
