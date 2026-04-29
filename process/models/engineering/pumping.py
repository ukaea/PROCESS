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
