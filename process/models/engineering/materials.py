"""Materials module for engineering calculations and properties."""

import logging

import numba
import numpy as np

logger = logging.getLogger(__name__)

poisson_steel: float = 0.3
"""Steel Poisson's ratio, Source : https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html"""


def eurofer97_thermal_conductivity(temp: float, fw_th_conductivity: float) -> float:
    """Calculates the thermal conductivity of the first wall material (Eurofer97).

    Parameters
    ----------
    temp:
        Property temperature in Kelvin (K).
    fw_th_conductivity:
        thermal conductivity of first wall material at 293 K (W/m/K)

    Returns
    -------
    :
        Thermal conductivity of Eurofer97 in W/m/K.

    Notes
    -----
    Valid up to about 800 K

    References
    ----------
    - A. A. Tavassoli et al., “Materials design data for reduced activation martensitic
      steel type EUROFER,”
      Journal of Nuclear Materials, vol. 329-333, pp. 257-262, Aug. 2004,
      doi: https://doi.org/10.1016/j.jnucmat.2004.04.020.

    - Tavassoli, F. "Fusion Demo Interim Structural Design Criteria (DISDC)/Appendix A
      Material Design Limit Data/A3. S18E Eurofer Steel."
      CEA, EFDA_TASK_TW4-TTMS-005-D01 (2004)
    """
    # temp in Kelvin
    return (
        (5.4308 + 0.13565 * temp - 0.00023862 * temp**2 + 1.3393e-7 * temp**3)
        * fw_th_conductivity
        / 28.34
    )


@numba.njit(cache=True)
def calculate_tresca_stress(
    stress_x: float | np.ndarray,
    stress_y: float | np.ndarray,
    stress_z: float | np.ndarray,
) -> float | np.ndarray:
    """Calculates the Tresca (maximum shear stress) criterion from three principal
    stress components.

    Parameters
    ----------
    stress_x:
        First principal stress in Pa.
    stress_y:
        Second principal stress in Pa.
    stress_z:
        Third principal stress in Pa.

    Returns
    -------
    :
        Tresca stress (maximum shear stress criterion) in Pa, defined as the
        maximum of |stress_x - stress_y|, |stress_y - stress_z|, |stress_x - stress_z|.
    """
    return np.maximum(
        np.maximum(
            np.abs(stress_x - stress_y),
            np.abs(stress_y - stress_z),
        ),
        np.abs(stress_x - stress_z),
    )


@numba.njit(cache=True)
def calculate_von_mises_stress(
    stress_x: float | np.ndarray,
    stress_y: float | np.ndarray,
    stress_z: float | np.ndarray,
    stress_shear_xy: float | np.ndarray,
    stress_shear_yz: float | np.ndarray,
    stress_shear_zx: float | np.ndarray,
) -> float | np.ndarray:
    """Calculates the von Mises stress criterion from three principal stress components.

    Parameters
    ----------
    stress_x:
        First principal stress in Pa.
    stress_y:
        Second principal stress in Pa.
    stress_z:
        Third principal stress in Pa.
    stress_shear_xy:
        Shear stress in the xy-plane in Pa.
    stress_shear_yz:
        Shear stress in the yz-plane in Pa.
    stress_shear_zx:
        Shear stress in the zx-plane in Pa.

    Returns
    -------
    :
        Von Mises stress in Pa, defined as:
        √(0.5 * ((sx - sy)² + (sy - sz)² + (sx - sz)²) + 6 * (τ_xy² + τ_yz² + τ_zx²))

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    """
    return np.sqrt(
        0.5
        * (
            (stress_x - stress_y) ** 2
            + (stress_y - stress_z) ** 2
            + (stress_z - stress_x) ** 2
            + 6 * (stress_shear_xy**2 + stress_shear_yz**2 + stress_shear_zx**2)
        )
    )
