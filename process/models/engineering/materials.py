"""Materials module for engineering calculations and properties."""

import logging

from process.data_structure import (
    fwbs_variables,
)

logger = logging.getLogger(__name__)


def eurofer97_thermal_conductivity(temp: float) -> float:
    """Calculates the thermal conductivity of the first wall material (Eurofer97).

    Parameters
    ----------
    temp:
        Property temperature in Kelvin (K).

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
        * fwbs_variables.fw_th_conductivity
        / 28.34
    )
