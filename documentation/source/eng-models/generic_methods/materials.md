# Materials Methods

## Eurofer97 thermal conductivity | `eurofer97_thermal_conductivity()`

The thermal conductivity of the first wall is assumed to be that of Eurofer97 using the relation below[^1] [^2]:
 
$$
K_{\text{Eurofer97}} = \frac{\texttt{fwbs\_variables.fw\_th\_conductivity}}{28.34}\left(5.4308 + 0.13565T - 0.00023862T^2 + 1.3393 \times 10^{-7} T^3\right)
$$

!!! warning Thermal conductivity validity bounds

    The sources for the stated thermal conductivity relation above state that the relation is only valid up to 800K [^1] [^2].

[^1]: A. A. Tavassoli et al., “Materials design data for reduced activation martensitic steel type EUROFER,” Journal of Nuclear Materials, vol. 329–333, pp. 257–262, Aug. 2004, doi: https://doi.org/10.1016/j.jnucmat.2004.04.020.

[^2]: Tavassoli, F. "Fusion Demo Interim Structural Design Criteria (DISDC)/Appendix A Material Design Limit Data/A3. S18E Eurofer Steel." CEA, EFDA_TASK_TW4-TTMS-005-D01 (2004).