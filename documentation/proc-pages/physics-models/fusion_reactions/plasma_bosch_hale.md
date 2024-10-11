# Bosch-Hale Methods

These methods are still kept in `physics_functions.py` but outside the [`FusionReactionRate`](plasma_reactions.md) class.

## Bosch-Hale Constants | `BoschHaleConstants`

DataClass which holds the constants required for the Bosch Hale calculation
for a given fusion reaction.

## Volumetric Fusion Rate | `bosch_hale_reactivity()`

This function calcualtes the relative velocity fusion reactivity $\langle \sigma v \rangle$ for each point in the plasma profile based on the temperature.

 |  Input Variable             |    |
    |----------------------------------|-----------|
    | Array of temperature values for the plasma profile            | `temperature_profile`  |
    | Bosch-Hale constants for the specific reaction                   | `reaction_constants`  |

$$
\theta = \frac{\text{T}}{\left[1-\frac{\text{T(C2+T(C4+TC6))}}{1+\text{T(C3+T(C5+TC7))}}  \right]}
$$

$$
\xi = \left(\frac{\text{B}_\text{G}^2}{4\theta}\right)^{\frac{1}{3}}
$$

$$
\langle \sigma v \rangle = \text{C1} \times \theta \times \sqrt{\frac{\xi}{m_{\text{r}}\text{c}^2\text{T}^3}} \times e^{-3\xi}
$$

This will output a numpy.array for of the relative velocity fusion reactivity $\langle \sigma v \rangle$ for each point in the temperature profile in units of $[\text{m}^3\text{s}^{-1}]$ After calculation each value is multiplied by $10^{-6}$ as the original Bosch-Hale calculation[^1] give the output in $[\text{cm}^3\text{s}^{-1}]$







## Fusion Rate Integral | `fusion_rate_integral()`

|  Input Variable             |    |
    |----------------------------------|-----------|
    | PlasmaProfile object            | `plasma_profile`  |
    | Bosch-Hale constants for the specific reaction                   | `reaction_constants`  |

Assign the ion temperature profile by scaling the electron temperature profile by the ratio of the volume averaged ion and electron temperatures.

[^1]: H.-S. Bosch and G. M. Hale, “Improved formulas for fusion       cross-sections and thermal reactivities,”Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,doi: https://doi.org/10.1088/0029-5515/32/4/i07.
