# Bosch-Hale Methods

These methods are still kept in `physics_functions.py` but outside the [`FusionReactionRate`](plasma_reactions.md) class.

## Bosch-Hale Constants | `BoschHaleConstants`

The `BoschHaleConstants` class is a data structure designed to hold the constants required for the Bosch-Hale calculation for a given fusion reaction.

### Attributes

- **bg (float)**: Represents the Gamow energy parameter.
- **mrc2 (float)**: Represents the reduced mass energy term.
- **cc1 (float)**: Coefficient for the first term in the Bosch-Hale polynomial.
- **cc2 (float)**: Coefficient for the second term in the Bosch-Hale polynomial.
- **cc3 (float)**: Coefficient for the third term in the Bosch-Hale polynomial.
- **cc4 (float)**: Coefficient for the fourth term in the Bosch-Hale polynomial.
- **cc5 (float)**: Coefficient for the fifth term in the Bosch-Hale polynomial.
- **cc6 (float)**: Coefficient for the sixth term in the Bosch-Hale polynomial.
- **cc7 (float)**: Coefficient for the seventh term in the Bosch-Hale polynomial.



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

This function calculates the integrand for the fusion power integration by evaluating the number of fusion reactions per unit volume per particle volume density (m^3/s). It scales the ion temperature profile by the ratio of the volume-averaged ion to electron temperature and normalizes the density profile by the volume-averaged density. The resulting integrand is used to compute the volume-averaged fusion reaction rate, which can be scaled with the volume-averaged ion density.

1. **Scale Ion Temperature Profile**: 
    - Scale the ion temperature profile by the ratio of the volume-averaged ion to electron temperature.
    ```python
    ion_temperature_profile = (
         physics_variables.ti / physics_variables.te
    ) * plasma_profile.teprofile.profile_y
    ```

2. **Calculate Fusion Reactivity**:
    - Calculate the number of fusion reactions per unit volume per particle volume density using the `bosch_hale_reactivity` function.
    ```python
    sigv = bosch_hale_reactivity(ion_temperature_profile, reaction_constants)
    ```

3. **Normalize Density Profile**:
    - Normalize the density profile by the volume-averaged density.
    ```python
    density_profile_normalised = (
         1.0 / physics_variables.dene
    ) * plasma_profile.neprofile.profile_y
    ```

4. **Compute Fusion Integral**:
    - Calculate the volume-averaged fusion reaction integral.
    ```python
    fusion_integral = (
         2.0
         * plasma_profile.teprofile.profile_x
         * sigv
         * density_profile_normalised**2
    )
    ```

5. **Return Result**:
    - Return the computed fusion integral.
    ```python
    return fusion_integral
    ```

[^1]: H.-S. Bosch and G. M. Hale, “Improved formulas for fusion       cross-sections and thermal reactivities,”Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,doi: https://doi.org/10.1088/0029-5515/32/4/i07.