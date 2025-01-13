# Plasma Radiation


## Impurity Radiation Class | `ImpurityRadiation`

### Initialization | `__init__()`

Initialize the FusionReactionRate class with the given plasma profile.

#### Parameters:
- `plasma_profile (PlasmaProfile)`: The parameterized temperature and density profiles of the plasma. Taken from the plasma_profile object.

#### Attributes:
- `plasma_profile` (`PlasmaProfile`): Plasma profile instance.

- `rho` (`numpy.ndarray`): Density profile along the x-axis.

- `rhodx` (`numpy.ndarray`): Density profile step size along the x-axis.

- `imp` (`numpy.ndarray`): Indices of impurities with a fraction greater than 1.0e-30.

- `pimp_profile` (`numpy.ndarray`): Impurity profile array initialized to zeros.

- `radtot_profile` (`numpy.ndarray`): Total radiation profile array initialized to zeros.

- `radcore_profile` (`numpy.ndarray`): Core radiation profile array initialized to zeros.

- `radtot` (`float`): Total radiation initialized to 0.0.

- `radcore` (`float`): Core radiation initialized to 0.0.


## Synchrotron radiation | `psync_albajar_fidone()`

The synchrotron radiation power[^1] [^2] is assumed to originate from the 
plasma core. 

$$
\begin{aligned}
P_{s y n, r}(\mathrm{MW})= & 3.84 \times 10^{-8}(1-r)^{1 / 2} R a^{1.38} \kappa^{0.79} \\
& \times B_t^{2.62} n_{e_0(20)}^{0.38} T_{e_0}\left(16+T_{e_0}\right)^{2.61} \\
& \times\left(1+0.12 \frac{T_{e_0}}{p_{a_0}^{0.41}}\right)^{-1.51} \\
& \times K\left(\alpha_n, \alpha_T, \beta_T\right) G(A)
\end{aligned}
$$

$$
K(\alpha_n, \alpha_T, \beta_T) = \left(\alpha_n +3.87\alpha_T +1.46\right)^{-0.79} \\
\times \left(1.98+\alpha_T\right)^{1.36}\beta_T^{2.14} \\
\times \left(\beta_T^{1.53}+1.87\alpha_T-0.16\right)^{-1.33}
$$

$$
G\left(A\right) = 0.93\left[1+0.85 e^{-0.82A}\right]
$$

$$
p_{a_0} = 6.04 \times 10^3 \frac{a n_{e_0(20)}}{B_T}
$$



--------------------

By changing the input parameter `coreradius`, the user may set the normalised 
radius defining the 'core' region. Only the impurity and synchrotron radiation 
from this affects the confinement scaling. Figure 1 below shows the
radiation power contributions.

The wall reflection factor `f_sync_reflect` may be set by the user.


![Schematic diagram of radiation power contributions](../images/radiation.png "Schematic diagram of radiation power contributions")
*Figure 1: Schematic diagram of the radiation power contributions and how they are split between core and edge radiation*

------------------

## Key Constraints

### Radiation power density upper limit

This constraint can be activated by stating `icc = 17` in the input file.

$$
\frac{P_{\text{rad}}}{V_{\text{plasma}}} < \frac{P_{\text{inj}} + P_{\alpha}f_{\alpha, \text{plasma}} + P_{\text{charged}} + P_{\text{ohmic}} }{V_{\text{plasma}}}
$$

Ensures that the calculated total radiation power density does not exceed the total
power available that can be converted to radiation (i.e. the sum of the fusion
alpha power coupled to the plasma, other charged particle fusion power, auxiliary injected power and
the ohmic heating power). $f_{\alpha, \text{plasma}}$ is the fraction of alpha power that is coupled to the plasma (`f_alpha_plasma`).

The scaling value `fradpwr` can be varied also.

**It is recommended to have this constraint on as it is a plasma stability model**

----------------

### Radiation wall load upper limit

This constraint can be activated by stating `icc = 67` in the input file.

The limiting value of $q_{\text{fw,rad}}$ in $\mathrm {MWm^{-2}}$ is be set using input parameter `maxradwallload`.

The scaling value `fradwall` can be varied also.


[^1]: F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron radiation losses in realistic tokamak plasmas,” Nuclear Fusion, vol. 41, no. 6, pp. 665–678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.
‌
[^2]: I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755–1758, Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
‌