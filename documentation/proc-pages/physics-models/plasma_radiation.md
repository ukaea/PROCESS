# Plasma Radiation

The radiation per unit volume is determined using loss functions computed by the `ADAS405` code [^1].
The effective collisional–radiative coefficients necessary to determine the ionization state and radiative losses of each ionic species, 
assuming equilibrium ionization balance in an optically thin plasma, were sourced from `ADF11`-derived data files [^2]. 
These coefficients utilize the generalized collisional-radiative approach [^3] for elements such as H, He, Be, C, N, O, Ne, and Si. 

For Ni, the data is based on [^4], for Fe it is derived from [^5]; and for W, the data is obtained from [^6].
The Ni and Fe rates incorporate a density dependence as described in [^7]. For Kr and Xe, data is taken from the ADAS baseline.

The computed loss functions exhibit a weak dependence on density but are evaluated at a fixed electron density of $10^{19} \text{m}^{-3}$.
This differs from strict coronal equilibrium, which assumes density independence. 
In practice, non-local effects arising from density and temperature gradients are significant but are not considered here. 
The loss functions account for Bremsstrahlung, line radiation, and recombination radiation, represented by:

$$
P_i = n_i n_e L_Z (Z_i, T)
$$

where \(P_i\) is the radiation per unit volume (excluding synchrotron radiation),
\(L_Z (Z_i, T)\) is the loss function for ion species \(i\) at temperature \(T\), and \(n_i\) is the density of ion species \(i\).

The radiation emission is numerically integrated over the plasma profile, 
using the corresponding temperature and density distributions. Emission is only considered from within the separatrix, 
as the PROCESS model does not account for the scrape-off layer. 
The plasma inside the separatrix is divided into two regions: the “core” and the “edge,” separated by a normalized minor radius defined by the user. 
Radiation is calculated separately for the core and edge regions, except for synchrotron radiation, which is assumed to originate solely from the core.

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

The synchrotron radiation power[^8][^9] is assumed to originate from the plasma core. 

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

The limiting value of $q_{\text{fw,rad}}$ in $\mathrm {MWm^{-2}}$ is be set using input parameter `pflux_fw_rad_max`.

The scaling value `fradwall` can be varied also.

[^1]: “ADAS: Docmentation,” Adas.ac.uk, 2024. https://www.adas.ac.uk/manual.php
[^2]: “OPEN-ADAS,” Adas.ac.uk, 2025. https://open.adas.ac.uk/adf11 (accessed Jan. 15, 2025).
[^3]: H. P. Summers et al., “Ionization state, excited populations and emission of impurities in dynamic finite density plasmas: I. The generalized collisional–radiative model for light elements,” Plasma Physics and Controlled Fusion, vol. 48, no. 2, pp. 263–293, Jan. 2006, doi: https://doi.org/10.1088/0741-3335/48/2/007.
[^4]: M. Arnaud and R. Rothenflug, “An updated evaluation of recombination and ionization rates,” Astronomy & Astrophysics Supplement Series, vol. 60, no. 3, pp. 425–457, Jun. 1985.
[^5]: M. Arnaud, R. Rothenflug, An updated evaluation of recombination and ionization rates, Astron. & Astrophys. Supp. Ser. 60 (1985) 425–457.
[^6]: T. Pütterich, R. Neu, R. Dux, A. D. Whiteford, and M. G. O’Mullane, “Modelling of measured tungsten spectra from ASDEX Upgrade and predictions for ITER,” Plasma Physics and Controlled Fusion, vol. 50, no. 8, p. 085016, Jun. 2008, doi: https://doi.org/10.1088/0741-3335/50/8/085016.
[^7]: Summers, H. P. (1974). Tables and Graphs of Collisional Dielectronic Recombination and Ionisation Coefficients and Ionisation Equilibria of H-like to A-like Ions of Elements. Appleton Laboratory.
[^8]: F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron radiation losses in realistic tokamak plasmas,"Nuclear Fusion, vol. 41, no. 6, pp. 665–678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.
[^9]: I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755–1758, Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
