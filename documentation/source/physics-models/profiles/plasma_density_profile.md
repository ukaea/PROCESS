# Density Profile | `NeProfile(Profile)`

The density profile class is organised around a central runner function that is called each time the plasma is parameterised by the parent [`PlasmaProfile()`](./plasma_profiles.md) class. It is called by [`pedestal_parameterisation()`](plasma_profiles.md#pedestal_parameterisation) and [`parabolic parameterisation()`](./plasma_profiles.md#parabolic_paramterisation). The sequence of the runner function can be seen below along with explanation of the following calculations.

## Runner function | `run()`

1. Firstly the profile x-dimension is normalised in [`normalise_profile_x()`](./plasma_profiles_abstract_class.md/#normalise-the-profile-in-x--normalise_profile_x) by simply dividing the profile size by its max value.

2. The steps between the normalised points is then done by [`calculate_profile_dx()`](./plasma_profiles_abstract_class.md#calculate-the-profile-steps-in-x--calculate_profile_dx) which divides the max x-dimension by the number of points.

3. The core electron and ion temperatures are claculated via [`set_physics_variables()`]()

### Calculate core values | `set_physics_variables()`

The core electron density is calculated using the [`ncore`](plasma_density_profile.md#electron-core-density-of-a-pedestalised-profile--ncore) method.
The core ion density is then set from $n_{\text{i}}$ (`nd_plasma_ions_total_vol_avg`) which is the total ion density such as:

$$
n_{\text{i0}} = \left(\frac{n_\text{i}}{n_\text{e}}\right)n_{\text{e0}}
$$

#### Electron core density of a pedestalised profile | `ncore()`

This function calculates the core electron density for a pedestalsied profile (`i_plasma_pedestal == 1`). It takes in values of

| Profile parameter / Input               | Density   |
|----------------------------------|-----------|
| Pedestal radius (r/a)            | `radius_plasma_pedestal_density_norm`, $\rho_{\text{ped,n}}$ |
| Pedestal value                   | `nd_plasma_pedestal_electron`, $n_{\text{ped}}$ |
| Separatrix value                 | `nd_plasma_separatrix_electron`, $n_{\text{sep}}$ |
| Average density             | `nd_plasma_electrons_vol_avg`, $\langle n \rangle$ |
| Profile index/ peaking parameter | `alphan`, $\alpha_n$ |

$$
n_0  =  \frac{1}{3\rho_{\text{ped,n}}^2}\left[3\langle n_{\text{e}} \rangle (1+\alpha_n)
+ n_{\text{sep}} (1+\alpha_n) \left(-2 + \rho_{\text{ped,n}} + \rho_{\text{ped,n}}^2\right) \\
- n_{\text{ped}}\left( (1 + \alpha_n)(1+ \rho_{\text{ped,n}}) + (\alpha_n -2)
\rho_{\text{ped,n}}^2\right)\right]
$$

If `ncore` is returned as being less than 0, it is forced into a state of `ncore = 1E-6` in order to help convergence. This will also give a warning to the user to raise the lower bound of the average electron density `nd_plasma_electrons_vol_avg`.

##### Derivation

We calculate the volume integrated profile and then divide by the volume of integration to get the volume average density $\langle n_{\text{e}} \rangle$. If we assume the plasma to be a torus of circular cross-section then we can use spherical coordinates. We can simplify the problem by representing the torus as a cylinder of height equal to the circumference of the torus equal to $2\pi R$ where $R$ is the major radius of the torus, and $a$ is the plasma minor radius in the poloidal plane.

The cylindrical volume element is given by:

$$
V = \int \int \int dV = \int^{2\pi R}_0 \int^{2\pi}_0 \int^a_0 r \ dr \ d\theta \ dz
$$

Inserting our density function to be integrated over we get in normalised radial coordinates ($\rho$) we get:

$$
\int^{2\pi R}_0 \int^{2\pi}_0 \int^{1}_0       \rho \ n_{\text{e}}(\rho) \ d\rho \ d\theta \ dz
$$

Since our density function is only a function of $\rho$, and the torus is symmetric around its center, the integration simplifies to integrating over $\rho$ and the $d\theta ,\ dz$ integrals are solved to give values for the full poloidal angle and cylindrical height / torus length, leading to:

$$
4\pi^2R \int^{1}_0     \rho \ n_{\text{e}}(\rho) \ d\rho  
$$

This is the general form for the full profile width without expansion. Separating out the density function into its separate functions for the core and pedestal region we get the fully expanded integration form.

$$
4\pi^2R\left[ \int^{\rho_{\text{ped,n}}}_0     \rho\left(n_{\text{ped}} + (n_0 - n_{\text{ped}}) \left( 1 -
\frac{\rho^2}{\rho_{\text{ped},n}^2}\right)^{\alpha_n}\right) \ d\rho \\
+\int^1_{\rho_{\text{ped,n}}}     \rho\left(n_{\text{sep}} + (n_{\text{ped}} - n_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},n}}\right)\right)\right] \ d\rho
$$

Integrating each part within its bounds:

$$
4\pi^2R\left[ \frac{\left(n_{\text{ped}} {\alpha}_{n} + n_{0}\right) {\rho}_{\text{ped,n}}^{2}}{2{\alpha}_{n} + 2} \\
+\frac{\left(1-{\rho}_{\text{ped,n}}\right) \left(\left(n_{\text{sep}} + 2n_{\text{ped}}\right) {\rho}_{\text{ped,n}} + 2n_{\text{sep}} + n_{\text{ped}}\right)}{6}\right]
$$

In the form of volume average density where the volume integrated density function has to be divided by the volume of the cylinder / torus, within the volume bounded by that pedestal position we get:

$$
\langle n_{\text{e}} \rangle = 4\pi^2R\left[ \frac{\frac{\left(n_{\text{ped}} {\alpha}_{n} + n_{0}\right) {\rho}_{\text{ped,n}}^{2}}{2{\alpha}_{n} + 2}
+\frac{\left(1-{\rho}_{\text{ped,n}}\right) \left(\left(n_{\text{sep}} + 2n_{\text{ped}}\right) {\rho}_{\text{ped,n}} + 2n_{\text{sep}} + n_{\text{ped}}\right)}{6}}{2\pi^2 R \rho^2}\right]
$$

In this case, the value of $\rho$ is equal to 1 as we integrated over the full profile.

$$
\langle n_{\text{e}} \rangle = 2\left[\frac{\left(n_{\text{ped}} {\alpha}_{n} + n_{0}\right) {\rho}_{\text{ped,n}}^{2}}{2{\alpha}_{n} + 2} \\
+\frac{\left(1-{\rho}_{\text{ped,n}}\right) \left(\left(n_{\text{sep}} + 2n_{\text{ped}}\right) {\rho}_{\text{ped,n}} + 2n_{\text{sep}} + n_{\text{ped}}\right)}{6}\right]
$$

$$
\langle n_{\text{e}} \rangle = \frac{(n_0+n_{\text{ped}}\alpha_\text{n})\rho_{\text{ped,n}}^2}{1+\alpha_n}+\frac{1}{3}(1-\rho_{\text{ped}})(n_{\text{ped}}+2n_{\text{ped}}\rho_{\text{ped}}+n_{\text{sep}}(2+\rho_{\text{ped,n}}))
$$

The above is then rearranged to get a function for $n_0$

$$
n_0  =  \frac{1}{3\rho_{\text{ped,n}}^2}\left[3\langle n_{\text{e}} \rangle (1+\alpha_n)
+ n_{\text{sep}} (1+\alpha_n) \left(-2 + \rho_{\text{ped,n}} + \rho_{\text{ped,n}}^2\right) \\
- n_{\text{ped}}\left( (1 + \alpha_n)(1+ \rho_{\text{ped,n}}) + (\alpha_n -2)
\rho_{\text{ped,n}}^2\right)\right]
$$

$\blacksquare$

------

4. The y profile is then calculated using [`calculate_profile_y()`](plasma_density_profile.md#calculate-density-at-each-radius-position-calculate_profile_y). This routine calculates the density at each normalised minor radius position $\rho$ for a HELIOS-type density pedestal profile[^1]

### Calculate density at each radius position | `calculate_profile_y()`

A table of the input variables can be found below

| Profile parameter / Input               | Density   |
|----------------------------------|-----------|
| Normalized plasma radii            | `profile_x` |
| Pedestal radius (r/a)            | `radius_plasma_pedestal_density_norm`, $\rho_{\text{ped,n}}$ |
| Core density                | `nd_plasma_electron_on_axis`, $n_{\text{e0}}$ |
| Pedestal value                   | `nd_plasma_pedestal_electron`, $n_{\text{ped}}$ |
| Separatrix value                 | `nd_plasma_separatrix_electron`, $n_{\text{sep}}$ |
| Profile index/ peaking parameter | `alphan`, $\alpha_n$ |

If `i_plasma_pedestal == 0` then the original parabolic profile form is used

$$
n(\rho) = n_0(1 - \rho^2)^{\alpha_n}
$$

The central density ($n_0$) is then checked to make sure it is not less than the pedestal density, $n_{\text{ped}}$.
If it is less than a logger warning is pushed to the terminal at runtime.

Values of the profile density are then assigned based on the density function below across bounds from 0 to `radius_plasma_pedestal_density_norm` and `radius_plasma_pedestal_density_norm` to 1.  

$$\begin{aligned}
\mbox{Density:} \ n(\rho) = \left\{
\begin{aligned}
    & n_{\text{ped}} + (n_0 - n_{\text{ped}}) \left( 1 -
    \frac{\rho^2}{\rho_{\text{ped,n}}^2}\right)^{\alpha_n}
& \ 0 \leq \rho \leq \rho_{\text{ped,n}} \\
& n_{\text{sep}} + (n_{\text{ped}} - n_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped,n}}}\right)
& \ \rho_{\text{ped,n}} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$

5. Profile is then integrated with `integrate_profile_y()` using Simpsons integration from the profile abstract base class

[^1]: Jean, J. (2011). *HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies*. Fusion Science and Technology, 59(2), 308–349. <https://doi.org/10.13182/FST11-A11650>
