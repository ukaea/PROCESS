# Temperature Profile | `TeProfile(Profile)`

The temperature profile class is organised around a central runner function that is called each time the plasma is parameterised by the parent [`PlasmaProfile()`](plasma_profiles.md#plasma-profile-class-plasmaprofile) class. It is called by [`pedestal_parameterisation()`](plasma_profiles.md#pedestal_parameterisation) and [`parabolic parameterisation()`](plasma_profiles.md#parabolic_paramterisation). The sequence of the runner function can be seen below along with explanation of the following calculations.

## Runner function | `run()`

1. Firstly the profile x-dimension is normalised in [`normalise_profile_x()`](./plasma_profiles_abstract_class.md/#normalise-the-profile-in-x--normalise_profile_x) by simply dividing the profile size by its max value.

2. The steps between the normalized points is then done by [`calculate_profile_dx()`](./plasma_profiles_abstract_class.md#calculate-the-profile-steps-in-x--calculate_profile_dx) which divides the max x-dimension by the number of points.

3. The core electron and ion temperatures are claculated via [`set_physics_variables()`]()

### Calculate core values | `set_physics_variables()`

The core electron temperature is calculated using the [`tcore`](plasma_temperature_profile.md#electron-core-density-of-a-pedestalised-profile--tcore) method.

#### Electron core density of a pedestalised profile | `tcore()`

This function calculates the core electron density for a pedestalsied profile in $\text{keV}$. The inclusion of a new $\beta_T$ exponent term allows a more accurate description of temperature profiles with a triangular shape or a strong gradient near the pedestal (characteristic of regimes with an [internal transport barrier](https://wiki.fusion.ciemat.es/wiki/Internal_Transport_Barrier)).

A list of input parameters for calculating the core plasma temperature can be found below.

| Profile parameter / Input               | Temperature   |
|----------------------------------|-----------|
| Pedestal radius (r/a)            | `radius_plasma_pedestal_density_norm`, $\rho_{\text{ped,T}}$ |
| Pedestal value                   | `nd_plasma_pedestal_electron`, $T_{\text{ped}}$ |
| Separatrix value                 | `nd_plasma_separatrix_electron`, $T_{\text{sep}}$ |
| Average temperature             | `nd_plasma_electrons_vol_avg`, $\langle T_\text{e} \rangle$ |
| Profile index/ peaking parameter | `alphan`, $\alpha_T$ |
| Profile index/ peaking parameter | `tbeta`, $\beta_T$ |

$$
T_0 = T_{\text{ped}}+\frac{\beta_T(3\langle T_{\text{e}} \rangle +T_{\text{sep}}(-2+\rho_{\text{ped}} +\rho_{\text{ped}}^2)-T_{\text{ped}}(1+\rho_{\text{ped}}+\rho_{\text{ped}}^2))}{6\rho_{\text{ped}}^2\text{B}\left[1+\alpha_T,\frac{2}{\beta_T}\right]}
$$

Where $\text{B}$ is the [Beta function](https://en.wikipedia.org/wiki/Beta_function)

---------

##### Derivation

We calculate the volume integrated profile and then divide by the volume of integration to get the volume average density $\langle T_{\text{e}} \rangle$. If we assume the plasma to be a torus of circular cross-section then we can use spherical coordinates. We can simplify the problem by representing the torus as a cylinder of height equal to the circumference of the torus equal to $2\pi R$ where $R$ is the major radius of the torus, and $a$ is the plasma minor radius in the poloidal plane.

The cylindrical volume element is given by:

$$
V = \int \int \int dV = \int^{2\pi R}_0 \int^{2\pi}_0 \int^a_0 r \ dr \ d\theta \ dz
$$

Inserting our temperature function to be integrated over we get in normalized radial cordinates ($\rho$) we get:

$$
\int^{2\pi R}_0 \int^{2\pi}_0 \int^{1}_0       \rho \ T_{\text{e}}(\rho) \ d\rho \ d\theta \ dz
$$

Since our temperature function is only a function of $\rho$, and the torus is symmetric around its center, the integration simplifies to integrating over $\rho$ and the $d\theta ,\ dz$ integrals are solved to give values for the full poloidal angle and cylindrical height / torus length, leading to:

$$
4\pi^2R \int^{1}_0     \rho \ T_{\text{e}}(\rho) \ d\rho  
$$

This is the general form for the full profile width without expansion. Separating out the temperature function into its separate functions for the core and pedestal region we get the fully expanded integration form.

$$
4\pi^2R\left[ \int^{\rho_{\text{ped,T}}}_0     \rho\left(T_{\text{ped}} + (T_0 - T_{\text{ped}}) \left( 1 -
\frac{\rho^{\beta_T}}{\rho_{\text{ped},T}^{\beta_T}}\right)^{\alpha_T}\right) \ d\rho \\
+\int^1_{\rho_{\text{ped,T}}}     \rho\left(T_{\text{sep}} + (T_{\text{ped}} - T_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},T}}\right)\right)\right] \ d\rho
$$

In the form of volume average temperature where the volume integrated temperature function has to be divided by the volume of the cylinder / torus, within the volume bounded by that pedestal position we get:

$$
\langle T_{\text{e}} \rangle = 4\pi^2R\left[ \frac{\frac{\left(T_{\text{ped}}\beta_T+(2T_0-2T_{\text{ped}})B\left(\alpha_T+1,\frac{2}{\beta_T}\right)\right)\rho_{\text{ped},T}^2}{2\beta_T}+\frac{(1-\rho_{\text{ped},T})\left((T_{\text{sep}}+2T_{\text{ped}})\rho_\text{ped}+2T_{\text{sep}}+T_{\text{ped}}\right)}{6}}{2\pi^2 R \rho^2}\right]
$$

In this case, the value of $\rho$ is equal to 1 as we integrated over the full profile.

$$
\langle T_{\text{e}} \rangle = 2\left[ \frac{\left(T_{\text{ped}}\beta_T+(2T_0-2T_{\text{ped}})\text{B}\left(\alpha_T+1,\frac{2}{\beta_T}\right)\right)\rho_{\text{ped},T}^2}{2\beta_T} \\
+\frac{(1-\rho_{\text{ped},T})\left((T_{\text{sep}}+2T_{\text{ped}})\rho_\text{ped}+2T_{\text{sep}}+T_{\text{ped}}\right)}{6}\right]
$$

$$
\langle T_{\text{e}} \rangle =  \frac{\left(T_{\text{ped}}\beta_T+(2T_0-2T_{\text{ped}})\text{B}\left[\alpha_T+1,\frac{2}{\beta_T}\right]\right)\rho_{\text{ped},T}^2}{\beta_T} \\
+\frac{(1-\rho_{\text{ped},T})\left((T_{\text{sep}}+2T_{\text{ped}})\rho_\text{ped}+2T_{\text{sep}}+T_{\text{ped}}\right)}{3}
$$

Where $\text{B}$ is the [Beta function](https://en.wikipedia.org/wiki/Beta_function)

Re-arranging to get $T_0$ we get:

$$
T_0 = T_{\text{ped}}+\frac{\beta_T(3\langle T_{\text{e}} \rangle +T_{\text{sep}}(-2+\rho_{\text{ped}} +\rho_{\text{ped}}^2)-T_{\text{ped}}(1+\rho_{\text{ped}}+\rho_{\text{ped}}^2))}{6\rho_{\text{ped}}^2\text{B}\left[1+\alpha_T,\frac{2}{\beta_T}\right]}
$$

$\blacksquare$

-----

The core ion temperature is then set such as:

$$
T_{\text{i0}} = \left(\frac{T_{\text{i}}}{T_{\text{e}}}\right)T_{\text{e0}}
$$

4. The y profile is then calculated using [`calculate_profile_y()`](plasma_temperature_profile.md#calculate-temperature-at-each-radius-position-calculate_profile_y). This routine calculates the temperature at each normalised minor radius position $\rho$ for a HELIOS-type temperature pedestal profile[^1]

-------

### Calculate temperature at each radius position | `calculate_profile_y()`

A table of the input variables can be found below

| Profile parameter / Input               | Density   |
|----------------------------------|-----------|
| Normalised plasma radii            | `profile_x` |
| Pedestal radius (r/a)            | `radius_plasma_pedestal_temp_norm`, $\rho_{\text{ped,T}}$ |
| Core density                | `temp_plasma_electron_on_axis_kev`, $T_{\text{e0}}$ |
| Pedestal value                   | `temp_plasma_pedestal_kev`, $T_{\text{ped}}$ |
| Separatrix value                 | `temp_plasma_separatrix_kev`, $T_{\text{sep}}$ |
| Profile index/ peaking parameter | `alphat`, $\alpha_T$ |
| 2nd profile index/ peaking parameter | `tbeta`, $\beta_T$ |

If `i_plasma_pedestal == 0` then the original parabolic profile form is used

$$
T(\rho) = T_0(1 - \rho^2)^{\alpha_T}
$$

The central temperature ($T_0$) is then checked to make sure it is not less than the pedestal temperature, $T_{\text{ped}}$.
If it is less than a logger warning is pushed to the terminal at runtime.

Values of the profile temperature are then assigned based on the density function below across bounds from 0 to `radius_plasma_pedestal_density_norm` and `radius_plasma_pedestal_density_norm` to 1.  

$$\begin{aligned}
\mbox{Temperature:} \ \ T(\rho) = \left\{
\begin{aligned}
& \text{T}_{\text{ped}} + (T_0 - \text{T}_{\text{ped}}) \left( 1 - \frac{\rho^{\beta_T}}
{\rho_{\text{ped},T}^{\beta_T}}\right)^{\alpha_T}  \  0 \leq \rho \leq \rho_{\text{ped},T} \\
& \text{T}_{\text{sep}} + (\text{T}_{\text{ped}} - \text{T}_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},T}}\right)
\ \rho_{\text{ped},T} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$

5. Profile is then integrated with `integrate_profile_y()` using Simpsons integration from the profile abstract base class

[^1]: Jean, J. (2011). *HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies*. Fusion Science and Technology, 59(2), 308–349. <https://doi.org/10.13182/FST11-A11650>
