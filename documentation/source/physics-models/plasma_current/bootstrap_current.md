## Overview

Bootstrap current in tokamaks originates from the pressure gradients within the plasma and the resulting collisions between particles. As the plasma pressure varies radially, it creates a differential in particle velocities, leading to a net drift of charged particles. This drift generates a toroidal current, known as the bootstrap current, which flows parallel to the magnetic field lines. The phenomenon is a consequence of the neoclassical transport theory, where the collisional processes in a magnetically confined plasma lead to a self-sustaining current. This current is particularly advantageous as it reduces the need for external current drive systems, thereby enhancing the efficiency and stability of the tokamak operation. The bootstrap current is proportional to the pressure gradient and the collisionality of the plasma, making it a critical factor in the design and operation of advanced tokamak reactors aiming for steady-state fusion.

Some more info can be found [here](https://wiki.fusion.ciemat.es/wiki/Bootstrap_current)
## Selection

The fraction of the plasma current provided by the bootstrap effect
can be either input into the code directly, or calculated using one of eleven
methods, as summarized below:

--------------

### ITER IPDG89 Scaling | `bootstrap_fraction_iter89()`
 
Original empirical ITER bootstrap scaling from the 1989 Physics Design Guidelines[^0]
Is selected by setting `i_bootstrap_current = 1`

Empirical fit for the bootstrap current fraction as:

$$
\frac{I_{\text{BS}}}{I} = C_{\text{BS}}\left(\epsilon^{0.5}\beta_{\text{pa}}\right)^{1.3}
$$

where:

$$
C_{\text{BS}} = 1.32 - 0.235\left(\frac{q_{95}}{q_0}\right) + 0.0185\left(\frac{q_{95}}{q_0}\right)^2
$$

$$
\beta_{\text{pa}} = \frac{\langle p \rangle}{\frac{B_{\text{pa}^2}}{2\mu_0}} = \beta_{\text{tot}}\left(\frac{B_0}{B_{\text{pa}}}\right)^2
$$

$$
B_{\text{pa}} = \frac{I}{5\langle a \rangle}
$$

$$
\langle a \rangle = \left(\frac{V}{2\pi^2 R_0}\right)^{0.5}
$$

Here, $\beta_{\text{tot}}$ is the average total plasma (toroidal) beta. $I$ is given in $\text{MA}$ and $B_0$ is the on-axis toroidal field in Tesla.

------------

### Nevins Scaling | `bootstrap_fraction_nevins()`

The general Nevins scaling is normally cited from a ITER specialists meeting in 1989[^1] which is not publicly accessible.
However it can be found in the appendix [here](https://doi.org/10.1016/j.fusengdes.2014.07.009)[^2].

Is selected by setting `i_bootstrap_current = 2`

$$
f_{\text{BS}} = \frac{2.5\beta_{e0}R_pB_{T}q_{95}}{I_{\text{P}}}\int^1_0 B_{\text{int}}\ dy
$$

$$
\beta_e = \frac{1.6022 \times 10^{-16}n_{\text{e}}T_{\text{e}}}{B_{T}^2 / 2\mu_0}
$$

$$
\beta_{e0} = \frac{1.6022 \times 10^{-16}n_{\text{e0}}T_{\text{e0}}}{B_{T}^2 / 2\mu_0}
$$

$$
\Delta = \epsilon y^{0.5}
$$

$$
x = \frac{1.46 \sqrt{\Delta}+2.4\Delta}{\left(1-\Delta\right)^{1.5}}
$$

$$
d = \sqrt{2}Z_{\text{eff}} + Z_{\text{eff}}^2 + x\left(0.75+2.657Z_{\text{eff}}+2Z_{\text{eff}}^2\right) \\
+x^2\left(0.348+1.243Z_{\text{eff}}+Z_{\text{eff}}^2\right)
$$

$$
A_1 = \left(\alpha_{\text{n}}+\alpha_{\text{T}}\right)\left(1-y\right)^{\alpha_{\text{n}}+\alpha_{\text{T}}-1}
$$

$$
A_2 = \alpha_{\text{T}}\left(1-y\right)^{\alpha_{\text{n}}+\alpha_{\text{T}}-1}
$$

$$
A_{l1} = \frac{x}{d}\left(0.754+2.21Z_{\text{eff}}+Z_{\text{eff}}^2+x\left(0.348+1.243Z_{\text{eff}}+Z_{\text{eff}}^2\right)\right)
$$

$$
A_{l2}=-x\frac{0.884+2.0742Z_{\text{eff}}}{d}
$$

$$
\alpha_i = -\frac{1.172}{1+0.462x}
$$

$$
q = q_0+\left(q_{95}-q_0\right)\frac{y+y^2+y^3}{3}
$$

$$
\beta_{\text{tot}} = \beta_{\text{T}}\frac{B_{\text{tot}}^2}{B_{\text{T}}^2} = \frac{\beta_{T}\left(B_{T}^2+\beta_{\text{P}}\right)}{B_{\text{T}}^2}
$$

$$
P_{\text{ratio}} = \frac{\beta_{\text{T}}-\beta_e}{\beta_e}
$$

$$
B_{\text{int}} = \frac{q}{q_{95}}\left[A_{l1}\left(A_1+P_{\text{ratio}}\left(A_1+\alpha_i A_2\right)\right)+A_{l2}A_2\right]
$$



----------------


### Wilson Scaling | `bootstrap_fraction_wilson()`

Wilson gives an empirical formula[^3] [^7] as a function of the pressure, 
temperature and total current profiles as well as the poloidal beta and aspect ratio of a tokamak. This empirical formula was compared with an expression obtained by the ITER group; also compared with an analytical result (valid at large aspect ratio). It is found that the determined empirical result agreed well with the large aspect ratio result, but not so well with the empirical formula of the ITER group[^0]

Is selected by setting `i_bootstrap_current = 3`

Data is fitted to 3000 equilibria evenly distributed across the parameter range below:

| Variable | $R$ | $A$ | $B_{\text{T}}$ | $\delta$ | $\kappa$ | $\alpha_{\text{P}}$ | $\alpha_{\text{T}}$ | $\alpha_{\text{J}}$ | $Z$ |
|----------|-----|-----|----------------|-----------|-----------|-------------------|-------------------|-------------------|-----|
| Value    | 5.31 | 1.1-5 | 6.2 | 0.2 | 2.0 | 1-3 | 0.1-$\alpha_{\text{P}}$ | 0.5-2.0 | 1-3 |


Using the relation:

$$
\frac{I_{\text{b}}}{I_{\text{P}}} = \beta_{\text{P}}\epsilon_0^{\frac{1}{2}} \sum_{i=1}^{12}a_i(\alpha_{\text{J}}, Z)b_i
$$

In the paper definition $\epsilon_0$ is not the standard inverse aspect ration but is defined as:

$$
\epsilon_0 = \frac{R_2-R_1}{R_2+R_1}
$$

Where $R_2$ and $R_1$ are the maximum and minimum radii of the plasma.

The poloidal beta term is defined as:

$$
\beta_{\text{P}} = \frac{2\mu_0\langle p \rangle}{\langle \langle B_{\text{P}} \rangle \rangle^2}
$$

Where $\langle p \rangle$ is the volume averaged pressure and $\langle \langle B_{\text{P}} \rangle \rangle$  is the plasma surface average of the polidal field.

The Wilson method extrapolates on the analytically derived expression for the bootstrap current fraction in the limit of large aspect ratio with circular cross-section. This large aspect expression assumed coincidence of a constant $r$ surface with a constant flux surface. The allowed the treatment of the temperature, pressure and current density as functions of $r$ only with the simple parabolic profile form. 
For the Wilson method this is expanded to the arbitrary aspect ratio case with D-shaped plasmas with a given triangularity and elongation. The profiles for pressure and temperature are now of the form:

$$
P=P_0\psi^{\alpha_{\text{P}}}
$$

Where $\psi$ is the flux function. For the current density we cannot model it in the form above as it varies across a flux surface. Instead it is convenient to consider a flux surface average of the total parallel current density which is taken to vary as:

$$
\frac{\langle j \cdot B \rangle}{\langle B^2 \rangle^{\frac{1}{2}}} = J_0 \psi^{\alpha_{\text{J}}}
$$

Describing this flux surface average by using a standard **parabolic profile** for the safety factor $(q)$.

$$
q(r) = q_0 +\left(q_{95}-q_0\right)\left(\frac{r}{a}\right)^2
$$

To relate the temperature, pressure and current density profiles across a flux surface average of the flux function we re-arrange their standard parabolic forms into one that can be integrated into the $(q)$ profile function. Using the temperature profile as an example:

$$
T(r) = T_0 \left(1-\left(\frac{r}{a}\right)^2\right)^{\alpha_{\text{T}}}
$$

$$
\frac{T(r)}{T_0} = \left(1-\left(\frac{r}{a}\right)^2\right)^{\alpha_{\text{T}}}
$$

$$
\left(\frac{r}{a}\right)^2 = 1-\left(\frac{T(r)}{T_0}\right)^{\frac{1}{\alpha_{\text{T}}}}
$$

Substituting the above into the $q$ profile and taking a median profile value of half the core value we get:

$$
q_0+(q_{95}-q_0)\times (1-0.5^{\frac{1}{\alpha_{\text{T}}}})
$$

To find the new flux surface average where the $q$ profile shares points where the temperature, pressure and current density are half of their core values, we relate the function above of the calculated $q$ value to the core $q_0$ value. Taking the natural log of the $q$ profile point at 50% of the core value for the temperature, pressure and current density we relate these two ratios compared to $q_{95}$. Dividing $\left(\ln 0.5\right)$ by this result provides new values of the profile exponents $\alpha_{\text{T}}, \alpha_{\text{P}} ,\alpha_{\text{J}}$ which are now averaged across a shared flux function. 


$$
\alpha = \frac{\ln{0.5}}{\ln{\frac{\ln \left(\frac{q_0+(q_{95}-q_0)\times (1-0.5^{\frac{1}{\alpha}})}{q95}\right)}{\ln{\frac{q_0}{q_{95}}}}}}
$$

These new calculated profile exponents $(\alpha)$ are what is used in the matrix values below and **not** the standard parabolic profile indexes.

$$
b_1 = 1 \ \ b_2 = \alpha_{\text{P}} \ \ b_3 = \alpha_{\text{T}} \ \ b_4 = \alpha_{\text{P}}\alpha_{\text{T}} \\
b_5 = \epsilon_0^{\frac{1}{2}} \ \ b_6 = \alpha_{\text{P}}\epsilon_0^{\frac{1}{2}}  \ \ b_7 = \alpha_{\text{T}}\epsilon_0^{\frac{1}{2}} \ \ b_8 = \alpha_{\text{P}}\alpha_{\text{T}}\epsilon_0^{\frac{1}{2}} \\
b_9 = \epsilon_0 \ \ b_{10} = \alpha_{\text{P}}\epsilon_0 \ \ b_{11} = \alpha_{\text{T}}\epsilon_0 \ \ b_{12} = \alpha_{\text{P}}\alpha_{\text{T}}\epsilon_0
$$

$$
a_1 = 1.41\left(1.0 - 0.28 \alpha_{\text{J}}^{\frac{1}{2}}\right)\left(1.0 + \frac{0.12}{Z}\right) \\
a_2 = 0.36  \left(1.0 - 0.59 \alpha_{\text{J}}^{\frac{1}{2}}\right) \left(1.0 + \frac{0.8}{Z}\right) \\
a_3 = -0.27  \left(1.0 - 0.47 \alpha_{\text{J}}^{\frac{1}{2}}\right)  \left(1.0 + \frac{3.0}{Z}\right) \\
a_4 = 0.0053  \left(1.0 + \frac{5.0}{Z}\right) \\
a_5 = -0.93  \left(1.0 - 0.34 \alpha_{\text{J}}^{\frac{1}{2}}\right)  \left(1.0 + \frac{0.15}{Z}\right) \\
a_6 = -0.26  \left(1.0 - 0.57 \alpha_{\text{J}}^{\frac{1}{2}}\right)  \left(1.0 - 0.27 Z\right) \\
a_7 = 0.064  \left(1.0 - 0.6 \alpha_{\text{J}} + 0.15 \alpha_{\text{J}}^2\right)  \left(1.0 + \frac{7.6}{Z}\right) \\
a_8 = -0.0011  \left(1.0 + \frac{9.0}{Z}\right) \\
a_9 = -0.33  \left(1.0 - \alpha_{\text{J}} + 0.33 \alpha_{\text{J}}^2\right) \\
a_{10} = -0.26  \left(1.0 - \frac{0.87}{\alpha_{\text{J}}^{\frac{1}{2}}} - 0.16 \alpha_{\text{J}}\right) \\
a_{11} = -0.14  \left(1.0 - \frac{1.14}{\alpha_{\text{J}}^{\frac{1}{2}}} - 0.45 \alpha_{\text{J}}^{\frac{1}{2}}\right) \\
a_{12} = -0.0069 \\
$$

Coefficients are obtained by a least squares fit to the 3000 equilibria numerical solutions. Error distribution shows an average error of 3.6% and a maximum error of 20%. These larger errors appear to come from the cases where the temperature profile is approximately equal to the pressure profile. For these cases the density profile is very flat over much of the plasma radius and (as it is forced to be zero at the plasma edge) this 
means that it falls off sharply at the plasma edge.

!!! quote "Excerpt from Wilson[^3]"

    *"The ITER group gives an expression for the bootstrap current[^0] which differs significantly from that 
    presented here. Their relatively simple expression,
    shows no variation with density and temperature profiles. It also does not reproduce the large aspect ratio 
    scaling with $\beta_{\text{P}}$, and $\epsilon$ in this limit."*   

!!! warning "Flux average profile indexes"

    The implementation of calculating the new profile indexes to be used in the $a$ & $b$ coefficients is not explicitly given in  Wilson[^3]. It is an adaption to make the method more aligned with what `PROCESS` calculates.    


--------------------

### Sauter Scaling | `bootstrap_fraction_sauter()`

Sauter et al.[^4] [^5] provides a formula using the exact Fokker–Planck operator and
without any approximation on the plasma geometry or collisionality. In this way we have been able to accurately determine the neoclassical resistivity and the coefficients for the
bootstrap current which allows one to calculate the bootstrap fraction.

Is selected by setting `i_bootstrap_current = 4`

$$
\left\langle j_{\|} B\right\rangle= \sigma_{\text {neo }}\left\langle E_{\|} B\right\rangle-I(\psi) p(\psi)\left[\mathcal{L}_{31} \frac{\partial \ln n_e}{\partial \psi} \\
+R_{p e}\left(\mathcal{L}_{31}+\mathcal{L}_{32}\right) \frac{\partial \ln T_e}{\partial \psi}+\left(1-R_{p e}\right) \times\left(1+\frac{\mathcal{L}_{34}}{\mathcal{L}_{31}} \alpha\right) \mathcal{L}_{31} \frac{\partial \ln T_i}{\partial \psi}\right]
$$

Note that the above $\left\langle j_{\|} B\right\rangle$ given by Sauter et.al.[^4]   gives the component of the current density in the direction of the field – not the toroidal component of the current.  The error is second order in the pitch angle, but can be important, especially in the outboard region of a low aspect ratio tokamak, where the pitch angle is large.  Moreover this error will always have the effect of overestimating the current. This is accounted for via poloidal correction function implemented as: [`beta_poloidal_sauter()`](#calculate-electron-only-poloidal-beta-correction-beta_poloidal_sauter) and[`beta_poloidal_total_sauter()`](#calculate-ion-and-electron-poloidal-beta-correction-beta_poloidal_total_sauter)

The correction is of the form:

$$
I_{\phi}^{\text{bs}} = 2\pi \int d\psi \frac{q(\psi)}{\langle B^{2}\rangle} \left\langle j_{\|} B\right\rangle
$$

where $q(\psi)$ is the safety factor.

The reconstructed implementation in PROCESS looks as such:

$$
\frac{I_{\text{b}}}{I_{\text{P}}} = \sum_2^{\rho_{\text{max}}} \left(2\pi\left[\rho\right]_{-1} \times \left(\left[\rho\right] -\left[\rho\right]_{-1}\right)  \right) \times  \\
\left(0.5 \times \left[\mathcal{L}_{31} \frac{\partial \ln n_e}{\partial \psi} 
+\left(\mathcal{L}_{31} + \mathcal{L}_{32}\right) \frac{\partial \ln T_e}{\partial \psi} 
+ \left(1 + \frac{\mathcal{L}_{34}}{\mathcal{L}_{31}} \alpha\right) \mathcal{L}_{31} \frac{\partial \ln T_i}{\partial \psi}\right] \times \\
 1 \times 10^6 \times \frac{-B_{0}\left[\rho\right]_{-1}\left[\frac{1}{q}\right]_{-1}}{0.2\pi R_0}\right)
$$

In this case square brackets denote array variables equal in length to $\rho_{\text{max}}$ representing the normalised radius elements across the profile. The $-1$ subscript denotes the previous array element in the summation.

It is not known fully if the $\left(\sigma_{\text {neo }}\left\langle E_{\|} B\right\rangle-I(\psi) p(\psi)\right)$ term is properly implemented into the `PROCESS` version. The $R_{pe}$ value is assumingly taken to be 0.5 as it is stated to approximately to be in Sauter et.al[^4].

!!! warning "Validity of the Sauter Bootstrap Scaling"

    In its current state the several base functions called by the Sauter scaling have no reference and cannot be verified. The ad-hoc adaption of the Sauter scaling for use in `PROCESS` is done knowing that `PROCESS` does not calculate flux surfaces across the plasma.

-----------

#### Calculate the trapped particle fraction | `_trapped_particle_fraction_sauter()`

This function calculates the trapped particle fraction $\left(f_t\right)$ used within other key internal Sauter scaling functions.

$$
f_t = \frac{1.0 - (1.0-\epsilon_{-1})\sqrt{1.0-\epsilon_{-1}}}{\left(1.0+1.46 \sqrt{\epsilon_{-1}}\right)}
$$

$\epsilon$ in this case is the local aspect ratio at that normalised radial point in the profile given by $\epsilon = \rho \left(\frac{a}{R}\right)$. The value of $\rho$ varies from 0 to 1 across the profile.

The $-1$ subscript in this case refers to the value of the variable in the previous array index value.

-------------

#### Calculate electron density coefficient | `_calculate_l31_coefficient()`

This function calculates and returns the $\mathcal{L}_{31}$ coefficient value for $\frac{\partial \ln n_e}{\partial \psi}$

$$
\mathcal{L}_{31}= F_{31}\left(X=f_{\text {teff }}^{31}\right) \equiv\left(1+\frac{1.4}{Z+1}\right) X-\frac{1.9}{Z+1} X^2+\frac{0.3}{Z+1} X^3 +\frac{0.2}{Z+1} X^4 \\
f_{\text {teff }}^{31}\left(\nu_{e *}\right)= \frac{f_t}{1+\left(1-0.1 f_t\right) \sqrt{\nu_{e *}}+0.5\left(1-f_t\right) \nu_{e *} / Z}
$$

The returned value is $\mathcal{L}_{31} \times$ [`_beta_poloidal_total_sauter()`](#calculate-ion-and-electron-poloidal-beta-correction-beta_poloidal_total_sauter)  

---------------

#### Calculate electron temperature coefficient | `_calculate_l31_32_coefficient()`

This function calculates and returns the $\left(\mathcal{L}_{31}+\mathcal{L}_{32}\right)$ coefficient value for $\frac{\partial \ln T_{\text{e}}}{\partial \psi}$.

$$
\begin{align*}
\mathcal{L}_{32} &= F_{32 \_e e}\left(X=f_{\text {teff }}^{32 \_e e}\right)+F_{32 \_e i}\left(Y=f_{\text {teff }}^{32 \_e i}\right) \\ \\
F_{32 \_e e}(X) &= \frac{0.05+0.62 Z}{Z(1+0.44 Z)}\left(X-X^4\right)+\frac{1}{1+0.22 Z}\left[X^2-X^4\right. \\
& \quad \left.-1.2\left(X^3-X^4\right)\right]+\frac{1.2}{1+0.5 Z} X^4 \\
f_{\text {teff }}^{32_{-e e}}\left(\nu_{e *}\right) &= \frac{f_t}{1+0.26\left(1-f_t\right) \sqrt{\nu_{e *}}+0.18\left(1-0.37 f_t\right) \frac{\nu_{e *}}{\sqrt{Z}}}
\end{align*}
$$

$$
F_{32 \_e i}(Y) = -\frac{0.56+1.93 Z}{Z(1+0.44 Z)}\left(Y-Y^4\right)+\frac{4.95}{1+2.48 Z}\left[Y^2-Y^4\right] \\
-0.55\left(Y^3-Y^4\right)-\frac{1.2}{1+0.5 Z} Y^4 \\
f_{\text {teff }}^{32 \_e i}\left(\nu_{e *}\right) = \frac{f_t}{1+\left(1+0.6 f_t\right) \sqrt{\nu_{e *}}+0.85\left(1-0.37 f_t\right) \nu_{e *}(1+Z)}
$$

The above is added to a call of [`_calculate_l31_coefficient()`](#calculate-electron-density-coefficient-calculate_l31_coefficient). This is then multiplied by [`_beta_poloidal_sauter()`](#calculate-electron-only-poloidal-beta-correction-beta_poloidal_sauter). 

This product above is then multiplied by ([`_beta_poloidal_sauter()`](#calculate-electron-only-poloidal-beta-correction-beta_poloidal_sauter) divided by [`_beta_poloidal_total_sauter()`](#calculate-ion-and-electron-poloidal-beta-correction-beta_poloidal_total_sauter))

---------------

#### Calculate ion temperature coefficient | `_calculate_l34_alpha_31_coefficient()`

This function calculates and returns the $\left(1+\frac{\mathcal{L}_{34}}{\mathcal{L}_{31}}\alpha\right)\mathcal{L}_{31}$ coefficient value for $\frac{\partial \ln T_{\text{i}}}{\partial \psi}$.

$$
\mathcal{L}_{34}=F_{31}\left(X=f_{\text {teff }}^{34}\right) 
$$

$$
f_{\text {teff }}^{34}\left(\nu_{e *}\right)=\frac{f_t}{1+\left(1-0.1 f_t\right) \sqrt{\nu_{e *}}+0.5\left(1-0.5 f_t\right) \nu_{e *} / Z} 
$$

$$
\alpha_0=-\frac{1.17\left(1-f_t\right)}{1-0.22 f_t-0.19 f_t^2} 
$$

$$
\alpha\left(\nu_{i *}\right)=\left[\frac{\alpha_0+0.25\left(1-f_t^2\right) \sqrt{\nu_{i *}}}{1+0.5 \sqrt{\nu_{i *}}} \\
+0.315 \nu_{i *}^2 f_t^6\right] \frac{1}{1+0.15 \nu_{i *}^2 f_t^6}
$$

The definition of $\alpha\left(\nu_{i *}\right)$ is that found in the erratum paper which changes the value of $-0.315\nu_{i *}^2 f_t^6$ to positive.[^5]

The return sequence is ([`_beta_poloidal_total_sauter()`](#calculate-ion-and-electron-poloidal-beta-correction-beta_poloidal_total_sauter) - [`_beta_poloidal_sauter()`](#calculate-electron-only-poloidal-beta-correction-beta_poloidal_sauter)) $\times (\mathcal{L}_{34} + \alpha)$ + [`_calculate_l31_coefficient()`](#calculate-electron-density-coefficient-calculate_l31_coefficient) $\times$ (1.0 -  [`_beta_poloidal_sauter()`](#calculate-electron-only-poloidal-beta-correction-beta_poloidal_sauter) divided by [`_beta_poloidal_total_sauter()`](#calculate-ion-and-electron-poloidal-beta-correction-beta_poloidal_total_sauter))



-------------

#### Calculate the Coulomb logarithm | `_coulomb_logarithm_sauter()`

$$
\ln \Lambda = 15.9 -0.5 \times \ln{n_{\text{e}}}+\ln{T_{\text{e}}}
$$

-----------

#### Calculate frequency of electron collisions | `_electron_collisions_sauter()`

Using the Coulomb logarithm ($\ln \Lambda$) calculated from [`_coulomb_logarithm_sauter()`](#calculate-the-coulomb-logarithm--_coulomb_logarithm_sauter) we get:

$$
\nu_{\text{e}} = 670 \times \frac{\ln \Lambda \times n_{\text{e}}}{T_{\text{e}}^{3/2}}
$$

------------

#### Calculate electron collisionality | `_electron_collisionality_sauter()`

The origins of the coefficients values are not known, but thought to be derived from a condition of the [Bohm diffusion coefficient](https://en.wikipedia.org/wiki/Bohm_diffusion)

Using the electron collision frequency ($\nu_{\text{e}}$) calculated from [`_electron_collisions_sauter()`](#calculate-frequency-of-electron-collisions--_electron_collisions_sauter) we get:
$$
\nu_{\text{e*}} = \frac{1.4 \ R \ \nu_{\text{e}}  \ Z_{\text{eff}}}{\left|\frac{1}{q}\epsilon^{3/2}\sqrt{T_{\text{e}}}\times 1.875\times10^7\right|}
$$

-------------

#### Calculate frequency of ion collisions | `_ion_collisions_sauter()`


$$
\nu_{\text{i}} = 320 \times \frac{Z_{\text{eff}}^4n_{\text{i}}}{T_{\text{i}}^{3/2}\sqrt{a_{\text{i}}}}
$$

-----

#### Calculate ion collisionality | `_ion_collisionality_sauter()`

The origins of the coefficients values are not known, but thought to be derived from a condition of the [Bohm diffusion coefficient](https://en.wikipedia.org/wiki/Bohm_diffusion)

Using the ion collision frequency ($\nu_{\text{i}}$) calculated from [`_ion_collisions_sauter()`](#calculate-frequency-of-ion-collisions--_ion_collisions_sauter) we get:

$$
\nu_{\text{e*}} = \frac{3.2\times10^{-6} \nu_{\text{i}} R}{\left|\left(\frac{1}{q}+0.0001\right)\epsilon^{3/2} \sqrt{\frac{T_{\text{i}}}{a_{\text{i}}}} \right|}
$$

----------------

#### Calculate electron only poloidal beta correction | `_beta_poloidal_sauter()`

This function returns an electron only local poloidal beta correction dependant on the array index of the profile.

If the current index is not equal to the size value (or end of the array) then the following is returned:

$$
\frac{1.6\times 10^{-4}\pi R \left(n_{\text{e}}+n_{\text{e-1}}\right)\times \left(T_{\text{e}}+T_{\text{e-1}}\right)}{\left(B_{\text{T}}\rho_{-1}\left|\left(\frac{1}{q}\right)_{-1}+1\times 10^{-4}\right|\right)^2}
$$

Otherwise the following is returned:

$$
\frac{6.4\times 10^{-4}\pi R \left(n_{\text{e-1}}T_{\text{e-1}}\right)}{\left(B_{\text{T}}\rho_{-1}\left|\left(\frac{1}{q}\right)_{-1}+1\times 10^{-4}\right|\right)^2}
$$

The $-1$ subscript in this case refers to the value of the variable in the previous array index value



---------------

#### Calculate ion and electron poloidal beta correction | `_beta_poloidal_total_sauter()`

This function returns the local poloidal beta correction with both electron and ion pressure dependant on the array index of the profile.

If the current index is not equal to the size value (or end of the array) then the following is returned:

$$
\frac{1.6\times 10^{-4}\pi R \left[\left(\left(n_{\text{e}}+n_{\text{e-1}}\right)\times \left(T_{\text{e}}+T_{\text{e-1}}\right)\right)+\left(\left(n_{\text{i}}+n_{\text{i-1}}\right)\times \left(T_{\text{i}}+T_{\text{i-1}}\right)\right)\right]}{\left(B_{\text{T}}\rho_{-1}\left|\left(\frac{1}{q}\right)_{-1}+1\times 10^{-4}\right|\right)^2}
$$

Otherwise the following is returned:

$$
\frac{6.4\times 10^{-4}\pi R \left[\left(n_{\text{e-1}}T_{\text{e-1}}\right)+\left(n_{\text{i-1}}T_{\text{i-1}}\right)\right]}{\left(B_{\text{T}}\rho_{-1}\left|\left(\frac{1}{q}\right)_{-1}+1\times 10^{-4}\right|\right)^2}
$$

The $-1$ subscript in this case refers to the value of the variable in the previous array index value


------------------

### Sakai Scaling | `bootstrap_fraction_sakai()`

Is selected by setting `i_bootstrap_current = 5`[^6]

$$
f_{\text{BS}} = 10^{0.951 \epsilon - 0.948} \cdot \beta_p^{1.226 \epsilon + 1.584} \cdot l_i^{-0.184\epsilon - 0.282} \cdot \left(\frac{q_{95}}{q_0}\right)^{-0.042 \epsilon - 0.02} \\
\cdot \alpha_n^{0.13 \epsilon + 0.05} \cdot \alpha_t^{0.502 \epsilon - 0.273}
$$

The model includes the toroidal diamagnetic current in the calculation due to the dataset, so `i_diamagnetic_current = 0` can only be used with it

-------------------

### ARIES Scaling | `bootstrap_fraction_aries()`

Is selected by setting `i_bootstrap_current = 6`[^8]

The source reference[^8] does not provide any info about the derivation of the formula. It is only stated like that shown below.

$$
a_1 = 1.10-1.165l_{\text{i}}+0.47l_{\text{i}}^2
$$

$$
b_1 = 0.806-0.885l_{\text{i}}+0.297l_{\text{i}}^2
$$

$$
C_{\text{BS}} = a_1+b_1\left(\frac{n(0)}{\langle n \rangle}\right)
$$

$$
f_{\text{BS}} = C_{\text{BS}} \sqrt{\epsilon}\beta_{\text{p}}
$$

---------------------

### Andrade Scaling 

Is selected by setting `i_bootstrap_current = 7`[^9]

Based off of 350 plasma profiles from Experimento Tokamak Esferico (ETE) spherical tokamak.
Profiles were taken to be Gaussian shaped functions.

The range of parameters from the discharges in ETE can be found in the table below:

| Parameter             | Value     | Parameter             | Value     |
|-----------------------|-----------|-----------------------|-----------|
| Aspect ratio, $A$                   | 1.5       | Electron temperature profile index, $\alpha_{\text{T}_\text{e}}$ | 0.02      |
| Major radius, $R_0$                 | 0.3 $\text{[m]}$     | Ion temperature profile index, $\alpha_{\text{T}_\text{i}}$ | 2         |
| Core plasma pressure, $p(0)$                | 15 $\text{[kPa]}$    | Elongation, $\kappa(a)$           | 2         |
| Core electron/ion temperature, $T_{\text{e,i}}(0)$   | 1 $\text{[keV]}$     | Triangularity, $\delta$              | 0.3       |
| Edge electron/ion temperature $T_{\text{e,i}}(a)$   | 0.1 $\text{[keV]}$   | Plasma Current, $I_{\text{p}}$        | 200 $\text{[kA]}$    |
| Pressure profile index, $\alpha_{\text{p}}$   | 3         | Toroidal field on-axis, $B_0$                 | 0.4 $\text{[T]}$     |
| Plasma total beta, $\beta$               | 4-10%     | Plasma effective charge, $Z_{\text{eff}}$      | 1         |

Errors mostly up to the order of 10% are obtained when both expressions are compared with the equilibrium estimates for the bootstrap current in ETE.

#### Scaling 1

!!! note "Applicability of 1st scaling"

    Andrade et.al[^9] actually presents two scalings, the first is:

    $$
    \frac{I_{\text{BS}}}{I_\text{p}}5C_{\text{bs}}c_{\text{p}}^{\lambda}\frac{\beta_{\text{N}}q_{\text{cyl}}}{\epsilon^{1/2}l_i^{\gamma}}\frac{R_0}{R_{\text{m}}}
    $$

    $R_{\text{m}}$ in this case is the major radius of the magnetic axis which `PROCESS` does not have as it does not currently calculate the [Shafranov shift](https://wiki.fusion.ciemat.es/wiki/Shafranov_shift). If this is implemented then the first Andrade scaling can be properly implemented.

#### Scaling 2 | `bootstrap_fraction_andrade()`

This form of the Andrade scaling in terms of $\beta_{\text{p}}$ found using a least-square fit to 347 of the equilibria points.

$$
C_{\text{BS}} = 0.2340 \pm 0.0007
$$

$$
c_p = \frac{p_0}{\langle p \rangle}
$$

$$
f_{\text{BS}} = C_{\text{BS}} \sqrt{\epsilon}\beta_{\text{p}}c_p^{0.8}
$$

---------------------

### Hoang Scaling | `bootstrap_fraction_hoang()`

Is selected by setting `i_bootstrap_current = 8`[^10]

This scaling is based off of 170 discharges from TFTR with the neoclassical bootstrap current being calculated by the TRANSP plasma analysis code.
The plasma parameters of the discharges can be seen in the table below:

| Parameter             | Value     |
|-----------------------|-----------|
| Plasma Current, $I_{\text{p}}$        | 0.6 - 2.7 $\text{[MA]}$    |  
| Edge safety factor, $q(a)$                 | 2.8 - 11.0       |
| Injected NBI power, $P_{\text{NBI}}$                | 2.0 - 35.0 $\text{[MW]}$    |
| Injected ICRH power, $P_{\text{ICRH}}$   | 1.5 - 6.0 $\text{[MW]}$     |
| Toroidal field on-axis, $B_{\text{T}}$   | 1.9 - 5.7 $\text{[T]}$   |
| Core electron density, $n_{\text{e,0}}$   | 0.2 - 1.2 $[10^{20} \text{m}^{-3}]$         |

A wide variety of discharge regimes are included, such as: L-mode  supershots, discharges
with reversed shear (RS) and enhanced reversed shear (ERS), and discharges with increased-$l_i$.
Discharges with both monotonic $q$ profiles and with reversed shear are included in the dataset.

For an example ERS discharge the bootstrap current is driven by when the thermal particles surpasses 1 MA ($f_{\text{boot}}$ = 63%). Some ERS discharges in TFTR achieved $f_{\text{boot}}$ greater than 100% transiently.


!!! note "Change of profile index definition"

    Hoang et.al uses a different definition for the profile indexes such that
    $\alpha_{\text{p}}$ is defined as the ratio of the central and the volume-averaged values, and the peakedness of the density of the total plasma current (defined as ratio of the central value and $I_{\text{p}}$), $\alpha_{\text{J}}$ 

    Assuming that the pressure and current profile is parabolic we can represent these ratios as $\frac{p_0}{\langle p \rangle}= \alpha_{\text{p}}+1$

    **This could lead to large changes in the value depending on interpretation of the profile index**

$$
C_{\text{BS}} = \sqrt{\frac{\alpha_{\text{p}}+1}{\alpha_{\text{j}}+1}}
$$

$$
f_{\text{BS}} = 0.4C_{\text{BS}} \sqrt{\epsilon}\beta_{\text{p}}^{0.9}
$$

---------------------

### Wong Scaling | `bootstrap_fraction_wong()`

Is selected by setting `i_bootstrap_current = 9`[^11]

This scaling data is based off of equilibria from Miller et.al.[^12].
The equilibria from Miller et.al. are in the range of $A$ =  1.2 - 3 that are stable to infinite $n$ ballooning and low $n$ kink modes at a bootstrap fraction of 99% for $\kappa$ = 2, 2.5, 3.0. The results were parameterized as a function of aspect ratio and elongation.

The parametric dependency of $\beta_{\text{p}}$ and $\beta_{\text{T}}$ are based on the fitting of the DIII-D high equivalent DT yield results.

$$
\beta_{\text{N}}=\frac{\left(3.09+\frac{3.35}{A}+\frac{3.87}{A^{0.5}}\right)\left(\frac{\kappa}{3}\right)^{0.5}}{f_{\text{peak}}^{0.5}}
$$

$$
\beta_{\text{T}}=\frac{25}{\beta_p}\left(\frac{1+\kappa^2}{2}\right)\left(\frac{\beta_N}{100}\right)^2
$$

Here $\beta_{\text{p}}$ is given by

$$
\beta_p=f_{b s} \frac{\sqrt{A}}{C_{b s} f_{\text{peak}}^{0.25}}
$$

$$
C_{\text{BS}} = 0.773+0.019\kappa
$$

Parabolic profiles should be used for best results as the pressure peaking value is calculated as the product of a parabolic temperature and density profile.

$$
f_{\text{peak}} = \left(\int_0^1 \left(1-\rho^2 \right)^{\alpha_{\text{T}}} \left(1-\rho^2 \right)^{\alpha_{\text{n}}} \ \ \mathrm{d\rho}\right)^{-1}
$$

The integral above is set by the definite solution below

$$
\frac{\operatorname{B}\left(\frac{1}{2},{\alpha}_{n} + {\alpha}_{T} + 1\right)}{2}
$$

Assuming that $\alpha_{\text{n}} + \alpha_{\text{T}} > 0,   \alpha_{\text{n}} + \alpha_{\text{T}} + 1 > 0$

$$
f_{\text{BS}} = C_{\text{BS}} f_{\text{peak}}^{0.25}\beta_{\text{p}}\sqrt{\epsilon}
$$

---------------------

### Gi Scaling's

This scaling is found by solving the Hirshman-Sigmar bootstrap current model using the matrix inversion method to create bootstrap current scalings with variables given explicitly in the TPC systems code[^13].
A 8800 point database for the bootstrap current fraction using the bootstrap current density calculation module in the ACCOME code is used, using the variable ranges in the table below:

The fitting of the variable exponents is done using the least squares method with a $R^2$ value of > 0.98 for [scaling one](#scaling-1--bootstrap_fraction_gi_i) and > 0.96 for [scaling two](#scaling-2--bootstrap_fraction_gi_ii) compared to the ACCOME data.


| Parameter                  | Range          | Points |
|----------------------------|----------------|--------|
| Major radius, $R$  | 5.0 $[\mathrm{m}]$   | 1      |
| Aspect ratio, $A$        | 1.3, 1.5, 1.7, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 5.0 | 10     |
| Elongation, $\kappa$   | $\sim$ 2 | 1      |
| Triangularity, $\delta$   | $\sim$ 0.3 | 1      |
| Density profile index, $a_{\text{n}}$      | 0.1-0.8 | 8      |
| Temperature profile index, $a_{\text{T}}$  | 1.0-3.0 | 11     |
| Effective charge, $Z_{\text{eff}}$ | 1.2-3.0 | 10     |

The plasma parameters for each point in the aspect ratio scan can be seen in the table below:

| Aspect ratio A | 1.3 | 1.5 | 1.7 | 2.0 | 2.2 | 2.5 | 3.0 | 3.5 | 4.0 | 5.0 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Electron density at axis, $n_{\text{e0}}\left[10^{20} \mathrm{~m}^{-3}\right]$ | 1.0 | 1.0 | 1.5 | 1.5 | 1.5 | 2.0 | 2.0 | 2.0 | 2.0 | 1.0 |
| Electron temperature at axis, $T_{\text{e0}}[\mathrm{keV}]$ | 40 | 20 | 30 | 30 | 30 | 40 | 20 | 40 | 20 | 30 |
| Plasma current, $I_p$ $[\mathrm{MA}]$ | 20 | 15 | 20 | 15 | 15 | 15 | 10 | 15 | 10 | 5 |
| Toroidal magnetic field at axis, $B_{\text{T}}$ $[\mathrm{T}]$ | 3.0 | 2.0 | 3.0 | 2.0 | 4.0 | 2.0 | 2.0 | 6.0 | 2.0 | 5.0 |
| Poloidal beta, $\beta_{\text{p}}$ | 0.9-2.6 | 0.6-1.8 | 0.6-1.8 | 0.8-1.9 | 0.7-2.0 | 0.9-2.7 | 0.7-2.2 | 0.5-1.4 | 0.4-1.2 | 0.8-2.3 |



#### Scaling 1 | `bootstrap_fraction_gi_I()`

Is selected by setting `i_bootstrap_current = 10`

Scaling 1 has better accuracy than Scaling 2. However, Scaling 1 overestimated the $f_{\text{BS}}$ value for reversed shear equilibrium. Although Scaling 2 does not have an internal current profile term, it can predict the $f_{\text{BS}}$ values to a certain extent for the high-$f_{\text{BS}}$ equilibria which are expected for next fusion devices.

$$
C_{\text{BS}} = 0.474 \epsilon^{-0.1} \alpha_{\text{p}}^{0.974} \alpha_{\text{T}}^{-0.416} Z_{\text{eff}}^{0.178} \left(\frac{q_{95}}{q_0}\right)^{-0.133}
$$

$$
f_{\text{BS}} = C_{\text{BS}} \beta_{\text{p}}\sqrt{\epsilon}
$$

#### Scaling 2 | `bootstrap_fraction_gi_II()`

Is selected by setting `i_bootstrap_current = 11`

This scaling has the $q$ profile dependance removed to obtain a scaling formula with much more flexible variables than that by a single profile factor for internal current profile.

$$
C_{\text{BS}} = 0.382 \epsilon^{-0.242} \alpha_{\text{p}}^{0.974} \alpha_{\text{T}}^{-0.416} Z_{\text{eff}}^{0.178}
$$

$$
f_{\text{BS}} = C_{\text{BS}} \beta_{\text{p}}\sqrt{\epsilon}
$$

---------------------

### Sugiyama Scaling's

This scaling is found by solving the Hirshman-Sigmar bootstrap current model using the matrix inversion method to create bootstrap current scalings with variables given explicitly in the TPC systems code[^14]. The databases are constructed with the ACCOME code, with the ohmic and externally driven currents not considered and the total plasma current is achieved by
introducing a virtual additional current of the shape:

$$
j_{\text{add}}(\rho) \propto \left(1-\rho^{b_{\text{j}}}\right)^{a_{\text{j}}}
$$

The parameters for the ACCOME database are as follows, with the scan range in the table below:

- $R$ = 6 m
- $\kappa \sim 1.8$ 
- $\delta \sim 0.4$
- $Z_{\text{eff}}$ has a uniform profile, and only fully stripped carbon is considered as an impurity.
- $T_{\text{e}} = T_{\text{i}}$

| Parameter                              | Range                                                                 |
|----------------------------------------|----------------------------------------------------------------------|
| Aspect ratio, $A$                      | 1.5, 1.7, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0                    |
| Plasma Current, $I_{\text{p}} [\mathrm{MA}]$ | 30.0, 27.0, 22.0, 20.0, 17.5, 13.0, 11.5, 10.0, 8.8, 8.0             |
| Toroidal magnetic field on axis, $B_{\text{T}} [\mathrm{T}]$ | 3.5, 3.5, 3.5, 4.1, 4.5, 4.5, 5.2, 5.67, 6.8, 8.5                  |

---------------

#### L-mode Scaling | `bootstrap_fraction_sugiyama_l_mode()`

Is selected by setting `i_bootstrap_current = 12`


| Parameter                              | Range                                                                 |
|----------------------------------------|----------------------------------------------------------------------|
| Aspect ratio, $A$                      | 1.5, 1.7, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0                    |
| $\alpha_{\text{n}}$ | 0.1, 0.3, 0.5, 0.7, 0.3, 1.1            |
| $\alpha_{\text{T}}$ | 1.0, 1.4, 1.8, 2.2, 2.6, 3.0            |
| $Z_{\text{eff}}$ | 1.2, 1.5, 2.0, 2.5, 3.0           |
| $n_{\text{e,0}} \  [10^{20} \text{m}^{-3}]$  | 0.5, 1.0, 2.0                  |
| $T_{\text{e,0}} \  [\text{keV}]$  | 20, 30, 30                  |
| $a_{\text{j}}$  | 0.2, 0.5, 1.0, 1.1, 2.0, 2.0, 4.0, 10.0, 50.0                 |
| $b_{\text{j}}$  | 2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 4.0, 4.0, 4.0                  |

The scan range for the L-mode scaling is seen in the table above and for the variables in the top table which are shared between both scalings.  For all combinations in both tables there are 48,600 cases of which 47,652 converged. The mean error of the scaling with the ACCOME data is, $\text{ME} = 5.92 \times 10^{-4}$, with the root mean squared error being,  $\text{RMSE} = 0.0236$.

$$
f_{\text{BS}}^{\text{L}} = 0.740 \epsilon^{0.418} \beta_{\text{p}}^{0.904} \alpha_{\text{n}}^{0.06} \alpha_{\text{T}}^{-0.138} Z_{\text{eff}}^{0.230} \left(\frac{q_{95}}{q_0}\right)^{-0.142}
$$

------------

#### H-mode Scaling | `bootstrap_fraction_sugiyama_h_mode()`

Is selected by setting `i_bootstrap_current = 13`

| Parameter                              | Range                                                                 |
|----------------------------------------|----------------------------------------------------------------------|
| Aspect ratio, $A$                      | 1.5, 1.7, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0                    |
| $\alpha_{\text{n}}$ | 0.8, 1.4, 2.2, 3.0           |
| $\alpha_{\text{T}}$ | 1.2, 1.5, 3.0, 6.0            |
| $b_{\text{T}}$ | 1.2, 1.5, 3.0, 6.0            |
| $Z_{\text{eff}}$ | 1.2, 1.5, 2.0, 3.0           |
| $n_{\text{e,0}} \  [10^{20} \text{m}^{-3}]$ | 1.5 |
| $T_{\text{e,0}} \  [\text{keV}]$  | 30                  |
| $n_{\text{e,ped}} \  [10^{20} \text{m}^{-3}]$  | 0.6, 1.05, 1.5                 |
| $T_{\text{ped}} \ [\text{keV}]$  | 0.5, 2.0, 4.0, 8.0 |
| $\rho_{\text{ped}}$  | 0.85, 0.92, 0.99 |
| $a_{\text{j}}$  | 0.6, 1.2, 1.5, 1.1, 10.0, 50.0                 |
| $b_{\text{j}}$  | 2.0, 1.5, 1.0, 4.0, 4.0                |

The scan range for the H-mode scaling is seen in the table above and for the variables in the top table which are shared between both scalings.  For all combinations in both tables there are 460,800 cases of which 330,149 converged. The mean error of the scaling with the ACCOME data is, $\text{ME} = 1.17 \times 10^{-4}$, with the root mean squared error being,  $\text{RMSE} = 0.0324$.

It is also assumed that:

- $n_{\text{sep}} = 0 \ \text{m}^{-3}$
- $T_{\text{sep}} = 0 \ \text{keV}$
- $\rho_{\text{ped}} = \rho_{\text{ped,n}} = \rho_{\text{ped,T}}$

$$
f_{\text{BS}}^{\text{H}} = 0.789 \epsilon^{0.606} \beta_{\text{p}}^{0.960} \alpha_{\text{n}}^{0.0319} \alpha_{\text{T}}^{0.00822} b_{\text{T}}^{-0.0783} Z_{\text{eff}}^{0.241} \\ \times \left(\frac{q_{95}}{q_0}\right)^{-0.103} \rho_{\text{ped}}^{0.367} \left(\frac{n_{\text{e,ped}}}{n_{\text{GW}}}\right)^{-0.174} T_{\text{ped,keV}}^{0.0552}
$$

---------------------

## Setting of maximum desirable bootstrap current fraction

The variable `f_c_plasma_bootstrap_max` can be set to the value of maximum desirable bootstrap current fraction for a specific design. When optimising if the value of the calculated `f_c_plasma_bootstrap` for the model selected with `i_bootstrap_current` exceeds this value, then `f_c_plasma_bootstrap` is set to the value of `f_c_plasma_bootstrap_max`.

An error is also raised to the user in the terminal output at the end of the run saying "Bootstrap fraction upper limit enforced".

## Fixing the bootstrap current fraction

If the user wants to set the value of the bootstrap current fraction directly then the value can be set by setting `i_bootstrap_current = 0` and then writing the value directly with 


```txt
>>> IN.DAT

# Setting a fixed bootstrap current fraction of 80%

i_bootstrap_current = 0
f_c_plasma_bootstrap = 0.8
```


[^0]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
[^1]: Nevins, W. M. "Summary report: ITER specialists’ meeting on heating and current drive." ITER-TN-PH-8-4, June 1988. 1988. 
[^2]: Keii Gi, Makoto Nakamura, Kenji Tobita, Yasushi Ono, Bootstrap current fraction scaling for a tokamak reactor design study,
Fusion Engineering and Design, Volume 89, Issue 11, 2014, Pages 2709-2715, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2014.07.009.
[^3]: Wilson, H.R. (1992). Bootstrap current scaling in tokamaks. Nuclear Fusion, 32(2), pp.257–263. doi:https://doi.org/10.1088/0029-5515/32/2/i05.
[^4]: O. Sauter, C. Angioni, Y. R. Lin-Liu; Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime. Phys. Plasmas 1 July 1999; 6 (7): 2834–2839. https://doi.org/10.1063/1.873240 
[^5]: O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum: “Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime” [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140. https://doi.org/10.1063/1.1517052  
[^6]: Ryosuke Sakai, Takaaki Fujita, Atsushi Okamoto, Derivation of bootstrap current fraction scaling formula for 0-D system code analysis, Fusion Engineering and Design, Volume 149, 2019, 111322, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2019.111322.
[^7]: T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992 
[^8]: Zoran Dragojlovic et al., “An advanced computational algorithm for systems analysis of tokamak power plants, ”Fusion Engineering and Design, vol. 85, no. 2, pp. 243–265, Apr. 2010, doi: https://doi.org/10.1016/j.fusengdes.2010.02.015.
[^9]: M. C. R. Andrade and G. O. Ludwig, “Scaling of bootstrap current on equilibrium and plasma profile parameters in tokamak plasmas,” Plasma Physics and Controlled Fusion, vol. 50, no. 6, pp. 065001–065001, Apr. 2008, doi: https://doi.org/10.1088/0741-3335/50/6/065001.
[^10]: G. T. Hoang and R. V. Budny, “The bootstrap fraction in TFTR,” AIP conference proceedings, Jan. 1997, doi: https://doi.org/10.1063/1.53414.
[^11]: C.-P. Wong, J. C. Wesley, R. D. Stambaugh, and E. T. Cheng, “Toroidal reactor designs as a function of aspect ratio and elongation,” vol. 42, no. 5, pp. 547–556, May 2002, doi: https://doi.org/10.1088/0029-5515/42/5/307.
[^12]: Miller, R L, "Stable bootstrap-current driven equilibria for low aspect ratio tokamaks". Switzerland: N. p., 1996. Web.https://fusion.gat.com/pubs-ext/MISCONF96/A22433.pdf
[^13]: K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, “Bootstrap current fraction scaling for a tokamak reactor design study,” Fusion Engineering and Design, vol. 89, no. 11, pp. 2709–2715, Aug. 2014, doi: https://doi.org/10.1016/j.fusengdes.2014.07.009.
[^14]: S. Sugiyama, T. Goto, H. Utoh, and Y. Sakamoto, “Improvement of core plasma power and current balance models for tokamak systems code considering H-mode plasma profiles,” Fusion Engineering and Design, vol. 216, p. 115022, Jul. 2025, doi: https://doi.org/10.1016/j.fusengdes.2025.115022.
‌ 