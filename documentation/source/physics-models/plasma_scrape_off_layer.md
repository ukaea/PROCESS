# Scrape off Layer | `ScrapeOffLayer`


## Parallel energy flux density

The peak unmitigated parallel upstream energy flux density, $q_{\parallel,u}$ $[\text{W} \text{m}^{-2}]$ is given by[^stangeby_boundary] [^henderson_step]:

$$
q_{\parallel,u} = \frac{P_{\text{up}}}{2\pi\lambda_{\text{q,u}}R_{\text{u}}}\frac{B_{\text{Tot,u}}}{B_{\text{p,u}}}
$$

$$
q_{\parallel,u} = \frac{P_{\text{up}}}{A_{\parallel,u}}
$$

$P_{\text{up}}$ is defined as $P_{\text{sep}}f_{\text{div}}$, where $P_{\text{sep}}$ is the plasma separatrix power and $f_{\text{div}}$ is the fraction of the power directed to the divertor. $\lambda_{\text{q,u}}$ is the [upstream power decay length](#power-decay-lengths), $R_{\text{u}}$ is the midplane separatrix radius, $B_{\text{Tot,u}}$ is the total magnetic field strength at the midplane separatrix and $B_{\text{p,u}}$ is the poloidal magnetic field strength at the midplane separatrix.

------------------

### Parallel area | `calculate_upstream_sol_outboard_parallel_area()`

The outboard parallel SOL flux area $A_{\parallel,u}$ [$\text{m}^2$] can be envisaged as a ribbon going around the outboard circumference of the plasma of width $\lambda_{\text{q,u}}$. It is then weighted by the field ratio $\frac{B_{\text{p,u}}}{B_{\text{Tot,u}}}$ to account for the field line angle at the outboard midplane:

$$
A_{\parallel,u} = 2\pi\lambda_{\text{q,u}}R_{\text{u}}\frac{B_{\text{p,u}}}{B_{\text{Tot,u}}}
$$

---------------

## Upstream radial decay | `calculate_outboard_midplane_near_sol_radial_profile()`

The radial decay length $\lambda_{\text{q}}$ of the scrape-off layer (SOL) at the outer midplane of a tokamak is defined as the e-folding distance over which plasma heat and particle fluxes decay exponentially outside the last closed flux surface. Therefore the decay of the total heat flux outside the separatrix towards the vessel walls can be modelled as [^eich_2011] [^eich_2013]:

$$
q_{\text{u}}(r) = q_{\parallel,\text{u}}e^{\frac{-r}{\lambda_{\text{q}}}}
$$

where $r = R - R_{\text{sep}}$, $R_{\text{sep}}$ being the major radius of the separatrix, $\lambda_{\text{q}}$ the [power decay length](#power-decay-lengths) and $q_{\parallel}$ the [upstream energy flux density](#upstream-radial-decay--calculate_outboard_midplane_near_sol_radial_profile)

----------------

## Eich parallel flux at target | `calculate_eich_target_heat_flux_profile()`

The Eich formula (often called the standard SOL heat flux profile) is the primary mathematical model used to describe the distribution of heat target loads on tokamak divertor plates. It convolutionally connects the physics of the plasma edge at the outer midplane with the geometric projection of the heat hitting the divertor surface [^eich_2011] [^eich_2013].

Heat transport into the private flux region is modeled by convolving the power profile $q_{\text{u}}(r)$ with a Gaussian function of width $S$ known as the [spreading parameter](#spreading-parameter). 

$$
q_{\parallel,t} = \frac{q_0}{2}\times \exp\left(\left(\frac{S}{2\lambda_{\text{q}}}\right)- \frac{\overline{s}}{\lambda_q f_x}\right) \times \operatorname{erfc}\left(\frac{S}{2\lambda_{\text{q}}}- \frac{\overline{s}}{S f_{x}}\right) + q_{\text{BG}}
$$

where $\overline{s} = s- s_0 = (R_{\text{sep}} - R) \times f_x $. $\operatorname{erfc}$ is the complementary error function, $q_{\text{BG}}$ is the background heat flux, $\lambda_{\text{q}}$ is the [power decay length](#power-decay-lengths), $f_x$ is the effective flux expansion in the region,

------------------

## Power Decay Lengths

In tokamak scrape-off layer (SOL) physics, $\lambda_q$ is the radial heat flux e-folding length, representing the narrow width over which exhaust power drops exponentially. It measures how tightly thermal energy exiting the core plasma is compressed onto open magnetic field lines.

-----------

### Eich 2013 Model | `calculate_eich2013_sol_power_decay_length()`

The power decay length in metres is given by[^eich_2013]:

$$
\lambda_q = 1.35 \times 10^{-3} P_{\text{sep}}^{-0.02}R_0^{0.04}B_{\text{p}}(a)^{-0.92}\epsilon^{0.42}
$$

Here $P_{\text{sep}}$ is the plasma separatrix power in $\text{MW}$, $R_0$ is the plasma major radius in $\text{m}$, $B_{\text{p}}(a)$ is the plasma poloidal field at the outboard mid-plane measured in $\text{T}$ and $\epsilon$ is the plasma inverse aspect ratio.

This can be found in Table 3 from Eich et.al [^eich_2013]

The $R^2$ value for this fit is 0.88

----------

### MAST 2014 I | `calculate_mast2014_sol_power_decay_length_1()`

The power decay length in metres is given by[^mast_2014]:

$$
\lambda_q = 1.84(\pm0.48) \times 10^{-3} P_{\text{sep}}^{0.18(\pm0.07)}B_{\text{p}}(a)^{-0.68(\pm0.14)}
$$

Here $P_{\text{sep}}$ is the plasma separatrix power in $\text{MW}$, and $B_{\text{p}}(a)$ is the plasma poloidal field at the outboard mid-plane measured in $\text{T}$.

This can be found in Table 2 and Equation 3 from Thornton et.al [^mast_2014]

The $R^2$ value for this fit is 0.56

----------

### MAST 2014 II | `calculate_mast2014_sol_power_decay_length_2()`

The power decay length in metres is given by[^mast_2014]:

$$
\lambda_q = 4.57(\pm0.54) \times 10^{-3} P_{\text{sep}}^{0.22(\pm0.08)}I_{\text{p}}^{-0.64(\pm0.15)}
$$

Here $P_{\text{sep}}$ is the plasma separatrix power in $\text{MW}$, and $I_{\text{p}}$ is the total plasma current measured in $\text{MA}$.

This can be found in Table 2 and Equation 4 from Thornton et.al [^mast_2014]

The $R^2$ value for this fit is 0.55

--------

## Spreading Parameter

The scrape-off layer (SOL) spreading parameter $S$ represents a Gaussian width that quantifies additional perpendicular heat spreading in the divertor leg. It works alongside the upstream heat flux decay length $\lambda_{q}$ to determine total target heat loads on the divertor.

Unlike $\lambda_{q}$, which is governed by robust upstream parallel and perpendicular transport physics at the plasma midplane, $S$ is inherently a "local" divertor parameter. Deriving a single, absolute multi-machine formula for $S$ is incredibly difficult due to several overlapping regional variables:

- Divertor Geometry: The path length from the X-point to the target tile heavily impacts how much the heat spreads radially.

- Plasma Recycling Regimes: Low-recycling, high-recycling, and detached plasma conditions completely alter the cross-field diffusion rates.

- Localized Radiation: Impurity seeding and neutral gas interactions dissipate power unevenly along the divertor leg, altering the effective Gaussian profile width.

-----------

### Scarabosio 2015 | `calculate_scarabosio2015_power_spreading_factor()`

The H-mode SOL spreading factor, $S$ is given in $\text{m}$ by[^scarabosio_2015]:

$$
S = (0.12(\pm0.07)\times 10^{-3}) P_{\text{sep}}^{0.21(\pm0.11)}R_0^{0.71(\pm0.5)}B_{\text{p}}(a)^{-0.82(\pm0.27)}n_{\text{sep}}^{0.71(\pm0.5)}
$$

- This was fitted from ASDEX Upgrade and JET outer target data
- The $R^2$ value of the regression fit was 0.65

------------

[^eich_2013]: T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9 p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

[^mast_2014]: A. J. Thornton and A. Kirk, “Scaling of the scrape-off layer width during inter-ELM H modes on MAST as measured by infrared thermography,”
Plasma Physics and Controlled Fusion, vol. 56, no. 5, p. 055008, Apr. 2014, doi: 10.1088/0741-3335/56/5/055008.

[^stangeby_boundary]: P. C. Stangeby, “The Plasma Boundary of Magnetic Fusion Devices,” Jan. 2000, doi: 10.1201/9780367801489.

[^henderson_step]: S. S. Henderson et al., “An overview of the STEP divertor design and the simple models driving the plasma exhaust scenario,” Nuclear Fusion, vol. 65, no. 1, pp. 016033–016033, Nov. 2024, doi: 10.1088/1741-4326/ad93e7.

[^scarabosio_2015]: A. Scarabosio et al., “Scaling of the divertor power spreading (S-factor) in open and closed divertor operation in JET and ASDEX Upgrade,” Journal of Nuclear Materials, vol. 463, pp. 49-54, Aug. 2015, doi: 10.1016/j.jnucmat.2014.11.076.

[^eich_2011]: T. Eich, B. Sieglin, A. Scarabosio, W. Fundamenski, R. J. Goldston, and A. Herrmann, “Inter-ELM Power Decay Length for JET and ASDEX Upgrade: Measurement and Comparison with Heuristic Drift-Based Model,” Physical Review Letters, vol. 107, no. 21, Nov. 2011, doi: https://doi.org/10.1103/PhysRevLett.107.215001