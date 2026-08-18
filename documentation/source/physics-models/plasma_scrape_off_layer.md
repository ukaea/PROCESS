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

### Eich 2011 JET Model | `calculate_eich2011_jet_sol_power_decay_length()`

The power decay length in metres is given by[^eich_2011]:

$$
\lambda_q = 0.7(\pm0.26) \times 10^{-3} B_{\text{T,0}}^{-0.84(\pm0.26)} q_{\text{cyl}}^{1.23(\pm0.26)} P_{\text{sep}}^{0.14(\pm0.14)}
$$

Here $B_{\text{T,0}}$ is the on axis toroidal magnetic field, $q_{\text{cyl}}$ is the cylindrical safety factor and $P_{\text{sep}}$ is the plasma separatrix power in $\text{MW}$.

This can be found in Table 2 from Eich et.al [^eich_2011]


----------

### Eich 2011 JET+ASDEX Model | `calculate_eich2011_jet_asdex_sol_power_decay_length()`

The power decay length in metres is given by[^eich_2011]:

$$
\lambda_q = 0.73(\pm0.38) \times 10^{-3} B_{\text{T,0}}^{-0.78(\pm0.25)} q_{\text{cyl}}^{1.2(\pm0.27)} P_{\text{sep}}^{0.1(\pm0.11)}R_0^{0.02(\pm0.2)}
$$

Here $B_{\text{T,0}}$ is the on axis toroidal magnetic field, $q_{\text{cyl}}$ is the cylindrical safety factor, $P_{\text{sep}}$ is the plasma separatrix power in $\text{MW}$ and $R_0$ is the plasma major radius.

This can be found in Table 2 from Eich et.al [^eich_2011]

------------------

[^eich_2013]: T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9 p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

[^mast_2014]: A. J. Thornton and A. Kirk, “Scaling of the scrape-off layer width during inter-ELM H modes on MAST as measured by infrared thermography,”
Plasma Physics and Controlled Fusion, vol. 56, no. 5, p. 055008, Apr. 2014, doi: 10.1088/0741-3335/56/5/055008.

[^stangeby_boundary]: P. C. Stangeby, “The Plasma Boundary of Magnetic Fusion Devices,” Jan. 2000, doi: 10.1201/9780367801489.

[^henderson_step]: S. S. Henderson et al., “An overview of the STEP divertor design and the simple models driving the plasma exhaust scenario,” Nuclear Fusion, vol. 65, no. 1, pp. 016033–016033, Nov. 2024, doi: 10.1088/1741-4326/ad93e7.

[^eich_2011]: T. Eich, B. Sieglin, A. Scarabosio, W. Fundamenski, Robert James Goldston, and A. Herrmann, “Inter-ELM Power Decay Length for JET and ASDEX Upgrade: Measurement and Comparison with Heuristic Drift-Based Model,” Physical Review Letters, vol. 107, no. 21, Nov. 2011, doi: https://doi.org/10.1103/PhysRevLett.107.215001. 