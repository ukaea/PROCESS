# Scrape off Layer | `ScrapeOffLayer`

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

[^eich_2013]: T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9 p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

[^mast_2014]: A. J. Thornton and A. Kirk, “Scaling of the scrape-off layer width during inter-ELM H modes on MAST as measured by infrared thermography,”
Plasma Physics and Controlled Fusion, vol. 56, no. 5, p. 055008, Apr. 2014, doi: 10.1088/0741-3335/56/5/055008.