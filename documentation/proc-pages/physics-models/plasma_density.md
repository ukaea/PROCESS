# Density Limit

Several density limit models[^1][^2] are available in PROCESS. These are
calculated in routine `calculate_density_limit`, which is called by `physics`. To enforce any of 
these limits, turn on constraint equation no. 5 with iteration variable no. 9 
(`fdene`). In addition, switch `i_density_limit` must be set to the relevant value, as 
follows:

For the `i_density_limit = 1-5` scalings we scale the function output by the separatrix to volume averaged electron density so that we can set the limit on the volume averaged. Therefore it is recommended to only use these scalings with an H-mode profile (`ipedestal == 1`) otherwise the separatrix density (`nesep`) will not be calculated.

For the models below $P_{\perp}$ is the mean heat flux density across the separatrix ($\mathrm{MW}/\mathrm{m^2}$), which we take as the divertor power divided by the plasma surface area.

-----------------

## ASDEX model

Switch value: `i_density_limit = 1`[^1][^2]

$$
n_{\text{b}}^{\text{crit}} = 1.54 \frac{P_{\perp}^{0.43}B_{\text{T}}^{0.31}}{\left(q_{95}R\right)^{0.45}}
$$

-----------------

## Borrass model for ITER, I

Switch value: `i_density_limit = 2` [^1]

$$
n_{\text{b}}^{\text{crit}} = C \frac{P_{\perp}^{0.53}B_{\text{T}}^{0.31}}{\left(q_{95}R\right)^{0.22}}
$$

$C \approx$  1.8 for ITER-like conditions.

-----------------

## Borrass model for ITER, II 

Switch value: `i_density_limit = 3` [^1]

$$
n_{\text{b}}^{\text{crit}} = 0.5 \frac{P_{\perp}^{0.57}B_{\text{T}}^{0.31}}{\left(q_{95}R\right)^{0.09}}
$$

-----------------

## JET edge radiation model

Switch value: `i_density_limit = 4` [^1]

$$
n_{\text{b}}^{\text{crit}} = P_{\text{in}}^{0.5} \frac{1}{\left[\left(Z_{\text{eff}}-1\right)\left(1-\frac{4}{3q_{\text{c}}}\right)\right]^{0.5}}
$$

-----------------

## JET simplified model

Switch value: `i_density_limit = 5` [^1]

$$
n_{\text{b}}^{\text{crit}} = 0.237 P^{0.5}
$$

For a radiation from a shell thickness $\Delta$, this may be written as:

$$
n_{\text{b}}^{\text{crit}} = 0.147 \frac{P^{0.5}}{\left[Ra\Delta \sqrt{\frac{\left(1+\kappa \right)}{2}}\right]^{0.5}}
$$

where $\kappa \approx 1.5, \Delta \approx 0.1a$ has been taken from JET.

-----------------

## Hugill-Murakami model

Switch value: `i_density_limit = 6` [^2]

$$
\langle n_{^{\text{crit}}} \rangle \approx \frac{3.0 B_{\text{T}}}{R_0 q_{\text{cyl}}}
$$


-----------------

## Greenwald model

Switch value: `i_density_limit = 7`

$$
\overline{n}_{\text{e}}^{\text{ crit}} = 1.0 \times 10^{14} \frac{I_\text{p}}{\pi a^2}
$$

For the Greenwald model the limit applies to the line-averaged electron density, not the volume-averaged density.

-----------------

## Key Constraints

[^1]: T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',