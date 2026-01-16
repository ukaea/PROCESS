# Density Limit

Several density limit models are available in PROCESS. These are
calculated in routine `calculate_density_limit()`, which is called by `physics`.

This constraint can be activated by stating `icc = 5` in the input file.

The value of `i_density_limit` can be set to apply the relevant limit. The variable `fdene` can be set to scale the constraint bound: `nd_plasma_electrons_vol_avg` / `nd_plasma_electrons_max` <= `fdene` (or `nd_plasma_electron_line` / `nd_plasma_electrons_max` <= `fdene` if `i_density_limit = 7`).

For the `i_density_limit = 1-5,8` scalings we scale the function output by the separatrix to volume averaged electron density so that we can set the limit on the volume averaged. **Therefore it is recommended to only use these scalings with an H-mode profile (`i_plasma_pedestal == 1`) otherwise the separatrix density (`nd_plasma_separatrix_electron`) will not be calculated.**

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

Switch value: `i_density_limit = 7` [^3][^4]

$$
\overline{n}_{\text{e}}^{\text{ crit}} = 1.0 \times 10^{14} \frac{I_\text{p}}{\pi a^2}
$$

For the Greenwald model the limit applies to the line-averaged electron density, not the volume-averaged density. The plasma current term is given in $[\mathrm{A}]$ and the minor radius in $[\mathrm{m}]$

---------------------

## ASDEX New model

Switch value: `i_density_limit = 8` [^5][^6]

$$
\overline{n}_{\text{sep}}^{\text{ crit}} = 1.0 \times 10^{20} \times 0.506 \pm 0.192 \frac{P_\text{heat}^{0.396\pm0.13} I_{\text{p}}^{0.265\pm 0.14}}{q_{95}^{0.323 \pm 0.14}}
$$

-----------------

[^1]: T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',

[^3]: M. Greenwald et al., “A new look at density limits in tokamaks,” Nuclear Fusion, vol. 28, no. 12, pp. 2199–2207, Dec. 1988, doi: https://doi.org/10.1088/0029-5515/28/12/009.

[^4]: M. Greenwald, “Density limits in toroidal plasmas,” Plasma Physics and Controlled Fusion, vol. 44, no. 8, pp. R27–R53, Jul. 2002, doi: https://doi.org/10.1088/0741-3335/44/8/201.

[^5]: J. W. Berkery et al., “Density limits as disruption forecasters for spherical tokamaks,” Plasma Physics and Controlled Fusion, vol. 65, no. 9, pp. 095003–095003, Jul. 2023, doi: https://doi.org/10.1088/1361-6587/ace476.

[^6]: M. Bernert et al., “The H-mode density limit in the full tungsten ASDEX Upgrade tokamak,” vol. 57, no. 1, pp. 014038–014038, Nov. 2014, doi: https://doi.org/10.1088/0741-3335/57/1/014038.
‌
‌