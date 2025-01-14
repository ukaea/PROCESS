# Divertor

The divertor provides a means of removing plasma reaching the scrape-off layer. 
The principal outputs from the code are the divertor heat load, used to 
determine its lifetime, and its peak temperature. The divertor is cooled either 
by gaseous helium or by pressurised water.

Switch `snull` controls the overall plasma configuration. Setting `snull = 0` 
corresponds to an up-down symmetric, double null configuration, while 
`snull = 1` (the default) assumes a single null plasma with the divertor at the 
bottom of the machine. The vertical build and PF coil current scaling 
algorithms take the value of this switch into account, although not the plasma 
geometry at present.

The Harrison-Kukushkin-Hotston divertor model[^1] developed for ITER is available, but is unlikely to be relevant for a reactor.

The divertor heat flux `hldiv` can be calculated or it can be input by the user. Options are selected using the switch `i_hldiv`:

| `i_hldiv` | Description |
| :-: | - |
| 0 | the user inputs the value for `hldiv` |
| 1 | the chamber model (`divtart`) is called to calculate `hldiv` |
| 2 | the Wade heat flux model (`divwade`) is called to calculate `hldiv` |

## Chamber model

!!! Note ""
    `i_hldiv == 1`

The tight aspect ratio tokamak divertor model (`divtart`) calculates the divertor heat flux by 
assuming that the power is evenly spread around the divertor chamber by the action of a gaseous 
target. Each divertor is assumed to be approximately triangular in the R,Z plane.

## Wade Heat Flux Model

!!! Note ""
    `i_hldiv == 2`

A divertor heat flux model is provided in Appendix A.II. of [^2].  This uses the Eich scaling 
[^3] and S-factor [^4] to calculate the SOL width at the outboard divertor, mapped to the midplane:

$$
\lambda_{int} = \lambda_{q,Eich} + 1.64S
$$

where

$$
\lambda_{q,Eich} = 1.35 \, P_{\mathrm{SOL}}^{-0.02} \, R_{o}^{0.04} \, B_{p}^{-0.92} \, \epsilon^{0.42}
$$

$$
S = 0.12(n_{e,mid}/10^{19})^{-0.02} \, P_{\mathrm{SOL}}^{-0.21} \, R_{o}^{0.71} \, B_{p}^{-0.82}.
$$

This is then used to calculate the wetted area in the divertor

$$
A_{wetted} = 2\pi N_{div} R \lambda_{int} F_{exp} \sin(\theta_{div})
$$

where $N_{div}$ is the number of divertors (1 or 2), $F_{exp}$ is the relevant flux expansion, and 
$\theta_{div}$ is the tilt of the separatrix relative to the target in the poloidal plane, and has the form

$$
\theta_{div} = \sin^{-1} [(1+1/\alpha_{div}^{2})\sin\beta_{div}],
$$

where

$$
\alpha_{div} = F_{exp}\alpha_{mid}
$$

$$
\alpha_{mid} = \tan^{-1}\frac{B_{p,mid}}{B_{T,mid}}
$$

where $B_{p,mid}$ and $B_{T,mid}$ are the poloidal and toroidal fields on the outer midplane. The 
parameter $\beta_{div}$ is the angle of incidence between the field line and the target.

The divertor heat flux in $\mathrm{MW}/\mathrm{m^{2}}$ is then 

$$
q_{div} = P_{\mathrm{SOL}}(1-f_{rad,div})/A_{wetted}
$$

where $f_{rad,div}$ is the SOL radiative fraction.

For the purposes of this model, the following are inputs:

- Flux expansion $F_{exp}$  (`flux_exp`, default = 2)  
- Field line angle with respect to divertor target plate (degrees) $\beta_{div}$ (`beta_div`), also 
  available as an iteration variable (170)  
- SOL radiative fraction, $f_{rad,div}$ (`rad_fraction_sol`).

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)

[^2]: M.R. Wade & J.A. Leuer, 'Cost Drivers for a Tokamak-Based Compact Pilot Plant, Fusion Science and Technology, 77:2, 119-143 (2021)

[^3]: T. Eich et al, 'Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER', Nucl. Fusion 53 093031 (2013)

[^4]: A. Scarabosio et al, 'Scaling of the divertor power spreading (S-factor) in open and closed divertor operation in JET and ASDEX Upgrade, Journal of Nuclear Materials, Vol. 463, 49-54 (2015)
