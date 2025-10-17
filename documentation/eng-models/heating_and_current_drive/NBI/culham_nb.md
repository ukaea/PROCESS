# Culham Neutral Beam Model | `culnbi()`

- `i_hcd_primary/i_hcd_secondary` = 8 



This routine calculates Neutral Beam current drive parameters
using the corrections outlined in AEA FUS 172 to the ITER method.
The result cannot be guaranteed for devices with aspect ratios far
from that of ITER (approx. 2.8).

| Output | Description |
|----------|-------------|
| $\mathtt{effnbss}$  | Neutral beam current drive efficiency in Amperes per Watt |
| $\mathtt{f_p_beam_injected_ions}$    | Fraction of NB power given to ions |
| $\mathtt{fshine}$   | Shine-through fraction of the beam |

$$
\mathtt{f_radius_beam_tangency_rmajor} = \frac{R_{\text{tan}}}{R_0}
$$

Where $R_{\text{tan}}$ is major radius at which the centre-line of the beam is tangential to the toroidal direction. This can be user defined

$$
\left(1+ \frac{1}{A}\right) < \mathtt{f_radius_beam_tangency_rmajor}
$$

A quick sanity check is done to make sure no negative roots are formed when calculating $\mathtt{dpath}$ this prevents setups where the NBI beam would miss the plasma


$$
\mathtt{dpath} = R_0 \sqrt{\left(1+\frac{1}{A}\right)^2-\mathtt{f_radius_beam_tangency_rmajor}^2}
$$

Beams topping cross section is calculated via $\mathtt{sigbeam}$ found [here](../NBI/nbi_overview.md/#beam-stopping-cross-section-sigbeam). This produces $\mathtt{sigstop}$

Calculate number of decay lengths to centre

$$
\mathtt{n_beam_decay_lengths_core} = \mathtt{dpath} \times n_{\text{e,0}} \times \mathtt{sigstop}
$$

Calculate the shine through fraction of the beam

$$
\mathtt{fshine} = e^{\left(-2 \times \mathtt{dpath} \times n_{\text{e,0}} \times \mathtt{sigstop}\right)}
$$

Deuterium and tritium beam densities

$$
\mathtt{dend} = n_{\text{ion}} \times (1-\mathtt{f_beam_tritium})
$$

$$
\mathtt{dent} = n_{\text{ion}} \times \mathtt{f_beam_tritium}
$$

Power split to the ions and electrons is calculated with the $\mathtt{cfnbi()}$ method found [here](../NBI/nbi_overview.md/#ion-coupled-power-cfnbi) and outputs $\mathtt{f_p_beam_injected_ions}$

## Current drive efficiency | `etanb2()`

This routine calculates the current drive efficiency in A/W of
a neutral beam system, based on the 1990 ITER model,
plus correction terms outlined in Culham Report AEA FUS 172.

| Input       | Description                          |
| :---------- | :----------------------------------- |
| $\mathtt{m_beam_amu}$      | beam ion mass (amu)   |
| $\mathtt{alphan}$, $\alpha_n$       | density profile factor   |
| $\mathtt{alphat}$, $\alpha_T$       | temperature profile factor  |
|  $\mathtt{aspect}$, $A$      |   aspect ratio                            |
|  $\mathtt{nd_plasma_electrons_vol_avg}$, $n_{\text{e}}$     |    volume averaged electron density $(\text{m}^{-3})$                           |
|  $\mathtt{nd_plasma_electron_line}$, $n_{\text{e,0}}$      |    line averaged electron density $(\text{m}^{-3})$                           |
|  $\mathtt{e_beam_kev}$      |  neutral beam energy $(\text{keV})$                             |
|  $\mathtt{f_radius_beam_tangency_rmajor}$      |   R_tangent / R_major for neutral beam injection                            |
|  $\mathtt{fshine}$      |  shine-through fraction of beam                             |
|  $\mathtt{rmajor}$, $R$      |  plasma major radius $(\text{m})$                              |
|  $\mathtt{rminor}$, $a$      |  plasma minor radius $(\text{m})$                             |
|  $\mathtt{temp_plasma_electron_density_weighted_kev}$      |    density weighted average electron temperature $(\text{keV})$                             |
|  $\mathtt{zeff}$, $Z_{\text{eff}}$     |   plasma effective charge                            |


Charge of beam ions
$$
\mathtt{zbeam} = 1.0
$$

Fitting factor (IPDG89)

$$
\mathtt{bbd} = 1.0
$$

Volume averaged electron density ($10^{20} \text{m}^{-3}$)

$$
\mathtt{dene20} = n_{\text{e,20}}
$$

Line averaged electron density ($10^{20} \text{m}^{-3}$)

$$    
\mathtt{dnla20} = n_{\text{(e,0) 20}} 
$$

Critical energy ($\text{MeV}$) (power to electrons = power to ions) (IPDG89)
N.B. temp_plasma_electron_density_weighted_kev is in keV

$$
\mathtt{ecrit} = 0.01 \times \mathtt{m_beam_amu} \times \mathtt{temp_plasma_electron_density_weighted_kev}
$$

Beam energy in MeV

$$
\mathtt{ebmev} = \frac{\mathtt{e_beam_kev}}{10^3}
$$

x and y coefficients of function J0(x,y) (IPDG89)

$$    
\mathtt{xjs} = \frac{\mathtt{ebmev}}{\mathtt{bbd}\times \mathtt{ecrit}}  
$$

$$
\mathtt{xj} = \sqrt{\mathtt{xjs}}
$$

$$
\mathtt{yj} = \frac{0.8 \times Z_{\text{eff}}}{\mathtt{m_beam_amu}}
$$

Fitting function J0(x,y)

$$
\mathtt{j0} = \frac{xjs}{(4.0 + 3.0 \times \mathtt{yj} + \mathtt{xjs} \times (\mathtt{xj} + 1.39 + 0.61 \times yj^{0.7}))}
$$

Effective inverse aspect ratio, with a limit on its maximum value

$$    
 \mathtt{epseff} = \text{min}(0.2, (0.5 / A))
$$   

Reduction in the reverse electron current
due to neoclassical effects

$$
\mathtt{gfac} = (1.55 + 0.85 / Z_{\text{eff}}) \times \sqrt{\mathtt{epseff}} - (0.2 + 1.55 / Z_{\text{eff}}) \times \mathtt{epseff}
$$

Reduction in the net beam driven current
due to the reverse electron current

$$
\mathtt{ffac} = 1.0 - \frac{\mathtt{zbeam}}{Z_{\text{eff}}} \times (1.0 - \mathtt{gfac})
$$

Normalisation to allow results to be valid for
non-ITER plasma size and density:

Line averaged electron density ($10^{20} \text{m}^{-3}$) normalised to ITER

$$    
\mathtt{nnorm} = 1.0
$$

Distance along beam to plasma centre

$$
\mathtt{r} = \text{max}(R, R \times \mathtt{f_radius_beam_tangency_rmajor})
$$

$$
\mathtt{eps1} = a / \mathtt{r}
$$


$$
\mathtt{d} = R \times \sqrt{((1.0 + \mathtt{eps1})^2 - \mathtt{f_radius_beam_tangency_rmajor}^2)}
$$

Distance along beam to plasma centre for ITER
assuming a tangency radius equal to the major radius
    
$$
\mathtt{epsitr} = 2.15 / 6.0
$$  
  
$$  
\mathtt{dnorm} = 6.0 \times \sqrt{(2.0 \times \mathtt{epsitr} + \mathtt{epsitr}^2)}
$$

Normalisation to beam energy (assumes a simplified formula for
the beam stopping cross-section)

$$
    \mathtt{ebnorm} = \mathtt{ebmev} \times ((\mathtt{nnorm} \times \mathtt{dnorm}) / (n_{\text{e,0}} \times \mathtt{d})) ^{1.0 / 0.78)}
$$

A_bd fitting coefficient, after normalisation with ebnorm

$$   
\mathtt{abd} = (
    0.107
    \times (1.0 - 0.35 \times \alpha_n + 0.14 \times \alpha_n^2)
    \times (1.0 - 0.21 \times \alpha_T)
    \times (1.0 - 0.2 \times \mathtt{ebnorm} + 0.09 \times \mathtt{ebnorm}^2)
    )
$$

Normalised current drive efficiency ($\text{A/W} \text{m}^{2}$) (IPDG89)

$$
\mathtt{gamnb} = 5.0 \times \mathtt{abd} \times 0.1 \times \mathtt{temp_plasma_electron_density_weighted_kev} \times (1.0 - \mathtt{fshine}) \times \mathtt{f_radius_beam_tangency_rmajor} \times \frac{\mathtt{j0}}{0.2} \times \mathtt{ffac}
$$

Current drive efficiency (A/W)

$$
\text{Current drive efficiency [A/W]} = \frac{\mathtt{gamnb}}{\mathtt{dene20}\times R}
$$
