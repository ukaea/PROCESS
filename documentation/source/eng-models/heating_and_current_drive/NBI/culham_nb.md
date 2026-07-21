# Culham Neutral Beam Model | `culnbi()`

- `i_hcd_primary/i_hcd_secondary` = 8 



This routine calculates Neutral Beam current drive parameters
using the corrections outlined in AEA FUS 172 to the ITER method.
The result cannot be guaranteed for devices with aspect ratios far
from that of ITER (approx. 2.8).

| Output | Description |
|----------|-------------|
| $\texttt{effnbss}$  | Neutral beam current drive efficiency in Amperes per Watt |
| $\texttt{f_p_beam_injected_ions}$    | Fraction of NB power given to ions |
| $\texttt{fshine}$   | Shine-through fraction of the beam |

$$
\texttt{f_radius_beam_tangency_rmajor} = \frac{R_{\text{tan}}}{R_0}
$$

Where $R_{\text{tan}}$ is major radius at which the centre-line of the beam is tangential to the toroidal direction. This can be user defined

$$
\left(1+ \frac{1}{A}\right) < \texttt{f_radius_beam_tangency_rmajor}
$$

A quick sanity check is done to make sure no negative roots are formed when calculating $\texttt{dpath}$ this prevents setups where the NBI beam would miss the plasma


$$
\texttt{dpath} = R_0 \sqrt{\left(1+\frac{1}{A}\right)^2-\texttt{f_radius_beam_tangency_rmajor}^2}
$$

Beams topping cross section is calculated via $\texttt{sigbeam}$ found [here](../NBI/nbi_overview.md/#beam-stopping-cross-section-sigbeam). This produces $\texttt{sigstop}$

Calculate number of decay lengths to centre

$$
\texttt{n_beam_decay_lengths_core} = \texttt{dpath} \times n_{\text{e,0}} \times \texttt{sigstop}
$$

Calculate the shine through fraction of the beam

$$
\texttt{fshine} = e^{\left(-2 \times \texttt{dpath} \times n_{\text{e,0}} \times \texttt{sigstop}\right)}
$$

Deuterium and tritium beam densities

$$
\texttt{dend} = n_{\text{ion}} \times (1-\texttt{f_beam_tritium})
$$

$$
\texttt{dent} = n_{\text{ion}} \times \texttt{f_beam_tritium}
$$

Power split to the ions and electrons is calculated with the $\texttt{cfnbi()}$ method found [here](../NBI/nbi_overview.md/#ion-coupled-power-cfnbi) and outputs $\texttt{f_p_beam_injected_ions}$

## Current drive efficiency | `etanb2()`

This routine calculates the current drive efficiency in A/W of
a neutral beam system, based on the 1990 ITER model,
plus correction terms outlined in Culham Report AEA FUS 172.

| Input       | Description                          |
| :---------- | :----------------------------------- |
| $\texttt{m_beam_amu}$      | beam ion mass (amu)   |
| $\texttt{alphan}$, $\alpha_n$       | density profile factor   |
| $\texttt{alphat}$, $\alpha_T$       | temperature profile factor  |
|  $\texttt{aspect}$, $A$      |   aspect ratio                            |
|  $\texttt{nd_plasma_electrons_vol_avg}$, $n_{\text{e}}$     |    volume averaged electron density $(\text{m}^{-3})$                           |
|  $\texttt{nd_plasma_electron_line}$, $n_{\text{e,0}}$      |    line averaged electron density $(\text{m}^{-3})$                           |
|  $\texttt{e_beam_kev}$      |  neutral beam energy $(\text{keV})$                             |
|  $\texttt{f_radius_beam_tangency_rmajor}$      |   R_tangent / R_major for neutral beam injection                            |
|  $\texttt{fshine}$      |  shine-through fraction of beam                             |
|  $\texttt{rmajor}$, $R$      |  plasma major radius $(\text{m})$                              |
|  $\texttt{rminor}$, $a$      |  plasma minor radius $(\text{m})$                             |
|  $\texttt{temp_plasma_electron_density_weighted_kev}$      |    density weighted average electron temperature $(\text{keV})$                             |
|  $\texttt{zeff}$, $Z_{\text{eff}}$     |   plasma effective charge                            |


Charge of beam ions
$$
\texttt{zbeam} = 1.0
$$

Fitting factor (IPDG89)

$$
\texttt{bbd} = 1.0
$$

Volume averaged electron density ($10^{20} \text{m}^{-3}$)

$$
\texttt{dene20} = n_{\text{e,20}}
$$

Line averaged electron density ($10^{20} \text{m}^{-3}$)

$$    
\texttt{dnla20} = n_{\text{(e,0) 20}} 
$$

Critical energy ($\text{MeV}$) (power to electrons = power to ions) (IPDG89)
N.B. temp_plasma_electron_density_weighted_kev is in keV

$$
\texttt{ecrit} = 0.01 \times \texttt{m_beam_amu} \times \texttt{temp_plasma_electron_density_weighted_kev}
$$

Beam energy in MeV

$$
\texttt{ebmev} = \frac{\texttt{e_beam_kev}}{10^3}
$$

x and y coefficients of function J0(x,y) (IPDG89)

$$    
\texttt{xjs} = \frac{\texttt{ebmev}}{\texttt{bbd}\times \texttt{ecrit}}  
$$

$$
\texttt{xj} = \sqrt{\texttt{xjs}}
$$

$$
\texttt{yj} = \frac{0.8 \times Z_{\text{eff}}}{\texttt{m_beam_amu}}
$$

Fitting function J0(x,y)

$$
\texttt{j0} = \frac{xjs}{(4.0 + 3.0 \times \texttt{yj} + \texttt{xjs} \times (\texttt{xj} + 1.39 + 0.61 \times yj^{0.7}))}
$$

Effective inverse aspect ratio, with a limit on its maximum value

$$    
 \texttt{epseff} = \text{min}(0.2, (0.5 / A))
$$   

Reduction in the reverse electron current
due to neoclassical effects

$$
\texttt{gfac} = (1.55 + 0.85 / Z_{\text{eff}}) \times \sqrt{\texttt{epseff}} - (0.2 + 1.55 / Z_{\text{eff}}) \times \texttt{epseff}
$$

Reduction in the net beam driven current
due to the reverse electron current

$$
\texttt{ffac} = 1.0 - \frac{\texttt{zbeam}}{Z_{\text{eff}}} \times (1.0 - \texttt{gfac})
$$

Normalisation to allow results to be valid for
non-ITER plasma size and density:

Line averaged electron density ($10^{20} \text{m}^{-3}$) normalised to ITER

$$    
\texttt{nnorm} = 1.0
$$

Distance along beam to plasma centre

$$
\texttt{r} = \text{max}(R, R \times \texttt{f_radius_beam_tangency_rmajor})
$$

$$
\texttt{eps1} = a / \texttt{r}
$$


$$
\texttt{d} = R \times \sqrt{((1.0 + \texttt{eps1})^2 - \texttt{f_radius_beam_tangency_rmajor}^2)}
$$

Distance along beam to plasma centre for ITER
assuming a tangency radius equal to the major radius
    
$$
\texttt{epsitr} = 2.15 / 6.0
$$  
  
$$  
\texttt{dnorm} = 6.0 \times \sqrt{(2.0 \times \texttt{epsitr} + \texttt{epsitr}^2)}
$$

Normalisation to beam energy (assumes a simplified formula for
the beam stopping cross-section)

$$
    \texttt{ebnorm} = \texttt{ebmev} \times ((\texttt{nnorm} \times \texttt{dnorm}) / (n_{\text{e,0}} \times \texttt{d})) ^{1.0 / 0.78)}
$$

A_bd fitting coefficient, after normalisation with ebnorm

$$   
\texttt{abd} = (
    0.107
    \times (1.0 - 0.35 \times \alpha_n + 0.14 \times \alpha_n^2)
    \times (1.0 - 0.21 \times \alpha_T)
    \times (1.0 - 0.2 \times \texttt{ebnorm} + 0.09 \times \texttt{ebnorm}^2)
    )
$$

Normalised current drive efficiency ($\text{A/W} \text{m}^{2}$) (IPDG89)

$$
\texttt{gamnb} = 5.0 \times \texttt{abd} \times 0.1 \times \texttt{temp_plasma_electron_density_weighted_kev} \times (1.0 - \texttt{fshine}) \times \texttt{f_radius_beam_tangency_rmajor} \times \frac{\texttt{j0}}{0.2} \times \texttt{ffac}
$$

Current drive efficiency (A/W)

$$
\text{Current drive efficiency [A/W]} = \frac{\texttt{gamnb}}{\texttt{dene20}\times R}
$$
