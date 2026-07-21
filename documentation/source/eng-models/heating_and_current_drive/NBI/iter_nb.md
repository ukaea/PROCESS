# ITER Neutral Beam Model | `iternb()`

- `i_hcd_primary/i_hcd_secondary` = 5

| Output | Description |
|----------|-------------|
| $\texttt{effnbss}$  | Neutral beam current drive efficiency in $\text{A/W}$ |
| $\texttt{f_p_beam_injected_ions}$    | Fraction of NB power given to ions |
| $\texttt{fshine}$   | Shine-through fraction of the beam |

This model calculates the current drive parameters for a neutral beam system, based on the 1990 ITER model.[^1]


Firstly the beam access is checked for such that
$$
\bigg(1+ \frac{1}{A}\bigg) > (R_{\text{tangential}}/R_0)
$$

The beam path length to centre is calculated:

$$
\underbrace{\texttt{dpath}}_{\text{Path length to centre}} = R_0 \sqrt{\left(\left(1+\frac{1}{A}\right)^2-\texttt{f_radius_beam_tangency_rmajor}^2\right)}
$$


Beam stopping cross-section ($\sigma_{\text{beam}}$) is calculated using the `sigbeam` method described [here](nbi_overview.md) :


Calculate number of electron decay lengths to centre

$$
\tau_{\text{beam}} = \texttt{dpath}\times n_e \sigma_{\text{beam}}
$$

Shine-through fraction of beam:
$$
f_{\text{shine}} = e^{(-2.0 \times  \texttt{dpath} \times  n_e  \sigma_{\text{beam}})} \\
$$

Deuterium and tritium beam densities:
$$
n_D = n_i * (1.0 - \texttt{f_beam_tritium}) 
$$

$$
n_T = n_i * \texttt{f_beam_tritium}
$$

Power split to ions / electrons is calculated via the the `cfnbi` method described [here](nbi_overview.md)


## Current drive efficiency | `etanb()`

This routine calculates the current drive efficiency of
a neutral beam system, based on the 1990 ITER model.
AEA FUS 251: A User's Guide to the PROCESS Systems Code
ITER Physics Design Guidelines: 1989 IPDG89, N. A. Uckan et al,
ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990

| Input   | Description                                               |
|---------|-----------------------------------------------------------|
| $\texttt{m_beam_amu}$, $m_{\text{u,ion}}$   | Beam ion mass ($\text{amu}$)                                       |
| $\texttt{alphan}$  | Density profile factor                                    |
| $\texttt{alphat}$  | Temperature profile factor                                |
| $\texttt{aspect}$, $A$  | Aspect ratio                                              |
| $\texttt{dene20}$, $n_{\text{e,20}}$    | Volume averaged electron density ($10^{20} \text{m}^{-3}$)                  |
| $\texttt{ebeam}$   | Neutral beam energy ($\text{keV}$)                                 |
| $\texttt{rmajor}$, R  | Plasma major radius ($\text{m}$)                                   |
| $\texttt{temp_plasma_electron_density_weighted_kev}$     | Density weighted average electron temperature ($\text{keV}$)       |
| $\texttt{zeff}$, $Z_{\text{eff}}$    | Plasma effective charge                                   |

| Output  | Description                                               |
|---------|-----------------------------------------------------------|
| $\texttt{etanb}$   | Neutral beam current drive efficiency in $\text{A/W}$ |




$$
\texttt{zbeam} = 1.0
$$

$$
\texttt{bbd} = 1.0
$$


Ratio of E_beam/E_crit

$$
\texttt{xjs} = \frac{\texttt{ebeam}}{10 \  m_{\text{u,ion}} \ T_e} 
$$

$$
\texttt{xj} = \sqrt{\texttt{xjs}}
$$
    
$$
\texttt{yj} = 0.8 \frac{Z_{\text{eff}}}{m_{\text{u,ion}}}
$$
        
$$
\texttt{rjfunc} = \frac{\texttt{xjs}}{((4.0 + 3.0\texttt{yj} + \texttt{xjs}(\texttt{xj} + 1.39 + 0.61\texttt{yj}^{0.7})))}
$$

$$
\texttt{epseff} = \frac{0.5}{A}
$$

$$
\texttt{gfac} = \left(1.55 + \frac{0.85}{Z_{\text{eff}}}\right)\left(\sqrt{\texttt{epseff}}-\left(0.2+\frac{1.55}{Z_{\text{eff}}}\right)\texttt{epseff}\right)
$$

$$
\texttt{ffac} = \frac{1}{\texttt{zbeam}} - \frac{(1-\texttt{gfac})}{Z_{\text{eff}}}
$$
    
$$
\texttt{abd} = 0.107 (1-0.35 \ \texttt{alphan}+0.14 \ \texttt{alphan}^2)(1-0.21 \ \texttt{alphat})(1-0.2\times 10^{-3}\texttt{ebeam}+0.09\times 10^{-6} \ \texttt{ebeam}^2)
$$

$$
\text{Current drive efficiency [A/W]} = \texttt{abd} \times\frac{5}{R_0} \times0.1\frac{\texttt{temp_plasma_electron_density_weighted_kev}}{n_{\text{e},20}} \times \frac{\texttt{rjfunc}}{0.2}\texttt{ffac}
$$

        

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)