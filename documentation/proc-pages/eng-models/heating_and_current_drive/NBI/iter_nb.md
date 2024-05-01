# ITER Neutral Beam Model

This model calculates the current drive parameters for a neutral beam system, based on the 1990 ITER model.[^1]


Firstly the beam access is checked for such that
$$
\bigg(1+ \frac{1}{A}\bigg) > (R_{\text{tangential}}/R_0)
$$

The beam path length to centre is calculated:

$$
\underbrace{\mathtt{dpath}}_{\text{Path length to centre}} = R_0 \sqrt{\left(\left(1+\frac{1}{A}\right)^2-\mathtt{frbeam}^2\right)}
$$


Beam stopping cross-section ($\sigma_{\text{beam}}$) is calculated using the `sigbeam` method described [here](nbi_overview.md) :


Calculate number of electron decay lengths to centre

$$
\tau_{\text{beam}} = \mathtt{dpath}\times n_e \sigma_{\text{beam}}
$$

Shine-through fraction of beam:
$$
f_{\text{shine}} = e^{(-2.0 \times  \mathtt{dpath} \times  n_e  \sigma_{\text{beam}})} \\
$$

Deuterium and tritium beam densities:
$$
n_D = n_i * (1.0 - \mathtt{ftritbm}) 
$$

$$
n_T = n_i * \mathtt{ftritbm}
$$

Power split to ions / electrons is calculated via the the `cfnbi` method described [here](nbi_overview.md)


## Current drive efficiency (`etanb`)
Uses the `etanb` method.

$$
\mathtt{zbeam} = 1.0
$$

$$
\mathtt{bbd} = 1.0
$$


Ratio of E_beam/E_crit

$$
\mathtt{xjs} = \frac{\mathtt{ebeam}}{10 \  m_{\text{u,ion}} \ T_e} 
$$

$$
\mathtt{xj} = \sqrt{\mathtt{xjs}}
$$
    
$$
\mathtt{yj} = 0.8 \frac{Z_{\text{eff}}}{m_{\text{u,beam}}}
$$
        
$$
\mathtt{rjfunc} = \frac{\mathtt{xjs}}{((4.0 + 3.0\mathtt{yj} + \mathtt{xjs}(\mathtt{xj} + 1.39 + 0.61\mathtt{yj}^{0.7})))}
$$

$$
\mathtt{epseff} = \frac{0.5}{A}
$$

$$
\mathtt{gfac} = \left(1.55 + \frac{0.85}{Z_{\text{eff}}}\right)\left(\sqrt{\mathtt{epseff}}-\left(0.2+\frac{1.55}{Z_{\text{eff}}}\right)\mathtt{epseff}\right)
$$

$$
\mathtt{ffac} = \frac{1}{\mathtt{zbeam}} - \frac{(1-\mathtt{gfac})}{Z_{\text{eff}}}
$$
    
$$
\mathtt{abd} = 0.107 (1-0.35 \ \mathtt{alphan}+0.14 \ \mathtt{alphan}^2)(1-0.21 \ \mathtt{alphat})(1-0.2\times 10^{-3}\mathtt{ebeam}+0.09\times 10^{-6} \ \mathtt{ebeam}^2)
$$

$$
\text{Current drive efficiency [A/W]} = \mathtt{abd} \times\frac{5}{R_0} \times0.1\frac{\mathtt{ten}}{n_{e,20}} \times \frac{\mathtt{rjfunc}}{0.2}\mathtt{ffac}
$$

        

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)