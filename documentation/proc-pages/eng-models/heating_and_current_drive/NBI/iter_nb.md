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


## Current drive efficiency
Uses the `etanb` method.

effnbss = current_drive_variables.frbeam * self.etanb(
    physics_variables.abeam,
    physics_variables.alphan,
    physics_variables.alphat,
    physics_variables.aspect,
    physics_variables.dene,
    current_drive_variables.enbeam,
    physics_variables.rmajor,
    physics_variables.ten,
    physics_variables.zeff,
)

zbeam = 1.0
        bbd = 1.0

        dene20 = 1e-20 * dene

Ratio of E_beam/E_crit

$$
\mathtt{xjs} = \frac{ebeam}{10 m_{\text{u,ion}}T_e} 
$$

$$
\mathtt{xj} = \sqrt{\mathtt{xjs}}
$$
    
$$
\mathtt{yj} = 0.8 \frac{Z_{\text{eff}}}{m_{\text{u,beam}}}
$$
        
$$
rjfunc = \frac{xjs}{((4.0 + 3.0\mathtt{yj} + \mathtt{xjs}(\mathtt{xj} + 1.39 + 0.61\mathtt{yj}^{0.7})))}
$$

$$
\mathtt{epseff} = \frac{0.5}{A}
$$

$$
\mathtt{gfac} = \left(1.55 + \frac{0.85}{Z_{\text{eff}}}\right)\left(\sqrt{\mathtt{epseff}}-(0.2+\frac{1.55}{Z_{\text{eff}}})\mathtt{epseff}\right)
$$

$$
\mathtt{ffac} = \frac{1}{\mathtt{zbeam}} - \frac{(1-\mathtt{gfac})}{Z_{\text{eff}}}
$$
    
$$
\mathtt{abd} = 0.107 (1-0.35\mathtt{alphan}+0.14\mathtt{alphan}^2)(1-0.21\mathtt{alphat})(1-0.2E-3\mathtt{ebeam}+0.09E-6\mathtt{ebeam}^2)
$$

$$
\text{Current drive efficiency [A/W]} = \mathtt{abd}(\frac{5}{R_0})(0.1\frac{\mathtt{ten}}{n_{e,20}})\frac{\mathtt{rjfunc}}{0.2}\mathtt{ffac}
$$

        

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)