# ECRH with Cutoff

- `iefrf/iefrffix` = 13

| Input | Description |
|-------|-------------|
| `dene`, $n_{\text{e}}$ | Avergae electron temperature $\left[10^{19}\text{m}^{-3}\right]$ |
| `te`, $T_{\text{e}}$ | Avergae electron temperature $\left[\text{keV}\right]$ |
| `rmajor`, $R_0$ | Major radius $\left[\text{m}\right]$ |
| `bt`, $B_{\text{T}}$ | Toroidal magnetic field $\left[\text{T}\right]$ |
| `zeff`, $Z_{\text{eff}}$ | Effective charge |
| `harnum` | Harmonic number |
| `mode` | RF mode |

----



$$
\mathtt{fc} = \frac{1}{2\pi}eB_{\text{T}}m_{\text{e}}
$$

$$
\mathtt{fp} = \frac{1}{2\pi}\sqrt{\frac{n_{\text{e,19}}e^2}{m_{\text{e}}\epsilon_0}}
$$

Apply effective charge correction from GRAY study

$$
\mathtt{xi_{CD}} = 0.18\left(\frac{4.8}{2+Z_{\text{eff}}}\right)
$$

$$
\mathtt{effrfss} = \frac{\mathtt{xi_{CD}}T_{\text{e}}}{3.27R_0n_{\text{e,19}}}
$$
                
For the O-mode case:

$$
\mathtt{f_{cutoff}} = \mathtt{fp}
$$

For the X-mode case:

$$
\mathtt{f_{cutoff}} = 0.5\left(\mathtt{fc}+\sqrt{\mathtt{harnum}\times\mathtt{fc}^2+4\mathtt{fp}^2}\right)
$$
               
Plasma coupling only occurs if the plasma cut-off is below the cyclotron harmonic
(a = 0.1).  This controls how sharply the transition is reached
                
$$
\mathtt{cutoff_{factor}} = 0.5\left(1+\tanh\left({\left(\frac{2}{a}\right)((\mathtt{harnum}\times \mathtt{fc} -\mathtt{f_cutoff})/\mathtt{fp -a })}\right)\right)
$$

$$
\text{Current drive efficiency [A/W]} = \mathtt{effrfss} \times \mathtt{cutoff_{factor}}
$$

<figure markdown>
![ECRH Cutoff](../images/ecrh_cutoff.png){ width = "100"}
<figcaption>Figure 1: The variation in current drive efficiency as a function of toroidal magnetic field at different harmonics and modes </figcaption>
</figure>