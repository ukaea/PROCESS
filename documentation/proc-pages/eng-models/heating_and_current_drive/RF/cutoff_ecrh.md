# Cuttoff

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
