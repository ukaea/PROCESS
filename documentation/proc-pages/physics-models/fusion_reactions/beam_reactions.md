## NBI Fusion

### Beam fusion cross section | `_sigbmfus()`

### Beam fusion reaction rate | `sgvhot()`

### Beam fusion reaction rate integrand | `_hot_beam_fusion_reaction_rate_integrand()`

### Beam fusion alpha power | `alpha_power_beam()`

### Beam energy given to ions | `beam_energy_to_ions()`

#### Derivation

The beam particle is born at energy $\mathrm{E}_0$, which is much greater than either $T_{\text{i}}$ or $T_{\text{e}}$, and gives up its energy to the ions and electrons in the process of slowing down. The power transferred to the ions at time $t$ is:

$$
P=A Z^2 \sqrt{\frac{M}{E(t)}} .
$$

The total energy transferred to the ions is

$$
W_1=\int A Z^2 \sqrt{\frac{M}{E(t)}} \  \mathrm{dt} \\
=-\int_{E_f}^{E_0} A Z^2 \sqrt{\frac{\bar{M}}{E}}\left(\frac{\mathrm{dt}}{\mathrm{dE}}\right) \mathrm{dE},
$$

where $E_{\mathrm{f}}$ is the final energy. Since $E_{\mathrm{f}} \ll \mathrm{E}_0$, its value is not important and we can take $E_{\mathrm{f}} = 0$ without introducing a significant error. The fraction $f_{\mathrm{i}}$ of the initial energy given to the ions can be written as

$$
f_1=\frac{W_i}{E_0}=\frac{2 v_c^2}{v_0^2} \int_0^{X_c}  \frac{x}{1+x^3} \  \mathrm{dx}
$$

where

$$
E_{\text{c}}=\frac{M v_{\text{c}}^2}{2},  \quad E_0=\frac{M v_0^2}{2}, \quad X_c= \frac{v_0}{v_c}
$$

The integral can be evaluated analytically to yield:

$$
f_i=\frac{2}{X_c^2}\left\{\frac{1}{\sqrt{3}} \tan ^{-1}\left(\frac{2 X_c-1}{\sqrt{3}}\right)+\frac{1}{\sqrt{3}} \tan ^{-1}\left(\frac{1}{\sqrt{3}}\right)-\frac{1}{6} \ln \left[\frac{X_c^2+2 X_c+1}{X_c^2-X_c^2+1}\right]\right\} .
$$

The fraction $f_{\text{e}}$ of the initial energy transferred to the electrons is simply $1-f_1$.

------------------

### Neutral beam alpha power and ion energy | `beamcalc()`

Set the fraction of the beam current for deuterium and tritium

calculate the velocity at which the deuterium ions are going at the critical speed


### Beam slowing down properties | `beamfus()`

#### Derivation of beam slowing down rate

The rate of slowing down of a test particle of mass $M$, charge $Z\text{e}$ and energy $E$, due to Coulomb collisions with a background species off mass $m_{\text{j}}$, charge $Z_{\text{j}}\text{e}$, density $n_{\text{j}}$ and temperature $T_{\text{j}}$ is given by:

$$
\frac{\mathrm{dE}}{\mathrm{dt}}=\left[-\Phi\left(x_{\text{j}}\right)+x_{\text{j}}\left(1+\frac{m_{\text{j}}}{M} \Phi^{\prime}\left(x_{\text{j}}\right)\right)\right] \frac{4 \pi n_{\text{j}}}{m_{\text{j}} V}\left(\frac{Z z_{\text{j}} Z \text{e}^2}{4 \pi \varepsilon_0}\right)^2 \ln \Lambda_{\text{j}},
$$

where $\Phi(x)$ is the error function,

$$
\Phi^{\prime}(x) = \frac{\mathrm{d\Phi}}{\mathrm{dx}}, \ \ V = \sqrt{\frac{2E}{M}}, \ \ V_{\text{j}} = \sqrt{\frac{2kT_{\text{j}}}{m_{\text{j}}}}, \ \ x_{\text{j}} = \frac{V}{V_{\text{j}}}
$$

and $\ln \Lambda_{\text{j}}$ is the usual Coulomb logarithm for the test particle and background species interactions.

For the contribution to the slowing down due to interaction with the background ions, we can use a large argument expansion for the error function. This is because the fusion born ions have a velocity $V$, much greater than the thermal velocity $V_{\text{j}}$, of the background ions. The velocity of the fast ions is much less than the thermal velocity of electrons, however. For the electron contribution to the slowing down we use the small argument expansion for the error function. The net slowing down rate is then given by
$$
\frac{\mathrm{d E}}{\mathrm{d t}}=-\frac{A Z^2 \sqrt{M}}{\sqrt{E}}-\frac{B Z^2 E}{M}
$$
where the coefficients $A$ and $B$ are given by
$$
\begin{aligned}
& A=\frac{4 \pi}{\sqrt{2}}\left(\frac{\text{e}^2}{4 \pi \varepsilon_0}\right)^2 \sum_j\left(\frac{n_{\text{j}} Z_{\text{j}}^2}{m_{\text{j}}} \ln \Lambda_j\right) \\
& B=\frac{16 \sqrt{\pi}}{3 k T_{\text{e}}} \sqrt{\frac{m_{\text{e}}}{2 k T_{\text{e}}}}\left(\frac{\text{e}^2}{4 \pi \varepsilon_0}\right)^2 n_{\text{e}} \ln \Lambda_{\text{e}}
\end{aligned}
$$

The sum over $\text{j}$ in $A$ is over the various ionic species. The quantities with subscript $\text{e}$ refer to the electrons.

$\frac{\mathrm{d E}}{\mathrm{d t}}$ can be rewritten in the form,

$$
\frac{d E}{d t}=-\frac{2 E}{\tau_s}\left[1+\left(\frac{E}{E}\right)^{3 / 2}\right]
$$

The critical energy $E_{\text{c}}$ is given by

$$
E_c=\left[\frac{3 \sqrt{\pi}}{4} \frac{M^{3 / 2}}{n_{\text{e}} \sqrt{m}_{\text{e}}} \sum_j\left(\frac{n_{\text{j}} z_{\text{j}}^2}{m_{\text{j}}} \ln \Lambda_{\text{j}}\right) \frac{1}{\ln \Lambda_{\text{e}}}\right]^{2 / 3} k T_{\text{e}}
$$

Some authors take $\ln \Lambda_i=\ln \Lambda_{\text{e}}$ in the expression for $E_C$, but this is not a good approximation for fast ion slowing down in fusion plasmas since $\ln \Lambda_{\text{e}} \approx 17$, while $\ln \Lambda_i \approx 22$. The slowing down time, $\tau_s$ is given by

$$
\tau_s=\frac{3(k T_{\text{e}})^{3 / 2}}{4 \sqrt{2 \pi m_{\text{e}}} Z^2}\left(\frac{4 \pi \varepsilon_0}{\text{e}^2}\right)^2 \frac{M}{n_{\text{e}} \ln \Lambda_{\text{e}}}
$$

When the particle energy $E$ is above $E_C$ the contribution of the electrons to the slowing down is larger than that of the ions. The slowing down time $\tau_s$ is actually the time scale for $V$ to decrease due to electron drag, i.e. $\tau_{s}=-V /(d V / d t) e^*$

$\blacksquare$

-------------------

Calculate the beam ion slowing down time given by:

$$
\tau_{\text{beam-slow}} = 1.99 \times 10^{19}\left(A_{\text{D}}\left(1.0-f_{\text{tritium-beam}}\right)+(A_{\text{T}}f_{\text{tritium-beam}})\right)\frac{\langle T_{\text{e}}\rangle^{3/2}}{\langle n_{\text{e}} \rangle \Lambda_{\text{ie}}}
$$

The alpha particles are born with an energy of 3.5 MeV and initially slow down mainly by collisions with electrons. At a critical energy $E_{\text{crit}}$ the rate of loss to the ions becomes equal to that to the electrons, and at lower
energies the loss to the ions predominates.[^1]

$$
E_{\text{crit}} = 14.8 A T_{\text{e}}\left[\frac{1}{n_{\text{e} \ln{ \Lambda_{\text{e}}}}}\left[\Sigma \frac{n_j Z_j^2\ln{\Lambda_{\text{i}}}}{A_j}\right]\right]^{2/3} \ [\text{eV}]
$$

This can be approximated to:

$$
E_{\text{crit}} \approx 0.1AT_{10} \ \ [\text{MeV}]
$$

Set the plasma tritium and ion densities.

Run `beamcalc()`

------------------------


[^1]: J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3, pp. 1015–1103,Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.