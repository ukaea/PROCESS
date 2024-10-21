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

### Beam slowing down properties | `beamfus()`

------------------------