# Resistive Plasma Heating

The ohmic component of the plasma heating is given by that from the ITER 1989 Physics Design Guidelines[^1]

-------------------

## Plasma Ohmic Heating Power | `plasma_ohmic_heating()`

Using the resistive loop voltage for a reference profile of parabolic shape with:

$$
\alpha_n \approx 0.5, \alpha_T \approx 1.0,  \alpha_J \approx 1.5
$$

We calculate the plasma resistance as:

$$
\mathtt{res\_plasma} = \Omega_{\text{plasma}} \approx 2.15 \times 10^{-3} Z_{\text{eff}}\langle \gamma_{\text{NC}} \rangle \frac{R_0}{\kappa a^2} \frac{1}{T_{10}^{1.5}}
$$

where $Z_{\text{eff}}$ is the plasma effective charge and $T_{10}$ is the density-weighted temperature in units of 10 keV.

The neoclassical (average) resistivity enhancement factor $\left(\langle \gamma_{\text{NC}} \rangle \right)$ is given by an empirical fit:

$$
\langle \gamma_{\text{NC}} \rangle = 4.3 -0.6A
$$

where $A$ is valid in the range of 2.5 - 4.0. If $A < 2.5$ then $\langle \gamma_{\text{NC}}\rangle$ is et equal to 1.0

---------------

The ohmic heating power in MW is then simply found using Joules law:

$$
\mathtt{p\_plasma\_ohmic\_mw} = P_{\text{OH}} = 1\times 10^{-6} \left[f_{\text{ind}}I_{\text{p}}^2\Omega_{\text{plasma}}\right]
$$

where $f_{\text{ind}}$ is the fraction of plasma current driven by inductive means.

Likewise, the ohmic heating per unit volume simply as:

$$
\mathtt{pden\_plasma\_ohmic\_mw} = \frac{\mathtt{p\_plasma\_ohmic\_mw}}{V_{\text{p}}}
$$


[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',