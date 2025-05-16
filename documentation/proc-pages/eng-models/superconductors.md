

`PROCESS` offers a variety of different superconducting coil types, both low-temperature superconductors (LTS) and high-temperature superconductors (HTS).


---------------

## Low-temperature superconductors

### Nb<sub>3</sub>Sn

#### Western Superconducting technologies Nb<sub>3</sub>Sn

The following parameters [^2] following the [Bottura scaling](#bottura-scaling--bottura_scaling) are shown below:

- **\( C_{\text{a1}} \)**: 50.06
- **\( C_{\text{a2}} \)**: 0
- **\( \epsilon_{\text{0,a}} \)**: 0.312%
- **\( \epsilon_{\text{m}} \)**: -0.059%
- **\( B_{C20 \text{max}}^* \)**: 33.24 $\text{T}$
- **\( T_{C0 \text{max}}^* \)**: 16.34 $\text{K}$
- **\( C \)**: 83075 $\text{AT/mm}^2$
- **\( p \)**: 0.593
- **\( q \)**: 2.156

The critical surface at zero strain looks as follows: 

<figure markdown>
![Wester Superconducting Ltd. Nb3Sn](./images/western_superconducting_Nb3Sn_zero_strain.png){ width = "100"}
<figcaption>Figure 1: Critical current density surface for the Wester Superconducting Ltd. Nb<sub>3</sub>Sn superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/western_superconducting_Nb3Sn_zero_strain.html)** :bar_chart:

!!! quote "Excerpt"
    *"TF Coil: The initial choice of EUTF4 right leg (OST) and its parametrization can be replaced by the actual WST strand (ø = 1.5 mm) used and tested in the DEMO TF prototypes manufactured in 2014 and measured in 2015 and 2016. The EUTF4 OST strand parameterisation can still be employed in the conductor design, if the designers prefer. The parameters are measured in the scope of the 2014 activities and reported below. Compared to the EUTF4 right leg (OST), the performance of the WST at low strain is superior by about 10%."*[^2]

----------------

##### Bottura scaling | `bottura_scaling()`

This is a generic scaling proposed for the characterization and production of 
ITER Nb<sub>3</sub>Sn strands. This is also known as the "ITER-2008 parametrization."[^1]

The datasets span 10 years of R&D and production, with 100 to 500 data points each, covering strain (1.5% to 0.4%), temperature (2.35 to 16 K), and field (0.5 to 19 T). The ITER-2008 parameterization achieves an average accuracy error of 3.8 Amps, with the best at 1.5 Amps and the worst at 7.5 Amps. Compared to Durham University's full parameterization, ITER-2008 yields root mean squared errors 1.5 times larger, which is significant but not dramatic.

The strain function is suitable only in the moderate strain region, down to 0.8%. Beyond this, the measured behavior shows an inflection and curvature change not captured in the equations. Fitting strain data beyond an intrinsic compressive strain of 0.5% impacts the moderate strain regime and reduces fit accuracy.


$$
J_{\text{C}} = \frac{C}{B}s(\epsilon)(1-t^{1.52})(1-t^2)b^p (1-b)^q
$$

$$
B_{C2}^*(T,\epsilon) = B_{C20 \text{max}}^* s(\epsilon)(1-t^{1.52})
$$

$$
T_{C}^*(T,\epsilon) = T_{C0 \text{max}}^* \left[s(\epsilon)\right]^{\frac{1}{3}}\left(1-\frac{B}{B_{C2}^*(0,\epsilon)}\right)^{\frac{1}{1.52}}
$$

$$
s(\epsilon) = 1 + \frac{C_{a1}\left(\sqrt{\epsilon_{sh}^2+\epsilon_{0,a}^2}- \sqrt{(\epsilon -\epsilon_{sh})2+\epsilon_{0,a}^2}\right)- C_{a2}\epsilon}{1- C_{a1} \epsilon_{0,a}}
$$

$$
\epsilon_{sh} = \frac{C_{a2}\epsilon_{0,a}}{\sqrt{C_{a1}^2-C_{a2}^2}}
$$

- **\( C \)**: Scaling constant
- **\( J_{\text{C}} \)**: Critical current density, representing the maximum current density a superconductor can carry without resistance.
- **\( B \)**: Magnetic field strength at the conductor
- **\( s(\epsilon) \)**: Strain scaling function, accounting for the effect of strain on superconducting properties.
- **\( t \)**: Reduced temperature, defined as \( T / T_{C0 \text{max}}^* \), where \( T \) is the temperature at the conductor
- **\( b \)**: Reduced magnetic field, defined as \( B / B_{C2}^* \), where \( B_{C2}^* \) is the upper critical field.
- **\( p \)**: Low field exponent of the pinning force ($\approx$ 0.5)
- **\( q \)**: High field exponent of the pinning force ($\approx$ 2)
- **\( B_{C2}^*(T,\epsilon) \)**: Upper critical field as a function of temperature and strain.
- **\( B_{C20 \text{max}}^* \)**: Maximum upper critical field at zero temperature and strain.
- **\( T_{C2}^*(T,\epsilon) \)**: Critical temperature as a function of temperature and strain.
- **\( T_{C0 \text{max}}^* \)**: Maximum critical temperature at zero field and strain.
- **\( \epsilon \)**: Strain applied to the superconductor.
- **\( \epsilon_{sh} \)**: Strain shift parameter, related to the intrinsic strain behavior.
- **\( \epsilon_{0,a} \)**: Residual strain component
- **\( C_{a1}, C_{a2} \)**: Material-specific strain fitting constants used in the strain scaling function.


--------------------

### NbTi

#### Old empirical | `jcrit_nbti()`

Parameters following the old empirical Lubell[^4] scaling are shown below:

- $B_{C20 \text{max}}^*: 15.0 \ \text{T}$
- $T_{C0 \text{max}}^*: 9.3 \ \text{K}$
- $C_0: 1 \times 10^{10}$

$$
B = \frac{B_{\text{conductor}}}{B_{C20 \text{max}}^*}
$$

$$
T_{\text{critical}} = T_{C0 \text{max}}^* \times (1-B)^{0.59}
$$

$$
T = 1- \frac{T_{\text{conductor}}}{T_{\text{critical}}}
$$

$$
J_{C} = C_0 \times (1-B)\times T
$$

The critical surface at zero strain looks as follows: 

<figure markdown>
![Old NbTi](./images/old_NbTi_zero_strain.png){ width = "100"}
<figcaption>Figure 2: Critical current density surface for old NbTi superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/old_NbTi_zero_strain.html)** :bar_chart:

!!! quote "Simon Chislett-McDonald and Prof. Damian P. Hampshire"
    *"The existing Nb-Ti critical surface model in PROCESS is based on an antiquated scaling law that has not been used in the superconductivity field for at least three decades. It is somewhat simplistic and linear in $J_C (B)$ (for a fixed temperature), lacking accuracy in both the high and low field regions: it under predicts J_c at low fields and over predicts both $B_{C20 \text{max}}^*$ and $J_C$ at high fields. The Nb-Ti strands that the law was fit to were also made to a different specification to those used in ITER. As such, both the law itself and the fit have been revised."*

------------


#### Durham model | `gl_nbti()`

Critical current density of the superconductor in an ITER
Nb-Ti strand based on the Ginzburg-Landau theory of superconductivity

This law is fitted to $J_c (B,T)$ measurements of five ITER specification Nb-Ti strands measured at Durham University. It accurately describes the behaviour of  $J_c (B,T)$ for all fields and temperatures. Additionally it incorporates the strain dependence of Nb-Ti’s $J_c$, as fit to a single filament of Nb-Ti. 

##### Derivation

The Durham scaling law is derived from the well-known equation for the volume pinning force $F_p$:

$$
F_p =  J_CB = A \frac{\left[B_{C2 }^*(T,\epsilon)\right]^n}{\left(\frac{2\pi h}{2\text{e}}\right)^{0.5}\mu_0 \left[\kappa_1^*(T,\epsilon)\right]^2} b^p(1-b)^q
$$

$\kappa_1^*$ is the effective Ginzburg-Landau parameter. The effective upper critical field, $B_{C2 }^*$ can be written in terms of $\kappa_1^*$ as:

$$
B_{C2 }^*(T,\epsilon) = \sqrt{2}\kappa_1^*(T,\epsilon)B_c(T,\epsilon)
$$

where $B_c$ is the thermodynamic critical field. From the two fluid model it is known:

$$
B_c(T) = B_c(0)(1-t^2)
$$

from the BCS equation $B_c (0,\epsilon) \propto T_c$ and extensive measurements have yielded:

$$
B_{c2}^*(T,\epsilon) = B_{c2}^*(0,\epsilon)(1-t^v)
$$

Substituting these equations into the $F_p$ equation gives:

$$
J_{\text{c,eng}}(B,T,\epsilon_I) = A^*(\epsilon_I)[T_c^*(\epsilon_I)(1-t^2)]^2[B_{c2}^*(\epsilon_I)(1-t^v)]^{n-3}b^{p-1}(1-b)^q
$$

where the introduced intrinsic strain $\epsilon_I$ and $A^*(\epsilon_I)$. The strain dependencies are related through:

$$
\frac{B_{c2}^*(0,\epsilon_I)}{B_{c2}^*(0,0)} = \left(\frac{T_c^*(\epsilon_I)}{T_c^*(0)}\right)^w = \left(\frac{A^*(\epsilon_I)}{A^*(0)}\right)^{\frac{w}{u}}
$$

By definition $\epsilon_a$ is the applied strain and $\epsilon_I = 0$ when J_c is maximum.

$$
\epsilon_I = \epsilon_a - \epsilon_m
$$

where $\epsilon_m$ is the applied strain at which $J_c$ is maximum. When $=0$, $A^*(\epsilon_I)$ is constant. The strain dependence of $B_{c2}^*$ can be specified as a polynomial given by:

$$
\frac{B_{c2}^*(0,\epsilon_I)}{B_{c2}^*(0,0)} = s(\epsilon_I) = 1+c_2\epsilon_I^2+c_3\epsilon_I^3+c_4\epsilon_I^4
$$

------------------

The materials data is from five ITER PF coil specification Nb-Ti strands. The Nb-Ti strands were produced by Chapetskiy Mechanical Plant (Glasov, Russia) for ITER PF6. Magnetisation measurements were performed on samples at 4.2 K. Transport measurements were performed at temperatures of 3.5 K, 4.0 K, 4.2 K, 5.0 K, 5.5 K, 6.0 K, 7.0 K and 8.0 K. 

The best fits of the Durham scaling law to the data calculated using the Python (version 3.7.4) `scipy.optimize.curve_fit` function. The strands’ diameters are 0.730 ± 0.005 mm, and they have a copper volume fraction of 69%. 

This produces the recommended variable values of:

- $A_{\text{non-Cu}}^*(0): 1102 \  \text{A}\text{T}^{3-n}\text{K}^{-2}$
- $B_{c2}^*(0,0): 14.9 \  \text{T}$
- $T_c^*(0): 9.0 \ \text{K}$
- $p: 0.49$
- $q: 0.56$
- $n: 1.83$
- $v: 1.42$
- $u: 0.0$
- $w: 2.2$

The stress dependance is derived by observing bespoke single-filament wires. In a multifilamentary Nb-Ti wire considerable strain would be applied when winding and mounting it to a Walter's spring. In such a sample, the filaments would be significantly compressed on the inboard side of the wire and significantly tensioned on the outboard side while soldering the wire to the spring. Such an effect is minimised by choice to measure a small single filament wire because the filament lies on the neutral axis. The critical current densities were measured at 4.2 K, at fields of 7 T to 10 T and between intrinsic strains of -1.03 % and +1.26%.  The strain fit parameters are shown in Table 5. 

The best fit parameters calculated below were done with Python (version 3.7.4) `curve_fit` function:

- $c_2 : -0.0025$
- $c_3 : -0.0003$
- $c_4 : -0.0001$
- $\epsilon_m : -2\times10^{-5}$

This gives a RMS error of 0.067 Amperes.

<figure markdown>
![Durham NbTi](./images/Durham_NbTi_zero_strain.png){ width = "100"}
<figcaption>Figure 3: Critical current density surface for Durham model NbTi superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/Durham_NbTi_zero_strain.html)** :bar_chart:

------------



## High-temperature superconductors

### Bi-2212 | `bi2212()`

Bi-2212 has demonstrated tolerating magnetic fields up to 34 T, surpassing the 24 T limit of Nb<sub>3</sub>Sn, can be achieved using a multifilament round wire conductor made of the high-temperature superconductor Bi₂Sr₂CaCu₂O₈₋ₓ (Bi-2212). Despite having many high-angle grain boundaries and lacking macroscopic texture, Bi-2212 attains high superconducting critical current densities ($J_C$) of 2500 A/mm² at 20 T and 4.2 K. Unlike REBa₂Cu₃O₇₋ₓ (REBCO) conductors, Bi-2212 does not require extreme texture, has a lower aspect ratio, and is less sensitive to defects.[^3]

This calculated critical current density and the temperature margin
for Bi-2212 superconductor is done using a fit by M. Kovari to measurements
described in the reference[^3], specifically from the points shown in Figure 6.

The validity of the fit below is for:

- $T < 20 \  \text{K}$
- $B > 6 \  \text{T}$

$$
B = \frac{B_{\text{conductor}}}{e^{-0.168\left(T_{\text{conductor}}-4.2\right)}}
$$

$$
J_C = f_{\text{strain}} \times \left(1.175\times 10^9 \times e^\left({-0.02115 \times B}\right)-1.288 \times 10^8\right)
$$

$B_{\text{conductor}}$ & $T_{\text{conductor}}$ are the field and temperature at the superconductor respectively. $f_{\text{strain}}$ can be used to account for strain, radiation damage, fatigue, or AC losses.



The critical surface at zero strain looks as follows: 

<figure markdown>
![Bi-2212](./images/Bi_2212_zero_strain.png){ width = "100"}
<figcaption>Figure 2: Critical current density surface for Bi-2212 superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/Bi_2212_zero_strain.html)** :bar_chart:


-----------

### REBCO

[^1]: L. Bottura and B. Bordini, “$J_{C}(B,T,\varepsilon)$ Parameterization for the ITER ${\rm Nb}_{3}{\rm Sn}$ Production,” IEEE Transactions on Applied Superconductivity, vol. 19, no. 3, pp. 1521-1524, Jun. 2009, doi: https://doi.org/10.1109/tasc.2009.2018278.
[^2]: V. Corato, “EUROFUSION WPMAG-REP(16) 16565 Common operating values for DEMO magnets design for 2016 REPORT.”Accessed: May 12, 2025. [Online].
Available: https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPMAGREP16_16565_submitted.pdf
[^3]: D. C. Larbalestier, J. Jiang, U. P. Trociewitz, F. Kametani, and E. E. Hellstrom, “A transformative superconducting magnet technology for fields well above 30 T using isotropic round wire multifilament Bi2Sr2CaCu2O8-x conductor,” May 06, 2013. https://www.researchgate.net/publication/236627864_A_transformative_superconducting_magnet_technology_for_fields_well_above_30_T_using_isotropic_round_wire_multifilament_Bi2Sr2CaCu2O8-x_conductor
[^4]: M. Lubell, “Empirical scaling formulas for critical current and critical field for commercial NbTi,” IEEE Transactions on Magnetics, vol. 19, no. 3, pp. 754-757, May 1983, doi: https://doi.org/10.1109/tmag.1983.1062311.