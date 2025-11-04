

`PROCESS` offers a variety of different superconducting coil types, both low-temperature superconductors (LTS) and high-temperature superconductors (HTS).


---------------

## Low-temperature superconductors

### Nb<sub>3</sub>Sn

#### Western Superconducting technologies Nb<sub>3</sub>Sn | `western_superconducting_nb3sn()`

The following parameters [^2] following the [Bottura scaling](#bottura-scaling-bottura_scaling) are shown below:

- $C_{\text{a1}}: 50.06$
- $C_{\text{a2}}: 0$
- $\epsilon_{\text{0,a}}: 0.312%$
- $\epsilon_{\text{m}}: -0.059%$
- $B_{C20 \text{max}}^* : 33.24 \ \text{T}$
- $T_{C0 \text{max}}^*:  16.34 \ \text{K}$
- $C: 83075 \ \text{AT/mm}^2$
- $p: 0.593$
- $q: 2.156$

The critical surface at zero strain looks as follows: 

<figure markdown>
![Wester Superconducting Ltd. Nb3Sn](./images/western_superconducting_Nb3Sn_zero_strain.png){ width = "100"}
<figcaption>Figure 1: Critical current density surface for the Wester Superconducting Ltd. Nb<sub>3</sub>Sn superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/western_superconducting_Nb3Sn_zero_strain.html)** :bar_chart:

!!! quote "Excerpt"
    *"TF Coil: The initial choice of EUTF4 right leg (OST) and its parametrization can be replaced by the actual WST strand (ø = 1.5 mm) used and tested in the DEMO TF prototypes manufactured in 2014 and measured in 2015 and 2016. The EUTF4 OST strand parameterisation can still be employed in the conductor design, if the designers prefer. The parameters are measured in the scope of the 2014 activities and reported below. Compared to the EUTF4 right leg (OST), the performance of the WST at low strain is superior by about 10%."*[^2]

----------------

#### EUTF4 Nb<sub>3</sub>Sn | `itersc()`



The following parameters [^2] following the [Bottura scaling](#bottura-scaling-bottura_scaling) are shown below:

- $C_{\text{a1}}: 44.48$
- $C_{\text{a2}}: 0$
- $\epsilon_{\text{0,a}}: 0.256%$
- $\epsilon_{\text{m}}: -0.110%$
- $B_{C20 \text{max}}^* : 32.97 \ \text{T}$
- $T_{C0 \text{max}}^*:  16.06 \ \text{K}$
- $C: 19922 \ \text{AT/mm}^2$
- $p: 0.63$
- $q: 2.10$

The critical surface at zero strain looks as follows: 

<figure markdown>
![EUTF4 Nb3Sn](./images/eutf4_Nb3Sn_zero_strain.png){ width = "100"}
<figcaption>Figure 2: Critical current density surface for the EUTF4 Nb<sub>3</sub>Sn superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/eutf4_Nb3Sn_zero_strain.html)** :bar_chart:

---------------

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
Nb-Ti strand based on the Ginzburg-Landau theory of superconductivity.

The model derivation can be found [here](#durham-ginzburg-landau-model-derivation)

This law is fitted to $J_c (B,T)$ measurements of five ITER specification Nb-Ti strands measured at Durham University. It accurately describes the behaviour of  $J_c (B,T)$ for all fields and temperatures. Additionally it incorporates the strain dependence of Nb-Ti’s $J_c$, as fit to a single filament of Nb-Ti. 


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

#### Bi-2212 | `bi2212()`

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


---------------


#### 2nd Generation REBCO | `jcrit_rebco()`

The following parameters following the scaling are shown below:

- $B_{C20 \text{max}}^*: 132.5 \ \text{T}$
- $T_{C0 \text{max}}^*: 90.0 \ \text{K}$
- $C: 1.82962 \times 10^8$
- $p: 0.5875$
- $q: 1.7$
- $\alpha: 1.54121$
- $\beta: 1.96679$

$$
B = B_{C20 \text{max}}^* \times \left(1-\frac{T_{\text{conductor}}}{T_{\text{C0max}}^*}\right)^{\alpha}
$$

$$
J_c = \frac{C}{B_{\text{conductor}}}B^{\beta} \times \left(\frac{B_{\text{conductor}}}{B}\right)^p\left(1-\frac{B_{\text{conductor}}}{B}\right)^q
$$

<figure markdown>
![2nd Generation REBCO](./images/2nd_gen_rebco_zero_strain.png){ width = "100"}
<figcaption>Figure 3: Critical current density surface for 2nd generation REBCO superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/2nd_gen_rebco_zero_strain.html)** :bar_chart:

-------------

#### Durham REBCO | `gl_rebco()`



The model derivation[^6] can be found [here](#durham-ginzburg-landau-model-derivation)


Values of scaling parameters $T_c(\epsilon=0), n, s, w, u, c_1, c_2, c_3, c_4$ and $\epsilon_m$ were taken from measurements undertaken in Durham[^7]:

- $T_c^*(0): 185 \ \text{K}$
- $n: 3.33$
- $u: 0.0$
- $w: 2.2$
- $s: 5.27$
- $c_2 : -0.0191$
- $c_3 : -0.0039$
- $c_4 : -0.0010$
- $\epsilon_m : 0.058$


Values of $A_{\text{Total-tape}}^*(0)$, $B_{c2}^*(0,0)$, $p$ and $q$ were calculated by fitting to state-of-the-art measurements shown by the National MagLab at Tallahassee [^5]:

- $A_{\text{Ttotal-tape}}^*(0): 295.0 \  \text{A} \text{m}^{-2}\text{T}^{3-n}\text{K}^{-2}$
- $B_{c2}^*(0,0): 429 \  \text{T}$
- $p: 0.32$
- $q: 2.50$


!!! warning "Model applicability"
    It should be noted that the values of $B_{c2}^*(0,0)$ and $T_c^*(0)$ here are not representative of the true REBCO material properties - rather these values yielded the best fit to the data within the data range $2 - 14 \ \text{T}$ and $4.2 - 60 \ \text{K}$ in Durham[^7], and $1 - 30 \ \text{T}$ at $4.2 \ \text{K}$ for the data published by Tallahassee[^5]. 

The strain fit is a polynomial, and as such produces un-physical results for large strains (both positive and compressive). Typically REBCO samples do not survive strains > + 0.7 %.

<figure markdown>
![Durham REBCO](./images/Durham_REBCO_zero_strain.png){ width = "100"}
<figcaption>Figure 4: Critical current density surface for Durham model REBCO superconductor as a function of magnetic field and temperature at the conductor.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/Durham_REBCO_zero_strain.html)** :bar_chart:

-------------

#### Hazelton-Zhain REBCO | `hijc_rebco()`

This model is based on the parametrization presented in Wolf et.al.[^8]

$$
I_c(B,T) = A\frac{B_0(T)^{\beta}}{B}\left(\frac{B}{B_0(T)}\right)^p\left(1- \frac{B}{B_0(T)}\right)^q
$$

$$
B_0(T) = B_0\left(1- \frac{T}{T_0}\right)
$$

The value of $A$ is transformed into a function $A(T)$ based on a Newton polynomial fit
considering $A(4.2 \ \text{K}) = 2.2 \times 10^8$, $A(20 \ \text{K}) = 2.3 \times 10^8$ and $A(65 \  \text{K}) = 3.5 \times 10^8.$

These values were selected manually. A good fit to the published data can be seen in the $4-10 \ \text{T}$ range but the fit deviates at very low or very high field.
Keep in mind that ITER's coils are measured against a critical current criterion of  $1\times 10^{-5} \  \text{V/m}$, while REBCO is measured against $1\times 10^{-4} \  \text{V/m}$.

$$
A(T) = A_0 + (uT^2)+vT
$$

The fitted coefficients are taken and fitted to both Hazelton et.al[^9] and the high $I_c$ parametrization fit adapted from Zhai et al[^10].

The fit values are:

- $T_c^*(0): 92 \ \text{K}$
- $B_{c2}^*(0,0): 138 \  \text{T}$
- $a: 1.4$
- $b: 2.005$
- $A_0: 2,2 \times 10^8$
- $p: 0.39$
- $q: 0.9$
- $u: 33450.0$
- $v: -176577.0$

The output of the model only gives the critical current ($I_c$) and not the critical current density ($J_c$). Therefore the value is multiplied by the tape paramters to find the critical current density in the REBCO strand

$$
J_c = I_c \frac{w_{\text{tape}} \times \Delta x_{\text{REBCO}}}{w_{\text{tape}} \times \Delta x_{\text{tape}} }
$$

where $w_{\text{tape}}$ is the width of the superconducting tape, $\Delta x_{\text{REBCO}}$ is the thickness of the REBCO layer and $\Delta x_{\text{tape}}$ is the full thickness of the superconducting tape.

The critical current density for a single REBCO strand can be seen below:


<figure markdown>
![Hazelton-Zhai REBCO](./images/Hazelton_Zhai_REBCO_zero_strain.png){ width = "100"}
<figcaption>Figure 4: Critical current density surface for Hazelton-Zhai model REBCO superconductor as a function of magnetic field and temperature at the conductor. The critical current is assumed only for the REBCO tape strand.</figcaption>
</figure>

:bar_chart: **An interactive version of the critical surface graph above can be found [here](./images/Hazelton_Zhai_REBCO_zero_strain.html)** :bar_chart:

-------------------------


### CroCo Cable Geometry | `calculate_croco_cable_geometry()`

The geometry of a single CroCo cable is calculated as follows:

1. The diameter of the circular internal tape region is given by the outer copper diameter minus its thickness:

    $$
    \overbrace{D_{\text{cable,internal}}}^{\texttt{dia_croco_strand_tape_region}} = \overbrace{D_{\text{cable}}}^{\texttt{dia_croco_strand}} - \overbrace{dx_{\text{cable,copper}}}^{\texttt{dx_croco_strand_copper}}
    $$

2. The total thickness of the HTS tape is found:

    $$
    \overbrace{dx_{\text{tape}}}^{\texttt{dx_hts_tape_total}} = \overbrace{dx_{\text{tape,REBCO}}}^{\texttt{dx_hts_tape_rebco}} + \overbrace{dx_{\text{tape,copper}}}^{\texttt{dx_hts_tape_copper}} + \overbrace{dx_{\text{tape,Hastelloy}}}^{\texttt{dx_hts_tape_hastelloy}}
    $$

3. The width of the tape is scaled to be:

    $$
    \overbrace{dr_{\text{tape}}^2}^{\texttt{dr_hts_tape}} = \frac{D_{\text{cable,internal}}\times 0.00375}{0.0054}
    $$

4. The height of the tape stack at the centre of the cable is found using the tape width and diameter of the region:

    $$
    \overbrace{dx_{\text{tape,stack}}}^{\texttt{dx_croco_strand_tape_stack}} = \sqrt{D_{\text{cable,internal}}^2-\overbrace{dr_{\text{tape}}^2}^{\texttt{dr_hts_tape}}}
    $$

5. The number of tapes in the stack is thus:

    $$
    \overbrace{N_{\text{cable,tapes}}}^{\texttt{n_croco_strand_hts_tapes}} = \frac{dx_{\text{tape,stack}}}{dx_{\text{tape}}} 
    $$

6. The total copper area in the strand (from the copper sheath and inside the tapes) is:

    $$
    \overbrace{A_{\text{cable,copper}}}^{\texttt{a_croco_strand_copper_total}} = \left(\left(\pi D_{\text{cable}} dx_{\text{cable,copper}}\right) - \pi dx_{\text{cable,copper}}^2\right)\\
     + \left(dx_{\text{tape,copper}}dr_{\text{tape}}N_{\text{cable,tapes}}\right)
    $$

7. The total Hastelloy area in the strand is:

    $$
    \overbrace{A_{\text{cable,Hastelloy}}}^{\texttt{a_croco_strand_hastelloy}} = dx_{\text{tape,Hastelloy}}dr_{\text{tape}}N_{\text{cable,tapes}}
    $$

8. The area of the solder surrounding the tape stack is:

    $$
    \overbrace{A_{\text{cable,solder}}}^{\texttt{a_croco_strand_solder}} = \frac{\pi}{4}D_{\text{cable,internal}}^2 - dx_{\text{tape,stack}}dr_{\text{tape}}
    $$

9. The total area of REBCO in the tape stack is:

    $$
    \overbrace{A_{\text{cable,REBCO}}}^{\texttt{a_croco_strand_hts_tapes}} = dx_{\text{tape,REBCO}}dr_{\text{tape}}N_{\text{cable,tapes}}
    $$

10. The total area of the cable is thus:

    $$
    \overbrace{A_{\text{cable}}}^{\texttt{a_croco_strand}} = \frac{\pi}{4}D_{\text{cable}}^2
    $$

------------------------

### Durham Ginzburg-Landau Model Derivation

The Durham scaling law [^6] is derived from the well-known equation for the volume pinning force $F_p$:

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

with the introduced intrinsic strain $\epsilon_I$ and $A^*(\epsilon_I)$. The strain dependencies are related through:

$$
\frac{B_{c2}^*(0,\epsilon_I)}{B_{c2}^*(0,0)} = \left(\frac{T_c^*(\epsilon_I)}{T_c^*(0)}\right)^w = \left(\frac{A^*(\epsilon_I)}{A^*(0)}\right)^{\frac{w}{u}}
$$

By definition $\epsilon_a$ is the applied strain and $\epsilon_I = 0$ when $J_c$ is maximum.

$$
\epsilon_I = \epsilon_a - \epsilon_m
$$

where $\epsilon_m$ is the applied strain at which $J_c$ is maximum. When $=0$, $A^*(\epsilon_I)$ is constant. The strain dependence of $B_{c2}^*$ can be specified as a polynomial given by:

$$
\frac{B_{c2}^*(0,\epsilon_I)}{B_{c2}^*(0,0)} = s(\epsilon_I) = 1+c_2\epsilon_I^2+c_3\epsilon_I^3+c_4\epsilon_I^4
$$



[^1]: L. Bottura and B. Bordini, “$J_{C}(B,T,\varepsilon)$ Parameterization for the ITER ${\rm Nb}_{3}{\rm Sn}$ Production,” IEEE Transactions on Applied Superconductivity, vol. 19, no. 3, pp. 1521-1524, Jun. 2009, doi: https://doi.org/10.1109/tasc.2009.2018278.
[^2]: V. Corato, “EUROFUSION WPMAG-REP(16) 16565 Common operating values for DEMO magnets design for 2016 REPORT.”Accessed: May 12, 2025. [Online].
Available: https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPMAGREP16_16565_submitted.pdf
[^3]: D. C. Larbalestier, J. Jiang, U. P. Trociewitz, F. Kametani, and E. E. Hellstrom, “A transformative superconducting magnet technology for fields well above 30 T using isotropic round wire multifilament Bi2Sr2CaCu2O8-x conductor,” May 06, 2013. https://www.researchgate.net/publication/236627864_A_transformative_superconducting_magnet_technology_for_fields_well_above_30_T_using_isotropic_round_wire_multifilament_Bi2Sr2CaCu2O8-x_conductor
[^4]: M. Lubell, “Empirical scaling formulas for critical current and critical field for commercial NbTi,” IEEE Transactions on Magnetics, vol. 19, no. 3, pp. 754-757, May 1983, doi: https://doi.org/10.1109/tmag.1983.1062311.
[^5]: N. High, “Plots - MagLab,” Nationalmaglab.org, 2018. https://nationalmaglab.org/magnet-development/applied-superconductivity-center/plots (accessed May 19, 2025).
[^6]: S B L Chislett-Mcdonald, Y. Tsui, E. Surrey, M. Kovari, and D. P. Hampshire, “The magnetic field, temperature, strain and angular dependence of the critical current density for Nb-Ti,” Journal of Physics Conference Series, vol. 1559, no. 1, pp. 012063–012063, Jun. 2020, doi: https://doi.org/10.1088/1742-6596/1559/1/012063.
[^7]: P. Branch, K. Osamura, and D. Hampshire, “Weak emergence in the angular dependence of the critical current density of the high temperature superconductor coated conductor REBCO,” Superconductor Science and Technology, vol. 33, no. 10, p. 104006, Sep. 2020, doi: https://doi.org/10.1088/1361-6668/abaebe.
[^8]: M. J. Wolf, Nadezda Bagrets, W. H. Fietz, C. Lange, and K.-P. Weiss, “Critical Current Densities of 482 A/mm2 in HTS CrossConductors at 4.2 K and 12 T,” IEEE Transactions on Applied Superconductivity, vol. 28, no. 4, pp. 1–4, Jun. 2018, doi: https://doi.org/10.1109/tasc.2018.2815767.
[^9]: D. W. Hazelton, “4th Workshop on Accelerator Magnets in HTS (WAMHTS-4) | 2G HTS Wire Development at SuperPower,”
Indico, 2017. https://indico.cern.ch/event/588810/contributions/2473740/ (accessed May 20, 2025).
[^10]: Y. Zhai, D. van der Laan, P. Connolly, and C. Kessel, “Conceptual design of HTS magnets for fusion nuclear science facility,” Fusion Engineering and Design, vol. 168, p. 112611, Jul. 2021, doi: https://doi.org/10.1016/j.fusengdes.2021.112611.
‌
‌
‌