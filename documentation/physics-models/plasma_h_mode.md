## L-H Power Threshold Scalings

Transitions from a standard confinement mode (L-mode) to an improved
confinement regime (H-mode), called L-H transitions, are observed in most
tokamaks.

A range of scaling laws are available that provide estimates of the
power terms required to initiate these transitions, via extrapolations
from present-day devices. PROCESS calculates these power threshold values
for the scaling laws listed [below](#l-h-scaling-options), in routine `l_h_threshold_power()`.

Depending on the value of the chosen scaling by setting `i_l_h_threshold`, a different L-H threshold power is set to the `p_l_h_threshold_mw` variable.

We define the net power across the seperatrix for the scaling as `p_plasma_separatrix_mw` below. This is equal to the net heating power of the plasma with radiation losses removed. This is then treated as the excess heating power for the plasma that is given to the divertors.

$$
\mathtt{p_plasma_separatrix_mw} = \frac{\mathrm{d}W}{\mathrm{d}t} =  \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}} - P_{\text{rad}}
$$

There are two separate constraint equations for enforcing the L-H threshold.

----------------

### Use the full divertor power

There are two constraints that can be used to enforce L-mode or H-mode.

Constraint 15 (`icc = 15`) is used to enforce H-mode by mandating that the power transported through the separatrix is greater than or equal to the L-H threshold power, a necessary condition for H-mode. `h_mode_threshold_margin >= 1.0` can be used to ensure the separatrix power exceeds the threshold by some margin. 

$$
\mathtt{p\_plasma\_separatrix\_mw} \ge \mathtt{h\_mode\_threshold\_margin} \times \underbrace{\mathtt{p\_l\_h\_threshold\_mw}}_{\text{Power from scaling}}
$$

For example, `h_mode_threshold_margin = 1.2` ensures that `p_plasma_separatrix_mw` is at least $1.2\times$ greater than the threshold power `p_l_h_threshold_mw`.


Constraint 22 (`icc = 22`) is the opposite of constraint 15 and ensures that the power transported through the separatrix is less than or equal to the L-H threshold power. `l_mode_threshold_margin >= 1.0` can be used to ensure that the threshold power is greater than the separatrix power by some margin.

$$
\underbrace{\mathtt{p\_l\_h\_threshold\_mw}}_{\text{Power from scaling}}  \ge \mathtt{l\_mode\_threshold\_margin} \times \mathtt{p\_plasma\_separatrix\_mw} 
$$

For example, `l_mode_threshold_margin = 1.2` ensures that `p_l_h_threshold_mw` is at least $1.2\times$ greater than the separatrix power `p_plasma_separatrix_mw`.


**It is recommended to use the H-mode constraint `icc = 15` at all times unless explicitly trying to enforce L-mode in which case `icc = 22` is used.**

-------

### Use the injected power reduced divertor power.

This constraint can be activated by stating `icc = 73` in the input file.

$$
1.0 - \frac{\mathtt{p_plasma_separatrix_mw}}{
\underbrace{\mathtt{p\_l\_h\_threshold\_mw}}_{\text{Power from scaling}}+ P_{\text{HCD}}}
$$


--------------------

## L-H scaling options

----------------

### ITER-1996 Scalings

The general form is:

$$
P_{\text{L-H}} = 0.45 (0.6\bar{n}_{e,20}R^2)^{\alpha} \bar{n}^{0.75}_{\text{e},20}B_{\text{T}}R^2
$$

where $\alpha$ lies in the range of $-0.25 \le \alpha \le 0.25$,  $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla and $R$ is the plasma major radius in metres.



------------------

#### ITER-1996 Nominal Scaling

Is selected with `i_l_h_threshold = 1` [^1] [^2]

$$
P_{\text{L-H}} = 0.45 \times \bar{n}^{0.75}_{\text{e},20}B_{\text{T}}R^2
$$

---------------

#### ITER-1996 Upper Scaling

Is selected with `i_l_h_threshold = 2` [^1] [^2]

$$
P_{\text{L-H}} = 0.3960502816 \times \bar{n}_{\text{e},20}B_{\text{T}}R^{2.5}
$$

---------------

#### ITER-1996 Lower Scaling

Is selected with `i_l_h_threshold = 3` [^1] [^2]

$$
P_{\text{L-H}} = 0.5112987149 \times \bar{n}_{\text{e},20}^{0.5}B_{\text{T}}R^{1.5}
$$

---------------


### Snipes 1997 ITER Scaling I

Is selected with `i_l_h_threshold = 4` [^3]

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$

$$
P_{\text{L-H}} = 0.65 \bar{n}_{\text{e},20}^{0.93} B_{\text{T}}^{0.86} R^{2.15}
$$

---------------

###  Snipes 1997 ITER Scaling II

Is selected with `i_l_h_threshold = 5` [^3]

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$


$$
P_{\text{L-H}} = 0.42 \bar{n}_{\text{e},20}^{0.8} B_{\text{T}}^{0.9} R^{1.99} \kappa^{0.76}
$$

---------------

### Martin 2008 Scalings

The general form is:

$$
P_{\text{L-H}} = 0.0488 e^{\pm0.057} \bar{n}_{\text{e},20}^{0.717 \pm 0.035} B_{\text{T}}^{0.803 \pm 0.032} S_{\text{p}}^{0.941 \pm 0.019}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla and $S_{\text{p}}$ is the plasma surface area in $\text{m}^2$.

We apply a mass-correction term to the scaling, stated by Martin et.al [^4] as per:

!!! quote "Mass dependence on threshold"

    "*It is also found in JET that the 
    threshold power in tritium discharges becomes further lower. The dependence of the threshold power 
    on the ion mass number M was roughly given by $P_{\text{L-H}} \propto \frac{1}{M}$ [^9]. When this mass dependence is 
    applied to the deuterium-tritium discharges for ITER, the above predicted values of PThresh can be 
    reduced by ~ 20%.*[^4]"

We thus apply a factor of $\left(\frac{2}{M_{\text{i}}}\right)$ to the end of the scalings, where $M_{\text{i}}$ is the average atomic mass of all ions. Therefore for a pure 50:50 D-T plasma giving a $M_{\text{i}} = 2.5$ the value of $P_{\text{L-H}}$ is dropped by 20%.

------------------

#### Martin 2008 Nominal Scaling

Is selected with `i_l_h_threshold = 6` [^4]

$$
P_{\text{L-H}} = 0.0488 \bar{n}_{\text{e},20}^{0.717} B_{\text{T}}^{0.803} S_{\text{p}}^{0.941}\left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Martin 2008 Upper Scaling

Is selected with `i_l_h_threshold = 7` [^4]

$$
P_{\text{L-H}} = 0.05166240355 \times \bar{n}_{\text{e},20}^{0.752} B_{\text{T}}^{0.835} S_{\text{p}}^{0.96}\left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Martin 2008 Lower Scaling

Is selected with `i_l_h_threshold = 8` [^4]

$$
P_{\text{L-H}} = 0.04609619059 \times \bar{n}_{\text{e},20}^{0.682} B_{\text{T}}^{0.771} S_{\text{p}}^{0.922}\left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

### Snipes 2000 Scalings

!!! quote "Excerpt fom Snipes et.al"
    *"A number of regression fits were performed on the full data set of all 10 tokamaks for all
    L–H transition points that fit the standard SELDB2 criteria in deuterium plasmas ($N$ = 702).
    A small improvement in the RMSE was obtained by correcting the total ICRF power from
    Alcator C-Mod with density- and toroidal-field-dependent corrections for the absorbed ICRF
    power based on experimental measurements. For hydrogen minority heating between 5 and
    6 $T$, $P_{\text{abs}} = 0.9\bar{n}_{\text{e}}^{-0.6} P_{\text{ICRH}}$ while for He3 minority heating above 6 T, the assumption is made
    that $P_{\text{abs}} = 0.75 P_{\text{ICRH}}$."*

- This scaling has a RMSE of 26.8%

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$


The general form is:

$$
P_{\text{L-H}} = 1.42\pm 0.127 \times \bar{n}_{\text{e},20}^{0.58 \pm 0.035} B_{\text{T}}^{0.82 \pm 0.031} R^{1.00 \pm 0.089} a^{0.81 \pm 0.066}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla, $R$ is the plasma major radius in metres and $a$ is the plasma major radius in metres.

We apply the same mass-correction factor given in the [Martin 2008 scalings above](#martin-2008-scalings).[^4] [^9] Snipes et.al provides the same justification [^9].

We thus apply a factor of $\left(\frac{2}{M_{\text{i}}}\right)$ to the end of the scalings, where $M_{\text{i}}$ is the average atomic mass of all ions. Therefore for a pure 50:50 D-T plasma giving a $M_{\text{i}} = 2.5$ the value of $P_{\text{L-H}}$ is dropped by 20%.

------------------

#### Snipes 2000 Nominal Scaling

Is selected with `i_l_h_threshold = 9` [^5]

$$
P_{\text{L-H}} = 1.42 \times \bar{n}_{\text{e},20}^{0.58} B_{\text{T}}^{0.82} R^{1.00} a^{0.81} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Snipes 2000 Upper Scaling

Is selected with `i_l_h_threshold = 10`  [^5]

$$
P_{\text{L-H}} = 1.547 \times \bar{n}_{\text{e},20}^{0.615} B_{\text{T}}^{0.851} R^{1.089} a^{0.876} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Snipes 2000 Lower Scaling

Is selected with `i_l_h_threshold = 11`  [^5]

$$
P_{\text{L-H}} = 1.293 \times \bar{n}_{\text{e},20}^{0.545} B_{\text{T}}^{0.789} R^{0.911} a^{0.744} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

### Snipes 2000 Closed Divertor Scalings

!!! quote "Excerpt fom Snipes et.al"
    *"Several machines have reported changes in the H-mode threshold due to divertor geometry modifications. ASDEX-Upgrade saw a 15% increase with a more-closed divertor, attributed to higher edge densities. Conversely, JET and JT-60U experienced a 20% decrease after installing more-closed divertors. Alcator C-Mod, with its inherently closed divertor, showed no change in threshold despite further closure of bypass gaps. Preliminary experiments with C-Mod's divertor bypass indicate no threshold change, even with a significant drop in divertor neutral pressure. These variations contribute to the scatter in H-mode threshold regression fits."*

    *"Only four tokamaks have closed-divertor data in the database (Alcator C-Mod, ASDEX-Upgrade, JET, and JT-60U). Although the resulting data set is limited ($N$ = 169)"*

- This scaling has a RMSE of 22.0%

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$


The general form is:

$$
P_{\text{L-H}} = 0.8\pm 0.067 \times \bar{n}_{\text{e},20}^{0.50 \pm 0.061} B_{\text{T}}^{0.53 \pm 0.058} R^{1.51 \pm 0.077}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla and $R$ is the plasma major radius in metres

We apply the same mass-correction factor given in the [Martin 2008 scalings above](#martin-2008-scalings).[^4] [^9] Snipes et.al provides the same justification [^9].

We thus apply a factor of $\left(\frac{2}{M_{\text{i}}}\right)$ to the end of the scalings, where $M_{\text{i}}$ is the average atomic mass of all ions. Therefore for a pure 50:50 D-T plasma giving a $M_{\text{i}} = 2.5$ the value of $P_{\text{L-H}}$ is dropped by 20%.

------------------

#### Snipes 2000 Closed Divertor Nominal Scaling

Is selected with `i_l_h_threshold = 12`  [^5]

$$
P_{\text{L-H}} = 0.8 \times \bar{n}_{\text{e},20}^{0.50} B_{\text{T}}^{0.53} R^{1.51} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Snipes 2000 Closed Divertor Upper Scaling

Is selected with `i_l_h_threshold = 13`  [^5]

$$
P_{\text{L-H}} = 0.867 \times \bar{n}_{\text{e},20}^{0.561} B_{\text{T}}^{0.588} R^{1.587} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

#### Snipes 2000 Closed Divertor Lower Scaling

Is selected with `i_l_h_threshold = 14`  [^5]

$$
P_{\text{L-H}} = 0.733 \times \bar{n}_{\text{e},20}^{0.439} B_{\text{T}}^{0.472} R^{1.433} \left(\frac{2}{M_{\text{i}}}\right)
$$

---------------

###  Hubbard 2012 L-I Scalings

The general form is:

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.94\pm 0.24}\bar{n}_{\text{e},20}^{0.65\pm 0.18}
$$

where $I_{\text{p}}$ is the plasma current in $\text{MA}$ and $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20}$.

------------------

#### Hubbard 2012 L-I Nominal Scaling

Is selected with `i_l_h_threshold = 15` [^6]

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.94}\bar{n}_{\text{e},20}^{0.65}
$$

---------------

#### Hubbard 2012 L-I Lower Scaling

Is selected with `i_l_h_threshold = 16` [^6]

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.7}\bar{n}_{\text{e},20}^{0.47}
$$

---------------

#### Hubbard 2012 L-I Upper Scaling

Is selected with `i_l_h_threshold = 17` [^6]

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{1.18}\bar{n}_{\text{e},20}^{0.83}
$$

---------------

###  Hubbard 2017 L-I Scaling

Is selected with `i_l_h_threshold = 18` [^7]

$$
P_{\text{L-H}} = 0.162 \times B_{\text{T}}^{0.26}\bar{n}_{\text{e},20} S_{\text{p}}
$$

where $B_{\text{T}}$ is the toroidal magnetic filed in $\text{T}$, $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20}$ and $S_{\text{p}}$ is the plasma surface area in $\text{m}^2$.

------------------

### Martin 2008 Aspect ratio corrected scalings

The general form is the same as the original [Martin 2008](#martin-2008-aspect-ratio-corrected-scalings) scaling with an aspect ratio correction factor from T. Takizuka et.al [^8]:

$$
P_{\text{L-H}} = 0.0488 e^{\pm0.057} \bar{n}_{\text{e},20}^{0.717 \pm 0.035} B_{\text{T}}^{0.803 \pm 0.032} S_{\text{p}}^{0.941 \pm 0.019} \\ 
\times  \left[0.098 \times \frac{A}{1.0 - \left(\frac{2.0}{(1.0 + A)}\right)^{0.5}}\right] \text{for} \  A \le 2.7
$$



where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla, $S_{\text{p}}$ is the plasma surface area in $\text{m}^2$ and $A$ is the plasma aspect ratio

We apply the same mass-correction done for the original scaling [discussed above](#martin-2008-scalings).

------------------

#### Martin 2008 Aspect Corrected Nominal Scaling

Is selected with `i_l_h_threshold = 19` [^4] [^8]

$$
P_{\text{L-H}} = 0.0488 \bar{n}_{\text{e},20}^{0.717} B_{\text{T}}^{0.803} S_{\text{p}}^{0.941}\left(\frac{2}{M_{\text{i}}}\right) \\ 
\times  \left[0.098 \times \frac{A}{1.0 - \left(\frac{2.0}{(1.0 + A)}\right)^{0.5}}\right] \text{for} \  A \le 2.7
$$

---------------

#### Martin 2008 Aspect Corrected Upper Scaling

Is selected with `i_l_h_threshold = 20` [^4] [^8]

$$
P_{\text{L-H}} = 0.05166240355 \times \bar{n}_{\text{e},20}^{0.752} B_{\text{T}}^{0.835} S_{\text{p}}^{0.96}\left(\frac{2}{M_{\text{i}}}\right) \\ 
\times  \left[0.098 \times \frac{A}{1.0 - \left(\frac{2.0}{(1.0 + A)}\right)^{0.5}}\right] \text{for} \  A \le 2.7
$$

---------------

#### Martin 2008 Aspect Corrected Lower Scaling

Is selected with `i_l_h_threshold = 21` [^4] [^8]

$$
P_{\text{L-H}} = 0.04609619059 \times \bar{n}_{\text{e},20}^{0.682} B_{\text{T}}^{0.771} S_{\text{p}}^{0.922}\left(\frac{2}{M_{\text{i}}}\right) \\ 
\times  \left[0.098 \times \frac{A}{1.0 - \left(\frac{2.0}{(1.0 + A)}\right)^{0.5}}\right] \text{for} \  A \le 2.7
$$


[^1]: T. Takizuka and International Atomic Energy Agency, Vienna (Austria),"Threshold power and energy confinement for ITER". 1996.
[^2]: J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,” Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997, doi: https://doi.org/10.1063/1.872406.
[^3]: J. A. Snipes and the ITER H-mode Threshold Database Working Group, "An Analysis of the H-mode Threshold in ITER," Controlled Fusion and Plasma Physics, 24th EPS Conference, Berchtesgaden, June 9th-13th 1997, vol.21A, part III, p.961. url:https://library.ipp.mpg.de/EPS_24_Vol3_1997.pdf.
[^4]: Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,” Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008, doi: https://doi.org/10.1088/1742-6596/123/1/012033.
[^5]: J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,” Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.
[^6]: A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012, doi: https://doi.org/10.1088/0029-5515/52/11/114009.
[^7]: A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,” Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.
[^8]: T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of the H-mode power threshold in tokamaks of the ITPA database,” Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233, Apr. 2004, doi: https://doi.org/10.1088/0741-3335/46/5a/024.
[^9]: E. Righi et al., “Isotope scaling of the H mode power threshold on JET,” Nuclear Fusion, vol. 39, no. 3, pp. 309–319, Mar. 1999, doi: https://doi.org/10.1088/0029-5515/39/3/302.
‌