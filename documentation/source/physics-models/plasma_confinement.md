# Plasma confinement time

## Overview

Confinement time scalings are empirical relationships derived from experimental data across various fusion machines. These scalings help predict how changes in tokamak parameters (like size, magnetic field strength, and plasma density) will affect the confinement time and overall performance.

The energy confinement time $\tau_E$ is calculated using one of a choice of empirical scalings. ($\tau_E$ is defined below.)

Normally most confinement scalings will be of the form:

$$
\tau_{\text{E}} =  C I_{\text{p}}^{\alpha_{I}} B_{\text{T}}^{\alpha_{B}} \overline{n}_e^{\alpha_{n}} P_{\text{L}}^{\alpha_{P}} R^{\alpha_{R}} \kappa^{\alpha_{\kappa}} \epsilon^{\alpha_{\epsilon}} M^{\alpha_{M}}
$$

Where $\tau_{\text{E}}$ is the confinement time in seconds, $C$ is a coefficient , $I_{\text{p}}$ [MA] is the plasma current, $B_{\text{T}}$ [T] is the toroidal magnetic field, $\overline{n}_e$ [$10^{19} \text{m}^{-3}$] is the electron central line averaged density, $P_{\text{L}}$ [MW] is the loss power, $R$ [m] is the major radius, $\kappa$ is the plasma elongation,
$\epsilon$ is the inverse aspect ratio and $M$ is the average atomic mass of the plasma.

Classically the loss power, $P_{\text{L}}$ is defined as:

$$
P_{\text{L}} = \frac{W}{\tau_{\text{E}}}
$$

where $W$ is the total thermal energy of the plasma. We can look at it mainly as the difference in heating and loss powers in the plasma, as such we interpret it as power transported out from the
“core” by charged particles. This leads to the classic definition of loss power for the scaling:

$$
P_{\text{L}} = \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}}
$$

where $f_{\alpha}$ is the [fraction of alpha power that is coupled to the plasma](../physics-models/fusion_reactions/plasma_reactions.md#coupled-alpha-particle-power), $P_{\alpha}$ is the alpha power, $P_{\text{c}}$ is the charged particle power, $P_{\text{OH}}$ is the ohmic heating power, $P_{\text{HCD}}$ is the plasma heating done by the external heating & current drive systems.

----------

## Calculating plasma confinement time | `calculate_confinement_time()`

The correspoding plasma confinement time is calculated by the `calculate_confinement_time()` function in `physics.py` with scalings taken from `confinement_time.py`.

A key definition of elongation is defined here and is used mainly in the ITER physics basis scalings [^12]:

$$
\kappa_{\text{IPB}} = \frac{V_{\text{p}}}{2\pi R}\frac{1}{\pi a^2}
$$

where $V_{\text{p}}$ is the plasma volume, $R$ is the plasma major radius and $a$ is the plasma minor radius.

The loss power $P_{\text{L}}$ [$\mathtt{p\_plasma\_loss\_mw}$] is calculated from above but may have a separate radiation term depending on the condition of `i_rad_loss` switch below.


-------------

### Effect of radiation on energy confinement

Published confinement scalings are all based on low radiation pulses. A power
plant will certainly be a high radiation machine, both in the core, due to
bremsstrahlung and synchrotron radiation, and in the edge due to impurity
seeding. The scaling data does not predict this radiation [^24] [^25], that needs to be
done by the radiation model. However, if the transport is very "stiff", as
predicted by some models, then the additional radiation causes an almost equal
drop in power transported by ions and electrons, leaving the confinement
nearly unchanged.

To allow for these uncertainties, three options are available, using the switch
`i_rad_loss`.

- For `i_rad_loss = 0` the total plasma radiation is taken from the loss power.

$$
P_{\text{L}} = \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}} - P_{\text{rad}}
$$

- For `i_rad_loss = 1` the plasma radiation only from the "core" region is taken from the loss power.

$$
P_{\text{L}} = \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}} - P_{\text{rad,core}}
$$

- For `i_rad_loss = 2` the plasma radiation is not taken from the loss power

$$
P_{\text{L}} = \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}}
$$

----------

### Ignition

Switch `i_plasma_ignited` can be used to denote whether the plasma is ignited, i.e. fully self-sustaining 
without the need for any injected auxiliary power during the burn. If `i_plasma_ignited = 1`, the heating and current drive power $P_{\text{HCD}}$, does not contribute to the loss power term. 

If `i_plasma_ignited = 0`, the plasma is not ignited, and the heating and current drive power $P_{\text{HCD}}$, does contribute to the loss power term.
phase. An ignited plasma will be difficult to control and is unlikely to be practical. This
option is not recommended.

----------

### Available confinement time scalings

Many energy confinement time scaling laws are available within PROCESS, for conventional aspect ratio tokamaks, spherical tokamaks, and stellarators.
The value of `i_confinement_time` determines which of the scalings is used in the plasma energy balance calculation.

The scaling chosen with `i_confinement_time` is then calculated and multiplied with the $H$-factor [$\mathtt{hfact}$]. $\mathtt{hfact}$ can be set as an interation variable by setting `ixc = 10` in the `IN.DAT` input file.

---------------

#### 0: User input confinement time

Is selected with `i_confinement_time = 0`

$$
\tau_{\text{E}} = \mathtt{t\_electron\_confinement\_in}
$$

------------

#### 1: Nec-Alcator scaling (Ohmic)

Is selected with `i_confinement_time = 1`[^1]

$$
\tau_{\text{E}} = 0.07 n_{20}aRq_{\text{cyl}}
$$

------------

#### 2: Mirnov scaling (H-mode)

Is selected with `i_confinement_time = 2`[^1]

$$
\tau_{\text{E}} = 0.2 a \sqrt{\kappa_{95}}I_{\text{p}}
$$

------------

#### 3: Merezhkin-Mukhovatov scaling (Ohmic / L-mode)

Is selected with `i_confinement_time = 3`[^1]

$$
\tau_{\text{E}} = 0.0035 \overline{n}_{20}a^{0.25}R^{2.75}q_{\text{cyl}}\kappa_{95}^{0.125}M_i^{0.5}T_{10}^{0.5}
$$


---------------

#### 4: Shimomura scaling (H-mode)

Is selected with `i_confinement_time = 4`[^1]

$$
\tau_{\text{E}} = 0.045 Ra B_{\text{T}}\sqrt{\kappa_{95}}\sqrt{M_{\text{i}}}
$$

----------------

#### 5: Kaye-Goldston scaling (L-mode)

Is selected with `i_confinement_time = 5`[^1]

$$
\tau_{\text{E}} = 0.055 I_{\text{p}}^{1.24}P_{\text{L}}^{-0.58}R^{1.65}a^{-0.49}\kappa_{95}^{0.28}n_{20}^{0.26}B_{\text{T}}^{-0.09}\left(\frac{M_{\text{i}}}{1.5}\right)^{0.5}
$$

----------------

#### 6: ITER 89-P scaling (L-mode)

Is selected with `i_confinement_time = 6`[^1] [^2]

$$
\tau_{\text{E}} = 0.048 I_{\text{p}}^{0.85}R^{1.2}a^{0.3}\kappa^{0.5}\overline{n}_{20}^{0.1}B_{\text{T}}^{0.2}M_{\text{i}}^{0.5} P_{\text{L}}^{-0.5}
$$

----------------

#### 7: ITER 89-0 scaling (L-mode)

Is selected with `i_confinement_time = 7` [^2]

$$
\begin{aligned}
\tau_E= & 0.04 I_{\text{p}}^{0.5} R^{0.3} a^{0.8} \kappa^{0.6} M_i^{0.5} \\
& +0.064 I_{\text{p}}^{0.8} R^{1.6} a^{0.6} \kappa^{0.2} \bar{n}_{20}^{0.6} B_0^{0.35} M_i^{0.2} / P_{\text{L}}
\end{aligned}
$$

----------------

#### 8: Rebut-Lallia scaling (L-mode)

Is selected with `i_confinement_time = 8` [^2]

$$
\begin{aligned}
\tau_E= & 1.65\left[1.2 \times 10^{-5} I_{\text{p}} \ell^{1.5} Z_{e f f}^{-0.5}\right. \\
& \left.+0.146 \bar{n}_{20}^{0.75} I_{\text{p}}^{0.5} B_0^{0.5} \ell^{2.75} Z_{e f f}^{0.25} / P_{\text{L}}\right]\left(A_i / 2\right)^{0.5}
\end{aligned}
$$

where $\ell = \left(a^2R\kappa\right)^{\frac{1}{3}}$

----------------

#### 9: Goldston scaling (L-mode)

Is selected with `i_confinement_time = 9` [^1]

$$
\tau_{\text{E}} = 0.037 I_{\text{p}} P_{\text{L}}^{-0.5} R^{1.75}a^{-0.37}\kappa_{95}^{0.5} \left(\frac{M_i}{1.5}\right)^{0.5}
$$

----------------

#### 10: T-10 scaling (L-mode)

Is selected with `i_confinement_time = 10` [^1]


$$
\tau_{\text{E}} =  0.095 a R B_{\text{T}} \kappa_{95}^{0.5} \frac{\overline{n}_{20}}{\overline{n}_{20*}}P_{\text{L}}^{-0.4} \left[\frac{Z_{\text{eff}}^2 I_{\text{p}}^4}{aRq_{\text{cyl}}^3\kappa_{95}^{1.5}}  \right]^{0.08}
$$

where $\overline{n}_{20*} = 1.3\left(\frac{B_{\text{T}}}{Rq_{\text{cyl}}}\right)$ and $\frac{\overline{n}_{20}}{\overline{n}_{20*}} \le 1$

----------------

#### 11: JAERI / Odajima-Shimomura scaling (L-mode)

Is selected with `i_confinement_time = 11` [^1]


$$
\tau_{\text{E}} =  \left[\frac{0.085\kappa_{95}a^2+0.069In_{20}^{0.6}B_{\text{T}}^{0.2}R^{1.6} a^{0.4} \kappa_{95}^{0.2} G\left(q_{\text{cyl}},Z_{\text{eff}}\right)}{P_{\text{L}}}\right]M_{\text{i}}^{0.5}
$$

where $G\left(q_{\text{cyl}},Z_{\text{eff}}\right) = Z_{\text{eff}}^{0.4}\left[\frac{\left(15 - Z_{\text{eff}}\right)}{20}\right]^{0.6}\left[3q_{\text{cyl}}\frac{q_{\text{cyl}}+5}{(q_{\text{cyl}}+2)(q_{\text{cyl}}+7)}\right]^{0.6}$


----------------

#### 12: Kaye "big" scaling (L-mode)

Is selected with `i_confinement_time = 12` [^1]

$$
\tau_{\text{E}} =  0.1051 I_{\text{p}}^{0.85} P_{\text{L}}^{-0.5} R^{0.5} a^{0.3} \kappa^{0.25} n_{20}^{0.1}B_{\text{T}}^{0.3}M_{\text{i}}^{0.5}
$$

-------------------------

#### 13: ITER H90-P scaling (H-mode)

Is selected with `i_confinement_time = 13` [^2]

$$
\tau_{\text{E}} =  0.064 I_{\text{p}}^{0.87} R^{1.82} a^{-0.12} \kappa_{95}^{0.35} \overline{n}_{20}^{0.09} B_{\text{T}}^{0.15} M_{\text{i}}^{0.5} P_{\text{L}}^{-0.5}
$$

-------------------------

#### 14: Minimum of ITER 89-P and ITER 89-O

Is selected with `i_confinement_time = 14` [^1] [^2]

Will return the value of [ITER 89-P](#6-iter-89-p-l-mode-scaling) or [ITER 89-O](#7-iter-89-0-l-mode-scaling), whichever is smaller.

-------------------------

#### 15: Riedel scaling (L-mode)

Is selected with `i_confinement_time = 15` [^2]

$$
\tau_{\text{E}} =  0.044 I_{\text{p}}^{0.93} R^{1.37} a^{-0.049} \kappa_{95}^{0.588} \overline{n}_{20}^{0.078} B_{\text{T}}^{0.152} P_{\text{L}}^{-0.537}
$$

-------------------------

#### 16: Christiansen scaling (L-mode)

Is selected with `i_confinement_time = 16` [^2]

$$
\tau_{\text{E}} =  0.24 I_{\text{p}}^{0.79} R^{0.56} a^{1.46} \kappa_{95}^{0.73} \overline{n}_{20}^{0.41} B_{\text{T}}^{0.29} P_{\text{L}}^{-0.79} M_{\text{i}}^{-0.02}
$$

-------------------------

#### 17: Lackner-Gottardi scaling (L-mode)

Is selected with `i_confinement_time = 17` [^2]

$$
\tau_{\text{E}} =  0.12 I_{\text{p}}^{0.8} R^{1.8} a^{0.4} \left(\frac{\kappa_{95}}{\left(1+\kappa_{95}\right)^{0.8}}\right) \overline{n}_{20}^{0.6} \hat{q}^{0.4} P_{\text{L}}^{-0.6}
$$

where $\hat{q} = \frac{(1+\kappa_{95}a^2B_{\text{T}})}{0.4 I_{\text{p}} R}$

-------------------------

#### 18: Neo-Kaye scaling (L-mode)

Is selected with `i_confinement_time = 18` [^2]

$$
\tau_{\text{E}} =  0.063 I_{\text{p}}^{1.12} R^{1.3} a^{-0.04} \kappa_{95}^{0.28} \overline{n}_{20}^{0.14} B_{\text{T}}^{0.04} P_{\text{L}}^{-0.59}
$$

-------------------------

#### 19: Riedel scaling (H-mode)

Is selected with `i_confinement_time = 19` [^2]

$$
\tau_{\text{E}} =  0.1 M_{\text{i}}^{0.5} I_{\text{p}}^{0.884} R^{1.24} a^{-0.23} \kappa_{95}^{0.317} \overline{n}_{20}^{0.105} B_{\text{T}}^{0.207} P_{\text{L}}^{-0.486}
$$

-------------------------

#### 20: Amended ITER H90-P scaling (H-mode)

Is selected with `i_confinement_time = 20` [^3]

$$
\tau_{\text{E}} =  0.082 M_{\text{i}}^{0.5} I_{\text{p}}^{1.02} R^{1.6}  \kappa_{95}^{-0.19}  B_{\text{T}}^{0.15} P_{\text{L}}^{-0.47}
$$

-------------------------

#### 21: Sudo et al. scaling (Stellarator)

Is selected with `i_confinement_time = 21` [^4]

$$
\tau_{\text{E}} =  0.17 P_{\text{L}}^{-0.58} \overline{n}_{20}^{0.69} B^{0.84} a^{2.0} R^{0.75}
$$

-------------------------

#### 22: Gyro reduced Bohm scaling (Stellarator)

Is selected with `i_confinement_time = 22` [^5]

$$
\tau_{\text{E}} =  0.25 P_{\text{L}}^{-0.6} \overline{n}_{20}^{0.6} B_{\text{T}}^{0.8} a^{2.4} R^{0.6}
$$

-------------------------

#### 23: Lackner-Gottardi scaling (Stellarator)

Is selected with `i_confinement_time = 23` [^6]

$$
\tau_{\text{E}} =  0.17 P_{\text{L}}^{-0.6} \overline{n}_{20}^{0.6} B_{\text{T}}^{0.8} a^{2.0} R q_{95}^{0.4}
$$

-------------------------

#### 24: ITER H93 ELM-free scaling (H-mode)

Is selected with `i_confinement_time = 24` [^7]

$$
\tau_{\text{E}} =  0.036  I_{\text{p}}^{1.06} B_{\text{T}}^{0.32} P_{\text{L}}^{-0.67} R^{1.79} \epsilon^{-0.11} \kappa^{0.66} \overline{n}_{20}^{0.17} M_{\text{i}}^{0.41}
$$

-------------------------

#### 25: TITAN Reversed-Field_Pinch scaling

Is selected with `i_confinement_time = 25`

!!! warning 
    This scaling has been removed

-------------------------

#### 26: ITER H-97P ELM-free scaling (H-mode)

Is selected with `i_confinement_time = 26` [^8]

$$
\tau_{\text{E}} =  0.031 M_{\text{i}}^{0.42} I_{\text{p}}^{0.95} R^{1.92} \epsilon^{0.08} \kappa_{95}^{0.63} \overline{n}_{19}^{0.35} B_{\text{T}}^{0.25} P_{\text{L}}^{-0.67}
$$

-------------------------

#### 27: ITER H-97P ELMy scaling (H-mode)

Is selected with `i_confinement_time = 27` [^8] [^9]

$$
\tau_{\text{E}} =  0.029 M_{\text{i}}^{0.2} I_{\text{p}}^{0.9} R^{2.03} \epsilon^{-0.19} \kappa_{95}^{0.92} \overline{n}_{19}^{0.4} B_{\text{T}}^{0.20} P_{\text{L}}^{-0.66}
$$

-------------------------

#### 28: ITER-96P (ITER-97L) scaling (L-mode)

Is selected with `i_confinement_time = 28` [^10]

$$
\tau_{\text{E}} =  0.023 M_{\text{i}}^{0.2} I_{\text{p}}^{0.96} R^{1.83} \epsilon^{-0.06} \kappa_{95}^{0.64} \overline{n}_{19}^{0.4} B_{\text{T}}^{0.03} P_{\text{L}}^{-0.73}
$$

-------------------------

#### 29: Valovic modified ELMy scaling (H-mode)

Is selected with `i_confinement_time = 29`

$$
\tau_{\text{E}} =  0.067 M_{\text{i}}^{0.05} I_{\text{p}}^{0.9} R^{1.31} \kappa^{0.56} \overline{n}_{19}^{0.45} B_{\text{T}}^{0.17} P_{\text{L}}^{-0.68} a^{0.79}
$$

!!! warning
    The origin, name and values of this scaling cannot be confirmed.

-------------------------

#### 30: Kaye 98 modified scaling (L-mode)

Is selected with `i_confinement_time = 30`

$$
\tau_{\text{E}} =  0.021 M_{\text{i}}^{0.25} I_{\text{p}}^{0.81} R^{2.01} \kappa^{0.7} \overline{n}_{19}^{0.47} B_{\text{T}}^{0.14} P_{\text{L}}^{-0.73} \epsilon^{0.18}
$$

!!! warning
    The origin, name and values of this scaling cannot be confirmed.

-------------------------

#### 31: ITERH-PB98P(y) scaling (H-mode)

Is selected with `i_confinement_time = 31` 

$$
\tau_{\text{E}} =  0.0615 M^{0.2} I_{\text{p}}^{0.9} R^{2.0} \kappa_{\text{IPB}}^{0.75} \overline{n}_{19}^{0.4} B_{\text{T}}^{0.1} P_{\text{L}}^{-0.66} \epsilon^{0.66}
$$

!!! warning
    The origin, name and values of this scaling cannot be confirmed.

-------------------------

#### 32: IPB98(y) ELMy scaling (H-mode)

Is selected with `i_confinement_time = 32` [^11] [^12]

$$
\tau_{\text{E}} =  0.0365 I_{\text{p}}^{0.97} B_{\text{T}}^{0.08} \overline{n}_{19}^{0.41} P_{\text{L}}^{-0.63} R^{1.93} \kappa^{0.67} \epsilon^{0.23} M^{0.2}
$$

-------------------------

#### 33: IPB98(y,1) ELMy scaling (H-mode)

Is selected with `i_confinement_time = 33` [^11] [^12]

$$
\tau_{\text{E}} =  0.0503 I_{\text{p}}^{0.91} B_{\text{T}}^{0.15} \overline{n}_{19}^{0.44} P_{\text{L}}^{-0.65} R^{2.05} \kappa_{\text{IPB}}^{0.72} \epsilon^{0.57} M^{0.13}
$$

-------------------------

#### 34: IPB98(y,2) ELMy scaling (H-mode)

Is selected with `i_confinement_time = 34` [^11] [^12]

$$
\tau_{\text{E}} =  0.0562 I_{\text{p}}^{0.93} B_{\text{T}}^{0.15} \overline{n}_{19}^{0.41} P_{\text{L}}^{-0.69} R^{1.97} \kappa_{\text{IPB}}^{0.78} \epsilon^{0.58} M^{0.19}
$$

-------------------------

#### 35: IPB98(y,3) ELMy scaling (H-mode)

Is selected with `i_confinement_time = 35` [^11] [^12]

$$
\tau_{\text{E}} =  0.0564 I_{\text{p}}^{0.88} B_{\text{T}}^{0.07} \overline{n}_{19}^{0.4} P_{\text{L}}^{-0.69} R^{2.15} \kappa_{\text{IPB}}^{0.78} \epsilon^{0.64} M^{0.2}
$$

-------------------------

#### 36: IPB98(y,4) ELMy scaling (H-mode)

Is selected with `i_confinement_time = 36` [^11] [^12]

$$
\tau_{\text{E}} =  0.0587 I_{\text{p}}^{0.85} B_{\text{T}}^{0.29} \overline{n}_{19}^{0.39} P_{\text{L}}^{-0.7} R^{2.08} \kappa_{\text{IPB}}^{0.76} \epsilon^{0.69} M^{0.17}
$$

-------------------------


#### 37: ISS95 scaling (Stellarator)

Is selected with `i_confinement_time = 37` [^13]

$$
\tau_{\text{E}} =  0.079 a^{2.21} R^{0.65} P_{\text{L}}^{-0.59} \overline{n}_{19}^{0.51} B_{\text{T}}^{0.83} \iota_{2/3}^{0.4}
$$

-------------------------


#### 38: ISS04 scaling (Stellarator)

Is selected with `i_confinement_time = 38` [^14]

$$
\tau_{\text{E}} =  0.134 a^{2.28} R^{0.64} P_{\text{L}}^{-0.61} \overline{n}_{19}^{0.54} B_{\text{T}}^{0.84} \iota_{2/3}^{0.41}
$$

-------------------------

#### 39: DS03 beta-independent scaling (H-mode)

Is selected with `i_confinement_time = 39` [^15]

$$
\tau_{\text{E}} =  0.028 I_{\text{p}}^{0.83} B_{\text{T}}^{0.07} \overline{n}_{19}^{0.49} P_{\text{L}}^{-0.55} R^{2.11} \kappa_{95}^{0.75} \epsilon^{0.3} M^{0.14}
$$

-------------------------

#### 40: Murari "Non-power law" scaling (H-mode)

Is selected with `i_confinement_time = 40` [^16]

$$
\tau_{\text{E}} =  0.0367 I_{\text{p}}^{1.006} R^{1.731} \kappa_{\text{IPB}}^{1.45} P_{\text{L}}^{-0.735} \\
\times \frac{\overline{n}_{19}^{0.49}}{1+e^\left({-9.403\left(\frac{\overline{n}_{19}^{0.49}}{B_{\text{T}}}\right)^{-1.365}}\right)}
$$

-------------------------

#### 41: Petty08 scaling (H-mode)

Is selected with `i_confinement_time = 41` [^17]

$$
\tau_{\text{E}} =  0.052 I_{\text{p}}^{0.75} B_{\text{T}}^{0.3} \overline{n}_{19}^{0.32} P_{\text{L}}^{-0.47} R^{2.09} \kappa_{\text{IPB}}^{0.88} \epsilon^{0.84}
$$

-------------------------

#### 42: Lang high density scaling (H-mode)

Is selected with `i_confinement_time = 42` [^18]

$$
\tau_{\text{E}} =  6.94\times 10^{-7} M^{0.2} \kappa_{\text{IPB}}^{0.37} \left(\frac{q_{95}}{q_{\text{cyl}}}\right)^{0.77} \\
\times A^{2.48205} \frac{I_{\text{p}}^{1.3678} B_{\text{T}}^{0.12} R^{1.2345} \overline{n}^{0.032236}}{A^{0.9\ln{A}}P_{\text{L}}^{0.74}} \left(\frac{\overline{n}_{e}}{n_{\text{GW}}}\right)^{-0.22 \ln{\left(\frac{\overline{n}_e}{n_{\text{GW}}}\right)}}
$$

-------------------------

#### 43: Hubbard nominal scaling (I-mode)

Is selected with `i_confinement_time = 43` [^19]

$$
\tau_{\text{E}} =  0.014 I_{\text{p}}^{0.68} B_{\text{T}}^{0.77} \overline{n}_{20}^{0.02} P_{\text{L}}^{-0.29}
$$

-------------------------

#### 44: Hubbard lower scaling (I-mode)

Is selected with `i_confinement_time = 44` [^19]

$$
\tau_{\text{E}} =  0.014 I_{\text{p}}^{0.6} B_{\text{T}}^{0.7} \overline{n}_{20}^{-0.03} P_{\text{L}}^{-0.33}
$$

-------------------------

#### 45: Hubbard upper scaling (I-mode)

Is selected with `i_confinement_time = 45` [^19]

$$
\tau_{\text{E}} =  0.014 I_{\text{p}}^{0.76} B_{\text{T}}^{0.84} \overline{n}_{20}^{-0.07} P_{\text{L}}^{-0.25}
$$

-------------------------


#### 46: Menard NSTX scaling (H-mode)

Is selected with `i_confinement_time = 46` [^20]

$$
\tau_{\text{E}} =  0.095 I_{\text{p}}^{0.75} B_{\text{T}}^{1.08} \overline{n}_{19}^{0.44} P_{\text{L}}^{-0.73} R^{1.97} \kappa_{\text{IPB}}^{0.78} \epsilon^{0.58} M^{0.19}
$$

-------------------------

#### 47: Menard NSTX-Petty08 hybrid scaling

Is selected with `i_confinement_time = 47` [^20]

- If $\epsilon \le 0.4 \  (A \ge 2.5)$ apply the [Petty08 scaling](#41-petty-h-mode-scaling)
- If $\epsilon \ge 0.6 \ (A \le 1.7)$ apply the [Menard NSTX scaling](#46-menard-nstx-h-mode-scaling)

Otherwise:

$$
\tau_{\text{E}} =  \frac{\epsilon - 0.4}{0.2}\tau_{\text{E,NSTX}}+ \frac{0.6-\epsilon}{0.2}\tau_{\text{E,Petty08}}
$$

-------------------------

#### 48: Buxton NSTX Gyro-Bohm scaling (H-mode)

Is selected with `i_confinement_time = 48` [^21]

$$
\tau_{\text{E}} =  0.21 I_{\text{p}}^{0.54} B_{\text{T}}^{0.91} \overline{n}_{20}^{-0.05} P_{\text{L}}^{-0.38} R^{2.14}
$$

-------------------------

#### 49: ITPA20 scaling (H-mode)

Is selected with `i_confinement_time = 49` [^22]

$$
\tau_{\text{E}} =  0.053 I_{\text{p}}^{0.98} B_{\text{T}}^{0.22} \overline{n}_{19}^{0.24} P_{\text{L}}^{-0.669} R^{1.71} \left(1+\delta \right)^{0.36}  \kappa_{\text{IPB}}^{0.8} \epsilon^{0.35} M^{0.2}
$$

-------------------------

#### 50: ITPA20-IL scaling (H-mode)

Is selected with `i_confinement_time = 50` [^23]

$$
\tau_{\text{E}} =  0.067 I_{\text{p}}^{1.29} B_{\text{T}}^{-0.13} P_{\text{L}}^{-0.644} \overline{n}_{19}^{0.15} M^{0.3} R^{1.19} \left(1+\delta \right)^{0.56} \kappa_{\text{IPB}}^{0.67} 
$$

-------------------------

### Transport Powers

After the confinement time scaling with $H$-factor correction has been calculated, the ion and electron transport power densities are found. `PROCESS` assumes the scaling confinement time to be equal to the ion and electron energy confinement time.

This is simply the volume averaged thermal energy of the electron and ions divided by the $H$-factor corrected confinement time from the chosen scaling.

$$
\mathtt{pden\_ion\_transport\_loss\_mw} = \frac{3}{2}\frac{n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}}{\tau_{\text{E}}}
$$

$$
\mathtt{pden\_electron\_transport\_loss\_mw} = \frac{3}{2}\frac{n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{n}}}{\tau_{\text{E}}}
$$

Here $\langle T_{\text{i}} \rangle$ and $\langle T_{\text{e}} \rangle$ are the ion and electron density weighted temperatures respectively.

Calculate the density and density weighted ratio:

$$
\frac{n_{\text{i}}}{n_{\text{e}}}\frac{\langle T_{\text{i}} \rangle_{\text{n}}}{\langle T_{\text{e}} \rangle_{\text{n}}}
$$

The density weighted global energy confinement time is then found in terms of this ratio:

$$
\tau_{\text{E}} = \frac{\frac{n_{\text{i}}}{n_{\text{e}}}\frac{\langle T_{\text{i}} \rangle_{\text{n}}}{\langle T_{\text{e}} \rangle_{\text{n}}} + 1}{\left(\frac{\frac{n_{\text{i}}}{n_{\text{e}}}\frac{\langle T_{\text{i}} \rangle_{\text{n}}}{\langle T_{\text{e}} \rangle_{\text{n}}}}{\tau_{\text{i}}}+\frac{1}{\tau_{\text{e}}}\right)}
$$

----------

## Key Constraints

### Global plasma power balance

This constraint can be activated by stating `icc = 2` in the input file.

This constraint ensures self consistency between the the transport loss power used for the confinement scalings and the calculated confinement time in relation to the plasmas total thermal energy:

$$
P_{\text{L}} = \frac{W}{\tau_{\text{E}}}
$$

$$
\underbrace{\frac{3}{2}\frac{n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}}{\tau_{\text{E}}} + \frac{3}{2}\frac{n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{n}}}{\tau_{\text{E}}}}_{\frac{W}{\tau_{\text{E}}}}  = \underbrace{\frac{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}{V_{\text{P}}} - \frac{P_{\text{rad}}}{V_{\text{p}}}}_{P_{\text{L}}}
$$

The $\frac{3}{2}n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}$ value is simply the volume averaged ion thermal energy density where $\langle T_{\text{i}} \rangle_{\text{n}}$ is the density weighted temperature. The same goes for the $\frac{3}{2}n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{e}}$ electron thermal energy density term. $\tau_{\text{E}}$ is the confinement time calculated from the chosen confinement scaling via `i_confinement_time`. 

The constraint uses the loss power and thermal densities hence the inclusion of the $V_{\text{p}}$ plasma volume term. The constraint is adapted depending on the condition of `i_rad_loss` which governs the radiation contribution to the loss power definition, see the [radiation and energy confinement section](#effect-of-radiation-on-energy-confinement) for more info. The injected heating and current drive contribution $P_{\text{HCD}}$ is also included or excluded depending if the plasma is deemed to be ignited with the `i_plasma_ignited` switch.

**It is highly recommended to always have this constraint on as it is a global consistency checker**

----------

### Lower limit on alpha particle confinement time ratio

This constraint can be activated by stating `icc = 62` in the input file.

The value of `f_alpha_energy_confinement_min` can be set to the desired minimum total ratio between the alpha confinement and energy confinement times.

[^1]: N. A. Uckan, International Atomic Energy Agency, Vienna (Austria) and ITER Physics Group, "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
[^2]: T.C. Hender et al., 'Physics Assessment of the European Reactor Study', AEA FUS 172, 1992.
[^3]: J. P. Christiansen et al., “Global energy confinement H-mode database for ITER,” Nuclear Fusion, vol. 32, no. 2, pp. 291-338, Feb. 1992, doi: https://doi.org/10.1088/0029-5515/32/2/i11.
[^4]: S. Sudo et al., “Scalings of energy confinement and density limit in stellarator/heliotron devices,” Nuclear Fusion, vol. 30, no. 1, pp. 11-21, Jan. 1990, doi: https://doi.org/10.1088/0029-5515/30/1/002.
[^5]: Goldston, R. J., H. Biglari, and G. W. Hammett. "E x B/B² vs. µ B/B as the Cause of Transport in Tokamaks." Bull. Am. Phys. Soc 34 (1989): 1964.
[^6]: K. Lackner and N. A. O. Gottardi, “Tokamak confinement in relation to plateau scaling,” Nuclear Fusion, vol. 30, no. 4, pp. 767-770, Apr. 1990, doi: https://doi.org/10.1088/0029-5515/30/4/018.
[^7]: K. Thomsen et al., “ITER H mode confinement database update,” vol. 34, no. 1, pp. 131-167, Jan. 1994, doi: https://doi.org/10.1088/0029-5515/34/1/i10.
[^8]: I. C. Database and M. W. G. (presented Cordey), “Energy confinement scaling and the extrapolation to ITER,” Plasma Physics and Controlled Fusion, vol. 39, no. 12B, pp. B115-B127, Dec. 1997, doi: https://doi.org/10.1088/0741-3335/39/12b/009.
[^9]: International Atomic Energy Agency, Vienna (Austria), "Technical basis for the ITER final design report, cost review and safety analysis (FDR)", no. 16. Dec. 1998.
[^10]: S. B. Kaye et al., “ITER L mode confinement database,” Nuclear Fusion, vol. 37, no. 9, pp. 1303-1328, Sep. 1997, doi: https://doi.org/10.1088/0029-5515/37/9/i10.
[^11]: I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,” Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.
[^12]: O. Kardaun, N. K. Thomsen, and A. Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,” Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
[^13]: U. Stroth et al., “Energy confinement scaling from the international stellarator database,” vol. 36, no. 8, pp. 1063-1077, Aug. 1996, doi: https://doi.org/10.1088/0029-5515/36/8/i11.
[^14]: H. Yamada et al., “Characterization of energy confinement in net-current free plasmas using the extended International Stellarator Database,” vol. 45, no. 12, pp. 1684-1693, Nov. 2005, doi: https://doi.org/10.1088/0029-5515/45/12/024.
[^15]: T. C. Luce, C. C. Petty, and J. G. Cordey, “Application of dimensionless parameter scaling techniques to the design and interpretation of magnetic fusion experiments,” Plasma Physics and Controlled Fusion, vol. 50, no. 4, p. 043001, Mar. 2008, doi: https://doi.org/10.1088/0741-3335/50/4/043001.
[^16]: A. Murari, E. Peluso, M. Gelfusa, I. Lupelli, and P. Gaudio, “A new approach to the formulation and validation of scaling expressions for plasma confinement in tokamaks,” Nuclear Fusion, vol. 55, no. 7, pp. 073009-073009, Jun. 2015, doi: https://doi.org/10.1088/0029-5515/55/7/073009.
[^17]: C. C. Petty, “Sizing up plasmas using dimensionless parameters,” Physics of Plasmas, vol. 15, no. 8, Aug. 2008, doi: https://doi.org/10.1063/1.2961043.
[^18]: P. T. Lang, C. Angioni, R. M. M. Dermott, R. Fischer, and H. Zohm, “Pellet Induced High Density Phases during ELM Suppression in ASDEX Upgrade,” 24th IAEA Conference Fusion Energy, 2012, Oct. 2012, Available: https://www.researchgate.net/publication/274456104_Pellet_Induced_High_Density_Phases_during_ELM_Suppression_in_ASDEX_Upgrade.
[^19]: A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,” Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.
[^20]: J. E. Menard, “Compact steady-state tokamak performance dependence on magnet and core physics limits,” Philosophical Transactions of the Royal Society A, vol. 377, no. 2141, pp. 20170440-20170440, Feb. 2019, doi: https://doi.org/10.1098/rsta.2017.0440.
[^21]: P. F. Buxton, L. Connor, A. E. Costley, M. Gryaznevich, and S. McNamara, “On the energy confinement time in spherical tokamaks: implications for the design of pilot plants and fusion reactors,” vol. 61, no. 3, pp. 035006-035006, Jan. 2019, doi: https://doi.org/10.1088/1361-6587/aaf7e5.
[^22]: G. Verdoolaege et al., “The updated ITPA global H-mode confinement database: description and analysis,” Nuclear Fusion, vol. 61, no. 7, pp. 076006-076006, Jan. 2021, doi: https://doi.org/10.1088/1741-4326/abdb91.
[^23]: T. Luda et al., “Validation of a full-plasma integrated modeling approach on ASDEX Upgrade,” Nuclear Fusion, vol. 61, no. 12, pp. 126048-126048, Nov. 2021, doi: https://doi.org/10.1088/1741-4326/ac3293.
[^24]: H. Lux, R. Kemp, E. Fable, and R. Wenninger, “Radiation and confinement in 0D fusion systems codes,” Plasma Physics and Controlled Fusion, vol. 58, no. 7, pp. 075001–075001, May 2016, doi: https://doi.org/10.1088/0741-3335/58/7/075001.
[^25]: H. Lux, R. Kemp, D. J. Ward, and M. Sertoli, “Impurity radiation in DEMO systems modelling,” Fusion Engineering and Design, vol. 101, pp. 42–51, Dec. 2015, doi: https://doi.org/10.1016/j.fusengdes.2015.10.002.
‌