# Confinement Time Scaling Laws

The energy confinement time $\tau_E$ is calculated using one of a choice of empirical scalings. ($\tau_E$ is defined below.)

Many energy confinement time scaling laws are available within PROCESS, for
tokamaks, RFPs and stellarators. These are calculated in routine `pcond`. The 
value of `i_confinement_time` determines which of the scalings is used in the plasma energy 
balance calculation. The table below summarises the available scaling laws. The 
most commonly used is the so-called IPB98(y,2) scaling.

## Effect of radiation on energy confinement

Published confinement scalings are all based on low radiation pulses. A power
plant will certainly be a high radiation machine --- both in the core, due to
bremsstrahlung and synchrotron radiation, and in the edge due to impurity
seeding. The scaling data do not predict this radiation --- that needs to be
done by the radiation model. However, if the transport is very "stiff", as
predicted by some models, then the additional radiation causes an almost equal
drop in power transported by ions and electrons, leaving the confinement
nearly unchanged.

To allow for these uncertainties, three options are available, using the switch
`i_rad_loss`. In each case, the particle transport loss power `pscaling` is
derived directly from the energy confinement scaling law.

`i_rad_loss = 0` -- Total power lost is scaling power plus radiation:

`pscaling + pradpv = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/plasma_volume`


`i_rad_loss = 1` -- Total power lost is scaling power plus radiation from a region defined as the "core":
  
`pscaling + pcoreradpv = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/plasma_volume`

`i_rad_loss = 2` -- Total power lost is scaling power only, with no additional 
allowance for radiation. This is not recommended for power plant models.

`pscaling = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/plasma_volume`

## Available confinement time scalings

### 0: User input confinement time

Is selected with `i_confinement_time = 0`

$$
\tau_{\text{E}} = \mathtt{tauee\_in}
$$

------------

### 1: Nec-Alcator (NA) OH scaling

Is selected with `i_confinement_time = 1`

$$
\tau_{\text{E}} = 0.07 n_{20}aRq_{\text{cyl}}
$$

------------

### 2: Mirnov scaling (H-mode)

Is selected with `i_confinement_time = 2`

$$
\tau_{\text{E}} = 0.2 a \sqrt{\kappa_{95}}I
$$

------------

### 3: Merezhkin-Mukhovatov  OH/L-mode scaling

Is selected with `i_confinement_time = 3`

$$
\tau_{\text{E}} = 0.0035 \overline{n}_{20}a^{0.25}R^{2.75}q_{\text{cyl}}\kappa_{95}^{0.125}M_i^{0.5}T_{10}^{0.5}
$$


---------------

### 4: Shimomura optimized H-mode scaling

Is selected with `i_confinement_time = 4`

$$
\tau_{\text{E}} = 0.045 Ra B_{\text{T}}\sqrt{\kappa_{95}}\sqrt{M_{\text{i}}}
$$

----------------

### 5: Kaye-Goldston L-mode scaling

Is selected with `i_confinement_time = 5`

$$
\tau_{\text{E}} = 0.055 I^{1.24}P^{-0.58}R^{1.65}a^{-0.49}\kappa_{95}^{0.28}n_{20}^{0.26}B_{\text{T}}^{-0.09}\left(\frac{M_{\text{i}}}{1.5}\right)^{0.5}
$$

----------------

### 6: ITER 89-P L-mode scaling

Is selected with `i_confinement_time = 6`

$$
\tau_{\text{E}} = 0.048 I^{0.85}R^{1.2}a^{0.3}\kappa^{0.5}\overline{n}_{20}^{0.1}B_{\text{T}}^{0.2}M_{\text{i}}^{0.5} P^{-0.5}
$$

----------------

### 7: ITER 89-0 L-mode scaling

Is selected with `i_confinement_time = 7`

$$
\begin{aligned}
\tau_E= & 0.04 I^{0.5} R^{0.3} a^{0.8} \kappa^{0.6} M_i^{0.5} \\
& +0.064 I^{0.8} R^{1.6} a^{0.6} \kappa^{0.2} \bar{n}_{20}^{0.6} B_0^{0.35} M_i^{0.2} / P
\end{aligned}
$$

----------------

### 8: Rebut-Lallia L-mode scaling

Is selected with `i_confinement_time = 8`

$$
\begin{aligned}
\tau_E= & 1.65\left[1.2 \times 10^{-5} I \ell^{1.5} Z_{e f f}^{-0.5}\right. \\
& \left.+0.146 \bar{n}_{20}^{0.75} I^{0.5} B_0^{0.5} \ell^{2.75} Z_{e f f}^{0.25} / P\right]\left(A_i / 2\right)^{0.5}
\end{aligned}
$$

where $\ell = \left(a^2R\kappa\right)^{\frac{1}{3}}$

----------------


### 9: Goldston L-mode scaling

Is selected with `i_confinement_time = 9`

$$
\tau_{\text{E}} = 0.037 I P^{-0.5} R^{1.75}a^{-0.37}\kappa_{95}^{0.5} \left(\frac{M_i}{1.5}\right)^{0.5}
$$

----------------

### 10: T-10 L-mode scaling

Is selected with `i_confinement_time = 10`


$$
\tau_{\text{E}} =  0.095 a R B_{\text{T}} \kappa_{95}^{0.5} \frac{\overline{n}_{20}}{\overline{n}_{20*}}P^{-0.4} \left[\frac{Z_{\text{eff}}^2 I^4}{aRq_{\text{cyl}}^3\kappa_{95}^{1.5}}  \right]^{0.08}
$$

where $\overline{n}_{20*} = 1.3\left(\frac{B_{\text{T}}}{Rq_{\text{cyl}}}\right)$ and $\frac{\overline{n}_{20}}{\overline{n}_{20*}} \le 1$

----------------

### 11: JAERI / Odajima-Shimomura L-mode scaling

Is selected with `i_confinement_time = 11`


$$
\tau_{\text{E}} =  \left[\frac{0.085\kappa_{95}a^2+0.069In_{20}^{0.6}B_{\text{T}}^{0.2}R^{1.6} a^{0.4} \kappa_{95}^{0.2} G\left(q_{\text{cyl}},Z_{\text{eff}}\right)}{P}\right]M_{\text{i}}^{0.5}
$$

where $G\left(q_{\text{cyl}},Z_{\text{eff}}\right) = Z_{\text{eff}}^{0.4}\left[\frac{\left(15 - Z_{\text{eff}}\right)}{20}\right]^{0.6}\left[3q_{\text{cyl}}\frac{q_{\text{cyl}}+5}{(q_{\text{cyl}}+2)(q_{\text{cyl}}+7)}\right]^{0.6}$


----------------

### 12: Kaye "big" L-mode scaling

Is selected with `i_confinement_time = 12`

$$
\tau_{\text{E}} =  0.1051 I^{0.85} P^{-0.5} R^{0.5} a^{0.3} \kappa^{0.25} n_{20}^{0.1}B_{\text{T}}^{0.3}M_{\text{i}}^{0.5}
$$

-------------------------

### 13: ITER H90-P H-mode scaling

Is selected with `i_confinement_time = 13`

$$
\tau_{\text{E}} =  0.064 I^{0.87} R^{1.82} a^{-0.12} \kappa_{95}^{0.35} \overline{n}_{20}^{0.09} B_{\text{T}}^{0.15} M_{\text{i}}^{0.5} P^{-0.5}
$$

-------------------------

### 14: Minimum of ITER 89-P and ITER 89-O

Is selected with `i_confinement_time = 14`

Will return the value of [ITER 89-P](#6-iter-89-p-l-mode-scaling) or [ITER 89-O](#7-iter-89-0-l-mode-scaling), whichever is smaller.

-------------------------

### 15: Riedel L-mode scaling

Is selected with `i_confinement_time = 15`

$$
\tau_{\text{E}} =  0.044 I^{0.93} R^{1.37} a^{-0.049} \kappa_{95}^{0.588} \overline{n}_{20}^{0.078} B_{\text{T}}^{0.152} P^{-0.537}
$$

-------------------------

### 16: Christiansen L-mode scaling

Is selected with `i_confinement_time = 16`

$$
\tau_{\text{E}} =  0.24 I^{0.79} R^{0.56} a^{1.46} \kappa_{95}^{0.73} \overline{n}_{20}^{0.41} B_{\text{T}}^{0.29} P^{-0.79} M_{\text{i}}^{-0.02}
$$

-------------------------

### 17: Lackner-Gottardi L-mode scaling

Is selected with `i_confinement_time = 17`

$$
\tau_{\text{E}} =  0.12 I^{0.8} R^{1.8} a^{0.4} \left(\frac{\kappa_{95}}{\left(1+\kappa_{95}\right)^{0.8}}\right) \overline{n}_{20}^{0.6} \hat{q}^{0.4} P^{-0.6}
$$

where $\hat{q} = \frac{(1+\kappa_{95}a^2B_{\text{T}})}{0.4 I R}$

-------------------------

### 18: Neo-Kaye L-mode scaling

Is selected with `i_confinement_time = 18`

$$
\tau_{\text{E}} =  0.063 I^{1.12} R^{1.3} a^{-0.04} \kappa_{95}^{0.28} \overline{n}_{20}^{0.14} B_{\text{T}}^{0.04} P^{-0.59}
$$

-------------------------

### 19: Riedel H-mode scaling

Is selected with `i_confinement_time = 19`

$$
\tau_{\text{E}} =  0.1 M_{\text{i}}^{0.5} I^{0.884} R^{1.24} a^{-0.23} \kappa_{95}^{0.317} \overline{n}_{20}^{0.105} B_{\text{T}}^{0.207} P^{-0.486}
$$

-------------------------

### 20: Amended ITER H90-P H-mode scaling 

Is selected with `i_confinement_time = 20`

$$
\tau_{\text{E}} =  0.082 M_{\text{i}}^{0.5} I^{1.02} R^{1.6}  \kappa_{95}^{-0.19}  B_{\text{T}}^{0.15} P^{-0.47}
$$

-------------------------

### 21: Sudo et al. stellarators/heliotron scaling

Is selected with `i_confinement_time = 21`

$$
\tau_{\text{E}} =  0.17 P^{-0.58} \overline{n}_{20}^{0.69} B^{0.84} a^{2.0} R^{0.75}
$$

-------------------------

### 22: Gyro reduced Bohm scaling

Is selected with `i_confinement_time = 22`

$$
\tau_{\text{E}} =  0.25 P^{-0.6} \overline{n}_{20}^{0.6} B_{\text{T}}^{0.8} a^{2.4} R^{0.6}
$$

-------------------------

### 23: Lackner-Gottardi Stellerator scaling

Is selected with `i_confinement_time = 23`

$$
\tau_{\text{E}} =  0.17 P^{-0.6} \overline{n}_{20}^{0.6} B_{\text{T}}^{0.8} a^{2.0} R q_{95}^{0.4}
$$

-------------------------

### 26: ELM-free: ITERH-97P scaling

Is selected with `i_confinement_time = 26`

$$
\tau_{\text{E}} =  0.031 M_{\text{i}}^{0.42} I^{0.95} R^{1.92} \epsilon^{0.08} \kappa_{95}^{0.63} \overline{n}_{19}^{0.35} B_{\text{T}}^{0.25} P^{-0.67}
$$

-------------------------

### 27: ELMy ITER H-H97-P(y) scaling

Is selected with `i_confinement_time = 27`

$$
\tau_{\text{E}} =  0.029 M_{\text{i}}^{0.2} I^{0.9} R^{2.03} \epsilon^{-0.19} \kappa_{95}^{0.92} \overline{n}_{19}^{0.4} B_{\text{T}}^{0.20} P^{-0.66}
$$

-------------------------

### 28: ITER-96P (ITER-97L) L-mode scaling

Is selected with `i_confinement_time = 28`

$$
\tau_{\text{E}} =  0.023 M_{\text{i}}^{0.2} I^{0.96} R^{1.83} \epsilon^{-0.06} \kappa_{95}^{0.64} \overline{n}_{19}^{0.4} B_{\text{T}}^{0.03} P^{-0.73}
$$

-------------------------

### 32: IPB98(y) ELMy H-mode scaling

Is selected with `i_confinement_time = 32`

$$
\tau_{\text{E}} =  0.0365 I^{0.97} B_{\text{T}}^{0.08} \overline{n}_{19}^{0.41} P^{-0.63} R^{1.93} \kappa^{0.67} \epsilon^{0.23} M^{0.2}
$$

-------------------------


### 33: IPB98(y,1) ELMy H-mode scaling

Is selected with `i_confinement_time = 33`

$$
\tau_{\text{E}} =  0.0503 I^{0.91} B_{\text{T}}^{0.15} \overline{n}_{19}^{0.44} P^{-0.65} R^{2.05} \kappa_{\text{IPB}}^{0.72} \epsilon^{0.57} M^{0.13}
$$

-------------------------

### 34: IPB98(y,2) ELMy H-mode scaling

Is selected with `i_confinement_time = 34`

$$
\tau_{\text{E}} =  0.0562 I^{0.93} B_{\text{T}}^{0.15} \overline{n}_{19}^{0.41} P^{-0.69} R^{1.97} \kappa_{\text{IPB}}^{0.78} \epsilon^{0.58} M^{0.19}
$$

-------------------------

### 35: IPB98(y,3) ELMy H-mode scaling

Is selected with `i_confinement_time = 35`

$$
\tau_{\text{E}} =  0.0564 I^{0.88} B_{\text{T}}^{0.07} \overline{n}_{19}^{0.4} P^{-0.69} R^{2.15} \kappa_{\text{IPB}}^{0.78} \epsilon^{0.64} M^{0.2}
$$

-------------------------

### 36: IPB98(y,4) ELMy H-mode scaling

Is selected with `i_confinement_time = 36`

$$
\tau_{\text{E}} =  0.0587 I^{0.85} B_{\text{T}}^{0.29} \overline{n}_{19}^{0.39} P^{-0.7} R^{2.08} \kappa_{\text{IPB}}^{0.76} \epsilon^{0.69} M^{0.17}
$$

-------------------------

| `i_confinement_time` | scaling law | reference |
| :-: | - | - |
| 1 | Neo-Alcator (ohmic) | [^1] | 
| 2 | Mirnov (H-mode) | [^1] |
| 3 | Merezhkin-Muhkovatov (L-mode) | [^1] |
| 4 | Shimomura (H-mode) | JAERI-M 87-080 (1987) |
| 5 | Kaye-Goldston (L-mode) | Nuclear Fusion **25** (1985) p.65 |
| 6 | ITER 89-P (L-mode) | Nuclear Fusion **30** (1990) p.1999 |
| 7 | ITER 89-O (L-mode) | [^2] |
| 8 | Rebut-Lallia (L-mode) | Plasma Physics and Controlled Nuclear Fusion Research **2** (1987) p. 187 |
| 9  | Goldston (L-mode)| Plas.\ Phys.\ Controlled Fusion **26** (1984) p.87 |
| 10 | T10 (L-mode) | [^2] |
| 11 | JAERI-88 (L-mode) | JAERI-M 88-068 (1988) |
| 12 | Kaye-Big Complex (L-mode) | Phys.\ Fluids B **2** (1990) p.2926 |
| 13 | ITER H90-P (H-mode) |  |
| 14 | ITER Mix (minimum of 6 and 7) |  |
| 15 | Riedel (L-mode) |  |
| 16 | Christiansen et al. (L-mode) | JET Report JET-P (1991) 03 |
| 17 | Lackner-Gottardi (L-mode) | Nuclear Fusion **30** (1990) p.767  |
| 18 | Neo-Kaye (L-mode) | [^2] |
| 19 | Riedel (H-mode) |  |
| 20 | ITER H90-P (amended) | Nuclear Fusion **32** (1992) p.318 |
| 21 | Large Helical Device (stellarator) | Nuclear Fusion **30** (1990) |
| 22 | Gyro-reduced Bohm (stellarator) | Bull. Am. Phys. Society, **34** (1989) p.1964 |
| 23 | Lackner-Gottardi (stellarator) | Nuclear Fusion **30** (1990) p.767 |
| 24 | ITER-93H (H-mode) | PPCF, Proc. 15th Int. Conf.Seville, 1994 IAEA-CN-60/E-P-3 |
| 25 | TITAN (RFP) | TITAN RFP Fusion Reactor Study, Scoping Phase Report, UCLA-PPG-1100, page 5--9, Jan 1987 |
| 26 | ITER H-97P ELM-free (H-mode) | J. G. Cordey et al., EPS Berchtesgaden, 1997 |
| 27 | ITER H-97P ELMy (H-mode) | J. G. Cordey et al., EPS Berchtesgaden, 1997 |
| 28 | ITER-96P (= ITER97-L) (L-mode) | Nuclear Fusion **37** (1997) p.1303 |
| 29 | Valovic modified ELMy (H-mode) | |
| 30 | Kaye PPPL April 98 (L-mode) |  |
| 31 | ITERH-PB98P(y) (H-mode) | |
| 32 | IPB98(y) (H-mode)   | Nuclear Fusion **39** (1999) p.2175, Table 5,  |
| 33 | IPB98(y,1) (H-mode) | Nuclear Fusion **39** (1999) p.2175, Table 5, full data |
| 34 | IPB98(y,2) (H-mode) | Nuclear Fusion **39** (1999) p.2175, Table 5, NBI only |
| 35 | IPB98(y,3) (H-mode) | Nuclear Fusion **39** (1999) p.2175, Table 5, NBI only, no C-Mod |
| 36 | IPB98(y,4) (H-mode) | Nuclear Fusion **39** (1999) p.2175, Table 5, NBI only ITER like |
| 37 | ISS95 (stellarator) | Nuclear Fusion **36** (1996) p.1063  |
| 38 | ISS04 (stellarator) | Nuclear Fusion **45** (2005) p.1684  |
| 39 | DS03 (H-mode) | Plasma Phys. Control. Fusion **50** (2008) 043001, equation 4.13  |
| 40 | Non-power law (H-mode) | A. Murari et al 2015 Nucl. Fusion 55 073009, Table 4.  |
| 41 | Petty 2008 (H-mode) |  C.C. Petty 2008 Phys. Plasmas **15** 080501, equation 36 |
| 42 | Lang 2012 (H-mode) | P.T. Lang et al. 2012 IAEA conference proceeding EX/P4-01 |
| 43 | Hubbard 2017 -- nominal (I-mode) | A.E. Hubbard et al. 2017, Nuclear Fusion **57** 126039 |
| 44 | Hubbard 2017 -- lower (I-mode) | A.E. Hubbard et al. 2017, Nuclear Fusion **57** 126039 |
| 45 | Hubbard 2017 -- upper (I-mode) | A.E. Hubbard et al. 2017, Nuclear Fusion **57** 126039 |
| 46 | NSTX (H-mode; spherical tokamak) | J. Menard 2019, Phil. Trans. R. Soc. A 377:201704401 |
| 47 | NSTX-Petty08 Hybrid (H-mode) | J. Menard 2019, Phil. Trans. R. Soc. A 377:201704401 |
| 48 | NSTX gyro-Bohm (Buxton) (H-mode; spherical tokamak) | P. Buxton et al. 2019 Plasma Phys. Control. Fusion 61 035006 |
| 49 | Use input `tauee_in` |  |
| 50 | ITPA20 (H-mode) | G. Verdoolaege et al 2021 Nucl. Fusion 61 076006 |

### Effect of radiation on energy confinement

Published confinement scalings are all based on low radiation pulses. A power
plant will certainly be a high radiation machine --- both in the core, due to
bremsstrahlung and synchrotron radiation, and in the edge due to impurity
seeding. The scaling data do not predict this radiation --- that needs to be
done by the radiation model. However, if the transport is very "stiff", as
predicted by some models, then the additional radiation causes an almost equal
drop in power transported by ions and electrons, leaving the confinement
nearly unchanged.

To allow for these uncertainties, three options are available, using the switch
`iradloss`. In each case, the particle transport loss power `pscaling` is
derived directly from the energy confinement scaling law.

`iradloss = 0` -- Total power lost is scaling power plus radiation:

`pscaling + pden_plasma_rad_mw = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/vol_plasma`


`iradloss = 1` -- Total power lost is scaling power plus radiation from a region defined as the "core":
  
`pscaling + pcoreradpv = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/vol_plasma`

`iradloss = 2` -- Total power lost is scaling power only, with no additional 
allowance for radiation. This is not recommended for power plant models.

`pscaling = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + pinjmw/vol_plasma`


## Ignition

Switch `ignite` can be used to denote whether the plasma is ignited, i.e. fully self-sustaining 
without the need for any injected auxiliary power during the burn. If `ignite` = 1, the calculated 
injected power does not contribute to the plasma power balance, although the cost of the auxiliary 
power system is taken into account (the system is then assumed to be required to provide heating 
and/or current drive during the plasma start-up phase only). If `ignite` = 0, the plasma is not 
ignited, and the auxiliary power is taken into account in the plasma power balance during the burn 
phase. An ignited plasma will be difficult to control and is unlikely to be practical. This 
option is not recommended.

[^1]: T. C. Hender et al., 'Physics Assessment for the European Reactor Study',
AEA Fusion Report AEA FUS 172 (1992)
[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)