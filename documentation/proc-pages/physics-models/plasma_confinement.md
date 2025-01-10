# Confinement Time Scaling Laws

The energy confinement time $\tau_E$ is calculated using one of a choice of empirical scalings. ($\tau_E$ is defined below.)

Many energy confinement time scaling laws are available within PROCESS, for
tokamaks, RFPs and stellarators. These are calculated in routine `pcond`. The 
value of `isc` determines which of the scalings is used in the plasma energy 
balance calculation. The table below summarises the available scaling laws. The 
most commonly used is the so-called IPB98(y,2) scaling.       

| `isc` | scaling law | reference |
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

## L-H Power Threshold Scalings

Transitions from a standard confinement mode (L-mode) to an improved
confinement regime (H-mode), called L-H transitions, are observed in most
tokamaks. A range of scaling laws are available that provide estimates of the
heating power required to initiate these transitions, via extrapolations
from present-day devices. PROCESS calculates these power threshold values
for the scaling laws listed in the table below, in routine `pthresh`.

For an H-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with
iteration variable no. 103 (`flhthresh`). By default, this will ensure
that the power reaching the divertor is at least equal to the threshold power
calculated for the chosen scaling, which is a necessary condition for
H-mode. 

For an L-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with 
iteration variable no. 103 (`flhthresh`). Set lower and upper bounds for 
the f-value `boundl(103) = 0.001` and `boundu(103) = 1.0` 
to ensure that the power does not exceed the calculated threshold, 
and therefore the machine remains in L-mode.


| `ilhthresh` | Name | Reference |
| :-: | - | - |
| 1 | ITER 1996 nominal | ITER Physics Design Description Document |
| 2 | ITER 1996 upper bound | D. Boucher, p.2-2 |
| 3 | ITER 1996 lower bound | 
| 4 | ITER 1997 excluding elongation | J. A. Snipes, ITER H-mode Threshold Database |
| 5 | ITER 1997 including elongation |  Working Group, Controlled Fusion and Plasma Physics, 24th EPS conference, Berchtesgaden, June 1997, vol.21A, part III, p.961 |
| 6 | Martin 2008 nominal | Martin et al, 11th IAEA Tech. Meeting |
| 7 | Martin 2008 95% upper bound |  H-mode Physics and Transport Barriers, Journal |
| 8 | Martin 2008 95% lower bound |  of Physics: Conference Series **123**, 2008 |
| 9 | Snipes 2000 nominal | J. A. Snipes and the International H-mode |
| 10| Snipes 2000 upper bound | Threshold Database Working Group |
| 11| Snipes 2000 lower bound |  2000, Plasma Phys. Control. Fusion, 42, A299 |
| 12| Snipes 2000 (closed divertor): nominal | 
| 13| Snipes 2000 (closed divertor): upper bound | 
| 14| Snipes 2000 (closed divertor): lower bound | 
| 15| Hubbard 2012 L-I threshold scaling: nominal | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009) |
| 16| Hubbard 2012 L-I threshold scaling: lower bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 17| Hubbard 2012 L-I threshold scaling: upper bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 18| Hubbard 2017 L-I threshold scaling | [Hubbard et al. (2017; Nucl. Fusion 57 126039)](https://iopscience.iop.org/article/10.1088/1741-4326/aa8570) |
| 19 | Martin 2008 aspect ratio corrected nominal | Martin et al (2008; J Phys Conf, 123, 012033) |
| 20 | Martin 2008 aspect ratio corrected 95% upper bound | [Takizuka et al. (2004; Plasma Phys. Contol. Fusion, 46, A227)](https://iopscience.iop.org/article/10.1088/0741-3335/46/5A/024)  |
| 21 | Martin 2008 aspect ratio corrected 95% lower bound |  

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