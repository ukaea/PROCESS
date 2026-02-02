"""Module containing tokamak plasma physics routines
author: P J Knight, CCFE, Culham Science Centre
N/A
This module contains all the primary plasma physics routines
for a tokamak device.

author: J. Morris (UKAEA)
Module containing global variables relating to the plasma physics
"""

import numpy as np

# From physics.f90:
iscz: int = None


err242: int = None


err243: int = None


f_p_plasma_separatrix_rad: float = None
"""Separatrix radiation fraction"""


e_plasma_beta: float = None
"""[J]"""


p_plasma_heating_total_mw: float = None
"""[W]"""


t_energy_confinement_beta: float = None
"""[s]"""


ptarmw: float = None


lambdaio: float = None


drsep: float = None


fio: float = None


fli: float = None


flo: float = None


fui: float = None


fuo: float = None


plimw: float = None


plomw: float = None


puimw: float = None


puomw: float = None


rho_star: float = None

nu_star: float = None

beta_mcdonald: float = None

itart_r: float = None

# Var in subroutine plasma_composition which requires re-initialisation on
# each new run:
first_call: int = None


# From physics_variables.f90:
N_CONFINEMENT_SCALINGS: int = 51
"""number of energy confinement time scaling laws"""

LABELS_CONFINEMENT_SCALINGS: list[str] = [
    "User input electron confinement   ",
    "Neo-Alcator                (Ohmic)",
    "Mirnov                         (H)",
    "Merezkhin-Muhkovatov    (Ohmic)(L)",
    "Shimomura                      (H)",
    "Kaye-Goldston                  (L)",
    "ITER 89-P                      (L)",
    "ITER 89-O                      (L)",
    "Rebut-Lallia                   (L)",
    "Goldston                       (L)",
    "T10                            (L)",
    "JAERI / Odajima-Shimomura      (L)",
    "Kaye-Big Complex               (L)",
    "ITER H90-P                     (H)",
    "ITER 89-P & 89-O min           (L)",
    "Riedel                         (L)",
    "Christiansen                   (L)",
    "Lackner-Gottardi               (L)",
    "Neo-Kaye                       (L)",
    "Riedel                         (H)",
    "ITER H90-P amended             (H)",
    "LHD                        (Stell)",
    "Gyro-reduced Bohm          (Stell)",
    "Lackner-Gottardi           (Stell)",
    "ITER-93H  ELM-free             (H)",
    "TITAN RFP OBSOLETE                ",
    "ITER H-97P ELM-free            (H)",
    "ITER H-97P ELMy                (H)",
    "ITER-96P (ITER-97L)            (L)",
    "Valovic modified ELMy          (H)",
    "Kaye 98 modified               (L)",
    "ITERH-PB98P(y)                 (H)",
    "IPB98(y)                       (H)",
    "IPB98(y,1)                     (H)",
    "IPB98(y,2)                     (H)",
    "IPB98(y,3)                     (H)",
    "IPB98(y,4)                     (H)",
    "ISS95                      (Stell)",
    "ISS04                      (Stell)",
    "DS03 beta-independent          (H)",
    'Murari "Non-power law"         (H)',
    "Petty 2008                 (ST)(H)",
    "Lang high density              (H)",
    "Hubbard 2017 - nominal         (I)",
    "Hubbard 2017 - lower           (I)",
    "Hubbard 2017 - upper           (I)",
    "Menard NSTX                (ST)(H)",
    "Menard NSTX-Petty08 hybrid (ST)(H)",
    "Buxton NSTX gyro-Bohm      (ST)(H)",
    "ITPA20                         (H)",
    "ITPA20-IL                      (H)",
]
"""labels describing energy confinement scaling laws"""


m_beam_amu: float = None
"""beam ion mass (amu)"""


m_fuel_amu: float = None
"""average mass of fuel portion of ions (amu)"""


m_ions_total_amu: float = None
"""average mass of all ions (amu)"""


m_plasma_fuel_ions: float = None
"""Mass of the plasma fuel ions (kg)"""


m_plasma_ions_total: float = None
"""Mass of all ions in plasma (kg)"""


m_plasma_alpha: float = None
"""Mass of the alpha particles in the plasma (kg)"""


m_plasma_electron: float = None
"""Mass of the electrons in the plasma (kg)"""


m_plasma: float = None
"""Total mass of the plasma (kg)"""


alphaj: float = None
"""current profile index"""


alphaj_wesson: float = None
"""Wesson-like current profile index"""


alphan: float = None
"""density profile index"""


alphap: float = None
"""pressure profile index"""


fusden_alpha_total: float = None
"""Alpha particle production rate per unit volume, from plasma and beams [particles/m3/sec]"""


fusden_plasma_alpha: float = None
"""Alpha particle production rate per unit volume, just from plasma [particles/m3/sec]"""


alphat: float = None
"""temperature profile index"""


aspect: float = None
"""aspect ratio (`iteration variable 1`)"""


beamfus0: float = None
"""multiplier for beam-background fusion calculation"""


beta_total_vol_avg: float = None
"""Volume averaged total plasma beta (`iteration variable 5`) (calculated if stellarator)"""


beta_fast_alpha: float = None
"""fast alpha beta component"""


beta_vol_avg_max: float = None
"""Max allowable volume averaged beta"""


beta_vol_avg_min: float = None
"""Minimum allowable volume averaged beta"""


beta_beam: float = None
"""neutral beam beta component"""


beta_poloidal_vol_avg: float = None
"""poloidal beta"""


beta_poloidal_eps: float = None
"""Poloidal beta and inverse aspcet ratio product"""


beta_toroidal_vol_avg: float = None
"""Plasma volume averaged toroidal beta"""

beta_thermal_toroidal_profile: list[float] = None
"""toroidal beta profile"""


beta_thermal_vol_avg: float = None
"""Plasma volume averaged thermal beta"""


beta_thermal_poloidal_vol_avg: float = None
"""Plasma volume averaged poloidal thermal beta"""


beta_thermal_toroidal_vol_avg: float = None
"""Plasma volume averaged toloidal thermal beta"""


beta_norm_total: float = None
"""normaised total beta"""


beta_norm_thermal: float = None
"""normaised thermal beta"""


beta_norm_toroidal: float = None
"""normaised toroidal beta"""


beta_norm_poloidal: float = None
"""normaised poloidal beta"""


e_plasma_beta_thermal: float = None
"""Plasma thermal energy derived from thermal beta"""


betbm0: float = None
"""leading coefficient for NB beta fraction"""


b_plasma_poloidal_average: float = None
"""Plasma average poloidal field (T)"""


b_plasma_toroidal_on_axis: float = None
"""Plasma toroidal field on axis (T) (`iteration variable 2`)"""

b_plasma_inboard_toroidal: float = None
"""Plasma inboard toroidal field (T)"""

b_plasma_outboard_toroidal: float = None
"""Plasma outboard toroidal field (T)"""

b_plasma_toroidal_profile: list[float] = None
"""toroidal field profile in plasma (T)"""


b_plasma_total: float = None
"""Sum of plasma total toroidal + poloidal field (T)"""

e_plasma_magnetic_stored: float = None
"""Plasma stored magnetic energy (J)"""


burnup: float = None
"""fractional plasma burnup"""


burnup_in: float = None
"""fractional plasma burnup user input"""


b_plasma_vertical_required: float = None
"""Vertical field needed for plasma equilibrium (T)"""


c_beta: float = None
"""Destabalisation parameter for i_beta_norm_max=4 beta limit"""


csawth: float = None
"""coeff. for sawteeth effects on burn V-s requirement"""


f_vol_plasma: float = None
"""multiplying factor for the plasma volume (normally=1)"""


f_r_conducting_wall: float = None
"""maximum ratio of conducting wall distance to plasma minor radius for
vertical stability (`constraint equation 23`)
"""


nd_plasma_electrons_vol_avg: float = None
"""Plasma volume averaged electron density (/m3) (`iteration variable 6`)"""


nd_plasma_fuel_ions_vol_avg: float = None
"""Plasma volume averaged fuel ion density (/m3)"""


dlamee: float = None
"""electron-electron coulomb logarithm"""


dlamie: float = None
"""ion-electron coulomb logarithm"""


nd_plasma_electron_max_array: list[float] = None
"""Array of plasma electron density upper limits values (/m3)"""


nd_plasma_alphas_vol_avg: float = None
"""Plasma volume averaged thermal alpha density (/m3)"""


nd_beam_ions: float = None
"""hot beam ion density, variable (/m3)"""


nd_beam_ions_out: float = None
"""hot beam ion density from calculation (/m3)"""


beta_norm_max: float = None
"""Troyon-like coefficient for beta scaling"""


beta_norm_max_wesson: float = None
"""Wesson-like coefficient for beta scaling"""


beta_norm_max_menard: float = None
"""Menard-like coefficient for beta scaling"""


beta_norm_max_original_scaling: float = None
"""Original scaling coefficient for beta scaling"""


beta_norm_max_tholerus: float = None
"""Tholerus-like coefficient for beta scaling"""


beta_norm_max_stambaugh: float = None
"""Stambaugh-like coefficient for beta scaling"""


nd_plasma_electrons_max: float = None
"""Plasma electron max density limit (/m3)"""


nd_plasma_ions_total_vol_avg: float = None
"""Plasma volume averaged total ion density (/m3)"""


nd_plasma_electron_line: float = None
"""Plasma line averaged electron density (/m3)"""


nd_plasma_protons_vol_avg: float = None
"""Plasma volume averaged proton ash density (/m3)"""


ntau: float = None
"""Fusion double product (s/m3)"""


nTtau: float = None
"""Lawson triple product [keV s / m3]"""


nd_plasma_impurities_vol_avg: float = None
"""Plasma volume averaged impurity (Z > 2) ion density (/m3)"""


gradient_length_ne: float = None
"""Max. normalized gradient length in el. density (i_plasma_pedestal==0 only)"""


gradient_length_te: float = None
"""Max. normalized gradient length in el. temperature (i_plasma_pedestal==0 only)"""


beta_poloidal_eps_max: float = None
"""maximum (eps*beta_poloidal) (`constraint equation 6`). Note: revised issue #346
"Operation at the tokamak equilibrium poloidal beta-limit in TFTR", 1992 Nucl. Fusion 32 1468
"""


eps: float = None
"""inverse aspect ratio"""


f_c_plasma_auxiliary: float = None
"""fraction of plasma current produced by auxiliary current drive"""


f_c_plasma_inductive: float = None
"""fraction of plasma current produced inductively"""


f_alpha_electron: float = None
"""fraction of alpha energy to electrons"""


f_p_alpha_plasma_deposited: float = None
"""Fraction of alpha power deposited in plasma. Default of 0.95 taken from https://doi.org/10.1088/0029-5515/39/12/305."""


f_alpha_ion: float = None
"""fraction of alpha power to ions"""


f_plasma_fuel_deuterium: float = None
"""Plasma deuterium fuel fraction"""


f_p_div_lower: float = None
"""fraction of power to the lower divertor in double null configuration
(`i_single_null = 0` only) (default assumes SN)
"""


ffwal: float = None
"""factor to convert plasma surface area to first wall area in neutron wall
load calculation (`i_pflux_fw_neutron=1`)
"""


f_nd_plasma_pedestal_greenwald: float = None
"""fraction of Greenwald density to set as pedestal-top density. If `<0`, pedestal-top
density set manually using nd_plasma_pedestal_electron (`i_plasma_pedestal==1`).
(`iteration variable 145`)
"""


f_nd_plasma_separatrix_greenwald: float = None
"""fraction of Greenwald density to set as separatrix density. If `<0`, separatrix
density set manually using nd_plasma_separatrix_electron (`i_plasma_pedestal==1`).
(`iteration variable 152`)
"""


f_plasma_fuel_helium3: float = None
"""Plasma Helium-3 fuel fraction"""


figmer: float = None
"""physics figure of merit (= plasma_current*aspect**sbar, where `sbar=1`)"""


fkzohm: float = None
"""Zohm elongation scaling adjustment factor (`i_plasma_geometry=2, 3`)"""


fplhsep: float = None
"""F-value for Psep >= Plh + Paux (`constraint equation 73`)"""


fp_plasma_separatrix_min_mw: float = None
"""F-value for minimum p_plasma_separatrix_mw (`constraint equation 80`)"""


fne0: float = None
"""f-value for the constraint ne(0) > ne(ped) (`constraint equation 81`)
(`Iteration variable 154`)
"""


f_plasma_fuel_tritium: float = None
"""Plasma tritium fuel fraction"""


fusden_total: float = None
"""fusion reaction rate density, from beams and plasma (reactions/m3/sec)"""


fusrat_total: float = None
"""fusion reaction rate, from beams and plasma (reactions/sec)"""

fusrat_plasma_dt_profile: list[float] = None
"""Profile of D-T fusion reaction rate in plasma, (reactions/sec)"""

fusrat_plasma_dd_triton_profile: list[float] = None
"""Profile of D-D fusion reaction rate (tritium branch) in plasma, (reactions/sec)"""

fusrat_plasma_dd_helion_profile: list[float] = None
"""Profile of D-D fusion reaction rate (helium branch) in plasma, (reactions/sec)"""

fusrat_plasma_dhe3_profile: list[float] = None
"""Profile of D-3He fusion reaction rate in plasma, (reactions/sec)"""

fusden_plasma: float = None
"""fusion reaction rate, just from plasma (reactions/m3/sec)"""


f_c_plasma_non_inductive: float = None
"""fraction of the plasma current produced by non-inductive means (`iteration variable 44`)"""


ejima_coeff: float = None
"""Ejima coefficient for resistive startup V-s formula"""


f_beta_alpha_beam_thermal: float = None
"""ratio of (fast alpha + neutral beam beta) to thermal beta"""


hfac: list[float] = None
"""H factors for an ignited plasma for each energy confinement time scaling law"""


hfact: float = None
"""H factor on energy confinement times, radiation corrected (`iteration variable 10`)."""


t_plasma_energy_confinement_max: float = None
"""Maximum allowed energy confinement time (s)"""


i_bootstrap_current: int = None
"""switch for bootstrap current scaling
- =1 ITER 1989 bootstrap scaling (high R/a only)
- =2 for Nevins et al general scaling
- =3 for Wilson et al numerical scaling
- =4 for Sauter et al scaling
- =5 for Sakai et al scaling
- =6 for ARIES scaling
- =7 for Andrade et al scaling
- =8 for Hoang et al scaling
- =9 for Wong et al scaling
- =10 for Gi-I et al scaling
- =11 for Gi-II et al scaling
- =12 for Sugiyama (L-mode) et al scaling
- =13 for Sugiyama (H-mode) et al scaling
"""


i_beta_component: int = None
"""switch for beta limit scaling (`constraint equation 24`)
- =0 apply limit to total beta
- =1 apply limit to thermal beta
- =2 apply limit to thermal + neutral beam beta
- =3 apply limit to toroidal beta
"""


i_plasma_current: int = None
"""switch for plasma current scaling to use
- =1 Peng analytic fit
- =2 Peng double null divertor scaling (ST)
- =3 simple ITER scaling (k = 2.2, d = 0.6)
- =4 later ITER scaling, a la Uckan
- =5 Todd empirical scaling I
- =6 Todd empirical scaling II
- =7 Connor-Hastie model
- =8 Sauter scaling allowing negative triangularity
- =9 FIESTA ST fit
"""


i_diamagnetic_current: int = None
"""switch for diamagnetic current scaling
- =0 Do not calculate
- =1 Use original TART scaling
- =2 Use SCENE scaling
"""


i_density_limit: int = None
"""switch for density limit to enforce (`constraint equation 5`)
- =1 old ASDEX
- =2 Borrass model for ITER (I)
- =3 Borrass model for ITER (II)
- =4 JET edge radiation
- =5 JET simplified
- =6 Hugill-Murakami Mq limit
- =7 Greenwald limit
- =8 ASDEX New
"""


n_divertors: int = None
"""number of divertors (calculated from `i_single_null`)"""


i_beta_fast_alpha: int = None
"""switch for fast alpha pressure calculation
- =0 ITER physics rules (Uckan) fit
- =1 Modified fit (D. Ward) - better at high temperature
"""


i_plasma_ignited: int = None
"""switch for ignition assumption. Obviously, i_plasma_ignited must be zero if current drive
is required. If i_plasma_ignited is 1, any auxiliary power is assumed to be used only during
plasma start-up, and is excluded from all steady-state power balance calculations.
- =0 do not assume plasma ignition
- =1 assume ignited (but include auxiliary power in costs)</UL
"""


i_plasma_pedestal: int = None
"""switch for pedestal profiles:
- =0 use original parabolic profiles
- =1 use pedestal profile
"""


i_pfirsch_schluter_current: int = None
"""switch for Pfirsch-Schlüter current scaling (issue #413):
- =0 Do not calculate
- =1 Use SCENE scaling
"""


nd_plasma_pedestal_electron: float = None
"""electron density of pedestal [m-3] (`i_plasma_pedestal==1)"""


nd_plasma_separatrix_electron: float = None
"""electron density at separatrix [m-3] (`i_plasma_pedestal==1)"""


alpha_crit: float = None
"""critical ballooning parameter value"""


nd_plasma_separatrix_electron_eich_max: float = None
"""Eich critical electron density at separatrix [m-3]"""


plasma_res_factor: float = None
"""plasma resistivity pre-factor"""


radius_plasma_pedestal_density_norm: float = None
"""Normalised radius of density pedestal (`i_plasma_pedestal==1`)"""


radius_plasma_pedestal_temp_norm: float = None
"""Normalised radius of temperature pedestal (`i_plasma_pedestal==1`)"""


rho_te_max: float = None
"""r/a where the temperature gradient is largest (`i_plasma_pedestal==0`)"""


rho_ne_max: float = None
"""r/a where the density gradient is largest (`i_plasma_pedestal==0`)"""


tbeta: float = None
"""temperature profile index beta  (`i_plasma_pedestal==1)"""


temp_plasma_pedestal_kev: float = None
"""Plasma electron temperature of pedestal (keV) (`i_plasma_pedestal==1`)"""


temp_plasma_separatrix_kev: float = None
"""Plasma electron temperature at separatrix (keV) (`i_plasma_pedestal==1`) calculated if reinke
criterion is used (`icc=78`)
"""


i_beta_norm_max: int = None
"""Switch for maximum normalised beta scaling:"""


i_ind_plasma_internal_norm: int = None
"""Switch for plasma normalised internal inductance scaling:"""


i_alphaj: int = None
"""Switch for current profile index scaling:"""


i_rad_loss: int = None
"""switch for radiation loss term usage in power balance (see User Guide):
- =0 total power lost is scaling power plus radiation
- =1 total power lost is scaling power plus core radiation only
- =2 total power lost is scaling power only, with no additional
allowance for radiation. This is not recommended for power plant models.
"""


i_confinement_time: int = None
"""switch for energy confinement time scaling law (see description in `LABELS_CONFINEMENT_SCALINGS`)"""


i_plasma_wall_gap: int = None
"""Switch for plasma-first wall clearances at the mid-plane:
- =0 use 10% of plasma minor radius
- =1 use input (`dr_fw_plasma_gap_inboard` and `dr_fw_plasma_gap_outboard`)
"""


i_plasma_geometry: int = None
"""switch for plasma elongation and triangularity calculations:
- =0 use input kappa, triang to calculate 95% values
- =1 scale q95_min, kappa, triang with aspect ratio (ST)
- =2 set kappa to the natural elongation value (Zohm ITER scaling), triang input
- =3 set kappa to the natural elongation value (Zohm ITER scaling), triang95 input
- =4 use input kappa95, triang95 to calculate separatrix values
- =5 use input kappa95, triang95 to calculate separatrix values based on MAST scaling (ST)
- =6 use input kappa, triang to calculate 95% values based on MAST scaling (ST)
- =7 use input kappa95, triang95 to calculate separatrix values based on fit to FIESTA (ST)
- =8 use input kappa, triang to calculate 95% values based on fit to FIESTA (ST)
- =9 set kappa to the natural elongation value, triang input
- =10 set kappa to maximum stable value at a given aspect ratio (2.6<A<3.6)), triang input (#1399)
- =11 set kappa Menard 2016 aspect-ratio-dependent scaling, triang input (#1439)
"""


i_plasma_shape: int = None
"""switch for plasma boundary shape:
- =0 use original PROCESS 2-arcs model
- =1 use the Sauter model
"""


itart: int = None
"""switch for spherical tokamak (ST) models:
- =0 use conventional aspect ratio models
- =1 use spherical tokamak models
"""


itartpf: int = None
"""switch for Spherical Tokamak PF models:
- =0 use Peng and Strickler (1986) model
- =1 use conventional aspect ratio model
"""


i_pflux_fw_neutron: int = None
"""switch for neutron wall load calculation:
- =1 use scaled plasma surface area
- =2 use first wall area directly
"""


plasma_square: float = None
"""plasma squareness used by Sauter plasma shape"""


kappa: float = None
"""plasma separatrix elongation (calculated if `i_plasma_geometry = 1-5, 7 or 9-10`)"""


kappa95: float = None
"""plasma elongation at 95% surface (calculated if `i_plasma_geometry = 0-3, 6, or 8-10`)"""


kappa_ipb: float = None
"""Separatrix elongation calculated for IPB scalings"""


nd_plasma_electron_on_axis: float = None
"""central electron density (/m3)"""


nd_plasma_ions_on_axis: float = None
"""central ion density (/m3)"""


m_s_limit: float = None
"""margin to vertical stability"""


pres_plasma_thermal_on_axis: float = None
"""Plasma central thermal pressure (no fast ions or beam pressure) (Pa)"""

pres_plasma_thermal_total_profile: list[float] = None
"""Profile of total pressure in plasma (Pa)"""

pres_plasma_electron_profile: list[float] = None
"""Profile of electron pressure in plasma (Pa)"""

pres_plasma_ion_total_profile: list[float] = None
"""Profile of ion pressure in plasma (Pa)"""

pres_plasma_fuel_profile: list[float] = None
"""Profile of fuel pressure in plasma (Pa)"""

j_plasma_on_axis: float = None
"""Central plasma current density (A/m2)"""

n_plasma_profile_elements: int = None
"""Number of elements in plasma profile"""

pres_plasma_thermal_vol_avg: float = None
"""Volume averaged thermal plasma pressure (no fast ions or beam pressure) (Pa)"""


f_dd_branching_trit: float = None
"""branching ratio for DD -> T"""


pden_plasma_alpha_mw: float = None
"""Alpha power per volume just from plasma [MW/m3]"""


pden_alpha_total_mw: float = None
"""Alpha power per volume from plasma and beams [MW/m3]"""


f_pden_alpha_electron_mw: float = None
"""Alpha power per volume to electrons [MW/m3]"""


p_fw_alpha_mw: float = None
"""alpha power escaping plasma and reaching first wall (MW)"""


f_pden_alpha_ions_mw: float = None
"""alpha power per volume to ions (MW/m3)"""


p_plasma_alpha_mw: float = None
"""Alpha power from only the plasma (MW)"""


p_alpha_total_mw: float = None
"""Total alpha power from plasma and beams (MW)"""


p_beam_alpha_mw: float = None
"""alpha power from hot neutral beam ions (MW)"""


p_beam_neutron_mw: float = None
"""neutron power from hot neutral beam ions (MW)"""


p_beam_dt_mw: float = None
"""D-T fusion power from hot neutral beam ions (MW)"""


p_non_alpha_charged_mw: float = None
"""non-alpha charged particle fusion power (MW)"""


p_charged_particle_mw: float = None
"""Total charged particle fusion power [MW]"""


pden_non_alpha_charged_mw: float = None
"""Non-alpha charged particle fusion power per volume [MW/m3]"""


f_temp_plasma_electron_density_vol_avg: float = None
"""Ratio of density weighted plasma electron tempertaurature to volume averaged (Profile Factor)"""


p_plasma_inner_rad_mw: float = None
"""radiation power from inner zone (MW)"""


pden_plasma_core_rad_mw: float = None
"""total core radiation power per volume (MW/m3)"""


p_dd_total_mw: float = None
"""deuterium-deuterium fusion power (MW)"""


p_dhe3_total_mw: float = None
"""deuterium-helium3 fusion power (MW)"""


p_plasma_separatrix_mw: float = None
"""power to conducted to the divertor region (MW)"""


p_div_lower_separatrix_mw: float = None
"""Separatrix power conducted to the lower divertor region (calculated if `i_single_null = 0`) (MW)"""


p_div_upper_separatrix_mw: float = None
"""Separatrix power conducted to the upper divertor region (calculated if `i_single_null = 0`) (MW)"""


p_div_separatrix_max_mw: float = None
"""Separatrix power conducted to the divertor with most load (calculated if `i_single_null = 0`) (MW)"""


p_dt_total_mw: float = None
"""Total deuterium-tritium fusion power, from plasma and beams [MW]"""


p_plasma_dt_mw: float = None
"""Deuterium-tritium fusion power, just from plasma [MW]"""


p_plasma_outer_rad_mw: float = None
"""radiation power from outer zone (MW)"""


pden_plasma_outer_rad_mw: float = None
"""edge radiation power per volume (MW/m3)"""


vs_plasma_internal: float = None
"""internal plasma V-s"""


pflux_fw_rad_mw: float = None
"""Nominal mean radiation load on inside surface of reactor (MW/m2)"""


pden_ion_electron_equilibration_mw: float = None
"""ion/electron equilibration power per volume (MW/m3)"""


plasma_current: float = None
"""plasma current (A)"""


p_plasma_neutron_mw: float = None
"""Neutron fusion power from just the plasma [MW]"""


p_neutron_total_mw: float = None
"""Total neutron fusion power from plasma and beams [MW]"""


pden_neutron_total_mw: float = None
"""neutron fusion power per volume from beams and plasma (MW/m3)"""


pden_plasma_neutron_mw: float = None
"""neutron fusion power per volume just from plasma (MW/m3)"""


p_plasma_ohmic_mw: float = None
"""ohmic heating power (MW)"""


pden_plasma_ohmic_mw: float = None
"""ohmic heating power per volume (MW/m3)"""


p_plasma_loss_mw: float = None
"""heating power (= transport loss power) (MW) used in confinement time calculation"""


p_fusion_total_mw: float = None
"""fusion power (MW)"""


len_plasma_poloidal: float = None
"""plasma poloidal perimeter (m)"""


p_plasma_rad_mw: float = None
"""total radiation power from inside LCFS (MW)"""


pden_plasma_rad_mw: float = None
"""total radiation power per volume (MW/m3)"""


pradsolmw: float = None
"""radiation power from SoL (MW)"""


proton_rate_density: float = None
"""Proton production rate [particles/m3/sec]"""


psolradmw: float = None
"""SOL radiation power (MW) (`stellarator only`)"""


pden_plasma_sync_mw: float = None
"""synchrotron radiation power per volume (MW/m3)"""


p_plasma_sync_mw: float = None
"""Total synchrotron radiation power from plasma (MW)"""


i_l_h_threshold: int = None
"""switch for L-H mode power threshold scaling to use (see l_h_threshold_powers for list)"""


p_l_h_threshold_mw: float = None
"""L-H mode power threshold (MW) (chosen via i_l_h_threshold, and enforced if
constraint equation 15 is on)
"""


l_h_threshold_powers: list[float] = None
"""L-H power threshold for various scalings (MW)
- =1 ITER 1996 scaling: nominal
- =2 ITER 1996 scaling: upper bound
- =3 ITER 1996 scaling: lower bound
- =4 ITER 1997 scaling: excluding elongation
- =5 ITER 1997 scaling: including elongation
- =6 Martin 2008 scaling: nominal
- =7 Martin 2008 scaling: 95% upper bound
- =8 Martin 2008 scaling: 95% lower bound
- =9 Snipes 2000 scaling: nominal
- =10 Snipes 2000 scaling: upper bound
- =11 Snipes 2000 scaling: lower bound
- =12 Snipes 2000 scaling (closed divertor): nominal
- =13 Snipes 2000 scaling (closed divertor): upper bound
- =14 Snipes 2000 scaling (closed divertor): lower bound
- =15 Hubbard et al. 2012 L-I threshold scaling: nominal
- =16 Hubbard et al. 2012 L-I threshold scaling: lower bound
- =17 Hubbard et al. 2012 L-I threshold scaling: upper bound
- =18 Hubbard et al. 2017 L-I threshold scaling
- =19 Martin 2008 aspect ratio corrected scaling: nominal
- =20 Martin 2008 aspect ratio corrected scaling: 95% upper bound
- =21 Martin 2008 aspect ratio corrected scaling: 95% lower bound
"""


p_electron_transport_loss_mw: float = None
"""electron transport power (MW)"""


pden_electron_transport_loss_mw: float = None
"""electron transport power per volume (MW/m3)"""


p_ion_transport_loss_mw: float = None
"""ion transport power (MW)"""


pscalingmw: float = None
"""Total transport power from scaling law (MW)"""


pden_ion_transport_loss_mw: float = None
"""ion transport power per volume (MW/m3)"""


q0: float = None
"""Safety factor on axis"""


q95: float = None
"""Safety factor at 95% flux surface (iteration variable 18) (unless icurr=2 (ST current scaling),
in which case q95 = mean edge safety factor qbar)
"""


molflow_plasma_fuelling_required: float = None
"""plasma fuelling rate (nucleus-pairs/s)"""


tauratio: float = None
"""tauratio /1.0/ : ratio of He and pellet particle confinement times"""


q95_min: float = None
"""lower limit for edge safety factor"""


qstar: float = None
"""cylindrical safety factor"""


rad_fraction_sol: float = None
"""SoL radiation fraction"""


rad_fraction_total: float = None
"""Radiation fraction total = SoL + LCFS radiation / total power deposited in plasma"""


f_nd_alpha_electron: float = None
"""thermal alpha density/electron density (`iteration variable 109`)"""


f_nd_protium_electrons: float = None
"""Seeded f_nd_protium_electrons density / electron density."""


ind_plasma_internal_norm: float = None
"""Plasma normalised internal inductance"""


ind_plasma_internal_norm_wesson: float = None
"""Wesson-like plasma normalised internal inductance"""


ind_plasma_internal_menard: float = None
"""Menard-like plasma normalised internal inductance"""


ind_plasma: float = None
"""plasma inductance (H)"""


rmajor: float = None
"""plasma major radius (m) (`iteration variable 3`)"""


rminor: float = None
"""plasma minor radius (m)"""


f_nd_beam_electron: float = None
"""hot beam density / n_e (`iteration variable 7`)"""


f_nd_plasma_carbon_electron: float = None
"""n_carbon / n_e"""


rndfuel: float = None
"""fuel burnup rate (reactions/second)"""


f_nd_plasma_iron_argon_electron: float = None
"""n_highZ / n_e"""


f_nd_plasma_oxygen_electron: float = None
"""n_oxygen / n_e"""


f_res_plasma_neo: float = None
"""neo-classical correction factor to res_plasma"""


res_plasma: float = None
"""plasma resistance (ohm)"""


t_plasma_res_diffusion: float = None
"""plasma current resistive diffusion time (s)"""


a_plasma_surface: float = None
"""plasma surface area"""


a_plasma_surface_outboard: float = None
"""outboard plasma surface area"""


i_single_null: int = None
"""switch for single null / double null plasma:
- =0 for double null
- =1 for single null (diverted side down)
"""


f_sync_reflect: float = None
"""synchrotron wall reflectivity factor"""


t_electron_energy_confinement: float = None
"""electron energy confinement time (sec)"""


tauee_in: float = None
"""Input electron energy confinement time (sec) (`i_confinement_time=48 only`)"""


t_energy_confinement: float = None
"""global thermal energy confinement time (sec)"""


t_ion_energy_confinement: float = None
"""ion energy confinement time (sec)"""


t_alpha_confinement: float = None
"""alpha particle confinement time (sec)"""


f_alpha_energy_confinement: float = None
"""alpha particle to energy confinement time ratio"""


temp_plasma_electron_vol_avg_kev: float = None
"""volume averaged electron temperature (keV) (`iteration variable 4`)"""


temp_plasma_electron_on_axis_kev: float = None
"""central electron temperature (keV)"""


temp_plasma_electron_density_weighted_kev: float = None
"""density weighted average electron temperature (keV)"""


temp_plasma_ion_vol_avg_kev: float = None
"""volume averaged ion temperature (keV). N.B. calculated from temp_plasma_electron_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`"""


temp_plasma_ion_on_axis_kev: float = None
"""central ion temperature (keV)"""


temp_plasma_ion_density_weighted_kev: float = None
"""density weighted average ion temperature (keV)"""


f_temp_plasma_ion_electron: float = None
"""ion temperature / electron temperature(used to calculate temp_plasma_ion_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`"""


triang: float = None
"""plasma separatrix triangularity (calculated if `i_plasma_geometry = 1, 3-5 or 7`)"""


triang95: float = None
"""plasma triangularity at 95% surface (calculated if `i_plasma_geometry = 0-2, 6, 8 or 9`)"""


vol_plasma: float = None
"""plasma volume (m3)"""


vs_plasma_burn_required: float = None
"""V-s needed during flat-top (heat + burn times) (Wb)"""


vs_plasma_ramp_required: float = None
"""V-s needed during ramp-up (Wb)"""


v_plasma_loop_burn: float = None
"""Plasma loop voltage during flat-top (V)"""


vs_plasma_ind_ramp: float = None
"""Total plasma inductive flux consumption for plasma current ramp-up (Vs)(Wb)"""


vs_plasma_res_ramp: float = None
"""Plasma resistive flux consumption for plasma current ramp-up (Vs)(Wb)"""


vs_plasma_total_required: float = None
"""total V-s needed (Wb)"""


pflux_fw_neutron_mw: float = None
"""average neutron wall load (MW/m2)"""


wtgpd: float = None
"""mass of fuel used per day (g)"""


a_plasma_poloidal: float = None
"""plasma poloidal cross-sectional area [m^2]"""


zeff: float = None
"""plasma effective charge"""


zeffai: float = None
"""mass weighted plasma effective charge"""


def init_physics_module():
    """Initialise the physics module"""
    global first_call
    global iscz
    global err242
    global err243
    global f_p_plasma_separatrix_rad
    global e_plasma_beta
    global p_plasma_heating_total_mw
    global t_energy_confinement_beta
    global ptarmw
    global lambdaio
    global drsep
    global fio
    global fli
    global flo
    global fui
    global fuo
    global plimw
    global plomw
    global puimw
    global puomw
    global rho_star
    global rho_ne_max
    global rho_te_max
    global nu_star
    global beta_mcdonald
    global itart_r

    first_call = 1
    iscz = 0
    err242 = 0
    err243 = 0
    f_p_plasma_separatrix_rad = 0.0
    e_plasma_beta = 0.0
    p_plasma_heating_total_mw = 0.0
    t_energy_confinement_beta = 0.0
    ptarmw = 0.0
    lambdaio = 0.0
    drsep = 0.0
    fio = 0.0
    fli = 0.0
    flo = 0.0
    fui = 0.0
    fuo = 0.0
    plimw = 0.0
    plomw = 0.0
    puimw = 0.0
    puomw = 0.0
    rho_ne_max = 0.0
    rho_te_max = 0.0
    rho_star = 0.0
    nu_star = 0.0
    beta_mcdonald = 0.0
    itart_r = 0.0


def init_physics_variables():
    global m_beam_amu
    global m_fuel_amu
    global m_ions_total_amu
    global m_plasma_fuel_ions
    global m_plasma_ions_total
    global m_plasma_alpha
    global m_plasma_electron
    global m_plasma
    global alphaj
    global i_alphaj
    global alphan
    global alphap
    global fusden_alpha_total
    global fusden_plasma_alpha
    global alphat
    global aspect
    global beamfus0
    global beta_total_vol_avg
    global beta_fast_alpha
    global beta_vol_avg_max
    global beta_vol_avg_min
    global beta_beam
    global beta_poloidal_vol_avg
    global beta_poloidal_eps
    global beta_toroidal_vol_avg
    global beta_thermal_toroidal_profile
    global beta_thermal_vol_avg
    global beta_thermal_poloidal_vol_avg
    global beta_thermal_toroidal_vol_avg
    global beta_norm_total
    global beta_norm_thermal
    global beta_norm_poloidal
    global e_plasma_beta_thermal
    global beta_norm_toroidal
    global betbm0
    global b_plasma_poloidal_average
    global b_plasma_toroidal_on_axis
    global b_plasma_toroidal_inboard
    global b_plasma_toroidal_outboard
    global b_plasma_toroidal_profile
    global b_plasma_total
    global e_plasma_magnetic_stored
    global burnup
    global burnup_in
    global b_plasma_vertical_required
    global c_beta
    global csawth
    global f_vol_plasma
    global f_r_conducting_wall
    global nd_plasma_electrons_vol_avg
    global nd_plasma_fuel_ions_vol_avg
    global dlamee
    global dlamie
    global nd_plasma_electron_max_array
    global nd_plasma_alphas_vol_avg
    global nd_beam_ions
    global nd_beam_ions_out
    global beta_norm_max
    global beta_norm_max_wesson
    global beta_norm_max_menard
    global beta_norm_max_original_scaling
    global beta_norm_max_tholerus
    global beta_norm_max_stambaugh
    global nd_plasma_electrons_max
    global nd_plasma_ions_total_vol_avg
    global nd_plasma_electron_line
    global nd_plasma_protons_vol_avg
    global ntau
    global nTtau
    global nd_plasma_impurities_vol_avg
    global beta_poloidal_eps_max
    global eps
    global f_c_plasma_auxiliary
    global f_c_plasma_inductive
    global f_alpha_electron
    global f_p_alpha_plasma_deposited
    global f_alpha_ion
    global f_plasma_fuel_deuterium
    global f_p_div_lower
    global ffwal
    global f_nd_plasma_pedestal_greenwald
    global f_nd_plasma_separatrix_greenwald
    global f_plasma_fuel_helium3
    global figmer
    global fkzohm
    global fplhsep
    global fp_plasma_separatrix_min_mw
    global fne0
    global f_plasma_fuel_tritium
    global fusden_total
    global fusrat_total
    global fusrat_plasma_dt_profile
    global fusrat_plasma_dd_triton_profile
    global fusrat_plasma_dd_helion_profile
    global fusrat_plasma_dhe3_profile
    global fusden_plasma
    global f_c_plasma_non_inductive
    global ejima_coeff
    global f_beta_alpha_beam_thermal
    global hfac
    global hfact
    global t_plasma_energy_confinement_max
    global i_bootstrap_current
    global i_beta_component
    global i_plasma_current
    global i_diamagnetic_current
    global i_density_limit
    global n_divertors
    global i_beta_fast_alpha
    global i_plasma_ignited
    global i_plasma_pedestal
    global i_pfirsch_schluter_current
    global nd_plasma_pedestal_electron
    global nd_plasma_separatrix_electron
    global alpha_crit
    global nd_plasma_separatrix_electron_eich_max
    global plasma_res_factor
    global radius_plasma_pedestal_density_norm
    global radius_plasma_pedestal_temp_norm
    global tbeta
    global temp_plasma_pedestal_kev
    global temp_plasma_separatrix_kev
    global i_beta_norm_max
    global i_rad_loss
    global i_confinement_time
    global i_plasma_wall_gap
    global i_plasma_geometry
    global i_plasma_shape
    global itart
    global itartpf
    global i_pflux_fw_neutron
    global plasma_square
    global kappa
    global kappa95
    global kappa_ipb
    global nd_plasma_electron_on_axis
    global nd_plasma_ions_on_axis
    global m_s_limit
    global pres_plasma_thermal_on_axis
    global pres_plasma_thermal_total_profile
    global pres_plasma_electron_profile
    global pres_plasma_ion_total_profile
    global pres_plasma_fuel_profile
    global j_plasma_on_axis
    global n_plasma_profile_elements
    global f_dd_branching_trit
    global pden_plasma_alpha_mw
    global pden_alpha_total_mw
    global f_pden_alpha_electron_mw
    global p_fw_alpha_mw
    global f_pden_alpha_ions_mw
    global p_alpha_total_mw
    global p_plasma_alpha_mw
    global p_beam_alpha_mw
    global p_beam_neutron_mw
    global p_beam_dt_mw
    global p_non_alpha_charged_mw
    global pden_non_alpha_charged_mw
    global f_temp_plasma_electron_density_vol_avg
    global p_plasma_inner_rad_mw
    global pden_plasma_core_rad_mw
    global p_dd_total_mw
    global p_dhe3_total_mw
    global p_plasma_separatrix_mw
    global p_div_lower_separatrix_mw
    global p_div_upper_separatrix_mw
    global p_div_separatrix_max_mw
    global p_dt_total_mw
    global p_plasma_dt_mw
    global p_plasma_outer_rad_mw
    global pden_plasma_outer_rad_mw
    global p_charged_particle_mw
    global vs_plasma_internal
    global pflux_fw_rad_mw
    global pden_ion_electron_equilibration_mw
    global plasma_current
    global p_plasma_neutron_mw
    global p_neutron_total_mw
    global pden_neutron_total_mw
    global pden_plasma_neutron_mw
    global p_plasma_ohmic_mw
    global pden_plasma_ohmic_mw
    global p_plasma_loss_mw
    global p_fusion_total_mw
    global len_plasma_poloidal
    global p_plasma_rad_mw
    global pden_plasma_rad_mw
    global pradsolmw
    global proton_rate_density
    global psolradmw
    global pden_plasma_sync_mw
    global p_plasma_sync_mw
    global i_l_h_threshold
    global p_l_h_threshold_mw
    global l_h_threshold_powers
    global p_electron_transport_loss_mw
    global pden_electron_transport_loss_mw
    global p_ion_transport_loss_mw
    global pscalingmw
    global pden_ion_transport_loss_mw
    global q0
    global q95
    global molflow_plasma_fuelling_required
    global tauratio
    global q95_min
    global qstar
    global rad_fraction_sol
    global rad_fraction_total
    global f_nd_alpha_electron
    global f_nd_protium_electrons
    global ind_plasma_internal_norm
    global ind_plasma_internal_norm_wesson
    global ind_plasma_internal_norm_menard
    global i_ind_plasma_internal_norm
    global ind_plasma
    global rmajor
    global rminor
    global f_nd_beam_electron
    global f_nd_plasma_carbon_electron
    global rndfuel
    global f_nd_plasma_iron_argon_electron
    global f_nd_plasma_oxygen_electron
    global f_res_plasma_neo
    global res_plasma
    global t_plasma_res_diffusion
    global a_plasma_surface
    global a_plasma_surface_outboard
    global i_single_null
    global f_sync_reflect
    global t_electron_energy_confinement
    global tauee_in
    global t_energy_confinement
    global t_ion_energy_confinement
    global t_alpha_confinement
    global f_alpha_energy_confinement
    global temp_plasma_electron_vol_avg_kev
    global temp_plasma_electron_on_axis_kev
    global temp_plasma_electron_density_weighted_kev
    global temp_plasma_ion_vol_avg_kev
    global temp_plasma_ion_on_axis_kev
    global temp_plasma_ion_density_weighted_kev
    global f_temp_plasma_ion_electron
    global triang
    global triang95
    global vol_plasma
    global vs_plasma_burn_required
    global vs_plasma_ramp_required
    global v_plasma_loop_burn
    global vs_plasma_ind_ramp
    global vs_plasma_res_ramp
    global vs_plasma_total_required
    global pflux_fw_neutron_mw
    global wtgpd
    global a_plasma_poloidal
    global zeff
    global zeffai

    m_beam_amu = 0.0
    m_fuel_amu = 0.0
    m_ions_total_amu = 0.0
    m_plasma_fuel_ions = 0.0
    m_plasma_ions_total = 0.0
    m_plasma_alpha = 0.0
    m_plasma_electron = 0.0
    m_plasma = 0.0
    alphaj = 1.0
    i_alphaj = 0
    alphan = 0.25
    alphap = 0.0
    fusden_alpha_total = 0.0
    fusden_plasma_alpha = 0.0
    alphat = 0.5
    aspect = 2.907
    beamfus0 = 1.0
    beta_total_vol_avg = 0.042
    beta_fast_alpha = 0.0
    beta_vol_avg_max = 0.0
    beta_vol_avg_min = 0.0
    beta_beam = 0.0
    beta_poloidal_vol_avg = 0.0
    beta_poloidal_eps = 0.0
    beta_toroidal_vol_avg = 0.0
    beta_thermal_toroidal_profile = []
    beta_thermal_vol_avg = 0.0
    beta_thermal_poloidal_vol_avg = 0.0
    beta_thermal_toroidal_vol_avg = 0.0
    beta_norm_total = 0.0
    beta_norm_thermal = 0.0
    beta_norm_poloidal = 0.0
    e_plasma_beta_thermal = 0.0
    beta_norm_toroidal = 0.0
    betbm0 = 1.5
    b_plasma_poloidal_average = 0.0
    b_plasma_toroidal_on_axis = 5.68
    b_plasma_toroidal_inboard = 0.0
    b_plasma_toroidal_outboard = 0.0
    b_plasma_toroidal_profile = []
    b_plasma_total = 0.0
    e_plasma_magnetic_stored = 0.0
    burnup = 0.0
    burnup_in = 0.0
    b_plasma_vertical_required = 0.0
    c_beta = 0.5
    csawth = 1.0
    f_vol_plasma = 1.0
    f_r_conducting_wall = 1.35
    nd_plasma_electrons_vol_avg = 9.8e19
    nd_plasma_fuel_ions_vol_avg = 0.0
    dlamee = 0.0
    dlamie = 0.0
    nd_plasma_electron_max_array = np.zeros(8, dtype=np.float64)
    nd_plasma_alphas_vol_avg = 0.0
    nd_beam_ions = 0.0
    nd_beam_ions_out = 0.0
    beta_norm_max = 3.5
    beta_norm_max_wesson = 0.0
    beta_norm_max_menard = 0.0
    beta_norm_max_original_scaling = 0.0
    beta_norm_max_tholerus = 0.0
    beta_norm_max_stambaugh = 0.0
    nd_plasma_electrons_max = 0.0
    nd_plasma_ions_total_vol_avg = 0.0
    nd_plasma_electron_line = 0.0
    nd_plasma_protons_vol_avg = 0.0
    ntau = 0.0
    nTtau = 0.0
    nd_plasma_impurities_vol_avg = 0.0
    beta_poloidal_eps_max = 1.38
    eps = 0.34399724802
    f_c_plasma_auxiliary = 0.0
    f_c_plasma_inductive = 0.0
    f_alpha_electron = 0.0
    f_p_alpha_plasma_deposited = 0.95
    f_alpha_ion = 0.0
    f_plasma_fuel_deuterium = 0.5
    f_p_div_lower = 1.0
    ffwal = 0.92
    f_nd_plasma_pedestal_greenwald = 0.85
    f_nd_plasma_separatrix_greenwald = 0.50
    f_plasma_fuel_helium3 = 0.0
    figmer = 0.0
    fkzohm = 1.0
    fplhsep = 1.0
    fp_plasma_separatrix_min_mw = 1.0
    fne0 = 1.0
    f_plasma_fuel_tritium = 0.5
    fusden_total = 0.0
    fusrat_total = 0.0
    fusrat_plasma_dt_profile = []
    fusrat_plasma_dd_triton_profile = []
    fusrat_plasma_dd_helion_profile = []
    fusrat_plasma_dhe3_profile = []
    fusden_plasma = 0.0
    f_c_plasma_non_inductive = 1.0
    ejima_coeff = 0.4
    f_beta_alpha_beam_thermal = 0.0
    hfac = np.zeros(N_CONFINEMENT_SCALINGS, dtype=np.float64)
    hfact = 1.0
    t_plasma_energy_confinement_max = 10.0
    i_bootstrap_current = 3
    i_beta_component = 0
    i_plasma_current = 4
    i_diamagnetic_current = 0
    i_density_limit = 8
    n_divertors = 2
    i_beta_fast_alpha = 1
    i_plasma_ignited = 0
    i_plasma_pedestal = 1
    i_pfirsch_schluter_current = 0
    nd_plasma_pedestal_electron = 4.0e19
    nd_plasma_separatrix_electron = 3.0e19
    alpha_crit = 0.0
    nd_plasma_separatrix_electron_eich_max = 0.0
    plasma_res_factor = 1.0
    radius_plasma_pedestal_density_norm = 1.0
    radius_plasma_pedestal_temp_norm = 1.0
    tbeta = 2.0
    temp_plasma_pedestal_kev = 1.0
    temp_plasma_separatrix_kev = 0.1
    i_beta_norm_max = 1
    i_rad_loss = 1
    i_confinement_time = 34
    i_plasma_wall_gap = 1
    i_plasma_geometry = 0
    i_plasma_shape = 0
    itart = 0
    itartpf = 0
    i_pflux_fw_neutron = 1
    plasma_square = 0.0
    kappa = 1.792
    kappa95 = 1.6
    kappa_ipb = 0.0
    nd_plasma_electron_on_axis = 0.0
    nd_plasma_ions_on_axis = 0.0
    m_s_limit = 0.3
    pres_plasma_thermal_on_axis = 0.0
    pres_plasma_thermal_total_profile = []
    pres_plasma_electron_profile = []
    pres_plasma_ion_total_profile = []
    pres_plasma_fuel_profile = []
    j_plasma_on_axis = 0.0
    n_plasma_profile_elements = 500
    f_dd_branching_trit = 0.0
    pden_plasma_alpha_mw = 0.0
    pden_alpha_total_mw = 0.0
    f_pden_alpha_electron_mw = 0.0
    p_fw_alpha_mw = 0.0
    f_pden_alpha_ions_mw = 0.0
    p_alpha_total_mw = 0.0
    p_plasma_alpha_mw = 0.0
    p_beam_alpha_mw = 0.0
    p_beam_neutron_mw = 0.0
    p_beam_dt_mw = 0.0
    p_non_alpha_charged_mw = 0.0
    pden_non_alpha_charged_mw = 0.0
    f_temp_plasma_electron_density_vol_avg = 0.0
    p_plasma_inner_rad_mw = 0.0
    pden_plasma_core_rad_mw = 0.0
    p_dd_total_mw = 0.0
    p_dhe3_total_mw = 0.0
    p_plasma_separatrix_mw = 0.0
    p_div_lower_separatrix_mw = 0.0
    p_div_upper_separatrix_mw = 0.0
    p_div_separatrix_max_mw = 0.0
    p_dt_total_mw = 0.0
    p_plasma_dt_mw = 0.0
    p_plasma_outer_rad_mw = 0.0
    pden_plasma_outer_rad_mw = 0.0
    p_charged_particle_mw = 0.0
    vs_plasma_internal = 0.0
    pflux_fw_rad_mw = 0.0
    pden_ion_electron_equilibration_mw = 0.0
    plasma_current = 0.0
    p_plasma_neutron_mw = 0.0
    p_neutron_total_mw = 0.0
    pden_neutron_total_mw = 0.0
    pden_plasma_neutron_mw = 0.0
    p_plasma_ohmic_mw = 0.0
    pden_plasma_ohmic_mw = 0.0
    p_plasma_loss_mw = 0.0
    p_fusion_total_mw = 0.0
    len_plasma_poloidal = 0.0
    p_plasma_rad_mw = 0.0
    pden_plasma_rad_mw = 0.0
    pradsolmw = 0.0
    proton_rate_density = 0.0
    psolradmw = 0.0
    pden_plasma_sync_mw = 0.0
    p_plasma_sync_mw = 0.0
    i_l_h_threshold = 19
    p_l_h_threshold_mw = 0.0
    l_h_threshold_powers = np.zeros(21, dtype=np.float64)
    p_electron_transport_loss_mw = 0.0
    pden_electron_transport_loss_mw = 0.0
    p_ion_transport_loss_mw = 0.0
    pscalingmw = 0.0
    pden_ion_transport_loss_mw = 0.0
    q0 = 1.0
    q95 = 0.0
    molflow_plasma_fuelling_required = 0.0
    tauratio = 1.0
    q95_min = 0.0
    qstar = 0.0
    rad_fraction_sol = 0.8
    rad_fraction_total = 0.0
    f_nd_alpha_electron = 0.1
    f_nd_protium_electrons = 0.0
    ind_plasma_internal_norm = 0.9
    ind_plasma_internal_norm_wesson = 0.0
    ind_plasma_internal_norm_menard = 0.0
    i_ind_plasma_internal_norm = 0
    ind_plasma = 0.0
    rmajor = 8.14
    rminor = 0.0
    f_nd_beam_electron = 0.005
    f_nd_plasma_carbon_electron = 0.0
    rndfuel = 0.0
    f_nd_plasma_iron_argon_electron = 0.0
    f_nd_plasma_oxygen_electron = 0.0
    f_res_plasma_neo = 0.0
    res_plasma = 0.0
    t_plasma_res_diffusion = 0.0
    a_plasma_surface = 0.0
    a_plasma_surface_outboard = 0.0
    i_single_null = 1
    f_sync_reflect = 0.6
    t_electron_energy_confinement = 0.0
    tauee_in = 0.0
    t_energy_confinement = 0.0
    t_ion_energy_confinement = 0.0
    t_alpha_confinement = 0.0
    f_alpha_energy_confinement = 0.0
    temp_plasma_electron_vol_avg_kev = 12.9
    temp_plasma_electron_on_axis_kev = 0.0
    temp_plasma_electron_density_weighted_kev = 0.0
    temp_plasma_ion_vol_avg_kev = 12.9
    temp_plasma_ion_on_axis_kev = 0.0
    temp_plasma_ion_density_weighted_kev = 0.0
    f_temp_plasma_ion_electron = 1.0
    triang = 0.36
    triang95 = 0.24
    vol_plasma = 0.0
    vs_plasma_burn_required = 0.0
    vs_plasma_ramp_required = 0.0
    v_plasma_loop_burn = 0.0
    vs_plasma_ind_ramp = 0.0
    vs_plasma_res_ramp = 0.0
    vs_plasma_total_required = 0.0
    pflux_fw_neutron_mw = 0.0
    wtgpd = 0.0
    a_plasma_poloidal = 0.0
    zeff = 0.0
    zeffai = 0.0
