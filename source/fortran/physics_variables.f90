module physics_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the plasma physics
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  integer, parameter :: n_confinement_scalings = 51
  !! number of energy confinement time scaling laws

  real(dp) :: m_beam_amu
  !! beam ion mass (amu)

  real(dp) :: m_fuel_amu
  !! average mass of fuel portion of ions (amu)

  real(dp) :: m_ions_total_amu
  !! average mass of all ions (amu)

  real(dp) :: m_plasma_fuel_ions
  !! Mass of the plasma fuel ions (kg)

  real(dp) :: m_plasma_ions_total
  !! Mass of all ions in plasma (kg)

  real(dp) :: m_plasma_alpha
  !! Mass of the alpha particles in the plasma (kg)

  real(dp) :: m_plasma_electron
  !! Mass of the electrons in the plasma (kg)

  real(dp) :: m_plasma
  !! Total mass of the plasma (kg)

  real(dp) :: alphaj
  !! current profile index

  real(dp) :: alphaj_wesson
  !! Wesson-like current profile index

  real(dp) :: alphan
  !! density profile index

  real(dp) :: alphap
  !! pressure profile index

  real(dp) :: fusden_alpha_total
  !! Alpha particle production rate per unit volume, from plasma and beams [particles/m3/sec]

  real(dp) :: fusden_plasma_alpha
  !! Alpha particle production rate per unit volume, just from plasma [particles/m3/sec]

  real(dp) :: alphat
  !! temperature profile index

  real(dp) :: aspect
  !! aspect ratio (`iteration variable 1`)

  real(dp) :: beamfus0
  !! multiplier for beam-background fusion calculation

  real(dp) :: beta
  !! total plasma beta (`iteration variable 5`) (calculated if stellarator)

  real(dp) :: beta_fast_alpha
  !! fast alpha beta component

  real(dp) :: beta_max
  !! Max allowable beta

  real(dp) :: beta_min
  !! allowable lower beta

  real(dp) :: beta_beam
  !! neutral beam beta component

  real(dp) :: beta_poloidal
  !! poloidal beta

  real(dp) :: beta_poloidal_eps
  !! Poloidal beta and inverse aspcet ratio product

  real(dp) :: beta_toroidal
  !! toroidal beta

  real(dp) :: beta_thermal
  !! thermal beta

  real(dp) :: beta_thermal_poloidal
  !! poloidal thermal beta

  real(dp) :: beta_thermal_toroidal
  !! poloidal thermal beta

  real(dp) :: beta_norm_total
  !! normaised total beta

  real(dp) :: beta_norm_thermal
  !! normaised thermal beta

  real(dp) :: beta_norm_toroidal
  !! normaised toroidal beta

  real(dp) :: beta_norm_poloidal
  !! normaised poloidal beta

  real(dp) :: e_plasma_beta_thermal
  !! Plasma thermal energy derived from thermal beta

  real(dp) :: betbm0
  !! leading coefficient for NB beta fraction

  real(dp) :: bp
  !! poloidal field (T)

  real(dp) :: bt
  !! toroidal field on axis (T) (`iteration variable 2`)

  real(dp) :: btot
  !! total toroidal + poloidal field (T)

  real(dp) :: burnup
  !! fractional plasma burnup

  real(dp) :: burnup_in
  !! fractional plasma burnup user input

  real(dp) :: bvert
  !! vertical field at plasma (T)

  real(dp) :: c_beta
  !! Destabalisation parameter for i_beta_norm_max=4 beta limit

  real(dp) :: csawth
  !! coeff. for sawteeth effects on burn V-s requirement

  real(dp) :: f_vol_plasma
  !! multiplying factor for the plasma volume (normally=1)

  real(dp) :: f_r_conducting_wall
  !! maximum ratio of conducting wall distance to plasma minor radius for
  !! vertical stability (`constraint equation 23`)

  real(dp) :: dene
  !! electron density (/m3) (`iteration variable 6`)

  real(dp) :: nd_fuel_ions
  !! fuel ion density (/m3)

  real(dp) :: dlamee
  !! electron-electron coulomb logarithm

  real(dp) :: dlamie
  !! ion-electron coulomb logarithm

  real(dp), dimension(8) :: dlimit
  !! density limit (/m3) as calculated using various models

  real(dp) :: nd_alphas
  !! thermal alpha density (/m3)

  real(dp) :: nd_beam_ions
  !! hot beam ion density, variable (/m3)

  real(dp) :: beam_density_out
  !! hot beam ion density from calculation (/m3)

  real(dp) :: beta_norm_max
  !! Troyon-like coefficient for beta scaling

  real(dp) :: beta_norm_max_wesson
  !! Wesson-like coefficient for beta scaling

  real(dp) :: beta_norm_max_menard
  !! Menard-like coefficient for beta scaling

  real(dp) :: beta_norm_max_original_scaling
  !! Original scaling coefficient for beta scaling

  real(dp) :: beta_norm_max_tholerus
  !! Tholerus-like coefficient for beta scaling

  real(dp) :: beta_norm_max_stambaugh
  !! Stambaugh-like coefficient for beta scaling

  real(dp) :: dnelimt
  !! density limit (/m3)

  real(dp) :: nd_ions_total
  !! total ion density (/m3)

  real(dp) :: dnla
  !! line averaged electron density (/m3)

  real(dp) :: nd_protons
  !! proton ash density (/m3)

  real(dp) :: ntau
  !! Fusion double product (s/m3)

  real(dp) :: nTtau
  !! Lawson triple product [keV s / m3]

  real(dp) :: nd_impurities
  !! high Z ion density (/m3)

  real(dp) :: gradient_length_ne
  !! Max. normalized gradient length in el. density (ipedestal==0 only)

  real(dp) :: gradient_length_te
  !! Max. normalized gradient length in el. temperature (ipedestal==0 only)

  real(dp) :: beta_poloidal_eps_max
  !! maximum (eps*beta_poloidal) (`constraint equation 6`). Note: revised issue #346
  !! "Operation at the tokamak equilibrium poloidal beta-limit in TFTR", 1992 Nucl. Fusion 32 1468

  real(dp) :: eps
  !! inverse aspect ratio

  real(dp) :: f_c_plasma_auxiliary
  !! fraction of plasma current produced by auxiliary current drive

  real(dp) :: f_c_plasma_inductive
  !! fraction of plasma current produced inductively

  real(dp) :: f_alpha_electron
  !! fraction of alpha energy to electrons

  real(dp) :: f_alpha_plasma
  !! Fraction of alpha power deposited in plasma. Default of 0.95 taken from https://doi.org/10.1088/0029-5515/39/12/305.

  real(dp) :: f_alpha_ion
  !! fraction of alpha power to ions

  real(dp) :: f_deuterium
  !! deuterium fuel fraction

  real(dp) :: f_p_div_lower
  !! fraction of power to the lower divertor in double null configuration
  !! (`i_single_null = 0` only) (default assumes SN)

  real(dp) :: ffwal
  !! factor to convert plasma surface area to first wall area in neutron wall
  !! load calculation (`iwalld=1`)

  real(dp) :: fgwped
  !! fraction of Greenwald density to set as pedestal-top density. If `<0`, pedestal-top
  !! density set manually using neped (`ipedestal==1`).
  !! (`iteration variable 145`)

  real(dp) :: fgwsep
  !! fraction of Greenwald density to set as separatrix density. If `<0`, separatrix
  !! density set manually using nesep (`ipedestal==1`).
  !! (`iteration variable 152`)

  real(dp) :: f_helium3
  !! helium-3 fuel fraction

  real(dp) :: figmer
  !! physics figure of merit (= plasma_current*aspect**sbar, where `sbar=1`)

  real(dp) :: fkzohm
  !! Zohm elongation scaling adjustment factor (`i_plasma_geometry=2, 3`)

  real(dp) :: fplhsep
  !! F-value for Psep >= Plh + Paux (`constraint equation 73`)

  real(dp) :: fp_plasma_separatrix_min_mw
  !! F-value for minimum p_plasma_separatrix_mw (`constraint equation 80`)

  real(dp) :: fne0
  !! f-value for the constraint ne(0) > ne(ped) (`constraint equation 81`)
  !! (`Iteration variable 154`)

  real(dp) :: f_tritium
  !! tritium fuel fraction

  real(dp) :: fusden_total
  !! fusion reaction rate density, from beams and plasma (reactions/m3/sec)

  real(dp) :: fusrat_total
  !! fusion reaction rate, from beams and plasma (reactions/sec)

  real(dp) :: fusden_plasma
  !! fusion reaction rate, just from plasma (reactions/m3/sec)

  real(dp) :: f_c_plasma_non_inductive
  !! fraction of the plasma current produced by non-inductive means (`iteration variable 44`)

  real(dp) :: ejima_coeff
  !! Ejima coefficient for resistive startup V-s formula

  real(dp) :: f_beta_alpha_beam_thermal
  !! ratio of (fast alpha + neutral beam beta) to thermal beta

  real(dp), dimension(n_confinement_scalings) :: hfac
  !! H factors for an ignited plasma for each energy confinement time scaling law

  real(dp) :: hfact
  !! H factor on energy confinement times, radiation corrected (`iteration variable 10`).

  real(dp) :: taumax
  !! Maximum allowed energy confinement time (s)

  integer :: i_bootstrap_current
  !! switch for bootstrap current scaling
  !!
  !! - =1 ITER 1989 bootstrap scaling (high R/a only)
  !! - =2 for Nevins et al general scaling
  !! - =3 for Wilson et al numerical scaling
  !! - =4 for Sauter et al scaling
  !! - =5 for Sakai et al scaling
  !! - =6 for ARIES scaling
  !! - =7 for Andrade et al scaling
  !! - =8 for Hoang et al scaling
  !! - =9 for Wong et al scaling
  !! - =10 for Gi-I et al scaling
  !! - =11 for Gi-II et al scaling
  !! - =12 for Sugiyama (L-mode) et al scaling
  !! - =13 for Sugiyama (H-mode) et al scaling

  integer :: i_beta_component
  !! switch for beta limit scaling (`constraint equation 24`)
  !!
  !! - =0 apply limit to total beta
  !! - =1 apply limit to thermal beta
  !! - =2 apply limit to thermal + neutral beam beta
  !! - =3 apply limit to toroidal beta

  integer :: i_plasma_current
  !! switch for plasma current scaling to use
  !!
  !! - =1 Peng analytic fit
  !! - =2 Peng double null divertor scaling (ST)
  !! - =3 simple ITER scaling (k = 2.2, d = 0.6)
  !! - =4 later ITER scaling, a la Uckan
  !! - =5 Todd empirical scaling I
  !! - =6 Todd empirical scaling II
  !! - =7 Connor-Hastie model
  !! - =8 Sauter scaling allowing negative triangularity
  !! - =9 FIESTA ST fit

  integer :: i_diamagnetic_current
  !! switch for diamagnetic current scaling
  !!
  !! - =0 Do not calculate
  !! - =1 Use original TART scaling
  !! - =2 Use SCENE scaling

  integer :: i_density_limit
  !! switch for density limit to enforce (`constraint equation 5`)
  !!
  !! - =1 old ASDEX
  !! - =2 Borrass model for ITER (I)
  !! - =3 Borrass model for ITER (II)
  !! - =4 JET edge radiation
  !! - =5 JET simplified
  !! - =6 Hugill-Murakami Mq limit
  !! - =7 Greenwald limit
  !! - =8 ASDEX New

  integer :: n_divertors
  !! number of divertors (calculated from `i_single_null`)

  integer :: i_beta_fast_alpha
  !! switch for fast alpha pressure calculation
  !!
  !! - =0 ITER physics rules (Uckan) fit
  !! - =1 Modified fit (D. Ward) - better at high temperature

  integer :: i_plasma_ignited
  !! switch for ignition assumption. Obviously, i_plasma_ignited must be zero if current drive
  !! is required. If i_plasma_ignited is 1, any auxiliary power is assumed to be used only during
  !! plasma start-up, and is excluded from all steady-state power balance calculations.
  !!
  !! - =0 do not assume plasma ignition
  !! - =1 assume ignited (but include auxiliary power in costs)</UL

  integer :: ipedestal
  !! switch for pedestal profiles:
  !!
  !! - =0 use original parabolic profiles
  !! - =1 use pedestal profile

  integer :: i_pfirsch_schluter_current
  !! switch for Pfirsch-Schlüter current scaling (issue #413):
  !!
  !! - =0 Do not calculate
  !! - =1 Use SCENE scaling

  real(dp) :: neped
  !! electron density of pedestal [m-3] (`ipedestal==1)

  real(dp) :: nesep
  !! electron density at separatrix [m-3] (`ipedestal==1)

  real(dp) :: alpha_crit
  !! critical ballooning parameter value

  real(dp) :: nesep_crit
  !! critical electron density at separatrix [m-3]

  real(dp) :: plasma_res_factor
  !! plasma resistivity pre-factor

  real(dp) :: rhopedn
  !! r/a of density pedestal (`ipedestal==1`)

  real(dp) :: rhopedt
  !! r/a of temperature pedestal (`ipedestal==1`)

  real(dp) :: rho_te_max
  !! r/a where the temperature gradient is largest (`ipedestal==0`)

  real(dp) :: rho_ne_max
  !! r/a where the density gradient is largest (`ipedestal==0`)

  real(dp) :: tbeta
  !! temperature profile index beta  (`ipedestal==1)

  real(dp) :: teped
  !! electron temperature of pedestal (keV) (`ipedestal==1`)

  real(dp) :: tesep
  !! electron temperature at separatrix (keV) (`ipedestal==1`) calculated if reinke
  !! criterion is used (`icc=78`)

  integer :: i_beta_norm_max
  !! Switch for maximum normalised beta scaling:

  integer :: i_ind_plasma_internal_norm
  !! Switch for plasma normalised internal inductance scaling:

  integer :: i_alphaj
  !! Switch for current profile index scaling:

  integer :: i_rad_loss
  !! switch for radiation loss term usage in power balance (see User Guide):
  !!
  !! - =0 total power lost is scaling power plus radiation
  !! - =1 total power lost is scaling power plus core radiation only
  !! - =2 total power lost is scaling power only, with no additional
  !!   allowance for radiation. This is not recommended for power plant models.

  integer :: i_confinement_time
  !! switch for energy confinement time scaling law (see description in `labels_confinement_scalings`)

  !! labels_confinement_scalings(n_confinement_scalings) : labels describing energy confinement scaling laws
  character*34, parameter, dimension(n_confinement_scalings) :: labels_confinement_scalings = (/  &
    'User input electron confinement   ', &
    'Neo-Alcator                (Ohmic)', &
    'Mirnov                         (H)', &
    'Merezkhin-Muhkovatov    (Ohmic)(L)', &
    'Shimomura                      (H)', &
    'Kaye-Goldston                  (L)', &
    'ITER 89-P                      (L)', &
    'ITER 89-O                      (L)', &
    'Rebut-Lallia                   (L)', &
    'Goldston                       (L)', &
    'T10                            (L)', &
    'JAERI / Odajima-Shimomura      (L)', &
    'Kaye-Big Complex               (L)', &
    'ITER H90-P                     (H)', &
    'ITER 89-P & 89-O min           (L)', &
    'Riedel                         (L)', &
    'Christiansen                   (L)', &
    'Lackner-Gottardi               (L)', &
    'Neo-Kaye                       (L)', &
    'Riedel                         (H)', &
    'ITER H90-P amended             (H)', &
    'LHD                        (Stell)', &
    'Gyro-reduced Bohm          (Stell)', &
    'Lackner-Gottardi           (Stell)', &
    'ITER-93H  ELM-free             (H)', &
    'TITAN RFP OBSOLETE                ', &
    'ITER H-97P ELM-free            (H)', &
    'ITER H-97P ELMy                (H)', &
    'ITER-96P (ITER-97L)            (L)', &
    'Valovic modified ELMy          (H)', &
    'Kaye 98 modified               (L)', &
    'ITERH-PB98P(y)                 (H)', &
    'IPB98(y)                       (H)', &
    'IPB98(y,1)                     (H)', &
    'IPB98(y,2)                     (H)', &
    'IPB98(y,3)                     (H)', &
    'IPB98(y,4)                     (H)', &
    'ISS95                      (Stell)', &
    'ISS04                      (Stell)', &
    'DS03 beta-independent          (H)', &
    'Murari "Non-power law"         (H)', &
    'Petty 2008                 (ST)(H)', &
    'Lang high density              (H)', &
    'Hubbard 2017 - nominal         (I)', &
    'Hubbard 2017 - lower           (I)', &
    'Hubbard 2017 - upper           (I)', &
    'Menard NSTX                (ST)(H)', &
    'Menard NSTX-Petty08 hybrid (ST)(H)', &
    'Buxton NSTX gyro-Bohm      (ST)(H)', &
    'ITPA20                         (H)', &
    'ITPA20-IL                      (H)' /)

  integer :: i_plasma_wall_gap
  !! Switch for plasma-first wall clearances at the mid-plane:
  !!
  !! - =0 use 10% of plasma minor radius
  !! - =1 use input (`dr_fw_plasma_gap_inboard` and `dr_fw_plasma_gap_outboard`)

  integer :: i_plasma_geometry
  !! switch for plasma elongation and triangularity calculations:
  !!
  !! - =0 use input kappa, triang to calculate 95% values
  !! - =1 scale q95_min, kappa, triang with aspect ratio (ST)
  !! - =2 set kappa to the natural elongation value (Zohm ITER scaling), triang input
  !! - =3 set kappa to the natural elongation value (Zohm ITER scaling), triang95 input
  !! - =4 use input kappa95, triang95 to calculate separatrix values
  !! - =5 use input kappa95, triang95 to calculate separatrix values based on MAST scaling (ST)
  !! - =6 use input kappa, triang to calculate 95% values based on MAST scaling (ST)
  !! - =7 use input kappa95, triang95 to calculate separatrix values based on fit to FIESTA (ST)
  !! - =8 use input kappa, triang to calculate 95% values based on fit to FIESTA (ST)
  !! - =9 set kappa to the natural elongation value, triang input
  !! - =10 set kappa to maximum stable value at a given aspect ratio (2.6<A<3.6)), triang input (#1399)
  !! - =11 set kappa Menard 2016 aspect-ratio-dependent scaling, triang input (#1439)

  integer :: i_plasma_shape
  !! switch for plasma boundary shape:
  !!
  !! - =0 use original PROCESS 2-arcs model
  !! - =1 use the Sauter model

  integer :: itart
  !! switch for spherical tokamak (ST) models:
  !!
  !! - =0 use conventional aspect ratio models
  !! - =1 use spherical tokamak models

  integer :: itartpf
  !! switch for Spherical Tokamak PF models:
  !!
  !! - =0 use Peng and Strickler (1986) model
  !! - =1 use conventional aspect ratio model

  integer :: iwalld
  !! switch for neutron wall load calculation:
  !!
  !! - =1 use scaled plasma surface area
  !! - =2 use first wall area directly

  real(dp) :: plasma_square
  !! plasma squareness used by Sauter plasma shape

  real(dp) :: kappa
  !! plasma separatrix elongation (calculated if `i_plasma_geometry = 1-5, 7 or 9-10`)

  real(dp) :: kappa95
  !! plasma elongation at 95% surface (calculated if `i_plasma_geometry = 0-3, 6, or 8-10`)


  real(dp) :: kappa_ipb
  !! Separatrix elongation calculated for IPB scalings

  real(dp) :: ne0
  !! central electron density (/m3)

  real(dp) :: ni0
  !! central ion density (/m3)

  real(dp) :: m_s_limit
  !! margin to vertical stability

  real(dp) :: p0
  !! central total plasma pressure (Pa)

  real(dp) :: j_plasma_0
  !! Central plasma current density (A/m2)

  real(dp) :: vol_avg_pressure
  !! Volume averaged plasma pressure (Pa)

  real(dp) :: f_dd_branching_trit
  !! branching ratio for DD -> T

  real(dp) :: pden_plasma_alpha_mw
  !! Alpha power per volume just from plasma [MW/m3]

  real(dp) :: pden_alpha_total_mw
  !! Alpha power per volume from plasma and beams [MW/m3]

  real(dp) :: f_pden_alpha_electron_mw
  !! Alpha power per volume to electrons [MW/m3]

  real(dp) :: p_fw_alpha_mw
  !! alpha power escaping plasma and reaching first wall (MW)

  real(dp) :: f_pden_alpha_ions_mw
  !! alpha power per volume to ions (MW/m3)

  real(dp) :: p_plasma_alpha_mw
  !! Alpha power from only the plasma (MW)

  real(dp) :: p_alpha_total_mw
  !! Total alpha power from plasma and beams (MW)

  real(dp) :: p_beam_alpha_mw
  !! alpha power from hot neutral beam ions (MW)

  real(dp) :: p_beam_neutron_mw
  !! neutron power from hot neutral beam ions (MW)

  real(dp) :: p_beam_dt_mw
  !! D-T fusion power from hot neutral beam ions (MW)

  real(dp) :: p_non_alpha_charged_mw
  !! non-alpha charged particle fusion power (MW)

  real(dp) :: p_charged_particle_mw
  !! Total charged particle fusion power [MW]

  real(dp) :: pden_non_alpha_charged_mw
  !! Non-alpha charged particle fusion power per volume [MW/m3]

  real(dp) :: pcoef
  !! profile factor (= n-weighted T / average T)

  real(dp) :: p_plasma_inner_rad_mw
  !! radiation power from inner zone (MW)

  real(dp) :: pden_plasma_core_rad_mw
  !! total core radiation power per volume (MW/m3)

  real(dp) :: p_dd_total_mw
  !! deuterium-deuterium fusion power (MW)

  real(dp) :: p_dhe3_total_mw
  !! deuterium-helium3 fusion power (MW)

  real(dp) :: p_plasma_separatrix_mw
  !! power to conducted to the divertor region (MW)

  real(dp) :: pdivl
  !! power conducted to the lower divertor region (calculated if `i_single_null = 0`) (MW)

  real(dp) :: pdivu
  !! power conducted to the upper divertor region (calculated if `i_single_null = 0`) (MW)

  real(dp) :: pdivmax
  !! power conducted to the divertor with most load (calculated if `i_single_null = 0`) (MW)

  real(dp) :: p_dt_total_mw
  !!  Total deuterium-tritium fusion power, from plasma and beams [MW]

  real(dp) :: p_plasma_dt_mw
  !!  Deuterium-tritium fusion power, just from plasma [MW]

  real(dp) :: p_plasma_outer_rad_mw
  !! radiation power from outer zone (MW)

  real(dp) :: pden_plasma_outer_rad_mw
  !! edge radiation power per volume (MW/m3)

  real(dp) :: vs_plasma_internal
  !! internal plasma V-s

  real(dp) :: pflux_fw_rad_mw
  !! Nominal mean radiation load on inside surface of reactor (MW/m2)

  real(dp) :: pden_ion_electron_equilibration_mw
  !! ion/electron equilibration power per volume (MW/m3)

  real(dp) :: plasma_current
  !! plasma current (A)

  real(dp) :: p_plasma_neutron_mw
  !! Neutron fusion power from just the plasma [MW]

  real(dp) :: p_neutron_total_mw
  !! Total neutron fusion power from plasma and beams [MW]

  real(dp) :: pden_neutron_total_mw
  !! neutron fusion power per volume from beams and plasma (MW/m3)

  real(dp) :: pden_plasma_neutron_mw
  !! neutron fusion power per volume just from plasma (MW/m3)

  real(dp) :: p_plasma_ohmic_mw
  !! ohmic heating power (MW)

  real(dp) :: pden_plasma_ohmic_mw
  !! ohmic heating power per volume (MW/m3)

  real(dp) :: p_plasma_loss_mw
  !! heating power (= transport loss power) (MW) used in confinement time calculation

  real(dp) :: p_fusion_total_mw
  !! fusion power (MW)

  real(dp) :: len_plasma_poloidal
  !! plasma poloidal perimeter (m)

  real(dp) :: p_plasma_rad_mw
  !! total radiation power from inside LCFS (MW)

  real(dp) :: pden_plasma_rad_mw
  !! total radiation power per volume (MW/m3)

  real(dp) :: pradsolmw
  !! radiation power from SoL (MW)

  real(dp) :: proton_rate_density
  !! Proton production rate [particles/m3/sec]

  real(dp) :: psolradmw
  !! SOL radiation power (MW) (`stellarator only`)

  real(dp) :: pden_plasma_sync_mw
  !! synchrotron radiation power per volume (MW/m3)

  real(dp) :: p_plasma_sync_mw
  !! Total synchrotron radiation power from plasma (MW)

  integer :: i_l_h_threshold
  !! switch for L-H mode power threshold scaling to use (see l_h_threshold_powers for list)

  real(dp) :: p_l_h_threshold_mw
  !! L-H mode power threshold (MW) (chosen via i_l_h_threshold, and enforced if
  !! constraint equation 15 is on)

  real(dp), dimension(21) :: l_h_threshold_powers
  !! L-H power threshold for various scalings (MW)
  !!
  !! - =1 ITER 1996 scaling: nominal
  !! - =2 ITER 1996 scaling: upper bound
  !! - =3 ITER 1996 scaling: lower bound
  !! - =4 ITER 1997 scaling: excluding elongation
  !! - =5 ITER 1997 scaling: including elongation
  !! - =6 Martin 2008 scaling: nominal
  !! - =7 Martin 2008 scaling: 95% upper bound
  !! - =8 Martin 2008 scaling: 95% lower bound
  !! - =9 Snipes 2000 scaling: nominal
  !! - =10 Snipes 2000 scaling: upper bound
  !! - =11 Snipes 2000 scaling: lower bound
  !! - =12 Snipes 2000 scaling (closed divertor): nominal
  !! - =13 Snipes 2000 scaling (closed divertor): upper bound
  !! - =14 Snipes 2000 scaling (closed divertor): lower bound
  !! - =15 Hubbard et al. 2012 L-I threshold scaling: nominal
  !! - =16 Hubbard et al. 2012 L-I threshold scaling: lower bound
  !! - =17 Hubbard et al. 2012 L-I threshold scaling: upper bound
  !! - =18 Hubbard et al. 2017 L-I threshold scaling
  !! - =19 Martin 2008 aspect ratio corrected scaling: nominal
  !! - =20 Martin 2008 aspect ratio corrected scaling: 95% upper bound
  !! - =21 Martin 2008 aspect ratio corrected scaling: 95% lower bound

  real(dp) :: p_electron_transport_loss_mw
  !! electron transport power (MW)

  real(dp) :: pden_electron_transport_loss_mw
  !! electron transport power per volume (MW/m3)

  real(dp) :: p_ion_transport_loss_mw
  !! ion transport power (MW)

  real(dp) :: pscalingmw
  !! Total transport power from scaling law (MW)

  real(dp) :: pden_ion_transport_loss_mw
  !! ion transport power per volume (MW/m3)

  real(dp) :: q0
  !! Safety factor on axis

  real(dp) :: q95
  !! Safety factor at 95% flux surface (iteration variable 18) (unless icurr=2 (ST current scaling),
  !! in which case q95 = mean edge safety factor qbar)

  real(dp) :: qfuel
  !! plasma fuelling rate (nucleus-pairs/s)

  real(dp) :: tauratio
  !! tauratio /1.0/ : ratio of He and pellet particle confinement times

  real(dp) :: q95_min
  !! lower limit for edge safety factor

  real(dp) :: qstar
  !! cylindrical safety factor

  real(dp) :: rad_fraction_sol
  !! SoL radiation fraction

  real(dp) :: rad_fraction_total
  !! Radiation fraction total = SoL + LCFS radiation / total power deposited in plasma

  real(dp) :: f_nd_alpha_electron
  !! thermal alpha density/electron density (`iteration variable 109`)

  real(dp) :: f_nd_protium_electrons
  !! Seeded f_nd_protium_electrons density / electron density.

  real(dp) :: ind_plasma_internal_norm
  !! Plasma normalised internal inductance

  real(dp) :: ind_plasma_internal_norm_wesson
  !! Wesson-like plasma normalised internal inductance

  real(dp) :: ind_plasma_internal_menard
  !! Menard-like plasma normalised internal inductance

  real(dp) :: ind_plasma
  !! plasma inductance (H)

  real(dp) :: rmajor
  !! plasma major radius (m) (`iteration variable 3`)

  real(dp) :: rminor
  !! plasma minor radius (m)

  real(dp) :: f_nd_beam_electron
  !! hot beam density / n_e (`iteration variable 7`)

  real(dp) :: rncne
  !! n_carbon / n_e

  real(dp) :: rndfuel
  !! fuel burnup rate (reactions/second)

  real(dp) :: rnfene
  !! n_highZ / n_e

  real(dp) :: rnone
  !! n_oxygen / n_e

  real(dp) :: f_res_plasma_neo
  !! neo-classical correction factor to res_plasma

  real(dp) :: res_plasma
  !! plasma resistance (ohm)

  real(dp) :: t_plasma_res_diffusion
  !! plasma current resistive diffusion time (s)

  real(dp) :: a_plasma_surface
  !! plasma surface area

  real(dp) :: a_plasma_surface_outboard
  !! outboard plasma surface area

  integer :: i_single_null
  !! switch for single null / double null plasma:
  !!
  !! - =0 for double null
  !! - =1 for single null (diverted side down)

  real(dp) :: f_sync_reflect
  !! synchrotron wall reflectivity factor

  real(dp) :: t_electron_energy_confinement
  !! electron energy confinement time (sec)

  real(dp) :: tauee_in
  !! Input electron energy confinement time (sec) (`i_confinement_time=48 only`)

  real(dp) :: t_energy_confinement
  !! global thermal energy confinement time (sec)

  real(dp) :: t_ion_energy_confinement
  !! ion energy confinement time (sec)

  real(dp) :: t_alpha_confinement
  !! alpha particle confinement time (sec)

  real(dp) :: f_alpha_energy_confinement
  !! alpha particle to energy confinement time ratio

  real(dp) :: te
  !! volume averaged electron temperature (keV) (`iteration variable 4`)

  real(dp) :: te0
  !! central electron temperature (keV)

  real(dp) :: ten
  !! density weighted average electron temperature (keV)

  real(dp) :: ti
  !! volume averaged ion temperature (keV). N.B. calculated from te if `tratio > 0.0`

  real(dp) :: ti0
  !! central ion temperature (keV)

  real(dp) :: tin
  !! density weighted average ion temperature (keV)

  real(dp) :: tratio
  !! ion temperature / electron temperature(used to calculate ti if `tratio > 0.0`

  real(dp) :: triang
  !! plasma separatrix triangularity (calculated if `i_plasma_geometry = 1, 3-5 or 7`)

  real(dp) :: triang95
  !! plasma triangularity at 95% surface (calculated if `i_plasma_geometry = 0-2, 6, 8 or 9`)

  real(dp) :: vol_plasma
  !! plasma volume (m3)

  real(dp) :: vs_plasma_burn_required
  !! V-s needed during flat-top (heat + burn times) (Wb)

  real(dp) :: vs_plasma_ramp_required
  !! V-s needed during ramp-up (Wb)

  real(dp) :: v_plasma_loop_burn
  !! Plasma loop voltage during flat-top (V)

  real(dp) :: vshift
  !! plasma/device midplane vertical shift - single null

  real(dp) :: vs_plasma_ind_ramp
  !! Total plasma inductive flux consumption for plasma current ramp-up (Vs)(Wb)

  real(dp) :: vs_plasma_res_ramp
  !! Plasma resistive flux consumption for plasma current ramp-up (Vs)(Wb)

  real(dp) :: vs_plasma_total_required
  !! total V-s needed (Wb)

  real(dp) :: pflux_fw_neutron_mw
  !! average neutron wall load (MW/m2)

  real(dp) :: wtgpd
  !! mass of fuel used per day (g)

  real(dp) :: a_plasma_poloidal
  !! plasma poloidal cross-sectional area [m^2]

  real(dp) :: zeff
  !! plasma effective charge

  real(dp) :: zeffai
  !! mass weighted plasma effective charge
end module physics_variables
