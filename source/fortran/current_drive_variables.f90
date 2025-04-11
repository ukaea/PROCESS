module current_drive_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the current drive system
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: beamwd
  !! width of neutral beam duct where it passes between the TF coils (m)
  !! T Inoue et al, Design of neutral beam system for ITER-FEAT,
  !! <A HREF=http://dx.doi.org/10.1016/S0920-3796(01)00339-8>
  !! Fusion Engineering and Design, Volumes 56-57, October 2001, Pages 517-521</A>)

  real(dp) :: bigq
  !! Fusion gain; P_fusion / (P_injection + P_ohmic)

  real(dp) :: f_c_plasma_bootstrap
  !! bootstrap current fraction (enforced; see i_bootstrap_current)

  real(dp) :: f_c_plasma_bootstrap_max
  !! maximum fraction of plasma current from bootstrap; if `f_c_plasma_bootstrap_max < 0`,
  !! bootstrap fraction = abs(f_c_plasma_bootstrap_max)

  real(dp) :: f_c_plasma_bootstrap_iter89
  !! bootstrap current fraction, ITER 1989 model

  real(dp) :: f_c_plasma_bootstrap_nevins
  !! bootstrap current fraction, Nevins et al model

  real(dp) :: f_c_plasma_bootstrap_sauter
  !! bootstrap current fraction, Sauter et al model

  real(dp) :: f_c_plasma_bootstrap_wilson
  !! bootstrap current fraction, Wilson et al model

  real(dp) :: f_c_plasma_bootstrap_sakai
  !! Bootstrap current fraction, Sakai et al model

  real(dp) :: f_c_plasma_bootstrap_aries
  !! Bootstrap current fraction, ARIES model

  real(dp) :: f_c_plasma_bootstrap_andrade
  !! Bootstrap current fraction, Andrade et al model

  real(dp) :: f_c_plasma_bootstrap_hoang
  !! Bootstrap current fraction, Hoang et al model

  real(dp) :: f_c_plasma_bootstrap_wong
  !! Bootstrap current fraction, Wong et al model

  real(dp) :: bscf_gi_I
  !! Bootstrap current fraction, first Gi et al model

  real(dp) :: bscf_gi_II
  !! Bootstrap current fraction, second Gi et al model

  real(dp) :: cboot
  !! bootstrap current fraction multiplier

  real(dp) :: c_beam_total
  !! neutral beam current (A)

  real(dp) :: f_c_plasma_diamagnetic_hender
  !! diamagnetic current fraction, Hender fit

  real(dp) :: f_c_plasma_diamagnetic_scene
  !! diamagnetic current fraction, SCENE fit

  real(dp) :: f_c_plasma_diamagnetic
  !! diamagnetic current fraction

  real(dp) :: p_ecrh_injected_mw
  !! ECH power (MW)

  real(dp) :: echwpow
  !! ECH wall plug power (MW)

  real(dp) :: eta_cd_hcd_primary
  !! current drive efficiency (A/W)

  real(dp) :: n_ecrh_harmonic
  !! cyclotron harmonic frequency number, used in cut-off function

  integer :: i_ecrh_wave_mode
  !! Switch for ECRH wave mode :
  !!
  !!  - =0 O-mode
  !!  - =1 X-mode

  real(dp) :: e_beam_kev
  !! neutral beam energy (keV) (`iteration variable 19`)

  real(dp) :: eta_hcd_primary_injector_wall_plug
  !! auxiliary power wall plug to injector efficiency

  real(dp) :: eta_hcd_secondary_injector_wall_plug
  !! secondary auxiliary power wall plug to injector efficiency

  real(dp) :: eta_ecrh_injector_wall_plug
  !! ECH wall plug to injector efficiency

  real(dp) :: eta_lowhyb_injector_wall_plug
  !! lower hybrid wall plug to injector efficiency

  real(dp) :: eta_beam_injector_wall_plug
  !! neutral beam wall plug to injector efficiency

  real(dp) :: f_p_beam_injected_ions
  !! fraction of beam energy to ions

  real(dp) :: p_beam_injected_mw
  !! neutral beam power entering vacuum vessel

  real(dp) :: f_c_plasma_pfirsch_schluter_scene
  !! Pfirsch-Schlüter current fraction, SCENE fit

  real(dp) :: p_beam_shine_through_mw
  !! neutral beam shine-through power

  real(dp) :: feffcd
  !! current drive efficiency fudge factor (`iteration variable 47`)

  real(dp) :: f_p_beam_orbit_loss
  !! fraction of neutral beam power lost after ionisation but before
  !! thermalisation (orbit loss fraction)

  real(dp) :: frbeam
  !! R_tangential / R_major for neutral beam injection

  real(dp) :: f_beam_tritium
  !! fraction of beam that is tritium

  real(dp) :: eta_cd_norm_hcd_primary
  !! Normalised current drive efficiency for primary HCD system [(1.0e20 A)/(W m^2)]

  real(dp) :: eta_cd_norm_hcd_secondary
  !! Normalised current drive efficiency for secondary HCD system [(1.0e20 A)/(W m^2)]

  real(dp) :: eta_cd_norm_ecrh
  !! User input ECRH gamma (1.0e20 A/(W m^2))

  real(dp) :: xi_ebw
  !! User scaling input for EBW plasma heating. Default 0.43

  integer :: i_hcd_primary
  !! Switch for current drive efficiency model:
  !!
  !!  - =1 Fenstermacher Lower Hybrid
  !!  - =2 Ion Cyclotron current drive
  !!  - =3 Fenstermacher ECH
  !!  - =4 Ehst Lower Hybrid
  !!  - =5 ITER Neutral Beam
  !!  - =6 new Culham Lower Hybrid model
  !!  - =7 new Culham ECCD model
  !!  - =8 new Culham Neutral Beam model
  !!  - =9 RFP option removed in PROCESS (issue #508)
  !!  - =10 ECRH user input gamma
  !!  - =11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.
  !!  - =12 EBW user scaling input. Scaling (S. Freethy)

  integer :: i_hcd_secondary
  !! Switch for 2nd current drive efficiency model:
  !!
  !! - =0 No fixed current drive
  !! - =1 Fenstermacher Lower Hybrid
  !! - =2 Ion Cyclotron current drive
  !! - =3 Fenstermacher ECH
  !! - =4 Ehst Lower Hybrid
  !! - =5 ITER Neutral Beam
  !! - =6 new Culham Lower Hybrid model
  !! - =7 new Culham ECCD model
  !! - =8 new Culham Neutral Beam model
  !! - =9 RFP option removed in PROCESS (issue #508)
  !! - =10 ECRH user input gamma
  !! - =11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.
  !! - =12 EBW user scaling input. Scaling (S. Freethy)

  integer :: i_hcd_calculations
  !! Switch for current drive calculation:
  !!
  !! - =0 turned off
  !! - =1 turned on

  real(dp) :: f_p_beam_shine_through
  !! neutral beam shine-through fraction

  real(dp) :: dx_beam_shield
  !! neutral beam duct shielding thickness (m)

  real(dp) :: p_hcd_primary_extra_heat_mw
  !! heating power not used for current drive (MW) (`iteration variable 11`)

  real(dp) :: p_hcd_secondary_extra_heat_mw
  !! secondary fixed heating power not used for current drive (MW)

  real(dp) :: p_hcd_injected_max
  !! maximum allowable value for injected power (MW) (`constraint equation 30`)

  real(dp) :: p_hcd_injected_electrons_mw
  !! auxiliary injected power to electrons (MW)

  real(dp) :: p_hcd_injected_ions_mw
  !! auxiliary injected power to ions (MW)

  real(dp) :: p_hcd_injected_total_mw
  !! total auxiliary injected power (MW)

  real(dp)  :: p_hcd_secondary_injected_mw
  !! secondary total fixed auxiliary injected power (MW)

  real(dp) :: f_c_plasma_internal
  !! plasma current fraction driven internally (Bootstrap + Diamagnetic + PS)

  real(dp) :: plhybd
  !! lower hybrid injection power (MW)

  real(dp) :: pnbeam
  !! neutral beam injection power (MW)

  real(dp) :: p_beam_orbit_loss_mw
  !! neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)

  real(dp) :: f_c_plasma_pfirsch_schluter
  !! Pfirsch-Schlüter current fraction

  real(dp) :: pwplh
  !! lower hybrid wall plug power (MW)

  real(dp) :: pwpnb
  !! neutral beam wall plug power (MW)

  real(dp) :: rtanbeam
  !! neutral beam centreline tangency radius (m)

  real(dp) :: rtanmax
  !! maximum tangency radius for centreline of beam (m)

  real(dp) :: n_beam_decay_lengths_core
  !! neutral beam e-decay lengths to plasma centre

  real(dp) :: tbeamin
  !! permitted neutral beam e-decay lengths to plasma centre
end module current_drive_variables
