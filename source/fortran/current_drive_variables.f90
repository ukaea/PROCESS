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

  real(dp) :: bootstrap_current_fraction
  !! bootstrap current fraction (enforced; see i_bootstrap_current)

  real(dp) :: bootstrap_current_fraction_max
  !! maximum fraction of plasma current from bootstrap; if `bootstrap_current_fraction_max < 0`,
  !! bootstrap fraction = abs(bootstrap_current_fraction_max)

  real(dp) :: bscf_iter89
  !! bootstrap current fraction, ITER 1989 model

  real(dp) :: bscf_nevins
  !! bootstrap current fraction, Nevins et al model

  real(dp) :: bscf_sauter
  !! bootstrap current fraction, Sauter et al model

  real(dp) :: bscf_wilson
  !! bootstrap current fraction, Wilson et al model

  real(dp) :: bscf_sakai
  !! Bootstrap current fraction, Sakai et al model

  real(dp) :: bscf_aries
  !! Bootstrap current fraction, ARIES model

  real(dp) :: bscf_andrade
  !! Bootstrap current fraction, Andrade et al model

  real(dp) :: bscf_hoang
  !! Bootstrap current fraction, Hoang et al model

  real(dp) :: bscf_wong
  !! Bootstrap current fraction, Wong et al model

  real(dp) :: bscf_gi_I
  !! Bootstrap current fraction, first Gi et al model

  real(dp) :: bscf_gi_II
  !! Bootstrap current fraction, second Gi et al model

  real(dp) :: cboot
  !! bootstrap current fraction multiplier

  real(dp) :: beam_current
  !! neutral beam current (A)

  real(dp) :: diacf_hender
  !! diamagnetic current fraction, Hender fit

  real(dp) :: diacf_scene
  !! diamagnetic current fraction, SCENE fit

  real(dp) :: diamagnetic_current_fraction
  !! diamagnetic current fraction

  real(dp) :: echpwr
  !! ECH power (MW)

  real(dp) :: echwpow
  !! ECH wall plug power (MW)

  real(dp) :: effcd
  !! current drive efficiency (A/W)

  real(dp) :: harnum
  !! cyclotron harmonic frequency number, used in cut-off function

  integer :: wave_mode
  !! Switch for ECRH wave mode :
  !!
  !!  - =0 O-mode
  !!  - =1 X-mode

  real(dp) :: beam_energy
  !! neutral beam energy (keV) (`iteration variable 19`)

  real(dp) :: etacd
  !! auxiliary power wall plug to injector efficiency

  real(dp) :: etacdfix
  !! secondary auxiliary power wall plug to injector efficiency

  real(dp) :: etaech
  !! ECH wall plug to injector efficiency

  real(dp) :: etalh
  !! lower hybrid wall plug to injector efficiency

  real(dp) :: etanbi
  !! neutral beam wall plug to injector efficiency

  real(dp) :: fpion
  !! fraction of beam energy to ions

  real(dp) :: pnbitot
  !! neutral beam power entering vacuum vessel

  real(dp) :: pscf_scene
  !! Pfirsch-Schlüter current fraction, SCENE fit

  real(dp) :: nbshinemw
  !! neutral beam shine-through power

  real(dp) :: feffcd
  !! current drive efficiency fudge factor (`iteration variable 47`)

  real(dp) :: forbitloss
  !! fraction of neutral beam power lost after ionisation but before
  !! thermalisation (orbit loss fraction)

  real(dp) :: frbeam
  !! R_tangential / R_major for neutral beam injection

  real(dp) :: f_tritium_beam
  !! fraction of beam that is tritium

  real(dp) :: gamcd
  !! normalised current drive efficiency (1.0e20 A/(W m^2))

  real(dp) :: gamma_ecrh
  !! User input ECRH gamma (1.0e20 A/(W m^2))

  real(dp) :: xi_ebw
  !! User scaling input for EBW plasma heating. Default 0.43

  integer :: iefrf
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

  integer :: iefrffix
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

  integer :: irfcd
  !! Switch for current drive calculation:
  !!
  !! - =0 turned off
  !! - =1 turned on

  real(dp) :: nbshinef
  !! neutral beam shine-through fraction

  real(dp) :: nbshield
  !! neutral beam duct shielding thickness (m)

  real(dp) :: pheat
  !! heating power not used for current drive (MW) (`iteration variable 11`)

  real(dp) :: pheatfix
  !! secondary fixed heating power not used for current drive (MW)

  real(dp) :: pinjalw
  !! maximum allowable value for injected power (MW) (`constraint equation 30`)

  real(dp) :: pinjemw
  !! auxiliary injected power to electrons (MW)

  real(dp) :: pinjimw
  !! auxiliary injected power to ions (MW)

  real(dp) :: pinjmw
  !! total auxiliary injected power (MW)

  real(dp)  :: pinjfixmw
  !! secondary total fixed auxiliary injected power (MW)

  real(dp) :: plasma_current_internal_fraction
  !! plasma current fraction driven internally (Bootstrap + Diamagnetic + PS)

  real(dp) :: plhybd
  !! lower hybrid injection power (MW)

  real(dp) :: pnbeam
  !! neutral beam injection power (MW)

  real(dp) :: porbitlossmw
  !! neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)

  real(dp) :: ps_current_fraction
  !! Pfirsch-Schlüter current fraction

  real(dp) :: pwplh
  !! lower hybrid wall plug power (MW)

  real(dp) :: pwpnb
  !! neutral beam wall plug power (MW)

  real(dp) :: rtanbeam
  !! neutral beam centreline tangency radius (m)

  real(dp) :: rtanmax
  !! maximum tangency radius for centreline of beam (m)

  real(dp) :: taubeam
  !! neutral beam e-decay lengths to plasma centre

  real(dp) :: tbeamin
  !! permitted neutral beam e-decay lengths to plasma centre
end module current_drive_variables
