module current_drive_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the current drive system
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

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

  real(dp) :: bootipf
  !! bootstrap current fraction (enforced; see ibss)

  real(dp) :: bscfmax
  !! maximum fraction of plasma current from bootstrap; if `bscfmax < 0`,
  !! bootstrap fraction = abs(bscfmax)

  real(dp) :: bscf_iter89
  !! bootstrap current fraction, ITER 1989 model

  real(dp) :: bscf_nevins
  !! bootstrap current fraction, Nevins et al model

  real(dp) :: bscf_sauter
  !! bootstrap current fraction, Sauter et al model

  real(dp) :: bscf_wilson
  !! bootstrap current fraction, Wilson et al model

  real(dp) :: cboot
  !! bootstrap current fraction multiplier (`ibss=1`)

  real(dp) :: cnbeam
  !! neutral beam current (A)

  real(dp) :: diacf_hender
  !! diamagnetic current fraction, Hender fit

  real(dp) :: diacf_scene
  !! diamagnetic current fraction, SCENE fit

  real(dp) :: diaipf
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
  
  real(dp) :: enbeam
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

  real(dp) :: ftritbm
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

  real(dp) :: plasipf
  !! plasma driven current fraction (Bootstrap + Diamagnetic + PS)

  real(dp) :: plhybd
  !! lower hybrid injection power (MW)

  real(dp) :: pnbeam
  !! neutral beam injection power (MW)

  real(dp) :: porbitlossmw
  !! neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)

  real(dp) :: psipf
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

  contains

  subroutine init_current_drive_variables
    !! Initialise module variables
    implicit none

    beamwd = 0.58D0
    bigq = 0.0D0
    bootipf = 0.0D0
    bscfmax = 0.9D0
    bscf_iter89 = 0.0D0
    bscf_nevins = 0.0D0
    bscf_sauter = 0.0D0
    bscf_wilson = 0.0D0
    cboot = 1.0D0
    cnbeam = 0.0D0
    diacf_hender = 0.0D0
    diacf_scene = 0.0D0
    diaipf = 0.0D0
    echpwr = 0.0D0
    echwpow = 0.0D0
    effcd = 0.0D0
    harnum = 2.0
    wave_mode = 0
    enbeam = 1.0D3
    etacd = 0.0D0
    etacdfix = 0.0D0
    etaech = 0.3D0
    etalh = 0.3D0
    etanbi = 0.3D0
    fpion = 0.5D0
    pnbitot = 0.0D0
    pscf_scene = 0.0D0
    nbshinemw = 0.0D0
    feffcd = 1.0D0
    forbitloss = 0.0D0
    frbeam = 1.05D0
    ftritbm = 1.0D-6
    gamcd = 0.0D0
    gamma_ecrh = 0.35D0
    xi_ebw = 0.8D0
    iefrf = 5
    iefrffix = 0
    irfcd = 1
    nbshinef = 0.0D0
    nbshield = 0.5D0
    pheat = 0.0D0
    pheatfix = 0.0D0
    pinjalw = 150.0D0
    pinjemw = 0.0D0
    pinjimw = 0.0D0
    pinjmw = 0.0D0
    pinjfixmw = 0.0D0
    plasipf = 0.0D0
    plhybd = 0.0D0
    pnbeam = 0.0D0
    porbitlossmw = 0.0D0
    psipf = 0.0D0
    pwplh = 0.0D0
    pwpnb = 0.0D0
    rtanbeam = 0.0D0
    rtanmax = 0.0D0
    taubeam = 0.0D0
    tbeamin = 3.0D0
  end subroutine init_current_drive_variables
end module current_drive_variables
