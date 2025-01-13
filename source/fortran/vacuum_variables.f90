module vacuum_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the vacuum system
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public
  character(len=6) :: vacuum_model
  !! switch for vacuum pumping model:
  !!
  !! - ='old' for old detailed ETR model
  !! - ='simple' for simple steady-state model with comparison to ITER cryopumps
  !#TODO: old and simple not suitable names.

  real(dp) :: niterpump
  !! number of high vacuum pumps (real number), each with the throughput of one
  !! ITER cryopump (50 Pa m3 s-1), all operating at the same time (`vacuum_model='simple'`)

  integer :: ntype
  !! switch for vacuum pump type:
  !!
  !! - =0 - for turbomolecular pump (magnetic bearing) with speed of 2.0 m3/s
  !!   (1.95 for N2, 1.8 for He, 1.8 for DT)
  !! - =1 - for compound cryopump with nominal speed of 10.0 m3/s
  !!   (9.0 for N2, 5.0 for He and 25.0 for DT)

  integer :: nvduct
  !! number of ducts (torus to pumps)

  real(dp) :: dlscal
  !! vacuum system duct length scaling

  real(dp) :: pbase
  !! base pressure during dwell before gas pre-fill(Pa)

  real(dp) :: prdiv
  !! divertor chamber pressure during burn (Pa)

  real(dp) :: pumptp
  !! Pump throughput (molecules/s) (default is ITER value)

  real(dp) :: rat
  !! plasma chamber wall outgassing rate (Pa-m/s)

  real(dp) :: tn
  !! neutral gas temperature in chamber (K)

  real(dp) :: vacdshm
  !! mass of vacuum duct shield (kg)

  real(dp) :: vcdimax
  !! diameter of duct passage (m)

  integer :: vpumpn
  !! number of high vacuum pumps

  integer :: dwell_pump
  !! switch for dwell pumping options:
  !!
  !! - =0 pumping only during t_between_pulse
  !! - =1 pumping only during t_precharge
  !! - =2 pumping during t_between_pulse + t_precharge

  ! The following are used in the Battes, Day and Rohde pump-down model
  ! See "Basic considerations on the pump-down time in the dwell phase of a pulsed fusion DEMO"
  ! http://dx.doi.org/10.1016/j.fusengdes.2015.07.011)(vacuum_model=simple')

  real(dp) :: pumpareafraction
  !! area of one pumping port as a fraction of plasma surface area

  real(dp) :: pumpspeedmax
  !! maximum pumping speed per unit area for deuterium & tritium, molecular flow

  real(dp) :: pumpspeedfactor
  !! effective pumping speed reduction factor due to duct impedance

  real(dp) :: initialpressure
  !! initial neutral pressure at the beginning of the dwell phase (Pa)

  real(dp) :: outgasindex
  !! outgassing decay index

  real(dp) :: outgasfactor
  !! outgassing prefactor kw: outgassing rate at 1 s per unit area (Pa m s-1)

  contains

  subroutine init_vacuum_variables
    !! Initialise module variables
    implicit none

    vacuum_model = 'old'
    niterpump = 0.0D0
    ntype = 1
    nvduct = 0
    dlscal = 0.0D0
    pbase = 5.0D-4
    prdiv = 0.36D0
    pumptp = 1.2155D22
    rat = 1.3D-8
    tn = 300.0D0
    vacdshm = 0.0D0
    vcdimax = 0.0D0
    vpumpn = 0
    dwell_pump = 0
    pumpareafraction = 0.0203D0
    pumpspeedmax = 27.3D0
    pumpspeedfactor = 0.167D0
    initialpressure = 1.0D0
    outgasindex = 1.0D0
    outgasfactor = 0.0235D0
  end subroutine init_vacuum_variables
end module vacuum_variables
