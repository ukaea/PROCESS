module pulse_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the pulsed reactor model
  !!
  !!### References
  !!
  !! - Work File Notes in F/MPE/MOD/CAG/PROCESS/PULSE
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: bctmp
  !! first wall bulk coolant temperature (C)

  real(dp) :: dtstor
  !! maximum allowable temperature change in stainless steel thermal storage block (K) (`istore=3`)

  integer :: istore
  !! Switch for thermal storage method:
  !!
  !! - =1 option 1 of Electrowatt report, AEA FUS 205
  !! - =2 option 2 of Electrowatt report, AEA FUS 205
  !! - =3 stainless steel block

  integer :: itcycl
  !! Switch for first wall axial stress model:
  !!
  !! - =1 total axial constraint, no bending
  !! - =2 no axial constraint, no bending
  !! - =3 no axial constraint, bending

  integer :: lpulse
  !! Switch for reactor model:
  !!
  !! - =0 continuous operation
  !! - =1 pulsed operation

  contains

  subroutine init_pulse_variables
    !! Initialise module variables
    implicit none

    bctmp = 320.0D0
    dtstor = 300.0D0
    istore = 1
    itcycl = 1
    lpulse = 0
  end subroutine init_pulse_variables
end module pulse_variables
