module startup_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the plasma start-up model
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

  real(dp) :: ftaue
  !! factor in energy confinement time formula

  real(dp) :: gtaue
  !! offset term in energy confinement time scaling

  real(dp) :: nign
  !! electron density at ignition (start-up) (/m3)

  real(dp) :: ptaue
  !! exponent for density term in energy confinement time formula

  real(dp) :: qtaue
  !! exponent for temperature term in energy confinement time formula

  real(dp) :: rtaue
  !! exponent for power term in energy confinement time formula

  real(dp) :: tign
  !! electron temperature at ignition (start-up) (keV)

  contains

  subroutine init_startup_variables
    !! Initialise module variables
    implicit none

    ftaue = 0.0D0
    gtaue  = 0.0D0
    nign  = 0.0D0
    ptaue  = 0.0D0
    qtaue  = 0.0D0
    rtaue  = 0.0D0
    tign  = 0.0D0
  end subroutine init_startup_variables
end module startup_variables
