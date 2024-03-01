 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_module

  !! Module containing tokamak plasma physics routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains all the primary plasma physics routines
  !! for a tokamak device.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  !  Module-level variables

  integer :: iscz
  integer :: err242, err243
  real(dp) :: rad_fraction_LCFS
  real(dp) :: total_plasma_internal_energy  ! [J]
  real(dp) :: total_loss_power        ! [W]
  real(dp) :: total_energy_conf_time  ! [s]
  real(dp) :: ptarmw, lambdaio, drsep
  real(dp) :: fio, fLI, fLO, fUI, fUO, pLImw, pLOmw, pUImw, pUOmw
  real(dp) :: rho_star
  real(dp) :: nu_star
  real(dp) :: beta_mcdonald
  real(dp) :: itart_r

  ! Var in subroutine plasma_composition which requires re-initialisation on
  ! each new run
  integer :: first_call

  contains

  subroutine init_physics_module
    !! Initialise module variables
    implicit none

    first_call = 1
    iscz = 0
    err242 = 0
    err243 = 0
    rad_fraction_LCFS = 0.0D0
    total_plasma_internal_energy = 0.0D0
    total_loss_power = 0.0D0
    total_energy_conf_time = 0.0D0
    ptarmw = 0.0D0
    lambdaio = 0.0D0
    drsep = 0.0D0
    fio = 0.0D0
    fLI = 0.0D0
    fLO = 0.0D0
    fUI = 0.0D0
    fUO = 0.0D0
    pLImw = 0.0D0
    pLOmw = 0.0D0
    pUImw = 0.0D0
    pUOmw = 0.0D0
    rho_star   = 0.0D0
    nu_star   = 0.0D0
    beta_mcdonald = 0.0D0
    itart_r = 0.0D0
  end subroutine init_physics_module

end module physics_module
