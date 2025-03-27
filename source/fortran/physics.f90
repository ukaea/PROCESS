 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_module

  !! Module containing tokamak plasma physics routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains all the primary plasma physics routines
  !! for a tokamak device.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  !  Module-level variables

  integer :: iscz
  integer :: err242, err243
  real(dp) :: rad_fraction_LCFS
  real(dp) :: e_plasma_beta  ! [J]
  real(dp) :: total_loss_power        ! [W]
  real(dp) :: t_energy_confinement_beta  ! [s]
  real(dp) :: ptarmw, lambdaio, drsep
  real(dp) :: fio, fLI, fLO, fUI, fUO, pLImw, pLOmw, pUImw, pUOmw
  real(dp) :: rho_star
  real(dp) :: nu_star
  real(dp) :: beta_mcdonald
  real(dp) :: itart_r

  ! Var in subroutine plasma_composition which requires re-initialisation on
  ! each new run
  integer :: first_call
end module physics_module
