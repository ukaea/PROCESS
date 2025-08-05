! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pfcoil_module
  !! Module containing PF coil routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! author: R Kemp, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of the PF coil systems for a fusion power plant.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
   use pfcoil_variables, only: nfixmx, n_pf_groups_max, n_pf_coils_in_group_max, ngc2
   implicit none

   public

   !  Local variables

   integer :: nef,nfxf
   real(dp) :: ricpf, ssq0, sig_axial, sig_hoop
   real(dp) :: axial_force
   ! Private module variable arrays have variable dimensions; can't be wrapped
   ! with f2py if made public
   ! #TODO Temporarily hardcode dimensions in order to make public and wrap
   ! real(dp), dimension(nfixmx), private :: rfxf,zfxf,cfxf,xind
   ! real(dp), dimension(n_pf_groups_max,n_pf_coils_in_group_max), private :: rcls,zcls
   ! real(dp), dimension(n_pf_groups_max), private :: ccls,ccl0
   ! real(dp), dimension(ngc2), private :: bpf2
   ! real(dp), dimension(ngc2,3), private :: vsdum
   real(dp), dimension(64) :: rfxf,zfxf,cfxf,xind
   real(dp), dimension(10,2) :: rcls,zcls
   real(dp), dimension(10) :: ccls,ccl0
   real(dp), dimension(22) :: bpf2
   real(dp), dimension(22,3) :: vsdum

   ! pfcoil subroutine var requiring re-initialisation before each new run
   logical :: first_call
   ! outpf subroutine var requiring re-initialisation before each new run
   logical :: CSlimit

 end module pfcoil_module
