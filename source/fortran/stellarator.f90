module stellarator_module

  !! Module containing stellarator routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of the first wall, blanket and shield components
  !! of a fusion power plant.

  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
   use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  ! scaling parameters to reference point.
  real(dp) :: f_n, f_r, f_aspect, f_b, f_i, f_a
  real(dp) :: f_coil_aspect, r_coil_major, r_coil_minor, f_coil_shape


  logical :: first_call = .true.
  logical :: first_call_stfwbs = .true.

end module stellarator_module
