module structure_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the support structure
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: aintmass
  !! intercoil structure mass (kg)

  real(dp) :: clgsmass
  !! gravity support structure for TF coil, PF coil and intercoil support systems (kg)

  real(dp) :: coldmass
  !! total mass of components at cryogenic temperatures (kg)

  real(dp) :: fncmass
  !! PF coil outer support fence mass (kg)

  real(dp) :: gsmass
  !! reactor core gravity support mass (kg)
end module structure_variables
