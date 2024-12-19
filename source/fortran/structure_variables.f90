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

  real(dp) :: m_components_cryo_cooled
  !! total mass of components at cryogenic temperatures (kg)

  real(dp) :: fncmass
  !! PF coil outer support fence mass (kg)

  real(dp) :: gsmass
  !! reactor core gravity support mass (kg)

  contains

  subroutine init_structure_variables
    !! Initialise module variables
    implicit none

    aintmass = 0.0D0
    clgsmass = 0.0D0
    m_components_cryo_cooled = 0.0D0
    fncmass = 0.0D0
    gsmass = 0.0D0
  end subroutine init_structure_variables
end module structure_variables
