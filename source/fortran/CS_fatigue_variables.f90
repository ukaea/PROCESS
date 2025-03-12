module CS_fatigue_variables

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: residual_sig_hoop
  !! residual hoop stress in strucutal material (Pa)

  real(dp) :: n_cycle
  !! Allowable number of cycles for CS stress model

  real(dp) :: n_cycle_min
  !! Minimum llowable number of cycles for CS stress model

  real(dp) :: t_crack_radial
  !! Initial depth of crack in thickness of conduit (m)

  real(dp) :: t_crack_vertical
  !! Inital vertical crack size (m)

  real(dp) :: t_structural_radial
  !! Thickness of CS conductor conduit (m)

  real(dp) :: t_structural_vertical
  !! Vertical thickness of CS conductor conduit (m)

  real(dp) :: bkt_life_csf
  !! Switch to pass bkt_life cycles to n_cycle_min

  real(dp) :: sf_vertical_crack
  !! Safety factor for vertical crack size (-)

  real(dp) :: sf_radial_crack
  !! Safety factor for radial crack size (-)

  real(dp) :: sf_fast_fracture
  !! safety factor for stress intensity factor (-)

  real(dp) :: paris_coefficient
  !! Paris equation material coefficient (-)

  real(dp) :: paris_power_law
  !! Paris equation material power law (-)

  real(dp) :: walker_coefficient
  !! walker coefficent (-)

  real(dp) :: fracture_toughness
  !! fracture toughness (MPa m^1/2)
end module CS_fatigue_variables
