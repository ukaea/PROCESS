module resistive_materials
  !! author: M. Kovari
  !!
  !! Variables relating to resistive materials in superconducting conductors

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  type resistive_material
    real(dp) :: cp
    !! Specific heat capacity J/(K kg)

    real(dp) :: rrr
    !! Residual resistivity ratio

    real(dp) :: resistivity
    !! Resistivity [ohm.m]

    real(dp) :: density
    !! kg/m3

    real(dp) :: cp_density
    !! Cp x density J/K/m3
  end type resistive_material

  type supercon_strand
      real(dp) :: area
      !! Superconducting strand area [m2]

      real(dp) :: critical_current
      !! Superconducting strand critical current [A]
  end type supercon_strand

  !#TODO: variables need descriptions
  type volume_fractions
    real(dp) :: copper_area,  copper_fraction
    real(dp) :: copper_bar_area
    real(dp) :: hastelloy_area, hastelloy_fraction
    real(dp) :: helium_area, helium_fraction
    real(dp) :: solder_area, solder_fraction
    real(dp) :: jacket_area, jacket_fraction
    real(dp) :: rebco_area,  rebco_fraction
    real(dp) :: critical_current
    real(dp) :: acs
    !! Area of cable space inside jacket
    real(dp) :: area
  end type volume_fractions
end module resistive_materials
