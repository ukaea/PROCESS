module divertor_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the tokamak divertor components
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: anginc
  !! angle of incidence of field line on plate (rad)

  real(dp) :: deg_div_field_plate
  !! field line angle wrt divertor target plate (degrees)

  real(dp) :: betai
  !! poloidal plane angle between divertor plate and leg, inboard (rad)

  real(dp) :: betao
  !! poloidal plane angle between divertor plate and leg, outboard (rad)

  real(dp) :: f_vol_div_coolant
  !! divertor coolant fraction

  real(dp) :: den_div_structure
  !! divertor structure density (kg/m3)

  real(dp) :: dz_divertor
  !! divertor structure vertical thickness (m)

  real(dp) :: m_div_plate
  !! divertor plate mass (kg)

  real(dp) :: dx_div_plate
  !! divertor plate thickness (m) (from Spears, Sept 1990)

  real(dp) :: a_div_surface_total
  !! divertor surface area (m2)

  real(dp) :: fdiva
  !! divertor area fudge factor (for ITER, Sept 1990)

  real(dp) :: f_div_flux_expansion
  !! The plasma flux expansion in the divertor (default 2; Wade 2020)

  real(dp) :: pflux_div_heat_load_mw
  !! divertor heat load (MW/m2)

  integer :: i_div_heat_load
  !! switch for user input pflux_div_heat_load_mw:
  !!
  !! - = 0: divtart model turned off and user inputs pflux_div_heat_load_mw
  !! - = 1: divtart model calculates pflux_div_heat_load_mw
  !! - = 2: divwade model calculates pflux_div_heat_load_mw

  real(dp) :: pflux_div_heat_load_max_mw
  !! heat load limit (MW/m2)

  real(dp) :: prn1
  !! n-scrape-off / n-average plasma; (input for `ipedestal=0`, = nesep/dene if `ipedestal>=1`)


  real(dp) :: rsrd
  !! effective separatrix/divertor radius ratio

  real(dp) :: tconl
  !! main plasma connection length (m)

  real(dp) :: tdiv
  !! temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)

  real(dp) :: xpertin
  !! perpendicular heat transport coefficient (m2/s)

  real(dp) :: zeffdiv
  !! Zeff in the divertor region (if `divdum/=0`)

end module divertor_variables
