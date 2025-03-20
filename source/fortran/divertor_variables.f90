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

  real(dp) :: adas
  !! area divertor / area main plasma (along separatrix)

  real(dp) :: anginc
  !! angle of incidence of field line on plate (rad)

  real(dp) :: beta_div
  !! field line angle wrt divertor target plate (degrees)

  real(dp) :: betai
  !! poloidal plane angle between divertor plate and leg, inboard (rad)

  real(dp) :: betao
  !! poloidal plane angle between divertor plate and leg, outboard (rad)

  real(dp) :: bpsout
  !! reference B_p at outboard divertor strike point (T)

  real(dp) :: c1div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: c2div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: c3div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: c4div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: c5div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: c6div
  !! fitting coefficient to adjust ptpdiv, ppdiv

  real(dp) :: delld
  !! coeff for power distribution along main plasma

  real(dp) :: dendiv
  !! plasma density at divertor (10**20 /m3)

  real(dp) :: densin
  !! density at plate (on separatrix) (10**20 /m3)

  real(dp) :: divclfr
  !! divertor coolant fraction

  real(dp) :: divdens
  !! divertor structure density (kg/m3)

  integer :: divdum
  !! switch for divertor Zeff model:
  !!
  !! - =0 calc
  !! - =1 input
  !#TODO: switch name should be changed to i_<something>

  real(dp) :: dz_divertor
  !! divertor structure vertical thickness (m)

  real(dp) :: divmas
  !! divertor plate mass (kg)

  real(dp) :: divplt
  !! divertor plate thickness (m) (from Spears, Sept 1990)

  real(dp) :: divsur
  !! divertor surface area (m2)

  real(dp) :: fdfs
  !! radial gradient ratio

  real(dp) :: fdiva
  !! divertor area fudge factor (for ITER, Sept 1990)

  real(dp) :: fhout
  !! fraction of power to outboard divertor (for single null)

  real(dp) :: fififi
  !! coefficient for gamdiv
  !#TODO: what the hell is this variable name...

  real(dp) :: flux_exp
  !! The plasma flux expansion in the divertor (default 2; Wade 2020)

  real(dp) :: frrp
  !! fraction of radiated power to plate

  real(dp) :: hldiv
  !! divertor heat load (MW/m2)

  integer :: i_hldiv
  !! switch for user input hldiv:
  !!
  !! - = 0: divtart model turned off and user inputs hldiv
  !! - = 1: divtart model calculates hldiv
  !! - = 2: divwade model calculates hldiv

  real(dp) :: hldivlim
  !! heat load limit (MW/m2)

  real(dp) :: ksic
  !! power fraction for outboard double-null scrape-off plasma

  real(dp) :: lamp
  !! power flow width (m)

  real(dp) :: minstang
  !! minimum strike angle for heat flux calculation

  real(dp) :: omegan
  !! pressure ratio (nT)_plasma / (nT)_scrape-off

  real(dp) :: omlarg
  !! power spillage to private flux factor

  real(dp) :: ppdivr
  !! peak heat load at plate (with radiation) (MW/m2)

  real(dp) :: prn1
  !! n-scrape-off / n-average plasma; (input for `ipedestal=0`, = nesep/dene if `ipedestal>=1`)

  real(dp) :: ptpdiv
  !! peak temperature at the plate (eV)

  real(dp) :: rconl
  !! connection length ratio, outboard side

  real(dp) :: rlclolcn
  !! ratio of collision length / connection length

  real(dp) :: rlenmax
  !! maximum value for length ratio (rlclolcn) (`constraintg eqn 22`)

  real(dp) :: rsrd
  !! effective separatrix/divertor radius ratio

  real(dp) :: tconl
  !! main plasma connection length (m)

  real(dp) :: tdiv
  !! temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)

  real(dp) :: tsep
  !! temperature at the separatrix (eV)

  real(dp) :: xparain
  !! parallel heat transport coefficient (m2/s)

  real(dp) :: xpertin
  !! perpendicular heat transport coefficient (m2/s)

  real(dp) :: zeffdiv
  !! Zeff in the divertor region (if `divdum/=0`)

end module divertor_variables
