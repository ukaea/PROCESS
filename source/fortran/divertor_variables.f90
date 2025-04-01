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

  real(dp) :: fdiva
  !! divertor area fudge factor (for ITER, Sept 1990)

  real(dp) :: fhout
  !! fraction of power to outboard divertor (for single null)

  real(dp) :: flux_exp
  !! The plasma flux expansion in the divertor (default 2; Wade 2020)

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


  real(dp) :: minstang
  !! minimum strike angle for heat flux calculation

  real(dp) :: omlarg
  !! power spillage to private flux factor


  real(dp) :: prn1
  !! n-scrape-off / n-average plasma; (input for `ipedestal=0`, = nesep/dene if `ipedestal>=1`)


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


  real(dp) :: xparain
  !! parallel heat transport coefficient (m2/s)

  real(dp) :: xpertin
  !! perpendicular heat transport coefficient (m2/s)

  real(dp) :: zeffdiv
  !! Zeff in the divertor region (if `divdum/=0`)

  contains

  subroutine init_divertor_variables
    !! Initialise divertor_variables
    implicit none

    adas = 0.0D0
    anginc = 0.262D0
    beta_div = 1.0D0
    betai = 1.0D0
    betao = 1.0D0
    bpsout = 0.60D0
    divclfr = 0.3D0
    divdens = 1.0D4
    divdum = 0
    dz_divertor = 0.2D0
    divmas = 0.0D0
    divplt = 0.035D0
    divsur = 0.0D0
    fdiva = 1.11D0
    fhout = 0.0D0
    flux_exp = 2.0D0
    hldiv = 0.0D0
    i_hldiv = 0
    hldivlim = 5.0D0
    ksic = 0.8D0
    minstang = 0.0D0
    omlarg = 0.0D0
    prn1 = 0.285D0
    rconl = 0.0D0
    rlclolcn = 0.0D0
    rlenmax = 0.5D0
    rsrd = 0.0D0
    tconl = 0.0D0
    tdiv = 2.0D0
    xparain = 2.1D3
    xpertin = 2.0D0
    zeffdiv = 1.0D0
  end subroutine init_divertor_variables
end module divertor_variables
