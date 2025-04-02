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

  real(dp) :: beta_div
  !! field line angle wrt divertor target plate (degrees)

  real(dp) :: betai
  !! poloidal plane angle between divertor plate and leg, inboard (rad)

  real(dp) :: betao
  !! poloidal plane angle between divertor plate and leg, outboard (rad)

  real(dp) :: divclfr
  !! divertor coolant fraction

  real(dp) :: divdens
  !! divertor structure density (kg/m3)

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

  real(dp) :: prn1
  !! n-scrape-off / n-average plasma; (input for `ipedestal=0`, = nesep/dene if `ipedestal>=1`)


  real(dp) :: rconl
  !! connection length ratio, outboard side

  real(dp) :: rsrd
  !! effective separatrix/divertor radius ratio

  real(dp) :: tconl
  !! main plasma connection length (m)

  real(dp) :: tdiv
  !! temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)

  real(dp) :: xpertin
  !! perpendicular heat transport coefficient (m2/s)

  contains

  subroutine init_divertor_variables
    !! Initialise divertor_variables
    implicit none

    anginc = 0.262D0
    beta_div = 1.0D0
    betai = 1.0D0
    betao = 1.0D0
    divclfr = 0.3D0
    divdens = 1.0D4
    dz_divertor = 0.2D0
    divmas = 0.0D0
    divplt = 0.035D0
    divsur = 0.0D0
    fdiva = 1.11D0
    flux_exp = 2.0D0
    hldiv = 0.0D0
    i_hldiv = 2
    hldivlim = 5.0D0
    prn1 = 0.285D0
    rconl = 0.0D0
    rsrd = 0.0D0
    tconl = 0.0D0
    tdiv = 2.0D0
    xpertin = 2.0D0
  end subroutine init_divertor_variables
end module divertor_variables
