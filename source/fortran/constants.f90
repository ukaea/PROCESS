module constants
  !! author: J. Morris (UAKEA)
  !!
  !! Module containing miscellaneous numerical and physical constants
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  ! File output indexes
  integer, parameter :: iotty    = 6
  !! Standard output unit identifier

  integer, parameter :: nout     = 11
  !! Output file unit identifier

  integer, parameter :: nplot    = 12
  !! Plot data file unit identifier

  integer, parameter :: mfile    = 13
  !! Machine-optimised output file unit

  integer, parameter :: vfile    = 14
  !! Verbose diagnostics file

  integer, parameter :: opt_file = 15
  !! Optimisation information output file number

  integer, parameter :: sig_file = 16
  !! TF inboard stress radial distributions file number

  real(dp), parameter :: degrad = 0.01745329251D0
  !! degrees to radians, = pi/180

  real(dp), parameter :: echarge = 1.60217733D-19
  !! electron charge [C]

  real(dp), parameter :: emass = 9.10938370D-31
  !! electron mass [kg]

  real(dp), parameter :: mproton = 1.6726231D-27
  !! proton mass [kg]

  real(dp), parameter :: pi = 3.1415926535897932D0
  !! pi

  real(dp), parameter :: rmu0 = 1.256637062D-6
  !! permeability of free space  [H/m]

  real(dp), parameter :: twopi = 6.2831853071795862D0
  !! 2 pi

  real(dp), parameter :: umass = 1.660538921D-27
  !! unified atomic mass unit [kg

  real(dp), parameter :: epsilon0 = 8.85418781D-12
  !! permittivity of free space [Farad/m]

  real(dp), parameter :: cph2o = 4180.0D0
  !! specific heat capacity of water (J/kg/K)

  real(dp) :: dcopper
  !! density of copper (kg/m3)

  real(dp) :: dalu
  !! density of aluminium (kg/m3)

  real(dp), parameter :: denh2o = 985.0D0
  !! density of water (kg/m3)

  real(dp), parameter :: k_copper = 330.0D0
  !! Copper thermal conductivity (W/m/K)

  real(dp), parameter :: kh2o = 0.651D0
  !! thermal conductivity of water (W/m/K)

  real(dp), parameter :: muh2o = 4.71D-4
  !! water dynamic viscosity (kg/m/s)

  real(dp), parameter :: n_day_year = 365.2425D0
  !! Average number of days in a year

  contains

  subroutine init_constants
    !! Initialise module variables
    implicit none

    dcopper = 8900.0D0
    dalu = 2700.0D0
  end subroutine init_constants
end module constants
