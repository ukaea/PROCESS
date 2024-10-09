module constants
  !! author: J. Morris (UAKEA)
  !!
  !! Module containing miscellaneous numerical and physical constants
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  ! Standard output unit identifier
  integer, parameter :: iotty    = 6

  ! Output file unit identifier
  integer, parameter :: nout     = 11

  ! Plot data file unit identifier
  integer, parameter :: nplot    = 12

  ! Machine-optimised output file unit
  integer, parameter :: mfile    = 13

  ! Verbose diagnostics file
  integer, parameter :: vfile    = 14

  ! Optimisation information output file number
  integer, parameter :: opt_file = 15

  ! TF inboard stress radial distributions file number
  integer, parameter :: sig_file = 16

  ! degrees to radians, = pi/180
  real(dp), parameter :: degrad = 0.01745329251D0

  ! electron charge [C]
  real(dp), parameter :: echarge = 1.60217733D-19

  ! Electron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=ELECTRON+MASS
  real(dp), parameter :: ELECTRON_MASS = 9.1093837139D-31

  ! Proton mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mp|search_for=PROTON+MASS
  real(dp), parameter :: PROTON_MASS = 1.67262192595D-27

  ! Deuteron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?md|search_for=DEUTERON+MASS
  real(dp), parameter :: DEUTERON_MASS = 3.3435837768D-27

  ! Triton mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mt|search_for=TRITON+MASS
  real(dp), parameter :: TRITON_MASS = 5.0073567512D-27

  ! Neutron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mn|search_for=NEUTRON+MASS
  real(dp), parameter :: NEUTRON_MASS = 1.67492750056D-27

  ! Alpha particle mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mal|search_for=alpha+mass
  real(dp), parameter :: ALPHA_MASS = 6.6446573450D-27

  ! Helion (3He) mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mh|search_for=HELION
  real(dp), parameter :: HELION_MASS = 5.0064127862D-27


  ! pi
  real(dp), parameter :: pi = 3.1415926535897932D0

  ! permeability of free space  [H/m]
  real(dp), parameter :: rmu0 = 1.256637062D-6

  ! 2 pi
  real(dp), parameter :: twopi = 6.2831853071795862D0

  ! unified atomic mass unit [kg]
  real(dp), parameter :: umass = 1.660538921D-27

  ! permittivity of free space [Farad/m]
  real(dp), parameter :: epsilon0 = 8.85418781D-12

  ! specific heat capacity of water (J/kg/K)
  real(dp), parameter :: cph2o = 4180.0D0

  ! density of copper (kg/m3)
  real(dp) :: dcopper

  ! density of aluminium (kg/m3)
  real(dp) :: dalu

  ! density of water (kg/m3)
  real(dp), parameter :: denh2o = 985.0D0

  ! Copper thermal conductivity (W/m/K)
  real(dp), parameter :: k_copper = 330.0D0

  ! thermal conductivity of water (W/m/K)
  real(dp), parameter :: kh2o = 0.651D0

  ! water dynamic viscosity (kg/m/s)
  real(dp), parameter :: muh2o = 4.71D-4

  ! Average number of days in a year
  real(dp), parameter :: n_day_year = 365.2425D0

  contains

  subroutine init_constants
    !! Initialise module variables
    implicit none

    dcopper = 8900.0D0
    dalu = 2700.0D0
  end subroutine init_constants
end module constants
