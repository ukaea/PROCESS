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

  ! Electron / elementary charge [C]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?e|search_for=electron+charge
  real(dp), parameter :: electron_charge = 1.602176634D-19

  ! While the electron charge is a fundamental constant, the electron volt is a derived unit and 
  ! is added here for convenience. This allows better syntax and is more readable than using the electron
  ! charge constant directly when working with units of energy.

  ! Electron volt [J]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?evj|search_for=electron+volt
  real(dp), parameter :: electron_volt = 1.602176634D-19

  ! Electron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=ELECTRON+MASS
  real(dp), parameter :: electron_mass = 9.1093837139D-31

  ! Proton mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mp|search_for=PROTON+MASS
  real(dp), parameter :: proton_mass = 1.67262192595D-27

  ! Deuteron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?md|search_for=DEUTERON+MASS
  real(dp), parameter :: deuteron_mass = 3.3435837768D-27

  ! Triton mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mt|search_for=TRITON+MASS
  real(dp), parameter :: triton_mass = 5.0073567512D-27

  ! Neutron mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mn|search_for=NEUTRON+MASS
  real(dp), parameter :: neutron_mass = 1.67492750056D-27

  ! Alpha particle mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mal|search_for=alpha+mass
  real(dp), parameter :: alpha_mass = 6.6446573450D-27

  ! Helion (3He) mass [kg]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?mh|search_for=HELION
  real(dp), parameter :: helion_mass = 5.0064127862D-27

  ! Speed of light in vacuum (c) [m/s]
  ! Reference: National Institute of Standards and Technology (NIST)
  !            https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=light
  real(dp), parameter :: speed_light = 299792458D0

  ! Deuterium - Tritium reaction energy [J]
  ! Find the mass difference in the reactancts and products of the D-T reaction
  ! Multiply by the speed of light squared to get the energy released
  real(dp), parameter :: d_t_energy = (((deuteron_mass+triton_mass)-(alpha_mass+neutron_mass))*speed_light**2)

  ! Deuterium - Helion (3He) reaction energy [J]
  ! Find the mass difference in the reactancts and products of the D-3He reaction
  ! Multiply by the speed of light squared to get the energy released
  real(dp), parameter :: d_helium_energy = (((deuteron_mass+helion_mass)-(alpha_mass+proton_mass))*speed_light**2)

  ! Deuterium - Deuterium (3He producing) reaction energy [J]
  ! Find the mass difference in the reactancts and products of the D-D reaction
  ! Multiply by the speed of light squared to get the energy released
  real(dp), parameter :: dd_helium_energy = (((deuteron_mass+deuteron_mass)-(helion_mass+neutron_mass))*speed_light**2)

  ! Deuterium - Deuterium (Triton producing) reaction energy [J]
  ! Find the mass difference in the reactancts and products of the D-D reaction
  ! Multiply by the speed of light squared to get the energy released
  real(dp), parameter :: dd_triton_energy = (((deuteron_mass+deuteron_mass)-(triton_mass+proton_mass))*speed_light**2)

  ! Deuterium - Tritium reaction energy fraction carried by neutron
  ! Assuming centre of mass frame as the momenta of the fusion products exceed 
  ! those of the fusion reagents by many orders of magnitude. Assumed to be non-relativistic.
  ! Roughly 79.867% of the energy is carried by the neutron
  real(dp), parameter :: dt_neutron_energy_fraction = (alpha_mass/(neutron_mass+alpha_mass))

  ! Deuterium - Deuterium (3He producing) reaction energy fraction carried by neutron
  ! Assuming centre of mass frame as the momenta of the fusion products exceed 
  ! those of the fusion reagents by many orders of magnitude. Assumed to be non-relativistic.
  ! Roughly 74.935% of the energy is carried by the neutron
  real(dp), parameter :: dd_neutron_energy_fraction = (helion_mass/(neutron_mass+helion_mass))

  ! Deuterium - Deuterium (Triton producing) reaction energy fraction carried by proton
  ! Assuming centre of mass frame as the momenta of the fusion products exceed 
  ! those of the fusion reagents by many orders of magnitude. Assumed to be non-relativistic.
  ! Roughly 74.960% of the energy is carried by the proton
  real(dp), parameter :: dd_proton_energy_fraction = (triton_mass/(proton_mass+triton_mass))

  ! Deuterium - Helion (3He) reaction energy fraction carried by proton
  ! Assuming centre of mass frame as the momenta of the fusion products exceed 
  ! those of the fusion reagents by many orders of magnitude. Assumed to be non-relativistic.
  ! Roughly 74.957% of the energy is carried by the proton
  real(dp), parameter :: dhelium_proton_energy_fraction = (helion_mass/(proton_mass+helion_mass))

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
