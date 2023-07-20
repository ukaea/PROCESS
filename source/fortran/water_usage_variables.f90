module water_usage_variables
  !! author: R Chapman (UKAEA)
  !!
  !! Module containing global variables relating to the water usage
  !!
  !!### References
  !!
  !! https://www.usgs.gov/special-topic/water-science-school/science/water-density
  !! https://www.thermal-engineering.org/what-is-latent-heat-of-vaporization-definition/
  !!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: airtemp
  !! ambient air temperature (degrees Celsius)

  real(dp) :: watertemp
  !! water temperature (degrees Celsius)

  real(dp) :: windspeed
  !! wind speed (m/s)

  real(dp) :: waterdens
  !! density of water (kg/m3)
  !!   for simplicity, set to static value applicable to water at 21 degC

  real(dp) :: latentheat
  !! latent heat of vaporization (J/kg)
  !!   for simplicity, set to static value applicable at 1 atm (100 kPa) air pressure

  real(dp) :: volheat
  !! volumetric heat of vaporization (J/m3)

  real(dp) :: evapratio
  !! evaporation ratio: ratio of the heat used to evaporate water
  !!   to the total heat discharged through the tower

  real(dp) :: evapvol
  !! evaporated volume of water (m3)

  real(dp) :: energypervol
  !! input waste (heat) energy cooled per evaporated volume (J/m3)

  real(dp) :: volperenergy
  !! volume evaporated by units of heat energy (m3/MJ)

  real(dp) :: waterusetower
  !! total volume of water used in cooling tower (m3)

  real(dp) :: wateruserecirc
  !! total volume of water used in recirculating system (m3)

  real(dp) :: wateruseonethru
  !! total volume of water used in once-through system (m3)

  contains

  subroutine init_watuse_variables
    !! Initialise module variables
    implicit none

    airtemp = 15.0D0
    watertemp = 5.0D0
    windspeed = 4.0D0
    waterdens = 998.02D0
    latentheat = 2257000.0D0
    volheat = 0.0D0
    evapratio = 0.0D0
    evapvol = 0.0D0
    energypervol = 0.0D0
    volperenergy = 0.0D0
    waterusetower = 0.0D0
    wateruserecirc = 0.0D0
    wateruseonethru = 0.0D0

  end subroutine init_watuse_variables

end module water_usage_variables
