! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fw_module
  !! Module containing first wall model
  !! author: J Morris, CCFE, Culham Science Centre
  !! N/A
  !! This module contains the PROCESS first wall model
  !! PROCESS Engineering paper (M. Kovari et al.)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Modules to import
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  private
  public :: friction, heat_transfer, fw_thermal_conductivity

contains

  subroutine friction(reynolds, darcy_friction)
    !! Calculate Darcy friction factor, using Haaland equation
    !! author: M Kovari, CCFE, Culham Science Centre
    !! reynolds : input real : Reynolds number
    !! darcy_friction : output real : Darcy friction factor
    !! Darcy friction factor, using Haaland equation, an approximation to the
    !! implicit Colebrook–White equationGnielinski correlation.
    !! https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use fwbs_variables, only: roughness, afw

    implicit none

    real(dp), intent(in) :: reynolds
    real(dp), intent(out) :: darcy_friction

    ! Local variables !
    ! !!!!!!!!!!!!!!!!!!

    real(dp) :: bracket

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Bracketed term in Haaland equation
    bracket = (roughness/afw/3.7)**1.11d0 + 6.9/reynolds

    ! Calculate Darcy friction factor
    darcy_friction = (1.8d0 * log10(bracket))**(-2)

  end subroutine friction

  function heat_transfer(masflx, rhof, radius, cf, viscf, kf)
    !! Calculate heat transfer coefficient using Gnielinski correlation
    !! author: M Kovari, CCFE, Culham Science Centre
    !! masflx : input real : coolant mass flux in a single channel (kg/m2/s)
    !! rhof : input real : coolant density (average of inlet and outlet) (kg/m3)
    !! radius : input real : coolant pipe radius (m)
    !! cf : input real : coolant specific heat capacity (average of inlet and outlet) (J/K)
    !! viscf : input real : coolant viscosity (average of inlet and outlet) (Pa.s)
    !! kf : input real : thermal conductivity of coolant (average of inlet and outlet) (W/m.K)
    !! Gnielinski correlation. Ignore the distinction between wall and
    !! bulk temperatures. Valid for:3000 < Re < 5e6, 0.5 < Pr < 2000
    !! https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, report_error

    implicit none

    ! Arguments

    ! Function output: Heat transfer coefficient (W/m2K)
    real(dp) :: heat_transfer

    ! Coolant mass flux in a single channel (kg/m2/s)
    real(dp) :: masflx

    ! Coolant density (average of inlet and outlet) (kg/m3)
    real(dp) :: rhof

    ! Coolant pipe radius (m)
    real(dp) :: radius

    ! Coolant specific heat capacity (average of inlet and outlet) (J/K)
    real(dp) :: cf

    ! Coolant viscosity (average of inlet and outlet) (Pa.s)
    real(dp) :: viscf

    ! Thermal conductivity of coolant (average of inlet and outlet) (W/m.K)
    real(dp) :: kf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Local variables

    ! Calculate flow velocity (m/s)
    real(dp) :: velocity

    ! Reynolds number
    real(dp) :: reynolds

    ! Prandtl number
    real(dp) :: pr

    ! Darcy friction factor
    real(dp) :: f

    ! Nusselt number
    real(dp) :: nusselt

    ! Pipe diameter (m)
    real(dp) ::diameter

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calculate pipe diameter (m)
    diameter = 2*radius

    ! Calculate flow velocity (m/s)
    velocity = masflx / rhof

    ! Calculate Reynolds number
    reynolds = rhof * velocity * diameter / viscf

    ! Calculate Prandtl number
    pr = cf * viscf / kf

    ! Calculate Darcy friction factor, using Haaland equation
    call friction(reynolds, f)

    ! Calculate the Nusselt number
    nusselt = (f/8.0d0)*(reynolds-1000.0d0)*pr / (1+12.7*sqrt(f/8.0d0)*(pr**0.6667-1.0d0))

    ! Calculate the heat transfer coefficient (W/m2K)
    heat_transfer = nusselt * kf / (2.0d0*radius)

    ! Check that Reynolds number is in valid range for the Gnielinski correlation
    if ( ( reynolds <= 3000.0d0 ).or.( reynolds > 5.0d6 ) ) then
      fdiags(1) = reynolds
      call report_error(225)
    end if

    ! Check that Prandtl number is in valid range for the Gnielinski correlation
    if ( ( pr < 0.5d0 ).or.( pr > 2000.0d0) ) then
      fdiags(1) = pr
      call report_error(226)
    end if

    ! Check that the Darcy friction factor is in valid range for the Gnielinski correlation
    if ( f <= 0.0d0 ) call report_error(227)

  end function heat_transfer

  function fw_thermal_conductivity(t)
    !! Calculates the thermal conductivity of the first wall
    !! t : input real : property temperature (K)
    !! Calculates the thermal conductivity of Eurofer (W/m/K).
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use fwbs_variables, only: fw_th_conductivity

    implicit none

    ! Function return value: thermal conductivity of first wall (W/m.K)
    real(dp) :: fw_thermal_conductivity

    ! Arguments !
    ! !!!!!!!!!!!!

    real(dp), intent(in) :: t

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Eurofer correlation, from "Fusion Demo Interim Structural Design Criteria -
    ! Appendix A Material Design Limit Data", F. Tavassoli, TW4-TTMS-005-D01, 2004
    ! t in Kelvin
    fw_thermal_conductivity = (5.4308D0 + 0.13565D0*t - 0.00023862D0*t*t + &
      1.3393D-7*t*t*t)*fw_th_conductivity/28.34D0

  end function fw_thermal_conductivity

end module
