module primary_pumping_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to `the primary_pumping=3` option.
  !! (Mechanical pumping power is calculated using specified pressure drop)
  !!
  !!### References
  !!
  !! - issue #503

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: gamma_he
  !! ratio of specific heats for helium (`primary_pumping=3`)

  ! if cp_he is required place here specific heat capacity at constant pressure: helium (`primary_pumping=3`) [J/(kg.K)]
  ! cp_he is only used in private routine hcll.f90 at the moment

  real(dp) :: t_in_bb
  !! temperature in FW and blanket coolant at blanket entrance (`primary_pumping=3`) [K]

  real(dp) :: t_out_bb
  !! temperature in FW and blanket coolant at blanket exit (`primary_pumping=3`) [K]

  real(dp) :: p_he
  !! pressure in FW and blanket coolant at pump exit (`primary_pumping=3`) [Pa]

  real(dp) :: dp_he
  !! pressure drop in FW and blanket coolant including heat exchanger and pipes (`primary_pumping=3`) [Pa]

  real(dp) :: dp_fw_blkt
  !! pressure drop in FW and blanket coolant including heat exchanger and pipes (`primary_pumping=3`) [Pa]

  real(dp) :: dp_fw
  !! pressure drop in FW coolant including heat exchanger and pipes (`primary_pumping=3`) [Pa]

  real(dp) :: dp_blkt
  !! pressure drop in blanket coolant including heat exchanger and pipes (`primary_pumping=3`) [Pa]

  real(dp) :: dp_liq
  !! pressure drop in liquid metal blanket coolant including heat exchanger and pipes (`primary_pumping=3`) [Pa]

  real(dp) :: p_fw_blanket_pumping_mw
  !! mechanical pumping power for FW and blanket including heat exchanger and
  !! pipes (`primary_pumping=3`) [MW]

  contains

  subroutine init_primary_pumping_variables
    !! Initialise module variables
    implicit none

    !! initialise variables with default values in the absence of a value in the input file.
    gamma_he = 1.667D0 ! Ratio of specific heats  Helium
    !! cp_he = 5.195.0D3 ! Specific Heat J/kg K (disabled at the moment)
    t_in_bb = 573.13D0 ! K
    t_out_bb = 773.13D0 ! K
    p_he = 8.0D6 ! Pa
    dp_he = 5.5D5 ! Pa
    dp_fw_blkt = 1.5D5 ! Pa
    dp_fw = 1.5D5 ! Pa
    dp_blkt = 3.5D3 ! Pa
    dp_liq = 1.0D7 ! Pa
    p_fw_blanket_pumping_mw = 0.0D0

  end subroutine init_primary_pumping_variables

end module primary_pumping_variables
