module primary_pumping_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to `the i_coolant_pumping=3` option.
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
  !! ratio of specific heats for helium (`i_coolant_pumping=3`)

  ! if cp_he is required place here specific heat capacity at constant pressure: helium (`i_coolant_pumping=3`) [J/(kg.K)]
  ! cp_he is only used in private routine hcll.f90 at the moment

  real(dp) :: t_in_bb
  !! temperature in FW and blanket coolant at blanket entrance (`i_coolant_pumping=3`) [K]

  real(dp) :: t_out_bb
  !! temperature in FW and blanket coolant at blanket exit (`i_coolant_pumping=3`) [K]

  real(dp) :: p_he
  !! pressure in FW and blanket coolant at pump exit (`i_coolant_pumping=3`) [Pa]

  real(dp) :: dp_he
  !! pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_coolant_pumping=3`) [Pa]

  real(dp) :: dp_fw_blkt
  !! pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_coolant_pumping=3`) [Pa]

  real(dp) :: dp_fw
  !! pressure drop in FW coolant including heat exchanger and pipes (`i_coolant_pumping=3`) [Pa]

  real(dp) :: dp_blkt
  !! pressure drop in blanket coolant including heat exchanger and pipes (`i_coolant_pumping=3`) [Pa]

  real(dp) :: dp_liq
  !! pressure drop in liquid metal blanket coolant including heat exchanger and pipes (`i_coolant_pumping=3`) [Pa]

  real(dp) :: p_fw_blkt_coolant_pump_mw
  !! mechanical pumping power for FW and blanket including heat exchanger and
  !! pipes (`i_coolant_pumping=3`) [MW]

  real(dp) :: f_p_fw_blkt_pump
  !! Pumping power for FW and Blanket multiplier factor
end module primary_pumping_variables
