module times_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the plasma pulse timings
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: pulsetimings
  !! Switch for pulse timings (if i_pulsed_plant=1):
  !!
  !! - =0, t_current_ramp_up = Ip(MA)/0.1 t_precharge, t_ramp_down = input
  !! - =1, t_current_ramp_up = iteration var or input. t_precharge/t_ramp_down max of input or t_current_ramp_up

  real(dp) :: t_burn
  !! flat-top duration (s) (calculated if `i_pulsed_plant=1`)

  real(dp) :: t_burn_0
  !! burn time (s) - used for internal consistency

  real(dp) :: t_cycle
  !! full cycle time (s)

  real(dp) :: tdown
  !! down time (s)

  real(dp) :: t_between_pulse
  !! time between pulses in a pulsed reactor (s) (`iteration variable 17`)

  real(dp) :: t_fusion_ramp
  !! time for plasma temperature and density rise to full values (s)

  real(dp), dimension(6) :: tim
  !! array of time points during plasma pulse (s)

  character*11, dimension(6) :: timelabel
  !! array of time labels during plasma pulse (s)

  character*19, dimension(5) :: intervallabel
  !! time intervals - as strings (s)

  real(dp) :: t_current_ramp_up
  !! time for plasma current to ramp up to approx. full value (s) (calculated if `i_pulsed_plant=0`)
  !! (`iteration variable 65`)

  integer :: i_t_current_ramp_up
  !! Switch for plasma current ramp-up time (if i_pulsed_plant=0):
  !!
  !! - = 0, t_current_ramp_up = t_precharge = t_ramp_down = Ip(MA)/0.5
  !! - = 1, t_current_ramp_up, t_precharge, t_ramp_down are input

  real(dp) :: t_pulse_repetition
  !! pulse length = t_current_ramp_up + t_fusion_ramp + t_burn + t_ramp_down

  real(dp) :: t_ramp_down
  !! time for plasma current, density, and temperature to ramp down to zero, simultaneously (s); if pulsed, = t_current_ramp_up
  !! the CS and PF coil currents also ramp to zero at the same time

  real(dp) :: t_precharge
  !! the time for the central solenoid and PF coils to ramp from zero to max current (s); if pulsed, = t_current_ramp_up
end module times_variables
