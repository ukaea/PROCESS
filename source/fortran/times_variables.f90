module times_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the plasma pulse timings
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: pulsetimings
  !! Switch for pulse timings (if lpulse=1):
  !!
  !! - =0, tohs = Ip(MA)/0.1 tramp, tqnch = input
  !! - =1, tohs = iteration var or input. tramp/tqnch max of input or tohs

  real(dp) :: tburn
  !! burn time (s) (calculated if `lpulse=1`)

  real(dp) :: tburn0
  !! burn time (s) - used for internal consistency

  real(dp) :: tcycle
  !! full cycle time (s)

  real(dp) :: tdown
  !! down time (s)

  real(dp) :: tdwell
  !! time between pulses in a pulsed reactor (s) (`iteration variable 17`)

  real(dp) :: t_fusion_ramp
  !! heating time, after current ramp up (s)

  real(dp), dimension(6) :: tim
  !! array of time points during plasma pulse (s)

  character*11, dimension(6) :: timelabel
  !! array of time labels during plasma pulse (s)

  character*11, dimension(5) :: intervallabel
  !! time intervals - as strings (s)

  real(dp) :: tohs
  !! plasma current ramp-up time for current initiation (s) (calculated if `lpulse=0`)
  !! (`iteration variable 65`)

  real(dp) :: tohsin
  !! Switch for plasma current ramp-up time (if lpulse=0):
  !!
  !! - = 0, tohs = tramp = tqnch = Ip(MA)/0.5
  !! - <>0, tohs = tohsin; tramp, tqnch are input

  real(dp) :: tpulse
  !! pulse length = tohs + t_fusion_ramp + tburn + tqnch

  real(dp) :: tqnch
  !! shut down time for PF coils (s); if pulsed, = tohs

  real(dp) :: tramp
  !! initial PF coil charge time (s); if pulsed, = tohs

  contains

  subroutine init_times_variables
    !! Initialise module variables
    implicit none

    pulsetimings = 1.0D0
    tburn = 1000.0D0
    tburn0 = 0.0D0
    tcycle = 0.0D0
    tdown = 0.0D0
    tdwell = 1800.0D0
    t_fusion_ramp = 10.0D0
    tim = 0.0D0
    timelabel = (/ 'Start', &
      'BOP  ', &
      'EOR  ', &
      'BOF  ', &
      'EOF  ', &
      'EOP  ' /)
    intervallabel = (/ 'tramp        ', &
      'tohs         ', &
      't_fusion_ramp', &
      'tburn        ', &
      'tqnch        ' /)
    tohs = 30.0D0
    tohsin = 0.0D0
    tpulse = 0.0D0
    tqnch = 15.0D0
    tramp = 15.0D0
  end subroutine init_times_variables
end module times_variables
