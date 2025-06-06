module pf_power_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the PF coil power conversion system
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: acptmax
  !! average of currents in PF circuits (kA)

  real(dp) :: ensxpfm
  !! maximum stored energy in the PF circuits (MJ)

  integer :: i_pf_energy_storage_source
  !! Switch for PF coil energy storage option:
  !!
  !! - =1 all power from MGF (motor-generator flywheel) units
  !! - =2 all pulsed power from line
  !! - =3 PF power from MGF, heating from line
  !   (In fact, options 1 and 3 are not treated differently)

  real(dp) :: pfckts
  !! number of PF coil circuits

  real(dp) :: spfbusl
  !! total PF coil circuit bus length (m)

  real(dp) :: spsmva
  !! sum of PF power supply ratings (MVA)

  real(dp) :: srcktpm
  !! sum of resistive PF coil power (kW)

  real(dp) :: vpfskv
  !! PF coil voltage (kV)

  real(dp) :: peakpoloidalpower
  !! Peak absolute rate of change of stored energy in poloidal field (MW)

  real(dp) :: maxpoloidalpower
  !! Maximum permitted absolute rate of change of stored energy in poloidal field (MW)

  real(dp), dimension(5) :: poloidalpower
  !! Poloidal power usage at time t (MW)
end module pf_power_variables
