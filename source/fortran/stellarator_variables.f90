module stellarator_variables
  !! author: S. Muldrew (UKAEA), F. Warmer, J. Lion (IPP Greifswald)
  !!
  !! Module containing global variables relating to the stellarator model
  !!
  !!### References
  !!
  !! - A general stellarator version of the Systems code PROCESS, J. Lion et al, Nucl. Fus. (2021), https://doi.org/10.1088/1741-4326/ac2dbf
  !! - Stellarator Divertor Model for the Systems Code PROCESS, F. Warmer, 21/06/2013

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  integer :: istell
  !! Switch for stellarator option (set via `device.dat`):
  !!
  !! - =0 use tokamak model
  !! - =1 use stellarator model: Helias5
  !! - =2 use stellarator model: Helias4
  !! - =3 use stellarator model: Helias3
  !! - =4 use stellarator model: Wendelstein 7-X with 50 Coils
  !! - =5 use stellarator model: Wendelstein 7-X with 30 Coils
  !! - =6 use stellarator model: Use stella_conf.json file (any modulear stellarator, see documentation)

  real(dp) :: bmn
  !! relative radial field perturbation

  real(dp) :: f_asym
  !! divertor heat load peaking factor

  real(dp) :: f_rad
  !! radiated power fraction in SOL

  real(dp) :: f_w
  !! island size fraction factor

  real(dp) :: fdivwet
  !! wetted fraction of the divertor area

  real(dp) :: flpitch
  !! field line pitch (rad)

  real(dp) :: hportamax
  !! maximum available area for horizontal ports (m2)

  real(dp) :: hportpmax
  !! maximum available poloidal extent for horizontal ports (m)

  real(dp) :: hporttmax
  !! maximum available toroidal extent for horizontal ports (m)

  real(dp) :: iotabar
  !! rotational transform (reciprocal of tokamak q) for stellarator confinement time scaling laws

  integer :: isthtr
  !! Switch for stellarator auxiliary heating method:
  !!
  !! - = 1electron cyclotron resonance heating
  !! - = 2lower hybrid heating
  !! - = 3neutral beam injection

  integer :: m_res
  !! poloidal resonance number (1)

  real(dp) :: max_gyrotron_frequency
  !! Maximal available gyrotron frequency (input parameter) (Hz)

  integer :: n_res
  !! toroidal resonance number (1)

  real(dp) :: shear
  !! magnetic shear, derivative of iotabar (1)

  real(dp) :: te0_ecrh_achievable
  !! maximal central electron temperature as achievable by the ECRH, input. (keV)

  real(dp) :: vportamax
  !! maximum available area for vertical ports (m2)

  real(dp) :: vportpmax
  !! maximum available poloidal extent for vertical ports (m)

  real(dp) :: vporttmax
  !! maximum available toroidal extent for vertical ports (m)

  real(dp) :: powerht_constraint
  real(dp) :: powerscaling_constraint

  real(dp) :: f_st_coil_aspect
  !! aspect ratio for coils

end module stellarator_variables
