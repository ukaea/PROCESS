! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module impurity_radiation_module

  !! Module for new impurity radiation calculations
  !! author: H Lux, CCFE, Culham Science Centre
  !! author: R Kemp, CCFE, Culham Science Centre
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! bremsstrahlung and line radiation of impurities
  !! including H  and He, assuming a coronal equilibrium.
  !! <P>The model is only valid for T &gt; 30 eV. For some impurity
  !! species there is also an upper temperature limit of T &lt; 40 keV.
  !! Johner, Fusion Science and Technology 59 (2011), pp 308-349
  !! Sertoli, private communication
  !! Kallenbach et al., Plasma Phys. Control. Fus. 55(2013) 124041
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  !! (It is recommended to turn on
  !! constraint eqn.17 with iteration variable 28: fradpwr.)

  integer, parameter :: n_impurities = 14
  !! n_impurities /14/ FIX : number of ion species in impurity radiation model

  real(dp) :: radius_plasma_core_norm
  !! radius_plasma_core_norm /0.6/ : normalised radius defining the 'core' region

  real(dp) :: coreradiationfraction
  !! coreradiationfraction /1.0/ : fraction of radiation from 'core' region that is subtracted from the loss power

  !! fimp(n_impurities) /1.0,0.1,0.02,0.0,0.0,0.0,0.0,0.0,0.0016,0.0,0.0,0.0,0.0,0.0/ :
  !!        impurity number density fractions relative to electron density
  !!
  real(dp), dimension(n_impurities) :: fimp

  character*2, dimension(n_impurities) :: imp_label
  !! imp_label(n_impurities) : impurity ion species names:<UL>
  !! <LI> ( 1)  Hydrogen  (fraction calculated by code)
  !! <LI> ( 2)  Helium
  !! <LI> ( 3)  Beryllium
  !! <LI> ( 4)  Carbon
  !! <LI> ( 5)  Nitrogen
  !! <LI> ( 6)  Oxygen
  !! <LI> ( 7)  Neon
  !! <LI> ( 8)  Silicon
  !! <LI> ( 9)  Argon
  !! <LI> (10)  Iron
  !! <LI> (11)  Nickel
  !! <LI> (12)  Krypton
  !! <LI> (13)  Xenon
  !! <LI> (14)  Tungsten</UL>

  !  Declare impurity data type

!   type :: imp_dat

!      character(len=2)  :: Label    !  Element name
!      integer           :: Z        !  Charge number
!      real(dp) :: amass    !  Atomic mass
!      real(dp) :: frac     !  Number density fraction (relative to ne)
!      integer           :: len_tab  !  Length of temperature vs. Lz table
!      !  Table of temperature values
!      real(dp), allocatable, dimension(:) :: Temp_keV
!      !  Table of corresponding Lz values
!      real(dp), allocatable, dimension(:) :: Lz_Wm3
!      !  Table of corresponding average atomic charge values
!      real(dp), allocatable, dimension(:) :: Zav

!   end type imp_dat

  ! derived type imp_dat (and hence impurity_arr) were
  ! incompatible with f2py and have been replaced with
  ! a less moder, but supported way of achieveing the
  ! same results

  integer, parameter :: all_array_hotfix_len = 200
  ! maximum length of the second dimensions of
  ! Temp_keV, Lz_Wm3, Zav
  ! since these can no longer be allocatable

  character*2, dimension(n_impurities) :: impurity_arr_Label
  integer, dimension(n_impurities) :: impurity_arr_Z
  real(dp), dimension(n_impurities) :: impurity_arr_amass
  real(dp), dimension(n_impurities) :: impurity_arr_frac
  integer, dimension(n_impurities) :: impurity_arr_len_tab
  real(dp), dimension(n_impurities, all_array_hotfix_len) :: impurity_arr_Temp_keV
  real(dp), dimension(n_impurities, all_array_hotfix_len) :: impurity_arr_Lz_Wm3
  real(dp), dimension(n_impurities, all_array_hotfix_len) :: impurity_arr_Zav

end module impurity_radiation_module
