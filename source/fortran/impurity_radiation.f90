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

  private
  public :: Zav_of_te, init_impurity_radiation_module, element2index
  public :: impurity_arr_amass, impurity_arr_frac, impurity_arr_Label, &
   impurity_arr_len_tab, impurity_arr_Lz_Wm3, impurity_arr_Temp_keV, &
   impurity_arr_Z, impurity_arr_Zav


  !! (It is recommended to turn on
  !! constraint eqn.17 with iteration variable 28: fradpwr.)

  integer, public, parameter :: nimp = 14
  !! nimp /14/ FIX : number of ion species in impurity radiation model

  real(dp), public :: coreradius
  !! coreradius /0.6/ : normalised radius defining the 'core' region

  real(dp), public :: coreradiationfraction
  !! coreradiationfraction /1.0/ : fraction of radiation from 'core' region that is subtracted from the loss power

  !! fimp(nimp) /1.0,0.1,0.02,0.0,0.0,0.0,0.0,0.0,0.0016,0.0,0.0,0.0,0.0,0.0/ :
  !!        impurity number density fractions relative to electron density
  !!
  real(dp), public, dimension(nimp) :: fimp

  character*2, public, dimension(nimp) :: imp_label
  !! imp_label(nimp) : impurity ion species names:<UL>
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

  character*2, dimension(nimp) :: impurity_arr_Label
  integer, dimension(nimp) :: impurity_arr_Z
  real(dp), dimension(nimp) :: impurity_arr_amass
  real(dp), dimension(nimp) :: impurity_arr_frac
  integer, dimension(nimp) :: impurity_arr_len_tab
  real(dp), dimension(nimp, all_array_hotfix_len) :: impurity_arr_Temp_keV
  real(dp), dimension(nimp, all_array_hotfix_len) :: impurity_arr_Lz_Wm3
  real(dp), dimension(nimp, all_array_hotfix_len) :: impurity_arr_Zav


!   type(imp_dat),  dimension(nimp), save, public :: impurity_arr

  logical, public :: toolow
  !! Used for reporting error in function pimpden

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_impurity_radiation_module
    !! Initialise module variables
    implicit none

    coreradius = 0.6D0
    coreradiationfraction = 1.0D0
    fimp = (/ 1.0D0, 0.1D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
      0.0D0, 0.00D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
    imp_label = (/ &
      'H_', &
      'He', &
      'Be', &
      'C_', &
      'N_', &
      'O_', &
      'Ne', &
      'Si', &
      'Ar', &
      'Fe', &
      'Ni', &
      'Kr', &
      'Xe', &
      'W_'/)
      toolow = .false.
      impurity_arr_Label = "  "
      impurity_arr_Z = 0
      impurity_arr_amass = 0.0D0
      impurity_arr_len_tab = 0.0D0
      impurity_arr_Temp_keV = 0.0D0
      impurity_arr_Lz_Wm3 = 0.0D0
      impurity_arr_Zav = 0.0D0
      ! Re-initialise entire array
  end subroutine init_impurity_radiation_module

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function element2index(element_label)

    !! Returns the index of the element in the impurity array with
    !! a given name
    !! author: P J Knight, CCFE, Culham Science Centre
    !! element_label : input string : impurity name
    !! This function returns the index of the element
    !! in the impurity array with the corresponding name.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: report_error
    implicit none

    integer :: element2index

    !  Arguments

    character(len=*), intent(in) :: element_label

    !  Local variables

    integer :: i

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, size(impurity_arr_Label)
       if (element_label == impurity_arr_Label(i)) then
          element2index = i
          return
       end if
    end do

    !  Should only get here if there is a problem

    call report_error(34)

  end function element2index

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function Zav_of_te(imp_element_index,te)

    !! Electron temperature dependent average atomic number
    !! author: H Lux, CCFE, Culham Science Centre
    !! imp_element : input imp_dat : impurity element
    !! te : input real : electron temperature (keV)
    !! This routine returns the interpolated average atomic
    !! charge for a given electron temperature.
    !! <P>The Zav versus temperature data is interpolated from
    !! lookup tables from the ADAS data base provided by Martin O'Mullane.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none

    real(dp) :: Zav_of_te

    !  Arguments

    integer,  intent(in) :: imp_element_index
    real(dp), intent(in) :: te

    !  Local variables

    integer :: i
    real(dp) :: xi, yi, c

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Interpolate tabulated data

    !  Temperatures lower than the minimum, or higher than the maximum,
    !  are dealt with by taking the Zav value at the nearest tabulated
    !  temperature point

    if (te <= impurity_arr_Temp_keV(imp_element_index, 1)) then
       ! This should not be too unreasonable.
       Zav_of_te = impurity_arr_Zav(imp_element_index,1)

    else if (te >= impurity_arr_Temp_keV(imp_element_index, impurity_arr_len_tab(imp_element_index))) then
       !  This should be okay, as most elements are fully ionised by now.
       Zav_of_te = impurity_arr_Zav(imp_element_index, impurity_arr_len_tab(imp_element_index))

    else

       do i = 1, impurity_arr_len_tab(imp_element_index)-1

          !  Linear interpolation in log-lin space

          if ( (te > impurity_arr_Temp_keV(imp_element_index, i)) .and. &
               (te <= impurity_arr_Temp_keV(imp_element_index, i+1)) ) then

             yi = impurity_arr_Zav(imp_element_index, i)
             xi = log(impurity_arr_Temp_keV(imp_element_index, i))
             c  = (impurity_arr_Zav(imp_element_index, i+1) - yi) / &
                  (log(impurity_arr_Temp_keV(imp_element_index,i+1)) - xi)
             Zav_of_te = yi + c * (log(te) - xi)
             exit
          end if

       end do

    end if

  end function Zav_of_te


end module impurity_radiation_module
