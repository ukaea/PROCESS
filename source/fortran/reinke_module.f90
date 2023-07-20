module reinke_module
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
   implicit none

  !private

  !  Module-level variables

  !integer ::
  real(dp) :: vcritx

contains

  subroutine init_reinke_module
    !! Initialise module variables
    implicit none

    vcritx = 0.0D0
  end subroutine init_reinke_module

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function reinke_tsep(bt, flh, qstar, rmajor, eps, fgw, kappa, lhat)

    !! Function for calculating upstream temperature(keV) in Reinke model
    !! author: H Lux, CCFE/UKAEA
    !! bt      : input real : toroidal field on axis (T)
    !! flh     : input real : fraction of Psep/P_LH
    !! qstar   : input real : safety factor similar to q95 (see #707)
    !! rmajor  : input real : major radius (m)
    !! eps     : input real : inverse aspect ratio
    !! fgw     : input real : ratio of volume averaged density to n_GW
    !! kappa   : input real : elongation
    !! lhat    : input real : connection length factor
    !! This function calculates the upstream temperature in the
    !! divertor/SoL model used for the Reinke citerion.
    !! Issue #707
    !! M.L. Reinke 2017 Nucl. Fusion 57 034004

    implicit none
    real(dp) :: reinke_tsep
    real(dp) :: bt, flh, qstar, rmajor, eps, fgw, kappa, lhat
    real(dp), parameter :: kappa_0 = 2D3 !Stangeby W/m/eV^(7/2)

    reinke_tsep = bt**0.72 * flh**0.29 * fgw**0.21 * qstar**0.08 * rmajor**0.33
    !reinke_tsep = bt**0.72 * flh**0.2857 * fgw**0.2057 * qstar**0.08 * rmajor**0.3314

    reinke_tsep = reinke_tsep * eps**0.15 * (1.D0 + kappa**2.)**0.34
    !reinke_tsep = reinke_tsep * eps**0.1486 * (1.D0 + kappa**2.)**0.34

    reinke_tsep = reinke_tsep * lhat**0.29 * kappa_0 **(-0.29) * 0.285D0
    !reinke_tsep = reinke_tsep * lhat**0.2857 * kappa_0 **(-0.2857) * 0.285D0

  end function reinke_tsep

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Disabled for ease of #1542 - Tim
  ! subroutine test_reinke()
  !   use impurity_radiation_module, only: imp_label
  !   implicit none

  !   real(dp) :: testResult_fZ_DEMOBASE, testResult_fZ_ASDEXBASE, testInput_tsep
  !   integer :: i, j
  !   real(dp), parameter :: test_Bt = 5.8547
  !   type(imp_dat),  dimension(14), save :: test_imp_arr
  !   real(dp), dimension(14) :: impurity_enrichment

  !   do i=1,14
  !      test_imp_arr(i)%frac = 0.0d0
  !      test_imp_arr(i)%label = imp_label(i)
  !      impurity_enrichment(i) = 1.0d0
  !   end do
  !   test_imp_arr(1)%frac = 1.0d0
  !   test_imp_arr(2)%frac = 0.1d0
  !   test_imp_arr(9)%frac = 0.001d0
  !   test_imp_arr(13)%frac = 4.4d-04
  !   test_imp_arr(14)%frac = 5d-05


  !   !testResult_tsep = reinke_tsep(test_Bt, 1.0d0, 3.0d0, 1.65d0, 0.33d0, 0.8d0, 1.7d0, 4.33d0)
  !   !                              bt, flh, qstar, rmajor, eps, fgw, kappa, lhat
  !   !write(*,*) 'reinke_tsep = ', testResult_tsep

  !  ! Open file to output the data for the test of the fzmin function.
  !   open(unit=1,file='FZMIN_TEST.DAT')

  !  ! We will output fzmin for a range of tesep values, using both a DEMO baseline and a ASDEX like parameters
  !  ! The ASDEX like fZ should be compariable with data in M.L. Reinke et al, 2017
  !   do j = 1,39
  !     testInput_tsep =0.02d0+0.02d0*j
  !     testResult_fZ_DEMOBASE = reinke_fzmin(test_Bt, 1.4114d0, 3.0d0, 9.0019d0, 0.3225d0, 0.5651d0, 0.8d0, 1.848d0, &
  !        4.33d0, 0.1d0, testInput_tsep, 9, test_imp_arr%frac, impurity_enrichment)
  !     testResult_fZ_ASDEXBASE = reinke_fzmin(test_Bt, 1.0d0, 3.0d0, 1.65d0, 0.33d0, 1.0d0, 0.8d0, 1.7d0, &
  !        4.33d0, 0.1d0, testInput_tsep, 9, test_imp_arr%frac, impurity_enrichment)
  !     write(1,*) testInput_tsep, testResult_fZ_DEMOBASE, testResult_fZ_ASDEXBASE, test_imp_arr(9)%frac
  !   end do

  !  close(unit=1)
  !  write(*,*) 'fz minium data written to FZMIN_TEST.DAT'

  !   !testResult_fZ = reinke_fzmin(test_Bt, 1.4114, 3.0d0, 9.0019d0, 0.3225d0, 1.0d0, 0.8d0, 1.848d0, &
  !   !      4.33d0, 0.1d0, testInput_tsep, 9, test_imp_arr%frac, impurity_enrichment)
  !                                 !bt, flh, qstar, rmajor, eps, fsep, fgw, kappa,
  !        !lhat, netau, impvardiv, impurity_arr%frac, impurity_enrichment
  !   !write(*,*) 'reinke_fzmin = ', testResult_fZ

  !   !if(testResult_fZ /= 1.4) then
  !     ! call report_error(217)
  !   !end if

  !   !if(testResult_tsep /= 1.4) then
  !      !call report_error(217)
  !   !end if

  ! end subroutine test_reinke



end module reinke_module
