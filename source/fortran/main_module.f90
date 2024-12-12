module main_module

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

contains

subroutine herror(ifail)

  !! Routine to print out relevant messages in the case of an
  !! unfeasible result from a HYBRD (non-optimisation) run
  !! author: P J Knight, CCFE, Culham Science Centre
  !! ifail  : input integer : error flag
  !! This routine prints out relevant messages in the case of
  !! an unfeasible result from a HYBRD (non-optimisation) run.
  !! <P>The messages are written to units NOUT and IOTTY, which are
  !! by default the output file and screen, respectively.
  !! <P>If <CODE>IFAIL=1</CODE> then a feasible solution has been
  !! found and therefore no error message is required.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: nout, iotty
  use process_output, only: oblnkl, ocmmnt
  implicit none

  !  Arguments
  integer, intent(in) :: ifail

  !  Local variables

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case (ifail)

  case (:-1)
     call ocmmnt(nout, 'User-terminated execution of HYBRD.')
     call ocmmnt(iotty,'User-terminated execution of HYBRD.')

  case (0)
     call ocmmnt(nout, 'Improper input parameters to the HYBRD routine.')
     call ocmmnt(nout, 'PROCESS coding must be checked.')

     call ocmmnt(iotty,'Improper input parameters to the HYBRD routine.')
     call ocmmnt(iotty,'PROCESS coding must be checked.')

  case (1)
     continue

  case (2)
     call ocmmnt(nout,'The maximum number of calls has been reached without solution,')
     call ocmmnt(nout,'suggesting that the iteration is not making good progress.')
     call ocmmnt(nout,'Try changing the variables in IXC.')

     call ocmmnt(iotty,'The maximum number of calls has been reached without solution,')
     call ocmmnt(iotty,'suggesting that the iteration is not making good progress.')
     call ocmmnt(iotty,'Try changing the variables in IXC.')

  case (3)
     call ocmmnt(nout,'The tolerance is too small: No further improvement in the approximate solution is possible.')
     call ocmmnt(nout,'Try raising the value of FTOL.')

     call ocmmnt(iotty, 'The tolerance is too small: No further improvement in the approximate solution is possible.')
     call ocmmnt(iotty,'in the approximate solution is possible.')
     call ocmmnt(iotty,'Try raising the value of FTOL.')

  case (4)
     call ocmmnt(nout,'The iteration is not making good progress.')
     call ocmmnt(nout,'The code may be stuck in a minimum in the residual')
     call ocmmnt(nout,'space that is significantly above zero.')
     call oblnkl(nout)
     call ocmmnt(nout,'There is either no solution possible, or the code')
     call ocmmnt(nout,'is failing to escape from a deep local minimum.')
     call ocmmnt(nout,'Try changing the variables in IXC, or modify their initial values.')

     call ocmmnt(iotty,'The iteration is not making good progress.')
     call ocmmnt(iotty,'The code may be stuck in a minimum in the residual')
     call ocmmnt(iotty,'space that is significantly above zero.')
     call oblnkl(iotty)
     call ocmmnt(iotty,'There is either no solution possible, or the code')
     call ocmmnt(iotty,'is failing to escape from a deep local minimum.')
     call ocmmnt(iotty,'Try changing the variables in IXC, or modify their initial values.')

  case default
     call ocmmnt(nout,'This value of IFAIL should not be possible...')
     call ocmmnt(nout,'See source code for details.')

     call ocmmnt(iotty,'This value of IFAIL should not be possible...')
     call ocmmnt(iotty,'See source code for details.')

  end select

end subroutine herror

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine verror(ifail)

  !! Routine to print out relevant messages in the case of an
  !! unfeasible result from a VMCON (optimisation) run
  !! author: P J Knight, CCFE, Culham Science Centre
  !! ifail  : input integer : error flag
  !! This routine prints out relevant messages in the case of
  !! an unfeasible result from a VMCON (optimisation) run.
  !! <P>The messages are written to units NOUT and IOTTY, which are
  !! by default the output file and screen, respectively.
  !! <P>If <CODE>IFAIL=1</CODE> then a feasible solution has been
  !! found and therefore no error message is required.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: nout, iotty
  use process_output, only: ocmmnt, oblnkl
  implicit none

  !  Arguments
  integer, intent(in) :: ifail

  !  Local variables

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case (ifail)

  case (:-1)
     call ocmmnt(nout, 'User-terminated execution of VMCON.')
     call ocmmnt(iotty,'User-terminated execution of VMCON.')

  case (0)
     call ocmmnt(nout, 'Improper input parameters to the VMCON routine.')
     call ocmmnt(nout, 'PROCESS coding must be checked.')

     call ocmmnt(iotty,'Improper input parameters to the VMCON routine.')
     call ocmmnt(iotty,'PROCESS coding must be checked.')

  case (1)
     continue

  case (2)
     call ocmmnt(nout,'The maximum number of calls has been reached without solution.')
     call ocmmnt(nout,'The code may be stuck in a minimum in the residual space that is significantly above zero.')
     call oblnkl(nout)
     call ocmmnt(nout,'There is either no solution possible, or the code')
     call ocmmnt(nout,'is failing to escape from a deep local minimum.')
     call ocmmnt(nout,'Try changing the variables in IXC, or modify their initial values.')

     call ocmmnt(iotty,'The maximum number of calls has been reached without solution.')
     call ocmmnt(iotty,'The code may be stuck in a minimum in the residual space that is significantly above zero.')
     call oblnkl(nout)
     call oblnkl(iotty)
     call ocmmnt(iotty,'There is either no solution possible, or the code')
     call ocmmnt(iotty,'is failing to escape from a deep local minimum.')
     call ocmmnt(iotty,'Try changing the variables in IXC, or modify their initial values.')

  case (3)
     call ocmmnt(nout,'The line search required the maximum of 10 calls.')
     call ocmmnt(nout,'A feasible solution may be difficult to achieve.')
     call ocmmnt(nout,'Try changing or adding variables to IXC.')

     call ocmmnt(iotty,'The line search required the maximum of 10 calls.')
     call ocmmnt(iotty,'A feasible solution may be difficult to achieve.')
     call ocmmnt(iotty,'Try changing or adding variables to IXC.')

  case (4)
     call ocmmnt(nout,'An uphill search direction was found.')
     call ocmmnt(nout,'Try changing the equations in ICC, or')
     call ocmmnt(nout,'adding new variables to IXC.')

     call ocmmnt(iotty,'An uphill search direction was found.')
     call ocmmnt(iotty,'Try changing the equations in ICC, or')
     call ocmmnt(iotty,'adding new variables to IXC.')

  case (5)
     call ocmmnt(nout, &
          'The quadratic programming technique was unable to')
     call ocmmnt(nout,'find a feasible point.')
     call oblnkl(nout)
     call ocmmnt(nout,'Try changing or adding variables to IXC, or modify')
     call ocmmnt(nout,'their initial values (especially if only 1 optimisation')
     call ocmmnt(nout,'iteration was performed).')

     call ocmmnt(iotty, &
          'The quadratic programming technique was unable to')
     call ocmmnt(iotty,'find a feasible point.')
     call oblnkl(iotty)
     call ocmmnt(iotty,'Try changing or adding variables to IXC, or modify')
     call ocmmnt(iotty,'their initial values (especially if only 1 optimisation')
     call ocmmnt(iotty,'iteration was performed).')

  case (6)
     call ocmmnt(nout, &
          'The quadratic programming technique was restricted')
     call ocmmnt(nout, &
          'by an artificial bound, or failed due to a singular')
     call ocmmnt(nout,'matrix.')
     call ocmmnt(nout,'Try changing the equations in ICC, or')
     call ocmmnt(nout,'adding new variables to IXC.')

     call ocmmnt(iotty, &
          'The quadratic programming technique was restricted')
     call ocmmnt(iotty, &
          'by an artificial bound, or failed due to a singular')
     call ocmmnt(iotty,'matrix.')
     call ocmmnt(iotty,'Try changing the equations in ICC, or')
     call ocmmnt(iotty,'adding new variables to IXC.')

  case default
     call ocmmnt(nout,'This value of IFAIL should not be possible...')
     call ocmmnt(nout,'See source code for details.')

     call ocmmnt(iotty,'This value of IFAIL should not be possible...')
     call ocmmnt(iotty,'See source code for details.')

  end select

end subroutine verror

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine runtests
   ! These tests should gradually be moved to pytest
  use constants, only: nout
  use maths_library, only: binomial
  use process_output, only: ocmmnt, ovarre
!   use pfcoil_module, only: brookscoil
!   use reinke_module, only: test_reinke
  implicit none
  real(dp) :: fshift, xf, enpa,ftherm,fpp,cdeff, ampperwatt
  logical :: Temperature_capped
  call ovarre(nout,'Binomial coefficients C(5,0): 1', '(binomial(5,0))', binomial(5,0))
  call ovarre(nout,'Binomial coefficients C(5,1): 5', '(binomial(5,1))', binomial(5,1))
  call ovarre(nout,'Binomial coefficients C(5,2): 10', '(binomial(5,2))', binomial(5,2))
  call ovarre(nout,'Binomial coefficients C(5,3): 10', '(binomial(5,3))', binomial(5,3))
  call ovarre(nout,'Binomial coefficients C(5,4): 5', '(binomial(5,4))', binomial(5,4))
  call ovarre(nout,'Binomial coefficients C(5,5): 1', '(binomial(5,5))', binomial(5,5))

   !   call brookscoil(nout) Moved to pytest
  ! Disabled for ease of #1542 - Tim
!   call test_reinke()
end subroutine runtests


subroutine get_DDMonYYTimeZone(dt_time)
  !! Routine to get date, time and timezone
  !! author: M Kumar, CCFE, Culham Science Centre
  !! dt_time : output string  : String containing formatted time and date
  !! This routine calls the intrinsic DATE_AND_TIME subroutine
  !! and format the output in
  !! DD Mon YYYY hr:minute:second time difference from UTC.

! Arguments
    CHARACTER(len = *), INTENT(OUT) :: dt_time
! Local variables
    INTEGER :: values(8)
    CHARACTER(len = 1), parameter :: tspt = ":"
    CHARACTER(len = 1), parameter :: spspt = " "

    CHARACTER(len = 2)  :: dd
    CHARACTER(len = 5)  :: mons(12)
    CHARACTER(len = 4)  :: yyyy
    CHARACTER(len = 2)  :: hr    ! Hour of the day
    CHARACTER(len = 2)  :: mnt   ! Minute of the hour
    CHARACTER(len = 2)  :: scnd  ! The seconds of the minute
    CHARACTER(len = 5)  :: zn    ! In form (+-)hhmm, representing the difference with respect to Coordinated Universal Time (UTC).
    CHARACTER(len = 20) :: znfrmt

    mons = [' Jan ',' Feb ',' Mar ',' Apr ',' May ',' Jun ',&
      ' Jul ',' Aug ',' Sep ',' Oct ',' Nov ',' Dec ']

    CALL DATE_AND_TIME(ZONE = zn, VALUES = values)
    znfrmt = zn(1:3)//":"//zn(4:5)//"(hh:mm) UTC"
    znfrmt = trim(znfrmt)
    WRITE(  dd,'(i2)') values(3)
    WRITE(yyyy,'(i4)') values(1)
    write(hr, '(i2)')  values(5)
    write(mnt, '(i2)')  values(6)
    write(scnd, '(i2)')  values(7)
    if(mnt(1:1) == " ")   mnt(1:1) = "0"
    if(scnd(1:1) == " ") scnd(1:1) = "0"

    dt_time = dd//mons(values(2))//yyyy//spspt// &
             hr//tspt//mnt//tspt//scnd//spspt//znfrmt
    dt_time = trim(dt_time)

  END subroutine get_DDMonYYTimeZone

end module main_module
