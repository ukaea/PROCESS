module main_module

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

contains

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
end module main_module
