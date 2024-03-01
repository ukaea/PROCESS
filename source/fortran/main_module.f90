#ifndef COMMSG
#error COMMSG not defined!
#endif

#ifndef tagno
#error tagno not defined!
#endif

#ifndef branch_name
#error branch_name not defined!
#endif

#ifndef untracked
#error untracked not defined!
#endif

module main_module

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

contains

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inform(progid)

  !! Routine to obtain information about the program being executed
  !! author: P J Knight, CCFE, Culham Science Centre
  !! progid(0:10) : output string array : Strings containing useful info
  !! This subroutine uses system calls to identify the user, date,
  !! machine etc. for the present run, and stores the information
  !! in a character string array.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: nout
  implicit none

  !  Arguments
  character(len=110), dimension(0:10) :: progid

  !  Local variables
  character(len=10) :: progname
  character(len=98) :: executable
  character(len=*), parameter :: progver = &  !  Beware: keep exactly same format...
       '3.0.2   Release Date :: 2024-01-25'
  character(len = 50) :: dt_time
  character(len=72), dimension(10) :: id

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Program name
  progname = 'PROCESS'
  call get_command_argument(0, executable)
  call get_DDMonYYTimeZone(dt_time)
  id(1) = trim(dt_time) ! values(3)//"/"// values(2)//"/"// values(1)  ! 5 6 7!date
  call getlog(id(2))    ! Get user ID
  call hostnm(id(3))    ! Get host name
  call getcwd(id(4))    ! Get current working directory

  !  Annotate information and store in PROGID character array
  !  for use in other program units via the routine argument

  progid(1) = '  Program : ' // executable
  progid(2) = '  Version : ' // progver
  progid(3) = 'Date/time : ' // id(1)
  progid(4) = '     User : ' // id(2)
  progid(5) = ' Computer : ' // id(3)
  progid(6) = 'Directory : ' // id(4)

  !  Summarise most useful data, and store in progid(0)
  progid(0) = trim(progname) // ' ' // trim(progver(:7)) // &
       ' : Run on ' // trim(id(1)) // ' by ' // trim(id(3))

end subroutine inform

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_summary

  !! Routine to print out a summary header
  !! author: P J Knight, CCFE, Culham Science Centre
  !! None
  !! This routine prints out a header summarising the program
  !! execution details, plus a list of the active iteration
  !! variables and constraint equations for the run.
  !! A User's Guide to the PROCESS Systems Code, P. J. Knight,
  !! AEA Fusion Report AEA FUS 251, 1993
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: nout, mfile, iotty, mfile
  use maths_library, only: integer2string, integer3string
  use global_variables, only: maxcal, fileprefix, icase, runtitle
  use numerics, only: nvar, neqns, ioptimz, nineqns, epsvmc, minmax, icc, &
    lablcc, lablmm
  use process_output, only: ocentr, oblnkl, ocmmnt, ostars, ovarst, ovarin
  use physics_variables, only: te
  implicit none

  !  Local variables
  integer, parameter :: width = 110
  integer :: lap, ii, outfile
  character(len = 110) :: progid(0:10)
  character(len = 9)   :: vstring
  character(len = 8)   :: date
  character(len = 10)  :: time
  character(len = 12)  :: dstring
  character(len = 7)   :: tstring
  character(len = 10)  :: ustring
  character(len = 100) :: rstring
  character(len = 60)  :: fom_string
  character(len = 14)  :: minmax_string
  character(len = 10)  :: eps_string
  character :: minmax_sign

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Obtain execution details for this run
  call inform(progid)

  !  Print code banner + run details to screen and output file
  do lap = 1,2
     if (lap == 1) then
        outfile = iotty
     else
        outfile = nout
     end if

     ! PROCESS code header
     call oblnkl(outfile)
     call ostars(outfile, width)
     call ocentr(outfile,'PROCESS', width)
     call ocentr(outfile,'Power Reactor Optimisation Code', width)
     call ostars(outfile, width)
     call oblnkl(outfile)

     !  Run execution details
     call ocmmnt(outfile, progid(1))  !  program name
     call ocmmnt(outfile, progid(2))  !  version
     if (untracked > 0) then  ! tag number
       call ocmmnt(outfile, '  Tag No. : '//tagno//' code contains untracked changes')
     else
       call ocmmnt(outfile, '  Tag No. : '//tagno)
     end if
     call ocmmnt(outfile, '   Branch : '//branch_name)
     call ocmmnt(outfile, '  Git log : '// &
     COMMSG)  !  Last git com message
     call ocmmnt(outfile, progid(3))  !  date/time
     call ocmmnt(outfile, progid(4))  !  user
     call ocmmnt(outfile, progid(5))  !  computer
     call ocmmnt(outfile, progid(6))  !  directory
     if (trim(fileprefix) == "") then
       call ocmmnt(outfile, '    Input : IN.DAT')  !  input file name
     else
       call ocmmnt(outfile, '    Input : '//trim(fileprefix))  !  input file name
     end if
     call ocmmnt(outfile, 'Run title : '//trim(runtitle))   ! run title
     call ocmmnt(outfile, ' Run type : Reactor concept design: '// trim(icase) // ', (c) CCFE')

     call oblnkl(outfile)
     call ostars(outfile, width)
     call oblnkl(outfile)
     call ocmmnt(outfile, '  Equality constraints : '//integer2string(neqns))  !  Number of equality constraints
     call ocmmnt(outfile, 'Inequality constraints : '//integer2string(nineqns))  !  Number of inequality constraints
     call ocmmnt(outfile, '     Total constraints : '//integer2string(neqns+nineqns))  !  Number of constraints
     call ocmmnt(outfile, '   Iteration variables : '//integer2string(nvar))  !  Number of iteration variables
     call ocmmnt(outfile, '        Max iterations : '//integer3string(maxcal))  !  Max number of iterations

     if (minmax > 0) then
      minmax_string = '  -- minimise '
      minmax_sign = "+"
     else
      minmax_string = '  -- maximise '
      minmax_sign = "-"
     end if
     fom_string = lablmm(abs(minmax))
     call ocmmnt(outfile, '      Figure of merit  : '//minmax_sign//integer2string(abs(minmax))//minmax_string//fom_string) ! Figure of merit

     write(eps_string, '(ES8.2)') epsvmc
     call ocmmnt(outfile, ' Convergence parameter : '//eps_string)  !  Convergence parameter
     call oblnkl(outfile)
     call ostars(outfile, width)
  end do

  call oblnkl(outfile)
  call ocmmnt(nout,'(Please include this header in any models, presentations and papers based on these results)')
  call oblnkl(nout)
  call ostars(nout, width)
  ! Issue #270
  call oblnkl(outfile)
  call ocmmnt(nout,'Quantities listed in standard row format are labelled as follows in columns 112-114:')
  call ocmmnt(nout,'ITV : Active iteration variable (in any output blocks)')
  call ocmmnt(nout,'OP  : Calculated output quantity')
  call ocmmnt(nout,'Unlabelled quantities in standard row format are generally inputs')
  call ocmmnt(nout,'Note that calculated quantities may be trivially rescaled from inputs, or equal to bounds which are input.')
  ! MDK Note that the label must be exactly three characters or none - I don't know how to fix this.

  !  Beware of possible future changes to the progid(...) layouts

  !  Relies on an internal read statement
  vstring = progid(2)(13:21)
  call ovarst(mfile,'PROCESS version number','(procver)','"'//vstring//'"')

  call date_and_time(date=date, time=time)

  !  Date output in the form "DD/MM/YYYY" (including quotes)
  dstring = '"'//date(7:8)//'/'//date(5:6)//'/'//date(1:4)//'"'
  call ovarst(mfile,'Date of run','(date)',dstring)

  !  Time output in the form "hh:mm" (including quotes)
  tstring = '"'//time(1:2)//':'//time(3:4)//'"'
  call ovarst(mfile,'Time of run','(time)',tstring)

  ustring = '"'//trim(progid(4)(13:20))//'"'
  call ovarst(mfile,'User','(username)',ustring)

  rstring = '"'//runtitle//'"'
  call ovarst(mfile,'PROCESS run title','(runtitle)',rstring)

  rstring = '"'//tagno//'"'
  call ovarst(mfile,'PROCESS tag number','(tagno)',rstring)

  rstring = '"'//branch_name//'"'
  call ovarst(mfile,'PROCESS git branch name','(branch_name)',rstring)

  rstring = '"'//COMMSG//'"'
  call ovarst(mfile,'PROCESS last commit message','(commsg)',rstring)

  call ovarst(mfile,'Input filename','(fileprefix)','"'//trim(fileprefix//'"'))

  if (ioptimz == -2) then
     call ovarin(mfile,'Optimisation switch','(ioptimz)',ioptimz)
     call ovarin(mfile,'Figure of merit switch','(minmax)',minmax)
  end if

#ifndef unit_test

  ! MDK Only print out the constraints here for HYBRD.
  ! For VMCON they are printed out later with residues.
  call oblnkl(nout)
  if (ioptimz == -1) then
      call ocmmnt(nout, 'The following constraint equations have been imposed,')
      call ocmmnt(nout, 'but limits will not be enforced by the code :')
      write(nout,30)
30    format(t10,'icc',t25,'label')
      call oblnkl(nout)
      write(nout,40) (ii,icc(ii),lablcc(icc(ii)), ii=1,neqns+nineqns)
40    format(t1,i3,t10,i3,t18,a33)
  end if

#endif

end subroutine run_summary

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Minpack Copyright Notice (1999) University of Chicago.  All rights reserved

! Redistribution and use in source and binary forms, with or
! without modification, are permitted provided that the
! following conditions are met:

! 1. Redistributions of source code must retain the above
! copyright notice, this list of conditions and the following
! disclaimer.

! 2. Redistributions in binary form must reproduce the above
! copyright notice, this list of conditions and the following
! disclaimer in the documentation and/or other materials
! provided with the distribution.

! 3. The end-user documentation included with the
! redistribution, if any, must include the following
! acknowledgment:

!    "This product includes software developed by the
!    University of Chicago, as Operator of Argonne National
!    Laboratory.

! Alternately, this acknowledgment may appear in the software
! itself, if and wherever such third-party acknowledgments
! normally appear.

! 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
! WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
! UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
! THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
! OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
! OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
! USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
! THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
! DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
! UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
! BE CORRECTED.

! 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
! HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
! ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
! INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
! ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
! PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
! SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
! (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
! EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
! POSSIBILITY OF SUCH LOSS OR DAMAGES.

! eqslv() has been temporarily commented out. Please see the comment in
! function_evaluator.fcnhyb() for an explanation.

! subroutine eqslv(ifail)

!   !! Routine to call the non-optimising equation solver
!   !! author: P J Knight, CCFE, Culham Science Centre
!   !! ifail   : output integer : error flag
!   !! This routine calls the non-optimising equation solver.
!   !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
!   !
!   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   use constants, only: nout, mfile, iotty
!   use constraints, only: constraint_eqns
!   use function_evaluator, only: fcnhyb
!   use error_handling, only: idiags, fdiags, errors_on, report_error
!   use numerics, only: ipeqns, epsfcn, factor, ftol, iptnt, ncalls, &
!     neqns, nfev1, nfev2, sqsumsq, xcm, rcm, xcs, resdl, scafc, ixc, lablxc, &
!     icc, lablcc, eqsolv
!   use process_output, only: ovarin, oblnkl, ocmmnt, oheadr, osubhd, &
!     ovarre, int_to_string3
!   use physics_variables, only: bt, aspect, rmajor, powfmw, wallmw
!   use define_iteration_variables, only: loadxc
!   implicit none

!   !  Arguments
!   integer, intent(out) :: ifail

!   !  Local variables
!   integer :: inn,nprint,nx
!   real(dp) :: sumsq
! !  real(dp), dimension(iptnt) :: wa
!   real(dp) :: wa(iptnt)
!   real(dp), dimension(ipeqns) :: con1, con2, err
!   character(len = 1), dimension(ipeqns) :: sym
!   character(len = 10), dimension(ipeqns) :: lab

!   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   ncalls = 0
!   nfev1 = 0
!   nfev2 = 0
!   nprint = 0

!   !  Use HYBRD to find a starting point
!   call loadxc
!   call eqsolv(fcnhyb,neqns,xcm,rcm,ftol,epsfcn,factor,nprint,ifail, &
!        wa,iptnt,resdl,nfev1)

!   !  Turn on error reporting
!   errors_on = .true.

!   !  Print out information on solution
!   call oheadr(nout,'Numerics')
!   call ocmmnt(nout, &
!        'PROCESS has performed a HYBRD (non-optimisation) run,')

!   if (ifail /= 1) then
!      call ocmmnt(nout,'but could NOT find a feasible set of parameters.')
!      call oblnkl(nout)
!      call ovarin(nout,'Number of iteration variables and constraints','(neqns)',neqns)
!      call ovarin(nout,'HYBRD error flag','(ifail)',ifail)

!      call oheadr(iotty,'PROCESS COULD NOT FIND A FEASIBLE SOLUTION')
!      call ovarin(iotty,'HYBRD error flag (ifail)','',ifail)
!      call oblnkl(iotty)

!      idiags(1) = ifail ; call report_error(131)

!   else
!      call ocmmnt(nout,'and found a feasible set of parameters.')
!      call oblnkl(nout)
!      call ovarin(nout,'HYBRD error flag','(ifail)',ifail)
!      call oblnkl(nout)
!      call oheadr(iotty,'PROCESS found a feasible solution')
!   end if

!   !  Sum the square of the residuals
!   sumsq = 0.0D0
!   do nx = 1,neqns
!      sumsq = sumsq + rcm(nx)**2
!   end do
!   sqsumsq = sqrt(sumsq)

!   call ovarre(nout,'Square root of the sum of squares of the constraint residuals','(sqsumsq)',sqsumsq, 'OP ')

!   !  If necessary, write out a relevant error message
!   if (ifail /= 1) then
!      call oblnkl(nout)
!      call herror(ifail)
!      call oblnkl(iotty)
!   else
!      !  Show a warning if the constraints appear high even if allegedly converged
!      if (sqsumsq >= 1.0D-2) then
!         call oblnkl(nout)
!         call ocmmnt(nout,'WARNING: Constraint residues are HIGH; consider re-running')
!         call ocmmnt(nout,'   with lower values of FTOL to confirm convergence...')
!         call ocmmnt(nout,'   (should be able to get down to about 1.0E-8 okay)')

!         call ocmmnt(iotty,'WARNING: Constraint residues are HIGH; consider re-running')
!         call ocmmnt(iotty,'   with lower values of FTOL to confirm convergence...')
!         call ocmmnt(iotty,'   (should be able to get down to about 1.0E-8 okay)')
!         call oblnkl(iotty)

!         fdiags(1) = sqsumsq ; call report_error(133)

!      end if
!   end if

!   call osubhd(nout,'The solution vector is comprised as follows :')

!   write(nout,10)
! 10 format(t5,'i',t23,'final',t33,'fractional',t46,'residue')

!   write(nout,20)
! 20 format(t23,'value',t35,'change')

!   call oblnkl(nout)

!   do inn = 1,neqns
!      xcs(inn) = xcm(inn)*scafc(inn)
!      write(nout,30) inn,lablxc(ixc(inn)),xcs(inn),xcm(inn),resdl(inn)
!      call ovarre(mfile,lablxc(ixc(inn)),'(itvar'//int_to_string3(inn)//')',xcs(inn))
!   end do
! !30 format(t2,i4,t8,a9,t19,1pe12.4,1pe12.4,1pe12.4)
! ! Make lablxc longer
! 30 format(t2,i4,t8,a30,t39,1pe12.4,1pe12.4,1pe12.4)

!   call osubhd(nout, &
!        'The following constraint residues should be close to zero :')

!   call constraint_eqns(neqns,con1,-1,con2,err,sym,lab)
!   write(nout,40)
! 40 format(t48,'physical',t73,'constraint',t100,'normalised')
!   write(nout,50)
! 50 format(t47,'constraint',t74,'residue',t101,'residue')
!   call oblnkl(nout)
!   do inn = 1,neqns
!      write(nout,60) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
!           lab(inn),err(inn),lab(inn),con1(inn)
!      call ovarre(mfile,lablcc(icc(inn))//' normalised residue', &
!           '(normres'//int_to_string3(inn)//')',con1(inn))
!   end do
! 60 format(t2,i4,t8,a33,t46,a1,t47,1pe12.4,t60,a10,t71,1pe12.4,t84,a10,t98,1pe12.4)

! end subroutine eqslv

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
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
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
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
  use maths_library, only: nearly_equal, binomial, test_secant_solve
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
  call test_secant_solve()
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
