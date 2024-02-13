module final_module
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  contains

  subroutine final_header(ifail)
    use process_output, only: oheadr
    use constants, only: nout
    implicit none

    integer, intent(in) :: ifail

    if (ifail == 1) then
      call oheadr(nout,'Final Feasible Point')
    else
      call oheadr(nout,'Final UNFEASIBLE Point')
    end if
  end subroutine final_header

  subroutine no_optimisation()
    use constants, only: mfile, nout
    use numerics, only: neqns, nineqns, ipeqns, icc, lablcc, rcm, norm_objf
    use process_output, only: int_to_string3, ovarre, ocmmnt, oblnkl, osubhd, &
      oheadr
    use constraints, only: constraint_eqns
    use function_evaluator, only: funfom

    implicit none

    integer :: inn
    real(dp), dimension(ipeqns) :: con1, con2, err
    character(len=1), dimension(ipeqns) :: sym
    character(len=10), dimension(ipeqns) :: lab


    call oheadr(nout,'Numerics')
    call ocmmnt(nout,'PROCESS has performed a run witout optimisation.')
    call oblnkl(nout)

    ! Evaluate objective function
    call funfom(norm_objf)
    call ovarre(mfile,'Normalised objective function','(norm_objf)',norm_objf)

    ! Print the residuals of the constraint equations
    call constraint_eqns(neqns+nineqns,-1,con1,con2,err,sym,lab)
    write(nout,120)
    120 format(t48,'physical',t73,'constraint',t100,'normalised')
    write(nout,130)
    130 format(t47,'constraint',t74,'residue',t101,'residue')
    call oblnkl(nout)
    do inn = 1,neqns
      write(nout,140) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
      lab(inn),err(inn),lab(inn),con1(inn)
      call ovarre(mfile,lablcc(icc(inn))//' normalised residue', &
      '(eq_con'//int_to_string3(icc(inn))//')',con1(inn))
    end do

    140 format(t2,i4,t8,a33,t46,a1,t47,1pe12.4,t60,a10,t71,1pe12.4,t84,a10,t98,1pe12.4)

    if (nineqns > 0) then
      call osubhd(nout, &
      'The following inequality constraint residues should be greater than or approximately equal to zero :')

      do inn = neqns+1,neqns+nineqns
        write(nout,140) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
        lab(inn), err(inn), lab(inn)
        call ovarre(mfile,lablcc(icc(inn)),'(ineq_con'//int_to_string3(icc(inn))//')',rcm(inn))
      end do
    end if
  end subroutine no_optimisation

  subroutine final_output()
    use constants, only: iotty
    use numerics, only: nfev1, ncalls, xcm, ioptimz, nviter, &
      nvar
    use define_iteration_variables, only: loadxc
    implicit none

    if (nfev1 == 0) then  !  no HYBRD call
      !if (nviter == 1) then
      !    write(iotty,10) nviter,ncalls
      !else
      !    write(iotty,20) nviter,ncalls
      !end if
    else if (nviter == 0) then  !  no VMCON call
      if (nfev1 == 1) then
        write(iotty,30) nfev1,ncalls
      else
        write(iotty,40) nfev1,ncalls
      end if
    else if (nfev1 == 1) then ! (unlikely that nviter is also 1...)
      write(iotty,50) nfev1,nviter,ncalls
    else if (nviter == 1) then ! (unlikely that nfev1 is also 1...)
        write(iotty,60) nfev1,nviter,ncalls
      else
        write(iotty,70) nfev1,nviter,ncalls
      end if

    30 format( &
          t2,'The HYBRD point required ',i5,' iteration',/, &
          t2,'There were ',i6,' function calls')
    40 format( &
          t2,'The HYBRD point required ',i5,' iterations',/, &
          t2,'There were ',i6,' function calls')
    50 format( &
          t2,'The HYBRD point required ',i5,' iteration',/, &
          t2,'The optimisation required ',i5,' iterations',/, &
          t2,'There were ',i6,' function calls')
    60 format( &
          t2,'The HYBRD point required ',i5,' iterations',/, &
          t2,'The optimisation required ',i5,' iteration',/, &
          t2,'There were ',i6,' function calls')
    70 format( &
          t2,'The HYBRD point required ',i5,' iterations',/, &
          t2,'The optimisation required ',i5,' iterations',/, &
          t2,'There were ',i6,' function calls')
  end subroutine final_output
end module final_module
