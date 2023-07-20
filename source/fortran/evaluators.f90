! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module function_evaluator

  !! Module containing function evaluators for HYBRD and VMCON
  !! solvers
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains the function evaluators required
  !! by the two equation solvers in the code.

  !! None
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fcnhyb() is commented out temporarily. It calls the caller() subroutine, which
! is being moved to Python. fcnhyb() is passed into maths_library.hybrd() as an
! external subroutine argument, which calls fcnhyb() within a goto block. This
! is difficult to unravel when converting to Python (so that caller() is called
! from Python), so it has been decided to comment out hybrd() and fcnhyb()
! temporarily, disabling the non-optimising solver, to allow Python conversion
! work using the optimising solver (vmcon()) to continue.

!   subroutine fcnhyb(n,xc,rc,iflag)

!     !! Function evaluator for EQSOLV
!     !! author: P J Knight, CCFE, Culham Science Centre
!     !! n : input integer : Number of equations and unknowns
!     !! xc(n) : input/output real array : On input XC must contain
!     !! an initial estimate of the solution vector. On output XC
!     !! contains the final estimate of the solution vector.
!     !! rc(n) : output real array : Functions evaluated at the output XC
!     !! iflag : input/output integer : Terminate execution of EQSOLV
!     !! by setting IFLAG to a negative integer.
!     !! This subroutine is the function evaluator for
!     !! <A HREF="eqsolv.html">EQSOLV</A> (q.v.).
!     !! None
!     !
!     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     use constraints, only: constraint_eqns
!     use numerics, only: neqns
!     use caller_module, only: caller
!     implicit none

!     !  Arguments

!     integer, intent(in) :: n
!     real(dp), dimension(n), intent(inout) :: xc
!     real(dp), dimension(n), intent(out) :: rc
!     integer, intent(inout) :: iflag

!     !  Local variables

!     integer :: ncon, nvars

!     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     nvars = neqns
!     ncon = neqns

!     call caller(xc,nvars)
!     call constraint_eqns(ncon,rc,-1)

!     !  Set iflag < 0 if program is to be terminated here.

!     iflag = 1 * iflag

!   end subroutine fcnhyb

  subroutine funfom(fc)

    !! Objective function evaluator for VMCON
    !! author: P J Knight, CCFE, Culham Science Centre
    !! fc : output real : value of objective function at the output point
    !! This routine evaluates the value of the objective function
    !! i.e. the (normalised) figure-of-merit, at the nvar-dimensional
    !! point of interest.
    !! <P>Each equation for <CODE>fc<CODE> gives a value of the
    !! order of unity for the sake of the numerics.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use global_variables, only: xlabel, iscan_global
		use constants, only: nout, iotty, mfile
		use constraints, only: constraint_eqns
		use cost_variables, only: concost, cfactr, cdirt, ireactor, iavail, coe
		use current_drive_variables, only: bigq, porbitlossmw, pinjmw
		use divertor_variables, only: hldiv
		use error_handling, only: idiags, fdiags, errors_on, report_error
		use heat_transport_variables, only: pnetelmw
    use numerics, only: minmax
		use physics_variables, only: powfmw, bt, rmajor, wallmw, aspect, pohmmw
		use pf_power_variables, only: srcktpm
		use process_output, only: int_to_string3
		use tfcoil_variables, only: tfcmw
		use times_variables, only: tburn
    implicit none

    !  Arguments

    real(dp), intent(out) :: fc
!    real(c_double), intent(out) :: fc
    !  Local variables

    integer :: iab
    real(dp) :: sgn

!        write(*,*) 'Figure of merit 2 (fusion power / input power) is not used.'
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iab = abs(minmax)
    sgn = sign(1.0D0, real(minmax, kind(1.0D0)))

    !  If sgn is -1 the value of fc will be maximised
    !  If sgn is +1 the value of fc will be minimised

    select case (iab)

    case (1)  !  major radius
       fc = sgn * 0.2D0 * rmajor

    case (2)  !  fusion power / input power
        write(*,*) 'Figure of merit 2 (fusion power / input power) is not used.'
        write(*,*) 'Figure of merit 5 (fusion gain Q) is available.'
        stop 1
       ! fc = sgn * powfmw / (pinjmw + porbitlossmw + tfcpmw + ppump/1.0D6)

    case (3)  !  neutron wall load
       fc = sgn * wallmw

    case (4)  !  TF coil + PF coil power
       fc = sgn * (tfcmw + 1.0D-3*srcktpm)/10.0D0

   case (5)  !  Q = fusion gain  Issue #540
       fc = sgn * powfmw / (pinjmw + porbitlossmw + pohmmw)
       !fc = sgn * powfmw / pinjmw

    case (6)  !  cost of electricity
       fc = sgn * coe/100.0D0

    case (7)  !  direct/constructed/capital cost
       if (ireactor == 0) then
          fc = sgn * cdirt/1.0D3
       else
          fc = sgn * concost/1.0D4
       end if

    case (8)  !  aspect ratio
       fc = sgn * aspect

    case (9)  !  divertor heat load
       fc = sgn * hldiv

    case (10)  !  toroidal field on axis
       fc = sgn * bt

    case (11)  !  injection power
       fc = sgn * pinjmw

    case (12)  !  hydrogen production capital cost
       ! #506 OBSOLETE
       write(*,*) 'Figure of Merit 13 (Hydrogen production) is no longer supported.'
       stop 1
    case (13)  !  hydrogen production rate
       ! #506 OBSOLETE
       write(*,*) 'Figure of Merit 13 (Hydrogen production) is no longer supported.'
       stop 1

    case (14)  !  pulse length
       fc = sgn * tburn / 2.0D4

    case (15)  !  plant availability factor (N.B. requires iavail = 1)

       if (iavail /= 1) call report_error(23)

       fc = sgn * cfactr

    case (16)  !  major radius/burn time
       fc = sgn * ( 0.95d0 * (rmajor/9.0d0) - 0.05d0 * (tburn/7200.d0) )

    case (17)  !  net electrical output
       fc = sgn * pnetelmw / 500.0d0

   case (18)  !  Null figure of merit
      fc = 1d0

   case (19)  !  major radius/burn time
      fc = sgn * ( -0.5d0 * (bigq/20.0D0) - 0.5d0 * (tburn/7200.d0) )

    case default
       idiags(1) = iab ; call report_error(24)

    end select

    !  Crude method of catching NaN errors

    if ((abs(fc) > 9.99D99).or.(fc /= fc)) then
       idiags(1) = iab ; call report_error(25)
    end if

  end subroutine funfom

end module function_evaluator
