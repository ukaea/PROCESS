module optimiz_module
  implicit none

  contains

  subroutine write_out(n, m)
    !! Write out optimiser information to the OPT.DAT file
    ! #TODO Move to Python once file IO done in Python
    use constants, only: opt_file
    use numerics, only: icc, ixc
    implicit none

    integer, intent(in) :: n
    !! Number of variables
    integer, intent(in) :: m
    !! Number of constraints

    ! Array defined for optimizer data output only
    integer, dimension(n) :: ixc_opt_out
    integer, dimension(m) :: icc_opt_out
    integer :: i

    ! Write the VMCON setup in OPT.DAT
    do i = 1, m
      icc_opt_out(i) = icc(i)
    end do
    do i = 1, n
      ixc_opt_out(i) = ixc(i)
    end do

    write(opt_file, *) ' number of constrains'
    write(opt_file, '(I4)') m
    write(opt_file, *) ' '
    write(opt_file, *) ' Constrains selection'
    write(opt_file, '(I3,*(I4))') icc_opt_out
    write(opt_file, *) ' '
    write(opt_file, *) ' number of variables'
    write(opt_file, '(I4)') n
    write(opt_file, *) ' '
    write(opt_file, *) ' Variables selection'
    write(opt_file, '(I3,*(I4))') ixc_opt_out
    write(opt_file, *) ' '
    write(opt_file, *) ' '
    write(opt_file, *) ' n VMCOM iter | Figure of merit | VMCON conv      | constrains quad sum |   residual,   input values &
                    &and  FoM input gradients'
    write(opt_file, '(A,*(I18))') '  niter          abs(objf)         sum                sqsumsq ', icc_opt_out, ixc_opt_out&
                    &, ixc_opt_out
  end subroutine write_out
end module optimiz_module
