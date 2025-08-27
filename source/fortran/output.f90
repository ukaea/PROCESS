! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module process_output_fortran

  !! Module containing routines to produce a uniform output style
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains a number of routines that allow the
  !! program to write output to a file unit in a uniform style.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public



contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write(file, string)
    !! Write a string to file.
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    implicit none

    !  Arguments
    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    write(file,*) trim(string)
  end subroutine write

  subroutine dblcol(file, desc, val1, val2)
    !! Write a description and 2 columns of values to 2dp in standard notation.
    !! file : input integer : Fortran output unit identifier
    !! desc : input character string : Character string to be used
    !! val1 : input real : Value of the left variable
    !! val2 : input real : Value of the right variable
    implicit none

    !  Arguments
    integer, intent(in) :: file
    character(len=70), intent(in) :: desc
    real(8), intent(in) :: val1, val2

    write(file,10) desc, val1, val2
    10 format(1x,a,t75,f10.2,t100,f10.2)
  end subroutine dblcol


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ovarst(file,descr,varnam,value)

    !! Routine to print out the details of a character variable
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the variable
    !! varnam : input character string : Name of the variable
    !! value : input character string : Value of the variable
    !! This routine writes out the description, name and value of a
    !! character string variable.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: sqsumsq
		use constants, only: mfile
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr, varnam, value

    !  Local variables

    character(len=72) :: dum72
    character(len=30) :: dum20

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr and varnam with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum72 = descr
    dum20 = varnam

    if (file /= mfile) then
       write(file,10) dum72, dum20, value
    end if

    call underscore(dum72)
    call underscore(dum20)
    write(mfile,10) dum72, dum20, value

10  format(1x,a,t75,a,t110,a)

  end subroutine ovarst

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine obuild(file,descr,thick,total,variable_name)

    !! Routine to print out a description, the thickness and
    !! summed build of a component of the radial or vertical build
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the component
    !! thick : input real : Thickness of the component (m)
    !! total : input real : Total build, including this component (m)
    !! This routine writes out a description, the thickness and
    !! summed build of a component of the radial or vertical build.
    !!     !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: boundl, boundu
		use constants, only: electron_charge
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr
    character(len=*), optional :: variable_name
    real(8), intent(in) :: thick, total

    !  Local variables

    character(len=50) :: dum30
    character(len=40) :: dum20

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr with dummy string of the correct length.
    !  This counters problems that would occur if the original string
    !  was the wrong length.

    dum30 = descr
    if (present(variable_name)) then
        dum20 = variable_name
    else
        dum20 = ''
    end if

    write(file,10) dum30, thick, total, dum20
10  format(1x,a,t42,f10.3,t58,f10.3,t71,a)

  end subroutine obuild

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine underscore(string)

    !! Routine that converts spaces in a string to underscores
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : input/output string : character string of interest
    !! This routine converts any space characters in the string
    !! to underscore characters.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: active_constraints, boundu, boundl
		use constants, only: electron_charge
    implicit none

    !  Arguments

    character(len=*), intent(inout) :: string

    !  Local variables

    integer :: loop, i

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    i = index(string, ' ')
    if (i > 0) then
       do loop = i, len(string)
          if (string(loop:loop) == ' ') string(loop:loop) = '_'
       end do
    end if

  end subroutine underscore

end module process_output_fortran
