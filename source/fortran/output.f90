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

end module process_output_fortran
