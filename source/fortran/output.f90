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

end module process_output_fortran
