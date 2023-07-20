module utilities

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

contains
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine upper_case(string,start,finish)

    !! Routine that converts a (sub-)string to uppercase
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : input string   : character string of interest
    !! start  : optional input integer  : starting character for conversion
    !! finish : optional input integer  : final character for conversion
    !! This routine converts the specified section of a string
    !! to uppercase. By default, the whole string will be converted.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(inout) :: string
    integer, optional, intent(in) :: start,finish

    !  Local variables

    character(len=1) :: letter
    character(len=27), parameter :: upptab = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_'
    integer :: loop, i

    integer :: first, last

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(start)) then
       first = start
    else
       first = 1
    end if

    if (present(finish)) then
       last = finish
    else
       last = len(string)
    end if

    if (first <= last) then
       do loop = first,last
          letter = string(loop:loop)
          i = index('abcdefghijklmnopqrstuvwxyz_',letter)
          if (i > 0) string(loop:loop) = upptab(i:i)
       end do
    end if

  end subroutine upper_case

end module utilities
