module maths_library

  !! Library of mathematical and numerical routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains a large number of routines to enable
  !! PROCESS to perform a variety of numerical procedures, including
  !! linear algebra, zero finding, integration and minimisation.
  !! The module is an amalgamation of the contents of several
  !! different pre-existing PROCESS source files, which themselves
  !! were derived from a number of different numerical libraries
  !! including BLAS and MINPAC.
  !! http://en.wikipedia.org/wiki/Gamma_function
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  !use process_output

  implicit none

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ------------------------------------------------------------------------
  pure function variable_error(variable)
      real(dp), intent(in) ::variable
      logical::variable_error

      if((variable/=variable).or.(variable<-9.99D99).or.(variable>9.99D99))then
          variable_error = .TRUE.
      else
          variable_error = .FALSE.
      end if

  end function variable_error

end module maths_library
