module global_variables
  !! author: J. Morris (UKAEA)
  !!
  !! This module contains miscellaneous global variables not well-suited to any
  !! of the other 'variables' modules.
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public

  character(len=48) :: icase
  !! power plant type

  character(len=180) :: runtitle
  !! short descriptive title for the run

  integer :: verbose
  !! switch for turning on/off diagnostic messages
  !!
  !! - =0 turn off diagnostics
  !! - =1 turn on diagnostics

  integer :: run_tests
  !! turns on built-in tests if set to 1

  integer :: maxcal
  !! maximum number of VMCON iterations

  character(len=400) :: fileprefix
  !! input file prefix

  character(len=400) :: output_prefix
  !! output file prefix

  character(len=25) :: xlabel
  !! scan parameter description label

  character(len=25) :: vlabel
  !! scan value name label

  character(len=25) :: xlabel_2
  !! scan parameter description label (2nd dimension)

  character(len=25) :: vlabel_2
  !! scan value name label (2nd dimension)

  integer :: iscan_global
  !! Makes iscan available globally.

  real(dp) :: convergence_parameter
  !! VMCON convergence parameter "sum"

  contains

  subroutine init_global_variables
    !! Initialise global variables
    implicit none

    icase = 'Steady-state tokamak model'
    runtitle = "Run Title (change this line using input variable 'runtitle')"
    verbose = 0
    run_tests = 0
    maxcal = 200
    fileprefix = ""
    output_prefix = ""
    xlabel = ""
    vlabel = ""
    xlabel_2 = ""
    vlabel_2 = ""
    iscan_global = 0
    convergence_parameter = 0.0D0
  end subroutine init_global_variables
end module global_variables
