module stellarator_module

  !! Module containing stellarator routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of the first wall, blanket and shield components
  !! of a fusion power plant.

  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
   use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  ! scaling parameters to reference point.
  real(dp) :: f_n, f_r, f_aspect, f_b, f_i, f_a


  logical :: first_call = .true.
  logical :: first_call_stfwbs = .true.

contains

  subroutine init_stellarator_module
    !! Initialise module variables
    implicit none

    first_call = .true.
    first_call_stfwbs = .true.
    f_n = 0.0D0
    f_r = 0.0D0
    f_a = 0.0D0
    f_b = 0.0D0
    f_i = 0.0D0
  end subroutine init_stellarator_module

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine stinit

    !! Routine to initialise the variables relevant to stellarators
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: F Warmer, IPP Greifswald
    !! None
    !! This routine initialises the variables relevant to stellarators.
    !! Many of these may override the values set in routine
    !! <A HREF="initial.html">initial</A>.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use build_variables, only: gapoh, iohcl, ohcth, tfootfi
    use current_drive_variables, only: irfcd
    use pfcoil_variables, only: ohhghf
    use physics_variables, only: aspect, dnbeta, kappa, kappa95, q, rmajor, &
      triang, hfac, tauscl
    use numerics, only: boundl, boundu
    use stellarator_variables, only: istell
    use tfcoil_variables, only: n_tf
    use times_variables, only: tburn, tcycle, tdown, tdwell, t_fusion_ramp, tohs, &
      tpulse, tqnch, tramp
		use global_variables, only: icase
		use constants, only: pi, rmu0, nout
    implicit none

    !  Arguments

    !  Local variables

    !real(dp) :: fsum

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! This routine is called before (!!!) the input file. put everything that depends on the input file in stcaller
    if (istell == 0) return

    boundu(1) = 40.0D0 ! allow higher aspect ratio
    !  Numerics quantities

    !boundl(1) = 5.0D0


    !boundu(3) = 30.0D0
    !boundu(29) = 20.0D0

    !  These lines switch off tokamak specifics (solenoid, pf coils, pulses etc.).
    !  Are they still up to date? (26/07/22 JL)

    !  Build quantities

    ohcth = 0.0D0
    iohcl = 0
    ohhghf = 0.0D0
    gapoh = 0.0D0
    tfootfi = 1.0D0

    !  Physics quantities

    dnbeta = 0.0D0
    kappa95 = 1.0D0
    triang = 0.0D0
    q = 1.03D0

    !  Turn off current drive

    irfcd = 0

    !  Times for different phases

    tramp = 0.0D0
    tohs = 0.0D0
    tburn = 3.15576D7  !  one year
    tqnch = 0.0D0
    tpulse = tohs + t_fusion_ramp + tburn + tqnch
    tdown  = tramp + tohs + tqnch + tdwell
    tcycle = tramp + tohs + t_fusion_ramp + tburn + tqnch + tdwell

  end subroutine stinit
end module stellarator_module
