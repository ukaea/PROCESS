module stellarator_module

  !! Module containing stellarator routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of the first wall, blanket and shield components
  !! of a fusion power plant.

  !!   !
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

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine stinit

    !! Routine to initialise the variables relevant to stellarators
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: F Warmer, IPP Greifswald
    !! None
    !! This routine initialises the variables relevant to stellarators.
    !! Many of these may override the values set in routine
    !! <A HREF="initial.html">initial</A>.
    !!     !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use build_variables, only: dr_cs_tf_gap, iohcl, dr_cs, tfootfi
    use current_drive_variables, only: i_hcd_calculations
    use pfcoil_variables, only: f_z_cs_tf_internal
    use physics_variables, only: aspect, beta_norm_max, kappa, kappa95, q95, rmajor, &
      triang, hfac, labels_confinement_scalings
    use numerics, only: boundl, boundu
    use stellarator_variables, only: istell
    use tfcoil_variables, only: n_tf_coils
    use times_variables, only: t_burn, t_cycle, tdown, t_between_pulse, t_fusion_ramp, t_current_ramp_up, &
      t_pulse_repetition, t_ramp_down, t_precharge
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

    dr_cs = 0.0D0
    iohcl = 0
    f_z_cs_tf_internal = 0.0D0
    dr_cs_tf_gap = 0.0D0
    tfootfi = 1.0D0

    !  Physics quantities

    beta_norm_max = 0.0D0
    kappa95 = 1.0D0
    triang = 0.0D0
    q95 = 1.03D0

    !  Turn off current drive

    i_hcd_calculations = 0

    !  Times for different phases

    t_precharge = 0.0D0
    t_current_ramp_up = 0.0D0
    t_burn = 3.15576D7  !  one year
    t_ramp_down = 0.0D0
    t_pulse_repetition = t_current_ramp_up + t_fusion_ramp + t_burn + t_ramp_down
    tdown  = t_precharge + t_current_ramp_up + t_ramp_down + t_between_pulse
    t_cycle = t_precharge + t_current_ramp_up + t_fusion_ramp + t_burn + t_ramp_down + t_between_pulse

  end subroutine stinit
end module stellarator_module
