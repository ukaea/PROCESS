! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constraints
  !! author: J Morris
  !!
  !! Module defining the constraint equations and the routine that evaluates the
  !! constraint equations.

  ! Import modules
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  use error_handling, only: report_error, idiags, fdiags

  implicit none


contains

   subroutine constraint_eqn_015(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for L-H power threshold limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for L-H power threshold limit
      !! #=# physics
      !! #=#=# fl_h_threshold, p_l_h_threshold_mw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fl_h_threshold : input real : f-value for L-H power threshold
      !! p_l_h_threshold_mw : input real : L-H mode power threshold (MW)
      !! p_plasma_separatrix_mw : input real : power to conducted to the divertor region (MW)
      use constraint_variables, only: fl_h_threshold
      use physics_variables, only: p_l_h_threshold_mw, p_plasma_separatrix_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fl_h_threshold * p_plasma_separatrix_mw / p_l_h_threshold_mw
      tmp_con = p_l_h_threshold_mw
      tmp_err = p_l_h_threshold_mw - p_plasma_separatrix_mw / fl_h_threshold
      if (fl_h_threshold > 1.0D0) then
         tmp_symbol = '>'
      else
         tmp_symbol = '<'
      end if
      tmp_units = 'MW'

   end subroutine constraint_eqn_015

   subroutine constraint_eqn_016(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for net electric power lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for net electric power lower limit
      !! #=# heat_transport
      !! #=#=# fpnetel, pnetelin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpnetel : input real : f-value for net electric power
      !! pnetelmw : input real : net electric power (MW)
      !! pnetelin : input real : required net electric power (MW)
      use constraint_variables, only: fpnetel, pnetelin
      use heat_transport_variables, only: pnetelmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fpnetel * pnetelmw / pnetelin
      tmp_con = pnetelin
      tmp_err = pnetelmw - pnetelin / fpnetel
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_016

   subroutine constraint_eqn_021(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for minor radius lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for minor radius lower limit
      !! #=# physics
      !! #=#=# frminor, aplasmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! frminor : input real : f-value for minor radius limit
      !! rminor : input real : plasma minor radius (m)
      !! aplasmin : input real : minimum minor radius (m)
      use constraint_variables, only: frminor
      use physics_variables, only: rminor
      use build_variables, only: aplasmin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - frminor * rminor/aplasmin
      tmp_con = aplasmin * (1.0D0 - tmp_cc)
      tmp_err = aplasmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_021

   subroutine constraint_eqn_028(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fusion gain (big Q) lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fusion gain (big Q) lower limit
      !! #=# physics
      !! #=#=# fqval, bigqmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fqval : input real : pf-value for Q
      !! bigq : input real : Fusion gain; P_fusion / (P_injection + P_ohmic)
      !! bigqmin : input real : minimum fusion gain Q
      !! ignite : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! Obviously, ignite must be zero if current drive is required.
      !! If ignite=1, any auxiliary power is assumed to be used only
      !! during plasma start-up, and is excluded from all steady-state
      !! power balance calculations.
      use constraint_variables, only: fqval, bigqmin
      use current_drive_variables, only: bigq
      use physics_variables, only: ignite
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if plasma is not ignited ...
      if (ignite == 0) then
         tmp_cc =  1.0D0 - fqval * bigq/bigqmin
         tmp_con = bigqmin * (1.0D0 - tmp_cc)
         tmp_err = bigqmin * tmp_cc
         tmp_symbol = '>'
         tmp_units = ''
      ! if plasma is ignited report error
      else
         tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''
         call report_error(4)
      end if

   end subroutine constraint_eqn_028

   subroutine constraint_eqn_030(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for injection power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for injection power upper limit
      !! #=# current_drive
      !! #=#=# fpinj, p_hcd_injected_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! p_hcd_injected_total_mw : input real : total auxiliary injected power (MW)
      !! fpinj : input real : f-value for injection power
      !! p_hcd_injected_max : input real : Maximum allowable value for injected power (MW)
      use current_drive_variables, only: p_hcd_injected_total_mw, p_hcd_injected_max
      use constraint_variables, only: fpinj
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  p_hcd_injected_total_mw/p_hcd_injected_max - 1.0D0 * fpinj
      tmp_con = p_hcd_injected_max
      tmp_err = p_hcd_injected_max  - p_hcd_injected_total_mw / fpinj
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_030

   subroutine constraint_eqn_031(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil case stress upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil case stress upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fstrcase, sig_tf_case_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstrcase : input real : f-value for TF coil case stress
      !! sig_tf_case_max : input real : Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)
      !! sig_tf_case : input real : Constrained stress in TF coil case (Pa)
      use constraint_variables, only: fstrcase
      use tfcoil_variables, only: sig_tf_case_max, sig_tf_case
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  sig_tf_case/sig_tf_case_max - 1.0D0 * fstrcase
      tmp_con = sig_tf_case_max
      tmp_err = sig_tf_case_max - sig_tf_case / fstrcase
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_031

   subroutine constraint_eqn_032(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil conduit stress upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil conduit stress upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fstrcond, sig_tf_wp_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstrcond : input real : f-value for TF coil conduit stress
      !! sig_tf_wp_max : input real : Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)
      !! sig_tf_wp : input real : Constrained stress in TF conductor conduit (Pa)
      use constraint_variables, only: fstrcond
      use tfcoil_variables, only: sig_tf_wp_max, sig_tf_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  sig_tf_wp/sig_tf_wp_max - 1.0D0 * fstrcond
      tmp_con = sig_tf_wp_max
      tmp_err = sig_tf_wp_max - sig_tf_wp / fstrcond
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_032

   subroutine constraint_eqn_033(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil operating/critical J upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil operating/critical J upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fiooic, jwdgcrt
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fiooic : input real : f-value for TF coil operating current / critical
      !! jwdgcrt : input real : critical current density for winding pack (A/m2)
      !! j_tf_wp : input real : winding pack current density (A/m2)
      use constraint_variables, only: fiooic
      use tfcoil_variables, only: jwdgcrt, j_tf_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      if (fiooic > 0.7D0) call report_error(285)
      tmp_cc =  j_tf_wp/jwdgcrt - 1.0D0 * fiooic
      tmp_con = jwdgcrt * (1.0D0 - tmp_cc)
      tmp_err = j_tf_wp * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_033

   subroutine constraint_eqn_034(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil dump voltage upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil dump voltage upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fvdump, vdalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fvdump : input real : f-value for dump voltage
      !! vdalw : input real : max voltage across TF coil during quench (kV)
      !! vtfskv : input real : voltage across a TF coil during quench (kV)
      use constraint_variables, only: fvdump
      use tfcoil_variables, only: vdalw, vtfskv
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  vtfskv/vdalw - 1.0D0 * fvdump
      tmp_con = vdalw
      tmp_err = vdalw - vtfskv
      tmp_symbol = '<'
      tmp_units = 'V'

   end subroutine constraint_eqn_034

   subroutine constraint_eqn_035(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil J_wp/J_prot upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil J_wp/J_prot upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fjprot, jwdgpro
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjprot : input real : f-value for TF coil winding pack current density
      !! jwdgpro : input real : allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)
      !! j_tf_wp : input real : winding pack current density (A/m2)
      use constraint_variables, only: fjprot
      use tfcoil_variables, only: jwdgpro, j_tf_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  j_tf_wp/jwdgpro - 1.0D0 * fjprot
      tmp_con = jwdgpro
      tmp_err =  j_tf_wp - jwdgpro
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_035

   subroutine constraint_eqn_036(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil s/c temperature margin lower limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil s/c temperature margin lower limit (SCTF)
      !! #=# tfcoil
      !! #=#=# ftmargtf, tmargmin_tf
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftmargtf : input real : f-value for TF coil temperature margin
      !! tmargtf : input real : TF coil temperature margin (K)
      !! tmargmin_tf : input real : minimum allowable temperature margin : TF coils (K)
      use constraint_variables, only: ftmargtf
      use tfcoil_variables, only: tmargtf, tmargmin_tf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftmargtf * tmargtf/tmargmin_tf
      tmp_con = tmargmin_tf
      tmp_err = tmargmin_tf - tmargtf
      tmp_symbol = '>'
      tmp_units = 'K'

   end subroutine constraint_eqn_036

   subroutine constraint_eqn_037(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for current drive gamma upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for current drive gamma upper limit
      !! #=# current_drive
      !! #=#=# fgamcd, gammax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fgamcd : input real : f-value for current drive gamma
      !! gammax : input real : maximum current drive gamma
      !! eta_cd_norm_hcd_primary : input real : normalised current drive efficiency (1.0e20 A/W-m2)
      use constraint_variables, only: fgamcd, gammax
      use current_drive_variables, only: eta_cd_norm_hcd_primary
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  eta_cd_norm_hcd_primary/gammax - 1.0D0 * fgamcd
      tmp_con = gammax * (1.0D0 - tmp_cc)
      tmp_err = eta_cd_norm_hcd_primary * tmp_cc
      tmp_symbol = '<'
      tmp_units = '1E20 A/Wm2'

   end subroutine constraint_eqn_037

   subroutine constraint_eqn_038(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      implicit none
	   ! Dummy formal arguments, for compliance with interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_038

   subroutine constraint_eqn_039(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for first wall temperature upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for first wall temperature upper limit
      !! #=# fwbs
      !! #=#=# ftpeak, temp_fw_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftpeak : input real : f-value for first wall peak temperature
      !! temp_fw_max : input real : maximum temperature of first wall material (K) (i_thermal_electric_conversion>1)
      !! temp_fw_peak : input real : peak first wall temperature (K)
      use constraint_variables, only: ftpeak
      use fwbs_variables, only: temp_fw_max, temp_fw_peak
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! If the temperature peak == 0 then report an error
      if (temp_fw_peak < 1.0D0) call report_error(5)
      tmp_cc =  temp_fw_peak/temp_fw_max - 1.0D0 * ftpeak
      tmp_con = temp_fw_max * (1.0D0 - tmp_cc)
      tmp_err = temp_fw_peak * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'K'

   end subroutine constraint_eqn_039

   subroutine constraint_eqn_040(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for auxiliary power lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for auxiliary power lower limit
      !! #=# current_drive
      !! #=#=# fauxmn, auxmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fauxmn : input real : f-value for minimum auxiliary power
      !! p_hcd_injected_total_mw : input real : total auxiliary injected power (MW)
      !! auxmin : input real : minimum auxiliary power (MW)
      use constraint_variables, only: fauxmn, auxmin
      use current_drive_variables, only: p_hcd_injected_total_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fauxmn * p_hcd_injected_total_mw/auxmin
      tmp_con = auxmin * (1.0D0 - tmp_cc)
      tmp_err = auxmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_040

   subroutine constraint_eqn_041(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for plasma current ramp-up time lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for plasma current ramp-up time lower limit
      !! #=# times
      !! #=#=# ft_current_ramp_up, t_current_ramp_up_min
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ft_current_ramp_up : input real : f-value for plasma current ramp-up time
      !! t_current_ramp_up : input real : plasma current ramp-up time for current initiation (s)
      !! t_current_ramp_up_min : input real : minimum plasma current ramp-up time (s)
      use constraint_variables, only: ft_current_ramp_up, t_current_ramp_up_min
      use times_variables, only: t_current_ramp_up
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ft_current_ramp_up * t_current_ramp_up/t_current_ramp_up_min
      tmp_con = t_current_ramp_up_min * (1.0D0 - tmp_cc)
      tmp_err = t_current_ramp_up_min * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'sec'

   end subroutine constraint_eqn_041

   subroutine constraint_eqn_042(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for cycle time lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for cycle time lower limit
      !! #=# times
      !! #=#=# ftcycl, tcycmn
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftcycl : input real : f-value for cycle time
      !! t_cycle : input real : full cycle time (s)
      !! tcycmn : input real : minimum cycle time (s)
      use constraint_variables, only: ftcycl, tcycmn
      use times_variables, only: t_cycle
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the minimum cycle time == 0 report an error
      if (tcycmn < 1.0D0) call report_error(6)
      tmp_cc =  1.0D0 - ftcycl * t_cycle/tcycmn
      tmp_con = tcycmn * (1.0D0 - tmp_cc)
      tmp_err = tcycmn * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'sec'

   end subroutine constraint_eqn_042

   subroutine constraint_eqn_044(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for centrepost temperature upper limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for centrepost temperature upper limit (TART)
      !! #=# tfcoil
      !! #=#=# fptemp, ptempalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fptemp : input real : f-value for peak centrepost temperature
      !! ptempalw : input real : maximum peak centrepost temperature (K)
      !! tcpmax : input real :  peak centrepost temperature (K)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use constraint_variables, only: fptemp
      use tfcoil_variables, only: ptempalw, tcpmax
      use physics_variables, only: itart
      use tfcoil_variables, only:  i_tf_sup
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(8)

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         ptempalw = ptempalw - 273.15D0
         tcpmax = tcpmax - 273.15D0
      end if

      tmp_cc =   tcpmax/ptempalw - 1.0D0 * fptemp
      tmp_con = ptempalw * (1.0D0 - tmp_cc)
      tmp_err = tcpmax * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'deg C'

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         ptempalw = ptempalw + 273.15D0
         tcpmax = tcpmax + 273.15D0
      end if

   end subroutine constraint_eqn_044

   subroutine constraint_eqn_045(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for edge safety factor lower limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for edge safety factor lower limit (TART)
      !! #=# tfcoil
      !! #=#=# fq, q95_min
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fq : input real : f-value for edge safety factor
      !! q95 : safety factor 'near' plasma edge
      !! (unless i_plasma_current = 2 (ST current scaling), in which case q = mean edge safety factor qbar)
      !! q95_min : input real :  lower limit for edge safety factor
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use constraint_variables, only: fq
      use physics_variables, only: q95, q95_min, itart
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(9)
      tmp_cc =   1.0D0 - fq * q95/q95_min
      tmp_con = q95_min * (1.0D0 - tmp_cc)
      tmp_err = q95_min * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_045

   subroutine constraint_eqn_046(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Ip/Irod upper limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Ip/Irod upper limit (TART)
      !! #=# tfcoil
      !! #=#=# fipir, cratmx
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! eps : input real :  inverse aspect ratio
      !! fipir : input real : f-value for Ip/Irod upper limit
      !! c_tf_total : input real : total (summed) current in TF coils (A)
      !! plasma_current : input real :  plasma current (A)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use physics_variables, only: eps, plasma_current, itart
      use constraint_variables, only: fipir
      use tfcoil_variables, only: c_tf_total
      implicit none
      ! cratmx : local real : maximum ratio of plasma current to centrepost current
      real(dp) :: cratmx
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(10)
      cratmx = 1.0D0 + 4.91D0*(eps-0.62D0)
      tmp_cc =  (plasma_current / c_tf_total) / cratmx - 1.0D0 * fipir
      tmp_con = cratmx * (1.0D0 - tmp_cc)
      tmp_err = plasma_current/c_tf_total * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_046

   subroutine constraint_eqn_047(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Issue #508 Remove RFP option: Relevant only to reversed field pinch devices
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Issue #508 Remove RFP option: Relevant only to reversed field pinch devices
      !! Equation for TF coil toroidal thickness upper limit
      !! #=# empty
      !! #=#=# empty
      implicit none
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_047

   subroutine constraint_eqn_048(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for poloidal beta upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for poloidal beta upper limit
      !! #=# physics
      !! #=#=# fbeta_poloidal, beta_poloidal_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbeta_poloidal : input real : rf-value for poloidal beta
      !! beta_poloidal_max : input real :  maximum poloidal beta
      !! beta_poloidal : input real :  poloidal beta
      use constraint_variables, only: fbeta_poloidal, beta_poloidal_max
      use physics_variables, only: beta_poloidal
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  beta_poloidal/beta_poloidal_max - 1.0D0 * fbeta_poloidal
      tmp_con = beta_poloidal_max * (1.0D0 - tmp_cc)
      tmp_err = beta_poloidal * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_048

   subroutine constraint_eqn_049(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Issue #508 Remove IFE option: Equation for repetition rate upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Issue #508 Remove IFE option: Equation for repetition rate upper limit
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_049

   subroutine constraint_eqn_050(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! IFE option: Equation for repetition rate upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! author: S I Muldrew, CCFE, Culham Science Centre
      !! IFE option: Equation for repetition rate upper limit
      !! #=# IFE
      !! #=#=# frrmax, rrmax
      use ife_variables, only: frrmax, ife, rrmax, reprat
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      if (ife /= 1) then
         call report_error(12)
      end if

      tmp_cc =  reprat/rrmax - 1.0D0 * frrmax
      tmp_con = rrmax * (1.0D0 - tmp_cc)
      tmp_err = reprat * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'Hz'

   end subroutine constraint_eqn_050

   subroutine constraint_eqn_052(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for tritium breeding ratio lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for tritium breeding ratio lower limit
      !! #=# fwbs
      !! #=#=# ftbr, tbrmin
      !! ? TODO should this only be for certain blanket models ?
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftbr : input real : f-value for minimum tritium breeding ratio
      !! tbr : input real :  tritium breeding ratio (i_blanket_type=2,3 (KIT HCPB/HCLL))
      !! tbrmin : input real :  minimum tritium breeding ratio (If i_blanket_type=1, tbrmin=minimum 5-year time-averaged tritium breeding ratio)
      use constraint_variables, only: ftbr, tbrmin
      use fwbs_variables, only: tbr
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftbr * tbr/tbrmin
      tmp_con = tbrmin * (1.0D0 - tmp_cc)
      tmp_err = tbrmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_052

   subroutine constraint_eqn_053(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fast neutron fluence on TF coil upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fast neutron fluence on TF coil upper limit
      !! #=# fwbs
      !! #=#=# fflutf, nflutfmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fflutf : input real : f-value for maximum TF coil nuclear heating
      !! nflutfmax : input real :  max fast neutron fluence on TF coil (n/m2)
      !! nflutf : input real :  peak fast neutron fluence on TF coil superconductor (n/m2)
      use constraint_variables, only: fflutf, nflutfmax
      use fwbs_variables, only: nflutf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  nflutf/nflutfmax - 1.0D0 * fflutf
      tmp_con = nflutfmax * (1.0D0 - tmp_cc)
      tmp_err = nflutf * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'neutron/m2'

   end subroutine constraint_eqn_053

   subroutine constraint_eqn_054(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for peak TF coil nuclear heating upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for peak TF coil nuclear heating upper limit
      !! #=# fwbs
      !! #=#=# fptfnuc, ptfnucmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fptfnuc : input real : f-value for maximum TF coil nuclear heating
      !! ptfnucmax : input real :  maximum nuclear heating in TF coil (MW/m3)
      !! ptfnucpm3 : input real :  nuclear heating in the TF coil (MW/m3) (blktmodel>0)
      use constraint_variables, only: fptfnuc, ptfnucmax
      use fwbs_variables, only: ptfnucpm3
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = ptfnucpm3/ptfnucmax - 1.0D0 * fptfnuc
      tmp_con = ptfnucmax * (1.0D0 - tmp_cc)
      tmp_err = ptfnucpm3 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m3'

   end subroutine constraint_eqn_054

   subroutine constraint_eqn_055(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! vvhemax is no longer calculated in PROCESS and this constraint is disabled
      implicit none

      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      call report_error(173)
   end subroutine constraint_eqn_055

   subroutine constraint_eqn_056(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for power through separatrix / major radius upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for power through separatrix / major radius upper limit
      !! #=# current_drive
      !! #=#=# fnbshinef, f_p_beam_shine_through_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpsepr : input real : f-value for maximum Psep/R limit
      !! pseprmax : input real :  maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
      !! p_plasma_separatrix_mw : input real :  power to be conducted to the divertor region (MW)
      !! rmajor : input real :  plasma major radius (m)
      use constraint_variables, only: fpsepr, pseprmax
      use physics_variables, only: p_plasma_separatrix_mw, rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = (p_plasma_separatrix_mw/rmajor)/pseprmax - 1.0D0 * fpsepr
      tmp_con = pseprmax * (1.0D0 - tmp_cc)
      tmp_err = (p_plasma_separatrix_mw/rmajor) * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m'

   end subroutine constraint_eqn_056

   subroutine constraint_eqn_057(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_057

   subroutine constraint_eqn_058(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_058

   subroutine constraint_eqn_059(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutral beam shine-through fraction upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutral beam shine-through fraction upper limit
      !! #=# current_drive
      !! #=#=# fnbshinef, f_p_beam_shine_through_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fnbshinef : input real : f-value for maximum neutral beam shine-through fraction
      !! f_p_beam_shine_through_max : input real :  maximum neutral beam shine-through fraction
      !! f_p_beam_shine_through : input real :  neutral beam shine-through fraction
      use constraint_variables, only: fnbshinef, f_p_beam_shine_through_max
      use current_drive_variables, only: f_p_beam_shine_through
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc = f_p_beam_shine_through/f_p_beam_shine_through_max - 1.0D0 * fnbshinef
      tmp_con = f_p_beam_shine_through_max * (1.0D0 - tmp_cc)
      tmp_err = f_p_beam_shine_through * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''
   end subroutine constraint_eqn_059

   subroutine constraint_eqn_060(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid s/c temperature margin lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid s/c temperature margin lower limit
      !! #=# tfcoil
      !! #=#=# ftmargoh, tmargmin_cs
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftmargoh : input real :  f-value for central solenoid temperature margin
      !! temp_cs_margin : input real :  Central solenoid temperature margin (K)
      !! tmargmin_cs : input real :  Minimum allowable temperature margin : CS (K)
      use constraint_variables, only: ftmargoh
      use pfcoil_variables, only: temp_cs_margin
      use tfcoil_variables, only: tmargmin_cs
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - ftmargoh * temp_cs_margin/tmargmin_cs
      tmp_con = tmargmin_cs
      tmp_err = tmargmin_cs - temp_cs_margin
      tmp_symbol = '>'
      tmp_units = 'K'

   end subroutine constraint_eqn_060

   subroutine constraint_eqn_061(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for availability lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for availability limit
      !! #=# cost
      !! #=#=# favail, avail_min
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! favail : input real : F-value for minimum availability
      !! cfactr : input real : Total plant availability fraction
      !! avail_min : input real : Minimum availability
      use cost_variables, only: favail, cfactr, avail_min
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - favail * cfactr / avail_min
      tmp_con = avail_min * (1.0D0 - tmp_cc)
      tmp_err = cfactr * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_061

   subroutine constraint_eqn_062(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
      !! #=# physics
      !! #=#=# falpha_energy_confinement, f_alpha_energy_confinement_min
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! falpha_energy_confinement : input real : f-value for lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement
      !! t_alpha_confinement : input real : alpha particle confinement time (s)
      !! t_energy_confinement : input real : global thermal energy confinement time (sec)
      !! f_alpha_energy_confinement_min : input real : Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
      use constraint_variables, only: falpha_energy_confinement, f_alpha_energy_confinement_min
      use physics_variables, only: t_alpha_confinement, t_energy_confinement
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - falpha_energy_confinement * (t_alpha_confinement / t_energy_confinement) / f_alpha_energy_confinement_min
      tmp_con = f_alpha_energy_confinement_min
      tmp_err = (t_alpha_confinement / t_energy_confinement) * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_062

   subroutine constraint_eqn_063(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on niterpump (vacuum_model = simple)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Upper limit on niterpump (vacuum_model = simple)
      !! #=# vacuum
      !! #=#=# fniterpump, tfno
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fniterpump : input real : f-value for constraint that number of pumps < tfno
      !! tfno : input real : number of TF coils (default = 50 for stellarators)
      !! niterpump : input real : number of high vacuum pumps (real number), each with the throughput
      use constraint_variables, only: fniterpump
      use tfcoil_variables, only: n_tf_coils
      use vacuum_variables, only: niterpump
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = niterpump/n_tf_coils - 1.0D0 * fniterpump
      tmp_con = n_tf_coils
      tmp_err = n_tf_coils * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_063

   subroutine constraint_eqn_064(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on Zeff
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Upper limit on Zeff
      !! #=# physics
      !! #=#=# fzeffmax, zeffmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fzeffmax : input real : f-value for maximum zeff
      !! zeffmax : input real : maximum value for Zeff
      !! zeff : input real : plasma effective charge
      use constraint_variables, only: fzeffmax, zeffmax
      use physics_variables, only: zeff
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = zeff/zeffmax - 1.0D0 * fzeffmax
      tmp_con = zeffmax
      tmp_err = zeffmax * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_064

   subroutine constraint_eqn_065(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on stress of the vacuum vessel that occurs when the TF coil quenches.
      !! author: Timothy Nunn, UKAEA
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fmaxvvstress : input real : f-value for constraint on maximum VV stress
      !! max_vv_stress : input real : Maximum permitted stress of the VV (Pa)
      !! vv_stress_quench : input real : Stress of the VV (Pa)
      use constraint_variables, only: fmaxvvstress
      use tfcoil_variables, only: max_vv_stress
      use sctfcoil_module, only: vv_stress_quench
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  vv_stress_quench / max_vv_stress - 1.0d0 * fmaxvvstress
      tmp_con = max_vv_stress
      tmp_err = max_vv_stress * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_065

   subroutine constraint_eqn_066(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on rate of change of energy in poloidal field
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Limit on rate of change of energy in poloidal field
      !! #=# pfcoil
      !! #=#=# fpoloidalpower, maxpoloidalpower
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpoloidalpower : input real : f-value for constraint on rate of change of energy in poloidal field
      !! maxpoloidalpower : input real : Maximum permitted absolute rate of change of stored energy in poloidal field (MW)
      !! peakpoloidalpower : input real : Peak absolute rate of change of stored energy in poloidal field (MW) (11/01/16)
      use constraint_variables, only: fpoloidalpower
      use pf_power_variables, only: maxpoloidalpower, peakpoloidalpower
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = peakpoloidalpower / maxpoloidalpower - 1.0d0 * fpoloidalpower
      tmp_con = maxpoloidalpower
      tmp_err = maxpoloidalpower * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_066

   subroutine constraint_eqn_067(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Simple upper limit on radiation wall load
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Simple upper limit on radiation wall load
      !! #=# physics
      !! #=#=# fradwall, pflux_fw_rad_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fradwall : input real : f-value for upper limit on radiation wall load
      !! pflux_fw_rad_max : input real : Maximum permitted radiation wall load (MW/m^2)
      !! pflux_fw_rad_max_mw : input real : Peak radiation wall load (MW/m^2)
      use constraint_variables, only: fradwall, pflux_fw_rad_max, pflux_fw_rad_max_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = pflux_fw_rad_max_mw / pflux_fw_rad_max - 1.0d0 * fradwall
      tmp_con = pflux_fw_rad_max
      tmp_err =  pflux_fw_rad_max * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m^2'

   end subroutine constraint_eqn_067

   subroutine constraint_eqn_068(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on Psep scaling (PsepB/qAR)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! New Psep scaling (PsepB/qAR)
      !! Issue #464
      !! #=# physics
      !! #=#=# fpsepbqar, psepbqarmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpsepbqar : input real : f-value for upper limit on psepbqar, maximum Psep*Bt/qAR limit
      !! psepbqarmax : input real : maximum permitted value of ratio of Psep*Bt/qAR (MWT/m)
      !! p_plasma_separatrix_mw : input real : Power to conducted to the divertor region (MW)
      !! bt : input real : toroidal field on axis (T) (iteration variable 2)
      !! q95 : input real : safety factor q at 95% flux surface
      !! aspect : input real : aspect ratio (iteration variable 1)
      !! rmajor : input real : plasma major radius (m) (iteration variable 3)
      !! i_q95_fixed : input int : Switch that allows for fixing q95 only in this constraint.
      !! q95_fixed : input real : fixed safety factor q at 95% flux surface
      use constraint_variables, only: fpsepbqar, psepbqarmax, i_q95_fixed, q95_fixed
      use physics_variables, only: p_plasma_separatrix_mw, bt, q95, aspect, rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      if (i_q95_fixed == 1) then
         tmp_cc = ((p_plasma_separatrix_mw*bt)/(q95_fixed*aspect*rmajor)) / psepbqarmax - 1.0d0 * fpsepbqar
         tmp_err = (p_plasma_separatrix_mw*bt)/(q95_fixed*aspect*rmajor) - psepbqarmax
      else
         tmp_cc = ((p_plasma_separatrix_mw*bt)/(q95*aspect*rmajor)) / psepbqarmax - 1.0d0 * fpsepbqar
         tmp_err = (p_plasma_separatrix_mw*bt)/(q95*aspect*rmajor) - psepbqarmax
      end if
      tmp_con = psepbqarmax
      tmp_symbol = '<'
      tmp_units = 'MWT/m'

   end subroutine constraint_eqn_068

   subroutine constraint_eqn_069(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Ensure separatrix power is less than value from Kallenbach divertor
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure separatrix power is less than value from Kallenbach divertor
      !! #=# divertor_kallenbach
      !! #=#=# consistency, psep_kallenbach
      !! fpsep has been removed from the equation.
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! psep_kallenbach : input real : Power conducted through the separatrix, as calculated by the divertor model [W]
      !! p_plasma_separatrix_mw : input real :  power to conducted to the divertor region (MW)
      ! use div_kal_vars, only: psep_kallenbach
      use physics_variables, only: p_plasma_separatrix_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_069

   subroutine constraint_eqn_070(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Separatrix density consistency
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Separatrix density consistency
      !! #=# divertor_kallenbach
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! teomp : input real : Separatrix temperature calculated by the Kallenbach divertor model [eV]
      !! tesep : input real : Electron temperature at separatrix [keV]
      ! use div_kal_vars, only: teomp
      use  physics_variables, only: tesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_070

   subroutine constraint_eqn_071(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! TODO Remove
      !! Separatrix density consistency
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Separatrix density consistency
      !! #=# divertor_kallenbach
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! neomp : input real : Mean SOL density at OMP calculated by the Kallenbach divertor model [m-3]
      !! nesep : input real :  electron density at separatrix [m-3] (ipedestal=1,2, calculated if 3)
      !! neratio : input real : Ratio of mean SOL density at OMP to separatrix density at OMP (iteration variable 121)
      ! use div_kal_vars, only: neomp, neratio
      use physics_variables, only: nesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_071

   subroutine constraint_eqn_072(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on central Solenoid Tresca yield stress
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Central Solenoid Tresca yield criterion
      !! #=# pfcoil
      !! #=#=# foh_stress, alstroh
      !! In the case if the bucked and wedged option ( i_tf_bucking >= 2 ) the constrained
      !! stress is the largest the largest stress of the
      !!  - CS stress at maximum current (conservative as the TF inward pressure is not taken
      !!    into account)
      !!  - CS stress at flux swing (no current in CS) from the TF inward pressure
      !! This allow to cover the 2 worst stress scenario in the bucked and wedged design
      !! Otherwise (free standing TF), the stress limits are only set by the CS stress at max current
      !! Reverse the sign so it works as an inequality constraint (tmp_cc > 0)
      !! This will have no effect if it is used as an equality constraint because it will be squared.
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! foh_stress : input real : f-value for Tresca yield criterion in Central Solenoid
      !! alstroh : input real :  allowable hoop stress in Central Solenoid structural material (Pa)
      !! s_shear_cs_peak : input real : Maximum shear stress coils/central solenoid (Pa)
      !! sig_tf_cs_bucked : input real : Maximum shear stress in CS case at flux swing (no current in CS)
      !!                       can be significant for the bucked and weged design
      !! i_tf_bucking : input integer : switch for TF structure design
      use constraint_variables, only: foh_stress
      use pfcoil_variables, only: alstroh, s_shear_cs_peak
      use tfcoil_variables, only: sig_tf_cs_bucked, i_tf_bucking
      use build_variables, only: i_tf_inside_cs
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! bucked and wedged desing (see subroutine comment)
      if ( i_tf_bucking >= 2 .and. i_tf_inside_cs == 0 ) then
         tmp_cc = max(s_shear_cs_peak, sig_tf_cs_bucked) / alstroh - 1.0d0 * foh_stress
         tmp_err = alstroh - max(s_shear_cs_peak, sig_tf_cs_bucked)
      ! Free standing CS
      else
         tmp_cc = s_shear_cs_peak / alstroh - 1.0d0 * foh_stress
         tmp_err = alstroh - s_shear_cs_peak
      end if

      tmp_con = alstroh
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_072

   subroutine constraint_eqn_073(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit to ensure separatrix power is greater than the L-H power + auxiliary power
      !! Related to constraint 15
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure separatrix power is greater than the L-H power + auxiliary power
      !! #=# physics
      !! #=#=# fplhsep, p_plasma_separatrix_mw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fplhsep : input real : F-value for Psep >= Plh + Paux : for consistency of two values of separatrix power
      !! p_l_h_threshold_mw : input real : L-H mode power threshold (MW)
      !! p_plasma_separatrix_mw : input real : power to be conducted to the divertor region (MW)
      !! p_hcd_injected_total_mw : inout real : total auxiliary injected power (MW)
      use physics_variables, only: fplhsep, p_l_h_threshold_mw, p_plasma_separatrix_mw
      use current_drive_variables, only: p_hcd_injected_total_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fplhsep * p_plasma_separatrix_mw / (p_l_h_threshold_mw+p_hcd_injected_total_mw)
      tmp_con = p_plasma_separatrix_mw
      tmp_err = p_plasma_separatrix_mw * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_073

   subroutine constraint_eqn_074(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit to ensure TF coil quench temperature < tmax_croco
      !! ONLY used for croco HTS coil
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure TF coil quench temperature < tmax_croco ONLY used for croco HTS coil
      !! #=# physics
      !! #=#=# fcqt, tmax_croco
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fcqt : input real : f-value: TF coil quench temparature remains below tmax_croco
      !! croco_quench_temperature : input real : CroCo strand: Actual temp reached during a quench (K)
      !! tmax_croco : input real : CroCo strand: maximum permitted temp during a quench (K)
      use constraint_variables, only: fcqt
      use tfcoil_variables, only: croco_quench_temperature, tmax_croco
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = croco_quench_temperature / tmax_croco - 1.0d0 * fcqt
      tmp_con = croco_quench_temperature
      tmp_err = croco_quench_temperature * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'K'

   end subroutine constraint_eqn_074

   subroutine constraint_eqn_075(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit to ensure that TF coil current / copper area < Maximum value
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure that TF coil current / copper area < Maximum value
      !! ONLY used for croco HTS coil
      !! #=# physics
      !! #=#=# f_coppera_m2, copperA_m2_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! copperA_m2 : input real : TF coil current / copper area (A/m2)
      !! copperA_m2_max : input real : Maximum TF coil current / copper area (A/m2)
      !! f_coppera_m2 : input real : f-value for TF coil current / copper area < copperA_m2_max
      use rebco_variables, only: copperA_m2, copperA_m2_max, f_coppera_m2
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = copperA_m2 / copperA_m2_max - 1.0d0 * f_coppera_m2
      tmp_con = copperA_m2
      tmp_err = copperA_m2 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_075

   subroutine constraint_eqn_076(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit for Eich critical separatrix density model: Added for issue 558
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Eich critical separatrix density model
      !! Added for issue 558 with ref to http://iopscience.iop.org/article/10.1088/1741-4326/aaa340/pdf
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! alpha_crit : output real : critical ballooning parameter value
      !! nesep_crit : output real : critical electron density at separatrix [m-3]
      !! kappa : input real : plasma separatrix elongation (calculated if i_plasma_geometry = 1-5, 7 or 9)
      !! triang : input real : plasma separatrix triangularity (calculated if i_plasma_geometry = 1, 3-5 or 7)
      !! aspect : input real : aspect ratio (iteration variable 1)
      !! p_plasma_separatrix_mw : input real : power to conducted to the divertor region (MW)
      !! dlimit(7) : input real array : density limit (/m3) as calculated using various models
      !! fnesep : input real : f-value for Eich critical separatrix density
      use physics_variables, only: alpha_crit, nesep_crit, kappa, triang, &
                                   aspect, p_plasma_separatrix_mw, dlimit, nesep
      use constraint_variables, only: fnesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      alpha_crit = (kappa ** 1.2D0) * (1.0D0 + 1.5D0 * triang)
      nesep_crit = 5.9D0 * alpha_crit * (aspect ** (-2.0D0/7.0D0)) * &
                (((1.0D0 + (kappa ** 2.0D0)) / 2.0D0) ** (-6.0D0/7.0D0)) &
                * ((p_plasma_separatrix_mw* 1.0D6) ** (-11.0D0/70.0D0)) * dlimit(7)
      tmp_cc = nesep / nesep_crit - 1.0D0 * fnesep
      tmp_con = nesep
      tmp_err = nesep * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm-3'

   end subroutine constraint_eqn_076

   subroutine constraint_eqn_077(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for maximum TF current per turn upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for maximum TF current per turn upper limit
      !! #=# tfcoil
      !! #=#=# fcpttf, cpttf, cpttf_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fcpttf : input : f-value for TF coil current per turn
      !! cpttf_max  : input : allowable TF coil current per turn [A/turn]
      !! cpttf  : input : TF coil current per turn [A/turn]
      use constraint_variables, only: fcpttf
      use tfcoil_variables, only: cpttf_max, cpttf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = cpttf / cpttf_max - 1.0D0 * fcpttf
      tmp_con = cpttf_max
      tmp_err = cpttf_max * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/turn'

   end subroutine constraint_eqn_077

   subroutine constraint_eqn_078(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Reinke criterion, divertor impurity fraction lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for Reinke criterion, divertor impurity fraction lower limit
      !! #=# divertor
      !! #=#=# freinke, fzactual, fzmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present;
      !! and con will be printed out only if present. Thesw conditions were missing.
      !! freinke : input : f-value for Reinke criterion (itv 147)
      !! fzmin : input : minimum impurity fraction from Reinke model
      !! fzactual : input : actual impurity fraction
      use constraint_variables, only: freinke
      use reinke_variables, only: fzactual, fzmin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! write(*,*) 'freinke, fzact, fzmin = ', freinke, ', ', fzactual, ', ', fzmin
      !            1.0,    0.0,   value
      tmp_cc = 1.0D0 - freinke *  fzactual/fzmin
      !The following two pre-existing lines are not understood:
      !KE note - cc is always 1, code never enters IF statement...
      tmp_con = fzmin * (1.0D0 - tmp_cc)
      tmp_err = fzmin * tmp_cc
      tmp_symbol = '>'
      tmp_units  = ''
      ! write(*,*) 'cc, con = ', tmp_cc, ', ', tmp_con
      ! write(*,*) 'freinke, fzactual, fzmin = ', freinke, ', ', fzactual, ', ', fzmin

   end subroutine constraint_eqn_078

   subroutine constraint_eqn_079(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for maximum CS field
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for maximum CS field
      !! #=# pfcoil
      !! #=#=# fb_cs_limit_max, b_cs_peak_flat_top_end, b_cs_peak_pulse_start, b_cs_limit_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fb_cs_limit_max : input : F-value for CS mmax field (cons. 79, itvar 149)
      !! b_cs_limit_max : input : Central solenoid max field limit [T]
      !! b_cs_peak_pulse_start : input : maximum field in central solenoid at beginning of pulse (T)
      !! b_cs_peak_flat_top_end : input real : maximum field in central solenoid at end of flat-top (EoF) (T)
      !! (Note: original code has "b_cs_peak_flat_top_end/b_cs_peak_pulse_start |  peak CS field [T]".)
      use pfcoil_variables, only: fb_cs_limit_max, b_cs_limit_max, b_cs_peak_pulse_start, b_cs_peak_flat_top_end
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc     = max(b_cs_peak_flat_top_end, b_cs_peak_pulse_start) / b_cs_limit_max - 1.0D0 * fb_cs_limit_max
      tmp_con    = b_cs_limit_max
      tmp_err    = max(b_cs_peak_flat_top_end, b_cs_peak_pulse_start) * tmp_cc
      tmp_symbol = '<'
      tmp_units  = 'A/turn'

   end subroutine constraint_eqn_079

   subroutine constraint_eqn_080(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for p_plasma_separatrix_mw lower limit
      !! author: J Morris, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units;
      !! output string; units string
      !! Lower limit p_plasma_separatrix_mw
      !! #=# physics
      !! #=#=# fp_plasma_separatrix_min_mw, p_plasma_separatrix_mw
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fp_plasma_separatrix_min_mw : input : F-value for lower limit on p_plasma_separatrix_mw (cons. 80, itvar 153)
      !! p_plasma_separatrix_min_mw : input : Minimum power crossing separatrix p_plasma_separatrix_mw [MW]
      !! p_plasma_separatrix_mw : input : Power crossing separatrix [MW]
      use physics_variables, only: fp_plasma_separatrix_min_mw, p_plasma_separatrix_mw
      use constraint_variables, only : p_plasma_separatrix_min_mw
      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc     = 1.0D0 - fp_plasma_separatrix_min_mw * p_plasma_separatrix_mw / p_plasma_separatrix_min_mw
      tmp_con    = p_plasma_separatrix_min_mw
      tmp_err    = p_plasma_separatrix_mw * tmp_cc
      tmp_symbol = '>'
      tmp_units  = 'MW'

   end subroutine constraint_eqn_080

   subroutine constraint_eqn_081(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit to ensure central density is larger that the pedestal one
      !! author: S Kahn, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Lower limit ne0 > neped
      !! !#=# physics
      !! !#=#=# ne0, neped
      !! Logic change during pre-factoring: err, symbol, units will be
      !! assigned only if present.
      !! fne0  : input : F-value for constraint on ne0 > neped
      !! ne0   : input : Central electron density [m-3]
      !! neped : input : Electron density at pedestal [m-3]
      use physics_variables, only: ne0, fne0, neped
      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc     = 1.0D0 - fne0 * ne0/neped
      tmp_con    = fne0
      tmp_err    = fne0 * tmp_cc
      tmp_symbol = '>'
      tmp_units  = '/m3'

   end subroutine constraint_eqn_081

   subroutine constraint_eqn_082(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for toroidal consistency of stellarator build
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! toroidalgap > dx_tf_inboard_out_toroidal
      !! #=# tfcoil
      !! #=#=# dx_tf_inboard_out_toroidal, ftoroidalgap
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftoroidalgap : input real : f-value for constraint toroidalgap > dx_tf_inboard_out_toroidal
      !! toroidalgap : input real :  minimal gap between two stellarator coils
      !! dx_tf_inboard_out_toroidal : input real :  total toroidal width of a tf coil
      use tfcoil_variables, only: dx_tf_inboard_out_toroidal,ftoroidalgap,toroidalgap
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftoroidalgap * toroidalgap/dx_tf_inboard_out_toroidal
      tmp_con = toroidalgap
      tmp_err = toroidalgap - dx_tf_inboard_out_toroidal/ftoroidalgap
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_082

   subroutine constraint_eqn_083(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radial consistency of stellarator build
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! available_radial_space > required_radial_space
      !! #=# build
      !! #=#=# required_radial_space, f_avspace
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! f_avspace : input real : f-value for constraint available_radial_space > required_radial_space
      !! available_radial_space : input real :  avaible space in radial direction as given by each s.-configuration
      !! required_radial_space : input real :  required space in radial direction
      use build_variables, only: available_radial_space, required_radial_space, f_avspace
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - f_avspace  * available_radial_space/required_radial_space
      tmp_con = available_radial_space * (1.0D0 - tmp_cc)
      tmp_err = required_radial_space * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'
   end subroutine constraint_eqn_083

   subroutine constraint_eqn_084(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for the lower limit of beta
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !!  (beta-beta_fast_alpha) > beta_min
      !! #=# physics
      !! #=#=# beta_fast_alpha, beta, fbeta_min
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbeta_min : input real : f-value for constraint beta-beta_fast_alpha > beta_min
      !! beta_min : input real :  Lower limit for beta
      !! beta : input real :  plasma beta

      use physics_variables, only: beta_min, beta
      use constraint_variables, only: fbeta_min
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units


      tmp_cc = 1.0D0 - fbeta_min * (beta)/beta_min
      tmp_con = beta_min * (1.0D0 - tmp_cc)
      tmp_err = (beta) * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''


   end subroutine constraint_eqn_084

   subroutine constraint_eqn_086(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on the turn edge length in the TF winding pack
      !! Author : S Kahn
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units;
      !! t_turn_tf : input real : TF coil turn edge length including turn insulation [m]
      !! f_t_turn_tf : input real : f-value for TF turn edge length constraint
      !! t_turn_tf_max : input real : TF turn edge length including turn insulation upper limit [m]
      use tfcoil_variables, only : t_turn_tf, f_t_turn_tf, t_turn_tf_max

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      !! Constraints output
      tmp_cc = t_turn_tf / t_turn_tf_max - 1.0D0 * f_t_turn_tf
      tmp_con = t_turn_tf_max * (1.0D0 - tmp_cc)
      tmp_err = t_turn_tf_max * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_086


   subroutine constraint_eqn_087(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil cryogenic power upper limit
      !! author: S. Kahn, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! crypmw : input real : cryogenic plant power (MW)
      !! f_crypmw : input real : f-value for maximum cryogenic plant power
      !! crypmw_max : input real : Maximum cryogenic plant power (MW)

      use heat_transport_variables, only: crypmw, crypmw_max, f_crypmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  crypmw / crypmw_max - 1.0D0 * f_crypmw
      tmp_con = crypmw_max * (1.0D0 - tmp_cc)
      tmp_err = crypmw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'
   end subroutine constraint_eqn_087

   subroutine constraint_eqn_088(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil vertical strain upper limit (absolute value)
      !! author: CPS Swanson, PPPL, USA
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil vertical strain upper limit (absolute value)
      !! #=# tfcoil
      !! #=#=# fstr_wp, str_wp_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstr_wp : input real : f-value for TF coil strain
      !! str_wp_max : input real : Allowable maximum TF coil vertical strain
      !! str_wp : input real : Constrained TF coil vertical strain
      use constraint_variables, only: fstr_wp
      use tfcoil_variables, only: str_wp_max, str_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  abs(str_wp) / str_wp_max - 1.0D0 * fstr_wp
      tmp_con = str_wp_max
      tmp_err = str_wp_max - abs(str_wp) / fstr_wp
      tmp_symbol = '<'
      tmp_units = ''
   end subroutine constraint_eqn_088

   subroutine constraint_eqn_089(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit to ensure that the Central Solenoid [OH] coil current / copper area < Maximum value
      !! author: G Turkington, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! #=# physics
      !! #=#=# f_copperaoh_m2, copperaoh_m2_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! copperaoh_m2 : input real : CS coil current at EOF / copper area [A/m2]
      !! copperaoh_m2_max : input real : maximum coil current / copper area [A/m2]
      !! f_copperaoh_m2 : input real : f-value for CS coil current / copper area
      use rebco_variables, only: copperaoh_m2, copperaoh_m2_max, f_copperaoh_m2
      implicit none
                  real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = copperaoh_m2 / copperaoh_m2_max - 1.0d0 * f_copperaoh_m2
      tmp_con = copperaoh_m2
      tmp_err = copperaoh_m2 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_089

   subroutine constraint_eqn_090(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit for CS coil stress load cycles
      !! author: A. Pearce, G Turkington CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! fncycle : input real : f-value for constraint n_cycle > n_cycle_min
      !! n_cycle : input real : Allowable number of cycles for CS
      !! n_cycle_min : input real :  Minimum required cycles for CS
      use CS_fatigue_variables, only: n_cycle, n_cycle_min, bkt_life_csf
      use constraint_variables, only: fncycle
      use cost_variables, only: ibkt_life, bktcycles
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Switch to relay the calculated fw/blanket lifetime cycles as the minimum required CS stress cycles.
      !! bkt_life_cycle = 1 turns on the relay. Otherwise the models run independently.
      if (ibkt_life == 1 .and. bkt_life_csf == 1 ) then
         n_cycle_min = bktcycles
      end if

      tmp_cc =  1.0D0 - fncycle * n_cycle / n_cycle_min
      tmp_con = n_cycle_min * (1.0D0 - tmp_cc)
      tmp_err = n_cycle * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_090

   subroutine constraint_eqn_091(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit to ensure ECRH te is greater than required te for ignition
      !! at lower values for n and B. Or if the design point is ECRH heatable (if ignite==0)
      !! stellarators only (but in principle usable also for tokamaks).
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !!  powerht_local > powerscaling
      !! #=# physics
      !! #=#=# fecrh_ignition, powerht_local, powerscaling
      !! fecrh_ignition : input real : f-value for constraint powerht_local > powerscaling
      !! max_gyrotron_frequency : input real :  Max. av. gyrotron frequency
      !! te0_ecrh_achievable : input real : Max. achievable electron temperature at ignition point
      use constraint_variables, only: fecrh_ignition
      use stellarator_variables, only: max_gyrotron_frequency, te0_ecrh_achievable, powerscaling_constraint, powerht_constraint
      use physics_variables, only: ignite
      use current_drive_variables, only: p_hcd_primary_extra_heat_mw
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! Achievable ECRH te needs to be larger than needed te for igntion
      if(ignite==0) then
         tmp_cc = 1.0D0 - fecrh_ignition* (powerht_constraint+p_hcd_primary_extra_heat_mw)/powerscaling_constraint
      else
         tmp_cc = 1.0D0 - fecrh_ignition* powerht_constraint/powerscaling_constraint
      endif

      tmp_con = powerscaling_constraint * (1.0D0 - tmp_cc)
      tmp_err = powerht_constraint * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'
   end subroutine constraint_eqn_091

end module constraints
