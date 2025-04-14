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
      !! i_plasma_ignited : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! Obviously, i_plasma_ignited must be zero if current drive is required.
      !! If i_plasma_ignited=1, any auxiliary power is assumed to be used only
      !! during plasma start-up, and is excluded from all steady-state
      !! power balance calculations.
      use constraint_variables, only: fqval, bigqmin
      use current_drive_variables, only: bigq
      use physics_variables, only: i_plasma_ignited
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if plasma is not ignited ...
      if (i_plasma_ignited == 0) then
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
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_045


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
      tmp_symbol = '>'
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
      tmp_symbol = '>'
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
      !! at lower values for n and B. Or if the design point is ECRH heatable (if i_plasma_ignited==0)
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
      use physics_variables, only: i_plasma_ignited
      use current_drive_variables, only: p_hcd_primary_extra_heat_mw
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! Achievable ECRH te needs to be larger than needed te for igntion
      if(i_plasma_ignited==0) then
         tmp_cc = 1.0D0 - fecrh_ignition* (powerht_constraint+p_hcd_primary_extra_heat_mw)/powerscaling_constraint
      else
         tmp_cc = 1.0D0 - fecrh_ignition* powerht_constraint/powerscaling_constraint
      endif

      tmp_con = powerscaling_constraint * (1.0D0 - tmp_cc)
      tmp_err = powerht_constraint * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'MW'
   end subroutine constraint_eqn_091

end module constraints
