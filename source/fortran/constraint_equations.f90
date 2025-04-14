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

   subroutine constraint_eqn_002(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
    !! author: J. Morris
    !! category: equality constraint
    !!
    !! Global plasma power balance equation
    !!
    !! \begin{equation} c_i =
    !! \end{equation}
    !!
    !! i_rad_loss : input integer : switch for radiation loss term usage in power balance (see User Guide):<UL>
    !! <LI> = 0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
    !! <LI> = 1 total power lost is scaling power plus core radiation only
    !! <LI> = 2 total power lost is scaling power only, with no additional
    !! allowance for radiation. This is not recommended for power plant models.</UL>
    !! i_plasma_ignited : input integer : switch for ignition assumption:<UL>
    !! <LI> = 0 do not assume plasma ignition;
    !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
    !! pden_electron_transport_loss_mw : input real : electron transport power per volume (MW/m3)
    !! pden_ion_transport_loss_mw : input real :  ion transport power per volume (MW/m3)
    !! pden_plasma_rad_mw : input real : total radiation power per volume (MW/m3)
    !! pden_plasma_core_rad_mw : input real : total core radiation power per volume (MW/m3)
    !! f_alpha_plasma : input real : fraction of alpha power deposited in plasma
    !! alpha_power_density_total : input real : alpha power per volume (MW/m3)
    !! charged_power_density : input real : non-alpha charged particle fusion power per volume (MW/m3)
    !! pden_plasma_ohmic_mw : input real : ohmic heating power per volume (MW/m3)
    !! p_hcd_injected_total_mw : input real : total auxiliary injected power (MW)
    !! vol_plasma : input real : plasma volume (m3)

    use physics_variables, only: i_rad_loss, i_plasma_ignited, pden_electron_transport_loss_mw, pden_ion_transport_loss_mw, pden_plasma_rad_mw, &
                                  pden_plasma_core_rad_mw, f_alpha_plasma, alpha_power_density_total, charged_power_density, &
                                  pden_plasma_ohmic_mw, vol_plasma
    use current_drive_variables, only: p_hcd_injected_total_mw

    implicit none

          real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
    !! constraint derived type

    ! pscaling : Local real : total transport power per volume (MW/m3)
    real(dp) :: pscaling
    real(dp) :: pnumerator, pdenom
    pscaling = pden_electron_transport_loss_mw + pden_ion_transport_loss_mw
    ! Total power lost is scaling power plus radiation:
    if (i_rad_loss == 0) then
        pnumerator = pscaling + pden_plasma_rad_mw
    else if (i_rad_loss == 1) then
        pnumerator = pscaling + pden_plasma_core_rad_mw
    else
        pnumerator = pscaling
    end if

    ! if plasma not ignited include injected power
    if (i_plasma_ignited == 0) then
      pdenom = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw + p_hcd_injected_total_mw/vol_plasma
    else
      ! if plasma ignited
      pdenom = f_alpha_plasma*alpha_power_density_total + charged_power_density + pden_plasma_ohmic_mw
    end if

    tmp_cc = 1.0D0 - pnumerator / pdenom
    tmp_con = pdenom * (1.0D0 - tmp_cc)
    tmp_err = pdenom * tmp_cc
    tmp_symbol = '='
    tmp_units = 'MW/m3'

   end subroutine constraint_eqn_002

   subroutine constraint_eqn_003(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Global power balance equation for ions
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Global power balance equation for ions
      !! This is a consistency equation (NBI)
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! i_plasma_ignited : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! pden_ion_transport_loss_mw : input real :  ion transport power per volume (MW/m3)
      !! piepv : input real : ion/electron equilibration power per volume (MW/m3)
      !! f_alpha_plasma : input real : fraction of alpha power deposited in plasma
      !! alpha_power_ions_density : input real : alpha power per volume to ions (MW/m3)
      !! p_hcd_injected_ions_mw : input real : auxiliary injected power to ions (MW)
      !! vol_plasma : input real : plasma volume (m3)
      use physics_variables, only: i_plasma_ignited, pden_ion_transport_loss_mw, piepv, f_alpha_plasma, alpha_power_ions_density, vol_plasma
      use current_drive_variables, only: p_hcd_injected_ions_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

	   ! No assume plasma ignition:
      if (i_plasma_ignited == 0) then
         tmp_cc     = 1.0D0 - (pden_ion_transport_loss_mw + piepv) / (f_alpha_plasma*alpha_power_ions_density + p_hcd_injected_ions_mw/vol_plasma)
         tmp_con    = (f_alpha_plasma*alpha_power_ions_density + p_hcd_injected_ions_mw/vol_plasma) * (1.0D0 - tmp_cc)
         tmp_err    = (f_alpha_plasma*alpha_power_ions_density + p_hcd_injected_ions_mw/vol_plasma) * tmp_cc
         tmp_symbol = '='
         tmp_units  = 'MW/m3'
	   ! Plasma ignited:
      else
         tmp_cc     = 1.0D0 - (pden_ion_transport_loss_mw+piepv) / (f_alpha_plasma*alpha_power_ions_density)
         tmp_con    = (f_alpha_plasma*alpha_power_ions_density) * (1.0D0 - tmp_cc)
         tmp_err    = (f_alpha_plasma*alpha_power_ions_density) * tmp_cc
         tmp_symbol = '='
         tmp_units  = 'MW/m3'
      end if

   end subroutine constraint_eqn_003

   subroutine constraint_eqn_004(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Global power balance equation for electrons
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Global power balance equation for electrons
      !! This is a consistency equation
      !! N.B. This constraint is currently NOT RECOMMENDED for use.
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! i_rad_loss : input integer : switch for radiation loss term usage in power balance (see User Guide):<UL>
      !! <LI> = 0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
      !! <LI> = 1 total power lost is scaling power plus core radiation only
      !! <LI> = 2 total power lost is scaling power only, with no additional
      !! allowance for radiation. This is not recommended for power plant models.</UL>
      !! i_plasma_ignited : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! pden_electron_transport_loss_mw : input real : electron transport power per volume (MW/m3)
      !! pden_plasma_rad_mw : input real : total radiation power per volume (MW/m3)
      !! pden_plasma_core_rad_mw : input real : total core radiation power per volume (MW/m3)
      !! f_alpha_plasma : input real : fraction of alpha power deposited in plasma
      !! alpha_power_electron_density : input real : alpha power per volume to electrons (MW/m3)
      !! piepv : input real : ion/electron equilibration power per volume (MW/m3)
      !! p_hcd_injected_electrons_mw : input real : auxiliary injected power to electrons (MW)
      !! vol_plasma : input real : plasma volume (m3)
      use physics_variables, only: i_rad_loss, i_plasma_ignited, pden_electron_transport_loss_mw, pden_plasma_core_rad_mw, f_alpha_plasma, &
                                 alpha_power_electron_density, piepv, vol_plasma, pden_plasma_rad_mw
      use current_drive_variables, only: p_hcd_injected_electrons_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! pscaling : Local real : total transport power per volume (MW/m3)
      real(dp) :: pscaling
      real(dp) :: pnumerator, pdenom
      pscaling = pden_electron_transport_loss_mw
	   ! Total power lost is scaling power plus radiation:
      if (i_rad_loss == 0) then
         pnumerator = pscaling + pden_plasma_rad_mw
      else if (i_rad_loss == 1) then
         pnumerator = pscaling + pden_plasma_core_rad_mw
      else
         pnumerator = pscaling
      end if

      ! if plasma not ignited include injected power
      if (i_plasma_ignited == 0) then
         pdenom = f_alpha_plasma*alpha_power_electron_density + piepv + p_hcd_injected_electrons_mw/vol_plasma
      else
      ! if plasma ignited
         pdenom = f_alpha_plasma*alpha_power_electron_density + piepv
      end if

      tmp_cc     = 1.0D0 - pnumerator / pdenom
      tmp_con    = pdenom * (1.0D0 - tmp_cc)
      tmp_err    = pdenom * tmp_cc
      tmp_symbol = '='
      tmp_units  = 'MW/m3'

   end subroutine constraint_eqn_004

   subroutine constraint_eqn_006(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for epsilon beta-poloidal upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for epsilon beta-poloidal upper limit
      !! #=# physics
      !! #=#=# fbeta_poloidal_eps, beta_poloidal_eps_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbeta_poloidal_eps : input real : f-value for epsilon beta-poloidal
      !! beta_poloidal_eps_max : input real : maximum (eps*beta_poloidal)
      !! eps : input real : inverse aspect ratio
      !! beta_poloidal : input real : poloidal beta
      use physics_variables, only: beta_poloidal_eps_max, eps, beta_poloidal
      use constraint_variables, only: fbeta_poloidal_eps
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  (eps*beta_poloidal)/beta_poloidal_eps_max - 1.0D0 * fbeta_poloidal_eps
      tmp_con = beta_poloidal_eps_max * (1.0D0 - tmp_cc)
      tmp_err = (eps*beta_poloidal) * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_006

   subroutine constraint_eqn_007(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for hot beam ion density
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for hot beam ion density
      !! This is a consistency equation (NBI)
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! i_plasma_ignited : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! Obviously, i_plasma_ignited must be zero if current drive is required.
      !! If i_plasma_ignited=1, any auxiliary power is assumed to be used only
      !! during plasma start-up, and is excluded from all steady-state
      !! power balance calculations.
      !! beam_density_out : input real :  hot beam ion density from calculation (/m3)
      !! nd_beam_ions : input real : hot beam ion density, variable (/m3)
      use physics_variables, only: i_plasma_ignited, beam_density_out, nd_beam_ions
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

	   ! Do not assume plasma ignition:
      if (i_plasma_ignited == 0) then
         tmp_cc     = 1.0D0 - beam_density_out/nd_beam_ions
         tmp_con    = nd_beam_ions * (1.0D0 - tmp_cc)
         tmp_err    = nd_beam_ions * tmp_cc
         tmp_symbol = '='
         tmp_units  = '/m3'
      else
         tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''
         call report_error(1)
      end if

   end subroutine constraint_eqn_007

   subroutine constraint_eqn_008(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutron wall load upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutron wall load upper limit
      !! #=# physics
      !! #=#=# fwalld, walalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fwalld : input real : f-value for maximum wall load
      !! walalw : input real : allowable wall-load (MW/m2)
      !! pflux_fw_neutron_mw : input real : average neutron wall load (MW/m2)
      use constraint_variables, only: fwalld, walalw
      use physics_variables, only: pflux_fw_neutron_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  pflux_fw_neutron_mw/walalw - 1.0D0 * fwalld
      tmp_con = fwalld * walalw
      tmp_err = fwalld * walalw - pflux_fw_neutron_mw
      tmp_symbol = '<'
      tmp_units = 'MW/m2'

   end subroutine constraint_eqn_008

   subroutine constraint_eqn_009(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fusion power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fusion power upper limit
      !! #=# physics
      !! #=#=# ffuspow, powfmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ffuspow : input real : f-value for maximum fusion power
      !! powfmax : input real : maximum fusion power (MW)
      !! fusion_power : input real : fusion power (MW)
      use constraint_variables, only: ffuspow, powfmax
      use physics_variables, only: fusion_power
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  fusion_power/powfmax - 1.0D0 * ffuspow
      tmp_con = powfmax * (1.0D0 - tmp_cc)
      tmp_err = fusion_power * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_009

   subroutine constraint_eqn_010(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Equation for field at TF coil
      !! This is a consistency equation
      !! (do not use for stellarators)
      !! #=# tfcoil
      !! #=#=# consistency
      !! rmajor |  plasma major radius (m)
      !! bt     |  toroidal field on axis (T)
      !! r_b_tf_inboard_peak  |  radius of maximum toroidal field (m)
      !! b_tf_inboard_peak |  peak field at toroidal field coil (T)

      !! This constraint is depreciated

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Constraints output

      ! This constraint is depreciated
      call report_error(236)

      tmp_con = 1.0D0
      tmp_err = 0.0D0
      tmp_symbol = '='
      tmp_units = ''

   end subroutine constraint_eqn_010

   subroutine constraint_eqn_011(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radial build
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for radial build
      !! (This is a consistency equation.)
      !! #=# build
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rbld : input real : sum of thicknesses to the major radius (m)
      !! rmajor : input real : plasma major radius (m)
      use build_variables, only: rbld
      use physics_variables, only: rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - rbld/rmajor
      tmp_con = rmajor * (1.0D0 - tmp_cc)
      tmp_err = rmajor * tmp_cc
      tmp_symbol = '='
      tmp_units = 'm'

   end subroutine constraint_eqn_011

   subroutine constraint_eqn_014(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation to fix number of NBI decay lengths to plasma centre
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation to fix number of NBI decay lengths to plasma centre
      !! This is a consistency equation
      !! #=# current_drive
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! n_beam_decay_lengths_core : input real : neutral beam e-decay lengths to plasma centre
      !! tbeamin : input real : permitted neutral beam e-decay lengths to plasma centre
      use current_drive_variables, only: n_beam_decay_lengths_core, tbeamin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - n_beam_decay_lengths_core/tbeamin
      tmp_con = tbeamin * (1.0D0 - tmp_cc)
      tmp_err = tbeamin * tmp_cc
      tmp_symbol = '='
      tmp_units = ''

   end subroutine constraint_eqn_014

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

   subroutine constraint_eqn_017(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radiation power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for radiation power upper limit
      !! #=# physics
      !! #=#=# fradpwr, pradmaxpv
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! f_alpha_plasma : input real : fraction of alpha power deposited in plasma
      !! p_hcd_injected_total_mw : input real : total auxiliary injected power (MW)
      !! vol_plasma : input real : plasma volume (m3)
      !! alpha_power_density_total : input real : alpha power per volume (MW/m3)
      !! charged_power_density :  input real : non-alpha charged particle fusion power per volume (MW/m3)
      !! pden_plasma_ohmic_mw : input real : ohmic heating power per volume (MW/m3)
      !! fradpwr : input real : f-value for core radiation power limit
      !! pden_plasma_rad_mw : input real : total radiation power per volume (MW/m3)
      use physics_variables, only: f_alpha_plasma, vol_plasma, alpha_power_density_total, charged_power_density, pden_plasma_ohmic_mw, pden_plasma_rad_mw
      use current_drive_variables, only: p_hcd_injected_total_mw
      use constraint_variables, only: fradpwr
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      real(dp) :: pradmaxpv
      !! Maximum possible power/vol_plasma that can be radiated (local)

      pradmaxpv = p_hcd_injected_total_mw/vol_plasma + alpha_power_density_total*f_alpha_plasma + charged_power_density + pden_plasma_ohmic_mw
      tmp_cc =  pden_plasma_rad_mw/pradmaxpv - 1.0D0 * fradpwr
      tmp_con = pradmaxpv * (1.0D0 - tmp_cc)
      tmp_err = pden_plasma_rad_mw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m3'

   end subroutine constraint_eqn_017

   subroutine constraint_eqn_018(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for divertor heat load upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for divertor heat load upper limit
      !! #=# divertor
      !! #=#=# fpflux_div_heat_load_mw, pflux_div_heat_load_max_mw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpflux_div_heat_load_mw : input real : f-value for divertor heat load
      !! pflux_div_heat_load_max_mw : input real : heat load limit (MW/m2)
      !! pflux_div_heat_load_mw : input real : divertor heat load (MW/m2)
      use constraint_variables, only: fpflux_div_heat_load_mw
      use divertor_variables, only: pflux_div_heat_load_max_mw, pflux_div_heat_load_mw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  pflux_div_heat_load_mw/pflux_div_heat_load_max_mw - 1.0D0 * fpflux_div_heat_load_mw
      tmp_con = pflux_div_heat_load_max_mw * (1.0D0 - tmp_cc)
      tmp_err = pflux_div_heat_load_mw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m2'

   end subroutine constraint_eqn_018

   subroutine constraint_eqn_019(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for MVA (power) upper limit: resistive TF coil set
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for MVA upper limit
      !! #=# tfcoil
      !! #=#=# fmva, mvalim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! tfcpmw : input real : peak resistive TF coil inboard leg power (total) (MW)
      !! tflegmw : input real : TF coil outboard leg resistive power (total) (MW)
      !! fmva : input real : f-value for maximum MVA
      !! mvalim : input real : MVA limit for resistive TF coil set (total) (MW)
      use tfcoil_variables, only: tfcpmw, tflegmw
      use constraint_variables, only: fmva, mvalim
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      ! totmva : local real : total MVA in TF coil (MW)
      real(dp) :: totmva

      totmva = tfcpmw + tflegmw
      tmp_cc =  totmva/mvalim - 1.0D0 * fmva
      tmp_con = mvalim * (1.0D0 - tmp_cc)
      tmp_err = totmva * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MVA'

   end subroutine constraint_eqn_019

   subroutine constraint_eqn_020(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutral beam tangency radius upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutral beam tangency radius upper limit
      !! #=# current_drive
      !! #=#=# fportsz, rtanmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fportsz : input real : f-value for neutral beam tangency radius limit
      !! rtanmax : input real : maximum tangency radius for centreline of beam (m)
      !! rtanbeam : input real : neutral beam centreline tangency radius (m)
      use constraint_variables, only: fportsz
      use current_drive_variables, only: rtanmax, rtanbeam
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  rtanbeam/rtanmax - 1.0D0 * fportsz
      tmp_con = rtanmax * (1.0D0 - tmp_cc)
      tmp_err = rtanbeam * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_020

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

   subroutine constraint_eqn_022(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for divertor collision/connection length ratio upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for divertor collision/connection length ratio upper limit
      !! #=# divertor
      !! #=#=# fdivcol, rlenmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fdivcol : input real : f-value for divertor collisionality
      !! rlenmax : input real : maximum value for length ratio (rlclolcn)
      !! rlclolcn : input real : ratio of collision length / connection length

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Constraints output

      ! This constraint is depreciated
      call report_error(289)

      tmp_con = 1.0D0
      tmp_err = 0.0D0
      tmp_symbol = '='
      tmp_units = ''

   end subroutine constraint_eqn_022

   subroutine constraint_eqn_023(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for conducting shell radius / rminor upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for conducting shell radius / rminor upper limit
      !! #=# physics
      !! #=#=#fr_conducting_wall, f_r_conducting_wall
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rminor : input real : plasma minor radius (m)
      !! dr_fw_plasma_gap_outboard : input real : gap between plasma and first wall, outboard side (m)
      !! dr_fw_outboard : input real : outboard first wall thickness, initial estimate (m)
      !! dr_blkt_outboard : input real : outboard blanket thickness (m)
      !!fr_conducting_wall : input real : f-value for conducting wall radius / rminor limit
      !! f_r_conducting_wall : input real : maximum ratio of conducting wall distance to plasma minor radius for vertical stability
      use physics_variables, only: rminor, f_r_conducting_wall
      use build_variables, only: dr_fw_plasma_gap_outboard, dr_fw_outboard, dr_blkt_outboard
      use constraint_variables, only:fr_conducting_wall
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      ! rcw : local real : conducting shell radius (m)
      real(dp) :: rcw

      rcw = rminor + dr_fw_plasma_gap_outboard + dr_fw_outboard + dr_blkt_outboard
      tmp_cc =  rcw / (f_r_conducting_wall * rminor) - 1.0D0 * fr_conducting_wall
      tmp_con = f_r_conducting_wall*rminor * (1.0D0 - tmp_cc)
      tmp_err = rcw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_023

   subroutine constraint_eqn_024(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for beta upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for beta upper limit
      !! #=# physics
      !! #=#=# fbeta_max, beta_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! i_beta_component : input integer : switch for beta limit scaling (constraint equation  24):<UL>
      !! <LI> = 0 apply limit to total beta;
      !! <LI> = 1 apply limit to thermal beta;
      !! <LI> = 2 apply limit to thermal + neutral beam beta
      !! <LI> = 3 apply limit to toroidal beta </UL>
      !! istell : input integer : switch for stellarator option (set via <CODE>device.dat</CODE>):<UL>
      !! <LI> = 0 use tokamak model;
      !! <LI> = 1 use stellarator model</UL>
      !! fbeta_max : input real : f-value for beta limit
      !! beta_max : input real : allowable beta
      !! beta : input real : total plasma beta (calculated if ipedestal =3)
      !! beta_fast_alpha : input real : fast alpha beta component
      !! beta_beam : input real : neutral beam beta component
      !! bt : input real : toroidal field
      !! btot : input real : total field
      use physics_variables, only: i_beta_component, beta_max, beta, beta_beam, beta_fast_alpha, bt, btot
      use stellarator_variables, only: istell
      use constraint_variables, only: fbeta_max
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! Include all beta components: relevant for both tokamaks and stellarators
      if ((i_beta_component == 0).or.(istell /= 0)) then
         tmp_cc =  beta/beta_max - 1.0D0 * fbeta_max
         tmp_con = beta_max
         tmp_err = beta_max - beta / fbeta_max
         tmp_symbol = '<'
         tmp_units = ''
      ! Here, the beta limit applies to only the thermal component, not the fast alpha or neutral beam parts
      else if (i_beta_component == 1) then
         tmp_cc = (beta-beta_fast_alpha-beta_beam)/beta_max - 1.0D0 * fbeta_max
         tmp_con = beta_max
         tmp_err = beta_max - (beta-beta_fast_alpha-beta_beam) / fbeta_max
         tmp_symbol = '<'
         tmp_units = ''
      ! Beta limit applies to thermal + neutral beam: components of the total beta, i.e. excludes alphas
      else if (i_beta_component == 2) then
         tmp_cc = (beta-beta_fast_alpha)/beta_max - 1.0D0 * fbeta_max
         tmp_con = beta_max * (1.0D0 - tmp_cc)
         tmp_err = (beta-beta_fast_alpha) * tmp_cc
         tmp_symbol = '<'
         tmp_units = ''
      ! Beta limit applies to toroidal beta
      else if (i_beta_component == 3) then
         tmp_cc =  (beta*(btot/bt)**2)/beta_max - 1.0D0 * fbeta_max
         tmp_con = beta_max
         tmp_err = beta_max - (beta*(btot/bt)**2) / fbeta_max
         tmp_symbol = '<'
         tmp_units = ''
      end if

   end subroutine constraint_eqn_024

   subroutine constraint_eqn_025(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for peak toroidal field upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for peak toroidal field upper limit
      !! #=# tfcoil
      !! #=#=# fpeakb, bmxlim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpeakb : input real : f-value for maximum toroidal field
      !! bmxlim : input real : maximum peak toroidal field (T)
      !! b_tf_inboard_peak : input real : mean peak field at TF coil (T)
      use constraint_variables, only: fpeakb, bmxlim
      use tfcoil_variables, only: b_tf_inboard_peak
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  b_tf_inboard_peak/bmxlim - 1.0D0 * fpeakb
      tmp_con = bmxlim * (1.0D0 - tmp_cc)
      tmp_err = b_tf_inboard_peak * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'T'

   end subroutine constraint_eqn_025

   subroutine constraint_eqn_026(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid current density upper limit at EOF
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid current density upper limit at EOF
      !! #=# pfcoil
      !! #=#=# fjohc, j_cs_critical_flat_top_end
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjohc : input real : f-value for central solenoid current at end-of-flattop
      !! j_cs_critical_flat_top_end : input real : allowable central solenoid current density at end of flat-top (A/m2)
      !! j_cs_flat_top_end : input real : central solenoid overall current density at end of flat-top (A/m2)
      use constraint_variables, only: fjohc
      use pfcoil_variables, only: j_cs_critical_flat_top_end, j_cs_flat_top_end
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  j_cs_flat_top_end/j_cs_critical_flat_top_end - 1.0D0 * fjohc
      tmp_con = j_cs_critical_flat_top_end
      tmp_err = j_cs_critical_flat_top_end - j_cs_flat_top_end / fjohc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_026

   subroutine constraint_eqn_027(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid current density upper limit at BOP
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid current density upper limit at BOP
      !! #=# pfcoil
      !! #=#=# fjohc0, j_cs_critical_pulse_start
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjohc0 : input real : f-value for central solenoid current at beginning of pulse
      !! j_cs_critical_pulse_start : input real : allowable central solenoid current density at beginning of pulse (A/m2)
      !! j_cs_pulse_start : input real : central solenoid overall current density at beginning of pulse (A/m2)
      use constraint_variables, only: fjohc0
      use pfcoil_variables, only: j_cs_critical_pulse_start, j_cs_pulse_start
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  j_cs_pulse_start/j_cs_critical_pulse_start - 1.0D0 * fjohc0
      tmp_con = j_cs_critical_pulse_start
      tmp_err = j_cs_critical_pulse_start - j_cs_pulse_start / fjohc0
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_027

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

   subroutine constraint_eqn_029(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for inboard major radius: This is a consistency equation
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for inboard major radius: This is a consistency equation
      !! #=# build
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rmajor : input real : plasma major radius (m) (iteration variable 3)
      !! rminor : input real : plasma minor radius (m)
      !! rinboard : input real : plasma inboard radius (m)
      use physics_variables, only: rmajor, rminor
      use build_variables, only: rinboard
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - (rmajor - rminor) / rinboard
      tmp_con = rinboard * (1.0D0 - tmp_cc)
      tmp_err = rinboard * tmp_cc
      tmp_symbol = '='
      tmp_units = 'm'

   end subroutine constraint_eqn_029

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

   subroutine constraint_eqn_043(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for average centrepost temperature: This is a consistency equation (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for average centrepost temperature: This is a consistency equation (TART)
      !! #=# tfcoil
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! temp_cp_average : input real : average temp of TF coil inboard leg conductor (C)e
      !! tcpav2 : input real : centrepost average temperature (C) (for consistency)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use tfcoil_variables, only: temp_cp_average, tcpav2
      use physics_variables, only: itart
      use tfcoil_variables, only:  i_tf_sup

      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(7)

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         temp_cp_average = temp_cp_average - 273.15D0
         tcpav2 = tcpav2 - 273.15D0
      end if

      tmp_cc =   1.0D0 - temp_cp_average/tcpav2
      tmp_con = tcpav2 * (1.0D0 - tmp_cc)
      tmp_err = tcpav2 * tmp_cc
      tmp_symbol = '='
      tmp_units = 'deg C'

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         temp_cp_average = temp_cp_average + 273.15D0
         tcpav2 = tcpav2 + 273.15D0
      end if



   end subroutine constraint_eqn_043

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

   subroutine constraint_eqn_051(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation to enforce startup flux = available startup flux
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation to enforce startup flux = available startup flux
      !! #=# pfcoil
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! vs_plasma_res_ramp : input real : resistive losses in startup V-s (Wb)
      !! vs_plasma_ind_ramp : input real :  internal and external plasma inductance V-s (Wb))
      !! vs_cs_pf_total_ramp : input real :  total flux swing for startup (Wb)
      use physics_variables, only: vs_plasma_res_ramp, vs_plasma_ind_ramp
      use pfcoil_variables, only: vs_cs_pf_total_ramp, fvs_cs_pf_total_ramp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fvs_cs_pf_total_ramp * abs((vs_plasma_res_ramp+vs_plasma_ind_ramp) / vs_cs_pf_total_ramp)
      tmp_con = vs_cs_pf_total_ramp * (1.0D0 - tmp_cc)
      tmp_err = vs_cs_pf_total_ramp * tmp_cc
      tmp_symbol = '='
      tmp_units = 'V.s'

   end subroutine constraint_eqn_051

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

   subroutine constraint_eqn_085(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equality constraint for the centerpost (CP) lifetime
      !! Author : S Kahn
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Depending on the chosen option : i_cp_lifetime
      !!  - 0 : The CP full power year lifelime is set by the user (cplife_input)
      !!  - 1 : The CP lifelime is equal to the divertor one
      !!  - 2 : The CP lifetime is equal to the breeding blankets one
      !!  - 3 : The CP lifetime is equal to the plant one
      !! #=# availability
      !! #=#=# consistency
      !! Logic change during pre-factoring: err, symbol, units will be assigned
      !! only if present.
      !! cplife : input real : calculated CP full power year lifetime (years)
      !! life_blkt_fpy : input real : calculated first wall/blanket power year lifetime (years)
      !! divlife : input real : calculated divertor  power year lifetime (years)
      !! i_cp_lifetime : input integer : switch chosing which plant element the CP
      !!                                 the CP lifetime must equate
      use cost_variables, only : cplife, divlife, cplife_input, &
         tlife, i_cp_lifetime
      use fwbs_variables, only : life_blkt_fpy

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Constraints output


      ! The CP lifetime is equal to the the divertor one
      if  ( i_cp_lifetime == 0 ) then
         tmp_cc = 1.0D0 - cplife/cplife_input

      else if ( i_cp_lifetime == 1 ) then
         tmp_cc = 1.0D0 - cplife/divlife

      ! The CP lifetime is equal to the tritium breeding blankets / FW one
      else if ( i_cp_lifetime == 2 ) then
         tmp_cc = 1.0D0 - cplife/life_blkt_fpy

      ! The CP lifetime is equal to the
      else if ( i_cp_lifetime == 3 ) then
         tmp_cc = 1.0D0 - cplife/tlife
      end if

      tmp_con = divlife * (1.0D0 - tmp_cc)
      tmp_err = divlife * tmp_cc
      tmp_symbol = '='
      tmp_units = 'years'

   end subroutine constraint_eqn_085

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
      tmp_symbol = '<'
      tmp_units = 'MW'
   end subroutine constraint_eqn_091

   subroutine constraint_eqn_092(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for checking is D/T ratio is consistent, and sums to 1.
      !! author: G Turkington, UKAEA
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! f_deuterium : input : fraction of deuterium ions
      !! f_tritium  : input : fraction of tritium ions
      !! f_helium3  : input : fraction of helium-3 ions
      use physics_variables, only: f_deuterium, f_tritium, f_helium3
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units


      ! Iterate over f_tritium and calculate f_deuterium
      f_deuterium = 1.0D0 - (f_tritium + f_helium3)
      tmp_cc = 1.0D0 - (f_deuterium + f_tritium + f_helium3)
      tmp_con = 1.0D0
      tmp_err = tmp_con * tmp_cc
      tmp_symbol = '='
      tmp_units = 'fraction'

   end subroutine constraint_eqn_092


end module constraints
