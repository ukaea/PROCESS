! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scan_module

  !! Module containing routines to perform a parameter scan
  !! author: P J Knight, CCFE, Culham Science Centre
  !! None
  !! This module contains routines to perform a parameter scan
  !! over a range of values of a particular scanning variable.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public

  integer, parameter :: ipnscns = 1000
  !! Maximum number of scan points

  integer, parameter :: ipnscnv = 79
  !! Number of available scan variables

  integer, parameter :: noutvars = 84
  integer, parameter :: width = 110

  integer :: scan_dim
  !! 1-D or 2-D scan switch (1=1D, 2=2D)

  integer :: isweep
  !! Number of scan points to calculate

  integer :: isweep_2
  !! Number of 2D scan points to calculate

  integer :: nsweep
  !! Switch denoting quantity to scan:<UL>
  !!         <LI> 1  aspect
  !!         <LI> 2  hldivlim
  !!         <LI> 3  pnetelin
  !!         <LI> 4  hfact
  !!         <LI> 5  oacdcp
  !!         <LI> 6  walalw
  !!         <LI> 7  beamfus0
  !!         <LI> 8  fqval
  !!         <LI> 9  te
  !!         <LI> 10 boundu(15: fvs)
  !!         <LI> 11 dnbeta
  !!         <LI> 12 bscfmax (use negative values only)
  !!         <LI> 13 boundu(10: hfact)
  !!         <LI> 14 fiooic
  !!         <LI> 15 fjprot
  !!         <LI> 16 rmajor
  !!         <LI> 17 bmxlim
  !!         <LI> 18 gammax
  !!         <LI> 19 boundl(16: ohcth)
  !!         <LI> 20 tbrnmn
  !!         <LI> 21 not used
  !!         <LI> 22 cfactr (N.B. requires iavail=0)
  !!         <LI> 23 boundu(72: fipir)
  !!         <LI> 24 powfmax
  !!         <LI> 25 kappa
  !!         <LI> 26 triang
  !!         <LI> 27 tbrmin (for blktmodel > 0 only)
  !!         <LI> 28 bt
  !!         <LI> 29 coreradius
  !!         <LI> 30 fimpvar # OBSOLETE
  !!         <LI> 31 taulimit
  !!         <LI> 32 epsvmc
  !!         <LI> 33 ttarget
  !!         <LI> 34 qtargettotal
  !!         <LI> 35 lambda_q_omp
  !!         <LI> 36 lambda_target
  !!         <LI> 37 lcon_factor
  !!         <LI> 38 Neon upper limit
  !!         <LI> 39 Argon upper limit
  !!         <LI> 40 Xenon upper limit
  !!         <LI> 41 blnkoth
  !!         <LI> 42 Argon fraction fimp(9)
  !!         <LI> 43 normalised minor radius at which electron cyclotron current drive is maximum
  !!         <LI> 44 Allowable maximum shear stress (Tresca) in tf coil structural material
  !!         <LI> 45 Minimum allowable temperature margin ; tf coils
  !!         <LI> 46 boundu(150) fgwsep
  !!         <LI> 47 impurity_enrichment(9) Argon impurity enrichment
  !!         <LI> 48 TF coil - n_pancake (integer turn winding pack)
  !!         <LI> 49 TF coil - n_layer (integer turn winding pack)
  !!         <LI> 50 Xenon fraction fimp(13)
  !!         <LI> 51 Power fraction to lower DN Divertor ftar
  !!         <LI> 52 SoL radiation fraction
  !!         <LI> 54 GL_nbti upper critical field at 0 Kelvin
  !!         <LI> 55 `shldith` : Inboard neutron shield thickness
  !!         <LI> 56 crypmw_max: Maximum cryogenic power (ixx=164, ixc=87)
  !!         <LI> 57 `bt` lower boundary
  !!         <LI> 58 `scrapli` : Inboard plasma-first wall gap
  !!         <LI> 59 `scraplo` : Outboard plasma-first wall gap
  !!         <LI> 60 sig_tf_wp_max: Allowable stress in TF Coil conduit (Tresca)
  !!         <LI> 61 copperaoh_m2_max : CS coil current / copper area
  !!         <LI> 62 coheof : CS coil current density at EOF
  !!         <LI> 63 ohcth : CS thickness (m)
  !!         <LI> 64 ohhghf : CS height (m)
  !!         <LI> 65 n_cycle_min : Minimum cycles for CS stress model constraint 90
  !!         <LI> 66 oh_steel_frac: Steel fraction in CS coil
  !!         <LI> 67 t_crack_vertical: Initial crack vertical dimension (m) </UL>
  !!         <LI> 68 `inlet_temp_liq' : Inlet temperature of blanket liquid metal coolant/breeder (K)
  !!         <LI> 69 `outlet_temp_liq' : Outlet temperature of blanket liquid metal coolant/breeder (K)
  !!         <LI> 70 `blpressure_liq' : Blanket liquid metal breeder/coolant pressure (Pa)
  !!         <LI> 71 `n_liq_recirc' : Selected number of liquid metal breeder recirculations per day
  !!         <LI> 72 `bz_channel_conduct_liq' : Conductance of liquid metal breeder duct walls (A V-1 m-1)
  !!         <LI> 73 `pnuc_fw_ratio_dcll' : Ratio of FW nuclear power as fraction of total (FW+BB)
  !!         <LI> 74 `f_nuc_pow_bz_struct' : Fraction of BZ power cooled by primary coolant for dual-coolant balnket
  !!         <LI> 75 pitch : pitch of first wall cooling channels (m)
  !!         <LI> 76 etath : Thermal conversion eff.
  !!         <LI> 77 startupratio : Gyrotron redundancy
  !!         <LI> 78 fkind : Multiplier for Nth of a kind costs
  !!         <LI> 79 etaech : ECH wall plug to injector efficiency

  integer :: nsweep_2
  !! nsweep_2 /3/ : switch denoting quantity to scan for 2D scan:

  real(dp), dimension(ipnscns) :: sweep
  !! sweep(ipnscns) /../: actual values to use in scan

  real(dp), dimension(ipnscns) :: sweep_2
  !! sweep_2(ipnscns) /../: actual values to use in 2D scan

  ! Vars in subroutines scan_1d and scan_2d requiring re-initialising before
  ! each new run
  logical :: first_call_1d
  logical :: first_call_2d

contains

  subroutine init_scan_module
    !! Initialise module variables
    implicit none

    scan_dim = 1
    isweep = 0
    isweep_2 = 0
    nsweep = 1
    nsweep_2 = 3
    sweep = 0.0D0
    sweep_2 = 0.0D0
    first_call_1d = .true.
    first_call_2d = .true.
  end subroutine init_scan_module

  subroutine scan_1d_write_point_header(iscan)
    use global_variables, only: iscan_global, xlabel, vlabel
    use constants, only: mfile, nout
    use process_output, only: ovarin, ostars, oblnkl
    implicit none
    integer, intent(in) :: iscan
    !! Scan point number

    ! Makes iscan available globally (read-only)
    iscan_global = iscan

    call scan_select(nsweep, sweep, iscan, vlabel, xlabel)

    ! Write banner to output file
    call oblnkl(nout)
    call ostars(nout,width)
    write(nout,10) ' ***** Scan point ', iscan,' of ',isweep,': &
        ',trim(xlabel),', ',trim(vlabel),' = ',sweep(iscan),' *****'
10     format(a,i2,a,i2,5a,1pe10.3,a)
    call ostars(nout,width)

    ! Write additional information to mfile
    call oblnkl(mfile)
    call ovarin(mfile,'Scan point number','(iscan)',iscan)

    ! Call the optimization routine VMCON at this scan point
    write(*,20)'Starting scan point ',iscan, ': ', trim(xlabel),', &
        ',trim(vlabel),' = ',sweep(iscan)
20     format(a,i2,a,4a,1pe10.3)
  end subroutine scan_1d_write_point_header

  subroutine scan_1d_store_output(iscan, ifail, noutvars_, ipnscns_, outvar)
    use constraint_variables, only: taulimit
    use cost_variables, only: cdirt, coe, coeoam, coefuelt, c222, ireactor, &
      capcost, coecap, c221
    use current_drive_variables, only: pheat, pinjmw, bootipf, enbeam, bigq
    use divertor_variables, only: hldiv
    use error_handling, only: errors_on
    use heat_transport_variables, only: pgrossmw, pinjwp, pnetelmw
    use impurity_radiation_module, only: fimp
    use pfcoil_variables, only: whtpf
    use pf_power_variables, only: srcktpm
    use process_output, only: oblnkl
    use numerics, only: sqsumsq
    use tfcoil_variables, only: tfareain, wwp2, sig_tf_wp, tfcmw, tcpmax, oacdcp, &
      tfcpmw, fcutfsu, acond, fcoolcp, rcool, whttf, ppump, vcool, wwp1, n_tf, &
      dr_tf_wp, b_crit_upper_nbti
    use fwbs_variables, only: tpeak
    use physics_variables, only: q, aspect, pradmw, dene, powfmw, btot, tesep, &
      pdivt, ralpne, ten, betap, hfac, teped, palpnb, qlim, rmajor, wallmw, &
      beta, betalim, bt, plascur
    use global_variables, only: verbose, maxcal, runtitle, run_tests
    use constants, only: nout
    implicit none

    integer, intent(in) :: iscan
    integer, intent(in) :: ifail
    ! outvar
    integer, intent(in) :: noutvars_, ipnscns_
    real(dp), dimension(noutvars_,ipnscns_), intent(out) :: outvar

    ! Turn off error reporting (until next output)
    errors_on = .false.

    ! Store values for MFILE.DAT output
    outvar( 1,iscan) = dble(ifail)
    outvar( 2,iscan) = sqsumsq
    outvar( 3,iscan) = coe
    outvar( 4,iscan) = coecap
    outvar( 5,iscan) = coefuelt
    outvar( 6,iscan) = coeoam
    outvar( 7,iscan) = capcost
    outvar( 8,iscan) = c221 + c222
    outvar( 9,iscan) = cdirt / 1.0D3
    outvar(10,iscan) = rmajor
    outvar(11,iscan) = aspect
    outvar(12,iscan) = 1.0D-6 * plascur
    outvar(13,iscan) = bt
    outvar(14,iscan) = btot
    outvar(15,iscan) = q
    outvar(16,iscan) = qlim
    outvar(17,iscan) = beta
    outvar(18,iscan) = betalim
    outvar(19,iscan) = betap / aspect
    outvar(20,iscan) = ten/10.0D0
    outvar(21,iscan) = dene/1.0D20
    outvar(22,iscan) = hfac(6)
    outvar(23,iscan) = hfac(7)
    outvar(24,iscan) = powfmw
    outvar(25,iscan) = palpnb * 5.0D0
    outvar(26,iscan) = wallmw
    outvar(27,iscan) = pinjmw
    outvar(28,iscan) = pinjwp
    outvar(29,iscan) = pheat
    outvar(30,iscan) = pinjmw - pheat
    outvar(31,iscan) = bigq
    outvar(32,iscan) = bootipf
    outvar(33,iscan) = enbeam/1.0D3
    outvar(34,iscan) = hldiv
    outvar(35,iscan) = tfcmw
    outvar(36,iscan) = whttf
    outvar(37,iscan) = sig_tf_wp
    outvar(38,iscan) = oacdcp/1.0D6
    outvar(39,iscan) = tcpmax
    outvar(40,iscan) = tfcpmw
    outvar(41,iscan) = fcoolcp
    outvar(42,iscan) = rcool
    outvar(43,iscan) = vcool
    outvar(44,iscan) = ppump/1.0D6
    outvar(45,iscan) = 1.0D-3 * srcktpm
    outvar(46,iscan) = whtpf
    outvar(47,iscan) = pgrossmw
    outvar(48,iscan) = pnetelmw
    if (ireactor == 1) then
        outvar(49,iscan) = (pgrossmw-pnetelmw) / pgrossmw
    else
        outvar(49,iscan) = 0.0D0
    end if
    outvar(50,iscan) = pdivt/rmajor
    !outvar(51,iscan) = fimpvar #OBSOLETE
    outvar(51,iscan) = 0.0d0
    outvar(52,iscan) = pradmw
    outvar(53,iscan) = tpeak
    outvar(54,iscan) = fcutfsu
    outvar(55,iscan) = (wwp1+wwp2)*dr_tf_wp
    outvar(56,iscan) = acond
    outvar(57,iscan) = tfareain/n_tf
    outvar(58,iscan) = taulimit
    outvar(66,iscan) = ralpne
    outvar(69,iscan) = fimp(1)
    outvar(70,iscan) = fimp(2)
    outvar(71,iscan) = fimp(3)
    outvar(72,iscan) = fimp(4)
    outvar(73,iscan) = fimp(5)
    outvar(74,iscan) = fimp(6)
    outvar(75,iscan) = fimp(7)
    outvar(76,iscan) = fimp(8)
    outvar(77,iscan) = fimp(9)
    outvar(78,iscan) = fimp(10)
    outvar(79,iscan) = fimp(11)
    outvar(80,iscan) = fimp(12)
    outvar(81,iscan) = fimp(13)
    outvar(82,iscan) = fimp(14)
    outvar(83,iscan) = teped
  end subroutine scan_1d_store_output

  subroutine scan_1d_write_plot(iscan, outvar)
    use global_variables, only: icase, xlabel
    use constants, only: nplot, mfile
    use process_output, only: ovarin
    implicit none

    integer, intent(inout) :: iscan
    real(dp), dimension(:,:), intent(in) :: outvar

    character(len=48) :: tlabel
    integer :: ivar
    character(len=25), dimension(noutvars), save :: plabel

    tlabel = icase

    !  Set up labels for plotting output
    !  Use underscores instead of spaces
    if (first_call_1d) then
        plabel( 1) = 'Ifail____________________'
        plabel( 2) = 'Sqsumsq__________________'
        plabel( 3) = 'Electric_cost_(mil/kwh)__'
        plabel( 4) = 'Capital_cost_(mil/kwh)___'
        plabel( 5) = 'Fuel_cost_(mil/kwh)______'
        plabel( 6) = 'Operations_cost_(mil/kwh)'
        plabel( 7) = 'Capital_cost_(millions)__'
        plabel( 8) = 'Core_costs_(millions)____'
        plabel( 9) = 'Direct_cost_(billions)___'
        plabel(10) = 'Major_Radius_(m)_________'
        plabel(11) = 'Aspect_Ratio_____________'
        plabel(12) = 'Plasma_Current_(MA)______'
        plabel(13) = 'B_Toroidal_Axis_(T)______'
        plabel(14) = 'B_total_on_axis_(T)______'
        plabel(15) = 'Safety_Factor____________'
        plabel(16) = 'qlim_(zero_if_ishape=0)__'
        plabel(17) = 'Beta_____________________'
        plabel(18) = 'Beta_Limit_______________'
        plabel(19) = 'Epsilon_Beta_Poloidal____'
        plabel(20) = 'Dens.weight_Te_(10keV)___'
        plabel(21) = 'Average_Dens_(10^20/m^3)_'
        plabel(22) = 'H-fact_Iter_Power________'
        plabel(23) = 'H-fact_Iter_Offset_______'
        plabel(24) = 'Fusion_Power_(MW)________'
        plabel(25) = 'nb_Fusion_Power_(MW)_____'
        plabel(26) = 'Wall_Load_(MW/m^2)_______'
        plabel(27) = 'Injection_Power_(MW)_____'
        plabel(28) = 'Inject_Pwr_Wall_Plug_(MW)'
        plabel(29) = 'Heating_Power_(MW)_______'
        plabel(30) = 'Current_Drive_(MW)_______'
        plabel(31) = 'Big_Q____________________'
        plabel(32) = 'Bootstrap_Fraction_______'
        plabel(33) = 'Neutral_Beam_Energy_(MeV)'
        plabel(34) = 'Divertor_Heat_(MW/m^2)___'
        plabel(35) = 'TF_coil_Power_(MW)_______'
        plabel(36) = 'TF_coil_weight_(kg)______'
        plabel(37) = 'vM_stress_in_TF_case_(Pa)'
        plabel(38) = 'J_TF_inboard_leg_(MA/m^2)'
        plabel(39) = 'Centrepost_max_T_(TART)__'
        plabel(40) = 'Res_TF_inbrd_leg_Pwr_(MW)'
        plabel(41) = 'Coolant_Fraction_Ctr.____'
        plabel(42) = 'C/P_coolant_radius_(m)___'
        plabel(43) = 'C/P_coolant_velocity(m/s)'
        plabel(44) = 'C/P_pump_power_(MW)______'
        plabel(45) = 'PF_coil_Power_(MW)_______'
        plabel(46) = 'PF_coil_weight_(kg)______'
        plabel(47) = 'Gross_Elect_Pwr_(MW)_____'
        plabel(48) = 'Net_electric_Pwr_(MW)____'
        plabel(49) = 'Recirculating_Fraction___'
        plabel(50) = 'Psep/R___________________'
        plabel(51) = '' !OBSOLETE
        plabel(52) = 'Tot._radiation_power_(MW)'
        plabel(53) = 'First_wall_peak_temp_(K)_'
        plabel(54) = 'Cu_frac_TFC_conductor____'
        plabel(55) = 'Winding_pack_area_TFC(m2)'
        plabel(56) = 'Conductor_area_TFC_(m2)__'
        plabel(57) = 'Area_TF_inboard_leg_(m2)_'
        plabel(58) = 'Taup/taueff_lower_limit__'
        plabel(59) = 'Plasma_temp_at_sep_[keV]_'
        plabel(60) = 'SOL_density_at_OMP_______'
        plabel(61) = 'Power_through__separatrix'
        plabel(62) = 'neomp/nesep______________'
        plabel(63) = 'qtargettotal_____________'
        plabel(64) = 'Total_pressure_at_target_'
        plabel(65) = 'Temperature_at_target____'
        plabel(66) = 'Helium_fraction__________'
        plabel(67) = 'Momentum_loss_factor_____'
        plabel(68) = 'totalpowerlost_[W]_______'
        plabel(69) = 'H__concentration_________'
        plabel(70) = 'He_concentration_________'
        plabel(71) = 'Be_concentration_________'
        plabel(72) = 'C__concentration_________'
        plabel(73) = 'N__concentration_________'
        plabel(74) = 'O__concentration_________'
        plabel(75) = 'Ne_concentration_________'
        plabel(76) = 'Si_concentration_________'
        plabel(77) = 'Ar_concentration_________'
        plabel(78) = 'Fe_concentration_________'
        plabel(79) = 'Ni_concentration_________'
        plabel(80) = 'Kr_concentration_________'
        plabel(81) = 'Xe_concentration_________'
        plabel(82) = 'W__concentration_________'
        plabel(83) = 'teped____________________'
        plabel(84) = 'Max_field_on_TF_coil_____'
        call ovarin(mfile,'Number of scan points','(isweep)',isweep)
        call ovarin(mfile,'Scanning variable number','(nsweep)',nsweep)

        first_call_1d = .false.
     end if

  end subroutine scan_1d_write_plot

  subroutine scan_2d_init
    !! Routine to call 2-D scan
    !! author: J Morris, UKAEA, Culham Science Centre
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: mfile
    use process_output, only: ovarin
    implicit none

    !  Set up labels for plotting output
    !  Use underscores instead of spaces
    call ovarin(mfile,'Number of first variable scan points','(isweep)',isweep)
    call ovarin(mfile,'Number of second variable scan points','(isweep_2)',isweep_2)
    call ovarin(mfile,'Scanning first variable number','(nsweep)',nsweep)
    call ovarin(mfile,'Scanning second variable number','(nsweep_2)',nsweep_2)
    call ovarin(mfile,'Scanning second variable number','(nsweep_2)',nsweep_2)
    call ovarin(mfile,'Scanning second variable number','(nsweep_2)',nsweep_2)
  end subroutine scan_2d_init

  subroutine scan_2d_write_point_header(iscan, iscan_1, iscan_2, iscan_R)
    use process_output, only: oblnkl, ostars, ovarin
    use global_variables, only: vlabel, vlabel_2, xlabel, xlabel_2, iscan_global
    use constants, only: nout, mfile
    implicit none

    integer, intent(in) :: iscan
    integer, intent(in) :: iscan_1
    integer, intent(in) :: iscan_2
    integer, intent(out) :: iscan_R

    integer :: ifail

    if (mod(iscan_1,2)==0) then
        iscan_R = isweep_2 - iscan_2 + 1
    else
        iscan_R = iscan_2
    end if
    ! Makes iscan available globally (read-only)
    iscan_global = iscan

    call scan_select(nsweep, sweep, iscan_1, vlabel, xlabel)
    call scan_select(nsweep_2, sweep_2, iscan_R, vlabel_2, xlabel_2)

    ! Write banner to output file
    call oblnkl(nout)
    call ostars(nout,width)
    write(nout,10) iscan, isweep*isweep_2, trim(vlabel), &
        sweep(iscan_1), trim(vlabel_2), sweep_2(iscan_R)
! 10    format(a, i2, a, i2, 5a, 1pe10.3, a)
10  format(' ***** 2D scan point ', i3, ' of ', i3, ' : ', a, ' = ', &
            1pe10.3, ' and ', a, ' = ', 1pe10.3, ' *****')
    call ostars(nout,width)

    ! Write additional information to mfile
    call oblnkl(mfile)
    call ovarin(mfile,'Scan point number','(iscan)',iscan)

    ! Call the optimization routine VMCON at this scan point
    write(*,20) iscan, trim(xlabel), trim(vlabel), sweep(iscan_1), &
        trim(xlabel_2), trim(vlabel_2), sweep_2(iscan_R)
    ! 20     format(a,i2,a,4a,1pe10.3)
20  format('Starting scan point ', i3, ': ', a, ', ', a, ' = ', &
            1pe10.3, ' and ', a, ', ', a, ' = ', 1pe10.3)
  end subroutine scan_2d_write_point_header

  subroutine scan_2d_store_output(ifail, iscan_1, iscan_R, iscan, noutvars_, ipnscns_, outvar, &
    sweep_1_vals, sweep_2_vals)
    implicit none

    integer, intent(in) :: ifail
    integer, intent(in) :: iscan_1
    integer, intent(in) :: iscan_R
    integer, intent(in) :: iscan
    integer, intent(in) :: noutvars_, ipnscns_
    ! Required for shape of intent(out) arrays
    real(dp), dimension(noutvars_,ipnscns_), intent(out) :: outvar
    real(dp), dimension(ipnscns_), intent(out) :: sweep_1_vals, sweep_2_vals

    call scan_1d_store_output(iscan, ifail, noutvars_, ipnscns_, outvar)

    sweep_1_vals(iscan) = sweep(iscan_1)
    sweep_2_vals(iscan) = sweep_2(iscan_R)
  end subroutine scan_2d_store_output

  subroutine scan_2d_write_plot(iscan, outvar, sweep_1_vals, sweep_2_vals)
    use constants, only: nplot
    use global_variables, only: icase, xlabel, xlabel_2
    implicit none

    integer, intent(inout) :: iscan
    real(dp), dimension(:,:), intent(in) :: outvar
    real(dp), dimension(:), intent(in) :: sweep_1_vals, sweep_2_vals

    integer :: ivar
    character(len=48) :: tlabel
    character(len=25), dimension(noutvars), save :: plabel

    plabel( 1) = 'Ifail____________________'
    plabel( 2) = 'Sqsumsq__________________'
    plabel( 3) = 'Electric_cost_(mil/kwh)__'
    plabel( 4) = 'Capital_cost_(mil/kwh)___'
    plabel( 5) = 'Fuel_cost_(mil/kwh)______'
    plabel( 6) = 'Operations_cost_(mil/kwh)'
    plabel( 7) = 'Capital_cost_(millions)__'
    plabel( 8) = 'Core_costs_(millions)____'
    plabel( 9) = 'Direct_cost_(billions)___'
    plabel(10) = 'Major_Radius_(m)_________'
    plabel(11) = 'Aspect_Ratio_____________'
    plabel(12) = 'Plasma_Current_(MA)______'
    plabel(13) = 'B_Toroidal_Axis_(T)______'
    plabel(14) = 'B_total_on_axis_(T)______'
    plabel(15) = 'Safety_Factor____________'
    plabel(16) = 'qlim_(zero_if_ishape=0)__'
    plabel(17) = 'Beta_____________________'
    plabel(18) = 'Beta_Limit_______________'
    plabel(19) = 'Epsilon_Beta_Poloidal____'
    plabel(20) = 'Dens.weight_Te_(10keV)___'
    plabel(21) = 'Average_Dens_(10^20/m^3)_'
    plabel(22) = 'H-fact_Iter_Power________'
    plabel(23) = 'H-fact_Iter_Offset_______'
    plabel(24) = 'Fusion_Power_(MW)________'
    plabel(25) = 'nb_Fusion_Power_(MW)_____'
    plabel(26) = 'Wall_Load_(MW/m^2)_______'
    plabel(27) = 'Injection_Power_(MW)_____'
    plabel(28) = 'Inject_Pwr_Wall_Plug_(MW)'
    plabel(29) = 'Heating_Power_(MW)_______'
    plabel(30) = 'Current_Drive_(MW)_______'
    plabel(31) = 'Big_Q____________________'
    plabel(32) = 'Bootstrap_Fraction_______'
    plabel(33) = 'Neutral_Beam_Energy_(MeV)'
    plabel(34) = 'Divertor_Heat_(MW/m^2)___'
    plabel(35) = 'TF_coil_Power_(MW)_______'
    plabel(36) = 'TF_coil_weight_(kg)______'
    plabel(37) = 'vM_stress_in_TF_cond_(Pa)'
    plabel(38) = 'J_TF_inboard_leg_(MA/m^2)'
    plabel(39) = 'Centrepost_max_T_(TART)__'
    plabel(40) = 'Res_TF_inbrd_leg_Pwr_(MW)'
    plabel(41) = 'Coolant_Fraction_Ctr.____'
    plabel(42) = 'C/P_coolant_radius_(m)___'
    plabel(43) = 'C/P_coolant_velocity(m/s)'
    plabel(44) = 'C/P_pump_power_(MW)______'
    plabel(45) = 'PF_coil_Power_(MW)_______'
    plabel(46) = 'PF_coil_weight_(kg)______'
    plabel(47) = 'Gross_Elect_Pwr_(MW)_____'
    plabel(48) = 'Net_electric_Pwr_(MW)____'
    plabel(49) = 'Recirculating_Fraction___'
    plabel(50) = 'Psep/R___________________'
    plabel(51) = '' !OBSOLETE
    plabel(52) = 'Tot._radiation_power_(MW)'
    plabel(53) = 'First_wall_peak_temp_(K)_'
    plabel(54) = 'Cu_frac_TFC_conductor____'
    plabel(55) = 'Winding_pack_area_TFC(m2)'
    plabel(56) = 'Conductor_area_TFC_(m2)__'
    plabel(57) = 'Area_TF_inboard_leg_(m2)_'
    plabel(58) = 'Taup/taueff_lower_limit__'
    plabel(59) = 'Plasma_temp_at_sep_[keV]_'
    plabel(60) = 'SOL_density_at_OMP_______'
    plabel(61) = 'Power_through__separatrix'
    plabel(62) = 'neomp/nesep______________'
    plabel(63) = 'qtargettotal_____________'
    plabel(64) = 'Total_pressure_at_target_'
    plabel(65) = 'Temperature_at_target____'
    plabel(66) = 'Helium_fraction__________'
    plabel(67) = 'Momentum_loss_factor_____'
    plabel(68) = 'totalpowerlost_[W]_______'
    plabel(69) = 'H__concentration_________'
    plabel(70) = 'He_concentration_________'
    plabel(71) = 'Be_concentration_________'
    plabel(72) = 'C__concentration_________'
    plabel(73) = 'N__concentration_________'
    plabel(74) = 'O__concentration_________'
    plabel(75) = 'Ne_concentration_________'
    plabel(76) = 'Si_concentration_________'
    plabel(77) = 'Ar_concentration_________'
    plabel(78) = 'Fe_concentration_________'
    plabel(79) = 'Ni_concentration_________'
    plabel(80) = 'Kr_concentration_________'
    plabel(81) = 'Xe_concentration_________'
    plabel(82) = 'W__concentration_________'
    plabel(83) = 'teped____________________'
    plabel(84) = 'Max_field_on_TF_coil_____'

    tlabel = icase

  end subroutine scan_2d_write_plot

  subroutine scan_select(nwp, swp, iscn, vlab, xlab)
    !! Routine to select first scan case
    !! author: J Morris, UKAEA, Culham Science Centre
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	use build_variables, only: blnkoth, shldith, scrapli, scraplo, ohcth
    use constraint_variables, only: fiooic, walalw, bmxlim, fqval, taulimit, &
        gammax, tbrnmn, tbrmin, fjprot, pnetelin, powfmax
	use cost_variables, only: cfactr, iavail, fkind, startupratio
	use current_drive_variables, only: bscfmax, etaech
	use divertor_variables, only: hldivlim
	use error_handling, only: idiags, report_error
    use fwbs_variables, only: inlet_temp_liq, outlet_temp_liq, blpressure_liq, &
        n_liq_recirc, bz_channel_conduct_liq, pnuc_fw_ratio_dcll, f_nuc_pow_bz_struct, pitch
	use impurity_radiation_module, only: fimp, coreradius, impurity_arr_frac
    use physics_variables, only: kappa, dnbeta, te, aspect, ftar, bt, &
        rad_fraction_sol, triang, rmajor, beamfus0, hfact
    use numerics, only: epsvmc, boundu, boundl
    use tfcoil_variables, only: tmargmin_tf, sig_tf_case_max, n_pancake, oacdcp, &
      n_layer, b_crit_upper_nbti, sig_tf_wp_max
    use heat_transport_variables, only: crypmw_max, etath
    use rebco_variables, only: copperaoh_m2_max
    use pfcoil_variables, only: coheof, ohhghf, oh_steel_frac
    use CS_fatigue_variables, only: n_cycle_min, t_crack_vertical
    implicit none

    ! Arguments
    integer, intent(in) :: nwp, iscn
    real(dp), intent(in), dimension(:) :: swp
    character(len=25), intent(out) :: vlab, xlab

    select case (nwp)
        ! Use underscores instead of spaces in xlabel
        ! MDK Remove the "=" from vlabel, to make it easier to compare with
        ! list of iteration variables

        case (1)
            aspect = swp(iscn)
            vlab = 'aspect' ; xlab = 'Aspect_ratio'
        case (2)
            hldivlim = swp(iscn)
            vlab = 'hldivlim' ; xlab = 'Div_heat_limit_(MW/m2)'
        case (3)
            pnetelin = swp(iscn)
            vlab = 'pnetelin' ; xlab = 'Net_electric_power_(MW)'
        case (4)
            hfact = swp(iscn)
            vlab = 'hfact' ; xlab = 'Confinement_H_factor'
        case (5)
            oacdcp = swp(iscn)
            vlab = 'oacdcp' ; xlab = 'TF_inboard_leg_J_(MA/m2)'
        case (6)
            walalw = swp(iscn)
            vlab = 'walalw' ; xlab = 'Allow._wall_load_(MW/m2)'
        case (7)
            beamfus0 = swp(iscn)
            vlab = 'beamfus0' ; xlab = 'Beam_bkgrd_multiplier'
        case (8)
            fqval = swp(iscn)
            vlab = 'fqval' ; xlab = 'Big_Q_f-value'
        case (9)
            te = swp(iscn)
            vlab = 'te' ; xlab = 'Electron_temperature_keV'
        case (10)
            boundu(15) = swp(iscn)
            vlab = 'boundu(15)' ; xlab = 'Volt-second_upper_bound'
        case (11)
            dnbeta = swp(iscn)
            vlab = 'dnbeta' ; xlab = 'Beta_coefficient'
        case (12)
            bscfmax = swp(iscn)
            vlab = 'bscfmax' ; xlab = 'Bootstrap_fraction'
        case (13)
            boundu(10) = swp(iscn)
            vlab = 'boundu(10)' ; xlab = 'H_factor_upper_bound'
        case (14)
            fiooic = swp(iscn)
            vlab = 'fiooic' ; xlab = 'TFC_Iop_/_Icrit_f-value'
        case (15)
            fjprot = swp(iscn)
            vlab = 'fjprot' ; xlab = 'TFC_Jprot_limit_f-value'
        case (16)
            rmajor = swp(iscn)
            vlab = 'rmajor' ; xlab = 'Plasma_major_radius_(m)'
        case (17)
            bmxlim = swp(iscn)
            vlab = 'bmxlim' ; xlab = 'Max_toroidal_field_(T)'
        case (18)
            gammax = swp(iscn)
            vlab = 'gammax' ; xlab = 'Maximum_CD_gamma'
        case (19)
            boundl(16) = swp(iscn)
            vlab = 'boundl(16)' ; xlab = 'CS_thickness_lower_bound'
        case (20)
            tbrnmn = swp(iscn)
            vlab = 'tbrnmn' ; xlab = 'Minimum_burn_time_(s)'
        case (21)
            ! sigpfalw = swp(iscn)
            vlab = 'obsolete' ; xlab = 'obsolete'
        case (22)
            if (iavail == 1) call report_error(95)
            cfactr = swp(iscn)
            vlab = 'cfactr' ; xlab = 'Plant_availability_factor'
        case (23)
            boundu(72) = swp(iscn)
            vlab = 'boundu(72)' ; xlab = 'Ip/Irod_upper_bound'
        case (24)
            powfmax = swp(iscn)
            vlab = 'powfmax' ; xlab = 'Fusion_power_limit_(MW)'
        case (25)
            kappa = swp(iscn)
            vlab = 'kappa' ; xlab = 'Plasma_elongation'
        case (26)
            triang = swp(iscn)
            vlab = 'triang' ; xlab = 'Plasma_triangularity'
        case (27)
            tbrmin = swp(iscn)
            vlab = 'tbrmin' ; xlab = 'Min_tritium_breed._ratio'
        case (28)
            bt = swp(iscn)
            vlab = 'bt' ; xlab = 'Tor._field_on_axis_(T)'
        case (29)
            coreradius = swp(iscn)
            vlab = 'coreradius' ; xlab = 'Core_radius'
        case (30)
            !fimpvar = swp(iscn)
            vlab = 'OBSOLETE' ; xlab = 'OBSOLETE'
        case (31)
            taulimit = swp(iscn)
            vlab = 'taulimit' ; xlab = 'Taup/taueff_lower_limit'
        case (32)
            epsvmc = swp(iscn)
            vlab = 'epsvmc' ; xlab = 'VMCON error tolerance'
        case (33, 34, 35, 36, 37, 47)
            write(*,*) 'Kallenbach model has been removed, remove the kallenbach scan variables'
            stop 1
        case (38)
            boundu(129) = swp(iscn)
            vlab = 'boundu(129)' ; xlab = ' Neon upper limit'
        case (39)
            boundu(131) = swp(iscn)
            vlab = 'boundu(131)' ; xlab = ' Argon upper limit'
        case (40)
            boundu(135) = swp(iscn)
            vlab = 'boundu(135)' ; xlab = ' Xenon upper limit'
        case (41)
            blnkoth = swp(iscn)
            vlab = 'blnkoth' ; xlab = 'Outboard blanket thick.'
        case (42)
            fimp(9) = swp(iscn)
            impurity_arr_frac(9) = fimp(9)
            vlab = 'fimp(9)' ; xlab = 'Argon fraction'
        case (43)
            ! rho_ecrh = swp(iscn)
            vlab = 'obsolete' ; xlab = 'obsolete'
        case (44)
            sig_tf_case_max = swp(iscn)
            vlab = 'sig_tf_case_max' ; xlab = 'Allowable_stress_in_tf_coil_case_Tresca_(pa)'
        case (45)
            tmargmin_tf = swp(iscn)
            vlab = 'tmargmin_tf' ; xlab = 'Minimum_allowable_temperature_margin'
        case (46)
            boundu(152) = swp(iscn)
            vlab = 'boundu(152)' ; xlab = 'Max allowable fgwsep'
        case (48)
            n_pancake = int(swp(iscn))
            vlab = 'n_pancake' ; xlab = 'TF Coil - n_pancake'
        case (49)
            n_layer = int(swp(iscn))
            vlab = 'n_layer' ; xlab = 'TF Coil - n_layer'
        case (50)
            fimp(13) = swp(iscn)
            impurity_arr_frac(13) = fimp(13)
            vlab = 'fimp(13)' ; xlab = 'Xenon fraction'
        case (51)
            ftar = swp(iscn)
            vlab = 'ftar' ; xlab = 'lower_divertor_power_fraction'
        case (52)
            rad_fraction_sol = swp(iscn)
            vlab = 'rad_fraction_sol' ; xlab = 'SoL radiation fraction'
        case (53)
            boundu(157) = swp(iscn)
            vlab = 'boundu(157)' ; xlab = 'Max allowable fvssu'
        case (54)
            b_crit_upper_nbti = swp(iscn)
            vlab = 'Bc2(0K)' ; xlab = 'GL_NbTi Bc2(0K)'
        case(55)
            shldith = swp(iscn)
            vlab = 'shldith' ; xlab = 'Inboard neutronic shield'
        case(56)
            crypmw_max = swp(iscn)
            vlab = 'crypmw_max' ; xlab = 'max allowable crypmw'
        case(57)
            boundl(2) = swp(iscn)
            vlab = 'boundl(2)' ; xlab = 'bt minimum'
        case(58)
            scrapli = swp(iscn)
            vlab = 'scrapli' ; xlab = 'Inboard FW-plasma sep gap'
        case(59)
            scraplo = swp(iscn)
            vlab = 'scraplo' ; xlab = 'Outboard FW-plasma sep gap'
        case (60)
            sig_tf_wp_max = swp(iscn)
            vlab = 'sig_tf_wp_max' ; xlab = 'Allowable_stress_in_tf_coil_conduit_Tresca_(pa)'
        case (61)
            copperaoh_m2_max = swp(iscn)
            vlab = 'copperaoh_m2_max' ; xlab = 'Max CS coil current / copper area'
        case (62)
            coheof = swp(iscn)
            vlab = 'coheof' ; xlab = 'CS coil current density at EOF (A/m2)'
        case (63)
            ohcth = swp(iscn)
            vlab = 'ohcth' ; xlab = 'CS coil thickness (m)'
        case (64)
            ohhghf = swp(iscn)
            vlab = 'ohhghf' ; xlab = 'CS height (m)'
        case (65)
            n_cycle_min = swp(iscn)
            vlab = 'n_cycle_min' ; xlab = 'CS stress cycles min'
        case (66)
            oh_steel_frac = swp(iscn)
            vlab = 'oh_steel_frac' ; xlab = 'CS steel fraction'
        case (67)
            t_crack_vertical = swp(iscn)
            vlab = 't_crack_vertical' ; xlab = 'Initial crack vertical size (m)'
        case (68)
            inlet_temp_liq = swp(iscn)
            vlab = 'inlet_temp_liq' ; xlab = 'Inlet Temperature Liquid Metal Breeder/Coolant (K)'
        case (69)
            outlet_temp_liq = swp(iscn)
            vlab = 'outlet_temp_liq' ; xlab = 'Outlet Temperature Liquid Metal Breeder/Coolant (K)'
        case(70)
            blpressure_liq = swp(iscn)
            vlab = 'blpressure_liq' ; xlab = 'Blanket liquid metal breeder/coolant pressure (Pa)'
        case(71)
            n_liq_recirc = swp(iscn)
            vlab = 'n_liq_recirc' ; xlab = 'Selected number of liquid metal breeder recirculations per day'
        case(72)
            bz_channel_conduct_liq = swp(iscn)
            vlab = 'bz_channel_conduct_liq'  ; xlab = 'Conductance of liquid metal breeder duct walls (A V-1 m-1)'
        case(73)
            pnuc_fw_ratio_dcll = swp(iscn)
            vlab = 'pnuc_fw_ratio_dcll' ; xlab = 'Ratio of FW nuclear power as fraction of total (FW+BB)'
        case(74)
            f_nuc_pow_bz_struct = swp(iscn)
            vlab = 'f_nuc_pow_bz_struct' ; xlab = 'Fraction of BZ power cooled by primary coolant for dual-coolant balnket'
        case(75)
            pitch = swp(iscn)
            vlab = 'pitch' ; xlab = 'pitch of first wall cooling channels (m)'
        case (76)
            etath = swp(iscn)
              vlab = 'etath' ; xlab = 'Thermal conversion eff.'
        case (77)
            startupratio = swp(iscn)
              vlab = 'startupratio' ; xlab = 'Gyrotron redundancy'
        case (78)
            fkind = swp(iscn)
              vlab = 'fkind' ; xlab = 'Multiplier for Nth of a kind costs'
        case (79)
            etaech = swp(iscn)
              vlab = 'etaech' ; xlab = 'ECH wall plug to injector efficiency'
        case default
            idiags(1) = nwp ; call report_error(96)

    end select

  end subroutine scan_select

  subroutine post_optimise(ifail)
  !! Called after calling the optimising equation solver from Python.
  !! author: P J Knight, CCFE, Culham Science Centre
  !! ifail   : input integer : error flag
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code

  use constraints
  use error_handling
  use function_evaluator
  use numerics
  use process_output
  use utilities, only:upper_case
  use main_module, only:verror
  ! for ipedestal = 2 option
  use global_variables, only: convergence_parameter
  use constants, only: iotty, nout, mfile
  use physics_variables, only: ipedestal
  use define_iteration_variables, only: boundxc, loadxc
  implicit none

  !  Arguments
  integer, intent(in) :: ifail

  !  Local variables
  integer :: ii,inn,iflag
  real(dp) :: summ,xcval,xmaxx,xminn,f,xnorm
  real(dp), dimension(ipeqns) :: con1, con2, err
  character(len=1), dimension(ipeqns) :: sym
  character(len=10), dimension(ipeqns) :: lab
  character(len=30) :: strfom
  character(len=60) :: string1, string2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Check on accuracy of solution by summing the
  !  squares of the residuals of the equality constraints
  summ = 0.0D0
  do ii = 1,neqns
     summ = summ + rcm(ii)*rcm(ii)
  end do
  sqsumsq = sqrt(summ)

  !  Turn on error reporting
  errors_on = .true.

  !  Print out information on solution
  call oheadr(nout,'Numerics')
  call ocmmnt(nout,'PROCESS has performed a VMCON (optimisation) run.')
  if (ifail /= 1) then
     !call ocmmnt(nout,'but could not find a feasible set of parameters.')
    !  call oheadr(nout,'PROCESS COULD NOT FIND A FEASIBLE SOLUTION')
    !  call ovarin(iotty,'VMCON error flag (ifail)','',ifail)
     call ovarin(nout,'VMCON error flag','(ifail)',ifail)
     call oheadr(iotty,'PROCESS COULD NOT FIND A FEASIBLE SOLUTION')
     call oblnkl(iotty)

     idiags(1) = ifail ; call report_error(132)

  else
     call ocmmnt(nout,'and found a feasible set of parameters.')
     call oblnkl(nout)
     call ovarin(nout,'VMCON error flag','(ifail)',ifail)
     call oheadr(iotty,'PROCESS found a feasible solution')
  end if

  !call oblnkl(nout)

  !  If necessary, write out a relevant error message
  if (ifail /= 1) then
     call verror(ifail)
     call oblnkl(nout)
     call oblnkl(iotty)
  else
     !  Show a warning if the constraints appear high even if allegedly converged
     if (sqsumsq >= 1.0D-2) then
        call oblnkl(nout)
        call ocmmnt(nout,'WARNING: Constraint residues are HIGH; consider re-running')
        call ocmmnt(nout,'   with lower values of EPSVMC to confirm convergence...')
        call ocmmnt(nout,'   (should be able to get down to about 1.0E-8 okay)')
        call oblnkl(nout)
        call ocmmnt(iotty,'WARNING: Constraint residues are HIGH; consider re-running')
        call ocmmnt(iotty,'   with lower values of EPSVMC to confirm convergence...')
        call ocmmnt(iotty,'   (should be able to get down to about 1.0E-8 okay)')
        call oblnkl(iotty)

        fdiags(1) = sqsumsq ; call report_error(134)

     end if
  end if

  call ovarin(nout,'Number of iteration variables','(nvar)',nvar)
  call ovarin(nout,'Number of constraints (total)','(neqns+nineqns)',neqns+nineqns)
  call ovarin(nout,'Optimisation switch','(ioptimz)',ioptimz)
  call ovarin(nout,'Figure of merit switch','(minmax)',minmax)
!   if (ifail /= 1) then
!      call ovarin(nout,'VMCON error flag','(ifail)',ifail)
!   end if

  call ovarre(nout,'Square root of the sum of squares of the constraint residuals','(sqsumsq)',sqsumsq, 'OP ')
  call ovarre(nout,'VMCON convergence parameter','(convergence_parameter)',convergence_parameter, 'OP ')
  call ovarre(nout,'Normalised objective function','(norm_objf)',norm_objf, 'OP ')
  call ovarin(nout,'Number of VMCON iterations','(nviter)',nviter, 'OP ')
  call oblnkl(nout)

  if (ifail == 1) then
     string1 = 'PROCESS has successfully optimised the iteration variables'
  else
     string1 = 'PROCESS has tried to optimise the iteration variables'
  end if

  if (minmax > 0) then
     string2 = ' to minimise the figure of merit: '
  else
     string2 = ' to maximise the figure of merit: '
  end if

  strfom = lablmm(abs(minmax))
  call upper_case(strfom)
  write(nout,10) trim(string1) // trim(string2),  trim(strfom)
10 format(a90, t92, a22)

  call oblnkl(nout)

  !  Check which variables are at bounds
  iflag = 0
  do ii = 1,nvar
     xminn = 1.01D0*bondl(ii)
     xmaxx = 0.99D0*bondu(ii)

     if (xcm(ii) < xminn) then
        if (iflag == 0) then
           call ocmmnt(nout, &
                'Certain operating limits have been reached,')
           call ocmmnt(nout, &
                'as shown by the following iteration variables that are')
           call ocmmnt(nout, &
                'at or near to the edge of their prescribed range :')
           call oblnkl(nout)
           iflag = 1
        end if
        xcval = xcm(ii)*scafc(ii)
        !write(nout,30) ii,lablxc(ixc(ii)),xcval,bondl(ii)*scafc(ii)
        write(nout,30) lablxc(ixc(ii)),xcval,bondl(ii)*scafc(ii)
     end if

     if (xcm(ii) > xmaxx) then
        if (iflag == 0) then
           call ocmmnt(nout, &
                'Certain operating limits have been reached,')
           call ocmmnt(nout, &
                'as shown by the following iteration variables that are')
           call ocmmnt(nout, &
                'at or near to the edge of their prescribed range :')
           call oblnkl(nout)
           iflag = 1
        end if
        xcval = xcm(ii)*scafc(ii)
        write(nout,40) lablxc(ixc(ii)),xcval,bondu(ii)*scafc(ii)
     end if
  end do

!30 format(t4,'Variable ',i3,' (',a9, &
!        ',',1pe12.4,') is at or below its lower bound:',1pe12.4)
30 format(t4, a30, '=',1pe12.4,' is at or below its lower bound:',1pe12.4)
40 format(t4, a30, '=',1pe12.4,' is at or above its upper bound:',1pe12.4)
!40 format(t4,'Variable ',i3,' (',a9, &
!        ',',1pe12.4,') is at or above its upper bound:',1pe12.4)

  !  Print out information on numerics
  call osubhd(nout,'The solution vector is comprised as follows :')
!  write(nout,50)
! Remove Lagrange multipliers as no-one understands them.
! MFILE not changed
!50 format(t47,'lower',t59,'upper')

  write(nout,60)
!60 format(t23,'final',t33,'fractional',t46,'Lagrange',t58,'Lagrange')
60 format(t43,'final',t55,'final /')


  write(nout,70)
!70 format(t5,'i',t23,'value',t35,'change',t45,'multiplier', &
!        t57,'multiplier')
70 format(t5,'i',t43,'value',t55,'initial')

  call oblnkl(nout)

  do inn = 1,nvar
     xcs(inn) = xcm(inn)*scafc(inn)
!     write(nout,80) inn,lablxc(ixc(inn)),xcs(inn),xcm(inn), &
!          vlam(neqns+nineqns+inn), vlam(neqns+nineqns+1+inn+nvar)
     write(nout,80) inn,lablxc(ixc(inn)),xcs(inn),xcm(inn)
!80 format(t2,i4,t8,a9,t19,4(1pe12.4))
!80 format(t2,i4,t8,a30,t39,2(1pe12.4))
80 format(t2,i4,t8,a30,t39,1pe12.4, t52, 0pf10.4)
! MDK The 0p is needed because of a bizarre "feature"/bug in fortran:
! the 1p in the previous format continues until changed.
     call ovarre(mfile,lablxc(ixc(inn)),'(itvar'//int_to_string3(inn)//')',xcs(inn))

     !  'Range-normalised' iteration variable values for MFILE:
     !  0.0 (at lower bound) to 1.0 (at upper bound)
     if (bondl(inn) == bondu(inn)) then
        xnorm = 1.0D0
     else
        xnorm = (xcm(inn) - bondl(inn)) / (bondu(inn) - bondl(inn))
        xnorm = max(xnorm, 0.0D0)
        xnorm = min(xnorm, 1.0D0)
     end if
     ! Added ratio final/initial to MFILE
     call ovarre(mfile,trim(lablxc(ixc(inn)))//' (final value/initial value)', &
          '(xcm'//int_to_string3(inn)//')',xcm(inn))
     call ovarre(mfile,trim(lablxc(ixc(inn)))//' (range normalised)', &
          '(nitvar'//int_to_string3(inn)//')',xnorm)
  end do


  call osubhd(nout, &
       'The following equality constraint residues should be close to zero :')

  call constraint_eqns(neqns+nineqns,-1,con1,con2,err,sym,lab)
  write(nout,90)
90 format(t48,'physical',t73,'constraint',t100,'normalised')
  write(nout,100)
100 format(t47,'constraint',t74,'residue',t101,'residue')
  call oblnkl(nout)
  do inn = 1,neqns
     write(nout,110) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
          lab(inn),err(inn),lab(inn),con1(inn)
     call ovarre(mfile,lablcc(icc(inn))//' normalised residue', &
          '(eq_con'//int_to_string3(icc(inn))//')',con1(inn))
  end do
110 format(t2,i4,t8,a33,t46,a1,t47,1pe12.4,t60,a10,t71,1pe12.4,t84,a10,t98,1pe12.4)

  if (nineqns > 0) then
     call osubhd(nout, &
          'The following inequality constraint residues should be greater than or approximately equal to zero :')

     do inn = neqns+1,neqns+nineqns
        !write(nout,120) inn,lablcc(icc(inn)),rcm(inn),vlam(inn)
        write(nout,110) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
                        lab(inn), err(inn), lab(inn)
        call ovarre(mfile,lablcc(icc(inn)),'(ineq_con'//int_to_string3(icc(inn))//')',rcm(inn))
     end do
  end if

! 120 format(t2,i4,t8,a33,t45,1pe12.4,1pe12.4)

end subroutine post_optimise

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module scan_module
