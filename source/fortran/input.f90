 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Uncomment #define line below to perform unit testing
!  Compile using pre-processor, e.g. ifort -cpp input.f90
!#define unit_test

module process_input

  !! Module containing the routines that perform the actual reading
  !! and parsing of the input file
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module provides a set of routines to read in data from the
  !! main PROCESS input file (IN.DAT). The format of the file is
  !! similar to the F90 NAMELIST structure, but with a few
  !! additional features:
  !! <OL>
  !! <P><LI>Comments can be read in that are copied to the standard
  !! output channel - these are lines with five (or more)
  !! consecutive '*' characters at the start.
  !! <P><LI>Other lines within the file can contain simple comments
  !! for the user - these are not copied to the standard output
  !! channel. They start with one to four '*' characters.
  !! </OL>
  !! <P>Character strings, integers and double precision values can
  !! be read in.
  !! <P>The following rules must be obeyed when writing an input
  !! file:
  !! <UL>
  !! <P><LI>Each variable must be on a separate line.
  !! <P><LI>Leading spaces are ignored.
  !! <P><LI>Variable names can be upper case, lower case, or a
  !! mixture of both.
  !! <P><LI>Spaces may not appear within a variable name or data
  !! value.
  !! <P><LI>Other spaces within a line, and trailing spaces, are
  !! ignored.
  !! <P><LI>Commas are not necessary between variables.
  !! <P><LI>Data can extend over more than one line.
  !! <P><LI>One-dimensional arrays can be explicitly subscripted, or
  !! unscripted, in which case the following element order is
  !! assumed: A(1), A(2), A(3), ...
  !! <P><LI>At present, multiple dimension arrays can only be
  !! handled without reference to explicit subscripts, in which case
  !! the following element order is assumed: B(1,1), B(2,1), B(3,1),
  !! etc. The use of the input file to specify multiple dimension
  !! array elements is prone to error.
  !! <P><LI>Unscripted array elements must be separated by commas.
  !! <P><LI>Blank lines are allowed anywhere in the input file.
  !! <P><LI>Lines starting with a * are assumed to be comments.
  !! <P><LI>Comment lines starting with five or more asterisks
  !! (i.e. *****) are reproduced verbatim in the output file. These
  !! should be used copiously to give a great deal of information
  !! about the run being performed, and should be updated before
  !! every single run of the code, as it is very easy to lose track
  !! of what is being attempted.
  !! </UL>
  !! A User's Guide to the PROCESS Systems Code, P. J. Knight,
  !! AEA Fusion Report AEA FUS 251, 1993
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  integer, parameter :: nin = 10

  integer, parameter :: maxlen = 2000  !  maximum line length
  character(len=maxlen) :: line  !  current line of text from input file
  integer :: linelen, lineno  !  current line length, line number
  integer :: iptr             !  current position on line
  integer :: infile, outfile, report_changes, icode
  logical :: subscript_present
  logical :: error
  character(len=78) :: error_message

  ! Vars for subroutine input() requiring re-initialisation before each new run
  integer :: show_changes
  logical :: constraints_exist

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_input
    !! Initialise module variables
    implicit none

    error = .False.
    show_changes = 0
    constraints_exist = .false.
    line = ""
    linelen = 0
    lineno = 0
    iptr = 0
    infile = 0
    outfile = 0
    report_changes = 0
    icode = 0
    subscript_present = .false.
    error_message = ""
  end subroutine init_input

  subroutine parse_input_file(in_file,out_file,show_changes)

    !! Routine that parses the contents of the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: J Morris, CCFE, Culham Science Centre
    !! author: F Warmer, IPP Greifswald
    !! in_file  : input integer : Fortran input unit identifier
    !! out_file : input integer : Fortran output unit identifier
    !! show_changes : input integer : switch to turn on (1) or off (0)
    !! reporting of changed values
    !! This routine reads the data from the PROCESS input file (IN.DAT),
    !! dealing with comments or blank lines correctly, and sets the
    !! value of any variables found in the file. Any changes
    !! from the default values may be reported if required.
    !! <P>Each possible variable in this block is dealt with
    !! individually. (To add additional input variables, simply copy
    !! and edit one of the similar existing examples.)
    !! The routine also does the extremely useful function of checking
    !! that the given value for a variable lies within a sensible
    !! predefined range, and stops the program if apparently
    !! nonsensical values are attempted.
    !! A User's Guide to the PROCESS Systems Code, P. J. Knight,
    !! AEA Fusion Report AEA FUS 251, 1993
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: dcopper, dalu
    use global_variables, only: run_tests, verbose, maxcal
    use build_variables, only: tf_in_cs, blbmoth, blbuith, dr_shld_outboard, &
      shldtth, shldlth, vgap_vv_thermalshield, plleni, dr_fw_outboard, dr_shld_blkt_gap, &
      dr_shld_thermal_inboard, dr_shld_thermal_outboard, thshield_vb, i_cs_precomp, &
      blbpith, aplasmin, blbuoth, dr_tf_inboard, &
      iohcl, dr_tf_shld_gap, f_z_cryostat, dr_bore, plleno, dr_fw_plasma_gap_inboard, gapomin, dr_cryostat, &
      rinboard, dr_blkt_outboard, fseppc, plsepo, dr_blkt_inboard, &
      dr_cs, plsepi, blbmith, dr_cs_tf_gap, fcspc, dr_fw_plasma_gap_outboard, vgaptop, &
      blbpoth, dr_shld_vv_gap_inboard, dr_fw_inboard, vgap_xpoint_divertor, dr_shld_inboard, sigallpc, tfootfi, f_avspace,&
      r_cp_top, dr_vv_inboard, dr_vv_outboard, d_vv_top, d_vv_bot, f_r_cp, i_r_cp_top
    use buildings_variables, only: hcwt, conv, wgt, trcl, rbwt, &
      esbldgm3, fndt, row, wgt2, pibv, dz_tf_cryostat, stcl, clh2, &
      tfcbv, hccl, rbrt, triv, shov, admv, i_bldgs_v, i_bldgs_size, &
      mbvfac, pfbldgm3, wsvfac, rbvfac, rxcl, shmf, &
      aux_build_l, aux_build_w, aux_build_h, auxcool_l, auxcool_w, auxcool_h, &
      bioshld_thk, chemlab_l, chemlab_w, chemlab_h, control_buildings_l, &
      control_buildings_w, control_buildings_h, crane_arm_h, crane_clrnc_h, &
      crane_clrnc_v, cryomag_l, cryomag_w, cryomag_h, cryostore_l, &
      cryostore_w, cryostore_h, cryostat_clrnc, elecdist_l, elecdist_w, &
      elecdist_h, elecload_l, elecload_w, elecload_h, elecstore_l, &
      elecstore_w, elecstore_h, fc_building_l, fc_building_w, &
      gas_buildings_l, gas_buildings_w, gas_buildings_h, ground_clrnc, &
      hcd_building_l, hcd_building_w, hcd_building_h, hw_storage_l, &
      hw_storage_w, hw_storage_h, heat_sink_l, heat_sink_w, heat_sink_h, &
      hot_sepdist, hotcell_h, ilw_smelter_l, ilw_smelter_w, ilw_smelter_h, &
      ilw_storage_l, ilw_storage_w, ilw_storage_h, llw_storage_l, &
      llw_storage_w, llw_storage_h, magnet_pulse_l, magnet_pulse_w, &
      magnet_pulse_h, magnet_trains_l, magnet_trains_w, magnet_trains_h, &
      maint_cont_l, maint_cont_w, maint_cont_h, nbi_sys_l, nbi_sys_w, &
      qnty_sfty_fac, reactor_clrnc, reactor_fndtn_thk, reactor_hall_l, &
      reactor_hall_w, reactor_hall_h, reactor_roof_thk, reactor_wall_thk, &
      robotics_l, robotics_w, robotics_h, sec_buildings_l, sec_buildings_w, &
      sec_buildings_h, staff_buildings_h, staff_buildings_area, &
      transp_clrnc, turbine_hall_l, turbine_hall_w, turbine_hall_h, &
      tw_storage_l, tw_storage_w, tw_storage_h, warm_shop_l, warm_shop_w, &
      warm_shop_h, water_buildings_l, water_buildings_w, water_buildings_h, &
      workshop_l, workshop_w, workshop_h
    use constraint_variables, only: fl_h_threshold, fpeakb, fpsep, fdivcol, ftcycl, &
      beta_poloidal_max, fpsepbqar, ftmargtf, fradwall, fptfnuc, fnesep, fportsz, tbrmin, &
      maxradwallload, pseprmax, fdene, fniterpump, fpinj, pnetelin, powfmax, &
      fgamcd, ftbr, mvalim, f_alpha_energy_confinement_min, walalw, fmva, fradpwr, nflutfmax, fipir, &
      fauxmn, fiooic,fr_conducting_wall, fjohc0, frminor, psepbqarmax, ftpeak, bigqmin, &
      fstrcond, fptemp, ftmargoh, fvs, fbeta_max, vvhealw, fpnetel, ft_burn, &
      ffuspow, fpsepr, ptfnucmax, fvdump, pdivtlim, falpha_energy_confinement, nbshinefmax, &
      fcqt, fzeffmax, fstrcase, fhldiv, foh_stress, fwalld, gammax, fjprot, &
      ft_current_ramp_up, tcycmn, auxmin, zeffmax, f_fw_rad_max, fdtmp, fpoloidalpower, &
      fnbshinef, freinke, fvvhe, fqval, fq, fmaxvvstress, fbeta_poloidal, fbeta_poloidal_eps, fjohc, &
      fflutf, bmxlim, t_burn_min, fbeta_min, fecrh_ignition, fstr_wp, fncycle
    use cost_variables, only: ucich, uctfsw, dintrt, ucblbe, uubop, dtlife, &
      cost_factor_vv, cfind, uccry, fcap0cp, uccase, uuves, cconshtf, conf_mag, &
      ucbllipb, ucfuel, uumag, ucpfbs, ireactor, uucd, div_umain_time, div_nu, &
      maintenance_gen, uctfps, uufw, tbktrepl, cost_factor_fwbs, decomf, &
      cconshpf, uche3, ucpfdr1, ucech, uudiv, cost_model, adivflnc, &
      cost_factor_rh, cost_factor_bop, ifueltyp, fcontng, fwbs_nref, &
      cost_factor_buildings, favail, cconfix, ucblli2o, abktflnc, ucf1, ucfnc, &
      ucpfps, iavail, ibkt_life, life_dpa, ucpfbk, cost_factor_tf_coils, costexp_pebbles, &
      ucmisc, cpstflnc, uccryo, costexp, fwbs_nu, ucpfic, ucblbreed, tcomrepl, uufuel, &
      ucdiv, uccpcl1, discount_rate, uctfbr, uccpclb, ucoam, div_prob_fail, ucnbi, &
      uccu, ucwst, cfactr, div_nref, amortization, ucwindtf, ucme, csi, cowner, &
      cost_factor_misc, fcr0, cturbb, lsa, fcap0, output_costs, &
      cost_factor_land, redun_vacp, ucrb, uctfbus, num_rh_systems, fkind, &
      fwbs_umain_time, uchrs, avail_min, uciac, ucshld, tdivrepl, &
      ucblli, ucpfcb, tlife, ipnet, fcdfuel, ucbus, ucpfb, uchts, &
      maintenance_fwbs, fwbs_prob_fail, uclh, ucblss, ucblvd, ucsc, ucturb, &
      ucpens, cland, ucwindpf, i_cp_lifetime, cplife_input, &
      startupratio, tmain, u_unplanned_cp, supercond_cost_model
    use current_drive_variables, only: pinjfixmw, etaech, pinjalw, etanbi, &
      f_tritium_beam, gamma_ecrh, pheat, beamwd, beam_energy, pheatfix, bootstrap_current_fraction_max, &
      forbitloss, nbshield, tbeamin, feffcd, iefrf, iefrffix, irfcd, cboot, &
      etalh, frbeam, harnum, xi_ebw, wave_mode
    use divertor_variables, only: fdfs, anginc, divdens, divclfr, c4div, &
      c5div, ksic, fififi, flux_exp, divplt, delld, c2div, beta_div, betao, divdum, tdiv, c6div, &
      omegan, prn1, frrp, xpertin, c1div, betai, bpsout, xparain, fdiva, &
      zeffdiv, hldivlim, rlenmax, divfix, c3div, &
      hldiv, i_hldiv
    use fwbs_variables, only: fblhebpo, vfblkt, fdiv, fvolso, i_fw_coolant_type, &
      dx_fw_module, i_blanket_type, blktmodel, afwi, fblli2o, nphcdin, breeder_multiplier, &
      fw_armour_thickness, roughness, fwclfr, breedmat, fblli, fblvd, &
      iblanket_thickness, vfcblkt, breeder_f, fbllipb, fhcd, vfshld, fblhebmi, &
      f_neut_shield, fw_th_conductivity, nblktmodti, dr_fw_wall, afwo, &
      fvolsi, etahtp, nblktmodpo, pres_fw_coolant, emult, temp_fw_coolant_out, nblktmodpi, &
      fblhebpi, fblss, inlet_temp, outlet_temp, fblbreed, qnuc, blpressure, &
      blpressure_liq, n_liq_recirc, pnuc_fw_ratio_dcll, f_nuc_pow_bz_struct, &
      declblkt, fblhebmo, blkttype, radius_fw_channel, inuclear, declshld, hcdportsize, &
      npdiv, f_fw_peak, primary_pumping, dr_pf_cryostat, secondary_cycle, secondary_cycle_liq, &
      denstl, declfw, nphcdout, i_blkt_inboard, vfpblkt, temp_fw_coolant_in, wallpf, fblbe, &
      fhole, fwbsshape, coolp, temp_fw_max, irefprop, len_fw_channel, &
      li6enrich, etaiso, nblktmodto, fvoldw, i_shield_mat, i_bb_liq, &
      icooldual, ifci, inlet_temp_liq, outlet_temp_liq, bz_channel_conduct_liq, ipump, ims, &
      coolwh, emult
    use heat_transport_variables, only: htpmw_fw, baseel, fmgdmw, htpmw_div, &
      pwpm2, etath, vachtmw, iprimshld, fpumpdiv, pinjmax, htpmw_blkt, etatf, &
      htpmw_min, fpumpblkt, ipowerflow, htpmw_shld, fpumpshld, trithtmw, &
      fpumpfw, crypmw_max, f_crypmw
    use ife_variables, only: bldzu, etali, sombdr, gainve, cdriv0, v1dzl, &
      bldrc, fauxbop, pfusife, dcdrv0, fwdr, pdrive, mcdriv, ucconc, shdr, &
      v3dzu, bldzl, rrin, maxmat, shmatf, fwmatf, drveff, flirad, shdzu, v2dzu, &
      pifecr, ifedrv, v2dr, chmatf, v1dr, v1matf, dcdrv1, chdzu, dcdrv2, &
      ifetyp, fwdzl, htpmw_ife, uccarb, v3matf, fbreed, edrive, ptargf, cdriv2, &
      fburn, fwdzu, etave, v3dr, uctarg, shdzl, ucflib, v3dzl, v1dzu, v2dzl, &
      chdzl, chrad, cdriv1, tgain, somtdr, v2matf, rrmax, bldr, frrmax, &
      blmatf, ife
    use impurity_radiation_module, only: coreradius, n_impurities, &
      coreradiationfraction, fimp
    use numerics, only: factor, minmax, neqns, nvar, epsfcn, ixc, &
      ftol, ipnvars, nineqns, ipeqns, boundu, icc, ipnfoms, name_xc
    use pfcoil_variables, only: rhopfbus, rjconpf, zref, fcuohsu, oh_steel_frac, vf, &
      coheof, sigpfcalw, alstroh, ipfres, fcupfsu, fvssu, etapsu, i_cs_stress, &
      fbmaxcs, ngc, rpf2, fcohbop, ohhghf, vfohc, isumatoh, ngrpmx, ngc2, rpf1, &
      ngrp, isumatpf, nfxfh, alfapf, routr, sigpfcf, pfclres, bmaxcs_lim, &
      ncls, nfixmx, cptdin, ipfloc, i_sup_pf_shape, rref, i_pf_current, &
      ccl0_ma, ccls_ma, ld_ratio_cst
    use physics_variables, only: ipedestal, taumax, i_single_null, fvsbrnni, &
      rhopedt, f_vol_plasma, f_deuterium, ffwal, i_beta_component, itartpf, i_l_h_threshold, &
      fpdivlim, beta_poloidal_eps_max, i_confinement_time, kappa95, aspect, f_r_conducting_wall, nesep, c_beta, csawth, dene, &
      ftar, plasma_res_factor, f_sync_reflect, f_nd_beam_electron, beta, neped, hfact, beta_norm_max, &
      fgwsep, rhopedn, tratio, q0, i_plasma_geometry, i_plasma_shape, fne0, ignite, f_tritium, &
      i_beta_fast_alpha, tauee_in, alphaj, alphat, i_plasma_current, q95, ti, tesep, ind_plasma_internal_norm, triang, &
      itart, f_nd_alpha_electron, iprofile, triang95, rad_fraction_sol, betbm0, f_nd_protium_electrons, &
      teped, f_helium3, iwalld, ejima_coeff, f_alpha_plasma, fgwped, tbeta, i_bootstrap_current, &
      i_rad_loss, te, alphan, rmajor, plasma_square, kappa, fkzohm, beamfus0, &
      tauratio, i_density_limit, bt, i_plasma_wall_gap, n_confinement_scalings, beta_max, beta_min, &
      i_diamagnetic_current, i_pfirsch_schluter_current, m_s_limit, burnup_in
    use pf_power_variables, only: iscenr, maxpoloidalpower
    use pulse_variables, only: i_pulsed_plant, dtstor, itcycl, istore, bctmp

    use primary_pumping_variables, only: t_in_bb, t_out_bb, dp_he, p_he, gamma_he, &
      dp_fw_blkt, dp_fw, dp_blkt, dp_liq

    use scan_module, only: isweep_2, nsweep, isweep, scan_dim, nsweep_2, &
      sweep_2, sweep, ipnscns, ipnscnv
    use stellarator_variables, only: f_asym, isthtr, n_res, iotabar, fdivwet, &
      f_w, bmn, shear, m_res, f_rad, flpitch, istell, max_gyrotron_frequency, &
      te0_ecrh_achievable
    use tfcoil_variables, only: fcoolcp, tfinsgap, vftf, &
      fhts, dr_tf_wp, rcool, rho_tf_leg, thkcas, &
      casthi, n_pancake, bcritsc, i_tf_sup, str_pf_con_res, thwcndut, &
      thicndut, tftmp, oacdcp, tmax_croco, ptempalw, tmargmin_tf, tmpcry, &
      sig_tf_case_max, dztop, dcond, str_cs_con_res, etapump, drtop, vcool, dcondins, &
      i_tf_tresca, dhecoil, tmaxpro, n_tf_coils, temp_cp_average, fcutfsu, j_tf_bus, &
      casthi_fraction, tmargmin_cs, vdalw, dcase, t_turn_tf,&
      cpttf_max, tdmptf, casths, i_tf_turns_integer, quench_model, &
      tcritsc, layer_ins, tinstf, n_layer, tcoolin, ripmax, frhocp, &
      cpttf, tmargmin, casths_fraction, eff_tf_cryo, eyoung_ins, &
      eyoung_steel, eyoung_res_tf_buck, eyoung_cond_axial, f_vforce_inboard, &
      f_a_tf_cool_outboard, frholeg, ftoroidalgap, i_tf_sc_mat, i_tf_shape, i_tf_bucking, &
      n_tf_graded_layers, n_tf_joints, n_tf_joints_contact, poisson_al, &
      poisson_copper, poisson_steel, rho_tf_joints, rho_tf_bus, th_joint_contact,&
      i_tf_stress_model, eyoung_al, i_tf_wp_geom, i_tf_case_geom, &
      i_tf_turns_integer, n_rad_per_layer, b_crit_upper_nbti, t_crit_nbti, &
      i_cp_joints, n_tf_turn, f_t_turn_tf, t_turn_tf_max, t_cable_tf, &
      sig_tf_wp_max, eyoung_cond_trans, i_tf_cond_eyoung_axial, i_tf_cond_eyoung_trans, &
      str_wp_max, str_tf_con_res, i_str_wp, max_vv_stress, theta1_coil, theta1_vv, &
      len_tf_bus

    use times_variables, only: t_current_ramp_up, pulsetimings, t_ramp_down, t_fusion_ramp, t_precharge, t_burn, &
      t_between_pulse, tohsin
    use vacuum_variables, only: dwell_pump, pbase, tn, pumpspeedfactor, &
      initialpressure, outgasfactor, prdiv, pumpspeedmax, rat, outgasindex, &
      pumpareafraction, ntype, vacuum_model, pumptp
    use rebco_variables, only: hastelloy_thickness, f_coppera_m2, &
      rebco_thickness, tape_thickness, tape_width, &
      copper_rrr, coppera_m2_max, croco_thick, copper_thick, f_copperaoh_m2, &
      copperaoh_m2, copperaoh_m2_max
    use reinke_variables, only: reinke_mode, fzactual, impvardiv, lhat
    use water_usage_variables, only: airtemp, watertemp, windspeed
    use CS_fatigue_variables, only: residual_sig_hoop, t_crack_radial, t_crack_vertical, &
      t_structural_vertical, t_structural_radial, n_cycle_min, bkt_life_csf, sf_vertical_crack, &
      sf_radial_crack, sf_fast_fracture ,paris_coefficient, paris_power_law, walker_coefficient, &
      fracture_toughness

    implicit none

    !  Arguments

    integer, intent(in) :: in_file, out_file, show_changes

    !  Local variables

    integer :: iost
    integer :: isub1,isub2,varlen
    integer :: no_constraints
    integer :: no_iteration
    integer :: foundAst

    character(len=32) :: varnam

    logical :: obsolete_var

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialise local variables
    no_constraints = 0
    no_iteration = 0
    obsolete_var = .false.

    !  Initialise module-wide variables

    infile = in_file
    outfile = out_file
    report_changes = show_changes

    icode = 0
    lineno = 0

    !  Main loop

    loop_over_lines: do

       subscript_present = .FALSE.

       read(infile,'(A)',iostat=iost) line

       !  On error or end, return
       if (iost /= 0) exit loop_over_lines

       lineno = lineno + 1

       line = adjustl(line)  !  rotate any leading blanks to the end
       linelen = len_trim(line)


20     continue

       !  Ignore blank lines

       if (line == ' ') cycle

       !  Ignore comments, unless they start with '*****',
       !  in which case print them.

       if (line(1:5) == '*****') write(outfile,*) line(1:76)
       if (line(1:1) == '*') cycle
       if (line(1:1) == '$') cycle  !  in case block delimiters are still present

       iptr = 1

       !Ignore input comments denoted by asterisk, before assigning variables

       if (index(line,'*') > 0) then
          foundAst = index(line,'*') - 1
          linelen = min(linelen, foundAst)
          line = line(:linelen)
       end if

       !  This must be an assignment line, so get the variable name

       call get_variable_name(varnam,varlen,isub1,isub2)
       if (isub1 /= 0) subscript_present = .TRUE.
       if (varlen == 0) then
          write(*,*) 'Error in IN.DAT at line ', lineno
          write(*,*) line
          error = .True.
       end if

       !  Read the associated data

       variable: select case (varnam(1:varlen))

       case ('dcond')
          call parse_real_array('dcond', dcond, isub1, 9, &
               'TF/PF coil superconductor density (kg/m3)', icode)




       case ('quench_model')
          call parse_string_variable('quench_model', quench_model, &
          'Switch for TF coil quench model (Only applies to REBCO magnet at present)')

       case ('cptdin')
          call parse_real_array('cptdin', cptdin, isub1, ngc2, &
               'Current per turn for PF coil', icode)

       case ('ipfloc')
          call parse_int_array('ipfloc', ipfloc, isub1, ngrpmx, &
               'PF coil location', icode)

       ! New optional argument startindex used MDK 3/3/17
       ! Allows simplified IN.DAT format for icc and ixc.
       case ('icc')
          no_constraints = no_constraints + 1
          call parse_int_array('icc', icc, isub1, ipeqns, &
               'Constraint equation', icode,no_constraints)
          no_constraints = isub1
      case ('ixc')
          no_iteration = no_iteration + 1
          call parse_int_array('ixc', ixc, isub1, ipnvars, &
                   'Iteration variable', icode,no_iteration)
          no_iteration = isub1

       case ('ncls')
          call parse_int_array('ncls', ncls, isub1, ngrpmx+2, &
               'No of coils in PF group', icode)

       case ('nfxfh')
          call parse_int_variable('nfxfh', nfxfh, 1, nfixmx/2, &
               'Central Solenoid splitting parameter')
       case ('ngrp')
          call parse_int_variable('ngrp', ngrp, 0, ngrpmx, &
               'No of groups of PF coils')

       case ('rref')
          call parse_real_array('rref', rref, isub1, ngrpmx, &
               'radius of location 4 coil groups, minor radii from major radius', icode)
       case ('ccl0_ma')
          call parse_real_array('ccl0_ma', ccl0_ma, isub1, ngrpmx, &
               'Flux-swing cancel current of PF coil groups, MA', icode)
       case ('ccls_ma')
          call parse_real_array('ccls_ma', ccls_ma, isub1, ngrpmx, &
               'Equilibrium current of PF coil groups, MA', icode)
       case ('vf')
          call parse_real_array('vf', vf, isub1, ngc2, &
               'Void fraction of PF coil', icode)

       case ('zref')
          call parse_real_array('zref', zref, isub1, ngrpmx, &
               'height of location 3 and 4 coil groups / minor radius', icode)


       case ('i_fw_coolant_type')
          call parse_string_variable('i_fw_coolant_type', i_fw_coolant_type, 'first wall coolant')
          call lower_case(i_fw_coolant_type)



       case ('i_blanket_type')
          call parse_int_variable('i_blanket_type', i_blanket_type, 1, 5, 'Switch for blanket model')
          if ((i_blanket_type == 2).or.(i_blanket_type == 4)) then
            write(outfile,*) ' '
            write(outfile,*) '**********'
            write(outfile,*) 'i_blanket_type = 2/4, KIT HCPB/HCLL model has been removed -'
            write(outfile,*) 'please select a different blanket model.'
            write(outfile,*) '**********'
            write(outfile,*) ' '
            obsolete_var = .true.
          endif

          if (i_blanket_type == 3) then
              dr_fw_inboard = 0.03D0
              dr_fw_outboard = 0.03D0
              fw_armour_thickness = 0.003D0
          end if


       case ('iblanket_thickness')
          call parse_int_variable('iblanket_thickness', iblanket_thickness, 1, 3, &
               'Blanket thickness switch')
          if (iblanket_thickness == 1) then
            dr_blkt_inboard = 0.53D0
            dr_blkt_outboard = 0.91D0
          else if (iblanket_thickness == 2) then
            dr_blkt_inboard = 0.64D0
            dr_blkt_outboard = 1.11D0
          else if (iblanket_thickness == 3) then
            dr_blkt_inboard = 0.75D0
            dr_blkt_outboard = 1.30D0
          end if

       case ('cfind')
          call parse_real_array('cfind', cfind, isub1, 4, &
               'Indirect cost factor vs LSA', icode)





       case ('uchts')
          call parse_real_array('uchts', uchts, isub1, 2, &
               'Cost of heat transp system equip per loop ($/W)', icode)

       case ('ucoam')
          call parse_real_array('ucoam', ucoam, isub1, 4, &
               'Annual cost of operation and maintenance', icode)

       case ('ucsc')
          call parse_real_array('ucsc', ucsc, isub1, 9, &
               'Cost of superconductor ($/kg)', icode)
       case ('ucturb')
          call parse_real_array('ucturb', ucturb, isub1, 2, &
               'Cost of turbine plant equipment ($)', icode)
       case ('ucwst')
          call parse_real_array('ucwst', ucwst, isub1, 4, &
               'cost of waste disposal (M$/yr)', icode)

       case ('isweep')
          call parse_int_variable('isweep', isweep, 0, ipnscns, &
               'Number of scans to perform')
       case ('nsweep')
          call parse_int_variable('nsweep', nsweep, 1, ipnscnv, &
               'Variable used in scan')
       case ('sweep')
          call parse_real_array('sweep', sweep, isub1, ipnscns, &
               'Actual values to use in scan', icode)

        case ('isweep_2')
          call parse_int_variable('isweep_2', isweep_2, 0, ipnscns, &
               'Number of 2D scans to perform')
        case ('nsweep_2')
          call parse_int_variable('nsweep_2', nsweep_2, 1, ipnscnv, &
               'Second variable used in 2D scan')
        case ('sweep_2')
          call parse_real_array('sweep_2', sweep_2, isub1, ipnscns, &
               'Actual values to use in 2D scan', icode)



       case ('vacuum_model')
          call parse_string_variable('vacuum_model', vacuum_model, 'vacuum_model')



       case ('impvardiv')
          call parse_int_variable('impvardiv', impvardiv, 3, n_impurities, &
               'Index of impurity to be iterated for Reike criterion')



       case ('blmatf')  !  N.B. actually a 2-D array
          call parse_real_array('blmatf', blmatf, isub1, 3*(maxmat+1), &
                    'IFE blanket material fraction', icode)



       case ('chmatf')
          call parse_real_array('chmatf', chmatf, isub1, maxmat+1, &
                    'IFE chamber material fraction', icode)

       case ('etave')
          call parse_real_array('etave', etave, isub1, 10, &
                    'IFE driver efficiency vs driver energy', icode)



       case ('fwmatf')  !  N.B. actually a 2-D array
          call parse_real_array('fwmatf', fwmatf, isub1, 3*(maxmat+1), &
                    'IFE first wall material fraction', icode)
       case ('gainve')
          call parse_real_array('gainve', gainve, isub1, 10, &
                    'IFE target gain vs driver energy', icode)




       case ('shmatf')  !  N.B. actually a 2-D array
          call parse_real_array('shmatf', shmatf, isub1, 3*(maxmat+1), &
                    'IFE shield material fraction', icode)



       case ('v1matf')  !  N.B. actually a 2-D array
          call parse_real_array('v1matf', v1matf, isub1, 3*(maxmat+1), &
                    'IFE void 1 material fraction', icode)
       case ('v2matf')  !  N.B. actually a 2-D array
          call parse_real_array('v2matf', v2matf, isub1, 3*(maxmat+1), &
                    'IFE void 2 material fraction', icode)
       case ('v3matf')  !  N.B. actually a 2-D array
          call parse_real_array('v3matf', v3matf, isub1, 3*(maxmat+1), &
                    'IFE void 3 material fraction', icode)

       case default
          error_message = 'Unknown variable in input file: '//varnam(1:varlen)
          write(*,*) error_message
          write(*,*) 'Error occurred at this line in the IN.DAT file:', lineno
          write(*,*) line
          error = .True.

       end select variable

       !  Uncomment the following to abort the code if an obsolete variable name
       !  has been found in the input file

       if (obsolete_var) then
          error_message = 'Obsolete variable specified'
          write(*,*) error_message
          write(*,*) 'Error occurred at this line in the IN.DAT file: ', lineno
          write(*,*) line
          obsolete_var = .false.
          error = .True.
       end if

       !  If we have just read in an array, a different loop-back is needed

       if (icode == -1) goto 20

       cycle

    end do loop_over_lines

    if(neqns == 0) then
        ! The value of neqns has not been set in the input file.  Default = 0.
        neqns = no_constraints - nineqns
    else
        ! The value of neqns has been set in the input file.
        nineqns = no_constraints - neqns
    end if

    nvar = no_iteration

    if (error .eqv. .True.) stop 1

    ! MDK Try allocating here
    if (allocated(name_xc)) deallocate(name_xc)
    allocate(name_xc(nvar))
    ! Ensure array is initialised
    name_xc = ""

  end subroutine parse_input_file

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parse_real_variable(varnam,varval,vmin,vmax,description)

    !! Routine that obtains the value of a real variable from the input
    !! file and checks that it lies within the expected range
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : input string : name of the variable
    !! varval : input/output real : value of the variable
    !! vmin : input real : minimum allowed value for the variable
    !! vmax : input real : maximum allowed value for the variable
    !! description : input string : brief description of the variable
    !! This routine parses a line containing a 'name = value' pair
    !! for a real variable, extracting the value from the line
    !! and checking whether it lies between user-defined lower and
    !! upper limits.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: varnam, description
    real(dp), intent(inout) :: varval
    real(dp), intent(in) :: vmin, vmax

    !  Local variables

    real(dp) :: oldval

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check whether a subscript was found by the preceding call to GET_VARIABLE_NAME
    !  and stop if this is the case

    if (subscript_present) then
       write(*,*) 'Unexpected subscript found at line ', lineno
       write(*,*) 'Variable name and description:'
       write(*,*) varnam, ', ', description
          error = .True.
    end if

    !  Obtain the new value for the variable

    oldval = varval

    call get_value_real(varval,icode)

    if (icode /= 0) then
       write(*,*) 'Error whilst reading input file.  Variable name and description:'
       write(*,*) varnam, ', ', description
       write(*,*) 'Comments should be indicated by an asterisk'
       error = .True.
    end if

    !  Check variable lies within range

    call check_range_real(varnam,varval,vmin,vmax)

    if ((report_changes == 1).and.(varval /= oldval)) then
       write(outfile,*) trim(description),', ',trim(varnam),' = ',varval
    end if

  end subroutine parse_real_variable

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parse_int_variable(varnam,varval,vmin,vmax,description)

    !! Routine that obtains the value of an integer variable from the
    !! input file and checks that it lies within the expected range
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : input string : name of the variable
    !! varval : input/output integer : value of the variable
    !! vmin : input integer : minimum allowed value for the variable
    !! vmax : input integer : maximum allowed value for the variable
    !! description : input string : brief description of the variable
    !! This routine parses a line containing a 'name = value' pair
    !! for an integer variable, extracting the value from the line
    !! and checking whether it lies between user-defined lower and
    !! upper limits.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: nout
    implicit none

    !  Arguments

    character(len=*), intent(in) :: varnam, description
    integer, intent(inout) :: varval
    integer, intent(in) :: vmin, vmax

    !  Local variables

    integer :: oldval

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check whether a subscript was found by the preceding call to GET_VARIABLE_NAME
    !  and stop if this is the case

    if (subscript_present) then
       write(*,*) 'Unexpected subscript found in IN.DAT at line number: ', lineno
       write(*,*) 'Name and description of variable: '
       write(*,*) varnam, description
       error = .True.
    end if

    !  Obtain the new value for the variable

    oldval = varval
    call get_value_int(varval,icode)
    if (icode /= 0) then
       write(*,*) 'Error found in input file, check line ',lineno
       write(*,*) 'Variable name, description:'
       write(*,*) varnam, ', ', description
          error = .True.
    end if

    !  Check variable lies within range

    call check_range_int(varnam,varval,vmin,vmax)

    if ((report_changes == 1).and.(varval /= oldval)) then
       write(outfile,*) trim(description),', ',trim(varnam),' = ',varval
    end if

  end subroutine parse_int_variable

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parse_string_variable(varnam,varval,description)

    !! Routine that obtains the value of a string variable from the
    !! input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : input string : name of the variable
    !! varval : input/output string : value of the variable
    !! description : input string : brief description of the variable
    !! This routine parses a line containing a 'name = value' pair
    !! for a string variable, extracting the value from the line.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: varnam, description
    character(len=*), intent(inout) :: varval

    !  Local variables

    character(len=maxlen) :: oldval

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check whether a subscript was found by the preceding call to GET_VARIABLE_NAME
    !  and stop if this is the case

    if (subscript_present) then
       write(*,*) 'Unexpected subscript found in IN.DAT at line number: ', lineno
       write(*,*) 'Name and description of variable: '
       write(*,*) varnam, description
       error = .True.
       !stop 1
    end if

    !  Obtain the new value for the variable

    oldval = varval
    call get_substring(varval,icode)
    if (icode /= 0) then
       write(*,*) 'Error in IN.DAT found at line ',lineno
       write(*,*) 'Variable name, description:'
       write(*,*) varnam, ', ', description
       error = .True.
    end if

    if ((report_changes == 1).and.(trim(varval) /= trim(oldval))) then
       write(outfile,*) trim(description),', ',trim(varnam),' = ',varval
    end if

  end subroutine parse_string_variable

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parse_real_array(varnam,varval,isub1,n,description,icode)

    !! Routine that obtains the values of a real array from the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : input string : name of the variable
    !! varval(n) : input/output real array : value of the variable
    !! isub1 : input integer : array element pointer
    !! n : input integer : size of varval array
    !! icode : output integer : diagnostic flag
    !! description : input string : brief description of the variable
    !! This routine parses a line in one of the two following forms:
    !! <PRE>
    !! name = v1[, v2, ...]
    !! name(element) = v
    !! </PRE>
    !! to read in and extract one or more values for a real 1-D array.
    !! <P>N.B. No array bounds or value range checking is performed.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: varnam, description
    integer, intent(inout) :: isub1
    integer, intent(in) :: n
    integer, intent(out) :: icode
    real(dp), dimension(n), intent(inout) :: varval

    !  Local variables

    real(dp) :: oldval, val

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check whether a subscript was found by the preceding call to GET_VARIABLE_NAME

    if (subscript_present) then

       oldval = varval(isub1)
       call get_value_real(val,icode)

       if (icode /= 0) then
          write(*,*) 'Error in IN.DAT found at line ',lineno
          write(*,*) 'Variable name, description:'
          write(*,*) varnam, ', ', description
          error = .True.
       end if

       varval(isub1) = val
       if ((report_changes == 1).and.(varval(isub1) /= oldval)) then
          write(outfile,10) trim(description),', ', &
               trim(varnam),'(',isub1,') = ',varval(isub1)
       end if

    else

       isub1 = 1
       do
          call get_value_real(val,icode)
          !  icode == 1 denotes an error
          !  icode == -1 denotes end of line, so the next line needs to be read in
          !  (hence the 'goto 20' in the calling routine)
          if (icode /= 0) return

          oldval = varval(isub1)
          varval(isub1) = val
          if ((report_changes == 1).and.(varval(isub1) /= oldval)) then
             write(outfile,10) trim(description),', ', &
                  trim(varnam),'(',isub1,') = ',varval(isub1)
          end if
          isub1 = isub1 + 1
       end do
    end if

10  format(a,a,a,a1,i3,a,e14.6e2)

  end subroutine parse_real_array

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parse_int_array(varnam,varval,isub1,n,description,icode,startindex)

    !! Routine that obtains the values of an integer array
    !! from the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : input string : name of the variable
    !! varval(n) : input/output integer array : value of the variable
    !! isub1 : input integer : array element pointer
    !! n : input integer : size of varval array
    !! icode : output integer : diagnostic flag
    !! description : input string : brief description of the variable
    !! This routine parses a line in one of the two following forms:
    !! <PRE>
    !! name = v1[, v2, ...]
    !! name(element) = v
    !! </PRE>
    !! to read in and extract one or more values for an integer 1-D array.
    !! <P>N.B. No array bounds or value range checking is performed.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments
    character(len=*), intent(in) :: varnam, description
    integer, intent(inout) :: isub1
    integer, intent(in) :: n
    integer, intent(out) :: icode
    integer, dimension(n), intent(inout) :: varval
    integer, intent(in), optional :: startindex

    !  Local variables
    integer :: oldval, val
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check whether a subscript was found by the preceding call to GET_VARIABLE_NAME

    if (subscript_present) then

       oldval = varval(isub1)
       call get_value_int(val,icode)

       if (icode /= 0) then
          write(*,*) 'Error in IN.DAT found at line ',lineno
          write(*,*) 'Variable name, description:'
          write(*,*) varnam, ', ', description
          error = .True.
       end if

       varval(isub1) = val
       if ((report_changes == 1).and.(varval(isub1) /= oldval)) then
          write(outfile,10) trim(description),', ', &
               trim(varnam),'(',isub1,') = ',varval(isub1)
       end if

   else  ! subscript is not present

       isub1 = 1
       if(present(startindex))isub1 = startindex
       do
          call get_value_int(val,icode)

          !  icode == 1 denotes an error
          !  icode == -1 denotes end of line
          if (icode /= 0) then
              ! Make sure isub1 = the last array index
              isub1 = isub1 - 1
              return
          end if

          oldval = varval(isub1)
          varval(isub1) = val
          if ((report_changes == 1).and.(varval(isub1) /= oldval)) then
             write(outfile,10) trim(description),', ', &
                  trim(varnam),'(',isub1,') = ',varval(isub1)
          end if
          isub1 = isub1 + 1
       end do

    end if

10  format(a,a,a,a1,i3,a,i12)

  end subroutine parse_int_array

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine string_to_int(string,length,ivar,icode)

    !! Routine that converts the ASCII digits in a string to
    !! an integer
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : input string : contains digits of the number
    !! length : input integer : useful length of character string
    !! ivar : output integer : value stored in the string
    !! icode : output integer : diagnostic flag
    !! This routine converts the ASCII digits in string(1:length)
    !! to the integer ivar. It is equivalent to doing
    !! 'READ(STRING(1:LENGTH),I) IVAR' but this routine conforms
    !! to the ANSI standard.
    !! Each digit is parsed in turn, the current total is multiplied
    !! by ten and the new digit is added.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: string
    integer, intent(in) :: length
    integer, intent(out) :: ivar, icode

    !  Local variables

    character(len=maxlen) :: xstr
    integer :: iptr,izero,xlen
    logical :: negate

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ivar = 0
    icode = 0

    if (length <= 0) goto 1000

    negate = .false.
    izero = ichar('0')
    iptr = 1
    xstr = string(1:length)

    ! *** Ignore trailing spaces

    xlen = len_trim(xstr)
    if (xlen <= 0) goto 1000

    ! *** Ignore leading spaces

10  continue
    if (xstr(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr > xlen) goto 1000
       goto 10
    end if

    ! *** Check for leading + or -

    if (xstr(iptr:iptr) == '+') then
       iptr = iptr + 1
       if (iptr > xlen) goto 1000
    else if (xstr(iptr:iptr) == '-') then
       negate = .true.
       iptr = iptr + 1
       if (iptr > xlen) goto 1000
    else
       continue
    end if

    ! *** Ignore leading zeros

20  continue
    if (xstr(iptr:iptr) == '0') then
       iptr = iptr + 1
       if (iptr > xlen) goto 1000
       goto 20
    end if

    ! *** Check for number too large

    if ((xlen-iptr+1) > 10) then
       if (negate) then
          ivar = -1234567890
       else
          ivar = 1234567890
          write(*,*) '1 Problem with IN file, please check line'
          write(*,*) xstr
          error = .True.
       end if
       icode = 1
       goto 1000
    else if ((xlen-iptr+1) == 10) then
       if (xstr(iptr:xlen) > '2147483647') then
          if (negate) then
             ivar = -1234567890
          else
             ivar = 1234567890
          end if
          icode = 1
          goto 1000
       end if
    else
       continue
    end if

    ! *** Parse the digits

30  continue
    if ((xstr(iptr:iptr) >= '0').and.(xstr(iptr:iptr) <= '9')) then
       ivar = (ivar * 10) + (ichar(xstr(iptr:iptr))-izero)
       iptr = iptr + 1
       if (iptr <= xlen) goto 30

       ! *** This is the normal exit path...

       if (negate) ivar = -ivar

    else
       if(ivar /= 0) then
          write(*,*) 'Problem with IN file, please check line'
          write(*,*) xstr
          write(*,*) 'Comments should be indicated by an asterisk (*)'
          error = .True.
       end if
       icode = 1
    end if

1000 continue

  end subroutine string_to_int

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine string_to_real(string,length,rval,icode)

    !! Routine that converts the ASCII digits in a string to
    !! a real value
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : input string : contains digits of the number
    !! length : input integer : useful length of character string
    !! rvar : output real : value stored in the string
    !! icode : output integer : diagnostic flag
    !! This routine converts the ASCII digits in string(1:length)
    !! to the real variable rvar.
    !! The string is parsed one character at a time, from the left,
    !! handling the mantissa, and all other components of the real
    !! number separately, combining them at the end.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: string
    integer, intent(in) :: length
    real(dp), intent(out) :: rval
    integer, intent(out) :: icode

    !  Local variables

    real(dp) :: valbdp,valadp,xexp
    integer :: iptr,izero,iexpon
    logical :: negatm,negate

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iptr = 1

    ! *** Ignore leading spaces

10  continue
    if (string(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= length) goto 10
    end if

    ! *** Initialise real value

    rval = 0.0D0

    ! *** ASCII '0'

    izero = ichar('0')

    ! *** If negative mantissa

    negatm = .false.

    ! *** If negative exponent

    negate = .false.

    ! *** Value before decimal point

    valbdp = 0.0D0

    ! *** Value after decimal point

    valadp = 0.0D0

    ! *** Exponent

    iexpon = 0

    ! *** First character can be +, -, ., or <digit>

    if (string(iptr:iptr) == '+') then
       iptr = iptr + 1
       if (iptr > length) goto 50
    else if (string(iptr:iptr) == '-') then
       iptr = iptr + 1
       if (iptr > length) goto 50
       negatm = .true.
    else
       continue
    end if

    ! *** Parse the mantissa - before the decimal point

    valbdp = 0.0D0
    xexp = -1.0D0
20  continue
    if ((string(iptr:iptr) >= '0').and.(string(iptr:iptr) <= '9')) then
       valbdp = (valbdp * 10.0D0) + dble(ichar(string(iptr:iptr))-izero)
       iptr = iptr + 1
       if (iptr > length) goto 50
       goto 20
    end if

    ! *** After the mantissa, we expect '.' or 'd' or 'e'

    if (string(iptr:iptr) == '.') then
       iptr = iptr + 1
       if (iptr > length) goto 50
    end if

    ! *** Parse the mantissa - after the decimal point

    valadp = 0.0D0
30  continue
    if ((string(iptr:iptr) >= '0').and.(string(iptr:iptr) <= '9')) then
       valadp = valadp + (dble(ichar(string(iptr:iptr))-izero)*(10.0D0 ** xexp))
       xexp = xexp - 1.0D0
       iptr = iptr + 1
       if (iptr > length) goto 50
       goto 30
    end if

    ! *** Now we expect the exponent

    if ( (string(iptr:iptr) == 'D').or. &
         (string(iptr:iptr) == 'E').or. &
         (string(iptr:iptr) == 'd').or. &
         (string(iptr:iptr) == 'e')) then
       iptr = iptr + 1
       if (iptr > length) goto 50

       ! *** First character can be +, -, ., or <digit>

       if (string(iptr:iptr) == '+') then
          iptr = iptr + 1
          if (iptr > length) goto 50
       else if (string(iptr:iptr) == '-') then
          iptr = iptr + 1
          if (iptr > length) goto 50
          negate = .true.
       else
          continue
       end if

       ! *** Parse the exponent

40     continue
       if ((string(iptr:iptr) >= '0').and.(string(iptr:iptr) <= '9')) then
          iexpon = (iexpon * 10) + (ichar(string(iptr:iptr))-izero)
          iptr = iptr + 1
          if (iptr <= length) goto 40
       else
          goto 60
       end if
    else
       goto 60
    end if

50  continue

    ! *** Negative exponent?

    if (negate) iexpon = -iexpon

    ! *** Build the number at last

    if (iexpon == 0) then
       rval = (valbdp + valadp)
    else
       rval = (valbdp + valadp) * (10.0D0 ** iexpon)
    end if

    ! *** Negative mantissa?

    if (negatm) rval = -rval

    ! *** All OK

    icode = 0
    goto 1000

    ! *** Errors

60  continue

    write(*,*) 'Problem with IN file, please check line'
    write(*,*) string
    write(*,*) 'Comments should be indicated by an asterisk (*)'
    error = .True.

    icode = 1

1000 continue

  end subroutine string_to_real

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_value_int(ival,icode)

    !! Routine that extracts an integer value from a line of the
    !! input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! ival   : output integer : extracted integer value
    !! icode  : output integer : diagnostic flag
    !! This routine extracts an integer value from the current line of
    !! the input file, i.e. the value of an integer variable as
    !! specified by the user.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(out) :: ival, icode

    !  Local variables

    character(len=maxlen) :: varval
    integer :: varlen
    integer :: foundComma, foundAst, foundPoint

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *** Ignore leading spaces

10  continue
    if (iptr <= linelen) then
       if (line(iptr:iptr) == ' ') then
          iptr = iptr + 1
          goto 10
       end if
    end if

    if (iptr > linelen) goto 60

! 40     continue !KE I guess I can remove this too

    ! *** Put rest of line into varval (makes it easier to parse)

    varval = line(iptr:)

    ! *** Exclude any input after * or , - these denote an input comment

    varlen = len_trim(varval)
    foundComma = varlen
    foundAst = varlen
    foundPoint = 0

    if (index(varval,',') > 0) then
       foundComma = index(varval,',') - 1
    end if
    if (index(varval,'*') > 0) then
       foundAst = index(varval,'*') - 1
    end if
    varlen = min(varlen, foundComma, foundAst)

    if (varlen <= 0) varlen = index(varval,' ') - 1
    if (varlen <= 0) varlen = iptr

    varval = varval(:varlen)

    varlen = len_trim(varval)

    foundPoint = index(varval,'.') - 1
    if (foundPoint > 0) then
       varlen = foundPoint
       write(*,*) 'Integer value expected in following input line...'
       write(*,*) '   ',line(1:50),'...'
       error = .True.
    end if

    ! *** Update pointer

    iptr = iptr + varlen

    ! *** Ignore trailing spaces

50  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= linelen) goto 50
    end if

    ! *** Ignore comma, if present

    if (iptr <= linelen) then
       if (line(iptr:iptr) == ',') iptr = iptr + 1
    end if

    ! *** Convert the ASCII text into an integer value

    call string_to_int(varval,varlen,ival,icode)

    goto 1000

60  continue
    icode = 1

1000 continue

  end subroutine get_value_int

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_value_real(rval,icode)

    !! Routine that extracts a real value from a line of the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rval   : output real : extracted real value
    !! icode  : output integer : diagnostic flag
    !! This routine extracts a real value from the current line of
    !! the input file, i.e. the value of a real variable as specified
    !! by the user.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(out) :: icode
    real(dp), intent(out) :: rval

    !  Local variables

    character(len=maxlen) :: varval
    integer :: varlen
    integer :: foundComma, foundAst

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *** Ignore leading spaces

10  continue
    if (iptr <= linelen) then
       if (line(iptr:iptr) == ' ') then
          iptr = iptr + 1
          goto 10
       end if
    end if
    if (iptr > linelen) goto 60

    ! *** Put rest of line into varval (makes it easier to parse)

    varval = line(iptr:)

    ! *** Exclude any input after * or , - these denote an input comment

    varlen = len_trim(varval)
    foundComma = varlen
    foundAst = varlen

    if (index(varval,',') > 0) then
       foundComma = index(varval,',') - 1
    end if
    if (index(varval,'*') > 0) then
       foundAst = index(varval,'*') - 1
    end if
    varlen = min(varlen, foundComma, foundAst)

    if (varlen <= 0) varlen = index(varval,' ') - 1
    if (varlen <= 0) varlen = iptr

    varval = varval(:varlen)

    varlen = len_trim(varval)

    ! *** Update pointer

    iptr = iptr + varlen

    ! *** Ignore trailing spaces

50  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= linelen) goto 50
    end if

    ! *** Ignore comma, if present

    if (iptr <= linelen) then
       if (line(iptr:iptr) == ',') iptr = iptr + 1
    end if

    ! *** Convert the ASCII text into a real value

    call string_to_real(varval,varlen,rval,icode)

    goto 1000

60  continue
    icode = 1

1000 continue

  end subroutine get_value_real

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_substring(string,icode)

    !! Routine that extracts a substring from a line of the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : output string : extracted string
    !! icode  : output integer : diagnostic flag
    !! This routine extracts a string from the current line of
    !! the input file, i.e. the value of a string variable as specified
    !! by the user. Unlike routine
    !! <A HREF="get_substring_trim.html">get_substring_trim</A>,
    !! this routine does not truncate the string found at its first
    !! non-leading blank.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(out) :: icode
    character(len=*), intent(out) :: string

    !  Local variables

    character(len=maxlen) :: varval
    integer :: varlen

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *** Ignore leading spaces

10  continue
    if (iptr <= linelen) then
       if (line(iptr:iptr) == ' ') then
          iptr = iptr + 1
          goto 10
       end if
    end if
    if (iptr > linelen) goto 60

    ! *** Put rest of line into varval (makes it easier to parse)

    varval = line(iptr:)
    varlen = len_trim(varval)

    if (varlen <= 0) varlen = iptr

    ! *** Update pointer

    iptr = iptr + varlen

    ! *** Ignore trailing spaces

50  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= linelen) goto 50
    end if

    ! *** Ignore comma, if present

    if (iptr <= linelen) then
       if (line(iptr:iptr) == ',') iptr = iptr + 1
    end if

    ! *** Write the text into the variable

    string = varval

    goto 1000

60  continue
    icode = 1

1000 continue

  end subroutine get_substring

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_subscript(isub1,isub2,icode)

    !! Routine that extracts any subscripts present in a line of
    !! the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! isub1  : output integer : first subscript found
    !! isub2  : output integer : second subscript found
    !! icode  : output integer : diagnostic flag
    !! This routine extracts any subscripts from the current line of
    !! the input file, i.e. if any array elements are specified
    !! by the user. It looks at the next non-space character in the
    !! line, and if it is a left bracket, it assumes that at
    !! least one subscript is to follow and extracts it/them.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(out) :: isub1, isub2, icode

    !  Local variables

    integer :: izero
    logical :: negate

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *** Initial values

    isub1 = 0
    isub2 = 0

    ! *** First character should be '('

    if (line(iptr:iptr) /= '(') goto 70
    iptr = iptr + 1
    if (iptr > linelen) goto 80

    ! *** Parse the first subscript
    ! *** Ignore leading spaces

10  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr > linelen) goto 80
       goto 10
    end if

    izero = ichar('0')
    negate = .false.

    ! *** Extract and evaluate the first subscript
    ! *** Subscript may be prefaced by '+' or '-'

    if (line(iptr:iptr) == '+') then
       iptr = iptr + 1
       if (iptr > linelen) goto 80
    else if (line(iptr:iptr) == '-') then
       negate = .true.
       iptr = iptr + 1
       if (iptr > linelen) goto 80
    else
       continue
    end if

20  continue

    if ((line(iptr:iptr) >= '0').and.(line(iptr:iptr) <= '9')) then
       isub1 = isub1 * 10 + ichar(line(iptr:iptr)) - izero
       iptr = iptr + 1
       if (iptr > linelen) goto 80
       goto 20
    end if
    if (negate) isub1 = -isub1

    ! *** Ignore trailing spaces of first subscript

30  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr > linelen) goto 70
       goto 30
    end if

    ! *** Is there a second subscript?

    if (line(iptr:iptr) == ',') then
       iptr = iptr + 1
       if (iptr > linelen) goto 80

       ! *** Ignore leading spaces of second subscript

40     continue
       if (line(iptr:iptr) == ' ') then
          iptr = iptr + 1
          if (iptr > linelen) goto 80
          goto 40
       end if

       ! *** Extract and evaluate the second subscript

       negate = .false.

       ! *** Subscript may be prefaced by '+' or '-'

       if (line(iptr:iptr) == '+') then
          iptr = iptr + 1
          if (iptr > linelen) goto 80
       else if (line(iptr:iptr) == '-') then
          negate = .true.
          iptr = iptr + 1
          if (iptr > linelen) goto 80
       else
          continue
       end if
50     continue
       if ((line(iptr:iptr) >= '0').and.(line(iptr:iptr) <= '9')) then
          isub2 = isub2 * 10 + ichar(line(iptr:iptr)) - izero
          iptr = iptr + 1
          if (iptr > linelen) goto 80
          goto 50
       end if

       ! *** Is it a negative subscript?

       if (negate) isub2 = -isub2

       ! *** Ignore trailing spaces of second subscript

60     continue
       if (line(iptr:iptr) == ' ') then
          iptr = iptr + 1
          if (iptr <= linelen) goto 60
       end if

    end if

    ! *** Must end with ')'

    if (line(iptr:iptr) /= ')') goto 80
    iptr = iptr + 1

70  continue
    icode = 0
    goto 1000

80  continue
    icode = 1

1000 continue

  end subroutine get_subscript

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_variable_name(varnam,varlen,isub1,isub2)

    !! Routine that extracts a variable name from a line of
    !! the input file
    !! author: P J Knight, CCFE, Culham Science Centre
    !! varnam : output string  : extracted variable name
    !! varlen : output integer : length of variable name
    !! isub1  : output integer : first subscript found
    !! isub2  : output integer : second subscript found
    !! This routine extracts a variable name from the current line of
    !! the input file. It also extracts any subscripts present.
    !! On exit, the counter <CODE>iptr</CODE> points to the first
    !! character of the value to be assigned to the variable.
    !! If the routine finds an error a value of 0 is returned in
    !! variable <CODE>varlen</CODE>.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(out) :: varlen, isub1, isub2
    character(len=*), intent(out) :: varnam

    !  Local variables

    character(len=maxlen) :: line1
    integer :: ifrom,ito,icode

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! *** Store LINE in local variable

    line1 = line

    ! *** Convert string to lower case

    call lower_case(line)

    varlen = 0
    ifrom = iptr

    ! *** First character must be alphabetic

    if ((line(iptr:iptr) < 'a').or.(line(iptr:iptr) > 'z')) goto 1000
    iptr = iptr + 1
    if (iptr > linelen) goto 1000

    ! *** Now parse the rest of the letters (must be alphanumeric or _ )

10  continue
    if ( ((line(iptr:iptr) >= 'a').and.(line(iptr:iptr) <= 'z')).or. &
         ((line(iptr:iptr) == '_')).or. &
         ((line(iptr:iptr) >= '0').and.(line(iptr:iptr) <= '9')) ) then
       iptr = iptr + 1
       if (iptr <= linelen) goto 10
    end if

    ! *** Extract variable name

    ito = iptr - 1
    varlen = ito - ifrom + 1
    if (varlen > 0) varnam = line(ifrom:ito)

    ! *** Ignore intervening spaces

20  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= linelen) goto 20
    end if

    ! *** Now extract any subscript

    call get_subscript(isub1,isub2,icode)
    if (icode /= 0) then
       varlen = 0
       goto 1000
    end if

    ! *** Ignore intervening spaces

30  continue
    if (line(iptr:iptr) == ' ') then
       iptr = iptr + 1
       if (iptr <= linelen) goto 30
    end if

    ! *** We now expect '='

    if (line(iptr:iptr) == '=') then
       iptr = iptr + 1

       ! *** Restore original string's upper/lower case after '=' sign

       line(iptr:linelen) = line1(iptr:linelen)

    else
       varlen = 0
    end if

1000 continue

  end subroutine get_variable_name

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_range_int(cvar,varval,min_value,max_value)

    !! Routine that checks whether an integer variable lies within
    !! the desired range
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer  : Fortran output unit identifier
    !! cvar    : input string   : name of variable
    !! varval  : input integer  : value of variable
    !! min_value : input integer : minimum allowed value of variable
    !! max_value : input integer : maximum allowed value of variable
    !! This routine checks whether an integer variable lies within
    !! the range predetermined by the user, and reports an error
    !! and stops if it doesn't.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: cvar
    integer, intent(in) :: varval,min_value,max_value

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (min_value > max_value) then
       write(outfile,*) 'Illegal relative values of min_value and max_value'
       write(outfile,*) 'for variable ',cvar

       write(*,*) 'Illegal relative values of min_value and max_value'
       write(*,*) 'for variable ',cvar
          error = .True.
    end if

    if ((varval < min_value).or.(varval > max_value)) then
       write(outfile,*) cvar,' lies outside its allowed range :'
       write(outfile,*) 'Minimum value = ',min_value
       write(outfile,*) 'Maximum value = ',max_value
       write(outfile,*) ' Actual value = ',varval

       write(*,*) cvar,' lies outside its allowed range :'
       write(*,*) 'Minimum value = ',min_value
       write(*,*) 'Maximum value = ',max_value
       write(*,*) ' Actual value = ',varval
          error = .True.
    end if

  end subroutine check_range_int

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_range_real(cvar,varval,min_value,max_value)

    !! Routine that checks whether a real variable lies within
    !! the desired range
    !! author: P J Knight, CCFE, Culham Science Centre
    !! cvar    : input string   : name of variable
    !! varval  : input real     : value of variable
    !! min_value : input real   : minimum allowed value of variable
    !! max_value : input real   : maximum allowed value of variable
    !! This routine checks whether a real variable lies within
    !! the range predetermined by the user, and reports an error
    !! and stops if it doesn't.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    character(len=*), intent(in) :: cvar
    real(dp), intent(in) :: varval,min_value,max_value

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (min_value > max_value) then
       write(outfile,*) 'Illegal relative values of min_value and max_value'
       write(outfile,*) 'for variable ',cvar

       write(*,*) 'Illegal relative values of min_value and max_value'
       write(*,*) 'for variable ',cvar
          error = .True.
    end if

    if ((varval < min_value).or.(varval > max_value)) then
       write(outfile,*) cvar,' lies outside its allowed range :'
       write(outfile,*) 'Minimum value = ',min_value
       write(outfile,*) 'Maximum value = ',max_value
       write(outfile,*) ' Actual value = ',varval

       write(*,*) cvar,' lies outside its allowed range :'
       write(*,*) 'Minimum value = ',min_value
       write(*,*) 'Maximum value = ',max_value
       write(*,*) ' Actual value = ',varval
          error = .True.
    end if

  end subroutine check_range_real


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lower_case(string,start,finish)

   !! Routine that converts a (sub-)string to lowercase
   !! author: P J Knight, CCFE, Culham Science Centre
   !! string : input string   : character string of interest
   !! start  : optional input integer  : starting character for conversion
   !! finish : optional input integer  : final character for conversion
   !! This routine converts the specified section of a string
   !! to lowercase. By default, the whole string will be converted.
   !! None
   !
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   !  Arguments

   character(len=*), intent(inout) :: string
   integer, optional, intent(in) :: start,finish

   !  Local variables

   character(len=1) :: letter
   character(len=27), parameter :: lowtab = 'abcdefghijklmnopqrstuvwxyz_'
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
         i = index('ABCDEFGHIJKLMNOPQRSTUVWXYZ_',letter)
         if (i > 0) string(loop:loop) = lowtab(i:i)
      end do
   end if

 end subroutine lower_case

end module process_input
