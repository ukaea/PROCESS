! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial

    !! Routine to initialise
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !!     !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use define_iteration_variables, only: init_itv_1, init_itv_2, init_itv_3, &
        init_itv_4, init_itv_5, init_itv_6, init_itv_7, init_itv_8, init_itv_9, &
        init_itv_10, init_itv_11, init_itv_12, init_itv_13, init_itv_14, init_itv_15, &
        init_itv_16, init_itv_17, init_itv_18, init_itv_19, init_itv_20, init_itv_21, &
        init_itv_23, init_itv_25, init_itv_26, init_itv_27, init_itv_28, init_itv_29, &
        init_itv_30, init_itv_31, init_itv_32, init_itv_33, init_itv_34, init_itv_35, &
        init_itv_36, init_itv_37, init_itv_38, init_itv_39, init_itv_40, init_itv_41, &
        init_itv_42, init_itv_44, init_itv_45, init_itv_46, init_itv_47, init_itv_48, &
        init_itv_49, init_itv_50, init_itv_51, init_itv_53, init_itv_54, &
        init_itv_56, init_itv_57, init_itv_58, init_itv_59, init_itv_60, init_itv_61, &
        init_itv_62, init_itv_63, init_itv_64, init_itv_65, init_itv_66, init_itv_67, &
        init_itv_68, init_itv_69, init_itv_70, init_itv_71, init_itv_72, init_itv_73, &
        init_itv_74, init_itv_75, init_itv_79, init_itv_81, init_itv_82, init_itv_83, &
        init_itv_84, init_itv_85, init_itv_86, init_itv_89, init_itv_90, init_itv_91, &
        init_itv_92, init_itv_93, init_itv_94, init_itv_95, init_itv_96, init_itv_97, &
        init_itv_98, init_itv_103, init_itv_104, init_itv_105, &
        init_itv_106, init_itv_107, init_itv_108, init_itv_109, init_itv_110, &
        init_itv_111, init_itv_112, init_itv_113, init_itv_114, init_itv_115, &
        init_itv_116, init_itv_117, init_itv_118, init_itv_119, init_itv_120, &
        init_itv_121, init_itv_122, init_itv_123, init_itv_124, init_itv_125, &
        init_itv_126, init_itv_127, init_itv_128, init_itv_129, init_itv_130, &
        init_itv_131, init_itv_132, init_itv_133, init_itv_134, init_itv_135, &
        init_itv_136, init_itv_137, init_itv_138, init_itv_139, init_itv_140, &
        init_itv_141, init_itv_142, init_itv_143, init_itv_144, init_itv_145, &
        init_itv_146, init_itv_147, init_itv_148, init_itv_149, &
        init_itv_152, init_itv_153, init_itv_154, init_itv_155, &
        init_itv_156, init_itv_157, init_itv_158, init_itv_159, init_itv_160, &
        init_itv_161, init_itv_162, init_itv_163, init_itv_164, init_itv_165, &
        init_itv_166, init_itv_167, init_itv_168, init_itv_169, init_itv_170, &
        init_itv_171, init_itv_172, init_itv_173, init_itv_174, init_itv_175
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    !  Arguments

    !  Local variables

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! boundl(ipnvars) /../ : lower bounds on iteration variables
    !! boundu(ipnvars) /../ : upper bounds on iteration variables

    ! Issue #287  The initialization subroutines for the iteration variables are called
    call init_itv_1
    call init_itv_2
    call init_itv_3
    call init_itv_4
    call init_itv_5
    call init_itv_6
    call init_itv_7
    call init_itv_8
    call init_itv_9
    call init_itv_10
    call init_itv_11
    call init_itv_12
    call init_itv_13
    call init_itv_14
    call init_itv_15
    call init_itv_16
    call init_itv_17
    call init_itv_18
    call init_itv_19
    call init_itv_20
    call init_itv_21

    call init_itv_23

    call init_itv_25
    call init_itv_26
    call init_itv_27
    call init_itv_28
    call init_itv_29
    call init_itv_30
    call init_itv_31
    call init_itv_32
    call init_itv_33
    call init_itv_34
    call init_itv_35
    call init_itv_36
    call init_itv_37
    call init_itv_38
    call init_itv_39
    call init_itv_40
    call init_itv_41
    call init_itv_42

    call init_itv_44
    call init_itv_45
    call init_itv_46
    call init_itv_47
    call init_itv_48
    call init_itv_49
    call init_itv_50
    call init_itv_51

    call init_itv_53
    call init_itv_54

    call init_itv_56
    call init_itv_57
    call init_itv_58
    call init_itv_59
    call init_itv_60
    call init_itv_61
    call init_itv_62
    call init_itv_63
    call init_itv_64
    call init_itv_65
    call init_itv_66
    call init_itv_67
    call init_itv_68
    call init_itv_69
    call init_itv_70
    call init_itv_71
    call init_itv_72
    call init_itv_73
    call init_itv_74
    call init_itv_75




    call init_itv_79









    call init_itv_89
    call init_itv_90
    call init_itv_91
    call init_itv_92
    call init_itv_93
    call init_itv_94
    call init_itv_95
    call init_itv_96
    call init_itv_97
    call init_itv_98
    !Not used
    call init_itv_103
    call init_itv_104
    call init_itv_105
    call init_itv_106
    call init_itv_107
    call init_itv_108
    call init_itv_109
    call init_itv_110
    call init_itv_111
    call init_itv_112
    call init_itv_113
    call init_itv_114
    call init_itv_115
    call init_itv_116
    call init_itv_117
    call init_itv_118
    call init_itv_119
    call init_itv_120
    call init_itv_121
    call init_itv_122
    call init_itv_123
    call init_itv_124
    call init_itv_125
    call init_itv_126
    call init_itv_127
    call init_itv_128
    call init_itv_129
    call init_itv_130
    call init_itv_131
    call init_itv_132
    call init_itv_133
    call init_itv_134
    call init_itv_135
    call init_itv_136
    call init_itv_137
    call init_itv_138
    call init_itv_139
    call init_itv_140
    call init_itv_141
    call init_itv_142
    call init_itv_143
    call init_itv_144
    call init_itv_145
    call init_itv_146
    call init_itv_147
    call init_itv_148
    call init_itv_149
    call init_itv_152
    call init_itv_153
    call init_itv_154
    call init_itv_155
    call init_itv_156
    call init_itv_157
    call init_itv_158
    call init_itv_159
    call init_itv_160
    call init_itv_161
    call init_itv_162
    call init_itv_163
    call init_itv_164
    call init_itv_165
    call init_itv_166
    call init_itv_167
    call init_itv_168
    call init_itv_169
    call init_itv_170
    call init_itv_171
    call init_itv_172
    call init_itv_173
    call init_itv_174
    call init_itv_175


end subroutine initial

subroutine check

    !! Routine to reset specific variables if certain options are
    !! being used
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !! This routine performs a sanity check of the input variables
    !! and ensures other dependent variables are given suitable values.

    !!     !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use build_variables, only: blnkith, bore, gapoh, ohcth, precomp, iprecomp, &
        i_r_cp_top, r_cp_top, vgaptop, vgap, shldtth, shldlth, d_vv_top, d_vv_bot, tf_in_cs
    use buildings_variables, only: esbldgm3, triv
    use current_drive_variables, only: gamcd, iefrf, irfcd
    use error_handling, only: errors_on, idiags, fdiags, report_error
    use fwbs_variables, only: breeder_multiplier, iblanket, vfcblkt, vfpblkt, &
        iblnkith
    use global_variables, only: icase
    use heat_transport_variables, only: trithtmw
    use ife_variables, only: ife
    use impurity_radiation_module, only: nimp, impurity_arr_frac, fimp
    use numerics, only: ixc, icc, ioptimz, neqns, nineqns, nvar, boundl, &
        boundu
    use pfcoil_variables, only: ipfres, ngrp, pfclres, ipfloc, ncls, isumatoh
    use physics_variables, only: aspect, fdeut, fgwped, fhe3, &
        fgwsep, ftrit, ibss, i_single_null, icurr, idivrt, ishape, &
        iradloss, isc, ipedestal, ilhthresh, itart, nesep, rhopedn, rhopedt, &
        rnbeam, neped, te, tauee_in, tesep, teped, itartpf, ftar, idia
    use pulse_variables, only: lpulse
    use reinke_variables, only: fzactual, impvardiv
    use tfcoil_variables, only: casthi, casthi_is_fraction, casths, i_tf_sup, &
        tcoolin, tcpav, tfc_sidewall_is_fraction, tmargmin, tmargmin_cs, &
        tmargmin_tf, eff_tf_cryo, eyoung_ins, i_tf_bucking, i_tf_shape, &
        n_tf_graded_layers, n_tf_stress_layers, tlegav,  i_tf_stress_model, &
        i_tf_sc_mat, i_tf_wp_geom, i_tf_turns_integer, tinstf, thwcndut, &
        tfinsgap, rcool, dhecoil, thicndut, i_cp_joints, t_turn_tf_is_input, &
        t_turn_tf, tftmp, t_cable_tf, t_cable_tf_is_input, tftmp, tmpcry, &
        i_tf_cond_eyoung_axial, eyoung_cond_axial, eyoung_cond_trans, &
        i_tf_cond_eyoung_trans, i_str_wp
    use stellarator_variables, only: istell
    use sctfcoil_module, only: initialise_cables
    use vacuum_variables, only: vacuum_model
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none

    !  Local variables

    integer :: i,j,k,imp
    real(dp) :: fsum

    real(dp) :: dr_tf_wp_min
    !! Minimal WP or conductor layer thickness [m]
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    errors_on = .true.

    !  Check that there are sufficient iteration variables
    if (nvar < neqns) then
        idiags(1) = nvar ; idiags(2) = neqns
        call report_error(137)
    end if

    !  Check that sufficient elements of ixc and icc have been specified
    if ( any(ixc(1:nvar) == 0) ) then
        idiags(1) = nvar
        call report_error(139)
    end if


    if ( any(icc(1:neqns+nineqns) == 0) ) then
        idiags(1) = neqns ; idiags(2) = nineqns
        call report_error(140)
    end if

    !  Deprecate constraints 3 and 4
    if ( any(icc(1:neqns+nineqns) == 3) ) then
        call report_error(162)
        write(*,*) 'PROCESS stopping'
        stop 1
    end if

    if ( any(icc(1:neqns+nineqns) == 4) ) then
        call report_error(163)
        write(*,*) 'PROCESS stopping'
        stop 1
    end if


    ! MDK Report error if constraint 63 is used with old vacuum model
    if (any(icc(1:neqns+nineqns) == 63).and.(vacuum_model.ne.'simple') ) then
        write(*,*) 'Constraint 63 is requested without the correct vacuum model ("simple").'
        write(*,*) 'vacuum_model = ', vacuum_model
        write(*,*) 'PROCESS stopping'
        stop 1
    end if

    if ( any(icc(1:neqns+nineqns) == 74) ) then
        write(*,*)'Constraint 74 (TF coil quench temperature for Croco HTS conductor) is not yet implemented'
        write(*,*) 'PROCESS stopping'
        stop 1
    end if

    !  Fuel ion fractions must add up to 1.0
    if (abs(1.0D0 - fdeut - ftrit - fhe3) > 1.0D-6) then
        fdiags(1) = fdeut; fdiags(2) = ftrit ; fdiags(3) = fhe3
        call report_error(36)
    end if

    if (ftrit < 1.0D-3) then  !  tritium fraction is negligible
        triv = 0.0D0
        trithtmw = 0.0D0
    end if

    if (fimp(2) .ne. 0.1D0) then
        write(*,*)'The thermal alpha/electron density ratio should be controlled using ralpne (itv 109) and not fimp(2).'
        write(*,*)'fimp(2) should be removed from the input file, or set to the default value 0.1D0.'
        stop 1
    end if

    !  Impurity fractions
    do imp = 1,nimp
        impurity_arr_frac(imp) = fimp(imp)
    end do

    ! The 1/R B field dependency constraint variable is being depreciated
    ! Stop the run if the constraint 10 is used
    if ( any( icc == 10 ) ) then
        call report_error(236)
        stop 1
    end if

    ! Stop the run if oacdcp is used as an optimisation variable
    ! As the current density is now calculated from bt without constraint 10
    if ( any( ixc == 12 ) ) then
        call report_error(236)
        stop 1
    end if

    !  Warn if ion power balance equation is being used with the new radiation model
    if (any(icc == 3)) then
        call report_error(138)
    end if

    !  Plasma profile consistency checks
    if (ife /= 1) then
        if (ipedestal == 1) then

            !  Temperature checks
            if (teped < tesep) then
                fdiags(1) = teped ; fdiags(2) = tesep
                call report_error(146)
            end if

            if ((abs(rhopedt-1.0D0) <= 1.0D-7).and.((teped-tesep) >= 1.0D-7)) then
                fdiags(1) = rhopedt ; fdiags(2) = teped ; fdiags(3) = tesep
                call report_error(147)
            end if

            !  Core temperature should always be calculated (later) as being
            !  higher than the pedestal temperature, if and only if the
            !  volume-averaged temperature never drops below the pedestal
            !  temperature. Prevent this by adjusting te, and its lower bound
            !  (which will only have an effect if this is an optimisation run)
            if (te <= teped) then
                fdiags(1) = te ; fdiags(2) = teped
                te = teped*1.001D0
                call report_error(149)
            end if

            if ((ioptimz >= 0).and.(any(ixc == 4)).and.(boundl(4) < teped*1.001D0)) then
                call report_error(150)
                boundl(4) = teped*1.001D0
                boundu(4) = max(boundu(4), boundl(4))
            end if

             !  Density checks
             !  Case where pedestal density is set manually
             ! ---------------
             if ( (fgwped < 0) .or. (.not.any(ixc==145)) ) then

                 ! Issue #589 Pedestal density is set manually using neped but it is less than nesep.
                 if ( neped < nesep ) then
                     fdiags(1) = neped ; fdiags(2) = nesep
                     call report_error(151)
                 end if

                 ! Issue #589 Pedestal density is set manually using neped,
                 ! but pedestal width = 0.
                 if ( (abs(rhopedn-1.0D0) <= 1.0D-7).and.((neped-nesep) >= 1.0D-7) ) then
                     fdiags(1) = rhopedn ; fdiags(2) = neped ; fdiags(3) = nesep
                     call report_error(152)
                 end if
             end if

             ! Issue #862 : Variable ne0/neped ratio without constraint eq 81 (ne0>neped)
             !  -> Potential hollowed density profile
             if ( (ioptimz >= 0) .and. (.not.any(icc==81)) ) then
                 if ( any(ixc == 145 )) call report_error(154)
                 if ( any(ixc ==   6 )) call report_error(155)
             end if
         end if
     end if
     ! ---------------


     ! Cannot use Psep/R and PsepB/qAR limits at the same time
     if(any(icc == 68) .and. any(icc == 56)) then
        call report_error(178)
     endif

     if ((any(ixc==145)) .and. (boundl(145) < fgwsep)) then  !if lower bound of fgwped < fgwsep
        fdiags(1) = boundl(145); fdiags(2) = fgwsep
        call report_error(186)
     end if

     if (any(icc == 78)) then

        !If Reinke criterion is used tesep is calculated and cannot be an
        !iteration variable
        if (any(ixc == 119)) then
           call report_error(219)
        endif

        !If Reinke criterion is used need to enforce LH-threshold
        !using Martin scaling for consistency
        if (.not. ilhthresh == 6) then
           call report_error(218)
        endif
        if  (.not. any(icc==15) .and. (ipedestal .ne. 3)) then
           call report_error(218)
        endif


     endif

     if (any(icc == 78)) then

        !If Reinke criterion is used tesep is calculated and cannot be an
        !iteration variable
        if (any(ixc == 119)) then
           call report_error(219)
        endif

        !If Reinke criterion is used need to enforce LH-threshold
        !using Martin scaling for consistency
        if (.not. ilhthresh == 6) then
           call report_error(218)
        endif
        if  (.not. any(icc==15) .and. (ipedestal .ne. 3)) then
           call report_error(218)
        endif


     endif

     if (i_single_null == 0) then
         idivrt = 2
         vgaptop = vgap
         shldtth = shldlth
         d_vv_top = d_vv_bot
         call report_error(272)
     else  !  i_single_null == 1
         idivrt = 1
     end if


    !  Tight aspect ratio options (ST)
    ! --------------------------------
    if ( itart == 1 ) then

        icase  = 'Tight aspect ratio tokamak model'

        ! Disabled Forcing that no inboard breeding blanket is used
        ! Disabled iblnkith = 0

        ! Check if the choice of plasma current is addapted for ST
        ! 2 : Peng Ip scaling (See STAR code documentation)
        ! 9 : Fiesta Ip scaling
        if (icurr /= 2 .and. icurr /= 9) then
            idiags(1) = icurr ; call report_error(37)
        end if

        !! If using Peng and Strickler (1986) model (itartpf == 0)
        ! Overwrite the location of the TF coils
        ! 2 : PF coil on top of TF coil
        ! 3 : PF coil outside of TF coil
        if (itartpf == 0) then
          ipfloc(1) = 2
          ipfloc(2) = 3
          ipfloc(3) = 3
        end if

        ! Water cooled copper magnets initalisation / checks
        if ( i_tf_sup == 0 ) then
            ! Check if the initial centrepost coolant loop adapted to the magnet technology
            ! Ice cannot flow so tcoolin > 273.15 K
            if ( tcoolin < 273.15D0 ) call report_error(234)

            ! Temperature of the TF legs cannot be cooled down
            if ( abs(tlegav+1.0D0) > epsilon(tlegav) .and. tlegav < 273.15D0 ) call report_error(239)

            ! Check if conductor upper limit is properly set to 50 K or below
            if ( any(ixc == 20 ) .and. boundu(20) < 273.15D0 ) call report_error(241)

        ! Call a lvl 3 error if superconductor magnets are used
        else if ( i_tf_sup == 1 ) then
            call report_error(233)

        ! Aluminium magnets initalisation / checks
        ! Initialize the CP conductor temperature to cryogenic temperature for cryo-al magnets (20 K)
        else  if ( i_tf_sup == 2 ) then

            ! Call a lvl 3 error if the inlet coolant temperature is too large
            ! Motivation : ill-defined aluminium resistivity fit for T > 40-50 K
            if ( tcoolin > 40.0D0 ) call report_error(235)

            ! Check if the leg average temperature is low enough for the resisitivity fit
            if ( tlegav > 50.0D0 ) call report_error(238)

            ! Check if conductor upper limit is properly set to 50 K or below
            if ( any(ixc == 20 ) .and. boundu(20) > 50.0D0 ) call report_error(240)

            ! Otherwise intitialise the average conductor temperature at
            tcpav = tcoolin

        end if

        ! Check if the boostrap current selection is addapted to ST
        if (ibss  == 1) call report_error(38)

        ! Check if a single null divertor is used in double null machine
        if (i_single_null == 0 .and. (ftar == 1.0 .or. ftar == 0.0)) then
            call report_error(39)
        end if

        ! Set the TF coil shape to picture frame (if default value)
        if ( i_tf_shape == 0 ) i_tf_shape = 2

        ! Warning stating that the CP fast neutron fluence calculation
        ! is not addapted for cryoaluminium calculations yet
        if ( i_tf_sup == 2 .and. any( icc == 85 ) .and. itart == 1 ) then
            call report_error(260)
        end if

        ! Setting the CP joints default options :
        !  0 : No joints for superconducting magents (i_tf_sup = 1)
        !  1 : Sliding joints for resistive magnets (i_tf_sup = 0, 2)
        if ( i_cp_joints == -1 ) then
            if ( i_tf_sup == 1 ) then
                i_cp_joints = 0
            else
                i_cp_joints = 1
            end if
        end if

        ! Checking the CP TF top radius
        if ( ( abs(r_cp_top) > epsilon(r_cp_top) .or. any(ixc(1:nvar) == 174) ) &
            .and. i_r_cp_top /= 1 ) then
             call report_error(267)
        end if
    ! --------------------------------


    ! Conventionnal aspect ratios specific
    ! ------------------------------------
    else

        if (icurr == 2 .or. icurr == 9) call report_error(40)

        ! Set the TF coil shape to PROCESS D-shape (if default value)
        if ( i_tf_shape == 0 ) i_tf_shape = 1

        !  Check PF coil configurations
        j = 0 ; k = 0
        do i = 1, ngrp
            if ((ipfloc(i) /= 2).and.(ncls(i) /= 2)) then
                idiags(1) = i ; idiags(2) = ncls(i)
                call report_error(41)
            end if

            if (ipfloc(i) == 2) then
                j = j + 1
                k = k + ncls(i)
            end if
        end do

        if (k == 1) call report_error(42)
        if (k > 2) call report_error(43)
        if ((i_single_null == 1).and.(j < 2)) call report_error(44)

        ! Constraint 10 is dedicated to ST designs with demountable joints
        if ( any(icc(1:neqns+nineqns) == 10 ) ) call report_error(259)

    end if
    ! ------------------------------------

    !  Pulsed power plant model
    if (lpulse == 1) then
        icase = 'Pulsed tokamak model'
    else
        esbldgm3 = 0.0D0
    end if

    !  Ensure minimum cycle time constraint is turned off
    !  (not currently available, as routine thrmal has been commented out)
    if ( any(icc == 42) ) then
        call report_error(164)
    end if



    ! TF coil
    ! -------
    ! TF stress model not defined of r_tf_inboard = 0
    ! Unless i_tf_stress_model == 2
    ! -> If bore + gapoh + ohcth = 0 and fixed and stress constraint is used
    !    Generate a lvl 3 error proposing not to use any stress constraints
    if (       ( .not. ( any(ixc == 16 ) .or. any(ixc == 29 ) .or. any(ixc == 42 ) ) ) & ! No bore,gapoh, ohcth iteration
         .and. ( abs(bore + gapoh + ohcth + precomp) < epsilon(bore) )                 & ! bore + gapoh + ohcth = 0
         .and. ( any(icc == 31) .or. any(icc == 32) )                                  & ! Stress constraints (31 or 32) is used
         .and. ( i_tf_stress_model /= 2 ) ) then                                         ! TF stress model can't handle no bore

        call report_error(246)
        stop 1
    end if

    ! Make sure that plane stress model is not used for resistive magnets
    if ( i_tf_stress_model == 1 .and. i_tf_sup /= 1 ) call report_error(253)

    ! bucking cylinder default option setting
    !  - bucking (casing) for SC i_tf_bucking ( i_tf_bucking = 1 )
    !  - No bucking for copper magnets ( i_tf_bucking = 0 )
    !  - Bucking for aluminium magnets ( i_tf_bucking = 1 )
    if ( i_tf_bucking == -1 ) then
        if ( i_tf_sup == 0 ) then
            i_tf_bucking = 0
        else
            i_tf_bucking = 1
        end if
    end if

    ! Ensure that the TF isnt placed against the
    ! CS which is now outside it
    if ( i_tf_bucking >= 2 .and. tf_in_cs == 1 ) then
        call report_error(281)
    end if
    ! Ensure that no pre-compression structure
    ! is used for bucked and wedged design
    if ( i_tf_bucking >= 2 .and. iprecomp == 1 ) then
        call report_error(252)
    end if

    ! Number of stress calculation layers
    ! +1 to add in the inboard TF coil case on the plasma side, per Issue #1509
    n_tf_stress_layers = i_tf_bucking + n_tf_graded_layers + 1

    ! If TFC sidewall has not been set by user
    if ( casths < 0.1d-10 ) tfc_sidewall_is_fraction = .true.

    ! If inboard TF coil case plasma side thickness has not been set by user
    if( casthi < 0.1d-10 ) casthi_is_fraction = .true.

    ! Setting the default cryo-plants efficiencies
    !-!
    if ( abs(eff_tf_cryo + 1.0D0) < epsilon(eff_tf_cryo) ) then

        ! The ITER cyoplant efficiency is used for SC
        if ( i_tf_sup == 1 ) then
            eff_tf_cryo = 0.13D0

        ! Strawbrige plot extrapolation is used for Cryo-Al
        else if ( i_tf_sup == 2 ) then
            eff_tf_cryo = 0.40D0
        end if

    ! Cryo-plane efficiency must be in [0-1.0]
    else if ( eff_tf_cryo >  1.0D0 .or. eff_tf_cryo < 0.0D0 ) then
        call report_error(248)
        stop 1
    end if
    !-!

    ! Integer turns option not yet available for REBCO taped turns
    !-!
    if ( i_tf_sc_mat == 6 .and. i_tf_turns_integer == 1 ) then
        call report_error(254)
        stop 1
    end if
    !-!


    ! Setting up insulation layer young modulae default values [Pa]
    !-!
    if ( abs(eyoung_ins - 1.0D8 ) < epsilon(eyoung_ins) ) then

        ! Copper magnets, no insulation material defined
        ! But use the ITER design by default
        if ( i_tf_sup == 0 ) then
            eyoung_ins = 20.0D9

        ! SC magnets
        ! Value from DDD11-2 v2 2 (2009)
        else if ( i_tf_sup == 1 ) then
            eyoung_ins = 20.0D9

        ! Cryo-aluminum magnets (Kapton polymer)
        else if ( i_tf_sup == 2 ) then
            eyoung_ins = 2.5D9
        end if
    end if
    !-!

    !-! Setting the default WP geometry
    !-!
    if ( i_tf_wp_geom == -1 ) then
        if ( i_tf_turns_integer == 0 ) i_tf_wp_geom = 1
        if ( i_tf_turns_integer == 1 ) i_tf_wp_geom = 0
    end if
    !-!

    !-! Setting the TF coil conductor elastic properties
    !-!
    if ( i_tf_cond_eyoung_axial == 0 ) then
        ! Conductor stiffness is not considered
        eyoung_cond_axial = 0
        eyoung_cond_trans = 0
    else if ( i_tf_cond_eyoung_axial == 2 ) then
        ! Select sensible defaults from the literature
        select case (i_tf_sc_mat)
            case (1,4,5)
                ! Nb3Sn: Nyilas, A et. al, Superconductor Science and Technology 16, no. 9 (2003): 1036–42. https://doi.org/10.1088/0953-2048/16/9/313.
                eyoung_cond_axial = 32D9
            case (2)
                ! Bi-2212: Brown, M. et al, IOP Conference Series: Materials Science and Engineering 279 (2017): 012022. https://doi.org/10.1088/1757-899X/279/1/012022.
                eyoung_cond_axial = 80D9
            case (3,7)
                ! NbTi: Vedrine, P. et. al, IEEE Transactions on Applied Superconductivity 9, no. 2 (1999): 236–39. https://doi.org/10.1109/77.783280.
                eyoung_cond_axial = 6.8D9
            case (6,8,9)
                ! REBCO: Fujishiro, H. et. al, Physica C: Superconductivity, 426–431 (2005): 699–704. https://doi.org/10.1016/j.physc.2005.01.045.
                eyoung_cond_axial = 145D9
        end select

        if ( i_tf_cond_eyoung_trans == 0) then
            ! Transverse stiffness is not considered
            eyoung_cond_trans = 0
        else
            ! Transverse stiffness is significant
            eyoung_cond_trans = eyoung_cond_axial
        end if
    end if
    !-!

    ! Check if the WP/conductor radial thickness (dr_tf_wp) is large enough
    ! To contains the insulation, cooling and the support structure
    ! Rem : Only verified if the WP thickness is used
    if ( any(ixc(1:nvar) == 140) ) then

        ! Minimal WP thickness
        if ( i_tf_sup == 1 ) then
            dr_tf_wp_min = 2.0D0 * ( tinstf + tfinsgap + thicndut + dhecoil )

            ! Steel conduit thickness (can be an iteration variable)
            if ( any(ixc(1:nvar) == 58 ) ) then
                dr_tf_wp_min = dr_tf_wp_min + 2.0D0 * boundl(58)
            else
                dr_tf_wp_min = dr_tf_wp_min + 2.0D0 * thwcndut
            end if

        ! Minimal conductor layer thickness
        else if ( i_tf_sup == 0 .or. i_tf_sup == 2 ) then
            dr_tf_wp_min = 2.0D0 * ( thicndut + tinstf ) + 4.0D0 * rcool
        end if

        if ( boundl(140) < dr_tf_wp_min ) then
            fdiags(1) = dr_tf_wp_min
            call report_error(255)
        end if
    end if

    ! Setting t_turn_tf_is_input to true if t_turn_tf is an input
    if ( abs(t_turn_tf) < epsilon(t_turn_tf) ) then
        t_turn_tf_is_input = .false.
    else
        t_turn_tf_is_input = .true.
    end if

    ! Impossible to set the turn size of integer turn option
    if ( t_turn_tf_is_input .and. i_tf_turns_integer == 1 ) then
        call report_error(269)
    end if

    if ( i_tf_wp_geom /= 0  .and. i_tf_turns_integer == 1 ) then
        call report_error(283)
    end if

    if ( ibss == 5  .and. idia /= 0 ) then
        call report_error(284)
    end if

    ! Setting t_cable_tf_is_input to true if t_cable_tf is an input
    if ( abs(t_cable_tf) < epsilon(t_cable_tf) ) then
        t_cable_tf_is_input = .false.
    else
        t_cable_tf_is_input = .true.
    end if

    ! Impossible to set the cable size of integer turn option
    if ( t_cable_tf_is_input .and. i_tf_turns_integer == 1 ) then
        call report_error(269)
    end if

    ! Impossible to set both the TF coil turn and the cable dimension
    if ( t_turn_tf_is_input .and. t_cable_tf_is_input ) then
        call report_error(271)
    end if

    ! Checking the SC temperature for LTS
    if ( ( i_tf_sc_mat == 1 .or. &
           i_tf_sc_mat == 3 .or. &
           i_tf_sc_mat == 4 .or. &
           i_tf_sc_mat == 5 ) .and. tftmp > 10.0D0 ) then
        call report_error(270)
    end if
    ! -------



    !  PF coil resistivity is zero if superconducting
    if (ipfres == 0) pfclres = 0.0D0

    !  If there is no NBI, then hot beam density should be zero
    if (irfcd == 1) then
        if ((iefrf /= 5).and.(iefrf /= 8)) rnbeam = 0.0D0
    else
        rnbeam = 0.0D0
    end if

    ! Set inboard blanket thickness to zero if no inboard blanket switch
    ! used (Issue #732)
    if (iblnkith == 0) blnkith = 0.0D0

    !  Solid breeder assumed if ipowerflow=0

    !if (ipowerflow == 0) blkttype = 3

    !  Set coolant fluid type

    !if ((blkttype == 1).or.(blkttype == 2)) then
    !   coolwh = 2  !  water
    !else
    !   coolwh = 1  !  helium
    !end if

    !  But... set coolant to water if blktmodel > 0
    !  Although the *blanket* is by definition helium-cooled in this case,
    !  the shield etc. are assumed to be water-cooled, and since water is
    !  heavier (and the unit cost of pumping it is higher), the calculation
    !  for coolmass is better done with coolwh=2 if blktmodel > 0 to give
    !  slightly pessimistic results.

    !if (blktmodel > 0) then
    !   secondary_cycle = 0
    !   blkttype = 3  !  HCPB
    !   coolwh = 2
    !end if

    !  Ensure that blanket material fractions allow non-zero space for steel
    !  CCFE HCPB Model

    if (istell == 0) then
        if ((iblanket == 1).or.(iblanket == 3)) then
            fsum = breeder_multiplier + vfcblkt + vfpblkt
            if (fsum >= 1.0D0) then
                idiags(1) = iblanket
                fdiags(2) = breeder_multiplier
                fdiags(3) = vfcblkt
                fdiags(4) = vfpblkt
                fdiags(5) = fsum
                call report_error(165)
            end if
        end if
    end if

    ! Initialise superconductor cable parameters
    if(i_tf_sup==1)then
        call initialise_cables()
    end if

    ! Check that the temperature margins are not overdetermined
    if(tmargmin>0.0001d0)then
        ! This limit has been input and will be applied to both TFC and CS
        if(tmargmin_tf>0.0001d0)then
            write(*,*)'tmargmin_tf and tmargmin should not both be specified in IN.DAT.'
            write(*,*)'tmargmin_tf has been ignored.'
        end if
        if(tmargmin_cs>0.0001d0)then
            write(*,*)'tmargmin_cs and tmargmin should not both be specified in IN.DAT.'
            write(*,*)'tmargmin_cs has been ignored.'
        end if
        tmargmin_tf = tmargmin
        tmargmin_cs = tmargmin
     end if

     if (tauee_in.ge.1.0D-10.and.isc.ne.48) then
        ! Report error if confinement time is in the input
        ! but the scaling to use it is not selected.
        call report_error(220)
     end if

     if (aspect.gt.1.7D0.and.isc.eq.46) then
        ! NSTX scaling is for A<1.7
        call report_error(221)
     end if

    if (icurr.eq.2.and.isc.eq.42) then
        call report_error(222)
    end if

    ! Cannot use temperature margin constraint with REBCO TF coils
    if(any(icc == 36) .and. ((i_tf_sc_mat == 8).or.(i_tf_sc_mat == 9))) then
        call report_error(265)
    endif

    ! Cannot use temperature margin constraint with REBCO CS coils
    if(any(icc == 60) .and. (isumatoh == 8)) then
        call report_error(264)
    endif

    ! Cold end of the cryocooler should be colder than the TF
    if(tmpcry > tftmp) then
        call report_error(273)
    endif

    ! Cannot use TF coil strain limit if i_str_wp is off:
    if(any(icc == 88) .and. (i_str_wp == 0)) then
        call report_error(275)
    endif

    errors_on = .false.

    ! Disable error logging only after all checks have been performed.
    ! (CPSS #1582: Why is error logging disabled at all?)
    errors_on = .false.


end subroutine check
