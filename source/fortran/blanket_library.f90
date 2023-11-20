module blanket_library

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! This library contains routines that can be shared by the blanket modules used in PROCESS.

    !!! These include:
    !!! - component_volumes
    !!! - component_masses
    !!! - thermo_hydraulic_model

    !!! author: G Graham, CCFE, Culham Science Centre

    !! Acronyms for this module:
    !!
    !!      BB          Breeding Blanket
    !!      FW          First Wall
    !!      BZ          Breeder Zone
    !!      MF/BSS      Manifold/Back Supporting Structure
    !!      LT          Low Temperature
    !!      HT          High Temperature
    !!      MMS         Multi Module Segment
    !!      SMS         Single Modle Segment
    !!      IB          Inboard
    !!      OB          Outboard
    !!      HCD         Heating & Current Drive
    !!      FCI         Flow Channel Insert

    !!! Any changes within a subroutine or function code will have a comment explaining the change

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
     use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    real(dp) :: volshldi, volshldo
    !! Volume of inboard and outboard shield (m3)

    real(dp) :: volvvi, volvvo
    !! Volume of inboard and outboard Vacuum Vessel (m3)

    real(dp) :: hcryopf
    !! Clearance between uppermost PF coil and cryostat lid (m)

    real(dp) :: vfblkti, vfblkto
    !! Inboard/outboard void fraction of blanket

    real(dp) :: bldepti, bldepto
    !! Inboard/outboard blanket coolant channel length (radial direction) (m)

    real(dp) :: blwidti, blwidto
    !! Inboard/outboard blanket mid-plan toroidal circumference for segment (m)

    real(dp) :: bllengi, bllengo
    !! Inboard/outboard blanket segment poloidal length (m)

    real(dp) :: bzfllengi, bzfllengo
    !! Inboard/outboard primary blanket flow lengths (m)

    real(dp) :: bzfllengi_liq, bzfllengo_liq
    !! Inboard/outboard secondary blanket flow lengths (m)

    real(dp) :: pnucfwi, pnucfwo
    !! Inboard/outboard first wall nuclear heating (MW)

    real(dp) :: tpeakfwi, tpeakfwo
    !! Inboard/outboard first wall peak temperature (K)

    real(dp) :: mffwi, mffwo, mffw
    !! Inboard/outboard total mass flow rate to remove inboard FW power (kg/s)

    real(dp) :: npfwi, npfwo
    !! Inboard/utboard total number of pipes

    real(dp) :: mffwpi, mffwpo
    !! Inboard/outboard mass flow rate per coolant pipe (kg/s)

    real(dp) :: pnucblkti, pnucblkto
    !! Neutron power deposited inboard/outboard blanket blanket (MW)

    real(dp) :: mfblkti, mfblkto, mfblkt
    !! Inboard/outboard blanket mass flow rate for coolant (kg/s)

    real(dp):: mfblkti_liq, mfblkto_liq, mfblkt_liq
    !! Inboard/outboard blanket mass flow rate for liquid breeder (kg/s)

    real(dp) :: mftotal
    !! Total mass flow rate for coolant (kg/s)

    real(dp) :: npblkti, npblkto
    !! Inboard/outboard total num of pipes

    real(dp) :: mfblktpi, mfblktpo
    !! Inboard/outboard mass flow rate per coolant pipe (kg/s)

    real(dp) :: velblkti, velblkto
    !! Inboard/outboard coolant velocity in blanket (m/s)

    real(dp) :: htpmw_fwi, htpmw_fwo
    !! Inboard/outboard first wall pumping power (MW)

    real(dp) :: htpmw_blkti, htpmw_blkto
    !! Inboard/outboard blanket pumping power (MW)

    real(dp) :: htpmw_fw_blkti, htpmw_fw_blkto
    !! Inboard/outboard fw and blanket pumping power (MW)

    real(dp) :: hblnkt
    !! Blanket internal half-height (m)

    real(dp) :: hshld
    !! Shield internal half-height (m)

    real(dp) :: hvv
    !! Vacuum vessel internal half-height (m)

    integer :: icomponent
    !! Switch used to specify selected component: blanket=0, shield=1, vacuum vessel=2

contains

    subroutine init_blanket_library
        !! Initialise module variables

        implicit none

        hblnkt = 0.0D0
        hshld = 0.0D0
        hcryopf = 0.0D0
        hvv = 0.0D0
        volshldi = 0.0D0
        volshldo = 0.0D0
        volvvi = 0.0D0
        volvvo = 0.0D0
        bldepti = 0.0D0
        bldepto = 0.0D0
        blwidti = 0.0D0
        blwidto = 0.0D0
        bllengi = 0.0D0
        bllengo = 0.0D0
        bzfllengi = 0.0D0
        bzfllengi_liq = 0.0D0
        bzfllengo_liq = 0.0D0
        bzfllengo = 0.0D0
        pnucfwi = 0.0D0
        pnucfwo = 0.0D0
        tpeakfwi = 0.0D0
        tpeakfwo = 0.0D0
        mffwi = 0.0D0
        mffwo = 0.0D0
        mffw = 0.0D0
        npfwi = 0.0D0
        npfwo = 0.0D0
        mffwpi = 0.0D0
        mffwpo = 0.0D0
        pnucblkti = 0.0D0
        pnucblkto = 0.0D0
        mfblkti = 0.0D0
        mfblkto = 0.0D0
        mfblkti_liq = 0.0D0
        mfblkto_liq = 0.0D0
        mfblkt_liq = 0.0D0
        mfblkt = 0.0D0
        mftotal = 0.0D0
        npblkti = 0.0D0
        npblkto = 0.0D0
        mfblktpi = 0.0D0
        mfblktpo = 0.0D0
        velblkti = 0.0D0
        velblkto = 0.0D0
        htpmw_fwi = 0.0D0
        htpmw_fwo = 0.0D0
        htpmw_blkti = 0.0D0
        htpmw_blkto = 0.0D0
        vfblkti = 0.0D0
        vfblkto = 0.0D0

    end subroutine init_blanket_library

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! VOLUME CALCULATIONS
    !!! component_volumes:
    !!!     - component_half_height
    !!!     - dshaped_component
    !!!     - elliptical_component
    !!!     - apply_coverage_factors
    !!!     - external_cryo_geometry
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine component_volumes
        !! Calculate the blanket, shield, vacuum vessel and cryostat volumes
        !! author: J. Morris, CCFE, Culham Science Centre
        !! Calculate the blanket, shield, vacuum vessel and cryostat volumes
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: fwbsshape
        use physics_variables, only: itart

        implicit none

        ! N.B. icomponent is a switch used to specify selected component: blanket=0, sheild=1, vacuum vessel=2
        ! Replaced seperate subroutines for blnkt, shld and vv with fuction/subroutine with icomponent switch.

        ! Calculate half-height
        ! Blanket
        hblnkt = component_half_height(icomponent=0)
        ! Shield
        hshld = component_half_height(icomponent=1)
        ! Vacuum Vessel
        hvv = component_half_height(icomponent=2)

        ! D-shaped blanket and shield
        if ((itart == 1).or.(fwbsshape == 1)) then

            do icomponent = 0,2
                call dshaped_component(icomponent)
            enddo

        ! Elliptical blanket and shield
        else

             do icomponent = 0,2
                 call elliptical_component(icomponent)
             enddo

            ! This will fail the hts_REBCO and 2D_scan regression tests,
            ! the number of VMCON iterations (nviter) is different.
            ! Seems to be because in the blanket calculations (icomponent=0):
            ! r2 = 1.3836567143743970 rather than old value of r2 = 1.3836567143743972,
            ! r3 = 3.7009701431231936 rather than r3 = 3.7009701431231923.

        end if

        ! Apply coverage factors to volumes and surface areas
        call apply_coverage_factors

        ! Calculate cryostat geometry
        call external_cryo_geometry

    end subroutine component_volumes

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function component_half_height(icomponent)
        !! Calculate the blanket, shield or vacuum vessel half-height
        !! Based on blanket_half_height, shield_half_height, vv_half_height
        !! original author: J. Morris, CCFE, Culham Science Centre
        !! author: G. Graham, CCFE, Culham Science Centre
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use build_variables, only: hmax, vgap, vgap2, blnktth, shldtth, scrapli, scraplo, &
        fwith, fwoth, d_vv_bot, d_vv_top
        use physics_variables, only: rminor, kappa, idivrt
        use divertor_variables, only: divfix

        implicit none

        ! Input variables
        integer, intent(in) :: icomponent
        ! icomponent - blnkt=0, shld=1, vv=2

        ! Local variables
        real(dp) :: hbot, htop

        ! Return variable
        real(dp) :: component_half_height

        ! Calculate component internal lower half-height (m)
        ! Blanket
        if (icomponent==0) hbot = rminor*kappa + vgap + divfix - blnktth
        ! Sheild
        if (icomponent==1) hbot = rminor*kappa + vgap + divfix
        ! Vacuum vessel
        if (icomponent==2) hbot = hmax - vgap2 - d_vv_bot

        ! Calculate component internal upper half-height (m)
         ! If a double null machine then symmetric
        if (idivrt == 2) then
            htop = hbot
        else
            ! Blanket
            htop = rminor*kappa + 0.5D0*(scrapli+scraplo + fwith+fwoth)
            ! Shield
            if (icomponent==1) htop = htop + blnktth
            ! Vacuum Vessel
            if (icomponent==2) htop = htop + blnktth + shldtth
        end if

        ! Average of top and bottom (m)
        component_half_height = 0.5D0*(htop + hbot)

    end function component_half_height

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dshaped_component(icomponent)

        !! Calculate component surface area and volume using dshaped scheme
        !! Based on dshaped_blanket, dshaped_shield, dshaped_vv
        !! original author: J. Morris, CCFE, Culham Science Centre
        !! author: G. Graham, CCFE, Culham Science Centre
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use maths_library, only: dshellarea, dshellvol
        use build_variables, only: rsldi, shldith, blnkith, fwith, scrapli, scraplo, &
        fwoth, blareaib, blareaob, blarea, blnkoth, blnktth, &
        shareaib, shareaob, sharea, shldoth, shldtth, &
        rsldo, d_vv_in, d_vv_out, d_vv_top, d_vv_bot
        use fwbs_variables, only: volblkti, volblkto, volblkt, volshld, vdewin
        use physics_variables, only: rminor

        implicit none

        ! Input variables
        integer, intent(in) :: icomponent
        ! icomponent - blnkt=0, shld=1, vv=2

        ! Local variables
        real(dp) :: r1, r2

        ! Calculate major radius to outer edge of inboard ...
        ! ... section (m)
        r1 = rsldi
        ! ... shield (m)
        if (icomponent==1) r1 = r1 + shldith
        ! ... blanket (m)
        if (icomponent==0) r1 = r1 + shldith + blnkith

        ! Horizontal distance between inside edges (m)
        ! i.e. outer radius of inboard part to inner radius of outboard part
        ! Blanket
        r2 = fwith + scrapli + 2.0D0*rminor + scraplo + fwoth
        ! Sheild
        if (icomponent==1) r2 = blnkith + r2 + blnkoth
        ! Vaccum Vessel
        if (icomponent==2) r2 = rsldo - r1

        ! Calculate surface area, assuming 100% coverage
        if (icomponent==0) call dshellarea(r1, r2, hblnkt, blareaib, blareaob, blarea)
        if (icomponent==1) call dshellarea(r1, r2, hshld, shareaib, shareaob, sharea)

        ! Calculate volumes, assuming 100% coverage
        if (icomponent==0) call dshellvol(r1, r2, hblnkt, blnkith, blnkoth, blnktth, volblkti, volblkto, volblkt)
        if (icomponent==1) call dshellvol(r1, r2, hshld, shldith, shldoth, shldtth, volshldi, volshldo, volshld)
        if (icomponent==2) call dshellvol(r1, r2, hvv, d_vv_in, d_vv_out, &
                    (d_vv_top+d_vv_bot)/2, volvvi, volvvo, vdewin)

    end subroutine dshaped_component

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine elliptical_component(icomponent)

        !! Calculate component surface area and volume using elliptical scheme
        !! Based on elliptical_blanket, elliptical_shield, elliptical_vv
        !! original author: J. Morris, CCFE, Culham Science Centre
        !! author: G. Graham, CCFE, Culham Science Centre
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use maths_library, only: eshellarea, eshellvol
        use build_variables, only: rsldi, shldith, blnkith, rsldo, shldoth, blnkoth, &
        blareaib, blareaob, blarea, blnktth, &
        shareaib, shareaob, sharea, shldtth, &
        d_vv_in, d_vv_out, d_vv_top, d_vv_bot
        use fwbs_variables, only: volblkti, volblkto, volblkt, volshld, vdewin
        use physics_variables, only: rmajor, rminor, triang

        implicit none

        ! Input variables
        integer, intent(in) :: icomponent
        ! icomponent - blnkt=0, shld=1, vv=2

        ! Local variables
        real(dp) :: r1, r2, r3

        ! Major radius to centre of inboard and outboard ellipses (m)
        ! (coincident in radius with top of plasma)
        r1 = rmajor - rminor*triang

        ! Calculate distance between r1 and outer edge of inboard ...
        ! ... section (m)
        r2 = r1 - rsldi
        ! ... shield (m)
        if (icomponent==1) r2 = r2 - shldith
        ! ... blanket (m)
        if (icomponent==0) r2 = r2 - shldith - blnkith

        ! Calculate distance between r1 and inner edge of outboard ...
        ! ... section (m)
        r3 = rsldo - r1
        ! ... shield (m)
        if (icomponent==1) r3 = r3 - shldoth
        ! ... blanket (m)
        if (icomponent==0) r3 = r3 - shldoth - blnkoth

        ! Calculate surface area, assuming 100% coverage
        if (icomponent==0) call eshellarea(r1, r2, r3, hblnkt, blareaib, blareaob, blarea)
        if (icomponent==1) call eshellarea(r1, r2, r3, hshld, shareaib, shareaob, sharea)

        ! Calculate volumes, assuming 100% coverage
        if (icomponent==0) call eshellvol(r1, r2, r3, hblnkt, blnkith, blnkoth, blnktth, &
                    volblkti, volblkto, volblkt)
        if (icomponent==1) call eshellvol(r1, r2, r3, hshld, shldith, shldoth, shldtth, &
                    volshldi, volshldo, volshld)
        if (icomponent==2) call eshellvol(r1, r2, r3, hvv, d_vv_in, d_vv_out, &
                    (d_vv_top+d_vv_bot)/2, volvvi, volvvo, vdewin)

    end subroutine elliptical_component

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine apply_coverage_factors
        !! Apply coverage factors to volumes
        !! author: J. Morris, CCFE, Culham Science Centre
        !! Apply coverage factors to volumes
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use build_variables, only: blareaob, blarea, blareaib, shareaib, shareaob, &
          sharea
        use fwbs_variables, only: fdiv, fhcd, volblkto, volblkti, volblkt, fvolsi, &
          fvolso, volshld, vdewin, fvoldw
        use physics_variables, only: idivrt

        implicit none

        ! Apply blanket coverage factors
        if (idivrt == 2) then
          ! double null configuration
          blareaob = blarea*(1.0D0-2.0D0*fdiv-fhcd) - blareaib
        else
          ! single null configuration
          blareaob = blarea*(1.0D0-fdiv-fhcd) - blareaib
        end if

        blarea = blareaib + blareaob

        volblkto = volblkt*(1.0D0-fdiv-fhcd) - volblkti
        volblkt = volblkti + volblkto

        ! Apply shield coverage factors
        shareaib = fvolsi*shareaib
        shareaob = fvolso*shareaob
        sharea = shareaib + shareaob

        volshldi = fvolsi*volshldi
        volshldo = fvolso*volshldo
        volshld = volshldi + volshldo

        ! Apply vacuum vessel coverage factor
        ! moved from dshaped_* and elliptical_* to keep coverage factor
        ! changes in the same location.
        vdewin = fvoldw*vdewin

    end subroutine apply_coverage_factors

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine external_cryo_geometry
        !! Calculate cryostat geometry
        !! author: J. Morris, CCFE, Culham Science Centre
        !! Calculate cryostat geometry
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use constants, only: pi
        use build_variables, only: clhsf, hmax, tfcth, ddwex
        use fwbs_variables, only: rdewex, rpf2dewar, zdewex, vdewex, vvmass, vdewin, &
          denstl, dewmkg
        use pfcoil_variables, only: rb, zh
        use buildings_variables, only: clh1

        implicit none

        ! cryostat radius (m)
        ! ISSUE #508 Remove RFP option
        ! rb(i) = outer radius of PF coil i (tokamaks)
        rdewex = maxval(rb) + rpf2dewar

        ! Clearance between uppermost PF coil and cryostat lid (m).
        ! Scaling from ITER by M. Kovari
        hcryopf = clhsf * (2.0D0*rdewex)/28.440D0

        ! Half-height of cryostat (m)
        ! ISSUE #508 Remove RFP option
        zdewex = maxval(zh) + hcryopf

        ! Vertical clearance between TF coil and cryostat (m)
        clh1 = zdewex - (hmax + tfcth)

        ! cryostat volume (m3)
        vdewex = ( (2.0D0*pi*rdewex) * 2.0D0*zdewex + (2.0D0*pi*rdewex**2) ) * ddwex

        ! Vacuum vessel mass (kg)
        vvmass = vdewin * denstl

        ! Sum of internal vacuum vessel and cryostat masses (kg)
        dewmkg = (vdewin + vdewex) * denstl

    end subroutine external_cryo_geometry

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! MASS CALCULATIONS
    !!! component_masses
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine component_masses

        !! Calculations for component masses
        !! author: J. Morris, CCFE, Culham Science Centre
        !! Calculations for component masses
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! variables used in ccfe hcpb but not kit hcpb:
        use constants, only: pi
        use build_variables, only: fwareaib, fwith, fwareaob, fwoth, fwarea
        use divertor_variables, only: divsur, divclfr, divplt, fdiva, divmas, divdens
        use physics_variables, only: rminor, rmajor, idivrt, sarea

        ! variables below used in kit hcpb but not ccfe hcpb:
        use build_variables, only: blnkith, blbuith, blbmith, blbpith, blnkoth, &
            blbuoth, blbmoth, blbpoth

        ! fwbs_variables:
        use fwbs_variables, only: volblkt, vfblkt, & !! CCFE and KIT HCPB --------
        whtblbe, whtblss, denstl, whtblkt, &
        volshld, vfshld, coolmass, fwclfr, & !! CCFE HCPB only -------------------
        breeder_f, breeder_multiplier, whtbltibe12, whtblli4sio4, wtblli2o, &
        vfcblkt, vfpblkt, whtshld, wpenshld, fwmass, fw_armour_vol, &
        fw_armour_thickness, fw_armour_mass, armour_fw_bl_mass, &
        volblkti, volblkto, iblnkith, fblhebmi, & !! KIT HCPB only ---------------
        fblhebpi,fblhebmo, fblhebpo, fblss, fblbe, &
        whtblbreed, densbreed, fblbreed, &
        iblanket, denw, vffwi, vffwo, volfw, & !! added --------------------
        fblss_ccfe, fblli2sio4, fbltibe12

        implicit none

        ! Only ccfe hcpb has local varibles i.e. coolvol
        ! Local variables

        ! Coolant volume (m3)
        real(dp) :: coolvol

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! CCFE HCPB modal calculates the coolant mass,
        ! have added an if staement using the iblanket switch for this.
        ! N.B. iblanket=1 for CCFE HCPB and iblanket=3 for the same with TBR using Shimwell.

        if ((iblanket==1).or.(iblanket==3)) then

            ! Start adding components of the coolant mass:
            ! Divertor coolant volume (m3)
            coolvol = divsur * divclfr * divplt

            ! Blanket coolant volume (m3)
            coolvol = coolvol + volblkt*vfblkt

            ! Shield coolant volume (m3)
            coolvol = coolvol + volshld*vfshld

            ! First wall coolant volume (m3)
            coolvol = coolvol + fwareaib*fwith*vffwi + fwareaob*fwoth*vffwo

            ! Mass of He coolant = volume * density at typical coolant temperatures and pressures (kg)
            coolmass = coolvol*1.517D0

            ! Average first wall coolant fraction, only used by old routines in fispact.f90, safety.f90
            fwclfr = (fwareaib*fwith*vffwi + fwareaob*fwoth*vffwo) / (fwarea*0.5D0*(fwith+fwoth))

        endif

        ! CCFE HCPB calculates the mass of the divertor, blanket (including seprate masses for each material),
        ! shield, FW and FW armour.
        ! KIT HCPB calculates the mass of the blanket (including seprate masses for each material)
        ! and the void fraction for the blanket.
        ! N.B. iblanket=1 for CCFE HCPB and iblanket=3 for the same with TBR using Shimwell.

        if ((iblanket==1).or.(iblanket==3)) then

            ! Component masses

            ! Divertor mass (kg)
            divsur = fdiva * 2.0D0 * pi * rmajor * rminor
            if (idivrt == 2) divsur = divsur * 2.0D0
            divmas = divsur * divdens * (1.0D0 - divclfr) * divplt

            ! Shield mass (kg)
            whtshld = volshld * denstl * (1.0D0 - vfshld)

            ! Penetration shield mass (set = internal shield) (kg)
            wpenshld = whtshld

            ! First wall volume (m^3)
            volfw = (fwareaib*fwith*(1.0D0-vffwi) + fwareaob*fwoth*(1.0D0-vffwo))

            ! First wall mass, excluding armour (kg)
            fwmass = denstl * volfw

            ! First wall armour volume (m^3)
            fw_armour_vol = sarea*fw_armour_thickness

            ! First wall armour mass (kg)
            fw_armour_mass = fw_armour_vol*denw

        endif

        if ((iblanket==1).or.(iblanket==3)) then

            if (breeder_f < 1.0D-10) breeder_f = 1.0D-10
            if (breeder_f > 1.0D0  ) breeder_f = 1.0D0

            ! fbltibe12 = fblli2sio4 * (1 - breeder_f)/breeder_f
            ! New combined variable breeder_multiplier
            ! Lithium orthosilicate fraction:
            fblli2sio4 = breeder_f * breeder_multiplier

            ! Titanium beryllide fraction, and mass (kg):
            fbltibe12  = breeder_multiplier - fblli2sio4
            whtbltibe12 = volblkt * fbltibe12 * 2260.0D0

            ! Blanket Lithium orthosilicate mass (kg)
            ! Ref: www.rockwoodlithium.com...
            whtblli4sio4 = volblkt * fblli2sio4 * 2400.0D0

            ! TODO sort this out so that costs model uses new variables.
            ! #327 For backwards compatibility, set the old blanket masses the same:
            whtblbe = whtbltibe12
            wtblli2o = whtblli4sio4

            ! Steel fraction by volume is the remainder:
            fblss_ccfe = 1.0D0 - fblli2sio4 - fbltibe12 - vfcblkt - vfpblkt

            ! Steel mass (kg)
            whtblss = volblkt * fblss_ccfe * denstl

            ! Total blanket mass (kg)
            whtblkt = whtbltibe12 + whtblli4sio4 + whtblss

            ! Total mass of first wall and blanket
            armour_fw_bl_mass = fw_armour_mass + fwmass + whtblkt

        endif

        if (iblanket==2) then

            ! Mass of steel in blanket (kg)
            if (iblnkith==1) then
                whtblss = denstl * ( volblkti/blnkith * ( blbuith * fblss + blbmith * (1.0D0-fblhebmi) + &
                blbpith * (1.0D0-fblhebpi) ) + volblkto/blnkoth * ( blbuoth * fblss + &
                blbmoth * (1.0D0-fblhebmo) + blbpoth * (1.0D0-fblhebpo) ) )
            else
                whtblss = denstl * ( volblkto/blnkoth * ( blbuoth * fblss + &
                blbmoth * (1.0D0-fblhebmo) + blbpoth * (1.0D0-fblhebpo) ) )
            end if

            ! Mass of beryllium in blanket (kg)
            if (iblnkith==1) then
                whtblbe = 1850.0D0 * fblbe * ( (volblkti * blbuith/blnkith) + &
                    (volblkto * blbuoth/blnkoth) )
            else
                whtblbe = 1850.0D0 * fblbe * (volblkto * blbuoth/blnkoth)
            end if

            ! Mass of breeder material in blanket (kg)
            if (iblnkith==1) then
                whtblbreed = densbreed * fblbreed * ( (volblkti * blbuith/blnkith) + &
                    (volblkto * blbuoth/blnkoth) )
            else
                whtblbreed = densbreed * fblbreed * (volblkto * blbuoth/blnkoth)
            end if

            ! Mass of blanket (kg)
            whtblkt = whtblss + whtblbe + whtblbreed

            ! Void fraction of blanket inboard portion
            if (iblnkith==1) then
                vfblkti = volblkti/volblkt * ( (blbuith/blnkith) * (1.0D0 - fblbe - fblbreed - fblss) &
                    + (blbmith/blnkith) * fblhebmi + (blbpith/blnkith) * fblhebpi )
            else
                vfblkti = 0.0D0
            end if

            ! Void fraction of blanket outboard portion
            vfblkto = volblkto/volblkt * ( (blbuoth/blnkoth) * (1.0D0 - fblbe - fblbreed - fblss) &
                + (blbmoth/blnkoth) * fblhebmo + (blbpoth/blnkoth) * fblhebpo )

            ! Void fraction of blanket
            vfblkt = vfblkti + vfblkto

        endif

    end subroutine component_masses

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! THERMOHYDRAUIC CALCULATIONS
    !!! thermo_hydraulic_model:
    !!!     - blanket_mod_pol_height
    !!!     - flow_velocity
    !!!     - liquid_breeder_properties
    !!!     - primary_coolant_properties
    !!!     - hydraulic_diameter
    !!!     - elbow_coeff
    !!!     - pressure_drop
    !!!     - liquid_breeder_pressure_drop_mhd
    !!!     - deltap_tot
    !!!     - pumppower
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine blanket_mod_pol_height

        !! Calculations for blanket module poloidal height
        !! author: J. Morris, CCFE, Culham Science Centre
        !! Calculations for blanket module poloidal height for D shaped and elliptical machines
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use constants, only: pi
        use build_variables, only: scrapli, scraplo
        use fwbs_variables, only: fwbsshape, nblktmodpi, fdiv, nblktmodpo, &
        iblanket
        use physics_variables, only: itart, rminor, idivrt, rmajor, triang

        implicit none

        ! Local variables

        ! Mid-plane distance from inboard to outboard side (m)
        real(dp) :: a

        ! Internal half-height of blanket (m)
        real(dp) :: b

        ! Calculate ellipse circumference using Ramanujan approximation (m)
        real(dp) :: ptor

        ! Major radius where half-ellipses 'meet' (m)
        real(dp) :: r1

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((itart == 1).or.(fwbsshape == 1)) then  ! D-shaped machine

            ! Segment vertical inboard surface (m)
            bllengi = (2.0D0*hblnkt) / nblktmodpi

            ! Calculate perimeter of ellipse that defines the internal
            ! surface of the outboard first wall / blanket

            ! Mid-plane distance from inboard to outboard side (m)
            a = scrapli + 2.0D0*rminor + scraplo

            ! Internal half-height of blanket (m)
            b = hblnkt

            ! Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = pi * ( 3.0D0*(a+b) - sqrt( (3.0D0*a + b)*(a + 3.0D0*b) ) )

            ! Calculate blanket poloidal length and segment, subtracting divertor length (m)
            ! kit hcll version only had the single null option
            if (idivrt == 2) then
                ! Double null configuration
                bllengo = 0.5D0*ptor * (1.0D0 - 2.0D0*fdiv) / nblktmodpo
            else
                ! single null configuration
                bllengo = 0.5D0*ptor * (1.0D0 - fdiv) / nblktmodpo
            end if

        ! shape defined by two half-ellipses
        else

            ! Major radius where half-ellipses 'meet' (m)
            r1 = rmajor - rminor*triang

            ! Internal half-height of blanket (m)
            b = hblnkt

            ! Distance between r1 and nearest edge of inboard first wall / blanket (m)
            a = r1 - (rmajor - rminor - scrapli)

            ! Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = pi * ( 3.0D0*(a+b) - sqrt( (3.0D0*a + b)*(a + 3.0D0*b) ) )

            ! Calculate inboard blanket poloidal length and segment, subtracting divertor length (m)
            ! Assume divertor lies between the two ellipses, so fraction fdiv still applies

            ! kit hcll version only had the single null option
            if (idivrt == 2) then
                ! Double null configuration
                bllengi = 0.5D0*ptor * (1.0D0 - 2.0D0*fdiv) / nblktmodpi
            else
                ! single null configuration
                bllengi = 0.5D0*ptor * (1.0D0 - fdiv) / nblktmodpi
            end if

            ! Distance between r1 and inner edge of outboard first wall / blanket (m)
            a = rmajor + rminor + scraplo - r1

            ! Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = pi * ( 3.0D0*(a+b) - sqrt( (3.0D0*a + b)*(a + 3.0D0*b) ) )

            ! kit hcll version only had the single null option
            ! Calculate outboard blanket poloidal length and segment, subtracting divertor length (m)
            if (idivrt == 2) then
                ! Double null configuration
                bllengo = 0.5D0*ptor * (1.0D0 - 2.0D0*fdiv) / nblktmodpo
            else
                ! single null configuration
                bllengo = 0.5D0*ptor * (1.0D0 - fdiv) / nblktmodpo
            end if

        end if

    end subroutine blanket_mod_pol_height

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function flow_velocity(i_channel_shape, mass_flow_rate, flow_density)

        !! Calculate the coolant flow velocity (m/s) for given pipe mass flow rate and pipe size/shape.
        !! N.B. Assumed that primary BB and FW coolants have same pipe radius (= afw_outboard).
        !! author: G. Graham, CCFE
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: afw_outboard, a_bz_liq, b_bz_liq
        use constants, only: pi

        implicit none

        !! Function return parameter !!!!!

        real(dp) :: flow_velocity

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        !! Swicth for circular or rectangular channel crossection.
        !! Shape depends on whether primary or secondary coolant.
        !!  - =1   circle (primary)
        !!  - =2   rectangle (secondary)
        integer, intent(in) :: i_channel_shape

        !! Coolant mass flow rate per pipe (kg/s)
        real(dp), intent(in) :: mass_flow_rate

        !! Coolant density
        real(dp), intent(in) :: flow_density

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! If primary coolant then circular channels assumed
        if (i_channel_shape==1) flow_velocity = mass_flow_rate / (flow_density*pi*afw_outboard*afw_outboard)

        !! If secondary coolant then rectangular channels assumed
        if (i_channel_shape==2) flow_velocity = mass_flow_rate / (flow_density * a_bz_liq * b_bz_liq)

    end function flow_velocity

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine liquid_breeder_properties(ip, ofile)

        !! Calculates the fluid properties of the Liquid Metal Breeder/Coolant in the Blanket BZ
        !! Uses middle value of input and output temperatures of Liquid Metal Breeder/Coolant
        !! Curently have PbLi but can expand with e.g., Lithium
        !!
        !! author: G Graham, CCFE
        !!
        !! References:
        !!
        !!      [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
        !!                  lead-lithium alloy, two candidate liquid metal breeder materials
        !!                  for self-cooled blankets, Fusion Engineering and Design 27, 399-406.
        !!
        !!      [Mas2008]   Mas de les Valles et al. (2008), Lead-lithium material database for
        !!                  nuclear fusion technology, Journal of Nuclear Materials, Vol. 376(6).
        !!
        !!      [Mar2019]   Martelli et al. (2019), Literature review of lead-lithium
        !!                  thermophysical properties, Fusion Engineering and Design, 138, 183-195.
        !!
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: inlet_temp_liq, outlet_temp_liq, a_bz_liq, b_bz_liq, den_liq, &
        specific_heat_liq, thermal_conductivity_liq, dynamic_viscosity_liq, electrical_conductivity_liq, &
        i_bb_liq, hartmann_liq, b_mag_blkt, iblnkith
        use physics_variables, only: bt, aspect, rmajor
        use build_variables, only: blnkith, blnkoth
        use error_handling, only: report_error

        implicit none

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        integer, intent(in) :: ip, ofile

        !! Local variables !!!!!!!!!!!!!!!

        !! Gas constant (J K-1 mol-1)
        real(dp)  :: r

        !! mid temp of liquid metal (K)
        real(dp)  :: mid_temp_liq

        !! Ratio of conductivity to dynamic viscosity
        real(dp)  :: con_vsc_rat

        !! Array of valid temperature ranges for breeder propertites
        real(dp)  :: t_ranges(5,2)

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Gas constant
        r = 8.314

        !! Use mid temp
        if (inlet_temp_liq==outlet_temp_liq) then
            mid_temp_liq = outlet_temp_liq
        else
            mid_temp_liq = (inlet_temp_liq + outlet_temp_liq)*0.5
        endif


        !! If the liquid metal is PbLi...
        if (i_bb_liq==0) then

            !! PbLi from [Mar2019]
            !! Constant pressure ~ 17 atmospheres ~ 1.7D6 Pa
            !! Li content is ~ 17%
            !!
            !! density                      kg m-3          T in Kelvin     range = 508-880 K
            !!
            !! specific_heat                J kg-1 K-1      T in Kelvin     range = 508-880 K
            !!
            !! thermal_conductivity         W m-1 K-1       T in Celcius    range = 508-773 K
            !!
            !! dynamic_viscosity            Pa s            T in Celcius    range = 508-873 K
            !!
            !! electrical_conductivity      A V-1 m-1       T in Kelvin     range = 600-800 K

            !! Caculate properties
            den_liq = 1.052D4*(1 - mid_temp_liq*1.13D-4)

            specific_heat_liq = 1.95D2 - mid_temp_liq*9.116D-3

            thermal_conductivity_liq = 1.95 + (mid_temp_liq - 273.15)*1.96D-2

            dynamic_viscosity_liq = 6.11D-3 -(2.257D-5 * (mid_temp_liq-273.15)) &
            + (3.766D-8 * (mid_temp_liq-273.15)**2) - (2.289D-11 * (mid_temp_liq-273.15)**3)

            t_ranges(1:4,1) = 508.0D0
            t_ranges(1:4,2) = 880.0D0

            electrical_conductivity_liq = 1.0D0/(1.03D-6 - (6.75D-11 * mid_temp_liq) + &
            (4.18D-13 * mid_temp_liq**2))

            t_ranges(5,1) = 600.0D0
            t_ranges(5,2) = 800.0D0

        !! If the liquid metal is Li...
        else if (i_bb_liq==1) then

            !! Temporary - should be updated with information from Li reviews conducted at CCFE once completed
            !! Li Properties from [Mal1995] at 300 Celcius
            !! den_liq = 505                           !! kg/m3
            !! specific_heat_liq = 4260                !! J kg-1 K-1
            !! thermal_conductivity_liq = 46           !! W m-1 K-1
            !! dynamic_viscosity_liq = 1.0D-6          !! m2 s-1
            !! electrical_conductivity_liq = 3.03D6    !! A V-1 m-1

            !! New from 'Application of lithium in systems of fusion reactors. 1. Physical and chemical properties of lithium'
            !! Lyublinski et al., 2009, Plasma Devicec and Operations
            den_liq = 504.43D0 - (0.2729D0 * mid_temp_liq) - (8.0035D-5 * mid_temp_liq**2) + (3.799D-8 * mid_temp_liq**3)
            specific_heat_liq = 31.227 + (0.205D6 * mid_temp_liq**(-2)) - (5.265D-3 * mid_temp_liq) + (2.628D6 * mid_temp_liq**(-2))
            !! thermal_conductivity_liq also in paper
            dynamic_viscosity_liq = exp(-4.16D0 - (0.64D0 * log(mid_temp_liq)) + (262.1/mid_temp_liq))
            electrical_conductivity_liq = (0.9249D9 * mid_temp_liq) + 2.3167D6 - (0.7131D3 * mid_temp_liq)

        endif

        !! Magnetic feild strength in T for Hartmann calculation
        !! IB
        if (iblnkith==1) b_mag_blkt(1) = bt * rmajor/(rmajor - (rmajor/aspect) - (blnkith/2))
        !! We do not use this if there is no IB blanket, but will use edge as fill value
        if (iblnkith==0) b_mag_blkt(1) = bt * rmajor/(rmajor - (rmajor/aspect))
        !! OB
        b_mag_blkt(2) = bt * rmajor/(rmajor + (rmajor/aspect) + (blnkoth/2))

        !! Calculate Hartmann number
        con_vsc_rat = electrical_conductivity_liq/dynamic_viscosity_liq
        !! Use toroidal width of the rectangular cooling channel as characteristic length scale
        hartmann_liq = b_mag_blkt * a_bz_liq/2.0D0 * sqrt(con_vsc_rat)

        !! Error for temperature range of breeder property realtions
        if (i_bb_liq == 0) then
           if ((any(t_ranges(:,1) > mid_temp_liq)).or.(any(t_ranges(:,2) < mid_temp_liq))) call report_error(280)
        endif

        !! Output !!!!!!!!!!!!!!!!!!!!!!!!
        if (ip == 0) return
        call write_output_liquid_breeder_properties

    contains

        subroutine write_output_liquid_breeder_properties

            use process_output, only: oheadr, osubhd, ovarrf, ovarre, &
            ocmmnt, ovarin, ovarst

            use fwbs_variables, only: i_bb_liq, den_liq, dynamic_viscosity_liq, &
            electrical_conductivity_liq, icooldual

            implicit none

            !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call oheadr(ofile,'Blanket : Liquid Breeder Properties')

            if (icooldual == 1) call ocmmnt(ofile, 'Single coolant: liquid metal circulted for tritium extraction.')
            if (icooldual == 2) call ocmmnt(ofile, 'Dual coolant: self-cooled liquid metal breeder.')

            if (i_bb_liq == 0) call ocmmnt(ofile, 'Blanket breeder type (i_bb_liq=0), PbLi (~ 17% Li)')
            if (i_bb_liq == 1) call ocmmnt(ofile, 'Blanket breeder type (i_bb_liq=1), Li')

            call ovarrf(ofile, 'Density (kg m-3)', '(den_liq)', den_liq, 'OP ')
            call ovarrf(ofile, 'Viscosity (Pa s)', '(dynamic_viscosity_liq)', dynamic_viscosity_liq, 'OP ')
            call ovarrf(ofile, 'Electrical Conductivity (A V-1 m-1)', '(electrical_conductivity_liq)', electrical_conductivity_liq, 'OP ')
            call ovarrf(ofile, 'Hartmann Number IB', '(hartmann_liq)', hartmann_liq(1), 'OP ')
            call ovarrf(ofile, 'Hartmann Number OB', '(hartmann_liq)', hartmann_liq(2), 'OP ')

            call ovarre(ofile, 'Inlet Temperature (Celcius)', '(inlet_temp_liq)', inlet_temp_liq, 'OP ')
            call ovarre(ofile, 'Outlet Temperature (Celcius)', '(outlet_temp_liq)', outlet_temp_liq, 'OP ')

            if (i_bb_liq == 0) then
                if ((any(t_ranges(:,1) > mid_temp_liq)).or.(any(t_ranges(:,2) < mid_temp_liq))) then
                    call ocmmnt(ofile, 'Outside temperature limit for one or more liquid metal breeder properties.')
                    call ovarrf(ofile, 'Liquid metal temperature (K)', '(mid_temp_liq)', mid_temp_liq, 'OP ')
                    call ocmmnt(ofile, 'Density: Max T = 880 K, Min T = 508 K')
                    call ocmmnt(ofile, 'Specific heat: Max T = 880 K, Min T = 508 K')
                    call ocmmnt(ofile, 'Thermal conductivity: Max T = 880 K, Min T = 508 K')
                    call ocmmnt(ofile, 'Dynamic viscosity : Max T = 880 K, Min T = 508 K')
                    call ocmmnt(ofile, 'Electrical conductivity: Max T = 800 K, Min T = 600 K')
                endif
            endif

        end subroutine write_output_liquid_breeder_properties

    end subroutine liquid_breeder_properties

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function pressure_drop(ip, ofile, i_ps, num_90, num_180, l_pipe, den, vsc, vv, label)

        !! Pressure drops are calculated for a pipe with a number of 90
        !! and 180 degree bends. The pressure drop due to frictional forces along
        !! the total straight length of the pipe is calculated, then the pressure
        !! drop due to the bends is calculated. The total pressure drop is the sum
        !! of all contributions.
        !!
        !! original author: P. J. Knight, CCFE
        !! moved from previous version of pumppower function by: G Graham, CCFE
        !!
        !! N.B Darcy-Weisbach Equation (straight pipe):
        !!
        !!  kstrght = lambda * L/D
        !!
        !!  pressure drop = kstrght * (rho*V^2)/2
        !!
        !!  lambda - Darcy friction factor, L - pipe length, D - hydraulic diameter,
        !!  rho - fluid density, V - fluid flow average velocity
        !!
        !! This function also calculates pressure drop equations for elbow bends,
        !! with modified coefficients.
        !!
        !! N.B. Darcy friction factor is estimated from the Haaland approximation.
        !!
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: afw_outboard, a_bz_liq, b_bz_liq
        use fw_module, only: friction

        implicit none

        !! Function return parameter !!!!!

        !! Pressure drop along the pipe (Pa)
        real(dp) :: pressure_drop

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        integer, intent(in) :: ip, ofile

        !! Swicth for primary or secondary coolant
        !!  - =1   primary
        !!  - =2   secondary
        integer, intent(in) :: i_ps

        !! Number of 90 and 180 degree bends in pipe
        integer, intent(in) :: num_90, num_180

        !! Total flow length along pipe (m)
        real(dp), intent(in) :: l_pipe

        !!  Coolant density (kg/m3)
        real(dp), intent(in) :: den

        !!  Coolant viscosity (Pa s)
        real(dp), intent(in) :: vsc

        !! Coolant flow velocity (m/s)
        real(dp), intent(in) :: vv

        !! Description of this calculation
        character(len=*), intent(in) :: label

        !! Local variables !!!!!!!!!!!!!!!

        !! Hydraulic diameter of coolant flow channels (m)
        real(dp)  :: dh

        !! Reynolds number
        real(dp)  :: reyn

        !! Darcy friction factor
        real(dp)  :: lambda

        !! Pressure drops for straight setions, 90 bends and 180 bends (Pa)
        real(dp)  :: pdropstraight, pdrop90, pdrop180

        !! Elbow radius for 90 and 180 deg brends (m)
        real(dp)  :: elbow_radius

        !! Pressure drop coefficients
        real(dp)  :: kelbwn, kelbwt, kstrght

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Calculate hydraulic dimater for round or retancular pipe (m)
        dh = hydraulic_diameter(i_ps)

        !! Reynolds number
        reyn = den * vv * dh / vsc

        !! Calculate Darcy friction factor
        !! N.B. friction function Uses Haaland approx. which assumes a filled circular pipe.
        !! Use dh which allows us to do fluid calculations for non-cicular tubes
        !! (dh is estimate appropriate for fully developed flow).
        call friction(reyn,lambda)

        !! Pressure drop coefficient !!!!!

        !! Straight section
        kstrght = lambda * l_pipe/dh

        !! In preveious version of pumppower:
        !! - elbow radius assumed = 0.018m for 90 degree elbow, from WCLL
        !! - elbow radius assumed half that of 90 deg case for 180 deg elbow
        !! Intialised value for afw_outboard is 0.006m, so elbow radius = 3 * afw_outboard,
        !! aka 1.5 * pipe diameter, which seems to be engineering standard for
        !! a steel pipe long-radius elbow (short-radius elbow = 2 * afw_outboard).

        !! If primary coolant...
        if (i_ps==1) then
            elbow_radius = 3 * afw_outboard
        !! If secondary coolant...
        else
            !! See DCLL
            elbow_radius = b_bz_liq
        endif

        !! 90 degree elbow pressure drop coefficient
        kelbwn = elbow_coeff(elbow_radius, 90.0D0, lambda, dh)

        !! 180 degree elbow pressure drop coefficient
        kelbwt = elbow_coeff(elbow_radius/2, 180.0D0, lambda, dh)

        !! Total (Pa)
        pdropstraight = kstrght * 0.5D0*den*vv*vv
        pdrop90 = num_90*kelbwn * 0.5D0*den*vv*vv
        pdrop180 = num_180*kelbwt * 0.5D0*den*vv*vv
        pressure_drop = pdropstraight + pdrop90 + pdrop180

        if (ip==0) return
        call write_output_pressure_drop

    contains

        subroutine write_output_pressure_drop

            use process_output, only: oheadr, osubhd, ovarrf, ovarre, &
            ocmmnt, ovarin, ovarst, oblnkl

            use global_variables, only: verbose

            implicit none

            !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call osubhd(ofile, 'Pressure drop (friction) for ' // label)

            call ovarre(ofile, 'Reynolds number', '(reyn)', reyn, 'OP ')
            call ovarre(ofile, 'Darcy friction factor', '(lambda)', lambda, 'OP ')
            call ovarre(ofile, 'Pressure drop (Pa)', '(pressure_drop)', pressure_drop, 'OP ')
            call ocmmnt(ofile, 'This is the sum of the following:')
            call ovarre(ofile, '            Straight sections (Pa)', '(pdropstraight)', &
                pdropstraight, 'OP ')
            call ovarre(ofile, '            90 degree bends (Pa)', '(pdrop90)', pdrop90, 'OP ')
            call ovarre(ofile, '            180 degree bends (Pa)', '(pdrop180)', pdrop180, 'OP ')
            call ocmmnt(ofile, 'Additional information is printed when verbose = 1')
            if (verbose==1) then
                call oblnkl(ofile)
                call ovarre(ofile, 'Straight section pressure drop coefficient', &
                    '(kstrght)', kstrght, 'OP ')
                call ovarre(ofile, '90 degree elbow coefficient', '(kelbwn)', kelbwn, 'OP ')
                call ovarre(ofile, '180 degree elbow coefficient coefficient', '(kelbwt)', kelbwt, 'OP ')
            end if

        end subroutine write_output_pressure_drop

    end function pressure_drop

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function hydraulic_diameter(i_channel_shape)

        !! Caculate the hydraulic diameter (m) for a given coolant pipe size/shape.
        !! author: G. Graham
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: afw_outboard, a_bz_liq, b_bz_liq

        implicit none

        !! Function return parameter !!!!!

        real(dp) :: hydraulic_diameter

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        !! Swicth for circular or rectangular channel crossection
        !! Shape depends on whether primary or secondary coolant
        !!  - =1   circle (primary)
        !!  - =2   rectangle (secondary)
        integer, intent(in) :: i_channel_shape

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! If primary coolant then circular channels assumed
        if (i_channel_shape==1) hydraulic_diameter = 2.0D0*afw_outboard

        !! If secondary coolant then rectangular channels assumed
        if (i_channel_shape==2) hydraulic_diameter = 2*a_bz_liq*b_bz_liq/(a_bz_liq+b_bz_liq)

    end function hydraulic_diameter

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function elbow_coeff(r_elbow, ang_elbow, lambda, dh)

        !! Function calculates elbow bends coefficients for pressure drop
        !! calculations.
        !!
        !! author: G. Graham, CCFE
        !!
        !! References:
        !!
        !!      [Ide1969]   Idel'Cik, I. E. (1969), Memento des pertes de charge,
        !!                  Collection de la Direction des Etudes et Recherches d'Electricité de France.
        !!
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use constants, only: pi

        implicit none

        !! Function return parameter !!!!!

        !! Elbow coefficient for pressure drop calculation
        real(dp) :: elbow_coeff

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        !! Pipe elbow radius (m)
        real(dp), intent(in) :: r_elbow

        !! Pipe elbow angle (degrees)
        real(dp), intent(in) :: ang_elbow

        !! Darcy Friction Factor
        real(dp), intent(in) :: lambda

        !! Hydraulic Diameter (m)
        real(dp), intent(in) :: dh

        !! Local variables !!!!!!!!!!!!!!!

        !! A and B from singularity coefficient equation
        real(dp)  :: a, b

        !! Singularity and friction coefficients
        real(dp)  :: xifn, xift, ximn, ximt

        !! Ratio of elbow radius to Darcy friction factor
        real(dp)  :: r_ratio

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Pressure loss for bends calculated using equivilent lengths.
        !! coeff = lambda * L/D --> get equiv. of L

        if (ang_elbow==90) then
            a = 1.0D0
        else if (ang_elbow<70) then
            a = 0.9D0 * sin(ang_elbow * pi/180.0D0)
        else if (ang_elbow>100) then
            a = 0.7D0 + (0.35D0 * sin((ang_elbow/90.0D0) * (pi/180.0D0)))
        else
            write(*,*) 'No formula for 70 <= elbow angle(deg) <= 100, only 90 deg option available in this range.'
            stop 1
        endif

        r_ratio = r_elbow/dh
        if (r_ratio > 1) b = 0.21D0/r_ratio**0.5
        if (r_ratio < 1) b = 0.21D0/r_ratio**2.5
        if (r_ratio == 1) b = 0.21D0

        !! Singularity
        ximt = a * b

        !! Friction
        xift = 0.0175D0 * lambda * (r_elbow/dh) * ang_elbow

        !! Elbow Coefficient
        elbow_coeff = ximt + xift

    end function elbow_coeff

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function liquid_breeder_pressure_drop_mhd(ip, ofile, vel, vsc, conduct_liq, l_channel, num_pol, label)

        !! Calculates the pressure drop in a liquid metal flow channel due to MHD effects. The total pressure
        !! drop is the sum of contributions. This is only used for secondary coolant/breeder so rectangular flow
        !! channels are assumed.
        !!
        !! author: G Graham, CCFE
        !!
        !! References:
        !!
        !!      [Miy1986]   Miyazaki et al. (1986), Magneto-Hydro-Dynamic Pressure Drop of Lithium
        !!                  Flow in Rectangular Ducts, Fusion Technology, 10:3P2A, 830-836, DOI: 10.13182/FST10-830
        !!
        !!      [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
        !!                  lead-lithium alloy, two candidate liquid metal breeder materials
        !!                  for self-cooled blankets, Fusion Engineering and Design 27, 399-406
        !!
        !!      [Iba2013]   Ibano et al (2013), Nutronics and pumping power analysis on the
        !!                  Tokamak reactor for the fusion-biomass hybrid concept,
        !!                  Fusion Engineering and Design, 88
        !!
        !!      [Sho2018]   Shoki et al (2018), MHD pressure drop measurement of PbLi flow
        !!                  in double-bended pipe, Fusion Engineering and Design, 136, 17-23
        !!
        !!      [Klu2019]   Kluber et al. (2019), Numerical simulations of 3D magnetohydrodynamic
        !!                  flows in dual-coolant lead lithium blankets, Fusion Engineering and Design,
        !!                  146, 684-687
        !!
        !!      [Sua2021]   MHD effects in geometrical sigularities on high velocity breeding
        !!                  blanket designs. Part II, ENR-PRD.BB-T007-D002, EFDA_D_2PDT9U.
        !!                  Also, see asssociated paper: Suarez et al. (2021), On the use of CFD
        !!                  to obtain head loss coefficients in hydraulic systems and it's appliaction
        !!                  to liquid metal flows in nuclear fusion reactor blankets, Plasma. Phys.
        !!                  Control fusion, 63, 124002
        !!
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fwbs_variables, only: ifci, a_bz_liq, b_bz_liq, hartmann_liq, b_mag_blkt, &
        bz_channel_conduct_liq, th_wall_secondary
        use constraint_variables, only: maxradwallload
        use physics_variables, only: btot

        implicit none

        !! Function return parameter !!!!!

        real(dp) :: liquid_breeder_pressure_drop_mhd

        !! Arguments !!!!!!!!!!!!!!!!!!!!!

        integer, intent(in) :: ip, ofile

        !! Liquid metal coolant/breeder  flow velocity (m/s)
        real(dp), intent(in) :: vel

        !! Liquid metal visosity
        real(dp), intent(in) :: vsc

        !! Liquid metal conductivity
        real(dp), intent(in) :: conduct_liq

        !! Length long poloidal sections of channel
        real(dp), intent(in) :: l_channel

        !! Number long poloidal sections of channel
        integer, intent(in) :: num_pol

        !! Description of this calculation
        character(len=*), intent(in) :: label

        !! Local Variables !!!!!!!!!!!!!!!

        !! Half-widths of channel (m)
        real(dp)  :: half_wth_a, half_wth_b

        !! MHD pressure drop for single channel
        real(dp)  :: mhd_pressure_drop

        !! Magnetic field strenght (T)
        real(dp) :: b_mag

        !! Internal resistance of fluid
        real(dp) :: r_i

        !! Resistance of side walls
        real(dp) :: r_w

        !! Resistance of electrode walls
        real(dp) :: r_e

        !! Wall-to-fluid conductance ratio
        real(dp) :: big_c

        !!  Term calculated for [Miy1986] pressure drop
        real(dp) :: kp

        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Magnetic feild strength in IB or OB blanket
        if (label=='Inboard blanket breeder liquid') b_mag = b_mag_blkt(1) !! IB
        if (label=='Outboard blanket breeder liquid') b_mag = b_mag_blkt(2) !! OB

        !! Half-widths
        !! N.B. a_bz_liq (width in the toroidal direction) is in B direction
        half_wth_a = a_bz_liq * 0.5
        half_wth_b = b_bz_liq *0.5

        !! If have thin conducting walls...
        if (ifci/=1) then

            !! Caculate resistances of fluid and walls
            r_i = half_wth_b/(conduct_liq * half_wth_a)
            r_w = half_wth_b/(bz_channel_conduct_liq * th_wall_secondary)
            r_e = half_wth_a/(bz_channel_conduct_liq * th_wall_secondary)
            big_c = r_i/r_w
            !!  Calculate pressure drop for conducting wall [Miy1986]
            kp = big_c/(1 + half_wth_a/(3 * half_wth_b) + big_c)
            mhd_pressure_drop = kp * conduct_liq * vel * (b_mag**2) * l_channel

        !! If have perfcetly insulating FCIs...
        else

            !! Calculate pressure drop for (perfectly) insulating FCI [Mal1995]
            mhd_pressure_drop =  vel * b_mag * l_channel * sqrt(conduct_liq * vsc / half_wth_a)

        endif

        !! Total (Pa)
        liquid_breeder_pressure_drop_mhd = num_pol * mhd_pressure_drop

        if (ip==0) return
        call write_output_liquid_breeder_pressure_drop_mhd

        contains

        subroutine write_output_liquid_breeder_pressure_drop_mhd

            use process_output, only: oheadr, osubhd, ovarrf, ovarre, &
            ocmmnt, ovarin, ovarst

            implicit none

            !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call osubhd(ofile, 'Liquid metal breeder/coolant MHD pressure drop for ' // label)

            if (ifci==0) then
                call ocmmnt(ofile, 'Flow channels have thin conducting walls (ifci==0)')
                call ovarre(ofile, 'Wall conductance (A V-1 m-1)', '(bz_channel_conduct_liq)', bz_channel_conduct_liq, 'OP ')
            else if (ifci==2) then
                call ocmmnt(ofile, 'Flow Channel Inserts (FCIs) used (ifci==2)')
                call ovarre(ofile, 'FCI conductance (A V-1 m-1)', '(bz_channel_conduct_liq)', bz_channel_conduct_liq, 'OP ')
            else
                call ocmmnt(ofile, 'Flow Channel Inserts - assumed perfect insulator (ifci==1)')
            endif
            call ovarre(ofile, 'Length of long poloidal secion of channel (m)', '(l_channel)', l_channel, 'OP ')
            call ovarin(ofile, 'Number of long poloidal secions of channel', '(num_pol)', num_pol, 'OP ')
            call ovarre(ofile, 'MHD pressure drop (Pa)', '(liquid_breeder_pressure_drop_mhd)', liquid_breeder_pressure_drop_mhd, 'OP ')

        end subroutine write_output_liquid_breeder_pressure_drop_mhd

    end function liquid_breeder_pressure_drop_mhd

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module blanket_library
