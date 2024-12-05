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

end module blanket_library
