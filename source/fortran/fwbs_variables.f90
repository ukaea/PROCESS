module fwbs_variables
  !! author: J. Morris (UKAEA), M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the first wall, blanket and
  !! shield components
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: bktlife
  !! Full power blanket lifetime (years)

  real(dp) :: coolmass
  !! mass of water coolant (in shield, blanket, first wall, divertor) [kg]

  real(dp) :: vvmass
  !! vacuum vessel mass [kg]

  real(dp) :: denstl
  !! density of steel [kg m^-3]
  !#TODO: should this be in constants. Is currently an input. Should be a list of preapproved options?

  real(dp) :: denw
  !! density of tungsten [kg m^-3]
  !#TODO: same as above with steel?

  real(dp) :: denwc
  !! density of tungsten carbide [kg m^-3]

  real(dp) :: dewmkg
  !! total mass of vacuum vessel + cryostat [kg] (calculated if blktmodel>0)
  !# TODO: blktmodel needs consolidating with iblanket

  real(dp) :: emult
  !! energy multiplication in blanket and shield

  real(dp) :: emultmw
  !! power due to energy multiplication in blanket and shield [MW]

  real(dp) :: fblss
  !! KIT blanket model: steel fraction of breeding zone

  real(dp) :: fdiv
  !! Solid angle fraction taken by one divertor

  real(dp) :: fhcd
  !! area fraction covered by heating/current drive apparatus plus diagnostics

  real(dp) :: fhole
  !! area fraction taken up by other holes (IFE)

  integer :: fwbsshape
  !! switch for first wall, blanket, shield and vacuum vessel shape:
  !!
  !! - =1 D-shaped (cylinder inboard + ellipse outboard)
  !! - =2 defined by two ellipses
  !#TODO: change to adopt switch naming convention

  real(dp) :: fwlife
  !! first wall full-power year lifetime (y)

  real(dp) :: fwmass
  !! first wall mass [kg]

  real(dp) :: fw_armour_mass
  !! first wall armour mass [kg]

  real(dp) :: fw_armour_thickness
  !! first wall armour thickness [m]

  real(dp) :: fw_armour_vol
  !! first wall armour volume [m^3]

  integer :: iblanket
  !! switch for blanket model:
  !!
  !! - =1 CCFE HCPB model
  !! - =2 KIT HCPB model
  !! - =3 CCFE HCPB model with Tritium Breeding Ratio calculation
  !! - =4 KIT HCLL model
  !! - =5 DCLL model -  no nutronics model included (in development) please check/choose values for
  !!                      'dual-coolant blanket' fractions (provided in this file).
  !!                 -  please use primary_pumping = 0 or 1.

  integer :: iblnkith
  !! switch for inboard blanket:
  !!
  !! - =0 No inboard blanket (blnkith=0.0)
  !! - =1 Inboard blanket present

  integer :: inuclear
  !! switch for nuclear heating in the coils:
  !!
  !! - =0 Frances Fox model (default)
  !! - =1 Fixed by user (qnuc)

  real(dp) :: qnuc
  !! nuclear heating in the coils (W) (`inuclear=1`)

  real(dp) :: li6enrich
  !! lithium-6 enrichment of breeding material (%)

  real(dp) :: pnucblkt
  !! nuclear heating in the blanket [MW]

  real(dp) :: pnuc_cp
  !! Total nuclear heating in the ST centrepost [MW]

  real(dp) :: pnuc_cp_sh
  !! Neutronic shield nuclear heating in the ST centrepost [MW]

  real(dp) :: pnuc_cp_tf
  !! TF neutronic nuclear heating in the ST centrepost [MW]

  real(dp) :: pnucdiv
  !! nuclear heating in the divertor [MW]

  real(dp) :: pnucfw
  !! nuclear heating in the first wall [MW]

  real(dp) :: pnuchcd
  !! nuclear heating in the HCD apparatus and diagnostics [MW]

  real(dp) :: pnucloss
  !! nuclear heating lost via holes [MW]

  real(dp) :: pnucvvplus
  !! nuclear heating to vacuum vessel and beyond [MW]

  real(dp) :: pnucshld
  !! nuclear heating in the shield [MW]

  real(dp) :: whtblkt
  !! mass of blanket [kg]

  real(dp) :: whtblss
  !! mass of blanket - steel part [kg]

  real(dp) :: armour_fw_bl_mass
  !! Total mass of armour, first wall and blanket [kg]

  ! CCFE HCPB Blanket Model (with or without TBR calculation) iblanket=1,3
  ! ----------

  real(dp) :: breeder_f
  !! Volume ratio: Li4SiO4/(Be12Ti+Li4SiO4) (`iteration variable 108`)

  real(dp) :: breeder_multiplier
  !! combined breeder/multipler fraction of blanket by volume

  real(dp) :: vfcblkt
  !! He coolant fraction of blanket by volume (`iblanket= 1,3` (CCFE HCPB))

  real(dp) :: vfpblkt
  !! He purge gas fraction of blanket by volume (`iblanket= 1,3` (CCFE HCPB))

  real(dp) :: whtblli4sio4
  !! mass of lithium orthosilicate in blanket [kg] (`iblanket=1,3` (CCFE HCPB))

  real(dp) :: whtbltibe12
  !! mass of titanium beryllide in blanket [kg] (`iblanket=1,3` (CCFE HCPB))

  real(dp) :: neut_flux_cp
  !! Centrepost TF fast neutron flux (E > 0.1 MeV) [m^(-2).^(-1)]
  !! This variable is only calculated for superconducting (i_tf_sup = 1 )
  !! spherical tokamal magnet designs (itart = 0)

  real(dp) :: f_neut_shield
  !! Fraction of nuclear power shielded before the CP magnet (ST)
  !! ( neut_absorb = -1 --> a fit on simplified MCNP neutronic
  !! calculation is used assuming water cooled (13%) tungesten carbyde )

  real(dp) :: vffwi, vffwo
  !! Inboard/outboard FW coolant void fraction

  real(dp) :: psurffwi, psurffwo
  !! Surface heat flux on first wall [MW] (sum = pradfw)

  real(dp) :: volfw
  !! First wall volume [m3]

  real(dp) :: fblss_ccfe, fblli2sio4, fbltibe12
  !! Fractions of blanket by volume: steel, lithium orthosilicate, titanium beryllide

  !  KIT HCPB blanket model (iblanket = 2)
  ! ----------

  integer :: breedmat
  !! breeder material switch (iblanket=2 (KIT HCPB)):
  !!
  !! - =1 Lithium orthosilicate
  !! - =2 Lithium methatitanate
  !! - =3 Lithium zirconate

  real(dp) :: densbreed
  !! density of breeder material [kg m^-3] (`iblanket=2` (KIT HCPB))

  real(dp) :: fblbe
  !! beryllium fraction of blanket by volume (if `iblanket=2`, is Be fraction of breeding zone)

  real(dp) :: fblbreed
  !! breeder fraction of blanket breeding zone by volume (`iblanket=2` (KIT HCPB))

  real(dp) :: fblhebmi
  !! helium fraction of inboard blanket box manifold by volume (`iblanket=2` (KIT HCPB))

  real(dp) :: fblhebmo
  !! helium fraction of outboard blanket box manifold by volume (`iblanket=2` (KIT HCPB))

  real(dp) :: fblhebpi
  !! helium fraction of inboard blanket back plate by volume (`iblanket=2` (KIT HCPB))

  real(dp) :: fblhebpo
  !! helium fraction of outboard blanket back plate by volume (`iblanket=2` (KIT HCPB))

  integer :: hcdportsize
  !! switch for size of heating/current drive ports (`iblanket=2` (KIT HCPB)):
  !!
  !! - =1 'small'
  !! - =2 'large'
  !#TODO: switch name and also large and small not descriptive enough

  real(dp) :: nflutf
  !! peak fast neutron fluence on TF coil superconductor [n m^-2] (`iblanket=2` (KIT HCPB))

  integer :: npdiv
  !! number of divertor ports (`iblanket=2` (KIT HCPB))

  integer :: nphcdin
  !! number of inboard ports for heating/current drive (`iblanket=2` (KIT HCPB))

  integer :: nphcdout
  !! number of outboard ports for heating/current drive (`iblanket=2` (KIT HCPB))

  real(dp) :: tbr
  !! tritium breeding ratio (`iblanket=2,3` (KIT HCPB/HCLL))

  real(dp) :: tritprate
  !! tritium production rate [g day^-1] (`iblanket=2` (KIT HCPB))

  real(dp) :: wallpf
  !! neutron wall load peaking factor (`iblanket=2` (KIT HCPB))

  real(dp) :: whtblbreed
  !! mass of blanket - breeder part [kg] (`iblanket=2` (KIT HCPB))

  real(dp) :: whtblbe
  !! mass of blanket - beryllium part [kg]

  ! CCFE HCPB model with Tritium Breeding Ratio calculation (iblanket=3)
  ! ---------------

  integer :: iblanket_thickness
  !! Blanket thickness switch (Do not set blnkith, blnkoth, fwith or fwoth when `iblanket=3`):
  !!
  !! - =1 thin    0.53 m inboard, 0.91 m outboard
  !! - =2 medium  0.64 m inboard, 1.11 m outboard
  !! - =3 thick   0.75 m inboard, 1.30 m outboard

  integer :: primary_pumping
  !! Switch for pumping power for primary coolant (mechanical power only and peak first wall
  !! temperature is only calculated if `primary_pumping=2`):
  !!
  !! - =0 User sets pump power directly (htpmw_blkt, htpmw_fw, htpmw_div, htpmw_shld)
  !! - =1 User sets pump power as a fraction of thermal power (fpumpblkt, fpumpfw, fpumpdiv, fpumpshld)
  !! - =2 Mechanical pumping power is calculated
  !! - =3 Mechanical pumping power is calculated using specified pressure drop

  integer :: i_shield_mat
  !! Switch for shield material - *currently only applied in costing routines* `cost_model = 2`
  !!
  !! - =0 Tungsten (default)
  !! - =1 Tungsten carbide

  integer :: secondary_cycle
  !! Switch for power conversion cycle:
  !!
  !! - =0 Set efficiency for chosen blanket, from detailed models (divertor heat not used)
  !! - =1 Set efficiency for chosen blanket, from detailed models (divertor heat used)
  !! - =2 user input thermal-electric efficiency (etath)
  !! - =3 steam Rankine cycle
  !! - =4 supercritical CO2 cycle

  integer :: secondary_cycle_liq
  !! Switch for power conversion cycle for the liquid breeder component of the blanket:
  !!
  !! - =2 user input thermal-electric efficiency (etath)
  !! - =4 supercritical CO2 cycle

  integer :: coolwh
  !! Switch for blanket coolant (set via blkttype):
  !!
  !! - =1 helium
  !! - =2 pressurized water
  !#TODO: change switch name to satisfy convention

  real(dp) :: afwi
  !! inner radius of inboard first wall/blanket coolant channels (stellarator only) [m]
  !#TODO move to stellarator?

  real(dp) :: afwo
  !! inner radius of outboard first wall/blanket coolant channels (stellarator only) [m]
  !#TODO move to stellarator?

  character(len=6) :: fwcoolant
  !! switch for first wall coolant (can be different from blanket coolant):
  !!
  !! - 'helium'
  !! - 'water'

  real(dp) :: fw_wall
  !! wall thickness of first wall coolant channels [m]

  real(dp) :: afw
  !! radius of first wall cooling channels [m]

  real(dp) :: pitch
  !! pitch of first wall cooling channels [m]

  real(dp) :: fwinlet
  !! inlet temperature of first wall coolant [K]

  real(dp) :: fwoutlet
  !! outlet temperature of first wall coolant [K]

  real(dp) :: fwpressure
  !! first wall coolant pressure [Pa] (`secondary_cycle>1`)

  real(dp) :: tpeak
  !! peak first wall temperature [K]

  real(dp) :: roughness
  !! first wall channel roughness epsilon [m]

  real(dp) :: fw_channel_length
  !! Length of a single first wall channel (all in parallel) [m]
  !! (`iteration variable 114`, useful for `constraint equation 39`)

  real(dp) :: peaking_factor
  !! peaking factor for first wall heat loads. (Applied separately to inboard and outboard loads.
  !! Applies to both neutron and surface loads. Only used to calculate peak temperature - not
  !! the coolant flow rate.)

  real(dp) :: blpressure
  !! blanket coolant pressure [Pa] (`secondary_cycle>1`)

  real(dp) :: inlet_temp
  !! inlet temperature of blanket coolant  [K] (`secondary_cycle>1`)

  real(dp) :: outlet_temp
  !! Outlet temperature of blanket coolant [K] (`secondary_cycle>1`)
  !!
  !! - input if `coolwh=1` (helium)
  !! - calculated if `coolwh=2` (water)

  real(dp) :: coolp
  !! blanket coolant pressure [Pa] (stellarator only)

  integer :: nblktmodpo
  !! number of outboard blanket modules in poloidal direction (`secondary_cycle>1`)

  integer :: nblktmodpi
  !! number of inboard blanket modules in poloidal direction (`secondary_cycle>1`)

  integer :: nblktmodto
  !! number of outboard blanket modules in toroidal direction (`secondary_cycle>1`)

  integer :: nblktmodti
  !! number of inboard blanket modules in toroidal direction (`secondary_cycle>1`)

  real(dp) :: tfwmatmax
  !! maximum temperature of first wall material [K] (`secondary_cycle>1`)

  real(dp) :: fw_th_conductivity
  !! thermal conductivity of first wall material at 293 K (W/m/K) (Temperature dependence
  !! is as for unirradiated Eurofer)

  real(dp) :: fvoldw
  !! area coverage factor for vacuum vessel volume

  real(dp) :: fvolsi
  !! area coverage factor for inboard shield volume

  real(dp) :: fvolso
  !! area coverage factor for outboard shield volume

  real(dp) :: fwclfr
  !! first wall coolant fraction (calculated if `lpulse=1` or `ipowerflow=1`)

  real(dp) :: praddiv
  !! Radiation power incident on the divertor (MW)

  real(dp) :: pradfw
  !! Radiation power incident on the first wall (MW)

  real(dp) :: pradhcd
  !! Radiation power incident on the heating and current drive system (MW)

  real(dp) :: pradloss
  !! Radiation power lost through holes (eventually hits shield) (MW)
  !! Only used for stellarator

  real(dp) :: ptfnuc
  !! nuclear heating in the TF coil (MW)

  real(dp) :: ptfnucpm3
  !! nuclear heating in the TF coil (MW/m3) (`blktmodel>0`)
  !#TODO: check usage of old blktmodel. Update to iblanket

  real(dp) :: rdewex
  !! cryostat radius [m]

  real(dp) :: zdewex
  !! cryostat height [m]

  real(dp) :: rpf2dewar
  !! radial distance between outer edge of largest (`ipfloc=3`) PF coil (or stellarator
  !! modular coil) and cryostat [m]

  real(dp) :: vdewex
  !! cryostat volume [m^3]

  real(dp) :: vdewin
  !! vacuum vessel volume [m^3]

  real(dp) :: vfshld
  !! coolant void fraction in shield

  real(dp) :: volblkt
  !! volume of blanket [m^3]

  real(dp) :: volblkti
  !! volume of inboard blanket [m^3]

  real(dp) :: volblkto
  !! volume of outboard blanket [m^3]

  real(dp) :: volshld
  !! volume of shield [m^3]

  real(dp) :: whtshld
  !! mass of shield [kg]

  real(dp) :: wpenshld
  !! mass of the penetration shield [kg]

  real(dp) :: wtshldi
  !! mass of inboard shield [kg]

  real(dp) :: wtshldo
  !! mass of outboard shield [kg]

  integer :: irefprop
  !! Switch to use REFPROP routines (stellarator only)
  !#TODO: number of stellarator only items here. Also appear in fispact. Tidy needed

  real(dp) :: fblli
  !! lithium fraction of blanket by volume (stellarator only)

  real(dp) :: fblli2o
  !! lithium oxide fraction of blanket by volume (stellarator only)

  real(dp) :: fbllipb
  !! lithium lead fraction of blanket by volume (stellarator only)

  real(dp) :: fblvd
  !! vanadium fraction of blanket by volume (stellarator only)

  real(dp) :: wtblli2o
  !! mass of blanket - Li_2O part [kg]

  real(dp) :: wtbllipb
  !! mass of blanket - Li-Pb part [kg]

  real(dp) :: whtblvd
  !! mass of blanket - vanadium part [kg]

  real(dp) :: whtblli
  !! mass of blanket - lithium part [kg]

  real(dp) :: vfblkt
  !! coolant void fraction in blanket (`blktmodel=0`), (calculated if `blktmodel > 0`)

  integer :: blktmodel
  !! switch for blanket/tritium breeding model (see iblanket):
  !!
  !! - =0 original simple model
  !! - =1 KIT model based on a helium-cooled pebble-bed blanket (HCPB) reference design
  !#TODO: this needs investigating and removing after any required functionality is in iblanket

  real(dp) :: declblkt
  !! neutron power deposition decay length of blanket structural material [m] (stellarators only)

  real(dp) :: declfw
  !! neutron power deposition decay length of first wall structural material [m] (stellarators only)

  real(dp) :: declshld
  !! neutron power deposition decay length of shield structural material [m] (stellarators only)

  integer :: blkttype
  !! Switch for blanket type:
  !!
  !! - =1 WCLL;
  !! - =2 HCLL; efficiency taken from M. Kovari 2016
  !! "PROCESS": A systems code for fusion power plants - Part 2: Engineering
  !! https://www.sciencedirect.com/science/article/pii/S0920379616300072
  !! Feedheat & reheat cycle assumed
  !! - =3 HCPB; efficiency taken from M. Kovari 2016
  !! "PROCESS": A systems code for fusion power plants - Part 2: Engineering
  !! https://www.sciencedirect.com/science/article/pii/S0920379616300072
  !! Feedheat & reheat cycle assumed
  !#TODO: this needs to be merged into iblanket and then removed.

  real(dp) :: etaiso
  !! isentropic efficiency of FW and blanket coolant pumps

  real(dp) :: etahtp
  !! electrical efficiency of primary coolant pumps

  !! -----------------------------------------------------
  !! BLANKET REFACTOR
  !! For DCLL, but to be used by all mods that share blanket library after testing.
  !! Thermodynamic Model for primary_pumping == 2
  !! -----------------------------------------------------

  !! Switch for whether the FW and BB are on the same pump system
  !! i.e. do they have the same primary coolant or not
  !!  - =0    FW and BB have the same primary coolant, flow = FWin->FWout->BBin->BBout
  !!  - =1    FW and BB have the different primary coolant and are on different pump systems
  integer :: ipump

  !! Switch for Liquid Metal Breeder Material
  !!  - =0   PbLi
  !!  - =1   Li
  integer :: i_bb_liq

  !! Switch to specify whether breeding blanket is single-cooled or dual-coolant.
  !!  - =0    Single coolant used for FW and Blanket (H2O or He). Solid Breeder.
  !!  - =1    Single coolant used for FW and Blanket (H2O or He). Liquid metal breeder
  !!          circulted for tritium extraction.
  !!  - =2    Dual coolant: primary coolant (H2O or He) for FW and blanket structure;
  !!          secondary coolant is self-cooled liquid metal breeder.
  integer :: icooldual

  !! Switch for Flow Channel Insert (FCI) type if liquid metal breeder blanket.
  !!  - =0    Thin conducting walls, default electrical conductivity (bz_channel_conduct_liq) is Eurofer
  !!  - =1    Insulating Material, assumed perfect electrical insulator, default density (den_ceramic) is for SiC
  !!  - =2    Insulating Material, electrical conductivity (bz_channel_conduct_liq) is input (default Eurofer), default density (den_ceramic) is for SiC
  integer :: ifci

  !! Switch for Multi Module Segment (MMS) or Single Modle Segment (SMS)
  !!  - =0    MMS
  !!  - =1    SMS
  integer :: ims

  !! Number of liquid metal breeder recirculations per day, for use with icooldual=1
  integer :: n_liq_recirc

  !! Radial fraction of BZ liquid channels
  real(dp) :: r_f_liq_ib, r_f_liq_ob

  !! Toroidal fraction of BZ liquid channels
  real(dp) :: w_f_liq_ib, w_f_liq_ob

  !! FCI material density
  real(dp) :: den_ceramic

  !! Liquid metal coolant/breeder wall thickness thin conductor or FCI [m]
  real(dp) :: th_wall_secondary

  !! Liquid metal coolant/breeder thin conductor or FCI wall conductance [A V^-1 m^-1]
  real(dp) :: bz_channel_conduct_liq

  !! Toroidal width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone
  real(dp) :: a_bz_liq

  !! Radial width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone
  real(dp) :: b_bz_liq

  !! Number of poloidal sections in a liquid metal breeder/coolant channel for module/segment
  integer :: nopol

  !! Number of Liquid metal breeder/coolant channels per module/segment
  integer :: nopipes

  !! Liquid metal breeder/coolant density [kg m^-3]
  real(dp) :: den_liq

  !! Liquid metal
  real(dp) ::wht_liq, wht_liq_ib, wht_liq_ob

  !! Liquid metal breeder/coolant specific heat [J kg^-1 K^-1]
  real(dp) :: specific_heat_liq

  !! Liquid metal breeder/coolant thermal conductivity [W m^-1 K^-1]
  real(dp) :: thermal_conductivity_liq

  !! Liquid metal breeder/coolant dynamic viscosity [Pa s]
  real(dp) :: dynamic_viscosity_liq

  !! Liquid metal breeder/coolant electrical conductivity [Ohm m]
  real(dp) :: electrical_conductivity_liq

  !! Hartmann number
  real(dp) :: hartmann_liq(2)

  !! Toroidal Magnetic feild strength for IB/OB blanket [T]
  real(dp) :: b_mag_blkt(2)

  !! Isentropic efficiency of blanket liquid breeder/coolant pumps
  real(dp) :: etaiso_liq

  !! blanket liquid metal breeder/coolant pressure [Pa]
  real(dp) :: blpressure_liq

  !! Inlet (scan var 68) and Outlet (scan var 69) temperature of the liquid breeder/coolant [K]
  real(dp) :: inlet_temp_liq, outlet_temp_liq

  !! Density of the FW primary coolant
  real(dp) :: rhof_fw

  !! Viscosity of the FW primary coolant
  real(dp) :: visc_fw

  !! Density of the blanket primary coolant
  real(dp) :: rhof_bl

  !! Viscosity of the blanket primary coolant
  real(dp) :: visc_bl

  !! Spesific heat for FW and blanket primary coolant(s)
  real(dp) :: cp_fw, cv_fw, cp_bl, cv_bl

  !! For a dual-coolant blanket, fraction of BZ power cooled by primary coolant
  real(dp) :: f_nuc_pow_bz_struct

  !! For a dual-coolant blanket, fraction of BZ self-cooled power (secondary coolant)
  real(dp) :: f_nuc_pow_bz_liq

  !! For a dual-coolant blanket, ratio of FW/Blanket nuclear power as fraction of total
  real(dp) :: pnuc_fw_ratio_dcll, pnuc_blkt_ratio_dcll

  !! Number of radial and poloidal sections that make up the total primary coolant flow
  !! length in a blanket module (IB and OB)
  integer :: bzfllengi_n_rad, bzfllengi_n_pol, bzfllengo_n_rad, bzfllengo_n_pol

  !! Number of radial and poloidal sections that make up the total secondary coolant/breeder
  !! flow length in a blanket module (IB and OB)
  integer :: bzfllengi_n_rad_liq, bzfllengi_n_pol_liq, bzfllengo_n_rad_liq, bzfllengo_n_pol_liq

  contains

  subroutine init_fwbs_variables
    !! Initialise fwbs variables
    implicit none

    bktlife = 0.0D0
    coolmass = 0.0D0
    vvmass = 0.0D0
    denstl = 7800.0D0
    denw = 19250.0D0
    denwc = 15630.0D0
    dewmkg = 0.0D0
    emult = 1.269D0
    emultmw = 0.0D0
    fblss = 0.09705D0
    fdiv = 0.115D0
    fhcd = 0.0D0
    fhole = 0.0D0
    fwbsshape = 2
    fwlife = 0.0D0
    fwmass = 0.0D0
    fw_armour_mass = 0.0D0
    fw_armour_thickness = 0.005D0
    fw_armour_vol = 0.0D0
    iblanket = 1
    iblnkith = 1
    inuclear = 0
    qnuc = 0.0D0
    li6enrich = 30.0D0
    pnucblkt = 0.0D0
    pnucdiv = 0.0D0
    pnucfw = 0.0D0
    pnuchcd = 0.0D0
    pnucloss = 0.0D0
    pnucvvplus = 0.0D0
    pnucshld = 0.0D0
    whtblkt = 0.0D0
    whtblss = 0.0D0
    armour_fw_bl_mass = 0.0D0
    breeder_f = 0.5D0
    breeder_multiplier = 0.75D0
    vfcblkt = 0.05295D0
    vfpblkt = 0.1D0
    whtblli4sio4 = 0.0D0
    whtbltibe12 = 0.0D0
    f_neut_shield = -1.0D0
    vffwi = 0.0D0
    vffwo = 0.0D0
    psurffwi = 0.0D0
    psurffwo = 0.0D0
    volfw = 0.0D0
    fblss_ccfe = 0.0D0
    fblli2sio4 = 0.0D0
    fbltibe12 = 0.0D0
    breedmat = 1
    densbreed = 0.0D0
    fblbe = 0.6D0
    fblbreed = 0.154D0
    fblhebmi = 0.4D0
    fblhebmo = 0.4D0
    fblhebpi = 0.6595D0
    fblhebpo = 0.6713D0
    hcdportsize = 1
    nflutf = 0.0D0
    npdiv = 2
    nphcdin = 2
    nphcdout = 2
    tbr = 0.0D0
    tritprate = 0.0D0
    wallpf = 1.21D0
    whtblbreed = 0.0D0
    whtblbe = 0.0D0
    iblanket_thickness = 2
    primary_pumping = 2
    i_shield_mat = 0
    secondary_cycle = 0
    secondary_cycle_liq = 4
    coolwh = 1
    afwi = 0.008D0
    afwo = 0.008D0
    fwcoolant = 'helium'
    fw_wall = 0.003D0
    afw = 0.006D0
    pitch = 0.02D0
    fwinlet = 573.0D0
    fwoutlet = 823.0D0
    fwpressure = 15.5D6
    tpeak = 873.0D0
    roughness = 1.0D-6
    fw_channel_length = 4.0D0
    peaking_factor = 1.0D0
    blpressure = 15.50D6
    inlet_temp = 573.0D0
    outlet_temp = 823.0D0
    coolp = 15.5D6
    nblktmodpo = 8
    nblktmodpi = 7
    nblktmodto = 48
    nblktmodti = 32
    tfwmatmax = 823.0D0
    fw_th_conductivity = 28.34D0
    fvoldw = 1.74D0
    fvolsi = 1.0D0
    fvolso = 0.64D0
    fwclfr = 0.15D0
    praddiv = 0.0D0
    pradfw = 0.0D0
    pradhcd = 0.0D0
    pradloss = 0.0D0
    ptfnuc = 0.0D0
    ptfnucpm3 = 0.0D0
    rdewex = 0.0D0
    zdewex = 0.0D0
    rpf2dewar = 0.5D0
    vdewex = 0.0D0
    vdewin = 0.0D0
    vfshld = 0.25D0
    volblkt = 0.0D0
    volblkti = 0.0D0
    volblkto = 0.0D0
    volshld = 0.0D0
    whtshld = 0.0D0
    wpenshld = 0.0D0
    wtshldi = 0.0D0
    wtshldo = 0.0D0
    irefprop = 1
    fblli = 0.0D0
    fblli2o = 0.08D0
    fbllipb = 0.68D0
    fblvd = 0.0D0
    wtblli2o = 0.0D0
    wtbllipb = 0.0D0
    whtblvd = 0.0D0
    whtblli = 0.0D0
    vfblkt = 0.25D0
    blktmodel = 0
    declblkt = 0.075D0
    declfw = 0.075D0
    declshld = 0.075D0
    blkttype = 3
    etaiso = 0.85D0
    etahtp = 0.95D0
    pnuc_cp = 0.0D0
    pnuc_cp_sh = 0.0D0
    pnuc_cp_tf = 0.0D0
    neut_flux_cp = 0.0D0
    ipump = 0
    i_bb_liq = 0
    icooldual = 0
    ifci = 0
    ims = 0
    n_liq_recirc = 10
    r_f_liq_ib=0.5
    r_f_liq_ob=0.5
    w_f_liq_ib=0.5
    w_f_liq_ob=0.5
    den_ceramic = 3.21D3
    th_wall_secondary = 1.25D-2
    bz_channel_conduct_liq = 8.33D5
    a_bz_liq = 0.2D0
    b_bz_liq = 0.2D0
    nopol = 2
    nopipes = 4
    den_liq = 9.5D3
    specific_heat_liq = 1.9D2
    thermal_conductivity_liq = 30.0
    wht_liq = 0.0D0
    wht_liq_ib = 0.0D0
    wht_liq_ob = 0.0D0
    dynamic_viscosity_liq = 0.0D0
    electrical_conductivity_liq = 0.0D0
    hartmann_liq = (0.0D0, 0.0D0)
    b_mag_blkt = (5.0D0, 5.0D0)
    etaiso_liq = 0.85D0
    blpressure_liq = 1.7D6
    inlet_temp_liq = 570.0D0
    outlet_temp_liq = 720.0D0
    rhof_fw = 0.0D0
    visc_fw = 0.0D0
    rhof_bl = 0.0D0
    visc_bl = 0.0D0
    cp_fw = 0.0D0
    cv_fw = 0.0D0
    cp_bl = 0.0D0
    cv_bl = 0.0D0
    f_nuc_pow_bz_struct = 0.34D0
    f_nuc_pow_bz_liq = 0.66D0
    pnuc_fw_ratio_dcll = 0.14D0
    pnuc_blkt_ratio_dcll = 0.86D0
    bzfllengi_n_rad = 4
    bzfllengi_n_pol = 2
    bzfllengo_n_rad = 4
    bzfllengo_n_pol = 2
    bzfllengi_n_rad_liq = 2
    bzfllengi_n_pol_liq = 2
    bzfllengo_n_rad_liq = 2
    bzfllengo_n_pol_liq = 2

  end subroutine init_fwbs_variables
end module fwbs_variables
