module fwbs_variables
  !! author: J. Morris (UKAEA), M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the first wall, blanket and
  !! shield components
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: life_blkt_fpy
  !! Full power blanket lifetime (years)

  real(dp) :: life_blkt
  !! Calendar year blanket lifetime (years)

  real(dp) :: m_fw_blkt_div_coolant_total
  !! mass of water coolant (in shield, blanket, first wall, divertor) [kg]

  real(dp) :: m_vv
  !! vacuum vessel mass [kg]

  real(dp) :: denstl
  !! density of steel [kg m^-3]
  !#TODO: should this be in constants. Is currently an input. Should be a list of preapproved options?

  real(dp) :: denwc
  !! density of tungsten carbide [kg m^-3]

  real(dp) :: dewmkg
  !! total mass of vacuum vessel + cryostat [kg] (calculated if blktmodel>0)
  !# TODO: blktmodel needs consolidating with i_blanket_type

  real(dp) :: f_p_blkt_multiplication
  !! energy multiplication in blanket and shield

  real(dp) :: p_blkt_multiplication_mw
  !! power due to energy multiplication in blanket and shield [MW]

  real(dp) :: fblss
  !! KIT blanket model: steel fraction of breeding zone

  real(dp) :: f_ster_div_single
  !! Solid angle fraction taken by one divertor

  real(dp) :: f_a_fw_outboard_hcd
  !! area fraction of first wall covered by heating/current drive apparatus plus diagnostics

  real(dp) :: fhole
  !! area fraction taken up by other holes (IFE)

  integer :: i_fw_blkt_vv_shape
  !! switch for first wall, blanket, shield and vacuum vessel shape:
  !!
  !! - =1 D-shaped (cylinder inboard + ellipse outboard)
  !! - =2 defined by two ellipses

  real(dp) :: life_fw_fpy
  !! first wall full-power year lifetime (y)

  real(dp) :: m_fw_total
  !! first wall mass [kg]

  real(dp) :: fw_armour_mass
  !! first wall armour mass [kg]

  real(dp) :: fw_armour_thickness
  !! first wall armour thickness [m]

  real(dp) :: fw_armour_vol
  !! first wall armour volume [m^3]

  integer :: i_blanket_type
  !! switch for blanket model:
  !!
  !! - =1 CCFE HCPB model
  !! - =2 KIT HCPB model  # REMOVED, no longer usable
  !! - =3 CCFE HCPB model with Tritium Breeding Ratio calculation
  !! - =4 KIT HCLL model  # REMOVED, no longer usable
  !! - =5 DCLL model -  no nutronics model included (in development) please check/choose values for
  !!                      'dual-coolant blanket' fractions (provided in this file).
  !!                 -  please use i_coolant_pumping = 0 or 1.

  integer :: i_blkt_inboard
  !! switch for inboard blanket:
  !!
  !! - =0 No inboard blanket (dr_blkt_inboard=0.0)
  !! - =1 Inboard blanket present

  integer :: inuclear
  !! switch for nuclear heating in the coils:
  !!
  !! - =0 Frances Fox model (default)
  !! - =1 Fixed by user (qnuc)

  real(dp) :: qnuc
  !! nuclear heating in the coils (W) (`inuclear=1`)

  real(dp) :: f_blkt_li6_enrichment
  !! lithium-6 enrichment of breeding material (%)

  real(dp) :: p_blkt_nuclear_heat_total_mw
  !! nuclear heating in the blanket [MW]

  real(dp) :: pnuc_cp
  !! Total nuclear heating in the ST centrepost [MW]

  real(dp) :: p_cp_shield_nuclear_heat_mw
  !! Neutronic shield nuclear heating in the ST centrepost [MW]

  real(dp) :: pnuc_cp_tf
  !! TF neutronic nuclear heating in the ST centrepost [MW]

  real(dp) :: p_div_nuclear_heat_total_mw
  !! nuclear heating in the divertor [MW]

  real(dp) :: p_fw_nuclear_heat_total_mw
  !! nuclear heating in the first wall [MW]

  real(dp) :: p_fw_hcd_nuclear_heat_mw
  !! Nuclear heating in the HCD apparatus and diagnostics on the first wall [MW]

  real(dp) :: pnucloss
  !! nuclear heating lost via holes [MW]

  real(dp) :: pnucvvplus
  !! nuclear heating to vacuum vessel and beyond [MW]

  real(dp) :: p_shld_nuclear_heat_mw
  !! nuclear heating in the shield [MW]

  real(dp) :: m_blkt_total
  !! mass of blanket [kg]

  real(dp) :: m_blkt_steel_total
  !! mass of blanket - steel part [kg]

  real(dp) :: armour_fw_bl_mass
  !! Total mass of armour, first wall and blanket [kg]

  ! CCFE HCPB Blanket Model (with or without TBR calculation) i_blanket_type=1,3
  ! ----------

  real(dp) :: breeder_f
  !! Volume ratio: Li4SiO4/(Be12Ti+Li4SiO4) (`iteration variable 108`)

  real(dp) :: breeder_multiplier
  !! combined breeder/multipler fraction of blanket by volume

  real(dp) :: vfcblkt
  !! He coolant fraction of blanket by volume (`i_blanket_type= 1,3` (CCFE HCPB))

  real(dp) :: vfpblkt
  !! He purge gas fraction of blanket by volume (`i_blanket_type= 1,3` (CCFE HCPB))

  real(dp) :: m_blkt_li4sio4
  !! mass of lithium orthosilicate in blanket [kg] (`i_blanket_type=1,3` (CCFE HCPB))

  real(dp) :: m_blkt_tibe12
  !! mass of titanium beryllide in blanket [kg] (`i_blanket_type=1,3` (CCFE HCPB))

  real(dp) :: neut_flux_cp
  !! Centrepost TF fast neutron flux (E > 0.1 MeV) [m^(-2).^(-1)]
  !! This variable is only calculated for superconducting (i_tf_sup = 1 )
  !! spherical tokamal magnet designs (itart = 0)

  real(dp) :: f_neut_shield
  !! Fraction of nuclear power shielded before the CP magnet (ST)
  !! ( neut_absorb = -1 --> a fit on simplified MCNP neutronic
  !! calculation is used assuming water cooled (13%) tungesten carbyde )

  real(dp) :: f_a_fw_coolant_inboard, f_a_fw_coolant_outboard
  !! Inboard/outboard FW coolant cross-sectional area void fraction

  real(dp) :: psurffwi, psurffwo
  !! Surface heat flux on first wall [MW] (sum = p_fw_rad_total_mw)

  real(dp) :: vol_fw_total
  !! First wall volume [m3]

  real(dp) :: f_vol_blkt_steel, f_vol_blkt_li4sio4, f_vol_blkt_tibe12
  !! Fractions of blanket by volume: steel, lithium orthosilicate, titanium beryllide

  !  KIT HCPB blanket model (i_blanket_type = 2)
  ! ----------

  integer :: breedmat
  !! breeder material switch (i_blanket_type=2 (KIT HCPB)):
  !!
  !! - =1 Lithium orthosilicate
  !! - =2 Lithium methatitanate
  !! - =3 Lithium zirconate

  real(dp) :: densbreed
  !! density of breeder material [kg m^-3] (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: fblbe
  !! beryllium fraction of blanket by volume (if `i_blanket_type=2`, is Be fraction of breeding zone)

  real(dp) :: fblbreed
  !! breeder fraction of blanket breeding zone by volume (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: fblhebmi
  !! helium fraction of inboard blanket box manifold by volume (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: fblhebmo
  !! helium fraction of outboard blanket box manifold by volume (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: fblhebpi
  !! helium fraction of inboard blanket back plate by volume (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: fblhebpo
  !! helium fraction of outboard blanket back plate by volume (`i_blanket_type=2` (KIT HCPB))

  integer :: hcdportsize
  !! switch for size of heating/current drive ports (`i_blanket_type=2` (KIT HCPB)):
  !!
  !! - =1 'small'
  !! - =2 'large'
  !#TODO: switch name and also large and small not descriptive enough

  real(dp) :: nflutf
  !! peak fast neutron fluence on TF coil superconductor [n m^-2] (`i_blanket_type=2` (KIT HCPB))

  integer :: npdiv
  !! number of divertor ports (`i_blanket_type=2` (KIT HCPB))

  integer :: nphcdin
  !! number of inboard ports for heating/current drive (`i_blanket_type=2` (KIT HCPB))

  integer :: nphcdout
  !! number of outboard ports for heating/current drive (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: tbr
  !! tritium breeding ratio (`i_blanket_type=2,3` (KIT HCPB/HCLL))

  real(dp) :: tritprate
  !! tritium production rate [g day^-1] (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: wallpf
  !! neutron wall load peaking factor (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: whtblbreed
  !! mass of blanket - breeder part [kg] (`i_blanket_type=2` (KIT HCPB))

  real(dp) :: m_blkt_beryllium
  !! mass of blanket - beryllium part [kg]

  integer :: i_coolant_pumping
  !! Switch for pumping power for primary coolant (mechanical power only and peak first wall
  !! temperature is only calculated if `i_coolant_pumping=2`):
  !!
  !! - =0 User sets pump power directly (p_blkt_coolant_pump_mw, p_fw_coolant_pump_mw, p_div_coolant_pump_mw, p_shld_coolant_pump_mw)
  !! - =1 User sets pump power as a fraction of thermal power (fpumpblkt, fpumpfw, fpumpdiv, fpumpshld)
  !! - =2 Mechanical pumping power is calculated
  !! - =3 Mechanical pumping power is calculated using specified pressure drop

  integer :: i_shield_mat
  !! Switch for shield material - *currently only applied in costing routines* `cost_model = 2`
  !!
  !! - =0 Tungsten (default)
  !! - =1 Tungsten carbide

  integer :: i_thermal_electric_conversion
  !! Switch for power conversion cycle:
  !!
  !! - =0 Set efficiency for chosen blanket, from detailed models (divertor heat not used)
  !! - =1 Set efficiency for chosen blanket, from detailed models (divertor heat used)
  !! - =2 user input thermal-electric efficiency (eta_turbine)
  !! - =3 steam Rankine cycle
  !! - =4 supercritical CO2 cycle

  integer :: secondary_cycle_liq
  !! Switch for power conversion cycle for the liquid breeder component of the blanket:
  !!
  !! - =2 user input thermal-electric efficiency (eta_turbine)
  !! - =4 supercritical CO2 cycle

  integer :: i_blkt_coolant_type
  !! Switch for blanket coolant (set via blkttype):
  !!
  !! - =1 helium
  !! - =2 pressurized water
  !#TODO: change switch name to satisfy convention

  character(len=6) :: i_fw_coolant_type
  !! switch for first wall coolant (can be different from blanket coolant):
  !!
  !! - 'helium'
  !! - 'water'

  real(dp) :: dr_fw_wall
  !! wall thickness of first wall coolant channels [m]

  real(dp) :: radius_fw_channel
  !! radius of first wall cooling channels [m]

  real(dp) :: dx_fw_module
  !! Width of a FW module containing a cooling channel [m]

  real(dp) :: temp_fw_coolant_in
  !! inlet temperature of first wall coolant [K]

  real(dp) :: temp_fw_coolant_out
  !! outlet temperature of first wall coolant [K]

  real(dp) :: pres_fw_coolant
  !! first wall coolant pressure [Pa] (`i_thermal_electric_conversion>1`)

  real(dp) :: temp_fw_peak
  !! peak first wall temperature [K]

  real(dp) :: roughness_fw_channel
  !! first wall channel roughness epsilon [m]

  real(dp) :: len_fw_channel
  !! Length of a single first wall channel (all in parallel) [m]
  !! (`iteration variable 114`, useful for `constraint equation 39`)

  real(dp) :: f_fw_peak
  !! peaking factor for first wall heat loads. (Applied separately to inboard and outboard loads.
  !! Applies to both neutron and surface loads. Only used to calculate peak temperature - not
  !! the coolant flow rate.)

  real(dp) :: pres_blkt_coolant
  !! blanket coolant pressure [Pa] (`i_thermal_electric_conversion>1`)

  real(dp) :: temp_blkt_coolant_in
  !! inlet temperature of blanket coolant  [K] (`i_thermal_electric_conversion>1`)

  real(dp) :: temp_blkt_coolant_out
  !! Outlet temperature of blanket coolant [K] (`i_thermal_electric_conversion>1`)
  !!
  !! - input if `i_blkt_coolant_type=1` (helium)
  !! - calculated if `i_blkt_coolant_type=2` (water)

  real(dp) :: coolp
  !! blanket coolant pressure [Pa] (stellarator only)

  integer :: n_blkt_outboard_modules_poloidal
  !! number of outboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)

  integer :: n_blkt_inboard_modules_poloidal
  !! number of inboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)

  integer :: n_blkt_outboard_modules_toroidal
  !! number of outboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)

  integer :: n_blkt_inboard_modules_toroidal
  !! number of inboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)

  real(dp) :: temp_fw_max
  !! maximum temperature of first wall material [K] (`i_thermal_electric_conversion>1`)

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
  !! first wall coolant fraction (calculated if `i_pulsed_plant=1` or `ipowerflow=1`)

  real(dp) :: p_div_rad_total_mw
  !! Total radiation power incident on the divertor(s) (MW)

  real(dp) :: p_fw_rad_total_mw
  !! Radiation power incident on the first wall (MW)

  real(dp) :: p_fw_hcd_rad_total_mw
  !! Radiation power incident on the heating and current drive systems on the first wall (MW)

  real(dp) :: pradloss
  !! Radiation power lost through holes (eventually hits shield) (MW)
  !! Only used for stellarator

  real(dp) :: p_tf_nuclear_heat_mw
  !! nuclear heating in the TF coil (MW)

  real(dp) :: ptfnucpm3
  !! nuclear heating in the TF coil (MW/m3) (`blktmodel>0`)
  !#TODO: check usage of old blktmodel. Update to i_blanket_type

  real(dp) :: r_cryostat_inboard
  !! cryostat radius [m]

  real(dp) :: z_cryostat_half_inside
  !! cryostat height [m]

  real(dp) :: dr_pf_cryostat
  !! Radial distance between outer edge of furthest away PF coil (or stellarator
  !! modular coil) and cryostat [m]

  real(dp) :: vol_cryostat
  !! Cryostat structure volume [m^3]

  real(dp) :: vol_cryostat_internal
  !! Internal volume of the cryostat [m^3]

  real(dp) :: vol_vv
  !! vacuum vessel volume [m^3]

  real(dp) :: vfshld
  !! coolant void fraction in shield

  real(dp) :: vol_blkt_total
  !! volume of blanket [m^3]

  real(dp) :: vol_blkt_inboard
  !! volume of inboard blanket [m^3]

  real(dp) :: vol_blkt_outboard
  !! volume of outboard blanket [m^3]

  real(dp) :: vol_shld_total
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

  real(dp) :: m_blkt_li2o
  !! mass of blanket - Li_2O part [kg]

  real(dp) :: wtbllipb
  !! mass of blanket - Li-Pb part [kg]

  real(dp) :: m_blkt_vanadium
  !! mass of blanket - vanadium part [kg]

  real(dp) :: m_blkt_lithium
  !! mass of blanket - lithium part [kg]

  real(dp) :: vfblkt
  !! coolant void fraction in blanket.

  integer :: blktmodel
  !! switch for blanket/tritium breeding model (see i_blanket_type):
  !!
  !! - =0 original simple model
  !! - =1 KIT model based on a helium-cooled pebble-bed blanket (HCPB) reference design
  !#TODO: this needs investigating and removing after any required functionality is in i_blanket_type

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
  !#TODO: this needs to be merged into i_blanket_type and then removed.

  real(dp) :: etaiso
  !! isentropic efficiency of FW and blanket coolant pumps

  real(dp) :: eta_coolant_pump_electric
  !! electrical efficiency of primary coolant pumps

  !! -----------------------------------------------------
  !! BLANKET REFACTOR
  !! For DCLL, but to be used by all mods that share blanket library after testing.
  !! Thermodynamic Model for i_coolant_pumping == 2
  !! -----------------------------------------------------

  integer :: i_fw_blkt_shared_coolant
  !! Switch for whether the FW and BB are on the same pump system
  !! i.e. do they have the same primary coolant or not
  !!  - =0    FW and BB have the same primary coolant, flow = FWin->FWout->BBin->BBout
  !!  - =1    FW and BB have the different primary coolant and are on different pump systems

  integer :: i_blkt_liquid_breeder_type
  !! Switch for Liquid Metal Breeder Material
  !!  - =0   PbLi
  !!  - =1   Li

  integer :: i_blkt_dual_coolant
  !! Switch to specify whether breeding blanket is single-cooled or dual-coolant.
  !!  - =0    Single coolant used for FW and Blanket (H2O or He). Solid Breeder.
  !!  - =1    Single coolant used for FW and Blanket (H2O or He). Liquid metal breeder
  !!          circulted for tritium extraction.
  !!  - =2    Dual coolant: primary coolant (H2O or He) for FW and blanket structure;
  !!          secondary coolant is self-cooled liquid metal breeder.

  integer :: i_blkt_liquid_breeder_channel_type
  !! Switch for Flow Channel Insert (FCI) type if liquid metal breeder blanket.
  !!  - =0    Thin conducting walls, default electrical conductivity (bz_channel_conduct_liq) is Eurofer
  !!  - =1    Insulating Material, assumed perfect electrical insulator, default density (den_ceramic) is for SiC
  !!  - =2    Insulating Material, electrical conductivity (bz_channel_conduct_liq) is input (default Eurofer), default density (den_ceramic) is for SiC

  integer :: i_blkt_module_segmentation
  !! Switch for Multi Module Segment (MMS) or Single Modle Segment (SMS)
  !!  - =0    MMS
  !!  - =1    SMS

  integer :: n_liq_recirc
  !! Number of liquid metal breeder recirculations per day, for use with i_blkt_dual_coolant=1

  real(dp) :: r_f_liq_ib, r_f_liq_ob
  !! Radial fraction of BZ liquid channels

  real(dp) :: w_f_liq_ib, w_f_liq_ob
  !! Toroidal fraction of BZ liquid channels

  real(dp) :: den_ceramic
  !! FCI material density

  real(dp) :: th_wall_secondary
  !! Liquid metal coolant/breeder wall thickness thin conductor or FCI [m]

  real(dp) :: bz_channel_conduct_liq
  !! Liquid metal coolant/breeder thin conductor or FCI wall conductance [A V^-1 m^-1]

  real(dp) :: a_bz_liq
  !! Toroidal width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone

  real(dp) :: b_bz_liq
  !! Radial width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone

  integer :: nopol
  !! Number of poloidal sections in a liquid metal breeder/coolant channel for module/segment

  integer :: nopipes
  !! Number of Liquid metal breeder/coolant channels per module/segment

  real(dp) :: den_liq
  !! Liquid metal breeder/coolant density [kg m^-3]

  real(dp) ::wht_liq, wht_liq_ib, wht_liq_ob
  !! Liquid metal

  real(dp) :: specific_heat_liq
  !! Liquid metal breeder/coolant specific heat [J kg^-1 K^-1]

  real(dp) :: thermal_conductivity_liq
  !! Liquid metal breeder/coolant thermal conductivity [W m^-1 K^-1]

  real(dp) :: dynamic_viscosity_liq
  !! Liquid metal breeder/coolant dynamic viscosity [Pa s]

  real(dp) :: electrical_conductivity_liq
  !! Liquid metal breeder/coolant electrical conductivity [Ohm m]

  real(dp) :: hartmann_liq(2)
  !! Hartmann number

  real(dp) :: b_mag_blkt(2)
  !! Toroidal Magnetic feild strength for IB/OB blanket [T]

  real(dp) :: etaiso_liq
  !! Isentropic efficiency of blanket liquid breeder/coolant pumps

  real(dp) :: blpressure_liq
  !! blanket liquid metal breeder/coolant pressure [Pa]

  real(dp) :: inlet_temp_liq, outlet_temp_liq
  !! Inlet (scan var 68) and Outlet (scan var 69) temperature of the liquid breeder/coolant [K]

  real(dp) :: den_fw_coolant
  !! Density of the FW primary coolant

  real(dp) :: visc_fw_coolant
  !! Viscosity of the FW primary coolant

  real(dp) :: den_blkt_coolant
  !! Density of the blanket primary coolant

  real(dp) :: visc_blkt_coolant
  !! Viscosity of the blanket primary coolant

  real(dp) :: cp_fw, cv_fw, cp_bl, cv_bl
  !! Spesific heat for FW and blanket primary coolant(s)

  real(dp) :: f_nuc_pow_bz_struct
  !! For a dual-coolant blanket, fraction of BZ power cooled by primary coolant

  real(dp) :: f_nuc_pow_bz_liq
  !! For a dual-coolant blanket, fraction of BZ self-cooled power (secondary coolant)

  real(dp) :: pnuc_fw_ratio_dcll, pnuc_blkt_ratio_dcll
  !! For a dual-coolant blanket, ratio of FW/Blanket nuclear power as fraction of total

  integer :: bzfllengi_n_rad, bzfllengi_n_pol, bzfllengo_n_rad, bzfllengo_n_pol
  !! Number of radial and poloidal sections that make up the total primary coolant flow
  !! length in a blanket module (IB and OB)

  integer :: bzfllengi_n_rad_liq, bzfllengi_n_pol_liq, bzfllengo_n_rad_liq, bzfllengo_n_pol_liq
  !! Number of radial and poloidal sections that make up the total secondary coolant/breeder
  !! flow length in a blanket module (IB and OB)
end module fwbs_variables
