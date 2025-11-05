"""author: J. Morris (UKAEA), M. Kovari (UKAEA)
Module containing global variables relating to the first wall, blanket and
shield components
### References
-
"""

life_blkt_fpy: float = None
"""Full power blanket lifetime (years)"""


life_blkt: float = None
"""Calendar year blanket lifetime (years)"""


m_fw_blkt_div_coolant_total: float = None
"""mass of water coolant (in shield, blanket, first wall, divertor) [kg]"""


m_vv: float = None
"""vacuum vessel mass [kg]"""


den_steel: float = None
"""density of steel [kg m^-3]"""


denwc: float = None
"""density of tungsten carbide [kg m^-3]"""


dewmkg: float = None
"""total mass of vacuum vessel + cryostat [kg] (calculated if blktmodel>0)"""


f_p_blkt_multiplication: float = None
"""energy multiplication in blanket and shield"""


p_blkt_multiplication_mw: float = None
"""power due to energy multiplication in blanket and shield [MW]"""


fblss: float = None
"""KIT blanket model: steel fraction of breeding zone"""


f_ster_div_single: float = None
"""Solid angle fraction taken by one divertor"""


f_a_fw_outboard_hcd: float = None
"""area fraction of first wall covered by heating/current drive apparatus plus diagnostics"""


fhole: float = None
"""area fraction taken up by other holes (IFE)"""


i_fw_blkt_vv_shape: int = None
"""switch for first wall, blanket, shield and vacuum vessel shape:
- =1 D-shaped (cylinder inboard + ellipse outboard)
- =2 defined by two ellipses
"""


life_fw_fpy: float = None
"""first wall full-power year lifetime (y)"""


m_fw_total: float = None
"""first wall mass [kg]"""


fw_armour_mass: float = None
"""first wall armour mass [kg]"""


fw_armour_thickness: float = None
"""first wall armour thickness [m]"""


fw_armour_vol: float = None
"""first wall armour volume [m^3]"""


i_blanket_type: int = None
"""switch for blanket model:
- =1 CCFE HCPB model
- =2 KIT HCPB model  # REMOVED, no longer usable
- =3 CCFE HCPB model with Tritium Breeding Ratio calculation
- =4 KIT HCLL model  # REMOVED, no longer usable
- =5 DCLL model -  no nutronics model included (in development) please check/choose values for
'dual-coolant blanket' fractions (provided in this file).
-  please use i_p_coolant_pumping = 0 or 1.
"""


i_blkt_inboard: int = None
"""switch for inboard blanket:
- =0 No inboard blanket (dr_blkt_inboard=0.0)
- =1 Inboard blanket present
"""


inuclear: int = None
"""switch for nuclear heating in the coils:
- =0 Frances Fox model (default)
- =1 Fixed by user (qnuc)
"""


qnuc: float = None
"""nuclear heating in the coils (W) (`inuclear=1`)"""


f_blkt_li6_enrichment: float = None
"""lithium-6 enrichment of breeding material (%)"""


p_blkt_nuclear_heat_total_mw: float = None
"""nuclear heating in the blanket [MW]"""


pnuc_cp: float = None
"""Total nuclear heating in the ST centrepost [MW]"""


p_cp_shield_nuclear_heat_mw: float = None
"""Neutronic shield nuclear heating in the ST centrepost [MW]"""


pnuc_cp_tf: float = None
"""TF neutronic nuclear heating in the ST centrepost [MW]"""


p_div_nuclear_heat_total_mw: float = None
"""nuclear heating in the divertor [MW]"""


p_fw_nuclear_heat_total_mw: float = None
"""nuclear heating in the first wall [MW]"""


p_fw_hcd_nuclear_heat_mw: float = None
"""Nuclear heating in the HCD apparatus and diagnostics on the first wall [MW]"""


pnucloss: float = None
"""nuclear heating lost via holes [MW]"""


pnucvvplus: float = None
"""nuclear heating to vacuum vessel and beyond [MW]"""


p_shld_nuclear_heat_mw: float = None
"""nuclear heating in the shield [MW]"""


m_blkt_total: float = None
"""mass of blanket [kg]"""


m_blkt_steel_total: float = None
"""mass of blanket - steel part [kg]"""


armour_fw_bl_mass: float = None
"""Total mass of armour, first wall and blanket [kg]"""


# CCFE HCPB Blanket Model i_blanket_type=1


breeder_f: float = None
"""Volume ratio: Li4SiO4/(Be12Ti+Li4SiO4) (`iteration variable 108`)"""


breeder_multiplier: float = None
"""combined breeder/multipler fraction of blanket by volume"""


vfcblkt: float = None
"""He coolant fraction of blanket by volume (`i_blanket_type= 1,3` (CCFE HCPB))"""


vfpblkt: float = None
"""He purge gas fraction of blanket by volume (`i_blanket_type= 1,3` (CCFE HCPB))"""


m_blkt_li4sio4: float = None
"""mass of lithium orthosilicate in blanket [kg] (`i_blanket_type=1,3` (CCFE HCPB))"""


m_blkt_tibe12: float = None
"""mass of titanium beryllide in blanket [kg] (`i_blanket_type=1,3` (CCFE HCPB))"""


neut_flux_cp: float = None
"""Centrepost TF fast neutron flux (E > 0.1 MeV) [m^(-2).^(-1)]
This variable is only calculated for superconducting (i_tf_sup = 1 )
spherical tokamal magnet designs (itart = 0)
"""


f_neut_shield: float = None
"""Fraction of nuclear power shielded before the CP magnet (ST)
( neut_absorb = -1 --> a fit on simplified MCNP neutronic
calculation is used assuming water cooled (13%) tungesten carbyde )
"""


f_a_fw_coolant_inboard: float = None
"""Inboard FW coolant cross-sectional area void fraction"""


f_a_fw_coolant_outboard: float = None
"""Outboard FW coolant cross-sectional area void fraction"""


psurffwi: float = None
"""Surface heat flux on first wall [MW] (sum = p_fw_rad_total_mw)"""


psurffwo: float = None
"""Surface heat flux on first wall [MW] (sum = p_fw_rad_total_mw)"""


vol_fw_total: float = None
"""First wall volume [m3]"""


f_vol_blkt_steel: float = None
"""Fractions of blanket by volume: steel"""


f_vol_blkt_li4sio4: float = None
"""Fractions of blanket by volume: lithium orthosilicate"""


f_vol_blkt_tibe12: float = None
"""Fractions of blanket by volume: titanium beryllide"""


breedmat: int = None
"""breeder material switch (i_blanket_type=2 (KIT HCPB)):
- =1 Lithium orthosilicate
- =2 Lithium methatitanate
- =3 Lithium zirconate
"""


densbreed: float = None
"""density of breeder material [kg m^-3] (`i_blanket_type=2` (KIT HCPB))"""


fblbe: float = None
"""beryllium fraction of blanket by volume (if `i_blanket_type=2`, is Be fraction of breeding zone)"""


fblbreed: float = None
"""breeder fraction of blanket breeding zone by volume (`i_blanket_type=2` (KIT HCPB))"""


fblhebmi: float = None
"""helium fraction of inboard blanket box manifold by volume (`i_blanket_type=2` (KIT HCPB))"""


fblhebmo: float = None
"""helium fraction of outboard blanket box manifold by volume (`i_blanket_type=2` (KIT HCPB))"""


fblhebpi: float = None
"""helium fraction of inboard blanket back plate by volume (`i_blanket_type=2` (KIT HCPB))"""


fblhebpo: float = None
"""helium fraction of outboard blanket back plate by volume (`i_blanket_type=2` (KIT HCPB))"""


hcdportsize: int = None
"""switch for size of heating/current drive ports (`i_blanket_type=2` (KIT HCPB)):
- =1 'small'
- =2 'large'
"""


nflutf: float = None
"""peak fast neutron fluence on TF coil superconductor [n m^-2] (`i_blanket_type=2` (KIT HCPB))"""


npdiv: int = None
"""number of divertor ports (`i_blanket_type=2` (KIT HCPB))"""


nphcdin: int = None
"""number of inboard ports for heating/current drive (`i_blanket_type=2` (KIT HCPB))"""


nphcdout: int = None
"""number of outboard ports for heating/current drive (`i_blanket_type=2` (KIT HCPB))"""


tbr: float = None
"""tritium breeding ratio (`i_blanket_type=2,3` (KIT HCPB/HCLL))"""


tritprate: float = None
"""tritium production rate [g day^-1] (`i_blanket_type=2` (KIT HCPB))"""


wallpf: float = None
"""neutron wall load peaking factor (`i_blanket_type=2` (KIT HCPB))"""


whtblbreed: float = None
"""mass of blanket - breeder part [kg] (`i_blanket_type=2` (KIT HCPB))"""


m_blkt_beryllium: float = None
"""mass of blanket - beryllium part [kg]"""


i_p_coolant_pumping: int = None
"""Switch for pumping power for primary coolant (mechanical power only and peak first wall
temperature is only calculated if `i_p_coolant_pumping=2`):
- =0 User sets pump power directly (p_blkt_coolant_pump_mw, p_fw_coolant_pump_mw, p_div_coolant_pump_mw, p_shld_coolant_pump_mw)
- =1 User sets pump power as a fraction of thermal power (f_p_blkt_coolant_pump_total_heat, f_p_fw_coolant_pump_total_heat, f_p_div_coolant_pump_total_heat, f_p_shld_coolant_pump_total_heat)
- =2 Mechanical pumping power is calculated
- =3 Mechanical pumping power is calculated using specified pressure drop
"""


i_shield_mat: int = None
"""Switch for shield material - *currently only applied in costing routines* `cost_model = 2`
- =0 Tungsten (default)
- =1 Tungsten carbide
"""


i_thermal_electric_conversion: int = None
"""Switch for power conversion cycle:
- =0 Set efficiency for chosen blanket, from detailed models (divertor heat not used)
- =1 Set efficiency for chosen blanket, from detailed models (divertor heat used)
- =2 user input thermal-electric efficiency (eta_turbine)
- =3 steam Rankine cycle
- =4 supercritical CO2 cycle
"""


secondary_cycle_liq: int = None
"""Switch for power conversion cycle for the liquid breeder component of the blanket:
- =2 user input thermal-electric efficiency (eta_turbine)
- =4 supercritical CO2 cycle
"""


i_blkt_coolant_type: int = None
"""Switch for blanket coolant (set via blkttype):
- =1 helium
- =2 pressurized water
"""


i_fw_coolant_type: str = None
"""switch for first wall coolant (can be different from blanket coolant):
- 'helium'
- 'water'
"""

i_fw_inboard_coolant_type: str = None
"""switch for inboard first wall coolant (can be different from blanket coolant):
- 'helium'
- 'water'
"""


dr_fw_wall: float = None
"""wall thickness of first wall coolant channels [m]"""

dr_fw_inboard_wall: float = None
"""wall thickness of inboard first wall coolant channels [m]"""


radius_fw_channel: float = None
"""radius of first wall cooling channels [m]"""

radius_fw_inboard_channel: float = None
"""radius of inboard first wall cooling channels [m]"""


dx_fw_module: float = None
"""Width of a FW module containing a cooling channel [m]"""

dx_fw_inboard_module: float = None
"""Width of the inboard FW module containing a cooling channel [m]"""


temp_fw_coolant_in: float = None
"""inlet temperature of first wall coolant [K]"""

temp_fw_inboard_coolant_in: float = None
"""inlet temperature of inboard first wall coolant [K]"""


temp_fw_coolant_out: float = None
"""outlet temperature of first wall coolant [K]"""

temp_fw_inboard_coolant_out: float = None
"""outlet temperature of inboard first wall coolant [K]"""

pres_fw_coolant: float = None
"""first wall coolant pressure [Pa] (`i_thermal_electric_conversion>1`)"""

pres_fw_inboard_coolant: float = None
"""inboard first wall coolant pressure [Pa] (`i_thermal_electric_conversion>1`)"""


temp_fw_peak: float = None
"""peak first wall temperature [K]"""

temp_fw_inboard_peak: float = None
"""peak inboard first wall temperature [K]"""

roughness_fw_channel: float = None
"""first wall channel roughness epsilon [m]"""

roughness_fw_inboard_channel: float = None
"""inboard first wall channel roughness epsilon [m]"""

len_fw_channel: float = None
"""Length of a single first wall channel (all in parallel) [m]
(`iteration variable 114`, useful for `constraint equation 39`)
"""

len_fw_inboard_channel: float = None
"""Length of a single inboard first wall channel (all in parallel) [m]
(`iteration variable 115`, useful for `constraint equation 40`)
"""


f_fw_peak: float = None
"""peaking factor for first wall heat loads. (Applied separately to inboard and outboard loads.
Applies to both neutron and surface loads. Only used to calculate peak temperature - not
the coolant flow rate.)
"""


pres_blkt_coolant: float = None
"""blanket coolant pressure [Pa] (`i_thermal_electric_conversion>1`)"""


temp_blkt_coolant_in: float = None
"""inlet temperature of blanket coolant  [K] (`i_thermal_electric_conversion>1`)"""


temp_blkt_coolant_out: float = None
"""Outlet temperature of blanket coolant [K] (`i_thermal_electric_conversion>1`)
- input if `i_blkt_coolant_type=1` (helium)
- calculated if `i_blkt_coolant_type=2` (water)
"""


coolp: float = None
"""blanket coolant pressure [Pa] (stellarator only)"""


n_blkt_outboard_modules_poloidal: int = None
"""number of outboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)"""


n_blkt_inboard_modules_poloidal: int = None
"""number of inboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)"""


n_blkt_outboard_modules_toroidal: int = None
"""number of outboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)"""


n_blkt_inboard_modules_toroidal: int = None
"""number of inboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)"""


temp_fw_max: float = None
"""maximum temperature of first wall material [K] (`i_thermal_electric_conversion>1`)"""


fw_th_conductivity: float = None
"""thermal conductivity of first wall material at 293 K (W/m/K) (Temperature dependence
is as for unirradiated Eurofer)
"""


fvoldw: float = None
"""area coverage factor for vacuum vessel volume"""


fvolsi: float = None
"""area coverage factor for inboard shield volume"""


fvolso: float = None
"""area coverage factor for outboard shield volume"""


fwclfr: float = None
"""first wall coolant fraction (calculated if `i_pulsed_plant=1` or `ipowerflow=1`)"""


p_div_rad_total_mw: float = None
"""Total radiation power incident on the divertor(s) (MW)"""


p_fw_rad_total_mw: float = None
"""Radiation power incident on the first wall (MW)"""


p_fw_hcd_rad_total_mw: float = None
"""Radiation power incident on the heating and current drive systems on the first wall (MW)"""


pradloss: float = None
"""Radiation power lost through holes (eventually hits shield) (MW)
Only used for stellarator
"""


p_tf_nuclear_heat_mw: float = None
"""nuclear heating in the TF coil (MW)"""


ptfnucpm3: float = None
"""nuclear heating in the TF coil (MW/m3) (`blktmodel>0`)"""


r_cryostat_inboard: float = None
"""cryostat radius [m]"""


z_cryostat_half_inside: float = None
"""cryostat height [m]"""


dr_pf_cryostat: float = None
"""Radial distance between outer edge of furthest away PF coil (or stellarator
modular coil) and cryostat [m]
"""


vol_cryostat: float = None
"""Cryostat structure volume [m^3]"""


vol_cryostat_internal: float = None
"""Internal volume of the cryostat [m^3]"""


vol_vv: float = None
"""vacuum vessel volume [m^3]"""


vfshld: float = None
"""coolant void fraction in shield"""


vol_blkt_total: float = None
"""volume of blanket [m^3]"""


vol_blkt_inboard: float = None
"""volume of inboard blanket [m^3]"""


vol_blkt_outboard: float = None
"""volume of outboard blanket [m^3]"""


vol_shld_total: float = None
"""volume of shield [m^3]"""


whtshld: float = None
"""mass of shield [kg]"""


wpenshld: float = None
"""mass of the penetration shield [kg]"""


wtshldi: float = None
"""mass of inboard shield [kg]"""


wtshldo: float = None
"""mass of outboard shield [kg]"""


irefprop: int = None
"""Switch to use REFPROP routines (stellarator only)"""


fblli: float = None
"""lithium fraction of blanket by volume (stellarator only)"""


fblli2o: float = None
"""lithium oxide fraction of blanket by volume (stellarator only)"""


fbllipb: float = None
"""lithium lead fraction of blanket by volume (stellarator only)"""


fblvd: float = None
"""vanadium fraction of blanket by volume (stellarator only)"""


m_blkt_li2o: float = None
"""mass of blanket - Li_2O part [kg]"""


wtbllipb: float = None
"""mass of blanket - Li-Pb part [kg]"""


m_blkt_vanadium: float = None
"""mass of blanket - vanadium part [kg]"""


m_blkt_lithium: float = None
"""mass of blanket - lithium part [kg]"""


f_a_blkt_cooling_channels: float = None
"""coolant void fraction in blanket."""


blktmodel: int = None
"""switch for blanket/tritium breeding model (see i_blanket_type):
- =0 original simple model
- =1 KIT model based on a helium-cooled pebble-bed blanket (HCPB) reference design
"""


declblkt: float = None
"""neutron power deposition decay length of blanket structural material [m] (stellarators only)"""


declfw: float = None
"""neutron power deposition decay length of first wall structural material [m] (stellarators only)"""


declshld: float = None
"""neutron power deposition decay length of shield structural material [m] (stellarators only)"""


blkttype: int = None
"""Switch for blanket type:
- =1 WCLL;
- =2 HCLL; efficiency taken from M. Kovari 2016
"PROCESS": A systems code for fusion power plants - Part 2: Engineering
https://www.sciencedirect.com/science/article/pii/S0920379616300072
Feedheat & reheat cycle assumed
- =3 HCPB; efficiency taken from M. Kovari 2016
"PROCESS": A systems code for fusion power plants - Part 2: Engineering
https://www.sciencedirect.com/science/article/pii/S0920379616300072
Feedheat & reheat cycle assumed
"""


etaiso: float = None
"""isentropic efficiency of FW and blanket coolant pumps"""


eta_coolant_pump_electric: float = None
"""electrical efficiency of primary coolant pumps"""


i_fw_blkt_shared_coolant: int = None
"""Switch for whether the FW and BB are on the same pump system
i.e. do they have the same primary coolant or not
- =0    FW and BB have the same primary coolant, flow = FWin->FWout->BBin->BBout
- =1    FW and BB have the different primary coolant and are on different pump systems
"""


i_blkt_liquid_breeder_type: int = None
"""Switch for Liquid Metal Breeder Material
- =0   PbLi
- =1   Li
"""


i_blkt_dual_coolant: int = None
"""Switch to specify whether breeding blanket is single-cooled or dual-coolant.
- =0    Single coolant used for FW and Blanket (H2O or He). Solid Breeder.
- =1    Single coolant used for FW and Blanket (H2O or He). Liquid metal breeder
circulted for tritium extraction.
- =2    Dual coolant: primary coolant (H2O or He) for FW and blanket structure;
secondary coolant is self-cooled liquid metal breeder.
"""


i_blkt_liquid_breeder_channel_type: int = None
"""Switch for Flow Channel Insert (FCI) type if liquid metal breeder blanket.
- =0    Thin conducting walls, default electrical conductivity (bz_channel_conduct_liq) is Eurofer
- =1    Insulating Material, assumed perfect electrical insulator, default density (den_ceramic) is for SiC
- =2    Insulating Material, electrical conductivity (bz_channel_conduct_liq) is input (default Eurofer), default density (den_ceramic) is for SiC
"""


i_blkt_module_segmentation: int = None
"""Switch for Multi Module Segment (MMS) or Single Modle Segment (SMS)
- =0    MMS
- =1    SMS
"""


n_liq_recirc: int = None
"""Number of liquid metal breeder recirculations per day, for use with i_blkt_dual_coolant=1"""


r_f_liq_ib: float = None
"""Radial fraction of BZ liquid channels"""


r_f_liq_ob: float = None
"""Radial fraction of BZ liquid channels"""


w_f_liq_ib: float = None
"""Toroidal fraction of BZ liquid channels"""


w_f_liq_ob: float = None
"""Toroidal fraction of BZ liquid channels"""


den_ceramic: float = None
"""FCI material density"""


th_wall_secondary: float = None
"""Liquid metal coolant/breeder wall thickness thin conductor or FCI [m]"""


bz_channel_conduct_liq: float = None
"""Liquid metal coolant/breeder thin conductor or FCI wall conductance [A V^-1 m^-1]"""


a_bz_liq: float = None
"""Toroidal width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone"""


b_bz_liq: float = None
"""Radial width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone"""


nopol: int = None
"""Number of poloidal sections in a liquid metal breeder/coolant channel for module/segment"""


nopipes: int = None
"""Number of Liquid metal breeder/coolant channels per module/segment"""


den_liq: float = None
"""Liquid metal breeder/coolant density [kg m^-3]"""


wht_liq: float = None
"""Liquid metal"""


wht_liq_ib: float = None
"""Liquid metal"""


wht_liq_ob: float = None
"""Liquid metal"""


specific_heat_liq: float = None
"""Liquid metal breeder/coolant specific heat [J kg^-1 K^-1]"""


thermal_conductivity_liq: float = None
"""Liquid metal breeder/coolant thermal conductivity [W m^-1 K^-1]"""


dynamic_viscosity_liq: float = None
"""Liquid metal breeder/coolant dynamic viscosity [Pa s]"""


electrical_conductivity_liq: float = None
"""Liquid metal breeder/coolant electrical conductivity [Ohm m]"""


hartmann_liq: list[float] = None
"""Hartmann number"""


b_mag_blkt: list[float] = None
"""Toroidal Magnetic feild strength for IB/OB blanket [T]"""


etaiso_liq: float = None
"""Isentropic efficiency of blanket liquid breeder/coolant pumps"""


blpressure_liq: float = None
"""blanket liquid metal breeder/coolant pressure [Pa]"""


inlet_temp_liq: float = None
"""Inlet (scan var 68) temperature of the liquid breeder/coolant [K]"""


outlet_temp_liq: float = None
"""Outlet (scan var 69) temperature of the liquid breeder/coolant [K]"""


den_fw_coolant: float = None
"""Density of the FW primary coolant"""


visc_fw_coolant: float = None
"""Viscosity of the FW primary coolant"""


den_blkt_coolant: float = None
"""Density of the blanket primary coolant"""


visc_blkt_coolant: float = None
"""Viscosity of the blanket primary coolant"""


cp_fw: float = None
"""Spesific heat for FW and blanket primary coolant(s)"""


cv_fw: float = None
"""Spesific heat for FW and blanket primary coolant(s)"""


cp_bl: float = None
"""Spesific heat for FW and blanket primary coolant(s)"""


cv_bl: float = None
"""Spesific heat for FW and blanket primary coolant(s)"""


f_nuc_pow_bz_struct: float = None
"""For a dual-coolant blanket, fraction of BZ power cooled by primary coolant"""


f_nuc_pow_bz_liq: float = None
"""For a dual-coolant blanket, fraction of BZ self-cooled power (secondary coolant)"""


pnuc_fw_ratio_dcll: float = None
"""For a dual-coolant blanket, ratio of FW nuclear power as fraction of total"""


pnuc_blkt_ratio_dcll: float = None
"""For a dual-coolant blanket, ratio of Blanket nuclear power as fraction of total"""


n_blkt_inboard_module_coolant_sections_radial: int = None
"""Number of radial and poloidal sections that make up the total primary coolant flow
length in a blanket module (IB and OB)
"""


n_blkt_inboard_module_coolant_sections_poloidal: int = None
"""Number of radial and poloidal sections that make up the total primary coolant flow
length in a blanket module (IB and OB)
"""


n_blkt_outboard_module_coolant_sections_radial: int = None
"""Number of radial and poloidal sections that make up the total primary coolant flow
length in a blanket module (IB and OB)
"""


n_blkt_outboard_module_coolant_sections_poloidal: int = None
"""Number of radial and poloidal sections that make up the total primary coolant flow
length in a blanket module (IB and OB)
"""

bzfllengi_n_rad_liq: int = None
"""Number of radial and poloidal sections that make up the total secondary coolant/breeder
flow length in a blanket module (IB and OB)
"""


bzfllengi_n_pol_liq: int = None
"""Number of radial and poloidal sections that make up the total secondary coolant/breeder
flow length in a blanket module (IB and OB)
"""


bzfllengo_n_rad_liq: int = None
"""Number of radial and poloidal sections that make up the total secondary coolant/breeder
flow length in a blanket module (IB and OB)
"""


bzfllengo_n_pol_liq: int = None
"""Number of radial and poloidal sections that make up the total secondary coolant/breeder
flow length in a blanket module (IB and OB)
"""

radius_blkt_channel: float = None
"""Radius of blanket cooling channels [m]"""

radius_blkt_channel_90_bend: float = None
"""Radius of blanket cooling channel 90° bend [m]"""

radius_fw_channel_90_bend: float = None
"""Radius of first wall cooling channel 90° bend [m]"""

radius_fw_channel_180_bend: float = None
"""Radius of first wall cooling channel 180° bend [m]"""

radius_blkt_channel_180_bend: float = None
"""Radius of blanket cooling channel 180° bend [m]"""


def init_fwbs_variables():
    """Initialise FWBS variables"""
    global life_blkt_fpy
    global life_blkt
    global m_fw_blkt_div_coolant_total
    global m_vv
    global den_steel
    global denwc
    global dewmkg
    global f_p_blkt_multiplication
    global p_blkt_multiplication_mw
    global fblss
    global f_ster_div_single
    global f_a_fw_outboard_hcd
    global fhole
    global i_fw_blkt_vv_shape
    global life_fw_fpy
    global m_fw_total
    global fw_armour_mass
    global fw_armour_thickness
    global fw_armour_vol
    global i_blanket_type
    global i_blkt_inboard
    global inuclear
    global qnuc
    global f_blkt_li6_enrichment
    global p_blkt_nuclear_heat_total_mw
    global pnuc_cp
    global p_cp_shield_nuclear_heat_mw
    global pnuc_cp_tf
    global p_div_nuclear_heat_total_mw
    global p_fw_nuclear_heat_total_mw
    global p_fw_hcd_nuclear_heat_mw
    global pnucloss
    global pnucvvplus
    global p_shld_nuclear_heat_mw
    global m_blkt_total
    global m_blkt_steel_total
    global armour_fw_bl_mass
    global breeder_f
    global breeder_multiplier
    global vfcblkt
    global vfpblkt
    global m_blkt_li4sio4
    global m_blkt_tibe12
    global neut_flux_cp
    global f_neut_shield
    global f_a_fw_coolant_inboard
    global f_a_fw_coolant_outboard
    global psurffwi
    global psurffwo
    global vol_fw_total
    global f_vol_blkt_steel
    global f_vol_blkt_li4sio4
    global f_vol_blkt_tibe12
    global breedmat
    global densbreed
    global fblbe
    global fblbreed
    global fblhebmi
    global fblhebmo
    global fblhebpi
    global fblhebpo
    global hcdportsize
    global nflutf
    global npdiv
    global nphcdin
    global nphcdout
    global tbr
    global tritprate
    global wallpf
    global whtblbreed
    global m_blkt_beryllium
    global i_p_coolant_pumping
    global i_shield_mat
    global i_thermal_electric_conversion
    global secondary_cycle_liq
    global i_blkt_coolant_type
    global i_fw_coolant_type
    global i_fw_inboard_coolant_type
    global dr_fw_wall
    global dr_fw_inboard_wall
    global radius_fw_channel
    global radius_fw_inboard_channel
    global dx_fw_module
    global dx_fw_inboard_module
    global temp_fw_coolant_in
    global temp_fw_inboard_coolant_in
    global temp_fw_coolant_out
    global temp_fw_inboard_coolant_out
    global pres_fw_coolant
    global pres_fw_inboard_coolant
    global temp_fw_peak
    global temp_fw_inboard_peak
    global roughness_fw_channel
    global roughness_fw_inboard_channel
    global len_fw_channel
    global len_fw_inboard_channel
    global f_fw_peak
    global pres_blkt_coolant
    global temp_blkt_coolant_in
    global temp_blkt_coolant_out
    global coolp
    global n_blkt_outboard_modules_poloidal
    global n_blkt_inboard_modules_poloidal
    global n_blkt_outboard_modules_toroidal
    global n_blkt_inboard_modules_toroidal
    global temp_fw_max
    global fw_th_conductivity
    global fvoldw
    global fvolsi
    global fvolso
    global fwclfr
    global p_div_rad_total_mw
    global p_fw_rad_total_mw
    global p_fw_hcd_rad_total_mw
    global pradloss
    global p_tf_nuclear_heat_mw
    global ptfnucpm3
    global r_cryostat_inboard
    global z_cryostat_half_inside
    global dr_pf_cryostat
    global vol_cryostat
    global vol_cryostat_internal
    global vol_vv
    global vfshld
    global vol_blkt_total
    global vol_blkt_inboard
    global vol_blkt_outboard
    global vol_shld_total
    global whtshld
    global wpenshld
    global wtshldi
    global wtshldo
    global irefprop
    global fblli
    global fblli2o
    global fbllipb
    global fblvd
    global m_blkt_li2o
    global wtbllipb
    global m_blkt_vanadium
    global m_blkt_lithium
    global f_a_blkt_cooling_channels
    global blktmodel
    global declblkt
    global declfw
    global declshld
    global blkttype
    global etaiso
    global eta_coolant_pump_electric
    global i_fw_blkt_shared_coolant
    global i_blkt_liquid_breeder_type
    global i_blkt_dual_coolant
    global i_blkt_liquid_breeder_channel_type
    global i_blkt_module_segmentation
    global n_liq_recirc
    global r_f_liq_ib
    global r_f_liq_ob
    global w_f_liq_ib
    global w_f_liq_ob
    global den_ceramic
    global th_wall_secondary
    global bz_channel_conduct_liq
    global a_bz_liq
    global b_bz_liq
    global nopol
    global nopipes
    global den_liq
    global wht_liq
    global wht_liq_ib
    global wht_liq_ob
    global specific_heat_liq
    global thermal_conductivity_liq
    global dynamic_viscosity_liq
    global electrical_conductivity_liq
    global hartmann_liq
    global b_mag_blkt
    global etaiso_liq
    global blpressure_liq
    global inlet_temp_liq
    global outlet_temp_liq
    global den_fw_coolant
    global visc_fw_coolant
    global den_blkt_coolant
    global visc_blkt_coolant
    global cp_fw
    global cv_fw
    global cp_bl
    global cv_bl
    global f_nuc_pow_bz_struct
    global f_nuc_pow_bz_liq
    global pnuc_fw_ratio_dcll
    global pnuc_blkt_ratio_dcll
    global n_blkt_inboard_module_coolant_sections_radial
    global n_blkt_inboard_module_coolant_sections_poloidal
    global n_blkt_outboard_module_coolant_sections_radial
    global n_blkt_outboard_module_coolant_sections_poloidal
    global bzfllengi_n_rad_liq
    global bzfllengi_n_pol_liq
    global bzfllengo_n_rad_liq
    global bzfllengo_n_pol_liq
    global radius_blkt_channel
    global radius_fw_channel_90_bend
    global radius_fw_channel_180_bend
    global radius_blkt_channel_90_bend
    global radius_blkt_channel_180_bend

    life_blkt_fpy = 0.0
    life_blkt = 0.0
    m_fw_blkt_div_coolant_total = 0.0
    m_vv = 0.0
    den_steel = 7800.0
    denwc = 15630.0
    dewmkg = 0.0
    f_p_blkt_multiplication = 1.269
    p_blkt_multiplication_mw = 0.0
    fblss = 0.09705
    f_ster_div_single = 0.115
    f_a_fw_outboard_hcd = 0.0
    fhole = 0.0
    i_fw_blkt_vv_shape = 2
    life_fw_fpy = 0.0
    m_fw_total = 0.0
    fw_armour_mass = 0.0
    fw_armour_thickness = 0.005
    fw_armour_vol = 0.0
    i_blanket_type = 1
    i_blkt_inboard = 1
    inuclear = 0
    qnuc = 0.0
    f_blkt_li6_enrichment = 30.0
    p_blkt_nuclear_heat_total_mw = 0.0
    p_div_nuclear_heat_total_mw = 0.0
    p_fw_nuclear_heat_total_mw = 0.0
    p_fw_hcd_nuclear_heat_mw = 0.0
    pnucloss = 0.0
    pnucvvplus = 0.0
    p_shld_nuclear_heat_mw = 0.0
    m_blkt_total = 0.0
    m_blkt_steel_total = 0.0
    armour_fw_bl_mass = 0.0
    breeder_f = 0.5
    breeder_multiplier = 0.75
    vfcblkt = 0.05295
    vfpblkt = 0.1
    m_blkt_li4sio4 = 0.0
    m_blkt_tibe12 = 0.0
    f_neut_shield = -1.0
    f_a_fw_coolant_inboard = 0.0
    f_a_fw_coolant_outboard = 0.0
    psurffwi = 0.0
    psurffwo = 0.0
    vol_fw_total = 0.0
    f_vol_blkt_steel = 0.0
    f_vol_blkt_li4sio4 = 0.0
    f_vol_blkt_tibe12 = 0.0
    breedmat = 1
    densbreed = 0.0
    fblbe = 0.6
    fblbreed = 0.154
    fblhebmi = 0.4
    fblhebmo = 0.4
    fblhebpi = 0.6595
    fblhebpo = 0.6713
    hcdportsize = 1
    nflutf = 0.0
    npdiv = 2
    nphcdin = 2
    nphcdout = 2
    tbr = 0.0
    tritprate = 0.0
    wallpf = 1.21
    whtblbreed = 0.0
    m_blkt_beryllium = 0.0
    i_p_coolant_pumping = 2
    i_shield_mat = 0
    i_thermal_electric_conversion = 0
    secondary_cycle_liq = 4
    i_blkt_coolant_type = 1
    i_fw_coolant_type = "helium"
    i_fw_inboard_coolant_type = "helium"
    dr_fw_wall = 0.003
    dr_fw_inboard_wall = 0.003
    radius_fw_channel = 0.006
    radius_fw_inboard_channel = 0.006
    dx_fw_module = 0.02
    dx_fw_inboard_module = 0.02
    temp_fw_coolant_in = 573.0
    temp_fw_inboard_coolant_in = 573.0
    temp_fw_coolant_out = 823.0
    temp_fw_inboard_coolant_out = 823.0
    pres_fw_coolant = 15.5e6
    pres_fw_inboard_coolant = 15.5e6
    temp_fw_peak = 873.0
    temp_fw_inboard_peak = 873.0
    roughness_fw_channel = 1.0e-6
    roughness_fw_inboard_channel = 1.0e-6
    len_fw_channel = 4.0
    len_fw_inboard_channel = 4.0
    f_fw_peak = 1.0
    pres_blkt_coolant = 15.50e6
    temp_blkt_coolant_in = 573.0
    temp_blkt_coolant_out = 823.0
    coolp = 15.5e6
    n_blkt_outboard_modules_poloidal = 8
    n_blkt_inboard_modules_poloidal = 7
    n_blkt_outboard_modules_toroidal = 48
    n_blkt_inboard_modules_toroidal = 32
    temp_fw_max = 823.0
    fw_th_conductivity = 28.34
    fvoldw = 1.74
    fvolsi = 1.0
    fvolso = 0.64
    fwclfr = 0.15
    p_div_rad_total_mw = 0.0
    p_fw_rad_total_mw = 0.0
    p_fw_hcd_rad_total_mw = 0.0
    pradloss = 0.0
    p_tf_nuclear_heat_mw = 0.0
    ptfnucpm3 = 0.0
    r_cryostat_inboard = 0.0
    z_cryostat_half_inside = 0.0
    dr_pf_cryostat = 0.5
    vol_cryostat = 0.0
    vol_cryostat_internal = 0.0
    vol_vv = 0.0
    vfshld = 0.25
    vol_blkt_total = 0.0
    vol_blkt_inboard = 0.0
    vol_blkt_outboard = 0.0
    vol_shld_total = 0.0
    whtshld = 0.0
    wpenshld = 0.0
    wtshldi = 0.0
    wtshldo = 0.0
    irefprop = 1
    fblli = 0.0
    fblli2o = 0.08
    fbllipb = 0.68
    fblvd = 0.0
    m_blkt_li2o = 0.0
    wtbllipb = 0.0
    m_blkt_vanadium = 0.0
    m_blkt_lithium = 0.0
    f_a_blkt_cooling_channels = 0.25
    blktmodel = 0
    declblkt = 0.075
    declfw = 0.075
    declshld = 0.075
    blkttype = 3
    etaiso = 0.85
    eta_coolant_pump_electric = 0.95
    pnuc_cp = 0.0
    p_cp_shield_nuclear_heat_mw = 0.0
    pnuc_cp_tf = 0.0
    neut_flux_cp = 0.0
    i_fw_blkt_shared_coolant = 0
    i_blkt_liquid_breeder_type = 0
    i_blkt_dual_coolant = 0
    i_blkt_liquid_breeder_channel_type = 0
    i_blkt_module_segmentation = 0
    n_liq_recirc = 10
    r_f_liq_ib = 0.5
    r_f_liq_ob = 0.5
    w_f_liq_ib = 0.5
    w_f_liq_ob = 0.5
    den_ceramic = 3.21e3
    th_wall_secondary = 1.25e-2
    bz_channel_conduct_liq = 8.33e5
    a_bz_liq = 0.2
    b_bz_liq = 0.2
    nopol = 2
    nopipes = 4
    den_liq = 9.5e3
    specific_heat_liq = 1.9e2
    thermal_conductivity_liq = 30.0
    wht_liq = 0.0
    wht_liq_ib = 0.0
    wht_liq_ob = 0.0
    dynamic_viscosity_liq = 0.0
    electrical_conductivity_liq = 0.0
    hartmann_liq = [0.0, 0.0]
    b_mag_blkt = [5.0, 5.0]
    etaiso_liq = 0.85
    blpressure_liq = 1.7e6
    inlet_temp_liq = 570.0
    outlet_temp_liq = 720.0
    den_fw_coolant = 0.0
    visc_fw_coolant = 0.0
    den_blkt_coolant = 0.0
    visc_blkt_coolant = 0.0
    cp_fw = 0.0
    cv_fw = 0.0
    cp_bl = 0.0
    cv_bl = 0.0
    f_nuc_pow_bz_struct = 0.34
    f_nuc_pow_bz_liq = 0.66
    pnuc_fw_ratio_dcll = 0.14
    pnuc_blkt_ratio_dcll = 0.86
    n_blkt_inboard_module_coolant_sections_radial = 4
    n_blkt_inboard_module_coolant_sections_poloidal = 2
    n_blkt_outboard_module_coolant_sections_radial = 4
    n_blkt_outboard_module_coolant_sections_poloidal = 2
    bzfllengi_n_rad_liq = 2
    bzfllengi_n_pol_liq = 2
    bzfllengo_n_rad_liq = 2
    bzfllengo_n_pol_liq = 2
    radius_blkt_channel = 0.0
    radius_fw_channel_90_bend = 0.0
    radius_fw_channel_180_bend = 0.0
    radius_blkt_channel_90_bend = 0.0
    radius_blkt_channel_180_bend = 0.0
