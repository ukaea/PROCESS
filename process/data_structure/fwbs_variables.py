"""
Module containing global variables relating to the first wall, blanket and
shield components
### References
-
"""

from dataclasses import dataclass, field


@dataclass
class FWBSData:
    life_blkt_fpy: float = 0.0
    """Full power blanket lifetime (years)"""

    life_blkt: float = 0.0
    """Calendar year blanket lifetime (years)"""

    m_fw_blkt_div_coolant_total: float = 0.0
    """mass of water coolant (in shield, blanket, first wall, divertor) [kg]"""

    m_vv: float = 0.0
    """vacuum vessel mass [kg]"""

    den_steel: float = 7800.0
    """density of steel [kg m^-3]"""

    denwc: float = 15630.0
    """density of tungsten carbide [kg m^-3]"""

    dewmkg: float = 0.0
    """total mass of vacuum vessel + cryostat [kg] (calculated if blktmodel>0)"""

    f_p_blkt_multiplication: float = 1.269
    """energy multiplication in blanket and shield"""

    p_blkt_multiplication_mw: float = 0.0
    """power due to energy multiplication in blanket and shield [MW]"""

    fblss: float = 0.09705
    """KIT blanket model: steel fraction of breeding zone"""

    f_ster_div_single: float = 0.115
    """Solid angle fraction taken by one divertor"""

    f_a_fw_outboard_hcd: float = 0.0
    """area fraction of first wall covered by heating/current drive apparatus plus diagnostics"""

    fhole: float = 0.0
    """area fraction taken up by other holes (IFE)"""

    i_fw_blkt_vv_shape: int = 2
    """switch for first wall, blanket, shield and vacuum vessel shape:
    - =1 D-shaped (cylinder inboard + ellipse outboard)
    - =2 defined by two ellipses
    """

    life_fw_fpy: float = 0.0
    """first wall full-power year lifetime (y)"""

    m_fw_total: float = 0.0
    """first wall mass [kg]"""

    fw_armour_mass: float = 0.0
    """first wall armour mass [kg]"""

    fw_armour_thickness: float = 0.005
    """first wall armour thickness [m]"""

    fw_armour_vol: float = 0.0
    """first wall armour volume [m^3]"""

    i_blanket_type: int = 1
    """switch for blanket model:
    - =1 CCFE HCPB model
    - =2 KIT HCPB model  # REMOVED, no longer usable
    - =3 CCFE HCPB model with Tritium Breeding Ratio calculation # REMOVED, no longer usable
    - =4 KIT HCLL model  # REMOVED, no longer usable
    - =5 DCLL model -  no nutronics model included (in development) please check/choose values for
    'dual-coolant blanket' fractions (provided in this file).
    -  please use i_p_coolant_pumping = 0 or 1.
    """

    i_blkt_inboard: int = 1
    """switch for inboard blanket:
    - =0 No inboard blanket (dr_blkt_inboard=0.0)
    - =1 Inboard blanket present
    """

    inuclear: int = 0
    """switch for nuclear heating in the coils:
    - =0 Frances Fox model (default)
    - =1 Fixed by user (qnuc)
    """

    qnuc: float = 0.0
    """nuclear heating in the coils (W) (`inuclear=1`)"""

    f_blkt_li6_enrichment: float = 30.0
    """lithium-6 enrichment of breeding material (%)"""

    p_blkt_nuclear_heat_total_mw: float = 0.0
    """nuclear heating in the blanket [MW]"""

    pnuc_cp: float = 0.0
    """Total nuclear heating in the ST centrepost [MW]"""

    p_cp_shield_nuclear_heat_mw: float = 0.0
    """Neutronic shield nuclear heating in the ST centrepost [MW]"""

    pnuc_cp_tf: float = 0.0
    """TF neutronic nuclear heating in the ST centrepost [MW]"""

    p_div_nuclear_heat_total_mw: float = 0.0
    """nuclear heating in the divertor [MW]"""

    p_fw_nuclear_heat_total_mw: float = 0.0
    """nuclear heating in the first wall [MW]"""

    p_fw_hcd_nuclear_heat_mw: float = 0.0
    """Nuclear heating in the HCD apparatus and diagnostics on the first wall [MW]"""

    pnucloss: float = 0.0
    """nuclear heating lost via holes [MW]"""

    pnucvvplus: float = 0.0
    """nuclear heating to vacuum vessel and beyond [MW]"""

    p_shld_nuclear_heat_mw: float = 0.0
    """nuclear heating in the shield [MW]"""

    m_blkt_total: float = 0.0
    """mass of blanket [kg]"""

    m_blkt_steel_total: float = 0.0
    """mass of blanket - steel part [kg]"""

    armour_fw_bl_mass: float = 0.0
    """Total mass of armour, first wall and blanket [kg]"""

    # CCFE HCPB Blanket Model i_blanket_type=1

    breeder_f: float = 0.5
    """Volume ratio: Li4SiO4/(Be12Ti+Li4SiO4) (`iteration variable 108`)"""

    breeder_multiplier: float = 0.75
    """combined breeder/multipler fraction of blanket by volume"""

    vfcblkt: float = 0.05295
    """He coolant fraction of blanket by volume (`i_blanket_type= 1` (CCFE HCPB))"""

    vfpblkt: float = 0.1
    """He purge gas fraction of blanket by volume (`i_blanket_type= 1` (CCFE HCPB))"""

    m_blkt_li4sio4: float = 0.0
    """mass of lithium orthosilicate in blanket [kg] (`i_blanket_type=1` (CCFE HCPB))"""

    m_blkt_tibe12: float = 0.0
    """mass of titanium beryllide in blanket [kg] (`i_blanket_type=1` (CCFE HCPB))"""

    neut_flux_cp: float = 0.0
    """Centrepost TF fast neutron flux (E > 0.1 MeV) [m^(-2).^(-1)]
    This variable is only calculated for superconducting (i_tf_sup = 1 )
    spherical tokamal magnet designs (itart = 0)
    """

    f_neut_shield: float = -1.0
    """Fraction of nuclear power shielded before the CP magnet (ST)
    ( neut_absorb = -1 --> a fit on simplified MCNP neutronic
    calculation is used assuming water cooled (13%) tungesten carbyde )
    """

    f_a_fw_coolant_inboard: float = 0.0
    """Inboard FW coolant cross-sectional area void fraction"""

    f_a_fw_coolant_outboard: float = 0.0
    """Outboard FW coolant cross-sectional area void fraction"""

    psurffwi: float = 0.0
    """Surface heat flux on first wall [MW] (sum = p_fw_rad_total_mw)"""

    psurffwo: float = 0.0
    """Surface heat flux on first wall [MW] (sum = p_fw_rad_total_mw)"""

    vol_fw_total: float = 0.0
    """First wall volume [m3]"""

    f_vol_blkt_steel: float = 0.0
    """Fractions of blanket by volume: steel"""

    f_vol_blkt_li4sio4: float = 0.0
    """Fractions of blanket by volume: lithium orthosilicate"""

    f_vol_blkt_tibe12: float = 0.0
    """Fractions of blanket by volume: titanium beryllide"""

    breedmat: int = 1
    """breeder material switch:
    - =1 Lithium orthosilicate
    - =2 Lithium methatitanate
    - =3 Lithium zirconate
    """

    densbreed: float = 0.0
    """density of breeder material [kg m^-3]"""

    fblbe: float = 0.6
    """beryllium fraction of blanket by volume"""

    fblbreed: float = 0.154
    """breeder fraction of blanket breeding zone by volume"""

    fblhebmi: float = 0.4
    """helium fraction of inboard blanket box manifold by volume"""

    fblhebmo: float = 0.4
    """helium fraction of outboard blanket box manifold by volume """

    fblhebpi: float = 0.6595
    """helium fraction of inboard blanket back plate by volume """

    fblhebpo: float = 0.6713
    """helium fraction of outboard blanket back plate by volume """

    hcdportsize: int = 1
    """switch for size of heating/current drive ports :
    - =1 'small'
    - =2 'large'
    """

    nflutf: float = 0.0
    """peak fast neutron fluence on TF coil superconductor [n m^-2] """

    npdiv: int = 2
    """number of divertor ports """

    nphcdin: int = 2
    """number of inboard ports for heating/current drive """

    nphcdout: int = 2
    """number of outboard ports for heating/current drive """

    tbr: float = 0.0
    """tritium breeding ratio"""

    tritprate: float = 0.0
    """tritium production rate [g day^-1] """

    wallpf: float = 1.21
    """neutron wall load peaking factor """

    whtblbreed: float = 0.0
    """mass of blanket - breeder part [kg] """

    m_blkt_beryllium: float = 0.0
    """mass of blanket - beryllium part [kg]"""

    i_p_coolant_pumping: int = 2
    """Switch for pumping power for primary coolant (mechanical power only and peak first wall
    temperature is only calculated if `i_p_coolant_pumping=2`):
    - =0 User sets pump power directly (p_blkt_coolant_pump_mw, p_fw_coolant_pump_mw, p_div_coolant_pump_mw, p_shld_coolant_pump_mw)
    - =1 User sets pump power as a fraction of thermal power (f_p_blkt_coolant_pump_total_heat, f_p_fw_coolant_pump_total_heat, f_p_div_coolant_pump_total_heat, f_p_shld_coolant_pump_total_heat)
    - =2 Mechanical pumping power is calculated
    - =3 Mechanical pumping power is calculated using specified pressure drop
    """

    i_shield_mat: int = 0
    """Switch for shield material - *currently only applied in costing routines* `cost_model = 2`
    - =0 Tungsten (default)
    - =1 Tungsten carbide
    """

    i_thermal_electric_conversion: int = 0
    """Switch for power conversion cycle:
    - =0 Set efficiency for chosen blanket, from detailed models (divertor heat not used)
    - =1 Set efficiency for chosen blanket, from detailed models (divertor heat used)
    - =2 user input thermal-electric efficiency (eta_turbine)
    - =3 steam Rankine cycle
    - =4 supercritical CO2 cycle
    """

    secondary_cycle_liq: int = 4
    """Switch for power conversion cycle for the liquid breeder component of the blanket:
    - =2 user input thermal-electric efficiency (eta_turbine)
    - =4 supercritical CO2 cycle
    """

    i_blkt_coolant_type: int = 1
    """Switch for blanket coolant (set via blkttype):
    - =1 helium
    - =2 pressurized water
    """

    i_fw_coolant_type: str = "helium"
    """switch for first wall coolant (can be different from blanket coolant):
    - 'helium'
    - 'water'
    """

    dr_fw_wall: float = 0.003
    """wall thickness of first wall coolant channels [m]"""

    radius_fw_channel: float = 0.006
    """radius of first wall cooling channels [m]"""

    dx_fw_module: float = 0.02
    """Width of a FW module containing a cooling channel [m]"""

    temp_fw_coolant_in: float = 573.0
    """inlet temperature of first wall coolant [K]"""

    temp_fw_coolant_out: float = 823.0
    """outlet temperature of first wall coolant [K]"""

    pres_fw_coolant: float = 15.5e6
    """first wall coolant pressure [Pa] (`i_thermal_electric_conversion>1`)"""

    temp_fw_peak: float = 873.0
    """peak first wall temperature [K]"""

    roughness_fw_channel: float = 1.0e-6
    """first wall channel roughness epsilon [m]"""

    len_fw_channel: float = 4.0
    """Length of a single first wall channel (all in parallel) [m]
    (`iteration variable 114`, useful for `constraint equation 39`)
    """

    f_fw_peak: float = 1.0
    """peaking factor for first wall heat loads. (Applied separately to inboard and outboard loads.
    Applies to both neutron and surface loads. Only used to calculate peak temperature - not
    the coolant flow rate.)
    """

    pres_blkt_coolant: float = 15.50e6
    """blanket coolant pressure [Pa] (`i_thermal_electric_conversion>1`)"""

    temp_blkt_coolant_in: float = 573.0
    """inlet temperature of blanket coolant  [K] (`i_thermal_electric_conversion>1`)"""

    temp_blkt_coolant_out: float = 823.0
    """Outlet temperature of blanket coolant [K] (`i_thermal_electric_conversion>1`)
    - input if `i_blkt_coolant_type=1` (helium)
    - calculated if `i_blkt_coolant_type=2` (water)
    """

    coolp: float = 15.5e6
    """blanket coolant pressure [Pa] (stellarator only)"""

    n_blkt_outboard_modules_poloidal: int = 8
    """number of outboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)"""

    n_blkt_inboard_modules_poloidal: int = 7
    """number of inboard blanket modules in poloidal direction (`i_thermal_electric_conversion>1`)"""

    n_blkt_outboard_modules_toroidal: int = 48
    """number of outboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)"""

    n_blkt_inboard_modules_toroidal: int = 32
    """number of inboard blanket modules in toroidal direction (`i_thermal_electric_conversion>1`)"""

    temp_fw_max: float = 823.0
    """maximum temperature of first wall material [K] (`i_thermal_electric_conversion>1`)"""

    fw_th_conductivity: float = 28.34
    """thermal conductivity of first wall material at 293 K (W/m/K) (Temperature dependence
    is as for unirradiated Eurofer)
    """

    fvoldw: float = 1.74
    """area coverage factor for vacuum vessel volume"""

    fvolsi: float = 1.0
    """area coverage factor for inboard shield volume"""

    fvolso: float = 0.64
    """area coverage factor for outboard shield volume"""

    fwclfr: float = 0.15
    """first wall coolant fraction (calculated if `i_pulsed_plant=1` or `ipowerflow=1`)"""

    p_div_rad_total_mw: float = 0.0
    """Total radiation power incident on the divertor(s) (MW)"""

    p_fw_rad_total_mw: float = 0.0
    """Radiation power incident on the first wall (MW)"""

    p_fw_hcd_rad_total_mw: float = 0.0
    """Radiation power incident on the heating and current drive systems on the first wall (MW)"""

    pradloss: float = 0.0
    """Radiation power lost through holes (eventually hits shield) (MW)
    Only used for stellarator
    """

    p_tf_nuclear_heat_mw: float = 0.0
    """nuclear heating in the TF coil (MW)"""

    ptfnucpm3: float = 0.0
    """nuclear heating in the TF coil (MW/m3) (`blktmodel>0`)"""

    r_cryostat_inboard: float = 0.0
    """cryostat radius [m]"""

    z_cryostat_half_inside: float = 0.0
    """cryostat height [m]"""

    dr_pf_cryostat: float = 0.5
    """Radial distance between outer edge of furthest away PF coil (or stellarator
    modular coil) and cryostat [m]
    """

    vol_cryostat: float = 0.0
    """Cryostat structure volume [m^3]"""

    vol_cryostat_internal: float = 0.0
    """Internal volume of the cryostat [m^3]"""

    vol_vv: float = 0.0
    """vacuum vessel volume [m^3]"""

    vfshld: float = 0.25
    """coolant void fraction in shield"""

    vol_blkt_total: float = 0.0
    """volume of blanket [m^3]"""

    vol_blkt_total_full_coverage: float = 0.0
    """Volume of blanket with no holes or ports (toroidally continuous) [m³]"""

    vol_blkt_inboard: float = 0.0
    """volume of inboard blanket [m^3]"""

    vol_blkt_inboard_full_coverage: float = 0.0
    """Volume of inboard blanket with no holes or ports (toroidally continuous) [m³]"""

    vol_blkt_outboard: float = 0.0
    """volume of outboard blanket [m^3]"""

    vol_blkt_outboard_full_coverage: float = 0.0
    """Volume of outboard blanket with no holes or ports (toroidally continuous) [m³]"""

    vol_shld_total: float = 0.0
    """volume of shield [m^3]"""

    whtshld: float = 0.0
    """mass of shield [kg]"""

    wpenshld: float = 0.0
    """mass of the penetration shield [kg]"""

    wtshldi: float = 0.0
    """mass of inboard shield [kg]"""

    wtshldo: float = 0.0
    """mass of outboard shield [kg]"""

    irefprop: int = 1
    """Switch to use REFPROP routines (stellarator only)"""

    fblli: float = 0.0
    """lithium fraction of blanket by volume (stellarator only)"""

    fblli2o: float = 0.08
    """lithium oxide fraction of blanket by volume (stellarator only)"""

    fbllipb: float = 0.68
    """lithium lead fraction of blanket by volume (stellarator only)"""

    fblvd: float = 0.0
    """vanadium fraction of blanket by volume (stellarator only)"""

    m_blkt_li2o: float = 0.0
    """mass of blanket - Li_2O part [kg]"""

    wtbllipb: float = 0.0
    """mass of blanket - Li-Pb part [kg]"""

    m_blkt_vanadium: float = 0.0
    """mass of blanket - vanadium part [kg]"""

    m_blkt_lithium: float = 0.0
    """mass of blanket - lithium part [kg]"""

    f_a_blkt_cooling_channels: float = 0.25
    """coolant void fraction in blanket."""

    blktmodel: int = 0
    """switch for blanket/tritium breeding model (see i_blanket_type):
    - =0 original simple model
    - =1 KIT model based on a helium-cooled pebble-bed blanket (HCPB) reference design
    """

    declblkt: float = 0.075
    """neutron power deposition decay length of blanket structural material [m] (stellarators only)"""

    declfw: float = 0.075
    """neutron power deposition decay length of first wall structural material [m] (stellarators only)"""

    declshld: float = 0.075
    """neutron power deposition decay length of shield structural material [m] (stellarators only)"""

    blkttype: int = 3
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

    etaiso: float = 0.85
    """isentropic efficiency of FW and blanket coolant pumps"""

    eta_coolant_pump_electric: float = 0.95
    """electrical efficiency of primary coolant pumps"""

    i_fw_blkt_shared_coolant: int = 0
    """Switch for whether the FW and BB are on the same pump system
    i.e. do they have the same primary coolant or not
    - =0    FW and BB have the same primary coolant, flow = FWin->FWout->BBin->BBout
    - =1    FW and BB have the different primary coolant and are on different pump systems
    """

    i_blkt_liquid_breeder_type: int = 0
    """Switch for Liquid Metal Breeder Material
    - =0   PbLi
    - =1   Li
    """

    i_blkt_dual_coolant: int = 0
    """Switch to specify whether breeding blanket is single-cooled or dual-coolant.
    - =0    Single coolant used for FW and Blanket (H2O or He). Solid Breeder.
    - =1    Single coolant used for FW and Blanket (H2O or He). Liquid metal breeder
    circulted for tritium extraction.
    - =2    Dual coolant: primary coolant (H2O or He) for FW and blanket structure;
    secondary coolant is self-cooled liquid metal breeder.
    """

    i_blkt_liquid_breeder_channel_type: int = 0
    """Switch for Flow Channel Insert (FCI) type if liquid metal breeder blanket.
    - =0    Thin conducting walls, default electrical conductivity (bz_channel_conduct_liq) is Eurofer
    - =1    Insulating Material, assumed perfect electrical insulator, default density (den_ceramic) is for SiC
    - =2    Insulating Material, electrical conductivity (bz_channel_conduct_liq) is input (default Eurofer), default density (den_ceramic) is for SiC
    """

    i_blkt_module_segmentation: int = 0
    """Switch for Multi Module Segment (MMS) or Single Modle Segment (SMS)
    - =0    MMS
    - =1    SMS
    """

    n_liq_recirc: int = 10
    """Number of liquid metal breeder recirculations per day, for use with i_blkt_dual_coolant=1"""

    r_f_liq_ib: float = 0.5
    """Radial fraction of BZ liquid channels"""

    r_f_liq_ob: float = 0.5
    """Radial fraction of BZ liquid channels"""

    w_f_liq_ib: float = 0.5
    """Toroidal fraction of BZ liquid channels"""

    w_f_liq_ob: float = 0.5
    """Toroidal fraction of BZ liquid channels"""

    den_ceramic: float = 3.21e3
    """FCI material density"""

    th_wall_secondary: float = 1.25e-2
    """Liquid metal coolant/breeder wall thickness thin conductor or FCI [m]"""

    bz_channel_conduct_liq: float = 8.33e5
    """Liquid metal coolant/breeder thin conductor or FCI wall conductance [A V^-1 m^-1]"""

    a_bz_liq: float = 0.2
    """Toroidal width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone"""

    b_bz_liq: float = 0.2
    """Radial width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone"""

    nopol: int = 2
    """Number of poloidal sections in a liquid metal breeder/coolant channel for module/segment"""

    nopipes: int = 4
    """Number of Liquid metal breeder/coolant channels per module/segment"""

    den_liq: float = 9.5e3
    """Liquid metal breeder/coolant density [kg m^-3]"""

    wht_liq: float = 0.0
    """Liquid metal"""

    wht_liq_ib: float = 0.0
    """Liquid metal"""

    wht_liq_ob: float = 0.0
    """Liquid metal"""

    specific_heat_liq: float = 1.9e2
    """Liquid metal breeder/coolant specific heat [J kg^-1 K^-1]"""

    thermal_conductivity_liq: float = 30.0
    """Liquid metal breeder/coolant thermal conductivity [W m^-1 K^-1]"""

    dynamic_viscosity_liq: float = 0.0
    """Liquid metal breeder/coolant dynamic viscosity [Pa s]"""

    electrical_conductivity_liq: float = 0.0
    """Liquid metal breeder/coolant electrical conductivity [Ohm m]"""

    hartmann_liq: list[float] = field(default_factory=lambda: [0.0, 0.0])
    """Hartmann number"""

    b_mag_blkt: list[float] = field(default_factory=lambda: [5.0, 5.0])
    """Toroidal Magnetic feild strength for IB/OB blanket [T]"""

    etaiso_liq: float = 0.85
    """Isentropic efficiency of blanket liquid breeder/coolant pumps"""

    blpressure_liq: float = 1.7e6
    """blanket liquid metal breeder/coolant pressure [Pa]"""

    inlet_temp_liq: float = 570.0
    """Inlet (scan var 68) temperature of the liquid breeder/coolant [K]"""

    outlet_temp_liq: float = 720.0
    """Outlet (scan var 69) temperature of the liquid breeder/coolant [K]"""

    den_fw_coolant: float = 0.0
    """Density of the FW primary coolant"""

    visc_fw_coolant: float = 0.0
    """Viscosity of the FW primary coolant"""

    den_blkt_coolant: float = 0.0
    """Density of the blanket primary coolant"""

    visc_blkt_coolant: float = 0.0
    """Viscosity of the blanket primary coolant"""

    cp_fw: float = 0.0
    """Spesific heat for FW and blanket primary coolant(s)"""

    cv_fw: float = 0.0
    """Spesific heat for FW and blanket primary coolant(s)"""

    cp_bl: float = 0.0
    """Spesific heat for FW and blanket primary coolant(s)"""

    cv_bl: float = 0.0
    """Spesific heat for FW and blanket primary coolant(s)"""

    f_nuc_pow_bz_struct: float = 0.34
    """For a dual-coolant blanket, fraction of BZ power cooled by primary coolant"""

    f_nuc_pow_bz_liq: float = 0.66
    """For a dual-coolant blanket, fraction of BZ self-cooled power (secondary coolant)"""

    pnuc_fw_ratio_dcll: float = 0.14
    """For a dual-coolant blanket, ratio of FW nuclear power as fraction of total"""

    pnuc_blkt_ratio_dcll: float = 0.86
    """For a dual-coolant blanket, ratio of Blanket nuclear power as fraction of total"""

    n_blkt_inboard_module_coolant_sections_radial: int = 4
    """Number of radial and poloidal sections that make up the total primary coolant flow
    length in a blanket module (IB and OB)
    """

    n_blkt_inboard_module_coolant_sections_poloidal: int = 2
    """Number of radial and poloidal sections that make up the total primary coolant flow
    length in a blanket module (IB and OB)
    """

    n_blkt_outboard_module_coolant_sections_radial: int = 4
    """Number of radial and poloidal sections that make up the total primary coolant flow
    length in a blanket module (IB and OB)
    """

    n_blkt_outboard_module_coolant_sections_poloidal: int = 2
    """Number of radial and poloidal sections that make up the total primary coolant flow
    length in a blanket module (IB and OB)
    """

    bzfllengi_n_rad_liq: int = 2
    """Number of radial and poloidal sections that make up the total secondary coolant/breeder
    flow length in a blanket module (IB and OB)
    """

    bzfllengi_n_pol_liq: int = 2
    """Number of radial and poloidal sections that make up the total secondary coolant/breeder
    flow length in a blanket module (IB and OB)
    """

    bzfllengo_n_rad_liq: int = 2
    """Number of radial and poloidal sections that make up the total secondary coolant/breeder
    flow length in a blanket module (IB and OB)
    """

    bzfllengo_n_pol_liq: int = 2
    """Number of radial and poloidal sections that make up the total secondary coolant/breeder
    flow length in a blanket module (IB and OB)
    """

    radius_blkt_channel: float = 0.0
    """Radius of blanket cooling channels [m]"""

    radius_blkt_channel_90_bend: float = 0.0
    """Radius of blanket cooling channel 90° bend [m]"""

    radius_fw_channel_90_bend: float = 0.0
    """Radius of first wall cooling channel 90° bend [m]"""

    radius_fw_channel_180_bend: float = 0.0
    """Radius of first wall cooling channel 180° bend [m]"""

    radius_blkt_channel_180_bend: float = 0.0
    """Radius of blanket cooling channel 180° bend [m]"""

    dz_fw_half: float = 0.0
    """Half-height of first wall structure [m]"""


CREATE_DICTS_FROM_DATACLASS = FWBSData
