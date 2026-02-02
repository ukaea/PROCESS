"""This library contains routines that can be shared by the blanket modules used in PROCESS.

These include:
- component_volumes
- component_masses
- thermo_hydraulic_model

author: G Graham, CCFE, Culham Science Centre

Acronyms for this module:

     BB          Breeding Blanket
     FW          First Wall
     BZ          Breeder Zone
     MF/BSS      Manifold/Back Supporting Structure
     LT          Low Temperature
     HT          High Temperature
     MMS         Multi Module Segment
     SMS         Single Modle Segment
     IB          Inboard
     OB          Outboard
     HCD         Heating & Current Drive
     FCI         Flow Channel Insert

Any changes within a subroutine or function code will have a comment explaining the change
"""

vol_shld_inboard: float = None
"""Volume of inboard shield (m3)"""


vol_shld_outboard: float = None
"""Volume of outboard shield (m3)"""


vol_vv_inboard: float = None
"""Volume of inboard Vacuum Vessel (m3)"""


vol_vv_outboard: float = None
"""Volume of outboard Vacuum Vessel (m3)"""


dz_pf_cryostat: float = None
"""Clearance between uppermost PF coil and cryostat lid (m)"""


vfblkti: float = None
"""Inboard void fraction of blanket"""


vfblkto: float = None
"""Outboard void fraction of blanket"""


len_blkt_inboard_coolant_channel_radial: float = None
"""Inboard blanket coolant channel length (radial direction) (m)"""


len_blkt_outboard_coolant_channel_radial: float = None
"""Outboard blanket coolant channel length (radial direction) (m)"""


len_blkt_inboard_segment_toroidal: float = None
"""Inboard blanket mid-plane toroidal circumference for segment (m)"""


len_blkt_outboard_segment_toroidal: float = None
"""Outboard blanket mid-plane toroidal circumference for segment (m)"""


len_blkt_inboard_segment_poloidal: float = None
"""Inboard blanket segment poloidal length (m)"""


len_blkt_outboard_segment_poloidal: float = None
"""Outboard blanket segment poloidal length (m)"""


len_blkt_inboard_channel_total: float = None
"""Inboard primary blanket flow lengths (m)"""


len_blkt_outboard_channel_total: float = None
"""Outboard primary blanket flow lengths (m)"""


bzfllengi_liq: float = None
"""Inboard secondary blanket flow lengths (m)"""


bzfllengo_liq: float = None
"""Outboard secondary blanket flow lengths (m)"""


p_fw_inboard_nuclear_heat_mw: float = None
"""Inboard first wall nuclear heating (MW)"""


p_fw_outboard_nuclear_heat_mw: float = None
"""Outboard first wall nuclear heating (MW)"""


temp_fw_inboard_peak: float = None
"""Inboard first wall peak temperature (K)"""


temp_fw_outboard_peak: float = None
"""Outboard first wall peak temperature (K)"""


mflow_fw_inboard_coolant_total: float = None
"""Inboard mass flow rate to remove inboard FW power (kg/s)"""


mflow_fw_outboard_coolant_total: float = None
"""Outboard mass flow rate to remove inboard FW power (kg/s)"""


mflow_fw_coolant_total: float = None
"""Total mass flow rate to remove inboard FW power (kg/s)"""


n_fw_inboard_channels: float = None
"""Inboard total number of first wall coolant channels"""


n_fw_outboard_channels: float = None
"""Outboard total number of first wall coolant channels"""


mflow_fw_inboard_coolant_channel: float = None
"""Inboard mass flow rate per coolant pipe (kg/s)"""


mflow_fw_outboard_coolant_channel: float = None
"""Outboard mass flow rate per coolant pipe (kg/s)"""


p_blkt_nuclear_heat_inboard_mw: float = None
"""Neutron power deposited inboard blanket blanket (MW)"""


p_blkt_nuclear_heat_outboard_mw: float = None
"""Neutron power deposited outboard blanket blanket (MW)"""


mflow_blkt_inboard_coolant: float = None
"""Inboard blanket mass flow rate for coolant (kg/s)"""


mflow_blkt_outboard_coolant: float = None
"""Outboard blanket mass flow rate for coolant (kg/s)"""


mflow_blkt_coolant_total: float = None
"""Total blanket mass flow rate for coolant (kg/s)"""


mfblkti_liq: float = None
"""Inboard blanket mass flow rate for liquid breeder (kg/s)"""


mfblkto_liq: float = None
"""Outboard blanket mass flow rate for liquid breeder (kg/s)"""


mfblkt_liq: float = None
"""Blanket mass flow rate for liquid breeder (kg/s)"""


mftotal: float = None
"""Total mass flow rate for coolant (kg/s)"""


n_blkt_inboard_channels: float = None
"""Inboard total number of blanket coolant pipes"""


n_blkt_outboard_channels: float = None
"""Outboard total number of blanket coolant pipes"""


mfblktpi: float = None
"""Inboard mass flow rate per coolant pipe (kg/s)"""


mfblktpo: float = None
"""Outboard mass flow rate per coolant pipe (kg/s)"""


vel_blkt_inboard_coolant: float = None
"""Inboard coolant velocity in blanket (m/s)"""


vel_blkt_outboard_coolant: float = None
"""Outboard coolant velocity in blanket (m/s)"""


htpmw_fwi: float = None
"""Inboard first wall pumping power (MW)"""


htpmw_fwo: float = None
"""Outboard first wall pumping power (MW)"""


htpmw_blkti: float = None
"""Inboard blanket pumping power (MW)"""


htpmw_blkto: float = None
"""Outboard blanket pumping power (MW)"""


htpmw_fw_blkti: float = None
"""Inboard fw and blanket pumping power (MW)"""


htpmw_fw_blkto: float = None
"""Outboard fw and blanket pumping power (MW)"""


dz_blkt_half: float = None
"""Blanket internal half-height (m)"""


dz_shld_half: float = None
"""Shield internal half-height (m)"""


dz_vv_half: float = None
"""Vacuum vessel internal half-height (m)"""


icomponent: int = None
"""Switch used to specify selected component: blanket=0, shield=1, vacuum vessel=2"""


def init_blanket_library():
    global dz_blkt_half
    global dz_shld_half
    global dz_pf_cryostat
    global dz_vv_half
    global vol_shld_inboard
    global vol_shld_outboard
    global vol_vv_inboard
    global vol_vv_outboard
    global len_blkt_inboard_coolant_channel_radial
    global len_blkt_outboard_coolant_channel_radial
    global len_blkt_inboard_segment_toroidal
    global len_blkt_outboard_segment_toroidal
    global len_blkt_inboard_segment_poloidal
    global len_blkt_outboard_segment_poloidal
    global len_blkt_inboard_channel_total
    global bzfllengi_liq
    global bzfllengo_liq
    global len_blkt_outboard_channel_total
    global p_fw_inboard_nuclear_heat_mw
    global p_fw_outboard_nuclear_heat_mw
    global temp_fw_inboard_peak
    global temp_fw_outboard_peak
    global mflow_fw_inboard_coolant_total
    global mflow_fw_outboard_coolant_total
    global mflow_fw_coolant_total
    global n_fw_inboard_channels
    global n_fw_outboard_channels
    global mflow_fw_inboard_coolant_channel
    global mflow_fw_outboard_coolant_channel
    global p_blkt_nuclear_heat_inboard_mw
    global p_blkt_nuclear_heat_outboard_mw
    global mflow_blkt_inboard_coolant
    global mflow_blkt_outboard_coolant
    global mfblkti_liq
    global mfblkto_liq
    global mfblkt_liq
    global mflow_blkt_coolant_total
    global mftotal
    global n_blkt_inboard_channels
    global n_blkt_outboard_channels
    global mfblktpi
    global mfblktpo
    global vel_blkt_inboard_coolant
    global vel_blkt_outboard_coolant
    global htpmw_fwi
    global htpmw_fwo
    global htpmw_blkti
    global htpmw_blkto
    global vfblkti
    global vfblkto

    dz_blkt_half = 0.0
    dz_shld_half = 0.0
    dz_pf_cryostat = 0.0
    dz_vv_half = 0.0
    vol_shld_inboard = 0.0
    vol_shld_outboard = 0.0
    vol_vv_inboard = 0.0
    vol_vv_outboard = 0.0
    len_blkt_inboard_coolant_channel_radial = 0.0
    len_blkt_outboard_coolant_channel_radial = 0.0
    len_blkt_inboard_segment_toroidal = 0.0
    len_blkt_outboard_segment_toroidal = 0.0
    len_blkt_inboard_segment_poloidal = 0.0
    len_blkt_outboard_segment_poloidal = 0.0
    len_blkt_inboard_channel_total = 0.0
    bzfllengi_liq = 0.0
    bzfllengo_liq = 0.0
    len_blkt_outboard_channel_total = 0.0
    p_fw_inboard_nuclear_heat_mw = 0.0
    p_fw_outboard_nuclear_heat_mw = 0.0
    temp_fw_inboard_peak = 0.0
    temp_fw_outboard_peak = 0.0
    mflow_fw_inboard_coolant_total = 0.0
    mflow_fw_outboard_coolant_total = 0.0
    mflow_fw_coolant_total = 0.0
    n_fw_inboard_channels = 0.0
    n_fw_outboard_channels = 0.0
    mflow_fw_inboard_coolant_channel = 0.0
    mflow_fw_outboard_coolant_channel = 0.0
    p_blkt_nuclear_heat_inboard_mw = 0.0
    p_blkt_nuclear_heat_outboard_mw = 0.0
    mflow_blkt_inboard_coolant = 0.0
    mflow_blkt_outboard_coolant = 0.0
    mfblkti_liq = 0.0
    mfblkto_liq = 0.0
    mfblkt_liq = 0.0
    mflow_blkt_coolant_total = 0.0
    mftotal = 0.0
    n_blkt_inboard_channels = 0.0
    n_blkt_outboard_channels = 0.0
    mfblktpi = 0.0
    mfblktpo = 0.0
    vel_blkt_inboard_coolant = 0.0
    vel_blkt_outboard_coolant = 0.0
    htpmw_fwi = 0.0
    htpmw_fwo = 0.0
    htpmw_blkti = 0.0
    htpmw_blkto = 0.0
    vfblkti = 0.0
    vfblkto = 0.0
