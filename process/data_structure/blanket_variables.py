"""This library contains routines that can be shared by the blanket modules used in PROCESS.

These include:
- component_volumes
- component_masses
- thermo_hydraulic_model



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

from dataclasses import dataclass


@dataclass
class BlanketData:
    vol_shld_inboard: float = 0.0
    """Volume of inboard shield (m3)"""

    vol_shld_outboard: float = 0.0
    """Volume of outboard shield (m3)"""

    vol_vv_inboard: float = 0.0
    """Volume of inboard Vacuum Vessel (m3)"""

    vol_vv_outboard: float = 0.0
    """Volume of outboard Vacuum Vessel (m3)"""

    dz_pf_cryostat: float = 0.0
    """Clearance between uppermost PF coil and cryostat lid (m)"""

    vfblkti: float = 0.0
    """Inboard void fraction of blanket"""

    vfblkto: float = 0.0
    """Outboard void fraction of blanket"""

    len_blkt_inboard_coolant_channel_radial: float = 0.0
    """Inboard blanket coolant channel length (radial direction) (m)"""

    len_blkt_outboard_coolant_channel_radial: float = 0.0
    """Outboard blanket coolant channel length (radial direction) (m)"""

    len_blkt_inboard_segment_toroidal: float = 0.0
    """Inboard blanket mid-plane toroidal circumference for segment (m)"""

    len_blkt_outboard_segment_toroidal: float = 0.0
    """Outboard blanket mid-plane toroidal circumference for segment (m)"""

    len_blkt_inboard_segment_poloidal: float = 0.0
    """Inboard blanket segment poloidal length (m)"""

    len_blkt_outboard_segment_poloidal: float = 0.0
    """Outboard blanket segment poloidal length (m)"""

    len_blkt_inboard_channel_total: float = 0.0
    """Inboard primary blanket flow lengths (m)"""

    len_blkt_outboard_channel_total: float = 0.0
    """Outboard primary blanket flow lengths (m)"""

    bzfllengi_liq: float = 0.0
    """Inboard secondary blanket flow lengths (m)"""

    bzfllengo_liq: float = 0.0
    """Outboard secondary blanket flow lengths (m)"""

    p_fw_inboard_nuclear_heat_mw: float = 0.0
    """Inboard first wall nuclear heating (MW)"""

    p_fw_outboard_nuclear_heat_mw: float = 0.0
    """Outboard first wall nuclear heating (MW)"""

    temp_fw_inboard_peak: float = 0.0
    """Inboard first wall peak temperature (K)"""

    temp_fw_outboard_peak: float = 0.0
    """Outboard first wall peak temperature (K)"""

    mflow_fw_inboard_coolant_total: float = 0.0
    """Inboard mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_outboard_coolant_total: float = 0.0
    """Outboard mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_coolant_total: float = 0.0
    """Total mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_inboard_coolant_channel: float = 0.0
    """Inboard mass flow rate per coolant pipe (kg/s)"""

    mflow_fw_outboard_coolant_channel: float = 0.0
    """Outboard mass flow rate per coolant pipe (kg/s)"""

    n_fw_inboard_channels: float = 0.0
    """Inboard total number of first wall coolant channels"""

    n_fw_outboard_channels: float = 0.0
    """Outboard total number of first wall coolant channels"""

    p_blkt_nuclear_heat_inboard_mw: float = 0.0
    """Neutron power deposited inboard blanket blanket (MW)"""

    p_blkt_nuclear_heat_outboard_mw: float = 0.0
    """Neutron power deposited outboard blanket blanket (MW)"""

    mflow_blkt_inboard_coolant: float = 0.0
    """Inboard blanket mass flow rate for coolant (kg/s)"""

    mflow_blkt_outboard_coolant: float = 0.0
    """Outboard blanket mass flow rate for coolant (kg/s)"""

    mflow_blkt_coolant_total: float = 0.0
    """Total blanket mass flow rate for coolant (kg/s)"""

    mfblkti_liq: float = 0.0
    """Inboard blanket mass flow rate for liquid breeder (kg/s)"""

    mfblkto_liq: float = 0.0
    """Outboard blanket mass flow rate for liquid breeder (kg/s)"""

    mfblkt_liq: float = 0.0
    """Blanket mass flow rate for liquid breeder (kg/s)"""

    mftotal: float = 0.0
    """Total mass flow rate for coolant (kg/s)"""

    n_blkt_inboard_channels: float = 0.0
    """Inboard total number of blanket coolant pipes"""

    n_blkt_outboard_channels: float = 0.0
    """Outboard total number of blanket coolant pipes"""

    mfblktpi: float = 0.0
    """Inboard mass flow rate per coolant pipe (kg/s)"""

    mfblktpo: float = 0.0
    """Outboard mass flow rate per coolant pipe (kg/s)"""

    vel_blkt_inboard_coolant: float = 0.0
    """Inboard coolant velocity in blanket (m/s)"""

    vel_blkt_outboard_coolant: float = 0.0
    """Outboard coolant velocity in blanket (m/s)"""

    htpmw_fwi: float = 0.0
    """Inboard first wall pumping power (MW)"""

    htpmw_fwo: float = 0.0
    """Outboard first wall pumping power (MW)"""

    htpmw_blkti: float = 0.0
    """Inboard blanket pumping power (MW)"""

    htpmw_blkto: float = 0.0
    """Outboard blanket pumping power (MW)"""

    htpmw_fw_blkti: float = None
    """Inboard fw and blanket pumping power (MW)"""

    htpmw_fw_blkto: float = None
    """Outboard fw and blanket pumping power (MW)"""

    dz_blkt_half: float = 0.0
    """Blanket internal half-height (m)"""

    dz_shld_half: float = 0.0
    """Shield internal half-height (m)"""

    dz_vv_half: float = 0.0
    """Vacuum vessel internal half-height (m)"""

    deg_blkt_outboard_poloidal_plasma: float = 0.0
    """Outboard blanket poloidal angle subtended by plasma (degrees)"""

    f_deg_blkt_outboard_poloidal_plasma: float = 0.0
    """Fraction of outboard blanket poloidal angle subtended by plasma (degrees)"""

    deg_blkt_inboard_poloidal_plasma: float = 0.0
    """Inboard blanket poloidal angle subtended by plasma (degrees)"""

    f_deg_blkt_inboard_poloidal_plasma: float = 0.0
    """Fraction of inboard blanket poloidal angle subtended by plasma (degrees)"""


CREATE_DICTS_FROM_DATACLASS = BlanketData
