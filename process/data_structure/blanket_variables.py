"""Module containing variables for the blanket models

Acronyms for this module:

     BB          Breeding Blanket
     FW          First Wall
     BZ          Breeder Zone
     MF/BSS      Manifold/Back Supporting Structure
     LT          Low Temperature
     HT          High Temperature
     MMS         Multi Module Segment
     SMS         Single Module Segment
     IB          Inboard
     OB          Outboard
     HCD         Heating & Current Drive
     FCI         Flow Channel Insert
"""

from dataclasses import dataclass
from enum import IntEnum, unique

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@unique
class BlktModelTypes(IntEnum):
    """Enum for blanket model types. `i_blanket_type`"""

    CCFE_HCPB = 1
    DCLL = 5


@dataclass(slots=True)
class BlanketData(PROCESSModelData):
    """Dataclass holding blanket variables"""

    vol_shld_inboard: Parameter[float] = 0.0
    """Volume of inboard shield (m3)"""

    vol_shld_outboard: Parameter[float] = 0.0
    """Volume of outboard shield (m3)"""

    vol_vv_inboard: Parameter[float] = 0.0
    """Volume of inboard Vacuum Vessel (m3)"""

    vol_vv_outboard: Parameter[float] = 0.0
    """Volume of outboard Vacuum Vessel (m3)"""

    dz_pf_cryostat: Parameter[float] = 0.0
    """Clearance between uppermost PF coil and cryostat lid (m)"""

    vfblkti: Parameter[float] = 0.0
    """Inboard void fraction of blanket"""

    vfblkto: Parameter[float] = 0.0
    """Outboard void fraction of blanket"""

    len_blkt_inboard_coolant_channel_radial: Parameter[float] = 0.0
    """Inboard blanket coolant channel length (radial direction) (m)"""

    len_blkt_outboard_coolant_channel_radial: Parameter[float] = 0.0
    """Outboard blanket coolant channel length (radial direction) (m)"""

    len_blkt_inboard_segment_toroidal: Parameter[float] = 0.0
    """Inboard blanket mid-plane toroidal circumference for segment (m)"""

    len_blkt_outboard_segment_toroidal: Parameter[float] = 0.0
    """Outboard blanket mid-plane toroidal circumference for segment (m)"""

    len_blkt_inboard_segment_poloidal: Parameter[float] = 0.0
    """Inboard blanket segment poloidal length (m)"""

    len_blkt_outboard_segment_poloidal: Parameter[float] = 0.0
    """Outboard blanket segment poloidal length (m)"""

    len_blkt_inboard_channel_total: Parameter[float] = 0.0
    """Inboard primary blanket flow lengths (m)"""

    len_blkt_outboard_channel_total: Parameter[float] = 0.0
    """Outboard primary blanket flow lengths (m)"""

    bzfllengi_liq: Parameter[float] = 0.0
    """Inboard secondary blanket flow lengths (m)"""

    bzfllengo_liq: Parameter[float] = 0.0
    """Outboard secondary blanket flow lengths (m)"""

    p_fw_inboard_nuclear_heat_mw: Parameter[float] = 0.0
    """Inboard first wall nuclear heating (MW)"""

    p_fw_outboard_nuclear_heat_mw: Parameter[float] = 0.0
    """Outboard first wall nuclear heating (MW)"""

    temp_fw_inboard_peak: Parameter[float] = 0.0
    """Inboard first wall peak temperature (K)"""

    temp_fw_outboard_peak: Parameter[float] = 0.0
    """Outboard first wall peak temperature (K)"""

    mflow_fw_inboard_coolant_total: Parameter[float] = 0.0
    """Inboard mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_outboard_coolant_total: Parameter[float] = 0.0
    """Outboard mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_coolant_total: Parameter[float] = 0.0
    """Total mass flow rate to remove inboard FW power (kg/s)"""

    mflow_fw_inboard_coolant_channel: Parameter[float] = 0.0
    """Inboard mass flow rate per coolant pipe (kg/s)"""

    mflow_fw_outboard_coolant_channel: Parameter[float] = 0.0
    """Outboard mass flow rate per coolant pipe (kg/s)"""

    n_fw_inboard_channels: Parameter[float] = 0.0
    """Inboard total number of first wall coolant channels"""

    n_fw_outboard_channels: Parameter[float] = 0.0
    """Outboard total number of first wall coolant channels"""

    p_blkt_nuclear_heat_inboard_mw: Parameter[float] = 0.0
    """Neutron power deposited inboard blanket blanket (MW)"""

    p_blkt_nuclear_heat_outboard_mw: Parameter[float] = 0.0
    """Neutron power deposited outboard blanket blanket (MW)"""

    mflow_blkt_inboard_coolant: Parameter[float] = 0.0
    """Inboard blanket mass flow rate for coolant (kg/s)"""

    mflow_blkt_outboard_coolant: Parameter[float] = 0.0
    """Outboard blanket mass flow rate for coolant (kg/s)"""

    mflow_blkt_coolant_total: Parameter[float] = 0.0
    """Total blanket mass flow rate for coolant (kg/s)"""

    mfblkti_liq: Parameter[float] = 0.0
    """Inboard blanket mass flow rate for liquid breeder (kg/s)"""

    mfblkto_liq: Parameter[float] = 0.0
    """Outboard blanket mass flow rate for liquid breeder (kg/s)"""

    mfblkt_liq: Parameter[float] = 0.0
    """Blanket mass flow rate for liquid breeder (kg/s)"""

    mftotal: Parameter[float] = 0.0
    """Total mass flow rate for coolant (kg/s)"""

    n_blkt_inboard_channels: Parameter[float] = 0.0
    """Inboard total number of blanket coolant pipes"""

    n_blkt_outboard_channels: Parameter[float] = 0.0
    """Outboard total number of blanket coolant pipes"""

    mfblktpi: Parameter[float] = 0.0
    """Inboard mass flow rate per coolant pipe (kg/s)"""

    mfblktpo: Parameter[float] = 0.0
    """Outboard mass flow rate per coolant pipe (kg/s)"""

    vel_blkt_inboard_coolant: Parameter[float] = 0.0
    """Inboard coolant velocity in blanket (m/s)"""

    vel_blkt_outboard_coolant: Parameter[float] = 0.0
    """Outboard coolant velocity in blanket (m/s)"""

    htpmw_fwi: Parameter[float] = 0.0
    """Inboard first wall pumping power (MW)"""

    htpmw_fwo: Parameter[float] = 0.0
    """Outboard first wall pumping power (MW)"""

    htpmw_blkti: Parameter[float] = 0.0
    """Inboard blanket pumping power (MW)"""

    htpmw_blkto: Parameter[float] = 0.0
    """Outboard blanket pumping power (MW)"""

    htpmw_fw_blkti: Parameter[float] = None
    """Inboard fw and blanket pumping power (MW)"""

    htpmw_fw_blkto: Parameter[float] = None
    """Outboard fw and blanket pumping power (MW)"""

    dz_blkt_half: Parameter[float] = 0.0
    """Blanket internal half-height (m)"""

    dz_shld_half: Parameter[float] = 0.0
    """Shield internal half-height (m)"""

    dz_vv_half: Parameter[float] = 0.0
    """Vacuum vessel internal half-height (m)"""

    deg_blkt_outboard_poloidal_plasma: Parameter[float] = 0.0
    """Outboard blanket poloidal angle subtended by plasma (degrees)"""

    f_deg_blkt_outboard_poloidal_plasma: Parameter[float] = 0.0
    """Fraction of outboard blanket poloidal angle subtended by plasma (degrees)"""

    deg_blkt_inboard_poloidal_plasma: Parameter[float] = 0.0
    """Inboard blanket poloidal angle subtended by plasma (degrees)"""

    f_deg_blkt_inboard_poloidal_plasma: Parameter[float] = 0.0
    """Fraction of inboard blanket poloidal angle subtended by plasma (degrees)"""


CREATE_DICTS_FROM_DATACLASS = BlanketData
