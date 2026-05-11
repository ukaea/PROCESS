"""
This module contains the PROCESS CCFE HCPB blanket model
based on CCFE HCPB model from the PROCESS engineering paper
PROCESS Engineering paper (M. Kovari et al.)
### References
- Kovari et al., Fusion Engineering and Design 104 (2016) 9-20
"""

from dataclasses import dataclass


@dataclass
class CCFEHCPBData:
    armour_density: float = 0.0
    """FW armour density [kg/m3]"""

    fw_density: float = 0.0
    """FW density [kg/m3]"""

    blanket_density: float = 0.0
    """Blanket density [kg/m3]"""

    shield_density: float = 0.0
    """Shield density [kg/m3]"""

    vv_density: float = 0.0
    """Vacuum vessel density [kg/m3]"""

    x_blanket: float = 0.0
    """Blanket exponent (tonne/m2)"""

    x_shield: float = 0.0
    """Shield exponent (tonne/m2)"""

    tfc_nuc_heating: float = 0.0
    """Unit nuclear heating in TF coil (W per W of fusion power)"""

    fw_armour_u_nuc_heating: float = 6.25e-7
    """Unit heating of FW and armour in FW armour (W/kg per W of fusion power)"""

    shld_u_nuc_heating: float = 0.0
    """Unit nuclear heating in shield (W per W of fusion power)"""

    pnuc_tot_blk_sector: float = None
    """Total nuclear power deposited in blanket covered sector (FW, BLKT, SHLD, TF) (MW)"""

    exp_blanket: float = 0.0
    """Exponential factors in nuclear heating calcs"""

    exp_shield1: float = 0.0
    """Exponential factors in nuclear heating calcs"""

    exp_shield2: float = 0.0
    """Exponential factors in nuclear heating calcs"""


CREATE_DICTS_FROM_DATACLASS = CCFEHCPBData
