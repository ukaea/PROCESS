"""author: J Morris (UKAEA)
This module contains the PROCESS CCFE HCPB blanket model
based on CCFE HCPB model from the PROCESS engineering paper
PROCESS Engineering paper (M. Kovari et al.)
### References
- Kovari et al., Fusion Engineering and Design 104 (2016) 9-20
"""

# Smeared densities of build sections

armour_density: float = None
"""FW armour density [kg/m3]"""


fw_density: float = None
"""FW density [kg/m3]"""


blanket_density: float = None
"""Blanket density [kg/m3]"""


shield_density: float = None
"""Shield density [kg/m3]"""


vv_density: float = None
"""Vacuum vessel density [kg/m3]"""


x_blanket: float = None
"""Blanket exponent (tonne/m2)"""


x_shield: float = None
"""Shield exponent (tonne/m2)"""


tfc_nuc_heating: float = None
"""Unit nuclear heating in TF coil (W per W of fusion power)"""


fw_armour_u_nuc_heating: float = None
"""Unit heating of FW and armour in FW armour (W/kg per W of fusion power)"""


shld_u_nuc_heating: float = None
"""Unit nuclear heating in shield (W per W of fusion power)"""


pnuc_tot_blk_sector: float = None
"""Total nuclear power deposited in blanket covered sector (FW, BLKT, SHLD, TF) (MW)"""


exp_blanket: float = None
"""Exponential factors in nuclear heating calcs"""


exp_shield1: float = None
"""Exponential factors in nuclear heating calcs"""


exp_shield2: float = None
"""Exponential factors in nuclear heating calcs"""


def init_ccfe_hcpb_module():
    global iso_fortran_env
    global armour_density
    global fw_density
    global blanket_density
    global shield_density
    global vv_density
    global x_blanket
    global x_shield
    global tfc_nuc_heating
    global fw_armour_u_nuc_heating
    global shld_u_nuc_heating
    global pnuc_tot_blk_sector
    global exp_blanket
    global exp_shield1
    global exp_shield2

    armour_density = 0.0
    fw_density = 0.0
    blanket_density = 0.0
    shield_density = 0.0
    vv_density = 0.0
    x_blanket = 0.0
    x_shield = 0.0
    tfc_nuc_heating = 0.0
    fw_armour_u_nuc_heating = 6.25e-7
    shld_u_nuc_heating = 0.0
    exp_blanket = 0.0
    exp_shield1 = 0.0
    exp_shield2 = 0.0
