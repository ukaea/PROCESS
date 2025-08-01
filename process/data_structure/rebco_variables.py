rebco_thickness: float = None
"""thickness of REBCO layer in tape (m) (`iteration variable 138`)"""

copper_thick: float = None
"""thickness of copper layer in tape (m) (`iteration variable 139`)"""

hastelloy_thickness: float = None
"""thickness of Hastelloy layer in tape (m)"""

tape_width: float = None
"""Mean width of tape (m)"""

tape_thickness: float = None
"""thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)"""

croco_od: float = None
"""Outer diameter of CroCo strand (m)"""

croco_id: float = None
"""Inner diameter of CroCo copper tube (m)"""

croco_thick: float = None
"""Thickness of CroCo copper tube (m) (`iteration variable 158`)"""

copper_rrr: float = None
"""residual resistivity ratio copper in TF superconducting cable"""

copperA_m2: float = None  # noqa: N816
"""TF coil current / copper area (A/m2)"""

coppera_m2_max: float = None
"""Maximum TF coil current / copper area (A/m2)"""

f_coppera_m2: float = None
"""f-value for constraint 75: TF coil current / copper area < copperA_m2_max"""

copperaoh_m2: float = None
"""CS coil current / copper area (A/m2) (`sweep variable 61`)"""

copperaoh_m2_max: float = None
"""Maximum CS coil current / copper area (A/m2)"""

f_copperaoh_m2: float = None
"""f-value for constraint 88: CS coil current / copper area < copperA_m2_max"""

stack_thickness: float = None

tapes: float = None

rebco_area: float = None

copper_area: float = None

hastelloy_area: float = None

solder_area: float = None

croco_area: float = None


def init_rebco_variables():
    """Initialise the REBCO variables"""
    global rebco_thickness
    global copper_thick
    global hastelloy_thickness
    global tape_width
    global croco_od
    global croco_id
    global croco_thick
    global copper_rrr
    global coppera_m2_max
    global f_coppera_m2
    global tape_thickness
    global stack_thickness
    global tapes
    global rebco_area
    global copper_area
    global hastelloy_area
    global solder_area
    global croco_area
    global copperA_m2
    global copperaoh_m2_max
    global f_copperaoh_m2
    global copperaoh_m2

    rebco_thickness = 1.0e-6
    copper_thick = 100.0e-6
    hastelloy_thickness = 50.0e-6
    tape_width = 4.0e-3
    croco_od = 0.0
    croco_id = 0.0
    croco_thick = 2.5e-3
    copper_rrr = 100.0
    coppera_m2_max = 1.0e8
    f_coppera_m2 = 1.0
    tape_thickness = 6.5e-5
    stack_thickness = 0.0
    tapes = 0.0
    rebco_area = 0.0
    copper_area = 0.0
    hastelloy_area = 0.0
    solder_area = 0.0
    croco_area = 0.0
    copperA_m2 = 0.0
    copperaoh_m2_max = 1.0e8
    f_copperaoh_m2 = 1.0
    copperaoh_m2 = 0.0
