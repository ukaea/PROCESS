dx_hts_tape_rebco: float = None
"""thickness of REBCO layer in tape (m) (`iteration variable 138`)"""

dx_hts_tape_copper: float = None
"""thickness of copper layer in tape (m) (`iteration variable 139`)"""

dx_hts_tape_hastelloy: float = None
"""thickness of Hastelloy layer in tape (m)"""

dr_hts_tape: float = None
"""Mean width of tape (m)"""

dx_hts_tape_total: float = None
"""thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)"""

dia_croco_strand: float = None
"""Outer diameter of CroCo strand (m)"""

croco_id: float = None
"""Inner diameter of CroCo copper tube (m)"""

dx_croco_strand_copper: float = None
"""Thickness of CroCo strand copper tube (m) (`iteration variable 158`)"""

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

dx_croco_strand_tape_stack: float = None
"Width / thickness of tape stack in CroCo strand (m)"

n_croco_strand_hts_tapes: float = None
"Number of HTS tapes in CroCo strand"

a_croco_strand_hts_tapes: float = None
"Area of HTS tapes in CroCo strand (m2)"

copper_area: float = None

a_croco_strand_hastelloy: float = None
"Area of Hastelloy in CroCo strand (m2)"

a_croco_strand_solder: float = None
"Area of solder in CroCo strand (m2)"

croco_area: float = None


def init_rebco_variables():
    """Initialise the REBCO variables"""
    global dx_hts_tape_rebco
    global dx_hts_tape_copper
    global dx_hts_tape_hastelloy
    global dr_hts_tape
    global dia_croco_strand
    global croco_id
    global dx_croco_strand_copper
    global copper_rrr
    global coppera_m2_max
    global f_coppera_m2
    global dx_hts_tape_total
    global dx_croco_strand_tape_stack
    global n_croco_strand_hts_tapes
    global a_croco_strand_hts_tapes
    global copper_area
    global a_croco_strand_hastelloy
    global a_croco_strand_solder
    global croco_area
    global copperA_m2
    global copperaoh_m2_max
    global f_copperaoh_m2
    global copperaoh_m2

    dx_hts_tape_rebco = 1.0e-6
    dx_hts_tape_copper = 100.0e-6
    dx_hts_tape_hastelloy = 50.0e-6
    dr_hts_tape = 4.0e-3
    dia_croco_strand = 0.0
    croco_id = 0.0
    dx_croco_strand_copper = 2.5e-3
    copper_rrr = 100.0
    coppera_m2_max = 1.0e8
    f_coppera_m2 = 1.0
    dx_hts_tape_total = 6.5e-5
    dx_croco_strand_tape_stack = 0.0
    n_croco_strand_hts_tapes = 0.0
    a_croco_strand_hts_tapes = 0.0
    copper_area = 0.0
    a_croco_strand_hastelloy = 0.0
    a_croco_strand_solder = 0.0
    croco_area = 0.0
    copperA_m2 = 0.0
    copperaoh_m2_max = 1.0e8
    f_copperaoh_m2 = 1.0
    copperaoh_m2 = 0.0
