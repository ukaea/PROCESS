from dataclasses import dataclass


@dataclass(slots=True)
class RebcoData:
    dx_tf_hts_tape_rebco: float = 1.0e-6
    """thickness of REBCO layer in tape (m) (`iteration variable 138`)"""

    dx_hts_tape_copper: float = 100.0e-6
    """thickness of copper layer in tape (m) (`iteration variable 139`)"""

    dx_hts_tape_hastelloy: float = 50.0e-6
    """thickness of Hastelloy layer in tape (m)"""

    dr_hts_tape: float = 4.0e-3
    """Mean width of tape (m)"""

    dx_hts_tape_total: float = 6.5e-5
    """thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)"""

    dia_croco_strand_tape_region: float = 0.0
    """Inner diameter of CroCo strand tape region (m)"""

    dx_croco_strand_copper: float = 2.5e-3
    """Thickness of CroCo strand copper tube (m) (`iteration variable 158`)"""

    copper_rrr: float = 100.0
    """residual resistivity ratio copper in TF superconducting cable"""

    coppera_m2: float = 0.0
    """TF coil current / copper area (A/m2)"""

    coppera_m2_max: float = 1.0e8
    """Maximum TF coil current / copper area (A/m2)"""

    copperaoh_m2: float = 0.0
    """CS coil current / copper area (A/m2) (`sweep variable 61`)"""

    copperaoh_m2_max: float = 1.0e8
    """Maximum CS coil current / copper area (A/m2)"""

    dx_croco_strand_tape_stack: float = 0.0
    """Width / thickness of tape stack in CroCo strand (m)"""

    n_croco_strand_hts_tapes: float = 0.0
    """Number of HTS tapes in CroCo strand"""

    a_croco_strand_rebco: float = 0.0
    """Area of REBCO in CroCo strand (m2)"""

    a_croco_strand_copper_total: float = 0.0
    """Area of copper in CroCo strand (includes tapes and outer tube) (m2)"""

    a_croco_strand_hastelloy: float = 0.0
    """Area of Hastelloy in CroCo strand (m2)"""

    a_croco_strand_solder: float = 0.0
    """Area of solder in CroCo strand (m2)"""

    a_croco_strand: float = 0.0
    """Total area of a CroCo strand (m2)"""


CREATE_DICTS_FROM_DATACLASS = RebcoData
