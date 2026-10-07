"""Module containing variables for the REBCO models"""

from dataclasses import dataclass

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class RebcoData(PROCESSModelData):
    """Dataclass holding REBCO variables"""

    dx_hts_tape_rebco: Parameter[float] = 1.0e-6
    """thickness of REBCO layer in tape (m) (`iteration variable 138`)"""

    dx_hts_tape_copper: Parameter[float] = 100.0e-6
    """thickness of copper layer in tape (m) (`iteration variable 139`)"""

    dx_hts_tape_hastelloy: Parameter[float] = 50.0e-6
    """thickness of Hastelloy layer in tape (m)"""

    dr_hts_tape: Parameter[float] = 4.0e-3
    """Mean width of tape (m)"""

    dx_hts_tape_total: Parameter[float] = 6.5e-5
    """thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)"""

    dia_croco_strand_tape_region: Parameter[float] = 0.0
    """Inner diameter of CroCo strand tape region (m)"""

    dx_croco_strand_copper: Parameter[float] = 2.5e-3
    """Thickness of CroCo strand copper tube (m) (`iteration variable 158`)"""

    copper_rrr: Parameter[float] = 100.0
    """residual resistivity ratio copper in TF superconducting cable"""

    coppera_m2: Parameter[float] = 0.0
    """TF coil current / copper area (A/m2)"""

    coppera_m2_max: Parameter[float] = 1.0e8
    """Maximum TF coil current / copper area (A/m2)"""

    copperaoh_m2: Parameter[float] = 0.0
    """CS coil current / copper area (A/m2) (`sweep variable 61`)"""

    copperaoh_m2_max: Parameter[float] = 1.0e8
    """Maximum CS coil current / copper area (A/m2)"""

    dx_croco_strand_tape_stack: Parameter[float] = 0.0
    """Width / thickness of tape stack in CroCo strand (m)"""

    n_croco_strand_hts_tapes: Parameter[float] = 0.0
    """Number of HTS tapes in CroCo strand"""

    a_croco_strand_rebco: Parameter[float] = 0.0
    """Area of REBCO in CroCo strand (m2)"""

    a_croco_strand_copper_total: Parameter[float] = 0.0
    """Area of copper in CroCo strand (includes tapes and outer tube) (m2)"""

    a_croco_strand_hastelloy: Parameter[float] = 0.0
    """Area of Hastelloy in CroCo strand (m2)"""

    a_croco_strand_solder: Parameter[float] = 0.0
    """Area of solder in CroCo strand (m2)"""

    a_croco_strand: Parameter[float] = 0.0
    """Total area of a CroCo strand (m2)"""


CREATE_DICTS_FROM_DATACLASS = RebcoData
