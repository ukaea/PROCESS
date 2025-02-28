from dataclasses import dataclass
from typing import Any

import process.fortran as fortran


@dataclass
class IterationVariable:
    name: str
    """The name of the variable"""
    module: Any
    """The Fortran module that this variable should be set on."""
    lower_bound: float
    """The default lower bound of the iteration variable"""
    upper_bound: float
    """The default upper bound of the iteration variable"""


ITERATION_VARIABLES = {
    1: IterationVariable("aspect", fortran.constants, 1.100e0, 10.00e0),
    2: IterationVariable("bt", fortran.build_variables, 0.010e0, 30.00e0),
    3: IterationVariable("rmajor", fortran.constraint_variables, 0.100e0, 50.00e0),
    4: IterationVariable("te", fortran.constraint_variables, 5.000e0, 150.0e0),
    5: IterationVariable("beta", fortran.tfcoil_variables, 0.001e0, 1.000e0),
    6: IterationVariable("dene", fortran.tfcoil_variables, 2.00e19, 1.00e21),
    7: IterationVariable(
        "f_nd_beam_electron", fortran.constraint_variables, 1.00e-6, 1.000e0
    ),
    8: IterationVariable(
        "fbeta_poloidal_eps", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    9: IterationVariable("fdene", fortran.fwbs_variables, 0.001e0, 1.000e0),
    10: IterationVariable("hfact", fortran.physics_variables, 0.100e0, 3.000e0),
    11: IterationVariable("pheat", fortran.physics_variables, 1.00e-3, 1.000e3),
    12: IterationVariable(
        "oacdcp", fortran.impurity_radiation_module, 1.000e5, 1.500e8
    ),
    13: IterationVariable("dr_tf_inboard", fortran.rebco_variables, 0.100e0, 5.000e0),
    14: IterationVariable("fwalld", fortran.pfcoil_variables, 0.001e0, 1.000e0),
    15: IterationVariable("fvs", fortran.tfcoil_variables, 0.001e0, 10.000),
    16: IterationVariable("dr_cs", fortran.stellarator_variables, 0.010e0, 10.00e0),
    17: IterationVariable("t_between_pulse", fortran.constants, 0.100e0, 1.000e8),
    18: IterationVariable("q95", fortran.physics_variables, 2.000e0, 50.00e0),
    19: IterationVariable(
        "beam_energy", fortran.current_drive_variables, 1.000e0, 1.000e6
    ),
    20: IterationVariable(
        "temp_cp_average", fortran.tfcoil_variables, 40.00e0, 3.000e2
    ),
    21: IterationVariable("ft_burn", fortran.constraint_variables, 0.001e0, 1.000e0),
    23: IterationVariable("fcoolcp", fortran.tfcoil_variables, 0.100e0, 0.500e0),
    25: IterationVariable("fpnetel", fortran.constraint_variables, 0.001e0, 1.000e0),
    26: IterationVariable("ffuspow", fortran.constraint_variables, 0.001e0, 1.000e0),
    27: IterationVariable("fhldiv", fortran.constraint_variables, 0.001e0, 1.000e0),
    28: IterationVariable("fradpwr", fortran.constraint_variables, 0.001e0, 0.990e0),
    29: IterationVariable("dr_bore", fortran.build_variables, 0.100e0, 10.00e0),
    30: IterationVariable("fmva", fortran.constraint_variables, 0.010e0, 1.000e0),
    31: IterationVariable("gapomin", fortran.build_variables, 0.001e0, 1.000e1),
    32: IterationVariable("frminor", fortran.constraint_variables, 0.001e0, 1.000e0),
    33: IterationVariable("fportsz", fortran.constraint_variables, 0.001e0, 1.000e0),
    34: IterationVariable("fdivcol", fortran.constraint_variables, 0.001e0, 1.000e0),
    35: IterationVariable("fpeakb", fortran.constraint_variables, 0.001e0, 1.000e0),
    36: IterationVariable("fbeta_max", fortran.constraint_variables, 0.001e0, 1.000e0),
    37: IterationVariable("coheof", fortran.pfcoil_variables, 1.000e5, 1.000e8),
    38: IterationVariable("fjohc", fortran.constraint_variables, 0.010e0, 1.000e0),
    39: IterationVariable("fjohc0", fortran.constraint_variables, 0.001e0, 1.000e0),
    40: IterationVariable("fgamcd", fortran.constraint_variables, 0.001e0, 1.000e0),
    41: IterationVariable("fcohbop", fortran.pfcoil_variables, 0.001e0, 1.000e0),
    42: IterationVariable("dr_cs_tf_gap", fortran.build_variables, 0.001e0, 10.00e0),
    44: IterationVariable("fvsbrnni", fortran.physics_variables, 0.001e0, 1.000e0),
    45: IterationVariable("fqval", fortran.constraint_variables, 0.001e0, 1.000e0),
    46: IterationVariable("fpinj", fortran.constraint_variables, 0.001e0, 1.000e0),
    47: IterationVariable("feffcd", fortran.current_drive_variables, 0.001e0, 1.000e0),
    48: IterationVariable("fstrcase", fortran.constraint_variables, 0.001e0, 1.000e0),
    49: IterationVariable("fstrcond", fortran.constraint_variables, 0.001e0, 1.000e0),
    50: IterationVariable("fiooic", fortran.constraint_variables, 0.001e0, 1.000e0),
    51: IterationVariable("fvdump", fortran.constraint_variables, 0.001e0, 1.000e0),
    53: IterationVariable("fjprot", fortran.constraint_variables, 0.001e0, 1.000e0),
    54: IterationVariable("ftmargtf", fortran.constraint_variables, 0.001e0, 1.000e0),
    56: IterationVariable("tdmptf", fortran.tfcoil_variables, 0.100e0, 100.0e0),
    57: IterationVariable("thkcas", fortran.error_handling, 0.050e0, 1.000e0),
    58: IterationVariable("thwcndut", fortran.tfcoil_variables, 0.001e0, 0.100e0),
    59: IterationVariable("fcutfsu", fortran.tfcoil_variables, 0.001e0, 1.000e0),
    60: IterationVariable("cpttf", fortran.error_handling, 0.001e0, 4.000e4),
    61: IterationVariable(
        "dr_shld_vv_gap_inboard", fortran.build_variables, 0.001e0, 10.00e0
    ),
    62: IterationVariable("fdtmp", fortran.constraint_variables, 0.001e0, 1.000e0),
    63: IterationVariable("ftpeak", fortran.constraint_variables, 0.001e0, 1.000e0),
    64: IterationVariable("fauxmn", fortran.constraint_variables, 0.001e0, 1.000e0),
    65: IterationVariable(
        "t_current_ramp_up", fortran.error_handling, 0.100e0, 1.000e3
    ),
    66: IterationVariable(
        "ft_current_ramp_up", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    67: IterationVariable("ftcycl", fortran.constraint_variables, 0.001e0, 1.000e0),
    68: IterationVariable("fptemp", fortran.constraint_variables, 0.001e0, 1.000e0),
    69: IterationVariable("rcool", fortran.tfcoil_variables, 0.001e0, 0.010e0),
    70: IterationVariable("vcool", fortran.tfcoil_variables, 1.000e0, 1.000e2),
    71: IterationVariable("fq", fortran.constraint_variables, 0.001e0, 1.000e0),
    72: IterationVariable("fipir", fortran.constraint_variables, 0.001e0, 1.000e0),
    73: IterationVariable(
        "dr_fw_plasma_gap_inboard", fortran.build_variables, 0.001e0, 10.00e0
    ),
    74: IterationVariable(
        "dr_fw_plasma_gap_outboard", fortran.build_variables, 0.001e0, 10.00e0
    ),
    75: IterationVariable("tfootfi", fortran.build_variables, 0.200e0, 5.000e0),
    79: IterationVariable(
        "fbeta_poloidal", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    81: IterationVariable("edrive", fortran.ife_variables, 1.000e5, 5.000e7),
    82: IterationVariable("drveff", fortran.ife_variables, 0.010e0, 1.000e0),
    83: IterationVariable("tgain", fortran.ife_variables, 1.000e0, 500.0e0),
    84: IterationVariable("chrad", fortran.ife_variables, 0.100e0, 20.00e0),
    85: IterationVariable("pdrive", fortran.ife_variables, 1.000e6, 200.0e6),
    86: IterationVariable("frrmax", fortran.ife_variables, 0.001e0, 1.000e0),
    89: IterationVariable("ftbr", fortran.constraint_variables, 0.001e0, 1.000e0),
    90: IterationVariable("blbuith", fortran.build_variables, 0.001e0, 2.000e0),
    91: IterationVariable("blbuoth", fortran.build_variables, 0.001e0, 2.000e0),
    92: IterationVariable("fflutf", fortran.constraint_variables, 0.001e0, 1.000e0),
    93: IterationVariable("dr_shld_inboard", fortran.build_variables, 0.001e0, 10.00e0),
    94: IterationVariable(
        "dr_shld_outboard", fortran.build_variables, 0.001e0, 10.00e0
    ),
    95: IterationVariable("fptfnuc", fortran.constraint_variables, 0.001e0, 1.000e0),
    96: IterationVariable("fvvhe", fortran.constraint_variables, 0.001e0, 1.000e0),
    97: IterationVariable("fpsepr", fortran.constraint_variables, 0.001e0, 1.000e0),
    98: IterationVariable("li6enrich", fortran.fwbs_variables, 10.00e0, 100.0e0),
    103: IterationVariable(
        "fl_h_threshold", fortran.constraint_variables, 1.000e0, 1.000e6
    ),
    104: IterationVariable("fcwr", fortran.constraint_variables, 0.001e0, 1.000e0),
    105: IterationVariable("fnbshinef", fortran.constraint_variables, 0.001e0, 1.000e0),
    106: IterationVariable("ftmargoh", fortran.constraint_variables, 0.001e0, 1.000e0),
    107: IterationVariable("favail", fortran.cost_variables, 0.001e0, 1.000e0),
    108: IterationVariable("breeder_f", fortran.fwbs_variables, 0.060e0, 1.000e0),
    109: IterationVariable(
        "f_nd_alpha_electron", fortran.physics_variables, 0.050e0, 0.150e0
    ),
    110: IterationVariable(
        "falpha_energy_confinement", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    111: IterationVariable(
        "fniterpump", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    112: IterationVariable("fzeffmax", fortran.constraint_variables, 0.001e0, 1.000e0),
    113: IterationVariable(
        "fmaxvvstress", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    114: IterationVariable("len_fw_channel", fortran.fwbs_variables, 0.001e0, 1.000e3),
    115: IterationVariable(
        "fpoloidalpower", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    116: IterationVariable("fradwall", fortran.constraint_variables, 0.001e0, 1.000e0),
    117: IterationVariable("fpsepbqar", fortran.constraint_variables, 0.001e0, 1.000e0),
    118: IterationVariable("fpsep", fortran.constraint_variables, 0.001e0, 1.000e0),
    119: IterationVariable("tesep", fortran.physics_variables, 0.000e0, 1.000e1),
    120: IterationVariable("ttarget", fortran.numerics, 1.000e0, 1.000e4),
    121: IterationVariable("neratio", fortran.numerics, 0.001e0, 1.000e0),
    122: IterationVariable("oh_steel_frac", fortran.pfcoil_variables, 0.001e0, 0.950e0),
    123: IterationVariable(
        "foh_stress", fortran.constraint_variables, 0.001e0, 1.000e0
    ),
    124: IterationVariable("qtargettotal", fortran.numerics, 0.001e0, 1.000e7),
    137: IterationVariable("fplhsep", fortran.physics_variables, 0.001e0, 1.000e0),
    138: IterationVariable(
        "rebco_thickness", fortran.physics_variables, 0.01e-6, 100.0e-6
    ),
    139: IterationVariable("copper_thick", fortran.rebco_variables, 1.00e-6, 1.00e-3),
    140: IterationVariable("dr_tf_wp", fortran.tfcoil_variables, 0.001e0, 2.000e0),
    141: IterationVariable("fcqt", fortran.constraint_variables, 0.001e0, 1.000e0),
    142: IterationVariable("nesep", fortran.physics_variables, 1.00e17, 1.00e20),
    143: IterationVariable("f_coppera_m2", fortran.rebco_variables, 0.001e0, 1.000e0),
    144: IterationVariable("fnesep", fortran.constraint_variables, 0.001e0, 1.000e0),
    145: IterationVariable("fgwped", fortran.physics_variables, 0.100e0, 0.9e0),
    146: IterationVariable("fcpttf", fortran.constraint_variables, 0.001e0, 1.000e0),
    147: IterationVariable("freinke", fortran.constraint_variables, 0.001e0, 1.000e0),
    148: IterationVariable("fzactual", fortran.numerics, 1.00e-8, 1.000e0),
    149: IterationVariable("fbmaxcs", fortran.pfcoil_variables, 0.001e0, 1.000e0),
    152: IterationVariable("fgwsep", fortran.physics_variables, 0.001e0, 0.5e0),
    153: IterationVariable("fpdivlim", fortran.physics_variables, 0.001e0, 1.000e0),
    154: IterationVariable("fne0", fortran.physics_variables, 0.001e0, 1.000e0),
    155: IterationVariable("pfusife", fortran.ife_variables, 5.000e2, 3.000e3),
    156: IterationVariable("rrin", fortran.ife_variables, 1.000e0, 1.000e1),
    157: IterationVariable("fvssu", fortran.pfcoil_variables, 1.00e-3, 1.000e1),
    158: IterationVariable("croco_thick", fortran.rebco_variables, 1.0e-3, 1.0e-1),
    159: IterationVariable("ftoroidalgap", fortran.tfcoil_variables, 1.0e-4, 1.0e0),
    160: IterationVariable("f_avspace", fortran.build_variables, 0.010e0, 1.000e0),
    161: IterationVariable("fbeta_min", fortran.constraint_variables, 0.010e0, 1.000e0),
    162: IterationVariable("r_cp_top", fortran.build_variables, 0.0010e0, 10.000e0),
    163: IterationVariable("f_t_turn_tf", fortran.tfcoil_variables, 0.0010e0, 1000.0e0),
    164: IterationVariable(
        "f_crypmw", fortran.heat_transport_variables, 0.001e0, 1.000e0
    ),
    165: IterationVariable("fstr_wp", fortran.constraint_variables, 1.0e-9, 1.0e0),
    166: IterationVariable("f_copperaoh_m2", fortran.rebco_variables, 0.001e0, 1.000e0),
    167: IterationVariable("fncycle", fortran.constraint_variables, 1.0e-8, 1.0e0),
    168: IterationVariable(
        "fecrh_ignition", fortran.constraint_variables, 0.010e0, 2.000e0
    ),
    169: IterationVariable(
        "te0_ecrh_achievable", fortran.stellarator_variables, 1.0e0, 1.0e3
    ),
    170: IterationVariable("beta_div", fortran.divertor_variables, 0.49, 5.01),
    171: IterationVariable("casths_fraction", fortran.tfcoil_variables, 0.01, 0.99),
    172: IterationVariable("casths", fortran.tfcoil_variables, 0.001, 1.0),
    173: IterationVariable("f_tritium", fortran.physics_variables, 0.000, 1.000),
}

# TODO: 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136


def loadxc():
    """Loads the physics and engineering variables into the optimisation variable array."""

    for i in range(fortran.numerics.nvar):
        variable_index = fortran.ixc[i].item()
        iteration_variable = ITERATION_VARIABLES[variable_index]
        fortran.numerics.xcm[i] = getattr(
            iteration_variable.module, iteration_variable.name
        )
