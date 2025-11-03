from copy import deepcopy
from dataclasses import dataclass
from typing import Any
from warnings import warn

import numpy as np

import process.data_structure as data_structure
import process.fortran as fortran
from process.exceptions import ProcessValueError
from process.utilities.f2py_string_patch import (
    f2py_compatible_to_string,
    string_to_f2py_compatible,
)


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
    target_name: str | None = None
    """If specified, the iteration variable is set as `module.target_name`
    rather than `module.name`.
    """
    array_index: int | None = None
    """If `module.name` is an array, the iteration variable can only modify
    `array_index` of that array.

    NOTE: The indexes start at 0 (despite indexing Fortran arrays).
    """


ITERATION_VARIABLES = {
    1: IterationVariable("aspect", fortran.physics_variables, 1.1, 10.00),
    2: IterationVariable("bt", fortran.physics_variables, 0.010, 30.00),
    3: IterationVariable("rmajor", fortran.physics_variables, 0.1, 50.00),
    4: IterationVariable("te", fortran.physics_variables, 5.0, 150.0),
    5: IterationVariable("beta", fortran.physics_variables, 0.001, 1.0),
    6: IterationVariable("dene", fortran.physics_variables, 2.0e19, 1.0e21),
    7: IterationVariable("f_nd_beam_electron", fortran.physics_variables, 1.0e-6, 1.0),
    8: IterationVariable(
        "fbeta_poloidal_eps", fortran.constraint_variables, 0.001, 1.0
    ),
    9: IterationVariable("fdene", fortran.constraint_variables, 0.001, 1.0),
    10: IterationVariable("hfact", fortran.physics_variables, 0.1, 3.0),
    11: IterationVariable(
        "p_hcd_primary_extra_heat_mw", fortran.current_drive_variables, 1.0e-3, 1.0e3
    ),
    12: IterationVariable("oacdcp", fortran.tfcoil_variables, 1.0e5, 1.50e8),
    13: IterationVariable("dr_tf_inboard", fortran.build_variables, 0.1, 5.0),
    14: IterationVariable(
        "fpflux_fw_neutron_max_mw", fortran.constraint_variables, 0.001, 1.0
    ),
    15: IterationVariable(
        "fvs_plasma_total_required", fortran.constraint_variables, 0.001, 10.000
    ),
    16: IterationVariable("dr_cs", fortran.build_variables, 0.01, 10.00),
    17: IterationVariable(
        "t_between_pulse", data_structure.times_variables, 0.1, 1.0e8
    ),
    18: IterationVariable("q95", fortran.physics_variables, 2.0, 50.00),
    19: IterationVariable("e_beam_kev", fortran.current_drive_variables, 1.0, 1.0e6),
    20: IterationVariable("temp_cp_average", fortran.tfcoil_variables, 40.00, 3.0e2),
    21: IterationVariable("ft_burn_min", fortran.constraint_variables, 0.001, 1.0),
    23: IterationVariable("fcoolcp", fortran.tfcoil_variables, 0.1, 0.50),
    25: IterationVariable(
        "fp_plant_electric_net_required_mw", fortran.constraint_variables, 0.001, 1.0
    ),
    26: IterationVariable(
        "fp_fusion_total_max_mw", fortran.constraint_variables, 0.001, 1.0
    ),
    27: IterationVariable(
        "fpflux_div_heat_load_mw", fortran.constraint_variables, 0.001, 1.0
    ),
    28: IterationVariable("fradpwr", fortran.constraint_variables, 0.001, 0.99),
    29: IterationVariable("dr_bore", fortran.build_variables, 0.1, 10.00),
    30: IterationVariable("fmva", fortran.constraint_variables, 0.010, 1.0),
    31: IterationVariable("gapomin", fortran.build_variables, 0.001, 1.0e1),
    32: IterationVariable("frminor", fortran.constraint_variables, 0.001, 1.0),
    33: IterationVariable(
        "fradius_beam_tangency", fortran.constraint_variables, 0.001, 1.0
    ),
    35: IterationVariable(
        "fb_tf_inboard_max", fortran.constraint_variables, 0.001, 1.0
    ),
    36: IterationVariable("fbeta_max", fortran.constraint_variables, 0.001, 1.0),
    37: IterationVariable("j_cs_flat_top_end", fortran.pfcoil_variables, 1.0e5, 1.0e8),
    38: IterationVariable("fjohc", fortran.constraint_variables, 0.010, 1.0),
    39: IterationVariable("fjohc0", fortran.constraint_variables, 0.001, 1.0),
    40: IterationVariable(
        "feta_cd_norm_hcd_primary_max", fortran.constraint_variables, 0.001, 1.0
    ),
    41: IterationVariable(
        "f_j_cs_start_pulse_end_flat_top", fortran.pfcoil_variables, 0.001, 1.0
    ),
    42: IterationVariable("dr_cs_tf_gap", fortran.build_variables, 0.001, 10.00),
    44: IterationVariable(
        "f_c_plasma_non_inductive", fortran.physics_variables, 0.001, 1.0
    ),
    45: IterationVariable("fqval", fortran.constraint_variables, 0.001, 1.0),
    46: IterationVariable(
        "fp_hcd_injected_max", fortran.constraint_variables, 0.001, 1.0
    ),
    47: IterationVariable("feffcd", fortran.current_drive_variables, 0.001, 1.0),
    48: IterationVariable("fstrcase", fortran.constraint_variables, 0.001, 1.0),
    49: IterationVariable("fstrcond", fortran.constraint_variables, 0.001, 1.0),
    50: IterationVariable("fiooic", fortran.constraint_variables, 0.001, 1.0),
    51: IterationVariable("fvdump", fortran.constraint_variables, 0.001, 1.0),
    53: IterationVariable("fjprot", fortran.constraint_variables, 0.001, 1.0),
    54: IterationVariable("ftmargtf", fortran.constraint_variables, 0.001, 1.0),
    56: IterationVariable("tdmptf", fortran.tfcoil_variables, 0.1, 100.0),
    57: IterationVariable("dr_tf_nose_case", fortran.tfcoil_variables, 0.05, 1.0),
    58: IterationVariable("dx_tf_turn_steel", fortran.tfcoil_variables, 0.001, 0.1),
    59: IterationVariable("fcutfsu", fortran.tfcoil_variables, 0.001, 1.0),
    60: IterationVariable("c_tf_turn", fortran.tfcoil_variables, 0.001, 4.0e4),
    61: IterationVariable(
        "dr_shld_vv_gap_inboard", fortran.build_variables, 0.001, 10.00
    ),
    62: IterationVariable("fdtmp", fortran.constraint_variables, 0.001, 1.0),
    63: IterationVariable("ftemp_fw_max", fortran.constraint_variables, 0.001, 1.0),
    64: IterationVariable("fauxmn", fortran.constraint_variables, 0.001, 1.0),
    65: IterationVariable(
        "t_current_ramp_up", data_structure.times_variables, 0.1, 1.0e3
    ),
    66: IterationVariable(
        "ft_current_ramp_up", fortran.constraint_variables, 0.001, 1.0
    ),
    67: IterationVariable("ft_cycle_min", fortran.constraint_variables, 0.001, 1.0),
    68: IterationVariable("fptemp", fortran.constraint_variables, 0.001, 1.0),
    69: IterationVariable("rcool", fortran.tfcoil_variables, 0.001, 0.010),
    70: IterationVariable("vcool", fortran.tfcoil_variables, 1.0, 1.0e2),
    71: IterationVariable("fq", fortran.constraint_variables, 0.001, 1.0),
    72: IterationVariable("fipir", fortran.constraint_variables, 0.001, 1.0),
    73: IterationVariable(
        "dr_fw_plasma_gap_inboard", fortran.build_variables, 0.001, 10.00
    ),
    74: IterationVariable(
        "dr_fw_plasma_gap_outboard", fortran.build_variables, 0.001, 10.00
    ),
    75: IterationVariable("tfootfi", fortran.build_variables, 0.200, 5.0),
    79: IterationVariable("fbeta_poloidal", fortran.constraint_variables, 0.001, 1.0),
    81: IterationVariable("edrive", fortran.ife_variables, 1.0e5, 5.0e7),
    82: IterationVariable("drveff", fortran.ife_variables, 0.010, 1.0),
    83: IterationVariable("tgain", fortran.ife_variables, 1.0, 500.0),
    84: IterationVariable("chrad", fortran.ife_variables, 0.1, 20.00),
    85: IterationVariable("pdrive", fortran.ife_variables, 1.0e6, 200.0e6),
    86: IterationVariable("frrmax", fortran.ife_variables, 0.001, 1.0),
    89: IterationVariable("ftbr", fortran.constraint_variables, 0.001, 1.0),
    90: IterationVariable("blbuith", fortran.build_variables, 0.001, 2.0),
    91: IterationVariable("blbuoth", fortran.build_variables, 0.001, 2.0),
    92: IterationVariable("fflutf", fortran.constraint_variables, 0.001, 1.0),
    93: IterationVariable("dr_shld_inboard", fortran.build_variables, 0.001, 10.00),
    94: IterationVariable("dr_shld_outboard", fortran.build_variables, 0.001, 10.00),
    95: IterationVariable("fptfnuc", fortran.constraint_variables, 0.001, 1.0),
    96: IterationVariable("fvvhe", fortran.constraint_variables, 0.001, 1.0),
    97: IterationVariable("fpsepr", fortran.constraint_variables, 0.001, 1.0),
    98: IterationVariable(
        "f_blkt_li6_enrichment", fortran.fwbs_variables, 10.00, 100.0
    ),
    103: IterationVariable("fl_h_threshold", fortran.constraint_variables, 0.001, 1.0),
    104: IterationVariable("fcwr", fortran.constraint_variables, 0.001, 1.0),
    105: IterationVariable("fnbshinef", fortran.constraint_variables, 0.001, 1.0),
    106: IterationVariable("ftmargoh", fortran.constraint_variables, 0.001, 1.0),
    107: IterationVariable("favail", data_structure.cost_variables, 0.001, 1.0),
    108: IterationVariable("breeder_f", fortran.fwbs_variables, 0.060, 1.0),
    109: IterationVariable(
        "f_nd_alpha_electron", fortran.physics_variables, 0.05, 0.15
    ),
    110: IterationVariable(
        "falpha_energy_confinement", fortran.constraint_variables, 0.001, 1.0
    ),
    111: IterationVariable("fniterpump", fortran.constraint_variables, 0.001, 1.0),
    112: IterationVariable("fzeffmax", fortran.constraint_variables, 0.001, 1.0),
    113: IterationVariable("fmaxvvstress", fortran.constraint_variables, 0.001, 1.0),
    114: IterationVariable("len_fw_channel", fortran.fwbs_variables, 0.001, 1.0e3),
    115: IterationVariable("fpoloidalpower", fortran.constraint_variables, 0.001, 1.0),
    116: IterationVariable(
        "fpflux_fw_rad_max", fortran.constraint_variables, 0.001, 1.0
    ),
    117: IterationVariable("fpsepbqar", fortran.constraint_variables, 0.001, 1.0),
    118: IterationVariable("fpsep", fortran.constraint_variables, 0.001, 1.0),
    119: IterationVariable("tesep", fortran.physics_variables, 0.0, 1.0e1),
    122: IterationVariable("f_a_cs_steel", fortran.pfcoil_variables, 0.001, 0.950),
    123: IterationVariable("foh_stress", fortran.constraint_variables, 0.001, 1.0),
    125: IterationVariable(
        "fimp(03)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=2,
    ),
    126: IterationVariable(
        "fimp(04)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=3,
    ),
    127: IterationVariable(
        "fimp(05)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=4,
    ),
    128: IterationVariable(
        "fimp(06)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=5,
    ),
    129: IterationVariable(
        "fimp(07)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=6,
    ),
    130: IterationVariable(
        "fimp(08)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=7,
    ),
    131: IterationVariable(
        "fimp(09)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=8,
    ),
    132: IterationVariable(
        "fimp(10)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=9,
    ),
    133: IterationVariable(
        "fimp(11)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=10,
    ),
    134: IterationVariable(
        "fimp(12)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=11,
    ),
    135: IterationVariable(
        "fimp(13)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=12,
    ),
    136: IterationVariable(
        "fimp(14)",
        fortran.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="impurity_arr_frac",
        array_index=13,
    ),
    137: IterationVariable("fplhsep", fortran.physics_variables, 0.001, 1.0),
    138: IterationVariable(
        "rebco_thickness", fortran.physics_variables, 0.01e-6, 100.0e-6
    ),
    139: IterationVariable(
        "copper_thick", data_structure.rebco_variables, 1.0e-6, 1.0e-3
    ),
    140: IterationVariable(
        "dr_tf_wp_with_insulation", fortran.tfcoil_variables, 0.001, 2.0
    ),
    141: IterationVariable("fcqt", fortran.constraint_variables, 0.001, 1.0),
    142: IterationVariable("nesep", fortran.physics_variables, 1.0e17, 1.0e20),
    143: IterationVariable("f_coppera_m2", data_structure.rebco_variables, 0.001, 1.0),
    144: IterationVariable("fnesep", fortran.constraint_variables, 0.001, 1.0),
    145: IterationVariable("fgwped", fortran.physics_variables, 0.1, 0.9),
    146: IterationVariable("fc_tf_turn_max", fortran.constraint_variables, 0.001, 1.0),
    147: IterationVariable("freinke", fortran.constraint_variables, 0.001, 1.0),
    149: IterationVariable("fb_cs_limit_max", fortran.pfcoil_variables, 0.001, 1.0),
    152: IterationVariable("fgwsep", fortran.physics_variables, 0.001, 0.5),
    153: IterationVariable("fpdivlim", fortran.physics_variables, 0.001, 1.0),
    154: IterationVariable("fne0", fortran.physics_variables, 0.001, 1.0),
    155: IterationVariable("pfusife", fortran.ife_variables, 5.0e2, 3.0e3),
    156: IterationVariable("rrin", fortran.ife_variables, 1.0, 1.0e1),
    157: IterationVariable(
        "fvs_cs_pf_total_ramp", fortran.pfcoil_variables, 1.0e-3, 1.0e1
    ),
    158: IterationVariable(
        "croco_thick", data_structure.rebco_variables, 1.0e-3, 1.0e-1
    ),
    159: IterationVariable("ftoroidalgap", fortran.tfcoil_variables, 1.0e-4, 1.0),
    160: IterationVariable("f_avspace", fortran.build_variables, 0.010, 1.0),
    161: IterationVariable("fbeta_min", fortran.constraint_variables, 0.010, 1.0),
    162: IterationVariable("r_cp_top", fortran.build_variables, 0.0010, 10.0),
    163: IterationVariable("f_t_turn_tf", fortran.tfcoil_variables, 0.0010, 1000.0),
    164: IterationVariable("f_crypmw", fortran.heat_transport_variables, 0.001, 1.0),
    165: IterationVariable("fstr_wp", fortran.constraint_variables, 1.0e-9, 1.0),
    166: IterationVariable(
        "f_copperaoh_m2", data_structure.rebco_variables, 0.001, 1.0
    ),
    167: IterationVariable("fncycle", fortran.constraint_variables, 1.0e-8, 1.0),
    168: IterationVariable("fecrh_ignition", fortran.constraint_variables, 0.010, 2.0),
    169: IterationVariable(
        "te0_ecrh_achievable", fortran.stellarator_variables, 1.0, 1.0e3
    ),
    170: IterationVariable(
        "deg_div_field_plate", data_structure.divertor_variables, 0.49, 5.01
    ),
    171: IterationVariable("casths_fraction", fortran.tfcoil_variables, 0.01, 0.99),
    172: IterationVariable("dx_tf_side_case_min", fortran.tfcoil_variables, 0.001, 1.0),
    173: IterationVariable("f_tritium", fortran.physics_variables, 0.000, 1.000),
    174: IterationVariable("triang", fortran.physics_variables, 0.00, 1.00),
    175: IterationVariable("kappa", fortran.physics_variables, 0.00, 10.00),
    176: IterationVariable("f_st_coil_aspect", fortran.stellarator_variables, 0.80, 1.20),
}


def check_iteration_variable(iteration_variable_value, name: str = ""):
    """Check that the iteration variable value is valid (not a weird number or too small).

    Raises an error upon encountering an invalid value, otherwise does nothing.
    """
    if abs(iteration_variable_value) <= 1e-12:
        error_msg = f"Iteration variable {name} is 0 (or very close)"
        raise ProcessValueError(
            error_msg, iteration_variable_value=iteration_variable_value
        )

    if np.isnan(iteration_variable_value) or np.isinf(iteration_variable_value):
        error_msg = f"Iteration variable {name} invalid number"
        raise ProcessValueError(
            error_msg, iteration_variable_value=iteration_variable_value
        )


def load_iteration_variables():
    """Loads the physics and engineering variables into the optimisation variable array."""

    for i in range(fortran.numerics.nvar):
        variable_index = fortran.numerics.ixc[i].item()
        iteration_variable = ITERATION_VARIABLES[variable_index]

        # use ... as the default return value because None might be a valid return from Fortran?
        iteration_variable_value = getattr(
            iteration_variable.module,
            iteration_variable.target_name or iteration_variable.name,
            ...,
        )

        if iteration_variable_value is ...:
            error_msg = (
                f"Could not get the value for iteration variable {variable_index} "
                f"({iteration_variable.name})"
            )
            raise ProcessValueError(error_msg)

        # if an array index is specified, iteration_variable_value is currently
        # the whole array and not just the element we are interested in. Lets extract
        # the correct element
        if iteration_variable.array_index is not None:
            iteration_variable_value = iteration_variable_value[
                iteration_variable.array_index
            ]

        fortran.numerics.xcm[i] = iteration_variable_value
        fortran.numerics.name_xc[i] = string_to_f2py_compatible(
            fortran.numerics.name_xc[i], iteration_variable.name
        )

        # warn of the iteration variable is also a scan variable because this will cause
        # the optimiser and scan to overwrite the same variable and conflict
        if iteration_variable.name in (
            f2py_compatible_to_string(fortran.global_variables.vlabel),
            f2py_compatible_to_string(fortran.global_variables.vlabel_2),
        ):
            warn(
                (
                    "The sweep variable is also an iteration variable and will be "
                    "overwritten by the optimiser"
                ),
                stacklevel=3,
            )

        # check that the iteration variable is valid (not 0, NaN, inf, or very large)

        check_iteration_variable(
            iteration_variable_value,
            name=f"{variable_index} ({iteration_variable.name})",
        )

        fortran.numerics.scale[i] = 1.0 / iteration_variable_value
        fortran.numerics.scafc[i] = 1.0 / fortran.numerics.scale[i]

        fortran.numerics.xcm[i] = iteration_variable_value * fortran.numerics.scale[i]


def set_scaled_iteration_variable(xc, nn: int):
    """Converts scaled iteration variables back to their real values and sets them in the code.

    :param xc: scaled iteration variable values
    :param nn: number of iteration variables
    """

    for i in range(nn):
        # there is less error handling here than in load_iteration_variables
        # because many errors will be caught in load_iteration_variables which is
        # run first. This verifies the variables exist and the module target is correct.
        variable_index = fortran.numerics.ixc[i].item()
        iteration_variable = ITERATION_VARIABLES[variable_index]

        ratio = xc[i] / fortran.numerics.scale[i]

        if iteration_variable.array_index is None:
            setattr(
                iteration_variable.module,
                iteration_variable.target_name or iteration_variable.name,
                ratio,
            )
        else:
            current_array = getattr(
                iteration_variable.module,
                iteration_variable.target_name or iteration_variable.name,
            )
            new_array = deepcopy(current_array)
            new_array[iteration_variable.array_index] = ratio
            setattr(
                iteration_variable.module,
                iteration_variable.target_name or iteration_variable.name,
                new_array,
            )

        check_iteration_variable(
            ratio, name=f"{variable_index} ({iteration_variable.name})"
        )


def load_scaled_bounds():
    """Sets the scaled bounds of the iteration variables."""

    for i in range(fortran.numerics.nvar.item()):
        variable_index = fortran.numerics.ixc[i].item() - 1
        fortran.numerics.itv_scaled_lower_bounds[i] = (
            fortran.numerics.boundl[variable_index] * fortran.numerics.scale[i]
        )
        fortran.numerics.itv_scaled_upper_bounds[i] = (
            fortran.numerics.boundu[variable_index] * fortran.numerics.scale[i]
        )


def initialise_iteration_variables():
    """Initialise the iteration variables (label and default bounds)"""
    for itv_index, itv in ITERATION_VARIABLES.items():
        fortran.numerics.lablxc[itv_index - 1] = string_to_f2py_compatible(
            fortran.numerics.lablxc[itv_index - 1], itv.name
        )

        fortran.numerics.boundl[itv_index - 1] = itv.lower_bound
        fortran.numerics.boundu[itv_index - 1] = itv.upper_bound
