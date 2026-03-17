from copy import deepcopy
from dataclasses import dataclass
from typing import Any
from warnings import warn

import numpy as np

from process import data_structure
from process.core.exceptions import ProcessValueError
from process.core.model import DataStructure


@dataclass
class IterationVariable:
    name: str
    """The name of the variable"""
    module: str | Any
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
    1: IterationVariable("aspect", "physics", 1.1, 10.00),
    2: IterationVariable("b_plasma_toroidal_on_axis", "physics", 0.010, 30.00),
    3: IterationVariable("rmajor", "physics", 0.1, 50.00),
    4: IterationVariable("temp_plasma_electron_vol_avg_kev", "physics", 5.0, 150.0),
    5: IterationVariable("beta_total_vol_avg", "physics", 0.001, 1.0),
    6: IterationVariable("nd_plasma_electrons_vol_avg", "physics", 2.0e19, 1.0e21),
    7: IterationVariable("f_nd_beam_electron", "physics", 1.0e-6, 1.0),
    10: IterationVariable("hfact", "physics", 0.1, 3.0),
    11: IterationVariable(
        "p_hcd_primary_extra_heat_mw",
        "current_drive",
        1.0e-3,
        1.0e3,
    ),
    12: IterationVariable("j_tf_coil_full_area", "tfcoil", 1.0e5, 1.50e8),
    13: IterationVariable("dr_tf_inboard", "build", 0.1, 5.0),
    16: IterationVariable("dr_cs", "build", 0.01, 10.00),
    17: IterationVariable("t_plant_pulse_dwell", "times", 0.1, 1.0e8),
    18: IterationVariable("q95", "physics", 2.0, 50.00),
    19: IterationVariable("e_beam_kev", "current_drive", 1.0, 1.0e6),
    20: IterationVariable("temp_cp_average", "tfcoil", 40.00, 573.0),
    23: IterationVariable("fcoolcp", "tfcoil", 0.1, 0.50),
    29: IterationVariable("dr_bore", "build", 0.1, 10.00),
    31: IterationVariable("gapomin", "build", 0.001, 1.0e1),
    37: IterationVariable("j_cs_flat_top_end", "pf_coil", 1.0e5, 1.0e8),
    41: IterationVariable("f_j_cs_start_pulse_end_flat_top", "pf_coil", 0.001, 1.0),
    42: IterationVariable("dr_cs_tf_gap", "build", 0.001, 10.00),
    44: IterationVariable("f_c_plasma_non_inductive", "physics", 0.001, 1.0),
    47: IterationVariable("feffcd", "current_drive", 0.001, 1.0),
    56: IterationVariable("t_tf_superconductor_quench", "tfcoil", 0.1, 100.0),
    57: IterationVariable("dr_tf_nose_case", "tfcoil", 0.05, 1.0),
    58: IterationVariable("dx_tf_turn_steel", "tfcoil", 0.001, 0.1),
    59: IterationVariable("f_a_tf_turn_cable_copper", "tfcoil", 0.001, 1.0),
    60: IterationVariable("c_tf_turn", "tfcoil", 0.001, 4.0e4),
    61: IterationVariable("dr_shld_vv_gap_inboard", "build", 0.001, 10.00),
    65: IterationVariable(
        "t_plant_pulse_plasma_current_ramp_up",
        "times",
        0.1,
        1.0e3,
    ),
    69: IterationVariable("radius_cp_coolant_channel", "tfcoil", 0.001, 0.010),
    70: IterationVariable("vel_cp_coolant_midplane", "tfcoil", 1.0, 1.0e2),
    73: IterationVariable("dr_fw_plasma_gap_inboard", "build", 0.001, 10.00),
    74: IterationVariable("dr_fw_plasma_gap_outboard", "build", 0.001, 10.00),
    75: IterationVariable("f_dr_tf_outboard_inboard", "build", 0.200, 5.0),
    81: IterationVariable("edrive", "ife", 1.0e5, 5.0e7),
    82: IterationVariable("drveff", "ife", 0.010, 1.0),
    83: IterationVariable("tgain", "ife", 1.0, 500.0),
    84: IterationVariable("chrad", "ife", 0.1, 20.00),
    85: IterationVariable("pdrive", "ife", 1.0e6, 200.0e6),
    90: IterationVariable("blbuith", "build", 0.001, 2.0),
    91: IterationVariable("blbuoth", "build", 0.001, 2.0),
    93: IterationVariable("dr_shld_inboard", "build", 0.001, 10.00),
    94: IterationVariable("dr_shld_outboard", "build", 0.001, 10.00),
    98: IterationVariable("f_blkt_li6_enrichment", "fwbs", 10.00, 100.0),
    104: IterationVariable("fcwr", "constraints", 0.001, 1.0),
    108: IterationVariable("breeder_f", "fwbs", 0.060, 1.0),
    109: IterationVariable("f_nd_alpha_electron", "physics", 0.05, 0.15),
    114: IterationVariable("len_fw_channel", "fwbs", 0.001, 1.0e3),
    119: IterationVariable("temp_plasma_separatrix_kev", "physics", 0.0, 1.0e1),
    122: IterationVariable("f_a_cs_turn_steel", "pf_coil", 0.001, 0.950),
    125: IterationVariable(
        "f_nd_impurity_electrons(03)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=2,
    ),
    126: IterationVariable(
        "f_nd_impurity_electrons(04)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=3,
    ),
    127: IterationVariable(
        "f_nd_impurity_electrons(05)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=4,
    ),
    128: IterationVariable(
        "f_nd_impurity_electrons(06)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=5,
    ),
    129: IterationVariable(
        "f_nd_impurity_electrons(07)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=6,
    ),
    130: IterationVariable(
        "f_nd_impurity_electrons(08)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=7,
    ),
    131: IterationVariable(
        "f_nd_impurity_electrons(09)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=8,
    ),
    132: IterationVariable(
        "f_nd_impurity_electrons(10)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=9,
    ),
    133: IterationVariable(
        "f_nd_impurity_electrons(11)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=10,
    ),
    134: IterationVariable(
        "f_nd_impurity_electrons(12)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=11,
    ),
    135: IterationVariable(
        "f_nd_impurity_electrons(13)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=12,
    ),
    136: IterationVariable(
        "f_nd_impurity_electrons(14)",
        "impurity_radiation",
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=13,
    ),
    138: IterationVariable(
        "dx_tf_hts_tape_rebco",
        "superconducting_tfcoil",
        0.01e-6,
        100.0e-6,
    ),
    139: IterationVariable(
        "dx_tf_hts_tape_copper",
        "superconducting_tfcoil",
        1.0e-6,
        1.0e-3,
    ),
    140: IterationVariable("dr_tf_wp_with_insulation", "tfcoil", 0.001, 2.0),
    142: IterationVariable(
        "nd_plasma_separatrix_electron",
        "physics",
        1.0e17,
        1.0e20,
    ),
    145: IterationVariable("f_nd_plasma_pedestal_greenwald", "physics", 0.1, 1.5),
    152: IterationVariable("f_nd_plasma_separatrix_greenwald", "physics", 0.001, 0.9),
    155: IterationVariable("pfusife", "ife", 5.0e2, 3.0e3),
    156: IterationVariable("rrin", "ife", 1.0, 1.0e1),
    158: IterationVariable(
        "dx_tf_croco_strand_copper",
        "superconducting_tfcoil",
        1.0e-3,
        1.0e-1,
    ),
    162: IterationVariable("r_cp_top", "build", 0.0010, 10.0),
    169: IterationVariable("te0_ecrh_achievable", "stellarator", 1.0, 1.0e3),
    170: IterationVariable("deg_div_field_plate", "divertor", 0.49, 5.01),
    171: IterationVariable("casths_fraction", "tfcoil", 0.01, 0.99),
    172: IterationVariable("dx_tf_side_case_min", "tfcoil", 0.001, 1.0),
    173: IterationVariable("f_plasma_fuel_tritium", "physics", 0.000, 1.000),
    174: IterationVariable("triang", "physics", 0.00, 1.00),
    175: IterationVariable("kappa", "physics", 0.00, 10.00),
    176: IterationVariable("f_st_coil_aspect", "stellarator", 0.70, 1.30),
    177: IterationVariable(
        "f_plasma_particles_lcfs_recycled", "physics", 0.01, 1.0
    ),
    178: IterationVariable(
        "eta_plasma_fuelling", "physics", 0.01, 1.0
    ),
    179: IterationVariable(
        "molflow_plasma_fuelling_vv_injected",
        "physics",
        1.0,
        1e24,
    ),
}


def check_iteration_variable(iteration_variable_value, name: str = ""):
    """Check that the iteration variable value is valid (not a weird number or too small).

    Raises an error upon encountering an invalid value, otherwise does nothing.

    Parameters
    ----------
    iteration_variable_value :

    name: str :
         (Default value = "")
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


def load_iteration_variables(data):
    """Loads the physics and engineering variables into the optimisation variable array."""
    for i in range(data_structure.numerics.nvar):
        variable_index = data_structure.numerics.ixc[i]
        iteration_variable = ITERATION_VARIABLES[variable_index]

        # use ... as the default return value because None might be a valid return from Fortran?

        module = (
            getattr(data, iteration_variable.module)
            if isinstance(iteration_variable.module, str)
            else iteration_variable.module
        )

        iteration_variable_value = getattr(
            module,
            iteration_variable.target_name or iteration_variable.name,
            ...,
        )

        # If iteration variable is missing
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

        data_structure.numerics.xcm[i] = iteration_variable_value
        data_structure.numerics.name_xc[i] = (
            data_structure.numerics.name_xc[i],
            iteration_variable.name,
        )

        # warn of the iteration variable is also a scan variable because this will cause
        # the optimiser and scan to overwrite the same variable and conflict
        if iteration_variable.name in {
            data.globals.vlabel,
            data.globals.vlabel_2,
        }:
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

        data_structure.numerics.scale[i] = 1.0 / iteration_variable_value
        data_structure.numerics.scafc[i] = 1.0 / data_structure.numerics.scale[i]

        data_structure.numerics.xcm[i] = (
            iteration_variable_value * data_structure.numerics.scale[i]
        )


def set_scaled_iteration_variable(xc, nn: int, data: DataStructure):
    """Converts scaled iteration variables back to their real values and sets them in the code.

    Parameters
    ----------
    xc :
        scaled iteration variable values
    nn :
        number of iteration variables
    data: DataStructure
        data structure
    """
    for i in range(nn):
        # there is less error handling here than in load_iteration_variables
        # because many errors will be caught in load_iteration_variables which is
        # run first. This verifies the variables exist and the module target is correct.
        variable_index = data_structure.numerics.ixc[i]
        iteration_variable = ITERATION_VARIABLES[variable_index]

        ratio = xc[i] / data_structure.numerics.scale[i]

        module = (
            getattr(data, iteration_variable.module)
            if isinstance(iteration_variable.module, str)
            else iteration_variable.module
        )

        if iteration_variable.array_index is None:
            setattr(
                module,
                iteration_variable.target_name or iteration_variable.name,
                ratio,
            )

        else:
            current_array = getattr(
                module,
                iteration_variable.target_name or iteration_variable.name,
            )
            new_array = deepcopy(current_array)
            new_array[iteration_variable.array_index] = ratio
            setattr(
                module,
                iteration_variable.target_name or iteration_variable.name,
                new_array,
            )

        check_iteration_variable(
            ratio, name=f"{variable_index} ({iteration_variable.name})"
        )


def load_scaled_bounds():
    """Sets the scaled bounds of the iteration variables."""
    for i in range(data_structure.numerics.nvar):
        variable_index = data_structure.numerics.ixc[i] - 1
        data_structure.numerics.itv_scaled_lower_bounds[i] = (
            data_structure.numerics.boundl[variable_index]
            * data_structure.numerics.scale[i]
        )
        data_structure.numerics.itv_scaled_upper_bounds[i] = (
            data_structure.numerics.boundu[variable_index]
            * data_structure.numerics.scale[i]
        )


def initialise_iteration_variables():
    """Initialise the iteration variables (label and default bounds)"""
    for itv_index, itv in ITERATION_VARIABLES.items():
        data_structure.numerics.lablxc[itv_index - 1] = itv.name

        data_structure.numerics.boundl[itv_index - 1] = itv.lower_bound
        data_structure.numerics.boundu[itv_index - 1] = itv.upper_bound
