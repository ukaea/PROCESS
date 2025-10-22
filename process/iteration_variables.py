from copy import deepcopy
from dataclasses import dataclass
from typing import Any
from warnings import warn

import numpy as np

import process.data_structure as data_structure
from process.exceptions import ProcessValueError


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
    1: IterationVariable("aspect", data_structure.physics_variables, 1.1, 10.00),
    2: IterationVariable(
        "b_plasma_toroidal_on_axis", data_structure.physics_variables, 0.010, 30.00
    ),
    3: IterationVariable("rmajor", data_structure.physics_variables, 0.1, 50.00),
    4: IterationVariable(
        "temp_plasma_electron_vol_avg_kev", data_structure.physics_variables, 5.0, 150.0
    ),
    5: IterationVariable(
        "beta_total_vol_avg", data_structure.physics_variables, 0.001, 1.0
    ),
    6: IterationVariable(
        "nd_plasma_electrons_vol_avg", data_structure.physics_variables, 2.0e19, 1.0e21
    ),
    7: IterationVariable(
        "f_nd_beam_electron", data_structure.physics_variables, 1.0e-6, 1.0
    ),
    10: IterationVariable("hfact", data_structure.physics_variables, 0.1, 3.0),
    11: IterationVariable(
        "p_hcd_primary_extra_heat_mw",
        data_structure.current_drive_variables,
        1.0e-3,
        1.0e3,
    ),
    12: IterationVariable("oacdcp", data_structure.tfcoil_variables, 1.0e5, 1.50e8),
    13: IterationVariable("dr_tf_inboard", data_structure.build_variables, 0.1, 5.0),
    16: IterationVariable("dr_cs", data_structure.build_variables, 0.01, 10.00),
    17: IterationVariable(
        "t_plant_pulse_dwell", data_structure.times_variables, 0.1, 1.0e8
    ),
    18: IterationVariable("q95", data_structure.physics_variables, 2.0, 50.00),
    19: IterationVariable(
        "e_beam_kev", data_structure.current_drive_variables, 1.0, 1.0e6
    ),
    20: IterationVariable(
        "temp_cp_average", data_structure.tfcoil_variables, 40.00, 573.0
    ),
    23: IterationVariable("fcoolcp", data_structure.tfcoil_variables, 0.1, 0.50),
    29: IterationVariable("dr_bore", data_structure.build_variables, 0.1, 10.00),
    31: IterationVariable("gapomin", data_structure.build_variables, 0.001, 1.0e1),
    37: IterationVariable(
        "j_cs_flat_top_end", data_structure.pfcoil_variables, 1.0e5, 1.0e8
    ),
    41: IterationVariable(
        "f_j_cs_start_pulse_end_flat_top", data_structure.pfcoil_variables, 0.001, 1.0
    ),
    42: IterationVariable("dr_cs_tf_gap", data_structure.build_variables, 0.001, 10.00),
    44: IterationVariable(
        "f_c_plasma_non_inductive", data_structure.physics_variables, 0.001, 1.0
    ),
    47: IterationVariable("feffcd", data_structure.current_drive_variables, 0.001, 1.0),
    56: IterationVariable(
        "t_tf_superconductor_quench", data_structure.tfcoil_variables, 0.1, 100.0
    ),
    57: IterationVariable("dr_tf_nose_case", data_structure.tfcoil_variables, 0.05, 1.0),
    58: IterationVariable(
        "dx_tf_turn_steel", data_structure.tfcoil_variables, 0.001, 0.1
    ),
    59: IterationVariable(
        "f_a_tf_turn_cable_copper", data_structure.tfcoil_variables, 0.001, 1.0
    ),
    60: IterationVariable("c_tf_turn", data_structure.tfcoil_variables, 0.001, 4.0e4),
    61: IterationVariable(
        "dr_shld_vv_gap_inboard", data_structure.build_variables, 0.001, 10.00
    ),
    65: IterationVariable(
        "t_plant_pulse_plasma_current_ramp_up",
        data_structure.times_variables,
        0.1,
        1.0e3,
    ),
    69: IterationVariable(
        "radius_cp_coolant_channel", data_structure.tfcoil_variables, 0.001, 0.010
    ),
    70: IterationVariable(
        "vel_cp_coolant_midplane", data_structure.tfcoil_variables, 1.0, 1.0e2
    ),
    73: IterationVariable(
        "dr_fw_plasma_gap_inboard", data_structure.build_variables, 0.001, 10.00
    ),
    74: IterationVariable(
        "dr_fw_plasma_gap_outboard", data_structure.build_variables, 0.001, 10.00
    ),
    75: IterationVariable(
        "f_dr_tf_outboard_inboard", data_structure.build_variables, 0.200, 5.0
    ),
    81: IterationVariable("edrive", data_structure.ife_variables, 1.0e5, 5.0e7),
    82: IterationVariable("drveff", data_structure.ife_variables, 0.010, 1.0),
    83: IterationVariable("tgain", data_structure.ife_variables, 1.0, 500.0),
    84: IterationVariable("chrad", data_structure.ife_variables, 0.1, 20.00),
    85: IterationVariable("pdrive", data_structure.ife_variables, 1.0e6, 200.0e6),
    90: IterationVariable("blbuith", data_structure.build_variables, 0.001, 2.0),
    91: IterationVariable("blbuoth", data_structure.build_variables, 0.001, 2.0),
    93: IterationVariable(
        "dr_shld_inboard", data_structure.build_variables, 0.001, 10.00
    ),
    94: IterationVariable(
        "dr_shld_outboard", data_structure.build_variables, 0.001, 10.00
    ),
    98: IterationVariable(
        "f_blkt_li6_enrichment", data_structure.fwbs_variables, 10.00, 100.0
    ),
    104: IterationVariable("fcwr", data_structure.constraint_variables, 0.001, 1.0),
    108: IterationVariable("breeder_f", data_structure.fwbs_variables, 0.060, 1.0),
    109: IterationVariable(
        "f_nd_alpha_electron", data_structure.physics_variables, 0.05, 0.15
    ),
    114: IterationVariable(
        "len_fw_channel", data_structure.fwbs_variables, 0.001, 1.0e3
    ),
    119: IterationVariable(
        "temp_plasma_separatrix_kev", data_structure.physics_variables, 0.0, 1.0e1
    ),
    122: IterationVariable(
        "f_a_cs_turn_steel", data_structure.pfcoil_variables, 0.001, 0.950
    ),
    125: IterationVariable(
        "f_nd_impurity_electrons(03)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=2,
    ),
    126: IterationVariable(
        "f_nd_impurity_electrons(04)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=3,
    ),
    127: IterationVariable(
        "f_nd_impurity_electrons(05)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=4,
    ),
    128: IterationVariable(
        "f_nd_impurity_electrons(06)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=5,
    ),
    129: IterationVariable(
        "f_nd_impurity_electrons(07)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=6,
    ),
    130: IterationVariable(
        "f_nd_impurity_electrons(08)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=7,
    ),
    131: IterationVariable(
        "f_nd_impurity_electrons(09)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=8,
    ),
    132: IterationVariable(
        "f_nd_impurity_electrons(10)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=9,
    ),
    133: IterationVariable(
        "f_nd_impurity_electrons(11)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=10,
    ),
    134: IterationVariable(
        "f_nd_impurity_electrons(12)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=11,
    ),
    135: IterationVariable(
        "f_nd_impurity_electrons(13)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=12,
    ),
    136: IterationVariable(
        "f_nd_impurity_electrons(14)",
        data_structure.impurity_radiation_module,
        1e-8,
        0.01,
        target_name="f_nd_impurity_electron_array",
        array_index=13,
    ),
    138: IterationVariable(
        "dx_hts_tape_rebco", data_structure.physics_variables, 0.01e-6, 100.0e-6
    ),
    139: IterationVariable(
        "dx_hts_tape_copper", data_structure.rebco_variables, 1.0e-6, 1.0e-3
    ),
    140: IterationVariable(
        "dr_tf_wp_with_insulation", data_structure.tfcoil_variables, 0.001, 2.0
    ),
    142: IterationVariable(
        "nd_plasma_separatrix_electron",
        data_structure.physics_variables,
        1.0e17,
        1.0e20,
    ),
    145: IterationVariable(
        "f_nd_plasma_pedestal_greenwald", data_structure.physics_variables, 0.1, 0.9
    ),
    152: IterationVariable(
        "f_nd_plasma_separatrix_greenwald", data_structure.physics_variables, 0.001, 0.5
    ),
    155: IterationVariable("pfusife", data_structure.ife_variables, 5.0e2, 3.0e3),
    156: IterationVariable("rrin", data_structure.ife_variables, 1.0, 1.0e1),
    158: IterationVariable(
        "dx_croco_strand_copper", data_structure.rebco_variables, 1.0e-3, 1.0e-1
    ),
    162: IterationVariable("r_cp_top", data_structure.build_variables, 0.0010, 10.0),
    169: IterationVariable(
        "te0_ecrh_achievable", data_structure.stellarator_variables, 1.0, 1.0e3
    ),
    170: IterationVariable(
        "deg_div_field_plate", data_structure.divertor_variables, 0.49, 5.01
    ),
    171: IterationVariable(
        "casths_fraction", data_structure.tfcoil_variables, 0.01, 0.99
    ),
    172: IterationVariable(
        "dx_tf_side_case_min", data_structure.tfcoil_variables, 0.001, 1.0
    ),
    173: IterationVariable(
        "f_plasma_fuel_tritium", data_structure.physics_variables, 0.000, 1.000
    ),
    174: IterationVariable("triang", data_structure.physics_variables, 0.00, 1.00),
    175: IterationVariable("kappa", data_structure.physics_variables, 0.00, 10.00),
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

    for i in range(data_structure.numerics.nvar):
        variable_index = data_structure.numerics.ixc[i]
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

        data_structure.numerics.xcm[i] = iteration_variable_value
        data_structure.numerics.name_xc[i] = (
            data_structure.numerics.name_xc[i],
            iteration_variable.name,
        )

        # warn of the iteration variable is also a scan variable because this will cause
        # the optimiser and scan to overwrite the same variable and conflict
        if iteration_variable.name in (
            data_structure.global_variables.vlabel,
            data_structure.global_variables.vlabel_2,
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

        data_structure.numerics.scale[i] = 1.0 / iteration_variable_value
        data_structure.numerics.scafc[i] = 1.0 / data_structure.numerics.scale[i]

        data_structure.numerics.xcm[i] = (
            iteration_variable_value * data_structure.numerics.scale[i]
        )


def set_scaled_iteration_variable(xc, nn: int):
    """Converts scaled iteration variables back to their real values and sets them in the code.

    :param xc: scaled iteration variable values
    :param nn: number of iteration variables
    """

    for i in range(nn):
        # there is less error handling here than in load_iteration_variables
        # because many errors will be caught in load_iteration_variables which is
        # run first. This verifies the variables exist and the module target is correct.
        variable_index = data_structure.numerics.ixc[i]
        iteration_variable = ITERATION_VARIABLES[variable_index]

        ratio = xc[i] / data_structure.numerics.scale[i]

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
