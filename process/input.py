"""Handle parsing, validation, and actioning of a PROCESS input file (*IN.DAT)."""

import copy
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from warnings import warn

import process
import process.data_structure as data_structure
from process.constraints import ConstraintManager
from process.exceptions import ProcessValidationError, ProcessValueError

NumberType = int | float
ValidInputTypes = NumberType | str

DataTypes = (int, float, str)
"""Valid input variable types alias"""


def _ixc_additional_actions(_name, value: int, _array_index, _config):
    data_structure.numerics.ixc[data_structure.numerics.nvar] = value
    data_structure.numerics.nvar += 1


def _icc_additional_actions(_name, value: int, _array_index, _config):
    data_structure.numerics.icc[data_structure.numerics.n_constraints] = value
    data_structure.numerics.n_constraints += 1


@dataclass
class InputVariable:
    """A variable to be parsed from the input file."""

    module: Any
    """The Fortran module that this variable should be set on."""
    type: type
    """The expected type of the variable"""
    range: tuple[NumberType, NumberType] | None = None
    """The valid range of values (int and float only) inclusive of the endpoints."""
    choices: list[ValidInputTypes] | None = None
    """The valid values of this input."""
    array: bool = False
    """Is this input assigning values to an array?"""
    additional_validation: (
        Callable[[str, ValidInputTypes, int | None, "InputVariable"], ValidInputTypes]
        | None
    ) = None
    """A function that takes the input variable: name, value, array index, and config (this dataclass)
    as input and returns a cleaned version of the input value. May raise
    a ProcessValidationError.

    NOTE: The value passed in has already been cleaned in the default way and has
    been cast to the specified `type`.
    """
    additional_actions: (
        Callable[[str, ValidInputTypes, int | None, "InputVariable"], None] | None
    ) = None
    """A function that takes the input variable: name, value, array index, and config (this dataclass)
    as input and performs some additional action in addition to the default actions prescribed by the variables
    config. May raise a ProcessValidationError.

    NOTE: The value passed in has already been cleaned in the default
    way and has been cast to the specified `type`.

    NOTE: the input name is the cleaned name as in the input file NOT the `target_variable`.
    """
    target_name: str | None = None
    """Indicates the parsed variable name is different to its target on the `module`."""
    set_variable: bool = True
    """Do not attempt to set the variable on the module. Instead only do any `additional_actions`."""

    def __post_init__(self):
        if self.type not in DataTypes:
            error_msg = (
                f"Cannot parse variable of type {self.type}, must be one of {DataTypes}"
            )
            raise ProcessValueError(error_msg)

        if self.type not in (int, float) and self.range is not None:
            error_msg = "Can only apply a range to integer or float variables"
            raise ProcessValueError(error_msg)

        if not self.set_variable and self.additional_actions is None:
            warn(
                "Not setting variable or performing an additional action. Why are you parsing this variable?",
                stacklevel=2,
            )


INPUT_VARIABLES = {
    "runtitle": InputVariable(data_structure.global_variables, str),
    "verbose": InputVariable(data_structure.global_variables, int, choices=[0, 1]),
    "run_tests": InputVariable(data_structure.global_variables, int, choices=[0, 1]),
    "ioptimz": InputVariable(data_structure.numerics, int, choices=[1, -2]),
    "epsvmc": InputVariable(data_structure.numerics, float, range=(0.0, 1.0)),
    "boundl": InputVariable(data_structure.numerics, float, array=True),
    "boundu": InputVariable(data_structure.numerics, float, array=True),
    "epsfcn": InputVariable(data_structure.numerics, float, range=(0.0, 1.0)),
    "maxcal": InputVariable(data_structure.global_variables, int, range=(0, 10000)),
    "minmax": InputVariable(data_structure.numerics, int),
    "neqns": InputVariable(
        data_structure.numerics, int, range=(0, ConstraintManager.num_constraints())
    ),
    "nineqns": InputVariable(
        data_structure.numerics, int, range=(0, ConstraintManager.num_constraints())
    ),
    "alphaj": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "alphan": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "alphat": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "aspect": InputVariable(
        data_structure.physics_variables, float, range=(1.001, 40.0)
    ),
    "beamfus0": InputVariable(
        data_structure.physics_variables, float, range=(0.01, 10.0)
    ),
    "beta_total_vol_avg": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "beta_vol_avg_max": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "beta_vol_avg_min": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "betbm0": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "b_plasma_toroidal_on_axis": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 30.0)
    ),
    "burnup_in": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "radius_plasma_core_norm": InputVariable(
        data_structure.impurity_radiation_module, float, range=(0.0, 1.0)
    ),
    "f_p_plasma_core_rad_reduction": InputVariable(
        data_structure.impurity_radiation_module, float, range=(0.0, 1.0)
    ),
    "c_beta": InputVariable(data_structure.physics_variables, float, range=(0.0, 1.0)),
    "csawth": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "f_vol_plasma": InputVariable(
        data_structure.physics_variables, float, range=(0.001, 10.0)
    ),
    "f_r_conducting_wall": InputVariable(
        data_structure.physics_variables, float, range=(1.0, 3.0)
    ),
    "nd_plasma_electrons_vol_avg": InputVariable(
        data_structure.physics_variables, float, range=(1.0e18, 1.0e22)
    ),
    "beta_norm_max": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 20.0)
    ),
    "beta_poloidal_eps_max": InputVariable(
        data_structure.physics_variables, float, range=(0.01, 10.0)
    ),
    "f_p_alpha_plasma_deposited": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_p_div_lower": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_plasma_fuel_deuterium": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "ffwal": InputVariable(data_structure.physics_variables, float, range=(0.0, 10.0)),
    "f_nd_plasma_pedestal_greenwald": InputVariable(
        data_structure.physics_variables, float, range=(-1.0, 5.0)
    ),
    "f_nd_plasma_separatrix_greenwald": InputVariable(
        data_structure.physics_variables, float, range=(-1.0, 5.0)
    ),
    "f_plasma_fuel_helium3": InputVariable(
        data_structure.physics_variables, float, range=(-1.0, 5.0)
    ),
    # TODO: does f_nd_impurity_electrons require an additional range?
    "f_nd_impurity_electrons": InputVariable(
        data_structure.impurity_radiation_module, float, array=True
    ),
    "fkzohm": InputVariable(data_structure.physics_variables, float, range=(0.5, 2.0)),
    "abktflnc": InputVariable(data_structure.cost_variables, float, range=(0.1, 100.0)),
    "adivflnc": InputVariable(data_structure.cost_variables, float, range=(0.1, 100.0)),
    "admv": InputVariable(
        data_structure.buildings_variables, float, range=(1.0e4, 1.0e6)
    ),
    "airtemp": InputVariable(
        data_structure.water_usage_variables, float, range=(-15.0, 40.0)
    ),
    "alfapf": InputVariable(data_structure.pfcoil_variables, float, range=(1e-12, 1.0)),
    "alstroh": InputVariable(
        data_structure.pfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "amortization": InputVariable(
        data_structure.cost_variables, float, range=(1.0, 50.0)
    ),
    "anginc": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "aplasmin": InputVariable(data_structure.build_variables, float, range=(0.01, 10.0)),
    "aux_build_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "aux_build_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "aux_build_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "auxcool_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "auxcool_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "auxcool_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "p_hcd_injected_min_mw": InputVariable(
        data_structure.constraint_variables, float, range=(0.01, 100.0)
    ),
    "avail_min": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "b_crit_upper_nbti": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 30.0)
    ),
    "p_plant_electric_base": InputVariable(
        data_structure.heat_transport_variables, float, range=(1000000.0, 10000000000.0)
    ),
    "bcritsc": InputVariable(data_structure.tfcoil_variables, float, range=(10.0, 50.0)),
    "bctmp": InputVariable(data_structure.pulse_variables, float, range=(1.0, 800.0)),
    "e_beam_kev": InputVariable(
        data_structure.current_drive_variables, float, range=(1.0, 1000000.0)
    ),
    "dx_beam_duct": InputVariable(
        data_structure.current_drive_variables, float, range=(0.001, 5.0)
    ),
    "deg_div_field_plate": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 360.0)
    ),
    "beta_poloidal_max": InputVariable(
        data_structure.constraint_variables, float, range=(0.01, 2.0)
    ),
    "betai": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "betao": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "big_q_plasma_min": InputVariable(
        data_structure.constraint_variables, float, range=(0.01, 100.0)
    ),
    "bioshld_thk": InputVariable(
        data_structure.buildings_variables, float, range=(0.25, 25.0)
    ),
    "bkt_life_csf": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.0, 1.0)
    ),
    "blbmith": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "blbmoth": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "blbpith": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "blbpoth": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "blbuith": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "blbuoth": InputVariable(data_structure.build_variables, float, range=(0.0, 2.0)),
    "bldr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "bldrc": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "bldzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "bldzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "pres_blkt_coolant": InputVariable(
        data_structure.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "blpressure_liq": InputVariable(
        data_structure.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "b_cs_limit_max": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.01, 100.0)
    ),
    "bmn": InputVariable(
        data_structure.stellarator_variables, float, range=(0.0001, 0.01)
    ),
    "b_tf_inboard_max": InputVariable(
        data_structure.constraint_variables, float, range=(0.1, 50.0)
    ),
    "f_c_plasma_bootstrap_max": InputVariable(
        data_structure.current_drive_variables, float, range=(-0.999, 0.999)
    ),
    "f_c_plasma_bootstrap": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "breeder_f": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "breeder_multiplier": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "bz_channel_conduct_liq": InputVariable(
        data_structure.fwbs_variables, float, range=(1e-06, 1000000.0)
    ),
    "dr_tf_plasma_case": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_dr_tf_plasma_case": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dx_tf_side_case_min": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "casths_fraction": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "cboot": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 10.0)
    ),
    "cconfix": InputVariable(data_structure.cost_variables, float, range=(50.0, 200.0)),
    "cconshpf": InputVariable(data_structure.cost_variables, float, range=(50.0, 200.0)),
    "cconshtf": InputVariable(data_structure.cost_variables, float, range=(50.0, 200.0)),
    "cdriv0": InputVariable(data_structure.ife_variables, float, range=(50.0, 500.0)),
    "cdriv1": InputVariable(data_structure.ife_variables, float, range=(50.0, 500.0)),
    "cdriv2": InputVariable(data_structure.ife_variables, float, range=(50.0, 500.0)),
    "f_t_plant_available": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "chdzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "chdzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "chemlab_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "chemlab_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "chemlab_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "chrad": InputVariable(data_structure.ife_variables, float, range=(0.1, 20.0)),
    "cland": InputVariable(data_structure.cost_variables, float, range=(10.0, 100.0)),
    "clh2": InputVariable(data_structure.buildings_variables, float, range=(0.0, 30.0)),
    "j_cs_flat_top_end": InputVariable(
        data_structure.pfcoil_variables, float, range=(10000.0, 500000000.0)
    ),
    "conf_mag": InputVariable(data_structure.cost_variables, float, range=(0.9, 1.0)),
    "control_buildings_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "control_buildings_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "control_buildings_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "conv": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "coolp": InputVariable(
        data_structure.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "copper_rrr": InputVariable(
        data_structure.rebco_variables, float, range=(1.0, 10000.0)
    ),
    "dx_hts_tape_copper": InputVariable(
        data_structure.rebco_variables, float, range=(0.0, 0.001)
    ),
    "copperaoh_m2": InputVariable(
        data_structure.rebco_variables, float, range=(1.0, 10000000000.0)
    ),
    "copperaoh_m2_max": InputVariable(
        data_structure.rebco_variables, float, range=(10000.0, 10000000000.0)
    ),
    "cost_factor_bop": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_buildings": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_fwbs": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_land": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_misc": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_rh": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_tf_coils": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "cost_factor_vv": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 10.0)
    ),
    "costexp": InputVariable(data_structure.cost_variables, float, range=(0.01, 5.0)),
    "costexp_pebbles": InputVariable(
        data_structure.cost_variables, float, range=(0.01, 5.0)
    ),
    "cowner": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "cplife_input": InputVariable(
        data_structure.cost_variables, float, range=(0.001, 50.0)
    ),
    "cpstflnc": InputVariable(data_structure.cost_variables, float, range=(0.01, 30.0)),
    "c_tf_turn": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.001, 1000000.0)
    ),
    "c_tf_turn_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(1.0, 1000000.0)
    ),
    "crane_arm_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "crane_clrnc_h": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "crane_clrnc_v": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "dx_croco_strand_copper": InputVariable(
        data_structure.rebco_variables, float, range=(0.001, 0.1)
    ),
    "cryomag_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "cryomag_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryomag_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryostat_clrnc": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "cryostore_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "cryostore_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryostore_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "p_cryo_plant_electric_max_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.01, 200.0)
    ),
    "csi": InputVariable(data_structure.cost_variables, float, range=(1.0, 100.0)),
    "cturbb": InputVariable(data_structure.cost_variables, float, range=(100.0, 1000.0)),
    "dz_vv_lower": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dz_vv_upper": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "den_aluminium": InputVariable(process.constants, float, range=(2500.0, 30000.0)),
    "den_tf_coil_case": InputVariable(
        data_structure.tfcoil_variables, float, range=(1000.0, 100000.0)
    ),
    "dcdrv0": InputVariable(data_structure.ife_variables, float, range=(0.0, 200.0)),
    "dcdrv1": InputVariable(data_structure.ife_variables, float, range=(0.0, 200.0)),
    "dcdrv2": InputVariable(data_structure.ife_variables, float, range=(0.0, 200.0)),
    "den_tf_wp_turn_insulation": InputVariable(
        data_structure.tfcoil_variables, float, range=(500.0, 10000.0)
    ),
    "den_copper": InputVariable(process.constants, float, range=(8000.0, 10000.0)),
    "declblkt": InputVariable(data_structure.fwbs_variables, float, range=(0.01, 0.2)),
    "declfw": InputVariable(data_structure.fwbs_variables, float, range=(0.01, 0.2)),
    "declshld": InputVariable(data_structure.fwbs_variables, float, range=(0.01, 0.2)),
    "decomf": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "den_steel": InputVariable(
        data_structure.fwbs_variables, float, range=(5000.0, 10000.0)
    ),
    "dia_tf_turn_coolant_channel": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "dia_tf_turn_superconducting_cable": InputVariable(
        data_structure.superconducting_tf_coil_variables, float, range=(0.0001, 0.01)
    ),
    "dintrt": InputVariable(data_structure.cost_variables, float, range=(0.0, 0.1)),
    "discount_rate": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 0.5)
    ),
    "div_nref": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000000.0)
    ),
    "f_vol_div_coolant": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.0)
    ),
    "den_div_structure": InputVariable(
        data_structure.divertor_variables, float, range=(0.1, 100000.0)
    ),
    "dz_divertor": InputVariable(
        data_structure.divertor_variables, float, range=(0.1, 5.0)
    ),
    "dx_div_plate": InputVariable(
        data_structure.divertor_variables, float, range=(0.01, 1.0)
    ),
    "dp_blkt": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10000000.0)
    ),
    "dp_fw": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10000000.0)
    ),
    "dp_fw_blkt": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10000000.0)
    ),
    "dp_he": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10000000.0)
    ),
    "dp_liq": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10000000.0)
    ),
    "f_p_fw_blkt_pump": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 10.0)
    ),
    "dr_blkt_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_blkt_outboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_bore": InputVariable(data_structure.build_variables, float, range=(0.0, 50.0)),
    "dr_cryostat": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_cs": InputVariable(data_structure.build_variables, float, range=(0.0, 10.0)),
    "dr_cs_tf_gap": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_fw_plasma_gap_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_fw_plasma_gap_outboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_fw_wall": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0005, 0.1)
    ),
    "dr_pf_cryostat": InputVariable(
        data_structure.fwbs_variables, float, range=(0.1, 5.0)
    ),
    "dr_shld_blkt_gap": InputVariable(
        data_structure.build_variables, float, range=(0.0, 5.0)
    ),
    "dr_shld_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_outboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_thermal_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_thermal_outboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_vv_gap_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_tf_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_tf_shld_gap": InputVariable(
        data_structure.build_variables, float, range=(0.0, 5.0)
    ),
    "dr_tf_wp_with_insulation": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 10.0)
    ),
    "dr_vv_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_vv_outboard": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "drtop": InputVariable(data_structure.tfcoil_variables, float, range=(-1.5, 1.5)),
    "drveff": InputVariable(data_structure.ife_variables, float, range=(0.01, 1.0)),
    "dtlife": InputVariable(data_structure.cost_variables, float, range=(0.0, 15.0)),
    "dtstor": InputVariable(data_structure.pulse_variables, float, range=(50.0, 500.0)),
    "dx_fw_module": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0005, 0.1)
    ),
    "dz_tf_cryostat": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 20.0)
    ),
    "dztop": InputVariable(data_structure.tfcoil_variables, float, range=(-0.5, 0.5)),
    "edrive": InputVariable(
        data_structure.ife_variables, float, range=(100000.0, 5000000000.0)
    ),
    "eff_tf_cryo": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "ejima_coeff": InputVariable(
        data_structure.physics_variables, float, range=(0.1, 1.0)
    ),
    "elecdist_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "elecdist_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecdist_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecload_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "elecload_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecload_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecstore_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "elecstore_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecstore_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "f_p_blkt_multiplication": InputVariable(
        data_structure.fwbs_variables, float, range=(1.0, 2.0)
    ),
    "esbldgm3": InputVariable(
        data_structure.buildings_variables, float, range=(1000.0, 1000000.0)
    ),
    "eta_ecrh_injector_wall_plug": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "eta_icrh_injector_wall_plug": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "eta_ebw_injector_wall_plug": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "eta_coolant_pump_electric": InputVariable(
        data_structure.fwbs_variables, float, range=(0.1, 1.0)
    ),
    "etaiso": InputVariable(data_structure.fwbs_variables, float, range=(0.1, 1.0)),
    "eta_lowhyb_injector_wall_plug": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "etali": InputVariable(data_structure.ife_variables, float, range=(0.0, 1.0)),
    "eta_beam_injector_wall_plug": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "etapsu": InputVariable(data_structure.pfcoil_variables, float, range=(0.0, 1.0)),
    "etapump": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 1.0)),
    "etatf": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1.0)
    ),
    "eta_turbine": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1.0)
    ),
    "eyoung_al": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 1.0)),
    "eyoung_cond_axial": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 10000000000000.0)
    ),
    "eyoung_cond_trans": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 10000000000000.0)
    ),
    "eyoung_ins": InputVariable(
        data_structure.tfcoil_variables, float, range=(100000000.0, 10000000000000.0)
    ),
    "eyoung_res_tf_buck": InputVariable(
        data_structure.tfcoil_variables, float, range=(1e-10, 1000000000000.0)
    ),
    "eyoung_steel": InputVariable(
        data_structure.tfcoil_variables, float, range=(100000000.0, 10000000000000.0)
    ),
    "f_a_tf_cool_outboard": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_alpha_energy_confinement_min": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 100.0)
    ),
    "f_asym": InputVariable(
        data_structure.stellarator_variables, float, range=(0.9, 2.0)
    ),
    "f_fw_peak": InputVariable(data_structure.fwbs_variables, float, range=(1.0, 100.0)),
    "f_fw_rad_max": InputVariable(
        data_structure.constraint_variables, float, range=(0.1, 10)
    ),
    "f_nd_alpha_electron": InputVariable(
        data_structure.physics_variables, float, range=(1e-12, 1.0)
    ),
    "f_nd_beam_electron": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_nd_protium_electrons": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_neut_shield": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "f_nuc_pow_bz_struct": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "f_r_cp": InputVariable(data_structure.build_variables, float, range=(1.0, 100.0)),
    "f_rad": InputVariable(
        data_structure.stellarator_variables, float, range=(0.0, 1.0)
    ),
    "f_sync_reflect": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_plasma_fuel_tritium": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_beam_tritium": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "f_vforce_inboard": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_w": InputVariable(data_structure.stellarator_variables, float, range=(0.1, 1.0)),
    "f_z_cryostat": InputVariable(
        data_structure.build_variables, float, range=(2.0, 10.0)
    ),
    "fauxbop": InputVariable(data_structure.ife_variables, float, range=(0.0, 1.0)),
    "fblbe": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblbreed": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebmi": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebmo": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebpi": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebpo": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblli": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblli2o": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fbllipb": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblss": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fblvd": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fbreed": InputVariable(data_structure.ife_variables, float, range=(0.0, 0.999)),
    "fburn": InputVariable(data_structure.ife_variables, float, range=(0.01, 1.0)),
    "fc_building_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "fc_building_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "fcap0": InputVariable(data_structure.cost_variables, float, range=(1.0, 1.5)),
    "fcap0cp": InputVariable(data_structure.cost_variables, float, range=(1.0, 1.5)),
    "fcdfuel": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "f_j_cs_start_pulse_end_flat_top": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 1.0)
    ),
    "fcontng": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "fcoolcp": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 1.0)),
    "fcr0": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "fcspc": InputVariable(data_structure.build_variables, float, range=(0.0, 1.0)),
    "fcuohsu": InputVariable(data_structure.pfcoil_variables, float, range=(0.0, 1.0)),
    "fcupfsu": InputVariable(data_structure.pfcoil_variables, float, range=(0.0, 1.0)),
    "f_a_tf_turn_cable_copper": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "fdene": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fiooic": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1.0)
    ),
    "fjohc": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1.0)
    ),
    "fjohc0": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1.0)
    ),
    "f_ster_div_single": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "fdiva": InputVariable(data_structure.divertor_variables, float, range=(0.1, 2.0)),
    "fdivwet": InputVariable(
        data_structure.stellarator_variables, float, range=(0.01, 1.0)
    ),
    "feffcd": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 20.0)
    ),
    "f_a_fw_outboard_hcd": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "fhole": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fhts": InputVariable(data_structure.tfcoil_variables, float, range=(0.01, 1.0)),
    "fkind": InputVariable(data_structure.cost_variables, float, range=(0.5, 1.0)),
    "h_mode_threshold_margin": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1000000.0)
    ),
    "l_mode_threshold_margin": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1000000.0)
    ),
    "flirad": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "flpitch": InputVariable(
        data_structure.stellarator_variables, float, range=(0.0001, 0.01)
    ),
    "f_div_flux_expansion": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "fmgdmw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "fndt": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "f_p_beam_orbit_loss": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 0.999)
    ),
    "f_p_blkt_coolant_pump_total_heat": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "f_p_div_coolant_pump_total_heat": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "f_p_fw_coolant_pump_total_heat": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "f_p_shld_coolant_pump_total_heat": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "fracture_toughness": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.1, 100000000.0)
    ),
    "fradpwr": InputVariable(
        data_structure.constraint_variables, float, range=(0.0, 1.0)
    ),
    "f_radius_beam_tangency_rmajor": InputVariable(
        data_structure.current_drive_variables, float, range=(0.5, 2.0)
    ),
    "frhocp": InputVariable(data_structure.tfcoil_variables, float, range=(0.01, 5.0)),
    "frholeg": InputVariable(data_structure.tfcoil_variables, float, range=(0.01, 5.0)),
    "fseppc": InputVariable(
        data_structure.build_variables, float, range=(1000000.0, 1000000000.0)
    ),
    "fvoldw": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 10.0)),
    "fvolsi": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 10.0)),
    "fvolso": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 10.0)),
    "f_c_plasma_non_inductive": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "fw_armour_thickness": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "fw_th_conductivity": InputVariable(
        data_structure.fwbs_variables, float, range=(1.0, 100.0)
    ),
    "fwbs_nref": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000000.0)
    ),
    "fwbs_nu": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000000.0)
    ),
    "fwbs_prob_fail": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "fwbs_umain_time": InputVariable(
        data_structure.cost_variables, float, range=(0.1, 2.0)
    ),
    "fwclfr": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "fwdr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "fwdzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "fwdzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "fzactual": InputVariable(data_structure.reinke_variables, float, range=(0.0, 1.0)),
    "eta_cd_norm_ecrh": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "gamma_he": InputVariable(
        data_structure.primary_pumping_variables, float, range=(1.0, 2.0)
    ),
    "eta_cd_norm_hcd_primary_max": InputVariable(
        data_structure.constraint_variables, float, range=(0.01, 10.0)
    ),
    "gapomin": InputVariable(data_structure.build_variables, float, range=(0.0, 10.0)),
    "gas_buildings_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "gas_buildings_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "gas_buildings_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ground_clrnc": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "n_ecrh_harmonic": InputVariable(
        data_structure.current_drive_variables, float, range=(1.0, 10.0)
    ),
    "dx_hts_tape_hastelloy": InputVariable(
        data_structure.rebco_variables, float, range=(1e-08, 0.001)
    ),
    "hccl": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "hcd_building_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "hcd_building_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hcd_building_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hcwt": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "heat_sink_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "heat_sink_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "heat_sink_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hfact": InputVariable(data_structure.physics_variables, float, range=(0.01, 10.0)),
    "pflux_div_heat_load_mw": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "pflux_div_heat_load_max_mw": InputVariable(
        data_structure.divertor_variables, float, range=(0.1, 20.0)
    ),
    "hot_sepdist": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "hotcell_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "p_blkt_coolant_pump_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "p_div_coolant_pump_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "p_fw_coolant_pump_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "htpmw_ife": InputVariable(data_structure.ife_variables, float, range=(0.0, 1000.0)),
    "p_shld_coolant_pump_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "hw_storage_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "hw_storage_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hw_storage_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_smelter_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "ilw_smelter_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_smelter_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_storage_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "ilw_storage_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_storage_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ind_plasma_internal_norm": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 10.0)
    ),
    "pres_vv_chamber_dwell_start": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 10000.0)
    ),
    "temp_blkt_coolant_in": InputVariable(
        data_structure.fwbs_variables, float, range=(200.0, 600.0)
    ),
    "inlet_temp_liq": InputVariable(
        data_structure.fwbs_variables, float, range=(508.0, 1500.0)
    ),
    "iotabar": InputVariable(
        data_structure.stellarator_variables, float, range=(0.1, 10.0)
    ),
    "j_tf_bus": InputVariable(
        data_structure.tfcoil_variables, float, range=(10000.0, 100000000.0)
    ),
    "kappa": InputVariable(data_structure.physics_variables, float, range=(0.99, 5.0)),
    "kappa95": InputVariable(data_structure.physics_variables, float, range=(0.99, 5.0)),
    "layer_ins": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 0.1)),
    "f_dr_dz_cs_turn": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 5.0)
    ),
    "len_fw_channel": InputVariable(
        data_structure.fwbs_variables, float, range=(0.001, 1000.0)
    ),
    "len_tf_bus": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.01, 1000.0)
    ),
    "lhat": InputVariable(data_structure.reinke_variables, float, range=(1.0, 15.0)),
    "f_blkt_li6_enrichment": InputVariable(
        data_structure.fwbs_variables, float, range=(7.4, 100.0)
    ),
    "life_dpa": InputVariable(data_structure.cost_variables, float, range=(10.0, 100.0)),
    "llw_storage_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "llw_storage_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "llw_storage_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "m_s_limit": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "magnet_pulse_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "magnet_pulse_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_pulse_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_trains_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "magnet_trains_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_trains_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maint_cont_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "maint_cont_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maint_cont_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maintenance_fwbs": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "maintenance_gen": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "max_gyrotron_frequency": InputVariable(
        data_structure.stellarator_variables,
        float,
        range=(1000000000.0, 100000000000000.0),
    ),
    "max_vv_stress": InputVariable(
        data_structure.tfcoil_variables, float, range=(100000.0, 500000000.0)
    ),
    "maxpoloidalpower": InputVariable(
        data_structure.pf_power_variables, float, range=(0.0, 2000.0)
    ),
    "pflux_fw_rad_max": InputVariable(
        data_structure.constraint_variables, float, range=(0.1, 10.0)
    ),
    "mbvfac": InputVariable(data_structure.buildings_variables, float, range=(0.9, 3.0)),
    "mcdriv": InputVariable(data_structure.ife_variables, float, range=(0.1, 10.0)),
    "mvalim": InputVariable(
        data_structure.constraint_variables, float, range=(0.0, 1000.0)
    ),
    "n_cycle_min": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.0, 100000000.0)
    ),
    "n_tf_coils": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "n_tf_coil_turns": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "nbi_sys_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "nbi_sys_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "dx_beam_shield": InputVariable(
        data_structure.current_drive_variables, float, range=(0.01, 0.5)
    ),
    "f_p_beam_shine_through_max": InputVariable(
        data_structure.constraint_variables, float, range=(1e-20, 0.1)
    ),
    "nd_plasma_pedestal_electron": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1e21)
    ),
    "nd_plasma_separatrix_electron": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1e21)
    ),
    "nflutfmax": InputVariable(
        data_structure.constraint_variables, float, range=(1e20, 1e24)
    ),
    "oacdcp": InputVariable(
        data_structure.tfcoil_variables, float, range=(10000.0, 1000000000.0)
    ),
    "f_a_cs_turn_steel": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.001, 0.999)
    ),
    "f_z_cs_tf_internal": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 2.0)
    ),
    "outgasfactor": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "outgasindex": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "temp_blkt_coolant_out": InputVariable(
        data_structure.fwbs_variables, float, range=(450.0, 900.0)
    ),
    "outlet_temp_liq": InputVariable(
        data_structure.fwbs_variables, float, range=(508.0, 1500.0)
    ),
    "p_he": InputVariable(
        data_structure.primary_pumping_variables, float, range=(0.0, 100000000.0)
    ),
    "paris_coefficient": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1e-20, 10.0)
    ),
    "paris_power_law": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1.0, 10.0)
    ),
    "pres_vv_chamber_base": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-08, 0.001)
    ),
    "p_plasma_separatrix_min_mw": InputVariable(
        data_structure.constraint_variables, float, range=(0.1, 1000.0)
    ),
    "pdrive": InputVariable(
        data_structure.ife_variables, float, range=(1000000.0, 200000000.0)
    ),
    "pfbldgm3": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "rho_pf_coil": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 0.0001)
    ),
    "pfusife": InputVariable(data_structure.ife_variables, float, range=(0.0, 10000.0)),
    "p_hcd_primary_extra_heat_mw": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "p_hcd_secondary_extra_heat_mw": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "pibv": InputVariable(
        data_structure.buildings_variables, float, range=(1000.0, 100000.0)
    ),
    "pifecr": InputVariable(data_structure.ife_variables, float, range=(0.0, 100.0)),
    "p_hcd_injected_max": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "p_hcd_secondary_injected_mw": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "plasma_res_factor": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "plasma_square": InputVariable(
        data_structure.physics_variables, float, range=(-5.0, 5.0)
    ),
    "plleni": InputVariable(data_structure.build_variables, float, range=(0.1, 10.0)),
    "plleno": InputVariable(data_structure.build_variables, float, range=(0.1, 10.0)),
    "plsepi": InputVariable(data_structure.build_variables, float, range=(0.1, 10.0)),
    "plsepo": InputVariable(data_structure.build_variables, float, range=(0.1, 10.0)),
    "p_plant_electric_net_required_mw": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 10000.0)
    ),
    "pnuc_fw_ratio_dcll": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "poisson_al": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "poisson_copper": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "poisson_steel": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "p_fusion_total_max_mw": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 10000.0)
    ),
    "pres_div_chamber_burn": InputVariable(
        data_structure.vacuum_variables, float, range=(0.0, 10.0)
    ),
    "pres_fw_coolant": InputVariable(
        data_structure.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "prn1": InputVariable(data_structure.divertor_variables, float, range=(0.0, 1.0)),
    "psepbqarmax": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 50.0)
    ),
    "pseprmax": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 60.0)
    ),
    "ptargf": InputVariable(data_structure.ife_variables, float, range=(0.1, 100.0)),
    "temp_cp_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(4.0, 573.15)
    ),
    "ptfnucmax": InputVariable(
        data_structure.constraint_variables, float, range=(1e-06, 1.0)
    ),
    "pulsetimings": InputVariable(
        data_structure.times_variables, float, range=(0.0, 1.0)
    ),
    "f_a_vac_pump_port_plasma_surface": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1.0)
    ),
    "f_volflow_vac_pumps_impedance": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1.0)
    ),
    "volflow_vac_pumps_max": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "molflow_vac_pumps": InputVariable(
        data_structure.vacuum_variables, float, range=(0.0, 1e30)
    ),
    "pflux_plant_floor_electric": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "q0": InputVariable(data_structure.physics_variables, float, range=(0.01, 20.0)),
    "q95": InputVariable(data_structure.physics_variables, float, range=(1.0, 50.0)),
    "q95_fixed": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 50.0)
    ),
    "qnty_sfty_fac": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "qnuc": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1000000.0)),
    "r_cp_top": InputVariable(
        data_structure.build_variables, float, range=(0.001, 10.0)
    ),
    "rad_fraction_sol": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 1.0)
    ),
    "radius_fw_channel": InputVariable(
        data_structure.fwbs_variables, float, range=(0.001, 0.5)
    ),
    "outgrat_fw": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-10, 1e-06)
    ),
    "rbrt": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "rbvfac": InputVariable(data_structure.buildings_variables, float, range=(0.9, 3.0)),
    "rbwt": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "radius_cp_coolant_channel": InputVariable(
        data_structure.tfcoil_variables, float, range=(1e-06, 1.0)
    ),
    "reactor_clrnc": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "reactor_fndtn_thk": InputVariable(
        data_structure.buildings_variables, float, range=(0.25, 25.0)
    ),
    "reactor_hall_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "reactor_hall_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "reactor_hall_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "reactor_roof_thk": InputVariable(
        data_structure.buildings_variables, float, range=(0.25, 25.0)
    ),
    "reactor_wall_thk": InputVariable(
        data_structure.buildings_variables, float, range=(0.25, 25.0)
    ),
    "dx_hts_tape_rebco": InputVariable(
        data_structure.rebco_variables,
        float,
        range=(1e-08, 0.0001),
        additional_actions=lambda _n, rt, _i, _c: rt <= 1e-6
        or warn(
            (
                "the relationship between REBCO layer thickness and current density is not linear."
                "REBCO layer thicknesses > 1um should be considered an aggressive extrapolation of"
                "current HTS technology and any results must be considered speculative."
            ),
            stacklevel=1,
        ),
    ),
    "redun_vacp": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 100.0)
    ),
    "residual_sig_hoop": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.0, 1000000000.0)
    ),
    "rho_tf_bus": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1e-05)
    ),
    "rho_tf_joints": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.01)
    ),
    "radius_plasma_pedestal_density_norm": InputVariable(
        data_structure.physics_variables, float, range=(0.01, 1.0)
    ),
    "radius_plasma_pedestal_temp_norm": InputVariable(
        data_structure.physics_variables, float, range=(0.01, 1.0)
    ),
    "rhopfbus": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 1e-05)
    ),
    "rinboard": InputVariable(data_structure.build_variables, float, range=(0.1, 10.0)),
    "ripple_b_tf_plasma_edge_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.1, 100.0)
    ),
    "rmajor": InputVariable(data_structure.physics_variables, float, range=(0.1, 50.0)),
    "robotics_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "robotics_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "robotics_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "roughness_fw_channel": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 0.01)
    ),
    "dr_pf_tf_outboard_out_offset": InputVariable(
        data_structure.pfcoil_variables, float, range=(-3.0, 3.0)
    ),
    "row": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "dr_pf_cs_middle_offset": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 3.0)
    ),
    "rpf2": InputVariable(data_structure.pfcoil_variables, float, range=(-3.0, 3.0)),
    "rrin": InputVariable(data_structure.ife_variables, float, range=(0.1, 50.0)),
    "rrmax": InputVariable(data_structure.ife_variables, float, range=(1.0, 50.0)),
    "rrr_tf_cu": InputVariable(
        data_structure.tfcoil_variables, float, range=(1.0, 1000.0)
    ),
    "rxcl": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "sec_buildings_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "sec_buildings_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "sec_buildings_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "sf_fast_fracture": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1.0, 10.0)
    ),
    "sf_radial_crack": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1.0, 10.0)
    ),
    "sf_vertical_crack": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1.0, 10.0)
    ),
    "shdr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "shdzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "shdzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "shear": InputVariable(
        data_structure.stellarator_variables, float, range=(0.1, 10.0)
    ),
    "dz_shld_lower": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dz_shld_upper": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "shmf": InputVariable(data_structure.buildings_variables, float, range=(0.0, 1.0)),
    "shov": InputVariable(
        data_structure.buildings_variables, float, range=(1000.0, 1000000.0)
    ),
    "sig_tf_case_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "sig_tf_wp_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "sigallpc": InputVariable(
        data_structure.build_variables, float, range=(0.0, 1000000000.0)
    ),
    "sigpfcalw": InputVariable(
        data_structure.pfcoil_variables, float, range=(1.0, 1000.0)
    ),
    "sigpfcf": InputVariable(data_structure.pfcoil_variables, float, range=(0.1, 1.0)),
    "sombdr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "somtdr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "staff_buildings_area": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "staff_buildings_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "startupratio": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 10.0)
    ),
    "stcl": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "str_cs_con_res": InputVariable(
        data_structure.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_pf_con_res": InputVariable(
        data_structure.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_tf_con_res": InputVariable(
        data_structure.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_wp_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.3)
    ),
    "t_plant_pulse_dwell": InputVariable(
        data_structure.times_variables, float, range=(0.0, 100000000.0)
    ),
    "t_plant_pulse_burn": InputVariable(
        data_structure.times_variables, float, range=(0.0, 100000000.0)
    ),
    "t_burn_min": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 1000000.0)
    ),
    "dx_tf_turn_cable_space_general": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "t_crack_radial": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1e-05, 1.0)
    ),
    "t_crack_vertical": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1e-05, 1.0)
    ),
    "t_crit_nbti": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 15.0)
    ),
    "t_plant_pulse_plasma_current_ramp_up": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_plant_pulse_fusion_ramp": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_in_bb": InputVariable(
        data_structure.primary_pumping_variables, float, range=(200.0, 1000.0)
    ),
    "t_tf_quench_detection": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "t_out_bb": InputVariable(
        data_structure.primary_pumping_variables, float, range=(200.0, 1000.0)
    ),
    "t_plant_pulse_coil_precharge": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_plant_pulse_plasma_current_ramp_down": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "dr_cs_turn_conduit": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.001, 1.0)
    ),
    "dz_cs_turn_conduit": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.001, 1.0)
    ),
    "dx_tf_turn_general": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "t_turn_tf_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dx_hts_tape_total": InputVariable(
        data_structure.rebco_variables, float, range=(0.0, 0.1)
    ),
    "dr_hts_tape": InputVariable(
        data_structure.rebco_variables, float, range=(0.0, 0.1)
    ),
    "tauee_in": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 100.0)
    ),
    "t_plasma_energy_confinement_max": InputVariable(
        data_structure.physics_variables, float, range=(0.1, 100.0)
    ),
    "tauratio": InputVariable(
        data_structure.physics_variables, float, range=(0.1, 100.0)
    ),
    "n_beam_decay_lengths_core_required": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 10.0)
    ),
    "tbeta": InputVariable(data_structure.physics_variables, float, range=(0.0, 4.0)),
    "t_blkt_replace_yrs": InputVariable(
        data_structure.cost_variables, float, range=(0.01, 2.0)
    ),
    "tbrmin": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 2.0)
    ),
    "tcomrepl": InputVariable(data_structure.cost_variables, float, range=(0.01, 2.0)),
    "temp_cp_coolant_inlet": InputVariable(
        data_structure.tfcoil_variables, float, range=(4.0, 373.15)
    ),
    "tcritsc": InputVariable(data_structure.tfcoil_variables, float, range=(1.0, 300.0)),
    "t_cycle_min": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 2000000.0)
    ),
    "tdiv": InputVariable(data_structure.divertor_variables, float, range=(0.1, 100.0)),
    "t_div_replace_yrs": InputVariable(
        data_structure.cost_variables, float, range=(0.01, 2.0)
    ),
    "t_tf_superconductor_quench": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.1, 100.0)
    ),
    "temp_plasma_electron_vol_avg_kev": InputVariable(
        data_structure.physics_variables, float, range=(1.0, 200.0)
    ),
    "te0_ecrh_achievable": InputVariable(
        data_structure.stellarator_variables, float, range=(1.0, 1000.0)
    ),
    "temp_cp_average": InputVariable(
        data_structure.tfcoil_variables, float, range=(4.0, 573.15)
    ),
    "temp_fw_coolant_in": InputVariable(
        data_structure.fwbs_variables, float, range=(300.0, 1500.0)
    ),
    "temp_fw_coolant_out": InputVariable(
        data_structure.fwbs_variables, float, range=(300.0, 1500.0)
    ),
    "temp_fw_max": InputVariable(
        data_structure.fwbs_variables, float, range=(500.0, 2000.0)
    ),
    "temp_plasma_pedestal_kev": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 20.0)
    ),
    "temp_plasma_separatrix_kev": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 20.0)
    ),
    "tfcbv": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "dx_tf_wp_insertion_gap": InputVariable(
        data_structure.tfcoil_variables, float, range=(1e-10, 0.1)
    ),
    "f_dr_tf_outboard_inboard": InputVariable(
        data_structure.build_variables, float, range=(0.2, 5.0)
    ),
    "tftmp": InputVariable(data_structure.tfcoil_variables, float, range=(0.01, 293.0)),
    "tgain": InputVariable(data_structure.ife_variables, float, range=(1.0, 500.0)),
    "th_joint_contact": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "theta1_coil": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.1, 60.0)
    ),
    "theta1_vv": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.1, 60.0)
    ),
    "dx_tf_turn_insulation": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "dr_tf_nose_case": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dz_shld_thermal": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dx_tf_turn_steel": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "temp_plasma_ion_vol_avg_kev": InputVariable(
        data_structure.physics_variables, float, range=(5.0, 50.0)
    ),
    "dx_tf_wp_insulation": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "life_plant": InputVariable(
        data_structure.cost_variables, float, range=(1.0, 100.0)
    ),
    "tmain": InputVariable(data_structure.cost_variables, float, range=(0.0, 100.0)),
    "tmargmin": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 20.0)),
    "temp_cs_superconductor_margin_min": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 20.0)
    ),
    "temp_tf_superconductor_margin_min": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 20.0)
    ),
    "temp_croco_quench_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(4.0, 1000.0)
    ),
    "temp_tf_conductor_quench_max": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1000.0)
    ),
    "temp_tf_cryo": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.01, 293.0)
    ),
    "temp_vv_chamber_gas_burn_end": InputVariable(
        data_structure.vacuum_variables, float, range=(1.0, 1000.0)
    ),
    "i_t_current_ramp_up": InputVariable(
        data_structure.times_variables, int, choices=[0, 1]
    ),
    "transp_clrnc": InputVariable(
        data_structure.buildings_variables, float, range=(0.0, 10.0)
    ),
    "f_temp_plasma_ion_electron": InputVariable(
        data_structure.physics_variables, float, range=(0.0, 2.0)
    ),
    "trcl": InputVariable(data_structure.buildings_variables, float, range=(0.0, 10.0)),
    "triang": InputVariable(data_structure.physics_variables, float, range=(-1.0, 1.0)),
    "triang95": InputVariable(data_structure.physics_variables, float, range=(0.0, 1.0)),
    "p_tritium_plant_electric_mw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "triv": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "turbine_hall_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "turbine_hall_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "turbine_hall_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "tw_storage_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "tw_storage_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "tw_storage_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "u_unplanned_cp": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "ucblbe": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucblbreed": InputVariable(
        data_structure.cost_variables, float, range=(1.0, 1000.0)
    ),
    "ucblli": InputVariable(data_structure.cost_variables, float, range=(10.0, 10000.0)),
    "ucblli2o": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 10000.0)
    ),
    "ucbllipb": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 10000.0)
    ),
    "ucblss": InputVariable(data_structure.cost_variables, float, range=(10.0, 1000.0)),
    "ucblvd": InputVariable(data_structure.cost_variables, float, range=(100.0, 1000.0)),
    "ucbus": InputVariable(data_structure.cost_variables, float, range=(0.01, 10.0)),
    "uccarb": InputVariable(data_structure.ife_variables, float, range=(10.0, 1000.0)),
    "uccase": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucconc": InputVariable(data_structure.ife_variables, float, range=(0.1, 1000.0)),
    "uccpcl1": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "uccpclb": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "uccry": InputVariable(
        data_structure.cost_variables, float, range=(10000.0, 1000000.0)
    ),
    "uccryo": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "uccu": InputVariable(data_structure.cost_variables, float, range=(10.0, 100.0)),
    "ucdiv": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 10000000.0)
    ),
    "ucech": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "ucf1": InputVariable(
        data_structure.cost_variables, float, range=(1000000.0, 50000000.0)
    ),
    "ucflib": InputVariable(data_structure.ife_variables, float, range=(10.0, 1000.0)),
    "ucfnc": InputVariable(data_structure.cost_variables, float, range=(10.0, 100.0)),
    "ucfuel": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "uche3": InputVariable(
        data_structure.cost_variables, float, range=(100000.0, 10000000.0)
    ),
    "uchrs": InputVariable(
        data_structure.cost_variables, float, range=(10000000.0, 500000000.0)
    ),
    "uciac": InputVariable(
        data_structure.cost_variables, float, range=(10000000.0, 1000000000.0)
    ),
    "ucich": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "uclh": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "ucme": InputVariable(
        data_structure.cost_variables, float, range=(10000000.0, 1000000000.0)
    ),
    "ucmisc": InputVariable(
        data_structure.cost_variables, float, range=(10000000.0, 50000000.0)
    ),
    "ucnbi": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "ucpens": InputVariable(data_structure.cost_variables, float, range=(1.0, 100.0)),
    "ucpfb": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucpfbk": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000.0)
    ),
    "ucpfbs": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 10000.0)
    ),
    "ucpfcb": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000.0)
    ),
    "ucpfdr1": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucpfic": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000.0)
    ),
    "ucpfps": InputVariable(
        data_structure.cost_variables, float, range=(1000.0, 100000.0)
    ),
    "ucrb": InputVariable(data_structure.cost_variables, float, range=(100.0, 1000.0)),
    "ucshld": InputVariable(data_structure.cost_variables, float, range=(1.0, 100.0)),
    "uctarg": InputVariable(data_structure.ife_variables, float, range=(0.1, 1000.0)),
    "uctfbr": InputVariable(data_structure.cost_variables, float, range=(1.0, 10.0)),
    "uctfbus": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "uctfps": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "uctfsw": InputVariable(data_structure.cost_variables, float, range=(0.1, 10.0)),
    "ucwindpf": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 1000.0)
    ),
    "ucwindtf": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 1000.0)
    ),
    "uubop": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uucd": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uudiv": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uufuel": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uufw": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uumag": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "uuves": InputVariable(data_structure.cost_variables, float, range=(0.005, 0.1)),
    "v1dr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v1dzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v1dzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v2dr": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v2dzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v2dzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 10.0)),
    "v3dr": InputVariable(data_structure.ife_variables, float, range=(0.0, 50.0)),
    "v3dzl": InputVariable(data_structure.ife_variables, float, range=(0.0, 30.0)),
    "v3dzu": InputVariable(data_structure.ife_variables, float, range=(0.0, 30.0)),
    "vachtmw": InputVariable(
        data_structure.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "vel_cp_coolant_midplane": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.001, 100.0)
    ),
    "v_tf_coil_dump_quench_max_kv": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "f_a_blkt_cooling_channels": InputVariable(
        data_structure.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "vfcblkt": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "f_a_cs_void": InputVariable(
        data_structure.pfcoil_variables, float, range=(0.0, 1.0)
    ),
    "vfpblkt": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "vfshld": InputVariable(data_structure.fwbs_variables, float, range=(0.0, 1.0)),
    "f_a_tf_turn_cable_space_extra_void": InputVariable(
        data_structure.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dz_shld_vv_gap": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dz_xpoint_divertor": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "dz_fw_plasma_gap": InputVariable(
        data_structure.build_variables, float, range=(0.0, 10.0)
    ),
    "pflux_fw_neutron_max_mw": InputVariable(
        data_structure.constraint_variables, float, range=(0.001, 50.0)
    ),
    "walker_coefficient": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.1, 10.0)
    ),
    "wallpf": InputVariable(data_structure.fwbs_variables, float, range=(1.0, 2.0)),
    "warm_shop_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "warm_shop_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "warm_shop_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "water_buildings_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "water_buildings_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "water_buildings_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "watertemp": InputVariable(
        data_structure.water_usage_variables, float, range=(0.0, 25.0)
    ),
    "wgt": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "wgt2": InputVariable(
        data_structure.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "windspeed": InputVariable(
        data_structure.water_usage_variables, float, range=(0.0, 10.0)
    ),
    "workshop_h": InputVariable(
        data_structure.buildings_variables, float, range=(1.0, 100.0)
    ),
    "workshop_l": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "workshop_w": InputVariable(
        data_structure.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "wsvfac": InputVariable(data_structure.buildings_variables, float, range=(0.9, 3.0)),
    "xi_ebw": InputVariable(
        data_structure.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "xpertin": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "zeff_max": InputVariable(
        data_structure.constraint_variables, float, range=(1.0, 10.0)
    ),
    "blktmodel": InputVariable(data_structure.fwbs_variables, int, choices=[0, 1]),
    "blkttype": InputVariable(data_structure.fwbs_variables, int, choices=[1, 2, 3]),
    "breedmat": InputVariable(data_structure.fwbs_variables, int, choices=[1, 2, 3]),
    "ccl0_ma": InputVariable(data_structure.pfcoil_variables, float, array=True),
    "ccls_ma": InputVariable(data_structure.pfcoil_variables, float, array=True),
    "cfind": InputVariable(data_structure.cost_variables, float, array=True),
    "i_blkt_coolant_type": InputVariable(
        data_structure.fwbs_variables, int, choices=[1, 2]
    ),
    "coppera_m2_max": InputVariable(
        data_structure.rebco_variables, float, range=(1.0e6, 1.0e10)
    ),
    "cost_model": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "i_vac_pump_dwell": InputVariable(
        data_structure.vacuum_variables, int, choices=[0, 1, 2]
    ),
    "i_fw_blkt_vv_shape": InputVariable(
        data_structure.fwbs_variables, int, range=(1, 2)
    ),
    "hcdportsize": InputVariable(data_structure.fwbs_variables, int, range=(1, 2)),
    "i_blkt_liquid_breeder_type": InputVariable(
        data_structure.fwbs_variables, int, choices=[0, 1]
    ),
    "i_beta_component": InputVariable(
        data_structure.physics_variables, int, range=(0, 3)
    ),
    "i_beta_fast_alpha": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_blanket_type": InputVariable(
        data_structure.fwbs_variables, int, choices=[1, 3, 5]
    ),
    "i_bldgs_size": InputVariable(
        data_structure.buildings_variables, int, choices=[0, 1]
    ),
    "i_bldgs_v": InputVariable(data_structure.buildings_variables, int, choices=[0, 1]),
    "i_blkt_inboard": InputVariable(data_structure.fwbs_variables, int, choices=[0, 1]),
    "i_bootstrap_current": InputVariable(
        data_structure.physics_variables, int, range=(0, 13)
    ),
    "i_cp_joints": InputVariable(data_structure.tfcoil_variables, int, choices=[0, 1]),
    "i_cp_lifetime": InputVariable(data_structure.cost_variables, int, range=(0, 3)),
    "i_cs_precomp": InputVariable(data_structure.build_variables, int, choices=[0, 1]),
    "i_cs_stress": InputVariable(data_structure.pfcoil_variables, int, choices=[0, 1]),
    "i_density_limit": InputVariable(
        data_structure.physics_variables, int, range=(1, 8)
    ),
    "i_diamagnetic_current": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1, 2]
    ),
    "i_div_heat_load": InputVariable(
        data_structure.divertor_variables, int, choices=[0, 1, 2]
    ),
    "i_l_h_threshold": InputVariable(
        data_structure.physics_variables, int, range=(1, 21)
    ),
    "i_pf_current": InputVariable(
        data_structure.pfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_pfirsch_schluter_current": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_plasma_current": InputVariable(
        data_structure.physics_variables, int, range=(1, 9)
    ),
    "i_plasma_geometry": InputVariable(
        data_structure.physics_variables, int, range=(0, 11)
    ),
    "i_plasma_shape": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_plasma_wall_gap": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_pulsed_plant": InputVariable(data_structure.pulse_variables, int, choices=[0, 1]),
    "i_q95_fixed": InputVariable(
        data_structure.constraint_variables, int, choices=[0, 1]
    ),
    "i_r_cp_top": InputVariable(data_structure.build_variables, int, choices=[0, 1, 2]),
    "i_rad_loss": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1, 2]
    ),
    "i_shield_mat": InputVariable(data_structure.fwbs_variables, int, choices=[0, 1]),
    "i_single_null": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_str_wp": InputVariable(data_structure.tfcoil_variables, int, choices=[0, 1]),
    "i_r_pf_outside_tf_placement": InputVariable(
        data_structure.pfcoil_variables, int, choices=[0, 1]
    ),
    "i_tf_bucking": InputVariable(data_structure.tfcoil_variables, int, range=(0, 3)),
    "i_tf_case_geom": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1]
    ),
    "i_tf_cond_eyoung_axial": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_tf_cond_eyoung_trans": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1]
    ),
    "i_tf_sc_mat": InputVariable(data_structure.tfcoil_variables, int, range=(1, 9)),
    "i_tf_shape": InputVariable(data_structure.tfcoil_variables, int, choices=[0, 1, 2]),
    "i_tf_stress_model": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_tf_sup": InputVariable(data_structure.tfcoil_variables, int, choices=[0, 1, 2]),
    "i_tf_tresca": InputVariable(data_structure.tfcoil_variables, int, choices=[0, 1]),
    "i_tf_turns_integer": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1]
    ),
    "i_tf_wp_geom": InputVariable(
        data_structure.tfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_plant_availability": InputVariable(
        data_structure.cost_variables, int, range=(0, 3)
    ),
    "ibkt_life": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "i_blkt_dual_coolant": InputVariable(
        data_structure.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "i_hcd_primary": InputVariable(
        data_structure.current_drive_variables, int, range=(1, 13)
    ),
    "i_hcd_secondary": InputVariable(
        data_structure.current_drive_variables, int, range=(0, 13)
    ),
    "i_blkt_liquid_breeder_channel_type": InputVariable(
        data_structure.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "ife": InputVariable(data_structure.ife_variables, int, choices=[0, 1]),
    "ifedrv": InputVariable(data_structure.ife_variables, int, range=(-1, 3)),
    "ifetyp": InputVariable(data_structure.ife_variables, int, range=(0, 4)),
    "ifueltyp": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "i_plasma_ignited": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_blkt_module_segmentation": InputVariable(
        data_structure.fwbs_variables, int, choices=[0, 1]
    ),
    "inuclear": InputVariable(data_structure.fwbs_variables, int, choices=[0, 1]),
    "iohcl": InputVariable(data_structure.build_variables, int, choices=[0, 1]),
    "i_plasma_pedestal": InputVariable(
        data_structure.physics_variables, int, choices=[0, 1]
    ),
    "i_pf_conductor": InputVariable(
        data_structure.pfcoil_variables, int, choices=[0, 1]
    ),
    "ipnet": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "ipowerflow": InputVariable(
        data_structure.heat_transport_variables, int, choices=[0, 1]
    ),
    "i_shld_primary_heat": InputVariable(
        data_structure.heat_transport_variables, int, choices=[0, 1]
    ),
    "i_beta_norm_max": InputVariable(
        data_structure.physics_variables, int, range=(0, 5)
    ),
    "i_ind_plasma_internal_norm": InputVariable(
        data_structure.physics_variables, int, range=(0, 2)
    ),
    "i_alphaj": InputVariable(data_structure.physics_variables, int, range=(0, 1)),
    "i_fw_blkt_shared_coolant": InputVariable(
        data_structure.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "ireactor": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "irefprop": InputVariable(data_structure.fwbs_variables, int, choices=[0, 1]),
    "i_hcd_calculations": InputVariable(
        data_structure.current_drive_variables, int, choices=[0, 1]
    ),
    "i_pf_energy_storage_source": InputVariable(
        data_structure.pf_power_variables, int, range=(1, 3)
    ),
    "istell": InputVariable(data_structure.stellarator_variables, int, range=(0, 6)),
    "isthtr": InputVariable(data_structure.stellarator_variables, int, range=(1, 3)),
    "istore": InputVariable(data_structure.pulse_variables, int, range=(1, 3)),
    "i_cs_superconductor": InputVariable(
        data_structure.pfcoil_variables, int, range=(1, 9)
    ),
    "i_pf_superconductor": InputVariable(
        data_structure.pfcoil_variables, int, range=(1, 9)
    ),
    "itart": InputVariable(data_structure.physics_variables, int, choices=[0, 1]),
    "itartpf": InputVariable(data_structure.physics_variables, int, choices=[0, 1]),
    "itcycl": InputVariable(data_structure.pulse_variables, int, range=(1, 3)),
    "i_pflux_fw_neutron": InputVariable(
        data_structure.physics_variables, int, range=(1, 2)
    ),
    "lsa": InputVariable(data_structure.cost_variables, int, range=(1, 4)),
    "m_res": InputVariable(data_structure.stellarator_variables, int, range=(1, 10)),
    "n_tf_wp_layers": InputVariable(
        data_structure.tfcoil_variables, int, range=(1, 100)
    ),
    "n_liq_recirc": InputVariable(data_structure.fwbs_variables, int, range=(1, 50)),
    "n_tf_wp_pancakes": InputVariable(
        data_structure.tfcoil_variables, int, range=(1, 100)
    ),
    "n_rad_per_layer": InputVariable(
        data_structure.tfcoil_variables, int, range=(1, 500)
    ),
    "n_res": InputVariable(data_structure.stellarator_variables, int, range=(3, 6)),
    "n_tf_graded_layers": InputVariable(
        data_structure.tfcoil_variables, int, range=(1, 20)
    ),
    "n_tf_joints": InputVariable(data_structure.tfcoil_variables, int, range=(1, 50)),
    "n_tf_joints_contact": InputVariable(
        data_structure.tfcoil_variables, int, range=(1, 50)
    ),
    "n_blkt_inboard_modules_poloidal": InputVariable(
        data_structure.fwbs_variables, int, range=(1, 16)
    ),
    "n_blkt_outboard_modules_poloidal": InputVariable(
        data_structure.fwbs_variables, int, range=(1, 16)
    ),
    "n_blkt_inboard_modules_toroidal": InputVariable(
        data_structure.fwbs_variables, int, range=(8, 96)
    ),
    "n_blkt_outboard_modules_toroidal": InputVariable(
        data_structure.fwbs_variables, int, range=(8, 96)
    ),
    "npdiv": InputVariable(data_structure.fwbs_variables, int, range=(0, 4)),
    "nphcdin": InputVariable(data_structure.fwbs_variables, int, range=(0, 4)),
    "nphcdout": InputVariable(data_structure.fwbs_variables, int, range=(0, 4)),
    "i_vacuum_pump_type": InputVariable(
        data_structure.vacuum_variables, int, choices=[0, 1]
    ),
    "num_rh_systems": InputVariable(data_structure.cost_variables, int, range=(1, 10)),
    "output_costs": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "i_p_coolant_pumping": InputVariable(
        data_structure.fwbs_variables, int, range=(0, 3)
    ),
    "reinke_mode": InputVariable(data_structure.reinke_variables, int, choices=[0, 1]),
    "scan_dim": InputVariable(data_structure.scan_variables, int, range=(1, 2)),
    "i_thermal_electric_conversion": InputVariable(
        data_structure.fwbs_variables, int, range=(0, 4)
    ),
    "secondary_cycle_liq": InputVariable(
        data_structure.fwbs_variables, int, range=(2, 4)
    ),
    "supercond_cost_model": InputVariable(
        data_structure.cost_variables, int, choices=[0, 1]
    ),
    "i_tf_inside_cs": InputVariable(data_structure.build_variables, int, choices=[0, 1]),
    "i_ecrh_wave_mode": InputVariable(
        data_structure.current_drive_variables, int, choices=[0, 1]
    ),
    "i_confinement_time": InputVariable(
        data_structure.physics_variables,
        int,
        choices=list(range(data_structure.physics_variables.N_CONFINEMENT_SCALINGS)),
    ),
    "quench_model": InputVariable(
        data_structure.tfcoil_variables, str, choices=["exponential", "linear"]
    ),
    "i_fw_coolant_type": InputVariable(
        data_structure.fwbs_variables, str, choices=["helium", "water"]
    ),
    "i_vacuum_pumping": InputVariable(
        data_structure.vacuum_variables, str, choices=["old", "simple"]
    ),
    "dcond": InputVariable(data_structure.tfcoil_variables, float, array=True),
    "c_pf_coil_turn_peak_input": InputVariable(
        data_structure.pfcoil_variables, float, array=True
    ),
    "i_pf_location": InputVariable(data_structure.pfcoil_variables, int, array=True),
    "n_pf_coils_in_group": InputVariable(
        data_structure.pfcoil_variables, int, array=True
    ),
    "n_cs_current_filaments": InputVariable(
        data_structure.pfcoil_variables, int, array=True
    ),
    "n_pf_coil_groups": InputVariable(
        data_structure.pfcoil_variables,
        int,
        range=(0, data_structure.pfcoil_variables.N_PF_GROUPS_MAX),
    ),
    "rref": InputVariable(data_structure.pfcoil_variables, float, array=True),
    "f_a_pf_coil_void": InputVariable(
        data_structure.pfcoil_variables, float, array=True
    ),
    "zref": InputVariable(data_structure.pfcoil_variables, float, array=True),
    "uchts": InputVariable(data_structure.cost_variables, float, array=True),
    "ucoam": InputVariable(data_structure.cost_variables, float, array=True),
    "ucsc": InputVariable(data_structure.cost_variables, float, array=True),
    "ucturb": InputVariable(data_structure.cost_variables, float, array=True),
    "ucwst": InputVariable(data_structure.cost_variables, float, array=True),
    "blmatf": InputVariable(data_structure.ife_variables, float, array=True),
    "chmatf": InputVariable(data_structure.ife_variables, float, array=True),
    "etave": InputVariable(data_structure.ife_variables, float, array=True),
    "fwmatf": InputVariable(data_structure.ife_variables, float, array=True),
    "gainve": InputVariable(data_structure.ife_variables, float, array=True),
    "shmatf": InputVariable(data_structure.ife_variables, float, array=True),
    "v1matf": InputVariable(data_structure.ife_variables, float, array=True),
    "v2matf": InputVariable(data_structure.ife_variables, float, array=True),
    "v3matf": InputVariable(data_structure.ife_variables, float, array=True),
    "isweep": InputVariable(
        data_structure.scan_variables,
        int,
        choices=range(data_structure.scan_variables.IPNSCNS + 1),
    ),
    "nsweep": InputVariable(
        data_structure.scan_variables,
        int,
        choices=range(1, data_structure.scan_variables.IPNSCNV + 1),
    ),
    "isweep_2": InputVariable(
        data_structure.scan_variables,
        int,
        choices=range(data_structure.scan_variables.IPNSCNS + 1),
    ),
    "nsweep_2": InputVariable(
        data_structure.scan_variables,
        int,
        choices=range(1, data_structure.scan_variables.IPNSCNV + 1),
    ),
    "sweep": InputVariable(data_structure.scan_variables, float, array=True),
    "sweep_2": InputVariable(data_structure.scan_variables, float, array=True),
    "impvardiv": InputVariable(
        data_structure.reinke_variables,
        int,
        choices=range(3, data_structure.impurity_radiation_module.N_IMPURITIES + 1),
    ),
    "j_pf_coil_wp_peak": InputVariable(
        data_structure.pfcoil_variables, float, array=True
    ),
    "ixc": InputVariable(
        None,
        int,
        range=(1, data_structure.numerics.ipnvars),
        additional_actions=_ixc_additional_actions,
        set_variable=False,
    ),
    "icc": InputVariable(
        None,
        int,
        range=(1, data_structure.numerics.ipeqns),
        additional_actions=_icc_additional_actions,
        set_variable=False,
    ),
    "force_vmcon_inequality_satisfication": InputVariable(
        data_structure.numerics,
        int,
        choices=(0, 1),
    ),
    "force_vmcon_inequality_tolerance": InputVariable(
        data_structure.numerics, float, range=(0.0, 1e10)
    ),
}


def parse_input_file():
    input_file = data_structure.global_variables.fileprefix

    input_file_path = Path("IN.DAT")
    if input_file != "":
        input_file_path = Path(input_file)

    with input_file_path.open("r") as f:
        lines = f.readlines()

    variables = {}

    for line_no, line in enumerate(lines, start=1):
        stripped_line = line.strip()

        # don't bother trying to process blank lines
        # or comment lines
        if stripped_line == "" or stripped_line[0] == "*":
            continue

        # matches (variable name, array index, value)
        # NOTE: array index is Fortran-based hence starts at 1.
        line_match = re.match(
            r"([a-zA-Z0-9_]+)(?:\(([0-9]+)\))?[ ]*=[ ]*([ +\-a-zA-Z0-9.,]+).*",
            stripped_line,
        )

        if line_match is None:
            error_msg = (
                f"Unable to parse line {line_no} of the input file ({stripped_line})"
            )
            raise ProcessValueError(error_msg)

        variable_name, array_index, variable_value = line_match.groups()
        variable_name = variable_name.lower()

        variable_config = INPUT_VARIABLES.get(variable_name)

        if variable_config is None:
            error_msg = (
                f"Unrecognised input '{variable_name}' at line {line_no} of input file."
            )
            raise ProcessValidationError(error_msg)

        # Validate the variable and also clean it (cast to correct type)
        # If the variable value (after the = sign) contains a ',' then it
        # defines the whole array so needs to be split down into its elements
        # and the parsed like an array defined as 'my_array(<index>) = <value>'
        if "," in variable_value:
            getattr(variable_config.module, variable_name)[:] = 0.0
            clean_variable_value = [
                validate_variable(
                    variable_name,
                    v.strip(),
                    str(i + 1),
                    variable_config,
                    line_no,
                )
                for i, v in enumerate(variable_value.split(","), start=1)
            ]
        else:
            clean_variable_value = validate_variable(
                variable_name,
                variable_value.strip(),
                array_index,
                variable_config,
                line_no,
            )

        # check if the target name in the module is different to the variable name
        # in the input file
        variable_name_in_module = variable_config.target_name or variable_name

        array_index_clean = None if array_index is None else int(array_index)

        if variable_config.set_variable:
            if isinstance(clean_variable_value, list):
                for idx, value in enumerate(clean_variable_value, start=1):
                    set_array_variable(
                        variable_name_in_module,
                        value,
                        idx,
                        variable_config,
                    )
            elif variable_config.array:
                set_array_variable(
                    variable_name_in_module,
                    clean_variable_value,
                    array_index_clean,
                    variable_config,
                )
            else:
                set_scalar_variable(
                    variable_name_in_module, clean_variable_value, variable_config
                )

        if variable_config.additional_actions is not None:
            # intentionally passing the variable name as in the input file,
            # not the module target name.
            variable_config.additional_actions(
                variable_name,
                clean_variable_value,
                array_index_clean,
                variable_config,
            )

        # add the variable to a dictionary indexed by the variable name (in the input file)
        variables[variable_name] = {
            "value": clean_variable_value,
            "index": array_index_clean,
            "config": variable_config,
        }

    return variables


def validate_variable(
    name: str,
    value: str,
    array_index: int | None,
    config: InputVariable,
    line_number: int,
) -> ValidInputTypes:
    """Validate an input.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param array_index: the array index of the variable in the input file.
    :param config: the config of the variable that describes how to validate and process it.
    :param line_number: line number of current line being parsed for error reporting.

    :returns: the input value with the correct type.
    """
    # check that if the variable should be an array, then an array index is provided
    # EXCEPT for if check_array is False. This should only be the case when parsing
    # entire arrays (e.g. my_array = 1,2,2,4,5) where there will be no array index.
    if array_index is None and config.array:
        error_msg = f"Expected '{name}' at line {line_number} to be an array."
        raise ProcessValidationError(error_msg)

    if array_index is not None and not config.array:
        error_msg = f"Not expecting '{name}' at line {line_number} to be an array."
        raise ProcessValidationError(error_msg)

    if config.type in [float, int]:
        value = value.lower().replace("d", "e")

    try:
        clean_value = config.type(value)
    except ValueError as e:
        error_msg = f"Cannot cast variable name '{name}' at line {line_number} to a {config.type} (value = {value})"
        raise ProcessValidationError(error_msg) from e

    if (
        config.range is not None
        and not config.range[0] <= clean_value <= config.range[1]
    ):
        error_msg = (
            f"Variable '{name}' at line {line_number} is not on the prescribed range"
            f" {config.range} (value = {value})"
        )
        raise ProcessValidationError(error_msg)

    if config.choices is not None and clean_value not in config.choices:
        error_msg = f"Variable '{name}' at line {line_number} is not one of {config.choices} (value = {value})"
        raise ProcessValidationError(error_msg)

    if config.additional_validation is not None:
        clean_value = config.additional_validation(
            name, clean_value, int(array_index), config
        )

    return clean_value


def set_scalar_variable(name: str, value: ValidInputTypes, config: InputVariable):
    """Set a scalar (not part of an array) variable in the `config.module`.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param config: the config of the variable that describes how to validate and process it.
    """
    current_value = getattr(config.module, name, ...)

    # use ... sentinel because None is probably a valid return from Fortran
    # and definately will be when moving to a Python data structure
    if current_value is ...:
        error_msg = (
            f"Fortran module '{config.module}' does not have a variable '{name}'."
        )
        raise ProcessValueError(error_msg)

    setattr(config.module, name, value)


def set_array_variable(name: str, value: str, array_index: int, config: InputVariable):
    """Set an array variable in the `config.module`.

    The way PROCESS input files are structured, each element of the array is provided on one line
    so this function just needs to set the `value` at `array_index` (-1) because of Fortran-based indexing.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param array_index: the array index of the variable in the input file.
    :param config: the config of the variable that describes how to validate and process it.
    """
    current_array = getattr(config.module, name, ...)
    shape = current_array.shape

    # use ... sentinel for same reason as above
    if current_array is ...:
        error_msg = f"Data structure '{config.module}' does not have an array '{name}'."
        raise ProcessValueError(error_msg)

    if config.type is str:
        raise NotImplementedError

    new_array = copy.deepcopy(current_array)

    if len(shape) > 1:
        new_array = new_array.T.ravel()

    new_array[array_index - 1] = value

    if len(shape) > 1:
        new_array = new_array.reshape(shape, order="F")

    setattr(config.module, name, new_array)
