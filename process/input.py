"""Handle parsing, validation, and actioning of a PROCESS input file (*IN.DAT)."""

import copy
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from warnings import warn

import process.data_structure as data_structure
import process.fortran as fortran
from process.constraints import ConstraintManager
from process.exceptions import ProcessValidationError, ProcessValueError
from process.utilities.f2py_string_patch import (
    f2py_compatible_to_string,
    string_to_f2py_compatible,
)

NumberType = int | float
ValidInputTypes = NumberType | str

DataTypes = (int, float, str)
"""Valid input variable types alias"""


def _ixc_additional_actions(_name, value: int, _array_index, _config):
    fortran.numerics.ixc[fortran.numerics.nvar.item()] = value
    fortran.numerics.nvar += 1


def _icc_additional_actions(_name, value: int, _array_index, _config):
    fortran.numerics.icc[fortran.numerics.n_constraints.item()] = value
    fortran.numerics.n_constraints += 1


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
    "runtitle": InputVariable(fortran.global_variables, str),
    "verbose": InputVariable(fortran.global_variables, int, choices=[0, 1]),
    "run_tests": InputVariable(fortran.global_variables, int, choices=[0, 1]),
    "ioptimz": InputVariable(fortran.numerics, int, choices=[1, -2]),
    "epsvmc": InputVariable(fortran.numerics, float, range=(0.0, 1.0)),
    "boundl": InputVariable(fortran.numerics, float, array=True),
    "boundu": InputVariable(fortran.numerics, float, array=True),
    "epsfcn": InputVariable(fortran.numerics, float, range=(0.0, 1.0)),
    "maxcal": InputVariable(fortran.global_variables, int, range=(0, 10000)),
    "minmax": InputVariable(fortran.numerics, int),
    "neqns": InputVariable(
        fortran.numerics, int, range=(1, ConstraintManager.num_constraints())
    ),
    "nineqns": InputVariable(
        fortran.numerics, int, range=(1, ConstraintManager.num_constraints())
    ),
    "alphaj": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "alphan": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "alphat": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "aspect": InputVariable(fortran.physics_variables, float, range=(1.001, 40.0)),
    "beamfus0": InputVariable(fortran.physics_variables, float, range=(0.01, 10.0)),
    "beta": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "beta_max": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "beta_min": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "betbm0": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "bt": InputVariable(fortran.physics_variables, float, range=(0.0, 30.0)),
    "burnup_in": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "radius_plasma_core_norm": InputVariable(
        fortran.impurity_radiation_module, float, range=(0.0, 1.0)
    ),
    "coreradiationfraction": InputVariable(
        fortran.impurity_radiation_module, float, range=(0.0, 1.0)
    ),
    "c_beta": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "csawth": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "f_vol_plasma": InputVariable(
        fortran.physics_variables, float, range=(0.001, 10.0)
    ),
    "f_r_conducting_wall": InputVariable(
        fortran.physics_variables, float, range=(1.0, 3.0)
    ),
    "dene": InputVariable(fortran.physics_variables, float, range=(1.0e18, 1.0e22)),
    "beta_norm_max": InputVariable(fortran.physics_variables, float, range=(0.0, 20.0)),
    "beta_poloidal_eps_max": InputVariable(
        fortran.physics_variables, float, range=(0.01, 10.0)
    ),
    "f_p_alpha_plasma_deposited": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_p_div_lower": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "f_deuterium": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "ffwal": InputVariable(fortran.physics_variables, float, range=(0.0, 10.0)),
    "fgwped": InputVariable(fortran.physics_variables, float, range=(-1.0, 5.0)),
    "fgwsep": InputVariable(fortran.physics_variables, float, range=(-1.0, 5.0)),
    "f_helium3": InputVariable(fortran.physics_variables, float, range=(-1.0, 5.0)),
    # TODO: does fimp require an additional range?
    "fimp": InputVariable(fortran.impurity_radiation_module, float, array=True),
    "fkzohm": InputVariable(fortran.physics_variables, float, range=(0.5, 2.0)),
    "fnesep": InputVariable(fortran.constraint_variables, float, range=(0.1, 20.0)),
    "abktflnc": InputVariable(data_structure.cost_variables, float, range=(0.1, 100.0)),
    "adivflnc": InputVariable(data_structure.cost_variables, float, range=(0.1, 100.0)),
    "admv": InputVariable(fortran.buildings_variables, float, range=(1.0e4, 1.0e6)),
    "airtemp": InputVariable(
        data_structure.water_usage_variables, float, range=(-15.0, 40.0)
    ),
    "alfapf": InputVariable(fortran.pfcoil_variables, float, range=(1e-12, 1.0)),
    "alstroh": InputVariable(
        fortran.pfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "amortization": InputVariable(
        data_structure.cost_variables, float, range=(1.0, 50.0)
    ),
    "anginc": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "aplasmin": InputVariable(fortran.build_variables, float, range=(0.01, 10.0)),
    "aux_build_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "aux_build_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "aux_build_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "auxcool_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "auxcool_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "auxcool_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "p_hcd_injected_min_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.01, 100.0)
    ),
    "avail_min": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "b_crit_upper_nbti": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 30.0)
    ),
    "p_plant_electric_base": InputVariable(
        fortran.heat_transport_variables, float, range=(1000000.0, 10000000000.0)
    ),
    "bcritsc": InputVariable(fortran.tfcoil_variables, float, range=(10.0, 50.0)),
    "bctmp": InputVariable(data_structure.pulse_variables, float, range=(1.0, 800.0)),
    "e_beam_kev": InputVariable(
        fortran.current_drive_variables, float, range=(1.0, 1000000.0)
    ),
    "dx_beam_duct": InputVariable(
        fortran.current_drive_variables, float, range=(0.001, 5.0)
    ),
    "deg_div_field_plate": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 360.0)
    ),
    "beta_poloidal_max": InputVariable(
        fortran.constraint_variables, float, range=(0.01, 2.0)
    ),
    "betai": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "betao": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 1.5707)
    ),
    "bigqmin": InputVariable(fortran.constraint_variables, float, range=(0.01, 100.0)),
    "bioshld_thk": InputVariable(
        fortran.buildings_variables, float, range=(0.25, 25.0)
    ),
    "bkt_life_csf": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.0, 1.0)
    ),
    "blbmith": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "blbmoth": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "blbpith": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "blbpoth": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "blbuith": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "blbuoth": InputVariable(fortran.build_variables, float, range=(0.0, 2.0)),
    "bldr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "bldrc": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "bldzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "bldzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "pres_blkt_coolant": InputVariable(
        fortran.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "blpressure_liq": InputVariable(
        fortran.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "b_cs_limit_max": InputVariable(
        fortran.pfcoil_variables, float, range=(0.01, 100.0)
    ),
    "bmn": InputVariable(fortran.stellarator_variables, float, range=(0.0001, 0.01)),
    "b_tf_inboard_max": InputVariable(
        fortran.constraint_variables, float, range=(0.1, 50.0)
    ),
    "f_c_plasma_bootstrap_max": InputVariable(
        fortran.current_drive_variables, float, range=(-0.999, 0.999)
    ),
    "f_c_plasma_bootstrap": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "breeder_f": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "breeder_multiplier": InputVariable(
        fortran.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "bz_channel_conduct_liq": InputVariable(
        fortran.fwbs_variables, float, range=(1e-06, 1000000.0)
    ),
    "dr_tf_plasma_case": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_dr_tf_plasma_case": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dx_tf_side_case_min": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "casths_fraction": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "cboot": InputVariable(fortran.current_drive_variables, float, range=(0.0, 10.0)),
    "cconfix": InputVariable(data_structure.cost_variables, float, range=(50.0, 200.0)),
    "cconshpf": InputVariable(
        data_structure.cost_variables, float, range=(50.0, 200.0)
    ),
    "cconshtf": InputVariable(
        data_structure.cost_variables, float, range=(50.0, 200.0)
    ),
    "cdriv0": InputVariable(fortran.ife_variables, float, range=(50.0, 500.0)),
    "cdriv1": InputVariable(fortran.ife_variables, float, range=(50.0, 500.0)),
    "cdriv2": InputVariable(fortran.ife_variables, float, range=(50.0, 500.0)),
    "cfactr": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "chdzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "chdzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "chemlab_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "chemlab_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "chemlab_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "chrad": InputVariable(fortran.ife_variables, float, range=(0.1, 20.0)),
    "cland": InputVariable(data_structure.cost_variables, float, range=(10.0, 100.0)),
    "clh2": InputVariable(fortran.buildings_variables, float, range=(0.0, 30.0)),
    "j_cs_flat_top_end": InputVariable(
        fortran.pfcoil_variables, float, range=(10000.0, 500000000.0)
    ),
    "conf_mag": InputVariable(data_structure.cost_variables, float, range=(0.9, 1.0)),
    "control_buildings_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "control_buildings_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "control_buildings_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "conv": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "coolp": InputVariable(
        fortran.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "copper_rrr": InputVariable(
        data_structure.rebco_variables, float, range=(1.0, 10000.0)
    ),
    "copper_thick": InputVariable(
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
        fortran.tfcoil_variables, float, range=(0.001, 1000000.0)
    ),
    "c_tf_turn_max": InputVariable(
        fortran.tfcoil_variables, float, range=(1.0, 1000000.0)
    ),
    "crane_arm_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "crane_clrnc_h": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "crane_clrnc_v": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "croco_thick": InputVariable(
        data_structure.rebco_variables, float, range=(0.001, 0.1)
    ),
    "cryomag_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "cryomag_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryomag_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryostat_clrnc": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "cryostore_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "cryostore_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "cryostore_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "p_cryo_plant_electric_max_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.01, 200.0)
    ),
    "csi": InputVariable(data_structure.cost_variables, float, range=(1.0, 100.0)),
    "cturbb": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 1000.0)
    ),
    "dz_vv_lower": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dz_vv_upper": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dalu": InputVariable(fortran.constants, float, range=(2500.0, 30000.0)),
    "dcase": InputVariable(fortran.tfcoil_variables, float, range=(1000.0, 100000.0)),
    "dcdrv0": InputVariable(fortran.ife_variables, float, range=(0.0, 200.0)),
    "dcdrv1": InputVariable(fortran.ife_variables, float, range=(0.0, 200.0)),
    "dcdrv2": InputVariable(fortran.ife_variables, float, range=(0.0, 200.0)),
    "dcondins": InputVariable(fortran.tfcoil_variables, float, range=(500.0, 10000.0)),
    "dcopper": InputVariable(fortran.constants, float, range=(8000.0, 10000.0)),
    "declblkt": InputVariable(fortran.fwbs_variables, float, range=(0.01, 0.2)),
    "declfw": InputVariable(fortran.fwbs_variables, float, range=(0.01, 0.2)),
    "declshld": InputVariable(fortran.fwbs_variables, float, range=(0.01, 0.2)),
    "decomf": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "denstl": InputVariable(fortran.fwbs_variables, float, range=(5000.0, 10000.0)),
    "dia_tf_turn_coolant_channel": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 0.1)
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
    "dr_blkt_inboard": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_blkt_outboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_bore": InputVariable(fortran.build_variables, float, range=(0.0, 50.0)),
    "dr_cryostat": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_cs": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_cs_tf_gap": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_fw_plasma_gap_inboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_fw_plasma_gap_outboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_fw_wall": InputVariable(fortran.fwbs_variables, float, range=(0.0005, 0.1)),
    "dr_pf_cryostat": InputVariable(fortran.fwbs_variables, float, range=(0.1, 5.0)),
    "dr_shld_blkt_gap": InputVariable(fortran.build_variables, float, range=(0.0, 5.0)),
    "dr_shld_inboard": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_shld_outboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_thermal_inboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_thermal_outboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_shld_vv_gap_inboard": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dr_tf_inboard": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_tf_shld_gap": InputVariable(fortran.build_variables, float, range=(0.0, 5.0)),
    "dr_tf_wp_with_insulation": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 10.0)
    ),
    "dr_vv_inboard": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dr_vv_outboard": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "drtop": InputVariable(fortran.tfcoil_variables, float, range=(-1.5, 1.5)),
    "drveff": InputVariable(fortran.ife_variables, float, range=(0.01, 1.0)),
    "dtlife": InputVariable(data_structure.cost_variables, float, range=(0.0, 15.0)),
    "dtstor": InputVariable(data_structure.pulse_variables, float, range=(50.0, 500.0)),
    "dx_fw_module": InputVariable(fortran.fwbs_variables, float, range=(0.0005, 0.1)),
    "dz_tf_cryostat": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 20.0)
    ),
    "dztop": InputVariable(fortran.tfcoil_variables, float, range=(-0.5, 0.5)),
    "edrive": InputVariable(
        fortran.ife_variables, float, range=(100000.0, 5000000000.0)
    ),
    "eff_tf_cryo": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "ejima_coeff": InputVariable(fortran.physics_variables, float, range=(0.1, 1.0)),
    "elecdist_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "elecdist_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecdist_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecload_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "elecload_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecload_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecstore_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "elecstore_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "elecstore_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "f_p_blkt_multiplication": InputVariable(
        fortran.fwbs_variables, float, range=(1.0, 2.0)
    ),
    "esbldgm3": InputVariable(
        fortran.buildings_variables, float, range=(1000.0, 1000000.0)
    ),
    "eta_ecrh_injector_wall_plug": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "eta_coolant_pump_electric": InputVariable(
        fortran.fwbs_variables, float, range=(0.1, 1.0)
    ),
    "etaiso": InputVariable(fortran.fwbs_variables, float, range=(0.1, 1.0)),
    "eta_lowhyb_injector_wall_plug": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "etali": InputVariable(fortran.ife_variables, float, range=(0.0, 1.0)),
    "eta_beam_injector_wall_plug": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "etapsu": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 1.0)),
    "etapump": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "etatf": InputVariable(fortran.heat_transport_variables, float, range=(0.0, 1.0)),
    "eta_turbine": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1.0)
    ),
    "eyoung_al": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "eyoung_cond_axial": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 10000000000000.0)
    ),
    "eyoung_cond_trans": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 10000000000000.0)
    ),
    "eyoung_ins": InputVariable(
        fortran.tfcoil_variables, float, range=(100000000.0, 10000000000000.0)
    ),
    "eyoung_res_tf_buck": InputVariable(
        fortran.tfcoil_variables, float, range=(1e-10, 1000000000000.0)
    ),
    "eyoung_steel": InputVariable(
        fortran.tfcoil_variables, float, range=(100000000.0, 10000000000000.0)
    ),
    "f_a_tf_cool_outboard": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_alpha_energy_confinement_min": InputVariable(
        fortran.constraint_variables, float, range=(1.0, 100.0)
    ),
    "f_asym": InputVariable(fortran.stellarator_variables, float, range=(0.9, 2.0)),
    "f_avspace": InputVariable(fortran.build_variables, float, range=(0.001, 10.0)),
    "f_coppera_m2": InputVariable(
        data_structure.rebco_variables, float, range=(0.001, 10.0)
    ),
    "f_copperaoh_m2": InputVariable(
        data_structure.rebco_variables, float, range=(0.001, 1.0)
    ),
    "f_crypmw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "f_fw_peak": InputVariable(fortran.fwbs_variables, float, range=(1.0, 100.0)),
    "f_fw_rad_max": InputVariable(fortran.constraint_variables, float, range=(0.1, 10)),
    "f_nd_alpha_electron": InputVariable(
        fortran.physics_variables, float, range=(1e-12, 1.0)
    ),
    "f_nd_beam_electron": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_nd_protium_electrons": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "f_neut_shield": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "f_nuc_pow_bz_struct": InputVariable(
        fortran.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "f_r_cp": InputVariable(fortran.build_variables, float, range=(1.0, 100.0)),
    "f_rad": InputVariable(fortran.stellarator_variables, float, range=(0.0, 1.0)),
    "f_sync_reflect": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "f_t_turn_tf": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "f_tritium": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "f_beam_tritium": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "f_vforce_inboard": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "f_w": InputVariable(fortran.stellarator_variables, float, range=(0.1, 1.0)),
    "f_st_coil_aspect": InputVariable(fortran.stellarator_variables, float, range=(0.1, 10.0)),
    "f_z_cryostat": InputVariable(fortran.build_variables, float, range=(2.0, 10.0)),
    "falpha_energy_confinement": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1.0)
    ),
    "fauxbop": InputVariable(fortran.ife_variables, float, range=(0.0, 1.0)),
    "fauxmn": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "favail": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "fbeta_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fbeta_min": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fbeta_poloidal": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fbeta_poloidal_eps": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fblbe": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblbreed": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebmi": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebmo": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebpi": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblhebpo": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblli": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblli2o": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fbllipb": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblss": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fblvd": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fb_cs_limit_max": InputVariable(
        fortran.pfcoil_variables, float, range=(0.01, 1.0)
    ),
    "fbreed": InputVariable(fortran.ife_variables, float, range=(0.0, 0.999)),
    "fburn": InputVariable(fortran.ife_variables, float, range=(0.01, 1.0)),
    "fc_building_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "fc_building_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "fcap0": InputVariable(data_structure.cost_variables, float, range=(1.0, 1.5)),
    "fcap0cp": InputVariable(data_structure.cost_variables, float, range=(1.0, 1.5)),
    "fcdfuel": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "f_j_cs_start_pulse_end_flat_top": InputVariable(
        fortran.pfcoil_variables, float, range=(0.0, 1.0)
    ),
    "fcontng": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "fcoolcp": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "fcqt": InputVariable(fortran.constraint_variables, float, range=(0.001, 1.0)),
    "fcr0": InputVariable(data_structure.cost_variables, float, range=(0.0, 1.0)),
    "fcspc": InputVariable(fortran.build_variables, float, range=(0.0, 1.0)),
    "fcuohsu": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 1.0)),
    "fcupfsu": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 1.0)),
    "fcutfsu": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "fdene": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "f_ster_div_single": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fdiva": InputVariable(data_structure.divertor_variables, float, range=(0.1, 2.0)),
    "fdivwet": InputVariable(fortran.stellarator_variables, float, range=(0.01, 1.0)),
    "fdtmp": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fecrh_ignition": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "feffcd": InputVariable(fortran.current_drive_variables, float, range=(0.0, 20.0)),
    "fflutf": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fp_fusion_total_max_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "feta_cd_norm_hcd_primary_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "f_a_fw_hcd": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fpflux_div_heat_load_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fhole": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fhts": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 1.0)),
    "fiooic": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fipir": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fjohc": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fjohc0": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fjprot": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fkind": InputVariable(data_structure.cost_variables, float, range=(0.5, 1.0)),
    "fl_h_threshold": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1000000.0)
    ),
    "flirad": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "flpitch": InputVariable(
        fortran.stellarator_variables, float, range=(0.0001, 0.01)
    ),
    "f_div_flux_expansion": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "fmaxvvstress": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1.0)
    ),
    "fmgdmw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "fmva": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fnbshinef": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fncycle": InputVariable(fortran.constraint_variables, float, range=(1e-08, 1.0)),
    "fndt": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "fne0": InputVariable(fortran.physics_variables, float, range=(0.001, 1.0)),
    "fniterpump": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "foh_stress": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1.0)
    ),
    "f_p_beam_orbit_loss": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 0.999)
    ),
    "fp_plasma_separatrix_min_mw": InputVariable(
        fortran.physics_variables, float, range=(0.001, 1.0)
    ),
    "fb_tf_inboard_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fp_hcd_injected_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fp_plant_electric_net_required_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fpoloidalpower": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1.0)
    ),
    "fradius_beam_tangency": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fpsep": InputVariable(fortran.constraint_variables, float, range=(0.001, 1.0)),
    "fpsepbqar": InputVariable(fortran.constraint_variables, float, range=(0.001, 1.0)),
    "fpsepr": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fptemp": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fptfnuc": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fpumpblkt": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "fpumpdiv": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "fpumpfw": InputVariable(fortran.heat_transport_variables, float, range=(0.0, 0.2)),
    "fpumpshld": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 0.2)
    ),
    "fq": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fqval": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fr_conducting_wall": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fracture_toughness": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.1, 100000000.0)
    ),
    "fradpwr": InputVariable(fortran.constraint_variables, float, range=(0.0, 1.0)),
    "fpflux_fw_rad_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1.0)
    ),
    "f_radius_beam_tangency_rmajor": InputVariable(
        fortran.current_drive_variables, float, range=(0.5, 2.0)
    ),
    "freinke": InputVariable(fortran.constraint_variables, float, range=(0.001, 1.0)),
    "frhocp": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 5.0)),
    "frholeg": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 5.0)),
    "frminor": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "frrmax": InputVariable(fortran.ife_variables, float, range=(1e-06, 1.0)),
    "fseppc": InputVariable(
        fortran.build_variables, float, range=(1000000.0, 1000000000.0)
    ),
    "fstr_wp": InputVariable(fortran.constraint_variables, float, range=(1e-09, 10.0)),
    "fstrcase": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fstrcond": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "ft_burn_min": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "ft_current_ramp_up": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "ftbr": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "ft_cycle_min": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "ftmargoh": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "ftmargtf": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "ftoroidalgap": InputVariable(fortran.tfcoil_variables, float, range=(0.001, 10.0)),
    "ftemp_fw_max": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "fvdump": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fvoldw": InputVariable(fortran.fwbs_variables, float, range=(0.0, 10.0)),
    "fvolsi": InputVariable(fortran.fwbs_variables, float, range=(0.0, 10.0)),
    "fvolso": InputVariable(fortran.fwbs_variables, float, range=(0.0, 10.0)),
    "fvs_plasma_total_required": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
    ),
    "f_c_plasma_non_inductive": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "fvs_cs_pf_total_ramp": InputVariable(
        fortran.pfcoil_variables, float, range=(0.001, 10.0)
    ),
    "fvvhe": InputVariable(fortran.constraint_variables, float, range=(0.001, 10.0)),
    "fw_armour_thickness": InputVariable(
        fortran.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "fw_th_conductivity": InputVariable(
        fortran.fwbs_variables, float, range=(1.0, 100.0)
    ),
    "fpflux_fw_neutron_max_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 10.0)
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
    "fwclfr": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "fwdr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "fwdzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "fwdzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "fzactual": InputVariable(data_structure.reinke_variables, float, range=(0.0, 1.0)),
    "fzeffmax": InputVariable(fortran.constraint_variables, float, range=(0.001, 1.0)),
    "eta_cd_norm_ecrh": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1.0)
    ),
    "gamma_he": InputVariable(
        data_structure.primary_pumping_variables, float, range=(1.0, 2.0)
    ),
    "eta_cd_norm_hcd_primary_max": InputVariable(
        fortran.constraint_variables, float, range=(0.01, 10.0)
    ),
    "gapomin": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "gas_buildings_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "gas_buildings_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "gas_buildings_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ground_clrnc": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "n_ecrh_harmonic": InputVariable(
        fortran.current_drive_variables, float, range=(1.0, 10.0)
    ),
    "hastelloy_thickness": InputVariable(
        data_structure.rebco_variables, float, range=(1e-08, 0.001)
    ),
    "hccl": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "hcd_building_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "hcd_building_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hcd_building_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hcwt": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "heat_sink_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "heat_sink_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "heat_sink_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hfact": InputVariable(fortran.physics_variables, float, range=(0.01, 10.0)),
    "pflux_div_heat_load_mw": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "pflux_div_heat_load_max_mw": InputVariable(
        data_structure.divertor_variables, float, range=(0.1, 20.0)
    ),
    "hot_sepdist": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "hotcell_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "p_blkt_coolant_pump_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "p_div_coolant_pump_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "p_fw_coolant_pump_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "htpmw_ife": InputVariable(fortran.ife_variables, float, range=(0.0, 1000.0)),
    "p_shld_coolant_pump_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "hw_storage_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "hw_storage_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "hw_storage_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_smelter_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "ilw_smelter_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_smelter_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_storage_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "ilw_storage_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ilw_storage_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "ind_plasma_internal_norm": InputVariable(
        fortran.physics_variables, float, range=(0.0, 10.0)
    ),
    "initialpressure": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 10000.0)
    ),
    "temp_blkt_coolant_in": InputVariable(
        fortran.fwbs_variables, float, range=(200.0, 600.0)
    ),
    "inlet_temp_liq": InputVariable(
        fortran.fwbs_variables, float, range=(508.0, 1500.0)
    ),
    "iotabar": InputVariable(fortran.stellarator_variables, float, range=(0.1, 10.0)),
    "j_tf_bus": InputVariable(
        fortran.tfcoil_variables, float, range=(10000.0, 100000000.0)
    ),
    "kappa": InputVariable(fortran.physics_variables, float, range=(0.99, 5.0)),
    "kappa95": InputVariable(fortran.physics_variables, float, range=(0.99, 5.0)),
    "layer_ins": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 0.1)),
    "ld_ratio_cst": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 5.0)),
    "len_fw_channel": InputVariable(
        fortran.fwbs_variables, float, range=(0.001, 1000.0)
    ),
    "len_tf_bus": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 1000.0)),
    "lhat": InputVariable(data_structure.reinke_variables, float, range=(1.0, 15.0)),
    "f_blkt_li6_enrichment": InputVariable(
        fortran.fwbs_variables, float, range=(7.4, 100.0)
    ),
    "life_dpa": InputVariable(
        data_structure.cost_variables, float, range=(10.0, 100.0)
    ),
    "llw_storage_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "llw_storage_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "llw_storage_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "m_s_limit": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "magnet_pulse_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "magnet_pulse_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_pulse_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_trains_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "magnet_trains_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "magnet_trains_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maint_cont_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "maint_cont_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maint_cont_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "maintenance_fwbs": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "maintenance_gen": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "max_gyrotron_frequency": InputVariable(
        fortran.stellarator_variables, float, range=(1000000000.0, 100000000000000.0)
    ),
    "max_vv_stress": InputVariable(
        fortran.tfcoil_variables, float, range=(100000.0, 500000000.0)
    ),
    "maxpoloidalpower": InputVariable(
        fortran.pf_power_variables, float, range=(0.0, 2000.0)
    ),
    "pflux_fw_rad_max": InputVariable(
        fortran.constraint_variables, float, range=(0.1, 10.0)
    ),
    "mbvfac": InputVariable(fortran.buildings_variables, float, range=(0.9, 3.0)),
    "mcdriv": InputVariable(fortran.ife_variables, float, range=(0.1, 10.0)),
    "mvalim": InputVariable(fortran.constraint_variables, float, range=(0.0, 1000.0)),
    "n_cycle_min": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.0, 100000000.0)
    ),
    "n_tf_coils": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 100.0)),
    "n_tf_coil_turns": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "nbi_sys_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "nbi_sys_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "dx_beam_shield": InputVariable(
        fortran.current_drive_variables, float, range=(0.01, 0.5)
    ),
    "f_p_beam_shine_through_max": InputVariable(
        fortran.constraint_variables, float, range=(1e-20, 0.1)
    ),
    "neped": InputVariable(fortran.physics_variables, float, range=(0.0, 1e21)),
    "nesep": InputVariable(fortran.physics_variables, float, range=(0.0, 1e21)),
    "nflutfmax": InputVariable(fortran.constraint_variables, float, range=(1e20, 1e24)),
    "oacdcp": InputVariable(
        fortran.tfcoil_variables, float, range=(10000.0, 1000000000.0)
    ),
    "f_a_cs_steel": InputVariable(
        fortran.pfcoil_variables, float, range=(0.001, 0.999)
    ),
    "f_z_cs_tf_internal": InputVariable(
        fortran.pfcoil_variables, float, range=(0.0, 2.0)
    ),
    "outgasfactor": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "outgasindex": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "temp_blkt_coolant_out": InputVariable(
        fortran.fwbs_variables, float, range=(450.0, 900.0)
    ),
    "outlet_temp_liq": InputVariable(
        fortran.fwbs_variables, float, range=(508.0, 1500.0)
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
    "pbase": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-08, 0.001)
    ),
    "p_plasma_separatrix_min_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.1, 1000.0)
    ),
    "pdrive": InputVariable(
        fortran.ife_variables, float, range=(1000000.0, 200000000.0)
    ),
    "pfbldgm3": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "rho_pf_coil": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 0.0001)),
    "pfusife": InputVariable(fortran.ife_variables, float, range=(0.0, 10000.0)),
    "p_hcd_primary_extra_heat_mw": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "p_hcd_secondary_extra_heat_mw": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "pibv": InputVariable(fortran.buildings_variables, float, range=(1000.0, 100000.0)),
    "pifecr": InputVariable(fortran.ife_variables, float, range=(0.0, 100.0)),
    "p_hcd_injected_max": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "p_hcd_secondary_injected_mw": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 1000.0)
    ),
    "plasma_res_factor": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "plasma_square": InputVariable(fortran.physics_variables, float, range=(-5.0, 5.0)),
    "plleni": InputVariable(fortran.build_variables, float, range=(0.1, 10.0)),
    "plleno": InputVariable(fortran.build_variables, float, range=(0.1, 10.0)),
    "plsepi": InputVariable(fortran.build_variables, float, range=(0.1, 10.0)),
    "plsepo": InputVariable(fortran.build_variables, float, range=(0.1, 10.0)),
    "p_plant_electric_net_required_mw": InputVariable(
        fortran.constraint_variables, float, range=(1.0, 10000.0)
    ),
    "pnuc_fw_ratio_dcll": InputVariable(
        fortran.fwbs_variables, float, range=(0.0, 1.0)
    ),
    "poisson_al": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "poisson_copper": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "poisson_steel": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "p_fusion_total_max_mw": InputVariable(
        fortran.constraint_variables, float, range=(1.0, 10000.0)
    ),
    "prdiv": InputVariable(data_structure.vacuum_variables, float, range=(0.0, 10.0)),
    "pres_fw_coolant": InputVariable(
        fortran.fwbs_variables, float, range=(100000.0, 100000000.0)
    ),
    "prn1": InputVariable(data_structure.divertor_variables, float, range=(0.0, 1.0)),
    "psepbqarmax": InputVariable(
        fortran.constraint_variables, float, range=(1.0, 50.0)
    ),
    "pseprmax": InputVariable(fortran.constraint_variables, float, range=(1.0, 60.0)),
    "ptargf": InputVariable(fortran.ife_variables, float, range=(0.1, 100.0)),
    "ptempalw": InputVariable(fortran.tfcoil_variables, float, range=(4.0, 573.15)),
    "ptfnucmax": InputVariable(fortran.constraint_variables, float, range=(1e-06, 1.0)),
    "pulsetimings": InputVariable(
        data_structure.times_variables, float, range=(0.0, 1.0)
    ),
    "pumpareafraction": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1.0)
    ),
    "pumpspeedfactor": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1.0)
    ),
    "pumpspeedmax": InputVariable(
        data_structure.vacuum_variables, float, range=(1e-06, 1000.0)
    ),
    "pumptp": InputVariable(data_structure.vacuum_variables, float, range=(0.0, 1e30)),
    "pflux_plant_floor_electric": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 1000.0)
    ),
    "q0": InputVariable(fortran.physics_variables, float, range=(0.01, 20.0)),
    "q95": InputVariable(fortran.physics_variables, float, range=(1.0, 50.0)),
    "q95_fixed": InputVariable(fortran.constraint_variables, float, range=(1.0, 50.0)),
    "qnty_sfty_fac": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "qnuc": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1000000.0)),
    "r_cp_top": InputVariable(fortran.build_variables, float, range=(0.001, 10.0)),
    "rad_fraction_sol": InputVariable(
        fortran.physics_variables, float, range=(0.0, 1.0)
    ),
    "radius_fw_channel": InputVariable(
        fortran.fwbs_variables, float, range=(0.001, 0.5)
    ),
    "rat": InputVariable(data_structure.vacuum_variables, float, range=(1e-10, 1e-06)),
    "rbrt": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "rbvfac": InputVariable(fortran.buildings_variables, float, range=(0.9, 3.0)),
    "rbwt": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "rcool": InputVariable(fortran.tfcoil_variables, float, range=(1e-06, 1.0)),
    "reactor_clrnc": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "reactor_fndtn_thk": InputVariable(
        fortran.buildings_variables, float, range=(0.25, 25.0)
    ),
    "reactor_hall_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "reactor_hall_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "reactor_hall_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "reactor_roof_thk": InputVariable(
        fortran.buildings_variables, float, range=(0.25, 25.0)
    ),
    "reactor_wall_thk": InputVariable(
        fortran.buildings_variables, float, range=(0.25, 25.0)
    ),
    "rebco_thickness": InputVariable(
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
    "rho_tf_bus": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1e-05)),
    "rho_tf_joints": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 0.01)),
    "rhopedn": InputVariable(fortran.physics_variables, float, range=(0.01, 1.0)),
    "rhopedt": InputVariable(fortran.physics_variables, float, range=(0.01, 1.0)),
    "rhopfbus": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 1e-05)),
    "rinboard": InputVariable(fortran.build_variables, float, range=(0.1, 10.0)),
    "ripmax": InputVariable(fortran.tfcoil_variables, float, range=(0.1, 100.0)),
    "rmajor": InputVariable(fortran.physics_variables, float, range=(0.1, 50.0)),
    "robotics_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "robotics_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "robotics_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "roughness_fw_channel": InputVariable(
        fortran.fwbs_variables, float, range=(0.0, 0.01)
    ),
    "routr": InputVariable(fortran.pfcoil_variables, float, range=(-3.0, 3.0)),
    "row": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "rpf1": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 3.0)),
    "rpf2": InputVariable(fortran.pfcoil_variables, float, range=(-3.0, 3.0)),
    "rrin": InputVariable(fortran.ife_variables, float, range=(0.1, 50.0)),
    "rrmax": InputVariable(fortran.ife_variables, float, range=(1.0, 50.0)),
    "rrr_tf_cu": InputVariable(fortran.tfcoil_variables, float, range=(1.0, 1000.0)),
    "rxcl": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "sec_buildings_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "sec_buildings_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "sec_buildings_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
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
    "shdr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "shdzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "shdzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "shear": InputVariable(fortran.stellarator_variables, float, range=(0.1, 10.0)),
    "dz_shld_lower": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dz_shld_upper": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "shmf": InputVariable(fortran.buildings_variables, float, range=(0.0, 1.0)),
    "shov": InputVariable(
        fortran.buildings_variables, float, range=(1000.0, 1000000.0)
    ),
    "sig_tf_case_max": InputVariable(
        fortran.tfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "sig_tf_wp_max": InputVariable(
        fortran.tfcoil_variables, float, range=(1000000.0, 100000000000.0)
    ),
    "sigallpc": InputVariable(
        fortran.build_variables, float, range=(0.0, 1000000000.0)
    ),
    "sigpfcalw": InputVariable(fortran.pfcoil_variables, float, range=(1.0, 1000.0)),
    "sigpfcf": InputVariable(fortran.pfcoil_variables, float, range=(0.1, 1.0)),
    "sombdr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "somtdr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "staff_buildings_area": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "staff_buildings_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "startupratio": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 10.0)
    ),
    "stcl": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "str_cs_con_res": InputVariable(
        fortran.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_pf_con_res": InputVariable(
        fortran.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_tf_con_res": InputVariable(
        fortran.tfcoil_variables, float, range=(-0.02, 0.02)
    ),
    "str_wp_max": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 0.3)),
    "t_between_pulse": InputVariable(
        data_structure.times_variables, float, range=(0.0, 100000000.0)
    ),
    "t_burn": InputVariable(
        data_structure.times_variables, float, range=(0.0, 100000000.0)
    ),
    "t_burn_min": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 1000000.0)
    ),
    "t_cable_tf": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 0.1)),
    "t_crack_radial": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1e-05, 1.0)
    ),
    "t_crack_vertical": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(1e-05, 1.0)
    ),
    "t_crit_nbti": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 15.0)),
    "t_current_ramp_up": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_fusion_ramp": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_in_bb": InputVariable(
        data_structure.primary_pumping_variables, float, range=(200.0, 1000.0)
    ),
    "t_tf_quench_detection": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 100.0)
    ),
    "t_out_bb": InputVariable(
        data_structure.primary_pumping_variables, float, range=(200.0, 1000.0)
    ),
    "t_precharge": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_ramp_down": InputVariable(
        data_structure.times_variables, float, range=(0.0, 10000.0)
    ),
    "t_structural_radial": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.001, 1.0)
    ),
    "t_structural_vertical": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.001, 1.0)
    ),
    "t_turn_tf": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 0.1)),
    "t_turn_tf_max": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "tape_thickness": InputVariable(
        data_structure.rebco_variables, float, range=(0.0, 0.1)
    ),
    "tape_width": InputVariable(
        data_structure.rebco_variables, float, range=(0.0, 0.1)
    ),
    "tauee_in": InputVariable(fortran.physics_variables, float, range=(0.0, 100.0)),
    "taumax": InputVariable(fortran.physics_variables, float, range=(0.1, 100.0)),
    "tauratio": InputVariable(fortran.physics_variables, float, range=(0.1, 100.0)),
    "n_beam_decay_lengths_core_required": InputVariable(
        fortran.current_drive_variables, float, range=(0.0, 10.0)
    ),
    "tbeta": InputVariable(fortran.physics_variables, float, range=(0.0, 4.0)),
    "tbktrepl": InputVariable(data_structure.cost_variables, float, range=(0.01, 2.0)),
    "tbrmin": InputVariable(fortran.constraint_variables, float, range=(0.001, 2.0)),
    "tcomrepl": InputVariable(data_structure.cost_variables, float, range=(0.01, 2.0)),
    "tcoolin": InputVariable(fortran.tfcoil_variables, float, range=(4.0, 373.15)),
    "tcritsc": InputVariable(fortran.tfcoil_variables, float, range=(1.0, 300.0)),
    "t_cycle_min": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 2000000.0)
    ),
    "tdiv": InputVariable(data_structure.divertor_variables, float, range=(0.1, 100.0)),
    "tdivrepl": InputVariable(data_structure.cost_variables, float, range=(0.01, 2.0)),
    "tdmptf": InputVariable(fortran.tfcoil_variables, float, range=(0.1, 100.0)),
    "te": InputVariable(fortran.physics_variables, float, range=(1.0, 200.0)),
    "te0_ecrh_achievable": InputVariable(
        fortran.stellarator_variables, float, range=(1.0, 1000.0)
    ),
    "temp_cp_average": InputVariable(
        fortran.tfcoil_variables, float, range=(4.0, 573.15)
    ),
    "temp_fw_coolant_in": InputVariable(
        fortran.fwbs_variables, float, range=(300.0, 1500.0)
    ),
    "temp_fw_coolant_out": InputVariable(
        fortran.fwbs_variables, float, range=(300.0, 1500.0)
    ),
    "temp_fw_max": InputVariable(fortran.fwbs_variables, float, range=(500.0, 2000.0)),
    "teped": InputVariable(fortran.physics_variables, float, range=(0.0, 20.0)),
    "tesep": InputVariable(fortran.physics_variables, float, range=(0.0, 20.0)),
    "tfcbv": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "dx_tf_wp_insertion_gap": InputVariable(
        fortran.tfcoil_variables, float, range=(1e-10, 0.1)
    ),
    "tfootfi": InputVariable(fortran.build_variables, float, range=(0.2, 5.0)),
    "tftmp": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 293.0)),
    "tgain": InputVariable(fortran.ife_variables, float, range=(1.0, 500.0)),
    "th_joint_contact": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "theta1_coil": InputVariable(fortran.tfcoil_variables, float, range=(0.1, 60.0)),
    "theta1_vv": InputVariable(fortran.tfcoil_variables, float, range=(0.1, 60.0)),
    "dx_tf_turn_insulation": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "dr_tf_nose_case": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1.0)),
    "dz_shld_thermal": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dx_tf_turn_steel": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "ti": InputVariable(fortran.physics_variables, float, range=(5.0, 50.0)),
    "dx_tf_wp_insulation": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 0.1)
    ),
    "tlife": InputVariable(data_structure.cost_variables, float, range=(1.0, 100.0)),
    "tmain": InputVariable(data_structure.cost_variables, float, range=(0.0, 100.0)),
    "tmargmin": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 20.0)),
    "tmargmin_cs": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 20.0)),
    "tmargmin_tf": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 20.0)),
    "tmax_croco": InputVariable(fortran.tfcoil_variables, float, range=(4.0, 1000.0)),
    "tmaxpro": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 1000.0)),
    "temp_tf_cryo": InputVariable(fortran.tfcoil_variables, float, range=(0.01, 293.0)),
    "tn": InputVariable(data_structure.vacuum_variables, float, range=(1.0, 1000.0)),
    "i_t_current_ramp_up": InputVariable(
        data_structure.times_variables, int, choices=[0, 1]
    ),
    "transp_clrnc": InputVariable(
        fortran.buildings_variables, float, range=(0.0, 10.0)
    ),
    "tratio": InputVariable(fortran.physics_variables, float, range=(0.0, 2.0)),
    "trcl": InputVariable(fortran.buildings_variables, float, range=(0.0, 10.0)),
    "triang": InputVariable(fortran.physics_variables, float, range=(-1.0, 1.0)),
    "triang95": InputVariable(fortran.physics_variables, float, range=(0.0, 1.0)),
    "p_tritium_plant_electric_mw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "triv": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "turbine_hall_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "turbine_hall_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "turbine_hall_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "tw_storage_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "tw_storage_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "tw_storage_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "u_unplanned_cp": InputVariable(
        data_structure.cost_variables, float, range=(0.0, 1.0)
    ),
    "ucblbe": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucblbreed": InputVariable(
        data_structure.cost_variables, float, range=(1.0, 1000.0)
    ),
    "ucblli": InputVariable(
        data_structure.cost_variables, float, range=(10.0, 10000.0)
    ),
    "ucblli2o": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 10000.0)
    ),
    "ucbllipb": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 10000.0)
    ),
    "ucblss": InputVariable(data_structure.cost_variables, float, range=(10.0, 1000.0)),
    "ucblvd": InputVariable(
        data_structure.cost_variables, float, range=(100.0, 1000.0)
    ),
    "ucbus": InputVariable(data_structure.cost_variables, float, range=(0.01, 10.0)),
    "uccarb": InputVariable(fortran.ife_variables, float, range=(10.0, 1000.0)),
    "uccase": InputVariable(data_structure.cost_variables, float, range=(1.0, 1000.0)),
    "ucconc": InputVariable(fortran.ife_variables, float, range=(0.1, 1000.0)),
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
    "ucflib": InputVariable(fortran.ife_variables, float, range=(10.0, 1000.0)),
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
    "uctarg": InputVariable(fortran.ife_variables, float, range=(0.1, 1000.0)),
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
    "v1dr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v1dzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v1dzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v2dr": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v2dzl": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v2dzu": InputVariable(fortran.ife_variables, float, range=(0.0, 10.0)),
    "v3dr": InputVariable(fortran.ife_variables, float, range=(0.0, 50.0)),
    "v3dzl": InputVariable(fortran.ife_variables, float, range=(0.0, 30.0)),
    "v3dzu": InputVariable(fortran.ife_variables, float, range=(0.0, 30.0)),
    "vachtmw": InputVariable(
        fortran.heat_transport_variables, float, range=(0.0, 100.0)
    ),
    "vcool": InputVariable(fortran.tfcoil_variables, float, range=(0.001, 100.0)),
    "vdalw": InputVariable(fortran.tfcoil_variables, float, range=(0.0, 100.0)),
    "vfblkt": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "vfcblkt": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "f_a_cs_void": InputVariable(fortran.pfcoil_variables, float, range=(0.0, 1.0)),
    "vfpblkt": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "vfshld": InputVariable(fortran.fwbs_variables, float, range=(0.0, 1.0)),
    "f_a_tf_turn_cable_space_extra_void": InputVariable(
        fortran.tfcoil_variables, float, range=(0.0, 1.0)
    ),
    "dz_shld_vv_gap": InputVariable(fortran.build_variables, float, range=(0.0, 10.0)),
    "dz_xpoint_divertor": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "dz_fw_plasma_gap": InputVariable(
        fortran.build_variables, float, range=(0.0, 10.0)
    ),
    "vvhealw": InputVariable(fortran.constraint_variables, float, range=(0.01, 10.0)),
    "pflux_fw_neutron_max_mw": InputVariable(
        fortran.constraint_variables, float, range=(0.001, 50.0)
    ),
    "walker_coefficient": InputVariable(
        data_structure.cs_fatigue_variables, float, range=(0.1, 10.0)
    ),
    "wallpf": InputVariable(fortran.fwbs_variables, float, range=(1.0, 2.0)),
    "warm_shop_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "warm_shop_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "warm_shop_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "water_buildings_h": InputVariable(
        fortran.buildings_variables, float, range=(1.0, 100.0)
    ),
    "water_buildings_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "water_buildings_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "watertemp": InputVariable(
        data_structure.water_usage_variables, float, range=(0.0, 25.0)
    ),
    "wgt": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "wgt2": InputVariable(
        fortran.buildings_variables, float, range=(10000.0, 1000000.0)
    ),
    "windspeed": InputVariable(
        data_structure.water_usage_variables, float, range=(0.0, 10.0)
    ),
    "workshop_h": InputVariable(fortran.buildings_variables, float, range=(1.0, 100.0)),
    "workshop_l": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "workshop_w": InputVariable(
        fortran.buildings_variables, float, range=(10.0, 1000.0)
    ),
    "wsvfac": InputVariable(fortran.buildings_variables, float, range=(0.9, 3.0)),
    "xi_ebw": InputVariable(fortran.current_drive_variables, float, range=(0.0, 1.0)),
    "xpertin": InputVariable(
        data_structure.divertor_variables, float, range=(0.0, 10.0)
    ),
    "zeffmax": InputVariable(fortran.constraint_variables, float, range=(1.0, 10.0)),
    "blktmodel": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "blkttype": InputVariable(fortran.fwbs_variables, int, choices=[1, 2, 3]),
    "breedmat": InputVariable(fortran.fwbs_variables, int, choices=[1, 2, 3]),
    "ccl0_ma": InputVariable(fortran.pfcoil_variables, float, array=True),
    "ccls_ma": InputVariable(fortran.pfcoil_variables, float, array=True),
    "cfind": InputVariable(data_structure.cost_variables, float, array=True),
    "i_blkt_coolant_type": InputVariable(fortran.fwbs_variables, int, choices=[1, 2]),
    "coppera_m2_max": InputVariable(
        data_structure.rebco_variables, float, range=(1.0e6, 1.0e10)
    ),
    "cost_model": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "dwell_pump": InputVariable(
        data_structure.vacuum_variables, int, choices=[0, 1, 2]
    ),
    "i_fw_blkt_vv_shape": InputVariable(fortran.fwbs_variables, int, range=(1, 2)),
    "hcdportsize": InputVariable(fortran.fwbs_variables, int, range=(1, 2)),
    "i_blkt_liquid_breeder_type": InputVariable(
        fortran.fwbs_variables, int, choices=[0, 1]
    ),
    "i_beta_component": InputVariable(fortran.physics_variables, int, range=(0, 3)),
    "i_beta_fast_alpha": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "i_blanket_type": InputVariable(fortran.fwbs_variables, int, choices=[1, 3, 5]),
    "i_bldgs_size": InputVariable(fortran.buildings_variables, int, choices=[0, 1]),
    "i_bldgs_v": InputVariable(fortran.buildings_variables, int, choices=[0, 1]),
    "i_blkt_inboard": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "i_bootstrap_current": InputVariable(fortran.physics_variables, int, range=(0, 13)),
    "i_cp_joints": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1]),
    "i_cp_lifetime": InputVariable(data_structure.cost_variables, int, range=(0, 3)),
    "i_cs_precomp": InputVariable(fortran.build_variables, int, choices=[0, 1]),
    "i_cs_stress": InputVariable(fortran.pfcoil_variables, int, choices=[0, 1]),
    "i_density_limit": InputVariable(fortran.physics_variables, int, range=(1, 8)),
    "i_diamagnetic_current": InputVariable(
        fortran.physics_variables, int, choices=[0, 1, 2]
    ),
    "i_div_heat_load": InputVariable(
        data_structure.divertor_variables, int, choices=[0, 1, 2]
    ),
    "i_l_h_threshold": InputVariable(fortran.physics_variables, int, range=(1, 21)),
    "i_pf_current": InputVariable(fortran.pfcoil_variables, int, choices=[0, 1, 2]),
    "i_pfirsch_schluter_current": InputVariable(
        fortran.physics_variables, int, choices=[0, 1]
    ),
    "i_plasma_current": InputVariable(fortran.physics_variables, int, range=(1, 9)),
    "i_plasma_geometry": InputVariable(fortran.physics_variables, int, range=(0, 11)),
    "i_plasma_shape": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "i_plasma_wall_gap": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "i_pulsed_plant": InputVariable(
        data_structure.pulse_variables, int, choices=[0, 1]
    ),
    "i_q95_fixed": InputVariable(fortran.constraint_variables, int, choices=[0, 1]),
    "i_r_cp_top": InputVariable(fortran.build_variables, int, choices=[0, 1, 2]),
    "i_rad_loss": InputVariable(fortran.physics_variables, int, choices=[0, 1, 2]),
    "i_shield_mat": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "i_single_null": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "i_str_wp": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1]),
    "i_sup_pf_shape": InputVariable(fortran.pfcoil_variables, int, choices=[0, 1]),
    "i_tf_bucking": InputVariable(fortran.tfcoil_variables, int, range=(0, 3)),
    "i_tf_case_geom": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1]),
    "i_tf_cond_eyoung_axial": InputVariable(
        fortran.tfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_tf_cond_eyoung_trans": InputVariable(
        fortran.tfcoil_variables, int, choices=[0, 1]
    ),
    "i_tf_sc_mat": InputVariable(fortran.tfcoil_variables, int, range=(1, 9)),
    "i_tf_shape": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1, 2]),
    "i_tf_stress_model": InputVariable(
        fortran.tfcoil_variables, int, choices=[0, 1, 2]
    ),
    "i_tf_sup": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1, 2]),
    "i_tf_tresca": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1]),
    "i_tf_turns_integer": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1]),
    "i_tf_wp_geom": InputVariable(fortran.tfcoil_variables, int, choices=[0, 1, 2]),
    "iavail": InputVariable(data_structure.cost_variables, int, range=(0, 3)),
    "ibkt_life": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "i_blkt_dual_coolant": InputVariable(
        fortran.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "i_hcd_primary": InputVariable(fortran.current_drive_variables, int, range=(1, 13)),
    "i_hcd_secondary": InputVariable(
        fortran.current_drive_variables, int, range=(0, 13)
    ),
    "i_blkt_liquid_breeder_channel_type": InputVariable(
        fortran.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "ife": InputVariable(fortran.ife_variables, int, choices=[0, 1]),
    "ifedrv": InputVariable(fortran.ife_variables, int, range=(-1, 3)),
    "ifetyp": InputVariable(fortran.ife_variables, int, range=(0, 4)),
    "ifueltyp": InputVariable(data_structure.cost_variables, int, choices=[0, 1, 2]),
    "i_plasma_ignited": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "ims": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "inuclear": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "iohcl": InputVariable(fortran.build_variables, int, choices=[0, 1]),
    "ipedestal": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "i_pf_conductor": InputVariable(fortran.pfcoil_variables, int, choices=[0, 1]),
    "ipnet": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "ipowerflow": InputVariable(fortran.heat_transport_variables, int, choices=[0, 1]),
    "i_shld_primary_heat": InputVariable(
        fortran.heat_transport_variables, int, choices=[0, 1]
    ),
    "i_beta_norm_max": InputVariable(fortran.physics_variables, int, range=(0, 5)),
    "i_ind_plasma_internal_norm": InputVariable(
        fortran.physics_variables, int, range=(0, 2)
    ),
    "i_alphaj": InputVariable(fortran.physics_variables, int, range=(0, 1)),
    "i_fw_blkt_shared_coolant": InputVariable(
        fortran.fwbs_variables, int, choices=[0, 1, 2]
    ),
    "ireactor": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "irefprop": InputVariable(fortran.fwbs_variables, int, choices=[0, 1]),
    "i_hcd_calculations": InputVariable(
        fortran.current_drive_variables, int, choices=[0, 1]
    ),
    "i_pf_energy_storage_source": InputVariable(
        fortran.pf_power_variables, int, range=(1, 3)
    ),
    "istell": InputVariable(fortran.stellarator_variables, int, range=(0, 6)),
    "isthtr": InputVariable(fortran.stellarator_variables, int, range=(1, 3)),
    "istore": InputVariable(data_structure.pulse_variables, int, range=(1, 3)),
    "i_cs_superconductor": InputVariable(fortran.pfcoil_variables, int, range=(1, 9)),
    "i_pf_superconductor": InputVariable(fortran.pfcoil_variables, int, range=(1, 9)),
    "itart": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "itartpf": InputVariable(fortran.physics_variables, int, choices=[0, 1]),
    "itcycl": InputVariable(data_structure.pulse_variables, int, range=(1, 3)),
    "iwalld": InputVariable(fortran.physics_variables, int, range=(1, 2)),
    "lsa": InputVariable(data_structure.cost_variables, int, range=(1, 4)),
    "m_res": InputVariable(fortran.stellarator_variables, int, range=(1, 10)),
    "n_layer": InputVariable(fortran.tfcoil_variables, int, range=(1, 100)),
    "n_liq_recirc": InputVariable(fortran.fwbs_variables, int, range=(1, 50)),
    "n_pancake": InputVariable(fortran.tfcoil_variables, int, range=(1, 100)),
    "n_rad_per_layer": InputVariable(fortran.tfcoil_variables, int, range=(1, 500)),
    "n_res": InputVariable(fortran.stellarator_variables, int, range=(3, 6)),
    "n_tf_graded_layers": InputVariable(fortran.tfcoil_variables, int, range=(1, 20)),
    "n_tf_joints": InputVariable(fortran.tfcoil_variables, int, range=(1, 50)),
    "n_tf_joints_contact": InputVariable(fortran.tfcoil_variables, int, range=(1, 50)),
    "n_blkt_inboard_modules_poloidal": InputVariable(
        fortran.fwbs_variables, int, range=(1, 16)
    ),
    "n_blkt_outboard_modules_poloidal": InputVariable(
        fortran.fwbs_variables, int, range=(1, 16)
    ),
    "n_blkt_inboard_modules_toroidal": InputVariable(
        fortran.fwbs_variables, int, range=(8, 96)
    ),
    "n_blkt_outboard_modules_toroidal": InputVariable(
        fortran.fwbs_variables, int, range=(8, 96)
    ),
    "npdiv": InputVariable(fortran.fwbs_variables, int, range=(0, 4)),
    "nphcdin": InputVariable(fortran.fwbs_variables, int, range=(0, 4)),
    "nphcdout": InputVariable(fortran.fwbs_variables, int, range=(0, 4)),
    "ntype": InputVariable(data_structure.vacuum_variables, int, choices=[0, 1]),
    "num_rh_systems": InputVariable(data_structure.cost_variables, int, range=(1, 10)),
    "output_costs": InputVariable(data_structure.cost_variables, int, choices=[0, 1]),
    "i_coolant_pumping": InputVariable(fortran.fwbs_variables, int, range=(0, 3)),
    "reinke_mode": InputVariable(data_structure.reinke_variables, int, choices=[0, 1]),
    "scan_dim": InputVariable(fortran.scan_module, int, range=(1, 2)),
    "i_thermal_electric_conversion": InputVariable(
        fortran.fwbs_variables, int, range=(0, 4)
    ),
    "secondary_cycle_liq": InputVariable(fortran.fwbs_variables, int, range=(2, 4)),
    "supercond_cost_model": InputVariable(
        data_structure.cost_variables, int, choices=[0, 1]
    ),
    "i_tf_inside_cs": InputVariable(fortran.build_variables, int, choices=[0, 1]),
    "i_ecrh_wave_mode": InputVariable(
        fortran.current_drive_variables, int, choices=[0, 1]
    ),
    "i_confinement_time": InputVariable(
        fortran.physics_variables,
        int,
        choices=list(range(fortran.physics_variables.n_confinement_scalings)),
    ),
    "quench_model": InputVariable(
        fortran.tfcoil_variables, str, choices=["exponential", "linear"]
    ),
    "i_fw_coolant_type": InputVariable(
        fortran.fwbs_variables, str, choices=["helium", "water"]
    ),
    "vacuum_model": InputVariable(
        data_structure.vacuum_variables, str, choices=["old", "simple"]
    ),
    "dcond": InputVariable(fortran.tfcoil_variables, float, array=True),
    "c_pf_coil_turn_peak_input": InputVariable(
        fortran.pfcoil_variables, float, array=True
    ),
    "i_pf_location": InputVariable(fortran.pfcoil_variables, int, array=True),
    "n_pf_coils_in_group": InputVariable(fortran.pfcoil_variables, int, array=True),
    "nfxfh": InputVariable(fortran.pfcoil_variables, int, array=True),
    "n_pf_coil_groups": InputVariable(
        fortran.pfcoil_variables,
        int,
        range=(0, fortran.pfcoil_variables.n_pf_groups_max.item()),
    ),
    "rref": InputVariable(fortran.pfcoil_variables, float, array=True),
    "f_a_pf_coil_void": InputVariable(fortran.pfcoil_variables, float, array=True),
    "zref": InputVariable(fortran.pfcoil_variables, float, array=True),
    "uchts": InputVariable(data_structure.cost_variables, float, array=True),
    "ucoam": InputVariable(data_structure.cost_variables, float, array=True),
    "ucsc": InputVariable(data_structure.cost_variables, float, array=True),
    "ucturb": InputVariable(data_structure.cost_variables, float, array=True),
    "ucwst": InputVariable(data_structure.cost_variables, float, array=True),
    "blmatf": InputVariable(fortran.ife_variables, float, array=True),
    "chmatf": InputVariable(fortran.ife_variables, float, array=True),
    "etave": InputVariable(fortran.ife_variables, float, array=True),
    "fwmatf": InputVariable(fortran.ife_variables, float, array=True),
    "gainve": InputVariable(fortran.ife_variables, float, array=True),
    "shmatf": InputVariable(fortran.ife_variables, float, array=True),
    "v1matf": InputVariable(fortran.ife_variables, float, array=True),
    "v2matf": InputVariable(fortran.ife_variables, float, array=True),
    "v3matf": InputVariable(fortran.ife_variables, float, array=True),
    "isweep": InputVariable(
        fortran.scan_module, int, choices=range(fortran.scan_module.ipnscns.item() + 1)
    ),
    "nsweep": InputVariable(
        fortran.scan_module,
        int,
        choices=range(1, fortran.scan_module.ipnscnv.item() + 1),
    ),
    "isweep_2": InputVariable(
        fortran.scan_module, int, choices=range(fortran.scan_module.ipnscns.item() + 1)
    ),
    "nsweep_2": InputVariable(
        fortran.scan_module,
        int,
        choices=range(1, fortran.scan_module.ipnscnv.item() + 1),
    ),
    "sweep": InputVariable(fortran.scan_module, float, array=True),
    "sweep_2": InputVariable(fortran.scan_module, float, array=True),
    "impvardiv": InputVariable(
        data_structure.reinke_variables,
        int,
        choices=range(3, fortran.impurity_radiation_module.n_impurities.item() + 1),
    ),
    "j_pf_coil_wp_peak": InputVariable(fortran.pfcoil_variables, float, array=True),
    "ixc": InputVariable(
        None,
        int,
        range=(1, fortran.numerics.ipnvars.item()),
        additional_actions=_ixc_additional_actions,
        set_variable=False,
    ),
    "icc": InputVariable(
        None,
        int,
        range=(1, fortran.numerics.ipeqns.item()),
        additional_actions=_icc_additional_actions,
        set_variable=False,
    ),
}


def parse_input_file():
    input_file = f2py_compatible_to_string(fortran.global_variables.fileprefix)

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

    if config.type is str:
        setattr(config.module, name, string_to_f2py_compatible(current_value, value))
    else:
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
        error_msg = f"Fortran module '{config.module}' does not have an array '{name}'."
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
