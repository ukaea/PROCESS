"""Module containing superconducting TF coil routines
author: P J Knight, CCFE, Culham Science Centre
author: J Morris, CCFE, Culham Science Centre
author: S Kahn, CCFE, Culham Science Centre
N/A
This module contains routines for calculating the
parameters of a superconducting TF coil system for a
fusion power plant.
PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
"""

tf_fit_t: float = None
"""Dimensionless winding pack width"""


tf_fit_z: float = None
"""Dimensionless winding pack radial thickness"""


f_b_tf_inboard_peak_ripple_symmetric: float = None
"""Ratio of peak field with ripple to nominal axisymmetric peak field"""


c_tf_coil: float = None
"""Current in each TF coil"""


a_tf_wp_with_insulation: float = None
"""Total cross-sectional area of winding pack including
GW insulation and insertion gap [m2]
"""


a_tf_wp_no_insulation: float = None
"""Total cross-sectional area of winding pack without
ground insulation and insertion gap [m2]
"""


a_tf_coil_inboard_steel: float = None
"""Inboard coil steel coil cross-sectional area [m2]"""


a_tf_coil_inboard_insulation: float = None
"""Inboard coil insulation cross-section per coil [m2]"""


f_a_tf_coil_inboard_steel: float = None
"""Inboard coil steel fraction [-]"""


f_a_tf_coil_inboard_insulation: float = None
"""Inboard coil insulation fraction [-]"""


z_cp_top: float = None
"""Vertical distance from the midplane to the top of the tapered section [m]"""


r_tf_outboard_in: float = None
"""Radial position of plasma-facing edge of TF coil outboard leg [m]"""


r_tf_outboard_out: float = None
"""Radial position of outer edge of TF coil inboard leg [m]"""


r_tf_wp_inboard_inner: float = None
"""Radial position of inner edge and centre of winding pack [m]"""


r_tf_wp_inboard_outer: float = None
"""Radial position of outer edge and centre of winding pack [m]"""


r_tf_wp_inboard_centre: float = None
"""Radial position of centre and centre of winding pack [m]"""


dr_tf_wp_top: float = None
"""Conductor layer radial thickness at centercollumn top [m]
Ground insulation layer included, only defined for itart = 1
"""


vol_ins_cp: float = None
"""CP turn insulation volume [m3]"""


vol_gr_ins_cp: float = None
"""CP ground insulation volume [m3]"""


vol_case_cp: float = None
"""Volume of the CP outer casing cylinder"""


dx_tf_wp_toroidal_min: float = None
"""Minimal toroidal thickness of of winding pack [m]"""


dx_tf_wp_toroidal_average: float = None
"""Averaged toroidal thickness of of winding pack [m]"""


dx_tf_side_case_average: float = None
"""Average lateral casing thickness [m]"""


a_tf_plasma_case: float = None
"""Front casing area [m2]"""


a_tf_coil_nose_case: float = None
"""Nose casing area [m2]"""


a_tf_wp_ground_insulation: float = None
"""Inboard mid-plane cross-section area of the WP ground insulation [m2]"""


a_leg_ins: float = None
"""TF ouboard leg turn insulation area per coil [m2]"""


a_leg_gr_ins: float = None
"""TF outboard leg ground insulation area per coil [m2]"""


a_leg_cond: float = None
"""Exact TF ouboard leg conductor area [m2]"""


rad_tf_coil_inboard_toroidal_half: float = None
"""Half toroidal angular extent of a single TF coil inboard leg"""


tan_theta_coil: float = None
"""Tan half toroidal angular extent of a single TF coil inboard leg"""


t_conductor_radial: float = None
"""Conductor area radial and toroidal dimension (integer turn only) [m]"""


t_conductor_toroidal: float = None
"""Conductor area radial and toroidal dimension (integer turn only) [m]"""


dr_tf_turn_cable_space: float = None
"""Cable area radial and toroidal dimension (integer turn only) [m]"""


dx_tf_turn_cable_space: float = None
"""Cable area radial and toroidal dimension (integer turn only) [m]"""


dr_tf_turn: float = None
"""Turn radial and toroidal dimension (integer turn only) [m]"""


dx_tf_turn: float = None
"""Turn radial and toroidal dimension (integer turn only) [m]"""


dx_tf_turn_cable_space_average: float = None
"""Cable area averaged dimension (square shape) [m]"""


radius_tf_turn_cable_space_corners: float = None
"""Radius of the corners of the cable space in the TF turn [m]"""


a_tf_turn_cable_space_effective: float = None
"""True cable area of WP turn. This includes the removal of the cooling pipe [m^2] """

vforce_inboard_tot: float = None
"""Total inboard vertical tension (all coils) [N]"""

dr_tf_wp_no_insulation: float = None
"""Radial thickness of winding pack without insulation [m]"""

dia_tf_turn_superconducting_cable: float = None
"""Diameter of the superconducting cable in the TF turn [m]"""

j_tf_superconductor_critical: float = None
"""Critical current density of the superconducting cable [A/m^2]"""

f_c_tf_turn_operating_critical: float = None
"""Ratio of the TF operating current to the critical current"""

j_tf_coil_turn: float = None
"""Current density in the TF coil turn [A/m^2]"""

b_tf_superconductor_critical_zero_temp_strain: float = None
"""Critical magnetic field of the superconducting cable at zero temperature and strain [T]"""

temp_tf_superconductor_critical_zero_field_strain: float = None
"""Critical temperature of the superconducting cable at zero magnetic field and strain [K]"""

# Vacuum Vessel stress on TF coil quench

vv_stress_quench: float = None
"""The Tresca stress experienced by the Vacuum Vessel when the SCTF coil quenches [Pa]"""


# croco_strand

croco_strand_area: float = None
croco_strand_critical_current: float = None


# conductor

conductor_copper_area: float = None
conductor_copper_fraction: float = None
conductor_copper_bar_area: float = None
conductor_hastelloy_area: float = None
conductor_hastelloy_fraction: float = None
conductor_helium_area: float = None
conductor_helium_fraction: float = None
conductor_solder_area: float = None
conductor_solder_fraction: float = None
conductor_jacket_area: float = None
conductor_jacket_fraction: float = None
conductor_rebco_area: float = None
conductor_rebco_fraction: float = None
conductor_critical_current: float = None
conductor_acs: float = None
conductor_area: float = None
"""Area of cable space inside jacket"""


t1: float = None

time2: float = None

tau2: float = None

e_tf_magnetic_stored_total: float = None

is_leg_cp_temp_same: int = None


def init_superconducting_tf_coil_variables():
    global is_leg_cp_temp_same
    global tf_fit_t
    global tf_fit_z
    global f_b_tf_inboard_peak_ripple_symmetric
    global c_tf_coil
    global a_tf_wp_with_insulation
    global a_tf_wp_no_insulation
    global a_tf_coil_inboard_steel
    global a_tf_coil_inboard_insulation
    global f_a_tf_coil_inboard_steel
    global f_a_tf_coil_inboard_insulation
    global z_cp_top
    global r_tf_outboard_in
    global r_tf_outboard_out
    global r_tf_wp_inboard_inner
    global r_tf_wp_inboard_outer
    global r_tf_wp_inboard_centre
    global dr_tf_wp_top
    global vol_ins_cp
    global vol_gr_ins_cp
    global vol_case_cp
    global dx_tf_wp_toroidal_min
    global dx_tf_wp_toroidal_average
    global dx_tf_side_case_average
    global a_tf_plasma_case
    global a_tf_coil_nose_case
    global a_tf_wp_ground_insulation
    global a_leg_ins
    global a_leg_gr_ins
    global a_leg_cond
    global rad_tf_coil_inboard_toroidal_half
    global tan_theta_coil
    global t_conductor_radial
    global t_conductor_toroidal
    global dr_tf_turn_cable_space
    global dx_tf_turn_cable_space
    global dr_tf_turn
    global dx_tf_turn
    global dx_tf_turn_cable_space_average
    global vforce_inboard_tot
    global t1
    global time2
    global tau2
    global e_tf_magnetic_stored_total
    global radius_tf_turn_cable_space_corners
    global a_tf_turn_cable_space_effective
    global dr_tf_wp_no_insulation
    global dia_tf_turn_superconducting_cable
    global j_tf_superconductor_critical
    global f_c_tf_turn_operating_critical
    global j_tf_coil_turn
    global b_tf_superconductor_critical_zero_temp_strain
    global temp_tf_superconductor_critical_zero_field_strain

    is_leg_cp_temp_same = 0
    tf_fit_t = 0.0
    tf_fit_z = 0.0
    f_b_tf_inboard_peak_ripple_symmetric = 0.0
    c_tf_coil = 0.0
    a_tf_wp_with_insulation = 0.0
    a_tf_wp_no_insulation = 0.0
    a_tf_coil_inboard_steel = 0.0
    a_tf_coil_inboard_insulation = 0.0
    f_a_tf_coil_inboard_steel = 0.0
    f_a_tf_coil_inboard_insulation = 0.0
    z_cp_top = 0.0
    r_tf_outboard_in = 0.0
    r_tf_outboard_out = 0.0
    r_tf_wp_inboard_inner = 0.0
    r_tf_wp_inboard_outer = 0.0
    r_tf_wp_inboard_centre = 0.0
    dr_tf_wp_top = 0.0
    vol_ins_cp = 0.0
    vol_gr_ins_cp = 0.0
    vol_case_cp = 0.0
    dx_tf_wp_toroidal_min = 0.0
    dx_tf_wp_toroidal_average = 0.0
    dx_tf_side_case_average = 0.0
    a_tf_plasma_case = 0.0
    a_tf_coil_nose_case = 0.0
    a_tf_wp_ground_insulation = 0.0
    a_leg_ins = 0.0
    a_leg_gr_ins = 0.0
    a_leg_cond = 0.0
    rad_tf_coil_inboard_toroidal_half = 0.0
    tan_theta_coil = 0.0
    t_conductor_radial = 0.0
    t_conductor_toroidal = 0.0
    dr_tf_turn_cable_space = 0.0
    dx_tf_turn_cable_space = 0.0
    dr_tf_turn = 0.0
    dx_tf_turn = 0.0
    dx_tf_turn_cable_space_average = 0.0
    vforce_inboard_tot = 0.0
    t1 = 0.0
    time2 = 0.0
    tau2 = 0.0
    e_tf_magnetic_stored_total = 0.0
    radius_tf_turn_cable_space_corners = 0.0
    a_tf_turn_cable_space_effective = 0.0
    dr_tf_wp_no_insulation = 0.0
    dia_tf_turn_superconducting_cable = 0.0
    j_tf_superconductor_critical = 0.0
    f_c_tf_turn_operating_critical = 0.0
    j_tf_coil_turn = 0.0
    b_tf_superconductor_critical_zero_temp_strain = 0.0
    temp_tf_superconductor_critical_zero_field_strain = 0.0
