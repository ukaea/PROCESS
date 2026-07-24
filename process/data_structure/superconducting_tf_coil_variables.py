"""Module containing superconducting TF coil routines



N/A
This module contains routines for calculating the
parameters of a superconducting TF coil system for a
fusion power plant.
PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
"""

from dataclasses import dataclass


@dataclass(slots=True)
class SuperconductingTFData:
    tf_fit_t: float = 0.0
    """Dimensionless winding pack width"""

    tf_fit_z: float = 0.0
    """Dimensionless winding pack radial thickness"""

    f_b_tf_inboard_peak_ripple_symmetric: float = 0.0
    """Ratio of peak field with ripple to nominal axisymmetric peak field"""

    c_tf_coil: float = 0.0
    """Current in each TF coil"""

    a_tf_wp_with_insulation: float = 0.0
    """Total cross-sectional area of winding pack including
    GW insulation and insertion gap [m²]
    """

    a_tf_wp_no_insulation: float = 0.0
    """Total cross-sectional area of winding pack without
    ground insulation and insertion gap [m²]
    """

    a_tf_coil_inboard_steel: float = 0.0
    """Inboard coil steel coil cross-sectional area [m²]"""

    a_tf_coil_inboard_insulation: float = 0.0
    """Inboard coil insulation cross-section per coil [m²]"""

    f_a_tf_coil_inboard_steel: float = 0.0
    """Inboard coil steel fraction [-]"""

    f_a_tf_coil_inboard_insulation: float = 0.0
    """Inboard coil insulation fraction [-]"""

    z_cp_top: float = 0.0
    """Vertical distance from the midplane to the top of the tapered section [m]"""

    r_tf_outboard_in: float = 0.0
    """Radial position of plasma-facing edge of TF coil outboard leg [m]"""

    r_tf_outboard_out: float = 0.0
    """Radial position of outer edge of TF coil inboard leg [m]"""

    r_tf_wp_inboard_inner: float = 0.0
    """Radial position of inner edge and centre of winding pack [m]"""

    r_tf_wp_inboard_outer: float = 0.0
    """Radial position of outer edge and centre of winding pack [m]"""

    r_tf_wp_inboard_centre: float = 0.0
    """Radial position of centre and centre of winding pack [m]"""

    dr_tf_wp_top: float = 0.0
    """Conductor layer radial thickness at centre column top [m]
    Ground insulation layer included, only defined for itart = 1
    """

    vol_ins_cp: float = 0.0
    """CP turn insulation volume [m3]"""

    vol_gr_ins_cp: float = 0.0
    """CP ground insulation volume [m3]"""

    vol_case_cp: float = 0.0
    """Volume of the CP outer casing cylinder"""

    dx_tf_wp_toroidal_min: float = 0.0
    """Minimal toroidal thickness of of winding pack [m]"""

    dx_tf_wp_toroidal_average: float = 0.0
    """Averaged toroidal thickness of of winding pack [m]"""

    dx_tf_side_case_average: float = 0.0
    """Average lateral casing thickness [m]"""

    a_tf_plasma_case: float = 0.0
    """Front casing area [m²]"""

    a_tf_coil_nose_case: float = 0.0
    """Nose casing area [m²]"""

    a_tf_wp_ground_insulation: float = 0.0
    """Inboard mid-plane cross-section area of the WP ground insulation [m²]"""

    a_leg_ins: float = 0.0
    """TF outboard leg turn insulation area per coil [m²]"""

    a_leg_gr_ins: float = 0.0
    """TF outboard leg ground insulation area per coil [m²]"""

    a_leg_cond: float = 0.0
    """Exact TF outboard leg conductor area [m²]"""

    rad_tf_coil_inboard_toroidal_half: float = 0.0
    """Half toroidal angular extent of a single TF coil inboard leg"""

    tan_theta_coil: float = 0.0
    """Tan half toroidal angular extent of a single TF coil inboard leg"""

    t_conductor_radial: float = 0.0
    """Conductor area radial and toroidal dimension (integer turn only) [m]"""

    t_conductor_toroidal: float = 0.0
    """Conductor area radial and toroidal dimension (integer turn only) [m]"""

    dr_tf_turn_cable_space: float = 0.0
    """Cable area radial and toroidal dimension (integer turn only) [m]"""

    dx_tf_turn_cable_space: float = 0.0
    """Cable area radial and toroidal dimension (integer turn only) [m]"""

    dr_tf_turn: float = 0.0
    """Turn radial and toroidal dimension (integer turn only) [m]"""

    dx_tf_turn: float = 0.0
    """Turn radial and toroidal dimension (integer turn only) [m]"""

    dx_tf_turn_cable_space_average: float = 0.0
    """Cable area averaged dimension (square shape) [m]"""

    radius_tf_turn_cable_space_corners: float = 0.0
    """Radius of the corners of the cable space in the TF turn [m]"""

    a_tf_turn_cable_space_effective: float = 0.0
    """True cable area of WP turn. This includes the removal of the cooling pipe [m²] """

    vforce_inboard_tot: float = 0.0
    """Total inboard vertical tension (all coils) [N]"""

    dr_tf_wp_no_insulation: float = 0.0
    """Radial thickness of winding pack without insulation [m]"""

    dia_tf_turn_superconducting_cable: float = 0.00073
    """Diameter of the superconducting cable in the TF turn [m]"""

    dia_tf_turn_croco_cable: float = 0.0
    """Diameter of the Croco cable in the TF turn [m]"""

    n_tf_turn_superconducting_cables: int = 0
    """Number of superconducting cables in the TF turn"""

    len_tf_coil_superconductor: float = 0.0
    """Length of superconducting cable in one TF coil [m]"""

    len_tf_superconductor_total: float = 0.0
    """Total length of superconducting cable in all TF coils [m]"""

    j_tf_superconductor_critical: float = 0.0
    """Critical current density of the superconducting cable [A/m²]"""

    f_c_tf_turn_operating_critical: float = 0.0
    """Ratio of the TF operating current to the critical current"""

    j_tf_coil_turn: float = 0.0
    """Current density in the TF coil turn [A/m²]"""

    b_tf_superconductor_critical_zero_temp_strain: float = 0.0
    """Critical magnetic field of the superconducting cable at zero temperature and strain [T]"""

    temp_tf_superconductor_critical_zero_field_strain: float = 0.0
    """Critical temperature of the superconducting cable at zero magnetic field and strain [K]"""

    f_a_tf_turn_cable_space_cooling: float = 0.0
    """Fraction of usable turn cable space area used for cooling"""

    c_tf_turn_cables_critical: float = 0.0
    """Critical current density in the turn cables [A/m²]"""

    j_tf_superconductor: float = 0.0
    """Current density in the superconducting cable/tape/strand [A/m²]"""

    i_tf_turn_type: int = 1
    """Switch for TF turn geometry type"""

    # Vacuum Vessel stress on TF coil quench

    vv_stress_quench: float = 0.0
    """The Tresca stress experienced by the Vacuum Vessel when the SCTF coil quenches [Pa]"""

    # REBCO tape and CroCo strand variables

    dx_tf_hts_tape_rebco: float = 1.0e-6
    """thickness of REBCO layer in tape (m) (`iteration variable 138`)"""

    dx_tf_hts_tape_copper: float = 100.0e-6
    """thickness of copper layer in tape (m) (`iteration variable 139`)"""

    dx_tf_hts_tape_hastelloy: float = 50.0e-6
    """thickness of Hastelloy layer in tape (m)"""

    dr_tf_hts_tape: float = 4.0e-3
    """Mean width of tape (m)"""

    dx_tf_hts_tape_total: float = 6.5e-5
    """thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)"""

    dia_tf_croco_strand_tape_region: float = 0.0
    """Inner diameter of CroCo strand tape region (m)"""

    dx_tf_croco_strand_copper: float = 2.5e-3
    """Thickness of CroCo strand copper tube (m) (`iteration variable 158`)"""

    copper_rrr: float = 100.0
    """residual resistivity ratio copper in TF superconducting cable"""

    tf_coppera_m2: float = 0.0
    """TF coil current / copper area (A/m²)"""

    tf_coppera_m2_max: float = 1.0e8
    """Maximum TF coil current / copper area (A/m²)"""

    dx_tf_croco_strand_tape_stack: float = 0.0
    """Width / thickness of tape stack in CroCo strand (m)"""

    n_tf_croco_strand_hts_tapes: float = 0.0
    """Number of HTS tapes in CroCo strand"""

    a_tf_croco_strand_rebco: float = 0.0
    """Area of REBCO in CroCo strand (m²)"""

    a_tf_croco_strand_copper_total: float = 0.0
    """Area of copper in CroCo strand (includes tapes and outer tube) (m²)"""

    a_tf_croco_strand_hastelloy: float = 0.0
    """Area of Hastelloy in CroCo strand (m²)"""

    a_tf_croco_strand_solder: float = 0.0
    """Area of solder in CroCo strand (m²)"""

    a_tf_croco_strand: float = 0.0
    """Total area of a CroCo strand (m²)"""

    # croco_strand

    tf_croco_strand_area: float = 0.0
    cur_tf_turn_croco_strand_critical: float = 0.0
    """Critical current in the TF turn CroCo strand (A)"""

    # conductor

    a_tf_turn_croco_cable_space_copper: float = None
    conductor_copper_fraction: float = None
    a_tf_turn_croco_copper_bar: float = None
    """Area of the central copper strand in the CroCo TF turn [m²]"""
    a_tf_turn_croco_hastelloy: float = None
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
    conductor_area: float = None
    """Area of cable space inside jacket"""

    t1: float = 0.0

    time2: float = 0.0

    tau2: float = 0.0

    is_leg_cp_temp_same: int = 0


CREATE_DICTS_FROM_DATACLASS = SuperconductingTFData
