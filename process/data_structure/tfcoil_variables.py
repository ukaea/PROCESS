"""
Module containing global variables relating to the toroidal field coil systems
"""

from dataclasses import dataclass, field

import numpy as np

N_RADIAL_ARRAY = 50
"""Size of the radial distribution arrays per layers
used for stress, strain and displacement distibution
"""


@dataclass(slots=True)
class TFData:
    a_tf_coil_inboard_case: float = 0.0
    """external case area per coil (inboard leg) (m2)"""

    a_tf_coil_outboard_case: float = 0.0
    """external case area per coil (outboard leg) (m2)"""

    a_tf_turn_steel: float = 0.0
    """area of the cable conduit (m2)"""

    a_tf_wp_conductor: float = 0.0
    """Winding pack conductor area [m2]
    Does not include the area of voids and central helium channel
    """

    a_res_tf_coil_conductor: float = 0.0
    """Area of resistive conductor in resistive TF coil [m2]"""

    a_tf_turn_cable_space_no_void: float = 0.0
    """Cable space area (per turn)  [m2]
    Includes the area of voids and central helium channel
    """

    a_tf_turn_insulation: float = 0.0
    """single turn insulation area (m2)"""

    a_tf_coil_wp_turn_insulation: float = 0.0
    """winding pack turn insulation area per coil (m2)"""

    sig_tf_case_max: float = 6.0e8
    """Allowable maximum shear stress (Tresca criterion) in TF coil case (Pa)"""

    sig_tf_wp_max: float = 6.0e8
    """Allowable maximum shear stress (Tresca criterion) in TF coil conduit (Pa)"""

    a_tf_leg_outboard: float = 0.0
    """outboard TF leg area (m2)"""

    a_tf_wp_steel: float = 0.0
    """Total area of all winding pack steel (sum of each conduit steel from each turn) (m2)"""

    a_tf_wp_extra_void: float = 0.0
    """winding pack void (He coolant) area (m2)"""

    a_tf_wp_coolant_channels: float = 0.0
    """winding pack He coil area (m2)"""

    bcritsc: float = 24.0
    """upper critical field (T) for Nb3Sn superconductor at zero temperature and
    strain (`i_tf_sc_mat=4, =bc20m`)
    """

    b_tf_inboard_peak_symmetric: float = 0.0
    """mean peak field at TF coil (T)"""

    b_tf_inboard_peak_with_ripple: float = 0.0
    """peak field at TF conductor with ripple (T)"""

    casestr: float = 0.0
    """case strain"""

    dr_tf_plasma_case: float = 0.0
    """inboard TF coil case plasma side thickness (m) (calculated for stellarators)"""

    f_dr_tf_plasma_case: float = 0.05
    """inboard TF coil case plasma side thickness as a fraction of dr_tf_inboard"""

    i_f_dr_tf_plasma_case: bool = False
    """logical switch to make dr_tf_plasma_case a fraction of TF coil thickness (`f_dr_tf_plasma_case`)"""

    dx_tf_side_case_min: float = 0.0
    """inboard TF coil minimum sidewall case thickness (m) (calculated for stellarators)"""

    dx_tf_side_case_peak: float = 0.0
    """inboard TF coil peak sidewall case thickness (m) (calculated for stellarators)"""

    casths_fraction: float = 0.06
    """inboard TF coil sidewall case thickness as a fraction of dx_tf_inboard_out_toroidal"""

    tfc_sidewall_is_fraction: bool = False
    """logical switch to make dx_tf_side_case_min a fraction of TF coil thickness (`casths_fraction`)"""

    t_conductor: float = 0.0
    """Conductor (cable + steel conduit) area averaged dimension [m]"""

    dx_tf_turn_general: float = 0.0
    """TF coil turn edge length including turn insulation [m]
    If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
    equivalent size is use to calculated this quantity
    If the dx_tf_turn_general is non zero, c_tf_turn is calculated
    """

    i_dx_tf_turn_general_input: bool = False
    """Boolean switch to activated when the user set the TF coil turn dimensions
    Not an input
    """

    t_turn_tf_max: float = 0.05
    """TF turn edge length including turn insulation upper limit [m]
    If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
    equivalent size is use for this constraint
    constraint equation icc = 86
    """

    dx_tf_turn_cable_space_general: float = 0.0
    """TF coil superconducting cable squared/rounded dimensions [m]
    If the turn is not a square (i_tf_turns_integer = 1) a squared cable of
    equivalent size is use to calculated this quantity
    If the dx_tf_turn_cable_space_general is non zero, c_tf_turn is calculated
    """

    i_dx_tf_turn_cable_space_general_input: bool = False
    """Boolean switch to activated when the user set the TF coil cable dimensions
    Not an input
    """

    acs: float = 0.0
    """Area of space inside conductor (m2)"""

    cdtfleg: float = 0.0
    """TF outboard leg current density (A/m2) (resistive coils only)"""

    cforce: float = 0.0
    """centering force on inboard leg (per coil) (N/m)"""

    cplen: float = 0.0
    """length of TF coil inboard leg ('centrepost') (`i_tf_sup = 1`)"""

    c_tf_turn: float = 7.0e4
    """TF coil current per turn (A). (calculated for stellarators) (calculated for
    integer-turn TF coils `i_tf_turns_integer=1`) (`iteration variable 60`)
    """

    c_tf_turn_max: float = 9.0e4
    """Max TF coil current per turn [A]. (for stellarators and `i_tf_turns_integer=1`)
    (`constraint equation 77`)
    """

    den_tf_coil_case: float = 8000.0
    """density of coil case (kg/m3)"""

    dcond: list[float] = field(
        default_factory=lambda: np.array([
            6080.0,
            6080.0,
            6070.0,
            6080.0,
            6080.0,
            8500.0,
            6070.0,
            8500.0,
            8500.0,
        ])
    )
    """density of superconductor type given by i_tf_sc_mat/i_cs_superconductor/i_pf_superconductor (kg/m3)"""

    den_tf_wp_turn_insulation: float = 1800.0
    """density of conduit + ground-wall insulation (kg/m3)"""

    dia_tf_turn_coolant_channel: float = 0.005
    """diameter of central helium channel in TF winding (m)"""

    e_tf_magnetic_stored_total: float = 0.0
    """Total magnetic stored energy in the toroidal field coils (J)"""

    e_tf_magnetic_stored_total_gj: float = 0.0
    """Total magnetic stored energy in the toroidal field coils (GJ)"""

    e_tf_coil_magnetic_stored: float = 0.0
    """Stored magnetic energy in a single TF coil (J)"""

    b_crit_upper_nbti: float = 14.86
    """upper critical field of GL_nbti"""

    t_crit_nbti: float = 9.04
    """critical temperature of GL_nbti"""

    max_force_density: float = 0.0
    """Maximal (WP averaged) force density in TF coils at 1 point. (MN/m3)"""

    f_a_tf_turn_cable_copper: float = 0.69
    """copper fraction of cable conductor (TF coils)
    (iteration variable 59)
    """

    fhts: float = 0.5
    """technology adjustment factor for critical current density fit for isumat..=2
    Bi-2212 superconductor, to describe the level of technology assumed (i.e. to
    account for stress, fatigue, radiation, AC losses, joints or manufacturing
    variations; 1.0 would be very optimistic)
    """

    insstrain: float = 0.0
    """Radial strain in insulator"""

    i_tf_stress_model: int = 1
    """Switch for the TF coil stress model
    0 : Generalized plane strain formulation, Issues #977 and #991, O(n^3)
    1 : Old plane stress model (only for SC)
    2 : Axisymmetric extended plane strain, Issues #1414 and #998, O(n)
    """

    i_tf_tresca: int = 0
    """Switch for TF coil conduit Tresca stress criterion:
    0 : Tresca (no adjustment);
    1 : Tresca with CEA adjustment factors (radial+2%, vertical+60%) </UL>
    """

    i_tf_wp_geom: int = -1
    """Switch for TF WP geometry selection
    0 : Rectangular geometry
    1 : Double rectangular geometry
    2 : Trapezoidal geometry (constant lateral casing thickness)
    Default setting for backward compatibility
    if i_tf_turns_integer = 0 : Double rectangular
    if i_tf_turns_integer = 1 : Rectangular
    """

    i_tf_case_geom: int = 0
    """Switch for TF case geometry selection
    0 : Circular front case (ITER design)
    1 : Straight front case
    """

    i_tf_turns_integer: int = 0
    """Switch for TF coil integer/non-integer turns:
    0 : non-integer turns
    1 : integer turns
    """

    i_tf_sc_mat: int = 1
    """Switch for superconductor material in TF coils:
    - =1 ITER Nb3Sn critical surface model with standard
    ITER parameters
    - =2 Bi-2212 high temperature superconductor (range of
    validity T < 20K, adjusted field b < 104 T, B > 6 T)
    - =3 NbTi
    - =4 ITER Nb3Sn model with user-specified parameters
    - =5 WST Nb3Sn parameterisation
    - =6 REBCO HTS tape in CroCo strand
    - =7 Durham Ginzburg-Landau critical surface model for Nb-Ti
    - =8 Durham Ginzburg-Landau critical surface model for REBCO
    - =9 Hazelton experimental data + Zhai conceptual model for REBCO
    """

    i_tf_sup: int = 1
    """Switch for TF coil conductor model:
    - =0 copper
    - =1 superconductor
    - =2 Cryogenic aluminium
    """

    i_tf_shape: int = 0
    """Switch for TF coil toroidal shape:
    - =0  Default value : Picture frame coil for TART / PROCESS D-shape for non itart
    - =1  PROCESS D-shape : parametrise with 2 arcs
    - =2  Picture frame coils
    """

    i_tf_cond_eyoung_axial: int = 0
    """Switch for the behavior of the TF coil conductor elastic axial properties
    - =0  Young's modulus is set to zero, and the conductor is not considered
    in the stress calculation. This corresponds to the case that the
    conductor is much less stiff than the conduit, or the case that the
    conductor is prevented (isolated) from taking axial loads.
    - =1  Elastic properties are set by user input, using the variable
    `eyoung_cond_axial`
    - =2  Elastic properties are set to reasonable defaults taking into
    account the superconducting material `i_tf_sc_mat`
    """

    i_tf_cond_eyoung_trans: int = 1
    """Switch for the behavior of the elastic properties of the TF coil
    conductorin the transverse direction. Only active if
    `i_tf_cond_eyoung_axial == 2`
    - =0  Cable not potted in solder. Transverse Young's modulus set to zero.
    - =1  Cable potted in solder. If `i_tf_cond_eyoung_axial == 2`, the
    transverse Young's modulus of the conductor is equal to the axial,
    which is set to a sensible material-dependent default.
    """

    n_tf_wp_pancakes: int = 10
    """Number of pancakes in TF coil. Only used if `i_tf_turns_integer=1`"""

    n_tf_wp_layers: int = 20
    """Number of layers in TF coil. Only used if `i_tf_turns_integer=1`"""

    n_rad_per_layer: int = 100
    """Size of the arrays per layers storing the radial dependent stress
    quantities (stresses, strain displacement etc..)
    """

    i_tf_bucking: int = -1
    """Switch for TF inboard support structure design:
    Default setting for backward compatibility
    - if copper resistive TF (i_tf_sup = 0) : Free standing TF without bucking structure
    - if Superconducting TF  (i_tf_sup = 1) : Free standing TF with a steel casing
    - if aluminium  TF       (i_tf_sup = 2) : Free standing TF with a bucking structure
    Rem : the case is a bucking structure
    - =0 : Free standing TF without case/bucking cyliner (only a conductor layer)
    - =1 : Free standing TF with a case/bucking cylinder made of
    - if copper resistive     TF (i_tf_sup = 0) : used defined bucking cylinder
    - if Superconducting      TF (i_tf_sup = 1) : Steel casing
    - if aluminium resistive TF (i_tf_sup = 2) : used defined bucking cylinder
    - =2 : The TF is in contact with the CS : "bucked and wedged design"
    Fast version : thin TF-CS interface neglected in the stress calculations (3 layers)
    The CS is frictionally decoupled from the TF, does not carry axial tension
    - =3 : The TF is in contact with the CS : "bucked and wedged design"
    Full version : thin TF-CS Kapton interface introduced in the stress calculations (4 layers)
    The CS and kaptop are frictionally decoupled from the TF, do not carry
    axial tension
    """

    n_tf_graded_layers: int = 1
    """Number of layers of different stress properties in the WP. If `n_tf_graded_layers > 1`,
    a graded coil is condidered
    """

    n_tf_stress_layers: int = 0
    """Number of layers considered for the inboard TF stress calculations
    """

    n_tf_wp_stress_layers: int = 5
    """Maximum number of layers that can be considered in the TF coil composited/smeared
    stress analysis. This is the layers of one turn, not the entire WP.
    Default: 5. void, conductor, copper, conduit, insulation.
    """

    j_tf_bus: float = 1.25e6
    """bussing current density (A/m2)"""

    j_crit_str_tf: float = 0.0
    """j_crit_str : superconductor strand critical current density under operating
    conditions (A/m2). Necessary for the cost calculation in $/kAm
    """

    j_crit_str_0: list[float] = field(
        default_factory=lambda: np.array([
            596905475.80390120,
            1925501534.8512938,
            724544682.96063495,
            549858624.45072436,
            669284509.85818779,
            0.0,
            898964415.36996782,
            1158752995.2559297,
            865652122.9071957,
        ])
    )
    """j_crit_str_pf_0 : superconductor strand critical current density at 6 T and 4.2 K (A/m2)
    Necessary for the cost calculation in $/kAm
    """

    j_tf_wp_critical: float = 0.0
    """critical current density for winding pack (A/m2)"""

    j_tf_wp_quench_heat_max: float = 0.0
    """allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)"""

    j_tf_wp: float = 0.0
    """winding pack engineering current density (A/m2)"""

    j_tf_coil_full_area: float = 0.0
    """Inboard leg mid-plane full coil area current density (A/m²)
    """

    eyoung_ins: float = 1.0e8
    """Insulator Young's modulus [Pa]. Default value (1.0D8) setup the following values
    - SC TF, eyoung_ins = 20 Gpa (default value from DDD11-2 v2 2 (2009))
    - Al TF, eyoung_ins = 2.5 GPa (Kapton polymer)
    """

    eyoung_steel: float = 2.05e11
    """Steel case Young's modulus (Pa) (default value from DDD11-2 v2 2 (2009))"""

    eyoung_cond_axial: float = 6.6e8
    """SC TF coil conductor Young's modulus in the parallel (along the wire/tape)
    direction [Pa]
    Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
    set by the behavior of that switch.
    """

    eyoung_cond_trans: float = 0.0
    """SC TF coil conductor Young's modulus in the transverse direction [Pa]
    Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
    set by the behavior of that switch.
    """

    eyoung_res_tf_buck: float = 150.0e9
    """Resistive TF magnets bucking cylinder young modulus (Pa)"""

    eyoung_copper: float = 117.0e9
    """Copper young modulus. Default value taken from wikipedia"""

    eyoung_al: float = 69.0e9
    """Aluminium young modulus.  Default value taken from wikipedia"""

    poisson_steel: float = 0.3
    """Steel Poisson's ratio, Source : https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html"""

    poisson_copper: float = 0.35
    """Copper Poisson's ratio. Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html"""

    poisson_al: float = 0.35
    """Aluminium Poisson's ratio.
    Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
    """

    poisson_ins: float = 0.34
    """Insulation Poisson's ratio. Default: Kapton.
    Source : DuPont™ Kapton® HN datasheet.
    """

    poisson_cond_axial: float = 0.3
    """SC TF coil conductor Poisson's ratio in the parallel-transverse direction"""

    poisson_cond_trans: float = 0.3
    """SC TF coil conductor Poisson's ratio in the transverse-transverse direction"""

    r_b_tf_inboard_peak: float = 0.0
    """Radius of maximum TF B-field (m)"""

    res_tf_leg: float = 0.0
    """TF coil leg resistance (ohm)"""

    toroidalgap: float = 1.0
    """Minimal distance between two toroidal coils. (m)"""

    ripple_b_tf_plasma_edge_max: float = 1.0
    """maximum allowable toroidal field ripple amplitude at plasma edge (%)"""

    ripple_b_tf_plasma_edge: float = 0.0
    """peak/average toroidal field ripple at plasma edge (%)"""

    c_tf_total: float = 0.0
    """total (summed) current in TF coils (A)"""

    radial_array: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """Array refining the radii of the stress calculations arrays"""

    sig_tf_r: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """TF Inboard leg radial stress in steel r distribution at mid-plane [Pa]"""

    sig_tf_t: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """TF Inboard leg tangential stress in steel r distribution at mid-plane [Pa]"""

    deflect: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """TF coil radial deflection (displacement) radial distribution [m]"""

    sig_tf_z: float = 0.0
    """TF Inboard leg vertical tensile stress in steel at mid-plane [Pa]"""

    sig_tf_vmises: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """TF Inboard leg Von-Mises stress in steel r distribution at mid-plane [Pa]"""

    s_shear_tf: list[float] = field(default_factory=lambda: np.zeros(N_RADIAL_ARRAY))
    """TF Inboard leg maximum shear stress (Tresca criterion) in steel r distribution at mid-plane [Pa]"""

    sig_tf_cs_bucked: float = 0.0

    sig_tf_case: float = 0.0
    """Maximum shear stress (Tresca criterion) in TF casing steel structures (Pa)"""

    sig_tf_wp: float = 0.0

    str_cs_con_res: float = -0.005
    """Residual manufacturing strain in CS superconductor material"""

    str_pf_con_res: float = -0.005
    """Residual manufacturing strain in PF superconductor material"""

    str_tf_con_res: float = -0.005
    """Residual manufacturing strain in TF superconductor material
    If `i_str_wp == 0`, used to compute the critical surface.
    Otherwise, the self-consistent winding pack `str_wp` is used.
    """

    str_wp: float = 0.0
    """Axial (vertical) strain in the TF coil winding pack found by
    self-consistent stress/strain calculation.
    if `i_str_wp == 1`, used to compute the critical surface.
    Otherwise, the input value `str_tf_con_res` is used.
    Constrain the absolute value using `constraint equation 88`
    You can't have constraint 88 and i_str_wp = 0 at the same time
    """

    str_wp_max: float = 0.7e-2
    """Maximum allowed absolute value of the strain in the TF coil
    (`Constraint equation 88`)
    """

    i_str_wp: int = 1
    """Switch for the behavior of the TF strain used to compute
    the strain-dependent critical surface:
    - =0  str_tf_con_res is used
    - =1  str_wp is used
    """

    quench_model: str = "exponential"
    """switch for TF coil quench model (Only applies to REBCO magnet at present, issue #522):
    - ='exponential' exponential quench with constant discharge resistor
    - ='linear' quench with constant voltage
    """

    time1: float = 0
    """Time at which TF quench is detected (s)"""

    tcritsc: float = 16.0
    """critical temperature (K) for superconductor at zero field and strain (`i_tf_sc_mat=4, =tc0m`)"""

    t_tf_superconductor_quench: float = 10.0
    """fast discharge time for TF coil in event of quench (s) (`iteration variable 56`)
    For REBCO model, meaning depends on quench_model:
    - exponential quench : e-folding time (s)`
    - linear quench : discharge time (s)
    """

    a_tf_inboard_total: float = 0.0
    """Total inboard area of all TF coils (m²)"""

    len_tf_bus: float = 300.0
    """TF coil bus length (m)"""

    m_tf_bus: float = 0.0
    """TF coil bus mass (kg)"""

    tfckw: float = 0.0
    """available DC power for charging the TF coils (kW)"""

    tfcmw: float = 0.0
    """Peak power per TF power supply (MW)"""

    p_cp_resistive_mw: float = 0.0
    """Peak resistive TF coil inboard leg power (MW)"""

    p_tf_joints_resistive_mw: float = 0.0
    """TF joints resistive power losses (MW)"""

    tfcryoarea: float = 0.0
    """surface area of toroidal shells covering TF coils (m2)"""

    tficrn: float = 0.0
    """TF coil half-width - inner dr_bore (m)"""

    ind_tf_coil: float = 0.0
    """TF coil inductance (H)"""

    dx_tf_wp_insertion_gap: float = 0.01
    """TF coil WP insertion gap (m)"""

    p_tf_leg_resistive_mw: float = 0.0
    """TF coil outboard leg resistive power (MW)"""

    rho_cp: float = 0.0
    """TF coil inboard leg resistivity [Ohm-m]. If `itart=0`, this variable is the
    average resistivity over the whole magnet
    """

    rho_tf_leg: float = 0.0
    """Resistivity of a TF coil leg (Ohm-m)"""

    rho_tf_bus: float = 1.86e-8
    """Resistivity of a TF coil bus (Ohm-m). Default values is for that of GLIDCOP AL-15 (C15715) at 293K"""

    frhocp: float = 1.0
    """Centrepost resistivity enhancement factor. For `itart=0`, this factor
    is used for the whole magnet
    """

    frholeg: float = 1.0
    """Outboard legs resistivity enhancement factor. Only used for `itart=1`."""

    i_cp_joints: int = -1
    """Switch for CP demoutable joints type
    -= 0 : Clampled joints
    -= 1 : Sliding joints
    Default value (-1) choses :
    Sliding joints for resistive magnets (i_tf_sup = 0, 2)
    Clampled joints for superconducting magnets (i_tf_sup = 1)
    """

    rho_tf_joints: float = 2.5e-10
    """TF joints surfacic resistivity [ohm.m]. Feldmetal joints assumed."""

    n_tf_joints_contact: int = 6
    """Number of contact per turn"""

    n_tf_joints: int = 4
    """Number of joints
    Ex: n_tf_joints = 2 for top and bottom CP joints
    """

    th_joint_contact: float = 0.03
    """TF sliding joints contact pad width [m]"""

    p_tf_joints_resistive: float = 0.0
    """Calculated TF joints resistive power losses [W]"""

    len_tf_coil: float = 0.0
    """TF coil circumference (m)"""

    eff_tf_cryo: float = -1.0
    """TF cryoplant efficiency (compared to pefect Carnot cycle).
    Using -1 set the default value depending on magnet technology:
    - i_tf_sup = 1 : SC magnet, eff_tf_cryo = 0.13 (ITER design)
    - i_tf_sup = 2 : Cryo-aluminium, eff_tf_cryo = 0.4
    """

    n_tf_coils: float = 16.0
    """Number of TF coils (default = 50 for stellarators). Number of TF coils outer legs for ST"""

    tfocrn: float = 0.0
    """TF coil half-width - outer dr_bore (m)"""

    tfsai: float = 0.0
    """area of the inboard TF coil legs (m2)"""

    tfsao: float = 0.0
    """area of the outboard TF coil legs (m2)"""

    tftmp: float = 4.5
    """peak helium coolant temperature in TF coils and PF coils (K)"""

    dx_tf_inboard_out_toroidal: float = 1.0
    """Inboard leg toroidal thickness at outer edge (m)"""

    dx_tf_turn_insulation: float = 8e-4
    """conduit insulation thickness (m)"""

    layer_ins: float = 0.0
    """Additional insulation thickness between layers (m)"""

    dr_tf_nose_case: float = 0.3
    """inboard TF coil case outer (non-plasma side) thickness (m) (`iteration variable 57`)
    (calculated for stellarators)
    """

    dr_tf_full_midplane: float = 0.0
    """Full radial thickness of TF coil at midplane (m)"""

    dr_tf_internal_midplane: float = 0.0
    """Internal radial thickness of TF coil at midplane (m)"""

    dr_tf_wp_with_insulation: float = 0.0
    """radial thickness of winding pack (m) (`iteration variable 140`) (issue #514)"""

    dx_tf_turn_steel: float = 8e-3
    """TF coil turn steel conduit case thickness (m) (`iteration variable 58`)"""

    dx_tf_wp_insulation: float = 0.018
    """Thickness of the ground insulation layer surrounding (m)
    - Superconductor TF (`i_tf_sup == 1`) : The TF coil Winding packs
    - Resistive magnets (`i_tf_sup /= 1`) : The TF coil wedges
    Rem : Thickness calculated for stellarators.
    """

    temp_tf_superconductor_margin_min: float = 0.0
    """minimum allowable temperature margin : TF coils (K)"""

    temp_cs_superconductor_margin_min: float = 0.0
    """minimum allowable temperature margin : CS (K)"""

    tmargmin: float = 0.0
    """minimum allowable temperature margin : TFC AND CS (K)"""

    temp_margin: float = 0.0
    """temperature margin (K)"""

    temp_tf_superconductor_margin: float = 0.0
    """TF coil superconductor temperature margin (K)"""

    temp_tf_conductor_quench_max: float = 150.0
    """maximum temp during a quench for protection (K)"""

    temp_croco_quench_max: float = 200.0
    """CroCo strand: maximum permitted temp during a quench (K)"""

    temp_croco_quench: float = 0.0
    """CroCo strand: Actual temp reached during a quench (K)"""

    temp_tf_cryo: float = 4.5
    """coil temperature for cryogenic plant power calculation (K)"""

    n_tf_coil_turns: float = 0.0
    """number of turns per TF coil"""

    v_tf_coil_dump_quench_max_kv: float = 20.0
    """max voltage across TF coil during quench (kV) (`iteration variable 52`)"""

    vforce: float = 0.0
    """vertical tension on inboard leg/coil (N)"""

    f_vforce_inboard: float = 0.5
    """Fraction of the total vertical force taken by the TF inboard leg tension
    Not used for resistive `itart=1` (sliding joints)
    """

    vforce_outboard: float = 0.0
    """Vertical tension on outboard leg/coil (N)"""

    f_a_tf_turn_cable_space_extra_void: float = 0.4
    """coolant fraction of TFC 'cable' (`i_tf_sup=1`), or of TFC leg (`i_tf_ssup=0`)"""

    voltfleg: float = 0.0
    """volume of each TF coil outboard leg (m3)"""

    vtfkv: float = 0.0
    """TF coil voltage for resistive coil including bus (kV)"""

    v_tf_coil_dump_quench_kv: float = 0.0
    """voltage across a TF coil during quench (kV)"""

    m_tf_coil_case: float = 0.0
    """mass per coil of external case (kg)"""

    m_tf_coil_conductor: float = 0.0
    """TF coil conductor mass per coil (kg/coil).
    For `itart=1`, coil is return limb plus centrepost/n_tf_coils
    """

    m_tf_coil_copper: float = 0.0
    """copper mass in TF coil conductor (kg/coil).
    For `itart=1`, coil is return limb plus centrepost/n_tf_coils
    """

    whtconal: float = 0.0
    """Aluminium mass in TF coil conductor (kg/coil).
    For `itart=1`, coil is return limb plus centrepost/n_tf_coils
    """

    m_tf_coil_wp_turn_insulation: float = 0.0
    """conduit insulation mass in TF coil conductor (kg/coil)"""

    m_tf_coil_superconductor: float = 0.0
    """superconductor mass in TF coil cable (kg/coil)"""

    m_tf_wp_steel_conduit: float = 0.0
    """steel conduit mass in TF coil conductor (kg/coil)"""

    m_tf_coil_wp_insulation: float = 0.0
    """mass of ground-wall insulation layer per coil (kg/coil)"""

    m_tf_coils_total: float = 0.0
    """total mass of the TF coils (kg)"""

    dx_tf_wp_primary_toroidal: float = 0.0
    """width of first step of winding pack (m)"""

    dx_tf_wp_secondary_toroidal: float = 0.0
    """width of second step of winding pack (m)"""

    # Superconducting TF coil shape parameters;
    # the TF inner surface top half is approximated by four circular arcs.
    # Arc 1 goes through points 1 and 2 on the inner surface. Arc 2
    # goes through points 2 and 3, etc.

    dthet: list[float] = field(default_factory=lambda: np.zeros(4))
    """angle of arc i (rad)"""

    radctf: list[float] = field(default_factory=lambda: np.zeros(4))
    """radius of arc i (m)"""

    r_tf_arc: list[float] = field(default_factory=lambda: np.zeros(5))
    """x location of arc point i on surface (m)"""

    xctfc: list[float] = field(default_factory=lambda: np.zeros(4))
    """x location of arc centre i (m)"""

    z_tf_arc: list[float] = field(default_factory=lambda: np.zeros(5))
    """y location of arc point i on surface (m)"""

    yctfc: list[float] = field(default_factory=lambda: np.zeros(4))
    """y location of arc centre i (m)"""

    # New TF shape:  Horizontal and vertical radii of inside edge of TF coil
    # Arcs are numbered clockwise:
    # 1=upper inboard, 2=upper outboard, 3=lower ouboard, 4=lower inboard

    tfa: list[float] = field(default_factory=lambda: np.zeros(4))
    """Horizontal radius of inside edge of TF coil (m)"""

    tfb: list[float] = field(default_factory=lambda: np.zeros(4))
    """Vertical radius of inside edge of TF coil (m)"""

    # Quantities relating to the spherical tokamak model (itart=1)
    # (and in some cases, also to resistive TF coils, i_tf_sup=0):

    drtop: float = 0.0
    """centrepost taper maximum radius adjustment (m)"""

    dztop: float = 0.0
    """centrepost taper height adjustment (m)"""

    etapump: float = 0.8
    """centrepost coolant pump efficiency"""

    fcoolcp: float = 0.3
    """coolant fraction of TF coil inboard legs (`iteration variable 23`)"""

    f_a_tf_cool_outboard: float = 0.2
    """coolant fraction of TF coil outboard legs"""

    a_cp_cool: float = 0.0
    """Centrepost cooling area toroidal cross-section (constant over the whole CP)"""

    n_cp_coolant_channels_total: float = 0.0
    """number of centrepost coolant tubes"""

    p_cp_coolant_pump_elec: float = 0.0
    """centrepost coolant pump power (W)"""

    p_cp_resistive: float = 0.0
    """resistive power in the centrepost (itart=1) [W].
    If `itart=0`, this variable is the ressitive power on the whole magnet
    """

    p_tf_leg_resistive: float = 0.0
    """Summed resistive power in the TF coil legs [W]. Remain 0 if `itart=0`."""

    temp_cp_max: float = 473.15  # 200 C
    """maximum peak centrepost temperature (K) (`constraint equation 44`)"""

    radius_cp_coolant_channel: float = 0.005
    """average radius of coolant channel (m) (`iteration variable 69`)"""

    temp_cp_coolant_inlet: float = 313.15  # 40 C
    """centrepost coolant inlet temperature (K)"""

    dtemp_cp_coolant: float = 0.0
    """inlet / outlet TF coil coolant temperature rise (K)"""

    temp_cp_average: float = 373.15  # 100 C
    """Average temperature of centrepost called CP (K). Only used for resistive coils
    to compute the resisitive heating. Must be an iteration variable for
    ST (`itart=1`) (`iteration variable 20`)
    """

    tcpav2: float = 0.0
    """Computed centrepost average temperature (K) (for consistency)"""

    temp_tf_legs_outboard: float = -1.0
    """Average temperature of the TF outboard legs [K]. If `temp_tf_legs_outboard=-1.0`, the ouboard
    legs and CP temperatures are the same. Fixed for now, should use a contraints eq like temp_cp_average
    """

    temp_cp_peak: float = 0.0
    """peak centrepost temperature (K)"""

    vel_cp_coolant_midplane: float = 20.0
    """inlet centrepost coolant flow speed at midplane (m/s) (`iteration variable 70`)"""

    vol_cond_cp: float = 0.0
    """Exact conductor volume in the centrepost (m3)"""

    whtcp: float = 0.0
    """mass of TF coil inboard legs (kg)"""

    whttflgs: float = 0.0
    """mass of the TF coil legs (kg)"""

    cryo_cool_req: float = 0.0
    """Cryo cooling requirement at helium temp 4.5K (kW)"""

    theta1_coil: float = 45.0
    """The angle of the outboard arc forming the TF coil current center line [deg]"""

    theta1_vv: float = 1.0
    """The angle of the outboard arc forming the Vacuum Vessel current center line [deg]"""

    max_vv_stress: float = 143.0e6
    """The allowable peak maximum shear stress in the vacuum vessel due to quench and fast discharge of the TF coils [Pa]"""

    t_tf_quench_detection: float = 3.0
    """TF coil quench detection time (s). Only used for TF coil quench protection."""

    rrr_tf_cu: float = 100.0
    """TF coil copper residual-resistance-ratio (RRR). Only used for quench protection."""

    a_tf_turn: float = 0.0
    """TF coil turn area (m²)"""

    drarea: float = 0.0

    r_b_tf_inboard_peak_symmetric: float = 0.0


CREATE_DICTS_FROM_DATACLASS = TFData
