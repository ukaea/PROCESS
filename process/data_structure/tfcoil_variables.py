"""author: J. Morris, M. Kovari, S. Kahn (UKAEA)
Module containing global variables relating to the toroidal field coil systems
### References
- ITER Magnets design description document DDD11-2 v2 2 (2009)
"""

import numpy as np

N_RADIAL_ARRAY = 50
"""Size of the radial distribution arrays per layers
used for stress, strain and displacement distibution
"""


a_tf_coil_inboard_case: float = None
"""external case area per coil (inboard leg) (m2)"""


a_tf_coil_outboard_case: float = None
"""external case area per coil (outboard leg) (m2)"""


a_tf_turn_steel: float = None
"""area of the cable conduit (m2)"""


a_tf_wp_conductor: float = None
"""Winding pack conductor area [m2]
Does not include the area of voids and central helium channel
"""


a_res_tf_coil_conductor: float = None
"""Area of resistive conductor in resistive TF coil [m2]"""


a_tf_turn_cable_space_no_void: float = None
"""Cable space area (per turn)  [m2]
Includes the area of voids and central helium channel
"""


a_tf_turn_insulation: float = None
"""single turn insulation area (m2)"""


a_tf_coil_wp_turn_insulation: float = None
"""winding pack turn insulation area per coil (m2)"""


sig_tf_case_max: float = None
"""Allowable maximum shear stress (Tresca criterion) in TF coil case (Pa)"""


sig_tf_wp_max: float = None
"""Allowable maximum shear stress (Tresca criterion) in TF coil conduit (Pa)"""


a_tf_leg_outboard: float = None
"""outboard TF leg area (m2)"""


a_tf_wp_steel: float = None
"""Total area of all winding pack steel (sum of each conduit steel from each turn) (m2)"""


a_tf_wp_extra_void: float = None
"""winding pack void (He coolant) area (m2)"""


a_tf_wp_coolant_channels: float = None
"""winding pack He coil area (m2)"""


bcritsc: float = None
"""upper critical field (T) for Nb3Sn superconductor at zero temperature and
strain (`i_tf_sc_mat=4, =bc20m`)
"""


b_tf_inboard_peak_symmetric: float = None
"""mean peak field at TF coil (T)"""


b_tf_inboard_peak_with_ripple: float = None
"""peak field at TF conductor with ripple (T)"""


casestr: float = None
"""case strain"""


dr_tf_plasma_case: float = None
"""inboard TF coil case plasma side thickness (m) (calculated for stellarators)"""


f_dr_tf_plasma_case: float = None
"""inboard TF coil case plasma side thickness as a fraction of dr_tf_inboard"""


i_f_dr_tf_plasma_case: bool = None
"""logical switch to make dr_tf_plasma_case a fraction of TF coil thickness (`f_dr_tf_plasma_case`)"""


dx_tf_side_case_min: float = None
"""inboard TF coil minimum sidewall case thickness (m) (calculated for stellarators)"""

dx_tf_side_case_peak: float = None
"""inboard TF coil peak sidewall case thickness (m) (calculated for stellarators)"""


casths_fraction: float = None
"""inboard TF coil sidewall case thickness as a fraction of dx_tf_inboard_out_toroidal"""


tfc_sidewall_is_fraction: bool = None
"""logical switch to make dx_tf_side_case_min a fraction of TF coil thickness (`casths_fraction`)"""


t_conductor: float = None
"""Conductor (cable + steel conduit) area averaged dimension [m]"""


dx_tf_turn_general: float = None
"""TF coil turn edge length including turn insulation [m]
If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
equivelent size is use to calculated this quantity
If the dx_tf_turn_general is non zero, c_tf_turn is calculated
"""


i_dx_tf_turn_general_input: bool = None
"""Boolean switch to activated when the user set the TF coil turn dimensions
Not an input
"""


f_t_turn_tf: float = None
"""f-value for TF turn edge length constraint
If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
equivelent size is use for this constraint
iteration variable ixc = 175
constraint equation icc = 86
"""


t_turn_tf_max: float = None
"""TF turn edge length including turn insulation upper limit [m]
If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
equivelent size is use for this constraint
constraint equation icc = 86
"""


dx_tf_turn_cable_space_general: float = None
"""TF coil superconducting cable squared/rounded dimensions [m]
If the turn is not a square (i_tf_turns_integer = 1) a squared cable of
equivelent size is use to calculated this quantity
If the dx_tf_turn_cable_space_general is non zero, c_tf_turn is calculated
"""


i_dx_tf_turn_cable_space_general_input: bool = None
"""Boolean switch to activated when the user set the TF coil cable dimensions
Not an input
"""


acs: float = None
"""Area of space inside conductor (m2)"""


cdtfleg: float = None
"""TF outboard leg current density (A/m2) (resistive coils only)"""


cforce: float = None
"""centering force on inboard leg (per coil) (N/m)"""


cplen: float = None
"""length of TF coil inboard leg ('centrepost') (`i_tf_sup = 1`)"""


c_tf_turn: float = None
"""TF coil current per turn (A). (calculated for stellarators) (calculated for
integer-turn TF coils `i_tf_turns_integer=1`) (`iteration variable 60`)
"""


c_tf_turn_max: float = None
"""Max TF coil current per turn [A]. (for stellarators and `i_tf_turns_integer=1`)
(`constraint equation 77`)
"""


den_tf_coil_case: float = None
"""density of coil case (kg/m3)"""


dcond: list[float] = None
"""density of superconductor type given by i_tf_sc_mat/i_cs_superconductor/i_pf_superconductor (kg/m3)"""


den_tf_wp_turn_insulation: float = None
"""density of conduit + ground-wall insulation (kg/m3)"""


dia_tf_turn_coolant_channel: float = None
"""diameter of central helium channel in TF winding (m)"""


e_tf_magnetic_stored_total_gj: float = None
"""total magnetic stored energy in the toroidal field coils (GJ)"""

e_tf_coil_magnetic_stored: float = None
"""Stored magnetic energy in a single TF coil (J)"""


b_crit_upper_nbti: float = None
"""upper critical field of GL_nbti"""


t_crit_nbti: float = None
"""critical temperature of GL_nbti"""


max_force_density: float = None
"""Maximal (WP averaged) force density in TF coils at 1 point. (MN/m3)"""


f_a_tf_turn_cable_copper: float = None
"""copper fraction of cable conductor (TF coils)
(iteration variable 59)
"""


fhts: float = None
"""technology adjustment factor for critical current density fit for isumat..=2
Bi-2212 superconductor, to describe the level of technology assumed (i.e. to
account for stress, fatigue, radiation, AC losses, joints or manufacturing
variations; 1.0 would be very optimistic)
"""


insstrain: float = None
"""Radial strain in insulator"""


i_tf_stress_model: int = None
"""Switch for the TF coil stress model
0 : Generalized plane strain formulation, Issues #977 and #991, O(n^3)
1 : Old plane stress model (only for SC)
2 : Axisymmetric extended plane strain, Issues #1414 and #998, O(n)
"""


i_tf_tresca: int = None
"""Switch for TF coil conduit Tresca stress criterion:
0 : Tresca (no adjustment);
1 : Tresca with CEA adjustment factors (radial+2%, vertical+60%) </UL>
"""


i_tf_wp_geom: int = None
"""Switch for TF WP geometry selection
0 : Rectangular geometry
1 : Double rectangular geometry
2 : Trapezoidal geometry (constant lateral casing thickness)
Default setting for backward compatibility
if i_tf_turns_integer = 0 : Double rectangular
if i_tf_turns_integer = 1 : Rectangular
"""


i_tf_case_geom: int = None
"""Switch for TF case geometry selection
0 : Circular front case (ITER design)
1 : Straight front case
"""


i_tf_turns_integer: int = None
"""Switch for TF coil integer/non-integer turns:
0 : non-integer turns
1 : integer turns
"""


i_tf_sc_mat: int = None
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


i_tf_sup: int = None
"""Switch for TF coil conductor model:
- =0 copper
- =1 superconductor
- =2 Cryogenic aluminium
"""


i_tf_shape: int = None
"""Switch for TF coil toroidal shape:
- =0  Default value : Picture frame coil for TART / PROCESS D-shape for non itart
- =1  PROCESS D-shape : parametrise with 2 arcs
- =2  Picture frame coils
"""


i_tf_cond_eyoung_axial: int = None
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


i_tf_cond_eyoung_trans: int = None
"""Switch for the behavior of the elastic properties of the TF coil
conductorin the transverse direction. Only active if
`i_tf_cond_eyoung_axial == 2`
- =0  Cable not potted in solder. Transverse Young's modulus set to zero.
- =1  Cable potted in solder. If `i_tf_cond_eyoung_axial == 2`, the
transverse Young's modulus of the conductor is equal to the axial,
which is set to a sensible material-dependent default.
"""


n_tf_wp_pancakes: int = None
"""Number of pancakes in TF coil. Only used if `i_tf_turns_integer=1`"""


n_tf_wp_layers: int = None
"""Number of layers in TF coil. Only used if `i_tf_turns_integer=1`"""


n_rad_per_layer: int = None
"""Size of the arrays per layers storing the radial dependent stress
quantities (stresses, strain displacement etc..)
"""


i_tf_bucking: int = None
"""Switch for TF inboard suport structure design:
Default setting for backward compatibility
- if copper resistive TF (i_tf_sup = 0) : Free standing TF without bucking structure
- if Superconducting TF  (i_tf_sup = 1) : Free standing TF with a steel casing
- if aluminium  TF       (i_tf_sup = 2) : Free standing TF with a bucking structure
Rem : the case is a bucking structure
- =0 : Free standing TF without case/bucking cyliner (only a conductor layer)
- =1 : Free standing TF with a case/bucking cylinder made of
- if copper resistive     TF (i_tf_sup = 0) : used defined bucking cylinder
- if Superconducting      TF (i_tf_sup = 1) : Steel casing
- if aluminium resisitive TF (i_tf_sup = 2) : used defined bucking cylinder
- =2 : The TF is in contact with the CS : "bucked and wedged design"
Fast version : thin TF-CS interface neglected in the stress calculations (3 layers)
The CS is frictionally decoupled from the TF, does not carry axial tension
- =3 : The TF is in contact with the CS : "bucked and wedged design"
Full version : thin TF-CS Kapton interface introduced in the stress calculations (4 layers)
The CS and kaptop are frictionally decoupled from the TF, do not carry
axial tension
"""


n_tf_graded_layers: int = None
"""Number of layers of different stress properties in the WP. If `n_tf_graded_layers > 1`,
a graded coil is condidered
"""


n_tf_stress_layers: int = None
"""Number of layers considered for the inboard TF stress calculations
set in initial.f90 from i_tf_bucking and n_tf_graded_layers
"""


n_tf_wp_stress_layers: int = None
"""Maximum number of layers that can be considered in the TF coil composited/smeared
stress analysis. This is the layers of one turn, not the entire WP.
Default: 5. void, conductor, copper, conduit, insulation.
"""


j_tf_bus: float = None
"""bussing current density (A/m2)"""


j_crit_str_tf: float = None
"""j_crit_str : superconductor strand critical current density under operating
conditions (A/m2). Necessary for the cost calculation in $/kAm
"""


j_crit_str_0: list[float] = None
"""j_crit_str_pf_0 : superconductor strand critical current density at 6 T and 4.2 K (A/m2)
Necessary for the cost calculation in $/kAm
"""


j_tf_wp_critical: float = None
"""critical current density for winding pack (A/m2)"""


j_tf_wp_quench_heat_max: float = None
"""allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)"""


j_tf_wp: float = None
"""winding pack engineering current density (A/m2)"""


oacdcp: float = None
"""Overall current density in TF coil inboard legs midplane (A/m2)
Rem SK : Not used in tfcoil to set the current any more. Should not be used as
iteration variable 12 any more. It is now calculated.
"""


eyoung_ins: float = None
"""Insulator Young's modulus [Pa]. Default value (1.0D8) setup the following values
- SC TF, eyoung_ins = 20 Gpa (default value from DDD11-2 v2 2 (2009))
- Al TF, eyoung_ins = 2.5 GPa (Kapton polymer)
"""


eyoung_steel: float = None
"""Steel case Young's modulus (Pa) (default value from DDD11-2 v2 2 (2009))"""


eyoung_cond_axial: float = None
"""SC TF coil conductor Young's modulus in the parallel (along the wire/tape)
direction [Pa]
Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
set by the behavior of that switch.
"""


eyoung_cond_trans: float = None
"""SC TF coil conductor Young's modulus in the transverse direction [Pa]
Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
set by the behavior of that switch.
"""


eyoung_res_tf_buck: float = None
"""Resistive TF magnets bucking cylinder young modulus (Pa)"""


eyoung_copper: float = None
"""Copper young modulus. Default value taken from wikipedia"""


eyoung_al: float = None
"""Aluminium young modulus.  Default value taken from wikipedia"""


poisson_steel: float = None
"""Steel Poisson's ratio, Source : https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html"""


poisson_copper: float = None
"""Copper Poisson's ratio. Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html"""


poisson_al: float = None
"""Aluminium Poisson's ratio.
Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
"""


poisson_ins: float = None
"""Insulation Poisson's ratio. Default: Kapton.
Source : DuPont™ Kapton® HN datasheet.
"""


poisson_cond_axial: float = None
"""SC TF coil conductor Poisson's ratio in the parallel-transverse direction"""


poisson_cond_trans: float = None
"""SC TF coil conductor Poisson's ratio in the transverse-transverse direction"""


r_b_tf_inboard_peak: float = None
"""Radius of maximum TF B-field (m)"""


res_tf_leg: float = None
"""TF coil leg resistance (ohm)"""


toroidalgap: float = None
"""Minimal distance between two toroidal coils. (m)"""


ftoroidalgap: float = None
"""F-value for minimum dx_tf_inboard_out_toroidal (`constraint equation 82`)"""


ripple_b_tf_plasma_edge_max: float = None
"""maximum allowable toroidal field ripple amplitude at plasma edge (%)"""


ripple_b_tf_plasma_edge: float = None
"""peak/average toroidal field ripple at plasma edge (%)"""


c_tf_total: float = None
"""total (summed) current in TF coils (A)"""


radial_array: list[float] = None
"""Array refining the radii of the stress calculations arrays"""


sig_tf_r: list[float] = None
"""TF Inboard leg radial stress in steel r distribution at mid-plane [Pa]"""


sig_tf_t: list[float] = None
"""TF Inboard leg tangential stress in steel r distribution at mid-plane [Pa]"""


deflect: list[float] = None
"""TF coil radial deflection (displacement) radial distribution [m]"""


sig_tf_z: float = None
"""TF Inboard leg vertical tensile stress in steel at mid-plane [Pa]"""


sig_tf_vmises: list[float] = None
"""TF Inboard leg Von-Mises stress in steel r distribution at mid-plane [Pa]"""


s_shear_tf: list[float] = None
"""TF Inboard leg maximum shear stress (Tresca criterion) in steel r distribution at mid-plane [Pa]"""


sig_tf_cs_bucked: float = None


sig_tf_case: float = None
"""Maximum shear stress (Tresca criterion) in TF casing steel structures (Pa)"""


sig_tf_wp: float = None


str_cs_con_res: float = None
"""Residual manufacturing strain in CS superconductor material"""


str_pf_con_res: float = None
"""Residual manufacturing strain in PF superconductor material"""


str_tf_con_res: float = None
"""Residual manufacturing strain in TF superconductor material
If `i_str_wp == 0`, used to compute the critical surface.
Otherwise, the self-consistent winding pack `str_wp` is used.
"""


str_wp: float = None
"""Axial (vertical) strain in the TF coil winding pack found by
self-consistent stress/strain calculation.
if `i_str_wp == 1`, used to compute the critical surface.
Otherwise, the input value `str_tf_con_res` is used.
Constrain the absolute value using `constraint equation 88`
You can't have constraint 88 and i_str_wp = 0 at the same time
"""


str_wp_max: float = None
"""Maximum allowed absolute value of the strain in the TF coil
(`Constraint equation 88`)
"""


i_str_wp: int = None
"""Switch for the behavior of the TF strain used to compute
the strain-dependent critical surface:
- =0  str_tf_con_res is used
- =1  str_wp is used
"""


quench_model: str = None
"""switch for TF coil quench model (Only applies to REBCO magnet at present, issue #522):
- ='exponential' exponential quench with constant discharge resistor
- ='linear' quench with constant voltage
"""


time1: float = None
"""Time at which TF quench is detected (s)"""


tcritsc: float = None
"""critical temperature (K) for superconductor at zero field and strain (`i_tf_sc_mat=4, =tc0m`)"""


t_tf_superconductor_quench: float = None
"""fast discharge time for TF coil in event of quench (s) (`iteration variable 56`)
For REBCO model, meaning depends on quench_model:
- exponential quench : e-folding time (s)`
- linear quench : discharge time (s)
"""


a_tf_inboard_total: float = None
"""Area of inboard midplane TF legs (m2)"""


len_tf_bus: float = None
"""TF coil bus length (m)"""


m_tf_bus: float = None
"""TF coil bus mass (kg)"""


tfckw: float = None
"""available DC power for charging the TF coils (kW)"""


tfcmw: float = None
"""Peak power per TF power supply (MW)"""


p_cp_resistive_mw: float = None
"""Peak resistive TF coil inboard leg power (MW)"""


p_tf_joints_resistive_mw: float = None
"""TF joints resistive power losses (MW)"""


tfcryoarea: float = None
"""surface area of toroidal shells covering TF coils (m2)"""


tficrn: float = None
"""TF coil half-width - inner dr_bore (m)"""


ind_tf_coil: float = None
"""TF coil inductance (H)"""


dx_tf_wp_insertion_gap: float = None
"""TF coil WP insertion gap (m)"""


p_tf_leg_resistive_mw: float = None
"""TF coil outboard leg resistive power (MW)"""


rho_cp: float = None
"""TF coil inboard leg resistivity [Ohm-m]. If `itart=0`, this variable is the
average resistivity over the whole magnet
"""


rho_tf_leg: float = None
"""Resistivity of a TF coil leg (Ohm-m)"""


rho_tf_bus: float = None
"""Resistivity of a TF coil bus (Ohm-m). Default values is for that of GLIDCOP AL-15 (C15715) at 293K"""


frhocp: float = None
"""Centrepost resistivity enhancement factor. For `itart=0`, this factor
is used for the whole magnet
"""


frholeg: float = None
"""Ouboard legs resistivity enhancement factor. Only used for `itart=1`."""


i_cp_joints: int = None
"""Switch for CP demoutable joints type
-= 0 : Clampled joints
-= 1 : Sliding joints
Default value (-1) choses :
Sliding joints for resistive magnets (i_tf_sup = 0, 2)
Clampled joints for superconducting magents (i_tf_sup = 1)
"""


rho_tf_joints: float = None
"""TF joints surfacic resistivity [ohm.m]. Feldmetal joints assumed."""


n_tf_joints_contact: int = None
"""Number of contact per turn"""


n_tf_joints: int = None
"""Number of joints
Ex: n_tf_joints = 2 for top and bottom CP joints
"""


th_joint_contact: float = None
"""TF sliding joints contact pad width [m]"""


p_tf_joints_resistive: float = None
"""Calculated TF joints resistive power losses [W]"""


len_tf_coil: float = None
"""TF coil circumference (m)"""


eff_tf_cryo: float = None
"""TF cryoplant efficiency (compared to pefect Carnot cycle).
Using -1 set the default value depending on magnet technology:
- i_tf_sup = 1 : SC magnet, eff_tf_cryo = 0.13 (ITER design)
- i_tf_sup = 2 : Cryo-aluminium, eff_tf_cryo = 0.4
"""


n_tf_coils: float = None
"""Number of TF coils (default = 50 for stellarators). Number of TF coils outer legs for ST"""


tfocrn: float = None
"""TF coil half-width - outer dr_bore (m)"""


tfsai: float = None
"""area of the inboard TF coil legs (m2)"""


tfsao: float = None
"""area of the outboard TF coil legs (m2)"""


tftmp: float = None
"""peak helium coolant temperature in TF coils and PF coils (K)"""


dx_tf_inboard_out_toroidal: float = None
"""TF coil toroidal thickness (m)"""


dx_tf_turn_insulation: float = None
"""conduit insulation thickness (m)"""


layer_ins: float = None
"""Additional insulation thickness between layers (m)"""


dr_tf_nose_case: float = None
"""inboard TF coil case outer (non-plasma side) thickness (m) (`iteration variable 57`)
(calculated for stellarators)
"""


dr_tf_wp_with_insulation: float = None
"""radial thickness of winding pack (m) (`iteration variable 140`) (issue #514)"""


dx_tf_turn_steel: float = None
"""TF coil turn steel conduit case thickness (m) (`iteration variable 58`)"""


dx_tf_wp_insulation: float = None
"""Thickness of the ground insulation layer surrounding (m)
- Superconductor TF (`i_tf_sup == 1`) : The TF coil Winding packs
- Resistive magnets (`i_tf_sup /= 1`) : The TF coil wedges
Rem : Thickness calculated for stellarators.
"""


temp_tf_superconductor_margin_min: float = None
"""minimum allowable temperature margin : TF coils (K)"""


temp_cs_superconductor_margin_min: float = None
"""minimum allowable temperature margin : CS (K)"""


tmargmin: float = None
"""minimum allowable temperature margin : TFC AND CS (K)"""


temp_margin: float = None
"""temperature margin (K)"""


temp_tf_superconductor_margin: float = None
"""TF coil superconductor temperature margin (K)"""


temp_tf_conductor_quench_max: float = None
"""maximum temp during a quench for protection (K)"""


temp_croco_quench_max: float = None
"""CroCo strand: maximum permitted temp during a quench (K)"""


temp_croco_quench: float = None
"""CroCo strand: Actual temp reached during a quench (K)"""


temp_tf_cryo: float = None
"""coil temperature for cryogenic plant power calculation (K)"""


n_tf_coil_turns: float = None
"""number of turns per TF coil"""


v_tf_coil_dump_quench_max_kv: float = None
"""max voltage across TF coil during quench (kV) (`iteration variable 52`)"""


vforce: float = None
"""vertical tension on inboard leg/coil (N)"""


f_vforce_inboard: float = None
"""Fraction of the total vertical force taken by the TF inboard leg tension
Not used for resistive `itart=1` (sliding joints)
"""


vforce_outboard: float = None
"""Vertical tension on outboard leg/coil (N)"""


f_a_tf_turn_cable_space_extra_void: float = None
"""coolant fraction of TFC 'cable' (`i_tf_sup=1`), or of TFC leg (`i_tf_ssup=0`)"""


voltfleg: float = None
"""volume of each TF coil outboard leg (m3)"""


vtfkv: float = None
"""TF coil voltage for resistive coil including bus (kV)"""


v_tf_coil_dump_quench_kv: float = None
"""voltage across a TF coil during quench (kV)"""


m_tf_coil_case: float = None
"""mass per coil of external case (kg)"""


m_tf_coil_conductor: float = None
"""TF coil conductor mass per coil (kg/coil).
For `itart=1`, coil is return limb plus centrepost/n_tf_coils
"""


m_tf_coil_copper: float = None
"""copper mass in TF coil conductor (kg/coil).
For `itart=1`, coil is return limb plus centrepost/n_tf_coils
"""


whtconal: float = None
"""Aluminium mass in TF coil conductor (kg/coil).
For `itart=1`, coil is return limb plus centrepost/n_tf_coils
"""


m_tf_coil_wp_turn_insulation: float = None
"""conduit insulation mass in TF coil conductor (kg/coil)"""


m_tf_coil_superconductor: float = None
"""superconductor mass in TF coil cable (kg/coil)"""


m_tf_wp_steel_conduit: float = None
"""steel conduit mass in TF coil conductor (kg/coil)"""


m_tf_coil_wp_insulation: float = None
"""mass of ground-wall insulation layer per coil (kg/coil)"""


m_tf_coils_total: float = None
"""total mass of the TF coils (kg)"""


dx_tf_wp_primary_toroidal: float = None
"""width of first step of winding pack (m)"""


dx_tf_wp_secondary_toroidal: float = None
"""width of second step of winding pack (m)"""


# Superconducting TF coil shape parameters;
# the TF inner surface top half is approximated by four circular arcs.
# Arc 1 goes through points 1 and 2 on the inner surface. Arc 2
# goes through points 2 and 3, etc.

dthet: list[float] = None
"""angle of arc i (rad)"""


radctf: list[float] = None
"""radius of arc i (m)"""


r_tf_arc: list[float] = None
"""x location of arc point i on surface (m)"""


xctfc: list[float] = None
"""x location of arc centre i (m)"""


z_tf_arc: list[float] = None
"""y location of arc point i on surface (m)"""


yctfc: list[float] = None
"""y location of arc centre i (m)"""


# New TF shape:  Horizontal and vertical radii of inside edge of TF coil
# Arcs are numbered clockwise:
# 1=upper inboard, 2=upper outboard, 3=lower ouboard, 4=lower inboard


tfa: list[float] = None
"""Horizontal radius of inside edge of TF coil (m)"""


tfb: list[float] = None
"""Vertical radius of inside edge of TF coil (m)"""


# Quantities relating to the spherical tokamak model (itart=1)
# (and in some cases, also to resistive TF coils, i_tf_sup=0):


drtop: float = None
"""centrepost taper maximum radius adjustment (m)"""


dztop: float = None
"""centrepost taper height adjustment (m)"""


etapump: float = None
"""centrepost coolant pump efficiency"""


fcoolcp: float = None
"""coolant fraction of TF coil inboard legs (`iteration variable 23`)"""


f_a_tf_cool_outboard: float = None
"""coolant fraction of TF coil outboard legs"""


a_cp_cool: float = None
"""Centrepost cooling area toroidal cross-section (constant over the whole CP)"""


ncool: float = None
"""number of centrepost coolant tubes"""


p_cp_coolant_pump_elec: float = None
"""centrepost coolant pump power (W)"""


p_cp_resistive: float = None
"""resistive power in the centrepost (itart=1) [W].
If `itart=0`, this variable is the ressitive power on the whole magnet
"""


p_tf_leg_resistive: float = None
"""Summed resistive power in the TF coil legs [W]. Remain 0 if `itart=0`."""


temp_cp_max: float = None
"""maximum peak centrepost temperature (K) (`constraint equation 44`)"""


rcool: float = None
"""average radius of coolant channel (m) (`iteration variable 69`)"""


tcoolin: float = None
"""centrepost coolant inlet temperature (K)"""


dtiocool: float = None
"""inlet / outlet TF coil coolant temperature rise (K)"""


temp_cp_average: float = None
"""Average temperature of centrepost called CP (K). Only used for resistive coils
to compute the resisitive heating. Must be an iteration variable for
ST (`itart=1`) (`iteration variable 20`)
"""


tcpav2: float = None
"""Computed centrepost average temperature (K) (for consistency)"""


temp_tf_legs_outboard: float = None
"""Average temperature of the TF outboard legs [K]. If `temp_tf_legs_outboard=-1.0`, the ouboard
legs and CP temperatures are the same. Fixed for now, should use a contraints eq like temp_cp_average
"""


temp_cp_peak: float = None
"""peak centrepost temperature (K)"""


vcool: float = None
"""inlet centrepost coolant flow speed at midplane (m/s) (`iteration variable 70`)"""


vol_cond_cp: float = None
"""Exact conductor volume in the centrepost (m3)"""


whtcp: float = None
"""mass of TF coil inboard legs (kg)"""


whttflgs: float = None
"""mass of the TF coil legs (kg)"""


cryo_cool_req: float = None
"""Cryo cooling requirement at helium temp 4.5K (kW)"""


theta1_coil: float = None
"""The angle of the outboard arc forming the TF coil current center line [deg]"""


theta1_vv: float = None
"""The angle of the outboard arc forming the Vacuum Vessel current center line [deg]"""


max_vv_stress: float = None
"""The allowable peak maximum shear stress in the vacuum vessel due to quench and fast discharge of the TF coils [Pa]"""


t_tf_quench_detection: float = None
"""TF coil quench detection time (s). Only used for TF coil quench protection."""


rrr_tf_cu: float = None
"""TF coil copper residual-resistance-ratio (RRR). Only used for quench protection."""


def init_tfcoil_variables():
    global \
        a_tf_coil_inboard_case, \
        a_tf_coil_outboard_case, \
        a_tf_turn_steel, \
        a_tf_wp_conductor, \
        a_res_tf_coil_conductor, \
        a_tf_turn_cable_space_no_void, \
        a_tf_turn_insulation, \
        a_tf_coil_wp_turn_insulation, \
        sig_tf_case_max, \
        sig_tf_wp_max, \
        a_tf_leg_outboard, \
        a_tf_wp_steel, \
        a_tf_wp_extra_void, \
        a_tf_wp_coolant_channels, \
        bcritsc, \
        b_tf_inboard_peak_symmetric, \
        b_tf_inboard_peak_with_ripple, \
        casestr, \
        dr_tf_plasma_case, \
        f_dr_tf_plasma_case, \
        i_f_dr_tf_plasma_case, \
        dx_tf_side_case_min, \
        dx_tf_side_case_peak, \
        casths_fraction, \
        tfc_sidewall_is_fraction, \
        t_conductor, \
        dx_tf_turn_general, \
        i_dx_tf_turn_general_input, \
        f_t_turn_tf, \
        t_turn_tf_max, \
        dx_tf_turn_cable_space_general, \
        i_dx_tf_turn_cable_space_general_input, \
        acs, \
        cdtfleg, \
        cforce, \
        cplen, \
        c_tf_turn, \
        c_tf_turn_max, \
        den_tf_coil_case, \
        dcond, \
        den_tf_wp_turn_insulation, \
        dia_tf_turn_coolant_channel, \
        e_tf_magnetic_stored_total_gj, \
        e_tf_coil_magnetic_stored, \
        b_crit_upper_nbti, \
        t_crit_nbti, \
        max_force_density, \
        f_a_tf_turn_cable_copper, \
        fhts, \
        insstrain, \
        i_tf_stress_model, \
        i_tf_tresca, \
        i_tf_wp_geom, \
        i_tf_case_geom, \
        i_tf_turns_integer, \
        i_tf_sc_mat, \
        i_tf_sup, \
        i_tf_shape, \
        i_tf_cond_eyoung_axial, \
        i_tf_cond_eyoung_trans, \
        n_tf_wp_pancakes, \
        n_tf_wp_layers, \
        n_rad_per_layer, \
        i_tf_bucking, \
        n_tf_graded_layers, \
        n_tf_stress_layers, \
        n_tf_wp_stress_layers, \
        j_tf_bus, \
        j_crit_str_tf, \
        j_crit_str_0, \
        j_tf_wp_critical, \
        j_tf_wp_quench_heat_max, \
        j_tf_wp, \
        oacdcp, \
        eyoung_ins, \
        eyoung_steel, \
        eyoung_cond_axial, \
        eyoung_cond_trans, \
        eyoung_res_tf_buck, \
        eyoung_copper, \
        eyoung_al, \
        poisson_steel, \
        poisson_copper, \
        poisson_al, \
        poisson_ins, \
        poisson_cond_axial, \
        poisson_cond_trans, \
        r_b_tf_inboard_peak, \
        res_tf_leg, \
        toroidalgap, \
        ftoroidalgap, \
        ripple_b_tf_plasma_edge_max, \
        ripple_b_tf_plasma_edge, \
        c_tf_total, \
        radial_array, \
        sig_tf_r, \
        sig_tf_t, \
        deflect, \
        sig_tf_z, \
        sig_tf_vmises, \
        s_shear_tf, \
        sig_tf_cs_bucked, \
        sig_tf_case, \
        sig_tf_wp, \
        str_cs_con_res, \
        str_pf_con_res, \
        str_tf_con_res, \
        str_wp, \
        str_wp_max, \
        i_str_wp, \
        quench_model, \
        time1, \
        tcritsc, \
        t_tf_superconductor_quench, \
        a_tf_inboard_total, \
        len_tf_bus, \
        m_tf_bus, \
        tfckw, \
        tfcmw, \
        p_cp_resistive_mw, \
        p_tf_joints_resistive_mw, \
        tfcryoarea, \
        tficrn, \
        ind_tf_coil, \
        dx_tf_wp_insertion_gap, \
        p_tf_leg_resistive_mw, \
        rho_cp, \
        rho_tf_leg, \
        rho_tf_bus, \
        frhocp, \
        frholeg, \
        i_cp_joints, \
        rho_tf_joints, \
        n_tf_joints_contact, \
        n_tf_joints, \
        th_joint_contact, \
        p_tf_joints_resistive, \
        len_tf_coil, \
        eff_tf_cryo, \
        n_tf_coils, \
        tfocrn, \
        tfsai, \
        tfsao, \
        tftmp, \
        dx_tf_inboard_out_toroidal, \
        dx_tf_turn_insulation, \
        layer_ins, \
        dr_tf_nose_case, \
        dr_tf_wp_with_insulation, \
        dx_tf_turn_steel, \
        dx_tf_wp_insulation, \
        temp_tf_superconductor_margin_min, \
        temp_cs_superconductor_margin_min, \
        tmargmin, \
        temp_margin, \
        temp_tf_superconductor_margin, \
        temp_tf_conductor_quench_max, \
        temp_croco_quench_max, \
        temp_croco_quench, \
        temp_tf_cryo, \
        n_tf_coil_turns, \
        v_tf_coil_dump_quench_max_kv, \
        vforce, \
        f_vforce_inboard, \
        vforce_outboard, \
        f_a_tf_turn_cable_space_extra_void, \
        voltfleg, \
        vtfkv, \
        v_tf_coil_dump_quench_kv, \
        m_tf_coil_case, \
        m_tf_coil_conductor, \
        m_tf_coil_copper, \
        whtconal, \
        m_tf_coil_wp_turn_insulation, \
        m_tf_coil_superconductor, \
        m_tf_wp_steel_conduit, \
        m_tf_coil_wp_insulation, \
        m_tf_coils_total, \
        dx_tf_wp_primary_toroidal, \
        dx_tf_wp_secondary_toroidal, \
        dthet, \
        radctf, \
        r_tf_arc, \
        xctfc, \
        z_tf_arc, \
        yctfc, \
        tfa, \
        tfb, \
        drtop, \
        dztop, \
        etapump, \
        fcoolcp, \
        f_a_tf_cool_outboard, \
        a_cp_cool, \
        ncool, \
        p_cp_coolant_pump_elec, \
        p_cp_resistive, \
        p_tf_leg_resistive, \
        temp_cp_max, \
        rcool, \
        tcoolin, \
        dtiocool, \
        temp_cp_average, \
        tcpav2, \
        temp_tf_legs_outboard, \
        temp_cp_peak, \
        vcool, \
        vol_cond_cp, \
        whtcp, \
        whttflgs, \
        cryo_cool_req, \
        theta1_coil, \
        theta1_vv, \
        max_vv_stress, \
        t_tf_quench_detection, \
        rrr_tf_cu

    a_tf_coil_inboard_case = 0.0
    a_tf_coil_outboard_case = 0.0
    a_tf_turn_steel = 0.0
    a_tf_wp_conductor = 0.0
    a_res_tf_coil_conductor = 0.0
    a_tf_turn_cable_space_no_void = 0.0
    a_tf_turn_insulation = 0.0
    a_tf_coil_wp_turn_insulation = 0.0
    sig_tf_case_max = 6.0e8
    sig_tf_wp_max = 6.0e8
    a_tf_leg_outboard = 0.0
    a_tf_wp_steel = 0.0
    a_tf_wp_extra_void = 0.0
    a_tf_wp_coolant_channels = 0.0
    bcritsc = 24.0
    b_tf_inboard_peak_symmetric = 0.0
    b_tf_inboard_peak_with_ripple = 0.0
    casestr = 0.0
    dr_tf_plasma_case = 0.0
    f_dr_tf_plasma_case = 0.05
    i_f_dr_tf_plasma_case = False
    dx_tf_side_case_min = 0.0
    dx_tf_side_case_peak = 0.0
    casths_fraction = 0.06
    t_conductor = 0.0
    dx_tf_turn_cable_space_general = 0.0
    i_dx_tf_turn_cable_space_general_input = False
    dx_tf_turn_general = 0.0
    i_dx_tf_turn_general_input = False
    f_t_turn_tf = 1.0
    t_turn_tf_max = 0.05
    acs = 0.0
    cdtfleg = 0.0
    cforce = 0.0
    cplen = 0.0
    c_tf_turn = 7.0e4
    c_tf_turn_max = 9.0e4
    den_tf_coil_case = 8000.0
    dcond = np.array([
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
    den_tf_wp_turn_insulation = 1800.0
    dia_tf_turn_coolant_channel = 0.005
    e_tf_magnetic_stored_total_gj = 0.0
    e_tf_coil_magnetic_stored = 0.0
    b_crit_upper_nbti = 14.86
    t_crit_nbti = 9.04
    max_force_density = 0.0
    f_a_tf_turn_cable_copper = 0.69
    fhts = 0.5
    insstrain = 0.0
    i_tf_stress_model = 1
    i_tf_tresca = 0
    i_tf_wp_geom = -1
    i_tf_case_geom = 0
    i_tf_turns_integer = 0
    i_tf_sc_mat = 1
    i_tf_sup = 1
    i_tf_shape = 0
    i_tf_cond_eyoung_axial = 0
    i_tf_cond_eyoung_trans = 1
    n_tf_wp_pancakes = 10
    n_tf_wp_layers = 20
    n_rad_per_layer = 100
    i_tf_bucking = -1
    n_tf_graded_layers = 1
    n_tf_stress_layers = 0
    n_tf_wp_stress_layers = 5
    j_tf_bus = 1.25e6
    j_crit_str_tf = 0.0
    j_crit_str_0 = np.array([
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
    j_tf_wp_critical = 0.0
    j_tf_wp_quench_heat_max = 0.0
    j_tf_wp = 0.0
    oacdcp = 0.0
    eyoung_ins = 1.0e8
    eyoung_steel = 2.05e11
    eyoung_cond_axial = 6.6e8
    eyoung_cond_trans = 0.0
    eyoung_res_tf_buck = 150.0e9
    eyoung_copper = 117.0e9
    eyoung_al = 69.0e9
    poisson_steel = 0.3
    poisson_copper = 0.35
    poisson_al = 0.35
    poisson_ins = 0.34
    poisson_cond_axial = 0.3
    poisson_cond_trans = 0.3
    r_b_tf_inboard_peak = 0.0
    res_tf_leg = 0.0
    toroidalgap = 1.0  # [m]
    ftoroidalgap = 1.0
    ripple_b_tf_plasma_edge_max = 1.0
    ripple_b_tf_plasma_edge = 0.0
    c_tf_total = 0.0
    radial_array = np.zeros(N_RADIAL_ARRAY)
    sig_tf_r = np.zeros(N_RADIAL_ARRAY)
    sig_tf_t = np.zeros(N_RADIAL_ARRAY)
    deflect = np.zeros(N_RADIAL_ARRAY)
    sig_tf_z = 0.0
    sig_tf_vmises = np.zeros(N_RADIAL_ARRAY)
    s_shear_tf = np.zeros(N_RADIAL_ARRAY)
    sig_tf_cs_bucked = 0.0
    sig_tf_case = 0.0
    sig_tf_wp = 0.0
    str_cs_con_res = -0.005
    str_pf_con_res = -0.005
    str_tf_con_res = -0.005
    str_wp = 0.0
    str_wp_max = 0.7e-2
    i_str_wp = 1
    quench_model = "exponential"
    time1 = 0
    tcritsc = 16.0
    t_tf_superconductor_quench = 10.0
    a_tf_inboard_total = 0.0
    len_tf_bus = 300.0
    m_tf_bus = 0.0
    tfckw = 0.0
    tfcmw = 0.0
    p_cp_resistive_mw = 0.0
    p_tf_joints_resistive_mw = 0.0
    tfcryoarea = 0.0
    tficrn = 0.0
    ind_tf_coil = 0.0
    dx_tf_wp_insertion_gap = 0.01
    p_tf_leg_resistive_mw = 0.0
    rho_cp = 0.0
    rho_tf_leg = 0.0
    rho_tf_bus = 1.86e-8
    frhocp = 1.0
    frholeg = 1.0
    rho_tf_joints = 2.5e-10
    n_tf_joints_contact = 6
    n_tf_joints = 4
    th_joint_contact = 0.03
    p_tf_joints_resistive = 0.0
    len_tf_coil = 0.0
    eff_tf_cryo = -1.0
    n_tf_coils = 16.0
    tfocrn = 0.0
    tfsai = 0.0
    tfsao = 0.0
    tftmp = 4.5
    dx_tf_inboard_out_toroidal = 1.0
    dx_tf_turn_insulation = 8e-4
    layer_ins = 0.0
    dr_tf_nose_case = 0.3
    dr_tf_wp_with_insulation = 0.0
    dx_tf_turn_steel = 8e-3
    dx_tf_wp_insulation = 0.018
    temp_tf_superconductor_margin_min = 0.0
    temp_cs_superconductor_margin_min = 0.0
    tmargmin = 0.0
    temp_margin = 0.0
    temp_tf_superconductor_margin = 0.0
    temp_tf_conductor_quench_max = 150.0
    temp_croco_quench_max = 200.0
    temp_croco_quench = 0.0
    temp_tf_cryo = 4.5
    n_tf_coil_turns = 0.0
    v_tf_coil_dump_quench_max_kv = 20.0
    vforce = 0.0
    f_vforce_inboard = 0.5
    vforce_outboard = 0.0
    f_a_tf_turn_cable_space_extra_void = 0.4
    voltfleg = 0.0
    vtfkv = 0.0
    v_tf_coil_dump_quench_kv = 0.0
    m_tf_coil_case = 0.0
    m_tf_coil_conductor = 0.0
    m_tf_coil_copper = 0.0
    whtconal = 0.0
    m_tf_coil_wp_turn_insulation = 0.0
    m_tf_coil_superconductor = 0.0
    m_tf_wp_steel_conduit = 0.0
    m_tf_coil_wp_insulation = 0.0
    m_tf_coils_total = 0.0
    dx_tf_wp_primary_toroidal = 0.0
    dx_tf_wp_secondary_toroidal = 0.0
    dthet = np.zeros(4)
    radctf = np.zeros(4)
    r_tf_arc = np.zeros(5)
    xctfc = np.zeros(4)
    z_tf_arc = np.zeros(5)
    yctfc = np.zeros(4)
    tfa = np.zeros(4)
    tfb = np.zeros(4)
    drtop = 0.0
    dztop = 0.0
    etapump = 0.8
    fcoolcp = 0.3
    f_a_tf_cool_outboard = 0.2
    a_cp_cool = 0.0
    ncool = 0.0
    p_cp_coolant_pump_elec = 0.0
    p_cp_resistive = 0.0
    p_tf_leg_resistive = 0.0
    temp_cp_max = 473.15  # 200 C
    rcool = 0.005
    tcoolin = 313.15  # 40 C
    dtiocool = 0.0
    temp_cp_average = 373.15  # 100 C
    tcpav2 = 0.0
    temp_tf_legs_outboard = -1.0
    temp_cp_peak = 0.0
    vcool = 20.0
    vol_cond_cp = 0.0
    whtcp = 0.0
    whttflgs = 0.0
    tfc_sidewall_is_fraction = False
    i_cp_joints = -1
    cryo_cool_req = 0.0
    theta1_coil = 45.0
    theta1_vv = 1.0  # 1 Deg
    max_vv_stress = 143.0e6
    rrr_tf_cu = 100.0
    t_tf_quench_detection = 3.0
