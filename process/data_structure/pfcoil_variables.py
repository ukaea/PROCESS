import numpy as np

N_PF_GROUPS_MAX = 10
"""maximum number of groups of PF coils"""

N_PF_COILS_IN_GROUP_MAX = 2
"""maximum number of PF coils in a given group"""

NPTSMX = 32
"""maximum number of points across the midplane of the plasma at which the field from
the PF coils is fixed
"""

NFIXMX = 64
"""maximum number of fixed current PF coils"""

NGC = N_PF_GROUPS_MAX * N_PF_COILS_IN_GROUP_MAX
"""maximum total number of coils across all groups"""

NGC2 = NGC + 2
"""new variable to include 2 additional circuits: plasma and central solenoid"""

# PF coil module variables
nef: int = None

nfxf: int = None

ricpf: float = None

ssq0: float = None

sig_axial: float = None

sig_hoop: float = None

axial_force: float = None

r_pf_cs_current_filaments: list[float] = None
"""array of radial positions of current filaments in central solenoid"""

z_pf_cs_current_filaments: list[float] = None
"""array of vertical positions of current filaments in central solenoid"""

c_pf_cs_current_filaments: list[float] = None
"""array of current in filaments in central solenoid"""

xind: list[float] = None

r_pf_coil_middle_group_array: list[float] = None
"""2D array of PF coil middle radii, indexed by group and coil in group"""

z_pf_coil_middle_group_array: list[float] = None
"""2D array of PF coil middle heights, indexed by group and coil in group"""

ccls: list[float] = None

ccl0: list[float] = None

bpf2: list[float] = None

vsdum: list[float] = None

first_call: bool = None

cslimit: bool = None

# PF coil variables

alfapf: float = None
"""smoothing parameter used in PF coil current calculation at the beginning of pulse (BoP)"""


alstroh: float = None
"""allowable hoop stress in Central Solenoid structural material (Pa)"""


i_cs_stress: int = None
"""Switch for CS stress calculation:
- =0 Hoop stress only
- =1 Hoop + Axial stress
"""


a_cs_poloidal: float = None
"""Central solenoid vertical cross-sectional area (m2)"""


a_cs_turn: float = None
"""Central solenoid (OH) trun cross-sectional area (m2)"""


awpoh: float = None
"""central solenoid conductor+void area with area of steel subtracted (m2)"""


b_cs_peak_flat_top_end: float = None
"""maximum field in central solenoid at end of flat-top (EoF) (T)"""


b_cs_peak_pulse_start: float = None
"""maximum field in central solenoid at beginning of pulse (T)"""


b_pf_coil_peak: list[float] = None
"""peak field at coil i (T)"""


ccl0_ma: list[float] = None
"""PF group current array, flux-swing cancellation current (MA)
Input if i_pf_current=0, computed otherwise
"""


ccls_ma: list[float] = None
"""PF group current array, equilibrium current (MA)
Input if i_pf_current=0, computed otherwise
"""


j_cs_pulse_start: float = None
"""Central solenoid overall current density at beginning of pulse (A/m2)"""


j_cs_flat_top_end: float = None
"""Central solenoid overall current density at end of flat-top (A/m2) (`iteration variable 37`) (`sweep variable 62`)"""


c_pf_coil_turn: list[float] = None
"""current per turn in coil i at time j (A)"""


c_pf_coil_turn_peak_input: list[float] = None
"""peak current per turn input for PF coil i (A)"""


c_pf_cs_coil_pulse_start_ma: list[float] = None
"""PF coil current array, at beginning of pulse (MA)
Indexed by coil number, not group number
"""


c_pf_cs_coil_flat_top_ma: list[float] = None
"""PF coil current array, at flat top (MA)
Indexed by coil number, not group number
"""


c_pf_cs_coil_pulse_end_ma: list[float] = None
"""PF coil current array, at end of pulse (MA)
Indexed by coil number, not group number
"""


etapsu: float = None
"""Efficiency of transfer of PF stored energy into or out of storage."""


f_j_cs_start_end_flat_top: float = None
"""ratio of central solenoid overall current density at beginning of flat-top / end of flat-top"""


f_j_cs_start_pulse_end_flat_top: float = None
"""ratio of central solenoid overall current density at beginning of pulse / end of flat-top
(`iteration variable 41`)
"""


fcuohsu: float = None
"""copper fraction of strand in central solenoid"""


fcupfsu: float = None
"""copper fraction of cable conductor (PF coils)"""


fvs_cs_pf_total_ramp: float = None
"""F-value for `constraint equation 51`"""


i_pf_location: list[int] = None
"""Switch for location of PF coil group i:
- =1 PF coil on top of central solenoid (flux ramp only)
- =2 PF coil on top of TF coil (flux ramp only)
- =3 PF coil outside of TF coil (equilibrium coil)
- =4 PF coil, general location (equilibrium coil)
"""


i_pf_conductor: int = None
"""switch for PF & CS coil conductor type:
- =0 superconducting PF coils
- =1 resistive PF coils
"""


itr_sum: float = None
"""total sum of I x turns x radius for all PF coils and CS (Am)"""


i_cs_superconductor: int = None
"""switch for superconductor material in central solenoid:
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


i_pf_superconductor: int = None
"""switch for superconductor material in PF coils:
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


j_crit_str_cs: float = None
"""superconductor strand critical current density under operating
conditions in central solenoid (A/m2). Necessary for the cost calculation in $/kA m
"""


j_crit_str_pf: float = None
"""superconductor strand critical current density under operating
conditions in PF coils (A/m2). Necessary for the cost calculation in $/kA m
"""


i_pf_current: int = None
"""Switch for controlling the current of the PF coils:
- =0 Input via the variables c_pf_cs_coil_pulse_start_ma, c_pf_cs_coil_flat_top_ma, c_pf_cs_coil_pulse_end_ma
- =1 SVD targets zero field across midplane (flux swing
coils) and the correct vertical field at the plasma
center (equilibrium coils)
"""


i_r_pf_outside_tf_placement: int = None
"""Switch for the placement of Location 3 (outboard) PF coils
- =0 (Default) Outboard PF coils follow TF shape
in an ellipsoidal winding surface
- =1 Outboard PF coils all have same radius, cylindrical
winding surface
"""


j_cs_conductor_critical_pulse_start: float = None
"""central solenoid superconductor critical current density (A/m2) at beginning-of-pulse"""


j_cs_conductor_critical_flat_top_end: float = None
"""central solenoid superconductor critical current density (A/m2) at end-of-flattop"""


jcableoh_bop: float = None
"""central solenoid cable critical current density (A/m2) at beginning-of-pulse"""


jcableoh_eof: float = None
"""central solenoid cable critical current density (A/m2) at end-of-flattop"""


n_pf_cs_plasma_circuits: int = None
"""number of PF circuits (including central solenoid and plasma)"""


n_pf_coils_in_group: list[int] = None
"""number of PF coils in group j"""


n_cs_current_filaments: int = None
"""number of filaments the top and bottom of the central solenoid should be broken
into during scaling (5 - 10 is good)
"""


n_pf_coil_groups: int = None
"""number of groups of PF coils. Symmetric coil pairs should all be in the same group"""


n_cs_pf_coils: int = None
"""number of PF coils (excluding the central solenoid) + 1"""


f_z_cs_tf_internal: float = None
"""Central solenoid height / TF coil internal height"""


f_a_cs_turn_steel: float = None
"""Fraction of CS turn poloidal area that is steel (`iteration variable 122`)"""


pf_current_safety_factor: float = None
"""Ratio of permissible PF coil conductor current density to critical conductor
current density based on short-sample DC measurements
"""


pfcaseth: list[float] = None
"""steel case thickness for PF coil i (m)"""


rho_pf_coil: float = None
"""PF coil resistivity (if i_pf_conductor=1) (Ohm-m)"""


rhopfbus: float = None
"""Resistivity of CS and PF coil bus bars (irrespective of
whether the coils themselves are superconducting or resistive) (Ohm-m)
"""


m_pf_coil_max: float = None
"""mass of heaviest PF coil (tonnes)"""


r_pf_coil_outer_max: float = None
"""radius of largest PF coil (m)"""


p_pf_electric_supplies_mw: float = None
"""Total mean wall plug power dissipated in PFC and CS power supplies (MW) (issue #713)"""


p_cs_resistive_flat_top: float = None
"""central solenoid resistive power during flattop (W)"""


p_pf_coil_resistive_total_flat_top: float = None
"""total PF coil resistive losses during flattop (W)"""


r_pf_coil_inner: list[float] = None
"""inner radius of coil i (m)"""


r_pf_coil_outer: list[float] = None
"""outer radius of coil i (m)"""


c_pf_cs_coils_peak_ma: list[float] = None
"""peak current in coil i (MA-turns)"""


j_pf_coil_wp_peak: list[float] = None
"""average winding pack current density of PF coil i (A/m2) at time of peak
current in that coil (calculated for `i_pf_location=1` coils)
"""


j_cs_critical_flat_top_end: float = None
"""allowable central solenoid current density at end of flat-top (A/m2)"""


j_cs_critical_pulse_start: float = None
"""allowable central solenoid current density at beginning of pulse (A/m2)"""


j_pf_wp_critical: list[float] = None
"""allowable winding pack current density of PF coil i (A/m2)"""


r_cs_middle: float = None
"""radius to the centre of the central solenoid (m)"""

dz_cs_full: float = None
"""Full height of the central solenoid (m)"""


dr_pf_tf_outboard_out_offset: float = None
"""radial distance (m) from outboard TF coil leg to centre of `i_pf_location=3` PF coils"""


r_pf_coil_middle: list[float] = None
"""radius of PF coil i (m)"""


dr_pf_cs_middle_offset: float = None
"""offset (m) of radial position of `i_pf_location=1` PF coils from being directly above
the central solenoid
"""


rpf2: float = None
"""offset (m) of radial position of `i_pf_location=2` PF coils from being at
rmajor (offset = rpf2*triang*rminor)
"""


rref: list[float] = None
"""PF coil radial positioning adjuster:
- for groups j with i_pf_location(j) = 1; rref(j) is ignored
- for groups j with i_pf_location(j) = 2; rref(j) is ignored
- for groups j with i_pf_location(j) = 3; rref(j) is ignored
- for groups j with i_pf_location(j) = 4; rref(j) is radius of
the coil in units of minor radii from the major radius
(r = rmajor + rref*rminor)
"""


s_shear_cs_peak: float = None
"""Maximum shear stress (Tresca criterion) coils/central solenoid [MPa]"""


sigpfcalw: float = None
"""maximum permissible tensile stress (MPa) in steel coil cases for superconducting
PF coils (`i_pf_conductor=0`)
"""


sigpfcf: float = None
"""fraction of JxB hoop force supported by steel case for superconducting PF coils (`i_pf_conductor=0`)"""


ind_pf_cs_plasma_mutual: list[float] = None
"""mutual inductance matrix (H)"""


temp_cs_superconductor_margin: float = None
"""Central solenoid temperature margin (K)"""


n_pf_coil_turns: list[float] = None
"""number of turns in PF coil i"""


f_a_pf_coil_void: list[float] = None
"""winding pack void fraction of PF coil i for coolant"""


f_a_cs_void: float = None
"""void fraction of central solenoid conductor for coolant"""


vs_cs_pf_total_burn: float = None
"""total flux swing available for burn (Wb)"""


vs_pf_coils_total_burn: float = None
"""flux swing from PF coils for burn (Wb)"""


vs_pf_coils_total_ramp: float = None
"""flux swing from PF coils for startup (Wb)"""


vs_pf_coils_total_pulse: float = None
"""total flux swing from PF coils (Wb)"""


vs_cs_total_pulse: float = None
"""total flux swing from the central solenoid (Wb)"""


vs_cs_burn: float = None
"""central solenoid flux swing for burn (Wb)"""


vs_cs_ramp: float = None
"""central solenoid flux swing for startup (Wb)"""


vs_cs_pf_total_ramp: float = None
"""total flux swing for startup (`constraint eqn 51` to enforce vs_cs_pf_total_ramp=vs_plasma_res_ramp+vs_plasma_ind_ramp) (Wb)"""


vs_cs_pf_total_pulse: float = None
"""total flux swing for pulse (Wb)"""


f_c_pf_cs_peak_time_array: list[float] = None
"""PF, CS coil current relative to peak current at time points 1 to 6"""


m_pf_coil_conductor_total: float = None
"""total mass of the PF coil conductor (kg)"""


m_pf_coil_structure_total: float = None
"""total mass of the PF coil structure (kg)"""


m_pf_coil_conductor: list[float] = None
"""conductor mass for PF coil i (kg)"""


m_pf_coil_structure: list[float] = None
"""structure mass for PF coil i (kg)"""


z_pf_coil_upper: list[float] = None
"""upper point of PF coil i (m)"""


z_pf_coil_lower: list[float] = None
"""lower point of PF coil i (m)"""


z_pf_coil_middle: list[float] = None
"""z (height) location of PF coil i (m)"""


zref: list[float] = None
"""PF coil vertical positioning adjuster:
- for groups j with i_pf_location(j) = 1; zref(j) is ignored
- for groups j with i_pf_location(j) = 2 AND itart=1 (only);
zref(j) is distance of centre of PF coil from inside
edge of TF coil (remember that PF coils for STs lie
within the TF coil)
- for groups j with i_pf_location(j) = 3; zref(j) = ratio of
height of coil group j to plasma minor radius</UL>
- for groups j with i_pf_location(j) = 4; zref(j) = ratio of
height of coil group j to plasma minor radius</UL>
"""


b_cs_limit_max: float = None
"""Central solenoid max field limit [T]"""


fb_cs_limit_max: float = None
"""F-value for CS mmax field (`cons. 79`, `itvar 149`)"""


f_dr_dz_cs_turn: float = None
"""Ratio of CS coil turn conduit length to depth"""


dr_cs_turn: float = None
"""Length of CS of CS coil turn conduit"""

dr_cs_full: float = None
"""Full radial thickness of the central solenoid (m)"""


dz_cs_turn: float = None
"""Depth/width of CS of CS coil turn conduit"""


radius_cs_turn_corners: float = None
"""Radius of curvature of CS coil turn corners (m)"""


radius_cs_turn_cable_space: float = None
"""Length of CS of CS coil turn conduit length"""


def init_pfcoil_module():
    global first_call
    global cslimit
    global nef
    global nfxf
    global ricpf
    global ssq0
    global sig_axial
    global sig_hoop
    global axial_force
    global r_pf_cs_current_filaments
    global z_pf_cs_current_filaments
    global c_pf_cs_current_filaments
    global xind
    global r_pf_coil_middle_group_array
    global z_pf_coil_middle_group_array
    global ccls
    global ccl0
    global bpf2
    global vsdum

    first_call = True
    cslimit = False
    nef = 0
    nfxf = 0
    ricpf = 0.0
    ssq0 = 0.0
    sig_axial = 0.0
    sig_hoop = 0.0
    axial_force = 0.0

    r_pf_cs_current_filaments = np.zeros(NFIXMX)
    z_pf_cs_current_filaments = np.zeros(NFIXMX)
    c_pf_cs_current_filaments = np.zeros(NFIXMX)
    xind = np.zeros(NFIXMX)
    r_pf_coil_middle_group_array = np.zeros((N_PF_GROUPS_MAX, N_PF_COILS_IN_GROUP_MAX))
    z_pf_coil_middle_group_array = np.zeros((N_PF_GROUPS_MAX, N_PF_COILS_IN_GROUP_MAX))
    ccls = np.zeros(N_PF_GROUPS_MAX)
    ccl0 = np.zeros(N_PF_GROUPS_MAX)
    bpf2 = np.zeros(NGC2)
    vsdum = np.zeros((NGC2, 3))


def init_pfcoil_variables():
    """Initialise the PF coil variables"""
    global alfapf
    global alstroh
    global i_cs_stress
    global a_cs_poloidal
    global a_cs_turn
    global awpoh
    global b_cs_peak_flat_top_end
    global b_cs_peak_pulse_start
    global b_pf_coil_peak
    global ccl0_ma
    global ccls_ma
    global j_cs_pulse_start
    global j_cs_flat_top_end
    global c_pf_coil_turn
    global c_pf_coil_turn_peak_input
    global c_pf_cs_coil_pulse_start_ma
    global c_pf_cs_coil_flat_top_ma
    global c_pf_cs_coil_pulse_end_ma
    global etapsu
    global f_j_cs_start_end_flat_top
    global f_j_cs_start_pulse_end_flat_top
    global fcuohsu
    global fcupfsu
    global fvs_cs_pf_total_ramp
    global i_pf_location
    global i_pf_conductor
    global itr_sum
    global i_cs_superconductor
    global i_pf_superconductor
    global j_crit_str_cs
    global j_crit_str_pf
    global i_pf_current
    global i_r_pf_outside_tf_placement
    global j_cs_conductor_critical_pulse_start
    global j_cs_conductor_critical_flat_top_end
    global jcableoh_bop
    global jcableoh_eof
    global n_pf_cs_plasma_circuits
    global n_pf_coils_in_group
    global n_cs_current_filaments
    global n_pf_coil_groups
    global n_cs_pf_coils
    global f_z_cs_tf_internal
    global f_a_cs_turn_steel
    global pf_current_safety_factor
    global pfcaseth
    global rho_pf_coil
    global rhopfbus
    global m_pf_coil_max
    global r_pf_coil_outer_max
    global p_pf_electric_supplies_mw
    global p_cs_resistive_flat_top
    global p_pf_coil_resistive_total_flat_top
    global r_pf_coil_inner
    global r_pf_coil_outer
    global c_pf_cs_coils_peak_ma
    global j_pf_coil_wp_peak
    global j_cs_critical_flat_top_end
    global j_cs_critical_pulse_start
    global j_pf_wp_critical
    global r_cs_middle
    global dz_cs_full
    global dr_pf_tf_outboard_out_offset
    global r_pf_coil_middle
    global dr_pf_cs_middle_offset
    global rpf2
    global rref
    global s_shear_cs_peak
    global sigpfcalw
    global sigpfcf
    global ind_pf_cs_plasma_mutual
    global temp_cs_superconductor_margin
    global n_pf_coil_turns
    global f_a_pf_coil_void
    global f_a_cs_void
    global vs_cs_pf_total_burn
    global vs_pf_coils_total_burn
    global vs_pf_coils_total_ramp
    global vs_pf_coils_total_pulse
    global vs_cs_total_pulse
    global vs_cs_burn
    global vs_cs_ramp
    global vs_cs_pf_total_ramp
    global vs_cs_pf_total_pulse
    global f_c_pf_cs_peak_time_array
    global m_pf_coil_conductor_total
    global m_pf_coil_structure_total
    global m_pf_coil_conductor
    global m_pf_coil_structure
    global z_pf_coil_upper
    global z_pf_coil_lower
    global z_pf_coil_middle
    global zref
    global b_cs_limit_max
    global fb_cs_limit_max
    global f_dr_dz_cs_turn
    global dr_cs_turn
    global dr_cs_full
    global dz_cs_turn
    global radius_cs_turn_corners
    global radius_cs_turn_cable_space

    alfapf = 5e-10
    alstroh = 4.0e8
    i_cs_stress = 0
    a_cs_poloidal = 0.0
    a_cs_turn = 0.0
    awpoh = 0.0
    b_cs_peak_flat_top_end = 0.0
    b_cs_peak_pulse_start = 0.0
    b_pf_coil_peak = np.zeros(NGC2)
    ccl0_ma = np.zeros(N_PF_GROUPS_MAX)
    ccls_ma = np.zeros(N_PF_GROUPS_MAX)
    j_cs_pulse_start = 0.0
    j_cs_flat_top_end = 1.85e7
    c_pf_coil_turn = np.zeros((NGC2, 6))
    c_pf_coil_turn_peak_input = np.full(NGC2, 4.0e4)
    c_pf_cs_coil_pulse_start_ma = np.zeros(NGC2)
    c_pf_cs_coil_flat_top_ma = np.zeros(NGC2)
    c_pf_cs_coil_pulse_end_ma = np.zeros(NGC2)
    etapsu = 0.9
    f_j_cs_start_end_flat_top = 0.0
    f_j_cs_start_pulse_end_flat_top = 0.9
    fcuohsu = 0.7
    fcupfsu = 0.69
    fvs_cs_pf_total_ramp = 1.0
    i_pf_location = np.array([2, 2, 3, 0, 0, 0, 0, 0, 0, 0])
    i_pf_conductor = 0
    itr_sum = 0.0
    i_cs_superconductor = 1
    i_pf_superconductor = 1
    j_crit_str_cs = 0.0
    j_crit_str_pf = 0.0
    i_pf_current = 1
    i_r_pf_outside_tf_placement = 0
    j_cs_conductor_critical_pulse_start = 0.0
    j_cs_conductor_critical_flat_top_end = 0.0
    jcableoh_bop = 0.0
    jcableoh_eof = 0.0
    n_pf_cs_plasma_circuits = 0
    n_pf_coils_in_group = np.array([1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    n_cs_current_filaments = 7
    n_pf_coil_groups = 3
    n_cs_pf_coils = 0
    f_z_cs_tf_internal = 0.71
    f_a_cs_turn_steel = 0.5
    pf_current_safety_factor = 1.0
    pfcaseth = np.zeros(NGC2)
    rho_pf_coil = 2.5e-8
    rhopfbus = 3.93e-8
    m_pf_coil_max = 0.0
    r_pf_coil_outer_max = 0.0
    p_pf_electric_supplies_mw = 0.0
    p_cs_resistive_flat_top = 0.0
    p_pf_coil_resistive_total_flat_top = 0.0
    r_pf_coil_inner = np.zeros(NGC2)
    r_pf_coil_outer = np.zeros(NGC2)
    c_pf_cs_coils_peak_ma = np.zeros(NGC2)
    j_pf_coil_wp_peak = np.full(NGC2, 3.0e7)
    j_cs_critical_flat_top_end = 0.0
    j_cs_critical_pulse_start = 0.0
    j_pf_wp_critical = np.zeros(NGC2)
    r_cs_middle = 0.0
    dz_cs_full = 0.0
    dr_pf_tf_outboard_out_offset = 1.5
    r_pf_coil_middle = np.zeros(NGC2)
    dr_pf_cs_middle_offset = 0.0
    rpf2 = -1.63
    rref = np.full(N_PF_GROUPS_MAX, 7.0)
    s_shear_cs_peak = 0.0
    sigpfcalw = 500.0
    sigpfcf = 1.0
    ind_pf_cs_plasma_mutual = np.zeros((NGC2, NGC2))
    temp_cs_superconductor_margin = 0.0
    n_pf_coil_turns = np.zeros(NGC2)
    f_a_pf_coil_void = np.full(NGC2, 0.3)
    f_a_cs_void = 0.3
    vs_cs_pf_total_burn = 0.0
    vs_pf_coils_total_burn = 0.0
    vs_pf_coils_total_ramp = 0.0
    vs_pf_coils_total_pulse = 0.0
    vs_cs_total_pulse = 0.0
    vs_cs_burn = 0.0
    vs_cs_ramp = 0.0
    vs_cs_pf_total_ramp = 0.0
    vs_cs_pf_total_pulse = 0.0
    f_c_pf_cs_peak_time_array = np.zeros((NGC2, 6))
    m_pf_coil_conductor_total = 0.0
    m_pf_coil_structure_total = 0.0
    m_pf_coil_conductor = np.zeros(NGC2)
    m_pf_coil_structure = np.zeros(NGC2)
    z_pf_coil_upper = np.zeros(NGC2)
    z_pf_coil_lower = np.zeros(NGC2)
    z_pf_coil_middle = np.zeros(NGC2)
    zref = np.array([3.6, 1.2, 2.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    b_cs_limit_max = 13.0
    fb_cs_limit_max = 1.0
    f_dr_dz_cs_turn = 70.0 / 22.0
    dr_cs_turn = 0.0
    dr_cs_full = 0.0
    dz_cs_turn = 0.0
    radius_cs_turn_cable_space = 0.0
    radius_cs_turn_corners = 3.0e-3
