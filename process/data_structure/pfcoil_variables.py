from dataclasses import dataclass, field

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


@dataclass
class PFCoilData:
    # PF coil module variables
    nef: int = 0

    nfxf: int = 0

    ricpf: float = 0.0

    ssq0: float = 0.0

    stress_z_cs_self_peak_midplane: float = 0.0
    """Peak axial stress (z) in central solenoid at midplane due to its own field (when at peak current) (Pa)"""

    sig_hoop: float = 0.0

    forc_z_cs_self_peak_midplane: float = 0.0
    """Axial force (z) on central solenoid at midplane due to its own field (when at peak current) (N)"""

    r_pf_cs_current_filaments: list[float] = field(
        default_factory=lambda: np.zeros(NFIXMX)
    )
    """array of radial positions of current filaments in central solenoid"""

    z_pf_cs_current_filaments: list[float] = field(
        default_factory=lambda: np.zeros(NFIXMX)
    )
    """array of vertical positions of current filaments in central solenoid"""

    c_pf_cs_current_filaments: list[float] = field(
        default_factory=lambda: np.zeros(NFIXMX)
    )
    """array of current in filaments in central solenoid"""

    xind: list[float] = field(default_factory=lambda: np.zeros(NFIXMX))

    r_pf_coil_middle_group_array: list[float] = field(
        default_factory=lambda: np.zeros((N_PF_GROUPS_MAX, N_PF_COILS_IN_GROUP_MAX))
    )
    """2D array of PF coil middle radii, indexed by group and coil in group"""

    z_pf_coil_middle_group_array: list[float] = field(
        default_factory=lambda: np.zeros((N_PF_GROUPS_MAX, N_PF_COILS_IN_GROUP_MAX))
    )
    """2D array of PF coil middle heights, indexed by group and coil in group"""

    ccls: list[float] = field(default_factory=lambda: np.zeros(N_PF_GROUPS_MAX))

    ccl0: list[float] = field(default_factory=lambda: np.zeros(N_PF_GROUPS_MAX))

    bpf2: list[float] = field(default_factory=lambda: np.zeros(NGC2))

    vsdum: list[float] = field(default_factory=lambda: np.zeros((NGC2, 3)))

    first_call: bool = True

    cslimit: bool = False

    # PF coil variables

    alfapf: float = 5e-10
    """smoothing parameter used in PF coil current calculation at the beginning of pulse (BoP)"""

    alstroh: float = 4.0e8
    """allowable hoop stress in Central Solenoid structural material (Pa)"""

    i_cs_stress: int = 1
    """Switch for CS stress calculation:
    - =0 Hoop stress only
    - =1 Hoop + Axial stress
    """

    a_cs_poloidal: float = 0.0
    """Central solenoid vertical cross-sectional area (m2)"""

    a_cs_turn: float = 0.0
    """Central solenoid (OH) trun cross-sectional area (m2)"""

    awpoh: float = 0.0
    """central solenoid conductor+void area with area of steel subtracted (m2)"""

    b_cs_peak_flat_top_end: float = 0.0
    """maximum field in central solenoid at end of flat-top (EoF) (T)"""

    b_cs_peak_pulse_start: float = 0.0
    """maximum field in central solenoid at beginning of pulse (T)"""

    b_pf_coil_peak: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """peak field at coil i (T)"""

    ccl0_ma: list[float] = field(default_factory=lambda: np.zeros(N_PF_GROUPS_MAX))
    """PF group current array, flux-swing cancellation current (MA)
    Input if i_pf_current=0, computed otherwise
    """

    ccls_ma: list[float] = field(default_factory=lambda: np.zeros(N_PF_GROUPS_MAX))
    """PF group current array, equilibrium current (MA)
    Input if i_pf_current=0, computed otherwise
    """

    j_cs_pulse_start: float = 0.0
    """Central solenoid overall current density at beginning of pulse (A/m2)"""

    j_cs_flat_top_end: float = 1.85e7
    """Central solenoid overall current density at end of flat-top (A/m2) (`iteration variable 37`) (`sweep variable 62`)"""

    c_pf_coil_turn: list[float] = field(default_factory=lambda: np.zeros((NGC2, 6)))
    """current per turn in coil i at time j (A)"""

    c_pf_coil_turn_peak_input: list[float] = field(
        default_factory=lambda: np.full(NGC2, 4.0e4)
    )
    """peak current per turn input for PF coil i (A)"""

    c_pf_cs_coil_pulse_start_ma: list[float] = field(
        default_factory=lambda: np.zeros(NGC2)
    )
    """PF coil current array, at beginning of pulse (MA)
    Indexed by coil number, not group number
    """

    c_pf_cs_coil_flat_top_ma: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """PF coil current array, at flat top (MA)
    Indexed by coil number, not group number
    """

    c_pf_cs_coil_pulse_end_ma: list[float] = field(
        default_factory=lambda: np.zeros(NGC2)
    )
    """PF coil current array, at end of pulse (MA)
    Indexed by coil number, not group number
    """

    etapsu: float = 0.9
    """Efficiency of transfer of PF stored energy into or out of storage."""

    f_j_cs_start_end_flat_top: float = 0.0
    """ratio of central solenoid overall current density at beginning of flat-top / end of flat-top"""

    f_j_cs_start_pulse_end_flat_top: float = 0.9
    """ratio of central solenoid overall current density at beginning of pulse / end of flat-top
    (`iteration variable 41`)
    """

    fcuohsu: float = 0.7
    """copper fraction of strand in central solenoid"""

    fcupfsu: float = 0.69
    """copper fraction of cable conductor (PF coils)"""

    i_pf_location: list[int] = field(
        default_factory=lambda: np.array([2, 2, 3, 0, 0, 0, 0, 0, 0, 0])
    )
    """Switch for location of PF coil group i:
    - =1 PF coil on top of central solenoid (flux ramp only)
    - =2 PF coil on top of TF coil (flux ramp only)
    - =3 PF coil outside of TF coil (equilibrium coil)
    - =4 PF coil, general location (equilibrium coil)
    """

    i_pf_conductor: int = 0
    """switch for PF & CS coil conductor type:
    - =0 superconducting PF coils
    - =1 resistive PF coils
    """

    itr_sum: float = 0.0
    """total sum of I x turns x radius for all PF coils and CS (Am)"""

    i_cs_superconductor: int = 1
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

    i_pf_superconductor: int = 1
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

    j_crit_str_cs: float = 0.0
    """superconductor strand critical current density under operating
    conditions in central solenoid (A/m2). Necessary for the cost calculation in $/kA m
    """

    j_crit_str_pf: float = 0.0
    """superconductor strand critical current density under operating
    conditions in PF coils (A/m2). Necessary for the cost calculation in $/kA m
    """

    i_pf_current: int = 1
    """Switch for controlling the current of the PF coils:
    - =0 Input via the variables c_pf_cs_coil_pulse_start_ma, c_pf_cs_coil_flat_top_ma, c_pf_cs_coil_pulse_end_ma
    - =1 SVD targets zero field across midplane (flux swing
    coils) and the correct vertical field at the plasma
    center (equilibrium coils)
    """

    i_r_pf_outside_tf_placement: int = 0
    """Switch for the placement of Location 3 (outboard) PF coils
    - =0 (Default) Outboard PF coils follow TF shape
    in an ellipsoidal winding surface
    - =1 Outboard PF coils all have same radius, cylindrical
    winding surface
    """

    j_cs_conductor_critical_pulse_start: float = 0.0
    """central solenoid superconductor critical current density (A/m2) at beginning-of-pulse"""

    j_cs_conductor_critical_flat_top_end: float = 0.0
    """central solenoid superconductor critical current density (A/m2) at end-of-flattop"""

    jcableoh_bop: float = 0.0
    """central solenoid cable critical current density (A/m2) at beginning-of-pulse"""

    jcableoh_eof: float = 0.0
    """central solenoid cable critical current density (A/m2) at end-of-flattop"""

    n_pf_cs_plasma_circuits: int = 0
    """number of PF circuits (including central solenoid and plasma)"""

    n_pf_coils_in_group: list[int] = field(
        default_factory=lambda: np.array([1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    )
    """number of PF coils in group j"""

    n_cs_current_filaments: int = 7
    """number of filaments the top and bottom of the central solenoid should be broken
    into during scaling (5 - 10 is good)
    """

    n_pf_coil_groups: int = 3
    """number of groups of PF coils. Symmetric coil pairs should all be in the same group"""

    n_cs_pf_coils: int = 0
    """number of PF coils (excluding the central solenoid) + 1"""

    f_z_cs_tf_internal: float = 0.71
    """Central solenoid height / TF coil internal height"""

    f_a_cs_turn_steel: float = 0.5
    """Fraction of CS turn poloidal area that is steel (`iteration variable 122`)"""

    pf_current_safety_factor: float = 1.0
    """Ratio of permissible PF coil conductor current density to critical conductor
    current density based on short-sample DC measurements
    """

    pfcaseth: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """steel case thickness for PF coil i (m)"""

    rho_pf_coil: float = 2.5e-8
    """PF coil resistivity (if i_pf_conductor=1) (Ohm-m)"""

    rhopfbus: float = 3.93e-8
    """Resistivity of CS and PF coil bus bars (irrespective of
    whether the coils themselves are superconducting or resistive) (Ohm-m)
    """

    m_pf_coil_max: float = 0.0
    """mass of heaviest PF coil (tonnes)"""

    r_pf_coil_outer_max: float = 0.0
    """radius of largest PF coil (m)"""

    p_pf_electric_supplies_mw: float = 0.0
    """Total mean wall plug power dissipated in PFC and CS power supplies (MW) (issue #713)"""

    p_cs_resistive_flat_top: float = 0.0
    """central solenoid resistive power during flattop (W)"""

    p_pf_coil_resistive_total_flat_top: float = 0.0
    """total PF coil resistive losses during flattop (W)"""

    r_pf_coil_inner: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """inner radius of coil i (m)"""

    r_pf_coil_outer: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """outer radius of coil i (m)"""

    c_pf_cs_coils_peak_ma: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """peak current in coil i (MA-turns)"""

    j_pf_coil_wp_peak: list[float] = field(default_factory=lambda: np.full(NGC2, 3.0e7))
    """average winding pack current density of PF coil i (A/m2) at time of peak
    current in that coil (calculated for `i_pf_location=1` coils)
    """

    j_cs_critical_flat_top_end: float = 0.0
    """allowable central solenoid current density at end of flat-top (A/m2)"""

    j_cs_critical_pulse_start: float = 0.0
    """allowable central solenoid current density at beginning of pulse (A/m2)"""

    j_pf_wp_critical: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """allowable winding pack current density of PF coil i (A/m2)"""

    r_cs_middle: float = 0.0
    """radius to the centre of the central solenoid (m)"""

    dz_cs_full: float = 0.0
    """Full height of the central solenoid (m)"""

    dr_pf_tf_outboard_out_offset: float = 1.5
    """radial distance (m) from outboard TF coil leg to centre of `i_pf_location=3` PF coils"""

    r_pf_coil_middle: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """radius of PF coil i (m)"""

    dr_pf_cs_middle_offset: float = 0.0
    """offset (m) of radial position of `i_pf_location=1` PF coils from being directly above
    the central solenoid
    """

    rpf2: float = -1.63
    """offset (m) of radial position of `i_pf_location=2` PF coils from being at
    rmajor (offset = rpf2*triang*rminor)
    """

    rref: list[float] = field(default_factory=lambda: np.full(N_PF_GROUPS_MAX, 7.0))
    """PF coil radial positioning adjuster:
    - for groups j with i_pf_location(j) = 1; rref(j) is ignored
    - for groups j with i_pf_location(j) = 2; rref(j) is ignored
    - for groups j with i_pf_location(j) = 3; rref(j) is ignored
    - for groups j with i_pf_location(j) = 4; rref(j) is radius of
    the coil in units of minor radii from the major radius
    (r = rmajor + rref*rminor)
    """

    s_shear_cs_peak: float = 0.0
    """Maximum shear stress (Tresca criterion) coils/central solenoid [MPa]"""

    sigpfcalw: float = 500.0
    """maximum permissible tensile stress (MPa) in steel coil cases for superconducting
    PF coils (`i_pf_conductor=0`)
    """

    sigpfcf: float = 1.0
    """fraction of JxB hoop force supported by steel case for superconducting PF coils (`i_pf_conductor=0`)"""

    ind_pf_cs_plasma_mutual: list[float] = field(
        default_factory=lambda: np.zeros((NGC2, NGC2))
    )
    """mutual inductance matrix (H)"""

    temp_cs_superconductor_margin: float = 0.0
    """Central solenoid temperature margin (K)"""

    n_pf_coil_turns: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """number of turns in PF coil i"""

    f_a_pf_coil_void: list[float] = field(default_factory=lambda: np.full(NGC2, 0.3))
    """winding pack void fraction of PF coil i for coolant"""

    f_a_cs_void: float = 0.3
    """void fraction of central solenoid conductor for coolant"""

    vs_cs_pf_total_burn: float = 0.0
    """total flux swing available for burn (Wb)"""

    vs_pf_coils_total_burn: float = 0.0
    """flux swing from PF coils for burn (Wb)"""

    vs_pf_coils_total_ramp: float = 0.0
    """flux swing from PF coils for startup (Wb)"""

    vs_pf_coils_total_pulse: float = 0.0
    """total flux swing from PF coils (Wb)"""

    vs_cs_total_pulse: float = 0.0
    """total flux swing from the central solenoid (Wb)"""

    vs_cs_burn: float = 0.0
    """central solenoid flux swing for burn (Wb)"""

    vs_cs_ramp: float = 0.0
    """central solenoid flux swing for startup (Wb)"""

    vs_cs_pf_total_ramp: float = 0.0
    """total flux swing for startup (`constraint eqn 51` to enforce vs_cs_pf_total_ramp=vs_plasma_res_ramp+vs_plasma_ind_ramp) (Wb)"""

    vs_cs_pf_total_pulse: float = 0.0
    """total flux swing for pulse (Wb)"""

    f_c_pf_cs_peak_time_array: list[float] = field(
        default_factory=lambda: np.zeros((NGC2, 6))
    )
    """PF, CS coil current relative to peak current at time points 1 to 6"""

    m_pf_coil_conductor_total: float = 0.0
    """total mass of the PF coil conductor (kg)"""

    m_pf_coil_structure_total: float = 0.0
    """total mass of the PF coil structure (kg)"""

    m_pf_coil_conductor: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """conductor mass for PF coil i (kg)"""

    m_pf_coil_structure: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """structure mass for PF coil i (kg)"""

    z_pf_coil_upper: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """upper point of PF coil i (m)"""

    z_pf_coil_lower: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """lower point of PF coil i (m)"""

    z_pf_coil_middle: list[float] = field(default_factory=lambda: np.zeros(NGC2))
    """z (height) location of PF coil i (m)"""

    zref: list[float] = field(
        default_factory=lambda: np.array([
            3.6,
            1.2,
            2.5,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ])
    )
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

    b_cs_limit_max: float = 13.0
    """Central solenoid max field limit [T]"""

    f_dr_dz_cs_turn: float = 70.0 / 22.0
    """Ratio of CS coil turn conduit length to depth"""

    dr_cs_turn: float = 0.0
    """Length of CS of CS coil turn conduit"""

    dr_cs_full: float = 0.0
    """Full radial thickness of the central solenoid (m)"""

    dz_cs_turn: float = 0.0
    """Depth/width of CS of CS coil turn conduit"""

    radius_cs_turn_corners: float = 3.0e-3
    """Radius of curvature of CS coil turn corners (m)"""

    radius_cs_turn_cable_space: float = 0.0
    """Length of CS of CS coil turn conduit length"""


CREATE_DICTS_FROM_DATACLASS = PFCoilData
