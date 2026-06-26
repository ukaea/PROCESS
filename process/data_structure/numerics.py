from dataclasses import dataclass, field
from enum import IntEnum
from types import DynamicClassAttribute

import numpy as np


class PROCESSRunMode(IntEnum):
    """Enumeration of the available PROCESS run modes, which determine the behaviour
    of the code in various places. This is controlled by the `ioptimz` variable
    """

    EVALUATION = (-2, "Evaluation mode (no optimisation)")
    """In this mode, the code will not perform any optimisation, and will instead
    simply evaluate the constraints for the given input parameters, which is useful
    for testing and for evaluating the performance of a given design point without
    trying to optimise it. Internally, PROCESS uses `fsolve` (a Newton-Krylov/hybrd
    root-finding method from `scipy.optimize`) to seek a *consistent* solution by
    varying a subset of the iteration variables until the consistency constraints
    (equality constraints whose residuals must be driven to zero) are simultaneously
    satisfied; no figure-of-merit is optimised, and the solver simply tries to find
    a root of the constraint-residual vector.
    """
    OPTIMISATION = (1, "Optimisation mode (e.g. via VMCON)")
    """In this mode, the code will perform optimisation using the VMCON solver
    (or a custom solver if specified) to try to find a design point that optimises
    the figure of merit while satisfying the constraints.  This is the default mode
    of operation for PROCESS.
    """

    def __new__(cls, value: int, description: str):
        """Create a new PROCESSRunMode enum member with description.

        Args:
            value: The integer value of the enum member.
            description: The description for this run mode.

        Returns
        -------
            The new enum member with attached description.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._description_ = description
        return obj

    @DynamicClassAttribute
    def description(self):
        """The description for this run mode."""
        return self._description_


class FiguresOfMerit(IntEnum):
    """Enumeration of the available figures of merit (FoM) that can be used as
    objective functions for optimisation in PROCESS.
    """

    MAJOR_RADIUS = (1, "Plasma major radius (R₀)")
    NEUTRON_WALL_LOAD = (3, "Neutron wall load")
    P_TF_PLUS_P_PF = (4, "TF & PF coil power")
    FUSION_GAIN_Q = (5, "Fusion gain (Qₚₗₐₛₘₐ)")
    COST_OF_ELECTRICITY = (6, "Cost of electricity")
    CAPITAL_COST = (7, "Plant capital cost")
    ASPECT_RATIO = (8, "Plasma aspect ratio")
    DIVERTOR_HEAT_LOAD = (9, "Divertor heat load")
    TOROIDAL_FIELD = (10, "Plasma toroidal field on axis (B₀)")
    TOTAL_INJECTED_POWER = (11, "Plasma total injected power (Pᵢₙⱼ)")
    PULSE_LENGTH = (14, "Pulse length")
    PLANT_AVAILABILITY_FACTOR = (15, "Plant availability factor")
    MIN_R0_MAX_TAU_BURN = (
        16,
        "Linear combination of major radius (minimised) and pulse length (maximised)",
    )
    NET_ELECTRICAL_OUTPUT = (17, "Plant net electrical output")
    NULL_FIGURE_OF_MERIT = (18, "Null Figure of Merit")
    MAX_Q_MAX_T_PLANT_PULSE_BURN = (
        19,
        "Linear combination of big Q and pulse length (maximised)",
    )

    def __new__(cls, value: int, description: str):
        """Create a new FiguresOfMerit enum member with description.

        Args:
            value: The integer value of the enum member.
            description: The description for this figure of merit.

        Returns
        -------
            The new enum member with attached description.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._description_ = description
        return obj

    @DynamicClassAttribute
    def description(self):
        """The description for this figure of merit."""
        return self._description_


IPNVARS = 177
"""total number of variables available for iteration"""

IPEQNS = 92
"""number of constraint equations available"""

IPNFOMS = 19
"""number of available figures of merit"""


@dataclass(slots=True)
class NumericsData:
    ioptimz: int = 1
    """Code operation switch:
    * -2 for evaluation mode (i.e. no optimisation)
    * 1 for optimisation mode (e.g. via VMCON)
    """

    minmax: int = 7
    """
    Switch for figure-of-merit (see `FiguresOfMerit` for descriptions)
    negative => maximise, positive => minimise
    """

    n_constraints: int = 0
    """Total number of constraints (neqns + nineqns)"""

    ncalls: int = 0
    """number of function calls during solution"""

    neqns: int = -1
    """number of equality constraints to be satisfied"""

    nfev1: int = 0
    """number of calls to FCNHYB (HYBRD function caller) made"""

    nfev2: int = 0
    """number of calls to FCNVMC1 (VMCON function caller) made"""

    nineqns: int = 0
    """number of inequality constraints VMCON must satisfy
    (leave at zero for now)
    """

    nvar: int = 0
    """number of iteration variables to use"""

    nviter: int = 0
    """number of optimisation iterations performed"""

    icc: list[int] = field(default_factory=lambda: np.array([0] * IPEQNS))

    active_constraints: list[bool] = field(default_factory=lambda: [False] * IPEQNS)
    """Logical array showing which constraints are active"""

    lablcc: list[str] = field(
        default_factory=lambda: [
            "⟨β⟩ consistency                   ",
            "Global power balance consistency ",
            "Ion power balance                ",
            "Electron power balance           ",
            "Electron density upper limit (nₑ<) ",
            "(βₚε) upper limit ",
            "Beam ion density consistency     ",
            "Neutron wall load upper limit    ",
            "Fusion power upper limit         ",
            "NOT USED",
            "Radial build consistency         ",
            "CS & PF system whole pulse Vs lower limit          ",
            "Burn time lower limit            ",
            "NBI decay lengths consistency    ",
            "Pₛₑₚ > Pₗₕ consistency    ",
            "Net electric power lower limit   ",
            "Radiation fraction (fᵧ) upper limit   ",
            "Divertor heat load upper limit   ",
            "Resistive TF MVA upper limit                  ",
            "Beam tangency radius upper limit ",
            "Plasma minor radius (a) lower limit  ",
            "Divertor collisionality upper limit",
            "Conducting shell radius upper limit",
            "⟨β⟩ upper limit                 ",
            "TF peak toroidal field upper limit  ",
            "CS coil EOF current density limit",
            "CS coil BOP current density limit",
            "Fusion gain (Qₚₗₐₛₘₐ) lower limit        ",
            "Inboard radial build consistency ",
            "Plasma injected power (Pₐᵤₓ) upper limit      ",
            "TF coil case stress upper limit  ",
            "TF coil conduit stress upper limit ",
            "TF coil superconductor critical current density upper limit  ",
            "TF quench dump voltage upper limit         ",
            "TF quench hotspot current density upper limit         ",
            "TF coil superconductor temperature margin lower limit ",
            "HCD normalised current drive efficiency (γ) upper limit        ",  # noqa: RUF001
            "FW coolant temperature rise upper limit ",
            "FW peak temperature limit",
            "Plasma injected power lower limit  ",
            "Plasma current (Iₚ) ramp time lower limit ",
            "Pulse cycle time lower limit           ",
            "Average CP temperature consistency  ",
            "Peak CP temperature upper limit",
            "Edge safety factor (q₉₅) lower limit   ",
            "Iₚ/I_rod upper limit              ",
            "NOT USED",
            "⟨βₚ⟩ upper limit        ",
            "NOT USED",
            "IFE repetition rate upper limit  ",
            "CS & PF system ramp-up Vs consistency ",
            "Tritium breeding ratio lower limit ",
            "TF fast neutron fluence upper limit ",
            "TF peak nuclear heating upper limit ",
            "NOT USED",
            "Pₛₑₚ / R₀ upper limit             ",
            "NOT USED",
            "NOT USED",
            "NB shine-through fraction upper limit",
            "CS superconductor temperature margin lower limit",
            "Plant availability lower limit",
            "Alpha to energy confinement ratio (τ_α/τₑ) lower limit       ",  # noqa: RUF001
            "ITER-like vacuum pump number upper limit",
            "Plasma volume averaged effective charge (⟨Zₑ⟩) upper limit                       ",
            "VV stress during TF quench upper limit     ",
            "Rate of change of energy in PF system upper limit",
            "FW radiation wall load upper limit",
            "(PₛₑₚBₜ / q₉₅AR₀) upper limit      ",
            "NOT USED",
            "NOT USED",
            "NOT USED",
            "CS Tresca yield criterion upper limit     ",
            "Pₛₑₚ > Pₗₕ + Pₐᵤₓ consistency           ",
            "TF quench temperature < temp_croco_quench_max",
            "TF current/copper area < Max    ",
            "Eich critical separatrix density upper limit ",
            "TF current per turn upper limit ",
            "Reinke criterion divertor impurity fraction lower limit",
            "Peak CS field upper limit        ",
            "Pₛₑₚ lower limit                ",
            "nₑ₀ > nₑ_pedestal constraint     ",
            "toroidalgap > dx_tf_inboard_out_t",  # Stellarator constraint
            "available_space > required_space ",  # Stellarator constraint
            "⟨β⟩ lower limit                   ",
            "CP lifetime consistency                  ",
            "TF turn dimension upper limit              ",
            "Cryogenic plant power upper limit           ",
            "TF WP vertical strain upper limit    ",
            "CS current to copper area upper limit   ",
            "CS achievable stress load cycles lower limit           ",
            "ECRH ignitability                ",  # Stellarator constraint
            "Fuel composition consistency     ",
        ]
    )
    """Labels describing constraint equations (corresponding itvs)<UL>
    * ( 1) Beta (consistency equation) (itv 5)
    * ( 2) Global power balance (consistency equation) (itv 10,1,2,3,4,6,11)
    * ( 3) Ion power balance DEPRECATED (itv 10,1,2,3,4,6,11)
    * ( 4) Electron power balance DEPRECATED (itv 10,1,2,3,4,6,11)
    * ( 5) Density upper limit (itv 9,1,2,3,4,5,6)
    * ( 6) (Epsilon x beta poloidal) upper limit (itv 8,1,2,3,4,6)
    * ( 7) Beam ion density (NBI) (consistency equation) (itv 7)
    * ( 8) Neutron wall load upper limit (itv 14,1,2,3,4,6)
    * ( 9) Fusion power upper limit (itv 26,1,2,3,4,6)
    * (10) Toroidal field 1/R (consistency equation) (itv 12,1,2,3,13 )
    * (11) Radial build (consistency equation) (itv 3,1,13,16,29,42,61)
    * (12) Volt second lower limit (STEADY STATE) (itv 15,1,2,3)
    * (13) Burn time lower limit (PULSE) (itv 21,1,16,17,29,42,44,61)
    (itv 19,1,2,3,6)
    * (14) Neutral beam decay lengths to plasma centre (NBI) (consistency equation)
    * (15) LH power threshold limit (itv 103)
    * (16) Net electric power lower limit (itv 25,1,2,3)
    * (17) Radiation fraction upper limit (itv 28)
    * (18) Divertor heat load upper limit (itv 27)
    * (19) MVA upper limit (itv 30)
    * (20) Neutral beam tangency radius upper limit (NBI) (itv 33,31,3,13)
    * (21) Plasma minor radius lower limit (itv 32)
    * (22) Divertor collisionality upper limit (itv 34,43)
    * (23) Conducting shell to plasma minor radius ratio upper limit
    (itv 104,1,74)
    * (24) Beta upper limit (itv 36,1,2,3,4,6,18)
    * (25) Peak toroidal field upper limit (itv 35,3,13,29)
    * (26) Central solenoid EOF current density upper limit (i_pf_conductor=0)
    (itv 38,37,41,12)
    * (27) Central solenoid BOP current density upper limit (i_pf_conductor=0)
    (itv 39,37,41,12)
    * (28) Fusion gain Q lower limit (itv 45,47,40)
    * (29) Inboard radial build consistency (itv 3,1,13,16,29,42,61)
    * (30) Injection power upper limit (itv 46,47,11)
    * (31) TF coil case stress upper limit (SCTF) (itv 48,56,57,58,59,60,24)
    * (32) TF coil conduit stress upper limit (SCTF) (itv 49,56,57,58,59,60,24)
    * (33) I_op / I_critical (TF coil) (SCTF) (itv 50,56,57,58,59,60,24)
    * (34) Dump voltage upper limit (SCTF) (itv 51,52,56,57,58,59,60,24)
    * (35) TF Quench Hotspot J limit (SCTF) (itv 53,56,57,58,59,60,24)
    * (36) TF coil temperature margin lower limit (SCTF) (itv 54,55,56,57,58,59,60,24)
    * (37) Current drive gamma upper limit (itv 40,47)
    * (38) First wall coolant temperature rise upper limit (itv 62)
    * (39) First wall peak temperature upper limit (itv 63)
    * (40) Start-up injection power lower limit (PULSE) (itv 64)
    * (41) Plasma current ramp-up time lower limit (PULSE) (itv  66,65)
    * (42) Cycle time lower limit (PULSE) (itv 17,67,65)
    * (43) Average centrepost temperature
    (TART) (consistency equation) (itv 13,20,69,70)
    * (44) Peak centrepost temperature upper limit (TART) (itv 68,69,70)
    * (45) Edge safety factor lower limit (TART) (itv 71,1,2,3)
    * (46) Equation for Ip/Irod upper limit (TART) (itv 72,2,60)
    * (47) NOT USED
    * (48) Poloidal beta upper limit (itv 79,2,3,18)
    * (49) NOT USED
    * (50) IFE repetition rate upper limit (IFE)
    * (51) Startup volt-seconds consistency (PULSE) (itv 16,29,3,1)
    * (52) Tritium breeding ratio lower limit (itv 89,90,91)
    * (53) Neutron fluence on TF coil upper limit (itv 92,93,94)
    * (54) Peak TF coil nuclear heating upper limit (itv 95,93,94)
    * (55) NOT USED
    * (56) Pseparatrix/Rmajor upper limit (itv 97,1,3)
    * (57) NOT USED
    * (58) NOT USED
    * (59) Neutral beam shine-through fraction upper limit (NBI) (itv 105,6,19,4 )
    * (60) Central solenoid temperature margin lower limit (SCTF) (itv 106)
    * (61) Minimum availability value (itv 107)
    * (62) f_t_alpha_energy_confinement the ratio of particle to energy confinement times (itv 110)
    * (63) The number of ITER-like vacuum pumps n_iter_vacuum_pumps < tfno (itv 111)
    * (64) Zeff less than or equal to n_charge_plasma_effective_vol_avg_max (itv 112)
    * (65) Dump time set by VV loads (itv 56, 113)
    * (66) Limit on rate of change of energy in poloidal field
    (Use iteration variable 65(t_plant_pulse_plasma_current_ramp_up), 115)
    * (67) Simple Radiation Wall load limit (itv 116, 4,6)
    * (68) Psep * Bt / qAR upper limit (itv 117)
    * (69) ensure separatrix power = the value from Kallenbach divertor (itv 118)
    * (70) ensure that teomp = separatrix temperature in the pedestal profile,
    (itv 119 (temp_plasma_separatrix_kev))
    * (71) ensure that neomp = separatrix density (nd_plasma_separatrix_electron) x neratio
    * (72) central solenoid shear stress limit (Tresca yield criterion)
    * (73) Psep >= Plh + Paux
    * (74) TFC quench < temp_croco_quench_max
    * (75) TFC current/copper area < Maximum
    * (76) Eich critical separatrix density
    * (77) TF coil current per turn upper limit
    * (78) Reinke criterion impurity fraction lower limit
    * (79) Peak CS field upper limit
    * (80) Divertor power lower limit p_plasma_separatrix_mw
    * (81) Ne(0) > ne(ped) constraint
    * (82) toroidalgap >  dx_tf_inboard_out_toroidal constraint
    * (83) Radial build consistency for stellarators
    * (84) Lower limit for beta (itv 173 fbeta_min)
    * (85) Constraint for CP lifetime
    * (86) Constraint for TF coil turn dimension
    * (87) Constraint for cryogenic power
    * (88) Constraint for TF coil strain absolute value
    * (89) Constraint for CS coil quench protection
    * (90) Lower Limit on number of stress load cycles for CS
    * (91) Checking if the design point is ECRH ignitable
    * (92) D/T/He3 ratio in fuel sums to 1
    """

    ixc: list[int] = field(default_factory=lambda: np.array([0] * IPNVARS))
    """Array defining which iteration variables to activate
    (see lablxc for descriptions)
    """

    lablxc: list[str] = field(default_factory=lambda: [""] * IPNVARS)
    """Labels describing iteration variables<UL>
    * ( 1) aspect
    * ( 2) b_plasma_toroidal_on_axis
    * ( 3) rmajor
    * ( 4) temp_plasma_electron_vol_avg_kev
    * ( 5) beta_total_vol_avg
    * ( 6) nd_plasma_electrons_vol_avg
    * ( 7) f_nd_beam_electron
    * ( 8) NOT USED
    * ( 9) NOT USED
    * (10) hfact
    * (11) p_hcd_primary_extra_heat_mw
    * (12) j_tf_coil_full_area
    * (13) dr_tf_inboard (NOT RECOMMENDED)
    * (14) NOT USED
    * (15) NOT USED
    * (16) dr_cs
    * (17) t_plant_pulse_dwell
    * (18) q
    * (19) e_beam_kev
    * (20) temp_cp_average
    * (21) NOT USED
    * (22) NOT USED
    * (23) fcoolcp
    * (24) NOT USED
    * (25) NOT USED
    * (26) NOT USED
    * (27) NOT USED
    * (28) NOT USED
    * (29) dr_bore
    * (30) NOT USED
    * (31) gapomin
    * (32) NOT USED
    * (33) NOT USED
    * (34) NOT USED
    * (35) NOT USED
    * (36) NOT USED
    * (37) j_cs_flat_top_end
    * (38) NOT USED
    * (39) NOT USED
    * (40) NOT USED
    * (41) f_j_cs_start_pulse_end_flat_top
    * (42) dr_cs_tf_gap
    * (43) NOT USED
    * (44) f_c_plasma_non_inductive
    * (45) NOT USED
    * (46) NOT USED
    * (47) feffcd
    * (48) NOT USED
    * (49) NOT USED
    * (50) NOT USED
    * (51) NOT USED
    * (52) NOT USED
    * (53) NOT USED
    * (54) NOT USED
    * (55) NOT USED
    * (56) t_tf_superconductor_quench
    * (57) dr_tf_nose_case
    * (58) dx_tf_turn_steel
    * (59) f_a_tf_turn_cable_copper
    * (60) c_tf_turn
    * (61) dr_shld_vv_gap_inboard
    * (62) NOT USED
    * (63) NOT USED
    * (64) NOT USED
    * (65) t_plant_pulse_plasma_current_ramp_up
    * (66) NOT USED
    * (67) NOT USED
    * (68) NOT USED
    * (69) radius_cp_coolant_channel
    * (70) vel_cp_coolant_midplane
    * (71) NOT USED
    * (72) NOT USED
    * (73) dr_fw_plasma_gap_inboard
    * (74) dr_fw_plasma_gap_outboard
    * (75) f_dr_tf_outboard_inboard
    * (76) NOT USED
    * (77) NOT USED
    * (78) NOT USED
    * (79) NOT USED
    * (80) NOT USED
    * (81) edrive
    * (82) drveff
    * (83) tgain
    * (84) chrad
    * (85) pdrive
    * (86) NOT USED
    * (87) NOT USED
    * (88) NOT USED
    * (89) NOT USED
    * (90) blbuith
    * (91) blbuoth
    * (92) NOT USED
    * (93) dr_shld_inboard
    * (94) dr_shld_outboard
    * (95) NOT USED
    * (96) NOT USED
    * (97) NOT USED
    * (98) f_blkt_li6_enrichment
    * (99) NOT USED
    * (100) NOT USED
    * (101) NOT USED
    * (102) f_nd_impurity_electronsvar # OBSOLETE
    * (103) NOT USED
    * (104) NOT USED
    * (105) NOT USED
    * (106) NOT USED
    * (107) NOT USED
    * (108) breeder_f: Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
    * (109) f_nd_alpha_electron: thermal alpha density / electron density
    * (110) NOT USED
    * (111) NOT USED
    * (112) NOT USED
    * (113) NOT USED
    * (114) len_fw_channel: Length of a single first wall channel
    * (115) NOT USED
    * (116) NOT USED
    * (117) NOT USED
    * (119) temp_plasma_separatrix_kev:  separatrix temperature calculated by the Kallenbach divertor model
    * (120) ttarget: Plasma temperature adjacent to divertor sheath [eV]
    * (121) neratio: ratio of mean SOL density at OMP to separatrix density at OMP
    * (122) f_a_cs_turn_steel : streel fraction of Central Solenoid
    * (123) NOT USED
    * (124) qtargettotal : Power density on target including surface recombination [W/m2]
    * (125) f_nd_impurity_electrons(3) :  Beryllium density fraction relative to electron density
    * (126) f_nd_impurity_electrons(4) :  Carbon density fraction relative to electron density
    * (127) f_nd_impurity_electrons(5) :  Nitrogen fraction relative to electron density
    * (128) f_nd_impurity_electrons(6) :  Oxygen density fraction relative to electron density
    * (129) f_nd_impurity_electrons(7) :  Neon density fraction relative to electron density
    * (130) f_nd_impurity_electrons(8) :  Silicon density fraction relative to electron density
    * (131) f_nd_impurity_electrons(9) :  Argon density fraction relative to electron density
    * (132) f_nd_impurity_electrons(10) :  Iron density fraction relative to electron density
    * (133) f_nd_impurity_electrons(11) :  Nickel density fraction relative to electron density
    * (134) f_nd_impurity_electrons(12) :  Krypton density fraction relative to electron density
    * (135) f_nd_impurity_electrons(13) :  Xenon density fraction relative to electron density
    * (136) f_nd_impurity_electrons(14) :  Tungsten density fraction relative to electron density
    * (137) NOT USED
    * (138) dx_tf_hts_tape_rebco : thickness of REBCO layer in tape (m)
    * (139) dx_tf_hts_tape_copper : thickness of copper layer in tape (m)
    * (140) dr_tf_wp_with_insulation : radial thickness of TFC winding pack (m)
    * (141) NOT USED
    * (142) nd_plasma_separatrix_electron : electron density at separatrix [m-3]
    * (143) f_copperA_m2 : TF coil current / copper area < Maximum value
    * (144) NOT USED
    * (145) f_nd_plasma_pedestal_greenwald :  fraction of Greenwald density to set as pedestal-top density
    * (146) NOT USED
    * (147) NOT USED
    * (148) fzactual : fraction of impurity at SOL with Reinke detachment criterion
    * (149) NOT USED
    * (150) NOT USED
    * (151) NOT USED
    * (152) f_nd_plasma_separatrix_greenwald : Ratio of separatrix density to Greenwald density
    * (153) NOT USED
    * (154) NOT USED
    * (155) pfusife : IFE input fusion power (MW) (ifedrv=3 only)
    * (156) rrin : Input IFE repetition rate (Hz) (ifedrv=3 only)
    * (157) NOT USED
    * (158) dx_tf_croco_strand_copper : Thickness of CroCo copper tube (m)
    * (159) NOT USED
    * (160) NOT USED
    * (161) NOT USED
    * (162) r_cp_top : Top outer radius of the centropost (ST only) (m)
    * (163) NOT USED
    * (164) NOT USED
    * (165) NOT USED
    * (166) NOT USED
    * (167) NOT USED
    * (168) NOT USED
    * (169) te0_ecrh_achievable: Max. achievable electron temperature at ignition point
    * (170) deg_div_field_plate : field line angle wrt divertor target plate (degrees)
    * (171) casths_fraction : TF side case thickness as fraction of toridal case thickness
    * (172) dx_tf_side_case_min : TF side case thickness [m]
    * (173) f_plasma_fuel_deuterium : Deuterium fraction in fuel
    * (174) NOT USED
    * (175) NOT USED
    """
    # WARNING These labels are used as variable names by new_indat(), and possibly
    # other python utilities, so they cannot easily be changed.

    name_xc: list[str] = field(default_factory=lambda: [""] * IPNVARS)

    sqsumsq: float = 0.0
    """sqrt of the sum of the square of the constraint residuals"""

    objf_name: str = ""
    """Description of the objective function"""

    norm_objf: float = 0.0
    """Normalised objective function (figure of merit)"""

    epsfcn: float = 1.0e-3
    """Finite difference step length for calculating derivatives"""

    epsvmc: float = 1.0e-6
    """Error tolerance for optimiser"""

    boundl: list[float] = field(default_factory=lambda: np.array([9.0e-99] * IPNVARS))
    """Lower bounds used on ixc variables during
    optimisation runs
    """

    boundu: list[float] = field(default_factory=lambda: np.array([9.0e99] * IPNVARS))
    """Upper bounds used on ixc variables"""

    itv_scaled_lower_bounds: list[float] = field(
        default_factory=lambda: np.array([0.0] * IPNVARS)
    )
    """Lower bound of the ixc variables scaled to (divided by)
    the initial value of the corresponding ixc
    """

    itv_scaled_upper_bounds: list[float] = field(
        default_factory=lambda: np.array([0.0] * IPNVARS)
    )
    """Upper bound of the ixc variables scaled to (divided by)
    the initial value of the corresponding ixc
    """

    rcm: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))

    resdl: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))

    scafc: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))
    """The initial value of each ixc variable"""

    scale: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))
    """The reciprocal of the initial value of each ixc variable"""

    xcm: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))

    xcs: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))

    vlam: list[float] = field(default_factory=lambda: np.array([0.0] * IPNVARS))

    force_vmcon_inequality_satisfication: int = 1
    """If 1, adds an additional convergence criteria to the VMCON solver
    that enforces a margin on the inequality constraints.
    I.e. VMCON cannot converge until all inequality constraints are satisfied
    to within a tolerance of `force_vmcon_inequality_tolerance`.

    Default is 1 (enabled).

    NOTE: this only affects the VMCON solver.
    """

    force_vmcon_inequality_tolerance: float = 1e-8
    """The relative tolerance for the additional VMCON convergence criteria
    that forces inequality constraints to be satisfied within a tolerance.

    Default is 1e-8.

    NOTE: has no effect if `force_vmcon_inequality_satisfication` is 0
    NOTE: this only affects the VMCON solver.
    """

    factor: float = 0.1e0
    ftol: float = 1.0e-4


CREATE_DICTS_FROM_DATACLASS = NumericsData
