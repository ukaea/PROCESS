import numpy as np

ipnvars: int = 177
"""total number of variables available for iteration"""

ipeqns: int = 92
"""number of constraint equations available"""

ipnfoms: int = 19
"""number of available figures of merit"""

ipvlam: int = ipeqns + 2 * ipnvars + 1
iptnt: int = (ipeqns * (3 * ipeqns + 13)) / 2
ipvp1: int = ipnvars + 1

ioptimz: int = None
"""Code operation switch:
* -2 for evaluation mode (i.e. no optimisation)
* 1 for optimisation mode (e.g. via VMCON)
"""

minmax: int = None
"""
Switch for figure-of-merit (see lablmm for descriptions)
negative => maximise, positive => minimise
"""

lablmm: list[str] = None
"""lablmm(ipnfoms) : labels describing figures of merit:<UL>
* ( 1) major radius
* ( 2) not used
* ( 3) neutron wall load
* ( 4) P_tf + P_pf
* ( 5) fusion gain Q
* ( 6) cost of electricity
* ( 7) capital cost (direct cost if ireactor=0,
constructed cost otherwise)
* ( 8) aspect ratio
* ( 9) divertor heat load
* (10) toroidal field
* (11) total injected power
* (12) hydrogen plant capital cost OBSOLETE
* (13) hydrogen production rate OBSOLETE
* (14) pulse length
* (15) plant availability factor (N.B. requires
iavail=1 to be set)
* (16) linear combination of major radius (minimised) and pulse length (maximised)
note: FoM should be minimised only!
* (17) net electrical output
* (18) Null Figure of Merit
* (19) linear combination of big Q and pulse length (maximised)
note: FoM should be minimised only!
"""

n_constraints: int = None
"""Total number of constraints (neqns + nineqns)"""


ncalls: int = None
"""number of function calls during solution"""

neqns: int = None
"""number of equality constraints to be satisfied"""

nfev1: int = None
"""number of calls to FCNHYB (HYBRD function caller) made"""

nfev2: int = None
"""number of calls to FCNVMC1 (VMCON function caller) made"""

nineqns: int = None
"""number of inequality constraints VMCON must satisfy
(leave at zero for now)
"""

nvar: int = None
"""number of iteration variables to use"""

nviter: int = None
"""number of optimisation iterations performed"""

icc: list[int] = None

active_constraints: list[bool] = True
"""Logical array showing which constraints are active"""

# TODO Do not change the comments for lablcc: they are used to create the
# Python-Fortran dictionaries. This must be improved on.

lablcc: list[str] = None
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
* (35) J_winding pack/J_protection upper limit (SCTF) (itv 53,56,57,58,59,60,24)
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
* (55) Vacuum vessel helium concentration upper limit i_blanket_type =2 (itv 96,93,94)
* (56) Pseparatrix/Rmajor upper limit (itv 97,1,3)
* (57) NOT USED
* (58) NOT USED
* (59) Neutral beam shine-through fraction upper limit (NBI) (itv 105,6,19,4 )
* (60) Central solenoid temperature margin lower limit (SCTF) (itv 106)
* (61) Minimum availability value (itv 107)
* (62) f_alpha_energy_confinement the ratio of particle to energy confinement times (itv 110)
* (63) The number of ITER-like vacuum pumps n_iter_vacuum_pumps < tfno (itv 111)
* (64) Zeff less than or equal to zeff_max (itv 112)
* (65) Dump time set by VV loads (itv 56, 113)
* (66) Limit on rate of change of energy in poloidal field
(Use iteration variable 65(t_plant_pulse_plasma_current_ramp_up), 115)
* (67) Simple Radiation Wall load limit (itv 116, 4,6)
* (68) Psep * Bt / qAR upper limit (itv 117)
* (69) ensure separatrix power = the value from Kallenbach divertor (itv 118)
* (70) ensure that teomp = separatrix temperature in the pedestal profile,
(itv 119 (temp_plasma_separatrix_kev))
* (71) ensure that neomp = separatrix density (nd_plasma_separatrix_electron) x neratio
* (72) central solenoid shear stress limit (Tresca yield criterion) (itv 123 foh_stress)
* (73) Psep >= Plh + Paux (itv 137 (fplhsep))
* (74) TFC quench < temp_croco_quench_max (itv 141 (ftemp_croco_quench_max))
* (75) TFC current/copper area < Maximum (itv 143 f_coppera_m2)
* (76) Eich critical separatrix density
* (77) TF coil current per turn upper limit
* (78) Reinke criterion impurity fraction lower limit (itv  147 freinke)
* (79) Peak CS field upper limit (itv  149 fb_cs_limit_max)
* (80) Divertor power lower limit p_plasma_separatrix_mw (itv  153 fp_plasma_separatrix_min_mw)
* (81) Ne(0) > ne(ped) constraint (itv  154 fne0)
* (82) toroidalgap >  dx_tf_inboard_out_toroidal constraint (itv  171 ftoroidalgap)
* (83) Radial build consistency for stellarators (itv 172 f_avspace)
* (84) Lower limit for beta (itv 173 fbeta_min)
* (85) Constraint for CP lifetime
* (86) Constraint for TF coil turn dimension
* (87) Constraint for cryogenic power
* (88) Constraint for TF coil strain absolute value
* (89) Constraint for CS coil quench protection
* (90) Lower Limit on number of stress load cycles for CS (itr 167 fncycle)
* (91) Checking if the design point is ECRH ignitable (itv 168 fecrh_ignition)
* (92) D/T/He3 ratio in fuel sums to 1
"""

ixc: list[int] = None
"""Array defining which iteration variables to activate
(see lablxc for descriptions)
"""

lablxc: list[str] = None
"""Labels describing iteration variables<UL>
* ( 1) aspect
* ( 2) b_plasma_toroidal_on_axis
* ( 3) rmajor
* ( 4) temp_plasma_electron_vol_avg_kev
* ( 5) beta_total_vol_avg
* ( 6) nd_plasma_electrons_vol_avg
* ( 7) f_nd_beam_electron
* ( 8) fbeta_poloidal_eps (f-value for equation 6)
* ( 9) fdene (f-value for equation 5)
* (10) hfact
* (11) p_hcd_primary_extra_heat_mw
* (12) oacdcp
* (13) dr_tf_inboard (NOT RECOMMENDED)
* (14) fpflux_fw_neutron_max_mw (f-value for equation 8)
* (15) fvs_plasma_total_required (f-value for equation 12)
* (16) dr_cs
* (17) t_plant_pulse_dwell
* (18) q
* (19) e_beam_kev
* (20) temp_cp_average
* (21) ft_burn_min (f-value for equation 13)
* (22) NOT USED
* (23) fcoolcp
* (24) NOT USED
* (25) fp_plant_electric_net_required_mw (f-value for equation 16)
* (26) fp_fusion_total_max_mw (f-value for equation 9)
* (27) fpflux_div_heat_load_mw (f-value for equation 18)
* (28) fradpwr (f-value for equation 17), total radiation fraction
* (29) dr_bore
* (30) fmva (f-value for equation 19)
* (31) gapomin
* (32) frminor (f-value for equation 21)
* (33) fradius_beam_tangency (f-value for equation 20)
* (34) NOT USED
* (35) fb_tf_inboard_max (f-value for equation 25)
* (36) fbeta_max (f-value for equation 24)
* (37) j_cs_flat_top_end
* (38) fjohc (f-value for equation 26)
* (39) fjohc0 (f-value for equation 27)
* (40) feta_cd_norm_hcd_primary_max (f-value for equation 37)
* (41) f_j_cs_start_pulse_end_flat_top
* (42) dr_cs_tf_gap
* (43) NOT USED
* (44) f_c_plasma_non_inductive
* (45) fbig_q_plasma_min (f-value for equation 28)
* (46) fp_hcd_injected_max (f-value for equation 30)
* (47) feffcd
* (48) fstrcase (f-value for equation 31)
* (49) fstrcond (f-value for equation 32)
* (50) fiooic (f-value for equation 33)
* (51) fvdump (f-value for equation 34)
* (52) NOT USED
* (53) fjprot (f-value for equation 35)
* (54) ftmargtf (f-value for equation 36)
* (55) NOT USED
* (56) t_tf_superconductor_quench
* (57) dr_tf_nose_case
* (58) dx_tf_turn_steel
* (59) f_a_tf_turn_cable_copper
* (60) c_tf_turn
* (61) dr_shld_vv_gap_inboard
* (62) fdtmp (f-value for equation 38)
* (63) ftemp_fw_max (f-value for equation 39)
* (64) fp_hcd_injected_min_mw (f-value for equation 40)
* (65) t_plant_pulse_plasma_current_ramp_up
* (66) ft_current_ramp_up (f-value for equation 41)
* (67) ft_cycle_min (f-value for equation 42)
* (68) fptemp (f-value for equation 44)
* (69) rcool
* (70) vcool
* (71) fq95_min (f-value for equation 45)
* (72) fipir (f-value for equation 46)
* (73) dr_fw_plasma_gap_inboard
* (74) dr_fw_plasma_gap_outboard
* (75) f_dr_tf_outboard_inboard
* (76) NOT USED
* (77) NOT USED
* (78) NOT USED
* (79) fbeta_poloidal (f-value for equation 48)
* (80) NOT USED
* (81) edrive
* (82) drveff
* (83) tgain
* (84) chrad
* (85) pdrive
* (86) frrmax (f-value for equation 50)
* (87) NOT USED
* (88) NOT USED
* (89) ftbr (f-value for equation 52)
* (90) blbuith
* (91) blbuoth
* (92) fflutf (f-value for equation 53)
* (93) dr_shld_inboard
* (94) dr_shld_outboard
* (95) fptfnuc (f-value for equation 54)
* (96) fvvhe (f-value for equation 55)
* (97) fpsepr (f-value for equation 56)
* (98) f_blkt_li6_enrichment
* (99) NOT USED
* (100) NOT USED
* (101) NOT USED
* (102) f_nd_impurity_electronsvar # OBSOLETE
* (103) fl_h_threshold (f-value for equation 15)
* (104)fr_conducting_wall (f-value for equation 23)
* (105) fnbshinef (f-value for equation 59)
* (106) ftmargoh (f-value for equation 60)
* (107) favail (f-value for equation 61)
* (108) breeder_f: Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
* (109) f_nd_alpha_electron: thermal alpha density / electron density
* (110) falpha_energy_confinement: Lower limit on f_alpha_energy_confinement the ratio of alpha
* (111) fniterpump: f-value for constraint that number
* (112) fzeff_max: f-value for max Zeff (f-value for equation 64)
* (113) ftaucq: f-value for minimum quench time (f-value for equation 65)
* (114) len_fw_channel: Length of a single first wall channel
* (115) fpoloidalpower: f-value for max rate of change of
* (116) fpflux_fw_rad_max: f-value for radiation wall load limit (eq. 67)
* (117) fpsepbqar: f-value for  Psep*Bt/qar upper limit (eq. 68)
* (119) temp_plasma_separatrix_kev:  separatrix temperature calculated by the Kallenbach divertor model
* (120) ttarget: Plasma temperature adjacent to divertor sheath [eV]
* (121) neratio: ratio of mean SOL density at OMP to separatrix density at OMP
* (122) f_a_cs_turn_steel : streel fraction of Central Solenoid
* (123) foh_stress : f-value for CS coil Tresca yield criterion (f-value for eq. 72)
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
* (137) fplhsep (f-value for equation 73)
* (138) dx_hts_tape_rebco : thickness of REBCO layer in tape (m)
* (139) dx_hts_tape_copper : thickness of copper layer in tape (m)
* (140) dr_tf_wp_with_insulation : radial thickness of TFC winding pack (m)
* (141) ftemp_croco_quench_max : TF coil quench temperature < temp_croco_quench_max (f-value for equation 74)
* (142) nd_plasma_separatrix_electron : electron density at separatrix [m-3]
* (143) f_copperA_m2 : TF coil current / copper area < Maximum value
* (144) fnesep : Eich critical electron density at separatrix
* (145) f_nd_plasma_pedestal_greenwald :  fraction of Greenwald density to set as pedestal-top density
* (146) fc_tf_turn_max : F-value for TF coil current per turn limit (constraint equation 77)
* (147) freinke : F-value for Reinke detachment criterion (constraint equation 78)
* (148) fzactual : fraction of impurity at SOL with Reinke detachment criterion
* (149) fb_cs_limit_max : F-value for max peak CS field (con. 79, itvar 149)
* (150) REMOVED
* (151) REMOVED
* (152) f_nd_plasma_separatrix_greenwald : Ratio of separatrix density to Greenwald density
* (153) fp_plasma_separatrix_min_mw : F-value for minimum p_plasma_separatrix_mw (con. 80)
* (154) fne0 : F-value for ne(0) > ne(ped) (con. 81)
* (155) pfusife : IFE input fusion power (MW) (ifedrv=3 only)
* (156) rrin : Input IFE repetition rate (Hz) (ifedrv=3 only)
* (157) fvs_cs_pf_total_ramp : F-value for available to required start up flux (con. 51)
* (158) dx_croco_strand_copper : Thickness of CroCo copper tube (m)
* (159) ftoroidalgap : F-value for toroidalgap >  dx_tf_inboard_out_toroidal constraint (con. 82)
* (160) f_avspace (f-value for equation 83)
* (161) fbeta_min (f-value for equation 84)
* (162) r_cp_top : Top outer radius of the centropost (ST only) (m)
* (163) f_t_turn_tf : f-value for TF coils WP trurn squared dimension constraint
* (164) f_crypmw : f-value for cryogenic plant power
* (165) fstr_wp : f-value for TF coil strain absolute value
* (166) f_copperaoh_m2 : CS coil current /copper area < Maximum value
* (167) fncycle : f-value for minimum CS coil stress load cycles
* (168) fecrh_ignition: f-value for equation 91
* (169) te0_ecrh_achievable: Max. achievable electron temperature at ignition point
* (170) deg_div_field_plate : field line angle wrt divertor target plate (degrees)
* (171) casths_fraction : TF side case thickness as fraction of toridal case thickness
* (172) dx_tf_side_case_min : TF side case thickness [m]
* (173) f_plasma_fuel_deuterium : Deuterium fraction in fuel
* (174) EMPTY : Description
* (175) EMPTY : Description
"""
# Issue 287 iteration variables are now defined in module define_iteration_variables in iteration variables.f90

name_xc: list[str] = None

sqsumsq: float = None
"""sqrt of the sum of the square of the constraint residuals"""

objf_name: str = None
"""Description of the objective function"""

norm_objf: float = None
"""Normalised objective function (figure of merit)"""

epsfcn: float = None
"""Finite difference step length for calculating derivatives"""

epsvmc: float = None
"""Error tolerance for optimiser"""

boundl: list[float] = None
"""Lower bounds used on ixc variables during
optimisation runs
"""

boundu: list[float] = None
"""Upper bounds used on ixc variables"""

itv_scaled_lower_bounds: list[float] = None
"""Lower bound of the ixc variables scaled to (divided by)
the initial value of the corresponding ixc
"""

itv_scaled_upper_bounds: list[float] = None
"""Upper bound of the ixc variables scaled to (divided by)
the initial value of the corresponding ixc
"""

rcm: list[float] = None

resdl: list[float] = None

scafc: list[float] = None
"""The initial value of each ixc variable"""

scale: list[float] = None
"""The reciprocal of the initial value of each ixc variable"""

xcm: list[float] = None

xcs: list[float] = None

vlam: list[float] = None

force_vmcon_inequality_satisfication: int = None
"""If 1, adds an additional convergence criteria to the VMCON solver
that enforces a margin on the inequality constraints.
I.e. VMCON cannot converge until all inequality constraints are satisfied
to within a tolerance of `force_vmcon_inequality_tolerance`.

Default is 1 (enabled).

NOTE: this only affects the VMCON solver.
"""

force_vmcon_inequality_tolerance: float = None
"""The relative tolerance for the additional VMCON convergence criteria
that forces inequality constraints to be satisfied within a tolerance.

Default is 1e-8.

NOTE: has no effect if `force_vmcon_inequality_satisfication` is 0
NOTE: this only affects the VMCON solver.
"""


def init_numerics():
    global ipnvars
    global ipeqns
    global ipnfoms
    global ipvlam
    global iptnt
    global ipvp1
    global ioptimz
    global minmax
    global lablmm
    global n_constraints
    global ncalls
    global neqns
    global nfev1
    global nfev2
    global nineqns
    global nvar
    global nviter
    global icc
    global active_constraints
    global lablcc
    global ixc
    global lablxc
    global name_xc
    global sqsumsq
    global objf_name
    global norm_objf
    global epsfcn
    global epsvmc
    global factor
    global ftol
    global boundl
    global boundu
    global itv_scaled_lower_bounds
    global itv_scaled_upper_bounds
    global rcm
    global resdl
    global scafc
    global scale
    global xcm
    global xcs
    global vlam
    global force_vmcon_inequality_satisfication
    global force_vmcon_inequality_tolerance

    """Initialise module variables"""
    ioptimz = 1
    minmax = 7
    lablmm = [
        "major radius          ",
        "not used              ",
        "neutron wall load     ",
        "P_tf + P_pf           ",
        "fusion gain           ",
        "cost of electricity   ",
        "capital cost          ",
        "aspect ratio          ",
        "divertor heat load    ",
        "toroidal field        ",
        "total injected power  ",
        "H plant capital cost  ",
        "H production rate     ",
        "pulse length          ",
        "plant availability    ",
        "min R0, max tau_burn  ",
        "net electrical output ",
        "Null figure of merit  ",
        "max Q, max t_plant_pulse_burn     ",
    ]

    ncalls = 0
    neqns = 0
    nfev1 = 0
    nfev2 = 0
    nineqns = 0
    nvar = 0
    n_constraints = 0
    nviter = 0
    icc = np.array([0] * ipeqns)
    active_constraints = [False] * ipeqns

    lablcc = [
        "Beta consistency                 ",
        "Global power balance consistency ",
        "Ion power balance                ",
        "Electron power balance           ",
        "Density upper limit              ",
        "(Epsilon x beta-pol) upper limit ",
        "Beam ion density consistency     ",
        "Neutron wall load upper limit    ",
        "Fusion power upper limit         ",
        "Toroidal field 1/R consistency   ",
        "Radial build consistency         ",
        "Volt second lower limit          ",
        "Burn time lower limit            ",
        "NBI decay lengths consistency    ",
        "L-H power threshold limit        ",
        "Net electric power lower limit   ",
        "Radiation fraction upper limit   ",
        "Divertor heat load upper limit   ",
        "MVA upper limit                  ",
        "Beam tangency radius upper limit ",
        "Plasma minor radius lower limit  ",
        "Divertor collisionality upper lim",
        "Conducting shell radius upper lim",
        "Beta upper limit                 ",
        "Peak toroidal field upper limit  ",
        "CS coil EOF current density limit",
        "CS coil BOP current density limit",
        "Fusion gain Q lower limit        ",
        "Inboard radial build consistency ",
        "Injection power upper limit      ",
        "TF coil case stress upper limit  ",
        "TF coil conduit stress upper lim ",
        "I_op / I_critical (TF coil)      ",
        "Dump voltage upper limit         ",
        "J_winding pack/J_protection limit",
        "TF coil temp. margin lower limit ",
        "Current drive gamma limit        ",
        "1st wall coolant temp rise limit ",
        "First wall peak temperature limit",
        "Start-up inj. power lower limit  ",
        "Plasma curr. ramp time lower lim ",
        "Cycle time lower limit           ",
        "Average centrepost temperature   ",
        "Peak centrepost temp. upper limit",
        "Edge safety factor lower limit   ",
        "Ip/Irod upper limit              ",
        "TF coil tor. thickness upper lim ",
        "Poloidal beta upper limit        ",
        "RFP reversal parameter < 0       ",
        "IFE repetition rate upper limit  ",
        "Startup volt-seconds consistency ",
        "Tritium breeding ratio lower lim ",
        "Neutron fluence on TF coil limit ",
        "Peak TF coil nucl. heating limit ",
        "Vessel helium concentration limit",
        "Psep / R upper limit             ",
        "TF coil leg rad width lower limit",
        "TF coil leg rad width lower limit",
        "NB shine-through frac upper limit",
        "CS temperature margin lower limit",
        "Minimum availability value       ",
        "f_alpha_energy_confinement       ",
        "number of ITER-like vacuum pumps ",
        "Zeff limit                       ",
        "Dump time set by VV stress       ",
        "Rate of change of energy in field",
        "Upper Lim. on Radiation Wall load",
        "Upper Lim. on Psep * Bt / q A R  ",
        "p_sep < psep_kallenbach divertor ",
        "Separatrix temp consistency      ",
        "Separatrix density consistency   ",
        "CS Tresca yield criterion        ",
        "Psep >= Plh + Paux               ",
        "TFC quench <temp_croco_quench_max",
        "TFC current/copper area < Max    ",
        "Eich critical separatrix density ",
        "TFC current per turn upper limit ",
        "Reinke criterion fZ lower limit  ",
        "Peak CS field upper limit        ",
        "p_sep lower limit                ",
        "nd_plasma_electron_on_axis > nd_plasma_pedestal_electron                      ",
        "toroidalgap > dx_tf_inboard_out_t",
        "available_space > required_space ",
        "beta > beta_vol_avg_min                  ",
        "CP lifetime                      ",
        "TFC turn dimension               ",
        "Cryogenic plant power            ",
        "TF coil strain absolute value    ",
        "CS current/copper area < Max     ",
        "CS stress load cycles            ",
        "ECRH ignitability                ",
        "Fuel composition consistency     ",
    ]

    ixc = np.array([0] * ipnvars)

    # WARNING These labels are used as variable names by write_new_in_dat.py, and possibly
    # other python utilities, so they cannot easily be changed.
    lablxc = [""] * ipnvars

    sqsumsq = 0.0
    objf_name = ""
    norm_objf = 0.0
    epsfcn = 1.0e-3
    epsvmc = 1.0e-6
    factor = 0.1e0
    ftol = 1.0e-4

    boundl = np.array([9.0e-99] * ipnvars)
    boundu = np.array([9.0e99] * ipnvars)

    itv_scaled_lower_bounds = np.array([0.0] * ipnvars)
    itv_scaled_upper_bounds = np.array([0.0] * ipnvars)
    rcm = np.array([0.0] * ipnvars)
    resdl = np.array([0.0] * ipnvars)
    scafc = np.array([0.0] * ipnvars)
    scale = np.array([0.0] * ipnvars)
    xcm = np.array([0.0] * ipnvars)
    xcs = np.array([0.0] * ipnvars)
    vlam = np.array([0.0] * ipnvars)
    name_xc = [""] * ipnvars
    force_vmcon_inequality_satisfication = 1
    force_vmcon_inequality_tolerance = 1e-8
