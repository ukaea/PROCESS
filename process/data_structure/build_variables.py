aplasmin: float = None
"""minimum minor radius (m)"""


available_radial_space: float = None
"""Minimal radial space between plasma and coils (m)"""


a_blkt_total_surface: float = None
"""blanket total surface area (m2)"""


a_blkt_inboard_surface: float = None
"""inboard blanket surface area (m2)"""


a_blkt_outboard_surface: float = None
"""outboard blanket surface area (m2)"""


blbmith: float = None
"""inboard blanket box manifold thickness (m) (`blktmodel>0`)"""


blbmoth: float = None
"""outboard blanket box manifold thickness (m) (`blktmodel>0`)"""


blbpith: float = None
"""inboard blanket base plate thickness (m) (`blktmodel>0`)"""


blbpoth: float = None
"""outboard blanket base plate thickness (m) (`blktmodel>0`)"""


blbuith: float = None
"""inboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 90`)"""


blbuoth: float = None
"""outboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 91`)"""


dr_blkt_inboard: float = None
"""inboard blanket thickness (m); (calculated if `blktmodel>0`) (=0.0 if `i_blkt_inboard=0`)"""


dr_blkt_outboard: float = None
"""outboard blanket thickness (m); calculated if `blktmodel>0`"""


dz_blkt_upper: float = None
"""top blanket thickness (m), = mean of inboard and outboard blanket thicknesses"""


dz_fw_upper: float = None
"""upper first wall thickness (m)"""


dr_bore: float = None
"""central solenoid inboard radius (m) (`iteration variable 29`)"""


f_z_cryostat: float = None
"""cryostat lid height scaling factor (tokamaks)"""


dr_cryostat: float = None
"""cryostat thickness (m)"""


dr_vv_inboard: float = None
"""vacuum vessel inboard thickness (TF coil / shield) (m)"""


dr_vv_outboard: float = None
"""vacuum vessel outboard thickness (TF coil / shield) (m)"""


dz_vv_upper: float = None
"""vacuum vessel topside thickness (TF coil / shield) (m) (= dz_vv_lower if double-null)"""


dz_vv_lower: float = None
"""vacuum vessel underside thickness (TF coil / shield) (m)"""


dr_vv_shells: float = None
"""vacuum vessel double walled shell thicknesses (m)"""


fcspc: float = None
"""Fraction of space occupied by CS pre-compression structure"""


fseppc: float = None
"""Separation force in CS coil pre-compression structure"""


a_fw_total_full_coverage: float = None
"""First wall total surface area with no holes or ports [m^2]"""


a_fw_inboard_full_coverage: float = None
"""Inboard first wall surface area with no holes or ports [m^2]"""


a_fw_outboard_full_coverage: float = None
"""Outboard first wall surface area with no holes or ports [m^2]"""

a_fw_total: float = None
"""First wall total surface area [m^2]"""


a_fw_inboard: float = None
"""Inboard first wall surface area [m^2]"""


a_fw_outboard: float = None
"""Outboard first wall surface area [m^2]"""


dr_fw_inboard: float = None
"""inboard first wall thickness, initial estimate as calculated (m)"""


dr_fw_outboard: float = None
"""outboard first wall thickness, initial estimate as calculated (m)"""


dr_shld_vv_gap_inboard: float = None
"""gap between inboard vacuum vessel and thermal shield (m) (`iteration variable 61`)"""


dr_cs_tf_gap: float = None
"""gap between central solenoid and TF coil (m) (`iteration variable 42`)"""


gapomin: float = None
"""minimum gap between outboard vacuum vessel and TF coil (m) (`iteration variable 31`)"""


dr_shld_vv_gap_outboard: float = None
"""gap between outboard vacuum vessel and TF coil (m)"""


z_tf_inside_half: float = None
"""maximum (half-)height of TF coil (inside edge) (m)"""


dz_tf_upper_lower_midplane: float = None
"""difference in distance from midplane of upper and lower portions of TF
legs (non-zero for single-null devices) (m)
"""


z_tf_top: float = None
"""height to top of (upper) TF coil leg (m)"""


hr1: float = None
"""half-height of TF coil inboard leg straight section (m)"""


iohcl: int = None
"""Switch for existence of central solenoid:
- =0 central solenoid not present
- =1 central solenoid exists
"""


i_cs_precomp: int = None
"""Switch for existence of central solenoid pre-compression structure:
- =0 no pre-compression structure
- =1 calculated pre-compression structure
"""


i_tf_inside_cs: int = None
"""Switch for placing the TF coil inside the CS
- = 0 TF coil is outside the CS (default)
- = 1 TF coil is inside the CS
"""


dr_cs: float = None
"""Central solenoid thickness (m) (`iteration variable 16`)"""


dr_cs_precomp: float = None
"""CS coil precompression structure thickness (m)"""


rbld: float = None
"""sum of thicknesses to the major radius (m)"""


required_radial_space: float = None
"""Required space between coil and plasma for blanket shield wall etc (m)"""


rinboard: float = None
"""plasma inboard radius (m) (`consistency equation 29`)"""


r_shld_inboard_inner: float = None
"""radius to inboard shield (inside point) (m)"""


r_shld_outboard_outer: float = None
"""radius to outboard shield (outside point) (m)"""


r_vv_inboard_out: float = None
"""Radial plasma facing side position of inboard vacuum vessel [m]"""


r_sh_inboard_in: float = None
"""Radial inner side position of inboard neutronic shield [m]"""


r_sh_inboard_out: float = None
"""Radial plasma facing side position of inboard neutronic shield [m]"""


r_tf_inboard_in: float = None
"""Mid-plane inboard TF coil leg radius at the centre-machine side [m]"""


r_tf_inboard_mid: float = None
"""Mid-plane inboard TF coil leg radius at middle of the coil [m]"""


r_tf_inboard_out: float = None
"""Mid-plane inboard TF coil leg radius at the plasma side [m]"""


r_tf_outboard_mid: float = None
"""Mid-plane outboard TF coil leg radius at the middle of the coil [m]"""


i_r_cp_top: int = None
"""Switch selecting the he parametrization of the outer radius of the top of the CP part of the TF coil
0 : `r_cp_top` is set by the plasma shape
1 : `r_cp_top` is a user input
2 : `r_cp_top` is set using the CP top and midplane CP radius ratio
"""


r_cp_top: float = None
"""Top outer radius of the centropost (ST only) (m)"""


f_r_cp: float = None
"""Ratio between the top and the midplane TF CP outer radius [-]
Not used by default (-1) must be larger than 1 otherwise
"""


dr_tf_inner_bore: float = None
"""TF coil horizontal inner dr_bore (m)"""


dh_tf_inner_bore: float = None
"""TF coil vertical inner dr_bore (m)"""


dr_fw_plasma_gap_inboard: float = None
"""Gap between plasma and first wall, inboard side (m) (if `i_plasma_wall_gap=1`)
Iteration variable: ixc = 73
Scan variable: nsweep = 58
"""


dr_fw_plasma_gap_outboard: float = None
"""Gap between plasma and first wall, outboard side (m) (if `i_plasma_wall_gap=1`)
Iteration variable: ixc = 74
Scan variable: nsweep = 59
"""


a_shld_total_surface: float = None
"""shield total surface area (m2)"""


a_shld_inboard_surface: float = None
"""inboard shield surface area (m2)"""


a_shld_outboard_surface: float = None
"""outboard shield surface area (m2)"""


dr_shld_inboard: float = None
"""inboard shield thickness (m) (`iteration variable 93`)"""


dz_shld_lower: float = None
"""lower (under divertor) shield thickness (m)"""


dr_shld_outboard: float = None
"""outboard shield thickness (m) (`iteration variable 94`)"""


dz_shld_upper: float = None
"""upper/lower shield thickness (m); calculated if `blktmodel > 0` (= dz_shld_lower if double-null)"""


sigallpc: float = None
"""allowable stress in CSpre-compression structure (Pa)"""


# TODO: Issue #514 Make dr_tf_inboard an output not an iteration variable
dr_tf_inboard: float = None
"""inboard TF coil thickness, (centrepost for ST) (m)
(input, calculated or `iteration variable 13`)
"""


dz_tf_plasma_centre_offset: float = None
"""Vertical distance between centre of TF coils and centre of plasma (m)"""


f_dr_tf_outboard_inboard: float = None
"""TF coil outboard leg / inboard leg radial thickness
ratio (`i_tf_sup=0` only) (`iteration variable 75`)
"""


dr_tf_outboard: float = None
"""Outboard TF coil thickness (m)"""


dr_tf_shld_gap: float = None
"""Minimum metal-to-metal gap between TF coil and thermal shield (m)"""


dr_shld_thermal_inboard: float = None
"""TF-VV thermal shield thickness, inboard (m)"""


dr_shld_thermal_outboard: float = None
"""TF-VV thermal shield thickness, outboard (m)"""


dz_shld_thermal: float = None
"""TF-VV thermal shield thickness, vertical build (m)"""


dz_shld_vv_gap: float = None
"""vertical gap between vacuum vessel and thermal shields (m)"""


dz_xpoint_divertor: float = None
"""vertical gap between x-point and divertor (m) (if = 0, it is calculated)"""


dz_fw_plasma_gap: float = None
"""vertical gap between top of plasma and first wall (m) (= dz_xpoint_divertor if double-null)"""


dr_shld_blkt_gap: float = None
"""gap between vacuum vessel and blanket (m)"""


plleni: float = None
"""length of inboard divertor plate (m)"""


plleno: float = None
"""length of outboard divertor plate (m)"""


plsepi: float = None
"""poloidal length, x-point to inboard strike point (m)"""


plsepo: float = None
"""poloidal length, x-point to outboard strike point (m)"""


rspo: float = None
"""outboard strike point radius (m)"""


z_plasma_xpoint_upper: float = None
"""Vertical height of the upper plasma x-point (m)"""


z_plasma_xpoint_lower: float = None
"""Vertical height of the lower plasma x-point (m)"""

ripflag: int = None
"""1 if the fitted range of applicability is exceeded for
the ripple calculation, else 0
"""


def init_build_variables():
    global \
        ripflag, \
        aplasmin, \
        available_radial_space, \
        a_blkt_total_surface, \
        a_blkt_inboard_surface, \
        a_blkt_outboard_surface, \
        blbmith, \
        blbmoth, \
        blbpith, \
        blbpoth, \
        blbuith, \
        blbuoth, \
        dr_blkt_inboard, \
        dr_blkt_outboard, \
        dz_blkt_upper, \
        dz_fw_upper, \
        dr_bore, \
        f_z_cryostat, \
        dr_cryostat, \
        dr_vv_inboard, \
        dr_vv_outboard, \
        dz_vv_upper, \
        dz_vv_lower, \
        dr_vv_shells, \
        fcspc, \
        fseppc, \
        a_fw_total_full_coverage, \
        a_fw_inboard_full_coverage, \
        a_fw_outboard_full_coverage, \
        a_fw_total, \
        a_fw_inboard, \
        a_fw_outboard, \
        dr_fw_inboard, \
        dr_fw_outboard, \
        dr_shld_vv_gap_inboard, \
        dr_cs_tf_gap, \
        gapomin, \
        dr_shld_vv_gap_outboard, \
        z_tf_inside_half, \
        dz_tf_upper_lower_midplane, \
        z_tf_top, \
        hr1, \
        iohcl, \
        i_cs_precomp, \
        i_tf_inside_cs, \
        dr_cs, \
        dr_cs_precomp, \
        rbld, \
        required_radial_space, \
        rinboard, \
        r_shld_inboard_inner, \
        r_shld_outboard_outer, \
        r_vv_inboard_out, \
        r_sh_inboard_in, \
        r_sh_inboard_out, \
        r_tf_inboard_in, \
        r_tf_inboard_mid, \
        r_tf_inboard_out, \
        r_tf_outboard_mid, \
        i_r_cp_top, \
        r_cp_top, \
        f_r_cp, \
        dr_tf_inner_bore, \
        dh_tf_inner_bore, \
        dr_fw_plasma_gap_inboard, \
        dr_fw_plasma_gap_outboard, \
        a_shld_total_surface, \
        a_shld_inboard_surface, \
        a_shld_outboard_surface, \
        dr_shld_inboard, \
        dz_shld_lower, \
        dr_shld_outboard, \
        dz_shld_upper, \
        sigallpc, \
        dr_tf_inboard, \
        dz_tf_plasma_centre_offset, \
        f_dr_tf_outboard_inboard, \
        dr_tf_outboard, \
        dr_tf_shld_gap, \
        dr_shld_thermal_inboard, \
        dr_shld_thermal_outboard, \
        dz_shld_thermal, \
        dz_shld_vv_gap, \
        dz_xpoint_divertor, \
        dz_fw_plasma_gap, \
        dr_shld_blkt_gap, \
        plleni, \
        plleno, \
        plsepi, \
        plsepo, \
        rspo, \
        z_plasma_xpoint_upper, \
        z_plasma_xpoint_lower

    ripflag = 0
    aplasmin = 0.25
    available_radial_space = 0.0
    a_blkt_total_surface = 0.0
    a_blkt_inboard_surface = 0.0
    a_blkt_outboard_surface = 0.0
    blbmith = 0.17
    blbmoth = 0.27
    blbpith = 0.30
    blbpoth = 0.35
    blbuith = 0.365
    blbuoth = 0.465
    dr_blkt_inboard = 0.115
    dr_blkt_outboard = 0.235
    dz_blkt_upper = 0.0
    dz_fw_upper = 0.0
    dr_bore = 1.42
    f_z_cryostat = 4.268
    dr_cryostat = 0.07
    dr_vv_inboard = 0.07
    dr_vv_outboard = 0.07
    dz_vv_upper = 0.07
    dz_vv_lower = 0.07
    dr_vv_shells = 0.12
    fcspc = 0.6
    fseppc = 3.5e8
    a_fw_total_full_coverage = 0.0
    a_fw_inboard_full_coverage = 0.0
    a_fw_outboard_full_coverage = 0.0
    a_fw_total = 0.0
    a_fw_inboard = 0.0
    a_fw_outboard = 0.0
    dr_fw_inboard = 0.0
    dr_fw_outboard = 0.0
    dr_shld_vv_gap_inboard = 0.155
    dr_cs_tf_gap = 0.08
    gapomin = 0.234
    dr_shld_vv_gap_outboard = 0.0
    z_tf_inside_half = 0.0
    dz_tf_upper_lower_midplane = 0.0
    z_tf_top = 0.0
    hr1 = 0.0
    iohcl = 1
    i_cs_precomp = 1
    i_tf_inside_cs = 0
    dr_cs = 0.811
    dr_cs_precomp = 0.0
    rbld = 0.0
    required_radial_space = 0.0
    rinboard = 0.651
    r_shld_inboard_inner = 0.0
    r_shld_outboard_outer = 0.0
    r_vv_inboard_out = 0.0
    r_sh_inboard_out = 0.0
    r_tf_inboard_in = 0.0
    r_tf_inboard_mid = 0.0
    r_tf_inboard_out = 0.0
    r_tf_outboard_mid = 0.0
    i_r_cp_top = 0
    r_cp_top = 0.0
    f_r_cp = 1.4
    dr_tf_inner_bore = 0.0
    dh_tf_inner_bore = 0.0
    dr_fw_plasma_gap_inboard = 0.14
    dr_fw_plasma_gap_outboard = 0.15
    a_shld_total_surface = 0.0
    a_shld_inboard_surface = 0.0
    a_shld_outboard_surface = 0.0
    dr_shld_inboard = 0.69
    dz_shld_lower = 0.7
    dr_shld_outboard = 1.05
    dz_shld_upper = 0.6
    sigallpc = 3.0e8
    dr_tf_inboard = 0.0
    dz_tf_plasma_centre_offset = 0.0
    f_dr_tf_outboard_inboard = 1.19
    dr_tf_outboard = 0.0
    dr_tf_shld_gap = 0.05
    dr_shld_thermal_inboard = 0.05
    dr_shld_thermal_outboard = 0.05
    dz_shld_thermal = 0.05
    dz_shld_vv_gap = 0.163
    dz_xpoint_divertor = 0.0
    dz_fw_plasma_gap = 0.60
    dr_shld_blkt_gap = 0.05
    plleni = 1.0
    plleno = 1.0
    plsepi = 1.0
    plsepo = 1.5
    rspo = 0.0
    r_sh_inboard_in = 0.0
    z_plasma_xpoint_upper = 0.0
    z_plasma_xpoint_lower = 0.0
