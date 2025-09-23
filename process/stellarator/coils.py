import numpy as np

from process.data_structure import rebco_variables
import process.superconductors as superconductors
from process import process_output as po
from process.exceptions import ProcessValueError

from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    error_handling,
    fwbs_variables,
    physics_variables,
    sctfcoil_module,
    stellarator_configuration,
    stellarator_variables,
    tfcoil_variables,
    stellarator_module as st,
)


def st_coil(stellarator, output: bool):
    """Routine that performs the calculations for stellarator coils
    author: J Lion, IPP Greifswald
    outfile : input integer : output file unit
    iprint : input integer : switch for writing to output file (1=yes)
    This routine calculates the properties of the coils for
    a stellarator device.
    <P>Some precalculated effective parameters for a stellarator power
    plant design are used as the basis for the calculations. The coils
    are assumed to be a fixed shape, but are scaled in size
    appropriately for the machine being modelled.
    """
    r_coil_major = st.r_coil_major
    r_coil_minor = st.r_coil_minor

    ########################################################################################
    # Winding Pack Geometry: for one conductor
    #
    # This one conductor will just be multiplied later to fit the winding pack size.
    #
    # [m] Dimension of square cable space inside insulation
    #     and case of the conduit of each turn
    dx_tf_turn_cable_space_average = tfcoil_variables.t_turn_tf - 2.0e0 * (
        tfcoil_variables.dx_tf_turn_steel + tfcoil_variables.dx_tf_turn_insulation
    )  # dx_tf_turn_cable_space_average = t_w
    if dx_tf_turn_cable_space_average < 0:
        print(
            "dx_tf_turn_cable_space_average is negative. Check t_turn, tfcoil_variables.dx_tf_turn_steel and dx_tf_turn_insulation."
        )
    # [m^2] Cross-sectional area of cable space per turn
    tfcoil_variables.a_tf_turn_cable_space_no_void = (
        0.9e0 * dx_tf_turn_cable_space_average**2
    )  # 0.9 to include some rounded corners. (tfcoil_variables.a_tf_turn_cable_space_no_void = pi (dx_tf_turn_cable_space_average/2)**2 = pi/4 *dx_tf_turn_cable_space_average**2 for perfect round conductor). This factor depends on how round the corners are.
    # [m^2] Cross-sectional area of conduit case per turn
    tfcoil_variables.a_tf_turn_steel = (
        dx_tf_turn_cable_space_average + 2.0e0 * tfcoil_variables.dx_tf_turn_steel
    ) ** 2 - tfcoil_variables.a_tf_turn_cable_space_no_void
    #######################################################################################

    #######################################################################################
    # Winding Pack total size:
    #
    # Total coil current (MA)
    coilcurrent = (
        st.f_b * stellarator_configuration.stella_config_i0 * st.f_r / st.f_coil_aspect / st.f_n
    )
    st.f_i = coilcurrent / stellarator_configuration.stella_config_i0

    n_it = 200  # number of iterations

    rhs = np.zeros((n_it,))
    lhs = np.zeros((n_it,))
    jcrit_vector = np.zeros((n_it,))
    wp_width_r = np.zeros((n_it,))
    b_max_k = np.zeros((n_it,))

    for k in range(n_it):
        # Sample coil winding pack
        wp_width_r[k] = (r_coil_minor / 40.0e0) + (k / (n_it - 1e0)) * (
            r_coil_minor / 1.0e0 - r_coil_minor / 40.0e0
        )
        if tfcoil_variables.i_tf_sc_mat == 6:
            wp_width_r[k] = (r_coil_minor / 150.0e0) + (k / (n_it - 1e0)) * (
                r_coil_minor / 1.0e0 - r_coil_minor / 150.0e0
            )

        #  B-field calculation
        b_max_k[k] = bmax_from_awp(
            wp_width_r[k],
            coilcurrent,
            tfcoil_variables.n_tf_coils,
            r_coil_major,
            r_coil_minor,
        )
        # Two margins can be applied for jcrit: direct or by temperature margin.
        # Temperature margin is implemented in the jcrit_vector definition, 
        # direct margin is implemented after jcrit is defined (equation below)
        # jcrit for this bmax:
        jcrit_vector[k] = jcrit_from_material(
            b_max_k[k],
            tfcoil_variables.tftmp + tfcoil_variables.tmargmin,
            tfcoil_variables.i_tf_sc_mat,
            tfcoil_variables.b_crit_upper_nbti,
            tfcoil_variables.bcritsc,
            tfcoil_variables.fcutfsu,
            tfcoil_variables.fhts,
            tfcoil_variables.t_crit_nbti,
            tfcoil_variables.tcritsc,
            tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
            tfcoil_variables.j_tf_wp,
        )  # Get here a temperature margin from tfcoil_variables.tmargtf.

    # The operation current density weighted with the global iop/icrit fraction
    lhs[:] = constraint_variables.fiooic * jcrit_vector

    # Superconductor fraction in wp
    f_a_scu_of_wp = (
        (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        )
        * (1.0e0 - tfcoil_variables.fcutfsu)
        / (tfcoil_variables.t_turn_tf**2)

    )
    # print *, "f_a_scu_of_wp. ",f_a_scu_of_wp,"Awp min: ",Awp(1)

    rhs[:] = coilcurrent / (
        wp_width_r**2 / stellarator_configuration.stella_config_wp_ratio * f_a_scu_of_wp
    )  # f_a_scu_of_wp should be the fraction of the sc that is in the winding pack.

    wp_width_r_min = (
        r_coil_minor / 10.0e0
    ) ** 2  # Initial guess for intersection routine
    if tfcoil_variables.i_tf_sc_mat == 6:
        wp_width_r_min = (
            r_coil_minor / 20.0e0
        ) ** 2  # If REBCO, : start at smaller winding pack ratios

    # Find the intersection between LHS and RHS (or: how much awp do I need to get to the desired coil current)
    wp_width_r_min = intersect(
        wp_width_r, lhs, wp_width_r, rhs, wp_width_r_min
    )

    # Maximum field at superconductor surface (T)
    wp_width_r_min = max(tfcoil_variables.t_turn_tf**2, wp_width_r_min)

    # Recalculate tfcoil_variables.b_tf_inboard_peak at the found awp_min:
    tfcoil_variables.b_tf_inboard_peak = bmax_from_awp(
        wp_width_r_min,
        coilcurrent,
        tfcoil_variables.n_tf_coils,
        r_coil_major,
        r_coil_minor,
    )

    # Winding pack toroidal, radial cross-sections (m)
    awp_tor = (
        wp_width_r_min / stellarator_configuration.stella_config_wp_ratio
    )  # Toroidal dimension
    awp_rad = wp_width_r_min  # Radial dimension

    tfcoil_variables.dx_tf_wp_primary_toroidal = (
        awp_tor  # [m] toroidal thickness of winding pack
    )
    tfcoil_variables.dx_tf_wp_secondary_toroidal = (
        awp_tor  # [m] toroidal thickness of winding pack (region in front)
    )
    tfcoil_variables.dr_tf_wp_with_insulation = (
        awp_rad  # [m] radial thickness of winding pack
    )

    #  [m^2] winding-pack cross sectional area including insulation (not global)
    a_tf_wp_with_insulation = (
        tfcoil_variables.dr_tf_wp_with_insulation
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    ) * (
        tfcoil_variables.dx_tf_wp_primary_toroidal
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )

    a_tf_wp_no_insulation = (
        awp_tor * awp_rad
    )  # [m^2] winding-pack cross sectional area
    tfcoil_variables.j_tf_wp = (
        coilcurrent * 1.0e6 / a_tf_wp_no_insulation
    )  # [A/m^2] winding pack current density
    tfcoil_variables.n_tf_coil_turns = (
        a_tf_wp_no_insulation / (tfcoil_variables.t_turn_tf**2)
    )  # estimated number of turns for a given turn size (not global). Take at least 1.
    tfcoil_variables.c_tf_turn = (
        coilcurrent * 1.0e6 / tfcoil_variables.n_tf_coil_turns
    )  # [A] current per turn - estimation
    # [m^2] Total conductor cross-sectional area, taking account of void area
    tfcoil_variables.a_tf_wp_conductor = (
        tfcoil_variables.a_tf_turn_cable_space_no_void
        * tfcoil_variables.n_tf_coil_turns
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
    )
    # [m^2] Void area in cable, for He
    tfcoil_variables.a_tf_wp_extra_void = (
        tfcoil_variables.a_tf_turn_cable_space_no_void
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.f_a_tf_turn_cable_space_extra_void
    )
    # [m^2] Insulation area (not including ground-wall)
    tfcoil_variables.a_tf_coil_wp_turn_insulation = (
        tfcoil_variables.n_tf_coil_turns
        * (
            tfcoil_variables.t_turn_tf**2
            - tfcoil_variables.a_tf_turn_steel
            - tfcoil_variables.a_tf_turn_cable_space_no_void
        )
    )
    # [m^2] Structure area for cable
    tfcoil_variables.a_tf_wp_steel = (
        tfcoil_variables.n_tf_coil_turns * tfcoil_variables.a_tf_turn_steel
    )
    # End of winding pack calculations
    #######################################################################################

    #######################################################################################
    #  Casing calculations
    #
    # Coil case thickness (m). Here assumed to be constant
    # until something better comes up.
    # case_thickness_constant = tfcoil_variables.dr_tf_nose_case #0.2e0 # #? Leave this constant for now... Check this## Should be scaled with forces I think.
    #  For now assumed to be constant in a bolted plate model.
    #
    tfcoil_variables.dr_tf_plasma_case = (
        tfcoil_variables.dr_tf_nose_case
    )  # [m] coil case thickness outboard distance (radial)
    # dr_tf_nose_case = case_thickness_constant/2.0e0 # [m] coil case thickness inboard distance  (radial).
    tfcoil_variables.dx_tf_side_case_min = (
        tfcoil_variables.dr_tf_nose_case
    )  # [m] coil case thickness toroidal distance (toroidal)

    # End of casing calculations
    #######################################################################################

    #######################################################################################
    #  Port calculations
    #
    #  Maximal toroidal port size (vertical ports) (m)
    #  The maximal distance is correct but the vertical extension of this port is not clear#
    #  This is simplified for now and can be made more accurate in the future#
    stellarator_variables.vporttmax = (
        0.4e0
        * stellarator_configuration.stella_config_max_portsize_width
        * st.f_r
        / st.f_n
    )  # This is not accurate yet. Needs more insight#

    #  Maximal poloidal port size (vertical ports) (m)
    stellarator_variables.vportpmax = (
        2.0 * stellarator_variables.vporttmax
    )  # Simple approximation

    #  Maximal vertical port clearance area (m2)
    stellarator_variables.vportamax = (
        stellarator_variables.vporttmax * stellarator_variables.vportpmax
    )

    #  Horizontal ports
    #  Maximal toroidal port size (horizontal ports) (m)
    stellarator_variables.hporttmax = (
        0.8e0
        * stellarator_configuration.stella_config_max_portsize_width
        * st.f_r
        / st.f_n
    )  # Factor 0.8 to take the variation with height into account

    #  Maximal poloidal port size (horizontal ports) (m)
    stellarator_variables.hportpmax = (
        2.0e0 * stellarator_variables.hporttmax
    )  # Simple approximation

    #  Maximal horizontal port clearance area (m2)
    stellarator_variables.hportamax = (
        stellarator_variables.hporttmax * stellarator_variables.hportpmax
    )
    # End of port calculations
    #######################################################################################

    #######################################################################################
    #  General Coil Geometry values
    #
    tfcoil_variables.dx_tf_inboard_out_toroidal = (
        tfcoil_variables.dx_tf_wp_primary_toroidal
        + 2.0e0 * tfcoil_variables.dx_tf_side_case_min
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )  # [m] Thickness of inboard leg in toroidal direction

    build_variables.dr_tf_inboard = (
        tfcoil_variables.dr_tf_nose_case
        + tfcoil_variables.dr_tf_wp_with_insulation
        + tfcoil_variables.dr_tf_plasma_case
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )  # [m] Thickness of inboard leg in radial direction
    build_variables.dr_tf_outboard = (
        tfcoil_variables.dr_tf_nose_case
        + tfcoil_variables.dr_tf_wp_with_insulation
        + tfcoil_variables.dr_tf_plasma_case
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )  # [m] Thickness of outboard leg in radial direction (same as inboard)
    tfcoil_variables.a_tf_leg_outboard = (
        build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
    )  # [m^2] overall coil cross-sectional area (assuming inboard and
    #       outboard leg are the same)
    tfcoil_variables.a_tf_coil_inboard_case = (
        build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
    ) - a_tf_wp_with_insulation  # [m^2] Cross-sectional area of surrounding case

    tfcoil_variables.tfocrn = (
        0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal
    )  # [m] Half-width of side of coil nearest torus centreline
    tfcoil_variables.tficrn = (
        0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal
    )  # [m] Half-width of side of coil nearest plasma

    # [m^2] Total surface area of coil side facing plasma: inboard region
    tfcoil_variables.tfsai = (
        tfcoil_variables.n_tf_coils
        * tfcoil_variables.dx_tf_inboard_out_toroidal
        * 0.5e0
        * tfcoil_variables.len_tf_coil
    )
    # [m^2] Total surface area of coil side facing plasma: outboard region
    tfcoil_variables.tfsao = (
        tfcoil_variables.tfsai
    )  # depends, how 'inboard' and 'outboard' are defined

    # [m] Minimal distance in toroidal direction between two stellarator coils (from mid to mid)
    # Consistency with coil width is checked in constraint equation 82
    tfcoil_variables.toroidalgap = (
        stellarator_configuration.stella_config_dmin
        * (r_coil_major - r_coil_minor)
        / (
            stellarator_configuration.stella_config_coil_rmajor
            - stellarator_configuration.stella_config_coil_rminor
        )
    )
    # Left-Over coil gap between two coils (m)
    coilcoilgap = (
        tfcoil_variables.toroidalgap - tfcoil_variables.dx_tf_inboard_out_toroidal
    )

    #  Variables for ALL coils.
    tfcoil_variables.a_tf_inboard_total = (
        tfcoil_variables.n_tf_coils * tfcoil_variables.a_tf_leg_outboard
    )  # [m^2] Total area of all coil legs (midplane)
    tfcoil_variables.c_tf_total = (
        tfcoil_variables.n_tf_coils * coilcurrent * 1.0e6
    )  # [A] Total current in ALL coils
    tfcoil_variables.oacdcp = (
        tfcoil_variables.c_tf_total / tfcoil_variables.a_tf_inboard_total
    )  # [A / m^2] overall current density
    tfcoil_variables.r_b_tf_inboard_peak = (
        r_coil_major - r_coil_minor + awp_rad
    )  # [m] radius of peak field occurrence, average
    # jlion: not sure what this will be used for. Not very
    # useful for stellarators

    # This uses the reference value for the inductance and scales it with a^2/R (toroid inductance scaling)
    inductance = (
        stellarator_configuration.stella_config_inductance
        / st.f_r
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
        * st.f_n**2
    )
    tfcoil_variables.e_tf_magnetic_stored_total_gj = (
        0.5e0
        * (
            stellarator_configuration.stella_config_inductance
            / st.f_r
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
            ** 2
            * st.f_n**2
        )
        * (tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils) ** 2
        * 1.0e-9
    )  # [GJ] Total magnetic energy

    #  Coil dimensions
    build_variables.z_tf_inside_half = (
        0.5e0
        * stellarator_configuration.stella_config_maximal_coil_height
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
    )  # [m] maximum half-height of coil
    r_tf_inleg_mid = (
        r_coil_major - r_coil_minor
    )  # This is not very well defined for a stellarator.
    # Though, this is taken as an average value.
    tf_total_h_width = (
        r_coil_minor  # ? not really sure what this is supposed to be. Estimated as
    )
    # the average minor coil radius

    tfborev = (
        2.0e0 * build_variables.z_tf_inside_half
    )  # [m] estimated vertical coil dr_bore

    tfcoil_variables.len_tf_coil = (
        stellarator_configuration.stella_config_coillength
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
        / tfcoil_variables.n_tf_coils
    )  # [m] estimated average length of a coil

    # [m^2] Total surface area of toroidal shells covering coils
    tfcoil_variables.tfcryoarea = (
        stellarator_configuration.stella_config_coilsurface * st.f_r
        * (st.r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
        * 1.1e0
    )
    # 1.1 to scale it out a bit, as the shell must be bigger than WP


    # Minimal bending radius:
    min_bending_radius = (
        stellarator_configuration.stella_config_min_bend_radius
        * st.f_r
        * 1.0
        / (1.0 - tfcoil_variables.dr_tf_wp_with_insulation / (2.0 * r_coil_minor))
    )

    # End of general coil geometry values
    #######################################################################################

    #######################################################################################
    #  Masses of conductor constituents
    #
    # [kg] Mass of case
    #  (no need for correction factors as is the case for tokamaks)
    # This is only correct if the winding pack is 'thin' (len_tf_coil>>sqrt(tfcoil_variables.a_tf_coil_inboard_case)).
    tfcoil_variables.whtcas = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.a_tf_coil_inboard_case
        * tfcoil_variables.dcase
    )
    # Mass of ground-wall insulation [kg]
    # (assumed to be same density/material as conduit insulation)
    tfcoil_variables.whtgw = (
        tfcoil_variables.len_tf_coil
        * (a_tf_wp_with_insulation - a_tf_wp_no_insulation)
        * tfcoil_variables.dcondins
    )
    # [kg] mass of Superconductor
    tfcoil_variables.whtconsc = (
        (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coil_turns
            * tfcoil_variables.a_tf_turn_cable_space_no_void
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
            * (1.0e0 - tfcoil_variables.fcutfsu)
            - tfcoil_variables.len_tf_coil
            * tfcoil_variables.a_tf_wp_coolant_channels
        )
        * tfcoil_variables.dcond[tfcoil_variables.i_tf_sc_mat - 1]
    )  # a_tf_wp_coolant_channels is 0 for a stellarator. but keep this term for now.
    # [kg] mass of Copper in conductor
    tfcoil_variables.whtconcu = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.a_tf_turn_cable_space_no_void
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        * tfcoil_variables.fcutfsu
        - tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_wp_coolant_channels
    ) * constants.dcopper
    # [kg] mass of Steel conduit (sheath)
    tfcoil_variables.m_tf_turn_steel_conduit = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.a_tf_turn_steel
        * fwbs_variables.denstl
    )
    # if (i_tf_sc_mat==6)   tfcoil_variables.m_tf_turn_steel_conduit = fcondsteel * a_tf_wp_no_insulation *tfcoil_variables.len_tf_coil* fwbs_variables.denstl
    # Conduit insulation mass [kg]
    # (tfcoil_variables.a_tf_coil_wp_turn_insulation already contains tfcoil_variables.n_tf_coil_turns)
    tfcoil_variables.whtconin = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.a_tf_coil_wp_turn_insulation
        * tfcoil_variables.dcondins
    )
    # [kg] Total conductor mass
    tfcoil_variables.whtcon = (
        tfcoil_variables.whtconsc
        + tfcoil_variables.whtconcu
        + tfcoil_variables.m_tf_turn_steel_conduit
        + tfcoil_variables.whtconin
    )
    # [kg] Total coil mass
    tfcoil_variables.m_tf_coils_total = (
        tfcoil_variables.whtcas + tfcoil_variables.whtcon + tfcoil_variables.whtgw
    ) * tfcoil_variables.n_tf_coils
    # End of general coil geometry values
    #######################################################################################

    #######################################################################################
    # Quench protection:
    #
    # This copied from the tokamak module:
    # Radial position of vacuum vessel [m]
    rad_vv_in = (
        physics_variables.rmajor
        - physics_variables.rminor
        - build_variables.dr_fw_plasma_gap_inboard
        - build_variables.dr_fw_inboard
        - build_variables.dr_blkt_inboard
        - build_variables.dr_shld_blkt_gap
        - build_variables.dr_shld_inboard
    )
    rad_vv_out = (
        physics_variables.rmajor
        + physics_variables.rminor
        + build_variables.dr_fw_plasma_gap_outboard
        + build_variables.dr_fw_outboard
        + build_variables.dr_blkt_outboard
        + build_variables.dr_shld_blkt_gap
        + build_variables.dr_shld_outboard
    )

    # Stellarator version is working on the W7-X scaling, so we should actual use vv r_major
    # plasma r_major is just an approximation, but exact calculations require 3D geometry
    # Maybe it can be added to the stella_config file in the future
    rad_vv = physics_variables.rmajor

    # Actual VV force density
    # Based on reference values from W-7X:
    # Bref = 3;
    # Iref = 1.3*50;
    # aref = 0.92;
    # \[Tau]ref = 3.;
    # Rref = 5.2;
    # dref = 14*10^-3;

    # MN/m^3
    f_vv_actual = (
        2.54
        * (3e0 / physics_variables.bt
            * 1.3e6 * 50e0 / tfcoil_variables.c_tf_total
            * 0.92e0**2e0 / physics_variables.rminor**2
            ) **(-1)
        * (
                3e0 / tfcoil_variables.tdmptf
                * 5.2e0 / rad_vv
                * 0.014e0 / ((build_variables.dr_vv_inboard + build_variables.dr_vv_outboard) / 2)
            )
    )

    # This is not correct - it gives pressure on the vv wall, not stress
    # N/m^2
    # is the vv width the correct length to multiply by to turn the
    # force density into a stress?
    # sctfcoil_module.vv_stress_quench = (
    #     f_vv_actual
    #     * 1e6
    #     * ((build_variables.dr_vv_inboard + build_variables.dr_vv_outboard) / 2)
    # )

    # This approach merge stress model from tokamaks with induced force calculated from W7-X scaling
    a_vv = (rad_vv_out + rad_vv_in) / (rad_vv_out - rad_vv_in)
    zeta = 1 + ((a_vv - 1) * np.log((a_vv + 1) / (a_vv - 1)) / (2 * a_vv))

    sctfcoil_module.vv_stress_quench =  zeta * f_vv_actual * 1e6 * rad_vv_in

    # the conductor fraction is meant of the cable space#
    # This is the old routine which is being replaced for now by the new one below
    #    protect(aio,  tfes,               acs,       aturn,   tdump,  fcond,  fcu,   tba,  tmax   ,ajwpro, vd)
    # call protect(c_tf_turn,e_tf_magnetic_stored_total_gj/tfcoil_variables.n_tf_coils*1.0e9,a_tf_turn_cable_space_no_void,
    #    tfcoil_variables.t_turn_tf**2   ,tdmptf,1-f_a_tf_turn_cable_space_extra_void,fcutfsu,tftmp,tmaxpro,jwdgpro2,vd)

    vd = u_max_protect_v(
        tfcoil_variables.e_tf_magnetic_stored_total_gj
        / tfcoil_variables.n_tf_coils
        * 1.0e9,
        tfcoil_variables.tdmptf,
        tfcoil_variables.c_tf_turn,
    )

    # comparison
    # the new quench protection routine, see #1047
    tfcoil_variables.jwdgpro = calculate_quench_protection_current_density(
        tau_quench=tfcoil_variables.tdmptf,
        t_detect=tfcoil_variables.t_tf_quench_detection,
        fcu=tfcoil_variables.fcutfsu,
        fcond=1 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
        temp=tfcoil_variables.tftmp,
        acs=tfcoil_variables.a_tf_turn_cable_space_no_void,
        aturn=tfcoil_variables.t_turn_tf**2,
    )

    # Also give the copper area for REBCO quench calculations:
    rebco_variables.coppera_m2 = (
        coilcurrent
        * 1.0e6
        / (tfcoil_variables.a_tf_wp_conductor * tfcoil_variables.fcutfsu)
    )
    tfcoil_variables.vtfskv = vd / 1.0e3  # Dump voltage
    #
    #######################################################################################

    # Forces scaling #
    tfcoil_variables.max_force_density = (
        stellarator_configuration.stella_config_max_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )

    # Approximate, very simple maxiumum stress: (needed for limitation of icc 32)
    tfcoil_variables.sig_tf_wp = (
        tfcoil_variables.max_force_density
        * tfcoil_variables.dr_tf_wp_with_insulation
        * 1.0e6
    )  # in Pa

    # Units: MN/m
    max_force_density_mnm = (
        stellarator_configuration.stella_config_max_force_density_mnm
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
    )
    #
    max_lateral_force_density = (
        stellarator_configuration.stella_config_max_lateral_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )
    max_radial_force_density = (
        stellarator_configuration.stella_config_max_radial_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )
    #
    # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
    centering_force_max_mn = (
        stellarator_configuration.stella_config_centering_force_max_mn
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_coillength
        / tfcoil_variables.n_tf_coils
        / tfcoil_variables.len_tf_coil
    )
    centering_force_min_mn = (
        stellarator_configuration.stella_config_centering_force_min_mn
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_coillength
        / tfcoil_variables.n_tf_coils
        / tfcoil_variables.len_tf_coil
    )
    centering_force_avg_mn = (
        stellarator_configuration.stella_config_centering_force_avg_mn
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_coillength
        / tfcoil_variables.n_tf_coils
        / tfcoil_variables.len_tf_coil
    )
    #
    ####################################

    if output:
        print_output(
            stellarator,
            a_tf_wp_no_insulation,
            centering_force_avg_mn,
            centering_force_max_mn,
            centering_force_min_mn,
            coilcoilgap,
            rebco_variables.coppera_m2,
            rebco_variables.coppera_m2_max,
            f_a_scu_of_wp,
            f_vv_actual,
            constraint_variables.fiooic,
            inductance,
            tfcoil_variables.max_force_density,
            max_force_density_mnm,
            max_lateral_force_density,
            max_radial_force_density,
            min_bending_radius,
            r_coil_major,
            r_coil_minor,
            r_tf_inleg_mid,
            tfcoil_variables.sig_tf_wp,
            tfcoil_variables.t_turn_tf,
            tfcoil_variables.tdmptf,
            tf_total_h_width,
            tfborev,
            tfcoil_variables.toroidalgap,
            tfcoil_variables.vdalw,
            tfcoil_variables.vtfskv,
        )


def u_max_protect_v(tfes, tdump, aio):
    """tfes : input real : Energy stored in one TF coil (J)
    tdump : input real : Dump time (sec)
    aio : input real : Operating current (A)
    """
    return 2 * tfes / (tdump * aio)


def calculate_quench_protection_current_density(tau_quench, t_detect, fcu, fcond, temp, acs, aturn):
    """
    Calculates the current density limited by the protection limit.

    Simplified 0-D adiabatic heat balance "hotspot criterion" model.

    This is slightly diffrent that tokamak version (also diffrent from the stellarator paper). 
    We skip the superconduc6tor contribution (this should be more conservative in theory). 
    """
    temp_k = [4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 104, 114, 124]
    q_cu_array_sa2m4 = [
        1.08514e17,
        1.12043e17,
        1.12406e17,
        1.05940e17,
        9.49741e16,
        8.43757e16,
        7.56346e16,
        6.85924e16,
        6.28575e16,
        5.81004e16,
        5.40838e16,
        5.06414e16,
        4.76531e16,
    ]
    q_he_array_sa2m4 = [
        3.44562e16,
        9.92398e15,
        4.90462e15,
        2.41524e15,
        1.26368e15,
        7.51617e14,
        5.01632e14,
        3.63641e14,
        2.79164e14,
        2.23193e14,
        1.83832e14,
        1.54863e14,
        1.32773e14,
    ]

    q_he = np.interp(temp, temp_k, q_he_array_sa2m4)
    q_cu = np.interp(temp, temp_k, q_cu_array_sa2m4)

    # This leaves out the contribution from the superconductor fraction for now
    return (acs / aturn) * np.sqrt(
        1
        / (0.5 * tau_quench + t_detect)
        * (fcu**2 * fcond**2 * q_cu + fcu * fcond * (1 - fcond) * q_he)
    )


def jcrit_from_material(
    bmax,
    thelium,
    i_tf_sc_mat,
    b_crit_upper_nbti,
    bcritsc,
    fcutfsu,
    fhts,
    t_crit_nbti,
    tcritsc,
    f_a_tf_turn_cable_space_extra_void,
    jwp,
):
    strain = -0.005  # for now a small value
    f_he = f_a_tf_turn_cable_space_extra_void  # this is helium fraction in the superconductor (set it to the fixed global variable here)

    f_tf_conductor_copper = fcutfsu  # fcutfsu is a global variable. Is the copper fraction
    # of a cable conductor.

    if i_tf_sc_mat == 1:  # ITER Nb3Sn critical surface parameterization
        bc20m = 32.97  # these are values taken from sctfcoil.f90
        tc0m = 16.06

        #  j_crit_sc returned by itersc is the critical current density in the
        #  superconductor - not the whole strand, which contains copper
        if bmax > bc20m:
            j_crit_sc = 1.0e-9  # Set to a small nonzero value
        else:
            (
                j_crit_sc,
                bcrit,
                tcrit,
            ) = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)

        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0 - f_tf_conductor_copper) * (1.0e0 - f_he)

        # This is needed right now. Can we change it later?
        j_crit_sc = max(1.0e-9, j_crit_sc)
        j_crit_cable = max(1.0e-9, j_crit_cable)

    elif i_tf_sc_mat == 2:
        # Bi-2212 high temperature superconductor parameterization
        #  Current density in a strand of Bi-2212 conductor
        #  N.B. jcrit returned by bi2212 is the critical current density
        #  in the strand, not just the superconducting portion.
        #  The parameterization for j_crit_cable assumes a particular strand
        #  composition that does not require a user-defined copper fraction,
        #  so this is irrelevant in this model

        jstrand = jwp / (1 - f_he)
        #  jstrand = 0  # as far as I can tell this will always be 0
        #  because jwp was never set in fortran (so 0)

        j_crit_cable, tmarg = superconductors.bi2212(
            bmax, jstrand, thelium, fhts
        )  # bi2212 outputs j_crit_cable
        j_crit_sc = j_crit_cable / (1 - f_tf_conductor_copper)
        tcrit = thelium + tmarg
    elif i_tf_sc_mat == 3:  # NbTi data
        bc20m = 15.0
        tc0m = 9.3
        c0 = 1.0

        if bmax > bc20m:
            j_crit_sc = 1.0e-9  # Set to a small nonzero value
        else:
            j_crit_sc, tcrit = superconductors.jcrit_nbti(
                thelium,
                bmax,
                c0,
                bc20m,
                tc0m,
            )
            # I dont need tcrit here so dont use it.

        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)

        # This is needed right now. Can we change it later?
        j_crit_sc = max(1.0e-9, j_crit_sc)
        j_crit_cable = max(1.0e-9, j_crit_cable)
    elif i_tf_sc_mat == 4:  # As (1), but user-defined parameters
        bc20m = bcritsc
        tc0m = tcritsc
        j_crit_sc, bcrit, tcrit = superconductors.itersc(
            thelium, bmax, strain, bc20m, tc0m
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif i_tf_sc_mat == 5:  # WST Nb3Sn parameterisation
        bc20m = 32.97
        tc0m = 16.06

        #  j_crit_sc returned by itersc is the critical current density in the
        #  superconductor - not the whole strand, which contains copper

        j_crit_sc, bcrit, tcrit = superconductors.western_superconducting_nb3sn(
            thelium,
            bmax,
            strain,
            bc20m,
            tc0m,
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif (
        i_tf_sc_mat == 6
    ):  # ! "REBCO" 2nd generation HTS superconductor in CrCo strand
        j_crit_sc, validity = superconductors.jcrit_rebco(thelium, bmax, 0)
        j_crit_sc = max(1.0e-9, j_crit_sc)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)

    elif i_tf_sc_mat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
        bc20m = b_crit_upper_nbti
        tc0m = t_crit_nbti
        j_crit_sc, bcrit, tcrit = superconductors.gl_nbti(
            thelium, bmax, strain, bc20m, tc0m
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif i_tf_sc_mat == 8:
        bc20m = 429
        tc0m = 185
        j_crit_sc, bcrit, tcrit = superconductors.gl_rebco(
            thelium, bmax, strain, bc20m, tc0m
        )
        # A0 calculated for tape cross section already
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    else:
        raise ProcessValueError(
            "Illegal value for i_pf_superconductor", i_tf_sc_mat=i_tf_sc_mat
        )

    return j_crit_sc * 1e-6


def intersect(x1, y1, x2, y2, xin):
        """Routine to find the x (abscissa) intersection point of two curves
        each defined by tabulated (x,y) values
        author: P J Knight, CCFE, Culham Science Centre
        x1(1:n1) : input real array : x values for first curve
        y1(1:n1) : input real array : y values for first curve
        n1       : input integer : length of arrays x1, y1
        x2(1:n2) : input real array : x values for first curve
        y2(1:n2) : input real array : y values for first curve
        n2       : input integer : length of arrays x2, y2
        x        : input/output real : initial x value guess on entry;
        x value at point of intersection on exit
        This routine estimates the x point (abscissa) at which two curves
        defined by tabulated (x,y) values intersect, using simple
        linear interpolation and the Newton-Raphson method.
        The routine will stop with an error message if no crossing point
        is found within the x ranges of the two curves.
        None
        """
        x = xin
        n1 = len(x1)
        n2 = len(x2)

        xmin = max(np.amin(x1), np.amin(x2))
        xmax = min(np.max(x1), np.amax(x2))

        if xmin >= xmax:
            error_handling.fdiags[0] = np.amin(x1)
            error_handling.fdiags[1] = np.amin(x2)
            error_handling.fdiags[2] = np.amax(x1)
            error_handling.fdiags[3] = np.amax(x2)
            error_handling.report_error(111)

        #  Ensure input guess for x is within this range

        if x < xmin:
            x = xmin
        elif x > xmax:
            x = xmax

        #  Find overall y range, and set tolerance
        #  in final difference in y values

        ymin = min(np.amin(y1), np.amin(y2))
        ymax = max(np.max(y1), np.max(y2))

        epsy = 1.0e-6 * (ymax - ymin)

        #  Finite difference dx

        dx = 0.01e0 / max(n1, n2) * (xmax - xmin)

        for _i in range(100):
            #  Find difference in y values at x

            y01 = np.interp(x, x1, y1)
            y02 = np.interp(x, x2, y2)
            y = y01 - y02

            if abs(y) < epsy:
                break

            #  Find difference in y values at x+dx

            y01 = np.interp(x + dx, x1, y1)
            y02 = np.interp(x + dx, x2, y2)
            yright = y01 - y02

            #  Find difference in y values at x-dx

            y01 = np.interp(x - dx, x1, y1)
            y02 = np.interp(x - dx, x2, y2)
            yleft = y01 - y02

            #  Adjust x using Newton-Raphson method

            x = x - 2.0e0 * dx * y / (yright - yleft)

            if x < xmin:
                error_handling.fdiags[0] = x
                error_handling.fdiags[1] = xmin
                error_handling.report_error(112)
                x = xmin
                break

            if x > xmax:
                error_handling.fdiags[0] = x
                error_handling.fdiags[1] = xmax
                error_handling.report_error(113)
                x = xmax
                break
        else:
            error_handling.report_error(114)

        return x


def bmax_from_awp(
    wp_width_radial, current, n_tf_coils, r_coil_major, r_coil_minor
):
    """Returns a fitted function for bmax for stellarators

    author: J Lion, IPP Greifswald
    Returns a fitted function for bmax in dependece
    of the winding pack. The stellarator type config
    is taken from the parent scope.
    """

    return (
        2e-1
        * current
        * n_tf_coils
        / (r_coil_major - r_coil_minor)
        * (
            stellarator_configuration.stella_config_a1
            + stellarator_configuration.stella_config_a2
            * r_coil_major
            / wp_width_radial
        )
    )


def print_output(
        stellarator,
        a_tf_wp_no_insulation,
        centering_force_avg_mn,
        centering_force_max_mn,
        centering_force_min_mn,
        coilcoilgap,
        coppera_m2,
        coppera_m2_max,
        f_a_scu_of_wp,
        f_vv_actual,
        fiooic,
        inductance,
        max_force_density,
        max_force_density_mnm,
        max_lateral_force_density,
        max_radial_force_density,
        min_bending_radius,
        r_coil_major,
        r_coil_minor,
        r_tf_inleg_mid,
        sig_tf_wp,
        t_turn_tf,
        tdmptf,
        tf_total_h_width,
        tfborev,
        toroidalgap,
        vdalw,
        vtfskv,
    ):
        """Writes stellarator modular coil output to file
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        This routine writes the stellarator modular coil results
        to the output file.
        None
        """
        po.oheadr(stellarator.outfile, "Modular Coils")

        po.osubhd(stellarator.outfile, "General Coil Parameters :")

        po.ovarre(
            stellarator.outfile,
            "Number of modular coils",
            "(n_tf_coils)",
            tfcoil_variables.n_tf_coils,
        )
        po.ovarre(stellarator.outfile, "Av. coil major radius", "(coil_r)", r_coil_major)
        po.ovarre(stellarator.outfile, "Av. coil minor radius", "(coil_a)", r_coil_minor)
        po.ovarre(
            stellarator.outfile,
            "Av. coil aspect ratio",
            "(coil_aspect)",
            r_coil_major / r_coil_minor,
        )

        po.ovarre(
            stellarator.outfile,
            "Cross-sectional area per coil (m2)",
            "(tfarea/n_tf_coils)",
            tfcoil_variables.a_tf_inboard_total / tfcoil_variables.n_tf_coils,
        )
        po.ovarre(
            stellarator.outfile,
            "Total inboard leg radial thickness (m)",
            "(dr_tf_inboard)",
            build_variables.dr_tf_inboard,
        )
        po.ovarre(
            stellarator.outfile,
            "Total outboard leg radial thickness (m)",
            "(dr_tf_outboard)",
            build_variables.dr_tf_outboard,
        )
        po.ovarre(
            stellarator.outfile,
            "Inboard leg outboard half-width (m)",
            "(tficrn)",
            tfcoil_variables.tficrn,
        )
        po.ovarre(
            stellarator.outfile,
            "Inboard leg inboard half-width (m)",
            "(tfocrn)",
            tfcoil_variables.tfocrn,
        )
        po.ovarre(
            stellarator.outfile,
            "Outboard leg toroidal thickness (m)",
            "(dx_tf_inboard_out_toroidal)",
            tfcoil_variables.dx_tf_inboard_out_toroidal,
        )
        po.ovarre(
            stellarator.outfile, "Minimum coil distance (m)", "(toroidalgap)", toroidalgap
        )
        po.ovarre(
            stellarator.outfile,
            "Minimal left gap between coils (m)",
            "(coilcoilgap)",
            coilcoilgap,
        )
        po.ovarre(
            stellarator.outfile,
            "Minimum coil bending radius (m)",
            "(min_bend_radius)",
            min_bending_radius,
        )
        po.ovarre(
            stellarator.outfile,
            "Mean coil circumference (m)",
            "(len_tf_coil)",
            tfcoil_variables.len_tf_coil,
        )
        po.ovarre(
            stellarator.outfile,
            "Total current (MA)",
            "(c_tf_total)",
            1.0e-6 * tfcoil_variables.c_tf_total,
        )
        po.ovarre(
            stellarator.outfile,
            "Current per coil(MA)",
            "(c_tf_total/n_tf_coils)",
            1.0e-6 * tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils,
        )
        po.ovarre(
            stellarator.outfile,
            "Winding pack current density (A/m2)",
            "(j_tf_wp)",
            tfcoil_variables.j_tf_wp,
        )
        po.ovarre(
            stellarator.outfile,
            "Max allowable current density as restricted by quench (A/m2)",
            "(jwdgpro)",
            tfcoil_variables.jwdgpro,
        )
        po.ovarre(
            stellarator.outfile,
            "Overall current density (A/m2)",
            "(oacdcp)",
            tfcoil_variables.oacdcp,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximum field on superconductor (T)",
            "(b_tf_inboard_peak)",
            tfcoil_variables.b_tf_inboard_peak,
        )
        po.ovarre(
            stellarator.outfile,
            "Total Stored energy (GJ)",
            "(e_tf_magnetic_stored_total_gj)",
            tfcoil_variables.e_tf_magnetic_stored_total_gj,
        )
        po.ovarre(
            stellarator.outfile, "Inductance of TF Coils (H)", "(inductance)", inductance
        )
        po.ovarre(
            stellarator.outfile,
            "Total mass of coils (kg)",
            "(m_tf_coils_total)",
            tfcoil_variables.m_tf_coils_total,
        )

        po.osubhd(stellarator.outfile, "Coil Geometry :")
        po.ovarre(
            stellarator.outfile,
            "Inboard leg centre radius (m)",
            "(r_tf_inleg_mid)",
            r_tf_inleg_mid,
        )
        po.ovarre(
            stellarator.outfile,
            "Outboard leg centre radius (m)",
            "(r_tf_outboard_mid)",
            build_variables.r_tf_outboard_mid,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximum inboard edge height (m)",
            "(z_tf_inside_half)",
            build_variables.z_tf_inside_half,
        )
        po.ovarre(
            stellarator.outfile,
            "Clear horizontal dr_bore (m)",
            "(tf_total_h_width)",
            tf_total_h_width,
        )
        po.ovarre(stellarator.outfile, "Clear vertical dr_bore (m)", "(tfborev)", tfborev)

        po.osubhd(stellarator.outfile, "Conductor Information :")
        po.ovarre(
            stellarator.outfile,
            "Superconductor mass per coil (kg)",
            "(whtconsc)",
            tfcoil_variables.whtconsc,
        )
        po.ovarre(
            stellarator.outfile,
            "Copper mass per coil (kg)",
            "(whtconcu)",
            tfcoil_variables.whtconcu,
        )
        po.ovarre(
            stellarator.outfile,
            "Steel conduit mass per coil (kg)",
            "(m_tf_turn_steel_conduit)",
            tfcoil_variables.m_tf_turn_steel_conduit,
        )
        po.ovarre(
            stellarator.outfile,
            "Total conductor cable mass per coil (kg)",
            "(whtcon)",
            tfcoil_variables.whtcon,
        )
        po.ovarre(
            stellarator.outfile,
            "Cable conductor + void area (m2)",
            "(a_tf_turn_cable_space_no_void)",
            tfcoil_variables.a_tf_turn_cable_space_no_void,
        )
        po.ovarre(
            stellarator.outfile,
            "Cable space coolant fraction",
            "(f_a_tf_turn_cable_space_extra_void)",
            tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
        )
        po.ovarre(
            stellarator.outfile,
            "Conduit case thickness (m)",
            "(dx_tf_turn_steel)",
            tfcoil_variables.dx_tf_turn_steel,
        )
        po.ovarre(
            stellarator.outfile,
            "Cable insulation thickness (m)",
            "(dx_tf_turn_insulation)",
            tfcoil_variables.dx_tf_turn_insulation,
        )

        ap = a_tf_wp_no_insulation
        po.osubhd(stellarator.outfile, "Winding Pack Information :")
        po.ovarre(stellarator.outfile, "Winding pack area", "(ap)", ap)
        po.ovarre(
            stellarator.outfile,
            "Conductor fraction of winding pack",
            "(a_tf_wp_conductor/ap)",
            tfcoil_variables.a_tf_wp_conductor / ap,
        )
        po.ovarre(
            stellarator.outfile,
            "Copper fraction of conductor",
            "(fcutfsu)",
            tfcoil_variables.fcutfsu,
        )
        po.ovarre(
            stellarator.outfile,
            "Structure fraction of winding pack",
            "(a_tf_wp_steel/ap)",
            tfcoil_variables.a_tf_wp_steel / ap,
        )
        po.ovarre(
            stellarator.outfile,
            "Insulator fraction of winding pack",
            "(a_tf_coil_wp_turn_insulation/ap)",
            tfcoil_variables.a_tf_coil_wp_turn_insulation / ap,
        )
        po.ovarre(
            stellarator.outfile,
            "Helium fraction of winding pack",
            "(a_tf_wp_extra_void/ap)",
            tfcoil_variables.a_tf_wp_extra_void / ap,
        )
        po.ovarre(
            stellarator.outfile,
            "Winding radial thickness (m)",
            "(dr_tf_wp_with_insulation)",
            tfcoil_variables.dr_tf_wp_with_insulation,
        )
        po.ovarre(
            stellarator.outfile,
            "Winding toroidal thickness (m)",
            "(dx_tf_wp_primary_toroidal)",
            tfcoil_variables.dx_tf_wp_primary_toroidal,
        )
        po.ovarre(
            stellarator.outfile,
            "Ground wall insulation thickness (m)",
            "(dx_tf_wp_insulation)",
            tfcoil_variables.dx_tf_wp_insulation,
        )
        po.ovarre(
            stellarator.outfile,
            "Number of turns per coil",
            "(n_tf_coil_turns)",
            tfcoil_variables.n_tf_coil_turns,
        )
        po.ovarre(
            stellarator.outfile,
            "Width of each turn (incl. insulation) (m)",
            "(t_turn_tf)",
            t_turn_tf,
        )
        po.ovarre(
            stellarator.outfile,
            "Current per turn (A)",
            "(c_tf_turn)",
            tfcoil_variables.c_tf_turn,
        )
        po.ovarre(stellarator.outfile, "jop/jcrit", "(fiooic)", fiooic)
        po.ovarre(
            stellarator.outfile,
            "Current density in conductor area (A/m2)",
            "(c_tf_total/a_tf_wp_conductor)",
            1.0e-6
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.a_tf_wp_conductor,
        )
        po.ovarre(
            stellarator.outfile,
            "Current density in SC area (A/m2)",
            "(c_tf_total/a_tf_wp_conductor/f_a_scu_of_wp)",
            1.0e-6
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
            / ap
            / f_a_scu_of_wp,
        )
        po.ovarre(stellarator.outfile, "Superconductor faction of WP (1)", "(f_a_scu_of_wp)", f_a_scu_of_wp)

        po.osubhd(stellarator.outfile, "Forces and Stress :")
        po.ovarre(
            stellarator.outfile,
            "Maximal toroidally and radially av. force density (MN/m3)",
            "(max_force_density)",
            max_force_density,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximal force density (MN/m)",
            "(max_force_density_Mnm)",
            max_force_density_mnm,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximal stress (approx.) (MPa)",
            "(sig_tf_wp)",
            sig_tf_wp * 1.0e-6,
        )

        po.ovarre(
            stellarator.outfile,
            "Maximal lateral force density (MN/m3)",
            "(max_lateral_force_density)",
            max_lateral_force_density,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximal radial force density (MN/m3)",
            "(max_radial_force_density)",
            max_radial_force_density,
        )

        po.ovarre(
            stellarator.outfile,
            "Max. centering force (coil) (MN)",
            "(centering_force_max_MN)",
            centering_force_max_mn,
        )
        po.ovarre(
            stellarator.outfile,
            "Min. centering force (coil) (MN)",
            "(centering_force_min_MN)",
            centering_force_min_mn,
        )
        po.ovarre(
            stellarator.outfile,
            "Avg. centering force per coil (MN)",
            "(centering_force_avg_MN)",
            centering_force_avg_mn,
        )

        po.osubhd(stellarator.outfile, "Quench Restrictions :")
        po.ovarre(
            stellarator.outfile,
            "Actual quench time (or time constant) (s)",
            "(tdmptf)",
            tdmptf,
        )
        po.ovarre(
            stellarator.outfile,
            "Actual quench vaccuum vessel force density (MN/m^3)",
            "(f_vv_actual)",
            f_vv_actual,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximum allowed voltage during quench due to insulation (kV)",
            "(vdalw)",
            vdalw,
        )
        po.ovarre(stellarator.outfile, "Actual quench voltage (kV)", "(vtfskv)", vtfskv, "OP ")
        po.ovarre(
            stellarator.outfile,
            "Current (A) per mm^2 copper (A/mm2)",
            "(coppera_m2)",
            coppera_m2 * 1.0e-6,
        )
        po.ovarre(
            stellarator.outfile,
            "Max Copper current fraction:",
            "(coppera_m2/coppera_m2_max)",
            coppera_m2 / coppera_m2_max,
        )

        po.osubhd(stellarator.outfile, "External Case Information :")

        po.ovarre(
            stellarator.outfile,
            "Case thickness, plasma side (m)",
            "(dr_tf_plasma_case)",
            tfcoil_variables.dr_tf_plasma_case,
        )
        po.ovarre(
            stellarator.outfile,
            "Case thickness, outer side (m)",
            "(dr_tf_nose_case)",
            tfcoil_variables.dr_tf_nose_case,
        )
        po.ovarre(
            stellarator.outfile,
            "Case toroidal thickness (m)",
            "(dx_tf_side_case_min)",
            tfcoil_variables.dx_tf_side_case_min,
        )
        po.ovarre(
            stellarator.outfile,
            "Case area per coil (m2)",
            "(a_tf_coil_inboard_case)",
            tfcoil_variables.a_tf_coil_inboard_case,
        )
        po.ovarre(
            stellarator.outfile,
            "External case mass per coil (kg)",
            "(whtcas)",
            tfcoil_variables.whtcas,
        )

        po.osubhd(stellarator.outfile, "Available Space for Ports :")

        po.ovarre(
            stellarator.outfile,
            "Max toroidal size of vertical ports (m)",
            "(vporttmax)",
            stellarator_variables.vporttmax,
        )
        po.ovarre(
            stellarator.outfile,
            "Max poloidal size of vertical ports (m)",
            "(vportpmax)",
            stellarator_variables.vportpmax,
        )
        po.ovarre(
            stellarator.outfile,
            "Max area of vertical ports (m2)",
            "(vportamax)",
            stellarator_variables.vportamax,
        )
        po.ovarre(
            stellarator.outfile,
            "Max toroidal size of horizontal ports (m)",
            "(hporttmax)",
            stellarator_variables.hporttmax,
        )
        po.ovarre(
            stellarator.outfile,
            "Max poloidal size of horizontal ports (m)",
            "(hportpmax)",
            stellarator_variables.hportpmax,
        )
        po.ovarre(
            stellarator.outfile,
            "Max area of horizontal ports (m2)",
            "(hportamax)",
            stellarator_variables.hportamax,
        )