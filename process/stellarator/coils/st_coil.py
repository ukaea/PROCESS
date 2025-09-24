from process.data_structure import rebco_variables
import process.stellarator.coils.forces as forces
from process.stellarator.coils.coils import bmax_from_awp, calculate_quench_protection_current_density, intersect, jcrit_from_material, max_dump_voltage

from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    fwbs_variables,
    sctfcoil_module,
    physics_variables,
    stellarator_configuration,
    stellarator_variables,
    tfcoil_variables,
)
from process.fortran import (
    stellarator_module as st,
)

import numpy as np

from process.stellarator.coils.output import write


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

    #######################################################################################
    winding_pack_geometry()

    coilcurrent, awp_rad, a_tf_wp_no_insulation, a_tf_wp_with_insulation, f_a_scu_of_wp = winding_pack_total_size(r_coil_major, r_coil_minor)
    # End of winding pack calculations

    #######################################################################################
    #  Casing calculations
    calculate_casing()

    #######################################################################################
    #  Port calculations
    vertical_ports()
    horizontal_ports()

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


    # comparison
    # the new quench protection routine, see #1047
    tfcoil_variables.jwdgpro = calculate_quench_protection_current_density(
        tau_quench=tfcoil_variables.tdmptf,
        t_detect=tfcoil_variables.t_tf_quench_detection,
        f_cu=tfcoil_variables.fcutfsu,
        f_cond=1 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
        temp=tfcoil_variables.tftmp,
        a_cable=tfcoil_variables.a_tf_turn_cable_space_no_void,
        a_turn=tfcoil_variables.t_turn_tf**2,
    )

    # Also give the copper area for REBCO quench calculations:
    rebco_variables.coppera_m2 = (
        coilcurrent
        * 1.0e6
        / (tfcoil_variables.a_tf_wp_conductor * tfcoil_variables.fcutfsu)
    )

    # Max volatage during fast discharge of TF coil (V)
    # (note that tf_coil_variable is in kV, while calculation is in V)
    tfcoil_variables.vtfskv = max_dump_voltage(
        tfcoil_variables.e_tf_magnetic_stored_total_gj
        / tfcoil_variables.n_tf_coils
        * 1.0e9,
        tfcoil_variables.tdmptf,
        tfcoil_variables.c_tf_turn,
    ) / 1.0e3  
    #
    #######################################################################################

    # Forces scaling #
    forces.calculate_max_force_density(a_tf_wp_no_insulation)
    forces.calculate_maximum_stress()

    # Units: MN/m
    max_force_density_mnm = (
        stellarator_configuration.stella_config_max_force_density_mnm
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
    )
    #
    max_lateral_force_density = forces.calculate_max_lateral_force_density(a_tf_wp_no_insulation)
    max_radial_force_density = forces.calculate_max_radial_force_density(a_tf_wp_no_insulation)
    #
    # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
    centering_force_max_mn = forces.calculate_centering_force_max_mn()
    centering_force_min_mn = forces.calculate_centering_force_min_mn()
    centering_force_avg_mn = forces.calculate_centering_force_avg_mn()
    #
    ####################################

    if output:
        write(
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


def winding_pack_geometry():
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


def winding_pack_total_size(r_coil_major, r_coil_minor):
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

    return coilcurrent, awp_rad, a_tf_wp_no_insulation, a_tf_wp_with_insulation, f_a_scu_of_wp


def calculate_casing():
    """    
    Coil case thickness (m). Here assumed to be constant until something better comes up.
    case_thickness_constant = tfcoil_variables.dr_tf_nose_case #0.2e0 # ?
    # Leave this constant for now... Check this## Should be scaled with forces I think.
    For now assumed to be constant in a bolted plate model.
    """
    # [m] coil case thickness outboard distance (radial)
    tfcoil_variables.dr_tf_plasma_case = (
        tfcoil_variables.dr_tf_nose_case
    )  
    # dr_tf_nose_case = case_thickness_constant/2.0e0 # [m] coil case thickness inboard distance  (radial).

    # [m] coil case thickness toroidal distance (toroidal)
    tfcoil_variables.dx_tf_side_case_min = (
        tfcoil_variables.dr_tf_nose_case
    )  


def vertical_ports():
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

def horizontal_ports():
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