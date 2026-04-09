import logging

import numpy as np

import process.models.stellarator.coils.forces as forces
from process.data_structure import (
    build_variables,
    constraint_variables,
    rebco_variables,
    stellarator_configuration,
    stellarator_variables,
    tfcoil_variables,
)
from process.models.stellarator.coils.coils import (
    bmax_from_awp,
    intersect,
    jcrit_from_material,
)
from process.models.stellarator.coils.mass import calculate_coils_mass
from process.models.stellarator.coils.output import write
from process.models.stellarator.coils.quench import calculate_quench_protection

logger = logging.getLogger(__name__)


def st_coil(stellarator, output: bool):
    """This routine calculates the properties of the coils for
    a stellarator device.

    Some precalculated effective parameters for a stellarator power
    plant design are used as the basis for the calculations. The coils
    are assumed to be a fixed shape, but are scaled in size
    appropriately for the machine being modelled.

    Parameters
    ----------
    stellarator :

    output:


    """
    r_coil_major = stellarator_variables.r_coil_major
    r_coil_minor = stellarator_variables.r_coil_minor

    #######################################################################################
    calculate_winding_pack_geometry()

    # Total coil current (MA)
    coilcurrent = calculate_current()

    awp_rad, a_tf_wp_no_insulation, a_tf_wp_with_insulation, f_a_scu_of_wp = (
        winding_pack_total_size(r_coil_major, r_coil_minor, coilcurrent)
    )

    #######################################################################################
    #  Casing calculations
    calculate_casing()

    #######################################################################################
    #  Port calculations
    calculate_vertical_ports()
    calculate_horizontal_ports()

    #######################################################################################
    #  General Coil Geometry values
    #
    calculate_coil_toroidal_thickness()
    calculate_coil_radial_thickness()

    calculate_coil_cross_sectional_area(a_tf_wp_with_insulation)

    calculate_coil_half_widths()

    calculate_plasma_facing_coil_area()

    coil_coil_gap, _ = calculate_coil_coil_toroidal_gap(r_coil_major, r_coil_minor)

    calculate_coils_summary_variables(coilcurrent, r_coil_major, r_coil_minor, awp_rad)

    inductance = calculate_inductnace(r_coil_minor)
    calculate_stored_magnetic_energy(r_coil_minor)

    #  Coil dimensions
    build_variables.z_tf_inside_half = (
        0.5e0
        * stellarator_configuration.stella_config_maximal_coil_height
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
    )  # [m] maximum half-height of coil

    # [m] estimated average length of a coil
    tfcoil_variables.len_tf_coil = (
        stellarator_configuration.stella_config_coillength
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
        / tfcoil_variables.n_tf_coils
    )

    # [m^2] Total surface area of toroidal shells covering coils
    tfcoil_variables.tfcryoarea = (
        stellarator_configuration.stella_config_coilsurface
        * stellarator_variables.f_st_rmajor
        * (
            stellarator_variables.r_coil_minor
            / stellarator_configuration.stella_config_coil_rminor
        )
        * 1.1e0
    )
    # 1.1 to scale it out a bit, as the shell must be bigger than WP

    # Minimal bending radius:
    min_bending_radius = (
        stellarator_configuration.stella_config_min_bend_radius
        * stellarator_variables.f_st_rmajor
        / (1.0 - tfcoil_variables.dr_tf_wp_with_insulation / (2.0 * r_coil_minor))
    )

    # End of general coil geometry values

    #######################################################################################
    #  Coil_mases calculations
    calculate_coils_mass(a_tf_wp_with_insulation, a_tf_wp_no_insulation)

    #######################################################################################
    # Quench protection:
    f_vv_actual = calculate_quench_protection(coilcurrent)

    #
    #######################################################################################
    # Forces scaling #
    forces.calculate_max_force_density(a_tf_wp_no_insulation)
    forces.calculate_maximum_stress()

    # Units: MN/m
    max_force_density_mnm = forces.calculate_max_force_density_mnm()
    #
    max_lateral_force_density = forces.calculate_max_lateral_force_density(
        a_tf_wp_no_insulation
    )
    max_radial_force_density = forces.calculate_max_radial_force_density(
        a_tf_wp_no_insulation
    )
    #
    # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
    centering_force_max_mn = forces.calculate_centering_force_max_mn()
    centering_force_min_mn = forces.calculate_centering_force_min_mn()
    centering_force_avg_mn = forces.calculate_centering_force_avg_mn()
    #
    ####################################

    if output:
        write(
            stellarator=stellarator,
            a_tf_wp_no_insulation=a_tf_wp_no_insulation,
            centering_force_avg_mn=centering_force_avg_mn,
            centering_force_max_mn=centering_force_max_mn,
            centering_force_min_mn=centering_force_min_mn,
            coilcoilgap=coil_coil_gap,
            coppera_m2=rebco_variables.coppera_m2,
            coppera_m2_max=rebco_variables.coppera_m2_max,
            f_a_scu_of_wp=f_a_scu_of_wp,
            f_vv_actual=f_vv_actual,
            fiooic=constraint_variables.fiooic,
            inductance=inductance,
            max_force_density=tfcoil_variables.max_force_density,
            max_force_density_mnm=max_force_density_mnm,
            max_lateral_force_density=max_lateral_force_density,
            max_radial_force_density=max_radial_force_density,
            min_bending_radius=min_bending_radius,
            r_coil_major=r_coil_major,
            r_coil_minor=r_coil_minor,
            sig_tf_wp=tfcoil_variables.sig_tf_wp,
            dx_tf_turn_general=tfcoil_variables.dx_tf_turn_general,
            t_tf_superconductor_quench=tfcoil_variables.t_tf_superconductor_quench,
            toroidalgap=tfcoil_variables.toroidalgap,
            allowed_quench_voltage=tfcoil_variables.v_tf_coil_dump_quench_max_kv,
            quench_voltage=tfcoil_variables.v_tf_coil_dump_quench_kv,
        )


def calculate_coil_toroidal_thickness():
    tfcoil_variables.dx_tf_inboard_out_toroidal = (
        tfcoil_variables.dx_tf_wp_primary_toroidal
        + 2.0e0 * tfcoil_variables.dx_tf_side_case_min
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )  # [m] Thickness of inboard leg in toroidal direction


def calculate_coil_radial_thickness():
    """Thickness of inboard and outboard leg in radial direction"""
    # [m] Thickness of inboard leg in radial direction
    build_variables.dr_tf_inboard = (
        tfcoil_variables.dr_tf_nose_case
        + tfcoil_variables.dr_tf_wp_with_insulation
        + tfcoil_variables.dr_tf_plasma_case
        + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
    )
    # [m] Thickness of outboard leg in radial direction (same as inboard)
    build_variables.dr_tf_outboard = build_variables.dr_tf_inboard


def calculate_coil_cross_sectional_area(a_tf_wp_with_insulation):
    # [m^2] overall coil cross-sectional area
    # (assuming inboard and outboard leg are the same)
    tfcoil_variables.a_tf_leg_outboard = (
        build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
    )
    # [m^2] Cross-sectional area of surrounding case
    tfcoil_variables.a_tf_coil_inboard_case = (
        build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
    ) - a_tf_wp_with_insulation


def calculate_coil_half_widths():
    # [m] Half-width of side of coil nearest torus centreline
    tfcoil_variables.tfocrn = 0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal
    # [m] Half-width of side of coil nearest plasma
    tfcoil_variables.tficrn = 0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal


def calculate_plasma_facing_coil_area():
    # [m^2] Total surface area of coil side facing plasma: inboard region
    tfcoil_variables.tfsai = (
        tfcoil_variables.n_tf_coils
        * tfcoil_variables.dx_tf_inboard_out_toroidal
        * 0.5e0
        * tfcoil_variables.len_tf_coil
    )
    # [m^2] Total surface area of coil side facing plasma: outboard region (same as inboard)
    tfcoil_variables.tfsao = tfcoil_variables.tfsai


def calculate_coil_coil_toroidal_gap(r_coil_major, r_coil_minor):
    """[m] Minimal distance in toroidal direction between two stellarator coils
    Consistency with coil width is checked in constraint equation 82

    Parameters
    ----------
    r_coil_major :

    r_coil_minor :

    """
    # [m] Toroidal gap between two coil filaments
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
    return coilcoilgap, tfcoil_variables.toroidalgap


def calculate_coils_summary_variables(coilcurrent, r_coil_major, r_coil_minor, awp_rad):
    """Variables for ALL coils.

    Parameters
    ----------
    coilcurrent :

    r_coil_major :

    r_coil_minor :

    awp_rad :

    """
    # [m^2] Total area of all coil legs (midplane)
    tfcoil_variables.a_tf_inboard_total = (
        tfcoil_variables.n_tf_coils * tfcoil_variables.a_tf_leg_outboard
    )
    # [A] Total current in ALL coils
    tfcoil_variables.c_tf_total = tfcoil_variables.n_tf_coils * coilcurrent * 1.0e6
    # [A / m^2] overall current density
    tfcoil_variables.oacdcp = (
        tfcoil_variables.c_tf_total / tfcoil_variables.a_tf_inboard_total
    )
    # [m] radius of peak field occurrence, average
    tfcoil_variables.r_b_tf_inboard_peak_symmetric = (
        r_coil_major - r_coil_minor + awp_rad
    )
    # jlion: not sure what this will be used for. Not very
    # useful for stellarators


def calculate_inductnace(r_coil_minor):
    """This uses the reference value for the inductance and scales it with a^2/R (toroid inductance scaling)

    Parameters
    ----------
    r_coil_minor :

    """
    return (
        stellarator_configuration.stella_config_inductance
        / stellarator_variables.f_st_rmajor
        * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
        * stellarator_variables.f_st_n_coils**2
    )


def calculate_stored_magnetic_energy(r_coil_minor):
    """[GJ] Total magnetic energy

    Parameters
    ----------
    r_coil_minor :

    """
    tfcoil_variables.e_tf_magnetic_stored_total_gj = (
        0.5e0
        * (
            stellarator_configuration.stella_config_inductance
            / stellarator_variables.f_st_rmajor
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
            * stellarator_variables.f_st_n_coils**2
        )
        * (tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils) ** 2
        * 1.0e-9
    )


def calculate_winding_pack_geometry():
    """Winding Pack Geometry: for one conductor
    This one conductor will just be multiplied later to fit the winding pack size.

    """
    # [m] Dimension of square cable space inside insulation
    #     and case of the conduit of each turn
    dx_tf_turn_cable_space_average = tfcoil_variables.dx_tf_turn_general - 2.0e0 * (
        tfcoil_variables.dx_tf_turn_steel + tfcoil_variables.dx_tf_turn_insulation
    )  # dx_tf_turn_cable_space_average = t_w
    if dx_tf_turn_cable_space_average < 0:
        logger.warning(
            "Warning: Negative cable space dimension in TF coil winding pack. Check input parameters."
        )
        logger.info(
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


def calculate_current():
    """Recalculate the coil current from global stellarator configuration and variables:
    coilcurrent = f_b * stella_config_i0 * f_r / f_n
    Update stellarator_variables.f_i

    """
    coilcurrent = (
        stellarator_variables.f_st_b
        * stellarator_configuration.stella_config_i0
        * stellarator_variables.f_st_rmajor
        / stellarator_variables.f_st_n_coils
    )
    stellarator_variables.f_st_i_total = (
        coilcurrent / stellarator_configuration.stella_config_i0
    )
    return coilcurrent


def winding_pack_total_size(
    r_coil_major: float, r_coil_minor: float, coilcurrent: float
):
    # Winding Pack total size:

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
            tfcoil_variables.f_a_tf_turn_cable_copper,
            tfcoil_variables.fhts,
            tfcoil_variables.t_crit_nbti,
            tfcoil_variables.tcritsc,
            tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
            tfcoil_variables.j_tf_wp,
        )  # Get here a temperature margin from tfcoil_variables.tmargtf.

    # The operation current density weighted with the global iop/icrit fraction
    lhs[:] = constraint_variables.fiooic * jcrit_vector

    # Superconductor fraction in wp
    fraction_area_superconductor_of_wp = (
        (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        )
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_copper)
        / (tfcoil_variables.dx_tf_turn_general**2)
    )
    # print *, "f_a_scu_of_wp. ",f_a_scu_of_wp,"Awp min: ",Awp(1)

    rhs[:] = coilcurrent / (
        wp_width_r**2
        / stellarator_configuration.stella_config_wp_ratio
        * fraction_area_superconductor_of_wp
    )  # f_a_scu_of_wp should be the fraction of the sc that is in the winding pack.

    wp_width_r_min = (
        r_coil_minor / 10.0e0
    ) ** 2  # Initial guess for intersection routine
    if tfcoil_variables.i_tf_sc_mat == 6:
        wp_width_r_min = (
            r_coil_minor / 20.0e0
        ) ** 2  # If REBCO, : start at smaller winding pack ratios

    # Find the intersection between LHS and RHS (or: how much awp do I need to get to the desired coil current)
    wp_width_r_min = intersect(wp_width_r, lhs, wp_width_r, rhs, wp_width_r_min)

    # Maximum field at superconductor surface (T)
    wp_width_r_min = max(tfcoil_variables.dx_tf_turn_general**2, wp_width_r_min)

    # Recalculate tfcoil_variables.b_tf_inboard_peak_symmetric at the found awp_min:
    tfcoil_variables.b_tf_inboard_peak_symmetric = bmax_from_awp(
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

    a_tf_wp_no_insulation = awp_tor * awp_rad  # [m^2] winding-pack cross sectional area
    tfcoil_variables.j_tf_wp = (
        coilcurrent * 1.0e6 / a_tf_wp_no_insulation
    )  # [A/m^2] winding pack current density
    tfcoil_variables.n_tf_coil_turns = a_tf_wp_no_insulation / (
        tfcoil_variables.dx_tf_turn_general**2
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
    tfcoil_variables.a_tf_coil_wp_turn_insulation = tfcoil_variables.n_tf_coil_turns * (
        tfcoil_variables.dx_tf_turn_general**2
        - tfcoil_variables.a_tf_turn_steel
        - tfcoil_variables.a_tf_turn_cable_space_no_void
    )
    # [m^2] Structure area for cable
    tfcoil_variables.a_tf_wp_steel = (
        tfcoil_variables.n_tf_coil_turns * tfcoil_variables.a_tf_turn_steel
    )

    return (
        awp_rad,
        a_tf_wp_no_insulation,
        a_tf_wp_with_insulation,
        fraction_area_superconductor_of_wp,
    )


def calculate_casing():
    """Coil case thickness (m).

    Here assumed to be constant until something better comes up.
    For now assumed to be constant in a bolted plate model.

    case_thickness_constant = tfcoil_variables.dr_tf_nose_case

    #0.2e0 ?
    Leave this constant for now... Check this

    ## Should be scaled with forces I think.
    """
    # [m] coil case thickness outboard distance (radial)
    tfcoil_variables.dr_tf_plasma_case = tfcoil_variables.dr_tf_nose_case
    # dr_tf_nose_case = case_thickness_constant/2.0e0 # [m] coil case thickness inboard distance  (radial).

    # [m] coil case thickness toroidal distance (toroidal)
    tfcoil_variables.dx_tf_side_case_min = tfcoil_variables.dr_tf_nose_case


def calculate_vertical_ports():
    #  Maximal toroidal port size (vertical ports) (m)
    #  The maximal distance is correct but the vertical extension of this port is not clear#
    #  This is simplified for now and can be made more accurate in the future#
    stellarator_variables.vporttmax = (
        0.4e0
        * stellarator_configuration.stella_config_max_portsize_width
        * stellarator_variables.f_st_rmajor
        / stellarator_variables.f_st_n_coils
    )  # This is not accurate yet. Needs more insight#

    #  Maximal poloidal port size (vertical ports) (m)
    stellarator_variables.vportpmax = (
        2.0 * stellarator_variables.vporttmax
    )  # Simple approximation

    #  Maximal vertical port clearance area (m2)
    stellarator_variables.vportamax = (
        stellarator_variables.vporttmax * stellarator_variables.vportpmax
    )


def calculate_horizontal_ports():
    #  Maximal toroidal port size (horizontal ports) (m)
    stellarator_variables.hporttmax = (
        0.8e0
        * stellarator_configuration.stella_config_max_portsize_width
        * stellarator_variables.f_st_rmajor
        / stellarator_variables.f_st_n_coils
    )  # Factor 0.8 to take the variation with height into account

    #  Maximal poloidal port size (horizontal ports) (m)
    stellarator_variables.hportpmax = (
        2.0e0 * stellarator_variables.hporttmax
    )  # Simple approximation

    #  Maximal horizontal port clearance area (m2)
    stellarator_variables.hportamax = (
        stellarator_variables.hporttmax * stellarator_variables.hportpmax
    )
