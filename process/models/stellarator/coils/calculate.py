import logging

import numpy as np

from process.core.model import DataStructure
from process.models.stellarator.coils import forces
from process.models.stellarator.coils.coils import (
    bmax_from_awp,
    intersect,
    jcrit_from_material,
)
from process.models.stellarator.coils.mass import calculate_coils_mass
from process.models.stellarator.coils.output import write
from process.models.stellarator.coils.quench import calculate_quench_protection

logger = logging.getLogger(__name__)


def st_coil(stellarator, output: bool, data: DataStructure):
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

    data: DataStructure
        data structure object to provide model data

    """
    r_coil_major = data.stellarator.r_coil_major
    r_coil_minor = data.stellarator.r_coil_minor

    calculate_winding_pack_geometry(data)

    # Total coil current (MA)
    coilcurrent = calculate_current(data)

    awp_rad, a_tf_wp_no_insulation, a_tf_wp_with_insulation, f_a_scu_of_wp = (
        winding_pack_total_size(r_coil_major, r_coil_minor, coilcurrent, data)
    )

    # Casing calculations
    calculate_casing(data)

    # Port calculations
    calculate_vertical_ports(data)
    calculate_horizontal_ports(data)

    # General Coil Geometry values
    #
    calculate_coil_toroidal_thickness(data)
    calculate_coil_radial_thickness(data)

    calculate_coil_cross_sectional_area(a_tf_wp_with_insulation, data)

    calculate_coil_half_widths(data)

    calculate_plasma_facing_coil_area(data)

    coil_coil_gap, _ = calculate_coil_coil_toroidal_gap(r_coil_major, r_coil_minor, data)

    calculate_coils_summary_variables(
        coilcurrent, r_coil_major, r_coil_minor, awp_rad, data
    )

    inductance = calculate_inductance(r_coil_minor, data)
    calculate_stored_magnetic_energy(r_coil_minor, data)

    # Coil dimensions
    data.build.z_tf_inside_half = (
        0.5e0
        * data.stellarator_config.stella_config_maximal_coil_height
        * (r_coil_minor / data.stellarator_config.stella_config_coil_rminor)
    )  # [m] maximum half-height of coil

    # [m] estimated average length of a coil
    data.tfcoil.len_tf_coil = (
        data.stellarator_config.stella_config_coillength
        * (r_coil_minor / data.stellarator_config.stella_config_coil_rminor)
        / data.tfcoil.n_tf_coils
    )

    # [m^2] Total surface area of toroidal shells covering coils
    data.tfcoil.tfcryoarea = (
        data.stellarator_config.stella_config_coilsurface
        * data.stellarator.f_st_rmajor
        * (
            data.stellarator.r_coil_minor
            / data.stellarator_config.stella_config_coil_rminor
        )
        * 1.1e0
    )
    # 1.1 to scale it out a bit, as the shell must be bigger than WP

    # Minimal bending radius:
    min_bending_radius = (
        data.stellarator_config.stella_config_min_bend_radius
        * data.stellarator.f_st_rmajor
        / (1.0 - data.tfcoil.dr_tf_wp_with_insulation / (2.0 * r_coil_minor))
    )

    # End of general coil geometry values

    # Coil_mases calculations
    calculate_coils_mass(a_tf_wp_with_insulation, a_tf_wp_no_insulation, data)

    # Quench protection:
    f_vv_actual = calculate_quench_protection(coilcurrent, data)

    #
    # Forces scaling #
    forces.calculate_max_force_density(a_tf_wp_no_insulation, data)
    forces.calculate_maximum_stress(data)

    # Units: MN/m
    max_force_density_mnm = forces.calculate_max_force_density_mnm(data)
    max_lateral_force_density = forces.calculate_max_lateral_force_density(
        a_tf_wp_no_insulation, data
    )
    max_radial_force_density = forces.calculate_max_radial_force_density(
        a_tf_wp_no_insulation, data
    )
    #
    # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
    centering_force_max_mn = forces.calculate_centering_force_max_mn(data)
    centering_force_min_mn = forces.calculate_centering_force_min_mn(data)
    centering_force_avg_mn = forces.calculate_centering_force_avg_mn(data)

    if output:
        write(
            stellarator=stellarator,
            a_tf_wp_no_insulation=a_tf_wp_no_insulation,
            centering_force_avg_mn=centering_force_avg_mn,
            centering_force_max_mn=centering_force_max_mn,
            centering_force_min_mn=centering_force_min_mn,
            coilcoilgap=coil_coil_gap,
            coppera_m2=data.rebco.coppera_m2,
            coppera_m2_max=data.rebco.coppera_m2_max,
            f_a_scu_of_wp=f_a_scu_of_wp,
            f_vv_actual=f_vv_actual,
            f_j_tf_wp_critical_max=data.constraints.f_j_tf_wp_critical_max,
            inductance=inductance,
            max_force_density=data.tfcoil.max_force_density,
            max_force_density_mnm=max_force_density_mnm,
            max_lateral_force_density=max_lateral_force_density,
            max_radial_force_density=max_radial_force_density,
            min_bending_radius=min_bending_radius,
            r_coil_major=r_coil_major,
            r_coil_minor=r_coil_minor,
            sig_tf_wp=data.tfcoil.sig_tf_wp,
            dx_tf_turn_general=data.tfcoil.dx_tf_turn_general,
            t_tf_superconductor_quench=data.tfcoil.t_tf_superconductor_quench,
            toroidalgap=data.tfcoil.toroidalgap,
            allowed_quench_voltage=data.tfcoil.v_tf_coil_dump_quench_max_kv,
            quench_voltage=data.tfcoil.v_tf_coil_dump_quench_kv,
            data=data,
        )


def calculate_coil_toroidal_thickness(data: DataStructure):
    data.tfcoil.dx_tf_inboard_out_toroidal = (
        data.tfcoil.dx_tf_wp_primary_toroidal
        + 2.0e0 * data.tfcoil.dx_tf_side_case_min
        + 2.0e0 * data.tfcoil.dx_tf_wp_insulation
    )  # [m] Thickness of inboard leg in toroidal direction


def calculate_coil_radial_thickness(data):
    """Thickness of inboard and outboard leg in radial direction"""
    # [m] Thickness of inboard leg in radial direction
    data.build.dr_tf_inboard = (
        data.tfcoil.dr_tf_nose_case
        + data.tfcoil.dr_tf_wp_with_insulation
        + data.tfcoil.dr_tf_plasma_case
        + 2.0e0 * data.tfcoil.dx_tf_wp_insulation
    )
    # [m] Thickness of outboard leg in radial direction (same as inboard)
    data.build.dr_tf_outboard = data.build.dr_tf_inboard


def calculate_coil_cross_sectional_area(a_tf_wp_with_insulation, data):
    # [m^2] overall coil cross-sectional area
    # (assuming inboard and outboard leg are the same)
    data.tfcoil.a_tf_leg_outboard = (
        data.build.dr_tf_inboard * data.tfcoil.dx_tf_inboard_out_toroidal
    )
    # [m^2] Cross-sectional area of surrounding case
    data.tfcoil.a_tf_coil_inboard_case = (
        data.build.dr_tf_inboard * data.tfcoil.dx_tf_inboard_out_toroidal
    ) - a_tf_wp_with_insulation


def calculate_coil_half_widths(data: DataStructure):
    # [m] Half-width of side of coil nearest torus centreline
    data.tfcoil.tfocrn = 0.5e0 * data.tfcoil.dx_tf_inboard_out_toroidal
    # [m] Half-width of side of coil nearest plasma
    data.tfcoil.tficrn = 0.5e0 * data.tfcoil.dx_tf_inboard_out_toroidal


def calculate_plasma_facing_coil_area(data: DataStructure):
    # [m^2] Total surface area of coil side facing plasma: inboard region
    data.tfcoil.tfsai = (
        data.tfcoil.n_tf_coils
        * data.tfcoil.dx_tf_inboard_out_toroidal
        * 0.5e0
        * data.tfcoil.len_tf_coil
    )
    # [m^2] Total surface area of coil side facing plasma: outboard region
    # (same as inboard)
    data.tfcoil.tfsao = data.tfcoil.tfsai


def calculate_coil_coil_toroidal_gap(r_coil_major, r_coil_minor, data: DataStructure):
    """[m] Minimal distance in toroidal direction between two stellarator coils
    Consistency with coil width is checked in constraint equation 82

    Parameters
    ----------
    r_coil_major :

    r_coil_minor :

    data: DataStructure
        data structure object

    """
    # [m] Toroidal gap between two coil filaments
    data.tfcoil.toroidalgap = (
        data.stellarator_config.stella_config_dmin
        * (r_coil_major - r_coil_minor)
        / (
            data.stellarator_config.stella_config_coil_rmajor
            - data.stellarator_config.stella_config_coil_rminor
        )
    )
    # Left-Over coil gap between two coils (m)
    coilcoilgap = data.tfcoil.toroidalgap - data.tfcoil.dx_tf_inboard_out_toroidal
    return coilcoilgap, data.tfcoil.toroidalgap


def calculate_coils_summary_variables(
    coilcurrent, r_coil_major, r_coil_minor, awp_rad, data: DataStructure
):
    """Variables for ALL coils.

    Parameters
    ----------
    coilcurrent :

    r_coil_major :

    r_coil_minor :

    awp_rad :

    """
    # [m^2] Total area of all coil legs (midplane)
    data.tfcoil.a_tf_inboard_total = (
        data.tfcoil.n_tf_coils * data.tfcoil.a_tf_leg_outboard
    )
    # [A] Total current in ALL coils
    data.tfcoil.c_tf_total = data.tfcoil.n_tf_coils * coilcurrent * 1.0e6
    # [A / m^2] overall current density
    data.tfcoil.j_tf_coil_full_area = (
        data.tfcoil.c_tf_total / data.tfcoil.a_tf_inboard_total
    )
    # [m] radius of peak field occurrence, average
    data.tfcoil.r_b_tf_inboard_peak_symmetric = r_coil_major - r_coil_minor + awp_rad
    # jlion: not sure what this will be used for. Not very
    # useful for stellarators


def calculate_inductance(r_coil_minor, data: DataStructure):
    """This uses the reference value for the inductance and scales it with a^2/R
    (toroid inductance scaling)

    Parameters
    ----------
    r_coil_minor :

    data: DataStructure
        data structure object
    """
    return (
        data.stellarator_config.stella_config_inductance
        / data.stellarator.f_st_rmajor
        * (r_coil_minor / data.stellarator_config.stella_config_coil_rminor) ** 2
        * data.stellarator.f_st_n_coils**2
    )


def calculate_stored_magnetic_energy(r_coil_minor, data: DataStructure):
    """[GJ] Total magnetic energy

    Parameters
    ----------
    r_coil_minor :

    data: DataStructure
        data structure object
    """
    data.tfcoil.e_tf_magnetic_stored_total_gj = (
        0.5e0
        * (
            data.stellarator_config.stella_config_inductance
            / data.stellarator.f_st_rmajor
            * (r_coil_minor / data.stellarator_config.stella_config_coil_rminor) ** 2
            * data.stellarator.f_st_n_coils**2
        )
        * (data.tfcoil.c_tf_total / data.tfcoil.n_tf_coils) ** 2
        * 1.0e-9
    )


def calculate_winding_pack_geometry(data: DataStructure):
    """Winding Pack Geometry: for one conductor
    This one conductor will just be multiplied later to fit the winding pack size.

    """
    # [m] Dimension of square cable space inside insulation
    #    and case of the conduit of each turn
    dx_tf_turn_cable_space_average = data.tfcoil.dx_tf_turn_general - 2.0e0 * (
        data.tfcoil.dx_tf_turn_steel + data.tfcoil.dx_tf_turn_insulation
    )  # dx_tf_turn_cable_space_average = t_w
    if dx_tf_turn_cable_space_average < 0:
        logger.warning(
            "Warning: Negative cable space dimension in TF coil winding pack. "
            "Check input parameters."
        )
        logger.info(
            "dx_tf_turn_cable_space_average is negative. "
            "Check t_turn, data.tfcoil.dx_tf_turn_steel and dx_tf_turn_insulation."
        )
    # [m^2] Cross-sectional area of cable space per turn
    # 0.9 to include some rounded corners.
    # (data.tfcoil.a_tf_turn_cable_space_no_void =
    # pi (dx_tf_turn_cable_space_average/2)**2 =
    # pi/4 *dx_tf_turn_cable_space_average**2 for perfect round conductor).
    # This factor depends on how round the corners are.
    data.tfcoil.a_tf_turn_cable_space_no_void = 0.9e0 * dx_tf_turn_cable_space_average**2
    # [m^2] Cross-sectional area of conduit case per turn
    data.tfcoil.a_tf_turn_steel = (
        dx_tf_turn_cable_space_average + 2.0e0 * data.tfcoil.dx_tf_turn_steel
    ) ** 2 - data.tfcoil.a_tf_turn_cable_space_no_void


def calculate_current(data: DataStructure):
    """Recalculate the coil current from global stellarator configuration and variables:
    coilcurrent = f_b * stella_config_i0 * f_r / f_n
    Update data.stellarator.f_i

    """
    coilcurrent = (
        data.stellarator.f_st_b
        * data.stellarator_config.stella_config_i0
        * data.stellarator.f_st_rmajor
        / data.stellarator.f_st_n_coils
    )
    data.stellarator.f_st_i_total = (
        coilcurrent / data.stellarator_config.stella_config_i0
    )
    return coilcurrent


def winding_pack_total_size(
    r_coil_major: float, r_coil_minor: float, coilcurrent: float, data: DataStructure
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
        if data.tfcoil.i_tf_sc_mat == 6:
            wp_width_r[k] = (r_coil_minor / 150.0e0) + (k / (n_it - 1e0)) * (
                r_coil_minor / 1.0e0 - r_coil_minor / 150.0e0
            )

        # B-field calculation
        b_max_k[k] = bmax_from_awp(
            wp_width_r[k],
            coilcurrent,
            data.tfcoil.n_tf_coils,
            r_coil_major,
            r_coil_minor,
            data,
        )
        # Two margins can be applied for jcrit: direct or by temperature margin.
        # Temperature margin is implemented in the jcrit_vector definition,
        # direct margin is implemented after jcrit is defined (equation below)
        # jcrit for this bmax:
        jcrit_vector[k] = jcrit_from_material(
            b_max_k[k],
            data.tfcoil.tftmp + data.tfcoil.tmargmin,
            data.tfcoil.i_tf_sc_mat,
            data.tfcoil.b_crit_upper_nbti,
            data.tfcoil.bcritsc,
            data.tfcoil.f_a_tf_turn_cable_copper,
            data.tfcoil.fhts,
            data.tfcoil.t_crit_nbti,
            data.tfcoil.tcritsc,
            data.tfcoil.f_a_tf_turn_cable_space_extra_void,
            data.tfcoil.j_tf_wp,
        )  # Get here a temperature margin from data.tfcoil.tmargtf.

    # The operation current density weighted with the global iop/icrit fraction
    lhs[:] = data.constraints.f_j_tf_wp_critical_max * jcrit_vector

    # Superconductor fraction in wp
    fraction_area_superconductor_of_wp = (
        (
            data.tfcoil.a_tf_turn_cable_space_no_void
            * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_space_extra_void)
        )
        * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_copper)
        / (data.tfcoil.dx_tf_turn_general**2)
    )
    # print *, "f_a_scu_of_wp. ",f_a_scu_of_wp,"Awp min: ",Awp(1)

    rhs[:] = coilcurrent / (
        wp_width_r**2
        / data.stellarator_config.stella_config_wp_ratio
        * fraction_area_superconductor_of_wp
    )  # f_a_scu_of_wp should be the fraction of the sc that is in the winding pack.

    wp_width_r_min = (
        r_coil_minor / 10.0e0
    ) ** 2  # Initial guess for intersection routine
    if data.tfcoil.i_tf_sc_mat == 6:
        wp_width_r_min = (
            r_coil_minor / 20.0e0
        ) ** 2  # If REBCO, : start at smaller winding pack ratios

    # Find the intersection between LHS and RHS
    # (or: how much awp do I need to get to the desired coil current)
    wp_width_r_min = intersect(wp_width_r, lhs, wp_width_r, rhs, wp_width_r_min)

    # Maximum field at superconductor surface (T)
    wp_width_r_min = max(data.tfcoil.dx_tf_turn_general**2, wp_width_r_min)

    # Recalculate data.tfcoil.b_tf_inboard_peak_symmetric at the found awp_min:
    data.tfcoil.b_tf_inboard_peak_symmetric = bmax_from_awp(
        wp_width_r_min,
        coilcurrent,
        data.tfcoil.n_tf_coils,
        r_coil_major,
        r_coil_minor,
        data,
    )

    # Winding pack toroidal, radial cross-sections (m)
    awp_tor = (
        wp_width_r_min / data.stellarator_config.stella_config_wp_ratio
    )  # Toroidal dimension
    awp_rad = wp_width_r_min  # Radial dimension

    data.tfcoil.dx_tf_wp_primary_toroidal = (
        awp_tor  # [m] toroidal thickness of winding pack
    )
    data.tfcoil.dx_tf_wp_secondary_toroidal = (
        awp_tor  # [m] toroidal thickness of winding pack (region in front)
    )
    data.tfcoil.dr_tf_wp_with_insulation = (
        awp_rad  # [m] radial thickness of winding pack
    )

    # [m^2] winding-pack cross sectional area including insulation (not global)
    a_tf_wp_with_insulation = (
        data.tfcoil.dr_tf_wp_with_insulation + 2.0e0 * data.tfcoil.dx_tf_wp_insulation
    ) * (data.tfcoil.dx_tf_wp_primary_toroidal + 2.0e0 * data.tfcoil.dx_tf_wp_insulation)

    a_tf_wp_no_insulation = awp_tor * awp_rad  # [m^2] winding-pack cross sectional area
    data.tfcoil.j_tf_wp = (
        coilcurrent * 1.0e6 / a_tf_wp_no_insulation
    )  # [A/m^2] winding pack current density
    data.tfcoil.n_tf_coil_turns = a_tf_wp_no_insulation / (
        data.tfcoil.dx_tf_turn_general**2
    )  # estimated number of turns for a given turn size (not global). Take at least 1.
    data.tfcoil.c_tf_turn = (
        coilcurrent * 1.0e6 / data.tfcoil.n_tf_coil_turns
    )  # [A] current per turn - estimation
    # [m^2] Total conductor cross-sectional area, taking account of void area
    data.tfcoil.a_tf_wp_conductor = (
        data.tfcoil.a_tf_turn_cable_space_no_void
        * data.tfcoil.n_tf_coil_turns
        * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_space_extra_void)
    )
    # [m^2] Void area in cable, for He
    data.tfcoil.a_tf_wp_extra_void = (
        data.tfcoil.a_tf_turn_cable_space_no_void
        * data.tfcoil.n_tf_coil_turns
        * data.tfcoil.f_a_tf_turn_cable_space_extra_void
    )
    # [m^2] Insulation area (not including ground-wall)
    data.tfcoil.a_tf_coil_wp_turn_insulation = data.tfcoil.n_tf_coil_turns * (
        data.tfcoil.dx_tf_turn_general**2
        - data.tfcoil.a_tf_turn_steel
        - data.tfcoil.a_tf_turn_cable_space_no_void
    )
    # [m^2] Structure area for cable
    data.tfcoil.a_tf_wp_steel = data.tfcoil.n_tf_coil_turns * data.tfcoil.a_tf_turn_steel

    return (
        awp_rad,
        a_tf_wp_no_insulation,
        a_tf_wp_with_insulation,
        fraction_area_superconductor_of_wp,
    )


def calculate_casing(data: DataStructure):
    """Coil case thickness (m).

    Here assumed to be constant until something better comes up.
    For now assumed to be constant in a bolted plate model.

    case_thickness_constant = data.tfcoil.dr_tf_nose_case

    #0.2e0 ?
    Leave this constant for now... Check this

    ## Should be scaled with forces I think.
    """
    # [m] coil case thickness outboard distance (radial)
    data.tfcoil.dr_tf_plasma_case = data.tfcoil.dr_tf_nose_case

    # [m] coil case thickness toroidal distance (toroidal)
    data.tfcoil.dx_tf_side_case_min = data.tfcoil.dr_tf_nose_case


def calculate_vertical_ports(data: DataStructure):
    # Maximal toroidal port size (vertical ports) (m)
    # The maximal distance is correct
    # but the vertical extension of this port is not clear
    # This is simplified for now and can be made more accurate in the future#
    data.stellarator.vporttmax = (
        0.4e0
        * data.stellarator_config.stella_config_max_portsize_width
        * data.stellarator.f_st_rmajor
        / data.stellarator.f_st_n_coils
    )  # This is not accurate yet. Needs more insight#

    # Maximal poloidal port size (vertical ports) (m)
    data.stellarator.vportpmax = 2.0 * data.stellarator.vporttmax  # Simple approximation

    # Maximal vertical port clearance area (m2)
    data.stellarator.vportamax = data.stellarator.vporttmax * data.stellarator.vportpmax


def calculate_horizontal_ports(data: DataStructure):
    # Maximal toroidal port size (horizontal ports) (m)
    data.stellarator.hporttmax = (
        0.8e0
        * data.stellarator_config.stella_config_max_portsize_width
        * data.stellarator.f_st_rmajor
        / data.stellarator.f_st_n_coils
    )  # Factor 0.8 to take the variation with height into account

    # Maximal poloidal port size (horizontal ports) (m)
    data.stellarator.hportpmax = (
        2.0e0 * data.stellarator.hporttmax
    )  # Simple approximation

    # Maximal horizontal port clearance area (m2)
    data.stellarator.hportamax = data.stellarator.hporttmax * data.stellarator.hportpmax
