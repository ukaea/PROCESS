from process import process_output as po

from process.fortran import (
    build_variables,
    stellarator_variables,
    tfcoil_variables,
)

def write(
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
        sig_tf_wp,
        dx_tf_turn_general,
        t_tf_superconductor_quench,
        toroidalgap,
        allowed_quench_voltage,
        quench_voltage,
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
            "(j_tf_wp_quench_heat_max)",
            tfcoil_variables.j_tf_wp_quench_heat_max,
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
            "(b_tf_inboard_peak_symmetric)",
            tfcoil_variables.b_tf_inboard_peak_symmetric,
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

        po.osubhd(stellarator.outfile, "Conductor Information :")
        po.ovarre(
            stellarator.outfile,
            "Superconductor mass per coil (kg)",
            "(m_tf_coil_superconductor)",
            tfcoil_variables.m_tf_coil_superconductor,
        )
        po.ovarre(
            stellarator.outfile,
            "Copper mass per coil (kg)",
            "(m_tf_coil_copper)",
            tfcoil_variables.m_tf_coil_copper,
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
            "(m_tf_coil_conductor)",
            tfcoil_variables.m_tf_coil_conductor,
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
            "(f_a_tf_turn_cable_copper)",
            tfcoil_variables.f_a_tf_turn_cable_copper,
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
            "(dx_tf_turn_general)",
            dx_tf_turn_general,
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
            "(t_tf_superconductor_quench)",
            t_tf_superconductor_quench,
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
            "(v_tf_coil_dump_quench_max_kv)",
            allowed_quench_voltage,
        )
        po.ovarre(
                stellarator.outfile,
                "Actual quench voltage (kV)", 
                "(v_tf_coil_dump_quench_kv)", 
                quench_voltage, 
                "OP "
        )
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
            "(m_tf_coil_case)",
            tfcoil_variables.m_tf_coil_case,
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