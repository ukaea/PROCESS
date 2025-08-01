import copy
import logging

import numpy as np
from scipy import optimize

import process.superconductors as superconductors
from process import process_output as po
from process.data_structure import divertor_variables, rebco_variables
from process.exceptions import ProcessValueError
from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    error_handling,
    global_variables,
    pfcoil_variables,
    physics_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.quench import calculate_quench_protection_current_density
from process.tf_coil import TFCoil
from process.utilities.f2py_string_patch import f2py_compatible_to_string

logger = logging.getLogger(__name__)

SUPERCONDUCTING_TF_TYPES = {
    1: "Nb3Sn ITER",
    2: "Bi-2212",
    3: "NbTi",
    4: "Nb3Sn user",
    5: "Nb3Sn WST",
    6: "REBCO Croco",
    7: "NbTi Ginzburg-Landau",
    8: "REBCO Ginzburg-Landau",
    9: "REBCO Hazelton-Zhai",
}


RMU0 = constants.rmu0
EPS = np.finfo(1.0).eps


class SuperconductingTFCoil(TFCoil):
    def __init__(self):
        self.outfile = constants.nout

    def run(self, output: bool):
        """
        Routine to call the superconductor module for the TF coils
        """
        self.iprint = 0
        (
            sctfcoil_module.rad_tf_coil_inboard_toroidal_half,
            sctfcoil_module.tan_theta_coil,
            tfcoil_variables.a_tf_inboard_total,
            sctfcoil_module.r_tf_outboard_in,
            sctfcoil_module.r_tf_outboard_out,
            tfcoil_variables.dx_tf_inboard_out_toroidal,
            tfcoil_variables.a_tf_leg_outboard,
            tfcoil_variables.dr_tf_plasma_case,
            tfcoil_variables.dx_tf_side_case_min,
        ) = super().tf_global_geometry(
            i_tf_case_geom=tfcoil_variables.i_tf_case_geom,
            i_f_dr_tf_plasma_case=tfcoil_variables.i_f_dr_tf_plasma_case,
            f_dr_tf_plasma_case=tfcoil_variables.f_dr_tf_plasma_case,
            tfc_sidewall_is_fraction=tfcoil_variables.tfc_sidewall_is_fraction,
            casths_fraction=tfcoil_variables.casths_fraction,
            n_tf_coils=tfcoil_variables.n_tf_coils,
            dr_tf_inboard=build_variables.dr_tf_inboard,
            dr_tf_nose_case=tfcoil_variables.dr_tf_nose_case,
            r_tf_inboard_out=build_variables.r_tf_inboard_out,
            r_tf_inboard_in=build_variables.r_tf_inboard_in,
            r_tf_outboard_mid=build_variables.r_tf_outboard_mid,
            dr_tf_outboard=build_variables.dr_tf_outboard,
        )

        # Radial position of peak toroidal field [m]

        tfcoil_variables.r_b_tf_inboard_peak = (
            build_variables.r_tf_inboard_out
            - tfcoil_variables.dr_tf_plasma_case
            - tfcoil_variables.dx_tf_wp_insulation
            - tfcoil_variables.dx_tf_wp_insertion_gap
        )

        (
            tfcoil_variables.b_tf_inboard_peak,
            tfcoil_variables.c_tf_total,
            sctfcoil_module.c_tf_coil,
            tfcoil_variables.oacdcp,
        ) = super().tf_current(
            n_tf_coils=tfcoil_variables.n_tf_coils,
            bt=physics_variables.bt,
            rmajor=physics_variables.rmajor,
            r_b_tf_inboard_peak=tfcoil_variables.r_b_tf_inboard_peak,
            a_tf_inboard_total=tfcoil_variables.a_tf_inboard_total,
        )

        (
            tfcoil_variables.len_tf_coil,
            tfcoil_variables.tfa,
            tfcoil_variables.tfb,
            tfcoil_variables.r_tf_arc,
            tfcoil_variables.z_tf_arc,
        ) = super().tf_coil_shape_inner(
            i_tf_shape=tfcoil_variables.i_tf_shape,
            itart=physics_variables.itart,
            i_single_null=physics_variables.i_single_null,
            r_tf_inboard_out=build_variables.r_tf_inboard_out,
            r_cp_top=build_variables.r_cp_top,
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            r_tf_outboard_in=sctfcoil_module.r_tf_outboard_in,
            z_tf_inside_half=build_variables.z_tf_inside_half,
            z_tf_top=build_variables.z_tf_top,
            dr_tf_inboard=build_variables.dr_tf_inboard,
            dr_tf_outboard=build_variables.dr_tf_outboard,
            r_tf_outboard_mid=build_variables.r_tf_outboard_mid,
            r_tf_inboard_mid=build_variables.r_tf_inboard_mid,
        )

        self.sc_tf_internal_geom(
            tfcoil_variables.i_tf_wp_geom,
            tfcoil_variables.i_tf_case_geom,
            tfcoil_variables.i_tf_turns_integer,
        )

        tfcoil_variables.ind_tf_coil = super().tf_coil_self_inductance(
            dr_tf_inboard=build_variables.dr_tf_inboard,
            r_tf_arc=tfcoil_variables.r_tf_arc,
            z_tf_arc=tfcoil_variables.z_tf_arc,
            itart=physics_variables.itart,
            i_tf_shape=tfcoil_variables.i_tf_shape,
            z_tf_inside_half=build_variables.z_tf_inside_half,
            dr_tf_outboard=build_variables.dr_tf_outboard,
            r_tf_outboard_mid=build_variables.r_tf_outboard_mid,
            r_tf_inboard_mid=build_variables.r_tf_inboard_mid,
        )

        # Total TF coil stored magnetic energy [J]
        sctfcoil_module.e_tf_magnetic_stored_total = (
            0.5e0 * tfcoil_variables.ind_tf_coil * tfcoil_variables.c_tf_total**2
        )

        # Total TF coil stored magnetic energy [Gigajoule]
        tfcoil_variables.e_tf_magnetic_stored_total_gj = (
            1.0e-9 * sctfcoil_module.e_tf_magnetic_stored_total
        )

        self.tf_field_and_force()

        # Calculate TF coil areas and masses
        self.tf_coil_area_and_masses()

        # Do stress calculations (writes the stress output)
        if output:
            tfcoil_variables.n_rad_per_layer = 500

        try:
            (
                sig_tf_r_max,
                sig_tf_t_max,
                sig_tf_z_max,
                sig_tf_vmises_max,
                s_shear_tf_peak,
                deflect,
                eyoung_axial,
                eyoung_trans,
                eyoung_wp_axial,
                eyoung_wp_trans,
                poisson_wp_trans,
                radial_array,
                s_shear_cea_tf_cond,
                poisson_wp_axial,
                sig_tf_r,
                sig_tf_smeared_r,
                sig_tf_smeared_t,
                sig_tf_smeared_z,
                sig_tf_t,
                s_shear_tf,
                sig_tf_vmises,
                sig_tf_z,
                str_tf_r,
                str_tf_t,
                str_tf_z,
                n_radial_array,
                n_tf_bucking,
                tfcoil_variables.sig_tf_wp,
                sig_tf_case,
                sig_tf_cs_bucked,
                str_wp,
                casestr,
                insstrain,
                sig_tf_wp_av_z,
            ) = self.stresscl(
                int(tfcoil_variables.n_tf_stress_layers),
                int(tfcoil_variables.n_rad_per_layer),
                int(tfcoil_variables.n_tf_wp_layers),
                int(tfcoil_variables.i_tf_bucking),
                float(build_variables.r_tf_inboard_in),
                build_variables.dr_bore,
                build_variables.z_tf_inside_half,
                pfcoil_variables.f_z_cs_tf_internal,
                build_variables.dr_cs,
                build_variables.i_tf_inside_cs,
                build_variables.dr_tf_inboard,
                build_variables.dr_cs_tf_gap,
                pfcoil_variables.i_pf_conductor,
                pfcoil_variables.j_cs_flat_top_end,
                pfcoil_variables.j_cs_pulse_start,
                pfcoil_variables.c_pf_coil_turn_peak_input,
                pfcoil_variables.n_pf_coils_in_group,
                pfcoil_variables.ld_ratio_cst,
                pfcoil_variables.r_out_cst,
                pfcoil_variables.f_a_cs_steel,
                tfcoil_variables.eyoung_steel,
                tfcoil_variables.poisson_steel,
                tfcoil_variables.eyoung_cond_axial,
                tfcoil_variables.poisson_cond_axial,
                tfcoil_variables.eyoung_cond_trans,
                tfcoil_variables.poisson_cond_trans,
                tfcoil_variables.eyoung_ins,
                tfcoil_variables.poisson_ins,
                tfcoil_variables.dx_tf_turn_insulation,
                tfcoil_variables.eyoung_copper,
                tfcoil_variables.poisson_copper,
                tfcoil_variables.i_tf_sup,
                tfcoil_variables.eyoung_res_tf_buck,
                sctfcoil_module.r_tf_wp_inboard_inner,
                sctfcoil_module.tan_theta_coil,
                sctfcoil_module.rad_tf_coil_inboard_toroidal_half,
                sctfcoil_module.r_tf_wp_inboard_outer,
                sctfcoil_module.a_tf_coil_inboard_steel,
                sctfcoil_module.a_tf_plasma_case,
                sctfcoil_module.a_tf_coil_nose_case,
                tfcoil_variables.dx_tf_wp_insertion_gap,
                tfcoil_variables.dx_tf_wp_insulation,
                tfcoil_variables.n_tf_coil_turns,
                int(tfcoil_variables.i_tf_turns_integer),
                sctfcoil_module.dx_tf_turn_cable_space_average,
                sctfcoil_module.dr_tf_turn_cable_space,
                tfcoil_variables.dia_tf_turn_coolant_channel,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.dx_tf_turn_steel,
                sctfcoil_module.dx_tf_side_case_average,
                sctfcoil_module.dx_tf_wp_toroidal_average,
                sctfcoil_module.a_tf_coil_inboard_insulation,
                tfcoil_variables.a_tf_wp_steel,
                tfcoil_variables.a_tf_wp_conductor,
                sctfcoil_module.a_tf_wp_with_insulation,
                tfcoil_variables.eyoung_al,
                tfcoil_variables.poisson_al,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf_graded_layers,
                tfcoil_variables.c_tf_total,
                tfcoil_variables.dr_tf_plasma_case,
                tfcoil_variables.i_tf_stress_model,
                sctfcoil_module.vforce_inboard_tot,
                tfcoil_variables.i_tf_tresca,
                tfcoil_variables.a_tf_coil_inboard_case,
                tfcoil_variables.vforce,
                tfcoil_variables.a_tf_turn_steel,
            )

            tfcoil_variables.sig_tf_case = (
                tfcoil_variables.sig_tf_case
                if tfcoil_variables.sig_tf_case is None
                else sig_tf_case
            )

            tfcoil_variables.sig_tf_cs_bucked = (
                tfcoil_variables.sig_tf_cs_bucked
                if tfcoil_variables.sig_tf_cs_bucked is None
                else sig_tf_cs_bucked
            )

            tfcoil_variables.str_wp = (
                tfcoil_variables.str_wp if tfcoil_variables.str_wp is None else str_wp
            )

            tfcoil_variables.casestr = (
                tfcoil_variables.casestr
                if tfcoil_variables.casestr is None
                else casestr
            )

            tfcoil_variables.insstrain = (
                tfcoil_variables.insstrain
                if tfcoil_variables.insstrain is None
                else insstrain
            )

            if output:
                self.out_stress(
                    sig_tf_r_max,
                    sig_tf_t_max,
                    sig_tf_z_max,
                    sig_tf_vmises_max,
                    s_shear_tf_peak,
                    deflect,
                    eyoung_axial,
                    eyoung_trans,
                    eyoung_wp_axial,
                    eyoung_wp_trans,
                    poisson_wp_trans,
                    radial_array,
                    s_shear_cea_tf_cond,
                    poisson_wp_axial,
                    sig_tf_r,
                    sig_tf_smeared_r,
                    sig_tf_smeared_t,
                    sig_tf_smeared_z,
                    sig_tf_t,
                    s_shear_tf,
                    sig_tf_vmises,
                    sig_tf_z,
                    str_tf_r,
                    str_tf_t,
                    str_tf_z,
                    n_radial_array,
                    n_tf_bucking,
                    sig_tf_wp_av_z,
                )
        except ValueError as e:
            if e.args[1] == 245 and e.args[2] == 0:
                error_handling.report_error(245)
                tfcoil_variables.sig_tf_case = 0.0e0
                tfcoil_variables.sig_tf_wp = 0.0e0
        peaktfflag = 0
        self.vv_stress_on_quench()

        # Peak field including ripple
        # Rem : as resistive magnets are axisymmetric, no inboard ripple is present
        tfcoil_variables.bmaxtfrp, peaktfflag = self.peak_tf_with_ripple(
            tfcoil_variables.n_tf_coils,
            tfcoil_variables.dx_tf_wp_primary_toroidal,
            tfcoil_variables.dr_tf_wp_with_insulation
            - 2.0e0
            * (
                tfcoil_variables.dx_tf_wp_insulation
                + tfcoil_variables.dx_tf_wp_insertion_gap
            ),
            sctfcoil_module.r_tf_wp_inboard_centre,
            tfcoil_variables.b_tf_inboard_peak,
        )

        tfes = sctfcoil_module.e_tf_magnetic_stored_total / tfcoil_variables.n_tf_coils
        # Cross-sectional area per turn
        a_tf_turn = tfcoil_variables.c_tf_total / (
            tfcoil_variables.j_tf_wp
            * tfcoil_variables.n_tf_coils
            * tfcoil_variables.n_tf_coil_turns
        )

        if tfcoil_variables.i_tf_sc_mat == 6:
            (tfcoil_variables.j_tf_wp_critical, tfcoil_variables.tmargtf) = (
                self.supercon_croco(
                    a_tf_turn,
                    tfcoil_variables.bmaxtfrp,
                    tfcoil_variables.c_tf_turn,
                    tfcoil_variables.tftmp,
                    output=output,
                )
            )

            tfcoil_variables.vtfskv = (
                self.croco_voltage() / 1.0e3
            )  # TFC Quench voltage in kV

        else:
            (
                tfcoil_variables.j_tf_wp_critical,
                vdump,
                tfcoil_variables.tmargtf,
            ) = self.supercon(
                tfcoil_variables.a_tf_turn_cable_space_no_void,
                a_tf_turn,
                tfcoil_variables.bmaxtfrp,
                tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.c_tf_turn,
                tfcoil_variables.j_tf_wp,
                tfcoil_variables.i_tf_sc_mat,
                tfcoil_variables.fhts,
                tfcoil_variables.tdmptf,
                tfes,
                tfcoil_variables.tftmp,
                tfcoil_variables.tmaxpro,
                tfcoil_variables.bcritsc,
                tfcoil_variables.tcritsc,
                output=output,
            )

            tfcoil_variables.vtfskv = vdump / 1.0e3  # TFC Quench voltage in kV

            if output:
                self.outtf(peaktfflag)

    def croco_voltage(self) -> float:
        if f2py_compatible_to_string(tfcoil_variables.quench_model) == "linear":
            sctfcoil_module.time2 = tfcoil_variables.tdmptf
            croco_voltage = (
                2.0e0
                / sctfcoil_module.time2
                * (
                    sctfcoil_module.e_tf_magnetic_stored_total
                    / tfcoil_variables.n_tf_coils
                )
                / tfcoil_variables.c_tf_turn
            )
        elif f2py_compatible_to_string(tfcoil_variables.quench_model) == "exponential":
            sctfcoil_module.tau2 = tfcoil_variables.tdmptf
            croco_voltage = (
                2.0e0
                / sctfcoil_module.tau2
                * (
                    sctfcoil_module.e_tf_magnetic_stored_total
                    / tfcoil_variables.n_tf_coils
                )
                / tfcoil_variables.c_tf_turn
            )
        else:
            return 0.0

        return croco_voltage

    def supercon_croco(self, a_tf_turn, b_tf_inboard_peak, iop, thelium, output: bool):
        """TF superconducting CroCo conductor using REBCO tape
        author: M Kovari, CCFE, Culham Science Centre
        b_tf_inboard_peak : input real : Peak field at conductor (T)
        iop : input real : Operating current per turn (A)
        thelium : input real : He temperature at peak field point (K)
        iprint : input integer : Switch for printing (1 = yes, 0 = no)
        outfile : input integer : Fortran output unit identifier
        j_tf_wp_critical : output real : Critical winding pack current density (A/m2)
        tmarg : output real : Temperature margin (K)
        """

        j_crit_sc: float = 0.0
        #  Find critical current density in superconducting cable, j_crit_cable
        j_crit_sc, _ = superconductors.jcrit_rebco(thelium, b_tf_inboard_peak)
        # tfcoil_variables.a_tf_turn_cable_space_no_void : Cable space - inside area (m2)
        # Set new rebco_variables.croco_od
        # allowing for scaling of rebco_variables.croco_od
        rebco_variables.croco_od = (
            tfcoil_variables.t_conductor / 3.0e0
            - tfcoil_variables.dx_tf_turn_steel * (2.0e0 / 3.0e0)
        )
        sctfcoil_module.conductor_acs = (
            9.0e0 / 4.0e0 * np.pi * rebco_variables.croco_od**2
        )
        tfcoil_variables.a_tf_turn_cable_space_no_void = sctfcoil_module.conductor_acs
        sctfcoil_module.conductor_area = (
            tfcoil_variables.t_conductor**2
        )  # does this not assume it's a sqaure???

        sctfcoil_module.conductor_jacket_area = (
            sctfcoil_module.conductor_area - sctfcoil_module.conductor_acs
        )
        tfcoil_variables.a_tf_turn_steel = sctfcoil_module.conductor_jacket_area

        sctfcoil_module.conductor_jacket_fraction = (
            sctfcoil_module.conductor_jacket_area / sctfcoil_module.conductor_area
        )
        (
            sctfcoil_module.croco_strand_area,
            sctfcoil_module.croco_strand_critical_current,
            sctfcoil_module.conductor_copper_area,
            sctfcoil_module.conductor_copper_fraction,
            sctfcoil_module.conductor_copper_bar_area,
            sctfcoil_module.conductor_hastelloy_area,
            sctfcoil_module.conductor_hastelloy_fraction,
            sctfcoil_module.conductor_helium_area,
            sctfcoil_module.conductor_helium_fraction,
            sctfcoil_module.conductor_solder_area,
            sctfcoil_module.conductor_solder_fraction,
            sctfcoil_module.conductor_rebco_area,
            sctfcoil_module.conductor_rebco_fraction,
            sctfcoil_module.conductor_critical_current,
        ) = superconductors.croco(
            j_crit_sc,
            sctfcoil_module.conductor_area,
            rebco_variables.croco_od,
            rebco_variables.croco_thick,
        )

        rebco_variables.coppera_m2 = iop / sctfcoil_module.conductor_copper_area

        icrit = sctfcoil_module.conductor_critical_current
        j_crit_cable = (
            sctfcoil_module.croco_strand_critical_current
            / sctfcoil_module.croco_strand_area
        )

        # Critical current density in winding pack
        # a_tf_turn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        j_tf_wp_critical = icrit / a_tf_turn
        #  Ratio of operating / critical current
        iooic = iop / icrit
        #  Operating current density
        jwdgop = iop / a_tf_turn
        #  Actual current density in superconductor,
        # which should be equal to jcrit(thelium+tmarg)

        #  when we have found the desired value of tmarg
        jsc = iooic * j_crit_sc

        # Temperature margin
        current_sharing_t = superconductors.current_sharing_rebco(
            b_tf_inboard_peak, jsc
        )
        tmarg = current_sharing_t - thelium
        tfcoil_variables.temp_margin = (
            tmarg  # Only used in the availabilty routine - see comment to Issue #526
        )

        if output:  # Output ----------------------------------
            total = (
                sctfcoil_module.conductor_copper_area
                + sctfcoil_module.conductor_hastelloy_area
                + sctfcoil_module.conductor_solder_area
                + sctfcoil_module.conductor_jacket_area
                + sctfcoil_module.conductor_helium_area
                + sctfcoil_module.conductor_rebco_area
            )

            if tfcoil_variables.temp_margin <= 0.0e0:
                logger.warning(
                    f"""Negative TFC temperature margin
                temp_margin: {tfcoil_variables.temp_margin}
                b_tf_inboard_peak: {b_tf_inboard_peak}"""
                )

            po.oheadr(self.outfile, "Superconducting TF Coils")
            po.ovarin(self.outfile, "Superconductor switch", "(isumat)", 6)
            po.ocmmnt(
                self.outfile, "Superconductor used: REBCO HTS tape in CroCo strand"
            )

            po.ovarre(
                self.outfile,
                "Thickness of REBCO layer in tape (m)",
                "(rebco_thickness)",
                rebco_variables.rebco_thickness,
            )
            po.ovarre(
                self.outfile,
                "Thickness of copper layer in tape (m)",
                "(copper_thick  )",
                rebco_variables.copper_thick,
            )
            po.ovarre(
                self.outfile,
                "Thickness of Hastelloy layer in tape (m) ",
                "(hastelloy_thickness)",
                rebco_variables.hastelloy_thickness,
            )

            po.ovarre(
                self.outfile,
                "Mean width of tape (m)",
                "(tape_width)",
                rebco_variables.tape_width,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Outer diameter of CroCo copper tube (m) ",
                "(croco_od)",
                rebco_variables.croco_od,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Inner diameter of CroCo copper tube (m) ",
                "(croco_id)",
                rebco_variables.croco_id,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Thickness of CroCo copper tube (m) ",
                "(croco_thick)",
                rebco_variables.croco_thick,
            )

            po.ovarre(
                self.outfile,
                "Thickness of each HTS tape ",
                "(tape_thickness)",
                rebco_variables.tape_thickness,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Thickness of stack of rebco_variables.tapes (m) ",
                "(stack_thickness)",
                rebco_variables.stack_thickness,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Number of rebco_variables.tapes in strand",
                "(tapes)",
                rebco_variables.tapes,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Area of REBCO in strand (m2)",
                "(rebco_area)",
                rebco_variables.rebco_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of copper in strand (m2)",
                "(copper_area)",
                rebco_variables.copper_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of hastelloy substrate in strand (m2) ",
                "(hastelloy_area)",
                rebco_variables.hastelloy_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of solder in strand (m2)  ",
                "(solder_area)",
                rebco_variables.solder_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total: area of CroCo strand (m2)  ",
                "(croco_strand_area)",
                sctfcoil_module.croco_strand_area,
                "OP ",
            )
            if (
                abs(
                    sctfcoil_module.croco_strand_area
                    - (
                        rebco_variables.rebco_area
                        + rebco_variables.copper_area
                        + rebco_variables.hastelloy_area
                        + rebco_variables.solder_area
                    )
                )
                > 1e-6
            ):
                po.ocmmnt(self.outfile, "ERROR: Areas in CroCo strand do not add up")
                logger.warning("Areas in CroCo strand do not add up - see OUT.DAT")

            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Cable information")
            po.ovarin(
                self.outfile,
                "Number of CroCo strands in the cable (fixed) ",
                "",
                6,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total area of cable space (m2)",
                "(a_tf_turn_cable_space_no_void)",
                tfcoil_variables.a_tf_turn_cable_space_no_void,
                "OP ",
            )

            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile,
                "Conductor information (includes jacket, not including insulation)",
            )
            po.ovarre(
                self.outfile,
                "Width of square conductor (cable + steel jacket) (m)",
                "(t_conductor)",
                tfcoil_variables.t_conductor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of conductor (m2)",
                "(area)",
                sctfcoil_module.conductor_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "REBCO area of conductor (mm2)",
                "(rebco_area)",
                sctfcoil_module.conductor_rebco_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of central copper bar (mm2)",
                "(copper_bar_area)",
                sctfcoil_module.conductor_copper_bar_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total copper area of conductor, total (mm2)",
                "(copper_area)",
                sctfcoil_module.conductor_copper_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hastelloy area of conductor (mm2)",
                "(hastelloy_area)",
                sctfcoil_module.conductor_hastelloy_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Solder area of conductor (mm2)",
                "(solder_area)",
                sctfcoil_module.conductor_solder_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Jacket area of conductor (mm2)",
                "(jacket_area)",
                sctfcoil_module.conductor_jacket_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Helium area of conductor (mm2)",
                "(helium_area)",
                sctfcoil_module.conductor_helium_area,
                "OP ",
            )
            if abs(total - sctfcoil_module.conductor_area) > 1e-8:
                po.ovarre(
                    self.outfile,
                    "ERROR: conductor areas do not add up:",
                    "(total)",
                    total,
                    "OP ",
                )
                logger.warning(f"conductor areas do not add up. total: {total}")

            po.ovarre(
                self.outfile,
                "Critical current of CroCo strand (A)",
                "(croco_strand_critical_current)",
                sctfcoil_module.croco_strand_critical_current,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current of conductor (A) ",
                "(conductor_critical_current)",
                sctfcoil_module.conductor_critical_current,
                "OP ",
            )

            if global_variables.run_tests == 1:
                po.oblnkl(self.outfile)
                po.ocmmnt(
                    self.outfile,
                    "PROCESS TF Coil peak field fit. Values for t, z and y:",
                )
                po.oblnkl(self.outfile)
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack width",
                    "(tf_fit_t)",
                    sctfcoil_module.tf_fit_t,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack radial thickness",
                    "(tf_fit_z)",
                    sctfcoil_module.tf_fit_z,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Ratio of actual peak field to nominal axisymmetric peak field",
                    "(tf_fit_y)",
                    sctfcoil_module.tf_fit_y,
                    "OP ",
                )

            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Helium temperature at peak field (= superconductor temperature) (K)",
                "(thelium)",
                thelium,
            )
            po.ovarre(
                self.outfile,
                "Critical current density in superconductor (A/m2)",
                "(j_crit_sc)",
                j_crit_sc,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in cable (A/m2)",
                "(j_crit_cable)",
                j_crit_cable,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in winding pack (A/m2)",
                "(j_tf_wp_critical)",
                j_tf_wp_critical,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current density in winding pack (A/m2)",
                "(jwdgop)",
                jwdgop,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Minimum allowed temperature margin in superconductor (K)",
                "(tmargmin_tf)",
                tfcoil_variables.tmargmin_tf,
            )

            po.ovarre(
                self.outfile,
                "Actual temperature margin in superconductor (K)",
                "(tmarg)",
                tmarg,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Current sharing temperature (K)",
                "(current_sharing_t)",
                current_sharing_t,
                "OP ",
            )
            po.ovarre(self.outfile, "Critical current (A)", "(icrit)", icrit, "OP ")
            po.ovarre(
                self.outfile,
                "Actual current (A)",
                "(c_tf_turn)",
                tfcoil_variables.c_tf_turn,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return j_tf_wp_critical, tmarg

    def supercon(
        self,
        a_tf_turn_cable_space: float,
        a_tf_turn: float,
        b_tf_inboard_peak: float,
        f_a_tf_turn_cooling_extra: float,
        f_a_tf_turn_cable_copper: float,
        c_tf_turn: float,
        j_tf_wp: float,
        i_tf_superconductor: int,
        f_strain_scale: float,
        t_tf_quench_dump: float,
        e_tf_coil_magnetic_stored: float,
        temp_tf_coolant_peak_field: float,
        temp_tf_conductor_peak_quench: float,
        bcritsc: float,
        tcritsc: float,
        output: bool,
    ) -> tuple[float, float, float]:
        """
        Calculates the properties of the TF superconducting conductor.

        :param float a_tf_turn_cable_space:
            Cable space - inside area (m²).
        :param float a_tf_turn:
            Area per turn (i.e. entire jacketed conductor) (m²).
        :param float b_tf_inboard_peak:
            Peak field at conductor (T).
        :param float f_a_tf_turn_cooling_extra:
            Additional fraction of turn cable space reserved for cooling.
        :param float f_a_tf_turn_cable_copper:
            Fraction of conductor that is copper.
        :param float c_tf_turn:
            Operating current per turn (A).
        :param float j_tf_wp:
            Actual winding pack current density (A/m²).
        :param int i_tf_superconductor:
            Switch for conductor type:
            - 1: ITER Nb3Sn, standard parameters
            - 2: Bi-2212 High Temperature Superconductor
            - 3: NbTi
            - 4: ITER Nb3Sn, user-defined parameters
            - 5: WST Nb3Sn parameterisation
            - 7: Durham Ginzburg-Landau Nb-Ti parameterisation
            - 8: Durham Ginzburg-Landau critical surface model for REBCO
            - 9: Hazelton experimental data + Zhai conceptual model for REBCO
        :param float f_strain_scale:
            Adjustment factor (<= 1) to account for strain, radiation damage, fatigue or AC losses.
        :param float t_tf_quench_dump:
            Dump time (s).
        :param float e_tf_coil_magnetic_stored:
            Energy stored in one TF coil (J).
        :param float temp_tf_coolant_peak_field:
            He temperature at peak field point (K).
        :param float temp_tf_conductor_peak_quench:
            Max conductor temperature during quench (K).
        :param float bcritsc:
            Critical field at zero temperature and strain (T) (used only if i_tf_superconductor=4).
        :param float tcritsc:
            Critical temperature at zero field and strain (K) (used only if i_tf_superconductor=4).
        :param bool output:
            Switch for printing output.

        :returns: tuple (j_tf_wp_critical, vd, tmarg)
            - j_tf_wp_critical (float): Critical winding pack current density (A/m²).
            - vd (float): Discharge voltage imposed on a TF coil (V).
            - tmarg (float): Temperature margin (K).

        :notes:
            This routine calculates the superconductor properties for the TF coils.
            It was originally programmed by J. Galambos (1991), from algorithms provided by J. Miller.
            The routine calculates the critical current density (winding pack) and also the protection
            information (for a quench). Not used for the CroCo conductor.

            The critical current density for a superconductor (``j_crit_sc``) is for the superconducting
            strands/tape, not including copper. The critical current density for a cable (``j_crit_cable``)
            accounts for both the fraction of the cable taken up by helium coolant channels, and the cable
            conductor copper fraction (i.e., the copper in the superconducting strands and any additional
            copper, such as REBCO tape support).
        """
        tdump = t_tf_quench_dump

        # Helium channel
        f_a_tf_turn_cable_space_cooling = (
            f_a_tf_turn_cooling_extra
            # Area of inner cooling channel
            + (
                (
                    (np.pi / 4.0e0)
                    * tfcoil_variables.dia_tf_turn_coolant_channel
                    * tfcoil_variables.dia_tf_turn_coolant_channel
                )
                / a_tf_turn_cable_space
            )
        )

        # Guard against negative conductor fraction f_a_tf_turn_cable_space_conductor
        # Kludge to allow solver to continue and hopefully be constrained away
        # from this point
        if f_a_tf_turn_cable_space_cooling > 0.99:
            f_a_tf_turn_cable_space_cooling = 0.99

        #  Conductor fraction (including central helium channel)
        f_a_tf_turn_cable_space_conductor = 1.0e0 - f_a_tf_turn_cable_space_cooling

        if tfcoil_variables.i_str_wp == 0:
            strain = tfcoil_variables.str_tf_con_res
        else:
            strain = tfcoil_variables.str_wp

        # Find critical current density in the superconducter (j_crit_sc)
        # and the superconducting cable (j_crit_cable)

        # =================================================================

        if i_tf_superconductor == 1:  # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            #  j_crit_sc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            j_crit_sc, _, _ = superconductors.itersc(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calculation for costing in $/kAm
            # = Superconducting filaments jc * (1 - strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        # =================================================================

        elif (
            i_tf_superconductor == 2
        ):  # Bi-2212 high temperature superconductor parameterization
            #  Current density in a strand of Bi-2212 conductor
            #  N.B. jcrit returned by superconductors.bi2212 is the critical current density
            #  in the strand, not just the superconducting portion.
            #  The parameterization for j_crit_cable assumes a particular strand
            #  composition that does not require a user-defined copper fraction,
            #  so this is irrelevant in this model
            jstrand = (
                j_tf_wp
                * a_tf_turn
                / (a_tf_turn_cable_space * f_a_tf_turn_cable_space_conductor)
            )

            j_crit_cable, tmarg = superconductors.bi2212(
                b_tf_inboard_peak, jstrand, temp_tf_coolant_peak_field, f_strain_scale
            )
            j_crit_sc = j_crit_cable / (1.0e0 - f_a_tf_turn_cable_copper)
            #  Critical current in cable
            icrit = (
                j_crit_cable * a_tf_turn_cable_space * f_a_tf_turn_cable_space_conductor
            )

            # Strand critical current calulation for costing in $ / kAm
            # Copper in the strand is already accounted for
            tfcoil_variables.j_crit_str_tf = j_crit_sc
        # =================================================================

        elif i_tf_superconductor == 3:  # NbTi data
            bc20m = 15.0e0
            tc0m = 9.3e0
            c0 = 1.0e10
            j_crit_sc, _ = superconductors.jcrit_nbti(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, c0, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        # =================================================================

        elif (
            i_tf_superconductor == 4
        ):  # ITER Nb3Sn parameterization, but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            j_crit_sc, _, _ = superconductors.itersc(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        # =================================================================

        elif i_tf_superconductor == 5:  # WST Nb3Sn parameterisation
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            #  j_crit_sc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            j_crit_sc, _, _ = superconductors.western_superconducting_nb3sn(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        # =================================================================

        elif (
            i_tf_superconductor == 6
        ):  # "REBCO" 2nd generation HTS superconductor in CrCo strand
            raise ProcessValueError(
                "sctfcoil.supercon has been called but tfcoil_variables.i_tf_sc_mat=6"
            )

        elif i_tf_superconductor == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
            bc20m = tfcoil_variables.b_crit_upper_nbti
            tc0m = tfcoil_variables.t_crit_nbti
            j_crit_sc, _, _ = superconductors.gl_nbti(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        # =================================================================

        elif (
            i_tf_superconductor == 8
        ):  # Durham Ginzburg-Landau critical surface model for REBCO
            bc20m = 430
            tc0m = 185
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.7e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.7e-2

            j_crit_sc, _, _ = superconductors.gl_rebco(
                temp_tf_coolant_peak_field, b_tf_inboard_peak, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # Already includes buffer and support layers so no need to include f_a_tf_turn_cable_copper here
            tfcoil_variables.j_crit_str_tf = j_crit_sc

        # =================================================================

        elif (
            i_tf_superconductor == 9
        ):  # Hazelton experimental data + Zhai conceptual model for REBCO
            bc20m = 138
            tc0m = 92
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.7e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.7e-2

            # 'high current density' as per parameterisation described in Wolf,
            #  and based on Hazelton experimental data and Zhai conceptual model;
            #  see subroutine for full references
            j_crit_sc, _, _ = superconductors.hijc_rebco(
                temp_tf_coolant_peak_field,
                b_tf_inboard_peak,
                bc20m,
                tc0m,
                rebco_variables.tape_width,
                rebco_variables.rebco_thickness,
                rebco_variables.tape_thickness,
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = (
                j_crit_sc
                * (1.0e0 - f_a_tf_turn_cable_copper)
                * f_a_tf_turn_cable_space_conductor
            )
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = j_crit_cable * a_tf_turn_cable_space

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (
                1.0e0 - f_a_tf_turn_cable_copper
            )

        else:
            raise ProcessValueError(
                "Illegal value for i_tf_sc_mat", i_tf_superconductor=i_tf_superconductor
            )

        # Critical current density in winding pack
        # a_tf_turn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        j_tf_wp_critical = icrit / a_tf_turn
        #  Ratio of operating / critical current
        iooic = c_tf_turn / icrit
        #  Operating current density
        j_tf_coil_turn = c_tf_turn / a_tf_turn
        #  Actual current density in superconductor, which should be equal to jcrit(temp_tf_coolant_peak_field+tmarg)
        #  when we have found the desired value of tmarg
        jsc = iooic * j_crit_sc

        if iooic <= 0e0:
            logger.warning(
                f"""Negative Iop/Icrit for TF coil
            jsc: {jsc}
            iooic: {iooic}
            j_crit_sc: {j_crit_sc}
            Check conductor dimensions. Cable space area a_tf_turn_cable_space likely gone negative. a_tf_turn_cable_space: {a_tf_turn_cable_space}
            This is likely because dr_tf_turn_cable_space or dx_tf_turn_cable_space has gone negative:
            dr_tf_turn_cable_space: {sctfcoil_module.dr_tf_turn_cable_space}
            dx_tf_turn_cable_space: {sctfcoil_module.dx_tf_turn_cable_space}
            """
            )

        # REBCO measurements from 2 T to 14 T, extrapolating outside this
        if (i_tf_superconductor == 8) and (tfcoil_variables.bmaxtfrp >= 14):
            error_handling.report_error(266)

        #  Temperature margin (already calculated in superconductors.bi2212 for i_tf_superconductor=2)

        if i_tf_superconductor in (
            1,
            3,
            4,
            5,
            7,
            8,
            9,
        ):  # Find temperature at which current density margin = 0
            if i_tf_superconductor == 3:
                arguments = (
                    i_tf_superconductor,
                    jsc,
                    b_tf_inboard_peak,
                    strain,
                    bc20m,
                    tc0m,
                    c0,
                )
            else:
                arguments = (
                    i_tf_superconductor,
                    jsc,
                    b_tf_inboard_peak,
                    strain,
                    bc20m,
                    tc0m,
                )

            another_estimate = 2 * temp_tf_coolant_peak_field
            t_zero_margin, root_result = optimize.newton(
                superconductors.current_density_margin,
                temp_tf_coolant_peak_field,
                fprime=None,
                args=arguments,
                # args=(i_tf_superconductor, jsc, b_tf_inboard_peak, strain, bc20m, tc0m,),
                tol=1.0e-06,
                maxiter=50,
                fprime2=None,
                x1=another_estimate,
                rtol=1.0e-6,
                full_output=True,
                disp=True,
            )
            # print(root_result)  # Diagnostic for newton method
            tmarg = t_zero_margin - temp_tf_coolant_peak_field
            tfcoil_variables.temp_margin = tmarg

        # Find the current density limited by the protection limit
        # At present only valid for LTS windings (Nb3Sn properties assumed)
        tfcoil_variables.jwdgpro, vd = self.protect(
            c_tf_turn,
            e_tf_coil_magnetic_stored,
            a_tf_turn_cable_space,
            a_tf_turn,
            tdump,
            f_a_tf_turn_cable_space_conductor,
            f_a_tf_turn_cable_copper,
            temp_tf_coolant_peak_field,
            temp_tf_conductor_peak_quench,
            b_tf_inboard_peak,
            tfcoil_variables.rrr_tf_cu,
            tfcoil_variables.t_tf_quench_detection,
            constraint_variables.nflutfmax,
        )

        if output:  # Output --------------------------
            if tmarg <= 0.0e0:
                logger.warning(
                    """Negative TFC temperature margin
                tmarg: {tmarg}
                b_tf_inboard_peak: {b_tf_inboard_peak}
                jcrit0: {jcrit0}
                jsc: {jsc}
                """
                )

            po.oheadr(self.outfile, "Superconducting TF Coils")
            po.ovarin(
                self.outfile,
                "Superconductor switch",
                "(i_tf_superconductor)",
                i_tf_superconductor,
            )

            if i_tf_superconductor == 1:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, "  (ITER Jcrit model, standard parameters)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 2:
                po.ocmmnt(self.outfile, "Superconductor used: Bi-2212 HTS")
            if i_tf_superconductor == 3:
                po.ocmmnt(self.outfile, "Superconductor used: NbTi")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 4:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, "  (ITER Jcrit model, user-defined parameters)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 5:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 7:
                po.ocmmnt(self.outfile, "Superconductor used: Nb-Ti")
                po.ocmmnt(
                    self.outfile, " (Durham Ginzburg-Landau critical surface model)"
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 8:
                po.ocmmnt(self.outfile, "Superconductor used: REBCO")
                po.ocmmnt(
                    self.outfile, " (Durham Ginzburg-Landau critical surface model)"
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if i_tf_superconductor == 9:
                po.ocmmnt(self.outfile, "Superconductor used: REBCO")
                po.ocmmnt(
                    self.outfile,
                    " (Hazelton experimental data + Zhai conceptual model)",
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )

            if global_variables.run_tests == 1:
                po.oblnkl(self.outfile)
                po.ocmmnt(
                    self.outfile,
                    "PROCESS TF Coil peak field fit. Values for t, z and y:",
                )
                po.oblnkl(self.outfile)
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack width",
                    "(tf_fit_t)",
                    sctfcoil_module.tf_fit_t,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack radial thickness",
                    "(tf_fit_z)",
                    sctfcoil_module.tf_fit_z,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Ratio of peak field with ripple to nominal axisymmetric peak field",
                    "(tf_fit_y)",
                    sctfcoil_module.tf_fit_y,
                    "OP ",
                )

            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Helium temperature at peak field (= superconductor temperature) (K)",
                "(temp_tf_coolant_peak_field)",
                temp_tf_coolant_peak_field,
            )
            po.ovarre(
                self.outfile,
                "Total helium fraction inside cable space",
                "(f_a_tf_turn_cable_space_cooling)",
                f_a_tf_turn_cable_space_cooling,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Copper fraction of conductor",
                "(fcutfsu)",
                f_a_tf_turn_cable_copper,
            )
            po.ovarre(
                self.outfile,
                "Residual manufacturing strain on superconductor",
                "(str_tf_con_res)",
                tfcoil_variables.str_tf_con_res,
            )
            po.ovarre(
                self.outfile,
                "Self-consistent strain on superconductor",
                "(str_wp)",
                tfcoil_variables.str_wp,
            )
            po.ovarre(
                self.outfile,
                "Critical current density in superconductor (A/m2)",
                "(j_crit_sc)",
                j_crit_sc,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in cable (A/m2)",
                "(j_crit_cable)",
                j_crit_cable,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in winding pack (A/m2)",
                "(j_tf_wp_critical)",
                j_tf_wp_critical,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current density in winding pack (A/m2)",
                "(j_tf_coil_turn)",
                j_tf_coil_turn,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Minimum allowed temperature margin in superconductor (K)",
                "(tmargmin_tf)",
                tfcoil_variables.tmargmin_tf,
            )
            po.ovarre(
                self.outfile,
                "Actual temperature margin in superconductor (K)",
                "(tmarg)",
                tmarg,
                "OP ",
            )
            po.ovarre(self.outfile, "Critical current (A)", "(icrit)", icrit, "OP ")
            po.ovarre(
                self.outfile,
                "Actual current (A)",
                "(c_tf_turn)",
                tfcoil_variables.c_tf_turn,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return j_tf_wp_critical, vd, tmarg

    def protect(
        self,
        aio,
        tfes,
        acs,
        aturn,
        tdump,
        fcond,
        fcu,
        tba,
        tmax,
        peak_field,
        cu_rrr,
        detection_time,
        fluence,
    ):
        """
        Calculates the maximum conductor current density limited by the protection limit,
        and the discharge voltage for a TF coil.

        :param float aio: Operating current (A)
        :param float tfes: Energy stored in one TF coil (J)
        :param float acs: Cable space - inside area (m²)
        :param float aturn: Area per turn (i.e. entire cable) (m²)
        :param float tdump: Dump time (s)
        :param float fcond: Fraction of cable space containing conductor
        :param float fcu: Fraction of conductor that is copper
        :param float tba: Helium temperature at peak field point (K)
        :param float tmax: Maximum conductor temperature during quench (K)
        :param float peak_field: Maximum conductor temperature during quench (T)
        :param float cur_rrr: Copper residual-resistance-ratio
        :param float detection_time: Quench detection time (s)
        :param float fluence: End-of-life neutron fluence in the copper (1/m²)

        :return float ajwpro: Winding pack current density from temperature rise protection (A/m²)
        :return float vd: Discharge voltage imposed on a TF coil (V)
        """
        #  Dump voltage
        vd = 2.0e0 * tfes / (tdump * aio)
        ajwpro = (
            acs
            / aturn
            * calculate_quench_protection_current_density(
                tdump,
                peak_field,
                fcu,
                1.0 - fcond,
                tba,
                tmax,
                cu_rrr,
                detection_time,
                fluence,
            )
        )

        return ajwpro, vd

    def vv_stress_on_quench(self):
        """Calculate the Tresca stress [Pa] of the Vacuum Vessel (VV)
        experienced when the TF coil quenches.

        Author: Timothy Nunn, UKAEA

        Assumes the current center line (CCL) of the TF coil is the
        middle of the coil.

        We assume vertical symmetry which is only true for double null
        machines.
        """
        H_coil = build_variables.z_tf_inside_half + (build_variables.dr_tf_inboard / 2)
        ri_coil = build_variables.r_tf_inboard_mid
        ro_coil = build_variables.r_tf_outboard_mid
        # NOTE: rm is measured from the outside edge of the coil because thats where
        # the radius of the first ellipse is measured from
        rm_coil = build_variables.r_tf_inboard_out + tfcoil_variables.tfa[0]

        H_vv = (
            build_variables.z_plasma_xpoint_upper
            + build_variables.dz_xpoint_divertor
            + divertor_variables.dz_divertor
            + build_variables.dz_shld_upper
            + (build_variables.dz_vv_upper / 2)
        )
        # ri and ro for VV dont consider the shield widths
        # because it is assumed the shield is on the plasma side
        # of the VV
        ri_vv = build_variables.r_vv_inboard_out - (build_variables.dr_vv_outboard / 2)
        ro_vv = (
            build_variables.r_tf_outboard_mid
            - (build_variables.dr_tf_outboard / 2)
            - build_variables.dr_tf_shld_gap
            - build_variables.dr_shld_thermal_outboard
            - build_variables.dr_shld_vv_gap_outboard
            - (build_variables.dr_vv_outboard / 2)
        )

        # Assume the radius of the first ellipse of the VV is in the same proportion to
        # that of the plasma facing radii of the two structures
        tf_vv_frac = build_variables.r_tf_inboard_out / build_variables.r_vv_inboard_out
        rm_vv = build_variables.r_vv_inboard_out + (
            tfcoil_variables.tfa[0] * tf_vv_frac
        )

        sctfcoil_module.vv_stress_quench = vv_stress_on_quench(
            # TF shape
            H_coil=H_coil,
            ri_coil=ri_coil,
            ro_coil=ro_coil,
            rm_coil=rm_coil,
            ccl_length_coil=tfcoil_variables.len_tf_coil,
            theta1_coil=tfcoil_variables.theta1_coil,
            # VV shape
            H_vv=H_vv,
            ri_vv=ri_vv,
            ro_vv=ro_vv,
            rm_vv=rm_vv,
            theta1_vv=tfcoil_variables.theta1_vv,
            # TF properties
            n_tf_coils=tfcoil_variables.n_tf_coils,
            n_tf_coil_turns=tfcoil_variables.n_tf_coil_turns,
            # Area of the radial plate taken to be the area of steel in the WP
            # TODO: value clipped due to #1883
            s_rp=np.clip(sctfcoil_module.a_tf_coil_inboard_steel, 0, None),
            s_cc=sctfcoil_module.a_tf_plasma_case
            + sctfcoil_module.a_tf_coil_nose_case
            + 2.0 * sctfcoil_module.dx_tf_side_case_average,
            taud=tfcoil_variables.tdmptf,
            # TODO: is this the correct current?
            i_op=sctfcoil_module.c_tf_coil / tfcoil_variables.n_tf_coil_turns,
            # VV properties
            d_vv=build_variables.dr_vv_shells,
        )

    def peak_tf_with_ripple(
        self,
        n_tf_coils,
        dx_tf_wp_primary_toroidal,
        dr_tf_wp_with_insulation,
        tfin,
        b_tf_inboard_peak,
    ):
        """Peak toroidal field on the conductor
        author: P J Knight, CCFE, Culham Science Centre
        This subroutine calculates the peak toroidal field at the
        outboard edge of the inboard TF coil winding pack, including
        the effects of ripple.
        <P>For 16, 18 or 20 coils, the calculation uses fitting formulae
        derived by M. Kovari using MAGINT calculations on coil sets based
        on a DEMO1 case.
        <P>For other numbers of coils, the original estimate using a 9%
        increase due to ripple from the axisymmetric calculation is used.
        M. Kovari, Toroidal Field Coils - Maximum Field and Ripple -
        Parametric Calculation, July 2014

        :param n_tf_coils: number of TF coils
        :type n_tf_coils: float
        :param dx_tf_wp_primary_toroidal: width of plasma-facing face of winding pack (m)
        :type dx_tf_wp_primary_toroidal: float
        :param dr_tf_wp_with_insulation: radial thickness of winding pack (m)
        :type dr_tf_wp_with_insulation: float
        :param tfin: major radius of centre of winding pack (m)
        :type tfin: float
        :param b_tf_inboard_peak: nominal (axisymmetric) peak toroidal field (T)
        :type b_tf_inboard_peak: float

        :returns: (bmaxtfrp, flag)
        * bmaxtfrp: peak toroidal field including ripple (T)
        * flag: flag warning of applicability problems

        :rtype: Tuple[float, int]

        """
        a = np.zeros((4,))
        flag = 0

        #  Set fitting coefficients for different numbers of TF coils

        int_n_tf = np.round(n_tf_coils)

        if int_n_tf == 16:
            a[0] = 0.28101e0
            a[1] = 1.8481e0
            a[2] = -0.88159e0
            a[3] = 0.93834e0
        elif int_n_tf == 18:
            a[0] = 0.29153e0
            a[1] = 1.81600e0
            a[2] = -0.84178e0
            a[3] = 0.90426e0
        elif int_n_tf == 20:
            a[0] = 0.29853e0
            a[1] = 1.82130e0
            a[2] = -0.85031e0
            a[3] = 0.89808e0

        else:
            bmaxtfrp = 1.09e0 * b_tf_inboard_peak
            return bmaxtfrp, flag

        #  Maximum winding pack width before adjacent packs touch
        #  (ignoring the external case and ground wall thicknesses)

        wmax = (2.0e0 * tfin + dr_tf_wp_with_insulation) * np.tan(np.pi / n_tf_coils)

        #  Dimensionless winding pack width

        sctfcoil_module.tf_fit_t = dx_tf_wp_primary_toroidal / wmax
        if (sctfcoil_module.tf_fit_t < 0.3e0) or (sctfcoil_module.tf_fit_t > 1.1e0):
            # write(*,*) 'PEAK_TF_WITH_RIPPLE: fitting problem; t = ',t
            flag = 1

        #  Dimensionless winding pack radial thickness

        sctfcoil_module.tf_fit_z = dr_tf_wp_with_insulation / wmax
        if (sctfcoil_module.tf_fit_z < 0.26e0) or (sctfcoil_module.tf_fit_z > 0.7e0):
            # write(*,*) 'PEAK_TF_WITH_RIPPLE: fitting problem; z = ',z
            flag = 2

        #  Ratio of peak field with ripple to nominal axisymmetric peak field

        sctfcoil_module.tf_fit_y = (
            a[0]
            + a[1] * np.exp(-sctfcoil_module.tf_fit_t)
            + a[2] * sctfcoil_module.tf_fit_z
            + a[3] * sctfcoil_module.tf_fit_z * sctfcoil_module.tf_fit_t
        )

        bmaxtfrp = sctfcoil_module.tf_fit_y * b_tf_inboard_peak

        return bmaxtfrp, flag

    def sc_tf_internal_geom(self, i_tf_wp_geom, i_tf_case_geom, i_tf_turns_integer):
        """
        Author : S. Kahn, CCFE
        Seting the WP, case and turns geometry for SC magnets
        """

        # Calculating the WP / ground insulation areas
        (
            sctfcoil_module.r_tf_wp_inboard_inner,
            sctfcoil_module.r_tf_wp_inboard_outer,
            sctfcoil_module.r_tf_wp_inboard_centre,
            sctfcoil_module.dx_tf_wp_toroidal_min,
            tfcoil_variables.dx_tf_wp_primary_toroidal,
            tfcoil_variables.dx_tf_wp_secondary_toroidal,
            sctfcoil_module.dx_tf_wp_toroidal_average,
            sctfcoil_module.a_tf_wp_with_insulation,
            sctfcoil_module.a_tf_wp_no_insulation,
            sctfcoil_module.a_tf_wp_ground_insulation,
        ) = self.superconducting_tf_wp_geometry(
            i_tf_wp_geom=i_tf_wp_geom,
            r_tf_inboard_in=build_variables.r_tf_inboard_in,
            dr_tf_nose_case=tfcoil_variables.dr_tf_nose_case,
            dr_tf_wp_with_insulation=tfcoil_variables.dr_tf_wp_with_insulation,
            tan_theta_coil=sctfcoil_module.tan_theta_coil,
            dx_tf_side_case_min=tfcoil_variables.dx_tf_side_case_min,
            dx_tf_wp_insulation=tfcoil_variables.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=tfcoil_variables.dx_tf_wp_insertion_gap,
        )

        # Calculating the TF steel casing areas
        (
            tfcoil_variables.a_tf_coil_inboard_case,
            tfcoil_variables.a_tf_coil_outboard_case,
            sctfcoil_module.a_tf_plasma_case,
            sctfcoil_module.a_tf_coil_nose_case,
            sctfcoil_module.dx_tf_side_case_average,
        ) = self.superconducting_tf_case_geometry(
            i_tf_case_geom=i_tf_case_geom,
            i_tf_wp_geom=i_tf_wp_geom,
            a_tf_inboard_total=tfcoil_variables.a_tf_inboard_total,
            n_tf_coils=tfcoil_variables.n_tf_coils,
            a_tf_wp_with_insulation=sctfcoil_module.a_tf_wp_with_insulation,
            a_tf_leg_outboard=tfcoil_variables.a_tf_leg_outboard,
            rad_tf_coil_inboard_toroidal_half=sctfcoil_module.rad_tf_coil_inboard_toroidal_half,
            r_tf_inboard_out=build_variables.r_tf_inboard_out,
            tan_theta_coil=sctfcoil_module.tan_theta_coil,
            r_tf_wp_inboard_outer=sctfcoil_module.r_tf_wp_inboard_outer,
            dr_tf_plasma_case=tfcoil_variables.dr_tf_plasma_case,
            r_tf_wp_inboard_inner=sctfcoil_module.r_tf_wp_inboard_inner,
            r_tf_inboard_in=build_variables.r_tf_inboard_in,
            dx_tf_side_case_min=tfcoil_variables.dx_tf_side_case_min,
            dr_tf_wp_with_insulation=tfcoil_variables.dr_tf_wp_with_insulation,
        )

        # WP/trun currents
        self.tf_wp_currents()

        # Setting the WP turn geometry / areas
        if i_tf_turns_integer == 0:
            # Non-ingeger number of turns
            (
                tfcoil_variables.a_tf_turn_cable_space_no_void,
                tfcoil_variables.a_tf_turn_steel,
                tfcoil_variables.a_tf_turn_insulation,
                tfcoil_variables.n_tf_coil_turns,
            ) = self.tf_averaged_turn_geom(
                tfcoil_variables.j_tf_wp,
                tfcoil_variables.dx_tf_turn_steel,
                tfcoil_variables.dx_tf_turn_insulation,
                tfcoil_variables.i_tf_sc_mat,
            )

        else:
            # Integer number of turns
            (
                tfcoil_variables.a_tf_turn_cable_space_no_void,
                tfcoil_variables.a_tf_turn_steel,
                tfcoil_variables.a_tf_turn_insulation,
                tfcoil_variables.c_tf_turn,
                tfcoil_variables.n_tf_coil_turns,
            ) = self.tf_integer_turn_geom(
                tfcoil_variables.n_layer,
                tfcoil_variables.n_pancake,
                tfcoil_variables.dx_tf_turn_steel,
                tfcoil_variables.dx_tf_turn_insulation,
            )

        # Areas and fractions
        # -------------------
        # Central helium channel down the conductor core [m2]
        tfcoil_variables.a_tf_wp_coolant_channels = (
            0.25e0
            * tfcoil_variables.n_tf_coil_turns
            * np.pi
            * tfcoil_variables.dia_tf_turn_coolant_channel**2
        )

        # Total conductor cross-sectional area, taking account of void area
        # and central helium channel [m2]
        tfcoil_variables.a_tf_wp_conductor = (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * tfcoil_variables.n_tf_coil_turns
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
            - tfcoil_variables.a_tf_wp_coolant_channels
        )

        # Void area in conductor for He, not including central channel [m2]
        tfcoil_variables.a_tf_wp_extra_void = (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * tfcoil_variables.n_tf_coil_turns
            * tfcoil_variables.f_a_tf_turn_cable_space_extra_void
        )

        # Area of inter-turn insulation: total [m2]
        tfcoil_variables.a_tf_coil_wp_turn_insulation = (
            tfcoil_variables.n_tf_coil_turns * tfcoil_variables.a_tf_turn_insulation
        )

        # Area of steel structure in winding pack [m2]
        tfcoil_variables.a_tf_wp_steel = (
            tfcoil_variables.n_tf_coil_turns * tfcoil_variables.a_tf_turn_steel
        )

        # Inboard coil steel area [m2]
        sctfcoil_module.a_tf_coil_inboard_steel = (
            tfcoil_variables.a_tf_coil_inboard_case + tfcoil_variables.a_tf_wp_steel
        )

        # Inboard coil steel fraction [-]
        sctfcoil_module.f_a_tf_coil_inboard_steel = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_coil_inboard_steel
            / tfcoil_variables.a_tf_inboard_total
        )

        # Inboard coil insulation cross-section [m2]
        sctfcoil_module.a_tf_coil_inboard_insulation = (
            tfcoil_variables.a_tf_coil_wp_turn_insulation
            + sctfcoil_module.a_tf_wp_ground_insulation
        )

        #  Inboard coil insulation fraction [-]
        sctfcoil_module.f_a_tf_coil_inboard_insulation = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_coil_inboard_insulation
            / tfcoil_variables.a_tf_inboard_total
        )

        # Negative areas or fractions error reporting
        if (
            tfcoil_variables.a_tf_wp_conductor <= 0.0e0
            or tfcoil_variables.a_tf_wp_extra_void <= 0.0e0
            or tfcoil_variables.a_tf_coil_wp_turn_insulation <= 0.0e0
            or tfcoil_variables.a_tf_wp_steel <= 0.0e0
            or sctfcoil_module.a_tf_coil_inboard_steel <= 0.0e0
            or sctfcoil_module.f_a_tf_coil_inboard_steel <= 0.0e0
            or sctfcoil_module.a_tf_coil_inboard_insulation <= 0.0e0
            or sctfcoil_module.f_a_tf_coil_inboard_insulation <= 0.0e0
        ):
            error_handling.fdiags[0] = tfcoil_variables.a_tf_wp_conductor
            error_handling.fdiags[1] = tfcoil_variables.a_tf_wp_extra_void
            error_handling.fdiags[2] = tfcoil_variables.a_tf_coil_wp_turn_insulation
            error_handling.fdiags[3] = tfcoil_variables.a_tf_wp_steel
            error_handling.fdiags[4] = sctfcoil_module.a_tf_coil_inboard_steel
            error_handling.fdiags[5] = sctfcoil_module.f_a_tf_coil_inboard_steel
            error_handling.fdiags[6] = sctfcoil_module.a_tf_coil_inboard_insulation
            error_handling.fdiags[7] = sctfcoil_module.f_a_tf_coil_inboard_insulation
            error_handling.report_error(276)

    def superconducting_tf_wp_geometry(
        self,
        i_tf_wp_geom: int,
        r_tf_inboard_in: float,
        dr_tf_nose_case: float,
        dr_tf_wp_with_insulation: float,
        tan_theta_coil: float,
        dx_tf_side_case_min: float,
        dx_tf_wp_insulation: float,
        dx_tf_wp_insertion_gap: float,
    ) -> tuple[
        float,  # r_tf_wp_inboard_inner
        float,  # r_tf_wp_inboard_outer
        float,  # r_tf_wp_inboard_centre
        float,  # dx_tf_wp_toroidal_min
        float,  # dx_tf_wp_primary_toroidal
        float,  # dx_tf_wp_secondary_toroidal
        float,  # dx_tf_wp_toroidal_average
        float,  # a_tf_wp_with_insulation
        float,  # a_tf_wp_no_insulation
        float,  # a_tf_wp_ground_insulation
    ]:
        """
        Calculates the winding pack (WP) geometry and cross-sectional areas for superconducting toroidal field (TF) coils.

        :param i_tf_wp_geom:
            - 0: Rectangular
            - 1: Double rectangular
            - 2: Trapezoidal
        :type i_tf_wp_geom: int
        :param r_tf_inboard_in: Inboard inner radius [m].
        :type r_tf_inboard_in: float
        :param dr_tf_nose_case: Radial thickness of nose case [m].
        :type dr_tf_nose_case: float
        :param dr_tf_wp_with_insulation: Radial thickness of winding pack including insulation [m].
        :type dr_tf_wp_with_insulation: float
        :param tan_theta_coil: Tangent of coil half angle [-].
        :type tan_theta_coil: float
        :param dx_tf_side_case_min: Side case thickness [m].
        :type dx_tf_side_case_min: float
        :param dx_tf_wp_insulation: Insulation thickness [m].
        :type dx_tf_wp_insulation: float
        :param dx_tf_wp_insertion_gap: Insertion gap thickness [m].
        :type dx_tf_wp_insertion_gap: float

        :returns:
            Tuple containing:
            - r_tf_wp_inboard_inner (float): WP inboard inner radius [m]
            - r_tf_wp_inboard_outer (float): WP inboard outer radius [m]
            - r_tf_wp_inboard_centre (float): WP inboard centre radius [m]
            - dx_tf_wp_toroidal_min (float): Minimal toroidal thickness of WP [m]
            - dx_tf_wp_primary_toroidal (float): Primary toroidal thickness [m]
            - dx_tf_wp_secondary_toroidal (float): Secondary toroidal thickness [m]
            - dx_tf_wp_toroidal_average (float): Averaged toroidal thickness [m]
            - a_tf_wp_with_insulation (float): WP cross-sectional area with insulation [m²]
            - a_tf_wp_no_insulation (float): WP cross-sectional area without insulation [m²]
            - a_tf_wp_ground_insulation (float): WP ground insulation cross-sectional area [m²]
        :rtype: tuple[float, float, float, float, float, float, float, float, float, float]

        :raises ValueError: If calculated winding pack area (with or without insulation) is non-positive.
        """

        r_tf_wp_inboard_inner = r_tf_inboard_in + dr_tf_nose_case

        # Radial position of outer edge of winding pack [m]
        r_tf_wp_inboard_outer = r_tf_wp_inboard_inner + dr_tf_wp_with_insulation

        # Radius of geometrical centre of winding pack [m]
        r_tf_wp_inboard_centre = 0.5e0 * (r_tf_wp_inboard_inner + r_tf_wp_inboard_outer)

        # TF toroidal thickness at the WP inner radius [m]
        dx_tf_wp_inner_toroidal = 2.0e0 * r_tf_wp_inboard_inner * tan_theta_coil

        # Minimal toroidal thickness of winding pack [m]
        dx_tf_wp_toroidal_min = dx_tf_wp_inner_toroidal - 2.0e0 * dx_tf_side_case_min

        # Rectangular WP
        # --------------
        if i_tf_wp_geom == 0:
            # Outer WP layer toroidal thickness [m]
            dx_tf_wp_primary_toroidal = dx_tf_wp_toroidal_min

            # No secondary WP here but will set for consistency
            dx_tf_wp_secondary_toroidal = dx_tf_wp_toroidal_min

            # Averaged toroidal thickness of of winding pack [m]
            dx_tf_wp_toroidal_average = dx_tf_wp_toroidal_min

            # Total cross-sectional area of winding pack [m²]
            a_tf_wp_with_insulation = (
                dr_tf_wp_with_insulation * dx_tf_wp_primary_toroidal
            )

            # WP cross-section without insertion gap and ground insulation [m²]
            a_tf_wp_no_insulation = (
                dr_tf_wp_with_insulation
                - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
            ) * (
                dx_tf_wp_primary_toroidal
                - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
            )

            # Cross-section area of the WP ground insulation [m²]
            a_tf_wp_ground_insulation = (
                dr_tf_wp_with_insulation - 2.0e0 * dx_tf_wp_insertion_gap
            ) * (
                dx_tf_wp_primary_toroidal - 2.0e0 * dx_tf_wp_insertion_gap
            ) - a_tf_wp_no_insulation

        # Double rectangular WP
        # ---------------------
        elif i_tf_wp_geom == 1:
            # Thickness of winding pack section at R > sctfcoil_module.r_tf_wp_inboard_centre [m]
            dx_tf_wp_primary_toroidal = 2.0e0 * (
                r_tf_wp_inboard_centre * tan_theta_coil - dx_tf_side_case_min
            )

            # Thickness of winding pack section at R < sctfcoil_module.r_tf_wp_inboard_centre [m]
            dx_tf_wp_secondary_toroidal = 2.0e0 * (
                r_tf_wp_inboard_inner * tan_theta_coil - dx_tf_side_case_min
            )

            # Averaged toroidal thickness of of winding pack [m]
            dx_tf_wp_toroidal_average = 0.5e0 * (
                dx_tf_wp_primary_toroidal + dx_tf_wp_secondary_toroidal
            )

            # Total cross-sectional area of winding pack [m²]
            # Including ground insulation and insertion gap
            a_tf_wp_with_insulation = (
                dr_tf_wp_with_insulation * dx_tf_wp_toroidal_average
            )

            # WP cross-section without insertion gap and ground insulation [m²]
            a_tf_wp_no_insulation = (
                0.5e0
                * (
                    dr_tf_wp_with_insulation
                    - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
                )
                * (
                    dx_tf_wp_primary_toroidal
                    + dx_tf_wp_secondary_toroidal
                    - 4.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
                )
            )

            # Cross-section area of the WP ground insulation [m²]
            a_tf_wp_ground_insulation = (
                0.5e0
                * (dr_tf_wp_with_insulation - 2.0e0 * dx_tf_wp_insertion_gap)
                * (
                    dx_tf_wp_primary_toroidal
                    + dx_tf_wp_secondary_toroidal
                    - 4.0e0 * dx_tf_wp_insertion_gap
                )
                - a_tf_wp_no_insulation
            )

        # Trapezoidal WP
        # --------------
        else:
            # Thickness of winding pack section at r_tf_wp_inboard_outer [m]
            dx_tf_wp_primary_toroidal = 2.0e0 * (
                r_tf_wp_inboard_outer * tan_theta_coil - dx_tf_side_case_min
            )

            # Thickness of winding pack section at r_tf_wp_inboard_inner [m]
            dx_tf_wp_secondary_toroidal = 2.0e0 * (
                r_tf_wp_inboard_inner * tan_theta_coil - dx_tf_side_case_min
            )

            # Averaged toroidal thickness of of winding pack [m]
            dx_tf_wp_toroidal_average = 0.5e0 * (
                dx_tf_wp_primary_toroidal + dx_tf_wp_secondary_toroidal
            )

            # Total cross-sectional area of winding pack [m²]
            # Including ground insulation and insertion gap
            a_tf_wp_with_insulation = (
                dr_tf_wp_with_insulation
                * 0.5
                * (dx_tf_wp_primary_toroidal + dx_tf_wp_secondary_toroidal)
            )

            # WP cross-section without insertion gap and ground insulation [m²]
            a_tf_wp_no_insulation = (
                (
                    dr_tf_wp_with_insulation
                    - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
                )
                * (
                    (
                        dx_tf_wp_secondary_toroidal
                        - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
                    )
                    + (
                        dx_tf_wp_primary_toroidal
                        - 2.0e0 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)
                    )
                )
                / 2
            )

            # Cross-section area of the WP ground insulation [m²]
            a_tf_wp_ground_insulation = (
                dr_tf_wp_with_insulation - 2.0e0 * dx_tf_wp_insertion_gap
            ) * (
                (
                    (dx_tf_wp_primary_toroidal - 2.0e0 * dx_tf_wp_insertion_gap)
                    + (dx_tf_wp_secondary_toroidal - 2.0e0 * dx_tf_wp_insertion_gap)
                )
                / 2
            ) - a_tf_wp_no_insulation

        # --------------
        # Negative WP area error reporting
        if a_tf_wp_no_insulation <= 0.0e0 or a_tf_wp_with_insulation <= 0.0e0:
            error_handling.fdiags[0] = a_tf_wp_no_insulation
            error_handling.fdiags[1] = a_tf_wp_with_insulation
            error_handling.report_error(99)

        return (
            r_tf_wp_inboard_inner,
            r_tf_wp_inboard_outer,
            r_tf_wp_inboard_centre,
            dx_tf_wp_toroidal_min,
            dx_tf_wp_primary_toroidal,
            dx_tf_wp_secondary_toroidal,
            dx_tf_wp_toroidal_average,
            a_tf_wp_with_insulation,
            a_tf_wp_no_insulation,
            a_tf_wp_ground_insulation,
        )

    def superconducting_tf_case_geometry(
        self,
        i_tf_wp_geom: int,
        i_tf_case_geom: int,
        a_tf_inboard_total: float,
        n_tf_coils: float,
        a_tf_wp_with_insulation: float,
        a_tf_leg_outboard: float,
        rad_tf_coil_inboard_toroidal_half: float,
        r_tf_inboard_out: float,
        tan_theta_coil: float,
        r_tf_wp_inboard_outer: float,
        dr_tf_plasma_case: float,
        r_tf_wp_inboard_inner: float,
        r_tf_inboard_in: float,
        dx_tf_side_case_min: float,
        dr_tf_wp_with_insulation: float,
    ) -> tuple[float, float, float, float, float]:
        """
        Setting the case geometry and area for SC magnets

        :param i_tf_wp_geom: Index specifying winding pack geometry (0: rectangular, 1: double rectangular, else: trapezoidal).
        :type i_tf_wp_geom: int
        :param i_tf_case_geom: Index specifying case geometry (0: circular, else: straight).
        :type i_tf_case_geom: int
        :param a_tf_inboard_total: Total inboard area for TF coils [m²].
        :type a_tf_inboard_total: float
        :param n_tf_coils: Number of TF coils.
        :type n_tf_coils: float
        :param a_tf_wp_with_insulation: Area of winding pack with insulation [m²].
        :type a_tf_wp_with_insulation: float
        :param a_tf_leg_outboard: Outboard leg cross-sectional area [m²].
        :type a_tf_leg_outboard: float
        :param rad_tf_coil_inboard_toroidal_half: Half toroidal radius of inboard coil [m].
        :type rad_tf_coil_inboard_toroidal_half: float
        :param r_tf_inboard_out: Outer radius of inboard TF coil [m].
        :type r_tf_inboard_out: float
        :param tan_theta_coil: Tangent of coil angle theta.
        :type tan_theta_coil: float
        :param r_tf_wp_inboard_outer: Outer radius of inboard winding pack [m].
        :type r_tf_wp_inboard_outer: float
        :param dr_tf_plasma_case: Radial thickness of plasma case [m].
        :type dr_tf_plasma_case: float
        :param r_tf_wp_inboard_inner: Inner radius of inboard winding pack [m].
        :type r_tf_wp_inboard_inner: float
        :param r_tf_inboard_in: Inner radius of inboard TF coil [m].
        :type r_tf_inboard_in: float
        :param dx_tf_side_case_min: Minimum lateral casing thickness [m].
        :type dx_tf_side_case_min: float
        :param dr_tf_wp_with_insulation: Radial thickness of winding pack with insulation [m].
        :type dr_tf_wp_with_insulation: float

        :returns: Tuple containing:
            - a_tf_coil_inboard_case (float): Inboard case area [m²].
            - a_tf_coil_outboard_case (float): Outboard case area [m²].
            - a_tf_plasma_case (float): Front casing area [m²].
            - a_tf_coil_nose_case (float): Nose casing area [m²].
            - dx_tf_side_case_average (float): Average lateral casing thickness [m].
        :rtype: tuple[float, float, float, float, float]

        :raises: Reports error if calculated casing areas are negative.
        """

        # Total area of inboard TF coil case [m²]
        a_tf_coil_inboard_case = (
            a_tf_inboard_total / n_tf_coils
        ) - a_tf_wp_with_insulation

        # Outboard leg cross-sectional area of surrounding case [m²]
        a_tf_coil_outboard_case = a_tf_leg_outboard - a_tf_wp_with_insulation

        # Front casing area [m²]
        if i_tf_case_geom == 0:
            # Circular front case
            a_tf_plasma_case = (
                rad_tf_coil_inboard_toroidal_half * r_tf_inboard_out**2
            ) - (tan_theta_coil * r_tf_wp_inboard_outer**2)
        else:
            # Straight front case [m²]
            a_tf_plasma_case = (
                (r_tf_wp_inboard_outer + dr_tf_plasma_case) ** 2
                - r_tf_wp_inboard_outer**2
            ) * tan_theta_coil

        # Nose casing area [m²]
        a_tf_coil_nose_case = (
            tan_theta_coil * r_tf_wp_inboard_inner**2
            - rad_tf_coil_inboard_toroidal_half * r_tf_inboard_in**2
        )

        # Report error if the casing area is negative
        if a_tf_coil_inboard_case <= 0.0e0 or a_tf_coil_outboard_case <= 0.0e0:
            error_handling.fdiags[0] = a_tf_coil_inboard_case
            error_handling.fdiags[1] = a_tf_coil_outboard_case
            error_handling.report_error(99)

        # Average lateral casing thickness [m]
        # --------------
        # Rectangular casing
        if i_tf_wp_geom == 0:
            dx_tf_side_case_average = (
                dx_tf_side_case_min + 0.5e0 * tan_theta_coil * dr_tf_wp_with_insulation
            )

        # Double rectangular WP
        elif i_tf_wp_geom == 1:
            dx_tf_side_case_average = (
                dx_tf_side_case_min + 0.25e0 * tan_theta_coil * dr_tf_wp_with_insulation
            )

        # Trapezoidal WP
        else:
            dx_tf_side_case_average = dx_tf_side_case_min

        return (
            a_tf_coil_inboard_case,
            a_tf_coil_outboard_case,
            a_tf_plasma_case,
            a_tf_coil_nose_case,
            dx_tf_side_case_average,
        )

    def tf_integer_turn_geom(
        self, n_layer, n_pancake, dx_tf_turn_steel, dx_tf_turn_insulation
    ):
        """
        Authors: J Morris & S Khan
        Setting the TF WP turn geometry for SC magnets from the number
        of turns rows in the radial direction. The turns can have any
        rectangular shapes.
        This calculation has two purposes, first to check if a turn can exist
        (positive cable space) and the second to provide its dimenesions,
        areas and the its associated current

        """
        sctfcoil_module.rbcndut = dx_tf_turn_steel * 0.75e0

        # Radial turn dimension [m]
        sctfcoil_module.dr_tf_turn = (
            tfcoil_variables.dr_tf_wp_with_insulation
            - 2.0e0
            * (
                tfcoil_variables.dx_tf_wp_insulation
                + tfcoil_variables.dx_tf_wp_insertion_gap
            )
        ) / n_layer

        if sctfcoil_module.dr_tf_turn <= (
            2.0e0 * dx_tf_turn_insulation + 2.0e0 * dx_tf_turn_steel
        ):
            error_handling.fdiags[0] = sctfcoil_module.dr_tf_turn
            error_handling.fdiags[1] = dx_tf_turn_insulation
            error_handling.fdiags[2] = dx_tf_turn_steel
            error_handling.report_error(100)

        # Toroidal turn dimension [m]
        sctfcoil_module.dx_tf_turn = (
            sctfcoil_module.dx_tf_wp_toroidal_min
            - 2.0e0
            * (
                tfcoil_variables.dx_tf_wp_insulation
                + tfcoil_variables.dx_tf_wp_insertion_gap
            )
        ) / n_pancake

        if sctfcoil_module.dx_tf_turn <= (
            2.0e0 * dx_tf_turn_insulation + 2.0e0 * dx_tf_turn_steel
        ):
            error_handling.fdiags[0] = sctfcoil_module.dx_tf_turn
            error_handling.fdiags[1] = dx_tf_turn_insulation
            error_handling.fdiags[2] = dx_tf_turn_steel
            error_handling.report_error(100)

        tfcoil_variables.t_turn_tf = np.sqrt(
            sctfcoil_module.dr_tf_turn * sctfcoil_module.dx_tf_turn
        )

        # Number of TF turns
        n_tf_coil_turns = np.double(n_layer * n_pancake)

        # Current per turn [A/turn]
        c_tf_turn = sctfcoil_module.c_tf_coil / n_tf_coil_turns

        # Radial and toroidal dimension of conductor [m]
        sctfcoil_module.t_conductor_radial = (
            sctfcoil_module.dr_tf_turn - 2.0e0 * dx_tf_turn_insulation
        )
        sctfcoil_module.t_conductor_toroidal = (
            sctfcoil_module.dx_tf_turn - 2.0e0 * dx_tf_turn_insulation
        )
        tfcoil_variables.t_conductor = np.sqrt(
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
        )

        # Dimension of square cable space inside conduit [m]
        sctfcoil_module.dr_tf_turn_cable_space = (
            sctfcoil_module.t_conductor_radial - 2.0e0 * dx_tf_turn_steel
        )
        sctfcoil_module.dx_tf_turn_cable_space = (
            sctfcoil_module.t_conductor_toroidal - 2.0e0 * dx_tf_turn_steel
        )
        sctfcoil_module.dx_tf_turn_cable_space_average = np.sqrt(
            sctfcoil_module.dr_tf_turn_cable_space
            * sctfcoil_module.dx_tf_turn_cable_space
        )

        # Cross-sectional area of cable space per turn
        # taking account of rounded inside corners [m2]
        a_tf_turn_cable_space_no_void = (
            sctfcoil_module.dr_tf_turn_cable_space
            * sctfcoil_module.dx_tf_turn_cable_space
        ) - (4.0e0 - np.pi) * sctfcoil_module.rbcndut**2

        if a_tf_turn_cable_space_no_void <= 0.0e0:
            if (sctfcoil_module.dr_tf_turn_cable_space < 0.0e0) or (
                sctfcoil_module.dx_tf_turn_cable_space < 0.0e0
            ):
                error_handling.fdiags[0] = a_tf_turn_cable_space_no_void
                error_handling.fdiags[1] = sctfcoil_module.dr_tf_turn_cable_space
                error_handling.fdiags[2] = sctfcoil_module.dx_tf_turn_cable_space
                error_handling.report_error(101)
            else:
                error_handling.fdiags[0] = a_tf_turn_cable_space_no_void
                error_handling.fdiags[1] = sctfcoil_module.dr_tf_turn_cable_space
                error_handling.fdiags[1] = sctfcoil_module.dx_tf_turn_cable_space
                error_handling.report_error(102)
                sctfcoil_module.rbcndut = 0.0e0
                a_tf_turn_cable_space_no_void = (
                    sctfcoil_module.dr_tf_turn_cable_space
                    * sctfcoil_module.dx_tf_turn_cable_space
                )

        # Cross-sectional area of conduit jacket per turn [m2]
        a_tf_turn_steel = (
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
            - a_tf_turn_cable_space_no_void
        )

        # Area of inter-turn insulation: single turn [m2]
        a_tf_turn_insulation = (
            sctfcoil_module.dr_tf_turn * sctfcoil_module.dx_tf_turn
            - a_tf_turn_steel
            - a_tf_turn_cable_space_no_void
        )
        return (
            a_tf_turn_cable_space_no_void,
            a_tf_turn_steel,
            a_tf_turn_insulation,
            c_tf_turn,
            n_tf_coil_turns,
        )

        # -------------

    def tf_wp_currents(self):
        """
        Author : S. Kahn, CCFE
        Turn engineering turn currents/densities
        """
        tfcoil_variables.j_tf_wp = max(
            1.0e0,
            tfcoil_variables.c_tf_total
            / (tfcoil_variables.n_tf_coils * sctfcoil_module.a_tf_wp_no_insulation),
        )

    def tf_averaged_turn_geom(
        self, j_tf_wp, dx_tf_turn_steel, dx_tf_turn_insulation, i_tf_sc_mat
    ):
        """
        subroutine straight from Python, see comments in tf_averaged_turn_geom_wrapper
        Authors : J. Morris, CCFE
        Authors : S. Kahn, CCFE
        Setting the TF WP turn geometry for SC magnets from the number
        the current per turn.
        This calculation has two purposes, first to check if a turn can exist
        (positive cable space) and the second to provide its dimensions,
        areas and the (float) number of turns
        """

        # Turn dimension is a an input
        if tfcoil_variables.t_turn_tf_is_input:
            # Turn area [m2]
            a_turn = tfcoil_variables.t_turn_tf**2

            # Current per turn [A]
            tfcoil_variables.c_tf_turn = a_turn * j_tf_wp

        # Turn cable dimension is an input
        elif tfcoil_variables.t_cable_tf_is_input:
            # Turn squared dimension [m]
            tfcoil_variables.t_turn_tf = tfcoil_variables.t_cable_tf + 2.0e0 * (
                dx_tf_turn_insulation + dx_tf_turn_steel
            )

            # Turn area [m2]
            a_turn = tfcoil_variables.t_turn_tf**2

            # Current per turn [A]
            tfcoil_variables.c_tf_turn = a_turn * j_tf_wp

        # Current per turn is an input
        else:
            # Turn area [m2]
            # Allow for additional inter-layer insulation MDK 13/11/18
            # Area of turn including conduit and inter-layer insulation
            a_turn = tfcoil_variables.c_tf_turn / j_tf_wp

            # Dimension of square cross-section of each turn including inter-turn insulation [m]
            tfcoil_variables.t_turn_tf = np.sqrt(a_turn)

        # Square turn assumption
        sctfcoil_module.dr_tf_turn = tfcoil_variables.t_turn_tf
        sctfcoil_module.dx_tf_turn = tfcoil_variables.t_turn_tf

        # See derivation in the following document
        # k:\power plant physics and technology\process\hts\hts coil module for process.docx
        tfcoil_variables.t_conductor = (
            -tfcoil_variables.layer_ins
            + np.sqrt(tfcoil_variables.layer_ins**2 + 4.0e00 * a_turn)
        ) / 2 - 2.0e0 * dx_tf_turn_insulation

        # Total number of turns per TF coil (not required to be an integer)
        n_tf_coil_turns = sctfcoil_module.a_tf_wp_no_insulation / a_turn

        # Area of inter-turn insulation: single turn [m2]
        a_tf_turn_insulation = a_turn - tfcoil_variables.t_conductor**2

        # NOTE: Fortran has a_tf_turn_cable_space_no_void as an intent(out) variable that was outputting
        # into tfcoil_variables.a_tf_turn_cable_space_no_void. The local variable, however, appears to
        # initially hold the value of tfcoil_variables.a_tf_turn_cable_space_no_void despite not being
        # intent(in). I have replicated this behaviour in Python for now.
        a_tf_turn_cable_space_no_void = copy.copy(
            tfcoil_variables.a_tf_turn_cable_space_no_void
        )

        # ITER like turn structure
        if i_tf_sc_mat != 6:
            # Radius of rounded corners of cable space inside conduit [m]
            rbcndut = dx_tf_turn_steel * 0.75e0

            # Dimension of square cable space inside conduit [m]
            sctfcoil_module.dx_tf_turn_cable_space_average = (
                tfcoil_variables.t_conductor - 2.0e0 * dx_tf_turn_steel
            )

            # Cross-sectional area of cable space per turn
            # taking account of rounded inside corners [m2]
            a_tf_turn_cable_space_no_void = (
                sctfcoil_module.dx_tf_turn_cable_space_average**2
                - (4.0e0 - np.pi) * rbcndut**2
            )

            if a_tf_turn_cable_space_no_void <= 0.0e0:
                if tfcoil_variables.t_conductor < 0.0e0:
                    error_handling.fdiags[0] = a_tf_turn_cable_space_no_void
                    error_handling.fdiags[1] = (
                        sctfcoil_module.dx_tf_turn_cable_space_average
                    )
                    error_handling.report_error(101)
                else:
                    error_handling.fdiags[0] = a_tf_turn_cable_space_no_void
                    error_handling.fdiags[1] = (
                        sctfcoil_module.dx_tf_turn_cable_space_average
                    )
                    error_handling.report_error(102)
                    rbcndut = 0.0e0
                    a_tf_turn_cable_space_no_void = (
                        sctfcoil_module.dx_tf_turn_cable_space_average**2
                    )

            # Cross-sectional area of conduit jacket per turn [m2]
            a_tf_turn_steel = (
                tfcoil_variables.t_conductor**2 - a_tf_turn_cable_space_no_void
            )

        # REBCO turn structure
        elif i_tf_sc_mat == 6:
            # Diameter of circular cable space inside conduit [m]
            sctfcoil_module.dx_tf_turn_cable_space_average = (
                tfcoil_variables.t_conductor - 2.0e0 * dx_tf_turn_steel
            )

            # Cross-sectional area of conduit jacket per turn [m2]
            a_tf_turn_steel = (
                tfcoil_variables.t_conductor**2 - a_tf_turn_cable_space_no_void
            )

        return (
            a_tf_turn_cable_space_no_void,
            a_tf_turn_steel,
            a_tf_turn_insulation,
            n_tf_coil_turns,
        )


@staticmethod
def lambda_term(tau: float, omega: float) -> float:
    """
    The lambda function used inegral in inductance calcuation found
    in Y. Itoh et al. The full form of the functions are given in appendix A.

    Author: Alexander Pearce, UKAEA

    :param tau: tau_{s,k} = (R_{s,k} - R_{c,k}) / R_k
    :param omega: omega_k = R_{c,k}/R_k
    :returns: integral
    """
    p = 1.0 - omega**2.0

    if p < 0:
        integral = (1.0 / np.sqrt(np.abs(p))) * np.arcsin(
            (1.0 + omega * tau) / (tau + omega)
        )
    else:
        integral = (1.0 / np.sqrt(np.abs(p))) * np.log(
            (2.0 * (1.0 + tau * omega - np.sqrt(p * (1 - tau**2.0)))) / (tau + omega)
        )

    return integral


@staticmethod
def _theta_factor_integral(
    ro_vv: float, ri_vv: float, rm_vv: float, h_vv: float, theta1_vv: float
) -> float:
    """
    The calcuation of the theta factor found in Eq 4 of Y. Itoh et al. The
    full form of the integral is given in appendix A.

    Author: Alexander Pearce, UKAEA

    :param ro_vv: the radius of the outboard edge of the VV CCL
    :param ri_vv: the radius of the inboard edge of the VV CCL
    :param rm_vv: the radius where the maximum height of the VV CCL occurs
    :param h_vv: the maximum height of the VV CCL
    :param theta1_vv: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the VV CCL
    :returns: theta factor.
    """
    theta2 = np.pi / 2.0 + theta1_vv
    a = (ro_vv - ri_vv) / 2.0
    rbar = (ro_vv + ri_vv) / 2.0
    delta = (rbar - rm_vv) / a
    kappa = h_vv / a
    iota = (1.0 + delta) / kappa

    denom = np.cos(theta1_vv) + np.sin(theta1_vv) - 1.0

    r1 = h_vv * ((np.cos(theta1_vv) + iota * (np.sin(theta1_vv) - 1.0)) / denom)
    r2 = h_vv * ((np.cos(theta1_vv) - 1.0 + iota * np.sin(theta1_vv)) / denom)
    r3 = h_vv * (1 - delta) / kappa

    rc1 = (h_vv / kappa) * ((rbar / a) + 1.0) - r1
    rc2 = rc1 + (r1 - r2) * np.cos(theta1_vv)
    rc3 = rc2
    zc2 = (r1 - r2) * np.sin(theta1_vv)
    zc3 = zc2 + r2 - r3

    tau = np.array([
        [np.cos(theta1_vv), np.cos(theta1_vv + theta2), -1.0],
        [1.0, np.cos(theta1_vv), np.cos(theta1_vv + theta2)],
    ])

    omega = np.array([rc1 / r1, rc2 / r2, rc3 / r3])

    # Assume up down symmetry and let Zc6 = - Zc3
    chi1 = (zc3 + np.abs(-zc3)) / ri_vv
    chi2 = 0.0

    for k in range(len(omega)):
        chi2 = chi2 + np.abs(
            lambda_term(tau[1, k], omega[k]) - lambda_term(tau[0, k], omega[k])
        )

    return (chi1 + 2.0 * chi2) / (2.0 * np.pi)


def vv_stress_on_quench(
    # TF shape
    H_coil: float,
    ri_coil: float,
    ro_coil: float,
    rm_coil: float,
    ccl_length_coil: float,
    theta1_coil: float,
    # VV shape
    H_vv: float,
    ri_vv: float,
    ro_vv: float,
    rm_vv: float,
    theta1_vv: float,
    # TF properties
    n_tf_coils: float,
    n_tf_coil_turns: float,
    s_rp: float,
    s_cc: float,
    taud: float,
    i_op: float,
    # VV properties
    d_vv: float,
) -> float:
    """Generic model to calculate the Tresca stress of the
    Vacuum Vessel (VV), experienced when the TF coil quenches.

    Author: Timothy Nunn, UKAEA

    The current center line (CCL) of a structure is an appoximation
    of where the poloidal current of a structure acts. This model
    considers how the current (self-)induced in the following structures
    affects the Tresca stress on the VV (corresponding to the structure
    index used in [1]):

    0. the TF coil conductors
    1. the TF coil steel structure
    2. the VV

    :param H_coil: the maximum height of the TF coil CCL
    :param ri_coil: the radius of the inboard edge of the TF coil CCL
    :param ro_coil: the radius of the outboard edge of the TF coil CCL
    :param rm_coil: the radius where the maximum height of the TF coil CCL occurs
    :param ccl_length_coil: the length of the TF coil CCL
    :param theta1_coil: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the coil CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).

    :param H_vv: the maximum height of the VV CCL
    :param ri_vv: the radius of the inboard edge of the VV CCL
    :param ro_vv: the radius of the outboard edge of the VV CCL
    :param rm_vv: the radius where the maximum height of the VV CCL occurs
    :param theta1_vv: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the VV CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).

    :param n_tf_coils: the number of TF coils
    :param n_tf_coil_turns: the number of turns per TF coil
    :param s_rp: the cross-sectional area of the radial plates of the TF coil
    :param s_cc: the cross-sectional area of the TF coil case
    :param taud: the discharge time of the TF coil when quench occurs
    :param i_op: the 'normal' operating current of the TF coil

    :param d_vv: the thickness of the vacuum vessel shell

    :returns: the maximum stress experienced by the vacuum vessel

    Notes
    -----
    The theta1 quantity for the TF coil and VV is not very meaningful. The
    impact of it of the inductance is rather small. Generally, the paper seems to
    suggest the TF coil is between 40 and 60, as this is the range they calculate
    the surrogates over. The thickness of the VV considers an ITER like design and
    only the outer and inner shells that which act of conductuve structural material.

    References
    ----------
    1. ITOH, Yasuyuki & Utoh, Hiroyasu & SAKAMOTO, Yoshiteru & Hiwatari, Ryoji. (2020).
    Empirical Formulas for Estimating Self and Mutual Inductances of Toroidal Field Coils and Structures.
    Plasma and Fusion Research. 15. 1405078-1405078. 10.1585/pfr.15.1405078.
    """
    # Convert angles into radians
    theta1_vv_rad = np.pi * (theta1_vv / 180.0)

    # Poloidal loop resistance (PLR) in ohms
    theta_vv = _theta_factor_integral(ro_vv, ri_vv, rm_vv, H_vv, theta1_vv_rad)
    plr_coil = ((0.5 * ccl_length_coil) / (n_tf_coils * (s_cc + s_rp))) * 1e-6
    plr_vv = ((0.84 / d_vv) * theta_vv) * 1e-6

    # relevant self-inductances in henry (H)
    coil_structure_self_inductance = (
        (constants.rmu0 / np.pi)
        * H_coil
        * _inductance_factor(H_coil, ri_coil, ro_coil, rm_coil, theta1_coil)
    )
    vv_self_inductance = (
        (constants.rmu0 / np.pi)
        * H_vv
        * _inductance_factor(H_vv, ri_vv, ro_vv, rm_vv, theta1_vv)
    )

    # s^-1
    lambda0 = 1 / taud
    lambda1 = (plr_coil) / coil_structure_self_inductance
    lambda2 = (plr_vv) / vv_self_inductance

    # approximate time at which the maximum force (and stress) will occur on the VV
    tmaxforce = np.log((lambda0 + lambda1) / (2 * lambda0)) / (lambda1 - lambda0)

    i0 = i_op * np.exp(-lambda0 * tmaxforce)
    i1 = (
        lambda0
        * n_tf_coils
        * n_tf_coil_turns
        * i_op
        * (
            (np.exp(-lambda1 * tmaxforce) - np.exp(-lambda0 * tmaxforce))
            / (lambda0 - lambda1)
        )
    )
    i2 = (lambda1 / lambda2) * i1

    a_vv = (ro_vv + ri_vv) / (ro_vv - ri_vv)
    b_vvi = (constants.rmu0 * (n_tf_coils * n_tf_coil_turns * i0 + i1 + (i2 / 2))) / (
        2 * np.pi * ri_vv
    )
    j_vvi = i2 / (2 * np.pi * d_vv * ri_vv)

    zeta = 1 + ((a_vv - 1) * np.log((a_vv + 1) / (a_vv - 1)) / (2 * a_vv))

    return zeta * b_vvi * j_vvi * ri_vv


def _inductance_factor(
    H: float, ri: float, ro: float, rm: float, theta1: float
) -> float:
    """Calculates the inductance factor for a toroidal structure
    using surrogate 2 #1866.

    Author: Timothy Nunn, UKAEA

    :param H: the maximum height of the structure
    :param ri: the radius of the inboard side of the structure CCL
    :param ro: the radius of the outboard side of the structure CCL
    :param rm: the radius at which `H` occurs
    :param theta1: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the structure CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).
    """
    # NOTE: the following properties are not those of the plasma but of
    # the VV/coil structures
    major_radius = (ro + ri) / 2
    minor_radius = (ro - ri) / 2

    aspect_ratio = major_radius / minor_radius
    triangularity = (major_radius - rm) / minor_radius
    elongation = H / minor_radius

    return (
        4.933
        + 0.03728 * elongation
        + 0.06980 * triangularity
        - 3.551 * aspect_ratio
        + 0.7629 * aspect_ratio**2
        - 0.06298 * (theta1 / 90)
    )


def init_sctfcoil_module():
    sctfcoil_module.is_leg_cp_temp_same = 0
    sctfcoil_module.tf_fit_t = 0.0
    sctfcoil_module.tf_fit_z = 0.0
    sctfcoil_module.tf_fit_y = 0.0
    sctfcoil_module.c_tf_coil = 0.0
    sctfcoil_module.a_tf_wp_with_insulation = 0.0
    sctfcoil_module.a_tf_wp_no_insulation = 0.0
    sctfcoil_module.a_tf_coil_inboard_steel = 0.0
    sctfcoil_module.a_tf_coil_inboard_insulation = 0.0
    sctfcoil_module.f_a_tf_coil_inboard_steel = 0.0
    sctfcoil_module.f_a_tf_coil_inboard_insulation = 0.0
    sctfcoil_module.z_cp_top = 0.0
    sctfcoil_module.r_tf_outboard_in = 0.0
    sctfcoil_module.r_tf_outboard_out = 0.0
    sctfcoil_module.r_tf_wp_inboard_inner = 0.0
    sctfcoil_module.r_tf_wp_inboard_outer = 0.0
    sctfcoil_module.r_tf_wp_inboard_centre = 0.0
    sctfcoil_module.dr_tf_wp_top = 0.0
    sctfcoil_module.vol_ins_cp = 0.0
    sctfcoil_module.vol_gr_ins_cp = 0.0
    sctfcoil_module.vol_case_cp = 0.0
    sctfcoil_module.dx_tf_wp_toroidal_min = 0.0
    sctfcoil_module.dx_tf_wp_toroidal_average = 0.0
    sctfcoil_module.dx_tf_side_case_average = 0.0
    sctfcoil_module.a_tf_plasma_case = 0.0
    sctfcoil_module.a_tf_coil_nose_case = 0.0
    sctfcoil_module.a_tf_wp_ground_insulation = 0.0
    sctfcoil_module.a_leg_ins = 0.0
    sctfcoil_module.a_leg_gr_ins = 0.0
    sctfcoil_module.a_leg_cond = 0.0
    sctfcoil_module.rad_tf_coil_inboard_toroidal_half = 0.0
    sctfcoil_module.tan_theta_coil = 0.0
    sctfcoil_module.t_conductor_radial = 0.0
    sctfcoil_module.t_conductor_toroidal = 0.0
    sctfcoil_module.dr_tf_turn_cable_space = 0.0
    sctfcoil_module.dx_tf_turn_cable_space = 0.0
    sctfcoil_module.dr_tf_turn = 0.0
    sctfcoil_module.dx_tf_turn = 0.0
    sctfcoil_module.dx_tf_turn_cable_space_average = 0.0
    sctfcoil_module.vforce_inboard_tot = 0.0
    sctfcoil_module.t1 = 0.0
    sctfcoil_module.time2 = 0.0
    sctfcoil_module.tau2 = 0.0
    sctfcoil_module.e_tf_magnetic_stored_total = 0.0
