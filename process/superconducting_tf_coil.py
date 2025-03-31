import copy
import logging

import numpy as np
from scipy import optimize

import process.superconductors as superconductors
from process import process_output as po
from process.fortran import (
    build_variables,
    constants,
    divertor_variables,
    error_handling,
    global_variables,
    numerics,
    pfcoil_variables,
    physics_variables,
    rebco_variables,
    sctfcoil_module,
    tfcoil_variables,
)
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
        self.tf_global_geometry()
        self.tf_current()
        self.tf_coil_shape_inner()
        self.sc_tf_internal_geom(
            tfcoil_variables.i_tf_wp_geom,
            tfcoil_variables.i_tf_case_geom,
            tfcoil_variables.i_tf_turns_integer,
        )

        if physics_variables.itart == 0 and tfcoil_variables.i_tf_shape == 1:
            tfcoil_variables.ind_tf_coil = self.tfcind(
                build_variables.dr_tf_inboard,
                tfcoil_variables.xarc,
                tfcoil_variables.yarc,
            )
        else:
            tfcoil_variables.ind_tf_coil = (
                (build_variables.hmax + build_variables.dr_tf_outboard)
                * RMU0
                / constants.pi
                * np.log(
                    build_variables.r_tf_outboard_mid / build_variables.r_tf_inboard_mid
                )
            )

        # Total TF coil stored magnetic energy [J]
        sctfcoil_module.estotft = (
            0.5e0 * tfcoil_variables.ind_tf_coil * tfcoil_variables.c_tf_total**2
        )

        # Total TF coil stored magnetic energy [Gigajoule]
        tfcoil_variables.estotftgj = 1.0e-9 * sctfcoil_module.estotft

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
                build_variables.hmax,
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
                tfcoil_variables.thicndut,
                tfcoil_variables.eyoung_copper,
                tfcoil_variables.poisson_copper,
                tfcoil_variables.i_tf_sup,
                tfcoil_variables.eyoung_res_tf_buck,
                sctfcoil_module.r_wp_inner,
                sctfcoil_module.tan_theta_coil,
                sctfcoil_module.rad_tf_coil_toroidal,
                sctfcoil_module.r_wp_outer,
                sctfcoil_module.a_tf_steel,
                sctfcoil_module.a_case_front,
                sctfcoil_module.a_case_nose,
                tfcoil_variables.tfinsgap,
                tfcoil_variables.tinstf,
                tfcoil_variables.n_tf_turn,
                int(tfcoil_variables.i_tf_turns_integer),
                sctfcoil_module.t_cable,
                sctfcoil_module.t_cable_radial,
                tfcoil_variables.dhecoil,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.thwcndut,
                sctfcoil_module.t_lat_case_av,
                sctfcoil_module.t_wp_toroidal_av,
                sctfcoil_module.a_tf_ins,
                tfcoil_variables.aswp,
                tfcoil_variables.acond,
                sctfcoil_module.awpc,
                tfcoil_variables.eyoung_al,
                tfcoil_variables.poisson_al,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf_graded_layers,
                tfcoil_variables.c_tf_total,
                tfcoil_variables.casthi,
                tfcoil_variables.i_tf_stress_model,
                sctfcoil_module.vforce_inboard_tot,
                tfcoil_variables.i_tf_tresca,
                tfcoil_variables.acasetf,
                tfcoil_variables.vforce,
                tfcoil_variables.acndttf,
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
            tfcoil_variables.wwp1,
            tfcoil_variables.dr_tf_wp
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap),
            sctfcoil_module.r_wp_centre,
            tfcoil_variables.b_tf_inboard_peak,
        )

        tfes = sctfcoil_module.estotft / tfcoil_variables.n_tf_coils
        # Cross-sectional area per turn
        aturn = tfcoil_variables.c_tf_total / (
            tfcoil_variables.j_tf_wp
            * tfcoil_variables.n_tf_coils
            * tfcoil_variables.n_tf_turn
        )

        if tfcoil_variables.i_tf_sc_mat == 6:
            (tfcoil_variables.jwdgcrt, tfcoil_variables.tmargtf) = self.supercon_croco(
                aturn,
                tfcoil_variables.bmaxtfrp,
                tfcoil_variables.cpttf,
                tfcoil_variables.tftmp,
                output=output,
            )

            tfcoil_variables.vtfskv = (
                self.croco_voltage() / 1.0e3
            )  # TFC Quench voltage in kV

        else:
            (
                tfcoil_variables.jwdgcrt,
                vdump,
                tfcoil_variables.tmargtf,
            ) = self.supercon(
                tfcoil_variables.acstf,
                aturn,
                tfcoil_variables.bmaxtfrp,
                tfcoil_variables.vftf,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.cpttf,
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
                * (sctfcoil_module.estotft / tfcoil_variables.n_tf_coils)
                / tfcoil_variables.cpttf
            )
        elif f2py_compatible_to_string(tfcoil_variables.quench_model) == "exponential":
            sctfcoil_module.tau2 = tfcoil_variables.tdmptf
            croco_voltage = (
                2.0e0
                / sctfcoil_module.tau2
                * (sctfcoil_module.estotft / tfcoil_variables.n_tf_coils)
                / tfcoil_variables.cpttf
            )
        else:
            return 0.0

        return croco_voltage

    def supercon_croco(self, aturn, bmax, iop, thelium, output: bool):
        """TF superconducting CroCo conductor using REBCO tape
        author: M Kovari, CCFE, Culham Science Centre
        bmax : input real : Peak field at conductor (T)
        iop : input real : Operating current per turn (A)
        thelium : input real : He temperature at peak field point (K)
        iprint : input integer : Switch for printing (1 = yes, 0 = no)
        outfile : input integer : Fortran output unit identifier
        jwdgcrt : output real : Critical winding pack current density (A/m2)
        tmarg : output real : Temperature margin (K)
        """

        j_crit_sc: float = 0.0
        #  Find critical current density in superconducting cable, j_crit_cable
        j_crit_sc, _ = superconductors.jcrit_rebco(thelium, bmax)
        # tfcoil_variables.acstf : Cable space - inside area (m2)
        # Set new rebco_variables.croco_od
        # allowing for scaling of rebco_variables.croco_od
        rebco_variables.croco_od = (
            tfcoil_variables.t_conductor / 3.0e0
            - tfcoil_variables.thwcndut * (2.0e0 / 3.0e0)
        )
        sctfcoil_module.conductor_acs = (
            9.0e0 / 4.0e0 * np.pi * rebco_variables.croco_od**2
        )
        tfcoil_variables.acstf = sctfcoil_module.conductor_acs
        sctfcoil_module.conductor_area = (
            tfcoil_variables.t_conductor**2
        )  # does this not assume it's a sqaure???

        sctfcoil_module.conductor_jacket_area = (
            sctfcoil_module.conductor_area - sctfcoil_module.conductor_acs
        )
        tfcoil_variables.acndttf = sctfcoil_module.conductor_jacket_area

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
        # aturn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        jwdgcrt = icrit / aturn
        #  Ratio of operating / critical current
        iooic = iop / icrit
        #  Operating current density
        jwdgop = iop / aturn
        #  Actual current density in superconductor,
        # which should be equal to jcrit(thelium+tmarg)

        #  when we have found the desired value of tmarg
        jsc = iooic * j_crit_sc

        # Temperature margin
        current_sharing_t = superconductors.current_sharing_rebco(bmax, jsc)
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
                bmax: {bmax}"""
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
                "(acstf)",
                tfcoil_variables.acstf,
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
                "(jwdgcrt)",
                jwdgcrt,
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
                "(cpttf)",
                tfcoil_variables.cpttf,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return jwdgcrt, tmarg

    def supercon(
        self,
        acs,
        aturn,
        bmax,
        fhe,
        fcu,
        iop,
        jwp,
        isumat,
        fhts,
        tdmptf,
        tfes,
        thelium,
        tmax,
        bcritsc,
        tcritsc,
        output: bool,
    ):
        """Routine to calculate the TF superconducting conductor  properties
        author: P J Knight, CCFE, Culham Science Centre
        author: J Galambos, ORNL
        author: R Kemp, CCFE, Culham Science Centre
        author: M Kovari, CCFE, Culham Science Centre
        author: J Miller, ORNL
        acs : input real : Cable space - inside area (m2)
        aturn : input real : Area per turn (i.e. entire jacketed conductor) (m2)
        bmax : input real : Peak field at conductor (T)
        fhe : input real : Fraction of cable space that is for He cooling
        fcu : input real : Fraction of conductor that is copper
        iop : input real : Operating current per turn (A)
        jwp : input real : Actual winding pack current density (A/m2)
        isumat : input integer : Switch for conductor type:
        1 = ITER Nb3Sn, standard parameters,
        2 = Bi-2212 High Temperature Superconductor,
        3 = NbTi,
        4 = ITER Nb3Sn, user-defined parameters
        5 = WST Nb3Sn parameterisation
        7 = Durham Ginzburg-Landau Nb-Ti parameterisation
        8 = Durham Ginzburg-Landau critical surface model for REBCO
        9 = Hazelton experimental data + Zhai conceptual model for REBCO
        fhts    : input real : Adjustment factor (<= 1) to account for strain,
        radiation damage, fatigue or AC losses
        tdmptf : input real : Dump time (sec)
        tfes : input real : Energy stored in one TF coil (J)
        thelium : input real : He temperature at peak field point (K)
        tmax : input real : Max conductor temperature during quench (K)
        bcritsc : input real : Critical field at zero temperature and strain (T) (isumat=4 only)
        tcritsc : input real : Critical temperature at zero field and strain (K) (isumat=4 only)
        iprint : input integer : Switch for printing (1 = yes, 0 = no)
        outfile : input integer : Fortran output unit identifier
        jwdgpro : output real : Winding pack current density from temperature
        rise protection (A/m2)
        jwdgcrt : output real : Critical winding pack current density (A/m2)
        vd : output real : Discharge voltage imposed on a TF coil (V)
        tmarg : output real : Temperature margin (K)
        This routine calculates the superconductor properties for the TF coils.
        It was originally programmed by J. Galambos 1991, from algorithms provided
        by J. Miller.
        <P>The routine calculates the critical current density (winding pack)
        and also the protection information (for a quench).
        NOT used for the Croco conductor

        N.B. critical current density for a super conductor (j_crit_sc)
        is for the superconducting strands/tape, not including copper.
        Critical current density for a cable (j_crit_cable) acounts for
        both the fraction of the cable taken up by helium coolant channels,
        and the cable conductor copper fraction - i.e., the copper in the
        superconducting strands AND any addtional copper, such as REBCO
        tape support.
        """
        tdump = tdmptf

        # Helium channel
        fhetot = (
            fhe
            + (np.pi / 4.0e0)
            * tfcoil_variables.dhecoil
            * tfcoil_variables.dhecoil
            / acs
        )

        # Guard against negative conductor fraction fcond
        # Kludge to allow solver to continue and hopefully be constrained away
        # from this point
        if fhetot > 0.99:
            fhetot = 0.99

        #  Conductor fraction (including central helium channel)
        fcond = 1.0e0 - fhetot

        if tfcoil_variables.i_str_wp == 0:
            strain = tfcoil_variables.str_tf_con_res
        else:
            strain = tfcoil_variables.str_wp

        # Find critical current density in the superconducter (j_crit_sc)
        # and the superconducting cable (j_crit_cable)
        if isumat == 1:  # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            #  j_crit_sc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            j_crit_sc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable
            icrit = j_crit_cable * acs

            # Strand critical current calculation for costing in $/kAm
            # = Superconducting filaments jc * (1 - strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        elif isumat == 2:  # Bi-2212 high temperature superconductor parameterization
            #  Current density in a strand of Bi-2212 conductor
            #  N.B. jcrit returned by superconductors.bi2212 is the critical current density
            #  in the strand, not just the superconducting portion.
            #  The parameterization for j_crit_cable assumes a particular strand
            #  composition that does not require a user-defined copper fraction,
            #  so this is irrelevant in this model
            jstrand = jwp * aturn / (acs * fcond)

            j_crit_cable, tmarg = superconductors.bi2212(bmax, jstrand, thelium, fhts)
            j_crit_sc = j_crit_cable / (1.0e0 - fcu)
            #  Critical current in cable
            icrit = j_crit_cable * acs * fcond

            # Strand critical current calulation for costing in $ / kAm
            # Copper in the strand is already accounted for
            tfcoil_variables.j_crit_str_tf = j_crit_sc

        elif isumat == 3:  # NbTi data
            bc20m = 15.0e0
            tc0m = 9.3e0
            c0 = 1.0e10
            j_crit_sc, _ = superconductors.jcrit_nbti(thelium, bmax, c0, bc20m, tc0m)
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        elif isumat == 4:  # ITER Nb3Sn parameterization, but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            j_crit_sc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        elif isumat == 5:  # WST Nb3Sn parameterisation
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.5e-2

            #  j_crit_sc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            j_crit_sc, _, _ = superconductors.wstsc(thelium, bmax, strain, bc20m, tc0m)
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        elif isumat == 6:  # "REBCO" 2nd generation HTS superconductor in CrCo strand
            raise ValueError(
                "sctfcoil.supercon has been called but tfcoil_variables.i_tf_sc_mat=6"
            )

        elif isumat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
            bc20m = tfcoil_variables.b_crit_upper_nbti
            tc0m = tfcoil_variables.t_crit_nbti
            j_crit_sc, _, _ = superconductors.gl_nbti(
                thelium, bmax, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        elif isumat == 8:  # Durham Ginzburg-Landau critical surface model for REBCO
            bc20m = 430
            tc0m = 185
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.7e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = np.sign(strain) * 0.7e-2

            j_crit_sc, _, _ = superconductors.gl_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # Already includes buffer and support layers so no need to include fcu here
            tfcoil_variables.j_crit_str_tf = j_crit_sc

        elif (
            isumat == 9
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
                thelium, bmax, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * fcond
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = j_crit_cable * acs

            # Strand critical current calulation for costing in $ / kAm
            # = superconducting filaments jc * (1 -strand copper fraction)
            tfcoil_variables.j_crit_str_tf = j_crit_sc * (1.0e0 - fcu)

        else:
            error_handling.idiags[0] = isumat
            error_handling.report_error(105)

        # Critical current density in winding pack
        # aturn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        jwdgcrt = icrit / aturn
        #  Ratio of operating / critical current
        iooic = iop / icrit
        #  Operating current density
        jwdgop = iop / aturn
        #  Actual current density in superconductor, which should be equal to jcrit(thelium+tmarg)
        #  when we have found the desired value of tmarg
        jsc = iooic * j_crit_sc

        if iooic <= 0e0:
            logger.warning(
                f"""Negative Iop/Icrit for TF coil
            jsc: {jsc}
            iooic: {iooic}
            j_crit_sc: {j_crit_sc}
            Check conductor dimensions. Cable space area acs likely gone negative. acs: {acs}
            This is likely because t_cable_radial or t_cable_toroidal has gone negative:
            t_cable_radial: {sctfcoil_module.t_cable_radial}
            t_cable_toroidal: {sctfcoil_module.t_cable_toroidal}
            """
            )

        # REBCO measurements from 2 T to 14 T, extrapolating outside this
        if (isumat == 8) and (tfcoil_variables.bmaxtfrp >= 14):
            error_handling.report_error(266)

        #  Temperature margin (already calculated in superconductors.bi2212 for isumat=2)

        if (
            (isumat == 1)
            or (isumat == 3)
            or (isumat == 4)
            or (isumat == 5)
            or (isumat == 7)
            or (isumat == 8)
            or (isumat == 9)
        ):  # Find temperature at which current density margin = 0
            if isumat == 3:
                arguments = (isumat, jsc, bmax, strain, bc20m, tc0m, c0)
            else:
                arguments = (isumat, jsc, bmax, strain, bc20m, tc0m)

            another_estimate = 2 * thelium
            t_zero_margin, root_result = optimize.newton(
                superconductors.current_density_margin,
                thelium,
                fprime=None,
                args=arguments,
                # args=(isumat, jsc, bmax, strain, bc20m, tc0m,),
                tol=1.0e-06,
                maxiter=50,
                fprime2=None,
                x1=another_estimate,
                rtol=1.0e-6,
                full_output=True,
                disp=True,
            )
            # print(root_result)  # Diagnostic for newton method
            tmarg = t_zero_margin - thelium
            tfcoil_variables.temp_margin = tmarg

        #  Find the current density limited by the protection limit
        #  (N.B. Unclear of this routine's relevance for Bi-2212 (isumat=2), due
        #  to presence of fcu argument, which is not used for this model above)

        tfcoil_variables.jwdgpro, vd = self.protect(
            iop, tfes, acs, aturn, tdump, fcond, fcu, thelium, tmax
        )

        if output:  # Output --------------------------
            if tmarg <= 0.0e0:
                logger.warning(
                    """Negative TFC temperature margin
                tmarg: {tmarg}
                bmax: {bmax}
                jcrit0: {jcrit0}
                jsc: {jsc}
                """
                )

            po.oheadr(self.outfile, "Superconducting TF Coils")
            po.ovarin(self.outfile, "Superconductor switch", "(isumat)", isumat)

            if isumat == 1:
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
            if isumat == 2:
                po.ocmmnt(self.outfile, "Superconductor used: Bi-2212 HTS")
            if isumat == 3:
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
            if isumat == 4:
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
            if isumat == 5:
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
            if isumat == 7:
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
            if isumat == 8:
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
            if isumat == 9:
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
                "(thelium)",
                thelium,
            )
            po.ovarre(
                self.outfile,
                "Total helium fraction inside cable space",
                "(fhetot)",
                fhetot,
                "OP ",
            )
            po.ovarre(self.outfile, "Copper fraction of conductor", "(fcutfsu)", fcu)
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
                "(jwdgcrt)",
                jwdgcrt,
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
            po.ovarre(self.outfile, "Critical current (A)", "(icrit)", icrit, "OP ")
            po.ovarre(
                self.outfile,
                "Actual current (A)",
                "(cpttf)",
                tfcoil_variables.cpttf,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return jwdgcrt, vd, tmarg

    def protect(self, aio, tfes, acs, aturn, tdump, fcond, fcu, tba, tmax):
        """Finds the current density limited by the protection limit
        author: P J Knight, CCFE, Culham Science Centre
        author: J Miller, ORNL
        aio : input real : Operating current (A)
        tfes : input real : Energy stored in one TF coil (J)
        acs : input real : Cable space - inside area (m2)
        aturn : input real : Area per turn (i.e.  entire cable) (m2)
        tdump : input real : Dump time (sec)
        fcond : input real : Fraction of cable space containing conductor
        fcu : input real : Fraction of conductor that is copper
        tba : input real : He temperature at peak field point (K)
        tmax : input real : Max conductor temperature during quench (K)
        ajwpro : output real :  Winding pack current density from temperature
        rise protection (A/m2)
        vd : output real :  Discharge voltage imposed on a TF coil (V)
        This routine calculates maximum conductor current density which
        limits the peak temperature in the winding to a given limit (tmax).
        It also finds the dump voltage.
        <P>These calculations are based on Miller's formulations.
        """
        # Integration coefficients p1,p2,p3
        p1 = (
            0.0e0,
            0.8e0,
            1.75e0,
            2.4e0,
            2.7e0,
            2.95e0,
            3.1e0,
            3.2e0,
            3.3e0,
            3.4e0,
            3.5e0,
        )
        p2 = (
            0.0e0,
            0.05e0,
            0.5e0,
            1.4e0,
            2.6e0,
            3.7e0,
            4.6e0,
            5.3e0,
            5.95e0,
            6.55e0,
            7.1e0,
        )
        p3 = (
            0.0e0,
            0.05e0,
            0.5e0,
            1.4e0,
            2.6e0,
            3.7e0,
            4.6e0,
            5.4e0,
            6.05e0,
            6.8e0,
            7.2e0,
        )

        #  Dump voltage

        vd = 2.0e0 * tfes / (tdump * aio)

        #  Current density limited by temperature rise during quench

        tav = 1.0e0 + (tmax - tba) / 20.0e0
        n_o = int(tav)
        n_p = n_o + 1
        n_p = min(n_p, 11)

        ai1 = 1.0e16 * (p1[n_o - 1] + (p1[n_p - 1] - p1[n_o - 1]) * (tav - n_o))
        ai2 = 1.0e16 * (p2[n_o - 1] + (p2[n_p - 1] - p2[n_o - 1]) * (tav - n_o))
        ai3 = 1.0e16 * (p3[n_o - 1] + (p3[n_p - 1] - p3[n_o - 1]) * (tav - n_o))

        aa = vd * aio / tfes
        bb = (1.0e0 - fcond) * fcond * fcu * ai1
        cc = (fcu * fcond) ** 2 * ai2
        dd = (1.0e0 - fcu) * fcu * fcond**2 * ai3
        ajcp = np.sqrt(aa * (bb + cc + dd))
        ajwpro = ajcp * (acs / aturn)

        return ajwpro, vd

        # ---

    def vv_stress_on_quench(self):
        """Calculate the Tresca stress [Pa] of the Vacuum Vessel (VV)
        experienced when the TF coil quenches.

        Author: Timothy Nunn, UKAEA

        Assumes the current center line (CCL) of the TF coil is the
        middle of the coil.

        We assume vertical symmetry which is only true for double null
        machines.
        """
        H_coil = build_variables.hmax + (build_variables.dr_tf_inboard / 2)
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
            n_tf_turn=tfcoil_variables.n_tf_turn,
            # Area of the radial plate taken to be the area of steel in the WP
            # TODO: value clipped due to #1883
            s_rp=np.clip(sctfcoil_module.a_tf_steel, 0, None),
            s_cc=sctfcoil_module.a_case_front
            + sctfcoil_module.a_case_nose
            + 2.0 * sctfcoil_module.t_lat_case_av,
            taud=tfcoil_variables.tdmptf,
            # TODO: is this the correct current?
            i_op=sctfcoil_module.c_tf_coil / tfcoil_variables.n_tf_turn,
            # VV properties
            d_vv=build_variables.dr_vv_shells,
        )

    def peak_tf_with_ripple(self, n_tf_coils, wwp1, dr_tf_wp, tfin, b_tf_inboard_peak):
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
        :param wwp1: width of plasma-facing face of winding pack (m)
        :type wwp1: float
        :param dr_tf_wp: radial thickness of winding pack (m)
        :type dr_tf_wp: float
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

        wmax = (2.0e0 * tfin + dr_tf_wp) * np.tan(np.pi / n_tf_coils)

        #  Dimensionless winding pack width

        sctfcoil_module.tf_fit_t = wwp1 / wmax
        if (sctfcoil_module.tf_fit_t < 0.3e0) or (sctfcoil_module.tf_fit_t > 1.1e0):
            # write(*,*) 'PEAK_TF_WITH_RIPPLE: fitting problem; t = ',t
            flag = 1

        #  Dimensionless winding pack radial thickness

        sctfcoil_module.tf_fit_z = dr_tf_wp / wmax
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
        self.tf_wp_geom(i_tf_wp_geom)

        # Calculating the TF steel casing areas
        self.tf_case_geom(i_tf_wp_geom, i_tf_case_geom)

        # WP/trun currents
        self.tf_wp_currents()

        # sctfcoil_module.sc_tf_internal_geom(
        #     i_tf_wp_geom, i_tf_case_geom, i_tf_turns_integer
        # )

        # Setting the WP turn geometry / areas
        if i_tf_turns_integer == 0:
            # Non-ingeger number of turns
            (
                tfcoil_variables.acstf,
                tfcoil_variables.acndttf,
                tfcoil_variables.insulation_area,
                tfcoil_variables.n_tf_turn,
            ) = self.tf_averaged_turn_geom(
                tfcoil_variables.j_tf_wp,
                tfcoil_variables.thwcndut,
                tfcoil_variables.thicndut,
                tfcoil_variables.i_tf_sc_mat,
            )

        else:
            # Integer number of turns
            (
                tfcoil_variables.acstf,
                tfcoil_variables.acndttf,
                tfcoil_variables.insulation_area,
                tfcoil_variables.cpttf,
                tfcoil_variables.n_tf_turn,
            ) = self.tf_integer_turn_geom(
                tfcoil_variables.n_layer,
                tfcoil_variables.n_pancake,
                tfcoil_variables.thwcndut,
                tfcoil_variables.thicndut,
            )

        # Areas and fractions
        # -------------------
        # Central helium channel down the conductor core [m2]
        tfcoil_variables.awphec = (
            0.25e0 * tfcoil_variables.n_tf_turn * np.pi * tfcoil_variables.dhecoil**2
        )

        # Total conductor cross-sectional area, taking account of void area
        # and central helium channel [m2]
        tfcoil_variables.acond = (
            tfcoil_variables.acstf
            * tfcoil_variables.n_tf_turn
            * (1.0e0 - tfcoil_variables.vftf)
            - tfcoil_variables.awphec
        )

        # Void area in conductor for He, not including central channel [m2]
        tfcoil_variables.avwp = (
            tfcoil_variables.acstf * tfcoil_variables.n_tf_turn * tfcoil_variables.vftf
        )

        # Area of inter-turn insulation: total [m2]
        tfcoil_variables.a_tf_coil_wp_turn_insulation = (
            tfcoil_variables.n_tf_turn * tfcoil_variables.insulation_area
        )

        # Area of steel structure in winding pack [m2]
        tfcoil_variables.aswp = tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf

        # Inboard coil steel area [m2]
        sctfcoil_module.a_tf_steel = tfcoil_variables.acasetf + tfcoil_variables.aswp

        # Inboard coil steel fraction [-]
        sctfcoil_module.f_tf_steel = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_steel
            / tfcoil_variables.a_tf_coil_inboard
        )

        # Inboard coil insulation cross-section [m2]
        sctfcoil_module.a_tf_ins = (
            tfcoil_variables.a_tf_coil_wp_turn_insulation + sctfcoil_module.a_ground_ins
        )

        #  Inboard coil insulation fraction [-]
        sctfcoil_module.f_tf_ins = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_ins
            / tfcoil_variables.a_tf_coil_inboard
        )

        # Negative areas or fractions error reporting
        if (
            tfcoil_variables.acond <= 0.0e0
            or tfcoil_variables.avwp <= 0.0e0
            or tfcoil_variables.a_tf_coil_wp_turn_insulation <= 0.0e0
            or tfcoil_variables.aswp <= 0.0e0
            or sctfcoil_module.a_tf_steel <= 0.0e0
            or sctfcoil_module.f_tf_steel <= 0.0e0
            or sctfcoil_module.a_tf_ins <= 0.0e0
            or sctfcoil_module.f_tf_ins <= 0.0e0
        ):
            error_handling.fdiags[0] = tfcoil_variables.acond
            error_handling.fdiags[1] = tfcoil_variables.avwp
            error_handling.fdiags[2] = tfcoil_variables.a_tf_coil_wp_turn_insulation
            error_handling.fdiags[3] = tfcoil_variables.aswp
            error_handling.fdiags[4] = sctfcoil_module.a_tf_steel
            error_handling.fdiags[5] = sctfcoil_module.f_tf_steel
            error_handling.fdiags[6] = sctfcoil_module.a_tf_ins
            error_handling.fdiags[7] = sctfcoil_module.f_tf_ins
            error_handling.report_error(276)

    def tf_wp_geom(self, i_tf_wp_geom):
        """
         Author : S. Kahn, CCFE
        Seting the WP geometry and area for SC magnets
        """
        sctfcoil_module.r_wp_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.dr_tf_nose_case
        )

        # Radial position of outer edge of winding pack [m]
        sctfcoil_module.r_wp_outer = (
            sctfcoil_module.r_wp_inner + tfcoil_variables.dr_tf_wp
        )

        # Radius of geometrical centre of winding pack [m]
        sctfcoil_module.r_wp_centre = 0.5e0 * (
            sctfcoil_module.r_wp_inner + sctfcoil_module.r_wp_outer
        )

        # TF toroidal thickness at the WP inner radius [m]
        t_tf_at_wp = 2.0e0 * sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil

        # Minimal toroidal thickness of winding pack [m]
        sctfcoil_module.t_wp_toroidal = (
            t_tf_at_wp - 2.0e0 * tfcoil_variables.dx_tf_side_case
        )

        # Rectangular WP
        # --------------
        if i_tf_wp_geom == 0:
            # Outer WP layer toroidal thickness [m]
            tfcoil_variables.wwp1 = sctfcoil_module.t_wp_toroidal

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = sctfcoil_module.t_wp_toroidal

            # Total cross-sectional area of winding pack [m2]
            sctfcoil_module.awpc = (
                tfcoil_variables.dr_tf_wp * sctfcoil_module.t_wp_toroidal
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                tfcoil_variables.dr_tf_wp
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            ) * (
                sctfcoil_module.t_wp_toroidal
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap
            ) * (
                sctfcoil_module.t_wp_toroidal - 2.0e0 * tfcoil_variables.tfinsgap
            ) - sctfcoil_module.awptf

        # Double rectangular WP
        # ---------------------
        elif i_tf_wp_geom == 1:
            # Thickness of winding pack section at R > sctfcoil_module.r_wp_centre [m]
            tfcoil_variables.wwp1 = 2.0e0 * (
                sctfcoil_module.r_wp_centre * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.dx_tf_side_case
            )

            # Thickness of winding pack section at R < sctfcoil_module.r_wp_centre [m]
            tfcoil_variables.wwp2 = 2.0e0 * (
                sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.dx_tf_side_case
            )

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = 0.5e0 * (
                tfcoil_variables.wwp1 + tfcoil_variables.wwp2
            )

            # Total cross-sectional area of winding pack [m2]
            # Including ground insulation and insertion gap
            sctfcoil_module.awpc = (
                tfcoil_variables.dr_tf_wp * sctfcoil_module.t_wp_toroidal_av
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                0.5e0
                * (
                    tfcoil_variables.dr_tf_wp
                    - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                )
                * (
                    tfcoil_variables.wwp1
                    + tfcoil_variables.wwp2
                    - 4.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                )
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                0.5e0
                * (tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap)
                * (
                    tfcoil_variables.wwp1
                    + tfcoil_variables.wwp2
                    - 4.0e0 * tfcoil_variables.tfinsgap
                )
                - sctfcoil_module.awptf
            )

        # Trapezoidal WP
        # --------------
        else:
            # Thickness of winding pack section at sctfcoil_module.r_wp_outer [m]
            tfcoil_variables.wwp1 = 2.0e0 * (
                sctfcoil_module.r_wp_outer * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.dx_tf_side_case
            )

            # Thickness of winding pack section at sctfcoil_module.r_wp_inner [m]
            tfcoil_variables.wwp2 = 2.0e0 * (
                sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.dx_tf_side_case
            )

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = 0.5e0 * (
                tfcoil_variables.wwp1 + tfcoil_variables.wwp2
            )

            # Total cross-sectional area of winding pack [m2]
            # Including ground insulation and insertion gap
            sctfcoil_module.awpc = tfcoil_variables.dr_tf_wp * (
                tfcoil_variables.wwp2
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                tfcoil_variables.dr_tf_wp
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            ) * (
                tfcoil_variables.wwp2
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap
            ) * (
                tfcoil_variables.wwp2
                - 2.0e0 * tfcoil_variables.tfinsgap
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            ) - sctfcoil_module.awptf

        # --------------

        # Negative WP area error reporting
        if sctfcoil_module.awptf <= 0.0e0 or sctfcoil_module.awpc <= 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awptf
            error_handling.fdiags[1] = sctfcoil_module.awpc
            error_handling.report_error(99)

    def tf_case_geom(self, i_tf_wp_geom, i_tf_case_geom):
        """
        Author : S. Kahn, CCFE
        Setting the case geometry and area for SC magnets
        """
        tfcoil_variables.acasetf = (
            tfcoil_variables.a_tf_coil_inboard / tfcoil_variables.n_tf_coils
        ) - sctfcoil_module.awpc

        # Outboard leg cross-sectional area of surrounding case [m2]
        tfcoil_variables.acasetfo = (
            tfcoil_variables.a_tf_leg_outboard - sctfcoil_module.awpc
        )

        # Front casing area [m2]
        if i_tf_case_geom == 0:
            # Circular front case
            sctfcoil_module.a_case_front = (
                sctfcoil_module.rad_tf_coil_toroidal
                * build_variables.r_tf_inboard_out**2
                - sctfcoil_module.tan_theta_coil * sctfcoil_module.r_wp_outer**2
            )
        else:
            # Straight front case
            sctfcoil_module.a_case_front = (
                (sctfcoil_module.r_wp_outer + tfcoil_variables.casthi) ** 2
                - sctfcoil_module.r_wp_outer**2
            ) * sctfcoil_module.tan_theta_coil

        # Nose casing area [m2]
        sctfcoil_module.a_case_nose = (
            sctfcoil_module.tan_theta_coil * sctfcoil_module.r_wp_inner**2
            - sctfcoil_module.rad_tf_coil_toroidal * build_variables.r_tf_inboard_in**2
        )

        # Report error if the casing area is negative
        if tfcoil_variables.acasetf <= 0.0e0 or tfcoil_variables.acasetfo <= 0.0e0:
            error_handling.fdiags[0] = tfcoil_variables.acasetf
            error_handling.fdiags[1] = tfcoil_variables.acasetfo
            error_handling.report_error(99)

        # Average lateral casing thickness
        # --------------
        # Rectangular casing
        if i_tf_wp_geom == 0:
            sctfcoil_module.t_lat_case_av = (
                tfcoil_variables.dx_tf_side_case
                + 0.5e0 * sctfcoil_module.tan_theta_coil * tfcoil_variables.dr_tf_wp
            )

        # Double rectangular WP
        elif i_tf_wp_geom == 1:
            sctfcoil_module.t_lat_case_av = (
                tfcoil_variables.dx_tf_side_case
                + 0.25e0 * sctfcoil_module.tan_theta_coil * tfcoil_variables.dr_tf_wp
            )

        # Trapezoidal WP
        else:
            sctfcoil_module.t_lat_case_av = tfcoil_variables.dx_tf_side_case

        # --------------

    def tf_integer_turn_geom(self, n_layer, n_pancake, thwcndut, thicndut):
        """
        Authors: J Morris & S Khan
        Setting the TF WP turn geometry for SC magnets from the number
        of turns rows in the radial direction. The turns can have any
        rectangular shapes.
        This calculation has two purposes, first to check if a turn can exist
        (positive cable space) and the second to provide its dimenesions,
        areas and the its associated current

        """
        sctfcoil_module.rbcndut = thwcndut * 0.75e0

        # Radial turn dimension [m]
        sctfcoil_module.t_turn_radial = (
            tfcoil_variables.dr_tf_wp
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
        ) / n_layer

        if sctfcoil_module.t_turn_radial <= (2.0e0 * thicndut + 2.0e0 * thwcndut):
            error_handling.fdiags[0] = sctfcoil_module.t_turn_radial
            error_handling.fdiags[1] = thicndut
            error_handling.fdiags[2] = thwcndut
            error_handling.report_error(100)

        # Toroidal turn dimension [m]
        sctfcoil_module.t_turn_toroidal = (
            sctfcoil_module.t_wp_toroidal
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
        ) / n_pancake

        if sctfcoil_module.t_turn_toroidal <= (2.0e0 * thicndut + 2.0e0 * thwcndut):
            error_handling.fdiags[0] = sctfcoil_module.t_turn_toroidal
            error_handling.fdiags[1] = thicndut
            error_handling.fdiags[2] = thwcndut
            error_handling.report_error(100)

        tfcoil_variables.t_turn_tf = np.sqrt(
            sctfcoil_module.t_turn_radial * sctfcoil_module.t_turn_toroidal
        )

        # Number of TF turns
        n_tf_turn = np.double(n_layer * n_pancake)

        # Current per turn [A/turn]
        cpttf = sctfcoil_module.c_tf_coil / n_tf_turn

        # Radial and toroidal dimension of conductor [m]
        sctfcoil_module.t_conductor_radial = (
            sctfcoil_module.t_turn_radial - 2.0e0 * thicndut
        )
        sctfcoil_module.t_conductor_toroidal = (
            sctfcoil_module.t_turn_toroidal - 2.0e0 * thicndut
        )
        tfcoil_variables.t_conductor = np.sqrt(
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
        )

        # Dimension of square cable space inside conduit [m]
        sctfcoil_module.t_cable_radial = (
            sctfcoil_module.t_conductor_radial - 2.0e0 * thwcndut
        )
        sctfcoil_module.t_cable_toroidal = (
            sctfcoil_module.t_conductor_toroidal - 2.0e0 * thwcndut
        )
        sctfcoil_module.t_cable = np.sqrt(
            sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal
        )

        # Cross-sectional area of cable space per turn
        # taking account of rounded inside corners [m2]
        acstf = (sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal) - (
            4.0e0 - np.pi
        ) * sctfcoil_module.rbcndut**2

        if acstf <= 0.0e0:
            if (sctfcoil_module.t_cable_radial < 0.0e0) or (
                sctfcoil_module.t_cable_toroidal < 0.0e0
            ):
                error_handling.fdiags[0] = acstf
                error_handling.fdiags[1] = sctfcoil_module.t_cable_radial
                error_handling.fdiags[2] = sctfcoil_module.t_cable_toroidal
                error_handling.report_error(101)
            else:
                error_handling.fdiags[0] = acstf
                error_handling.fdiags[1] = sctfcoil_module.t_cable_radial
                error_handling.fdiags[1] = sctfcoil_module.t_cable_toroidal
                error_handling.report_error(102)
                sctfcoil_module.rbcndut = 0.0e0
                acstf = (
                    sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal
                )

        # Cross-sectional area of conduit jacket per turn [m2]
        acndttf = (
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
            - acstf
        )

        # Area of inter-turn insulation: single turn [m2]
        insulation_area = (
            sctfcoil_module.t_turn_radial * sctfcoil_module.t_turn_toroidal
            - acndttf
            - acstf
        )
        return acstf, acndttf, insulation_area, cpttf, n_tf_turn

        # -------------

    def tf_wp_currents(self):
        """
        Author : S. Kahn, CCFE
        Turn engineering turn currents/densities
        """
        tfcoil_variables.j_tf_wp = max(
            1.0e0,
            tfcoil_variables.c_tf_total
            / (tfcoil_variables.n_tf_coils * sctfcoil_module.awptf),
        )

    def outtf(self, peaktfflag):
        """Writes superconducting TF coil output to file
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        peaktfflag : input integer : warning flag from peak TF calculation
        This routine writes the superconducting TF coil results
        to the output file.
        PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
        """

        # General coil parameters
        po.osubhd(self.outfile, "TF design")
        po.ovarin(
            self.outfile,
            "Conductor technology",
            "(i_tf_sup)",
            tfcoil_variables.i_tf_sup,
        )

        if tfcoil_variables.i_tf_sup == 0:
            po.ocmmnt(
                self.outfile, "  -> Resitive coil : Water cooled copper (GLIDCOP AL-15)"
            )
        elif tfcoil_variables.i_tf_sup == 1:
            po.ocmmnt(self.outfile, "  -> Superconducting coil (SC)")
        elif tfcoil_variables.i_tf_sup == 2:
            po.ocmmnt(self.outfile, "  -> Reisitive coil : Helium cooled aluminium")

        # SC material scaling
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarin(
                self.outfile,
                "Superconductor material",
                "(i_tf_sc_mat)",
                tfcoil_variables.i_tf_sc_mat,
            )

            if tfcoil_variables.i_tf_sc_mat == 1:
                po.ocmmnt(self.outfile, "  -> ITER Nb3Sn critical surface model")
            elif tfcoil_variables.i_tf_sc_mat == 2:
                po.ocmmnt(self.outfile, "  -> Bi-2212 high temperature superconductor")
            elif tfcoil_variables.i_tf_sc_mat == 3:
                po.ocmmnt(self.outfile, "  -> NbTi")
            elif tfcoil_variables.i_tf_sc_mat == 4:
                po.ocmmnt(
                    self.outfile,
                    "  -> ITER Nb3Sn critical surface model, user-defined parameters",
                )
            elif tfcoil_variables.i_tf_sc_mat == 5:
                po.ocmmnt(self.outfile, "  -> WST Nb3Sn")
            elif tfcoil_variables.i_tf_sc_mat == 6:
                po.ocmmnt(
                    self.outfile,
                    "  -> High temperature superconductor: REBCO HTS tape in CroCo strand",
                )
            elif tfcoil_variables.i_tf_sc_mat == 7:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Durham Ginzburg-Landau critical surface model for Nb-Ti",
                )
            elif tfcoil_variables.i_tf_sc_mat == 8:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Durham Ginzburg-Landau critical surface model for REBCO",
                )
            elif tfcoil_variables.i_tf_sc_mat == 9:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Hazelton experimental data + Zhai conceptual model for REBCO",
                )

        # Joints strategy
        po.ovarin(
            self.outfile,
            "Presence of TF demountable joints",
            "(itart)",
            physics_variables.itart,
        )
        if physics_variables.itart == 1:
            po.ocmmnt(
                self.outfile, "  -> TF coil made of a Centerpost (CP) and outer legs"
            )
            po.ocmmnt(self.outfile, "     interfaced with demountable joints")
        else:
            po.ocmmnt(self.outfile, "  -> Coils without demountable joints")

        # Centring forces support strategy
        po.ovarin(
            self.outfile,
            "TF inboard leg support strategy",
            "(i_tf_bucking)",
            tfcoil_variables.i_tf_bucking,
        )

        if tfcoil_variables.i_tf_bucking == 0:
            po.ocmmnt(self.outfile, "  -> No support structure")
        elif tfcoil_variables.i_tf_bucking == 1:
            if tfcoil_variables.i_tf_sup == 1:
                po.ocmmnt(self.outfile, "  -> Steel casing")
            elif (
                abs(tfcoil_variables.eyoung_res_tf_buck - 205.0e9)
                < np.finfo(float(tfcoil_variables.eyoung_res_tf_buck)).eps
            ):
                po.ocmmnt(self.outfile, "  -> Steel bucking cylinder")
            else:
                po.ocmmnt(self.outfile, "  -> Bucking cylinder")

        elif (
            tfcoil_variables.i_tf_bucking in (2, 3)
            and build_variables.i_tf_inside_cs == 1
        ):
            po.ocmmnt(
                self.outfile,
                "  -> TF in contact with dr_bore filler support (bucked and weged design)",
            )

        elif (
            tfcoil_variables.i_tf_bucking in (2, 3)
            and build_variables.i_tf_inside_cs == 0
        ):
            po.ocmmnt(
                self.outfile, "  -> TF in contact with CS (bucked and weged design)"
            )

        # TF coil geometry
        po.osubhd(self.outfile, "TF coil Geometry :")
        po.ovarin(
            self.outfile,
            "Number of TF coils",
            "(n_tf_coils)",
            int(tfcoil_variables.n_tf_coils),
        )
        po.ovarre(
            self.outfile,
            "Inboard leg centre radius (m)",
            "(r_tf_inboard_mid)",
            build_variables.r_tf_inboard_mid,
            "OP ",
        )
        po.ovarre(
            constants.mfile,
            "Inboard leg inner radius (m)",
            "(r_tf_inboard_in)",
            build_variables.r_tf_inboard_in,
            "OP ",
        )
        po.ovarre(
            constants.mfile,
            "Inboard leg outer radius (m)",
            "(r_tf_inboard_out)",
            build_variables.r_tf_inboard_out,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "WP shape selection switch",
            "(i_tf_wp_geom)",
            tfcoil_variables.i_tf_wp_geom,
        )
        po.ovarre(
            constants.mfile,
            "Radial position of inner edge and centre of winding pack (m)",
            "(r_wp_inner)",
            sctfcoil_module.r_wp_inner,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Outboard leg centre radius (m)",
            "(r_tf_outboard_mid)",
            build_variables.r_tf_outboard_mid,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "Outboard leg nose case type",
            "(i_tf_case_geom)",
            tfcoil_variables.i_tf_case_geom,
        )
        po.ovarre(
            self.outfile,
            "Total inboard leg radial thickness (m)",
            "(dr_tf_inboard)",
            build_variables.dr_tf_inboard,
        )
        po.ovarre(
            self.outfile,
            "Total outboard leg radial thickness (m)",
            "(dr_tf_outboard)",
            build_variables.dr_tf_outboard,
        )
        po.ovarre(
            self.outfile,
            "Outboard leg toroidal thickness (m)",
            "(dx_tf_inboard_out_toroidal)",
            tfcoil_variables.dx_tf_inboard_out_toroidal,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum inboard edge height (m)",
            "(hmax)",
            build_variables.hmax,
            "OP ",
        )
        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Mean coil circumference (inboard leg not included) (m)",
                "(len_tf_coil)",
                tfcoil_variables.len_tf_coil,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Length of the inboard segment (m)",
                "(cplen)",
                tfcoil_variables.cplen,
                "OP ",
            )
        else:
            po.ovarre(
                self.outfile,
                "Mean coil circumference (including inboard leg length) (m)",
                "(len_tf_coil)",
                tfcoil_variables.len_tf_coil,
                "OP ",
            )

        # Vertical shape
        po.ovarin(
            self.outfile,
            "Vertical TF shape",
            "(i_tf_shape)",
            tfcoil_variables.i_tf_shape,
        )
        if tfcoil_variables.i_tf_shape == 1:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "D-shape coil, inner surface shape approximated by")
            po.ocmmnt(
                self.outfile,
                "by a straight segment and elliptical arcs between the following points:",
            )
            po.oblnkl(self.outfile)
        elif tfcoil_variables.i_tf_shape == 2:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Picture frame coil, inner surface approximated by")
            po.ocmmnt(
                self.outfile, "by a straight segment between the following points:"
            )
            po.oblnkl(self.outfile)

        po.write(self.outfile, "  point              x(m)              y(m)")
        for ii in range(5):
            po.write(
                self.outfile,
                f"  {ii}              {tfcoil_variables.xarc[ii]}              {tfcoil_variables.yarc[ii]}",
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} R (m)",
                f"(xarc({ii + 1}))",
                tfcoil_variables.xarc[ii],
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} Z (m)",
                f"(yarc({ii + 1}))",
                tfcoil_variables.yarc[ii],
            )

        # CP tapering geometry
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Tapered Centrepost TF coil Dimensions:")
            po.ovarre(
                self.outfile,
                "TF coil centrepost outer radius at midplane (m)",
                "(r_tf_inboard_out)",
                build_variables.r_tf_inboard_out,
            )
            po.ovarre(
                self.outfile,
                "TF coil centrepost outer radius at its top (m)",
                "(r_cp_top)",
                build_variables.r_cp_top,
            )
            po.ovarre(
                self.outfile,
                "Top/miplane TF CP radius ratio (-)",
                "(f_r_cp)",
                build_variables.f_r_cp,
            )
            po.ovarre(
                self.outfile,
                "Distance from the midplane to the top of the tapered section (m)",
                "(z_cp_top)",
                sctfcoil_module.z_cp_top,
            )
            po.ovarre(
                self.outfile,
                "Distance from the midplane to the top of the centrepost (m)",
                "(hmax + dr_tf_outboard)",
                build_variables.hmax + build_variables.dr_tf_outboard,
            )

        # Turn/WP gemoetry
        if tfcoil_variables.i_tf_sup == 1:
            # Total material fraction
            po.osubhd(self.outfile, "Global material area/fractions:")
            po.ovarre(
                self.outfile,
                "TF cross-section (total) (m2)",
                "(a_tf_coil_inboard)",
                tfcoil_variables.a_tf_coil_inboard,
            )
            po.ovarre(
                self.outfile,
                "Total steel cross-section (m2)",
                "(a_tf_steel*n_tf_coils)",
                sctfcoil_module.a_tf_steel * tfcoil_variables.n_tf_coils,
            )
            po.ovarre(
                self.outfile,
                "Total steel TF fraction",
                "(f_tf_steel)",
                sctfcoil_module.f_tf_steel,
            )
            po.ovarre(
                self.outfile,
                "Total Insulation cross-section (total) (m2)",
                "(a_tf_ins*n_tf_coils)",
                sctfcoil_module.a_tf_ins * tfcoil_variables.n_tf_coils,
            )
            po.ovarre(
                self.outfile,
                "Total Insulation fraction",
                "(f_tf_ins)",
                sctfcoil_module.f_tf_ins,
            )

            # External casing
            po.osubhd(self.outfile, "External steel Case Information :")
            po.ovarre(
                self.outfile,
                "Casing cross section area (per leg) (m2)",
                "(acasetf)",
                tfcoil_variables.acasetf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case plasma side wall thickness (m)",
                "(casthi)",
                tfcoil_variables.casthi,
            )
            po.ovarre(
                self.outfile,
                'Inboard leg case inboard "nose" thickness (m)',
                "(dr_tf_nose_case)",
                tfcoil_variables.dr_tf_nose_case,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case sidewall thickness at its narrowest point (m)",
                "(dx_tf_side_case)",
                tfcoil_variables.dx_tf_side_case,
            )
            po.ovarre(
                self.outfile,
                "External case mass per coil (kg)",
                "(whtcas)",
                tfcoil_variables.whtcas,
                "OP ",
            )

            # Winding pack structure
            po.osubhd(self.outfile, "TF winding pack (WP) geometry:")
            po.ovarre(
                self.outfile,
                "WP cross section area with insulation and insertion (per coil) (m2)",
                "(awpc)",
                sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "WP cross section area (per coil) (m2)",
                "(aswp)",
                sctfcoil_module.awptf,
            )
            po.ovarre(
                self.outfile,
                "Winding pack radial thickness (m)",
                "(dr_tf_wp)",
                tfcoil_variables.dr_tf_wp,
                "OP ",
            )
            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width (m)",
                    "(wwp1)",
                    tfcoil_variables.wwp1,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width 1 (m)",
                    "(wwp1)",
                    tfcoil_variables.wwp1,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width 2 (m)",
                    "(wwp2)",
                    tfcoil_variables.wwp2,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Ground wall insulation thickness (m)",
                "(tinstf)",
                tfcoil_variables.tinstf,
            )
            po.ovarre(
                self.outfile,
                "Winding pack insertion gap (m)",
                "(tfinsgap)",
                tfcoil_variables.tfinsgap,
            )

            # WP material fraction
            po.osubhd(self.outfile, "TF winding pack (WP) material area/fractions:")
            po.ovarre(
                self.outfile,
                "Steel WP cross-section (total) (m2)",
                "(aswp*n_tf_coils)",
                tfcoil_variables.aswp * tfcoil_variables.n_tf_coils,
            )
            po.ovarre(
                self.outfile,
                "Steel WP fraction",
                "(aswp/awpc)",
                tfcoil_variables.aswp / sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Insulation WP fraction",
                "(a_tf_coil_wp_turn_insulation/awpc)",
                tfcoil_variables.a_tf_coil_wp_turn_insulation / sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Cable WP fraction",
                "((awpc-aswp-a_tf_coil_wp_turn_insulation)/awpc)",
                (
                    sctfcoil_module.awpc
                    - tfcoil_variables.aswp
                    - tfcoil_variables.a_tf_coil_wp_turn_insulation
                )
                / sctfcoil_module.awpc,
            )

            # Number of turns
            po.osubhd(self.outfile, "WP turn information:")
            po.ovarin(
                self.outfile,
                "Turn parametrisation",
                "(i_tf_turns_integer)",
                tfcoil_variables.i_tf_turns_integer,
            )
            if tfcoil_variables.i_tf_turns_integer == 0:
                po.ocmmnt(self.outfile, "  Non-integer number of turns")
            else:
                po.ocmmnt(self.outfile, "  Integer number of turns")

            po.ovarre(
                self.outfile,
                "Number of turns per TF coil",
                "(n_tf_turn)",
                tfcoil_variables.n_tf_turn,
                "OP ",
            )
            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarin(
                    self.outfile,
                    "Number of TF pancakes",
                    "(n_pancake)",
                    tfcoil_variables.n_pancake,
                )
                po.ovarin(
                    self.outfile,
                    "Number of TF layers",
                    "(n_layer)",
                    tfcoil_variables.n_layer,
                )

            po.oblnkl(self.outfile)

            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarre(
                    self.outfile,
                    "Radial width of turn (m)",
                    "(t_turn_radial)",
                    sctfcoil_module.t_turn_radial,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of turn (m)",
                    "(t_turn_toroidal)",
                    sctfcoil_module.t_turn_toroidal,
                )
                po.ovarre(
                    self.outfile,
                    "Radial width of conductor (m)",
                    "(t_conductor_radial)",
                    sctfcoil_module.t_conductor_radial,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of conductor (m)",
                    "(t_conductor_toroidal)",
                    sctfcoil_module.t_conductor_toroidal,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Radial width of cable space",
                    "(t_cable_radial)",
                    sctfcoil_module.t_cable_radial,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of cable space",
                    "(t_cable_toroidal)",
                    sctfcoil_module.t_cable_toroidal,
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Width of turn including inter-turn insulation (m)",
                    "(t_turn_tf)",
                    tfcoil_variables.t_turn_tf,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Width of conductor (square) (m)",
                    "(t_conductor)",
                    tfcoil_variables.t_conductor,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Width of space inside conductor (m)",
                    "(t_cable)",
                    sctfcoil_module.t_cable,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Steel conduit thickness (m)",
                "(thwcndut)",
                tfcoil_variables.thwcndut,
            )
            po.ovarre(
                self.outfile,
                "Inter-turn insulation thickness (m)",
                "(thicndut)",
                tfcoil_variables.thicndut,
            )

            if tfcoil_variables.i_tf_sc_mat in (1, 2, 3, 4, 5, 7, 8, 9):
                po.osubhd(self.outfile, "Conductor information:")
                po.ovarre(
                    self.outfile,
                    "Diameter of central helium channel in cable",
                    "(dhecoil)",
                    tfcoil_variables.dhecoil,
                )
                po.ocmmnt(self.outfile, "Fractions by area")
                po.ovarre(
                    self.outfile,
                    "internal area of the cable space",
                    "(acstf)",
                    tfcoil_variables.acstf,
                )
                po.ovarre(
                    self.outfile,
                    "Coolant fraction in conductor excluding central channel",
                    "(vftf)",
                    tfcoil_variables.vftf,
                )
                po.ovarre(
                    self.outfile,
                    "Copper fraction of conductor",
                    "(fcutfsu)",
                    tfcoil_variables.fcutfsu,
                )
                po.ovarre(
                    self.outfile,
                    "Superconductor fraction of conductor",
                    "(1-fcutfsu)",
                    1 - tfcoil_variables.fcutfsu,
                )
                # TODO
                # po.ovarre(self.outfile,'Conductor fraction of winding pack','(tfcoil_variables.acond/ap)',acond/ap, 'OP ')
                # po.ovarre(self.outfile,'Conduit fraction of winding pack','(tfcoil_variables.n_tf_turn*tfcoil_variables.acndttf/ap)',n_tf_turn*tfcoil_variables.acndttf/ap, 'OP ')
                # po.ovarre(self.outfile,'Insulator fraction of winding pack','(tfcoil_variables.a_tf_coil_wp_turn_insulation/ap)',a_tf_coil_wp_turn_insulation/ap, 'OP ')
                # po.ovarre(self.outfile,'Helium area fraction of winding pack excluding central channel','(tfcoil_variables.avwp/ap)',avwp/ap, 'OP ')
                # po.ovarre(self.outfile,'Central helium channel area as fraction of winding pack','(tfcoil_variables.awphec/ap)',awphec/ap, 'OP ')
                ap = (
                    tfcoil_variables.acond
                    + tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf
                    + tfcoil_variables.a_tf_coil_wp_turn_insulation
                    + tfcoil_variables.avwp
                    + tfcoil_variables.awphec
                )
                po.ovarrf(
                    self.outfile,
                    "Check total area fractions in winding pack = 1",
                    "",
                    (
                        tfcoil_variables.acond
                        + tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf
                        + tfcoil_variables.a_tf_coil_wp_turn_insulation
                        + tfcoil_variables.avwp
                        + tfcoil_variables.awphec
                    )
                    / ap,
                )
                po.ovarrf(
                    self.outfile,
                    "minimum TF conductor temperature margin  (K)",
                    "(tmargmin_tf)",
                    tfcoil_variables.tmargmin_tf,
                )
                po.ovarrf(
                    self.outfile,
                    "TF conductor temperature margin (K)",
                    "(tmargtf)",
                    tfcoil_variables.tmargtf,
                )

                po.ovarin(
                    self.outfile,
                    "Elastic properties behavior",
                    "(i_tf_cond_eyoung_axial)",
                    tfcoil_variables.i_tf_cond_eyoung_axial,
                )
                if tfcoil_variables.i_tf_cond_eyoung_axial == 0:
                    po.ocmmnt(self.outfile, "  Conductor stiffness neglected")
                elif tfcoil_variables.i_tf_cond_eyoung_axial == 1:
                    po.ocmmnt(self.outfile, "  Conductor stiffness is user-input")
                elif tfcoil_variables.i_tf_cond_eyoung_axial == 2:
                    po.ocmmnt(
                        self.outfile,
                        "  Conductor stiffness is set by material-specific default",
                    )

                po.ovarre(
                    self.outfile,
                    "Conductor axial Youngs modulus",
                    "(eyoung_cond_axial)",
                    tfcoil_variables.eyoung_cond_axial,
                )
                po.ovarre(
                    self.outfile,
                    "Conductor transverse Youngs modulus",
                    "(eyoung_cond_trans)",
                    tfcoil_variables.eyoung_cond_trans,
                )
        else:
            # External casing
            po.osubhd(self.outfile, "Bucking cylinder information:")
            po.ovarre(
                self.outfile,
                "Casing cross section area (per leg) (m2)",
                "(acasetf)",
                tfcoil_variables.acasetf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case plasma side wall thickness (m)",
                "(casthi)",
                tfcoil_variables.casthi,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg bucking cylinder thickness (m)",
                "(dr_tf_nose_case)",
                tfcoil_variables.dr_tf_nose_case,
            )

            # Conductor layer geometry
            po.osubhd(self.outfile, "Inboard TFC conductor sector geometry:")
            po.ovarre(
                self.outfile,
                "Inboard TFC conductor sector area with gr insulation (per leg) (m2)",
                "(awpc)",
                sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Inboard TFC conductor sector area, NO ground & gap (per leg) (m2)",
                "(awptf)",
                sctfcoil_module.awptf,
            )
            po.ovarre(
                self.outfile,
                "Inboard conductor sector radial thickness (m)",
                "(dr_tf_wp)",
                tfcoil_variables.dr_tf_wp,
            )
            if physics_variables.itart == 1:
                po.ovarre(
                    self.outfile,
                    "Central collumn top conductor sector radial thickness (m)",
                    "(dr_tf_wp_top)",
                    sctfcoil_module.dr_tf_wp_top,
                )

            po.ovarre(
                self.outfile,
                "Ground wall insulation thickness (m)",
                "(tinstf)",
                tfcoil_variables.tinstf,
            )
            # Turn info
            po.osubhd(self.outfile, "Coil turn information:")
            po.ovarre(
                self.outfile,
                "Number of turns per TF leg",
                "(n_tf_turn)",
                tfcoil_variables.n_tf_turn,
            )
            po.ovarre(
                self.outfile,
                "Turn insulation thickness",
                "(thicndut)",
                tfcoil_variables.thicndut,
            )
            po.ovarre(
                self.outfile,
                "Mid-plane CP cooling fraction",
                "(fcoolcp)",
                tfcoil_variables.fcoolcp,
            )
            po.ovarre(
                self.outfile,
                "Outboard leg current per turn (A)",
                "(cpttf)",
                tfcoil_variables.cpttf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg conductor volume (m3)",
                "(vol_cond_cp)",
                tfcoil_variables.vol_cond_cp,
            )
            po.ovarre(
                self.outfile,
                "Outboard leg volume per coil (m3)",
                "(voltfleg)",
                tfcoil_variables.voltfleg,
            )

        # Coil masses
        po.osubhd(self.outfile, "TF coil mass:")
        po.ovarre(
            self.outfile,
            "Superconductor mass per coil (kg)",
            "(whtconsc)",
            tfcoil_variables.whtconsc,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Copper mass per coil (kg)",
            "(whtconcu)",
            tfcoil_variables.whtconcu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Steel conduit mass per coil (kg)",
            "(m_tf_turn_steel_conduit)",
            tfcoil_variables.m_tf_turn_steel_conduit,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Conduit insulation mass per coil (kg)",
            "(whtconin)",
            tfcoil_variables.whtconin,
            "OP ",
        )
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarre(
                self.outfile,
                "Total conduit mass per coil (kg)",
                "(whtcon)",
                tfcoil_variables.whtcon,
                "OP ",
            )

        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Mass of inboard legs (kg)",
                "(whtcp)",
                tfcoil_variables.whtcp,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mass of outboard legs (kg)",
                "(whttflgs)",
                tfcoil_variables.whttflgs,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Mass of each TF coil (kg)",
            "(m_tf_coils_total/n_tf_coils)",
            tfcoil_variables.m_tf_coils_total / tfcoil_variables.n_tf_coils,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total TF coil mass (kg)",
            "(m_tf_coils_total)",
            tfcoil_variables.m_tf_coils_total,
            "OP ",
        )

        # TF current and field
        po.osubhd(self.outfile, "Maximum B field and currents:")
        po.ovarre(
            self.outfile,
            "Nominal peak field assuming toroidal symmetry (T)",
            "(b_tf_inboard_peak)",
            tfcoil_variables.b_tf_inboard_peak,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total current in all TF coils (MA)",
            "(c_tf_total/1.D6)",
            1.0e-6 * tfcoil_variables.c_tf_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "TF coil current (summed over all coils) (A)",
            "(c_tf_total)",
            tfcoil_variables.c_tf_total,
        )
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarre(
                self.outfile,
                "Actual peak field at discrete conductor (T)",
                "(bmaxtfrp)",
                tfcoil_variables.bmaxtfrp,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Winding pack current density (A/m2)",
                "(j_tf_wp)",
                tfcoil_variables.j_tf_wp,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Inboard leg mid-plane conductor current density (A/m2)",
            "(oacdcp)",
            tfcoil_variables.oacdcp,
        )
        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Outboard leg conductor current density (A/m2)",
                "(cdtfleg)",
                tfcoil_variables.cdtfleg,
            )

        po.ovarre(
            self.outfile,
            "Total stored energy in TF coils (GJ)",
            "(estotftgj)",
            tfcoil_variables.estotftgj,
            "OP ",
        )

        # TF forces
        po.osubhd(self.outfile, "TF Forces:")
        po.ovarre(
            self.outfile,
            "Inboard vertical tension per coil (N)",
            "(vforce)",
            tfcoil_variables.vforce,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Outboard vertical tension per coil (N)",
            "(vforce_outboard)",
            tfcoil_variables.vforce_outboard,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "inboard vertical tension fraction (-)",
            "(f_vforce_inboard)",
            tfcoil_variables.f_vforce_inboard,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Centring force per coil (N/m)",
            "(cforce)",
            tfcoil_variables.cforce,
            "OP ",
        )

        # Resistive coil parameters
        if tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Resitive loss parameters:")
            if tfcoil_variables.i_tf_sup == 0:
                po.ocmmnt(
                    self.outfile,
                    "Resistive Material : GLIDCOP AL-15 - Dispersion Strengthened Copper",
                )
            elif tfcoil_variables.i_tf_sup == 2:
                po.ocmmnt(
                    self.outfile, "Resistive Material : Pure Aluminium (99.999+ %)"
                )

            if physics_variables.itart == 1:
                po.ovarre(
                    self.outfile,
                    "CP resistivity (ohm.m)",
                    "(rho_cp)",
                    tfcoil_variables.rho_cp,
                )
                po.ovarre(
                    self.outfile,
                    "Leg resistivity (ohm.m)",
                    "(rho_tf_leg)",
                    tfcoil_variables.rho_tf_leg,
                )
                po.ovarre(
                    self.outfile,
                    "CP resistive power loss (W)",
                    "(p_cp_resistive)",
                    tfcoil_variables.p_cp_resistive,
                )
                po.ovarre(
                    self.outfile,
                    "Total legs resitive power loss, (W)",
                    "(p_tf_leg_resistive)",
                    tfcoil_variables.p_tf_leg_resistive,
                )
                po.ovarre(
                    self.outfile,
                    "joints resistive power loss (W)",
                    "(pres_joints)",
                    tfcoil_variables.pres_joints,
                )
                po.ovarre(
                    self.outfile,
                    "Outboard leg resistance per coil (ohm)",
                    "(res_tf_leg)",
                    tfcoil_variables.res_tf_leg,
                )
                po.ovarre(
                    self.outfile,
                    "Average CP temperature (K)",
                    "(temp_cp_average)",
                    tfcoil_variables.temp_cp_average,
                )
                po.ovarre(
                    self.outfile,
                    "Average leg temperature (K)",
                    "(temp_tf_legs_outboard)",
                    tfcoil_variables.temp_tf_legs_outboard,
                )

            else:
                po.ovarre(
                    self.outfile,
                    "TF resistivity (ohm.m)",
                    "(p_cp_resistive)",
                    tfcoil_variables.rho_cp,
                )
                po.ovarre(
                    self.outfile,
                    "TF coil resistive power less (total) (ohm.m)",
                    "(p_cp_resistive)",
                    tfcoil_variables.p_cp_resistive,
                )
                po.ovarre(
                    self.outfile,
                    "Average coil temperature (K)",
                    "(temp_cp_average)",
                    tfcoil_variables.temp_cp_average,
                )

        # Ripple calculations
        po.osubhd(self.outfile, "Ripple information:")
        if tfcoil_variables.i_tf_shape == 1:
            if peaktfflag == 1:
                error_handling.report_error(144)
            elif peaktfflag == 2:
                error_handling.report_error(145)

            po.ovarre(
                self.outfile,
                "Max allowed tfcoil_variables.ripple amplitude at plasma outboard midplane (%)",
                "(ripmax)",
                tfcoil_variables.ripmax,
            )
            po.ovarre(
                self.outfile,
                "Ripple amplitude at plasma outboard midplane (%)",
                "(ripple)",
                tfcoil_variables.ripple,
                "OP ",
            )
        else:
            po.ovarre(
                self.outfile,
                "Max allowed tfcoil_variables.ripple amplitude at plasma (%)",
                "(ripmax)",
                tfcoil_variables.ripmax,
            )
            po.ovarre(
                self.outfile,
                "Ripple at plasma edge (%)",
                "(ripple)",
                tfcoil_variables.ripple,
            )
            po.ocmmnt(
                self.outfile,
                "  Ripple calculation to be re-defined for picure frame coils",
            )

        # Quench information
        if tfcoil_variables.i_tf_sup == 1:
            po.osubhd(self.outfile, "Quench information :")
            po.ovarre(
                self.outfile,
                "Actual quench time (or time constant) (s)",
                "(tdmptf)",
                tfcoil_variables.tdmptf,
            )
            po.ovarre(
                self.outfile,
                "Vacuum Vessel stress on quench (Pa)",
                "(vv_stress_quench)",
                sctfcoil_module.vv_stress_quench,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Maximum allowed voltage during quench due to insulation (kV)",
                "(vdalw)",
                tfcoil_variables.vdalw,
            )
            po.ovarre(
                self.outfile,
                "Actual quench voltage (kV)",
                "(vtfskv)",
                tfcoil_variables.vtfskv,
                "OP ",
            )

            if tfcoil_variables.i_tf_sc_mat in (1, 2, 3, 4, 5):
                po.ovarre(
                    self.outfile,
                    "Maximum allowed temp rise during a quench (K)",
                    "(tmaxpro)",
                    tfcoil_variables.tmaxpro,
                )
            elif tfcoil_variables == 6:
                po.ocmmnt(self.outfile, "CroCo cable with jacket: ")

                if 75 in numerics.icc:
                    po.ovarre(
                        self.outfile,
                        "Maximum permitted TF coil current / copper area (A/m2)",
                        "(copperA_m2_max)",
                        rebco_variables.copperA_m2_max,
                    )

                po.ovarre(
                    self.outfile,
                    "Actual TF coil current / copper area (A/m2)",
                    "(copperA_m2)",
                    rebco_variables.copperA_m2,
                )

        # TF coil radial build
        po.osubhd(self.outfile, "Radial build of TF coil centre-line :")
        # po.write(self.outfile,5)
        # 5   format(t43,'Thickness (m)',t60,'Outer radius (m)')

        radius = build_variables.r_tf_inboard_in
        po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

        # Radial build for SC TF coils
        if tfcoil_variables.i_tf_sup == 1:
            radius = radius + tfcoil_variables.dr_tf_nose_case
            po.obuild(
                self.outfile,
                'Coil case ("nose")',
                tfcoil_variables.dr_tf_nose_case,
                radius,
                "(dr_tf_nose_case)",
            )

            radius = radius + tfcoil_variables.tfinsgap
            po.obuild(
                self.outfile,
                "Insertion gap for winding pack",
                tfcoil_variables.tfinsgap,
                radius,
                "(tfinsgap)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Winding pack ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius
                + 0.5e0 * tfcoil_variables.dr_tf_wp
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            po.obuild(
                self.outfile,
                "Winding - first half",
                tfcoil_variables.dr_tf_wp / 2e0
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap,
                radius,
                "(dr_tf_wp/2-tinstf-tfinsgap)",
            )

            radius = (
                radius
                + 0.5e0 * tfcoil_variables.dr_tf_wp
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            po.obuild(
                self.outfile,
                "Winding - second half",
                tfcoil_variables.dr_tf_wp / 2e0
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap,
                radius,
                "(dr_tf_wp/2-tinstf-tfinsgap)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Winding pack insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.tfinsgap
            po.obuild(
                self.outfile,
                "Insertion gap for winding pack",
                tfcoil_variables.tfinsgap,
                radius,
                "(tfinsgap)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side case min radius",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

            radius = radius / np.cos(np.pi / tfcoil_variables.n_tf_coils)
            po.obuild(
                self.outfile,
                "Plasma side case max radius",
                build_variables.r_tf_inboard_out,
                radius,
                "(r_tf_inboard_out)",
            )

        # Radial build for restive coil
        else:
            radius = radius + tfcoil_variables.dr_tf_nose_case
            po.obuild(
                self.outfile,
                "Coil bucking cylindre",
                tfcoil_variables.dr_tf_nose_case,
                radius,
                "(dr_tf_nose_case)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius + 0.5e0 * tfcoil_variables.dr_tf_wp - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - first half",
                tfcoil_variables.dr_tf_wp / 2e0 - tfcoil_variables.tinstf,
                radius,
                "(tfcoil_variables.dr_tf_wp/2-tfcoil_variables.tinstf)",
            )

            radius = (
                radius + 0.5e0 * tfcoil_variables.dr_tf_wp - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - second half",
                tfcoil_variables.dr_tf_wp / 2e0 - tfcoil_variables.tinstf,
                radius,
                "(tfcoil_variables.dr_tf_wp/2-tfcoil_variables.tinstf)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

        # Radial build consistency check
        if (
            abs(
                radius - build_variables.r_tf_inboard_in - build_variables.dr_tf_inboard
            )
            < 10.0e0 * np.finfo(float(radius)).eps
        ):
            po.ocmmnt(self.outfile, "TF coil dimensions are consistent")
        else:
            po.ocmmnt(self.outfile, "ERROR: TF coil dimensions are NOT consistent:")
            po.ovarre(
                self.outfile,
                "Radius of plasma-facing side of inner leg SHOULD BE [m]",
                "",
                build_variables.r_tf_inboard_in + build_variables.dr_tf_inboard,
            )
            po.ovarre(
                self.outfile,
                "Inboard TF coil radial thickness [m]",
                "(dr_tf_inboard)",
                build_variables.dr_tf_inboard,
            )
            po.oblnkl(self.outfile)

        tf_total_height = (
            build_variables.dh_tf_inner_bore + 2 * build_variables.dr_tf_inboard
        )
        tf_total_width = (
            build_variables.dr_tf_inner_bore
            + build_variables.dr_tf_inboard
            + build_variables.dr_tf_outboard
        )
        po.oblnkl(self.outfile)
        po.obuild(
            self.outfile,
            "Total height and width of TFC [m]",
            tf_total_height,
            tf_total_width,
        )

        # Top section TF coil radial build (physics_variables.itart = 1 only)
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Radial build of TF coil at central collumn top :")
            # write(self.outfile,5)

            # Restart the radial build at bucking cylindre inner radius
            radius = build_variables.r_tf_inboard_in
            po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

            radius = radius + tfcoil_variables.dr_tf_nose_case
            po.obuild(
                self.outfile,
                "Coil bucking cylindre",
                tfcoil_variables.dr_tf_nose_case,
                radius,
                "(dr_tf_nose_case)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius + 0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - first half",
                0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf,
                radius,
                "(dr_tf_wp_top/2-tinstf)",
            )

            radius = (
                radius + 0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - second half",
                0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf,
                radius,
                "(dr_tf_wp_top/2-tinstf)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

            # Consistency check
            if abs(radius - build_variables.r_cp_top) < np.finfo(float(radius)).eps:
                po.ocmmnt(self.outfile, "Top TF coil dimensions are consistent")
            else:
                po.ocmmnt(self.outfile, "ERROR: TF coil dimensions are NOT consistent:")
                po.ovarre(
                    self.outfile,
                    "Radius of plasma-facing side of inner leg SHOULD BE [m]",
                    "",
                    build_variables.r_cp_top,
                )
                po.oblnkl(self.outfile)

    def tf_averaged_turn_geom(self, j_tf_wp, thwcndut, thicndut, i_tf_sc_mat):
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
            tfcoil_variables.cpttf = a_turn * j_tf_wp

        # Turn cable dimension is an input
        elif tfcoil_variables.t_cable_tf_is_input:
            # Turn squared dimension [m]
            tfcoil_variables.t_turn_tf = tfcoil_variables.t_cable_tf + 2.0e0 * (
                thicndut + thwcndut
            )

            # Turn area [m2]
            a_turn = tfcoil_variables.t_turn_tf**2

            # Current per turn [A]
            tfcoil_variables.cpttf = a_turn * j_tf_wp

        # Current per turn is an input
        else:
            # Turn area [m2]
            # Allow for additional inter-layer insulation MDK 13/11/18
            # Area of turn including conduit and inter-layer insulation
            a_turn = tfcoil_variables.cpttf / j_tf_wp

            # Dimension of square cross-section of each turn including inter-turn insulation [m]
            tfcoil_variables.t_turn_tf = np.sqrt(a_turn)

        # Square turn assumption
        sctfcoil_module.t_turn_radial = tfcoil_variables.t_turn_tf
        sctfcoil_module.t_turn_toroidal = tfcoil_variables.t_turn_tf

        # See derivation in the following document
        # k:\power plant physics and technology\process\hts\hts coil module for process.docx
        tfcoil_variables.t_conductor = (
            -tfcoil_variables.layer_ins
            + np.sqrt(tfcoil_variables.layer_ins**2 + 4.0e00 * a_turn)
        ) / 2 - 2.0e0 * thicndut

        # Total number of turns per TF coil (not required to be an integer)
        n_tf_turn = sctfcoil_module.awptf / a_turn

        # Area of inter-turn insulation: single turn [m2]
        insulation_area = a_turn - tfcoil_variables.t_conductor**2

        # NOTE: Fortran has acstf as an intent(out) variable that was outputting
        # into tfcoil_variables.acstf. The local variable, however, appears to
        # initially hold the value of tfcoil_variables.acstf despite not being
        # intent(in). I have replicated this behaviour in Python for now.
        acstf = copy.copy(tfcoil_variables.acstf)

        # ITER like turn structure
        if i_tf_sc_mat != 6:
            # Radius of rounded corners of cable space inside conduit [m]
            rbcndut = thwcndut * 0.75e0

            # Dimension of square cable space inside conduit [m]
            sctfcoil_module.t_cable = tfcoil_variables.t_conductor - 2.0e0 * thwcndut

            # Cross-sectional area of cable space per turn
            # taking account of rounded inside corners [m2]
            acstf = sctfcoil_module.t_cable**2 - (4.0e0 - np.pi) * rbcndut**2

            if acstf <= 0.0e0:
                if tfcoil_variables.t_conductor < 0.0e0:
                    error_handling.fdiags[0] = acstf
                    error_handling.fdiags[1] = sctfcoil_module.t_cable
                    error_handling.report_error(101)
                else:
                    error_handling.fdiags[0] = acstf
                    error_handling.fdiags[1] = sctfcoil_module.t_cable
                    error_handling.report_error(102)
                    rbcndut = 0.0e0
                    acstf = sctfcoil_module.t_cable**2

            # Cross-sectional area of conduit jacket per turn [m2]
            acndttf = tfcoil_variables.t_conductor**2 - acstf

        # REBCO turn structure
        elif i_tf_sc_mat == 6:
            # Diameter of circular cable space inside conduit [m]
            sctfcoil_module.t_cable = tfcoil_variables.t_conductor - 2.0e0 * thwcndut

            # Cross-sectional area of conduit jacket per turn [m2]
            acndttf = tfcoil_variables.t_conductor**2 - acstf

        return acstf, acndttf, insulation_area, n_tf_turn


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
    n_tf_turn: float,
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
    :param n_tf_turn: the number of turns per TF coil
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
        * n_tf_turn
        * i_op
        * (
            (np.exp(-lambda1 * tmaxforce) - np.exp(-lambda0 * tmaxforce))
            / (lambda0 - lambda1)
        )
    )
    i2 = (lambda1 / lambda2) * i1

    a_vv = (ro_vv + ri_vv) / (ro_vv - ri_vv)
    b_vvi = (constants.rmu0 * (n_tf_coils * n_tf_turn * i0 + i1 + (i2 / 2))) / (
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
    sctfcoil_module.tfc_current = 0.0
    sctfcoil_module.awpc = 0.0
    sctfcoil_module.awptf = 0.0
    sctfcoil_module.a_tf_steel = 0.0
    sctfcoil_module.a_tf_ins = 0.0
    sctfcoil_module.f_tf_steel = 0.0
    sctfcoil_module.f_tf_ins = 0.0
    sctfcoil_module.h_cp_top = 0.0
    sctfcoil_module.r_tf_outboard_in = 0.0
    sctfcoil_module.r_tf_outboard_out = 0.0
    sctfcoil_module.r_wp_inner = 0.0
    sctfcoil_module.r_wp_outer = 0.0
    sctfcoil_module.r_wp_centre = 0.0
    sctfcoil_module.dr_tf_wp_top = 0.0
    sctfcoil_module.vol_ins_cp = 0.0
    sctfcoil_module.vol_gr_ins_cp = 0.0
    sctfcoil_module.vol_case_cp = 0.0
    sctfcoil_module.t_wp_toroidal = 0.0
    sctfcoil_module.t_wp_toroidal_av = 0.0
    sctfcoil_module.t_lat_case_av = 0.0
    sctfcoil_module.a_case_front = 0.0
    sctfcoil_module.a_case_nose = 0.0
    sctfcoil_module.a_ground_ins = 0.0
    sctfcoil_module.a_leg_ins = 0.0
    sctfcoil_module.a_leg_gr_ins = 0.0
    sctfcoil_module.a_leg_cond = 0.0
    sctfcoil_module.theta_coil = 0.0
    sctfcoil_module.tan_theta_coil = 0.0
    sctfcoil_module.t_conductor_radial = 0.0
    sctfcoil_module.t_conductor_toroidal = 0.0
    sctfcoil_module.t_cable_radial = 0.0
    sctfcoil_module.t_cable_toroidal = 0.0
    sctfcoil_module.t_turn_radial = 0.0
    sctfcoil_module.t_turn_toroidal = 0.0
    sctfcoil_module.t_cable = 0.0
    sctfcoil_module.vforce_inboard_tot = 0.0
    sctfcoil_module.t1 = 0.0
    sctfcoil_module.time2 = 0.0
    sctfcoil_module.tau2 = 0.0
    sctfcoil_module.estotft = 0.0


def init_rebco_variables():
    """Initialise the REBCO variables"""
    rebco_variables.rebco_thickness = 1.0e-6
    rebco_variables.copper_thick = 100.0e-6
    rebco_variables.hastelloy_thickness = 50.0e-6
    rebco_variables.tape_width = 4.0e-3
    rebco_variables.croco_od = 0.0
    rebco_variables.croco_id = 0.0
    rebco_variables.croco_thick = 2.5e-3
    rebco_variables.copper_rrr = 100.0
    rebco_variables.coppera_m2_max = 1.0e8
    rebco_variables.f_coppera_m2 = 1.0
    rebco_variables.tape_thickness = 6.5e-5
    rebco_variables.stack_thickness = 0.0
    rebco_variables.tapes = 0.0
    rebco_variables.rebco_area = 0.0
    rebco_variables.copper_area = 0.0
    rebco_variables.hastelloy_area = 0.0
    rebco_variables.solder_area = 0.0
    rebco_variables.croco_area = 0.0
    rebco_variables.copperA_m2 = 0.0
    rebco_variables.copperaoh_m2_max = 1.0e8
    rebco_variables.f_copperaoh_m2 = 1.0
    rebco_variables.copperaoh_m2 = 0.0
