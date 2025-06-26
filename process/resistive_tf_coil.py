import logging

import numba
import numpy as np

from process.fortran import (
    build_variables,
    constants,
    error_handling,
    pfcoil_variables,
    physics_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.tf_coil import TFCoil

logger = logging.getLogger(__name__)

EPS = np.finfo(1.0).eps
RMU0 = constants.rmu0


class ResistiveTFCoil(TFCoil):
    def __init__(self):
        self.outfile = constants.nout

    def run(self, output: bool):
        """Run main tfcoil subroutine without outputting."""
        self.iprint = 0
        (
            sctfcoil_module.rad_tf_coil_toroidal,
            sctfcoil_module.tan_theta_coil,
            tfcoil_variables.a_tf_coil_inboard,
            sctfcoil_module.r_tf_outboard_in,
            sctfcoil_module.r_tf_outboard_out,
            tfcoil_variables.dx_tf_inboard_out_toroidal,
            tfcoil_variables.a_tf_leg_outboard,
            tfcoil_variables.dr_tf_plasma_case,
            tfcoil_variables.dx_tf_side_case,
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

        # Resistive coils : No approx necessary as the symmetry is cylindrical
        # The turn insulation th (tfcoil_variables.dx_tf_turn_insulation) is also subtracted too here
        tfcoil_variables.r_b_tf_inboard_peak = (
            build_variables.r_tf_inboard_out
            - tfcoil_variables.dr_tf_plasma_case
            - tfcoil_variables.dx_tf_turn_insulation
            - tfcoil_variables.dx_tf_wp_insulation
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
            a_tf_coil_inboard=tfcoil_variables.a_tf_coil_inboard,
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
        self.res_tf_internal_geom()
        self.tf_res_heating()

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
        tfcoil_variables.estotftgj = 1.0e-9 * sctfcoil_module.e_tf_magnetic_stored_total

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
                sctfcoil_module.r_tf_wp_inner,
                sctfcoil_module.tan_theta_coil,
                sctfcoil_module.rad_tf_coil_toroidal,
                sctfcoil_module.r_tf_wp_outer,
                sctfcoil_module.a_tf_coil_inboard_steel,
                sctfcoil_module.a_case_front,
                sctfcoil_module.a_tf_coil_nose_case,
                tfcoil_variables.dx_tf_wp_insertion_gap,
                tfcoil_variables.dx_tf_wp_insulation,
                tfcoil_variables.n_tf_coil_turns,
                int(tfcoil_variables.i_tf_turns_integer),
                sctfcoil_module.t_cable,
                sctfcoil_module.dr_tf_turn_cable_space,
                tfcoil_variables.dia_tf_turn_coolant_channel,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.dx_tf_turn_steel,
                sctfcoil_module.t_lat_case_av,
                sctfcoil_module.t_wp_toroidal_av,
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
        if output:
            self.outtf(0)

    def res_tf_internal_geom(self):
        """
        Author : S. Kahn
        Resisitve TF turn geometry, equivalent to winding_pack subroutines

        """
        sctfcoil_module.r_tf_wp_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.dr_tf_nose_case
        )
        sctfcoil_module.r_tf_wp_outer = (
            build_variables.r_tf_inboard_out - tfcoil_variables.dr_tf_plasma_case
        )

        # Conductor layer radial thickness at centercollumn top [m]
        if physics_variables.itart == 1:
            sctfcoil_module.dr_tf_wp_top = (
                build_variables.r_cp_top
                - tfcoil_variables.dr_tf_plasma_case
                - tfcoil_variables.dr_tf_nose_case
                - build_variables.r_tf_inboard_in
            )

        # Number of turns
        # Set by user (no turn structure by default, i.e. tfcoil_variables.n_tf_coil_turns = 1 )
        if (
            abs(tfcoil_variables.n_tf_coil_turns)
            < np.finfo(float(tfcoil_variables.n_tf_coil_turns)).eps
        ):
            tfcoil_variables.n_tf_coil_turns = 1.0e0

        # Total mid-plane cross-sectional area of winding pack, [m2]
        # including the surrounding ground-wall insulation layer
        sctfcoil_module.a_tf_wp_with_insulation = (
            np.pi
            * (sctfcoil_module.r_tf_wp_outer**2 - sctfcoil_module.r_tf_wp_inner**2)
            / tfcoil_variables.n_tf_coils
        )

        # Area of the front case, the plasma-facing case of the inner TF coil [m2]
        sctfcoil_module.a_case_front = (
            np.pi
            * (
                (sctfcoil_module.r_tf_wp_outer + tfcoil_variables.dr_tf_plasma_case)
                ** 2
                - sctfcoil_module.r_tf_wp_outer**2
            )
            / tfcoil_variables.n_tf_coils
        )

        # WP mid-plane cross-section excluding ground insulation per coil [m2]
        sctfcoil_module.a_tf_wp_no_insulation = (
            np.pi
            * (
                (sctfcoil_module.r_tf_wp_outer - tfcoil_variables.dx_tf_wp_insulation)
                ** 2
                - (sctfcoil_module.r_tf_wp_inner + tfcoil_variables.dx_tf_wp_insulation)
                ** 2
            )
            / tfcoil_variables.n_tf_coils
            - 2.0e0
            * tfcoil_variables.dx_tf_wp_insulation
            * (
                tfcoil_variables.dr_tf_wp_with_insulation
                - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
            )
        )

        # Ground insulation cross-section area per coil [m2]
        sctfcoil_module.a_tf_wp_ground_insulation = (
            sctfcoil_module.a_tf_wp_with_insulation
            - sctfcoil_module.a_tf_wp_no_insulation
        )

        # Exact mid-plane cross-section area of the conductor per TF coil [m2]
        a_tf_cond = np.pi * (
            (
                sctfcoil_module.r_tf_wp_outer
                - tfcoil_variables.dx_tf_wp_insulation
                - tfcoil_variables.dx_tf_turn_insulation
            )
            ** 2
            - (
                sctfcoil_module.r_tf_wp_inner
                + tfcoil_variables.dx_tf_wp_insulation
                + tfcoil_variables.dx_tf_turn_insulation
            )
            ** 2
        ) / tfcoil_variables.n_tf_coils - (
            tfcoil_variables.dr_tf_wp_with_insulation
            - 2.0e0
            * (
                tfcoil_variables.dx_tf_wp_insulation
                + tfcoil_variables.dx_tf_turn_insulation
            )
        ) * 2.0e0 * (
            tfcoil_variables.dx_tf_wp_insulation
            + tfcoil_variables.dx_tf_turn_insulation * tfcoil_variables.n_tf_coil_turns
        )
        a_tf_cond = a_tf_cond * (1.0e0 - tfcoil_variables.fcoolcp)

        # Inter turn insulation area per coil [m2]
        tfcoil_variables.a_tf_coil_wp_turn_insulation = (
            sctfcoil_module.a_tf_wp_no_insulation
            - a_tf_cond / (1.0e0 - tfcoil_variables.fcoolcp)
        )

        # Total insulation cross-section per coil [m2]
        sctfcoil_module.a_tf_coil_inboard_insulation = (
            tfcoil_variables.a_tf_coil_wp_turn_insulation
            + sctfcoil_module.a_tf_wp_ground_insulation
        )

        # Insulation fraction [-]
        sctfcoil_module.f_tf_ins = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_coil_inboard_insulation
            / tfcoil_variables.a_tf_coil_inboard
        )

        # Total cross-sectional area of the bucking cylindre and the outer support
        # support structure per coil [m2]
        # physics_variables.itart = 1 : Only valid at mid-plane
        tfcoil_variables.a_tf_coil_inboard_case = (
            tfcoil_variables.a_tf_coil_inboard / tfcoil_variables.n_tf_coils
        ) - sctfcoil_module.a_tf_wp_with_insulation

        # Current per turn
        tfcoil_variables.c_tf_turn = tfcoil_variables.c_tf_total / (
            tfcoil_variables.n_tf_coil_turns * tfcoil_variables.n_tf_coils
        )

        # Exact current density on TF oubard legs
        tfcoil_variables.cdtfleg = tfcoil_variables.c_tf_total / (
            (1.0e0 - tfcoil_variables.fcoolcp)
            * (
                tfcoil_variables.dx_tf_inboard_out_toroidal
                - 2.0e0
                * (
                    tfcoil_variables.n_tf_coil_turns
                    * tfcoil_variables.dx_tf_turn_insulation
                    + tfcoil_variables.dx_tf_wp_insulation
                )
            )
            * (
                build_variables.dr_tf_outboard
                - 2.0e0
                * (
                    tfcoil_variables.dx_tf_turn_insulation
                    + tfcoil_variables.dx_tf_wp_insulation
                )
            )
        )

        # Reporting negative WP areas issues
        if sctfcoil_module.a_tf_wp_with_insulation < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.a_tf_wp_with_insulation
            error_handling.fdiags[0] = tfcoil_variables.dr_tf_wp_with_insulation
            error_handling.report_error(99)

        elif sctfcoil_module.a_tf_wp_no_insulation < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.a_tf_wp_no_insulation
            error_handling.report_error(101)

    def tf_res_heating(self) -> None:
        """
        Calculate resistive heating for resistive magnets.

        This method calculates the resistive heating for resistive magnets.
        It considers the following scenarios:
        - Clamped joints in superconductors might have resistive power losses on the joints.
        - Sliding joints might have a region of high resistivity.

        Notes:
        - The copper resisitivty is set to be that for GLIDCOP AL-15 at 20°C for copper (i_tf_sup = 0).
        - The coefficient of resistivity is set to be that of pure copper

        References:
            - https://www.spotweldingconsultants.com/GlidCop_AL_15.pdf

            - https://cirris.com/temperature-coefficient-of-copper/
        """

        # Resistivity of the Glidcop copper centerpost
        if tfcoil_variables.i_tf_sup == 0:
            tfcoil_variables.rho_cp = (
                # 1.86 is the resitivity at `20°C` for GLIDCOP AL-15
                # 0.00393 is the coefficient of resistivity for copper
                tfcoil_variables.frhocp
                * (1.86e0 + 0.00393e0 * (tfcoil_variables.temp_cp_average - 293.15e0))
                * 1.0e-8
            )

        # Resistivity of the aluminium centerpost
        if tfcoil_variables.i_tf_sup == 2:
            tfcoil_variables.rho_cp = tfcoil_variables.frhocp * (
                2.00016e-14 * tfcoil_variables.temp_cp_average**3
                - 6.75384e-13 * tfcoil_variables.temp_cp_average**2
                + 8.89159e-12 * tfcoil_variables.temp_cp_average
            )

        # Calculations dedicated for configurations with CP
        if physics_variables.itart == 1:
            # Tricky trick to make the leg / CP tempearture the same
            if (
                abs(tfcoil_variables.temp_tf_legs_outboard + 1.0e0)
                < np.finfo(float(tfcoil_variables.temp_tf_legs_outboard)).eps
            ):
                sctfcoil_module.is_leg_cp_temp_same = 1
                tfcoil_variables.temp_tf_legs_outboard = (
                    tfcoil_variables.temp_cp_average
                )

            # Leg resistivity (different leg temperature as separate cooling channels)
            if tfcoil_variables.i_tf_sup == 0:
                tfcoil_variables.rho_tf_leg = (
                    tfcoil_variables.frholeg
                    * (
                        1.86e0
                        + 0.00393e0
                        * (tfcoil_variables.temp_tf_legs_outboard - 293.15e0)
                    )
                    * 1.0e-8
                )
            elif tfcoil_variables.i_tf_sup == 2:
                tfcoil_variables.rho_tf_leg = tfcoil_variables.frholeg * (
                    2.00016e-14 * tfcoil_variables.temp_tf_legs_outboard**3
                    - 6.75384e-13 * tfcoil_variables.temp_tf_legs_outboard**2
                    + 8.89159e-12 * tfcoil_variables.temp_tf_legs_outboard
                )

            # Tricky trick to make the leg / CP tempearture the same
            if sctfcoil_module.is_leg_cp_temp_same == 1:
                tfcoil_variables.temp_tf_legs_outboard = -1.0e0

            # Centrepost resisitivity and conductor/insulation volume

            (
                tfcoil_variables.a_cp_cool,
                tfcoil_variables.vol_cond_cp,
                tfcoil_variables.p_cp_resistive,
                sctfcoil_module.vol_ins_cp,
                sctfcoil_module.vol_case_cp,
                sctfcoil_module.vol_gr_ins_cp,
            ) = self.cpost(
                build_variables.r_tf_inboard_in,
                build_variables.r_tf_inboard_out,
                build_variables.r_cp_top,
                sctfcoil_module.z_cp_top,
                build_variables.z_tf_inside_half + build_variables.dr_tf_outboard,
                tfcoil_variables.dr_tf_nose_case,
                tfcoil_variables.dr_tf_plasma_case,
                tfcoil_variables.dx_tf_wp_insulation,
                tfcoil_variables.dx_tf_turn_insulation,
                tfcoil_variables.n_tf_coil_turns,
                tfcoil_variables.c_tf_total,
                tfcoil_variables.rho_cp,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf_coils,
            )

        # Leg cross-section areas
        # Rem : For physics_variables.itart = 1, these quantitire corresponds to the outer leg only
        # ---
        # Leg ground insulation area per coil [m2]
        sctfcoil_module.a_leg_gr_ins = tfcoil_variables.a_tf_leg_outboard - (
            tfcoil_variables.dx_tf_inboard_out_toroidal
            - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        ) * (
            build_variables.dr_tf_outboard
            - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        )

        # Outboard leg turns insulation area per coil [m2]
        sctfcoil_module.a_leg_ins = (
            2.0e0
            * tfcoil_variables.dx_tf_turn_insulation
            * (
                tfcoil_variables.dx_tf_inboard_out_toroidal
                - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
            )
            + 2.0e0
            * tfcoil_variables.dx_tf_turn_insulation
            * tfcoil_variables.n_tf_coil_turns
            * (
                build_variables.dr_tf_outboard
                - 2.0e0
                * (
                    tfcoil_variables.dx_tf_turn_insulation
                    + tfcoil_variables.dx_tf_wp_insulation
                )
            )
        )  # toroidal direction + radial direction

        # Exact TF outboard leg conductor area per coil [m2]
        sctfcoil_module.a_leg_cond = (1.0e0 - tfcoil_variables.f_a_tf_cool_outboard) * (
            tfcoil_variables.a_tf_leg_outboard
            - sctfcoil_module.a_leg_gr_ins
            - sctfcoil_module.a_leg_ins
        )
        # ---

        if physics_variables.itart == 1:
            # Outer leg resistive power loss
            # ---
            # TF outboard leg's resistance calculation (per leg) [ohm]
            tfcoil_variables.res_tf_leg = (
                tfcoil_variables.rho_tf_leg
                * tfcoil_variables.len_tf_coil
                / sctfcoil_module.a_leg_cond
            )

            # TF outer leg resistive power (TOTAL) [W]
            tfcoil_variables.p_tf_leg_resistive = (
                tfcoil_variables.res_tf_leg
                * (tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils) ** 2
            ) * tfcoil_variables.n_tf_coils
            # ---

            # Sliding joints resistive heating
            # ---
            if tfcoil_variables.i_cp_joints != 0:
                # Number of contact area per joint (all legs)
                n_contact_tot = (
                    tfcoil_variables.n_tf_joints_contact
                    * np.round(tfcoil_variables.n_tf_coil_turns)
                    * np.round(tfcoil_variables.n_tf_coils)
                )

                # Area of joint contact (all legs)
                a_joints = (
                    build_variables.dr_tf_outboard
                    * tfcoil_variables.th_joint_contact
                    * n_contact_tot
                )

                # Total joints resistive power losses
                tfcoil_variables.p_tf_joints_resistive = (
                    tfcoil_variables.n_tf_joints
                    * tfcoil_variables.rho_tf_joints
                    * tfcoil_variables.c_tf_total**2
                    / a_joints
                )
            else:
                # Joints resistance to be evaluated for SC
                tfcoil_variables.p_tf_joints_resistive = 0.0e0

            # ---

        # Case of a resistive magnet without joints
        # ***
        else:
            # TF resistive powers
            tfcoil_variables.p_cp_resistive = (
                tfcoil_variables.rho_cp
                * tfcoil_variables.c_tf_total**2
                * tfcoil_variables.len_tf_coil
                / (sctfcoil_module.a_leg_cond * tfcoil_variables.n_tf_coils)
            )

            # tfcoil_variables.p_cp_resistive containts the the total resistive power losses
            tfcoil_variables.p_tf_leg_resistive = 0.0e0

            # No joints if physics_variables.itart = 0
            tfcoil_variables.p_tf_joints_resistive = 0.0e0

    @staticmethod
    @numba.njit(cache=True)
    def cpost(
        r_tf_inboard_in,
        r_tf_inboard_out,
        r_cp_top,
        ztop,
        hmaxi,
        cas_in_th,
        cas_out_th,
        gr_ins_th,
        ins_th,
        n_tf_coil_turns,
        curr,
        rho,
        fcool,
        n_tf_coils,
    ):
        """
        author: P J Knight, CCFE, Culham Science Centre
        Calculates the volume and resistive power losses of a TART centrepost
        This routine calculates the volume and resistive power losses
        of a TART centrepost. It is assumed to be tapered - narrowest at
        the midplane and reaching maximum thickness at the height of the
        plasma. Above/below the plasma, the centrepost is cylindrical.
        The shape of the taper is assumed to be an arc of a circle.
        P J Knight, CCFE, Culham Science Centre
        21/10/96 PJK Initial version
        08/05/12 PJK Initial F90 version
        16/10/12 PJK Added constants; removed argument pi
        26/06/14 PJK Added error handling
        12/11/19 SK Using fixed cooling cross-section area along the CP
        26/11/19 SK added the coolant area, the conuctor/isulator/outer casing volume
        30/11/20 SK added the ground outer ground insulation volume
        F/MI/PJK/LOGBOOK12, pp.33,34
        """
        yy_ins = np.zeros((101,))  # Exact conductor area (to be integrated)
        yy_cond = np.zeros((101,))  # Turn insulation area (to be integrated)
        yy_gr_ins = np.zeros((101,))  # Outter ground insulation area (to be integrated)
        yy_casout = np.zeros((101,))  # Outter case area (to be integrated)

        rtop = r_cp_top - cas_out_th - gr_ins_th

        # Conductor outer radius at CP mid-plane [m]
        rmid = r_tf_inboard_out - cas_out_th - gr_ins_th

        # Conductor inner radius [m]
        r_tfin_inleg = r_tf_inboard_in + cas_in_th + gr_ins_th
        # -#

        #  Error traps
        # ------------
        # if rtop <= 0.0e0:
        #     error_handling.fdiags[0] = rtop
        #     error_handling.report_error(115)

        # if ztop <= 0.0e0:
        #     error_handling.fdiags[0] = ztop
        #     error_handling.report_error(116)

        # if rmid <= 0.0e0:
        #     error_handling.fdiags[0] = rmid
        #     error_handling.report_error(117)

        # if build_variables.z_tf_inside_half <= 0.0e0:
        #     error_handling.fdiags[0] = build_variables.z_tf_inside_half
        #     error_handling.report_error(118)

        # if (fcool < 0.0e0) or (fcool > 1.0e0):
        #     error_handling.fdiags[0] = fcool
        #     error_handling.report_error(119)

        # if rtop < rmid:
        #     error_handling.fdiags[0] = rtop
        #     error_handling.fdiags[1] = rmid
        #     error_handling.report_error(120)

        # if build_variables.z_tf_inside_half < ztop:
        #     error_handling.fdiags[0] = build_variables.z_tf_inside_half
        #     error_handling.fdiags[1] = ztop
        #     error_handling.report_error(121)

        # ------------

        # Mid-plane area calculations
        # ---------------------------
        # Total number of CP turns
        n_turns_tot = n_tf_coils * n_tf_coil_turns

        # Area of the innner TF central hole [m2]
        a_tfin_hole = np.pi * r_tfin_inleg**2

        # Mid-plane outer casing cross-section area [m2]
        a_casout = np.pi * (
            (rmid + gr_ins_th + cas_out_th) ** 2 - (rmid + gr_ins_th) ** 2
        )

        # Mid-plane outter ground insulation thickness [m2]
        a_cp_gr_ins = (
            np.pi * ((rmid + gr_ins_th) ** 2 - rmid**2)
            + 2.0e0 * gr_ins_th * (rmid - r_tfin_inleg) * n_tf_coils
        )

        # Mid-plane turn layer cross-section area [m2]
        a_cp_ins = (
            np.pi
            * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)  # Inner layer volume
            + np.pi * (rmid**2 - (rmid - ins_th) ** 2)  # Outter layer volume
            + 2.0e0 * n_turns_tot * ins_th * (rmid - r_tfin_inleg - 2.0e0 * ins_th)
        )  # inter turn separtion

        # Cooling pipes cross-section per coil [m2]
        a_cp_cool = fcool * (
            (np.pi * rmid**2 - a_tfin_hole - a_cp_ins) / n_tf_coils
            - 2.0e0 * gr_ins_th * (rmid - r_tfin_inleg)
        )  # Wedge ground insulation
        # ---------------------------

        #  Trivial solutions
        # ------------------
        if np.abs(fcool) < EPS:
            vol_cond_cp = 0.0e0
            respow = 0.0e0
            vol_case_cp = 0.0e0
            vol_gr_ins_cp = 0.0e0
            vol_ins_cp = 0.0e0
            # error_handling.report_error(122)
            return (
                a_cp_cool,
                vol_cond_cp,
                respow,
                vol_ins_cp,
                vol_case_cp,
                vol_gr_ins_cp,
            )

        if np.abs(rmid - rtop) < EPS:
            # Exact conductor cross-section
            a_cond_midplane = (
                np.pi * rmid**2 - a_tfin_hole - n_tf_coils * a_cp_cool - a_cp_ins
            )

            # Volumes and resisitive losses calculations
            vol_cond_cp = 2.0e0 * hmaxi * a_cond_midplane
            vol_ins_cp = 2.0e0 * hmaxi * a_cp_ins
            vol_gr_ins_cp = 2.0e0 * hmaxi * a_cp_gr_ins
            respow = 2.0e0 * hmaxi * curr**2 * rho / a_cond_midplane
            vol_case_cp = 2.0e0 * hmaxi * a_casout

            return (
                a_cp_cool,
                vol_cond_cp,
                respow,
                vol_ins_cp,
                vol_case_cp,
                vol_gr_ins_cp,
            )

        # ------------------

        # Find centre of circle (RC,0) defining the taper's arc
        # (r1,z1) is midpoint of line joining (rmid,0) and (rtop,ztop)
        # Rem : The taper arc is defined using the outer radius of the
        #       conductor including turn unsulation
        # -------------------------------------------------------------
        r1 = 0.5e0 * (rmid + rtop)
        z1 = 0.5e0 * ztop

        x = (r1 - rmid) ** 2 + z1**2
        y = ztop**2 / ((rtop - rmid) ** 2 + ztop**2)

        rc = rmid + np.sqrt(x / (1.0e0 - y))
        # -------------------------------------------------------------

        #  Find volume of tapered section of centrepost, and the resistive
        #  power losses, by integrating along the centrepost from the midplane
        # --------------------------------------------------------------------
        #  Calculate centrepost radius and cross-sectional areas at each Z
        dz = 0.01e0 * ztop

        for ii in range(101):
            z = ii * dz
            z = np.fmin(np.array(z), ztop)

            r = rc - np.sqrt((rc - rmid) ** 2 - z * z)

            # if r <= 0.0e0:
            #     error_handling.fdiags[0] = r
            #     error_handling.fdiags[1] = rc
            #     error_handling.fdiags[2] = rmid
            #     error_handling.fdiags[3] = z

            #     error_handling.report_error(123)

            # Insulation cross-sectional area at z
            yy_ins[ii] = (
                np.pi * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)
                + np.pi * (r**2 - (r - ins_th) ** 2)  # Inner layer volume
                + 2.0e0  # Outter layer volume
                * ins_th
                * (r - r_tfin_inleg - 2.0e0 * ins_th)
                * n_turns_tot
            )  # inter turn layers

            #  Conductor cross-sectional area at z
            yy_cond[ii] = (
                np.pi * r**2
                - a_tfin_hole
                - n_tf_coils * a_cp_cool
                - yy_ins[ii]
                - 2.0e0 * n_tf_coils * gr_ins_th * (r - r_tfin_inleg)
            )  # Wedge ground insulation

            #  Outer ground insulation area at z
            yy_gr_ins[ii] = np.pi * (
                (r + gr_ins_th) ** 2 - r**2
            ) + 2.0e0 * n_tf_coils * gr_ins_th * (r - r_tfin_inleg)

            #  Outer casing Cross-sectional area at z
            yy_casout[ii] = np.pi * (
                (r + gr_ins_th + cas_out_th) ** 2 - (r + gr_ins_th) ** 2
            )

        #  Perform integrals using trapezium rule
        sum1 = 0.0e0
        sum2 = 0.0e0
        sum3 = 0.0e0
        sum4 = 0.0e0
        sum5 = 0.0e0
        for ii in range(1, 100):
            sum1 = sum1 + yy_cond[ii]
            sum2 = sum2 + 1.0e0 / yy_cond[ii]
            sum3 = sum3 + yy_ins[ii]
            sum4 = sum4 + yy_casout[ii]
            sum5 = sum5 + yy_gr_ins[ii]

        sum1 = 0.5e0 * dz * (yy_cond[0] + yy_cond[100] + 2.0e0 * sum1)
        sum2 = 0.5e0 * dz * (1.0e0 / yy_cond[0] + 1.0e0 / yy_cond[100] + 2.0e0 * sum2)
        sum3 = 0.5e0 * dz * (yy_ins[0] + yy_ins[100] + 2.0e0 * sum3)
        sum4 = 0.5e0 * dz * (yy_casout[0] + yy_casout[100] + 2.0e0 * sum4)
        sum5 = 0.5e0 * dz * (yy_gr_ins[0] + yy_gr_ins[100] + 2.0e0 * sum5)

        # Turn insulation layer cross section at CP top  [m2]
        a_cp_ins = (
            np.pi * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)
            + np.pi * (rtop**2 - (rtop - ins_th) ** 2)  # Inner layer volume
            + 2.0e0  # Outter layer volume
            * ins_th
            * (rtop - r_tfin_inleg - 2.0e0 * ins_th)
            * n_turns_tot
        )  # turn separtion layers

        # Ground insulation layer cross-section at CP top [m2]
        a_cp_gr_ins = (
            np.pi * ((rtop + gr_ins_th) ** 2 - rtop**2)
            + 2.0e0 * gr_ins_th * (rtop - r_tfin_inleg) * n_tf_coils
        )

        # Outer casing cross-section area at CP top [m2]
        a_casout = np.pi * (
            (rmid + gr_ins_th + cas_out_th) ** 2 - (rmid + gr_ins_th) ** 2
        )

        # Centrepost volume (ignoring coolant fraction) [m3]
        vol_cond_cp = 2.0e0 * sum1 + 2.0e0 * (  # Tapered section
            hmaxi - ztop
        ) * (  # Straight section vertical height
            np.pi * rtop**2
            - a_tfin_hole
            - a_cp_ins
            - n_tf_coils * a_cp_cool
            - 2.0e0 * n_tf_coils * gr_ins_th * (rtop - r_tfin_inleg)
        )  # subtracting ground insulation wedge separation

        # Resistive power losses in taped section (variable radius section) [W]
        res_taped = rho * curr**2 * sum2

        # Centrepost insulator volume [m3]
        vol_ins_cp = 2.0e0 * (sum3 + (hmaxi - ztop) * a_cp_ins)

        # Ground insulation volume [m3]
        vol_gr_ins_cp = 2.0e0 * (
            sum5
            + (hmaxi - ztop) * a_cp_gr_ins
            + hmaxi * np.pi * (r_tfin_inleg**2 - (r_tfin_inleg - gr_ins_th) ** 2)
        )

        # CP casing volume [m3]
        vol_case_cp = 2.0e0 * (
            sum4
            + (hmaxi - ztop) * a_casout
            + hmaxi
            * np.pi
            * (
                (r_tfin_inleg - gr_ins_th) ** 2
                - (r_tfin_inleg - gr_ins_th - cas_in_th) ** 2
            )
        )

        # Resistive power losses in cylindrical section (constant radius) [W]
        res_cyl = (
            rho
            * curr**2
            * (
                (hmaxi - ztop)
                / (
                    np.pi * rtop**2
                    - a_tfin_hole
                    - a_cp_ins
                    - n_tf_coils * a_cp_cool
                    - 2.0e0 * n_tf_coils * gr_ins_th * (rtop - r_tfin_inleg)
                )
            )
        )  # ground insulation separation

        # Total CP resistive power [W]
        respow = 2.0e0 * (res_cyl + res_taped)

        return (
            a_cp_cool,
            vol_cond_cp,
            respow,
            vol_ins_cp,
            vol_case_cp,
            vol_gr_ins_cp,
        )
