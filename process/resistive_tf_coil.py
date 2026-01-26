import logging

import numba
import numpy as np

from process import constants
from process.data_structure import (
    build_variables,
    fwbs_variables,
    pfcoil_variables,
    physics_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
)
from process.tf_coil import TFCoil

logger = logging.getLogger(__name__)

EPS = np.finfo(1.0).eps


class ResistiveTFCoil(TFCoil):
    def __init__(self):
        self.outfile = constants.NOUT

    def run(self, output: bool):
        """Run main tfcoil subroutine without outputting."""
        self.iprint = 0

        # Set up TF values share by all coil types
        self.run_base_tf()

        self.res_tf_internal_geom()
        self.tf_res_heating()

        tfcoil_variables.ind_tf_coil = self.tf_coil_self_inductance(
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

        (
            superconducting_tf_coil_variables.e_tf_magnetic_stored_total,
            tfcoil_variables.e_tf_magnetic_stored_total_gj,
            tfcoil_variables.e_tf_coil_magnetic_stored,
        ) = self.tf_stored_magnetic_energy(
            ind_tf_coil=tfcoil_variables.ind_tf_coil,
            c_tf_total=tfcoil_variables.c_tf_total,
            n_tf_coils=tfcoil_variables.n_tf_coils,
        )

        (
            tfcoil_variables.cforce,
            tfcoil_variables.vforce,
            tfcoil_variables.vforce_outboard,
            superconducting_tf_coil_variables.vforce_inboard_tot,
            tfcoil_variables.f_vforce_inboard,
        ) = self.tf_field_and_force(
            i_tf_sup=tfcoil_variables.i_tf_sup,
            r_tf_wp_inboard_outer=superconducting_tf_coil_variables.r_tf_wp_inboard_outer,
            r_tf_wp_inboard_inner=superconducting_tf_coil_variables.r_tf_wp_inboard_inner,
            r_tf_outboard_in=superconducting_tf_coil_variables.r_tf_outboard_in,
            dx_tf_wp_insulation=tfcoil_variables.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=tfcoil_variables.dx_tf_wp_insertion_gap,
            b_tf_inboard_peak_symmetric=tfcoil_variables.b_tf_inboard_peak_symmetric,
            c_tf_total=tfcoil_variables.c_tf_total,
            n_tf_coils=tfcoil_variables.n_tf_coils,
            dr_tf_plasma_case=tfcoil_variables.dr_tf_plasma_case,
            rmajor=physics_variables.rmajor,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            r_cp_top=build_variables.r_cp_top,
            itart=physics_variables.itart,
            i_cp_joints=tfcoil_variables.i_cp_joints,
            f_vforce_inboard=tfcoil_variables.f_vforce_inboard,
        )

        # Calculate TF coil areas and masses
        self.generic_tf_coil_area_and_masses()
        self.resistive_tf_coil_areas_and_masses()

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
                int(tfcoil_variables.n_tf_wp_stress_layers),
                int(tfcoil_variables.i_tf_bucking),
                float(build_variables.r_tf_inboard_in),
                build_variables.dr_bore,
                build_variables.dr_cs,
                build_variables.i_tf_inside_cs,
                build_variables.dr_tf_inboard,
                build_variables.dr_cs_tf_gap,
                pfcoil_variables.i_pf_conductor,
                pfcoil_variables.j_cs_flat_top_end,
                pfcoil_variables.j_cs_pulse_start,
                pfcoil_variables.c_pf_coil_turn_peak_input,
                pfcoil_variables.n_pf_coils_in_group,
                pfcoil_variables.f_dr_dz_cs_turn,
                pfcoil_variables.radius_cs_turn_corners,
                pfcoil_variables.f_a_cs_turn_steel,
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
                superconducting_tf_coil_variables.r_tf_wp_inboard_inner,
                superconducting_tf_coil_variables.tan_theta_coil,
                superconducting_tf_coil_variables.rad_tf_coil_inboard_toroidal_half,
                superconducting_tf_coil_variables.r_tf_wp_inboard_outer,
                superconducting_tf_coil_variables.a_tf_coil_inboard_steel,
                superconducting_tf_coil_variables.a_tf_plasma_case,
                superconducting_tf_coil_variables.a_tf_coil_nose_case,
                tfcoil_variables.dx_tf_wp_insertion_gap,
                tfcoil_variables.dx_tf_wp_insulation,
                tfcoil_variables.n_tf_coil_turns,
                int(tfcoil_variables.i_tf_turns_integer),
                superconducting_tf_coil_variables.dx_tf_turn_cable_space_average,
                superconducting_tf_coil_variables.dr_tf_turn_cable_space,
                tfcoil_variables.dia_tf_turn_coolant_channel,
                tfcoil_variables.f_a_tf_turn_cable_copper,
                tfcoil_variables.dx_tf_turn_steel,
                superconducting_tf_coil_variables.dx_tf_side_case_average,
                superconducting_tf_coil_variables.dx_tf_wp_toroidal_average,
                superconducting_tf_coil_variables.a_tf_coil_inboard_insulation,
                tfcoil_variables.a_tf_wp_steel,
                tfcoil_variables.a_tf_wp_conductor,
                superconducting_tf_coil_variables.a_tf_wp_with_insulation,
                tfcoil_variables.eyoung_al,
                tfcoil_variables.poisson_al,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf_graded_layers,
                tfcoil_variables.c_tf_total,
                tfcoil_variables.dr_tf_plasma_case,
                tfcoil_variables.i_tf_stress_model,
                superconducting_tf_coil_variables.vforce_inboard_tot,
                tfcoil_variables.i_tf_tresca,
                tfcoil_variables.a_tf_coil_inboard_case,
                tfcoil_variables.vforce,
                tfcoil_variables.a_tf_turn_steel,
                pfcoil_variables.a_cs_poloidal,
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
                tfcoil_variables.casestr if tfcoil_variables.casestr is None else casestr
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
                logger.error(
                    "Invalid stress model (r_tf_inboard = 0), stress constraint switched off"
                )
                tfcoil_variables.sig_tf_case = 0.0e0
                tfcoil_variables.sig_tf_wp = 0.0e0
        if output:
            self.outtf()

    def res_tf_internal_geom(self):
        """
        Author : S. Kahn
        Resisitve TF turn geometry, equivalent to winding_pack subroutines

        """
        superconducting_tf_coil_variables.r_tf_wp_inboard_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.dr_tf_nose_case
        )
        superconducting_tf_coil_variables.r_tf_wp_inboard_outer = (
            build_variables.r_tf_inboard_out - tfcoil_variables.dr_tf_plasma_case
        )

        # Conductor layer radial thickness at centercollumn top [m]
        if physics_variables.itart == 1:
            superconducting_tf_coil_variables.dr_tf_wp_top = (
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
        superconducting_tf_coil_variables.a_tf_wp_with_insulation = (
            np.pi
            * (
                superconducting_tf_coil_variables.r_tf_wp_inboard_outer**2
                - superconducting_tf_coil_variables.r_tf_wp_inboard_inner**2
            )
            / tfcoil_variables.n_tf_coils
        )

        # Area of the front case, the plasma-facing case of the inner TF coil [m2]
        superconducting_tf_coil_variables.a_tf_plasma_case = (
            np.pi
            * (
                (
                    superconducting_tf_coil_variables.r_tf_wp_inboard_outer
                    + tfcoil_variables.dr_tf_plasma_case
                )
                ** 2
                - superconducting_tf_coil_variables.r_tf_wp_inboard_outer**2
            )
            / tfcoil_variables.n_tf_coils
        )

        # WP mid-plane cross-section excluding ground insulation per coil [m2]
        superconducting_tf_coil_variables.a_tf_wp_no_insulation = (
            np.pi
            * (
                (
                    superconducting_tf_coil_variables.r_tf_wp_inboard_outer
                    - tfcoil_variables.dx_tf_wp_insulation
                )
                ** 2
                - (
                    superconducting_tf_coil_variables.r_tf_wp_inboard_inner
                    + tfcoil_variables.dx_tf_wp_insulation
                )
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
        superconducting_tf_coil_variables.a_tf_wp_ground_insulation = (
            superconducting_tf_coil_variables.a_tf_wp_with_insulation
            - superconducting_tf_coil_variables.a_tf_wp_no_insulation
        )

        # Exact mid-plane cross-section area of the conductor per TF coil [m2]
        tfcoil_variables.a_res_tf_coil_conductor = np.pi * (
            (
                superconducting_tf_coil_variables.r_tf_wp_inboard_outer
                - tfcoil_variables.dx_tf_wp_insulation
                - tfcoil_variables.dx_tf_turn_insulation
            )
            ** 2
            - (
                superconducting_tf_coil_variables.r_tf_wp_inboard_inner
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
        tfcoil_variables.a_res_tf_coil_conductor = (
            tfcoil_variables.a_res_tf_coil_conductor * (1.0e0 - tfcoil_variables.fcoolcp)
        )

        # Inter turn insulation area per coil [m2]
        tfcoil_variables.a_tf_coil_wp_turn_insulation = (
            superconducting_tf_coil_variables.a_tf_wp_no_insulation
            - tfcoil_variables.a_res_tf_coil_conductor
            / (1.0e0 - tfcoil_variables.fcoolcp)
        )

        # Total insulation cross-section per coil [m2]
        superconducting_tf_coil_variables.a_tf_coil_inboard_insulation = (
            tfcoil_variables.a_tf_coil_wp_turn_insulation
            + superconducting_tf_coil_variables.a_tf_wp_ground_insulation
        )

        # Insulation fraction [-]
        superconducting_tf_coil_variables.f_a_tf_coil_inboard_insulation = (
            tfcoil_variables.n_tf_coils
            * superconducting_tf_coil_variables.a_tf_coil_inboard_insulation
            / tfcoil_variables.a_tf_inboard_total
        )

        # Total cross-sectional area of the bucking cylindre and the outer support
        # support structure per coil [m2]
        # physics_variables.itart = 1 : Only valid at mid-plane
        tfcoil_variables.a_tf_coil_inboard_case = (
            tfcoil_variables.a_tf_inboard_total / tfcoil_variables.n_tf_coils
        ) - superconducting_tf_coil_variables.a_tf_wp_with_insulation

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
        if superconducting_tf_coil_variables.a_tf_wp_with_insulation < 0.0e0:
            logger.error(
                f"Winding pack cross-section problem... {superconducting_tf_coil_variables.a_tf_wp_with_insulation=} "
                f"{tfcoil_variables.dr_tf_wp_with_insulation=}"
            )

        elif superconducting_tf_coil_variables.a_tf_wp_no_insulation < 0.0e0:
            logger.error(
                f"Negative cable space dimension. {superconducting_tf_coil_variables.a_tf_wp_no_insulation=}"
            )

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
                superconducting_tf_coil_variables.is_leg_cp_temp_same = 1
                tfcoil_variables.temp_tf_legs_outboard = tfcoil_variables.temp_cp_average

            # Leg resistivity (different leg temperature as separate cooling channels)
            if tfcoil_variables.i_tf_sup == 0:
                tfcoil_variables.rho_tf_leg = (
                    tfcoil_variables.frholeg
                    * (
                        1.86e0
                        + 0.00393e0 * (tfcoil_variables.temp_tf_legs_outboard - 293.15e0)
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
            if superconducting_tf_coil_variables.is_leg_cp_temp_same == 1:
                tfcoil_variables.temp_tf_legs_outboard = -1.0e0

            # Centrepost resisitivity and conductor/insulation volume

            (
                tfcoil_variables.a_cp_cool,
                tfcoil_variables.vol_cond_cp,
                tfcoil_variables.p_cp_resistive,
                superconducting_tf_coil_variables.vol_ins_cp,
                superconducting_tf_coil_variables.vol_case_cp,
                superconducting_tf_coil_variables.vol_gr_ins_cp,
            ) = self.cpost(
                build_variables.r_tf_inboard_in,
                build_variables.r_tf_inboard_out,
                build_variables.r_cp_top,
                superconducting_tf_coil_variables.z_cp_top,
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
        superconducting_tf_coil_variables.a_leg_gr_ins = (
            tfcoil_variables.a_tf_leg_outboard
            - (
                tfcoil_variables.dx_tf_inboard_out_toroidal
                - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
            )
            * (
                build_variables.dr_tf_outboard
                - 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
            )
        )

        # Outboard leg turns insulation area per coil [m2]
        superconducting_tf_coil_variables.a_leg_ins = (
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
        superconducting_tf_coil_variables.a_leg_cond = (
            1.0e0 - tfcoil_variables.f_a_tf_cool_outboard
        ) * (
            tfcoil_variables.a_tf_leg_outboard
            - superconducting_tf_coil_variables.a_leg_gr_ins
            - superconducting_tf_coil_variables.a_leg_ins
        )
        # ---

        if physics_variables.itart == 1:
            # Outer leg resistive power loss
            # ---
            # TF outboard leg's resistance calculation (per leg) [ohm]
            tfcoil_variables.res_tf_leg = (
                tfcoil_variables.rho_tf_leg
                * tfcoil_variables.len_tf_coil
                / superconducting_tf_coil_variables.a_leg_cond
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
                / (
                    superconducting_tf_coil_variables.a_leg_cond
                    * tfcoil_variables.n_tf_coils
                )
            )

            # tfcoil_variables.p_cp_resistive containts the the total resistive power losses
            tfcoil_variables.p_tf_leg_resistive = 0.0e0

            # No joints if physics_variables.itart = 0
            tfcoil_variables.p_tf_joints_resistive = 0.0e0

    def resistive_tf_coil_areas_and_masses(self):
        """
        Calculate the areas and masses of the resistive TF coil
        """

        vol_case = 0.0e0  # Total TF case volume [m3]
        vol_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_gr_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_cond = 0.0e0  # Total conductor insulator volume [m3]
        vol_ins_leg = 0.0e0  # Outboard leg turn isulation volume [m3]
        vol_gr_ins_leg = 0.0e0  # Outboard leg turn insulation volume [m3]
        vol_cond_leg = 0.0e0  # Outboard leg conductor insulator volume [m3]

        # Volumes
        # -------
        # CP with joints
        # ---
        if physics_variables.itart == 1:
            # Total volume of one outerleg [m3]
            tfcoil_variables.voltfleg = (
                tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_leg_outboard
            )

            # Outboard leg TF conductor volume [m3]
            vol_cond_leg = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_cond
            )

            # Total TF conductor volume [m3]
            vol_cond = (
                tfcoil_variables.vol_cond_cp + tfcoil_variables.n_tf_coils * vol_cond_leg
            )

            # Outboard leg TF turn insulation layer volume (per leg) [m3]
            vol_ins_leg = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_ins
            )

            # Total turn insulation layer volume [m3]
            vol_ins = (
                superconducting_tf_coil_variables.vol_ins_cp
                + tfcoil_variables.n_tf_coils * vol_ins_leg
            )

            # Ouboard leg TF ground insulation layer volume (per leg) [m3]
            vol_gr_ins_leg = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_gr_ins
            )

            # Total ground insulation layer volume [m3]
            vol_gr_ins = (
                superconducting_tf_coil_variables.vol_gr_ins_cp
                + tfcoil_variables.n_tf_coils * vol_gr_ins_leg
            )

            # Total volume of the CP casing [m3]
            # Rem : no outer leg case
            vol_case = superconducting_tf_coil_variables.vol_case_cp

        # No joints
        # ---
        else:
            # Total TF outer leg conductor volume [m3]
            vol_cond = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_cond
                * tfcoil_variables.n_tf_coils
            )

            # Total turn insulation layer volume [m3]
            vol_ins = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_ins
                * tfcoil_variables.n_tf_coils
            )

            # Total ground insulation volume [m3]
            vol_gr_ins = (
                tfcoil_variables.len_tf_coil
                * superconducting_tf_coil_variables.a_leg_gr_ins
                * tfcoil_variables.n_tf_coils
            )

            # Total case volume [m3]
            vol_case = (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.a_tf_coil_inboard_case
                * tfcoil_variables.n_tf_coils
            )

        # ---
        # -------

        # Copper magnets casing/conductor weights per coil [kg]
        if tfcoil_variables.i_tf_sup == 0:
            tfcoil_variables.m_tf_coil_case = (
                fwbs_variables.den_steel * vol_case / tfcoil_variables.n_tf_coils
            )  # Per TF leg, no casing for outer leg
            tfcoil_variables.m_tf_coil_copper = (
                constants.den_copper * vol_cond / tfcoil_variables.n_tf_coils
            )
            tfcoil_variables.whtconal = 0.0e0

            # Outer legs/CP weights
            if physics_variables.itart == 1:
                # Weight of all the TF legs
                tfcoil_variables.whttflgs = tfcoil_variables.n_tf_coils * (
                    constants.den_copper * vol_cond_leg
                    + tfcoil_variables.den_tf_wp_turn_insulation
                    * (vol_ins_leg + vol_gr_ins_leg)
                )

                # CP weight
                tfcoil_variables.whtcp = (
                    constants.den_copper * tfcoil_variables.vol_cond_cp
                    + tfcoil_variables.den_tf_wp_turn_insulation
                    * (
                        superconducting_tf_coil_variables.vol_ins_cp
                        + superconducting_tf_coil_variables.vol_gr_ins_cp
                    )
                    + superconducting_tf_coil_variables.vol_case_cp
                    * fwbs_variables.den_steel
                )

        # Cryo-aluminium conductor weights
        # Casing made of re-inforced aluminium alloy
        elif tfcoil_variables.i_tf_sup == 2:
            # Casing weight (CP only if physics_variables.itart = 1)bper leg/coil
            tfcoil_variables.m_tf_coil_case = (
                constants.den_aluminium * vol_case / tfcoil_variables.n_tf_coils
            )
            tfcoil_variables.m_tf_coil_copper = 0.0e0
            tfcoil_variables.whtconal = (
                constants.den_aluminium * vol_cond / tfcoil_variables.n_tf_coils
            )

            # Outer legs/CP weights
            if physics_variables.itart == 1:
                # Weight of all the TF legs
                tfcoil_variables.whttflgs = tfcoil_variables.n_tf_coils * (
                    constants.den_aluminium * vol_cond_leg
                    + tfcoil_variables.den_tf_wp_turn_insulation
                    * (vol_ins_leg + vol_gr_ins_leg)
                )

                # CP weight
                tfcoil_variables.whtcp = (
                    constants.den_aluminium * tfcoil_variables.vol_cond_cp
                    + tfcoil_variables.den_tf_wp_turn_insulation
                    * (
                        superconducting_tf_coil_variables.vol_ins_cp
                        + superconducting_tf_coil_variables.vol_gr_ins_cp
                    )
                    + superconducting_tf_coil_variables.vol_case_cp
                    * fwbs_variables.den_steel
                )

        # Turn insulation mass [kg]
        tfcoil_variables.m_tf_coil_wp_turn_insulation = (
            tfcoil_variables.den_tf_wp_turn_insulation
            * vol_ins
            / tfcoil_variables.n_tf_coils
        )

        # Ground wall insulation layer weight
        tfcoil_variables.m_tf_coil_wp_insulation = (
            tfcoil_variables.den_tf_wp_turn_insulation
            * vol_gr_ins
            / tfcoil_variables.n_tf_coils
        )

        # Total weight
        tfcoil_variables.m_tf_coils_total = (
            tfcoil_variables.m_tf_coil_case
            + tfcoil_variables.m_tf_coil_copper
            + tfcoil_variables.whtconal
            + tfcoil_variables.m_tf_coil_wp_turn_insulation
            + tfcoil_variables.m_tf_coil_wp_insulation
        ) * tfcoil_variables.n_tf_coils

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


class CopperTFCoil(ResistiveTFCoil):
    """
    Copper TF coil class for resistive TF coil calculations.
    Inherits from ResistiveTFCoil and implements specific methods for copper TF coils.
    """


class AluminiumTFCoil(ResistiveTFCoil):
    """
    Aluminium TF coil class for resistive TF coil calculations.
    Inherits from ResistiveTFCoil and implements specific methods for aluminium TF coils.
    """
