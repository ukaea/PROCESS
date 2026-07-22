"""Resistive TF coil model."""

import logging

import numba
import numpy as np

from process.core import constants
from process.core import process_output as po
from process.models.tfcoil.base import TFCoil, TFConductorModel

logger = logging.getLogger(__name__)

EPS = np.finfo(1.0).eps


class ResistiveTFCoil(TFCoil):
    """Class for resistive TF coil calculations."""

    def __init__(self):
        self.outfile = constants.NOUT

    def output(self):
        """Run main tfcoil subroutine with outputting."""
        self.run(output=True)
        self.output_general_resistive_tf_info()

    def run(self, output: bool = False):
        """Run main tfcoil subroutine without outputting.

        Parameters
        ----------
        output: bool

        """
        # Set up TF values share by all coil types
        self.run_base_tf()

        # Set the peak inboard field with ripple to be the same as the symmetric peak
        # field as a default for resistive TF coils. The conductor segments are all
        # parallel and there is no ripple in the field.
        # This also allows the peak inboard field constraint to be used for resistive
        # TF coils, which is important for the design of the coil.
        self.data.tfcoil.b_tf_inboard_peak_with_ripple = (
            self.data.tfcoil.b_tf_inboard_peak_symmetric
        )

        self.res_tf_internal_geom()
        self.tf_res_heating()

        self.data.tfcoil.ind_tf_coil = self.tf_coil_self_inductance(
            dr_tf_inboard=self.data.build.dr_tf_inboard,
            r_tf_arc=self.data.tfcoil.r_tf_arc,
            z_tf_arc=self.data.tfcoil.z_tf_arc,
            itart=self.data.physics.itart,
            i_tf_shape=self.data.tfcoil.i_tf_shape,
            z_tf_inside_half=self.data.build.z_tf_inside_half,
            dr_tf_outboard=self.data.build.dr_tf_outboard,
            r_tf_outboard_mid=self.data.build.r_tf_outboard_mid,
            r_tf_inboard_mid=self.data.build.r_tf_inboard_mid,
        )

        (
            self.data.tfcoil.e_tf_magnetic_stored_total,
            self.data.tfcoil.e_tf_magnetic_stored_total_gj,
            self.data.tfcoil.e_tf_coil_magnetic_stored,
        ) = self.tf_stored_magnetic_energy(
            ind_tf_coil=self.data.tfcoil.ind_tf_coil,
            c_tf_total=self.data.tfcoil.c_tf_total,
            n_tf_coils=self.data.tfcoil.n_tf_coils,
        )

        (
            self.data.tfcoil.cforce,
            self.data.tfcoil.vforce,
            self.data.tfcoil.vforce_outboard,
            self.data.superconducting_tfcoil.vforce_inboard_tot,
            self.data.tfcoil.f_vforce_inboard,
        ) = self.tf_field_and_force(
            i_tf_sup=self.data.tfcoil.i_tf_sup,
            r_tf_wp_inboard_outer=self.data.superconducting_tfcoil.r_tf_wp_inboard_outer,
            r_tf_wp_inboard_inner=self.data.superconducting_tfcoil.r_tf_wp_inboard_inner,
            r_tf_outboard_in=self.data.superconducting_tfcoil.r_tf_outboard_in,
            dx_tf_wp_insulation=self.data.tfcoil.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=self.data.tfcoil.dx_tf_wp_insertion_gap,
            b_tf_inboard_peak_symmetric=self.data.tfcoil.b_tf_inboard_peak_symmetric,
            c_tf_total=self.data.tfcoil.c_tf_total,
            n_tf_coils=self.data.tfcoil.n_tf_coils,
            dr_tf_plasma_case=self.data.tfcoil.dr_tf_plasma_case,
            rmajor=self.data.physics.rmajor,
            b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
            r_cp_top=self.data.build.r_cp_top,
            itart=self.data.physics.itart,
            i_cp_joints=self.data.tfcoil.i_cp_joints,
            f_vforce_inboard=self.data.tfcoil.f_vforce_inboard,
        )

        # Calculate TF coil areas and masses
        self.generic_tf_coil_area_and_masses()
        self.resistive_tf_coil_areas_and_masses()

        # Do stress calculations (writes the stress output)
        if output:
            self.data.tfcoil.n_rad_per_layer = 500

        try:  # noqa: PLW0717
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
                self.data.tfcoil.sig_tf_wp,
                sig_tf_case,
                sig_tf_cs_bucked,
                str_wp,
                casestr,
                insstrain,
                sig_tf_wp_av_z,
            ) = self.stresscl(
                n_tf_layer=int(self.data.tfcoil.n_tf_stress_layers),
                n_radial_array=int(self.data.tfcoil.n_rad_per_layer),
                n_tf_wp_stress_layers=int(self.data.tfcoil.n_tf_wp_stress_layers),
                i_tf_bucking=int(self.data.tfcoil.i_tf_bucking),
                r_tf_inboard_in=float(self.data.build.r_tf_inboard_in),
                dr_bore=self.data.build.dr_bore,
                dr_cs=self.data.build.dr_cs,
                i_tf_inside_cs=self.data.build.i_tf_inside_cs,
                dr_tf_inboard=self.data.build.dr_tf_inboard,
                dr_cs_tf_gap=self.data.build.dr_cs_tf_gap,
                i_pf_conductor=self.data.pf_coil.i_pf_conductor,
                j_cs_flat_top_end=self.data.pf_coil.j_cs_flat_top_end,
                j_cs_pulse_start=self.data.pf_coil.j_cs_pulse_start,
                c_pf_coil_turn_peak_input=self.data.pf_coil.c_pf_coil_turn_peak_input,
                n_pf_coils_in_group=self.data.pf_coil.n_pf_coils_in_group,
                f_dr_dz_cs_turn=self.data.pf_coil.f_dr_dz_cs_turn,
                radius_cs_turn_corners=self.data.pf_coil.radius_cs_turn_corners,
                f_a_cs_turn_steel=self.data.pf_coil.f_a_cs_turn_steel,
                eyoung_steel=self.data.tfcoil.eyoung_steel,
                poisson_steel=self.data.tfcoil.poisson_steel,
                eyoung_cond_axial=self.data.tfcoil.eyoung_cond_axial,
                poisson_cond_axial=self.data.tfcoil.poisson_cond_axial,
                eyoung_cond_trans=self.data.tfcoil.eyoung_cond_trans,
                poisson_cond_trans=self.data.tfcoil.poisson_cond_trans,
                eyoung_ins=self.data.tfcoil.eyoung_ins,
                poisson_ins=self.data.tfcoil.poisson_ins,
                dx_tf_turn_insulation=self.data.tfcoil.dx_tf_turn_insulation,
                eyoung_copper=self.data.tfcoil.eyoung_copper,
                poisson_copper=self.data.tfcoil.poisson_copper,
                i_tf_sup=self.data.tfcoil.i_tf_sup,
                eyoung_res_tf_buck=self.data.tfcoil.eyoung_res_tf_buck,
                r_tf_wp_inboard_inner=self.data.superconducting_tfcoil.r_tf_wp_inboard_inner,
                tan_theta_coil=self.data.superconducting_tfcoil.tan_theta_coil,
                rad_tf_coil_inboard_toroidal_half=self.data.superconducting_tfcoil.rad_tf_coil_inboard_toroidal_half,
                r_tf_wp_inboard_outer=self.data.superconducting_tfcoil.r_tf_wp_inboard_outer,
                a_tf_coil_inboard_steel=self.data.superconducting_tfcoil.a_tf_coil_inboard_steel,
                a_tf_plasma_case=self.data.superconducting_tfcoil.a_tf_plasma_case,
                a_tf_coil_nose_case=self.data.superconducting_tfcoil.a_tf_coil_nose_case,
                dx_tf_wp_insertion_gap=self.data.tfcoil.dx_tf_wp_insertion_gap,
                dx_tf_wp_insulation=self.data.tfcoil.dx_tf_wp_insulation,
                n_tf_coil_turns=self.data.tfcoil.n_tf_coil_turns,
                i_tf_turns_integer=int(self.data.tfcoil.i_tf_turns_integer),
                dx_tf_turn_cable_space_average=self.data.superconducting_tfcoil.dx_tf_turn_cable_space_average,
                dr_tf_turn_cable_space=self.data.superconducting_tfcoil.dr_tf_turn_cable_space,
                dia_tf_turn_coolant_channel=self.data.tfcoil.dia_tf_turn_coolant_channel,
                f_a_tf_turn_cable_copper=self.data.tfcoil.f_a_tf_turn_cable_copper,
                dx_tf_turn_steel=self.data.tfcoil.dx_tf_turn_steel,
                dx_tf_side_case_average=self.data.superconducting_tfcoil.dx_tf_side_case_average,
                dx_tf_wp_toroidal_average=self.data.superconducting_tfcoil.dx_tf_wp_toroidal_average,
                a_tf_coil_inboard_insulation=self.data.superconducting_tfcoil.a_tf_coil_inboard_insulation,
                a_tf_wp_steel=self.data.tfcoil.a_tf_wp_steel,
                a_tf_wp_conductor=self.data.tfcoil.a_tf_wp_conductor,
                a_tf_wp_with_insulation=self.data.superconducting_tfcoil.a_tf_wp_with_insulation,
                eyoung_al=self.data.tfcoil.eyoung_al,
                poisson_al=self.data.tfcoil.poisson_al,
                fcoolcp=self.data.tfcoil.fcoolcp,
                n_tf_graded_layers=self.data.tfcoil.n_tf_graded_layers,
                c_tf_total=self.data.tfcoil.c_tf_total,
                dr_tf_plasma_case=self.data.tfcoil.dr_tf_plasma_case,
                i_tf_stress_model=self.data.tfcoil.i_tf_stress_model,
                vforce_inboard_tot=self.data.superconducting_tfcoil.vforce_inboard_tot,
                i_tf_tresca=self.data.tfcoil.i_tf_tresca,
                a_tf_coil_inboard_case=self.data.tfcoil.a_tf_coil_inboard_case,
                vforce=self.data.tfcoil.vforce,
                a_tf_turn_steel=self.data.tfcoil.a_tf_turn_steel,
                a_cs_poloidal=self.data.pf_coil.a_cs_poloidal,
            )

            self.data.tfcoil.sig_tf_case = (
                self.data.tfcoil.sig_tf_case
                if self.data.tfcoil.sig_tf_case is None
                else sig_tf_case
            )

            self.data.tfcoil.sig_tf_cs_bucked = (
                self.data.tfcoil.sig_tf_cs_bucked
                if self.data.tfcoil.sig_tf_cs_bucked is None
                else sig_tf_cs_bucked
            )

            self.data.tfcoil.str_wp = (
                self.data.tfcoil.str_wp if self.data.tfcoil.str_wp is None else str_wp
            )

            self.data.tfcoil.casestr = (
                self.data.tfcoil.casestr if self.data.tfcoil.casestr is None else casestr
            )

            self.data.tfcoil.insstrain = (
                self.data.tfcoil.insstrain
                if self.data.tfcoil.insstrain is None
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
                    "Invalid stress model (r_tf_inboard = 0), stress constraint "
                    "switched off"
                )
                self.data.tfcoil.sig_tf_case = 0.0e0
                self.data.tfcoil.sig_tf_wp = 0.0e0
        if output:
            self.output_general_tf_info()

    def res_tf_internal_geom(self):
        """
        Resistive TF turn geometry, equivalent to winding_pack subroutines
        """
        self.data.superconducting_tfcoil.r_tf_wp_inboard_inner = (
            self.data.build.r_tf_inboard_in + self.data.tfcoil.dr_tf_nose_case
        )
        self.data.superconducting_tfcoil.r_tf_wp_inboard_outer = (
            self.data.build.r_tf_inboard_out - self.data.tfcoil.dr_tf_plasma_case
        )

        # Conductor layer radial thickness at centercollumn top [m]
        if self.data.physics.itart == 1:
            self.data.superconducting_tfcoil.dr_tf_wp_top = (
                self.data.build.r_cp_top
                - self.data.tfcoil.dr_tf_plasma_case
                - self.data.tfcoil.dr_tf_nose_case
                - self.data.build.r_tf_inboard_in
            )

        # Number of turns
        # Set by user (no turn structure by default, i.e.
        # self.data.tfcoil.n_tf_coil_turns = 1 )
        if (
            abs(self.data.tfcoil.n_tf_coil_turns)
            < np.finfo(float(self.data.tfcoil.n_tf_coil_turns)).eps
        ):
            self.data.tfcoil.n_tf_coil_turns = 1.0e0

        # Total mid-plane cross-sectional area of winding pack, [m2]
        # including the surrounding ground-wall insulation layer
        self.data.superconducting_tfcoil.a_tf_wp_with_insulation = (
            np.pi
            * (
                self.data.superconducting_tfcoil.r_tf_wp_inboard_outer**2
                - self.data.superconducting_tfcoil.r_tf_wp_inboard_inner**2
            )
            / self.data.tfcoil.n_tf_coils
        )

        # Area of the front case, the plasma-facing case of the inner TF coil [m2]
        self.data.superconducting_tfcoil.a_tf_plasma_case = (
            np.pi
            * (
                (
                    self.data.superconducting_tfcoil.r_tf_wp_inboard_outer
                    + self.data.tfcoil.dr_tf_plasma_case
                )
                ** 2
                - self.data.superconducting_tfcoil.r_tf_wp_inboard_outer**2
            )
            / self.data.tfcoil.n_tf_coils
        )

        # WP mid-plane cross-section excluding ground insulation per coil [m2]
        self.data.superconducting_tfcoil.a_tf_wp_no_insulation = (
            np.pi
            * (
                (
                    self.data.superconducting_tfcoil.r_tf_wp_inboard_outer
                    - self.data.tfcoil.dx_tf_wp_insulation
                )
                ** 2
                - (
                    self.data.superconducting_tfcoil.r_tf_wp_inboard_inner
                    + self.data.tfcoil.dx_tf_wp_insulation
                )
                ** 2
            )
            / self.data.tfcoil.n_tf_coils
            - 2.0e0
            * self.data.tfcoil.dx_tf_wp_insulation
            * (
                self.data.tfcoil.dr_tf_wp_with_insulation
                - 2.0e0 * self.data.tfcoil.dx_tf_wp_insulation
            )
        )

        # Ground insulation cross-section area per coil [m2]
        self.data.superconducting_tfcoil.a_tf_wp_ground_insulation = (
            self.data.superconducting_tfcoil.a_tf_wp_with_insulation
            - self.data.superconducting_tfcoil.a_tf_wp_no_insulation
        )

        # Exact mid-plane cross-section area of the conductor per TF coil [m2]
        self.data.tfcoil.a_res_tf_coil_conductor = np.pi * (
            (
                self.data.superconducting_tfcoil.r_tf_wp_inboard_outer
                - self.data.tfcoil.dx_tf_wp_insulation
                - self.data.tfcoil.dx_tf_turn_insulation
            )
            ** 2
            - (
                self.data.superconducting_tfcoil.r_tf_wp_inboard_inner
                + self.data.tfcoil.dx_tf_wp_insulation
                + self.data.tfcoil.dx_tf_turn_insulation
            )
            ** 2
        ) / self.data.tfcoil.n_tf_coils - (
            self.data.tfcoil.dr_tf_wp_with_insulation
            - 2.0e0
            * (
                self.data.tfcoil.dx_tf_wp_insulation
                + self.data.tfcoil.dx_tf_turn_insulation
            )
        ) * 2.0e0 * (
            self.data.tfcoil.dx_tf_wp_insulation
            + self.data.tfcoil.dx_tf_turn_insulation * self.data.tfcoil.n_tf_coil_turns
        )
        self.data.tfcoil.a_res_tf_coil_conductor *= 1.0e0 - self.data.tfcoil.fcoolcp

        # Inter turn insulation area per coil [m2]
        self.data.tfcoil.a_tf_coil_wp_turn_insulation = (
            self.data.superconducting_tfcoil.a_tf_wp_no_insulation
            - self.data.tfcoil.a_res_tf_coil_conductor
            / (1.0e0 - self.data.tfcoil.fcoolcp)
        )

        # Total insulation cross-section per coil [m2]
        self.data.superconducting_tfcoil.a_tf_coil_inboard_insulation = (
            self.data.tfcoil.a_tf_coil_wp_turn_insulation
            + self.data.superconducting_tfcoil.a_tf_wp_ground_insulation
        )

        # Insulation fraction [-]
        self.data.superconducting_tfcoil.f_a_tf_coil_inboard_insulation = (
            self.data.tfcoil.n_tf_coils
            * self.data.superconducting_tfcoil.a_tf_coil_inboard_insulation
            / self.data.tfcoil.a_tf_inboard_total
        )

        # Total cross-sectional area of the bucking cylinder and the outer support
        # support structure per coil [m2]
        # self.data.physics.itart = 1 : Only valid at mid-plane
        self.data.tfcoil.a_tf_coil_inboard_case = (
            self.data.tfcoil.a_tf_inboard_total / self.data.tfcoil.n_tf_coils
        ) - self.data.superconducting_tfcoil.a_tf_wp_with_insulation

        # Current per turn
        self.data.tfcoil.c_tf_turn = self.data.tfcoil.c_tf_total / (
            self.data.tfcoil.n_tf_coil_turns * self.data.tfcoil.n_tf_coils
        )

        # Exact current density on TF outboard legs
        self.data.tfcoil.cdtfleg = self.data.tfcoil.c_tf_total / (
            (1.0e0 - self.data.tfcoil.fcoolcp)
            * (
                self.data.tfcoil.dx_tf_inboard_out_toroidal
                - 2.0e0
                * (
                    self.data.tfcoil.n_tf_coil_turns
                    * self.data.tfcoil.dx_tf_turn_insulation
                    + self.data.tfcoil.dx_tf_wp_insulation
                )
            )
            * (
                self.data.build.dr_tf_outboard
                - 2.0e0
                * (
                    self.data.tfcoil.dx_tf_turn_insulation
                    + self.data.tfcoil.dx_tf_wp_insulation
                )
            )
        )

        # Reporting negative WP areas issues
        if self.data.superconducting_tfcoil.a_tf_wp_with_insulation < 0.0e0:
            logger.error(
                f"Winding pack cross-section problem... "
                f"{self.data.superconducting_tfcoil.a_tf_wp_with_insulation=} "
                f"{self.data.tfcoil.dr_tf_wp_with_insulation=}"
            )

        elif self.data.superconducting_tfcoil.a_tf_wp_no_insulation < 0.0e0:
            logger.error(
                f"Negative cable space dimension. "
                f"{self.data.superconducting_tfcoil.a_tf_wp_no_insulation=}"
            )

    def tf_res_heating(self):
        """Calculate resistive heating for resistive magnets.

        This method calculates the resistive heating for resistive magnets.
        It considers the following scenarios:
        - Clamped joints in superconductors might have resistive power losses on
          the joints.
        - Sliding joints might have a region of high resistivity.

        Notes
        -----
        - The copper resistivity is set to be that for GLIDCOP AL-15 at 20°C for copper
          (i_tf_sup = 0).
        - The coefficient of resistivity is set to be that of pure copper

        References
        ----------
            - https://www.spotweldingconsultants.com/GlidCop_AL_15.pdf

            - https://cirris.com/temperature-coefficient-of-copper/
        """
        # Resistivity of the Glidcop copper centerpost
        if self.data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
            self.data.tfcoil.rho_cp = (
                # 1.86 is the resistivity at `20°C` for GLIDCOP AL-15
                # 0.00393 is the coefficient of resistivity for copper
                self.data.tfcoil.frhocp
                * (1.86e0 + 0.00393e0 * (self.data.tfcoil.temp_cp_average - 293.15e0))
                * 1.0e-8
            )

        # Resistivity of the aluminium centerpost
        if self.data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            self.data.tfcoil.rho_cp = self.data.tfcoil.frhocp * (
                2.00016e-14 * self.data.tfcoil.temp_cp_average**3
                - 6.75384e-13 * self.data.tfcoil.temp_cp_average**2
                + 8.89159e-12 * self.data.tfcoil.temp_cp_average
            )

        # Calculations dedicated for configurations with CP
        if self.data.physics.itart == 1:
            # Tricky trick to make the leg / CP temperatures the same
            if (
                abs(self.data.tfcoil.temp_tf_legs_outboard + 1.0e0)
                < np.finfo(float(self.data.tfcoil.temp_tf_legs_outboard)).eps
            ):
                self.data.superconducting_tfcoil.is_leg_cp_temp_same = 1
                self.data.tfcoil.temp_tf_legs_outboard = self.data.tfcoil.temp_cp_average

            # Leg resistivity (different leg temperature as separate cooling channels)
            if self.data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
                self.data.tfcoil.rho_tf_leg = (
                    self.data.tfcoil.frholeg
                    * (
                        1.86e0
                        + 0.00393e0 * (self.data.tfcoil.temp_tf_legs_outboard - 293.15e0)
                    )
                    * 1.0e-8
                )
            elif self.data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
                self.data.tfcoil.rho_tf_leg = self.data.tfcoil.frholeg * (
                    2.00016e-14 * self.data.tfcoil.temp_tf_legs_outboard**3
                    - 6.75384e-13 * self.data.tfcoil.temp_tf_legs_outboard**2
                    + 8.89159e-12 * self.data.tfcoil.temp_tf_legs_outboard
                )

            # Tricky trick to make the leg / CP temperatures the same
            if self.data.superconducting_tfcoil.is_leg_cp_temp_same == 1:
                self.data.tfcoil.temp_tf_legs_outboard = -1.0e0

            # Centrepost resisitivity and conductor/insulation volume

            (
                self.data.tfcoil.a_cp_cool,
                self.data.tfcoil.vol_cond_cp,
                self.data.tfcoil.p_cp_resistive,
                self.data.superconducting_tfcoil.vol_ins_cp,
                self.data.superconducting_tfcoil.vol_case_cp,
                self.data.superconducting_tfcoil.vol_gr_ins_cp,
            ) = self.cpost(
                self.data.build.r_tf_inboard_in,
                self.data.build.r_tf_inboard_out,
                self.data.build.r_cp_top,
                self.data.superconducting_tfcoil.z_cp_top,
                self.data.build.z_tf_inside_half + self.data.build.dr_tf_outboard,
                self.data.tfcoil.dr_tf_nose_case,
                self.data.tfcoil.dr_tf_plasma_case,
                self.data.tfcoil.dx_tf_wp_insulation,
                self.data.tfcoil.dx_tf_turn_insulation,
                self.data.tfcoil.n_tf_coil_turns,
                self.data.tfcoil.c_tf_total,
                self.data.tfcoil.rho_cp,
                self.data.tfcoil.fcoolcp,
                self.data.tfcoil.n_tf_coils,
            )

        # Leg cross-section areas
        # Rem : For self.data.physics.itart = 1, these quantities corresponds to
        # the outer leg only
        # ---
        # Leg ground insulation area per coil [m2]
        self.data.superconducting_tfcoil.a_leg_gr_ins = (
            self.data.tfcoil.a_tf_leg_outboard
            - (
                self.data.tfcoil.dx_tf_inboard_out_toroidal
                - 2.0e0 * self.data.tfcoil.dx_tf_wp_insulation
            )
            * (
                self.data.build.dr_tf_outboard
                - 2.0e0 * self.data.tfcoil.dx_tf_wp_insulation
            )
        )

        # Outboard leg turns insulation area per coil [m2]
        self.data.superconducting_tfcoil.a_leg_ins = (
            2.0e0
            * self.data.tfcoil.dx_tf_turn_insulation
            * (
                self.data.tfcoil.dx_tf_inboard_out_toroidal
                - 2.0e0 * self.data.tfcoil.dx_tf_wp_insulation
            )
            + 2.0e0
            * self.data.tfcoil.dx_tf_turn_insulation
            * self.data.tfcoil.n_tf_coil_turns
            * (
                self.data.build.dr_tf_outboard
                - 2.0e0
                * (
                    self.data.tfcoil.dx_tf_turn_insulation
                    + self.data.tfcoil.dx_tf_wp_insulation
                )
            )
        )  # toroidal direction + radial direction

        # Exact TF outboard leg conductor area per coil [m2]
        self.data.superconducting_tfcoil.a_leg_cond = (
            1.0e0 - self.data.tfcoil.f_a_tf_cool_outboard
        ) * (
            self.data.tfcoil.a_tf_leg_outboard
            - self.data.superconducting_tfcoil.a_leg_gr_ins
            - self.data.superconducting_tfcoil.a_leg_ins
        )
        # ---

        if self.data.physics.itart == 1:
            # Outer leg resistive power loss
            # ---
            # TF outboard leg's resistance calculation (per leg) [ohm]
            self.data.tfcoil.res_tf_leg = (
                self.data.tfcoil.rho_tf_leg
                * self.data.tfcoil.len_tf_coil
                / self.data.superconducting_tfcoil.a_leg_cond
            )

            # TF outer leg resistive power (TOTAL) [W]
            self.data.tfcoil.p_tf_leg_resistive = (
                self.data.tfcoil.res_tf_leg
                * (self.data.tfcoil.c_tf_total / self.data.tfcoil.n_tf_coils) ** 2
            ) * self.data.tfcoil.n_tf_coils
            # ---

            # Sliding joints resistive heating
            # ---
            if self.data.tfcoil.i_cp_joints != 0:
                # Number of contact area per joint (all legs)
                n_contact_tot = (
                    self.data.tfcoil.n_tf_joints_contact
                    * np.round(self.data.tfcoil.n_tf_coil_turns)
                    * np.round(self.data.tfcoil.n_tf_coils)
                )

                # Area of joint contact (all legs)
                a_joints = (
                    self.data.build.dr_tf_outboard
                    * self.data.tfcoil.th_joint_contact
                    * n_contact_tot
                )

                # Total joints resistive power losses
                self.data.tfcoil.p_tf_joints_resistive = (
                    self.data.tfcoil.n_tf_joints
                    * self.data.tfcoil.rho_tf_joints
                    * self.data.tfcoil.c_tf_total**2
                    / a_joints
                )
            else:
                # Joints resistance to be evaluated for SC
                self.data.tfcoil.p_tf_joints_resistive = 0.0e0

            # ---

        # Case of a resistive magnet without joints
        # ***
        else:
            # TF resistive powers
            self.data.tfcoil.p_cp_resistive = (
                self.data.tfcoil.rho_cp
                * self.data.tfcoil.c_tf_total**2
                * self.data.tfcoil.len_tf_coil
                / (
                    self.data.superconducting_tfcoil.a_leg_cond
                    * self.data.tfcoil.n_tf_coils
                )
            )

            # self.data.tfcoil.p_cp_resistive contains the the total resistive
            # power losses
            self.data.tfcoil.p_tf_leg_resistive = 0.0e0

            # No joints if self.data.physics.itart = 0
            self.data.tfcoil.p_tf_joints_resistive = 0.0e0

    def resistive_tf_coil_areas_and_masses(self):
        """Calculate the areas and masses of the resistive TF coil"""
        vol_case = 0.0e0  # Total TF case volume [m3]
        vol_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_gr_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_cond = 0.0e0  # Total conductor insulator volume [m3]
        vol_ins_leg = 0.0e0  # Outboard leg turn insulation volume [m3]
        vol_gr_ins_leg = 0.0e0  # Outboard leg turn insulation volume [m3]
        vol_cond_leg = 0.0e0  # Outboard leg conductor insulator volume [m3]

        # Volumes
        # -------
        # CP with joints
        # ---
        if self.data.physics.itart == 1:
            # Total volume of one outerleg [m3]
            self.data.tfcoil.voltfleg = (
                self.data.tfcoil.len_tf_coil * self.data.tfcoil.a_tf_leg_outboard
            )

            # Outboard leg TF conductor volume [m3]
            vol_cond_leg = (
                self.data.tfcoil.len_tf_coil
                * self.data.superconducting_tfcoil.a_leg_cond
            )

            # Total TF conductor volume [m3]
            vol_cond = (
                self.data.tfcoil.vol_cond_cp + self.data.tfcoil.n_tf_coils * vol_cond_leg
            )

            # Outboard leg TF turn insulation layer volume (per leg) [m3]
            vol_ins_leg = (
                self.data.tfcoil.len_tf_coil * self.data.superconducting_tfcoil.a_leg_ins
            )

            # Total turn insulation layer volume [m3]
            vol_ins = (
                self.data.superconducting_tfcoil.vol_ins_cp
                + self.data.tfcoil.n_tf_coils * vol_ins_leg
            )

            # Ouboard leg TF ground insulation layer volume (per leg) [m3]
            vol_gr_ins_leg = (
                self.data.tfcoil.len_tf_coil
                * self.data.superconducting_tfcoil.a_leg_gr_ins
            )

            # Total ground insulation layer volume [m3]
            vol_gr_ins = (
                self.data.superconducting_tfcoil.vol_gr_ins_cp
                + self.data.tfcoil.n_tf_coils * vol_gr_ins_leg
            )

            # Total volume of the CP casing [m3]
            # Rem : no outer leg case
            vol_case = self.data.superconducting_tfcoil.vol_case_cp

        # No joints
        # ---
        else:
            # Total TF outer leg conductor volume [m3]
            vol_cond = (
                self.data.tfcoil.len_tf_coil
                * self.data.superconducting_tfcoil.a_leg_cond
                * self.data.tfcoil.n_tf_coils
            )

            # Total turn insulation layer volume [m3]
            vol_ins = (
                self.data.tfcoil.len_tf_coil
                * self.data.superconducting_tfcoil.a_leg_ins
                * self.data.tfcoil.n_tf_coils
            )

            # Total ground insulation volume [m3]
            vol_gr_ins = (
                self.data.tfcoil.len_tf_coil
                * self.data.superconducting_tfcoil.a_leg_gr_ins
                * self.data.tfcoil.n_tf_coils
            )

            # Total case volume [m3]
            vol_case = (
                self.data.tfcoil.len_tf_coil
                * self.data.tfcoil.a_tf_coil_inboard_case
                * self.data.tfcoil.n_tf_coils
            )

        # Copper magnets casing/conductor weights per coil [kg]
        if self.data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
            self.data.tfcoil.m_tf_coil_case = (
                self.data.fwbs.den_steel * vol_case / self.data.tfcoil.n_tf_coils
            )  # Per TF leg, no casing for outer leg
            self.data.tfcoil.m_tf_coil_copper = (
                constants.DEN_COPPER * vol_cond / self.data.tfcoil.n_tf_coils
            )
            self.data.tfcoil.whtconal = 0.0e0

            # Outer legs/CP weights
            if self.data.physics.itart == 1:
                # Weight of all the TF legs
                self.data.tfcoil.whttflgs = self.data.tfcoil.n_tf_coils * (
                    constants.DEN_COPPER * vol_cond_leg
                    + self.data.tfcoil.den_tf_wp_turn_insulation
                    * (vol_ins_leg + vol_gr_ins_leg)
                )

                # CP weight
                self.data.tfcoil.whtcp = (
                    constants.DEN_COPPER * self.data.tfcoil.vol_cond_cp
                    + self.data.tfcoil.den_tf_wp_turn_insulation
                    * (
                        self.data.superconducting_tfcoil.vol_ins_cp
                        + self.data.superconducting_tfcoil.vol_gr_ins_cp
                    )
                    + self.data.superconducting_tfcoil.vol_case_cp
                    * self.data.fwbs.den_steel
                )

        # Cryo-aluminium conductor weights
        # Casing made of re-inforced aluminium alloy
        elif self.data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            # Casing weight (CP only if self.data.physics.itart = 1)bper leg/coil
            self.data.tfcoil.m_tf_coil_case = (
                constants.DEN_ALUMINIUM * vol_case / self.data.tfcoil.n_tf_coils
            )
            self.data.tfcoil.m_tf_coil_copper = 0.0e0
            self.data.tfcoil.whtconal = (
                constants.DEN_ALUMINIUM * vol_cond / self.data.tfcoil.n_tf_coils
            )

            # Outer legs/CP weights
            if self.data.physics.itart == 1:
                # Weight of all the TF legs
                self.data.tfcoil.whttflgs = self.data.tfcoil.n_tf_coils * (
                    constants.DEN_ALUMINIUM * vol_cond_leg
                    + self.data.tfcoil.den_tf_wp_turn_insulation
                    * (vol_ins_leg + vol_gr_ins_leg)
                )

                # CP weight
                self.data.tfcoil.whtcp = (
                    constants.DEN_ALUMINIUM * self.data.tfcoil.vol_cond_cp
                    + self.data.tfcoil.den_tf_wp_turn_insulation
                    * (
                        self.data.superconducting_tfcoil.vol_ins_cp
                        + self.data.superconducting_tfcoil.vol_gr_ins_cp
                    )
                    + self.data.superconducting_tfcoil.vol_case_cp
                    * self.data.fwbs.den_steel
                )

        # Turn insulation mass [kg]
        self.data.tfcoil.m_tf_coil_wp_turn_insulation = (
            self.data.tfcoil.den_tf_wp_turn_insulation
            * vol_ins
            / self.data.tfcoil.n_tf_coils
        )

        # Ground wall insulation layer weight
        self.data.tfcoil.m_tf_coil_wp_insulation = (
            self.data.tfcoil.den_tf_wp_turn_insulation
            * vol_gr_ins
            / self.data.tfcoil.n_tf_coils
        )

        # Total weight
        self.data.tfcoil.m_tf_coils_total = (
            self.data.tfcoil.m_tf_coil_case
            + self.data.tfcoil.m_tf_coil_copper
            + self.data.tfcoil.whtconal
            + self.data.tfcoil.m_tf_coil_wp_turn_insulation
            + self.data.tfcoil.m_tf_coil_wp_insulation
        ) * self.data.tfcoil.n_tf_coils

    def output_general_resistive_tf_info(self) -> None:
        """Output general information on the resistive TF coil
        Should only contain information found in the `ResistiveTFCoil` class
        """
        po.oheadr(self.outfile, "General Resistive TF Coil Parameters")

        # External casing
        po.osubhd(self.outfile, "Bucking cylinder information:")
        po.ovarre(
            self.outfile,
            "Casing cross section area (per leg) (m2)",
            "(a_tf_coil_inboard_case)",
            self.data.tfcoil.a_tf_coil_inboard_case,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg case plasma side wall thickness (m)",
            "(dr_tf_plasma_case)",
            self.data.tfcoil.dr_tf_plasma_case,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg plasma case area (m^2)",
            "(a_tf_plasma_case)",
            self.data.superconducting_tfcoil.a_tf_plasma_case,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg bucking cylinder thickness (m)",
            "(dr_tf_nose_case)",
            self.data.tfcoil.dr_tf_nose_case,
        )
        po.ovarre(
            self.outfile,
            'Inboard leg case inboard "nose" area (m^2)',
            "(a_tf_coil_nose_case)",
            self.data.superconducting_tfcoil.a_tf_coil_nose_case,
        )

        # Conductor layer geometry
        po.osubhd(self.outfile, "Inboard TFC conductor sector geometry:")
        po.ovarre(
            self.outfile,
            "Inboard TFC conductor sector area with gr insulation (per leg) (m2)",
            "(a_tf_wp_with_insulation)",
            self.data.superconducting_tfcoil.a_tf_wp_with_insulation,
        )
        po.ovarre(
            self.outfile,
            "Inboard TFC conductor sector area, NO ground & gap (per leg) (m2)",
            "(a_tf_wp_no_insulation)",
            self.data.superconducting_tfcoil.a_tf_wp_no_insulation,
        )
        po.ovarre(
            self.outfile,
            "Ground wall insulation area (m^2)",
            "(a_tf_wp_ground_insulation)",
            self.data.superconducting_tfcoil.a_tf_wp_ground_insulation,
        )
        po.ovarre(
            self.outfile,
            "Inboard conductor sector radial thickness (m)",
            "(dr_tf_wp_with_insulation)",
            self.data.tfcoil.dr_tf_wp_with_insulation,
        )
        if self.data.physics.itart == 1:
            po.ovarre(
                self.outfile,
                "Central collumn top conductor sector radial thickness (m)",
                "(dr_tf_wp_top)",
                self.data.superconducting_tfcoil.dr_tf_wp_top,
            )

        po.ovarre(
            self.outfile,
            "Ground wall insulation thickness (m)",
            "(dx_tf_wp_insulation)",
            self.data.tfcoil.dx_tf_wp_insulation,
        )
        # Turn info
        po.osubhd(self.outfile, "Coil turn information:")
        po.ovarre(
            self.outfile,
            "Number of turns per TF leg",
            "(n_tf_coil_turns)",
            self.data.tfcoil.n_tf_coil_turns,
        )
        po.ovarre(
            self.outfile,
            "Turn insulation thickness",
            "(dx_tf_turn_insulation)",
            self.data.tfcoil.dx_tf_turn_insulation,
        )
        po.ovarre(
            self.outfile,
            "Mid-plane CP cooling fraction",
            "(fcoolcp)",
            self.data.tfcoil.fcoolcp,
        )
        po.ovarre(
            self.outfile,
            "Area of resistive conductor per coil",
            "(a_res_tf_coil_conductor)",
            self.data.tfcoil.a_res_tf_coil_conductor,
        )
        po.ovarre(
            self.outfile,
            "Outboard leg current per turn (A)",
            "(c_tf_turn)",
            self.data.tfcoil.c_tf_turn,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg conductor volume (m3)",
            "(vol_cond_cp)",
            self.data.tfcoil.vol_cond_cp,
        )
        po.ovarre(
            self.outfile,
            "Outboard leg volume per coil (m3)",
            "(voltfleg)",
            self.data.tfcoil.voltfleg,
        )

        # TF coil radial build
        po.osubhd(self.outfile, "Radial build of TF coil centre-line :")

        radius = self.data.build.r_tf_inboard_in
        po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

        radius += self.data.tfcoil.dr_tf_nose_case
        po.obuild(
            self.outfile,
            "Coil bucking cylindre",
            self.data.tfcoil.dr_tf_nose_case,
            radius,
            "(dr_tf_nose_case)",
        )

        radius += self.data.tfcoil.dx_tf_wp_insulation
        po.obuild(
            self.outfile,
            "Conductor ground insulation",
            self.data.tfcoil.dx_tf_wp_insulation,
            radius,
            "(dx_tf_wp_insulation)",
        )

        radius = (
            radius
            + 0.5e0 * self.data.tfcoil.dr_tf_wp_with_insulation
            - self.data.tfcoil.dx_tf_wp_insulation
        )
        po.obuild(
            self.outfile,
            "Conductor - first half",
            self.data.tfcoil.dr_tf_wp_with_insulation / 2e0
            - self.data.tfcoil.dx_tf_wp_insulation,
            radius,
            "(self.data.tfcoil.dr_tf_wp_with_insulation/2-self.data.tfcoil.dx_tf_wp_insulation)",
        )

        radius = (
            radius
            + 0.5e0 * self.data.tfcoil.dr_tf_wp_with_insulation
            - self.data.tfcoil.dx_tf_wp_insulation
        )
        po.obuild(
            self.outfile,
            "Conductor - second half",
            self.data.tfcoil.dr_tf_wp_with_insulation / 2e0
            - self.data.tfcoil.dx_tf_wp_insulation,
            radius,
            "(self.data.tfcoil.dr_tf_wp_with_insulation/2-self.data.tfcoil.dx_tf_wp_insulation)",
        )

        radius += self.data.tfcoil.dx_tf_wp_insulation
        po.obuild(
            self.outfile,
            "Conductor ground insulation",
            self.data.tfcoil.dx_tf_wp_insulation,
            radius,
            "(dx_tf_wp_insulation)",
        )

        radius += self.data.tfcoil.dr_tf_plasma_case
        po.obuild(
            self.outfile,
            "Plasma side TF coil support",
            self.data.tfcoil.dr_tf_plasma_case,
            radius,
            "(dr_tf_plasma_case)",
        )

        # Radial build consistency check
        if not (
            abs(radius - self.data.build.r_tf_inboard_in - self.data.build.dr_tf_inboard)
            < 10.0e0 * np.finfo(float(radius)).eps
        ):
            logger.error(
                "TF coil dimensions are not consistent. "
                "Radius of plasma-facing side of inner leg should be "
                f"{self.data.build.r_tf_inboard_in + self.data.build.dr_tf_inboard}m"
            )

        tf_total_height = (
            self.data.build.dh_tf_inner_bore + 2 * self.data.build.dr_tf_inboard
        )
        tf_total_width = (
            self.data.build.dr_tf_inner_bore
            + self.data.build.dr_tf_inboard
            + self.data.build.dr_tf_outboard
        )
        po.oblnkl(self.outfile)
        po.obuild(
            self.outfile,
            "Total height and width of TFC [m]",
            tf_total_height,
            tf_total_width,
        )

        # Resistive coil parameters
        po.osubhd(self.outfile, "Resitive loss parameters:")
        if self.data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
            po.ocmmnt(
                self.outfile,
                "Resistive Material : GLIDCOP AL-15 - Dispersion Strengthened Copper",
            )
        elif self.data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            po.ocmmnt(self.outfile, "Resistive Material : Pure Aluminium (99.999+ %)")

        if self.data.physics.itart == 1:
            po.ovarre(
                self.outfile,
                "CP resistivity (ohm.m)",
                "(rho_cp)",
                self.data.tfcoil.rho_cp,
            )
            po.ovarre(
                self.outfile,
                "Leg resistivity (ohm.m)",
                "(rho_tf_leg)",
                self.data.tfcoil.rho_tf_leg,
            )
            po.ovarre(
                self.outfile,
                "CP resistive power loss (W)",
                "(p_cp_resistive)",
                self.data.tfcoil.p_cp_resistive,
            )
            po.ovarre(
                self.outfile,
                "Total legs resitive power loss, (W)",
                "(p_tf_leg_resistive)",
                self.data.tfcoil.p_tf_leg_resistive,
            )
            po.ovarre(
                self.outfile,
                "joints resistive power loss (W)",
                "(p_tf_joints_resistive)",
                self.data.tfcoil.p_tf_joints_resistive,
            )
            po.ovarre(
                self.outfile,
                "Outboard leg resistance per coil (ohm)",
                "(res_tf_leg)",
                self.data.tfcoil.res_tf_leg,
            )
            po.ovarre(
                self.outfile,
                "Average CP temperature (K)",
                "(temp_cp_average)",
                self.data.tfcoil.temp_cp_average,
            )
            po.ovarre(
                self.outfile,
                "Average leg temperature (K)",
                "(temp_tf_legs_outboard)",
                self.data.tfcoil.temp_tf_legs_outboard,
            )

        else:
            po.ovarre(
                self.outfile,
                "TF resistivity (ohm.m)",
                "(p_cp_resistive)",
                self.data.tfcoil.rho_cp,
            )
            po.ovarre(
                self.outfile,
                "TF coil resistive power less (total) (ohm.m)",
                "(p_cp_resistive)",
                self.data.tfcoil.p_cp_resistive,
            )
            po.ovarre(
                self.outfile,
                "Average coil temperature (K)",
                "(temp_cp_average)",
                self.data.tfcoil.temp_cp_average,
            )

        # Top section TF coil radial build (physics_variables.itart = 1 only)
        if self.data.physics.itart == 1:
            po.osubhd(self.outfile, "Radial build of TF coil at central collumn top :")
            # write(self.outfile,5)

            # Restart the radial build at bucking cylindre inner radius
            radius = self.data.build.r_tf_inboard_in
            po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

            radius += self.data.tfcoil.dr_tf_nose_case
            po.obuild(
                self.outfile,
                "Coil bucking cylindre",
                self.data.tfcoil.dr_tf_nose_case,
                radius,
                "(dr_tf_nose_case)",
            )

            radius += self.data.tfcoil.dx_tf_wp_insulation
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                self.data.tfcoil.dx_tf_wp_insulation,
                radius,
                "(dx_tf_wp_insulation)",
            )

            radius = (
                radius
                + 0.5e0 * self.data.superconducting_tfcoil.dr_tf_wp_top
                - self.data.tfcoil.dx_tf_wp_insulation
            )
            po.obuild(
                self.outfile,
                "Conductor - first half",
                0.5e0 * self.data.superconducting_tfcoil.dr_tf_wp_top
                - self.data.tfcoil.dx_tf_wp_insulation,
                radius,
                "(dr_tf_wp_top/2-dx_tf_wp_insulation)",
            )

            radius = (
                radius
                + 0.5e0 * self.data.superconducting_tfcoil.dr_tf_wp_top
                - self.data.tfcoil.dx_tf_wp_insulation
            )
            po.obuild(
                self.outfile,
                "Conductor - second half",
                0.5e0 * self.data.superconducting_tfcoil.dr_tf_wp_top
                - self.data.tfcoil.dx_tf_wp_insulation,
                radius,
                "(dr_tf_wp_top/2-dx_tf_wp_insulation)",
            )

            radius += self.data.tfcoil.dx_tf_wp_insulation
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                self.data.tfcoil.dx_tf_wp_insulation,
                radius,
                "(dx_tf_wp_insulation)",
            )

            radius += self.data.tfcoil.dr_tf_plasma_case
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                self.data.tfcoil.dr_tf_plasma_case,
                radius,
                "(dr_tf_plasma_case)",
            )

            # Consistency check
            if abs(radius - self.data.build.r_cp_top) < np.finfo(float(radius)).eps:
                po.ocmmnt(self.outfile, "Top TF coil dimensions are consistent")
            else:
                po.ocmmnt(self.outfile, "ERROR: TF coil dimensions are NOT consistent:")
                po.ovarre(
                    self.outfile,
                    "Radius of plasma-facing side of inner leg SHOULD BE [m]",
                    "",
                    self.data.build.r_cp_top,
                )
                po.oblnkl(self.outfile)

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
        Calculates the volume and resistive power losses of a TART centrepost
        This routine calculates the volume and resistive power losses
        of a TART centrepost. It is assumed to be tapered - narrowest at
        the midplane and reaching maximum thickness at the height of the
        plasma. Above/below the plasma, the centrepost is cylindrical.
        The shape of the taper is assumed to be an arc of a circle.

        F/MI/PJK/LOGBOOK12, pp.33,34

        Parameters
        ----------
        r_tf_inboard_in :

        r_tf_inboard_out :

        r_cp_top :

        ztop :

        hmaxi :

        cas_in_th :

        cas_out_th :

        gr_ins_th :

        ins_th :

        n_tf_coil_turns :

        curr :

        rho :

        fcool :

        n_tf_coils :

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

            # Volumes and resistive losses calculations
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
            sum1 += yy_cond[ii]
            sum2 += 1.0e0 / yy_cond[ii]
            sum3 += yy_ins[ii]
            sum4 += yy_casout[ii]
            sum5 += yy_gr_ins[ii]

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
    """Copper TF coil class for resistive TF coil calculations.
    Inherits from ResistiveTFCoil and implements specific methods for copper TF coils.
    """


class AluminiumTFCoil(ResistiveTFCoil):
    """Aluminium TF coil class for resistive TF coil calculations.
    Inherits from ResistiveTFCoil and implements specific methods for aluminium TF coils.
    """
