import copy
import json

import numba
import numpy as np

from process import fortran as ft
from process import process_output as po
from process.build import Build
from process.exceptions import ProcessValueError
from process.fortran import (
    build_variables,
    constants,
    error_handling,
    fwbs_variables,
    global_variables,
    numerics,
    pfcoil_variables,
    physics_variables,
    rebco_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.fortran import build_variables as bv
from process.fortran import error_handling as eh
from process.fortran import fwbs_variables as fwbsv
from process.fortran import tfcoil_variables as tfv
from process.utilities.f2py_string_patch import (
    f2py_compatible_to_string,
    string_to_f2py_compatible,
)

RMU0 = constants.rmu0


class TFCoil:
    """Calculates the parameters of a resistive TF coil system for a fusion power plant"""

    def __init__(self, build: Build):
        """Initialise Fortran module variables."""
        self.outfile = ft.constants.nout  # output file unit
        self.iprint = 0  # switch for writing to output file (1=yes)
        self.build = build
        self.a_tf_coil_inboard = tfcoil_variables.a_tf_coil_inboard

    def run(self, output):
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
        ) = self.tf_global_geometry(
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

        (
            tfcoil_variables.b_tf_inboard_peak,
            tfcoil_variables.c_tf_total,
            sctfcoil_module.c_tf_coil,
            tfcoil_variables.oacdcp,
        ) = self.tf_current(
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
        ) = self.tf_coil_shape_inner(
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

        if physics_variables.itart == 0 and tfcoil_variables.i_tf_shape == 1:
            tfcoil_variables.ind_tf_coil = self.tfcind(
                build_variables.dr_tf_inboard,
                tfcoil_variables.r_tf_arc,
                tfcoil_variables.z_tf_arc,
            )
        else:
            tfcoil_variables.ind_tf_coil = (
                (build_variables.z_tf_inside_half + build_variables.dr_tf_outboard)
                * RMU0
                / constants.pi
                * np.log(
                    build_variables.r_tf_outboard_mid / build_variables.r_tf_inboard_mid
                )
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
                sctfcoil_module.dr_tf_turn_cable_space,
                tfcoil_variables.dia_tf_turn_coolant_channel,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.dx_tf_turn_steel,
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
                tfcoil_variables.dr_tf_plasma_case,
                tfcoil_variables.i_tf_stress_model,
                sctfcoil_module.vforce_inboard_tot,
                tfcoil_variables.i_tf_tresca,
                tfcoil_variables.acasetf,
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

    def output(self):
        """Run main tfcoil subroutine and write output."""
        self.iprint = 1
        self.tfcoil(output=bool(self.iprint))

    def tf_global_geometry(
        self,
        i_tf_case_geom: int,
        i_f_dr_tf_plasma_case: bool,
        f_dr_tf_plasma_case: float,
        tfc_sidewall_is_fraction: bool,
        casths_fraction: float,
        n_tf_coils: int,
        dr_tf_inboard: float,
        dr_tf_nose_case: float,
        r_tf_inboard_out: float,
        r_tf_inboard_in: float,
        r_tf_outboard_mid: float,
        dr_tf_outboard: float,
    ) -> tuple[float, float, float, float, float, float, float, float, float]:
        """
        Calculate the global geometry of the Toroidal Field (TF) coil.

        This method computes the overall geometry of the TF coil, including:

        - Toroidal angular spacing and tangent of the angle.
        - Cross-sectional areas of the inboard and outboard legs.
        - Radii and widths of the inboard and outboard legs.
        - Case thicknesses for plasma-facing and sidewall regions.

        :param i_tf_case_geom:
            Geometry type of the TF coil case (e.g., circular or straight plasma-facing front case).
        :type i_tf_case_geom: int
        :param i_f_dr_tf_plasma_case:
            Whether the plasma-facing case thickness is specified as a fraction of the inboard thickness.
        :type i_f_dr_tf_plasma_case: bool
        :param f_dr_tf_plasma_case:
            Fraction of the inboard thickness used for the plasma-facing case thickness.
        :type f_dr_tf_plasma_case: float
        :param tfc_sidewall_is_fraction:
            Whether the sidewall case thickness is specified as a fraction of the inboard radius.
        :type tfc_sidewall_is_fraction: bool
        :param casths_fraction:
            Fraction of the inboard radius used for the sidewall case thickness.
        :type casths_fraction: float
        :param n_tf_coils:
            Number of TF coils.
        :type n_tf_coils: int
        :param dr_tf_inboard:
            Radial thickness of the inboard leg of the TF coil [m].
        :type dr_tf_inboard: float
        :param dr_tf_nose_case:
            Thickness of the inboard leg case at the nose region [m].
        :type dr_tf_nose_case: float
        :param r_tf_inboard_out:
            Outer radius of the inboard leg of the TF coil [m].
        :type r_tf_inboard_out: float
        :param r_tf_inboard_in:
            Inner radius of the inboard leg of the TF coil [m].
        :type r_tf_inboard_in: float
        :param r_tf_outboard_mid:
            Mid-plane radius of the outboard leg of the TF coil [m].
        :type r_tf_outboard_mid: float
        :param dr_tf_outboard:
            Radial thickness of the outboard leg of the TF coil [m].
        :type dr_tf_outboard: float

        :returns:
            A tuple containing:
            - **rad_tf_coil_toroidal** (*float*): Toroidal angular spacing of each TF coil [radians].
            - **tan_theta_coil** (*float*): Tangent of the toroidal angular spacing.
            - **a_tf_coil_inboard** (*float*): Cross-sectional area of the inboard leg of the TF coil [m²].
            - **r_tf_outboard_in** (*float*): Inner radius of the outboard leg of the TF coil [m].
            - **r_tf_outboard_out** (*float*): Outer radius of the outboard leg of the TF coil [m].
            - **dx_tf_inboard_out_toroidal** (*float*): Width of the inboard leg at the outer edge in the toroidal direction [m].
            - **a_tf_leg_outboard** (*float*): Cross-sectional area of the outboard leg of the TF coil [m²].
        :rtype: tuple
        """

        # The angular space of each TF coil in the toroidal direction [rad]
        # The angular space between each coil is the same as that of the coils

        rad_tf_coil_toroidal = np.pi / n_tf_coils
        tan_theta_coil = np.tan(rad_tf_coil_toroidal)

        # TF coil inboard legs total mid-plane cross-section area [m^2]
        if i_tf_case_geom == 0:
            # Circular plasma facing front case
            a_tf_coil_inboard = np.pi * (r_tf_inboard_out**2 - r_tf_inboard_in**2)
        else:
            # Straight plasma facing front case
            a_tf_coil_inboard = (
                n_tf_coils
                * np.sin(rad_tf_coil_toroidal)
                * np.cos(rad_tf_coil_toroidal)
                * r_tf_inboard_out**2
                - np.pi * r_tf_inboard_in**2
            )

        # TF coil width in toroidal direction at inboard leg outer edge [m]

        dx_tf_inboard_out_toroidal = (
            2.0e0 * r_tf_inboard_out * np.sin(rad_tf_coil_toroidal)
        )

        # Outer leg geometry
        # Mid-plane inner/out radial position of the TF coil outer leg [m]

        r_tf_outboard_in = r_tf_outboard_mid - (dr_tf_outboard * 0.5e0)
        r_tf_outboard_out = r_tf_outboard_mid + (dr_tf_outboard * 0.5e0)

        # Area of rectangular cross-section TF outboard leg [m^2]
        a_tf_leg_outboard = dx_tf_inboard_out_toroidal * dr_tf_outboard

        # Plasma facing front case thickness [m]
        if i_f_dr_tf_plasma_case:
            dr_tf_plasma_case = f_dr_tf_plasma_case * dr_tf_inboard
        else:
            dr_tf_plasma_case = tfcoil_variables.dr_tf_plasma_case

        # Case thickness of side wall [m]
        if tfc_sidewall_is_fraction:
            dx_tf_side_case = (
                casths_fraction
                * (r_tf_inboard_in + dr_tf_nose_case)
                * np.tan(np.pi / n_tf_coils)
            )
        else:
            dx_tf_side_case = tfcoil_variables.dx_tf_side_case

        return (
            rad_tf_coil_toroidal,
            tan_theta_coil,
            a_tf_coil_inboard,
            r_tf_outboard_in,
            r_tf_outboard_out,
            dx_tf_inboard_out_toroidal,
            a_tf_leg_outboard,
            dr_tf_plasma_case,
            dx_tf_side_case,
        )

    def tf_current(
        self,
        n_tf_coils: int,
        bt: float,
        rmajor: float,
        r_b_tf_inboard_peak: float,
        a_tf_coil_inboard: float,
    ) -> tuple[float, float, float, float]:
        """
        Calculate the maximum B field and the corresponding TF current.

        :param n_tf_coils: Number of TF coils.
        :type n_tf_coils: int
        :param bt: Toroidal magnetic field at the plasma center [T].
        :type bt: float
        :param rmajor: Major radius of the plasma [m].
        :type rmajor: float
        :param r_b_tf_inboard_peak: Radius at which the peak inboard B field occurs [m].
        :type r_b_tf_inboard_peak: float
        :param a_tf_coil_inboard: Cross-sectional area of the inboard leg of the TF coil [m²].
        :type a_tf_coil_inboard: float

        :returns: A tuple containing:
            - **b_tf_inboard_peak** (*float*): Maximum B field on the magnet [T].
            - **c_tf_total** (*float*): Total current in TF coils [A].
            - **c_tf_coil** (*float*): Current per TF coil [A].
            - **oacdcp** (*float*): Global inboard leg average current density in TF coils [A/m²].
        :rtype: tuple[float, float, float, float]
        """

        # Calculation of the maximum B field on the magnet [T]
        b_tf_inboard_peak = bt * rmajor / r_b_tf_inboard_peak

        # Total current in TF coils [A]
        c_tf_total = (
            b_tf_inboard_peak * r_b_tf_inboard_peak * (2 * np.pi / constants.rmu0)
        )

        # Current per TF coil [A]
        c_tf_coil = c_tf_total / n_tf_coils

        # Global inboard leg average current in TF coils [A/m2]
        oacdcp = c_tf_total / a_tf_coil_inboard

        return b_tf_inboard_peak, c_tf_total, c_tf_coil, oacdcp

    def tf_coil_shape_inner(
        self,
        i_tf_shape: int,
        itart: int,
        i_single_null: int,
        r_tf_inboard_out: float,
        r_cp_top: float,
        rmajor: float,
        rminor: float,
        r_tf_outboard_in: float,
        z_tf_inside_half: float,
        z_tf_top: float,
        dr_tf_inboard: float,
        dr_tf_outboard: float,
        r_tf_outboard_mid: float,
        r_tf_inboard_mid: float,
    ) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate the shape of the inside of the TF coil.

        This method approximates the TF coil by a straight inboard section and four elliptical arcs.
        The model is ad hoc and does not have a physics or engineering basis.

        :param i_tf_shape: TF coil shape selection switch.
        :type i_tf_shape: int
        :param itart: TART (tight aspect ratio tokamak) switch.
        :type itart: int
        :param i_single_null: Single null switch.
        :type i_single_null: int
        :param r_tf_inboard_out: Outer radius of the inboard leg of the TF coil [m].
        :type r_tf_inboard_out: float
        :param r_cp_top: Centrepost outer radius at the top [m].
        :type r_cp_top: float
        :param rmajor: Major radius of the plasma [m].
        :type rmajor: float
        :param rminor: Minor radius of the plasma [m].
        :type rminor: float
        :param r_tf_outboard_in: Inner radius of the outboard leg of the TF coil [m].
        :type r_tf_outboard_in: float
        :param z_tf_inside_half: Maximum inboard edge height [m].
        :type z_tf_inside_half: float
        :param z_tf_top: Vertical position of the top of the TF coil [m].
        :type z_tf_top: float
        :param dr_tf_inboard: Radial thickness of the inboard leg of the TF coil [m].
        :type dr_tf_inboard: float
        :param dr_tf_outboard: Radial thickness of the outboard leg of the TF coil [m].
        :type dr_tf_outboard: float
        :param r_tf_outboard_mid: Mid-plane radius of the outboard leg of the TF coil [m].
        :type r_tf_outboard_mid: float
        :param r_tf_inboard_mid: Mid-plane radius of the inboard leg of the TF coil [m].
        :type r_tf_inboard_mid: float

        :returns:
            - len_tf_coil (float): Total length of the TF coil.
            - tfa (np.ndarray): Arc segment lengths in R.
            - tfb (np.ndarray): Arc segment lengths in Z.
            - r_tf_arc (np.ndarray): R coordinates of arc points.
            - z_tf_arc (np.ndarray): Z coordinates of arc points.
        :rtype: tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        """
        FSTRAIGHT = 0.6
        len_tf_coil = 0.0
        r_tf_arc = np.zeros(5)
        z_tf_arc = np.zeros(5)
        tfa = np.zeros(4)
        tfb = np.zeros(4)

        if i_tf_shape == 1 and itart == 0:
            # PROCESS D-shape parametrisation
            r_tf_arc[0] = r_tf_inboard_out
            r_tf_arc[1] = rmajor - 0.2e0 * rminor
            r_tf_arc[2] = r_tf_outboard_in
            r_tf_arc[3] = r_tf_arc[1]
            r_tf_arc[4] = r_tf_arc[0]

            if i_single_null == 0:
                z_tf_arc[0] = FSTRAIGHT * (z_tf_top - dr_tf_inboard)
                z_tf_arc[1] = z_tf_top - dr_tf_inboard
                z_tf_arc[2] = 0
                z_tf_arc[3] = -z_tf_inside_half
                z_tf_arc[4] = -FSTRAIGHT * z_tf_inside_half
            else:
                z_tf_arc[0] = FSTRAIGHT * z_tf_inside_half
                z_tf_arc[1] = z_tf_inside_half
                z_tf_arc[2] = 0
                z_tf_arc[3] = -z_tf_inside_half
                z_tf_arc[4] = -FSTRAIGHT * z_tf_inside_half

            len_tf_coil = z_tf_arc[0] - z_tf_arc[4]

            for ii in range(4):
                tfa[ii] = abs(r_tf_arc[ii + 1] - r_tf_arc[ii])
                tfb[ii] = abs(z_tf_arc[ii + 1] - z_tf_arc[ii])
                aa = tfa[ii] + 0.5e0 * dr_tf_inboard
                bb = tfb[ii] + 0.5e0 * dr_tf_inboard
                len_tf_coil += 0.25e0 * self.circumference(aa, bb)

        elif i_tf_shape == 1 and itart == 1:
            # Centrepost with D-shaped
            r_tf_arc[0] = r_cp_top
            r_tf_arc[1] = rmajor - 0.2e0 * rminor
            r_tf_arc[2] = r_tf_outboard_in
            r_tf_arc[3] = r_tf_arc[1]
            r_tf_arc[4] = r_tf_arc[0]

            z_tf_arc[0] = z_tf_top - dr_tf_inboard
            z_tf_arc[1] = z_tf_top - dr_tf_inboard
            z_tf_arc[2] = 0
            z_tf_arc[3] = -z_tf_inside_half
            z_tf_arc[4] = -z_tf_inside_half

            len_tf_coil = 2 * (r_tf_arc[1] - r_tf_arc[0])

            for ii in range(1, 3):
                tfa[ii] = abs(r_tf_arc[ii + 1] - r_tf_arc[ii])
                tfb[ii] = abs(z_tf_arc[ii + 1] - z_tf_arc[ii])
                aa = tfa[ii] + 0.5e0 * dr_tf_outboard
                bb = tfb[ii] + 0.5e0 * dr_tf_outboard
                len_tf_coil += 0.25e0 * self.circumference(aa, bb)

        elif i_tf_shape == 2:
            # Picture frame coil
            if itart == 0:
                r_tf_arc[0] = r_tf_inboard_out
            if itart == 1:
                r_tf_arc[0] = r_cp_top
            r_tf_arc[1] = r_tf_outboard_in
            r_tf_arc[2] = r_tf_arc[1]
            r_tf_arc[3] = r_tf_arc[1]
            r_tf_arc[4] = r_tf_arc[0]

            z_tf_arc[0] = z_tf_top - dr_tf_inboard
            z_tf_arc[1] = z_tf_top - dr_tf_inboard
            z_tf_arc[2] = 0
            z_tf_arc[3] = -z_tf_inside_half
            z_tf_arc[4] = -z_tf_inside_half

            if itart == 0:
                len_tf_coil = 2.0e0 * (
                    2.0e0 * z_tf_inside_half
                    + dr_tf_inboard
                    + r_tf_outboard_mid
                    - r_tf_inboard_mid
                )
            elif itart == 1:
                len_tf_coil = (
                    z_tf_inside_half + z_tf_top + 2.0e0 * (r_tf_outboard_mid - r_cp_top)
                )

        return len_tf_coil, tfa, tfb, r_tf_arc, z_tf_arc

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
            "(z_tf_inside_half)",
            build_variables.z_tf_inside_half,
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
                f"  {ii}              {tfcoil_variables.r_tf_arc[ii]}              {tfcoil_variables.z_tf_arc[ii]}",
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} R (m)",
                f"(r_tf_arc({ii + 1}))",
                tfcoil_variables.r_tf_arc[ii],
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} Z (m)",
                f"(z_tf_arc({ii + 1}))",
                tfcoil_variables.z_tf_arc[ii],
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
                "(z_tf_inside_half + dr_tf_outboard)",
                build_variables.z_tf_inside_half + build_variables.dr_tf_outboard,
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
                "(dr_tf_plasma_case)",
                tfcoil_variables.dr_tf_plasma_case,
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
                    "(dr_tf_turn)",
                    sctfcoil_module.dr_tf_turn,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of turn (m)",
                    "(dx_tf_turn)",
                    sctfcoil_module.dx_tf_turn,
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
                    "(dr_tf_turn_cable_space)",
                    sctfcoil_module.dr_tf_turn_cable_space,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of cable space",
                    "(dx_tf_turn_cable_space)",
                    sctfcoil_module.dx_tf_turn_cable_space,
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
                "(dx_tf_turn_steel)",
                tfcoil_variables.dx_tf_turn_steel,
            )
            po.ovarre(
                self.outfile,
                "Inter-turn insulation thickness (m)",
                "(dx_tf_turn_insulation)",
                tfcoil_variables.dx_tf_turn_insulation,
            )

            if tfcoil_variables.i_tf_sc_mat in (1, 2, 3, 4, 5, 7, 8, 9):
                po.osubhd(self.outfile, "Conductor information:")
                po.ovarre(
                    self.outfile,
                    "Diameter of central helium channel in cable",
                    "(dia_tf_turn_coolant_channel)",
                    tfcoil_variables.dia_tf_turn_coolant_channel,
                )
                po.ocmmnt(self.outfile, "Fractions by area")
                po.ovarre(
                    self.outfile,
                    "internal area of the cable space",
                    "(a_tf_turn_cable_space)",
                    tfcoil_variables.a_tf_turn_cable_space,
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
                # po.ovarre(self.outfile,'Conduit fraction of winding pack','(tfcoil_variables.n_tf_turn*tfcoil_variables.a_tf_turn_steel/ap)',n_tf_turn*tfcoil_variables.a_tf_turn_steel/ap, 'OP ')
                # po.ovarre(self.outfile,'Insulator fraction of winding pack','(tfcoil_variables.a_tf_coil_wp_turn_insulation/ap)',a_tf_coil_wp_turn_insulation/ap, 'OP ')
                # po.ovarre(self.outfile,'Helium area fraction of winding pack excluding central channel','(tfcoil_variables.avwp/ap)',avwp/ap, 'OP ')
                # po.ovarre(self.outfile,'Central helium channel area as fraction of winding pack','(tfcoil_variables.a_tf_wp_coolant_channels/ap)',a_tf_wp_coolant_channels/ap, 'OP ')
                ap = (
                    tfcoil_variables.acond
                    + tfcoil_variables.n_tf_turn * tfcoil_variables.a_tf_turn_steel
                    + tfcoil_variables.a_tf_coil_wp_turn_insulation
                    + tfcoil_variables.avwp
                    + tfcoil_variables.a_tf_wp_coolant_channels
                )
                po.ovarrf(
                    self.outfile,
                    "Check total area fractions in winding pack = 1",
                    "",
                    (
                        tfcoil_variables.acond
                        + tfcoil_variables.n_tf_turn * tfcoil_variables.a_tf_turn_steel
                        + tfcoil_variables.a_tf_coil_wp_turn_insulation
                        + tfcoil_variables.avwp
                        + tfcoil_variables.a_tf_wp_coolant_channels
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
                "(dr_tf_plasma_case)",
                tfcoil_variables.dr_tf_plasma_case,
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
                "(dx_tf_turn_insulation)",
                tfcoil_variables.dx_tf_turn_insulation,
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
                "(c_tf_turn)",
                tfcoil_variables.c_tf_turn,
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
            "Radius of maximum field position on inboard TF coil (m)",
            "(r_b_tf_inboard_peak)",
            tfcoil_variables.r_b_tf_inboard_peak,
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
                    "(p_tf_joints_resistive)",
                    tfcoil_variables.p_tf_joints_resistive,
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

            radius = radius + tfcoil_variables.dr_tf_plasma_case
            po.obuild(
                self.outfile,
                "Plasma side case min radius",
                tfcoil_variables.dr_tf_plasma_case,
                radius,
                "(dr_tf_plasma_case)",
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

            radius = radius + tfcoil_variables.dr_tf_plasma_case
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.dr_tf_plasma_case,
                radius,
                "(dr_tf_plasma_case)",
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

            radius = radius + tfcoil_variables.dr_tf_plasma_case
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.dr_tf_plasma_case,
                radius,
                "(dr_tf_plasma_case)",
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

    @staticmethod
    def circumference(aaa, bbb):
        """Calculate ellipse arc circumference using Ramanujan approximation (m)
        See https://www.johndcook.com/blog/2013/05/05/ramanujan-circumference-ellipse/
        for a discussion of the precision of the formula

        An ellipse has the following formula: (x/a)² + (y/b)² = 1

        :param aaa: the value of a in the formula of the ellipse.
        :type aaa: float

        :param bbb: the value of b in the formula of the ellipse.
        :type bbb: float

        :returns: an approximation of the circumference of the ellipse
        :rtype: float
        """
        hh = (aaa - bbb) ** 2 / (aaa + bbb) ** 2
        return (
            np.pi
            * (aaa + bbb)
            * (1.0e0 + (3.0e0 * hh) / (10.0e0 + np.sqrt(4.0e0 - 3.0e0 * hh)))
        )

    def cntrpst(self):
        """
        Evaluates the properties of a TART centrepost
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This subroutine evaluates the parameters of the centrepost for a
        tight aspect ratio tokamak. The centrepost is assumed to be tapered,
        i.e. narrowest on the midplane (z=0).
        """

        # Vertical distance from the midplane to the top of the tapered section [m]
        if physics_variables.itart == 1:
            sctfcoil_module.z_cp_top = (
                build_variables.z_plasma_xpoint_upper + tfcoil_variables.dztop
            )

        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            tfcoil_variables.dx_tf_inboard_out_toroidal = (
                2.0e0
                * build_variables.r_cp_top
                * np.sin(sctfcoil_module.rad_tf_coil_toroidal)
            )

        # Temperature margin used in calculations (K)
        tmarg = 10.0e0

        # Number of integral step used for the coolant temperature rise
        n_tcool_it = 20

        # Coolant channels:
        acool = tfv.a_cp_cool * tfv.n_tf_coils  # Cooling cross-sectional area
        dcool = 2.0e0 * tfv.rcool  # Diameter
        lcool = 2.0e0 * (bv.z_tf_inside_half + bv.dr_tf_outboard)  # Length
        tfv.ncool = acool / (constants.pi * tfv.rcool**2)  # Number

        # Average conductor cross-sectional area to cool (with cooling area)
        acpav = (
            0.5e0 * tfv.vol_cond_cp / (bv.z_tf_inside_half + bv.dr_tf_outboard) + acool
        )
        ro = (acpav / (constants.pi * tfv.ncool)) ** 0.5

        # Inner legs total heating power (to be removed by coolant)
        ptot = tfv.p_cp_resistive + fwbsv.pnuc_cp_tf * 1.0e6

        # Temperature calculations
        # -------------------------
        # Temperature rise in coolant (inlet to outlet)
        # **********************************************
        # Water coollant
        # --------------
        if tfv.i_tf_sup == 0:
            # Water coolant physical properties
            coolant_density = constants.denh2o
            coolant_cp = constants.cph2o
            coolant_visco = constants.muh2o
            coolant_th_cond = constants.kh2o

            # Mass flow rate [kg/s]
            cool_mass_flow = acool * coolant_density * tfv.vcool

            # Water temperature rise
            tfv.dtiocool = ptot / (cool_mass_flow * coolant_cp)

            # Constant coolant velocity
            vcool_max = tfv.vcool
            # --------------

        # Helium coolant
        # --------------
        elif tfv.i_tf_sup == 2:
            # Inlet coolant density [kg/m3]
            coolant_density = self.he_density(tfv.tcoolin)

            # Mass flow rate [kg/s]
            cool_mass_flow = acool * coolant_density * tfv.vcool

            # Infinitesimal power deposition used in the integral
            dptot = ptot / n_tcool_it

            tcool_calc = copy.copy(tfv.tcoolin)  # K
            for _i in range(n_tcool_it):
                # Thermal capacity Cp
                coolant_cp = self.he_cp(tcool_calc)

                # Temperature infinitesimal increase
                tcool_calc += dptot / (cool_mass_flow * coolant_cp)

            # Outlet coolant density (minimal coolant density value)
            coolant_density = self.he_density(tcool_calc)

            # Maxium coolant velocity
            vcool_max = cool_mass_flow / (acool * coolant_density)

            # Getting the global in-outlet temperature increase
            tfv.dtiocool = tcool_calc - tfv.tcoolin
        # --------------

        # Average coolant temperature
        tcool_av = tfv.tcoolin + 0.5e0 * tfv.dtiocool
        # **********************************************

        # Film temperature rise
        # *********************
        # Rem : The helium cooling properties are calculated using the outlet ones
        # this is not an exact approximation for average temperature rise

        # Helium viscosity
        if tfv.i_tf_sup == 2:
            coolant_visco = self.he_visco(tcool_av)

        # Reynolds number
        reyn = coolant_density * tfv.vcool * dcool / coolant_visco

        # Helium thermal conductivity [W/(m.K)]
        if tfv.i_tf_sup == 2:
            coolant_th_cond = self.he_th_cond(tcool_av)

        # Prandlt number
        prndtl = coolant_cp * coolant_visco / coolant_th_cond

        # Film temperature difference calculations
        # Originally prandtl was prndtl**0.3e0 but this is incorrect as from
        # Dittus-Boelter correlation where the fluid is being heated it should be as below
        nuselt = 0.023e0 * reyn**0.8e0 * prndtl**0.4e0
        h = nuselt * coolant_th_cond / dcool
        dtfilmav = ptot / (h * 2.0e0 * constants.pi * tfv.rcool * tfv.ncool * lcool)

        # Average film temperature (in contact with te conductor)
        tcool_film = tcool_av + dtfilmav
        # *********************

        # Temperature rise in conductor
        # ------------------------------
        # Conductor thermal conductivity
        # ******
        # Copper conductor
        if tfv.i_tf_sup == 0:
            conductor_th_cond = constants.k_copper

        # Aluminium
        elif tfv.i_tf_sup == 2:
            conductor_th_cond = self.al_th_cond(tcool_film)
        # ******

        # Average temperature rise : To be changed with Garry Voss' better documented formula ?
        dtcncpav = (
            (ptot / tfv.vol_cond_cp)
            / (2.0e0 * conductor_th_cond * (ro**2 - tfv.rcool**2))
            * (
                ro**2 * tfv.rcool**2
                - 0.25e0 * tfv.rcool**4
                - 0.75e0 * ro**4
                + ro**4 * np.log(ro / tfv.rcool)
            )
        )

        # Peak temperature rise : To be changed with Garry Voss' better documented formula ?
        dtconcpmx = (
            (ptot / tfv.vol_cond_cp)
            / (2.0e0 * conductor_th_cond)
            * ((tfv.rcool**2 - ro**2) / 2.0e0 + ro**2 * np.log(ro / tfv.rcool))
        )

        # If the average conductor temperature difference is negative, set it to 0
        if dtcncpav < 0.0e0:
            eh.report_error(249)
            dtcncpav = 0.0e0

        # If the average conductor temperature difference is negative, set it to 0
        if dtconcpmx < 0.0e0:
            eh.report_error(250)
            dtconcpmx = 0.0e0

        # Average conductor temperature
        tfv.tcpav2 = tfv.tcoolin + dtcncpav + dtfilmav + 0.5e0 * tfv.dtiocool

        # Peak wall temperature
        tfv.tcpmax = tfv.tcoolin + tfv.dtiocool + dtfilmav + dtconcpmx
        tcoolmx = tfv.tcoolin + tfv.dtiocool + dtfilmav
        # -------------------------

        # Thermal hydraulics: friction factor from Z. Olujic, Chemical
        # Engineering, Dec. 1981, p. 91
        roughrat = 4.6e-5 / dcool
        fricfac = (
            1.0e0
            / (
                -2.0e0
                * np.log10(
                    roughrat / 3.7e0
                    - 5.02e0 / reyn * np.log10(roughrat / 3.7e0 + 14.5e0 / reyn)
                )
            )
            ** 2
        )

        # Pumping efficiency
        if tfv.i_tf_sup == 0:  # Water cooled
            tfv.etapump = 0.8e0
        elif tfv.i_tf_sup == 2:  # Cryogenic helium
            tfv.etapump = 0.6e0

        # Pressure drop calculation
        dpres = fricfac * (lcool / dcool) * coolant_density * 0.5e0 * tfv.vcool**2
        tfv.p_cp_coolant_pump_elec = dpres * acool * tfv.vcool / tfv.etapump

        # Critical pressure in saturation pressure calculations (Pa)
        pcrt = 2.24e7

        # Saturation pressure
        # Ref : Keenan, Keyes, Hill, Moore, steam tables, Wiley & Sons, 1969
        # Rem 1 : ONLY VALID FOR WATER !
        # Rem 2 : Not used anywhere else in the code ...
        tclmx = tcoolmx + tmarg
        tclmxs = min(tclmx, 374.0e0)
        fc = 0.65e0 - 0.01e0 * tclmxs
        sum_ = (
            -741.9242e0
            - 29.721e0 * fc
            - 11.55286e0 * fc**2
            - 0.8685635e0 * fc**3
            + 0.1094098e0 * fc**4
            + 0.439993e0 * fc**5
            + 0.2520658e0 * fc**6
            + 0.0518684e0 * fc**7
        )
        psat = pcrt * np.exp(0.01e0 / (tclmxs + 273.0e0) * (374.0e0 - tclmxs) * sum_)
        presin = psat + dpres

        # Output section
        if self.iprint == 1:
            po.oheadr(self.outfile, "Centrepost Coolant Parameters")
            po.ovarre(
                self.outfile, "Centrepost coolant fraction", "(fcoolcp)", tfv.fcoolcp
            )
            po.ovarre(
                self.outfile, "Average coolant channel diameter (m)", "(dcool)", dcool
            )
            po.ovarre(self.outfile, "Coolant channel length (m)", "(lcool)", lcool)
            po.ovarre(
                self.outfile, "Inlet coolant flow speed (m/s)", "(vcool)", tfv.vcool
            )
            po.ovarre(
                self.outfile,
                "Outlet coolant flow speed (m/s)",
                "(vcool_max)",
                vcool_max,
            )
            po.ovarre(
                self.outfile,
                "Coolant mass flow rate (kg/s)",
                "(cool_mass_flow)",
                cool_mass_flow,
            )
            po.ovarre(self.outfile, "Number of coolant tubes", "(ncool)", tfv.ncool)
            po.ovarre(self.outfile, "Reynolds number", "(reyn)", reyn)
            po.ovarre(self.outfile, "Prandtl number", "(prndtl)", prndtl)
            po.ovarre(self.outfile, "Nusselt number", "(nuselt)", nuselt)

            po.osubhd(self.outfile, "Resistive Heating :")
            po.ovarre(
                self.outfile,
                "Average conductor resistivity (ohm.m)",
                "(rho_cp)",
                tfv.rho_cp,
            )
            po.ovarre(
                self.outfile,
                "Resistive heating (MW)",
                "(p_cp_resistive/1.0e6)",
                tfv.p_cp_resistive / 1.0e6,
            )
            po.ovarre(
                self.outfile, "Nuclear heating (MW)", "(pnuc_cp_tf)", fwbsv.pnuc_cp_tf
            )
            po.ovarre(self.outfile, "Total heating (MW)", "(ptot/1.0e6)", ptot / 1.0e6)

            po.osubhd(self.outfile, "Temperatures :")
            po.ovarre(
                self.outfile,
                "Input coolant temperature (K)",
                "(tfv.tcoolin)",
                tfv.tcoolin,
            )
            po.ovarre(
                self.outfile,
                "Input-output coolant temperature rise (K)",
                "(dtiocool)",
                tfv.dtiocool,
            )
            po.ovarre(self.outfile, "Film temperature rise (K)", "(dtfilmav)", dtfilmav)
            po.ovarre(
                self.outfile,
                "Average temp gradient in conductor (K/m)",
                "(dtcncpav)",
                dtcncpav,
            )
            po.ovarre(
                self.outfile,
                "Average centrepost temperature (K)",
                "(tcpav2)",
                tfv.tcpav2,
            )
            po.ovarre(
                self.outfile, "Peak centrepost temperature (K)", "(tcpmax)", tfv.tcpmax
            )

            po.osubhd(self.outfile, "Pump Power :")
            po.ovarre(self.outfile, "Coolant pressure drop (Pa)", "(dpres)", dpres)
            if tfv.i_tf_sup == 0:  # Saturation pressure calculated with Water data ...
                po.ovarre(
                    self.outfile, "Coolant inlet pressure (Pa)", "(presin)", presin
                )

            po.ovarre(
                self.outfile,
                "Pump power (W)",
                "(p_cp_coolant_pump_elec)",
                tfv.p_cp_coolant_pump_elec,
            )

    def tf_field_and_force(self):
        """
        Calculate the TF coil field, force and VV quench consideration, and the resistive magnets resistance/volume
        """

        # Outer/inner WP radius removing the ground insulation layer and the insertion gap [m]
        if tfcoil_variables.i_tf_sup == 1:
            r_out_wp = (
                sctfcoil_module.r_wp_outer
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            r_in_wp = (
                sctfcoil_module.r_wp_inner
                + tfcoil_variables.tinstf
                + tfcoil_variables.tfinsgap
            )
        else:
            r_out_wp = sctfcoil_module.r_wp_outer - tfcoil_variables.tinstf
            r_in_wp = sctfcoil_module.r_wp_inner + tfcoil_variables.tinstf

        # Associated WP thickness
        dr_wp = r_out_wp - r_in_wp

        # In plane forces
        # ---
        # Centering force = net inwards radial force per meters per TF coil [N/m]
        tfcoil_variables.cforce = (
            0.5e0
            * tfcoil_variables.b_tf_inboard_peak
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
        )

        # Vertical force per coil [N]
        # ***
        # Rem : this force does not depends on the TF shape or the presence of
        #        sliding joints, the in/outboard vertical tension repartition is
        # -#
        # Ouboard leg WP plasma side radius without ground insulation/insertion gat [m]
        if tfcoil_variables.i_tf_sup == 1:
            r_in_outwp = (
                sctfcoil_module.r_tf_outboard_in
                + tfcoil_variables.dr_tf_plasma_case
                + tfcoil_variables.tinstf
                + tfcoil_variables.tfinsgap
            )
        else:
            r_in_outwp = sctfcoil_module.r_tf_outboard_in + tfcoil_variables.tinstf

        # If the TF coil has no dr_bore it would induce division by 0.
        # In this situation, the dr_bore radius is set to a very small value : 1.0e-9 m
        if abs(r_in_wp) < np.finfo(float(r_in_wp)).eps:
            r_in_wp = 1.0e-9

        # May the force be with you
        vforce_tot = (
            0.5e0
            * (
                physics_variables.bt
                * physics_variables.rmajor
                * tfcoil_variables.c_tf_total
            )
            / (tfcoil_variables.n_tf_coils * dr_wp**2)
            * (
                r_out_wp**2 * np.log(r_out_wp / r_in_wp)
                + r_in_outwp**2 * np.log((r_in_outwp + dr_wp) / r_in_outwp)
                + dr_wp**2 * np.log((r_in_outwp + dr_wp) / r_in_wp)
                - dr_wp * (r_out_wp + r_in_outwp)
                + 2.0e0
                * dr_wp
                * (
                    r_out_wp * np.log(r_in_wp / r_out_wp)
                    + r_in_outwp * np.log((r_in_outwp + dr_wp) / r_in_outwp)
                )
            )
        )

        # Case of a centrepost (physics_variables.itart == 1) with sliding joints (the CP vertical are separated from the leg ones)
        # Rem SK : casing/insulation thickness not subtracted as part of the CP is genuinely connected to the legs..
        if physics_variables.itart == 1 and tfcoil_variables.i_cp_joints == 1:
            # CP vertical tension [N]
            tfcoil_variables.vforce = (
                0.25e0
                * (
                    physics_variables.bt
                    * physics_variables.rmajor
                    * tfcoil_variables.c_tf_total
                )
                / (tfcoil_variables.n_tf_coils * dr_wp**2)
                * (
                    2.0e0 * r_out_wp**2 * np.log(r_out_wp / r_in_wp)
                    + 2.0e0 * dr_wp**2 * np.log(build_variables.r_cp_top / r_in_wp)
                    + 3.0e0 * dr_wp**2
                    - 2.0e0 * dr_wp * r_out_wp
                    + 4.0e0 * dr_wp * r_out_wp * np.log(r_in_wp / r_out_wp)
                )
            )

            # Vertical tension applied on the outer leg [N]
            tfcoil_variables.vforce_outboard = vforce_tot - tfcoil_variables.vforce

            # Inboard vertical tension fraction
            tfcoil_variables.f_vforce_inboard = tfcoil_variables.vforce / vforce_tot

        # Case of TF without joints or with clamped joints vertical tension
        else:
            # Inboard vertical tension [N]
            tfcoil_variables.vforce = tfcoil_variables.f_vforce_inboard * vforce_tot

            # Ouboard vertical tension [N]
            tfcoil_variables.vforce_outboard = tfcoil_variables.vforce * (
                (1.0e0 / tfcoil_variables.f_vforce_inboard) - 1.0e0
            )

        # ***

        # Total vertical force
        sctfcoil_module.vforce_inboard_tot = (
            tfcoil_variables.vforce * tfcoil_variables.n_tf_coils
        )

    @staticmethod
    def he_density(temp: float) -> float:
        """Author : S. Kahn
        Subroutine calculating temperature dependent helium density at 100 bar
        from fit using the following data, valid in [4-50] K
        Ref : R.D. McCarty, Adv. Cryo. Eng., 1990, 35, 1465-1475.

        :param temp: Helium temperature [K]
        :type temp: float

        :returns density: Heliyn density [kg/m3]
        :type density: float
        """

        # Fit range validation
        if temp < 4.0e0 or temp > 50.0e0:
            eh.fdiags[0] = temp
            eh.report_error(257)

        # Oder 3 polynomial fit
        if temp < 29.5e0:
            density = (
                217.753831e0
                - 1.66564525e0 * temp
                - 0.160654724e0 * temp**2
                + 0.00339003258e0 * temp**3
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif temp < 30.5e0:
            density = 231.40661479377616e0 - 3.917589985552496e0 * temp

        # Oder 2 polynomial fit
        else:
            density = 212.485251e0 - 4.18059786e0 * temp + 0.0289632937e0 * temp**2

        return density

    @staticmethod
    def he_cp(temp: float) -> float:
        """Author : S. Kahn
        Subroutine calculating temperature dependent thermal capacity at
        constant pressures at 100 Bar from fit using the following data
        valid in [4-50] K
        Ref : R.D. McCarty, Adv. Cryo. Eng., 1990, 35, 1465-1475.

        :param temp: Helium temperature [K]
        :type temp: float

        :return cp: Themal capacity at constant pressure [K/(kg.K)]
        :type cp: float
        """

        # Fit range validation
        if temp < 4.0e0 or temp > 50.0e0:
            eh.fdiags[0] = temp
            eh.report_error(257)

        # Order 3 polynomial fit in [4-30] K on the dimenion [K/(g.K)]
        if temp < 29.5e0:
            cp = (
                -0.834218557e0
                + 0.637079569e0 * temp
                - 0.0208839696e0 * temp**2
                + 0.000233433748e0 * temp**3
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif (
            temp < 30.5e0
        ):  # Linear interpolation between the fits to avoid discontinuity
            cp = 4.924018467550791e0 + 0.028953709588498633e0 * temp

        # Linear fit in [30-60] K on the dimenion [K/(g.K)]
        else:
            cp = 6.11883125e0 - 0.01022048e0 * temp

        # conversion to [K/(kg.K)] and return
        return cp * 1.0e3

    @staticmethod
    def he_visco(temp: float) -> float:
        """Author : S. Kahn
        Subroutine calculating temperature dependent He viscosity at 100 Bar
        from fit using the following data, valid in [4-50] K
        Ref : V.D. Arp,; R.D. McCarty ; Friend, D.G., Technical Note 1334, National
        Institute of Standards and Technology, Boulder, CO, 1998, 0.

        :param temp: Helium temperature [K]
        :type temp: float

        :return visco: Themal capacity at constant pressure [Pa.s]
        :type visco: float
        """

        if temp < 4.0e0 or temp > 50.0e0:
            eh.fdiags[0] = temp
            eh.report_error(257)

        # Order 4 polynomial exponential fit in [4-25] K
        if temp < 22.5e0:
            visco = np.exp(
                -9.19688182e0
                - 4.83007225e-1 * temp
                + 3.47720002e-2 * temp**2
                - 1.17501538e-3 * temp**3
                + 1.54218249e-5 * temp**4
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif temp < 27.5e0:
            visco = 6.708587487790973e-6 + 5.776427353055518e-9 * temp

        # Linear fit in [25-60] K
        else:
            visco = 5.41565319e-6 + 5.279222e-8 * temp

        return visco

    @staticmethod
    def he_th_cond(temp: float) -> float:
        """Author : S. Kahn
        Subroutine calculating temperature dependent He thermal conductivity
        at 100 Bar from fit using the following data, valid in [4-50] K
        Ref : B.A. Hands B.A., Cryogenics, 1981, 21, 12, 697-703.

        :param temp: Helium temperature [K]
        :type temp: float

        :return th_cond: Themal conductivity [W/(m.K)]
        :type th_cond: float
        """

        # Fit range validation
        if temp < 4.0e0 or temp > 50.0e0:
            eh.fdiags[0] = temp
            eh.report_error(257)

        # Order 4 polynomial fit
        if temp < 24.0e0:
            th_cond = (
                -7.56066334e-3
                + 1.62626819e-2 * temp
                - 1.3633619e-3 * temp**2
                + 4.84227752e-5 * temp**3
                - 6.31264281e-7 * temp**4
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif temp < 25.0e0:
            th_cond = 0.05858194642349288e0 - 5.706361831471496e-5 * temp

        # Order 2 polynomial fit
        elif temp < 50.0e0:
            th_cond = (
                0.0731268577e0
                - 0.0013826223e0 * temp
                + 3.55551245e-5 * temp**2
                - 2.32185411e-7 * temp**3
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif temp < 51.0e0:
            th_cond = 4.450475632499988e-2 + 3.871124250000024e-4 * temp

        # Linear fit
        else:
            th_cond = 0.04235676e0 + 0.00042923e0 * temp

        return th_cond

    @staticmethod
    def al_th_cond(temp: float) -> float:
        """Author : S. Kahn
        Subroutine calculating temperature dependent Al thermal conductivity

        :param temp: Helium temperature [K]
        :type temp: float

        :return th_cond: Themal conductivity [W/(m.K)]
        :type th_cond: float
        """

        # Fiting range verification
        if temp < 15.0e0 or temp > 150.0e0:
            eh.fdiags[0] = temp
            eh.report_error(258)

        # fit 15 < T < 60 K (order 3 poly)
        if temp < 60.0e0:
            th_cond = (
                16332.2073e0
                - 776.91775e0 * temp
                + 13.405688e0 * temp**2
                - 8.01e-02 * temp**3
            )

        # Linear interpolation between the fits to avoid discontinuity
        elif temp < 70.0e0:
            th_cond = 1587.9108966527328e0 - 15.19819661087886e0 * temp

        # fit 70 < T < 150 K (order 2 poly)
        elif temp < 150.0e0:
            th_cond = 1782.77406e0 - 24.7778504e0 * temp + 9.70842050e-2 * temp**2

        # constant value after that set with the fit upper limit to avoid discontinuities
        else:
            th_cond = 250.4911087866094e0

        return th_cond

    @staticmethod
    @numba.njit(cache=True)
    def tfcind(tfthk, r_tf_arc, z_tf_arc):
        """Calculates the self inductance of a TF coil
        This routine calculates the self inductance of a TF coil
        approximated by a straight inboard section and two elliptical arcs.
        The inductance of the TFC (considered as a single axisymmetric turn)
        is calculated by numerical integration over the cross-sectional area.
        The contribution from the cross-sectional area of the
        coil itself is calculated by taking the field as B(r)/2.
        The field in the dr_bore is calculated for unit current.
        Top/bottom symmetry is assumed.

        :param tfthk: TF coil thickness (m)
        :type tfthk: float
        """
        NINTERVALS = 100

        # Integrate over the whole TF area, including the coil thickness.
        x0 = r_tf_arc[1]
        y0 = z_tf_arc[1]

        # Minor and major radii of the inside and outside perimeters of the the
        # Inboard leg and arc.
        # Average the upper and lower halves, which are different in the
        # single null case
        ai = r_tf_arc[1] - r_tf_arc[0]
        bi = (z_tf_arc[1] - z_tf_arc[3]) / 2.0e0 - z_tf_arc[0]
        ao = ai + tfthk
        bo = bi + tfthk
        # Interval used for integration
        dr = ao / NINTERVALS
        # Start both integrals from the centre-point where the arcs join.
        # Initialise major radius
        r = x0 - dr / 2.0e0

        ind_tf_coil = 0

        for _ in range(NINTERVALS):
            # Field in the dr_bore for unit current
            b = RMU0 / (2.0e0 * np.pi * r)
            # Find out if there is a dr_bore
            if x0 - r < ai:
                h_bore = y0 + bi * np.sqrt(1 - ((r - x0) / ai) ** 2)
                h_thick = bo * np.sqrt(1 - ((r - x0) / ao) ** 2) - h_bore
            else:
                h_bore = 0.0e0
                # Include the contribution from the straight section
                h_thick = bo * np.sqrt(1 - ((r - x0) / ao) ** 2) + z_tf_arc[0]

            # Assume B in TF coil = 1/2  B in dr_bore
            # Multiply by 2 for upper and lower halves of coil
            ind_tf_coil += b * dr * (2.0e0 * h_bore + h_thick)
            r = r - dr

        # Outboard arc
        ai = r_tf_arc[2] - r_tf_arc[1]
        bi = (z_tf_arc[1] - z_tf_arc[3]) / 2.0e0
        ao = ai + tfthk
        bo = bi + tfthk
        dr = ao / NINTERVALS
        # Initialise major radius
        r = x0 + dr / 2.0e0

        for _ in range(NINTERVALS):
            # Field in the dr_bore for unit current
            b = RMU0 / (2.0e0 * np.pi * r)
            # Find out if there is a dr_bore
            if r - x0 < ai:
                h_bore = y0 + bi * np.sqrt(1 - ((r - x0) / ai) ** 2)
                h_thick = bo * np.sqrt(1 - ((r - x0) / ao) ** 2) - h_bore
            else:
                h_bore = 0.0e0
                h_thick = bo * np.sqrt(1 - ((r - x0) / ao) ** 2)

            # Assume B in TF coil = 1/2  B in dr_bore
            # Multiply by 2 for upper and lower halves of coil
            ind_tf_coil += b * dr * (2.0e0 * h_bore + h_thick)
            r = r + dr

        return ind_tf_coil

    def tf_coil_area_and_masses(self):
        """Subroutine to calculate the TF coil areas and masses"""
        vol_case = 0.0e0  # Total TF case volume [m3]
        vol_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_gr_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_cond = 0.0e0  # Total conductor insulator volume [m3]
        vol_ins_leg = 0.0e0  # Outboard leg turn isulation volume [m3]
        vol_gr_ins_leg = 0.0e0  # Outboard leg turn insulation volume [m3]
        vol_cond_leg = 0.0e0  # Outboard leg conductor insulator volume [m3]
        # ---

        # Surface areas (for cryo system) [m2]
        wbtf = (
            build_variables.r_tf_inboard_out
            * np.sin(sctfcoil_module.rad_tf_coil_toroidal)
            - build_variables.r_tf_inboard_in * sctfcoil_module.tan_theta_coil
        )
        tfcoil_variables.tfocrn = (
            build_variables.r_tf_inboard_in * sctfcoil_module.tan_theta_coil
        )
        tfcoil_variables.tficrn = tfcoil_variables.tfocrn + wbtf

        # Total surface area of two toroidal shells covering the TF coils [m2]
        # (inside and outside surfaces)
        # = 2 * centroid coil length * 2 pi R, where R is average of i/b and o/b centres
        tfcoil_variables.tfcryoarea = (
            2.0e0
            * tfcoil_variables.len_tf_coil
            * constants.twopi
            * 0.5e0
            * (build_variables.r_tf_inboard_mid + build_variables.r_tf_outboard_mid)
        )

        # Superconductor coil design specific calculation
        # ---
        if tfcoil_variables.i_tf_sup == 1:
            # Mass of case [kg]
            # ***

            # Mass of ground-wall insulation [kg]
            # (assumed to be same density/material as turn insulation)
            tfcoil_variables.whtgw = (
                tfcoil_variables.len_tf_coil
                * (sctfcoil_module.awpc - sctfcoil_module.awptf)
                * tfcoil_variables.dcondins
            )

            # The length of the vertical section is that of the first (inboard) segment
            # = height of TF coil inner edge + (2 * coil thickness)
            tfcoil_variables.cplen = (2.0e0 * build_variables.z_tf_inside_half) + (
                2.0e0 * build_variables.dr_tf_inboard
            )

            # The 2.2 factor is used as a scaling factor to fit
            # to the ITER-FDR value of 450 tonnes; see CCFE note T&M/PKNIGHT/PROCESS/026
            if physics_variables.itart == 1:
                # tfcoil_variables.len_tf_coil does not include inboard leg ('centrepost') length in TART
                tfcoil_variables.whtcas = (
                    2.2e0
                    * tfcoil_variables.dcase
                    * (
                        tfcoil_variables.cplen * tfcoil_variables.acasetf
                        + tfcoil_variables.len_tf_coil * tfcoil_variables.acasetfo
                    )
                )
            else:
                tfcoil_variables.whtcas = (
                    2.2e0
                    * tfcoil_variables.dcase
                    * (
                        tfcoil_variables.cplen * tfcoil_variables.acasetf
                        + (tfcoil_variables.len_tf_coil - tfcoil_variables.cplen)
                        * tfcoil_variables.acasetfo
                    )
                )

            # ***

            # Masses of conductor constituents
            # ---------------------------------
            # Superconductor mass [kg]
            # Includes space allowance for central helium channel, area tfcoil_variables.a_tf_wp_coolant_channels
            tfcoil_variables.whtconsc = (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.a_tf_turn_cable_space
                * (1.0e0 - tfcoil_variables.vftf)
                * (1.0e0 - tfcoil_variables.fcutfsu)
                - tfcoil_variables.len_tf_coil
                * tfcoil_variables.a_tf_wp_coolant_channels
            ) * tfcoil_variables.dcond[tfcoil_variables.i_tf_sc_mat - 1]

            # Copper mass [kg]
            tfcoil_variables.whtconcu = (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.a_tf_turn_cable_space
                * (1.0e0 - tfcoil_variables.vftf)
                * tfcoil_variables.fcutfsu
                - tfcoil_variables.len_tf_coil
                * tfcoil_variables.a_tf_wp_coolant_channels
            ) * constants.dcopper
            if tfcoil_variables.whtconcu <= 0.0e0:
                tfcoil_variables.whtconcu = 0.0e0

            # Steel conduit (sheath) mass [kg]
            tfcoil_variables.m_tf_turn_steel_conduit = (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.a_tf_turn_steel
                * fwbs_variables.denstl
            )

            # Conduit insulation mass [kg]
            # (tfcoil_variables.a_tf_coil_wp_turn_insulation already contains tfcoil_variables.n_tf_turn)
            tfcoil_variables.whtconin = (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.a_tf_coil_wp_turn_insulation
                * tfcoil_variables.dcondins
            )

            # Total conductor mass [kg]
            tfcoil_variables.whtcon = (
                tfcoil_variables.whtconsc
                + tfcoil_variables.whtconcu
                + tfcoil_variables.m_tf_turn_steel_conduit
                + tfcoil_variables.whtconin
            )
            # ---------------------------------

            # Total TF coil mass [kg] (all coils)
            tfcoil_variables.m_tf_coils_total = (
                tfcoil_variables.whtcas
                + tfcoil_variables.whtcon
                + tfcoil_variables.whtgw
            ) * tfcoil_variables.n_tf_coils

            # If spherical tokamak, distribute between centrepost and outboard legs
            # (in this case, total TF coil length = inboard `cplen` + outboard `len_tf_coil`)
            if physics_variables.itart == 1:
                tfleng_sph = tfcoil_variables.cplen + tfcoil_variables.len_tf_coil
                tfcoil_variables.whtcp = tfcoil_variables.m_tf_coils_total * (
                    tfcoil_variables.cplen / tfleng_sph
                )
                tfcoil_variables.whttflgs = tfcoil_variables.m_tf_coils_total * (
                    tfcoil_variables.len_tf_coil / tfleng_sph
                )

        # Resitivive magnets weights
        # ---
        # Rem SK : No casing for the outboard leg is considered for now #
        else:
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
                vol_cond_leg = tfcoil_variables.len_tf_coil * sctfcoil_module.a_leg_cond

                # Total TF conductor volume [m3]
                vol_cond = (
                    tfcoil_variables.vol_cond_cp
                    + tfcoil_variables.n_tf_coils * vol_cond_leg
                )

                # Outboard leg TF turn insulation layer volume (per leg) [m3]
                vol_ins_leg = tfcoil_variables.len_tf_coil * sctfcoil_module.a_leg_ins

                # Total turn insulation layer volume [m3]
                vol_ins = (
                    sctfcoil_module.vol_ins_cp
                    + tfcoil_variables.n_tf_coils * vol_ins_leg
                )

                # Ouboard leg TF ground insulation layer volume (per leg) [m3]
                vol_gr_ins_leg = (
                    tfcoil_variables.len_tf_coil * sctfcoil_module.a_leg_gr_ins
                )

                # Total ground insulation layer volume [m3]
                vol_gr_ins = (
                    sctfcoil_module.vol_gr_ins_cp
                    + tfcoil_variables.n_tf_coils * vol_gr_ins_leg
                )

                # Total volume of the CP casing [m3]
                # Rem : no outer leg case
                vol_case = sctfcoil_module.vol_case_cp

            # No joints
            # ---
            else:
                # Total TF outer leg conductor volume [m3]
                vol_cond = (
                    tfcoil_variables.len_tf_coil
                    * sctfcoil_module.a_leg_cond
                    * tfcoil_variables.n_tf_coils
                )

                # Total turn insulation layer volume [m3]
                vol_ins = (
                    tfcoil_variables.len_tf_coil
                    * sctfcoil_module.a_leg_ins
                    * tfcoil_variables.n_tf_coils
                )

                # Total ground insulation volume [m3]
                vol_gr_ins = (
                    tfcoil_variables.len_tf_coil
                    * sctfcoil_module.a_leg_gr_ins
                    * tfcoil_variables.n_tf_coils
                )

                # Total case volume [m3]
                vol_case = (
                    tfcoil_variables.len_tf_coil
                    * tfcoil_variables.acasetf
                    * tfcoil_variables.n_tf_coils
                )

            # ---
            # -------

            # Copper magnets casing/conductor weights per coil [kg]
            if tfcoil_variables.i_tf_sup == 0:
                tfcoil_variables.whtcas = (
                    fwbs_variables.denstl * vol_case / tfcoil_variables.n_tf_coils
                )  # Per TF leg, no casing for outer leg
                tfcoil_variables.whtconcu = (
                    constants.dcopper * vol_cond / tfcoil_variables.n_tf_coils
                )
                tfcoil_variables.whtconal = 0.0e0

                # Outer legs/CP weights
                if physics_variables.itart == 1:
                    # Weight of all the TF legs
                    tfcoil_variables.whttflgs = tfcoil_variables.n_tf_coils * (
                        constants.dcopper * vol_cond_leg
                        + tfcoil_variables.dcondins * (vol_ins_leg + vol_gr_ins_leg)
                    )

                    # CP weight
                    tfcoil_variables.whtcp = (
                        constants.dcopper * tfcoil_variables.vol_cond_cp
                        + tfcoil_variables.dcondins
                        * (sctfcoil_module.vol_ins_cp + sctfcoil_module.vol_gr_ins_cp)
                        + sctfcoil_module.vol_case_cp * fwbs_variables.denstl
                    )

            # Cryo-aluminium conductor weights
            # Casing made of re-inforced aluminium alloy
            elif tfcoil_variables.i_tf_sup == 2:
                # Casing weight (CP only if physics_variables.itart = 1)bper leg/coil
                tfcoil_variables.whtcas = (
                    constants.dalu * vol_case / tfcoil_variables.n_tf_coils
                )
                tfcoil_variables.whtconcu = 0.0e0
                tfcoil_variables.whtconal = (
                    constants.dalu * vol_cond / tfcoil_variables.n_tf_coils
                )

                # Outer legs/CP weights
                if physics_variables.itart == 1:
                    # Weight of all the TF legs
                    tfcoil_variables.whttflgs = tfcoil_variables.n_tf_coils * (
                        constants.dalu * vol_cond_leg
                        + tfcoil_variables.dcondins * (vol_ins_leg + vol_gr_ins_leg)
                    )

                    # CP weight
                    tfcoil_variables.whtcp = (
                        constants.dalu * tfcoil_variables.vol_cond_cp
                        + tfcoil_variables.dcondins
                        * (sctfcoil_module.vol_ins_cp + sctfcoil_module.vol_gr_ins_cp)
                        + sctfcoil_module.vol_case_cp * fwbs_variables.denstl
                    )

            # Turn insulation mass [kg]
            tfcoil_variables.whtconin = (
                tfcoil_variables.dcondins * vol_ins / tfcoil_variables.n_tf_coils
            )

            # Ground wall insulation layer weight
            tfcoil_variables.whtgw = (
                tfcoil_variables.dcondins * vol_gr_ins / tfcoil_variables.n_tf_coils
            )

            # Total weight
            tfcoil_variables.m_tf_coils_total = (
                tfcoil_variables.whtcas
                + tfcoil_variables.whtconcu
                + tfcoil_variables.whtconal
                + tfcoil_variables.whtconin
                + tfcoil_variables.whtgw
            ) * tfcoil_variables.n_tf_coils

    @staticmethod
    @numba.njit(cache=True)
    def stresscl(
        n_tf_layer,
        n_radial_array,
        n_tf_wp_layers,
        i_tf_bucking,
        r_tf_inboard_in,
        dr_bore,
        z_tf_inside_half,
        f_z_cs_tf_internal,
        dr_cs,
        i_tf_inside_cs,
        dr_tf_inboard,
        dr_cs_tf_gap,
        i_pf_conductor,
        j_cs_flat_top_end,
        j_cs_pulse_start,
        c_pf_coil_turn_peak_input,
        n_pf_coils_in_group,
        ld_ratio_cst,
        r_out_cst,
        f_a_cs_steel,
        eyoung_steel,
        poisson_steel,
        eyoung_cond_axial,
        poisson_cond_axial,
        eyoung_cond_trans,
        poisson_cond_trans,
        eyoung_ins,
        poisson_ins,
        dx_tf_turn_insulation,
        eyoung_copper,
        poisson_copper,
        i_tf_sup,
        eyoung_res_tf_buck,
        r_wp_inner,
        tan_theta_coil,
        rad_tf_coil_toroidal,
        r_wp_outer,
        a_tf_steel,
        a_case_front,
        a_case_nose,
        tfinsgap,
        tinstf,
        n_tf_turn,
        i_tf_turns_integer,
        t_cable,
        dr_tf_turn_cable_space,
        dia_tf_turn_coolant_channel,
        fcutfsu,
        dx_tf_turn_steel,
        t_lat_case_av,
        t_wp_toroidal_av,
        a_tf_ins,
        aswp,
        acond,
        awpc,
        eyoung_al,
        poisson_al,
        fcoolcp,
        n_tf_graded_layers,
        c_tf_total,
        dr_tf_plasma_case,
        i_tf_stress_model,
        vforce_inboard_tot,
        i_tf_tresca,
        acasetf,
        vforce,
        a_tf_turn_steel,
    ):
        """TF coil stress routine


        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        author: S Kahn, CCFE, Culham Science Centre
        author: J Galambos, FEDC/ORNL

        This subroutine sets up the stress calculations for the
        TF coil set.
        PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
        """
        jeff = np.zeros((n_tf_layer,))
        # Effective current density [A/m2]

        radtf = np.zeros((n_tf_layer + 1,))
        # Radii used to define the layers used in the stress models [m]
        # Layers are labelled from inboard to outbard

        eyoung_trans = np.zeros((n_tf_layer,))
        # Young's moduli (one per layer) of the TF coil in the
        # transverse (radial/toroidal) direction. Used in the stress
        # models [Pa]

        poisson_trans = np.zeros(
            (n_tf_layer,),
        )
        # Poisson's ratios (one per layer) of the TF coil between the
        # two transverse directions (radial and toroidal). Used in the
        # stress models.

        eyoung_member_array = np.zeros((n_tf_wp_layers,))
        # Array to store the Young's moduli of the members to composite into smeared
        # properties [Pa]

        poisson_member_array = np.zeros((n_tf_wp_layers,))
        # Array to store the Poisson's ratios of the members to composite into smeared
        # properti

        l_member_array = np.zeros((n_tf_wp_layers,))
        # Array to store the linear dimension (thickness) of the members to composite into smeared
        # properties [m]

        eyoung_axial = np.zeros((n_tf_layer,))
        # Young's moduli (one per layer) of the TF coil in the vertical
        # direction used in the stress models [Pa]

        poisson_axial = np.zeros((n_tf_layer,))
        # Poisson's ratios (one per layer) of the TF coil between the
        # vertical and transverse directions (in that order). Used in the
        # stress models. d(transverse strain)/d(vertical strain) with
        # only vertical stress.

        sig_tf_wp_av_z = np.zeros(((n_tf_layer - i_tf_bucking) * n_radial_array,))
        # TF Inboard leg WP smeared vertical stress r distribution at mid-plane [Pa]

        sig_tf_r_max = np.zeros((n_tf_layer,))
        # Radial stress of the point of maximum shear stress of each layer [Pa]

        sig_tf_t_max = np.zeros((n_tf_layer,))
        # Toroidal stress of the point of maximum shear stress of each layer [Pa]

        sig_tf_z_max = np.zeros((n_tf_layer,))
        # Vertical stress of the point of maximum shear stress of each layer [Pa]
        # Rem : Currently constant but will be r dependent in the future

        sig_tf_vmises_max = np.zeros((n_tf_layer,))
        # Von-Mises stress of the point of maximum shear stress of each layer [Pa]

        s_shear_tf_peak = np.zeros((n_tf_layer,))
        # Maximum shear stress, for the Tresca yield criterion of each layer [Pa]
        # If the CEA correction is addopted, the CEA corrected value is used

        sig_tf_z = np.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg vertical tensile stress at mid-plane [Pa]

        sig_tf_smeared_r = np.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg radial smeared stress r distribution at mid-plane [Pa]

        sig_tf_smeared_t = np.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg tangential smeared stress r distribution at mid-plane [Pa]

        sig_tf_smeared_z = np.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg vertical smeared stress r distribution at mid-plane [Pa]

        fac_sig_z_wp_av = 0.0
        # WP averaged vertical stress unsmearing factor

        str_tf_r = np.zeros((n_radial_array * n_tf_layer,))
        str_tf_t = np.zeros((n_radial_array * n_tf_layer,))
        str_tf_z = np.zeros((n_radial_array * n_tf_layer,))

        sig_tf_case = None
        sig_tf_cs_bucked = None
        str_wp = None
        casestr = None
        insstrain = None

        # New extended plane strain model can handle it
        if (
            abs(r_tf_inboard_in) < np.finfo(float(r_tf_inboard_in)).eps
            and i_tf_stress_model != 2
        ):
            raise ProcessValueError("r_tf_inboard_in is ~= 0", 245)

        # TODO: following is no longer used/needed?
        # if tfcoil_variables.a_tf_turn_cable_space >= 0.0e0:
        #     tcbs = numpy.sqrt(tfcoil_variables.a_tf_turn_cable_space)
        # else:
        #     tcbs = 0.0e0

        # LAYER ELASTIC PROPERTIES
        # ------------------------
        # Number of bucking layers
        n_tf_bucking = i_tf_bucking

        # CS properties (bucked and wedged)
        # ---
        if i_tf_bucking >= 2:
            # Calculation performed at CS flux swing (no current on CS)
            jeff[0] = 0.0e0

            # Inner radius of the CS
            if i_tf_inside_cs == 1:
                # CS not used as wedge support i_tf_inside_cs = 1
                radtf[0] = 0.001
            else:
                radtf[0] = dr_bore

            # Superconducting CS
            if i_pf_conductor == 0:
                # Getting the turn dimention from scratch
                # as the TF is called before CS in caller.f90
                # -#

                # CS vertical cross-section area [m2]
                if i_tf_inside_cs == 1:
                    a_oh = (
                        2.0e0
                        * z_tf_inside_half
                        * f_z_cs_tf_internal
                        * (dr_bore - dr_tf_inboard)
                    )
                else:
                    a_oh = 2.0e0 * z_tf_inside_half * f_z_cs_tf_internal * dr_cs

                # Maximum current in Central Solenoid, at either BOP or EOF [MA-turns]
                # Absolute value
                curr_oh_max = (
                    1.0e-6 * np.maximum(j_cs_flat_top_end, j_cs_pulse_start) * a_oh
                )

                #  Number of turns
                n_oh_turns = (
                    1.0e6
                    * curr_oh_max
                    / c_pf_coil_turn_peak_input[sum(n_pf_coils_in_group)]
                )

                # CS Turn vertical cross-sectionnal area
                a_cs_turn = a_oh / n_oh_turns

                # CS coil turn geometry calculation - stadium shape
                # Literature: https://doi.org/10.1016/j.fusengdes.2017.04.052
                dz_cs_turn = (
                    a_cs_turn / ld_ratio_cst
                ) ** 0.5  # width of cs turn conduit
                dr_cs_turn = ld_ratio_cst * dz_cs_turn  # length of cs turn conduit
                # Radius of turn space = radius_cs_turn_cable_space
                # Radius of curved outer corrner r_out_cst = 3mm from literature
                # ld_ratio_cst = 70 / 22 from literature
                p1 = ((dr_cs_turn - dz_cs_turn) / np.pi) ** 2
                p2 = (
                    (dr_cs_turn * dz_cs_turn)
                    - (4 - np.pi) * (r_out_cst**2)
                    - (a_cs_turn * f_a_cs_steel)
                ) / np.pi
                radius_cs_turn_cable_space = -(
                    (dr_cs_turn - dz_cs_turn) / np.pi
                ) + np.sqrt(p1 + p2)
                t_cond_oh = (
                    dz_cs_turn / 2
                ) - radius_cs_turn_cable_space  # thickness of steel conduit in cs turn

                # OH/CS conduit thickness calculated assuming square conduit [m]
                # The CS insulation layer is assumed to the same as the TF one

                # CS turn cable space thickness
                t_cable_oh = radius_cs_turn_cable_space * 2
                # -#

                # Smeared elastic properties of the CS
                # These smearing functions were written assuming transverse-
                # isotropic materials; that is not true of the CS, where the
                # stiffest dimension is toroidal and the radial and vertical
                # dimension are less stiff. Nevertheless this attempts to
                # hit the mark.
                # [EDIT: eyoung_cond is for the TF coil, not the CS coil]

                # Get transverse properties
                (eyoung_trans[0], a_working, poisson_trans[0]) = eyoung_parallel(
                    eyoung_steel,
                    f_a_cs_steel,
                    poisson_steel,
                    eyoung_cond_axial,
                    1e0 - f_a_cs_steel,
                    poisson_cond_axial,
                )

                # Get vertical properties
                # Split up into "members", concentric squares in cross section
                # (described in Figure 10 of the TF coil documentation)
                # Conductor
                eyoung_member_array[0] = eyoung_cond_trans
                poisson_member_array[0] = poisson_cond_trans
                l_member_array[0] = t_cable_oh
                # Steel conduit
                eyoung_member_array[1] = eyoung_steel
                poisson_member_array[1] = poisson_steel
                l_member_array[1] = 2 * t_cond_oh
                # Insulation
                eyoung_member_array[2] = eyoung_ins
                poisson_member_array[2] = poisson_ins
                l_member_array[2] = 2 * dx_tf_turn_insulation
                # [EDIT: Add central cooling channel? Would be new member #1]

                # Compute the composited (smeared) properties
                (
                    eyoung_axial[0],
                    a_working,
                    poisson_axial[0],
                    eyoung_cs_stiffest_leg,
                ) = eyoung_t_nested_squares(
                    3, eyoung_member_array, l_member_array, poisson_member_array
                )

            # resistive CS (copper)
            else:
                # Here is a rough approximation
                eyoung_trans[0] = eyoung_copper
                eyoung_axial[0] = eyoung_copper
                poisson_trans[0] = poisson_copper
                poisson_axial[0] = poisson_copper

        # ---

        # CS-TF interlayer properties
        # ---
        if i_tf_bucking == 3:
            # No current in this layer
            jeff[1] = 0.0e0

            # Outer radius of the CS
            if i_tf_inside_cs == 1:
                radtf[1] = dr_bore - dr_tf_inboard - dr_cs_tf_gap
            else:
                radtf[1] = dr_bore + dr_cs

            # Assumed to be Kapton for the moment
            # Ref : https://www.dupont.com/content/dam/dupont/products-and-services/membranes-and-films/polyimde-films/documents/DEC-Kapton-summary-of-properties.pdf
            eyoung_trans[1] = 2.5e9
            eyoung_axial[1] = 2.5e9
            poisson_trans[1] = 0.34e0  # Default value for young modulus
            poisson_axial[1] = 0.34e0  # Default value for young modulus

        # ---

        # bucking cylinder/casing properties
        # ---
        if i_tf_bucking >= 1:
            # No current in bucking cylinder/casing
            jeff[n_tf_bucking - 1] = 0.0e0

            if i_tf_sup == 1:
                eyoung_trans[n_tf_bucking - 1] = eyoung_steel
                eyoung_axial[n_tf_bucking - 1] = eyoung_steel
                poisson_trans[n_tf_bucking - 1] = poisson_steel
                poisson_axial[n_tf_bucking - 1] = poisson_steel

            # Bucking cylinder properties
            else:
                eyoung_trans[n_tf_bucking - 1] = eyoung_res_tf_buck
                eyoung_axial[n_tf_bucking - 1] = eyoung_res_tf_buck
                poisson_trans[n_tf_bucking - 1] = poisson_steel  # Seek better value #
                poisson_axial[n_tf_bucking - 1] = poisson_steel  # Seek better value #

            # Innernost TF casing radius
            radtf[n_tf_bucking - 1] = r_tf_inboard_in

        # ---

        # (Super)conductor layer properties
        # ---
        # SC coil
        if i_tf_sup == 1:
            # Inner/outer radii of the layer representing the WP in stress calculations [m]
            # These radii are chosen to preserve the true WP area; see Issue #1048
            r_wp_inner_eff = r_wp_inner * np.sqrt(tan_theta_coil / rad_tf_coil_toroidal)
            r_wp_outer_eff = r_wp_outer * np.sqrt(tan_theta_coil / rad_tf_coil_toroidal)

            # Area of the cylinder representing the WP in stress calculations [m2]
            a_wp_eff = (r_wp_outer_eff**2 - r_wp_inner_eff**2) * rad_tf_coil_toroidal

            # Steel cross-section under the area representing the WP in stress calculations [m2]
            a_wp_steel_eff = a_tf_steel - a_case_front - a_case_nose

            # WP effective insulation thickness (SC only) [m]
            # include groundwall insulation + insertion gap in tfcoil_variables.dx_tf_turn_insulation
            # inertion gap is tfcoil_variables.tfinsgap on 4 sides
            t_ins_eff = dx_tf_turn_insulation + (tfinsgap + tinstf) / n_tf_turn

            # Effective WP young modulus in the toroidal direction [Pa]
            # The toroidal property drives the stress calculation (J. Last report no 4)
            # Hence, the radial direction is relevant for the property smearing
            # Rem : This assumption might be re-defined for bucked and wedged design

            # Non-integer or interger number of turns
            t_cable_eyng = (
                t_cable if i_tf_turns_integer == 0 else dr_tf_turn_cable_space
            )

            # Average WP Young's modulus in the transverse
            # (radial and toroidal) direction
            # Split up into "members", concentric squares in cross section
            # (described in Figure 10 of the TF coil documentation)
            # Helium
            eyoung_member_array[0] = 0e0
            poisson_member_array[0] = poisson_steel
            l_member_array[0] = dia_tf_turn_coolant_channel
            # Conductor and co-wound copper
            (
                eyoung_member_array[1],
                l_member_array[1],
                poisson_member_array[1],
            ) = eyoung_series(
                np.double(eyoung_cond_trans),
                (t_cable_eyng - dia_tf_turn_coolant_channel) * (1.0e0 - fcutfsu),
                np.double(poisson_cond_trans),
                np.double(eyoung_copper),
                (t_cable_eyng - dia_tf_turn_coolant_channel) * fcutfsu,
                np.double(poisson_copper),
            )
            # Steel conduit
            eyoung_member_array[2] = eyoung_steel
            poisson_member_array[2] = poisson_steel
            l_member_array[2] = 2 * dx_tf_turn_steel
            # Insulation
            eyoung_member_array[3] = eyoung_ins
            poisson_member_array[3] = poisson_ins
            l_member_array[3] = 2 * t_ins_eff

            # Compute the composited (smeared) properties
            (
                eyoung_wp_trans,
                a_working,
                poisson_wp_trans,
                eyoung_wp_stiffest_leg,
            ) = eyoung_t_nested_squares(
                4,
                eyoung_member_array,
                l_member_array,
                poisson_member_array,
            )

            # Lateral casing correction (series-composition)
            (eyoung_wp_trans_eff, a_working, poisson_wp_trans_eff) = eyoung_series(
                eyoung_wp_trans,
                np.double(t_wp_toroidal_av),
                poisson_wp_trans,
                np.double(eyoung_steel),
                2.0e0 * t_lat_case_av,
                np.double(poisson_steel),
            )

            poisson_wp_trans_eff = np.double(poisson_wp_trans_eff)

            # Average WP Young's modulus in the vertical direction
            # Split up into "members", concentric squares in cross section
            # (described in Figure 10 of the TF coil documentation)
            # Steel conduit
            eyoung_member_array[0] = eyoung_steel
            poisson_member_array[0] = poisson_steel
            l_member_array[0] = aswp
            # Insulation
            eyoung_member_array[1] = eyoung_ins
            poisson_member_array[1] = poisson_ins
            l_member_array[1] = a_tf_ins
            # Copper
            eyoung_member_array[2] = eyoung_copper
            poisson_member_array[2] = poisson_copper
            l_member_array[2] = acond * fcutfsu
            # Conductor
            eyoung_member_array[3] = eyoung_cond_axial
            poisson_member_array[3] = poisson_cond_axial
            l_member_array[3] = acond * (1.0e0 - fcutfsu)
            # Helium and void
            eyoung_member_array[4] = 0e0
            poisson_member_array[4] = poisson_steel
            l_member_array[4] = awpc - acond - a_tf_ins - aswp
            # Compute the composite / smeared properties:
            (eyoung_wp_axial, a_working, poisson_wp_axial) = eyoung_parallel_array(
                5,
                eyoung_member_array,
                l_member_array,
                poisson_member_array,
            )

            # Average WP Young's modulus in the vertical direction, now including the lateral case
            # Parallel-composite the steel and insulation, now including the lateral case (sidewalls)
            (eyoung_wp_axial_eff, a_working, poisson_wp_axial_eff) = eyoung_parallel(
                eyoung_steel,
                a_wp_steel_eff - aswp,
                poisson_steel,
                eyoung_wp_axial,
                a_working,
                poisson_wp_axial,
            )

        # Resistive coil
        else:
            # Picking the conductor material Young's modulus
            if i_tf_sup == 0:
                eyoung_cond = eyoung_copper
                poisson_cond = poisson_copper
            elif i_tf_sup == 2:
                eyoung_cond = eyoung_al
                poisson_cond = poisson_al

            # Effective WP young modulus in the toroidal direction [Pa]
            # Rem : effect of cooling pipes and insulation not taken into account
            #       for now as it needs a radially dependent Young modulus
            eyoung_wp_trans_eff = eyoung_cond
            eyoung_wp_trans = np.double(eyoung_cond)
            poisson_wp_trans_eff = np.double(poisson_cond)
            poisson_wp_trans = np.double(poisson_cond)

            # WP area using the stress model circular geometry (per coil) [m2]
            a_wp_eff = (r_wp_outer**2 - r_wp_inner**2) * rad_tf_coil_toroidal

            # Effective conductor region young modulus in the vertical direction [Pa]
            # Parallel-composite conductor and insulator
            (eyoung_wp_axial, a_working, poisson_wp_axial) = eyoung_parallel(
                eyoung_cond,
                (a_wp_eff - a_tf_ins) * (1.0e0 - fcoolcp),
                poisson_cond,
                eyoung_ins,
                a_tf_ins,
                poisson_ins,
            )
            # Parallel-composite cooling pipes into that
            (eyoung_wp_axial, a_working, poisson_wp_axial) = eyoung_parallel(
                0e0,
                (a_wp_eff - a_tf_ins) * fcoolcp,
                poisson_cond,
                eyoung_wp_axial,
                a_working,
                poisson_wp_axial,
            )

            # Effective young modulus used in stress calculations
            eyoung_wp_axial_eff = eyoung_wp_axial
            poisson_wp_axial_eff = poisson_wp_axial

            # Effect conductor layer inner/outer radius
            r_wp_inner_eff = np.double(r_wp_inner)
            r_wp_outer_eff = np.double(r_wp_outer)

        # Thickness of the layer representing the WP in stress calcualtions [m]
        dr_tf_wp_eff = r_wp_outer_eff - r_wp_outer_eff

        # Thickness of WP with homogeneous stress property [m]
        dr_wp_layer = dr_tf_wp_eff / n_tf_graded_layers

        for ii in range(np.intc(n_tf_graded_layers)):
            # Homogeneous current in (super)conductor
            jeff[n_tf_bucking + ii] = c_tf_total / (
                np.pi * (r_wp_outer_eff**2 - r_wp_inner_eff**2)
            )

            # Same thickness for all WP layers in stress calculation
            radtf[n_tf_bucking + ii] = r_wp_inner_eff + ii * dr_wp_layer

            # Young modulus
            eyoung_trans[n_tf_bucking + ii] = eyoung_wp_trans_eff
            eyoung_axial[n_tf_bucking + ii] = eyoung_wp_axial_eff

            # Poisson's ratio
            poisson_trans[n_tf_bucking + ii] = poisson_wp_trans_eff
            poisson_axial[n_tf_bucking + ii] = poisson_wp_axial_eff

        # Steel case on the plasma side of the inboard TF coil
        # As per Issue #1509
        jeff[n_tf_layer - 1] = 0.0e0
        radtf[n_tf_layer - 1] = r_wp_outer_eff
        eyoung_trans[n_tf_layer - 1] = eyoung_steel
        eyoung_axial[n_tf_layer - 1] = eyoung_steel
        poisson_trans[n_tf_layer - 1] = poisson_steel
        poisson_axial[n_tf_layer - 1] = poisson_steel

        # last layer radius
        radtf[n_tf_layer] = r_wp_outer_eff + dr_tf_plasma_case

        # The ratio between the true cross sectional area of the
        # front case, and that considered by the plane strain solver
        f_tf_stress_front_case = (
            a_case_front
            / rad_tf_coil_toroidal
            / (radtf[n_tf_layer] ** 2 - radtf[n_tf_layer - 1] ** 2)
        )

        # Correct for the missing axial stiffness from the missing
        # outer case steel as per the updated description of
        # Issue #1509
        eyoung_axial[n_tf_layer - 1] = (
            eyoung_axial[n_tf_layer - 1] * f_tf_stress_front_case
        )

        # ---
        # ------------------------

        # RADIAL STRESS SUBROUTINES CALL
        # ------------------------------
        # Stress model not valid the TF does not contain any hole
        # (Except if i_tf_plane_stress == 2; see Issue 1414)
        # Current action : trigger and error and add a little hole
        #                  to allow stress calculations
        # Rem SK : Can be easily ameneded playing around the boundary conditions
        # New extended plane strain model can handle it
        if abs(radtf[0]) < np.finfo(float(radtf[0])).eps and i_tf_stress_model != 2:
            radtf[0] = 1.0e-9
        # ---

        # Old generalized plane stress model
        # ---
        if i_tf_stress_model == 1:
            # Plane stress calculation (SC) [Pa]

            (sig_tf_r, sig_tf_t, deflect, radial_array) = plane_stress(
                nu=poisson_trans,
                rad=radtf,
                ey=eyoung_trans,
                j=jeff,
                nlayers=int(n_tf_layer),
                n_radial_array=int(n_radial_array),
            )

            # Vertical stress [Pa]
            sig_tf_z[:] = vforce / (
                acasetf + a_tf_turn_steel * n_tf_turn
            )  # Array equation [EDIT: Are you sure? It doesn't look like one to me]

            # Strain in vertical direction on WP
            str_wp = sig_tf_z[n_tf_bucking] / eyoung_wp_axial_eff

            # Case strain
            casestr = sig_tf_z[n_tf_bucking - 1] / eyoung_steel

            # Radial strain in insulator
            insstrain = (
                sig_tf_r[n_radial_array - 1]
                * eyoung_wp_stiffest_leg
                / eyoung_wp_trans_eff
                / eyoung_ins
            )
        # ---

        # New generalized plane strain formulation
        # ---
        # See #1670 for why i_tf_stress_model=0 was removed
        elif i_tf_stress_model in (0, 2):
            # Extended plane strain calculation [Pa]
            # Issues #1414 and #998
            # Permits build_variables.dr_bore >= 0, O(n) in layers
            # If build_variables.dr_bore > 0, same result as generalized plane strain calculation

            (
                radial_array,
                sig_tf_r,
                sig_tf_t,
                sig_tf_z,
                str_tf_r,
                str_tf_t,
                str_tf_z,
                deflect,
            ) = extended_plane_strain(
                poisson_trans,
                poisson_axial,
                eyoung_trans,
                eyoung_axial,
                radtf,
                jeff,
                vforce_inboard_tot,
                int(n_tf_layer),
                int(n_radial_array),
                int(n_tf_bucking),
            )

            # Strain in TF conductor material
            str_wp = str_tf_z[n_tf_bucking * n_radial_array]

        # ---

        # Storing the smeared properties for output

        sig_tf_smeared_r[:] = sig_tf_r  # Array equation
        sig_tf_smeared_t[:] = sig_tf_t  # Array equation
        sig_tf_smeared_z[:] = sig_tf_z  # Array equation

        # ------------------------------

        # STRESS DISTRIBUTIONS CORRECTIONS
        # --------------------------------
        # SC central solenoid coil stress unsmearing (bucked and wedged only)
        # ---
        if i_tf_bucking >= 2 and i_pf_conductor == 0:
            # Central Solenoid (OH) steel conduit stress unsmearing factors
            for ii in range(n_radial_array):
                sig_tf_r[ii] = sig_tf_r[ii] * eyoung_cs_stiffest_leg / eyoung_axial[0]
                sig_tf_t[ii] = sig_tf_t[ii] * eyoung_steel / eyoung_trans[0]
                sig_tf_z[ii] = sig_tf_z[ii] * eyoung_cs_stiffest_leg / eyoung_axial[0]

        # ---

        # No TF vertical forces on CS and CS-TF layer (bucked and wedged only)
        # ---
        # This correction is only applied if the plane stress model is used
        # as the generalized plane strain calculates the vertical stress properly
        if i_tf_bucking >= 2 and i_tf_stress_model == 1:
            for ii in range((n_tf_bucking - 1) * n_radial_array):
                sig_tf_z[ii] = 0.0e0

        # ---

        # Toroidal coil unsmearing
        # ---
        # Copper magnets
        if i_tf_sup == 0:
            # Vertical force unsmearing factor
            fac_sig_z = eyoung_copper / eyoung_wp_axial_eff

            # Toroidal WP steel stress unsmearing factor
            fac_sig_t = 1.0e0
            fac_sig_r = 1.0e0

        elif i_tf_sup == 1:
            # Vertical WP steel stress unsmearing factor
            if i_tf_stress_model != 1:
                fac_sig_z = eyoung_steel / eyoung_wp_axial_eff
                fac_sig_z_wp_av = eyoung_wp_axial / eyoung_wp_axial_eff
            else:
                fac_sig_z = 1.0e0

            # Toroidal WP steel conduit stress unsmearing factor
            fac_sig_t = eyoung_wp_stiffest_leg / eyoung_wp_trans_eff

            # Radial WP steel conduit stress unsmearing factor
            fac_sig_r = eyoung_wp_stiffest_leg / eyoung_wp_trans_eff

        elif i_tf_sup == 2:
            # Vertical WP steel stress unsmearing factor
            fac_sig_z = eyoung_al / eyoung_wp_axial_eff

            # Toroidal WP steel stress unsmearing factor
            # NO CALCULTED FOR THE MOMENT (to be done later)
            fac_sig_t = 1.0e0
            fac_sig_r = 1.0e0

        # Application of the unsmearing to the WP layers
        # For each point within the winding pack / conductor, unsmear the
        # stress. This is n_radial_array test points within tfcoil_variables.n_tf_graded_layers
        # layers starting at n_tf_bucking + 1
        # GRADED MODIF : add another do loop to allow the graded properties
        #                to be taken into account
        for ii in range(
            n_tf_bucking * n_radial_array,
            ((n_tf_bucking + n_tf_graded_layers) * n_radial_array),
        ):
            sig_tf_wp_av_z[ii - n_tf_bucking * n_radial_array] = (
                sig_tf_z[ii] * fac_sig_z_wp_av
            )
            sig_tf_r[ii] = sig_tf_r[ii] * fac_sig_r
            sig_tf_t[ii] = sig_tf_t[ii] * fac_sig_t
            sig_tf_z[ii] = sig_tf_z[ii] * fac_sig_z

        # For each point within the front case,
        # remove the correction for the missing axial
        # stiffness as per the updated description of
        # Issue #1509
        for ii in range(
            (n_tf_bucking + n_tf_graded_layers) * n_radial_array,
            (n_tf_layer * n_radial_array),
        ):
            sig_tf_z[ii] = sig_tf_z[ii] / f_tf_stress_front_case
        # ---

        # Tresca / Von Mises yield criteria calculations
        # -----------------------------
        # Array equation
        s_shear_tf = np.maximum(
            np.absolute(sig_tf_r - sig_tf_z), np.absolute(sig_tf_z - sig_tf_t)
        )

        # Array equation

        sig_tf_vmises = np.sqrt(
            0.5e0
            * (
                (sig_tf_r - sig_tf_t) ** 2
                + (sig_tf_r - sig_tf_z) ** 2
                + (sig_tf_z - sig_tf_t) ** 2
            )
        )

        # Array equation
        s_shear_cea_tf_cond = s_shear_tf.copy()

        # SC conducting layer stress distribution corrections
        # ---
        if i_tf_sup == 1:
            # GRADED MODIF : add another do loop to allow the graded properties
            #                to be taken into account
            for ii in range(
                (n_tf_bucking * n_radial_array), n_tf_layer * n_radial_array
            ):
                # Addaped Von-mises stress calculation to WP strucure [Pa]

                svmxz = sigvm(0.0e0, sig_tf_t[ii], sig_tf_z[ii], 0.0e0, 0.0e0, 0.0e0)

                svmyz = sigvm(sig_tf_r[ii], 0.0e0, sig_tf_z[ii], 0.0e0, 0.0e0, 0.0e0)
                sig_tf_vmises[ii] = max(svmxz, svmyz)

                # Maximum shear stress for the Tresca yield criterion using CEA calculation [Pa]
                s_shear_cea_tf_cond[ii] = (
                    1.02e0 * abs(sig_tf_r[ii]) + 1.6e0 * sig_tf_z[ii]
                )

        # ---
        # -----------------------------

        # Output formating : Maximum shear stress of each layer for the Tresca yield criterion
        # ----------------
        for ii in range(n_tf_layer):
            sig_max = 0.0e0
            ii_max = 0

            for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
                # CEA out of plane approximation
                if i_tf_tresca == 1 and i_tf_sup == 1 and ii >= i_tf_bucking + 1:
                    if sig_max < s_shear_cea_tf_cond[jj]:
                        sig_max = s_shear_cea_tf_cond[jj]
                        ii_max = jj

                # Conventional Tresca
                else:
                    if sig_max < s_shear_tf[jj]:
                        sig_max = s_shear_tf[jj]
                        ii_max = jj

            # OUT.DAT output

            sig_tf_r_max[ii] = sig_tf_r[ii_max]
            sig_tf_t_max[ii] = sig_tf_t[ii_max]
            sig_tf_z_max[ii] = sig_tf_z[ii_max]
            sig_tf_vmises_max[ii] = sig_tf_vmises[ii_max]

            # Maximum shear stress for the Tresca yield criterion (or CEA OOP correction)

            if i_tf_tresca == 1 and i_tf_sup == 1 and ii >= i_tf_bucking + 1:
                s_shear_tf_peak[ii] = s_shear_cea_tf_cond[ii_max]
            else:
                s_shear_tf_peak[ii] = s_shear_tf[ii_max]

        # Constraint equation for the Tresca yield criterion

        sig_tf_wp = s_shear_tf_peak[n_tf_bucking]
        # Maximum assumed in the first graded layer

        if i_tf_bucking >= 1:
            sig_tf_case = s_shear_tf_peak[n_tf_bucking - 1]
        if i_tf_bucking >= 2:
            sig_tf_cs_bucked = s_shear_tf_peak[0]
        # ----------------

        return (
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
            sig_tf_wp,
            sig_tf_case,
            sig_tf_cs_bucked,
            str_wp,
            casestr,
            insstrain,
            sig_tf_wp_av_z,
        )

    def out_stress(
        self,
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
    ):
        """Subroutine showing the writing the TF midplane stress analysis
        in the output file and the stress distribution in the SIG_TF.json
        file used to plot stress distributions
        Author : S. Kahn
        """

        def table_format_arrays(a, mult=1, delim="\t\t"):
            return delim.join([str(i * mult) for i in a])

        # Stress output section
        po.oheadr(self.outfile, "TF coils ")
        po.osubhd(self.outfile, "TF Coil Stresses (CCFE model) :")

        if tfcoil_variables.i_tf_stress_model == 1:
            po.ocmmnt(self.outfile, "Plane stress model with smeared properties")
        else:
            po.ocmmnt(self.outfile, "Generalized plane strain model")

        po.ovarre(
            self.outfile,
            "Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)",
            "(sig_tf_case_max)",
            tfcoil_variables.sig_tf_case_max,
        )

        po.ovarre(
            self.outfile,
            "Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)",
            "(sig_tf_wp_max)",
            tfcoil_variables.sig_tf_wp_max,
        )
        if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
            po.ocmmnt(
                self.outfile,
                "WP conduit Tresca criterion corrected using CEA formula (i_tf_tresca = 1)",
            )

        if tfcoil_variables.i_tf_bucking >= 3:
            po.ocmmnt(
                self.outfile, "No stress limit imposed on the TF-CS interface layer"
            )
            po.ocmmnt(
                self.outfile, "  -> Too much unknow on it material choice/properties"
            )

        # OUT.DAT data on maximum shear stress values for the Tresca criterion
        po.ocmmnt(
            self.outfile,
            "Materal stress of the point of maximum shear stress (Tresca criterion) for each layer",
        )
        po.ocmmnt(
            self.outfile,
            "Please use utilities/plot_stress_tf.py for radial plots plots summary",
        )

        if tfcoil_variables.i_tf_bucking == 0:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(self.outfile, "  Layers \t\t\t\t WP \t\t Outer case")
            else:
                po.write(self.outfile, "  Layers \t\t\t\t conductor \t\t Outer case")

        elif tfcoil_variables.i_tf_bucking == 1:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile, "  Layers \t\t\t\t Steel case \t\t WP \t\t Outer case"
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t bucking \t\t conductor \t\t Outer case",
                )

        elif tfcoil_variables.i_tf_bucking == 2:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t Steel case \t\t WP \t\t Outer case",
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t bucking \t\t conductor \t\t Outer case",
                )

        elif tfcoil_variables.i_tf_bucking == 3:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t interface \t\t Steel case \t\t WP \t\t Outer case",
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t interface \t\t bucking \t\t conductor \t\t Outer case",
                )

        po.write(
            self.outfile,
            f"  Radial stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_r_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Toroidal stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_t_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Vertical stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_z_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Von-Mises stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_vmises_max, 1e-6)}",
        )

        if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
            po.write(
                self.outfile,
                f"  Shear (CEA Tresca) \t\t\t (MPa) \t\t {table_format_arrays(s_shear_tf_peak, 1e-6)}",
            )
        else:
            po.write(
                self.outfile,
                f"  Shear (Tresca) \t\t\t (MPa) \t\t {table_format_arrays(s_shear_tf_peak, 1e-6)}",
            )

        po.write(self.outfile, "")
        po.write(
            self.outfile,
            f"  Toroidal modulus \t\t\t (GPa) \t\t {table_format_arrays(eyoung_trans, 1e-9)}",
        )
        po.write(
            self.outfile,
            f"  Vertical modulus \t\t\t (GPa) \t\t {table_format_arrays(eyoung_axial, 1e-9)}",
        )
        po.write(self.outfile, "")
        po.ovarre(
            self.outfile,
            "WP transverse modulus (GPa)",
            "(eyoung_wp_trans*1.0d-9)",
            eyoung_wp_trans * 1.0e-9,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP vertical modulus (GPa)",
            "(eyoung_wp_axial*1.0d-9)",
            eyoung_wp_axial * 1.0e-9,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP transverse Poissons ratio",
            "(poisson_wp_trans)",
            poisson_wp_trans,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP vertical-transverse Pois. rat.",
            "(poisson_wp_axial)",
            poisson_wp_axial,
            "OP ",
        )

        # MFILE.DAT data
        for ii in range(n_tf_bucking + 2):
            po.ovarre(
                constants.mfile,
                f"Radial    stress at maximum shear of layer {ii + 1} (Pa)",
                f"(sig_tf_r_max({ii + 1}))",
                sig_tf_r_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"toroidal  stress at maximum shear of layer {ii + 1} (Pa)",
                f"(sig_tf_t_max({ii + 1}))",
                sig_tf_t_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"Vertical  stress at maximum shear of layer {ii + 1} (Pa)",
                f"(sig_tf_z_max({ii + 1}))",
                sig_tf_z_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"Von-Mises stress at maximum shear of layer {ii + 1} (Pa)",
                f"(sig_tf_vmises_max({ii + 1}))",
                sig_tf_vmises_max[ii],
            )
            if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
                po.ovarre(
                    constants.mfile,
                    f"Maximum shear stress for CEA Tresca yield criterion {ii + 1} (Pa)",
                    f"(s_shear_tf_peak({ii + 1}))",
                    s_shear_tf_peak[ii],
                )
            else:
                po.ovarre(
                    constants.mfile,
                    f"Maximum shear stress for the Tresca yield criterion {ii + 1} (Pa)",
                    f"(s_shear_tf_peak({ii + 1}))",
                    s_shear_tf_peak[ii],
                )

        # SIG_TF.json storage
        sig_file_data = {
            "Points per layers": n_radial_array,
            "Radius (m)": radial_array,
            "Radial stress (MPa)": sig_tf_r * 1e-6,
            "Toroidal stress (MPa)": sig_tf_t * 1e-6,
            "Vertical stress (MPa)": sig_tf_z * 1e-6,
            "Radial smear stress (MPa)": sig_tf_smeared_r * 1e-6,
            "Toroidal smear stress (MPa)": sig_tf_smeared_t * 1e-6,
            "Vertical smear stress (MPa)": sig_tf_smeared_z * 1e-6,
            "Von-Mises stress (MPa)": sig_tf_vmises * 1e-6,
            "CEA Tresca stress (MPa)": (
                s_shear_cea_tf_cond * 1e-6
                if tfcoil_variables.i_tf_sup == 1
                else s_shear_tf * 1e-6
            ),
            "rad. displacement (mm)": deflect * 1e3,
        }
        if tfcoil_variables.i_tf_stress_model != 1:
            sig_file_data = {
                **sig_file_data,
                "Radial strain": str_tf_r,
                "Toroidal strain": str_tf_t,
                "Vertical strain": str_tf_z,
            }

        if tfcoil_variables.i_tf_sup == 1:
            sig_file_data = {
                **sig_file_data,
                "WP smeared stress (MPa)": sig_tf_wp_av_z * 1.0e-6,
            }

        sig_file_data = {
            k: list(v) if isinstance(v, np.ndarray) else v
            for k, v in sig_file_data.items()
        }

        sig_filename = (
            f2py_compatible_to_string(global_variables.output_prefix) + "SIG_TF.json"
        )
        with open(sig_filename, "w") as f:
            json.dump(sig_file_data, f)

        # TODO sig_tf_wp_av_z is always undefined here. This needs correcting or removing
        # if ( tfcoil_variables.i_tf_sup == 1 ) :
        # write(constants.sig_file,'(t2, "WP"    ," smeared stress", t20, "(MPa)",t26, *(F11.3,3x))') sig_tf_wp_av_z*1.0e-6
        #

        # Quantities from the plane stress stress formulation (no resitive coil output)
        if tfcoil_variables.i_tf_stress_model == 1 and tfcoil_variables.i_tf_sup == 1:
            # Other quantities (displacement strain, etc..)
            po.ovarre(
                self.outfile,
                "Maximum radial deflection at midplane (m)",
                "(deflect)",
                deflect[n_radial_array - 1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Vertical strain on casing",
                "(casestr)",
                tfcoil_variables.casestr,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Radial strain on insulator",
                "(insstrain)",
                tfcoil_variables.insstrain,
                "OP ",
            )


@numba.njit(cache=True)
def eyoung_parallel(
    eyoung_j_1, a_1, poisson_j_perp_1, eyoung_j_2, a_2, poisson_j_perp_2
):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in parallel with each other.
    The force goes in direction j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    This is pretty easy because the smeared properties are simply
    the average weighted by the cross-sectional areas perpendicular
    to j.
    The assumption is that the strains in j are equal.
    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_parallel(triplet1, triplet2, tripletOUT)
    call eyoung_parallel(triplet3, tripletOUT, tripletOUT)
    call eyoung_parallel(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """
    poisson_j_perp_3 = (poisson_j_perp_1 * a_1 + poisson_j_perp_2 * a_2) / (a_1 + a_2)
    eyoung_j_3 = (eyoung_j_1 * a_1 + eyoung_j_2 * a_2) / (a_1 + a_2)
    a_3 = a_1 + a_2

    return eyoung_j_3, a_3, poisson_j_perp_3


@numba.njit(cache=True)
def sigvm(sx: float, sy: float, sz: float, txy: float, txz: float, tyz: float) -> float:
    """Calculates Von Mises stress in a TF coil
    author: P J Knight, CCFE, Culham Science Centre
    author: B Reimer, FEDC
    This routine calculates the Von Mises combination of
    stresses (Pa) in a TF coil.

    :param sx: In-plane stress in X direction [Pa]
    :param sy: In-plane stress in Y direction [Pa]
    :param sz: In-plane stress in Z direction [Pa]
    :param txy: Out-of-plane stress in X-Y plane [Pa]
    :param txz: Out-of-plane stress in X-Z plane [Pa]
    :param tyz: Out-of-plane stress in Y-Z plane [Pa]

    :returns: Von Mises combination of stresses (Pa) in a TF coil.
    """

    return np.sqrt(
        0.5
        * (
            (sx - sy) ** 2
            + (sx - sz) ** 2
            + (sz - sy) ** 2
            + 6 * (txy**2 + txz**2 + tyz**2)
        )
    )


@numba.njit(cache=True, error_model="numpy")
def extended_plane_strain(
    nu_t,
    nu_zt,
    ey_t,
    ey_z,
    rad,
    d_curr,
    v_force,
    nlayers,
    n_radial_array,
    i_tf_bucking,
):
    """Author : C. Swanson, PPPL and S. Kahn, CCFE
    September 2021
    There is a writeup of the derivation of this model on the gitlab server.
    https://git.ccfe.ac.uk/process/process/-/issues/1414
    This surboutine estimates the radial displacement, stresses, and strains of
    the inboard midplane of the TF. It assumes that structures are axisymmetric
    and long in the axial (z) direction, the "axisymmetric extended plane strain"
    problem. The TF is assumed to be constructed from some number of layers,
    within which materials properties and current densities are uniform.
    The 1D radially-resolved solution is reduced to a 0D matrix inversion problem
    using analytic solutions to Lame's thick cylinder problem. Materials may be
    transverse-isotropic in Young's modulus and Poisson's ratio. The consraints
    are: Either zero radial stress or zero radial displacement at the inner
    surface (depending on whether the inner radius is zero), zero radial stress
    at the outer surface, total axial force (tension) is equal to the input,
    and optionally the axial force of an inner subset of cylinders is zero (slip
    conditions between the inner and outer subset). The matrix inversion / linear
    solve is always 4x4, no matter how many layers there are.
    The problem is formulated around a solution vector (A,B,eps_z,1.0,eps_z_slip)
    where A and B are the parameters in Lame's solution where u = A*r + B/r, u
    is the radial displacement. eps_z is the axial strain on the outer, force-
    carrying layers. eps_z_slip is the axial strain on the inner, non-force-
    carrying layers (optionally). The solution vector is A,B at the outermost
    radius, and is transformed via matrix multiplication into those A,B
    values at other radii. The constraints are inner products with this vector,
    and so when stacked form a matrix to invert.
    """
    # outputs
    sigr = np.zeros((n_radial_array * nlayers,))
    # Stress distribution in the radial direction (r) [Pa]

    sigt = np.zeros((n_radial_array * nlayers,))
    # Stress distribution in the toroidal direction (t) [Pa]

    sigz = np.zeros((n_radial_array * nlayers,))
    # Stress distribution in the vertical direction (z)

    str_r = np.zeros((n_radial_array * nlayers,))
    # Strain distribution in the radial direction (r)

    str_t = np.zeros((n_radial_array * nlayers,))
    # Strain distribution in the toroidal direction (t)

    str_z = np.zeros((n_radial_array * nlayers,))
    # Uniform strain in the vertical direction (z)

    r_deflect = np.zeros((n_radial_array * nlayers,))
    # Radial displacement radial distribution [m]

    rradius = np.zeros((n_radial_array * nlayers,))
    # Radius array [m]

    # local arrays
    # Stiffness form of compliance tensor
    nu_tz = np.zeros((nlayers,))
    # Transverse-axial Poisson's ratio
    # (ratio of axial strain to transverse strain upon transverse stress)
    ey_bar_z = np.zeros((nlayers,))
    # Axial effective Young's modulus (zero cross-strains, not stresses) [Pa]
    ey_bar_t = np.zeros((nlayers,))
    # Transverse effective Young's modulus [Pa]
    nu_bar_t = np.zeros((nlayers,))
    # Transverse effective Poisson's ratio
    nu_bar_tz = np.zeros((nlayers,))
    # Transverse-axial effective Poisson's ratio
    nu_bar_zt = np.zeros((nlayers,))
    # Axial-transverse effective Poisson's ratio

    # Lorentz force parameters
    currents = np.zeros((nlayers,))
    # Currents in each layer [A]
    currents_enclosed = np.zeros((nlayers,))
    # Currents enclosed by inner radius of each layer [A]
    f_lin_fac = np.zeros((nlayers,))
    # Factor that multiplies r linearly in the force density
    f_rec_fac = np.zeros((nlayers,))
    # Factor that multiplies r reciprocally in the force density
    f_int_a = np.zeros((nlayers,))
    # Force density integral that adds to Lame parameter A
    f_int_b = np.zeros((nlayers,))
    # Force density integral that adds to Lame parameter B

    # Layer transfer matrices
    m_int = np.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # outer radius to the inner radius of each layer
    m_ext = np.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # inner radius of one layer to the outer radius of the
    # next inner.
    m_tot = np.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # outer radius of the outer layer to the inner radius of
    # each.

    # Axial force inner product
    v_force_row = np.zeros(
        (
            1,
            5,
        ),
    )
    # Row vector (matrix multiplication is inner product) to
    # obtain the axial force from the force-carrying layers
    v_force_row_slip = np.zeros(
        (
            1,
            5,
        ),
    )
    # Row vector (matrix multiplication is inner product) to
    # obtain the axial force inner slip layers (no net force)
    rad_row_helper = np.zeros(
        (
            1,
            5,
        ),
    )
    # A helper variable to store [radius, 1, 0, 0, 0] in row

    # Boundary condition matrix
    m_bc = np.zeros(
        (
            4,
            5,
        ),
    )
    # Boundary condition matrix. Multiply this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain a zero vector.
    m_toinv = np.zeros(
        (
            4,
            4,
        ),
    )
    # Matrix to invert to get the solution
    RHS_vec = np.zeros((4,))
    # Right-hand-side vector to divide m_toinv
    a_vec_solution = np.zeros((5,))
    # Solution vector, Lame parameters at outer radius, strain
    # of force-carrying layers, and strain of slip layers
    # (A,B,eps_z,1,eps_z_slip)

    # Constructing the solution everywhere
    a_vec_layer = np.zeros((5,))
    # Lame parameters and strains vector at outer radius
    # of each layer

    # The stress calcualtion differential equations is analytically sloved
    # The final solution is given by the layer boundary conditions on
    # radial stress and displacement between layers solved
    # The problem is set as aa.cc = bb, cc being the constant we search

    # Inner slip layers parameters
    # Section 15 in the writeup
    # Innermost layer that takes axial force. Layers inner of this
    # have zero net axial force, to include CS decoupling.
    # This configuration JUST HAPPENS to work out because of
    # the specific i_tf_bucking options; if those are changed,
    # will need a switch here.
    nonslip_layer = i_tf_bucking

    if nonslip_layer < 1:
        nonslip_layer = 1

    # Stiffness tensor factors
    # Section 3 in the writeup
    # With Section 12 anisotropic materials properties
    # ***
    # Dependent Poisson's ratio: nu-transverse-axial
    # from nu-axial-transverse and the Young's moduli
    nu_tz[:] = nu_zt * ey_t / ey_z

    # Effective Young's Moduli and Poisson's ratios
    # holding strain, not stress, cross-terms constant
    ey_bar_z[:] = ey_z * (1 - nu_t) / (1 - nu_t - 2 * nu_tz * nu_zt)
    ey_bar_t[:] = (
        ey_t * (1 - nu_tz * nu_zt) / (1 - nu_t - 2 * nu_tz * nu_zt) / (1 + nu_t)
    )

    nu_bar_t[:] = (nu_t + nu_tz * nu_zt) / (1 - nu_tz * nu_zt)
    nu_bar_tz[:] = nu_tz / (1 - nu_t)
    nu_bar_zt[:] = nu_zt * (1 + nu_t) / (1 - nu_tz * nu_zt)

    # Lorentz force parameters
    # Section 13 in the writeup
    # ***
    # Currents in each layer [A]
    currents[:] = np.pi * d_curr * (rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2)
    # Currents enclosed by inner radius of each layer [A]
    currents_enclosed[0] = 0.0e0

    for ii in range(1, nlayers):
        currents_enclosed[ii] = currents_enclosed[ii - 1] + currents[ii - 1]
    # Factor that multiplies r linearly in the force density
    f_lin_fac[:] = RMU0 / 2.0e0 * d_curr**2
    # Factor that multiplies r reciprocally in the force density
    f_rec_fac[:] = (
        RMU0
        / 2.0e0
        * (d_curr * currents_enclosed / np.pi - d_curr**2 * rad[:nlayers] ** 2)
    )
    # Force density integral that adds to Lame parameter A
    f_int_a[:] = 0.5e0 * f_lin_fac * (
        rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2
    ) + f_rec_fac * np.log(rad[1 : nlayers + 1] / rad[:nlayers])
    if f_rec_fac[0] == 0e0:
        f_int_a[0] = 0.5e0 * f_lin_fac[0] * (rad[1] ** 2 - rad[0] ** 2)

    # Force density integral that adds to Lame parameter B
    f_int_b[:] = 0.25e0 * f_lin_fac * (
        rad[1 : nlayers + 1] ** 4 - rad[:nlayers] ** 4
    ) + 0.5e0 * f_rec_fac * (rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2)

    # Transformation matrix from outer to inner Lame parameters
    # Section 5 in the writeup
    # With Section 12 anisotropic materials properties
    # ***
    # m_int[kk] multiplies Lame parameter vector of layer kk (A,B,eps_z,1.0,eps_z_slip)
    # and transforms the values at the outer radius to the values at the inner radius
    for kk in range(nlayers):
        m_int[0, 0, kk] = 1.0e0
        m_int[1, 1, kk] = 1.0e0
        m_int[2, 2, kk] = 1.0e0
        m_int[3, 3, kk] = 1.0e0
        m_int[4, 4, kk] = 1.0e0

        m_int[0, 3, kk] = -0.5e0 / ey_bar_t[kk] * f_int_a[kk]
        m_int[1, 3, kk] = 0.5e0 / ey_bar_t[kk] * f_int_b[kk]

    # Transformation matrix between layers
    # Section 6 in the writeup
    # With Section 12 anisotropic materials properties
    # With Section 15 inner slip-decoupled layers
    # ***
    # m_ext[kk] multiplies Lame parameter vector of layer kk (A,B,eps_z,1.0,eps_z_slip)
    # and transforms the values at the inner radius to the values at the outer radius
    # of layer kk-1
    for kk in range(1, nonslip_layer - 1):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        m_ext[0, 2, kk] = 0.0e0
        m_ext[0, 4, kk] = 0.5e0 * (ey_fac * nu_bar_zt[kk] - nu_bar_zt[kk - 1])

    if nonslip_layer > 1:
        ey_fac = ey_bar_t[nonslip_layer - 1] / ey_bar_t[nonslip_layer - 2]
        m_ext[0, 2, nonslip_layer - 1] = 0.5e0 * ey_fac * nu_bar_zt[nonslip_layer - 1]
        m_ext[0, 4, nonslip_layer - 1] = 0.5e0 * (-nu_bar_zt[nonslip_layer - 2])

    for kk in range(nonslip_layer, nlayers):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        m_ext[0, 2, kk] = 0.5e0 * (ey_fac * nu_bar_zt[kk] - nu_bar_zt[kk - 1])
        m_ext[0, 4, kk] = 0.0e0

    for kk in range(1, nlayers):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        m_ext[0, 0, kk] = 0.5e0 * (ey_fac * (1 + nu_bar_t[kk]) + 1 - nu_bar_t[kk - 1])
        if rad[kk] > 0e0:
            m_ext[0, 1, kk] = (
                0.5e0
                / rad[kk] ** 2
                * (1 - nu_bar_t[kk - 1] - ey_fac * (1 - nu_bar_t[kk]))
            )

        m_ext[1, 0, kk] = rad[kk] ** 2 * (1 - m_ext[0, 0, kk])
        m_ext[1, 1, kk] = 1 - rad[kk] ** 2 * m_ext[0, 1, kk]
        m_ext[1, 2, kk] = -(rad[kk] ** 2) * m_ext[0, 2, kk]
        m_ext[1, 4, kk] = -(rad[kk] ** 2) * m_ext[0, 4, kk]
        m_ext[2, 2, kk] = 1.0e0
        m_ext[3, 3, kk] = 1.0e0
        m_ext[4, 4, kk] = 1.0e0

    # Total transformation matrix, from Lame parmeters at outside to
    # Lame parameters at inside of each layer
    # Section 7 in the writeup
    # ***
    m_tot[:, :, nlayers - 1] = m_int[:, :, nlayers - 1]

    for kk in range(nlayers - 2, -1, -1):
        m_tot[:, :, kk] = np.ascontiguousarray(m_int[:, :, kk]) @ (
            np.ascontiguousarray(m_ext[:, :, kk + 1])
            @ np.ascontiguousarray(m_tot[:, :, kk + 1])
        )

    # Axial force inner product. Dot-product this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain the axial force.
    # Section 8 in the writeup
    # ***
    # Axial stiffness products
    ey_bar_z_area = np.pi * sum(
        ey_bar_z[nonslip_layer - 1 : nlayers]
        * (
            rad[nonslip_layer : nlayers + 1] ** 2
            - rad[nonslip_layer - 1 : nlayers] ** 2
        )
    )
    ey_bar_z_area_slip = np.pi * sum(
        ey_bar_z[: nonslip_layer - 1]
        * (rad[1:nonslip_layer] ** 2 - rad[: nonslip_layer - 1] ** 2)
    )

    # Axial stiffness inner product, for layers which carry axial force
    rad_row_helper[0, :] = [rad[nlayers] ** 2, 1e0, 0e0, 0e0, 0e0]
    v_force_row[:, :] = (
        2e0 * np.pi * ey_bar_z[nlayers - 1] * nu_bar_tz[nlayers - 1] * rad_row_helper
    )
    rad_row_helper[0, :] = [rad[nonslip_layer - 1] ** 2, 1e0, 0e0, 0e0, 0e0]
    v_force_row[:, :] = v_force_row - 2e0 * np.pi * ey_bar_z[
        nonslip_layer - 1
    ] * nu_bar_tz[nonslip_layer - 1] * (
        rad_row_helper @ np.ascontiguousarray(m_tot[:, :, nonslip_layer - 1])
    )
    for kk in range(nonslip_layer, nlayers):
        rad_row_helper[0, :] = [rad[kk] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row[:, :] = v_force_row + 2e0 * np.pi * (
            ey_bar_z[kk - 1] * nu_bar_tz[kk - 1] - ey_bar_z[kk] * nu_bar_tz[kk]
        ) * (rad_row_helper @ np.ascontiguousarray(m_tot[:, :, kk]))

    # Include the effect of axial stiffness
    v_force_row[0, 2] += ey_bar_z_area

    # Axial stiffness inner product, for layers which DON'T carry force
    if nonslip_layer > 1:
        rad_row_helper[0, :] = [rad[nonslip_layer - 1] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row_slip[:, :] = (
            2e0
            * np.pi
            * ey_bar_z[nonslip_layer - 2]
            * nu_bar_tz[nonslip_layer - 2]
            * (rad_row_helper @ np.ascontiguousarray(m_tot[:, :, nonslip_layer - 1]))
        )
        rad_row_helper[0, :] = [rad[0] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row_slip[:, :] -= (
            2e0
            * np.pi
            * ey_bar_z[0]
            * nu_bar_tz[0]
            * (rad_row_helper @ np.ascontiguousarray(m_tot[:, :, 0]))
        )
        for kk in range(1, nonslip_layer - 1):
            rad_row_helper[0, :] = [rad[kk] ** 2, 1e0, 0e0, 0e0, 0e0]
            v_force_row_slip[:, :] += (
                2e0
                * np.pi
                * (ey_bar_z[kk - 1] * nu_bar_tz[kk - 1] - ey_bar_z[kk] * nu_bar_tz[kk])
                * (rad_row_helper @ np.ascontiguousarray(m_tot[:, :, kk]))
            )
        # Include the effect of axial stiffness
        v_force_row_slip[0, 4] += ey_bar_z_area_slip
    else:
        # If there's no inner slip layer, still need a finite 5th
        # element to ensure no singular matrix
        v_force_row_slip[0, :] = [0e0, 0e0, 0e0, 0e0, 1e0]

    # Boundary condition matrix. Multiply this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain a zero vector.
    # Solved to get the Lame parameters.
    # Section 9 in the writeup
    # ***
    # Outer boundary condition row, zero radial stress
    m_bc[0, :] = [
        (1e0 + nu_bar_t[nlayers - 1]) * rad[nlayers] ** 2,
        -1e0 + nu_bar_t[nlayers - 1],
        nu_bar_zt[nlayers - 1] * rad[nlayers] ** 2,
        0e0,
        0e0,
    ]
    # Inner boundary condition row, zero radial stress
    # or zero displacement if rad(1)=0
    if nonslip_layer > 1:
        m_bc[1, :] = [
            (1e0 + nu_bar_t[0]) * rad[0] ** 2,
            -1e0 + nu_bar_t[0],
            0e0,
            0e0,
            nu_bar_zt[0] * rad[0] ** 2,
        ]
    else:
        m_bc[1, :] = [
            (1e0 + nu_bar_t[0]) * rad[0] ** 2,
            -1e0 + nu_bar_t[0],
            nu_bar_zt[0] * rad[0] ** 2,
            0e0,
            0e0,
        ]

    m_bc[1, :] = np.ascontiguousarray(m_bc[1, :]) @ np.ascontiguousarray(m_tot[:, :, 0])
    # Axial force boundary condition
    m_bc[2, :] = v_force_row[0, :]
    m_bc[2, 3] = m_bc[2, 3] - v_force
    # Axial force boundary condition of slip layers
    m_bc[3, :] = v_force_row_slip[0, :]

    # The solution, the outermost Lame parameters A,B
    # and the axial strains of the force-carrying and
    # slip layers eps_z and eps_z_slip.
    # Section 10 in the writeup
    # ***
    m_toinv[:, :3] = m_bc[
        :,
        :3,
    ]
    m_toinv[:, 3] = m_bc[:, 4]
    RHS_vec[:] = -m_bc[:, 3]

    a_vec_solution[:4] = np.linalg.solve(m_toinv, RHS_vec)

    a_vec_solution[4] = a_vec_solution[3]
    a_vec_solution[3] = 1

    # Radial/toroidal/vertical stress radial distribution
    # ------
    # Radial displacement, stress and strain distributions

    a_vec_layer[:] = a_vec_solution[:]
    for ii in range(nlayers - 1, -1, -1):
        a_layer = a_vec_layer[0]
        b_layer = a_vec_layer[1]

        dradius = (rad[ii + 1] - rad[ii]) / (n_radial_array - 1)

        for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
            rradius[jj] = rad[ii] + dradius * (jj - (n_radial_array * ii))

            f_int_a_plot = 0.5e0 * f_lin_fac[ii] * (
                rad[ii + 1] ** 2 - rradius[jj] ** 2
            ) + f_rec_fac[ii] * np.log(rad[ii + 1] / (rradius[jj]))
            f_int_b_plot = 0.25e0 * f_lin_fac[ii] * (
                rad[ii + 1] ** 4 - rradius[jj] ** 4
            ) + 0.5e0 * f_rec_fac[ii] * (rad[ii + 1] ** 2 - rradius[jj] ** 2)
            a_plot = a_layer - 0.5e0 / ey_bar_t[ii] * f_int_a_plot
            b_plot = b_layer + 0.5e0 / ey_bar_t[ii] * f_int_b_plot

            # Radial displacement
            r_deflect[jj] = a_plot * rradius[jj] + b_plot / rradius[jj]

            # Radial strain
            str_r[jj] = a_plot - b_plot / rradius[jj] ** 2
            # Azimuthal strain
            str_t[jj] = a_plot + b_plot / rradius[jj] ** 2
            # Axial strain
            if ii < nonslip_layer - 1:
                str_z[jj] = a_vec_solution[4]
            else:
                str_z[jj] = a_vec_solution[2]

            # Radial stress
            sigr[jj] = ey_bar_t[ii] * (
                str_r[jj] + (nu_bar_t[ii] * str_t[jj]) + (nu_bar_zt[ii] * str_z[jj])
            )
            # Aximuthal stress
            sigt[jj] = ey_bar_t[ii] * (
                str_t[jj] + (nu_bar_t[ii] * str_r[jj]) + (nu_bar_zt[ii] * str_z[jj])
            )
            # Axial stress
            sigz[jj] = ey_bar_z[ii] * (
                str_z[jj] + (nu_bar_tz[ii] * (str_r[jj] + str_t[jj]))
            )

        a_vec_layer = np.ascontiguousarray(m_tot[:, :, ii]) @ a_vec_solution
        a_vec_layer = np.ascontiguousarray(m_ext[:, :, ii]) @ a_vec_layer
    # ------

    return rradius, sigr, sigt, sigz, str_r, str_t, str_z, r_deflect


@numba.njit(cache=True)
def plane_stress(nu, rad, ey, j, nlayers, n_radial_array):
    """Calculates the stresses in a superconductor TF coil
    inboard leg at the midplane using the plain stress approximation
    author: P J Knight, CCFE, Culham Science Centre
    author: J Morris, CCFE, Culham Science Centre
    author: S Kahn, CCFE, Culham Science Centre
    This routine calculates the stresses in a superconductor TF coil
    inboard leg at midplane.
    <P>A 2 layer plane stress model developed by CCFE is used. The first layer
    is the steel case inboard of the winding pack, and the second
    layer is the winding pack itself.
    PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
    """
    alpha = np.zeros((nlayers,))
    beta = np.zeros((nlayers,))
    # Lorentz body force parametres

    area = np.zeros((nlayers,))
    # Layer area

    aa = np.zeros((
        2 * nlayers,
        2 * nlayers,
    ))
    # Matrix encoding the integration constant cc coeficients

    bb = np.zeros((2 * nlayers,))
    # Vector encoding the alpha/beta (lorentz forces) contribution

    cc = np.zeros((2 * nlayers,))
    c1 = np.zeros((nlayers,))
    c2 = np.zeros((nlayers,))
    # Integration constants vector (solution)

    rradius = np.zeros((nlayers * n_radial_array,))
    # Radius array [m]

    sigr = np.zeros((nlayers * n_radial_array,))
    # Radial stress radial distribution [Pa]

    sigt = np.zeros((nlayers * n_radial_array,))
    # Toroidal stress radial distribution [Pa]

    r_deflect = np.zeros((nlayers * n_radial_array,))
    # Radial deflection (displacement) radial distribution [m]

    kk = ey / (1 - nu**2)

    # Lorentz forces parametrisation coeficients (array equation)
    alpha = 0.5e0 * RMU0 * j**2 / kk

    inner_layer_curr = 0.0e0
    for ii in range(nlayers):
        beta[ii] = (
            0.5e0
            * RMU0
            * j[ii]
            * (inner_layer_curr - np.pi * j[ii] * rad[ii] ** 2)
            / (np.pi * kk[ii])
        )

        # Layer area
        area[ii] = np.pi * (rad[ii + 1] ** 2 - rad[ii] ** 2)

        # Total current carried by the inners layers
        inner_layer_curr = inner_layer_curr + area[ii] * j[ii]
    # ***

    # Null radial stress at R(1)
    aa[0, 0] = kk[0] * (1.0e0 + nu[0])
    aa[0, 1] = -kk[0] * (1.0e0 - nu[0]) / (rad[0] ** 2)

    # Inter-layer boundary conditions
    if nlayers != 1:
        for ii in range(nlayers - 1):
            # Continuous radial normal stress at R(ii+1)
            aa[2 * ii + 1, 2 * ii] = kk[ii] * (1.0e0 + nu[ii])
            aa[2 * ii + 1, 2 * ii + 1] = -kk[ii] * (1.0e0 - nu[ii]) / rad[ii + 1] ** 2
            aa[2 * ii + 1, 2 * ii + 2] = -kk[ii + 1] * (1.0e0 + nu[ii + 1])
            aa[2 * ii + 1, 2 * ii + 3] = (
                kk[ii + 1] * (1.0e0 - nu[ii + 1]) / rad[ii + 1] ** 2
            )

            # Continuous displacement at R(ii+1)
            aa[2 * ii + 2, 2 * ii] = rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 1] = 1.0e0 / rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 2] = -rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 3] = -1.0e0 / rad[ii + 1]

    # Radial stress = 0
    aa[2 * (nlayers - 1) + 1, 2 * (nlayers - 1)] = kk[nlayers - 1] * (
        1.0e0 + nu[nlayers - 1]
    )
    aa[2 * (nlayers - 1) + 1, 2 * (nlayers - 1) + 1] = (
        -kk[nlayers - 1] * (1.0e0 - nu[nlayers - 1]) / rad[nlayers] ** 2
    )
    # ***

    # Right hand side vector bb
    # ***
    # Null radial stress at R(1)
    bb[0] = -kk[0] * (
        0.125e0 * alpha[0] * (3.0e0 + nu[0]) * rad[0] ** 2
        + 0.5e0 * beta[0] * (1.0e0 + (1.0e0 + nu[0]) * np.log(rad[0]))
    )

    # Inter-layer boundary conditions
    if nlayers != 1:
        for ii in range(nlayers - 1):
            # Continuous radial normal stress at R[ii+1]
            bb[2 * ii + 1] = -kk[ii] * (
                0.125e0 * alpha[ii] * (3.0e0 + nu[ii]) * rad[ii + 1] ** 2
                + 0.5e0 * beta[ii] * (1.0e0 + (1.0e0 + nu[ii]) * np.log(rad[ii + 1]))
            ) + kk[ii + 1] * (
                0.125e0 * alpha[ii + 1] * (3.0e0 + nu[ii + 1]) * rad[ii + 1] ** 2
                + 0.5e0
                * beta[ii + 1]
                * (1.0e0 + (1.0e0 + nu[ii + 1]) * np.log(rad[ii + 1]))
            )

            # Continuous displacement at R(ii+1)
            bb[2 * ii + 2] = (
                -0.125e0 * alpha[ii] * rad[ii + 1] ** 3
                - 0.5e0 * beta[ii] * rad[ii + 1] * np.log(rad[ii + 1])
                + 0.125e0 * alpha[ii + 1] * rad[ii + 1] ** 3
                + 0.5e0 * beta[ii + 1] * rad[ii + 1] * np.log(rad[ii + 1])
            )

    # Null radial stress at R(nlayers+1)
    bb[2 * (nlayers - 1) + 1] = -kk[nlayers - 1] * (
        0.125e0 * alpha[nlayers - 1] * (3.0e0 + nu[nlayers - 1]) * rad[nlayers] ** 2
        + 0.5e0
        * beta[nlayers - 1]
        * (1.0e0 + (1.0e0 + nu[nlayers - 1]) * np.log(rad[nlayers]))
    )
    # ***

    #  Find solution vector cc
    # ***
    aa = np.asfortranarray(aa)
    cc = np.linalg.solve(aa, bb)

    #  Multiply c by (-1) (John Last, internal CCFE memorandum, 21/05/2013)
    for ii in range(nlayers):
        c1[ii] = cc[2 * ii]
        c2[ii] = cc[2 * ii + 1]
    # ***
    # ------

    # Radial/toroidal/vertical stress radial distribution
    # ------

    for ii in range(nlayers):
        dradius = (rad[ii + 1] - rad[ii]) / n_radial_array
        for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
            rad_c = rad[ii] + dradius * (jj - n_radial_array * ii)
            rradius[jj] = rad_c

            # Radial stress radial distribution [Pa]
            sigr[jj] = kk[ii] * (
                (1.0e0 + nu[ii]) * c1[ii]
                - ((1.0e0 - nu[ii]) * c2[ii]) / rad_c**2
                + 0.125e0 * (3.0e0 + nu[ii]) * alpha[ii] * rad_c**2
                + 0.5e0 * beta[ii] * (1.0e0 + (1.0e0 + nu[ii]) * np.log(rad_c))
            )

            # Radial stress radial distribution [Pa]
            sigt[jj] = kk[ii] * (
                (1.0e0 + nu[ii]) * c1[ii]
                + (1.0e0 - nu[ii]) * c2[ii] / rad_c**2
                + 0.125e0 * (1.0e0 + 3.0e0 * nu[ii]) * alpha[ii] * rad_c**2
                + 0.5e0 * beta[ii] * (nu[ii] + (1.0e0 + nu[ii]) * np.log(rad_c))
            )

            #  Deflection [m]
            r_deflect[jj] = (
                c1[ii] * rad_c
                + c2[ii] / rad_c
                + 0.125e0 * alpha[ii] * rad_c**3
                + 0.5e0 * beta[ii] * rad_c * np.log(rad_c)
            )

    return sigr, sigt, r_deflect, rradius


@numba.njit(cache=True)
def eyoung_parallel_array(n, eyoung_j_in, a_in, poisson_j_perp_in):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in parallel with each other.
    The force goes in direction j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    This is pretty easy because the smeared properties are simply
    the average weighted by the cross-sectional areas perpendicular
    to j.
    The assumption is that the strains in j are equal.
    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_parallel(triplet1, triplet2, tripletOUT)
    call eyoung_parallel(triplet3, tripletOUT, tripletOUT)
    call eyoung_parallel(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """
    eyoung_j_out = 0
    a_out = 0
    poisson_j_perp_out = 0

    # Parallel-composite them all together
    for ii in range(n):
        eyoung_j_out, a_out, poisson_j_perp_out = eyoung_parallel(
            eyoung_j_in[ii],
            a_in[ii],
            poisson_j_perp_in[ii],
            eyoung_j_out,
            a_out,
            poisson_j_perp_out,
        )

    return eyoung_j_out, a_out, poisson_j_perp_out


@numba.njit(cache=True)
def eyoung_t_nested_squares(n, eyoung_j_in, l_in, poisson_j_perp_in):
    """
    Author : C. Swanson, PPPL
    January 2022
    This subroutine gives the smeared transverse elastic
    properties of n members whose cross sectional areas are
    nested squares. It uses the subroutines eyoung_series and
    eyoung_parallel, above, so please be aware of the assumptions
    inherent in those subroutines.

    It assumes that each "leg" of the square cross section
    (vertical slice, as described in Figure 10 of the TF coil
    documentation) is composed of several layers under stress in
    series, and each leg is under stress in parallel with every
    other leg.
    """
    eyoung_j_working = np.zeros((n,))
    l_working = np.zeros((n,))
    poisson_j_perp_working = np.zeros((n,))

    # First member
    eyoung_j_working[0] = eyoung_j_in[0]
    l_working[0] = l_in[0]
    poisson_j_perp_working[0] = poisson_j_perp_in[0]

    for ii in range(1, n):
        # Initialize the leg of which this is the new member
        eyoung_j_working[ii] = eyoung_j_in[ii]
        l_working[ii] = l_working[ii - 1] + l_in[ii]
        poisson_j_perp_working[ii] = poisson_j_perp_in[ii]

        # Serial-composite the new layer of this member into the previous legs
        # changed from range(ii-1) because range(0) == []
        for jj in range(ii):
            (
                eyoung_j_working[jj],
                l_working[jj],
                poisson_j_perp_working[jj],
            ) = eyoung_series(
                eyoung_j_working[ii],
                l_in[ii],
                poisson_j_perp_working[ii],
                eyoung_j_working[jj],
                l_working[jj],
                poisson_j_perp_working[jj],
            )

    # Find stiffest leg
    eyoung_stiffest = max(eyoung_j_working)

    eyoung_j_out = 0
    l_out = 0
    poisson_j_perp_out = 0

    # Parallel-composite them all together
    for ii in range(n):
        eyoung_j_out, l_out, poisson_j_perp_out = eyoung_parallel(
            eyoung_j_working[ii],
            l_in[ii],
            poisson_j_perp_working[ii],
            eyoung_j_out,
            l_out,
            poisson_j_perp_out,
        )

    return eyoung_j_out, l_out, poisson_j_perp_out, eyoung_stiffest


@numba.njit(cache=True)
def eyoung_series(eyoung_j_1, l_1, poisson_j_perp_1, eyoung_j_2, l_2, poisson_j_perp_2):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in series with each other.
    The force goes in direction j.
    The assumption is that the stresses in j are equal.
    The smeared Young's modulus is the inverse of the average of
    the inverse of the Young's moduli, weighted by the length
    of the members in j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    The smeared Poisson's ratio is the averaged of the Poisson's
    ratios, weighted by the quantity (Young's modulus / length of
    the members in j).

    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_series(triplet1, triplet2, tripletOUT)
    call eyoung_series(triplet3, tripletOUT, tripletOUT)
    call eyoung_series(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """

    if eyoung_j_1 * eyoung_j_2 == 0:
        # poisson_j_perp_3 = 0
        poisson_j_perp_3 = poisson_j_perp_1 if eyoung_j_1 == 0 else poisson_j_perp_2

        eyoung_j_3 = 0.0
        l_3 = l_1 + l_2
    else:
        poisson_j_perp_3 = (
            poisson_j_perp_1 * l_1 / eyoung_j_1 + poisson_j_perp_2 * l_2 / eyoung_j_2
        ) / (l_1 / eyoung_j_1 + l_2 / eyoung_j_2)
        eyoung_j_3 = (l_1 + l_2) / (l_1 / eyoung_j_1 + l_2 / eyoung_j_2)
        l_3 = l_1 + l_2

    eyoung_j_3 = np.array(eyoung_j_3)
    poisson_j_perp_3 = np.array(poisson_j_perp_3)

    return eyoung_j_3, l_3, poisson_j_perp_3


def init_tfcoil_variables():
    tfv.acasetf = 0.0
    tfv.acasetfo = 0.0
    tfv.a_tf_turn_steel = 0.0
    tfv.acond = 0.0
    tfv.a_tf_turn_cable_space = 0.0
    tfv.a_tf_turn_insulation = 0.0
    tfv.a_tf_coil_wp_turn_insulation = 0.0
    tfv.sig_tf_case_max = 6.0e8
    tfv.sig_tf_wp_max = 6.0e8
    tfv.a_tf_leg_outboard = 0.0
    tfv.aswp = 0.0
    tfv.avwp = 0.0
    tfv.a_tf_wp_coolant_channels = 0.0
    tfv.bcritsc = 24.0
    tfv.b_tf_inboard_peak = 0.0
    tfv.bmaxtfrp = 0.0
    tfv.casestr = 0.0
    tfv.dr_tf_plasma_case = 0.0
    tfv.f_dr_tf_plasma_case = 0.05
    tfv.i_f_dr_tf_plasma_case = False
    tfv.dx_tf_side_case = 0.0
    tfv.casths_fraction = 0.06
    tfv.t_conductor = 0.0
    tfv.t_cable_tf = 0.0
    tfv.t_cable_tf_is_input = False
    tfv.t_turn_tf = 0.0
    tfv.t_turn_tf_is_input = False
    tfv.f_t_turn_tf = 1.0
    tfv.t_turn_tf_max = 0.05
    tfv.acs = 0.0
    tfv.cdtfleg = 0.0
    tfv.cforce = 0.0
    tfv.cplen = 0.0
    tfv.c_tf_turn = 7.0e4
    tfv.cpttf_max = 9.0e4
    tfv.dcase = 8000.0
    tfv.dcond = [6080.0, 6080.0, 6070.0, 6080.0, 6080.0, 8500.0, 6070.0, 8500.0, 8500.0]
    tfv.dcondins = 1800.0
    tfv.dia_tf_turn_coolant_channel = 0.005
    tfv.estotftgj = 0.0
    tfv.b_crit_upper_nbti = 14.86
    tfv.t_crit_nbti = 9.04
    tfv.max_force_density = 0.0
    tfv.fcutfsu = 0.69
    tfv.fhts = 0.5
    tfv.insstrain = 0.0
    tfv.i_tf_stress_model = 1
    tfv.i_tf_tresca = 0
    tfv.i_tf_wp_geom = -1
    tfv.i_tf_case_geom = 0
    tfv.i_tf_turns_integer = 0
    tfv.i_tf_sc_mat = 1
    tfv.i_tf_sup = 1
    tfv.i_tf_shape = 0
    tfv.i_tf_cond_eyoung_axial = 0
    tfv.i_tf_cond_eyoung_trans = 1
    tfv.n_pancake = 10
    tfv.n_layer = 20
    tfv.n_rad_per_layer = 100
    tfv.i_tf_bucking = -1
    tfv.n_tf_graded_layers = 1
    tfv.n_tf_stress_layers = 0
    tfv.n_tf_wp_layers = 5
    tfv.j_tf_bus = 1.25e6
    tfv.j_crit_str_tf = 0.0
    tfv.j_crit_str_0 = [
        596905475.80390120,
        1925501534.8512938,
        724544682.96063495,
        549858624.45072436,
        669284509.85818779,
        0.0,
        898964415.36996782,
        1158752995.2559297,
        865652122.9071957,
    ]
    tfv.j_tf_wp_critical = 0.0
    tfv.jwdgpro = 0.0
    tfv.j_tf_wp = 0.0
    tfv.oacdcp = 0.0
    tfv.eyoung_ins = 1.0e8
    tfv.eyoung_steel = 2.05e11
    tfv.eyoung_cond_axial = 6.6e8
    tfv.eyoung_cond_trans = 0.0
    tfv.eyoung_res_tf_buck = 150.0e9
    tfv.eyoung_copper = 117.0e9
    tfv.eyoung_al = 69.0e9
    tfv.poisson_steel = 0.3
    tfv.poisson_copper = 0.35
    tfv.poisson_al = 0.35
    tfv.poisson_ins = 0.34
    tfv.poisson_cond_axial = 0.3
    tfv.poisson_cond_trans = 0.3
    tfv.r_b_tf_inboard_peak = 0.0
    tfv.res_tf_leg = 0.0
    tfv.toroidalgap = 1.0  # [m]
    tfv.ftoroidalgap = 1.0
    tfv.ripmax = 1.0
    tfv.ripple = 0.0
    tfv.c_tf_total = 0.0
    tfv.radial_array[:] = 0.0
    tfv.sig_tf_r[:] = 0.0
    tfv.sig_tf_t[:] = 0.0
    tfv.deflect[:] = 0.0
    tfv.sig_tf_z = 0.0
    tfv.sig_tf_vmises[:] = 0.0
    tfv.s_shear_tf[:] = 0.0
    tfv.sig_tf_cs_bucked = 0.0
    tfv.sig_tf_case = 0.0
    tfv.sig_tf_wp = 0.0
    tfv.str_cs_con_res = -0.005
    tfv.str_pf_con_res = -0.005
    tfv.str_tf_con_res = -0.005
    tfv.str_wp = 0.0
    tfv.str_wp_max = 0.7e-2
    tfv.i_str_wp = 1
    tfv.quench_model = string_to_f2py_compatible(tfv.quench_model, "exponential")
    tfv.time1 = 0
    tfv.tcritsc = 16.0
    tfv.tdmptf = 10.0
    tfv.a_tf_coil_inboard = 0.0
    tfv.len_tf_bus = 300.0
    tfv.m_tf_bus = 0.0
    tfv.tfckw = 0.0
    tfv.tfcmw = 0.0
    tfv.p_cp_resistive_mw = 0.0
    tfv.p_tf_joints_resistive_mw = 0.0
    tfv.tfcryoarea = 0.0
    tfv.tficrn = 0.0
    tfv.ind_tf_coil = 0.0
    tfv.tfinsgap = 0.01
    tfv.p_tf_leg_resistive_mw = 0.0
    tfv.rho_cp = 0.0
    tfv.rho_tf_leg = 0.0
    tfv.rho_tf_bus = 1.86e-8
    tfv.frhocp = 1.0
    tfv.frholeg = 1.0
    tfv.rho_tf_joints = 2.5e-10
    tfv.n_tf_joints_contact = 6
    tfv.n_tf_joints = 4
    tfv.th_joint_contact = 0.03
    tfv.p_tf_joints_resistive = 0.0
    tfv.len_tf_coil = 0.0
    tfv.eff_tf_cryo = -1.0
    tfv.n_tf_coils = 16.0
    tfv.tfocrn = 0.0
    tfv.tfsai = 0.0
    tfv.tfsao = 0.0
    tfv.tftmp = 4.5
    tfv.dx_tf_inboard_out_toroidal = 1.0
    tfv.dx_tf_turn_insulation = 8e-4
    tfv.layer_ins = 0.0
    tfv.dr_tf_nose_case = 0.3
    tfv.dr_tf_wp = 0.0
    tfv.dx_tf_turn_steel = 8e-3
    tfv.tinstf = 0.018
    tfv.tmargmin_tf = 0.0
    tfv.tmargmin_cs = 0.0
    tfv.tmargmin = 0.0
    tfv.temp_margin = 0.0
    tfv.tmargtf = 0.0
    tfv.tmaxpro = 150.0
    tfv.tmax_croco = 200.0
    tfv.croco_quench_temperature = 0.0
    tfv.temp_tf_cryo = 4.5
    tfv.n_tf_turn = 0.0
    tfv.vdalw = 20.0
    tfv.vforce = 0.0
    tfv.f_vforce_inboard = 0.5
    tfv.vforce_outboard = 0.0
    tfv.vftf = 0.4
    tfv.voltfleg = 0.0
    tfv.vtfkv = 0.0
    tfv.vtfskv = 0.0
    tfv.whtcas = 0.0
    tfv.whtcon = 0.0
    tfv.whtconcu = 0.0
    tfv.whtconal = 0.0
    tfv.whtconin = 0.0
    tfv.whtconsc = 0.0
    tfv.m_tf_turn_steel = 0.0
    tfv.whtgw = 0.0
    tfv.m_tf_coils_total = 0.0
    tfv.wwp1 = 0.0
    tfv.wwp2 = 0.0
    tfv.dthet[:] = 0.0
    tfv.radctf[:] = 0.0
    tfv.r_tf_arc[:] = 0.0
    tfv.xctfc[:] = 0.0
    tfv.z_tf_arc[:] = 0.0
    tfv.yctfc[:] = 0.0
    tfv.tfa[:] = 0.0
    tfv.tfb[:] = 0.0
    tfv.drtop = 0.0
    tfv.dztop = 0.0
    tfv.etapump = 0.8
    tfv.fcoolcp = 0.3
    tfv.f_a_tf_cool_outboard = 0.2
    tfv.a_cp_cool = 0.0
    tfv.ncool = 0.0
    tfv.p_cp_coolant_pump_elec = 0.0
    tfv.p_cp_resistive = 0.0
    tfv.p_tf_leg_resistive = 0.0
    tfv.ptempalw = 473.15  # 200 C
    tfv.rcool = 0.005
    tfv.tcoolin = 313.15  # 40 C
    tfv.dtiocool = 0.0
    tfv.temp_cp_average = 373.15  # 100 C
    tfv.tcpav2 = 0.0
    tfv.temp_tf_legs_outboard = -1.0
    tfv.tcpmax = 0.0
    tfv.vcool = 20.0
    tfv.vol_cond_cp = 0.0
    tfv.whtcp = 0.0
    tfv.whttflgs = 0.0
    tfv.tfc_sidewall_is_fraction = False
    tfv.i_cp_joints = -1
    tfv.cryo_cool_req = 0.0
    tfv.theta1_coil = 45.0
    tfv.theta1_vv = 1.0  # 1 Deg
    tfv.max_vv_stress = 143.0e6
