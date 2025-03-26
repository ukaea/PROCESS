import logging

import numpy as np

from process.fortran import (
    build_variables,
    constants,
    error_handling,
    physics_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.tfcoil import TFCoil

logger = logging.getLogger(__name__)

logger = logging.getLogger(__name__)


class ResistiveTFCoil(TFCoil):
    def __init__(self):
        self.outfile = constants.nout

    def run(self):
        self.res_tf_internal_geom()

    def res_tf_internal_geom(self):
        """
        Author : S. Kahn
        Resisitve TF turn geometry, equivalent to winding_pack subroutines

        """
        sctfcoil_module.r_wp_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.thkcas
        )
        sctfcoil_module.r_wp_outer = (
            build_variables.r_tf_inboard_out - tfcoil_variables.casthi
        )

        # Conductor layer radial thickness at centercollumn top [m]
        if physics_variables.itart == 1:
            sctfcoil_module.dr_tf_wp_top = (
                build_variables.r_cp_top
                - tfcoil_variables.casthi
                - tfcoil_variables.thkcas
                - build_variables.r_tf_inboard_in
            )

        # Number of turns
        # Set by user (no turn structure by default, i.e. tfcoil_variables.n_tf_turn = 1 )
        if (
            abs(tfcoil_variables.n_tf_turn)
            < np.finfo(float(tfcoil_variables.n_tf_turn)).eps
        ):
            tfcoil_variables.n_tf_turn = 1.0e0

        # Total mid-plane cross-sectional area of winding pack, [m2]
        # including the surrounding ground-wall insulation layer
        sctfcoil_module.awpc = (
            np.pi
            * (sctfcoil_module.r_wp_outer**2 - sctfcoil_module.r_wp_inner**2)
            / tfcoil_variables.n_tf_coils
        )

        # Area of the front case, the plasma-facing case of the inner TF coil [m2]
        sctfcoil_module.a_case_front = (
            np.pi
            * (
                (sctfcoil_module.r_wp_outer + tfcoil_variables.casthi) ** 2
                - sctfcoil_module.r_wp_outer**2
            )
            / tfcoil_variables.n_tf_coils
        )

        # WP mid-plane cross-section excluding ground insulation per coil [m2]
        sctfcoil_module.awptf = np.pi * (
            (sctfcoil_module.r_wp_outer - tfcoil_variables.tinstf) ** 2
            - (sctfcoil_module.r_wp_inner + tfcoil_variables.tinstf) ** 2
        ) / tfcoil_variables.n_tf_coils - 2.0e0 * tfcoil_variables.tinstf * (
            tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tinstf
        )

        # Ground insulation cross-section area per coil [m2]
        sctfcoil_module.a_ground_ins = sctfcoil_module.awpc - sctfcoil_module.awptf

        # Exact mid-plane cross-section area of the conductor per TF coil [m2]
        a_tf_cond = np.pi * (
            (
                sctfcoil_module.r_wp_outer
                - tfcoil_variables.tinstf
                - tfcoil_variables.thicndut
            )
            ** 2
            - (
                sctfcoil_module.r_wp_inner
                + tfcoil_variables.tinstf
                + tfcoil_variables.thicndut
            )
            ** 2
        ) / tfcoil_variables.n_tf_coils - (
            tfcoil_variables.dr_tf_wp
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.thicndut)
        ) * 2.0e0 * (
            tfcoil_variables.tinstf
            + tfcoil_variables.thicndut * tfcoil_variables.n_tf_turn
        )
        a_tf_cond = a_tf_cond * (1.0e0 - tfcoil_variables.fcoolcp)

        # Inter turn insulation area per coil [m2]
        tfcoil_variables.aiwp = sctfcoil_module.awptf - a_tf_cond / (
            1.0e0 - tfcoil_variables.fcoolcp
        )

        # Total insulation cross-section per coil [m2]
        sctfcoil_module.a_tf_ins = tfcoil_variables.aiwp + sctfcoil_module.a_ground_ins

        # Insulation fraction [-]
        sctfcoil_module.f_tf_ins = (
            tfcoil_variables.n_tf_coils
            * sctfcoil_module.a_tf_ins
            / tfcoil_variables.tfareain
        )

        # Total cross-sectional area of the bucking cylindre and the outer support
        # support structure per coil [m2]
        # physics_variables.itart = 1 : Only valid at mid-plane
        tfcoil_variables.acasetf = (
            tfcoil_variables.tfareain / tfcoil_variables.n_tf_coils
        ) - sctfcoil_module.awpc

        # Current per turn
        tfcoil_variables.cpttf = tfcoil_variables.c_tf_total / (
            tfcoil_variables.n_tf_turn * tfcoil_variables.n_tf_coils
        )

        # Exact current density on TF oubard legs
        tfcoil_variables.cdtfleg = tfcoil_variables.c_tf_total / (
            (1.0e0 - tfcoil_variables.fcoolcp)
            * (
                tfcoil_variables.tftort
                - 2.0e0
                * (
                    tfcoil_variables.n_tf_turn * tfcoil_variables.thicndut
                    + tfcoil_variables.tinstf
                )
            )
            * (
                build_variables.dr_tf_outboard
                - 2.0e0 * (tfcoil_variables.thicndut + tfcoil_variables.tinstf)
            )
        )

        # Reporting negative WP areas issues
        if sctfcoil_module.awpc < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awpc
            error_handling.fdiags[0] = tfcoil_variables.dr_tf_wp
            error_handling.report_error(99)

        elif sctfcoil_module.awptf < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awptf
            error_handling.report_error(101)
