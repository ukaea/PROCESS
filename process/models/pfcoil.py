import logging
import math

import numba
import numpy as np
from scipy import optimize
from scipy.linalg import svd
from scipy.special import ellipe, ellipk

import process.models.superconductors as superconductors
from process import process_output as op
from process.core import constants
from process.core.exceptions import ProcessValueError
from process.data_structure import build_variables as bv
from process.data_structure import constraint_variables as ctv
from process.data_structure import cs_fatigue_variables as csfv
from process.data_structure import fwbs_variables as fwbsv
from process.data_structure import (
    pfcoil_variables,
    superconducting_tf_coil_variables,
)
from process.data_structure import physics_variables as pv
from process.data_structure import rebco_variables as rcv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv

logger = logging.getLogger(__name__)


class PFCoil:
    """Calculate poloidal field coil system parameters."""

    def __init__(self, cs_fatigue):
        """Initialise Fortran module variables."""
        self.outfile = constants.NOUT  # output file unit
        self.mfile = constants.MFILE  # mfile file unit
        pfcoil_variables.init_pfcoil_module()
        self.cs_fatigue = cs_fatigue
        self.cs_coil = CSCoil(cs_fatigue)

    def run(self):
        """Run the PF coil model."""
        self.pfcoil()

        # Poloidal field coil inductance calculation
        self.induct(False)

        # Volt-second capability of PF coil set
        self.vsec()

    def output(self):
        """Output results to output file."""
        self.cs_coil.output_cs_structure()
        self.outpf()
        self.outvolt()

    def output_induct(self):
        """Output poloidal field coil inductance calculation."""
        self.induct(True)

    def pfcoil(self):
        """Routine to perform calculations for the PF and Central Solenoid coils.

        This subroutine performs the calculations for the PF and
        Central Solenoid coils, to determine their size, location, current waveforms,
        stresses etc.
        """
        lrow1 = 2 * pfcoil_variables.NPTSMX + pfcoil_variables.N_PF_GROUPS_MAX
        lcol1 = pfcoil_variables.N_PF_GROUPS_MAX

        pcls0 = np.zeros(pfcoil_variables.N_PF_GROUPS_MAX, dtype=int)
        ncls0 = np.zeros(pfcoil_variables.N_PF_GROUPS_MAX + 2, dtype=int)

        pfcoil_variables.rcls0, pfcoil_variables.zcls0 = np.zeros(
            (
                2,
                pfcoil_variables.N_PF_GROUPS_MAX,
                pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
            ),
            order="F",
        )
        pfcoil_variables.ccls0 = np.zeros(int(pfcoil_variables.N_PF_GROUPS_MAX / 2))
        brin, bzin, rpts, zpts = np.zeros((4, pfcoil_variables.NPTSMX))
        bfix, bvec = np.zeros((2, lrow1))
        gmat, _, _ = np.zeros((3, lrow1, lcol1), order="F")
        signn = np.zeros(2)
        aturn = np.zeros(pfcoil_variables.NGC2)

        # Toggle switch for i_pf_location()=2 coils above/below midplane
        top_bottom = 1

        # Set up the number of PF coils including the Central Solenoid (n_cs_pf_coils),
        # and the number of PF circuits including the plasma (n_pf_cs_plasma_circuits)
        if pfcoil_variables.n_pf_coil_groups > pfcoil_variables.N_PF_GROUPS_MAX:
            raise ProcessValueError(
                "n_pf_coil_groups is larger than n_pf_groups_max",
                n_pf_coil_groups=pfcoil_variables.n_pf_coil_groups,
                n_pf_groups_max=pfcoil_variables.N_PF_GROUPS_MAX,
            )

        # Total the number of PF coils in all groups, and check that none
        # exceeds the limit
        pfcoil_variables.n_cs_pf_coils = 0
        for i in range(pfcoil_variables.n_pf_coil_groups):
            if (
                pfcoil_variables.n_pf_coils_in_group[i]
                > pfcoil_variables.N_PF_COILS_IN_GROUP_MAX
            ):
                raise ProcessValueError(
                    "PFCOIL: Too many coils in a PF coil group",
                    i=i,
                    n_pf_coils_in_group=pfcoil_variables.n_pf_coils_in_group[i],
                    n_pf_coils_in_group_max=pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
                )

            pfcoil_variables.n_cs_pf_coils = (
                pfcoil_variables.n_cs_pf_coils + pfcoil_variables.n_pf_coils_in_group[i]
            )

        # Add one if an Central Solenoid is present, and make an extra group
        if bv.iohcl != 0:
            pfcoil_variables.n_cs_pf_coils = pfcoil_variables.n_cs_pf_coils + 1
            pfcoil_variables.n_pf_coils_in_group[pfcoil_variables.n_pf_coil_groups] = 1

        # Add one for the plasma
        pfcoil_variables.n_pf_cs_plasma_circuits = pfcoil_variables.n_cs_pf_coils + 1

        # Overall current density in the Central Solenoid at beginning of pulse
        pfcoil_variables.j_cs_pulse_start = (
            pfcoil_variables.j_cs_flat_top_end
            * pfcoil_variables.f_j_cs_start_pulse_end_flat_top
        )

        # Set up array of times
        tv.t_pulse_cumulative[0] = 0.0e0
        tv.t_pulse_cumulative[1] = tv.t_plant_pulse_coil_precharge
        tv.t_pulse_cumulative[2] = (
            tv.t_pulse_cumulative[1] + tv.t_plant_pulse_plasma_current_ramp_up
        )
        tv.t_pulse_cumulative[3] = (
            tv.t_pulse_cumulative[2] + tv.t_plant_pulse_fusion_ramp
        )
        tv.t_pulse_cumulative[4] = tv.t_pulse_cumulative[3] + tv.t_plant_pulse_burn
        tv.t_pulse_cumulative[5] = (
            tv.t_pulse_cumulative[4] + tv.t_plant_pulse_plasma_current_ramp_down
        )

        # Set up call to MHD scaling routine for coil currents.
        # First break up Central Solenoid solenoid into 'filaments'

        (
            pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.z_pf_coil_lower[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_cs_middle,
            pfcoil_variables.z_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.a_cs_poloidal,
            pfcoil_variables.dz_cs_full,
            pfcoil_variables.dr_cs_full,
        ) = self.cs_coil.calculate_cs_geometry(
            z_tf_inside_half=bv.z_tf_inside_half,
            f_z_cs_tf_internal=pfcoil_variables.f_z_cs_tf_internal,
            dr_cs=bv.dr_cs,
            dr_bore=bv.dr_bore,
        )

        # nfxf is the total no of filaments into which the Central Solenoid is split,
        # if present
        if bv.iohcl == 0:
            pfcoil_variables.nfxf = 0
            c_cs_flat_top_end = 0.0e0
        else:
            pfcoil_variables.nfxf = 2 * pfcoil_variables.n_cs_current_filaments

            # total Central Solenoid current at EOF
            c_cs_flat_top_end = -(
                pfcoil_variables.a_cs_poloidal * pfcoil_variables.j_cs_flat_top_end
            )

            if pfcoil_variables.nfxf > pfcoil_variables.NFIXMX:
                raise ProcessValueError(
                    "Too many filaments nfxf repesenting the OH coil",
                    nfxf=pfcoil_variables.nfxf,
                    nfixmx=pfcoil_variables.NFIXMX,
                )

            # Symmetric up/down Central Solenoid : Find (R,Z) and current of each filament at BOP

            (
                pfcoil_variables.r_pf_cs_current_filaments,
                pfcoil_variables.z_pf_cs_current_filaments,
                pfcoil_variables.c_pf_cs_current_filaments,
            ) = self.cs_coil.place_cs_filaments(
                n_cs_current_filaments=pfcoil_variables.n_cs_current_filaments,
                r_cs_middle=pfcoil_variables.r_cs_middle,
                z_cs_inside_half=pfcoil_variables.dz_cs_full / 2,
                c_cs_flat_top_end=c_cs_flat_top_end,
                f_j_cs_start_pulse_end_flat_top=pfcoil_variables.f_j_cs_start_pulse_end_flat_top,
                nfxf=pfcoil_variables.nfxf,
            )

        # Scale PF coil locations
        signn[0] = 1.0e0
        signn[1] = -1.0e0
        pfcoil_variables.r_pf_outside_tf_midplane = (
            superconducting_tf_coil_variables.r_tf_outboard_out
            + pfcoil_variables.dr_pf_tf_outboard_out_offset
        )

        # Place the PF coils:

        # N.B. Problems here if coil=n_pf_coils_in_group(group) is greater than 2.
        for group in range(pfcoil_variables.n_pf_coil_groups):
            if pfcoil_variables.i_pf_location[group] == 1:
                # PF coil is stacked on top of the Central Solenoid
                # Use a helper function to compute r_pf_coil_middle_group_array and
                # z_pf_coil_middle_group_array arrays for this group

                r_pf_coil_middle_group_array, z_pf_coil_middle_group_array = (
                    self.place_pf_above_cs(
                        n_pf_coils_in_group=pfcoil_variables.n_pf_coils_in_group,
                        n_pf_group=group,
                        r_cs_middle=pfcoil_variables.r_cs_middle,
                        dr_pf_cs_middle_offset=pfcoil_variables.dr_pf_cs_middle_offset,
                        z_tf_inside_half=bv.z_tf_inside_half,
                        dr_tf_inboard=bv.dr_tf_inboard,
                        z_cs_coil_upper=pfcoil_variables.dz_cs_full / 2,
                    )
                )
                for coil in range(pfcoil_variables.n_pf_coils_in_group[group]):
                    pfcoil_variables.r_pf_coil_middle_group_array[group, coil] = (
                        r_pf_coil_middle_group_array[group, coil]
                    )
                    pfcoil_variables.z_pf_coil_middle_group_array[group, coil] = (
                        z_pf_coil_middle_group_array[group, coil]
                    )

            # =========================================================================

            elif pfcoil_variables.i_pf_location[group] == 2:
                # PF coil is on top of the TF coil
                (
                    r_pf_coil_middle_group_array,
                    z_pf_coil_middle_group_array,
                    top_bottom,
                ) = self.place_pf_above_tf(
                    n_pf_coils_in_group=pfcoil_variables.n_pf_coils_in_group,
                    n_pf_group=group,
                    rmajor=pv.rmajor,
                    triang=pv.triang,
                    rminor=pv.rminor,
                    itart=pv.itart,
                    itartpf=pv.itartpf,
                    z_tf_inside_half=bv.z_tf_inside_half,
                    dz_tf_upper_lower_midplane=bv.dz_tf_upper_lower_midplane,
                    z_tf_top=bv.z_tf_top,
                    top_bottom=top_bottom,
                    rpf2=pfcoil_variables.rpf2,
                    zref=pfcoil_variables.zref,
                )

                for coil in range(pfcoil_variables.n_pf_coils_in_group[group]):
                    pfcoil_variables.r_pf_coil_middle_group_array[group, coil] = (
                        r_pf_coil_middle_group_array[group, coil]
                    )
                    pfcoil_variables.z_pf_coil_middle_group_array[group, coil] = (
                        z_pf_coil_middle_group_array[group, coil]
                    )
            # =========================================================================

            elif pfcoil_variables.i_pf_location[group] == 3:
                # PF coil is radially outside the TF coil
                (
                    r_pf_coil_middle_group_array,
                    z_pf_coil_middle_group_array,
                ) = self.place_pf_outside_tf(
                    n_pf_coils_in_group=pfcoil_variables.n_pf_coils_in_group,
                    n_pf_group=group,
                    rminor=pv.rminor,
                    zref=pfcoil_variables.zref,
                    i_tf_shape=tfv.i_tf_shape,
                    i_r_pf_outside_tf_placement=pfcoil_variables.i_r_pf_outside_tf_placement,
                    r_pf_outside_tf_midplane=pfcoil_variables.r_pf_outside_tf_midplane,
                )

                for coil in range(pfcoil_variables.n_pf_coils_in_group[group]):
                    pfcoil_variables.r_pf_coil_middle_group_array[group, coil] = (
                        r_pf_coil_middle_group_array[group, coil]
                    )
                    pfcoil_variables.z_pf_coil_middle_group_array[group, coil] = (
                        z_pf_coil_middle_group_array[group, coil]
                    )

            # =========================================================================

            elif pfcoil_variables.i_pf_location[group] == 4:
                (
                    r_pf_coil_middle_group_array,
                    z_pf_coil_middle_group_array,
                ) = self.place_pf_generally(
                    n_pf_coils_in_group=pfcoil_variables.n_pf_coils_in_group,
                    n_pf_group=group,
                    rminor=pv.rminor,
                    rmajor=pv.rmajor,
                    zref=pfcoil_variables.zref,
                    rref=pfcoil_variables.rref,
                )

                for coil in range(pfcoil_variables.n_pf_coils_in_group[group]):
                    pfcoil_variables.r_pf_coil_middle_group_array[group, coil] = (
                        r_pf_coil_middle_group_array[group, coil]
                    )
                    pfcoil_variables.z_pf_coil_middle_group_array[group, coil] = (
                        z_pf_coil_middle_group_array[group, coil]
                    )

            else:
                raise ProcessValueError(
                    "Illegal i_pf_location value",
                    group=group,
                    i_pf_location=pfcoil_variables.i_pf_location[group],
                )
            # =========================================================================

        # Allocate current to the PF coils:
        # "Flux swing coils" participate in cancellation of the CS
        # field during a flux swing. "Equilibrium coils" are varied
        # to create the equilibrium field, targeting the correct
        # vertical field
        # As implemented, all coils are flux swing coils
        # As implemented, Location 3 and 4 coils are equilibrium
        # coils.

        # Flux swing coils:
        if pfcoil_variables.j_cs_pulse_start != 0.0e0:
            # Find currents for plasma initiation to null field across plasma
            npts = 32  # Number of test points across plasma midplane
            if npts > pfcoil_variables.NPTSMX:
                raise ProcessValueError(
                    "Too many test points npts across plasma midplane",
                    npts=npts,
                    nptsmx=pfcoil_variables.NPTSMX,
                )

            # Position and B-field at each test point
            drpt = 2.0e0 * pv.rminor / (npts - 1)
            rpt0 = pv.rmajor - pv.rminor

            for i in range(npts):
                rpts[i] = rpt0 + (i) * drpt
                zpts[i] = 0.0e0
                brin[i] = 0.0e0
                bzin[i] = 0.0e0

                # Calculate currents in coils to produce the given field
            pfcoil_variables.ssq0, pfcoil_variables.ccl0 = self.efc(
                npts,
                rpts,
                zpts,
                brin,
                bzin,
                pfcoil_variables.nfxf,
                pfcoil_variables.r_pf_cs_current_filaments,
                pfcoil_variables.z_pf_cs_current_filaments,
                pfcoil_variables.c_pf_cs_current_filaments,
                pfcoil_variables.n_pf_coil_groups,
                pfcoil_variables.n_pf_coils_in_group,
                pfcoil_variables.r_pf_coil_middle_group_array,
                pfcoil_variables.z_pf_coil_middle_group_array,
                pfcoil_variables.alfapf,
                bfix,
                gmat,
                bvec,
            )

        # Equilibrium coil currents determined by SVD targeting B
        if pfcoil_variables.i_pf_current == 1:
            # Simple coil current scaling for STs (good only for A < about 1.8)
            # Bypasses SVD solver
            if pv.itart == 1 and pv.itartpf == 0:
                for i in range(pfcoil_variables.n_pf_coil_groups):
                    if pfcoil_variables.i_pf_location[i] == 1:
                        # PF coil is stacked on top of the Central Solenoid
                        pfcoil_variables.ccls[i] = 0.0e0
                        raise ProcessValueError(
                            "i_pf_location(i) should not be 1 if itart=1", i=i
                        )

                    if pfcoil_variables.i_pf_location[i] == 2:
                        # PF coil is on top of the TF coil
                        pfcoil_variables.ccls[i] = (
                            0.3e0 * pv.aspect**1.6e0 * pv.plasma_current
                        )

                    elif pfcoil_variables.i_pf_location[i] == 3:
                        # PF coil is radially outside the TF coil
                        pfcoil_variables.ccls[i] = -0.4e0 * pv.plasma_current

                    else:
                        raise ProcessValueError(
                            "Illegal value of i_pf_location(i)",
                            i=i,
                            i_pf_location=pfcoil_variables.i_pf_location[i],
                        )

                # Vertical field (T)
                pv.b_plasma_vertical_required = (
                    -1.0e-7
                    * pv.plasma_current
                    / pv.rmajor
                    * (
                        math.log(8.0e0 * pv.aspect)
                        + pv.beta_poloidal_vol_avg
                        + (pv.ind_plasma_internal_norm / 2.0e0)
                        - 1.5e0
                    )
                )

            else:
                # Conventional aspect ratio scaling
                nfxf0 = 0
                ngrp0 = 0
                nocoil = 0
                for i in range(pfcoil_variables.n_pf_coil_groups):
                    if pfcoil_variables.i_pf_location[i] == 1:
                        # Do not allow if no central solenoid
                        if bv.iohcl == 0:
                            raise ProcessValueError(
                                "i_pf_location(i) should not be 1 if iohcl=0"
                            )
                        # PF coil is stacked on top of the Central Solenoid
                        # This coil is to balance Central Solenoid flux and should not be involved
                        # in equilibrium calculation -- RK 07/12
                        pfcoil_variables.ccls[i] = 0.0e0
                        nfxf0 = nfxf0 + pfcoil_variables.n_pf_coils_in_group[i]
                        for ccount in range(pfcoil_variables.n_pf_coils_in_group[i]):
                            pfcoil_variables.r_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.r_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.z_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.z_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.c_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.ccls[i]
                            )
                            nocoil = nocoil + 1

                    elif pfcoil_variables.i_pf_location[i] == 2:
                        # PF coil is on top of the TF coil; divertor coil
                        # This is a fixed current for this calculation -- RK 07/12

                        pfcoil_variables.ccls[i] = (
                            pv.plasma_current
                            * 2.0e0
                            * (
                                1.0e0
                                - (pv.kappa * pv.rminor)
                                / abs(
                                    pfcoil_variables.z_pf_coil_middle_group_array[i, 0]
                                )
                            )
                        )
                        nfxf0 = nfxf0 + pfcoil_variables.n_pf_coils_in_group[i]
                        for ccount in range(pfcoil_variables.n_pf_coils_in_group[i]):
                            pfcoil_variables.r_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.r_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.z_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.z_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.c_pf_cs_current_filaments[nocoil] = (
                                pfcoil_variables.ccls[i]
                            )
                            nocoil = nocoil + 1

                    elif pfcoil_variables.i_pf_location[i] == 3:
                        # PF coil is radially outside the TF coil
                        # This is an equilibrium coil, current must be solved for

                        pcls0[ngrp0] = i + 1
                        ngrp0 = ngrp0 + 1

                    elif pfcoil_variables.i_pf_location[i] == 4:
                        # PF coil is generally placed
                        # See issue 1418
                        # https://git.ccfe.ac.uk/process/process/-/issues/1418
                        # This is an equilibrium coil, current must be solved for

                        pcls0[ngrp0] = i + 1
                        ngrp0 = ngrp0 + 1

                    else:
                        raise ProcessValueError(
                            "Illegal value of i_pf_location(i)",
                            i=i,
                            i_pf_location=pfcoil_variables.i_pf_location[i],
                        )

                for ccount in range(ngrp0):
                    ncls0[ccount] = 2
                    pfcoil_variables.rcls0[ccount, 0] = (
                        pfcoil_variables.r_pf_coil_middle_group_array[
                            pcls0[ccount] - 1, 0
                        ]
                    )
                    pfcoil_variables.rcls0[ccount, 1] = (
                        pfcoil_variables.r_pf_coil_middle_group_array[
                            pcls0[ccount] - 1, 1
                        ]
                    )
                    pfcoil_variables.zcls0[ccount, 0] = (
                        pfcoil_variables.z_pf_coil_middle_group_array[
                            pcls0[ccount] - 1, 0
                        ]
                    )
                    pfcoil_variables.zcls0[ccount, 1] = (
                        pfcoil_variables.z_pf_coil_middle_group_array[
                            pcls0[ccount] - 1, 1
                        ]
                    )

                npts0 = 1
                rpts[0] = pv.rmajor
                zpts[0] = 0.0e0
                brin[0] = 0.0e0

                # Added pv.ind_plasma_internal_norm term correctly -- RK 07/12

                bzin[0] = (
                    -1.0e-7
                    * pv.plasma_current
                    / pv.rmajor
                    * (
                        math.log(8.0e0 * pv.aspect)
                        + pv.beta_poloidal_vol_avg
                        + (pv.ind_plasma_internal_norm / 2.0e0)
                        - 1.5e0
                    )
                )

                pv.b_plasma_vertical_required = bzin[0]

                _ssqef, pfcoil_variables.ccls0 = self.efc(
                    npts0,
                    rpts,
                    zpts,
                    brin,
                    bzin,
                    nfxf0,
                    pfcoil_variables.r_pf_cs_current_filaments,
                    pfcoil_variables.z_pf_cs_current_filaments,
                    pfcoil_variables.c_pf_cs_current_filaments,
                    ngrp0,
                    ncls0,
                    pfcoil_variables.rcls0,
                    pfcoil_variables.zcls0,
                    pfcoil_variables.alfapf,
                    bfix,
                    gmat,
                    bvec,
                )

                for ccount in range(ngrp0):
                    pfcoil_variables.ccls[pcls0[ccount] - 1] = pfcoil_variables.ccls0[
                        ccount
                    ]

        # Flux swing from vertical field

        # If this is the first visit to the routine the inductance matrix
        # ind_pf_cs_plasma_mutual and the turns array have not yet been calculated, so we set
        # them to (very) approximate values to avoid strange behaviour...
        if pfcoil_variables.first_call:
            pfcoil_variables.ind_pf_cs_plasma_mutual[:, :] = 1.0e0
            pfcoil_variables.n_pf_coil_turns[:] = 100.0e0
            pfcoil_variables.first_call = False

        pfflux = 0.0e0
        nocoil = 0
        for ccount in range(pfcoil_variables.n_pf_coil_groups):
            for _i in range(pfcoil_variables.n_pf_coils_in_group[ccount]):
                pfflux = pfflux + (
                    pfcoil_variables.ccls[ccount]
                    * pfcoil_variables.ind_pf_cs_plasma_mutual[
                        nocoil, pfcoil_variables.n_pf_cs_plasma_circuits - 1
                    ]
                    / pfcoil_variables.n_pf_coil_turns[nocoil]
                )
                nocoil = nocoil + 1

        # Flux swing required from CS coil
        csflux = -(pv.vs_plasma_res_ramp + pv.vs_plasma_ind_ramp) - pfflux

        if bv.iohcl == 1:
            # Required current change in CS coil

            # Proposed new calculation...
            # dics = csflux / ind_pf_cs_plasma_mutual(n_cs_pf_coils,n_pf_cs_plasma_circuits)
            # BUT... ind_pf_cs_plasma_mutual(n_cs_pf_coils,n_pf_cs_plasma_circuits) is around 2000 times ddics below...

            ddics = (
                4.0e-7
                * np.pi
                * np.pi
                * (
                    (bv.dr_bore * bv.dr_bore)
                    + (bv.dr_cs * bv.dr_cs) / 6.0e0
                    + (bv.dr_cs * bv.dr_bore) / 2.0e0
                )
                / (bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * 2.0e0)
            )
            dics = csflux / ddics

            pfcoil_variables.f_j_cs_start_end_flat_top = (
                (-c_cs_flat_top_end * pfcoil_variables.f_j_cs_start_pulse_end_flat_top)
                + dics
            ) / c_cs_flat_top_end
            if np.abs(pfcoil_variables.f_j_cs_start_end_flat_top) > 1.0:
                logger.warning(
                    "Ratio of central solenoid overall current density at "
                    "beginning of flat-top / end of flat-top > 1 (|f_j_cs_start_end_flat_top| > 1)"
                )
        else:
            dics = 0.0e0
            pfcoil_variables.f_j_cs_start_end_flat_top = 1.0e0
            logger.error("OH coil not present; check volt-second calculations...")

        # Split groups of coils into one set containing ncl coils
        ncl = 0
        for nng in range(pfcoil_variables.n_pf_coil_groups):
            for ng2 in range(pfcoil_variables.n_pf_coils_in_group[nng]):
                pfcoil_variables.r_pf_coil_middle[ncl] = (
                    pfcoil_variables.r_pf_coil_middle_group_array[nng, ng2]
                )
                pfcoil_variables.z_pf_coil_middle[ncl] = (
                    pfcoil_variables.z_pf_coil_middle_group_array[nng, ng2]
                )

                # Currents at different times:

                # If PF coil currents are computed, not input via ccl0_ma, ccls_ma:
                # Then set ccl0_ma,ccls_ma from the computed pfcoil_variables.ccl0,pfcoil_variables.ccls
                if pfcoil_variables.i_pf_current != 0:
                    pfcoil_variables.ccl0_ma[nng] = 1.0e-6 * pfcoil_variables.ccl0[nng]
                    pfcoil_variables.ccls_ma[nng] = 1.0e-6 * pfcoil_variables.ccls[nng]
                else:
                    # Otherwise set pfcoil_variables.ccl0,pfcoil_variables.ccls via the input ccl0_ma and ccls_ma
                    pfcoil_variables.ccl0[nng] = 1.0e6 * pfcoil_variables.ccl0_ma[nng]
                    pfcoil_variables.ccls[nng] = 1.0e6 * pfcoil_variables.ccls_ma[nng]

                # Beginning of pulse: t = tv.t_plant_pulse_coil_precharge
                pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ncl] = (
                    1.0e-6 * pfcoil_variables.ccl0[nng]
                )

                # Beginning of flat-top: t = tv.t_plant_pulse_coil_precharge+tv.t_plant_pulse_plasma_current_ramp_up
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ncl] = 1.0e-6 * (
                    pfcoil_variables.ccls[nng]
                    - (
                        pfcoil_variables.ccl0[nng]
                        * pfcoil_variables.f_j_cs_start_end_flat_top
                        / pfcoil_variables.f_j_cs_start_pulse_end_flat_top
                    )
                )

                # End of flat-top: t = tv.t_plant_pulse_coil_precharge+tv.t_plant_pulse_plasma_current_ramp_up+tv.t_plant_pulse_fusion_ramp+tv.t_plant_pulse_burn
                pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * (
                    pfcoil_variables.ccls[nng]
                    - (
                        pfcoil_variables.ccl0[nng]
                        * (1.0e0 / pfcoil_variables.f_j_cs_start_pulse_end_flat_top)
                    )
                )

                ncl = ncl + 1

        # Current in Central Solenoid as a function of time
        # N.B. If the Central Solenoid is not present then c_cs_flat_top_end is zero.
        pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ncl] = (
            -1.0e-6
            * c_cs_flat_top_end
            * pfcoil_variables.f_j_cs_start_pulse_end_flat_top
        )
        pfcoil_variables.c_pf_cs_coil_flat_top_ma[ncl] = (
            1.0e-6 * c_cs_flat_top_end * pfcoil_variables.f_j_cs_start_end_flat_top
        )
        pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * c_cs_flat_top_end

        # Set up coil current waveforms, normalised to the peak current in
        # each coil
        self.waveform()  # sets c_pf_cs_coils_peak_ma(), f_c_pf_cs_peak_time_array()

        # Calculate PF coil geometry, current and number of turns
        # Dimensions are those of the winding pack, and exclude
        # the steel supporting case
        i = 0
        pfcoil_variables.r_pf_coil_outer_max = 0.0e0

        dz = 0

        for ii in range(pfcoil_variables.n_pf_coil_groups):
            for _ij in range(pfcoil_variables.n_pf_coils_in_group[ii]):
                if pfcoil_variables.i_pf_location[ii] == 1:
                    # PF coil is stacked on top of the Central Solenoid
                    dx = 0.5e0 * bv.dr_cs
                    dz = 0.5e0 * (
                        bv.z_tf_inside_half
                        * (1.0e0 - pfcoil_variables.f_z_cs_tf_internal)
                        + bv.dr_tf_inboard
                        + 0.1e0
                    )  # ???
                    area = 4.0e0 * dx * dz * pfcoil_variables.pf_current_safety_factor

                    # Number of turns
                    # c_pf_coil_turn_peak_input[i] is the current per turn (input)
                    pfcoil_variables.n_pf_coil_turns[i] = abs(
                        (pfcoil_variables.c_pf_cs_coils_peak_ma[i] * 1.0e6)
                        / pfcoil_variables.c_pf_coil_turn_peak_input[i]
                    )
                    aturn[i] = area / pfcoil_variables.n_pf_coil_turns[i]

                    # Actual winding pack current density
                    pfcoil_variables.j_pf_coil_wp_peak[i] = (
                        1.0e6 * abs(pfcoil_variables.c_pf_cs_coils_peak_ma[i]) / area
                    )

                    # Location of edges of each coil:
                    # r_pf_coil_inner = inner radius, r_pf_coil_outer = outer radius
                    # z_pf_coil_lower = 'lower' edge z (i.e. edge nearer to midplane)
                    # z_pf_coil_upper = 'upper' edge z (i.e. edge further from midplane)
                    pfcoil_variables.r_pf_coil_inner[i] = (
                        pfcoil_variables.r_pf_coil_middle[i] - dx
                    )
                    pfcoil_variables.r_pf_coil_outer[i] = (
                        pfcoil_variables.r_pf_coil_middle[i] + dx
                    )

                    pfcoil_variables.z_pf_coil_lower[i] = (
                        pfcoil_variables.z_pf_coil_middle[i] - dz
                    )
                    if pfcoil_variables.z_pf_coil_middle[i] < 0.0e0:
                        pfcoil_variables.z_pf_coil_lower[i] = (
                            pfcoil_variables.z_pf_coil_middle[i] + dz
                        )

                    pfcoil_variables.z_pf_coil_upper[i] = (
                        pfcoil_variables.z_pf_coil_middle[i] + dz
                    )

                    if pfcoil_variables.z_pf_coil_middle[i] < 0.0e0:
                        pfcoil_variables.z_pf_coil_upper[i] = (
                            pfcoil_variables.z_pf_coil_middle[i] - dz
                        )

                else:
                    # Other coils. N.B. Current density j_pf_coil_wp_peak[i] is defined in
                    # routine INITIAL for these coils.
                    area = (
                        abs(
                            pfcoil_variables.c_pf_cs_coils_peak_ma[i]
                            * 1.0e6
                            / pfcoil_variables.j_pf_coil_wp_peak[i]
                        )
                        * pfcoil_variables.pf_current_safety_factor
                    )

                    pfcoil_variables.n_pf_coil_turns[i] = abs(
                        (pfcoil_variables.c_pf_cs_coils_peak_ma[i] * 1.0e6)
                        / pfcoil_variables.c_pf_coil_turn_peak_input[i]
                    )
                    aturn[i] = area / pfcoil_variables.n_pf_coil_turns[i]

                    dx = 0.5e0 * math.sqrt(area)  # square cross-section

                    pfcoil_variables.r_pf_coil_inner[i] = (
                        pfcoil_variables.r_pf_coil_middle[i] - dx
                    )
                    pfcoil_variables.r_pf_coil_outer[i] = (
                        pfcoil_variables.r_pf_coil_middle[i] + dx
                    )

                    pfcoil_variables.z_pf_coil_lower[i] = (
                        pfcoil_variables.z_pf_coil_middle[i] - dx
                    )
                    if pfcoil_variables.z_pf_coil_middle[i] < 0.0e0:
                        pfcoil_variables.z_pf_coil_lower[i] = (
                            pfcoil_variables.z_pf_coil_middle[i] + dx
                        )

                    pfcoil_variables.z_pf_coil_upper[i] = (
                        pfcoil_variables.z_pf_coil_middle[i] + dx
                    )
                    if pfcoil_variables.z_pf_coil_middle[i] < 0.0e0:
                        pfcoil_variables.z_pf_coil_upper[i] = (
                            pfcoil_variables.z_pf_coil_middle[i] - dx
                        )

                # Outside radius of largest PF coil (m)
                pfcoil_variables.r_pf_coil_outer_max = max(
                    pfcoil_variables.r_pf_coil_outer_max,
                    pfcoil_variables.r_pf_coil_outer[i],
                )

                i = i + 1

        # Calculate peak field, allowable current density, resistive
        # power losses and volumes and weights for each PF coil, index i
        i = 0
        it = 0
        pfcoil_variables.p_pf_coil_resistive_total_flat_top = 0.0e0
        pfcoil_variables.m_pf_coil_max = 0.0e0

        for ii in range(pfcoil_variables.n_pf_coil_groups):
            iii = ii
            for ij in range(pfcoil_variables.n_pf_coils_in_group[ii]):
                # Peak field

                if ij == 0:
                    # Index args +1ed
                    _bri, _bro, _bzi, _bzo = peak_b_field_at_pf_coil(
                        n_coil=i + 1, n_coil_group=iii + 1, t_b_field_peak=it
                    )  # returns b_pf_coil_peak, bpf2

                # Issue 1871.  MDK
                # Allowable current density (for superconducting coils) for each coil, index i
                if pfcoil_variables.i_pf_conductor == 0:
                    bmax = max(
                        abs(pfcoil_variables.b_pf_coil_peak[i]),
                        abs(pfcoil_variables.bpf2[i]),
                    )

                    pfcoil_variables.j_pf_wp_critical[i], _jstrand, jsc, _tmarg = (
                        superconpf(
                            bmax,
                            pfcoil_variables.f_a_pf_coil_void[i],
                            pfcoil_variables.fcupfsu,
                            pfcoil_variables.j_pf_coil_wp_peak[i],
                            pfcoil_variables.i_pf_superconductor,
                            tfv.fhts,
                            tfv.str_pf_con_res,
                            tfv.tftmp,
                            tfv.bcritsc,
                            tfv.tcritsc,
                        )
                    )

                    # Strand critical current calculation for costing in $/kAm
                    # = superconducting filaments jc * (1 - strand copper fraction)
                    if pfcoil_variables.i_cs_superconductor in {2, 6, 8}:
                        pfcoil_variables.j_crit_str_pf = jsc
                    else:
                        pfcoil_variables.j_crit_str_pf = jsc * (
                            1 - pfcoil_variables.fcupfsu
                        )

                # Length of conductor

                rll = (
                    2.0e0
                    * np.pi
                    * pfcoil_variables.r_pf_coil_middle[i]
                    * pfcoil_variables.n_pf_coil_turns[i]
                )

                # Resistive coils

                if pfcoil_variables.i_pf_conductor == 1:
                    # Coil resistance (f_a_pf_coil_void is the void fraction)

                    respf = (
                        pfcoil_variables.rho_pf_coil
                        * rll
                        / (aturn[i] * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i]))
                    )

                    # Sum resistive power losses

                    pfcoil_variables.p_pf_coil_resistive_total_flat_top = (
                        pfcoil_variables.p_pf_coil_resistive_total_flat_top
                        + respf
                        * (
                            1.0e6
                            * pfcoil_variables.c_pf_cs_coil_pulse_start_ma[i]
                            / pfcoil_variables.n_pf_coil_turns[i]
                        )
                        ** 2
                    )

                # Winding pack volume

                volpf = aturn[i] * rll

                # Conductor weight (f_a_pf_coil_void is the void fraction)

                if pfcoil_variables.i_pf_conductor == 0:
                    pfcoil_variables.m_pf_coil_conductor[i] = (
                        volpf
                        * tfv.dcond[pfcoil_variables.i_pf_superconductor - 1]
                        * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i])
                    )
                else:
                    pfcoil_variables.m_pf_coil_conductor[i] = (
                        volpf
                        * constants.den_copper
                        * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i])
                    )

                # (J x B) force on coil

                forcepf = (
                    0.5e6
                    * (pfcoil_variables.b_pf_coil_peak[i] + pfcoil_variables.bpf2[i])
                    * abs(pfcoil_variables.c_pf_cs_coils_peak_ma[i])
                    * pfcoil_variables.r_pf_coil_middle[i]
                )

                # Stress ==> cross-sectional area of supporting steel to use

                if pfcoil_variables.i_pf_conductor == 0:
                    # Superconducting coil
                    # Updated assumptions: 500 MPa stress limit with all of the force
                    # supported in the conduit (steel) case.
                    # Now, 500 MPa replaced by sigpfcalw, sigpfcf now defaultly set to 1

                    areaspf = (
                        pfcoil_variables.sigpfcf
                        * forcepf
                        / (pfcoil_variables.sigpfcalw * 1.0e6)
                    )

                    # Thickness of hypothetical steel casing assumed to encase the PF
                    # winding pack; in reality, the steel is distributed
                    # throughout the conductor. Issue #152
                    # Assume a case of uniform thickness around coil cross-section
                    # Thickness found via a simple quadratic equation

                    drpdz = (
                        pfcoil_variables.r_pf_coil_outer[i]
                        - pfcoil_variables.r_pf_coil_inner[i]
                        + abs(
                            pfcoil_variables.z_pf_coil_upper[i]
                            - pfcoil_variables.z_pf_coil_lower[i]
                        )
                    )  # dr + dz
                    pfcoil_variables.pfcaseth[i] = 0.25e0 * (
                        -drpdz + math.sqrt(drpdz * drpdz + 4.0e0 * areaspf)
                    )

                else:
                    areaspf = 0.0e0  # Resistive coil - no steel needed
                    pfcoil_variables.pfcaseth[i] = 0.0e0

                # Weight of steel case

                pfcoil_variables.m_pf_coil_structure[i] = (
                    areaspf
                    * 2.0e0
                    * np.pi
                    * pfcoil_variables.r_pf_coil_middle[i]
                    * fwbsv.den_steel
                )

                # Mass of heaviest PF coil (tonnes)

                pfcoil_variables.m_pf_coil_max = max(
                    pfcoil_variables.m_pf_coil_max,
                    (
                        1.0e-3
                        * (
                            pfcoil_variables.m_pf_coil_conductor[i]
                            + pfcoil_variables.m_pf_coil_structure[i]
                        )
                    ),
                )
                i = i + 1

        # Find sum of current x turns x radius for all coils for 2015 costs model
        c = 0
        pfcoil_variables.itr_sum = 0.0e0
        for m in range(pfcoil_variables.n_pf_coil_groups):
            for _n in range(pfcoil_variables.n_pf_coils_in_group[m]):
                pfcoil_variables.itr_sum = pfcoil_variables.itr_sum + (
                    pfcoil_variables.r_pf_coil_middle[c]
                    * pfcoil_variables.n_pf_coil_turns[c]
                    * pfcoil_variables.c_pf_coil_turn_peak_input[c]
                )
                c = c + 1

        pfcoil_variables.itr_sum = pfcoil_variables.itr_sum + (
            (bv.dr_bore + 0.5 * bv.dr_cs)
            * pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1]
            * pfcoil_variables.c_pf_coil_turn_peak_input[
                pfcoil_variables.n_cs_pf_coils - 1
            ]
        )

        # Find Central Solenoid information
        if bv.iohcl != 0:
            self.cs_coil.ohcalc()

        # Summation of weights and current
        pfcoil_variables.m_pf_coil_conductor_total = 0.0e0
        pfcoil_variables.m_pf_coil_structure_total = 0.0e0
        pfcoil_variables.ricpf = 0.0e0

        for i in range(pfcoil_variables.n_cs_pf_coils):
            pfcoil_variables.m_pf_coil_conductor_total = (
                pfcoil_variables.m_pf_coil_conductor_total
                + pfcoil_variables.m_pf_coil_conductor[i]
            )
            pfcoil_variables.m_pf_coil_structure_total = (
                pfcoil_variables.m_pf_coil_structure_total
                + pfcoil_variables.m_pf_coil_structure[i]
            )
            pfcoil_variables.ricpf = pfcoil_variables.ricpf + abs(
                pfcoil_variables.c_pf_cs_coils_peak_ma[i]
            )

        # Plasma size and shape
        pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils] = (
            pv.rminor * pv.kappa
        )
        pfcoil_variables.z_pf_coil_lower[pfcoil_variables.n_cs_pf_coils] = (
            -pv.rminor * pv.kappa
        )
        pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils] = (
            pv.rmajor - pv.rminor
        )
        pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils] = (
            pv.rmajor + pv.rminor
        )
        pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils] = 1.0e0

        # Generate coil currents as a function of time using
        # user-provided waveforms etc. (c_pf_coil_turn_peak_input, f_j_cs_start_pulse_end_flat_top, f_j_cs_start_end_flat_top)
        for k in range(6):  # time points
            for i in range(pfcoil_variables.n_pf_cs_plasma_circuits - 1):
                pfcoil_variables.c_pf_coil_turn[i, k] = (
                    pfcoil_variables.f_c_pf_cs_peak_time_array[i, k]
                    * math.copysign(
                        pfcoil_variables.c_pf_coil_turn_peak_input[i],
                        pfcoil_variables.c_pf_cs_coils_peak_ma[i],
                    )
                )

        # Plasma wave form
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 0
        ] = 0.0e0
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 1
        ] = 0.0e0
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 2
        ] = pv.plasma_current
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 3
        ] = pv.plasma_current
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 4
        ] = pv.plasma_current
        pfcoil_variables.c_pf_coil_turn[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1, 5
        ] = 0.0e0

    def place_pf_above_cs(
        self,
        n_pf_coils_in_group: np.ndarray,
        n_pf_group: int,
        r_cs_middle: float,
        dr_pf_cs_middle_offset: float,
        z_tf_inside_half: float,
        dr_tf_inboard: float,
        z_cs_coil_upper: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the placement of PF coils stacked above the Central Solenoid.

        Parameters
        ----------
        n_pf_coils_in_group : np.ndarray
            Array containing the number of coils in each PF group.
        n_pf_group : int
            Index of the PF coil group.
        r_cs_middle : float
            Radial coordinate of CS coil centre (m).
        dr_pf_cs_middle_offset : float
            Radial offset for PF coil placement (m).
        z_tf_inside_half : float
            Half-height of the TF bore (m).
        dr_tf_inboard : float
            Thickness of the TF inboard leg (m).
        z_cs_coil_upper : float
            Upper z coordinate of the CS coil (m).

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Tuple of arrays containing the radial and vertical coordinates of PF coils in the group.
        """

        # Initialise as empty arrays; will be resized in the loop
        r_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))
        z_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))

        for coil in range(n_pf_coils_in_group[n_pf_group]):
            # Positions PF coil directly above centre of CS with offset from dr_pf_cs_middle_offset
            r_pf_coil_middle_group_array[n_pf_group, coil] = (
                r_cs_middle + dr_pf_cs_middle_offset
            )

            # Z coordinate of coil enforced so as not
            # to occupy the same space as the Central Solenoid
            # Set sign: +1 for coil 0, -1 for coil 1
            sign = 1.0 if coil == 0 else -1.0
            z_pf_coil_middle_group_array[n_pf_group, coil] = sign * (
                z_cs_coil_upper
                + 0.1e0
                + 0.5e0 * ((z_tf_inside_half - z_cs_coil_upper) + dr_tf_inboard + 0.1e0)
            )
        return r_pf_coil_middle_group_array, z_pf_coil_middle_group_array

    def place_pf_above_tf(
        self,
        n_pf_coils_in_group: np.ndarray,
        n_pf_group: int,
        rmajor: float,
        triang: float,
        rminor: float,
        itart: int,
        itartpf: int,
        z_tf_inside_half: float,
        dz_tf_upper_lower_midplane: float,
        z_tf_top: float,
        top_bottom: int,
        rpf2: float,
        zref: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, int]:
        """Calculates and places poloidal field (PF) coils above the toroidal field (TF) coils for a given group.

        Parameters
        ----------
        n_pf_coils_in_group : np.ndarray
            Array containing the number of PF coils in each group.
        n_pf_group : int
            Index of the PF coil group to process.
        rmajor : float
            Major radius of the device.
        triang : float
            Triangularity parameter for coil placement.
        rminor : float
            Minor radius of the device.
        itart : int
            Flag indicating ST configuration.
        itartpf : int
            Flag indicating PF coil configuration for ST.
        z_tf_inside_half : float
            Half-height of the TF coil inside region.
        dz_tf_upper_lower_midplane : float
            Height difference parameter for PF coil placement.
        z_tf_top : float
            Top z-coordinate of the TF coil.
        top_bottom : int
            Indicator for coil placement above (+1) or below (-1) the midplane.
        rpf2 : float
            Radial offset parameter for PF coil placement.
        zref : np.ndarray
            Array of reference z-coordinates for PF coil placement.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, int]
            Tuple containing arrays of radial and vertical positions of PF coil middles for the specified group,
            and the updated top_bottom indicator.
        """
        # Initialise as empty arrays; will be resized in the loop
        r_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))
        z_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))

        for coil in range(n_pf_coils_in_group[n_pf_group]):
            # Place PF coils at radius determined by rmajor, triang and rminor
            r_pf_coil_middle_group_array[n_pf_group, coil] = (
                rmajor + rpf2 * triang * rminor
            )

            # Set sign: +1 for coil 0, -1 for coil 1
            sign = 1.0 if coil == 0 else -1.0
            if itart == 1 and itartpf == 0:
                z_pf_coil_middle_group_array[n_pf_group, coil] = (
                    z_tf_inside_half - zref[n_pf_group]
                ) * sign
            else:
                if top_bottom == 1:  # this coil is above midplane
                    z_pf_coil_middle_group_array[n_pf_group, coil] = z_tf_top + 0.86e0
                    top_bottom = -1
                else:  # this coil is below midplane
                    z_pf_coil_middle_group_array[n_pf_group, coil] = -1.0e0 * (
                        z_tf_top - dz_tf_upper_lower_midplane + 0.86e0
                    )
                    top_bottom = 1

        return r_pf_coil_middle_group_array, z_pf_coil_middle_group_array, top_bottom

    def place_pf_outside_tf(
        self,
        n_pf_coils_in_group: np.ndarray,
        n_pf_group: int,
        rminor: float,
        zref: np.ndarray,
        i_tf_shape: int,
        i_r_pf_outside_tf_placement: int,
        r_pf_outside_tf_midplane: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculates the radial and vertical positions of poloidal field (PF) coils placed outside the toroidal field (TF) coil.

        Parameters
        ----------
        n_pf_coils_in_group : np.ndarray
            Array containing the number of PF coils in each group.
        n_pf_group : int
            Index of the PF coil group to process.
        rminor : float
            Minor radius of the device.
        zref : np.ndarray
            Reference vertical positions for each PF coil group.
        i_tf_shape : int
            Integer flag indicating TF coil shape (2 for picture frame, others for D-shape).
        i_r_pf_outside_tf_placement : int
            Placement switch for PF coil radius (1 for constant/stacked, 0 for following TF curve).
        r_pf_outside_tf_midplane : float
            Radial position of PF coil at the midplane.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Tuple containing arrays of radial and vertical positions of PF coil centers for the specified group.
        """

        # Initialise as empty arrays; will be resized in the loop
        r_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))
        z_pf_coil_middle_group_array = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))

        # PF coil is radially outside the TF coil

        for coil in range(n_pf_coils_in_group[n_pf_group]):
            sign = 1.0 if coil == 0 else -1.0

            z_pf_coil_middle_group_array[n_pf_group, coil] = (
                rminor * zref[n_pf_group] * sign
            )
            # Coil radius is constant / stacked for picture frame TF or if placement switch is set
            if i_tf_shape == 2 or i_r_pf_outside_tf_placement == 1:
                r_pf_coil_middle_group_array[n_pf_group, coil] = r_pf_outside_tf_midplane
            else:
                # Coil radius follows TF coil curve for TF (D-shape)
                r_pf_coil_middle_group_array[n_pf_group, coil] = math.sqrt(
                    r_pf_outside_tf_midplane**2
                    - z_pf_coil_middle_group_array[n_pf_group, coil] ** 2
                )
                try:
                    assert r_pf_coil_middle_group_array[n_pf_group, coil] < np.inf
                except AssertionError:
                    logger.exception(
                        "Element of pfcoil_variables.r_pf_coil_middle_group_array is inf. Kludging to 1e10."
                    )
                    r_pf_coil_middle_group_array[n_pf_group, coil] = 1e10
        return (
            r_pf_coil_middle_group_array,
            z_pf_coil_middle_group_array,
        )

    def place_pf_generally(
        self,
        n_pf_coils_in_group: np.ndarray,
        n_pf_group: int,
        rminor: float,
        rmajor: float,
        zref: np.ndarray,
        rref: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculates the radial and vertical positions of poloidal field (PF) coils placed in a general location.

        Parameters
        ----------
        n_pf_coils_in_group : numpy.ndarray
            Array containing the number of PF coils in each group.
        n_pf_group : int
            Index of the PF coil group to process.
        rminor : float
            Minor radius of the device.
        rmajor : float
            Major radius of the device.
        zref : numpy.ndarray
            Reference vertical positions for each PF coil group.
        rref : numpy.ndarray
            Reference radial positions for each PF coil group.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Tuple containing arrays of radial and vertical positions of PF coil centers for the specified group.
        """
        r_pf_coil_middle_group_array: np.ndarray = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))
        z_pf_coil_middle_group_array: np.ndarray = np.zeros((
            pfcoil_variables.n_pf_coil_groups,
            pfcoil_variables.N_PF_COILS_IN_GROUP_MAX,
        ))

        for coil in range(n_pf_coils_in_group[n_pf_group]):
            sign: float = 1.0 if coil == 0 else -1.0

            # Place as mutiples of minor radius from the midplane
            z_pf_coil_middle_group_array[n_pf_group, coil] = (
                rminor * zref[n_pf_group] * sign
            )
            # Place as multiples of minor radius from the plasma centre
            r_pf_coil_middle_group_array[n_pf_group, coil] = (
                rminor * rref[n_pf_group] + rmajor
            )
        return (
            r_pf_coil_middle_group_array,
            z_pf_coil_middle_group_array,
        )

    def efc(
        self,
        npts,
        rpts,
        zpts,
        brin,
        bzin,
        nfix,
        rfix,
        zfix,
        cfix,
        n_pf_coil_groups,
        n_pf_coils_in_group,
        r_pf_coil_middle_group_array,
        z_pf_coil_middle_group_array,
        alfa,
        bfix,
        gmat,
        bvec,
    ):
        """Calculates field coil currents.

        This routine calculates the currents required in a group
        of ring coils to produce a fixed field at prescribed
        locations. Additional ring coils with fixed currents are
        also allowed.

        Parameters
        ----------
        npts : int
            number of data points at which field is to be fixed; should
            be <= nptsmx
        rpts : np.ndarray
            coords of data points (m)
        zpts : np.ndarray
            coords of data points (m)
        brin : np.ndarray
            field components at data points (T)
        bzin : np.ndarray
            field components at data points (T)
        nfix : int
            number of coils with fixed currents, <= nfixmx
        rfix : np.ndarray
            coordinates of coils with fixed currents (m)
        zfix : np.ndarray
            coordinates of coils with fixed currents (m)
        cfix : np.ndarray
            Fixed currents (A)
        n_pf_coil_groups : int
            number of coil groups, where all coils in a group have the
            same current, <= n_pf_groups_max
        n_pf_coils_in_group : np.ndarray
            number of coils in each group, each value <= n_pf_coils_in_group_max
        r_pf_coil_middle_group_array : np.ndarray
            coords R(i,j), Z(i,j) of coil j in group i (m)
        z_pf_coil_middle_group_array : np.ndarray
            coords R(i,j), Z(i,j) of coil j in group i (m)
        alfa : float
            smoothing parameter (0 = no smoothing, 1.0D-9 = large
            smoothing)
        bfix : np.ndarray
            work array
        gmat : np.ndarray
            work array
        bvec : np.ndarray
            work array

        Returns
        -------
        tuple[float, np.ndarray]
            sum of squares of elements of residual vector, solution vector
            of coil currents in each group (A)
        """
        lrow1 = bfix.shape[0]
        lcol1 = gmat.shape[1]
        bfix = fixb(lrow1, npts, rpts, zpts, int(nfix), rfix, zfix, cfix)

        # Set up matrix equation
        nrws, gmat, bvec = mtrx(
            lrow1,
            lcol1,
            npts,
            rpts,
            zpts,
            brin,
            bzin,
            int(n_pf_coil_groups),
            n_pf_coils_in_group,
            r_pf_coil_middle_group_array,
            z_pf_coil_middle_group_array,
            alfa,
            bfix,
            int(pfcoil_variables.N_PF_COILS_IN_GROUP_MAX),
        )

        # Solve matrix equation
        ccls = self.solv(
            pfcoil_variables.N_PF_GROUPS_MAX, n_pf_coil_groups, nrws, gmat, bvec
        )

        # Calculate the norm of the residual vectors
        _brssq, _brnrm, _bzssq, _bznrm, ssq = rsid(
            npts, brin, bzin, nfix, int(n_pf_coil_groups), ccls, bfix, gmat
        )

        return ssq, ccls

    def tf_pf_collision_detector(self):
        #  Collision test between TF and PF coils for picture frame TF
        #  See issue 1612
        #  https://git.ccfe.ac.uk/process/process/-/issues/1612

        if tfv.i_tf_shape == 2:
            pf_tf_collision = 0

            for i in range(pfcoil_variables.n_pf_coil_groups):
                for ii in range(pfcoil_variables.n_pf_coil_groups):
                    for ij in range(pfcoil_variables.n_pf_coils_in_group[ii]):
                        if pfcoil_variables.r_pf_coil_middle_group_array[
                            ii, ij
                        ] <= (  # Outboard TF coil collision
                            pfcoil_variables.r_pf_outside_tf_midplane
                            - pfcoil_variables.dr_pf_tf_outboard_out_offset
                            + pfcoil_variables.r_pf_coil_middle[i]
                        ) and pfcoil_variables.r_pf_coil_middle_group_array[ii, ij] >= (
                            bv.r_tf_outboard_mid
                            - (0.5 * bv.dr_tf_outboard)
                            - pfcoil_variables.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1
                        if pfcoil_variables.r_pf_coil_middle_group_array[
                            ii, ij
                        ] <= (  # Inboard TF coil collision
                            bv.dr_bore
                            + bv.dr_cs
                            + bv.dr_cs_precomp
                            + bv.dr_cs_tf_gap
                            + bv.dr_tf_inboard
                            + pfcoil_variables.r_pf_coil_middle[i]
                        ) and pfcoil_variables.r_pf_coil_middle_group_array[ii, ij] >= (
                            bv.dr_bore
                            + bv.dr_cs
                            + bv.dr_cs_precomp
                            + bv.dr_cs_tf_gap
                            - pfcoil_variables.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1
                        if (  # Vertical TF coil collision
                            abs(pfcoil_variables.z_pf_coil_middle_group_array[ii, ij])
                            <= bv.z_tf_top + pfcoil_variables.r_pf_coil_middle[i]
                            and abs(
                                pfcoil_variables.z_pf_coil_middle_group_array[ii, ij]
                            )
                            >= bv.z_tf_top
                            - (0.5 * bv.dr_tf_outboard)
                            - pfcoil_variables.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1

                        if pf_tf_collision >= 1:
                            logger.error(
                                "One or more collision between TF and PF coils. Check PF placement."
                            )

    def solv(self, n_pf_groups_max, n_pf_coil_groups, nrws, gmat, bvec):
        """Solve a matrix using singular value decomposition.

        This routine solves the matrix equation for calculating the
        currents in a group of ring coils.

        Parameters
        ----------
        n_pf_groups_max : int
            maximum number of PF coil groups
        n_pf_coil_groups : int
            number of coil groups, where all coils in a group have the
            same current, <= n_pf_groups_max
        nrws : int
            actual number of rows to use
        gmat : numpy.ndarray
            work array
        bvec : numpy.ndarray
            work array

        Returns
        -------
        :
            solution vector of coil currents
            in each group (A) (ccls), rest are work arrays
        """
        ccls = np.zeros(n_pf_groups_max)
        work2 = np.zeros(n_pf_groups_max)

        umat, sigma, vmat = svd(gmat)

        for i in range(n_pf_coil_groups):
            work2[i] = 0.0e0
            for j in range(nrws):
                work2[i] = work2[i] + umat[j, i] * bvec[j]

        # Compute currents
        for i in range(n_pf_coil_groups):
            zvec = 0.0e0
            for j in range(n_pf_coil_groups):
                if sigma[j] > 1.0e-10:
                    zvec = work2[j] / sigma[j]

                ccls[i] = ccls[i] + vmat[j, i] * zvec

        return ccls

    def vsec(self):
        """Calculation of volt-second capability of PF system.


        This routine calculates the volt-second capability of the PF
        coil system.
        """
        if bv.iohcl == 0:
            # No Central Solenoid
            pfcoil_variables.nef = pfcoil_variables.n_pf_cs_plasma_circuits - 1
        else:
            pfcoil_variables.nef = pfcoil_variables.n_pf_cs_plasma_circuits - 2

        pfcoil_variables.vs_pf_coils_total_ramp = 0.0e0

        for i in range(pfcoil_variables.nef):
            pfcoil_variables.vsdum[i, 0] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1, i
                ]
                * pfcoil_variables.c_pf_coil_turn[i, 1]
            )
            pfcoil_variables.vsdum[i, 1] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1, i
                ]
                * pfcoil_variables.c_pf_coil_turn[i, 2]
            )
            pfcoil_variables.vs_pf_coils_total_ramp = (
                pfcoil_variables.vs_pf_coils_total_ramp
                + (pfcoil_variables.vsdum[i, 1] - pfcoil_variables.vsdum[i, 0])
            )

        # Central Solenoid startup volt-seconds
        if bv.iohcl != 0:
            pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 0] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1,
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2,
                ]
                * pfcoil_variables.c_pf_coil_turn[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2, 1
                ]
            )
            pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 1] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1,
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2,
                ]
                * pfcoil_variables.c_pf_coil_turn[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2, 2
                ]
            )
            pfcoil_variables.vs_cs_ramp = (
                pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 1]
                - pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 0]
            )

        # Total available volt-seconds for start-up
        pfcoil_variables.vs_cs_pf_total_ramp = (
            pfcoil_variables.vs_cs_ramp + pfcoil_variables.vs_pf_coils_total_ramp
        )

        # Burn volt-seconds
        if bv.iohcl != 0:
            pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 2] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1,
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2,
                ]
                * pfcoil_variables.c_pf_coil_turn[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 2, 4
                ]
            )
            pfcoil_variables.vs_cs_burn = (
                pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 2]
                - pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 1]
            )

        # PF volt-seconds during burn
        pfcoil_variables.vs_pf_coils_total_burn = 0.0e0
        for i in range(pfcoil_variables.nef):
            pfcoil_variables.vsdum[i, 2] = (
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1, i
                ]
                * pfcoil_variables.c_pf_coil_turn[i, 4]
            )
            pfcoil_variables.vs_pf_coils_total_burn = (
                pfcoil_variables.vs_pf_coils_total_burn
                + (pfcoil_variables.vsdum[i, 2] - pfcoil_variables.vsdum[i, 1])
            )

        pfcoil_variables.vs_cs_pf_total_burn = (
            pfcoil_variables.vs_cs_burn + pfcoil_variables.vs_pf_coils_total_burn
        )

        pfcoil_variables.vs_cs_pf_total_pulse = (
            pfcoil_variables.vs_cs_pf_total_ramp + pfcoil_variables.vs_cs_pf_total_burn
        )
        pfcoil_variables.vs_pf_coils_total_pulse = (
            pfcoil_variables.vs_pf_coils_total_ramp
            + pfcoil_variables.vs_pf_coils_total_burn
        )
        pfcoil_variables.vs_cs_total_pulse = (
            pfcoil_variables.vs_cs_burn + pfcoil_variables.vs_cs_ramp
        )

    def induct(self, output):
        """Calculates PF coil set mutual inductance matrix.


        This routine calculates the mutual inductances between all the
        PF coils.

        Parameters
        ----------
        output : bool
            switch for writing to output file
        """
        nohmax = 200
        nplas = 1

        _br = 0.0
        _bz = 0.0
        _psi = 0.0
        rc = np.zeros(pfcoil_variables.NGC2 + nohmax)
        zc = np.zeros(pfcoil_variables.NGC2 + nohmax)
        xc = np.zeros(pfcoil_variables.NGC2 + nohmax)
        cc = np.zeros(pfcoil_variables.NGC2 + nohmax)
        xcin = np.zeros(pfcoil_variables.NGC2 + nohmax)
        xcout = np.zeros(pfcoil_variables.NGC2 + nohmax)
        rplasma = np.zeros(nplas)
        zplasma = np.zeros(nplas)

        pfcoil_variables.ind_pf_cs_plasma_mutual[:, :] = 0.0

        # Break Central Solenoid into noh segments
        #
        # Choose noh so that the radial thickness of the coil is not thinner
        # than each segment is tall, i.e. the segments are pancake-like,
        # for the benefit of the mutual inductance calculations later

        noh = math.ceil(
            2.0e0
            * pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]
            / (
                pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]
                - pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1]
            )
        )

        if noh > nohmax:
            logger.error(
                "Max no. of segments noh for OH coil > nohmax; increase dr_cs lower bound"
                f"{noh=} {nohmax=} {bv.dr_cs=}"
            )

        noh = min(noh, nohmax)

        # TODO In FNSF case, noh = -7! noh should always be positive. Fortran
        # array allocation with -ve bound previously coerced to 0
        if noh < 0:
            noh = 0

        roh = np.zeros(noh)
        zoh = np.zeros(noh)

        if bv.iohcl != 0:
            roh[:] = pfcoil_variables.r_cs_middle

            delzoh = (
                2.0e0
                * pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]
                / noh
            )  # z_pf_coil_upper(n_cs_pf_coils) is the half-height of the coil
            for i in range(noh):
                zoh[i] = pfcoil_variables.z_pf_coil_upper[
                    pfcoil_variables.n_cs_pf_coils - 1
                ] - delzoh * (0.5e0 + i)

        rplasma[0] = pv.rmajor  # assumes nplas==1
        zplasma[0] = 0.0

        # Central Solenoid / plasma mutual inductance
        #
        # Improved calculation: Each Central Solenoid segment is now split into two filaments,
        # of radius reqv+deltar and reqv-deltar, respectively. The mutual inductance
        # of the segment with a plasma circuit is the mean of that calculated
        # using the two equivalent filaments.
        # Formulas and tables for the calculation of mutual and self-inductance
        # [Revised], Rosa and Grover, Scientific papers of the Bureau of Standards,
        # No. 169, 3rd ed., 1916. page 33

        for i in range(nplas):
            rc[i] = rplasma[i]
            zc[i] = zplasma[i]

        if bv.iohcl != 0:
            xohpl = 0.0
            if bv.dr_cs >= delzoh:
                deltar = math.sqrt((bv.dr_cs**2 - delzoh**2) / 12.0e0)
            else:
                # Set deltar to something small and +ve instead; allows solver
                # to continue and hopefully be constrained away from this point
                deltar = 1.0e-6

            for i in range(noh):
                rp = roh[i]
                zp = zoh[i]

                reqv = rp * (1.0e0 + delzoh**2 / (24.0e0 * rp**2))

                xcin, _br, _bz, _psi = calculate_b_field_at_point(
                    r_current_loop=rc,
                    z_current_loop=zc,
                    c_current_loop=cc,
                    r_test_point=reqv - deltar,
                    z_test_point=zp,
                )
                xcout, _br, _bz, _psi = calculate_b_field_at_point(
                    r_current_loop=rc,
                    z_current_loop=zc,
                    c_current_loop=cc,
                    r_test_point=reqv + deltar,
                    z_test_point=zp,
                )

                for ii in range(nplas):
                    xc[ii] = 0.5e0 * (xcin[ii] + xcout[ii])
                    xohpl = xohpl + xc[ii]

            pfcoil_variables.ind_pf_cs_plasma_mutual[
                pfcoil_variables.n_pf_cs_plasma_circuits - 1,
                pfcoil_variables.n_cs_pf_coils - 1,
            ] = (
                xohpl
                / (nplas * noh)
                * pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1]
            )
            pfcoil_variables.ind_pf_cs_plasma_mutual[
                pfcoil_variables.n_cs_pf_coils - 1,
                pfcoil_variables.n_pf_cs_plasma_circuits - 1,
            ] = pfcoil_variables.ind_pf_cs_plasma_mutual[
                pfcoil_variables.n_pf_cs_plasma_circuits - 1,
                pfcoil_variables.n_cs_pf_coils - 1,
            ]

        # Plasma self inductance
        pfcoil_variables.ind_pf_cs_plasma_mutual[
            pfcoil_variables.n_pf_cs_plasma_circuits - 1,
            pfcoil_variables.n_pf_cs_plasma_circuits - 1,
        ] = pv.ind_plasma

        # PF coil / plasma mutual inductances
        ncoils = 0

        for i in range(pfcoil_variables.n_pf_coil_groups):
            xpfpl = 0.0
            ncoils = ncoils + pfcoil_variables.n_pf_coils_in_group[i]
            rp = pfcoil_variables.r_pf_coil_middle[ncoils - 1]
            zp = pfcoil_variables.z_pf_coil_middle[ncoils - 1]
            xc, _br, _bz, _psi = calculate_b_field_at_point(
                r_current_loop=rc,
                z_current_loop=zc,
                c_current_loop=cc,
                r_test_point=rp,
                z_test_point=zp,
            )
            for ii in range(nplas):
                xpfpl = xpfpl + xc[ii]

            for j in range(pfcoil_variables.n_pf_coils_in_group[i]):
                ncoilj = ncoils + 1 - (j + 1)
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    ncoilj - 1, pfcoil_variables.n_pf_cs_plasma_circuits - 1
                ] = xpfpl / nplas * pfcoil_variables.n_pf_coil_turns[ncoilj - 1]
                pfcoil_variables.ind_pf_cs_plasma_mutual[
                    pfcoil_variables.n_pf_cs_plasma_circuits - 1, ncoilj - 1
                ] = pfcoil_variables.ind_pf_cs_plasma_mutual[
                    ncoilj - 1, pfcoil_variables.n_pf_cs_plasma_circuits - 1
                ]

        if bv.iohcl != 0:
            # Central Solenoid self inductance
            a = pfcoil_variables.r_cs_middle  # mean radius of coil
            b = (
                2.0e0
                * pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]
            )  # length of coil
            c = (
                pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]
                - pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1]
            )  # radial winding thickness
            pfcoil_variables.ind_pf_cs_plasma_mutual[
                pfcoil_variables.n_cs_pf_coils - 1, pfcoil_variables.n_cs_pf_coils - 1
            ] = self.selfinductance(
                a,
                b,
                c,
                pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1],
            )

            # Central Solenoid / PF coil mutual inductances
            for i in range(noh):
                rc[i] = roh[i]
                zc[i] = zoh[i]

            ncoils = 0
            for i in range(pfcoil_variables.n_pf_coil_groups):
                xohpf = 0.0
                ncoils = ncoils + pfcoil_variables.n_pf_coils_in_group[i]
                rp = pfcoil_variables.r_pf_coil_middle[ncoils - 1]
                zp = pfcoil_variables.z_pf_coil_middle[ncoils - 1]
                xc, _br, _bz, _psi = calculate_b_field_at_point(
                    r_current_loop=rc,
                    z_current_loop=zc,
                    c_current_loop=cc,
                    r_test_point=rp,
                    z_test_point=zp,
                )
                for ii in range(noh):
                    xohpf = xohpf + xc[ii]

                for j in range(pfcoil_variables.n_pf_coils_in_group[i]):
                    ncoilj = ncoils + 1 - (j + 1)
                    pfcoil_variables.ind_pf_cs_plasma_mutual[
                        ncoilj - 1, pfcoil_variables.n_cs_pf_coils - 1
                    ] = (
                        xohpf
                        * pfcoil_variables.n_pf_coil_turns[ncoilj - 1]
                        * pfcoil_variables.n_pf_coil_turns[
                            pfcoil_variables.n_cs_pf_coils - 1
                        ]
                        / noh
                    )
                    pfcoil_variables.ind_pf_cs_plasma_mutual[
                        pfcoil_variables.n_cs_pf_coils - 1, ncoilj - 1
                    ] = pfcoil_variables.ind_pf_cs_plasma_mutual[
                        ncoilj - 1, pfcoil_variables.n_cs_pf_coils - 1
                    ]

        # PF coil - PF coil inductances
        if bv.iohcl == 0:
            pfcoil_variables.nef = pfcoil_variables.n_cs_pf_coils
        else:
            pfcoil_variables.nef = pfcoil_variables.n_cs_pf_coils - 1

        for i in range(pfcoil_variables.nef):
            for j in range(pfcoil_variables.nef - 1):
                jj = j + 1 + 1 if j >= i else j + 1

                zc[j] = pfcoil_variables.z_pf_coil_middle[jj - 1]
                rc[j] = pfcoil_variables.r_pf_coil_middle[jj - 1]

            rp = pfcoil_variables.r_pf_coil_middle[i]
            zp = pfcoil_variables.z_pf_coil_middle[i]
            xc, _br, _bz, _psi = calculate_b_field_at_point(
                r_current_loop=rc,
                z_current_loop=zc,
                c_current_loop=cc,
                r_test_point=rp,
                z_test_point=zp,
            )
            for k in range(pfcoil_variables.nef):
                if k < i:
                    pfcoil_variables.ind_pf_cs_plasma_mutual[i, k] = (
                        xc[k]
                        * pfcoil_variables.n_pf_coil_turns[k]
                        * pfcoil_variables.n_pf_coil_turns[i]
                    )
                elif k == i:
                    rl = abs(
                        pfcoil_variables.z_pf_coil_upper[k]
                        - pfcoil_variables.z_pf_coil_lower[k]
                    ) / math.sqrt(np.pi)
                    pfcoil_variables.ind_pf_cs_plasma_mutual[k, k] = (
                        constants.RMU0
                        * pfcoil_variables.n_pf_coil_turns[k] ** 2
                        * pfcoil_variables.r_pf_coil_middle[k]
                        * (
                            math.log(8.0e0 * pfcoil_variables.r_pf_coil_middle[k] / rl)
                            - 1.75e0
                        )
                    )
                else:
                    pfcoil_variables.ind_pf_cs_plasma_mutual[i, k] = (
                        xc[k - 1]
                        * pfcoil_variables.n_pf_coil_turns[k]
                        * pfcoil_variables.n_pf_coil_turns[i]
                    )

        # Output section
        if not output:
            return

        op.oheadr(self.outfile, "PF Coil Inductances")
        op.ocmmnt(self.outfile, "Inductance matrix [H]:")
        op.oblnkl(self.outfile)

        with np.printoptions(precision=1):
            for ig in range(pfcoil_variables.nef):
                op.write(
                    self.outfile,
                    f"{ig}\t{pfcoil_variables.ind_pf_cs_plasma_mutual[: pfcoil_variables.n_pf_cs_plasma_circuits, ig]}",
                )

            if bv.iohcl != 0:
                op.write(
                    self.outfile,
                    f"CS\t\t\t{pfcoil_variables.ind_pf_cs_plasma_mutual[: pfcoil_variables.n_pf_cs_plasma_circuits, pfcoil_variables.n_pf_cs_plasma_circuits - 2]}",
                )

            op.write(
                self.outfile,
                f"Plasma\t{pfcoil_variables.ind_pf_cs_plasma_mutual[: pfcoil_variables.n_pf_cs_plasma_circuits, pfcoil_variables.n_pf_cs_plasma_circuits - 1]}",
            )

    def outpf(self):
        """Routine to write output from PF coil module to file.


        This routine writes the PF coil information to the output file.
        """
        op.oheadr(self.outfile, "Central Solenoid and PF Coils")

        if bv.iohcl == 0:
            op.ocmmnt(self.outfile, "No central solenoid included")
            op.oblnkl(self.outfile)
            op.ovarin(self.mfile, "Existence_of_central_solenoid", "(iohcl)", bv.iohcl)
        else:
            if pfcoil_variables.i_pf_conductor == 0:
                op.ocmmnt(self.outfile, "Superconducting central solenoid")

                op.ovarin(
                    self.outfile,
                    "Central solenoid superconductor material",
                    "(i_cs_superconductor)",
                    pfcoil_variables.i_cs_superconductor,
                )

                if pfcoil_variables.i_cs_superconductor == 1:
                    op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
                elif pfcoil_variables.i_cs_superconductor == 2:
                    op.ocmmnt(
                        self.outfile, "  (Bi-2212 high temperature superconductor)"
                    )
                elif pfcoil_variables.i_cs_superconductor == 3:
                    op.ocmmnt(self.outfile, "  (NbTi)")
                elif pfcoil_variables.i_cs_superconductor == 4:
                    op.ocmmnt(
                        self.outfile,
                        "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                    )
                elif pfcoil_variables.i_cs_superconductor == 5:
                    op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
                elif pfcoil_variables.i_cs_superconductor == 6:
                    op.ocmmnt(self.outfile, " (REBCO HTS)")
                elif pfcoil_variables.i_cs_superconductor == 7:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                    )
                elif pfcoil_variables.i_cs_superconductor == 8:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for REBCO)",
                    )
                elif pfcoil_variables.i_cs_superconductor == 9:
                    op.ocmmnt(
                        self.outfile,
                        " (Hazelton experimental data + Zhai conceptual model for REBCO)",
                    )

                op.osubhd(self.outfile, "Central Solenoid Current Density Limits :")
                op.ovarre(
                    self.outfile,
                    "Maximum field at Beginning Of Pulse (T)",
                    "(b_cs_peak_pulse_start)",
                    pfcoil_variables.b_cs_peak_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at BOP (A/m2)",
                    "(j_cs_conductor_critical_pulse_start)",
                    pfcoil_variables.j_cs_conductor_critical_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical cable current density at BOP (A/m2)",
                    "(jcableoh_bop)",
                    pfcoil_variables.jcableoh_bop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at BOP (A/m2)",
                    "(j_cs_critical_pulse_start)",
                    pfcoil_variables.j_cs_critical_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at BOP (A/m2)",
                    "(j_cs_pulse_start)",
                    pfcoil_variables.j_cs_pulse_start,
                    "OP ",
                )
                op.oblnkl(self.outfile)
                op.ovarre(
                    self.outfile,
                    "Maximum field at End Of Flattop (T)",
                    "(b_cs_peak_flat_top_end)",
                    pfcoil_variables.b_cs_peak_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at EOF (A/m2)",
                    "(j_cs_conductor_critical_flat_top_end)",
                    pfcoil_variables.j_cs_conductor_critical_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical cable current density at EOF (A/m2)",
                    "(jcableoh_eof)",
                    pfcoil_variables.jcableoh_eof,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at EOF (A/m2)",
                    "(j_cs_critical_flat_top_end)",
                    pfcoil_variables.j_cs_critical_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at EOF (A/m2)",
                    "(j_cs_flat_top_end)",
                    pfcoil_variables.j_cs_flat_top_end,
                )
                for i in range(len(pfcoil_variables.r_pf_cs_current_filaments)):
                    op.ovarre(
                        self.mfile,
                        f"Radial position of CS filament {i}",
                        f"r_pf_cs_current_filaments{i}",
                        pfcoil_variables.r_pf_cs_current_filaments[i],
                    )
                for i in range(len(pfcoil_variables.z_pf_cs_current_filaments)):
                    op.ovarre(
                        self.mfile,
                        f"Vertical position of CS filament {i}",
                        f"z_pf_cs_current_filaments{i}",
                        pfcoil_variables.z_pf_cs_current_filaments[i],
                    )
                op.oblnkl(self.outfile)
                # MDK add bv.dr_cs, bv.dr_bore and bv.dr_cs_tf_gap as they can be iteration variables
                op.ovarre(self.outfile, "CS inside radius (m)", "(dr_bore)", bv.dr_bore)
                op.ovarre(self.outfile, "CS thickness (m)", "(dr_cs)", bv.dr_cs)
                op.ovarre(
                    self.outfile,
                    "Gap between central solenoid and TF coil (m)",
                    "(dr_cs_tf_gap)",
                    bv.dr_cs_tf_gap,
                )
                op.ovarre(
                    self.outfile,
                    "CS overall cross-sectional area (m2)",
                    "(a_cs_poloidal)",
                    pfcoil_variables.a_cs_poloidal,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS radial middle (m)",
                    "(r_cs_middle)",
                    pfcoil_variables.r_cs_middle,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS conductor+void cross-sectional area (m2)",
                    "(awpoh)",
                    pfcoil_variables.awpoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "   CS conductor cross-sectional area (m2)",
                    "(awpoh*(1-f_a_cs_void))",
                    pfcoil_variables.awpoh * (1.0e0 - pfcoil_variables.f_a_cs_void),
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "   CS void cross-sectional area (m2)",
                    "(awpoh*f_a_cs_void)",
                    pfcoil_variables.awpoh * pfcoil_variables.f_a_cs_void,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel cross-sectional area (m2)",
                    "(a_cs_poloidal-awpoh)",
                    pfcoil_variables.a_cs_poloidal - pfcoil_variables.awpoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel area fraction",
                    "(f_a_cs_turn_steel)",
                    pfcoil_variables.f_a_cs_turn_steel,
                )
                if pfcoil_variables.i_cs_stress == 1:
                    op.ocmmnt(self.outfile, "Hoop + axial stress considered")
                else:
                    op.ocmmnt(self.outfile, "Only hoop stress considered")

                op.ovarin(
                    self.outfile,
                    "Switch for CS stress calculation",
                    "(i_cs_stress)",
                    pfcoil_variables.i_cs_stress,
                )
                op.ovarre(
                    self.outfile,
                    "Allowable stress in CS steel (Pa)",
                    "(alstroh)",
                    pfcoil_variables.alstroh,
                )
                op.ovarre(
                    self.outfile,
                    "Hoop stress in CS steel (Pa)",
                    "(sig_hoop)",
                    pfcoil_variables.sig_hoop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Axial stress in CS steel (Pa)",
                    "(stress_z_cs_self_peak_midplane)",
                    pfcoil_variables.stress_z_cs_self_peak_midplane,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Maximum shear stress in CS steel for the Tresca criterion (Pa)",
                    "(s_shear_cs_peak)",
                    pfcoil_variables.s_shear_cs_peak,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Axial force in CS (N)",
                    "(forc_z_cs_self_peak_midplane)",
                    pfcoil_variables.forc_z_cs_self_peak_midplane,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Residual manufacturing strain in CS superconductor material",
                    "(str_cs_con_res)",
                    tfv.str_cs_con_res,
                )
                op.ovarre(
                    self.outfile,
                    "Copper fraction in strand",
                    "(fcuohsu)",
                    pfcoil_variables.fcuohsu,
                )
                # If REBCO material is used, print copperaoh_m2
                if (
                    pfcoil_variables.i_cs_superconductor == 6
                    or pfcoil_variables.i_cs_superconductor == 8
                    or pfcoil_variables.i_cs_superconductor == 9
                ):
                    op.ovarre(
                        self.outfile,
                        "CS current/copper area (A/m2)",
                        "(copperaoh_m2)",
                        rcv.copperaoh_m2,
                    )
                    op.ovarre(
                        self.outfile,
                        "Max CS current/copper area (A/m2)",
                        "(copperaoh_m2_max)",
                        rcv.copperaoh_m2_max,
                    )

                op.ovarre(
                    self.outfile,
                    "Void (coolant) fraction in conductor",
                    "(f_a_cs_void)",
                    pfcoil_variables.f_a_cs_void,
                )
                op.ovarre(
                    self.outfile,
                    "Helium coolant temperature (K)",
                    "(tftmp)",
                    tfv.tftmp,
                )
                op.ovarre(
                    self.outfile,
                    "CS temperature margin (K)",
                    "(temp_cs_superconductor_margin)",
                    pfcoil_variables.temp_cs_superconductor_margin,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Minimum permitted temperature margin (K)",
                    "(temp_cs_superconductor_margin_min)",
                    tfv.temp_cs_superconductor_margin_min,
                )
                # only output CS fatigue model for pulsed reactor design
                if pv.f_c_plasma_inductive > 0.0e-4:
                    op.ovarre(
                        self.outfile,
                        "Residual hoop stress in CS Steel (Pa)",
                        "(residual_sig_hoop)",
                        csfv.residual_sig_hoop,
                    )
                    op.ovarre(
                        self.outfile,
                        "Minimum burn time (s)",
                        "(t_burn_min)",
                        ctv.t_burn_min,
                    )
                    op.ovarre(
                        self.outfile,
                        "Initial vertical crack size (m)",
                        "(t_crack_vertical)",
                        csfv.t_crack_vertical,
                    )
                    op.ovarre(
                        self.outfile,
                        "Initial radial crack size (m)",
                        "(t_crack_radial)",
                        csfv.t_crack_radial,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn area (m)",
                        "(a_cs_turn)",
                        pfcoil_variables.a_cs_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn length (m)",
                        "(dr_cs_turn)",
                        pfcoil_variables.dr_cs_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn internal cable space radius (m)",
                        "(radius_cs_turn_cable_space)",
                        pfcoil_variables.radius_cs_turn_cable_space,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn width (m)",
                        "(dz_cs_turn)",
                        pfcoil_variables.dz_cs_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS structural vertical thickness (m)",
                        "(dz_cs_turn_conduit)",
                        csfv.dz_cs_turn_conduit,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS structural radial thickness (m)",
                        "(dr_cs_turn_conduit)",
                        csfv.dr_cs_turn_conduit,
                    )
                    op.ovarre(
                        self.outfile,
                        "Allowable number of cycles till CS fracture",
                        "(n_cycle)",
                        csfv.n_cycle,
                        "OP ",
                    )
                    op.ovarre(
                        self.outfile,
                        "Minimum number of cycles required till CS fracture",
                        "(n_cycle_min)",
                        csfv.n_cycle_min,
                        "OP ",
                    )
                # Check whether CS coil is hitting any limits
                if (
                    abs(pfcoil_variables.j_cs_flat_top_end)
                    > 0.99e0
                    * abs(ctv.fjohc * pfcoil_variables.j_cs_critical_flat_top_end)
                ) or (
                    abs(pfcoil_variables.j_cs_pulse_start)
                    > 0.99e0
                    * abs(ctv.fjohc0 * pfcoil_variables.j_cs_critical_pulse_start)
                ):
                    pfcoil_variables.cslimit = True

                if (
                    pfcoil_variables.j_cs_flat_top_end
                    / pfcoil_variables.j_cs_critical_flat_top_end
                    > 0.7
                ):
                    logger.error(
                        "j_cs_flat_top_end / j_cs_critical_flat_top_end shouldn't be above 0.7 "
                        "for engineering reliability"
                    )

                if (
                    pfcoil_variables.j_cs_pulse_start
                    / pfcoil_variables.j_cs_critical_pulse_start
                    > 0.7
                ):
                    logger.error(
                        "j_cs_pulse_start / j_cs_critical_pulse_start shouldn't be above 0.7 "
                        "for engineering reliability"
                    )

                if (
                    pfcoil_variables.temp_cs_superconductor_margin
                    < 1.01e0 * tfv.temp_cs_superconductor_margin_min
                ):
                    pfcoil_variables.cslimit = True
                if not pfcoil_variables.cslimit:
                    logger.warning(
                        "CS not using max current density: further optimisation may be possible"
                    )

                # REBCO fractures in strains above ~+/- 0.7%
                if (
                    pfcoil_variables.i_cs_superconductor == 6
                    or pfcoil_variables.i_cs_superconductor == 8
                    or pfcoil_variables.i_cs_superconductor == 9
                ) and abs(tfv.str_cs_con_res) > 0.7e-2:
                    logger.error(
                        "Non physical strain used in CS. Use superconductor strain < +/- 0.7%"
                    )

                if (
                    pfcoil_variables.i_pf_superconductor == 6
                    or pfcoil_variables.i_pf_superconductor == 8
                    or pfcoil_variables.i_pf_superconductor == 9
                ) and abs(tfv.str_pf_con_res) > 0.7e-2:
                    logger.error(
                        "Non physical strain used in PF. Use superconductor strain < +/- 0.7%"
                    )

            else:
                op.ocmmnt(self.outfile, "Resistive central solenoid")

        if pfcoil_variables.i_pf_conductor == 0:
            op.oblnkl(self.outfile)
            op.ocmmnt(self.outfile, "Superconducting PF coils")

            op.ovarin(
                self.outfile,
                "PF coil superconductor material",
                "(i_pf_superconductor)",
                pfcoil_variables.i_pf_superconductor,
            )

            if pfcoil_variables.i_pf_superconductor == 1:
                op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
            elif pfcoil_variables.i_pf_superconductor == 2:
                op.ocmmnt(self.outfile, "  (Bi-2212 high temperature superconductor)")
            elif pfcoil_variables.i_pf_superconductor == 3:
                op.ocmmnt(self.outfile, "  (NbTi)")
            elif pfcoil_variables.i_pf_superconductor == 4:
                op.ocmmnt(
                    self.outfile,
                    "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                )
            elif pfcoil_variables.i_pf_superconductor == 5:
                op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
            elif pfcoil_variables.i_pf_superconductor == 6:
                op.ocmmnt(
                    self.outfile,
                    " (REBCO 2nd generation HTS superconductor in CrCo strand)",
                )
            elif pfcoil_variables.i_pf_superconductor == 7:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                )
            elif pfcoil_variables.i_pf_superconductor == 8:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for REBCO)",
                )
            elif pfcoil_variables.i_pf_superconductor == 9:
                op.ocmmnt(
                    self.outfile,
                    " (Hazelton experimental data + Zhai conceptual model for REBCO)",
                )

            op.ovarre(
                self.outfile,
                "Copper fraction in conductor",
                "(fcupfsu)",
                pfcoil_variables.fcupfsu,
            )

            op.osubhd(self.outfile, "PF Coil Case Stress :")
            op.ovarre(
                self.outfile,
                "Maximum permissible tensile stress (MPa)",
                "(sigpfcalw)",
                pfcoil_variables.sigpfcalw,
            )
            op.ovarre(
                self.outfile,
                "JxB hoop force fraction supported by case",
                "(sigpfcf)",
                pfcoil_variables.sigpfcf,
            )

        else:
            op.oblnkl(self.outfile)
            op.ocmmnt(self.outfile, "Resistive PF coils")

            op.osubhd(self.outfile, "Resistive Power :")
            op.ovarre(
                self.outfile,
                "PF coil resistive power (W)",
                "(p_pf_coil_resistive_total_flat_top)",
                pfcoil_variables.p_pf_coil_resistive_total_flat_top,
                "OP ",
            )
            if bv.iohcl != 0:
                op.ovarre(
                    self.outfile,
                    "Central solenoid resistive power (W)",
                    "(p_cs_resistive_flat_top)",
                    pfcoil_variables.p_cs_resistive_flat_top,
                    "OP ",
                )

        # pfcoil_variables.nef is the number of coils excluding the Central Solenoid
        pfcoil_variables.nef = pfcoil_variables.n_cs_pf_coils
        if bv.iohcl != 0:
            pfcoil_variables.nef = pfcoil_variables.nef - 1

        op.osubhd(self.outfile, "Geometry of PF coils, central solenoid and plasma:")
        op.write(
            self.outfile,
            "coil\t\t\tR(m)\t\tZ(m)\t\tdR(m)\t\tdZ(m)\t\tturns",
        )
        op.oblnkl(self.outfile)

        # PF coils
        for k in range(pfcoil_variables.nef):
            op.write(
                self.outfile,
                f"PF {k}\t\t\t{pfcoil_variables.r_pf_coil_middle[k]:.2e}\t{pfcoil_variables.z_pf_coil_middle[k]:.2e}\t{pfcoil_variables.r_pf_coil_outer[k] - pfcoil_variables.r_pf_coil_inner[k]:.2e}\t{abs(pfcoil_variables.z_pf_coil_upper[k] - pfcoil_variables.z_pf_coil_lower[k]):.2e}\t{pfcoil_variables.n_pf_coil_turns[k]:.2e}",
            )

        for k in range(pfcoil_variables.nef):
            op.ovarre(
                self.mfile,
                f"PF coil {k} radius (m)",
                f"(r_pf_coil_middle[{k}])",
                pfcoil_variables.r_pf_coil_middle[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical position (m)",
                f"(z_pf_coil_middle[{k}])",
                pfcoil_variables.z_pf_coil_middle[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} radial thickness (m)",
                f"(pfdr({k}))",
                pfcoil_variables.r_pf_coil_outer[k]
                - pfcoil_variables.r_pf_coil_inner[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical thickness (m)",
                f"(pfdz({k}))",
                pfcoil_variables.z_pf_coil_upper[k]
                - pfcoil_variables.z_pf_coil_lower[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} turns",
                f"(n_pf_coil_turns[{k}])",
                pfcoil_variables.n_pf_coil_turns[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} current (MA)",
                f"(c_pf_cs_coils_peak_ma[{k}])",
                pfcoil_variables.c_pf_cs_coils_peak_ma[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} field (T)",
                f"(b_pf_coil_peak[{k}])",
                pfcoil_variables.b_pf_coil_peak[k],
            )
        self.tf_pf_collision_detector()

        # Central Solenoid, if present
        if bv.iohcl != 0:
            op.write(
                self.outfile,
                f"CS\t\t\t\t{pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.z_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1] - pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{abs(pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1] - pfcoil_variables.z_pf_coil_lower[pfcoil_variables.n_cs_pf_coils - 1]):.2e}\t{pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.pfcaseth[pfcoil_variables.n_cs_pf_coils - 1]:.2e}",
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radius (m)",
                "(r_pf_coil_middle[n_cs_pf_coils-1])",
                pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical position (m)",
                "(z_pf_coil_middle[n_cs_pf_coils-1])",
                pfcoil_variables.z_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radial thickness (m)",
                "(ohdr)",
                (
                    pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]
                    - pfcoil_variables.r_pf_coil_inner[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                ),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical thickness (m)",
                "(dz_cs_full)",
                (pfcoil_variables.dz_cs_full),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid full radial width (m)",
                "(dr_cs_full)",
                (pfcoil_variables.dr_cs_full),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid turns",
                "(n_pf_coil_turns[n_cs_pf_coils-1])",
                pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid current (MA)",
                "(c_pf_cs_coils_peak_ma[n_cs_pf_coils-1])",
                pfcoil_variables.c_pf_cs_coils_peak_ma[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid field (T)",
                "(b_pf_coil_peak[n_cs_pf_coils-1])",
                pfcoil_variables.b_pf_coil_peak[pfcoil_variables.n_cs_pf_coils - 1],
            )

        # Plasma
        op.write(
            self.outfile,
            f"Plasma\t\t\t{pv.rmajor:.2e}\t0.0e0\t\t{2.0e0 * pv.rminor:.2e}\t{2.0e0 * pv.rminor * pv.kappa:.2e}\t1.0e0",
        )

        op.osubhd(self.outfile, "PF Coil Information at Peak Current:")

        op.write(
            self.outfile,
            "coil\tcurrent\t\tallowed J\tactual J\tJ\t\tcond. mass\tsteel mass\tfield",
        )
        op.write(
            self.outfile,
            "\t(MA)\t\t(A/m2)\t\t(A/m2)\t\f_temp_plasma_ion_electron\t\t(kg)\t\t(kg)\t\t(T)",
        )

        op.oblnkl(self.outfile)

        # PF coils
        for k in range(pfcoil_variables.nef):
            if pfcoil_variables.i_pf_conductor == 0:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfcoil_variables.c_pf_cs_coils_peak_ma[k]:.2e}\t{pfcoil_variables.j_pf_wp_critical[k]:.2e}\t{pfcoil_variables.j_pf_coil_wp_peak[k]:.2e}\t{pfcoil_variables.j_pf_coil_wp_peak[k] / pfcoil_variables.j_pf_wp_critical[k]:.2e}\t{pfcoil_variables.m_pf_coil_conductor[k]:.2e}\t{pfcoil_variables.m_pf_coil_structure[k]:.2e}\t{pfcoil_variables.b_pf_coil_peak[k]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfcoil_variables.c_pf_cs_coils_peak_ma[k]:.2e}\t-1.0e0\t{pfcoil_variables.j_pf_coil_wp_peak[k]:.2e}\t1.0e0\t{pfcoil_variables.m_pf_coil_conductor[k]:.2e}\t{pfcoil_variables.m_pf_coil_structure[k]:.2e}\t{pfcoil_variables.b_pf_coil_peak[k]:.2e}\t",
                )

        # Central Solenoid, if present
        if bv.iohcl != 0:
            if pfcoil_variables.i_pf_conductor == 0:
                # Issue #328
                op.write(
                    self.outfile,
                    f"CS\t\t{pfcoil_variables.c_pf_cs_coils_peak_ma[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.j_pf_wp_critical[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{max(abs(pfcoil_variables.j_cs_pulse_start), abs(pfcoil_variables.j_cs_flat_top_end)):.2e}\t{max(abs(pfcoil_variables.j_cs_pulse_start), abs(pfcoil_variables.j_cs_flat_top_end)) / pfcoil_variables.j_pf_wp_critical[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.m_pf_coil_conductor[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.m_pf_coil_structure[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.b_pf_coil_peak[pfcoil_variables.n_cs_pf_coils - 1]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"CS\t\t{pfcoil_variables.c_pf_cs_coils_peak_ma[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t-1.0e0\t{max(abs(pfcoil_variables.j_cs_pulse_start)):.2e}\t{abs(pfcoil_variables.j_cs_flat_top_end):.2e}\t1.0e0\t{pfcoil_variables.m_pf_coil_conductor[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.m_pf_coil_structure[pfcoil_variables.n_cs_pf_coils - 1]:.2e}\t{pfcoil_variables.b_pf_coil_peak[pfcoil_variables.n_cs_pf_coils - 1]:.2e}",
                )

        # Miscellaneous totals
        op.write(
            self.outfile,
            "\t" * 1 + "------" + "\t" * 8 + "---------" + "\t" + "---------",
        )

        op.write(
            self.outfile,
            "\t" * 1
            + f"{pfcoil_variables.ricpf:.2e}"
            + "\t" * 7
            + f"{pfcoil_variables.m_pf_coil_conductor_total:.2e}\t{pfcoil_variables.m_pf_coil_structure_total:.2e}",
        )

        op.osubhd(self.outfile, "PF coil current scaling information :")
        op.ovarre(
            self.outfile,
            "Sum of squares of residuals ",
            "(ssq0)",
            pfcoil_variables.ssq0,
            "OP ",
        )
        op.ovarre(
            self.outfile, "Smoothing parameter ", "(alfapf)", pfcoil_variables.alfapf
        )

    def outvolt(self):
        """Writes volt-second information to output file.

        This routine writes the PF coil volt-second data to the
        output file.
        """
        op.oheadr(self.outfile, "Volt Second Consumption")

        op.write(self.outfile, "\t" * 3 + "volt-sec\t\t\tvolt-sec\t\tvolt-sec")
        op.write(self.outfile, "\t" * 3 + "start-up\t\t\t_burn\t\t\ttotal")
        op.write(
            self.outfile,
            f"PF coils:\t\t{pfcoil_variables.vs_pf_coils_total_ramp:.2f}\t\t\t\t{pfcoil_variables.vs_pf_coils_total_burn:.2f}\t\t\t{pfcoil_variables.vs_pf_coils_total_pulse:.2f}",
        )
        op.write(
            self.outfile,
            f"CS coil:\t\t{pfcoil_variables.vs_cs_ramp:.2f}\t\t\t\t{pfcoil_variables.vs_cs_burn:.2f}\t\t\t{pfcoil_variables.vs_cs_total_pulse:.2f}",
        )
        op.write(
            self.outfile, "\t" * 3 + "-" * 7 + "\t" * 4 + "-" * 7 + "\t" * 3 + "-" * 7
        )
        op.write(
            self.outfile,
            f"Total:\t\t\t{pfcoil_variables.vs_cs_pf_total_ramp:.2f}\t\t\t\t{pfcoil_variables.vs_cs_pf_total_burn:.2f}\t\t\t{pfcoil_variables.vs_cs_pf_total_pulse:.2f}",
        )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Total volt-second consumption by coils (Wb)",
            "(vs_cs_pf_total_pulse)",
            pfcoil_variables.vs_cs_pf_total_pulse,
            "OP",
        )
        op.ovarre(
            self.outfile,
            "Total volt-second available for burn phase (Wb)",
            "(vs_cs_pf_total_burn)",
            pfcoil_variables.vs_cs_pf_total_burn,
            "OP",
        )

        op.osubhd(self.outfile, "Summary of volt-second consumption by circuit (Wb):")

        op.write(self.outfile, "circuit\t\t\tBOP\t\t\tBOF\t\tEOF")
        op.oblnkl(self.outfile)

        for k in range(pfcoil_variables.nef):
            op.write(
                self.outfile,
                f"\t{k}\t\t\t{pfcoil_variables.vsdum[k, 0]:.3f}\t\t\t{pfcoil_variables.vsdum[k, 1]:.3f}\t\t{pfcoil_variables.vsdum[k, 2]:.3f}",
            )

        op.write(
            self.outfile,
            f"\tCS coil\t\t\t{pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 0]:.3f}\t\t\t{pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 1]:.3f}\t\t{pfcoil_variables.vsdum[pfcoil_variables.n_cs_pf_coils - 1, 2]:.3f}",
        )

        op.oshead(self.outfile, "Waveforms")
        op.ocmmnt(self.outfile, "Currents (Amps/coil) as a function of time:")
        op.oblnkl(self.outfile)

        op.write(self.outfile, "\t" * 8 + "time (sec)")
        line = "\t\t"
        for k in range(6):
            line += f"\t\t{tv.t_pulse_cumulative[k]:.2f}"
        op.write(self.outfile, line)

        line = "\t\t"
        for k in range(6):
            label = tv.timelabel[k]
            line += f"\t\t{label}"
        op.write(self.outfile, line)

        op.ocmmnt(self.outfile, "circuit")

        for k in range(pfcoil_variables.n_pf_cs_plasma_circuits - 1):
            line = f"\t{k}\t\t"
            for jj in range(6):
                line += f"\t{pfcoil_variables.c_pf_coil_turn[k, jj] * pfcoil_variables.n_pf_coil_turns[k]:.3e}"
            op.write(self.outfile, line)

        line = "Plasma (A)\t\t"
        for jj in range(6):
            line += f"\t{pfcoil_variables.c_pf_coil_turn[pfcoil_variables.n_pf_cs_plasma_circuits - 1, jj]:.3e}"

        op.write(self.outfile, line)

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "This consists of: CS coil field balancing:")
        for k in range(pfcoil_variables.n_pf_cs_plasma_circuits - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{pfcoil_variables.c_pf_coil_turn[k, 0] * pfcoil_variables.n_pf_coil_turns[k]:.3e}\t"
                    f"{pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.n_pf_coil_turns[k]:.3e}\t"
                    f"{-pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.n_pf_coil_turns[k] * (pfcoil_variables.f_j_cs_start_end_flat_top / pfcoil_variables.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{-pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.n_pf_coil_turns[k] * (pfcoil_variables.f_j_cs_start_end_flat_top / pfcoil_variables.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{-pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.n_pf_coil_turns[k] * (1.0e0 / pfcoil_variables.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{pfcoil_variables.c_pf_coil_turn[k, 5] * pfcoil_variables.n_pf_coil_turns[k]:.3e}"
                ),
            )

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "And: equilibrium field:")
        for k in range(pfcoil_variables.n_pf_cs_plasma_circuits - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{0.0:.3e}\t{0.0:.3e}\t"
                    f"{(pfcoil_variables.c_pf_coil_turn[k, 2] + pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.f_j_cs_start_end_flat_top / pfcoil_variables.f_j_cs_start_pulse_end_flat_top) * pfcoil_variables.n_pf_coil_turns[k]:.3e}\t"
                    f"{(pfcoil_variables.c_pf_coil_turn[k, 3] + pfcoil_variables.c_pf_coil_turn[k, 1] * pfcoil_variables.f_j_cs_start_end_flat_top / pfcoil_variables.f_j_cs_start_pulse_end_flat_top) * pfcoil_variables.n_pf_coil_turns[k]:.3e}\t"
                    f"{(pfcoil_variables.c_pf_coil_turn[k, 4] + pfcoil_variables.c_pf_coil_turn[k, 1] * 1.0e0 / pfcoil_variables.f_j_cs_start_pulse_end_flat_top) * pfcoil_variables.n_pf_coil_turns[k]:.3e}\t"
                    "0.0e0"
                ),
            )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Pulse / end of flat-top",
            "(f_j_cs_start_pulse_end_flat_top)",
            pfcoil_variables.f_j_cs_start_pulse_end_flat_top,
        )
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Flat-top / end of flat-top",
            "(f_j_cs_start_end_flat_top)",
            pfcoil_variables.f_j_cs_start_end_flat_top,
            "OP ",
        )

        op.oshead(self.outfile, "PF Circuit Waveform Data")
        op.ovarin(
            self.outfile,
            "Number of PF circuits including CS and plasma",
            "(n_pf_cs_plasma_circuits)",
            pfcoil_variables.n_pf_cs_plasma_circuits,
        )
        for k in range(pfcoil_variables.n_pf_cs_plasma_circuits):
            for jjj in range(6):
                if k == pfcoil_variables.n_pf_cs_plasma_circuits - 1:
                    circuit_name = f"Plasma Time point {jjj} (A)"
                    circuit_var_name = f"(plasmat{jjj})"
                elif k == pfcoil_variables.n_pf_cs_plasma_circuits - 2:
                    circuit_name = f"CS Circuit Time point {jjj} (A)"
                    circuit_var_name = f"(cs t{jjj})"
                else:
                    circuit_name = f"PF Circuit {k} Time point {jjj} (A)"
                    circuit_var_name = f"(pfc{k}t{jjj})"

                op.ovarre(
                    self.outfile,
                    circuit_name,
                    circuit_var_name,
                    pfcoil_variables.c_pf_coil_turn[k, jjj]
                    * pfcoil_variables.n_pf_coil_turns[k],
                )

    def selfinductance(self, a, b, c, n):
        """Calculates the selfinductance using Bunet's formula.


        This routine calculates the self inductance in Henries
        Radiotron Designers Handbook (4th Edition) chapter 10

        Parameters
        ----------
        a : float
            mean radius of coil (m)
        b : float
            length of coil (m) (given as l in the reference)
        c : float
            radial winding thickness (m)
        n : float
            number of turns


        Returns
        -------
        :
            the self inductance in Henries
        """
        return (
            (1.0e-6 / 0.0254e0)
            * a**2
            * n**2
            / (9.0e0 * a + 10.0e0 * b + 8.4e0 * c + 3.2e0 * c * b / a)
        )

    def waveform(self):
        """Sets up the PF coil waveforms.


        This routine sets up the PF coil current waveforms.
        f_c_pf_cs_peak_time_array[i,j] is the current in coil i, at time j,
        normalized to the peak current in that coil at any time.
        """
        nplas = pfcoil_variables.n_cs_pf_coils + 1
        for it in range(6):
            pfcoil_variables.f_c_pf_cs_peak_time_array[nplas - 1, it] = 1.0e0

        for ic in range(pfcoil_variables.n_cs_pf_coils):
            # Find where the peak current occurs
            # Beginning of pulse, t = t_plant_pulse_coil_precharge
            if (
                abs(pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic])
            ) and (
                abs(pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic])
            ):
                pfcoil_variables.c_pf_cs_coils_peak_ma[ic] = (
                    pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic]
                )

            # Beginning of flat-top, t = t_plant_pulse_coil_precharge + t_plant_pulse_plasma_current_ramp_up
            if (
                abs(pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic])
            ) and (
                abs(pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic])
            ):
                pfcoil_variables.c_pf_cs_coils_peak_ma[ic] = (
                    pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic]
                )

            # End of flat-top, t = t_plant_pulse_coil_precharge + t_plant_pulse_plasma_current_ramp_up + t_plant_pulse_fusion_ramp + t_plant_pulse_burn
            if (
                abs(pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic])
            ) and (
                abs(pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic])
                >= abs(pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic])
            ):
                pfcoil_variables.c_pf_cs_coils_peak_ma[ic] = (
                    pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic]
                )

            # Set normalized current waveforms
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 0] = 0.0e0
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 1] = (
                pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 2] = (
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 3] = (
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 4] = (
                pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.f_c_pf_cs_peak_time_array[ic, 5] = 0.0e0


class CSCoil:
    """Calculate central solenoid coil system parameters."""

    def __init__(self, cs_fatigue):
        """Initialise Fortran module variables."""
        self.outfile = constants.NOUT  # output file unit``
        self.mfile = constants.MFILE  # mfile file unit
        self.cs_fatigue = cs_fatigue

    def calculate_cs_geometry(
        self,
        z_tf_inside_half: float,
        f_z_cs_tf_internal: float,
        dr_cs: float,
        dr_bore: float,
    ) -> tuple[float, float, float, float, float, float, float, float, float]:
        """Calculate the geometry of the Central Solenoid (CS) coil.

        Parameters
        ----------
        z_tf_inside_half : float
            Half-height of the TF bore (m)
        f_z_cs_tf_internal : float
            Fractional height of CS relative to TF bore
        dr_cs : float
            Thickness of the CS coil (m)
        dr_bore : float
            Radius of the TF bore (m)

        Returns
        -------
        tuple[float, float, float, float, float, float, float, float]
            Tuple containing:
            - z_cs_coil_upper: Upper Z coordinate of CS coil (m)
            - z_cs_coil_lower: Lower Z coordinate of CS coil (m)
            - r_cs_coil_middle: Radial coordinate of CS coil centre (m)
            - z_cs_coil_middle: Z coordinate of CS coil centre (m)
            - r_cs_coil_outer: Outer radius of CS coil (m)
            - r_cs_coil_inner: Inner radius of CS coil (m)
            - a_cs_poloidal: Total poloidal cross-sectional area of CS coil (m²)
            - dz_cs_full: Full height of CS coil (m)
        """

        # Central Solenoid mean radius
        r_cs_middle = dr_bore + (0.5e0 * dr_cs)

        # Scale the CS height relative to TF bore height
        z_cs_inside_half = z_tf_inside_half * f_z_cs_tf_internal

        dz_cs_full = 2.0e0 * z_cs_inside_half  # Full height of CS coil

        # Z coordinates of CS coil edges
        z_cs_coil_upper = z_cs_inside_half
        z_cs_coil_lower = -z_cs_coil_upper

        # (R,Z) coordinates of coil centre
        r_cs_coil_middle = r_cs_middle
        z_cs_coil_middle = 0.0e0

        # Radius of outer edge
        r_cs_coil_outer = r_cs_middle + 0.5e0 * dr_cs

        # Radius of inner edge
        r_cs_coil_inner = r_cs_coil_outer - dr_cs

        # Full radial thickness of CS coil
        dr_cs_full = 2 * r_cs_coil_outer

        # Total cross-sectional area
        a_cs_poloidal = 2.0e0 * z_cs_inside_half * dr_cs

        return (
            z_cs_coil_upper,
            z_cs_coil_lower,
            r_cs_coil_middle,
            r_cs_middle,
            z_cs_coil_middle,
            r_cs_coil_outer,
            r_cs_coil_inner,
            a_cs_poloidal,
            dz_cs_full,
            dr_cs_full,
        )

    def calculate_cs_turn_geometry_eu_demo(
        self,
        a_cs_turn: float,
        f_dr_dz_cs_turn: float,
        radius_cs_turn_corners: float,
        f_a_cs_turn_steel: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the geometry of a CS (Central Solenoid) turn using the EU DEMO stadium-shaped model.

        Parameters
         ----------
         a_cs_turn : float
             Poloidal area of a CS turn (m^2)
         f_dr_dz_cs_turn : float
             Length-to-height ratio of the CS turn
         radius_cs_turn_corners : float
             Radius of curved outer corner (m)
         f_a_cs_turn_steel : float
             Fraction of steel area in the CS turn

        Returns
        -------
        :
            Tuple containing:
            - dz_cs_turn: Depth/width of CS turn conduit (m)
            - dr_cs_turn: Length of CS turn conduit (m)
            - radius_cs_turn_cable_space: Radius of CS turn cable space (m)
            - dr_cs_turn_conduit: Radial thickness of steel conduit (m)
            - dz_cs_turn_conduit: Vertical thickness of steel conduit (m)

        Notes
        -----
            - The calculation assumes a stadium-shaped cross-section for the CS turn.
            - If the calculated conduit thickness is negative or too small, it is set to a minimum value of 1 mm.

        References
        ----------
            - R. Wesche et al., “Central solenoid winding pack design for DEMO,”
            Fusion Engineering and Design, vol. 124, pp. 82-85, Apr. 2017,
            doi: https://doi.org/10.1016/j.fusengdes.2017.04.052.
        """
        # Vertical height of CS turn conduit/turn
        dz_cs_turn = (a_cs_turn / f_dr_dz_cs_turn) ** 0.5

        # Radial width of CS turn conduit/turn
        dr_cs_turn = f_dr_dz_cs_turn * dz_cs_turn

        # Calculate radius of cable space in CS turn
        radius_cs_turn_cable_space = -((dr_cs_turn - dz_cs_turn) / np.pi) + math.sqrt(
            (((dr_cs_turn - dz_cs_turn) / np.pi) ** 2)
            + (
                (
                    (dr_cs_turn * dz_cs_turn)
                    - (4 - np.pi) * (radius_cs_turn_corners**2)
                    - (a_cs_turn * f_a_cs_turn_steel)
                )
                / np.pi
            )
        )

        # Vertical thickness of steel conduit in CS turn
        dz_cs_turn_conduit = (dz_cs_turn / 2) - radius_cs_turn_cable_space

        # In this model the vertical and radial have the same thickness
        dr_cs_turn_conduit = dz_cs_turn_conduit
        # add a check for negative conduit thickness
        if dr_cs_turn_conduit < 1.0e-3:
            dr_cs_turn_conduit = 1.0e-3
            logger.error("CS turn conduit radial thickness < 1 mm, kludged to 1 mm")

        return (
            dz_cs_turn,
            dr_cs_turn,
            radius_cs_turn_cable_space,
            dr_cs_turn_conduit,
            dz_cs_turn_conduit,
        )

    def place_cs_filaments(
        self,
        n_cs_current_filaments: int,
        r_cs_middle: float,
        z_cs_inside_half: float,
        c_cs_flat_top_end: float,
        f_j_cs_start_pulse_end_flat_top: float,
        nfxf: int,
    ):
        """Places central solenoid (CS) filaments and assigns their positions and currents.

        This function calculates the radial (R) and vertical (Z) positions, as well as the current values,
        for a set of CS filaments based on the provided parameters. Each filament is placed symmetrically
        about the midplane, and currents are assigned according to the flat-top end current and scaling factors.

        Parameters
        ----------
        n_cs_current_filaments : int
            Number of CS current filaments to place (per side).
        r_cs_middle : float
            Radial coordinate of the middle of the CS.
        z_cs_inside_half : float
            Half-height of the CS in the vertical (Z) direction.
        c_cs_flat_top_end : float
            Flat-top end current for the CS.
        f_j_cs_start_pulse_end_flat_top : float
            Scaling factor for the CS current at the start of the pulse and flat-top end.
        nfxf : int
            Number of flux loops or scaling factor for current distribution.

        Returns
        -------
        tuple[list[float], list[float], list[float]]
            Tuple containing:
            - r_pf_cs_current_filaments (list of float): Radial positions of the CS filaments.
            - z_pf_cs_current_filaments (list of float): Vertical positions of the CS filaments.
            - c_pf_cs_current_filaments (list of float): Current values assigned to each CS filament.
        """
        r_pf_cs_current_filaments = np.zeros(pfcoil_variables.NFIXMX)
        z_pf_cs_current_filaments = np.zeros(pfcoil_variables.NFIXMX)
        c_pf_cs_current_filaments = np.zeros(pfcoil_variables.NFIXMX)

        for filament in range(n_cs_current_filaments):
            # Set the R coordinate of the filaments
            r_pf_cs_current_filaments[filament] = r_cs_middle
            r_pf_cs_current_filaments[filament + n_cs_current_filaments] = (
                r_pf_cs_current_filaments[filament]
            )

            # Set the Z cordinate of the filaments
            z_pf_cs_current_filaments[filament] = (
                z_cs_inside_half / n_cs_current_filaments * ((filament + 1) - 0.5e0)
            )
            z_pf_cs_current_filaments[
                filament + n_cs_current_filaments
            ] = -z_pf_cs_current_filaments[filament]

            # Assign currents to the filaments
            c_pf_cs_current_filaments[filament] = (
                -c_cs_flat_top_end / nfxf * f_j_cs_start_pulse_end_flat_top
            )
            c_pf_cs_current_filaments[filament + n_cs_current_filaments] = (
                c_pf_cs_current_filaments[filament]
            )

        return (
            r_pf_cs_current_filaments,
            z_pf_cs_current_filaments,
            c_pf_cs_current_filaments,
        )

    def ohcalc(self):
        """Routine to perform calculations for the Central Solenoid."""

        (
            pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.z_pf_coil_lower[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_cs_middle,
            pfcoil_variables.z_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.a_cs_poloidal,
            pfcoil_variables.dz_cs_full,
            pfcoil_variables.dr_cs_full,
        ) = self.calculate_cs_geometry(
            z_tf_inside_half=bv.z_tf_inside_half,
            f_z_cs_tf_internal=pfcoil_variables.f_z_cs_tf_internal,
            dr_cs=bv.dr_cs,
            dr_bore=bv.dr_bore,
        )

        # Maximum current (MA-turns) in central Solenoid, at either BOP or EOF
        if pfcoil_variables.j_cs_pulse_start > pfcoil_variables.j_cs_flat_top_end:
            sgn = 1.0e0
            pfcoil_variables.c_pf_cs_coils_peak_ma[
                pfcoil_variables.n_cs_pf_coils - 1
            ] = (
                sgn
                * 1.0e-6
                * pfcoil_variables.j_cs_pulse_start
                * pfcoil_variables.a_cs_poloidal
            )
        else:
            sgn = -1.0e0
            pfcoil_variables.c_pf_cs_coils_peak_ma[
                pfcoil_variables.n_cs_pf_coils - 1
            ] = (
                sgn
                * 1.0e-6
                * pfcoil_variables.j_cs_flat_top_end
                * pfcoil_variables.a_cs_poloidal
            )

        # Number of turns
        pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1] = (
            1.0e6
            * abs(
                pfcoil_variables.c_pf_cs_coils_peak_ma[
                    pfcoil_variables.n_cs_pf_coils - 1
                ]
            )
            / pfcoil_variables.c_pf_coil_turn_peak_input[
                pfcoil_variables.n_cs_pf_coils - 1
            ]
        )

        # Turn vertical cross-sectionnal area
        pfcoil_variables.a_cs_turn = (
            pfcoil_variables.a_cs_poloidal
            / pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1]
        )

        (
            pfcoil_variables.dz_cs_turn,
            pfcoil_variables.dr_cs_turn,
            pfcoil_variables.radius_cs_turn_cable_space,
            csfv.dr_cs_turn_conduit,
            csfv.dz_cs_turn_conduit,
        ) = self.calculate_cs_turn_geometry_eu_demo(
            a_cs_turn=pfcoil_variables.a_cs_turn,
            f_dr_dz_cs_turn=pfcoil_variables.f_dr_dz_cs_turn,
            radius_cs_turn_corners=pfcoil_variables.radius_cs_turn_corners,
            f_a_cs_turn_steel=pfcoil_variables.f_a_cs_turn_steel,
        )

        # Non-steel area void fraction for coolant
        pfcoil_variables.f_a_pf_coil_void[pfcoil_variables.n_cs_pf_coils - 1] = (
            pfcoil_variables.f_a_cs_void
        )

        # Peak field at the End-Of-Flattop (EOF)
        # Occurs at inner edge of coil; bmaxoh2 and bzi are of opposite sign at EOF

        # Peak field due to central Solenoid itself
        bmaxoh2 = self.calculate_cs_self_peak_magnetic_field(
            j_cs=pfcoil_variables.j_cs_flat_top_end,
            r_cs_inner=pfcoil_variables.r_pf_coil_inner[
                pfcoil_variables.n_cs_pf_coils - 1
            ],
            r_cs_outer=pfcoil_variables.r_pf_coil_outer[
                pfcoil_variables.n_cs_pf_coils - 1
            ],
            dz_cs_half=pfcoil_variables.z_pf_coil_upper[
                pfcoil_variables.n_cs_pf_coils - 1
            ],
        )

        # Peak field due to other PF coils plus plasma
        timepoint = 5
        _bri, _bro, bzi, bzo = peak_b_field_at_pf_coil(
            n_coil=pfcoil_variables.n_cs_pf_coils,
            n_coil_group=99,
            t_b_field_peak=timepoint,
        )

        pfcoil_variables.b_cs_peak_flat_top_end = abs(bzi - bmaxoh2)

        # Peak field on outboard side of central Solenoid
        # (self-field is assumed to be zero - long solenoid approximation)
        bohco = abs(bzo)

        # Peak field at the Beginning-Of-Pulse (BOP)
        # Occurs at inner edge of coil; b_cs_peak_pulse_start and bzi are of same sign at BOP
        pfcoil_variables.b_cs_peak_pulse_start = (
            self.calculate_cs_self_peak_magnetic_field(
                j_cs=pfcoil_variables.j_cs_pulse_start,
                r_cs_inner=pfcoil_variables.r_pf_coil_inner[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
                r_cs_outer=pfcoil_variables.r_pf_coil_outer[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
                dz_cs_half=pfcoil_variables.z_pf_coil_upper[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
            )
        )
        timepoint = 2
        _bri, _bro, bzi, bzo = peak_b_field_at_pf_coil(
            n_coil=pfcoil_variables.n_cs_pf_coils,
            n_coil_group=99,
            t_b_field_peak=timepoint,
        )

        pfcoil_variables.b_cs_peak_pulse_start = abs(
            pfcoil_variables.b_cs_peak_pulse_start + bzi
        )

        # Maximum field values
        pfcoil_variables.b_pf_coil_peak[pfcoil_variables.n_cs_pf_coils - 1] = max(
            pfcoil_variables.b_cs_peak_flat_top_end,
            abs(pfcoil_variables.b_cs_peak_pulse_start),
        )
        pfcoil_variables.bpf2[pfcoil_variables.n_cs_pf_coils - 1] = max(bohco, abs(bzo))

        # Stress ==> cross-sectional area of supporting steel to use
        if pfcoil_variables.i_pf_conductor == 0:
            # Superconducting coil

            # New calculation from M. N. Wilson for hoop stress
            pfcoil_variables.sig_hoop = self.hoop_stress(
                pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1]
            )

            # New calculation from Y. Iwasa for axial stress
            (
                pfcoil_variables.stress_z_cs_self_peak_midplane,
                pfcoil_variables.forc_z_cs_self_peak_midplane,
            ) = self.calculate_cs_self_peak_midplane_axial_stress(
                r_cs_outer=pfcoil_variables.r_pf_coil_outer[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
                r_cs_inner=pfcoil_variables.r_pf_coil_inner[
                    pfcoil_variables.n_cs_pf_coils - 1
                ],
                dz_cs_half=pfcoil_variables.dz_cs_full / 2.0,
                c_cs_peak=pfcoil_variables.c_pf_cs_coils_peak_ma[
                    pfcoil_variables.n_cs_pf_coils - 1
                ]
                * 1.0e6,
            )

            # Allowable (hoop) stress (Pa) alstroh
            # Now a user input
            # alstroh = min( (2.0e0*csytf/3.0e0), (0.5e0*csutf) )

            # Calculation of CS fatigue
            # this is only valid for pulsed reactor design
            if pv.f_c_plasma_inductive > 0.0e-4:
                csfv.n_cycle, csfv.t_crack_radial = self.cs_fatigue.ncycle(
                    pfcoil_variables.sig_hoop,
                    csfv.residual_sig_hoop,
                    csfv.t_crack_vertical,
                    csfv.dz_cs_turn_conduit,
                    csfv.dr_cs_turn_conduit,
                )

            # Now steel area fraction is iteration variable and constraint
            # equation is used for Central Solenoid stress

            # Area of steel in Central Solenoid
            areaspf = pfcoil_variables.f_a_cs_turn_steel * pfcoil_variables.a_cs_poloidal

            if pfcoil_variables.i_cs_stress == 1:
                pfcoil_variables.s_shear_cs_peak = max(
                    abs(
                        pfcoil_variables.sig_hoop
                        - pfcoil_variables.stress_z_cs_self_peak_midplane
                    ),
                    abs(pfcoil_variables.stress_z_cs_self_peak_midplane - 0.0e0),
                    abs(0.0e0 - pfcoil_variables.sig_hoop),
                )
            else:
                pfcoil_variables.s_shear_cs_peak = max(
                    abs(pfcoil_variables.sig_hoop - 0.0e0),
                    abs(0.0e0 - 0.0e0),
                    abs(0.0e0 - pfcoil_variables.sig_hoop),
                )

            # Thickness of hypothetical steel cylinders assumed to encase the CS along
            # its inside and outside edges; in reality, the steel is distributed
            # throughout the conductor
            pfcoil_variables.pfcaseth[pfcoil_variables.n_cs_pf_coils - 1] = (
                0.25e0
                * areaspf
                / pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]
            )

        else:
            areaspf = 0.0e0  # Resistive Central Solenoid - no steel needed
            pfcoil_variables.pfcaseth[pfcoil_variables.n_cs_pf_coils - 1] = 0.0e0

        # Weight of steel
        pfcoil_variables.m_pf_coil_structure[pfcoil_variables.n_cs_pf_coils - 1] = (
            areaspf
            * 2.0e0
            * np.pi
            * pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]
            * fwbsv.den_steel
        )

        # Non-steel cross-sectional area
        pfcoil_variables.awpoh = pfcoil_variables.a_cs_poloidal - areaspf

        # Issue #97. Fudge to ensure awpoh is positive; result is continuous, smooth and
        # monotonically decreases

        da = 0.0001e0  # 1 cm^2
        if pfcoil_variables.awpoh < da:
            pfcoil_variables.awpoh = da * da / (2.0e0 * da - pfcoil_variables.awpoh)

        # Weight of conductor in central Solenoid
        if pfcoil_variables.i_pf_conductor == 0:
            pfcoil_variables.m_pf_coil_conductor[pfcoil_variables.n_cs_pf_coils - 1] = (
                pfcoil_variables.awpoh
                * (1.0e0 - pfcoil_variables.f_a_cs_void)
                * 2.0e0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]
                * tfv.dcond[pfcoil_variables.i_cs_superconductor - 1]
            )
        else:
            pfcoil_variables.m_pf_coil_conductor[pfcoil_variables.n_cs_pf_coils - 1] = (
                pfcoil_variables.awpoh
                * (1.0e0 - pfcoil_variables.f_a_cs_void)
                * 2.0e0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]
                * constants.den_copper
            )

        if pfcoil_variables.i_pf_conductor == 0:
            # Allowable coil overall current density at EOF
            # (superconducting coils only)

            (
                jcritwp,
                pfcoil_variables.jcableoh_eof,
                pfcoil_variables.j_cs_conductor_critical_flat_top_end,
                tmarg1,
            ) = superconpf(
                pfcoil_variables.b_cs_peak_flat_top_end,
                pfcoil_variables.f_a_cs_void,
                pfcoil_variables.fcuohsu,
                (
                    abs(
                        pfcoil_variables.c_pf_cs_coils_peak_ma[
                            pfcoil_variables.n_cs_pf_coils - 1
                        ]
                    )
                    / pfcoil_variables.awpoh
                )
                * 1.0e6,
                pfcoil_variables.i_cs_superconductor,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )
            # Strand critical current calculation for costing in $/kAm
            # = superconducting filaments jc * (1 - strand copper fraction)
            if pfcoil_variables.i_cs_superconductor in {2, 6, 8}:
                pfcoil_variables.j_crit_str_cs = (
                    pfcoil_variables.j_cs_conductor_critical_flat_top_end
                )
            else:
                pfcoil_variables.j_crit_str_cs = (
                    pfcoil_variables.j_cs_conductor_critical_flat_top_end
                    * (1 - pfcoil_variables.fcuohsu)
                )

            pfcoil_variables.j_cs_critical_flat_top_end = (
                jcritwp * pfcoil_variables.awpoh / pfcoil_variables.a_cs_poloidal
            )

            # Allowable coil overall current density at BOP

            (
                jcritwp,
                pfcoil_variables.jcableoh_bop,
                pfcoil_variables.j_cs_conductor_critical_pulse_start,
                tmarg2,
            ) = superconpf(
                pfcoil_variables.b_cs_peak_pulse_start,
                pfcoil_variables.f_a_cs_void,
                pfcoil_variables.fcuohsu,
                (
                    abs(
                        pfcoil_variables.c_pf_cs_coils_peak_ma[
                            pfcoil_variables.n_cs_pf_coils - 1
                        ]
                    )
                    / pfcoil_variables.awpoh
                )
                * 1.0e6,
                pfcoil_variables.i_cs_superconductor,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )

            pfcoil_variables.j_pf_wp_critical[pfcoil_variables.n_cs_pf_coils - 1] = (
                jcritwp * pfcoil_variables.awpoh / pfcoil_variables.a_cs_poloidal
            )
            pfcoil_variables.j_cs_critical_pulse_start = (
                pfcoil_variables.j_pf_wp_critical[pfcoil_variables.n_cs_pf_coils - 1]
            )

            pfcoil_variables.temp_cs_superconductor_margin = min(tmarg1, tmarg2)

        else:
            # Resistive power losses (non-superconducting coil)

            pfcoil_variables.p_cs_resistive_flat_top = (
                2.0e0
                * np.pi
                * pfcoil_variables.r_cs_middle
                * pfcoil_variables.rho_pf_coil
                / (
                    pfcoil_variables.a_cs_poloidal
                    * (1.0e0 - pfcoil_variables.f_a_cs_void)
                )
                * (
                    1.0e6
                    * pfcoil_variables.c_pf_cs_coils_peak_ma[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                )
                ** 2
            )
            pfcoil_variables.p_pf_coil_resistive_total_flat_top = (
                pfcoil_variables.p_pf_coil_resistive_total_flat_top
                + pfcoil_variables.p_cs_resistive_flat_top
            )

    def calculate_cs_self_peak_magnetic_field(
        self,
        j_cs: float,
        r_cs_inner: float,
        r_cs_outer: float,
        dz_cs_half: float,
    ) -> float:
        """Calculates the maximum field of a solenoid of circular winding and rectangular cross-section.

        Parameters
        ----------
        j_cs : float
            Overall current density (A/m²)
        r_cs_inner : float
            Solenoid inner radius (m)
        r_cs_outer : float
            Solenoid outer radius (m)
        dz_cs_half : float
            Solenoid half height (m)

        Returns
        -------
        float
            Maximum field of solenoid (T)

        References
        ----------
            - Fits are taken from the figure on p.22 of M. Wilson's book
            "Superconducting Magnets", Clarendon Press, Oxford, N.Y., 1983,
            ISBN 13: 9780198548102
        """
        beta = dz_cs_half / r_cs_inner
        alpha = r_cs_outer / r_cs_inner

        # Field at the centre of the bore R=0, Z=0 of the solenoid
        b_cs_bore_centre = (
            j_cs
            * constants.RMU0
            * dz_cs_half
            * math.log(
                (alpha + math.sqrt(alpha**2 + beta**2))
                / (1.0 + math.sqrt(1.0 + beta**2))
            )
        )

        # Fits are for 1 < alpha < 2 , and 0.5 < beta < very large
        if beta > 3.0:
            b1 = constants.RMU0 * j_cs * (r_cs_outer - r_cs_inner)
            f = (3.0 / beta) ** 2
            b_cs_peak = (
                f * b_cs_bore_centre * (1.007 + (alpha - 1.0) * 0.0055) + (1.0 - f) * b1
            )

        elif beta > 2.0:
            rat = (1.025 - (beta - 2.0) * 0.018) + (alpha - 1.0) * (
                0.01 - (beta - 2.0) * 0.0045
            )
            b_cs_peak = rat * b_cs_bore_centre

        elif beta > 1.0:
            rat = (1.117 - (beta - 1.0) * 0.092) + (alpha - 1.0) * (beta - 1.0) * 0.01
            b_cs_peak = rat * b_cs_bore_centre

        elif beta > 0.75:
            rat = (1.30 - 0.732 * (beta - 0.75)) + (alpha - 1.0) * (
                0.2 * (beta - 0.75) - 0.05
            )
            b_cs_peak = rat * b_cs_bore_centre

        else:
            rat = (1.65 - 1.4 * (beta - 0.5)) + (alpha - 1.0) * (
                0.6 * (beta - 0.5) - 0.20
            )
            b_cs_peak = rat * b_cs_bore_centre

        return b_cs_peak

    def output_cs_structure(self):
        op.oheadr(self.outfile, "Central Solenoid Structure")

        op.ocmmnt(self.outfile, "CS turn structure")
        op.oblnkl(self.outfile)

        op.ovarre(
            self.outfile,
            "Poloidal area of a CS turn [m^2]",
            "(a_cs_turn)",
            pfcoil_variables.a_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Radial width a CS turn [m^2]",
            "(dz_cs_turn)",
            pfcoil_variables.dz_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Length of a CS turn [m]",
            "(dr_cs_turn)",
            pfcoil_variables.dr_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Length to diameter ratio of a CS turn",
            "(f_dr_dz_cs_turn)",
            pfcoil_variables.f_dr_dz_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Radius of CS turn cable space [m]",
            "(radius_cs_turn_cable_space)",
            pfcoil_variables.radius_cs_turn_cable_space,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Radial thickness of steel conduit to cable space [m]",
            "(dr_cs_turn_conduit)",
            csfv.dr_cs_turn_conduit,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Vertical thickness of steel conduit to cable space [m]",
            "(dz_cs_turn_conduit)",
            csfv.dz_cs_turn_conduit,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Corner radius of CS turn [m]",
            "(radius_cs_turn_corners)",
            pfcoil_variables.radius_cs_turn_corners,
            "OP ",
        )

    def calculate_cs_self_peak_midplane_axial_stress(
        self,
        r_cs_outer: float,
        r_cs_inner: float,
        dz_cs_half: float,
        c_cs_peak: float,
    ) -> tuple[float, float]:
        """Calculate axial stress and axial force for the central solenoid.

        Parameters
        ----------
        r_cs_outer:
            Outer radius of the central solenoid (m).
        r_cs_inner:
            Inner radius of the central solenoid (m).
        dz_cs_half:
            Half-height of the central solenoid (m).
        c_cs_peak:
            Peak CS coil current (A).

        Returns
        -------
        tuple(float, float)
            A tuple containing the unsmeared axial stress and the axial force.
                The first element is the unsmeared axial stress in MPa.
                The second element is the axial force in newtons (N).

        :note:
           The axial force is computed using elliptic-integral based terms and the
           unsmeared axial stress is obtained by dividing the axial force by
           the effective steel area associated with the CS turns.

        References
        ----------
            - Case Studies in Superconducting Magnets. Boston, MA: Springer US, 2009.
              doi: https://doi.org/10.1007/b112047.
        """

        # kb term for elliptical integrals
        # kb2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + hl**2))
        kb2 = (4.0e0 * r_cs_outer**2) / (4.0e0 * r_cs_outer**2 + dz_cs_half**2)

        # k2b term for elliptical integrals
        # k2b2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + 4.0e0*hl**2))
        k2b2 = (4.0e0 * r_cs_outer**2) / (4.0e0 * r_cs_outer**2 + 4.0e0 * dz_cs_half**2)

        # term 1
        axial_term_1 = (
            -(constants.RMU0 / 2.0e0) * (c_cs_peak / (2.0e0 * dz_cs_half)) ** 2
        )

        # term 2
        ekb2_1 = ellipk(kb2)
        ekb2_2 = ellipe(kb2)
        axial_term_2 = (
            2.0e0
            * dz_cs_half
            * (math.sqrt(4.0e0 * r_cs_outer**2 + dz_cs_half**2))
            * (ekb2_1 - ekb2_2)
        )

        # term 3
        ek2b2_1 = ellipk(k2b2)
        ek2b2_2 = ellipe(k2b2)
        axial_term_3 = (
            2.0e0
            * dz_cs_half
            * (math.sqrt(4.0e0 * r_cs_outer**2 + 4.0e0 * dz_cs_half**2))
            * (ek2b2_1 - ek2b2_2)
        )

        # calculate axial force [N]
        forc_z_cs_self_peak_midplane = axial_term_1 * (axial_term_2 - axial_term_3)

        # axial area [m2]
        area_ax = np.pi * (r_cs_outer**2 - r_cs_inner**2)

        # Calculate unsmeared axial stress
        # Average axial stress at the interface of each half of the coil
        s_axial = forc_z_cs_self_peak_midplane / (0.5 * area_ax)

        return s_axial, forc_z_cs_self_peak_midplane

    def hoop_stress(self, r):
        """Calculation of hoop stress of central solenoid.

        This routine calculates the hoop stress of the central solenoid
        from "Superconducting magnets", M. N. Wilson OUP

        Parameters
        ----------
        r : float
            radial position a < r < b

        Returns
        -------
        float
            hoop stress (MPa)
        """
        a = pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1]

        # Outer radius of central Solenoid [m]
        b = pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]

        # alpha
        alpha = b / a

        # epsilon
        epsilon = r / a

        # Field at inner radius of coil [T]
        b_a = pfcoil_variables.b_cs_peak_pulse_start

        # Field at outer radius of coil [T]
        # Assume to be 0 for now
        b_b = 0.0e0

        # current density [A/m^2]
        j = pfcoil_variables.j_cs_pulse_start

        # K term
        k = ((alpha * b_a - b_b) * j * a) / (alpha - 1.0e0)

        # M term
        m = ((b_a - b_b) * j * a) / (alpha - 1.0e0)

        # calculate hoop stress terms
        hp_term_1 = k * ((2.0e0 + tfv.poisson_steel) / (3.0e0 * (alpha + 1.0e0)))

        hp_term_2 = (
            alpha**2
            + alpha
            + 1.0e0
            + alpha**2 / epsilon**2
            - epsilon
            * (
                ((1.0e0 + 2.0e0 * tfv.poisson_steel) * (alpha + 1.0e0))
                / (2.0e0 + tfv.poisson_steel)
            )
        )

        hp_term_3 = m * ((3.0e0 + tfv.poisson_steel) / (8.0e0))

        hp_term_4 = (
            alpha**2
            + 1.0e0
            + alpha**2 / epsilon**2
            - epsilon**2
            * ((1.0e0 + 3.0e0 * tfv.poisson_steel) / (3.0e0 + tfv.poisson_steel))
        )

        s_hoop_nom = hp_term_1 * hp_term_2 - hp_term_3 * hp_term_4

        return s_hoop_nom / pfcoil_variables.f_a_cs_turn_steel


def peak_b_field_at_pf_coil(
    n_coil: int, n_coil_group: int, t_b_field_peak: int
) -> tuple[float, float, float, float]:
    """Calculates the peak magnetic field components at the inner and outer edges of a given PF coil.

    Parameters
    ----------
    n_coil : int
        Coil number (1-based index)
    n_coil_group : int
        Group number (1-based index)
    t_b_field_peak : int
        Time point at which the field is highest

    Returns
    -------
    tuple[float, float, float, float]
        Tuple containing:
        - b_pf_inner_radial (float): Radial field at inner edge (T)
        - b_pf_outer_radial (float): Radial field at outer edge (T)
        - b_pf_inner_vertical (float): Vertical field at inner edge (T)
        - b_pf_outer_vertical (float): Vertical field at outer edge (T)

    Notes
    -----
    This routine calculates the peak magnetic field components at the inner and outer edges of a given PF coil.
    The calculation includes the effects from all the coils and the plasma.
    """
    if bv.iohcl != 0 and n_coil == pfcoil_variables.n_cs_pf_coils:
        # Peak field is to be calculated at the Central Solenoid itself,
        # so exclude its own contribution; its self field is
        # dealt with externally using routine calculate_cs_self_peak_magnetic_field()
        kk = 0
    else:
        # Check different times for maximum current
        if (
            abs(
                pfcoil_variables.c_pf_cs_coil_pulse_start_ma[n_coil - 1]
                - pfcoil_variables.c_pf_cs_coils_peak_ma[n_coil - 1]
            )
            < 1.0e-12
        ):
            t_b_field_peak = 2
        elif (
            abs(
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[n_coil - 1]
                - pfcoil_variables.c_pf_cs_coils_peak_ma[n_coil - 1]
            )
            < 1.0e-12
        ):
            t_b_field_peak = 4
        elif (
            abs(
                pfcoil_variables.c_pf_cs_coil_pulse_end_ma[n_coil - 1]
                - pfcoil_variables.c_pf_cs_coils_peak_ma[n_coil - 1]
            )
            < 1.0e-12
        ):
            t_b_field_peak = 5
        else:
            raise ProcessValueError(
                "Illegal value of it; possible rounding error",
                t_b_field_peak=t_b_field_peak,
            )

        if bv.iohcl == 0:
            # No Central Solenoid
            kk = 0
        else:
            sgn = (
                1.0
                if pfcoil_variables.j_cs_pulse_start > pfcoil_variables.j_cs_flat_top_end
                else -1.0
            )

            # Current in each filament representing part of the Central Solenoid
            for iohc in range(pfcoil_variables.nfxf):
                pfcoil_variables.c_pf_cs_current_filaments[iohc] = (
                    pfcoil_variables.f_c_pf_cs_peak_time_array[
                        pfcoil_variables.n_cs_pf_coils - 1, t_b_field_peak - 1
                    ]
                    * pfcoil_variables.j_cs_flat_top_end
                    * sgn
                    * bv.dr_cs
                    * pfcoil_variables.f_z_cs_tf_internal
                    * bv.z_tf_inside_half
                    / pfcoil_variables.nfxf
                    * 2.0e0
                )

            kk = pfcoil_variables.nfxf

    # Non-Central Solenoid coils' contributions
    jj = 0
    for iii in range(pfcoil_variables.n_pf_coil_groups):
        for _jjj in range(pfcoil_variables.n_pf_coils_in_group[iii]):
            jj = jj + 1
            # Radius, z-coordinate and current for each coil
            if iii == n_coil_group - 1:
                # Self field from coil (Lyle's Method)
                kk = kk + 1

                dzpf = (
                    pfcoil_variables.z_pf_coil_upper[jj - 1]
                    - pfcoil_variables.z_pf_coil_lower[jj - 1]
                )
                pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.r_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.z_pf_coil_middle[jj - 1] + dzpf * 0.125e0
                )
                pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                    * pfcoil_variables.f_c_pf_cs_peak_time_array[
                        jj - 1, t_b_field_peak - 1
                    ]
                    * 0.25e6
                )
                kk = kk + 1
                pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.r_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.z_pf_coil_middle[jj - 1] + dzpf * 0.375e0
                )
                pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                    * pfcoil_variables.f_c_pf_cs_peak_time_array[
                        jj - 1, t_b_field_peak - 1
                    ]
                    * 0.25e6
                )
                kk = kk + 1
                pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.r_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.z_pf_coil_middle[jj - 1] - dzpf * 0.125e0
                )
                pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                    * pfcoil_variables.f_c_pf_cs_peak_time_array[
                        jj - 1, t_b_field_peak - 1
                    ]
                    * 0.25e6
                )
                kk = kk + 1
                pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.r_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.z_pf_coil_middle[jj - 1] - dzpf * 0.375e0
                )
                pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                    * pfcoil_variables.f_c_pf_cs_peak_time_array[
                        jj - 1, t_b_field_peak - 1
                    ]
                    * 0.25e6
                )

            else:
                # Field from different coil
                kk = kk + 1
                pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.r_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.z_pf_coil_middle[jj - 1]
                )
                pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = (
                    pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                    * pfcoil_variables.f_c_pf_cs_peak_time_array[
                        jj - 1, t_b_field_peak - 1
                    ]
                    * 1.0e6
                )

    # Plasma contribution
    if t_b_field_peak > 2:
        kk = kk + 1
        pfcoil_variables.r_pf_cs_current_filaments[kk - 1] = pv.rmajor
        pfcoil_variables.z_pf_cs_current_filaments[kk - 1] = 0.0e0
        pfcoil_variables.c_pf_cs_current_filaments[kk - 1] = pv.plasma_current

    # Calculate the field at the inner and outer edges
    # of the coil of interest
    pfcoil_variables.xind[:kk], b_pf_inner_radial, b_pf_inner_vertical, _psi = (
        calculate_b_field_at_point(
            r_current_loop=pfcoil_variables.r_pf_cs_current_filaments[:kk],
            z_current_loop=pfcoil_variables.z_pf_cs_current_filaments[:kk],
            c_current_loop=pfcoil_variables.c_pf_cs_current_filaments[:kk],
            r_test_point=pfcoil_variables.r_pf_coil_inner[n_coil - 1],
            z_test_point=pfcoil_variables.z_pf_coil_middle[n_coil - 1],
        )
    )
    pfcoil_variables.xind[:kk], b_pf_outer_radial, b_pf_outer_vertical, _psi = (
        calculate_b_field_at_point(
            r_current_loop=pfcoil_variables.r_pf_cs_current_filaments[:kk],
            z_current_loop=pfcoil_variables.z_pf_cs_current_filaments[:kk],
            c_current_loop=pfcoil_variables.c_pf_cs_current_filaments[:kk],
            r_test_point=pfcoil_variables.r_pf_coil_outer[n_coil - 1],
            z_test_point=pfcoil_variables.z_pf_coil_middle[n_coil - 1],
        )
    )

    # b_pf_coil_peak and bpf2 for the Central Solenoid are calculated in OHCALC
    if (bv.iohcl != 0) and (n_coil == pfcoil_variables.n_cs_pf_coils):
        return (
            b_pf_inner_radial,
            b_pf_outer_radial,
            b_pf_inner_vertical,
            b_pf_outer_vertical,
        )

    bpfin = math.sqrt(b_pf_inner_radial**2 + b_pf_inner_vertical**2)
    bpfout = math.sqrt(b_pf_outer_radial**2 + b_pf_outer_vertical**2)
    for n in range(pfcoil_variables.n_pf_coils_in_group[n_coil_group - 1]):
        pfcoil_variables.b_pf_coil_peak[n_coil - 1 + n] = bpfin
        pfcoil_variables.bpf2[n_coil - 1 + n] = bpfout

    return (
        b_pf_inner_radial,
        b_pf_outer_radial,
        b_pf_inner_vertical,
        b_pf_outer_vertical,
    )


def superconpf(bmax, fhe, fcu, jwp, isumat, fhts, strain, thelium, bcritsc, tcritsc):
    """Routine to calculate the PF coil superconductor properties.

    This routine calculates the superconductor critical winding pack
    current density for the PF coils, plus the temperature margin.
    It is based on the TF coil version, supercon.

    N.B. critical current density for a super conductor (j_crit_sc)
    is for the superconducting strands/tape, not including copper.
    Critical current density for a cable (j_crit_cable) acounts for
    both the fraction of the cable taken up by helium coolant channels,
    and the cable conductor copper fraction - i.e., the copper in the
    superconducting strands AND any addtional copper, such as REBCO
    tape support.

    Parameters
    ----------
    bmax : float
        peak field at conductor (T)
    fhe : float
        fraction of cable space that is for He cooling
    fcu : float
        fraction of cable conductor that is copper
    jwp : float
        actual winding pack current density (A/m2)
    isumat : int
        switch for conductor type
        1 = ITER Nb3Sn, standard parameters,
        2 = Bi-2212 High Temperature Superconductor,
        3 = NbTi,
        4 = ITER Nb3Sn, user-defined parameters
        5 = WST Nb3Sn parameterisation
        7 = Durham Ginzbug-Landau Nb-Ti parameterisation
    fhts : float
        Adjustment factor (<= 1) to account for strain,
        radiation damage, fatigue or AC losses
    strain : float
        Strain on superconductor at operation conditions
    thelium : float
        He temperature at peak field point (K)
    bcritsc : float
        Critical field at zero temperature and strain (T) (isumat=4 only)
    tcritsc : float
        Critical temperature at zero field and strain (K) (isumat=4 only)

    Returns
    -------
    tuple[float, float, float, float]
        Critical winding pack current density (A/m2) (j_crit_wp),
        Critical cable current density (A/m2) (j_crit_cable)
        Superconducting strand non-copper critical current density (A/m2) (j_crit_sc)
        Temperature margin (K) (tmarg)
    """

    # Find critical current density in superconducting strand, jcritstr
    if isumat == 1:
        # ITER Nb3Sn critical surface parameterization
        bc20m = 32.97e0
        tc0m = 16.06e0

        # j_crit_sc returned by superconductors.itersc is the critical current density in the
        # superconductor - not the whole strand, which contains copper

        j_crit_sc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 2:
        # Bi-2212 high temperature superconductor parameterization

        # Current density in a strand of Bi-2212 conductor
        # N.B. jcrit returned by superconductors.bi2212 is the critical current density
        # in the strand, not just the superconducting portion.
        # The parameterization for j_crit_cable assumes a particular strand
        # composition that does not require a user-defined copper fraction,
        # so this is irrelevant in this model

        #  jwp / conductor fraction of cable
        jstrand = jwp / (1.0e0 - fhe)
        j_crit_cable, tmarg = superconductors.bi2212(bmax, jstrand, thelium, fhts)
        #  j_crit_cable / non-copper fraction of conductor
        j_crit_sc = j_crit_cable / (1.0e0 - fcu)

    elif isumat == 3:
        # NbTi data
        bc20m = 15.0e0
        tc0m = 9.3e0
        c0 = 1.0e10
        j_crit_sc, _ = superconductors.jcrit_nbti(thelium, bmax, c0, bc20m, tc0m)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 4:
        # As (1), but user-defined parameters
        bc20m = bcritsc
        tc0m = tcritsc
        j_crit_sc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 5:
        # WST Nb3Sn parameterisation
        bc20m = 32.97e0
        tc0m = 16.06e0

        # j_crit_sc returned by superconductors.itersc is the critical current density in the
        # superconductor - not the whole strand, which contains copper

        j_crit_sc, _, _ = superconductors.western_superconducting_nb3sn(
            thelium, bmax, strain, bc20m, tc0m
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 6:
        # "REBCO" 2nd generation HTS superconductor in CrCo strand
        j_crit_sc, _ = superconductors.jcrit_rebco(thelium, bmax)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 7:
        # Durham Ginzburg-Landau critical surface model for Nb-Ti
        bc20m = tfv.b_crit_upper_nbti
        tc0m = tfv.t_crit_nbti
        j_crit_sc, _, _ = superconductors.gl_nbti(thelium, bmax, strain, bc20m, tc0m)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

        # The CS coil current at EOF
        # c_cs_flat_top_end = bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfcoil_variables.j_cs_flat_top_end

    elif isumat == 8:
        # Durham Ginzburg-Landau critical surface model for REBCO
        bc20m = 429e0
        tc0m = 185e0
        j_crit_sc, _, _ = superconductors.gl_rebco(thelium, bmax, strain, bc20m, tc0m)
        # A0 calculated for tape cross section already
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    elif isumat == 9:
        # Hazelton experimental data + Zhai conceptual model for REBCO
        bc20m = 138
        tc0m = 92
        j_crit_sc, _, _ = superconductors.hijc_rebco(
            thelium,
            bmax,
            bc20m,
            tc0m,
            rcv.dr_hts_tape,
            rcv.dx_hts_tape_rebco,
            rcv.dx_hts_tape_total,
        )
        # A0 calculated for tape cross section already
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

    else:
        # Error condition
        raise ProcessValueError("Illegal value for i_pf_superconductor", isumat=isumat)

    #  Critical current density in winding pack
    jcritwp = j_crit_cable
    #  jwp / conductor fraction of cable
    jstrand = jwp / (1.0e0 - fhe)
    #  jstrand / non-copper fraction of conductor
    jsc = jstrand / (1.0e0 - fcu)

    # Temperature margin (already calculated in superconductors.bi2212 for isumat=2)

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
        t_zero_margin, _root_result = optimize.newton(
            superconductors.superconductor_current_density_margin,
            thelium,
            fprime=None,
            args=arguments,
            tol=1.0e-06,
            maxiter=50,
            fprime2=None,
            x1=another_estimate,
            rtol=1.0e-6,
            full_output=True,
            disp=False,
        )
        tmarg = t_zero_margin - thelium

    return jcritwp, j_crit_cable, j_crit_sc, tmarg


@numba.njit(cache=True)
def calculate_b_field_at_point(
    r_current_loop: np.ndarray,
    z_current_loop: np.ndarray,
    c_current_loop: np.ndarray,
    r_test_point: float,
    z_test_point: float,
) -> tuple[np.ndarray, float, float, float]:
    """Calculate the magnetic field and mutual inductance at a point due to currents in circular poloidal conductor loops.

    Parameters
    ----------
    r_current_loop : np.ndarray
        Array of R coordinates of current loops (m)
    z_current_loop : np.ndarray
        Array of Z coordinates of current loops (m)
    c_current_loop : np.ndarray
        Array of currents in loops (A)
    r_test_point : float
        R coordinate of the test point (m)
    z_test_point : float
        Z coordinate of the test point (m)

    Returns
    -------
    tuple[np.ndarray, float, float, float]
        Tuple containing:
        - ind_mutual_array: Mutual inductances (H) between each loop and the test point
        - b_test_point_radial: Radial field component at the test point (T)
        - b_test_point_vertical: Vertical field component at the test point (T)
        - web_test_point_poloidal: Poloidal flux at the test point (Wb)

    Notes
    -----
    - This routine calculates the magnetic field components and the poloidal flux at a given (R, Z) point,
    given the locations and currents of a set of conductor loops. The mutual inductances between the loops
    and a poloidal filament at the (R, Z) point of interest are also computed.
    """

    #  Elliptic integral coefficients

    a0 = 1.38629436112
    a1 = 0.09666344259
    a2 = 0.03590092383
    a3 = 0.03742563713
    a4 = 0.01451196212
    b0 = 0.5
    b1 = 0.12498593597
    b2 = 0.06880248576
    b3 = 0.03328355346
    b4 = 0.00441787012
    c1 = 0.44325141463
    c2 = 0.06260601220
    c3 = 0.04757383546
    c4 = 0.01736506451
    d1 = 0.24998368310
    d2 = 0.09200180037
    d3 = 0.04069697526
    d4 = 0.00526449639

    n_current_loops = len(r_current_loop)

    ind_mutual_array = np.empty((n_current_loops,))
    b_test_point_radial = 0
    b_test_point_vertical = 0
    web_test_point_poloidal = 0

    for i in range(n_current_loops):
        d = (r_test_point + r_current_loop[i]) ** 2 + (
            z_test_point - z_current_loop[i]
        ) ** 2
        s = 4.0 * r_test_point * r_current_loop[i] / d

        # Kludge: avoid s >= 1.0, a goes inf
        if s > 0.999999:
            s = 0.999999

        t = 1.0 - s
        a = np.log(1.0 / t)

        dz = z_test_point - z_current_loop[i]
        zs = dz**2
        dr = r_test_point - r_current_loop[i]
        sd = np.sqrt(d)

        if dr == 0.0:
            # Kludge to avoid NaNs
            dr = 1e-6

        #  Evaluation of elliptic integrals

        xk = (
            a0
            + t * (a1 + t * (a2 + t * (a3 + a4 * t)))
            + a * (b0 + t * (b1 + t * (b2 + t * (b3 + b4 * t))))
        )
        xe = (
            1.0
            + t * (c1 + t * (c2 + t * (c3 + c4 * t)))
            + a * t * (d1 + t * (d2 + t * (d3 + d4 * t)))
        )

        #  Mutual inductances

        ind_mutual_array[i] = 0.5 * constants.RMU0 * sd * ((2.0 - s) * xk - 2.0 * xe)

        #  Radial, vertical fields

        brx = (
            constants.RMU0
            * c_current_loop[i]
            * dz
            / (2 * np.pi * r_test_point * sd)
            * (-xk + (r_current_loop[i] ** 2 + r_test_point**2 + zs) / (dr**2 + zs) * xe)
        )
        bzx = (
            constants.RMU0
            * c_current_loop[i]
            / (2 * np.pi * sd)
            * (xk + (r_current_loop[i] ** 2 - r_test_point**2 - zs) / (dr**2 + zs) * xe)
        )

        #  Sum fields, flux

        b_test_point_radial += brx
        b_test_point_vertical += bzx
        web_test_point_poloidal += ind_mutual_array[i] * c_current_loop[i]

    return (
        ind_mutual_array,
        b_test_point_radial,
        b_test_point_vertical,
        web_test_point_poloidal,
    )


@numba.njit(cache=True)
def rsid(npts, brin, bzin, nfix, n_pf_coil_groups, ccls, bfix, gmat):
    """Computes the norm of the residual vectors.

    This routine calculates the residuals from the matrix
    equation for calculation of the currents in a group of ring coils.

    Parameters
    ----------
    npts : int
        number of data points at which field is  to be fixed;
        should be <= nptsmx
    brin : numpy.ndarray
        field components at data points (T)
    bzin : numpy.ndarray
        field components at data points (T)
    nfix : int
        number of coils with fixed currents, <= nfixmx
    n_pf_coil_groups : int
        number of coil groups, where all coils in a group have the
        same current, <= n_pf_groups_max
    ccls : numpy.ndarray
        coil currents in each group (A)
    bfix : numpy.ndarray
        work array
    gmat : numpy.ndarray
        work array

    Returns
    -------
    tuple[float, float, float, float, float]
        sum of squares of radial field residues (brssq), radial field
        residue norm (brnrm), sum of squares of vertical field residues (bzssq),
        vertical field residue norm (bznrm), sum of squares of elements of
        residual vector (ssq)
    """
    brnrm = 0.0e0
    brssq = 0.0e0

    for i in range(npts):
        svec = 0.0e0
        if nfix > 0:
            svec = bfix[i]

        for j in range(n_pf_coil_groups):
            svec = svec + gmat[i, j] * ccls[j]

        rvec = svec - brin[i]
        brnrm = brnrm + brin[i] ** 2
        brssq = brssq + rvec**2

    bznrm = 0.0e0
    bzssq = 0.0e0

    for i in range(npts):
        svec = 0.0e0
        if nfix > 0:
            svec = bfix[i + npts]
        for j in range(n_pf_coil_groups):
            svec = svec + gmat[i + npts, j] * ccls[j]

        rvec = svec - bzin[i]
        bznrm = bznrm + bzin[i] ** 2
        bzssq = bzssq + rvec**2

    ssq = brssq / (1.0e0 + brnrm) + bzssq / (1.0e0 + bznrm)

    return brssq, brnrm, bzssq, bznrm, ssq


@numba.njit(cache=True)
def fixb(lrow1, npts, rpts, zpts, nfix, rfix, zfix, cfix):
    """Calculates the field from the fixed current loops.

    This routine calculates the fields at the points specified by
    (rpts,zpts) from the set of coils with fixed currents.

    Parameters
    ----------
    lrow1 : int
        row length of array bfix; should be >= nptsmx
    npts : int
        number of data points at which field is to be fixed;
        should be <= nptsmx
    rpts : numpy.ndarray
        coords of data points (m)
    zpts : numpy.ndarray
        coords of data points (m)
    nfix : int
        number of coils with fixed currents, <= nfixmx
    rfix : numpy.ndarray
        coordinates of coils with fixed currents (m)
    zfix : numpy.ndarray
        coordinates of coils with fixed currents (m)
    cfix : numpy.ndarray
        Fixed currents (A)

    Returns
    -------
    numpy.ndarray
        Fields at data points (T)
    """
    bfix = np.zeros(lrow1)

    if nfix <= 0:
        return bfix

    for i in range(npts):
        # calculate_b_field_at_point() only operates correctly on nfix slices of array
        # arguments, not entire arrays
        _, brw, bzw, _ = calculate_b_field_at_point(
            r_current_loop=rfix[:nfix],
            z_current_loop=zfix[:nfix],
            c_current_loop=cfix[:nfix],
            r_test_point=rpts[i],
            z_test_point=zpts[i],
        )
        bfix[i] = brw
        bfix[npts + i] = bzw

    return bfix


@numba.njit(cache=True)
def mtrx(
    lrow1,
    lcol1,
    npts,
    rpts,
    zpts,
    brin,
    bzin,
    n_pf_coil_groups,
    n_pf_coils_in_group,
    r_pf_coil_middle_group_array,
    z_pf_coil_middle_group_array,
    alfa,
    bfix,
    n_pf_coils_in_group_max,
):
    """Calculate the currents in a group of ring coils.

    Set up the matrix equation to calculate the currents in a group of ring
    coils.

    Parameters
    ----------
    lrow1 : int
        row length of arrays bfix, bvec, gmat, umat, vmat; should
        be >= (2*nptsmx + n_pf_groups_max)
    lcol1 : int
        column length of arrays gmat, umat, vmat; should be >=
        n_pf_groups_max
    npts : int
        number of data points at which field is to be fixed; should
        be <= nptsmx
    rpts : numpy.ndarray
        coords of data points (m)
    zpts : numpy.ndarray
        coords of data points (m)
    brin : numpy.ndarray
        field components at data points (T)
    bzin : numpy.ndarray
        field components at data points (T)
    n_pf_coil_groups : int
        number of coil groups, where all coils in a group have the
        same current, <= n_pf_groups_max
    n_pf_coils_in_group : numpy.ndarray
        number of coils in each group, each value <= n_pf_coils_in_group_max
    r_pf_coil_middle_group_array : numpy.ndarray
        coords R(i,j), Z(i,j) of coil j in group i (m)
    z_pf_coil_middle_group_array : numpy.ndarray
        coords R(i,j), Z(i,j) of coil j in group i (m)
    alfa : float
        smoothing parameter (0 = no smoothing, 1.0D-9 = large
        smoothing)
    bfix : numpy.ndarray
        Fields at data points (T)
    n_pf_coils_in_group_max :


    Returns
    -------
    :
        actual number of rows to use, work array, work array,
        Coordinates of conductor loops (m), Coordinates of conductor loops (m),
        Currents in conductor loops (A), Mutual inductances (H)
    """
    bvec = np.zeros(lrow1)
    gmat = np.zeros((lrow1, lcol1))
    cc = np.ones(n_pf_coils_in_group_max)

    for i in range(npts):
        bvec[i] = brin[i] - bfix[i]
        bvec[i + npts] = bzin[i] - bfix[i + npts]

        for j in range(n_pf_coil_groups):
            nc = n_pf_coils_in_group[j]

            _, gmat[i, j], gmat[i + npts, j], _ = calculate_b_field_at_point(
                r_current_loop=r_pf_coil_middle_group_array[j, :nc],
                z_current_loop=z_pf_coil_middle_group_array[j, :nc],
                c_current_loop=cc[:nc],
                r_test_point=rpts[i],
                z_test_point=zpts[i],
            )

    # Add constraint equations
    nrws = 2 * npts

    bvec[nrws : nrws + n_pf_coil_groups] = 0.0
    np.fill_diagonal(
        gmat[nrws : nrws + n_pf_coil_groups, :n_pf_coil_groups],
        n_pf_coils_in_group[:n_pf_coil_groups] * alfa,
    )

    nrws = 2 * npts + n_pf_coil_groups

    # numba doesnt like np.zeros(..., order="F") so this acts as a work
    # around to that missing signature
    gmat = np.asfortranarray(gmat)

    return nrws, gmat, bvec
