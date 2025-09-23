import logging
import math

import numba
import numpy as np
from scipy import optimize
from scipy.linalg import svd
from scipy.special import ellipe, ellipk

import process.superconductors as superconductors
from process import constants
from process import process_output as op
from process.data_structure import build_variables as bv
from process.data_structure import constraint_variables as ctv
from process.data_structure import cs_fatigue_variables as csfv
from process.data_structure import fwbs_variables as fwbsv
from process.data_structure import (
    numerics,
    pfcoil_variables,
    superconducting_tf_coil_variables,
)
from process.data_structure import physics_variables as pv
from process.data_structure import rebco_variables as rcv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class PFCoil:
    """Calculate poloidal field coil system parameters."""

    def __init__(self, cs_fatigue) -> None:
        """Initialise Fortran module variables."""
        self.outfile = constants.NOUT  # output file unit
        self.mfile = constants.MFILE  # mfile file unit
        pfcoil_variables.init_pfcoil_module()
        self.cs_fatigue = cs_fatigue

    def run(self):
        """Run the PF coil model."""
        self.pfcoil()

        # Poloidal field coil inductance calculation
        self.induct(False)

        # Volt-second capability of PF coil set
        self.vsec()

    def output(self):
        """Output results to output file."""
        self.output_cs_structure()
        self.outpf()
        self.outvolt()

    def output_induct(self):
        """Output poloidal field coil inductance calculation."""
        self.induct(True)

    def pfcoil(self):
        """Routine to perform calculations for the PF and Central Solenoid coils.
        author: P J Knight, CCFE, Culham Science Centre
        author: R Kemp, CCFE, Culham Science Centre
        None
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
        sigma, work2 = np.zeros((2, pfcoil_variables.N_PF_GROUPS_MAX))
        rc, zc, cc, xc = np.zeros((4, pfcoil_variables.N_PF_COILS_IN_GROUP_MAX))
        brin, bzin, rpts, zpts = np.zeros((4, pfcoil_variables.NPTSMX))
        bfix, bvec = np.zeros((2, lrow1))
        gmat, umat, vmat = np.zeros((3, lrow1, lcol1), order="F")
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
        tv.tim[0] = 0.0e0
        tv.tim[1] = tv.t_precharge
        tv.tim[2] = tv.tim[1] + tv.t_current_ramp_up
        tv.tim[3] = tv.tim[2] + tv.t_fusion_ramp
        tv.tim[4] = tv.tim[3] + tv.t_burn
        tv.tim[5] = tv.tim[4] + tv.t_ramp_down

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
        ) = self.calculate_cs_geometry(
            z_tf_inside_half=bv.z_tf_inside_half,
            f_z_cs_tf_internal=pfcoil_variables.f_z_cs_tf_internal,
            dr_cs=bv.dr_cs,
            dr_bore=bv.dr_bore,
        )

        # nfxf is the total no of filaments into which the Central Solenoid is split,
        # if present
        if bv.iohcl == 0:
            pfcoil_variables.nfxf = 0
            ioheof = 0.0e0
        else:
            pfcoil_variables.nfxf = 2 * pfcoil_variables.nfxfh

            # total Central Solenoid current at EOF
            ioheof = (
                -bv.z_tf_inside_half
                * pfcoil_variables.f_z_cs_tf_internal
                * bv.dr_cs
                * 2.0e0
                * pfcoil_variables.j_cs_flat_top_end
            )

            if pfcoil_variables.nfxf > pfcoil_variables.NFIXMX:
                raise ProcessValueError(
                    "Too many filaments nfxf repesenting the OH coil",
                    nfxf=pfcoil_variables.nfxf,
                    nfixmx=pfcoil_variables.NFIXMX,
                )

            # Symmetric up/down Central Solenoid : Find (R,Z) and current of each filament at BOP

            for nng in range(pfcoil_variables.nfxfh):
                pfcoil_variables.rfxf[nng] = pfcoil_variables.r_cs_middle
                pfcoil_variables.rfxf[nng + pfcoil_variables.nfxfh] = (
                    pfcoil_variables.rfxf[nng]
                )
                pfcoil_variables.zfxf[nng] = (
                    bv.z_tf_inside_half
                    * pfcoil_variables.f_z_cs_tf_internal
                    / pfcoil_variables.nfxfh
                    * ((nng + 1) - 0.5e0)
                )
                pfcoil_variables.zfxf[
                    nng + pfcoil_variables.nfxfh
                ] = -pfcoil_variables.zfxf[nng]
                pfcoil_variables.cfxf[nng] = (
                    -ioheof
                    / pfcoil_variables.nfxf
                    * pfcoil_variables.f_j_cs_start_pulse_end_flat_top
                )
                pfcoil_variables.cfxf[nng + pfcoil_variables.nfxfh] = (
                    pfcoil_variables.cfxf[nng]
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
                # PF coil is in general location
                # See issue 1418
                # https://git.ccfe.ac.uk/process/process/-/issues/1418
                for coil in range(pfcoil_variables.n_pf_coils_in_group[group]):
                    pfcoil_variables.z_pf_coil_middle_group_array[group, coil] = (
                        pv.rminor * pfcoil_variables.zref[group] * signn[coil]
                    )
                    pfcoil_variables.r_pf_coil_middle_group_array[group, coil] = (
                        pv.rminor * pfcoil_variables.rref[group] + pv.rmajor
                    )

            else:
                raise ProcessValueError(
                    "Illegal i_pf_location value",
                    group=group,
                    i_pf_location=pfcoil_variables.i_pf_location[group],
                )

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
                pfcoil_variables.rfxf,
                pfcoil_variables.zfxf,
                pfcoil_variables.cfxf,
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
                pv.bvert = (
                    -1.0e-7
                    * pv.plasma_current
                    / pv.rmajor
                    * (
                        math.log(8.0e0 * pv.aspect)
                        + pv.beta_poloidal
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
                            pfcoil_variables.rfxf[nocoil] = (
                                pfcoil_variables.r_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.zfxf[nocoil] = (
                                pfcoil_variables.z_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.cfxf[nocoil] = pfcoil_variables.ccls[i]
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
                            pfcoil_variables.rfxf[nocoil] = (
                                pfcoil_variables.r_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.zfxf[nocoil] = (
                                pfcoil_variables.z_pf_coil_middle_group_array[i, ccount]
                            )
                            pfcoil_variables.cfxf[nocoil] = pfcoil_variables.ccls[i]
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
                        + pv.beta_poloidal
                        + (pv.ind_plasma_internal_norm / 2.0e0)
                        - 1.5e0
                    )
                )

                pv.bvert = bzin[0]

                ssqef, pfcoil_variables.ccls0 = self.efc(
                    npts0,
                    rpts,
                    zpts,
                    brin,
                    bzin,
                    nfxf0,
                    pfcoil_variables.rfxf,
                    pfcoil_variables.zfxf,
                    pfcoil_variables.cfxf,
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
                * constants.PI
                * constants.PI
                * (
                    (bv.dr_bore * bv.dr_bore)
                    + (bv.dr_cs * bv.dr_cs) / 6.0e0
                    + (bv.dr_cs * bv.dr_bore) / 2.0e0
                )
                / (bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * 2.0e0)
            )
            dics = csflux / ddics

            pfcoil_variables.f_j_cs_start_end_flat_top = (
                (-ioheof * pfcoil_variables.f_j_cs_start_pulse_end_flat_top) + dics
            ) / ioheof
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

                # Beginning of pulse: t = tv.t_precharge
                pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ncl] = (
                    1.0e-6 * pfcoil_variables.ccl0[nng]
                )

                # Beginning of flat-top: t = tv.t_precharge+tv.t_current_ramp_up
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ncl] = 1.0e-6 * (
                    pfcoil_variables.ccls[nng]
                    - (
                        pfcoil_variables.ccl0[nng]
                        * pfcoil_variables.f_j_cs_start_end_flat_top
                        / pfcoil_variables.f_j_cs_start_pulse_end_flat_top
                    )
                )

                # End of flat-top: t = tv.t_precharge+tv.t_current_ramp_up+tv.t_fusion_ramp+tv.t_burn
                pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * (
                    pfcoil_variables.ccls[nng]
                    - (
                        pfcoil_variables.ccl0[nng]
                        * (1.0e0 / pfcoil_variables.f_j_cs_start_pulse_end_flat_top)
                    )
                )

                ncl = ncl + 1

        # Current in Central Solenoid as a function of time
        # N.B. If the Central Solenoid is not present then ioheof is zero.
        pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ncl] = (
            -1.0e-6 * ioheof * pfcoil_variables.f_j_cs_start_pulse_end_flat_top
        )
        pfcoil_variables.c_pf_cs_coil_flat_top_ma[ncl] = (
            1.0e-6 * ioheof * pfcoil_variables.f_j_cs_start_end_flat_top
        )
        pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * ioheof

        # Set up coil current waveforms, normalised to the peak current in
        # each coil
        self.waveform()  # sets c_pf_cs_coils_peak_ma(), waves()

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
                    bri, bro, bzi, bzo = self.peakb(
                        i + 1, iii + 1, it
                    )  # returns b_pf_coil_peak, bpf2

                # Issue 1871.  MDK
                # Allowable current density (for superconducting coils) for each coil, index i
                if pfcoil_variables.i_pf_conductor == 0:
                    bmax = max(
                        abs(pfcoil_variables.b_pf_coil_peak[i]),
                        abs(pfcoil_variables.bpf2[i]),
                    )

                    pfcoil_variables.j_pf_wp_critical[i], jstrand, jsc, tmarg = (
                        self.superconpf(
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
                    * constants.PI
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
                    * constants.PI
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
            self.ohcalc()

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
                pfcoil_variables.c_pf_coil_turn[i, k] = pfcoil_variables.waves[
                    i, k
                ] * math.copysign(
                    pfcoil_variables.c_pf_coil_turn_peak_input[i],
                    pfcoil_variables.c_pf_cs_coils_peak_ma[i],
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
        """
        Calculate the placement of PF coils stacked above the Central Solenoid.

        :param n_pf_coils_in_group: Array containing the number of coils in each PF group.
        :type n_pf_coils_in_group: np.ndarray
        :param n_pf_group: Index of the PF coil group.
        :type n_pf_group: int
        :param r_cs_middle: Radial coordinate of CS coil centre (m).
        :type r_cs_middle: float
        :param dr_pf_cs_middle_offset: Radial offset for PF coil placement (m).
        :type dr_pf_cs_middle_offset: float
        :param z_tf_inside_half: Half-height of the TF bore (m).
        :type z_tf_inside_half: float
        :param dr_tf_inboard: Thickness of the TF inboard leg (m).
        :type dr_tf_inboard: float
        :param z_cs_coil_upper: Upper z coordinate of the CS coil (m).
        :type z_cs_coil_upper: float
        :return: Tuple of arrays containing the radial and vertical coordinates of PF coils in the group.
        :rtype: tuple[np.ndarray, np.ndarray]
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
        """
        Calculates and places poloidal field (PF) coils above the toroidal field (TF) coils for a given group.

        :param n_pf_coils_in_group: Array containing the number of PF coils in each group.
        :type n_pf_coils_in_group: np.ndarray
        :param n_pf_group: Index of the PF coil group to process.
        :type n_pf_group: int
        :param rmajor: Major radius of the device.
        :type rmajor: float
        :param triang: Triangularity parameter for coil placement.
        :type triang: float
        :param rminor: Minor radius of the device.
        :type rminor: float
        :param itart: Flag indicating ST configuration.
        :type itart: int
        :param itartpf: Flag indicating PF coil configuration for ST.
        :type itartpf: int
        :param z_tf_inside_half: Half-height of the TF coil inside region.
        :type z_tf_inside_half: float
        :param dz_tf_upper_lower_midplane: Height difference parameter for PF coil placement.
        :type dz_tf_upper_lower_midplane: float
        :param z_tf_top: Top z-coordinate of the TF coil.
        :type z_tf_top: float
        :param top_bottom: Indicator for coil placement above (+1) or below (-1) the midplane.
        :type top_bottom: int
        :param rpf2: Radial offset parameter for PF coil placement.
        :type rpf2: float
        :param zref: Array of reference z-coordinates for PF coil placement.
        :type zref: np.ndarray

        :returns: Tuple containing arrays of radial and vertical positions of PF coil middles for the specified group,
                  and the updated top_bottom indicator.
        :rtype: tuple[np.ndarray, np.ndarray, int]
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
        """
        Calculates the radial and vertical positions of poloidal field (PF) coils placed outside the toroidal field (TF) coil.

        :param n_pf_coils_in_group: Array containing the number of PF coils in each group.
        :type n_pf_coils_in_group: np.ndarray
        :param n_pf_group: Index of the PF coil group to process.
        :type n_pf_group: int
        :param rminor: Minor radius of the device.
        :type rminor: float
        :param zref: Reference vertical positions for each PF coil group.
        :type zref: np.ndarray
        :param i_tf_shape: Integer flag indicating TF coil shape (2 for picture frame, others for D-shape).
        :type i_tf_shape: int
        :param i_r_pf_outside_tf_placement: Placement switch for PF coil radius (1 for constant/stacked, 0 for following TF curve).
        :type i_r_pf_outside_tf_placement: int
        :param r_pf_outside_tf_midplane: Radial position of PF coil at the midplane.
        :type r_pf_outside_tf_midplane: float

        :returns: Tuple containing arrays of radial and vertical positions of PF coil centers for the specified group.
        :rtype: tuple[np.ndarray, np.ndarray]
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
                r_pf_coil_middle_group_array[n_pf_group, coil] = (
                    r_pf_outside_tf_midplane
                )
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

        author: P J Knight, CCFE, Culham Science Centre
        author: D Strickler, ORNL
        author: J Galambos, ORNL
        author: P C Shipe, ORNL
        This routine calculates the currents required in a group
        of ring coils to produce a fixed field at prescribed
        locations. Additional ring coils with fixed currents are
        also allowed.

        :param npts: number of data points at which field is to be fixed; should
        be <= nptsmx
        :type npts: int
        :param rpts: coords of data points (m)
        :type rpts: np.ndarray
        :param zpts: coords of data points (m)
        :type zpts: np.ndarray
        :param brin: field components at data points (T)
        :type brin: np.ndarray
        :param bzin: field components at data points (T)
        :type bzin: np.ndarray
        :param nfix: number of coils with fixed currents, <= nfixmx
        :type nfix: int
        :param rfix: coordinates of coils with fixed currents (m)
        :type rfix: np.ndarray
        :param zfix: coordinates of coils with fixed currents (m)
        :type zfix: np.ndarray
        :param cfix: Fixed currents (A)
        :type cfix: np.ndarray
        :param n_pf_coil_groups: number of coil groups, where all coils in a group have the
        same current, <= n_pf_groups_max
        :type n_pf_coil_groups: int
        :param n_pf_coils_in_group: number of coils in each group, each value <= n_pf_coils_in_group_max
        :type n_pf_coils_in_group: np.ndarray
        :param r_pf_coil_middle_group_array: coords R(i,j), Z(i,j) of coil j in group i (m)
        :type r_pf_coil_middle_group_array: np.ndarray
        :param z_pf_coil_middle_group_array: coords R(i,j), Z(i,j) of coil j in group i (m)
        :type z_pf_coil_middle_group_array: np.ndarray
        :param alfa: smoothing parameter (0 = no smoothing, 1.0D-9 = large
        smoothing)
        :type alfa: float
        :param bfix: work array
        :type bfix: np.ndarray
        :param gmat: work array
        :type gmat: np.ndarray
        :param bvec: work array
        :type bvec: np.ndarray
        :param rc: work array
        :type rc: np.ndarray
        :param zc: work array
        :type zc: np.ndarray
        :param cc: work array
        :type cc: np.ndarray
        :param xc: work array
        :type xc: np.ndarray
        :param umat: work array
        :type umat: np.ndarray
        :param vmat: work array
        :type vmat: np.ndarray
        :param sigma: work array
        :type sigma: np.ndarray
        :param work2: work array
        :type work2: np.ndarray
        :return: sum of squares of elements of residual vector, solution vector
        of coil currents in each group (A)
        :rtype: tuple[float, np.ndarray]
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
        brssq, brnrm, bzssq, bznrm, ssq = rsid(
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
        author: P J Knight, CCFE, Culham Science Centre
        author: D Strickler, ORNL
        author: J Galambos, ORNL
        author: P C Shipe, ORNL

        :param n_pf_groups_max: maximum number of PF coil groups
        :type n_pf_groups_max: int
        :param n_pf_coil_groups: number of coil groups, where all coils in a group have the
        same current, <= n_pf_groups_max
        :type n_pf_coil_groups: int
        :param nrws: actual number of rows to use
        :type nrws: int
        :param gmat: work array
        :type gmat: numpy.ndarray
        :param bvec: work array
        :type bvec: numpy.ndarray
        :return: solution vector of coil currents
        in each group (A) (ccls), rest are work arrays
        :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray,
        numpy.ndarray, numpy.ndarray]
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

    def ohcalc(self):
        """Routine to perform calculations for the Central Solenoid.

        author: P J Knight, CCFE, Culham Science Centre
        """

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

        # Depth/width of cs turn conduit
        pfcoil_variables.dz_cs_turn = (
            pfcoil_variables.a_cs_turn / pfcoil_variables.ld_ratio_cst
        ) ** 0.5

        # length of cs turn conduit
        pfcoil_variables.dr_cs_turn = (
            pfcoil_variables.ld_ratio_cst * pfcoil_variables.dz_cs_turn
        )

        # Radius of turn space = pfcoil_variables.radius_cs_turn_cable_space
        # Radius of curved outer corrner pfcoil_variables.r_out_cst = 3mm from literature
        # pfcoil_variables.ld_ratio_cst = 70 / 22 from literature

        # CS coil turn geometry calculation - stadium shape
        # Literature: https://doi.org/10.1016/j.fusengdes.2017.04.052
        pfcoil_variables.radius_cs_turn_cable_space = -(
            (pfcoil_variables.dr_cs_turn - pfcoil_variables.dz_cs_turn) / constants.PI
        ) + math.sqrt(
            (
                (
                    (pfcoil_variables.dr_cs_turn - pfcoil_variables.dz_cs_turn)
                    / constants.PI
                )
                ** 2
            )
            + (
                (
                    (pfcoil_variables.dr_cs_turn * pfcoil_variables.dz_cs_turn)
                    - (4 - constants.PI) * (pfcoil_variables.r_out_cst**2)
                    - (pfcoil_variables.a_cs_turn * pfcoil_variables.f_a_cs_steel)
                )
                / constants.PI
            )
        )

        # Thickness of steel conduit in cs turn
        csfv.t_structural_vertical = (
            pfcoil_variables.dz_cs_turn / 2
        ) - pfcoil_variables.radius_cs_turn_cable_space
        # In this model the vertical and radial have the same thickness
        csfv.t_structural_radial = csfv.t_structural_vertical
        # add a check for negative conduit thickness
        if csfv.t_structural_radial < 1.0e-3:
            csfv.t_structural_radial = 1.0e-3

        # Non-steel area void fraction for coolant
        pfcoil_variables.f_a_pf_coil_void[pfcoil_variables.n_cs_pf_coils - 1] = (
            pfcoil_variables.f_a_cs_void
        )

        # Peak field at the End-Of-Flattop (EOF)
        # Occurs at inner edge of coil; bmaxoh2 and bzi are of opposite sign at EOF

        # Peak field due to central Solenoid itself
        bmaxoh2 = self.bfmax(
            pfcoil_variables.j_cs_flat_top_end,
            pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1],
        )

        # Peak field due to other PF coils plus plasma
        timepoint = 5
        bri, bro, bzi, bzo = self.peakb(pfcoil_variables.n_cs_pf_coils, 99, timepoint)

        pfcoil_variables.b_cs_peak_flat_top_end = abs(bzi - bmaxoh2)

        # Peak field on outboard side of central Solenoid
        # (self-field is assumed to be zero - long solenoid approximation)
        bohco = abs(bzo)

        # Peak field at the Beginning-Of-Pulse (BOP)
        # Occurs at inner edge of coil; b_cs_peak_pulse_start and bzi are of same sign at BOP
        pfcoil_variables.b_cs_peak_pulse_start = self.bfmax(
            pfcoil_variables.j_cs_pulse_start,
            pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1],
            pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1],
        )
        timepoint = 2
        bri, bro, bzi, bzo = self.peakb(pfcoil_variables.n_cs_pf_coils, 99, timepoint)

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
            pfcoil_variables.sig_axial, pfcoil_variables.axial_force = (
                self.axial_stress()
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
                    csfv.t_structural_vertical,
                    csfv.t_structural_radial,
                )

            # Now steel area fraction is iteration variable and constraint
            # equation is used for Central Solenoid stress

            # Area of steel in Central Solenoid
            areaspf = pfcoil_variables.f_a_cs_steel * pfcoil_variables.a_cs_poloidal

            if pfcoil_variables.i_cs_stress == 1:
                pfcoil_variables.s_shear_cs_peak = max(
                    abs(pfcoil_variables.sig_hoop - pfcoil_variables.sig_axial),
                    abs(pfcoil_variables.sig_axial - 0.0e0),
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
            * constants.PI
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
                * constants.PI
                * pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]
                * tfv.dcond[pfcoil_variables.i_cs_superconductor - 1]
            )
        else:
            pfcoil_variables.m_pf_coil_conductor[pfcoil_variables.n_cs_pf_coils - 1] = (
                pfcoil_variables.awpoh
                * (1.0e0 - pfcoil_variables.f_a_cs_void)
                * 2.0e0
                * constants.PI
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
            ) = self.superconpf(
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
            ) = self.superconpf(
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
                * constants.PI
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

    def calculate_cs_geometry(
        self,
        z_tf_inside_half: float,
        f_z_cs_tf_internal: float,
        dr_cs: float,
        dr_bore: float,
    ) -> tuple[float, float, float, float, float, float, float, float, float]:
        """
        Calculate the geometry of the Central Solenoid (CS) coil.

        :param z_tf_inside_half: Half-height of the TF bore (m)
        :type z_tf_inside_half: float
        :param f_z_cs_tf_internal: Fractional height of CS relative to TF bore
        :type f_z_cs_tf_internal: float
        :param dr_cs: Thickness of the CS coil (m)
        :type dr_cs: float
        :param dr_bore: Radius of the TF bore (m)
        :type dr_bore: float
        :return: Tuple containing:
            - z_cs_coil_upper: Upper Z coordinate of CS coil (m)
            - z_cs_coil_lower: Lower Z coordinate of CS coil (m)
            - r_cs_coil_middle: Radial coordinate of CS coil centre (m)
            - z_cs_coil_middle: Z coordinate of CS coil centre (m)
            - r_cs_coil_outer: Outer radius of CS coil (m)
            - r_cs_coil_inner: Inner radius of CS coil (m)
            - a_cs_poloidal: Total poloidal cross-sectional area of CS coil (m²)
            - dz_cs_full: Full height of CS coil (m)
        :rtype: tuple[float, float, float, float, float, float, float, float]
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

    def peakb(self, i, ii, it):
        """Calculates the peak field at a PF coil.

        author: P J Knight, CCFE, Culham Science Centre
        i : input integer : coil number
        ii : input integer : group number
        it : input integer : time point at which field is highest
        bri : output real : radial field at inner edge (T)
        bro : output real : radial field at outer edge (T)
        bzi : output real : vertical field at inner edge (T)
        bzo : output real : vertical field at outer edge (T)
        This routine calculates the peak magnetic field components
        at the inner and outer edges of a given PF coil.
        The calculation includes the effects from all the coils
        and the plasma.
        """
        if bv.iohcl != 0 and i == pfcoil_variables.n_cs_pf_coils:
            # Peak field is to be calculated at the Central Solenoid itself,
            # so exclude its own contribution; its self field is
            # dealt with externally using routine BFMAX
            kk = 0
        else:
            # Check different times for maximum current
            if (
                abs(
                    pfcoil_variables.c_pf_cs_coil_pulse_start_ma[i - 1]
                    - pfcoil_variables.c_pf_cs_coils_peak_ma[i - 1]
                )
                < 1.0e-12
            ):
                it = 2
            elif (
                abs(
                    pfcoil_variables.c_pf_cs_coil_flat_top_ma[i - 1]
                    - pfcoil_variables.c_pf_cs_coils_peak_ma[i - 1]
                )
                < 1.0e-12
            ):
                it = 4
            elif (
                abs(
                    pfcoil_variables.c_pf_cs_coil_pulse_end_ma[i - 1]
                    - pfcoil_variables.c_pf_cs_coils_peak_ma[i - 1]
                )
                < 1.0e-12
            ):
                it = 5
            else:
                raise ProcessValueError(
                    "Illegal value of it; possible rounding error", it=it
                )

            if bv.iohcl == 0:
                # No Central Solenoid
                kk = 0
            else:
                sgn = (
                    1.0
                    if pfcoil_variables.j_cs_pulse_start
                    > pfcoil_variables.j_cs_flat_top_end
                    else -1.0
                )

                # Current in each filament representing part of the Central Solenoid
                for iohc in range(pfcoil_variables.nfxf):
                    pfcoil_variables.cfxf[iohc] = (
                        pfcoil_variables.waves[
                            pfcoil_variables.n_cs_pf_coils - 1, it - 1
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
                if iii == ii - 1:
                    # Self field from coil (Lyle's Method)
                    kk = kk + 1

                    dzpf = (
                        pfcoil_variables.z_pf_coil_upper[jj - 1]
                        - pfcoil_variables.z_pf_coil_lower[jj - 1]
                    )
                    pfcoil_variables.rfxf[kk - 1] = pfcoil_variables.r_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.zfxf[kk - 1] = (
                        pfcoil_variables.z_pf_coil_middle[jj - 1] + dzpf * 0.125e0
                    )
                    pfcoil_variables.cfxf[kk - 1] = (
                        pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfcoil_variables.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pfcoil_variables.rfxf[kk - 1] = pfcoil_variables.r_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.zfxf[kk - 1] = (
                        pfcoil_variables.z_pf_coil_middle[jj - 1] + dzpf * 0.375e0
                    )
                    pfcoil_variables.cfxf[kk - 1] = (
                        pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfcoil_variables.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pfcoil_variables.rfxf[kk - 1] = pfcoil_variables.r_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.zfxf[kk - 1] = (
                        pfcoil_variables.z_pf_coil_middle[jj - 1] - dzpf * 0.125e0
                    )
                    pfcoil_variables.cfxf[kk - 1] = (
                        pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfcoil_variables.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pfcoil_variables.rfxf[kk - 1] = pfcoil_variables.r_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.zfxf[kk - 1] = (
                        pfcoil_variables.z_pf_coil_middle[jj - 1] - dzpf * 0.375e0
                    )
                    pfcoil_variables.cfxf[kk - 1] = (
                        pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfcoil_variables.waves[jj - 1, it - 1]
                        * 0.25e6
                    )

                else:
                    # Field from different coil
                    kk = kk + 1
                    pfcoil_variables.rfxf[kk - 1] = pfcoil_variables.r_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.zfxf[kk - 1] = pfcoil_variables.z_pf_coil_middle[
                        jj - 1
                    ]
                    pfcoil_variables.cfxf[kk - 1] = (
                        pfcoil_variables.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfcoil_variables.waves[jj - 1, it - 1]
                        * 1.0e6
                    )

        # Plasma contribution
        if it > 2:
            kk = kk + 1
            pfcoil_variables.rfxf[kk - 1] = pv.rmajor
            pfcoil_variables.zfxf[kk - 1] = 0.0e0
            pfcoil_variables.cfxf[kk - 1] = pv.plasma_current

        # Calculate the field at the inner and outer edges
        # of the coil of interest
        pfcoil_variables.xind[:kk], bri, bzi, psi = bfield(
            pfcoil_variables.rfxf[:kk],
            pfcoil_variables.zfxf[:kk],
            pfcoil_variables.cfxf[:kk],
            pfcoil_variables.r_pf_coil_inner[i - 1],
            pfcoil_variables.z_pf_coil_middle[i - 1],
        )
        pfcoil_variables.xind[:kk], bro, bzo, psi = bfield(
            pfcoil_variables.rfxf[:kk],
            pfcoil_variables.zfxf[:kk],
            pfcoil_variables.cfxf[:kk],
            pfcoil_variables.r_pf_coil_outer[i - 1],
            pfcoil_variables.z_pf_coil_middle[i - 1],
        )

        # b_pf_coil_peak and bpf2 for the Central Solenoid are calculated in OHCALC
        if (bv.iohcl != 0) and (i == pfcoil_variables.n_cs_pf_coils):
            return bri, bro, bzi, bzo

        bpfin = math.sqrt(bri**2 + bzi**2)
        bpfout = math.sqrt(bro**2 + bzo**2)
        for n in range(pfcoil_variables.n_pf_coils_in_group[ii - 1]):
            pfcoil_variables.b_pf_coil_peak[i - 1 + n] = bpfin
            pfcoil_variables.bpf2[i - 1 + n] = bpfout

        return bri, bro, bzi, bzo

    def bfmax(self, rj, a, b, h):
        """Calculates the maximum field of a solenoid.

        author: P J Knight, CCFE, Culham Science Centre
        This routine calculates the peak field (T) at a solenoid's
        inner radius, using fits taken from the figure
        on p.22 of M. Wilson's book Superconducting Magnets,
        Clarendon Press, Oxford, N.Y., 1983

        :param rj: overall current density (A/m2)
        :type rj: float
        :param a: solenoid inner radius (m)
        :type a: float
        :param b: solenoid outer radius (m)
        :type b: float
        :param h: solenoid half height (m)
        :type h: float
        :return bfmax: maximum field of solenoid
        :rtype: float
        """
        beta = h / a
        alpha = b / a

        # Fits are for 1 < alpha < 2 , and 0.5 < beta < very large
        b0 = (
            rj
            * constants.RMU0
            * h
            * math.log(
                (alpha + math.sqrt(alpha**2 + beta**2))
                / (1.0 + math.sqrt(1.0 + beta**2))
            )
        )

        if beta > 3.0:
            b1 = constants.RMU0 * rj * (b - a)
            f = (3.0 / beta) ** 2
            bfmax = f * b0 * (1.007 + (alpha - 1.0) * 0.0055) + (1.0 - f) * b1

        elif beta > 2.0:
            rat = (1.025 - (beta - 2.0) * 0.018) + (alpha - 1.0) * (
                0.01 - (beta - 2.0) * 0.0045
            )
            bfmax = rat * b0

        elif beta > 1.0:
            rat = (1.117 - (beta - 1.0) * 0.092) + (alpha - 1.0) * (beta - 1.0) * 0.01
            bfmax = rat * b0

        elif beta > 0.75:
            rat = (1.30 - 0.732 * (beta - 0.75)) + (alpha - 1.0) * (
                0.2 * (beta - 0.75) - 0.05
            )
            bfmax = rat * b0

        else:
            rat = (1.65 - 1.4 * (beta - 0.5)) + (alpha - 1.0) * (
                0.6 * (beta - 0.5) - 0.20
            )
            bfmax = rat * b0

        return bfmax

    def vsec(self):
        """Calculation of volt-second capability of PF system.

        author: P J Knight, CCFE, Culham Science Centre
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

    def hoop_stress(self, r):
        """Calculation of hoop stress of central solenoid.

        author: J Morris, CCFE, Culham Science Centre
        This routine calculates the hoop stress of the central solenoid
        from "Superconducting magnets", M. N. Wilson OUP

        :param r: radial position a < r < b
        :type r: float
        :return: hoop stress (MPa)
        :rtype: float
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

        return s_hoop_nom / pfcoil_variables.f_a_cs_steel

    def axial_stress(self):
        """Calculation of axial stress of central solenoid.

        author: J Morris, CCFE, Culham Science Centre
        This routine calculates the axial stress of the central solenoid
        from "Case studies in superconducting magnets", Y. Iwasa, Springer

        :return: unsmeared axial stress [MPa], axial force [N]
        :rtype: tuple[float, float]
        """
        b = pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]

        # Half height of central Solenoid [m]
        hl = pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]

        # Central Solenoid current [A]
        ni = (
            pfcoil_variables.c_pf_cs_coils_peak_ma[pfcoil_variables.n_cs_pf_coils - 1]
            * 1.0e6
        )

        # kb term for elliptical integrals
        # kb2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + hl**2))
        kb2 = (4.0e0 * b**2) / (4.0e0 * b**2 + hl**2)

        # k2b term for elliptical integrals
        # k2b2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + 4.0e0*hl**2))
        k2b2 = (4.0e0 * b**2) / (4.0e0 * b**2 + 4.0e0 * hl**2)

        # term 1
        axial_term_1 = -(constants.RMU0 / 2.0e0) * (ni / (2.0e0 * hl)) ** 2

        # term 2
        ekb2_1 = ellipk(kb2)
        ekb2_2 = ellipe(kb2)
        axial_term_2 = (
            2.0e0 * hl * (math.sqrt(4.0e0 * b**2 + hl**2)) * (ekb2_1 - ekb2_2)
        )

        # term 3
        ek2b2_1 = ellipk(k2b2)
        ek2b2_2 = ellipe(k2b2)
        axial_term_3 = (
            2.0e0 * hl * (math.sqrt(4.0e0 * b**2 + 4.0e0 * hl**2)) * (ek2b2_1 - ek2b2_2)
        )

        # calculate axial force [N]
        axial_force = axial_term_1 * (axial_term_2 - axial_term_3)

        # axial area [m2]
        area_ax = constants.PI * (
            pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1] ** 2
            - pfcoil_variables.r_pf_coil_inner[pfcoil_variables.n_cs_pf_coils - 1] ** 2
        )

        # calculate unsmeared axial stress [MPa]
        s_axial = axial_force / (pfcoil_variables.f_a_cs_steel * 0.5 * area_ax)

        return s_axial, axial_force

    def induct(self, output):
        """Calculates PF coil set mutual inductance matrix.

        author: P J Knight, CCFE, Culham Science Centre
        This routine calculates the mutual inductances between all the
        PF coils.

        :param output: switch for writing to output file
        :type output: bool
        """
        nohmax = 200
        nplas = 1

        br = 0.0
        bz = 0.0
        psi = 0.0
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

        noh = int(
            math.ceil(
                2.0e0
                * pfcoil_variables.z_pf_coil_upper[pfcoil_variables.n_cs_pf_coils - 1]
                / (
                    pfcoil_variables.r_pf_coil_outer[pfcoil_variables.n_cs_pf_coils - 1]
                    - pfcoil_variables.r_pf_coil_inner[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                )
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

                xcin, br, bz, psi = bfield(rc, zc, cc, reqv - deltar, zp)
                xcout, br, bz, psi = bfield(rc, zc, cc, reqv + deltar, zp)

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
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
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
                xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
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
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
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
                    ) / math.sqrt(constants.PI)
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
            "(ld_ratio_cst)",
            pfcoil_variables.ld_ratio_cst,
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
            "(t_structural_radial)",
            csfv.t_structural_radial,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Vertical thickness of steel conduit to cable space [m]",
            "(t_structural_vertical)",
            csfv.t_structural_vertical,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Corner radius of CS turn [m]",
            "(r_out_cst)",
            pfcoil_variables.r_out_cst,
            "OP ",
        )

    def outpf(self):
        """Routine to write output from PF coil module to file.

        author: P J Knight, CCFE, Culham Science Centre
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
                    "(f_a_cs_steel)",
                    pfcoil_variables.f_a_cs_steel,
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
                    "(sig_axial)",
                    pfcoil_variables.sig_axial,
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
                    "(axial_force)",
                    pfcoil_variables.axial_force,
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
                        "(t_structural_vertical)",
                        csfv.t_structural_vertical,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS structural radial thickness (m)",
                        "(t_structural_radial)",
                        csfv.t_structural_radial,
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
                # iteration variable (39) fjohc0
                # iteration variable(38) fjohc
                if (
                    abs(pfcoil_variables.j_cs_flat_top_end)
                    > 0.99e0
                    * abs(
                        numerics.boundu[37]
                        * pfcoil_variables.j_cs_critical_flat_top_end
                    )
                ) or (
                    abs(pfcoil_variables.j_cs_pulse_start)
                    > 0.99e0
                    * abs(
                        numerics.boundu[38] * pfcoil_variables.j_cs_critical_pulse_start
                    )
                ):
                    pfcoil_variables.cslimit = True
                if (
                    pfcoil_variables.temp_cs_superconductor_margin
                    < 1.01e0 * tfv.temp_cs_superconductor_margin_min
                ):
                    pfcoil_variables.cslimit = True
                if not pfcoil_variables.cslimit:
                    logger.warning(
                        "CS not using max current density: further optimisation may be possible"
                    )

                # Check whether CS coil currents are feasible from engineering POV
                if ctv.fjohc > 0.7:
                    logger.error(
                        "fjohc shouldn't be above 0.7 for engineering reliability"
                    )
                if ctv.fjohc0 > 0.7:
                    logger.error(
                        "fjohc0 shouldn't be above 0.7 for engineering reliability"
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
            "\t(MA)\t\t(A/m2)\t\t(A/m2)\t\tratio\t\t(kg)\t\t(kg)\t\t(T)",
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

        author: P J Knight, CCFE, Culham Science Centre
        author: R Kemp, CCFE, Culham Science Centre
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
            line += f"\t\t{tv.tim[k]:.2f}"
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

        author: M. Kovari, CCFE
        This routine calculates the self inductance in Henries
        Radiotron Designers Handbook (4th Edition) chapter 10

        :param a: mean radius of coil (m)
        :type a: float
        :param b: length of coil (m) (given as l in the reference)
        :type b: float
        :param c: radial winding thickness (m)
        :type c: float
        :param N: number of turns
        :type N: float
        :return selfinductance: the self inductance in Henries
        :rtype: float
        """
        return (
            (1.0e-6 / 0.0254e0)
            * a**2
            * n**2
            / (9.0e0 * a + 10.0e0 * b + 8.4e0 * c + 3.2e0 * c * b / a)
        )

    def waveform(self):
        """Sets up the PF coil waveforms.

        author: P J Knight, CCFE, Culham Science Centre
        This routine sets up the PF coil current waveforms.
        waves[i,j] is the current in coil i, at time j,
        normalized to the peak current in that coil at any time.
        """
        nplas = pfcoil_variables.n_cs_pf_coils + 1
        for it in range(6):
            pfcoil_variables.waves[nplas - 1, it] = 1.0e0

        for ic in range(pfcoil_variables.n_cs_pf_coils):
            # Find where the peak current occurs
            # Beginning of pulse, t = t_precharge
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

            # Beginning of flat-top, t = t_precharge + t_current_ramp_up
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

            # End of flat-top, t = t_precharge + t_current_ramp_up + t_fusion_ramp + t_burn
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
            pfcoil_variables.waves[ic, 0] = 0.0e0
            pfcoil_variables.waves[ic, 1] = (
                pfcoil_variables.c_pf_cs_coil_pulse_start_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.waves[ic, 2] = (
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.waves[ic, 3] = (
                pfcoil_variables.c_pf_cs_coil_flat_top_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.waves[ic, 4] = (
                pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            pfcoil_variables.waves[ic, 5] = 0.0e0

    def superconpf(
        self, bmax, fhe, fcu, jwp, isumat, fhts, strain, thelium, bcritsc, tcritsc
    ):
        """Routine to calculate the PF coil superconductor properties.

        This routine calculates the superconductor critical winding pack
        current density for the PF coils, plus the temperature margin.
        It is based on the TF coil version, supercon.
        author: P J Knight, CCFE, Culham Science Centre

        N.B. critical current density for a super conductor (j_crit_sc)
        is for the superconducting strands/tape, not including copper.
        Critical current density for a cable (j_crit_cable) acounts for
        both the fraction of the cable taken up by helium coolant channels,
        and the cable conductor copper fraction - i.e., the copper in the
        superconducting strands AND any addtional copper, such as REBCO
        tape support.

        :param bmax: peak field at conductor (T)
        :type bmax: float
        :param fhe: fraction of cable space that is for He cooling
        :type fhe: float
        :param fcu: fraction of cable conductor that is copper
        :type fcu: float
        :param jwp: actual winding pack current density (A/m2)
        :type jwp: float
        :param isumat: switch for conductor type
        1 = ITER Nb3Sn, standard parameters,
        2 = Bi-2212 High Temperature Superconductor,
        3 = NbTi,
        4 = ITER Nb3Sn, user-defined parameters
        5 = WST Nb3Sn parameterisation
        7 = Durham Ginzbug-Landau Nb-Ti parameterisation
        :type isumat: int
        :param fhts: Adjustment factor (<= 1) to account for strain,
        radiation damage, fatigue or AC losses
        :type fhts: float
        :param strain: Strain on superconductor at operation conditions
        :type strain: float
        :param thelium: He temperature at peak field point (K)
        :type thelium: float
        :param bcritsc: Critical field at zero temperature and strain (T) (isumat=4 only)
        :type bcritsc: float
        :param tcritsc: Critical temperature at zero field and strain (K) (isumat=4 only)
        :type tcritsc: float
        :return: Critical winding pack current density (A/m2) (j_crit_wp),
        Critical cable current density (A/m2) (j_crit_cable)
        Superconducting strand non-copper critical current density (A/m2) (j_crit_sc)
        Temperature margin (K) (tmarg)
        :rtype: tuple[float, float, float, float]
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

            # The CS coil current at EOF
            # ioheof = bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfcoil_variables.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # Copper area = (area of coil - area of steel)*(1- void fraction)*
            # (fraction of copper in strands)
            # rcv.copperaoh_m2 = ioheof / (pfcoil_variables.awpoh * (1.0 - pfcoil_variables.f_a_cs_void) * pfcoil_variables.fcuohsu)

        elif isumat == 7:
            # Durham Ginzburg-Landau critical surface model for Nb-Ti
            bc20m = tfv.b_crit_upper_nbti
            tc0m = tfv.t_crit_nbti
            j_crit_sc, _, _ = superconductors.gl_nbti(
                thelium, bmax, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

            # The CS coil current at EOF
            # ioheof = bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfcoil_variables.j_cs_flat_top_end

        elif isumat == 8:
            # Durham Ginzburg-Landau critical surface model for REBCO
            bc20m = 429e0
            tc0m = 185e0
            j_crit_sc, _, _ = superconductors.gl_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

            # The CS coil current at EOF
            # ioheof = bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfcoil_variables.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # rcv.copperaoh_m2 = ioheof / (pfcoil_variables.awpoh * (1.0 - pfcoil_variables.f_a_cs_void) * pfcoil_variables.fcuohsu)

        elif isumat == 9:
            # Hazelton experimental data + Zhai conceptual model for REBCO
            bc20m = 138
            tc0m = 92
            j_crit_sc, _, _ = superconductors.hijc_rebco(
                thelium,
                bmax,
                bc20m,
                tc0m,
                rcv.tape_width,
                rcv.rebco_thickness,
                rcv.tape_thickness,
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

            # The CS coil current at EOF
            # ioheof = bv.z_tf_inside_half * pfcoil_variables.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfcoil_variables.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # rcv.copperaoh_m2 = ioheof / (pfcoil_variables.awpoh * (1.0 - pfcoil_variables.f_a_cs_void) * pfcoil_variables.fcuohsu)

        else:
            # Error condition
            raise ProcessValueError(
                "Illegal value for i_pf_superconductor", isumat=isumat
            )

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
            t_zero_margin, root_result = optimize.newton(
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
def bfield(rc, zc, cc, rp, zp):
    """Calculate the field at a point due to currents in a number
    of circular poloidal conductor loops.
    author: P J Knight, CCFE, Culham Science Centre
    author: D Strickler, ORNL
    author: J Galambos, ORNL
    nc : input integer : number of loops
    rc(nc) : input real array : R coordinates of loops (m)
    zc(nc) : input real array : Z coordinates of loops (m)
    cc(nc) : input real array : Currents in loops (A)
    xc(nc) : output real array : Mutual inductances (H)
    rp, zp : input real : coordinates of point of interest (m)
    br : output real : radial field component at (rp,zp) (T)
    bz : output real : vertical field component at (rp,zp) (T)
    psi : output real : poloidal flux at (rp,zp) (Wb)
    This routine calculates the magnetic field components and
    the poloidal flux at an (R,Z) point, given the locations
    and currents of a set of conductor loops.
    <P>The mutual inductances between the loops and a poloidal
    filament at the (R,Z) point of interest is also found."""

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

    nc = len(rc)

    xc = np.empty((nc,))
    br = 0
    bz = 0
    psi = 0

    for i in range(nc):
        d = (rp + rc[i]) ** 2 + (zp - zc[i]) ** 2
        s = 4.0 * rp * rc[i] / d

        # Kludge: avoid s >= 1.0, a goes inf
        if s > 0.999999:
            s = 0.999999

        t = 1.0 - s
        a = np.log(1.0 / t)

        dz = zp - zc[i]
        zs = dz**2
        dr = rp - rc[i]
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

        xc[i] = 0.5 * constants.RMU0 * sd * ((2.0 - s) * xk - 2.0 * xe)

        #  Radial, vertical fields

        brx = (
            constants.RMU0
            * cc[i]
            * dz
            / (2 * np.pi * rp * sd)
            * (-xk + (rc[i] ** 2 + rp**2 + zs) / (dr**2 + zs) * xe)
        )
        bzx = (
            constants.RMU0
            * cc[i]
            / (2 * np.pi * sd)
            * (xk + (rc[i] ** 2 - rp**2 - zs) / (dr**2 + zs) * xe)
        )

        #  Sum fields, flux

        br += brx
        bz += bzx
        psi += xc[i] * cc[i]

    return xc, br, bz, psi


@numba.njit(cache=True)
def rsid(npts, brin, bzin, nfix, n_pf_coil_groups, ccls, bfix, gmat):
    """Computes the norm of the residual vectors.

    author: P J Knight, CCFE, Culham Science Centre
    author: D Strickler, ORNL
    author: J Galambos, ORNL
    author: P C Shipe, ORNL
    This routine calculates the residuals from the matrix
    equation for calculation of the currents in a group of ring coils.

    :param npts: number of data points at which field is  to be fixed;
    should be <= nptsmx
    :type npts: int
    :param brin: field components at data points (T)
    :type brin: numpy.ndarray
    :param bzin: field components at data points (T)
    :type bzin: numpy.ndarray
    :param nfix: number of coils with fixed currents, <= nfixmx
    :type nfix: int
    :param n_pf_coil_groups: number of coil groups, where all coils in a group have the
    same current, <= n_pf_groups_max
    :type n_pf_coil_groups: int
    :param ccls: coil currents in each group (A)
    :type ccls: numpy.ndarray
    :param bfix: work array
    :type bfix: numpy.ndarray
    :param gmat: work array
    :type gmat: numpy.ndarray
    :return: sum of squares of radial field residues (brssq), radial field
    residue norm (brnrm), sum of squares of vertical field residues (bzssq),
    vertical field residue norm (bznrm), sum of squares of elements of
    residual vector (ssq)
    :rtype: tuple[float, float, float, float, float]
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

    author: P J Knight, CCFE, Culham Science Centre
    author: D Strickler, ORNL
    author: J Galambos, ORNL
    This routine calculates the fields at the points specified by
    (rpts,zpts) from the set of coils with fixed currents.

    :param lrow1: row length of array bfix; should be >= nptsmx
    :type lrow1: int
    :param npts: number of data points at which field is to be fixed;
    should be <= nptsmx
    :type npts: int
    :param rpts: coords of data points (m)
    :type rpts: numpy.ndarray
    :param zpts: coords of data points (m)
    :type zpts: numpy.ndarray
    :param nfix: number of coils with fixed currents, <= nfixmx
    :type nfix: int
    :param rfix: coordinates of coils with fixed currents (m)
    :type rfix: numpy.ndarray
    :param zfix: coordinates of coils with fixed currents (m)
    :type zfix: numpy.ndarray
    :param cfix: Fixed currents (A)
    :type cfix: numpy.ndarray
    :return: Fields at data points (T)
    :rtype: numpy.ndarray
    """
    bfix = np.zeros(lrow1)

    if nfix <= 0:
        return bfix

    for i in range(npts):
        # bfield() only operates correctly on nfix slices of array
        # arguments, not entire arrays
        _, brw, bzw, _ = bfield(rfix[:nfix], zfix[:nfix], cfix[:nfix], rpts[i], zpts[i])
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
    author: P J Knight, CCFE, Culham Science Centre
    author: D Strickler, ORNL
    author: J Galambos, ORNL

    :param lrow1: row length of arrays bfix, bvec, gmat, umat, vmat; should
    be >= (2*nptsmx + n_pf_groups_max)
    :type lrow1: int
    :param lcol1: column length of arrays gmat, umat, vmat; should be >=
    n_pf_groups_max
    :type lcol1: int
    :param npts: number of data points at which field is to be fixed; should
    be <= nptsmx
    :type npts: int
    :param rpts: coords of data points (m)
    :type rpts: numpy.ndarray
    :param zpts: coords of data points (m)
    :type zpts: numpy.ndarray
    :param brin: field components at data points (T)
    :type brin: numpy.ndarray
    :param bzin: field components at data points (T)
    :type bzin: numpy.ndarray
    :param n_pf_coil_groups: number of coil groups, where all coils in a group have the
    same current, <= n_pf_groups_max
    :type n_pf_coil_groups: int
    :param n_pf_coils_in_group: number of coils in each group, each value <= n_pf_coils_in_group_max
    :type n_pf_coils_in_group: numpy.ndarray
    :param r_pf_coil_middle_group_array: coords R(i,j), Z(i,j) of coil j in group i (m)
    :type r_pf_coil_middle_group_array: numpy.ndarray
    :param z_pf_coil_middle_group_array: coords R(i,j), Z(i,j) of coil j in group i (m)
    :type z_pf_coil_middle_group_array: numpy.ndarray
    :param alfa: smoothing parameter (0 = no smoothing, 1.0D-9 = large
    smoothing)
    :type alfa: float
    :param bfix: Fields at data points (T)
    :type bfix: numpy.ndarray
    :return: actual number of rows to use, work array, work array,
    Coordinates of conductor loops (m), Coordinates of conductor loops (m),
    Currents in conductor loops (A), Mutual inductances (H)
    :rtype: tuple[int, numpy.ndarray, numpy.ndarray, numpy.ndarray
    numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    bvec = np.zeros(lrow1)
    gmat = np.zeros((lrow1, lcol1))
    cc = np.ones(n_pf_coils_in_group_max)

    for i in range(npts):
        bvec[i] = brin[i] - bfix[i]
        bvec[i + npts] = bzin[i] - bfix[i + npts]

        for j in range(n_pf_coil_groups):
            nc = n_pf_coils_in_group[j]

            _, gmat[i, j], gmat[i + npts, j], _ = bfield(
                r_pf_coil_middle_group_array[j, :nc],
                z_pf_coil_middle_group_array[j, :nc],
                cc[:nc],
                rpts[i],
                zpts[i],
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
