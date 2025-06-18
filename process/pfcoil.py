import logging
import math

import numba
import numpy as np
from scipy import optimize
from scipy.linalg import svd
from scipy.special import ellipe, ellipk

import process.superconductors as superconductors
from process import fortran as ft
from process import process_output as op
from process.exceptions import ProcessValueError
from process.fortran import build_variables as bv
from process.fortran import constants, numerics, rebco_variables
from process.fortran import constraint_variables as ctv
from process.fortran import cs_fatigue_variables as csfv
from process.fortran import error_handling as eh
from process.fortran import fwbs_variables as fwbsv
from process.fortran import pfcoil_module as pf
from process.fortran import pfcoil_variables as pfv
from process.fortran import physics_variables as pv
from process.fortran import rebco_variables as rcv
from process.fortran import tfcoil_variables as tfv
from process.fortran import times_variables as tv
from process.utilities.f2py_string_patch import f2py_compatible_to_string

logger = logging.getLogger(__name__)

RMU0 = ft.constants.rmu0


class PFCoil:
    """Calculate poloidal field coil system parameters."""

    def __init__(self, cs_fatigue) -> None:
        """Initialise Fortran module variables."""
        self.outfile = ft.constants.nout  # output file unit
        self.mfile = ft.constants.mfile  # mfile file unit
        init_pfcoil_module()
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
        lrow1 = 2 * pfv.nptsmx + pfv.n_pf_groups_max
        lcol1 = pfv.n_pf_groups_max

        pcls0 = np.zeros(pfv.n_pf_groups_max, dtype=int)
        ncls0 = np.zeros(pfv.n_pf_groups_max + 2, dtype=int)

        pf.rcls0, pf.zcls0 = np.zeros(
            (2, pfv.n_pf_groups_max, pfv.n_pf_coils_in_group_max), order="F"
        )
        pf.ccls0 = np.zeros(int(pfv.n_pf_groups_max / 2))
        sigma, work2 = np.zeros((2, pfv.n_pf_groups_max))
        rc, zc, cc, xc = np.zeros((4, pfv.n_pf_coils_in_group_max))
        brin, bzin, rpts, zpts = np.zeros((4, pfv.nptsmx))
        bfix, bvec = np.zeros((2, lrow1))
        gmat, umat, vmat = np.zeros((3, lrow1, lcol1), order="F")
        signn = np.zeros(2)
        aturn = np.zeros(pfv.ngc2)

        # Toggle switch for i_pf_location()=2 coils above/below midplane
        top_bottom = 1

        # Set up the number of PF coils including the Central Solenoid (n_cs_pf_coils),
        # and the number of PF circuits including the plasma (n_pf_cs_plasma_circuits)
        if pfv.n_pf_coil_groups > pfv.n_pf_groups_max:
            raise ProcessValueError(
                "n_pf_coil_groups is larger than n_pf_groups_max",
                n_pf_coil_groups=pfv.n_pf_coil_groups,
                n_pf_groups_max=pfv.n_pf_groups_max,
            )

        # Total the number of PF coils in all groups, and check that none
        # exceeds the limit
        pfv.n_cs_pf_coils = 0
        for i in range(pfv.n_pf_coil_groups):
            if pfv.n_pf_coils_in_group[i] > pfv.n_pf_coils_in_group_max:
                raise ProcessValueError(
                    "PFCOIL: Too many coils in a PF coil group",
                    i=i,
                    n_pf_coils_in_group=pfv.n_pf_coils_in_group[i],
                    n_pf_coils_in_group_max=pfv.n_pf_coils_in_group_max,
                )

            pfv.n_cs_pf_coils = pfv.n_cs_pf_coils + pfv.n_pf_coils_in_group[i]

        # Add one if an Central Solenoid is present, and make an extra group
        if bv.iohcl != 0:
            pfv.n_cs_pf_coils = pfv.n_cs_pf_coils + 1
            pfv.n_pf_coils_in_group[pfv.n_pf_coil_groups] = 1

        # Add one for the plasma
        pfv.n_pf_cs_plasma_circuits = pfv.n_cs_pf_coils + 1

        # Overall current density in the Central Solenoid at beginning of pulse
        pfv.j_cs_pulse_start = (
            pfv.j_cs_flat_top_end * pfv.f_j_cs_start_pulse_end_flat_top
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

        # Central Solenoid mean radius
        pfv.r_cs_middle = bv.dr_bore + 0.5e0 * bv.dr_cs

        # nfxf is the total no of filaments into which the Central Solenoid is split,
        # if present
        if bv.iohcl == 0:
            pf.nfxf = 0
            ioheof = 0.0e0
        else:
            pf.nfxf = 2 * pfv.nfxfh

            # total Central Solenoid current at EOF
            ioheof = (
                -bv.hmax
                * pfv.f_z_cs_tf_internal
                * bv.dr_cs
                * 2.0e0
                * pfv.j_cs_flat_top_end
            )

            if pf.nfxf > pfv.nfixmx:
                raise ProcessValueError(
                    "Too many filaments nfxf repesenting the OH coil",
                    nfxf=pf.nfxf,
                    nfixmx=pfv.nfixmx,
                )

            # Symmetric up/down Central Solenoid : Find (R,Z) and current of each filament at BOP

            for nng in range(pfv.nfxfh):
                pf.rfxf[nng] = pfv.r_cs_middle
                pf.rfxf[nng + pfv.nfxfh] = pf.rfxf[nng]
                pf.zfxf[nng] = (
                    bv.hmax * pfv.f_z_cs_tf_internal / pfv.nfxfh * ((nng + 1) - 0.5e0)
                )
                pf.zfxf[nng + pfv.nfxfh] = -pf.zfxf[nng]
                pf.cfxf[nng] = -ioheof / pf.nfxf * pfv.f_j_cs_start_pulse_end_flat_top
                pf.cfxf[nng + pfv.nfxfh] = pf.cfxf[nng]

        # Scale PF coil locations
        signn[0] = 1.0e0
        signn[1] = -1.0e0
        pf.rclsnorm = bv.r_tf_outboard_mid + 0.5e0 * bv.dr_tf_outboard + pfv.routr

        # Place the PF coils:

        # N.B. Problems here if k=n_pf_coils_in_group(group) is greater than 2.
        for j in range(pfv.n_pf_coil_groups):
            if pfv.i_pf_location[j] == 1:
                # PF coil is stacked on top of the Central Solenoid
                for k in range(pfv.n_pf_coils_in_group[j]):
                    pf.rcls[j, k] = pfv.r_cs_middle + pfv.rpf1

                    # Z coordinate of coil enforced so as not
                    # to occupy the same space as the Central Solenoid
                    pf.zcls[j, k] = signn[k] * (
                        bv.hmax * pfv.f_z_cs_tf_internal
                        + 0.1e0
                        + 0.5e0
                        * (
                            bv.hmax * (1.0e0 - pfv.f_z_cs_tf_internal)
                            + bv.dr_tf_inboard
                            + 0.1e0
                        )
                    )

            elif pfv.i_pf_location[j] == 2:
                # PF coil is on top of the TF coil
                for k in range(pfv.n_pf_coils_in_group[j]):
                    pf.rcls[j, k] = pv.rmajor + pfv.rpf2 * pv.triang * pv.rminor
                    if pv.itart == 1 and pv.itartpf == 0:
                        pf.zcls[j, k] = (bv.hmax - pfv.zref[j]) * signn[k]
                    else:
                        # pf.zcls(j,k) = (bv.hmax + bv.dr_tf_inboard + 0.86e0) * signn(k)
                        if top_bottom == 1:  # this coil is above midplane
                            pf.zcls[j, k] = bv.hpfu + 0.86e0
                            top_bottom = -1
                        else:  # this coil is below midplane
                            pf.zcls[j, k] = -1.0e0 * (
                                bv.hpfu - 2.0e0 * bv.hpfdif + 0.86e0
                            )
                            top_bottom = 1

            elif pfv.i_pf_location[j] == 3:
                # PF coil is radially outside the TF coil
                for k in range(pfv.n_pf_coils_in_group[j]):
                    pf.zcls[j, k] = pv.rminor * pfv.zref[j] * signn[k]
                    # Coil radius follows TF coil curve for SC TF (D-shape)
                    # otherwise stacked for resistive TF (rectangle-shape)
                    if tfv.i_tf_sup != 1 or pfv.i_sup_pf_shape == 1:
                        pf.rcls[j, k] = pf.rclsnorm
                    else:
                        pf.rcls[j, k] = math.sqrt(pf.rclsnorm**2 - pf.zcls[j, k] ** 2)
                        try:
                            assert pf.rcls[j, k] < np.inf
                        except AssertionError:
                            logger.exception(
                                "Element of pf.rcls is inf. Kludging to 1e10."
                            )
                            pf.rcls[j, k] = 1e10

            elif pfv.i_pf_location[j] == 4:
                # PF coil is in general location
                # See issue 1418
                # https://git.ccfe.ac.uk/process/process/-/issues/1418
                for k in range(pfv.n_pf_coils_in_group[j]):
                    pf.zcls[j, k] = pv.rminor * pfv.zref[j] * signn[k]
                    pf.rcls[j, k] = pv.rminor * pfv.rref[j] + pv.rmajor

            else:
                raise ProcessValueError(
                    "Illegal i_pf_location value",
                    j=j,
                    i_pf_location=pfv.i_pf_location[j],
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
        if pfv.j_cs_pulse_start != 0.0e0:
            # Find currents for plasma initiation to null field across plasma
            npts = 32  # Number of test points across plasma midplane
            if npts > pfv.nptsmx:
                raise ProcessValueError(
                    "Too many test points npts across plasma midplane",
                    npts=npts,
                    nptsmx=pfv.nptsmx,
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
            pf.ssq0, pf.ccl0 = self.efc(
                npts,
                rpts,
                zpts,
                brin,
                bzin,
                pf.nfxf,
                pf.rfxf,
                pf.zfxf,
                pf.cfxf,
                pfv.n_pf_coil_groups,
                pfv.n_pf_coils_in_group,
                pf.rcls,
                pf.zcls,
                pfv.alfapf,
                bfix,
                gmat,
                bvec,
            )

        # Equilibrium coil currents determined by SVD targeting B
        if pfv.i_pf_current == 1:
            # Simple coil current scaling for STs (good only for A < about 1.8)
            # Bypasses SVD solver
            if pv.itart == 1 and pv.itartpf == 0:
                for i in range(pfv.n_pf_coil_groups):
                    if pfv.i_pf_location[i] == 1:
                        # PF coil is stacked on top of the Central Solenoid
                        pf.ccls[i] = 0.0e0
                        raise ProcessValueError(
                            "i_pf_location(i) should not be 1 if itart=1", i=i
                        )

                    if pfv.i_pf_location[i] == 2:
                        # PF coil is on top of the TF coil
                        pf.ccls[i] = 0.3e0 * pv.aspect**1.6e0 * pv.plasma_current

                    elif pfv.i_pf_location[i] == 3:
                        # PF coil is radially outside the TF coil
                        pf.ccls[i] = -0.4e0 * pv.plasma_current

                    else:
                        raise ProcessValueError(
                            "Illegal value of i_pf_location(i)",
                            i=i,
                            i_pf_location=pfv.i_pf_location[i],
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
                for i in range(pfv.n_pf_coil_groups):
                    if pfv.i_pf_location[i] == 1:
                        # Do not allow if no central solenoid
                        if bv.iohcl == 0:
                            raise ProcessValueError(
                                "i_pf_location(i) should not be 1 if iohcl=0"
                            )
                        # PF coil is stacked on top of the Central Solenoid
                        # This coil is to balance Central Solenoid flux and should not be involved
                        # in equilibrium calculation -- RK 07/12
                        pf.ccls[i] = 0.0e0
                        nfxf0 = nfxf0 + pfv.n_pf_coils_in_group[i]
                        for ccount in range(pfv.n_pf_coils_in_group[i]):
                            pf.rfxf[nocoil] = pf.rcls[i, ccount]
                            pf.zfxf[nocoil] = pf.zcls[i, ccount]
                            pf.cfxf[nocoil] = pf.ccls[i]
                            nocoil = nocoil + 1

                    elif pfv.i_pf_location[i] == 2:
                        # PF coil is on top of the TF coil; divertor coil
                        # This is a fixed current for this calculation -- RK 07/12

                        pf.ccls[i] = (
                            pv.plasma_current
                            * 2.0e0
                            * (1.0e0 - (pv.kappa * pv.rminor) / abs(pf.zcls[i, 0]))
                        )
                        nfxf0 = nfxf0 + pfv.n_pf_coils_in_group[i]
                        for ccount in range(pfv.n_pf_coils_in_group[i]):
                            pf.rfxf[nocoil] = pf.rcls[i, ccount]
                            pf.zfxf[nocoil] = pf.zcls[i, ccount]
                            pf.cfxf[nocoil] = pf.ccls[i]
                            nocoil = nocoil + 1

                    elif pfv.i_pf_location[i] == 3:
                        # PF coil is radially outside the TF coil
                        # This is an equilibrium coil, current must be solved for

                        pcls0[ngrp0] = i + 1
                        ngrp0 = ngrp0 + 1

                    elif pfv.i_pf_location[i] == 4:
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
                            i_pf_location=pfv.i_pf_location[i],
                        )

                for ccount in range(ngrp0):
                    ncls0[ccount] = 2
                    pf.rcls0[ccount, 0] = pf.rcls[pcls0[ccount] - 1, 0]
                    pf.rcls0[ccount, 1] = pf.rcls[pcls0[ccount] - 1, 1]
                    pf.zcls0[ccount, 0] = pf.zcls[pcls0[ccount] - 1, 0]
                    pf.zcls0[ccount, 1] = pf.zcls[pcls0[ccount] - 1, 1]

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

                ssqef, pf.ccls0 = self.efc(
                    npts0,
                    rpts,
                    zpts,
                    brin,
                    bzin,
                    nfxf0,
                    pf.rfxf,
                    pf.zfxf,
                    pf.cfxf,
                    ngrp0,
                    ncls0,
                    pf.rcls0,
                    pf.zcls0,
                    pfv.alfapf,
                    bfix,
                    gmat,
                    bvec,
                )

                for ccount in range(ngrp0):
                    pf.ccls[pcls0[ccount] - 1] = pf.ccls0[ccount]

        # Flux swing from vertical field

        # If this is the first visit to the routine the inductance matrix
        # ind_pf_cs_plasma_mutual and the turns array have not yet been calculated, so we set
        # them to (very) approximate values to avoid strange behaviour...
        if pf.first_call:
            pfv.ind_pf_cs_plasma_mutual[:, :] = 1.0e0
            pfv.n_pf_coil_turns[:] = 100.0e0
            pf.first_call = False

        pfflux = 0.0e0
        nocoil = 0
        for ccount in range(pfv.n_pf_coil_groups):
            for _i in range(pfv.n_pf_coils_in_group[ccount]):
                pfflux = pfflux + (
                    pf.ccls[ccount]
                    * pfv.ind_pf_cs_plasma_mutual[
                        nocoil, pfv.n_pf_cs_plasma_circuits - 1
                    ]
                    / pfv.n_pf_coil_turns[nocoil]
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
                * constants.pi
                * constants.pi
                * (
                    (bv.dr_bore * bv.dr_bore)
                    + (bv.dr_cs * bv.dr_cs) / 6.0e0
                    + (bv.dr_cs * bv.dr_bore) / 2.0e0
                )
                / (bv.hmax * pfv.f_z_cs_tf_internal * 2.0e0)
            )
            dics = csflux / ddics

            pfv.f_j_cs_start_end_flat_top = (
                (-ioheof * pfv.f_j_cs_start_pulse_end_flat_top) + dics
            ) / ioheof
            if np.abs(pfv.f_j_cs_start_end_flat_top) > 1.0:
                logger.warning(
                    "Ratio of central solenoid overall current density at "
                    "beginning of flat-top / end of flat-top > 1 (|f_j_cs_start_end_flat_top| > 1)"
                )
        else:
            dics = 0.0e0
            pfv.f_j_cs_start_end_flat_top = 1.0e0
            eh.report_error(71)

        # Split groups of coils into one set containing ncl coils
        ncl = 0
        for nng in range(pfv.n_pf_coil_groups):
            for ng2 in range(pfv.n_pf_coils_in_group[nng]):
                pfv.r_pf_coil_middle[ncl] = pf.rcls[nng, ng2]
                pfv.z_pf_coil_middle[ncl] = pf.zcls[nng, ng2]

                # Currents at different times:

                # If PF coil currents are computed, not input via ccl0_ma, ccls_ma:
                # Then set ccl0_ma,ccls_ma from the computed pf.ccl0,pf.ccls
                if pfv.i_pf_current != 0:
                    pfv.ccl0_ma[nng] = 1.0e-6 * pf.ccl0[nng]
                    pfv.ccls_ma[nng] = 1.0e-6 * pf.ccls[nng]
                else:
                    # Otherwise set pf.ccl0,pf.ccls via the input ccl0_ma and ccls_ma
                    pf.ccl0[nng] = 1.0e6 * pfv.ccl0_ma[nng]
                    pf.ccls[nng] = 1.0e6 * pfv.ccls_ma[nng]

                # Beginning of pulse: t = tv.t_precharge
                pfv.c_pf_cs_coil_pulse_start_ma[ncl] = 1.0e-6 * pf.ccl0[nng]

                # Beginning of flat-top: t = tv.t_precharge+tv.t_current_ramp_up
                pfv.c_pf_cs_coil_flat_top_ma[ncl] = 1.0e-6 * (
                    pf.ccls[nng]
                    - (
                        pf.ccl0[nng]
                        * pfv.f_j_cs_start_end_flat_top
                        / pfv.f_j_cs_start_pulse_end_flat_top
                    )
                )

                # End of flat-top: t = tv.t_precharge+tv.t_current_ramp_up+tv.t_fusion_ramp+tv.t_burn
                pfv.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * (
                    pf.ccls[nng]
                    - (pf.ccl0[nng] * (1.0e0 / pfv.f_j_cs_start_pulse_end_flat_top))
                )

                ncl = ncl + 1

        # Current in Central Solenoid as a function of time
        # N.B. If the Central Solenoid is not present then ioheof is zero.
        pfv.c_pf_cs_coil_pulse_start_ma[ncl] = (
            -1.0e-6 * ioheof * pfv.f_j_cs_start_pulse_end_flat_top
        )
        pfv.c_pf_cs_coil_flat_top_ma[ncl] = (
            1.0e-6 * ioheof * pfv.f_j_cs_start_end_flat_top
        )
        pfv.c_pf_cs_coil_pulse_end_ma[ncl] = 1.0e-6 * ioheof

        # Set up coil current waveforms, normalised to the peak current in
        # each coil
        self.waveform()  # sets c_pf_cs_coils_peak_ma(), waves()

        # Calculate PF coil geometry, current and number of turns
        # Dimensions are those of the winding pack, and exclude
        # the steel supporting case
        i = 0
        pfv.r_pf_coil_outer_max = 0.0e0

        dz = 0

        for ii in range(pfv.n_pf_coil_groups):
            for _ij in range(pfv.n_pf_coils_in_group[ii]):
                if pfv.i_pf_location[ii] == 1:
                    # PF coil is stacked on top of the Central Solenoid
                    dx = 0.5e0 * bv.dr_cs
                    dz = 0.5e0 * (
                        bv.hmax * (1.0e0 - pfv.f_z_cs_tf_internal)
                        + bv.dr_tf_inboard
                        + 0.1e0
                    )  # ???
                    area = 4.0e0 * dx * dz * pfv.pf_current_safety_factor

                    # Number of turns
                    # c_pf_coil_turn_peak_input[i] is the current per turn (input)
                    pfv.n_pf_coil_turns[i] = abs(
                        (pfv.c_pf_cs_coils_peak_ma[i] * 1.0e6)
                        / pfv.c_pf_coil_turn_peak_input[i]
                    )
                    aturn[i] = area / pfv.n_pf_coil_turns[i]

                    # Actual winding pack current density
                    pfv.j_pf_coil_wp_peak[i] = (
                        1.0e6 * abs(pfv.c_pf_cs_coils_peak_ma[i]) / area
                    )

                    # Location of edges of each coil:
                    # r_pf_coil_inner = inner radius, r_pf_coil_outer = outer radius
                    # z_pf_coil_lower = 'lower' edge z (i.e. edge nearer to midplane)
                    # z_pf_coil_upper = 'upper' edge z (i.e. edge further from midplane)
                    pfv.r_pf_coil_inner[i] = pfv.r_pf_coil_middle[i] - dx
                    pfv.r_pf_coil_outer[i] = pfv.r_pf_coil_middle[i] + dx

                    pfv.z_pf_coil_lower[i] = pfv.z_pf_coil_middle[i] - dz
                    if pfv.z_pf_coil_middle[i] < 0.0e0:
                        pfv.z_pf_coil_lower[i] = pfv.z_pf_coil_middle[i] + dz

                    pfv.z_pf_coil_upper[i] = pfv.z_pf_coil_middle[i] + dz

                    if pfv.z_pf_coil_middle[i] < 0.0e0:
                        pfv.z_pf_coil_upper[i] = pfv.z_pf_coil_middle[i] - dz

                else:
                    # Other coils. N.B. Current density j_pf_coil_wp_peak[i] is defined in
                    # routine INITIAL for these coils.
                    area = (
                        abs(
                            pfv.c_pf_cs_coils_peak_ma[i]
                            * 1.0e6
                            / pfv.j_pf_coil_wp_peak[i]
                        )
                        * pfv.pf_current_safety_factor
                    )

                    pfv.n_pf_coil_turns[i] = abs(
                        (pfv.c_pf_cs_coils_peak_ma[i] * 1.0e6)
                        / pfv.c_pf_coil_turn_peak_input[i]
                    )
                    aturn[i] = area / pfv.n_pf_coil_turns[i]

                    dx = 0.5e0 * math.sqrt(area)  # square cross-section

                    pfv.r_pf_coil_inner[i] = pfv.r_pf_coil_middle[i] - dx
                    pfv.r_pf_coil_outer[i] = pfv.r_pf_coil_middle[i] + dx

                    pfv.z_pf_coil_lower[i] = pfv.z_pf_coil_middle[i] - dx
                    if pfv.z_pf_coil_middle[i] < 0.0e0:
                        pfv.z_pf_coil_lower[i] = pfv.z_pf_coil_middle[i] + dx

                    pfv.z_pf_coil_upper[i] = pfv.z_pf_coil_middle[i] + dx
                    if pfv.z_pf_coil_middle[i] < 0.0e0:
                        pfv.z_pf_coil_upper[i] = pfv.z_pf_coil_middle[i] - dx

                # Outside radius of largest PF coil (m)
                pfv.r_pf_coil_outer_max = max(
                    pfv.r_pf_coil_outer_max, pfv.r_pf_coil_outer[i]
                )

                i = i + 1

        # Calculate peak field, allowable current density, resistive
        # power losses and volumes and weights for each PF coil, index i
        i = 0
        it = 0
        pfv.p_pf_coil_resistive_total_flat_top = 0.0e0
        pfv.m_pf_coil_max = 0.0e0

        for ii in range(pfv.n_pf_coil_groups):
            iii = ii
            for ij in range(pfv.n_pf_coils_in_group[ii]):
                # Peak field

                if ij == 0:
                    # Index args +1ed
                    bri, bro, bzi, bzo = self.peakb(
                        i + 1, iii + 1, it
                    )  # returns b_pf_coil_peak, bpf2

                # Issue 1871.  MDK
                # Allowable current density (for superconducting coils) for each coil, index i
                if pfv.i_pf_conductor == 0:
                    bmax = max(abs(pfv.b_pf_coil_peak[i]), abs(pf.bpf2[i]))

                    pfv.j_pf_wp_critical[i], jstrand, jsc, tmarg = self.superconpf(
                        bmax,
                        pfv.f_a_pf_coil_void[i],
                        pfv.fcupfsu,
                        pfv.j_pf_coil_wp_peak[i],
                        pfv.i_pf_superconductor,
                        tfv.fhts,
                        tfv.str_pf_con_res,
                        tfv.tftmp,
                        tfv.bcritsc,
                        tfv.tcritsc,
                    )

                    # Strand critical current calculation for costing in $/kAm
                    # = superconducting filaments jc * (1 - strand copper fraction)
                    if pfv.i_cs_superconductor.item() in {2, 6, 8}:
                        pfv.j_crit_str_pf = jsc
                    else:
                        pfv.j_crit_str_pf = jsc * (1 - pfv.fcupfsu)

                # Length of conductor

                rll = (
                    2.0e0
                    * constants.pi
                    * pfv.r_pf_coil_middle[i]
                    * pfv.n_pf_coil_turns[i]
                )

                # Resistive coils

                if pfv.i_pf_conductor == 1:
                    # Coil resistance (f_a_pf_coil_void is the void fraction)

                    respf = (
                        pfv.rho_pf_coil
                        * rll
                        / (aturn[i] * (1.0e0 - pfv.f_a_pf_coil_void[i]))
                    )

                    # Sum resistive power losses

                    pfv.p_pf_coil_resistive_total_flat_top = (
                        pfv.p_pf_coil_resistive_total_flat_top
                        + respf
                        * (
                            1.0e6
                            * pfv.c_pf_cs_coil_pulse_start_ma[i]
                            / pfv.n_pf_coil_turns[i]
                        )
                        ** 2
                    )

                # Winding pack volume

                volpf = aturn[i] * rll

                # Conductor weight (f_a_pf_coil_void is the void fraction)

                if pfv.i_pf_conductor == 0:
                    pfv.m_pf_coil_conductor[i] = (
                        volpf
                        * tfv.dcond[pfv.i_pf_superconductor - 1]
                        * (1.0e0 - pfv.f_a_pf_coil_void[i])
                    )
                else:
                    pfv.m_pf_coil_conductor[i] = (
                        volpf * constants.dcopper * (1.0e0 - pfv.f_a_pf_coil_void[i])
                    )

                # (J x B) force on coil

                forcepf = (
                    0.5e6
                    * (pfv.b_pf_coil_peak[i] + pf.bpf2[i])
                    * abs(pfv.c_pf_cs_coils_peak_ma[i])
                    * pfv.r_pf_coil_middle[i]
                )

                # Stress ==> cross-sectional area of supporting steel to use

                if pfv.i_pf_conductor == 0:
                    # Superconducting coil
                    # Updated assumptions: 500 MPa stress limit with all of the force
                    # supported in the conduit (steel) case.
                    # Now, 500 MPa replaced by sigpfcalw, sigpfcf now defaultly set to 1

                    areaspf = pfv.sigpfcf * forcepf / (pfv.sigpfcalw * 1.0e6)

                    # Thickness of hypothetical steel casing assumed to encase the PF
                    # winding pack; in reality, the steel is distributed
                    # throughout the conductor. Issue #152
                    # Assume a case of uniform thickness around coil cross-section
                    # Thickness found via a simple quadratic equation

                    drpdz = (
                        pfv.r_pf_coil_outer[i]
                        - pfv.r_pf_coil_inner[i]
                        + abs(pfv.z_pf_coil_upper[i] - pfv.z_pf_coil_lower[i])
                    )  # dr + dz
                    pfv.pfcaseth[i] = 0.25e0 * (
                        -drpdz + math.sqrt(drpdz * drpdz + 4.0e0 * areaspf)
                    )

                else:
                    areaspf = 0.0e0  # Resistive coil - no steel needed
                    pfv.pfcaseth[i] = 0.0e0

                # Weight of steel case

                pfv.m_pf_coil_structure[i] = (
                    areaspf
                    * 2.0e0
                    * constants.pi
                    * pfv.r_pf_coil_middle[i]
                    * fwbsv.denstl
                )

                # Mass of heaviest PF coil (tonnes)

                pfv.m_pf_coil_max = max(
                    pfv.m_pf_coil_max,
                    (
                        1.0e-3
                        * (pfv.m_pf_coil_conductor[i] + pfv.m_pf_coil_structure[i])
                    ),
                )
                i = i + 1

        # Find sum of current x turns x radius for all coils for 2015 costs model
        c = 0
        pfv.itr_sum = 0.0e0
        for m in range(pfv.n_pf_coil_groups):
            for _n in range(pfv.n_pf_coils_in_group[m]):
                pfv.itr_sum = pfv.itr_sum + (
                    pfv.r_pf_coil_middle[c]
                    * pfv.n_pf_coil_turns[c]
                    * pfv.c_pf_coil_turn_peak_input[c]
                )
                c = c + 1

        pfv.itr_sum = pfv.itr_sum + (
            (bv.dr_bore + 0.5 * bv.dr_cs)
            * pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1]
            * pfv.c_pf_coil_turn_peak_input[pfv.n_cs_pf_coils - 1]
        )

        # Find Central Solenoid information
        if bv.iohcl != 0:
            self.ohcalc()

        # Summation of weights and current
        pfv.m_pf_coil_conductor_total = 0.0e0
        pfv.m_pf_coil_structure_total = 0.0e0
        pf.ricpf = 0.0e0

        for i in range(pfv.n_cs_pf_coils):
            pfv.m_pf_coil_conductor_total = (
                pfv.m_pf_coil_conductor_total + pfv.m_pf_coil_conductor[i]
            )
            pfv.m_pf_coil_structure_total = (
                pfv.m_pf_coil_structure_total + pfv.m_pf_coil_structure[i]
            )
            pf.ricpf = pf.ricpf + abs(pfv.c_pf_cs_coils_peak_ma[i])

        # Plasma size and shape
        pfv.z_pf_coil_upper[pfv.n_cs_pf_coils] = pv.rminor * pv.kappa
        pfv.z_pf_coil_lower[pfv.n_cs_pf_coils] = -pv.rminor * pv.kappa
        pfv.r_pf_coil_inner[pfv.n_cs_pf_coils] = pv.rmajor - pv.rminor
        pfv.r_pf_coil_outer[pfv.n_cs_pf_coils] = pv.rmajor + pv.rminor
        pfv.n_pf_coil_turns[pfv.n_cs_pf_coils] = 1.0e0

        # Generate coil currents as a function of time using
        # user-provided waveforms etc. (c_pf_coil_turn_peak_input, f_j_cs_start_pulse_end_flat_top, f_j_cs_start_end_flat_top)
        for k in range(6):  # time points
            for i in range(pfv.n_pf_cs_plasma_circuits - 1):
                pfv.c_pf_coil_turn[i, k] = pfv.waves[i, k] * math.copysign(
                    pfv.c_pf_coil_turn_peak_input[i], pfv.c_pf_cs_coils_peak_ma[i]
                )

        # Plasma wave form
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 0] = 0.0e0
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 1] = 0.0e0
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 2] = pv.plasma_current
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 3] = pv.plasma_current
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 4] = pv.plasma_current
        pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, 5] = 0.0e0

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
        rcls,
        zcls,
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
        :param rcls: coords R(i,j), Z(i,j) of coil j in group i (m)
        :type rcls: np.ndarray
        :param zcls: coords R(i,j), Z(i,j) of coil j in group i (m)
        :type zcls: np.ndarray
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
            rcls,
            zcls,
            alfa,
            bfix,
            int(pfv.n_pf_coils_in_group_max),
        )

        # Solve matrix equation
        ccls = self.solv(pfv.n_pf_groups_max, n_pf_coil_groups, nrws, gmat, bvec)

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

            for i in range(pfv.n_pf_coil_groups):
                for ii in range(pfv.n_pf_coil_groups):
                    for ij in range(pfv.n_pf_coils_in_group[ii]):
                        if pf.rcls[ii, ij] <= (  # Outboard TF coil collision
                            pf.rclsnorm - pfv.routr + pfv.r_pf_coil_middle[i]
                        ) and pf.rcls[ii, ij] >= (
                            bv.r_tf_outboard_mid
                            - (0.5 * bv.dr_tf_outboard)
                            - pfv.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1
                        if pf.rcls[ii, ij] <= (  # Inboard TF coil collision
                            bv.dr_bore
                            + bv.dr_cs
                            + bv.dr_cs_precomp
                            + bv.dr_cs_tf_gap
                            + bv.dr_tf_inboard
                            + pfv.r_pf_coil_middle[i]
                        ) and pf.rcls[ii, ij] >= (
                            bv.dr_bore
                            + bv.dr_cs
                            + bv.dr_cs_precomp
                            + bv.dr_cs_tf_gap
                            - pfv.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1
                        if (  # Vertical TF coil collision
                            abs(pf.zcls[ii, ij]) <= bv.hpfu + pfv.r_pf_coil_middle[i]
                            and abs(pf.zcls[ii, ij])
                            >= bv.hpfu
                            - (0.5 * bv.dr_tf_outboard)
                            - pfv.r_pf_coil_middle[i]
                        ):
                            pf_tf_collision += 1

                        if pf_tf_collision >= 1:
                            eh.report_error(277)

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
        hohc = bv.hmax * pfv.f_z_cs_tf_internal

        # Z coordinates of coil edges
        pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1] = hohc
        pfv.z_pf_coil_lower[pfv.n_cs_pf_coils - 1] = -pfv.z_pf_coil_upper[
            pfv.n_cs_pf_coils - 1
        ]

        # (R,Z) coordinates of coil centre
        pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1] = pfv.r_cs_middle
        pfv.z_pf_coil_middle[pfv.n_cs_pf_coils - 1] = 0.0e0

        # Radius of outer edge
        pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1] = pfv.r_cs_middle + 0.5e0 * bv.dr_cs

        # Radius of inner edge
        pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1] = (
            pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1] - bv.dr_cs
        )

        # Total cross-sectional area
        pfv.a_cs_poloidal = 2.0e0 * hohc * bv.dr_cs

        # Maximum current (MA-turns) in central Solenoid, at either BOP or EOF
        if pfv.j_cs_pulse_start > pfv.j_cs_flat_top_end:
            sgn = 1.0e0
            pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1] = (
                sgn * 1.0e-6 * pfv.j_cs_pulse_start * pfv.a_cs_poloidal
            )
        else:
            sgn = -1.0e0
            pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1] = (
                sgn * 1.0e-6 * pfv.j_cs_flat_top_end * pfv.a_cs_poloidal
            )

        # Number of turns
        pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1] = (
            1.0e6
            * abs(pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1])
            / pfv.c_pf_coil_turn_peak_input[pfv.n_cs_pf_coils - 1]
        )

        # Turn vertical cross-sectionnal area
        pfv.a_cs_turn = pfv.a_cs_poloidal / pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1]

        # Depth/width of cs turn conduit
        pfv.dz_cs_turn = (pfv.a_cs_turn / pfv.ld_ratio_cst) ** 0.5

        # length of cs turn conduit
        pfv.dr_cs_turn = pfv.ld_ratio_cst * pfv.dz_cs_turn

        # Radius of turn space = pfv.radius_cs_turn_cable_space
        # Radius of curved outer corrner pfv.r_out_cst = 3mm from literature
        # pfv.ld_ratio_cst = 70 / 22 from literature

        # CS coil turn geometry calculation - stadium shape
        # Literature: https://doi.org/10.1016/j.fusengdes.2017.04.052
        pfv.radius_cs_turn_cable_space = -(
            (pfv.dr_cs_turn - pfv.dz_cs_turn) / constants.pi
        ) + math.sqrt(
            (((pfv.dr_cs_turn - pfv.dz_cs_turn) / constants.pi) ** 2)
            + (
                (
                    (pfv.dr_cs_turn * pfv.dz_cs_turn)
                    - (4 - constants.pi) * (pfv.r_out_cst**2)
                    - (pfv.a_cs_turn * pfv.f_a_cs_steel)
                )
                / constants.pi
            )
        )

        # Thickness of steel conduit in cs turn
        csfv.t_structural_vertical = (
            pfv.dz_cs_turn / 2
        ) - pfv.radius_cs_turn_cable_space
        # In this model the vertical and radial have the same thickness
        csfv.t_structural_radial = csfv.t_structural_vertical
        # add a check for negative conduit thickness
        if csfv.t_structural_radial < 1.0e-3:
            csfv.t_structural_radial = 1.0e-3

        # Non-steel area void fraction for coolant
        pfv.f_a_pf_coil_void[pfv.n_cs_pf_coils - 1] = pfv.f_a_cs_void

        # Peak field at the End-Of-Flattop (EOF)
        # Occurs at inner edge of coil; bmaxoh2 and bzi are of opposite sign at EOF

        # Peak field due to central Solenoid itself
        bmaxoh2 = self.bfmax(
            pfv.j_cs_flat_top_end,
            pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1],
            pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1],
            hohc,
        )

        # Peak field due to other PF coils plus plasma
        timepoint = 5
        bri, bro, bzi, bzo = self.peakb(pfv.n_cs_pf_coils, 99, timepoint)

        pfv.b_cs_peak_flat_top_end = abs(bzi - bmaxoh2)

        # Peak field on outboard side of central Solenoid
        # (self-field is assumed to be zero - long solenoid approximation)
        bohco = abs(bzo)

        # Peak field at the Beginning-Of-Pulse (BOP)
        # Occurs at inner edge of coil; b_cs_peak_pulse_start and bzi are of same sign at BOP
        pfv.b_cs_peak_pulse_start = self.bfmax(
            pfv.j_cs_pulse_start,
            pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1],
            pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1],
            hohc,
        )
        timepoint = 2
        bri, bro, bzi, bzo = self.peakb(pfv.n_cs_pf_coils, 99, timepoint)

        pfv.b_cs_peak_pulse_start = abs(pfv.b_cs_peak_pulse_start + bzi)

        # Maximum field values
        pfv.b_pf_coil_peak[pfv.n_cs_pf_coils - 1] = max(
            pfv.b_cs_peak_flat_top_end, abs(pfv.b_cs_peak_pulse_start)
        )
        pf.bpf2[pfv.n_cs_pf_coils - 1] = max(bohco, abs(bzo))

        # Stress ==> cross-sectional area of supporting steel to use
        if pfv.i_pf_conductor == 0:
            # Superconducting coil

            # New calculation from M. N. Wilson for hoop stress
            pf.sig_hoop = self.hoop_stress(pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1])

            # New calculation from Y. Iwasa for axial stress
            pf.sig_axial, pf.axial_force = self.axial_stress()

            # Allowable (hoop) stress (Pa) alstroh
            # Now a user input
            # alstroh = min( (2.0e0*csytf/3.0e0), (0.5e0*csutf) )

            # Calculation of CS fatigue
            # this is only valid for pulsed reactor design
            if pv.f_c_plasma_inductive > 0.0e-4:
                csfv.n_cycle, csfv.t_crack_radial = self.cs_fatigue.ncycle(
                    pf.sig_hoop,
                    csfv.residual_sig_hoop,
                    csfv.t_crack_vertical,
                    csfv.t_structural_vertical,
                    csfv.t_structural_radial,
                )

            # Now steel area fraction is iteration variable and constraint
            # equation is used for Central Solenoid stress

            # Area of steel in Central Solenoid
            areaspf = pfv.f_a_cs_steel * pfv.a_cs_poloidal

            if pfv.i_cs_stress == 1:
                pfv.s_shear_cs_peak = max(
                    abs(pf.sig_hoop - pf.sig_axial),
                    abs(pf.sig_axial - 0.0e0),
                    abs(0.0e0 - pf.sig_hoop),
                )
            else:
                pfv.s_shear_cs_peak = max(
                    abs(pf.sig_hoop - 0.0e0),
                    abs(0.0e0 - 0.0e0),
                    abs(0.0e0 - pf.sig_hoop),
                )

            # Thickness of hypothetical steel cylinders assumed to encase the CS along
            # its inside and outside edges; in reality, the steel is distributed
            # throughout the conductor
            pfv.pfcaseth[pfv.n_cs_pf_coils - 1] = 0.25e0 * areaspf / hohc

        else:
            areaspf = 0.0e0  # Resistive Central Solenoid - no steel needed
            pfv.pfcaseth[pfv.n_cs_pf_coils - 1] = 0.0e0

        # Weight of steel
        pfv.m_pf_coil_structure[pfv.n_cs_pf_coils - 1] = (
            areaspf
            * 2.0e0
            * constants.pi
            * pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1]
            * fwbsv.denstl
        )

        # Non-steel cross-sectional area
        pfv.awpoh = pfv.a_cs_poloidal - areaspf

        # Issue #97. Fudge to ensure awpoh is positive; result is continuous, smooth and
        # monotonically decreases

        da = 0.0001e0  # 1 cm^2
        if pfv.awpoh < da:
            pfv.awpoh = da * da / (2.0e0 * da - pfv.awpoh)

        # Weight of conductor in central Solenoid
        if pfv.i_pf_conductor == 0:
            pfv.m_pf_coil_conductor[pfv.n_cs_pf_coils - 1] = (
                pfv.awpoh
                * (1.0e0 - pfv.f_a_cs_void)
                * 2.0e0
                * constants.pi
                * pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1]
                * tfv.dcond[pfv.i_cs_superconductor - 1]
            )
        else:
            pfv.m_pf_coil_conductor[pfv.n_cs_pf_coils - 1] = (
                pfv.awpoh
                * (1.0e0 - pfv.f_a_cs_void)
                * 2.0e0
                * constants.pi
                * pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1]
                * constants.dcopper
            )

        if pfv.i_pf_conductor == 0:
            # Allowable coil overall current density at EOF
            # (superconducting coils only)

            (
                jcritwp,
                pfv.jcableoh_eof,
                pfv.j_cs_conductor_critical_flat_top_end,
                tmarg1,
            ) = self.superconpf(
                pfv.b_cs_peak_flat_top_end,
                pfv.f_a_cs_void,
                pfv.fcuohsu,
                (abs(pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1]) / pfv.awpoh)
                * 1.0e6,
                pfv.i_cs_superconductor,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )
            # Strand critical current calculation for costing in $/kAm
            # = superconducting filaments jc * (1 - strand copper fraction)
            if pfv.i_cs_superconductor.item() in {2, 6, 8}:
                pfv.j_crit_str_cs = pfv.j_cs_conductor_critical_flat_top_end
            else:
                pfv.j_crit_str_cs = pfv.j_cs_conductor_critical_flat_top_end * (
                    1 - pfv.fcuohsu
                )

            pfv.j_cs_critical_flat_top_end = jcritwp * pfv.awpoh / pfv.a_cs_poloidal

            # Allowable coil overall current density at BOP

            (
                jcritwp,
                pfv.jcableoh_bop,
                pfv.j_cs_conductor_critical_pulse_start,
                tmarg2,
            ) = self.superconpf(
                pfv.b_cs_peak_pulse_start,
                pfv.f_a_cs_void,
                pfv.fcuohsu,
                (abs(pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1]) / pfv.awpoh)
                * 1.0e6,
                pfv.i_cs_superconductor,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )

            pfv.j_pf_wp_critical[pfv.n_cs_pf_coils - 1] = (
                jcritwp * pfv.awpoh / pfv.a_cs_poloidal
            )
            pfv.j_cs_critical_pulse_start = pfv.j_pf_wp_critical[pfv.n_cs_pf_coils - 1]

            pfv.temp_cs_margin = min(tmarg1, tmarg2)

        else:
            # Resistive power losses (non-superconducting coil)

            pfv.p_cs_resistive_flat_top = (
                2.0e0
                * constants.pi
                * pfv.r_cs_middle
                * pfv.rho_pf_coil
                / (pfv.a_cs_poloidal * (1.0e0 - pfv.f_a_cs_void))
                * (1.0e6 * pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1]) ** 2
            )
            pfv.p_pf_coil_resistive_total_flat_top = (
                pfv.p_pf_coil_resistive_total_flat_top + pfv.p_cs_resistive_flat_top
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
        if bv.iohcl != 0 and i == pfv.n_cs_pf_coils:
            # Peak field is to be calculated at the Central Solenoid itself,
            # so exclude its own contribution; its self field is
            # dealt with externally using routine BFMAX
            kk = 0
        else:
            # Check different times for maximum current
            if (
                abs(
                    pfv.c_pf_cs_coil_pulse_start_ma[i - 1]
                    - pfv.c_pf_cs_coils_peak_ma[i - 1]
                )
                < 1.0e-12
            ):
                it = 2
            elif (
                abs(
                    pfv.c_pf_cs_coil_flat_top_ma[i - 1]
                    - pfv.c_pf_cs_coils_peak_ma[i - 1]
                )
                < 1.0e-12
            ):
                it = 4
            elif (
                abs(
                    pfv.c_pf_cs_coil_pulse_end_ma[i - 1]
                    - pfv.c_pf_cs_coils_peak_ma[i - 1]
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
                sgn = 1.0 if pfv.j_cs_pulse_start > pfv.j_cs_flat_top_end else -1.0

                # Current in each filament representing part of the Central Solenoid
                for iohc in range(pf.nfxf):
                    pf.cfxf[iohc] = (
                        pfv.waves[pfv.n_cs_pf_coils - 1, it - 1]
                        * pfv.j_cs_flat_top_end
                        * sgn
                        * bv.dr_cs
                        * pfv.f_z_cs_tf_internal
                        * bv.hmax
                        / pf.nfxf
                        * 2.0e0
                    )

                kk = pf.nfxf

        # Non-Central Solenoid coils' contributions
        jj = 0
        for iii in range(pfv.n_pf_coil_groups):
            for _jjj in range(pfv.n_pf_coils_in_group[iii]):
                jj = jj + 1
                # Radius, z-coordinate and current for each coil
                if iii == ii - 1:
                    # Self field from coil (Lyle's Method)
                    kk = kk + 1

                    dzpf = pfv.z_pf_coil_upper[jj - 1] - pfv.z_pf_coil_lower[jj - 1]
                    pf.rfxf[kk - 1] = pfv.r_pf_coil_middle[jj - 1]
                    pf.zfxf[kk - 1] = pfv.z_pf_coil_middle[jj - 1] + dzpf * 0.125e0
                    pf.cfxf[kk - 1] = (
                        pfv.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfv.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.r_pf_coil_middle[jj - 1]
                    pf.zfxf[kk - 1] = pfv.z_pf_coil_middle[jj - 1] + dzpf * 0.375e0
                    pf.cfxf[kk - 1] = (
                        pfv.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfv.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.r_pf_coil_middle[jj - 1]
                    pf.zfxf[kk - 1] = pfv.z_pf_coil_middle[jj - 1] - dzpf * 0.125e0
                    pf.cfxf[kk - 1] = (
                        pfv.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfv.waves[jj - 1, it - 1]
                        * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.r_pf_coil_middle[jj - 1]
                    pf.zfxf[kk - 1] = pfv.z_pf_coil_middle[jj - 1] - dzpf * 0.375e0
                    pf.cfxf[kk - 1] = (
                        pfv.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfv.waves[jj - 1, it - 1]
                        * 0.25e6
                    )

                else:
                    # Field from different coil
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.r_pf_coil_middle[jj - 1]
                    pf.zfxf[kk - 1] = pfv.z_pf_coil_middle[jj - 1]
                    pf.cfxf[kk - 1] = (
                        pfv.c_pf_cs_coils_peak_ma[jj - 1]
                        * pfv.waves[jj - 1, it - 1]
                        * 1.0e6
                    )

        # Plasma contribution
        if it > 2:
            kk = kk + 1
            pf.rfxf[kk - 1] = pv.rmajor
            pf.zfxf[kk - 1] = 0.0e0
            pf.cfxf[kk - 1] = pv.plasma_current

        # Calculate the field at the inner and outer edges
        # of the coil of interest
        pf.xind[:kk], bri, bzi, psi = bfield(
            pf.rfxf[:kk],
            pf.zfxf[:kk],
            pf.cfxf[:kk],
            pfv.r_pf_coil_inner[i - 1],
            pfv.z_pf_coil_middle[i - 1],
        )
        pf.xind[:kk], bro, bzo, psi = bfield(
            pf.rfxf[:kk],
            pf.zfxf[:kk],
            pf.cfxf[:kk],
            pfv.r_pf_coil_outer[i - 1],
            pfv.z_pf_coil_middle[i - 1],
        )

        # b_pf_coil_peak and bpf2 for the Central Solenoid are calculated in OHCALC
        if (bv.iohcl != 0) and (i == pfv.n_cs_pf_coils):
            return bri, bro, bzi, bzo

        bpfin = math.sqrt(bri**2 + bzi**2)
        bpfout = math.sqrt(bro**2 + bzo**2)
        for n in range(pfv.n_pf_coils_in_group[ii - 1]):
            pfv.b_pf_coil_peak[i - 1 + n] = bpfin
            pf.bpf2[i - 1 + n] = bpfout

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
            * constants.rmu0
            * h
            * math.log(
                (alpha + math.sqrt(alpha**2 + beta**2))
                / (1.0 + math.sqrt(1.0 + beta**2))
            )
        )

        if beta > 3.0:
            b1 = constants.rmu0 * rj * (b - a)
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
            pf.nef = pfv.n_pf_cs_plasma_circuits - 1
        else:
            pf.nef = pfv.n_pf_cs_plasma_circuits - 2

        pfv.vs_pf_coils_total_ramp = 0.0e0

        for i in range(pf.nef):
            pf.vsdum[i, 0] = (
                pfv.ind_pf_cs_plasma_mutual[pfv.n_pf_cs_plasma_circuits - 1, i]
                * pfv.c_pf_coil_turn[i, 1]
            )
            pf.vsdum[i, 1] = (
                pfv.ind_pf_cs_plasma_mutual[pfv.n_pf_cs_plasma_circuits - 1, i]
                * pfv.c_pf_coil_turn[i, 2]
            )
            pfv.vs_pf_coils_total_ramp = pfv.vs_pf_coils_total_ramp + (
                pf.vsdum[i, 1] - pf.vsdum[i, 0]
            )

        # Central Solenoid startup volt-seconds
        if bv.iohcl != 0:
            pf.vsdum[pfv.n_cs_pf_coils - 1, 0] = (
                pfv.ind_pf_cs_plasma_mutual[
                    pfv.n_pf_cs_plasma_circuits - 1, pfv.n_pf_cs_plasma_circuits - 2
                ]
                * pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 2, 1]
            )
            pf.vsdum[pfv.n_cs_pf_coils - 1, 1] = (
                pfv.ind_pf_cs_plasma_mutual[
                    pfv.n_pf_cs_plasma_circuits - 1, pfv.n_pf_cs_plasma_circuits - 2
                ]
                * pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 2, 2]
            )
            pfv.vs_cs_ramp = (
                pf.vsdum[pfv.n_cs_pf_coils - 1, 1] - pf.vsdum[pfv.n_cs_pf_coils - 1, 0]
            )

        # Total available volt-seconds for start-up
        pfv.vs_cs_pf_total_ramp = pfv.vs_cs_ramp + pfv.vs_pf_coils_total_ramp

        # Burn volt-seconds
        if bv.iohcl != 0:
            pf.vsdum[pfv.n_cs_pf_coils - 1, 2] = (
                pfv.ind_pf_cs_plasma_mutual[
                    pfv.n_pf_cs_plasma_circuits - 1, pfv.n_pf_cs_plasma_circuits - 2
                ]
                * pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 2, 4]
            )
            pfv.vs_cs_burn = (
                pf.vsdum[pfv.n_cs_pf_coils - 1, 2] - pf.vsdum[pfv.n_cs_pf_coils - 1, 1]
            )

        # PF volt-seconds during burn
        pfv.vs_pf_coils_total_burn = 0.0e0
        for i in range(pf.nef):
            pf.vsdum[i, 2] = (
                pfv.ind_pf_cs_plasma_mutual[pfv.n_pf_cs_plasma_circuits - 1, i]
                * pfv.c_pf_coil_turn[i, 4]
            )
            pfv.vs_pf_coils_total_burn = pfv.vs_pf_coils_total_burn + (
                pf.vsdum[i, 2] - pf.vsdum[i, 1]
            )

        pfv.vs_cs_pf_total_burn = pfv.vs_cs_burn + pfv.vs_pf_coils_total_burn

        pfv.vs_cs_pf_total_pulse = pfv.vs_cs_pf_total_ramp + pfv.vs_cs_pf_total_burn
        pfv.vs_pf_coils_total_pulse = (
            pfv.vs_pf_coils_total_ramp + pfv.vs_pf_coils_total_burn
        )
        pfv.vs_cs_total_pulse = pfv.vs_cs_burn + pfv.vs_cs_ramp

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
        a = pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1]

        # Outer radius of central Solenoid [m]
        b = pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1]

        # alpha
        alpha = b / a

        # epsilon
        epsilon = r / a

        # Field at inner radius of coil [T]
        b_a = pfv.b_cs_peak_pulse_start

        # Field at outer radius of coil [T]
        # Assume to be 0 for now
        b_b = 0.0e0

        # current density [A/m^2]
        j = pfv.j_cs_pulse_start

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

        return s_hoop_nom / pfv.f_a_cs_steel

    def axial_stress(self):
        """Calculation of axial stress of central solenoid.

        author: J Morris, CCFE, Culham Science Centre
        This routine calculates the axial stress of the central solenoid
        from "Case studies in superconducting magnets", Y. Iwasa, Springer

        :return: unsmeared axial stress [MPa], axial force [N]
        :rtype: tuple[float, float]
        """
        b = pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1]

        # Half height of central Solenoid [m]
        hl = pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1]

        # Central Solenoid current [A]
        ni = pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1] * 1.0e6

        # kb term for elliptical integrals
        # kb2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + hl**2))
        kb2 = (4.0e0 * b**2) / (4.0e0 * b**2 + hl**2)

        # k2b term for elliptical integrals
        # k2b2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + 4.0e0*hl**2))
        k2b2 = (4.0e0 * b**2) / (4.0e0 * b**2 + 4.0e0 * hl**2)

        # term 1
        axial_term_1 = -(constants.rmu0 / 2.0e0) * (ni / (2.0e0 * hl)) ** 2

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
        area_ax = constants.pi * (
            pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1] ** 2
            - pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1] ** 2
        )

        # calculate unsmeared axial stress [MPa]
        s_axial = axial_force / (pfv.f_a_cs_steel * 0.5 * area_ax)

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
        rc = np.zeros(pfv.ngc2 + nohmax)
        zc = np.zeros(pfv.ngc2 + nohmax)
        xc = np.zeros(pfv.ngc2 + nohmax)
        cc = np.zeros(pfv.ngc2 + nohmax)
        xcin = np.zeros(pfv.ngc2 + nohmax)
        xcout = np.zeros(pfv.ngc2 + nohmax)
        rplasma = np.zeros(nplas)
        zplasma = np.zeros(nplas)

        pfv.ind_pf_cs_plasma_mutual[:, :] = 0.0

        # Break Central Solenoid into noh segments
        #
        # Choose noh so that the radial thickness of the coil is not thinner
        # than each segment is tall, i.e. the segments are pancake-like,
        # for the benefit of the mutual inductance calculations later

        noh = int(
            math.ceil(
                2.0e0
                * pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1]
                / (
                    pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1]
                    - pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1]
                )
            )
        )

        if noh > nohmax:
            eh.idiags[0] = noh
            eh.idiags[1] = nohmax
            eh.fdiags[0] = bv.dr_cs
            eh.report_error(73)

        noh = min(noh, nohmax)

        # TODO In FNSF case, noh = -7! noh should always be positive. Fortran
        # array allocation with -ve bound previously coerced to 0
        if noh < 0:
            noh = 0

        roh = np.zeros(noh)
        zoh = np.zeros(noh)

        if bv.iohcl != 0:
            roh[:] = pfv.r_cs_middle

            delzoh = (
                2.0e0 * pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1] / noh
            )  # z_pf_coil_upper(n_cs_pf_coils) is the half-height of the coil
            for i in range(noh):
                zoh[i] = pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1] - delzoh * (
                    0.5e0 + i
                )

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
                # eh.fdiags[0] = bv.dr_cs
                # eh.fdiags[1] = delzoh
                # eh.report_error(74)
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

            pfv.ind_pf_cs_plasma_mutual[
                pfv.n_pf_cs_plasma_circuits - 1, pfv.n_cs_pf_coils - 1
            ] = xohpl / (nplas * noh) * pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1]
            pfv.ind_pf_cs_plasma_mutual[
                pfv.n_cs_pf_coils - 1, pfv.n_pf_cs_plasma_circuits - 1
            ] = pfv.ind_pf_cs_plasma_mutual[
                pfv.n_pf_cs_plasma_circuits - 1, pfv.n_cs_pf_coils - 1
            ]

        # Plasma self inductance
        pfv.ind_pf_cs_plasma_mutual[
            pfv.n_pf_cs_plasma_circuits - 1, pfv.n_pf_cs_plasma_circuits - 1
        ] = pv.ind_plasma

        # PF coil / plasma mutual inductances
        ncoils = 0

        for i in range(pfv.n_pf_coil_groups):
            xpfpl = 0.0
            ncoils = ncoils + pfv.n_pf_coils_in_group[i]
            rp = pfv.r_pf_coil_middle[ncoils - 1]
            zp = pfv.z_pf_coil_middle[ncoils - 1]
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
            for ii in range(nplas):
                xpfpl = xpfpl + xc[ii]

            for j in range(pfv.n_pf_coils_in_group[i]):
                ncoilj = ncoils + 1 - (j + 1)
                pfv.ind_pf_cs_plasma_mutual[
                    ncoilj - 1, pfv.n_pf_cs_plasma_circuits - 1
                ] = xpfpl / nplas * pfv.n_pf_coil_turns[ncoilj - 1]
                pfv.ind_pf_cs_plasma_mutual[
                    pfv.n_pf_cs_plasma_circuits - 1, ncoilj - 1
                ] = pfv.ind_pf_cs_plasma_mutual[
                    ncoilj - 1, pfv.n_pf_cs_plasma_circuits - 1
                ]

        if bv.iohcl != 0:
            # Central Solenoid self inductance
            a = pfv.r_cs_middle  # mean radius of coil
            b = 2.0e0 * pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1]  # length of coil
            c = (
                pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1]
                - pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1]
            )  # radial winding thickness
            pfv.ind_pf_cs_plasma_mutual[
                pfv.n_cs_pf_coils - 1, pfv.n_cs_pf_coils - 1
            ] = self.selfinductance(a, b, c, pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1])

            # Central Solenoid / PF coil mutual inductances
            for i in range(noh):
                rc[i] = roh[i]
                zc[i] = zoh[i]

            ncoils = 0
            for i in range(pfv.n_pf_coil_groups):
                xohpf = 0.0
                ncoils = ncoils + pfv.n_pf_coils_in_group[i]
                rp = pfv.r_pf_coil_middle[ncoils - 1]
                zp = pfv.z_pf_coil_middle[ncoils - 1]
                xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
                for ii in range(noh):
                    xohpf = xohpf + xc[ii]

                for j in range(pfv.n_pf_coils_in_group[i]):
                    ncoilj = ncoils + 1 - (j + 1)
                    pfv.ind_pf_cs_plasma_mutual[ncoilj - 1, pfv.n_cs_pf_coils - 1] = (
                        xohpf
                        * pfv.n_pf_coil_turns[ncoilj - 1]
                        * pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1]
                        / noh
                    )
                    pfv.ind_pf_cs_plasma_mutual[pfv.n_cs_pf_coils - 1, ncoilj - 1] = (
                        pfv.ind_pf_cs_plasma_mutual[ncoilj - 1, pfv.n_cs_pf_coils - 1]
                    )

        # PF coil - PF coil inductances
        if bv.iohcl == 0:
            pf.nef = pfv.n_cs_pf_coils
        else:
            pf.nef = pfv.n_cs_pf_coils - 1

        for i in range(pf.nef):
            for j in range(pf.nef - 1):
                jj = j + 1 + 1 if j >= i else j + 1

                zc[j] = pfv.z_pf_coil_middle[jj - 1]
                rc[j] = pfv.r_pf_coil_middle[jj - 1]

            rp = pfv.r_pf_coil_middle[i]
            zp = pfv.z_pf_coil_middle[i]
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
            for k in range(pf.nef):
                if k < i:
                    pfv.ind_pf_cs_plasma_mutual[i, k] = (
                        xc[k] * pfv.n_pf_coil_turns[k] * pfv.n_pf_coil_turns[i]
                    )
                elif k == i:
                    rl = abs(
                        pfv.z_pf_coil_upper[k] - pfv.z_pf_coil_lower[k]
                    ) / math.sqrt(constants.pi)
                    pfv.ind_pf_cs_plasma_mutual[k, k] = (
                        constants.rmu0
                        * pfv.n_pf_coil_turns[k] ** 2
                        * pfv.r_pf_coil_middle[k]
                        * (math.log(8.0e0 * pfv.r_pf_coil_middle[k] / rl) - 1.75e0)
                    )
                else:
                    pfv.ind_pf_cs_plasma_mutual[i, k] = (
                        xc[k - 1] * pfv.n_pf_coil_turns[k] * pfv.n_pf_coil_turns[i]
                    )

        # Output section
        if not output:
            return

        op.oheadr(self.outfile, "PF Coil Inductances")
        op.ocmmnt(self.outfile, "Inductance matrix [H]:")
        op.oblnkl(self.outfile)

        with np.printoptions(precision=1):
            for ig in range(pf.nef):
                op.write(
                    self.outfile,
                    f"{ig}\t{pfv.ind_pf_cs_plasma_mutual[: pfv.n_pf_cs_plasma_circuits, ig]}",
                )

            if bv.iohcl != 0:
                op.write(
                    self.outfile,
                    f"CS\t\t\t{pfv.ind_pf_cs_plasma_mutual[: pfv.n_pf_cs_plasma_circuits, pfv.n_pf_cs_plasma_circuits - 2]}",
                )

            op.write(
                self.outfile,
                f"Plasma\t{pfv.ind_pf_cs_plasma_mutual[: pfv.n_pf_cs_plasma_circuits, pfv.n_pf_cs_plasma_circuits - 1]}",
            )

    def output_cs_structure(self):
        op.oheadr(self.outfile, "Central Solenoid Structure")

        op.ocmmnt(self.outfile, "CS turn structure")
        op.oblnkl(self.outfile)

        op.ovarre(
            self.outfile,
            "Poloidal area of a CS turn [m^2]",
            "(a_cs_turn)",
            pfv.a_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Radial width a CS turn [m^2]",
            "(dz_cs_turn)",
            pfv.dz_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Length of a CS turn [m]",
            "(dr_cs_turn)",
            pfv.dr_cs_turn,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Length to diameter ratio of a CS turn",
            "(ld_ratio_cst)",
            pfv.ld_ratio_cst,
            "OP ",
        )
        op.ovarre(
            self.outfile,
            "Radius of CS turn cable space [m]",
            "(radius_cs_turn_cable_space)",
            pfv.radius_cs_turn_cable_space,
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
            pfv.r_out_cst,
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
            if pfv.i_pf_conductor == 0:
                op.ocmmnt(self.outfile, "Superconducting central solenoid")

                op.ovarin(
                    self.outfile,
                    "Central solenoid superconductor material",
                    "(i_cs_superconductor)",
                    pfv.i_cs_superconductor,
                )

                if pfv.i_cs_superconductor == 1:
                    op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
                elif pfv.i_cs_superconductor == 2:
                    op.ocmmnt(
                        self.outfile, "  (Bi-2212 high temperature superconductor)"
                    )
                elif pfv.i_cs_superconductor == 3:
                    op.ocmmnt(self.outfile, "  (NbTi)")
                elif pfv.i_cs_superconductor == 4:
                    op.ocmmnt(
                        self.outfile,
                        "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                    )
                elif pfv.i_cs_superconductor == 5:
                    op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
                elif pfv.i_cs_superconductor == 6:
                    op.ocmmnt(self.outfile, " (REBCO HTS)")
                elif pfv.i_cs_superconductor == 7:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                    )
                elif pfv.i_cs_superconductor == 8:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for REBCO)",
                    )
                elif pfv.i_cs_superconductor == 9:
                    op.ocmmnt(
                        self.outfile,
                        " (Hazelton experimental data + Zhai conceptual model for REBCO)",
                    )

                op.osubhd(self.outfile, "Central Solenoid Current Density Limits :")
                op.ovarre(
                    self.outfile,
                    "Maximum field at Beginning Of Pulse (T)",
                    "(b_cs_peak_pulse_start)",
                    pfv.b_cs_peak_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at BOP (A/m2)",
                    "(j_cs_conductor_critical_pulse_start)",
                    pfv.j_cs_conductor_critical_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical cable current density at BOP (A/m2)",
                    "(jcableoh_bop)",
                    pfv.jcableoh_bop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at BOP (A/m2)",
                    "(j_cs_critical_pulse_start)",
                    pfv.j_cs_critical_pulse_start,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at BOP (A/m2)",
                    "(j_cs_pulse_start)",
                    pfv.j_cs_pulse_start,
                    "OP ",
                )
                op.oblnkl(self.outfile)
                op.ovarre(
                    self.outfile,
                    "Maximum field at End Of Flattop (T)",
                    "(b_cs_peak_flat_top_end)",
                    pfv.b_cs_peak_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at EOF (A/m2)",
                    "(j_cs_conductor_critical_flat_top_end)",
                    pfv.j_cs_conductor_critical_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical cable current density at EOF (A/m2)",
                    "(jcableoh_eof)",
                    pfv.jcableoh_eof,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at EOF (A/m2)",
                    "(j_cs_critical_flat_top_end)",
                    pfv.j_cs_critical_flat_top_end,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at EOF (A/m2)",
                    "(j_cs_flat_top_end)",
                    pfv.j_cs_flat_top_end,
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
                    pfv.a_cs_poloidal,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS conductor+void cross-sectional area (m2)",
                    "(awpoh)",
                    pfv.awpoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "   CS conductor cross-sectional area (m2)",
                    "(awpoh*(1-f_a_cs_void))",
                    pfv.awpoh * (1.0e0 - pfv.f_a_cs_void),
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "   CS void cross-sectional area (m2)",
                    "(awpoh*f_a_cs_void)",
                    pfv.awpoh * pfv.f_a_cs_void,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel cross-sectional area (m2)",
                    "(a_cs_poloidal-awpoh)",
                    pfv.a_cs_poloidal - pfv.awpoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel area fraction",
                    "(f_a_cs_steel)",
                    pfv.f_a_cs_steel,
                )
                if pfv.i_cs_stress == 1:
                    op.ocmmnt(self.outfile, "Hoop + axial stress considered")
                else:
                    op.ocmmnt(self.outfile, "Only hoop stress considered")

                op.ovarin(
                    self.outfile,
                    "Switch for CS stress calculation",
                    "(i_cs_stress)",
                    pfv.i_cs_stress,
                )
                op.ovarre(
                    self.outfile,
                    "Allowable stress in CS steel (Pa)",
                    "(alstroh)",
                    pfv.alstroh,
                )
                op.ovarre(
                    self.outfile,
                    "Hoop stress in CS steel (Pa)",
                    "(sig_hoop)",
                    pf.sig_hoop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Axial stress in CS steel (Pa)",
                    "(sig_axial)",
                    pf.sig_axial,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Maximum shear stress in CS steel for the Tresca criterion (Pa)",
                    "(s_shear_cs_peak)",
                    pfv.s_shear_cs_peak,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Axial force in CS (N)",
                    "(axial_force)",
                    pf.axial_force,
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
                    pfv.fcuohsu,
                )
                # If REBCO material is used, print copperaoh_m2
                if (
                    pfv.i_cs_superconductor == 6
                    or pfv.i_cs_superconductor == 8
                    or pfv.i_cs_superconductor == 9
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
                    pfv.f_a_cs_void,
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
                    "(temp_cs_margin)",
                    pfv.temp_cs_margin,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Minimum permitted temperature margin (K)",
                    "(tmargmin_cs)",
                    tfv.tmargmin_cs,
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
                        pfv.a_cs_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn length (m)",
                        "(dr_cs_turn)",
                        pfv.dr_cs_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn internal cable space radius (m)",
                        "(radius_cs_turn_cable_space)",
                        pfv.radius_cs_turn_cable_space,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn width (m)",
                        "(dz_cs_turn)",
                        pfv.dz_cs_turn,
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
                    abs(pfv.j_cs_flat_top_end)
                    > 0.99e0 * abs(numerics.boundu[37] * pfv.j_cs_critical_flat_top_end)
                ) or (
                    abs(pfv.j_cs_pulse_start)
                    > 0.99e0 * abs(numerics.boundu[38] * pfv.j_cs_critical_pulse_start)
                ):
                    pf.cslimit = True
                if pfv.temp_cs_margin < 1.01e0 * tfv.tmargmin_cs:
                    pf.cslimit = True
                if not pf.cslimit:
                    eh.report_error(135)

                # Check whether CS coil currents are feasible from engineering POV
                if ctv.fjohc > 0.7:
                    eh.report_error(286)
                if ctv.fjohc0 > 0.7:
                    eh.report_error(287)

                # REBCO fractures in strains above ~+/- 0.7%
                if (
                    pfv.i_cs_superconductor == 6
                    or pfv.i_cs_superconductor == 8
                    or pfv.i_cs_superconductor == 9
                ) and abs(tfv.str_cs_con_res) > 0.7e-2:
                    eh.report_error(262)

                if (
                    pfv.i_pf_superconductor == 6
                    or pfv.i_pf_superconductor == 8
                    or pfv.i_pf_superconductor == 9
                ) and abs(tfv.str_pf_con_res) > 0.7e-2:
                    eh.report_error(263)

            else:
                op.ocmmnt(self.outfile, "Resistive central solenoid")

        if pfv.i_pf_conductor == 0:
            op.oblnkl(self.outfile)
            op.ocmmnt(self.outfile, "Superconducting PF coils")

            op.ovarin(
                self.outfile,
                "PF coil superconductor material",
                "(i_pf_superconductor)",
                pfv.i_pf_superconductor,
            )

            if pfv.i_pf_superconductor == 1:
                op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
            elif pfv.i_pf_superconductor == 2:
                op.ocmmnt(self.outfile, "  (Bi-2212 high temperature superconductor)")
            elif pfv.i_pf_superconductor == 3:
                op.ocmmnt(self.outfile, "  (NbTi)")
            elif pfv.i_pf_superconductor == 4:
                op.ocmmnt(
                    self.outfile,
                    "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                )
            elif pfv.i_pf_superconductor == 5:
                op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
            elif pfv.i_pf_superconductor == 6:
                op.ocmmnt(
                    self.outfile,
                    " (REBCO 2nd generation HTS superconductor in CrCo strand)",
                )
            elif pfv.i_pf_superconductor == 7:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                )
            elif pfv.i_pf_superconductor == 8:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for REBCO)",
                )
            elif pfv.i_pf_superconductor == 9:
                op.ocmmnt(
                    self.outfile,
                    " (Hazelton experimental data + Zhai conceptual model for REBCO)",
                )

            op.ovarre(
                self.outfile,
                "Copper fraction in conductor",
                "(fcupfsu)",
                pfv.fcupfsu,
            )

            op.osubhd(self.outfile, "PF Coil Case Stress :")
            op.ovarre(
                self.outfile,
                "Maximum permissible tensile stress (MPa)",
                "(sigpfcalw)",
                pfv.sigpfcalw,
            )
            op.ovarre(
                self.outfile,
                "JxB hoop force fraction supported by case",
                "(sigpfcf)",
                pfv.sigpfcf,
            )

        else:
            op.oblnkl(self.outfile)
            op.ocmmnt(self.outfile, "Resistive PF coils")

            op.osubhd(self.outfile, "Resistive Power :")
            op.ovarre(
                self.outfile,
                "PF coil resistive power (W)",
                "(p_pf_coil_resistive_total_flat_top)",
                pfv.p_pf_coil_resistive_total_flat_top,
                "OP ",
            )
            if bv.iohcl != 0:
                op.ovarre(
                    self.outfile,
                    "Central solenoid resistive power (W)",
                    "(p_cs_resistive_flat_top)",
                    pfv.p_cs_resistive_flat_top,
                    "OP ",
                )

        # pf.nef is the number of coils excluding the Central Solenoid
        pf.nef = pfv.n_cs_pf_coils
        if bv.iohcl != 0:
            pf.nef = pf.nef - 1

        op.osubhd(self.outfile, "Geometry of PF coils, central solenoid and plasma:")
        op.write(
            self.outfile,
            "coil\t\t\tR(m)\t\tZ(m)\t\tdR(m)\t\tdZ(m)\t\tturns",
        )
        op.oblnkl(self.outfile)

        # PF coils
        for k in range(pf.nef):
            op.write(
                self.outfile,
                f"PF {k}\t\t\t{pfv.r_pf_coil_middle[k]:.2e}\t{pfv.z_pf_coil_middle[k]:.2e}\t{pfv.r_pf_coil_outer[k] - pfv.r_pf_coil_inner[k]:.2e}\t{abs(pfv.z_pf_coil_upper[k] - pfv.z_pf_coil_lower[k]):.2e}\t{pfv.n_pf_coil_turns[k]:.2e}",
            )

        for k in range(pf.nef):
            op.ovarre(
                self.mfile,
                f"PF coil {k} radius (m)",
                f"(r_pf_coil_middle[{k}])",
                pfv.r_pf_coil_middle[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical position (m)",
                f"(z_pf_coil_middle[{k}])",
                pfv.z_pf_coil_middle[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} radial thickness (m)",
                f"(pfdr({k}))",
                pfv.r_pf_coil_outer[k] - pfv.r_pf_coil_inner[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical thickness (m)",
                f"(pfdz({k}))",
                pfv.z_pf_coil_upper[k] - pfv.z_pf_coil_lower[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} turns",
                f"(n_pf_coil_turns[{k}])",
                pfv.n_pf_coil_turns[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} current (MA)",
                f"(c_pf_cs_coils_peak_ma[{k}])",
                pfv.c_pf_cs_coils_peak_ma[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} field (T)",
                f"(b_pf_coil_peak[{k}])",
                pfv.b_pf_coil_peak[k],
            )
        self.tf_pf_collision_detector()

        # Central Solenoid, if present
        if bv.iohcl != 0:
            op.write(
                self.outfile,
                f"CS\t\t\t\t{pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.z_pf_coil_middle[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1] - pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1]:.2e}\t{abs(pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1] - pfv.z_pf_coil_lower[pfv.n_cs_pf_coils - 1]):.2e}\t{pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.pfcaseth[pfv.n_cs_pf_coils - 1]:.2e}",
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radius (m)",
                "(r_pf_coil_middle[n_cs_pf_coils-1])",
                pfv.r_pf_coil_middle[pfv.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical position (m)",
                "(z_pf_coil_middle[n_cs_pf_coils-1])",
                pfv.z_pf_coil_middle[pfv.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radial thickness (m)",
                "(ohdr)",
                (
                    pfv.r_pf_coil_outer[pfv.n_cs_pf_coils - 1]
                    - pfv.r_pf_coil_inner[pfv.n_cs_pf_coils - 1]
                ),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical thickness (m)",
                "(ohdz)",
                (
                    pfv.z_pf_coil_upper[pfv.n_cs_pf_coils - 1]
                    - pfv.z_pf_coil_lower[pfv.n_cs_pf_coils - 1]
                ),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid turns",
                "(n_pf_coil_turns[n_cs_pf_coils-1])",
                pfv.n_pf_coil_turns[pfv.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid current (MA)",
                "(c_pf_cs_coils_peak_ma[n_cs_pf_coils-1])",
                pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid field (T)",
                "(b_pf_coil_peak[n_cs_pf_coils-1])",
                pfv.b_pf_coil_peak[pfv.n_cs_pf_coils - 1],
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
        for k in range(pf.nef):
            if pfv.i_pf_conductor == 0:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfv.c_pf_cs_coils_peak_ma[k]:.2e}\t{pfv.j_pf_wp_critical[k]:.2e}\t{pfv.j_pf_coil_wp_peak[k]:.2e}\t{pfv.j_pf_coil_wp_peak[k] / pfv.j_pf_wp_critical[k]:.2e}\t{pfv.m_pf_coil_conductor[k]:.2e}\t{pfv.m_pf_coil_structure[k]:.2e}\t{pfv.b_pf_coil_peak[k]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfv.c_pf_cs_coils_peak_ma[k]:.2e}\t-1.0e0\t{pfv.j_pf_coil_wp_peak[k]:.2e}\t1.0e0\t{pfv.m_pf_coil_conductor[k]:.2e}\t{pfv.m_pf_coil_structure[k]:.2e}\t{pfv.b_pf_coil_peak[k]:.2e}\t",
                )

        # Central Solenoid, if present
        if bv.iohcl != 0:
            if pfv.i_pf_conductor == 0:
                # Issue #328
                op.write(
                    self.outfile,
                    f"CS\t\t{pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.j_pf_wp_critical[pfv.n_cs_pf_coils - 1]:.2e}\t{max(abs(pfv.j_cs_pulse_start), abs(pfv.j_cs_flat_top_end)):.2e}\t{max(abs(pfv.j_cs_pulse_start), abs(pfv.j_cs_flat_top_end)) / pfv.j_pf_wp_critical[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.m_pf_coil_conductor[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.m_pf_coil_structure[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.b_pf_coil_peak[pfv.n_cs_pf_coils - 1]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"CS\t\t{pfv.c_pf_cs_coils_peak_ma[pfv.n_cs_pf_coils - 1]:.2e}\t-1.0e0\t{max(abs(pfv.j_cs_pulse_start)):.2e}\t{abs(pfv.j_cs_flat_top_end):.2e}\t1.0e0\t{pfv.m_pf_coil_conductor[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.m_pf_coil_structure[pfv.n_cs_pf_coils - 1]:.2e}\t{pfv.b_pf_coil_peak[pfv.n_cs_pf_coils - 1]:.2e}",
                )

        # Miscellaneous totals
        op.write(
            self.outfile,
            "\t" * 1 + "------" + "\t" * 8 + "---------" + "\t" + "---------",
        )

        op.write(
            self.outfile,
            "\t" * 1
            + f"{pf.ricpf:.2e}"
            + "\t" * 7
            + f"{pfv.m_pf_coil_conductor_total:.2e}\t{pfv.m_pf_coil_structure_total:.2e}",
        )

        op.osubhd(self.outfile, "PF coil current scaling information :")
        op.ovarre(
            self.outfile, "Sum of squares of residuals ", "(ssq0)", pf.ssq0, "OP "
        )
        op.ovarre(self.outfile, "Smoothing parameter ", "(alfapf)", pfv.alfapf)

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
            f"PF coils:\t\t{pfv.vs_pf_coils_total_ramp:.2f}\t\t\t\t{pfv.vs_pf_coils_total_burn:.2f}\t\t\t{pfv.vs_pf_coils_total_pulse:.2f}",
        )
        op.write(
            self.outfile,
            f"CS coil:\t\t{pfv.vs_cs_ramp:.2f}\t\t\t\t{pfv.vs_cs_burn:.2f}\t\t\t{pfv.vs_cs_total_pulse:.2f}",
        )
        op.write(
            self.outfile, "\t" * 3 + "-" * 7 + "\t" * 4 + "-" * 7 + "\t" * 3 + "-" * 7
        )
        op.write(
            self.outfile,
            f"Total:\t\t\t{pfv.vs_cs_pf_total_ramp:.2f}\t\t\t\t{pfv.vs_cs_pf_total_burn:.2f}\t\t\t{pfv.vs_cs_pf_total_pulse:.2f}",
        )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Total volt-second consumption by coils (Wb)",
            "(vs_cs_pf_total_pulse)",
            f"{pfv.vs_cs_pf_total_pulse:.2}",
            "OP",
        )

        op.osubhd(self.outfile, "Summary of volt-second consumption by circuit (Wb):")

        op.write(self.outfile, "circuit\t\t\tBOP\t\t\tBOF\t\tEOF")
        op.oblnkl(self.outfile)

        for k in range(pf.nef):
            op.write(
                self.outfile,
                f"\t{k}\t\t\t{pf.vsdum[k, 0]:.3f}\t\t\t{pf.vsdum[k, 1]:.3f}\t\t{pf.vsdum[k, 2]:.3f}",
            )

        op.write(
            self.outfile,
            f"\tCS coil\t\t\t{pf.vsdum[pfv.n_cs_pf_coils - 1, 0]:.3f}\t\t\t{pf.vsdum[pfv.n_cs_pf_coils - 1, 1]:.3f}\t\t{pf.vsdum[pfv.n_cs_pf_coils - 1, 2]:.3f}",
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
            label = f2py_compatible_to_string(tv.timelabel[k])
            line += f"\t\t{label}"
        op.write(self.outfile, line)

        op.ocmmnt(self.outfile, "circuit")

        for k in range(pfv.n_pf_cs_plasma_circuits - 1):
            line = f"\t{k}\t\t"
            for jj in range(6):
                line += f"\t{pfv.c_pf_coil_turn[k, jj] * pfv.n_pf_coil_turns[k]:.3e}"
            op.write(self.outfile, line)

        line = "Plasma (A)\t\t"
        for jj in range(6):
            line += f"\t{pfv.c_pf_coil_turn[pfv.n_pf_cs_plasma_circuits - 1, jj]:.3e}"

        op.write(self.outfile, line)

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "This consists of: CS coil field balancing:")
        for k in range(pfv.n_pf_cs_plasma_circuits - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{pfv.c_pf_coil_turn[k, 0] * pfv.n_pf_coil_turns[k]:.3e}\t"
                    f"{pfv.c_pf_coil_turn[k, 1] * pfv.n_pf_coil_turns[k]:.3e}\t"
                    f"{-pfv.c_pf_coil_turn[k, 1] * pfv.n_pf_coil_turns[k] * (pfv.f_j_cs_start_end_flat_top / pfv.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{-pfv.c_pf_coil_turn[k, 1] * pfv.n_pf_coil_turns[k] * (pfv.f_j_cs_start_end_flat_top / pfv.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{-pfv.c_pf_coil_turn[k, 1] * pfv.n_pf_coil_turns[k] * (1.0e0 / pfv.f_j_cs_start_pulse_end_flat_top):.3e}\t"
                    f"{pfv.c_pf_coil_turn[k, 5] * pfv.n_pf_coil_turns[k]:.3e}"
                ),
            )

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "And: equilibrium field:")
        for k in range(pfv.n_pf_cs_plasma_circuits - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{0.0:.3e}\t{0.0:.3e}\t"
                    f"{(pfv.c_pf_coil_turn[k, 2] + pfv.c_pf_coil_turn[k, 1] * pfv.f_j_cs_start_end_flat_top / pfv.f_j_cs_start_pulse_end_flat_top) * pfv.n_pf_coil_turns[k]:.3e}\t"
                    f"{(pfv.c_pf_coil_turn[k, 3] + pfv.c_pf_coil_turn[k, 1] * pfv.f_j_cs_start_end_flat_top / pfv.f_j_cs_start_pulse_end_flat_top) * pfv.n_pf_coil_turns[k]:.3e}\t"
                    f"{(pfv.c_pf_coil_turn[k, 4] + pfv.c_pf_coil_turn[k, 1] * 1.0e0 / pfv.f_j_cs_start_pulse_end_flat_top) * pfv.n_pf_coil_turns[k]:.3e}\t"
                    "0.0e0"
                ),
            )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Pulse / end of flat-top",
            "(f_j_cs_start_pulse_end_flat_top)",
            pfv.f_j_cs_start_pulse_end_flat_top,
        )
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Flat-top / end of flat-top",
            "(f_j_cs_start_end_flat_top)",
            pfv.f_j_cs_start_end_flat_top,
            "OP ",
        )

        op.oshead(self.outfile, "PF Circuit Waveform Data")
        op.ovarin(
            self.outfile,
            "Number of PF circuits including CS and plasma",
            "(n_pf_cs_plasma_circuits)",
            pfv.n_pf_cs_plasma_circuits,
        )
        for k in range(pfv.n_pf_cs_plasma_circuits):
            for jjj in range(6):
                if k == pfv.n_pf_cs_plasma_circuits - 1:
                    circuit_name = f"Plasma Time point {jjj} (A)"
                    circuit_var_name = f"(plasmat{jjj})"
                elif k == pfv.n_pf_cs_plasma_circuits - 2:
                    circuit_name = f"CS Circuit Time point {jjj} (A)"
                    circuit_var_name = f"(cs t{jjj})"
                else:
                    circuit_name = f"PF Circuit {k} Time point {jjj} (A)"
                    circuit_var_name = f"(pfc{k}t{jjj})"

                op.ovarre(
                    self.outfile,
                    circuit_name,
                    circuit_var_name,
                    pfv.c_pf_coil_turn[k, jjj] * pfv.n_pf_coil_turns[k],
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
        nplas = pfv.n_cs_pf_coils + 1
        for it in range(6):
            pfv.waves[nplas - 1, it] = 1.0e0

        for ic in range(pfv.n_cs_pf_coils):
            # Find where the peak current occurs
            # Beginning of pulse, t = t_precharge
            if (
                abs(pfv.c_pf_cs_coil_pulse_start_ma[ic])
                >= abs(pfv.c_pf_cs_coil_pulse_end_ma[ic])
            ) and (
                abs(pfv.c_pf_cs_coil_pulse_start_ma[ic])
                >= abs(pfv.c_pf_cs_coil_flat_top_ma[ic])
            ):
                pfv.c_pf_cs_coils_peak_ma[ic] = pfv.c_pf_cs_coil_pulse_start_ma[ic]

            # Beginning of flat-top, t = t_precharge + t_current_ramp_up
            if (
                abs(pfv.c_pf_cs_coil_flat_top_ma[ic])
                >= abs(pfv.c_pf_cs_coil_pulse_start_ma[ic])
            ) and (
                abs(pfv.c_pf_cs_coil_flat_top_ma[ic])
                >= abs(pfv.c_pf_cs_coil_pulse_end_ma[ic])
            ):
                pfv.c_pf_cs_coils_peak_ma[ic] = pfv.c_pf_cs_coil_flat_top_ma[ic]

            # End of flat-top, t = t_precharge + t_current_ramp_up + t_fusion_ramp + t_burn
            if (
                abs(pfv.c_pf_cs_coil_pulse_end_ma[ic])
                >= abs(pfv.c_pf_cs_coil_pulse_end_ma[ic])
            ) and (
                abs(pfv.c_pf_cs_coil_pulse_end_ma[ic])
                >= abs(pfv.c_pf_cs_coil_flat_top_ma[ic])
            ):
                pfv.c_pf_cs_coils_peak_ma[ic] = pfv.c_pf_cs_coil_pulse_end_ma[ic]

            # Set normalized current waveforms
            pfv.waves[ic, 0] = 0.0e0
            pfv.waves[ic, 1] = (
                pfv.c_pf_cs_coil_pulse_start_ma[ic] / pfv.c_pf_cs_coils_peak_ma[ic]
            )
            pfv.waves[ic, 2] = (
                pfv.c_pf_cs_coil_flat_top_ma[ic] / pfv.c_pf_cs_coils_peak_ma[ic]
            )
            pfv.waves[ic, 3] = (
                pfv.c_pf_cs_coil_flat_top_ma[ic] / pfv.c_pf_cs_coils_peak_ma[ic]
            )
            pfv.waves[ic, 4] = (
                pfv.c_pf_cs_coil_pulse_end_ma[ic] / pfv.c_pf_cs_coils_peak_ma[ic]
            )
            pfv.waves[ic, 5] = 0.0e0

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
            # ioheof = bv.hmax * pfv.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfv.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # Copper area = (area of coil - area of steel)*(1- void fraction)*
            # (fraction of copper in strands)
            # rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.f_a_cs_void) * pfv.fcuohsu)

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
            # ioheof = bv.hmax * pfv.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfv.j_cs_flat_top_end

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
            # ioheof = bv.hmax * pfv.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfv.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.f_a_cs_void) * pfv.fcuohsu)

        elif isumat == 9:
            # Hazelton experimental data + Zhai conceptual model for REBCO
            bc20m = 138
            tc0m = 92
            j_crit_sc, _, _ = superconductors.hijc_rebco(
                thelium,
                bmax,
                bc20m,
                tc0m,
                rebco_variables.tape_width,
                rebco_variables.rebco_thickness,
                rebco_variables.tape_thickness,
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0e0 - fcu) * (1.0e0 - fhe)

            # The CS coil current at EOF
            # ioheof = bv.hmax * pfv.f_z_cs_tf_internal * bv.dr_cs * 2.0 * pfv.j_cs_flat_top_end
            # The CS coil current/copper area calculation for quench protection
            # rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.f_a_cs_void) * pfv.fcuohsu)

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
                superconductors.current_density_margin,
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

        xc[i] = 0.5 * RMU0 * sd * ((2.0 - s) * xk - 2.0 * xe)

        #  Radial, vertical fields

        brx = (
            RMU0
            * cc[i]
            * dz
            / (2 * np.pi * rp * sd)
            * (-xk + (rc[i] ** 2 + rp**2 + zs) / (dr**2 + zs) * xe)
        )
        bzx = (
            RMU0
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
    rcls,
    zcls,
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
    :param rcls: coords R(i,j), Z(i,j) of coil j in group i (m)
    :type rcls: numpy.ndarray
    :param zcls: coords R(i,j), Z(i,j) of coil j in group i (m)
    :type zcls: numpy.ndarray
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
                rcls[j, :nc], zcls[j, :nc], cc[:nc], rpts[i], zpts[i]
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


def init_pfcoil_module():
    pf.first_call = True
    pf.cslimit = False
    pf.nef = 0.0
    pf.nfxf = 0.0
    pf.ricpf = 0.0
    pf.ssq0 = 0.0
    pf.sig_axial = 0.0
    pf.sig_hoop = 0.0
    pf.axial_force = 0
    pf.rfxf[:] = 0.0
    pf.zfxf[:] = 0.0
    pf.cfxf[:] = 0.0
    pf.xind[:] = 0.0
    pf.rcls[:] = 0.0
    pf.zcls[:] = 0.0
    pf.ccls[:] = 0.0
    pf.ccl0[:] = 0.0
    pf.bpf2[:] = 0.0
    pf.vsdum[:] = 0.0


def init_pfcoil_variables():
    """Initialise the PF coil variables"""
    pfv.alfapf = 5e-10
    pfv.alstroh = 4.0e8
    pfv.i_cs_stress = 0
    pfv.a_cs_poloidal = 0.0
    pfv.a_cs_turn = 0.0
    pfv.awpoh = 0.0
    pfv.b_cs_peak_flat_top_end = 0.0
    pfv.b_cs_peak_pulse_start = 0.0
    pfv.b_pf_coil_peak[:] = 0.0
    pfv.ccl0_ma[:] = 0.0
    pfv.ccls_ma[:] = 0.0
    pfv.j_cs_pulse_start = 0.0
    pfv.j_cs_flat_top_end = 1.85e7
    pfv.c_pf_coil_turn[:] = 0.0
    pfv.c_pf_coil_turn_peak_input[:] = 4.0e4
    pfv.c_pf_cs_coil_pulse_start_ma[:] = 0.0
    pfv.c_pf_cs_coil_flat_top_ma[:] = 0.0
    pfv.c_pf_cs_coil_pulse_end_ma[:] = 0.0
    pfv.etapsu = 0.9
    pfv.f_j_cs_start_end_flat_top = 0.0
    pfv.f_j_cs_start_pulse_end_flat_top = 0.9
    pfv.fcuohsu = 0.7
    pfv.fcupfsu = 0.69
    pfv.fvs_cs_pf_total_ramp = 1.0
    pfv.i_pf_location = [2, 2, 3, 0, 0, 0, 0, 0, 0, 0]
    pfv.i_pf_conductor = 0
    pfv.itr_sum = 0.0
    pfv.i_cs_superconductor = 1
    pfv.i_pf_superconductor = 1
    pfv.j_crit_str_cs = 0.0
    pfv.j_crit_str_pf = 0.0
    pfv.i_pf_current = 1
    pfv.i_sup_pf_shape = 0
    pfv.j_cs_conductor_critical_pulse_start = 0.0
    pfv.j_cs_conductor_critical_flat_top_end = 0.0
    pfv.jcableoh_bop = 0.0
    pfv.jcableoh_eof = 0.0
    pfv.n_pf_cs_plasma_circuits = 0
    pfv.n_pf_coils_in_group = [1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    pfv.nfxfh = 7
    pfv.n_pf_coil_groups = 3
    pfv.n_cs_pf_coils = 0
    pfv.f_z_cs_tf_internal = 0.71
    pfv.f_a_cs_steel = 0.5
    pfv.pf_current_safety_factor = 1.0
    pfv.pfcaseth[:] = 0.0
    pfv.rho_pf_coil = 2.5e-8
    pfv.rhopfbus = 3.93e-8
    pfv.m_pf_coil_max = 0.0
    pfv.r_pf_coil_outer_max = 0.0
    pfv.p_pf_electric_supplies_mw = 0.0
    pfv.p_cs_resistive_flat_top = 0.0
    pfv.p_pf_coil_resistive_total_flat_top = 0.0
    pfv.r_pf_coil_inner[:] = 0.0
    pfv.r_pf_coil_outer[:] = 0.0
    pfv.c_pf_cs_coils_peak_ma[:] = 0.0
    pfv.j_pf_coil_wp_peak[:] = 3.0e7
    pfv.j_cs_critical_flat_top_end = 0.0
    pfv.j_cs_critical_pulse_start = 0.0
    pfv.j_pf_wp_critical[:] = 0.0
    pfv.r_cs_middle = 0.0
    pfv.routr = 1.5
    pfv.r_pf_coil_middle[:] = 0.0
    pfv.rpf1 = 0.0
    pfv.rpf2 = -1.63
    pfv.rref[:] = 7.0
    pfv.s_shear_cs_peak = 0.0
    pfv.sigpfcalw = 500.0
    pfv.sigpfcf = 1.0
    pfv.ind_pf_cs_plasma_mutual[:] = 0.0
    pfv.temp_cs_margin = 0.0
    pfv.n_pf_coil_turns[:] = 0.0
    pfv.f_a_pf_coil_void[:] = 0.3
    pfv.f_a_cs_void = 0.3
    pfv.vs_cs_pf_total_burn = 0.0
    pfv.vs_pf_coils_total_burn = 0.0
    pfv.vs_pf_coils_total_ramp = 0.0
    pfv.vs_pf_coils_total_pulse = 0.0
    pfv.vs_cs_total_pulse = 0.0
    pfv.vs_cs_burn = 0.0
    pfv.vs_cs_ramp = 0.0
    pfv.vs_cs_pf_total_ramp = 0.0
    pfv.vs_cs_pf_total_pulse = 0.0
    pfv.waves[:] = 0.0
    pfv.m_pf_coil_conductor_total = 0.0
    pfv.m_pf_coil_structure_total = 0.0
    pfv.m_pf_coil_conductor[:] = 0.0
    pfv.m_pf_coil_structure[:] = 0.0
    pfv.z_pf_coil_upper[:] = 0.0
    pfv.z_pf_coil_lower[:] = 0.0
    pfv.z_pf_coil_middle[:] = 0.0
    pfv.zref = [3.6, 1.2, 2.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    pfv.b_cs_limit_max = 13.0
    pfv.fb_cs_limit_max = 1.0
    pfv.ld_ratio_cst = 70.0 / 22.0
    pfv.dr_cs_turn = 0.0
    pfv.dz_cs_turn = 0.0
    pfv.radius_cs_turn_cable_space = 0.0
    pfv.r_out_cst = 3.0e-3
