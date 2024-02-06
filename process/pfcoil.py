from process.fortran import pfcoil_module as pf
from process.fortran import pfcoil_variables as pfv
from process.fortran import times_variables as tv
from process.fortran import error_handling as eh
from process.fortran import build_variables as bv
from process.fortran import physics_variables as pv
from process.fortran import tfcoil_variables as tfv
from process.fortran import fwbs_variables as fwbsv
from process.fortran import constants
from process.fortran import cs_fatigue_variables as csfv
from process.fortran import maths_library as ml
from process.fortran import process_output as op
from process.fortran import numerics
from process.fortran import rebco_variables as rcv
from process.fortran import constraint_variables as ctv
from process import maths_library as pml
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from process import fortran as ft
import process.superconductors as superconductors
import math
import numpy as np
import numba
import logging

logger = logging.getLogger(__name__)

RMU0 = ft.constants.rmu0


class PFCoil:
    """Calculate poloidal field coil system parameters."""

    def __init__(self, cs_fatigue) -> None:
        """Initialise Fortran module variables."""
        self.outfile = ft.constants.nout  # output file unit
        self.mfile = ft.constants.mfile  # mfile file unit
        pf.init_pfcoil_module()
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        lrow1 = 2 * pfv.nptsmx + pfv.ngrpmx
        lcol1 = pfv.ngrpmx

        pcls0 = np.zeros(pfv.ngrpmx, dtype=int)
        ncls0 = np.zeros(pfv.ngrpmx + 2, dtype=int)

        pf.rcls0, pf.zcls0 = np.zeros((2, pfv.ngrpmx, pfv.nclsmx), order="F")
        pf.ccls0 = np.zeros(int(pfv.ngrpmx / 2))
        sigma, work2 = np.zeros((2, pfv.ngrpmx))
        rc, zc, cc, xc = np.zeros((4, pfv.nclsmx))
        brin, bzin, rpts, zpts = np.zeros((4, pfv.nptsmx))
        bfix, bvec = np.zeros((2, lrow1))
        gmat, umat, vmat = np.zeros((3, lrow1, lcol1), order="F")
        signn = np.zeros(2)
        aturn = np.zeros(pfv.ngc2)

        # Toggle switch for ipfloc()=2 coils above/below midplane
        top_bottom = 1

        # Set up the number of PF coils including the Central Solenoid (nohc),
        # and the number of PF circuits including the plasma (ncirt)
        if pfv.ngrp > pfv.ngrpmx:
            eh.idiags[0] = pfv.ngrp
            eh.idiags[1] = pfv.ngrpmx
            eh.report_error(64)

        # Total the number of PF coils in all groups, and check that none
        # exceeds the limit
        pfv.nohc = 0
        for i in range(pfv.ngrp):
            if pfv.ncls[i] > pfv.nclsmx:
                eh.idiags[0] = i
                eh.idiags[1] = pfv.ncls[i]
                eh.idiags[2] = pfv.nclsmx
                eh.report_error(65)

            pfv.nohc = pfv.nohc + pfv.ncls[i]

        # Add one if an Central Solenoid is present, and make an extra group
        if bv.iohcl != 0:
            pfv.nohc = pfv.nohc + 1
            pfv.ncls[pfv.ngrp] = 1

        # Add one for the plasma
        pfv.ncirt = pfv.nohc + 1

        # Overall current density in the Central Solenoid at beginning of pulse
        pfv.cohbop = pfv.coheof * pfv.fcohbop

        # Set up array of times
        tv.tim[0] = 0.0e0
        tv.tim[1] = tv.tramp
        tv.tim[2] = tv.tim[1] + tv.tohs
        tv.tim[3] = tv.tim[2] + tv.theat
        tv.tim[4] = tv.tim[3] + tv.tburn
        tv.tim[5] = tv.tim[4] + tv.tqnch

        # Set up call to MHD scaling routine for coil currents.
        # First break up Central Solenoid solenoid into 'filaments'

        # Central Solenoid radius
        pfv.rohc = bv.bore + 0.5e0 * bv.ohcth

        # nfxf is the total no of filaments into which the Central Solenoid is split,
        # if present
        if bv.iohcl == 0:
            pf.nfxf = 0
            ioheof = 0.0e0
        else:
            pf.nfxf = 2 * pfv.nfxfh

            # total Central Solenoid current at EOF
            ioheof = -bv.hmax * pfv.ohhghf * bv.ohcth * 2.0e0 * pfv.coheof

            if pf.nfxf > pfv.nfixmx:
                eh.idiags[0] = pf.nfxf
                eh.idiags[1] = pfv.nfixmx
                eh.report_error(66)

            # Symmetric up/down Central Solenoid : Find (R,Z) and current of each filament at BOP

            for nng in range(pfv.nfxfh):
                pf.rfxf[nng] = pfv.rohc
                pf.rfxf[nng + pfv.nfxfh] = pf.rfxf[nng]
                pf.zfxf[nng] = bv.hmax * pfv.ohhghf / pfv.nfxfh * ((nng + 1) - 0.5e0)
                pf.zfxf[nng + pfv.nfxfh] = -pf.zfxf[nng]
                pf.cfxf[nng] = -ioheof / pf.nfxf * pfv.fcohbop
                pf.cfxf[nng + pfv.nfxfh] = pf.cfxf[nng]

        # Scale PF coil locations
        signn[0] = 1.0e0
        signn[1] = -1.0e0
        pf.rclsnorm = bv.r_tf_outboard_mid + 0.5e0 * bv.tfthko + pfv.routr

        # Place the PF coils:

        # N.B. Problems here if k=ncls(group) is greater than 2.
        for j in range(pfv.ngrp):
            if pfv.ipfloc[j] == 1:
                # PF coil is stacked on top of the Central Solenoid
                for k in pfv.ncls[j]:
                    pf.rcls[j, k] = pfv.rohc + pfv.rpf1

                    # Z coordinate of coil enforced so as not
                    # to occupy the same space as the Central Solenoid
                    pf.zcls[j, k] = signn[k] * (
                        bv.hmax * pfv.ohhghf
                        + 0.1e0
                        + 0.5e0 * (bv.hmax * (1.0e0 - pfv.ohhghf) + bv.tfcth + 0.1e0)
                    )

            elif pfv.ipfloc[j] == 2:
                # PF coil is on top of the TF coil
                for k in range(pfv.ncls[j]):
                    pf.rcls[j, k] = pv.rmajor + pfv.rpf2 * pv.triang * pv.rminor
                    if pv.itart == 1 and pv.itartpf == 0:
                        pf.zcls[j, k] = (bv.hmax - pfv.zref[j]) * signn[k]
                    else:
                        # pf.zcls(j,k) = (bv.hmax + bv.tfcth + 0.86e0) * signn(k)
                        if top_bottom == 1:  # this coil is above midplane
                            pf.zcls[j, k] = bv.hpfu + 0.86e0
                            top_bottom = -1
                        else:  # this coil is below midplane
                            pf.zcls[j, k] = -1.0e0 * (
                                bv.hpfu - 2.0e0 * bv.hpfdif + 0.86e0
                            )
                            top_bottom = 1

            elif pfv.ipfloc[j] == 3:
                # PF coil is radially outside the TF coil
                for k in range(pfv.ncls[j]):
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

            elif pfv.ipfloc[j] == 4:
                # PF coil is in general location
                # See issue 1418
                # https://git.ccfe.ac.uk/process/process/-/issues/1418
                for k in range(pfv.ncls[j]):
                    pf.zcls[j, k] = pv.rminor * pfv.zref[j] * signn[k]
                    pf.rcls[j, k] = pv.rminor * pfv.rref[j] + pv.rmajor

            else:
                eh.idiags[0] = j
                eh.idiags[1] = pfv.ipfloc[j]
                eh.report_error(67)

        # Allocate current to the PF coils:
        # "Flux swing coils" participate in cancellation of the CS
        # field during a flux swing. "Equilibrium coils" are varied
        # to create the equilibrium field, targeting the correct
        # vertical field
        # As implemented, all coils are flux swing coils
        # As implemented, Location 3 and 4 coils are equilibrium
        # coils.

        # Flux swing coils:
        if pfv.cohbop != 0.0e0:
            # Find currents for plasma initiation to null field across plasma
            npts = 32  # Number of test points across plasma midplane
            if npts > pfv.nptsmx:
                eh.idiags[0] = npts
                eh.idiags[1] = pfv.nptsmx
                eh.report_error(68)

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
                pfv.ngrp,
                pfv.ncls,
                pf.rcls,
                pf.zcls,
                pfv.alfapf,
                bfix,
                gmat,
                bvec,
                rc,
                zc,
                cc,
                xc,
                umat,
                vmat,
                sigma,
                work2,
            )

        # Equilibrium coil currents determined by SVD targeting B
        if pfv.i_pf_current == 1:
            # Simple coil current scaling for STs (good only for A < about 1.8)
            # Bypasses SVD solver
            if pv.itart == 1 and pv.itartpf == 0:
                for i in range(pfv.ngrp):
                    if pfv.ipfloc[i] == 1:
                        # PF coil is stacked on top of the Central Solenoid
                        pf.ccls[i] = 0.0e0
                        eh.idiags[0] = i
                        eh.report_error(69)

                    elif pfv.ipfloc[i] == 2:
                        # PF coil is on top of the TF coil
                        pf.ccls[i] = 0.3e0 * pv.aspect**1.6e0 * pv.plascur

                    elif pfv.ipfloc[i] == 3:
                        # PF coil is radially outside the TF coil
                        pf.ccls[i] = -0.4e0 * pv.plascur

                    else:
                        eh.idiags[0] = i
                        eh.idiags[1] = pfv.ipfloc[i]
                        eh.report_error(70)

                # Vertical field (T)
                pv.bvert = (
                    -1.0e-7
                    * pv.plascur
                    / pv.rmajor
                    * (
                        math.log(8.0e0 * pv.aspect)
                        + pv.betap
                        + (pv.rli / 2.0e0)
                        - 1.5e0
                    )
                )

            else:
                # Conventional aspect ratio scaling
                nfxf0 = 0
                ngrp0 = 0
                nocoil = 0
                for i in range(pfv.ngrp):
                    if pfv.ipfloc[i] == 1:
                        # PF coil is stacked on top of the Central Solenoid
                        # This coil is to balance Central Solenoid flux and should not be involved
                        # in equilibrium calculation -- RK 07/12
                        pf.ccls[i] = 0.0e0
                        nfxf0 = nfxf0 + pfv.ncls[i]
                        for ccount in range(pfv.ncls[i]):
                            pf.rfxf[nocoil] = pf.rcls[i, ccount]
                            pf.zfxf[nocoil] = pf.zcls[i, ccount]
                            pf.cfxf[nocoil] = pf.ccls[i]
                            nocoil = nocoil + 1

                    elif pfv.ipfloc[i] == 2:
                        # PF coil is on top of the TF coil; divertor coil
                        # This is a fixed current for this calculation -- RK 07/12

                        pf.ccls[i] = (
                            pv.plascur
                            * 2.0e0
                            * (1.0e0 - (pv.kappa * pv.rminor) / abs(pf.zcls[i, 0]))
                        )
                        nfxf0 = nfxf0 + pfv.ncls[i]
                        for ccount in range(pfv.ncls[i]):
                            pf.rfxf[nocoil] = pf.rcls[i, ccount]
                            pf.zfxf[nocoil] = pf.zcls[i, ccount]
                            pf.cfxf[nocoil] = pf.ccls[i]
                            nocoil = nocoil + 1

                    elif pfv.ipfloc[i] == 3:
                        # PF coil is radially outside the TF coil
                        # This is an equilibrium coil, current must be solved for

                        pcls0[ngrp0] = i + 1
                        ngrp0 = ngrp0 + 1

                    elif pfv.ipfloc[i] == 4:
                        # PF coil is generally placed
                        # See issue 1418
                        # https://git.ccfe.ac.uk/process/process/-/issues/1418
                        # This is an equilibrium coil, current must be solved for

                        pcls0[ngrp0] = i + 1
                        ngrp0 = ngrp0 + 1

                    else:
                        eh.idiags[0] = i
                        eh.idiags[1] = pfv.ipfloc[i]
                        eh.report_error(70)

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

                # Added pv.rli term correctly -- RK 07/12

                bzin[0] = (
                    -1.0e-7
                    * pv.plascur
                    / pv.rmajor
                    * (
                        math.log(8.0e0 * pv.aspect)
                        + pv.betap
                        + (pv.rli / 2.0e0)
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
                    rc,
                    zc,
                    cc,
                    xc,
                    umat,
                    vmat,
                    sigma,
                    work2,
                )

                for ccount in range(ngrp0):
                    pf.ccls[pcls0[ccount] - 1] = pf.ccls0[ccount]

        # Flux swing from vertical field

        # If this is the first visit to the routine the inductance matrix
        # sxlg and the turns array have not yet been calculated, so we set
        # them to (very) approximate values to avoid strange behaviour...
        if pf.first_call:
            pfv.sxlg[:, :] = 1.0e0
            pfv.turns[:] = 100.0e0
            pf.first_call = False

        pfflux = 0.0e0
        nocoil = 0
        for ccount in range(pfv.ngrp):
            for i in range(pfv.ncls[ccount]):
                pfflux = pfflux + (
                    pf.ccls[ccount]
                    * pfv.sxlg[nocoil, pfv.ncirt - 1]
                    / pfv.turns[nocoil]
                )
                nocoil = nocoil + 1

        # Flux swing required from CS coil
        csflux = -(pv.vsres + pv.vsind) - pfflux

        if bv.iohcl == 1:
            # Required current change in CS coil

            # Proposed new calculation...
            # dics = csflux / sxlg(nohc,ncirt)
            # BUT... sxlg(nohc,ncirt) is around 2000 times ddics below...

            ddics = (
                4.0e-7
                * constants.pi
                * constants.pi
                * (
                    (bv.bore * bv.bore)
                    + (bv.ohcth * bv.ohcth) / 6.0e0
                    + (bv.ohcth * bv.bore) / 2.0e0
                )
                / (bv.hmax * pfv.ohhghf * 2.0e0)
            )
            dics = csflux / ddics

            pfv.fcohbof = ((-ioheof * pfv.fcohbop) + dics) / ioheof
            pfv.fcohbof = min(pfv.fcohbof, 1.0e0)  # constrains abs(fcohbof) <= 1.0;
            pfv.fcohbof = max(pfv.fcohbof, -1.0e0)  # probably un-necessary

        else:
            dics = 0.0e0
            pfv.fcohbof = 1.0e0
            eh.report_error(71)

        # Split groups of coils into one set containing ncl coils
        ncl = 0
        for nng in range(pfv.ngrp):
            for ng2 in range(pfv.ncls[nng]):
                pfv.rpf[ncl] = pf.rcls[nng, ng2]
                pfv.zpf[ncl] = pf.zcls[nng, ng2]

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

                # Beginning of pulse: t = tv.tramp
                pfv.curpfs[ncl] = 1.0e-6 * pf.ccl0[nng]

                # Beginning of flat-top: t = tv.tramp+tv.tohs
                pfv.curpff[ncl] = 1.0e-6 * (
                    pf.ccls[nng] - (pf.ccl0[nng] * pfv.fcohbof / pfv.fcohbop)
                )

                # End of flat-top: t = tv.tramp+tv.tohs+tv.theat+tv.tburn
                pfv.curpfb[ncl] = 1.0e-6 * (
                    pf.ccls[nng] - (pf.ccl0[nng] * (1.0e0 / pfv.fcohbop))
                )

                ncl = ncl + 1

        # Current in Central Solenoid as a function of time
        # N.B. If the Central Solenoid is not present then ioheof is zero.
        pfv.curpfs[ncl] = -1.0e-6 * ioheof * pfv.fcohbop
        pfv.curpff[ncl] = 1.0e-6 * ioheof * pfv.fcohbof
        pfv.curpfb[ncl] = 1.0e-6 * ioheof

        # Set up coil current waveforms, normalised to the peak current in
        # each coil
        self.waveform()  # sets ric(), waves()

        # Calculate PF coil geometry, current and number of turns
        # Dimensions are those of the winding pack, and exclude
        # the steel supporting case
        i = 0
        pfv.pfrmax = 0.0e0

        dz = 0

        for ii in range(pfv.ngrp):
            for ij in range(pfv.ncls[ii]):
                if pfv.ipfloc[ii] == 1:
                    # PF coil is stacked on top of the Central Solenoid
                    dx = 0.5e0 * bv.ohcth
                    dz = 0.5e0 * (
                        bv.hmax * (1.0e0 - pfv.ohhghf) + bv.tfcth + 0.1e0
                    )  # ???
                    area = 4.0e0 * dx * dz * pfv.pf_current_safety_factor

                    # Number of turns
                    # CPTDIN[i] is the current per turn (input)
                    pfv.turns[i] = abs((pfv.ric[i] * 1.0e6) / pfv.cptdin[i])
                    aturn[i] = area / pfv.turns[i]

                    # Actual winding pack current density
                    pfv.rjconpf[i] = 1.0e6 * abs(pfv.ric[i]) / area

                    # Location of edges of each coil:
                    # ra = inner radius, rb = outer radius
                    # zl = 'lower' edge z (i.e. edge nearer to midplane)
                    # zh = 'upper' edge z (i.e. edge further from midplane)
                    pfv.ra[i] = pfv.rpf[i] - dx
                    pfv.rb[i] = pfv.rpf[i] + dx

                    pfv.zl[i] = pfv.zpf[i] - dz
                    if pfv.zpf[i] < 0.0e0:
                        pfv.zl[i] = pfv.zpf[i] + dz

                    pfv.zh[i] = pfv.zpf[i] + dz

                    if pfv.zpf[i] < 0.0e0:
                        pfv.zh[i] = pfv.zpf[i] - dz

                else:
                    # Other coils. N.B. Current density RJCONPF[i] is defined in
                    # routine INITIAL for these coils.
                    area = (
                        abs(pfv.ric[i] * 1.0e6 / pfv.rjconpf[i])
                        * pfv.pf_current_safety_factor
                    )

                    pfv.turns[i] = abs((pfv.ric[i] * 1.0e6) / pfv.cptdin[i])
                    aturn[i] = area / pfv.turns[i]

                    dx = 0.5e0 * math.sqrt(area)  # square cross-section

                    pfv.ra[i] = pfv.rpf[i] - dx
                    pfv.rb[i] = pfv.rpf[i] + dx

                    pfv.zl[i] = pfv.zpf[i] - dx
                    if pfv.zpf[i] < 0.0e0:
                        pfv.zl[i] = pfv.zpf[i] + dx

                    pfv.zh[i] = pfv.zpf[i] + dx
                    if pfv.zpf[i] < 0.0e0:
                        pfv.zh[i] = pfv.zpf[i] - dx

                # Outside radius of largest PF coil (m)
                pfv.pfrmax = max(pfv.pfrmax, pfv.rb[i])

                i = i + 1

        # Calculate peak field, allowable current density, resistive
        # power losses and volumes and weights for each PF coil
        i = 0
        it = 0
        pfv.powpfres = 0.0e0
        pfv.pfmmax = 0.0e0

        for ii in range(pfv.ngrp):
            iii = ii
            for ij in range(pfv.ncls[ii]):
                # Peak field

                if ij == 0:
                    # Index args +1ed
                    bri, bro, bzi, bzo = self.peakb(
                        i + 1, iii + 1, it
                    )  # returns bpf, bpf2

                # Allowable current density (for superconducting coils)

                if pfv.ipfres == 0:
                    bmax = max(abs(pfv.bpf[i]), abs(pf.bpf2[i]))
                    pfv.rjpfalw[i], jstrand, jsc, tmarg = self.superconpf(
                        bmax,
                        pfv.vf[i],
                        pfv.fcupfsu,
                        pfv.rjconpf[i],
                        pfv.isumatpf,
                        tfv.fhts,
                        tfv.str_pf_con_res,
                        tfv.tftmp,
                        tfv.bcritsc,
                        tfv.tcritsc,
                    )

                # Length of conductor

                rll = 2.0e0 * constants.pi * pfv.rpf[i] * pfv.turns[i]

                # Resistive coils

                if pfv.ipfres == 1:
                    # Coil resistance (vf is the void fraction)

                    respf = pfv.pfclres * rll / (aturn[i] * (1.0e0 - pfv.vf[i]))

                    # Sum resistive power losses

                    pfv.powpfres = (
                        pfv.powpfres
                        + respf * (1.0e6 * pfv.curpfb[i] / pfv.turns[i]) ** 2
                    )

                # Winding pack volume

                volpf = aturn[i] * rll

                # Conductor weight (vf is the void fraction)

                if pfv.ipfres == 0:
                    pfv.wtc[i] = (
                        volpf * tfv.dcond[pfv.isumatpf - 1] * (1.0e0 - pfv.vf[i])
                    )
                else:
                    pfv.wtc[i] = volpf * constants.dcopper * (1.0e0 - pfv.vf[i])

                # (J x B) force on coil

                forcepf = (
                    0.5e6 * (pfv.bpf[i] + pf.bpf2[i]) * abs(pfv.ric[i]) * pfv.rpf[i]
                )

                # Stress ==> cross-sectional area of supporting steel to use

                if pfv.ipfres == 0:
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
                        pfv.rb[i] - pfv.ra[i] + abs(pfv.zh[i] - pfv.zl[i])
                    )  # dr + dz
                    pfv.pfcaseth[i] = 0.25e0 * (
                        -drpdz + math.sqrt(drpdz * drpdz + 4.0e0 * areaspf)
                    )

                else:
                    areaspf = 0.0e0  # Resistive coil - no steel needed
                    pfv.pfcaseth[i] = 0.0e0

                # Weight of steel case

                pfv.wts[i] = areaspf * 2.0e0 * constants.pi * pfv.rpf[i] * fwbsv.denstl

                # Mass of heaviest PF coil (tonnes)

                pfv.pfmmax = max(
                    pfv.pfmmax,
                    (1.0e-3 * (pfv.wtc[i] + pfv.wts[i])),
                )
                i = i + 1

        # Find sum of current x turns x radius for all coils for 2015 costs model
        c = 0
        pfv.itr_sum = 0.0e0
        for m in range(pfv.ngrp):
            for n in range(pfv.ncls[m]):
                pfv.itr_sum = pfv.itr_sum + (pfv.rpf[c] * pfv.turns[c] * pfv.cptdin[c])
                c = c + 1

        pfv.itr_sum = pfv.itr_sum + (
            (bv.bore + 0.5 * bv.ohcth)
            * pfv.turns[pfv.nohc - 1]
            * pfv.cptdin[pfv.nohc - 1]
        )

        # Find Central Solenoid information
        if bv.iohcl != 0:
            self.ohcalc()

        # Summation of weights and current
        pfv.whtpf = 0.0e0
        pfv.whtpfs = 0.0e0
        pf.ricpf = 0.0e0

        for i in range(pfv.nohc):
            pfv.whtpf = pfv.whtpf + pfv.wtc[i]
            pfv.whtpfs = pfv.whtpfs + pfv.wts[i]
            pf.ricpf = pf.ricpf + abs(pfv.ric[i])

        # Plasma size and shape
        pfv.zh[pfv.nohc] = pv.rminor * pv.kappa
        pfv.zl[pfv.nohc] = -pv.rminor * pv.kappa
        pfv.ra[pfv.nohc] = pv.rmajor - pv.rminor
        pfv.rb[pfv.nohc] = pv.rmajor + pv.rminor
        pfv.turns[pfv.nohc] = 1.0e0

        # Generate coil currents as a function of time using
        # user-provided waveforms etc. (cptdin, fcohbop, fcohbof)
        for k in range(6):  # time points
            for i in range(pfv.ncirt - 1):
                pfv.cpt[i, k] = pfv.waves[i, k] * math.copysign(
                    pfv.cptdin[i], pfv.ric[i]
                )

        # Plasma wave form
        pfv.cpt[pfv.ncirt - 1, 0] = 0.0e0
        pfv.cpt[pfv.ncirt - 1, 1] = 0.0e0
        pfv.cpt[pfv.ncirt - 1, 2] = pv.plascur
        pfv.cpt[pfv.ncirt - 1, 3] = pv.plascur
        pfv.cpt[pfv.ncirt - 1, 4] = pv.plascur
        pfv.cpt[pfv.ncirt - 1, 5] = 0.0e0

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
        ngrp,
        ncls,
        rcls,
        zcls,
        alfa,
        bfix,
        gmat,
        bvec,
        rc,
        zc,
        cc,
        xc,
        umat,
        vmat,
        sigma,
        work2,
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
        :param ngrp: number of coil groups, where all coils in a group have the
        same current, <= ngrpmx
        :type ngrp: int
        :param ncls: number of coils in each group, each value <= nclsmx
        :type ncls: np.ndarray
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
            int(ngrp),
            ncls,
            rcls,
            zcls,
            alfa,
            bfix,
            int(pfv.nclsmx),
        )

        # Solve matrix equation
        ccls, umat, vmat, sigma, work2 = self.solv(pfv.ngrpmx, ngrp, nrws, gmat, bvec)

        # Calculate the norm of the residual vectors
        brssq, brnrm, bzssq, bznrm, ssq = rsid(
            npts, brin, bzin, nfix, int(ngrp), ccls, bfix, gmat
        )

        return ssq, ccls

    def tf_pf_collision_detector(self):
        #  Collision test between TF and PF coils for picture frame TF
        #  See issue 1612
        #  https://git.ccfe.ac.uk/process/process/-/issues/1612

        if tfv.i_tf_shape == 2:
            pf_tf_collision = 0

            for i in range(pfv.ngrp):
                for ii in range(pfv.ngrp):
                    for ij in range(pfv.ncls[ii]):
                        if pf.rcls[ii, ij] <= (  # Outboard TF coil collision
                            pf.rclsnorm - pfv.routr + pfv.rpf[i]
                        ) and pf.rcls[ii, ij] >= (
                            bv.r_tf_outboard_mid - (0.5 * bv.tfthko) - pfv.rpf[i]
                        ):
                            pf_tf_collision += 1
                        if pf.rcls[ii, ij] <= (  # Inboard TF coil collision
                            bv.bore
                            + bv.ohcth
                            + bv.precomp
                            + bv.gapoh
                            + bv.tfcth
                            + pfv.rpf[i]
                        ) and pf.rcls[ii, ij] >= (
                            bv.bore + bv.ohcth + bv.precomp + bv.gapoh - pfv.rpf[i]
                        ):
                            pf_tf_collision += 1
                        if (  # Vertical TF coil collision
                            abs(pf.zcls[ii, ij]) <= bv.hpfu + pfv.rpf[i]
                            and abs(pf.zcls[ii, ij])
                            >= bv.hpfu - (0.5 * bv.tfthko) - pfv.rpf[i]
                        ):
                            pf_tf_collision += 1

                        if pf_tf_collision >= 1:
                            eh.report_error(277)

    def solv(self, ngrpmx, ngrp, nrws, gmat, bvec):
        """Solve a matrix using singular value decomposition.

        This routine solves the matrix equation for calculating the
        currents in a group of ring coils.
        author: P J Knight, CCFE, Culham Science Centre
        author: D Strickler, ORNL
        author: J Galambos, ORNL
        author: P C Shipe, ORNL

        :param ngrpmx: maximum number of PF coil groups
        :type ngrpmx: int
        :param ngrp: number of coil groups, where all coils in a group have the
        same current, <= ngrpmx
        :type ngrp: int
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
        truth = True
        eps = 1.0e-10
        ccls = np.zeros(ngrpmx)

        sigma, umat, vmat, ierr, work2 = ml.svd(
            nrws, np.asfortranarray(gmat), truth, truth
        )

        for i in range(ngrp):
            work2[i] = 0.0e0
            for j in range(nrws):
                work2[i] = work2[i] + umat[j, i] * bvec[j]

        # Compute currents
        for i in range(ngrp):
            ccls[i] = 0.0e0
            zvec = 0.0e0
            for j in range(ngrp):
                if sigma[j] > eps:
                    zvec = work2[j] / sigma[j]

                ccls[i] = ccls[i] + vmat[i, j] * zvec

        return ccls, umat, vmat, sigma, work2

    def ohcalc(self):
        """Routine to perform calculations for the Central Solenoid solenoid.

        author: P J Knight, CCFE, Culham Science Centre
        This subroutine performs the calculations for the
        Central Solenoid solenoid coil.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        hohc = bv.hmax * pfv.ohhghf

        # Z coordinates of coil edges
        pfv.zh[pfv.nohc - 1] = hohc
        pfv.zl[pfv.nohc - 1] = -pfv.zh[pfv.nohc - 1]

        # (R,Z) coordinates of coil centre
        pfv.rpf[pfv.nohc - 1] = pfv.rohc
        pfv.zpf[pfv.nohc - 1] = 0.0e0

        # Radius of outer edge
        pfv.rb[pfv.nohc - 1] = pfv.rohc + 0.5e0 * bv.ohcth

        # Radius of inner edge
        pfv.ra[pfv.nohc - 1] = pfv.rb[pfv.nohc - 1] - bv.ohcth

        # Total cross-sectional area
        pfv.areaoh = 2.0e0 * hohc * bv.ohcth

        # Maximum current (MA-turns) in central Solenoid, at either BOP or EOF
        if pfv.cohbop > pfv.coheof:
            sgn = 1.0e0
            pfv.ric[pfv.nohc - 1] = sgn * 1.0e-6 * pfv.cohbop * pfv.areaoh
        else:
            sgn = -1.0e0
            pfv.ric[pfv.nohc - 1] = sgn * 1.0e-6 * pfv.coheof * pfv.areaoh

        # Number of turns
        pfv.turns[pfv.nohc - 1] = (
            1.0e6 * abs(pfv.ric[pfv.nohc - 1]) / pfv.cptdin[pfv.nohc - 1]
        )

        # Turn vertical cross-sectionnal area
        pfv.a_oh_turn = pfv.areaoh / pfv.turns[pfv.nohc - 1]

        # Depth/width of cs turn conduit
        pfv.d_cond_cst = (pfv.a_oh_turn / pfv.ld_ratio_cst) ** 0.5
        # length of cs turn conduit
        pfv.l_cond_cst = pfv.ld_ratio_cst * pfv.d_cond_cst
        # Radius of turn space = pfv.r_in_cst
        # Radius of curved outer corrner pfv.r_out_cst = 3mm from literature
        # pfv.ld_ratio_cst = 70 / 22 from literature
        p1_cst = ((pfv.l_cond_cst - pfv.d_cond_cst) / constants.pi) ** 2
        p2_cst = (
            (pfv.l_cond_cst * pfv.d_cond_cst)
            - (4 - constants.pi) * (pfv.r_out_cst**2)
            - (pfv.a_oh_turn * pfv.oh_steel_frac)
        ) / constants.pi
        # CS coil turn geometry calculation - stadium shape
        # Literature: https://doi.org/10.1016/j.fusengdes.2017.04.052
        pfv.r_in_cst = -((pfv.l_cond_cst - pfv.d_cond_cst) / constants.pi) + math.sqrt(
            p1_cst + p2_cst
        )
        # Thickness of steel conduit in cs turn
        csfv.t_structural_radial = (pfv.d_cond_cst / 2) - pfv.r_in_cst
        # In this model the vertical and radial have the same thickness
        csfv.t_structural_vertical = csfv.t_structural_radial
        # add a check for negative conduit thickness
        if csfv.t_structural_radial < 1.0e-3:
            csfv.t_structural_radial = 1.0e-3

        # Non-steel area void fraction for coolant
        pfv.vf[pfv.nohc - 1] = pfv.vfohc

        # Peak field at the End-Of-Flattop (EOF)
        # Occurs at inner edge of coil; bmaxoh2 and bzi are of opposite sign at EOF

        # Peak field due to central Solenoid itself
        bmaxoh2 = self.bfmax(
            pfv.coheof,
            pfv.ra[pfv.nohc - 1],
            pfv.rb[pfv.nohc - 1],
            hohc,
        )

        # Peak field due to other PF coils plus plasma
        timepoint = 5
        bri, bro, bzi, bzo = self.peakb(pfv.nohc, 99, timepoint)

        pfv.bmaxoh = abs(bzi - bmaxoh2)

        # Peak field on outboard side of central Solenoid
        # (self-field is assumed to be zero - long solenoid approximation)
        bohco = abs(bzo)

        # Peak field at the Beginning-Of-Pulse (BOP)
        # Occurs at inner edge of coil; bmaxoh0 and bzi are of same sign at BOP
        pfv.bmaxoh0 = self.bfmax(
            pfv.cohbop,
            pfv.ra[pfv.nohc - 1],
            pfv.rb[pfv.nohc - 1],
            hohc,
        )
        timepoint = 2
        bri, bro, bzi, bzo = self.peakb(pfv.nohc, 99, timepoint)

        pfv.bmaxoh0 = abs(pfv.bmaxoh0 + bzi)

        # Maximum field values
        pfv.bpf[pfv.nohc - 1] = max(pfv.bmaxoh, abs(pfv.bmaxoh0))
        pf.bpf2[pfv.nohc - 1] = max(bohco, abs(bzo))

        # Stress ==> cross-sectional area of supporting steel to use
        if pfv.ipfres == 0:
            # Superconducting coil

            # New calculation from M. N. Wilson for hoop stress
            pf.sig_hoop = self.hoop_stress(pfv.ra[pfv.nohc - 1])

            # New calculation from Y. Iwasa for axial stress
            pf.sig_axial, pf.axial_force = self.axial_stress()

            # Allowable (hoop) stress (Pa) alstroh
            # Now a user input
            # alstroh = min( (2.0e0*csytf/3.0e0), (0.5e0*csutf) )

            # Calculation of CS fatigue
            # this is only valid for pulsed reactor design
            if pv.facoh > 0.0e-4:
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
            areaspf = pfv.oh_steel_frac * pfv.areaoh

            if pfv.i_cs_stress == 1:
                pfv.s_tresca_oh = max(
                    abs(pf.sig_hoop - pf.sig_axial),
                    abs(pf.sig_axial - 0.0e0),
                    abs(0.0e0 - pf.sig_hoop),
                )
            else:
                pfv.s_tresca_oh = max(
                    abs(pf.sig_hoop - 0.0e0),
                    abs(0.0e0 - 0.0e0),
                    abs(0.0e0 - pf.sig_hoop),
                )

            # Thickness of hypothetical steel cylinders assumed to encase the CS along
            # its inside and outside edges; in reality, the steel is distributed
            # throughout the conductor
            pfv.pfcaseth[pfv.nohc - 1] = 0.25e0 * areaspf / hohc

        else:
            areaspf = 0.0e0  # Resistive Central Solenoid - no steel needed
            pfv.pfcaseth[pfv.nohc - 1] = 0.0e0

        # Weight of steel
        pfv.wts[pfv.nohc - 1] = (
            areaspf * 2.0e0 * constants.pi * pfv.rpf[pfv.nohc - 1] * fwbsv.denstl
        )

        # Non-steel cross-sectional area
        pfv.awpoh = pfv.areaoh - areaspf

        # Issue #97. Fudge to ensure awpoh is positive; result is continuous, smooth and
        # monotonically decreases

        da = 0.0001e0  # 1 cm^2
        if pfv.awpoh < da:
            pfv.awpoh = da * da / (2.0e0 * da - pfv.awpoh)

        # Weight of conductor in central Solenoid
        if pfv.ipfres == 0:
            pfv.wtc[pfv.nohc - 1] = (
                pfv.awpoh
                * (1.0e0 - pfv.vfohc)
                * 2.0e0
                * constants.pi
                * pfv.rpf[pfv.nohc - 1]
                * tfv.dcond[pfv.isumatoh - 1]
            )
        else:
            pfv.wtc[pfv.nohc - 1] = (
                pfv.awpoh
                * (1.0e0 - pfv.vfohc)
                * 2.0e0
                * constants.pi
                * pfv.rpf[pfv.nohc - 1]
                * constants.dcopper
            )

        if pfv.ipfres == 0:
            # Allowable coil overall current density at EOF
            # (superconducting coils only)

            (jcritwp, pfv.jstrandoh_eof, pfv.jscoh_eof, tmarg1,) = self.superconpf(
                pfv.bmaxoh,
                pfv.vfohc,
                pfv.fcuohsu,
                (abs(pfv.ric[pfv.nohc - 1]) / pfv.awpoh) * 1.0e6,
                pfv.isumatoh,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )

            pfv.rjohc = jcritwp * pfv.awpoh / pfv.areaoh

            # Allowable coil overall current density at BOP

            (jcritwp, pfv.jstrandoh_bop, pfv.jscoh_bop, tmarg2,) = self.superconpf(
                pfv.bmaxoh0,
                pfv.vfohc,
                pfv.fcuohsu,
                (abs(pfv.ric[pfv.nohc - 1]) / pfv.awpoh) * 1.0e6,
                pfv.isumatoh,
                tfv.fhts,
                tfv.str_cs_con_res,
                tfv.tftmp,
                tfv.bcritsc,
                tfv.tcritsc,
            )

            pfv.rjpfalw[pfv.nohc - 1] = jcritwp * pfv.awpoh / pfv.areaoh
            pfv.rjohc0 = pfv.rjpfalw[pfv.nohc - 1]

            pfv.tmargoh = min(tmarg1, tmarg2)

        else:
            # Resistive power losses (non-superconducting coil)

            pfv.powohres = (
                2.0e0
                * constants.pi
                * pfv.rohc
                * pfv.pfclres
                / (pfv.areaoh * (1.0e0 - pfv.vfohc))
                * (1.0e6 * pfv.ric[pfv.nohc - 1]) ** 2
            )
            pfv.powpfres = pfv.powpfres + pfv.powohres

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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if bv.iohcl != 0 and i == pfv.nohc:
            # Peak field is to be calculated at the Central Solenoid itself,
            # so exclude its own contribution; its self field is
            # dealt with externally using routine BFMAX
            kk = 0
        else:
            # Check different times for maximum current
            if abs(pfv.curpfs[i - 1] - pfv.ric[i - 1]) < 1.0e-12:
                it = 2
            elif abs(pfv.curpff[i - 1] - pfv.ric[i - 1]) < 1.0e-12:
                it = 4
            elif abs(pfv.curpfb[i - 1] - pfv.ric[i - 1]) < 1.0e-12:
                it = 5
            else:
                eh.idiags[0] = it
                eh.report_error(72)

            if bv.iohcl == 0:
                # No Central Solenoid
                kk = 0
            else:
                if pfv.cohbop > pfv.coheof:
                    sgn = 1.0e0
                else:
                    sgn = -1.0e0

                # Current in each filament representing part of the Central Solenoid
                for iohc in range(pf.nfxf):
                    pf.cfxf[iohc] = (
                        pfv.waves[pfv.nohc - 1, it - 1]
                        * pfv.coheof
                        * sgn
                        * bv.ohcth
                        * pfv.ohhghf
                        * bv.hmax
                        / pf.nfxf
                        * 2.0e0
                    )

                kk = pf.nfxf

        # Non-Central Solenoid coils' contributions
        jj = 0
        for iii in range(pfv.ngrp):
            for jjj in range(pfv.ncls[iii]):
                jj = jj + 1
                # Radius, z-coordinate and current for each coil
                if iii == ii - 1:
                    # Self field from coil (Lyle's Method)
                    kk = kk + 1

                    dzpf = pfv.zh[jj - 1] - pfv.zl[jj - 1]
                    pf.rfxf[kk - 1] = pfv.rpf[jj - 1]
                    pf.zfxf[kk - 1] = pfv.zpf[jj - 1] + dzpf * 0.125e0
                    pf.cfxf[kk - 1] = (
                        pfv.ric[jj - 1] * pfv.waves[jj - 1, it - 1] * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.rpf[jj - 1]
                    pf.zfxf[kk - 1] = pfv.zpf[jj - 1] + dzpf * 0.375e0
                    pf.cfxf[kk - 1] = (
                        pfv.ric[jj - 1] * pfv.waves[jj - 1, it - 1] * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.rpf[jj - 1]
                    pf.zfxf[kk - 1] = pfv.zpf[jj - 1] - dzpf * 0.125e0
                    pf.cfxf[kk - 1] = (
                        pfv.ric[jj - 1] * pfv.waves[jj - 1, it - 1] * 0.25e6
                    )
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.rpf[jj - 1]
                    pf.zfxf[kk - 1] = pfv.zpf[jj - 1] - dzpf * 0.375e0
                    pf.cfxf[kk - 1] = (
                        pfv.ric[jj - 1] * pfv.waves[jj - 1, it - 1] * 0.25e6
                    )

                else:
                    # Field from different coil
                    kk = kk + 1
                    pf.rfxf[kk - 1] = pfv.rpf[jj - 1]
                    pf.zfxf[kk - 1] = pfv.zpf[jj - 1]
                    pf.cfxf[kk - 1] = (
                        pfv.ric[jj - 1] * pfv.waves[jj - 1, it - 1] * 1.0e6
                    )

        # Plasma contribution
        if it > 2:
            kk = kk + 1
            pf.rfxf[kk - 1] = pv.rmajor
            pf.zfxf[kk - 1] = 0.0e0
            pf.cfxf[kk - 1] = pv.plascur

        # Calculate the field at the inner and outer edges
        # of the coil of interest
        pf.xind[:kk], bri, bzi, psi = bfield(
            pf.rfxf[:kk],
            pf.zfxf[:kk],
            pf.cfxf[:kk],
            pfv.ra[i - 1],
            pfv.zpf[i - 1],
        )
        pf.xind[:kk], bro, bzo, psi = bfield(
            pf.rfxf[:kk],
            pf.zfxf[:kk],
            pf.cfxf[:kk],
            pfv.rb[i - 1],
            pfv.zpf[i - 1],
        )

        # bpf and bpf2 for the Central Solenoid are calculated in OHCALC
        if (bv.iohcl != 0) and (i == pfv.nohc):
            return bri, bro, bzi, bzo

        bpfin = math.sqrt(bri**2 + bzi**2)
        bpfout = math.sqrt(bro**2 + bzo**2)
        for n in range(pfv.ncls[ii - 1]):
            pfv.bpf[i - 1 + n] = bpfin
            pf.bpf2[i - 1 + n] = bpfout

        return bri, bro, bzi, bzo

    def bfmax(self, rj, a, b, h):
        """Calculates the maximum field of a solenoid.

        author: P J Knight, CCFE, Culham Science Centre
        This routine calculates the peak field (T) at a solenoid's
        inner radius, using fits taken from the figure
        on p.22 of M. Wilson's book Superconducting Magnets,
        Clarendon Press, Oxford, N.Y., 1983
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

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
            pf.nef = pfv.ncirt - 1
        else:
            pf.nef = pfv.ncirt - 2

        pfv.vsefsu = 0.0e0

        for i in range(pf.nef):
            pf.vsdum[i, 0] = pfv.sxlg[pfv.ncirt - 1, i] * pfv.cpt[i, 1]
            pf.vsdum[i, 1] = pfv.sxlg[pfv.ncirt - 1, i] * pfv.cpt[i, 2]
            pfv.vsefsu = pfv.vsefsu + (pf.vsdum[i, 1] - pf.vsdum[i, 0])

        # Central Solenoid startup volt-seconds
        if bv.iohcl != 0:
            pf.vsdum[pfv.nohc - 1, 0] = (
                pfv.sxlg[pfv.ncirt - 1, pfv.ncirt - 2] * pfv.cpt[pfv.ncirt - 2, 1]
            )
            pf.vsdum[pfv.nohc - 1, 1] = (
                pfv.sxlg[pfv.ncirt - 1, pfv.ncirt - 2] * pfv.cpt[pfv.ncirt - 2, 2]
            )
            pfv.vsohsu = pf.vsdum[pfv.nohc - 1, 1] - pf.vsdum[pfv.nohc - 1, 0]

        # Total available volt-seconds for start-up
        pfv.vssu = pfv.vsohsu + pfv.vsefsu

        # Burn volt-seconds
        if bv.iohcl != 0:
            pf.vsdum[pfv.nohc - 1, 2] = (
                pfv.sxlg[pfv.ncirt - 1, pfv.ncirt - 2] * pfv.cpt[pfv.ncirt - 2, 4]
            )
            pfv.vsohbn = pf.vsdum[pfv.nohc - 1, 2] - pf.vsdum[pfv.nohc - 1, 1]

        # PF volt-seconds during burn
        pfv.vsefbn = 0.0e0
        for i in range(pf.nef):
            pf.vsdum[i, 2] = pfv.sxlg[pfv.ncirt - 1, i] * pfv.cpt[i, 4]
            pfv.vsefbn = pfv.vsefbn + (pf.vsdum[i, 2] - pf.vsdum[i, 1])

        pfv.vsbn = pfv.vsohbn + pfv.vsefbn

        pfv.vstot = pfv.vssu + pfv.vsbn
        pfv.vseft = pfv.vsefsu + pfv.vsefbn
        pfv.vsoh = pfv.vsohbn + pfv.vsohsu

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
        a = pfv.ra[pfv.nohc - 1]

        # Outer radius of central Solenoid [m]
        b = pfv.rb[pfv.nohc - 1]

        # alpha
        alpha = b / a

        # epsilon
        epsilon = r / a

        # Field at inner radius of coil [T]
        B_a = pfv.bmaxoh0

        # Field at outer radius of coil [T]
        # Assume to be 0 for now
        B_b = 0.0e0

        # current density [A/m^2]
        j = pfv.cohbop

        # K term
        K = ((alpha * B_a - B_b) * j * a) / (alpha - 1.0e0)

        # M term
        M = ((B_a - B_b) * j * a) / (alpha - 1.0e0)

        # calculate hoop stress terms
        hp_term_1 = K * ((2.0e0 + tfv.poisson_steel) / (3.0e0 * (alpha + 1.0e0)))

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

        hp_term_3 = M * ((3.0e0 + tfv.poisson_steel) / (8.0e0))

        hp_term_4 = (
            alpha**2
            + 1.0e0
            + alpha**2 / epsilon**2
            - epsilon**2
            * ((1.0e0 + 3.0e0 * tfv.poisson_steel) / (3.0e0 + tfv.poisson_steel))
        )

        s_hoop_nom = hp_term_1 * hp_term_2 - hp_term_3 * hp_term_4

        s_hoop = s_hoop_nom / pfv.oh_steel_frac

        return s_hoop

    def axial_stress(self):
        """Calculation of axial stress of central solenoid.

        author: J Morris, CCFE, Culham Science Centre
        This routine calculates the axial stress of the central solenoid
        from "Case studies in superconducting magnets", Y. Iwasa, Springer

        :return: unsmeared axial stress [MPa], axial force [N]
        :rtype: tuple[float, float]
        """
        b = pfv.rb[pfv.nohc - 1]

        # Half height of central Solenoid [m]
        hl = pfv.zh[pfv.nohc - 1]

        # Central Solenoid current [A]
        ni = pfv.ric[pfv.nohc - 1] * 1.0e6

        # kb term for elliptical integrals
        # kb2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + hl**2))
        kb2 = (4.0e0 * b**2) / (4.0e0 * b**2 + hl**2)

        # k2b term for elliptical integrals
        # k2b2 = SQRT((4.0e0*b**2)/(4.0e0*b**2 + 4.0e0*hl**2))
        k2b2 = (4.0e0 * b**2) / (4.0e0 * b**2 + 4.0e0 * hl**2)

        # term 1
        axial_term_1 = -(constants.rmu0 / 2.0e0) * (ni / (2.0e0 * hl)) ** 2

        # term 2
        ekb2_1, ekb2_2 = ml.ellipke(kb2)
        axial_term_2 = (
            2.0e0 * hl * (math.sqrt(4.0e0 * b**2 + hl**2)) * (ekb2_1 - ekb2_2)
        )

        # term 3
        ek2b2_1, ek2b2_2 = ml.ellipke(k2b2)
        axial_term_3 = (
            2.0e0
            * hl
            * (math.sqrt(4.0e0 * b**2 + 4.0e0 * hl**2))
            * (ek2b2_1 - ek2b2_2)
        )

        # calculate axial force [N]
        axial_force = axial_term_1 * (axial_term_2 - axial_term_3)

        # axial area [m2]
        area_ax = constants.pi * (pfv.rb[pfv.nohc - 1] ** 2 - pfv.ra[pfv.nohc - 1] ** 2)

        # calculate unsmeared axial stress [MPa]
        s_axial = axial_force / (pfv.oh_steel_frac * 0.5 * area_ax)

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

        pfv.sxlg[:, :] = 0.0

        # Break Central Solenoid into noh segments
        #
        # Choose noh so that the radial thickness of the coil is not thinner
        # than each segment is tall, i.e. the segments are pancake-like,
        # for the benefit of the mutual inductance calculations later

        noh = int(
            math.ceil(
                2.0e0
                * pfv.zh[pfv.nohc - 1]
                / (pfv.rb[pfv.nohc - 1] - pfv.ra[pfv.nohc - 1])
            )
        )

        if noh > nohmax:
            eh.idiags[0] = noh
            eh.idiags[1] = nohmax
            eh.fdiags[0] = bv.ohcth
            eh.report_error(73)

        noh = min(noh, nohmax)

        # TODO In FNSF case, noh = -7! noh should always be positive. Fortran
        # array allocation with -ve bound previously coerced to 0
        if noh < 0:
            noh = 0

        roh = np.zeros(noh)
        zoh = np.zeros(noh)

        if bv.iohcl != 0:
            roh[:] = pfv.rohc

            delzoh = (
                2.0e0 * pfv.zh[pfv.nohc - 1] / noh
            )  # zh(nohc) is the half-height of the coil
            for i in range(noh):
                zoh[i] = pfv.zh[pfv.nohc - 1] - delzoh * (0.5e0 + i)

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
            if bv.ohcth >= delzoh:
                deltar = math.sqrt((bv.ohcth**2 - delzoh**2) / 12.0e0)
            else:
                # eh.fdiags[0] = bv.ohcth
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

            pfv.sxlg[pfv.ncirt - 1, pfv.nohc - 1] = (
                xohpl / (nplas * noh) * pfv.turns[pfv.nohc - 1]
            )
            pfv.sxlg[pfv.nohc - 1, pfv.ncirt - 1] = pfv.sxlg[
                pfv.ncirt - 1, pfv.nohc - 1
            ]

        # Plasma self inductance
        pfv.sxlg[pfv.ncirt - 1, pfv.ncirt - 1] = pv.rlp

        # PF coil / plasma mutual inductances
        ncoils = 0

        for i in range(pfv.ngrp):
            xpfpl = 0.0
            ncoils = ncoils + pfv.ncls[i]
            rp = pfv.rpf[ncoils - 1]
            zp = pfv.zpf[ncoils - 1]
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
            for ii in range(nplas):
                xpfpl = xpfpl + xc[ii]

            for j in range(pfv.ncls[i]):
                ncoilj = ncoils + 1 - (j + 1)
                pfv.sxlg[ncoilj - 1, pfv.ncirt - 1] = (
                    xpfpl / nplas * pfv.turns[ncoilj - 1]
                )
                pfv.sxlg[pfv.ncirt - 1, ncoilj - 1] = pfv.sxlg[
                    ncoilj - 1, pfv.ncirt - 1
                ]

        if bv.iohcl != 0:
            # Central Solenoid self inductance
            a = pfv.rohc  # mean radius of coil
            b = 2.0e0 * pfv.zh[pfv.nohc - 1]  # length of coil
            c = pfv.rb[pfv.nohc - 1] - pfv.ra[pfv.nohc - 1]  # radial winding thickness
            pfv.sxlg[pfv.nohc - 1, pfv.nohc - 1] = self.selfinductance(
                a, b, c, pfv.turns[pfv.nohc - 1]
            )

            # Central Solenoid / PF coil mutual inductances
            for i in range(noh):
                rc[i] = roh[i]
                zc[i] = zoh[i]

            ncoils = 0
            for i in range(pfv.ngrp):
                xohpf = 0.0
                ncoils = ncoils + pfv.ncls[i]
                rp = pfv.rpf[ncoils - 1]
                zp = pfv.zpf[ncoils - 1]
                xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
                for ii in range(noh):
                    xohpf = xohpf + xc[ii]

                for j in range(pfv.ncls[i]):
                    ncoilj = ncoils + 1 - (j + 1)
                    pfv.sxlg[ncoilj - 1, pfv.nohc - 1] = (
                        xohpf * pfv.turns[ncoilj - 1] * pfv.turns[pfv.nohc - 1] / noh
                    )
                    pfv.sxlg[pfv.nohc - 1, ncoilj - 1] = pfv.sxlg[
                        ncoilj - 1, pfv.nohc - 1
                    ]

        # PF coil - PF coil inductances
        if bv.iohcl == 0:
            pf.nef = pfv.nohc
        else:
            pf.nef = pfv.nohc - 1

        for i in range(pf.nef):
            for j in range(pf.nef - 1):
                if j >= i:
                    jj = j + 1 + 1
                else:
                    jj = j + 1

                zc[j] = pfv.zpf[jj - 1]
                rc[j] = pfv.rpf[jj - 1]

            rp = pfv.rpf[i]
            zp = pfv.zpf[i]
            xc, br, bz, psi = bfield(rc, zc, cc, rp, zp)
            for k in range(pf.nef):
                if k < i:
                    pfv.sxlg[i, k] = xc[k] * pfv.turns[k] * pfv.turns[i]
                elif k == i:
                    rl = abs(pfv.zh[k] - pfv.zl[k]) / math.sqrt(constants.pi)
                    pfv.sxlg[k, k] = (
                        constants.rmu0
                        * pfv.turns[k] ** 2
                        * pfv.rpf[k]
                        * (math.log(8.0e0 * pfv.rpf[k] / rl) - 1.75e0)
                    )
                else:
                    pfv.sxlg[i, k] = xc[k - 1] * pfv.turns[k] * pfv.turns[i]

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
                    f"{ig}\t{pfv.sxlg[:pfv.ncirt,ig]}",
                )

            if bv.iohcl != 0:
                op.write(
                    self.outfile,
                    f"CS\t\t\t{pfv.sxlg[:pfv.ncirt,pfv.ncirt-2]}",
                )

            op.write(
                self.outfile,
                f"Plasma\t{pfv.sxlg[:pfv.ncirt,pfv.ncirt-1]}",
            )

    def outpf(self):
        """Routine to write output from PF coil module to file.

        author: P J Knight, CCFE, Culham Science Centre
        This routine writes the PF coil information to the output file.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        op.oheadr(self.outfile, "Central Solenoid and PF Coils")

        if bv.iohcl == 0:
            op.ocmmnt(self.outfile, "No central solenoid included")
            op.oblnkl(self.outfile)
            op.ovarin(
                self.mfile, "Existence_of_central_solenoid", "(bv.iohcl)", bv.iohcl
            )
        else:
            if pfv.ipfres == 0:
                op.ocmmnt(self.outfile, "Superconducting central solenoid")

                op.ovarin(
                    self.outfile,
                    "Central solenoid superconductor material",
                    "(isumatoh)",
                    pfv.isumatoh,
                )

                if pfv.isumatoh == 1:
                    op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
                elif pfv.isumatoh == 2:
                    op.ocmmnt(
                        self.outfile, "  (Bi-2212 high temperature superconductor)"
                    )
                elif pfv.isumatoh == 3:
                    op.ocmmnt(self.outfile, "  (NbTi)")
                elif pfv.isumatoh == 4:
                    op.ocmmnt(
                        self.outfile,
                        "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                    )
                elif pfv.isumatoh == 5:
                    op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
                elif pfv.isumatoh == 6:
                    op.ocmmnt(self.outfile, " (REBCO HTS)")
                elif pfv.isumatoh == 7:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                    )
                elif pfv.isumatoh == 8:
                    op.ocmmnt(
                        self.outfile,
                        " (Durham Ginzburg-Landau critical surface model for REBCO)",
                    )
                elif pfv.isumatoh == 9:
                    op.ocmmnt(
                        self.outfile,
                        " (Hazelton experimental data + Zhai conceptual model for REBCO)",
                    )

                op.osubhd(self.outfile, "Central Solenoid Current Density Limits :")
                op.ovarre(
                    self.outfile,
                    "Maximum field at Beginning Of Pulse (T)",
                    "(bmaxoh0)",
                    pfv.bmaxoh0,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at BOP (A/m2)",
                    "(jscoh_bop)",
                    pfv.jscoh_bop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical strand current density at BOP (A/m2)",
                    "(jstrandoh_bop)",
                    pfv.jstrandoh_bop,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at BOP (A/m2)",
                    "(rjohc0)",
                    pfv.rjohc0,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at BOP (A/m2)",
                    "(cohbop)",
                    pfv.cohbop,
                    "OP ",
                )
                op.oblnkl(self.outfile)
                op.ovarre(
                    self.outfile,
                    "Maximum field at End Of Flattop (T)",
                    "(bmaxoh)",
                    pfv.bmaxoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical superconductor current density at EOF (A/m2)",
                    "(jscoh_eof)",
                    pfv.jscoh_eof,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Critical strand current density at EOF (A/m2)",
                    "(jstrandoh_eof)",
                    pfv.jstrandoh_eof,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Allowable overall current density at EOF (A/m2)",
                    "(rjohc)",
                    pfv.rjohc,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Actual overall current density at EOF (A/m2)",
                    "(coheof)",
                    pfv.coheof,
                )
                op.oblnkl(self.outfile)
                # MDK add bv.ohcth, bv.bore and bv.gapoh as they can be iteration variables
                op.ovarre(self.outfile, "CS inside radius (m)", "(bore)", bv.bore)
                op.ovarre(self.outfile, "CS thickness (m)", "(ohcth)", bv.ohcth)
                op.ovarre(
                    self.outfile,
                    "Gap between central solenoid and TF coil (m)",
                    "(gapoh)",
                    bv.gapoh,
                )
                op.ovarre(
                    self.outfile,
                    "CS overall cross-sectional area (m2)",
                    "(areaoh)",
                    pfv.areaoh,
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
                    "(awpoh*(1-vfohc))",
                    pfv.awpoh * (1.0e0 - pfv.vfohc),
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "   CS void cross-sectional area (m2)",
                    "(awpoh*vfohc)",
                    pfv.awpoh * pfv.vfohc,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel cross-sectional area (m2)",
                    "(areaoh-awpoh)",
                    pfv.areaoh - pfv.awpoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "CS steel area fraction",
                    "(oh_steel_frac)",
                    pfv.oh_steel_frac,
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
                    "(s_tresca_oh)",
                    pfv.s_tresca_oh,
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
                    "(tfcoil_variables.str_cs_con_res)",
                    tfv.str_cs_con_res,
                )
                op.ovarre(
                    self.outfile,
                    "Copper fraction in strand",
                    "(fcuohsu)",
                    pfv.fcuohsu,
                )
                # If REBCO material is used, print copperaoh_m2
                if pfv.isumatoh == 6 or pfv.isumatoh == 8 or pfv.isumatoh == 9:
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
                    "(vfohc)",
                    pfv.vfohc,
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
                    "(tmargoh)",
                    pfv.tmargoh,
                    "OP ",
                )
                op.ovarre(
                    self.outfile,
                    "Minimum permitted temperature margin (K)",
                    "(tmargmin_cs)",
                    tfv.tmargmin_cs,
                )
                # only output CS fatigue model for pulsed reactor design
                if pv.facoh > 0.0e-4:
                    op.ovarre(
                        self.outfile,
                        "Residual hoop stress in CS Steel (Pa)",
                        "(residual_sig_hoop)",
                        csfv.residual_sig_hoop,
                    )
                    op.ovarre(
                        self.outfile,
                        "Minimum burn time (s)",
                        "(tbrnmn)",
                        ctv.tbrnmn,
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
                        "(a_oh_turn)",
                        pfv.a_oh_turn,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn length (m)",
                        "(l_cond_cst)",
                        pfv.l_cond_cst,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn internal cable space radius (m)",
                        "(r_in_cst)",
                        pfv.r_in_cst,
                    )
                    op.ovarre(
                        self.outfile,
                        "CS turn width (m)",
                        "(d_cond_cst)",
                        pfv.d_cond_cst,
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
                    abs(pfv.coheof) > 0.99e0 * abs(numerics.boundu[37] * pfv.rjohc)
                ) or (abs(pfv.cohbop) > 0.99e0 * abs(numerics.boundu[38] * pfv.rjohc0)):
                    pf.cslimit = True
                if pfv.tmargoh < 1.01e0 * tfv.tmargmin_cs:
                    pf.cslimit = True
                if not pf.cslimit:
                    eh.report_error(135)

                # REBCO fractures in strains above ~+/- 0.7%
                if (
                    (pfv.isumatoh == 6 or pfv.isumatoh == 8 or pfv.isumatoh == 9)
                    and tfv.str_cs_con_res > 0.7e-2
                    or tfv.str_cs_con_res < -0.7e-2
                ):
                    eh.report_error(262)

                if (
                    (pfv.isumatpf == 6 or pfv.isumatpf == 8 or pfv.isumatpf == 9)
                    and tfv.str_pf_con_res > 0.7e-2
                    or tfv.str_pf_con_res < -0.7e-2
                ):
                    eh.report_error(263)

            else:
                op.ocmmnt(self.outfile, "Resistive central solenoid")

        if pfv.ipfres == 0:
            op.oblnkl(self.outfile)
            op.ocmmnt(self.outfile, "Superconducting PF coils")

            op.ovarin(
                self.outfile,
                "PF coil superconductor material",
                "(isumatpf)",
                pfv.isumatpf,
            )

            if pfv.isumatpf == 1:
                op.ocmmnt(self.outfile, "  (ITER Nb3Sn critical surface model)")
            elif pfv.isumatpf == 2:
                op.ocmmnt(self.outfile, "  (Bi-2212 high temperature superconductor)")
            elif pfv.isumatpf == 3:
                op.ocmmnt(self.outfile, "  (NbTi)")
            elif pfv.isumatpf == 4:
                op.ocmmnt(
                    self.outfile,
                    "  (ITER Nb3Sn critical surface model, user-defined parameters)",
                )
            elif pfv.isumatpf == 5:
                op.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
            elif pfv.isumatpf == 6:
                op.ocmmnt(
                    self.outfile,
                    " (REBCO 2nd generation HTS superconductor in CrCo strand)",
                )
            elif pfv.isumatpf == 7:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for Nb-Ti)",
                )
            elif pfv.isumatpf == 8:
                op.ocmmnt(
                    self.outfile,
                    " (Durham Ginzburg-Landau critical surface model for REBCO)",
                )
            elif pfv.isumatpf == 9:
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
                "(powpfres)",
                pfv.powpfres,
                "OP ",
            )
            if bv.iohcl != 0:
                op.ovarre(
                    self.outfile,
                    "Central solenoid resistive power (W)",
                    "(powohres)",
                    pfv.powohres,
                    "OP ",
                )

        # pf.nef is the number of coils excluding the Central Solenoid
        pf.nef = pfv.nohc
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
                f"PF {k}\t\t\t{pfv.rpf[k]:.2e}\t{pfv.zpf[k]:.2e}\t{pfv.rb[k]-pfv.ra[k]:.2e}\t{abs(pfv.zh[k]-pfv.zl[k]):.2e}\t{pfv.turns[k]:.2e}",
            )

        for k in range(pf.nef):
            op.ovarre(
                self.mfile,
                f"PF coil {k} radius (m)",
                f"(rpf[{k}]",
                pfv.rpf[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical position (m)",
                f"(zpf[{k}])",
                pfv.zpf[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} radial thickness (m)",
                f"(pfdr({k}))",
                pfv.rb[k] - pfv.ra[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} vertical thickness (m)",
                f"(pfdz({k}))",
                pfv.zh[k] - pfv.zl[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} turns",
                f"(turns[{k}])",
                pfv.turns[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} current (MA)",
                f"(ric[{k}])",
                pfv.ric[k],
            )
            op.ovarre(
                self.mfile,
                f"PF coil {k} field (T)",
                f"(bpf[{k}])",
                pfv.bpf[k],
            )
        self.tf_pf_collision_detector()

        # Central Solenoid, if present
        if bv.iohcl != 0:
            op.write(
                self.outfile,
                f"CS\t\t\t\t{pfv.rpf[pfv.nohc-1]:.2e}\t{pfv.zpf[pfv.nohc-1]:.2e}\t{pfv.rb[pfv.nohc-1]-pfv.ra[pfv.nohc-1]:.2e}\t{abs(pfv.zh[pfv.nohc-1]-pfv.zl[pfv.nohc-1]):.2e}\t{pfv.turns[pfv.nohc-1]:.2e}\t{pfv.pfcaseth[pfv.nohc-1]:.2e}",
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radius (m)",
                "(rpf[nohc-1])",
                pfv.rpf[pfv.nohc - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical position (m)",
                "(zpf[nohc-1])",
                pfv.zpf[pfv.nohc - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid radial thickness (m)",
                "(ohdr)",
                (pfv.rb[pfv.nohc - 1] - pfv.ra[pfv.nohc - 1]),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid vertical thickness (m)",
                "(ohdz)",
                (pfv.zh[pfv.nohc - 1] - pfv.zl[pfv.nohc - 1]),
            )
            op.ovarre(
                self.mfile,
                "Central solenoid turns",
                "(turns[nohc-1])",
                pfv.turns[pfv.nohc - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid current (MA)",
                "(ric[nohc-1])",
                pfv.ric[pfv.nohc - 1],
            )
            op.ovarre(
                self.mfile,
                "Central solenoid field (T)",
                "(bpf[nohc-1])",
                pfv.bpf[pfv.nohc - 1],
            )

        # Plasma
        op.write(
            self.outfile,
            f"Plasma\t\t\t{pv.rmajor:.2e}\t0.0e0\t\t{2.0e0*pv.rminor:.2e}\t{2.0e0*pv.rminor*pv.kappa:.2e}\t1.0e0",
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
            if pfv.ipfres == 0:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfv.ric[k]:.2e}\t{pfv.rjpfalw[k]:.2e}\t{pfv.rjconpf[k]:.2e}\t{pfv.rjconpf[k]/pfv.rjpfalw[k]:.2e}\t{pfv.wtc[k]:.2e}\t{pfv.wts[k]:.2e}\t{pfv.bpf[k]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"PF {k}\t{pfv.ric[k]:.2e}\t-1.0e0\t{pfv.rjconpf[k]:.2e}\t1.0e0\t{pfv.wtc[k]:.2e}\t{pfv.wts[k]:.2e}\t{pfv.bpf[k]:.2e}\t",
                )

        # Central Solenoid, if present
        if bv.iohcl != 0:
            if pfv.ipfres == 0:
                # Issue #328
                op.write(
                    self.outfile,
                    f"CS\t\t{pfv.ric[pfv.nohc-1]:.2e}\t{pfv.rjpfalw[pfv.nohc-1]:.2e}\t{max(abs(pfv.cohbop),abs(pfv.coheof)):.2e}\t{max(abs(pfv.cohbop),abs(pfv.coheof))/pfv.rjpfalw[pfv.nohc-1]:.2e}\t{pfv.wtc[pfv.nohc-1]:.2e}\t{pfv.wts[pfv.nohc-1]:.2e}\t{pfv.bpf[pfv.nohc-1]:.2e}",
                )
            else:
                op.write(
                    self.outfile,
                    f"CS\t\t{pfv.ric[pfv.nohc-1]:.2e}\t-1.0e0\t{max(abs(pfv.cohbop)):.2e}\t{abs(pfv.coheof):.2e}\t1.0e0\t{pfv.wtc[pfv.nohc-1]:.2e}\t{pfv.wts[pfv.nohc-1]:.2e}\t{pfv.bpf[pfv.nohc-1]:.2e}",
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
            + f"{pfv.whtpf:.2e}\t{pfv.whtpfs:.2e}",
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
        op.write(self.outfile, "\t" * 3 + "start-up\t\t\tburn\t\t\ttotal")
        op.write(
            self.outfile,
            f"PF coils:\t\t{pfv.vsefsu:.2f}\t\t\t\t{pfv.vsefbn:.2f}\t\t\t{pfv.vseft:.2f}",
        )
        op.write(
            self.outfile,
            f"CS coil:\t\t{pfv.vsohsu:.2f}\t\t\t\t{pfv.vsohbn:.2f}\t\t\t{pfv.vsoh:.2f}",
        )
        op.write(
            self.outfile, "\t" * 3 + "-" * 7 + "\t" * 4 + "-" * 7 + "\t" * 3 + "-" * 7
        )
        op.write(
            self.outfile,
            f"Total:\t\t\t{pfv.vssu:.2f}\t\t\t\t{pfv.vsbn:.2f}\t\t\t{pfv.vstot:.2f}",
        )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Total volt-second consumption by coils (Wb)",
            "(vstot)",
            f"{pfv.vstot:.2}",
            "OP",
        )

        op.osubhd(self.outfile, "Summary of volt-second consumption by circuit (Wb):")

        op.write(self.outfile, "circuit\t\t\tBOP\t\t\tBOF\t\tEOF")
        op.oblnkl(self.outfile)

        for k in range(pf.nef):
            op.write(
                self.outfile,
                f"\t{k}\t\t\t{pf.vsdum[k,0]:.3f}\t\t\t{pf.vsdum[k,1]:.3f}\t\t{pf.vsdum[k,2]:.3f}",
            )

        op.write(
            self.outfile,
            f"\tCS coil\t\t\t{pf.vsdum[pfv.nohc-1,0]:.3f}\t\t\t{pf.vsdum[pfv.nohc-1,1]:.3f}\t\t{pf.vsdum[pfv.nohc-1,2]:.3f}",
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

        for k in range(pfv.ncirt - 1):
            line = f"\t{k}\t\t"
            for jj in range(6):
                line += f"\t{pfv.cpt[k,jj]*pfv.turns[k]:.3e}"
            op.write(self.outfile, line)

        line = "Plasma (A)\t\t"
        for jj in range(6):
            line += f"\t{pfv.cpt[pfv.ncirt-1,jj]:.3e}"

        op.write(self.outfile, line)

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "This consists of: CS coil field balancing:")
        for k in range(pfv.ncirt - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{pfv.cpt[k,0]*pfv.turns[k]:.3e}\t"
                    f"{pfv.cpt[k,1]*pfv.turns[k]:.3e}\t"
                    f"{-pfv.cpt[k,1]*pfv.turns[k]*(pfv.fcohbof/pfv.fcohbop):.3e}\t"
                    f"{-pfv.cpt[k,1]*pfv.turns[k]*(pfv.fcohbof/pfv.fcohbop):.3e}\t"
                    f"{-pfv.cpt[k,1]*pfv.turns[k]*(1.0e0/pfv.fcohbop):.3e}\t"
                    f"{pfv.cpt[k,5]*pfv.turns[k]:.3e}"
                ),
            )

        op.oblnkl(self.outfile)
        op.ocmmnt(self.outfile, "And: equilibrium field:")
        for k in range(pfv.ncirt - 1):
            op.write(
                self.outfile,
                (
                    f"{k}\t\t\t{0.0:.3e}\t{0.0:.3e}\t"
                    f"{(pfv.cpt[k,2]+pfv.cpt[k,1]*pfv.fcohbof/pfv.fcohbop)*pfv.turns[k]:.3e}\t"
                    f"{(pfv.cpt[k,3]+pfv.cpt[k,1]*pfv.fcohbof/pfv.fcohbop)*pfv.turns[k]:.3e}\t"
                    f"{(pfv.cpt[k,4]+pfv.cpt[k,1]*1.0e0/pfv.fcohbop)*pfv.turns[k]:.3e}\t"
                    "0.0e0"
                ),
            )

        op.oblnkl(self.outfile)
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Pulse / end of flat-top",
            "(fcohbop)",
            pfv.fcohbop,
        )
        op.ovarre(
            self.outfile,
            "Ratio of central solenoid current at beginning of Flat-top / end of flat-top",
            "(fcohbof)",
            pfv.fcohbof,
            "OP ",
        )

        op.oshead(self.outfile, "PF Circuit Waveform Data")
        op.ovarin(
            self.outfile,
            "Number of PF circuits including CS and plasma",
            "(ncirt)",
            pfv.ncirt,
        )
        for k in range(pfv.ncirt):
            for jjj in range(6):
                if k == pfv.ncirt - 1:
                    circuit_name = f"Plasma Time point {jjj} (A)"
                    circuit_var_name = f"(plasmat{jjj})"
                elif k == pfv.ncirt - 2:
                    circuit_name = f"CS Circuit Time point {jjj} (A)"
                    circuit_var_name = f"(cs t{jjj})"
                else:
                    circuit_name = f"PF Circuit {k} Time point {jjj} (A)"
                    circuit_var_name = f"(pfc{k}t{jjj})"

                op.ovarre(
                    self.outfile,
                    circuit_name,
                    circuit_var_name,
                    pfv.cpt[k, jjj] * pfv.turns[k],
                )

    def selfinductance(self, a, b, c, N):
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
        selfinductance = (
            (1.0e-6 / 0.0254e0)
            * a**2
            * N**2
            / (9.0e0 * a + 10.0e0 * b + 8.4e0 * c + 3.2e0 * c * b / a)
        )
        return selfinductance

    def waveform(self):
        """Sets up the PF coil waveforms.

        author: P J Knight, CCFE, Culham Science Centre
        This routine sets up the PF coil current waveforms.
        waves[i,j] is the current in coil i, at time j,
        normalized to the peak current in that coil at any time.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        nplas = pfv.nohc + 1
        for it in range(6):
            pfv.waves[nplas - 1, it] = 1.0e0

        for ic in range(pfv.nohc):
            # Find where the peak current occurs
            # Beginning of pulse, t = tramp
            if (abs(pfv.curpfs[ic]) >= abs(pfv.curpfb[ic])) and (
                abs(pfv.curpfs[ic]) >= abs(pfv.curpff[ic])
            ):
                pfv.ric[ic] = pfv.curpfs[ic]

            # Beginning of flat-top, t = tramp + tohs
            if (abs(pfv.curpff[ic]) >= abs(pfv.curpfb[ic])) and (
                abs(pfv.curpff[ic]) >= abs(pfv.curpfs[ic])
            ):
                pfv.ric[ic] = pfv.curpff[ic]

            # End of flat-top, t = tramp + tohs + theat + tburn
            if (abs(pfv.curpfb[ic]) >= abs(pfv.curpfs[ic])) and (
                abs(pfv.curpfb[ic]) >= abs(pfv.curpff[ic])
            ):
                pfv.ric[ic] = pfv.curpfb[ic]

            # Set normalized current waveforms
            pfv.waves[ic, 0] = 0.0e0
            pfv.waves[ic, 1] = pfv.curpfs[ic] / pfv.ric[ic]
            pfv.waves[ic, 2] = pfv.curpff[ic] / pfv.ric[ic]
            pfv.waves[ic, 3] = pfv.curpff[ic] / pfv.ric[ic]
            pfv.waves[ic, 4] = pfv.curpfb[ic] / pfv.ric[ic]
            pfv.waves[ic, 5] = 0.0e0

    def superconpf(
        self, bmax, fhe, fcu, jwp, isumat, fhts, strain, thelium, bcritsc, tcritsc
    ):
        """Routine to calculate the PF coil superconductor properties.

        This routine calculates the superconductor critical winding pack
        current density for the PF coils, plus the temperature margin.
        It is based on the TF coil version, supercon.
        author: P J Knight, CCFE, Culham Science Centre

        :param bmax: peak field at conductor (T)
        :type bmax: float
        :param fhe: fraction of cable space that is for He cooling
        :type fhe: float
        :param fcu: fraction of strand that is copper
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
        :return: Critical winding pack current density (A/m2) (jcritwp),
        Critical strand current density (A/m2) (jcritstr)
        Critical superconductor current density (A/m2) (jcritsc)
        Temperature margin (K) (tmarg)
        :rtype: tuple[float, float, float, float]
        """

        """TODO maths_library.secant_solve() requires a function of one variable,
        e.g. f(x). However, this function can require other variables as arguments
        e.g. constants. Access to these variables (e.g. bmax, bc20m, tc0m) has
        previously been provided through nested functions with implicit access
        to the parent scope of superconpf() when the functions are passed to
        secant_solve(). Now in Python, it might be better to explcitly pass
        these variables as optional argument(s) to secant_solve() and remove
        this nested function requirement.
        """

        def deltaj_nbti(temperature):
            """Critical current density and current density difference in NbTi.

            :param temperature: temperature
            :type temperature: float
            :return: difference in current density
            :rtype: float
            """
            jcrit0, _ = superconductors.jcrit_nbti(temperature, bmax, c0, bc20m, tc0m)
            if ml.variable_error(jcrit0):  # superconductors.jcrit_nbti has failed.
                print(f"superconductors.jcrit_nbti: {bmax=} {temperature=} {jcrit0=}")

            deltaj_nbti = jcrit0 - jsc
            return deltaj_nbti

        def deltaj_wst(temperature):
            """Critical current density and current density difference for WST Nb3Sn.

            :param temperature: temperature
            :type temperature: float
            :return: difference in current density
            :rtype: float
            """
            jcrit0, _, _ = superconductors.wstsc(temperature, bmax, strain, bc20m, tc0m)
            if ml.variable_error(jcrit0):  # superconductors.wstsc has failed.
                print(f"deltaj_wst: {bmax=} {temperature=} {jcrit0=}")

            deltaj_wst = jcrit0 - jsc
            return deltaj_wst

        def deltaj_gl_nbti(temperature):
            """Critical current density and current density difference in GL NbTi.

            :param temperature: temperature
            :type temperature: float
            :return: difference in current density
            :rtype: float
            """
            jcrit0, _, _ = superconductors.gl_nbti(
                temperature, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(jcrit0):  # GL_Nbti has failed.
                print(f"deltaj_GL_nbti: {bmax=} {temperature=} {jcrit0=}")

            deltaj_gl_nbti = jcrit0 - jsc
            return deltaj_gl_nbti

        def deltaj_gl_rebco(temperature):
            """Critical current density and current density difference in GL REBCO.

            :param temperature: temperature
            :type temperature: float
            :return: difference in current density
            :rtype: float
            """
            jcrit0, _, _ = superconductors.gl_rebco(
                temperature, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(jcrit0):  # superconductors.GL_REBCO has failed.
                print(f"deltaj_gl_REBCO: {bmax=} {temperature=} {jcrit0=}")

            deltaj_gl_rebco = jcrit0 - jsc
            return deltaj_gl_rebco

        def deltaj_hijc_rebco(temperature):
            """Critical current density and current density difference in high current density REBCO.

            :param temperature: temperature
            :type temperature: float
            :return: difference in current density
            :rtype: float
            """
            jcrit0, _, _ = superconductors.hijc_rebco(
                temperature, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(jcrit0):  # superconductors.GL_REBCO has failed.
                print(f"deltaj_hijc_REBCO: {bmax=} {temperature=} {jcrit0=}")

            deltaj_hijc_rebco = jcrit0 - jsc
            return deltaj_hijc_rebco

        # Find critical current density in superconducting strand, jcritstr
        if isumat == 1:
            # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97e0
            tc0m = 16.06e0

            # jcritsc returned by superconductors.itersc is the critical current density in the
            # superconductor - not the whole strand, which contains copper

            jcritsc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)

        elif isumat == 2:
            # Bi-2212 high temperature superconductor parameterization

            # Current density in a strand of Bi-2212 conductor
            # N.B. jcrit returned by superconductors.bi2212 is the critical current density
            # in the strand, not just the superconducting portion.
            # The parameterization for jcritstr assumes a particular strand
            # composition that does not require a user-defined copper fraction,
            # so this is irrelevant in this model

            jstrand = jwp / (1.0e0 - fhe)

            jcritstr, tmarg = superconductors.bi2212(bmax, jstrand, thelium, fhts)
            jcritsc = jcritstr / (1.0e0 - fcu)

        elif isumat == 3:
            # NbTi data
            bc20m = 15.0e0
            tc0m = 9.3e0
            c0 = 1.0e10
            jcritsc, _ = superconductors.jcrit_nbti(thelium, bmax, c0, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)

        elif isumat == 4:
            # As (1), but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            jcritsc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)

        elif isumat == 5:
            # WST Nb3Sn parameterisation
            bc20m = 32.97e0
            tc0m = 16.06e0

            # jcritsc returned by superconductors.itersc is the critical current density in the
            # superconductor - not the whole strand, which contains copper

            jcritsc, _, _ = superconductors.wstsc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)

        elif isumat == 6:
            # "REBCO" 2nd generation HTS superconductor in CrCo strand
            jcritsc, _ = superconductors.jcrit_rebco(thelium, bmax)
            jcritstr = jcritsc * (1.0e0 - fcu)

            # The CS coil current at EOF
            ioheof = bv.hmax * pfv.ohhghf * bv.ohcth * 2.0 * pfv.coheof
            # The CS coil current/copper area calculation for quench protection
            # Copper area = (area of coil - area of steel)*(1- void fraction)*
            # (fraction of copper in strands)
            rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.vfohc) * pfv.fcuohsu)

        elif isumat == 7:
            # Durham Ginzburg-Landau critical surface model for Nb-Ti
            bc20m = tfv.b_crit_upper_nbti
            tc0m = tfv.t_crit_nbti
            jcritsc, _, _ = superconductors.gl_nbti(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)

            # The CS coil current at EOF
            ioheof = bv.hmax * pfv.ohhghf * bv.ohcth * 2.0 * pfv.coheof

        elif isumat == 8:
            # Durham Ginzburg-Landau critical surface model for REBCO
            bc20m = 429e0
            tc0m = 185e0
            jcritsc, _, _ = superconductors.gl_rebco(thelium, bmax, strain, bc20m, tc0m)
            # A0 calculated for tape cross section already
            jcritstr = jcritsc * (1.0e0 - fcu)

            # The CS coil current at EOF
            ioheof = bv.hmax * pfv.ohhghf * bv.ohcth * 2.0 * pfv.coheof
            # The CS coil current/copper area calculation for quench protection
            rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.vfohc) * pfv.fcuohsu)

        elif isumat == 9:
            # Hazelton experimental data + Zhai conceptual model for REBCO
            bc20m = 138
            tc0m = 92
            jcritsc, _, _ = superconductors.hijc_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            jcritstr = jcritsc * (1.0e0 - fcu)

            # The CS coil current at EOF
            ioheof = bv.hmax * pfv.ohhghf * bv.ohcth * 2.0 * pfv.coheof
            # The CS coil current/copper area calculation for quench protection
            rcv.copperaoh_m2 = ioheof / (pfv.awpoh * (1.0 - pfv.vfohc) * pfv.fcuohsu)

        else:
            # Error condition
            eh.idiag[0] = isumat
            eh.report_error(156)

        # Critical current density in winding pack
        jcritwp = jcritstr * (1.0e0 - fhe)
        jstrand = jwp / (1.0e0 - fhe)
        jsc = jstrand / (1.0e0 - fcu)

        # Temperature margin (already calculated in superconductors.bi2212 for isumat=2)
        if (isumat == 1) or (isumat == 4):
            # Newton-Raphson method; start at requested minimum temperature margin
            ttest = thelium + tfv.tmargmin_cs
            delt = 0.01e0
            jtol = 1.0e4

            # Actual current density in superconductor, which should be equal to jcrit(thelium+tmarg)
            # when we have found the desired value of tmarg
            for lap in range(100):
                if ttest <= 0.0:
                    eh.idiags[0] = lap
                    eh.fdiags[0] = ttest
                    eh.report_error(158)
                    break

                ttestm = ttest - delt
                ttestp = ttest + delt

                if isumat in [1, 4]:
                    jcrit0, _, _ = superconductors.itersc(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break

                    jcritm, _, _ = superconductors.itersc(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.itersc(
                        ttestp, bmax, strain, bc20m, tc0m
                    )

                # Kludge to avoid divide by 0
                if jcritm == jcritp:
                    jcritp = jcritp + (jcritp * 1e-6)

                ttest = ttest - 2.0e0 * delt * (jcrit0 - jsc) / (jcritp - jcritm)
            else:
                eh.idiags[0] = lap
                eh.fdiags[0] = ttest
                eh.report_error(158)

            tmarg = ttest - thelium

        # MDK 13/7/18 Use secant solver for NbTi.
        elif isumat == 3:
            x1 = 4e0  # Initial values of temperature
            x2 = 6e0
            # Solve for deltaj_nbti = 0
            current_sharing_t, error, residual = pml.secant_solve(
                deltaj_nbti, x1, x2, 100e0
            )
            tmarg = current_sharing_t - thelium
            jcrit0, _ = superconductors.jcrit_nbti(
                current_sharing_t, bmax, c0, bc20m, tc0m
            )
            if ml.variable_error(
                current_sharing_t
            ):  # current sharing secant solver has failed.
                print(
                    f"NbTi: {current_sharing_t=} {tmarg=} {jsc=} {jcrit0=} {residual=}"
                )

        # MDK 13/7/18 Use secant solver for WST.
        elif isumat == 5:
            # Current sharing temperature for WST Nb3Sn
            x1 = 4e0  # Initial values of temperature
            x2 = 6e0
            # Solve for deltaj_wst = 0
            current_sharing_t, error, residual = pml.secant_solve(
                deltaj_wst, x1, x2, 100e0
            )
            tmarg = current_sharing_t - thelium
            jcrit0, _, _ = superconductors.wstsc(
                current_sharing_t, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(
                current_sharing_t
            ):  # current sharing secant solver has failed.
                print(
                    f"WST: {current_sharing_t=} {tmarg=} {jsc=} {jcrit0=} {residual=}"
                )

        # Temperature margin: An alternative method using secant solver
        elif isumat == 6:
            current_sharing_t = superconductors.current_sharing_rebco(bmax, jsc)
            tmarg = current_sharing_t - thelium
            tfv.temp_margin = tmarg

        # SCM 16/03/20 Use secant solver for GL_nbti.
        elif isumat == 7:
            # Current sharing temperature for Durham Ginzburg-Landau Nb-Ti
            x1 = 4.0e0  # Initial values of temperature
            x2 = 6.0e0
            # Solve for deltaj_GL_nbti = 0
            current_sharing_t, error, residual = pml.secant_solve(
                deltaj_gl_nbti, x1, x2, 100e0
            )
            tmarg = current_sharing_t - thelium
            jcrit0, _, _ = superconductors.gl_nbti(
                current_sharing_t, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(
                current_sharing_t
            ):  # current sharing secant solver has failed.
                print(
                    f"Gl_nbti: {current_sharing_t=} {tmarg=} {jsc=} {jcrit0=} {residual=}"
                )

        # SCM 10/08/20 Use secant solver for superconductors.GL_REBCO.
        elif isumat == 8:
            # Current sharing temperature for Durham Ginzburg-Landau REBCO
            x1 = 4.0e0  # Initial values of temperature
            x2 = 6.0e0
            # Solve for deltaj_superconductors.GL_REBCO = 0
            current_sharing_t, error, residual = pml.secant_solve(
                deltaj_gl_rebco, x1, x2, 100e0
            )
            tmarg = current_sharing_t - thelium
            jcrit0, _, _ = superconductors.gl_rebco(
                current_sharing_t, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(
                current_sharing_t
            ):  # current sharing secant solver has failed.
                print(
                    f"Gl_REBCO: {current_sharing_t=} {tmarg=} {jsc=} {jcrit0=} {residual=}"
                )

        elif isumat == 9:
            # Current sharing temperature for Hazelton REBCO
            x1 = 19.0e0  # Initial values of temperature
            x2 = 21.0e0
            # Solve for deltaj_superconductors.HIJC_REBCO = 0
            current_sharing_t, error, residual = pml.secant_solve(
                deltaj_hijc_rebco, x1, x2, 100e0
            )
            tmarg = current_sharing_t - thelium
            jcrit0, _, _ = superconductors.hijc_rebco(
                current_sharing_t, bmax, strain, bc20m, tc0m
            )
            if ml.variable_error(
                current_sharing_t
            ):  # current sharing secant solver has failed.
                print(
                    f"HIJC_REBCO: {current_sharing_t=} {tmarg=} {jsc=} {jcrit0=} {residual=}"
                )

        return jcritwp, jcritstr, jcritsc, tmarg


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
def rsid(npts, brin, bzin, nfix, ngrp, ccls, bfix, gmat):
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
    :param ngrp: number of coil groups, where all coils in a group have the
    same current, <= ngrpmx
    :type ngrp: int
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

        for j in range(ngrp):
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
        for j in range(ngrp):
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
    ngrp,
    ncls,
    rcls,
    zcls,
    alfa,
    bfix,
    nclsmx,
):
    """Calculate the currents in a group of ring coils.

    Set up the matrix equation to calculate the currents in a group of ring
    coils.
    author: P J Knight, CCFE, Culham Science Centre
    author: D Strickler, ORNL
    author: J Galambos, ORNL

    :param lrow1: row length of arrays bfix, bvec, gmat, umat, vmat; should
    be >= (2*nptsmx + ngrpmx)
    :type lrow1: int
    :param lcol1: column length of arrays gmat, umat, vmat; should be >=
    ngrpmx
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
    :param ngrp: number of coil groups, where all coils in a group have the
    same current, <= ngrpmx
    :type ngrp: int
    :param ncls: number of coils in each group, each value <= nclsmx
    :type ncls: numpy.ndarray
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
    cc = np.ones(nclsmx)

    for i in range(npts):
        bvec[i] = brin[i] - bfix[i]
        bvec[i + npts] = bzin[i] - bfix[i + npts]

        for j in range(ngrp):
            nc = ncls[j]

            _, gmat[i, j], gmat[i + npts, j], _ = bfield(
                rcls[j, :nc], zcls[j, :nc], cc[:nc], rpts[i], zpts[i]
            )

    # Add constraint equations
    nrws = 2 * npts

    bvec[nrws : nrws + ngrp] = 0.0
    np.fill_diagonal(gmat[nrws : nrws + ngrp, :ngrp], ncls[:ngrp] * alfa)

    nrws = 2 * npts + ngrp

    # numba doesnt like np.zeros(..., order="F") so this acts as a work
    # around to that missing signature
    gmat = np.asfortranarray(gmat)

    return nrws, gmat, bvec
