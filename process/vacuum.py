import logging
import math
import numpy as np

from process.utilities.f2py_string_patch import f2py_compatible_to_string

from process.fortran import constants
from process.fortran import physics_variables as pv
from process.fortran import vacuum_variables as vacv
from process.fortran import build_variables as buv
from process.fortran import tfcoil_variables as tfv
from process.fortran import times_variables as tv
from process.fortran import process_output as po
from process.fortran import error_handling as eh

logger = logging.getLogger(__name__)


class Vacuum:
    """Module containing vacuum system routines
    author: P J Knight, CCFE, Culham Science Centre

    This module contains routines for calculating the
    parameters of the vacuum system for a fusion power plant.
    AEA FUS 251: A User's Guide to the PROCESS Systems Code
    """

    def __init__(self) -> None:
        self.outfile: int = constants.nout

    def run(self, output: bool) -> None:
        """Routine to call the vacuum module
        author: P J Knight, CCFE, Culham Science Centre

        This routine calls the main vacuum package.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        # (should be) NBI gas load (deuterons/second)

        qtorus = 0.0e0

        #  Total fuel gas load (kg/s)
        #  2 nuclei * nucleus-pairs/sec * mass/nucleus

        # MDK Check this!!
        gasld = 2.0e0 * pv.qfuel * pv.afuel * constants.umass

        self.vacuum_model = f2py_compatible_to_string(vacv.vacuum_model)

        # vacuum_model required to be compared to a b string
        # as this is what f2py returns
        if self.vacuum_model == "old":
            pumpn, vacv.nvduct, vacv.dlscal, vacv.vacdshm, vacv.vcdimax = self.vacuum(
                pv.powfmw,
                pv.rmajor,
                pv.rminor,
                0.5e0 * (buv.scrapli + buv.scraplo),
                pv.sarea,
                pv.vol,
                buv.shldoth,
                buv.shldith,
                buv.tfcth,
                buv.rsldi - buv.gapds - buv.d_vv_in,
                tfv.n_tf,
                tv.tdwell,
                pv.dene,
                pv.idivrt,
                qtorus,
                gasld,
                output=output,
            )
            # MDK pumpn is real: convert to integer by rounding.
            vacv.vpumpn = math.floor(pumpn + 0.5e0)
        elif self.vacuum_model == "simple":
            vacv.niterpump = self.vacuum_simple(output=output)
        else:
            logger.warning(f"vacuum_model seems to be invalid: {vacv.vacuum_model}")
            po.ocmmnt(
                self.outfile,
                f'ERROR "vacuum_model" seems to be invalid: {vacv.vacuum_model}',
            )

    def vacuum_simple(self, output) -> float:
        """Simple model of vacuum pumping system
        author: MD Kovari, CCFE, Culham Science Centre

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean

        :return npump: number of pumps for pumpdown and steady-state
        :rtype: float
        """

        # Steady-state model (super simple)
        # One ITER torus cryopump has a throughput of 50 Pa m3/s = 1.2155e+22 molecules/s
        # Issue #304
        niterpump = pv.qfuel / vacv.pumptp

        # Pump-down:
        # Pumping speed per pump m3/s
        pumpspeed = (
            vacv.pumpspeedmax
            * vacv.pumpareafraction
            * vacv.pumpspeedfactor
            * pv.sarea
            / tfv.n_tf
        )

        wallarea = (pv.sarea / 1084.0e0) * 2000.0e0
        # Required pumping speed for pump-down
        pumpdownspeed = (vacv.outgasfactor * wallarea / vacv.pbase) * tv.tdwell ** (
            -vacv.outgasindex
        )
        # Number of pumps required for pump-down
        npumpdown = pumpdownspeed / pumpspeed

        # Combine the two (somewhat inconsistent) models
        # Note that 'npump' can be constrained by constraint equation 63
        npump = max(niterpump, npumpdown)

        #  Output section
        if output:

            po.oheadr(self.outfile, "Vacuum System")
            po.ovarst(
                self.outfile,
                "Switch for vacuum pumping model",
                "(vacuum_model)",
                '"' + self.vacuum_model + '"',
            )
            po.ocmmnt(
                self.outfile,
                "Simple steady-state model with comparison to ITER cryopumps",
            )
            po.ovarre(
                self.outfile,
                "Plasma fuelling rate (nucleus-pairs/s)",
                "(qfuel)",
                pv.qfuel,
                "OP ",
            )
            po.ocmmnt(
                self.outfile, "Number of high vacuum pumps, each with the throughput"
            )
            po.ocmmnt(
                self.outfile,
                " of one ITER cryopump (50 Pa m3 s-1 = 1.2e+22 molecules/s),",
            )
            po.ovarre(
                self.outfile,
                " all operating at the same time",
                "(niterpump)",
                niterpump,
                "OP ",
            )

            po.ovarre(self.outfile, "Dwell time", "(tdwell)", tv.tdwell)
            po.ovarre(
                self.outfile,
                "Number of pumps required for pump-down",
                "(npumpdown)",
                npumpdown,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Number of pumps required overall",
                "(npump)",
                npump,
                "OP ",
            )

        return npump

    def vacuum(
        self,
        pfusmw,
        r0,
        aw,
        dsol,
        plasma_sarea,
        plasma_vol,
        thshldo,
        thshldi,
        thtf,
        ritf,
        n_tf,
        tdwell,
        nplasma,
        ndiv,
        qtorus,
        gasld,
        output,
    ):
        """Routine to calculate the parameters of the vacuum system
        author: P J Knight, CCFE, Culham Science Centre
        author: J Haines, FEDC (originator)
        author: P Dubois, LLNL
        author: J Galambos, ORNL
        author: P C Shipe, ORNL

        This routine calculates the parameters of the vacuum system.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param pfusmw: Fusion power (MW)
        :type pfusmw: float

        :param r0: Major radius (m)
        :type r0: float

        :param aw: Minor radius (m)
        :type aw: float

        :param dsol: Scrape-off layer average width (m)
        :type : float

        :param plasma_sarea: Plasma surface area (m2)
        :type : float

        :param plasma_vol: Plasma volume (m3)
        :type : float

        :param thshldo: Outboard shield thickness (m)
        :type : float

        :param thshldi: Inboard shield thickness (m)
        :type : float

        :param thtf:  TF coil thickness (m)
        :type : float

        :param ritf: Radius of inboard TF leg point nearest plasma (m)
        :type : float

        :param tfno:  Number of TF coils
        :type : int

        :param tdwell: Dwell time between pulses (s)
        :type : float

        :param nplasma: Plasma density (m**-3)
        :type : float

        :param ndiv: Number of divertors with pumping (single null = 1, double null = 2 if pumping provided at both locations)
        :type : int

        :param qtorus: Gas load  from NBI (deuterons/second)
        :type : float

        :param gasld: Total D-T gas load (kg/s)
        :type : float

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean

        :returns:
            - pumpn (`float`) - Number of high vacuum pumps
            - nduct (`int`) - Number of ducts
            - dlscalc (`float`) - Duct-length equivalent for costing purposes (m)
            - mvdsh (`float`) - Mass of a single vacuum duct shield (kg)
            - dimax (`float`) -  Diameter of passage from divertor to pumping ducts (m)
        """
        k = 1.38e-23  # Boltzmann's constant (J/K)
        densh = 7900.0e0  # Density of shielding material (kg/m2)
        fsolid = 0.9e0  # Fraction of duct shielding that is solid material

        #  Pump type;
        #    ntype = 0 for turbomolecular pump (mag. bearing) with a nominal
        #              speed of 2.0 m^3/s (1.95 for N2, 1.8 for He, 1.8 for DT)
        #    ntype = 1 for compound cryopump with nominal speed of 10 m^3/s
        #              (9.0 for N2, 5.0 for He and 25. for DT)

        pfus = pfusmw * 1.0e6  # Fusion power (W)
        ntf = int(n_tf)

        #  Feed rate (gas load) of D-T into chamber (pellets + gas puffing +
        #     NBI + ...) = load from fueller + load from NBI
        #  frate (kg/s) = gasld (kg/s) + qtorus (D2/s) * 6.64e-27 (kg/D2)

        frate = gasld + qtorus * 6.64e-27

        #  Set duct shield thickness to zero for no biological shielding
        #  instead of thshldo/3.0e0

        thdsh = 0.0e0

        #  Shielding (m) between duct and TF coils is scaled from inboard shield
        #  thickness

        thcsh = thshldi / 3.0e0

        #  Multiplier to convert conductance from gas species i to nitrogen
        xmult = [1.0e0, 0.423e0, 0.378e0, 0.423e0]
        # nitrogen, D-T, helium, D-T again

        nduct = ntf * ndiv

        #  Speed of high-vacuum pumps (m^3/s)

        if vacv.ntype == 0:
            sp = [1.95, 1.8, 1.8, 1.8]
            # nitrogen, DT, helium, DT again
        else:
            sp = [9.0, 25.0, 5.0, 25.0]
            # nitrogen, DT, helium, DT again

        #  Calculate required pumping speeds

        s = []

        #  Initial pumpdown based on outgassing
        #  s(1) = net pump speed (N2) required for pumpdown to base pressure (m^3/s)
        #  area = vacuum chamber/fw area (m^2)  ;  outgassing area = 10 x area
        #  rat = outgassing rate (effective for N2) of plasma chamber surface (Pa-m/s)
        #  pbase = base pressure (Pa)

        #  Old method: area = 4.0e0 * pi*pi * r0 * aw * sqrt(0.5e0*(1.0e0 + kappa*kappa))

        area = plasma_sarea * (aw + dsol) / aw

        ogas = vacv.rat * area * 10.0e0  # Outgassing rate (Pa-m^3/s)
        s.append(ogas / vacv.pbase)

        #  Pumpdown between burns
        #  s(2) = net pump speed (DT) required for pumpdown between burns (m^3/s)
        #  tn = temperature of neutral gas in chamber (K)
        #  tdwell = dwell time between burns (s)

        pend = (
            0.5e0 * nplasma * k * vacv.tn
        )  # pressure in plasma chamber after burn (Pa)
        pstart = 0.01e0 * pend  # pressure in chamber before start of burn (Pa)

        #  Chamber volume (m^3)

        #  Old method: volume = 2.0e0 * pi*pi * r0 * aw*aw * kappa

        volume = plasma_vol * (aw + dsol) * (aw + dsol) / (aw * aw)

        #  dwell pumping options
        if (vacv.dwell_pump == 1) or (tdwell == 0):
            tpump = tv.tramp
        elif vacv.dwell_pump == 2:
            tpump = tdwell + tv.tramp
        else:
            tpump = tdwell

        s.append(volume / tpump * math.log(pend / pstart))

        #  Helium ash removal
        #  s(3) = net pump speed (He) required for helium ash removal (m^3/s)
        #  source = alpha production rate (pa - m^3/s)
        #  fhe = fraction of neutral gas in divertor chamber that is helium
        #  prdiv = pressure in divertor chamber during burn (Pa)

        source = pfus * 1.47e-09
        fhe = source / (frate * 4.985e5)
        s.append(source / vacv.prdiv / fhe)

        #  Removal of dt on steady state basis
        #  s(4) = net speed (D-T) required to remove dt at fuelling rate (m^3/s)

        s.append((frate * 4.985e5 - source) / (vacv.prdiv * (1.0e0 - fhe)))

        #  Calculate conductance of a single duct

        imax = 1
        cmax = 0.01e0
        pumpn = 1.0e0
        nflag = 0  # Control option if ducts are too small in x-sectional area
        #  = 1 if problem is identified in output, but run continues
        #  = 0 otherwise

        l1 = thshldo + thtf  # Length of passage from divertor to ducts (m)
        l2 = thshldo + 4.0e0  # Length of ducts from divertor passage to elbow (m)
        l3 = 2.0e0  # Length of ducts from elbow to hi-vac pumps (m)
        ltot = l1 + l2 + l3

        # ceff and d require initialising to small positive values; they're not
        # always overwritten in the following loop and can cause div by 0 errors
        # otherwise
        ceff = np.full(4, 1e-6)
        d = np.full(4, 1e-6)

        for i in range(4):

            sss = nduct / (
                1.0e0 / sp[i] / pumpn + 1.0e0 / cmax * xmult[i] / xmult[imax]
            )
            if sss > s[i]:
                continue
            imax = i

            ccc = 2.0e0 * s[i] / nduct
            pumpn1 = 1.0e0 / (sp[i] * (nduct / s[i] - 1.0e0 / ccc))
            pumpn2 = 1.01e0 * s[i] / (sp[i] * nduct)
            pumpn = max(pumpn, pumpn1, pumpn2)
            ceff[i] = 1.0e0 / (nduct / s[i] - 1.0e0 / (sp[i] * pumpn))

            #  Newton's method solution for duct diameter
            while True:
                d[i] = 1.0e0

                for _ in range(100):
                    a1 = (
                        0.25e0 * math.pi * d[i] * d[i]
                    )  # Area of aperture and duct (m^2)
                    a2 = 1.44e0 * a1
                    a3 = a2
                    k1 = 4.0e0 / 3.0e0 * d[i] / (l1 + 4.0e0 / 3.0e0 * d[i])
                    k2 = (
                        4.0e0
                        / 3.0e0
                        * d[i]
                        * 1.2e0
                        / (l2 + 4.0e0 / 3.0e0 * d[i] * 1.2e0)
                    )
                    k3 = (
                        4.0e0
                        / 3.0e0
                        * d[i]
                        * 1.2e0
                        / (l3 + 4.0e0 / 3.0e0 * d[i] * 1.2e0)
                    )
                    cap = 119.0e0 * a1 / xmult[i]
                    dcap = 2.0e0 * cap / d[i]
                    c1 = 119.0e0 * a1 * k1 / xmult[i]
                    dc1 = c1 / d[i] * (3.0e0 - k1)
                    c2 = 119.0e0 * a2 * k2 / xmult[i]
                    dc2 = c2 / d[i] / 1.2e0 * (3.0e0 - k2)
                    c3 = 119.0e0 * a3 * k3 / xmult[i]
                    dc3 = c3 / d[i] / 1.2e0 * (3.0e0 - k3)
                    cnew = 1.0e0 / (1.0e0 / cap + 1.0e0 / c1 + 1.0e0 / c2 + 1.0e0 / c3)
                    y = -ceff[i] + cnew
                    dy = (
                        cnew
                        * cnew
                        * (
                            dcap / cap / cap
                            + dc1 / c1 / c1
                            + dc2 / c2 / c2
                            + dc3 / c3 / c3
                        )
                    )
                    dnew = d[i] - y / dy
                    dd = abs((d[i] - dnew) / d[i])
                    d[i] = dnew
                    if dd <= 0.01e0:
                        break

                else:
                    eh.fdiags[0] = pv.powfmw
                    eh.fdiags[1] = pv.te
                    eh.report_error(124)

                theta = math.pi / ntf

                #  Area between adjacent TF coils available for pump ducts
                #  ritf = outer radius of inboard leg of TF coil (m)

                a1max = (r0 + aw - ritf - thcsh / math.tan(theta)) ** 2 * math.tan(
                    theta
                )
                d1max = math.sqrt(4.0e0 * a1max / math.pi)  # Equivalent diameter
                if a1 < a1max:
                    break

                ceff[i] = 0.9e0 * ceff[i]
                if ceff[i] <= (1.1e0 * s[i]):
                    #  Ducts are not big enough. Flag and continue.
                    nflag = 1
                    break

            cmax = ceff[i]

        pumpn = pumpn * nduct

        #  d[imax]= diameter of passage from divertor to pumping ducts (m)
        #  dout    = diameter of ducts from passage to hi-vac pumps (m)

        dout = d[imax] * 1.2e0

        #  Net pumping speeds provided by vacuum pumping system
        #  snet(1) - net pump speed (N2) provided (m^3/s)
        #  snet(2) - net pump speed (D-T) provided (m^3/s)
        #  snet(3) - net pump speed (He) provided (m^3/s)
        #  snet(4) - snet(2)
        snet = []
        for i in range(4):
            ceff1 = ceff[imax] * nduct
            snet.append(
                1.0e0
                / (1.0e0 / (ceff1 * xmult[imax] / xmult[i]) + 1.0e0 / sp[i] / pumpn)
            )

        #  If cryopumps are used then an additional pump is required
        #  for continuous operation with regeneration.

        if vacv.ntype == 1:
            pumpn = pumpn * 2.0e0

        #  Information for costing routine

        dlscalc = l1 * d[imax] ** 1.4e0 + (ltot - l1) * (d[imax] * 1.2e0) ** 1.4e0

        #  Mass of duct shielding

        arsh = (
            0.25e0 * math.pi * ((d[imax] * 1.2e0 + thdsh) ** 2 - (d[imax] * 1.2e0) ** 2)
        )
        mvdsh = arsh * (ltot - l1) * densh * fsolid

        dimax = d[imax]

        if output:

            #  Output section

            po.oheadr(self.outfile, "Vacuum System")

            po.ocmmnt(self.outfile, "Pumpdown to Base Pressure :")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile, "First wall outgassing rate (Pa m/s)", "(rat)", vacv.rat
            )
            po.ovarre(
                self.outfile, "Total outgassing load (Pa m3/s)", "(ogas)", ogas, "OP "
            )
            po.ovarre(
                self.outfile, "Base pressure required (Pa)", "(pbase)", vacv.pbase
            )
            po.ovarre(
                self.outfile, "Required N2 pump speed (m3/s)", "(s(1))", s[0], "OP "
            )
            po.ovarre(
                self.outfile,
                "N2 pump speed provided (m3/s)",
                "(snet(1))",
                snet[0],
                "OP ",
            )

            po.osubhd(self.outfile, "Pumpdown between Burns :")
            po.ovarre(
                self.outfile, "Plasma chamber volume (m3)", "(volume)", volume, "OP "
            )
            po.ovarre(
                self.outfile, "Chamber pressure after burn (Pa)", "(pend)", pend, "OP "
            )
            po.ovarre(
                self.outfile, "Chamber pressure before burn (Pa)", "(pstart)", pstart
            )
            po.ovarin(
                self.outfile,
                "Allowable pumping time switch",
                "(dwell_pump)",
                vacv.dwell_pump,
            )
            po.ovarre(self.outfile, "Dwell time between burns (s)", "(tdwell.)", tdwell)
            po.ovarre(self.outfile, "CS ramp-up time burns (s)", "(tramp.)", tv.tramp)
            po.ovarre(
                self.outfile,
                "Allowable pumping time between burns (s)",
                "(tpump)",
                tpump,
            )
            po.ovarre(
                self.outfile, "Required D-T pump speed (m3/s)", "(s(2))", s[1], "OP "
            )
            po.ovarre(
                self.outfile,
                "D-T pump speed provided (m3/s)",
                "(snet(2))",
                snet[1],
                "OP ",
            )

            po.osubhd(self.outfile, "Helium Ash Removal :")
            po.ovarre(
                self.outfile,
                "Divertor chamber gas pressure (Pa)",
                "(prdiv)",
                vacv.prdiv,
            )
            po.ovarre(
                self.outfile,
                "Helium gas fraction in divertor chamber",
                "(fhe)",
                fhe,
                "OP ",
            )
            po.ovarre(
                self.outfile, "Required helium pump speed (m3/s)", "(s(3))", s[2], "OP "
            )
            po.ovarre(
                self.outfile,
                "Helium pump speed provided (m3/s)",
                "(snet(3))",
                snet[2],
                "OP ",
            )

            po.osubhd(self.outfile, "D-T Removal at Fuelling Rate :")
            po.ovarre(self.outfile, "D-T fuelling rate (kg/s)", "(frate)", frate, "OP ")
            po.ovarre(
                self.outfile, "Required D-T pump speed (m3/s)", "(s(4))", s[3], "OP "
            )
            po.ovarre(
                self.outfile,
                "D-T pump speed provided (m3/s)",
                "(snet(4))",
                snet[3],
                "OP ",
            )

            if nflag == 1:
                po.oblnkl(self.outfile)
                po.ocmmnt(self.outfile, "Vacuum pumping ducts are space limited.")
                po.ocmmnt(self.outfile, f"Maximum duct diameter is only {d1max} m")
                po.ocmmnt(self.outfile, "Conductance is inadequate.")
                po.oblnkl(self.outfile)

            if vacv.ntype == 1:
                ipump = "cryo "
            else:
                ipump = "turbo"

            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "The vacuum pumping system size is governed by the")

            if imax == 1:
                po.ocmmnt(self.outfile, "requirements for pumpdown to base pressure.")
            elif imax == 2:
                po.ocmmnt(self.outfile, "requirements for pumpdown between burns.")
            elif imax == 3:
                po.ocmmnt(self.outfile, "requirements for helium ash removal.")
            else:
                po.ocmmnt(
                    self.outfile, "requirements for D-T removal at fuelling rate."
                )

            po.oblnkl(self.outfile)
            po.ovarin(self.outfile, "Number of large pump ducts", "(nduct)", nduct)
            po.ovarre(
                self.outfile,
                "Passage diameter, divertor to ducts (m)",
                "(d(imax))",
                d[imax],
                "OP ",
            )
            po.ovarre(self.outfile, "Passage length (m)", "(l1)", l1, "OP ")
            po.ovarre(self.outfile, "Diameter of ducts (m)", "(dout)", dout, "OP ")

            po.ovarre(
                self.outfile, "Duct length, divertor to elbow (m)", "(l2)", l2, "OP "
            )
            po.ovarre(self.outfile, "Duct length, elbow to pumps (m)", "(l3)", l3)
            po.ovarre(self.outfile, "Number of pumps", "(pumpn)", pumpn, "OP ")
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, f"The vacuum system uses {ipump} pumps.")

        return pumpn, nduct, dlscalc, mvdsh, dimax
