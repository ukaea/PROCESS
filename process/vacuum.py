import logging
import math

import numpy as np

from process import constants
from process import process_output as po
from process.blanket_library import dshellvol, eshellvol
from process.data_structure import blanket_library as blanket_library
from process.data_structure import build_variables as buv
from process.data_structure import divertor_variables as dv
from process.data_structure import ccfe_hcpb_module as ccfe_hcpb_module
from process.data_structure import fwbs_variables as fwbs_variables
from process.data_structure import physics_variables as physics_variables
from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.data_structure import vacuum_variables as vacv

logger = logging.getLogger(__name__)


class Vacuum:
    """Module containing vacuum system routines
    author: P J Knight, CCFE, Culham Science Centre

    This module contains routines for calculating the
    parameters of the vacuum system for a fusion power plant.
    """

    def __init__(self) -> None:
        self.outfile: int = constants.NOUT

    def run(self, output: bool) -> None:
        """Routine to call the vacuum module
        author: P J Knight, CCFE, Culham Science Centre

        This routine calls the main vacuum package.

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        # (should be) NBI gas load (deuterons/second)

        qtorus = 0.0e0

        #  Total fuel gas load (kg/s)
        #  2 nuclei * nucleus-pairs/sec * mass/nucleus

        # MDK Check this!!
        gasld = (
            2.0e0 * pv.molflow_plasma_fuelling_required * pv.m_fuel_amu * constants.UMASS
        )

        self.i_vacuum_pumping = vacv.i_vacuum_pumping

        # i_vacuum_pumping required to be compared to a b string
        # as this is what f2py returns
        if self.i_vacuum_pumping == "old":
            (
                pumpn,
                vacv.n_vv_vacuum_ducts,
                vacv.dlscal,
                vacv.m_vv_vacuum_duct_shield,
                vacv.dia_vv_vacuum_ducts,
            ) = self.vacuum(
                pv.p_fusion_total_mw,
                pv.rmajor,
                pv.rminor,
                0.5e0 * (buv.dr_fw_plasma_gap_inboard + buv.dr_fw_plasma_gap_outboard),
                pv.a_plasma_surface,
                pv.vol_plasma,
                buv.dr_shld_outboard,
                buv.dr_shld_inboard,
                buv.dr_tf_inboard,
                buv.rsldi - buv.dr_shld_vv_gap_inboard - buv.dr_vv_inboard,
                tfv.n_tf_coils,
                tv.t_plant_pulse_dwell,
                pv.nd_plasma_electrons_vol_avg,
                dv.n_divertors,
                qtorus,
                gasld,
                output=output,
            )
            # MDK pumpn is real: convert to integer by rounding.
            vacv.n_vac_pumps_high = math.floor(pumpn + 0.5e0)
        elif self.i_vacuum_pumping == "simple":
            vacv.n_iter_vacuum_pumps = self.vacuum_simple(output=output)
        else:
            logger.error(f"i_vacuum_pumping is invalid: {vacv.i_vacuum_pumping}")

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
        n_iter_vacuum_pumps = (
            pv.molflow_plasma_fuelling_required / vacv.molflow_vac_pumps
        )

        # Pump-down:
        # Pumping speed per pump m3/s
        pumpspeed = (
            vacv.volflow_vac_pumps_max
            * vacv.f_a_vac_pump_port_plasma_surface
            * vacv.f_volflow_vac_pumps_impedance
            * pv.a_plasma_surface
            / tfv.n_tf_coils
        )

        wallarea = (pv.a_plasma_surface / 1084.0e0) * 2000.0e0
        # Required pumping speed for pump-down
        pumpdownspeed = (
            vacv.outgasfactor * wallarea / vacv.pres_vv_chamber_base
        ) * tv.t_plant_pulse_dwell ** (-vacv.outgasindex)
        # Number of pumps required for pump-down
        npumpdown = pumpdownspeed / pumpspeed

        # Combine the two (somewhat inconsistent) models
        # Note that 'npump' can be constrained by constraint equation 63
        npump = max(n_iter_vacuum_pumps, npumpdown)

        #  Output section
        if output:
            po.oheadr(self.outfile, "Vacuum System")
            po.ovarst(
                self.outfile,
                "Switch for vacuum pumping model",
                "(i_vacuum_pumping)",
                '"' + self.i_vacuum_pumping + '"',
            )
            po.ocmmnt(
                self.outfile,
                "Simple steady-state model with comparison to ITER cryopumps",
            )
            po.ovarre(
                self.outfile,
                "Plasma fuelling rate (nucleus-pairs/s)",
                "(molflow_plasma_fuelling_required)",
                pv.molflow_plasma_fuelling_required,
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
                "(n_iter_vacuum_pumps)",
                n_iter_vacuum_pumps,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Dwell time",
                "(t_plant_pulse_dwell)",
                tv.t_plant_pulse_dwell,
            )
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
        n_tf_coils,
        t_plant_pulse_dwell,
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

        :param t_plant_pulse_dwell: Dwell time between pulses (s)
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
        #    i_vacuum_pump_type = 0 for turbomolecular pump (mag. bearing) with a nominal
        #              speed of 2.0 m^3/s (1.95 for N2, 1.8 for He, 1.8 for DT)
        #    i_vacuum_pump_type = 1 for compound cryopump with nominal speed of 10 m^3/s
        #              (9.0 for N2, 5.0 for He and 25. for DT)

        pfus = pfusmw * 1.0e6  # Fusion power (W)
        ntf = int(n_tf_coils)

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

        # nitrogen, DT, helium, DT again
        sp = (
            [1.95, 1.8, 1.8, 1.8]
            if vacv.i_vacuum_pump_type == 0
            else [9.0, 25.0, 5.0, 25.0]
        )

        #  Calculate required pumping speeds

        s = []

        #  Initial pumpdown based on outgassing
        #  s(1) = net pump speed (N2) required for pumpdown to base pressure (m^3/s)
        #  area = vacuum chamber/fw area (m^2)  ;  outgassing area = 10 x area
        #  outgrat_fw = outgassing rate (effective for N2) of plasma chamber surface (Pa-m/s)
        #  pres_vv_chamber_base = base pressure (Pa)

        #  Old method: area = 4.0e0 * pi*pi * r0 * aw * sqrt(0.5e0*(1.0e0 + kappa*kappa))

        area = plasma_sarea * (aw + dsol) / aw

        ogas = vacv.outgrat_fw * area * 10.0e0  # Outgassing rate (Pa-m^3/s)
        s.append(ogas / vacv.pres_vv_chamber_base)

        #  Pumpdown between burns
        #  s(2) = net pump speed (DT) required for pumpdown between burns (m^3/s)
        #  temp_vv_chamber_gas_burn_end = temperature of neutral gas in chamber (K)
        #  t_plant_pulse_dwell = dwell time between burns (s)

        pend = (
            0.5e0 * nplasma * k * vacv.temp_vv_chamber_gas_burn_end
        )  # pressure in plasma chamber after burn (Pa)
        pstart = 0.01e0 * pend  # pressure in chamber before start of burn (Pa)

        #  Chamber volume (m^3)

        #  Old method: volume = 2.0e0 * pi*pi * r0 * aw*aw * kappa

        volume = plasma_vol * (aw + dsol) * (aw + dsol) / (aw * aw)

        #  dwell pumping options
        if (vacv.i_vac_pump_dwell == 1) or (t_plant_pulse_dwell == 0):
            tpump = tv.t_plant_pulse_coil_precharge
        elif vacv.i_vac_pump_dwell == 2:
            tpump = t_plant_pulse_dwell + tv.t_plant_pulse_coil_precharge
        else:
            tpump = t_plant_pulse_dwell

        s.append(volume / tpump * math.log(pend / pstart))

        #  Helium ash removal
        #  s(3) = net pump speed (He) required for helium ash removal (m^3/s)
        #  source = alpha production rate (pa - m^3/s)
        #  fhe = fraction of neutral gas in divertor chamber that is helium
        #  pres_div_chamber_burn = pressure in divertor chamber during burn (Pa)

        source = pfus * 1.47e-09
        fhe = source / (frate * 4.985e5)
        s.extend(
            (
                (source / vacv.pres_div_chamber_burn / fhe),
                #  Removal of dt on steady state basis
                #  s(4) = net speed (D-T) required to remove dt at fuelling rate (m^3/s)
                (
                    (frate * 4.985e5 - source)
                    / (vacv.pres_div_chamber_burn * (1.0e0 - fhe))
                ),
            ),
        )

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
            sss = nduct / (1.0e0 / sp[i] / pumpn + 1.0e0 / cmax * xmult[i] / xmult[imax])
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
                    logger.error(
                        f"Newton's method not converging; check fusion power, te {pv.p_fusion_total_mw=} {pv.temp_plasma_electron_vol_avg_kev=}"
                    )

                theta = math.pi / ntf

                #  Area between adjacent TF coils available for pump ducts
                #  ritf = outer radius of inboard leg of TF coil (m)

                a1max = (r0 + aw - ritf - thcsh / math.tan(theta)) ** 2 * math.tan(theta)
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

        if vacv.i_vacuum_pump_type == 1:
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
                self.outfile,
                "First wall outgassing rate (Pa m/s)",
                "(outgrat_fw)",
                vacv.outgrat_fw,
            )
            po.ovarre(
                self.outfile, "Total outgassing load (Pa m3/s)", "(ogas)", ogas, "OP "
            )
            po.ovarre(
                self.outfile,
                "Base pressure required (Pa)",
                "(pres_vv_chamber_base)",
                vacv.pres_vv_chamber_base,
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
                "(i_vac_pump_dwell)",
                vacv.i_vac_pump_dwell,
            )
            po.ovarre(
                self.outfile,
                "Dwell time between burns (s)",
                "(t_plant_pulse_dwell.)",
                t_plant_pulse_dwell,
            )
            po.ovarre(
                self.outfile,
                "CS ramp-up time burns (s)",
                "(t_plant_pulse_coil_precharge.)",
                tv.t_plant_pulse_coil_precharge,
            )
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
                "(pres_div_chamber_burn)",
                vacv.pres_div_chamber_burn,
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

            i_fw_blkt_shared_coolant = (
                "cryo " if vacv.i_vacuum_pump_type == 1 else "turbo"
            )

            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "The vacuum pumping system size is governed by the")

            if imax == 1:
                po.ocmmnt(self.outfile, "requirements for pumpdown to base pressure.")
            elif imax == 2:
                po.ocmmnt(self.outfile, "requirements for pumpdown between burns.")
            elif imax == 3:
                po.ocmmnt(self.outfile, "requirements for helium ash removal.")
            else:
                po.ocmmnt(self.outfile, "requirements for D-T removal at fuelling rate.")

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
            po.ocmmnt(
                self.outfile,
                f"The vacuum system uses {i_fw_blkt_shared_coolant} pumps.",
            )

        return pumpn, nduct, dlscalc, mvdsh, dimax


class VacuumVessel:
    """Class containing vacuum vessel routines"""

    def __init__(self) -> None:
        pass

    def run(self) -> None:
        blanket_library.dz_vv_half = self.calculate_vessel_half_height(
            z_tf_inside_half=buv.z_tf_inside_half,
            dz_shld_vv_gap=buv.dz_shld_vv_gap,
            dz_vv_lower=buv.dz_vv_lower,
            n_divertors=pv.n_divertors,
            dz_blkt_upper=buv.dz_blkt_upper,
            dz_shld_upper=buv.dz_shld_upper,
            z_plasma_xpoint_upper=buv.z_plasma_xpoint_upper,
            dr_fw_plasma_gap_inboard=buv.dr_fw_plasma_gap_inboard,
            dr_fw_plasma_gap_outboard=buv.dr_fw_plasma_gap_outboard,
            dr_fw_inboard=buv.dr_fw_inboard,
            dr_fw_outboard=buv.dr_fw_outboard,
        )
        # D-shaped blanket and shield
        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            (
                blanket_library.vol_vv_inboard,
                blanket_library.vol_vv_outboard,
                fwbs_variables.vol_vv,
            ) = self.calculate_dshaped_vessel_volumes(
                rsldi=buv.rsldi,
                rsldo=buv.rsldo,
                dz_vv_half=blanket_library.dz_vv_half,
                dr_vv_inboard=buv.dr_vv_inboard,
                dr_vv_outboard=buv.dr_vv_outboard,
                dz_vv_upper=buv.dz_vv_upper,
                dz_vv_lower=buv.dz_vv_lower,
            )
        else:
            (
                blanket_library.vol_vv_inboard,
                blanket_library.vol_vv_outboard,
                fwbs_variables.vol_vv,
            ) = self.calculate_elliptical_vessel_volumes(
                rmajor=pv.rmajor,
                rminor=pv.rminor,
                triang=pv.triang,
                rsldi=buv.rsldi,
                rsldo=buv.rsldo,
                dz_vv_half=blanket_library.dz_vv_half,
                dr_vv_inboard=buv.dr_vv_inboard,
                dr_vv_outboard=buv.dr_vv_outboard,
                dz_vv_upper=buv.dz_vv_upper,
                dz_vv_lower=buv.dz_vv_lower,
            )

        # Apply vacuum vessel coverage factor
        # moved from dshaped_* and elliptical_* to keep coverage factor
        # changes in the same location.
        fwbs_variables.vol_vv = fwbs_variables.fvoldw * fwbs_variables.vol_vv

        ccfe_hcpb_module.vv_density = fwbs_variables.m_vv / fwbs_variables.vol_vv

    def calculate_vessel_half_height(
        self,
        z_tf_inside_half: float,
        dz_shld_vv_gap: float,
        dz_vv_lower: float,
        n_divertors: int,
        dz_blkt_upper: float,
        dz_shld_upper: float,
        z_plasma_xpoint_upper: float,
        dr_fw_plasma_gap_inboard: float,
        dr_fw_plasma_gap_outboard: float,
        dr_fw_inboard: float,
        dr_fw_outboard: float,
    ) -> float:
        """Calculate vacuum vessel internal half-height (m)"""

        z_bottom = z_tf_inside_half - dz_shld_vv_gap - dz_vv_lower

        # Calculate component internal upper half-height (m)
        # If a double null machine then symmetric
        if n_divertors == 2:
            z_top = z_bottom
        else:
            z_top = z_plasma_xpoint_upper + 0.5 * (
                dr_fw_plasma_gap_inboard
                + dr_fw_plasma_gap_outboard
                + dr_fw_inboard
                + dr_fw_outboard
            )
        z_top = z_top + dz_blkt_upper + dz_shld_upper

        # Average of top and bottom (m)
        return 0.5 * (z_top + z_bottom)

    def calculate_dshaped_vessel_volumes(
        self,
        rsldi: float,
        rsldo: float,
        dz_vv_half: float,
        dr_vv_inboard: float,
        dr_vv_outboard: float,
        dz_vv_upper: float,
        dz_vv_lower: float,
    ) -> tuple[float, float, float]:
        """Calculate volumes of D-shaped vacuum vessel segments"""

        r_1 = rsldi
        r_2 = rsldo - r_1

        (
            vol_vv_inboard,
            vol_vv_outboard,
            vol_vv,
        ) = dshellvol(
            rmajor=r_1,
            rminor=r_2,
            zminor=dz_vv_half,
            drin=dr_vv_inboard,
            drout=dr_vv_outboard,
            dz=(dz_vv_upper + dz_vv_lower) / 2,
        )

        return vol_vv_inboard, vol_vv_outboard, vol_vv

    def calculate_elliptical_vessel_volumes(
        self,
        rmajor: float,
        rminor: float,
        triang: float,
        rsldi: float,
        rsldo: float,
        dz_vv_half: float,
        dr_vv_inboard: float,
        dr_vv_outboard: float,
        dz_vv_upper: float,
        dz_vv_lower: float,
    ) -> tuple[float, float, float]:
        # Major radius to centre of inboard and outboard ellipses (m)
        # (coincident in radius with top of plasma)
        r_1 = rmajor - rminor * triang

        # Calculate distance between r1 and outer edge of inboard ...
        # ... section (m)
        r_2 = r_1 - rsldi
        r_3 = rsldo - r_1

        (
            vol_vv_inboard,
            vol_vv_outboard,
            vol_vv,
        ) = eshellvol(
            r_1,
            r_2,
            r_3,
            dz_vv_half,
            dr_vv_inboard,
            dr_vv_outboard,
            (dz_vv_upper + dz_vv_lower) / 2,
        )
        return vol_vv_inboard, vol_vv_outboard, vol_vv
