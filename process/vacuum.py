import logging
import math

import numpy as np

from process import constants, process_output
from process import process_output as po
from process.blanket_library import dshellvol, eshellvol
from process.data_structure import (
    blanket_library,
    build_variables,
    divertor_variables,
    fwbs_variables,
    physics_variables,
    tfcoil_variables,
    times_variables,
    vacuum_variables,
)

logger = logging.getLogger(__name__)


class Vacuum:
    """Module containing vacuum system routines

    This module contains routines for calculating the
    parameters of the vacuum system for a fusion power plant.
    """

    def __init__(self):
        self.outfile: int = constants.NOUT

    def run(self, output: bool):
        """Routine to call the vacuum module
        This routine calls the main vacuum package.

        Parameters
        ----------
        output:
            indicate whether output should be written to the output file, or not

        """
        # (should be) NBI gas load (deuterons/second)

        qtorus = 0.0e0

        #  Total fuel gas load (kg/s)
        #  2 nuclei * nucleus-pairs/sec * mass/nucleus

        # MDK Check this!!
        gasld = (
            2.0e0
            * physics_variables.molflow_plasma_fuelling_required
            * physics_variables.m_fuel_amu
            * constants.UMASS
        )

        self.i_vacuum_pumping = vacuum_variables.i_vacuum_pumping

        # i_vacuum_pumping required to be compared to a b string
        # as this is what f2py returns
        if self.i_vacuum_pumping == "old":
            (
                pumpn,
                vacuum_variables.n_vv_vacuum_ducts,
                vacuum_variables.dlscal,
                vacuum_variables.m_vv_vacuum_duct_shield,
                vacuum_variables.dia_vv_vacuum_ducts,
            ) = self.vacuum(
                physics_variables.p_fusion_total_mw,
                physics_variables.rmajor,
                physics_variables.rminor,
                0.5e0
                * (
                    build_variables.dr_fw_plasma_gap_inboard
                    + build_variables.dr_fw_plasma_gap_outboard
                ),
                physics_variables.a_plasma_surface,
                physics_variables.vol_plasma,
                build_variables.dr_shld_outboard,
                build_variables.dr_shld_inboard,
                build_variables.dr_tf_inboard,
                build_variables.r_shld_inboard_inner
                - build_variables.dr_shld_vv_gap_inboard
                - build_variables.dr_vv_inboard,
                tfcoil_variables.n_tf_coils,
                times_variables.t_plant_pulse_dwell,
                physics_variables.nd_plasma_electrons_vol_avg,
                divertor_variables.n_divertors,
                qtorus,
                gasld,
                output=output,
            )
            # MDK pumpn is real: convert to integer by rounding.
            vacuum_variables.n_vac_pumps_high = math.floor(pumpn + 0.5e0)
        elif self.i_vacuum_pumping == "simple":
            vacuum_variables.n_iter_vacuum_pumps = self.vacuum_simple(output=output)
        else:
            logger.error(
                f"i_vacuum_pumping is invalid: {vacuum_variables.i_vacuum_pumping}"
            )

    def vacuum_simple(self, output) -> float:
        """Simple model of vacuum pumping system


        Parameters
        ----------
        output : boolean

        Returns
        -------
        npump:
            number of pumps for pumpdown and steady-state
            indicate whether output should be written to the output file, or not
        """

        # Steady-state model (super simple)
        # One ITER torus cryopump has a throughput of 50 Pa m3/s = 1.2155e+22 molecules/s
        # Issue #304
        n_iter_vacuum_pumps = (
            physics_variables.molflow_plasma_fuelling_required
            / vacuum_variables.molflow_vac_pumps
        )

        # Pump-down:
        # Pumping speed per pump m3/s
        pumpspeed = (
            vacuum_variables.volflow_vac_pumps_max
            * vacuum_variables.f_a_vac_pump_port_plasma_surface
            * vacuum_variables.f_volflow_vac_pumps_impedance
            * physics_variables.a_plasma_surface
            / tfcoil_variables.n_tf_coils
        )

        wallarea = (physics_variables.a_plasma_surface / 1084.0e0) * 2000.0e0
        # Required pumping speed for pump-down
        pumpdownspeed = (
            vacuum_variables.outgasfactor
            * wallarea
            / vacuum_variables.pres_vv_chamber_base
        ) * times_variables.t_plant_pulse_dwell ** (-vacuum_variables.outgasindex)
        # Number of pumps required for pump-down
        npumpdown = pumpdownspeed / pumpspeed

        # Combine the two (somewhat inconsistent) models
        # Note that 'npump' can be constrained by constraint equation 63
        npump = max(n_iter_vacuum_pumps, npumpdown)

        #  Output section
        if output:
            process_output.oheadr(self.outfile, "Vacuum System")
            process_output.ovarst(
                self.outfile,
                "Switch for vacuum pumping model",
                "(i_vacuum_pumping)",
                '"' + self.i_vacuum_pumping + '"',
            )
            process_output.ocmmnt(
                self.outfile,
                "Simple steady-state model with comparison to ITER cryopumps",
            )
            process_output.ovarre(
                self.outfile,
                "Plasma fuelling rate (nucleus-pairs/s)",
                "(molflow_plasma_fuelling_required)",
                physics_variables.molflow_plasma_fuelling_required,
                "OP ",
            )
            process_output.ocmmnt(
                self.outfile, "Number of high vacuum pumps, each with the throughput"
            )
            process_output.ocmmnt(
                self.outfile,
                " of one ITER cryopump (50 Pa m3 s-1 = 1.2e+22 molecules/s),",
            )
            process_output.ovarre(
                self.outfile,
                " all operating at the same time",
                "(n_iter_vacuum_pumps)",
                n_iter_vacuum_pumps,
                "OP ",
            )

            process_output.ovarre(
                self.outfile,
                "Dwell time",
                "(t_plant_pulse_dwell)",
                times_variables.t_plant_pulse_dwell,
            )
            process_output.ovarre(
                self.outfile,
                "Number of pumps required for pump-down",
                "(npumpdown)",
                npumpdown,
                "OP ",
            )
            process_output.ovarre(
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

        Parameters
        ----------
        pfusmw : float
            Fusion power (MW)
        r0 : float
            Major radius (m)
        aw : float
            Minor radius (m)
        dsol :
            Scrape-off layer average width (m)
        plasma_sarea :
            Plasma surface area (m2)
        plasma_vol :
            Plasma volume (m3)
        thshldo :
            Outboard shield thickness (m)
        thshldi :
            Inboard shield thickness (m)
        thtf :
            TF coil thickness (m)
        ritf :
            Radius of inboard TF leg point nearest plasma (m)
        n_tf_coils :
            Number of TF coils
        t_plant_pulse_dwell :
            Dwell time between pulses (s)
        nplasma :
            Plasma density (m**-3)
        ndiv :
            Number of divertors with pumping (single null = 1, double null = 2 if pumping provided at both locations)
        qtorus :
            Gas load  from NBI (deuterons/second)
        gasld :
            Total D-T gas load (kg/s)
        output : boolean
            indicate whether output should be written to the output file, or not


        Returns
        -------
        :
            pumpn (`float`) - Number of high vacuum pumps
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
            if vacuum_variables.i_vacuum_pump_type == 0
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

        ogas = vacuum_variables.outgrat_fw * area * 10.0e0  # Outgassing rate (Pa-m^3/s)
        s.append(ogas / vacuum_variables.pres_vv_chamber_base)

        #  Pumpdown between burns
        #  s(2) = net pump speed (DT) required for pumpdown between burns (m^3/s)
        #  temp_vv_chamber_gas_burn_end = temperature of neutral gas in chamber (K)
        #  t_plant_pulse_dwell = dwell time between burns (s)

        pend = (
            0.5e0 * nplasma * k * vacuum_variables.temp_vv_chamber_gas_burn_end
        )  # pressure in plasma chamber after burn (Pa)
        pstart = 0.01e0 * pend  # pressure in chamber before start of burn (Pa)

        #  Chamber volume (m^3)

        #  Old method: volume = 2.0e0 * pi*pi * r0 * aw*aw * kappa

        volume = plasma_vol * (aw + dsol) * (aw + dsol) / (aw * aw)

        #  dwell pumping options
        if (vacuum_variables.i_vac_pump_dwell == 1) or (t_plant_pulse_dwell == 0):
            tpump = times_variables.t_plant_pulse_coil_precharge
        elif vacuum_variables.i_vac_pump_dwell == 2:
            tpump = t_plant_pulse_dwell + times_variables.t_plant_pulse_coil_precharge
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
                (source / vacuum_variables.pres_div_chamber_burn / fhe),
                #  Removal of dt on steady state basis
                #  s(4) = net speed (D-T) required to remove dt at fuelling rate (m^3/s)
                (
                    (frate * 4.985e5 - source)
                    / (vacuum_variables.pres_div_chamber_burn * (1.0e0 - fhe))
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
                        f"Newton's method not converging; check fusion power, te {physics_variables.p_fusion_total_mw=} {physics_variables.temp_plasma_electron_vol_avg_kev=}"
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

        if vacuum_variables.i_vacuum_pump_type == 1:
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

            process_output.oheadr(self.outfile, "Vacuum System")

            process_output.ocmmnt(self.outfile, "Pumpdown to Base Pressure :")
            process_output.oblnkl(self.outfile)
            process_output.ovarre(
                self.outfile,
                "First wall outgassing rate (Pa m/s)",
                "(outgrat_fw)",
                vacuum_variables.outgrat_fw,
            )
            process_output.ovarre(
                self.outfile, "Total outgassing load (Pa m3/s)", "(ogas)", ogas, "OP "
            )
            process_output.ovarre(
                self.outfile,
                "Base pressure required (Pa)",
                "(pres_vv_chamber_base)",
                vacuum_variables.pres_vv_chamber_base,
            )
            process_output.ovarre(
                self.outfile, "Required N2 pump speed (m3/s)", "(s(1))", s[0], "OP "
            )
            process_output.ovarre(
                self.outfile,
                "N2 pump speed provided (m3/s)",
                "(snet(1))",
                snet[0],
                "OP ",
            )

            process_output.osubhd(self.outfile, "Pumpdown between Burns :")
            process_output.ovarre(
                self.outfile, "Plasma chamber volume (m3)", "(volume)", volume, "OP "
            )
            process_output.ovarre(
                self.outfile, "Chamber pressure after burn (Pa)", "(pend)", pend, "OP "
            )
            process_output.ovarre(
                self.outfile, "Chamber pressure before burn (Pa)", "(pstart)", pstart
            )
            process_output.ovarin(
                self.outfile,
                "Allowable pumping time switch",
                "(i_vac_pump_dwell)",
                vacuum_variables.i_vac_pump_dwell,
            )
            process_output.ovarre(
                self.outfile,
                "Dwell time between burns (s)",
                "(t_plant_pulse_dwell.)",
                t_plant_pulse_dwell,
            )
            process_output.ovarre(
                self.outfile,
                "CS ramp-up time burns (s)",
                "(t_plant_pulse_coil_precharge.)",
                times_variables.t_plant_pulse_coil_precharge,
            )
            process_output.ovarre(
                self.outfile,
                "Allowable pumping time between burns (s)",
                "(tpump)",
                tpump,
            )
            process_output.ovarre(
                self.outfile, "Required D-T pump speed (m3/s)", "(s(2))", s[1], "OP "
            )
            process_output.ovarre(
                self.outfile,
                "D-T pump speed provided (m3/s)",
                "(snet(2))",
                snet[1],
                "OP ",
            )

            process_output.osubhd(self.outfile, "Helium Ash Removal :")
            process_output.ovarre(
                self.outfile,
                "Divertor chamber gas pressure (Pa)",
                "(pres_div_chamber_burn)",
                vacuum_variables.pres_div_chamber_burn,
            )
            process_output.ovarre(
                self.outfile,
                "Helium gas fraction in divertor chamber",
                "(fhe)",
                fhe,
                "OP ",
            )
            process_output.ovarre(
                self.outfile, "Required helium pump speed (m3/s)", "(s(3))", s[2], "OP "
            )
            process_output.ovarre(
                self.outfile,
                "Helium pump speed provided (m3/s)",
                "(snet(3))",
                snet[2],
                "OP ",
            )

            process_output.osubhd(self.outfile, "D-T Removal at Fuelling Rate :")
            process_output.ovarre(
                self.outfile, "D-T fuelling rate (kg/s)", "(frate)", frate, "OP "
            )
            process_output.ovarre(
                self.outfile, "Required D-T pump speed (m3/s)", "(s(4))", s[3], "OP "
            )
            process_output.ovarre(
                self.outfile,
                "D-T pump speed provided (m3/s)",
                "(snet(4))",
                snet[3],
                "OP ",
            )

            if nflag == 1:
                process_output.oblnkl(self.outfile)
                process_output.ocmmnt(
                    self.outfile, "Vacuum pumping ducts are space limited."
                )
                process_output.ocmmnt(
                    self.outfile, f"Maximum duct diameter is only {d1max} m"
                )
                process_output.ocmmnt(self.outfile, "Conductance is inadequate.")
                process_output.oblnkl(self.outfile)

            i_fw_blkt_shared_coolant = (
                "cryo " if vacuum_variables.i_vacuum_pump_type == 1 else "turbo"
            )

            process_output.oblnkl(self.outfile)
            process_output.ocmmnt(
                self.outfile, "The vacuum pumping system size is governed by the"
            )

            if imax == 1:
                process_output.ocmmnt(
                    self.outfile, "requirements for pumpdown to base pressure."
                )
            elif imax == 2:
                process_output.ocmmnt(
                    self.outfile, "requirements for pumpdown between burns."
                )
            elif imax == 3:
                process_output.ocmmnt(
                    self.outfile, "requirements for helium ash removal."
                )
            else:
                process_output.ocmmnt(
                    self.outfile, "requirements for D-T removal at fuelling rate."
                )

            process_output.oblnkl(self.outfile)
            process_output.ovarin(
                self.outfile, "Number of large pump ducts", "(nduct)", nduct
            )
            process_output.ovarre(
                self.outfile,
                "Passage diameter, divertor to ducts (m)",
                "(d(imax))",
                d[imax],
                "OP ",
            )
            process_output.ovarre(self.outfile, "Passage length (m)", "(l1)", l1, "OP ")
            process_output.ovarre(
                self.outfile, "Diameter of ducts (m)", "(dout)", dout, "OP "
            )

            process_output.ovarre(
                self.outfile, "Duct length, divertor to elbow (m)", "(l2)", l2, "OP "
            )
            process_output.ovarre(
                self.outfile, "Duct length, elbow to pumps (m)", "(l3)", l3
            )
            process_output.ovarre(
                self.outfile, "Number of pumps", "(pumpn)", pumpn, "OP "
            )
            process_output.oblnkl(self.outfile)
            process_output.ocmmnt(
                self.outfile,
                f"The vacuum system uses {i_fw_blkt_shared_coolant} pumps.",
            )

        return pumpn, nduct, dlscalc, mvdsh, dimax


class VacuumVessel:
    """Class containing vacuum vessel routines"""

    def __init__(self):
        self.outfile = constants.NOUT

    def run(self):
        blanket_library.dz_vv_half = self.calculate_vessel_half_height(
            z_tf_inside_half=build_variables.z_tf_inside_half,
            dz_shld_vv_gap=build_variables.dz_shld_vv_gap,
            dz_vv_lower=build_variables.dz_vv_lower,
            n_divertors=divertor_variables.n_divertors,
            dz_blkt_upper=build_variables.dz_blkt_upper,
            dz_shld_upper=build_variables.dz_shld_upper,
            z_plasma_xpoint_upper=build_variables.z_plasma_xpoint_upper,
            dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
            dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            dr_fw_inboard=build_variables.dr_fw_inboard,
            dr_fw_outboard=build_variables.dr_fw_outboard,
        )
        # D-shaped blanket and shield
        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            (
                blanket_library.vol_vv_inboard,
                blanket_library.vol_vv_outboard,
                fwbs_variables.vol_vv,
            ) = self.calculate_dshaped_vessel_volumes(
                r_shld_inboard_inner=build_variables.r_shld_inboard_inner,
                r_shld_outboard_outer=build_variables.r_shld_outboard_outer,
                dz_vv_half=blanket_library.dz_vv_half,
                dr_vv_inboard=build_variables.dr_vv_inboard,
                dr_vv_outboard=build_variables.dr_vv_outboard,
                dz_vv_upper=build_variables.dz_vv_upper,
                dz_vv_lower=build_variables.dz_vv_lower,
            )
        else:
            (
                blanket_library.vol_vv_inboard,
                blanket_library.vol_vv_outboard,
                fwbs_variables.vol_vv,
            ) = self.calculate_elliptical_vessel_volumes(
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                triang=physics_variables.triang,
                r_shld_inboard_inner=build_variables.r_shld_inboard_inner,
                r_shld_outboard_outer=build_variables.r_shld_outboard_outer,
                dz_vv_half=blanket_library.dz_vv_half,
                dr_vv_inboard=build_variables.dr_vv_inboard,
                dr_vv_outboard=build_variables.dr_vv_outboard,
                dz_vv_upper=build_variables.dz_vv_upper,
                dz_vv_lower=build_variables.dz_vv_lower,
            )

        # Apply vacuum vessel coverage factor
        # moved from dshaped_* and elliptical_* to keep coverage factor
        # changes in the same location.
        fwbs_variables.vol_vv = fwbs_variables.fvoldw * fwbs_variables.vol_vv

        # Vacuum vessel mass (kg)
        fwbs_variables.m_vv = fwbs_variables.vol_vv * fwbs_variables.den_steel

    @staticmethod
    def calculate_vessel_half_height(
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
        """Calculate vacuum vessel internal half-height (m)

        Parameters
        ----------
        z_tf_inside_half:

        dz_shld_vv_gap:

        dz_vv_lower:

        n_divertors: int :

        dz_blkt_upper:

        dz_shld_upper:

        z_plasma_xpoint_upper:

        dr_fw_plasma_gap_inboard:

        dr_fw_plasma_gap_outboard:

        dr_fw_inboard:

        dr_fw_outboard:

        """

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

    @staticmethod
    def calculate_dshaped_vessel_volumes(
        r_shld_inboard_inner: float,
        r_shld_outboard_outer: float,
        dz_vv_half: float,
        dr_vv_inboard: float,
        dr_vv_outboard: float,
        dz_vv_upper: float,
        dz_vv_lower: float,
    ) -> tuple[float, float, float]:
        """Calculate volumes of D-shaped vacuum vessel segments

        Parameters
        ----------
        r_shld_inboard_inner:

        r_shld_outboard_outer:

        dz_vv_half:

        dr_vv_inboard:

        dr_vv_outboard:

        dz_vv_upper:

        dz_vv_lower:

        """

        r_1 = r_shld_inboard_inner
        r_2 = r_shld_outboard_outer - r_1

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

    @staticmethod
    def calculate_elliptical_vessel_volumes(
        rmajor: float,
        rminor: float,
        triang: float,
        r_shld_inboard_inner: float,
        r_shld_outboard_outer: float,
        dz_vv_half: float,
        dr_vv_inboard: float,
        dr_vv_outboard: float,
        dz_vv_upper: float,
        dz_vv_lower: float,
    ) -> tuple[float, float, float]:
        """Calculate volumes of elliptical vacuum vessel segments

        Parameters
        ----------
        rmajor:

        rminor:

        triang:

        r_shld_inboard_inner:

        r_shld_outboard_outer:

        dz_vv_half:

        dr_vv_inboard:

        dr_vv_outboard:

        dz_vv_upper:

        dz_vv_lower:

        """
        # Major radius to centre of inboard and outboard ellipses (m)
        # (coincident in radius with top of plasma)
        r_1 = rmajor - rminor * triang

        # Calculate distance between r1 and outer edge of inboard ...
        # ... section (m)
        r_2 = r_1 - r_shld_inboard_inner
        r_3 = r_shld_outboard_outer - r_1

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

    def output_vv_areas_and_volumes(self):
        """Output shield areas and volumes to log."""

        po.oheadr(self.outfile, "Vacuum Vessel Areas and Volumes")

        po.ovarrf(
            self.outfile,
            "Volume of inboard vacuum vessel (m^3)",
            "(vol_vv_inboard)",
            blanket_library.vol_vv_inboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Volume of outboard vacuum vessel (m^3)",
            "(vol_vv_outboard)",
            blanket_library.vol_vv_outboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Total volume of vacuum vessel (m^3)",
            "(vol_vv)",
            fwbs_variables.vol_vv,
            "OP ",
        )
