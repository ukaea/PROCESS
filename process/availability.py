import logging
import math

from scipy.special import comb as combinations

from process import constants
from process import process_output as po
from process.data_structure import constraint_variables as ctv
from process.data_structure import cost_variables as cv
from process.data_structure import divertor_variables as dv
from process.data_structure import fwbs_variables as fwbsv
from process.data_structure import ife_variables as ifev
from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.data_structure import vacuum_variables as vacv
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)

DAY_SECONDS = 60 * 60 * 24
# Number of seconds in a day [s]

DAYS_IN_YEAR = 365.25
# Number of days in a year

YEAR_SECONDS = DAY_SECONDS * DAYS_IN_YEAR
# Number of seconds in a year [s]

AVAILABILITY_MODELS = {
    0: "Input value for f_t_plant_available",
    1: "Ward and Taylor model (1999)",
    2: "Morris model (2015)",
    3: "ST model (2023)",
}


class Availability:
    """Module containing plant availability routines


    This module contains routines for calculating the
    plant availability and component lifetimes for a fusion power plant.
    """

    def __init__(self) -> None:
        self.outfile = constants.NOUT  # output file unit

    def run(self, output: bool = False):
        """Run appropriate availability model

        Availability switch values
        No.  |  model
        ---- | ------
        0    |  Input value for f_t_plant_available
        1    |  Ward and Taylor model (1999)
        2    |  Morris model (2015)
        3    |  ST model (2023)

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not (default = False)
        output: bool :
             (Default value = False)
        """

        if cv.i_plant_availability == 3:
            if pv.itart != 1:
                raise ProcessValueError(
                    f"{cv.i_plant_availability=} is for a Spherical Tokamak. Please set itart=1 to use this model."
                )
            self.avail_st(output)  # ST model (2023)
        elif cv.i_plant_availability == 2:
            self.avail_2(output)  # Morris model (2015)
        else:
            self.avail(output)  # Taylor and Ward model (1999)

    def avail(self, output: bool):
        """Routine to calculate component lifetimes and the overall plant availability


        This routine calculates the component lifetimes and the overall
        plant availability.
        F/PL/PJK/PROCESS/CODE/043

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :

        """

        # Full power lifetime (in years)
        if ifev.ife != 1:
            # Calculate DPA per FPY - based on neutronics-derived fusion power relation to DEMO blanket lifetime provided by Matti Coleman
            # Detailed and cited in T. Franke 2020, "The EU DEMO equatorial outboard limiter — Design and port integration concept"
            # https://www.sciencedirect.com/science/article/pii/S0920379620301952#bib0075
            # Scaling w.r.t. fusion power drops out a large number of factors relating to neutronics, such as:
            # - the actual neutron flux
            # - the geometry and material composition leading to the neutron flux at the EUROfer FW OMP
            # - the neutron energy spectrum
            # - all of the above and more leading to the dpa/fpy in EUROfer at the FW OMP
            # About a relatively "constant" reference point, we can reasonably assume they all equal to 1.0.
            ref_fusion_power = 2.0e3  # (MW) fusion power for EU-DEMO
            f_scale = pv.p_fusion_total_mw / ref_fusion_power
            ref_dpa_fpy = (
                10.0e0  # dpa per fpy from T. Franke 2020 states up to 10 dpa/FPY
            )
            dpa_fpy = f_scale * ref_dpa_fpy

            # First wall / blanket lifetime (years)
            # TODO MDK Do this calculation whatever the value of blktmodel (whatever that is)
            # For some reason life_fw_fpy is not always calculated, so ignore it if it is still zero.
            if fwbsv.life_fw_fpy < 0.0001e0:
                # Calculate blanket lifetime using neutron fluence model (ibkt_life=0)
                # or DEMO fusion power model (ibkt_life=1)
                if cv.ibkt_life == 0:
                    fwbsv.life_blkt_fpy = (
                        cv.life_plant
                        if pv.pflux_fw_neutron_mw == 0.0
                        else min(
                            (cv.abktflnc / pv.pflux_fw_neutron_mw),
                            cv.life_plant,
                        )
                    )
                else:
                    fwbsv.life_blkt_fpy = min(
                        cv.life_dpa / dpa_fpy, cv.life_plant
                    )  # DEMO
            else:
                if cv.ibkt_life == 0:
                    fwbsv.life_blkt_fpy = min(
                        fwbsv.life_fw_fpy,
                        cv.abktflnc / pv.pflux_fw_neutron_mw,
                        cv.life_plant,
                    )
                else:
                    fwbsv.life_blkt_fpy = min(
                        fwbsv.life_fw_fpy, cv.life_dpa / dpa_fpy, cv.life_plant
                    )  # DEMO

            # TODO Issue #834
            # Add a test for pflux_div_heat_load_mw=0
            if dv.pflux_div_heat_load_mw < 1.0e-10:
                dv.pflux_div_heat_load_mw = 1.0e-10

            # Divertor lifetime (years)
            cv.life_div_fpy = self.divertor_lifetime()

            # Centrepost lifetime (years) (ST machines only)
            if pv.itart == 1:
                cv.cplife = self.cp_lifetime()

        # Plant Availability (i_plant_availability=0,1)

        # Calculate the number of fusion cycles for a given blanket lifetime
        pulse_fpy = tv.t_plant_pulse_total / YEAR_SECONDS
        cv.bktcycles = (fwbsv.life_blkt_fpy / pulse_fpy) + 1

        # if i_plant_availability = 0 use input value for f_t_plant_available

        # Taylor and Ward 1999 model (i_plant_availability=1)
        if cv.i_plant_availability == 1:
            # Which component has the shorter life?
            if cv.life_div_fpy < fwbsv.life_blkt_fpy:
                ld = cv.life_div_fpy
                lb = fwbsv.life_blkt_fpy
                td = cv.t_div_replace_yrs
            else:
                ld = fwbsv.life_blkt_fpy
                lb = cv.life_div_fpy
                td = cv.t_blkt_replace_yrs

            # Number of outages between each combined outage
            n = math.ceil(lb / ld) - 1

            # Planned unavailability
            uplanned = (n * td + cv.tcomrepl) / ((n + 1) * ld + (n * td + cv.tcomrepl))

            # Unplanned unavailability
            # Rather than simply summing the individual terms, the following protects
            # against the total availability becoming zero or negative

            uutot = cv.uubop  # balance of plant
            uutot = uutot + (1.0e0 - uutot) * cv.uucd  # current drive
            uutot = uutot + (1.0e0 - uutot) * cv.uudiv  # divertor
            uutot = uutot + (1.0e0 - uutot) * cv.uufuel  # fuel system
            uutot = uutot + (1.0e0 - uutot) * cv.uufw  # first wall + blanket
            uutot = uutot + (1.0e0 - uutot) * cv.uumag  # magnets
            uutot = uutot + (1.0e0 - uutot) * cv.uuves  # vacuum vessel

            # Total availability
            cv.f_t_plant_available = 1.0e0 - (uplanned + uutot - (uplanned * uutot))

        # Capacity factor
        # Using the amount of time burning for a given pulse cycle
        cv.cpfact = cv.f_t_plant_available * (
            tv.t_plant_pulse_burn / tv.t_plant_pulse_total
        )

        # Modify lifetimes to take account of the availability
        if ifev.ife != 1:
            # First wall / blanket
            if fwbsv.life_blkt_fpy < cv.life_plant:
                fwbsv.life_blkt_fpy = min(
                    fwbsv.life_blkt_fpy / cv.f_t_plant_available, cv.life_plant
                )

            # Divertor
            if cv.life_div_fpy < cv.life_plant:
                cv.life_div_fpy = min(
                    cv.life_div_fpy / cv.f_t_plant_available, cv.life_plant
                )

            # Centrepost
            if pv.itart == 1 and cv.cplife < cv.life_plant:
                cv.cplife = min(cv.cplife / cv.f_t_plant_available, cv.life_plant)

        # Current drive system lifetime (assumed equal to first wall and blanket lifetime)
        cv.life_hcd_fpy = fwbsv.life_blkt_fpy

        # Output section
        if output:
            po.oheadr(self.outfile, "Plant Availability")
            if fwbsv.blktmodel == 0:
                po.ovarre(
                    self.outfile,
                    "Allowable blanket neutron fluence (MW-yr/m2)",
                    "(abktflnc)",
                    cv.abktflnc,
                )

            po.ovarre(
                self.outfile,
                "Allowable divertor heat fluence (MW-yr/m2)",
                "(adivflnc)",
                cv.adivflnc,
            )
            po.ovarre(
                self.outfile,
                "First wall / blanket lifetime (years)",
                "(life_blkt_fpy)",
                fwbsv.life_blkt_fpy,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Divertor lifetime (years)",
                "(life_div_fpy)",
                cv.life_div_fpy,
                "OP ",
            )

            if pv.itart == 1:
                po.ovarre(
                    self.outfile,
                    "Centrepost lifetime (years)",
                    "(cplife)",
                    cv.cplife,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Heating/CD system lifetime (years)",
                "(life_hcd_fpy)",
                cv.life_hcd_fpy,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total plant lifetime (years)",
                "(life_plant)",
                cv.life_plant,
            )

            if cv.i_plant_availability == 1:
                if cv.life_div_fpy < fwbsv.life_blkt_fpy:
                    po.ovarre(
                        self.outfile,
                        "Time needed to replace divertor (years)",
                        "(t_div_replace_yrs)",
                        cv.t_div_replace_yrs,
                    )
                else:
                    po.ovarre(
                        self.outfile,
                        "Time needed to replace blanket (years)",
                        "(t_blkt_replace_yrs)",
                        cv.t_blkt_replace_yrs,
                    )

                po.ovarre(
                    self.outfile,
                    "Time needed to replace blkt + div (years)",
                    "(tcomrepl)",
                    cv.tcomrepl,
                )
                po.ovarre(
                    self.outfile,
                    "Planned unavailability fraction",
                    "(uplanned)",
                    uplanned,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Unplanned unavailability fraction",
                    "(uutot)",
                    uutot,
                    "OP ",
                )

            if cv.i_plant_availability == 0:
                po.ovarre(
                    self.outfile,
                    "Total plant availability fraction",
                    "(f_t_plant_available)",
                    cv.f_t_plant_available,
                )
                po.ovarre(
                    self.outfile,
                    "Number of fusion cycles to reach allowable fw/blanket DPA",
                    "(bktcycles)",
                    cv.bktcycles,
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Total plant availability fraction",
                    "(f_t_plant_available)",
                    cv.f_t_plant_available,
                    "OP ",
                )

    def avail_2(self, output: bool):
        """Routine to calculate component lifetimes and the overall plant availability


        This routine calculates the component lifetimes and the overall
        plant availability using an updated model linked to the 2014 EUROfusion
        RAMI task
        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :

        """

        # Plant Availability

        # Planned unavailability

        u_planned = self.calc_u_planned(output)

        # Operational time (years)
        cv.t_plant_operational_total_yrs = cv.life_plant * (1.0e0 - u_planned)

        # Un-planned unavailability

        # Magnets
        u_unplanned_magnets = self.calc_u_unplanned_magnets(output)

        # Divertor
        u_unplanned_div = self.calc_u_unplanned_divertor(output)

        # First wall and blanket
        u_unplanned_fwbs = self.calc_u_unplanned_fwbs(output)

        # Balance of plant
        u_unplanned_bop = self.calc_u_unplanned_bop(output)

        # Heating and current drive
        u_unplanned_hcd = self.calc_u_unplanned_hcd()

        # Vacuum systems

        # Number of redundant pumps
        cv.redun_vac = math.floor(vacv.n_vac_pumps_high * cv.redun_vacp / 100.0 + 0.5e0)

        u_unplanned_vacuum = self.calc_u_unplanned_vacuum(output)

        # Total unplanned unavailability
        u_unplanned = min(
            1.0e0,
            u_unplanned_magnets
            + u_unplanned_div
            + u_unplanned_fwbs
            + u_unplanned_bop
            + u_unplanned_hcd
            + u_unplanned_vacuum,
        )

        # Total availability
        cv.f_t_plant_available = max(
            1.0e0 - (u_planned + u_unplanned + u_planned * u_unplanned), 0.0e0
        )

        # Modify lifetimes to take account of the availability
        if ifev.ife != 1:
            # First wall / blanket
            if fwbsv.life_blkt_fpy < cv.life_plant:
                fwbsv.life_blkt_fpy = min(
                    fwbsv.life_blkt_fpy / cv.f_t_plant_available, cv.life_plant
                )
                # Current drive system lifetime (assumed equal to first wall and blanket lifetime)
                cv.life_hcd_fpy = fwbsv.life_blkt_fpy

            # Divertor
            if cv.life_div_fpy < cv.life_plant:
                cv.life_div_fpy = min(
                    cv.life_div_fpy / cv.f_t_plant_available, cv.life_plant
                )

            # Centrepost
            if pv.itart == 1 and cv.cplife < cv.life_plant:
                cv.cplife = min(cv.cplife / cv.f_t_plant_available, cv.life_plant)

        # Capacity factor
        cv.cpfact = cv.f_t_plant_available * (
            tv.t_plant_pulse_burn / tv.t_plant_pulse_total
        )

        # Output
        if output:
            po.ovarre(
                self.outfile,
                "First wall / blanket lifetime (FPY)",
                "(life_blkt_fpy)",
                fwbsv.life_blkt_fpy,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Divertor lifetime (FPY)",
                "(life_div_fpy)",
                cv.life_div_fpy,
                "OP ",
            )
            if pv.itart == 1:
                po.ovarre(
                    self.outfile,
                    "Centrepost lifetime (FPY)",
                    "(cplife)",
                    cv.cplife,
                    "OP ",
                )
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Total unavailability:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Total planned unavailability",
                "(u_planned)",
                u_planned,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total unplanned unavailability",
                "(u_unplanned)",
                u_unplanned,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Total plant availability fraction",
                "(f_t_plant_available)",
                cv.f_t_plant_available,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total DT operational time (years)",
                "(t_plant_operational_total_yrs)",
                cv.t_plant_operational_total_yrs,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total plant lifetime (years)",
                "(life_plant)",
                cv.life_plant,
            )
            po.ovarre(
                self.outfile,
                "Capacity factor: total lifetime elec. energy output / output power",
                "(cpfact)",
                cv.cpfact,
                "OP ",
            )

    def calc_u_planned(self, output: bool) -> float:
        """Calculates the planned unavailability of the plant


        This routine calculates the planned unavailability of the
        plant, using the methodology outlined in the 2014 EUROfusion
        RAMI report.
        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            planned unavailability of plant
        """

        # Full power lifetimes (in years) !

        # Caculate DPA per FPY
        # Detailed and cited in T. Franke 2020, "The EU DEMO equatorial outboard limiter — Design and port integration concept"
        # https://www.sciencedirect.com/science/article/pii/S0920379620301952#bib0075
        # Scaling w.r.t. fusion power drops out a large number of factors relating to neutronics, such as:
        # - the actual neutron flux
        # - the geometry and material composition leading to the neutron flux at the EUROfer FW OMP
        # - the neutron energy spectrum
        # - all of the above and more leading to the dpa/fpy in EUROfer at the FW OMP
        # About a relatively "constant" reference point, we can reasonably assume they all equal to 1.0.
        ref_fusion_power = 2.0e3  # (MW) fusion power for EU-DEMO
        f_scale = pv.p_fusion_total_mw / ref_fusion_power
        ref_dpa_fpy = 10.0e0  # dpa per fpy from T. Franke 2020 states up to 10 dpa/FPY
        dpa_fpy = f_scale * ref_dpa_fpy

        # First wall / blanket lifetime (years)
        # Calculate blanket lifetime using neutron fluence model (ibkt_life=0)
        # or DEMO fusion power model (ibkt_life=1)
        if cv.ibkt_life == 0:
            fwbsv.life_blkt_fpy = min(
                cv.abktflnc / pv.pflux_fw_neutron_mw, cv.life_plant
            )
        else:
            fwbsv.life_blkt_fpy = min(cv.life_dpa / dpa_fpy, cv.life_plant)  # DEMO

        # Divertor lifetime (years)
        cv.life_div_fpy = self.divertor_lifetime()

        # Centrepost lifetime (years) (ST only)
        if pv.itart == 1:
            cv.cplife = self.cp_lifetime()

        # Current drive lifetime (assumed equal to first wall and blanket lifetime)
        cv.life_hcd_fpy = fwbsv.life_blkt_fpy

        # Calculate the blanket and divertor replacement times !

        # Blanket replacement time
        # ( Calculated using scaling from 2014 EUROfusion RAMI report )

        # Mean time to repair blanket is same as replacing both blanket and divertor.
        # The +2.0 at the end is for the 1 month cooldown and pump down at either end
        # of the maintenance period
        mttr_blanket = (21.0e0 * cv.num_rh_systems ** (-0.9e0) + 2.0e0) / 12.0e0

        # Mean time to repair divertor is 70% of time taken to replace blanket
        # This is taken from Oliver Crofts 2014 paper
        mttr_divertor = 0.7e0 * mttr_blanket

        #  Which component has the shorter life?
        if cv.life_div_fpy < fwbsv.life_blkt_fpy:
            lifetime_shortest = cv.life_div_fpy
            lifetime_longest = fwbsv.life_blkt_fpy
            mttr_shortest = mttr_divertor
        else:
            lifetime_shortest = fwbsv.life_blkt_fpy
            lifetime_longest = cv.life_div_fpy
            mttr_shortest = mttr_blanket

        # Number of outages between each combined outage
        n = math.ceil(lifetime_longest / lifetime_shortest) - 1

        # Planned unavailability
        u_planned = (n * mttr_shortest + mttr_blanket) / (
            (n + 1) * lifetime_shortest + (n * mttr_shortest + mttr_blanket)
        )

        # Output
        if output:
            po.oheadr(self.outfile, "Plant Availability (2014 Model)")

            po.ocmmnt(self.outfile, "Planned unavailability:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Allowable blanket neutron fluence (MW-yr/m2)",
                "(abktflnc)",
                cv.abktflnc,
            )
            po.ovarre(
                self.outfile,
                "Allowable divertor heat fluence (MW-yr/m2)",
                "(adivflnc)",
                cv.adivflnc,
            )
            po.ovarin(
                self.outfile,
                "Number of remote handling systems",
                "(num_rh_systems)",
                cv.num_rh_systems,
            )
            po.ovarre(
                self.outfile,
                "Time needed to replace divertor (yrs)",
                "(mttr_divertor)",
                mttr_divertor,
            )
            po.ovarre(
                self.outfile,
                "Time needed to replace blanket (yrs)",
                "(mttr_blanket)",
                mttr_blanket,
            )
            po.ovarre(
                self.outfile,
                "Time needed to replace blkt + div (yrs)",
                "(mttr_blanket)",
                mttr_blanket,
            )
            po.ovarre(
                self.outfile,
                "Total planned unavailability",
                "(uplanned)",
                u_planned,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_planned

    def calc_u_unplanned_magnets(self, output: bool) -> float:
        """Calculates the unplanned unavailability of the magnets


        This routine calculates the unplanned unavailability of the magnets,
        using the methodology outlined in the 2014 EUROfusion
        RAMI report.
        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            unplanned unavailability of magnets
        """

        # Magnet temperature margin limit (K)
        # Use the lower of the two values.  Issue #526
        tmargmin = min(
            tfv.temp_tf_superconductor_margin_min, tfv.temp_cs_superconductor_margin_min
        )
        mag_temp_marg_limit = tmargmin

        # Magnet maintenance time (years)
        mag_main_time = 0.5e0

        # Minimum unplanned unavailability
        mag_min_u_unplanned = mag_main_time / (
            cv.t_plant_operational_total_yrs + mag_main_time
        )

        # Point at which risk of unplanned unavailability increases
        # conf_mag is the c factor, which determines the temperature margin at which
        # lifetime starts to decline.
        start_of_risk = mag_temp_marg_limit / cv.conf_mag

        # Determine if temperature margin is in region with risk of unplanned unavailability
        if tfv.temp_margin >= start_of_risk:
            u_unplanned_magnets = mag_min_u_unplanned
        else:
            # Linear decrease in expected lifetime when approaching the limit
            t_life = max(
                0.0e0,
                (cv.t_plant_operational_total_yrs / (start_of_risk - tmargmin))
                * (tfv.temp_margin - tmargmin),
            )
            u_unplanned_magnets = mag_main_time / (t_life + mag_main_time)

        # Output !
        # !!!!!!!!!

        if output:
            po.ocmmnt(self.outfile, "Magnets:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile, "Minimum temperature margin (K)", "(tmargmin)", tmargmin
            )
            po.ovarre(
                self.outfile,
                "c parameter, determining the temp margin where lifetime declines",
                "(conf_mag)",
                cv.conf_mag,
            )
            po.ovarre(
                self.outfile,
                "Temperature Margin (K)",
                "(temp_margin)",
                tfv.temp_margin,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Magnets unplanned unavailability",
                "(u_unplanned_magnets)",
                u_unplanned_magnets,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_unplanned_magnets

    def calc_u_unplanned_divertor(self, output: bool) -> float:
        """Calculates the unplanned unavailability of the divertor


        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            unplanned unavailability of the divertor
        """

        # Calculate cycle limit in terms of days
        # Number of cycles between planned blanket replacements, N
        n = cv.life_div_fpy * YEAR_SECONDS / tv.t_plant_pulse_total

        # The probability of failure in one pulse cycle (before the reference cycle life)
        pf = (cv.div_prob_fail / DAY_SECONDS) * tv.t_plant_pulse_total
        a0 = 1.0e0 - pf * cv.div_umain_time * YEAR_SECONDS / tv.t_plant_pulse_total

        # Integrating the instantaneous availability gives the mean
        # availability over the planned cycle life N
        if cv.div_nu <= cv.div_nref:
            logger.error(
                """div_nu <= div_nref
            The cycle when the divertor fails with 100% probability <= & Reference value for cycle life of divertor
            """
            )
            po.ocmmnt(
                self.outfile,
                "ERROR: The cycle when the divertor fails with 100% probability & <= Reference value for cycle cycle life of divertor",
            )

        # Check number of cycles

        # Less than reference (availability is min availability)
        if n <= cv.div_nref:
            div_avail = a0

        # Greater than cycle number with 100% failure rate
        elif n >= cv.div_nu:
            div_avail = 0.0e0

        # Else number of cycles is inbetween and is given by formula below
        else:
            div_avail = (a0 / (cv.div_nu - cv.div_nref)) * (
                cv.div_nu - 0.5e0 * cv.div_nref**2.0e0 / n - 0.5e0 * n
            )

        # Unplanned unavailability for divertor
        u_unplanned_div = 1.0e0 - div_avail

        # Output
        if output:
            po.ocmmnt(self.outfile, "Divertor:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Probability of failure per operational day",
                "(div_prob_fail)",
                cv.div_prob_fail,
            )
            po.ovarre(
                self.outfile,
                "Repair time (years)",
                "(div_umain_time)",
                cv.div_umain_time,
            )
            po.ovarre(
                self.outfile,
                "Reference value for cycle life",
                "(div_nref)",
                cv.div_nref,
            )
            po.ovarre(
                self.outfile,
                "The cycle when failure is 100% certain",
                "(div_nu)",
                cv.div_nu,
            )
            po.ovarre(
                self.outfile, "Number of cycles between planned replacements", "(n)", n
            )
            po.ovarre(
                self.outfile,
                "Unplanned unavailability",
                "(u_unplanned_div)",
                u_unplanned_div,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_unplanned_div

    def calc_u_unplanned_fwbs(self, output: bool) -> float:
        """Calculates the unplanned unavailability of the first wall and blanket


        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            unplanned unavailability of first wall and blanket
        """

        # Calculate cycle limit in terms of days

        # Number of cycles between planned blanket replacements, N
        n = fwbsv.life_blkt_fpy * YEAR_SECONDS / tv.t_plant_pulse_total

        # The probability of failure in one pulse cycle
        # (before the reference cycle life)
        pf = (cv.fwbs_prob_fail / DAY_SECONDS) * tv.t_plant_pulse_total
        a0 = 1.0e0 - pf * cv.fwbs_umain_time * YEAR_SECONDS / tv.t_plant_pulse_total

        if cv.fwbs_nu <= cv.fwbs_nref:
            logger.error(
                """fwbs_nu <= fwbs_nref
            The cycle when the blanket fails with 100% probability <= &Reference value for cycle life of blanket
            """
            )
            po.ocmmnt(
                self.outfile,
                "EROROR: The cycle when the blanket fails with 100% probability& <= Reference value for cycle life of blanket",
            )

        # Integrating the instantaneous availability gives the mean
        # availability over the planned cycle life N
        if n <= cv.fwbs_nref:
            fwbs_avail = a0
        elif n >= cv.fwbs_nu:
            fwbs_avail = 0.0e0
        else:
            fwbs_avail = (a0 / (cv.fwbs_nu - cv.fwbs_nref)) * (
                cv.fwbs_nu - 0.5e0 * cv.fwbs_nref**2.0e0 / n - 0.5e0 * n
            )

        # First wall / blanket unplanned unavailability
        u_unplanned_fwbs = 1.0e0 - fwbs_avail

        # Output
        if output:
            po.ocmmnt(self.outfile, "First wall / Blanket:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Probability of failure per operational day",
                "(fwbs_prob_fail)",
                cv.fwbs_prob_fail,
            )
            po.ovarre(
                self.outfile,
                "Repair time (years)",
                "(fwbs_umain_time)",
                cv.fwbs_umain_time,
            )
            po.ovarre(
                self.outfile,
                "Reference value for cycle life",
                "(fwbs_nref)",
                cv.fwbs_nref,
            )
            po.ovarre(
                self.outfile,
                "The cycle when failure is 100% certain",
                "(fwbs_nu)",
                cv.fwbs_nu,
            )
            po.ovarre(
                self.outfile, "Number of cycles between planned replacements", "(n)", n
            )
            po.ovarre(
                self.outfile,
                "Unplanned unavailability",
                "(u_unplanned_fwbs)",
                u_unplanned_fwbs,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_unplanned_fwbs

    def calc_u_unplanned_bop(self, output: bool) -> float:
        """Calculates the unplanned unavailability of the balance of plant


        This routine calculates the unplanned unavailability of the balance of plant,
        using the methodology outlined in the 2014 EUROfusion
        RAMI report.
        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            unplanned unavailability of balance of plant
        """

        # Balance of plant failure rate (failures per hour)
        # ENEA study WP13-DTM02-T01
        bop_fail_rate = 9.39e-5

        # Number of balance of plant failures in plant operational lifetime
        bop_num_failures = math.ceil(
            bop_fail_rate * DAYS_IN_YEAR * 24.0e0 * cv.t_plant_operational_total_yrs
        )

        # Balance of plant mean time to repair (years)
        # ENEA study WP13-DTM02-T01
        bop_mttr = 96.0e0 / (24.0e0 * DAYS_IN_YEAR)

        # Unplanned downtime balance of plant
        u_unplanned_bop = (bop_mttr * bop_num_failures) / (
            cv.t_plant_operational_total_yrs
        )

        # Output
        if output:
            po.ocmmnt(self.outfile, "Balance of plant:")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile, "Failure rate (1/h)", "(bop_fail_rate)", bop_fail_rate
            )
            po.ovarin(
                self.outfile,
                "Number of failures in lifetime",
                "(bop_num_failures)",
                bop_num_failures,
                "OP ",
            )
            po.ovarre(self.outfile, "Balance of plant MTTR", "(bop_mttr)", bop_mttr)
            po.ovarre(
                self.outfile,
                "Balance of plant unplanned unavailability",
                "(u_unplanned_bop)",
                u_unplanned_bop,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_unplanned_bop

    def calc_u_unplanned_hcd(self) -> float:
        """Calculates the unplanned unavailability of the heating and current drive system


        This routine calculates the unplanned unavailability of the heating
        and current drive system, using the methodology outlined in the
        2014 EUROfusion RAMI report.
        2014 EUROfusion RAMI report, &quot;Availability in PROCESS&quot;

        Returns
        -------
        float
            unplanned unavailability of hcd
        """

        # Currently just a fixed value until more information available or Q.
        # Tran's response provides useful data.

        return 0.02e0

    def calc_u_unplanned_vacuum(self, output: bool) -> float:
        """Calculates the unplanned unavailability of the vacuum system


        This routine calculates the unplanned unavailability of the vacuum system,
        using the methodology outlined in the 2014 EUROfusion
        RAMI report.
        2014 EUROfusion RAMI report, &quot;Availability in
        PROCESS&quot;

        Parameters
        ----------
        output : boolean
            indicate whether output should be written to the output file, or not
        output: bool :


        Returns
        -------
        float
            unplanned unavailability of vacuum system
        """

        # Number of shutdowns
        n_shutdown: int = round(
            (cv.life_plant - cv.t_plant_operational_total_yrs)
            / ((21.0e0 * cv.num_rh_systems ** (-0.9e0) + 2.0e0) / 12.0e0)
        )

        # Operational time between shutdowns
        t_op_bt = cv.t_plant_operational_total_yrs / (n_shutdown + 1.0e0)

        # Cryopump maintenance time (y) = 2 months
        cryo_main_time = 1.0e0 / 6.0e0

        # Total pumps = pumps + redundant pumps
        total_pumps = vacv.n_vac_pumps_high + cv.redun_vac

        # Cryopump failure rate per machine operational period
        # From "Selected component failure rate values from fusion
        # safety assessment tasks", Cadwallader (1994)

        # probability of pump failure per operational period
        cryo_failure_rate = 2.0e-6 * DAYS_IN_YEAR * 24.0e0 * t_op_bt

        # probability of no pump failure per operational period
        cryo_nfailure_rate = 1.0e0 - cryo_failure_rate

        sum_prob = 0.0e0

        for n in range(cv.redun_vac + 1, total_pumps + 1):
            # Probability for n failures in the operational period, n > number of redundant pumps

            # calculate sum in formula for downtime
            sum_prob = sum_prob + combinations(total_pumps, n) * (
                cryo_nfailure_rate ** (total_pumps - n)
            ) * (cryo_failure_rate**n) * (n - cv.redun_vac)

        # Total down-time in reactor life
        t_down = (n_shutdown + 1.0e0) * cryo_main_time * sum_prob

        # Total vacuum unplanned unavailability
        u_unplanned_vacuum = max(
            0.005, t_down / (cv.t_plant_operational_total_yrs + t_down)
        )

        # Output
        if output:
            po.ocmmnt(self.outfile, "Vacuum:")
            po.oblnkl(self.outfile)
            po.ovarin(
                self.outfile,
                "Number of pumps (excluding redundant pumps)",
                "(n_vac_pumps_high)",
                vacv.n_vac_pumps_high,
                "OP ",
            )
            po.ovarin(
                self.outfile,
                "Number of redundant pumps",
                "(redun_vac)",
                cv.redun_vac,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total unplanned down-time due to pumps, excl fixed 0.5% (years)",
                "(t_down)",
                t_down,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Vacuum unplanned unavailability",
                "(u_unplanned_vacuum)",
                u_unplanned_vacuum,
                "OP ",
            )
            po.oblnkl(self.outfile)

        return u_unplanned_vacuum

    def avail_st(self, output: bool):
        """Calculate the availability for a plant with a Spherical Tokamak.

        This routine calculates the availability of a plant by considering various factors such as
        the lifetime of different components, planned and unplanned unavailability, and maintenance cycles.

        Parameters:
        -----------
        output : bool
            Indicates whether the output should be written to the output file or not.

        Detailed Description:
        ---------------------
        - The method calculates the Displacements Per Atom (DPA) per Full Power Year (FPY) based on the fusion power.
        - It determines the lifetime of the first wall and blanket, divertor, and current drive.
        - It calculates the time for a maintenance cycle and the number of maintenance cycles over the plant's lifetime.
        - It computes the planned and unplanned unavailability of various components such as magnets, divertor, first wall and blanket, balance of plant, heating and current drive, and vacuum systems.
        - The total availability of the plant is then calculated considering both planned and unplanned unavailability.
        - The method also adjusts the lifetimes of components based on the calculated availability.
        - Finally, it calculates the capacity factor and operational time of the plant.

        If `output` is True, the method writes detailed availability information to the output file.

        References:
        -----------
        - T. Franke 2020, "The EU DEMO equatorial outboard limiter — Design and port integration concept"
          https://www.sciencedirect.com/science/article/pii/S0920379620301952#bib0075

        Notes:
        ------
        - The method assumes certain constants and reference points for calculations.
        - The method modifies the lifetimes of components to account for the calculated availability.

        Parameters
        ----------
        output: bool :
        """

        ref_powfmw = 2.0e3  # (MW) fusion power for EU-DEMO
        f_scale = pv.p_fusion_total_mw / ref_powfmw
        ref_dpa_fpy = 10.0e0  # dpa per fpy from T. Franke 2020 states up to 10 dpa/FPY
        dpa_fpy = f_scale * ref_dpa_fpy

        if cv.ibkt_life == 0:
            fwbsv.life_blkt_fpy = min(
                cv.abktflnc / pv.pflux_fw_neutron_mw, cv.life_plant
            )
        else:
            fwbsv.life_blkt_fpy = min(cv.life_dpa / dpa_fpy, cv.life_plant)  # DEMO

        # Divertor lifetime (years)
        cv.life_div_fpy = self.divertor_lifetime()

        # CP lifetime (years)
        cv.cplife = self.cp_lifetime()

        # Current drive lifetime (assumed equal to first wall and blanket lifetime)
        cv.life_hcd_fpy = fwbsv.life_blkt_fpy

        # Time for a maintenance cycle (years)
        # Shortest component lifetime + time to replace
        shortest_lifetime = min(
            fwbsv.life_blkt_fpy,
            cv.life_div_fpy,
            cv.cplife,
            cv.life_hcd_fpy,
            cv.life_plant,
        )
        maint_cycle = shortest_lifetime + cv.tmain

        # Number of maintenance cycles over plant lifetime
        n_cycles_main = cv.life_plant / maint_cycle

        # Number of centre columns over plant lifetime
        n_centre_cols = math.ceil(n_cycles_main)

        # Planned unavailability
        u_planned = cv.tmain / maint_cycle

        # Operational time (years)
        cv.t_plant_operational_total_yrs = cv.life_plant * (1.0e0 - u_planned)

        if output:
            po.oheadr(self.outfile, "Plant Availability")

        # Un-planned unavailability

        # Magnets
        u_unplanned_magnets = self.calc_u_unplanned_magnets(output)

        # Divertor
        u_unplanned_div = self.calc_u_unplanned_divertor(output)

        # First wall and blanket
        u_unplanned_fwbs = self.calc_u_unplanned_fwbs(output)

        # Balance of plant
        u_unplanned_bop = self.calc_u_unplanned_bop(output)

        # Heating and current drive
        u_unplanned_hcd = self.calc_u_unplanned_hcd()

        # Vacuum systems

        # Number of redundant pumps
        cv.redun_vac = math.floor(vacv.n_vac_pumps_high * cv.redun_vacp / 100.0 + 0.5e0)

        u_unplanned_vacuum = self.calc_u_unplanned_vacuum(output)

        # Total unplanned unavailability
        u_unplanned = min(
            1.0e0,
            u_unplanned_magnets
            + u_unplanned_div
            + u_unplanned_fwbs
            + u_unplanned_bop
            + u_unplanned_hcd
            + u_unplanned_vacuum
            + cv.u_unplanned_cp,
        )

        # Total availability
        cv.f_t_plant_available = max(
            1.0e0 - (u_planned + u_unplanned + u_planned * u_unplanned), 0.0e0
        )

        # Modify lifetimes to take account of the availability
        if ifev.ife != 1:
            # First wall / blanket
            if fwbsv.life_blkt_fpy < cv.life_plant:
                fwbsv.life_blkt_fpy = min(
                    fwbsv.life_blkt_fpy / cv.f_t_plant_available, cv.life_plant
                )
                cv.life_hcd_fpy = fwbsv.life_blkt_fpy

            # Divertor
            if cv.life_div_fpy < cv.life_plant:
                cv.life_div_fpy = min(
                    cv.life_div_fpy / cv.f_t_plant_available, cv.life_plant
                )

            # Centrepost
            if pv.itart == 1 and cv.cplife < cv.life_plant:
                cv.cplife = min(cv.cplife / cv.f_t_plant_available, cv.life_plant)

        # Capacity factor
        cv.cpfact = cv.f_t_plant_available * (
            tv.t_plant_pulse_burn / tv.t_plant_pulse_total
        )

        if output:
            po.ocmmnt(self.outfile, "Plant Availability")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Allowable blanket neutron fluence (MW-yr/m2)",
                "(abktflnc)",
                cv.abktflnc,
            )
            po.ovarre(
                self.outfile,
                "Allowable divertor heat fluence (MW-yr/m2)",
                "(adivflnc)",
                cv.adivflnc,
            )
            po.ovarre(
                self.outfile,
                "First wall / blanket lifetime (FPY)",
                "(life_blkt_fpy)",
                fwbsv.life_blkt_fpy,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Divertor lifetime (FPY)",
                "(life_div_fpy)",
                cv.life_div_fpy,
                "OP ",
            )
            if tfv.i_tf_sup == 1:
                po.ovarre(
                    self.outfile,
                    "Max fast neutron fluence on TF coil (n/m2)",
                    "(nflutfmax)",
                    ctv.nflutfmax,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Centrepost TF fast neutron flux (E > 0.1 MeV) (m^(-2).^(-1))",
                    "(neut_flux_cp)",
                    fwbsv.neut_flux_cp,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Allowable ST centrepost neutron fluence (MW-yr/m2)",
                    "(cpstflnc)",
                    cv.cpstflnc,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Average neutron wall load (MW/m2)",
                    "(pflux_fw_neutron_mw)",
                    pv.pflux_fw_neutron_mw,
                    "OP ",
                )
            po.ovarre(
                self.outfile,
                "Centrepost lifetime (years)",
                "(cplife)",
                cv.cplife,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Maintenance time for replacing CP (years)",
                "(tmain)",
                cv.tmain,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Length of maintenance cycle (years)",
                "(maint_cycle)",
                maint_cycle,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Number of maintenance cycles over lifetime",
                "(n_cycles_main)",
                n_cycles_main,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Number of centre columns over lifetime",
                "(n_centre_cols)",
                n_centre_cols,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Total planned unavailability",
                "(u_planned)",
                u_planned,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total unplanned unavailability",
                "(u_unplanned)",
                u_unplanned,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total plant availability fraction",
                "(f_t_plant_available)",
                cv.f_t_plant_available,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Capacity factor: total lifetime elec. energy output / output power",
                "(cpfact)",
                cv.cpfact,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total DT operational time (years)",
                "(t_plant_operational_total_yrs)",
                cv.t_plant_operational_total_yrs,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total plant lifetime (years)",
                "(life_plant)",
                cv.life_plant,
                "OP",
            )

    @staticmethod
    def cp_lifetime():
        """Calculate Centrepost Lifetime

        This routine calculates the lifetime of the centrepost,
        either for superconducting or aluminium/resistive magnets.

        Returns
        -------
        float
            CP lifetime
        """
        # SC magnets CP lifetime
        # Rem : only the TF maximum fluence is considered for now
        if tfv.i_tf_sup == 1:
            cplife = (
                cv.life_plant
                if fwbsv.neut_flux_cp <= 0.0
                else min(
                    (ctv.nflutfmax / (fwbsv.neut_flux_cp * YEAR_SECONDS)),
                    cv.life_plant,
                )
            )

        # Aluminium/Copper magnets CP lifetime
        # For now, we keep the original def, developed for GLIDCOP magnets ...
        else:
            cplife = min(cv.cpstflnc / pv.pflux_fw_neutron_mw, cv.life_plant)

        return cplife

    @staticmethod
    def divertor_lifetime():
        """Calculate Divertor Lifetime

        This routine calculates the lifetime of the divertor based on the allowable divertor heat fluence.

        Returns
        -------
        float
            Divertor lifetime
        """
        # Divertor lifetime
        # Either 0.0, calculated from allowable divertor fluence and heat load, or lifetime of the plant
        return max(0.0, min(cv.adivflnc / dv.pflux_div_heat_load_mw, cv.life_plant))
