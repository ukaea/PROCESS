import copy

import numpy as np

from process import fortran as ft
from process import process_output as po
from process.build import Build
from process.fortran import build_variables as bv
from process.fortran import constants
from process.fortran import error_handling as eh
from process.fortran import fwbs_variables as fwbsv
from process.fortran import tfcoil_variables as tfv
from process.sctfcoil import Sctfcoil


class TFcoil:
    """Calculates the parameters of a resistive TF coil system for a fusion power plant"""

    def __init__(self, build: Build, sctfcoil: Sctfcoil):
        """Initialise Fortran module variables."""
        self.outfile = ft.constants.nout  # output file unit
        self.iprint = 0  # switch for writing to output file (1=yes)
        self.build = build
        self.sctfcoil = sctfcoil

    def run(self):
        """Run main tfcoil subroutine without outputting."""
        self.iprint = 0
        self.tfcoil()

    def output(self):
        """Run main tfcoil subroutine and write output."""
        self.iprint = 1
        self.tfcoil()

    def tfcoil(self):
        """TF coil module
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This subroutine calculates various parameters for the TF coil set.
        If the TF coils are superconducting the calculations are performed
        in routine <A HREF="sctfcoil.html">sctfcoil</A> instead.
        """

        # TF coil calculations
        self.sctfcoil.sctfcoil(output=bool(self.iprint))

        # Port size calculation
        self.build.portsz()

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

        # Temperature margin used in calculations (K)
        tmarg = 10.0e0

        # Number of integral step used for the coolant temperature rise
        n_tcool_it = 20

        # Coolant channels:
        acool = tfv.a_cp_cool * tfv.n_tf_coils  # Cooling cross-sectional area
        dcool = 2.0e0 * tfv.rcool  # Diameter
        lcool = 2.0e0 * (bv.hmax + bv.dr_tf_outboard)  # Length
        tfv.ncool = acool / (constants.pi * tfv.rcool**2)  # Number

        # Average conductor cross-sectional area to cool (with cooling area)
        acpav = 0.5e0 * tfv.vol_cond_cp / (bv.hmax + bv.dr_tf_outboard) + acool
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
        tfv.ppump = dpres * acool * tfv.vcool / tfv.etapump

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

            po.ovarre(self.outfile, "Pump power (W)", "(ppump)", tfv.ppump)

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
