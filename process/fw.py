import numpy as np

from process import process_output as po
from process.coolprop_interface import FluidProperties
from process.data_structure import build_variables
from process.fortran import (
    blanket_library,
    constants,
    error_handling,
    fwbs_variables,
)
from process.fortran import (
    error_handling as eh,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string


class Fw:
    def __init__(self) -> None:
        self.outfile = constants.nout

    def run(self):
        (
            blanket_library.n_fw_inboard_channels,
            blanket_library.n_fw_outboard_channels,
        ) = self.calculate_total_fw_channels(
            build_variables.a_fw_inboard,
            build_variables.a_fw_outboard,
            fwbs_variables.len_fw_channel,
            fwbs_variables.dx_fw_module,
        )

        self.set_fw_geometry()

    def set_fw_geometry(self):
        build_variables.dr_fw_inboard = (
            2 * fwbs_variables.radius_fw_channel + 2 * fwbs_variables.dr_fw_wall
        )
        build_variables.dr_fw_outboard = build_variables.dr_fw_inboard

    def fw_temp(
        self,
        output: bool,
        radius_fw_channel: float,
        dr_fw: float,
        a_fw: float,
        prad_incident: float,
        pnuc_deposited: float,
        label: str,
    ) -> tuple:
        """
        Thermo-hydraulic calculations for the first wall.

        :param output: Flag to indicate if output is required.
        :type output: bool
        :param radius_fw_channel: First wall coolant channel radius (m).
        :type radius_fw_channel: float
        :param dr_fw: First wall thickness (m).
        :type dr_fw: float
        :param a_fw: Area of first wall section under consideration (m^2).
        :type a_fw: float
        :param prad_incident: Radiation surface heat flux on first wall (MW).
        :type prad_incident: float
        :param pnuc_deposited: Nuclear power deposited in FW (MW).
        :type pnuc_deposited: float
        :param label: Information string.
        :type label: str

        :returns: Contains peak first wall temperature (K), coolant specific heat capacity at constant pressure (J/kg/K),
        :rtype: tuple

        Detailed thermal hydraulic model for the blanket (first wall + breeding zone).
        Given the heating incident on the first wall, and the coolant outlet temperature,
        the maximum temperature of the first wall is calculated to check it is below material limits.
        The routine is called separately for the inboard and outboard sides.
        The calculation of the maximum temperature is described by Gardner:
        "Temperature distribution in the first wall", K:\\Power Plant Physics and Technology\\ PROCESS\\PROCESS References & Systems Codes\\Pulsed option - Gardner.
        This is in turn taken from "Methods of First Wall Structural Analysis with Application to the Long Pulse Commercial Tokamak Reactor Design", R.J. LeClaire, MIT, PFC/RR-84-9.
        """

        # First wall volume (inboard or outboard depending on arguments) (m^3)
        vol_fw = a_fw * dr_fw

        # First wall channel area (m^2)
        a_fw_channel = np.pi * radius_fw_channel**2

        # Heat generation in the first wall due to neutron flux deposited in the material (W/m3)
        pden_fw_nuclear = 1e6 * pnuc_deposited / vol_fw

        # the nuclear heating in the coolant is small. (W/m2)
        # Note that the full first wall volume is used including coolant even though
        nuclear_heat_per_area = pden_fw_nuclear * dr_fw

        # Heat flux incident on the first wall surface (W/m2)
        pflux_fw_rad = 1e6 * prad_incident / a_fw

        # Calculate inlet coolant fluid properties (fixed pressure)
        inlet_coolant_properties = FluidProperties.of(
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type),
            temperature=fwbs_variables.temp_fw_coolant_in.item(),
            pressure=fwbs_variables.pres_fw_coolant.item(),
        )

        # Calculate outlet coolant fluid properties (fixed pressure)
        outlet_coolant_properties = FluidProperties.of(
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type),
            temperature=fwbs_variables.temp_fw_coolant_out.item(),
            pressure=fwbs_variables.pres_fw_coolant.item(),
        )

        # Mean properties (inlet + outlet)/2
        # Average coolant density (kg/m3)
        den_fw_coolant_average = (
            inlet_coolant_properties.density + outlet_coolant_properties.density
        ) / 2

        # Mean properties (inlet + outlet)/2
        # Average coolant specific heat capacity (J/K)
        heatcap_fw_coolant_average = (
            inlet_coolant_properties.specific_heat_const_p
            + outlet_coolant_properties.specific_heat_const_p
        ) / 2

        # Heat load per unit length of one first wall segment (W/m)
        # Nuclear particle and radiation heating
        load = (nuclear_heat_per_area + pflux_fw_rad) * fwbs_variables.dx_fw_module

        # Coolant mass flow rate (kg/s) (use mean properties)
        mflow_fw_coolant = (
            fwbs_variables.len_fw_channel
            * load
            / heatcap_fw_coolant_average
            / (fwbs_variables.temp_fw_coolant_out - fwbs_variables.temp_fw_coolant_in)
        )

        # Coolant mass flux in a single channel (kg/m2/s)
        mflux_fw_coolant = mflow_fw_coolant / a_fw_channel

        # Conditions at the outlet, where the temperature is highest
        # -----------------------------------------------------------

        # Coolant velocity (m/s)
        vel_fw_coolant_average = mflux_fw_coolant / outlet_coolant_properties.density

        # Mean temperature of the wall material on the plasma side of the coolant 'temp_fw_peak'
        # is the estimate from the previous iteration of the wall surface temperature
        # (underneath the armour)
        temp_k = (fwbs_variables.temp_fw_coolant_out + fwbs_variables.temp_fw_peak) / 2

        # Print debug info if temperature too low/high or NaN/Inf
        if np.isnan(temp_k):
            eh.report_error(223)
        elif (temp_k <= 100) or (temp_k > 1500):
            eh.fdiags[0] = temp_k
            eh.report_error(224)

        # Thermal conductivity of first wall material (W/m.K)
        tkfw = self.fw_thermal_conductivity(temp_k)

        # Heat transfer coefficient (W m^-2 K^-1)
        hcoeff = self.heat_transfer(
            mflux_fw_coolant,
            outlet_coolant_properties.density,
            radius_fw_channel,
            outlet_coolant_properties.specific_heat_const_p,
            outlet_coolant_properties.viscosity,
            outlet_coolant_properties.thermal_conductivity,
            fwbs_variables.roughness_fw_channel,
        )

        # Temperature drops between first-wall surface and bulk coolant !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Model B is given for comparison
        # Model C is used
        # Model A LeClaire formula for circular pipes (removed GitHub #389).

        # Model B: Simple 1-dimensional calculation !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # This is optimistic as it neglects the higher temperature midway between the channels.
        # I have included 1/4 of the volume load:
        # 1/2 is absorbed in the plasma-facing wall (A)
        # which on average has to pass through 1/2 the wall thickness.
        #  ______________
        #        A
        #  --------------
        #     'channel'
        #  --------------
        #  ______________

        # Worst case load (as above) per unit length in 1-D calculation (W/m)
        onedload = fwbs_variables.f_fw_peak * (
            pden_fw_nuclear * fwbs_variables.dx_fw_module * dr_fw / 4
            + pflux_fw_rad * fwbs_variables.dx_fw_module
        )

        # Effective area for heat transfer (m2)
        effective_area_for_heat_transfer = 2 * radius_fw_channel

        # Temperature drop in first-wall material (K)
        deltat_solid_1D = onedload * dr_fw / (tkfw * effective_area_for_heat_transfer)

        # Model C: A more realistic model !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Calculate maximum distance travelled by surface heat load (m)
        # dr_fw_wall | Minimum distance travelled by surface heat load (m)
        diagonal = np.sqrt(
            (fwbs_variables.dx_fw_module / 2 - radius_fw_channel) ** 2
            + (radius_fw_channel + dr_fw) ** 2
        )

        # Mean distance travelled by surface heat (m)
        mean_distance = (fwbs_variables.dr_fw_wall + diagonal) / 2

        # This heat starts off spread over width = 'dx_fw_module'.
        # It ends up spread over one half the circumference.
        # Use the mean of these values.
        mean_width = (
            fwbs_variables.dx_fw_module + np.pi * radius_fw_channel
        ) / 2  # (m)

        # As before, use a combined load 'onedload'
        # Temperature drop in first-wall material (K)
        deltat_solid = onedload * mean_distance / (tkfw * mean_width)

        # Temperature drop between channel inner wall and bulk coolant (K)
        deltat_coolant = load / (2 * np.pi * radius_fw_channel * hcoeff)

        # Peak first wall temperature (K)
        tpeakfw = fwbs_variables.temp_fw_coolant_out + deltat_solid + deltat_coolant

        if output:
            po.oheadr(
                self.outfile, "Heat transfer parameters at the coolant outlet: " + label
            )
            po.ovarre(
                self.outfile,
                "Radius of FW coolant channel (m)",
                "(radius_fw_channel)",
                radius_fw_channel,
            )
            po.ovarre(
                self.outfile,
                "Mean surface radiation flux on first wall (W/m2) ",
                "(pflux_fw_rad)",
                pflux_fw_rad,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mean nuclear power deposited in first wall per unit area (W/m2)",
                "",
                nuclear_heat_per_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Ratio of peak local heat load (surface and nuclear) to mean",
                "(f_fw_peak)",
                fwbs_variables.f_fw_peak,
            )
            po.ovarre(
                self.outfile,
                "Vertical length of a single coolant channel (all in parallel) (m)",
                "(len_fw_channel)",
                fwbs_variables.len_fw_channel,
            )
            po.ovarre(
                self.outfile,
                "Width of a FW module containing a cooling channel [m]",
                "(dx_fw_module)",
                fwbs_variables.dx_fw_module,
            )
            po.ovarre(
                self.outfile,
                "Thermal conductivity of first wall material (W/K/m)",
                "(tkfw)",
                tkfw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Coolant density (kg/m3)",
                "(rhofo)",
                outlet_coolant_properties.density,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Coolant mass flow rate in one channel (kg/s)",
                "(mflow_fw_coolant)",
                mflow_fw_coolant,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Coolant velocity (m/s)",
                "(vel_fw_coolant_average)",
                vel_fw_coolant_average,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Outlet temperature of first wall coolant (K)",
                "(temp_fw_coolant_out)",
                fwbs_variables.temp_fw_coolant_out,
            )
            po.ovarre(
                self.outfile, "Heat transfer coefficient", "(hcoeff)", hcoeff, "OP "
            )
            po.ovarre(
                self.outfile,
                "Temperature drop in the wall material (simple model)",
                "(deltat_solid)",
                deltat_solid,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Temperature drop in the coolant (wall to bulk)",
                "(deltat_coolant)",
                deltat_coolant,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "First wall temperature (excluding armour) (K)",
                "(tpeakfw)",
                tpeakfw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Temperature drop in the wall material: 1D estimate",
                "(deltat_solid_1D)",
                deltat_solid_1D,
                "OP ",
            )

        return (
            tpeakfw,
            heatcap_fw_coolant_average,
            den_fw_coolant_average,
            mflow_fw_coolant,
        )

    def fw_thermal_conductivity(self, temp: float) -> float:
        """
        Calculates the thermal conductivity of the first wall material (Eurofer97).

        :param temp: Property temperature in Kelvin (K).
        :type temp: float
        :return: Thermal conductivity of Eurofer97 in W/m/K.
        :rtype: float

        :notes:
            Valid up to about 800 K

        :references:
            - A. A. Tavassoli et al., “Materials design data for reduced activation martensitic steel type EUROFER,”
            Journal of Nuclear Materials, vol. 329-333, pp. 257-262, Aug. 2004,
            doi: https://doi.org/10.1016/j.jnucmat.2004.04.020.

            - Tavassoli, F. "Fusion Demo Interim Structural Design Criteria (DISDC)/Appendix A Material Design Limit Data/A3. S18E Eurofer Steel."
              CEA, EFDA_TASK_TW4-TTMS-005-D01 (2004)
        """

        # temp in Kelvin
        return (
            (5.4308 + 0.13565 * temp - 0.00023862 * temp**2 + 1.3393e-7 * temp**3)
            * fwbs_variables.fw_th_conductivity
            / 28.34
        )

    def heat_transfer(
        self,
        mflux_coolant: float,
        den_coolant: float,
        radius_channel: float,
        heatcap_coolant: float,
        visc_coolant: float,
        thermcond_coolant: float,
        roughness_fw_channel: float,
    ) -> float:
        """
        Calculate heat transfer coefficient using Gnielinski correlation.

        :param mflux_coolant: Coolant mass flux in a single channel (kg/m^2/s).
        :type mflux_coolant: float
        :param den_coolant: Coolant density (average of inlet and outlet) (kg/m^3).
        :type den_coolant: float
        :param radius_channel: Coolant pipe radius (m).
        :type radius_channel: float
        :param heatcap_coolant: Coolant specific heat capacity (average of inlet and outlet) (J/kg/K).
        :type heatcap_coolant: float
        :param visc_coolant: Coolant viscosity (average of inlet and outlet) (Pa.s).
        :type visc_coolant: float
        :param thermcond_coolant: Thermal conductivity of coolant (average of inlet and outlet) (W/m.K).
        :type thermcond_coolant: float
        :param roughness_fw_channel: Roughness of the first wall coolant channel (m).
        :type roughness_fw_channel: float
        :return: Heat transfer coefficient (W/m^2K).
        :rtype: float

        :notes:
            Gnielinski correlation. Ignore the distinction between wall and
            bulk temperatures. Valid for: 3000 < Re < 5e6, 0.5 < Pr < 2000

        :references:
            - https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation
        """
        # Calculate pipe diameter (m)
        diameter = 2 * radius_channel

        # Calculate flow velocity (m/s)
        velocity = mflux_coolant / den_coolant

        # Calculate Reynolds number
        reynolds = den_coolant * velocity * diameter / visc_coolant

        # Calculate Prandtl number
        pr = heatcap_coolant * visc_coolant / thermcond_coolant

        # Calculate Darcy friction factor, using Haaland equation
        f = self.darcy_friction_haaland(
            reynolds,
            roughness_fw_channel,
            radius_channel,
        )

        # Calculate the Nusselt number
        nusselt = (
            (f / 8.0)
            * (reynolds - 1000.0)
            * pr
            / (1 + 12.7 * np.sqrt(f / 8.0) * (pr ** (2 / 3) - 1.0))
        )

        # Calculate the heat transfer coefficient (W/m^2K)
        heat_transfer_coefficient = nusselt * thermcond_coolant / (2.0 * radius_channel)

        # Check that Reynolds number is in valid range for the Gnielinski correlation
        if (reynolds <= 3000.0) or (reynolds > 5.0e6):
            error_handling.fdiags[0] = reynolds
            error_handling.report_error(225)

        # Check that Prandtl number is in valid range for the Gnielinski correlation
        if (pr < 0.5) or (pr > 2000.0):
            error_handling.fdiags[0] = pr
            error_handling.report_error(226)

        # Check that the Darcy friction factor is in valid range for the Gnielinski correlation
        if f <= 0.0:
            error_handling.report_error(227)

        return heat_transfer_coefficient

    def darcy_friction_haaland(
        self, reynolds: float, roughness_fw_channel: float, radius_fw_channel: float
    ) -> float:
        """
        Calculate Darcy friction factor using the Haaland equation.

        :param reynolds: Reynolds number.
            :type reynolds: float
            :param roughness_fw_channel: Roughness of the first wall coolant channel (m).
            :type roughness_fw_channel: float
            :param radius_fw_channel: Radius of the first wall coolant channel (m).
            :type radius_fw_channel: float

            :return: Darcy friction factor.
            :rtype: float

        :Notes:
            The Haaland equation is an approximation to the implicit Colebrook-White equation.
            It is used to calculate the Darcy friction factor for turbulent flow in pipes.

        :References:
            - https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
        """

        # Bracketed term in Haaland equation
        bracket = (
            roughness_fw_channel / radius_fw_channel / 3.7
        ) ** 1.11 + 6.9 / reynolds

        # Calculate Darcy friction factor
        return (1.8 * np.log10(bracket)) ** (-2)

    @staticmethod
    def calculate_total_fw_channels(
        a_fw_inboard: float,
        a_fw_outboard: float,
        len_fw_channel: float,
        dx_fw_module: float,
    ) -> tuple[int, int]:
        """
        Calculate the total number of first wall channels for inboard and outboard sections.

        Args:
            a_fw_inboard (float): Area of the inboard first wall section (m^2).
            a_fw_outboard (float): Area of the outboard first wall section (m^2).
            len_fw_channel (float): Length of each first wall channel (m).
            dx_fw_module (float): Toroidal width of each first wall module (m).

        Returns:
            tuple: Number of inboard and outboard first wall channels.
        """
        n_fw_inboard_channels = a_fw_inboard / (len_fw_channel * dx_fw_module)
        n_fw_outboard_channels = a_fw_outboard / (len_fw_channel * dx_fw_module)
        return int(n_fw_inboard_channels), int(n_fw_outboard_channels)

    def output_fw_geometry(self):
        """
        Outputs the first wall geometry details to the output file.

        Returns:
            None
        """
        po.oheadr(self.outfile, "First wall build")

        po.ovarrf(
            self.outfile,
            "Radius of first wall cooling channels (m)",
            "(radius_fw_channel)",
            fwbs_variables.radius_fw_channel,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Radial wall thickness surrounding first wall coolant channel (m)",
            "(dr_fw_wall)",
            fwbs_variables.dr_fw_wall,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Toroidal width of each first wall module (m)",
            "(dx_fw_module)",
            fwbs_variables.dx_fw_module,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Length of each first wall channel (m)",
            "(len_fw_channel)",
            fwbs_variables.len_fw_channel,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Radial thickness off inboard first wall (m)",
            "(dr_fw_inboard)",
            build_variables.dr_fw_inboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Radial thickness off outboard first wall (m)",
            "(dr_fw_outboard)",
            build_variables.dr_fw_outboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Number of inboard first wall cooling channels",
            "(n_fw_inboard_channels)",
            blanket_library.n_fw_inboard_channels,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Number of outboard first wall cooling channels",
            "(n_fw_outboard_channels)",
            blanket_library.n_fw_outboard_channels,
            "OP ",
        )

    def output_fw_pumping(self):
        """
        Outputs the first wall pumping details to the output file.

        Returns:
            None
        """
        po.oheadr(self.outfile, "First wall pumping")

        po.ovarst(
            self.outfile,
            "First wall coolant type",
            "(i_fw_coolant_type)",
            f"'{f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type)}'",
        )
        po.ovarrf(
            self.outfile,
            "Outlet temperature of first wall coolant [K]",
            "(temp_fw_coolant_out)",
            fwbs_variables.temp_fw_coolant_out,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Inlet temperature of first wall coolant [K]",
            "(temp_fw_coolant_in)",
            fwbs_variables.temp_fw_coolant_in,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pressure of first wall coolant [Pa]",
            "(pres_fw_coolant)",
            fwbs_variables.pres_fw_coolant,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Peak temperature of first wall [K]",
            "(temp_fw_peak)",
            fwbs_variables.temp_fw_peak,
            "OP ",
        )


def init_fwbs_variables():
    """Initialise FWBS variables"""
    fwbs_variables.life_blkt_fpy = 0.0
    fwbs_variables.life_blkt = 0.0
    fwbs_variables.m_fw_blkt_div_coolant_total = 0.0
    fwbs_variables.m_vv = 0.0
    fwbs_variables.denstl = 7800.0
    fwbs_variables.denwc = 15630.0
    fwbs_variables.dewmkg = 0.0
    fwbs_variables.f_p_blkt_multiplication = 1.269
    fwbs_variables.p_blkt_multiplication_mw = 0.0
    fwbs_variables.fblss = 0.09705
    fwbs_variables.f_ster_div_single = 0.115
    fwbs_variables.f_a_fw_outboard_hcd = 0.0
    fwbs_variables.fhole = 0.0
    fwbs_variables.i_fw_blkt_vv_shape = 2
    fwbs_variables.life_fw_fpy = 0.0
    fwbs_variables.m_fw_total = 0.0
    fwbs_variables.fw_armour_mass = 0.0
    fwbs_variables.fw_armour_thickness = 0.005
    fwbs_variables.fw_armour_vol = 0.0
    fwbs_variables.i_blanket_type = 1
    fwbs_variables.i_blkt_inboard = 1
    fwbs_variables.inuclear = 0
    fwbs_variables.qnuc = 0.0
    fwbs_variables.f_blkt_li6_enrichment = 30.0
    fwbs_variables.p_blkt_nuclear_heat_total_mw = 0.0
    fwbs_variables.p_div_nuclear_heat_total_mw = 0.0
    fwbs_variables.p_fw_nuclear_heat_total_mw = 0.0
    fwbs_variables.p_fw_hcd_nuclear_heat_mw = 0.0
    fwbs_variables.pnucloss = 0.0
    fwbs_variables.pnucvvplus = 0.0
    fwbs_variables.p_shld_nuclear_heat_mw = 0.0
    fwbs_variables.m_blkt_total = 0.0
    fwbs_variables.m_blkt_steel_total = 0.0
    fwbs_variables.armour_fw_bl_mass = 0.0
    fwbs_variables.breeder_f = 0.5
    fwbs_variables.breeder_multiplier = 0.75
    fwbs_variables.vfcblkt = 0.05295
    fwbs_variables.vfpblkt = 0.1
    fwbs_variables.m_blkt_li4sio4 = 0.0
    fwbs_variables.m_blkt_tibe12 = 0.0
    fwbs_variables.f_neut_shield = -1.0
    fwbs_variables.f_a_fw_coolant_inboard = 0.0
    fwbs_variables.f_a_fw_coolant_outboard = 0.0
    fwbs_variables.psurffwi = 0.0
    fwbs_variables.psurffwo = 0.0
    fwbs_variables.vol_fw_total = 0.0
    fwbs_variables.f_vol_blkt_steel = 0.0
    fwbs_variables.f_vol_blkt_li4sio4 = 0.0
    fwbs_variables.f_vol_blkt_tibe12 = 0.0
    fwbs_variables.breedmat = 1
    fwbs_variables.densbreed = 0.0
    fwbs_variables.fblbe = 0.6
    fwbs_variables.fblbreed = 0.154
    fwbs_variables.fblhebmi = 0.4
    fwbs_variables.fblhebmo = 0.4
    fwbs_variables.fblhebpi = 0.6595
    fwbs_variables.fblhebpo = 0.6713
    fwbs_variables.hcdportsize = 1
    fwbs_variables.nflutf = 0.0
    fwbs_variables.npdiv = 2
    fwbs_variables.nphcdin = 2
    fwbs_variables.nphcdout = 2
    fwbs_variables.tbr = 0.0
    fwbs_variables.tritprate = 0.0
    fwbs_variables.wallpf = 1.21
    fwbs_variables.whtblbreed = 0.0
    fwbs_variables.m_blkt_beryllium = 0.0
    fwbs_variables.i_coolant_pumping = 2
    fwbs_variables.i_shield_mat = 0
    fwbs_variables.i_thermal_electric_conversion = 0
    fwbs_variables.secondary_cycle_liq = 4
    fwbs_variables.i_blkt_coolant_type = 1
    fwbs_variables.i_fw_coolant_type = "helium"
    fwbs_variables.dr_fw_wall = 0.003
    fwbs_variables.radius_fw_channel = 0.006
    fwbs_variables.dx_fw_module = 0.02
    fwbs_variables.temp_fw_coolant_in = 573.0
    fwbs_variables.temp_fw_coolant_out = 823.0
    fwbs_variables.pres_fw_coolant = 15.5e6
    fwbs_variables.temp_fw_peak = 873.0
    fwbs_variables.roughness_fw_channel = 1.0e-6
    fwbs_variables.len_fw_channel = 4.0
    fwbs_variables.f_fw_peak = 1.0
    fwbs_variables.pres_blkt_coolant = 15.50e6
    fwbs_variables.temp_blkt_coolant_in = 573.0
    fwbs_variables.temp_blkt_coolant_out = 823.0
    fwbs_variables.coolp = 15.5e6
    fwbs_variables.n_blkt_outboard_modules_poloidal = 8
    fwbs_variables.n_blkt_inboard_modules_poloidal = 7
    fwbs_variables.n_blkt_outboard_modules_toroidal = 48
    fwbs_variables.n_blkt_inboard_modules_toroidal = 32
    fwbs_variables.temp_fw_max = 823.0
    fwbs_variables.fw_th_conductivity = 28.34
    fwbs_variables.fvoldw = 1.74
    fwbs_variables.fvolsi = 1.0
    fwbs_variables.fvolso = 0.64
    fwbs_variables.fwclfr = 0.15
    fwbs_variables.p_div_rad_total_mw = 0.0
    fwbs_variables.p_fw_rad_total_mw = 0.0
    fwbs_variables.p_fw_hcd_rad_total_mw = 0.0
    fwbs_variables.pradloss = 0.0
    fwbs_variables.p_tf_nuclear_heat_mw = 0.0
    fwbs_variables.ptfnucpm3 = 0.0
    fwbs_variables.r_cryostat_inboard = 0.0
    fwbs_variables.z_cryostat_half_inside = 0.0
    fwbs_variables.dr_pf_cryostat = 0.5
    fwbs_variables.vol_cryostat = 0.0
    fwbs_variables.vol_cryostat_internal = 0.0
    fwbs_variables.vol_vv = 0.0
    fwbs_variables.vfshld = 0.25
    fwbs_variables.vol_blkt_total = 0.0
    fwbs_variables.vol_blkt_inboard = 0.0
    fwbs_variables.vol_blkt_outboard = 0.0
    fwbs_variables.vol_shld_total = 0.0
    fwbs_variables.whtshld = 0.0
    fwbs_variables.wpenshld = 0.0
    fwbs_variables.wtshldi = 0.0
    fwbs_variables.wtshldo = 0.0
    fwbs_variables.irefprop = 1
    fwbs_variables.fblli = 0.0
    fwbs_variables.fblli2o = 0.08
    fwbs_variables.fbllipb = 0.68
    fwbs_variables.fblvd = 0.0
    fwbs_variables.m_blkt_li2o = 0.0
    fwbs_variables.wtbllipb = 0.0
    fwbs_variables.m_blkt_vanadium = 0.0
    fwbs_variables.m_blkt_lithium = 0.0
    fwbs_variables.f_a_blkt_cooling_channels = 0.25
    fwbs_variables.blktmodel = 0
    fwbs_variables.declblkt = 0.075
    fwbs_variables.declfw = 0.075
    fwbs_variables.declshld = 0.075
    fwbs_variables.blkttype = 3
    fwbs_variables.etaiso = 0.85
    fwbs_variables.eta_coolant_pump_electric = 0.95
    fwbs_variables.pnuc_cp = 0.0
    fwbs_variables.p_cp_shield_nuclear_heat_mw = 0.0
    fwbs_variables.pnuc_cp_tf = 0.0
    fwbs_variables.neut_flux_cp = 0.0
    fwbs_variables.i_fw_blkt_shared_coolant = 0
    fwbs_variables.i_blkt_liquid_breeder_type = 0
    fwbs_variables.i_blkt_dual_coolant = 0
    fwbs_variables.i_blkt_liquid_breeder_channel_type = 0
    fwbs_variables.i_blkt_module_segmentation = 0
    fwbs_variables.n_liq_recirc = 10
    fwbs_variables.r_f_liq_ib = 0.5
    fwbs_variables.r_f_liq_ob = 0.5
    fwbs_variables.w_f_liq_ib = 0.5
    fwbs_variables.w_f_liq_ob = 0.5
    fwbs_variables.den_ceramic = 3.21e3
    fwbs_variables.th_wall_secondary = 1.25e-2
    fwbs_variables.bz_channel_conduct_liq = 8.33e5
    fwbs_variables.a_bz_liq = 0.2
    fwbs_variables.b_bz_liq = 0.2
    fwbs_variables.nopol = 2
    fwbs_variables.nopipes = 4
    fwbs_variables.den_liq = 9.5e3
    fwbs_variables.specific_heat_liq = 1.9e2
    fwbs_variables.thermal_conductivity_liq = 30.0
    fwbs_variables.wht_liq = 0.0
    fwbs_variables.wht_liq_ib = 0.0
    fwbs_variables.wht_liq_ob = 0.0
    fwbs_variables.dynamic_viscosity_liq = 0.0
    fwbs_variables.electrical_conductivity_liq = 0.0
    fwbs_variables.hartmann_liq = (0.0, 0.0)
    fwbs_variables.b_mag_blkt = (5.0, 5.0)
    fwbs_variables.etaiso_liq = 0.85
    fwbs_variables.blpressure_liq = 1.7e6
    fwbs_variables.inlet_temp_liq = 570.0
    fwbs_variables.outlet_temp_liq = 720.0
    fwbs_variables.den_fw_coolant = 0.0
    fwbs_variables.visc_fw_coolant = 0.0
    fwbs_variables.den_blkt_coolant = 0.0
    fwbs_variables.visc_blkt_coolant = 0.0
    fwbs_variables.cp_fw = 0.0
    fwbs_variables.cv_fw = 0.0
    fwbs_variables.cp_bl = 0.0
    fwbs_variables.cv_bl = 0.0
    fwbs_variables.f_nuc_pow_bz_struct = 0.34
    fwbs_variables.f_nuc_pow_bz_liq = 0.66
    fwbs_variables.pnuc_fw_ratio_dcll = 0.14
    fwbs_variables.pnuc_blkt_ratio_dcll = 0.86
    fwbs_variables.n_blkt_inboard_module_coolant_sections_radial = 4
    fwbs_variables.n_blkt_inboard_module_coolant_sections_poloidal = 2
    fwbs_variables.bzfllengo_n_rad = 4
    fwbs_variables.bzfllengo_n_pol = 2
    fwbs_variables.bzfllengi_n_rad_liq = 2
    fwbs_variables.bzfllengi_n_pol_liq = 2
    fwbs_variables.bzfllengo_n_rad_liq = 2
    fwbs_variables.bzfllengo_n_pol_liq = 2
