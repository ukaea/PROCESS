import logging

import numpy as np

from process import constants
from process import process_output as po
from process.blanket_library import BlanketLibrary, dshellarea, eshellarea
from process.coolprop_interface import FluidProperties
from process.data_structure import (
    blanket_library,
    build_variables,
    divertor_variables,
    first_wall_variables,
    fwbs_variables,
    physics_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class FirstWall:
    def __init__(self) -> None:
        self.outfile = constants.NOUT
        self.blanket_library = BlanketLibrary(fw=self)

    def run(self):
        fwbs_variables.dz_fw_half = self.calculate_first_wall_half_height(
            z_plasma_xpoint_lower=build_variables.z_plasma_xpoint_lower,
            dz_xpoint_divertor=build_variables.dz_xpoint_divertor,
            dz_divertor=divertor_variables.dz_divertor,
            dz_blkt_upper=build_variables.dz_blkt_upper,
            z_plasma_xpoint_upper=build_variables.z_plasma_xpoint_upper,
            dz_fw_plasma_gap=build_variables.dz_fw_plasma_gap,
            n_divertors=divertor_variables.n_divertors,
            dr_fw_inboard=build_variables.dr_fw_inboard,
            dr_fw_outboard=build_variables.dr_fw_outboard,
        )

        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            (
                first_wall_variables.a_fw_inboard_full_coverage,
                first_wall_variables.a_fw_outboard_full_coverage,
                first_wall_variables.a_fw_total_full_coverage,
            ) = self.calculate_dshaped_first_wall_areas(
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                dz_fw_half=fwbs_variables.dz_fw_half,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            )

        else:
            (
                first_wall_variables.a_fw_inboard_full_coverage,
                first_wall_variables.a_fw_outboard_full_coverage,
                first_wall_variables.a_fw_total_full_coverage,
            ) = self.calculate_elliptical_first_wall_areas(
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                triang=physics_variables.triang,
                dz_fw_half=fwbs_variables.dz_fw_half,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            )

        (
            first_wall_variables.a_fw_inboard,
            first_wall_variables.a_fw_outboard,
            first_wall_variables.a_fw_total,
        ) = self.apply_first_wall_coverage_factors(
            n_divertors=divertor_variables.n_divertors,
            f_ster_div_single=fwbs_variables.f_ster_div_single,
            f_a_fw_outboard_hcd=fwbs_variables.f_a_fw_outboard_hcd,
            a_fw_inboard_full_coverage=first_wall_variables.a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage=first_wall_variables.a_fw_outboard_full_coverage,
        )

        (
            blanket_library.n_fw_inboard_channels,
            blanket_library.n_fw_outboard_channels,
        ) = self.calculate_total_fw_channels(
            first_wall_variables.a_fw_inboard,
            first_wall_variables.a_fw_outboard,
            fwbs_variables.len_fw_channel,
            fwbs_variables.dx_fw_module,
        )

        self.set_fw_geometry()

        (
            fwbs_variables.radius_fw_channel_90_bend,
            fwbs_variables.radius_fw_channel_180_bend,
        ) = self.blanket_library.calculate_pipe_bend_radius(i_ps=1)

    @staticmethod
    def calculate_first_wall_half_height(
        z_plasma_xpoint_lower: float,
        dz_xpoint_divertor: float,
        dz_divertor: float,
        dz_blkt_upper: float,
        z_plasma_xpoint_upper: float,
        dz_fw_plasma_gap: float,
        n_divertors: int,
        dr_fw_inboard: float,
        dr_fw_outboard: float,
    ) -> float:
        """Calculate the half-height of the first wall."""

        #  Half-height of first wall (internal surface)
        z_bottom = (
            z_plasma_xpoint_lower
            + dz_xpoint_divertor
            + dz_divertor
            - dz_blkt_upper
            - 0.5e0 * (dr_fw_inboard + dr_fw_outboard)
        )
        if n_divertors == 2:
            z_top = z_bottom
        else:
            z_top = z_plasma_xpoint_upper + dz_fw_plasma_gap

        return 0.5e0 * (z_top + z_bottom)

    @staticmethod
    def calculate_dshaped_first_wall_areas(
        rmajor: float,
        rminor: float,
        dz_fw_half: float,
        dr_fw_plasma_gap_inboard: float,
        dr_fw_plasma_gap_outboard: float,
    ) -> tuple[float, float, float]:
        # D-shaped
        #  Major radius to outer edge of inboard section
        r1 = rmajor - rminor - dr_fw_plasma_gap_inboard

        #  Horizontal distance between inside edges,
        #  i.e. outer radius of inboard part to inner radius of outboard part

        r2 = (rmajor + rminor + dr_fw_plasma_gap_outboard) - r1
        #  Calculate surface area, assuming 100% coverage

        (
            a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage,
            a_fw_total_full_coverage,
        ) = dshellarea(rmajor=r1, rminor=r2, zminor=dz_fw_half)

        return (
            a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage,
            a_fw_total_full_coverage,
        )

    @staticmethod
    def calculate_elliptical_first_wall_areas(
        rmajor: float,
        rminor: float,
        triang: float,
        dz_fw_half: float,
        dr_fw_plasma_gap_inboard: float,
        dr_fw_plasma_gap_outboard: float,
    ) -> tuple[float, float, float]:
        """Calculate the first wall areas for an elliptical cross-section."""

        # Cross-section is assumed to be defined by two ellipses
        #  Major radius to centre of inboard and outboard ellipses
        #  (coincident in radius with top of plasma)

        r1 = rmajor - rminor * triang

        #  Distance between r1 and outer edge of inboard section

        r2 = r1 - (rmajor - rminor - dr_fw_plasma_gap_inboard)

        #  Distance between r1 and inner edge of outboard section

        r3 = (rmajor + rminor + dr_fw_plasma_gap_outboard) - r1

        #  Calculate surface area, assuming 100% coverage

        (
            a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage,
            a_fw_total_full_coverage,
        ) = eshellarea(rshell=r1, rmini=r2, rmino=r3, zminor=dz_fw_half)

        return (
            a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage,
            a_fw_total_full_coverage,
        )

    @staticmethod
    def apply_first_wall_coverage_factors(
        n_divertors: int,
        f_ster_div_single: float,
        f_a_fw_outboard_hcd: float,
        a_fw_inboard_full_coverage: float,
        a_fw_outboard_full_coverage: float,
    ) -> tuple[float, float, float]:
        """Apply first wall coverage factors to calculate actual first wall areas.

        Parameters
        ----------
        n_divertors : int
            Number of divertors (1 or 2).
        f_ster_div_single : float
            Fractional area of first wall sterically blocked by single divertor.
        f_a_fw_outboard_hcd : float
            Fractional area of outboard first wall covered by high heat flux components.
        a_fw_inboard_full_coverage : float
            First wall inboard area assuming 100% coverage (m^2).
        a_fw_outboard_full_coverage : float
            First wall outboard area assuming 100% coverage (m^2).
        n_divertors: int :

        f_ster_div_single: float :

        f_a_fw_outboard_hcd: float :

        a_fw_inboard_full_coverage: float :

        a_fw_outboard_full_coverage: float :


        Returns
        -------
        tuple[float, float, float]
            Contains first wall inboard area, outboard area, and total area (m^2).
        """
        if n_divertors == 2:
            # Double null configuration
            a_fw_outboard = a_fw_outboard_full_coverage * (
                1.0e0 - 2.0e0 * f_ster_div_single - f_a_fw_outboard_hcd
            )
            a_fw_inboard = a_fw_inboard_full_coverage * (
                1.0e0 - 2.0e0 * f_ster_div_single
            )
        else:
            # Single null configuration
            a_fw_outboard = a_fw_outboard_full_coverage * (
                1.0e0 - f_ster_div_single - f_a_fw_outboard_hcd
            )
            a_fw_inboard = a_fw_inboard_full_coverage * (1.0e0 - f_ster_div_single)

        a_fw_total = a_fw_inboard + a_fw_outboard

        if a_fw_outboard <= 0.0e0:
            raise ProcessValueError(
                "fhole+f_ster_div_single+f_a_fw_outboard_hcd is too high for a credible outboard wall area",
                f_ster_div_single=f_ster_div_single,
                f_a_fw_outboard_hcd=f_a_fw_outboard_hcd,
            )

        return a_fw_inboard, a_fw_outboard, a_fw_total

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
            fwbs_variables.i_fw_coolant_type,
            temperature=fwbs_variables.temp_fw_coolant_in,
            pressure=fwbs_variables.pres_fw_coolant,
        )

        # Calculate outlet coolant fluid properties (fixed pressure)
        outlet_coolant_properties = FluidProperties.of(
            fwbs_variables.i_fw_coolant_type,
            temperature=fwbs_variables.temp_fw_coolant_out,
            pressure=fwbs_variables.pres_fw_coolant,
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
            logger.error("NaN first wall temperature")
        elif (temp_k <= 100) or (temp_k > 1500):
            logger.error(
                f"First wall temperature (temp_k) out of range : [100-1500] K. {temp_k=}"
            )

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
        mean_width = (fwbs_variables.dx_fw_module + np.pi * radius_fw_channel) / 2  # (m)

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
            logger.error(f"Reynolds number out of range : [3e3-5000e3]. {reynolds=}")

        # Check that Prandtl number is in valid range for the Gnielinski correlation
        if (pr < 0.5) or (pr > 2000.0):
            logger.error(f"Prandtl number out of range : [0.5-2000]. {pr=}")

        # Check that the Darcy friction factor is in valid range for the Gnielinski correlation
        if f <= 0.0:
            logger.error(f"Negative Darcy friction factor (f). {f=}")

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
        po.ovarre(
            self.outfile,
            "Radius of 90 degree coolant channel bend (m)",
            "(radius_fw_channel_90_bend)",
            fwbs_variables.radius_fw_channel_90_bend,
        )
        po.ovarre(
            self.outfile,
            "Radius of 180 degree coolant channel bend (m)",
            "(radius_fw_channel_180_bend)",
            fwbs_variables.radius_fw_channel_180_bend,
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
            f"'{fwbs_variables.i_fw_coolant_type}'",
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
