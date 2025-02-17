import numpy as np

from process.coolprop_interface import FluidProperties
from process.fortran import (
    constants,
    error_handling,
    fwbs_variables,
)
from process.fortran import (
    error_handling as eh,
)
from process.fortran import (
    process_output as po,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string


class Fw:
    def __init__(self) -> None:
        self.outfile = constants.nout

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

        Args:
            output (bool): Flag to indicate if output is required.
            radius_fw_channel (float): First wall coolant channel radius (m).
            dr_fw (float): First wall thickness (m).
            a_fw (float): Area of first wall section under consideration (m^2).
            prad_incident (float): Surface heat flux on first wall (MW).
            pnuc_deposited (float): Nuclear power deposited in FW (MW).
            label (str): Information string.

        Returns:
            tuple: Contains peak first wall temperature (K), coolant specific heat capacity at constant pressure (J/kg/K),
               coolant density (kg/m^3), and coolant mass flow rate in a single channel (kg/s).

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

        # kfmean = (kfi + kfo) / 2  # coolant thermal conductivity (W/m.K)
        # viscfmean = (viscfi + viscfo) / 2  # coolant viscosity (Pa.s)

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

        # coolant velocity (m/s)
        vel_fw_coolant_average = mflux_fw_coolant / outlet_coolant_properties.density

        # Mean temperature of the wall material on the plasma side of the coolant 'tpeak'
        # is the estimate from the previous iteration of the wall surface temperature
        # (underneath the armour)
        temp_k = (fwbs_variables.temp_fw_coolant_out + fwbs_variables.tpeak) / 2

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
        onedload = fwbs_variables.peaking_factor * (
            pden_fw_nuclear * fwbs_variables.dx_fw_module * dr_fw / 4
            + pflux_fw_rad * fwbs_variables.dx_fw_module
        )

        # Note I do NOT assume that the channel covers the full width of the first wall:
        # Effective area for heat transfer (m2)
        effective_area_for_heat_transfer = 2 * radius_fw_channel

        # Temperature drop in first-wall material (K)
        deltat_solid_1D = (
            onedload
            * fwbs_variables.dr_fw_wall
            / (tkfw * effective_area_for_heat_transfer)
        )

        # Model C: A more realistic model !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Calculate maximum distance travelled by surface heat load (m)
        # dr_fw_wall | Minimum distance travelled by surface heat load (m)
        diagonal = np.sqrt(
            (fwbs_variables.dx_fw_module / 2 - radius_fw_channel) ** 2
            + (radius_fw_channel + fwbs_variables.dr_fw_wall) ** 2
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
                "(peaking_factor)",
                fwbs_variables.peaking_factor,
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

    def fw_thermal_conductivity(self, t):
        """Calculates the thermal conductivity of the first wall
        t : input real : property temperature (K)
        Calculates the thermal conductivity of Eurofer (W/m/K).
        """
        # Eurofer correlation, from "Fusion Demo Interim Structural Design Criteria -
        # Appendix A Material Design Limit Data", F. Tavassoli, TW4-TTMS-005-D01, 2004
        # t in Kelvin
        return (
            (5.4308 + 0.13565 * t - 0.00023862 * t * t + 1.3393e-7 * t * t * t)
            * fwbs_variables.fw_th_conductivity
            / 28.34
        )

    def heat_transfer(self, masflx, rhof, radius, cf, viscf, kf):
        """Calculate heat transfer coefficient using Gnielinski correlation
        author: M Kovari, CCFE, Culham Science Centre
        masflx : input real : coolant mass flux in a single channel (kg/m2/s)
        rhof : input real : coolant density (average of inlet and outlet) (kg/m3)
        radius : input real : coolant pipe radius (m)
        cf : input real : coolant specific heat capacity (average of inlet and outlet) (J/K)
        viscf : input real : coolant viscosity (average of inlet and outlet) (Pa.s)
        kf : input real : thermal conductivity of coolant (average of inlet and outlet) (W/m.K)
        Gnielinski correlation. Ignore the distinction between wall and
        bulk temperatures. Valid for:3000 < Re < 5e6, 0.5 < Pr < 2000
        https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation
        """
        # Calculate pipe diameter (m)
        diameter = 2 * radius

        # Calculate flow velocity (m/s)
        velocity = masflx / rhof

        # Calculate Reynolds number
        reynolds = rhof * velocity * diameter / viscf

        # Calculate Prandtl number
        pr = cf * viscf / kf

        # Calculate Darcy friction factor, using Haaland equation
        f = self.friction(reynolds)

        # Calculate the Nusselt number
        nusselt = (
            (f / 8.0)
            * (reynolds - 1000.0)
            * pr
            / (1 + 12.7 * np.sqrt(f / 8.0) * (pr**0.6667 - 1.0))
        )

        # Calculate the heat transfer coefficient (W/m2K)
        heat_transfer = nusselt * kf / (2.0 * radius)

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

        return heat_transfer

    def friction(self, reynolds):
        """Calculate Darcy friction factor, using Haaland equation
        author: M Kovari, CCFE, Culham Science Centre
        reynolds : input real : Reynolds number
        darcy_friction : output real : Darcy friction factor
        Darcy friction factor, using Haaland equation, an approximation to the
        implicit Colebrook-White equationGnielinski correlation.
        https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
        """

        # Bracketed term in Haaland equation
        bracket = (
            fwbs_variables.roughness / fwbs_variables.radius_fw_channel / 3.7
        ) ** 1.11 + 6.9 / reynolds

        # Calculate Darcy friction factor
        return (1.8 * np.log10(bracket)) ** (-2)
