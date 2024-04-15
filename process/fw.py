import numpy as np

from process.fortran import (
    constants,
    fwbs_variables,
    error_handling as eh,
    fw_module,
    process_output as po,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from process.coolprop_interface import FluidProperties


class Fw:
    def __init__(self) -> None:
        self.outfile = constants.nout

    def fw_temp(
        self,
        output: bool,
        afw,
        thickness,
        area,
        prad_incident,
        pnuc_deposited,
        label,
    ):
        """Thermo-hydraulic calculations for the first wall
        author: P J Knight, CCFE, Culham Science Centre
        afw : input real : first wall coolant channel radius (m)
        thickness : first wall thickness (fwith or fwoth) (m)
        area : input real : area of first wall section under consideration (m2)
        (i.e. area of inboard wall or outboard wall)
        prad_incident : input real : Surface heat flux on first wall (outboard and inboard) (MW)
        pnuc_deposited : input real : nuclear power deposited in FW (IB or OB) (MW)
        tpeakfw : output real : peak first wall temperature (K)
        cfmean : output real : coolant specific heat capacity at constant
        pressure (J/kg/K)
        rhofmean : output real : coolant density (kg/m3)
        massrate : output real : coolant mass flow rate in a single channel (kg/s)
        label : input string : information string
        Detailed thermal hydraulic model for the blanket (first wall +
        breeding zone).
        Given the heating incident on the first wall, and the coolant
        outlet temperature, the maximum temperature of the first wall is
        calculated to check it is below material limits (tfwmatmax).
        The routine is called separately for the inboard and outboard sides.
        The calculation of the maximum temperature is described by Gardner:
        "Temperature distribution in the first wall", K:\\Power Plant Physics and
        Technology\\ PROCESS\\PROCESS References & Systems Codes\\Pulsed option -
        Gardner.
        This is in turn taken from "Methods of First Wall Structural
        Analysis with Application to the Long Pulse Commercial Tokamak Reactor
        Design", R.J. LeClaire, MIT, PFC/RR-84-9
        """
        # First wall volume (inboard or outboard depending on arguments) (m3)
        fwvol = area * thickness

        # First wall channel area (m2)
        channel_area = np.pi * afw**2

        # Heat generation in the first wall due to neutron flux deposited in the material (W/m3)
        qppp = 1e6 * pnuc_deposited / fwvol

        # the nuclear heating in the coolant is small. (W/m2)
        # Note that the full first wall volume is used including coolant even though
        nuclear_heat_per_area = qppp * thickness

        # Heat flux incident on the first wall surface (W/m2)
        qpp = 1e6 * prad_incident / area

        # Calculate inlet coolant fluid properties (fixed pressure)
        ib_fluid_properties = FluidProperties.of(
            f2py_compatible_to_string(fwbs_variables.fwcoolant),
            temperature=fwbs_variables.fwinlet.item(),
            pressure=fwbs_variables.fwpressure.item(),
        )

        # Calculate outlet coolant fluid properties (fixed pressure)
        ob_fluid_properties = FluidProperties.of(
            f2py_compatible_to_string(fwbs_variables.fwcoolant),
            temperature=fwbs_variables.fwoutlet.item(),
            pressure=fwbs_variables.fwpressure.item(),
        )

        # Mean properties (inlet + outlet)/2
        rhofmean = (
            ib_fluid_properties.density + ob_fluid_properties.density
        ) / 2  # coolant density (kg/m3)
        # kfmean = (kfi + kfo) / 2  # coolant thermal conductivity (W/m.K)
        # viscfmean = (viscfi + viscfo) / 2  # coolant viscosity (Pa.s)
        cfmean = (
            ib_fluid_properties.specific_heat_const_p
            + ob_fluid_properties.specific_heat_const_p
        ) / 2  # coolant specific heat capacity (J/K)

        # Heat load per unit length of one first wall pipe (W/m)
        load = (nuclear_heat_per_area + qpp) * fwbs_variables.pitch

        # Coolant mass flow rate (kg/s) (use mean properties)
        massrate = (
            fwbs_variables.fw_channel_length
            * load
            / cfmean
            / (fwbs_variables.fwoutlet - fwbs_variables.fwinlet)
        )

        # Coolant mass flux in a single channel (kg/m2/s)
        masflx = massrate / channel_area

        # Conditions at the outlet, where the temperature is highest
        # -----------------------------------------------------------

        # Outlet coolant velocity (m/s)
        velocity = masflx / ob_fluid_properties.density

        # Mean temperature of the wall material on the plasma side of the coolant 'tpeak'
        # is the estimate from the previous iteration of the wall surface temperature
        # (underneath the armour)
        temp_k = (fwbs_variables.fwoutlet + fwbs_variables.tpeak) / 2  # (K)

        # Print debug info if temperature too low/high or NaN/Inf
        if np.isnan(temp_k):
            eh.report_error(223)
        elif (temp_k <= 100) or (temp_k > 1500):
            eh.fdiags[0] = temp_k
            eh.report_error(224)

        # Thermal conductivity of first wall material (W/m.K)
        tkfw = fw_module.fw_thermal_conductivity(temp_k)

        # Heat transfer coefficient (W/m2K)
        hcoeff = fw_module.heat_transfer(
            masflx,
            ob_fluid_properties.density,
            afw,
            ob_fluid_properties.specific_heat_const_p,
            ob_fluid_properties.viscosity,
            ob_fluid_properties.thermal_conductivity,
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
            qppp * fwbs_variables.pitch * thickness / 4 + qpp * fwbs_variables.pitch
        )

        # Note I do NOT assume that the channel covers the full width of the first wall:
        # Effective area for heat transfer (m2)
        effective_area_for_heat_transfer = 2 * afw

        # Temperature drop in first-wall material (K)
        deltat_solid_1D = (
            onedload
            * fwbs_variables.fw_wall
            / (tkfw * effective_area_for_heat_transfer)
        )

        # Model C: A more realistic model !
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Calculate maximum distance travelled by surface heat load (m)
        # fw_wall | Minimum distance travelled by surface heat load (m)
        diagonal = np.sqrt(
            (fwbs_variables.pitch / 2 - afw) ** 2 + (afw + fwbs_variables.fw_wall) ** 2
        )

        # Mean distance travelled by surface heat (m)
        mean_distance = (fwbs_variables.fw_wall + diagonal) / 2

        # This heat starts off spread over width = 'pitch'.
        # It ends up spread over one half the circumference.
        # Use the mean of these values.
        mean_width = (fwbs_variables.pitch + np.pi * afw) / 2  # (m)

        # As before, use a combined load 'onedload'
        # Temperature drop in first-wall material (K)
        deltat_solid = onedload * mean_distance / (tkfw * mean_width)

        # Temperature drop between channel inner wall and bulk coolant (K)
        deltat_coolant = load / (2 * np.pi * afw * hcoeff)

        # Peak first wall temperature (K)
        tpeakfw = fwbs_variables.fwoutlet + deltat_solid + deltat_coolant

        if output:
            po.oheadr(
                self.outfile, "Heat transfer parameters at the coolant outlet: " + label
            )
            po.ovarre(self.outfile, "Radius of coolant channel (m)", "(afw)", afw)
            po.ovarre(
                self.outfile,
                "Mean surface heat flux on first wall (W/m2) ",
                "(qpp)",
                qpp,
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
                "Length of a single coolant channel (all in parallel) (m)",
                "(fw_channel_length)",
                fwbs_variables.fw_channel_length,
            )
            po.ovarre(
                self.outfile,
                "Pitch of coolant channels (m)",
                "(pitch)",
                fwbs_variables.pitch,
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
                ob_fluid_properties.density,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Coolant mass flow rate in one channel (kg/s)",
                "(massrate)",
                massrate,
                "OP ",
            )
            po.ovarre(
                self.outfile, "Coolant velocity (m/s)", "(velocity)", velocity, "OP "
            )
            po.ovarre(
                self.outfile,
                "Outlet temperature of first wall coolant (K)",
                "(fwoutlet)",
                fwbs_variables.fwoutlet,
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

        return tpeakfw, cfmean, rhofmean, massrate
