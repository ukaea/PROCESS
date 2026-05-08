import logging

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.coolprop_interface import FluidProperties
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure import (
    blanket_library,
    build_variables,
    constraint_variables,
    divertor_variables,
    physics_variables,
)
from process.models.build import FwBlktVVShape
from process.models.engineering.ivc_functions import (
    calculate_pipe_bend_radius,
    dshellarea,
    eshellarea,
)
from process.models.engineering.materials import eurofer97_thermal_conductivity
from process.models.engineering.pumping import (
    gnielinski_heat_transfer_coefficient,
)

logger = logging.getLogger(__name__)


class FirstWall(Model):
    def __init__(self):
        self.outfile = constants.NOUT

    def output(self):
        # First wall geometry
        self.output_fw_geometry()

        # First wall surface loads
        self.output_fw_surface_loads()

        # First wall pumping
        self.output_fw_pumping()

    def run(self):
        self.data.fwbs.dz_fw_half = self.calculate_first_wall_half_height(
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

        if (
            physics_variables.itart == 1
            or self.data.fwbs.i_fw_blkt_vv_shape == FwBlktVVShape.D_SHAPED
        ):
            (
                self.data.first_wall.a_fw_inboard_full_coverage,
                self.data.first_wall.a_fw_outboard_full_coverage,
                self.data.first_wall.a_fw_total_full_coverage,
            ) = self.calculate_dshaped_first_wall_areas(
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                dz_fw_half=self.data.fwbs.dz_fw_half,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            )

        else:
            (
                self.data.first_wall.a_fw_inboard_full_coverage,
                self.data.first_wall.a_fw_outboard_full_coverage,
                self.data.first_wall.a_fw_total_full_coverage,
            ) = self.calculate_elliptical_first_wall_areas(
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                triang=physics_variables.triang,
                dz_fw_half=self.data.fwbs.dz_fw_half,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            )

        (
            self.data.first_wall.a_fw_inboard,
            self.data.first_wall.a_fw_outboard,
            self.data.first_wall.a_fw_total,
        ) = self.apply_first_wall_coverage_factors(
            n_divertors=divertor_variables.n_divertors,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            f_a_fw_outboard_hcd=self.data.fwbs.f_a_fw_outboard_hcd,
            a_fw_inboard_full_coverage=self.data.first_wall.a_fw_inboard_full_coverage,
            a_fw_outboard_full_coverage=self.data.first_wall.a_fw_outboard_full_coverage,
        )

        (
            blanket_library.n_fw_inboard_channels,
            blanket_library.n_fw_outboard_channels,
        ) = self.calculate_total_fw_channels(
            self.data.first_wall.a_fw_inboard,
            self.data.first_wall.a_fw_outboard,
            self.data.fwbs.len_fw_channel,
            self.data.fwbs.dx_fw_module,
        )

        self.set_fw_geometry()

        (
            self.data.fwbs.radius_fw_channel_90_bend,
            self.data.fwbs.radius_fw_channel_180_bend,
        ) = calculate_pipe_bend_radius(
            i_ps=1,
            radius_fw_channel=self.data.fwbs.radius_fw_channel,
            b_bz_liq=self.data.fwbs.b_bz_liq,
        )

        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.ffwal
                * physics_variables.pflux_plasma_surface_neutron_avg_mw
            )
        else:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.p_neutron_total_mw / self.data.first_wall.a_fw_total
            )

        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_rad_mw = (
                physics_variables.ffwal
                * physics_variables.p_plasma_rad_mw
                / physics_variables.a_plasma_surface
            )
        else:
            physics_variables.pflux_fw_rad_mw = (
                physics_variables.p_plasma_rad_mw / self.data.first_wall.a_fw_total
            )

        constraint_variables.pflux_fw_rad_max_mw = (
            physics_variables.pflux_fw_rad_mw * constraint_variables.f_fw_rad_max
        )

        # Power transported to the first wall by escaped alpha particles
        physics_variables.p_fw_alpha_mw = physics_variables.p_alpha_total_mw * (
            1.0e0 - physics_variables.f_p_alpha_plasma_deposited
        )

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
        """Calculate the half-height of the first wall.

        Parameters
        ----------
        z_plasma_xpoint_lower:

        dz_xpoint_divertor:

        dz_divertor:

        dz_blkt_upper:

        z_plasma_xpoint_upper:

        dz_fw_plasma_gap:

        n_divertors: int :

        dr_fw_inboard:

        dr_fw_outboard:

        """
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
        """Calculate the first wall areas for an elliptical cross-section.

        Parameters
        ----------
        rmajor:

        rminor:

        triang:

        dz_fw_half:

        dr_fw_plasma_gap_inboard:

        dr_fw_plasma_gap_outboard:

        """
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
            2 * self.data.fwbs.radius_fw_channel + 2 * self.data.fwbs.dr_fw_wall
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
        """Thermo-hydraulic calculations for the first wall.

        Parameters
        ----------
        output:
            Flag to indicate if output is required.
        radius_fw_channel:
            First wall coolant channel radius (m).
        dr_fw:
            First wall thickness (m).
        a_fw:
            Area of first wall section under consideration (m^2).
        prad_incident:
            Radiation surface heat flux on first wall (MW).
        pnuc_deposited:
            Nuclear power deposited in FW (MW).
        label:
            Information string.

        Returns
        -------
        tuple

        Detailed thermal hydraulic model for the blanket (first wall + breeding zone).
        Given the heating incident on the first wall, and the coolant outlet temperature,
        the maximum temperature of the first wall is calculated to check it is below material limits.
        The routine is called separately for the inboard and outboard sides.
        The calculation of the maximum temperature is described by Gardner:
        "Temperature distribution in the first wall", K:\\Power Plant Physics and Technology\\ PROCESS\\PROCESS References & Systems Codes\\Pulsed option - Gardner.
        This is in turn taken from "Methods of First Wall Structural Analysis with Application to the Long Pulse Commercial Tokamak Reactor Design", R.J. LeClaire, MIT, PFC/RR-84-9.
        Contains peak first wall temperature (K), coolant specific heat capacity at constant pressure (J/kg/K),
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
            self.data.fwbs.i_fw_coolant_type,
            temperature=self.data.fwbs.temp_fw_coolant_in,
            pressure=self.data.fwbs.pres_fw_coolant,
        )

        # Calculate outlet coolant fluid properties (fixed pressure)
        outlet_coolant_properties = FluidProperties.of(
            self.data.fwbs.i_fw_coolant_type,
            temperature=self.data.fwbs.temp_fw_coolant_out,
            pressure=self.data.fwbs.pres_fw_coolant,
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
        load = (nuclear_heat_per_area + pflux_fw_rad) * self.data.fwbs.dx_fw_module

        # Coolant mass flow rate (kg/s) (use mean properties)
        mflow_fw_coolant = (
            self.data.fwbs.len_fw_channel
            * load
            / heatcap_fw_coolant_average
            / (self.data.fwbs.temp_fw_coolant_out - self.data.fwbs.temp_fw_coolant_in)
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
        temp_k = (self.data.fwbs.temp_fw_coolant_out + self.data.fwbs.temp_fw_peak) / 2

        # Print debug info if temperature too low/high or NaN/Inf
        if np.isnan(temp_k):
            logger.error("NaN first wall temperature")
        elif (temp_k <= 100) or (temp_k > 1500):
            logger.error(
                "First wall temperature (temp_k) out of range : [100-1500] K. %s", temp_k
            )

        # Thermal conductivity of first wall material (W/m.K)
        tkfw = eurofer97_thermal_conductivity(
            temp=temp_k, fw_th_conductivity=self.data.fwbs.fw_th_conductivity
        )

        # Heat transfer coefficient (W m^-2 K^-1)
        hcoeff = gnielinski_heat_transfer_coefficient(
            mflux_coolant=mflux_fw_coolant,
            den_coolant=outlet_coolant_properties.density,
            radius_channel=radius_fw_channel,
            heatcap_coolant=outlet_coolant_properties.specific_heat_const_p,
            visc_coolant=outlet_coolant_properties.viscosity,
            thermcond_coolant=outlet_coolant_properties.thermal_conductivity,
            roughness_channel=self.data.fwbs.roughness_fw_channel,
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
        onedload = self.data.fwbs.f_fw_peak * (
            pden_fw_nuclear * self.data.fwbs.dx_fw_module * dr_fw / 4
            + pflux_fw_rad * self.data.fwbs.dx_fw_module
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
            (self.data.fwbs.dx_fw_module / 2 - radius_fw_channel) ** 2
            + (radius_fw_channel + dr_fw) ** 2
        )

        # Mean distance travelled by surface heat (m)
        mean_distance = (self.data.fwbs.dr_fw_wall + diagonal) / 2

        # This heat starts off spread over width = 'dx_fw_module'.
        # It ends up spread over one half the circumference.
        # Use the mean of these values.
        mean_width = (self.data.fwbs.dx_fw_module + np.pi * radius_fw_channel) / 2  # (m)

        # As before, use a combined load 'onedload'
        # Temperature drop in first-wall material (K)
        deltat_solid = onedload * mean_distance / (tkfw * mean_width)

        # Temperature drop between channel inner wall and bulk coolant (K)
        deltat_coolant = load / (2 * np.pi * radius_fw_channel * hcoeff)

        # Peak first wall temperature (K)
        tpeakfw = self.data.fwbs.temp_fw_coolant_out + deltat_solid + deltat_coolant

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
                self.data.fwbs.f_fw_peak,
            )
            po.ovarre(
                self.outfile,
                "Vertical length of a single coolant channel (all in parallel) (m)",
                "(len_fw_channel)",
                self.data.fwbs.len_fw_channel,
            )
            po.ovarre(
                self.outfile,
                "Width of a FW module containing a cooling channel [m]",
                "(dx_fw_module)",
                self.data.fwbs.dx_fw_module,
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
                self.data.fwbs.temp_fw_coolant_out,
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

    @staticmethod
    def calculate_total_fw_channels(
        a_fw_inboard: float,
        a_fw_outboard: float,
        len_fw_channel: float,
        dx_fw_module: float,
    ) -> tuple[int, int]:
        """Calculate the total number of first wall channels for inboard and outboard sections.

        Parameters
        ----------
        a_fw_inboard:
            Area of the inboard first wall section (m^2).
        a_fw_outboard:
            Area of the outboard first wall section (m^2).
        len_fw_channel:
            Length of each first wall channel (m).
        dx_fw_module:
            Toroidal width of each first wall module (m).

        Returns
        -------
        tuple
            Number of inboard and outboard first wall channels.
        """
        n_fw_inboard_channels = a_fw_inboard / (len_fw_channel * dx_fw_module)
        n_fw_outboard_channels = a_fw_outboard / (len_fw_channel * dx_fw_module)
        return int(n_fw_inboard_channels), int(n_fw_outboard_channels)

    def output_fw_geometry(self):
        """Outputs the first wall geometry details to the output file."""
        po.oheadr(self.outfile, "First wall build")

        po.ovarrf(
            self.outfile,
            "Radius of first wall cooling channels (m)",
            "(radius_fw_channel)",
            self.data.fwbs.radius_fw_channel,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Radius of 90 degree coolant channel bend (m)",
            "(radius_fw_channel_90_bend)",
            self.data.fwbs.radius_fw_channel_90_bend,
        )
        po.ovarre(
            self.outfile,
            "Radius of 180 degree coolant channel bend (m)",
            "(radius_fw_channel_180_bend)",
            self.data.fwbs.radius_fw_channel_180_bend,
        )
        po.ovarrf(
            self.outfile,
            "Radial wall thickness surrounding first wall coolant channel (m)",
            "(dr_fw_wall)",
            self.data.fwbs.dr_fw_wall,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Toroidal width of each first wall module (m)",
            "(dx_fw_module)",
            self.data.fwbs.dx_fw_module,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Length of each first wall channel (m)",
            "(len_fw_channel)",
            self.data.fwbs.len_fw_channel,
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
        """Outputs the first wall pumping details to the output file."""
        po.oheadr(self.outfile, "First wall pumping")

        po.ovarst(
            self.outfile,
            "First wall coolant type",
            "(i_fw_coolant_type)",
            f"'{self.data.fwbs.i_fw_coolant_type}'",
        )
        po.ovarrf(
            self.outfile,
            "Outlet temperature of first wall coolant [K]",
            "(temp_fw_coolant_out)",
            self.data.fwbs.temp_fw_coolant_out,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Inlet temperature of first wall coolant [K]",
            "(temp_fw_coolant_in)",
            self.data.fwbs.temp_fw_coolant_in,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pressure of first wall coolant [Pa]",
            "(pres_fw_coolant)",
            self.data.fwbs.pres_fw_coolant,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Peak temperature of first wall [K]",
            "(temp_fw_peak)",
            self.data.fwbs.temp_fw_peak,
            "OP ",
        )

    def output_fw_surface_loads(self):
        """Outputs the first wall surface load details to the output file."""
        po.oheadr(self.outfile, "First wall surface loads")

        po.ovarre(
            self.outfile,
            "Nominal mean radiation load on vessel first-wall (MW/m^2)",
            "(pflux_fw_rad_mw)",
            physics_variables.pflux_fw_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Peaking factor for radiation first-wall load",
            "(f_fw_rad_max)",
            constraint_variables.f_fw_rad_max,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum permitted radiation first-wall load (MW/m^2)",
            "(pflux_fw_rad_max)",
            constraint_variables.pflux_fw_rad_max,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Peak radiation wall load (MW/m^2)",
            "(pflux_fw_rad_max_mw)",
            constraint_variables.pflux_fw_rad_max_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fast alpha particle power incident on the first-wall (MW)",
            "(p_fw_alpha_mw)",
            physics_variables.p_fw_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean neutron load on vessel first-wall (MW/m^2)",
            "(pflux_fw_neutron_mw)",
            physics_variables.pflux_fw_neutron_mw,
            "OP ",
        )
