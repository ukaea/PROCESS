"""Module containing divertor routines"""

import math
from dataclasses import dataclass

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure.divertor_variables import DivertorHeatLoadModel
from process.data_structure.physics_variables import DivertorNumberModels


@dataclass
class WadeDivertorMetrics:
    """Dataclass to hold the output metrics for the Wade divertor model."""

    pflux_div_heat_load_mw: float
    """Divertor heat load [MW/m²]"""
    dx_div_lower_outboard_strike: float
    """Lower outboard divertor strike point width [m]"""
    a_div_lower_outboard_wetted: float
    """Lower divertor outboard wetted area [m²]"""
    f_div_lower_flux_expansion: float
    """Lower divertor flux expansion factor (fₓ)"""
    deg_b_div_lower_flux: float
    """Lower divertor flux angle [degrees]"""
    deg_div_lower_outboard_plate_separatrix_poloidal: float
    """Lower divertor outboard plate-separatrix poloidal angle [degrees]"""
    deg_b_div_lower_outboard_grazing: float
    """Lower divertor outboard grazing angle [degrees]"""


class Divertor(Model):
    """Module containing divertor routines

    This module contains routines relevant for calculating the
    divertor parameters for a fusion power plant.
    """

    def __init__(self):
        self.outfile = constants.NOUT  # output file unit

    def output(self):
        """Output divertor information"""
        self.run(output=True)

    def run(self, output: bool = False):
        """Routine to call the divertor model

        This subroutine calls the divertor routine. This routine scales
        dimensions, powers and field levels which are used as input to
        the Harrison divertor model.

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """
        self.data.divertor.deg_div_poloidal_plasma = self.single_divertor_angle
        self.data.fwbs.f_ster_div_single = (
            self.data.divertor.deg_div_poloidal_plasma / 360.0
        )

        self.data.fwbs.p_div_nuclear_heat_total_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=self.data.physics.p_plasma_neutron_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=self.data.divertor.n_divertors,
        )

        self.data.fwbs.p_div_rad_total_mw = self.incident_radiation_power(
            p_plasma_rad_mw=self.data.physics.p_plasma_rad_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=self.data.divertor.n_divertors,
        )

        if (
            DivertorHeatLoadModel(self.data.divertor.i_div_heat_load)
            == DivertorHeatLoadModel.USER_INPUT
            and output
        ):
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m²)",
                "(pflux_div_heat_load_mw)",
                self.data.divertor.pflux_div_heat_load_mw,
            )
            return
        if (
            DivertorHeatLoadModel(self.data.divertor.i_div_heat_load)
            == DivertorHeatLoadModel.PENG_CHAMBER
        ):
            self.divtart(
                self.data.physics.rmajor,
                self.data.physics.rminor,
                self.data.physics.triang,
                self.data.build.dr_fw_plasma_gap_inboard,
                self.data.build.dz_xpoint_divertor,
                self.data.physics.p_plasma_separatrix_mw,
                output=output,
                i_single_null=self.data.physics.i_single_null,
                dz_divertor=self.data.divertor.dz_divertor,
            )
            return
        if (
            DivertorHeatLoadModel(self.data.divertor.i_div_heat_load)
            == DivertorHeatLoadModel.WADE
        ):
            wade_divertor = self.divwade(
                rmajor=self.data.physics.rmajor,
                p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
                len_plasma_sol_power_decay=self.data.physics.len_plasma_sol_eich13_power_decay,
                len_plasma_sol_power_spreading=self.data.physics.len_plasma_sol_scrabosio14_power_spreading,
                f_div_flux_expansion=self.data.divertor.f_div_flux_expansion,
                deg_b_div_lower_outboard_grazing=self.data.divertor.deg_b_div_lower_outboard_grazing,
                rad_fraction_sol=self.data.physics.rad_fraction_sol,
                f_p_div_lower=self.data.physics.f_p_div_lower,
                deg_b_plasma_outboard_flux_midplane=self.data.physics.deg_b_plasma_outboard_flux_midplane,
            )

            self.data.divertor.pflux_div_heat_load_mw = (
                wade_divertor.pflux_div_heat_load_mw
            )
            self.data.divertor.dx_div_lower_outboard_strike = (
                wade_divertor.dx_div_lower_outboard_strike
            )
            self.data.divertor.a_div_lower_outboard_wetted = (
                wade_divertor.a_div_lower_outboard_wetted
            )
            self.data.divertor.f_div_lower_flux_expansion = (
                wade_divertor.f_div_lower_flux_expansion
            )
            self.data.divertor.deg_b_div_lower_flux = wade_divertor.deg_b_div_lower_flux
            self.data.divertor.deg_div_lower_outboard_plate_separatrix_poloidal = (
                wade_divertor.deg_div_lower_outboard_plate_separatrix_poloidal
            )
            self.data.divertor.deg_b_div_lower_outboard_grazing = (
                wade_divertor.deg_b_div_lower_outboard_grazing
            )

            if output:
                self.output_wade_divertor_model()

            return

    @property
    def single_divertor_angle(self):
        """
        Calculate the angle subtended by a single divertor.
        Angle is calculated as 180 degrees minus the inboard
        blanket poloidal angle, divided by 2 (for two divertors).
        """
        return (180.0 - self.data.blanket.deg_blkt_inboard_poloidal_plasma) / 2.0

    def divtart(
        self,
        rmajor: float,
        rminor: float,
        triang: float,
        dr_fw_plasma_gap_inboard: float,
        dz_xpoint_divertor: float,
        p_plasma_separatrix_mw: float,
        output: bool,
        i_single_null: int,
        dz_divertor: float,
    ) -> float:
        """Tight aspect ratio tokamak divertor model

        This method calculates the divertor heat load for a tight aspect
        ratio machine, assuming that the power is evenly distributed around the
        divertor chamber by the action of a gaseous target. Each divertor is
        modeled as approximately triangular in the R,Z plane.

        Parameters
        ----------
        rmajor : float
            Plasma major radius (m)
        rminor : float
            Plasma minor radius (m)
        triang : float
            Plasma triangularity
        dr_fw_plasma_gap_inboard : float
            Inboard scrape-off width (m)
        dz_xpoint_divertor : float
            Vertical distance from X-point to divertor (m)
        p_plasma_separatrix_mw : float
            Power to the divertor (MW)
        output : bool
            Indicates whether output should be written to the output file
        i_single_null : int
            1 for single null configuration, 0 for double null
        dz_divertor : float
            Vertical height of the divertor (m)

        Returns
        -------
        float
            Divertor heat load for a tight aspect ratio machine (MW/m2)

        Raises
        ------
        ProcessValueError
            If dz_xpoint_divertor is non-positive

        Notes
        -----
            - This model assumes a tight aspect ratio tokamak with a gaseous target
              divertor. The divertor chamber is modeled as triangular in the R,Z plane,
              and the heat load is calculated based on the total divertor surface area.
            - The method accounts for both single null and double null configurations.

        References
        ----------
            - Y.-K. M. Peng, J. B. Hicks, AEA Fusion, Culham (UK),
            "Engineering feasibility of tight aspect ratio Tokamak (spherical torus)
             reactors".
              1990. https://inis.iaea.org/records/ey2rf-dah04

            - Y.-K. M. Peng, J. B. Hicks,
            “Engineering feasibility of tight aspect ratio tokamak (spherical torus)
            reactors,”
              Osti.gov, 1991. https://www.osti.gov/biblio/1022679
              (accessed Mar. 24, 2025).
        """
        #  Thickness of centrepost + first wall at divertor height

        r1 = (
            rmajor
            - rminor * triang
            - 3.0e0 * dr_fw_plasma_gap_inboard
            + self.data.tfcoil.drtop
        )

        #  Outer radius of divertor region

        r2 = rmajor + rminor

        #  Angle of diagonal divertor plate from horizontal

        if dz_xpoint_divertor <= 0.0e0:
            raise ProcessValueError(
                "Non-positive dz_xpoint_divertor", dz_xpoint_divertor=dz_xpoint_divertor
            )

        theta = math.atan(dz_divertor / (r2 - r1))

        #  Vertical plate area

        a1 = 2.0e0 * np.pi * r1 * dz_divertor

        #  Horizontal plate area

        a2 = np.pi * (r2 * r2 - r1 * r1)

        #  Diagonal plate area

        a3 = a2 / (math.cos(theta) * math.cos(theta))

        #  Total divertor area

        # Single null case
        if i_single_null == DivertorNumberModels.SINGLE_NULL:
            areadv = a1 + a2 + a3
        # Double null case
        elif i_single_null == DivertorNumberModels.DOUBLE_NULL:
            areadv = 2.0 * (a1 + a2 + a3)

        if (
            DivertorHeatLoadModel(self.data.divertor.i_div_heat_load)
            == DivertorHeatLoadModel.PENG_CHAMBER
        ):
            self.data.divertor.pflux_div_heat_load_mw = p_plasma_separatrix_mw / areadv

        if (
            output
            and DivertorHeatLoadModel(self.data.divertor.i_div_heat_load)
            == DivertorHeatLoadModel.PENG_CHAMBER
        ):
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power to the divertor (MW)",
                "(p_plasma_separatrix_mw.)",
                p_plasma_separatrix_mw,
            )
            po.ovarre(self.outfile, "Divertor surface area (m²)", "(areadv)", areadv)
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m²)",
                "(pflux_div_heat_load_mw)",
                self.data.divertor.pflux_div_heat_load_mw,
            )

        elif output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ovarre(
                self.outfile,
                "Power to the divertor (MW)",
                "(p_plasma_separatrix_mw.)",
                p_plasma_separatrix_mw,
            )
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m²)",
                "(pflux_div_heat_load_mw)",
                self.data.divertor.pflux_div_heat_load_mw,
            )
        return self.data.divertor.pflux_div_heat_load_mw

    def divwade(
        self,
        rmajor: float,
        p_plasma_separatrix_mw: float,
        len_plasma_sol_power_decay: float,
        len_plasma_sol_power_spreading: float,
        f_div_flux_expansion: float,
        deg_b_div_lower_outboard_grazing: float,
        rad_fraction_sol: float,
        f_p_div_lower: float,
        deg_b_plasma_outboard_flux_midplane: float,
    ) -> WadeDivertorMetrics:
        """Divertor heat load model (Wade 2020)

        This subroutine calculates the divertor heat flux for any machine,
        with either a single null or double null configuration.
        It uses the Eich scaling (Eich et al. 2013) and spreading factor
        (Scarabosio et al. 2014)
        to calculate the SOL width.
        This is then used with a flux expansion factor to calculate
        the wetted area and then the heat flux.

        Parameters
        ----------
        rmajor : float
            plasma major radius (m)
        p_plasma_separatrix_mw : float
            power to divertor (MW)
        len_plasma_sol_power_decay : float
            SOL power decay length (λ_q) [m]
        len_plasma_sol_power_spreading : float
            SOL power spreading factor (S) [m]
        f_div_flux_expansion : float
            plasma flux expansion in divertor
        deg_b_div_lower_outboard_grazing : float
            field line angle wrt divertor target plate (degrees)
        rad_fraction_sol : float
            SOL radiation fraction
        f_p_div_lower : float
            fraction of power to the lower divertor in double null configuration

        Returns
        -------
        WadeDivertorMetrics


        References
        ----------
        [1] M. R. Wade and J. A. Leuer, “Cost Drivers for a Tokamak-Based Compact Pilot
        Plant,” Fusion Science and Technology, vol. 77, no. 2, pp. 119-143, Feb. 2021,
        doi: 10.1080/15361055.2020.1858670.
        """
        # SOL width at divertor plates (λ_int) [m]
        # λ_int = λ_q + 1.64 * S
        lambda_int = len_plasma_sol_power_decay + 1.64 * len_plasma_sol_power_spreading

        # Flux angle in the divertor as a pitch-like quantity
        alpha_div = f_div_flux_expansion * math.tan(
            math.radians(deg_b_plasma_outboard_flux_midplane)
        )

        if math.isclose(alpha_div, 0.0, abs_tol=1.0e-12):
            raise ProcessValueError(
                "Zero divertor flux angle",
                f_div_flux_expansion=f_div_flux_expansion,
                deg_b_plasma_outboard_flux_midplane=deg_b_plasma_outboard_flux_midplane,
            )

        # Flux angle in the divertor for output [degrees]
        deg_b_div_lower_flux = math.degrees(math.atan(alpha_div))

        # Tilt of the separatrix relative to the target in the poloidal plane
        sin_theta_div = (1.0 + 1.0 / alpha_div**2) * math.sin(
            math.radians(deg_b_div_lower_outboard_grazing)
        )
        sin_theta_div = max(-1.0, min(1.0, sin_theta_div))
        theta_div = math.asin(sin_theta_div)

        # Wetted area
        area_wetted = (
            2.0
            * np.pi
            * rmajor
            * lambda_int
            * f_div_flux_expansion
            * math.sin(theta_div)
        )

        if area_wetted <= 0.0:
            raise ProcessValueError("Non-positive wetted area", area_wetted=area_wetted)

        strike_width = (
            lambda_int * f_div_flux_expansion * math.sin(theta_div)
        )

        # Divertor heat load
        hldiv_base = p_plasma_separatrix_mw * (1.0 - rad_fraction_sol) / area_wetted

        # For double null, calculate heat loads to upper and lower divertors
        # and use the highest
        if self.data.divertor.n_divertors == 2:
            hldiv_lower = f_p_div_lower * hldiv_base
            hldiv_upper = (1.0 - f_p_div_lower) * hldiv_base
            pflux_div_heat_load_mw = max(hldiv_lower, hldiv_upper)
        else:
            pflux_div_heat_load_mw = hldiv_base

        return WadeDivertorMetrics(
            pflux_div_heat_load_mw=pflux_div_heat_load_mw,
            dx_div_lower_outboard_strike=strike_width,
            a_div_lower_outboard_wetted=area_wetted,
            f_div_lower_flux_expansion=f_div_flux_expansion,
            deg_b_div_lower_flux=deg_b_div_lower_flux,
            deg_div_lower_outboard_plate_separatrix_poloidal=math.degrees(theta_div),
            deg_b_div_lower_outboard_grazing=deg_b_div_lower_outboard_grazing,
        )

    def output_wade_divertor_model(self):

        po.oheadr(self.outfile, "Wade Divertor Heat Load Model")
        po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Flux expansion factor (fₓ)",
            "(f_div_flux_expansion)",
            self.data.divertor.f_div_lower_flux_expansion,
        )
        po.ovarre(
            self.outfile,
            "Field line angle wrt to target divertor plate [degrees]",
            "(deg_b_div_lower_outboard_grazing)",
            self.data.divertor.deg_b_div_lower_outboard_grazing,
        )
        po.ovarre(
            self.outfile,
            "Lower divertor flux angle [degrees]",
            "(deg_b_div_lower_flux)",
            self.data.divertor.deg_b_div_lower_flux,
        )
        po.ovarre(
            self.outfile,
            "Lower divertor outboard plate-separatrix poloidal angle [degrees]",
            "(deg_div_lower_outboard_plate_separatrix_poloidal)",
            self.data.divertor.deg_div_lower_outboard_plate_separatrix_poloidal,
        )

        po.ovarre(
            self.outfile,
            "Strike point width (m)",
            "(dx_div_lower_outboard_strike)",
            self.data.divertor.dx_div_lower_outboard_strike,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Divertor heat load [MW/m²]",
            "(pflux_div_heat_load_mw)",
            self.data.divertor.pflux_div_heat_load_mw,
        )

    @staticmethod
    def incident_radiation_power(
        p_plasma_rad_mw: float,
        f_ster_div_single: float,
        n_divertors: int,
    ) -> float:
        """Calculates the total incident radiation power on the divertor box.

        Parameters
        ----------
        p_plasma_rad_mw : float
            Total plasma radiated power in megawatts (MW).
        f_ster_div_single : float
            Fraction of the solid angle subtended by a single divertor.
        n_divertors : int
            Number of divertors.

        Returns
        -------
        float
            Total incident radiation power on the divertor box in megawatts (MW).
        """
        return p_plasma_rad_mw * f_ster_div_single * n_divertors

    @staticmethod
    def incident_neutron_power(
        p_plasma_neutron_mw: float,
        f_ster_div_single: float,
        n_divertors: int,
    ) -> float:
        """Calculates the total incident neutron power on the divertor box.

        Parameters
        ----------
        p_plasma_neutron_mw : float
            Total plasma neutron power in megawatts (MW).
        f_ster_div_single : float
            Fraction of the solid angle subtended by a single divertor.
        n_divertors : int
            Number of divertors.

        Returns
        -------
        float
            Total incident radiation power on the divertor box in megawatts (MW).
        """
        return p_plasma_neutron_mw * f_ster_div_single * n_divertors


class LowerDivertor(Divertor):
    """Module containing lower divertor routines"""

    def run(self, output: bool):
        """Run the LowerDivertor routines"""
        super().run(output=output)

        self.data.divertor.p_div_lower_nuclear_heat_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=self.data.physics.p_plasma_neutron_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=1,
        )

        self.data.divertor.p_div_lower_rad_mw = self.incident_radiation_power(
            p_plasma_rad_mw=self.data.physics.p_plasma_rad_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=1,
        )


class UpperDivertor(Divertor):
    """Module containing upper divertor routines"""

    def run(self, output: bool):
        """Run the UpperDivertor routine"""
        super().run(output=output)

        self.data.divertor.p_div_upper_nuclear_heat_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=self.data.physics.p_plasma_neutron_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=1,
        )

        self.data.divertor.p_div_upper_rad_mw = self.incident_radiation_power(
            p_plasma_rad_mw=self.data.physics.p_plasma_rad_mw,
            f_ster_div_single=self.data.fwbs.f_ster_div_single,
            n_divertors=1,
        )
