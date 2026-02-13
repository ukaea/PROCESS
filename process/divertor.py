import math

import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import build_variables as bv
from process.data_structure import divertor_variables as dv
from process.data_structure import fwbs_variables as fwbs
from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.exceptions import ProcessValueError


class Divertor:
    """Module containing divertor routines

    This module contains routines relevant for calculating the
    divertor parameters for a fusion power plant.
    """

    def __init__(self):
        self.outfile = constants.NOUT  # output file unit

    def run(self, output: bool):
        """Routine to call the divertor model

        This subroutine calls the divertor routine. This routine scales
        dimensions, powers and field levels which are used as input to
        the Harrison divertor model.

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """

        fwbs.p_div_nuclear_heat_total_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=pv.p_plasma_neutron_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=dv.n_divertors,
        )

        fwbs.p_div_rad_total_mw = self.incident_radiation_power(
            p_plasma_rad_mw=pv.p_plasma_rad_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=dv.n_divertors,
        )

        if dv.i_div_heat_load == 0 and output:
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m2)",
                "(pflux_div_heat_load_mw)",
                dv.pflux_div_heat_load_mw,
            )
            return
        if dv.i_div_heat_load == 1:
            self.divtart(
                pv.rmajor,
                pv.rminor,
                pv.triang,
                bv.dr_fw_plasma_gap_inboard,
                bv.dz_xpoint_divertor,
                pv.p_plasma_separatrix_mw,
                output=output,
                i_single_null=pv.i_single_null,
                dz_divertor=dv.dz_divertor,
            )
            return
        if dv.i_div_heat_load == 2:
            self.divwade(
                pv.rmajor,
                pv.rminor,
                pv.aspect,
                pv.b_plasma_toroidal_on_axis,
                pv.b_plasma_poloidal_average,
                pv.p_plasma_separatrix_mw,
                dv.f_div_flux_expansion,
                pv.nd_plasma_separatrix_electron,
                dv.deg_div_field_plate,
                pv.rad_fraction_sol,
                pv.f_p_div_lower,
                output=output,
            )
            return

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


        Notes
        -----
            - This model assumes a tight aspect ratio tokamak with a gaseous target
              divertor. The divertor chamber is modeled as triangular in the R,Z plane,
              and the heat load is calculated based on the total divertor surface area.
            - The method accounts for both single null and double null configurations.

        :references:
            - Y.-K. M. Peng, J. B. Hicks, AEA Fusion, Culham (UK), "Engineering feasibility of tight aspect ratio Tokamak (spherical torus) reactors".
              1990. https://inis.iaea.org/records/ey2rf-dah04

            - Y.-K. M. Peng, J. B. Hicks, “Engineering feasibility of tight aspect ratio tokamak (spherical torus) reactors,”
              Osti.gov, 1991. https://www.osti.gov/biblio/1022679 (accessed Mar. 24, 2025).
        """

        #  Thickness of centrepost + first wall at divertor height

        r1 = rmajor - rminor * triang - 3.0e0 * dr_fw_plasma_gap_inboard + tfv.drtop

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
        if i_single_null == 1:
            areadv = a1 + a2 + a3
        # Double null case
        elif i_single_null == 0:
            areadv = 2.0 * (a1 + a2 + a3)

        if dv.i_div_heat_load == 1:
            dv.pflux_div_heat_load_mw = p_plasma_separatrix_mw / areadv

        if output and dv.i_div_heat_load == 1:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power to the divertor (MW)",
                "(p_plasma_separatrix_mw.)",
                p_plasma_separatrix_mw,
            )
            po.ovarre(self.outfile, "Divertor surface area (m2)", "(areadv)", areadv)
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m2)",
                "(pflux_div_heat_load_mw)",
                dv.pflux_div_heat_load_mw,
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
                "Divertor heat load (MW/m2)",
                "(pflux_div_heat_load_mw)",
                dv.pflux_div_heat_load_mw,
            )
        return dv.pflux_div_heat_load_mw

    def divwade(
        self,
        rmajor: float,
        rminor: float,
        aspect: float,
        b_plasma_toroidal_on_axis: float,
        b_plasma_poloidal_average: float,
        p_plasma_separatrix_mw: float,
        f_div_flux_expansion: float,
        nd_plasma_separatrix_electron: float,
        deg_div_field_plate: float,
        rad_fraction_sol: float,
        f_p_div_lower: float,
        output: bool,
    ) -> float:
        """Divertor heat load model (Wade 2020)

        This subroutine calculates the divertor heat flux for any machine,
        with either a single null or double null configuration.
        It uses the Eich scaling (Eich et al. 2013) and spreading factor (Scarabosio et al. 2014)
        to calculate the SOL width. This is then used with a flux expansion factor to calculate
        the wetted area and then the heat flux.

        Parameters
        ----------
        rmajor : float
            plasma major radius (m)
        rminor : float
            plasma minor radius (m)
        aspect : float
            tokamak aspect ratio
        b_plasma_toroidal_on_axis : float
            toroidal field (T)
        b_plasma_poloidal_average : float
            poloidal field (T)
        p_plasma_separatrix_mw : float
            power to divertor (MW)
        f_div_flux_expansion : float
            plasma flux expansion in divertor
        nd_plasma_separatrix_electron : float
            electron density at separatrix (m-3)
        deg_div_field_plate : float
            field line angle wrt divertor target plate (degrees)
        rad_fraction_sol : float
            SOL radiation fraction
        f_p_div_lower : float
            fraction of power to the lower divertor in double null configuration

        Returns
        -------
        float
            divertor heat load for a tight aspect ratio machine
        """

        # Radius on midplane
        r_omp = rmajor + rminor

        # B fields on midplane
        Bp_omp = -b_plasma_poloidal_average * rmajor / r_omp

        Bt_omp = -b_plasma_toroidal_on_axis * rmajor / r_omp

        # Eich scaling for lambda_q
        lambda_eich = (
            1.35
            * p_plasma_separatrix_mw**-0.02
            * rmajor**0.04
            * b_plasma_poloidal_average**-0.92
            * aspect**0.42
        )

        # Spreading factor
        spread_fact = (
            0.12
            * (nd_plasma_separatrix_electron / 1e19) ** -0.02
            * p_plasma_separatrix_mw**-0.21
            * rmajor**0.71
            * b_plasma_poloidal_average**-0.82
        )

        # SOL width
        lambda_int = lambda_eich + 1.64 * spread_fact

        # Flux angle on midplane
        alpha_mid = math.degrees(math.atan(Bp_omp / Bt_omp))

        # Flux angle in the divertor
        alpha_div = f_div_flux_expansion * alpha_mid

        # Tilt of the separatrix relative to the target in the poloidal plane
        theta_div = math.asin(
            (1 + 1 / alpha_div**2) * math.sin(math.radians(deg_div_field_plate))
        )

        # Wetted area
        area_wetted = (
            2 * np.pi * rmajor * lambda_int * f_div_flux_expansion * math.sin(theta_div)
        )

        # Divertor heat load
        hldiv_base = p_plasma_separatrix_mw * (1 - rad_fraction_sol) / area_wetted

        # For double null, calculate heat loads to upper and lower divertors and use the highest
        if dv.n_divertors == 2:
            hldiv_lower = f_p_div_lower * hldiv_base
            hldiv_upper = (1.0 - f_p_div_lower) * hldiv_base
            dv.pflux_div_heat_load_mw = max(hldiv_lower, hldiv_upper)
        else:
            dv.pflux_div_heat_load_mw = hldiv_base

        if output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Flux expansion",
                "(f_div_flux_expansion)",
                f_div_flux_expansion,
            )
            po.ovarre(
                self.outfile,
                "Field line angle wrt to target divertor plate (degrees)",
                "(deg_div_field_plate)",
                deg_div_field_plate,
            )
            po.ovarre(
                self.outfile,
                "Divertor heat load (MW/m2)",
                "(pflux_div_heat_load_mw)",
                dv.pflux_div_heat_load_mw,
            )
        return dv.pflux_div_heat_load_mw

    def incident_radiation_power(
        self,
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

    def incident_neutron_power(
        self,
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
        super().run(output=output)

        dv.p_div_lower_nuclear_heat_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=pv.p_plasma_neutron_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=1,
        )

        dv.p_div_lower_rad_mw = self.incident_radiation_power(
            p_plasma_rad_mw=pv.p_plasma_rad_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=1,
        )


class UpperDivertor(Divertor):
    """Module containing upper divertor routines"""

    def run(self, output: bool):
        super().run(output=output)

        dv.p_div_upper_nuclear_heat_mw = self.incident_neutron_power(
            p_plasma_neutron_mw=pv.p_plasma_neutron_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=1,
        )

        dv.p_div_upper_rad_mw = self.incident_radiation_power(
            p_plasma_rad_mw=pv.p_plasma_rad_mw,
            f_ster_div_single=fwbs.f_ster_div_single,
            n_divertors=1,
        )
