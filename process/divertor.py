import math

from process import process_output as po
from process.fortran import build_variables as bv
from process.fortran import constants
from process.fortran import divertor_variables as dv
from process.fortran import error_handling as eh
from process.fortran import physics_variables as pv
from process.fortran import tfcoil_variables as tfv


class Divertor:
    """Module containing divertor routines
    author: P J Knight, CCFE, Culham Science Centre

    This module contains routines relevant for calculating the
    divertor parameters for a fusion power plant.
    """

    def __init__(self) -> None:
        self.outfile = constants.nout  # output file unit

    def run(self, output: bool) -> None:
        """Routine to call the divertor model
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine calls the divertor routine. This routine scales
        dimensions, powers and field levels which are used as input to
        the Harrison divertor model.

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        if dv.i_hldiv == 0 and output:
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)
            return
        if dv.i_hldiv == 1:
            self.divtart(
                pv.rmajor,
                pv.rminor,
                pv.triang,
                bv.dr_fw_plasma_gap_inboard,
                bv.dz_xpoint_divertor,
                pv.pdivt,
                output=output,
                i_single_null=pv.i_single_null,
            )
            return
        if dv.i_hldiv == 2:
            self.divwade(
                pv.rmajor,
                pv.rminor,
                pv.aspect,
                pv.bt,
                pv.bp,
                pv.pdivt,
                dv.flux_exp,
                pv.nesep,
                dv.beta_div,
                pv.rad_fraction_sol,
                pv.ftar,
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
        pdivt: float,
        output: bool,
        i_single_null: int,
    ) -> float:
        """Tight aspect ratio tokamak divertor model
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine calculates the divertor heat load for a tight aspect
        ratio machine, by assuming that the power is evenly spread around the
        divertor chamber by the action of a gaseous target. Each divertor is
        assumed to be approximately triangular in the R,Z plane.

        :param rmajor: plasma major radius (m)
        :type rmajor: float

        :param rminor: plasma minor radius (m)
        :type rminor: float

        :param triang: plasma triangularity
        :type triang: float

        :param dr_fw_plasma_gap_inboard: inboard scrape-off width (m)
        :type dr_fw_plasma_gap_inboard: float

        :param dz_xpoint_divertor: vertical distance from X-point to divertor (m)
        :type dz_xpoint_divertor: float

        :param pdivt: power to the divertor (MW)
        :type pdivt: float

        :param output: indicate whether output should be written to the output file, or not
        :type output: bool

        :param i_single_null: 1 for single null configuration, 0 for double null
        :type i_single_null: int

        :returns: divertor heat load for a tight aspect ratio machine (MW/m2)
        :rtype: float

        :notes:
            This model assumes a tight aspect ratio tokamak with a gaseous target
            divertor. The divertor chamber is modeled as triangular in the R,Z plane,
            and the heat load is calculated based on the total divertor surface area.

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
            eh.fdiags[0] = dz_xpoint_divertor
            eh.report_error(22)

        theta = math.atan(dz_xpoint_divertor / (r2 - r1))

        #  Vertical plate area

        a1 = 2.0e0 * constants.pi * r1 * dz_xpoint_divertor

        #  Horizontal plate area

        a2 = constants.pi * (r2 * r2 - r1 * r1)

        #  Diagonal plate area

        a3 = a2 / (math.cos(theta) * math.cos(theta))

        #  Total divertor area

        # Single null case
        if i_single_null == 1:
            areadv = a1 + a2 + a3
        # Double null case
        elif i_single_null == 0:
            areadv = 2.0 * (a1 + a2 + a3)

        if dv.i_hldiv == 1:
            dv.hldiv = pdivt / areadv

        if output and dv.i_hldiv == 1:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(self.outfile, "Power to the divertor (MW)", "(pdivt.)", pdivt)
            po.ovarre(self.outfile, "Divertor surface area (m2)", "(areadv)", areadv)
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)

        elif output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ovarre(self.outfile, "Power to the divertor (MW)", "(pdivt.)", pdivt)
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)
        return dv.hldiv

    def divwade(
        self,
        rmajor: float,
        rminor: float,
        aspect: float,
        bt: float,
        bp: float,
        pdivt: float,
        flux_exp: float,
        nesep: float,
        beta_div: float,
        rad_fraction_sol: float,
        ftar: float,
        output: bool,
    ) -> float:
        """Divertor heat load model (Wade 2020)
        author: J A Foster, CCFE, Culham Science Centre

        This subroutine calculates the divertor heat flux for any machine,
        with either a single null or double null configuration.
        It uses the Eich scaling (Eich et al. 2013) and spreading factor (Scarabosio et al. 2014)
        to calculate the SOL width. This is then used with a flux expansion factor to calculate
        the wetted area and then the heat flux.

        :param rmajor: plasma major radius (m)
        :type rmajor: float

        :param rminor: plasma minor radius (m)
        :type rminor: float

        :param aspect: tokamak aspect ratio
        :type aspect: float

        :param bt: toroidal field (T)
        :type bt: float

        :param bp: poloidal field (T)
        :type bp: float

        :param pdivt: power to divertor (MW)
        :type pdivt: float

        :param flux_exp: plasma flux expansion in divertor
        :type flux_exp: float

        :param nesep: electron density at separatrix (m-3)
        :type nesep: float

        :param beta_div: field line angle wrt divertor target plate (degrees)
        :type beta_div: float

        :param rad_fraction_sol: SOL radiation fraction
        :type rad_fraction_sol: float

        :param ftar: fraction of power to the lower divertor in double null configuration
        :type ftar: float

        :returns: divertor heat load for a tight aspect ratio machine
        :rtype: float

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean

        :notes:

        :references:




        """

        # Radius on midplane
        r_omp = rmajor + rminor

        # B fields on midplane
        Bp_omp = -bp * rmajor / r_omp

        Bt_omp = -bt * rmajor / r_omp

        # Eich scaling for lambda_q
        lambda_eich = 1.35 * pdivt**-0.02 * rmajor**0.04 * bp**-0.92 * aspect**0.42

        # Spreading factor
        spread_fact = (
            0.12 * (nesep / 1e19) ** -0.02 * pdivt**-0.21 * rmajor**0.71 * bp**-0.82
        )

        # SOL width
        lambda_int = lambda_eich + 1.64 * spread_fact

        # Flux angle on midplane
        alpha_mid = math.degrees(math.atan(Bp_omp / Bt_omp))

        # Flux angle in the divertor
        alpha_div = flux_exp * alpha_mid

        # Tilt of the separatrix relative to the target in the poloidal plane
        theta_div = math.asin((1 + 1 / alpha_div**2) * math.sin(math.radians(beta_div)))

        # Wetted area
        area_wetted = (
            2 * constants.pi * rmajor * lambda_int * flux_exp * math.sin(theta_div)
        )

        # Divertor heat load
        hldiv_base = pdivt * (1 - rad_fraction_sol) / area_wetted

        # For double null, calculate heat loads to upper and lower divertors and use the highest
        if pv.idivrt == 2:
            hldiv_lower = ftar * hldiv_base
            hldiv_upper = (1.0 - ftar) * hldiv_base
            dv.hldiv = max(hldiv_lower, hldiv_upper)
        else:
            dv.hldiv = hldiv_base

        if output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(self.outfile, "Flux expansion", "(flux_exp)", flux_exp)
            po.ovarre(
                self.outfile,
                "Field line angle wrt to target divertor plate (degrees)",
                "(beta_div)",
                beta_div,
            )
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)
        return dv.hldiv
