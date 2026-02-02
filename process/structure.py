import logging
import math

import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import build_variables as bv
from process.data_structure import divertor_variables as divv
from process.data_structure import fwbs_variables as fwbsv
from process.data_structure import pfcoil_variables as pfv
from process.data_structure import physics_variables as pv
from process.data_structure import structure_variables as stv
from process.data_structure import tfcoil_variables as tfv

logger = logging.getLogger(__name__)


class Structure:
    """Class containing support structure calculations
    author: P J Knight, CCFE, Culham Science Centre

    This class contains routines for calculating the
    parameters of the support structure for a
    fusion power plant.
    """

    def __init__(self) -> None:
        self.outfile = constants.NOUT  # output file unit

    def run(self, output: bool = False) -> None:
        """Structure calculation caller
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine calls the support structure mass calculations.

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        # Total weight of the PF coil conductor and its structure
        total_weight_pf = pfv.m_pf_coil_conductor_total + pfv.m_pf_coil_structure_total

        (
            stv.fncmass,
            stv.aintmass,
            stv.clgsmass,
            stv.coldmass,
            stv.gsmass,
        ) = self.structure(
            pv.plasma_current,
            pv.rmajor,
            pv.rminor,
            pv.kappa,
            pv.b_plasma_toroidal_on_axis,
            tfv.i_tf_sup,
            pfv.i_pf_conductor,
            bv.dr_tf_inner_bore + bv.dr_tf_outboard + bv.dr_tf_inboard,
            bv.z_tf_inside_half,
            fwbsv.whtshld,
            divv.m_div_plate,
            total_weight_pf,
            tfv.m_tf_coils_total,
            fwbsv.m_fw_total,
            fwbsv.m_blkt_total,
            fwbsv.m_fw_blkt_div_coolant_total,
            fwbsv.dewmkg,
            output=output,
        )

    def structure(
        self,
        ai,
        r0,
        a,
        akappa,
        b0,
        i_tf_sup,
        i_pf_conductor,
        tf_h_width,
        tfhmax,
        shldmass,
        dvrtmass,
        pfmass,
        tfmass,
        m_fw_total,
        blmass,
        m_fw_blkt_div_coolant_total,
        dewmass,
        output,
    ):
        """Method to calculate mass of support structure
        author: P J Knight, CCFE, Culham Science Centre
        author: J Galambos, ORNL
        Reprogrammed by J. Galambos to match the ITER (i.e. B. Spears) rules.

        :param ai: plasma current (max design value) (A)
        :type ai: float

        :param r0: plasma major radius (m)
        :type r0: float

        :param a: plasma minor radius (m)
        :type a: float

        :param akappa: plasma elongation
        :type akappa: float

        :param b0: axial B-field (T)
        :type b0: float

        :param itfsup: switch denoting whether TF coils are superconducting
        :type itfsup: integer

        :param i_pf_conductor: switch denoting whether PF & CS coils are resistive
        :type i_pf_conductor: integer

        :param tf_h_width: TF coil horizontal dr_bore (m)
        :type tf_h_width: float

        :param tfhmax: TF coil max height (m)
        :type tfhmax: float

        :param shldmass: total mass of shield (kg)
        :type shldmass: float

        :param dvrtmass: total mass of divertor and assoc. structure (kg)
        :type dvrtmass: float

        :param pfmass: total mass of PF coils plus cases (kg)
        :type pfmass: float

        :param tfmass: total mass of TF coils plus cases (kg)
        :type tfmass: float

        :param blmass: blanket mass (kg)
        :type blmass: float

        :param m_fw_total: first wall mass (kg)
        :type m_fw_total: float

        :param m_fw_blkt_div_coolant_total: total water coolant mass (kg)
        :type m_fw_blkt_div_coolant_total: float

        :param dewmass: vacuum vessel + cryostat mass (kg)
        :type dewmass: float

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean

        :returns:
            - fncmass (`float`) mass of outer pf coil support fence (kg)
            - aintmass (`float`) mass of intercoil support (kg)
            - clgsmass (`float`) coil gravity support mass (kg)
            - coldmass (`float`) total mass of cryogenic temp. stuff (kg)
            - gsm (`float`) gravity support for magnets, and shield/blanket (kg)
        """

        #  Outer PF coil fence (1990 ITER fit)
        fncmass = 2.1e-11 * ai * ai * r0 * akappa * a

        #  Intercoil support between TF coils to react overturning moment
        #  (scaled to 1990 ITER fit)
        aintmass = 1.4e6 * (ai / 2.2e7) * b0 / 4.85e0 * tf_h_width**2 / 50.0e0
        try:
            assert aintmass < np.inf
        except AssertionError:
            logger.exception("aintmass is inf. Kludging to 1e10.")
            aintmass = 1e10

        #  Total mass of coils plus support plus vacuum vessel + cryostat
        coilmass = tfmass + pfmass + aintmass + dewmass

        #  Total mass of cooled components
        coldmass = 0.0e0
        if i_tf_sup == 1:
            coldmass = coldmass + tfmass + aintmass + dewmass
        if i_pf_conductor != 1:
            coldmass = coldmass + pfmass

        #  Coil gravity support mass
        #  Set density (kg/m3) and allowable stress (Pa)

        dens = 7.8e3
        sigal = 2.5e7
        clgsmass = coilmass * (r0 / 6.0e0) * 9.1e0 * 9.807e0 * dens / sigal

        #  Gravity support masses scaled from Spears algorithms (9/90) :

        #  Torus leg support

        ws1 = m_fw_blkt_div_coolant_total + m_fw_total + blmass + shldmass + dvrtmass
        gsm1 = 5.0e0 * 9.807e0 * ws1 * dens / sigal

        #  Ring beam

        ws2 = ws1 + tfmass + pfmass + aintmass + clgsmass
        gsm2 = 1.0e-3 * 34.77e0 * (r0 + 1.0e0) * math.sqrt(0.001e0 * ws2) * dens

        #  Ring legs (JG: this term may be too big, need to check)

        gsm3 = 1.0e-6 * 0.3e0 * (tfhmax + 2.0e0) * ws2 * dens

        gsm = gsm1 + gsm2 + gsm3

        #  Output section

        if output:
            po.oheadr(self.outfile, "Support Structure")
            po.ovarre(
                self.outfile,
                "Outer PF coil fence mass (kg)",
                "(fncmass)",
                fncmass,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Intercoil support structure mass (kg)",
                "(aintmass)",
                aintmass,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mass of cooled components (kg)",
                "(coldmass)",
                coldmass,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Gravity support structure mass (kg)",
                "(clgsmass)",
                clgsmass,
                "OP ",
            )
            po.ovarre(
                self.outfile, "Torus leg support mass (kg)", "(gsm1)", gsm1, "OP "
            )
            po.ovarre(self.outfile, "Ring beam mass (kg)", "(gsm2)", gsm2, "OP ")
            po.ovarre(self.outfile, "Ring legs mass (kg)", "(gsm3)", gsm3, "OP ")

        return fncmass, aintmass, clgsmass, coldmass, gsm
