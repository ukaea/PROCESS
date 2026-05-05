import logging

from process.core import constants
from process.core import process_output as po
from process.data_structure import (
    divertor_variables,
    physics_variables,
)

logger = logging.getLogger(__name__)


class PlasmaExhaust:
    """Class to hold plasma exhaust calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def output(self):
        """Output plasma exhaust results to the output file."""
        po.oheadr(self.outfile, "Plasma Exhaust")
        po.ovarre(
            self.outfile,
            "Plasma separatrix power (Pₛₑₚ) (MW)",
            "(p_plasma_separatrix_mw)",
            physics_variables.p_plasma_separatrix_mw,
            "OP ",
        )

        if physics_variables.p_plasma_separatrix_mw <= 0.001e0:
            logger.error(
                "Possible problem with high radiation power, forcing p_plasma_separatrix_mw to odd values. "
                f"{physics_variables.p_plasma_separatrix_mw=}"
            )
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(self.outfile, "          Power into divertor zone is unrealistic;")
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if divertor_variables.n_divertors == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Plasma separatrix power over major radius (Pₛₑₚ / R₀) (MW/m) (On peak divertor)",
                "(p_plasma_separatrix_rmajor_mw)",
                physics_variables.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "EU-DEMO divertor protection re-attachment metric (PₛₑₚBₜ / q₉₅AR₀) (MWT/m) (On peak divertor)",
                "(p_div_bt_q_aspect_rmajor_mw)",
                physics_variables.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Plasma separatrix power over major radius (Pₛₑₚ / R₀) (MW/m)",
                "(p_plasma_separatrix_rmajor_mw)",
                physics_variables.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "EU-DEMO divertor protection re-attachment metric (PₛₑₚBₜ / q₉₅AR₀) (MWT/m)",
                "(p_div_bt_q_aspect_rmajor_mw)",
                physics_variables.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )

        po.oblnkl(self.outfile)

    @staticmethod
    def calculate_separatrix_power(
        f_p_alpha_plasma_deposited: float,
        p_alpha_total_mw: float,
        p_non_alpha_charged_mw: float,
        p_hcd_injected_total_mw: float,
        p_plasma_ohmic_mw: float,
        p_plasma_rad_mw: float,
    ) -> float:
        """
        Calculate the power crossing the separatrix (P_sep).

        Parameters
        ----------
        f_p_alpha_plasma_deposited : float
            Fraction of alpha power deposited in plasma.
        p_alpha_total_mw : float
            Total alpha power produced (MW).
        p_non_alpha_charged_mw : float
            Power from non-alpha charged particles (MW).
        p_hcd_injected_total_mw : float
            Total power injected by heating and current drive (MW).
        p_plasma_ohmic_mw : float
            Ohmic heating power (MW).
        p_plasma_rad_mw : float
            Radiated power from plasma (MW).

        Returns
        -------
        float
            Power crossing the separatrix (MW).
        """
        return (
            f_p_alpha_plasma_deposited * p_alpha_total_mw
            + p_non_alpha_charged_mw
            + p_hcd_injected_total_mw
            + p_plasma_ohmic_mw
            - p_plasma_rad_mw
        )

    @staticmethod
    def calculate_psep_over_r_metric(
        p_plasma_separatrix_mw: float, rmajor: float
    ) -> float:
        """
        Calculate the power crossing the separatrix per unit major radius (P_sep/R).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (MW).
        rmajor : float
            Plasma major radius (m).

        Returns
        -------
        float
            Power crossing the separatrix per unit major radius (MW/m).
        """
        return p_plasma_separatrix_mw / rmajor

    @staticmethod
    def calculate_eu_demo_re_attachment_metric(
        p_plasma_separatrix_mw: float,
        b_plasma_toroidal_on_axis: float,
        q95: float,
        aspect: float,
        rmajor: float,
    ) -> float:
        """Calculate the EU DEMO divertor protection re-attachment metric for plasma exhaust.

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (MW).
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field on axis (T).
        q95 : float
            Safety factor at 95% flux surface.
        aspect : float
            Aspect ratio of the plasma.
        rmajor : float
            Plasma major radius (m).

        Returns
        -------
        float
            EU DEMO re-attachment metric (MW T /m).

        References
        ----------
        - M. Siccinio, G. Federici, R. Kembleton, H. Lux, F. Maviglia, and J. Morris,
          "Figure of merit for divertor protection in the preliminary design of the EU-DEMO reactor,"
          Nuclear Fusion, vol. 59, no. 10, pp. 106026-106026, Jul. 2019,
          doi: https://doi.org/10.1088/1741-4326/ab3153.

        - H. Zohm et al.,
          "A stepladder approach to a tokamak fusion power plant,"
          Nuclear Fusion, vol. 57, no. 8, pp. 086002-086002, May 2017,
          doi: https://doi.org/10.1088/1741-4326/aa739e.
        """
        return (p_plasma_separatrix_mw * b_plasma_toroidal_on_axis) / (
            q95 * aspect * rmajor
        )
