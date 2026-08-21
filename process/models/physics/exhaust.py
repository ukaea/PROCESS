"""Module for plasma exhaust calculations and analysis."""

import logging
from dataclasses import dataclass

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure.physics_variables import DivertorNumberModels

logger = logging.getLogger(__name__)


@dataclass
class DivertorSeparatrixPowerSplits:
    """Dataclass to hold the power splits to the divertor targets."""

    f_p_div_inboard_separatrix: float = 0.0
    """Fraction of total separatrix power to inboard divertor targets"""

    f_p_div_outboard_separatrix: float = 0.0
    """Fraction of total separatrix power to outboard divertor targets"""

    f_p_div_inboard_lower_separatrix: float = 0.0
    """Fraction of total separatrix power to inboard lower divertor target"""

    f_p_div_inboard_upper_separatrix: float = 0.0
    """Fraction of total separatrix power to inboard upper divertor target"""

    f_p_div_outboard_lower_separatrix: float = 0.0
    """Fraction of total separatrix power to outboard lower divertor target"""

    f_p_div_outboard_upper_separatrix: float = 0.0
    """Fraction of total separatrix power to outboard upper divertor target"""


class PlasmaExhaust(Model):
    """Class to hold plasma exhaust calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """PlasmaExhaust model isn't run."""

    def output(self):
        """Output plasma exhaust results to the output file."""
        po.oheadr(self.outfile, "Plasma Exhaust")
        po.ovarre(
            self.outfile,
            "Plasma separatrix power (Pₛₑₚ) (MW)",
            "(p_plasma_separatrix_mw)",
            self.data.physics.p_plasma_separatrix_mw,
            "OP ",
        )

        if self.data.physics.p_plasma_separatrix_mw <= 0.001e0:
            logger.error(
                "Possible problem with high radiation power, forcing "
                "p_plasma_separatrix_mw to odd values. "
                f"{self.data.physics.p_plasma_separatrix_mw=}"
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

        if self.data.divertor.n_divertors == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Plasma separatrix power over major radius (Pₛₑₚ / R₀) [MW/m] "
                "(On peak divertor)",
                "(p_plasma_separatrix_rmajor_mw)",
                self.data.physics.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "EU-DEMO divertor protection re-attachment metric (PₛₑₚBₜ / q₉₅AR₀) "
                "[MWT/m] (On peak divertor)",
                "(p_div_bt_q_aspect_rmajor_mw)",
                self.data.physics.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Plasma separatrix power over major radius (Pₛₑₚ / R₀) [MW/m]",
                "(p_plasma_separatrix_rmajor_mw)",
                self.data.physics.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "EU-DEMO divertor protection re-attachment metric (PₛₑₚBₜ / q₉₅AR₀) "
                "[MWT/m]",
                "(p_div_bt_q_aspect_rmajor_mw)",
                self.data.physics.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        self.output_brunner_divertor_power_splits()

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
        Calculate the power crossing the separatrix (Pₛₑₚ).

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
        Calculate the power crossing the separatrix per unit major radius (Pₛₑₚ / R₀).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) [MW].
        rmajor : float
            Plasma major radius (R₀) [m].

        Returns
        -------
        float
            Power crossing the separatrix per unit major radius (Pₛₑₚ / R₀) [MW/m].
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
        """Calculate the EU DEMO divertor protection re-attachment metric for plasma
        exhaust (PₛₑₚBₜ / q₉₅AR₀).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) [MW].
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field on axis (Bₜ) [T].
        q95 : float
            Safety factor at 95% flux surface (q₉₅).
        aspect : float
            Aspect ratio of the plasma (A).
        rmajor : float
            Plasma major radius (R₀) [m].

        Returns
        -------
        float
            EU DEMO re-attachment metric (PₛₑₚBₜ / q₉₅AR₀) [MW T /m].

        References
        ----------
        [1] M. Siccinio, G. Federici, R. Kembleton, H. Lux, F. Maviglia, and J. Morris,
        "Figure of merit for divertor protection in the preliminary design of the
        EU-DEMO reactor," Nuclear Fusion, vol. 59, no. 10, pp. 106026-106026,
        Jul. 2019, doi: https://doi.org/10.1088/1741-4326/ab3153.

        [2] H. Zohm et al.,
        "A stepladder approach to a tokamak fusion power plant,"
        Nuclear Fusion, vol. 57, no. 8, pp. 086002-086002, May 2017,
        doi: https://doi.org/10.1088/1741-4326/aa739e.
        """
        return (p_plasma_separatrix_mw * b_plasma_toroidal_on_axis) / (
            q95 * aspect * rmajor
        )

    @staticmethod
    def calculate_radiation_fraction(
        p_plasma_rad_mw: float, p_plasma_heating_mw: float
    ) -> float:
        """
        Calculate the radiation fraction of the plasma.

        Parameters
        ----------
        p_plasma_rad_mw : float
            Radiated power from plasma (MW).
        p_plasma_heating_mw : float
            Total plasma heating power (MW).

        Returns
        -------
        float
            Radiation fraction of the plasma.
        """
        if p_plasma_heating_mw == 0:
            logger.warning(
                "Total plasma heating power is zero, "
                "cannot calculate radiation fraction."
            )
            return 0.0

        return p_plasma_rad_mw / p_plasma_heating_mw

    def output_brunner_divertor_power_splits(self):
        """Output the Brunner divertor power splits to the output file."""
        if self.data.stellarator.istell == 0:
            po.osubhd(self.outfile, "Brunner Divertor Power Splits:")

            for op in [
                (
                    (
                        "Requested fraction of power to the lower divertor in double "
                        "null configuration"
                    ),
                    "(f_p_div_lower_separatrix)",
                    self.data.physics.f_p_div_lower_separatrix,
                ),
                (
                    (
                        "Required distance between the first and second plasma "
                        "separatrixes at the outer midplane (δR_sep) [m]"
                    ),
                    "(dr_plasma_outboard_midplane_separatrix_separation)",
                    self.data.physics.dr_plasma_outboard_midplane_separatrix_separation,
                ),
                None,
                (
                    "Outboard side heat flux decay length (m)",
                    "(len_sol_outboard_power_decay)",
                    self.data.physics.len_sol_outboard_power_decay,
                ),
                None,
                (
                    "Fraction of separatrix power on the inner target(s)",
                    "(f_p_div_inboard_separatrix)",
                    self.data.physics.f_p_div_inboard_separatrix,
                ),
                (
                    "Fraction of separatrix power on the outer target(s)",
                    "(f_p_div_outboard_separatrix)",
                    self.data.physics.f_p_div_outboard_separatrix,
                ),
                None,
                (
                    "Fraction of separatrix power on the inner lower target",
                    "(f_p_div_lower_inboard_separatrix)",
                    self.data.physics.f_p_div_lower_inboard_separatrix,
                ),
                (
                    "Separatrix power on the inner lower target",
                    "(p_div_lower_inboard_separatrix_mw)",
                    self.data.physics.p_div_lower_inboard_separatrix_mw,
                ),
                None,
                (
                    "Fraction of separatrix power on the outer lower target",
                    "(f_p_div_lower_outboard_separatrix)",
                    self.data.physics.f_p_div_lower_outboard_separatrix,
                ),
                (
                    "Separatrix power on the outer lower target",
                    "(p_div_lower_outboard_separatrix_mw)",
                    self.data.physics.p_div_lower_outboard_separatrix_mw,
                ),
            ]:
                if op is None:
                    po.oblnkl(self.outfile)
                else:
                    desc, var, val = op
                    po.ovarre(self.outfile, desc, var, val, "OP ")

            if (
                DivertorNumberModels(self.data.physics.i_single_null)
                == DivertorNumberModels.DOUBLE_NULL
            ):
                po.oblnkl(self.outfile)
                po.ocmmnt(self.outfile, "----------------------------")
                po.oblnkl(self.outfile)

                for op in [
                    (
                        "Fraction of separatrix power on the inner upper target",
                        "(f_p_div_upper_inboard_separatrix)",
                        self.data.physics.f_p_div_upper_inboard_separatrix,
                    ),
                    (
                        "Separatrix power on the inner upper target",
                        "(p_div_upper_inboard_separatrix_mw)",
                        self.data.physics.p_div_upper_inboard_separatrix_mw,
                    ),
                    None,
                    (
                        "Fraction of separatrix power on the outer upper target",
                        "(f_p_div_upper_outboard_separatrix)",
                        self.data.physics.f_p_div_upper_outboard_separatrix,
                    ),
                    (
                        "Separatrix power on the outer upper target",
                        "(p_div_upper_outboard_separatrix_mw)",
                        self.data.physics.p_div_upper_outboard_separatrix_mw,
                    ),
                ]:
                    if op is None:
                        po.oblnkl(self.outfile)
                    else:
                        desc, var, val = op
                        po.ovarre(self.outfile, desc, var, val, "OP ")


def calculate_brunner_divertor_power_splits(
    dr_outboard_midplane_sep: float,
    len_plasma_sol_outboard_power_decay: float,
    len_plasma_sol_inboard_power_decay: float,
) -> DivertorSeparatrixPowerSplits:
    """
    Calculate the power splits to the divertor targets using Brunner's method.

    Parameters
    ----------
    dr_outboard_midplane_sep : float
        Radial separation of the plasma outboard midplane separatrix (δR_sep) [m].
    len_plasma_sol_outboard_power_decay : float
        Power decay length in the scrape-off layer (λ_q) [m].
    len_plasma_sol_inboard_power_decay : float
        Power decay length in the scrape-off layer for the inboard side (λᵢₙ_q) [m].

    Returns
    -------
    DivertorSeparatrixPowerSplits
        Dataclass containing the fractions of total separatrix power to each
        divertor target.

    Notes
    -----
    - The fitted value for `f_p_inner_sep_0` and `f_p_inner_sep_infinity` are taken
    from the fit given by Petrie et al. [2]

    References
    ----------
    [1] D. Brunner, A. Q. Kuang, B. LaBombard, and J. L. Terry, “The dependence of
    divertor power sharing on magnetic flux balance in near double-null
    configurations on Alcator C-Mod,” Nuclear Fusion, vol. 58, no. 7, p. 076010,
    May 2018, doi: https://doi.org/10.1088/1741-4326/aac006.

    [2] T. W. Petrie et al., “The effect of divertor magnetic balance on H-mode
    performance in DIII-D,” Journal of Nuclear Materials, vol. 290-293, pp. 935-939,
    Mar. 2001, doi: https://doi.org/10.1016/S0022-3115(00)00492-X
    """
    # Fraction of the power to the inner divertors at δR_sep = 0 and δR_sep → ∞
    f_p_inner_sep_0 = 0.16e0
    f_p_inner_sep_infinity = 0.41e0

    # Fractions of total outboard power going to each target
    # Outboard lower divertor
    f_p_outboard_lower = 1 / (
        1 + np.exp(dr_outboard_midplane_sep / len_plasma_sol_outboard_power_decay)
    )

    # Outboard upper divertor
    f_p_outboard_upper = 1 / (
        1 + np.exp(-dr_outboard_midplane_sep / len_plasma_sol_outboard_power_decay)
    )

    # Fractions of the total inboard power going to each target
    f_p_inboard_lower = 1 / (
        1 + np.exp(dr_outboard_midplane_sep / (len_plasma_sol_inboard_power_decay))
    )

    f_p_total_inboard = f_p_inner_sep_0 + (f_p_inner_sep_0 - f_p_inner_sep_infinity) * (
        1.0e0
        - (
            2.0e0
            / (
                1.0e0
                + np.exp(
                    -(
                        (dr_outboard_midplane_sep / len_plasma_sol_outboard_power_decay)
                        ** 2
                    )
                )
            )
        )
    )

    f_p_total_outboard = 1.0e0 - f_p_total_inboard

    f_p_inboard_upper = f_p_total_inboard * (1.0e0 - f_p_inboard_lower)
    f_p_outboard_upper = f_p_total_outboard * (1.0e0 - f_p_outboard_lower)
    f_p_inboard_lower = f_p_total_inboard * f_p_inboard_lower
    f_p_outboard_lower = f_p_total_outboard * f_p_outboard_lower

    return DivertorSeparatrixPowerSplits(
        f_p_div_inboard_separatrix=f_p_total_inboard,
        f_p_div_outboard_separatrix=f_p_total_outboard,
        f_p_div_inboard_lower_separatrix=f_p_inboard_lower,
        f_p_div_inboard_upper_separatrix=f_p_inboard_upper,
        f_p_div_outboard_lower_separatrix=f_p_outboard_lower,
        f_p_div_outboard_upper_separatrix=f_p_outboard_upper,
    )
