"""Module for calculating plasma scrape off layers physics"""

import logging

from process.core import constants
from process.core import process_output as po
from process.core.model import Model

logger = logging.getLogger(__name__)


class ScrapeOffLayer(Model):
    """Model for calculating plasma scrape off layers physics."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """Calculate the scrape off layer physics and update the physics variables."""
        self.data.physics.len_plasma_sol_eich13_power_decay = self.calculate_eich2013_sol_power_decay_length(  # noqa: E501
            p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
            rmajor=self.data.physics.rmajor,
            b_plasma_surface_poloidal_average=self.data.physics.b_plasma_surface_poloidal_average,
            aspect=self.data.physics.aspect,
        )

        self.data.physics.len_plasma_sol_mast14_power_decay_1 = self.calculate_mast2014_sol_power_decay_length_1(  # noqa: E501
            p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
            b_plasma_surface_poloidal_average=self.data.physics.b_plasma_surface_poloidal_average,
        )

        self.data.physics.len_plasma_sol_mast14_power_decay_2 = (
            self.calculate_mast2014_sol_power_decay_length_2(
                p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
                cur_plasma_ma=self.data.physics.plasma_current / 1e6,  # Convert A to MA
            )
        )

    def output(self) -> None:
        """Output plasma scrape off layer physics information."""
        po.oheadr(self.outfile, "Plasma Scrape Off Layer")

        po.osubhd(self.outfile, "Power Decay Lengths (λ_q)")

        po.ovarre(
            self.outfile,
            "Eich 2013 SOL power decay length (λ_q) [m]",
            "(len_plasma_sol_eich13_power_decay)",
            self.data.physics.len_plasma_sol_eich13_power_decay,
        )
        po.ovarre(
            self.outfile,
            "MAST 2014 SOL power decay length 1 (λ_q) [m]",
            "(len_plasma_sol_mast14_power_decay_1)",
            self.data.physics.len_plasma_sol_mast14_power_decay_1,
        )
        po.ovarre(
            self.outfile,
            "MAST 2014 SOL power decay length 2 (λ_q) [m]",
            "(len_plasma_sol_mast14_power_decay_2)",
            self.data.physics.len_plasma_sol_mast14_power_decay_2,
        )

    @staticmethod
    def calculate_eich2013_sol_power_decay_length(
        p_plasma_separatrix_mw: float,
        rmajor: float,
        b_plasma_surface_poloidal_average: float,
        aspect: float,
    ) -> float:
        """Calculate the Eich 2013 SOL power decay length (λ_q).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) (MW)
        rmajor : float
            Major radius of the plasma (R₀) [m]
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (⟨Bₚₒₗ(a)⟩)  [T]
        aspect : float
            Aspect ratio of the plasma  (A)

        Returns
        -------
        float
            Eich 2013 SOL power decay length (λ_q) [m]

        Notes
        -----
        - The paper states that the poloidal field terms is for the outer midplane
        Bₚₒₗ(a), we are using the outer surface average

        References
        ----------
        [1] T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode
        power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9,
        p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

        """
        return (
            1.35e-3
            * p_plasma_separatrix_mw**-0.02
            * rmajor**0.04
            * b_plasma_surface_poloidal_average**-0.92
            * aspect**-0.42
        )

    @staticmethod
    def calculate_mast2014_sol_power_decay_length_1(
        p_plasma_separatrix_mw: float,
        b_plasma_surface_poloidal_average: float,
    ) -> float:
        """Calculate the MAST 2014 SOL power decay length (λ_q).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) (MW)
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (⟨Bₚₒₗ(a)⟩)  [T]

        Returns
        -------
        float
            MAST 2014 SOL power decay length (λ_q) [m]

        Notes
        -----
        - The paper states that the poloidal field terms is for the outer midplane
        Bₚₒₗ(a), we are using the outer surface average

        References
        ----------
        [1] A. J. Thornton and A. Kirk, “Scaling of the scrape-off layer width during
        inter-ELM H modes on MAST as measured by infrared thermography,”
        Plasma Physics and Controlled Fusion, vol. 56, no. 5, p. 055008, Apr. 2014,
        doi: 10.1088/0741-3335/56/5/055008.

        """
        return (
            1.84e-3
            * p_plasma_separatrix_mw**0.18
            * b_plasma_surface_poloidal_average**-0.68
        )

    @staticmethod
    def calculate_mast2014_sol_power_decay_length_2(
        p_plasma_separatrix_mw: float,
        cur_plasma_ma: float,
    ) -> float:
        """Calculate the MAST 2014 SOL power decay length (λ_q).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) (MW)
        cur_plasma_ma : float
            Plasma current (Iₚ) [MA]

        Returns
        -------
        float
            MAST 2014 SOL power decay length (λ_q) [m]


        References
        ----------
        [1] A. J. Thornton and A. Kirk, “Scaling of the scrape-off layer width during
        inter-ELM H modes on MAST as measured by infrared thermography,”
        Plasma Physics and Controlled Fusion, vol. 56, no. 5, p. 055008, Apr. 2014,
        doi: 10.1088/0741-3335/56/5/055008.

        """
        return 4.57e-3 * p_plasma_separatrix_mw**0.22 * cur_plasma_ma**-0.64
