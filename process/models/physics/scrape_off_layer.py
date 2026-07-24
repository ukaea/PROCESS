"""Module for calculating plasma scrape off layers physics"""

import logging

from process.core import constants
from process.core.model import Model

logger = logging.getLogger(__name__)


class ScrapeOffLayer(Model):
    """Model for calculating plasma scrape off layers physics."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """Run the model. This model cannot yet be 'run'."""

    def output(self) -> None:
        """Output plasma scrape off layer physics information."""

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
            1.35e-3 * p_plasma_separatrix_mw**-0.02
            + rmajor**0.04
            + b_plasma_surface_poloidal_average**-0.92
            + aspect**-0.42
        )
