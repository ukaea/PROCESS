import logging

import numpy as np

from process.core import constants
from process.data_structure import (
    physics_variables,
)
from process.models.physics.plasma_current import PlasmaCurrent

logger = logging.getLogger(__name__)


class PlasmaFields:
    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.current = PlasmaCurrent()

    def calculate_surface_averaged_poloidal_field(
        self,
        i_plasma_current: int,
        ip: float,
        q95: float,
        aspect: float,
        eps: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
        perim: float,
    ) -> float:
        """Function to calculate poloidal field from the plasma current

        This function calculates the poloidal field from the plasma current in Tesla,
        using a simple calculation using Ampere's law for conventional
        tokamaks, or for TARTs, a scaling from Peng, Galambos and
        Shipe (1992).

        Parameters
        ----------
        i_plasma_current :
            current scaling model to use
        ip :
            plasma current (A)
        q95 :
            95% flux surface safety factor
        aspect :
            plasma aspect ratio
        eps :
            inverse aspect ratio
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        kappa :
            plasma elongation
        delta :
            plasma triangularity
        perim :
            plasma perimeter (m)

        Returns
        -------
        :
            poloidal field in Tesla


        References
        ----------
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971

        """
        # Use Ampere's law using the plasma poloidal cross-section
        if i_plasma_current != 2:
            return constants.RMU0 * ip / perim
        # Use the relation from Peng, Galambos and Shipe (1992) [STAR code] otherwise
        ff1, ff2, _, _ = self.current._plascar_bpol(aspect, eps, kappa, delta)

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        return b_plasma_toroidal_on_axis * (ff1 + ff2) / (2.0 * np.pi * qbar)

    