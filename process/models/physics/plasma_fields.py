import logging

import numba as nb
import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure import (
    physics_variables,
)
from process.models.physics.plasma_current import PlasmaCurrent

logger = logging.getLogger(__name__)


class PlasmaFields(Model):
    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.current = PlasmaCurrent()

    def run(self):
        """Run the model. This model cannot yet be 'run'."""

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
        """Function to calculate surface-averaged poloidal field (⟨Bₚ(a)⟩) from the plasma current

        This function calculates the surface-averaged poloidal field from the plasma current in Tesla,
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
            surface-averaged poloidal field in Tesla ⟨Bₚ(a)⟩


        References
        ----------
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971

        """
        # Use Ampere's law using the plasma poloidal cross-section this simply returns
        # ⟨Bₚ(a)⟩
        if i_plasma_current != 2:
            return constants.RMU0 * ip / perim
        # Use the relation from Peng, Galambos and Shipe (1992) [STAR code] otherwise
        ff1, ff2, _, _ = self.current.plascar_bpol(aspect, eps, kappa, delta)

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        return b_plasma_toroidal_on_axis * (ff1 + ff2) / (2.0 * np.pi * qbar)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_inboard_toroidal_field(
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
    ) -> float:
        """Calculate the toroidal field at the plasma inboard midplane (Bᴛ(R₀-a))

        Parameters
        ----------
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        rmajor :
            plasma major radius (m)
        rminor :
            plasma minor radius (m)

        Returns
        -------
        :
            toroidal field at the plasma inboard midplane (T)
        """

        return rmajor * b_plasma_toroidal_on_axis / (rmajor - rminor)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_outboard_toroidal_field(
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
    ) -> float:
        """Calculate the toroidal field at the plasma outboard midplane (Bᴛ(R₀+a))

        Parameters
        ----------
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        rmajor :
            plasma major radius (m)
        rminor :
            plasma minor radius (m)

        Returns
        -------
        :
            toroidal field at the plasma outboard midplane (T)
        """

        return rmajor * b_plasma_toroidal_on_axis / (rmajor + rminor)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_toroidal_field_profile(
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
        n_plasma_profile_elements: int,
    ) -> np.ndarray:
        """Calculate the toroidal field profile across the plasma midplane

        Parameters
        ----------
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        rmajor :
            plasma major radius (m)
        rminor :
            plasma minor radius (m)
        n_plasma_profile_elements :
            Number of elements to use in the plasma profile calculation
        """

        # Calculate the toroidal field across the plasma
        # Calculate the toroidal field profile across the plasma (1/R dependence)
        # Double element size to include both sides of the plasma
        rho = np.linspace(
            rmajor - rminor,
            rmajor + rminor,
            2 * n_plasma_profile_elements,
        )

        # Avoid division by zero at the magnetic axis
        rho = np.where(rho == 0, 1e-10, rho)
        return rmajor * b_plasma_toroidal_on_axis / rho

    @staticmethod
    @nb.njit(cache=True)
    def calculate_total_magnetic_field(
        b_plasma_toroidal: float, b_plasma_poloidal: float
    ) -> float:
        """Calculate the total magnetic field at the plasma edge

        Parameters
        ----------
        b_plasma_toroidal :
            toroidal field at point of interest (T)
        b_plasma_poloidal :
            poloidal field at point of interest (T)

        Returns
        -------
        :
            total magnetic field at the plasma edge (T)
        """
        return np.sqrt(b_plasma_toroidal**2 + b_plasma_poloidal**2)

    def output(self):
        po.oheadr(self.outfile, "Plasma magnetic fields")

        po.ovarrf(
            self.outfile,
            "Vertical field at plasma (T)",
            "(b_plasma_vertical_required)",
            physics_variables.b_plasma_vertical_required,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Vacuum toroidal field at R (T)",
            "(b_plasma_toroidal_on_axis)",
            physics_variables.b_plasma_toroidal_on_axis,
        )
        po.ovarrf(
            self.outfile,
            "Toroidal field at plasma inboard (T)",
            "(b_plasma_inboard_toroidal)",
            physics_variables.b_plasma_inboard_toroidal,
        )
        po.ovarrf(
            self.outfile,
            "Toroidal field at plasma outboard (T)",
            "(b_plasma_outboard_toroidal)",
            physics_variables.b_plasma_outboard_toroidal,
        )

        for i in range(len(physics_variables.b_plasma_toroidal_profile)):
            po.ovarre(
                self.mfile,
                f"Toroidal field in plasma at point {i}",
                f"b_plasma_toroidal_profile{i}",
                physics_variables.b_plasma_toroidal_profile[i],
            )
        po.ovarrf(
            self.outfile,
            "Plasma surface averaged poloidal field (T)",
            "(b_plasma_surface_poloidal_average)",
            physics_variables.b_plasma_surface_poloidal_average,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Total field (sqrt(b_plasma_surface_poloidal_average^2 + b_plasma_toroidal_on_axis^2)) (T)",
            "(b_plasma_total)",
            physics_variables.b_plasma_total,
            "OP ",
        )
