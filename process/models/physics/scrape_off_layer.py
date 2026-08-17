"""Module for calculating plasma scrape off layer physics"""

import logging

import numpy as np
import scipy

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure.physics_variables import OutbordSOLPowerDecayLengthModel

logger = logging.getLogger(__name__)


class ScrapeOffLayer(Model):
    """Model for calculating plasma scrape off layer physics."""

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

        # NOTE: converting plasma_current from A to MA
        self.data.physics.len_plasma_sol_mast14_power_decay_2 = (
            self.calculate_mast2014_sol_power_decay_length_2(
                p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
                cur_plasma_ma=self.data.physics.plasma_current / 1e6,
            )
        )

        # Set to user input if OutbordSOLPowerDecayLengthModel = 1/USER_INUT

        if (
            OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            )
            == OutbordSOLPowerDecayLengthModel.EICH_2013
        ):
            self.data.physics.len_sol_outboard_power_decay = (
                self.data.physics.len_plasma_sol_eich13_power_decay
            )
        elif (
            OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            )
            == OutbordSOLPowerDecayLengthModel.MAST_2014_1
        ):
            self.data.physics.len_sol_outboard_power_decay = (
                self.data.physics.len_plasma_sol_mast14_power_decay_1
            )
        elif (
            OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            )
            == OutbordSOLPowerDecayLengthModel.MAST_2014_2
        ):
            self.data.physics.len_sol_outboard_power_decay = (
                self.data.physics.len_plasma_sol_mast14_power_decay_2
            )

        self.data.physics.a_plasma_outboard_sol_parallel = self.calculate_upstream_sol_outboard_parallel_area(  # noqa: E501
            rmajor=self.data.physics.rmajor,
            rminor=self.data.physics.rminor,
            len_plasma_sol_power_decay=self.data.physics.len_sol_outboard_power_decay,
            b_plasma_outboard_total=self.data.physics.b_plasma_outboard_total,
            b_plasma_surface_poloidal_average=self.data.physics.b_plasma_surface_poloidal_average,
        )

        self.data.physics.a_plasma_outboard_sol_eich13_parallel = self.calculate_upstream_sol_outboard_parallel_area(  # noqa: E501
            rmajor=self.data.physics.rmajor,
            rminor=self.data.physics.rminor,
            len_plasma_sol_power_decay=self.data.physics.len_plasma_sol_eich13_power_decay,
            b_plasma_outboard_total=self.data.physics.b_plasma_outboard_total,
            b_plasma_surface_poloidal_average=self.data.physics.b_plasma_surface_poloidal_average,
        )

        self.data.physics.pflux_plasma_outboard_sol_parallel_mw = (
            self.data.physics.p_plasma_separatrix_mw
            / self.data.physics.a_plasma_outboard_sol_parallel
        )

        self.data.physics.pflux_plasma_outboard_sol_eich13_parallel_mw = (
            self.data.physics.p_plasma_separatrix_mw
            / self.data.physics.a_plasma_outboard_sol_eich13_parallel
        )

    def output(self) -> None:
        """Output plasma scrape off layer physics information."""
        po.oheadr(self.outfile, "Plasma Scrape Off Layer")

        po.osubhd(self.outfile, "Power Decay Lengths (λ_q):")

        po.ovarre(
            self.outfile,
            "Outboard SOL power decay length (λ_q) [m]",
            "(len_sol_outboard_power_decay)",
            self.data.physics.len_sol_outboard_power_decay,
        )
        po.ocmmnt(
            self.outfile,
            "-> "
            + OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            ).description
            + " ",
        )
        po.oblnkl(self.outfile)
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
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")

        po.osubhd(self.outfile, "Upstream Outboard SOL Parallel Area and Power Flux:")

        po.ovarre(
            self.outfile,
            "Plasma outboard midplane SOL parallel area (Aₗₗ,ᵤ) [m²]",
            "(a_plasma_outboard_sol_parallel)",
            self.data.physics.a_plasma_outboard_sol_parallel,
        )
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Plasma outboard midplane Eich 2013 SOL parallel area (Aₗₗ,ᵤ) [m²]",
            "(a_plasma_outboard_sol_eich13_parallel)",
            self.data.physics.a_plasma_outboard_sol_eich13_parallel,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma outboard midplane SOL parallel power flux (qₗₗ,ᵤ) [MW/m²]",
            "(pflux_plasma_outboard_sol_parallel_mw)",
            self.data.physics.pflux_plasma_outboard_sol_parallel_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma outboard midplane Eich 2013 SOL parallel power flux (qₗₗ,ᵤ) [MW/m²]",
            "(pflux_plasma_outboard_sol_eich13_parallel_mw)",
            self.data.physics.pflux_plasma_outboard_sol_eich13_parallel_mw,
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
            Power crossing the separatrix (Pₛₑₚ) [MW]
        rmajor : float
            Major radius of the plasma (R₀) [m]
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (⟨Bₚₒₗ(a)⟩) [T]
        aspect : float
            Aspect ratio of the plasma (A)

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
            Power crossing the separatrix (Pₛₑₚ) [MW]
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (⟨Bₚₒₗ(a)⟩) [T]

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
            Power crossing the separatrix (Pₛₑₚ) [MW]
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

    @staticmethod
    def calculate_upstream_sol_outboard_parallel_area(
        rmajor: float,
        rminor: float,
        len_plasma_sol_power_decay: float,
        b_plasma_outboard_total: float,
        b_plasma_surface_poloidal_average: float,
    ) -> float:
        """Calculate the outboard SOL upstream parallel area (Aₗₗ,ᵤ) [m²].

        Parameters
        ----------
        rmajor : float
            Major radius of the plasma (R₀) [m]
        rminor : float
            Minor radius of the plasma (a) [m]
        len_plasma_sol_power_decay : float
            Power decay length (λ_q) [m]
        b_plasma_outboard_total : float
            Total magnetic field at the plasma outboard (Bₜₒₜ(R₀+a)) [T]
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (⟨Bₚₒₗ(a)⟩) [T]

        Returns
        -------
        float
            Upstream outboard SOL parallel area (Aₗₗ,ᵤ) [m²]

        References
        ----------
        [1] P. C. Stangeby, “The Plasma Boundary of Magnetic Fusion Devices,” Jan. 2000,
        doi: 10.1201/9780367801489.

        [2] S. S. Henderson et al., “An overview of the STEP divertor design and the
        simple models driving the plasma exhaust scenario,” Nuclear Fusion, vol. 65,
        no. 1, pp. 016033-016033, Nov. 2024, doi: 10.1088/1741-4326/ad93e7.

        """
        return (
            (2 * np.pi * (rmajor + rminor))
            * len_plasma_sol_power_decay
            * (b_plasma_surface_poloidal_average / b_plasma_outboard_total)
        )

    @staticmethod
    def calculate_outboard_midplane_near_sol_radial_profile(
        rmajor: float,
        rminor: float,
        len_plasma_sol_power_decay: float,
        pflux_plasma_outboard_sol_parallel_mw: float,
        r: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the outboard midplane near SOL radial profile (qₗₗ(r)) [MW/m²].

        Parameters
        ----------
        rmajor : float
            Major radius of the plasma (R₀) [m]
        rminor : float
            Minor radius of the plasma (a) [m]
        len_plasma_sol_power_decay : float
            Power decay length (λ_q) [m]
        pflux_plasma_outboard_sol_parallel_mw : float
            Parallel power flux at the outboard midplane (qₗₗ,ᵤ) [MW/m²]
        r : float|np.ndarray
            Radial position(s) at which to calculate the SOL profile [m]

        Returns
        -------
        float|np.ndarray
            Outboard midplane SOL radial profile (qₗₗ(r)) [MW/m²]

        Notes
        -----
        - The exponential model is highly valid in the "near-SOL" (typically the first
        few millimeters to a centimeter outside the separatrix). In this region, parallel
        heat transport is dominated by classical electron heat conduction
        (Spitzer-Härm conductivity), which is vastly faster than perpendicular diffusion.
        This competition between fast parallel conduction and slow perpendicular
        diffusion naturally produces an exponential radial profile.

        - The midplane exponential assumes steady-state H-mode conditions without the
        massive, transient convective bursts caused by ELMs, which momentarily
        flatten the entire midplane profile.

        References
        ----------
        [1] T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode
        power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9,
        p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

        """
        if np.any(r < (rmajor + rminor)):
            raise ValueError(
                f"Radial position r={r} must be greater than or equal to the plasma "
                f"edge (rmajor + rminor)={rmajor + rminor}."
            )

        return pflux_plasma_outboard_sol_parallel_mw * np.exp(
            -(r - (rmajor + rminor)) / len_plasma_sol_power_decay
        )

    @staticmethod
    def calculate_eich_target_heat_flux_profile(
        rmajor: float,
        rminor: float,
        pflux_plasma_sol_parallel_mw: float,
        len_plasma_sol_power_decay: float,
        f_b_div_flux_expansion: float,
        len_plasma_sol_power_spreading: float,
        plux_target_background_heat_flux_mw: float,
        r: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the Eich target heat flux profile (qₜ(r)) [MW/m²].

        Parameters
        ----------
        rmajor : float
            Major radius of the plasma (R₀) [m]
        rminor : float
            Minor radius of the plasma (a) [m]
        pflux_plasma_sol_parallel_mw : float
            Parallel power flux at the outboard midplane (qₗₗ,ᵤ) [MW/m²]
        len_plasma_sol_power_decay : float
            Power decay length (λ_q) [m]
        f_b_div_flux_expansion : float
            Divertor flux expansion factor (fₓ) [-]
        len_plasma_sol_power_spreading : float
            Power spreading length in the divertor (S) [m]
        plux_target_background_heat_flux_mw : float
            Background heat flux at the divertor target [MW/m²]
        r : float|np.ndarray
            Radial position(s) at which to calculate the target heat flux profile [m]

        Returns
        -------
        float|np.ndarray
            Eich target heat flux profile (qₜ(r)) [MW/m²]

        Notes
        -----
        - The Eich target heat flux profile is derived from the midplane exponential
        profile, taking into account the magnetic geometry and flux expansion between
        the midplane and the divertor target. The profile is typically characterized by
        a combination of an exponential decay and a Gaussian spreading due to cross-field
        transport in the divertor leg.

        References
        ----------
        [1] T. Eich, B. Sieglin, A. Scarabosio, W. Fundamenski, R. J. Goldston, and
        A. Herrmann, “Inter-ELM Power Decay Length for JET and ASDEX Upgrade: Measurement
        and Comparison with Heuristic Drift-Based Model,” Physical Review Letters,
        vol. 107, no. 21, Nov. 2011, doi: https://doi.org/10.1103/PhysRevLett.107.215001

        [2] T. Eich et al., “Scaling of the tokamak near the scrape-off layer H-mode
        power width and implications for ITER,” Nuclear Fusion, vol. 53, no. 9,
        p. 093031, Aug. 2013, doi: 10.1088/0029-5515/53/9/093031.

        """
        return (pflux_plasma_sol_parallel_mw / 2) * np.exp(
            (
                (len_plasma_sol_power_spreading)
                / (2 * len_plasma_sol_power_decay * f_b_div_flux_expansion)
            )
            ** 2
            - (
                (r - (rmajor + rminor))
                * f_b_div_flux_expansion
                / (len_plasma_sol_power_spreading * f_b_div_flux_expansion)
            )
        ) * scipy.special.erfc(
            (
                len_plasma_sol_power_spreading
                / (2 * len_plasma_sol_power_decay * f_b_div_flux_expansion)
            )
            - (
                (r - (rmajor + rminor))
                * f_b_div_flux_expansion
                / (len_plasma_sol_power_spreading)
            )
        ) + plux_target_background_heat_flux_mw
