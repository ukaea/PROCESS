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

        self.data.physics.len_plasma_sol_eich11_jet_power_decay = (
            self.calculate_eich2011_jet_sol_power_decay_length(
                b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
                q_cyl=self.data.physics.qstar,
                p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
            )
        )

        self.data.physics.len_plasma_sol_eich11_jet_asdex_power_decay = (
            self.calculate_eich2011_jet_asdex_sol_power_decay_length(
                b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
                q_cyl=self.data.physics.qstar,
                p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
                rmajor=self.data.physics.rmajor,
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
        elif (
            OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            )
            == OutbordSOLPowerDecayLengthModel.EICH_2011_JET
        ):
            self.data.physics.len_sol_outboard_power_decay = (
                self.data.physics.len_plasma_sol_eich11_jet_power_decay
            )
        elif (
            OutbordSOLPowerDecayLengthModel(
                self.data.physics.i_len_sol_outboard_power_decay
            )
            == OutbordSOLPowerDecayLengthModel.EICH_2011_JET_ASDEX
        ):
            self.data.physics.len_sol_outboard_power_decay = (
                self.data.physics.len_plasma_sol_eich11_jet_asdex_power_decay
            )

        self.data.physics.len_sol_inboard_power_decay = (
            self.data.physics.f_len_sol_power_decay_inboard_outboard
            * self.data.physics.len_sol_outboard_power_decay
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
            * self.data.physics.f_p_div_outboard_separatrix
            / self.data.physics.a_plasma_outboard_sol_parallel
        )

        self.data.physics.pflux_plasma_outboard_sol_eich13_parallel_mw = (
            self.data.physics.p_plasma_separatrix_mw
            / self.data.physics.a_plasma_outboard_sol_eich13_parallel
        )

        self.data.physics.len_div_outboard_lower_scarabosio15_power_spreading = self.calculate_scarabosio2015_power_spreading_factor(  # noqa: E501
            p_plasma_separatrix_mw=self.data.physics.p_plasma_separatrix_mw,
            b_plasma_surface_poloidal_average=self.data.physics.b_plasma_surface_poloidal_average,
            nd_plasma_separatrix_electron_19=self.data.physics.nd_plasma_separatrix_electron
            / 1e19,
            rmajor=self.data.physics.rmajor,
        )

        self.data.physics.len_div_outboard_lower_power_spreading = (
            self.data.physics.len_div_outboard_lower_scarabosio15_power_spreading
        )

    def output(self) -> None:
        """Output plasma scrape off layer physics information."""
        po.oheadr(self.outfile, "Plasma Scrape Off Layer")

        po.osubhd(self.outfile, "Power Decay Lengths (λ_q):")

        po.ovarre(
            self.outfile,
            "Outboard SOL power decay length (λₒᵤₜ_q) [m]",
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
            "Inboard to outboard SOL power decay length ratio (λᵢₙ_q/λₒᵤₜ_q)",
            "(f_len_sol_power_decay_inboard_outboard)",
            self.data.physics.f_len_sol_power_decay_inboard_outboard,
        )
        po.ovarre(
            self.outfile,
            "Inboard SOL power decay length (λᵢₙ_q) [m]",
            "(len_sol_inboard_power_decay)",
            self.data.physics.len_sol_inboard_power_decay,
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
        po.ovarre(
            self.outfile,
            "Eich 2011 JET SOL power decay length (λ_q) [m]",
            "(len_plasma_sol_eich11_jet_power_decay)",
            self.data.physics.len_plasma_sol_eich11_jet_power_decay,
        )
        po.ovarre(
            self.outfile,
            "Eich 2011 JET + ASDEX Upgrade SOL power decay length (λ_q) [m]",
            "(len_plasma_sol_eich11_jet_asdex_power_decay)",
            self.data.physics.len_plasma_sol_eich11_jet_asdex_power_decay,
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
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Power Spreading Factors (S):")

        po.ovarre(
            self.outfile,
            "Outboard lower divertor power spreading factor (S) [m]",
            "(len_div_outboard_lower_power_spreading)",
            self.data.physics.len_div_outboard_lower_power_spreading,
        )
        po.ovarre(
            self.outfile,
            "Scarabosio 2015 H-mode power spreading factor (S) [m]",
            "(len_div_outboard_lower_scarabosio15_power_spreading)",
            self.data.physics.len_div_outboard_lower_scarabosio15_power_spreading,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Outboard lower divertor flux expansion factor for the divertor targets "
            "(fₓ)",
            "(f_b_div_outboard_lower_flux_expansion)",
            self.data.physics.f_b_div_outboard_lower_flux_expansion,
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
    def calculate_eich2011_jet_sol_power_decay_length(
        b_plasma_toroidal_on_axis: float,
        q_cyl: float,
        p_plasma_separatrix_mw: float,
    ) -> float:
        """Calculate the Eich 2011 JET SOL power decay length (λ_q).

        Parameters
        ----------
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field at the plasma axis (Bᴛ(R₀)) [T]
        q_cyl : float
            Cylindrical safety factor (q_cyl) [-]
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) [MW]

        Returns
        -------
        float
            Eich 2011 JET SOL power decay length (λ_q) [m]

        Notes
        -----
        - The fit values can be found in Table 2 of [1].
        - The scaling is done for type-I ELMy H-mode plasmas

        References
        ----------
        [1] T. Eich, B. Sieglin, A. Scarabosio, W. Fundamenski, Robert James Goldston,
        and A. Herrmann, “Inter-ELM Power Decay Length for JET and ASDEX Upgrade:
        Measurement and Comparison with Heuristic Drift-Based Model,”
        Physical Review Letters, vol. 107, no. 21, Nov. 2011,
        doi: https://doi.org/10.1103/PhysRevLett.107.215001.

        """
        return (
            0.7e-3
            * b_plasma_toroidal_on_axis**-0.84
            * q_cyl**1.23
            * p_plasma_separatrix_mw**0.14
        )

    @staticmethod
    def calculate_eich2011_jet_asdex_sol_power_decay_length(
        b_plasma_toroidal_on_axis: float,
        q_cyl: float,
        p_plasma_separatrix_mw: float,
        rmajor: float,
    ) -> float:
        """Calculate the Eich 2011 JET + ASDEX Upgrade SOL power decay length (λ_q).

        Parameters
        ----------
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field at the plasma axis (Bᴛ(R₀)) [T]
        q_cyl : float
            Cylindrical safety factor (q_cyl) [-]
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) [MW]
        rmajor : float
            Major radius of the plasma (R₀) [m]

        Returns
        -------
        float
            Eich 2011 JET + ASDEX Upgrade SOL power decay length (λ_q) [m]

        Notes
        -----
        - The fit values can be found in Table 2 of [1].
        - The scaling is done for type-I ELMy H-mode plasmas

        References
        ----------
        [1] T. Eich, B. Sieglin, A. Scarabosio, W. Fundamenski, Robert James Goldston,
        and A. Herrmann, “Inter-ELM Power Decay Length for JET and ASDEX Upgrade:
        Measurement and Comparison with Heuristic Drift-Based Model,”
        Physical Review Letters, vol. 107, no. 21, Nov. 2011,
        doi: https://doi.org/10.1103/PhysRevLett.107.215001.

        """
        return (
            0.73e-3
            * b_plasma_toroidal_on_axis**-0.78
            * q_cyl**1.2
            * p_plasma_separatrix_mw**0.1
            * rmajor**0.02
        )

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

        Raises
        ------
        ValueError
            If any radial position r is inside the plasma edge (r < rmajor + rminor)

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
        pflux_target_background_heat_flux_mw: float,
        r: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the Eich parallel target heat flux profile (qₗₗ,ₜ(r)) [MW/m²].

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
        pflux_target_background_heat_flux_mw : float
            Background heat flux at the divertor target [MW/m²]
        r : float|np.ndarray
            Radial position(s) at which to calculate the target heat flux profile [m]

        Returns
        -------
        float|np.ndarray
            Eich parallel target heat flux profile (qₗₗ,ₜ(r)) [MW/m²]

        Notes
        -----
        - The Eich parallel target heat flux profile is derived from the midplane
        exponential profile, taking into account the magnetic geometry and flux expansion
        between the midplane and the divertor target. The profile is typically
        characterized by a combination of an exponential decay and a Gaussian spreading
        due to cross-field transport in the divertor leg.

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
                / (len_plasma_sol_power_decay * f_b_div_flux_expansion)
            )
        ) * scipy.special.erfc(
            (
                len_plasma_sol_power_spreading
                / (2 * len_plasma_sol_power_decay * f_b_div_flux_expansion)
            )
            - ((r - (rmajor + rminor)) / (len_plasma_sol_power_spreading))
        ) + pflux_target_background_heat_flux_mw

    @staticmethod
    def calculate_scarabosio2015_power_spreading_factor(
        p_plasma_separatrix_mw: float,
        b_plasma_surface_poloidal_average: float,
        nd_plasma_separatrix_electron_19: float,
        rmajor: float,
    ) -> float:
        """Calculate the Scarabosio 2015 H-mode power spreading factor (S).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (Pₛₑₚ) [MW]
        b_plasma_surface_poloidal_average : float
            Poloidal magnetic field at the plasma surface (Bₚₒₗ(a))  [T]
        nd_plasma_separatrix_electron_19 : float
            Electron density at the separatrix (nₑ,ₛₑₚ) [10¹⁹ m⁻³]
        rmajor : float
            Major radius of the plasma (R₀) [m]

        Returns
        -------
        float
            Scarabosio 2015 H-mode power spreading factor (S) [m]

        Notes
        -----
        - The R² for the fit is 0.65

        References
        ----------
        [1] A. Scarabosio et al., “Scaling of the divertor power spreading (S-factor) in
        open and closed divertor operation in JET and ASDEX Upgrade,”
        Journal of Nuclear Materials, vol. 463, pp. 49-54, Aug. 2015,
        doi: 10.1016/j.jnucmat.2014.11.076.

        """
        return (
            0.12e-3
            * p_plasma_separatrix_mw**0.21
            * b_plasma_surface_poloidal_average**-0.82
            * nd_plasma_separatrix_electron_19**-0.02
            * rmajor**0.71
        )


class BasicTwoPointModel(Model):
    r"""
    Basic two-point model for scrape-off layer physics.

    This model provides a simplified representation of the plasma parameters
    along the scrape-off layer, assuming a two-point connection between the
    upstream (midplane) and downstream (target) conditions.

    Notes
    -----
    - Electron and ion temperatures are assumed to be equal.
    - Constant electron pressure, (n_tT_t=n_uT_u);
    - Braginskii parallel heat conduction (q_|| = -κ_e ∇_|| T_e)
    - Bohm speed (c_s=\sqrt{2eT_t/m_i});
    - Sheath heat transmission coefficient (\gamma).

    References
    ----------
    [1] Philippe Ghendrih, “The Plasma Boundary of Magnetic Fusion Devices,”
    Plasma Physics and Controlled Fusion, vol. 43, no. 2, pp. 223-224, Jan. 2001,
    doi: 10.1088/0741-3335/43/2/702.

    [2] P.C. Stangeby, “Basic physical processes and reduced models for plasma
    detachment,” vol. 60, no. 4, pp. 044022-044022, Mar. 2018,
    doi: 10.1088/1361-6587/aaacf6.
    """

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """Run the basic two-point model calculations."""

    def output(self):
        """Retrieve the output of the basic two-point model calculations."""

    @staticmethod
    def calculate_upstream_temperature(
        pflux_plasma_outboard_sol_parallel: float,
        len_connection: float,
        temp_target_ev: float = 0.0,
        electron_thermal_conductivity: float = 2000.0,
    ) -> float:
        """
        Calculate the upstream electron temperature (Tₑ,ᵤ) in the scrape-off layer.

        Parameters
        ----------
        pflux_plasma_outboard_sol_parallel : float
            Parallel heat flux along the scrape-off layer [W/m2]
        len_connection : float
            Connection length from the midplane to the target along the scrape-off layer [m]
        electron_thermal_conductivity : float, optional
            Electron thermal conductivity [W/m·eV^(-7/2)] (default is 2000.0)

        Returns
        -------
        float
            Upstream electron temperature [eV]

        Notes
        -----
        - The calculation assumes a simplified two-point model for the scrape-off layer.

        """
        return (
            temp_target_ev ** (7.0 / 2.0)
            + (7.0 / 2.0)
            * pflux_plasma_outboard_sol_parallel
            * len_connection
            / electron_thermal_conductivity
        ) ** (2.0 / 7.0)

    @staticmethod
    def calculate_connection_length(
        temp_electron_upstream_ev: float,
        temp_target_ev: float,
        pflux_plasma_outboard_sol_parallel: float,
        electron_thermal_conductivity: float = 2000.0,
    ) -> float:
        """Calculate the connection length for given upstream and target temperatures."""
        return (
            (2.0 / 7.0)
            * electron_thermal_conductivity
            / pflux_plasma_outboard_sol_parallel
            * (temp_electron_upstream_ev ** (7.0 / 2.0) - temp_target_ev ** (7.0 / 2.0))
        )

    def calculate_target_electron_temperature(
        self,
        m_ion_average: float,
        pflux_plasma_outboard_sol_parallel: float,
        nd_electron_upstream: float,
        temp_electron_upstream_ev: float,
        f_temp_ion_electron: float = 1.0,
        f_nd_electron_ion: float = 1.0,
        f_vel_ion_mach: float = 0.0,
        sheath_transmission_coefficient: float = 7.0,
    ) -> float:
        """
        Calculate the target electron temperature (Tₑ,ₜ) in the scrape-off layer.

        Parameters
        ----------
        m_ion_average : float
            Average ion mass [kg]
        pflux_plasma_outboard_sol_parallel : float
            Parallel heat flux along the scrape-off layer [W/m2]
        nd_electron_upstream : float
            Electron density at the upstream (midplane) [m^-3]
        temp_electron_upstream_ev : float
            Electron temperature at the upstream (midplane) [eV]
        sheath_transmission_coefficient : float, optional
            Sheath transmission coefficient (default is 7.0)

        Returns
        -------
        float
            Target electron temperature [eV]

        Notes
        -----
        - Taken from Equation 24 of Reference 2
        - The ion mass term is described as the fuel ion mass. Though we have used
        the total ion mass in the calculation, which includes impurities.

        """
        return (
            8
            * m_ion_average
            * pflux_plasma_outboard_sol_parallel**2
            / (
                sheath_transmission_coefficient**2
                * constants.ELECTRON_CHARGE
                * (
                    self.calculate_total_pressure(
                        nd_electron=nd_electron_upstream,
                        temp_electron_ev=temp_electron_upstream_ev,
                        f_temp_ion_electron=f_temp_ion_electron,
                        f_nd_electron_ion=f_nd_electron_ion,
                        f_vel_ion_mach=f_vel_ion_mach,
                    )
                )
                ** 2
            )
        )

    @staticmethod
    def calculate_total_pressure(
        nd_electron: float,
        temp_electron_ev: float,
        f_temp_ion_electron: float,
        f_nd_electron_ion: float,
        f_vel_ion_mach: float,
    ) -> float:
        """
        Calculate the total pressure in the scrape-off layer.

        Parameters
        ----------
        nd_electron : float
            Electron density [m^-3]
        temp_electron_ev : float
            Electron temperature [eV]
        f_temp_ion_electron : float
            Ratio of ion temperature to electron temperature
        f_nd_electron_ion : float
            Ratio of electron density to ion density
        f_vel_ion_mach : float
            Ion velocity in terms of Mach number

        Returns
        -------
        float
            Total pressure [Pa]

        Notes
        -----
        - Taken from Equation 20 of Reference 2
        - The total pressure is calculated considering both electron and ion
          contributions.
        - The ion contribution is scaled by the Mach number squared and the
          temperature and density ratios.

        """
        return (
            (1 + f_vel_ion_mach**2)
            * nd_electron
            * temp_electron_ev
            * (1 + f_temp_ion_electron / f_nd_electron_ion)
            * constants.ELECTRON_CHARGE
        )

    def solve_basic_two_point_model(
        self,
        len_connection: float,
        nd_electron_upstream: float,
        q_parallel: float,
        m_ion_average: float = 1.6726219e-27,
        electron_thermal_conductivity: float = 2000.0,
        sheath_transmission_coefficient: float = 7.0,
    ) -> tuple[float, float]:
        """Solve the coupled, no-loss basic two-point model in eV."""
        if len_connection <= 0.0:
            raise ValueError("Connection length must be positive.")
        if nd_electron_upstream <= 0.0:
            raise ValueError("Upstream density must be positive.")
        if q_parallel <= 0.0:
            raise ValueError("Parallel heat flux must be positive.")
        if m_ion_average <= 0.0:
            raise ValueError("Ion mass must be positive.")
        if electron_thermal_conductivity <= 0.0:
            raise ValueError("Thermal conductivity must be positive.")
        if sheath_transmission_coefficient <= 0.0:
            raise ValueError("Sheath transmission coefficient must be positive.")

        def residual(temp_electron_upstream_ev: float) -> float:
            temp_target_ev = self.calculate_target_electron_temperature(
                m_ion_average=m_ion_average,
                pflux_plasma_outboard_sol_parallel=q_parallel,
                nd_electron_upstream=nd_electron_upstream,
                temp_electron_upstream_ev=temp_electron_upstream_ev,
                sheath_transmission_coefficient=sheath_transmission_coefficient,
            )
            return temp_electron_upstream_ev - self.calculate_upstream_temperature(
                pflux_plasma_outboard_sol_parallel=q_parallel,
                len_connection=len_connection,
                temp_target_ev=temp_target_ev,
                electron_thermal_conductivity=electron_thermal_conductivity,
            )

        lower_bound = self.calculate_upstream_temperature(
            pflux_plasma_outboard_sol_parallel=q_parallel,
            len_connection=len_connection,
            electron_thermal_conductivity=electron_thermal_conductivity,
        )
        upper_bound = 2.0 * lower_bound
        while residual(upper_bound) < 0.0:
            upper_bound *= 2.0

        temp_electron_upstream_ev = scipy.optimize.brentq(
            residual, lower_bound, upper_bound
        )
        temp_target_ev = self.calculate_target_electron_temperature(
            m_ion_average=m_ion_average,
            pflux_plasma_outboard_sol_parallel=q_parallel,
            nd_electron_upstream=nd_electron_upstream,
            temp_electron_upstream_ev=temp_electron_upstream_ev,
            sheath_transmission_coefficient=sheath_transmission_coefficient,
        )
        return temp_electron_upstream_ev, temp_target_ev


    def calculate_temperature_profile(
        self,
        len_connection: float,
        nd_electron_upstream: float,
        q_parallel: float,
        m_ion_average: float = 1.6726219e-27,
        electron_thermal_conductivity: float = 2000.0,
        sheath_transmission_coefficient: float = 7.0,
        number_of_points: int = 100,
    ) -> np.ndarray:
        """Calculate the temperature profile using the two-point model.

        Parameters
        ----------
        len_connection : float
            Connection length along the magnetic field [m]
        nd_electron_upstream : float
            Electron density at the upstream (midplane) [m^-3]
        q_parallel : float
            Parallel heat flux, taken as an input from the mfile [W/m^2]
        m_ion_average : float, optional
            Average ion mass (default is proton mass)
        electron_thermal_conductivity : float, optional
            Electron thermal conductivity (default is 2000.0)
        sheath_transmission_coefficient : float, optional
            Sheath transmission coefficient (default is 7.0)
        number_of_points : int, optional
            Number of points in the temperature profile (default is 100)



        Raises
        ------
        ValueError
            If any of the input parameters are non-positive.

        """
        temp_upstream_target_ev, temp_target_ev = self.solve_basic_two_point_model(
            len_connection=len_connection,
            nd_electron_upstream=nd_electron_upstream,
            q_parallel=q_parallel,
            m_ion_average=m_ion_average,
            electron_thermal_conductivity=electron_thermal_conductivity,
            sheath_transmission_coefficient=sheath_transmission_coefficient,
        )

        temp_upstream_ev_pow = temp_upstream_target_ev ** (7.0 / 2.0)
        temp_target_ev_pow = temp_target_ev ** (7.0 / 2.0)

        profile_points = np.linspace(0.0, len_connection, number_of_points)

        # temp_upstream_ev_pow is, by construction, always >= temp_target_ev_pow
        # so the term being raised to (2/7) below is always non-negative and the
        # result is always real. Clip to zero to guard against small negative
        # values arising from floating-point rounding.
        profile_pow = temp_upstream_ev_pow - (
            temp_upstream_ev_pow - temp_target_ev_pow
        ) * (profile_points / len_connection)

        return np.clip(profile_pow, 0.0, None) ** (2.0 / 7.0)
