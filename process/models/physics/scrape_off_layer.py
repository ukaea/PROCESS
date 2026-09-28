"""Module for calculating plasma scrape off layer physics"""

import logging

import numpy as np

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


class Zhang0DBoxModel(Model):
    """Zhang 0D box model for scrape-off-layer (SOL) point physics.

    This is an extension of the conduction-limited two-point model.

    References
    ----------
    [1] X. Zhang, F. M. Poli, E. D. Emdee, and M. Podesta, "Reduced physics model
    of the tokamak Scrape-Off-Layer for pulse design," Nuclear Materials and Energy,
    vol. 34, p. 101354, Mar. 2023, doi: 10.1016/j.nme.2022.101354.
    """

    @staticmethod
    def calculate_parallel_electron_conductivity(
        ion_charge_number: float, coulomb_logarithm: float
    ) -> float:
        """Calculate the electron conductivity coefficient [W m^-1 eV^-7/2]."""
        return 30692.0 / (ion_charge_number * coulomb_logarithm)

    @staticmethod
    def calculate_parallel_ion_conductivity(
        ion_charge_number: float,
        ion_mass_amu: float,
        coulomb_logarithm: float,
    ) -> float:
        """Calculate the ion conductivity coefficient [W m^-1 eV^-7/2]."""
        return 1249.0 / (ion_charge_number**4 * ion_mass_amu**0.5 * coulomb_logarithm)

    @staticmethod
    def calculate_electron_ion_coulomb_logarithm(
        electron_density: float,
        electron_temperature: float,
        ion_density: float,
        ion_temperature: float,
        ion_mass_amu: float,
        ion_charge_number: float,
    ) -> float:
        """Calculate the electron-ion Coulomb logarithm for density in m^-3 and eV."""
        ion_mass = ion_mass_amu * constants.ATOMIC_MASS_UNIT
        low_temperature_limit = ion_temperature * constants.ELECTRON_MASS / ion_mass
        high_temperature_limit = 10.0 * ion_charge_number**2

        if electron_temperature < low_temperature_limit:
            return 16.0 - np.log(
                (ion_density * 1.0e-6) ** 0.5
                / ion_temperature**1.5
                * ion_charge_number**2
                / (ion_mass / constants.PROTON_MASS)
            )
        if electron_temperature > high_temperature_limit:
            return 24.0 - np.log(
                (electron_density * 1.0e-6) ** 0.5 / electron_temperature
            )
        return 23.0 - np.log(
            (electron_density * 1.0e-6) ** 0.5
            * ion_charge_number
            / electron_temperature**1.5
        )

    @staticmethod
    def calculate_upstream_temperature(
        temp_target_ev: float,
        f_conduction: float,
        pflux_upstream_parallel: float,
        cond_parallel: float,
        len_sol_connection: float,
        flux_expansion: float,
    ) -> float:
        """Calculate the upstream temperature (Tᵤ) [eV] based on the target
        temperature and parallel heat flux. This works for any electron and ion
        conduction scenario.

        Parameters
        ----------
        temp_target_ev :
            Target temperature (Tₜ) [eV].
        f_conduction :
            Ratio of the parallel heat flux carried by thermal conduction to
            the total parallel heat flux
        pflux_upstream_parallel :
            Parallel particle heat flux upstream (qₗₗ,ᵤ) [W/m²].
        cond_parallel :
            Parallel particle thermal conductivity along the magnetic field lines
            (κₗₗ) [W/m eV⁻⁷⸍²].
        len_sol_connection :
            Connection length of the SOL along the magnetic field lines (Lₗₗ) [m].
        flux_expansion :
            Flux expansion factor of the magnetic field lines in the SOL (fₓ).

        Returns
        -------
        :
            Upstream particle species temperature [eV].

        Notes
        -----
        - This equation is given by Equation 10 in Reference [1] for the class
        """
        return (
            (temp_target_ev) ** (7 / 2)
            + (
                f_conduction
                * (7 / 4)
                * pflux_upstream_parallel
                / cond_parallel
                * len_sol_connection
                * np.log(flux_expansion)
                / (flux_expansion - 1)
            )
        ) ** (2 / 7)

    @staticmethod
    def calculate_target_density(
        flux_expansion: float,
        particle_flux_density: float,
        vel_target_sound: float,
        f_plasma_particles_lcfs_recycled: float,
    ) -> float:
        """Calculate the target density (nₜ) [m⁻³] for a given particle species

        Parameters
        ----------
        flux_expansion :
            Flux expansion factor of the magnetic field lines in the SOL (fₓ).
        particle_flux_density :
            Particle flux density of the given species (Γ) [particles/(s m²)].
        vel_target_sound :
            Sound speed at the target (cₛ,ₜ) [m/s].
        f_plasma_particles_lcfs_recycled :
            Fraction of plasma particles at the LCFS/Divertor that are recycled back
            into the SOL (R).

        Returns
        -------
        :
            Target density (nₜ) [m⁻³].

        Notes
        -----
        - This equation is given by Equation 16 in Reference [1] for the class
        """
        return (
            ((flux_expansion + 1) / (2 * flux_expansion))
            * particle_flux_density
            / (vel_target_sound * (1.0 - f_plasma_particles_lcfs_recycled))
        )

    @staticmethod
    def calculate_target_sound_speed(
        temp_target_ion_ev: float, temp_target_electron_ev: float, m_ion: float
    ) -> float:
        """Calculate the sound speed at the target (cₛ,ₜ) [m/s] for a given ion species.

        Parameters
        ----------
        temp_target_ion_ev :
            Ion temperature at the target (Tᵢ,ₜ) [eV].
        temp_target_electron_ev :
            Electron temperature at the target (Tₑ,ₜ) [eV].
        m_ion :
            Mass of the ion species (mᵢ) [kg].

        Returns
        -------
        :
            Sound speed at the target (cₛ,ₜ) [m/s].

        Notes
        -----
        - This equation is given by Equation 5 in Reference [1] for the class

        """
        return np.sqrt(
            (temp_target_ion_ev + temp_target_electron_ev)
            * constants.ELEMENTARY_CHARGE
            / m_ion
        )

    @classmethod
    def solve_at_point(
        cls,
        pflux_electron_parallel_w: float,
        particle_flux_electron_parallel: float,
        pflux_ion_parallel_w: float,
        particle_flux_ion_parallel: float,
        ion_mass_amu: float,
        ion_charge_number: float,
        len_parallel: float,
        flux_expansion: float,
        recycling_fraction: float,
        electron_power_loss_fraction: float,
        ion_power_loss_fraction: float,
        momentum_loss_factor: float,
        conduction_loss_fraction: float,
        tolerance: float = 1.0e-4,
        max_iterations: int = 1000,
    ) -> dict[str, float]:
        """Solve target and upstream SOL conditions at a single outboard point.

        Parameters use PROCESS conventions: powers are fluxes in W/m^2, particle
        fluxes are in particles/(s m^2), lengths are in m, temperatures returned
        in eV, and densities returned in m^-3.
        """
        if (
            pflux_electron_parallel_w <= 0.0
            or particle_flux_electron_parallel <= 0.0
            or pflux_ion_parallel_w <= 0.0
            or particle_flux_ion_parallel <= 0.0
            or ion_mass_amu <= 0.0
            or ion_charge_number <= 0.0
            or len_parallel <= 0.0
            or flux_expansion < 1.0
            or not 0.0 <= recycling_fraction < 1.0
            or not 0.0 <= electron_power_loss_fraction <= 1.0
            or not 0.0 <= ion_power_loss_fraction <= 1.0
            or momentum_loss_factor <= 0.0
            or not 0.0 <= conduction_loss_fraction <= 1.0
            or tolerance <= 0.0
            or max_iterations < 1
        ):
            raise ValueError("Invalid Zhang 0D box model input.")

        ion_mass = ion_mass_amu * constants.ATOMIC_MASS_UNIT
        particle_loss_factor = 1.0 - recycling_fraction
        heat_transmission_factor = 2.5 + 3.0 * particle_loss_factor
        target_electron_temperature = 100.0
        target_ion_temperature = 100.0
        target_electron_density = 1.0e19
        target_ion_density = 1.0e19

        for _ in range(max_iterations):
            coulomb_logarithm = cls.calculate_electron_ion_coulomb_logarithm(
                electron_density=target_electron_density,
                electron_temperature=target_electron_temperature,
                ion_density=target_ion_density,
                ion_temperature=target_ion_temperature,
                ion_mass_amu=ion_mass_amu,
                ion_charge_number=ion_charge_number,
            )
            electron_conductivity = cls.calculate_parallel_electron_conductivity(
                ion_charge_number, coulomb_logarithm
            )
            ion_conductivity = cls.calculate_parallel_ion_conductivity(
                ion_charge_number, ion_mass_amu, coulomb_logarithm
            )

            target_electron_energy = (
                2.0
                * particle_loss_factor
                / heat_transmission_factor
                * (1.0 - electron_power_loss_fraction)
                / (flux_expansion + 1.0)
                * pflux_electron_parallel_w
                / particle_flux_electron_parallel
            )
            target_ion_energy = (
                2.0
                * particle_loss_factor
                / heat_transmission_factor
                * (1.0 - ion_power_loss_fraction)
                / (flux_expansion + 1.0)
                * pflux_ion_parallel_w
                / particle_flux_ion_parallel
            )
            new_target_electron_temperature = (
                target_electron_energy / constants.ELECTRON_VOLT
            )
            new_target_ion_temperature = target_ion_energy / constants.ELECTRON_VOLT
            target_sound_speed = np.sqrt(
                (target_electron_energy + target_ion_energy) / ion_mass
            )
            new_target_electron_density = (
                (flux_expansion + 1.0)
                / (2.0 * flux_expansion)
                * particle_flux_electron_parallel
                / (particle_loss_factor * target_sound_speed)
            )
            new_target_ion_density = (
                (flux_expansion + 1.0)
                / (2.0 * flux_expansion)
                * particle_flux_ion_parallel
                / (particle_loss_factor * target_sound_speed)
            )
            flux_expansion_logarithm = (
                1.0
                if flux_expansion == 1.0
                else np.log(flux_expansion) / (flux_expansion - 1.0)
            )
            upstream_electron_temperature = (
                new_target_electron_temperature**3.5
                + conduction_loss_fraction
                * 7.0
                / 4.0
                * pflux_electron_parallel_w
                / electron_conductivity
                * len_parallel
                * flux_expansion_logarithm
            ) ** (2.0 / 7.0)
            upstream_ion_temperature = (
                new_target_ion_temperature**3.5
                + conduction_loss_fraction
                * 7.0
                / 4.0
                * pflux_ion_parallel_w
                / ion_conductivity
                * len_parallel
                * flux_expansion_logarithm
            ) ** (2.0 / 7.0)
            upstream_density_factor = (
                (flux_expansion + 1.0)
                / (flux_expansion * momentum_loss_factor)
                * (new_target_ion_temperature + new_target_electron_temperature)
                / (upstream_ion_temperature + upstream_electron_temperature)
                / (particle_loss_factor * target_sound_speed)
            )
            upstream_electron_density = (
                upstream_density_factor * particle_flux_electron_parallel
            )
            upstream_ion_density = upstream_density_factor * particle_flux_ion_parallel

            relative_error = 0.5 * (
                abs(new_target_electron_temperature - target_electron_temperature)
                / new_target_electron_temperature
                + abs(new_target_ion_temperature - target_ion_temperature)
                / new_target_ion_temperature
            )
            target_electron_temperature = new_target_electron_temperature
            target_ion_temperature = new_target_ion_temperature
            target_electron_density = new_target_electron_density
            target_ion_density = new_target_ion_density
            if relative_error <= tolerance:
                return {
                    "temp_plasma_divertor_electron_ev": target_electron_temperature,
                    "temp_plasma_divertor_ion_ev": target_ion_temperature,
                    "nd_plasma_divertor_electron": target_electron_density,
                    "nd_plasma_divertor_ion": target_ion_density,
                    "temp_plasma_outboard_midplane_electron_ev": upstream_electron_temperature,
                    "temp_plasma_outboard_midplane_ion_ev": upstream_ion_temperature,
                    "nd_plasma_outboard_midplane_electron": upstream_electron_density,
                    "nd_plasma_outboard_midplane_ion": upstream_ion_density,
                }

        raise RuntimeError("Zhang 0D box model did not converge.")
