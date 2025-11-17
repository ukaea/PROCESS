import logging
from dataclasses import dataclass

import numpy as np
from scipy import integrate

from process import constants
from process.data_structure import physics_variables
from process.plasma_profiles import PlasmaProfile

logger = logging.getLogger(__name__)


REACTION_CONSTANTS_DT = {
    "bg": 34.3827,
    "mrc2": 1.124656e6,
    "cc1": 1.17302e-9,
    "cc2": 1.51361e-2,
    "cc3": 7.51886e-2,
    "cc4": 4.60643e-3,
    "cc5": 1.35000e-2,
    "cc6": -1.06750e-4,
    "cc7": 1.36600e-5,
}

REACTION_CONSTANTS_DHE3 = {
    "bg": 68.7508,
    "mrc2": 1.124572e6,
    "cc1": 5.51036e-10,
    "cc2": 6.41918e-3,
    "cc3": -2.02896e-3,
    "cc4": -1.91080e-5,
    "cc5": 1.35776e-4,
    "cc6": 0.0,
    "cc7": 0.0,
}

REACTION_CONSTANTS_DD1 = {
    "bg": 31.3970,
    "mrc2": 0.937814e6,
    "cc1": 5.43360e-12,
    "cc2": 5.85778e-3,
    "cc3": 7.68222e-3,
    "cc4": 0.0,
    "cc5": -2.96400e-6,
    "cc6": 0.0,
    "cc7": 0.0,
}

REACTION_CONSTANTS_DD2 = {
    "bg": 31.3970,
    "mrc2": 0.937814e6,
    "cc1": 5.65718e-12,
    "cc2": 3.41267e-3,
    "cc3": 1.99167e-3,
    "cc4": 0.0,
    "cc5": 1.05060e-5,
    "cc6": 0.0,
    "cc7": 0.0,
}


class FusionReactionRate:
    """
    Calculate the fusion reaction rate for each reaction case (DT, DHE3, DD1, DD2).

    This class provides methods to numerically integrate over the plasma cross-section
    to find the core plasma fusion power for different fusion reactions. The reactions
    considered are:
        - Deuterium-Tritium (D-T)
        - Deuterium-Helium-3 (D-3He)
        - Deuterium-Deuterium (D-D) first branch
        - Deuterium-Deuterium (D-D) second branch

    The class uses the Bosch-Hale parametrization to compute the volumetric fusion reaction
    rate <sigma v> and integrates over the plasma cross-section to find the core plasma
    fusion power.

    Attributes:
        plasma_profile (PlasmaProfile): The parameterized temperature and density profiles of the plasma.
        sigmav_dt_average (float): Average fusion reaction rate <sigma v> for D-T.
        dhe3_power_density (float): Fusion power density produced by the D-3He reaction.
        dd_power_density (float): Fusion power density produced by the D-D reactions.
        dt_power_density (float): Fusion power density produced by the D-T reaction.
        alpha_power_density (float): Power density of alpha particles produced.
        pden_non_alpha_charged_mw (float): Power density of charged particles produced.
        neutron_power_density (float): Power density of neutrons produced.
        fusion_rate_density (float): Fusion reaction rate density.
        alpha_rate_density (float): Alpha particle production rate density.
        proton_rate_density (float): Proton production rate density.
        f_dd_branching_trit (float): The rate of tritium producing D-D reactions to 3He ones.

    Methods:
        deuterium_branching(ion_temperature: float) -> float:
            Calculate the relative rate of tritium producing D-D reactions to 3He ones based on the volume averaged ion temperature.

        dt_reaction() -> None:
            Calculate the fusion reaction rate and power density for the deuterium-tritium (D-T) fusion reaction.

        dhe3_reaction() -> None:
            Calculate the fusion reaction rate and power density for the deuterium-helium-3 (D-3He) fusion reaction.

        dd_helion_reaction() -> None:
            Calculate the fusion reaction rate and power density for the deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces helium-3 (3He) and a neutron (n).

        dd_triton_reaction() -> None:
            Calculate the fusion reaction rate and power density for the deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces tritium (T) and a proton (p).

        sum_fusion_rates(alpha_power_add: float, charged_power_add: float, neutron_power_add: float, fusion_rate_add: float, alpha_rate_add: float, proton_rate_add: float) -> None:
            Sum the fusion rate at the end of each reaction.

        calculate_fusion_rates() -> None:
            Initiate all the fusion rate calculations.

        set_physics_variables() -> None:
            Set the required physics variables in the physics_variables and physics_module modules.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.

    """

    def __init__(self, plasma_profile: PlasmaProfile) -> None:
        """
        Initialize the FusionReactionRate class with the given plasma profile.

        Parameters:
            plasma_profile (PlasmaProfile): The parameterized temperature and density profiles of the plasma.

        Attributes:
            plasma_profile (PlasmaProfile): The parameterized temperature and density profiles of the plasma.
            sigmav_dt_average (float): Average fusion reaction rate <sigma v> for D-T.
            dhe3_power_density (float): Fusion power density produced by the D-3He reaction.
            dd_power_density (float): Fusion power density produced by the D-D reactions.
            dt_power_density (float): Fusion power density produced by the D-T reaction.
            alpha_power_density (float): Power density of alpha particles produced.
            pden_non_alpha_charged_mw (float): Power density of charged particles produced.
            neutron_power_density (float): Power density of neutrons produced.
            fusion_rate_density (float): Fusion reaction rate density.
            alpha_rate_density (float): Alpha particle production rate density.
            proton_rate_density (float): Proton production rate density.
            f_dd_branching_trit (float): The rate of tritium producing D-D reactions to 3He ones.
        """
        self.plasma_profile = plasma_profile
        self.sigmav_dt_average = 0.0
        self.dhe3_power_density = 0.0
        self.dd_power_density = 0.0
        self.dt_power_density = 0.0
        self.alpha_power_density = 0.0
        self.pden_non_alpha_charged_mw = 0.0
        self.neutron_power_density = 0.0
        self.fusion_rate_density = 0.0
        self.alpha_rate_density = 0.0
        self.proton_rate_density = 0.0
        self.f_dd_branching_trit = 0.0

    def deuterium_branching(self, ion_temperature: float) -> float:
        """
        Calculate the relative rate of tritium producing D-D reactions to 3He ones based on the volume averaged ion temperature

        Parameters:
            ion_temperature (float): Volume averaged ion temperature in keV

        The method updates the following attributes:
            -f_dd_branching_trit: The rate of tritium producing D-D reactions to 3He ones

        Notes:
            - For ion temperatures between 0.5 keV and 200 keV.
            - The deviation of the fit from the R-matrix branching ratio is always smaller than 0.5%.

        References:
            - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
              Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
              doi: https://doi.org/10.1088/0029-5515/32/4/i07.
        ‌
        """
        # Divide by 2 to get the branching ratio for the D-D reaction that produces tritium as the output
        # is just the ratio of the two normalized cross sections
        self.f_dd_branching_trit = (
            1.02934
            - 8.3264e-3 * ion_temperature
            + 1.7631e-4 * ion_temperature**2
            - 1.8201e-6 * ion_temperature**3
            + 6.9855e-9 * ion_temperature**4
        ) / 2.0

    def dt_reaction(self) -> None:
        """D + T --> 4He + n reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-tritium (D-T) fusion reaction. It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate <sigma v> and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.sigmav_dt_average: Average fusion reaction rate <sigma v> for D-T.
            - self.dt_power_density: Fusion power density produced by the D-T reaction.
            - self.alpha_power_density: Power density of alpha particles produced.
            - self.pden_non_alpha_charged_mw: Power density of charged particles produced.
            - self.neutron_power_density: Power density of neutrons produced.
            - self.fusion_rate_density: Fusion reaction rate density.
            - self.alpha_rate_density: Alpha particle production rate density.
            - self.proton_rate_density: Proton production rate density.

        Returns:
            None
        """
        # Initialize Bosch-Hale constants for the D-T reaction
        dt = BoschHaleConstants(**REACTION_CONSTANTS_DT)

        physics_variables.fusrat_plasma_dt_profile = (
            bosch_hale_reactivity(
                (
                    physics_variables.temp_plasma_ion_vol_avg_kev
                    / physics_variables.temp_plasma_electron_vol_avg_kev
                )
                * self.plasma_profile.teprofile.profile_y,
                dt,
            )
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.f_plasma_fuel_tritium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Calculate the fusion reaction rate integral using Simpson's rule
        sigmav = integrate.simpson(
            fusion_rate_integral(self.plasma_profile, dt),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Store the average fusion reaction rate
        self.sigmav_dt_average = sigmav

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.D_T_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            * (
                physics_variables.f_plasma_fuel_tritium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
        )

        # Power densities for different particles [MW/m^3]
        # Alpha particle gets approximately 20% of the fusion power
        alpha_power_density = (
            1.0 - constants.DT_NEUTRON_ENERGY_FRACTION
        ) * fusion_power_density
        pden_non_alpha_charged_mw = 0.0
        neutron_power_density = (
            constants.DT_NEUTRON_ENERGY_FRACTION * fusion_power_density
        )

        # Calculate the fusion rate density [reactions/m^3/second]
        fusion_rate_density = fusion_power_density / reaction_energy
        alpha_rate_density = fusion_rate_density
        proton_rate_density = 0.0

        # Update the cumulative D-T power density
        self.dt_power_density = fusion_power_density

        # Sum the fusion rates for all particles
        self.sum_fusion_rates(
            alpha_power_density,
            pden_non_alpha_charged_mw,
            neutron_power_density,
            fusion_rate_density,
            alpha_rate_density,
            proton_rate_density,
        )

    def dhe3_reaction(self) -> None:
        """D + 3He --> 4He + p reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-helium-3 (D-3He) fusion reaction. It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate <sigma v> and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.dhe3_power_density: Fusion power density produced by the D-3He reaction.
            - self.alpha_power_density: Power density of alpha particles produced.
            - self.pden_non_alpha_charged_mw: Power density of charged particles produced.
            - self.neutron_power_density: Power density of neutrons produced.
            - self.fusion_rate_density: Fusion reaction rate density.
            - self.alpha_rate_density: Alpha particle production rate density.
            - self.proton_rate_density: Proton production rate density.

        Returns:
            None
        """
        # Initialize Bosch-Hale constants for the D-3He reaction
        dhe3 = BoschHaleConstants(**REACTION_CONSTANTS_DHE3)

        # Calculate the fusion reaction rate integral using Simpson's rule
        sigmav = integrate.simpson(
            fusion_rate_integral(self.plasma_profile, dhe3),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )

        physics_variables.fusrat_plasma_dhe3_profile = (
            bosch_hale_reactivity(
                (
                    physics_variables.temp_plasma_ion_vol_avg_kev
                    / physics_variables.temp_plasma_electron_vol_avg_kev
                )
                * self.plasma_profile.teprofile.profile_y,
                dhe3,
            )
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.f_plasma_fuel_helium3
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.D_HELIUM_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            * (
                physics_variables.f_plasma_fuel_helium3
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
        )

        # Power densities for different particles [MW/m^3]
        # Alpha particle gets approximately 20% of the fusion power
        alpha_power_density = (
            1.0 - constants.DHELIUM_PROTON_ENERGY_FRACTION
        ) * fusion_power_density
        pden_non_alpha_charged_mw = (
            constants.DHELIUM_PROTON_ENERGY_FRACTION * fusion_power_density
        )
        neutron_power_density = 0.0

        # Calculate the fusion rate density [reactions/m^3/second]
        fusion_rate_density = fusion_power_density / reaction_energy
        alpha_rate_density = fusion_rate_density
        proton_rate_density = fusion_rate_density  # Proton production rate [m^3/second]

        # Update the cumulative D-3He power density
        self.dhe3_power_density = fusion_power_density

        # Sum the fusion rates for all particles
        self.sum_fusion_rates(
            alpha_power_density,
            pden_non_alpha_charged_mw,
            neutron_power_density,
            fusion_rate_density,
            alpha_rate_density,
            proton_rate_density,
        )

    def dd_helion_reaction(self) -> None:
        """D + D --> 3He + n reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces
        helium-3 (3He) and a neutron (n). It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate <sigma v> and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.dd_power_density: Fusion power density produced by the D-D reaction.
            - self.alpha_power_density: Power density of alpha particles produced.
            - self.pden_non_alpha_charged_mw: Power density of charged particles produced.
            - self.neutron_power_density: Power density of neutrons produced.
            - self.fusion_rate_density: Fusion reaction rate density.
            - self.alpha_rate_density: Alpha particle production rate density.
            - self.proton_rate_density: Proton production rate density.

        Returns:
            None
        """
        # Initialize Bosch-Hale constants for the D-D reaction
        dd1 = BoschHaleConstants(**REACTION_CONSTANTS_DD1)

        # Calculate the fusion reaction rate integral using Simpson's rule
        sigmav = integrate.simpson(
            fusion_rate_integral(self.plasma_profile, dd1),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )

        physics_variables.fusrat_plasma_dd_helion_profile = (
            bosch_hale_reactivity(
                (
                    physics_variables.temp_plasma_ion_vol_avg_kev
                    / physics_variables.temp_plasma_electron_vol_avg_kev
                )
                * self.plasma_profile.teprofile.profile_y,
                dd1,
            )
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.f_plasma_fuel_deuterium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.DD_HELIUM_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        # The power density is scaled by the branching ratio to simulate the different
        # product pathways
        fusion_power_density = (
            sigmav
            * reaction_energy
            * (1.0 - self.f_dd_branching_trit)
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
        )

        # Power densities for different particles [MW/m^3]
        # Neutron particle gets approximately 75% of the fusion power
        alpha_power_density = 0.0
        pden_non_alpha_charged_mw = (
            1.0 - constants.DD_NEUTRON_ENERGY_FRACTION
        ) * fusion_power_density
        neutron_power_density = (
            constants.DD_NEUTRON_ENERGY_FRACTION * fusion_power_density
        )

        # Calculate the fusion rate density [reactions/m^3/second]
        fusion_rate_density = fusion_power_density / reaction_energy
        alpha_rate_density = 0.0
        proton_rate_density = 0.0

        # Update the cumulative D-D power density
        self.dd_power_density += fusion_power_density

        # Sum the fusion rates for all particles
        self.sum_fusion_rates(
            alpha_power_density,
            pden_non_alpha_charged_mw,
            neutron_power_density,
            fusion_rate_density,
            alpha_rate_density,
            proton_rate_density,
        )

    def dd_triton_reaction(self) -> None:
        """D + D --> T + p reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces
        tritium (T) and a proton (p). It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate <sigma v> and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.dd_power_density: Fusion power density produced by the D-D reaction.
            - self.alpha_power_density: Power density of alpha particles produced.
            - self.pden_non_alpha_charged_mw: Power density of charged particles produced.
            - self.neutron_power_density: Power density of neutrons produced.
            - self.fusion_rate_density: Fusion reaction rate density.
            - self.alpha_rate_density: Alpha particle production rate density.
            - self.proton_rate_density: Proton production rate density.

        Returns:
            None
        """
        # Initialize Bosch-Hale constants for the D-D reaction
        dd2 = BoschHaleConstants(**REACTION_CONSTANTS_DD2)

        # Calculate the fusion reaction rate integral using Simpson's rule
        sigmav = integrate.simpson(
            fusion_rate_integral(self.plasma_profile, dd2),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )

        physics_variables.fusrat_plasma_dd_triton_profile = (
            bosch_hale_reactivity(
                (
                    physics_variables.temp_plasma_ion_vol_avg_kev
                    / physics_variables.temp_plasma_electron_vol_avg_kev
                )
                * self.plasma_profile.teprofile.profile_y,
                dd2,
            )
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.f_plasma_fuel_deuterium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.DD_TRITON_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        # The power density is scaled by the branching ratio to simulate the different
        # product pathways
        fusion_power_density = (
            sigmav
            * reaction_energy
            * self.f_dd_branching_trit
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            * (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
        )

        # Power densities for different particles [MW/m^3]
        alpha_power_density = 0.0
        pden_non_alpha_charged_mw = fusion_power_density
        neutron_power_density = 0.0

        # Calculate the fusion rate density [reactions/m^3/second]
        fusion_rate_density = fusion_power_density / reaction_energy
        alpha_rate_density = 0.0
        proton_rate_density = fusion_rate_density  # Proton production rate [m^3/second]

        # Update the cumulative D-D power density
        self.dd_power_density += fusion_power_density

        # Sum the fusion rates for all particles
        self.sum_fusion_rates(
            alpha_power_density,
            pden_non_alpha_charged_mw,
            neutron_power_density,
            fusion_rate_density,
            alpha_rate_density,
            proton_rate_density,
        )

    def sum_fusion_rates(
        self,
        alpha_power_add: float,
        charged_power_add: float,
        neutron_power_add: float,
        fusion_rate_add: float,
        alpha_rate_add: float,
        proton_rate_add: float,
    ) -> None:
        """Sum the fusion rate at the end of each reaction.

        This method updates the cumulative fusion power densities and reaction rates
        for alpha particles, charged particles, neutrons, and protons.

        Parameters:
            alpha_power_add (float): Alpha particle fusion power per unit volume [MW/m3].
            charged_power_add (float): Other charged particle fusion power per unit volume [MW/m3].
            neutron_power_add (float): Neutron fusion power per unit volume [MW/m3].
            fusion_rate_add (float): Fusion reaction rate per unit volume [reactions/m3/s].
            alpha_rate_add (float): Alpha particle production rate per unit volume [/m3/s].
            proton_rate_add (float): Proton production rate per unit volume [/m3/s].

        Returns:
            None
        """
        self.alpha_power_density += alpha_power_add
        self.pden_non_alpha_charged_mw += charged_power_add
        self.neutron_power_density += neutron_power_add
        self.fusion_rate_density += fusion_rate_add
        self.alpha_rate_density += alpha_rate_add
        self.proton_rate_density += proton_rate_add

    def calculate_fusion_rates(self) -> None:
        """
        Initiate all the fusion rate calculations.

        This method sequentially calculates the fusion reaction rates and power densities
        for the following reactions:
            - Deuterium-Tritium (D-T)
            - Deuterium-Helium-3 (D-3He)
            - Deuterium-Deuterium (D-D) first branch
            - Deuterium-Deuterium (D-D) second branch

        It updates the instance attributes for the cumulative power densities and reaction rates
        for alpha particles, charged particles, neutrons, and protons.

        Returns:
            None
        """
        self.dt_reaction()
        self.dhe3_reaction()
        self.dd_helion_reaction()
        self.dd_triton_reaction()

    def set_physics_variables(self) -> None:
        """
        Set the required physics variables in the physics_variables and physics_module modules.

        This method updates the global physics variables and module variables with the
        current instance's fusion power densities and reaction rates.

        Returns:
            None
        """
        physics_variables.pden_plasma_alpha_mw = self.alpha_power_density
        physics_variables.pden_non_alpha_charged_mw = self.pden_non_alpha_charged_mw
        physics_variables.pden_plasma_neutron_mw = self.neutron_power_density
        physics_variables.fusden_plasma = self.fusion_rate_density
        physics_variables.fusden_plasma_alpha = self.alpha_rate_density
        physics_variables.proton_rate_density = self.proton_rate_density
        physics_variables.sigmav_dt_average = self.sigmav_dt_average
        physics_variables.dt_power_density_plasma = self.dt_power_density
        physics_variables.dhe3_power_density = self.dhe3_power_density
        physics_variables.dd_power_density = self.dd_power_density
        physics_variables.f_dd_branching_trit = self.f_dd_branching_trit


@dataclass
class BoschHaleConstants:
    """DataClass which holds the constants required for the Bosch Hale calculation
    for a given fusion reaction.
    """

    bg: float
    mrc2: float
    cc1: float
    cc2: float
    cc3: float
    cc4: float
    cc5: float
    cc6: float
    cc7: float


def fusion_rate_integral(
    plasma_profile: PlasmaProfile, reaction_constants: BoschHaleConstants
) -> np.ndarray:
    """
    Evaluate the integrand for the fusion power integration.

    Parameters:
        plasma_profile (PlasmaProfile): Parameterised temperature and density profiles.
        reactionconstants (BoschHaleConstants): Bosch-Hale reaction constants.

    Returns:
        np.ndarray: Integrand for the fusion power.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """

    # Since the electron temperature profile is only calculated directly, we scale the ion temperature
    # profile by the ratio of the volume averaged ion to electron temperature
    ion_temperature_profile = (
        physics_variables.temp_plasma_ion_vol_avg_kev
        / physics_variables.temp_plasma_electron_vol_avg_kev
    ) * plasma_profile.teprofile.profile_y

    # Number of fusion reactions per unit volume per particle volume density (m^3/s)
    sigv = bosch_hale_reactivity(ion_temperature_profile, reaction_constants)

    # Integrand for the volume averaged fusion reaction rate sigmav:
    # sigmav = integral(2 rho (sigv(rho) ni(rho)^2) drho),
    # divided by the square of the volume-averaged ion density
    # to retain the dimensions m^3/s (this is multiplied back in later)

    # Set each point in the desnity profile as a fraction of the volume averaged desnity
    density_profile_normalised = (
        1.0 / physics_variables.nd_plasma_electrons_vol_avg
    ) * plasma_profile.neprofile.profile_y

    # Calculate a volume averaged fusion reaction integral that allows for fusion power to be scaled with
    # just the volume averged ion density.
    return (
        2.0 * plasma_profile.teprofile.profile_x * sigv * density_profile_normalised**2
    )


def bosch_hale_reactivity(
    ion_temperature_profile: np.ndarray, reaction_constants: BoschHaleConstants
) -> np.ndarray:
    """
    Calculate the volumetric fusion reaction rate 〈sigmav〉 (m^3/s) for one of four nuclear reactions using
    the Bosch-Hale parametrization.

    The valid range of the fit is 0.2 keV < t < 100 keV except for D-3He where it is 0.5 keV < t < 190 keV.

    Reactions:
        1. D-T reaction
        2. D-3He reaction
        3. D-D 1st reaction
        4. D-D 2nd reaction

    Parameters:
        ion_temperature_profile (np.ndarray): Plasma ion temperature profile in keV.
        reaction_constants (BoschHaleConstants): Bosch-Hale reaction constants.

    Returns:
        np.ndarray: Volumetric fusion reaction rate 〈sigmav〉 in m^3/s for each point in the ion temperature profile.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """
    theta1 = (
        ion_temperature_profile
        * (
            reaction_constants.cc2
            + ion_temperature_profile
            * (
                reaction_constants.cc4
                + ion_temperature_profile * reaction_constants.cc6
            )
        )
        / (
            1.0
            + ion_temperature_profile
            * (
                reaction_constants.cc3
                + ion_temperature_profile
                * (
                    reaction_constants.cc5
                    + ion_temperature_profile * reaction_constants.cc7
                )
            )
        )
    )
    theta = ion_temperature_profile / (1.0 - theta1)

    xi = ((reaction_constants.bg**2) / (4.0 * theta)) ** (1 / 3)

    # Volumetric reaction rate / reactivity 〈sigmav〉 (m^3/s)
    # Original form is in [cm^3/s], so multiply by 1.0e-6 to convert to [m^3/s]
    sigmav = (
        1.0e-6
        * reaction_constants.cc1
        * theta
        * np.sqrt(xi / (reaction_constants.mrc2 * ion_temperature_profile**3))
        * np.exp(-3.0 * xi)
    )

    # if t = 0, sigmav = 0. Use this mask to set sigmav to zero.
    t_mask = ion_temperature_profile == 0.0
    sigmav[t_mask] = 0.0

    # Return np.ndarray of sigmav for each point in the ion temperature profile
    return sigmav


def set_fusion_powers(
    f_alpha_electron: float,
    f_alpha_ion: float,
    p_beam_alpha_mw: float,
    pden_non_alpha_charged_mw: float,
    pden_plasma_neutron_mw: float,
    vol_plasma: float,
    pden_plasma_alpha_mw: float,
) -> tuple:
    """

    This function computes various fusion power metrics based on the provided plasma parameters.

    Parameters:
        f_alpha_electron (float): Fraction of alpha energy to electrons.
        f_alpha_ion (float): Fraction of alpha energy to ions.
        p_beam_alpha_mw (float): Alpha power from hot neutral beam ions (MW).
        pden_non_alpha_charged_mw (float): Other charged particle fusion power per unit volume (MW/m^3).
        pden_plasma_neutron_mw (float): Neutron fusion power per unit volume just from plasma (MW/m^3).
        vol_plasma (float): Plasma volume (m^3).
        pden_plasma_alpha_mw (float): Alpha power per unit volume just from plasma (MW/m^3).

    Returns:
        tuple: A tuple containing the following elements:
            - pden_neutron_total_mw (float): Neutron fusion power per unit volume from plasma and beams [MW/m^3].
            - p_plasma_alpha_mw (float): Alpha fusion power from only the plasma [MW].
            - p_alpha_total_mw (float): Total alpha fusion power from plasma and beams [MW].
            - p_plasma_neutron_mw (float): Neutron fusion power from only the plasma [MW].
            - p_neutron_total_mw (float): Total neutron fusion power from plasma and beams [MW].
            - p_non_alpha_charged_mw (float): Other total charged particle fusion power [MW].
            - pden_alpha_total_mw (float): Alpha power per unit volume, from beams and plasma [MW/m^3].
            - f_pden_alpha_electron_mw (float): Alpha power per unit volume to electrons [MW/m^3].
            - f_pden_alpha_ions_mw (float): Alpha power per unit volume to ions [MW/m^3].
            - p_charged_particle_mw (float): Charged particle fusion power [MW].
            - p_fusion_total_mw (float): Total fusion power [MW].

    References:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990

    """
    # Alpha power

    # Calculate alpha power produced just by the plasma
    p_plasma_alpha_mw = pden_plasma_alpha_mw * vol_plasma

    # Add neutral beam alpha power / volume
    pden_alpha_total_mw = pden_plasma_alpha_mw + (p_beam_alpha_mw / vol_plasma)

    # Total alpha power
    p_alpha_total_mw = pden_alpha_total_mw * vol_plasma

    # Neutron Power

    # Calculate neutron power produced just by the plasma
    p_plasma_neutron_mw = pden_plasma_neutron_mw * vol_plasma

    # Add extra neutron power from beams
    pden_neutron_total_mw = pden_plasma_neutron_mw + (
        (
            (
                constants.DT_NEUTRON_ENERGY_FRACTION
                / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION)
            )
            * p_beam_alpha_mw
        )
        / vol_plasma
    )

    # Total neutron power
    p_neutron_total_mw = pden_neutron_total_mw * vol_plasma

    # Charged particle power

    # Total non-alpha charged particle power
    p_non_alpha_charged_mw = pden_non_alpha_charged_mw * vol_plasma

    # Charged particle fusion power
    p_charged_particle_mw = p_alpha_total_mw + p_non_alpha_charged_mw

    # Total fusion power
    p_fusion_total_mw = p_alpha_total_mw + p_neutron_total_mw + p_non_alpha_charged_mw

    # Alpha power to electrons and ions (used with electron
    # and ion power balance equations only)
    # No consideration of pden_non_alpha_charged_mw here.
    f_pden_alpha_ions_mw = (
        physics_variables.f_p_alpha_plasma_deposited * pden_alpha_total_mw * f_alpha_ion
    )
    f_pden_alpha_electron_mw = (
        physics_variables.f_p_alpha_plasma_deposited
        * pden_alpha_total_mw
        * f_alpha_electron
    )

    return (
        pden_neutron_total_mw,
        p_plasma_alpha_mw,
        p_alpha_total_mw,
        p_plasma_neutron_mw,
        p_neutron_total_mw,
        p_non_alpha_charged_mw,
        pden_alpha_total_mw,
        f_pden_alpha_electron_mw,
        f_pden_alpha_ions_mw,
        p_charged_particle_mw,
        p_fusion_total_mw,
    )


def beam_fusion(
    beamfus0: float,
    betbm0: float,
    b_plasma_poloidal_average: float,
    b_plasma_toroidal_on_axis: float,
    c_beam_total: float,
    nd_plasma_electrons_vol_avg: float,
    nd_plasma_fuel_ions_vol_avg: float,
    ion_electron_coulomb_log: float,
    e_beam_kev: float,
    f_deuterium_plasma: float,
    f_tritium_plasma: float,
    f_beam_tritium: float,
    sigmav_dt_average: float,
    temp_plasma_electron_density_weighted_kev: float,
    temp_plasma_ion_density_weighted_kev: float,
    vol_plasma: float,
    n_charge_plasma_effective_mass_weighted_vol_avg: float,
) -> tuple:
    """
            Routine to calculate beam slowing down properties.

            This function computes the neutral beam beta component, hot beam ion density,
            and alpha power from hot neutral beam ions based on the provided plasma parameters.

            Parameters:
                beamfus0 (float): Multiplier for beam-background fusion calculation.
                betbm0 (float): Leading coefficient for neutral beam beta fraction.
                b_plasma_poloidal_average (float): Poloidal field (T).
                b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
                c_beam_total (float): Neutral beam current (A).
                nd_plasma_electrons_vol_avg (float): Electron density (m^-3).
                nd_plasma_fuel_ions_vol_avg (float): Fuel ion density (m^-3).
                ion_electron_coulomb_log (float): Ion-electron coulomb logarithm.
                e_beam_kev (float): Neutral beam energy (keV).
                f_deuterium_plasma (float): Deuterium fraction of main plasma.
                f_tritium_plasma (float): Tritium fraction of main plasma.
                f_beam_tritium (float): Tritium fraction of neutral beam.
                sigmav_dt_average (float): Profile averaged <sigma v> for D-T (m^3/s).
                temp_plasma_electron_density_weighted_kev (float): Density-weighted electron temperature (keV).
                temp_plasma_ion_density_weighted_kev (float): Density-weighted ion temperature (keV).
                vol_plasma (float): Plasma volume (m^3).
                n_charge_plasma_effective_mass_weighted_vol_avg (float): Mass weighted plasma effective charge.

            Returns:
                tuple: A tuple containing the following elements:
                    - beta_beam (float): Neutral beam beta component.
                    - nd_beam_ions_out (float): Hot beam ion density (m^-3).
                    - p_beam_alpha_mw (float): Alpha power from hot neutral beam ions (MW).

            Notes:
                - The function uses the Bosch-Hale parametrization to compute the reactivity.
                - The critical energy for electron/ion slowing down of the beam ion is calculated
                  for both deuterium and tritium neutral beams.
                - The function integrates the hot beam fusion reaction rate integrand over the
                  range of beam velocities up to the critical velocity.

             References:
                - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
                  Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
                  doi: https://doi.org/10.1088/0029-5515/32/4/i07.

                - J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3, pp. 1015-1103,
                  Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.

                - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,
                  vol. 9, no. 3, pp. 136-141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf
    ‌
    """

    # Beam ion slowing down time given by Deng Baiquan and G. A. Emmert 1987
    beam_slow_time = (
        1.99e19
        * (
            constants.M_DEUTERON_AMU * (1.0 - f_beam_tritium)
            + (constants.M_TRITON_AMU * f_beam_tritium)
        )
        * (temp_plasma_electron_density_weighted_kev**1.5 / nd_plasma_electrons_vol_avg)
        / ion_electron_coulomb_log
    )

    # Critical energy for electron/ion slowing down of the beam ion
    # (deuterium and tritium neutral beams, respectively) (keV)
    # Taken from J.W Sheffield, “The physics of magnetic fusion reactors,”
    critical_energy_deuterium = (
        14.8
        * constants.M_DEUTERON_AMU
        * temp_plasma_electron_density_weighted_kev
        * n_charge_plasma_effective_mass_weighted_vol_avg ** (2 / 3)
        * (ion_electron_coulomb_log + 4.0)
        / ion_electron_coulomb_log
    )
    critical_energy_tritium = critical_energy_deuterium * (
        constants.M_TRITON_AMU / constants.M_DEUTERON_AMU
    )

    # Deuterium and tritium ion densities
    deuterium_density = nd_plasma_fuel_ions_vol_avg * f_deuterium_plasma
    tritium_density = nd_plasma_fuel_ions_vol_avg * f_tritium_plasma

    (
        deuterium_beam_alpha_power,
        tritium_beam_alpha_power,
        hot_beam_density,
        beam_deposited_energy,
    ) = beamcalc(
        deuterium_density,
        tritium_density,
        e_beam_kev,
        critical_energy_deuterium,
        critical_energy_tritium,
        beam_slow_time,
        f_beam_tritium,
        c_beam_total,
        temp_plasma_ion_density_weighted_kev,
        vol_plasma,
        sigmav_dt_average,
    )

    # Neutral beam alpha power
    p_beam_alpha_mw = beamfus0 * (deuterium_beam_alpha_power + tritium_beam_alpha_power)

    # Neutral beam beta
    beta_beam = (
        betbm0
        * 4.03e-22
        * (2 / 3)
        * hot_beam_density
        * beam_deposited_energy
        / (b_plasma_toroidal_on_axis**2 + b_plasma_poloidal_average**2)
    )

    return beta_beam, hot_beam_density, p_beam_alpha_mw


def beamcalc(
    nd: float,
    nt: float,
    e_beam_kev: float,
    critical_energy_deuterium: float,
    critical_energy_tritium: float,
    beam_slow_time: float,
    f_beam_tritium: float,
    c_beam_total: float,
    temp_plasma_ion_vol_avg_kev: float,
    vol_plasma: float,
    svdt: float,
) -> tuple[float, float, float, float]:
    """
    Calculate neutral beam alpha power and ion energy.

    This function computes the alpha power generated from the interaction between
    hot beam ions and thermal ions in the plasma, as well as the hot beam ion density
    and average hot beam ion energy.

    Parameters:
        nd (float): Thermal deuterium density (m^-3).
        nt (float): Thermal tritium density (m^-3).
        e_beam_kev (float): Beam energy (keV).
        critical_energy_deuterium (float): Critical energy for electron/ion slowing down of the beam ion (deuterium neutral beam) (keV).
        critical_energy_tritium (float): Critical energy for beam slowing down (tritium neutral beam) (keV).
        beam_slow_time (float): Beam ion slowing down time on electrons (s).
        f_beam_tritium (float): Beam tritium fraction (0.0 = deuterium beam).
        c_beam_total (float): Beam current (A).
        temp_plasma_ion_vol_avg_kev (float): Thermal ion temperature (keV).
        vol_plasma (float): Plasma volume (m^3).
        svdt (float): Profile averaged <sigma v> for D-T (m^3/s).

    Returns:
        tuple[float, float, float, float]: A tuple containing the following elements:
            - Alpha power from deuterium beam-background fusion (MW).
            - Alpha power from tritium beam-background fusion (MW).
            - Hot beam ion density (m^-3).
            - Average hot beam ion energy (keV).

    Notes:
        - The function uses the Bosch-Hale parametrization to compute the reactivity.
        - The critical energy for electron/ion slowing down of the beam ion is calculated
          for both deuterium and tritium neutral beams.
        - The function integrates the hot beam fusion reaction rate integrand over the
          range of beam velocities up to the critical velocity.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.

        - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,
          vol. 9, no. 3, pp. 136-141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

        - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
          International Series of Monographs on Physics, Volume 149.

        - J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3, pp. 1015-1103,
          Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.

    """

    # D and T beam current fractions
    beam_current_deuterium = c_beam_total * (1.0 - f_beam_tritium)
    beam_current_tritium = c_beam_total * f_beam_tritium

    # At a critical energy the rate of loss to the ions becomes equal to that to the electrons,
    # and at lower energies the loss to the ions predominates.

    # Ratio of beam energy to critical energy for deuterium
    beam_energy_ratio_deuterium = e_beam_kev / critical_energy_deuterium

    # Calculate the characterstic time for the deuterium ions to slow down to the thermal energy, eg E = 0.
    characteristic_deuterium_beam_slow_time = (
        beam_slow_time / 3.0 * np.log(1.0 + (beam_energy_ratio_deuterium) ** 1.5)
    )

    deuterium_beam_density = (
        beam_current_deuterium
        * characteristic_deuterium_beam_slow_time
        / (constants.ELECTRON_CHARGE * vol_plasma)
    )

    # Ratio of beam energy to critical energy for tritium
    beam_energy_ratio_tritium = e_beam_kev / critical_energy_tritium

    # Calculate the characterstic time for the tritium to slow down to the thermal energy, eg E = 0.
    # Wesson, J. (2011) Tokamaks.
    characteristic_tritium_beam_slow_time = (
        beam_slow_time / 3.0 * np.log(1.0 + (beam_energy_ratio_tritium) ** 1.5)
    )

    tritium_beam_density = (
        beam_current_tritium
        * characteristic_tritium_beam_slow_time
        / (constants.ELECTRON_CHARGE * vol_plasma)
    )

    hot_beam_density = deuterium_beam_density + tritium_beam_density

    # Find the speed of the deuterium particle when it has the critical energy.
    # Re-arrange kinetic energy equation to find speed. Non-relativistic.
    deuterium_critical_energy_speed = np.sqrt(
        2.0
        * constants.KILOELECTRON_VOLT
        * critical_energy_deuterium
        / (constants.ATOMIC_MASS_UNIT * constants.M_DEUTERON_AMU)
    )

    # Find the speed of the tritium particle when it has the critical energy.
    # Re-arrange kinetic energy equation to find speed. Non-relativistic.
    tritium_critical_energy_speed = np.sqrt(
        2.0
        * constants.KILOELECTRON_VOLT
        * critical_energy_tritium
        / (constants.ATOMIC_MASS_UNIT * constants.M_TRITON_AMU)
    )

    # Source term representing the number of ions born per unit time per unit volume.
    # D.Baiquan et.al.  “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,
    # vol. 9, no. 3, pp. 136-141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

    source_deuterium = beam_current_deuterium / (constants.ELECTRON_CHARGE * vol_plasma)

    source_tritium = beam_current_tritium / (constants.ELECTRON_CHARGE * vol_plasma)

    pressure_coeff_deuterium = (
        constants.M_DEUTERON_AMU
        * constants.ATOMIC_MASS_UNIT
        * beam_slow_time
        * deuterium_critical_energy_speed**2
        * source_deuterium
        / (constants.KILOELECTRON_VOLT * 3.0)
    )
    pressure_coeff_tritium = (
        constants.M_TRITON_AMU
        * constants.ATOMIC_MASS_UNIT
        * beam_slow_time
        * tritium_critical_energy_speed**2
        * source_tritium
        / (constants.KILOELECTRON_VOLT * 3.0)
    )

    # Fast Ion Pressure
    # This is the same form as the ideal gas law pressure, P=1/3 * nmv^2
    deuterium_pressure = pressure_coeff_deuterium * fast_ion_pressure_integral(
        e_beam_kev, critical_energy_deuterium
    )
    tritium_pressure = pressure_coeff_tritium * fast_ion_pressure_integral(
        e_beam_kev, critical_energy_tritium
    )

    # Beam deposited energy
    # Find the energy from the ideal gas pressure, P=1/3 * nmv^2 = 2/3 * n<E>
    deuterium_deposited_energy = 1.5 * deuterium_pressure / deuterium_beam_density
    tritium_deposited_energy = 1.5 * tritium_pressure / tritium_beam_density

    total_depsoited_energy = (
        (deuterium_beam_density * deuterium_deposited_energy)
        + (tritium_beam_density * tritium_deposited_energy)
    ) / hot_beam_density

    hot_deuterium_rate = 1e-4 * beam_reaction_rate(
        constants.M_DEUTERON_AMU, deuterium_critical_energy_speed, e_beam_kev
    )

    hot_tritium_rate = 1e-4 * beam_reaction_rate(
        constants.M_TRITON_AMU, tritium_critical_energy_speed, e_beam_kev
    )

    deuterium_beam_alpha_power = alpha_power_beam(
        deuterium_beam_density,
        nt,
        hot_deuterium_rate,
        vol_plasma,
        temp_plasma_ion_vol_avg_kev,
        svdt,
    )
    tritium_beam_alpha_power = alpha_power_beam(
        tritium_beam_density,
        nd,
        hot_tritium_rate,
        vol_plasma,
        temp_plasma_ion_vol_avg_kev,
        svdt,
    )

    return (
        deuterium_beam_alpha_power,
        tritium_beam_alpha_power,
        hot_beam_density,
        total_depsoited_energy,
    )


def fast_ion_pressure_integral(e_beam_kev: float, critical_energy: float) -> float:
    """
    Calculate the fraction of initial beam energy given to the ions.

    This function computes the fraction of initial beam energy given to the ions. based on the neutral beam energy
    and the critical energy for electron/ion slowing down of the beam ion.

    Parameters:
        e_beam_kev (float): Neutral beam energy (keV).
        critical_energy (float): Critical energy for electron/ion slowing down of the beam ion (keV).

    Returns:
        float: Fraction of initial beam energy given to the ions.

    Notes:
        - The function uses the ratio of the beam energy to the critical energy to compute
          the hot ion energy parameter.
        - The calculation involves logarithmic and arctangent functions to account for
          the energy distribution of the hot ions.

    References:
        - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,
          vol. 9, no. 3, pp. 136-141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

        - W.A Houlberg, “Thermalization of an Energetic Heavy Ion in a Multi-species Plasma,” University of Wisconsin Fusion Technology Institute,
          Report UWFDM-103 1974, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm103.pdf
    """

    xcs = e_beam_kev / critical_energy
    xc = np.sqrt(xcs)

    t1 = xcs / 2.0
    t2 = np.log((xcs + 2.0 * xc + 1.0) / (xcs - xc + 1.0)) / 6.0

    xarg = (2.0 * xc - 1.0) / np.sqrt(3.0)
    t3 = np.arctan(xarg) / np.sqrt(3.0)
    t4 = (1 / np.sqrt(3.0)) * np.arctan(1 / np.sqrt(3.0))

    return t1 + t2 - t3 - t4


def alpha_power_beam(
    beam_ion_desnity: float,
    plasma_ion_desnity: float,
    sigv: float,
    vol_plasma: float,
    temp_plasma_ion_vol_avg_kev: float,
    sigmav_dt: float,
) -> float:
    """
    Calculate alpha power from beam-background fusion.

    This function computes the alpha power generated from the interaction between
    hot beam ions and thermal ions in the plasma.

    Parameters:
        beam_ion_desnity (float): Hot beam ion density (m^-3).
        plasma_ion_desnity (float): Thermal ion density (m^-3).
        sigv (float): Hot beam fusion reaction rate (m^3/s).
        vol_plasma (float): Plasma volume (m^3).
        temp_plasma_ion_vol_avg_kev (float): Thermal ion temperature (keV).
        sigmav_dt (float): Profile averaged <sigma v> for D-T (m^3/s).

    Returns:
        float: Alpha power from beam-background fusion (MW).

    Notes:
        - The function uses the Bosch-Hale parametrization to compute the reactivity.
        - The ratio of the profile-averaged <sigma v> to the reactivity at the given
          thermal ion temperature is used to scale the alpha power.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """
    # Calculate the reactivity ratio
    ratio = (
        sigmav_dt
        / bosch_hale_reactivity(
            np.array([temp_plasma_ion_vol_avg_kev]),
            BoschHaleConstants(**REACTION_CONSTANTS_DT),
        ).item()
    )

    # Calculate and return the alpha power
    return (
        beam_ion_desnity
        * plasma_ion_desnity
        * sigv
        * (constants.DT_ALPHA_ENERGY / 1e6)
        * vol_plasma
        * ratio
    )


def beam_reaction_rate(
    relative_mass_ion: float, critical_velocity: float, beam_energy_keV: float
) -> float:
    """
    Calculate the hot beam fusion reaction rate.

    This function computes the fusion reaction rate for hot beam ions
    using the critical velocity for electron/ion slowing down and the
    neutral beam energy.

    Parameters:
        relative_mass_ion (float): Relative atomic mass of the ion (e.g., approx 2.0 for D, 3.0 for T).
        critical_velocity (float): Critical velocity for electron/ion slowing down of the beam ion [m/s].
        beam_energy_keV(float): Neutral beam energy [keV].

    Returns:
        float: Hot beam fusion reaction rate (m^3/s).

    Notes:
        - The function integrates the hot beam fusion reaction rate integrand
          over the range of beam velocities up to the critical velocity.
        - The integration is performed using the quad function from scipy.integrate.

    References:
        - P J Knight, CCFE, Culham Science Centre
    """

    # Find the speed of the beam particle when it has the critical energy.
    # Re-arrange kinetic energy equation to find speed. Non-relativistic.
    beam_velocity = np.sqrt(
        (beam_energy_keV * constants.KILOELECTRON_VOLT)
        * 2.0
        / (relative_mass_ion * constants.ATOMIC_MASS_UNIT)
    )

    relative_velocity = beam_velocity / critical_velocity
    integral_coefficient = (
        3.0 * critical_velocity / np.log(1.0 + (relative_velocity**3))
    )

    fusion_integral = integrate.quad(
        _hot_beam_fusion_reaction_rate_integrand,
        0.0,
        relative_velocity,
        args=(critical_velocity,),
    )[0]

    return integral_coefficient * fusion_integral


def _hot_beam_fusion_reaction_rate_integrand(
    velocity_ratio: float, critical_velocity: float
) -> float:
    """
    Integrand function for the hot beam fusion reaction rate.

    This function computes the integrand for the hot beam fusion reaction rate
    based on the ratio of beam velocity to the critical velocity and the critical
    velocity for electron/ion slowing down of the beam ion.

    Parameters:
        velocity_ratio (float): Ratio of beam velocity to the critical velocity.
        critical_velocity (float): Critical velocity for electron/ion slowing down of the beam ion (m/s).

    Returns:
        float: Value of the integrand for the hot beam fusion reaction rate.

    Notes:
        - The function uses the ratio of the beam velocity to the critical velocity
          to compute the integrand.
        - The integrand involves the fusion reaction cross-section and the critical
          velocity for electron/ion slowing down of the beam ion.

    References:
        - P J Knight, CCFE, Culham Science Centre
    """
    intgeral_term = (velocity_ratio**3) / (1.0 + velocity_ratio**3)

    # critical_velocity : critical velocity for electron/ion slowing down of beam ion (m/s)
    beam_velcoity = critical_velocity * velocity_ratio

    # Calculate the beam kinetic energy per amu and normalise to keV
    xvcs = beam_velcoity**2 * constants.ATOMIC_MASS_UNIT / (constants.KILOELECTRON_VOLT)

    # Calculate the fusion reaction cross-section from beam kinetic energy
    cross_section = _beam_fusion_cross_section(xvcs)

    return intgeral_term * cross_section


def _beam_fusion_cross_section(vrelsq: float) -> float:
    """
    Calculate the fusion reaction cross-section.

    This function computes the fusion reaction cross-section based on the
    square of the speed of the beam ion (keV/amu). The functional form of
    the cross-section is in terms of the equivalent deuterium energy, i.e.,
    for a tritium beam at 500 keV the energy used in the cross-section
    function is 333 keV.

    Parameters:
        vrelsq (float): Square of the speed of the beam ion (keV/amu).

    Returns:
        float: Fusion reaction cross-section (cm^2).

    Notes:
        - The cross-section is limited at low and high beam energies.
        - For beam kinetic energy less than 10 keV, the cross-section is set to 1.0e-27 cm^2.
        - For beam kinetic energy greater than 10,000 keV, the cross-section is set to 8.0e-26 cm^2.
        - The cross-section is calculated using a functional form with parameters a1 to a5.

    References:
        - None
    """
    a1 = 45.95
    a2 = 5.02e4
    a3 = 1.368e-2
    a4 = 1.076
    a5 = 4.09e2

    # Beam kinetic energy
    e_beam_kev = 0.5 * constants.M_DEUTERON_AMU * vrelsq

    # Set limits on cross-section at low and high beam energies
    if e_beam_kev < 10.0:
        return 1.0e-27
    if e_beam_kev > 1.0e4:
        return 8.0e-26
    t1 = a2 / (1.0 + (a3 * e_beam_kev - a4) ** 2) + a5
    t2 = e_beam_kev * (np.exp(a1 / np.sqrt(e_beam_kev)) - 1.0)
    return 1.0e-24 * t1 / t2
