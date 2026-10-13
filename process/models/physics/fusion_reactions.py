"""Fusion reaction physics models and constants."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Self

import numpy as np
from scipy import integrate

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.models.physics.profiles import calculate_vol_avg_of_profile
from process.data_structure.physics_variables import PhysicsData
from process.models.physics.plasma_equilibrium import PlasmaEquilibrium

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from process.models.physics.plasma_profiles import PlasmaProfile

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

@dataclass
class FusionYieldDensities:
    """Fusion power and reaction rate densities produced by a reaction (or their sum).

    Attributes
    ----------
    pden_plasma_alpha_vol_avg_mw :
        Alpha particle fusion power per unit volume [MW/m³].
    pden_non_alpha_charged_vol_avg_mw :
        Non-alpha charged particle fusion power per unit volume [MW/m³].
    pden_neutron_vol_avg_mw :
        Neutron fusion power per unit volume [MW/m³].
    fusden_vol_avg :
        Fusion reaction rate per unit volume [reactions/m³/s].
    fusden_plasma_alpha_vol_avg :
        Alpha particle production rate per unit volume [particles/m³/s].
    fusden_plasma_protons_vol_avg :
        Proton production rate per unit volume [particles/m³/s].


    """

    pden_plasma_alpha_vol_avg_mw: float = 0.0
    pden_non_alpha_charged_vol_avg_mw: float = 0.0
    pden_neutron_vol_avg_mw: float = 0.0
    fusden_vol_avg: float = 0.0
    fusden_plasma_alpha_vol_avg: float = 0.0
    fusden_plasma_protons_vol_avg: float = 0.0

    def __iadd__(self, other: FusionYieldDensities) -> Self:
        """In-place addition of fusion rate densities.

        Raises
        ------
        TypeError
            If an attempt is made to add an object that is not an instance of
            FusionYieldDensities.
        """
        if not isinstance(other, type(self)):
            raise TypeError(f"other must be an instance of {type(self).__name__}")
        self.pden_plasma_alpha_vol_avg_mw += other.pden_plasma_alpha_vol_avg_mw
        self.pden_non_alpha_charged_vol_avg_mw += other.pden_non_alpha_charged_vol_avg_mw
        self.pden_neutron_vol_avg_mw += other.pden_neutron_vol_avg_mw
        self.fusden_vol_avg += other.fusden_vol_avg
        self.fusden_plasma_alpha_vol_avg += other.fusden_plasma_alpha_vol_avg
        self.fusden_plasma_protons_vol_avg += other.fusden_plasma_protons_vol_avg
        return self


class PlasmaReactions(Model):
    """Calculate the fusion reaction rate for each reaction case (DT, DHE3, DD1, DD2).

    This class provides methods to numerically integrate over the plasma cross-section
    to find the core plasma fusion power for different fusion reactions. The reactions
    considered are:
        - Deuterium-Tritium (D-T)
        - Deuterium-Helium-3 (D-3He)
        - Deuterium-Deuterium (D-D) first branch
        - Deuterium-Deuterium (D-D) second branch

    The class uses the Bosch-Hale parametrization to compute the volumetric fusion
    reaction rate 〈σv〉 and integrates over the plasma cross-section to find the
    core plasma fusion power.

    Attributes
    ----------
    plasma_profile : PlasmaProfile
        Parameterized temperature and density profiles of the plasma.
    sigma_v_plasma_dt_vol_avg : float
        Volume-averaged D-T fusion reactivity, 〈σv〉ᵥ, for D-T.
    pden_dhe3_total_vol_avg_mw : float
        Fusion power density produced by the D-3He reaction.
    pden_dd_total_vol_avg_mw : float
        Fusion power density produced by the D-D reactions.
    dt_power_density : float
        Fusion power density produced by the D-T reaction.
    fusion_rates : FusionYieldDensities
        Cumulative fusion power and reaction rate densities, summed over all
        reactions.
    f_dd_branching_trit : float
        Ratio of tritium-producing D-D reactions to 3He-producing D-D reactions.

    References
    ----------
    [1] H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections
    and thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611-631,
    Apr. 1992, doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """  # noqa: RUF002

    def __init__(self, plasma_profile: PlasmaProfile):
        """
        Initialize the PlasmaReactions class with the given plasma profile.

        Parameters
        ----------
        plasma_profile:
            The parameterized temperature and density profiles of the plasma.


        """
        self.plasma_profile = plasma_profile
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.sigma_v_plasma_dt_vol_avg = 0.0
        self.sigma_v_plasma_pb_vol_avg = 0.0
        self.pden_dhe3_total_vol_avg_mw = 0.0
        self.pden_dd_total_vol_avg_mw = 0.0
        self.dt_power_density = 0.0
        self.pb_power_density = 0.0
        self.fusion_rates = FusionYieldDensities()
        self.f_dd_branching_trit = 0.0

    def run(self) -> None:
        """Run the plasma reactions calculations."""

    def deuterium_branching(self, ion_temperature: float) -> float:
        """Calculate the relative rate of tritium producing D-D reactions to 3He ones
        based on the volume averaged ion temperature

        Parameters
        ----------
        ion_temperature :
            float

        Notes
        -----
        - For ion temperatures between 0.5 keV and 200 keV.
        - The deviation of the fit from the R-matrix branching ratio is always smaller
        than 0.5%.

        References
        ----------
        [1] H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections
        and thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611-631,
        Apr. 1992, doi: https://doi.org/10.1088/0029-5515/32/4/i07.
        """
        # Divide by 2 to get the branching ratio for the D-D reaction that produces
        # tritium as the output is just the ratio of the two normalized cross sections
        self.f_dd_branching_trit = (
            1.02934
            - 8.3264e-3 * ion_temperature
            + 1.7631e-4 * ion_temperature**2
            - 1.8201e-6 * ion_temperature**3
            + 6.9855e-9 * ion_temperature**4
        ) / 2.0

    def dt_reaction(self):
        """D + T --> 4He + n reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-tritium (D-T) fusion reaction. It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate 〈σv〉 and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.sigma_v_plasma_dt_vol_avg: Volume averaged D-T fusion
              reactivity 〈σv〉ᵥ for D-T.
            - self.dt_power_density: Fusion power density produced by the D-T reaction.
            - self.fusion_rates: Cumulative fusion power and reaction rate densities.

        """  # noqa: RUF002
        # Initialize Bosch-Hale constants for the D-T reaction
        dt = BoschHaleConstants(**REACTION_CONSTANTS_DT)

        # Velocity-space average reactivity D-T profile [m³/s]
        sigma_v_profile = bosch_hale_reactivity(
            ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            reaction_constants=dt,
        )

        self.sigma_v_plasma_dt_vol_avg = calculate_vol_avg_of_profile(
            profile_x=self.plasma_profile.neprofile.profile_x,
            profile_y=sigma_v_profile,
            profile_dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Full fusion reaction rate density profile for the D-T reaction [reactions/m³/s]
        self.data.physics.fusden_plasma_dt_profile = fusden_plasma_dt_profile = (
            sigma_v_profile
            * self.data.physics.f_plasma_fuel_deuterium
            * self.data.physics.f_plasma_fuel_tritium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    self.data.physics.nd_plasma_fuel_ions_vol_avg
                    / self.data.physics.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Volume-averaged fusion reaction rate for the D-T reaction [reactions/m³/s]
        fusden_dt_vol_avg = calculate_vol_avg_of_profile(
            profile_x=self.plasma_profile.neprofile.profile_x,
            profile_y=fusden_plasma_dt_profile,
            profile_dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Reaction energy in MegaJoules [MJ]
        e_reaction_mj = constants.D_T_ENERGY / 1.0e6

        # Volume-averaged power density for the D-T reaction [MW/m³]
        pden_dt_vol_avg_mw = fusden_dt_vol_avg * e_reaction_mj

        # Power densities for different particles [MW/m³]
        # Alpha particle gets approximately 20% of the fusion power
        pden_alpha_vol_avg_mw = (
            1.0 - constants.DT_NEUTRON_ENERGY_FRACTION
        ) * pden_dt_vol_avg_mw
        pden_neutron_vol_avg_mw = (
            constants.DT_NEUTRON_ENERGY_FRACTION * pden_dt_vol_avg_mw
        )

        # Update the cumulative D-T power density
        self.dt_power_density = pden_dt_vol_avg_mw

        # Sum the fusion rates for all particles
        self.fusion_rates += FusionYieldDensities(
            pden_plasma_alpha_vol_avg_mw=pden_alpha_vol_avg_mw,
            pden_neutron_vol_avg_mw=pden_neutron_vol_avg_mw,
            fusden_vol_avg=fusden_dt_vol_avg,
            fusden_plasma_alpha_vol_avg=fusden_dt_vol_avg,
        )

    def dhe3_reaction(self):
        """D + 3He --> 4He + p reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-helium-3 (D-3He) fusion reaction. It uses the Bosch-Hale
        parametrization to compute the volumetric fusion reaction rate 〈σv〉 and
        integrates over the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.pden_dhe3_total_vol_avg_mw: Fusion power density produced by the D-3He
              reaction.
            - self.fusion_rates: Cumulative fusion power and reaction rate densities.
        """  # noqa: RUF002
        # Initialize Bosch-Hale constants for the D-3He reaction
        dhe3 = BoschHaleConstants(**REACTION_CONSTANTS_DHE3)

        # Velocity-space average reactivity D-3He profile [m³/s]
        sigma_v_profile = bosch_hale_reactivity(
            ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            reaction_constants=dhe3,
        )

        # Full fusion reaction rate profile for the D-3He reaction [reactions/m³/s]
        self.data.physics.fusden_plasma_dhe3_profile = fusden_plasma_dhe3_profile = (
            sigma_v_profile
            * self.data.physics.f_plasma_fuel_deuterium
            * self.data.physics.f_plasma_fuel_helium3
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    self.data.physics.nd_plasma_fuel_ions_vol_avg
                    / self.data.physics.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Volume-averaged fusion reaction rate for the D-3He reaction [reactions/m³/s]
        fusden_dhe3_vol_avg = calculate_vol_avg_of_profile(
            profile_x=self.plasma_profile.neprofile.profile_x,
            profile_y=fusden_plasma_dhe3_profile,
            profile_dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Reaction energy in MegaJoules [MJ]
        e_reaction_mj = constants.D_HELIUM_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m³]
        pden_dhe3_vol_avg_mw = fusden_dhe3_vol_avg * e_reaction_mj

        # Power densities for different particles [MW/m³]
        # Alpha particle gets approximately 20% of the fusion power
        pden_alpha_vol_avg_mw = (
            1.0 - constants.DHELIUM_PROTON_ENERGY_FRACTION
        ) * pden_dhe3_vol_avg_mw
        pden_non_alpha_charged_mw = (
            constants.DHELIUM_PROTON_ENERGY_FRACTION * pden_dhe3_vol_avg_mw
        )

        # Update the cumulative D-3He power density
        self.pden_dhe3_total_vol_avg_mw = pden_dhe3_vol_avg_mw

        # Sum the fusion rates for all particles
        self.fusion_rates += FusionYieldDensities(
            pden_plasma_alpha_vol_avg_mw=pden_alpha_vol_avg_mw,
            pden_non_alpha_charged_vol_avg_mw=pden_non_alpha_charged_mw,
            fusden_vol_avg=fusden_dhe3_vol_avg,
            fusden_plasma_alpha_vol_avg=fusden_dhe3_vol_avg,
            fusden_plasma_protons_vol_avg=fusden_dhe3_vol_avg,
        )

    def dd_helion_reaction(self):
        """D + D --> 3He + n reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces
        helium-3 (3He) and a neutron (n). It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate 〈σv〉 and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
        - self.pden_dd_total_vol_avg_mw: Fusion power density produced by the
          D-D reaction.
        - self.fusion_rates: Cumulative fusion power and reaction rate densities.
        """  # noqa: RUF002
        # Initialize Bosch-Hale constants for the D-D reaction
        dd1 = BoschHaleConstants(**REACTION_CONSTANTS_DD1)

        # Velocity-space average reactivity D-D -> 3He + n profile [m³/s]
        sigma_v_profile = bosch_hale_reactivity(
            ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            reaction_constants=dd1,
        )

        # Full fusion reaction rate density profile for the D-D -> 3He + n
        # reaction [reactions/m³/s]
        self.data.physics.fusden_plasma_dd_helion_profile = (
            fusden_plasma_dd_helion_profile
        ) = (
            sigma_v_profile
            * self.data.physics.f_plasma_fuel_deuterium
            * self.data.physics.f_plasma_fuel_deuterium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    self.data.physics.nd_plasma_fuel_ions_vol_avg
                    / self.data.physics.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Volume-averaged fusion reaction rate for the D-D -> 3He + n
        # reaction [reactions/m³/s]
        fusden_dd_helion_vol_avg = calculate_vol_avg_of_profile(
            profile_x=self.plasma_profile.neprofile.profile_x,
            profile_y=fusden_plasma_dd_helion_profile,
            profile_dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Reaction energy in MegaJoules [MJ]
        e_reaction_mj = constants.DD_HELIUM_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m³]
        # The power density is scaled by the branching ratio to simulate the different
        # product pathways
        pden_dd_helion_vol_avg_mw = (
            fusden_dd_helion_vol_avg * e_reaction_mj * (1.0 - self.f_dd_branching_trit)
        )

        # Power densities for different particles [MW/m³]
        # Neutron particle gets approximately 75% of the fusion power
        pden_non_alpha_charged_vol_avg_mw = (
            1.0 - constants.DD_NEUTRON_ENERGY_FRACTION
        ) * pden_dd_helion_vol_avg_mw
        pden_neutron_vol_avg_mw = (
            constants.DD_NEUTRON_ENERGY_FRACTION * pden_dd_helion_vol_avg_mw
        )

        # Update the cumulative D-D power density
        self.pden_dd_total_vol_avg_mw += pden_dd_helion_vol_avg_mw

        # Sum the fusion rates for all particles
        self.fusion_rates += FusionYieldDensities(
            pden_non_alpha_charged_vol_avg_mw=pden_non_alpha_charged_vol_avg_mw,
            pden_neutron_vol_avg_mw=pden_neutron_vol_avg_mw,
            fusden_vol_avg=pden_dd_helion_vol_avg_mw,
        )

    def dd_triton_reaction(self):
        """D + D --> T + p reaction

        This method calculates the fusion reaction rate and power density for the
        deuterium-deuterium (D-D) fusion reaction, specifically the branch that produces
        tritium (T) and a proton (p). It uses the Bosch-Hale parametrization
        to compute the volumetric fusion reaction rate 〈σv〉 and integrates over
        the plasma cross-section to find the core plasma fusion power.

        The method updates the following attributes:
            - self.pden_dd_total_vol_avg_mw: Fusion power density produced by the
              D-D reaction.
            - self.fusion_rates: Cumulative fusion power and reaction rate densities.
        """  # noqa: RUF002
        # Initialize Bosch-Hale constants for the D-D reaction
        dd2 = BoschHaleConstants(**REACTION_CONSTANTS_DD2)

        # Velocity-space average reactivity for the D-D -> T + p reaction
        sigma_v_profile = bosch_hale_reactivity(
            ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            reaction_constants=dd2,
        )

        # Full fusion reaction rate density profile for the D-D -> T + p reaction
        # [reactions/m³/s]
        self.data.physics.fusden_plasma_dd_triton_profile = (
            fusden_plasma_dd_triton_profile
        ) = (
            sigma_v_profile
            * self.data.physics.f_plasma_fuel_deuterium
            * self.data.physics.f_plasma_fuel_deuterium
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    self.data.physics.nd_plasma_fuel_ions_vol_avg
                    / self.data.physics.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Volume-averaged fusion reaction rate density for the D-D -> T + p
        # reaction [reactions/m³/s]
        fusden_dd_triton_vol_avg = calculate_vol_avg_of_profile(
            profile_x=self.plasma_profile.neprofile.profile_x,
            profile_y=fusden_plasma_dd_triton_profile,
            profile_dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Reaction energy in MegaJoules [MJ]
        e_reaction_mj = constants.DD_TRITON_ENERGY / 1.0e6

        # Calculate the fusion power density produced [MW/m³]
        # The power density is scaled by the branching ratio to simulate the different
        # product pathways
        pden_dd_triton_vol_avg_mw = (
            fusden_dd_triton_vol_avg * e_reaction_mj * self.f_dd_branching_trit
        )

        # Calculate the fusion rate density [reactions/m³/s]
        fusion_rate_density = pden_dd_triton_vol_avg_mw / e_reaction_mj

        # Proton production rate [particles/m³/s]
        fusden_plasma_protons_vol_avg = fusion_rate_density

        # Update the cumulative D-D power density
        self.pden_dd_total_vol_avg_mw += pden_dd_triton_vol_avg_mw

        # Sum the fusion rates for all particles
        self.fusion_rates += FusionYieldDensities(
            pden_non_alpha_charged_vol_avg_mw=pden_dd_triton_vol_avg_mw,
            fusden_vol_avg=fusion_rate_density,
            fusden_plasma_protons_vol_avg=fusden_plasma_protons_vol_avg,
        )

    def pb_reaction(self):
        """PB reaction
        
        This method calculates the fusion reaction rate and power density for the
        pb fusion reaction.
        
        """
        if self.data.physics.i_pb_reactivity_model == 1:
            sigma_v_profile = nevins_reactivity(
                ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            )
        elif self.data.physics.i_pb_reactivity_model == 2:
            sigma_v_profile = tentori_belloni_reactivity(
                ion_temperature_profile=self.plasma_profile.tiprofile.profile_y,
            )
        else:
            raise ValueError(f"Invalid PB reaction model: {self.data.physics.i_pb_reactivity_model}")

        self.data.physics.fusden_plasma_pb_profile = (
            sigma_v_profile
            * self.data.physics.f_plasma_fuel_proton
            * self.data.physics.f_plasma_fuel_boron11
            * (
                self.plasma_profile.neprofile.profile_y
                * (
                    self.data.physics.nd_plasma_fuel_ions_vol_avg
                    / self.data.physics.nd_plasma_electrons_vol_avg
                )
            )
            ** 2
        )

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.PB_ENERGY / 1.0e6

        if self.data.physics.i_equilibrium_solve == 1 and self.data.veqpy.equilibrium is not None:
            eq = self.data.veqpy.equilibrium
            self.sigma_v_plasma_pb_vol_avg = PlasmaEquilibrium.get_volume_average(
                eq,
                self.plasma_profile.neprofile.profile_x, 
                sigma_v_profile,
            )
            fusden_pb_vol_avg = PlasmaEquilibrium.get_volume_average(
                eq,
                self.plasma_profile.neprofile.profile_x, 
                self.data.physics.fusden_plasma_pb_profile,
            )
        else:
            self.sigma_v_plasma_pb_vol_avg = calculate_vol_avg_of_profile(
                profile_x=self.plasma_profile.neprofile.profile_x,
                profile_y=sigma_v_profile,
                profile_dx=self.plasma_profile.neprofile.profile_dx,
            )
            fusden_pb_vol_avg = calculate_vol_avg_of_profile(
                profile_x=self.plasma_profile.neprofile.profile_x,
                profile_y=self.data.physics.fusden_plasma_pb_profile,
                profile_dx=self.plasma_profile.neprofile.profile_dx,
            )

        pb_power_density = fusden_pb_vol_avg * reaction_energy
        self.pb_power_density = pb_power_density

        self.fusion_rates += FusionYieldDensities(
            pden_plasma_alpha_vol_avg_mw=pb_power_density,
            fusden_vol_avg=fusden_pb_vol_avg,
            fusden_plasma_alpha_vol_avg=fusden_pb_vol_avg * 3.0,
        )


    def calculate_fusion_rates(self):
        """Initiate all the fusion rate calculations.

        This method sequentially calculates the fusion reaction rates and power
        densities for the following reactions:
            - Deuterium-Tritium (D-T)
            - Deuterium-Helium-3 (D-3He)
            - Deuterium-Deuterium (D-D) first branch
            - Deuterium-Deuterium (D-D) second branch
            - Proton-Boron-11 (p-B11), when ``i_fusion_reactions == "p-b11"``

        It updates the instance attributes for the cumulative power densities and
        reaction rates
        for alpha particles, charged particles, neutrons, and protons.


        """
        # Reset cumulative attributes as this instance persists across calls
        self.fusion_rates = FusionYieldDensities()
        self.pden_dd_total_vol_avg_mw = 0.0

        self.dt_reaction()
        self.dhe3_reaction()
        self.dd_helion_reaction()
        self.dd_triton_reaction()
        if self.data.physics.i_fusion_reactions == "p-b11":
            self.pb_reaction()

    def set_physics_variables(self):
        """Set the required physics variables in the physics_variables and
        physics modules.

        This method updates the global physics variables and module variables with the
        current instance's fusion power densities and reaction rates.


        """
        self.data.physics.pden_plasma_alpha_vol_avg_mw = (
            self.fusion_rates.pden_plasma_alpha_vol_avg_mw
        )
        self.data.physics.pden_non_alpha_charged_mw = (
            self.fusion_rates.pden_non_alpha_charged_vol_avg_mw
        )
        self.data.physics.pden_plasma_neutron_vol_avg_mw = (
            self.fusion_rates.pden_neutron_vol_avg_mw
        )
        self.data.physics.fusden_plasma_vol_avg = self.fusion_rates.fusden_vol_avg
        self.data.physics.fusden_plasma_alpha_vol_avg = (
            self.fusion_rates.fusden_plasma_alpha_vol_avg
        )
        self.data.physics.fusden_plasma_protons_vol_avg = (
            self.fusion_rates.fusden_plasma_protons_vol_avg
        )
        self.data.physics.sigma_v_plasma_dt_vol_avg = self.sigma_v_plasma_dt_vol_avg
        self.data.physics.pden_plasma_dt_vol_avg_mw = self.dt_power_density
        self.data.physics.pden_dhe3_total_vol_avg_mw = self.pden_dhe3_total_vol_avg_mw
        self.data.physics.pden_dd_total_vol_avg_mw = self.pden_dd_total_vol_avg_mw
        self.data.physics.f_dd_branching_trit = self.f_dd_branching_trit
        self.data.physics.sigma_v_plasma_pb_vol_avg = self.sigma_v_plasma_pb_vol_avg
        self.data.physics.pden_plasma_pb_vol_avg_mw = self.pb_power_density

    def output(self):
        """Output fusion reaction information."""
        po.oheadr(self.outfile, "Plasma Reactions :")

        po.osubhd(self.outfile, "Fuel Constituents :")
        po.ovarre(
            self.outfile,
            "Deuterium fuel fraction",
            "(f_plasma_fuel_deuterium)",
            self.data.physics.f_plasma_fuel_deuterium,
        )
        po.ovarre(
            self.outfile,
            "Tritium fuel fraction",
            "(f_plasma_fuel_tritium)",
            self.data.physics.f_plasma_fuel_tritium,
        )
        po.ovarre(
            self.outfile,
            "3-Helium fuel fraction",
            "(f_plasma_fuel_helium3)",
            self.data.physics.f_plasma_fuel_helium3,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion rates :")
        po.ovarre(
            self.outfile,
            "Fusion rate: total [reactions/sec]",
            "(fusrat_total)",
            self.data.physics.fusrat_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged fusion rate density: total [reactions/m³/sec]",
            "(fusden_total_vol_avg)",
            self.data.physics.fusden_total_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged fusion rate density: plasma [reactions/m³/sec]",
            "(fusden_plasma_vol_avg)",
            self.data.physics.fusden_plasma_vol_avg,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion Powers :")
        po.ocmmnt(
            self.outfile,
            "Fusion power totals from the main plasma and beam-plasma interactions "
            "(if present)\n",
        )

        po.ovarre(
            self.outfile,
            "Total fusion power [MW]",
            "(p_fusion_total_mw)",
            self.data.physics.p_fusion_total_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "D-T fusion power: total [MW]",
            "(p_dt_total_mw)",
            self.data.physics.p_dt_total_mw,
            "OP ",
        )
        for i in range(len(self.data.physics.fusden_plasma_dt_profile)):
            po.ovarre(
                self.mfile,
                f"DT fusion rate density in plasma at point {i}",
                f"fusden_plasma_dt_profile{i}",
                self.data.physics.fusden_plasma_dt_profile[i],
            )

        for i in range(len(self.data.physics.fusden_plasma_dd_triton_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> T fusion rate density in plasma at point {i}",
                f"fusden_plasma_dd_triton_profile{i}",
                self.data.physics.fusden_plasma_dd_triton_profile[i],
            )
        for i in range(len(self.data.physics.fusden_plasma_dd_helion_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> 3He fusion rate density in plasma at point {i}",
                f"fusden_plasma_dd_helion_profile{i}",
                self.data.physics.fusden_plasma_dd_helion_profile[i],
            )
        for i in range(len(self.data.physics.fusden_plasma_dhe3_profile)):
            po.ovarre(
                self.mfile,
                f"D-3He fusion rate density in plasma at point {i}",
                f"fusden_plasma_dhe3_profile{i}",
                self.data.physics.fusden_plasma_dhe3_profile[i],
            )
        po.ovarre(
            self.outfile,
            "D-T fusion power: plasma [MW]",
            "(p_plasma_dt_mw)",
            self.data.physics.p_plasma_dt_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged D-T fusion power density: plasma [MW/m³]",
            "(pden_plasma_dt_vol_avg_mw)",
            self.data.physics.pden_plasma_dt_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged plasma D-T fusion reactivity (〈σv〉ᵥ) [m³/s]",  # noqa: RUF001
            "(sigma_v_plasma_dt_vol_avg)",
            self.data.physics.sigma_v_plasma_dt_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-T fusion power: beam [MW]",
            "(p_beam_dt_mw)",
            self.data.physics.p_beam_dt_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "D-D fusion power [MW]",
            "(p_dd_total_mw)",
            self.data.physics.p_dd_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged D-D fusion power density [MW/m³]",
            "(pden_dd_total_vol_avg_mw)",
            self.data.physics.pden_dd_total_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-D branching ratio for tritium producing reactions",
            "(f_dd_branching_trit)",
            self.data.physics.f_dd_branching_trit,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "D-He3 fusion power [MW]",
            "(p_dhe3_total_mw)",
            self.data.physics.p_dhe3_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged D-3He fusion power density [MW/m³]",
            "(pden_dhe3_total_vol_avg_mw)",
            self.data.physics.pden_dhe3_total_vol_avg_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Alpha Powers (α) :")  # noqa: RUF001
        po.ovarre(
            self.outfile,
            "Volume-averaged alpha rate density: total [particles/m³/sec]",
            "(fusden_alpha_total_vol_avg)",
            self.data.physics.fusden_alpha_total_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged alpha rate density: plasma [particles/m³/sec]",
            "(fusden_plasma_alpha_vol_avg)",
            self.data.physics.fusden_plasma_alpha_vol_avg,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Alpha power: total [MW]",
            "(p_alpha_total_mw)",
            self.data.physics.p_alpha_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged alpha power density: total [MW/m³]",
            "(pden_alpha_total_vol_avg_mw)",
            self.data.physics.pden_alpha_total_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: plasma only [MW]",
            "(p_plasma_alpha_mw)",
            self.data.physics.p_plasma_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged alpha power density: plasma [MW/m³]",
            "(pden_plasma_alpha_vol_avg_mw)",
            self.data.physics.pden_plasma_alpha_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: beam-plasma [MW]",
            "(p_beam_alpha_mw)",
            self.data.physics.p_beam_alpha_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to electrons [MW/m³]",
            "(f_pden_alpha_electron_mw)",
            self.data.physics.f_pden_alpha_electron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to ions [MW/m³]",
            "(f_pden_alpha_ions_mw)",
            self.data.physics.f_pden_alpha_ions_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Neutron Powers (n) :")
        po.ovarre(
            self.outfile,
            "Neutron power: total [MW]",
            "(p_neutron_total_mw)",
            self.data.physics.p_neutron_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged neutron power density: total [MW/m³]",
            "(pden_neutron_total_vol_avg_mw)",
            self.data.physics.pden_neutron_total_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: plasma only [MW]",
            "(p_plasma_neutron_mw)",
            self.data.physics.p_plasma_neutron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged neutron power density: plasma [MW/m³]",
            "(pden_plasma_neutron_vol_avg_mw)",
            self.data.physics.pden_plasma_neutron_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: beam-plasma [MW]",
            "(p_beam_neutron_mw)",
            self.data.physics.p_beam_neutron_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Average neutron flux at plasma surface [MW/m²]",
            "(pflux_plasma_surface_neutron_avg_mw)",
            self.data.physics.pflux_plasma_surface_neutron_avg_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")

        po.osubhd(self.outfile, "Charged Particle Powers :")

        po.ovarre(
            self.outfile,
            "Charged particle power (p, 3He, T) (excluding alphas) [MW]",
            "(p_non_alpha_charged_mw)",
            self.data.physics.p_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Charged particle power density (p, 3He, T) (excluding alphas) [MW]",
            "(pden_non_alpha_charged_mw)",
            self.data.physics.pden_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total charged particle power (including alphas) [MW]",
            "(p_charged_particle_mw)",
            self.data.physics.p_charged_particle_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Volume-averaged proton rate density [particles/m³/sec]",
            "(fusden_plasma_protons_vol_avg)",
            self.data.physics.fusden_plasma_protons_vol_avg,
            "OP ",
        )
        po.oblnkl(self.outfile)


class BeamReactions(Model):
    """Calculate neutral beam slowing-down and beam-target fusion power.

    This class calculates the beam-target fusion contribution from neutral beam
    injection, complementing the thermal plasma reactions calculated by
    PlasmaReactions.

    Attributes
    ----------
    beta_beam : float
        Neutral beam beta contribution.
    nd_beam_ions_out : float
        Hot beam ion density [m⁻³].
    p_beam_alpha_mw : float
        Alpha power from beam-target fusion [MW].
    p_beam_neutron_mw : float
        Neutron power from beam-target fusion [MW].
    p_beam_dt_mw : float
        D-T fusion power from beam-target fusion [MW].
    """

    def __init__(self):

        self.beta_beam = 0.0
        self.nd_beam_ions_out = 0.0
        self.p_beam_alpha_mw = 0.0
        self.p_beam_neutron_mw = 0.0
        self.p_beam_dt_mw = 0.0
        self.p_beam_pb_mw = 0.0

    def run(self) -> None:
        """Run the beam reactions calculations."""

    def output(self) -> None:
        """Output the beam reactions results."""

    def calculate_beam_fusion(self):
        """Calculate neutral beam slowing-down and beam-target fusion power.

        Neutral beam fusion is neglected, leaving beam attributes at zero, if
        there is no beam current or the plasma is ignited (beams cannot be
        present in an ignited plasma).

        The method updates the following attributes:
            - self.beta_beam: Neutral beam beta contribution.
            - self.nd_beam_ions_out: Hot beam ion density [m⁻³].
            - self.p_beam_alpha_mw: Alpha power from beam-target fusion [MW].
            - self.p_beam_neutron_mw: Neutron power from beam-target fusion [MW].
            - self.p_beam_dt_mw: D-T fusion power from beam-target fusion [MW].
        """
        # Reset as this instance persists across calls
        self.beta_beam = 0.0
        self.nd_beam_ions_out = 0.0
        self.p_beam_alpha_mw = 0.0
        self.p_beam_neutron_mw = 0.0
        self.p_beam_dt_mw = 0.0
        self.p_beam_pb_mw = 0.0

        self.beta_beam, self.nd_beam_ions_out, self.p_beam_alpha_mw = beam_fusion(
            self.data.physics.beamfus0,
            self.data.physics.betbm0,
            self.data.physics.b_plasma_total,
            self.data.current_drive.c_beam_total,
            self.data.physics.nd_plasma_electrons_vol_avg,
            self.data.physics.nd_plasma_fuel_ions_vol_avg,
            self.data.physics.dlamie,
            self.data.current_drive.e_beam_kev,
            self.data.physics.f_plasma_fuel_deuterium,
            self.data.physics.f_plasma_fuel_tritium,
            self.data.current_drive.f_beam_tritium,
            self.data.physics.temp_plasma_electron_density_weighted_kev,
            self.data.physics.vol_plasma,
            self.data.physics.n_charge_plasma_effective_mass_weighted_vol_avg,
        )

        self.p_beam_neutron_mw = self.p_beam_alpha_mw * (
            constants.DT_NEUTRON_ENERGY_FRACTION
            / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION)
        )
        self.p_beam_dt_mw = self.p_beam_alpha_mw * (
            1.0 / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION)
        )
 
        if self.data.physics.i_fusion_reactions == "p-b11":
            (
                self.beta_beam,
                self.nd_beam_ions_out,
                self.p_beam_alpha_mw,
            ) = beam_fusion_p_b11(
                sigmv_pb=self.data.physics.sigma_v_plasma_pb_vol_avg,
                i_pb_reactivity_model=self.data.physics.i_pb_reactivity_model,
                ti=self.data.physics.temp_plasma_ion_vol_avg_kev,
                te=self.data.physics.temp_plasma_electron_vol_avg_kev,
                vol_plasma=self.data.physics.vol_plasma,
                current_beam_total=self.data.current_drive.c_beam_total,
                b_total=self.data.physics.b_plasma_total,
                ne=self.data.physics.nd_plasma_electrons_vol_avg,
                nd_fuel_b11=self.data.physics.nd_plasma_fuel_ions_vol_avg * self.data.physics.f_plasma_fuel_boron11,
                charge_eff_mass_weighted=self.data.physics.n_charge_plasma_effective_mass_weighted_vol_avg,
                dlamie=self.data.physics.dlamie,
                dlamee=self.data.physics.dlamee,
                betbm0=self.data.physics.betbm0,
                beamfus0=self.data.physics.beamfus0,
                e_beam_kev=self.data.current_drive.e_beam_kev,
            )
            self.p_beam_neutron_mw = 0.0
            self.p_beam_dt_mw = 0.0
            self.p_beam_pb_mw = self.p_beam_alpha_mw

    def set_physics_variables(self):
        """Set the beam fusion physics variables in the physics module."""
        self.data.physics.beta_beam = self.beta_beam
        self.data.physics.nd_beam_ions_out = self.nd_beam_ions_out
        self.data.physics.p_beam_alpha_mw = self.p_beam_alpha_mw
        self.data.physics.p_beam_neutron_mw = self.p_beam_neutron_mw
        self.data.physics.p_beam_dt_mw = self.p_beam_dt_mw
        self.data.physics.p_beam_pb_mw = self.p_beam_pb_mw


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



def bosch_hale_reactivity(
    ion_temperature_profile: np.ndarray, reaction_constants: BoschHaleConstants
) -> np.ndarray:
    """Calculate the velocity-space average reactivity 〈σv〉 [m³/s] for one of four
    nuclear reactions using the Bosch-Hale parametrization.

    The valid range of the fit is 0.2 keV < t < 100 keV except for D-3He where it is
    0.5 keV < t < 190 keV.

    Reactions:
        1. D-T reaction
        2. D-3He reaction
        3. D-D 1st reaction
        4. D-D 2nd reaction

    Parameters
    ----------
    ion_temperature_profile :
        Plasma ion temperature profile in keV.
    reaction_constants :
        Bosch-Hale reaction constants.

    Returns
    -------
    :
        np.ndarray: Velocity-space average reactivity 〈σv〉 [m³/s] for each point in
        the ion temperature profile.

    References
    ----------
    [1] H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and
    thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
    doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """  # noqa: RUF002
    theta1 = (
        ion_temperature_profile
        * (
            reaction_constants.cc2
            + ion_temperature_profile
            * (reaction_constants.cc4 + ion_temperature_profile * reaction_constants.cc6)
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

    # Volumetric reaction rate / reactivity 〈σv〉 [m³/s] # noqa: RUF003
    # Original form is in [cm³/s], so multiply by 1.0e-6 to convert to [m³/s]
    sigmav = (
        1.0e-6
        * reaction_constants.cc1
        * theta
        * np.sqrt(xi / (reaction_constants.mrc2 * ion_temperature_profile**3))
        * np.exp(-3.0 * xi)
    )

    # if t = 0, sigmav = 0. Use this mask to set sigmav to zero.
    t_mask = ion_temperature_profile == 0.0  # noqa: RUF069
    sigmav[t_mask] = 0.0

    # Return np.ndarray of sigmav for each point in the ion temperature profile
    return sigmav

def fsgmpb(e):
    """
    Calculate fusion cross-section for p + B11 -> 3He4
    Ref: Nevins 2000
    """
    bgpb = np.sqrt(22.589e3)

    # ind1 = np.argmin(np.abs(e - 400))
    # ind2 = np.argmin(np.abs(e - 642))
    if isinstance(e, float):
        if e <= 400:
            spb = 1.97e5 + 0.24e3 * e + 2.31e-1 * e**2 + 1.82e7 / ((e - 148.0)**2 + 2.35**2)
        elif e <= 642:
            e1 = e / 100.0 - 4.0
            spb = 3.30e5 + 66.1e3 * e1 - 20.3e3 * e1**2 - 1.58e3 * e1**5
        else:
            spb = (4.38e3 + 2.57e9 / ((e - 581.3)**2 + 85.7**2) +
                    5.67e8 / ((e - 1083.0)**2 + 234.0**2) +
                    1.34e8 / ((e - 2405.0)**2 + 138.0**2) +
                    5.68e8 / ((e - 3344.0)**2 + 309.0**2))
    else:
        ind1 = np.where(e <= 400)[0]
        ind2 = np.where((e <= 642) & (e > 400))[0]
        ind3 = np.where(e > 642)[0]

        e1 = e[ind1]
        # mid-band polynomial uses ξ = E/100 − 4 (same as scalar branch / Nevins 2000)
        e2 = e[ind2] / 100.0 - 4.0
        e3 = e[ind3]

        spb1 = 1.97e5 + 0.24e3 * e1 + 2.31e-1 * e1**2 + 1.82e7 / ((e1 - 148.0)**2 + 2.35**2)
        spb2 = 3.30e5 + 66.1e3 * e2 - 20.3e3 * e2**2 - 1.58e3 * e2**5
        spb3 = (4.38e3 + 2.57e9 / ((e3 - 581.3)**2 + 85.7**2) +
                5.67e8 / ((e3 - 1083.0)**2 + 234.0**2) +
                1.34e8 / ((e3 - 2405.0)**2 + 138.0**2) +
                5.68e8 / ((e3 - 3344.0)**2 + 309.0**2))

        spb = np.concatenate([spb1, spb2, spb3])
    sigmapb = spb / e * np.exp(-bgpb / np.sqrt(e)) * 1e-28
    return sigmapb

def nevins_reactivity(
    ion_temperature_profile: np.ndarray,
) -> np.ndarray:
    """Calculate the volumetric fusion reaction rate 〈sigmav〉 (m³/s) for 
    p-B reaction using the Nevins parametrization.
    """
    mp = 1.007276466621 * 1.66053906660e-27
    mb = 11 * mp
    mr = (mp * mb) / (mp + mb) 

    e = 10 ** np.arange(-1, 3.6005, 0.0005)  # keV

    sigma = fsgmpb(e)

    # analytic integral
    sgmv = np.zeros_like(ion_temperature_profile)
    if isinstance(ion_temperature_profile, float):
        len_teff = 1
    else:
        len_teff = len(ion_temperature_profile)
    for j in range(len_teff):
        if len_teff == 1:
            teff_now = ion_temperature_profile
        else:
            teff_now = ion_temperature_profile[j]
        sgmv_now = (
            np.sqrt(constants.KILOELECTRON_VOLT * teff_now * 8 / (np.pi * mr))
            / (teff_now ** 2)
            * np.nansum(e[1:] * sigma[1:] * np.diff(e) * np.exp(-e[1:] / teff_now))
        )
        if len_teff == 1:
            sgmv = sgmv_now
        else:
            sgmv[j] = sgmv_now
    return sgmv

def fsgmpb_tentori_belloni(e):
    c0 = 197.0
    c1 = 0.269
    c2 = 2.54e-4
    d0 = 346.0
    d1 = 150.0
    d2 = -59.9
    d5 = -0.460
    a0 = 1.98e6
    a1 = 3.89e6
    a2 = 1.36e6
    a3 = 3.71e6
    e0 = 640.9
    e1 = 1211.0
    e2 = 2340.0
    e3 = 3294.0
    de0 = 85.5
    de1 = 414.0
    de2 = 221
    de3 = 351
    b = 0.381
    if e <= 0.400:
        se = c0 + c1 * e * 1.0e3 + c2 * (e * 1.0e3)**2
    elif e <= 0.668:
        e_eff = (e - 0.400) * 10.0
        se = d0 + d1 * e_eff + d2 * e_eff**2 + d5 * e_eff**5
    elif e <= 9.76:
        e_eff = e * 1.0e3
        se = (
            b + a0 / ((e_eff - e0)**2 + de0**2)
            + a1 / ((e_eff - e1)**2 + de1**2)
            + a2 / ((e_eff - e2)**2 + de2**2)
            + a3 / ((e_eff - e3)**2 + de3**2)
        )
    else:
        e_eff = 9.76 * 1.0e3
        se = (
            b + a0 / ((e_eff - e0)**2 + de0**2)
            + a1 / ((e_eff - e1)**2 + de1**2)
            + a2 / ((e_eff - e2)**2 + de2**2)
            + a3 / ((e_eff - e3)**2 + de3**2)
        )
    return se / e * np.exp(-np.sqrt(22.589 / e)) * 1.0e-28

def tentori_belloni_reactivity(
    ion_temperature_profile: np.ndarray,
) -> np.ndarray:
    """Calculate the volumetric fusion reaction rate 〈sigmav〉 (m³/s) for 
    p-B reaction using the Tentori-Belloni parametrization.
    """
    mp = 1.007276466621 * 1.66053906660e-27
    mb = 11 * mp
    mr = (mp * mb) / (mp + mb) 

    e = 10 ** np.arange(-2, np.log10(9.76), 0.0005)  # MeV

    sigma = np.array([fsgmpb_tentori_belloni(ee) for ee in e])

    # analytic integral
    sgmv = np.zeros_like(ion_temperature_profile)
    if isinstance(ion_temperature_profile, float):
        len_teff = 1
    else:
        len_teff = len(ion_temperature_profile)
    for j in range(len_teff):
        if len_teff == 1:
            teff_now = ion_temperature_profile
        else:
            teff_now = ion_temperature_profile[j]
        sgmv_now = (
            np.sqrt(constants.KILOELECTRON_VOLT * teff_now * 8 / (np.pi * mr))
            / (teff_now ** 2)
            * np.nansum(e[1:] * sigma[1:] * np.diff(e) * np.exp(-e[1:] * 1.0e3 / teff_now))
        ) * 1.0e6
        if len_teff == 1:
            sgmv = sgmv_now
        else:
            sgmv[j] = sgmv_now
    return sgmv

def set_fusion_powers(
    f_alpha_electron: float,
    f_alpha_ion: float,
    p_beam_alpha_mw: float,
    p_beam_neutron_mw: float,
    pden_non_alpha_charged_mw: float,
    pden_plasma_neutron_vol_avg_mw: float,
    vol_plasma: float,
    pden_plasma_alpha_vol_avg_mw: float,
    f_p_alpha_plasma_deposited: float,
) -> tuple:
    """Computes various fusion power metrics based on the provided
    plasma parameters.

    Parameters
    ----------
    f_alpha_electron :
        float
    f_alpha_ion :
        float
    p_beam_alpha_mw :
        float
    pden_non_alpha_charged_mw :
        float
    pden_plasma_neutron_vol_avg_mw :
        float
    vol_plasma :
        float
    pden_plasma_alpha_vol_avg_mw :
        float
    f_p_alpha_plasma_deposited: float
        Fraction of alpha power deposited in plasma

    Returns
    -------
    :
        tuple: A tuple containing the following elements:
        - pden_neutron_total_vol_avg_mw (float): Neutron fusion power per unit volume
          fromplasma and beams [MW/m³].
        - p_plasma_alpha_mw (float): Alpha fusion power from only the plasma [MW].
        - p_alpha_total_mw (float): Total alpha fusion power from plasma and beams [MW].
        - p_plasma_neutron_mw (float): Neutron fusion power from only the plasma [MW].
        - p_neutron_total_mw (float): Total neutron fusion power from plasma and
          beams [MW].
        - p_non_alpha_charged_mw (float): Other total charged particle fusion power [MW].
        - pden_alpha_total_vol_avg_mw (float): Alpha power per unit volume, from beams
          and plasma [MW/m³].
        - f_pden_alpha_electron_mw (float): Alpha power per unit volume to
          electrons [MW/m³].
        - f_pden_alpha_ions_mw (float): Alpha power per unit volume to ions [MW/m³].
        - p_charged_particle_mw (float): Charged particle fusion power [MW].
        - p_fusion_total_mw (float): Total fusion power [MW].

    References
    ----------
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    """
    # Alpha power

    # Calculate alpha power produced just by the plasma
    p_plasma_alpha_mw = pden_plasma_alpha_vol_avg_mw * vol_plasma

    # Add neutral beam alpha power / volume
    pden_alpha_total_vol_avg_mw = pden_plasma_alpha_vol_avg_mw + (
        p_beam_alpha_mw / vol_plasma
    )

    # Total alpha power
    p_alpha_total_mw = pden_alpha_total_vol_avg_mw * vol_plasma

    # Neutron Power

    # Calculate neutron power produced just by the plasma
    p_plasma_neutron_mw = pden_plasma_neutron_vol_avg_mw * vol_plasma

    # Add extra neutron power from beams
    pden_neutron_total_vol_avg_mw = pden_plasma_neutron_vol_avg_mw + (
        p_beam_neutron_mw
        / vol_plasma
    )

    # Total neutron power
    p_neutron_total_mw = pden_neutron_total_vol_avg_mw * vol_plasma

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
        f_p_alpha_plasma_deposited * pden_alpha_total_vol_avg_mw * f_alpha_ion
    )
    f_pden_alpha_electron_mw = (
        f_p_alpha_plasma_deposited * pden_alpha_total_vol_avg_mw * f_alpha_electron
    )

    return (
        pden_neutron_total_vol_avg_mw,
        p_plasma_alpha_mw,
        p_alpha_total_mw,
        p_plasma_neutron_mw,
        p_neutron_total_mw,
        p_non_alpha_charged_mw,
        pden_alpha_total_vol_avg_mw,
        f_pden_alpha_electron_mw,
        f_pden_alpha_ions_mw,
        p_charged_particle_mw,
        p_fusion_total_mw,
    )


def beam_fusion(
    beamfus0: float,
    betbm0: float,
    b_plasma_total: float,
    c_beam_total: float,
    nd_plasma_electrons_vol_avg: float,
    nd_plasma_fuel_ions_vol_avg: float,
    ion_electron_coulomb_log: float,
    e_beam_kev: float,
    f_deuterium_plasma: float,
    f_tritium_plasma: float,
    f_beam_tritium: float,
    temp_plasma_electron_density_weighted_kev: float,
    vol_plasma: float,
    n_charge_plasma_effective_mass_weighted_vol_avg: float,
) -> tuple[float, float, float]:
    """Calculate neutral beam slowing-down properties and beam-target fusion.

    This function computes the neutral beam beta contribution, hot beam ion density,
    and alpha power generated from beam-target fusion reactions in the plasma.

    Parameters
    ----------
    beamfus0 :
        Multiplier for beam target fusion calculation.
    betbm0 :
        Leading coefficient for neutral beam beta fraction.
    b_plasma_total :
        Total magnetic field (T).
    c_beam_total :
        Total neutral beam current (A).
    nd_plasma_electrons_vol_avg :
        Volume-averaged electron density (m^-3).
    nd_plasma_fuel_ions_vol_avg :
        Volume-averaged fuel ion density (m^-3).
    ion_electron_coulomb_log :
        Ion-electron Coulomb logarithm.
    e_beam_kev :
        Neutral beam injection energy (keV).
    f_deuterium_plasma :
        Deuterium fraction in the bulk plasma.
    f_tritium_plasma :
        Tritium fraction in the bulk plasma.
    f_beam_tritium :
        Tritium fraction in the injected beam.
    temp_plasma_electron_density_weighted_kev :
        Density weighted average electron temperature (keV).
    vol_plasma :
        Plasma volume (m^3).
    n_charge_plasma_effective_mass_weighted_vol_avg :
        Mass-weighted effective plasma charge.

    Returns
    -------
    beta_beam :
        Neutral beam beta contribution.
    n_beam :
        Hot beam ion density (m^-3).
    p_alpha_beam :
        Alpha power from beam-target fusion (MW).

    Notes
    -----
        - The function uses the Bosch-Hale parametrization to evaluate fusion reactivity.
        - Beam ion slowing-down properties are computed using Coulomb collision physics.
        - Fusion power is obtained by integrating the beam-target reaction rate over
        the fast-ion distribution.

    References
    ----------
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and
        thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611-631, Apr. 1992,
        doi: https://doi.org/10.1088/0029-5515/32/4/i07.

        - J. W. Sheffield, “The physics of magnetic fusion reactors,”
        Rev. Mod. Phys., vol. 66, no. 3, pp. 1015-1103, Jul. 1994, Jul. 1994,
        doi: https://doi.org/10.1103/revmodphys.66.1015.

        - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,”
        Nuclear Fusion and Plasma Physics, vol. 9, no. 3, pp. 136-141, 2022.
        Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf
    """
    # Beam ion slowing down time given by Deng Baiquan and G. A. Emmert 1987
    t_beam_slow = (
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
    nd_plasma_deuterium = nd_plasma_fuel_ions_vol_avg * f_deuterium_plasma
    nd_plasma_tritium = nd_plasma_fuel_ions_vol_avg * f_tritium_plasma

    beam_state = beam_slowing_down_state(
        e_beam_kev,
        critical_energy_deuterium,
        critical_energy_tritium,
        t_beam_slow,
        f_beam_tritium,
        c_beam_total,
        vol_plasma,
    )

    sigv_beam_deuterium = beam_reaction_rate_coefficient(
        constants.M_DEUTERON_AMU,
        beam_state.deuterium_critical_energy_speed,
        e_beam_kev,
    )

    sigv_beam_tritium = beam_reaction_rate_coefficient(
        constants.M_TRITON_AMU,
        beam_state.tritium_critical_energy_speed,
        e_beam_kev,
    )

    reaction_rate_beam_deuterium_dt = beam_target_reaction_rate(
        beam_state.deuterium_beam_density,
        nd_plasma_tritium,
        sigv_beam_deuterium,
        vol_plasma,
    )
    p_alpha_beam_deuterium_dt_mw = alpha_power_beam(
        reaction_rate_beam_deuterium_dt,
    )

    reaction_rate_beam_tritium_dt = beam_target_reaction_rate(
        beam_state.tritium_beam_density,
        nd_plasma_deuterium,
        sigv_beam_tritium,
        vol_plasma,
    )
    p_alpha_beam_tritium_dt_mw = alpha_power_beam(
        reaction_rate_beam_tritium_dt,
    )

    # Neutral beam alpha power
    p_beam_alpha_mw = beamfus0 * (
        p_alpha_beam_deuterium_dt_mw + p_alpha_beam_tritium_dt_mw
    )

    # Neutral beam beta
    beta_beam = (
        betbm0
        * 4.03e-22
        * (2 / 3)
        * beam_state.nd_beam_hot
        * beam_state.e_beam_deposited_kev
        / (b_plasma_total**2)
    )

    return beta_beam, beam_state.nd_beam_hot, p_beam_alpha_mw

def sigv_beam_target_integrand(e_rel, v_th, vb, m_fi, m):
    mu = m_fi * m / (m_fi + m) * constants.ATOMIC_MASS_UNIT # reduced mass (kg)        
    vr = np.sqrt(2.0 * e_rel * constants.KILOELECTRON_VOLT / mu)
    fexp = np.exp(-((vb - vr) / v_th)**2) - np.exp(-((vb + vr) / v_th)**2)
    sig = fsgmpb(e_rel)
    return sig * fexp * vr * constants.KILOELECTRON_VOLT / mu

def sigv_beam_target(eb, v_th, g, m_fi, m):
    vb = np.sqrt(2.0 * eb * constants.KILOELECTRON_VOLT / (m_fi * constants.ATOMIC_MASS_UNIT))
    # 被积函数在 vr≈vb 处呈热速度宽度的尖峰；内层积分值常 ~1e-9，
    # 默认 epsabs(~1e-8) 过大，会触发 IntegrationWarning(roundoff)。
    # 用相对容差，并把积分限收在峰附近若干 v_th。
    mu = m_fi * m / (m_fi + m) * constants.ATOMIC_MASS_UNIT
    vr_lo = max(vb - 8.0 * v_th, 0.0)
    vr_hi = vb + 8.0 * v_th
    e_lo = max(0.5 * mu * vr_lo**2 / constants.KILOELECTRON_VOLT, 1.0e-5)
    e_hi = max(0.5 * mu * vr_hi**2 / constants.KILOELECTRON_VOLT, e_lo * 10.0)
    sigv = integrate.quad(
        sigv_beam_target_integrand,
        e_lo,
        e_hi,
        args=(v_th, vb, m_fi, m),
        epsabs=0.0,
        epsrel=1.0e-8,
        limit=200,
    )[0]
    sigv += integrate.quad(sigv_beam_target_integrand, 1.0e-5, e_lo, args=(v_th, vb, m_fi, m))[0]
    if e_hi < 1000.0:
        sigv += integrate.quad(sigv_beam_target_integrand, e_hi, 1000.0, args=(v_th, vb, m_fi, m))[0]
    coef = 1.0 / (v_th * vb * np.sqrt(np.pi))
    return sigv * coef * g(eb)

def cal_beam_fusion_sigv(e_beam_kev, tau_slow_beam, e_crit_beam, s0_beam, m_fi, n_fi, m, temp):
    """
    积分外层函数得到反应率 <sigma*v>
    """

    dedt = lambda e, t, e_crit: -2.0 * e / t * (1.0 + (e_crit / e)**1.5)
    g_beam = lambda e, t, e_crit, s0: - s0 / dedt(e, t, e_crit) * (e_crit > e) # m^-3 keV^-1
    g = lambda e: g_beam(e, tau_slow_beam, e_crit_beam, s0_beam)

    v_th = np.sqrt(2.0 * temp * constants.KILOELECTRON_VOLT / (m * constants.ATOMIC_MASS_UNIT))
    # 嵌套积分：外层同样避免过紧的绝对容差
    sigv = integrate.quad(
        sigv_beam_target,
        1.0e-5,
        e_beam_kev,
        args=(v_th, g, m_fi, m),
        epsabs=0.0,
        epsrel=1.0e-6,
        limit=200,
    )[0] / n_fi
    return sigv

def beam_fusion_p_b11(
    sigmv_pb: float,
    i_pb_reactivity_model: int,
    ti: float,
    te: float,
    vol_plasma: float,
    current_beam_total: float,
    b_total: float,
    ne: float,
    nd_fuel_b11: float,
    charge_eff_mass_weighted: float,
    dlamie: float,
    dlamee: float,
    betbm0: float,
    beamfus0: float,
    e_beam_kev: float,
):

    if i_pb_reactivity_model == 1:
        f_avg = sigmv_pb / nevins_reactivity(ti)
    elif i_pb_reactivity_model == 2:
        f_avg = sigmv_pb / tentori_belloni_reactivity(ti)
    else:
        raise ValueError(f"Invalid PB reaction model: {i_pb_reactivity_model}")
    y = 8.68 * constants.ELECTRON_VOLT # MJ

    #
    const3 = 1000.0 * constants.ELECTRON_VOLT / constants.ATOMIC_MASS_UNIT
    inv_e_vol = 1.0 / (vol_plasma * constants.ELECTRON_VOLT)
    me = constants.ELECTRON_MASS


    # reference:
    # 1. the particle kinetics of the plasmas, page 307, subsection 7-11
    # 2. fast ion pressure in fusion plasmas

    # critical velocity
    # assume：
    # 0.proton beam, no boron
    # 1.background particle is of Maxwellian distribution
    # 2.x_i = v_fast_ion / v_ion > 1, x_e = v_fast_ion / v_electron < 1 => use asymptotic error function
    # Gi = 1, Ge = 4/(3 * sqrt(pi)) * xe**3
    # 3.ideal gas, E = 1/2 * m * v**2
    v_crit = (
        (
            3.0 * np.sqrt(np.pi) /4.0 * me / constants.ATOMIC_MASS_UNIT 
            / ne * (charge_eff_mass_weighted * dlamie) / dlamee
        )**(1.0/3.0) 
        * np.sqrt(2 * te / me * constants.KILOELECTRON_VOLT)
    ) # m/s
    e_crit_p = 0.5 * constants.M_PROTON_AMU * v_crit**2 / const3 # keV


    # slowing down time
    # × (dlamie / dlamee) not included in the formula in the reference
    tau_slow_p = (
        3.0 * (np.sqrt(2.0*np.pi))**3 * constants.EPSILON0**2 
        * (te * constants.KILOELECTRON_VOLT)**1.5 * constants.M_PROTON_AMU        
        / (np.sqrt(me) * constants.ELECTRON_CHARGE**4 * ne * dlamee) * (dlamie / dlamee)
    ) * constants.ATOMIC_MASS_UNIT

    # s0: unit time, unit volume of fast ion production
    s0_p = current_beam_total * inv_e_vol


    # fast ion density(fast ion density is the integral of the distribution function g from 0 to e_beam_kev)
    x_crit_p = np.sqrt(e_beam_kev / e_crit_p)
    nd_beam_p = tau_slow_p * s0_p * np.log(1.0 + x_crit_p**3) / 3.0
    nd_beam = nd_beam_p


    # fast ion pressure
    # injection speed
    v_p_i = np.sqrt(2 * e_beam_kev / constants.M_PROTON_AMU * const3)
    pressure = lambda x, m, t, s0: constants.ATOMIC_MASS_UNIT * m * s0 * t * v_crit**2 / 3.0 *(
        -np.atan(1/np.sqrt(3.0))/np.sqrt(3.0) + 0.5 * x**2 
        - 1/np.sqrt(3.0) * np.atan((2.0 * x - 1) / np.sqrt(3.0))
        + 0.5 * np.log(1 + x) - 1/6.0 * np.log(1 + x**3)
    ) # Pa
    pressure_p = pressure(v_p_i / v_crit, constants.M_PROTON_AMU, tau_slow_p, s0_p)

    sigv_p = cal_beam_fusion_sigv(
        e_beam_kev, tau_slow_p, e_crit_p, s0_p, constants.M_PROTON_AMU, nd_beam_p, constants.M_BORON11_AMU, ti
    )

    sigmv_p = sigv_p * f_avg

    reaction_rate_beam = sigmv_p * beamfus0 * nd_fuel_b11 * nd_beam_p

    p_beam_alpha_mw = reaction_rate_beam * y * vol_plasma * beamfus0
    
    beta_beam = betbm0 * 2 * constants.RMU0 * pressure_p / b_total**2

    return beta_beam, nd_beam, p_beam_alpha_mw

@dataclass(frozen=True)
class BeamSlowingDownState:
    """Beam slowing-down state and hot-ion properties.

    Attributes
    ----------
    deuterium_beam_density :
        Steady-state deuterium beam ion density (m^-3).
    tritium_beam_density :
        Steady-state tritium beam ion density (m^-3).
    deuterium_critical_energy_speed :
        Critical speed of deuterium beam ions corresponding to the
        critical energy (m/s).
    tritium_critical_energy_speed :
        Critical speed of tritium beam ions corresponding to the
        critical energy (m/s).
    nd_beam_hot :
        Total hot beam ion density (m^-3).
    e_beam_deposited_kev :
        Density-weighted average deposited beam ion energy (keV).
    """

    deuterium_beam_density: float
    tritium_beam_density: float
    deuterium_critical_energy_speed: float
    tritium_critical_energy_speed: float
    nd_beam_hot: float
    e_beam_deposited_kev: float


def beam_slowing_down_state(
    e_beam_kev: float,
    critical_energy_deuterium: float,
    critical_energy_tritium: float,
    t_beam_slow: float,
    f_beam_tritium: float,
    c_beam_total: float,
    vol_plasma: float,
) -> BeamSlowingDownState:
    """Calculate beam slowing-down state and hot ion properties.

    Estimates the steady-state hot D and T beam-ion densities from beam
    source rate times slowing-down residence time, computes species critical
    velocities from critical energies, evaluates a slowing-down-distribution
    fast-ion pressure for each species, and converts that pressure into a
    pressure-equivalent mean fast-ion energy.

    Parameters
    ----------
    e_beam_kev :
        Beam energy (keV).
    critical_energy_deuterium :
        Critical energy for electron/ion slowing down of deuterium beam ions (keV).
    critical_energy_tritium :
        Critical energy for electron/ion slowing down of tritium beam ions (keV).
    t_beam_slow :
        Beam ion slowing-down time on electrons (s).
    f_beam_tritium :
        Tritium fraction in the injected beam (0.0 = pure deuterium beam).
    c_beam_total :
        Total neutral beam current (A).
    vol_plasma :
        Plasma volume (m^3).

    Returns
    -------
    :
        Beam slowing-down state containing deuterium and tritium hot beam-ion
        densities, deuterium and tritium critical speeds, total hot beam-ion
        density, and average deposited beam-ion energy.

    Notes
    -----
        - The critical energy for electron/ion slowing down of the beam ion is calculated
        for both deuterium and tritium neutral beams.
        - The function integrates the hot beam fusion reaction rate integrand over the
        range of beam velocities up to the critical velocity.

    References
    ----------
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections
        and thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611-631,
        Apr. 1992, doi: https://doi.org/10.1088/0029-5515/32/4/i07.

        - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,”
        Nuclear Fusion and Plasma Physics, vol. 9, no. 3, pp. 136-141, 2022,
        Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

        - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
        International Series of Monographs on Physics, Volume 149.

        - J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3,
        pp. 1015-1103, Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.
    """
    # D and T beam current fractions
    beam_current_deuterium = c_beam_total * (1.0 - f_beam_tritium)
    beam_current_tritium = c_beam_total * f_beam_tritium

    # At a critical energy the rate of loss to the ions becomes equal to that to
    # the electrons, and at lower energies the loss to the ions predominates.

    # Ratio of beam energy to critical energy for deuterium
    beam_energy_ratio_deuterium = e_beam_kev / critical_energy_deuterium

    # Calculate the characterstic time for the deuterium ions to slow down to the
    # thermal energy, eg E = 0.
    characteristic_deuterium_t_beam_slow = (
        t_beam_slow / 3.0 * np.log(1.0 + (beam_energy_ratio_deuterium) ** 1.5)
    )

    deuterium_beam_density = (
        beam_current_deuterium
        * characteristic_deuterium_t_beam_slow
        / (constants.ELECTRON_CHARGE * vol_plasma)
    )

    # Ratio of beam energy to critical energy for tritium
    beam_energy_ratio_tritium = e_beam_kev / critical_energy_tritium

    # Calculate the characterstic time for the tritium to slow down to the thermal
    # energy, eg E = 0.
    # Wesson, J. (2011) Tokamaks.
    characteristic_tritium_t_beam_slow = (
        t_beam_slow / 3.0 * np.log(1.0 + (beam_energy_ratio_tritium) ** 1.5)
    )

    tritium_beam_density = (
        beam_current_tritium
        * characteristic_tritium_t_beam_slow
        / (constants.ELECTRON_CHARGE * vol_plasma)
    )

    nd_beam_hot = deuterium_beam_density + tritium_beam_density

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
    # D.Baiquan et.al.  “Fast ion pressure in fusion plasma,” Nuclear Fusion and
    # Plasma Physics, vol. 9, no. 3, pp. 136-141, 2022,
    # Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

    source_deuterium = beam_current_deuterium / (constants.ELECTRON_CHARGE * vol_plasma)

    source_tritium = beam_current_tritium / (constants.ELECTRON_CHARGE * vol_plasma)

    pressure_coeff_deuterium = (
        constants.M_DEUTERON_AMU
        * constants.ATOMIC_MASS_UNIT
        * t_beam_slow
        * deuterium_critical_energy_speed**2
        * source_deuterium
        / (constants.KILOELECTRON_VOLT * 3.0)
    )
    pressure_coeff_tritium = (
        constants.M_TRITON_AMU
        * constants.ATOMIC_MASS_UNIT
        * t_beam_slow
        * tritium_critical_energy_speed**2
        * source_tritium
        / (constants.KILOELECTRON_VOLT * 3.0)
    )

    # Fast Ion Pressure
    # This is the same form as the ideal gas law pressure, P=1/3 * nmv²
    deuterium_pressure = pressure_coeff_deuterium * fast_ion_pressure_integral(
        e_beam_kev, critical_energy_deuterium
    )
    tritium_pressure = pressure_coeff_tritium * fast_ion_pressure_integral(
        e_beam_kev, critical_energy_tritium
    )

    # Beam deposited energy
    # Find the energy from the ideal gas pressure, P=1/3 * nmv² = 2/3 * n<E>
    deuterium_deposited_energy = (
        1.5 * deuterium_pressure / deuterium_beam_density
        if deuterium_beam_density > 0.0
        else 0.0
    )
    tritium_deposited_energy = (
        1.5 * tritium_pressure / tritium_beam_density
        if tritium_beam_density > 0.0
        else 0.0
    )

    e_beam_deposited_kev = (
        (
            deuterium_beam_density * deuterium_deposited_energy
            + tritium_beam_density * tritium_deposited_energy
        )
        / nd_beam_hot
        if nd_beam_hot > 0.0
        else 0.0
    )

    return BeamSlowingDownState(
        deuterium_beam_density=deuterium_beam_density,
        tritium_beam_density=tritium_beam_density,
        deuterium_critical_energy_speed=deuterium_critical_energy_speed,
        tritium_critical_energy_speed=tritium_critical_energy_speed,
        nd_beam_hot=nd_beam_hot,
        e_beam_deposited_kev=e_beam_deposited_kev,
    )


def fast_ion_pressure_integral(e_beam_kev: float, critical_energy: float) -> float:
    """Evaluate the dimensionless fast-ion pressure integral.

    This function returns the closed-form dimensionless integral factor
    appearing in the analytic fast-ion pressure expression for a
    monoenergetic source with a classical slowing-down distribution.

    Parameters
    ----------
    e_beam_kev :
        Beam birth energy (keV).
    critical_energy :
        Critical energy for electron/ion slowing down (keV).

    Returns
    -------
    :
        Dimensionless pressure integral factor.

    Notes
    -----
        - The function uses the ratio of the beam energy to the critical energy to
          compute the hot ion energy parameter.
        - The calculation involves logarithmic and arctangent functions to account for
          the energy distribution of the hot ions.

    References
    ----------
        - Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,”
        Nuclear Fusion and Plasma Physics,
        vol. 9, no. 3, pp. 136-141, 2022,
        Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf

        - W.A Houlberg, “Thermalization of an Energetic Heavy Ion in a Multi-species
        Plasma,” University of Wisconsin Fusion Technology Institute,
        Report UWFDM-103 1974,
        Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm103.pdf
    """
    xcs = e_beam_kev / critical_energy
    xc = np.sqrt(xcs)

    t1 = xcs / 2.0
    t2 = np.log((xcs + 2.0 * xc + 1.0) / (xcs - xc + 1.0)) / 6.0

    xarg = (2.0 * xc - 1.0) / np.sqrt(3.0)
    t3 = np.arctan(xarg) / np.sqrt(3.0)
    t4 = (1 / np.sqrt(3.0)) * np.arctan(1 / np.sqrt(3.0))

    return t1 + t2 - t3 - t4


def beam_target_reaction_rate(
    nd_beam_ion: float,
    nd_target_ion: float,
    sigv_beam: float,
    vol_plasma: float,
) -> float:
    """Calculate the total beam-target fusion reaction rate.

    This function computes the total number of beam-target fusion reactions
    occurring per second in the plasma.

    The total beam-target reaction rate is given by:
        reaction_rate = nd_beam_ion * nd_target_ion * sigv_beam * vol_plasma

    Unit analysis:
        (m^-3) * (m^-3) * (m^3/s) * (m^3) = s^-1

    Parameters
    ----------
    nd_beam_ion :
        Hot beam ion density (m^-3).
    nd_target_ion :
        Thermal target ion density (m^-3).
    sigv_beam :
        Beam-target fusion reactivity (m^3/s).
    vol_plasma :
        Plasma volume (m^3).

    Returns
    -------
    :
        Total beam-target fusion reaction rate (s^-1).

    Notes
    -----
    Any spatial/profile effects should be introduced through the
    supplied nd_target_ion.
    """
    return nd_beam_ion * nd_target_ion * sigv_beam * vol_plasma


def alpha_power_beam(
    beam_target_reaction_rate: float,
) -> float:
    """Calculate alpha power from beam-target fusion.

    This function converts the total beam-target fusion reaction rate into
    alpha particle power.

    The alpha power is given by:
        P_alpha = beam_target_reaction_rate * DT_ALPHA_ENERGY

    Unit analysis:
        (s^-1) * (J) = W

    The result is converted from watts to megawatts (MW).

    Parameters
    ----------
    beam_target_reaction_rate :
        Total beam-target fusion reaction rate (s^-1).

    Returns
    -------
    :
        Alpha power from beam-target fusion (MW).
    """
    return beam_target_reaction_rate * constants.DT_ALPHA_ENERGY / 1.0e6


def beam_reaction_rate_coefficient(
    relative_mass_ion: float, critical_velocity: float, beam_energy_kev: float
) -> float:
    """Calculate an effective beam-target fusion rate coefficient.

    This function evaluates an effective fusion rate coefficient for a
    slowing-down fast-ion population interacting with a thermal target
    plasma. The beam ions are assumed to follow a classical slowing-down
    distribution characterised by the critical velocity, and the fusion
    contribution is obtained by integrating a beam-energy-dependent fusion
    cross-section over the slowing-down distribution.

    Parameters
    ----------
    relative_mass_ion :
        Relative atomic mass of the ion (e.g. approx 2.0 for D, 3.0 for T).
    critical_velocity :
        Critical velocity for electron/ion slowing down of the beam ion [m/s].
    beam_energy_kev :
        Neutral beam energy [keV].

    Returns
    -------
    :
        Hot beam fusion rate coefficient [m^3/s].

    Notes
    -----
    - The function returns a reactivity rate coefficient, not the reaction rate itself.
    - The integration variable is the beam speed normalised to the critical speed.
    - The underlying beam fusion cross-section is evaluated in cm^2 and converted
      internally to m^2 so that this function returns an SI rate coefficient.
    - The integration is performed using scipy.integrate.quad.

    References
    ----------
    - S. Niikura and M. Nagami, "Improvement of fusion reactivity and fusion
      power multiplication factor in the presence of fast ions," Fusion
      Engineering and Design, vol. 12, pp. 467-480, 1990.
    """
    # Find the speed of the beam particle when it has the critical energy.
    # Re-arrange kinetic energy equation to find speed. Non-relativistic.
    beam_velocity = np.sqrt(
        (beam_energy_kev * constants.KILOELECTRON_VOLT)
        * 2.0
        / (relative_mass_ion * constants.ATOMIC_MASS_UNIT)
    )

    relative_velocity = beam_velocity / critical_velocity
    integral_coefficient = 3.0 * critical_velocity / np.log(1.0 + relative_velocity**3)

    fusion_integral_cm2 = integrate.quad(
        _hot_beam_fusion_reaction_rate_integrand,
        0.0,
        relative_velocity,
        args=(critical_velocity,),
    )[0]

    # _hot_beam_fusion_cross_section returns cm^2, so the integrated quantity is in
    # cm^2. Convert to m^2 here so the final rate coefficient is returned in m^3/s.
    fusion_integral_m2 = fusion_integral_cm2 * 1.0e-4

    return integral_coefficient * fusion_integral_m2


def _hot_beam_fusion_reaction_rate_integrand(
    velocity_ratio: float, critical_velocity: float
) -> float:
    """Evaluate the integrand for the beam-target fusion rate coefficient.

    This function returns the slowing-down-weighted integrand used to
    compute the effective beam-target fusion rate coefficient for a fast-ion
    population. The integration variable is the beam speed normalised to the
    critical speed, and the integrand combines a slowing-down weighting with
    the beam fusion cross-section evaluated at the corresponding beam energy.

    Parameters
    ----------
    velocity_ratio :
        Ratio of beam velocity to the critical velocity.
    critical_velocity :
        Critical velocity for electron/ion slowing down of the beam ion (m/s).


    Returns
    -------
    :
        Value of the integrand for the hot beam fusion reaction rate.

    Notes
    -----
    - The factor `velocity_ratio**3 / (1.0 + velocity_ratio**3)` represents
      the slowing-down weighting in the normalised velocity variable.
    - The beam energy passed to `_beam_fusion_cross_section` is constructed
      from the instantaneous beam speed and expressed in keV/amu.
    - This function is an internal helper used by
      `beam_reaction_rate_coefficient`.

    References
    ----------
    - S. Niikura and M. Nagami, "Improvement of fusion reactivity and fusion
      power multiplication factor in the presence of fast ions," Fusion
      Engineering and Design, vol. 12, pp. 467-480, 1990.

    """
    intgeral_term = (velocity_ratio**3) / (1.0 + velocity_ratio**3)

    # critical_velocity : critical velocity for electron/ion slowing down of beam
    # ion (m/s)
    beam_velcoity = critical_velocity * velocity_ratio

    # Calculate the beam kinetic energy per amu and normalise to keV
    xvcs = beam_velcoity**2 * constants.ATOMIC_MASS_UNIT / (constants.KILOELECTRON_VOLT)

    # Calculate the fusion reaction cross-section from beam kinetic energy
    cross_section = _beam_fusion_cross_section(xvcs)

    return intgeral_term * cross_section


def _beam_fusion_cross_section(vrelsq: float) -> float:
    """Calculate the beam fusion cross-section from beam energy.

    This function computes the fusion reaction cross-section based on the
    square of the speed of the beam ion (keV/amu). The functional form of
    the cross-section is in terms of the equivalent deuterium energy, i.e.,
    for a tritium beam at 500 keV the energy used in the cross-section
    function is 333 keV.

    Parameters
    ----------
    vrelsq :
        Square of the speed of the beam ion (keV/amu).

    Returns
    -------
    :
        Fusion reaction cross-section (cm²).

    Notes
    -----
    - The cross-section is limited at low and high beam energies.
    - For beam kinetic energy less than 10 keV, the cross-section is set to 1.0e-27 cm².
    - For beam kinetic energy greater than 10,000 keV, the cross-section is set to
      8.0e-26 cm².
    - The cross-section is calculated using a functional form with parameters a1 to a5.
    - The exact provenance of this particular fit and its coefficients has
      not yet been traced in the present code documentation.
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
