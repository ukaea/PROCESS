import logging
import numpy
from scipy import integrate
from dataclasses import dataclass
from process.fortran import physics_variables, physics_module, constants
from process.plasma_profiles import PlasmaProfile
import process.impurity_radiation as impurity


logger = logging.getLogger(__name__)

ATOMIC_MASS_DEUTERIUM = 2.0
ATOMIC_MASS_TRITIUM = 3.0

REACTION_CONSTANTS_DT = dict(
    bg=34.3827,
    mrc2=1.124656e6,
    cc1=1.17302e-9,
    cc2=1.51361e-2,
    cc3=7.51886e-2,
    cc4=4.60643e-3,
    cc5=1.35000e-2,
    cc6=-1.06750e-4,
    cc7=1.36600e-5,
)

REACTION_CONSTANTS_DHE3 = dict(
    bg=68.7508,
    mrc2=1.124572e6,
    cc1=5.51036e-10,
    cc2=6.41918e-3,
    cc3=-2.02896e-3,
    cc4=-1.91080e-5,
    cc5=1.35776e-4,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD1 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.43360e-12,
    cc2=5.85778e-3,
    cc3=7.68222e-3,
    cc4=0.0,
    cc5=-2.96400e-6,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD2 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.65718e-12,
    cc2=3.41267e-3,
    cc3=1.99167e-3,
    cc4=0.0,
    cc5=1.05060e-5,
    cc6=0.0,
    cc7=0.0,
)

ATOMIC_MASS_DEUTERIUM = 2.0
ATOMIC_MASS_TRITIUM = 3.0

REACTION_CONSTANTS_DT = dict(
    bg=34.3827,
    mrc2=1.124656e6,
    cc1=1.17302e-9,
    cc2=1.51361e-2,
    cc3=7.51886e-2,
    cc4=4.60643e-3,
    cc5=1.35000e-2,
    cc6=-1.06750e-4,
    cc7=1.36600e-5,
)

REACTION_CONSTANTS_DHE3 = dict(
    bg=68.7508,
    mrc2=1.124572e6,
    cc1=5.51036e-10,
    cc2=6.41918e-3,
    cc3=-2.02896e-3,
    cc4=-1.91080e-5,
    cc5=1.35776e-4,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD1 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.43360e-12,
    cc2=5.85778e-3,
    cc3=7.68222e-3,
    cc4=0.0,
    cc5=-2.96400e-6,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD2 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.65718e-12,
    cc2=3.41267e-3,
    cc3=1.99167e-3,
    cc4=0.0,
    cc5=1.05060e-5,
    cc6=0.0,
    cc7=0.0,
)


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
        charged_power_density (float): Power density of charged particles produced.
        neutron_power_density (float): Power density of neutrons produced.
        fusion_rate_density (float): Fusion reaction rate density.
        alpha_rate_density (float): Alpha particle production rate density.
        proton_rate_density (float): Proton production rate density.

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
          Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,
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
            charged_power_density (float): Power density of charged particles produced.
            neutron_power_density (float): Power density of neutrons produced.
            fusion_rate_density (float): Fusion reaction rate density.
            alpha_rate_density (float): Alpha particle production rate density.
            proton_rate_density (float): Proton production rate density.
        """
        self.plasma_profile = plasma_profile
        self.sigmav_dt_average = 0.0
        self.dhe3_power_density = 0.0
        self.dd_power_density = 0.0
        self.dt_power_density = 0.0
        self.alpha_power_density = 0.0
        self.charged_power_density = 0.0
        self.neutron_power_density = 0.0
        self.fusion_rate_density = 0.0
        self.alpha_rate_density = 0.0
        self.proton_rate_density = 0.0

    def deuterium_branching(self, ion_temperature: float) -> float:
        """
        Calculate the relative rate of tritium producing D-D reactions to 3He ones based on the volume averaged ion temperature

        Parameters:
            ion_temperature (float): Volume averaged ion temperature in keV

         Returns:
            float: The rate of tritium producing D-D reactions to 3He ones

        Notes:
            - For ion temperatures between 0.5 keV and 200 keV.
            - The deviation of the fit from the R-matrix branching ratio is always smaller than 0.5%.

        References:
            - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
              Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,
              doi: https://doi.org/10.1088/0029-5515/32/4/i07.
        ‌
        """
        return (
            1.02934
            - 8.3264e-3 * ion_temperature
            + 1.7631e-4 * ion_temperature**2
            - 1.8201e-6 * ion_temperature**3
            + 6.9855e-9 * ion_temperature**4
        )

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
            - self.charged_power_density: Power density of charged particles produced.
            - self.neutron_power_density: Power density of neutrons produced.
            - self.fusion_rate_density: Fusion reaction rate density.
            - self.alpha_rate_density: Alpha particle production rate density.
            - self.proton_rate_density: Proton production rate density.

        Returns:
            None
        """
        # Initialize Bosch-Hale constants for the D-T reaction
        dt = BoschHaleConstants(**REACTION_CONSTANTS_DT)

        # Calculate the fusion reaction rate integral using Simpson's rule
        sigmav = integrate.simpson(
            fusion_rate_integral(self.plasma_profile, dt),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )

        # Store the average fusion reaction rate
        self.sigmav_dt_average = sigmav

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.d_t_energy / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * (physics_variables.fdeut * physics_variables.deni)
            * (physics_variables.ftrit * physics_variables.deni)
        )

        # Power densities for different particles [MW/m^3]
        # Alpha particle gets approximately 20% of the fusion power
        alpha_power_density = (
            1.0 - constants.dt_neutron_energy_fraction
        ) * fusion_power_density
        charged_power_density = 0.0
        neutron_power_density = (
            constants.dt_neutron_energy_fraction * fusion_power_density
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
            charged_power_density,
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
            - self.charged_power_density: Power density of charged particles produced.
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

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.d_helium_energy / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * (physics_variables.fdeut * physics_variables.deni)
            * (physics_variables.fhe3 * physics_variables.deni)
        )

        # Power densities for different particles [MW/m^3]
        # Alpha particle gets approximately 20% of the fusion power
        alpha_power_density = (
            1.0 - constants.dhelium_proton_energy_fraction
        ) * fusion_power_density
        charged_power_density = (
            constants.dhelium_proton_energy_fraction * fusion_power_density
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
            charged_power_density,
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
            - self.charged_power_density: Power density of charged particles produced.
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

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.dd_helium_energy / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * 0.5  # Factor for D-D reaction
            * (physics_variables.fdeut * physics_variables.deni)
            * (physics_variables.fdeut * physics_variables.deni)
        )

        # Power densities for different particles [MW/m^3]
        # Neutron particle gets approximately 75% of the fusion power
        alpha_power_density = 0.0
        charged_power_density = (
            1.0 - constants.dd_neutron_energy_fraction
        ) * fusion_power_density
        neutron_power_density = (
            constants.dd_neutron_energy_fraction * fusion_power_density
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
            charged_power_density,
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
            - self.charged_power_density: Power density of charged particles produced.
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

        # Reaction energy in MegaJoules [MJ]
        reaction_energy = constants.dd_triton_energy / 1.0e6

        # Calculate the fusion power density produced [MW/m^3]
        fusion_power_density = (
            sigmav
            * reaction_energy
            * 0.5  # Factor for D-D reaction
            * (physics_variables.fdeut * physics_variables.deni)
            * (physics_variables.fdeut * physics_variables.deni)
        )

        # Power densities for different particles [MW/m^3]
        alpha_power_density = 0.0
        charged_power_density = fusion_power_density
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
            charged_power_density,
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
        self.charged_power_density += charged_power_add
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
        physics_variables.alpha_power_density = self.alpha_power_density
        physics_variables.charged_power_density = self.charged_power_density
        physics_variables.neutron_power_density = self.neutron_power_density
        physics_variables.fusion_rate_density = self.fusion_rate_density
        physics_variables.alpha_rate_density = self.alpha_rate_density
        physics_variables.proton_rate_density = self.proton_rate_density
        physics_module.sigmav_dt_average = self.sigmav_dt_average
        physics_module.dt_power_density = self.dt_power_density
        physics_module.dhe3_power_density = self.dhe3_power_density
        physics_module.dd_power_density = self.dd_power_density


def radpwr(plasma_profile):
    """This routine finds the radiation powers in MW/m3 by calling
    relevant routines.
    Author:
        P J Knight, CCFE, Culham Science Centre

    :param plasma_profile:
    :type plasma_profile: PlasmaProfile
    :return: RadpwrData - dataclass returning radiation power
    :rtype: DataClass
    """
    imp_rad = impurity.ImpurityRadiation(plasma_profile)
    imp_rad.calculate_imprad()

    pedgeradpv = imp_rad.radtot - imp_rad.radcore

    # Synchrotron radiation power/volume; assumed to be from core only.
    # Synchrotron radiation power/volume; assumed to be from core only.
    psyncpv = psync_albajar_fidone()

    # Total core radiation power/volume.
    # Total core radiation power/volume.
    pcoreradpv = imp_rad.radcore + psyncpv

    # Total radiation power/volume.
    # Total radiation power/volume.
    pradpv = imp_rad.radtot + psyncpv

    return RadpwrData(psyncpv, pcoreradpv, pedgeradpv, pradpv)


def psync_albajar_fidone():
    """This routine finds the synchrotron radiation power in MW/m3,
    using the method of Albajar and Fidone.
    References:
        Albajar, Nuclear Fusion 41 (2001) 665
        Fidone, Giruzzi, Granata, Nuclear Fusion 41 (2001) 1755
    Authors:
        P J Knight, CCFE, Culham Science Centre
        R Kemp, CCFE, Culham Science Centre

    :return: psyncpv synchrotron radiation power/volume (MW/m3)
    :rtype: float
    """
    tbet = 2.0e0

    # rpow is the (1-Rsyn) power dependence based on plasma shape
    # (see Fidone)
    # rpow is the (1-Rsyn) power dependence based on plasma shape
    # (see Fidone)

    rpow = 0.62e0

    kap = 0.0
    de2o = 0.0
    pao = 0.0
    gfun = 0.0
    kfun = 0.0
    dum = 0.0
    psync = 0.0
    psyncpv = 0.0

    kap = physics_variables.plasma_volume / (
        2.0e0 * numpy.pi**2 * physics_variables.rmajor * physics_variables.rminor**2
    )

    # No account is taken of pedestal profiles here, other than use of
    # the correct physics_variables.ne0 and physics_variables.te0...
    # No account is taken of pedestal profiles here, other than use of
    # the correct physics_variables.ne0 and physics_variables.te0...

    de2o = 1.0e-20 * physics_variables.ne0
    pao = 6.04e3 * (physics_variables.rminor * de2o) / physics_variables.bt
    gfun = 0.93e0 * (
        1.0e0
        + 0.85e0
        * numpy.exp(-0.82e0 * physics_variables.rmajor / physics_variables.rminor)
    )
    kfun = (physics_variables.alphan + 3.87e0 * physics_variables.alphat + 1.46e0) ** (
        -0.79e0
    )
    kfun = kfun * (1.98e0 + physics_variables.alphat) ** 1.36e0 * tbet**2.14e0
    kfun = kfun * (tbet**1.53e0 + 1.87e0 * physics_variables.alphat - 0.16e0) ** (
        -1.33e0
    )
    dum = (
        1.0e0
        + 0.12e0
        * (physics_variables.te0 / (pao**0.41e0))
        * (1.0e0 - physics_variables.ssync) ** 0.41e0
    )

    # Very high T modification, from Fidone
    # Very high T modification, from Fidone

    dum = dum ** (-1.51e0)

    psync = (
        3.84e-8
        * (1.0e0 - physics_variables.ssync) ** rpow
        * physics_variables.rmajor
        * physics_variables.rminor**1.38e0
    )
    psync = psync * kap**0.79e0 * physics_variables.bt**2.62e0 * de2o**0.38e0
    psync = (
        psync
        * physics_variables.te0
        * (16.0e0 + physics_variables.te0) ** 2.61e0
        * dum
        * gfun
        * kfun
    )

    # psyncpv should be per unit volume; Albajar gives it as total
    # psyncpv should be per unit volume; Albajar gives it as total

    psyncpv = psync / physics_variables.plasma_volume

    return psyncpv


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
) -> numpy.ndarray:
    """
    Evaluate the integrand for the fusion power integration.

    Parameters:
        plasma_profile (PlasmaProfile): Parameterised temperature and density profiles.
        reactionconstants (BoschHaleConstants): Bosch-Hale reaction constants.

    Returns:
        numpy.ndarray: Integrand for the fusion power.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,
          doi: https://doi.org/10.1088/0029-5515/32/4/i07.
    """

    # Since the electron temperature profile is only calculated directly, we scale the ion temperature
    # profile by the ratio of the volume averaged ion to electron temperature
    ion_temperature_profile = (
        physics_variables.ti / physics_variables.te
    ) * plasma_profile.teprofile.profile_y

    # Number of fusion reactions per unit volume per particle volume density (m^3/s)
    sigv = bosch_hale_reactivity(ion_temperature_profile, reaction_constants)

    # Integrand for the volume averaged fusion reaction rate sigmav:
    # sigmav = integral(2 rho (sigv(rho) ni(rho)^2) drho),
    # divided by the square of the volume-averaged ion density
    # to retain the dimensions m^3/s (this is multiplied back in later)

    # Set each point in the desnity profile as a fraction of the volume averaged desnity
    density_profile_normalised = (
        1.0 / physics_variables.dene
    ) * plasma_profile.neprofile.profile_y

    # Calculate a volume averaged fusion reaction integral that allows for fusion power to be scaled with
    # just the volume averged ion density.
    fusion_integral = (
        2.0
        * plasma_profile.teprofile.profile_x
        * sigv
        * density_profile_normalised**2
    )

    return fusion_integral


def bosch_hale_reactivity(
    ion_temperature_profile: numpy.ndarray, reaction_constants: BoschHaleConstants
) -> numpy.ndarray:
    """
    Calculate the volumetric fusion reaction rate <sigma v> (m^3/s) for one of four nuclear reactions using
    the Bosch-Hale parametrization.

    The valid range of the fit is 0.2 keV < t < 100 keV except for D-3He where it is 0.5 keV < t < 190 keV.

    Reactions:
        1. D-T reaction
        2. D-3He reaction
        3. D-D 1st reaction
        4. D-D 2nd reaction

    Parameters:
        ion_temperature_profile (numpy.ndarray): Plasma ion temperature profile in keV.
        reaction_constants (BoschHaleConstants): Bosch-Hale reaction constants.

    Returns:
        numpy.ndarray: Volumetric fusion reaction rate <sigma v> in m^3/s for each point in the ion temperature profile.

    References:
        - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
          Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,
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

    # Volumetric reaction rate / reactivity <sigma v> (m^3/s)
    # Original form is in [cm^3/s], so multiply by 1.0e-6 to convert to [m^3/s]
    sigmav = (
        1.0e-6
        * reaction_constants.cc1
        * theta
        * numpy.sqrt(xi / (reaction_constants.mrc2 * ion_temperature_profile**3))
        * numpy.exp(-3.0 * xi)
    )

    # if t = 0, sigmav = 0. Use this mask to set sigmav to zero.
    t_mask = ion_temperature_profile == 0.0
    sigmav[t_mask] = 0.0

    # Return np.ndarray of sigmav for each point in the ion temperature profile
    return sigmav


@dataclass
class RadpwrData:
    """DataClass which holds the output of the function radpwr"""

    psyncpv: float
    pcoreradpv: float
    pedgeradpv: float
    pradpv: float


def palph2(
    bp: float,
    bt: float,
    dene: float,
    deni: float,
    dnitot: float,
    falpe: float,
    falpi: float,
    alpha_power_beams: float,
    charged_power_density: float,
    neutron_power_density: float,
    ten: float,
    tin: float,
    plasma_volume: float,
    alpha_power_density: float,
    ifalphap: int,
) -> tuple:
    """
    Calculate fusion power and fast alpha pressure.

    This function computes various fusion power metrics and the fast alpha pressure
    based on the provided plasma parameters.

    Parameters:
        bp (float): Poloidal field (T).
        bt (float): Toroidal field on axis (T).
        dene (float): Electron density (m^-3).
        deni (float): Fuel ion density (m^-3).
        dnitot (float): Total ion density (m^-3).
        falpe (float): Fraction of alpha energy to electrons.
        falpi (float): Fraction of alpha energy to ions.
        alpha_power_beams (float): Alpha power from hot neutral beam ions (MW).
        charged_power_density (float): Other charged particle fusion power per unit volume (MW/m^3).
        neutron_power_density (float): Neutron fusion power per unit volume (MW/m^3).
        ten (float): Density-weighted electron temperature (keV).
        tin (float): Density-weighted ion temperature (keV).
        plasma_volume (float): Plasma volume (m^3).
        alpha_power_density (float): Alpha power per unit volume (MW/m^3).
        ifalphap (int): Switch for fast alpha pressure method.

    Returns:
        tuple: A tuple containing the following elements:
            - neutron_power_density_out (float): Neutron fusion power per unit volume [MW/m^3].
            - alpha_power_plasma (float): Alpha fusion power from only the plasma [MW].
            - alpha_power_total (float): Total alpha fusion power from plasma and beams [MW].
            - neutron_power_plasma (float): Neutron fusion power from only the plasma [MW].
            - neutron_power_total (float): Total neutron fusion power from plasma and beams [MW].
            - non_alpha_charged_power (float): Other total charged particle fusion power [MW].
            - betaft (float): Fast alpha beta component.
            - alpha_power_density_out (float): Alpha power per unit volume [MW/m^3].
            - alpha_power_electron_density (float): Alpha power per unit volume to electrons [MW/m^3].
            - alpha_power_ions_density (float): Alpha power per unit volume to ions [MW/m^3].
            - charged_particle_power(float): Charged particle fusion power [MW].
            - fusion_power (float): Total fusion power [MW].

    References:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990

    """
    # Alpha power

    # Calculate alpha power produced just by the plasma
    alpha_power_plasma = alpha_power_density * plasma_volume

    # Add neutral beam alpha power / volume
    alpha_power_density_out = alpha_power_density + (alpha_power_beams / plasma_volume)

    # Total alpha power
    alpha_power_total = alpha_power_density_out * plasma_volume

    # Neutron Power

    # Calculate neutron power produced just by the plasma
    neutron_power_plasma = neutron_power_density * plasma_volume

    # Add extra neutron power from beams
    neutron_power_density_out = neutron_power_density + (
        (4.0 * alpha_power_beams) / plasma_volume
    )

    # Total neutron power
    neutron_power_total = neutron_power_density_out * plasma_volume

    # Charged particle power

    # Total non-alpha charged particle power
    non_alpha_charged_power = charged_power_density * plasma_volume

    # Charged particle fusion power
    charged_particle_power = alpha_power_total + non_alpha_charged_power

    # Total fusion power
    fusion_power = alpha_power_total + neutron_power_total + non_alpha_charged_power

    # Alpha power to electrons and ions (used with electron
    # and ion power balance equations only)
    # No consideration of charged_power_density here...
    alpha_power_ions_density = (
        physics_variables.f_alpha_plasma * alpha_power_density_out * falpi
    )
    alpha_power_electron_density = (
        physics_variables.f_alpha_plasma * alpha_power_density_out * falpe
    )

    # Determine average fast alpha density
    if physics_variables.fdeut < 1.0:

        betath = (
            2.0e3
            * constants.rmu0
            * constants.electron_charge
            * (dene * ten + dnitot * tin)
            / (bt**2 + bp**2)
        )

        # jlion: This "fact" model is heavily flawed for smaller temperatures! It is unphysical for a stellarator (high n low T)
        # IPDG89 fast alpha scaling
        if ifalphap == 0:
            fact = min(0.3, 0.29 * (deni / dene) ** 2 * ((ten + tin) / 20.0 - 0.37))

        # Modified scaling, D J Ward
        else:
            fact = min(
                0.30,
                0.26
                * (deni / dene) ** 2
                * numpy.sqrt(max(0.0, ((ten + tin) / 20.0 - 0.65))),
            )

        fact = max(fact, 0.0)
        fact2 = alpha_power_density_out / alpha_power_density
        betaft = betath * fact * fact2

    else:  # negligible alpha production, alpha_power_density = alpha_power_beams = 0
        betaft = 0.0

    return (
        neutron_power_density_out,
        alpha_power_plasma,
        alpha_power_total,
        neutron_power_plasma,
        neutron_power_total,
        non_alpha_charged_power,
        betaft,
        alpha_power_density_out,
        alpha_power_electron_density,
        alpha_power_ions_density,
        charged_particle_power,
        fusion_power,
    )


def beamfus(
    beamfus0,
    betbm0,
    bp,
    bt,
    cnbeam,
    dene,
    deni,
    dlamie,
    enbeam,
    fdeut,
    ftrit,
    ftritbm,
    sigmav_dt_average,
    ten,
    tin,
    plasma_volume,
    zeffai,
):
    """Routine to calculate beam slowing down properties
    author: P J Knight, CCFE, Culham Science Centre

    :param beamfus0: multiplier for beam-background fusion calculation
    :param betbm0: leading coefficient for neutral beam beta fraction
    :param bp: poloidal field (T)
    :param bt: toroidal field on axis (T)
    :param cnbeam: neutral beam current (A)
    :param dene: electron density (m^-3)
    :param deni: fuel ion density (m^-3)
    :param dlamie: ion-electron coulomb logarithm
    :param enbeam: neutral beam energy (keV)
    :param fdeut: deuterium fraction of main plasma
    :param ftrit: tritium fraction of main plasma
    :param ftritbm: tritium fraction of neutral beam
    :param sigmav_dt_average: profile averaged <sigma v> for D-T (m3/s)
    :param ten: density-weighted electron temperature (keV)
    :param tin: density-weighted ion temperature (keV)
    :param plasma_volume: plasma volume (m3)
    :param zeffai: mass weighted plasma effective charge

    :returns: neutral beam beta component, hot beam ion density (m^-3),
    alpha power from hot neutral beam ions (MW)
    """

    tausl = (
        1.99e19
        * (2.0 * (1.0 - ftritbm) + (3.0 * ftritbm))
        * (ten**1.5 / dene)
        / dlamie
    )

    # Critical energy for electron/ion slowing down of the beam ion
    # (deuterium and tritium neutral beams, respectively) (keV)

    ecritd = 14.8 * ten * 2.0 * zeffai**0.6666 * (dlamie + 4.0) / dlamie
    ecritt = ecritd * 1.5

    # Deuterium and tritium ion densities

    denid = deni * fdeut
    denit = deni * ftrit

    palpdb, palptb, dnbeam2, ehotnb = beamcalc(
        denid,
        denit,
        enbeam,
        ecritd,
        ecritt,
        tausl,
        ftritbm,
        cnbeam,
        tin,
        plasma_volume,
        sigmav_dt_average,
    )

    # Neutral beam alpha power

    alpha_power_beams = beamfus0 * (palpdb + palptb)

    # Neutral beam beta

    betanb = betbm0 * 4.03e-22 * 0.66666 * dnbeam2 * ehotnb / (bt**2 + bp**2)

    return betanb, dnbeam2, alpha_power_beams


def beamcalc(
    nd, nt, ebeam, ecritd, ecritt, tausbme, ftritbm, ibeam, ti, plasma_volume, svdt
):
    """Neutral beam alpha power and ion energy
    author: P J Knight, CCFE, Culham Science Centre

    :param nd: thermal deuterium density (/m3)
    :param nt: thermal tritium density   (/m3)
    :param ebeam: beam energy (keV)
    :param ecritd: critical energy for electron/ion slowing down of
    the beam ion (deuterium neutral beam) (keV)
    :param ecritt: critical energy for beam slowing down (tritium neutral beam) (keV)
    :param ftritbm: beam tritium fraction (0.0 = deuterium beam)
    :param ibeam: beam current (A)
    :param svdt: profile averaged <sigma v> for D-T (m3/s)
    :param tausbme: beam ion slowing down time on electrons (s)
    :param ti: thermal ion temperature (keV)
    :param plasma_volume: plasma volume (m3) (95% flux surface)

    :returns: alpha power from deut. beam-background fusion (MW),
    alpha power from trit. beam-background fusion (MW), hot beam ion density (m^-3),
    average hot beam ion energy (keV)
    """

    # D and T beam current fractions
    ifbmd = ibeam * (1.0 - ftritbm)
    ifbmt = ibeam * ftritbm

    ebmratd = ebeam / ecritd
    vcritd = numpy.sqrt(
        2.0
        * constants.electron_charge
        * 1000.0
        * ecritd
        / (constants.proton_mass * ATOMIC_MASS_DEUTERIUM)
    )
    tauseffd = tausbme / 3.0 * numpy.log(1.0 + (ebmratd) ** 1.5)
    nhotmsd = (
        (1.0 - ftritbm) * ibeam * tauseffd / (constants.electron_charge * plasma_volume)
    )

    ebmratt = ebeam / ecritt
    vcritt = numpy.sqrt(
        2.0
        * constants.electron_charge
        * 1000.0
        * ecritt
        / (constants.proton_mass * ATOMIC_MASS_TRITIUM)
    )
    tausefft = tausbme / 3.0 * numpy.log(1.0 + (ebmratt) ** 1.5)
    nhotmst = ftritbm * ibeam * tausefft / (constants.electron_charge * plasma_volume)

    nhot = nhotmsd + nhotmst
    ndhot = nhotmsd
    nthot = nhotmst

    # Average hot ion energy from Deng & Emmert, UWFDM-718, Jan 87
    vcds = (
        2.0
        * ecritd
        * constants.electron_charge
        * 1000.0
        / (2.0 * constants.proton_mass)
    )
    vcts = (
        2.0
        * ecritt
        * constants.electron_charge
        * 1000.0
        / (3.0 * constants.proton_mass)
    )

    s0d = ifbmd / (constants.electron_charge * plasma_volume)
    s0t = ifbmt / (constants.electron_charge * plasma_volume)

    xcoefd = (
        ATOMIC_MASS_DEUTERIUM
        * constants.proton_mass
        * tausbme
        * vcds
        * s0d
        / (constants.electron_charge * 1000.0 * 3.0)
    )
    xcoeft = (
        ATOMIC_MASS_TRITIUM
        * constants.proton_mass
        * tausbme
        * vcts
        * s0t
        / (constants.electron_charge * 1000.0 * 3.0)
    )

    presd = xcoefd * xbrak(ebeam, ecritd)
    prest = xcoeft * xbrak(ebeam, ecritt)

    ehotd = 1.5 * presd / ndhot
    ehott = 1.5 * prest / nthot
    ehot = (ndhot * ehotd + nthot * ehott) / nhot

    iabm = 2
    svdhotn = 1e-4 * sgvhot(iabm, vcritd, ebeam)
    iabm = 3
    svthotn = 1e-4 * sgvhot(iabm, vcritt, ebeam)

    palfdb = palphabm(ndhot, nt, svdhotn, plasma_volume, ti, svdt)
    palftb = palphabm(nthot, nd, svthotn, plasma_volume, ti, svdt)

    return palfdb, palftb, nhot, ehot


def xbrak(e0, ec):
    """Hot ion energy parameter
    author: P J Knight, CCFE, Culham Science Centre

    :param e0: neutral beam energy (keV)
    :param ec: critical energy for electron/ion slowing down of the beam ion (keV)
    """
    xcs = e0 / ec
    xc = numpy.sqrt(xcs)

    t1 = xcs / 2.0
    t2 = (numpy.log((xcs + 2.0 * xc + 1.0) / (xcs - xc + 1.0))) / 6.0

    xarg = (2.0 * xc - 1.0) / numpy.sqrt(3.0)
    t3 = (numpy.arctan(xarg)) / numpy.sqrt(3.0)
    t4 = 0.3022999

    return t1 + t2 - t3 - t4


def palphabm(nbm, nblk, sigv, plasma_volume, ti, svdt):
    """Alpha power from beam-background fusion
    author: P J Knight, CCFE, Culham Science Centre

    :param nblk: thermal ion density (/m3)
    :param nbm: hot beam ion density (/m3)
    :param sigv: hot beam fusion reaction rate (m3/s)
    :param plasma_volume: plasma volume (m3)
    :param ti: thermal ion temperature (keV)
    :param svdt: profile averaged <sigma v> for D-T (m3/s)
    """

    # [ti] because bosch_hale expects temperature profile
    # so we pass it a profile of length 1
    ratio = svdt / bosch_hale_reactivity(
        numpy.array([ti]), BoschHaleConstants(**REACTION_CONSTANTS_DT)
    )
    return nbm * nblk * sigv * constants.dt_alpha_energy * plasma_volume * ratio.item()


def sgvhot(rmass_ion, vcrx, ebeam):
    """Hot beam fusion reaction rate
    author: P J Knight, CCFE, Culham Science Centre

    :param rmass_ion: relative atomic mass of the ion (of D or T)
    :param vcrx: critical velocity for electron/ion slowing down of the beam ion (m/s)
    :param ebeam: neutral beam energy (keV)
    """
    # Beam velocity

    vbeams = (
        ebeam
        * constants.electron_charge
        * 1000.0
        * 2.0
        / (rmass_ion * constants.proton_mass)
    )
    vbeam = numpy.sqrt(vbeams)

    xv = vbeam / vcrx
    t1 = 3.0 * vcrx / numpy.log(1.0 + (xv**3))

    svint = integrate.quad(
        _hot_beam_fusion_reaction_rate_integrand, 0.0, xv, args=(vcrx,)
    )[0]

    return t1 * svint


def _hot_beam_fusion_reaction_rate_integrand(u, vcritx):
    """Integrand function for the hot beam fusion reaction rate
    author: P J Knight, CCFE, Culham Science Centre

    :param u: ratio of beam velocity to the critical velocity
    """
    t1 = (u**3) / (1.0 + u**3)

    # vcritx : critical velocity for electron/ion slowing down of beam ion (m/s)
    xvc = vcritx * u
    xvcs = xvc * xvc * constants.proton_mass / (constants.electron_charge * 1000.0)
    t2 = _sigbmfus(xvcs)

    return t1 * t2


def _sigbmfus(vrelsq):
    """Fusion reaction cross-section
    author: P J Knight, CCFE, Culham Science Centre

    The functional form of the cross-section is in terms of the equivalent
    deuterium energy, i.e. for a tritium beam at 500 keV the energy
    used in the cross-section function is 333 keV.

    :param vrelsq: square of the speed of the beam ion (keV/amu)
    """
    a1 = 45.95
    a2 = 5.02e4
    a3 = 1.368e-2
    a4 = 1.076
    a5 = 4.09e2

    # Beam kinetic energy

    ebm = 0.5 * ATOMIC_MASS_DEUTERIUM * vrelsq

    # Set limits on cross-section at low and high beam energies

    if ebm < 10.0:
        return 1.0e-27
    elif ebm > 1.0e4:
        return 8.0e-26
    else:
        t1 = a2 / (1.0e0 + (a3 * ebm - a4) ** 2) + a5
        t2 = ebm * (numpy.exp(a1 / numpy.sqrt(ebm)) - 1.0)
        return 1.0e-24 * t1 / t2
