"""Module for plasma profile definitions and utilities.

This module provides abstract base classes and implementations for creating
and managing plasma profiles such as temperature and density distributions.
"""

import logging
from abc import ABC, abstractmethod
from enum import IntEnum
from types import DynamicClassAttribute

import numpy as np
import scipy as sp

from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.models.physics.density_limit import PlasmaDensityLimit

logger = logging.getLogger(__name__)


class PlasmaProfileShapeType(IntEnum):
    """Enum for i_plasma_pedestal method types"""

    PARABOLIC_PROFILE = (0, "Parabolic Profile (L-mode)")
    PEDESTAL_PROFILE = (1, "Pedestal Profile (H-mode)")

    def __new__(cls, value: int, description: str):
        """Create a new PlasmaProfileShapeType instance."""
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._description_ = description
        return obj

    @DynamicClassAttribute
    def description(self):
        """The description of the plasma profile shape."""
        return self._description_


class Profile(Model, ABC):
    """Abstract base class used to create and hold profiles (temperature, density)"""

    def __init__(self):
        """
        Initialize a Profiles object.

        Attributes
        ----------
        - profile_size (int): The size of the profile.
        - profile_x (ndarray): An array of values ranging from 0 to profile_size-1.
        - profile_y (ndarray): An array of zeros with length profile_size.
        - profile_integ (int): The integral of the profile_y array.
        - profile_dx (int): The step size between consecutive values in profile_x.
        """
        self.profile_integ = 0
        self.profile_dx = 0

    def run(self):
        """Initialise profile_x and profile_y"""
        self.profile_x = np.arange(
            self.data.physics.n_plasma_profile_elements, dtype=float
        )
        self.profile_y = np.zeros(self.data.physics.n_plasma_profile_elements)

    def output(self):
        """Profile model doesn't have any output"""

    def normalise_profile_x(self):
        """Normalizes the x-dimension of the profile.

        This method divides the values in the `profile_x` attribute by the maximum value
        in the `profile_x` array, resulting in a normalized version of the x-dimension.

        Example:
            If `profile_x` is [1, 2, 3, 4, 5], after normalization it will become
            [0.2, 0.4, 0.6, 0.8, 1.0].

        Note:
            This method modifies the `profile_x` attribute in-place.


        """
        self.profile_x /= max(self.profile_x)

    def calculate_profile_dx(self):
        """Calculates the differential between points in the profile.

        This method calculates the differential between points in the profile by
        dividing the maximum x value in the profile by the difference in size between
        the points. The result is stored in the `profile_dx` attribute.
        """
        self.profile_dx = max(self.profile_x) / (
            self.data.physics.n_plasma_profile_elements - 1
        )

    @abstractmethod
    def calculate_profile_y(self):
        """Use a profile function to act on self.profile_x to calculate and set the
        values of self.profile_y.
        """

    def integrate_profile_y(self):
        """Integrate profile_y values using scipy.integrate.simpson() function.

        This method calculates the integral of the profile_y values using the Simpson's
        rule provided by the scipy.integrate.simpson() function. The integral is stored
        in the `profile_integ` attribute.
        """
        self.profile_integ = sp.integrate.simpson(
            self.profile_y, x=self.profile_x, dx=self.profile_dx
        )


class DensityProfilePedestalType(IntEnum):
    """Enum for i_nd_plasma_pedestal_separatrix types"""

    USER_INPUT = (0, "User input direct values")
    GREENWALD_FRACTION = (1, "Fractions of the Greenwald limit")

    def __new__(cls, value: int, description: str):
        """Create a new DensityProfilePedestalType instance.

        Parameters
        ----------
            value: Integer value for the enum member.
            description: Human-readable description for the enum member.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._description_ = description
        return obj

    @DynamicClassAttribute
    def description(self):
        """The description of the plasma profile shape."""
        return self._description_


class ElectronDensityProfile(Profile):
    """Electron density (nₑ) profile class. Contains a function to calculate the electron
    density profile and store the data.
    """

    def run(self):
        """Subroutine which calls profile functions and stores neprofile data."""
        super().run()
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            rho=self.profile_x,
            radius_plasma_pedestal_density_norm=self.data.physics.radius_plasma_pedestal_density_norm,
            nd_on_axis=self.data.physics.nd_plasma_electron_on_axis,
            nd_pedestal=self.data.physics.nd_plasma_pedestal_electron,
            nd_separatrix=self.data.physics.nd_plasma_separatrix_electron,
            alphan=self.data.physics.alphan,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        radius_plasma_pedestal_density_norm: float,
        nd_on_axis: float,
        nd_pedestal: float,
        nd_separatrix: float,
        alphan: float,
    ) -> None:
        """Calculates the number density at each normalised minor radius (ρ) position.

        Parameters
        ----------
        rho :
            Normalised minor radius (ρ) vector.
        radius_plasma_pedestal_density_norm :
            Normalised minor radius pedestal position (ρₙ,pedestal).
        nd_on_axis :
            Central number density (n₀) [m⁻³].
        nd_pedestal :
            Pedestal density (n_pedestal) [m⁻³].
        nd_separatrix :
            Separatrix density (n_sep) [m⁻³].
        alphan :
            Density peaking parameter (αₙ).
        """  # noqa: RUF002
        if (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PARABOLIC_PROFILE
        ):
            self.profile_y = nd_on_axis * (1 - rho**2) ** alphan

        # Input checks

        if nd_on_axis < nd_pedestal:
            logger.info(
                "NPROFILE: Pedestal density is higher than core density. %s, %s",
                nd_pedestal,
                nd_on_axis,
            )
        rho_index = rho <= radius_plasma_pedestal_density_norm
        self.profile_y[rho_index] = (
            nd_pedestal
            + (nd_on_axis - nd_pedestal)
            * (1 - (rho[rho_index] / radius_plasma_pedestal_density_norm) ** 2) ** alphan
        )
        # Invert the rho_index
        self.profile_y[~rho_index] = nd_separatrix + (nd_pedestal - nd_separatrix) * (
            1 - rho[~rho_index]
        ) / (1 - radius_plasma_pedestal_density_norm)

    @staticmethod
    def calculate_pedestal_profile_on_axis_density(
        radius_plasma_pedestal_density_norm: float,
        nd_pedestal: float,
        nd_separatrix: float,
        nd_vol_average: float,
        alphan: float,
    ) -> float:
        """Calculates the core density (n₀) of a pedestalised profile.

        Parameters
        ----------
        radius_plasma_pedestal_density_norm :
            Normalised minor radius pedestal position (ρₙ,pedestal).
        nd_pedestal: float,
            The pedestal density (n_pedestal) [m⁻³].
        nd_separatrix: float,
            The separatrix density (n_sep) [m⁻³].
        nd_vol_average: float,
            The volume averaged density (⟨n⟩) [m⁻³].
        alphan: float,
            The density peaking parameter (αₙ).

        Returns
        -------
        :
            The core on-axis density (n₀) [/m³].
        """
        nd_on_axis = (
            1
            / (3 * radius_plasma_pedestal_density_norm**2)
            * (
                3 * nd_vol_average * (1 + alphan)
                + nd_separatrix
                * (1 + alphan)
                * (
                    -2
                    + radius_plasma_pedestal_density_norm
                    + radius_plasma_pedestal_density_norm**2
                )
                - nd_pedestal
                * (
                    (1 + alphan) * (1 + radius_plasma_pedestal_density_norm)
                    + (alphan - 2) * radius_plasma_pedestal_density_norm**2
                )
            )
        )

        if nd_on_axis < 0.0:
            # Allows solver to continue and
            # warns the user to raise the lower bound on nd_plasma_electrons_vol_avg
            # if the run did not converge
            logger.error(
                "nd_on_axis is going negative when solving. Please raise the value of "
                "nd_plasma_electrons_vol_avg (⟨nₑ⟩) and or its lower limit."
            )
            nd_on_axis = 1.0e-6
        return nd_on_axis

    def set_pedestal_and_separatrix_values(self):
        """Sets the pedestal and separatrix density values based on the user input
        or greenwald fraction method.
        """
        i_nd_plasma_pedestal_separatrix = DensityProfilePedestalType(
            self.data.physics.i_nd_plasma_pedestal_separatrix
        )

        if i_nd_plasma_pedestal_separatrix == DensityProfilePedestalType.USER_INPUT:
            self.data.physics.f_nd_plasma_pedestal_greenwald = (
                self.data.physics.nd_plasma_pedestal_electron
                / (
                    PlasmaDensityLimit.calculate_greenwald_density_limit(
                        c_plasma=self.data.physics.plasma_current,
                        rminor=self.data.physics.rminor,
                    )
                )
            )

            self.data.physics.f_nd_plasma_separatrix_greenwald = (
                self.data.physics.nd_plasma_separatrix_electron
                / (
                    PlasmaDensityLimit.calculate_greenwald_density_limit(
                        c_plasma=self.data.physics.plasma_current,
                        rminor=self.data.physics.rminor,
                    )
                )
            )
        elif (
            i_nd_plasma_pedestal_separatrix
            == DensityProfilePedestalType.GREENWALD_FRACTION
        ):
            self.data.physics.nd_plasma_pedestal_electron = (
                self.data.physics.f_nd_plasma_pedestal_greenwald
                * PlasmaDensityLimit.calculate_greenwald_density_limit(
                    c_plasma=self.data.physics.plasma_current,
                    rminor=self.data.physics.rminor,
                )
            )
            self.data.physics.nd_plasma_separatrix_electron = (
                self.data.physics.f_nd_plasma_separatrix_greenwald
                * PlasmaDensityLimit.calculate_greenwald_density_limit(
                    c_plasma=self.data.physics.plasma_current,
                    rminor=self.data.physics.rminor,
                )
            )

    def set_physics_variables(self):
        """Calculates and sets physics variables required for the profile."""
        if (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PARABOLIC_PROFILE
        ):
            self.data.physics.nd_plasma_electron_on_axis = (
                self.data.physics.nd_plasma_electrons_vol_avg
                * (1.0 + self.data.physics.alphan)
            )
        elif (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PEDESTAL_PROFILE
        ):
            self.data.physics.nd_plasma_electron_on_axis = self.calculate_pedestal_profile_on_axis_density(
                radius_plasma_pedestal_density_norm=self.data.physics.radius_plasma_pedestal_density_norm,
                nd_pedestal=self.data.physics.nd_plasma_pedestal_electron,
                nd_separatrix=self.data.physics.nd_plasma_separatrix_electron,
                nd_vol_avg=self.data.physics.nd_plasma_electrons_vol_avg,
                alphan=self.data.physics.alphan,
            )
        self.data.physics.nd_plasma_ions_on_axis = (
            self.data.physics.nd_plasma_ions_total_vol_avg
            / self.data.physics.nd_plasma_electrons_vol_avg
            * self.data.physics.nd_plasma_electron_on_axis
        )


class ElectronTemperatureProfile(Profile):
    """Electron temperature (Tₑ) profile class. Contains a function to calculate the
    temperature profile and store the data.
    """

    def run(self):
        """Subroutine to initialise neprofile and execute calculations."""
        super().run()
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            rho=self.profile_x,
            radius_plasma_pedestal_temp_norm=self.data.physics.radius_plasma_pedestal_temp_norm,
            temp_on_axis_kev=self.data.physics.temp_plasma_electron_on_axis_kev,
            temp_pedestal_kev=self.data.physics.temp_plasma_pedestal_kev,
            temp_separatrix_kev=self.data.physics.temp_plasma_separatrix_kev,
            alphat=self.data.physics.alphat,
            tbeta=self.data.physics.tbeta,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        radius_plasma_pedestal_temp_norm: float,
        temp_on_axis_kev: float,
        temp_pedestal_kev: float,
        temp_separatrix_kev: float,
        alphat: float,
        tbeta: float,
    ) -> None:
        """Calculates the temperature at each normalised minor radius (ρ) position.

        Parameters
        ----------
        rho :
            Normalised minor radius (ρ) vector
        radius_plasma_pedestal_temp_norm :
            Normalised minor radius pedestal position (ρₜ,pedestal).
        temp_on_axis_kev :
            Central on-axis temperature (T₀) [keV].
        temp_pedestal_kev :
            Pedestal temperature (T_pedestal) [keV].
        temp_separatrix_kev :
            Separatrix temperature (T_separatrix) [keV].
        alphat :
            Temperature peaking parameter (αₜ).
        tbeta :
            Second temperature exponent (βₜ).

        Raises
        ------
        ProcessValueError
            If negative temperature in plasma profile

        """  # noqa: RUF002
        if (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PARABOLIC_PROFILE
        ):
            # profile values of 0 cause divide by 0 errors so ensure the profile value
            # is at least 1e-8
            # which is small enough that it won't make a difference to any calculations
            self.profile_y = np.maximum(temp_on_axis_kev * (1 - rho**2) ** alphat, 1e-8)
            return

        if temp_on_axis_kev < temp_pedestal_kev:
            logger.info(
                "TPROFILE: Pedestal temperature is higher than core temperature. %s, %s",
                temp_pedestal_kev,
                temp_on_axis_kev,
            )

        rho_index = rho <= radius_plasma_pedestal_temp_norm
        self.profile_y[rho_index] = (
            temp_pedestal_kev
            + (temp_on_axis_kev - temp_pedestal_kev)
            * (1 - (rho[rho_index] / radius_plasma_pedestal_temp_norm) ** tbeta)
            ** alphat
        )
        self.profile_y[~rho_index] = temp_separatrix_kev + (
            temp_pedestal_kev - temp_separatrix_kev
        ) * (1 - rho[~rho_index]) / (1 - radius_plasma_pedestal_temp_norm)

        # Check for any negative temperature in profile: always fatal in
        # later models eventually
        if (self.profile_y < 0).any():
            raise ProcessValueError("Negative temperature in plasma profile")

    @staticmethod
    def calculate_pedestal_profile_on_axis_temperature(
        radius_plasma_pedestal_temp_norm: float,
        temp_pedestal_kev: float,
        temp_separatrix_kev: float,
        temp_vol_avg_kev: float,
        alphat: float,
        tbeta: float,
    ) -> float:
        """Calculates the core density (T₀) of a pedestalised profile.

        Parameters
        ----------
        radius_plasma_pedestal_temp_norm :
            Normalised minor radius pedestal position (ρₜ,pedestal).
        temp_pedestal_kev :
            Pedestal temperature (T_pedestal) [keV].
        temp_separatrix_kev :
            Separatrix temperature (T_separatrix) [keV].
        temp_vol_avg_kev :
            Volume average temperature (⟨T⟩) [keV].
        alphat :
            Temperature peaking parameter (αₜ).
        tbeta :
            Second temperature exponent (βₜ).

        Returns
        -------
        :
            The core on-axis temperature (T₀) [keV]

        """
        #  Calculate core temperature

        return temp_pedestal_kev + (
            (
                tbeta
                * (
                    3 * temp_vol_avg_kev
                    + temp_separatrix_kev
                    * (
                        -2.0
                        + radius_plasma_pedestal_temp_norm
                        + radius_plasma_pedestal_temp_norm**2
                    )
                    - temp_pedestal_kev
                    * (
                        1
                        + radius_plasma_pedestal_temp_norm
                        + radius_plasma_pedestal_temp_norm**2
                    )
                )
            )
            / (
                6
                * radius_plasma_pedestal_temp_norm**2
                * sp.special.beta(1 + alphat, 2 / tbeta)
            )
        )

    def set_physics_variables(self):
        """Calculates and sets physics variables required for the temperature profile."""
        if (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PARABOLIC_PROFILE
        ):
            self.data.physics.temp_plasma_electron_on_axis_kev = (
                self.data.physics.temp_plasma_electron_vol_avg_kev
                * (1.0 + self.data.physics.alphat)
            )
        elif (
            PlasmaProfileShapeType(self.data.physics.i_plasma_pedestal)
            == PlasmaProfileShapeType.PEDESTAL_PROFILE
        ):
            self.data.physics.temp_plasma_electron_on_axis_kev = (
                self.calculate_pedestal_profile_on_axis_temperature(
                    self.data.physics.radius_plasma_pedestal_temp_norm,
                    self.data.physics.temp_plasma_pedestal_kev,
                    self.data.physics.temp_plasma_separatrix_kev,
                    self.data.physics.temp_plasma_electron_vol_avg_kev,
                    self.data.physics.alphat,
                    self.data.physics.tbeta,
                )
            )

        self.data.physics.temp_plasma_ion_on_axis_kev = (
            self.data.physics.temp_plasma_ion_vol_avg_kev
            / self.data.physics.temp_plasma_electron_vol_avg_kev
            * self.data.physics.temp_plasma_electron_on_axis_kev
        )
