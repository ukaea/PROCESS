import logging
from abc import ABC, abstractmethod

import numpy as np
import scipy as sp

from process.data_structure import physics_variables

logger = logging.getLogger(__name__)


class Profile(ABC):
    """Abstract base class used to create and hold profiles (temperature, density)"""

    def __init__(self, profile_size: int):
        """
        Initialize a Profiles object.

        Parameters:
        - profile_size (int): The size of the profile.

        Attributes:
        - profile_size (int): The size of the profile.
        - profile_x (ndarray): An array of values ranging from 0 to profile_size-1.
        - profile_y (ndarray): An array of zeros with length profile_size.
        - profile_integ (int): The integral of the profile_y array.
        - profile_dx (int): The step size between consecutive values in profile_x.
        """
        self.profile_size = profile_size
        self.profile_x = np.arange(self.profile_size)
        self.profile_y = np.zeros(self.profile_size)
        self.profile_integ = 0
        self.profile_dx = 0

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
        self.profile_x = self.profile_x / max(self.profile_x)

    def calculate_profile_dx(self):
        """Calculates the differential between points in the profile.

        This method calculates the differential between points in the profile by dividing the maximum x value in the profile
        by the difference in size between the points. The result is stored in the `profile_dx` attribute.
        """
        self.profile_dx = max(self.profile_x) / (self.profile_size - 1)

    @abstractmethod
    def calculate_profile_y(self):
        """Use a profile function to act on self.profile_x to calculate and set the
        values of self.profile_y.
        """

    def integrate_profile_y(self):
        """Integrate profile_y values using scipy.integrate.simpson() function.

        This method calculates the integral of the profile_y values using the Simpson's rule
        provided by the scipy.integrate.simpson() function. The integral is stored in the
        self.profile_integ attribute.
        """
        self.profile_integ = sp.integrate.simpson(
            self.profile_y, x=self.profile_x, dx=self.profile_dx
        )


class NeProfile(Profile):
    """Electron density profile class. Contains a function to calculate the electron density profile and
    store the data.
    """

    def run(self):
        """Subroutine which calls profile functions and stores neprofile data."""
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            self.profile_x,
            physics_variables.radius_plasma_pedestal_density_norm,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.nd_plasma_pedestal_electron,
            physics_variables.nd_plasma_separatrix_electron,
            physics_variables.alphan,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        radius_plasma_pedestal_density_norm: float,
        n0: float,
        nped: float,
        nsep: float,
        alphan: float,
    ):
        """This routine calculates the density at each normalised minor radius position
        rho for a HELIOS-type density pedestal profile (neprofile).

        Parameters
        ----------
        rho :
            Normalised minor radius vector.
        radius_plasma_pedestal_density_norm :
            Normalised minor radius pedestal position.
        n0 :
            Central density (/m3).
        nped :
            Pedestal density (/m3).
        nsep :
            Separatrix density (/m3)
        alphan :
            Density peaking parameter.
        """

        if physics_variables.i_plasma_pedestal == 0:
            self.profile_y = n0 * (1 - rho**2) ** alphan

        # Input checks

        if n0 < nped:
            logger.info(
                f"NPROFILE: density pedestal is higher than core density. {nped = }, {n0 = }"
            )
        rho_index = rho <= radius_plasma_pedestal_density_norm
        self.profile_y[rho_index] = (
            nped
            + (n0 - nped)
            * (1 - (rho[rho_index] / radius_plasma_pedestal_density_norm) ** 2) ** alphan
        )
        # Invert the rho_index
        self.profile_y[~rho_index] = nsep + (nped - nsep) * (1 - rho[~rho_index]) / (
            1 - radius_plasma_pedestal_density_norm
        )

    @staticmethod
    def ncore(
        radius_plasma_pedestal_density_norm: float,
        nped: float,
        nsep: float,
        nav: float,
        alphan: float,
    ) -> float:
        """This routine calculates the core density of a pedestalised profile.
        The solution comes from integrating and summing the two separate density profiles for the core
        and pedestal region within their bounds. This has to be multiplied by the torus volume element before integration which leads
        to an added rho term in each part of the profile. When dividing by the volume of integration to get the average density
        the simplification leads to a factor of 2 having to be multiplied on to each of the integration results.
        This function for the average density can then be re-arranged to calculate the central plasma density n_0 / ncore.

        Parameters
        ----------
        radius_plasma_pedestal_density_norm :
             The normalised minor radius pedestal position.
        nped :
            The pedestal density (/m3).
        nsep :
            The separatrix density (/m3).
        nav :
            The electron density (/m3).
        alphan :
            The density peaking parameter

        Returns
        -------
        :
            The core density.

        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650
        """

        ncore = (
            1
            / (3 * radius_plasma_pedestal_density_norm**2)
            * (
                3 * nav * (1 + alphan)
                + nsep
                * (1 + alphan)
                * (
                    -2
                    + radius_plasma_pedestal_density_norm
                    + radius_plasma_pedestal_density_norm**2
                )
                - nped
                * (
                    (1 + alphan) * (1 + radius_plasma_pedestal_density_norm)
                    + (alphan - 2) * radius_plasma_pedestal_density_norm**2
                )
            )
        )

        if ncore < 0.0:
            # Allows solver to continue and
            # warns the user to raise the lower bound on nd_plasma_electrons_vol_avg if the run did not converge
            logger.error(
                "ncore is going negative when solving. Please raise the value of nd_plasma_electrons_vol_avg and or its lower limit."
            )
            ncore = 1.0e-6
        return ncore

    def set_physics_variables(self):
        """Calculates and sets physics variables required for the profile."""

        if physics_variables.i_plasma_pedestal == 0:
            physics_variables.nd_plasma_electron_on_axis = (
                physics_variables.nd_plasma_electrons_vol_avg
                * (1.0 + physics_variables.alphan)
            )
        elif physics_variables.i_plasma_pedestal == 1:
            physics_variables.nd_plasma_electron_on_axis = self.ncore(
                physics_variables.radius_plasma_pedestal_density_norm,
                physics_variables.nd_plasma_pedestal_electron,
                physics_variables.nd_plasma_separatrix_electron,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.alphan,
            )
        physics_variables.nd_plasma_ions_on_axis = (
            physics_variables.nd_plasma_ions_total_vol_avg
            / physics_variables.nd_plasma_electrons_vol_avg
            * physics_variables.nd_plasma_electron_on_axis
        )


class TeProfile(Profile):
    """Electron temperature profile class. Contains a function to calculate the temperature profile and store the data."""

    def run(self):
        """Subroutine to initialise neprofile and execute calculations."""
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            self.profile_x,
            physics_variables.radius_plasma_pedestal_temp_norm,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.temp_plasma_pedestal_kev,
            physics_variables.temp_plasma_separatrix_kev,
            physics_variables.alphat,
            physics_variables.tbeta,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        radius_plasma_pedestal_temp_norm: float,
        t0: float,
        temp_plasma_pedestal_kev: float,
        temp_plasma_separatrix_kev: float,
        alphat: float,
        tbeta: float,
    ):
        """Calculates the temperature at a normalised minor radius position rho for a pedestalised profile (teprofile).
        If i_plasma_pedestal = 0 the original parabolic profile form is used instead.

        Parameters
        ----------
        rho : np.array
            Normalised minor radius.
        radius_plasma_pedestal_temp_norm : float
            Normalised minor radius pedestal position.
        t0 : float
            Central temperature (keV).
        temp_plasma_pedestal_kev : float
            Pedestal temperature (keV).
        temp_plasma_separatrix_kev : float
            Separatrix temperature (keV).
        alphat : float
            Temperature peaking parameter.
        tbeta : float
            Second temperature exponent.

        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650
        """
        if physics_variables.i_plasma_pedestal == 0:
            # profile values of 0 cause divide by 0 errors so ensure the profile value is at least 1e-8
            # which is small enough that it won't make a difference to any calculations
            self.profile_y = np.maximum(t0 * (1 - rho**2) ** alphat, 1e-8)
            return

        if t0 < temp_plasma_pedestal_kev:
            logger.info(
                f"TPROFILE: temperature pedestal is higher than core temperature. {temp_plasma_pedestal_kev = }, {t0 = }"
            )

        rho_index = rho <= radius_plasma_pedestal_temp_norm
        self.profile_y[rho_index] = (
            temp_plasma_pedestal_kev
            + (t0 - temp_plasma_pedestal_kev)
            * (1 - (rho[rho_index] / radius_plasma_pedestal_temp_norm) ** tbeta)
            ** alphat
        )
        self.profile_y[~rho_index] = temp_plasma_separatrix_kev + (
            temp_plasma_pedestal_kev - temp_plasma_separatrix_kev
        ) * (1 - rho[~rho_index]) / (1 - radius_plasma_pedestal_temp_norm)

    @staticmethod
    def tcore(
        radius_plasma_pedestal_temp_norm: float,
        temp_plasma_pedestal_kev: float,
        temp_plasma_separatrix_kev: float,
        tav: float,
        alphat: float,
        tbeta: float,
    ) -> float:
        """This routine calculates the core temperature (keV)
        of a pedestalised profile. The solution comes from integrating and summing the two seprate temperature profiles for the core
        and pedestal region within their bounds. This has to be multiplied by the torus volume element before integration which leads
        to an added rho term in each part of the profile. When dividing by the volume of integration to get the average temperature
        the simplification leads to a factor of 2 having to be multiplied on to each of the integration results.
        This function for the average temperature can then be re-arranged to calculate the central plasma temeprature T_0 / tcore.
        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650

        Parameters
        ----------
        radius_plasma_pedestal_temp_norm : float
            Normalised minor radius pedestal position.
        temp_plasma_pedestal_kev : float
            Pedestal temperature (keV).
        temp_plasma_separatrix_kev : float
            Separatrix temperature (keV).
        tav : float
            Volume average temperature (keV).
        alphat : float
            Temperature peaking parameter.
        tbeta : float
            Second temperature exponent.

        Returns
        -------
        float
            Core temperature.
        """
        #  Calculate core temperature

        return temp_plasma_pedestal_kev + (
            (
                tbeta
                * (
                    3 * tav
                    + temp_plasma_separatrix_kev
                    * (
                        -2.0
                        + radius_plasma_pedestal_temp_norm
                        + radius_plasma_pedestal_temp_norm**2
                    )
                    - temp_plasma_pedestal_kev
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
        if physics_variables.i_plasma_pedestal == 0:
            physics_variables.temp_plasma_electron_on_axis_kev = (
                physics_variables.temp_plasma_electron_vol_avg_kev
                * (1.0 + physics_variables.alphat)
            )
        elif physics_variables.i_plasma_pedestal == 1:
            physics_variables.temp_plasma_electron_on_axis_kev = self.tcore(
                physics_variables.radius_plasma_pedestal_temp_norm,
                physics_variables.temp_plasma_pedestal_kev,
                physics_variables.temp_plasma_separatrix_kev,
                physics_variables.temp_plasma_electron_vol_avg_kev,
                physics_variables.alphat,
                physics_variables.tbeta,
            )

        physics_variables.temp_plasma_ion_on_axis_kev = (
            physics_variables.temp_plasma_ion_vol_avg_kev
            / physics_variables.temp_plasma_electron_vol_avg_kev
            * physics_variables.temp_plasma_electron_on_axis_kev
        )
