import logging
from abc import ABC, abstractmethod

import numpy as np
import scipy as sp

from process.fortran import error_handling, physics_variables

logger = logging.getLogger(__name__)
# Logging handler for console output
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.ERROR)
logger.addHandler(s_handler)


class Profile(ABC):
    """Abstract base class used to create and hold profiles (temperature, density)"""

    def __init__(self, profile_size: int) -> None:
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

    def normalise_profile_x(self) -> None:
        """
        Normalizes the x-dimension of the profile.

        This method divides the values in the `profile_x` attribute by the maximum value
        in the `profile_x` array, resulting in a normalized version of the x-dimension.

        Example:
            If `profile_x` is [1, 2, 3, 4, 5], after normalization it will become
            [0.2, 0.4, 0.6, 0.8, 1.0].

        Note:
            This method modifies the `profile_x` attribute in-place.

        Returns:
            None
        """
        self.profile_x = self.profile_x / max(self.profile_x)

    def calculate_profile_dx(self) -> None:
        """Calculates the differential between points in the profile.

        This method calculates the differential between points in the profile by dividing the maximum x value in the profile
        by the difference in size between the points. The result is stored in the `profile_dx` attribute.

        """
        self.profile_dx = max(self.profile_x) / (self.profile_size - 1)

    @abstractmethod
    def calculate_profile_y(self) -> None:
        """Use a profile function to act on self.profile_x to calculate and set the
        values of self.profile_y.
        """

    def integrate_profile_y(self) -> None:
        """
        Integrate profile_y values using scipy.integrate.simpson() function.

        This method calculates the integral of the profile_y values using the Simpson's rule
        provided by the scipy.integrate.simpson() function. The integral is stored in the
        self.profile_integ attribute.

        Parameters:
        None

        Returns:
        None
        """
        self.profile_integ = sp.integrate.simpson(
            self.profile_y, x=self.profile_x, dx=self.profile_dx
        )


class NeProfile(Profile):
    """Electron density profile class. Contains a function to calculate the electron density profile and
    store the data.

    Attributes:
        Inherits attributes from the base class `Profile`.

    Methods:
        run(): Subroutine which calls functions and stores neprofile data.
        calculate_profile_y(rho, rhopedn, n0, nped, nsep, alphan): Calculates the density at each normalised minor radius position.
        ncore(rhopedn, nped, nsep, nav, alphan): Calculates the core density of a pedestalised profile.
        set_physics_variables(): Calculates and sets physics variables required for the profile.
    """

    def run(self) -> None:
        """Subroutine which calls profile functions and stores neprofile data."""
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            self.profile_x,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        rhopedn: float,
        n0: float,
        nped: float,
        nsep: float,
        alphan: float,
    ) -> None:
        """
        This routine calculates the density at each normalised minor radius position
        rho for a HELIOS-type density pedestal profile (neprofile).

        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre

        Parameters:
            - rho (np.array): Normalised minor radius vector.
            - rhopedn (float): Normalised minor radius pedestal position.
            - n0 (float): Central density (/m3).
            - nped (float): Pedestal density (/m3).
            - nsep (float): Separatrix density (/m3).
            - alphan (float): Density peaking parameter.

        Returns:
            None
        """

        if physics_variables.ipedestal == 0:
            self.profile_y = n0 * (1 - rho**2) ** alphan

        # Input checks

        if n0 < nped:
            logger.info(
                f"NPROFILE: density pedestal is higher than core density. {nped = }, {n0 = }"
            )
        rho_index = rho <= rhopedn
        self.profile_y[rho_index] = (
            nped + (n0 - nped) * (1 - (rho[rho_index] / rhopedn) ** 2) ** alphan
        )
        # Invert the rho_index
        self.profile_y[~rho_index] = nsep + (nped - nsep) * (1 - rho[~rho_index]) / (
            1 - rhopedn
        )

    @staticmethod
    def ncore(
        rhopedn: float, nped: float, nsep: float, nav: float, alphan: float
    ) -> float:
        """
        This routine calculates the core density of a pedestalised profile.
        The solution comes from integrating and summing the two separate density profiles for the core
        and pedestal region within their bounds. This has to be multiplied by the torus volume element before integration which leads
        to an added rho term in each part of the profile. When dividing by the volume of integration to get the average density
        the simplification leads to a factor of 2 having to be multiplied on to each of the integration results.
        This function for the average density can then be re-arranged to calculate the central plasma density n_0 / ncore.
        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650
        Authors:
            Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre
            C. Ashe, CCFE, Culham Science Centre

        Parameters:
        - rhopedn (float): The normalised minor radius pedestal position.
        - nped (float): The pedestal density (/m3).
        - nsep (float): The separatrix density (/m3).
        - nav (float): The electron density (/m3).
        - alphan (float): The density peaking parameter.

        Returns:
        - ncore (float): The core density.

        """

        ncore = (
            1
            / (3 * rhopedn**2)
            * (
                3 * nav * (1 + alphan)
                + nsep * (1 + alphan) * (-2 + rhopedn + rhopedn**2)
                - nped * ((1 + alphan) * (1 + rhopedn) + (alphan - 2) * rhopedn**2)
            )
        )

        if ncore < 0.0:
            # Allows solver to continue and
            # warns the user to raise the lower bound on dene if the run did not converge
            error_handling.report_error(282)
            ncore = 1.0e-6
        return ncore

    def set_physics_variables(self) -> None:
        """Calculates and sets physics variables required for the profile."""

        if physics_variables.ipedestal == 0:
            physics_variables.ne0 = physics_variables.dene * (
                1.0 + physics_variables.alphan
            )
        elif physics_variables.ipedestal == 1:
            physics_variables.ne0 = self.ncore(
                physics_variables.rhopedn,
                physics_variables.neped,
                physics_variables.nesep,
                physics_variables.dene,
                physics_variables.alphan,
            )
        physics_variables.ni0 = (
            physics_variables.nd_ions_total
            / physics_variables.dene
            * physics_variables.ne0
        )


class TeProfile(Profile):
    """Electron temperature profile class. Contains a function to calculate the temperature profile and store the data."""

    def run(self) -> None:
        """Subroutine to initialise neprofile and execute calculations."""
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.set_physics_variables()
        self.calculate_profile_y(
            self.profile_x,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )
        self.integrate_profile_y()

    def calculate_profile_y(
        self,
        rho: np.array,
        rhopedt: float,
        t0: float,
        teped: float,
        tesep: float,
        alphat: float,
        tbeta: float,
    ) -> None:
        """
        Calculates the temperature at a normalised minor radius position rho for a pedestalised profile (teprofile).
        If ipedestal = 0 the original parabolic profile form is used instead.
        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650
        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre

        Args:
            rho (np.array): Normalised minor radius.
            rhopedt (float): Normalised minor radius pedestal position.
            t0 (float): Central temperature (keV).
            teped (float): Pedestal temperature (keV).
            tesep (float): Separatrix temperature (keV).
            alphat (float): Temperature peaking parameter.
            tbeta (float): Second temperature exponent.
        """
        if physics_variables.ipedestal == 0:
            self.profile_y = t0 * (1 - rho**2) ** alphat

            if t0 < teped:
                logger.info(
                    f"TPROFILE: temperature pedestal is higher than core temperature. {teped = }, {t0 = }"
                )
        else:
            rho_index = rho <= rhopedt
            self.profile_y[rho_index] = (
                teped + (t0 - teped) * (1 - (rho[rho_index] / rhopedt) ** tbeta) ** alphat
            )
            self.profile_y[~rho_index] = tesep + (teped - tesep) * (1 - rho[~rho_index]) / (
                1 - rhopedt
            )

    @staticmethod
    def tcore(
        rhopedt: float,
        teped: float,
        tesep: float,
        tav: float,
        alphat: float,
        tbeta: float,
    ) -> float:
        """
        This routine calculates the core temperature (keV)
        of a pedestalised profile. The solution comes from integrating and summing the two seprate temperature profiles for the core
        and pedestal region within their bounds. This has to be multiplied by the torus volume element before integration which leads
        to an added rho term in each part of the profile. When dividing by the volume of integration to get the average temperature
        the simplification leads to a factor of 2 having to be multiplied on to each of the integration results.
        This function for the average temperature can then be re-arranged to calculate the central plasma temeprature T_0 / tcore.
        References:
            Jean, J. (2011). HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies. Fusion Science and Technology, 59(2), 308-349. https://doi.org/10.13182/FST11-A11650
        Authors:
            Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre
            C. Ashe, CCFE, Culham Science Centre

        Args:
            rhopedt (float): Normalised minor radius pedestal position.
            teped (float): Pedestal temperature (keV).
            tesep (float): Separatrix temperature (keV).
            tav (float): Volume average temperature (keV).
            alphat (float): Temperature peaking parameter.
            tbeta (float): Second temperature exponent.

        Returns:
            float: Core temperature.
        """
        #  Calculate core temperature

        return teped + (
            (
                tbeta
                * (
                    3 * tav
                    + tesep * (-2.0 + rhopedt + rhopedt**2)
                    - teped * (1 + rhopedt + rhopedt**2)
                )
            )
            / (6 * rhopedt**2 * sp.special.beta(1 + alphat, 2 / tbeta))
        )

    def set_physics_variables(self) -> None:
        """
        Calculates and sets physics variables required for the temperature profile.

        Args:
            None

        Returns:
            None
        """
        if physics_variables.ipedestal == 0:
            physics_variables.te0 = physics_variables.te * (
                1.0 + physics_variables.alphat
            )
        elif physics_variables.ipedestal == 1:
            physics_variables.te0 = self.tcore(
                physics_variables.rhopedt,
                physics_variables.teped,
                physics_variables.tesep,
                physics_variables.te,
                physics_variables.alphat,
                physics_variables.tbeta,
            )

        physics_variables.ti0 = (
            physics_variables.ti / physics_variables.te * physics_variables.te0
        )
