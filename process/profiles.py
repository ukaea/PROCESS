import numpy as np
import logging
from scipy import integrate
from abc import ABC, abstractmethod

from process.fortran import maths_library, physics_variables, error_handling

logger = logging.getLogger(__name__)
# Logging handler for console output
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.ERROR)
logger.addHandler(s_handler)


class Profile(ABC):
    """Abstract base class used to create and hold profiles (temperature, desnity,)

    :param ABC: Abstract Base Class
    :type ABC: ABC, optional
    """

    def __init__(self, profile_size):
        """
        :param profile_size: Number of points in the profile
        :type profile_size: int
        """
        self.profile_size = profile_size
        self.profile_x = np.arange(self.profile_size)
        self.profile_y = np.zeros(self.profile_size)
        self.profile_integ = 0
        self.profile_dx = 0

    def normalise_profile_x(self):
        """Normalise the profile x-dimension."""
        self.profile_x = self.profile_x / max(self.profile_x)

    def calculate_profile_dx(self):
        """Calculates the differential between points in the profile."""
        self.profile_dx = max(self.profile_x) / (self.profile_size - 1)

    @abstractmethod
    def calculate_profile_y(self):
        """Use a profile function to act on self.profile_x to calculate and set the
        values of self.profile_y.
        """
        pass

    def integrate_profile_y(self):
        """
        Integrate profile_y values using scipy.integrate.simpson() function.
        """
        self.profile_integ = integrate.simpson(
            self.profile_y, x=self.profile_x, dx=self.profile_dx
        )


class NProfile(Profile):
    """Density profile class.
    Contains a function to calculate the electron density profile and
    store the data.
    """

    def run(self):
        """_summary_
        Subroutine which calls functions and stores nprofile data.
        """
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

    def calculate_profile_y(self, rho, rhopedn, n0, nped, nsep, alphan):
        """This routine calculates the density at each normalised minor radius position
        rho for a ELIOS-type density pedestal profile (nprofile).
        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre
        References:
            J.Johner, Fusion Science and Technology 59 (2011), pp 308-349

        :param rho: normalised minor radius vector
        :type rho: float
        :param rhopedn: normalised minor radius pedestal position
        :type rhopedn: float
        :param n0: central density (/m3)
        :type n0: float
        :param nped: oedestal desnity (/m3)
        :type nped: float
        :param nsep: separatrix density (/m3)
        :type nsep: float
        :param alphan: density peaking parameter
        :type alphan: float
        """

        if physics_variables.ipedestal == 0:
            self.profile_y = n0 * (1 - rho**2) ** alphan

        #  Error trap; shouldn't happen unless volume-averaged density has
        #  been allowed to drop below nped. This may happen during a HYBRD case,
        #  but should have been prevented for optimisation runs.

        #  Input checks

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
        """This routine calculates the core denesity of a pedestalised profile.

        :param rhopedn: normalised minor radius pedestal position
        :type rhopedn: numpy.array
        :param nped: pedestal density (/m3)
        :type nped: float
        :param nsep: separatrix desnity (/m3)
        :type nsep: float
        :param nav: electron density (/m3)
        :type nav: float
        :param alphan: density peaking parameter
        :type alphan: float
        :return: Core density
        :type: float
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

    def set_physics_variables(self):
        """_summary_
        Calculates and sets physics vcariables required for the profile.
        Returns:
            physics_variables.ne0 (float) : central electron density (/m3)
            physics_variables.ni0 (float) : central ion density (/m3)
        """
        physics_variables.ne0 = self.ncore(
            physics_variables.rhopedn,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.dene,
            physics_variables.alphan,
        )
        physics_variables.ni0 = (
            physics_variables.dnitot / physics_variables.dene * physics_variables.ne0
        )


class TProfile(Profile):
    """
    Temperature profile class.
    Contains a function to calculate the temperature profile and store the data.
    """

    def run(self):
        """Subroutine to initialise nprofile and execute calculations."""
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

    def calculate_profile_y(self, rho, rhopedt, t0, teped, tesep, alphat, tbeta):
        """This routine calculates the temperature at a normalised minor
        radius position rho for a pedestalised profile (tprofile).
        If ipedestal = 0 the original parabolic profile form is used instead.
        References:
            J.Johner, Fusion Science and Technology 59 (2011), pp 308-349
        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre

        :param rho: normalised minor radius
        :type rho: numpy.array
        :param rhopedt:  normalised minor radius pedestal position
        :type rhopedt: numpy.array
        :param t0: central temperature (keV)
        :type t0: float
        :param teped: pedestal temperature (keV)
        :type teped: float
        :param tesep: separatrix temperature (keV)
        :type tesep: float
        :param alphat:temperature peaking parameter
        :type alphat: float
        :param tbeta: second temperature exponent
        :type tbeta: float
        """
        if physics_variables.ipedestal == 0:
            self.profile_y = t0 * (1 - rho**2) ** alphat

        #  Error trap; shouldn't happen unless volume-averaged temperature has
        #  been allowed to drop below teped. This may happen during a HYBRD case,
        #  but should have been prevented for optimisation runs.

        if t0 < teped:
            logger.info(
                f"TPROFILE: temperature pedestal is higher than core temperature. {teped = }, {t0 = }"
            )
        rho_index = rho <= rhopedt
        self.profile_y[rho_index] = (
            teped + (t0 - teped) * (1 - (rho[rho_index] / rhopedt) ** tbeta) ** alphat
        )
        self.profile_y[~rho_index] = tesep + (teped - tesep) * (1 - rho[~rho_index]) / (
            1 - rhopedt
        )

    @staticmethod
    def tcore(rhopedt, teped, tesep, tav, alphat, tbeta):
        """This routine calculates the core temperature (keV)
        of a pedestalised profile.
        References:
            J.Johner, Fusion Science and Technology 59 (2011), pp 308-349
        Authors:
            Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre

        :param rhopedt: normalised minor radius pedestal position
        :type rhopedt: numpy.array
        :param teped: pedestal temperature (keV)
        :type teped: float
        :param tesep: separatrix temperature (keV)
        :type tesep: float
        :param tav: volume average temperature (keV)
        :type tav: float
        :param alphat: temperature peaking parameter
        :type alphat: float
        :param tbeta: second temperature exponent
        :type tbeta: float
        :return: core temperature
        :rtype: numpy.array
        """
        #  For integer values of alphat, the limit of
        #  gamfun(-alphat)*sin(pi*alphat) needs to be calculated directly

        gamfac = (
            maths_library.gamfun(1 + alphat + 2 / tbeta)
            / maths_library.gamfun((2 + tbeta) / tbeta)
            / rhopedt**2
        )
        if abs(alphat - np.around(alphat)) <= 1e-7:
            gamfac = -gamfac / maths_library.gamfun(1 + alphat)
        else:
            gamfac = (
                gamfac * maths_library.gamfun(-alphat) * np.sin(np.pi * alphat) / np.pi
            )

        #  Calculate core temperature

        return teped + gamfac * (
            teped * rhopedt**2
            - tav
            + (1 - rhopedt) / 3 * ((1 + 2 * rhopedt) * teped + (2 + rhopedt) * tesep)
        )

    def set_physics_variables(self):
        """Calculates and sets physics vcariables required for the temperature profile."""
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
