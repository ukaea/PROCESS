import logging
import numpy as np
import scipy as sp
import process.profiles as profiles
import json
import matplotlib.pyplot as plt

from process.fortran import (
    constants,
    divertor_variables,
    physics_variables,
)

logger = logging.getLogger(__name__)


class PlasmaProfile:
    """
    Plasma profile class. Initiates the density and temperature profiles and
    handles the required physics variables.

    Attributes:
        profile_size (int): The size of the plasma profile.
        outfile (str): The output file path.
        neprofile (NProfile): An instance of the NProfile class.
        teprofile (TProfile): An instance of the TProfile class.

    Methods:
        run(): Subroutine to execute PlasmaProfile functions.
        parameterise_plasma(): Initializes the density and temperature profile averages and peak values.
        parabolic_paramterisation(): Parameterizes plasma profiles in the case where ipedestal=0.
        pedestal_parameterisation(): Instance temperature and density profiles then integrate them, setting physics variables ten and tin.
        calculate_profile_factors(): Calculate and set the central pressure (p0) using the ideal gas law and the pressure profile index (alphap).
        calculate_parabolic_profile_factors(): The gradient information for ipedestal = 0.
    """

    def __init__(self, custom_profile: dict = None) -> None:
        """
        Initialize the PlasmaProfile class.

        Args:
            profile_size (int): The size of the plasma profile.
            outfile (str): The output file path.
            neprofile (NProfile): An instance of the NProfile class.
            teprofile (TProfile): An instance of the TProfile class.
        """
        # Default profile_size = 501, but it's possible to experiment with this value.
        self.profile_size = 501
        self.outfile = constants.nout
        self.custom_profile = custom_profile

        if custom_profile:
            neprofile = custom_profile.get("neprofile", {})
            teprofile = custom_profile.get("teprofile", {})
            niprofile = custom_profile.get("niprofile", {})
            tiprofile = custom_profile.get("tiprofile", {})
            # Extract values, raise error if keys are not found
            try:
                neprofile_x = neprofile["profile_x"]  # Raises KeyError if not found
                neprofile_y = neprofile["profile_y"]
                teprofile_x = teprofile["profile_x"]
                teprofile_y = teprofile["profile_y"]
                niprofile_x = niprofile["profile_x"]  # Raises KeyError if not found
                niprofile_y = niprofile["profile_y"]
                tiprofile_x = tiprofile["profile_x"]
                tiprofile_y = tiprofile["profile_y"]
            except KeyError as e:
                raise KeyError(f"Missing required key in profile data: {e}")
            # Set correct profile size
            self.profile_size = len(neprofile_x)
        else:
            # If custom_profile is not provided, use default values
            neprofile_x = np.arange(self.profile_size)  # Default x range
            neprofile_y = np.zeros(
                self.profile_size
            )  # Default density profile as zeros
            teprofile_x = np.arange(self.profile_size)  # Default x range
            teprofile_y = np.zeros(
                self.profile_size
            )  # Default temperature profile as zeros
        # Initialize the density and temperature profiles
        self.neprofile = profiles.NProfile(self.profile_size, neprofile_x, neprofile_y)
        self.teprofile = profiles.TProfile(self.profile_size, teprofile_x, teprofile_y)

    def run(self) -> None:
        """
        Subroutine to execute PlasmaProfile functions.

        This method calls the parameterise_plasma() method to initialize the plasma profiles.

        Returns:
            None
        """
        self.parameterise_plasma()

    def parameterise_plasma(self) -> None:
        """
        This routine initializes the density and temperature
        profile averages and peak values, given the main
        parameters describing these profiles.

        Args:
            None

        Returns:
            None

        Authors:
            P J Knight, CCFE, Culham Science Centre

        References:
            T&M/PKNIGHT/LOGBOOK24, pp.4-7
        """

        #  Volume-averaged ion temperature
        #  (input value used directly if tratio=0.0)

        if physics_variables.tratio > 0.0e0:
            physics_variables.ti = physics_variables.tratio * physics_variables.te

        # Parabolic profile case
        if physics_variables.ipedestal == 0:
            if self.custom_profile is None:
                self.calculate_parabolic_profiles()
            else:
                self.calculate_custom_profiles()

            self.parabolic_paramterisation()
            self.calculate_profile_factors()
            self.calculate_parabolic_profile_factors()
        # Pedestal profile case
        else:
            if self.custom_profile is None:
                self.teprofile.run()
                self.neprofile.run()
            else:
                self.calculate_custom_profiles()

            self.pedestal_parameterisation()
            self.calculate_profile_factors()

    def calculate_custom_profiles(self) -> None:
        """
        Calculate the properties of custom input temperature and density profiles.

        This method normalizes the profile x-values, calculates the profile dx-values,
        and sets the physics variables for the density and temperature profiles.
        It also fits the custom profiles to the provided data and stores the original
        and fitted profiles in JSON format.

        Args:
            None

        Returns:
            None
        """
        # First fit the temp profile. The fit is better, get a better estimate of pedestal
        # position.
        self.teprofile.normalise_profile_x()
        self.teprofile.calculate_profile_dx()
        # Known values of t0 and tesep
        t0 = self.teprofile.profile_y[0]
        tesep = self.teprofile.profile_y[-1]
        # Initial guesses for rhopedt, teped, alphat, and tbeta
        initial_guess = [0.5, 2.8, 0.5, 1.5]

        def fit_temperature_profile(
            rho_data, temperature_data, t0, tesep, initial_guess
        ):
            # Fit the curve to the data
            popt, pcov = sp.optimize.curve_fit(
                lambda rho, rhopedt, teped, alphat, tbeta: self.teprofile.calculate_profile_y(
                    rho, rhopedt, t0, teped, tesep, alphat, tbeta
                ),
                rho_data,
                temperature_data,
                p0=initial_guess,
            )
            perr = np.sqrt(np.diag(pcov))
            # Extract the fitted parameters
            rhopedt_fitted, teped_fitted, alphat_fitted, tbeta_fitted = popt

            return (
                rhopedt_fitted,
                teped_fitted,
                alphat_fitted,
                tbeta_fitted,
                perr,
            )

        # Initial guesses for rhopedt, teped, alphat, and tbeta
        initial_guess = [0.96, 0.5, 0.5, 1.5]
        # Fit the profile and get the fitted parameters
        (
            rhopedt_fitted,
            teped_fitted,
            alphat_fitted,
            tbeta_fitted,
            teperr,
        ) = fit_temperature_profile(
            self.teprofile.profile_x,
            self.teprofile.profile_y,
            t0,
            tesep,
            initial_guess,
        )
        # Create the fitted profile using the fitted parameters
        fitted_profile = self.teprofile.calculate_profile_y(
            self.teprofile.profile_x,
            rhopedt_fitted,
            t0,
            teped_fitted,
            tesep,
            alphat_fitted,
            tbeta_fitted,
        )
        # Store the original and fitted profiles in JSON format
        profiles = {
            "original_profile": self.teprofile.profile_y.tolist(),
            "fitted_profile": fitted_profile.tolist(),
        }

        with open("temperature_profiles.json", "w") as f:
            json.dump(profiles, f)

        self.teprofile.set_physics_variables()
        self.teprofile.integrate_profile_y()

        self.neprofile.normalise_profile_x()
        self.neprofile.calculate_profile_dx()
        self.neprofile.set_physics_variables()
        # Initial guess for the parameters (nped, nsep, nav, alphan)
        nsep = self.neprofile.profile_y[-1]
        n0 = self.neprofile.profile_y[0]
        # Initial guesses for nped, and alphan
        initial_guess = [0.75e20, 0.24]
        self.neprofile.integrate_profile_y()

        # Fit the density profile.
        def fit_density_profile(
            rho_data, density_data, n0, nsep, rhopedn, initial_guess
        ):
            # Fit the curve to the data
            popt, pcov = sp.optimize.curve_fit(
                lambda rho, nped, alphan, rhopedn=rhopedn: self.neprofile.calculate_profile_y(
                    rho, rhopedn, n0, nped, nsep, alphan
                ),
                rho_data,
                density_data,
                p0=initial_guess,
            )
            # Calculate the standard deviation errors
            perr = np.sqrt(np.diag(pcov))
            # Extract the fitted parameters
            nped_fitted, alphan_fitted = popt

            return nped_fitted, alphan_fitted, perr

        # Fit the profile and get the fitted parameters
        nped_fitted, alphan_fitted, perr = fit_density_profile(
            rho_data=self.neprofile.profile_x,
            density_data=self.neprofile.profile_y,
            n0=n0,
            nsep=nsep,
            rhopedn=rhopedt_fitted,
            initial_guess=initial_guess,
        )
        # Create the fitted profile using the fitted parameters
        fitted_profile = self.neprofile.calculate_profile_y(
            self.neprofile.profile_x,
            rhopedt_fitted,
            n0,
            nped_fitted,
            nsep,
            alphan_fitted,
        )

        # Store the original and fitted profiles in JSON format
        profiles = {
            "original_profile": self.neprofile.profile_y.tolist(),
            "fitted_profile": fitted_profile.tolist(),
        }

        with open("profiles.json", "w") as f:
            json.dump(profiles, f)
        physics_variables.alphan = alphan_fitted
        physics_variables.alphat = alphat_fitted
        physics_variables.tbeta = tbeta_fitted

    def calculate_parabolic_profiles(self) -> None:
        """Reset the pedestal values and calculate the profiles."""
        # Reset pedestal values to agree with original parabolic profiles
        if (
            physics_variables.rhopedt != 1.0
            or physics_variables.rhopedn != 1.0
            or physics_variables.teped != 0.0
            or physics_variables.tesep != 0.0
            or physics_variables.neped != 0.0
            or physics_variables.nesep != 0.0
            or physics_variables.tbeta != 2.0
        ):
            logger.warning(
                "Parabolic plasma profiles is used for an L-Mode plasma, "
                "but the physics variables do not describe an L-Mode plasma. "
                "'rhopedt', 'rhopedn', 'teped', 'tesep', 'neped', 'nesep', "
                "and 'tbeta' have all been reset to L-Mode appropriate values"
            )

            physics_variables.rhopedt = 1.0e0
            physics_variables.rhopedn = 1.0e0
            physics_variables.teped = 0.0e0
            physics_variables.tesep = 0.0e0
            physics_variables.neped = 0.0e0
            physics_variables.nesep = 0.0e0
            physics_variables.tbeta = 2.0e0

        # Re-caluclate core and profile values
        self.teprofile.run()
        self.neprofile.run()

    def parabolic_paramterisation(self) -> None:
        """
        Parameterise plasma profiles in the case where ipedestal == 0.

        This routine calculates the parameterization of plasma profiles in the case where ipedestal=0.
        It sets the necessary physics variables for the parabolic profile case.

        Args:
            None

        Returns:
            None
        """
        #  Profile factor; ratio of density-weighted to volume-averaged
        #  temperature

        physics_variables.pcoef = (
            (1.0e0 + physics_variables.alphan)
            * (1.0e0 + physics_variables.alphat)
            / (1.0e0 + physics_variables.alphan + physics_variables.alphat)
        )

        # Line averaged electron density (IPDG89)
        # Taken by integrating the parabolic profile over rho in the bounds of 0 and 1 and dividng by the width of the integration bounds

        physics_variables.dnla = (
            physics_variables.dene
            * (1.0 + physics_variables.alphan)
            * (sp.special.gamma(0.5) / 2.0)
            * sp.special.gamma(physics_variables.alphan + 1.0)
            / sp.special.gamma(physics_variables.alphan + 1.5)
        )

        #  Density-weighted temperatures

        physics_variables.ten = physics_variables.te * physics_variables.pcoef
        physics_variables.tin = physics_variables.ti * physics_variables.pcoef

        #  Central values for temperature (keV) and density (m**-3)

        physics_variables.te0 = physics_variables.te * (1.0 + physics_variables.alphat)
        physics_variables.ti0 = physics_variables.ti * (1.0 + physics_variables.alphat)

        physics_variables.ne0 = physics_variables.dene * (
            1.0 + physics_variables.alphan
        )
        physics_variables.ni0 = physics_variables.dnitot * (
            1.0 + physics_variables.alphan
        )

    def pedestal_parameterisation(self) -> None:
        """
        Instance temperature and density profiles then integrate them, setting physics variables ten and tin.

        This routine instances temperature and density profiles and integrates them to calculate the values of the physics variables `ten` and `tin`.

        Args:
            None

        Returns:
            None
        """
        #  Perform integrations to calculate ratio of density-weighted
        #  to volume-averaged temperature, etc.
        #  Density-weighted temperature = integral(n.T dV) / integral(n dV)
        #  which is approximately equal to the ratio
        #  integral(rho.n(rho).T(rho) drho) / integral(rho.n(rho) drho)

        drho = self.neprofile.profile_dx
        dens = self.neprofile.profile_y
        temp = self.teprofile.profile_y
        rho = self.neprofile.profile_x

        arg1 = rho * dens * temp
        arg2 = rho * dens

        integ1 = sp.integrate.simpson(arg1, x=rho, dx=drho)
        integ2 = sp.integrate.simpson(arg2, x=rho, dx=drho)

        physics_variables.integ1 = integ1
        physics_variables.integ2 = integ2
        #  Density-weighted temperatures
        physics_variables.ten = integ1 / integ2
        physics_variables.tin = (
            physics_variables.ti / physics_variables.te * physics_variables.ten
        )
        #  Profile factor; ratio of density-weighted to volume-averaged
        #  temperature

        physics_variables.pcoef = physics_variables.ten / physics_variables.te

        #  Line-averaged electron density
        #  = integral(n(rho).drho)

        physics_variables.dnla = self.neprofile.profile_integ

        #  Scrape-off density / volume averaged density
        #  (Input value is used if ipedestal = 0)

        divertor_variables.prn1 = max(
            0.01e0, physics_variables.nesep / physics_variables.dene
        )  # Preventing division by zero later

    def calculate_profile_factors(self) -> None:
        """
        Calculate and set the central pressure (p0) using the ideal gas law and the pressure profile index (alphap).

        This method calculates the central pressure (p0) using the ideal gas law and the pressure profile index (alphap).
        It sets the value of the physics variable `p0`.

        Args:
            None

        Returns:
            None
        """

        #  Central pressure (Pa), from ideal gas law : p = nkT

        physics_variables.p0 = (
            (
                physics_variables.ne0 * physics_variables.te0
                + physics_variables.ni0 * physics_variables.ti0
            )
            * 1.0e3
            * constants.electron_charge
        )

        #  Pressure profile index (N.B. no pedestal effects included here)
        #  N.B. p0 is NOT equal to <p> * (1 + alphap), but p(rho) = n(rho)*T(rho)
        #  and <p> = <n>.T_n where <...> denotes volume-averages and T_n is the
        #  density-weighted temperature

        physics_variables.alphap = physics_variables.alphan + physics_variables.alphat

    @staticmethod
    def calculate_parabolic_profile_factors() -> None:
        """
        Calculate the gradient information for ipedestal = 0.

        This function calculates the gradient information for the plasma profiles at the pedestal region
        when the value of ipedestal is 0. It is used by the stellarator routines.

        The function uses analytical parametric formulas to calculate the gradient information.
        The maximum normalized radius (rho_max) is obtained by equating the second derivative to zero.

        Args:
            None

        Returns:
            None
        """
        if physics_variables.ipedestal == 0:
            if physics_variables.alphat > 1.0:
                # Rho (normalized radius), where temperature derivative is largest
                rho_te_max = 1.0 / np.sqrt(-1.0 + 2.0 * physics_variables.alphat)
                dtdrho_max = (
                    -(2.0**physics_variables.alphat)
                    * (-1.0 + physics_variables.alphat)
                    ** (-1.0 + physics_variables.alphat)
                    * physics_variables.alphat
                    * (-1.0 + 2.0 * physics_variables.alphat)
                    ** (0.5e0 - physics_variables.alphat)
                    * physics_variables.te0
                )
                te_max = (
                    physics_variables.te0
                    * (1 - rho_te_max**2) ** physics_variables.alphat
                )

            elif physics_variables.alphat <= 1.0 and physics_variables.alphat > 0.0:
                # This makes the profiles very 'boxy'
                # The gradient diverges here at the edge so define some 'wrong' value of 0.9
                # to approximate the gradient
                rho_te_max = 0.9
                dtdrho_max = (
                    -2.0
                    * physics_variables.alphat
                    * rho_te_max
                    * (1 - rho_te_max**2) ** (-1.0 + physics_variables.alphat)
                    * physics_variables.te0
                )
                te_max = (
                    physics_variables.te0
                    * (1 - rho_te_max**2) ** physics_variables.alphat
                )
            else:
                raise ValueError(f"alphat is negative: { physics_variables.alphat}")

            # Same for density
            if physics_variables.alphan > 1.0:
                rho_ne_max = 1.0 / np.sqrt(-1.0 + 2.0 * physics_variables.alphan)
                dndrho_max = (
                    -(2.0**physics_variables.alphan)
                    * (-1.0 + physics_variables.alphan)
                    ** (-1.0 + physics_variables.alphan)
                    * physics_variables.alphan
                    * (-1.0 + 2.0 * physics_variables.alphan)
                    ** (0.5 - physics_variables.alphan)
                    * physics_variables.ne0
                )
                ne_max = (
                    physics_variables.ne0
                    * (1e0 - rho_ne_max**2) ** physics_variables.alphan
                )
            elif physics_variables.alphan <= 1.0 and physics_variables.alphan > 0.0:
                # This makes the profiles very 'boxy'
                # The gradient diverges here at the edge so define some 'wrong' value of 0.9
                # to approximate the gradient
                rho_ne_max = 0.9
                dndrho_max = (
                    -2.0
                    * physics_variables.alphan
                    * rho_ne_max
                    * (1 - rho_ne_max**2) ** (-1.0 + physics_variables.alphan)
                    * physics_variables.ne0
                )
                ne_max = (
                    physics_variables.ne0
                    * (1 - rho_ne_max**2) ** physics_variables.alphan
                )
            else:
                raise ValueError(f"alphan is negative: { physics_variables.alphan}")

            # set normalized gradient length
            # te at rho_te_max
            physics_variables.gradient_length_te = (
                -dtdrho_max * physics_variables.rminor * rho_te_max / te_max
            )
            # same for density:
            physics_variables.gradient_length_ne = (
                -dndrho_max * physics_variables.rminor * rho_ne_max / ne_max
            )
