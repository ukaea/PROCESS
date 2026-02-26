import logging
from enum import IntEnum

import numpy as np
import scipy
import scipy.integrate as integrate

from process import constants
from process import process_output as po
from process.data_structure import (
    current_drive_variables,
    physics_variables,
)
from process.exceptions import ProcessValueError
from process.models.physics.plasma_profiles import PlasmaProfile

logger = logging.getLogger(__name__)


class BootstrapCurrentFractionModel(IntEnum):
    """Bootstrap plasma current fraction (f_BS) model types"""

    USER_INPUT = 0
    ITER_89 = 1
    NEVINS = 2
    WILSON = 3
    SAUTER = 4
    SAKAI = 5
    ARIES = 6
    ANDRADE = 7
    HOANG = 8
    WONG = 9
    GI_1 = 10
    GI_2 = 11
    SUGIYAMA_L_MODE = 12
    SUGIYAMA_H_MODE = 13


class PlasmaBootstrapCurrent:
    """Class to hold plasma bootstrap current for plasma processing."""

    def __init__(self, plasma_profile: PlasmaProfile) -> None:
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile
        self.sauter_bootstrap = SauterBootstrapCurrent()

    def run(self) -> None:
        # Calculate bootstrap current fraction using various models
        current_drive_variables.f_c_plasma_bootstrap_iter89 = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_iter89(
                physics_variables.aspect,
                physics_variables.beta_total_vol_avg,
                physics_variables.b_plasma_total,
                physics_variables.plasma_current,
                physics_variables.q95,
                physics_variables.q0,
                physics_variables.rmajor,
                physics_variables.vol_plasma,
            )
        )

        current_drive_variables.f_c_plasma_bootstrap_nevins = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_nevins(
                physics_variables.alphan,
                physics_variables.alphat,
                physics_variables.beta_toroidal_vol_avg,
                physics_variables.b_plasma_toroidal_on_axis,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.plasma_current,
                physics_variables.q95,
                physics_variables.q0,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.temp_plasma_electron_vol_avg_kev,
                physics_variables.n_charge_plasma_effective_vol_avg,
            )
        )

        # Wilson scaling uses thermal poloidal beta, not total
        current_drive_variables.f_c_plasma_bootstrap_wilson = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wilson(
                physics_variables.alphaj,
                physics_variables.alphap,
                physics_variables.alphat,
                physics_variables.beta_thermal_poloidal_vol_avg,
                physics_variables.q0,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.rminor,
            )
        )

        (
            current_drive_variables.f_c_plasma_bootstrap_sauter,
            physics_variables.j_plasma_bootstrap_sauter_profile,
        ) = self.bootstrap_fraction_sauter(self.plasma_profile)
        current_drive_variables.f_c_plasma_bootstrap_sauter *= (
            current_drive_variables.cboot
        )

        current_drive_variables.f_c_plasma_bootstrap_sakai = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sakai(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                eps=physics_variables.eps,
                ind_plasma_internal_norm=physics_variables.ind_plasma_internal_norm,
            )
        )

        current_drive_variables.f_c_plasma_bootstrap_aries = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_aries(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                ind_plasma_internal_norm=physics_variables.ind_plasma_internal_norm,
                core_density=physics_variables.nd_plasma_electron_on_axis,
                average_density=physics_variables.nd_plasma_electrons_vol_avg,
                inverse_aspect=physics_variables.eps,
            )
        )

        current_drive_variables.f_c_plasma_bootstrap_andrade = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_andrade(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                core_pressure=physics_variables.pres_plasma_thermal_on_axis,
                average_pressure=physics_variables.pres_plasma_thermal_vol_avg,
                inverse_aspect=physics_variables.eps,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_hoang = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_hoang(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                pressure_index=physics_variables.alphap,
                current_index=physics_variables.alphaj,
                inverse_aspect=physics_variables.eps,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_wong = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wong(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                density_index=physics_variables.alphan,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                elongation=physics_variables.kappa,
            )
        )
        current_drive_variables.bscf_gi_i = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_gi_I(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                pressure_index=physics_variables.alphap,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                effective_charge=physics_variables.n_charge_plasma_effective_vol_avg,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
            )
        )

        current_drive_variables.bscf_gi_ii = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_gi_II(
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                pressure_index=physics_variables.alphap,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                effective_charge=physics_variables.n_charge_plasma_effective_vol_avg,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_sugiyama_l = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sugiyama_l_mode(
                eps=physics_variables.eps,
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_sugiyama_h = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sugiyama_h_mode(
                eps=physics_variables.eps,
                beta_poloidal=physics_variables.beta_poloidal_vol_avg,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                tbeta=physics_variables.tbeta,
                zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
                radius_plasma_pedestal_density_norm=physics_variables.radius_plasma_pedestal_density_norm,
                nd_plasma_pedestal_electron=physics_variables.nd_plasma_pedestal_electron,
                n_greenwald=physics_variables.nd_plasma_electron_max_array[6],
                temp_plasma_pedestal_kev=physics_variables.temp_plasma_pedestal_kev,
            )
        )

        # Calculate beta_norm_max based on i_beta_norm_max
        try:
            model = BootstrapCurrentFractionModel(
                int(physics_variables.i_bootstrap_current)
            )
            current_drive_variables.f_c_plasma_bootstrap = (
                self.get_bootstrap_current_fraction_value(model)
            )
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_bootstrap_current",
                i_bootstrap_current=physics_variables.i_bootstrap_current,
            ) from None

    def get_bootstrap_current_fraction_value(
        self, model: BootstrapCurrentFractionModel
    ) -> float:
        """Get the plasma current bootstrap fraction (f_BS) for the specified model."""
        model_map = {
            BootstrapCurrentFractionModel.USER_INPUT: current_drive_variables.f_c_plasma_bootstrap,
            BootstrapCurrentFractionModel.ITER_89: current_drive_variables.f_c_plasma_bootstrap_iter89,
            BootstrapCurrentFractionModel.NEVINS: current_drive_variables.f_c_plasma_bootstrap_nevins,
            BootstrapCurrentFractionModel.WILSON: current_drive_variables.f_c_plasma_bootstrap_wilson,
            BootstrapCurrentFractionModel.SAUTER: current_drive_variables.f_c_plasma_bootstrap_sauter,
            BootstrapCurrentFractionModel.SAKAI: current_drive_variables.f_c_plasma_bootstrap_sakai,
            BootstrapCurrentFractionModel.ARIES: current_drive_variables.f_c_plasma_bootstrap_aries,
            BootstrapCurrentFractionModel.ANDRADE: current_drive_variables.f_c_plasma_bootstrap_andrade,
            BootstrapCurrentFractionModel.HOANG: current_drive_variables.f_c_plasma_bootstrap_hoang,
            BootstrapCurrentFractionModel.WONG: current_drive_variables.f_c_plasma_bootstrap_wong,
            BootstrapCurrentFractionModel.GI_1: current_drive_variables.bscf_gi_i,
            BootstrapCurrentFractionModel.GI_2: current_drive_variables.bscf_gi_ii,
            BootstrapCurrentFractionModel.SUGIYAMA_L_MODE: current_drive_variables.f_c_plasma_bootstrap_sugiyama_l,
            BootstrapCurrentFractionModel.SUGIYAMA_H_MODE: current_drive_variables.f_c_plasma_bootstrap_sugiyama_h,
        }
        return model_map[model]

    @staticmethod
    def bootstrap_fraction_iter89(
        aspect: float,
        beta: float,
        b_plasma_toroidal_on_axis: float,
        plasma_current: float,
        q95: float,
        q0: float,
        rmajor: float,
        vol_plasma: float,
    ) -> float:
        """
        Calculate the bootstrap-driven fraction of the plasma current.

        Args:
            aspect (float): Plasma aspect ratio.
            beta (float): Plasma total beta.
            b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
            plasma_current (float): Plasma current (A).
            q95 (float): Safety factor at 95% surface.
            q0 (float): Central safety factor.
            rmajor (float): Plasma major radius (m).
            vol_plasma (float): Plasma volume (m3).

        Returns:
            float: The bootstrap-driven fraction of the plasma current.

        This function performs the original ITER calculation of the plasma current bootstrap fraction.

        Reference:
        - ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        - ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """

        # Calculate the bootstrap current coefficient
        c_bs = 1.32 - 0.235 * (q95 / q0) + 0.0185 * (q95 / q0) ** 2

        # Calculate the average minor radius
        average_a = np.sqrt(vol_plasma / (2 * np.pi**2 * rmajor))

        b_pa = (plasma_current / 1e6) / (5 * average_a)

        # Calculate the poloidal beta for bootstrap current
        betapbs = beta * (b_plasma_toroidal_on_axis / b_pa) ** 2

        # Ensure betapbs is positive
        if betapbs <= 0.0:
            return 0.0

        # Calculate and return the bootstrap current fraction
        return c_bs * (betapbs / np.sqrt(aspect)) ** 1.3

    @staticmethod
    def bootstrap_fraction_wilson(
        alphaj: float,
        alphap: float,
        alphat: float,
        betpth: float,
        q0: float,
        q95: float,
        rmajor: float,
        rminor: float,
    ) -> float:
        """
        Bootstrap current fraction from Wilson et al scaling

        Args:
            alphaj (float): Current profile index.
            alphap (float): Pressure profile index.
            alphat (float): Temperature profile index.
            betpth (float): Thermal component of poloidal beta.
            q0 (float): Safety factor on axis.
            q95 (float): Edge safety factor.
            rmajor (float): Major radius (m).
            rminor (float): Minor radius (m).

        Returns:
            float: The bootstrap current fraction.

        This function calculates the bootstrap current fraction using the numerically fitted algorithm written by Howard Wilson.

        Reference:
            - AEA FUS 172: Physics Assessment for the European Reactor Study, 1989
            - H. R. Wilson, Nuclear Fusion 32 (1992) 257
        """

        term1 = np.log(0.5)
        term2 = np.log(q0 / q95)

        # Re-arranging of parabolic profile to be equal to (r/a)^2 where the profile value is half of the core value

        termp = 1.0 - 0.5 ** (1.0 / alphap)
        termt = 1.0 - 0.5 ** (1.0 / alphat)
        termj = 1.0 - 0.5 ** (1.0 / alphaj)

        # Assuming a parabolic safety factor profile of the form q = q0 + (q95 - q0) * (r/a)^2
        # Substitute (r/a)^2 term from temperature,pressure and current profiles into q profile when values is 50% of core value
        # Take natural log of q profile over q95 and q0 to get the profile index

        alfpnw = term1 / np.log(np.log((q0 + (q95 - q0) * termp) / q95) / term2)
        alftnw = term1 / np.log(np.log((q0 + (q95 - q0) * termt) / q95) / term2)
        aj = term1 / np.log(np.log((q0 + (q95 - q0) * termj) / q95) / term2)

        # Crude check for NaN errors or other illegal values.
        if np.isnan(aj) or np.isnan(alfpnw) or np.isnan(alftnw) or aj < 0:
            raise ProcessValueError(
                "Illegal profile value found", aj=aj, alfpnw=alfpnw, alftnw=alftnw
            )

        # Ratio of ionic charge to electron charge

        z = 1.0

        # Inverse aspect ratio: r2 = maximum plasma radius, r1 = minimum
        # This is the definition used in the paper
        r2 = rmajor + rminor
        r1 = rmajor - rminor
        eps1 = (r2 - r1) / (r2 + r1)

        # Coefficients fitted using least squares techniques

        # Square root of current profile index term
        saj = np.sqrt(aj)

        a = np.array([
            1.41 * (1.0 - 0.28 * saj) * (1.0 + 0.12 / z),
            0.36 * (1.0 - 0.59 * saj) * (1.0 + 0.8 / z),
            -0.27 * (1.0 - 0.47 * saj) * (1.0 + 3.0 / z),
            0.0053 * (1.0 + 5.0 / z),
            -0.93 * (1.0 - 0.34 * saj) * (1.0 + 0.15 / z),
            -0.26 * (1.0 - 0.57 * saj) * (1.0 - 0.27 * z),
            0.064 * (1.0 - 0.6 * aj + 0.15 * aj * aj) * (1.0 + 7.6 / z),
            -0.0011 * (1.0 + 9.0 / z),
            -0.33 * (1.0 - aj + 0.33 * aj * aj),
            -0.26 * (1.0 - 0.87 / saj - 0.16 * aj),
            -0.14 * (1.0 - 1.14 / saj - 0.45 * saj),
            -0.0069,
        ])

        seps1 = np.sqrt(eps1)

        b = np.array([
            1.0,
            alfpnw,
            alftnw,
            alfpnw * alftnw,
            seps1,
            alfpnw * seps1,
            alftnw * seps1,
            alfpnw * alftnw * seps1,
            eps1,
            alfpnw * eps1,
            alftnw * eps1,
            alfpnw * alftnw * eps1,
        ])

        # Empirical bootstrap current fraction
        return seps1 * betpth * (a * b).sum()

    @staticmethod
    def _nevins_integral(
        y: float,
        nd_plasma_electrons_vol_avg: float,
        te: float,
        b_plasma_toroidal_on_axis: float,
        rminor: float,
        rmajor: float,
        zeff: float,
        alphat: float,
        alphan: float,
        q0: float,
        q95: float,
        beta_toroidal: float,
    ) -> float:
        """Integrand function for Nevins et al bootstrap current scaling.

        This function calculates the integrand function for the Nevins et al bootstrap current scaling.

        Parameters
        ----------
        y :
            abscissa of integration, normalized minor radius
        nd_plasma_electrons_vol_avg :
            volume averaged electron density (/m^3)
        te :
            volume averaged electron temperature (keV)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        rminor :
            plasma minor radius (m)
        rmajor :
            plasma major radius (m)
        zeff :
            plasma effective charge
        alphat :
            temperature profile index
        alphan :
            density profile index
        q0 :
            normalized safety factor at the magnetic axis
        q95 :
            normalized safety factor at 95% of the plasma radius
        beta_toroidal :
            Toroidal plasma beta


        Returns
        -------
        type
            - float, the integrand value


        Reference: See appendix of:
            Keii Gi, Makoto Nakamura, Kenji Tobita, Yasushi Ono,
            Bootstrap current fraction scaling for a tokamak reactor design study,
            Fusion Engineering and Design, Volume 89, Issue 11, 2014, Pages 2709-2715,
            ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2014.07.009.

            Nevins, W. M. "Summary report: ITER specialists' meeting on heating and current drive."
            ITER-TN-PH-8-4, June 1988. 1988.

        """

        # Compute average electron beta
        betae = (
            nd_plasma_electrons_vol_avg
            * te
            * 1.0e3
            * constants.ELECTRON_CHARGE
            / (b_plasma_toroidal_on_axis**2 / (2.0 * constants.RMU0))
        )

        nabla = rminor * np.sqrt(y) / rmajor
        x = (1.46 * np.sqrt(nabla) + 2.4 * nabla) / (1.0 - nabla) ** 1.5
        z = zeff
        d = (
            1.414 * z
            + z**2
            + x * (0.754 + 2.657 * z + (2.0 * z**2))
            + (x**2 * (0.348 + 1.243 * z + z**2))
        )
        a1 = (alphan + alphat) * (1.0 - y) ** (alphan + alphat - 1.0)
        a2 = alphat * (1.0 - y) ** (alphan + alphat - 1.0)
        al1 = (x / d) * (0.754 + 2.21 * z + z**2 + x * (0.348 + 1.243 * z + z**2))
        al2 = -x * ((0.884 + 2.074 * z) / d)
        alphai = -1.172 / (1.0 + 0.462 * x)

        # q-profile
        q = q0 + (q95 - q0) * ((y + y**2 + y**3) / (3.0))

        pratio = (beta_toroidal - betae) / betae

        return (q / q95) * (al1 * (a1 + (pratio * (a1 + alphai * a2))) + al2 * a2)

    def bootstrap_fraction_sauter(
        self, plasma_profile: PlasmaProfile
    ) -> tuple[float, np.ndarray]:
        """Get the bootstrap current fraction using Sauter et al scaling.

        Args:
            plasma_profile: The plasma profile object.

        Returns:
            tuple: (bootstrap fraction, bootstrap current density profile)
        """
        return self.sauter_bootstrap.bootstrap_fraction_sauter(plasma_profile)

    def bootstrap_fraction_nevins(
        self,
        alphan: float,
        alphat: float,
        beta_toroidal: float,
        b_plasma_toroidal_on_axis: float,
        nd_plasma_electrons_vol_avg: float,
        plasma_current: float,
        q95: float,
        q0: float,
        rmajor: float,
        rminor: float,
        te: float,
        zeff: float,
    ) -> float:
        """
        Calculate the bootstrap current fraction from Nevins et al scaling.

        Args:
            alphan (float): Density profile index.
            alphat (float): Temperature profile index.
            beta_toroidal (float): Toroidal plasma beta.
            b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
            nd_plasma_electrons_vol_avg (float): Electron density (/m3).
            plasma_current (float): Plasma current (A).
            q0 (float): Central safety factor.
            q95 (float): Safety factor at 95% surface.
            rmajor (float): Plasma major radius (m).
            rminor (float): Plasma minor radius (m).
            te (float): Volume averaged plasma temperature (keV).
            zeff (float): Plasma effective charge.

        Returns:
            float: The bootstrap current fraction.

        This function calculates the bootstrap current fraction using the Nevins et al method.

        Reference: See appendix of:
        Keii Gi, Makoto Nakamura, Kenji Tobita, Yasushi Ono,
        Bootstrap current fraction scaling for a tokamak reactor design study,
        Fusion Engineering and Design, Volume 89, Issue 11, 2014, Pages 2709-2715,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2014.07.009.

        Nevins, W. M. "Summary report: ITER specialists' meeting on heating and current drive."
        ITER-TN-PH-8-4, June 1988. 1988.

        """
        # Calculate peak electron beta at plasma centre, this is not the form used in the paper
        # The paper assumes parabolic profiles for calculating core values with the profile indexes.
        # We instead use the directly calculated electron density and temperature values at the core.
        # So that it is compatible with all profiles

        betae0 = (
            physics_variables.nd_plasma_electron_on_axis
            * physics_variables.temp_plasma_electron_on_axis_kev
            * 1.0e3
            * constants.ELECTRON_CHARGE
            / (b_plasma_toroidal_on_axis**2 / (2.0 * constants.RMU0))
        )

        # Call integration routine using definite integral routine from scipy

        ainteg, _ = integrate.quad(
            lambda y: self._nevins_integral(
                y,
                nd_plasma_electrons_vol_avg,
                te,
                b_plasma_toroidal_on_axis,
                rminor,
                rmajor,
                zeff,
                alphat,
                alphan,
                q0,
                q95,
                beta_toroidal,
            ),
            0,  # Lower bound
            1.0,  # Upper bound
        )

        # Calculate bootstrap current and fraction

        aibs = 2.5 * betae0 * rmajor * b_plasma_toroidal_on_axis * q95 * ainteg
        return 1.0e6 * aibs / plasma_current

    @staticmethod
    def bootstrap_fraction_sakai(
        beta_poloidal: float,
        q95: float,
        q0: float,
        alphan: float,
        alphat: float,
        eps: float,
        ind_plasma_internal_norm: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the Sakai formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        q95 (float): Safety factor at 95% of the plasma radius.
        q0 (float): Safety factor at the magnetic axis.
        alphan (float): Density profile index
        alphat (float): Temperature profile index
        eps (float): Inverse aspect ratio.
        ind_plasma_internal_norm (float): Plasma normalised internal inductance

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
        The profile assumed for the alphan and alpat indexes is only a parabolic profile without a pedestal (L-mode).
        The Root Mean Squared Error for the fitting database of this formula was 0.025
        Concentrating on the positive shear plasmas using the ACCOME code equilibria with the fully non-inductively driven
        conditions with neutral beam (NB) injection only are calculated.
        The electron temperature and the ion temperature were assumed to be equal
        This can be used for all apsect ratios.
        The diamagnetic fraction is included in this formula.

        References:
        Ryosuke Sakai, Takaaki Fujita, Atsushi Okamoto, Derivation of bootstrap current fraction scaling formula for 0-D system code analysis,
        Fusion Engineering and Design, Volume 149, 2019, 111322, ISSN 0920-3796,
        https://doi.org/10.1016/j.fusengdes.2019.111322.
        """
        # Sakai states that the ACCOME dataset used has the toridal diamagnetic current included in the bootstrap current
        # So the diamganetic current should not be calculated with this. i_diamagnetic_current = 0
        return (
            10 ** (0.951 * eps - 0.948)
            * beta_poloidal ** (1.226 * eps + 1.584)
            * ind_plasma_internal_norm ** (-0.184 * eps - 0.282)
            * (q95 / q0) ** (-0.042 * eps - 0.02)
            * alphan ** (0.13 * eps + 0.05)
            * alphat ** (0.502 * eps - 0.273)
        )

    @staticmethod
    def bootstrap_fraction_aries(
        beta_poloidal: float,
        ind_plasma_internal_norm: float,
        core_density: float,
        average_density: float,
        inverse_aspect: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the ARIES formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        ind_plasma_internal_norm (float): Plasma normalized internal inductance.
        core_density (float): Core plasma density.
        average_density (float): Average plasma density.
        inverse_aspect (float): Inverse aspect ratio.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - The source reference does not provide any info about the derivation of the formula. It is only stated

        References:
            - Zoran Dragojlovic et al., “An advanced computational algorithm for systems analysis of tokamak power plants,”
              Fusion Engineering and Design, vol. 85, no. 2, pp. 243-265, Apr. 2010,
              doi: https://doi.org/10.1016/j.fusengdes.2010.02.015.

        """
        # Using the standard variable naming from the ARIES paper
        a_1 = (
            1.10 - 1.165 * ind_plasma_internal_norm + 0.47 * ind_plasma_internal_norm**2
        )
        b_1 = (
            0.806
            - 0.885 * ind_plasma_internal_norm
            + 0.297 * ind_plasma_internal_norm**2
        )

        c_bs = a_1 + b_1 * (core_density / average_density)

        return c_bs * np.sqrt(inverse_aspect) * beta_poloidal

    @staticmethod
    def bootstrap_fraction_andrade(
        beta_poloidal: float,
        core_pressure: float,
        average_pressure: float,
        inverse_aspect: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the Andrade et al formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        core_pressure (float): Core plasma pressure.
        average_pressure (float): Average plasma pressure.
        inverse_aspect (float): Inverse aspect ratio.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - Based off 350 plasma profiles from Experimento Tokamak Esferico (ETE) spherical tokamak
            - A = 1.5, R_0 = 0.3m, I_p = 200kA, B_0=0.4T, beta = 4-10%. Profiles taken as Gaussian shaped functions.
            - Errors mostly up to the order of 10% are obtained when both expressions are compared with the equilibrium estimates for the
              bootstrap current in ETE

        References:
            - M. C. R. Andrade and G. O. Ludwig, “Scaling of bootstrap current on equilibrium and plasma profile parameters in tokamak plasmas,”
              Plasma Physics and Controlled Fusion, vol. 50, no. 6, pp. 065001-065001, Apr. 2008,
              doi: https://doi.org/10.1088/0741-3335/50/6/065001.

        """

        # Using the standard variable naming from the Andrade et.al. paper
        c_p = core_pressure / average_pressure

        # Error +- 0.0007
        c_bs = 0.2340

        return c_bs * np.sqrt(inverse_aspect) * beta_poloidal * c_p**0.8

    @staticmethod
    def bootstrap_fraction_hoang(
        beta_poloidal: float,
        pressure_index: float,
        current_index: float,
        inverse_aspect: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the Hoang et al formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        pressure_index (float): Pressure profile index.
        current_index (float): Current profile index.
        inverse_aspect (float): Inverse aspect ratio.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - Based off of TFTR data calculated using the TRANSP plasma analysis code
            - 170 discharges which was assembled to  study the tritium influx and transport in discharges with D-only neutral beam
              injection (NBI)
            - Contains L-mode, supershots, reversed shear, enhanced reversed shear and increased li discharges
            - Discharges with monotonic flux profiles with reversed shear are also included
            - Is applied to circular cross-section plasmas

        References:
            - G. T. Hoang and R. V. Budny, “The bootstrap fraction in TFTR,” AIP conference proceedings,
              Jan. 1997, doi: https://doi.org/10.1063/1.53414.
        """

        # Using the standard variable naming from the Hoang et.al. paper
        # Hoang et.al uses a different definition for the profile indexes such that
        # alpha_p is defined as the ratio of the central and the volume-averaged values, and the peakedness of the density of the total plasma current
        # (defined as ratio of the central value and I_p), alpha_j$

        # We assume the pressure and current profile is parabolic and use the (profile_index +1) term in lieu
        # The term represents the ratio of the the core to volume averaged value

        # This could lead to large changes in the value depending on interpretation of the profile index

        c_bs = np.sqrt((pressure_index + 1) / (current_index + 1))

        return 0.4 * np.sqrt(inverse_aspect) * beta_poloidal**0.9 * c_bs

    @staticmethod
    def bootstrap_fraction_wong(
        beta_poloidal: float,
        density_index: float,
        temperature_index: float,
        inverse_aspect: float,
        elongation: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the Wong et al formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        density_index (float): Density profile index.
        temperature_index (float): Temperature profile index.
        inverse_aspect (float): Inverse aspect ratio.
        elongation (float): Plasma elongation.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - Data is based off of equilibria from Miller et al.
            - A: 1.2 - 3.0 and stable to n ballooning and low n kink modes at a bootstrap fraction of 99% for kappa = 2, 2.5 and 3
            - The results were parameterized as a function of aspect ratio and elongation
            - The parametric dependency of beta_p and beta_T are based on fitting of the DIII-D high equivalent DT yield results
            - Parabolic profiles should be used for best results as the pressure peaking value is calculated as the product of a parabolic
              temperature and density profile

        References:
            - C.-P. Wong, J. C. Wesley, R. D. Stambaugh, and E. T. Cheng, “Toroidal reactor designs as a function of aspect ratio and elongation,”
              vol. 42, no. 5, pp. 547-556, May 2002, doi: https://doi.org/10.1088/0029-5515/42/5/307.

            - Miller, R L, "Stable bootstrap-current driven equilibria for low aspect ratio tokamaks".
              Switzerland: N. p., 1996. Web.https://fusion.gat.com/pubs-ext/MISCONF96/A22433.pdf
        """
        # Using the standard variable naming from the Wong et.al. paper
        f_peak = 2.0 / scipy.special.beta(0.5, density_index + temperature_index + 1)

        c_bs = 0.773 + 0.019 * elongation

        return c_bs * f_peak**0.25 * beta_poloidal * np.sqrt(inverse_aspect)

    @staticmethod
    def bootstrap_fraction_gi_I(  # noqa: N802
        beta_poloidal: float,
        pressure_index: float,
        temperature_index: float,
        inverse_aspect: float,
        effective_charge: float,
        q95: float,
        q0: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the first scaling from the Gi et al formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        pressure_index (float): Pressure profile index.
        temperature_index (float): Temperature profile index.
        inverse_aspect (float): Inverse aspect ratio.
        effective_charge (float): Plasma effective charge.
        q95 (float): Safety factor at 95% of the plasma radius.
        q0 (float): Safety factor at the magnetic axis.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - Scaling found by solving the Hirshman-Sigmar bootstrap modelusing the matrix inversion method
            - Method was done to put the scaling into parameters compatible with the TPC systems code
            - Uses the ACCOME code to create bootstrap current fractions without using the itrative calculations of the
              curent drive and equilibrium models in the scan
            - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, traing = 0.3, alpha_n = 0.1 - 0.8, alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0
            - Uses parabolic plasma profiles only.
            - Scaling 1 has better accuracy than Scaling 2. However, Scaling 1 overestimated the f_BS value for reversed shear
              equilibrium.

        References:
            - K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, “Bootstrap current fraction scaling for a tokamak reactor design study,”
              Fusion Engineering and Design, vol. 89, no. 11, pp. 2709-2715, Aug. 2014,
              doi: https://doi.org/10.1016/j.fusengdes.2014.07.009.
        """

        # Using the standard variable naming from the Gi et.al. paper

        c_bs = (
            0.474
            * inverse_aspect**-0.1
            * pressure_index**0.974
            * temperature_index**-0.416
            * effective_charge**0.178
            * (q95 / q0) ** -0.133
        )

        return c_bs * np.sqrt(inverse_aspect) * beta_poloidal

    @staticmethod
    def bootstrap_fraction_gi_II(  # noqa: N802
        beta_poloidal: float,
        pressure_index: float,
        temperature_index: float,
        inverse_aspect: float,
        effective_charge: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the second scaling from the Gi et al formula.

        Parameters:
        beta_poloidal (float): Plasma poloidal beta.
        pressure_index (float): Pressure profile index.
        temperature_index (float): Temperature profile index.
        inverse_aspect (float): Inverse aspect ratio.
        effective_charge (float): Plasma effective charge.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
            - Scaling found by solving the Hirshman-Sigmar bootstrap modelusing the matrix inversion method
            - Method was done to put the scaling into parameters compatible with the TPC systems code
            - Uses the ACCOME code to create bootstrap current fractions without using the itrative calculations of the
              curent drive and equilibrium models in the scan
            - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, traing = 0.3, alpha_n = 0.1 - 0.8, alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0
            - Uses parabolic plasma profiles only.
            - This scaling has the q profile dependance removed to obtain a scaling formula with much more flexible variables than
              that by a single profile factor for internal current profile.

        References:
            - K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, “Bootstrap current fraction scaling for a tokamak reactor design study,”
              Fusion Engineering and Design, vol. 89, no. 11, pp. 2709-2715, Aug. 2014,
              doi: https://doi.org/10.1016/j.fusengdes.2014.07.009.
        """

        # Using the standard variable naming from the Gi et.al. paper

        c_bs = (
            0.382
            * inverse_aspect**-0.242
            * pressure_index**0.974
            * temperature_index**-0.416
            * effective_charge**0.178
        )

        return c_bs * np.sqrt(inverse_aspect) * beta_poloidal

    @staticmethod
    def bootstrap_fraction_sugiyama_l_mode(
        eps: float,
        beta_poloidal: float,
        alphan: float,
        alphat: float,
        zeff: float,
        q95: float,
        q0: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the L-mode scaling from the Sugiyama et al formula.

        :param eps: Inverse aspect ratio.
        :type eps: float
        :param beta_poloidal: Plasma poloidal beta.
        :type beta_poloidal: float
        :param alphan: Density profile index.
        :type alphan: float
        :param alphat: Temperature profile index.
        :type alphat: float
        :param zeff: Plasma effective charge.
        :type zeff: float
        :param q95: Safety factor at 95% of the plasma radius.
        :type q95: float
        :param q0: Safety factor at the magnetic axis.
        :type q0: float

        :returns: The calculated bootstrap fraction.
        :rtype: float

        :notes:
            - This scaling is derived for L-mode plasmas.
            - Ion and electron temperature are the same
            - Z_eff has a uniform profile, with only fully stripped carbon impurity

        :references:
            - S. Sugiyama, T. Goto, H. Utoh, and Y. Sakamoto, “Improvement of core plasma power and
              current balance models for tokamak systems code considering H-mode plasma profiles,”
              Fusion Engineering and Design, vol. 216, p. 115022, Jul. 2025, doi:
              https://doi.org/10.1016/j.fusengdes.2025.115022.
        """

        return (
            0.740
            * eps**0.418
            * beta_poloidal**0.904
            * alphan**0.06
            * alphat**-0.138
            * zeff**0.230
            * (q95 / q0) ** -0.142
        )

    @staticmethod
    def bootstrap_fraction_sugiyama_h_mode(
        eps: float,
        beta_poloidal: float,
        alphan: float,
        alphat: float,
        tbeta: float,
        zeff: float,
        q95: float,
        q0: float,
        radius_plasma_pedestal_density_norm: float,
        nd_plasma_pedestal_electron: float,
        n_greenwald: float,
        temp_plasma_pedestal_kev: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the H-mode scaling from the Sugiyama et al formula.

        :param eps: Inverse aspect ratio.
        :type eps: float
        :param beta_poloidal: Plasma poloidal beta.
        :type beta_poloidal: float
        :param alphan: Density profile index.
        :type alphan: float
        :param alphat: Temperature profile index.
        :type alphat: float
        :param tbeta: Second temperature profile index.
        :type tbeta: float
        :param zeff: Plasma effective charge.
        :type zeff: float
        :param q95: Safety factor at 95% of the plasma radius.
        :type q95: float
        :param q0: Safety factor at the magnetic axis.
        :type q0: float
        :param radius_plasma_pedestal_density_norm: Normalised plasma radius of density pedestal.
        :type radius_plasma_pedestal_density_norm: float
        :param nd_plasma_pedestal_electron: Electron number density at the pedestal [m^-3].
        :type nd_plasma_pedestal_electron: float
        :param n_greenwald: Greenwald density limit [m^-3].
        :type n_greenwald: float
        :param temp_plasma_pedestal_kev: Electron temperature at the pedestal [keV].
        :type temp_plasma_pedestal_kev: float

        :returns: The calculated bootstrap fraction.
        :rtype: float

        :notes:
            - This scaling is derived for H-mode plasmas.
            - The temperature and density pedestal positions are the same
            - Separatrix temperature and density are zero
            - Ion and electron temperature are the same
            - Z_eff has a uniform profile, with only fully stripped carbon impurity

        :references:
            - S. Sugiyama, T. Goto, H. Utoh, and Y. Sakamoto, “Improvement of core plasma power and
              current balance models for tokamak systems code considering H-mode plasma profiles,”
              Fusion Engineering and Design, vol. 216, p. 115022, Jul. 2025, doi:
              https://doi.org/10.1016/j.fusengdes.2025.115022.
        """

        return (
            0.789
            * eps**0.606
            * beta_poloidal**0.960
            * alphan**0.0319
            * alphat**0.00822
            * tbeta**-0.0783
            * zeff**0.241
            * (q95 / q0) ** -0.103
            * radius_plasma_pedestal_density_norm**0.367
            * (nd_plasma_pedestal_electron / n_greenwald) ** -0.174
            * temp_plasma_pedestal_kev**0.0552
        )

    def output_bootstrap_current_information(self):
        """Output the calculated bootstrap current information to the output file."""

        po.ovarrf(
            self.outfile,
            "Bootstrap current fraction multiplier",
            "(cboot)",
            current_drive_variables.cboot,
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (ITER 1989)",
            "(f_c_plasma_bootstrap_iter89)",
            current_drive_variables.f_c_plasma_bootstrap_iter89,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Sauter et al)",
            "(f_c_plasma_bootstrap_sauter)",
            current_drive_variables.f_c_plasma_bootstrap_sauter,
            "OP ",
        )
        for point in range(len(physics_variables.j_plasma_bootstrap_sauter_profile)):
            po.ovarrf(
                self.mfile,
                f"Sauter et al bootstrap current density profile at point {point}",
                f"(j_plasma_bootstrap_sauter_profile{point})",
                physics_variables.j_plasma_bootstrap_sauter_profile[point],
                "OP ",
            )

        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Nevins et al)",
            "(f_c_plasma_bootstrap_nevins)",
            current_drive_variables.f_c_plasma_bootstrap_nevins,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Wilson)",
            "(f_c_plasma_bootstrap_wilson)",
            current_drive_variables.f_c_plasma_bootstrap_wilson,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Sakai)",
            "(f_c_plasma_bootstrap_sakai)",
            current_drive_variables.f_c_plasma_bootstrap_sakai,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (ARIES)",
            "(f_c_plasma_bootstrap_aries)",
            current_drive_variables.f_c_plasma_bootstrap_aries,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Andrade)",
            "(f_c_plasma_bootstrap_andrade)",
            current_drive_variables.f_c_plasma_bootstrap_andrade,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Hoang)",
            "(f_c_plasma_bootstrap_hoang)",
            current_drive_variables.f_c_plasma_bootstrap_hoang,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Wong)",
            "(f_c_plasma_bootstrap_wong)",
            current_drive_variables.f_c_plasma_bootstrap_wong,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Gi I)",
            "(bscf_gi_i)",
            current_drive_variables.bscf_gi_i,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Gi II)",
            "(bscf_gi_ii)",
            current_drive_variables.bscf_gi_ii,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Sugiyama L-mode)",
            "(f_c_plasma_bootstrap_sugiyama_l)",
            current_drive_variables.f_c_plasma_bootstrap_sugiyama_l,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (Sugiyama H-mode)",
            "(f_c_plasma_bootstrap_sugiyama_h)",
            current_drive_variables.f_c_plasma_bootstrap_sugiyama_h,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction (Hender)",
            "(f_c_plasma_diamagnetic_hender)",
            current_drive_variables.f_c_plasma_diamagnetic_hender,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction (SCENE)",
            "(f_c_plasma_diamagnetic_scene)",
            current_drive_variables.f_c_plasma_diamagnetic_scene,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pfirsch-Schlueter fraction (SCENE)",
            "(f_c_plasma_pfirsch_schluter_scene)",
            current_drive_variables.f_c_plasma_pfirsch_schluter_scene,
            "OP ",
        )
        # Error to catch if bootstap fraction limit has been enforced
        if physics_variables.err242 == 1:
            logger.error("Bootstrap fraction upper limit enforced")

        # Error to catch if self-driven current fraction limit has been enforced
        if physics_variables.err243 == 1:
            logger.error(
                "Predicted plasma driven current is more than upper limit on non-inductive fraction"
            )

        if physics_variables.i_bootstrap_current == 0:
            po.ocmmnt(self.outfile, "  (User-specified bootstrap current fraction used)")
        elif physics_variables.i_bootstrap_current == 1:
            po.ocmmnt(
                self.outfile, "  (ITER 1989 bootstrap current fraction model used)"
            )
        elif physics_variables.i_bootstrap_current == 2:
            po.ocmmnt(
                self.outfile,
                "  (Nevins et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 3:
            po.ocmmnt(self.outfile, "  (Wilson bootstrap current fraction model used)")
        elif physics_variables.i_bootstrap_current == 4:
            po.ocmmnt(
                self.outfile,
                "  (Sauter et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 5:
            po.ocmmnt(
                self.outfile,
                "  (Sakai et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 6:
            po.ocmmnt(
                self.outfile,
                "  (ARIES bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 7:
            po.ocmmnt(
                self.outfile,
                "  (Andrade et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 8:
            po.ocmmnt(
                self.outfile,
                "  (Hoang et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 9:
            po.ocmmnt(
                self.outfile,
                "  (Wong et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 10:
            po.ocmmnt(
                self.outfile,
                "  (Gi-I et al bootstrap current fraction model used)",
            )
        elif physics_variables.i_bootstrap_current == 11:
            po.ocmmnt(
                self.outfile,
                "  (Gi-II et al bootstrap current fraction model used)",
            )

        if physics_variables.i_diamagnetic_current == 0:
            po.ocmmnt(self.outfile, "  (Diamagnetic current fraction not calculated)")
            # Error to show if diamagnetic current is above 1% but not used
            if current_drive_variables.f_c_plasma_diamagnetic_scene > 0.01e0:
                logger.error(
                    "Diamagnetic fraction is more than 1%, but not calculated. "
                    "Consider using i_diamagnetic_current=2 and i_pfirsch_schluter_current=1"
                )

        elif physics_variables.i_diamagnetic_current == 1:
            po.ocmmnt(
                self.outfile, "  (Hender diamagnetic current fraction scaling used)"
            )
        elif physics_variables.i_diamagnetic_current == 2:
            po.ocmmnt(
                self.outfile, "  (SCENE diamagnetic current fraction scaling used)"
            )

        if physics_variables.i_pfirsch_schluter_current == 0:
            po.ocmmnt(self.outfile, "  Pfirsch-Schluter current fraction not calculated")
        elif physics_variables.i_pfirsch_schluter_current == 1:
            po.ocmmnt(
                self.outfile,
                "  (SCENE Pfirsch-Schluter current fraction scaling used)",
            )

        po.ovarrf(
            self.outfile,
            "Bootstrap fraction (enforced)",
            "(f_c_plasma_bootstrap.)",
            current_drive_variables.f_c_plasma_bootstrap,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction (enforced)",
            "(f_c_plasma_diamagnetic.)",
            current_drive_variables.f_c_plasma_diamagnetic,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pfirsch-Schlueter fraction (enforced)",
            "(f_c_plasma_pfirsch_schluter.)",
            current_drive_variables.f_c_plasma_pfirsch_schluter,
            "OP ",
        )


class SauterBootstrapCurrent:
    """Class to calculate the bootstrap current using the Sauter et al formula."""

    def bootstrap_fraction_sauter(self, plasma_profile: float) -> float:
        """Calculate the bootstrap current fraction from the Sauter et al scaling.

        Args:
            plasma_profile (PlasmaProfile): The plasma profile object.

        Returns:
            float: The bootstrap current fraction.

        This function calculates the bootstrap current fraction using the Sauter, Angioni, and Lin-Liu scaling.

        Reference:
        - O. Sauter, C. Angioni, Y. R. Lin-Liu;
          Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime.
          Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        - O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
          Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime
          [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140. https://doi.org/10.1063/1.1517052

        Note:
        The code was supplied by Emiliano Fable, IPP Garching (private communication).
        """

        # Radial points from 0 to 1 seperated by 1/profile_size
        roa = plasma_profile.neprofile.profile_x

        # Local circularised minor radius
        rho = np.sqrt(physics_variables.a_plasma_poloidal / np.pi) * roa

        # Square root of local aspect ratio
        sqeps = np.sqrt(roa * (physics_variables.rminor / physics_variables.rmajor))

        # Calculate electron and ion density profiles
        ne = plasma_profile.neprofile.profile_y * 1e-19
        ni = (
            physics_variables.nd_plasma_ions_total_vol_avg
            / physics_variables.nd_plasma_electrons_vol_avg
        ) * ne

        # Calculate electron and ion temperature profiles
        tempe = plasma_profile.teprofile.profile_y
        tempi = (
            physics_variables.temp_plasma_ion_vol_avg_kev
            / physics_variables.temp_plasma_electron_vol_avg_kev
        ) * tempe

        # Flat Zeff profile assumed
        # Return tempi like array object filled with zeff
        zeff = np.full_like(tempi, physics_variables.n_charge_plasma_effective_vol_avg)

        # inverse_q = 1/safety factor
        # Parabolic q profile assumed
        inverse_q = 1 / (
            physics_variables.q0
            + (physics_variables.q95 - physics_variables.q0) * roa**2
        )
        # Create new array of average mass of fuel portion of ions
        amain = np.full_like(inverse_q, physics_variables.m_ions_total_amu)

        # Create new array of average main ion charge
        zmain = np.full_like(inverse_q, 1.0 + physics_variables.f_plasma_fuel_helium3)

        # Calculate total bootstrap current (MA) by summing along profiles
        # Looping from 2 because _calculate_l31_coefficient() etc should return 0 @ j == 1
        radial_elements = np.arange(2, plasma_profile.profile_size)

        # Change in localised minor radius to be used as delta term in derivative
        drho = rho[radial_elements] - rho[radial_elements - 1]

        # Area of annulus, assuming circular plasma cross-section
        da = 2 * np.pi * rho[radial_elements - 1] * drho  # area of annulus

        # Create the partial derivatives using numpy gradient (central differences)
        dlogte_drho = np.gradient(np.log(tempe), rho)[radial_elements - 1]
        dlogti_drho = np.gradient(np.log(tempi), rho)[radial_elements - 1]
        dlogne_drho = np.gradient(np.log(ne), rho)[radial_elements - 1]

        jboot = (
            0.5
            * (
                self._calculate_l31_coefficient(
                    radial_elements,
                    plasma_profile.profile_size,
                    physics_variables.rmajor,
                    physics_variables.b_plasma_toroidal_on_axis,
                    physics_variables.triang,
                    ne,
                    ni,
                    tempe,
                    tempi,
                    inverse_q,
                    rho,
                    zeff,
                    sqeps,
                )
                * dlogne_drho
                + self._calculate_l31_32_coefficient(
                    radial_elements,
                    plasma_profile.profile_size,
                    physics_variables.rmajor,
                    physics_variables.b_plasma_toroidal_on_axis,
                    physics_variables.triang,
                    ne,
                    ni,
                    tempe,
                    tempi,
                    inverse_q,
                    rho,
                    zeff,
                    sqeps,
                )
                * dlogte_drho
                + self._calculate_l34_alpha_31_coefficient(
                    radial_elements,
                    plasma_profile.profile_size,
                    physics_variables.rmajor,
                    physics_variables.b_plasma_toroidal_on_axis,
                    physics_variables.triang,
                    inverse_q,
                    sqeps,
                    tempi,
                    tempe,
                    amain,
                    zmain,
                    ni,
                    ne,
                    rho,
                    zeff,
                )
                * dlogti_drho
            )
            * 1.0e6
            * (
                -physics_variables.b_plasma_toroidal_on_axis
                / (0.2 * np.pi * physics_variables.rmajor)
                * rho[radial_elements - 1]
                * inverse_q[radial_elements - 1]
            )
        )  # A/m2

        return (np.sum(da * jboot, axis=0) / physics_variables.plasma_current), jboot

    @staticmethod
    def _coulomb_logarithm_sauter(
        radial_elements: int, tempe: np.ndarray, ne: np.ndarray
    ) -> np.ndarray:
        """Calculate the Coulomb logarithm used in the arrays for the Sauter bootstrap current scaling.

        This function calculates the Coulomb logarithm, which is valid for e-e collisions (T_e > 0.01 keV)
        and for e-i collisions (T_e > 0.01*Zeff^2) (Alexander, 9/5/1994).

        Parameters
        ----------
        radial_elements :
            the radial element indexes in the range 1 to nr
        tempe :
            the electron temperature array
        ne :
            the electron density array

        Returns
        -------
        type
            - np.ndarray, the Coulomb logarithm at each array point

        Reference:
            - C. A. Ordonez, M. I. Molina;
            Evaluation of the Coulomb logarithm using cutoff and screened Coulomb interaction potentials.
            Phys. Plasmas 1 August 1994; 1 (8): 2515-2518. https://doi.org/10.1063/1.870578
            - Y. R. Shen, “Recent advances in nonlinear optics,” Reviews of Modern Physics, vol. 48, no. 1,
            pp. 1-32, Jan. 1976, doi: https://doi.org/10.1103/revmodphys.48.1.

        """
        return (
            15.9
            - 0.5 * np.log(ne[radial_elements - 1])
            + np.log(tempe[radial_elements - 1])
        )

    def _electron_collisions_sauter(
        self, radial_elements: np.ndarray, tempe: np.ndarray, ne: np.ndarray
    ) -> np.ndarray:
        """Calculate the frequency of electron-electron collisions used in the arrays for the Sauter bootstrap current scaling.

        This function calculates the frequency of electron-electron collisions

        Parameters
        ----------
        radial_elements :
            the radial element indexes in the range 1 to nr
        tempe :
            the electron temperature array
        ne :
            the electron density array

        Returns
        -------
        :
            the frequency of electron-electron collisions (Hz)


        Reference:
            - Yushmanov, 25th April 1987 (?), updated by Pereverzev, 9th November 1994 (?)

        """
        return (
            670.0
            * self._coulomb_logarithm_sauter(radial_elements, tempe, ne)
            * ne[radial_elements - 1]
            / (tempe[radial_elements - 1] * np.sqrt(tempe[radial_elements - 1]))
        )

    def _electron_collisionality_sauter(
        self,
        radial_elements: np.ndarray,
        rmajor: float,
        zeff: np.ndarray,
        inverse_q: np.ndarray,
        sqeps: np.ndarray,
        tempe: np.ndarray,
        ne: np.ndarray,
    ) -> np.ndarray:
        """Calculate the electron collisionality used in the arrays for the Sauter bootstrap current scaling.

        Parameters
        ----------
        radial_elements :
            the radial element index in the range 1 to nr
        rmajor :
            the plasma major radius (m)
        zeff :
            the effective charge array
        inverse_q :
            inverse safety factor profile
        sqeps :
            the square root of the inverse aspect ratio array
        tempe :
            the electron temperature array
        ne :
            the electron density array

        Returns
        -------
        :
            the relative frequency of electron collisions

        Reference:
            - Yushmanov, 30th April 1987 (?)

        """
        return (
            self._electron_collisions_sauter(radial_elements, tempe, ne)
            * 1.4
            * zeff[radial_elements - 1]
            * rmajor
            / np.abs(
                inverse_q[radial_elements - 1]
                * (sqeps[radial_elements - 1] ** 3)
                * np.sqrt(tempe[radial_elements - 1])
                * 1.875e7
            )
        )

    @staticmethod
    def _ion_collisions_sauter(
        radial_elements: np.ndarray,
        zeff: np.ndarray,
        ni: np.ndarray,
        tempi: np.ndarray,
        amain: np.ndarray,
    ) -> np.ndarray:
        """Calculate the full frequency of ion collisions used in the arrays for the Sauter bootstrap current scaling.

        This function calculates the full frequency of ion collisions using the Coulomb logarithm of 15.

        Parameters
        ----------
        radial_elements :
            the radial element indexes in the range 1 to nr
        zeff :
            the effective charge array
        ni :
            the ion density array
        tempi :
            the ion temperature array
        amain :
            the atomic mass of the main ion species array

        Returns
        -------
        :
            the full frequency of ion collisions (Hz)
        """
        return (
            zeff[radial_elements - 1] ** 4
            * ni[radial_elements - 1]
            * 322.0
            / (
                tempi[radial_elements - 1]
                * np.sqrt(tempi[radial_elements - 1] * amain[radial_elements - 1])
            )
        )

    def _ion_collisionality_sauter(
        self,
        radial_elements: np.ndarray,
        rmajor: float,
        inverse_q: np.ndarray,
        sqeps: np.ndarray,
        tempi: np.ndarray,
        amain: np.ndarray,
        zeff: np.ndarray,
        ni: np.ndarray,
    ) -> float:
        """Calculate the ion collisionality to be used in the Sauter bootstrap current scaling.

        Parameters
        ----------
        radial_elements :
            the radial element indexes in the range 1 to nr
        rmajor :
            the plasma major radius (m)
        inverse_q :
            inverse safety factor profile
        sqeps :
            the square root of the inverse aspect ratio profile
        tempi :
            the ion temperature profile (keV)
        amain :
            the atomic mass of the main ion species profile
        zeff :
            the effective charge of the main ion species
        ni :
            the ion density profile (/m^3)

        Returns
        -------
        :
            the ion collisionality
        """
        return (
            3.2e-6
            * self._ion_collisions_sauter(radial_elements, zeff, ni, tempi, amain)
            * rmajor
            / (
                np.abs(inverse_q[radial_elements - 1] + 1.0e-4)
                * sqeps[radial_elements - 1] ** 3
                * np.sqrt(tempi[radial_elements - 1] / amain[radial_elements - 1])
            )
        )

    def _calculate_l31_coefficient(
        self,
        radial_elements: np.ndarray,
        number_of_elements: int,
        rmajor: float,
        b_plasma_toroidal_on_axis: float,
        triang: float,
        ne: np.ndarray,
        ni: np.ndarray,
        tempe: np.ndarray,
        tempi: np.ndarray,
        inverse_q: np.ndarray,
        rho: np.ndarray,
        zeff: np.ndarray,
        sqeps: np.ndarray,
    ) -> float:
        """L31 coefficient before Grad(ln(ne)) in the Sauter bootstrap scaling.


        Parameters
        ----------
        radial_elements :
            radial element indexes in range 2 to nr
        number_of_elements :
            maximum value of radial_elements
        rmajor :
            plasma major radius (m)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        triang :
            plasma triangularity
        ne :
            electron density profile (/m^3)
        ni :
            ion density profile (/m^3)
        tempe :
            electron temperature profile (keV)
        tempi :
            ion temperature profile (keV)
        inverse_q :
            inverse safety factor profile
        rho :
            normalized minor radius profile
        zeff :
            effective charge profile
        sqeps :
            square root of inverse aspect ratio profile

        Returns
        -------
        type
            - float, the coefficient scaling grad(ln(ne)) in the Sauter bootstrap current scaling.

            This function calculates the coefficient scaling grad(ln(ne)) in the Sauter bootstrap current scaling.

        Reference:
            - O. Sauter, C. Angioni, Y. R. Lin-Liu;
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime.
            Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
            - O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime
            [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140. https://doi.org/10.1063/1.1517052

        """
        # Prevents first element being 0
        charge_profile = zeff[radial_elements - 1]

        # Calculate trapped particle fraction
        f_trapped = self._trapped_particle_fraction_sauter(
            radial_elements, triang, sqeps
        )

        # Calculated electron collisionality; nu_e*
        electron_collisionality = self._electron_collisionality_sauter(
            radial_elements, rmajor, zeff, inverse_q, sqeps, tempe, ne
        )

        # $f^{31}_{teff}(\nu_{e*})$, Eq.14b
        f31_teff = f_trapped / (
            (1.0 + (1.0 - 0.1 * f_trapped) * np.sqrt(electron_collisionality))
            + (0.5 * (1.0 - f_trapped) * electron_collisionality) / charge_profile
        )

        l31_coefficient = (
            ((1.0 + 1.4 / (charge_profile + 1.0)) * f31_teff)
            - (1.9 / (charge_profile + 1.0) * f31_teff**2)
            + ((0.3 * f31_teff**3 + 0.2 * f31_teff**4) / (charge_profile + 1.0))
        )

        # Corrections suggested by Fable, 15/05/2015
        return l31_coefficient * self._beta_poloidal_total_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            ni,
            tempe,
            tempi,
            inverse_q,
            rho,
        )

    def _calculate_l31_32_coefficient(
        self,
        radial_elements: np.ndarray,
        number_of_elements: int,
        rmajor: float,
        b_plasma_toroidal_on_axis: float,
        triang: float,
        ne: np.ndarray,
        ni: np.ndarray,
        tempe: np.ndarray,
        tempi: np.ndarray,
        inverse_q: np.ndarray,
        rho: np.ndarray,
        zeff: np.ndarray,
        sqeps: np.ndarray,
    ) -> float:
        """L31 & L32 coefficient before Grad(ln(Te)) in the Sauter bootstrap scaling.

        This function calculates the coefficient scaling grad(ln(Te)) in the Sauter bootstrap current scaling.


        Parameters
        ----------
        radial_elements :
            radial element indexes in range 2 to nr
        number_of_elements :
            maximum value of radial_elements
        rmajor :
            plasma major radius (m)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        triang :
            plasma triangularity
        ne :
            electron density profile (/m^3)
        ni :
            ion density profile (/m^3)
        tempe :
            electron temperature profile (keV)
        tempi :
            ion temperature profile (keV)
        inverse_q :
            inverse safety factor profile
        rho :
            normalized minor radius profile
        zeff :
            effective charge profile
        sqeps :
            square root of inverse aspect ratio profile


        Returns
        -------
        :
            the L31 & L32 coefficient scaling grad(ln(Te)) in the Sauter bootstrap current scaling.


        Reference:
            - O. Sauter, C. Angioni, Y. R. Lin-Liu;
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime.
            Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
            - O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime
            [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140. https://doi.org/10.1063/1.1517052

        """

        # Prevents first element being 0
        charge_profile = zeff[radial_elements - 1]

        # Calculate trapped particle fraction
        f_trapped = self._trapped_particle_fraction_sauter(
            radial_elements, triang, sqeps
        )

        # Calculated electron collisionality; nu_e*
        electron_collisionality = self._electron_collisionality_sauter(
            radial_elements, rmajor, zeff, inverse_q, sqeps, tempe, ne
        )

        # $f^{32\_ee}_{teff}(\nu_{e*})$, Eq.15d
        f32ee_teff = f_trapped / (
            1.0
            + 0.26 * (1.0 - f_trapped) * np.sqrt(electron_collisionality)
            + (
                0.18
                * (1.0 - 0.37 * f_trapped)
                * electron_collisionality
                / np.sqrt(charge_profile)
            )
        )

        # $f^{32\_ei}_{teff}(\nu_{e*})$, Eq.15e
        f32ei_teff = f_trapped / (
            (1.0 + (1.0 + 0.6 * f_trapped) * np.sqrt(electron_collisionality))
            + (
                0.85
                * (1.0 - 0.37 * f_trapped)
                * electron_collisionality
                * (1.0 + charge_profile)
            )
        )

        # $F_{32\_ee}(X)$, Eq.15b
        big_f32ee_teff = (
            (
                (0.05 + 0.62 * charge_profile)
                / charge_profile
                / (1.0 + 0.44 * charge_profile)
                * (f32ee_teff - f32ee_teff**4)
            )
            + (
                (f32ee_teff**2 - f32ee_teff**4 - 1.2 * (f32ee_teff**3 - f32ee_teff**4))
                / (1.0 + 0.22 * charge_profile)
            )
            + (1.2 / (1.0 + 0.5 * charge_profile) * f32ee_teff**4)
        )

        # $F_{32\_ei}(Y)$, Eq.15c
        big_f32ei_teff = (
            (
                -(0.56 + 1.93 * charge_profile)
                / charge_profile
                / (1.0 + 0.44 * charge_profile)
                * (f32ei_teff - f32ei_teff**4)
            )
            + (
                4.95
                / (1.0 + 2.48 * charge_profile)
                * (
                    f32ei_teff**2
                    - f32ei_teff**4
                    - 0.55 * (f32ei_teff**3 - f32ei_teff**4)
                )
            )
            - (1.2 / (1.0 + 0.5 * charge_profile) * f32ei_teff**4)
        )

        # big_f32ee_teff + big_f32ei_teff = L32 coefficient

        # Corrections suggested by Fable, 15/05/2015
        return self._beta_poloidal_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            tempe,
            inverse_q,
            rho,
        ) * (big_f32ee_teff + big_f32ei_teff) + self._calculate_l31_coefficient(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            triang,
            ne,
            ni,
            tempe,
            tempi,
            inverse_q,
            rho,
            zeff,
            sqeps,
        ) * self._beta_poloidal_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            tempe,
            inverse_q,
            rho,
        ) / self._beta_poloidal_total_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            ni,
            tempe,
            tempi,
            inverse_q,
            rho,
        )

    def _calculate_l34_alpha_31_coefficient(
        self,
        radial_elements: np.ndarray,
        number_of_elements: int,
        rmajor: float,
        b_plasma_toroidal_on_axis: float,
        triang: float,
        inverse_q: np.ndarray,
        sqeps: np.ndarray,
        tempi: np.ndarray,
        tempe: np.ndarray,
        amain: float,
        zmain: float,
        ni: np.ndarray,
        ne: np.ndarray,
        rho: np.ndarray,
        zeff: np.ndarray,
    ) -> float:
        """L34, alpha and L31 coefficient before Grad(ln(Ti)) in the Sauter bootstrap scaling.

        This function calculates the coefficient scaling grad(ln(Ti)) in the Sauter bootstrap current scaling.


        Parameters
        ----------
        radial_elements :
            radial element indexes in range 2 to nr
        number_of_elements :
            maximum value of radial_elements
        rmajor :
            plasma major radius (m)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        triang :
            plasma triangularity
        inverse_q :
            inverse safety factor profile
        sqeps :
            square root of inverse aspect ratio profile
        tempi :
            ion temperature profile (keV)
        tempe :
            electron temperature profile (keV)
        amain :
            atomic mass of the main ion
        zmain :
            charge of the main ion
        ni :
            ion density profile (/m^3)
        ne :
            electron density profile (/m^3)
        rho :
            normalized minor radius profile
        zeff :
            effective charge profile


        Returns
        -------
        :
            the L34, alpha and L31 coefficient scaling grad(ln(Ti)) in the Sauter bootstrap current scaling.


        Reference:
            - O. Sauter, C. Angioni, Y. R. Lin-Liu;
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime.
            Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
            - O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime
            [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140. https://doi.org/10.1063/1.1517052

        """
        # Prevents first element being 0
        charge_profile = zeff[radial_elements - 1]

        # Calculate trapped particle fraction
        f_trapped = self._trapped_particle_fraction_sauter(
            radial_elements, triang, sqeps
        )

        # Calculated electron collisionality; nu_e*
        electron_collisionality = self._electron_collisionality_sauter(
            radial_elements, rmajor, zeff, inverse_q, sqeps, tempe, ne
        )

        # $f^{34}_{teff}(\nu_{e*})$, Eq.16b
        f34_teff = f_trapped / (
            (1.0 + (1.0 - 0.1 * f_trapped) * np.sqrt(electron_collisionality))
            + 0.5 * (1.0 - 0.5 * f_trapped) * electron_collisionality / charge_profile
        )

        # Eq.16a
        l34_coefficient = (
            ((1.0 + (1.4 / (charge_profile + 1.0))) * f34_teff)
            - ((1.9 / (charge_profile + 1.0)) * f34_teff**2)
            + ((0.3 / (charge_profile + 1.0)) * f34_teff**3)
            + ((0.2 / (charge_profile + 1.0)) * f34_teff**4)
        )

        # $\alpha_0$, Eq.17a
        alpha_0 = (-1.17 * (1.0 - f_trapped)) / (
            1.0 - (0.22 * f_trapped) - 0.19 * f_trapped**2
        )

        # Calculate the ion collisionality
        ion_collisionality = self._ion_collisionality_sauter(
            radial_elements, rmajor, inverse_q, sqeps, tempi, amain, zmain, ni
        )

        # $\alpha(\nu_{i*})$, Eq.17b
        alpha = (
            (alpha_0 + (0.25 * (1.0 - f_trapped**2)) * np.sqrt(ion_collisionality))
            / (1.0 + (0.5 * np.sqrt(ion_collisionality)))
            + (0.315 * ion_collisionality**2 * f_trapped**6)
        ) / (1.0 + (0.15 * ion_collisionality**2 * f_trapped**6))

        # Corrections suggested by Fable, 15/05/2015
        # Below calculates the L34 * alpha + L31 coefficient
        return (
            self._beta_poloidal_total_sauter(
                radial_elements,
                number_of_elements,
                rmajor,
                b_plasma_toroidal_on_axis,
                ne,
                ni,
                tempe,
                tempi,
                inverse_q,
                rho,
            )
            - self._beta_poloidal_sauter(
                radial_elements,
                number_of_elements,
                rmajor,
                b_plasma_toroidal_on_axis,
                ne,
                tempe,
                inverse_q,
                rho,
            )
        ) * (l34_coefficient * alpha) + self._calculate_l31_coefficient(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            triang,
            ne,
            ni,
            tempe,
            tempi,
            inverse_q,
            rho,
            zeff,
            sqeps,
        ) * (
            1.0
            - self._beta_poloidal_sauter(
                radial_elements,
                number_of_elements,
                rmajor,
                b_plasma_toroidal_on_axis,
                ne,
                tempe,
                inverse_q,
                rho,
            )
            / self._beta_poloidal_total_sauter(
                radial_elements,
                number_of_elements,
                rmajor,
                b_plasma_toroidal_on_axis,
                ne,
                ni,
                tempe,
                tempi,
                inverse_q,
                rho,
            )
        )

    @staticmethod
    def _beta_poloidal_sauter(
        radial_elements: np.ndarray,
        nr: int,
        rmajor: float,
        b_plasma_toroidal_on_axis: float,
        ne: np.ndarray,
        tempe: np.ndarray,
        inverse_q: np.ndarray,
        rho: np.ndarray,
    ) -> np.ndarray:
        """Calculate the local beta poloidal using only electron profiles for the Sauter bootstrap current scaling.

        Parameters
        ----------
        radial_elements :
            radial element indexes in range 1 to nr
        nr :
            maximum value of radial_elements
        rmajor :
            plasma major radius (m)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        ne :
            electron density profile (/m^3)
        tempe :
            electron temperature profile (keV)
        inverse_q :
            inverse safety factor profile
        rho :
            normalized minor radius profile


        Returns
        -------
        :
            the local beta poloidal
        """
        return (
            np.where(
                radial_elements != nr,
                1.6e-4
                * np.pi
                * (ne[radial_elements] + ne[radial_elements - 1])
                * (tempe[radial_elements] + tempe[radial_elements - 1]),
                6.4e-4 * np.pi * ne[radial_elements - 1] * tempe[radial_elements - 1],
            )
            * (
                rmajor
                / (
                    b_plasma_toroidal_on_axis
                    * rho[radial_elements - 1]
                    * np.abs(inverse_q[radial_elements - 1] + 1.0e-4)
                )
            )
            ** 2
        )

    @staticmethod
    def _beta_poloidal_total_sauter(
        radial_elements: np.ndarray,
        nr: int,
        rmajor: float,
        b_plasma_toroidal_on_axis: float,
        ne: np.ndarray,
        ni: np.ndarray,
        tempe: np.ndarray,
        tempi: np.ndarray,
        inverse_q: np.ndarray,
        rho: np.ndarray,
    ) -> np.ndarray:
        """Calculate the local beta poloidal including ion pressure for the Sauter bootstrap current scaling.

        Parameters
        ----------
        radial_elements :
            radial element indexes in range 1 to nr
        nr :
            maximum value of radial_elements
        rmajor :
            plasma major radius (m)
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        ne :
            electron density profile (/m^3)
        ni :
            ion density profile (/m^3)
        tempe :
            electron temperature profile (keV)
        tempi :
            ion temperature profile (keV)
        inverse_q :
            inverse safety factor profile
        rho :
            normalized minor radius profile

        Returns
        -------
        :
            the local total beta poloidal
        """
        return (
            np.where(
                radial_elements != nr,
                1.6e-4
                * np.pi
                * (
                    (
                        (ne[radial_elements] + ne[radial_elements - 1])
                        * (tempe[radial_elements] + tempe[radial_elements - 1])
                    )
                    + (
                        (ni[radial_elements] + ni[radial_elements - 1])
                        * (tempi[radial_elements] + tempi[radial_elements - 1])
                    )
                ),
                6.4e-4
                * np.pi
                * (
                    ne[radial_elements - 1] * tempe[radial_elements - 1]
                    + ni[radial_elements - 1] * tempi[radial_elements - 1]
                ),
            )
            * (
                rmajor
                / (
                    b_plasma_toroidal_on_axis
                    * rho[radial_elements - 1]
                    * np.abs(inverse_q[radial_elements - 1] + 1.0e-4)
                )
            )
            ** 2
        )

    @staticmethod
    def _trapped_particle_fraction_sauter(
        radial_elements: np.ndarray, triang: float, sqeps: np.ndarray, fit: int = 0
    ) -> np.ndarray:
        """Calculates the trapped particle fraction to be used in the Sauter bootstrap current scaling.
        Parameters
        ----------
        radial_elements :
            radial element index in range 1 to nr
        triang :
            plasma triangularity
        sqeps :
            square root of local aspect ratio
        fit :
            fit method (1 = ASTRA method, 2 = Equation from Sauter 2002, 3 = Equation from Sauter 2016)

        Returns
        -------
        type
            - list, trapped particle fraction

            This function calculates the trapped particle fraction at a given radius.

        References
        ----------
            Used in this paper:
            - O. Sauter, C. Angioni, Y. R. Lin-Liu;
            Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime.
            Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

            - O. Sauter, R. J. Buttery, R. Felton, T. C. Hender, D. F. Howell, and contributors to the E.-J. Workprogramme,
            “Marginal-limit for neoclassical tearing modes in JET H-mode discharges,”
            Plasma Physics and Controlled Fusion, vol. 44, no. 9, pp. 1999-2019, Aug. 2002,
            doi: https://doi.org/10.1088/0741-3335/44/9/315.

            - O. Sauter, Geometric formulas for system codes including the effect of negative triangularity,
            Fusion Engineering and Design, Volume 112, 2016, Pages 633-645, ISSN 0920-3796,
            https://doi.org/10.1016/j.fusengdes.2016.04.033.

        """
        # Prevent first element from being zero
        sqeps_reduced = sqeps[radial_elements - 1]
        eps = sqeps_reduced**2

        if fit == 0:
            # ASTRA method, from Emiliano Fable, private communication
            # (Excluding h term which dominates for inverse aspect ratios < 0.5,
            # and tends to take the trapped particle fraction to 1)

            zz = 1.0 - eps
            return 1.0 - zz * np.sqrt(zz) / (1.0 + 1.46 * sqeps_reduced)

        if fit == 1:
            # Equation 4 of Sauter 2002; https://doi.org/10.1088/0741-3335/44/9/315.
            # Similar to, but not quite identical to above

            return 1.0 - (
                ((1.0 - eps) ** 2)
                / ((1.0 + 1.46 * sqeps_reduced) * np.sqrt(1.0 - eps**2))
            )

        if fit == 2:
            # Sauter 2016; https://doi.org/10.1016/j.fusengdes.2016.04.033.
            # Includes correction for triangularity

            epseff = 0.67 * (1.0 - 1.4 * triang * np.abs(triang)) * eps

            return 1.0 - (
                ((1.0 - epseff) / (1.0 + 2.0 * np.sqrt(epseff)))
                * (np.sqrt((1.0 - eps) / (1.0 + eps)))
            )

        raise ProcessValueError(f"fit={fit} is not valid. Must be 1, 2, or 3")
