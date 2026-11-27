"""Bootstrap current models and calculations for plasma physics simulations."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import IntEnum, unique
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numba as nb
import numpy as np
import scipy
from scipy import integrate

from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import Model

if TYPE_CHECKING:
    from process.models.physics.plasma_profiles import PlasmaProfile

logger = logging.getLogger(__name__)


@unique
class BootstrapCurrentFractionModel(IntEnum):
    """Bootstrap plasma current fraction (f_BS) model types"""

    USER_INPUT = (0, "User Input")
    ITER_89 = (1, "ITER IPDG89 scaling")
    NEVINS = (2, "Nevins scaling")
    WILSON = (3, "Wilson scaling")
    SAUTER = (4, "Sauter scaling")
    SAKAI = (5, "Sakai scaling")
    ARIES = (6, "Aries scaling")
    ANDRADE = (7, "Andrade scaling")
    HOANG = (8, "Hoang scaling")
    WONG = (9, "Wong scaling")
    GI_1 = (10, "GI 1 scaling")
    GI_2 = (11, "GI 2 scaling")
    SUGIYAMA_L_MODE = (12, "Sugiyama L Mode scaling")
    SUGIYAMA_H_MODE = (13, "Sugiyama H Mode scaling")

    def __new__(cls, value: int, full_name: str):
        """Create a new BootstrapCurrentFractionModel enum instance.

        Parameters
        ----------
        value : int
            The integer value of the enum member.
        full_name : str
            The full name/description of the bootstrap current model.

        Returns
        -------
        BootstrapCurrentFractionModel
            A new enum instance with the specified value and full_name.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        """The full name of the bootstrap current model."""
        return self._full_name_


class PlasmaBootstrapCurrent(Model):
    """Class to hold plasma bootstrap current for plasma processing."""

    def __init__(
        self, plasma_profile: PlasmaProfile, sauter_bootstrap: SauterBootstrapCurrent
    ) -> None:
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile
        self.sauter_bootstrap = sauter_bootstrap

    def run(self) -> None:
        """Calculate bootstrap current fraction using various models.

        Raises
        ------
        ProcessValueError
            If an illegal value of i_bootstrap_current is provided.
        """
        # Calculate bootstrap current fraction using various models
        self.data.current_drive.f_c_plasma_bootstrap_iter89 = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_iter89(
                aspect=self.data.physics.aspect,
                beta=self.data.physics.beta_total_vol_avg,
                b_plasma_toroidal_on_axis=self.data.physics.b_plasma_total,
                plasma_current=self.data.physics.plasma_current,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
                rmajor=self.data.physics.rmajor,
                vol_plasma=self.data.physics.vol_plasma,
            )
        )

        self.data.current_drive.f_c_plasma_bootstrap_nevins = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_nevins(
                alphan=self.data.physics.alphan,
                alphat=self.data.physics.alphat,
                beta_toroidal=self.data.physics.beta_toroidal_vol_avg,
                b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
                nd_plasma_electrons_vol_avg=self.data.physics.nd_plasma_electrons_vol_avg,
                plasma_current=self.data.physics.plasma_current,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
                rmajor=self.data.physics.rmajor,
                rminor=self.data.physics.rminor,
                te=self.data.physics.temp_plasma_electron_vol_avg_kev,
                zeff=self.data.physics.n_charge_plasma_effective_vol_avg,
            )
        )

        # Wilson scaling uses thermal poloidal beta, not total
        self.data.current_drive.f_c_plasma_bootstrap_wilson = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_wilson(
                alphaj=self.data.physics.alphaj,
                alphap=self.data.physics.alphap,
                alphat=self.data.physics.alphat,
                betpth=self.data.physics.beta_thermal_poloidal_vol_avg,
                q0=self.data.physics.q0,
                q95=self.data.physics.q95,
                rmajor=self.data.physics.rmajor,
                rminor=self.data.physics.rminor,
            )
        )

        sauter_outputs = self.bootstrap_fraction_sauter()
        self.data.current_drive.f_c_plasma_bootstrap_sauter = (
            sauter_outputs.f_cur_plasma_bootstrap * self.data.current_drive.cboot
        )
        self.data.physics.j_plasma_bootstrap_sauter_profile = sauter_outputs.j_bootstrap

        self.data.current_drive.f_c_plasma_bootstrap_sakai = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_sakai(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
                alphan=self.data.physics.alphan,
                alphat=self.data.physics.alphat,
                eps=self.data.physics.eps,
                ind_plasma_internal_norm=self.data.physics.ind_plasma_internal_norm,
            )
        )

        self.data.current_drive.f_c_plasma_bootstrap_aries = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_aries(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                ind_plasma_internal_norm=self.data.physics.ind_plasma_internal_norm,
                core_density=self.data.physics.nd_plasma_electron_on_axis,
                average_density=self.data.physics.nd_plasma_electrons_vol_avg,
                inverse_aspect=self.data.physics.eps,
            )
        )

        self.data.current_drive.f_c_plasma_bootstrap_andrade = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_andrade(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                core_pressure=self.data.physics.pres_plasma_thermal_on_axis,
                average_pressure=self.data.physics.pres_plasma_thermal_vol_avg,
                inverse_aspect=self.data.physics.eps,
            )
        )
        self.data.current_drive.f_c_plasma_bootstrap_hoang = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_hoang(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                pressure_index=self.data.physics.alphap,
                current_index=self.data.physics.alphaj,
                inverse_aspect=self.data.physics.eps,
            )
        )
        self.data.current_drive.f_c_plasma_bootstrap_wong = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_wong(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                density_index=self.data.physics.alphan,
                temperature_index=self.data.physics.alphat,
                inverse_aspect=self.data.physics.eps,
                elongation=self.data.physics.kappa,
            )
        )
        self.data.current_drive.bscf_gi_i = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_gi_I(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                pressure_index=self.data.physics.alphap,
                temperature_index=self.data.physics.alphat,
                inverse_aspect=self.data.physics.eps,
                effective_charge=self.data.physics.n_charge_plasma_effective_vol_avg,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
            )
        )

        self.data.current_drive.bscf_gi_ii = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_gi_II(
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                pressure_index=self.data.physics.alphap,
                temperature_index=self.data.physics.alphat,
                inverse_aspect=self.data.physics.eps,
                effective_charge=self.data.physics.n_charge_plasma_effective_vol_avg,
            )
        )
        self.data.current_drive.f_c_plasma_bootstrap_sugiyama_l = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_sugiyama_l_mode(
                eps=self.data.physics.eps,
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                alphan=self.data.physics.alphan,
                alphat=self.data.physics.alphat,
                zeff=self.data.physics.n_charge_plasma_effective_vol_avg,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
            )
        )
        self.data.current_drive.f_c_plasma_bootstrap_sugiyama_h = (
            self.data.current_drive.cboot
            * self.bootstrap_fraction_sugiyama_h_mode(
                eps=self.data.physics.eps,
                beta_poloidal=self.data.physics.beta_poloidal_vol_avg,
                alphan=self.data.physics.alphan,
                alphat=self.data.physics.alphat,
                tbeta=self.data.physics.tbeta,
                zeff=self.data.physics.n_charge_plasma_effective_vol_avg,
                q95=self.data.physics.q95,
                q0=self.data.physics.q0,
                radius_plasma_pedestal_density_norm=self.data.physics.radius_plasma_pedestal_density_norm,
                nd_plasma_pedestal_electron=self.data.physics.nd_plasma_pedestal_electron,
                n_greenwald=self.data.physics.nd_plasma_electron_max_array[6],
                temp_plasma_pedestal_electron_kev=self.data.physics.temp_plasma_pedestal_electron_kev,
            )
        )

        # Calculate beta_norm_max based on i_beta_norm_max
        try:
            model = BootstrapCurrentFractionModel(
                int(self.data.physics.i_bootstrap_current)
            )
            self.data.current_drive.f_c_plasma_bootstrap = (
                self.get_bootstrap_current_fraction_value(model)
            )
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_bootstrap_current",
                i_bootstrap_current=self.data.physics.i_bootstrap_current,
            ) from None

    def get_bootstrap_current_fraction_value(
        self,
        model: BootstrapCurrentFractionModel,
    ) -> float:
        """Get the plasma current bootstrap fraction (f_BS) for the specified model.

        Parameters
        ----------
        model : BootstrapCurrentFractionModel
            The bootstrap current model type.

        Returns
        -------
        float
            The bootstrap current fraction value.
        """
        cd = self.data.current_drive
        bs_cf_m = BootstrapCurrentFractionModel
        model_map = {
            bs_cf_m.USER_INPUT: cd.f_c_plasma_bootstrap,
            bs_cf_m.ITER_89: cd.f_c_plasma_bootstrap_iter89,
            bs_cf_m.NEVINS: cd.f_c_plasma_bootstrap_nevins,
            bs_cf_m.WILSON: cd.f_c_plasma_bootstrap_wilson,
            bs_cf_m.SAUTER: cd.f_c_plasma_bootstrap_sauter,
            bs_cf_m.SAKAI: cd.f_c_plasma_bootstrap_sakai,
            bs_cf_m.ARIES: cd.f_c_plasma_bootstrap_aries,
            bs_cf_m.ANDRADE: cd.f_c_plasma_bootstrap_andrade,
            bs_cf_m.HOANG: cd.f_c_plasma_bootstrap_hoang,
            bs_cf_m.WONG: cd.f_c_plasma_bootstrap_wong,
            bs_cf_m.GI_1: cd.bscf_gi_i,
            bs_cf_m.GI_2: cd.bscf_gi_ii,
            bs_cf_m.SUGIYAMA_L_MODE: cd.f_c_plasma_bootstrap_sugiyama_l,
            bs_cf_m.SUGIYAMA_H_MODE: cd.f_c_plasma_bootstrap_sugiyama_h,
        }
        return model_map[model]

    @staticmethod
    @nb.njit(cache=True)
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
        """Calculate the bootstrap-driven fraction of the plasma current.

        Parameters
        ----------
        aspect : float
            Plasma aspect ratio.
        beta : float
            Plasma total beta.
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        plasma_current : float
            Plasma current (A).
        q95 : float
            Safety factor at 95% surface.
        q0 : float
            Central safety factor.
        rmajor : float
            Plasma major radius (m).
        vol_plasma : float
            Plasma volume (m³).

        Returns
        -------
        float
            The bootstrap-driven fraction of the plasma current.

        Notes
        -----
        This function performs the original ITER calculation of the plasma current
        bootstrap fraction.

        References
        ----------
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
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
    @nb.njit(cache=True)
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
        """Bootstrap current fraction from Wilson et al scaling.

        Parameters
        ----------
        alphaj : float
            Current profile index.
        alphap : float
            Pressure profile index.
        alphat : float
            Temperature profile index.
        betpth : float
            Thermal component of poloidal beta.
        q0 : float
            Safety factor on axis.
        q95 : float
            Edge safety factor.
        rmajor : float
            Major radius (m).
        rminor : float
            Minor radius (m).

        Returns
        -------
        float
            The bootstrap current fraction.

        Raises
        ------
        ProcessValueError
            If illegal profile values are found (NaN errors or aj < 0).

        Notes
        -----
        This function calculates the bootstrap current fraction using the numerically
        fitted algorithm written by Howard Wilson.

        References
        ----------
        AEA FUS 172: Physics Assessment for the European Reactor Study, 1989

        H. R. Wilson, Nuclear Fusion 32 (1992) 257
        """
        term1 = np.log(0.5)
        term2 = np.log(q0 / q95)

        # Re-arranging of parabolic profile to be equal to (r/a)^2 where the profile
        # value is half of the core value

        termp = 1.0 - 0.5 ** (1.0 / alphap)
        termt = 1.0 - 0.5 ** (1.0 / alphat)
        termj = 1.0 - 0.5 ** (1.0 / alphaj)

        # Assuming a parabolic safety factor profile of the form
        # q = q0 + (q95 - q0) * (r/a)^2
        # Substitute (r/a)^2 term from temperature,pressure and current
        # profiles into q profile when values is 50% of core value. Take natural log of
        # q profile over q95 and q0 to get the profile index

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
    @nb.njit(cache=True)
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

        Parameters
        ----------
        y : float
            Abscissa of integration, normalized minor radius.
        nd_plasma_electrons_vol_avg : float
            Volume averaged electron density (/m³).
        te : float
            Volume averaged electron temperature (keV).
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        rminor : float
            Plasma minor radius (m).
        rmajor : float
            Plasma major radius (m).
        zeff : float
            Plasma effective charge.
        alphat : float
            Temperature profile index.
        alphan : float
            Density profile index.
        q0 : float
            Normalized safety factor at the magnetic axis.
        q95 : float
            Normalized safety factor at 95% of the plasma radius.
        beta_toroidal : float
            Toroidal plasma beta.

        Returns
        -------
        float
            The integrand value.

        References
        ----------
        Keii Gi, Makoto Nakamura, Kenji Tobita, Yasushi Ono,
        "Bootstrap current fraction scaling for a tokamak reactor design study,"
        Fusion Engineering and Design, Volume 89, Issue 11, 2014, Pages 2709-2715,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2014.07.009.

        Nevins, W. M. "Summary report: ITER specialists' meeting on heating and current
        drive." ITER-TN-PH-8-4, June 1988.
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
        self,
    ) -> SauterBootstrapCurrentParams:
        """Get the bootstrap current fraction using Sauter et al scaling.

        Returns
        -------
        :
            The calculated bootstrap current parameters.
        """
        return self.sauter_bootstrap.bootstrap_fraction_sauter(
            a_plasma_poloidal=self.data.physics.a_plasma_poloidal,
            rminor=self.data.physics.rminor,
            rmajor=self.data.physics.rmajor,
            triang=self.data.physics.triang,
            n_charge_plasma_effective_vol_avg=self.data.physics.n_charge_plasma_effective_vol_avg,
            nd_plasma_electron_profile=self.plasma_profile.neprofile.profile_y,
            roa=self.plasma_profile.neprofile.profile_x,
            temp_plasma_electron_profile_kev=self.plasma_profile.teprofile.profile_y,
            temp_plasma_ion_profile_kev=self.plasma_profile.tiprofile.profile_y,
            nd_plasma_ions_total_vol_avg=self.data.physics.nd_plasma_ions_total_vol_avg,
            nd_plasma_electrons_vol_avg=self.data.physics.nd_plasma_electrons_vol_avg,
            q0=self.data.physics.q0,
            q95=self.data.physics.q95,
            m_ions_total_amu=self.data.physics.m_ions_total_amu,
            f_plasma_fuel_helium3=self.data.physics.f_plasma_fuel_helium3,
            b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
            plasma_current=self.data.physics.plasma_current,
        )

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
        """Calculate the bootstrap current fraction from Nevins et al scaling.

        Parameters
        ----------
        alphan : float
            Density profile index.
        alphat : float
            Temperature profile index.
        beta_toroidal : float
            Toroidal plasma beta.
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        nd_plasma_electrons_vol_avg : float
            Electron density (/m³).
        plasma_current : float
            Plasma current (A).
        q95 : float
            Safety factor at 95% surface.
        q0 : float
            Central safety factor.
        rmajor : float
            Plasma major radius (m).
        rminor : float
            Plasma minor radius (m).
        te : float
            Volume averaged plasma temperature (keV).
        zeff : float
            Plasma effective charge.

        Returns
        -------
        float
            The bootstrap current fraction.

        Notes
        -----
        This function calculates the bootstrap current fraction using the Nevins et al
        method.

        References
        ----------
        K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, "Bootstrap current fraction
        scaling for a tokamak reactor design study," Fusion Engineering and Design,
        vol. 89, no. 11, pp. 2709-2715, Aug. 2014,
        doi: https://doi.org/10.1016/j.fusengdes.2014.07.009.

        Nevins, W. M. "Summary report: ITER specialists' meeting on heating and current
        drive." ITER-TN-PH-8-4, June 1988.
        """
        # Calculate peak electron beta at plasma centre, this is not the form used in
        # the paper. The paper assumes parabolic profiles for calculating core values
        # with the profile indexes. We instead use the directly calculated electron
        # density and temperature values at the core. So that it is compatible with all
        # profiles.

        betae0 = (
            self.data.physics.nd_plasma_electron_on_axis
            * self.data.physics.temp_plasma_electron_on_axis_kev
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
    @nb.njit(cache=True)
    def bootstrap_fraction_sakai(
        beta_poloidal: float,
        q95: float,
        q0: float,
        alphan: float,
        alphat: float,
        eps: float,
        ind_plasma_internal_norm: float,
    ) -> float:
        """Calculate the bootstrap fraction using the Sakai formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        q95 : float
            Safety factor at 95% of the plasma radius.
        q0 : float
            Safety factor at the magnetic axis.
        alphan : float
            Density profile index.
        alphat : float
            Temperature profile index.
        eps : float
            Inverse aspect ratio.
        ind_plasma_internal_norm : float
            Plasma normalized internal inductance.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - The profile assumed for the alphan and alphat indexes is only a parabolic
          profile without a pedestal (L-mode).
        - The Root Mean Squared Error for the fitting database of this formula was
          0.025.
        - Concentrating on the positive shear plasmas using the ACCOME code equilibria
          with the fully non-inductively driven conditions with neutral beam (NB)
          injection only are calculated.
        - The electron temperature and the ion temperature were assumed to be equal.
        - This can be used for all aspect ratios.
        - The diamagnetic fraction is included in this formula.

        References
        ----------
        R. Sakai, T. Fujita, and A. Okamoto, "Derivation of bootstrap current fraction
        scaling formula for 0-D system code analysis," Fusion Engineering and Design,
        vol. 149, p. 111322, Dec. 2019,
        doi: https://doi.org/10.1016/j.fusengdes.2019.111322.
        """
        # Sakai states that the ACCOME dataset used has the toridal diamagnetic current
        # included in the bootstrap current. So the diamganetic current should not be
        # calculated with this. i_diamagnetic_current = 0
        return (
            10 ** (0.951 * eps - 0.948)
            * beta_poloidal ** (1.226 * eps + 1.584)
            * ind_plasma_internal_norm ** (-0.184 * eps - 0.282)
            * (q95 / q0) ** (-0.042 * eps - 0.02)
            * alphan ** (0.13 * eps + 0.05)
            * alphat ** (0.502 * eps - 0.273)
        )

    @staticmethod
    @nb.njit(cache=True)
    def bootstrap_fraction_aries(
        beta_poloidal: float,
        ind_plasma_internal_norm: float,
        core_density: float,
        average_density: float,
        inverse_aspect: float,
    ) -> float:
        """Calculate the bootstrap fraction using the ARIES formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        ind_plasma_internal_norm : float
            Plasma normalized internal inductance.
        core_density : float
            Core plasma density.
        average_density : float
            Average plasma density.
        inverse_aspect : float
            Inverse aspect ratio.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        The source reference does not provide any info about the derivation of the
        formula. It is only stated.

        References
        ----------
        Zoran Dragojlovic et al., "An advanced computational algorithm for systems
        analysis of tokamak power plants," Fusion Engineering and Design, vol. 85,
        no. 2, pp. 243-265, Apr. 2010,
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
    @nb.njit(cache=True)
    def bootstrap_fraction_andrade(
        beta_poloidal: float,
        core_pressure: float,
        average_pressure: float,
        inverse_aspect: float,
    ) -> float:
        """Calculate the bootstrap fraction using the Andrade et al formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        core_pressure : float
            Core plasma pressure.
        average_pressure : float
            Average plasma pressure.
        inverse_aspect : float
            Inverse aspect ratio.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - Based off 350 plasma profiles from Experimento Tokamak Esferico (ETE)
          spherical tokamak with A = 1.5, R_0 = 0.3m, I_p = 200kA, B_0=0.4T,
          beta = 4-10%.
        - Profiles taken as Gaussian shaped functions.
        - Errors mostly up to the order of 10% are obtained when both expressions are
          compared with the equilibrium estimates for the bootstrap current in ETE.

        References
        ----------
        M. C. R. Andrade and G. O. Ludwig, "Scaling of bootstrap current on equilibrium
        and plasma profile parameters in tokamak plasmas," Plasma Physics and Controlled
        Fusion, vol. 50, no. 6, pp. 065001-065001, Apr. 2008,
        doi: https://doi.org/10.1088/0741-3335/50/6/065001.
        """
        # Using the standard variable naming from the Andrade et.al. paper
        c_p = core_pressure / average_pressure

        # Error +- 0.0007
        c_bs = 0.2340

        return c_bs * np.sqrt(inverse_aspect) * beta_poloidal * c_p**0.8

    @staticmethod
    @nb.njit(cache=True)
    def bootstrap_fraction_hoang(
        beta_poloidal: float,
        pressure_index: float,
        current_index: float,
        inverse_aspect: float,
    ) -> float:
        """Calculate the bootstrap fraction using the Hoang et al formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        pressure_index : float
            Pressure profile index.
        current_index : float
            Current profile index.
        inverse_aspect : float
            Inverse aspect ratio.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - Based off of TFTR data calculated using the TRANSP plasma analysis code.
        - 170 discharges which was assembled to study the tritium influx and transport
          in discharges with D-only neutral beam injection (NBI).
        - Contains L-mode, supershots, reversed shear, enhanced reversed shear and
          increased li discharges.
        - Discharges with monotonic flux profiles with reversed shear are also included.
        - Is applied to circular cross-section plasmas.

        References
        ----------
        G. T. Hoang and R. V. Budny, "The bootstrap fraction in TFTR," AIP conference
        proceedings, Jan. 1997, doi: https://doi.org/10.1063/1.53414.
        """
        # Using the standard variable naming from the Hoang et.al. paper
        # Hoang et.al uses a different definition for the profile indexes such that
        # alpha_p is defined as the ratio of the central and the volume-averaged values,
        # and the peakedness of the density of the total plasma current
        # (defined as ratio of the central value and I_p), alpha_j$

        # We assume the pressure and current profile is parabolic and use the
        # (profile_index +1) term in lieu. The term represents the ratio of the the
        # core to volume averaged value

        # This could lead to large changes in the value depending on interpretation of
        # the profile index

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
        """Calculate the bootstrap fraction using the Wong et al formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        density_index : float
            Density profile index.
        temperature_index : float
            Temperature profile index.
        inverse_aspect : float
            Inverse aspect ratio.
        elongation : float
            Plasma elongation.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - Data is based off of equilibria from Miller et al.
        - A: 1.2 - 3.0 and stable to n ballooning and low n kink modes at a bootstrap
          fraction of 99% for kappa = 2, 2.5 and 3.
        - The results were parameterized as a function of aspect ratio and elongation.
        - The parametric dependency of beta_p and beta_T are based on fitting of the
          DIII-D high equivalent DT yield results.
        - Parabolic profiles should be used for best results as the pressure peaking
          value is calculated as the product of a parabolic temperature and density
          profile.

        References
        ----------
        C.-P. Wong, J. C. Wesley, R. D. Stambaugh, and E. T. Cheng, "Toroidal reactor
        designs as a function of aspect ratio and elongation," vol. 42, no. 5,
        pp. 547-556, May 2002, doi: https://doi.org/10.1088/0029-5515/42/5/307.

        Miller, R L, "Stable bootstrap-current driven equilibria for low aspect ratio
        tokamaks". Switzerland: N. p., 1996.
        Web.https://fusion.gat.com/pubs-ext/MISCONF96/A22433.pdf
        """
        # Using the standard variable naming from the Wong et.al. paper
        f_peak = 2.0 / scipy.special.beta(0.5, density_index + temperature_index + 1)

        c_bs = 0.773 + 0.019 * elongation

        return c_bs * f_peak**0.25 * beta_poloidal * np.sqrt(inverse_aspect)

    @staticmethod
    @nb.njit(cache=True)
    def bootstrap_fraction_gi_I(  # noqa: N802
        beta_poloidal: float,
        pressure_index: float,
        temperature_index: float,
        inverse_aspect: float,
        effective_charge: float,
        q95: float,
        q0: float,
    ) -> float:
        """Calculate the bootstrap fraction using the first scaling from the Gi et al
        formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        pressure_index : float
            Pressure profile index.
        temperature_index : float
            Temperature profile index.
        inverse_aspect : float
            Inverse aspect ratio.
        effective_charge : float
            Plasma effective charge.
        q95 : float
            Safety factor at 95% of the plasma radius.
        q0 : float
            Safety factor at the magnetic axis.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - Scaling found by solving the Hirshman-Sigmar bootstrap model using the matrix
          inversion method.
        - Method was done to put the scaling into parameters compatible with the TPC
          systems code.
        - Uses the ACCOME code to create bootstrap current fractions without using the
          iterative calculations of the current drive and equilibrium models in the
          scan.
        - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, triang = 0.3, alpha_n = 0.1 - 0.8,
          alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0.
        - Uses parabolic plasma profiles only.
        - Scaling 1 has better accuracy than Scaling 2. However, Scaling 1
          overestimated the f_BS value for reversed shear equilibrium.

        References
        ----------
        K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, "Bootstrap current fraction
        scaling for a tokamak reactor design study," Fusion Engineering and Design,
        vol. 89, no. 11, pp. 2709-2715, Aug. 2014,
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
    @nb.njit(cache=True)
    def bootstrap_fraction_gi_II(  # noqa: N802
        beta_poloidal: float,
        pressure_index: float,
        temperature_index: float,
        inverse_aspect: float,
        effective_charge: float,
    ) -> float:
        """Calculate the bootstrap fraction using the second scaling from the Gi et al
        formula.

        Parameters
        ----------
        beta_poloidal : float
            Plasma poloidal beta.
        pressure_index : float
            Pressure profile index.
        temperature_index : float
            Temperature profile index.
        inverse_aspect : float
            Inverse aspect ratio.
        effective_charge : float
            Plasma effective charge.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - Scaling found by solving the Hirshman-Sigmar bootstrap model using the matrix
          inversion method.
        - Method was done to put the scaling into parameters compatible with the TPC
          systems code.
        - Uses the ACCOME code to create bootstrap current fractions without using the
          iterative calculations of the current drive and equilibrium models in the
          scan.
        - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, triang = 0.3, alpha_n = 0.1 - 0.8,
          alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0.
        - Uses parabolic plasma profiles only.
        - This scaling has the q profile dependence removed to obtain a scaling formula
          with much more flexible variables than that by a single profile factor for
          internal current profile.

        References
        ----------
        K. Gi, M. Nakamura, Kenji Tobita, and Y. Ono, "Bootstrap current fraction
        scaling for a tokamak reactor design study," Fusion Engineering and Design,
        vol. 89, no. 11, pp. 2709-2715, Aug. 2014,
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
    @nb.njit(cache=True)
    def bootstrap_fraction_sugiyama_l_mode(
        eps: float,
        beta_poloidal: float,
        alphan: float,
        alphat: float,
        zeff: float,
        q95: float,
        q0: float,
    ) -> float:
        """Calculate the bootstrap fraction using the L-mode scaling from the Sugiyama
        et al formula.

        Parameters
        ----------
        eps : float
            Inverse aspect ratio.
        beta_poloidal : float
            Plasma poloidal beta.
        alphan : float
            Density profile index.
        alphat : float
            Temperature profile index.
        zeff : float
            Plasma effective charge.
        q95 : float
            Safety factor at 95% of the plasma radius.
        q0 : float
            Safety factor at the magnetic axis.

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - This scaling is derived for L-mode plasmas.
        - Ion and electron temperature are the same.
        - Z_eff has a uniform profile, with only fully stripped carbon impurity.

        References
        ----------
        S. Sugiyama, T. Goto, H. Utoh, and Y. Sakamoto, "Improvement of core plasma
        power and current balance models for tokamak systems code considering H-mode
        plasma profiles," Fusion Engineering and Design, vol. 216, p. 115022,
        Jul. 2025, doi:https://doi.org/10.1016/j.fusengdes.2025.115022.
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
        temp_plasma_pedestal_electron_kev: float,
    ) -> float:
        """Calculate the bootstrap fraction using the H-mode scaling from the Sugiyama
        et al formula.

        Parameters
        ----------
        eps : float
            Inverse aspect ratio.
        beta_poloidal : float
            Plasma poloidal beta.
        alphan : float
            Density profile index.
        alphat : float
            Temperature profile index.
        tbeta : float
            Second temperature profile index.
        zeff : float
            Plasma effective charge.
        q95 : float
            Safety factor at 95% of the plasma radius.
        q0 : float
            Safety factor at the magnetic axis.
        radius_plasma_pedestal_density_norm : float
            Normalised plasma radius of density pedestal.
        nd_plasma_pedestal_electron : float
            Electron number density at the pedestal [m⁻³].
        n_greenwald : float
            Greenwald density limit [m⁻³].
        temp_plasma_pedestal_electron_kev : float
            Electron temperature at the pedestal [keV].

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
        - This scaling is derived for H-mode plasmas.
        - The temperature and density pedestal positions are the same.
        - Separatrix temperature and density are zero.
        - Ion and electron temperature are the same.
        - Z_eff has a uniform profile, with only fully stripped carbon impurity.

        References
        ----------
        S. Sugiyama, T. Goto, H. Utoh, and Y. Sakamoto, "Improvement of core plasma
        power and current balance models for tokamak systems code considering H-mode
        plasma profiles," Fusion Engineering and Design, vol. 216, p. 115022,
        Jul. 2025, doi: https://doi.org/10.1016/j.fusengdes.2025.115022.
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
            * temp_plasma_pedestal_electron_kev**0.0552
        )

    def output(self):
        """Output the calculated bootstrap current information to the output file."""
        po.oheadr(self.outfile, "Plasma Bootstrap Current Fraction")
        po.ovarre(
            self.outfile,
            "Plasma bootstrap current fraction scaling used",
            "(i_bootstrap_current)",
            self.data.physics.i_bootstrap_current,
        )
        po.ocmmnt(
            self.outfile,
            "Bootstrap current fraction model selected: "
            f"{BootstrapCurrentFractionModel(self.data.physics.i_bootstrap_current).full_name} ",  # noqa: E501
        )
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Bootstrap current fraction multiplier",
            "(cboot)",
            self.data.current_drive.cboot,
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.ITER_89.full_name})",
            "(f_c_plasma_bootstrap_iter89)",
            self.data.current_drive.f_c_plasma_bootstrap_iter89,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.NEVINS.full_name})",
            "(f_c_plasma_bootstrap_nevins)",
            self.data.current_drive.f_c_plasma_bootstrap_nevins,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.WILSON.full_name})",
            "(f_c_plasma_bootstrap_wilson)",
            self.data.current_drive.f_c_plasma_bootstrap_wilson,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.SAUTER.full_name})",
            "(f_c_plasma_bootstrap_sauter)",
            self.data.current_drive.f_c_plasma_bootstrap_sauter,
            "OP ",
        )
        for point in range(len(self.data.physics.j_plasma_bootstrap_sauter_profile)):
            po.ovarre(
                self.mfile,
                f"Sauter et al bootstrap current density profile at point {point}",
                f"(j_plasma_bootstrap_sauter_profile{point})",
                self.data.physics.j_plasma_bootstrap_sauter_profile[point],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.SAKAI.full_name})",
            "(f_c_plasma_bootstrap_sakai)",
            self.data.current_drive.f_c_plasma_bootstrap_sakai,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.ARIES.full_name})",
            "(f_c_plasma_bootstrap_aries)",
            self.data.current_drive.f_c_plasma_bootstrap_aries,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.ANDRADE.full_name})",
            "(f_c_plasma_bootstrap_andrade)",
            self.data.current_drive.f_c_plasma_bootstrap_andrade,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.HOANG.full_name})",
            "(f_c_plasma_bootstrap_hoang)",
            self.data.current_drive.f_c_plasma_bootstrap_hoang,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.WONG.full_name})",
            "(f_c_plasma_bootstrap_wong)",
            self.data.current_drive.f_c_plasma_bootstrap_wong,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.GI_1.full_name})",
            "(bscf_gi_i)",
            self.data.current_drive.bscf_gi_i,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction ({BootstrapCurrentFractionModel.GI_2.full_name})",
            "(bscf_gi_ii)",
            self.data.current_drive.bscf_gi_ii,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction"
            f" ({BootstrapCurrentFractionModel.SUGIYAMA_L_MODE.full_name})",
            "(f_c_plasma_bootstrap_sugiyama_l)",
            self.data.current_drive.f_c_plasma_bootstrap_sugiyama_l,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"Bootstrap fraction"
            f" ({BootstrapCurrentFractionModel.SUGIYAMA_H_MODE.full_name})",
            "(f_c_plasma_bootstrap_sugiyama_h)",
            self.data.current_drive.f_c_plasma_bootstrap_sugiyama_h,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Pfirsch-Schlueter fraction (SCENE)",
            "(f_c_plasma_pfirsch_schluter_scene)",
            self.data.current_drive.f_c_plasma_pfirsch_schluter_scene,
            "OP ",
        )
        # Error to catch if bootstap fraction limit has been enforced
        if self.data.physics.err242 == 1:
            logger.error("Bootstrap fraction upper limit enforced")

        # Error to catch if self-driven current fraction limit has been enforced
        if self.data.physics.err243 == 1:
            logger.error(
                "Predicted plasma driven current is more than upper limit on "
                "non-inductive fraction"
            )

        if self.data.physics.i_pfirsch_schluter_current == 0:
            po.ocmmnt(self.outfile, "  Pfirsch-Schluter current fraction not calculated")
        elif self.data.physics.i_pfirsch_schluter_current == 1:
            po.ocmmnt(
                self.outfile,
                "  (SCENE Pfirsch-Schluter current fraction scaling used)",
            )

        po.ovarre(
            self.outfile,
            "Bootstrap fraction (enforced)",
            "(f_c_plasma_bootstrap.)",
            self.data.current_drive.f_c_plasma_bootstrap,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Diamagnetic fraction (enforced)",
            "(f_c_plasma_diamagnetic.)",
            self.data.current_drive.f_c_plasma_diamagnetic,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Pfirsch-Schlueter fraction (enforced)",
            "(f_c_plasma_pfirsch_schluter.)",
            self.data.current_drive.f_c_plasma_pfirsch_schluter,
            "OP ",
        )


@dataclass
class SauterBootstrapCurrentParams:
    """Output parameters for the Sauter bootstrap current calculation."""

    f_cur_plasma_bootstrap: float
    "Bootstrap current fraction [-]"
    j_bootstrap: float | np.ndarray
    "Bootstrap current density (j) [A/m²]"
    f_trapped_particles: float | np.ndarray
    "Trapped particle fraction [-]"
    normalised_electron_collisionality: float | np.ndarray
    "Normalised electron collisionality [-]"
    normalised_ion_collisionality: float | np.ndarray
    "Normalised ion collisionality [-]"
    l31_coefficient: np.ndarray
    "L₃₁ coefficient [-]"
    l32_coefficient: np.ndarray
    "L₃₂ coefficient [-]"
    l34_coefficient: np.ndarray
    "L₃₄ coefficient [-]"
    alpha: np.ndarray
    "α coefficient [-]"  # noqa: RUF001
    alpha_0: np.ndarray
    "α₀ coefficient [-]"
    q_profile: float | np.ndarray
    "Safety factor profile (q) [-]"
    inverse_aspect: float | np.ndarray
    "Inverse aspect ratio (ε) [-]"
    sqrt_inverse_aspect: float | np.ndarray
    "Square root of inverse aspect ratio (√ε) [-]"


class SauterBootstrapCurrent(Model):
    """Class to calculate the bootstrap current using the Sauter et al formula.


    References
    ----------
    [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
    Neoclassical conductivity and bootstrap current formulas for general
    axisymmetric equilibria and arbitrary collisionality regime. Phys. Plasmas
    1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

    [2] O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
    Neoclassical conductivity and bootstrap current formulas for general
    axisymmetric equilibria and arbitrary collisionality regime
    [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12):
    5140. https://doi.org/10.1063/1.1517052


    """

    def run(self):
        """SauterBootstrapCurrent model isn't run"""

    def output(self):
        """SauterBootstrapCurrent model doesn't have any output"""

    def bootstrap_fraction_sauter(
        self,
        a_plasma_poloidal: float,
        rminor: float,
        rmajor: float,
        triang: float,
        n_charge_plasma_effective_vol_avg: float,
        nd_plasma_electron_profile: np.ndarray,
        roa: np.ndarray,
        temp_plasma_electron_profile_kev: np.ndarray,
        temp_plasma_ion_profile_kev: np.ndarray,
        nd_plasma_ions_total_vol_avg: float,
        nd_plasma_electrons_vol_avg: float,
        q0: float,
        q95: float,
        m_ions_total_amu: float,  # noqa: ARG002
        f_plasma_fuel_helium3: float,
        b_plasma_toroidal_on_axis: float,
        plasma_current: float,
    ) -> SauterBootstrapCurrentParams:
        """Calculate the bootstrap current fraction from the Sauter et al scaling.

        Parameters
        ----------
        a_plasma_poloidal:
            The poloidal cross-sectional area of the plasma [m²].
        rminor:
            The minor radius of the plasma (a) [m].
        rmajor:
            The major radius of the plasma (R₀) [m].
        triang:
            The triangularity of the plasma [-].
        n_charge_plasma_effective_vol_avg:
            The volume-averaged effective charge of the plasma (Z_eff) [-].
        nd_plasma_electron_profile:
            The electron density profile of the plasma [m⁻³].
        roa:
            Normalised minor-radius profile r/a [-].
        temp_plasma_electron_profile_kev:
            The electron temperature profile of the plasma [keV].
        temp_plasma_ion_profile_kev:
            The ion temperature profile of the plasma [keV].
        nd_plasma_ions_total_vol_avg:
            The volume-averaged total ion density of the plasma [m⁻³].
        nd_plasma_electrons_vol_avg:
            The volume-averaged electron density of the plasma [m⁻³].
        q0:
            The central safety factor [-].
        q95:
            The edge safety factor at 95% of the minor radius [-].
        m_ions_total_amu:
            The total ion mass in atomic mass units (amu) [-]. Retained for
            compatibility; the Sauter collisionality formula does not depend on it.
        f_plasma_fuel_helium3:
            The fraction of helium-3 in the plasma fuel [-].
        b_plasma_toroidal_on_axis:
            The toroidal magnetic field on the plasma axis (B₀) [T].
        plasma_current:
            The total plasma current (I_p) [A].

        Returns
        -------
        :
            The calculated bootstrap current parameters.

        Raises
        ------
        ProcessValueError
            If the collisionality inputs are not positive and finite, or the
            safety factor or elliptical-surface geometry is invalid.

        Notes
        -----
        This function calculates the bootstrap current fraction using the Sauter,
        Angioni, and Lin-Liu scaling. The total ion density, including thermal
        impurities and hot beam ions, is treated as one ion population at the
        thermal ion temperature and with the effective main-ion charge. Fast-alpha
        pressure is not included. This is a multispecies approximation, not the
        single-ion model derived in the paper.

        """
        if not np.isfinite(a_plasma_poloidal) or a_plasma_poloidal <= 0:
            raise ProcessValueError(
                "Sauter elliptical surfaces require a finite, positive plasma area",
                a_plasma_poloidal=a_plasma_poloidal,
            )

        if not np.isfinite(rminor) or rminor <= 0:
            raise ProcessValueError(
                "Sauter elliptical surfaces require a finite, positive minor radius",
                rminor=rminor,
            )
        kappa_effective = a_plasma_poloidal / (np.pi * rminor**2)
        r_horizontal = rminor * roa
        # Area-equivalent radius keeps the existing annulus areas and gradients.
        r_area = np.sqrt(a_plasma_poloidal / np.pi) * roa

        if (
            not np.isfinite(rmajor)
            or rmajor <= 0
            or not np.isfinite(kappa_effective)
            or kappa_effective <= 0
            or not np.all(np.isfinite(r_horizontal))
            or np.any(r_horizontal < 0)
            or np.any(r_horizontal >= rmajor)
        ):
            raise ProcessValueError(
                "Sauter elliptical surfaces require 0 <= r_horizontal < rmajor "
                "and a finite, positive rmajor",
                rmajor=rmajor,
                rminor=rminor,
                a_plasma_poloidal=a_plasma_poloidal,
            )

        # Square root of local aspect ratio
        inverse_aspect = (rminor / rmajor) * roa
        sqrt_inverse_aspect = np.sqrt(inverse_aspect)

        # Prevent log(0) using the existing density floor in m⁻³.
        nd_plasma_electron_profile = np.maximum(nd_plasma_electron_profile, 1.0e13)
        nd_plasma_ions_total_profile = (
            nd_plasma_ions_total_vol_avg / nd_plasma_electrons_vol_avg
        ) * nd_plasma_electron_profile

        # Flat effective charge profile assumed
        n_charge_plasma_effective_profile = np.full_like(
            temp_plasma_ion_profile_kev, n_charge_plasma_effective_vol_avg, dtype=float
        )

        # Parabolic q profile assumed
        if not np.isfinite(q0) or q0 <= 0 or not np.isfinite(q95) or q95 <= 0:
            raise ProcessValueError(
                "Sauter safety factors q0 and q95 must be positive and finite",
                q0=q0,
                q95=q95,
            )
        q_profile = q0 + (q95 - q0) * roa**2
        if not np.all(np.isfinite(q_profile)) or np.any(q_profile <= 0):
            raise ProcessValueError(
                "Sauter safety factor profile must be positive and finite",
                q0=q0,
                q95=q95,
            )
        inverse_q = 1 / q_profile
        # Create new array of average main ion charge
        zmain = np.full_like(inverse_q, 1.0 + f_plasma_fuel_helium3)

        # Calculate total bootstrap current (MA) by summing along profiles
        # Looping from 2 because _calculate_l31_coefficient() etc should return 0
        # @ j == 1
        radial_elements = np.arange(2, len(nd_plasma_electron_profile))

        # Annulus area for self-similar ellipses.
        dr_area = r_area[radial_elements] - r_area[radial_elements - 1]
        da = 2 * np.pi * r_area[radial_elements - 1] * dr_area

        # Create the partial derivatives using numpy gradient (central differences)
        dlogte_dr_area = np.gradient(np.log(temp_plasma_electron_profile_kev), r_area)[
            radial_elements - 1
        ]
        dlogti_dr_area = np.gradient(np.log(temp_plasma_ion_profile_kev), r_area)[
            radial_elements - 1
        ]
        dlogne_dr_area = np.gradient(np.log(nd_plasma_electron_profile), r_area)[
            radial_elements - 1
        ]

        # Calculate trapped particle fraction
        f_trapped_particles = self._trapped_particle_fraction_sauter(
            radial_elements=radial_elements,
            triang=triang,
            sqrt_inverse_aspect=sqrt_inverse_aspect,
        )

        # Calculate normalised electron collisionality
        normalised_electron_collisionality = self._electron_collisionality_sauter(
            radial_elements=radial_elements,
            rmajor=rmajor,
            n_charge_plasma_effective=n_charge_plasma_effective_profile,
            inverse_q=inverse_q,
            sqrt_inverse_aspect=sqrt_inverse_aspect,
            temp_plasma_electron_profile_kev=temp_plasma_electron_profile_kev,
            nd_plasma_electron_profile=nd_plasma_electron_profile,
        )

        # Calculate normalised ion collisionality
        normalised_ion_collisionality = self._ion_collisionality_sauter(
            radial_elements=radial_elements,
            rmajor=rmajor,
            inverse_q=inverse_q,
            sqrt_inverse_aspect=sqrt_inverse_aspect,
            temp_plasma_ion_profile_kev=temp_plasma_ion_profile_kev,
            n_charge_main_ion=zmain,
            nd_plasma_ion_profile=nd_plasma_ions_total_profile,
        )

        l31_coefficient = self._calculate_l31_coefficient(
            radial_elements=radial_elements,
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            zeff=n_charge_plasma_effective_profile,
        )
        l32_coefficient = self._calculate_l32_coefficient(
            radial_elements=radial_elements,
            f_trapped_particles=f_trapped_particles,
            zeff=n_charge_plasma_effective_profile,
            normalised_electron_collisionality=normalised_electron_collisionality,
        )
        l34_coefficient = self._calculate_l34_coefficient(
            radial_elements=radial_elements,
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            zeff=n_charge_plasma_effective_profile,
        )
        alpha_0 = self._calculate_alpha_0_coefficient(
            f_trapped_particles=f_trapped_particles
        )
        alpha = self._calculate_alpha_coefficient(
            alpha_0=alpha_0,
            f_trapped_particles=f_trapped_particles,
            normalised_ion_collisionality=normalised_ion_collisionality,
        )

        previous = radial_elements - 1
        ne_mid = 0.5 * (
            nd_plasma_electron_profile[radial_elements]
            + nd_plasma_electron_profile[previous]
        )
        ni_mid = 0.5 * (
            nd_plasma_ions_total_profile[radial_elements]
            + nd_plasma_ions_total_profile[previous]
        )
        te_mid = 0.5 * (
            temp_plasma_electron_profile_kev[radial_elements]
            + temp_plasma_electron_profile_kev[previous]
        )
        ti_mid = 0.5 * (
            temp_plasma_ion_profile_kev[radial_elements]
            + temp_plasma_ion_profile_kev[previous]
        )
        # Convert density in m⁻³ and temperature in keV to pressure in Pa.
        pressure_electron = constants.KILOELECTRON_VOLT * ne_mid * te_mid
        pressure_ion = constants.KILOELECTRON_VOLT * ni_mid * ti_mid
        pressure_total = pressure_electron + pressure_ion

        # Effective field, not the local field on an elliptical flux surface.
        field_scale = b_plasma_toroidal_on_axis * r_area[previous] / rmajor
        poloidal_field_effective = field_scale * inverse_q[previous]

        # Ion density follows the electron density, so both share this gradient.
        density_gradient_term = pressure_total * l31_coefficient * dlogne_dr_area
        electron_temperature_term = (
            pressure_electron * (l31_coefficient + l32_coefficient) * dlogte_dr_area
        )
        ion_temperature_term = (
            pressure_ion * (l31_coefficient + l34_coefficient * alpha) * dlogti_dr_area
        )

        j_bootstrap_toroidal = (
            -(density_gradient_term + electron_temperature_term + ion_temperature_term)
            / poloidal_field_effective
        )  # A/m2

        # Elliptical surfaces with j_parallel / B constant on each surface.
        epsilon_horizontal = r_horizontal[previous] / rmajor
        geometry_root = np.sqrt(1.0 - epsilon_horizontal**2)
        field_pitch_scale = (
            kappa_effective * epsilon_horizontal * inverse_q[previous] / geometry_root
        )
        shape_average = (1.0 + geometry_root / kappa_effective**2) / (
            1.0 + geometry_root
        )
        j_bootstrap_toroidal *= geometry_root / (
            1.0 + field_pitch_scale**2 * shape_average
        )

        return SauterBootstrapCurrentParams(
            f_cur_plasma_bootstrap=np.sum(da * j_bootstrap_toroidal, axis=0)
            / plasma_current,
            j_bootstrap=j_bootstrap_toroidal,
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            normalised_ion_collisionality=normalised_ion_collisionality,
            l31_coefficient=l31_coefficient,
            l32_coefficient=l32_coefficient,
            l34_coefficient=l34_coefficient,
            alpha=alpha,
            alpha_0=alpha_0,
            q_profile=1 / inverse_q,
            inverse_aspect=inverse_aspect,
            sqrt_inverse_aspect=sqrt_inverse_aspect,
        )

    @staticmethod
    def _electron_collisionality_sauter(
        radial_elements: np.ndarray,
        rmajor: float,
        n_charge_plasma_effective: float | np.ndarray,
        inverse_q: np.ndarray,
        sqrt_inverse_aspect: float | np.ndarray,
        temp_plasma_electron_profile_kev: float | np.ndarray,
        nd_plasma_electron_profile: float | np.ndarray,
    ) -> np.ndarray:
        """Calculate normalized electron collisionality (νₑ⁎) using Sauter Eqs. 18b
        and 18d.

        Parameters
        ----------
        radial_elements:
            The radial element index in the range 1 to nr.
        rmajor:
            The plasma major radius (R) [m].
        n_charge_plasma_effective:
            Plasma effective charge (Z_eff) [-].
        inverse_q:
            Inverse safety factor profile [-].
        sqrt_inverse_aspect:
            The square root of the inverse aspect ratio array (√ε) [-].
        temp_plasma_electron_profile_kev:
            The electron temperature profile [keV].
        nd_plasma_electron_profile:
            The electron density profile [m⁻³].

        Returns
        -------
        np.ndarray
            The normalized electron collisionality (νₑ⁎).

        Raises
        ------
        ProcessValueError
            If an input to the selected points is non-positive or non-finite.

        Notes
        -----
        - The normalised electron collisionality (νₑ⁎) is calculated from the
        electron-electron coulomb logarithm given by Equation 18d in reference [1].
        This is then used in the full equation in Equation 18b.


        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        """  # noqa: RUF002
        selected = radial_elements - 1
        ne_si = nd_plasma_electron_profile[selected]
        temp_plasma_electron_profile_ev = (
            temp_plasma_electron_profile_kev[selected] * 1.0e3
        )
        inverse_q_selected = inverse_q[selected]
        if not all(
            np.all(np.isfinite(value) & (value > 0))
            for value in (
                ne_si,
                temp_plasma_electron_profile_ev,
                inverse_q_selected,
                n_charge_plasma_effective[selected],
                sqrt_inverse_aspect[selected],
                rmajor,
            )
        ):
            raise ProcessValueError(
                "Sauter electron collisionality inputs must be positive and finite."
            )
        # Electron-electron Coulomb logarithm given by Equation 18d in reference [1]
        log_lambda_e = 31.3 - np.log(np.sqrt(ne_si) / temp_plasma_electron_profile_ev)
        return (
            6.921e-18
            * rmajor
            * ne_si
            * n_charge_plasma_effective[selected]
            * log_lambda_e
            / (
                inverse_q_selected
                * temp_plasma_electron_profile_ev**2
                * sqrt_inverse_aspect[selected] ** 3
            )
        )

    @staticmethod
    def _ion_collisionality_sauter(
        radial_elements: np.ndarray,
        rmajor: float,
        inverse_q: np.ndarray,
        sqrt_inverse_aspect: np.ndarray,
        temp_plasma_ion_profile_kev: np.ndarray,
        n_charge_main_ion: np.ndarray,
        nd_plasma_ion_profile: np.ndarray,
    ) -> np.ndarray:
        """Calculate normalized ion collisionality (νᵢ⁎) using Sauter Eqs. 18c and 18e.

        Parameters
        ----------
        radial_elements:
            The radial element indexes in the range 1 to nr.
        rmajor:
            The plasma major radius (R₀) [m].
        inverse_q:
            Inverse safety factor profile [-].
        sqrt_inverse_aspect:
            The square root of the inverse aspect ratio profile [-].
        temp_plasma_ion_profile_kev:
            The ion temperature profile [keV].
        n_charge_main_ion : np.ndarray
            The charge of the main ion species [-].
        nd_plasma_ion_profile : np.ndarray
            The ion density profile [m⁻³].

        Returns
        -------
        np.ndarray
            The normalized ion collisionality (νᵢ⁎).

        Raises
        ------
        ProcessValueError
            If an input to the selected points is non-positive or non-finite.

        Notes
        -----
        - The normalised electron collisionality (νₑ⁎) is calculated from the
        electron-electron coulomb logarithm given by Equation 18e in reference [1].
        This is then used in the full equation in Equation 18c.


        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        """  # noqa: RUF002
        selected = radial_elements - 1
        ni_si = nd_plasma_ion_profile[selected]
        temp_plasma_ion_profile_ev = temp_plasma_ion_profile_kev[selected] * 1.0e3
        zion_selected = n_charge_main_ion[selected]
        inverse_q_selected = inverse_q[selected]
        if not all(
            np.all(np.isfinite(value) & (value > 0))
            for value in (
                ni_si,
                temp_plasma_ion_profile_ev,
                zion_selected,
                inverse_q_selected,
                sqrt_inverse_aspect[selected],
                rmajor,
            )
        ):
            raise ProcessValueError(
                "Sauter ion collisionality inputs must be positive and finite."
            )
        # Calculate the ion-ion Coulomb logarithm. Given by Equation 18e in
        # reference [1].
        log_lambda_ii = 30.0 - np.log(
            zion_selected**3 * np.sqrt(ni_si) / temp_plasma_ion_profile_ev**1.5
        )
        return (
            4.90e-18
            * rmajor
            * ni_si
            * zion_selected**4
            * log_lambda_ii
            / (
                inverse_q_selected
                * temp_plasma_ion_profile_ev**2
                * sqrt_inverse_aspect[selected] ** 3
            )
        )

    def _calculate_l31_coefficient(
        self,
        radial_elements: np.ndarray,
        f_trapped_particles: np.ndarray,
        normalised_electron_collisionality: np.ndarray,
        zeff: np.ndarray,
    ) -> np.ndarray:
        """Calculate the dimensionless L31 coefficient.

        Parameters
        ----------
        radial_elements:
            Radial element indexes.
        f_trapped_particles:
            Trapped particle fraction profile at the selected elements.
        normalised_electron_collisionality:
            Electron collisionality profile at the selected elements.
        zeff:
            Effective charge profile.

        Returns
        -------
        :
            Dimensionless L31 coefficient.
        """
        charge_profile = zeff[radial_elements - 1]

        # f³¹_{teff}(νₑ⁎), Eq.14b # noqa: RUF003
        f31_teff = self._calculate_f_teff_31_coefficient(
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            n_charge_plasma_effective=charge_profile,
        )

        # L₃₁, Eq.14a
        return self._calculate_f_31_factor(
            n_charge_plasma_effective=charge_profile,
            x=f31_teff,
        )

    def _calculate_l32_coefficient(
        self,
        radial_elements: np.ndarray,
        f_trapped_particles: np.ndarray,
        zeff: np.ndarray,
        normalised_electron_collisionality: np.ndarray,
    ) -> np.ndarray:
        """Calculate the dimensionless L32 coefficient.

        Parameters
        ----------
        radial_elements:
            Radial element indexes.
        f_trapped_particles:
            Trapped particle fraction profile at the selected elements.
        zeff:
            Effective charge profile.
        normalised_electron_collisionality:
            Electron collisionality profile at the selected elements.

        Returns
        -------
        :
            Dimensionless L32 coefficient.
        """
        charge_profile = zeff[radial_elements - 1]

        # f³²ᵉᵉ_{teff}(νₑ⁎), Eq.15d # noqa: RUF003
        f32ee_teff = self._calculate_f_teff_32ee_coefficient(
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            n_charge_plasma_effective=charge_profile,
        )

        # f³²ᵉⁱ_{teff}(νₑ⁎), Eq.15e # noqa: RUF003
        f32ei_teff = self._calculate_f_teff_32ei_coefficient(
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            n_charge_plasma_effective=charge_profile,
        )

        # F₃₂ₑₑ(X), Eq.15b
        big_f32ee_teff = self._calculate_f_32_ee_factor(
            n_charge_plasma_effective=charge_profile, x=f32ee_teff
        )

        # F₃₂ₑᵢ(Y), Eq.15c
        big_f32ei_teff = self._calculate_f_32_ei_factor(
            n_charge_plasma_effective=charge_profile, y=f32ei_teff
        )

        return big_f32ee_teff + big_f32ei_teff

    def _calculate_l34_coefficient(
        self,
        radial_elements: np.ndarray,
        f_trapped_particles: np.ndarray,
        normalised_electron_collisionality: np.ndarray,
        zeff: np.ndarray,
    ) -> np.ndarray:
        """Calculate the dimensionless L34 coefficient.

        Parameters
        ----------
        radial_elements:
            Radial element indexes.
        f_trapped_particles:
            Trapped particle fraction profile at the selected elements.
        normalised_electron_collisionality:
            Electron collisionality profile at the selected elements.
        zeff:
            Effective charge profile.

        Returns
        -------
        :
            Dimensionless L34 coefficient.
        """
        charge_profile = zeff[radial_elements - 1]

        # f³⁴_{teff}(νₑ⁎), Eq.16b, # noqa: RUF003
        f34_teff = self._calculate_f_teff_34_coefficient(
            f_trapped_particles=f_trapped_particles,
            normalised_electron_collisionality=normalised_electron_collisionality,
            n_charge_plasma_effective=charge_profile,
        )

        # L₃₄ Eq.16a
        return self._calculate_f_31_factor(
            n_charge_plasma_effective=charge_profile, x=f34_teff
        )

    @staticmethod
    def _calculate_f_teff_31_coefficient(
        f_trapped_particles: float | np.ndarray,
        normalised_electron_collisionality: float | np.ndarray,
        n_charge_plasma_effective: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the effective trapped particle fraction coefficient f³¹_{teff} as a
        function of electron collisionality and effective charge profile.

        Parameters
        ----------
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].
        normalised_electron_collisionality:
            Normalised electron collisionality (νₑ⁎) [-].
        n_charge_plasma_effective:
            Plasma effective charge (Z_eff) [-].

        Returns
        -------
        :
            The effective trapped particle fraction coefficient f³¹_{teff}.

        Notes
        -----
        - This equation is given by Eq.14b in reference [1]

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        """  # noqa: RUF002
        # f³¹_{teff}(νₑ⁎), Eq.14b # noqa: RUF003
        return f_trapped_particles / (
            (
                1.0
                + (1.0 - 0.1 * f_trapped_particles)
                * np.sqrt(normalised_electron_collisionality)
            )
            + (0.5 * (1.0 - f_trapped_particles) * normalised_electron_collisionality)
            / n_charge_plasma_effective
        )

    @staticmethod
    def _calculate_f_teff_32ee_coefficient(
        f_trapped_particles: float | np.ndarray,
        normalised_electron_collisionality: float | np.ndarray,
        n_charge_plasma_effective: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the effective trapped particle fraction coefficient f³²ᵉᵉ_{teff}
        as a function of normalised electron collisionality and effective charge
        profile.

        Parameters
        ----------
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].
        normalised_electron_collisionality:
            Normalised electron collisionality (νₑ⁎) [-].
        n_charge_plasma_effective:
            Plasma effective charge (Z_eff) [-].

        Returns
        -------
        :
            The effective trapped particle fraction coefficient f³²ᵉᵉ_{teff}(νₑ⁎).

        Notes
        -----
        - This equation is given by Eq.15d in reference [1]

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        """  # noqa: RUF002
        # f³²ᵉᵉ_{teff}(νₑ⁎), Eq.15d # noqa: RUF003
        return f_trapped_particles / (
            1.0
            + 0.26
            * (1.0 - f_trapped_particles)
            * np.sqrt(normalised_electron_collisionality)
            + (
                0.18
                * (1.0 - 0.37 * f_trapped_particles)
                * normalised_electron_collisionality
                / np.sqrt(n_charge_plasma_effective)
            )
        )

    @staticmethod
    def _calculate_f_teff_32ei_coefficient(
        f_trapped_particles: float | np.ndarray,
        normalised_electron_collisionality: float | np.ndarray,
        n_charge_plasma_effective: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the effective trapped particle fraction coefficient f³²ᵉⁱ_{teff}
        as a function of normalised electron collisionality and effective charge profile.

        Parameters
        ----------
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].
        normalised_electron_collisionality:
            Normalised electron collisionality (νₑ⁎) [-].
        n_charge_plasma_effective:
            Effective charge profile (Z_eff) [-].

        Returns
        -------
        :
            The effective trapped particle fraction coefficient f³²ᵉⁱ_{teff}(νₑ⁎).

        Notes
        -----
        - This equation is given by Eq.15e in reference [1]

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        """  # noqa: RUF002
        # f³²ᵉⁱ_{teff}(νₑ⁎), Eq.15e # noqa: RUF003
        return f_trapped_particles / (
            (
                1.0
                + (1.0 + 0.6 * f_trapped_particles)
                * np.sqrt(normalised_electron_collisionality)
            )
            + (
                0.85
                * (1.0 - 0.37 * f_trapped_particles)
                * normalised_electron_collisionality
                * (1.0 + n_charge_plasma_effective)
            )
        )

    @staticmethod
    def _calculate_f_teff_34_coefficient(
        f_trapped_particles: float | np.ndarray,
        normalised_electron_collisionality: float | np.ndarray,
        n_charge_plasma_effective: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the effective trapped particle fraction coefficient f³⁴_{teff}
        as a function of normalised electron collisionality and effective charge profile.

        Parameters
        ----------
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].
        normalised_electron_collisionality:
            Normalised electron collisionality (νₑ⁎) [-].
        n_charge_plasma_effective:
            Effective charge profile (Z_eff) [-].

        Returns
        -------
        :
            The effective trapped particle fraction coefficient f³⁴_{teff}(νₑ⁎).

        Notes
        -----
        - This equation is given by Eq.16b in reference [1]

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        """  # noqa: RUF002
        # f³⁴_{teff}(νₑ⁎), Eq.16b, # noqa: RUF003

        return f_trapped_particles / (
            (
                1.0
                + (1.0 - 0.1 * f_trapped_particles)
                * np.sqrt(normalised_electron_collisionality)
            )
            + 0.5
            * (1.0 - 0.5 * f_trapped_particles)
            * normalised_electron_collisionality
            / n_charge_plasma_effective
        )

    @staticmethod
    def _calculate_alpha_0_coefficient(
        f_trapped_particles: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the (α₀) coefficient. The collisionless (banana-regime) limit of
        the ion-temperature-gradient bootstrap coefficient.

        Parameters
        ----------
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].

        Returns
        -------
        :
            The (α₀) coefficient.

        Notes
        -----
        - This equation is given by Eq.17a in reference [1]

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        """
        # α₀, Eq.17a
        return (-1.17 * (1.0 - f_trapped_particles)) / (
            1.0 - (0.22 * f_trapped_particles) - (0.19 * f_trapped_particles**2)
        )

    @staticmethod
    def _calculate_alpha_coefficient(
        alpha_0: float | np.ndarray,
        f_trapped_particles: float | np.ndarray,
        normalised_ion_collisionality: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the (α) coefficient, the finite-collisionality correction to the
        ion-temperature-gradient contribution to the bootstrap current.

        Parameters
        ----------
        alpha_0:
            The (α₀) coefficient [-].
        f_trapped_particles:
            Trapped particle fraction (fₜ) [-].
        normalised_ion_collisionality:
            Normalised ion collisionality (νᵢ⁎) [-].

        Returns
        -------
        :
            The (α) coefficient.

        Notes
        -----
        - This equation is originally given by Eq.17b in reference [1], though an error
        concerning the sign of the 0.315 term is corrected to be positive in
        reference [2].

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        [2] O. Sauter, C. Angioni, Y. R. Lin-Liu; Erratum:
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime
        [Phys. Plasmas 6, 2834 (1999)]. Phys. Plasmas 1 December 2002; 9 (12): 5140.
        https://doi.org/10.1063/1.1517052
        """  # noqa: RUF002
        # α(νᵢ⁎), Eq.17b # noqa: RUF003
        return (
            (
                alpha_0
                + (0.25 * (1.0 - f_trapped_particles**2))
                * np.sqrt(normalised_ion_collisionality)
            )
            / (1.0 + (0.5 * np.sqrt(normalised_ion_collisionality)))
            + (0.315 * normalised_ion_collisionality**2 * f_trapped_particles**6)
        ) / (1.0 + (0.15 * normalised_ion_collisionality**2 * f_trapped_particles**6))

    @staticmethod
    def _calculate_f_32_ei_factor(
        n_charge_plasma_effective: float | np.ndarray,
        y: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the F₃₂ₑᵢ(Y) factor, the electron-ion friction correction factor
        for the electron temperature-gradient bootstrap coefficient L₃₂.

        Parameters
        ----------
        n_charge_plasma_effective:
            Effective plasma charge (Z_eff) [-].
        y:
            Normalized parameter Y [-].

        Returns
        -------
        :
            The F₃₂ₑᵢ(Y) factor.

        Notes
        -----
        - This equation is originally given by Eq.15c in reference [1].

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        """
        # F₃₂ₑᵢ(Y), Eq.15c
        return (
            (
                -(0.56 + 1.93 * n_charge_plasma_effective)
                / n_charge_plasma_effective
                / (1.0 + 0.44 * n_charge_plasma_effective)
                * (y - y**4)
            )
            + (
                4.95
                / (1.0 + 2.48 * n_charge_plasma_effective)
                * (y**2 - y**4 - 0.55 * (y**3 - y**4))
            )
            - (1.2 / (1.0 + 0.5 * n_charge_plasma_effective) * y**4)
        )

    @staticmethod
    def _calculate_f_32_ee_factor(
        n_charge_plasma_effective: float | np.ndarray, x: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the F₃₂ₑₑ(X) factor, the electron-electron friction correction
        factor for the electron temperature-gradient bootstrap coefficient L₃₂.

        Parameters
        ----------
        n_charge_plasma_effective:
            Effective plasma charge (Z_eff) [-].
        x:
            Normalized parameter X [-].

        Returns
        -------
        :
            The F₃₂ₑₑ(X) factor.

        Notes
        -----
        - This equation is originally given by Eq.15b in reference [1].

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        """
        # F₃₂ₑₑ(X), Eq.15b
        return (
            (
                (0.05 + 0.62 * n_charge_plasma_effective)
                / n_charge_plasma_effective
                / (1.0 + 0.44 * n_charge_plasma_effective)
                * (x - x**4)
            )
            + (
                (x**2 - x**4 - 1.2 * (x**3 - x**4))
                / (1.0 + 0.22 * n_charge_plasma_effective)
            )
            + (1.2 / (1.0 + 0.5 * n_charge_plasma_effective) * x**4)
        )

    @staticmethod
    def _calculate_f_31_factor(
        n_charge_plasma_effective: float | np.ndarray, x: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the F₃₁ factor. Used to find L₃₁ and L₃₄.

        Parameters
        ----------
        n_charge_plasma_effective:
            Effective plasma charge (Z_eff) [-].
        x:
            Normalized parameter X [-].

        Returns
        -------
        :
            The F₃₁ factor.

        Notes
        -----
        - This equation is originally given by Eq.14a in reference [1].

        References
        ----------
        [1] O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime.
        Phys. Plasmas 1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240
        """
        # F₃₁ Eq.14a
        return (
            ((1.0 + 1.4 / (n_charge_plasma_effective + 1.0)) * x)
            - (1.9 / (n_charge_plasma_effective + 1.0) * x**2)
            + ((0.3 * x**3 + 0.2 * x**4) / (n_charge_plasma_effective + 1.0))
        )

    @staticmethod
    @nb.njit(cache=True)
    def _trapped_particle_fraction_sauter(
        radial_elements: np.ndarray,
        triang: float,
        sqrt_inverse_aspect: np.ndarray,
        fit: int = 0,
    ) -> np.ndarray:
        """Calculates the trapped particle fraction to be used in the Sauter bootstrap
        current scaling.

        Parameters
        ----------
        radial_elements:
            Radial element index in range 1 to nr.
        triang:
            Plasma triangularity.
        sqrt_inverse_aspect:
            Square root of the inverse of the local aspect ratio.
        fit:
            Optional. Fit method (0 = ASTRA method, 1 = Sauter 2002, 2 = Sauter 2016).
            Default is 0.

        Returns
        -------
        :
            Trapped particle fraction.

        Raises
        ------
        ProcessValueError
            If fit is not 0, 1, or 2.

        Notes
        -----
        This function calculates the trapped particle fraction at a given radius.

        References
        ----------
        Used in this paper:

        O. Sauter, C. Angioni, Y. R. Lin-Liu;
        Neoclassical conductivity and bootstrap current formulas for general
        axisymmetric equilibria and arbitrary collisionality regime. Phys. Plasmas
        1 July 1999; 6 (7): 2834-2839. https://doi.org/10.1063/1.873240

        O. Sauter, R. J. Buttery, R. Felton, T. C. Hender, D. F. Howell, and
        contributors to the E.-J. Workprogramme, "Marginal-limit for neoclassical
        tearing modes in JET H-mode discharges," Plasma Physics and Controlled Fusion,
        vol. 44, no. 9, pp. 1999-2019, Aug. 2002,
        doi: https://doi.org/10.1088/0741-3335/44/9/315.

        O. Sauter, Geometric formulas for system codes including the effect of negative
        triangularity, Fusion Engineering and Design, Volume 112, 2016, Pages 633-645,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2016.04.033.
        """
        # Prevent first element from being zero
        sqeps_reduced = sqrt_inverse_aspect[radial_elements - 1]
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

        raise ProcessValueError(f"fit={fit} is not valid. Must be 0, 1, or 2")
