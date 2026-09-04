"""Module containing neoclassics routines"""

import logging
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.polynomial import polyval

from process.core import constants
from process.core.model import Model
from process.models.stellarator.stellarator import KEV

logger = logging.getLogger(__name__)


@dataclass
class NormalisedCollisionality:
    """Storage for normalised collisionality"""

    e: float
    """electron"""
    D: float
    """Deutrerium"""
    T: float
    """Tritium"""
    He: float
    """Helium"""


class Neoclassics(Model):
    """Module containing neoclassics routines"""

    @property
    def no_roots(self):
        """Obtain number of Gauss Laguerre roots"""
        return self.data.neoclassics.roots.shape[0]

    @property
    def mass(self):
        """Component mass array for electron, deuterium, tritium and Helium"""
        return np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])

    def output(self):
        """Neoclassics model doesn't have any output"""

    def run(self):
        """Neoclassics model doesn't need to be run"""

    def init_neoclassics(self, r_effin, eps_effin, iotain):
        """Constructor of the neoclassics object from the effective radius,
        epsilon effective and iota only.

        Parameters
        ----------
        r_effin :

        eps_effin :

        iotain :

        """
        (
            self.data.neoclassics.densities,
            self.data.neoclassics.temperatures,
            self.data.neoclassics.dr_densities,
            self.data.neoclassics.dr_temperatures,
        ) = self.init_profile_values_from_PROCESS(r_effin)
        self.data.neoclassics.roots = np.array([
            4.740718054080526184e-2,
            2.499239167531593919e-1,
            6.148334543927683749e-1,
            1.143195825666101451,
            1.836454554622572344,
            2.696521874557216147,
            3.725814507779509288,
            4.927293765849881879,
            6.304515590965073635,
            7.861693293370260349,
            9.603775985479263255,
            1.153654659795613924e1,
            1.366674469306423489e1,
            1.600222118898106771e1,
            1.855213484014315029e1,
            2.132720432178312819e1,
            2.434003576453269346e1,
            2.760555479678096091e1,
            3.114158670111123683e1,
            3.496965200824907072e1,
            3.911608494906788991e1,
            4.361365290848483056e1,
            4.850398616380419980e1,
            5.384138540650750571e1,
            5.969912185923549686e1,
            6.618061779443848991e1,
            7.344123859555988076e1,
            8.173681050672767867e1,
            9.155646652253683726e1,
            1.041575244310588886e2,
        ])
        self.data.neoclassics.weights = np.array([
            1.160440860204388913e-1,
            2.208511247506771413e-1,
            2.413998275878537214e-1,
            1.946367684464170855e-1,
            1.237284159668764899e-1,
            6.367878036898660943e-2,
            2.686047527337972682e-2,
            9.338070881603925677e-3,
            2.680696891336819664e-3,
            6.351291219408556439e-4,
            1.239074599068830081e-4,
            1.982878843895233056e-5,
            2.589350929131392509e-6,
            2.740942840536013206e-7,
            2.332831165025738197e-8,
            1.580745574778327984e-9,
            8.427479123056716393e-11,
            3.485161234907855443e-12,
            1.099018059753451500e-13,
            2.588312664959080167e-15,
            4.437838059840028968e-17,
            5.365918308212045344e-19,
            4.393946892291604451e-21,
            2.311409794388543236e-23,
            7.274588498292248063e-26,
            1.239149701448267877e-28,
            9.832375083105887477e-32,
            2.842323553402700938e-35,
            1.878608031749515392e-39,
            8.745980440465011553e-45,
        ])

        self.data.neoclassics.kt = self.neoclassics_calc_KT()
        self.data.neoclassics.nu = self.neoclassics_calc_nu()
        self.data.neoclassics.nu_star = self.neoclassics_calc_nu_star()
        self.data.neoclassics.nu_star_averaged = self.neoclassics_calc_nu_star_fromT(
            iotain
        )
        self.data.neoclassics.vd = self.neoclassics_calc_vd()

        self.data.neoclassics.d11_plateau = self.neoclassics_calc_D11_plateau()

        self.data.neoclassics.d11_mono = self.neoclassics_calc_d11_mono(
            eps_effin
        )  # for using epseff

        self.data.neoclassics.d111 = self.calc_integrated_radial_transport_coeffs(
            index=1
        )
        self.data.neoclassics.d112 = self.calc_integrated_radial_transport_coeffs(
            index=2
        )
        self.data.neoclassics.d113 = self.calc_integrated_radial_transport_coeffs(
            index=3
        )

        self.data.neoclassics.gamma_flux = self.neoclassics_calc_gamma_flux(
            self.data.neoclassics.densities,
            self.data.neoclassics.temperatures,
            self.data.neoclassics.dr_densities,
            self.data.neoclassics.dr_temperatures,
        )
        self.data.neoclassics.q_flux = self.neoclassics_calc_q_flux()

    def init_profile_values_from_PROCESS(self, rho) -> tuple[np.ndarray, ...]:
        """Initialises the profile_values object from PROCESS' parabolic profiles

        Parameters
        ----------
        rho :

        Returns
        -------
        dens:
            density of electron, deuterium, tritum and alpha
        temp:
            temperature of electron, deuterium, tritum and alpha
        dr_dens:
            derivative density of electron, deuterium, tritum and alpha
        dr_temp:
            derivative temperature of electron, deuterium, tritum and alpha
        """
        t_suffix = (1 - rho**2) ** self.data.physics.alphat * KEV

        dens_suffix = (1 - rho**2) ** self.data.physics.alphan

        deriv_prefix = -2.0 * 1.0 / self.data.physics.rminor * rho

        # array are set up as electron, deuterium, tritum, alpha
        density_array = np.array([
            self.data.physics.nd_plasma_electron_on_axis,
            self.data.physics.f_plasma_fuel_deuterium
            * self.data.physics.nd_plasma_ions_on_axis,
            (1 - self.data.physics.f_plasma_fuel_deuterium)
            * self.data.physics.nd_plasma_ions_on_axis,
            self.data.physics.nd_plasma_alphas_thermal_vol_avg
            * (1 + self.data.physics.alphan),
        ])

        temperature_array = np.array([
            self.data.physics.temp_plasma_electron_on_axis_kev,
            self.data.physics.temp_plasma_ion_on_axis_kev,
            self.data.physics.temp_plasma_ion_on_axis_kev,
            self.data.physics.temp_plasma_ion_on_axis_kev,
        ])

        # Derivatives in real space
        dr_dens = (
            deriv_prefix
            * density_array
            * (1.0 - rho**2) ** (self.data.physics.alphan - 1.0)
            * self.data.physics.alphan
        )
        dr_temp = (
            deriv_prefix
            * temperature_array
            * (1.0 - rho**2) ** (self.data.physics.alphat - 1.0)
            * self.data.physics.alphat
            * KEV
        )
        return (
            (density_array * dens_suffix),
            (temperature_array * t_suffix),
            dr_dens,
            dr_temp,
        )

    def calc_neoclassics(self):
        """Calculate neoclassics parameters"""
        if self.data.stellarator_config.stella_config_epseff < 0:
            logger.error(
                "epseff value lower than 0: "
                f"{self.data.stellarator_config.stella_config_epseff}"
            )
        self.init_neoclassics(
            0.6,
            self.data.stellarator_config.stella_config_epseff,
            self.data.stellarator.iotabar,
        )

        q_PROCESS_r1 = (
            (
                self.data.physics.f_p_alpha_plasma_deposited
                * self.data.physics.pden_alpha_total_mw
                - self.data.physics.pden_plasma_core_rad_mw
            )
            * self.data.physics.vol_plasma
            / self.data.physics.a_plasma_surface
        )

        q_PROCESS = q_PROCESS_r1 * self.data.impurity_radiation.radius_plasma_core_norm

        dndt_neo_e = self.data.neoclassics.gamma_flux[0]
        q_neo_e = self.data.neoclassics.q_flux[0] * 1e-6
        g_neo_e = 1e-6 * dndt_neo_e * self.data.neoclassics.temperatures[0]
        total_q_neo_e = 4 * (q_neo_e + g_neo_e)

        dmdt_neo_fuel_from_e = (
            4e6
            * dndt_neo_e
            * self.data.physics.a_plasma_surface
            * self.data.impurity_radiation.radius_plasma_core_norm
            * self.data.physics.m_fuel_amu
            * constants.PROTON_MASS
        )  # kg

        chi_neo_e = -(
            self.data.neoclassics.q_flux[0]
            + dndt_neo_e * self.data.neoclassics.temperatures[0]
        ) / (
            self.data.neoclassics.densities[0] * self.data.neoclassics.dr_temperatures[0]
            + self.data.neoclassics.temperatures[0]
            * self.data.neoclassics.dr_densities[0]
        )

        chi_PROCESS_e = self.st_calc_eff_chi()

        # Unused calculations
        # q_neo_sum = 1e-6 * sum(self.data.neoclassics.q_flux)
        # gamma_neo = 1e-6 * sum(
        #     self.data.neoclassics.gamma_flux * self.data.neoclassics.temperatures
        # )

        # total_q_neo = 1e-6 * sum(
        #     self.data.neoclassics.q_flux
        #     + self.data.neoclassics.gamma_flux * self.data.neoclassics.temperatures
        # )

        # dndt_neo_D = self.data.neoclassics.gamma_flux[1]
        # dndt_neo_a = self.data.neoclassics.gamma_flux[3]
        # dndt_neo_T = self.data.neoclassics.gamma_flux[2]

        # dndt_neo_fuel = (
        #     (dndt_neo_D + dndt_neo_T)
        #     * self.data.physics.a_plasma_surface
        #     * self.data.impurity_radiation.radius_plasma_core_norm
        # )
        # dmdt_neo_fuel = (
        #     dndt_neo_fuel
        #     * self.data.physics.m_fuel_amu
        #     * constants.PROTON_MASS
        #     * 1.0e6
        # )  # mg
        return (
            q_PROCESS,
            q_PROCESS_r1,
            total_q_neo_e,
            q_neo_e,
            g_neo_e,
            dndt_neo_e,
            dmdt_neo_fuel_from_e,
            chi_neo_e,
            chi_PROCESS_e,
            NormalisedCollisionality(
                e=self.data.neoclassics.nu_star_averaged[0],
                D=self.data.neoclassics.nu_star_averaged[1],
                T=self.data.neoclassics.nu_star_averaged[2],
                He=self.data.neoclassics.nu_star_averaged[3],
            ),
        )

    def neoclassics_calc_KT(self):
        """Calculates the energy on the given grid
        which is given by the gauss laguerre roots.
        """
        return (
            self.data.neoclassics.roots[None] / KEV
        ) * self.data.neoclassics.temperatures[:, None]

    @staticmethod
    def _nu_erfn(xk, expxk):
        """Error function"""
        # Rational approximation for erf
        # t term
        t = 1.0 / (1.0 + 0.3275911 * np.sqrt(xk))
        # Expand: t * (c0 + c1*t + c2*t^2 + ...)
        coeffs = np.array([
            0,
            0.254829592,
            -0.284496736,
            1.421413741,
            -1.453152027,
            1.061405429,
        ])
        return 1.0 - polyval(t, coeffs) * expxk

    def neoclassics_calc_nu(self):
        """Calculates the collision frequency"""
        mass = self.mass
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.
        # Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(self.data.neoclassics.densities[0])
            + 2.3
            * np.log10(self.data.neoclassics.temperatures[0] / constants.ELECTRON_CHARGE)
        )

        roots = self.data.neoclassics.roots
        temp = self.data.neoclassics.temperatures

        inv_mass = 1 / mass
        xk = np.einsum("k,j,j,k,r->jkr", mass, inv_mass, temp, 1 / temp, roots)

        expxk = np.exp(-xk)

        erfn = self._nu_erfn(xk, expxk)

        phixmgx = (1.0 - 0.5 * 1 / xk) * erfn + expxk / np.sqrt(np.pi * xk)

        v = np.sqrt(2.0 * np.einsum("r,j,j->jr", roots, temp, inv_mass))
        denom = 1 / (4 * np.pi * constants.EPSILON0**2 * mass[:, None] ** 2 * v**3)

        return lnlambda * np.einsum(
            "k,jk,jkr,jr->jr",
            self.data.neoclassics.densities,
            np.einsum("i,j->ij", z, z) ** 2,
            phixmgx,
            denom,
        )

    def neoclassics_calc_nu_star(self):
        """Calculates the normalized collision frequency"""
        kk = (
            self.data.neoclassics.roots[None]
            * self.data.neoclassics.temperatures[:, None]
        )

        mass = self.mass

        v = constants.SPEED_LIGHT * np.sqrt(
            1.0 - (kk / (mass[:, None] * constants.SPEED_LIGHT**2) + 1) ** -1
        )
        return (
            self.data.physics.rmajor
            * self.data.neoclassics.nu
            / (self.data.neoclassics.iota * v)
        )

    def neoclassics_calc_nu_star_fromT(self, iota):
        """Calculates the collision frequency

        Parameters
        ----------
        iota :


        """
        temp = KEV * np.array([
            self.data.physics.temp_plasma_electron_vol_avg_kev,
            self.data.physics.temp_plasma_ion_vol_avg_kev,
            self.data.physics.temp_plasma_ion_vol_avg_kev,
            self.data.physics.temp_plasma_ion_vol_avg_kev,
        ])

        density = np.array([
            self.data.physics.nd_plasma_electrons_vol_avg,
            self.data.physics.nd_plasma_fuel_ions_vol_avg
            * self.data.physics.f_plasma_fuel_deuterium,
            self.data.physics.nd_plasma_fuel_ions_vol_avg
            * (1 - self.data.physics.f_plasma_fuel_deuterium),
            self.data.physics.nd_plasma_alphas_thermal_vol_avg,
        ])

        mass = self.mass
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.
        # Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(density[0])
            + 2.3 * np.log10(temp[0] / constants.ELECTRON_CHARGE)
        )

        inv_mass = 1 / mass
        v = np.sqrt(2.0 * temp * inv_mass)
        xk = np.einsum("k,j,j,k->jk", mass, inv_mass, temp, 1 / temp)

        # exp(-xk), clipped at xk >= 200
        mask_xk_lt_200 = xk < 200.0
        expxk = np.zeros_like(xk)
        expxk[mask_xk_lt_200] = np.exp(-xk[mask_xk_lt_200])

        erfn = self._nu_erfn(xk, expxk)

        phixmgx = (1.0 - 0.5 * 1 / xk) * erfn + expxk / np.sqrt(np.pi * xk)

        denom = 1 / (4 * np.pi * constants.EPSILON0**2 * mass**2 * v**4)

        # sum over k dimension
        return (
            lnlambda
            * self.data.physics.rmajor
            / iota
            * np.einsum(
                "k,jk,jk,j->j", density, np.einsum("i,j->ij", z, z) ** 2, phixmgx, denom
            )
        )

    def neoclassics_calc_vd(self):
        """Calculates the drift velocity on GL roots"""
        # alpha denominator is 2*vd_suffix
        vd_suffix = (
            constants.ELECTRON_CHARGE
            * self.data.physics.rmajor
            * self.data.physics.b_plasma_toroidal_on_axis
            * np.array([1, 1, 1, 2])
        )
        return (
            self.data.neoclassics.roots[None]
            * self.data.neoclassics.temperatures[:, None]
            / vd_suffix[:, None]
        )

    def neoclassics_calc_D11_plateau(self):
        """Calculates the plateau transport coefficients (D11_star sometimes)"""
        mass = self.mass

        v = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (
                self.data.neoclassics.kt[None]
                / (mass[:, None] * constants.SPEED_LIGHT**2)
                + 1
            )
            ** (-1)
        )

        return (
            np.pi
            / 4.0
            * self.data.neoclassics.vd**2
            * self.data.physics.rmajor
            / self.data.neoclassics.iota
            / v
        )

    def neoclassics_calc_d11_mono(self, eps_eff):
        """Calculates the monoenergetic radial transport coefficients
        using epsilon effective

        Parameters
        ----------
        eps_eff :

        """
        return (
            4.0
            / (9.0 * np.pi)
            * (2.0 * eps_eff) ** (3.0 / 2.0)
            * self.data.neoclassics.vd**2
            / self.data.neoclassics.nu
        )

    def calc_integrated_radial_transport_coeffs(self, index: int):
        """Calculates the integrated radial transport coefficients (index `index`)
        It uses Gauss laguerre integration
        https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature

        Parameters
        ----------
        index:

        """
        return np.sum(
            2.0
            / np.sqrt(np.pi)
            * self.data.neoclassics.d11_mono
            * self.data.neoclassics.roots ** (index - 0.5)
            * self.data.neoclassics.weights,
            axis=1,
        )

    def neoclassics_calc_gamma_flux(
        self, densities, temperatures, dr_densities, dr_temperatures
    ):
        """Calculates the Energy flux by neoclassical particle transport

        Parameters
        ----------
        densities :

        temperatures :

        dr_densities :

        dr_temperatures :

        """
        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -densities
            * self.data.neoclassics.d111
            * (
                (dr_densities / densities - z * self.data.neoclassics.er / temperatures)
                + (self.data.neoclassics.d112 / self.data.neoclassics.d111 - 3.0 / 2.0)
                * dr_temperatures
                / temperatures
            )
        )

    def neoclassics_calc_q_flux(self):
        """Calculates the Energy flux by neoclassicsal energy transport"""
        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -self.data.neoclassics.densities
            * self.data.neoclassics.temperatures
            * self.data.neoclassics.d112
            * (
                (
                    self.data.neoclassics.dr_densities / self.data.neoclassics.densities
                    - z * self.data.neoclassics.er / self.data.neoclassics.temperatures
                )
                + (self.data.neoclassics.d113 / self.data.neoclassics.d112 - 3.0 / 2.0)
                * self.data.neoclassics.dr_temperatures
                / self.data.neoclassics.temperatures
            )
        )

    def st_calc_eff_chi(self):
        """Calculates the effective thermal diffusivity from the alpha
        heating in the core to the boundary
        """
        volscaling = (
            self.data.physics.vol_plasma
            * self.data.stellarator.f_st_rmajor
            * (
                self.data.impurity_radiation.radius_plasma_core_norm
                * self.data.physics.rminor
                / self.data.stellarator_config.stella_config_rminor_ref
            )
            ** 2
        )
        surfacescaling = (
            self.data.physics.a_plasma_surface
            * self.data.stellarator.f_st_rmajor
            * (
                self.data.impurity_radiation.radius_plasma_core_norm
                * self.data.physics.rminor
                / self.data.stellarator_config.stella_config_rminor_ref
            )
        )

        nominator = (
            self.data.physics.f_p_alpha_plasma_deposited
            * self.data.physics.pden_alpha_total_mw
            - self.data.physics.pden_plasma_core_rad_mw
        ) * volscaling

        # in fortran there was a 0*alphan term which I have removed for obvious reasons
        # the following comment seems to describe this?
        # "include alphan if chi should be incorporate density gradients too"
        # but the history can be consulted if required (23/11/22 TN)
        denominator = (
            (
                3
                * self.data.physics.nd_plasma_electron_on_axis
                * constants.ELECTRON_CHARGE
                * self.data.physics.temp_plasma_electron_on_axis_kev
                * 1e3
                * self.data.physics.alphat
                * self.data.impurity_radiation.radius_plasma_core_norm
                * (1 - self.data.impurity_radiation.radius_plasma_core_norm**2)
                ** (self.data.physics.alphan + self.data.physics.alphat - 1)
            )
            * surfacescaling
            * 1e-6
        )

        return nominator / denominator
