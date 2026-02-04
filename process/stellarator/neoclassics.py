import logging

import numpy as np

from process import constants
from process.data_structure import (
    impurity_radiation_module,
    neoclassics_variables,
    physics_variables,
    stellarator_configuration,
    stellarator_variables,
)
from process.stellarator.stellarator import KEV

logger = logging.getLogger(__name__)


class Neoclassics:
    @property
    def no_roots(self):
        return neoclassics_variables.roots.shape[0]

    def init_neoclassics(self, r_effin, eps_effin, iotain):
        """Constructor of the neoclassics object from the effective radius,
        epsilon effective and iota only.
        """
        (
            neoclassics_variables.densities,
            neoclassics_variables.temperatures,
            neoclassics_variables.dr_densities,
            neoclassics_variables.dr_temperatures,
        ) = self.init_profile_values_from_PROCESS(r_effin)
        neoclassics_variables.roots = np.array([
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
        neoclassics_variables.weights = np.array([
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

        neoclassics_variables.kt = self.neoclassics_calc_KT()
        neoclassics_variables.nu = self.neoclassics_calc_nu()
        neoclassics_variables.nu_star = self.neoclassics_calc_nu_star()
        neoclassics_variables.nu_star_averaged = self.neoclassics_calc_nu_star_fromT(
            iotain
        )
        neoclassics_variables.vd = self.neoclassics_calc_vd()

        neoclassics_variables.d11_plateau = self.neoclassics_calc_D11_plateau()

        neoclassics_variables.d11_mono = self.neoclassics_calc_d11_mono(
            eps_effin
        )  # for using epseff

        neoclassics_variables.d111 = self.calc_integrated_radial_transport_coeffs(
            index=1
        )
        neoclassics_variables.d112 = self.calc_integrated_radial_transport_coeffs(
            index=2
        )
        neoclassics_variables.d113 = self.calc_integrated_radial_transport_coeffs(
            index=3
        )

        neoclassics_variables.gamma_flux = self.neoclassics_calc_gamma_flux(
            neoclassics_variables.densities,
            neoclassics_variables.temperatures,
            neoclassics_variables.dr_densities,
            neoclassics_variables.dr_temperatures,
        )
        neoclassics_variables.q_flux = self.neoclassics_calc_q_flux()

    def init_profile_values_from_PROCESS(self, rho):
        """Initializes the profile_values object from PROCESS' parabolic profiles"""
        tempe = (
            physics_variables.temp_plasma_electron_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempT = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempD = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempa = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )

        dense = (
            physics_variables.nd_plasma_electron_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densT = (
            (1 - physics_variables.f_plasma_fuel_deuterium)
            * physics_variables.nd_plasma_ions_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densD = (
            physics_variables.f_plasma_fuel_deuterium
            * physics_variables.nd_plasma_ions_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densa = (
            physics_variables.nd_plasma_alphas_vol_avg
            * (1 + physics_variables.alphan)
            * (1 - rho**2) ** physics_variables.alphan
        )

        # Derivatives in real space
        dr_tempe = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_electron_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )

        dr_dense = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.nd_plasma_electron_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * (1 - physics_variables.f_plasma_fuel_deuterium)
            * physics_variables.nd_plasma_ions_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.nd_plasma_ions_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.nd_plasma_alphas_vol_avg
            * (1 + physics_variables.alphan)
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )

        dens = np.array([dense, densD, densT, densa])
        temp = np.array([tempe, tempD, tempT, tempa])
        dr_dens = np.array([dr_dense, dr_densD, dr_densT, dr_densa])
        dr_temp = np.array([dr_tempe, dr_tempD, dr_tempT, dr_tempa])

        return dens, temp, dr_dens, dr_temp

    def calc_neoclassics(self):
        if stellarator_configuration.stella_config_epseff < 0:
            logger.error(
                f"epseff value lower than 0:  {stellarator_configuration.stella_config_epseff}"
            )
        self.init_neoclassics(
            0.6,
            stellarator_configuration.stella_config_epseff,
            stellarator_variables.iotabar,
        )

        q_PROCESS = (
            (
                physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        q_PROCESS_r1 = (
            (
                physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
        )

        q_neo = sum(neoclassics_variables.q_flux * 1e-6)
        gamma_neo = sum(
            neoclassics_variables.gamma_flux * neoclassics_variables.temperatures * 1e-6
        )

        total_q_neo = sum(
            neoclassics_variables.q_flux * 1e-6
            + neoclassics_variables.gamma_flux
            * neoclassics_variables.temperatures
            * 1e-6
        )

        total_q_neo_e = (
            2
            * 2
            * (
                neoclassics_variables.q_flux[0] * 1e-6
                + neoclassics_variables.gamma_flux[0]
                * neoclassics_variables.temperatures[0]
                * 1e-6
            )
        )

        q_neo_e = neoclassics_variables.q_flux[0] * 1e-6
        q_neo_D = neoclassics_variables.q_flux[1] * 1e-6
        q_neo_a = neoclassics_variables.q_flux[3] * 1e-6
        q_neo_T = neoclassics_variables.q_flux[2] * 1e-6

        g_neo_e = (
            neoclassics_variables.gamma_flux[0]
            * 1e-6
            * neoclassics_variables.temperatures[0]
        )
        g_neo_D = (
            neoclassics_variables.gamma_flux[1]
            * 1e-6
            * neoclassics_variables.temperatures[1]
        )
        g_neo_a = (
            neoclassics_variables.gamma_flux[3]
            * 1e-6
            * neoclassics_variables.temperatures[3]
        )
        g_neo_T = (
            neoclassics_variables.gamma_flux[2]
            * 1e-6
            * neoclassics_variables.temperatures[2]
        )

        dndt_neo_e = neoclassics_variables.gamma_flux[0]
        dndt_neo_D = neoclassics_variables.gamma_flux[1]
        dndt_neo_a = neoclassics_variables.gamma_flux[3]
        dndt_neo_T = neoclassics_variables.gamma_flux[2]

        dndt_neo_fuel = (
            (dndt_neo_D + dndt_neo_T)
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        dmdt_neo_fuel = (
            dndt_neo_fuel * physics_variables.m_fuel_amu * constants.PROTON_MASS * 1.0e6
        )  # mg
        dmdt_neo_fuel_from_e = (
            4
            * dndt_neo_e
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
            * physics_variables.m_fuel_amu
            * constants.PROTON_MASS
            * 1.0e6
        )  # kg

        chi_neo_e = -(
            neoclassics_variables.q_flux[0]
            + neoclassics_variables.gamma_flux[0] * neoclassics_variables.temperatures[0]
        ) / (
            neoclassics_variables.densities[0] * neoclassics_variables.dr_temperatures[0]
            + neoclassics_variables.temperatures[0]
            * neoclassics_variables.dr_densities[0]
        )

        chi_PROCESS_e = self.st_calc_eff_chi()

        nu_star_e = neoclassics_variables.nu_star_averaged[0]
        nu_star_d = neoclassics_variables.nu_star_averaged[1]
        nu_star_T = neoclassics_variables.nu_star_averaged[2]
        nu_star_He = neoclassics_variables.nu_star_averaged[3]

        return (
            q_PROCESS,
            q_PROCESS_r1,
            q_neo,
            gamma_neo,
            total_q_neo,
            total_q_neo_e,
            q_neo_e,
            q_neo_D,
            q_neo_a,
            q_neo_T,
            g_neo_e,
            g_neo_D,
            g_neo_a,
            g_neo_T,
            dndt_neo_e,
            dndt_neo_D,
            dndt_neo_a,
            dndt_neo_T,
            dndt_neo_fuel,
            dmdt_neo_fuel,
            dmdt_neo_fuel_from_e,
            chi_neo_e,
            chi_PROCESS_e,
            nu_star_e,
            nu_star_d,
            nu_star_T,
            nu_star_He,
        )

    def neoclassics_calc_KT(self):
        """Calculates the energy on the given grid
        which is given by the gauss laguerre roots.
        """
        k = np.repeat((neoclassics_variables.roots / KEV)[:, np.newaxis], 4, axis=1)

        return (k * neoclassics_variables.temperatures).T

    def neoclassics_calc_nu(self):
        """Calculates the collision frequency"""
        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(neoclassics_variables.densities[0])
            + 2.3
            * np.log10(neoclassics_variables.temperatures[0] / constants.ELECTRON_CHARGE)
        )

        neoclassics_calc_nu = np.zeros((4, self.no_roots), order="F")

        for j in range(4):
            for i in range(self.no_roots):
                x = neoclassics_variables.roots[i]
                for k in range(4):
                    xk = (
                        (mass[k] / mass[j])
                        * (
                            neoclassics_variables.temperatures[j]
                            / neoclassics_variables.temperatures[k]
                        )
                        * x
                    )
                    expxk = np.exp(-xk)
                    t = 1.0 / (1.0 + 0.3275911 * np.sqrt(xk))
                    erfn = (
                        1.0
                        - t
                        * (
                            0.254829592
                            + t
                            * (
                                -0.284496736
                                + t
                                * (1.421413741 + t * (-1.453152027 + t * 1.061405429))
                            )
                        )
                        * expxk
                    )
                    phixmgx = (1.0 - 0.5 / xk) * erfn + expxk / np.sqrt(np.pi * xk)
                    v = np.sqrt(
                        2.0 * x * neoclassics_variables.temperatures[j] / mass[j]
                    )
                    neoclassics_calc_nu[j, i] = neoclassics_calc_nu[
                        j, i
                    ] + neoclassics_variables.densities[k] * (
                        z[j] * z[k]
                    ) ** 2 * lnlambda * phixmgx / (
                        4.0 * np.pi * constants.EPSILON0**2 * mass[j] ** 2 * v**3
                    )

        return neoclassics_calc_nu

    def neoclassics_calc_nu_star(self):
        """Calculates the normalized collision frequency"""

        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[0, :] / (mass[0] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[1, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[1, :] / (mass[1] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[2, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[2, :] / (mass[2] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[3, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[3, :] / (mass[3] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )

        return (
            physics_variables.rmajor
            * neoclassics_variables.nu
            / (neoclassics_variables.iota * v)
        )

    def neoclassics_calc_nu_star_fromT(self, iota):
        """Calculates the collision frequency"""
        temp = (
            np.array([
                physics_variables.temp_plasma_electron_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
            ])
            * KEV
        )
        density = np.array([
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * physics_variables.f_plasma_fuel_deuterium,
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * (1 - physics_variables.f_plasma_fuel_deuterium),
            physics_variables.nd_plasma_alphas_vol_avg,
        ])

        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(density[0])
            + 2.3 * np.log10(temp[0] / constants.ELECTRON_CHARGE)
        )

        neoclassics_calc_nu_star_fromT = np.zeros((4,))

        for j in range(4):
            v = np.sqrt(2.0 * temp[j] / mass[j])
            for k in range(4):
                xk = (mass[k] / mass[j]) * (temp[j] / temp[k])

                expxk = 0.0
                if xk < 200.0:
                    expxk = np.exp(-xk)

                t = 1.0 / (1.0 + 0.3275911 * np.sqrt(xk))
                erfn = (
                    1.0
                    - t
                    * (
                        0.254829592
                        + t
                        * (
                            -0.284496736
                            + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))
                        )
                    )
                    * expxk
                )
                phixmgx = (1.0 - 0.5 / xk) * erfn + expxk / np.sqrt(np.pi * xk)
                neoclassics_calc_nu_star_fromT[j] = (
                    neoclassics_calc_nu_star_fromT[j]
                    + density[k]
                    * (z[j] * z[k]) ** 2
                    * lnlambda
                    * phixmgx
                    / (4.0 * np.pi * constants.EPSILON0**2 * mass[j] ** 2 * v**4)
                    * physics_variables.rmajor
                    / iota
                )
        return neoclassics_calc_nu_star_fromT

    def neoclassics_calc_vd(self):
        vde = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[0]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vdD = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[1]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vdT = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[2]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vda = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[3]
            / (
                2.0
                * constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )

        vd = np.empty((4, self.no_roots))

        vd[0, :] = vde
        vd[1, :] = vdD
        vd[2, :] = vdT
        vd[3, :] = vda

        return vd

    def neoclassics_calc_D11_plateau(self):
        """Calculates the plateau transport coefficients (D11_star sometimes)"""
        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[0, :] / (mass[0] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[1, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[1, :] / (mass[1] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[2, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[2, :] / (mass[2] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )
        v[3, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (neoclassics_variables.kt[3, :] / (mass[3] * constants.SPEED_LIGHT**2) + 1)
            ** (-1)
        )

        return (
            np.pi
            / 4.0
            * neoclassics_variables.vd**2
            * physics_variables.rmajor
            / neoclassics_variables.iota
            / v
        )

    def neoclassics_calc_d11_mono(self, eps_eff):
        """Calculates the monoenergetic radial transport coefficients
        using epsilon effective
        """
        return (
            4.0
            / (9.0 * np.pi)
            * (2.0 * eps_eff) ** (3.0 / 2.0)
            * neoclassics_variables.vd**2
            / neoclassics_variables.nu
        )

    def calc_integrated_radial_transport_coeffs(self, index: int):
        """Calculates the integrated radial transport coefficients (index `index`)
        It uses Gauss laguerre integration
        https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        """
        return np.sum(
            2.0
            / np.sqrt(np.pi)
            * neoclassics_variables.d11_mono
            * neoclassics_variables.roots ** (index - 0.5)
            * neoclassics_variables.weights,
            axis=1,
        )

    def neoclassics_calc_gamma_flux(
        self, densities, temperatures, dr_densities, dr_temperatures
    ):
        """Calculates the Energy flux by neoclassical particle transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -densities
            * neoclassics_variables.d111
            * (
                (dr_densities / densities - z * neoclassics_variables.er / temperatures)
                + (neoclassics_variables.d112 / neoclassics_variables.d111 - 3.0 / 2.0)
                * dr_temperatures
                / temperatures
            )
        )

    def neoclassics_calc_q_flux(self):
        """Calculates the Energy flux by neoclassicsal energy transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -neoclassics_variables.densities
            * neoclassics_variables.temperatures
            * neoclassics_variables.d112
            * (
                (
                    neoclassics_variables.dr_densities / neoclassics_variables.densities
                    - z * neoclassics_variables.er / neoclassics_variables.temperatures
                )
                + (neoclassics_variables.d113 / neoclassics_variables.d112 - 3.0 / 2.0)
                * neoclassics_variables.dr_temperatures
                / neoclassics_variables.temperatures
            )
        )

    def st_calc_eff_chi(self):
        volscaling = (
            physics_variables.vol_plasma
            * stellarator_variables.f_st_rmajor
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
            ** 2
        )
        surfacescaling = (
            physics_variables.a_plasma_surface
            * stellarator_variables.f_st_rmajor
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
        )

        nominator = (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.pden_alpha_total_mw
            - physics_variables.pden_plasma_core_rad_mw
        ) * volscaling

        # in fortran there was a 0*alphan term which I have removed for obvious reasons
        # the following comment seems to describe this?
        # "include alphan if chi should be incorporate density gradients too"
        # but the history can be consulted if required (23/11/22 TN)
        denominator = (
            (
                3
                * physics_variables.nd_plasma_electron_on_axis
                * constants.ELECTRON_CHARGE
                * physics_variables.temp_plasma_electron_on_axis_kev
                * 1e3
                * physics_variables.alphat
                * impurity_radiation_module.radius_plasma_core_norm
                * (1 - impurity_radiation_module.radius_plasma_core_norm**2)
                ** (physics_variables.alphan + physics_variables.alphat - 1)
            )
            * surfacescaling
            * 1e-6
        )

        return nominator / denominator
