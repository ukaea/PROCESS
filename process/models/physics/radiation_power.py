"""Module for calculating radiation power densities in tokamak plasmas.

This module provides functions to compute synchrotron and impurity radiation
power densities for plasma systems.
"""

import logging
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

import process.models.physics.impurity_radiation as impurity
from process.core.data_structure.base import DataStructure
from process.models.physics.plasma_profiles import PlasmaProfile
import process.core.constants as constants
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData

logger = logging.getLogger(__name__)


@dataclass
class RadpwrData:
    """DataClass which holds the output of the function radpwr"""

    pden_plasma_sync_mw: float
    pden_plasma_core_rad_mw: float
    pden_plasma_outer_rad_mw: float
    pden_plasma_rad_mw: float


def calculate_radiation_powers(
    plasma_profile: PlasmaProfile,
    nd_plasma_electron_on_axis: float,
    rminor: float,
    b_plasma_toroidal_on_axis: float,
    aspect: float,
    alphan: float,
    alphat: float,
    tbeta: float,
    temp_plasma_electron_on_axis_kev: float,
    f_sync_reflect: float,
    rmajor: float,
    kappa: float,
    vol_plasma: float,
    data_structure: DataStructure,
) -> RadpwrData:
    """Calculate the radiation powers in MW/m^3 by calling relevant routines.

    This function computes the radiation power densities for the plasma, including
    impurity radiation and synchrotron radiation. It returns a dataclass containing
    the calculated radiation power densities.

    Parameters
    ----------
    plasma_profile : PlasmaProfile
        The parameterized temperature and density profiles of the plasma.
    nd_plasma_electron_on_axis : float
        Central electron density (m^-3).
    rminor : float
        Minor radius of the plasma (m).
    b_plasma_toroidal_on_axis : float
        Toroidal magnetic field (T).
    aspect : float
        Aspect ratio of the plasma.
    alphan : float
        Alpha parameter for density profile.
    alphat : float
        Alpha parameter for temperature profile.
    tbeta : float
        Beta parameter for temperature profile.
    temp_plasma_electron_on_axis_kev : float
        Central electron temperature (keV).
    f_sync_reflect : float
        Fraction of synchrotron radiation reflected.
    rmajor : float
        Major radius of the plasma (m).
    kappa : float
        Elongation of the plasma.
    vol_plasma : float
        Plasma volume (m^3).

    Returns
    -------
    RadpwrData
        A dataclass containing the following radiation power densities:
        - pden_plasma_sync_mw (float): Synchrotron radiation power per unit
          volume (MW/m^3).
        - pden_plasma_core_rad_mw (float): Total core radiation power per unit
          volume (MW/m^3).
        - pden_plasma_outer_rad_mw (float): Edge radiation power per unit
          volume (MW/m^3).
        - pden_plasma_rad_mw (float): Total radiation power per unit volume (MW/m^3).

    References
    ----------
        - F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron
          radiation losses in realistic tokamak plasmas,” Nuclear Fusion, vol. 41,
          no. 6, pp. 665-678, Jun. 2001,
          doi: https://doi.org/10.1088/0029-5515/41/6/301.

        - I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks
          of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758,
          Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
    """
    imp_rad = impurity.ImpurityRadiation(plasma_profile, data_structure)
    imp_rad.calculate_imprad()

    pden_plasma_outer_rad_mw = (
        imp_rad.pden_impurity_rad_total_mw - imp_rad.pden_impurity_core_rad_total_mw
    )

    # Synchrotron radiation power/volume; assumed to be from core only.
    pden_plasma_sync_mw = psync_albajar_fidone(
        nd_plasma_electron_on_axis,
        rminor,
        b_plasma_toroidal_on_axis,
        aspect,
        alphan,
        alphat,
        tbeta,
        temp_plasma_electron_on_axis_kev,
        f_sync_reflect,
        rmajor,
        kappa,
        vol_plasma,
    )

    # Total core radiation power/volume.
    pden_plasma_core_rad_mw = (
        imp_rad.pden_impurity_core_rad_total_mw + pden_plasma_sync_mw
    )

    # Total radiation power/volume.
    pden_plasma_rad_mw = imp_rad.pden_impurity_rad_total_mw + pden_plasma_sync_mw

    return RadpwrData(
        pden_plasma_sync_mw,
        pden_plasma_core_rad_mw,
        pden_plasma_outer_rad_mw,
        pden_plasma_rad_mw,
    )


def psync_albajar_fidone(
    nd_plasma_electron_on_axis: float,
    rminor: float,
    b_plasma_toroidal_on_axis: float,
    aspect: float,
    alphan: float,
    alphat: float,
    tbeta: float,
    temp_plasma_electron_on_axis_kev: float,
    f_sync_reflect: float,
    rmajor: float,
    kappa: float,
    vol_plasma: float,
) -> float:
    """Calculate the synchrotron radiation power in MW/m^3.

    This function computes the synchrotron radiation power density for the plasma based
    on the plasma shape, major and minor radii, electron density, and temperature
    profiles.

    Parameters
    ----------
    nd_plasma_electron_on_axis : float
        Central electron density (m^-3).
    rminor : float
        Minor radius of the plasma (m).
    b_plasma_toroidal_on_axis : float
        Toroidal magnetic field (T).
    aspect : float
        Aspect ratio of the plasma.
    alphan : float
        Alpha parameter for density profile.
    alphat : float
        Alpha parameter for temperature profile.
    tbeta : float
        Beta parameter for temperature profile.
    temp_plasma_electron_on_axis_kev : float
        Central electron temperature (keV).
    f_sync_reflect : float
        Fraction of synchrotron radiation reflected.
    rmajor : float
        Major radius of the plasma (m).
    kappa : float
        Elongation of the plasma.
    vol_plasma : float
        Plasma volume (m^3).

    Returns
    -------
    float
        Synchrotron radiation power per unit volume (MW/m^3).

    References
    ----------
    - F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron
      radiation losses in realistic tokamak plasmas,” Nuclear Fusion, vol. 41, no. 6,
      pp. 665-678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.

    - I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of
      arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758, Dec. 2001,
      doi: https://doi.org/10.1088/0029-5515/41/12/102.
    """
    # Variable names are created to closely match those from the reference papers.

    ne0_20 = 1.0e-20 * nd_plasma_electron_on_axis

    p_a0 = 6.04e3 * (rminor * ne0_20) / b_plasma_toroidal_on_axis

    g_function = 0.93 * (1.0 + 0.85 * np.exp(-0.82 * aspect))

    k_function = (
        (alphan + 3.87 * alphat + 1.46) ** -0.79
        * (1.98 + alphat) ** 1.36
        * tbeta**2.14
        * (tbeta**1.53 + 1.87 * alphat - 0.16) ** -1.33
    )

    dum = (
        1.0
        + 0.12
        * (temp_plasma_electron_on_axis_kev / p_a0**0.41)
        * (1.0 - f_sync_reflect) ** 0.41
    ) ** -1.51

    p_sync_mw = (
        3.84e-8
        * (1.0 - f_sync_reflect) ** 0.62
        * rmajor
        * rminor**1.38
        * kappa**0.79
        * b_plasma_toroidal_on_axis**2.62
        * ne0_20**0.38
        * temp_plasma_electron_on_axis_kev
        * (16.0 + temp_plasma_electron_on_axis_kev) ** 2.61
        * dum
        * g_function
        * k_function
    )

    # pden_plasma_sync_mw should be per unit volume; Albajar gives it as total

    return p_sync_mw / vol_plasma

def brem_gaunt_factor(te, z):
    mc2 = (
        constants.ELECTRON_MASS * constants.SPEED_LIGHT**2 
        / constants.KILOELECTRON_VOLT
    )
    t = te / mc2
    c0 = 2.0 * np.sqrt(3.0) / np.pi
    c_nr = [0.4302, 24.2255e-5, 0.7546e-5, 0.5282, 0.3301, 0.0911]
    c_r = [0.55467, 2.6346, -2.277595, 1.1480, -0.36465, 0.07451, -0.00975, 0.0007885, -3.5841e-5, 6.99834e-7]
    c_z = [5.760e4, 3.440, 16.80, 0.1333]
    f_nr = (
        c_nr[0] * (1.0 - np.exp(-(c_nr[1] * z**2 / t)**c_nr[3])) 
        - (c_nr[0] + c_nr[5]) * np.exp(-(t / (c_nr[2] * z**2))**c_nr[4])
    )
    f_r = 0.0
    for i in range(10):
        f_r += c0 * (c_r[i] * t**(i + 1))
    x = 100.0 * t * np.sqrt(10.0 / z)
    f_z = z * 1.0e-1 * c_z[0] * x**c_z[1] / (np.exp(c_z[2] * x**c_z[3]) - 1.0)
    largef_ee = 0.5 * (np.tanh(0.602 * (np.log10(t) + 5.06)) + 1.0)
    largef_nr = (
        3.0 * np.sqrt(3.0) / (np.sqrt(2.0) * np.pi) * t 
        * (np.tanh(-2.153 * np.log10(t / 0.43)) + 1.0)
    )
    if t <= 10.0:
        g_ee = (
            largef_ee * largef_nr 
            * (1.0 + 0.53 * t + 9.48 * t**2 - 0.67 * t**3 + 0.027 * t**4)
        )
        g_ei = c0 * (1.0 + f_nr - f_z) + f_r
    else:
        g_ee = (
            9.0 * np.sqrt(6.0) / (4.0 * np.sqrt(np.pi)) * np.sqrt(t) 
            * (np.log(2.0 * t) + 1.25 - 0.5772)
        )
        g_ei = (
            9.0 * np.sqrt(6.0) / (8.0 * np.sqrt(np.pi)) * np.sqrt(t) 
            * (np.log(2.0 * t) + 1.5 - 0.5772)
        )
    return g_ei, g_ee

def calculate_power_bremsstrahlung_radiation(
    f_density_h_isotope_electron: float,
    f_density_he_isotope_electron: float,
    ndensity_b11_fuel_vol_avg: float,
    impurity_radiation: ImpurityRadiationData,
    ne: float, 
    te: float,
    ):
    # Xie, 2024, Equation 61 and Section 5.
    # neglect line radiation
    c1 = (
        32.0 * np.pi * constants.ELECTRON_CHARGE**6 
        * ne**2 / (
            3.0 * (4.0 * np.pi * constants.EPSILON0)**3 
            * constants.PLANCK_CONSTANT * constants.ELECTRON_MASS * constants.SPEED_LIGHT**3
        )
    )
    c2 = np.sqrt(
        2.0 * np.pi * te * constants.KILOELECTRON_VOLT 
        / (3.0 * constants.ELECTRON_MASS)
    )
    g_h_isotope, g_ee = brem_gaunt_factor(te, 1.0) 
    g_h_isotope = g_h_isotope * 1.0 * f_density_h_isotope_electron
    g_he_isotope, _ = brem_gaunt_factor(te, 2.0) 
    g_he_isotope = g_he_isotope * 4.0 * f_density_he_isotope_electron
    g_b11, _ = brem_gaunt_factor(te, 5.0) 
    g_b11 = g_b11 * 25.0 * ndensity_b11_fuel_vol_avg / ne
    g_imp = 0.0
    for imp in range(15):
        z_imp = impurity_radiation.impurity_arr_z[imp]
        if z_imp > 2 and z_imp != 5:
            g_ei, _ = brem_gaunt_factor(te, z_imp) 
            g_ei = g_ei * z_imp**2 * impurity_radiation.f_nd_impurity_electron_array[imp]
            g_imp += g_ei
    g_total = g_h_isotope + g_he_isotope + g_b11 + g_imp
    pden = c1 * c2 * (g_total + g_ee) * 1.0e-6
    return pden

def calculate_power_radiation(
    nd_plasma_electron_on_axis: float,
    rminor: float,
    b_plasma_toroidal_on_axis: float,
    aspect: float,
    alphan: float,
    alphat: float,
    tbeta: float,
    temp_plasma_electron_on_axis_kev: float,
    f_sync_reflect: float,
    rmajor: float,
    kappa: float,
    vol_plasma: float,
    i_calculate_radiation: int,
    i_equilibrium_solve: int,
    f_density_h_isotope_electron: float,
    f_density_he_isotope_electron: float,
    ndensity_b11_fuel_vol_avg: float,
    impurity_radiation: ImpurityRadiationData,
    ne_vol_avg: float,
    temp_e_vol_avg: float,
    rho: np.array,
    ne_profile: np.array,
    te_profile: np.array,
    f_power_rad_brem_core_reduction: float,
    rho_plasma_core_norm: float,
    eq,
):

    # Synchrotron radiation power/volume; assumed to be from core only.
    pden_plasma_sync_mw = psync_albajar_fidone(
        nd_plasma_electron_on_axis,
        rminor,
        b_plasma_toroidal_on_axis,
        aspect,
        alphan,
        alphat,
        tbeta,
        temp_plasma_electron_on_axis_kev,
        f_sync_reflect,
        rmajor,
        kappa,
        vol_plasma,
    )

    if i_calculate_radiation == 1:
        # Xie, 2024, using average temperature and average density
        # assume all radiation power is from core
        # only consider bremsstrahlung radiation, no line radiation
        pden_rad_brem = calculate_power_bremsstrahlung_radiation(
            f_density_h_isotope_electron=f_density_h_isotope_electron,
            f_density_he_isotope_electron=f_density_he_isotope_electron,
            ndensity_b11_fuel_vol_avg=ndensity_b11_fuel_vol_avg,
            impurity_radiation=impurity_radiation,
            ne=ne_vol_avg, 
            te=temp_e_vol_avg,
        )
        pden_plasma_core_rad_mw = pden_rad_brem + pden_plasma_sync_mw
        pden_plasma_outer_rad_mw = 0.0
        pden_plasma_rad_mw = pden_rad_brem + pden_plasma_sync_mw
    
    elif i_calculate_radiation == 2:
        # Xie, 2024, simple integration/equilibrium integration
        # only consider bremsstrahlung radiation, no line radiation
        ne_fun = lambda rhox: np.interp(rhox, rho, ne_profile)
        te_fun = lambda rhox: np.interp(rhox, rho, te_profile)
        fun = lambda rhox: calculate_power_bremsstrahlung_radiation(
            f_density_h_isotope_electron=f_density_h_isotope_electron,
            f_density_he_isotope_electron=f_density_he_isotope_electron,
            ndensity_b11_fuel_vol_avg=ndensity_b11_fuel_vol_avg,
            impurity_radiation=impurity_radiation,
            ne=ne_fun(rhox), 
            te=te_fun(rhox),
        )        

        fun1 = lambda rhox: fun(rhox) * rhox
        if i_equilibrium_solve == 0:
            pden_rad_brem = quad(fun1, 0, 1)[0] * 2.0
            power_rad_brem_core = (
                2.0 * f_power_rad_brem_core_reduction 
                * quad(fun1, 0, rho_plasma_core_norm)[0] 
                * vol_plasma
            )
            power_rad_brem_edge = (
                2.0 
                * quad(fun1, rho_plasma_core_norm, 1)[0] 
                * vol_plasma
            )
        else:
            # pden(ρ) 在 (ρ,θ) 上；体积元 ∝ J R（与 get_average_quantity_2 一致）
            pden = np.array([fun(rhox) for rhox in eq.rho])[:, None]
            dV = eq.J * eq.R
            integrand = pden * dV
            vol = float(eq.grid.integrate(dV))
            pden_rad_brem = float(eq.grid.integrate(integrand)) / max(vol, 1.0e-30)

            # θ∈[0,2π]、ρ∈[0, ρ_core] 与 ρ∈[ρ_core,1] 分区积分
            rho_core = rho_plasma_core_norm
            mask_core = eq.rho <= rho_core
            integrand_core = np.zeros_like(integrand)
            integrand_edge = np.zeros_like(integrand)
            integrand_core[mask_core, :] = integrand[mask_core, :]
            integrand_edge[~mask_core, :] = integrand[~mask_core, :]
            power_rad_brem_core = (
                f_power_rad_brem_core_reduction
                * float(eq.grid.integrate(integrand_core))
            )
            power_rad_brem_edge = float(eq.grid.integrate(integrand_edge))
        
        pden_plasma_core_rad_mw = power_rad_brem_core / vol_plasma + pden_plasma_sync_mw
        pden_plasma_outer_rad_mw = power_rad_brem_edge / vol_plasma
        pden_plasma_rad_mw = pden_plasma_core_rad_mw + pden_plasma_outer_rad_mw

    return RadpwrData(
        pden_plasma_sync_mw,
        pden_plasma_core_rad_mw,
        pden_plasma_outer_rad_mw,
        pden_plasma_rad_mw,
    )
