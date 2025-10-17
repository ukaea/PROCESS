import logging
from dataclasses import dataclass

import numpy as np

import process.impurity_radiation as impurity
from process import constants
from process.data_structure import physics_variables
from process.plasma_profiles import PlasmaProfile

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
) -> RadpwrData:
    """
    Calculate the radiation powers in MW/m^3 by calling relevant routines.

    This function computes the radiation power densities for the plasma, including
    impurity radiation and synchrotron radiation. It returns a dataclass containing
    the calculated radiation power densities.

    :param plasma_profile: The parameterized temperature and density profiles of the plasma.
    :type plasma_profile: PlasmaProfile
    :param nd_plasma_electron_on_axis: Central electron density (m^-3).
    :type nd_plasma_electron_on_axis: float
    :param rminor: Minor radius of the plasma (m).
    :type rminor: float
    :param b_plasma_toroidal_on_axis: Toroidal magnetic field (T).
    :type b_plasma_toroidal_on_axis: float
    :param aspect: Aspect ratio of the plasma.
    :type aspect: float
    :param alphan: Alpha parameter for density profile.
    :type alphan: float
    :param alphat: Alpha parameter for temperature profile.
    :type alphat: float
    :param tbeta: Beta parameter for temperature profile.
    :type tbeta: float
    :param temp_plasma_electron_on_axis_kev: Central electron temperature (keV).
    :type temp_plasma_electron_on_axis_kev: float
    :param f_sync_reflect: Fraction of synchrotron radiation reflected.
    :type f_sync_reflect: float
    :param rmajor: Major radius of the plasma (m).
    :type rmajor: float
    :param kappa: Elongation of the plasma.
    :type kappa: float
    :param vol_plasma: Plasma volume (m^3).
    :type vol_plasma: float

    :returns: A dataclass containing the following radiation power densities:
        - pden_plasma_sync_mw (float): Synchrotron radiation power per unit volume (MW/m^3).
        - pden_plasma_core_rad_mw (float): Total core radiation power per unit volume (MW/m^3).
        - pden_plasma_outer_rad_mw (float): Edge radiation power per unit volume (MW/m^3).
        - pden_plasma_rad_mw (float): Total radiation power per unit volume (MW/m^3).
    :rtype: RadpwrData

    :references:
        - F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron radiation losses in realistic tokamak plasmas,”
          Nuclear Fusion, vol. 41, no. 6, pp. 665-678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.
        - I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of arbitrary geometry,”
          Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758, Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.

    """
    imp_rad = impurity.ImpurityRadiation(plasma_profile)
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
    """
    Calculate the synchrotron radiation power in MW/m^3.

    This function computes the synchrotron radiation power density for the plasma based on
    the plasma shape, major and minor radii, electron density, and temperature profiles.

    :param nd_plasma_electron_on_axis: Central electron density (m^-3).
    :type nd_plasma_electron_on_axis: float
    :param rminor: Minor radius of the plasma (m).
    :type rminor: float
    :param b_plasma_toroidal_on_axis: Toroidal magnetic field (T).
    :type b_plasma_toroidal_on_axis: float
    :param aspect: Aspect ratio of the plasma.
    :type aspect: float
    :param alphan: Alpha parameter for density profile.
    :type alphan: float
    :param alphat: Alpha parameter for temperature profile.
    :type alphat: float
    :param tbeta: Beta parameter for temperature profile.
    :type tbeta: float
    :param temp_plasma_electron_on_axis_kev: Central electron temperature (keV).
    :type temp_plasma_electron_on_axis_kev: float
    :param f_sync_reflect: Fraction of synchrotron radiation reflected.
    :type f_sync_reflect: float
    :param rmajor: Major radius of the plasma (m).
    :type rmajor: float
    :param kappa: Elongation of the plasma.
    :type kappa: float
    :param vol_plasma: Plasma volume (m^3).
    :type vol_plasma: float

    :returns: Synchrotron radiation power per unit volume (MW/m^3).
    :rtype: float

    :references:
        - F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron radiation losses in realistic tokamak plasmas,”
          Nuclear Fusion, vol. 41, no. 6, pp. 665-678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.

        - I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of arbitrary geometry,”
          Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758, Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
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


def fast_alpha_beta(
    b_plasma_poloidal_average: float,
    b_plasma_toroidal_on_axis: float,
    nd_plasma_electrons_vol_avg: float,
    nd_plasma_fuel_ions_vol_avg: float,
    nd_plasma_ions_total_vol_avg: float,
    temp_plasma_electron_density_weighted_kev: float,
    temp_plasma_ion_density_weighted_kev: float,
    pden_alpha_total_mw: float,
    pden_plasma_alpha_mw: float,
    i_beta_fast_alpha: int,
) -> float:
    """
    Calculate the fast alpha beta component.

    This function computes the fast alpha beta contribution based on the provided plasma parameters.

    Parameters:
        b_plasma_poloidal_average (float): Poloidal field (T).
        b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
        nd_plasma_electrons_vol_avg (float): Electron density (m^-3).
        nd_plasma_fuel_ions_vol_avg (float): Fuel ion density (m^-3).
        nd_plasma_ions_total_vol_avg (float): Total ion density (m^-3).
        temp_plasma_electron_density_weighted_kev (float): Density-weighted electron temperature (keV).
        temp_plasma_ion_density_weighted_kev (float): Density-weighted ion temperature (keV).
        pden_alpha_total_mw (float): Alpha power per unit volume, from beams and plasma (MW/m^3).
        pden_plasma_alpha_mw (float): Alpha power per unit volume just from plasma (MW/m^3).
        i_beta_fast_alpha (int): Switch for fast alpha pressure method.

    Returns:
        float: Fast alpha beta component.

    Notes:
        - For IPDG89 scaling applicability is Z_eff = 1.5, T_i/T_e = 1, <T> = 5-20 keV


    References:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
          https://inis.iaea.org/collection/NCLCollectionStore/_Public/21/068/21068960.pdf

        - Uckan, N. A., Tolliver, J. S., Houlberg, W. A., and Attenberger, S. E.
          Influence of fast alpha diffusion and thermal alpha buildup on tokamak reactor performance.
          United States: N. p., 1987. Web.https://www.osti.gov/servlets/purl/5611706

    """

    # Determine average fast alpha density
    if physics_variables.f_plasma_fuel_deuterium < 1.0:
        beta_thermal = (
            2.0
            * constants.RMU0
            * constants.KILOELECTRON_VOLT
            * (
                nd_plasma_electrons_vol_avg * temp_plasma_electron_density_weighted_kev
                + nd_plasma_ions_total_vol_avg * temp_plasma_ion_density_weighted_kev
            )
            / (b_plasma_toroidal_on_axis**2 + b_plasma_poloidal_average**2)
        )

        # jlion: This "fact" model is heavily flawed for smaller temperatures! It is unphysical for a stellarator (high n low T)
        # IPDG89 fast alpha scaling
        if i_beta_fast_alpha == 0:
            fact = min(
                0.3,
                0.29
                * (nd_plasma_fuel_ions_vol_avg / nd_plasma_electrons_vol_avg) ** 2
                * (
                    (
                        temp_plasma_electron_density_weighted_kev
                        + temp_plasma_ion_density_weighted_kev
                    )
                    / 20.0
                    - 0.37
                ),
            )

        # Modified scaling, D J Ward
        else:
            fact = min(
                0.30,
                0.26
                * (nd_plasma_fuel_ions_vol_avg / nd_plasma_electrons_vol_avg) ** 2
                * np.sqrt(
                    max(
                        0.0,
                        (
                            (
                                temp_plasma_electron_density_weighted_kev
                                + temp_plasma_ion_density_weighted_kev
                            )
                            / 20.0
                            - 0.65
                        ),
                    )
                ),
            )

        fact = max(fact, 0.0)
        fact2 = pden_alpha_total_mw / pden_plasma_alpha_mw
        beta_fast_alpha = beta_thermal * fact * fact2

    else:  # negligible alpha production, alpha_power_density = p_beam_alpha_mw = 0
        beta_fast_alpha = 0.0

    return beta_fast_alpha
