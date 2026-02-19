import logging
from dataclasses import dataclass

import numpy as np

import process.impurity_radiation as impurity
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
            - pden_plasma_sync_mw (float): Synchrotron radiation power per unit volume (MW/m^3).
            - pden_plasma_core_rad_mw (float): Total core radiation power per unit volume (MW/m^3).
            - pden_plasma_outer_rad_mw (float): Edge radiation power per unit volume (MW/m^3).
            - pden_plasma_rad_mw (float): Total radiation power per unit volume (MW/m^3).

    References
    ----------
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
    """Calculate the synchrotron radiation power in MW/m^3.

    This function computes the synchrotron radiation power density for the plasma based on
    the plasma shape, major and minor radii, electron density, and temperature profiles.

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
