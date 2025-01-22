import numpy as np


def neo_alcator_confinement_time(
    dene20: float, rminor: float, rmajor: float, qstar: float
) -> float:
    """
    Calculate the Nec-Alcator(NA) OH scaling confinement time

    Parameters:
    dene20 (float): Volume averaged electron density in units of 10**20 m**-3
    rminor (float): Plasma minor radius [m]
    rmajor (float): Plasma major radius [m]
    qstar (float): Equivalent cylindrical edge safety factor

    Returns:
    float: Neo-Alcator confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
         ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.

    """
    return 0.07e0 * dene20 * rminor * rmajor * rmajor * qstar


def mirnov_confinement_time(rminor: float, kappa95: float, pcur: float) -> float:
    """
    Calculate the Mirnov scaling (H-mode) confinement time

    Parameters:
    hfact (float): H-factor
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    pcur (float): Plasma current [MA]

    Returns:
    float: Mirnov scaling confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
         ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return 0.2e0 * rminor * np.sqrt(kappa95) * pcur


def merezhkin_muhkovatov_confinement_time(
    rmajor: float,
    rminor: float,
    kappa95: float,
    qstar: float,
    dnla20: float,
    afuel: float,
    ten: float,
) -> float:
    """
    Calculate the Merezhkin-Mukhovatov (MM) OH/L-mode scaling confinement time

    Parameters:
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    qstar (float): Equivalent cylindrical edge safety factor
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    afuel (float): Fuel atomic mass number
    ten (float): Electron temperature [keV]

    Returns:
    float: Merezhkin-Mukhovatov confinement time [s]

     Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
         ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return (
        3.5e-3
        * rmajor**2.75e0
        * rminor**0.25e0
        * kappa95**0.125e0
        * qstar
        * dnla20
        * np.sqrt(afuel)
        / np.sqrt(ten / 10.0e0)
    )


def shimomura_confinement_time(
    rmajor: float, rminor: float, bt: float, kappa95: float, afuel: float
) -> float:
    """
    Calculate the  Shimomura (S) optimized H-mode scaling confinement time

    Parameters:
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    bt (float): Toroidal magnetic field [T]
    kappa95 (float): Plasma elongation at 95% flux surface
    afuel (float): Fuel atomic mass number

    Returns:
    float: Shimomura confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
         ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return 0.045e0 * rmajor * rminor * bt * np.sqrt(kappa95) * np.sqrt(afuel)


def kaye_goldston_confinement_time(
    kappa95: float,
    pcur: float,
    n20: float,
    rmajor: float,
    afuel: float,
    bt: float,
    rminor: float,
    powerht: float,
) -> float:
    """
    Calculate the Kaye-Goldston (KG) L-mode scaling confinement time

    Parameters:
    hfact (float): H-factor
    kappa95 (float): Plasma elongation at 95% flux surface
    pcur (float): Plasma current [MA]
    n20 (float): Line averaged electron density in units of 10**20 m**-3
    rmajor (float): Plasma major radius [m]
    afuel (float): Fuel atomic mass number
    bt (float): Toroidal magnetic field [T]
    rminor (float): Plasma minor radius [m]
    powerht (float): Net Heating power [MW]

    Returns:
    float: Kaye-Goldston confinement time [s]

    Notes:
        - An isotope correction factor (M_i/1.5)^0.5 is added to the original scaling to reflect the fact
          that the empirical fits to the data were from experiments with H and D mixture, M_i = 1.5

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
         ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return (
        0.055e0
        * kappa95**0.28e0
        * pcur**1.24e0
        * n20**0.26e0
        * rmajor**1.65e0
        * np.sqrt(afuel / 1.5e0)
        / (bt**0.09e0 * rminor**0.49e0 * powerht**0.58e0)
    )


def iter_89P_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa: float,
    dnla20: float,
    bt: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the ITER Power scaling - ITER 89-P (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa (float): Plasma elongation
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: ITER 89-P confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
          ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return (
        0.048e0
        * pcur**0.85e0
        * rmajor**1.2e0
        * rminor**0.3e0
        * np.sqrt(kappa)
        * dnla20**0.1e0
        * bt**0.2e0
        * np.sqrt(afuel)
        / np.sqrt(powerht)
    )


def iter_89_0_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa: float,
    dnla20: float,
    bt: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the ITER Offset linear scaling - ITER 89-O (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa (float): Plasma elongation
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: ITER 89-O confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """
    term1 = (
        0.04e0
        * pcur**0.5e0
        * rmajor**0.3e0
        * rminor**0.8e0
        * kappa**0.6e0
        * afuel**0.5e0
    )
    term2 = (
        0.064e0
        * pcur**0.8e0
        * rmajor**1.6e0
        * rminor**0.6e0
        * kappa**0.5e0
        * dnla20**0.6e0
        * bt**0.35e0
        * afuel**0.2e0
        / powerht
    )
    return term1 + term2


def rebut_lallia_confinement_time(
    rminor: float,
    rmajor: float,
    kappa: float,
    afuel: float,
    pcur: float,
    zeff: float,
    dnla20: float,
    bt: float,
    powerht: float,
) -> float:
    """
    Calculate the Rebut-Lallia offset linear scaling (L-mode) confinement time

    Parameters:
    rminor (float): Plasma minor radius [m]
    rmajor (float): Plasma major radius [m]
    kappa (float): Plasma elongation at 95% flux surface
    afuel (float): Fuel atomic mass number
    pcur (float): Plasma current [MA]
    zeff (float): Effective charge
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]

    Returns:
    float: Rebut-Lallia confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    rll = (rminor**2 * rmajor * kappa) ** (1.0e0 / 3.0e0)
    term1 = 1.2e-2 * pcur * rll**1.5e0 / np.sqrt(zeff)
    term2 = (
        0.146e0
        * dnla20**0.75e0
        * np.sqrt(pcur)
        * np.sqrt(bt)
        * rll**2.75e0
        * zeff**0.25e0
        / powerht
    )
    return 1.65e0 * np.sqrt(afuel / 2.0e0) * (term1 + term2)


def goldston_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the Goldston scaling (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: Goldston confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
        ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.

    """
    return (
        0.037e0
        * pcur
        * rmajor**1.75e0
        * rminor ** (-0.37e0)
        * np.sqrt(kappa95)
        * np.sqrt(afuel / 1.5e0)
        / np.sqrt(powerht)
    )


def t10_confinement_time(
    dnla20: float,
    rmajor: float,
    qstar: float,
    bt: float,
    rminor: float,
    kappa95: float,
    powerht: float,
    zeff: float,
    pcur: float,
) -> float:
    """
    Calculate the T-10 scaling confinement time

    Parameters:
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    rmajor (float): Plasma major radius [m]
    qstar (float): Equivalent cylindrical edge safety factor
    bt (float): Toroidal magnetic field [T]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    powerht (float): Net Heating power [MW]
    zeff (float): Effective charge
    pcur (float): Plasma current [MA]

    Returns:
    float: T-10 confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    denfac = dnla20 * rmajor * qstar / (1.3e0 * bt)
    denfac = min(1.0e0, denfac)
    tauee = (
        0.095e0
        * rmajor
        * rminor
        * bt
        * np.sqrt(kappa95)
        * denfac
        / powerht**0.4e0
        * (zeff**2 * pcur**4 / (rmajor * rminor * qstar**3 * kappa95**1.5e0)) ** 0.08e0
    )
    return tauee


def jaeri_confinement_time(
    kappa95: float,
    rminor: float,
    afuel: float,
    n20: float,
    pcur: float,
    bt: float,
    rmajor: float,
    qstar: float,
    zeff: float,
    powerht: float,
) -> float:
    """
    Calculate the JAERI / Odajima-Shimomura L-mode scaling confinement time

    Parameters:
    kappa95 (float): Plasma elongation at 95% flux surface
    rminor (float): Plasma minor radius [m]
    afuel (float): Fuel atomic mass number
    n20 (float): Line averaged electron density in units of 10**20 m**-3
    pcur (float): Plasma current [MA]
    bt (float): Toroidal magnetic field [T]
    rmajor (float): Plasma major radius [m]
    qstar (float): Equivalent cylindrical edge safety factor
    zeff (float): Effective charge
    powerht (float): Net Heating power [MW]

    Returns:
    float: JAERI confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
          ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    gjaeri = (
        zeff**0.4e0
        * ((15.0e0 - zeff) / 20.0e0) ** 0.6e0
        * (3.0e0 * qstar * (qstar + 5.0e0) / ((qstar + 2.0e0) * (qstar + 7.0e0)))
        ** 0.6e0
    )
    return (
        0.085e0 * kappa95 * rminor**2 * np.sqrt(afuel)
        + 0.069e0
        * n20**0.6e0
        * pcur
        * bt**0.2e0
        * rminor**0.4e0
        * rmajor**1.6e0
        * np.sqrt(afuel)
        * gjaeri
        * kappa95**0.2e0
        / powerht
    )


def kaye_big_confinement_time(
    rmajor: float,
    rminor: float,
    bt: float,
    kappa95: float,
    pcur: float,
    n20: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the Kaye-Big scaling confinement time

    Parameters:
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    bt (float): Toroidal magnetic field [T]
    kappa95 (float): Plasma elongation at 95% flux surface
    pcur (float): Plasma current [MA]
    n20 (float): Line averaged electron density in units of 10**20 m**-3
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: Kaye-Big confinement time [s]

    Notes:

    References:
        - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            ‘ITER physics design guidelines: 1989’, no. No. 10. Feb. 1990.
    """
    return (
        0.105e0
        * np.sqrt(rmajor)
        * rminor**0.8e0
        * bt**0.3e0
        * kappa95**0.25e0
        * pcur**0.85e0
        * n20**0.1e0
        * np.sqrt(afuel)
        / np.sqrt(powerht)
    )


def iter_h90_p_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa: float,
    dnla20: float,
    bt: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the ITER H-mode scaling - ITER H90-P confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa (float): Plasma elongation
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: ITER H90-P confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    return (
        0.064e0
        * pcur**0.87e0
        * rmajor**1.82e0
        * rminor ** (-0.12e0)
        * kappa**0.35e0
        * dnla20**0.09e0
        * bt**0.15e0
        * np.sqrt(afuel)
        / np.sqrt(powerht)
    )


def riedel_l_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    dnla20: float,
    bt: float,
    powerht: float,
) -> float:
    """
    Calculate the Riedel scaling (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]

    Returns:
    float: Riedel confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    return (
        0.044e0
        * pcur**0.93e0
        * rmajor**1.37e0
        * rminor ** (-0.049e0)
        * kappa95**0.588e0
        * dnla20**0.078e0
        * bt**0.152e0
        / powerht**0.537e0
    )


def christiansen_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    dnla20: float,
    bt: float,
    powerht: float,
    afuel: float,
) -> float:
    """
    Calculate the Christiansen et al scaling (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]
    afuel (float): Fuel atomic mass number

    Returns:
    float: Christiansen confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    return (
        0.24e0
        * pcur**0.79e0
        * rmajor**0.56e0
        * rminor**1.46e0
        * kappa95**0.73e0
        * dnla20**0.41e0
        * bt**0.29e0
        / (powerht**0.79e0 * afuel**0.02e0)
    )


def lackner_gottardi_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    dnla20: float,
    bt: float,
    powerht: float,
) -> float:
    """
    Calculate the Lackner-Gottardi scaling (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]

    Returns:
    float: Lackner-Gottardi confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """
    qhat = (1.0e0 + kappa95**2) * rminor**2 * bt / (0.4e0 * pcur * rmajor)
    return (
        0.12e0
        * pcur**0.8e0
        * rmajor**1.8e0
        * rminor**0.4e0
        * kappa95
        * (1.0e0 + kappa95) ** (-0.8e0)
        * dnla20**0.6e0
        * qhat**0.4e0
        / powerht**0.6e0
    )


def neo_kaye_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    dnla20: float,
    bt: float,
    powerht: float,
) -> float:
    """
    Calculate the Neo-Kaye scaling (L-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]

    Returns:
    float: Neo-Kaye confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    return (
        0.063e0
        * pcur**1.12e0
        * rmajor**1.3e0
        * rminor ** (-0.04e0)
        * kappa95**0.28e0
        * dnla20**0.14e0
        * bt**0.04e0
        / powerht**0.59e0
    )


def riedel_h_confinement_time(
    pcur: float,
    rmajor: float,
    rminor: float,
    kappa95: float,
    dnla20: float,
    bt: float,
    afuel: float,
    powerht: float,
) -> float:
    """
    Calculate the Riedel scaling (H-mode) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    rmajor (float): Plasma major radius [m]
    rminor (float): Plasma minor radius [m]
    kappa95 (float): Plasma elongation at 95% flux surface
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    bt (float): Toroidal magnetic field [T]
    afuel (float): Fuel atomic mass number
    powerht (float): Net Heating power [MW]

    Returns:
    float: Riedel H-mode confinement time [s]

    Notes:

    References:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
    """
    return (
        0.1e0
        * np.sqrt(afuel)
        * pcur**0.884e0
        * rmajor**1.24e0
        * rminor ** (-0.23e0)
        * kappa95**0.317e0
        * bt**0.207e0
        * dnla20**0.105e0
        / powerht**0.486e0
    )


def iter_h90_p_amended_confinement_time(
    pcur: float,
    bt: float,
    afuel: float,
    rmajor: float,
    powerht: float,
    kappa: float,
) -> float:
    """
        Calculate the amended ITER H90-P confinement time

        Parameters:
        pcur (float): Plasma current [MA]
        bt (float): Toroidal magnetic field [T]
        afuel (float): Fuel atomic mass number
        rmajor (float): Plasma major radius [m]
        powerht (float): Net Heating power [MW]
        kappa (float): Plasma elongation

        Returns:
        float: Amended ITER H90-P confinement time [s]

        Notes:

        References:
            - J. P. Christiansen et al., “Global energy confinement H-mode database for ITER,”
            Nuclear Fusion, vol. 32, no. 2, pp. 291–338, Feb. 1992,
            doi: https://doi.org/10.1088/0029-5515/32/2/i11.
    ‌
    """
    return (
        0.082e0
        * pcur**1.02e0
        * bt**0.15e0
        * np.sqrt(afuel)
        * rmajor**1.60e0
        / (powerht**0.47e0 * kappa**0.19e0)
    )


def sudo_et_al_confinement_time(
    rmajor: float,
    rminor: float,
    dnla20: float,
    bt: float,
    powerht: float,
) -> float:
    """
        Calculate the Sudo et al. scaling confinement time

        Parameters:
        rmajor (float): Plasma major radius [m]
        rminor (float): Plasma minor radius [m]
        dnla20 (float): Line averaged electron density in units of 10**20 m**-3
        bt (float): Toroidal magnetic field [T]
        powerht (float): Net Heating power [MW]

        Returns:
        float: Sudo et al. confinement time [s]

        Notes:

        References:
            - S. Sudo et al., “Scalings of energy confinement and density limit in stellarator/heliotron devices,”
            Nuclear Fusion, vol. 30, no. 1, pp. 11–21, Jan. 1990,
            doi: https://doi.org/10.1088/0029-5515/30/1/002.
    ‌"""

    return (
        0.17e0
        * rmajor**0.75e0
        * rminor**2
        * dnla20**0.69e0
        * bt**0.84e0
        * powerht ** (-0.58e0)
    )


def gyro_reduced_bohm_confinement_time(
    bt: float,
    dnla20: float,
    powerht: float,
    rminor: float,
    rmajor: float,
) -> float:
    """
    Calculate the Gyro-reduced Bohm scaling confinement time

    Parameters:
    bt (float): Toroidal magnetic field [T]
    dnla20 (float): Line averaged electron density in units of 10**20 m**-3
    powerht (float): Net Heating power [MW]
    rminor (float): Plasma minor radius [m]
    rmajor (float): Plasma major radius [m]

    Returns:
    float: Gyro-reduced Bohm confinement time [s]

    Notes:

    References:
        - Goldston, R. J., H. Biglari, and G. W. Hammett. "E× B/B 2 vs. µ B/B as the Cause of Transport in Tokamaks."
          Bull. Am. Phys. Soc 34 (1989): 1964.
    """
    return (
        0.25e0
        * bt**0.8e0
        * dnla20**0.6e0
        * powerht ** (-0.6e0)
        * rminor**2.4e0
        * rmajor**0.6e0
    )


def lackner_gottardi_stellarator_confinement_time(
    rmajor: float,
    rminor: float,
    dnla20: float,
    bt: float,
    powerht: float,
    q: float,
) -> float:
    """
        Calculate the Lackner-Gottardi stellarator scaling confinement time

        Parameters:
        rmajor (float): Plasma major radius [m]
        rminor (float): Plasma minor radius [m]
        dnla20 (float): Line averaged electron density in units of 10**20 m**-3
        bt (float): Toroidal magnetic field [T]
        powerht (float): Net Heating power [MW]
        q (float): Edge safety factor

        Returns:
        float: Lackner-Gottardi stellarator confinement time [s]

        Notes:

        References:
            - K. Lackner and N. A. O. Gottardi, “Tokamak confinement in relation to plateau scaling,”
            Nuclear Fusion, vol. 30, no. 4, pp. 767–770, Apr. 1990,
            doi: https://doi.org/10.1088/0029-5515/30/4/018.
    ‌
    """
    return (
        0.17e0
        * rmajor
        * rminor**2
        * dnla20**0.6e0
        * bt**0.8e0
        * powerht ** (-0.6e0)
        * q**0.4e0
    )


def iter_h97p_confinement_time(
    pcur: float,
    bt: float,
    powerht: float,
    dnla19: float,
    rmajor: float,
    aspect: float,
    kappa: float,
    afuel: float,
) -> float:
    """
        Calculate the ELM-free ITER H-mode scaling - ITER H97-P confinement time

        Parameters:
        pcur (float): Plasma current [MA]
        bt (float): Toroidal magnetic field [T]
        powerht (float): Net Heating power [MW]
        dnla19 (float): Line averaged electron density in units of 10**19 m**-3
        rmajor (float): Plasma major radius [m]
        aspect (float): Aspect ratio
        kappa (float): Plasma elongation
        afuel (float): Fuel atomic mass number

        Returns:
        float: ITER H97-P confinement time [s]

        Notes:

        References:
            - I. C. Database and M. W. G. (presented Cordey), “Energy confinement scaling and the extrapolation to ITER,”
            Plasma Physics and Controlled Fusion, vol. 39, no. 12B, pp. B115–B127, Dec. 1997,
            doi: https://doi.org/10.1088/0741-3335/39/12b/009.
    ‌
    """
    return (
        0.031e0
        * pcur**0.95e0
        * bt**0.25e0
        * powerht ** (-0.67e0)
        * dnla19**0.35e0
        * rmajor**1.92e0
        * aspect ** (-0.08e0)
        * kappa**0.63e0
        * afuel**0.42e0
    )


def iter_h97p_elmy_confinement_time(
    pcur: float,
    bt: float,
    powerht: float,
    dnla19: float,
    rmajor: float,
    aspect: float,
    kappa: float,
    afuel: float,
) -> float:
    """
    Calculate the ELMy ITER H-mode scaling - ITER H97-P(y) confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    bt (float): Toroidal magnetic field [T]
    powerht (float): Net Heating power [MW]
    dnla19 (float): Line averaged electron density in units of 10**19 m**-3
    rmajor (float): Plasma major radius [m]
    aspect (float): Aspect ratio
    kappa (float): Plasma elongation
    afuel (float): Fuel atomic mass number

    Returns:
    float: ITER H97-P(y) confinement time [s]

    Notes:

    References:
        - I. C. Database and M. W. G. (presented Cordey), “Energy confinement scaling and the extrapolation to ITER,”
          Plasma Physics and Controlled Fusion, vol. 39, no. 12B, pp. B115–B127, Dec. 1997,
          doi: https://doi.org/10.1088/0741-3335/39/12b/009.

        - International Atomic Energy Agency, Vienna (Austria), ‘Technical basis for the ITER final design report, cost review and safety analysis (FDR)’,
          no. no.16. Dec. 1998.
    """
    return (
        0.029e0
        * pcur**0.90e0
        * bt**0.20e0
        * powerht ** (-0.66e0)
        * dnla19**0.40e0
        * rmajor**2.03e0
        * aspect ** (-0.19e0)
        * kappa**0.92e0
        * afuel**0.2e0
    )


def iter_96p_confinement_time(
    pcur: float,
    bt: float,
    kappa95: float,
    rmajor: float,
    aspect: float,
    dnla19: float,
    afuel: float,
    powerht: float,
) -> float:
    """
        Calculate the ITER-96P (= ITER-97L) L-mode scaling confinement time

        Parameters:
        pcur (float): Plasma current [MA]
        bt (float): Toroidal magnetic field [T]
        kappa95 (float): Plasma elongation at 95% flux surface
        rmajor (float): Plasma major radius [m]
        aspect (float): Aspect ratio
        dnla19 (float): Line averaged electron density in units of 10**19 m**-3
        afuel (float): Fuel atomic mass number
        powerht (float): Net Heating power [MW]

        Returns:
        float: ITER-96P confinement time [s]

        Notes:
            - The thermal energy confinement time is given below

        References:
            - S. B. Kaye et al., “ITER L mode confinement database,”
            Nuclear Fusion, vol. 37, no. 9, pp. 1303–1328, Sep. 1997,
            doi: https://doi.org/10.1088/0029-5515/37/9/i10.
    ‌
    """
    return (
        0.023e0
        * pcur**0.96e0
        * bt**0.03e0
        * kappa95**0.64e0
        * rmajor**1.83e0
        * aspect**0.06e0
        * dnla19**0.40e0
        * afuel**0.20e0
        * powerht ** (-0.73e0)
    )


def iter_ipb98y_confinement_time(
    pcur: float,
    bt: float,
    dnla19: float,
    powerht: float,
    rmajor: float,
    kappa: float,
    aspect: float,
    afuel: float,
) -> float:
    """
    Calculate the IPB98(y) ELMy H-mode scaling confinement time

    Parameters:
    pcur (float): Plasma current [MA]
    bt (float): Toroidal magnetic field [T]
    dnla19 (float): Line averaged electron density in units of 10**19 m**-3
    powerht (float): Net Heating power [MW]
    rmajor (float): Plasma major radius [m]
    kappa (float): Plasma separatrix elongation
    aspect (float): Aspect ratio
    afuel (float): Fuel atomic mass number

    Returns:
    float: IPB98(y) ELMy H-mode confinement time [s]

    Notes:
        - Unlike the other IPB98 scaling laws, the IPB98(y) scaling law uses the true separatrix elongation.
        - See correction paper below for more information

    References:
        - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
        Nuclear Fusion, vol. 39, no. 12, pp. 2175–2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

        - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
          Nuclear Fusion, vol. 48, no. 9, pp. 099801–099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.

    """
    return (
        0.0365e0
        * pcur**0.97e0
        * bt**0.08e0
        * dnla19**0.41e0
        * powerht ** (-0.63e0)
        * rmajor**1.93e0
        * kappa**0.67e0
        * aspect ** (-0.23e0)
        * afuel**0.2e0
    )


if __name__ == "__main__":
    pass
