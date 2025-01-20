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
