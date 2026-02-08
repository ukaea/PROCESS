import logging
import math
from enum import IntEnum

import numba as nb
import numpy as np
import scipy
import scipy.integrate as integrate
from scipy.optimize import root_scalar

import process.models.physics.confinement_time as confinement
import process.models.physics.fusion_reactions as reactions
import process.models.physics.impurity_radiation as impurity_radiation
import process.models.physics.l_h_transition as transition
import process.models.physics.radiation_power as physics_funcs
from process import constants
from process import process_output as po
from process.data_structure import (
    build_variables,
    constraint_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    impurity_radiation_module,
    numerics,
    physics_variables,
    pulse_variables,
    reinke_variables,
    stellarator_variables,
    times_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


@nb.jit(nopython=True, cache=True)
def rether(
    alphan,
    alphat,
    nd_plasma_electrons_vol_avg,
    dlamie,
    te,
    temp_plasma_ion_vol_avg_kev,
    n_charge_plasma_effective_mass_weighted_vol_avg,
):
    """Routine to find the equilibration power between the
    ions and electrons
    This routine calculates the equilibration power between the
    ions and electrons.
    Unknown origin

    Parameters
    ----------
    alphan :
        density profile index
    alphat :
        temperature profile index
    nd_plasma_electrons_vol_avg :
        electron density (/m3)
    dlamie :
        ion-electron coulomb logarithm
    te :
        electron temperature (keV)
    temp_plasma_ion_vol_avg_kev :
        ion temperature (keV)
    n_charge_plasma_effective_mass_weighted_vol_avg :
        mass weighted plasma effective charge

    Returns
    -------
    pden_ion_electron_equilibration_mw  :
        ion/electron equilibration power (MW/m3)

    """
    profie = (1.0 + alphan) ** 2 / (
        (2.0 * alphan - 0.5 * alphat + 1.0) * np.sqrt(1.0 + alphat)
    )
    conie = (
        2.42165e-41
        * dlamie
        * nd_plasma_electrons_vol_avg**2
        * n_charge_plasma_effective_mass_weighted_vol_avg
        * profie
    )

    return conie * (temp_plasma_ion_vol_avg_kev - te) / (te**1.5)


# -----------------------------------------------------
# Plasma Current & Poloidal Field Calculations
# -----------------------------------------------------


@nb.jit(nopython=True, cache=True)
def _plascar_bpol(
    aspect: float, eps: float, kappa: float, delta: float
) -> tuple[float, float, float, float]:
    """Calculate the poloidal field coefficients for determining the plasma current
    and poloidal field.


    This internal function calculates the poloidal field coefficients,
    which is used to calculate the poloidal field and the plasma current.

    Parameters
    ----------
    aspect :
        plasma aspect ratio
    eps :
        inverse aspect ratio
    kappa :
        plasma elongation
    delta :
        plasma triangularity

    Returns
    -------
    :
        coefficients ff1, ff2, d1, d2

    References
    ----------
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document

    """
    # Original coding, only suitable for TARTs [STAR Code]

    c1 = (kappa**2 / (1.0 + delta)) + delta
    c2 = (kappa**2 / (1.0 - delta)) - delta

    d1 = (kappa / (1.0 + delta)) ** 2 + 1.0
    d2 = (kappa / (1.0 - delta)) ** 2 + 1.0

    c1_aspect = ((c1 * eps) - 1.0) if aspect < c1 else (1.0 - (c1 * eps))

    y1 = np.sqrt(c1_aspect / (1.0 + eps)) * ((1.0 + delta) / kappa)
    y2 = np.sqrt((c2 * eps + 1.0) / (1.0 - eps)) * ((1.0 - delta) / kappa)

    h2 = (1.0 + (c2 - 1.0) * (eps / 2.0)) / np.sqrt((1.0 - eps) * (c2 * eps + 1.0))
    f2 = (d2 * (1.0 - delta) * eps) / ((1.0 - eps) * ((c2 * eps) + 1.0))
    g = (eps * kappa) / (1.0 - (eps * delta))
    ff2 = f2 * (g + 2.0 * h2 * np.arctan(y2))

    h1 = (1.0 + (1.0 - c1) * (eps / 2.0)) / np.sqrt((1.0 + eps) * c1_aspect)
    f1 = (d1 * (1.0 + delta) * eps) / ((1.0 + eps) * (c1 * eps - 1.0))

    if aspect < c1:
        ff1 = f1 * (g - h1 * np.log((1.0 + y1) / (1.0 - y1)))
    else:
        ff1 = -f1 * (-g + 2.0 * h1 * np.arctan(y1))

    return ff1, ff2, d1, d2


def calculate_poloidal_field(
    i_plasma_current: int,
    ip: float,
    q95: float,
    aspect: float,
    eps: float,
    b_plasma_toroidal_on_axis: float,
    kappa: float,
    delta: float,
    perim: float,
    rmu0: float,
) -> float:
    """Function to calculate poloidal field from the plasma current

    This function calculates the poloidal field from the plasma current in Tesla,
    using a simple calculation using Ampere's law for conventional
    tokamaks, or for TARTs, a scaling from Peng, Galambos and
    Shipe (1992).

    Parameters
    ----------
    i_plasma_current :
        current scaling model to use
    ip :
        plasma current (A)
    q95 :
        95% flux surface safety factor
    aspect :
        plasma aspect ratio
    eps :
        inverse aspect ratio
    b_plasma_toroidal_on_axis :
        toroidal field on axis (T)
    kappa :
        plasma elongation
    delta :
        plasma triangularity
    perim :
        plasma perimeter (m)
    rmu0 :
        vacuum permeability (H/m)

    Returns
    -------
    :
        poloidal field in Tesla


    References
    ----------
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971

    """
    # Use Ampere's law using the plasma poloidal cross-section
    if i_plasma_current != 2:
        return rmu0 * ip / perim
    # Use the relation from Peng, Galambos and Shipe (1992) [STAR code] otherwise
    ff1, ff2, _, _ = _plascar_bpol(aspect, eps, kappa, delta)

    # Transform q95 to qbar
    qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

    return b_plasma_toroidal_on_axis * (ff1 + ff2) / (2.0 * np.pi * qbar)


def calculate_current_coefficient_peng(
    eps: float, len_plasma_poloidal: float, rminor: float
) -> float:
    """Calculate the plasma current scaling coefficient for the Peng scaling from the STAR code.

    Parameters
    ----------
    eps : float
        Plasma inverse aspect ratio.
    len_plasma_poloidal : float
        Plasma poloidal perimeter length [m].
    rminor : float
        Plasma minor radius [m].

    Returns
    -------
    float
        The plasma current scaling coefficient.

    """

    return (
        (1.22 - 0.68 * eps)
        / ((1.0 - eps * eps) ** 2)
        * (len_plasma_poloidal / (2.0 * np.pi * rminor)) ** 2
    )


def calculate_plasma_current_peng(
    q95: float,
    aspect: float,
    eps: float,
    rminor: float,
    b_plasma_toroidal_on_axis: float,
    kappa: float,
    delta: float,
) -> float:
    """Function to calculate plasma current (Peng scaling from the STAR code)

    This function calculates the plasma current in MA,
    using a scaling from Peng, Galambos and Shipe (1992).
    It is primarily used for Tight Aspect Ratio Tokamaks and is
    selected via i_plasma_current=2.

    Parameters
    ----------
    q95 :
        95% flux surface safety factor
    aspect :
        plasma aspect ratio
    eps :
        inverse aspect ratio
    rminor :
        plasma minor radius (m)
    b_plasma_toroidal_on_axis :
        toroidal field on axis (T)
    kappa :
        plasma elongation
    delta :
        plasma triangularity


    Returns
    -------
    :
        plasma current in MA


    References
    ----------
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971

    """

    # Transform q95 to qbar
    qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

    ff1, ff2, d1, d2 = _plascar_bpol(aspect, eps, kappa, delta)

    e1 = (2.0 * kappa) / (d1 * (1.0 + delta))
    e2 = (2.0 * kappa) / (d2 * (1.0 - delta))

    return (
        rminor
        * b_plasma_toroidal_on_axis
        / qbar
        * 5.0
        * kappa
        / (2.0 * np.pi**2)
        * (np.arcsin(e1) / e1 + np.arcsin(e2) / e2)
        * (ff1 + ff2)
    )


@nb.jit(nopython=True, cache=True)
def calculate_current_coefficient_ipdg89(
    eps: float, kappa95: float, triang95: float
) -> float:
    """Calculate the fq coefficient from the IPDG89 guidlines used in the plasma current scaling.

    This function calculates the fq coefficient used in the IPDG89 plasma current scaling,
    based on the given plasma parameters.

    Parameters
    ----------
    eps :
        plasma inverse aspect ratio
    kappa95 :
        plasma elongation 95%
    triang95 :
        plasma triangularity 95%

    Returns
    -------
    :
        the fq plasma current coefficient


    References
    ----------
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """
    return (
        0.5
        * (1.17 - 0.65 * eps)
        / ((1.0 - eps * eps) ** 2)
        * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
    )


@nb.jit(nopython=True, cache=True)
def calculate_current_coefficient_todd(
    eps: float, kappa95: float, triang95: float, model: int
) -> float:
    """Calculate the fq coefficient used in the two Todd plasma current scalings.

    This function calculates the fq coefficient based on the given plasma parameters for the two Todd scalings.

    Parameters
    ----------
    eps :
        plasma inverse aspect ratio
    kappa95 :
        plasma elongation 95%
    triang95 :
        plasma triangularity 95%
    model:

    Returns
    -------
    :
        the fq plasma current coefficient

    References
    ----------
        - D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """
    # Calculate the Todd scaling based on the model
    base_scaling = (
        (1.0 + 2.0 * eps**2)
        * ((1.0 + kappa95**2) / 0.5)
        * (1.24 - 0.54 * kappa95 + 0.3 * (kappa95**2 + triang95**2) + 0.125 * triang95)
    )
    if model == 1:
        return base_scaling
    if model == 2:
        return base_scaling * (1.0 + (abs(kappa95 - 1.2)) ** 3)
    raise ProcessValueError(f"model = {model} is an invalid option")


@nb.jit(nopython=True, cache=True)
def calculate_current_coefficient_hastie(
    alphaj: float,
    alphap: float,
    b_plasma_toroidal_on_axis: float,
    delta95: float,
    eps: float,
    kappa95: float,
    pres_plasma_on_axis: float,
    rmu0: float,
) -> float:
    """Routine to calculate the f_q coefficient for the Connor-Hastie model used for scaling the plasma current.

    This routine calculates the f_q coefficient used for scaling the plasma current,
    using the Connor-Hastie scaling


    Parameters
    ----------
    alphaj :
        the current profile index
    alphap :
        the pressure profile index
    b_plasma_toroidal_on_axis :
        the toroidal field on axis (T)
    delta95 :
        the plasma triangularity 95%
    eps :
        the inverse aspect ratio
    kappa95 :
        the plasma elongation 95%
    pres_plasma_on_axis :
        the central plasma pressure (Pa)
    rmu0 :
        the vacuum permeability (H/m)

    Returns
    -------
    :
        the F coefficient


    Reference:
        - J.W.Connor and R.J.Hastie, Culham Lab Report CLM-M106 (1985).
        https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-M106-1.pdf
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """
    # Exponent in Connor-Hastie current profile
    lamda = alphaj

    # Exponent in Connor-Hastie pressure profile
    nu = alphap

    # Central plasma beta
    beta0 = 2.0 * rmu0 * pres_plasma_on_axis / (b_plasma_toroidal_on_axis**2)

    # Plasma internal inductance
    lamp1 = 1.0 + lamda
    li = lamp1 / lamda * (lamp1 / lamda * np.log(lamp1) - 1.0)

    # T/r in AEA FUS 172
    kap1 = kappa95 + 1.0
    tr = kappa95 * delta95 / kap1**2

    # E/r in AEA FUS 172
    er = (kappa95 - 1.0) / kap1

    # T primed in AEA FUS 172
    tprime = 2.0 * tr * lamp1 / (1.0 + 0.5 * lamda)

    # E primed in AEA FUS 172
    eprime = er * lamp1 / (1.0 + lamda / 3.0)

    # Delta primed in AEA FUS 172
    deltap = (0.5 * kap1 * eps * 0.5 * li) + (beta0 / (0.5 * kap1 * eps)) * lamp1**2 / (
        1.0 + nu
    )

    # Delta/R0 in AEA FUS 172
    deltar = beta0 / 6.0 * (1.0 + 5.0 * lamda / 6.0 + 0.25 * lamda**2) + (
        0.5 * kap1 * eps
    ) ** 2 * 0.125 * (1.0 - (lamda**2) / 3.0)

    # F coefficient
    return (0.5 * kap1) ** 2 * (
        1.0
        + eps**2 * (0.5 * kap1) ** 2
        + 0.5 * deltap**2
        + 2.0 * deltar
        + 0.5 * (eprime**2 + er**2)
        + 0.5 * (tprime**2 + 4.0 * tr**2)
    )


@nb.jit(nopython=True, cache=True)
def calculate_current_coefficient_sauter(
    eps: float,
    kappa: float,
    triang: float,
) -> float:
    """Routine to calculate the f_q coefficient for the Sauter model used for scaling the plasma current.

    Parameters
    ----------
    eps :
        inverse aspect ratio
    kappa :
        plasma elongation at the separatrix
    triang :
        plasma triangularity at the separatrix

    Returns
    -------
    :
        the fq coefficient

    Reference:
        - O. Sauter, Geometric formulas for system codes including the effect of negative triangularity,
        Fusion Engineering and Design, Volume 112, 2016, Pages 633-645,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2016.04.033.

    """
    w07 = 1.0  # zero squareness - can be modified later if required

    return (
        (4.1e6 / 5.0e6)
        * (1.0 + 1.2 * (kappa - 1.0) + 0.56 * (kappa - 1.0) ** 2)
        * (1.0 + 0.09 * triang + 0.16 * triang**2)
        * (1.0 + 0.45 * triang * eps)
        / (1.0 - 0.74 * eps)
        * (1.0 + 0.55 * (w07 - 1.0))
    )


@nb.jit(nopython=True, cache=True)
def calculate_current_coefficient_fiesta(
    eps: float, kappa: float, triang: float
) -> float:
    """Calculate the fq coefficient used in the FIESTA plasma current scaling.

    This function calculates the fq coefficient based on the given plasma parameters for the FIESTA scaling.

    Parameters
    ----------
    eps :
        plasma inverse aspect ratio
    kappa :
        plasma elongation at the separatrix
    triang :
        plasma triangularity at the separatrix

    Returns
    -------
    :
        the fq plasma current coefficient


    References
    ----------
        - S.Muldrew et.al,“PROCESS”: Systems studies of spherical tokamaks, Fusion Engineering and Design,
        Volume 154, 2020, 111530, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2020.111530.

    """
    return 0.538 * (1.0 + 2.440 * eps**2.736) * kappa**2.154 * triang**0.060


# --------------------------------
# Bootstrap Current Calculations
# --------------------------------


@nb.jit(nopython=True, cache=True)
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


# -----------------------------------------------------
# Diamagnetic and Pfirsch-Schlüter Current Calculations
# -----------------------------------------------------


@nb.jit(nopython=True, cache=True)
def diamagnetic_fraction_hender(beta: float) -> float:
    """Calculate the diamagnetic fraction based on the Hender fit.

    Parameters
    ----------
    beta :
        the plasma beta value

    Returns
    -------
    :
        the diamagnetic fraction

    """
    return beta / 2.8


@nb.jit(nopython=True, cache=True)
def diamagnetic_fraction_scene(beta: float, q95: float, q0: float) -> float:
    """Calculate the diamagnetic fraction based on the SCENE fit by Tim Hender.

    Parameters
    ----------
    beta :
        the plasma beta value
    q95 :
        the normalized safety factor at 95% of the plasma radius
    q0 :
        the normalized safety factor at the magnetic axis

    Returns
    -------
    :
        the diamagnetic fraction
    """
    return beta * (0.1 * q95 / q0 + 0.44) * 0.414


@nb.jit(nopython=True, cache=True)
def ps_fraction_scene(beta: float) -> float:
    """Calculate the Pfirsch-Schlüter fraction based on the SCENE fit by Tim Hender 2019.

    Parameters
    ----------
    beta :
        the plasma beta value


    Returns
    -------
    :
        the Pfirsch-Schlüter current fraction

    """
    return -9e-2 * beta


# --------------------------------------------------
# Sauter Bootstrap Current Scaling Functions
# --------------------------------------------------


@nb.jit(nopython=True, cache=True)
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
        15.9 - 0.5 * np.log(ne[radial_elements - 1]) + np.log(tempe[radial_elements - 1])
    )


@nb.jit(nopython=True, cache=True)
def _electron_collisions_sauter(
    radial_elements: np.ndarray, tempe: np.ndarray, ne: np.ndarray
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
        * _coulomb_logarithm_sauter(radial_elements, tempe, ne)
        * ne[radial_elements - 1]
        / (tempe[radial_elements - 1] * np.sqrt(tempe[radial_elements - 1]))
    )


@nb.jit(nopython=True, cache=True)
def _electron_collisionality_sauter(
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
        _electron_collisions_sauter(radial_elements, tempe, ne)
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


@nb.jit(nopython=True, cache=True)
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


@nb.jit(nopython=True, cache=True)
def _ion_collisionality_sauter(
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
        * _ion_collisions_sauter(radial_elements, zeff, ni, tempi, amain)
        * rmajor
        / (
            np.abs(inverse_q[radial_elements - 1] + 1.0e-4)
            * sqeps[radial_elements - 1] ** 3
            * np.sqrt(tempi[radial_elements - 1] / amain[radial_elements - 1])
        )
    )


@nb.jit(nopython=True, cache=True)
def _calculate_l31_coefficient(
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
    f_trapped = _trapped_particle_fraction_sauter(radial_elements, triang, sqeps)

    # Calculated electron collisionality; nu_e*
    electron_collisionality = _electron_collisionality_sauter(
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
    return l31_coefficient * _beta_poloidal_total_sauter(
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


@nb.jit(nopython=True, cache=True)
def _calculate_l31_32_coefficient(
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
    f_trapped = _trapped_particle_fraction_sauter(radial_elements, triang, sqeps)

    # Calculated electron collisionality; nu_e*
    electron_collisionality = _electron_collisionality_sauter(
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
            * (f32ei_teff**2 - f32ei_teff**4 - 0.55 * (f32ei_teff**3 - f32ei_teff**4))
        )
        - (1.2 / (1.0 + 0.5 * charge_profile) * f32ei_teff**4)
    )

    # big_f32ee_teff + big_f32ei_teff = L32 coefficient

    # Corrections suggested by Fable, 15/05/2015
    return _beta_poloidal_sauter(
        radial_elements,
        number_of_elements,
        rmajor,
        b_plasma_toroidal_on_axis,
        ne,
        tempe,
        inverse_q,
        rho,
    ) * (big_f32ee_teff + big_f32ei_teff) + _calculate_l31_coefficient(
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
    ) * _beta_poloidal_sauter(
        radial_elements,
        number_of_elements,
        rmajor,
        b_plasma_toroidal_on_axis,
        ne,
        tempe,
        inverse_q,
        rho,
    ) / _beta_poloidal_total_sauter(
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


@nb.jit(nopython=True, cache=True)
def _calculate_l34_alpha_31_coefficient(
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
    f_trapped = _trapped_particle_fraction_sauter(radial_elements, triang, sqeps)

    # Calculated electron collisionality; nu_e*
    electron_collisionality = _electron_collisionality_sauter(
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
    ion_collisionality = _ion_collisionality_sauter(
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
        _beta_poloidal_total_sauter(
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
        - _beta_poloidal_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            tempe,
            inverse_q,
            rho,
        )
    ) * (l34_coefficient * alpha) + _calculate_l31_coefficient(
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
        - _beta_poloidal_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            b_plasma_toroidal_on_axis,
            ne,
            tempe,
            inverse_q,
            rho,
        )
        / _beta_poloidal_total_sauter(
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


@nb.jit(nopython=True, cache=True)
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


@nb.jit(nopython=True, cache=True)
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


@nb.jit(nopython=True, cache=True)
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
            ((1.0 - eps) ** 2) / ((1.0 + 1.46 * sqeps_reduced) * np.sqrt(1.0 - eps**2))
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


class Physics:
    def __init__(
        self,
        plasma_profile,
        current_drive,
        plasma_beta,
        plasma_inductance,
        plasma_current,
    ):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile
        self.current_drive = current_drive
        self.beta = plasma_beta
        self.inductance = plasma_inductance
        self.current = plasma_current

    def physics(self):
        """Routine to calculate tokamak plasma physics information
        This routine calculates all the primary plasma physics parameters for a tokamak fusion reactor.

        References
        ----------
        - M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants - Part 1: Physics
          https://www.sciencedirect.com/science/article/pii/S0920379614005961
        - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
          https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
        - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
          https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
        """

        # Calculate plasma composition
        # Issue #261 Remove old radiation model (imprad_model=0)
        self.plasma_composition()

        (
            physics_variables.m_plasma_fuel_ions,
            physics_variables.m_plasma_ions_total,
            physics_variables.m_plasma_alpha,
            physics_variables.m_plasma_electron,
            physics_variables.m_plasma,
        ) = self.calculate_plasma_masses(
            physics_variables.m_fuel_amu,
            physics_variables.m_ions_total_amu,
            physics_variables.nd_plasma_ions_total_vol_avg,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.nd_plasma_alphas_vol_avg,
            physics_variables.vol_plasma,
            physics_variables.nd_plasma_electrons_vol_avg,
        )

        # Define coulomb logarithm
        # (collisions: ion-electron, electron-electron)
        physics_variables.dlamee = (
            31.0
            - (np.log(physics_variables.nd_plasma_electrons_vol_avg) / 2.0)
            + np.log(physics_variables.temp_plasma_electron_vol_avg_kev * 1000.0)
        )
        physics_variables.dlamie = (
            31.3
            - (np.log(physics_variables.nd_plasma_electrons_vol_avg) / 2.0)
            + np.log(physics_variables.temp_plasma_electron_vol_avg_kev * 1000.0)
        )

        # Calculate plasma current
        (
            physics_variables.b_plasma_poloidal_average,
            physics_variables.qstar,
            physics_variables.plasma_current,
        ) = self.current.calculate_plasma_current(
            physics_variables.alphaj,
            physics_variables.alphap,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.eps,
            physics_variables.i_plasma_current,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pres_plasma_thermal_on_axis,
            physics_variables.len_plasma_poloidal,
            physics_variables.q95,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.triang,
            physics_variables.triang95,
        )

        # -----------------------------------------------------
        # Plasma Current Profile
        # -----------------------------------------------------

        physics_variables.alphaj_wesson = self.calculate_current_profile_index_wesson(
            qstar=physics_variables.qstar, q0=physics_variables.q0
        )

        # Map calculation methods to a dictionary
        alphaj_calculations = {
            0: physics_variables.alphaj,
            1: physics_variables.alphaj_wesson,
        }

        # Calculate alphaj based on i_alphaj
        if int(physics_variables.i_alphaj) in alphaj_calculations:
            physics_variables.alphaj = alphaj_calculations[
                int(physics_variables.i_alphaj)
            ]
        else:
            raise ProcessValueError(
                "Illegal value of i_alphaj",
                i_alphaj=physics_variables.i_alphaj,
            )

        # ==================================================

        # -----------------------------------------------------
        # Plasma Normalised Internal Inductance
        # -----------------------------------------------------

        self.inductance.run()

        # ===================================================

        # Calculate density and temperature profile quantities
        # If physics_variables.i_plasma_pedestal = 1 then set pedestal density to
        #   physics_variables.f_nd_plasma_pedestal_greenwald * Greenwald density limit
        # Note: this used to be done before plasma current
        if (physics_variables.i_plasma_pedestal == 1) and (
            physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0
        ):
            physics_variables.nd_plasma_pedestal_electron = (
                physics_variables.f_nd_plasma_pedestal_greenwald
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        if (physics_variables.i_plasma_pedestal == 1) and (
            physics_variables.f_nd_plasma_separatrix_greenwald >= 0e0
        ):
            physics_variables.nd_plasma_separatrix_electron = (
                physics_variables.f_nd_plasma_separatrix_greenwald
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        self.plasma_profile.run()

        # Calculate total magnetic field [T]
        physics_variables.b_plasma_total = np.sqrt(
            physics_variables.b_plasma_toroidal_on_axis**2
            + physics_variables.b_plasma_poloidal_average**2
        )

        # Calculate the toroidal field across the plasma
        # Calculate the toroidal field profile across the plasma (1/R dependence)
        # Double element size to include both sides of the plasma
        rho = np.linspace(
            physics_variables.rmajor - physics_variables.rminor,
            physics_variables.rmajor + physics_variables.rminor,
            2 * physics_variables.n_plasma_profile_elements,
        )

        # Calculate the inboard and outboard toroidal field
        physics_variables.b_plasma_inboard_toroidal = (
            physics_variables.rmajor
            * physics_variables.b_plasma_toroidal_on_axis
            / (physics_variables.rmajor - physics_variables.rminor)
        )

        physics_variables.b_plasma_outboard_toroidal = (
            physics_variables.rmajor
            * physics_variables.b_plasma_toroidal_on_axis
            / (physics_variables.rmajor + physics_variables.rminor)
        )

        # Avoid division by zero at the magnetic axis
        rho = np.where(rho == 0, 1e-10, rho)
        physics_variables.b_plasma_toroidal_profile = (
            physics_variables.rmajor * physics_variables.b_plasma_toroidal_on_axis / rho
        )

        # ============================================

        # -----------------------------------------------------
        # Beta Components
        # -----------------------------------------------------

        self.beta.run()

        # =======================================================

        # Set PF coil ramp times
        if pulse_variables.i_pulsed_plant != 1:
            if times_variables.i_t_current_ramp_up == 0:
                times_variables.t_plant_pulse_plasma_current_ramp_up = (
                    physics_variables.plasma_current / 5.0e5
                )
                times_variables.t_plant_pulse_coil_precharge = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )
                times_variables.t_plant_pulse_plasma_current_ramp_down = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )

        else:
            if times_variables.pulsetimings == 0.0e0:
                # times_variables.t_plant_pulse_coil_precharge is input
                times_variables.t_plant_pulse_plasma_current_ramp_up = (
                    physics_variables.plasma_current / 1.0e5
                )
                times_variables.t_plant_pulse_plasma_current_ramp_down = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )

            else:
                # times_variables.t_plant_pulse_plasma_current_ramp_up is set either in INITIAL or INPUT, or by being
                # iterated using limit equation 41.
                times_variables.t_plant_pulse_coil_precharge = max(
                    times_variables.t_plant_pulse_coil_precharge,
                    times_variables.t_plant_pulse_plasma_current_ramp_up,
                )
                # t_plant_pulse_plasma_current_ramp_down = max(t_plant_pulse_plasma_current_ramp_down,t_plant_pulse_plasma_current_ramp_up)
                times_variables.t_plant_pulse_plasma_current_ramp_down = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )

        # Reset second times_variables.t_plant_pulse_burn value (times_variables.t_burn_0).
        # This is used to ensure that the burn time is used consistently;
        # see convergence loop in fcnvmc1, evaluators.f90
        times_variables.t_burn_0 = times_variables.t_plant_pulse_burn

        # Time during the pulse in which a plasma is present
        times_variables.t_plant_pulse_plasma_present = (
            times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_fusion_ramp
            + times_variables.t_plant_pulse_burn
            + times_variables.t_plant_pulse_plasma_current_ramp_down
        )
        times_variables.t_plant_pulse_no_burn = (
            times_variables.t_plant_pulse_coil_precharge
            + times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_plasma_current_ramp_down
            + times_variables.t_plant_pulse_dwell
            + times_variables.t_plant_pulse_fusion_ramp
        )

        # Total cycle time
        times_variables.t_plant_pulse_total = (
            times_variables.t_plant_pulse_coil_precharge
            + times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_fusion_ramp
            + times_variables.t_plant_pulse_burn
            + times_variables.t_plant_pulse_plasma_current_ramp_down
            + times_variables.t_plant_pulse_dwell
        )

        # ***************************** #
        #      DIAMAGNETIC CURRENT      #
        # ***************************** #

        # Hender scaling for diamagnetic current at tight physics_variables.aspect ratio
        current_drive_variables.f_c_plasma_diamagnetic_hender = (
            diamagnetic_fraction_hender(physics_variables.beta_total_vol_avg)
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.f_c_plasma_diamagnetic_scene = (
            diamagnetic_fraction_scene(
                physics_variables.beta_total_vol_avg,
                physics_variables.q95,
                physics_variables.q0,
            )
        )

        if physics_variables.i_diamagnetic_current == 1:
            current_drive_variables.f_c_plasma_diamagnetic = (
                current_drive_variables.f_c_plasma_diamagnetic_hender
            )
        elif physics_variables.i_diamagnetic_current == 2:
            current_drive_variables.f_c_plasma_diamagnetic = (
                current_drive_variables.f_c_plasma_diamagnetic_scene
            )

        # ***************************** #
        #    PFIRSCH-SCHLÜTER CURRENT   #
        # ***************************** #

        # Pfirsch-Schlüter scaling for diamagnetic current
        current_drive_variables.f_c_plasma_pfirsch_schluter_scene = ps_fraction_scene(
            physics_variables.beta_total_vol_avg
        )

        if physics_variables.i_pfirsch_schluter_current == 1:
            current_drive_variables.f_c_plasma_pfirsch_schluter = (
                current_drive_variables.f_c_plasma_pfirsch_schluter_scene
            )

        # ***************************** #
        #       BOOTSTRAP CURRENT       #
        # ***************************** #

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

        bootstrap_map = {
            0: current_drive_variables.f_c_plasma_bootstrap,
            1: current_drive_variables.f_c_plasma_bootstrap_iter89,
            2: current_drive_variables.f_c_plasma_bootstrap_nevins,
            3: current_drive_variables.f_c_plasma_bootstrap_wilson,
            4: current_drive_variables.f_c_plasma_bootstrap_sauter,
            5: current_drive_variables.f_c_plasma_bootstrap_sakai,
            6: current_drive_variables.f_c_plasma_bootstrap_aries,
            7: current_drive_variables.f_c_plasma_bootstrap_andrade,
            8: current_drive_variables.f_c_plasma_bootstrap_hoang,
            9: current_drive_variables.f_c_plasma_bootstrap_wong,
            10: current_drive_variables.bscf_gi_i,
            11: current_drive_variables.bscf_gi_ii,
            12: current_drive_variables.f_c_plasma_bootstrap_sugiyama_l,
            13: current_drive_variables.f_c_plasma_bootstrap_sugiyama_h,
        }
        if int(physics_variables.i_bootstrap_current) in bootstrap_map:
            current_drive_variables.f_c_plasma_bootstrap = bootstrap_map[
                int(physics_variables.i_bootstrap_current)
            ]
        else:
            raise ProcessValueError(
                "Illegal value of i_bootstrap_current",
                i_bootstrap_current=physics_variables.i_bootstrap_current,
            )

        physics_variables.err242 = 0
        if (
            current_drive_variables.f_c_plasma_bootstrap
            > current_drive_variables.f_c_plasma_bootstrap_max
        ) and physics_variables.i_bootstrap_current != 0:
            current_drive_variables.f_c_plasma_bootstrap = min(
                current_drive_variables.f_c_plasma_bootstrap,
                current_drive_variables.f_c_plasma_bootstrap_max,
            )
            physics_variables.err242 = 1

        current_drive_variables.f_c_plasma_internal = (
            current_drive_variables.f_c_plasma_bootstrap
            + current_drive_variables.f_c_plasma_diamagnetic
            + current_drive_variables.f_c_plasma_pfirsch_schluter
        )

        # Plasma driven current fraction (Bootstrap + Diamagnetic
        # + Pfirsch-Schlüter) constrained to be less than
        # or equal to the total fraction of the plasma current
        # produced by non-inductive means (which also includes
        # the current drive proportion)
        physics_variables.err243 = 0
        if (
            current_drive_variables.f_c_plasma_internal
            > physics_variables.f_c_plasma_non_inductive
        ):
            current_drive_variables.f_c_plasma_internal = min(
                current_drive_variables.f_c_plasma_internal,
                physics_variables.f_c_plasma_non_inductive,
            )
            physics_variables.err243 = 1

        # Fraction of plasma current produced by inductive means
        physics_variables.f_c_plasma_inductive = max(
            1.0e-10, (1.0e0 - physics_variables.f_c_plasma_non_inductive)
        )
        #  Fraction of plasma current produced by auxiliary current drive
        physics_variables.f_c_plasma_auxiliary = (
            physics_variables.f_c_plasma_non_inductive
            - current_drive_variables.f_c_plasma_internal
        )

        # Auxiliary current drive power calculations

        if current_drive_variables.i_hcd_calculations != 0:
            self.current_drive.cudriv()

        # ***************************** #
        #        FUSION REACTIONS       #
        # ***************************** #

        # Calculate fusion power + components

        fusion_reactions = reactions.FusionReactionRate(self.plasma_profile)
        fusion_reactions.deuterium_branching(
            physics_variables.temp_plasma_ion_vol_avg_kev
        )
        fusion_reactions.calculate_fusion_rates()
        fusion_reactions.set_physics_variables()

        # This neglects the power from the beam
        physics_variables.p_plasma_dt_mw = (
            physics_variables.dt_power_density_plasma * physics_variables.vol_plasma
        )
        physics_variables.p_dhe3_total_mw = (
            physics_variables.dhe3_power_density * physics_variables.vol_plasma
        )
        physics_variables.p_dd_total_mw = (
            physics_variables.dd_power_density * physics_variables.vol_plasma
        )

        # Calculate neutral beam slowing down effects
        # If ignited, then ignore beam fusion effects

        if (current_drive_variables.c_beam_total != 0.0e0) and (
            physics_variables.i_plasma_ignited == 0
        ):
            (
                physics_variables.beta_beam,
                physics_variables.nd_beam_ions_out,
                physics_variables.p_beam_alpha_mw,
            ) = reactions.beam_fusion(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.b_plasma_poloidal_average,
                physics_variables.b_plasma_toroidal_on_axis,
                current_drive_variables.c_beam_total,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_fuel_ions_vol_avg,
                physics_variables.dlamie,
                current_drive_variables.e_beam_kev,
                physics_variables.f_plasma_fuel_deuterium,
                physics_variables.f_plasma_fuel_tritium,
                current_drive_variables.f_beam_tritium,
                physics_variables.sigmav_dt_average,
                physics_variables.temp_plasma_electron_density_weighted_kev,
                physics_variables.temp_plasma_ion_density_weighted_kev,
                physics_variables.vol_plasma,
                physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            )
            physics_variables.fusden_total = (
                physics_variables.fusden_plasma
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.fusden_alpha_total = (
                physics_variables.fusden_plasma_alpha
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.p_dt_total_mw = (
                physics_variables.p_plasma_dt_mw
                + (1.0 / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION))
                * physics_variables.p_beam_alpha_mw
            )
            physics_variables.p_beam_neutron_mw = physics_variables.p_beam_alpha_mw * (
                constants.DT_NEUTRON_ENERGY_FRACTION
                / (1 - constants.DT_NEUTRON_ENERGY_FRACTION)
            )
            physics_variables.p_beam_dt_mw = physics_variables.p_beam_alpha_mw * (
                1 / (1 - constants.DT_NEUTRON_ENERGY_FRACTION)
            )
        else:
            # If no beams present then the total alpha rates and power are the same as the plasma values
            physics_variables.fusden_total = physics_variables.fusden_plasma
            physics_variables.fusden_alpha_total = physics_variables.fusden_plasma_alpha
            physics_variables.p_dt_total_mw = physics_variables.p_plasma_dt_mw

        physics_variables.fusrat_total = (
            physics_variables.fusden_total * physics_variables.vol_plasma
        )

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.pden_neutron_total_mw,
            physics_variables.p_plasma_alpha_mw,
            physics_variables.p_alpha_total_mw,
            physics_variables.p_plasma_neutron_mw,
            physics_variables.p_neutron_total_mw,
            physics_variables.p_non_alpha_charged_mw,
            physics_variables.pden_alpha_total_mw,
            physics_variables.f_pden_alpha_electron_mw,
            physics_variables.f_pden_alpha_ions_mw,
            physics_variables.p_charged_particle_mw,
            physics_variables.p_fusion_total_mw,
        ) = reactions.set_fusion_powers(
            physics_variables.f_alpha_electron,
            physics_variables.f_alpha_ion,
            physics_variables.p_beam_alpha_mw,
            physics_variables.pden_non_alpha_charged_mw,
            physics_variables.pden_plasma_neutron_mw,
            physics_variables.vol_plasma,
            physics_variables.pden_plasma_alpha_mw,
        )

        physics_variables.beta_fast_alpha = self.beta.fast_alpha_beta(
            b_plasma_poloidal_average=physics_variables.b_plasma_poloidal_average,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
            nd_plasma_fuel_ions_vol_avg=physics_variables.nd_plasma_fuel_ions_vol_avg,
            nd_plasma_ions_total_vol_avg=physics_variables.nd_plasma_ions_total_vol_avg,
            temp_plasma_electron_density_weighted_kev=physics_variables.temp_plasma_electron_density_weighted_kev,
            temp_plasma_ion_density_weighted_kev=physics_variables.temp_plasma_ion_density_weighted_kev,
            pden_alpha_total_mw=physics_variables.pden_alpha_total_mw,
            pden_plasma_alpha_mw=physics_variables.pden_plasma_alpha_mw,
            i_beta_fast_alpha=physics_variables.i_beta_fast_alpha,
        )

        # Nominal mean neutron wall load on entire first wall area including divertor and beam holes
        # Note that 'a_fw_total' excludes these, so they have been added back in.
        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.ffwal
                * physics_variables.p_neutron_total_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if divertor_variables.n_divertors == 2:
                # Double null configuration
                physics_variables.pflux_fw_neutron_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - 2.0e0 * fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_neutron_total_mw
                    / build_variables.a_fw_total
                )
            else:
                # Single null Configuration
                physics_variables.pflux_fw_neutron_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_neutron_total_mw
                    / build_variables.a_fw_total
                )

        # Calculate ion/electron equilibration power

        physics_variables.pden_ion_electron_equilibration_mw = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.dlamie,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.temp_plasma_ion_vol_avg_kev,
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
        )

        # Calculate radiation power

        radpwrdata = physics_funcs.calculate_radiation_powers(
            self.plasma_profile,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.rminor,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.aspect,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.tbeta,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.f_sync_reflect,
            physics_variables.rmajor,
            physics_variables.kappa,
            physics_variables.vol_plasma,
        )
        physics_variables.pden_plasma_sync_mw = radpwrdata.pden_plasma_sync_mw
        physics_variables.pden_plasma_core_rad_mw = radpwrdata.pden_plasma_core_rad_mw
        physics_variables.pden_plasma_outer_rad_mw = radpwrdata.pden_plasma_outer_rad_mw
        physics_variables.pden_plasma_rad_mw = radpwrdata.pden_plasma_rad_mw

        physics_variables.p_plasma_sync_mw = (
            physics_variables.pden_plasma_sync_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_inner_rad_mw = (
            physics_variables.pden_plasma_core_rad_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_outer_rad_mw = (
            physics_variables.pden_plasma_outer_rad_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_rad_mw = (
            physics_variables.pden_plasma_rad_mw * physics_variables.vol_plasma
        )

        # Calculate ohmic power
        (
            physics_variables.pden_plasma_ohmic_mw,
            physics_variables.p_plasma_ohmic_mw,
            physics_variables.f_res_plasma_neo,
            physics_variables.res_plasma,
        ) = self.plasma_ohmic_heating(
            physics_variables.f_c_plasma_inductive,
            physics_variables.kappa95,
            physics_variables.plasma_current,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.vol_plasma,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        # Calculate L- to H-mode power threshold for different scalings
        physics_variables.l_h_threshold_powers = l_h_threshold_power(
            physics_variables.nd_plasma_electron_line,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.kappa,
            physics_variables.a_plasma_surface,
            physics_variables.m_ions_total_amu,
            physics_variables.aspect,
            physics_variables.plasma_current,
        )

        # Enforced L-H power threshold value (if constraint 15 is turned on)
        physics_variables.p_l_h_threshold_mw = physics_variables.l_h_threshold_powers[
            physics_variables.i_l_h_threshold - 1
        ]

        # Power transported to the divertor by charged particles,
        # i.e. excludes neutrons and radiation, and also NBI orbit loss power,
        # which is assumed to be absorbed by the first wall
        pinj = (
            current_drive_variables.p_hcd_injected_total_mw
            if physics_variables.i_plasma_ignited == 0
            else 0.0
        )

        physics_variables.p_plasma_separatrix_mw = (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + pinj
            + physics_variables.p_plasma_ohmic_mw
            - physics_variables.p_plasma_rad_mw
        )

        physics_variables.plfux_plasma_surface_neutron_avg_mw = (
            physics_variables.p_neutron_total_mw / physics_variables.a_plasma_surface
        )

        # KLUDGE: Ensure p_plasma_separatrix_mw is continuously positive (physical, rather than
        # negative potential power), as required by other models (e.g.
        # Physics.calculate_density_limit())
        physics_variables.p_plasma_separatrix_mw = (
            physics_variables.p_plasma_separatrix_mw
            / (1 - np.exp(-physics_variables.p_plasma_separatrix_mw))
        )

        # if double null configuration share the power
        # over the upper and lower divertor, where physics_variables.f_p_div_lower gives
        # the factor of power conducted to the lower divertor
        if divertor_variables.n_divertors == 2:
            physics_variables.p_div_lower_separatrix_mw = (
                physics_variables.f_p_div_lower
                * physics_variables.p_plasma_separatrix_mw
            )
            physics_variables.p_div_upper_separatrix_mw = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_variables.p_plasma_separatrix_mw
            physics_variables.p_div_separatrix_max_mw = max(
                physics_variables.p_div_lower_separatrix_mw,
                physics_variables.p_div_upper_separatrix_mw,
            )

        # Resistive diffusion time = current penetration time ~ mu0.a^2/resistivity
        physics_variables.t_plasma_res_diffusion = res_diff_time(
            physics_variables.rmajor,
            physics_variables.res_plasma,
            physics_variables.kappa95,
        )

        # Power transported to the first wall by escaped alpha particles
        physics_variables.p_fw_alpha_mw = physics_variables.p_alpha_total_mw * (
            1.0e0 - physics_variables.f_p_alpha_plasma_deposited
        )

        # Density limit
        (
            physics_variables.nd_plasma_electron_max_array,
            physics_variables.nd_plasma_electrons_max,
        ) = self.calculate_density_limit(
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.i_density_limit,
            physics_variables.p_plasma_separatrix_mw,
            current_drive_variables.p_hcd_injected_total_mw,
            physics_variables.plasma_current,
            divertor_variables.prn1,
            physics_variables.qstar,
            physics_variables.q95,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.a_plasma_surface,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        # Calculate transport losses and energy confinement time using the
        # chosen scaling law
        (
            physics_variables.pden_electron_transport_loss_mw,
            physics_variables.pden_ion_transport_loss_mw,
            physics_variables.t_electron_energy_confinement,
            physics_variables.t_energy_confinement,
            physics_variables.t_ion_energy_confinement,
            physics_variables.p_plasma_loss_mw,
        ) = self.calculate_confinement_time(
            physics_variables.m_fuel_amu,
            physics_variables.p_alpha_total_mw,
            physics_variables.aspect,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.nd_plasma_ions_total_vol_avg,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_electron_line,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.i_confinement_time,
            physics_variables.i_plasma_ignited,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.p_non_alpha_charged_mw,
            current_drive_variables.p_hcd_injected_total_mw,
            physics_variables.plasma_current,
            physics_variables.pden_plasma_core_rad_mw,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.temp_plasma_ion_density_weighted_kev,
            physics_variables.q95,
            physics_variables.qstar,
            physics_variables.vol_plasma,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        # Total transport power from scaling law (MW)
        physics_variables.p_electron_transport_loss_mw = (
            physics_variables.pden_electron_transport_loss_mw
            * physics_variables.vol_plasma
        )
        physics_variables.p_ion_transport_loss_mw = (
            physics_variables.pden_ion_transport_loss_mw * physics_variables.vol_plasma
        )

        # Calculate Volt-second requirements
        (
            physics_variables.vs_plasma_internal,
            physics_variables.ind_plasma,
            physics_variables.vs_plasma_burn_required,
            physics_variables.vs_plasma_ramp_required,
            physics_variables.vs_plasma_ind_ramp,
            physics_variables.vs_plasma_res_ramp,
            physics_variables.vs_plasma_total_required,
            physics_variables.v_plasma_loop_burn,
        ) = self.inductance.calculate_volt_second_requirements(
            physics_variables.csawth,
            physics_variables.eps,
            physics_variables.f_c_plasma_inductive,
            physics_variables.ejima_coeff,
            physics_variables.kappa,
            physics_variables.rmajor,
            physics_variables.res_plasma,
            physics_variables.plasma_current,
            times_variables.t_plant_pulse_fusion_ramp,
            times_variables.t_plant_pulse_burn,
            physics_variables.ind_plasma_internal_norm,
        )

        physics_variables.e_plasma_magnetic_stored = (
            0.5e0 * physics_variables.ind_plasma * physics_variables.plasma_current**2
        )

        # Calculate auxiliary physics related information
        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.ntau,
            physics_variables.nTtau,
            physics_variables.figmer,
            physics_variables.fusrat,
            physics_variables.molflow_plasma_fuelling_required,
            physics_variables.rndfuel,
            physics_variables.t_alpha_confinement,
            physics_variables.f_alpha_energy_confinement,
        ) = self.phyaux(
            physics_variables.aspect,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.fusden_total,
            physics_variables.fusden_alpha_total,
            physics_variables.plasma_current,
            sbar,
            physics_variables.nd_plasma_alphas_vol_avg,
            physics_variables.t_energy_confinement,
            physics_variables.vol_plasma,
        )

        # Total transport power from scaling law (MW)
        physics_variables.pscalingmw = (
            physics_variables.p_electron_transport_loss_mw
            + physics_variables.p_ion_transport_loss_mw
        )

        # ============================================================

        # MDK
        # Nominal mean photon wall load on entire first wall area including divertor and beam holes
        # Note that 'a_fw_total' excludes these, so they have been added back in.
        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_rad_mw = (
                physics_variables.ffwal
                * physics_variables.p_plasma_rad_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if divertor_variables.n_divertors == 2:
                # Double Null configuration in - including SoL radiation
                physics_variables.pflux_fw_rad_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - 2.0e0 * fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / build_variables.a_fw_total
                    + (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - 2.0e0 * fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.rad_fraction_sol
                    * physics_variables.p_plasma_separatrix_mw
                    / (build_variables.a_fw_total)
                )
            else:
                # Single null configuration - including SoL radaition
                physics_variables.pflux_fw_rad_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / build_variables.a_fw_total
                    + (
                        1.0e0
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.rad_fraction_sol
                    * physics_variables.p_plasma_separatrix_mw
                    / build_variables.a_fw_total
                )

        constraint_variables.pflux_fw_rad_max_mw = (
            physics_variables.pflux_fw_rad_mw * constraint_variables.f_fw_rad_max
        )

        # Calculate the target imbalances
        # find the total power into the targets
        physics_variables.ptarmw = physics_variables.p_plasma_separatrix_mw * (
            1.0e0 - physics_variables.rad_fraction_sol
        )
        # use physics_variables.f_p_div_lower to find deltarsep
        # Parameters taken from double null machine
        # D. Brunner et al
        physics_variables.lambdaio = 1.57e-3

        # Issue #1559 Infinities in physics_module.drsep when running single null in a double null machine
        # C W Ashe
        if physics_variables.f_p_div_lower < 4.5e-5:
            physics_variables.drsep = 1.5e-2
        elif physics_variables.f_p_div_lower > (1.0e0 - 4.5e-5):
            physics_variables.drsep = -1.5e-2
        else:
            physics_variables.drsep = (
                -2.0e0
                * 1.5e-3
                * math.atanh(2.0e0 * (physics_variables.f_p_div_lower - 0.5e0))
            )
        # Model Taken from D3-D paper for conventional divertor
        # Journal of Nuclear Materials
        # Volumes 290-293, March 2001, Pages 935-939
        # Find the innner and outer lower target imbalance

        physics_variables.fio = 0.16e0 + (0.16e0 - 0.41e0) * (
            1.0e0
            - (
                2.0e0
                / (
                    1.0e0
                    + np.exp(
                        -((physics_variables.drsep / physics_variables.lambdaio) ** 2)
                    )
                )
            )
        )
        if divertor_variables.n_divertors == 2:
            # Double Null configuration
            # Find all the power fractions accross the targets
            # Taken from D3-D conventional divertor design
            physics_variables.fli = (
                physics_variables.f_p_div_lower * physics_variables.fio
            )
            physics_variables.flo = physics_variables.f_p_div_lower * (
                1.0e0 - physics_variables.fio
            )
            physics_variables.fui = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_variables.fio
            physics_variables.fuo = (1.0e0 - physics_variables.f_p_div_lower) * (
                1.0e0 - physics_variables.fio
            )
            # power into each target
            physics_variables.plimw = physics_variables.fli * physics_variables.ptarmw
            physics_variables.plomw = physics_variables.flo * physics_variables.ptarmw
            physics_variables.puimw = physics_variables.fui * physics_variables.ptarmw
            physics_variables.puomw = physics_variables.fuo * physics_variables.ptarmw
        else:
            # Single null configuration
            physics_variables.fli = physics_variables.fio
            physics_variables.flo = 1.0e0 - physics_variables.fio
            # power into each target
            physics_variables.plimw = physics_variables.fli * physics_variables.ptarmw
            physics_variables.plomw = physics_variables.flo * physics_variables.ptarmw

        # Calculate some derived quantities that may not have been defined earlier
        physics_variables.p_plasma_heating_total_mw = 1e6 * (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            + current_drive_variables.p_hcd_injected_total_mw
        )
        physics_variables.f_p_plasma_separatrix_rad = (
            1.0e6
            * physics_variables.p_plasma_rad_mw
            / physics_variables.p_plasma_heating_total_mw
        )
        physics_variables.rad_fraction_total = (
            physics_variables.f_p_plasma_separatrix_rad
            + (1.0e0 - physics_variables.f_p_plasma_separatrix_rad)
            * physics_variables.rad_fraction_sol
        )
        physics_variables.pradsolmw = (
            physics_variables.rad_fraction_sol * physics_variables.p_plasma_separatrix_mw
        )

        if 78 in numerics.icc:
            po.write(
                self.outfile,
                (
                    f"reinke t and fz, physics = {physics_variables.temp_plasma_separatrix_kev} , {reinke_variables.fzmin}"
                ),
            )
            fgw = (
                physics_variables.nd_plasma_electron_max_array[6]
                / physics_variables.nd_plasma_electrons_vol_avg
            )
            # calculate separatrix temperature, if Reinke criterion is used
            physics_variables.temp_plasma_separatrix_kev = reinke_tsep(
                physics_variables.b_plasma_toroidal_on_axis,
                physics_variables.p_plasma_separatrix_mw
                / physics_variables.p_l_h_threshold_mw,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.eps,
                fgw,
                physics_variables.kappa,
                reinke_variables.lhat,
            )

            if reinke_variables.fzmin >= 1.0e0:
                logger.error(
                    "REINKE IMPURITY MODEL: fzmin is greater than or equal to 1.0, this is at least notable"
                )

            po.write(
                self.outfile,
                (
                    f" 'fzactual, frac, reinke_variables.impvardiv = {reinke_variables.fzactual},"
                    f" {impurity_radiation_module.f_nd_impurity_electron_array(reinke_variables.impvardiv)},"
                    f" {reinke_variables.impvardiv}"
                ),
            )

    @staticmethod
    def calculate_current_profile_index_wesson(qstar: float, q0: float) -> float:
        """Calculate the Wesson current profile index.

        Parameters
        ----------
        qstar : float
            Cylindrical safety factor.
        q0 : float
            Safety factor on axis.

        Returns
        -------
        float
            The Wesson current profile index.

        Notes
        -----
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              normalised internal inductance.
            - This relation is only true for the cyclindrical plasma approximation.

        References
        ----------
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.

        """
        return qstar / q0 - 1.0

    @staticmethod
    def calculate_density_limit(
        b_plasma_toroidal_on_axis: float,
        i_density_limit: int,
        p_plasma_separatrix_mw: float,
        p_hcd_injected_total_mw: float,
        plasma_current: float,
        prn1: float,
        qcyl: float,
        q95: float,
        rmajor: float,
        rminor: float,
        a_plasma_surface: float,
        zeff: float,
    ) -> tuple[np.ndarray, float]:
        """Calculate the density limit using various models.

        Parameters
        ----------
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        i_density_limit : int
            Switch denoting which formula to enforce.
        p_plasma_separatrix_mw : float
            Power flowing to the edge plasma via charged particles (MW).
        p_hcd_injected_total_mw : float
            Power injected into the plasma (MW).
        plasma_current : float
            Plasma current (A).
        prn1 : float
            Edge density / average plasma density.
        qcyl : float
            Equivalent cylindrical safety factor (qstar).
        q95 : float
            Safety factor at 95% surface.
        rmajor : float
            Plasma major radius (m).
        rminor : float
            Plasma minor radius (m).
        a_plasma_surface : float
            Plasma surface area (m^2).
        zeff : float
            Plasma effective charge.

        Returns
        -------
        Tuple[np.ndarray, float]
            A tuple containing:
            - nd_plasma_electron_max_array (np.ndarray): Average plasma density limit using seven different models (m^-3).
            - nd_plasma_electrons_max (float): Enforced average plasma density limit (m^-3).

        Raises
        ------
        ValueError
            If i_density_limit is not between 1 and 7.
            Notes
            -----
        ValueError
            If i_density_limit is not between 1 and 7.
            Notes
            -----
            This routine calculates several different formulae for the density limit and enforces the one chosen by the user.
        ValueError
            If i_density_limit is not between 1 and 7.
            Notes
            -----
            This routine calculates several different formulae for the density limit and enforces the one chosen by the user.
            For i_density_limit = 1-5, 8, we scale the sepatrix density limit output by the ratio of the separatrix to volume averaged density

        References
        ----------
        - AEA FUS 172
            Physics Assessment for the European Reactor Study
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines
            1989
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines
            1989
            - M. Bernert et al., “The H-mode density limit in the full tungsten ASDEX Upgrade tokamak,”
            vol. 57, no. 1, pp. 014038-014038, Nov. 2014, doi: https://doi.org/10.1088/0741-3335/57/1/014038. ‌

        """

        if i_density_limit < 1 or i_density_limit > 7:
            raise ProcessValueError(
                "Illegal value for i_density_limit", i_density_limit=i_density_limit
            )

        nd_plasma_electron_max_array = np.empty((8,))

        # Power per unit area crossing the plasma edge
        # (excludes radiation and neutrons)

        p_perp = p_plasma_separatrix_mw / a_plasma_surface

        # Old ASDEX density limit formula
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        nd_plasma_electron_max_array[0] = (
            1.54e20
            * p_perp**0.43
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.45
        ) / prn1

        # Borrass density limit model for ITER (I)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # Borrass et al, ITER-TN-PH-9-6 (1989)

        nd_plasma_electron_max_array[1] = (
            1.8e20
            * p_perp**0.53
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.22
        ) / prn1

        # Borrass density limit model for ITER (II)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # This formula is (almost) identical to that in the original routine
        # denlim (now deleted).

        nd_plasma_electron_max_array[2] = (
            0.5e20
            * p_perp**0.57
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.09
        ) / prn1

        # JET edge radiation density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # qcyl=qstar here, but literature is not clear.

        denom = (zeff - 1.0) * (1.0 - 4.0 / (3.0 * qcyl))
        if denom <= 0.0:
            if i_density_limit == 4:
                logger.error(
                    f"qcyl < 4/3; nd_plasma_electron_max_array(4) set to zero; model 5 will be enforced instead. {denom=} {qcyl=}"
                )
                i_density_limit = 5

            nd_plasma_electron_max_array[3] = 0.0
        else:
            nd_plasma_electron_max_array[3] = (
                1.0e20 * np.sqrt(p_hcd_injected_total_mw / denom)
            ) / prn1

        # JET simplified density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        nd_plasma_electron_max_array[4] = (
            0.237e20
            * b_plasma_toroidal_on_axis
            * np.sqrt(p_plasma_separatrix_mw)
            / rmajor
        ) / prn1

        # Hugill-Murakami M.q limit
        # qcyl=qstar here, which is okay according to the literature

        nd_plasma_electron_max_array[5] = (
            3.0e20 * b_plasma_toroidal_on_axis / (rmajor * qcyl)
        )

        # Greenwald limit

        nd_plasma_electron_max_array[6] = (
            1.0e14 * plasma_current / (np.pi * rminor * rminor)
        )

        nd_plasma_electron_max_array[7] = (
            1.0e20
            * 0.506
            * (p_hcd_injected_total_mw**0.396 * (plasma_current / 1.0e6) ** 0.265)
            / (q95**0.323)
        ) / prn1

        # Enforce the chosen density limit

        return nd_plasma_electron_max_array, nd_plasma_electron_max_array[
            i_density_limit - 1
        ]

    @staticmethod
    def plasma_composition():
        """Calculates various plasma component fractional makeups.

        This subroutine determines the various plasma component fractional makeups.
        It is the replacement for the original routine betcom(), and is used in conjunction
        with the new impurity radiation model.

        This function performs the following calculations:
        - Determines the alpha ash portion.
        - Calculates the proton density.
        - Calculates the beam hot ion component.
        - Sums the ion densities for all impurity ions with charge greater than helium.
        - Ensures charge neutrality by adjusting the fuel portion.
        - Calculates the total ion density.
        - Sets the relative impurity densities for radiation calculations.
        - Calculates the effective charge.
        - Defines the Coulomb logarithm for ion-electron and electron-electron collisions.
        - Calculates the fraction of alpha energy to ions and electrons.
        - Calculates the average atomic masses of injected fuel species and neutral beams.
        - Calculates the density weighted mass and mass weighted plasma effective charge.
        """

        # Alpha ash portion
        physics_variables.nd_plasma_alphas_vol_avg = (
            physics_variables.nd_plasma_electrons_vol_avg
            * physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Protons
        # This calculation will be wrong on the first call as the particle
        # production rates are evaluated later in the calling sequence
        # Issue #557 Allow f_nd_protium_electrons impurity to be specified: 'f_nd_protium_electrons'
        # This will override the calculated value which is a minimum.
        if physics_variables.fusden_alpha_total < 1.0e-6:  # not calculated yet...
            physics_variables.nd_plasma_protons_vol_avg = max(
                physics_variables.f_nd_protium_electrons
                * physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_alphas_vol_avg
                * (physics_variables.f_plasma_fuel_helium3 + 1.0e-3),
            )  # rough estimate
        else:
            physics_variables.nd_plasma_protons_vol_avg = max(
                physics_variables.f_nd_protium_electrons
                * physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_alphas_vol_avg
                * physics_variables.proton_rate_density
                / physics_variables.fusden_alpha_total,
            )

        # ======================================================================

        # Beam hot ion component
        # If ignited, prevent beam fusion effects
        if physics_variables.i_plasma_ignited == 0:
            physics_variables.nd_beam_ions = (
                physics_variables.nd_plasma_electrons_vol_avg
                * physics_variables.f_nd_beam_electron
            )
        else:
            physics_variables.nd_beam_ions = 0.0

        # ======================================================================

        # Sum of Zi.ni for all impurity ions (those with charge > helium)
        znimp = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                znimp += impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.temp_plasma_electron_vol_avg_kev])
                ).squeeze() * (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * physics_variables.nd_plasma_electrons_vol_avg
                )

        # ======================================================================

        # Fuel portion - conserve charge neutrality
        # znfuel is the sum of Zi.ni for the three fuel ions
        znfuel = (
            physics_variables.nd_plasma_electrons_vol_avg
            - 2.0 * physics_variables.nd_plasma_alphas_vol_avg
            - physics_variables.nd_plasma_protons_vol_avg
            - physics_variables.nd_beam_ions
            - znimp
        )

        # ======================================================================

        # Fuel ion density, nd_plasma_fuel_ions_vol_avg
        # For D-T-He3 mix, nd_plasma_fuel_ions_vol_avg = nD + nT + nHe3, while znfuel = nD + nT + 2*nHe3
        # So nd_plasma_fuel_ions_vol_avg = znfuel - nHe3 = znfuel - f_plasma_fuel_helium3*nd_plasma_fuel_ions_vol_avg
        physics_variables.nd_plasma_fuel_ions_vol_avg = znfuel / (
            1.0 + physics_variables.f_plasma_fuel_helium3
        )

        # ======================================================================

        # Set hydrogen and helium relative impurity densities for
        # radiation calculations
        impurity_radiation_module.f_nd_impurity_electron_array[
            impurity_radiation.element2index("H_")
        ] = (
            physics_variables.nd_plasma_protons_vol_avg
            + (
                physics_variables.f_plasma_fuel_deuterium
                + physics_variables.f_plasma_fuel_tritium
            )
            * physics_variables.nd_plasma_fuel_ions_vol_avg
            + physics_variables.nd_beam_ions
        ) / physics_variables.nd_plasma_electrons_vol_avg

        impurity_radiation_module.f_nd_impurity_electron_array[
            impurity_radiation.element2index("He")
        ] = (
            physics_variables.f_plasma_fuel_helium3
            * physics_variables.nd_plasma_fuel_ions_vol_avg
            / physics_variables.nd_plasma_electrons_vol_avg
            + physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Total impurity density
        physics_variables.nd_plasma_impurities_vol_avg = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.nd_plasma_impurities_vol_avg += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * physics_variables.nd_plasma_electrons_vol_avg
                )

        # ======================================================================

        # Total ion density
        physics_variables.nd_plasma_ions_total_vol_avg = (
            physics_variables.nd_plasma_fuel_ions_vol_avg
            + physics_variables.nd_plasma_alphas_vol_avg
            + physics_variables.nd_plasma_protons_vol_avg
            + physics_variables.nd_beam_ions
            + physics_variables.nd_plasma_impurities_vol_avg
        )

        # ======================================================================

        # Set some relative impurity density variables
        # for the benefit of other routines
        physics_variables.f_nd_plasma_carbon_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("C_")
            ]
        )
        physics_variables.f_nd_plasma_oxygen_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("O_")
            ]
        )
        physics_variables.f_nd_plasma_iron_argon_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("Fe")
            ]
            + impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("Ar")
            ]
        )

        # ======================================================================

        # Effective charge
        # Calculation should be sum(ni.Zi^2) / sum(ni.Zi),
        # but ne = sum(ni.Zi) through quasineutrality
        physics_variables.n_charge_plasma_effective_vol_avg = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            physics_variables.n_charge_plasma_effective_vol_avg += (
                impurity_radiation_module.f_nd_impurity_electron_array[imp]
                * impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.temp_plasma_electron_vol_avg_kev])
                ).squeeze()
                ** 2
            )

        # ======================================================================

        # Fraction of alpha energy to ions and electrons
        # From Max Fenstermacher
        # (used with electron and ion power balance equations only)
        # No consideration of pden_non_alpha_charged_mw here...

        # f_temp_plasma_electron_density_vol_avg now calculated in plasma_profiles, after the very first
        # call of plasma_composition; use old parabolic profile estimate
        # in this case
        if physics_variables.first_call == 1:
            pc = (
                (1.0 + physics_variables.alphan)
                * (1.0 + physics_variables.alphat)
                / (1.0 + physics_variables.alphan + physics_variables.alphat)
            )
            physics_variables.first_call = 0
        else:
            pc = physics_variables.f_temp_plasma_electron_density_vol_avg

        physics_variables.f_alpha_electron = 0.88155 * np.exp(
            -physics_variables.temp_plasma_electron_vol_avg_kev * pc / 67.4036
        )
        physics_variables.f_alpha_ion = 1.0 - physics_variables.f_alpha_electron

        # ======================================================================

        # Average atomic masses of injected fuel species
        physics_variables.m_fuel_amu = (
            (constants.M_DEUTERON_AMU * physics_variables.f_plasma_fuel_deuterium)
            + (constants.M_TRITON_AMU * physics_variables.f_plasma_fuel_tritium)
            + (constants.M_HELION_AMU * physics_variables.f_plasma_fuel_helium3)
        )

        # ======================================================================

        # Average atomic masses of injected fuel species in the neutral beams
        # Only deuterium and tritium in the beams
        physics_variables.m_beam_amu = (
            constants.M_DEUTERON_AMU * (1.0 - current_drive_variables.f_beam_tritium)
        ) + (constants.M_TRITON_AMU * current_drive_variables.f_beam_tritium)

        # ======================================================================

        # Average mass of all ions
        physics_variables.m_ions_total_amu = (
            (
                physics_variables.m_fuel_amu
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            + (constants.M_ALPHA_AMU * physics_variables.nd_plasma_alphas_vol_avg)
            + (physics_variables.nd_plasma_protons_vol_avg * constants.M_PROTON_AMU)
            + (physics_variables.m_beam_amu * physics_variables.nd_beam_ions)
        )
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.m_ions_total_amu += (
                    physics_variables.nd_plasma_electrons_vol_avg
                    * impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation_module.m_impurity_amu_array[imp]
                )

        physics_variables.m_ions_total_amu = (
            physics_variables.m_ions_total_amu
            / physics_variables.nd_plasma_ions_total_vol_avg
        )

        # ======================================================================

        # Mass weighted plasma effective charge
        # Sum of (Zi^2*n_i) / m_i
        physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg = (
            (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_DEUTERON_AMU
            )
            + (
                physics_variables.f_plasma_fuel_tritium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_TRITON_AMU
            )
            + (
                4.0
                * physics_variables.f_plasma_fuel_helium3
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_HELION_AMU
            )
            + (4.0 * physics_variables.nd_plasma_alphas_vol_avg / constants.M_ALPHA_AMU)
            + (physics_variables.nd_plasma_protons_vol_avg / constants.M_PROTON_AMU)
            + (
                (1.0 - current_drive_variables.f_beam_tritium)
                * physics_variables.nd_beam_ions
                / constants.M_DEUTERON_AMU
            )
            + (
                current_drive_variables.f_beam_tritium
                * physics_variables.nd_beam_ions
                / constants.M_TRITON_AMU
            )
        ) / physics_variables.nd_plasma_electrons_vol_avg
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation.zav_of_te(
                        imp,
                        np.array([physics_variables.temp_plasma_electron_vol_avg_kev]),
                    ).squeeze()
                    ** 2
                    / impurity_radiation_module.m_impurity_amu_array[imp]
                )

        # ======================================================================

    @staticmethod
    def phyaux(
        aspect: float,
        nd_plasma_electrons_vol_avg: float,
        te: float,
        nd_plasma_fuel_ions_vol_avg: float,
        fusden_total: float,
        fusden_alpha_total: float,
        plasma_current: float,
        sbar: float,
        nd_plasma_alphas_vol_avg: float,
        t_energy_confinement: float,
        vol_plasma: float,
    ) -> tuple[float, float, float, float, float, float, float, float]:
        """Auxiliary physics quantities

        Parameters
        ----------
        aspect : float
            Plasma aspect ratio.
        nd_plasma_electrons_vol_avg : float
            Electron density (/m3).
        te : float
            Volume avergaed electron temperature (keV).
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (/m3).
        fusden_total : float
            Fusion reaction rate from plasma and beams (/m3/s).
        fusden_alpha_total : float
            Alpha particle production rate (/m3/s).
        plasma_current : float
            Plasma current (A).
        sbar : float
            Exponent for aspect ratio (normally 1).
        nd_plasma_alphas_vol_avg : float
            Alpha ash density (/m3).
        t_energy_confinement : float
            Global energy confinement time (s).
        vol_plasma : float
            Plasma volume (m3).

        Returns
        -------
        tuple
            A tuple containing:
            - burnup (float): Fractional plasma burnup.
            - ntau (float): Plasma average n-tau (s/m3).
            - nTtau (float): Plasma triple product nT-tau (s/m3).
            - figmer (float): Physics figure of merit.
            - fusrat (float): Number of fusion reactions per second.
            - molflow_plasma_fuelling_required (float): Fuelling rate for D-T (nucleus-pairs/sec).
            - rndfuel (float): Fuel burnup rate (reactions/s).
            - t_alpha_confinement (float): Alpha particle confinement time (s).
            - f_alpha_energy_confinement (float): Fraction of alpha energy confinement.
            This subroutine calculates extra physics related items needed by other parts of the code.

        """
        figmer = 1e-6 * plasma_current * aspect**sbar

        ntau = t_energy_confinement * nd_plasma_electrons_vol_avg
        nTtau = ntau * te

        # Fusion reactions per second
        fusrat = fusden_total * vol_plasma

        # Alpha particle confinement time (s)
        # Number of alphas / alpha production rate
        if fusden_alpha_total != 0.0:
            t_alpha_confinement = nd_plasma_alphas_vol_avg / fusden_alpha_total
        else:  # only likely if DD is only active fusion reaction
            t_alpha_confinement = 0.0

        # Fractional burnup
        # (Consider detailed model in: G. L. Jackson, V. S. Chan, R. D. Stambaugh,
        # Fusion Science and Technology, vol.64, no.1, July 2013, pp.8-12)
        # The ratio of ash to fuel particle confinement times is given by
        # tauratio
        # Possible logic...
        # burnup = fuel ion-pairs burned/m3 / initial fuel ion-pairs/m3;
        # fuel ion-pairs burned/m3 = alpha particles/m3 (for both D-T and D-He3 reactions)
        # initial fuel ion-pairs/m3 = burnt fuel ion-pairs/m3 + unburnt fuel-ion pairs/m3
        # Remember that unburnt fuel-ion pairs/m3 = 0.5 * unburnt fuel-ions/m3
        if physics_variables.burnup_in <= 1.0e-9:
            burnup = (
                nd_plasma_alphas_vol_avg
                / (nd_plasma_alphas_vol_avg + 0.5 * nd_plasma_fuel_ions_vol_avg)
                / physics_variables.tauratio
            )
        else:
            burnup = physics_variables.burnup_in

        # Fuel burnup rate (reactions/second) (previously Amps)
        rndfuel = fusrat

        # Required fuelling rate (fuel ion pairs/second) (previously Amps)
        molflow_plasma_fuelling_required = rndfuel / burnup

        f_alpha_energy_confinement = t_alpha_confinement / t_energy_confinement

        return (
            burnup,
            ntau,
            nTtau,
            figmer,
            fusrat,
            molflow_plasma_fuelling_required,
            rndfuel,
            t_alpha_confinement,
            f_alpha_energy_confinement,
        )

    @staticmethod
    def plasma_ohmic_heating(
        f_c_plasma_inductive: float,
        kappa95: float,
        plasma_current: float,
        rmajor: float,
        rminor: float,
        temp_plasma_electron_density_weighted_kev: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float]:
        """Calculate the ohmic heating power and related parameters.

        Parameters
        ----------
        f_c_plasma_inductive : float
            Fraction of plasma current driven inductively.
        kappa95 : float
            Plasma elongation at 95% surface.
        plasma_current : float
            Plasma current (A).
        rmajor : float
            Major radius (m).
        rminor : float
            Minor radius (m).
        temp_plasma_electron_density_weighted_kev : float
            Density weighted average electron temperature (keV).
        vol_plasma : float
            Plasma volume (m^3).
        zeff : float
            Plasma effective charge.

        Returns
        -------
        Tuple[float, float, float, float]
            Tuple containing:
            - pden_plasma_ohmic_mw (float): Ohmic heating power per unit volume (MW/m^3).
            - p_plasma_ohmic_mw (float): Total ohmic heating power (MW).
            - f_res_plasma_neo (float): Neo-classical resistivity enhancement factor.
            - res_plasma (float): Plasma resistance (ohm).

        References
        ----------
        - ITER Physics Design Guidelines
            1989 [IPDG89], N. A. Uckan et al,

        """
        # Density weighted electron temperature in 10 keV units
        t10 = temp_plasma_electron_density_weighted_kev / 10.0

        # Plasma resistance, from loop voltage calculation in ITER Physics Design Guidelines: 1989
        res_plasma = (
            physics_variables.plasma_res_factor
            * 2.15e-9
            * zeff
            * rmajor
            / (kappa95 * rminor**2 * t10**1.5)
        )

        # Neo-classical resistivity enhancement factor
        # Taken from ITER Physics Design Guidelines: 1989
        # The expression is valid for aspect ratios in the range 2.5 to 4.0

        f_res_plasma_neo = (
            1.0 if 2.5 >= rmajor / rminor <= 4.0 else 4.3 - 0.6 * rmajor / rminor
        )

        res_plasma = res_plasma * f_res_plasma_neo

        # Check to see if plasma resistance is negative
        # (possible if aspect ratio is too high)
        if res_plasma <= 0.0:
            logger.error(
                f"Negative plasma resistance res_plasma. {res_plasma=} {physics_variables.aspect=}"
            )

        # Ohmic heating power per unit volume
        # Corrected from: pden_plasma_ohmic_mw = (f_c_plasma_inductive*plasma_current)**2 * ...

        pden_plasma_ohmic_mw = (
            f_c_plasma_inductive * plasma_current**2 * res_plasma * 1.0e-6 / vol_plasma
        )

        # Total ohmic heating power
        p_plasma_ohmic_mw = pden_plasma_ohmic_mw * vol_plasma

        return pden_plasma_ohmic_mw, p_plasma_ohmic_mw, f_res_plasma_neo, res_plasma

    @staticmethod
    def calculate_plasma_current(
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        eps: float,
        i_plasma_current: int,
        kappa: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        len_plasma_poloidal: float,
        q95: float,
        rmajor: float,
        rminor: float,
        triang: float,
        triang95: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma current.

        This routine calculates the plasma current based on the edge safety factor q95.
        It will also make the current profile parameters consistent with the q-profile if required.

        Parameters
        ----------
        alphaj : float
            Current profile index.
        alphap : float
            Pressure profile index.
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        eps : float
            Inverse aspect ratio.
        i_plasma_current : int
            Current scaling model to use.
            1 = Peng analytic fit
            2 = Peng divertor scaling (TART,STAR)
            3 = Simple ITER scaling
            4 = IPDG89 scaling
            5 = Todd empirical scaling I
            6 = Todd empirical scaling II
            7 = Connor-Hastie model
            8 = Sauter scaling (allowing negative triangularity)
            9 = FIESTA ST scaling
        kappa : float
            Plasma elongation.
        kappa95 : float
            Plasma elongation at 95% surface.
        pres_plasma_on_axis : float
            Central plasma pressure (Pa).
        len_plasma_poloidal : float
            Plasma perimeter length (m).
        q95 : float
            Plasma safety factor at 95% flux (= q-bar for i_plasma_current=2).
        ind_plasma_internal_norm : float
            Plasma normalised internal inductance.
        rmajor : float
            Major radius (m).
        rminor : float
            Minor radius (m).
        triang : float
            Plasma triangularity.
        triang95 : float
            Plasma triangularity at 95% surface.

        Returns
        -------
        Tuple[float, float, float,]
            Tuple containing b_plasma_poloidal_average, qstar, plasma_current,

        Raises
        ------
        ValueError
            If invalid value for i_plasma_current is provided.
        ValueError
            If invalid value for i_plasma_current is provided.
        ValueError
            If invalid value for i_plasma_current is provided.

        References
        ----------
        - J D Galambos, STAR Code
            Spherical Tokamak Analysis and Reactor Code, unpublished internal Oak Ridge document
        - J D Galambos, STAR Code
            Spherical Tokamak Analysis and Reactor Code, unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971
        - ITER Physics Design Guidelines
            1989 [IPDG89], N. A. Uckan et al, ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        - M. Kovari et al, 2014, "PROCESS"
            A systems code for fusion power plants - Part 1: Physics
        - M. Kovari et al, 2014, "PROCESS"
            A systems code for fusion power plants - Part 1: Physics
            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
        - M. Kovari et al, 2014, "PROCESS"
            A systems code for fusion power plants - Part 1: Physics
            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
            - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
        - M. Kovari et al, 2014, "PROCESS"
            A systems code for fusion power plants - Part 1: Physics
            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
            - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
            - Sauter, Geometric formulas for systems codes..., FED 2016

        """
        # Aspect ratio
        aspect_ratio = 1.0 / eps

        # Only the Sauter scaling (i_plasma_current=8) is suitable for negative triangularity:
        if i_plasma_current != 8 and triang < 0.0:
            raise ProcessValueError(
                f"Triangularity is negative without i_plasma_current = 8 selected: {triang=}, {i_plasma_current=}"
            )

        # Calculate the function Fq that scales the edge q from the
        # circular cross-section cylindrical case

        # Peng analytical fit
        if i_plasma_current == 1:
            fq = calculate_current_coefficient_peng(eps, len_plasma_poloidal, rminor)

        # Peng scaling for double null divertor; TARTs [STAR Code]
        elif i_plasma_current == 2:
            plasma_current = 1.0e6 * calculate_plasma_current_peng(
                q95, aspect_ratio, eps, rminor, b_plasma_toroidal_on_axis, kappa, triang
            )

        # Simple ITER scaling (simply the cylindrical case)
        elif i_plasma_current == 3:
            fq = 1.0

        # ITER formula (IPDG89)
        elif i_plasma_current == 4:
            fq = calculate_current_coefficient_ipdg89(eps, kappa95, triang95)

        # Todd empirical scalings
        # D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        elif i_plasma_current in [5, 6]:
            fq = calculate_current_coefficient_todd(eps, kappa95, triang95, model=1)

            if i_plasma_current == 6:
                fq = calculate_current_coefficient_todd(eps, kappa95, triang95, model=2)

        # Connor-Hastie asymptotically-correct expression
        elif i_plasma_current == 7:
            fq = calculate_current_coefficient_hastie(
                alphaj,
                alphap,
                b_plasma_toroidal_on_axis,
                triang95,
                eps,
                kappa95,
                pres_plasma_on_axis,
                constants.RMU0,
            )

        # Sauter scaling allowing negative triangularity [FED May 2016]
        # https://doi.org/10.1016/j.fusengdes.2016.04.033.
        elif i_plasma_current == 8:
            # Assumes zero squareness, note takes kappa, delta at separatrix not _95
            fq = calculate_current_coefficient_sauter(eps, kappa, triang)

        # FIESTA ST scaling
        # https://doi.org/10.1016/j.fusengdes.2020.111530.
        elif i_plasma_current == 9:
            fq = calculate_current_coefficient_fiesta(eps, kappa, triang)
        else:
            raise ProcessValueError(f"Invalid value {i_plasma_current=}")

        # Main plasma current calculation using the fq value from the different settings
        if i_plasma_current != 2:
            plasma_current = (
                (2.0 * np.pi / constants.RMU0)
                * rminor**2
                / (rmajor * q95)
                * fq
                * b_plasma_toroidal_on_axis
            )
        # i_plasma_current == 2 case covered above

        # Calculate cyclindrical safety factor from IPDG89
        qstar = (
            5.0e6
            * rminor**2
            / (rmajor * plasma_current / b_plasma_toroidal_on_axis)
            * 0.5
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

        # Calculate the poloidal field generated by the plasma current
        b_plasma_poloidal_average = calculate_poloidal_field(
            i_plasma_current,
            plasma_current,
            q95,
            aspect_ratio,
            eps,
            b_plasma_toroidal_on_axis,
            kappa,
            triang,
            len_plasma_poloidal,
            constants.RMU0,
        )

        return b_plasma_poloidal_average, qstar, plasma_current

    def outtim(self):
        po.oheadr(self.outfile, "Times")

        po.ovarrf(
            self.outfile,
            "Initial charge time for CS from zero current (s)",
            "(t_plant_pulse_coil_precharge)",
            times_variables.t_plant_pulse_coil_precharge,
        )
        po.ovarrf(
            self.outfile,
            "Plasma current ramp-up time (s)",
            "(t_plant_pulse_plasma_current_ramp_up)",
            times_variables.t_plant_pulse_plasma_current_ramp_up,
        )
        po.ovarrf(
            self.outfile,
            "Heating time (s)",
            "(t_plant_pulse_fusion_ramp)",
            times_variables.t_plant_pulse_fusion_ramp,
        )
        po.ovarre(
            self.outfile,
            "Burn time (s)",
            "(t_plant_pulse_burn)",
            times_variables.t_plant_pulse_burn,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Reset time to zero current for CS (s)",
            "(t_plant_pulse_plasma_current_ramp_down)",
            times_variables.t_plant_pulse_plasma_current_ramp_down,
        )
        po.ovarrf(
            self.outfile,
            "Time between pulses (s)",
            "(t_plant_pulse_dwell)",
            times_variables.t_plant_pulse_dwell,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plant cycle time (s)",
            "(t_plant_pulse_total)",
            times_variables.t_plant_pulse_total,
            "OP ",
        )

    def calculate_effective_charge_ionisation_profiles(self):
        """Calculate the effective charge profiles for ionisation calculations."""

        # Calculate the effective charge (zeff) profile across the plasma
        # Returns an array of zeff at each radial point
        zeff_profile = np.zeros_like(self.plasma_profile.teprofile.profile_y)
        for i in range(len(zeff_profile)):
            zeff_profile[i] = 0.0
            for imp in range(impurity_radiation_module.N_IMPURITIES):
                zeff_profile[i] += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation.zav_of_te(
                        imp, np.array([self.plasma_profile.teprofile.profile_y[i]])
                    ).squeeze()
                    ** 2
                )
        physics_variables.n_charge_plasma_effective_profile = zeff_profile

        # Assign the charge profiles of each species
        n_impurities = impurity_radiation_module.N_IMPURITIES
        te_profile = self.plasma_profile.teprofile.profile_y
        n_points = len(te_profile)
        # Create a 2D array: (n_impurities, n_points)
        charge_profiles = np.zeros((n_impurities, n_points))
        for imp in range(n_impurities):
            for i in range(n_points):
                charge_profiles[imp, i] = impurity_radiation.zav_of_te(
                    imp, np.array([te_profile[i]])
                ).squeeze()
        impurity_radiation_module.n_charge_impurity_profile = charge_profiles

    def outplas(self):
        """Subroutine to output the plasma physics information
        self.outfile : input integer : Fortran output unit identifier
        This routine writes the plasma physics information
        to a file, in a tidy format.
        """
        # Dimensionless plasma parameters. See reference below.
        physics_variables.nu_star = (
            1
            / constants.RMU0
            * (15.0e0 * constants.ELECTRON_CHARGE**4 * physics_variables.dlamie)
            / (4.0e0 * np.pi**1.5e0 * constants.EPSILON0**2)
            * physics_variables.vol_plasma**2
            * physics_variables.rmajor**2
            * physics_variables.b_plasma_toroidal_on_axis
            * np.sqrt(physics_variables.eps)
            * physics_variables.nd_plasma_electron_line**3
            * physics_variables.kappa
            / (physics_variables.e_plasma_beta**2 * physics_variables.plasma_current)
        )

        physics_variables.rho_star = np.sqrt(
            2.0e0
            * constants.PROTON_MASS
            * physics_variables.m_ions_total_amu
            * physics_variables.e_plasma_beta
            / (
                3.0e0
                * physics_variables.vol_plasma
                * physics_variables.nd_plasma_electron_line
            )
        ) / (
            constants.ELECTRON_CHARGE
            * physics_variables.b_plasma_toroidal_on_axis
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_variables.beta_mcdonald = (
            4.0e0
            / 3.0e0
            * constants.RMU0
            * physics_variables.e_plasma_beta
            / (
                physics_variables.vol_plasma
                * physics_variables.b_plasma_toroidal_on_axis**2
            )
        )

        po.oheadr(self.outfile, "Plasma")

        if stellarator_variables.istell == 0:
            if divertor_variables.n_divertors == 0:
                po.ocmmnt(self.outfile, "Plasma configuration = limiter")
            elif divertor_variables.n_divertors == 1:
                po.ocmmnt(self.outfile, "Plasma configuration = single null divertor")
            elif divertor_variables.n_divertors == 2:
                po.ocmmnt(self.outfile, "Plasma configuration = double null divertor")
            else:
                raise ProcessValueError(
                    "Illegal value of n_divertors",
                    n_divertors=divertor_variables.n_divertors,
                )
        else:
            po.ocmmnt(self.outfile, "Plasma configuration = stellarator")

        if stellarator_variables.istell == 0:
            if physics_variables.itart == 0:
                physics_variables.itart_r = physics_variables.itart
                po.ovarin(
                    self.outfile,
                    "Tokamak aspect ratio = Conventional, itart = 0",
                    "(itart)",
                    physics_variables.itart_r,
                )
            elif physics_variables.itart == 1:
                physics_variables.itart_r = physics_variables.itart
                po.ovarin(
                    self.outfile,
                    "Tokamak aspect ratio = Spherical, itart = 1",
                    "(itart)",
                    physics_variables.itart_r,
                )

        po.osubhd(self.outfile, "Plasma Geometry :")
        po.ovarin(
            self.outfile,
            "Plasma shaping model",
            "(i_plasma_shape)",
            physics_variables.i_plasma_shape,
        )
        if physics_variables.i_plasma_shape == 0:
            po.osubhd(self.outfile, "Classic PROCESS plasma shape model is used :")
        elif physics_variables.i_plasma_shape == 1:
            po.osubhd(self.outfile, "Sauter plasma shape model is used :")
        po.ovarrf(self.outfile, "Major radius (m)", "(rmajor)", physics_variables.rmajor)
        po.ovarrf(
            self.outfile,
            "Minor radius (m)",
            "(rminor)",
            physics_variables.rminor,
            "OP ",
        )
        po.ovarrf(self.outfile, "Aspect ratio", "(aspect)", physics_variables.aspect)
        po.ovarrf(
            self.outfile,
            "Plasma squareness",
            "(plasma_square)",
            physics_variables.plasma_square,
            "IP",
        )
        if stellarator_variables.istell == 0:
            if physics_variables.i_plasma_geometry in [0, 6, 8]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (input value used)",
                    "(kappa)",
                    physics_variables.kappa,
                    "IP ",
                )
            elif physics_variables.i_plasma_geometry == 1:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (TART scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.i_plasma_geometry in [2, 3]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (Zohm scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Zohm scaling adjustment factor",
                    "(fkzohm)",
                    physics_variables.fkzohm,
                )
            elif physics_variables.i_plasma_geometry in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from kappa95)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.i_plasma_geometry == 9:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and li(3))",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.i_plasma_geometry == 10:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and stability margin)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.i_plasma_geometry == 11:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio via Menard 2016)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            else:
                raise ProcessValueError(
                    "Illegal value of i_plasma_geometry",
                    i_plasma_geometry=physics_variables.i_plasma_geometry,
                )

            if physics_variables.i_plasma_geometry in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (input value used)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (calculated from kappa)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "OP ",
                )

            if physics_variables.i_plasma_geometry in [0, 2, 6, 8, 9, 10, 11]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (input value used)",
                    "(triang)",
                    physics_variables.triang,
                    "IP ",
                )
            elif physics_variables.i_plasma_geometry == 1:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (TART scaling)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (calculated from triang95)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )

            if physics_variables.i_plasma_geometry in [3, 4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (input value used)",
                    "(triang95)",
                    physics_variables.triang95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (calculated from triang)",
                    "(triang95)",
                    physics_variables.triang95,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Plasma poloidal perimeter (m)",
                "(len_plasma_poloidal)",
                physics_variables.len_plasma_poloidal,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Plasma cross-sectional area (m2)",
                "(a_plasma_poloidal)",
                physics_variables.a_plasma_poloidal,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma surface area (m2)",
                "(a_plasma_surface)",
                physics_variables.a_plasma_surface,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma volume (m3)",
                "(vol_plasma)",
                physics_variables.vol_plasma,
                "OP ",
            )

            po.osubhd(self.outfile, "Current and Field :")

            if stellarator_variables.istell == 0:
                po.oblnkl(self.outfile)
                po.ovarin(
                    self.outfile,
                    "Plasma current scaling law used",
                    "(i_plasma_current)",
                    physics_variables.i_plasma_current,
                )

                po.ovarrf(
                    self.outfile,
                    "Plasma current (MA)",
                    "(plasma_current_MA)",
                    physics_variables.plasma_current / 1.0e6,
                    "OP ",
                )

                if physics_variables.i_alphaj == 1:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                        "OP ",
                    )
                else:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                    )
                po.ocmmnt(self.outfile, "Current profile index scalings:")
                po.oblnkl(self.outfile)

                po.ovarrf(
                    self.outfile,
                    "J. Wesson plasma current profile index",
                    "(alphaj_wesson)",
                    physics_variables.alphaj_wesson,
                    "OP ",
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "On-axis plasma current density (A/m2)",
                    "(j_plasma_on_axis)",
                    physics_variables.j_plasma_on_axis,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Vertical field at plasma (T)",
                    "(b_plasma_vertical_required)",
                    physics_variables.b_plasma_vertical_required,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Vacuum toroidal field at R (T)",
                "(b_plasma_toroidal_on_axis)",
                physics_variables.b_plasma_toroidal_on_axis,
            )
            po.ovarrf(
                self.outfile,
                "Toroidal field at plasma inboard (T)",
                "(b_plasma_inboard_toroidal)",
                physics_variables.b_plasma_inboard_toroidal,
            )
            po.ovarrf(
                self.outfile,
                "Toroidal field at plasma outboard (T)",
                "(b_plasma_outboard_toroidal)",
                physics_variables.b_plasma_outboard_toroidal,
            )

            for i in range(len(physics_variables.b_plasma_toroidal_profile)):
                po.ovarre(
                    self.mfile,
                    f"Toroidal field in plasma at point {i}",
                    f"b_plasma_toroidal_profile{i}",
                    physics_variables.b_plasma_toroidal_profile[i],
                )
            po.ovarrf(
                self.outfile,
                "Average poloidal field (T)",
                "(b_plasma_poloidal_average)",
                physics_variables.b_plasma_poloidal_average,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Total field (sqrt(b_plasma_poloidal_average^2 + b_plasma_toroidal_on_axis^2)) (T)",
                "(b_plasma_total)",
                physics_variables.b_plasma_total,
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.ovarrf(
                self.outfile, "Safety factor on axis", "(q0)", physics_variables.q0
            )

            if physics_variables.i_plasma_current == 2:
                po.ovarrf(
                    self.outfile,
                    "Mean edge safety factor",
                    "(q95)",
                    physics_variables.q95,
                )

            po.ovarrf(
                self.outfile,
                "Safety factor at 95% flux surface",
                "(q95)",
                physics_variables.q95,
            )

            po.ovarrf(
                self.outfile,
                "Cylindrical safety factor (qcyl)",
                "(qstar)",
                physics_variables.qstar,
                "OP ",
            )

            if physics_variables.i_plasma_geometry == 1:
                po.ovarrf(
                    self.outfile,
                    "Lower limit for edge safety factor q95",
                    "(q95_min)",
                    physics_variables.q95_min,
                    "OP ",
                )
            po.ovarrf(
                self.outfile,
                "Plasma normalised internal inductance",
                "(ind_plasma_internal_norm)",
                physics_variables.ind_plasma_internal_norm,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Plasma normalised internal inductance scalings:")
            po.oblnkl(self.outfile)
            po.ovarrf(
                self.outfile,
                "J. Wesson plasma normalised internal inductance",
                "(ind_plasma_internal_norm_wesson)",
                physics_variables.ind_plasma_internal_norm_wesson,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "J. Menard plasma normalised internal inductance",
                "(ind_plasma_internal_norm_menard)",
                physics_variables.ind_plasma_internal_norm_menard,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "ITER li(3) plasma normalised internal inductance",
                "(ind_plasma_internal_norm_iter_3)",
                physics_variables.ind_plasma_internal_norm_iter_3,
                "OP ",
            )

        else:
            po.ovarrf(
                self.outfile,
                "Rotational transform",
                "(iotabar)",
                stellarator_variables.iotabar,
            )
        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        # Output beta information
        self.beta.output_beta_information()

        po.osubhd(self.outfile, "Temperature and Density (volume averaged) :")
        po.ovarrf(
            self.outfile,
            "Number of radial points in plasma profiles",
            "(n_plasma_profile_elements)",
            physics_variables.n_plasma_profile_elements,
        )
        po.ovarrf(
            self.outfile,
            "Volume averaged electron temperature (keV)",
            "(temp_plasma_electron_vol_avg_kev)",
            physics_variables.temp_plasma_electron_vol_avg_kev,
        )
        po.ovarrf(
            self.outfile,
            "Ratio of ion to electron volume-averaged temperature",
            "(f_temp_plasma_ion_electron)",
            physics_variables.f_temp_plasma_ion_electron,
            "IP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temperature on axis (keV)",
            "(temp_plasma_electron_on_axis_kev)",
            physics_variables.temp_plasma_electron_on_axis_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ion temperature (keV)",
            "(temp_plasma_ion_vol_avg_kev)",
            physics_variables.temp_plasma_ion_vol_avg_kev,
        )
        po.ovarrf(
            self.outfile,
            "Ion temperature on axis (keV)",
            "(temp_plasma_ion_on_axis_kev)",
            physics_variables.temp_plasma_ion_on_axis_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temp., density weighted (keV)",
            "(temp_plasma_electron_density_weighted_kev)",
            physics_variables.temp_plasma_electron_density_weighted_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ratio of electron density weighted temp. to volume averaged temp.",
            "(f_temp_plasma_electron_density_vol_avg)",
            physics_variables.f_temp_plasma_electron_density_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged electron number density (/m3)",
            "(nd_plasma_electrons_vol_avg)",
            physics_variables.nd_plasma_electrons_vol_avg,
        )
        po.ovarre(
            self.outfile,
            "Electron number density on axis (/m3)",
            "(nd_plasma_electron_on_axis)",
            physics_variables.nd_plasma_electron_on_axis,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line-averaged electron number density (/m3)",
            "(nd_plasma_electron_line)",
            physics_variables.nd_plasma_electron_line,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma pressure on axis (Pa)",
            "(pres_plasma_thermal_on_axis)",
            physics_variables.pres_plasma_thermal_on_axis,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged plasma pressure (Pa)",
            "(pres_plasma_thermal_vol_avg)",
            physics_variables.pres_plasma_thermal_vol_avg,
            "OP ",
        )
        # As array output is not currently supported, each element is output as a float instance
        # Output plasma pressure profiles to mfile
        for i in range(len(physics_variables.pres_plasma_thermal_total_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma pressure at point {i}",
                f"(pres_plasma_thermal_total_profile{i})",
                physics_variables.pres_plasma_thermal_total_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma electron pressure at point {i}",
                f"(pres_plasma_electron_profile{i})",
                physics_variables.pres_plasma_electron_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_ion_total_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma ion pressure at point {i}",
                f"(pres_plasma_ion_total_profile{i})",
                physics_variables.pres_plasma_ion_total_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_fuel_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma fuel pressure at point {i}",
                f"(pres_plasma_fuel_profile{i})",
                physics_variables.pres_plasma_fuel_profile[i],
            )

        if stellarator_variables.istell == 0:
            po.ovarre(
                self.outfile,
                "Line-averaged electron density / Greenwald density",
                "(dnla_gw)",
                physics_variables.nd_plasma_electron_line
                / physics_variables.nd_plasma_electron_max_array[6],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total Ion number density (/m3)",
            "(nd_plasma_ions_total_vol_avg)",
            physics_variables.nd_plasma_ions_total_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel ion number density (/m3)",
            "(nd_plasma_fuel_ions_vol_avg)",
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total impurity number density with Z > 2 (no He) (/m3)",
            "(nd_plasma_impurities_vol_avg)",
            physics_variables.nd_plasma_impurities_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion number density (thermalised ions only) (/m3)",
            "(nd_plasma_alphas_vol_avg)",
            physics_variables.nd_plasma_alphas_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) / electron number density",
            "(f_nd_alpha_electron)",
            physics_variables.f_nd_alpha_electron,
        )
        po.ovarre(
            self.outfile,
            "Plasma volume averaged proton number density (/m3)",
            "(nd_plasma_protons_vol_avg)",
            physics_variables.nd_plasma_protons_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Proton number density / electron number density",
            "(f_nd_protium_electrons)",
            physics_variables.f_nd_protium_electrons,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Hot beam ion number density (/m3)",
            "(nd_beam_ions)",
            physics_variables.nd_beam_ions,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hot beam ion number density / electron density",
            "(f_nd_beam_electron)",
            physics_variables.f_nd_beam_electron,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Impurities:")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Plasma ion densities / electron density:")

        for imp in range(impurity_radiation_module.N_IMPURITIES):
            # MDK Update f_nd_impurity_electrons, as this will make the ITV output work correctly.
            impurity_radiation_module.f_nd_impurity_electrons[imp] = (
                impurity_radiation_module.f_nd_impurity_electron_array[imp]
            )
            str1 = (
                impurity_radiation_module.impurity_arr_label[imp].item()
                + " concentration"
            )
            str2 = f"(f_nd_impurity_electrons({imp + 1:02}))"
            # MDK Add output flag for H which is calculated.
            if imp == 0:
                po.ovarre(
                    self.outfile,
                    str1,
                    str2,
                    impurity_radiation_module.f_nd_impurity_electrons[imp],
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    str1,
                    str2,
                    impurity_radiation_module.f_nd_impurity_electrons[imp],
                )

        po.ovarre(
            self.outfile,
            "Average mass of all ions (amu)",
            "(m_ions_total_amu)",
            physics_variables.m_ions_total_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all ions in plasma (kg)",
            "(m_plasma_ions_total)",
            physics_variables.m_plasma_ions_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Average mass of all fuel ions (amu)",
            "(m_fuel_amu)",
            physics_variables.m_fuel_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all fuel ions in plasma (kg)",
            "(m_plasma_fuel_ions)",
            physics_variables.m_plasma_fuel_ions,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Average mass of all beam ions (amu)",
            "(m_beam_amu)",
            physics_variables.m_beam_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all alpha particles in plasma (kg)",
            "(m_plasma_alpha)",
            physics_variables.m_plasma_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all electrons in plasma (kg)",
            "(m_plasma_electron)",
            physics_variables.m_plasma_electron,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of the plasma (kg)",
            "(m_plasma)",
            physics_variables.m_plasma,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Volume averaged plasma effective charge",
            "(n_charge_plasma_effective_vol_avg)",
            physics_variables.n_charge_plasma_effective_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.n_charge_plasma_effective_profile)):
            po.ovarre(
                self.mfile,
                "Effective charge at point",
                f"(n_charge_plasma_effective_profile{i})",
                physics_variables.n_charge_plasma_effective_profile[i],
                "OP ",
            )

        for imp in range(impurity_radiation_module.N_IMPURITIES):
            for i in range(physics_variables.n_plasma_profile_elements):
                po.ovarre(
                    self.mfile,
                    "Impurity charge at point",
                    f"(n_charge_plasma_profile{imp}_{i})",
                    impurity_radiation_module.n_charge_impurity_profile[imp][i],
                    "OP ",
                )

        po.ovarrf(
            self.outfile,
            "Mass-weighted Effective charge",
            "(n_charge_plasma_effective_mass_weighted_vol_avg)",
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            "OP ",
        )

        po.ovarrf(
            self.outfile, "Density profile factor", "(alphan)", physics_variables.alphan
        )
        po.ovarin(
            self.outfile,
            "Plasma profile model",
            "(i_plasma_pedestal)",
            physics_variables.i_plasma_pedestal,
        )

        if physics_variables.i_plasma_pedestal >= 1:
            if (
                physics_variables.nd_plasma_electron_on_axis
                < physics_variables.nd_plasma_pedestal_electron
            ):
                logger.error("Central density is less than pedestal density")

            po.ocmmnt(self.outfile, "Pedestal profiles are used.")
            po.ovarrf(
                self.outfile,
                "Density pedestal r/a location",
                "(radius_plasma_pedestal_density_norm)",
                physics_variables.radius_plasma_pedestal_density_norm,
            )
            if physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(nd_plasma_pedestal_electron)",
                    physics_variables.nd_plasma_pedestal_electron,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(nd_plasma_pedestal_electron)",
                    physics_variables.nd_plasma_pedestal_electron,
                )

            # must be assigned to their exisiting values#
            fgwped_out = (
                physics_variables.nd_plasma_pedestal_electron
                / physics_variables.nd_plasma_electron_max_array[6]
            )
            fgwsep_out = (
                physics_variables.nd_plasma_separatrix_electron
                / physics_variables.nd_plasma_electron_max_array[6]
            )
            if physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0:
                physics_variables.f_nd_plasma_pedestal_greenwald = (
                    physics_variables.nd_plasma_pedestal_electron
                    / physics_variables.nd_plasma_electron_max_array[6]
                )
            if physics_variables.f_nd_plasma_separatrix_greenwald >= 0e0:
                physics_variables.f_nd_plasma_separatrix_greenwald = (
                    physics_variables.nd_plasma_separatrix_electron
                    / physics_variables.nd_plasma_electron_max_array[6]
                )

            po.ovarre(
                self.outfile,
                "Electron density at pedestal / nGW",
                "(fgwped_out)",
                fgwped_out,
            )
            po.ovarrf(
                self.outfile,
                "Temperature pedestal r/a location",
                "(radius_plasma_pedestal_temp_norm)",
                physics_variables.radius_plasma_pedestal_temp_norm,
            )

            po.ovarrf(
                self.outfile,
                "Electron temp. pedestal height (keV)",
                "(temp_plasma_pedestal_kev)",
                physics_variables.temp_plasma_pedestal_kev,
            )
            if 78 in numerics.icc:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(temp_plasma_separatrix_kev)",
                    physics_variables.temp_plasma_separatrix_kev,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(temp_plasma_separatrix_kev)",
                    physics_variables.temp_plasma_separatrix_kev,
                )

            po.ovarre(
                self.outfile,
                "Electron density at separatrix (/m3)",
                "(nd_plasma_separatrix_electron)",
                physics_variables.nd_plasma_separatrix_electron,
            )
            po.ovarre(
                self.outfile,
                "Electron density at separatrix / nGW",
                "(fgwsep_out)",
                fgwsep_out,
            )

        # Issue 558 - addition of constraint 76 to limit the value of nd_plasma_separatrix_electron, in proportion with the ballooning parameter and Greenwald density
        if 76 in numerics.icc:
            po.ovarre(
                self.outfile,
                "Critical ballooning parameter value",
                "(alpha_crit)",
                physics_variables.alpha_crit,
            )
            po.ovarre(
                self.outfile,
                "Critical electron density at separatrix (/m3)",
                "(nd_plasma_separatrix_electron_eich_max)",
                physics_variables.nd_plasma_separatrix_electron_eich_max,
            )

        po.ovarrf(
            self.outfile,
            "Temperature profile index",
            "(alphat)",
            physics_variables.alphat,
        )
        po.ovarrf(
            self.outfile,
            "Temperature profile index beta",
            "(tbeta)",
            physics_variables.tbeta,
        )
        po.ovarrf(
            self.outfile,
            "Pressure profile index",
            "(alphap)",
            physics_variables.alphap,
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "Density Limit using different models :")
            po.ovarre(
                self.outfile,
                "Old ASDEX model",
                "(nd_plasma_electron_max_array(1))",
                physics_variables.nd_plasma_electron_max_array[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model I",
                "(nd_plasma_electron_max_array(2))",
                physics_variables.nd_plasma_electron_max_array[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model II",
                "(nd_plasma_electron_max_array(3))",
                physics_variables.nd_plasma_electron_max_array[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET edge radiation model",
                "(nd_plasma_electron_max_array(4))",
                physics_variables.nd_plasma_electron_max_array[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET simplified model",
                "(nd_plasma_electron_max_array(5))",
                physics_variables.nd_plasma_electron_max_array[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hugill-Murakami Mq model",
                "(nd_plasma_electron_max_array(6))",
                physics_variables.nd_plasma_electron_max_array[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Greenwald model",
                "(nd_plasma_electron_max_array(7))",
                physics_variables.nd_plasma_electron_max_array[6],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ASDEX New",
                "(nd_plasma_electron_max_array(8))",
                physics_variables.nd_plasma_electron_max_array[7],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Density limit from scaling (/m3)",
                "(nd_plasma_electrons_max)",
                physics_variables.nd_plasma_electrons_max,
                "OP ",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.oheadr(self.outfile, "Plasma Reactions :")

        po.osubhd(self.outfile, "Fuel Constituents :")
        po.ovarrf(
            self.outfile,
            "Deuterium fuel fraction",
            "(f_plasma_fuel_deuterium)",
            physics_variables.f_plasma_fuel_deuterium,
        )
        po.ovarrf(
            self.outfile,
            "Tritium fuel fraction",
            "(f_plasma_fuel_tritium)",
            physics_variables.f_plasma_fuel_tritium,
        )
        po.ovarrf(
            self.outfile,
            "3-Helium fuel fraction",
            "(f_plasma_fuel_helium3)",
            physics_variables.f_plasma_fuel_helium3,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion rates :")
        po.ovarre(
            self.outfile,
            "Fusion rate: total (reactions/sec)",
            "(fusrat_total)",
            physics_variables.fusrat_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion rate density: total (reactions/m3/sec)",
            "(fusden_total)",
            physics_variables.fusden_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion rate density: plasma (reactions/m3/sec)",
            "(fusden_plasma)",
            physics_variables.fusden_plasma,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion Powers :")
        po.ocmmnt(
            self.outfile,
            "Fusion power totals from the main plasma and beam-plasma interactions (if present)\n",
        )

        po.ovarre(
            self.outfile,
            "Total fusion power (MW)",
            "(p_fusion_total_mw)",
            physics_variables.p_fusion_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-T fusion power: total (MW)",
            "(p_dt_total_mw)",
            physics_variables.p_dt_total_mw,
            "OP ",
        )
        for i in range(len(physics_variables.fusrat_plasma_dt_profile)):
            po.ovarre(
                self.mfile,
                f"DT fusion rate at point {i}",
                f"fusrat_plasma_dt_profile{i}",
                physics_variables.fusrat_plasma_dt_profile[i],
            )

        for i in range(len(physics_variables.fusrat_plasma_dd_triton_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> T fusion rate at point {i}",
                f"fusrat_plasma_dd_triton_profile{i}",
                physics_variables.fusrat_plasma_dd_triton_profile[i],
            )
        for i in range(len(physics_variables.fusrat_plasma_dd_helion_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> 3He fusion rate at point {i}",
                f"fusrat_plasma_dd_helion_profile{i}",
                physics_variables.fusrat_plasma_dd_helion_profile[i],
            )
        for i in range(len(physics_variables.fusrat_plasma_dhe3_profile)):
            po.ovarre(
                self.mfile,
                f"D-3He fusion rate at point {i}",
                f"fusrat_plasma_dhe3_profile{i}",
                physics_variables.fusrat_plasma_dhe3_profile[i],
            )
        po.ovarre(
            self.outfile,
            "D-T fusion power: plasma (MW)",
            "(p_plasma_dt_mw)",
            physics_variables.p_plasma_dt_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-T fusion power: beam (MW)",
            "(p_beam_dt_mw)",
            physics_variables.p_beam_dt_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-D fusion power (MW)",
            "(p_dd_total_mw)",
            physics_variables.p_dd_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-D branching ratio for tritium producing reactions",
            "(f_dd_branching_trit)",
            physics_variables.f_dd_branching_trit,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-He3 fusion power (MW)",
            "(p_dhe3_total_mw)",
            physics_variables.p_dhe3_total_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Alpha Powers :")
        po.ovarre(
            self.outfile,
            "Alpha rate density: total (particles/m3/sec)",
            "(fusden_alpha_total)",
            physics_variables.fusden_alpha_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha rate density: plasma (particles/m3/sec)",
            "(fusden_plasma_alpha)",
            physics_variables.fusden_plasma_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: total (MW)",
            "(p_alpha_total_mw)",
            physics_variables.p_alpha_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power density: total (MW/m^3)",
            "(pden_alpha_total_mw)",
            physics_variables.pden_alpha_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: plasma only (MW)",
            "(p_plasma_alpha_mw)",
            physics_variables.p_plasma_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power density: plasma (MW/m^3)",
            "(pden_plasma_alpha_mw)",
            physics_variables.pden_plasma_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: beam-plasma (MW)",
            "(p_beam_alpha_mw)",
            physics_variables.p_beam_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to electrons (MW/m3)",
            "(f_pden_alpha_electron_mw)",
            physics_variables.f_pden_alpha_electron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to ions (MW/m3)",
            "(f_pden_alpha_ions_mw)",
            physics_variables.f_pden_alpha_ions_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Neutron Powers :")
        po.ovarre(
            self.outfile,
            "Neutron power: total (MW)",
            "(p_neutron_total_mw)",
            physics_variables.p_neutron_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power density: total (MW/m^3)",
            "(pden_neutron_total_mw)",
            physics_variables.pden_neutron_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: plasma only (MW)",
            "(p_plasma_neutron_mw)",
            physics_variables.p_plasma_neutron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power density: plasma (MW/m^3)",
            "(pden_plasma_neutron_mw)",
            physics_variables.pden_plasma_neutron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: beam-plasma (MW)",
            "(p_beam_neutron_mw)",
            physics_variables.p_beam_neutron_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")

        po.osubhd(self.outfile, "Charged Particle Powers :")

        po.ovarre(
            self.outfile,
            "Charged particle power (p, 3He, T) (excluding alphas) (MW)",
            "(p_non_alpha_charged_mw)",
            physics_variables.p_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Charged particle power density (p, 3He, T) (excluding alphas) (MW)",
            "(pden_non_alpha_charged_mw)",
            physics_variables.pden_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total charged particle power (including alphas) (MW)",
            "(p_charged_particle_mw)",
            physics_variables.p_charged_particle_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.osubhd(self.outfile, "Plasma radiation powers (excluding SOL):")
        po.ovarre(
            self.outfile,
            "Plasma total synchrotron radiation power (MW)",
            "(p_plasma_sync_mw)",
            physics_variables.p_plasma_sync_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma total synchrotron radiation power density (MW/m^3)",
            "(pden_plasma_sync_mw)",
            physics_variables.pden_plasma_sync_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Synchrotron wall reflectivity factor",
            "(f_sync_reflect)",
            physics_variables.f_sync_reflect,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma normalised minor radius defining 'core' region",
            "(radius_plasma_core_norm)",
            impurity_radiation_module.radius_plasma_core_norm,
        )
        po.ovarre(
            self.outfile,
            "Fractional scaling of radiation power along core profile",
            "(f_p_plasma_core_rad_reduction)",
            impurity_radiation_module.f_p_plasma_core_rad_reduction,
        )
        po.ovarre(
            self.outfile,
            "Plasma total radiation power from core region (MW)",
            "(p_plasma_inner_rad_mw)",
            physics_variables.p_plasma_inner_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma total radiation power from edge region (MW)",
            "(p_plasma_outer_rad_mw)",
            physics_variables.p_plasma_outer_rad_mw,
            "OP ",
        )

        if stellarator_variables.istell != 0:
            po.ovarre(
                self.outfile,
                "SOL radiation power as imposed by f_rad (MW)",
                "(psolradmw)",
                physics_variables.psolradmw,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Plasma total radiation power from inside last closed flux surface (MW)",
            "(p_plasma_rad_mw)",
            physics_variables.p_plasma_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Separatrix radiation fraction (MW)",
            "(f_p_plasma_separatrix_rad)",
            physics_variables.f_p_plasma_separatrix_rad,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean radiation load on vessel first-wall (MW/m^2)",
            "(pflux_fw_rad_mw)",
            physics_variables.pflux_fw_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Average neutron flux at plasma surface (MW/m^2)",
            "(plfux_plasma_surface_neutron_avg_mw)",
            physics_variables.plfux_plasma_surface_neutron_avg_mw,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Peaking factor for radiation first-wall load",
            "(f_fw_rad_max)",
            constraint_variables.f_fw_rad_max,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum permitted radiation first-wall load (MW/m^2)",
            "(pflux_fw_rad_max)",
            constraint_variables.pflux_fw_rad_max,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Peak radiation wall load (MW/m^2)",
            "(pflux_fw_rad_max_mw)",
            constraint_variables.pflux_fw_rad_max_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fast alpha particle power incident on the first-wall (MW)",
            "(p_fw_alpha_mw)",
            physics_variables.p_fw_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean neutron load on vessel first-wall (MW/m^2)",
            "(pflux_fw_neutron_mw)",
            physics_variables.pflux_fw_neutron_mw,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power incident on the divertor targets (MW)",
                "(ptarmw)",
                physics_variables.ptarmw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power to the lower divertor",
                "(f_p_div_lower)",
                physics_variables.f_p_div_lower,
                "IP ",
            )
            po.ovarre(
                self.outfile,
                "Outboard side heat flux decay length (m)",
                "(lambdaio)",
                physics_variables.lambdaio,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Midplane seperation of the two magnetic closed flux surfaces (m)",
                    "(drsep)",
                    physics_variables.drsep,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Fraction of power on the inner targets",
                "(fio)",
                physics_variables.fio,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower inner target",
                "(fLI)",
                physics_variables.fli,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower outer target",
                "(fLO)",
                physics_variables.flo,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper inner target",
                    "(fUI)",
                    physics_variables.fui,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper outer target",
                    "(fUO)",
                    physics_variables.fuo,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Power incident on the lower inner target (MW)",
                "(pLImw)",
                physics_variables.plimw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Power incident on the lower outer target (MW)",
                "(pLOmw)",
                physics_variables.plomw,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper innner target (MW)",
                    "(pUImw)",
                    physics_variables.puimw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper outer target (MW)",
                    "(pUOmw)",
                    physics_variables.puomw,
                    "OP ",
                )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Ohmic heating power (MW)",
            "(p_plasma_ohmic_mw)",
            physics_variables.p_plasma_ohmic_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power deposited in plasma",
            "(f_p_alpha_plasma_deposited)",
            physics_variables.f_p_alpha_plasma_deposited,
            "IP",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to electrons",
            "(f_alpha_electron)",
            physics_variables.f_alpha_electron,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to ions",
            "(f_alpha_ion)",
            physics_variables.f_alpha_ion,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ion transport (MW)",
            "(p_ion_transport_loss_mw)",
            physics_variables.p_ion_transport_loss_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electron transport (MW)",
            "(p_electron_transport_loss_mw)",
            physics_variables.p_electron_transport_loss_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to ions (MW)",
            "(p_hcd_injected_ions_mw)",
            current_drive_variables.p_hcd_injected_ions_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to electrons (MW)",
            "(p_hcd_injected_electrons_mw)",
            current_drive_variables.p_hcd_injected_electrons_mw,
            "OP ",
        )
        if physics_variables.i_plasma_ignited == 1:
            po.ocmmnt(self.outfile, "  (Injected power only used for start-up phase)")

        po.ovarin(
            self.outfile,
            "Ignited plasma switch (0=not ignited, 1=ignited)",
            "(i_plasma_ignited)",
            physics_variables.i_plasma_ignited,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Power into divertor zone via charged particles (MW)",
            "(p_plasma_separatrix_mw)",
            physics_variables.p_plasma_separatrix_mw,
            "OP ",
        )

        if physics_variables.p_plasma_separatrix_mw <= 0.001e0:
            logger.error(
                "Possible problem with high radiation power, forcing p_plasma_separatrix_mw to odd values. "
                f"{physics_variables.p_plasma_separatrix_mw=}"
            )
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(self.outfile, "          Power into divertor zone is unrealistic;")
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if divertor_variables.n_divertors == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Pdivt / R ratio (MW/m) (On peak divertor)",
                "(p_div_separatrix_max_mw/physics_variables.rmajor)",
                physics_variables.p_div_separatrix_max_mw / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)",
                "(pdivmaxbt/qar)",
                (
                    (
                        physics_variables.p_div_separatrix_max_mw
                        * physics_variables.b_plasma_toroidal_on_axis
                    )
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Psep / R ratio (MW/m)",
                "(p_plasma_separatrix_mw/rmajor)",
                physics_variables.p_plasma_separatrix_mw / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Psep Bt / qAR ratio (MWT/m)",
                "(pdivtbt_over_qar)",
                (
                    (
                        physics_variables.p_plasma_separatrix_mw
                        * physics_variables.b_plasma_toroidal_on_axis
                    )
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "H-mode Power Threshold Scalings :")

            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: nominal (MW)",
                "(l_h_threshold_powers(1))",
                physics_variables.l_h_threshold_powers[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: upper bound (MW)",
                "(l_h_threshold_powers(2))",
                physics_variables.l_h_threshold_powers[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: lower bound (MW)",
                "(l_h_threshold_powers(3))",
                physics_variables.l_h_threshold_powers[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (1) (MW)",
                "(l_h_threshold_powers(4))",
                physics_variables.l_h_threshold_powers[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (2) (MW)",
                "(l_h_threshold_powers(5))",
                physics_variables.l_h_threshold_powers[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: nominal (MW)",
                "(l_h_threshold_powers(6))",
                physics_variables.l_h_threshold_powers[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% upper bound (MW)",
                "(l_h_threshold_powers(7))",
                physics_variables.l_h_threshold_powers[6],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% lower bound (MW)",
                "(l_h_threshold_powers(8))",
                physics_variables.l_h_threshold_powers[7],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: nominal (MW)",
                "(l_h_threshold_powers(9))",
                physics_variables.l_h_threshold_powers[8],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: upper bound (MW)",
                "(l_h_threshold_powers(10))",
                physics_variables.l_h_threshold_powers[9],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: lower bound (MW)",
                "(l_h_threshold_powers(11))",
                physics_variables.l_h_threshold_powers[10],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): nominal (MW)",
                "(l_h_threshold_powers(12))",
                physics_variables.l_h_threshold_powers[11],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): upper bound (MW)",
                "(l_h_threshold_powers(13))",
                physics_variables.l_h_threshold_powers[12],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): lower bound (MW)",
                "(l_h_threshold_powers(14))",
                physics_variables.l_h_threshold_powers[13],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - nominal (MW)",
                "(l_h_threshold_powers(15))",
                physics_variables.l_h_threshold_powers[14],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - lower bound (MW)",
                "(l_h_threshold_powers(16))",
                physics_variables.l_h_threshold_powers[15],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - upper bound (MW)",
                "(l_h_threshold_powers(17))",
                physics_variables.l_h_threshold_powers[16],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2017 L-I threshold",
                "(l_h_threshold_powers(18))",
                physics_variables.l_h_threshold_powers[17],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: nominal (MW)",
                "(l_h_threshold_powers(19))",
                physics_variables.l_h_threshold_powers[18],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% upper bound (MW)",
                "(l_h_threshold_powers(20))",
                physics_variables.l_h_threshold_powers[19],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% lower bound (MW)",
                "(l_h_threshold_powers(21))",
                physics_variables.l_h_threshold_powers[20],
                "OP ",
            )
            po.oblnkl(self.outfile)
            if physics_variables.i_l_h_threshold in [9, 10, 11]:
                if (physics_variables.b_plasma_toroidal_on_axis < 0.78e0) or (
                    physics_variables.b_plasma_toroidal_on_axis > 7.94e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(b_plasma_toroidal_on_axis outside Snipes 2000 fitted range)",
                    )
                    logger.warning(
                        "b_plasma_toroidal_on_axis outside Snipes 2000 fitted range"
                    )

                if (physics_variables.rminor < 0.15e0) or (
                    physics_variables.rminor > 1.15e0
                ):
                    po.ocmmnt(self.outfile, "(rminor outside Snipes 2000 fitted range)")
                    logger.warning("rminor outside Snipes 2000 fitted range")

                if (physics_variables.rmajor < 0.55e0) or (
                    physics_variables.rmajor > 3.37e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.rmajor outside Snipes 2000 fitted range)",
                    )
                    logger.warning("rmajor outside Snipes 2000 fitted range")

                if (physics_variables.nd_plasma_electron_line < 0.09e20) or (
                    physics_variables.nd_plasma_electron_line > 3.16e20
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.nd_plasma_electron_line outside Snipes 2000 fitted range)",
                    )
                    logger.warning(
                        "nd_plasma_electron_line outside Snipes 2000 fitted range"
                    )

                if (physics_variables.kappa < 1.0e0) or (
                    physics_variables.kappa > 2.04e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.kappa outside Snipes 2000 fitted range)",
                    )
                    logger.warning("kappa outside Snipes 2000 fitted range")

                if (physics_variables.triang < 0.07e0) or (
                    physics_variables.triang > 0.74e0
                ):
                    po.ocmmnt(self.outfile, "(triang outside Snipes 2000 fitted range)")
                    logger.warning("triang outside Snipes 2000 fitted range")

            po.oblnkl(self.outfile)

            if physics_variables.i_l_h_threshold in [12, 13, 14]:
                po.ocmmnt(
                    self.outfile,
                    "(L-H threshold for closed divertor only. Limited data used in Snipes fit)",
                )
                po.oblnkl(self.outfile)
                logger.warning("Closed divertor only. Limited data used in Snipes fit")

            if (numerics.ioptimz > 0) and (numerics.active_constraints[14]):
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (MW)",
                    "(p_l_h_threshold_mw)",
                    physics_variables.p_l_h_threshold_mw,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (NOT enforced) (MW)",
                    "(p_l_h_threshold_mw)",
                    physics_variables.p_l_h_threshold_mw,
                    "OP ",
                )

        po.osubhd(self.outfile, "Confinement :")

        if physics_variables.i_plasma_ignited == 1:
            po.ocmmnt(
                self.outfile,
                "Device is assumed to be ignited for the calculation of confinement time",
            )
            po.oblnkl(self.outfile)

        tauelaw = physics_variables.LABELS_CONFINEMENT_SCALINGS[
            physics_variables.i_confinement_time
        ]

        po.ocmmnt(
            self.outfile,
            f"Confinement scaling law: {tauelaw}",
        )

        po.ovarst(
            self.outfile,
            "Confinement scaling law",
            "(tauelaw)",
            f'"{tauelaw.strip().split(" ")[0]}"',
        )

        po.ovarrf(
            self.outfile, "Confinement H factor", "(hfact)", physics_variables.hfact
        )
        po.ovarrf(
            self.outfile,
            "Global thermal energy confinement time, from scaling (s)",
            "(t_energy_confinement)",
            physics_variables.t_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Directly calculated total energy confinement time (s)",
            "(t_energy_confinement_beta)",
            physics_variables.t_energy_confinement_beta,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "(Total thermal energy derived from total plasma beta / loss power)",
        )
        po.ovarrf(
            self.outfile,
            "Ion energy confinement time, from scaling (s)",
            "(t_ion_energy_confinement)",
            physics_variables.t_ion_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron energy confinement time, from scaling (s)",
            "(t_electron_energy_confinement)",
            physics_variables.t_electron_energy_confinement,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion double product (s/m3)",
            "(ntau)",
            physics_variables.ntau,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Lawson Triple product (keV s/m3)",
            "(nTtau)",
            physics_variables.nTtau,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Transport loss power assumed in scaling law (MW)",
            "(p_plasma_loss_mw)",
            physics_variables.p_plasma_loss_mw,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "Switch for radiation loss term usage in power balance",
            "(i_rad_loss)",
            physics_variables.i_rad_loss,
        )
        if physics_variables.i_rad_loss == 0:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.p_plasma_rad_mw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is total radiation power)")
        elif physics_variables.i_rad_loss == 1:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.p_plasma_inner_rad_mw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is core radiation power)")
        else:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                0.0e0,
            )
            po.ocmmnt(self.outfile, "  (No radiation correction applied)")
        po.ovarrf(
            self.outfile,
            "H* non-radiation corrected",
            "(hstar)",
            physics_variables.hstar,
            "OP",
        )
        po.ocmmnt(self.outfile, "  (H* assumes IPB98(y,2), ELMy H-mode scaling)")
        po.ovarrf(
            self.outfile,
            "Alpha particle confinement time (s)",
            "(t_alpha_confinement)",
            physics_variables.t_alpha_confinement,
            "OP ",
        )
        # Note alpha confinement time is no longer equal to fuel particle confinement time.
        po.ovarrf(
            self.outfile,
            "Alpha particle/energy confinement time ratio",
            "(f_alpha_energy_confinement)",
            physics_variables.f_alpha_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Lower limit on f_alpha_energy_confinement",
            "(f_alpha_energy_confinement_min)",
            constraint_variables.f_alpha_energy_confinement_min,
        )

        # Plot table of al the H-factor scalings and coparison values
        self.output_confinement_comparison(istell=stellarator_variables.istell)

        if stellarator_variables.istell == 0:
            # Issues 363 Output dimensionless plasma parameters MDK
            po.osubhd(self.outfile, "Dimensionless plasma parameters")
            po.ocmmnt(self.outfile, "For definitions see")
            po.ocmmnt(
                self.outfile,
                "Recent progress on the development and analysis of the ITPA global H-mode confinement database",
            )
            po.ocmmnt(
                self.outfile,
                "D.C. McDonald et al, 2007 Nuclear Fusion v47, 147. (nu_star missing 1/mu0)",
            )
            po.ovarre(
                self.outfile,
                "Normalized plasma pressure beta as defined by McDonald et al",
                "(beta_mcdonald)",
                physics_variables.beta_mcdonald,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized ion Larmor radius",
                "(rho_star)",
                physics_variables.rho_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized collisionality",
                "(nu_star)",
                physics_variables.nu_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER Physics Basis definition of elongation",
                "(kappa_ipb)",
                physics_variables.kappa_ipb,
                "OP ",
            )

            po.oblnkl(self.outfile)
            po.ostars(self.outfile, 110)
            po.oblnkl(self.outfile)

            self.inductance.output_volt_second_information()

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
                po.ocmmnt(
                    self.outfile, "  (User-specified bootstrap current fraction used)"
                )
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
                po.ocmmnt(
                    self.outfile, "  (Wilson bootstrap current fraction model used)"
                )
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
                po.ocmmnt(
                    self.outfile, "  (Diamagnetic current fraction not calculated)"
                )
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
                po.ocmmnt(
                    self.outfile, "  Pfirsch-Schluter current fraction not calculated"
                )
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

        po.osubhd(self.outfile, "Fuelling :")
        po.ovarre(
            self.outfile,
            "Ratio of He and pellet particle confinement times",
            "(tauratio)",
            physics_variables.tauratio,
        )
        po.ovarre(
            self.outfile,
            "Fuelling rate (nucleus-pairs/s)",
            "(molflow_plasma_fuelling_required)",
            physics_variables.molflow_plasma_fuelling_required,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel burn-up rate (reactions/s)",
            "(rndfuel)",
            physics_variables.rndfuel,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Burn-up fraction",
            "(burnup)",
            physics_variables.burnup,
            "OP ",
        )

        if 78 in numerics.icc:
            po.osubhd(self.outfile, "Reinke Criterion :")
            po.ovarin(
                self.outfile,
                "index of impurity to be iterated for divertor detachment",
                "(impvardiv)",
                reinke_variables.impvardiv,
            )
            po.ovarre(
                self.outfile,
                "Minimum Impurity fraction from Reinke",
                "(fzmin)",
                reinke_variables.fzmin,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual Impurity fraction",
                "(fzactual)",
                reinke_variables.fzactual,
            )

    def output_confinement_comparison(self, istell: int):
        """Routine to calculate ignition margin for different confinement scalings and equivalent confinement times for H=1.

        This routine calculates the ignition margin at the final point with different scalings and outputs the results to a file.

        The output includes:
        - Energy confinement times
        - Required H-factors for power balance

        The routine iterates over a range of confinement times, skipping the first user input and a specific index (25). For each confinement time, it calculates various parameters related to confinement and ignition using the `calculate_confinement_time` method. It then calculates the H-factor for when the plasma is ignited using the `find_other_h_factors` method and writes the results to the output file.

        Output format:
        - Header: "Energy confinement times, and required H-factors :"
        - Columns: "Scaling law", "Confinement time [s]", "H-factor for power balance"

        Methods used:
        - `calculate_confinement_time`: Calculates confinement-related parameters.
        - `find_other_h_factors`: Calculates the H-factor for a given confinement time.

        Parameters
        ----------
        istell :
            Indicator for stellarator (0 for tokamak, >=1 for stellarator).

        """

        po.oheadr(self.outfile, "Energy confinement times, and required H-factors :")
        po.ocmmnt(
            self.outfile,
            f"{'':>2}{'Scaling law':<27}{'Electron confinement time [s]':<32}Equivalent H-factor for",
        )
        po.ocmmnt(
            self.outfile,
            f"{'':>38}{'for H = 1':<23}same confinement time",
        )
        po.oblnkl(self.outfile)

        # List of key values for stellarator scalings
        stellarator_scalings = [21, 22, 23, 37, 38]

        # Plot all of the confinement scalings for comparison when H = 1
        # Start from range 1 as the first i_confinement_time is a user input
        # If stellarator, use the stellarator scalings
        for i_confinement_time in (
            range(1, physics_variables.N_CONFINEMENT_SCALINGS)
            if istell == 0
            else stellarator_scalings
        ):
            if i_confinement_time == 25:
                continue
            (
                _,
                _,
                taueez,
                _,
                _,
                _,
            ) = self.calculate_confinement_time(
                physics_variables.m_fuel_amu,
                physics_variables.p_alpha_total_mw,
                physics_variables.aspect,
                physics_variables.b_plasma_toroidal_on_axis,
                physics_variables.nd_plasma_ions_total_vol_avg,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_electron_line,
                physics_variables.eps,
                1.0,
                i_confinement_time,
                physics_variables.i_plasma_ignited,
                physics_variables.kappa,
                physics_variables.kappa95,
                physics_variables.p_non_alpha_charged_mw,
                current_drive_variables.p_hcd_injected_total_mw,
                physics_variables.plasma_current,
                physics_variables.pden_plasma_core_rad_mw,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.temp_plasma_electron_density_weighted_kev,
                physics_variables.temp_plasma_ion_density_weighted_kev,
                physics_variables.q95,
                physics_variables.qstar,
                physics_variables.vol_plasma,
                physics_variables.n_charge_plasma_effective_vol_avg,
            )

            try:
                # Calculate the H-factor for the same confinement time in other scalings
                physics_variables.hfac[i_confinement_time - 1] = (
                    self.find_other_h_factors(i_confinement_time)
                )
            except ValueError:
                # This is only used for a table in the OUT.DAT so if it fails
                # just write a NaN--its not worth crashing PROCESS over.
                physics_variables.hfac[i_confinement_time - 1] = np.nan

            po.ocmmnt(
                self.outfile,
                f"{'':>2}{physics_variables.LABELS_CONFINEMENT_SCALINGS[i_confinement_time]:<38}"
                f"{taueez:<28.3f}{physics_variables.hfac[i_confinement_time - 1]:.3f}",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)

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
        """Calculate the bootstrap-driven fraction of the plasma current.

        This function performs the original ITER calculation of the plasma current bootstrap fraction.

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
            Plasma volume (m3).

        Returns
        -------
        float
            The bootstrap-driven fraction of the plasma current.

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
        """Bootstrap current fraction from Wilson et al scaling

        This function calculates the bootstrap current fraction using the numerically fitted algorithm written by Howard Wilson.

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

        Reference:
        - AEA FUS 172
            Physics Assessment for the European Reactor Study, 1989
        - AEA FUS 172
            Physics Assessment for the European Reactor Study, 1989
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
    def bootstrap_fraction_nevins(
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

        This function calculates the bootstrap current fraction using the Nevins et al method.

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
            Electron density (/m3).
        plasma_current : float
            Plasma current (A).
        q0 : float
            Central safety factor.
        q95 : float
            Safety factor at 95% surface.
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
            lambda y: _nevins_integral(
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
    def bootstrap_fraction_sauter(plasma_profile: float) -> float:
        """Calculate the bootstrap current fraction from the Sauter et al scaling.

        This function calculates the bootstrap current fraction using the Sauter, Angioni, and Lin-Liu scaling.

        Parameters
        ----------
        plasma_profile : PlasmaProfile
            The plasma profile object.

        Returns
        -------
        float
            The bootstrap current fraction.

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
                _calculate_l31_coefficient(
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
                + _calculate_l31_32_coefficient(
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
                + _calculate_l34_alpha_31_coefficient(
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
        beta_poloidal :
            Plasma poloidal beta.
        q95 :
            Safety factor at 95% of the plasma radius.
        q0 :
            Safety factor at the magnetic axis.
        alphan :
            Density profile index
        alphat :
            Temperature profile index
        eps :
            Inverse aspect ratio.
        ind_plasma_internal_norm :
            Plasma normalised internal inductance

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            The profile assumed for the alphan and alpat indexes is only a parabolic profile without a pedestal (L-mode).
            The Root Mean Squared Error for the fitting database of this formula was 0.025
            Concentrating on the positive shear plasmas using the ACCOME code equilibria with the fully non-inductively driven
            conditions with neutral beam (NB) injection only are calculated.
            The electron temperature and the ion temperature were assumed to be equal
            This can be used for all apsect ratios.
            The diamagnetic fraction is included in this formula.

        References
        ----------
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
        """Calculate the bootstrap fraction using the ARIES formula.
        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        ind_plasma_internal_norm :
            Plasma normalized internal inductance.
        core_density :
            Core plasma density.
        average_density :
            Average plasma density.
        inverse_aspect :
            Inverse aspect ratio.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            - The source reference does not provide any info about the derivation of the formula. It is only stated

        References
        ----------
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
        """Calculate the bootstrap fraction using the Andrade et al formula.

        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        core_pressure :
            Core plasma pressure.
        average_pressure :
            Average plasma pressure.
        inverse_aspect :
            Inverse aspect ratio.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            - Based off 350 plasma profiles from Experimento Tokamak Esferico (ETE) spherical tokamak
            - A = 1.5, R_0 = 0.3m, I_p = 200kA, B_0=0.4T, beta = 4-10%. Profiles taken as Gaussian shaped functions.
            - Errors mostly up to the order of 10% are obtained when both expressions are compared with the equilibrium estimates for the
            bootstrap current in ETE

        References
        ----------
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
        """Calculate the bootstrap fraction using the Hoang et al formula.
        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        pressure_index :
            Pressure profile index.
        current_index :
            Current profile index.
        inverse_aspect :
            Inverse aspect ratio.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            - Based off of TFTR data calculated using the TRANSP plasma analysis code
            - 170 discharges which was assembled to  study the tritium influx and transport in discharges with D-only neutral beam
            injection (NBI)
            - Contains L-mode, supershots, reversed shear, enhanced reversed shear and increased li discharges
            - Discharges with monotonic flux profiles with reversed shear are also included
            - Is applied to circular cross-section plasmas

        References
        ----------
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
        """Calculate the bootstrap fraction using the Wong et al formula.
        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        density_index :
            Density profile index.
        temperature_index :
            Temperature profile index.
        inverse_aspect :
            Inverse aspect ratio.
        elongation :
            Plasma elongation.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            - Data is based off of equilibria from Miller et al.
            - A: 1.2 - 3.0 and stable to n ballooning and low n kink modes at a bootstrap fraction of 99% for kappa = 2, 2.5 and 3
            - The results were parameterized as a function of aspect ratio and elongation
            - The parametric dependency of beta_p and beta_T are based on fitting of the DIII-D high equivalent DT yield results
            - Parabolic profiles should be used for best results as the pressure peaking value is calculated as the product of a parabolic
            temperature and density profile

        References
        ----------
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
        """Calculate the bootstrap fraction using the first scaling from the Gi et al formula.

        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        pressure_index :
            Pressure profile index.
        temperature_index :
            Temperature profile index.
        inverse_aspect :
            Inverse aspect ratio.
        effective_charge :
            Plasma effective charge.
        q95 :
            Safety factor at 95% of the plasma radius.
        q0 :
            Safety factor at the magnetic axis.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
        - Scaling found by solving the Hirshman-Sigmar bootstrap modelusing the matrix inversion method
        - Method was done to put the scaling into parameters compatible with the TPC systems code
        - Uses the ACCOME code to create bootstrap current fractions without using the itrative calculations of the
        current drive and equilibrium models in the scan
        - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, traing = 0.3, alpha_n = 0.1 - 0.8, alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0
        - Uses parabolic plasma profiles only.
        - Scaling 1 has better accuracy than Scaling 2. However, Scaling 1 overestimated the f_BS value for reversed shear
        equilibrium.

        References
        ----------
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
        """Calculate the bootstrap fraction using the second scaling from the Gi et al formula.

        Parameters
        ----------
        beta_poloidal :
            Plasma poloidal beta.
        pressure_index :
            Pressure profile index.
        temperature_index :
            Temperature profile index.
        inverse_aspect :
            Inverse aspect ratio.
        effective_charge :
            Plasma effective charge.

        Returns
        -------
        :
            The calculated bootstrap fraction.

        Notes
        -----
            - Scaling found by solving the Hirshman-Sigmar bootstrap modelusing the matrix inversion method
            - Method was done to put the scaling into parameters compatible with the TPC systems code
            - Uses the ACCOME code to create bootstrap current fractions without using the itrative calculations of the
            curent drive and equilibrium models in the scan
            - R = 5.0 m, A = 1.3 - 5.0, kappa = 2, traing = 0.3, alpha_n = 0.1 - 0.8, alpha_t = 1.0 - 3.0, Z_eff = 1.2 - 3.0
            - Uses parabolic plasma profiles only.
            - This scaling has the q profile dependance removed to obtain a scaling formula with much more flexible variables than
            that by a single profile factor for internal current profile.

            References
            ----------
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
        """Calculate the bootstrap fraction using the L-mode scaling from the Sugiyama et al formula.

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
            - Ion and electron temperature are the same
            - Z_eff has a uniform profile, with only fully stripped carbon impurity

        References
        ----------
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
        """Calculate the bootstrap fraction using the H-mode scaling from the Sugiyama et al formula.

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
            Electron number density at the pedestal [m^-3].
        n_greenwald : float
            Greenwald density limit [m^-3].
        temp_plasma_pedestal_kev : float
            Electron temperature at the pedestal [keV].

        Returns
        -------
        float
            The calculated bootstrap fraction.

        Notes
        -----
            - This scaling is derived for H-mode plasmas.
            - The temperature and density pedestal positions are the same
            - Separatrix temperature and density are zero
            - Ion and electron temperature are the same
            - Z_eff has a uniform profile, with only fully stripped carbon impurity

        References
        ----------
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

    def find_other_h_factors(self, i_confinement_time: int) -> float:
        """Function to find H-factor for the equivalent confinement time in other scalings.

        Parameters
        ----------
        i_confinement_time : int
            Index of the confinement time scaling to use.

        Returns
        -------
        float
            The calculated H-factor.

        """

        def fhz(hfact: float) -> float:
            """Function used to find power balance.

            Parameters
            ----------
            hfact : float
                H-factor to be used in the calculation.
            hfact: float :


            Returns
            -------
            float
                The difference between the calculated power and the required power for balance.

            """
            (
                ptrez,
                ptriz,
                _,
                _,
                _,
                _,
            ) = self.calculate_confinement_time(
                physics_variables.m_fuel_amu,
                physics_variables.p_alpha_total_mw,
                physics_variables.aspect,
                physics_variables.b_plasma_toroidal_on_axis,
                physics_variables.nd_plasma_ions_total_vol_avg,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_electron_line,
                physics_variables.eps,
                hfact,
                i_confinement_time,
                physics_variables.i_plasma_ignited,
                physics_variables.kappa,
                physics_variables.kappa95,
                physics_variables.p_non_alpha_charged_mw,
                current_drive_variables.p_hcd_injected_total_mw,
                physics_variables.plasma_current,
                physics_variables.pden_plasma_core_rad_mw,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.temp_plasma_electron_density_weighted_kev,
                physics_variables.temp_plasma_ion_density_weighted_kev,
                physics_variables.q95,
                physics_variables.qstar,
                physics_variables.vol_plasma,
                physics_variables.n_charge_plasma_effective_vol_avg,
            )

            # At power balance, fhz is zero.
            fhz_value = (
                ptrez
                + ptriz
                - physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_non_alpha_charged_mw
                - physics_variables.pden_plasma_ohmic_mw
            )

            # Take into account whether injected power is included in tau_e calculation (i.e. whether device is ignited)
            if physics_variables.i_plasma_ignited == 0:
                fhz_value -= (
                    current_drive_variables.p_hcd_injected_total_mw
                    / physics_variables.vol_plasma
                )

            # Include the radiation power if requested
            if physics_variables.i_rad_loss == 0:
                fhz_value += physics_variables.pden_plasma_rad_mw
            elif physics_variables.i_rad_loss == 1:
                fhz_value += physics_variables.pden_plasma_core_rad_mw

            return fhz_value

        return root_scalar(fhz, bracket=(0.01, 150), xtol=0.001).root

    @staticmethod
    def calculate_confinement_time(
        m_fuel_amu: float,
        p_alpha_total_mw: float,
        aspect: float,
        b_plasma_toroidal_on_axis: float,
        nd_plasma_ions_total_vol_avg: float,
        nd_plasma_electrons_vol_avg: float,
        nd_plasma_electron_line: float,
        eps: float,
        hfact: float,
        i_confinement_time: int,
        i_plasma_ignited: int,
        kappa: float,
        kappa95: float,
        p_non_alpha_charged_mw: float,
        p_hcd_injected_total_mw: float,
        plasma_current: float,
        pden_plasma_core_rad_mw: float,
        rmajor: float,
        rminor: float,
        temp_plasma_electron_density_weighted_kev: float,
        temp_plasma_ion_density_weighted_kev: float,
        q95: float,
        qstar: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate the confinement times and the transport power loss terms.

        Parameters
        ----------
        m_fuel_amu :
            Average mass of fuel (amu)
        p_alpha_total_mw :
            Alpha particle power (MW)
        aspect :
            Aspect ratio
        b_plasma_toroidal_on_axis :
            Toroidal field on axis (T)
        nd_plasma_ions_total_vol_avg :
            Total ion density (/m3)
        nd_plasma_electrons_vol_avg :
            Volume averaged electron density (/m3)
        nd_plasma_electron_line :
            Line-averaged electron density (/m3)
        eps :
            Inverse aspect ratio
        hfact :
            H factor on energy confinement scalings
        i_confinement_time :
            Switch for energy confinement scaling to use
        i_plasma_ignited :
            Switch for ignited calculation
        kappa :
            Plasma elongation
        kappa95 :
            Plasma elongation at 95% surface
        p_non_alpha_charged_mw :
            Non-alpha charged particle fusion power (MW)
        p_hcd_injected_total_mw :
            Auxiliary power to ions and electrons (MW)
        plasma_current :
            Plasma current (A)
        pden_plasma_core_rad_mw :
            Total core radiation power (MW/m3)
        q95 :
            Edge safety factor (tokamaks), or rotational transform iotabar (stellarators)
        qstar :
            Equivalent cylindrical edge safety factor
        rmajor :
            Plasma major radius (m)
        rminor :
            Plasma minor radius (m)
        temp_plasma_electron_density_weighted_kev :
            Density weighted average electron temperature (keV)
        temp_plasma_ion_density_weighted_kev :
            Density weighted average ion temperature (keV)
        vol_plasma :
            Plasma volume (m3)
        a_plasma_poloidal :
            Plasma cross-sectional area (m2)
        zeff :
            Plasma effective charge

        Returns
        -------
        type
            Tuple containing:
            - pden_electron_transport_loss_mw (float): Electron transport power (MW/m3)
            - pden_ion_transport_loss_mw (float): Ion transport power (MW/m3)
            - t_electron_energy_confinement (float): Electron energy confinement time (s)
            - t_ion_energy_confinement (float): Ion energy confinement time (s)
            - t_energy_confinement (float): Global energy confinement time (s)
            - p_plasma_loss_mw (float): Heating power (MW) assumed in calculation

        """

        # ========================================================================

        # Calculate heating power (MW)
        p_plasma_loss_mw = (
            physics_variables.f_p_alpha_plasma_deposited * p_alpha_total_mw
            + p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
        )

        # If the device is not ignited, add the injected auxiliary power
        if i_plasma_ignited == 0:
            p_plasma_loss_mw = p_plasma_loss_mw + p_hcd_injected_total_mw

        # Include the radiation as a loss term if requested
        if physics_variables.i_rad_loss == 0:
            p_plasma_loss_mw = (
                p_plasma_loss_mw - physics_variables.pden_plasma_rad_mw * vol_plasma
            )
        elif physics_variables.i_rad_loss == 1:
            p_plasma_loss_mw = (
                p_plasma_loss_mw - pden_plasma_core_rad_mw * vol_plasma
            )  # shouldn't this be vol_core instead of vol_plasma?
        # else do not adjust p_plasma_loss_mw for radiation

        # Ensure heating power is positive (shouldn't be necessary)
        p_plasma_loss_mw = max(p_plasma_loss_mw, 1.0e-3)

        # ========================================================================

        # Line averaged electron density in scaled units
        dnla20 = nd_plasma_electron_line * 1.0e-20
        dnla19 = nd_plasma_electron_line * 1.0e-19

        # Volume averaged electron density in units of 10**20 m**-3
        n20 = nd_plasma_electrons_vol_avg / 1.0e20

        # Plasma current in MA
        pcur = plasma_current / 1.0e6

        # Separatrix kappa defined with plasma volume for IPB scalings
        # Updated version of kappa used by the IPB98 scalings correction in:

        # None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy,
        # “Corrections to a sequence of papers in Nuclear Fusion,” Nuclear Fusion,
        # vol. 48, no. 9, pp. 099801099801, Aug. 2008,
        # doi: https://doi.org/10.1088/0029-5515/48/9/099801.

        physics_variables.kappa_ipb = (vol_plasma / (2.0 * np.pi * rmajor)) / (
            np.pi * rminor**2
        )

        # Electron energy confinement times

        # ========================================================================

        # User defined confinement time
        if i_confinement_time == 0:  # t_electron_energy_confinement is an input
            t_electron_confinement = physics_variables.tauee_in

        # ========================================================================

        # Nec-Alcator(NA) OH scaling
        if i_confinement_time == 1:
            t_electron_confinement = confinement.neo_alcator_confinement_time(
                n20, rminor, rmajor, qstar
            )

        # ========================================================================

        # "Mirnov"-like scaling (H-mode)
        elif i_confinement_time == 2:  # Mirnov scaling (H-mode)
            t_electron_confinement = confinement.mirnov_confinement_time(
                rminor, kappa95, pcur
            )

        # ========================================================================

        # Merezhkin-Mukhovatov (MM) OH/L-mode scaling
        elif i_confinement_time == 3:
            t_electron_confinement = confinement.merezhkin_muhkovatov_confinement_time(
                rmajor,
                rminor,
                kappa95,
                qstar,
                dnla20,
                m_fuel_amu,
                temp_plasma_electron_density_weighted_kev,
            )

        # ========================================================================

        # Shimomura (S) optimized H-mode scaling
        elif i_confinement_time == 4:
            t_electron_confinement = confinement.shimomura_confinement_time(
                rmajor, rminor, b_plasma_toroidal_on_axis, kappa95, m_fuel_amu
            )

        # ========================================================================

        # Kaye-Goldston scaling (L-mode)
        elif i_confinement_time == 5:
            t_electron_confinement = confinement.kaye_goldston_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER Power scaling - ITER 89-P (L-mode)
        elif i_confinement_time == 6:
            t_electron_confinement = confinement.iter_89p_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER Offset linear scaling - ITER 89-O (L-mode)
        elif i_confinement_time == 7:
            t_electron_confinement = confinement.iter_89_0_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )
        # ========================================================================

        # Rebut-Lallia offset linear scaling (L-mode)
        elif i_confinement_time == 8:
            t_electron_confinement = confinement.rebut_lallia_confinement_time(
                rminor,
                rmajor,
                kappa,
                m_fuel_amu,
                pcur,
                zeff,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Goldston scaling (L-mode)
        elif i_confinement_time == 9:  # Goldston scaling (L-mode)
            t_electron_confinement = confinement.goldston_confinement_time(
                pcur, rmajor, rminor, kappa95, m_fuel_amu, p_plasma_loss_mw
            )

        # ========================================================================

        # T-10 scaling (L-mode)
        elif i_confinement_time == 10:
            t_electron_confinement = confinement.t10_confinement_time(
                dnla20,
                rmajor,
                qstar,
                b_plasma_toroidal_on_axis,
                rminor,
                kappa95,
                p_plasma_loss_mw,
                zeff,
                pcur,
            )

        # ========================================================================

        # JAERI / Odajima-Shimomura L-mode scaling
        elif i_confinement_time == 11:  # JAERI scaling
            t_electron_confinement = confinement.jaeri_confinement_time(
                kappa95,
                rminor,
                m_fuel_amu,
                n20,
                pcur,
                b_plasma_toroidal_on_axis,
                rmajor,
                qstar,
                zeff,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Kaye "big"  L-mode scaling (based only on big tokamak data)
        elif i_confinement_time == 12:
            t_electron_confinement = confinement.kaye_big_confinement_time(
                rmajor,
                rminor,
                b_plasma_toroidal_on_axis,
                kappa95,
                pcur,
                n20,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER H90-P H-mode scaling
        elif i_confinement_time == 13:
            t_electron_confinement = confinement.iter_h90_p_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Minimum of ITER 89-P and ITER 89-O
        elif i_confinement_time == 14:
            t_electron_confinement = min(
                confinement.iter_89p_confinement_time(
                    pcur,
                    rmajor,
                    rminor,
                    kappa,
                    dnla20,
                    b_plasma_toroidal_on_axis,
                    m_fuel_amu,
                    p_plasma_loss_mw,
                ),
                confinement.iter_89_0_confinement_time(
                    pcur,
                    rmajor,
                    rminor,
                    kappa,
                    dnla20,
                    b_plasma_toroidal_on_axis,
                    m_fuel_amu,
                    p_plasma_loss_mw,
                ),
            )

        # ========================================================================

        # Riedel scaling (L-mode)
        elif i_confinement_time == 15:
            t_electron_confinement = confinement.riedel_l_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Christiansen et al scaling (L-mode)
        elif i_confinement_time == 16:
            t_electron_confinement = confinement.christiansen_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                m_fuel_amu,
            )

        # ========================================================================

        # Lackner-Gottardi scaling (L-mode)
        elif i_confinement_time == 17:
            t_electron_confinement = confinement.lackner_gottardi_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Neo-Kaye scaling (L-mode)
        elif i_confinement_time == 18:
            t_electron_confinement = confinement.neo_kaye_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ======== ================================================================

        # Riedel scaling (H-mode)
        elif i_confinement_time == 19:
            t_electron_confinement = confinement.riedel_h_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Amended version of ITER H90-P law
        elif i_confinement_time == 20:
            t_electron_confinement = confinement.iter_h90_p_amended_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                rmajor,
                p_plasma_loss_mw,
                kappa,
            )

        # ==========================================================================

        # Sudo et al. scaling (stellarators/heliotron)
        elif i_confinement_time == 21:
            t_electron_confinement = confinement.sudo_et_al_confinement_time(
                rmajor,
                rminor,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Gyro-reduced Bohm scaling
        elif i_confinement_time == 22:
            t_electron_confinement = confinement.gyro_reduced_bohm_confinement_time(
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
                rminor,
                rmajor,
            )

        # ==========================================================================

        # Lackner-Gottardi stellarator scaling
        elif i_confinement_time == 23:
            t_electron_confinement = (
                confinement.lackner_gottardi_stellarator_confinement_time(
                    rmajor,
                    rminor,
                    dnla20,
                    b_plasma_toroidal_on_axis,
                    p_plasma_loss_mw,
                    q95,
                )
            )

        # ==========================================================================

        # ITER_93 ELM-free H-mode scaling
        elif i_confinement_time == 24:
            t_electron_confinement = confinement.iter_93h_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                m_fuel_amu,
                rmajor,
                dnla20,
                aspect,
                kappa,
            )

        # ==========================================================================
        # Scaling removed
        elif i_confinement_time == 25:
            raise ProcessValueError("Scaling removed")
        # ==========================================================================

        # ELM-free: ITERH-97P
        elif i_confinement_time == 26:
            t_electron_confinement = confinement.iter_h97p_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                rmajor,
                aspect,
                kappa,
                m_fuel_amu,
            )

        # ==========================================================================

        # ELMy: ITERH-97P(y)
        elif i_confinement_time == 27:
            t_electron_confinement = confinement.iter_h97p_elmy_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                rmajor,
                aspect,
                kappa,
                m_fuel_amu,
            )

        # ==========================================================================

        # ITER-96P (= ITER-97L) L-mode scaling
        elif i_confinement_time == 28:
            t_electron_confinement = confinement.iter_96p_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                kappa95,
                rmajor,
                aspect,
                dnla19,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Valovic modified ELMy-H mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif i_confinement_time == 29:
            t_electron_confinement = confinement.valovic_elmy_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                m_fuel_amu,
                rmajor,
                rminor,
                kappa,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Kaye PPPL Workshop April 1998 L-mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif i_confinement_time == 30:
            t_electron_confinement = confinement.kaye_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                kappa,
                rmajor,
                aspect,
                dnla19,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # ITERH-PB98P(y), ELMy H-mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif i_confinement_time == 31:
            t_electron_confinement = confinement.iter_pb98py_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y), ELMy H-mode scaling
        elif i_confinement_time == 32:
            t_electron_confinement = confinement.iter_ipb98y_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,1), ELMy H-mode scaling
        elif i_confinement_time == 33:
            t_electron_confinement = confinement.iter_ipb98y1_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,2), ELMy H-mode scaling
        elif i_confinement_time == 34:
            t_electron_confinement = confinement.iter_ipb98y2_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,3), ELMy H-mode scaling
        elif i_confinement_time == 35:
            t_electron_confinement = confinement.iter_ipb98y3_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,4), ELMy H-mode scaling
        elif i_confinement_time == 36:
            t_electron_confinement = confinement.iter_ipb98y4_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # ISS95 stellarator scaling
        elif i_confinement_time == 37:
            # dummy argument q95 is actual argument iotabar for stellarators
            iotabar = q95
            t_electron_confinement = confinement.iss95_stellarator_confinement_time(
                rminor,
                rmajor,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                iotabar,
            )

        # ==========================================================================

        # ISS04 stellarator scaling
        elif i_confinement_time == 38:
            # dummy argument q95 is actual argument iotabar for stellarators
            iotabar = q95
            t_electron_confinement = confinement.iss04_stellarator_confinement_time(
                rminor,
                rmajor,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                iotabar,
            )

        # ==========================================================================

        # DS03 beta-independent H-mode scaling
        elif i_confinement_time == 39:
            t_electron_confinement = confinement.ds03_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa95,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        #  Murari "Non-power law" scaling
        elif i_confinement_time == 40:
            t_electron_confinement = confinement.murari_confinement_time(
                pcur,
                rmajor,
                physics_variables.kappa_ipb,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Petty08, beta independent dimensionless scaling
        elif i_confinement_time == 41:
            t_electron_confinement = confinement.petty08_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
            )

        # ==========================================================================

        # Lang high density relevant confinement scaling
        elif i_confinement_time == 42:
            t_electron_confinement = confinement.lang_high_density_confinement_time(
                plasma_current,
                b_plasma_toroidal_on_axis,
                nd_plasma_electron_line,
                p_plasma_loss_mw,
                rmajor,
                rminor,
                q95,
                qstar,
                aspect,
                m_fuel_amu,
                physics_variables.kappa_ipb,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - nominal
        elif i_confinement_time == 43:
            t_electron_confinement = confinement.hubbard_nominal_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - lower
        elif i_confinement_time == 44:
            t_electron_confinement = confinement.hubbard_lower_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - upper
        elif i_confinement_time == 45:
            t_electron_confinement = confinement.hubbard_upper_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Menard NSTX, ELMy H-mode scaling
        elif i_confinement_time == 46:
            t_electron_confinement = confinement.menard_nstx_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # Menard NSTX-Petty08 Hybrid
        elif i_confinement_time == 47:
            t_electron_confinement = (
                confinement.menard_nstx_petty08_hybrid_confinement_time(
                    pcur,
                    b_plasma_toroidal_on_axis,
                    dnla19,
                    p_plasma_loss_mw,
                    rmajor,
                    physics_variables.kappa_ipb,
                    aspect,
                    m_fuel_amu,
                )
            )

        # ==========================================================================

        # NSTX gyro-Bohm (Buxton)
        elif i_confinement_time == 48:
            t_electron_confinement = confinement.nstx_gyro_bohm_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                rmajor,
                dnla20,
            )

        # ==========================================================================

        # ITPA20 H-mode scaling
        elif i_confinement_time == 49:
            t_electron_confinement = confinement.itpa20_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.triang,
                physics_variables.kappa_ipb,
                eps,
                physics_variables.m_ions_total_amu,
            )

        # ==========================================================================

        # ITPA20-IL confinement time scaling
        elif i_confinement_time == 50:
            t_electron_confinement = confinement.itpa20_il_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                physics_variables.m_ions_total_amu,
                rmajor,
                physics_variables.triang,
                physics_variables.kappa_ipb,
            )

        # ==========================================================================

        else:
            raise ProcessValueError(
                "Illegal value for i_confinement_time",
                i_confinement_time=i_confinement_time,
            )

        # Apply H-factor correction to chosen scaling
        t_electron_energy_confinement = hfact * t_electron_confinement

        # Ion energy confinement time
        t_ion_energy_confinement = t_electron_energy_confinement

        # Calculate H* non-radiation corrected H factor
        # Note: we will assume the IPB-98y2 scaling.
        if physics_variables.i_rad_loss == 1:
            physics_variables.hstar = (
                hfact
                * (
                    p_plasma_loss_mw
                    / (
                        p_plasma_loss_mw
                        + physics_variables.pden_plasma_sync_mw * vol_plasma
                        + physics_variables.p_plasma_inner_rad_mw
                    )
                )
                ** 0.31
            )
        elif physics_variables.i_rad_loss == 0:
            physics_variables.hstar = (
                hfact
                * (
                    p_plasma_loss_mw
                    / (
                        p_plasma_loss_mw
                        + physics_variables.pden_plasma_rad_mw * vol_plasma
                    )
                )
                ** 0.31
            )
        elif physics_variables.i_rad_loss == 2:
            physics_variables.hstar = hfact

        # Calculation of the transport power loss terms
        # Transport losses in Watts/m3 are 3/2 * n.e.T / tau , with T in eV
        # (here, temp_plasma_ion_density_weighted_kev and temp_plasma_electron_density_weighted_kev are in keV, and pden_electron_transport_loss_mw and pden_ion_transport_loss_mw are in MW/m3)

        # The transport losses is just the electron and ion thermal energies divided by the confinement time.
        pden_ion_transport_loss_mw = (
            (3 / 2)
            * (constants.ELECTRON_CHARGE / 1e3)
            * nd_plasma_ions_total_vol_avg
            * temp_plasma_ion_density_weighted_kev
            / t_ion_energy_confinement
        )
        pden_electron_transport_loss_mw = (
            (3 / 2)
            * (constants.ELECTRON_CHARGE / 1e3)
            * nd_plasma_electrons_vol_avg
            * temp_plasma_electron_density_weighted_kev
            / t_electron_energy_confinement
        )

        ratio = (nd_plasma_ions_total_vol_avg / nd_plasma_electrons_vol_avg) * (
            temp_plasma_ion_density_weighted_kev
            / temp_plasma_electron_density_weighted_kev
        )

        # Global energy confinement time

        t_energy_confinement = (ratio + 1.0e0) / (
            ratio / t_ion_energy_confinement + 1.0e0 / t_electron_energy_confinement
        )

        # For comparison directly calculate the confinement time from the stored energy calculated
        # from the total plasma beta and the loss power used above.
        physics_variables.t_energy_confinement_beta = (
            physics_variables.e_plasma_beta / 1e6
        ) / p_plasma_loss_mw

        return (
            pden_electron_transport_loss_mw,
            pden_ion_transport_loss_mw,
            t_electron_energy_confinement,
            t_ion_energy_confinement,
            t_energy_confinement,
            p_plasma_loss_mw,
        )

    @staticmethod
    def calculate_plasma_masses(
        m_fuel_amu: float,
        m_ions_total_amu: float,
        nd_plasma_ions_total_vol_avg: float,
        nd_plasma_fuel_ions_vol_avg: float,
        nd_plasma_alphas_vol_avg: float,
        vol_plasma: float,
        nd_plasma_electrons_vol_avg: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma masses.

        Parameters
        ----------
        m_fuel_amu : float
            Average mass of fuel (amu).
        m_ions_total_amu : float
            Average mass of all ions (amu).
        nd_plasma_ions_total_vol_avg : float
            Total ion density (/m3).
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (/m3).
        nd_plasma_alphas_vol_avg : float
            Alpha ash density (/m3).
        vol_plasma : float
            Plasma volume (m3).
        nd_plasma_electrons_vol_avg : float
            Volume averaged electron density (/m3).

        Returns
        -------
        tuple[float, float, float, float, float]
            A tuple containing:

        """

        # Calculate mass of fuel ions
        m_plasma_fuel_ions = (m_fuel_amu * constants.ATOMIC_MASS_UNIT) * (
            nd_plasma_fuel_ions_vol_avg * vol_plasma
        )

        m_plasma_ions_total = (m_ions_total_amu * constants.ATOMIC_MASS_UNIT) * (
            nd_plasma_ions_total_vol_avg * vol_plasma
        )

        m_plasma_alpha = (nd_plasma_alphas_vol_avg * vol_plasma) * constants.ALPHA_MASS

        m_plasma_electron = constants.ELECTRON_MASS * (
            nd_plasma_electrons_vol_avg * vol_plasma
        )

        m_plasma = m_plasma_electron + m_plasma_ions_total

        return (
            m_plasma_fuel_ions,
            m_plasma_ions_total,
            m_plasma_alpha,
            m_plasma_electron,
            m_plasma,
        )


def res_diff_time(rmajor, res_plasma, kappa95):
    """Calculates resistive diffusion time

    Parameters
    ----------
    rmajor :
        plasma major radius (m)
    res_plasma :
        plasma resistivity (Ohms)
    kappa95 :
        plasma elongation at 95% flux surface

    """

    return 2 * constants.RMU0 * rmajor / (res_plasma * kappa95)


def l_h_threshold_power(
    nd_plasma_electron_line: float,
    b_plasma_toroidal_on_axis: float,
    rmajor: float,
    rminor: float,
    kappa: float,
    a_plasma_surface: float,
    m_ions_total_amu: float,
    aspect: float,
    plasma_current: float,
) -> list[float]:
    """L-mode to H-mode power threshold calculation.

    Parameters
    ----------
    nd_plasma_electron_line : float
        Line-averaged electron density (/m3)
    b_plasma_toroidal_on_axis : float
        Toroidal field on axis (T)
    rmajor : float
        Plasma major radius (m)
    rminor : float
        Plasma minor radius (m)
    kappa : float
        Plasma elongation
    a_plasma_surface : float
        Plasma surface area (m2)
    m_ions_total_amu : float
        Average mass of all ions (amu)
    aspect : float
        Aspect ratio
    plasma_current : float
        Plasma current (A)

    Returns
    -------
    list[float]
        Array of power thresholds

    """

    dnla20 = 1e-20 * nd_plasma_electron_line

    # ========================================================================

    # ITER-1996 H-mode power threshold database
    # Fit to 1996 H-mode power threshold database: nominal

    # i_l_h_threshold = 1
    iterdd = transition.calculate_iter1996_nominal(
        dnla20, b_plasma_toroidal_on_axis, rmajor
    )

    # Fit to 1996 H-mode power threshold database: upper bound
    # i_l_h_threshold = 2
    iterdd_ub = transition.calculate_iter1996_upper(
        dnla20, b_plasma_toroidal_on_axis, rmajor
    )

    # Fit to 1996 H-mode power threshold database: lower bound
    # i_l_h_threshold = 3
    iterdd_lb = transition.calculate_iter1996_lower(
        dnla20, b_plasma_toroidal_on_axis, rmajor
    )

    # ========================================================================

    # Snipes 1997 ITER H-mode power threshold

    # i_l_h_threshold = 4
    snipes_1997 = transition.calculate_snipes1997_iter(
        dnla20, b_plasma_toroidal_on_axis, rmajor
    )

    # i_l_h_threshold = 5
    snipes_1997_kappa = transition.calculate_snipes1997_kappa(
        dnla20, b_plasma_toroidal_on_axis, rmajor, kappa
    )

    # ========================================================================

    # Martin et al (2008) for recent ITER scaling, with mass correction
    # and 95% confidence limits

    # i_l_h_threshold = 6
    martin_nominal = transition.calculate_martin08_nominal(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
    )

    # i_l_h_threshold = 7
    martin_ub = transition.calculate_martin08_upper(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
    )

    # i_l_h_threshold = 8
    martin_lb = transition.calculate_martin08_lower(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
    )

    # ========================================================================

    # Snipes et al (2000) scaling with mass correction
    # Nominal, upper and lower

    # i_l_h_threshold = 9
    snipes_2000 = transition.calculate_snipes2000_nominal(
        dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
    )

    # i_l_h_threshold = 10
    snipes_2000_ub = transition.calculate_snipes2000_upper(
        dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
    )

    # i_l_h_threshold = 11
    snipes_2000_lb = transition.calculate_snipes2000_lower(
        dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
    )

    # ========================================================================

    # Snipes et al (2000) scaling (closed divertor) with mass correction
    # Nominal, upper and lower

    # i_l_h_threshold = 12
    snipes_2000_cd = transition.calculate_snipes2000_closed_divertor_nominal(
        dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
    )

    # i_l_h_threshold = 13
    snipes_2000_cd_ub = transition.calculate_snipes2000_closed_divertor_upper(
        dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
    )

    # i_l_h_threshold = 14
    snipes_2000_cd_lb = transition.calculate_snipes2000_closed_divertor_lower(
        dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
    )

    # ========================================================================

    # Hubbard et al. 2012 L-I threshold scaling

    # i_l_h_threshold = 15
    hubbard_2012 = transition.calculate_hubbard2012_nominal(plasma_current, dnla20)

    # i_l_h_threshold = 16
    hubbard_2012_lb = transition.calculate_hubbard2012_lower(plasma_current, dnla20)

    # i_l_h_threshold = 17
    hubbard_2012_ub = transition.calculate_hubbard2012_upper(plasma_current, dnla20)

    # ========================================================================

    # Hubbard et al. 2017 L-I threshold scaling

    # i_l_h_threshold = 18
    hubbard_2017 = transition.calculate_hubbard2017(
        dnla20, a_plasma_surface, b_plasma_toroidal_on_axis
    )

    # ========================================================================

    # Aspect ratio corrected Martin et al (2008)

    # i_l_h_threshold = 19
    martin_nominal_aspect = transition.calculate_martin08_aspect_nominal(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
    )

    # i_l_h_threshold = 20
    martin_ub_aspect = transition.calculate_martin08_aspect_upper(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
    )

    # i_l_h_threshold = 21
    martin_lb_aspect = transition.calculate_martin08_aspect_lower(
        dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
    )

    # ========================================================================

    return [
        iterdd,
        iterdd_ub,
        iterdd_lb,
        snipes_1997,
        snipes_1997_kappa,
        martin_nominal,
        martin_ub,
        martin_lb,
        snipes_2000,
        snipes_2000_ub,
        snipes_2000_lb,
        snipes_2000_cd,
        snipes_2000_cd_ub,
        snipes_2000_cd_lb,
        hubbard_2012,
        hubbard_2012_lb,
        hubbard_2012_ub,
        hubbard_2017,
        martin_nominal_aspect,
        martin_ub_aspect,
        martin_lb_aspect,
    ]


def reinke_tsep(b_plasma_toroidal_on_axis, flh, qstar, rmajor, eps, fgw, kappa, lhat):
    """Function for calculating upstream temperature(keV) in Reinke model
    This function calculates the upstream temperature in the
    divertor/SoL model used for the Reinke citerion.
    Issue #707
    M.L. Reinke 2017 Nucl. Fusion 57 034004

    Parameters
    ---------_
    b_plasma_toroidal_on_axis :
        toroidal field on axis (T)
    flh :
        fraction of Psep/P_LH
    qstar :
        safety factor similar to q95 (see #707)
    rmajor :
        major radius (m)
    eps :
        inverse aspect ratio
    fgw :
        ratio of volume averaged density to n_GW
    kappa :
        elongation
    lhat :
        connection length factor
    """
    kappa_0 = 2.0e3  # Stangeby W/m/eV^(7/2)

    return (
        (
            b_plasma_toroidal_on_axis**0.72
            * flh**0.29
            * fgw**0.21
            * qstar**0.08
            * rmajor**0.33
        )
        * (eps**0.15 * (1.0 + kappa**2.0) ** 0.34)
        * (lhat**0.29 * kappa_0 ** (-0.29) * 0.285)
    )


class BetaNormMaxModel(IntEnum):
    """Beta norm max (β_N_max) model types"""

    USER_INPUT = 0
    WESSON = 1
    ORIGINAL_SCALING = 2
    MENARD = 3
    THLOREUS = 4
    STAMBAUGH = 5


class PlasmaBeta:
    """Class to hold plasma beta calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def get_beta_norm_max_value(self, model: BetaNormMaxModel) -> float:
        """Get the beta norm max value (β_N_max) for the specified model.

        Parameters
        ----------
        model: BetaNormMaxModel :
        """
        model_map = {
            BetaNormMaxModel.USER_INPUT: physics_variables.beta_norm_max,
            BetaNormMaxModel.WESSON: physics_variables.beta_norm_max_wesson,
            BetaNormMaxModel.ORIGINAL_SCALING: physics_variables.beta_norm_max_original_scaling,
            BetaNormMaxModel.MENARD: physics_variables.beta_norm_max_menard,
            BetaNormMaxModel.THLOREUS: physics_variables.beta_norm_max_thloreus,
            BetaNormMaxModel.STAMBAUGH: physics_variables.beta_norm_max_stambaugh,
        }
        return model_map[model]

    def run(self):
        """Calculate plasma beta values."""

        # -----------------------------------------------------
        # Normalised Beta Limit
        # -----------------------------------------------------

        # Normalised beta from Troyon beta limit
        physics_variables.beta_norm_total = self.calculate_normalised_beta(
            beta=physics_variables.beta_total_vol_avg,
            rminor=physics_variables.rminor,
            c_plasma=physics_variables.plasma_current,
            b_field=physics_variables.b_plasma_toroidal_on_axis,
        )

        # Define beta_norm_max calculations

        physics_variables.beta_norm_max_wesson = self.calculate_beta_norm_max_wesson(
            ind_plasma_internal_norm=physics_variables.ind_plasma_internal_norm
        )

        # Original scaling law
        physics_variables.beta_norm_max_original_scaling = (
            self.calculate_beta_norm_max_original(eps=physics_variables.eps)
        )

        # J. Menard scaling law
        physics_variables.beta_norm_max_menard = self.calculate_beta_norm_max_menard(
            eps=physics_variables.eps
        )

        # E. Tholerus scaling law
        physics_variables.beta_norm_max_thloreus = self.calculate_beta_norm_max_thloreus(
            c_beta=physics_variables.c_beta,
            pres_plasma_on_axis=physics_variables.pres_plasma_thermal_on_axis,
            pres_plasma_vol_avg=physics_variables.pres_plasma_thermal_vol_avg,
        )

        # R. D. Stambaugh scaling law
        physics_variables.beta_norm_max_stambaugh = (
            self.calculate_beta_norm_max_stambaugh(
                f_c_plasma_bootstrap=current_drive_variables.f_c_plasma_bootstrap,
                kappa=physics_variables.kappa,
                aspect=physics_variables.aspect,
            )
        )

        # Calculate beta_norm_max based on i_beta_norm_max
        try:
            model = BetaNormMaxModel(int(physics_variables.i_beta_norm_max))
            physics_variables.beta_norm_max = self.get_beta_norm_max_value(model)
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_beta_norm_max",
                i_beta_norm_max=physics_variables.i_beta_norm_max,
            ) from None

        # calculate_beta_limit() returns the beta_vol_avg_max for beta
        physics_variables.beta_vol_avg_max = self.calculate_beta_limit_from_norm(
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            beta_norm_max=physics_variables.beta_norm_max,
            plasma_current=physics_variables.plasma_current,
            rminor=physics_variables.rminor,
        )

        physics_variables.beta_toroidal_vol_avg = (
            physics_variables.beta_total_vol_avg
            * physics_variables.b_plasma_total**2
            / physics_variables.b_plasma_toroidal_on_axis**2
        )

        # Calculate physics_variables.beta poloidal [-]
        physics_variables.beta_poloidal_vol_avg = self.calculate_poloidal_beta(
            b_plasma_total=physics_variables.b_plasma_total,
            b_plasma_poloidal_average=physics_variables.b_plasma_poloidal_average,
            beta=physics_variables.beta_total_vol_avg,
        )

        physics_variables.beta_thermal_vol_avg = (
            physics_variables.beta_total_vol_avg
            - physics_variables.beta_fast_alpha
            - physics_variables.beta_beam
        )

        physics_variables.beta_poloidal_eps = (
            physics_variables.beta_poloidal_vol_avg * physics_variables.eps
        )

        physics_variables.beta_thermal_poloidal_vol_avg = (
            physics_variables.beta_thermal_vol_avg
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_poloidal_average
            )
            ** 2
        )
        physics_variables.beta_thermal_toroidal_vol_avg = (
            physics_variables.beta_thermal_vol_avg
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )

        # =======================================================

        # Mirror the pressure profiles to match the doubled toroidal field profile
        pres_profile_total = np.concatenate([
            physics_variables.pres_plasma_thermal_total_profile[::-1],
            physics_variables.pres_plasma_thermal_total_profile,
        ])

        physics_variables.beta_thermal_toroidal_profile = np.array([
            self.calculate_plasma_beta(
                pres_plasma=pres_profile_total[i],
                b_field=physics_variables.b_plasma_toroidal_profile[i],
            )
            for i in range(len(physics_variables.b_plasma_toroidal_profile))
        ])

        # =======================================================

        physics_variables.beta_norm_thermal = self.calculate_normalised_beta(
            beta=physics_variables.beta_thermal_vol_avg,
            rminor=physics_variables.rminor,
            c_plasma=physics_variables.plasma_current,
            b_field=physics_variables.b_plasma_toroidal_on_axis,
        )

        physics_variables.beta_norm_toroidal = (
            physics_variables.beta_norm_total
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )
        physics_variables.beta_norm_poloidal = (
            physics_variables.beta_norm_total
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_poloidal_average
            )
            ** 2
        )

        physics_variables.f_beta_alpha_beam_thermal = (
            physics_variables.beta_fast_alpha + physics_variables.beta_beam
        ) / physics_variables.beta_thermal_vol_avg

        # Plasma thermal energy derived from the thermal beta
        physics_variables.e_plasma_beta_thermal = self.calculate_plasma_energy_from_beta(
            beta=physics_variables.beta_thermal_vol_avg,
            b_field=physics_variables.b_plasma_total,
            vol_plasma=physics_variables.vol_plasma,
        )

        # Plasma thermal energy derived from the total beta
        physics_variables.e_plasma_beta = self.calculate_plasma_energy_from_beta(
            beta=physics_variables.beta_total_vol_avg,
            b_field=physics_variables.b_plasma_total,
            vol_plasma=physics_variables.vol_plasma,
        )

    @staticmethod
    def calculate_plasma_beta(
        pres_plasma: float | np.ndarray, b_field: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the plasma beta (β) for a given pressure and field.

        Plasma beta is the ratio of plasma pressure to magnetic pressure.

        Parameters
        ----------
        pres_plasma : float | np.ndarray
            Plasma pressure (in Pascals).
        b_field : float | np.ndarray
            Magnetic field strength (in Tesla).

        Returns
        -------
        float | np.ndarray
            The plasma beta (dimensionless).
        """

        return 2 * constants.RMU0 * pres_plasma / (b_field**2)

    @staticmethod
    def calculate_beta_norm_max_wesson(ind_plasma_internal_norm: float) -> float:
        """Calculate the Wesson normalsied beta upper limit.

        Parameters
        ----------
        ind_plasma_internal_norm : float
            Plasma normalised internal inductance


         Returns
        -------
        float
            The Wesson normalised beta upper limit.

         Notes
         -----
             - It is recommended to use this method with the other Wesson relations for normalsied internal
             inductance and current profile index.
             - This fit is derived from the DIII-D database for β_N >= 2.5

         References
         ----------
             - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
             International Series of Monographs on Physics, Volume 149.

             - T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,”
             Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).

        """
        return 4 * ind_plasma_internal_norm

    @staticmethod
    def calculate_beta_norm_max_original(eps: float) -> float:
        """Calculate the original scaling law normalsied beta upper limit.

        Parameters
        ----------
        eps : float
            Plasma normalised internal inductance

        Returns
        -------
        float

        References
        ----------
            The original scaling law normalised beta upper limit.

        """
        return 2.7 * (1.0 + 5.0 * eps**3.5)

    @staticmethod
    def calculate_beta_norm_max_menard(eps: float) -> float:
        """Calculate the Menard normalsied beta upper limit.

        Parameters
        ----------
        eps : float
            Plasma normalised internal inductance

        Returns
        -------
        float
            The Menard normalised beta upper limit.

        Notes
        -----
            - Found as a reasonable fit to the computed no wall limit at f_BS ≈ 50%
            - Uses maximum κ data from NSTX at A = 1.45, A = 1.75. Along with record
              β_T data from DIII-D at A = 2.9 and high κ.

        References
        ----------
            - # J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.

        """
        return 3.12 + 3.5 * eps**1.7

    @staticmethod
    def calculate_beta_norm_max_thloreus(
        c_beta: float, pres_plasma_on_axis: float, pres_plasma_vol_avg: float
    ) -> float:
        """Calculate the E. Tholerus normalized beta upper limit.

        Parameters
        ----------
        c_beta : float
            Pressure peaking factor coefficient.
        pres_plasma_on_axis : float
            Central plasma pressure (Pa).
        pres_plasma_vol_avg : float
            Volume-averaged plasma pressure (Pa).
        c_beta: float :

        pres_plasma_on_axis: float :

        pres_plasma_vol_avg: float :


        Returns
        -------
        float
            The E. Tholerus normalized beta upper limit.

        Notes
        -----
            - This method calculates the normalized beta upper limit based on the pressure peaking factor (Fp),
              which is defined as the ratio of the peak pressure to the average pressure.
            - The formula is derived from operational space studies of flat-top plasma in the STEP power plant.

        References
        ----------
            - E. Tholerus et al., “Flat-top plasma operational space of the STEP power plant,”
              Nuclear Fusion, Aug. 2024, doi: https://doi.org/10.1088/1741-4326/ad6ea2.

        """
        return 3.7 + (
            (c_beta / (pres_plasma_on_axis / pres_plasma_vol_avg))
            * (12.5 - 3.5 * (pres_plasma_on_axis / pres_plasma_vol_avg))
        )

    @staticmethod
    def calculate_beta_norm_max_stambaugh(
        f_c_plasma_bootstrap: float,
        kappa: float,
        aspect: float,
    ) -> float:
        """Calculate the Stambaugh normalized beta upper limit.

        Parameters
        ----------
        f_c_plasma_bootstrap : float
            Bootstrap current fraction.
        kappa : float
            Plasma separatrix elongation.
        aspect : float
            Plasma aspect ratio.


        Returns
        -------
        float
            The Stambaugh normalized beta upper limit.

        Notes
        -----
            - This method calculates the normalized beta upper limit based on the Stambaugh scaling.
            - The formula is derived from empirical fits to high-performance, steady-state tokamak equilibria.

        References
        ----------
            - R. D. Stambaugh et al., “Fusion Nuclear Science Facility Candidates,”
              Fusion Science and Technology, vol. 59, no. 2, pp. 279-307, Feb. 2011,
              doi: https://doi.org/10.13182/fst59-279.

            - Y. R. Lin-Liu and R. D. Stambaugh, “Optimum equilibria for high performance, steady state tokamaks,”
              Nuclear Fusion, vol. 44, no. 4, pp. 548-554, Mar. 2004,
              doi: https://doi.org/10.1088/0029-5515/44/4/009.

        """
        return (
            f_c_plasma_bootstrap
            * 10
            * (-0.7748 + (1.2869 * kappa) - (0.2921 * kappa**2) + (0.0197 * kappa**3))
            / (aspect**0.5523 * np.tanh((1.8524 + (0.2319 * kappa)) / aspect**0.6163))
        )

    @staticmethod
    def calculate_normalised_beta(
        beta: float, rminor: float, c_plasma: float, b_field: float
    ) -> float:
        """Calculate normalised beta (β_N).

        Parameters
        ----------
        beta : float
            Plasma beta (fraction).
        rminor : float
            Plasma minor radius (m).
        c_plasma : float
            Plasma current (A).
        b_field : float
            Magnetic field (T).

        Returns
        -------
        float
            Normalised beta.

        Notes
        -----
        - 1.0e8 is a conversion factor to get beta_N in standard units, as plasma current is normally in MA and
        beta is in percentage instead of fraction.

        """

        return 1.0e8 * (beta * rminor * b_field) / c_plasma

    @staticmethod
    def calculate_plasma_energy_from_beta(
        beta: float, b_field: float, vol_plasma: float
    ) -> float:
        """Calculate plasma thermal energy from beta.

        E_plasma = 1.5 * β * B² / (2 * μ_0) * V

        Parameters
        ----------
        beta : float
            Plasma beta (fraction).
        b_field : float
            Magnetic field (T).
        vol_plasma : float
            Plasma volume (m³).

        Returns
        -------
        float
            Plasma energy (J).
        """

        return (1.5e0 * beta * b_field**2) / (2.0e0 * constants.RMU0) * vol_plasma

    @staticmethod
    def calculate_beta_limit_from_norm(
        b_plasma_toroidal_on_axis: float,
        beta_norm_max: float,
        plasma_current: float,
        rminor: float,
    ) -> float:
        """Calculate the maximum allowed beta (β) from a given normalised (β_N).

        This subroutine calculates the beta limit using the algorithm documented in AEA FUS 172.
        The limit applies to beta defined with respect to the total B-field.
        Switch i_beta_component determines which components of beta to include.

        Parameters
        ----------
        b_plasma_toroidal_on_axis : float
            Toroidal B-field on plasma axis [T].
        beta_norm_max : float
            Troyon-like g coefficient.
        plasma_current : float
            Plasma current [A].
        rminor : float
            Plasma minor axis [m].

        Returns
        -------
        float
            Beta limit as defined below.

        Notes
        -----
            - If i_beta_component = 0, then the limit is applied to the total beta.
            - If i_beta_component = 1, then the limit is applied to the thermal beta only.
            - If i_beta_component = 2, then the limit is applied to the thermal + neutral beam beta components.
            - If i_beta_component = 3, then the limit is applied to the toroidal beta.

            - The default value for the g coefficient is beta_norm_max = 3.5.

        References
        ----------
            - F. Troyon et.al,  “Beta limit in tokamaks. Experimental and computational status,”
            Plasma Physics and Controlled Fusion, vol. 30, no. 11, pp. 1597-1609, Oct. 1988,
            doi: https://doi.org/10.1088/0741-3335/30/11/019.

            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """

        # Multiplied by 0.01 to convert from % to fraction
        return (
            0.01
            * beta_norm_max
            * (plasma_current / 1.0e6)
            / (rminor * b_plasma_toroidal_on_axis)
        )

    @staticmethod
    def calculate_poloidal_beta(
        b_plasma_total: float, b_plasma_poloidal_average: float, beta: float
    ) -> float:
        """Calculates total poloidal beta (β_p)

        Parameters
        ----------
        b_plasma_poloidal_average : float
            The average poloidal magnetic field of the plasma (in Tesla).
        beta : float
            The plasma beta, a dimensionless parameter representing the ratio of plasma pressure to magnetic pressure.

        Returns
        -------
        float
            The calculated total poloidal beta.

        References
        ----------
        - J.P. Freidberg, "Plasma physics and fusion energy", Cambridge University Press (2007)
        Page 270 ISBN 0521851076

        """
        return beta * (b_plasma_total / b_plasma_poloidal_average) ** 2

    @staticmethod
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
        """Calculate the fast alpha beta (β_fast_alpha) component.

        This function computes the fast alpha beta contribution based on the provided plasma parameters.

        Parameters
        ----------
        b_plasma_poloidal_average : float
            Poloidal field (T).
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        nd_plasma_electrons_vol_avg : float
            Electron density (m⁻³).
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (m⁻³).
        nd_plasma_ions_total_vol_avg : float
            Total ion density (m⁻³).
        temp_plasma_electron_density_weighted_kev : float
            Density-weighted electron temperature (keV).
        temp_plasma_ion_density_weighted_kev : float
            Density-weighted ion temperature (keV).
        pden_alpha_total_mw : float
            Alpha power per unit volume, from beams and plasma (MW/m³).
        pden_plasma_alpha_mw : float
            Alpha power per unit volume just from plasma (MW/m³).
        i_beta_fast_alpha : int
            Switch for fast alpha pressure method.

        Returns
        -------
        float
            Fast alpha beta component.

        Notes
        -----
            - For IPDG89 scaling applicability is Z_eff = 1.5, T_i/T_e = 1, 〈T〉 = 5-20 keV


        References
        ----------
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
                    nd_plasma_electrons_vol_avg
                    * temp_plasma_electron_density_weighted_kev
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

    def output_beta_information(self):
        """Output beta information to file."""

        po.osubhd(self.outfile, "Beta Information :")
        if physics_variables.i_beta_component == 0:
            po.ovarrf(
                self.outfile,
                "Upper limit on total beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )
        elif physics_variables.i_beta_component == 1:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal + NB beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total plasma beta",
            "(beta_total_vol_avg)",
            physics_variables.beta_total_vol_avg,
        )
        if physics_variables.i_beta_component == 0:
            po.ovarrf(
                self.outfile,
                "Lower limit on total beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        elif physics_variables.i_beta_component == 1:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal + NB beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        po.ovarre(
            self.outfile,
            "Upper limit on poloidal beta",
            "(beta_poloidal_max)",
            constraint_variables.beta_poloidal_max,
            "IP",
        )
        po.ovarre(
            self.outfile,
            "Total poloidal beta",
            "(beta_poloidal_vol_avg)",
            physics_variables.beta_poloidal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged toroidal beta",
            "(beta_toroidal_vol_avg)",
            physics_variables.beta_toroidal_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.beta_thermal_toroidal_profile)):
            po.ovarre(
                self.mfile,
                f"Beta toroidal profile at point {i}",
                f"beta_thermal_toroidal_profile{i}",
                physics_variables.beta_thermal_toroidal_profile[i],
            )

        po.ovarre(
            self.outfile,
            "Fast alpha beta",
            "(beta_fast_alpha)",
            physics_variables.beta_fast_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutral Beam ion beta",
            "(beta_beam)",
            physics_variables.beta_beam,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ratio of fast alpha and beam beta to thermal beta",
            "(f_beta_alpha_beam_thermal)",
            physics_variables.f_beta_alpha_beam_thermal,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Volume averaged thermal beta",
            "(beta_thermal_vol_avg)",
            physics_variables.beta_thermal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal poloidal beta",
            "(beta_thermal_poloidal_vol_avg)",
            physics_variables.beta_thermal_poloidal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal toroidal beta",
            "(beta_thermal_toroidal_vol_avg)",
            physics_variables.beta_thermal_toroidal_vol_avg,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Poloidal beta and inverse aspect ratio",
            "(beta_poloidal_eps)",
            physics_variables.beta_poloidal_eps,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Poloidal beta and inverse aspect ratio upper limit",
            "(beta_poloidal_eps_max)",
            physics_variables.beta_poloidal_eps_max,
        )
        po.osubhd(self.outfile, "Normalised Beta Information :")
        if stellarator_variables.istell == 0:
            if physics_variables.i_beta_norm_max != 0:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(beta_norm_max)",
                    physics_variables.beta_norm_max,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(beta_norm_max)",
                    physics_variables.beta_norm_max,
                )
            po.ovarrf(
                self.outfile,
                "Normalised total beta",
                "(beta_norm_total)",
                physics_variables.beta_norm_total,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Normalised thermal beta",
                "(beta_norm_thermal) ",
                physics_variables.beta_norm_thermal,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised toroidal beta",
                "(beta_norm_toroidal) ",
                physics_variables.beta_norm_toroidal,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised poloidal beta",
                "(beta_norm_poloidal) ",
                physics_variables.beta_norm_poloidal,
                "OP ",
            )

            po.osubhd(self.outfile, "Maximum normalised beta scalings :")
            po.ovarrf(
                self.outfile,
                "J. Wesson normalised beta upper limit",
                "(beta_norm_max_wesson) ",
                physics_variables.beta_norm_max_wesson,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Original normalsied beta upper limit",
                "(beta_norm_max_original_scaling) ",
                physics_variables.beta_norm_max_original_scaling,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "J. Menard normalised beta upper limit",
                "(beta_norm_max_menard) ",
                physics_variables.beta_norm_max_menard,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "E. Thloreus normalised beta upper limit",
                "(beta_norm_max_thloreus) ",
                physics_variables.beta_norm_max_thloreus,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "R. Stambaugh normalised beta upper limit",
                "(beta_norm_max_stambaugh) ",
                physics_variables.beta_norm_max_stambaugh,
                "OP ",
            )

        po.osubhd(self.outfile, "Plasma energies derived from beta :")
        po.ovarre(
            self.outfile,
            "Plasma thermal energy derived from thermal beta (J)",
            "(e_plasma_beta_thermal) ",
            physics_variables.e_plasma_beta_thermal,
            "OP",
        )

        po.ovarre(
            self.outfile,
            "Plasma thermal energy derived from the total beta (J)",
            "(e_plasma_beta)",
            physics_variables.e_plasma_beta,
            "OP",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)


class IndInternalNormModel(IntEnum):
    """Normalised internal inductance (l_i) model types"""

    USER_INPUT = 0
    WESSON = 1
    MENARD = 2


class PlasmaInductance:
    """Class to hold plasma inductance calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        physics_variables.ind_plasma_internal_norm_wesson = (
            self.calculate_internal_inductance_wesson(alphaj=physics_variables.alphaj)
        )

        # Spherical Tokamak relation for internal inductance
        # Menard et al. (2016), Nuclear Fusion, 56, 106023
        physics_variables.ind_plasma_internal_norm_menard = (
            self.calculate_internal_inductance_menard(kappa=physics_variables.kappa)
        )

        physics_variables.ind_plasma_internal_norm_iter_3 = (
            self.calculate_normalised_internal_inductance_iter_3(
                b_plasma_poloidal_vol_avg=physics_variables.b_plasma_poloidal_average,
                c_plasma=physics_variables.plasma_current,
                vol_plasma=physics_variables.vol_plasma,
                rmajor=physics_variables.rmajor,
            )
        )

        # Calculate ind_plasma_internal_norm based on i_ind_plasma_internal_norm
        try:
            model = IndInternalNormModel(
                int(physics_variables.i_ind_plasma_internal_norm)
            )
            physics_variables.ind_plasma_internal_norm = (
                self.get_ind_internal_norm_value(model)
            )
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_ind_plasma_internal_norm",
                i_ind_plasma_internal_norm=physics_variables.i_ind_plasma_internal_norm,
            ) from None

    def get_ind_internal_norm_value(self, model: IndInternalNormModel) -> float:
        """Get the normalised internal inductance (l_i) for the specified model.

        Parameters
        ----------
        model: IndInternalNormModel :
        """
        model_map = {
            IndInternalNormModel.USER_INPUT: physics_variables.ind_plasma_internal_norm,
            IndInternalNormModel.WESSON: physics_variables.ind_plasma_internal_norm_wesson,
            IndInternalNormModel.MENARD: physics_variables.ind_plasma_internal_norm_menard,
        }
        return model_map[model]

    @staticmethod
    def calculate_volt_second_requirements(
        csawth: float,
        eps: float,
        f_c_plasma_inductive: float,
        ejima_coeff: float,
        kappa: float,
        rmajor: float,
        res_plasma: float,
        plasma_current: float,
        t_plant_pulse_fusion_ramp: float,
        t_plant_pulse_burn: float,
        ind_plasma_internal_norm: float,
    ) -> tuple[float, float, float, float, float, float]:
        """Calculate the volt-second requirements and related parameters for plasma physics.

        Parameters
        ----------
        csawth : float
            Coefficient for sawteeth effects
        eps : float
            Inverse aspect ratio
        f_c_plasma_inductive : float
            Fraction of plasma current produced inductively
        ejima_coeff : float
            Ejima coefficient for resistive start-up V-s component
        kappa : float
            Plasma elongation
        rmajor : float
            Plasma major radius (m)
        res_plasma : float
            Plasma resistance (ohm)
        plasma_current : float
            Plasma current (A)
        t_plant_pulse_fusion_ramp : float
            Heating time (s)
        t_plant_pulse_burn : float
            Burn time (s)
        ind_plasma_internal_norm : float
            Plasma normalized internal inductance


        Returns
        -------
        tuple[float, float, float, float, float, float]
            A tuple containing:
            - vs_plasma_internal: Internal plasma volt-seconds (Wb)
            - ind_plasma_internal: Plasma inductance (H)
            - vs_plasma_burn_required: Volt-seconds needed during flat-top (heat+burn) (Wb)
            - vs_plasma_ramp_required: Volt-seconds needed during ramp-up (Wb)
            - ind_plasma_total,: Internal and external plasma inductance V-s (Wb)
            - vs_res_ramp: Resistive losses in start-up volt-seconds (Wb)
            - vs_plasma_total_required: Total volt-seconds needed (Wb)

        References
        ----------
            - S. Ejima, R. W. Callis, J. L. Luxon, R. D. Stambaugh, T. S. Taylor, and J. C. Wesley,
            “Volt-second analysis and consumption in Doublet III plasmas,”
            Nuclear Fusion, vol. 22, no. 10, pp. 1313-1319, Oct. 1982, doi:
            https://doi.org/10.1088/0029-5515/22/10/006.

            - S. C. Jardin, C. E. Kessel, and N Pomphrey,
            “Poloidal flux linkage requirements for the International Thermonuclear Experimental Reactor,”
            Nuclear Fusion, vol. 34, no. 8, pp. 1145-1160, Aug. 1994,
            doi: https://doi.org/10.1088/0029-5515/34/8/i07.

            - S. P. Hirshman and G. H. Neilson, “External inductance of an axisymmetric plasma,”
            The Physics of Fluids, vol. 29, no. 3, pp. 790-793, Mar. 1986,
            doi: https://doi.org/10.1063/1.865934.

        """
        # Plasma internal inductance

        ind_plasma_internal = constants.RMU0 * rmajor * ind_plasma_internal_norm / 2.0

        # Internal plasma flux (V-s) component
        vs_plasma_internal = ind_plasma_internal * plasma_current

        # Start-up resistive component
        # Uses ITER formula without the 10 V-s add-on

        vs_res_ramp = ejima_coeff * constants.RMU0 * plasma_current * rmajor

        # ======================================================================

        # Hirshman and Neilson fit for external inductance

        aeps = (1.0 + 1.81 * np.sqrt(eps) + 2.05 * eps) * np.log(8.0 / eps) - (
            2.0 + 9.25 * np.sqrt(eps) - 1.21 * eps
        )
        beps = 0.73 * np.sqrt(eps) * (1.0 + 2.0 * eps**4 - 6.0 * eps**5 + 3.7 * eps**6)

        ind_plasma_external = (
            rmajor * constants.RMU0 * aeps * (1.0 - eps) / (1.0 - eps + beps * kappa)
        )

        # ======================================================================

        ind_plasma_total = ind_plasma_external + ind_plasma_internal

        # Inductive V-s component

        vs_self_ind_ramp = ind_plasma_total * plasma_current
        vs_plasma_ramp_required = vs_res_ramp + vs_self_ind_ramp

        # Plasma loop voltage during flat-top
        # Include enhancement factor in flattop V-s requirement
        # to account for MHD sawtooth effects.

        v_plasma_loop_burn = plasma_current * res_plasma * f_c_plasma_inductive

        v_burn_resistive = v_plasma_loop_burn * csawth

        # N.B. t_plant_pulse_burn on first iteration will not be correct
        # if the pulsed reactor option is used, but the value
        # will be correct on subsequent calls.

        vs_plasma_burn_required = v_burn_resistive * (
            t_plant_pulse_fusion_ramp + t_plant_pulse_burn
        )
        vs_plasma_total_required = vs_plasma_ramp_required + vs_plasma_burn_required

        return (
            vs_plasma_internal,
            ind_plasma_total,
            vs_plasma_burn_required,
            vs_plasma_ramp_required,
            vs_self_ind_ramp,
            vs_res_ramp,
            vs_plasma_total_required,
            v_plasma_loop_burn,
        )

    @staticmethod
    def calculate_normalised_internal_inductance_iter_3(
        b_plasma_poloidal_vol_avg: float,
        c_plasma: float,
        vol_plasma: float,
        rmajor: float,
    ) -> float:
        """Calculate the normalised internal inductance using ITER-3 scaling li(3).

        Parameters
        ------------
        b_plasma_poloidal_vol_avg : float
            Volume-averaged poloidal magnetic field (T).
        c_plasma : float
            Plasma current (A).
        vol_plasma : float
            Plasma volume (m^3).
        rmajor : float
            Plasma major radius (m).

        Returns
        --------
        float
            The li(3) normalised internal inductance.

        References
        -----------
            - T. C. Luce, D. A. Humphreys, G. L. Jackson, and W. M. Solomon,
            “Inductive flux usage and its optimization in tokamak operation,”
            Nuclear Fusion, vol. 54, no. 9, p. 093005, Jul. 2014,
            doi: https://doi.org/10.1088/0029-5515/54/9/093005.

            - G. L. Jackson et al., “ITER startup studies in the DIII-D tokamak,”
            Nuclear Fusion, vol. 48, no. 12, p. 125002, Nov. 2008,
            doi: https://doi.org/10.1088/0029-5515/48/12/125002.

        """

        return (
            2
            * vol_plasma
            * b_plasma_poloidal_vol_avg**2
            / (constants.RMU0**2 * c_plasma**2 * rmajor)
        )

    @staticmethod
    def calculate_internal_inductance_menard(kappa: float) -> float:
        """Calculate the Menard plasma normalized internal inductance.

        Parameters
        ----------
        kappa : float
            Plasma separatrix elongation.

        Returns
        -------
        float
            The Menard plasma normalised internal inductance.

        Notes
        -----
            - This relation is based off of data from NSTX for l_i in the range of 0.4-0.85
            - This model is only recommended to be used for ST's with kappa > 2.5

        References
        ----------
            - J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.

        """
        return 3.4 - kappa

    @staticmethod
    def calculate_internal_inductance_wesson(alphaj: float) -> float:
        """Calculate the Wesson plasma normalized internal inductance.

        Parameters
        ----------
        alphaj : float
            Current profile index.

        Returns
        -------
        float
            The Wesson plasma normalised internal inductance.

        Notes
        -----
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              current profile index.
            - This relation is only true for the cyclindrical plasma approximation with parabolic profiles.

        References
        ----------
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.

        """
        return np.log(1.65 + 0.89 * alphaj)

    def output_volt_second_information(self):
        """Output volt-second information to file."""

        po.osubhd(self.outfile, "Plasma Volt-second Requirements :")
        po.ovarre(
            self.outfile,
            "Total plasma volt-seconds required for pulse (Wb)",
            "(vs_plasma_total_required)",
            physics_variables.vs_plasma_total_required,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plasma inductive flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_ind_ramp)",
            physics_variables.vs_plasma_ind_ramp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma resistive flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_res_ramp)",
            physics_variables.vs_plasma_res_ramp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_ramp_required)",
            physics_variables.vs_plasma_ramp_required,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ejima coefficient",
            "(ejima_coeff)",
            physics_variables.ejima_coeff,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Internal plasma V-s",
            "(vs_plasma_internal)",
            physics_variables.vs_plasma_internal,
        )

        po.ovarre(
            self.outfile,
            "Plasma volt-seconds needed for flat-top (heat + burn times) (Wb)",
            "(vs_plasma_burn_required)",
            physics_variables.vs_plasma_burn_required,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Plasma loop voltage during burn (V)",
            "(v_plasma_loop_burn)",
            physics_variables.v_plasma_loop_burn,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Coefficient for sawtooth effects on burn V-s requirement",
            "(csawth)",
            physics_variables.csawth,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma resistance (ohm)",
            "(res_plasma)",
            physics_variables.res_plasma,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Plasma resistive diffusion time (s)",
            "(t_plasma_res_diffusion)",
            physics_variables.t_plasma_res_diffusion,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma inductance (H)",
            "(ind_plasma)",
            physics_variables.ind_plasma,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma magnetic energy stored (J)",
            "(e_plasma_magnetic_stored)",
            physics_variables.e_plasma_magnetic_stored,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Plasma normalised internal inductance",
            "(ind_plasma_internal_norm)",
            physics_variables.ind_plasma_internal_norm,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)


class PlasmaCurrent:
    """Class to hold plasma current calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def calculate_plasma_current(
        self,
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        eps: float,
        i_plasma_current: int,
        kappa: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        len_plasma_poloidal: float,
        q95: float,
        rmajor: float,
        rminor: float,
        triang: float,
        triang95: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma current.

        Args:
            alphaj (float): Current profile index.
            alphap (float): Pressure profile index.
            b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
            eps (float): Inverse aspect ratio.
            i_plasma_current (int): Current scaling model to use.
                1 = Peng analytic fit
                2 = Peng divertor scaling (TART,STAR)
                3 = Simple ITER scaling
                4 = IPDG89 scaling
                5 = Todd empirical scaling I
                6 = Todd empirical scaling II
                7 = Connor-Hastie model
                8 = Sauter scaling (allowing negative triangularity)
                9 = FIESTA ST scaling
            kappa (float): Plasma elongation.
            kappa95 (float): Plasma elongation at 95% surface.
            pres_plasma_on_axis (float): Central plasma pressure (Pa).
            len_plasma_poloidal (float): Plasma perimeter length (m).
            q95 (float): Plasma safety factor at 95% flux (= q-bar for i_plasma_current=2).
            ind_plasma_internal_norm (float): Plasma normalised internal inductance.
            rmajor (float): Major radius (m).
            rminor (float): Minor radius (m).
            triang (float): Plasma triangularity.
            triang95 (float): Plasma triangularity at 95% surface.

        Returns:
            Tuple[float, float, float,]: Tuple containing b_plasma_poloidal_average, qstar, plasma_current,

        Raises:
            ValueError: If invalid value for i_plasma_current is provided.

        Notes:
            This routine calculates the plasma current based on the edge safety factor q95.
            It will also make the current profile parameters consistent with the q-profile if required.

        References:
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code, unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
              'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
              1729-1738. https://doi.org/10.13182/FST92-A29971
            - ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al, ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
            - M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants - Part 1: Physics
            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
            - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
            - Sauter, Geometric formulas for systems codes..., FED 2016
        """
        # Aspect ratio
        aspect_ratio = 1.0 / eps

        # Only the Sauter scaling (i_plasma_current=8) is suitable for negative triangularity:
        if i_plasma_current != 8 and triang < 0.0:
            raise ProcessValueError(
                f"Triangularity is negative without i_plasma_current = 8 selected: {triang=}, {i_plasma_current=}"
            )

        # Calculate the function Fq that scales the edge q from the
        # circular cross-section cylindrical case

        # Peng analytical fit
        if i_plasma_current == 1:
            fq = self.current.calculate_current_coefficient_peng(
                eps, len_plasma_poloidal, rminor
            )

        # Peng scaling for double null divertor; TARTs [STAR Code]
        elif i_plasma_current == 2:
            plasma_current = 1.0e6 * self.current.calculate_plasma_current_peng(
                q95, aspect_ratio, eps, rminor, b_plasma_toroidal_on_axis, kappa, triang
            )

        # Simple ITER scaling (simply the cylindrical case)
        elif i_plasma_current == 3:
            fq = 1.0

        # ITER formula (IPDG89)
        elif i_plasma_current == 4:
            fq = self.calculate_current_coefficient_ipdg89(eps, kappa95, triang95)

        # Todd empirical scalings
        # D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        elif i_plasma_current in [5, 6]:
            fq = self.calculate_current_coefficient_todd(eps, kappa95, triang95, model=1)

            if i_plasma_current == 6:
                fq = self.calculate_current_coefficient_todd(
                    eps, kappa95, triang95, model=2
                )

        # Connor-Hastie asymptotically-correct expression
        elif i_plasma_current == 7:
            fq = self.calculate_current_coefficient_hastie(
                alphaj,
                alphap,
                b_plasma_toroidal_on_axis,
                triang95,
                eps,
                kappa95,
                pres_plasma_on_axis,
                constants.RMU0,
            )

        # Sauter scaling allowing negative triangularity [FED May 2016]
        # https://doi.org/10.1016/j.fusengdes.2016.04.033.
        elif i_plasma_current == 8:
            # Assumes zero squareness, note takes kappa, delta at separatrix not _95
            fq = self.calculate_current_coefficient_sauter(eps, kappa, triang)

        # FIESTA ST scaling
        # https://doi.org/10.1016/j.fusengdes.2020.111530.
        elif i_plasma_current == 9:
            fq = self.calculate_current_coefficient_fiesta(eps, kappa, triang)
        else:
            raise ProcessValueError(f"Invalid value {i_plasma_current=}")

        # Main plasma current calculation using the fq value from the different settings
        if i_plasma_current != 2:
            plasma_current = (
                (2.0 * np.pi / constants.RMU0)
                * rminor**2
                / (rmajor * q95)
                * fq
                * b_plasma_toroidal_on_axis
            )
        # i_plasma_current == 2 case covered above

        # Calculate cyclindrical safety factor from IPDG89
        qstar = (
            5.0e6
            * rminor**2
            / (rmajor * plasma_current / b_plasma_toroidal_on_axis)
            * 0.5
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

        # Calculate the poloidal field generated by the plasma current
        b_plasma_poloidal_average = calculate_poloidal_field(
            i_plasma_current,
            plasma_current,
            q95,
            aspect_ratio,
            eps,
            b_plasma_toroidal_on_axis,
            kappa,
            triang,
            len_plasma_poloidal,
            constants.RMU0,
        )

        return b_plasma_poloidal_average, qstar, plasma_current

    @staticmethod
    def calculate_current_coefficient_peng(
        eps: float, len_plasma_poloidal: float, rminor: float
    ) -> float:
        """
        Calculate the plasma current scaling coefficient for the Peng scaling from the STAR code.

        :param eps: Plasma inverse aspect ratio.
        :type eps: float
        :param len_plasma_poloidal: Plasma poloidal perimeter length [m].
        :type len_plasma_poloidal: float
        :param rminor: Plasma minor radius [m].
        :type rminor: float

        :return: The plasma current scaling coefficient.
        :rtype: float

        :references: None
        """

        return (
            (1.22 - 0.68 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (len_plasma_poloidal / (2.0 * np.pi * rminor)) ** 2
        )

    @staticmethod
    def calculate_plasma_current_peng(
        q95: float,
        aspect: float,
        eps: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
    ) -> float:
        """
        Function to calculate plasma current (Peng scaling from the STAR code)

        Parameters:
        - q95: float, 95% flux surface safety factor
        - aspect: float, plasma aspect ratio
        - eps: float, inverse aspect ratio
        - rminor: float, plasma minor radius (m)
        - b_plasma_toroidal_on_axis: float, toroidal field on axis (T)
        - kappa: float, plasma elongation
        - delta: float, plasma triangularity

        Returns:
        - float, plasma current in MA

        This function calculates the plasma current in MA,
        using a scaling from Peng, Galambos and Shipe (1992).
        It is primarily used for Tight Aspect Ratio Tokamaks and is
        selected via i_plasma_current=2.

        References:
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971
        """

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        ff1, ff2, d1, d2 = _plascar_bpol(aspect, eps, kappa, delta)

        e1 = (2.0 * kappa) / (d1 * (1.0 + delta))
        e2 = (2.0 * kappa) / (d2 * (1.0 - delta))

        return (
            rminor
            * b_plasma_toroidal_on_axis
            / qbar
            * 5.0
            * kappa
            / (2.0 * np.pi**2)
            * (np.arcsin(e1) / e1 + np.arcsin(e2) / e2)
            * (ff1 + ff2)
        )

    @staticmethod
    def calculate_current_coefficient_ipdg89(
        eps: float, kappa95: float, triang95: float
    ) -> float:
        """
        Calculate the fq coefficient from the IPDG89 guidlines used in the plasma current scaling.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient used in the IPDG89 plasma current scaling,
        based on the given plasma parameters.

        References:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.5
            * (1.17 - 0.65 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

    @staticmethod
    def calculate_current_coefficient_todd(
        eps: float, kappa95: float, triang95: float, model: int
    ) -> float:
        """
        Calculate the fq coefficient used in the two Todd plasma current scalings.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient based on the given plasma parameters for the two Todd scalings.

        References:
        - D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        # Calculate the Todd scaling based on the model
        base_scaling = (
            (1.0 + 2.0 * eps**2)
            * ((1.0 + kappa95**2) / 0.5)
            * (
                1.24
                - 0.54 * kappa95
                + 0.3 * (kappa95**2 + triang95**2)
                + 0.125 * triang95
            )
        )
        if model == 1:
            return base_scaling
        if model == 2:
            return base_scaling * (1.0 + (abs(kappa95 - 1.2)) ** 3)
        raise ProcessValueError(f"model = {model} is an invalid option")

    @staticmethod
    def calculate_current_coefficient_hastie(
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        delta95: float,
        eps: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        rmu0: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Connor-Hastie model used for scaling the plasma current.

        Parameters:
        - alphaj: float, the current profile index
        - alphap: float, the pressure profile index
        - b_plasma_toroidal_on_axis: float, the toroidal field on axis (T)
        - delta95: float, the plasma triangularity 95%
        - eps: float, the inverse aspect ratio
        - kappa95: float, the plasma elongation 95%
        - pres_plasma_on_axis: float, the central plasma pressure (Pa)
        - rmu0: float, the vacuum permeability (H/m)

        Returns:
        - float, the F coefficient

        This routine calculates the f_q coefficient used for scaling the plasma current,
        using the Connor-Hastie scaling

        Reference:
        - J.W.Connor and R.J.Hastie, Culham Lab Report CLM-M106 (1985).
        https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-M106-1.pdf
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        # Exponent in Connor-Hastie current profile
        lamda = alphaj

        # Exponent in Connor-Hastie pressure profile
        nu = alphap

        # Central plasma beta
        beta0 = 2.0 * rmu0 * pres_plasma_on_axis / (b_plasma_toroidal_on_axis**2)

        # Plasma internal inductance
        lamp1 = 1.0 + lamda
        li = lamp1 / lamda * (lamp1 / lamda * np.log(lamp1) - 1.0)

        # T/r in AEA FUS 172
        kap1 = kappa95 + 1.0
        tr = kappa95 * delta95 / kap1**2

        # E/r in AEA FUS 172
        er = (kappa95 - 1.0) / kap1

        # T primed in AEA FUS 172
        tprime = 2.0 * tr * lamp1 / (1.0 + 0.5 * lamda)

        # E primed in AEA FUS 172
        eprime = er * lamp1 / (1.0 + lamda / 3.0)

        # Delta primed in AEA FUS 172
        deltap = (0.5 * kap1 * eps * 0.5 * li) + (
            beta0 / (0.5 * kap1 * eps)
        ) * lamp1**2 / (1.0 + nu)

        # Delta/R0 in AEA FUS 172
        deltar = beta0 / 6.0 * (1.0 + 5.0 * lamda / 6.0 + 0.25 * lamda**2) + (
            0.5 * kap1 * eps
        ) ** 2 * 0.125 * (1.0 - (lamda**2) / 3.0)

        # F coefficient
        return (0.5 * kap1) ** 2 * (
            1.0
            + eps**2 * (0.5 * kap1) ** 2
            + 0.5 * deltap**2
            + 2.0 * deltar
            + 0.5 * (eprime**2 + er**2)
            + 0.5 * (tprime**2 + 4.0 * tr**2)
        )

    @staticmethod
    def calculate_current_coefficient_sauter(
        eps: float,
        kappa: float,
        triang: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Sauter model used for scaling the plasma current.

        Parameters:
        - eps: float, inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns:
        - float, the fq coefficient

        Reference:
        - O. Sauter, Geometric formulas for system codes including the effect of negative triangularity,
        Fusion Engineering and Design, Volume 112, 2016, Pages 633-645,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2016.04.033.
        """
        w07 = 1.0  # zero squareness - can be modified later if required

        return (
            (4.1e6 / 5.0e6)
            * (1.0 + 1.2 * (kappa - 1.0) + 0.56 * (kappa - 1.0) ** 2)
            * (1.0 + 0.09 * triang + 0.16 * triang**2)
            * (1.0 + 0.45 * triang * eps)
            / (1.0 - 0.74 * eps)
            * (1.0 + 0.55 * (w07 - 1.0))
        )

    @staticmethod
    def calculate_current_coefficient_fiesta(
        eps: float, kappa: float, triang: float
    ) -> float:
        """
        Calculate the fq coefficient used in the FIESTA plasma current scaling.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient based on the given plasma parameters for the FIESTA scaling.

        References:
        - S.Muldrew et.al,“PROCESS”: Systems studies of spherical tokamaks, Fusion Engineering and Design,
        Volume 154, 2020, 111530, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2020.111530.
        """
        return 0.538 * (1.0 + 2.440 * eps**2.736) * kappa**2.154 * triang**0.060


class DetailedPhysics:
    """Class to hold detailed physics models for plasma processing."""

    def __init__(self, plasma_profile):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile

    def run(self):
        # ---------------------------
        #  Debye length calculation
        # ---------------------------

        physics_variables.len_plasma_debye_electron_vol_avg = self.calculate_debye_length(
            temp_plasma_species_kev=physics_variables.temp_plasma_electron_vol_avg_kev,
            nd_plasma_species=physics_variables.nd_plasma_electrons_vol_avg,
        )

        physics_variables.len_plasma_debye_electron_profile = (
            self.calculate_debye_length(
                temp_plasma_species_kev=self.plasma_profile.teprofile.profile_y,
                nd_plasma_species=self.plasma_profile.neprofile.profile_y,
            )
        )

        # ============================
        # Particle relativistic speeds
        # ============================

        physics_variables.vel_plasma_electron_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT,
                mass=constants.ELECTRON_MASS,
            )
        )

        physics_variables.vel_plasma_deuteron_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT
                * physics_variables.f_temp_plasma_ion_electron,
                mass=constants.DEUTERON_MASS,
            )
        )

        physics_variables.vel_plasma_triton_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT
                * physics_variables.f_temp_plasma_ion_electron,
                mass=constants.TRITON_MASS,
            )
        )

        # ============================
        # Plasma frequencies
        # ============================

        physics_variables.freq_plasma_electron_profile = self.calculate_plasma_frequency(
            nd_particle=self.plasma_profile.neprofile.profile_y,
            m_particle=constants.ELECTRON_MASS,
            z_particle=1.0,
        )

        # ============================
        # Larmor frequencies
        # ============================

        physics_variables.freq_plasma_larmor_toroidal_electron_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.ELECTRON_MASS,
                z_particle=1.0,
            )
        )

        physics_variables.freq_plasma_larmor_toroidal_deuteron_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.DEUTERON_MASS,
                z_particle=1.0,
            )
        )

        physics_variables.freq_plasma_larmor_toroidal_triton_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.TRITON_MASS,
                z_particle=1.0,
            )
        )

        # ============================
        # Coulomb logarithm
        # ============================

        physics_variables.plasma_coulomb_log_electron_electron_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_electron_deuteron_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.DEUTERON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_deuteron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.DEUTERON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_deuteron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_electron_triton_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

    @staticmethod
    def calculate_debye_length(
        temp_plasma_species_kev: float | np.ndarray,
        nd_plasma_species: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the Debye length for a plasma.

        Parameters
        ----------
        temp_plasma_species_kev : float | np.ndarray
            Species temperature in keV.
        nd_plasma_species : float | np.ndarray
            Species number density (/m^3).

        Returns
        -------
        float | np.ndarray
            Debye length in meters.

        """
        return (
            (constants.EPSILON0 * temp_plasma_species_kev * constants.KILOELECTRON_VOLT)
            / (nd_plasma_species * constants.ELECTRON_CHARGE**2)
        ) ** 0.5

    @staticmethod
    def calculate_lorentz_factor(velocity: float | np.ndarray) -> float | np.ndarray:
        """Calculate the Lorentz factor for a given velocity.

        Parameters
        ----------
        velocity : float | np.ndarray
            Velocity in m/s.

        Returns
        -------
        float | np.ndarray
            Lorentz factor (dimensionless).

        """
        return 1 / (1 - (velocity / constants.SPEED_LIGHT) ** 2) ** 0.5

    @staticmethod
    def calculate_relativistic_particle_speed(
        e_kinetic: float | np.ndarray, mass: float
    ) -> float | np.ndarray:
        """Calculate the speed of a particle given its kinetic energy and mass using relativistic mechanics.

        Parameters
        ----------
        e_kinetic : float | np.ndarray
            Kinetic energy in Joules.
        mass : float
            Mass of the particle in kg.

        Returns
        -------
        float | np.ndarray
            Speed of the particle in m/s.

        """
        return (
            constants.SPEED_LIGHT
            * (1 - (1 / ((e_kinetic / (mass * constants.SPEED_LIGHT**2)) + 1) ** 2))
            ** 0.5
        )

    def calculate_coulomb_log_from_impact(
        self, impact_param_max: float, impact_param_min: float
    ) -> float:
        """Calculate the Coulomb logarithm from maximum and minimum impact parameters.

        Parameters
        ----------
        impact_param_max : float
            Maximum impact parameter in meters.
        impact_param_min : float
            Minimum impact parameter in meters.

        Returns
        -------
        float
            Coulomb logarithm (dimensionless).

        """
        return np.log(impact_param_max / impact_param_min)

    @staticmethod
    def calculate_classical_distance_of_closest_approach(
        charge1: float,
        charge2: float,
        m_reduced: float,
        vel_relative: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the classical distance of closest approach for two charged particles.

        Parameters
        ----------
        charge1 :
            Charge of particle 1 in units of elementary charge.
        charge2 :
            Charge of particle 2 in units of elementary charge.
        m_reduced:
            Reduced mass of the two-particle system in kg.
        vel_relative:
            Relative velocity of the two particles in m/s.

        Returns
        -------
        float | np.ndarray
            Distance of closest approach in meters.
        """

        return (charge1 * charge2 * constants.ELECTRON_CHARGE**2) / (
            2 * np.pi * constants.EPSILON0 * m_reduced * vel_relative**2
        )

    @staticmethod
    def calculate_debroglie_wavelength(
        mass: float, velocity: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the de Broglie wavelength of a particle.

        Parameters
        ----------
        mass : float
            Mass of the particle in kg.
        velocity : float | np.ndarray
            Velocity of the particle in m/s.

        Returns
        -------
        float | np.ndarray
            de Broglie wavelength in meters.

        :note: Reduced Planck constant (h-bar) is used in the calculation as this is for scattering.

        """
        return (constants.PLANCK_CONSTANT / (2 * np.pi)) / (mass * velocity)

    @staticmethod
    def calculate_plasma_frequency(
        nd_particle: float | np.ndarray, m_particle: float, z_particle: float
    ) -> float | np.ndarray:
        """Calculate the plasma frequency for a particle species.

        Parameters
        ----------
        nd_particle : float | np.ndarray
            Number density of the particle species (/m^3).
        m_particle : float
            Mass of the particle species (kg).
        Z_particle : float
            Charge state of the particle species (dimensionless).

        Returns
        -------
        float | np.ndarray
            Plasma frequency in Hz.

        """
        return (
            (
                (nd_particle * z_particle**2 * constants.ELECTRON_CHARGE**2)
                / (m_particle * constants.EPSILON0)
            )
            ** 0.5
        ) / (2 * np.pi)

    @staticmethod
    def calculate_larmor_frequency(
        b_field: float | np.ndarray, m_particle: float, z_particle: float
    ) -> float | np.ndarray:
        """Calculate the Larmor frequency for a particle species.

        Parameters
        ----------
        b_field : float | np.ndarray
            Magnetic field strength (T).
        m_particle : float
            Mass of the particle species (kg).
        Z_particle : float
            Charge state of the particle species (dimensionless).

        Returns
        -------
        float | np.ndarray
            Larmor frequency in Hz.

        """
        return (z_particle * constants.ELECTRON_CHARGE * b_field) / (
            2 * np.pi * m_particle
        )

    @staticmethod
    def calculate_reduced_mass(mass1: float, mass2: float) -> float:
        """
        Calculate the reduced mass of two particles.

        Parameters
        ----------
        mass1:
            Mass of particle 1 (kg).
        mass2:
            Mass of particle 2 (kg).

        Returns
        -------
            Reduced mass (kg).
        """
        return (mass1 * mass2) / (mass1 + mass2)

    @staticmethod
    def calculate_average_relative_velocity(
        velocity_1: float | np.ndarray, velocity_2: float | np.ndarray
    ) -> float | np.ndarray:
        """
        Calculate the average relative velocity between two particles.

        Parameters
        ----------
        velocity_1:
            Velocity of particle 1 (m/s).
        velocity_2:
            Velocity of particle 2 (m/s).

        Returns
        -------
            Average relative velocity (m/s).
        """
        return np.sqrt(velocity_1**2 + velocity_2**2)

    def output_detailed_physics(self):
        """Outputs detailed physics variables to file."""

        po.oheadr(self.outfile, "Detailed Plasma")

        po.osubhd(self.outfile, "Debye lengths:")

        po.ovarrf(
            self.outfile,
            "Plasma volume averaged electron Debye length (m)",
            "(len_plasma_debye_electron_vol_avg)",
            physics_variables.len_plasma_debye_electron_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.len_plasma_debye_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron Debye length at point {i}",
                f"(len_plasma_debye_electron_profile{i})",
                physics_variables.len_plasma_debye_electron_profile[i],
            )

        po.osubhd(self.outfile, "Velocities:")

        for i in range(len(physics_variables.vel_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron thermal velocity at point {i}",
                f"(vel_plasma_electron_profile{i})",
                physics_variables.vel_plasma_electron_profile[i],
            )
        for i in range(len(physics_variables.vel_plasma_deuteron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma deuteron thermal velocity at point {i}",
                f"(vel_plasma_deuteron_profile{i})",
                physics_variables.vel_plasma_deuteron_profile[i],
            )
        for i in range(len(physics_variables.vel_plasma_triton_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma triton thermal velocity at point {i}",
                f"(vel_plasma_triton_profile{i})",
                physics_variables.vel_plasma_triton_profile[i],
            )

        po.osubhd(self.outfile, "Frequencies:")

        for i in range(len(physics_variables.freq_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron frequency at point {i}",
                f"(freq_plasma_electron_profile{i})",
                physics_variables.freq_plasma_electron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_electron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma electron toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_electron_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_electron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_deuteron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma deuteron toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_deuteron_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_deuteron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_triton_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma triton toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_triton_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_triton_profile[i],
            )

        po.osubhd(self.outfile, "Coulomb Logarithms:")

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_electron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-electron Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_electron_profile{i})",
                physics_variables.plasma_coulomb_log_electron_electron_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_deuteron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-deuteron Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_deuteron_profile{i})",
                physics_variables.plasma_coulomb_log_electron_deuteron_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_triton_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-triton Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_triton_profile{i})",
                physics_variables.plasma_coulomb_log_electron_triton_profile[i],
            )
