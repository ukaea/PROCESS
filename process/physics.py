import math

import numba as nb
import numpy as np
import scipy
import scipy.integrate as integrate
from scipy.optimize import root_scalar

import process.confinement_time as confinement
import process.fusion_reactions as reactions
import process.impurity_radiation as impurity_radiation
import process.l_h_transition as transition
import process.physics_functions as physics_funcs
from process import (
    process_output as po,
)
from process.data_structure import (
    divertor_variables,
    pulse_variables,
    reinke_variables,
    times_variables,
)
from process.exceptions import ProcessValueError
from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    current_drive_variables,
    error_handling,
    fwbs_variables,
    impurity_radiation_module,
    numerics,
    physics_module,
    physics_variables,
    stellarator_variables,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string

ELECTRON_CHARGE: float = constants.electron_charge.item()
RMU0: float = constants.rmu0.item()


@nb.jit(nopython=True, cache=True)
def rether(alphan, alphat, dene, dlamie, te, ti, zeffai):
    """Routine to find the equilibration power between the
    ions and electrons
    author: P J Knight, CCFE, Culham Science Centre
    alphan : input real :  density profile index
    alphat : input real :  temperature profile index
    dene   : input real :  electron density (/m3)
    dlamie : input real :  ion-electron coulomb logarithm
    te     : input real :  electron temperature (keV)
    ti     : input real :  ion temperature (keV)
    zeffai : input real :  mass weighted plasma effective charge
    pden_ion_electron_equilibration_mw  : output real : ion/electron equilibration power (MW/m3)
    This routine calculates the equilibration power between the
    ions and electrons.
    Unknown origin
    """
    profie = (1.0 + alphan) ** 2 / (
        (2.0 * alphan - 0.5 * alphat + 1.0) * np.sqrt(1.0 + alphat)
    )
    conie = 2.42165e-41 * dlamie * dene**2 * zeffai * profie

    return conie * (ti - te) / (te**1.5)


@nb.jit(nopython=True, cache=True)
def calculate_volt_second_requirements(
    csawth: float,
    eps: float,
    f_c_plasma_inductive: float,
    ejima_coeff: float,
    kappa: float,
    rmajor: float,
    res_plasma: float,
    plasma_current: float,
    t_fusion_ramp: float,
    t_burn: float,
    ind_plasma_internal_norm: float,
) -> tuple[float, float, float, float, float, float]:
    """Calculate the volt-second requirements and related parameters for plasma physics.

    :param csawth: Coefficient for sawteeth effects
    :type csawth: float
    :param eps: Inverse aspect ratio
    :type eps: float
    :param f_c_plasma_inductive: Fraction of plasma current produced inductively
    :type f_c_plasma_inductive: float
    :param ejima_coeff: Ejima coefficient for resistive start-up V-s component
    :type ejima_coeff: float
    :param kappa: Plasma elongation
    :type kappa: float
    :param rmajor: Plasma major radius (m)
    :type rmajor: float
    :param res_plasma: Plasma resistance (ohm)
    :type res_plasma: float
    :param plasma_current: Plasma current (A)
    :type plasma_current: float
    :param t_fusion_ramp: Heating time (s)
    :type t_fusion_ramp: float
    :param t_burn: Burn time (s)
    :type t_burn: float
    :param ind_plasma_internal_norm: Plasma normalized internal inductance
    :type ind_plasma_internal_norm: float

    :return: A tuple containing:
        - vs_plasma_internal: Internal plasma volt-seconds (Wb)
        - ind_plasma_internal: Plasma inductance (H)
        - vs_plasma_burn_required: Volt-seconds needed during flat-top (heat+burn) (Wb)
        - vs_plasma_ramp_required: Volt-seconds needed during ramp-up (Wb)
        - ind_plasma_total,: Internal and external plasma inductance V-s (Wb)
        - vs_res_ramp: Resistive losses in start-up volt-seconds (Wb)
        - vs_plasma_total_required: Total volt-seconds needed (Wb)
    :rtype: tuple[float, float, float, float, float, float]

    :notes:

    :references:
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

    ind_plasma_internal = RMU0 * rmajor * ind_plasma_internal_norm / 2.0

    # Internal plasma flux (V-s) component
    vs_plasma_internal = ind_plasma_internal * plasma_current

    # Start-up resistive component
    # Uses ITER formula without the 10 V-s add-on

    vs_res_ramp = ejima_coeff * RMU0 * plasma_current * rmajor

    # ======================================================================

    # Hirshman and Neilson fit for external inductance

    aeps = (1.0 + 1.81 * np.sqrt(eps) + 2.05 * eps) * np.log(8.0 / eps) - (
        2.0 + 9.25 * np.sqrt(eps) - 1.21 * eps
    )
    beps = 0.73 * np.sqrt(eps) * (1.0 + 2.0 * eps**4 - 6.0 * eps**5 + 3.7 * eps**6)

    ind_plasma_external = (
        rmajor * RMU0 * aeps * (1.0 - eps) / (1.0 - eps + beps * kappa)
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

    # N.B. t_burn on first iteration will not be correct
    # if the pulsed reactor option is used, but the value
    # will be correct on subsequent calls.

    vs_plasma_burn_required = v_burn_resistive * (t_fusion_ramp + t_burn)
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


@nb.jit(nopython=True, cache=True)
def calculate_beta_limit(
    bt: float, beta_norm_max: float, plasma_current: float, rminor: float
) -> float:
    """
    Calculate the beta scaling limit.

    :param bt: Toroidal B-field on plasma axis [T].
    :type bt: float
    :param beta_norm_max: Troyon-like g coefficient.
    :type beta_norm_max: float
    :param plasma_current: Plasma current [A].
    :type plasma_current: float
    :param rminor: Plasma minor axis [m].
    :type rminor: float
    :return: Beta limit as defined below.
    :rtype: float

    This subroutine calculates the beta limit using the algorithm documented in AEA FUS 172.
    The limit applies to beta defined with respect to the total B-field.
    Switch i_beta_component determines which components of beta to include.

    Notes:
        - If i_beta_component = 0, then the limit is applied to the total beta.
        - If i_beta_component = 1, then the limit is applied to the thermal beta only.
        - If i_beta_component = 2, then the limit is applied to the thermal + neutral beam beta components.
        - If i_beta_component = 3, then the limit is applied to the toroidal beta.

        - The default value for the g coefficient is beta_norm_max = 3.5.

    References:
        - F. Troyon et.al,  “Beta limit in tokamaks. Experimental and computational status,”
          Plasma Physics and Controlled Fusion, vol. 30, no. 11, pp. 1597-1609, Oct. 1988,
          doi: https://doi.org/10.1088/0741-3335/30/11/019.

        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

    """

    # Multiplied by 0.01 to convert from % to fraction
    return 0.01 * beta_norm_max * (plasma_current / 1.0e6) / (rminor * bt)


# -----------------------------------------------------
# Plasma Current & Poloidal Field Calculations
# -----------------------------------------------------


@nb.jit(nopython=True, cache=True)
def _plascar_bpol(
    aspect: float, eps: float, kappa: float, delta: float
) -> tuple[float, float, float, float]:
    """
    Calculate the poloidal field coefficients for determining the plasma current
    and poloidal field.

    Parameters:
    - aspect: float, plasma aspect ratio
    - eps: float, inverse aspect ratio
    - kappa: float, plasma elongation
    - delta: float, plasma triangularity

    Returns:
    - Tuple[float, float, float, float], coefficients ff1, ff2, d1, d2

    This internal function calculates the poloidal field coefficients,
    which is used to calculate the poloidal field and the plasma current.

    References:
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
    bt: float,
    kappa: float,
    delta: float,
    perim: float,
    rmu0: float,
) -> float:
    """
    Function to calculate poloidal field from the plasma current

    Parameters:
    - i_plasma_current: int, current scaling model to use
    - ip: float, plasma current (A)
    - q95: float, 95% flux surface safety factor
    - aspect: float, plasma aspect ratio
    - eps: float, inverse aspect ratio
    - bt: float, toroidal field on axis (T)
    - kappa: float, plasma elongation
    - delta: float, plasma triangularity
    - perim: float, plasma perimeter (m)
    - rmu0: float, vacuum permeability (H/m)

    Returns:
    - float, poloidal field in Tesla

    This function calculates the poloidal field from the plasma current in Tesla,
    using a simple calculation using Ampere's law for conventional
    tokamaks, or for TARTs, a scaling from Peng, Galambos and
    Shipe (1992).

    References:
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

    return bt * (ff1 + ff2) / (2.0 * np.pi * qbar)


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


def calculate_plasma_current_peng(
    q95: float,
    aspect: float,
    eps: float,
    rminor: float,
    bt: float,
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
    - bt: float, toroidal field on axis (T)
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
        * bt
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


@nb.jit(nopython=True, cache=True)
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
    bt: float,
    delta95: float,
    eps: float,
    kappa95: float,
    p0: float,
    rmu0: float,
) -> float:
    """
    Routine to calculate the f_q coefficient for the Connor-Hastie model used for scaling the plasma current.

    Parameters:
    - alphaj: float, the current profile index
    - alphap: float, the pressure profile index
    - bt: float, the toroidal field on axis (T)
    - delta95: float, the plasma triangularity 95%
    - eps: float, the inverse aspect ratio
    - kappa95: float, the plasma elongation 95%
    - p0: float, the central plasma pressure (Pa)
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
    beta0 = 2.0 * rmu0 * p0 / (bt**2)

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


@nb.jit(nopython=True, cache=True)
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


# --------------------------------
# Bootstrap Current Calculations
# --------------------------------


@nb.jit(nopython=True, cache=True)
def _nevins_integral(
    y: float,
    dene: float,
    te: float,
    bt: float,
    rminor: float,
    rmajor: float,
    zeff: float,
    alphat: float,
    alphan: float,
    q0: float,
    q95: float,
    beta_toroidal: float,
) -> float:
    """
    Integrand function for Nevins et al bootstrap current scaling.

    Parameters:
    - y: float, abscissa of integration, normalized minor radius
    - dene: float, volume averaged electron density (/m^3)
    - te: float, volume averaged electron temperature (keV)
    - bt: float, toroidal field on axis (T)
    - rminor: float, plasma minor radius (m)
    - rmajor: float, plasma major radius (m)
    - zeff: float, plasma effective charge
    - alphat: float, temperature profile index
    - alphan: float, density profile index
    - q0: float, normalized safety factor at the magnetic axis
    - q95: float, normalized safety factor at 95% of the plasma radius
    - beta_toroidal: float, Toroidal plasma beta

    Returns:
    - float, the integrand value

    This function calculates the integrand function for the Nevins et al bootstrap current scaling.

    Reference: See appendix of:
        Keii Gi, Makoto Nakamura, Kenji Tobita, Yasushi Ono,
        Bootstrap current fraction scaling for a tokamak reactor design study,
        Fusion Engineering and Design, Volume 89, Issue 11, 2014, Pages 2709-2715,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2014.07.009.

        Nevins, W. M. "Summary report: ITER specialists' meeting on heating and current drive."
        ITER-TN-PH-8-4, June 1988. 1988.

    """

    # Compute average electron beta
    betae = dene * te * 1.0e3 * ELECTRON_CHARGE / (bt**2 / (2.0 * RMU0))

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
    """
    Calculate the diamagnetic fraction based on the Hender fit.

    Parameters:
    - beta: float, the plasma beta value

    Returns:
    - float, the diamagnetic fraction
    """
    return beta / 2.8


@nb.jit(nopython=True, cache=True)
def diamagnetic_fraction_scene(beta: float, q95: float, q0: float) -> float:
    """
    Calculate the diamagnetic fraction based on the SCENE fit by Tim Hender.

    Parameters:
    - beta: float, the plasma beta value
    - q95: float, the normalized safety factor at 95% of the plasma radius
    - q0: float, the normalized safety factor at the magnetic axis

    Returns:
    - float, the diamagnetic fraction
    """
    return beta * (0.1 * q95 / q0 + 0.44) * 0.414


@nb.jit(nopython=True, cache=True)
def ps_fraction_scene(beta: float) -> float:
    """
    Calculate the Pfirsch-Schlüter fraction based on the SCENE fit by Tim Hender 2019.

    Parameters:
    - beta: float, the plasma beta value

    Returns:
    - float, the Pfirsch-Schlüter current fraction
    """
    return -9e-2 * beta


# --------------------------------------------------
# Sauter Bootstrap Current Scaling Functions
# --------------------------------------------------


@nb.jit(nopython=True, cache=True)
def _coulomb_logarithm_sauter(
    radial_elements: int, tempe: np.ndarray, ne: np.ndarray
) -> np.ndarray:
    """
    Calculate the Coulomb logarithm used in the arrays for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, the radial element indexes in the range 1 to nr
    - tempe: np.ndarray, the electron temperature array
    - ne: np.ndarray, the electron density array

    Returns:
    - np.ndarray, the Coulomb logarithm at each array point

    This function calculates the Coulomb logarithm, which is valid for e-e collisions (T_e > 0.01 keV)
    and for e-i collisions (T_e > 0.01*Zeff^2) (Alexander, 9/5/1994).

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


@nb.jit(nopython=True, cache=True)
def _electron_collisions_sauter(
    radial_elements: np.ndarray, tempe: np.ndarray, ne: np.ndarray
) -> np.ndarray:
    """
    Calculate the frequency of electron-electron collisions used in the arrays for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, the radial element index in the range 1 to nr
    - tempe: np.ndarray, the electron temperature array
    - ne: np.ndarray, the electron density array

    Returns:
    - float, the frequency of electron-electron collisions (Hz)

    This function calculates the frequency of electron-electron collisions

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
    """
    Calculate the electron collisionality used in the arrays for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, the radial element index in the range 1 to nr
    - rmajor: float, the plasma major radius (m)
    - zeff: np.ndarray, the effective charge array
    - inverse_q: np.ndarray, inverse safety factor profile
    - sqeps: np.ndarray, the square root of the inverse aspect ratio array
    - tempe: np.ndarray, the electron temperature array
    - ne: np.ndarray, the electron density array

    Returns:
    - np.ndarray, the relative frequency of electron collisions

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
    """
    Calculate the full frequency of ion collisions used in the arrays for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, the radial element indexes in the range 1 to nr
    - zeff: np.ndarray, the effective charge array
    - ni: np.ndarray, the ion density array
    - tempi: np.ndarray, the ion temperature array
    - amain: np.ndarray, the atomic mass of the main ion species array

    Returns:
    - np.ndarray, the full frequency of ion collisions (Hz)

    This function calculates the full frequency of ion collisions using the Coulomb logarithm of 15.

    Reference:
    - None
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
    """
    Calculate the ion collisionality to be used in the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, the radial element indexes in the range 1 to nr
    - rmajor: float, the plasma major radius (m)
    - inverse_q: np.ndarray, inverse safety factor profile
    - sqeps: np.ndarray, the square root of the inverse aspect ratio profile
    - tempi: np.ndarray, the ion temperature profile (keV)
    - amain: np.ndarray, the atomic mass of the main ion species profile
    - zeff: np.ndarray, the effective charge of the main ion species
    - ni: np.ndarray, the ion density profile (/m^3)

    Returns:
    - float, the ion collisionality

    Reference:
    - None
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
    bt: float,
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
    """
    L31 coefficient before Grad(ln(ne)) in the Sauter bootstrap scaling.

    Parameters:
    - radial_elements: np.ndarray, radial element indexes in range 2 to nr
    - number_of_elements: int, maximum value of radial_elements
    - rmajor: float, plasma major radius (m)
    - bt: float, toroidal field on axis (T)
    - triang: float, plasma triangularity
    - ne: np.ndarray, electron density profile (/m^3)
    - ni: np.ndarray, ion density profile (/m^3)
    - tempe: np.ndarray, electron temperature profile (keV)
    - tempi: np.ndarray, ion temperature profile (keV)
    - inverse_q: np.ndarray, inverse safety factor profile
    - rho: np.ndarray, normalized minor radius profile
    - zeff: np.ndarray, effective charge profile
    - sqeps: np.ndarray, square root of inverse aspect ratio profile

    Returns:
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
        bt,
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
    bt: float,
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
    """
    L31 & L32 coefficient before Grad(ln(Te)) in the Sauter bootstrap scaling.

    Parameters:
    - radial_elements: np.ndarray, radial element indexes in range 2 to nr
    - number_of_elements: int, maximum value of radial_elements
    - rmajor: float, plasma major radius (m)
    - bt: float, toroidal field on axis (T)
    - triang: float, plasma triangularity
    - ne: np.ndarray, electron density profile (/m^3)
    - ni: np.ndarray, ion density profile (/m^3)
    - tempe: np.ndarray, electron temperature profile (keV)
    - tempi: np.ndarray, ion temperature profile (keV)
    - inverse_q: np.ndarray, inverse safety factor profile
    - rho: np.ndarray, normalized minor radius profile
    - zeff: np.ndarray, effective charge profile
    - sqeps: np.ndarray, square root of inverse aspect ratio profile

    Returns:
    - float, the L31 & L32 coefficient scaling grad(ln(Te)) in the Sauter bootstrap current scaling.

    This function calculates the coefficient scaling grad(ln(Te)) in the Sauter bootstrap current scaling.

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
        radial_elements, number_of_elements, rmajor, bt, ne, tempe, inverse_q, rho
    ) * (big_f32ee_teff + big_f32ei_teff) + _calculate_l31_coefficient(
        radial_elements,
        number_of_elements,
        rmajor,
        bt,
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
        radial_elements, number_of_elements, rmajor, bt, ne, tempe, inverse_q, rho
    ) / _beta_poloidal_total_sauter(
        radial_elements,
        number_of_elements,
        rmajor,
        bt,
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
    bt: float,
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
    """
    L34, alpha and L31 coefficient before Grad(ln(Ti)) in the Sauter bootstrap scaling.

    Parameters:
    - radial_elements: np.ndarray, radial element indexes in range 2 to nr
    - number_of_elements: int, maximum value of radial_elements
    - rmajor: float, plasma major radius (m)
    - bt: float, toroidal field on axis (T)
    - triang: float, plasma triangularity
    - inverse_q: np.ndarray, inverse safety factor profile
    - sqeps: np.ndarray, square root of inverse aspect ratio profile
    - tempi: np.ndarray, ion temperature profile (keV)
    - tempe: np.ndarray, electron temperature profile (keV)
    - amain: float, atomic mass of the main ion
    - zmain: float, charge of the main ion
    - ni: np.ndarray, ion density profile (/m^3)
    - ne: np.ndarray, electron density profile (/m^3)
    - rho: np.ndarray, normalized minor radius profile
    - zeff: np.ndarray, effective charge profile

    Returns:
    - float, the the L34, alpha and L31 coefficient scaling grad(ln(Ti)) in the Sauter bootstrap current scaling.

    This function calculates the coefficient scaling grad(ln(Ti)) in the Sauter bootstrap current scaling.

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
            bt,
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
            bt,
            ne,
            tempe,
            inverse_q,
            rho,
        )
    ) * (l34_coefficient * alpha) + _calculate_l31_coefficient(
        radial_elements,
        number_of_elements,
        rmajor,
        bt,
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
            bt,
            ne,
            tempe,
            inverse_q,
            rho,
        )
        / _beta_poloidal_total_sauter(
            radial_elements,
            number_of_elements,
            rmajor,
            bt,
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
    bt: float,
    ne: np.ndarray,
    tempe: np.ndarray,
    inverse_q: np.ndarray,
    rho: np.ndarray,
) -> np.ndarray:
    """
    Calculate the local beta poloidal using only electron profiles for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, radial element indexes in range 1 to nr
    - nr: int, maximum value of radial_elements
    - rmajor: float, plasma major radius (m)
    - bt: float, toroidal field on axis (T)
    - ne: np.ndarray, electron density profile (/m^3)
    - tempe: np.ndarray, electron temperature profile (keV)
    - inverse_q: np.ndarray, inverse safety factor profile
    - rho: np.ndarray, normalized minor radius profile

    Returns:
    - np.ndarray, the local beta poloidal

    Reference:
    - None
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
                bt
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
    bt: float,
    ne: np.ndarray,
    ni: np.ndarray,
    tempe: np.ndarray,
    tempi: np.ndarray,
    inverse_q: np.ndarray,
    rho: np.ndarray,
) -> np.ndarray:
    """
    Calculate the local beta poloidal including ion pressure for the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: np.ndarray, radial element indexes in range 1 to nr
    - nr: int, maximum value of radial_elements
    - rmajor: float, plasma major radius (m)
    - bt: float, toroidal field on axis (T)
    - ne: np.ndarray, electron density profile (/m^3)
    - ni: np.ndarray, ion density profile (/m^3)
    - tempe: np.ndarray, electron temperature profile (keV)
    - tempi: np.ndarray, ion temperature profile (keV)
    - inverse_q: np.ndarray, inverse safety factor profile
    - rho: np.ndarray, normalized minor radius profile

    Returns:
    - np.ndarray, the local total beta poloidal

    Reference:
    - None
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
                bt
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
    """
    Calculates the trapped particle fraction to be used in the Sauter bootstrap current scaling.

    Parameters:
    - radial_elements: list, radial element index in range 1 to nr
    - triang: float, plasma triangularity
    - sqeps: list, square root of local aspect ratio
    - fit: int, fit method (1 = ASTRA method, 2 = Equation from Sauter 2002, 3 = Equation from Sauter 2016)

    Returns:
    - list, trapped particle fraction

    This function calculates the trapped particle fraction at a given radius.

    References:
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
    def __init__(self, plasma_profile, current_drive):
        self.outfile = constants.nout
        self.plasma_profile = plasma_profile
        self.current_drive = current_drive

    def physics(self):
        """
        Routine to calculate tokamak plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine calculates all the primary plasma physics parameters for a tokamak fusion reactor.

        :References:
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
            physics_variables.nd_ions_total,
            physics_variables.nd_fuel_ions,
            physics_variables.nd_alphas,
            physics_variables.vol_plasma,
            physics_variables.dene,
        )

        # Define coulomb logarithm
        # (collisions: ion-electron, electron-electron)
        physics_variables.dlamee = (
            31.0
            - (np.log(physics_variables.dene) / 2.0)
            + np.log(physics_variables.te * 1000.0)
        )
        physics_variables.dlamie = (
            31.3
            - (np.log(physics_variables.dene) / 2.0)
            + np.log(physics_variables.te * 1000.0)
        )

        # Calculate plasma current
        (
            physics_variables.bp,
            physics_variables.qstar,
            physics_variables.plasma_current,
        ) = self.calculate_plasma_current(
            physics_variables.alphaj,
            physics_variables.alphap,
            physics_variables.bt,
            physics_variables.eps,
            physics_variables.i_plasma_current,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.p0,
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

        physics_variables.alphaj_wesson = (
            self.calculate_current_profile_index_wesson(
                qstar=physics_variables.qstar, q0=physics_variables.q0
            ),
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

        physics_variables.ind_plasma_internal_norm_wesson = (
            self.calculate_internal_inductance_wesson(alphaj=physics_variables.alphaj)
        )

        # Spherical Tokamak relation for internal inductance
        # Menard et al. (2016), Nuclear Fusion, 56, 106023
        physics_variables.ind_plasma_internal_norm_menard = (
            self.calculate_internal_inductance_menard(kappa=physics_variables.kappa)
        )

        # Map calculation methods to a dictionary
        ind_plasma_internal_norm_calculations = {
            0: physics_variables.ind_plasma_internal_norm,
            1: physics_variables.ind_plasma_internal_norm_wesson,
            2: physics_variables.ind_plasma_internal_norm_menard,
        }

        # Calculate ind_plasma_internal_normx based on i_ind_plasma_internal_norm
        if (
            int(physics_variables.i_ind_plasma_internal_norm)
            in ind_plasma_internal_norm_calculations
        ):
            physics_variables.ind_plasma_internal_norm = (
                ind_plasma_internal_norm_calculations[
                    int(physics_variables.i_ind_plasma_internal_norm)
                ]
            )
        else:
            raise ProcessValueError(
                "Illegal value of i_ind_plasma_internal_norm",
                i_ind_plasma_internal_norm=physics_variables.i_ind_plasma_internal_norm,
            )

        # ===================================================

        # Calculate density and temperature profile quantities
        # If physics_variables.ipedestal = 1 then set pedestal density to
        #   physics_variables.fgwped * Greenwald density limit
        # Note: this used to be done before plasma current
        if (physics_variables.ipedestal == 1) and (physics_variables.fgwped >= 0e0):
            physics_variables.neped = (
                physics_variables.fgwped
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        if (physics_variables.ipedestal == 1) and (physics_variables.fgwsep >= 0e0):
            physics_variables.nesep = (
                physics_variables.fgwsep
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        self.plasma_profile.run()

        # Calculate total magnetic field [T]
        physics_variables.btot = np.sqrt(
            physics_variables.bt**2 + physics_variables.bp**2
        )

        # ============================================

        # -----------------------------------------------------
        # Beta Components
        # -----------------------------------------------------

        physics_variables.beta_toroidal = (
            physics_variables.beta * physics_variables.btot**2 / physics_variables.bt**2
        )

        # Calculate physics_variables.beta poloidal [-]
        physics_variables.beta_poloidal = calculate_poloidal_beta(
            physics_variables.btot, physics_variables.bp, physics_variables.beta
        )

        physics_variables.beta_thermal = (
            physics_variables.beta
            - physics_variables.beta_fast_alpha
            - physics_variables.beta_beam
        )

        physics_variables.beta_poloidal_eps = (
            physics_variables.beta_poloidal * physics_variables.eps
        )

        physics_variables.beta_thermal_poloidal = (
            physics_variables.beta_thermal
            * (physics_variables.btot / physics_variables.bp) ** 2
        )
        physics_variables.beta_thermal_toroidal = (
            physics_variables.beta_thermal
            * (physics_variables.btot / physics_variables.bt) ** 2
        )
        physics_variables.beta_norm_thermal = (
            1.0e8
            * physics_variables.beta_thermal
            * physics_variables.rminor
            * physics_variables.bt
            / physics_variables.plasma_current,
        )
        physics_variables.beta_norm_toroidal = (
            physics_variables.beta_norm_total
            * (physics_variables.btot / physics_variables.bt) ** 2
        )
        physics_variables.beta_norm_poloidal = (
            physics_variables.beta_norm_total
            * (physics_variables.btot / physics_variables.bp) ** 2
        )

        physics_variables.f_beta_alpha_beam_thermal = (
            physics_variables.beta_fast_alpha + physics_variables.beta_beam
        ) / physics_variables.beta_thermal

        # Plasma thermal energy derived from the thermal beta
        physics_variables.e_plasma_beta_thermal = (
            1.5e0
            * physics_variables.beta_thermal
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol_plasma
        )

        # Plasma thermal energy derived from the total beta
        physics_module.e_plasma_beta = (
            1.5e0
            * physics_variables.beta
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol_plasma
        )

        # =======================================================

        # Set PF coil ramp times
        if pulse_variables.i_pulsed_plant != 1:
            if times_variables.i_t_current_ramp_up == 0:
                times_variables.t_current_ramp_up = (
                    physics_variables.plasma_current / 5.0e5
                )
                times_variables.t_precharge = times_variables.t_current_ramp_up
                times_variables.t_ramp_down = times_variables.t_current_ramp_up

        else:
            if times_variables.pulsetimings == 0.0e0:
                # times_variables.t_precharge is input
                times_variables.t_current_ramp_up = (
                    physics_variables.plasma_current / 1.0e5
                )
                times_variables.t_ramp_down = times_variables.t_current_ramp_up

            else:
                # times_variables.t_current_ramp_up is set either in INITIAL or INPUT, or by being
                # iterated using limit equation 41.
                times_variables.t_precharge = max(
                    times_variables.t_precharge, times_variables.t_current_ramp_up
                )
                # t_ramp_down = max(t_ramp_down,t_current_ramp_up)
                times_variables.t_ramp_down = times_variables.t_current_ramp_up

        # Reset second times_variables.t_burn value (times_variables.t_burn_0).
        # This is used to ensure that the burn time is used consistently;
        # see convergence loop in fcnvmc1, evaluators.f90
        times_variables.t_burn_0 = times_variables.t_burn

        # Pulse and down times : The reactor is assumed to be 'down'
        # at all times outside of the plasma current flat-top period.
        # The pulse length is the duration of non-zero plasma current
        times_variables.t_pulse_repetition = (
            times_variables.t_current_ramp_up
            + times_variables.t_fusion_ramp
            + times_variables.t_burn
            + times_variables.t_ramp_down
        )
        times_variables.tdown = (
            times_variables.t_precharge
            + times_variables.t_current_ramp_up
            + times_variables.t_ramp_down
            + times_variables.t_between_pulse
        )

        # Total cycle time
        times_variables.t_cycle = (
            times_variables.t_precharge
            + times_variables.t_current_ramp_up
            + times_variables.t_fusion_ramp
            + times_variables.t_burn
            + times_variables.t_ramp_down
            + times_variables.t_between_pulse
        )

        # ***************************** #
        #      DIAMAGNETIC CURRENT      #
        # ***************************** #

        # Hender scaling for diamagnetic current at tight physics_variables.aspect ratio
        current_drive_variables.f_c_plasma_diamagnetic_hender = (
            diamagnetic_fraction_hender(physics_variables.beta)
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.f_c_plasma_diamagnetic_scene = (
            diamagnetic_fraction_scene(
                physics_variables.beta, physics_variables.q95, physics_variables.q0
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
            physics_variables.beta
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
                physics_variables.beta,
                physics_variables.btot,
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
                physics_variables.beta_toroidal,
                physics_variables.bt,
                physics_variables.dene,
                physics_variables.plasma_current,
                physics_variables.q95,
                physics_variables.q0,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.te,
                physics_variables.zeff,
            )
        )

        # Wilson scaling uses thermal poloidal beta, not total
        current_drive_variables.f_c_plasma_bootstrap_wilson = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wilson(
                physics_variables.alphaj,
                physics_variables.alphap,
                physics_variables.alphat,
                physics_variables.beta_thermal_poloidal,
                physics_variables.q0,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.rminor,
            )
        )

        current_drive_variables.f_c_plasma_bootstrap_sauter = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sauter(self.plasma_profile)
        )

        current_drive_variables.f_c_plasma_bootstrap_sakai = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sakai(
                beta_poloidal=physics_variables.beta_poloidal,
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
                beta_poloidal=physics_variables.beta_poloidal,
                ind_plasma_internal_norm=physics_variables.ind_plasma_internal_norm,
                core_density=physics_variables.ne0,
                average_density=physics_variables.dene,
                inverse_aspect=physics_variables.eps,
            )
        )

        current_drive_variables.f_c_plasma_bootstrap_andrade = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_andrade(
                beta_poloidal=physics_variables.beta_poloidal,
                core_pressure=physics_variables.p0,
                average_pressure=physics_variables.vol_avg_pressure,
                inverse_aspect=physics_variables.eps,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_hoang = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_hoang(
                beta_poloidal=physics_variables.beta_poloidal,
                pressure_index=physics_variables.alphap,
                current_index=physics_variables.alphaj,
                inverse_aspect=physics_variables.eps,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_wong = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wong(
                beta_poloidal=physics_variables.beta_poloidal,
                density_index=physics_variables.alphan,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                elongation=physics_variables.kappa,
            )
        )
        current_drive_variables.bscf_gi_I = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_gi_I(
                beta_poloidal=physics_variables.beta_poloidal,
                pressure_index=physics_variables.alphap,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                effective_charge=physics_variables.zeff,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
            )
        )

        current_drive_variables.bscf_gi_II = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_gi_II(
                beta_poloidal=physics_variables.beta_poloidal,
                pressure_index=physics_variables.alphap,
                temperature_index=physics_variables.alphat,
                inverse_aspect=physics_variables.eps,
                effective_charge=physics_variables.zeff,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_sugiyama_l = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sugiyama_l_mode(
                eps=physics_variables.eps,
                beta_poloidal=physics_variables.beta_poloidal,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                zeff=physics_variables.zeff,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
            )
        )
        current_drive_variables.f_c_plasma_bootstrap_sugiyama_h = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sugiyama_h_mode(
                eps=physics_variables.eps,
                beta_poloidal=physics_variables.beta_poloidal,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                tbeta=physics_variables.tbeta,
                zeff=physics_variables.zeff,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
                rhopedn=physics_variables.rhopedn,
                neped=physics_variables.neped,
                n_greenwald=physics_variables.dlimit[6],
                teped=physics_variables.teped,
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
            10: current_drive_variables.bscf_gi_I,
            11: current_drive_variables.bscf_gi_II,
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

        physics_module.err242 = 0
        if (
            current_drive_variables.f_c_plasma_bootstrap
            > current_drive_variables.f_c_plasma_bootstrap_max
        ) and physics_variables.i_bootstrap_current != 0:
            current_drive_variables.f_c_plasma_bootstrap = min(
                current_drive_variables.f_c_plasma_bootstrap,
                current_drive_variables.f_c_plasma_bootstrap_max,
            )
            physics_module.err242 = 1

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
        physics_module.err243 = 0
        if (
            current_drive_variables.f_c_plasma_internal
            > physics_variables.f_c_plasma_non_inductive
        ):
            current_drive_variables.f_c_plasma_internal = min(
                current_drive_variables.f_c_plasma_internal,
                physics_variables.f_c_plasma_non_inductive,
            )
            physics_module.err243 = 1

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
        fusion_reactions.deuterium_branching(physics_variables.ti)
        fusion_reactions.calculate_fusion_rates()
        fusion_reactions.set_physics_variables()

        # This neglects the power from the beam
        physics_variables.p_plasma_dt_mw = (
            physics_module.dt_power_density_plasma * physics_variables.vol_plasma
        )
        physics_variables.p_dhe3_total_mw = (
            physics_module.dhe3_power_density * physics_variables.vol_plasma
        )
        physics_variables.p_dd_total_mw = (
            physics_module.dd_power_density * physics_variables.vol_plasma
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
                physics_variables.bp,
                physics_variables.bt,
                current_drive_variables.c_beam_total,
                physics_variables.dene,
                physics_variables.nd_fuel_ions,
                physics_variables.dlamie,
                current_drive_variables.e_beam_kev,
                physics_variables.f_deuterium,
                physics_variables.f_tritium,
                current_drive_variables.f_beam_tritium,
                physics_module.sigmav_dt_average,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.vol_plasma,
                physics_variables.zeffai,
            )
            physics_variables.fusden_total = (
                physics_variables.fusden_plasma
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.dt_alpha_energy)
                / physics_variables.vol_plasma
            )
            physics_variables.fusden_alpha_total = (
                physics_variables.fusden_plasma_alpha
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.dt_alpha_energy)
                / physics_variables.vol_plasma
            )
            physics_variables.p_dt_total_mw = (
                physics_variables.p_plasma_dt_mw
                + (1.0 / (1.0 - constants.dt_neutron_energy_fraction))
                * physics_variables.p_beam_alpha_mw
            )
            physics_variables.p_beam_neutron_mw = physics_variables.p_beam_alpha_mw * (
                constants.dt_neutron_energy_fraction
                / (1 - constants.dt_neutron_energy_fraction)
            )
            physics_variables.p_beam_dt_mw = physics_variables.p_beam_alpha_mw * (
                1 / (1 - constants.dt_neutron_energy_fraction)
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

        physics_variables.beta_fast_alpha = physics_funcs.fast_alpha_beta(
            physics_variables.bp,
            physics_variables.bt,
            physics_variables.dene,
            physics_variables.nd_fuel_ions,
            physics_variables.nd_ions_total,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.pden_alpha_total_mw,
            physics_variables.pden_plasma_alpha_mw,
            physics_variables.i_beta_fast_alpha,
        )

        # Nominal mean neutron wall load on entire first wall area including divertor and beam holes
        # Note that 'a_fw_total' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.ffwal
                * physics_variables.p_neutron_total_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if physics_variables.n_divertors == 2:
                # Double null configuration
                physics_variables.pflux_fw_neutron_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_hcd
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
                        - fwbs_variables.f_a_fw_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_neutron_total_mw
                    / build_variables.a_fw_total
                )

        # Calculate ion/electron equilibration power

        physics_variables.pden_ion_electron_equilibration_mw = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.dene,
            physics_variables.dlamie,
            physics_variables.te,
            physics_variables.ti,
            physics_variables.zeffai,
        )

        # Calculate radiation power

        radpwrdata = physics_funcs.calculate_radiation_powers(
            self.plasma_profile,
            physics_variables.ne0,
            physics_variables.rminor,
            physics_variables.bt,
            physics_variables.aspect,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.tbeta,
            physics_variables.te0,
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
            physics_variables.ten,
            physics_variables.vol_plasma,
            physics_variables.zeff,
        )

        # Calculate L- to H-mode power threshold for different scalings
        physics_variables.l_h_threshold_powers = l_h_threshold_power(
            physics_variables.nd_electron_line,
            physics_variables.bt,
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
        if physics_variables.n_divertors == 2:
            physics_variables.pdivl = (
                physics_variables.f_p_div_lower
                * physics_variables.p_plasma_separatrix_mw
            )
            physics_variables.pdivu = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_variables.p_plasma_separatrix_mw
            physics_variables.pdivmax = max(
                physics_variables.pdivl, physics_variables.pdivu
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
            physics_variables.dlimit,
            physics_variables.dnelimt,
        ) = self.calculate_density_limit(
            physics_variables.bt,
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
            physics_variables.zeff,
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
            physics_variables.bt,
            physics_variables.nd_ions_total,
            physics_variables.dene,
            physics_variables.nd_electron_line,
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
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.q95,
            physics_variables.qstar,
            physics_variables.vol_plasma,
            physics_variables.zeff,
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
        ) = calculate_volt_second_requirements(
            physics_variables.csawth,
            physics_variables.eps,
            physics_variables.f_c_plasma_inductive,
            physics_variables.ejima_coeff,
            physics_variables.kappa,
            physics_variables.rmajor,
            physics_variables.res_plasma,
            physics_variables.plasma_current,
            times_variables.t_fusion_ramp,
            times_variables.t_burn,
            physics_variables.ind_plasma_internal_norm,
        )

        # Calculate auxiliary physics related information
        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.ntau,
            physics_variables.nTtau,
            physics_variables.figmer,
            physics_module.fusrat,
            physics_variables.qfuel,
            physics_variables.rndfuel,
            physics_variables.t_alpha_confinement,
            physics_variables.f_alpha_energy_confinement,
        ) = self.phyaux(
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.te,
            physics_variables.nd_fuel_ions,
            physics_variables.fusden_total,
            physics_variables.fusden_alpha_total,
            physics_variables.plasma_current,
            sbar,
            physics_variables.nd_alphas,
            physics_variables.t_energy_confinement,
            physics_variables.vol_plasma,
        )

        # Total transport power from scaling law (MW)
        physics_variables.pscalingmw = (
            physics_variables.p_electron_transport_loss_mw
            + physics_variables.p_ion_transport_loss_mw
        )

        # ============================================================

        # -----------------------------------------------------
        # Normalised Beta Limit
        # -----------------------------------------------------

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
        physics_variables.beta_norm_max_thloreus = (
            self.calculate_beta_norm_max_thloreus(
                c_beta=physics_variables.c_beta,
                p0=physics_variables.p0,
                vol_avg_pressure=physics_variables.vol_avg_pressure,
            )
        )

        # R. D. Stambaugh scaling law
        physics_variables.beta_norm_max_stambaugh = (
            self.calculate_beta_norm_max_stambaugh(
                f_c_plasma_bootstrap=current_drive_variables.f_c_plasma_bootstrap,
                kappa=physics_variables.kappa,
                aspect=physics_variables.aspect,
            )
        )

        # Map calculation methods to a dictionary
        beta_norm_max_calculations = {
            0: physics_variables.beta_norm_max,
            1: physics_variables.beta_norm_max_wesson,
            2: physics_variables.beta_norm_max_original_scaling,
            3: physics_variables.beta_norm_max_menard,
            4: physics_variables.beta_norm_max_thloreus,
            5: physics_variables.beta_norm_max_stambaugh,
        }

        # Calculate beta_norm_max based on i_beta_norm_max
        if int(physics_variables.i_beta_norm_max) in beta_norm_max_calculations:
            physics_variables.beta_norm_max = beta_norm_max_calculations[
                int(physics_variables.i_beta_norm_max)
            ]
        else:
            raise ProcessValueError(
                "Illegal value of i_beta_norm_max",
                i_beta_norm_max=physics_variables.i_beta_norm_max,
            )

        # calculate_beta_limit() returns the beta_max for beta
        physics_variables.beta_max = calculate_beta_limit(
            physics_variables.bt,
            physics_variables.beta_norm_max,
            physics_variables.plasma_current,
            physics_variables.rminor,
        )

        # ============================================================

        # MDK
        # Nominal mean photon wall load on entire first wall area including divertor and beam holes
        # Note that 'a_fw_total' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.pflux_fw_rad_mw = (
                physics_variables.ffwal
                * physics_variables.p_plasma_rad_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if physics_variables.n_divertors == 2:
                # Double Null configuration in - including SoL radiation
                physics_variables.pflux_fw_rad_mw = (
                    (
                        1.0e0
                        - fwbs_variables.f_a_fw_hcd
                        - 2.0e0 * fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / build_variables.a_fw_total
                    + (
                        1.0e0
                        - fwbs_variables.f_a_fw_hcd
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
                        - fwbs_variables.f_a_fw_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / build_variables.a_fw_total
                    + (
                        1.0e0
                        - fwbs_variables.f_a_fw_hcd
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
        physics_module.ptarmw = physics_variables.p_plasma_separatrix_mw * (
            1.0e0 - physics_variables.rad_fraction_sol
        )
        # use physics_variables.f_p_div_lower to find deltarsep
        # Parameters taken from double null machine
        # D. Brunner et al
        physics_module.lambdaio = 1.57e-3

        # Issue #1559 Infinities in physics_module.drsep when running single null in a double null machine
        # C W Ashe
        if physics_variables.f_p_div_lower < 4.5e-5:
            physics_module.drsep = 1.5e-2
        elif physics_variables.f_p_div_lower > (1.0e0 - 4.5e-5):
            physics_module.drsep = -1.5e-2
        else:
            physics_module.drsep = (
                -2.0e0
                * 1.5e-3
                * math.atanh(2.0e0 * (physics_variables.f_p_div_lower - 0.5e0))
            )
        # Model Taken from D3-D paper for conventional divertor
        # Journal of Nuclear Materials
        # Volumes 290-293, March 2001, Pages 935-939
        # Find the innner and outer lower target imbalance

        physics_module.fio = 0.16e0 + (0.16e0 - 0.41e0) * (
            1.0e0
            - (
                2.0e0
                / (
                    1.0e0
                    + np.exp(-((physics_module.drsep / physics_module.lambdaio) ** 2))
                )
            )
        )
        if physics_variables.n_divertors == 2:
            # Double Null configuration
            # Find all the power fractions accross the targets
            # Taken from D3-D conventional divertor design
            physics_module.fli = physics_variables.f_p_div_lower * physics_module.fio
            physics_module.flo = physics_variables.f_p_div_lower * (
                1.0e0 - physics_module.fio
            )
            physics_module.fui = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_module.fio
            physics_module.fuo = (1.0e0 - physics_variables.f_p_div_lower) * (
                1.0e0 - physics_module.fio
            )
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw
            physics_module.puimw = physics_module.fui * physics_module.ptarmw
            physics_module.puomw = physics_module.fuo * physics_module.ptarmw
        else:
            # Single null configuration
            physics_module.fli = physics_module.fio
            physics_module.flo = 1.0e0 - physics_module.fio
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw

        # Calculate some derived quantities that may not have been defined earlier
        physics_module.total_loss_power = 1e6 * (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            + current_drive_variables.p_hcd_injected_total_mw
        )
        physics_module.rad_fraction_lcfs = (
            1.0e6 * physics_variables.p_plasma_rad_mw / physics_module.total_loss_power
        )
        physics_variables.rad_fraction_total = (
            physics_module.rad_fraction_lcfs
            + (1.0e0 - physics_module.rad_fraction_lcfs)
            * physics_variables.rad_fraction_sol
        )
        physics_variables.pradsolmw = (
            physics_variables.rad_fraction_sol
            * physics_variables.p_plasma_separatrix_mw
        )

        if any(numerics.icc == 78):
            po.write(
                self.outfile,
                (
                    f"reinke t and fz, physics = {physics_variables.tesep} , {reinke_variables.fzmin}"
                ),
            )
            # fsep = physics_variables.nesep / physics_variables.dene
            fgw = physics_variables.dlimit(7) / physics_variables.dene
            # calculate separatrix temperature, if Reinke criterion is used
            physics_variables.tesep = reinke_tsep(
                physics_variables.bt,
                constraint_variables.fl_h_threshold,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.eps,
                fgw,
                physics_variables.kappa,
                reinke_variables.lhat,
            )

            if reinke_variables.fzmin >= 1.0e0:
                error_handling.report_error(217)

            po.write(
                self.outfile,
                (
                    f" 'fzactual, frac, reinke_variables.impvardiv = {reinke_variables.fzactual},"
                    f" {impurity_radiation_module.impurity_arr_frac(reinke_variables.impvardiv)},"
                    f" {reinke_variables.impvardiv}"
                ),
            )

    @staticmethod
    def calculate_current_profile_index_wesson(qstar: float, q0: float) -> float:
        """
        Calculate the Wesson current profile index.

        :param qstar: Cylindrical safety factor.
        :type qstar: float
        :param q0: Safety factor on axis.
        :type q0: float


        :return: The Wesson current profile index.
        :rtype: float

        :Notes:
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              normalised internal inductance.
            - This relation is only true for the cyclindrical plasma approximation.

        :References:
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.
        """
        return qstar / q0 - 1.0

    @staticmethod
    def calculate_internal_inductance_wesson(alphaj: float) -> float:
        """
        Calculate the Wesson plasma normalized internal inductance.

        :param alphaj: Current profile index.
        :type alphaj: float

        :return: The Wesson plasma normalised internal inductance.
        :rtype: float

        :Notes:
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              current profile index.
            - This relation is only true for the cyclindrical plasma approximation with parabolic profiles.

        :References:
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.
        """
        return np.log(1.65 + 0.89 * alphaj)

    @staticmethod
    def calculate_beta_norm_max_wesson(ind_plasma_internal_norm: float) -> float:
        """
        Calculate the Wesson normalsied beta upper limit.

        :param ind_plasma_internal_norm: Plasma normalised internal inductance
        :type ind_plasma_internal_norm: float

        :return: The Wesson normalised beta upper limit.
        :rtype: float

        :Notes:
            - It is recommended to use this method with the other Wesson relations for normalsied internal
            inductance and current profile index.
            - This fit is derived from the DIII-D database for β_N >= 2.5

        :References:
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.

            - T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,”
            Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).
        """
        return 4 * ind_plasma_internal_norm

    @staticmethod
    def calculate_beta_norm_max_original(eps: float) -> float:
        """
        Calculate the original scaling law normalsied beta upper limit.

        :param eps: Plasma normalised internal inductance
        :type eps: float

        :return: The original scaling law normalised beta upper limit.
        :rtype: float

        :Notes:

        :References:

        """
        return 2.7 * (1.0 + 5.0 * eps**3.5)

    @staticmethod
    def calculate_beta_norm_max_menard(eps: float) -> float:
        """
        Calculate the Menard normalsied beta upper limit.

        :param eps: Plasma normalised internal inductance
        :type eps: float

        :return: The Menard normalised beta upper limit.
        :rtype: float

        :Notes:
            - Found as a reasonable fit to the computed no wall limit at f_BS ≈ 50%
            - Uses maximum κ data from NSTX at A = 1.45, A = 1.75. Along with record
              β_T data from DIII-D at A = 2.9 and high κ.

        :References:
            - # J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.

        """
        return 3.12 + 3.5 * eps**1.7

    @staticmethod
    def calculate_beta_norm_max_thloreus(
        c_beta: float, p0: float, vol_avg_pressure: float
    ) -> float:
        """
        Calculate the E. Tholerus normalized beta upper limit.

        :param c_beta: Pressure peaking factor coefficient.
        :type c_beta: float
        :param p0: Central plasma pressure (Pa).
        :type p0: float
        :param vol_avg_pressure: Volume-averaged plasma pressure (Pa).
        :type vol_avg_pressure: float

        :return: The E. Tholerus normalized beta upper limit.
        :rtype: float

        :Notes:
            - This method calculates the normalized beta upper limit based on the pressure peaking factor (Fp),
              which is defined as the ratio of the peak pressure to the average pressure.
            - The formula is derived from operational space studies of flat-top plasma in the STEP power plant.

        :References:
            - E. Tholerus et al., “Flat-top plasma operational space of the STEP power plant,”
              Nuclear Fusion, Aug. 2024, doi: https://doi.org/10.1088/1741-4326/ad6ea2.
        """
        return 3.7 + (
            (c_beta / (p0 / vol_avg_pressure)) * (12.5 - 3.5 * (p0 / vol_avg_pressure))
        )

    @staticmethod
    def calculate_beta_norm_max_stambaugh(
        f_c_plasma_bootstrap: float,
        kappa: float,
        aspect: float,
    ) -> float:
        """
        Calculate the Stambaugh normalized beta upper limit.

        :param f_c_plasma_bootstrap: Bootstrap current fraction.
        :type f_c_plasma_bootstrap: float
        :param kappa: Plasma separatrix elongation.
        :type kappa: float
        :param aspect: Plasma aspect ratio.
        :type aspect: float

        :return: The Stambaugh normalized beta upper limit.
        :rtype: float

        :Notes:
            - This method calculates the normalized beta upper limit based on the Stambaugh scaling.
            - The formula is derived from empirical fits to high-performance, steady-state tokamak equilibria.

        :References:
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
    def calculate_internal_inductance_menard(kappa: float) -> float:
        """
        Calculate the Menard plasma normalized internal inductance.

        :param kappa: Plasma separatrix elongation.
        :type kappa: float

        :return: The Menard plasma normalised internal inductance.
        :rtype: float

        :Notes:
            - This relation is based off of data from NSTX for l_i in the range of 0.4-0.85
            - This model is only recommneded to be used for ST's with kappa > 2.5

        :References:
            - J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.
        """
        return 3.4 - kappa

    @staticmethod
    def calculate_density_limit(
        bt: float,
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
        """
        Calculate the density limit using various models.

        Args:
            bt (float): Toroidal field on axis (T).
            i_density_limit (int): Switch denoting which formula to enforce.
            p_plasma_separatrix_mw (float): Power flowing to the edge plasma via charged particles (MW).
            p_hcd_injected_total_mw (float): Power injected into the plasma (MW).
            plasma_current (float): Plasma current (A).
            prn1 (float): Edge density / average plasma density.
            qcyl (float): Equivalent cylindrical safety factor (qstar).
            q95 (float): Safety factor at 95% surface.
            rmajor (float): Plasma major radius (m).
            rminor (float): Plasma minor radius (m).
            a_plasma_surface (float): Plasma surface area (m^2).
            zeff (float): Plasma effective charge.

        Returns:
            Tuple[np.ndarray, float]: A tuple containing:
                - dlimit (np.ndarray): Average plasma density limit using seven different models (m^-3).
                - dnelimt (float): Enforced average plasma density limit (m^-3).

        Raises:
            ValueError: If i_density_limit is not between 1 and 7.

        Notes:
            This routine calculates several different formulae for the density limit and enforces the one chosen by the user.
            For i_density_limit = 1-5, 8, we scale the sepatrix density limit output by the ratio of the separatrix to volume averaged density

        References:
            - AEA FUS 172: Physics Assessment for the European Reactor Study

            - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989

            - M. Bernert et al., “The H-mode density limit in the full tungsten ASDEX Upgrade tokamak,”
              vol. 57, no. 1, pp. 014038-014038, Nov. 2014, doi: https://doi.org/10.1088/0741-3335/57/1/014038. ‌
        """

        if i_density_limit < 1 or i_density_limit > 7:
            raise ProcessValueError(
                "Illegal value for i_density_limit", i_density_limit=i_density_limit
            )

        dlimit = np.empty((8,))

        # Power per unit area crossing the plasma edge
        # (excludes radiation and neutrons)

        p_perp = p_plasma_separatrix_mw / a_plasma_surface

        # Old ASDEX density limit formula
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        dlimit[0] = (1.54e20 * p_perp**0.43 * bt**0.31 / (q95 * rmajor) ** 0.45) / prn1

        # Borrass density limit model for ITER (I)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # Borrass et al, ITER-TN-PH-9-6 (1989)

        dlimit[1] = (1.8e20 * p_perp**0.53 * bt**0.31 / (q95 * rmajor) ** 0.22) / prn1

        # Borrass density limit model for ITER (II)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # This formula is (almost) identical to that in the original routine
        # denlim (now deleted).

        dlimit[2] = (0.5e20 * p_perp**0.57 * bt**0.31 / (q95 * rmajor) ** 0.09) / prn1

        # JET edge radiation density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # qcyl=qstar here, but literature is not clear.

        denom = (zeff - 1.0) * (1.0 - 4.0 / (3.0 * qcyl))
        if denom <= 0.0:
            if i_density_limit == 4:
                error_handling.fdiags[0] = denom
                error_handling.fdiags[1] = qcyl
                error_handling.report_error(80)
                i_density_limit = 5

            dlimit[3] = 0.0
        else:
            dlimit[3] = (1.0e20 * np.sqrt(p_hcd_injected_total_mw / denom)) / prn1

        # JET simplified density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        dlimit[4] = (0.237e20 * bt * np.sqrt(p_plasma_separatrix_mw) / rmajor) / prn1

        # Hugill-Murakami M.q limit
        # qcyl=qstar here, which is okay according to the literature

        dlimit[5] = 3.0e20 * bt / (rmajor * qcyl)

        # Greenwald limit

        dlimit[6] = 1.0e14 * plasma_current / (np.pi * rminor * rminor)

        dlimit[7] = (
            1.0e20
            * 0.506
            * (p_hcd_injected_total_mw**0.396 * (plasma_current / 1.0e6) ** 0.265)
            / (q95**0.323)
        ) / prn1

        # Enforce the chosen density limit

        return dlimit, dlimit[i_density_limit - 1]

    @staticmethod
    def plasma_composition() -> None:
        """
        Calculates various plasma component fractional makeups.

        Author: P J Knight, CCFE, Culham Science Centre

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
        - Sets impurity fractions for radiation calculations.
        - Calculates the effective charge.
        - Defines the Coulomb logarithm for ion-electron and electron-electron collisions.
        - Calculates the fraction of alpha energy to ions and electrons.
        - Calculates the average atomic masses of injected fuel species and neutral beams.
        - Calculates the density weighted mass and mass weighted plasma effective charge.

        Notes:


        References:
        """

        # Alpha ash portion
        physics_variables.nd_alphas = (
            physics_variables.dene * physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Protons
        # This calculation will be wrong on the first call as the particle
        # production rates are evaluated later in the calling sequence
        # Issue #557 Allow f_nd_protium_electrons impurity to be specified: 'f_nd_protium_electrons'
        # This will override the calculated value which is a minimum.
        if physics_variables.fusden_alpha_total < 1.0e-6:  # not calculated yet...
            physics_variables.nd_protons = max(
                physics_variables.f_nd_protium_electrons * physics_variables.dene,
                physics_variables.nd_alphas * (physics_variables.f_helium3 + 1.0e-3),
            )  # rough estimate
        else:
            physics_variables.nd_protons = max(
                physics_variables.f_nd_protium_electrons * physics_variables.dene,
                physics_variables.nd_alphas
                * physics_variables.proton_rate_density
                / physics_variables.fusden_alpha_total,
            )

        # ======================================================================

        # Beam hot ion component
        # If ignited, prevent beam fusion effects
        if physics_variables.i_plasma_ignited == 0:
            physics_variables.nd_beam_ions = (
                physics_variables.dene * physics_variables.f_nd_beam_electron
            )
        else:
            physics_variables.nd_beam_ions = 0.0

        # ======================================================================

        # Sum of Zi.ni for all impurity ions (those with charge > helium)
        znimp = 0.0
        for imp in range(impurity_radiation_module.n_impurities):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                znimp += impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.te])
                ).squeeze() * (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * physics_variables.dene
                )

        # ======================================================================

        # Fuel portion - conserve charge neutrality
        # znfuel is the sum of Zi.ni for the three fuel ions
        znfuel = (
            physics_variables.dene
            - 2.0 * physics_variables.nd_alphas
            - physics_variables.nd_protons
            - physics_variables.nd_beam_ions
            - znimp
        )

        # ======================================================================

        # Fuel ion density, nd_fuel_ions
        # For D-T-He3 mix, nd_fuel_ions = nD + nT + nHe3, while znfuel = nD + nT + 2*nHe3
        # So nd_fuel_ions = znfuel - nHe3 = znfuel - f_helium3*nd_fuel_ions
        physics_variables.nd_fuel_ions = znfuel / (1.0 + physics_variables.f_helium3)

        # ======================================================================

        # Set hydrogen and helium impurity fractions for
        # radiation calculations
        impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("H_")
        ] = (
            physics_variables.nd_protons
            + (physics_variables.f_deuterium + physics_variables.f_tritium)
            * physics_variables.nd_fuel_ions
            + physics_variables.nd_beam_ions
        ) / physics_variables.dene

        impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("He")
        ] = (
            physics_variables.f_helium3
            * physics_variables.nd_fuel_ions
            / physics_variables.dene
            + physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Total impurity density
        physics_variables.nd_impurities = 0.0
        for imp in range(impurity_radiation_module.n_impurities):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.nd_impurities += (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * physics_variables.dene
                )

        # ======================================================================

        # Total ion density
        physics_variables.nd_ions_total = (
            physics_variables.nd_fuel_ions
            + physics_variables.nd_alphas
            + physics_variables.nd_protons
            + physics_variables.nd_beam_ions
            + physics_variables.nd_impurities
        )

        # ======================================================================

        # Set some impurity fraction variables
        # for the benefit of other routines
        physics_variables.rncne = impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("C_")
        ]
        physics_variables.rnone = impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("O_")
        ]
        physics_variables.rnfene = (
            impurity_radiation_module.impurity_arr_frac[
                impurity_radiation.element2index("Fe")
            ]
            + impurity_radiation_module.impurity_arr_frac[
                impurity_radiation.element2index("Ar")
            ]
        )

        # ======================================================================

        # Effective charge
        # Calculation should be sum(ni.Zi^2) / sum(ni.Zi),
        # but ne = sum(ni.Zi) through quasineutrality
        physics_variables.zeff = 0.0
        for imp in range(impurity_radiation_module.n_impurities):
            physics_variables.zeff += (
                impurity_radiation_module.impurity_arr_frac[imp]
                * impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.te])
                ).squeeze()
                ** 2
            )

        # ======================================================================

        # Fraction of alpha energy to ions and electrons
        # From Max Fenstermacher
        # (used with electron and ion power balance equations only)
        # No consideration of pden_non_alpha_charged_mw here...

        # pcoef now calculated in plasma_profiles, after the very first
        # call of plasma_composition; use old parabolic profile estimate
        # in this case
        if physics_module.first_call == 1:
            pc = (
                (1.0 + physics_variables.alphan)
                * (1.0 + physics_variables.alphat)
                / (1.0 + physics_variables.alphan + physics_variables.alphat)
            )
            physics_module.first_call = 0
        else:
            pc = physics_variables.pcoef

        physics_variables.f_alpha_electron = 0.88155 * np.exp(
            -physics_variables.te * pc / 67.4036
        )
        physics_variables.f_alpha_ion = 1.0 - physics_variables.f_alpha_electron

        # ======================================================================

        # Average atomic masses of injected fuel species
        physics_variables.m_fuel_amu = (
            (constants.m_deuteron_amu * physics_variables.f_deuterium)
            + (constants.m_triton_amu * physics_variables.f_tritium)
            + (constants.m_helion_amu * physics_variables.f_helium3)
        )

        # ======================================================================

        # Average atomic masses of injected fuel species in the neutral beams
        # Only deuterium and tritium in the beams
        physics_variables.m_beam_amu = (
            constants.m_deuteron_amu * (1.0 - current_drive_variables.f_beam_tritium)
        ) + (constants.m_triton_amu * current_drive_variables.f_beam_tritium)

        # ======================================================================

        # Average mass of all ions
        physics_variables.m_ions_total_amu = (
            (physics_variables.m_fuel_amu * physics_variables.nd_fuel_ions)
            + (constants.m_alpha_amu * physics_variables.nd_alphas)
            + (physics_variables.nd_protons * constants.m_proton_amu)
            + (physics_variables.m_beam_amu * physics_variables.nd_beam_ions)
        )
        for imp in range(impurity_radiation_module.n_impurities):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.m_ions_total_amu += (
                    physics_variables.dene
                    * impurity_radiation_module.impurity_arr_frac[imp]
                    * impurity_radiation_module.impurity_arr_amass[imp]
                )

        physics_variables.m_ions_total_amu = (
            physics_variables.m_ions_total_amu / physics_variables.nd_ions_total
        )

        # ======================================================================

        # Mass weighted plasma effective charge
        # Sum of (Zi^2*n_i) / m_i
        physics_variables.zeffai = (
            (
                physics_variables.f_deuterium
                * physics_variables.nd_fuel_ions
                / constants.m_deuteron_amu
            )
            + (
                physics_variables.f_tritium
                * physics_variables.nd_fuel_ions
                / constants.m_triton_amu
            )
            + (
                4.0
                * physics_variables.f_helium3
                * physics_variables.nd_fuel_ions
                / constants.m_helion_amu
            )
            + (4.0 * physics_variables.nd_alphas / constants.m_alpha_amu)
            + (physics_variables.nd_protons / constants.m_proton_amu)
            + (
                (1.0 - current_drive_variables.f_beam_tritium)
                * physics_variables.nd_beam_ions
                / constants.m_deuteron_amu
            )
            + (
                current_drive_variables.f_beam_tritium
                * physics_variables.nd_beam_ions
                / constants.m_triton_amu
            )
        ) / physics_variables.dene
        for imp in range(impurity_radiation_module.n_impurities):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.zeffai += (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * impurity_radiation.zav_of_te(
                        imp, np.array([physics_variables.te])
                    ).squeeze()
                    ** 2
                    / impurity_radiation_module.impurity_arr_amass[imp]
                )

        # ======================================================================

    @staticmethod
    def phyaux(
        aspect: float,
        dene: float,
        te: float,
        nd_fuel_ions: float,
        fusden_total: float,
        fusden_alpha_total: float,
        plasma_current: float,
        sbar: float,
        nd_alphas: float,
        t_energy_confinement: float,
        vol_plasma: float,
    ) -> tuple[float, float, float, float, float, float, float, float]:
        """Auxiliary physics quantities

        Args:
            aspect (float): Plasma aspect ratio.
            dene (float): Electron density (/m3).
            te (float): Volume avergaed electron temperature (keV).
            nd_fuel_ions (float): Fuel ion density (/m3).
            fusden_total (float): Fusion reaction rate from plasma and beams (/m3/s).
            fusden_alpha_total (float): Alpha particle production rate (/m3/s).
            plasma_current (float): Plasma current (A).
            sbar (float): Exponent for aspect ratio (normally 1).
            nd_alphas (float): Alpha ash density (/m3).
            t_energy_confinement (float): Global energy confinement time (s).
            vol_plasma (float): Plasma volume (m3).

        Returns:
            tuple: A tuple containing:
                - burnup (float): Fractional plasma burnup.
                - ntau (float): Plasma average n-tau (s/m3).
                - nTtau (float): Plasma triple product nT-tau (s/m3).
                - figmer (float): Physics figure of merit.
                - fusrat (float): Number of fusion reactions per second.
                - qfuel (float): Fuelling rate for D-T (nucleus-pairs/sec).
                - rndfuel (float): Fuel burnup rate (reactions/s).
                - t_alpha_confinement (float): Alpha particle confinement time (s).
                - f_alpha_energy_confinement (float): Fraction of alpha energy confinement.

        This subroutine calculates extra physics related items needed by other parts of the code.
        """
        figmer = 1e-6 * plasma_current * aspect**sbar

        ntau = t_energy_confinement * dene
        nTtau = ntau * te

        # Fusion reactions per second
        fusrat = fusden_total * vol_plasma

        # Alpha particle confinement time (s)
        # Number of alphas / alpha production rate
        if fusden_alpha_total != 0.0:
            t_alpha_confinement = nd_alphas / fusden_alpha_total
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
                nd_alphas
                / (nd_alphas + 0.5 * nd_fuel_ions)
                / physics_variables.tauratio
            )
        else:
            burnup = physics_variables.burnup_in

        # Fuel burnup rate (reactions/second) (previously Amps)
        rndfuel = fusrat

        # Required fuelling rate (fuel ion pairs/second) (previously Amps)
        qfuel = rndfuel / burnup

        f_alpha_energy_confinement = t_alpha_confinement / t_energy_confinement

        return (
            burnup,
            ntau,
            nTtau,
            figmer,
            fusrat,
            qfuel,
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
        ten: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float]:
        """
        Calculate the ohmic heating power and related parameters.

        Args:
            f_c_plasma_inductive (float): Fraction of plasma current driven inductively.
            kappa95 (float): Plasma elongation at 95% surface.
            plasma_current (float): Plasma current (A).
            rmajor (float): Major radius (m).
            rminor (float): Minor radius (m).
            ten (float): Density weighted average electron temperature (keV).
            vol_plasma (float): Plasma volume (m^3).
            zeff (float): Plasma effective charge.

        Returns:
            Tuple[float, float, float, float]: Tuple containing:
                - pden_plasma_ohmic_mw (float): Ohmic heating power per unit volume (MW/m^3).
                - p_plasma_ohmic_mw (float): Total ohmic heating power (MW).
                - f_res_plasma_neo (float): Neo-classical resistivity enhancement factor.
                - res_plasma (float): Plasma resistance (ohm).

        Notes:

        References:
            - ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,

        """
        # Density weighted electron temperature in 10 keV units
        t10 = ten / 10.0

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
            error_handling.fdiags[0] = res_plasma
            error_handling.fdiags[1] = physics_variables.aspect
            error_handling.report_error(83)

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
        bt: float,
        eps: float,
        i_plasma_current: int,
        kappa: float,
        kappa95: float,
        p0: float,
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
            bt (float): Toroidal field on axis (T).
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
            p0 (float): Central plasma pressure (Pa).
            len_plasma_poloidal (float): Plasma perimeter length (m).
            q95 (float): Plasma safety factor at 95% flux (= q-bar for i_plasma_current=2).
            ind_plasma_internal_norm (float): Plasma normalised internal inductance.
            rmajor (float): Major radius (m).
            rminor (float): Minor radius (m).
            triang (float): Plasma triangularity.
            triang95 (float): Plasma triangularity at 95% surface.

        Returns:
            Tuple[float, float, float,]: Tuple containing bp, qstar, plasma_current,

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
            fq = calculate_current_coefficient_peng(eps, len_plasma_poloidal, rminor)

        # Peng scaling for double null divertor; TARTs [STAR Code]
        elif i_plasma_current == 2:
            plasma_current = 1.0e6 * calculate_plasma_current_peng(
                q95, aspect_ratio, eps, rminor, bt, kappa, triang
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
                alphaj, alphap, bt, triang95, eps, kappa95, p0, constants.rmu0
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
                (constants.twopi / constants.rmu0)
                * rminor**2
                / (rmajor * q95)
                * fq
                * bt
            )
        # i_plasma_current == 2 case covered above

        # Calculate cyclindrical safety factor from IPDG89
        qstar = (
            5.0e6
            * rminor**2
            / (rmajor * plasma_current / bt)
            * 0.5
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

        # Normalised beta from Troyon beta limit
        physics_variables.beta_norm_total = (
            1.0e8 * physics_variables.beta * rminor * bt / plasma_current
        )

        # Calculate the poloidal field generated by the plasma current
        bp = calculate_poloidal_field(
            i_plasma_current,
            plasma_current,
            q95,
            aspect_ratio,
            eps,
            bt,
            kappa,
            triang,
            len_plasma_poloidal,
            constants.rmu0,
        )

        return bp, qstar, plasma_current

    def outtim(self):
        po.oheadr(self.outfile, "Times")

        po.ovarrf(
            self.outfile,
            "Initial charge time for CS from zero current (s)",
            "(t_precharge)",
            times_variables.t_precharge,
        )
        po.ovarrf(
            self.outfile,
            "Plasma current ramp-up time (s)",
            "(t_current_ramp_up)",
            times_variables.t_current_ramp_up,
        )
        po.ovarrf(
            self.outfile,
            "Heating time (s)",
            "(t_fusion_ramp)",
            times_variables.t_fusion_ramp,
        )
        po.ovarre(
            self.outfile, "Burn time (s)", "(t_burn)", times_variables.t_burn, "OP "
        )
        po.ovarrf(
            self.outfile,
            "Reset time to zero current for CS (s)",
            "(t_ramp_down)",
            times_variables.t_ramp_down,
        )
        po.ovarrf(
            self.outfile,
            "Time between pulses (s)",
            "(t_between_pulse)",
            times_variables.t_between_pulse,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plant cycle time (s)",
            "(t_cycle)",
            times_variables.t_cycle,
            "OP ",
        )

    def outplas(self):
        """Subroutine to output the plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        self.outfile : input integer : Fortran output unit identifier
        This routine writes the plasma physics information
        to a file, in a tidy format.
        """

        # ###############################################
        # Dimensionless plasma parameters. See reference below.
        physics_module.nu_star = (
            1
            / constants.rmu0
            * (15.0e0 * constants.electron_charge**4 * physics_variables.dlamie)
            / (4.0e0 * np.pi**1.5e0 * constants.epsilon0**2)
            * physics_variables.vol_plasma**2
            * physics_variables.rmajor**2
            * physics_variables.bt
            * np.sqrt(physics_variables.eps)
            * physics_variables.nd_electron_line**3
            * physics_variables.kappa
            / (physics_module.e_plasma_beta**2 * physics_variables.plasma_current)
        )

        physics_module.rho_star = np.sqrt(
            2.0e0
            * constants.proton_mass
            * physics_variables.m_ions_total_amu
            * physics_module.e_plasma_beta
            / (
                3.0e0
                * physics_variables.vol_plasma
                * physics_variables.nd_electron_line
            )
        ) / (
            constants.electron_charge
            * physics_variables.bt
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_module.beta_mcdonald = (
            4.0e0
            / 3.0e0
            * constants.rmu0
            * physics_module.e_plasma_beta
            / (physics_variables.vol_plasma * physics_variables.bt**2)
        )

        po.oheadr(self.outfile, "Plasma")

        if stellarator_variables.istell == 0:
            if physics_variables.n_divertors == 0:
                po.ocmmnt(self.outfile, "Plasma configuration = limiter")
            elif physics_variables.n_divertors == 1:
                po.ocmmnt(self.outfile, "Plasma configuration = single null divertor")
            elif physics_variables.n_divertors == 2:
                po.ocmmnt(self.outfile, "Plasma configuration = double null divertor")
            else:
                raise ProcessValueError(
                    "Illegal value of n_divertors",
                    n_divertors=physics_variables.n_divertors,
                )
        else:
            po.ocmmnt(self.outfile, "Plasma configuration = stellarator")

        if stellarator_variables.istell == 0:
            if physics_variables.itart == 0:
                physics_module.itart_r = physics_variables.itart
                po.ovarin(
                    self.outfile,
                    "Tokamak aspect ratio = Conventional, itart = 0",
                    "(itart)",
                    physics_module.itart_r,
                )
            elif physics_variables.itart == 1:
                physics_module.itart_r = physics_variables.itart
                po.ovarin(
                    self.outfile,
                    "Tokamak aspect ratio = Spherical, itart = 1",
                    "(itart)",
                    physics_module.itart_r,
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
        po.ovarrf(
            self.outfile, "Major radius (m)", "(rmajor)", physics_variables.rmajor
        )
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
                    "(j_plasma_0)",
                    physics_variables.j_plasma_0,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Vertical field at plasma (T)",
                    "(bvert)",
                    physics_variables.bvert,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Vacuum toroidal field at R (T)",
                "(bt)",
                physics_variables.bt,
            )
            po.ovarrf(
                self.outfile,
                "Average poloidal field (T)",
                "(bp)",
                physics_variables.bp,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Total field (sqrt(bp^2 + bt^2)) (T)",
                "(btot)",
                physics_variables.btot,
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
        po.osubhd(self.outfile, "Beta Information :")
        if physics_variables.i_beta_component == 0:
            po.ovarrf(
                self.outfile,
                "Upper limit on total beta",
                "(beta_max)",
                physics_variables.beta_max,
                "OP ",
            )
        elif physics_variables.i_beta_component == 1:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal beta",
                "(beta_max)",
                physics_variables.beta_max,
                "OP ",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal + NB beta",
                "(beta_max)",
                physics_variables.beta_max,
                "OP ",
            )

        po.ovarre(self.outfile, "Total plasma beta", "(beta)", physics_variables.beta)
        if physics_variables.i_beta_component == 0:
            po.ovarrf(
                self.outfile,
                "Lower limit on total beta",
                "(beta_min)",
                physics_variables.beta_min,
                "IP",
            )
        elif physics_variables.i_beta_component == 1:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal beta",
                "(beta_min)",
                physics_variables.beta_min,
                "IP",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal + NB beta",
                "(beta_min)",
                physics_variables.beta_min,
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
            "(beta_poloidal)",
            physics_variables.beta_poloidal,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total toroidal beta",
            "(beta_toroidal)",
            physics_variables.beta_toroidal,
            "OP ",
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
            "Thermal beta",
            "(beta_thermal)",
            physics_variables.beta_thermal,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal poloidal beta",
            "(beta_thermal_poloidal)",
            physics_variables.beta_thermal_poloidal,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal toroidal beta",
            "(beta_thermal_toroidal)",
            physics_variables.beta_thermal_toroidal,
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
            physics_module.e_plasma_beta,
            "OP",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.osubhd(self.outfile, "Temperature and Density (volume averaged) :")
        po.ovarrf(
            self.outfile,
            "Volume averaged electron temperature (keV)",
            "(te)",
            physics_variables.te,
        )
        po.ovarrf(
            self.outfile,
            "Ratio of ion to electron volume-averaged temperature",
            "(tratio)",
            physics_variables.tratio,
            "IP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temperature on axis (keV)",
            "(te0)",
            physics_variables.te0,
            "OP ",
        )
        po.ovarrf(self.outfile, "Ion temperature (keV)", "(ti)", physics_variables.ti)
        po.ovarrf(
            self.outfile,
            "Ion temperature on axis (keV)",
            "(ti0)",
            physics_variables.ti0,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temp., density weighted (keV)",
            "(ten)",
            physics_variables.ten,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged electron number density (/m3)",
            "(dene)",
            physics_variables.dene,
        )
        po.ovarre(
            self.outfile,
            "Electron number density on axis (/m3)",
            "(ne0)",
            physics_variables.ne0,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line-averaged electron number density (/m3)",
            "(nd_electron_line)",
            physics_variables.nd_electron_line,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma pressure on axis (Pa)",
            "(p0)",
            physics_variables.p0,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged plasma pressure (Pa)",
            "(vol_avg_pressure)",
            physics_variables.vol_avg_pressure,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.ovarre(
                self.outfile,
                "Line-averaged electron density / Greenwald density",
                "(dnla_gw)",
                physics_variables.nd_electron_line / physics_variables.dlimit[6],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total Ion number density (/m3)",
            "(nd_ions_total)",
            physics_variables.nd_ions_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel ion number density (/m3)",
            "(nd_fuel_ions)",
            physics_variables.nd_fuel_ions,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total impurity number density with Z > 2 (no He) (/m3)",
            "(nd_impurities)",
            physics_variables.nd_impurities,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion number density (thermalised ions only) (/m3)",
            "(nd_alphas)",
            physics_variables.nd_alphas,
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
            "Proton number density (/m3)",
            "(nd_protons)",
            physics_variables.nd_protons,
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

        for imp in range(impurity_radiation_module.n_impurities):
            # MDK Update fimp, as this will make the ITV output work correctly.
            impurity_radiation_module.fimp[imp] = (
                impurity_radiation_module.impurity_arr_frac[imp]
            )
            str1 = (
                f2py_compatible_to_string(
                    impurity_radiation_module.impurity_arr_label[imp]
                )
                + " concentration"
            )
            str2 = f"(fimp({imp + 1:02}))"
            # MDK Add output flag for H which is calculated.
            if imp == 0:
                po.ovarre(
                    self.outfile, str1, str2, impurity_radiation_module.fimp[imp], "OP "
                )
            else:
                po.ovarre(self.outfile, str1, str2, impurity_radiation_module.fimp[imp])

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
            self.outfile, "Effective charge", "(zeff)", physics_variables.zeff, "OP "
        )

        po.ovarrf(
            self.outfile,
            "Mass-weighted Effective charge",
            "(zeffai)",
            physics_variables.zeffai,
            "OP ",
        )

        po.ovarrf(
            self.outfile, "Density profile factor", "(alphan)", physics_variables.alphan
        )
        po.ovarin(
            self.outfile,
            "Plasma profile model",
            "(ipedestal)",
            physics_variables.ipedestal,
        )

        if physics_variables.ipedestal >= 1:
            if physics_variables.ne0 < physics_variables.neped:
                error_handling.report_error(213)

            po.ocmmnt(self.outfile, "Pedestal profiles are used.")
            po.ovarrf(
                self.outfile,
                "Density pedestal r/a location",
                "(rhopedn)",
                physics_variables.rhopedn,
            )
            if physics_variables.fgwped >= 0e0:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                )

            # This code is ODD# Don't change it# No explanation why fgwped and physics_variables.fgwsep
            # must be assigned to their exisiting values#
            fgwped_out = physics_variables.neped / physics_variables.dlimit[6]
            fgwsep_out = physics_variables.nesep / physics_variables.dlimit[6]
            if physics_variables.fgwped >= 0e0:
                physics_variables.fgwped = (
                    physics_variables.neped / physics_variables.dlimit[6]
                )
            if physics_variables.fgwsep >= 0e0:
                physics_variables.fgwsep = (
                    physics_variables.nesep / physics_variables.dlimit[6]
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
                "(rhopedt)",
                physics_variables.rhopedt,
            )

            po.ovarrf(
                self.outfile,
                "Electron temp. pedestal height (keV)",
                "(teped)",
                physics_variables.teped,
            )
            if any(numerics.icc == 78):
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                )

            po.ovarre(
                self.outfile,
                "Electron density at separatrix (/m3)",
                "(nesep)",
                physics_variables.nesep,
            )
            po.ovarre(
                self.outfile,
                "Electron density at separatrix / nGW",
                "(fgwsep_out)",
                fgwsep_out,
            )

        # Issue 558 - addition of constraint 76 to limit the value of nesep, in proportion with the ballooning parameter and Greenwald density
        if any(numerics.icc == 76):
            po.ovarre(
                self.outfile,
                "Critical ballooning parameter value",
                "(alpha_crit)",
                physics_variables.alpha_crit,
            )
            po.ovarre(
                self.outfile,
                "Critical electron density at separatrix (/m3)",
                "(nesep_crit)",
                physics_variables.nesep_crit,
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
                "(dlimit(1))",
                physics_variables.dlimit[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model I",
                "(dlimit(2))",
                physics_variables.dlimit[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model II",
                "(dlimit(3))",
                physics_variables.dlimit[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET edge radiation model",
                "(dlimit(4))",
                physics_variables.dlimit[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET simplified model",
                "(dlimit(5))",
                physics_variables.dlimit[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hugill-Murakami Mq model",
                "(dlimit(6))",
                physics_variables.dlimit[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Greenwald model",
                "(dlimit(7))",
                physics_variables.dlimit[6],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ASDEX New",
                "(dlimit(8))",
                physics_variables.dlimit[7],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Density limit from scaling (/m3)",
                "(dnelimt)",
                physics_variables.dnelimt,
                "OP ",
            )
            if (numerics.ioptimz > 0) and (numerics.active_constraints[4]):
                po.ovarre(
                    self.outfile,
                    "Density limit (enforced) (/m3)",
                    "(boundu(9)*dnelimt)",
                    numerics.boundu[8] * physics_variables.dnelimt,
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
            "(f_deuterium)",
            physics_variables.f_deuterium,
        )
        po.ovarrf(
            self.outfile,
            "Tritium fuel fraction",
            "(f_tritium)",
            physics_variables.f_tritium,
        )
        po.ovarrf(
            self.outfile,
            "3-Helium fuel fraction",
            "(f_helium3)",
            physics_variables.f_helium3,
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
            "Fraction of core radiation subtracted from P_L",
            "(coreradiationfraction)",
            impurity_radiation_module.coreradiationfraction,
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
            "LCFS radiation fraction = total radiation in LCFS / total power deposited in plasma",
            "(rad_fraction_LCFS)",
            physics_module.rad_fraction_lcfs,
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
                physics_module.ptarmw,
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
                physics_module.lambdaio,
                "OP ",
            )
            if physics_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Midplane seperation of the two magnetic closed flux surfaces (m)",
                    "(drsep)",
                    physics_module.drsep,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Fraction of power on the inner targets",
                "(fio)",
                physics_module.fio,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower inner target",
                "(fLI)",
                physics_module.fli,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower outer target",
                "(fLO)",
                physics_module.flo,
                "OP ",
            )
            if physics_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper inner target",
                    "(fUI)",
                    physics_module.fui,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper outer target",
                    "(fUO)",
                    physics_module.fuo,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Power incident on the lower inner target (MW)",
                "(pLImw)",
                physics_module.plimw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Power incident on the lower outer target (MW)",
                "(pLOmw)",
                physics_module.plomw,
                "OP ",
            )
            if physics_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper innner target (MW)",
                    "(pUImw)",
                    physics_module.puimw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper outer target (MW)",
                    "(pUOmw)",
                    physics_module.puomw,
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
            error_handling.fdiags[0] = physics_variables.p_plasma_separatrix_mw
            error_handling.report_error(87)
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(
                self.outfile, "          Power into divertor zone is unrealistic;"
            )
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if physics_variables.n_divertors == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Pdivt / R ratio (MW/m) (On peak divertor)",
                "(pdivmax/physics_variables.rmajor)",
                physics_variables.pdivmax / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)",
                "(pdivmaxbt/qar)",
                (
                    (physics_variables.pdivmax * physics_variables.bt)
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
                    (physics_variables.p_plasma_separatrix_mw * physics_variables.bt)
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
                if (physics_variables.bt < 0.78e0) or (physics_variables.bt > 7.94e0):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.bt outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(201)

                if (physics_variables.rminor < 0.15e0) or (
                    physics_variables.rminor > 1.15e0
                ):
                    po.ocmmnt(self.outfile, "(rminor outside Snipes 2000 fitted range)")
                    error_handling.report_error(202)

                if (physics_variables.rmajor < 0.55e0) or (
                    physics_variables.rmajor > 3.37e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.rmajor outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(203)

                if (physics_variables.nd_electron_line < 0.09e20) or (
                    physics_variables.nd_electron_line > 3.16e20
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.nd_electron_line outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(204)

                if (physics_variables.kappa < 1.0e0) or (
                    physics_variables.kappa > 2.04e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.kappa outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(205)

                if (physics_variables.triang < 0.07e0) or (
                    physics_variables.triang > 0.74e0
                ):
                    po.ocmmnt(self.outfile, "(triang outside Snipes 2000 fitted range)")
                    error_handling.report_error(206)

            po.oblnkl(self.outfile)

            if physics_variables.i_l_h_threshold in [12, 13, 14]:
                po.ocmmnt(
                    self.outfile,
                    "(L-H threshold for closed divertor only. Limited data used in Snipes fit)",
                )
                po.oblnkl(self.outfile)
                error_handling.report_error(207)

            if (numerics.ioptimz > 0) and (numerics.active_constraints[14]):
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (enforced) (MW)",
                    "(boundl(103)*p_l_h_threshold_mw)",
                    numerics.boundl[102] * physics_variables.p_l_h_threshold_mw,
                    "OP ",
                )
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

        tauelaw = f2py_compatible_to_string(
            physics_variables.labels_confinement_scalings[
                physics_variables.i_confinement_time
            ]
        )

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
            physics_module.t_energy_confinement_beta,
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
        if physics_variables.i_rad_loss == 1:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact
                * (
                    physics_variables.p_plasma_loss_mw
                    / (
                        physics_variables.p_plasma_loss_mw
                        + physics_variables.pden_plasma_sync_mw
                        + physics_variables.p_plasma_inner_rad_mw
                    )
                )
                ** 0.31,
                "OP",
            )
        elif physics_variables.i_rad_loss == 0:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact
                * (
                    physics_variables.p_plasma_loss_mw
                    / (
                        physics_variables.p_plasma_loss_mw
                        + physics_variables.pden_plasma_rad_mw
                        * physics_variables.vol_plasma
                    )
                )
                ** 0.31,
                "OP",
            )
        elif physics_variables.i_rad_loss == 2:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact,
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
                physics_module.beta_mcdonald,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized ion Larmor radius",
                "(rho_star)",
                physics_module.rho_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized collisionality",
                "(nu_star)",
                physics_module.nu_star,
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
                "(bscf_gi_I)",
                current_drive_variables.bscf_gi_I,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Gi II)",
                "(bscf_gi_II)",
                current_drive_variables.bscf_gi_II,
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
            if physics_module.err242 == 1:
                error_handling.report_error(242)

            # Error to catch if self-driven current fraction limit has been enforced
            if physics_module.err243 == 1:
                error_handling.report_error(243)

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
                    error_handling.report_error(244)

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
            "(qfuel)",
            physics_variables.qfuel,
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

        if any(numerics.icc == 78):
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

    def output_confinement_comparison(self, istell: int) -> None:
        """
        Routine to calculate ignition margin for different confinement scalings and equivalent confinement times for H=1.

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

        Parameters:
        - istell (int): Indicator for stellarator (0 for tokamak, >=1 for stellarator).

        Returns:
        - None
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
            range(1, physics_variables.n_confinement_scalings)
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
                physics_variables.bt,
                physics_variables.nd_ions_total,
                physics_variables.dene,
                physics_variables.nd_electron_line,
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
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.q95,
                physics_variables.qstar,
                physics_variables.vol_plasma,
                physics_variables.zeff,
            )

            # Calculate the H-factor for the same confinement time in other scalings
            physics_variables.hfac[i_confinement_time - 1] = self.find_other_h_factors(
                i_confinement_time
            )

            po.ocmmnt(
                self.outfile,
                f"{'':>2}{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[i_confinement_time]):<38}"
                f"{taueez:<28.3f}{physics_variables.hfac[i_confinement_time - 1]:.3f}",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)

    @staticmethod
    def bootstrap_fraction_iter89(
        aspect: float,
        beta: float,
        bt: float,
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
            bt (float): Toroidal field on axis (T).
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
        betapbs = beta * (bt / b_pa) ** 2

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
    def bootstrap_fraction_nevins(
        alphan: float,
        alphat: float,
        beta_toroidal: float,
        bt: float,
        dene: float,
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
            bt (float): Toroidal field on axis (T).
            dene (float): Electron density (/m3).
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
            physics_variables.ne0
            * physics_variables.te0
            * 1.0e3
            * constants.electron_charge
            / (bt**2 / (2.0 * constants.rmu0))
        )

        # Call integration routine using definite integral routine from scipy

        ainteg, _ = integrate.quad(
            lambda y: _nevins_integral(
                y,
                dene,
                te,
                bt,
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

        aibs = 2.5 * betae0 * rmajor * bt * q95 * ainteg
        return 1.0e6 * aibs / plasma_current

    @staticmethod
    def bootstrap_fraction_sauter(plasma_profile: float) -> float:
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

        # Number of radial data points along the profile
        nr = plasma_profile.profile_size

        # Radial points from 0 to 1 seperated by 1/profile_size
        roa = plasma_profile.neprofile.profile_x

        # Local circularised minor radius
        rho = np.sqrt(physics_variables.a_plasma_poloidal / np.pi) * roa

        # Square root of local aspect ratio
        sqeps = np.sqrt(roa * (physics_variables.rminor / physics_variables.rmajor))

        # Calculate electron and ion density profiles
        ne = plasma_profile.neprofile.profile_y * 1e-19
        ni = (physics_variables.nd_ions_total / physics_variables.dene) * ne

        # Calculate electron and ion temperature profiles
        tempe = plasma_profile.teprofile.profile_y
        tempi = (physics_variables.ti / physics_variables.te) * tempe

        # Flat Zeff profile assumed
        # Return tempi like array object filled with zeff
        zeff = np.full_like(tempi, physics_variables.zeff)

        # inverse_q = 1/safety factor
        # Parabolic q profile assumed
        inverse_q = 1 / (
            physics_variables.q0
            + (physics_variables.q95 - physics_variables.q0) * roa**2
        )
        # Create new array of average mass of fuel portion of ions
        amain = np.full_like(inverse_q, physics_variables.m_fuel_amu)

        # Create new array of average main ion charge
        zmain = np.full_like(inverse_q, 1.0 + physics_variables.f_helium3)

        # Prevent division by zero
        if ne[nr - 1] == 0.0:
            ne[nr - 1] = 1e-4 * ne[nr - 2]
            ni[nr - 1] = 1e-4 * ni[nr - 2]

        # Prevent division by zero
        if tempe[nr - 1] == 0.0:
            tempe[nr - 1] = 1e-4 * tempe[nr - 2]
            tempi[nr - 1] = 1e-4 * tempi[nr - 2]

        # Calculate total bootstrap current (MA) by summing along profiles
        # Looping from 2 because _calculate_l31_coefficient() etc should return 0 @ j == 1
        radial_elements = np.arange(2, nr)

        # Change in localised minor radius to be used as delta term in derivative
        drho = rho[radial_elements] - rho[radial_elements - 1]

        # Area of annulus, assuming circular plasma cross-section
        da = 2 * np.pi * rho[radial_elements - 1] * drho  # area of annulus

        # Create the partial derivatives
        dlogte_drho = (
            np.log(tempe[radial_elements]) - np.log(tempe[radial_elements - 1])
        ) / drho
        dlogti_drho = (
            np.log(tempi[radial_elements]) - np.log(tempi[radial_elements - 1])
        ) / drho
        dlogne_drho = (
            np.log(ne[radial_elements]) - np.log(ne[radial_elements - 1])
        ) / drho

        jboot = (
            0.5
            * (
                _calculate_l31_coefficient(
                    radial_elements,
                    nr,
                    physics_variables.rmajor,
                    physics_variables.bt,
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
                    nr,
                    physics_variables.rmajor,
                    physics_variables.bt,
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
                    nr,
                    physics_variables.rmajor,
                    physics_variables.bt,
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
                -physics_variables.bt
                / (0.2 * np.pi * physics_variables.rmajor)
                * rho[radial_elements - 1]
                * inverse_q[radial_elements - 1]
            )
        )  # A/m2

        return np.sum(da * jboot, axis=0) / physics_variables.plasma_current

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
        rhopedn: float,
        neped: float,
        n_greenwald: float,
        teped: float,
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
        :param rhopedn: Normalised plasma radius of density pedestal.
        :type rhopedn: float
        :param neped: Electron number density at the pedestal [m^-3].
        :type neped: float
        :param n_greenwald: Greenwald density limit [m^-3].
        :type n_greenwald: float
        :param teped: Electron temperature at the pedestal [keV].
        :type teped: float

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
            * rhopedn**0.367
            * (neped / n_greenwald) ** -0.174
            * teped**0.0552
        )

    def find_other_h_factors(self, i_confinement_time: int) -> float:
        """
        Function to find H-factor for the equivalent confinement time in other scalings.

        Args:
            i_confinement_time (int): Index of the confinement time scaling to use.

        Returns:
            float: The calculated H-factor.
        """

        def fhz(hfact: float) -> float:
            """
            Function used to find power balance.

            Args:
                hfact (float): H-factor to be used in the calculation.

            Returns:
                float: The difference between the calculated power and the required power for balance.
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
                physics_variables.bt,
                physics_variables.nd_ions_total,
                physics_variables.dene,
                physics_variables.nd_electron_line,
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
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.q95,
                physics_variables.qstar,
                physics_variables.vol_plasma,
                physics_variables.zeff,
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
        bt: float,
        nd_ions_total: float,
        dene: float,
        nd_electron_line: float,
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
        ten: float,
        tin: float,
        q95: float,
        qstar: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float, float, float, float]:
        """
        Calculate the confinement times and the transport power loss terms.

        :param m_fuel_amu: Average mass of fuel (amu)
        :param p_alpha_total_mw: Alpha particle power (MW)
        :param aspect: Aspect ratio
        :param bt: Toroidal field on axis (T)
        :param nd_ions_total: Total ion density (/m3)
        :param dene: Volume averaged electron density (/m3)
        :param nd_electron_line: Line-averaged electron density (/m3)
        :param eps: Inverse aspect ratio
        :param hfact: H factor on energy confinement scalings
        :param i_confinement_time: Switch for energy confinement scaling to use
        :param i_plasma_ignited: Switch for ignited calculation
        :param kappa: Plasma elongation
        :param kappa95: Plasma elongation at 95% surface
        :param p_non_alpha_charged_mw: Non-alpha charged particle fusion power (MW)
        :param p_hcd_injected_total_mw: Auxiliary power to ions and electrons (MW)
        :param plasma_current: Plasma current (A)
        :param pden_plasma_core_rad_mw: Total core radiation power (MW/m3)
        :param q95: Edge safety factor (tokamaks), or rotational transform iotabar (stellarators)
        :param qstar: Equivalent cylindrical edge safety factor
        :param rmajor: Plasma major radius (m)
        :param rminor: Plasma minor radius (m)
        :param ten: Density weighted average electron temperature (keV)
        :param tin: Density weighted average ion temperature (keV)
        :param vol_plasma: Plasma volume (m3)
        :param a_plasma_poloidal: Plasma cross-sectional area (m2)
        :param zeff: Plasma effective charge

        :return: Tuple containing:
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
        dnla20 = nd_electron_line * 1.0e-20
        dnla19 = nd_electron_line * 1.0e-19

        # Volume averaged electron density in units of 10**20 m**-3
        n20 = dene / 1.0e20

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
                rmajor, rminor, kappa95, qstar, dnla20, m_fuel_amu, ten
            )

        # ========================================================================

        # Shimomura (S) optimized H-mode scaling
        elif i_confinement_time == 4:
            t_electron_confinement = confinement.shimomura_confinement_time(
                rmajor, rminor, bt, kappa95, m_fuel_amu
            )

        # ========================================================================

        # Kaye-Goldston scaling (L-mode)
        elif i_confinement_time == 5:
            t_electron_confinement = confinement.kaye_goldston_confinement_time(
                pcur, rmajor, rminor, kappa, dnla20, bt, m_fuel_amu, p_plasma_loss_mw
            )

        # ========================================================================

        # ITER Power scaling - ITER 89-P (L-mode)
        elif i_confinement_time == 6:
            t_electron_confinement = confinement.iter_89p_confinement_time(
                pcur, rmajor, rminor, kappa, dnla20, bt, m_fuel_amu, p_plasma_loss_mw
            )

        # ========================================================================

        # ITER Offset linear scaling - ITER 89-O (L-mode)
        elif i_confinement_time == 7:
            t_electron_confinement = confinement.iter_89_0_confinement_time(
                pcur, rmajor, rminor, kappa, dnla20, bt, m_fuel_amu, p_plasma_loss_mw
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
                bt,
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
                dnla20, rmajor, qstar, bt, rminor, kappa95, p_plasma_loss_mw, zeff, pcur
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
                bt,
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
                bt,
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
                bt,
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
                    bt,
                    m_fuel_amu,
                    p_plasma_loss_mw,
                ),
                confinement.iter_89_0_confinement_time(
                    pcur,
                    rmajor,
                    rminor,
                    kappa,
                    dnla20,
                    bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Amended version of ITER H90-P law
        elif i_confinement_time == 20:
            t_electron_confinement = confinement.iter_h90_p_amended_confinement_time(
                pcur,
                bt,
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
                bt,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Gyro-reduced Bohm scaling
        elif i_confinement_time == 22:
            t_electron_confinement = confinement.gyro_reduced_bohm_confinement_time(
                bt,
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
                    bt,
                    p_plasma_loss_mw,
                    q95,
                )
            )

        # ==========================================================================

        # ITER_93 ELM-free H-mode scaling
        elif i_confinement_time == 24:
            t_electron_confinement = confinement.iter_93h_confinement_time(
                pcur,
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
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
                bt,
                p_plasma_loss_mw,
                iotabar,
            )

        # ==========================================================================

        # DS03 beta-independent H-mode scaling
        elif i_confinement_time == 39:
            t_electron_confinement = confinement.ds03_confinement_time(
                pcur,
                bt,
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
                bt,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Petty08, beta independent dimensionless scaling
        elif i_confinement_time == 41:
            t_electron_confinement = confinement.petty08_confinement_time(
                pcur,
                bt,
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
                bt,
                nd_electron_line,
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
                bt,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - lower
        elif i_confinement_time == 44:
            t_electron_confinement = confinement.hubbard_lower_confinement_time(
                pcur,
                bt,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - upper
        elif i_confinement_time == 45:
            t_electron_confinement = confinement.hubbard_upper_confinement_time(
                pcur,
                bt,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Menard NSTX, ELMy H-mode scaling
        elif i_confinement_time == 46:
            t_electron_confinement = confinement.menard_nstx_confinement_time(
                pcur,
                bt,
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
                    bt,
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
                bt,
                p_plasma_loss_mw,
                rmajor,
                dnla20,
            )

        # ==========================================================================

        # ITPA20 H-mode scaling
        elif i_confinement_time == 49:
            t_electron_confinement = confinement.itpa20_confinement_time(
                pcur,
                bt,
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
                bt,
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

        # Calculation of the transport power loss terms
        # Transport losses in Watts/m3 are 3/2 * n.e.T / tau , with T in eV
        # (here, tin and ten are in keV, and pden_electron_transport_loss_mw and pden_ion_transport_loss_mw are in MW/m3)

        # The transport losses is just the electron and ion thermal energies divided by the confinement time.
        pden_ion_transport_loss_mw = (
            (3 / 2)
            * (constants.electron_charge / 1e3)
            * nd_ions_total
            * tin
            / t_ion_energy_confinement
        )
        pden_electron_transport_loss_mw = (
            (3 / 2)
            * (constants.electron_charge / 1e3)
            * dene
            * ten
            / t_electron_energy_confinement
        )

        ratio = (nd_ions_total / dene) * (tin / ten)

        # Global energy confinement time

        t_energy_confinement = (ratio + 1.0e0) / (
            ratio / t_ion_energy_confinement + 1.0e0 / t_electron_energy_confinement
        )

        # For comparison directly calculate the confinement time from the stored energy calculated
        # from the total plasma beta and the loss power used above.
        physics_module.t_energy_confinement_beta = (
            physics_module.e_plasma_beta / 1e6
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
        nd_ions_total: float,
        nd_fuel_ions: float,
        nd_alphas: float,
        vol_plasma: float,
        dene: float,
    ) -> tuple[float, float, float, float, float]:
        """
        Calculate the plasma masses.

        :param m_fuel_amu: Average mass of fuel (amu).
        :type m_fuel_amu: float
        :param m_ions_total_amu: Average mass of all ions (amu).
        :type m_ions_total_amu: float
        :param nd_ions_total: Total ion density (/m3).
        :type nd_ions_total: float
        :param nd_fuel_ions: Fuel ion density (/m3).
        :type nd_fuel_ions: float
        :param nd_alphas: Alpha ash density (/m3).
        :type nd_alphas: float
        :param vol_plasma: Plasma volume (m3).
        :type vol_plasma: float
        :param dene: Volume averaged electron density (/m3).
        :type dene: float

        :returns: A tuple containing:
        :rtype: tuple[float, float, float, float, float]
        """

        # Calculate mass of fuel ions
        m_plasma_fuel_ions = (m_fuel_amu * constants.atomic_mass_unit) * (
            nd_fuel_ions * vol_plasma
        )

        m_plasma_ions_total = (m_ions_total_amu * constants.atomic_mass_unit) * (
            nd_ions_total * vol_plasma
        )

        m_plasma_alpha = (nd_alphas * vol_plasma) * constants.alpha_mass

        m_plasma_electron = constants.electron_mass * (dene * vol_plasma)

        m_plasma = m_plasma_electron + m_plasma_ions_total

        return (
            m_plasma_fuel_ions,
            m_plasma_ions_total,
            m_plasma_alpha,
            m_plasma_electron,
            m_plasma,
        )


def calculate_poloidal_beta(btot, bp, beta):
    """Calculates total poloidal beta

    Author: James Morris (UKAEA)

    J.P. Freidberg, "Plasma physics and fusion energy", Cambridge University Press (2007)
    Page 270 ISBN 0521851076

    :param btot: sum of the toroidal and poloidal fields (T)
    :param bp: poloidal field (T)
    :param beta: total plasma beta
    """
    return beta * (btot / bp) ** 2


def res_diff_time(rmajor, res_plasma, kappa95):
    """Calculates resistive diffusion time

    Author: James Morris (UKAEA)

    :param rmajor: plasma major radius (m)
    :param res_plasma: plasma resistivity (Ohms)
    :param kappa95: plasma elongation at 95% flux surface
    """

    return 2 * constants.rmu0 * rmajor / (res_plasma * kappa95)


def l_h_threshold_power(
    nd_electron_line: float,
    bt: float,
    rmajor: float,
    rminor: float,
    kappa: float,
    a_plasma_surface: float,
    m_ions_total_amu: float,
    aspect: float,
    plasma_current: float,
) -> list[float]:
    """
    L-mode to H-mode power threshold calculation.

    :param nd_electron_line: Line-averaged electron density (/m3)
    :type nd_electron_line: float
    :param bt: Toroidal field on axis (T)
    :type bt: float
    :param rmajor: Plasma major radius (m)
    :type rmajor: float
    :param rminor: Plasma minor radius (m)
    :type rminor: float
    :param kappa: Plasma elongation
    :type kappa: float
    :param a_plasma_surface: Plasma surface area (m2)
    :type a_plasma_surface: float
    :param m_ions_total_amu: Average mass of all ions (amu)
    :type m_ions_total_amu: float
    :param aspect: Aspect ratio
    :type aspect: float
    :param plasma_current: Plasma current (A)
    :type plasma_current: float

    :returns: Array of power thresholds
    :rtype: list[float]
    """

    dnla20 = 1e-20 * nd_electron_line

    # ========================================================================

    # ITER-1996 H-mode power threshold database
    # Fit to 1996 H-mode power threshold database: nominal

    # i_l_h_threshold = 1
    iterdd = transition.calculate_iter1996_nominal(dnla20, bt, rmajor)

    # Fit to 1996 H-mode power threshold database: upper bound
    # i_l_h_threshold = 2
    iterdd_ub = transition.calculate_iter1996_upper(dnla20, bt, rmajor)

    # Fit to 1996 H-mode power threshold database: lower bound
    # i_l_h_threshold = 3
    iterdd_lb = transition.calculate_iter1996_lower(dnla20, bt, rmajor)

    # ========================================================================

    # Snipes 1997 ITER H-mode power threshold

    # i_l_h_threshold = 4
    snipes_1997 = transition.calculate_snipes1997_iter(dnla20, bt, rmajor)

    # i_l_h_threshold = 5
    snipes_1997_kappa = transition.calculate_snipes1997_kappa(dnla20, bt, rmajor, kappa)

    # ========================================================================

    # Martin et al (2008) for recent ITER scaling, with mass correction
    # and 95% confidence limits

    # i_l_h_threshold = 6
    martin_nominal = transition.calculate_martin08_nominal(
        dnla20, bt, a_plasma_surface, m_ions_total_amu
    )

    # i_l_h_threshold = 7
    martin_ub = transition.calculate_martin08_upper(
        dnla20, bt, a_plasma_surface, m_ions_total_amu
    )

    # i_l_h_threshold = 8
    martin_lb = transition.calculate_martin08_lower(
        dnla20, bt, a_plasma_surface, m_ions_total_amu
    )

    # ========================================================================

    # Snipes et al (2000) scaling with mass correction
    # Nominal, upper and lower

    # i_l_h_threshold = 9
    snipes_2000 = transition.calculate_snipes2000_nominal(
        dnla20, bt, rmajor, rminor, m_ions_total_amu
    )

    # i_l_h_threshold = 10
    snipes_2000_ub = transition.calculate_snipes2000_upper(
        dnla20, bt, rmajor, rminor, m_ions_total_amu
    )

    # i_l_h_threshold = 11
    snipes_2000_lb = transition.calculate_snipes2000_lower(
        dnla20, bt, rmajor, rminor, m_ions_total_amu
    )

    # ========================================================================

    # Snipes et al (2000) scaling (closed divertor) with mass correction
    # Nominal, upper and lower

    # i_l_h_threshold = 12
    snipes_2000_cd = transition.calculate_snipes2000_closed_divertor_nominal(
        dnla20, bt, rmajor, m_ions_total_amu
    )

    # i_l_h_threshold = 13
    snipes_2000_cd_ub = transition.calculate_snipes2000_closed_divertor_upper(
        dnla20, bt, rmajor, m_ions_total_amu
    )

    # i_l_h_threshold = 14
    snipes_2000_cd_lb = transition.calculate_snipes2000_closed_divertor_lower(
        dnla20, bt, rmajor, m_ions_total_amu
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
    hubbard_2017 = transition.calculate_hubbard2017(dnla20, a_plasma_surface, bt)

    # ========================================================================

    # Aspect ratio corrected Martin et al (2008)

    # i_l_h_threshold = 19
    martin_nominal_aspect = transition.calculate_martin08_aspect_nominal(
        dnla20, bt, a_plasma_surface, m_ions_total_amu, aspect
    )

    # i_l_h_threshold = 20
    martin_ub_aspect = transition.calculate_martin08_aspect_upper(
        dnla20, bt, a_plasma_surface, m_ions_total_amu, aspect
    )

    # i_l_h_threshold = 21
    martin_lb_aspect = transition.calculate_martin08_aspect_lower(
        dnla20, bt, a_plasma_surface, m_ions_total_amu, aspect
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


def reinke_tsep(bt, flh, qstar, rmajor, eps, fgw, kappa, lhat):
    """Function for calculating upstream temperature(keV) in Reinke model
    author: H Lux, CCFE/UKAEA
    bt      : input real : toroidal field on axis (T)
    flh     : input real : fraction of Psep/P_LH
    qstar   : input real : safety factor similar to q95 (see #707)
    rmajor  : input real : major radius (m)
    eps     : input real : inverse aspect ratio
    fgw     : input real : ratio of volume averaged density to n_GW
    kappa   : input real : elongation
    lhat    : input real : connection length factor
    This function calculates the upstream temperature in the
    divertor/SoL model used for the Reinke citerion.
    Issue #707
    M.L. Reinke 2017 Nucl. Fusion 57 034004
    """
    kappa_0 = 2.0e3  # Stangeby W/m/eV^(7/2)

    return (
        (bt**0.72 * flh**0.29 * fgw**0.21 * qstar**0.08 * rmajor**0.33)
        * (eps**0.15 * (1.0 + kappa**2.0) ** 0.34)
        * (lhat**0.29 * kappa_0 ** (-0.29) * 0.285)
    )


def init_physics_variables():
    physics_variables.m_beam_amu = 0.0
    physics_variables.m_fuel_amu = 0.0
    physics_variables.m_ions_total_amu = 0.0
    physics_variables.m_plasma_fuel_ions = 0.0
    physics_variables.m_plasma_ions_total = 0.0
    physics_variables.m_plasma_alpha = 0.0
    physics_variables.m_plasma_electron = 0.0
    physics_variables.m_plasma = 0.0
    physics_variables.alphaj = 1.0
    physics_variables.i_alphaj = 0
    physics_variables.alphan = 0.25
    physics_variables.alphap = 0.0
    physics_variables.fusden_alpha_total = 0.0
    physics_variables.fusden_plasma_alpha = 0.0
    physics_variables.alphat = 0.5
    physics_variables.aspect = 2.907
    physics_variables.beamfus0 = 1.0
    physics_variables.beta = 0.042
    physics_variables.beta_fast_alpha = 0.0
    physics_variables.beta_max = 0.0
    physics_variables.beta_min = 0.0
    physics_variables.beta_beam = 0.0
    physics_variables.beta_poloidal = 0.0
    physics_variables.beta_poloidal_eps = 0.0
    physics_variables.beta_toroidal = 0.0
    physics_variables.beta_thermal = 0.0
    physics_variables.beta_thermal_poloidal = 0.0
    physics_variables.beta_thermal_toroidal = 0.0
    physics_variables.beta_norm_total = 0.0
    physics_variables.beta_norm_thermal = 0.0
    physics_variables.beta_norm_poloidal = 0.0
    physics_variables.e_plasma_beta_thermal = 0.0
    physics_variables.beta_norm_toroidal = 0.0
    physics_variables.betbm0 = 1.5
    physics_variables.bp = 0.0
    physics_variables.bt = 5.68
    physics_variables.btot = 0.0
    physics_variables.burnup = 0.0
    physics_variables.burnup_in = 0.0
    physics_variables.bvert = 0.0
    physics_variables.c_beta = 0.5
    physics_variables.csawth = 1.0
    physics_variables.f_vol_plasma = 1.0
    physics_variables.f_r_conducting_wall = 1.35
    physics_variables.dene = 9.8e19
    physics_variables.nd_fuel_ions = 0.0
    physics_variables.dlamee = 0.0
    physics_variables.dlamie = 0.0
    physics_variables.dlimit[:] = 0.0
    physics_variables.nd_alphas = 0.0
    physics_variables.nd_beam_ions = 0.0
    physics_variables.nd_beam_ions_out = 0.0
    physics_variables.beta_norm_max = 3.5
    physics_variables.beta_norm_max_wesson = 0.0
    physics_variables.beta_norm_max_menard = 0.0
    physics_variables.beta_norm_max_original_scaling = 0.0
    physics_variables.beta_norm_max_tholerus = 0.0
    physics_variables.beta_norm_max_stambaugh = 0.0
    physics_variables.dnelimt = 0.0
    physics_variables.nd_ions_total = 0.0
    physics_variables.nd_electron_line = 0.0
    physics_variables.nd_protons = 0.0
    physics_variables.ntau = 0.0
    physics_variables.nTtau = 0.0
    physics_variables.nd_impurities = 0.0
    physics_variables.beta_poloidal_eps_max = 1.38
    physics_variables.eps = 0.34399724802
    physics_variables.f_c_plasma_auxiliary = 0.0
    physics_variables.f_c_plasma_inductive = 0.0
    physics_variables.f_alpha_electron = 0.0
    physics_variables.f_p_alpha_plasma_deposited = 0.95
    physics_variables.f_alpha_ion = 0.0
    physics_variables.f_deuterium = 0.5
    physics_variables.f_p_div_lower = 1.0
    physics_variables.ffwal = 0.92
    physics_variables.fgwped = 0.85
    physics_variables.fgwsep = 0.50
    physics_variables.f_helium3 = 0.0
    physics_variables.figmer = 0.0
    physics_variables.fkzohm = 1.0
    physics_variables.fplhsep = 1.0
    physics_variables.fpdivlim = 1.0
    physics_variables.fne0 = 1.0
    physics_variables.f_tritium = 0.5
    physics_variables.fusden_total = 0.0
    physics_variables.fusrat_total = 0.0
    physics_variables.fusden_plasma = 0.0
    physics_variables.f_c_plasma_non_inductive = 1.0
    physics_variables.ejima_coeff = 0.4
    physics_variables.f_beta_alpha_beam_thermal = 0.0
    physics_variables.hfac[:] = 0.0
    physics_variables.hfact = 1.0
    physics_variables.taumax = 10.0
    physics_variables.i_bootstrap_current = 3
    physics_variables.i_beta_component = 0
    physics_variables.i_plasma_current = 4
    physics_variables.i_diamagnetic_current = 0
    physics_variables.i_density_limit = 8
    physics_variables.n_divertors = 2
    physics_variables.i_beta_fast_alpha = 1
    physics_variables.i_plasma_ignited = 0
    physics_variables.ipedestal = 1
    physics_variables.i_pfirsch_schluter_current = 0
    physics_variables.neped = 4.0e19
    physics_variables.nesep = 3.0e19
    physics_variables.alpha_crit = 0.0
    physics_variables.nesep_crit = 0.0
    physics_variables.plasma_res_factor = 1.0
    physics_variables.rhopedn = 1.0
    physics_variables.rhopedt = 1.0
    physics_variables.tbeta = 2.0
    physics_variables.teped = 1.0
    physics_variables.tesep = 0.1
    physics_variables.i_beta_norm_max = 1
    physics_variables.i_rad_loss = 1
    physics_variables.i_confinement_time = 34
    physics_variables.i_plasma_wall_gap = 1
    physics_variables.i_plasma_geometry = 0
    physics_variables.i_plasma_shape = 0
    physics_variables.itart = 0
    physics_variables.itartpf = 0
    physics_variables.iwalld = 1
    physics_variables.plasma_square = 0.0
    physics_variables.kappa = 1.792
    physics_variables.kappa95 = 1.6
    physics_variables.kappa_ipb = 0.0
    physics_variables.ne0 = 0.0
    physics_variables.ni0 = 0.0
    physics_variables.m_s_limit = 0.3
    physics_variables.p0 = 0.0
    physics_variables.j_plasma_0 = 0.0
    physics_variables.f_dd_branching_trit = 0.0
    physics_variables.pden_plasma_alpha_mw = 0.0
    physics_variables.pden_alpha_total_mw = 0.0
    physics_variables.f_pden_alpha_electron_mw = 0.0
    physics_variables.p_fw_alpha_mw = 0.0
    physics_variables.f_pden_alpha_ions_mw = 0.0
    physics_variables.p_alpha_total_mw = 0.0
    physics_variables.p_plasma_alpha_mw = 0.0
    physics_variables.p_beam_alpha_mw = 0.0
    physics_variables.p_beam_neutron_mw = 0.0
    physics_variables.p_beam_dt_mw = 0.0
    physics_variables.p_non_alpha_charged_mw = 0.0
    physics_variables.pden_non_alpha_charged_mw = 0.0
    physics_variables.pcoef = 0.0
    physics_variables.p_plasma_inner_rad_mw = 0.0
    physics_variables.pden_plasma_core_rad_mw = 0.0
    physics_variables.p_dd_total_mw = 0.0
    physics_variables.p_dhe3_total_mw = 0.0
    physics_variables.p_plasma_separatrix_mw = 0.0
    physics_variables.pdivl = 0.0
    physics_variables.pdivu = 0.0
    physics_variables.pdivmax = 0.0
    physics_variables.p_dt_total_mw = 0.0
    physics_variables.p_plasma_dt_mw = 0.0
    physics_variables.p_plasma_outer_rad_mw = 0.0
    physics_variables.pden_plasma_outer_rad_mw = 0.0
    physics_variables.p_charged_particle_mw = 0.0
    physics_variables.vs_plasma_internal = 0.0
    physics_variables.pflux_fw_rad_mw = 0.0
    physics_variables.pden_ion_electron_equilibration_mw = 0.0
    physics_variables.plasma_current = 0.0
    physics_variables.p_plasma_neutron_mw = 0.0
    physics_variables.p_neutron_total_mw = 0.0
    physics_variables.pden_neutron_total_mw = 0.0
    physics_variables.pden_plasma_neutron_mw = 0.0
    physics_variables.p_plasma_ohmic_mw = 0.0
    physics_variables.pden_plasma_ohmic_mw = 0.0
    physics_variables.p_plasma_loss_mw = 0.0
    physics_variables.p_fusion_total_mw = 0.0
    physics_variables.len_plasma_poloidal = 0.0
    physics_variables.p_plasma_rad_mw = 0.0
    physics_variables.pden_plasma_rad_mw = 0.0
    physics_variables.pradsolmw = 0.0
    physics_variables.proton_rate_density = 0.0
    physics_variables.psolradmw = 0.0
    physics_variables.pden_plasma_sync_mw = 0.0
    physics_variables.p_plasma_sync_mw = 0.0
    physics_variables.i_l_h_threshold = 19
    physics_variables.p_l_h_threshold_mw = 0.0
    physics_variables.l_h_threshold_powers[:] = 0.0
    physics_variables.p_electron_transport_loss_mw = 0.0
    physics_variables.pden_electron_transport_loss_mw = 0.0
    physics_variables.p_ion_transport_loss_mw = 0.0
    physics_variables.pscalingmw = 0.0
    physics_variables.pden_ion_transport_loss_mw = 0.0
    physics_variables.q0 = 1.0
    physics_variables.q95 = 0.0
    physics_variables.qfuel = 0.0
    physics_variables.tauratio = 1.0
    physics_variables.q95_min = 0.0
    physics_variables.qstar = 0.0
    physics_variables.rad_fraction_sol = 0.8
    physics_variables.rad_fraction_total = 0.0
    physics_variables.f_nd_alpha_electron = 0.1
    physics_variables.f_nd_protium_electrons = 0.0
    physics_variables.ind_plasma_internal_norm = 0.9
    physics_variables.ind_plasma_internal_norm_wesson = 0.0
    physics_variables.ind_plasma_internal_norm_menard = 0.0
    physics_variables.i_ind_plasma_internal_norm = 0
    physics_variables.ind_plasma = 0.0
    physics_variables.rmajor = 8.14
    physics_variables.rminor = 0.0
    physics_variables.f_nd_beam_electron = 0.005
    physics_variables.rncne = 0.0
    physics_variables.rndfuel = 0.0
    physics_variables.rnfene = 0.0
    physics_variables.rnone = 0.0
    physics_variables.f_res_plasma_neo = 0.0
    physics_variables.res_plasma = 0.0
    physics_variables.t_plasma_res_diffusion = 0.0
    physics_variables.a_plasma_surface = 0.0
    physics_variables.a_plasma_surface_outboard = 0.0
    physics_variables.i_single_null = 1
    physics_variables.f_sync_reflect = 0.6
    physics_variables.t_electron_energy_confinement = 0.0
    physics_variables.tauee_in = 0.0
    physics_variables.t_energy_confinement = 0.0
    physics_variables.t_ion_energy_confinement = 0.0
    physics_variables.t_alpha_confinement = 0.0
    physics_variables.f_alpha_energy_confinement = 0.0
    physics_variables.te = 12.9
    physics_variables.te0 = 0.0
    physics_variables.ten = 0.0
    physics_variables.ti = 12.9
    physics_variables.ti0 = 0.0
    physics_variables.tin = 0.0
    physics_variables.tratio = 1.0
    physics_variables.triang = 0.36
    physics_variables.triang95 = 0.24
    physics_variables.vol_plasma = 0.0
    physics_variables.vs_plasma_burn_required = 0.0
    physics_variables.vs_plasma_ramp_required = 0.0
    physics_variables.v_plasma_loop_burn = 0.0
    physics_variables.vshift = 0.0
    physics_variables.vs_plasma_ind_ramp = 0.0
    physics_variables.vs_plasma_res_ramp = 0.0
    physics_variables.vs_plasma_total_required = 0.0
    physics_variables.pflux_fw_neutron_mw = 0.0
    physics_variables.wtgpd = 0.0
    physics_variables.a_plasma_poloidal = 0.0
    physics_variables.zeff = 0.0
    physics_variables.zeffai = 0.0


def init_physics_module():
    """Initialise the physics module"""
    physics_module.first_call = 1
    physics_module.iscz = 0
    physics_module.err242 = 0
    physics_module.err243 = 0
    physics_module.rad_fraction_lcfs = 0.0
    physics_module.e_plasma_beta = 0.0
    physics_module.total_loss_power = 0.0
    physics_module.t_energy_confinement_beta = 0.0
    physics_module.ptarmw = 0.0
    physics_module.lambdaio = 0.0
    physics_module.drsep = 0.0
    physics_module.fio = 0.0
    physics_module.fli = 0.0
    physics_module.flo = 0.0
    physics_module.fui = 0.0
    physics_module.fuo = 0.0
    physics_module.plimw = 0.0
    physics_module.plomw = 0.0
    physics_module.puimw = 0.0
    physics_module.puomw = 0.0
    physics_module.rho_star = 0.0
    physics_module.nu_star = 0.0
    physics_module.beta_mcdonald = 0.0
    physics_module.itart_r = 0.0
