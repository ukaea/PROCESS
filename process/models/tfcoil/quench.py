"""Module for quench protection calculations and material properties for TF coils."""

import logging
from typing import Final

import numpy as np

# TODO: Use of CoolProp prevents nb.jit at present...
from process.core.coolprop_interface import FluidProperties

__all__ = ["calculate_quench_protection_current_density"]


logger = logging.getLogger(__name__)

# Material property parameterisations

COPPER_DENSITY = 8960.0  # [kg/m³]
NB3SN_DENSITY = 8040.0  # [kg/m³]


def _copper_specific_heat_capacity(temperature: float) -> float:
    """Calculates the specific heat capacity of cryogenic copper at a given
    temperature [J/(kg·K)].

    Parameters
    ----------
    temperature:
        Temperature [K].

    Returns
    -------
    float
        Specific heat capacity of copper at the given temperature [J/(kg·K)].

    Notes
    -----
    Assumes high-purity copper with negligible impurity effects.

    References
    ----------
    [1] J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of
        Copper and Copper Alloys at Cryogenic Temperatures",
        U.S. Government Printing Office, February 1992.
        https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nistmonograph177.pdf
        Equation 7-1
    """
    poly_coeffs: Final[list[float]] = [1.131, -9.454, 12.99, -5.501, 0.7637]
    logt = np.log10(temperature)
    logcp = sum(c * logt**i for i, c in enumerate(poly_coeffs))
    return 10**logcp


def _copper_rrr_resistivity(temperature: float, rrr: float) -> float:
    """Calculates the electrical resistivity of cryogenic copper with temperature
    and RRR dependence  [Ω·m].

    Parameters
    ----------
    temperature:
        Operating temperature [K].
    rrr:
        Residual resistivity ratio (dimensionless).

    Returns
    -------
    float
        Electrical resistivity of copper [Ω·m].

    References
    ----------
    [1] J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of
    Copper and Copper Alloys at Cryogenic Temperatures", U.S. Government Printing Office,
    February 1992. https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nistmonograph177.pdf
    Equation 8-1

    [2] J. G. Hust, A. B. Lankford, NBSIR 84-3007 "THERMAL CONDUCTIVITY OF ALUMINUM,
    COPPER, IRON, AND TUNGSTEN FOR TEMPERATURES FROM  1 K TO THE MELTING POINT", 1984
    https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir84-3007.pdf
    """
    p1: Final[float] = 1.171e-17
    p2: Final[float] = 4.49
    p3: Final[float] = 3.841e10
    p4: Final[
        float
    ] = -1.14  # Appears to be a typo in the original papers... given as 1.14
    p5: Final[float] = 50.0
    p6: Final[float] = 6.428
    p7: Final[float] = 0.4531
    # Experimentally determined deviation from Mattheisen... (not given anywhere)
    rho_c: Final[float] = 0.0
    p9: Final[float] = 1.553e-8

    t = temperature

    # Compute rho_o
    rho_o = p9 / rrr
    # Compute rho_i
    numerator = p1 * t**p2
    denominator = 1.0 + p1 * p3 * t ** (p2 + p4) * np.exp(-((p5 / t) ** p6))
    rho_i = numerator / denominator + rho_c
    rho_io = p7 * rho_i * rho_o / (rho_i + rho_o)

    # Compute rho
    return rho_o + rho_i + rho_io


def _copper_irradiation_resistivity(fluence: float) -> float:
    """Calculates the radiation-induced electrical resistivity of copper.

    Estimates the increase in copper's electrical resistivity [Ω·m] due to neutron
    irradiation, as a function of total neutron fluence.


    Parameters
    ----------
    fluence:
        Total neutron fluence [n/m²].

    Returns
    -------
    :
        Radiation-induced resistivity of copper [Ω·m].

    Notes
    -----
    - Fit to data at low-temperature conditions (around 4.6 K).
    - Fit to data in fast neutron spectrum (E > 0.1 MeV).
    - Damage and transmutation effects therefore both included, but transmutation
      effects may be underestimated.
    - This is an additive contribution to the base residual resistivity of copper.

    References
    ----------
    [1] M. Kovari, 09/11/2012, internal notes (Excel / Mathcad), Technology Program,
    WP12, PEX, Super-X Divertor for DEMO.

    [2] M. Nakagawa et al., "High-dose neutron-irradiation effects in fcc metals at
    4.6 K", *Phys. Rev. B*, 16, 5285 (1977). https://doi.org/10.1103/PhysRevB.16.5285
    Figure 6
    """
    c1: Final[float] = 0.00283
    c2: Final[float] = -0.0711
    c3: Final[float] = 0.77982
    res_scale: Final[float] = 1e-9
    flu_scale: Final[float] = 1e-22

    fluence_norm = flu_scale * fluence

    return res_scale * (c1 * fluence_norm**3 + c2 * fluence_norm**2 + c3 * fluence_norm)


def _copper_magneto_resistivity(resistivity: float, field: float) -> float:
    """Calculates the electrical resistivity of cryogenic copper due to magnetoresistive
    effects [Ω·m].

    Parameters
    ----------
    resistivity:
        Operating resistivity [K].
    field:
        Operating magnetic field [T].

    Returns
    -------
    :
        Electrical resistivity of copper [Ω·m].

    Notes
    -----
        - Resistivity increases with magnetic field due to magnetoresistance effects.

    References
    ----------
    [1] J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of
    Copper and Copper Alloys at Cryogenic Temperatures", U.S. Government Printing Office,
    February 1992. https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nistmonograph177.pdf
    Equation 8-7
    """
    p9: Final[float] = 1.553e-8

    # TODO: This feels strange, but cut-off necessary, and B < 1.0 is "possible" and
    # well-behaved at low T
    if field > 1e-2:
        poly_coeffs: Final[list[float]] = [-2.662, 0.3168, 0.6229, -0.1839, 0.01827]

        x = np.log10(p9 * field / resistivity)
        a = sum(c * x**i for i, c in enumerate(poly_coeffs))
        return resistivity * (1.0 + 10**a)

    return resistivity


def _copper_electrical_resistivity(
    temperature: float, field: float, rrr: float, fluence: float
) -> float:
    """Calculates the electrical resistivity of cryogenic copper with temperature, RRR,
    magnetic field, and fluence dependence  [Ω·m].

    Parameters
    ----------
    temperature:
        Operating temperature [K].
    rrr:
        Residual resistivity ratio (dimensionless).
    field:
        Operating magnetic field [T].
    fluence:
        Total end-of-life neutron fluence [n/m²].

    Returns
    -------
    float
        Electrical resistivity of copper [Ω·m].

    Notes
    -----
    - Resistivity increases with magnetic field due to magnetoresistance effects.

    References
    ----------
    [1] J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of
    Copper and Copper Alloys at Cryogenic Temperatures", U.S. Government Printing Office,
    February 1992. https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/nistmonograph177.pdf
    Equation 8-1

    [2] J. G. Hust, A. B. Lankford, NBSIR 84-3007 "THERMAL CONDUCTIVITY OF  ALUMINUM,
    COPPER, IRON, AND TUNGSTEN FOR TEMPERATURES FROM  1 K TO THE MELTING POINT", 1984
    https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir84-3007.pdf

    [3] M. Kovari, 09/11/2012, internal notes (Excel / Mathcad), Technology Program,
    WP12, PEX, Super-X Divertor for DEMO.

    [4] M. Nakagawa et al., "High-dose neutron-irradiation effects in fcc metals at
    4.6 K", *Phys. Rev. B*, 16, 5285 (1977). https://doi.org/10.1103/PhysRevB.16.5285
    Figure 6
    """
    rho_rrr = _copper_rrr_resistivity(temperature, rrr)
    rho_irr = _copper_irradiation_resistivity(fluence)
    return _copper_magneto_resistivity(rho_rrr + rho_irr, field)


def _nb3sn_specific_heat_capacity(temperature: float) -> float:
    """Calculates the specific heat capacity of Nb₃Sn as a function of temperature.

    Provides the temperature-dependent specific heat capacity [J/(kg·K)] of the A15
    superconductor Nb₃Sn.

    Parameters
    ----------
    temperature:
        Temperature [K].

    Returns
    -------
    float
        Specific heat capacity of Nb₃Sn at the given temperature [J/(kg·K)].

    Notes
    -----
    - The superconducting part is ignored, which is typical in thermal quench
      calculations.

    References
    ----------
    [1] EFDA Material Data Compilation for Superconductor Simulation, P. Bauer,
    H. Rajainmaki, E. Salpietro, EFDA CSU, Garching, 04/18/07.

    [2] ITER DRG1 Annex, Superconducting Material Database, Article 5,
    N 11 FDR 42 01-07-05 R 0.1.

    [3] V.D. Arp, Stability and Thermal Quenches in Force-Cooled Superconducting Cables,
    Superconducting MHD Magnet Design Conf., MIT, pp 142-157, 1980.

    [4] G.S. Knapp, S.D. Bader, Z. Fisk, Phonon properties of A-15 superconductors
    obtained from heat capacity measurements, Phys. Rev. B, 13(9), pp 3783-3789, 1976.
    https://doi.org/10.1103/PhysRevB.13.3783
    """
    gamma: Final[float] = 0.1  # [J/K²/kg] (Grueneisen)
    beta: Final[float] = 0.001  # [J/K⁴/kg] (Debye)
    cp_300: Final[float] = 210.0  # [J/K/kg] Room-temperature specific heat

    # Normally conducting specific heat capacity (i.e. ignoring transition from
    # superconducting state)
    cp_low = beta * temperature**3 + gamma * temperature
    return 1.0 / (1.0 / cp_300 + 1.0 / cp_low)


# Quench model
# Gauss-Legendre quadrature nodes and weights for numerical integration
# 75 points is a good compromise between speed and accuracy for the integrals used
# in this model.
GAUSS_LEG_NODES, GAUSS_LEG_WEIGHTS = np.polynomial.legendre.leggauss(75)


def _quench_integrand_at_temperature(
    temperature: float,
    field: float,
    rrr: float,
    fluence: float,
    pressure: float = 6e5,
) -> tuple[float, float, float]:
    """Calculates the per-unit-temperature integrands for quench protection.

    Returns the volumetric heat capacity divided by copper resistivity for each
    material component at a given temperature.

    Parameters
    ----------
    temperature:
        Temperature [K].
    field:
        Magnetic field [T].
    rrr:
        Residual resistivity ratio of copper.
    fluence:
        Neutron fluence [n/m²] (for irradiation effects).
    pressure:
        Helium pressure [Pa]. Defaults to 6e5 (ITER TF coolant pressure).

    Returns
    -------
    tuple[float, float, float]
        Tuple of (helium, copper, superconductor) integrands [J·s/(m³·Ω·m·K)]
        i.e. (rho·cp / rho_cu) for each material at the given temperature.

    Notes
    -----
    - Helium is assumed to be at constant pressure.
    """
    nu_cu = _copper_electrical_resistivity(temperature, field, rrr, fluence)
    he_properties = FluidProperties.of("He", temperature=temperature, pressure=pressure)

    ihe_integrand = he_properties.specific_heat_const_p * he_properties.density / nu_cu
    icu_integrand = _copper_specific_heat_capacity(temperature) * COPPER_DENSITY / nu_cu
    isc_integrand = _nb3sn_specific_heat_capacity(temperature) * NB3SN_DENSITY / nu_cu

    return ihe_integrand, icu_integrand, isc_integrand


def _quench_integrals(
    t_he_peak: float,
    temp_quench_max: float,
    field: float,
    rrr: float,
    fluence: float,
) -> tuple[float, float, float]:
    """Calculates the material property integrals for quench protection.

    Evaluates the integrals over temperature for helium, copper, and superconductor
    contributions. These integrals are used in determining current density limits
    during a quench.

    Parameters
    ----------
    t_he_peak:
        t_he_peak: Lower temperature bound of integration [K].
    temp_quench_max:
        temp_quench_max: Upper temperature bound of integration [K].
    field:
        field: Magnetic field [T].
    rrr:
        rrr: Residual resistivity ratio of copper.
    fluence:
        fluence: Neutron fluence [n/cm²] (for irradiation effects).

    Returns
    -------
    Tuple[float, float, float]
        Tuple of integrals for helium, copper, and superconductor contributions
        (I_He, I_Cu, I_sc).

    Notes
    -----
    - Integrals assume temperature-dependent material models are defined for the entire
      range [t_he_peak, temp_quench_max].
    - Helium is assumed to be at constant pressure throughout the quench
      (i.e. some PRV in the cooling system)
    """
    # Helium pressure [Pa] (assumed to be constant throughout quench) - no plans to
    # make input
    pressure = 6e5  # ITER TF coolant pressure

    ihe = 0.0
    icu = 0.0
    isc = 0.0

    for xi, wi in zip(GAUSS_LEG_NODES, GAUSS_LEG_WEIGHTS, strict=False):
        ti = 0.5 * (xi + 1.0) * (temp_quench_max - t_he_peak) + t_he_peak
        dti = 0.5 * wi * (temp_quench_max - t_he_peak)

        ihe_i, icu_i, isc_i = _quench_integrand_at_temperature(
            ti, field, rrr, fluence, pressure
        )

        ihe += dti * ihe_i
        icu += dti * icu_i
        isc += dti * isc_i

    return ihe, icu, isc


def _build_cumulative_quench_integral(
    temp_he_peak: float,
    temp_quench_max: float,
    field: float,
    rrr: float,
    fluence: float,
    f_a_cable_space_helium: float,
    f_cu_cable: float,
    f_sc_cable: float,
    n_temp: int = 300,
) -> tuple[np.ndarray, np.ndarray]:
    """Builds a cumulative temperature integral lookup for hotspot temperature inversion.

    Computes the cumulative integral of the combined volumetric heat capacity
    (weighted by material fractions) divided by copper resistivity, from
    `temp_he_peak` to `temp_quench_max`. This lookup is used to map
    accumulated (∫J² dt) to hotspot temperature.

    Parameters
    ----------
    temp_he_peak:
        Lower temperature bound [K].
    temp_quench_max:
        Upper temperature bound [K].
    field:
        Magnetic field [T].
    rrr:
        Residual resistivity ratio of copper.
    fluence:
        Neutron fluence [n/m²].
    f_a_cable_space_helium:
        Fraction of cable space occupied by helium.
    f_cu_cable:
        Effective copper fraction of cable cross-section.
    f_sc_cable:
        Effective superconductor fraction of cable cross-section.
    n_temp:
        Number of temperature points in the lookup array.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        `(temp_array, cum_integral)` where `cum_integral[k]` is the
        cumulative integral from `temp_he_peak` to `temp_array[k]`.
    """
    pressure = 6e5  # ITER TF coolant pressure

    # Evaluate the integrand on a fine uniform grid and integrate with the
    # trapezoidal rule.  This is much cheaper than nested Gauss-Legendre
    # quadrature because each temperature point is evaluated exactly once
    # (avoiding repeated CoolProp calls) and the result is accumulated with
    # np.cumsum - which is important when this function is called multiple
    # times (e.g. from plot_quench_time_evolution).
    temp_array = np.linspace(temp_he_peak, temp_quench_max, n_temp)
    integrand = np.empty(n_temp)
    for k, t in enumerate(temp_array):
        ihe_i, icu_i, isc_i = _quench_integrand_at_temperature(
            t, field, rrr, fluence, pressure
        )
        integrand[k] = (
            f_a_cable_space_helium * ihe_i + f_cu_cable * icu_i + f_sc_cable * isc_i
        )

    # Cumulative trapezoid integration
    dt = temp_array[1] - temp_array[0]
    cum_integral = np.zeros(n_temp)
    cum_integral[1:] = np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * dt)

    return temp_array, cum_integral


def calculate_quench_protection_current_density(
    tau_discharge: float,
    b_peak: float,
    f_a_cable_copper: float,
    f_a_cable_space_helium: float,
    temp_he_peak: float,
    temp_quench_max: float,
    cu_rrr: float,
    t_quench_detection: float,
    fluence: float,
) -> float:
    """Calculates the current density limited by the protection limit.

    Simplified 0-D adiabatic heat balance "hotspot criterion" model.

    Uses temperature- and magnetic field-dependent material properties to determine the
    maximum allowable winding pack current density during a quench. Assumes that the
    magnetic field does not decay over time and accounts for contributions from copper,
    helium, and superconductor materials using temperature integrals.


    Parameters
    ----------
    tau_discharge:
        tau_discharge: Quench discharge time constant [s].
    b_peak:
        b_peak: Magnetic field at the peak point [T].
    f_a_cable_copper:
        f_a_cable_copper: Fraction of cable cross-section that is copper.
    f_a_cable_space_helium:
        f_a_cable_space_helium: Fraction of cable space occupied by helium.
    temp_he_peak:
        temp_he_peak: Peak helium temperature at quench initiation [K].
    temp_quench_max:
        temp_quench_max: Maximum allowed conductor temperature during quench [K].
    cu_rrr:
        cu_rrr: Residual resistivity ratio of copper.
    t_quench_detection:
        t_quench_detection: Detection time delay [s].
    fluence:
        fluence: Neutron fluence [n/m²].

    Returns
    -------
    float
        Maximum allowable winding pack current density [A/m²].

    Notes
    -----
    - Assumes constant magnetic field over the duration of the quench.
    - Assumes the dump resistor has a constant resistance much higher
    than that of the TF coil.
    - Operates on the current-carring cross-section of a conductor. The
    jacket and insulation are ignored. The actual allowable WP current
    density must be weighted with the ratio of current-carrying cross-section
    vs. total WP cross-section (including jacket and insulation).
    - Presently only applicable to LTS TF coil winding packs (Nb3Sn assumed)
    """
    # Default fluence is too high for this model
    if (fluence < 0.0) | (fluence > 1.5e23):
        logger.warning(
            "Fluence values out of range [0.0, 1.5e23]; kludging.", stacklevel=2
        )
        fluence = np.clip(fluence, 0.0, 1.5e23)

    i_he, i_cu, i_sc = _quench_integrals(
        temp_he_peak, temp_quench_max, b_peak, cu_rrr, fluence
    )

    f_cu_cable = (1.0 - f_a_cable_space_helium) * f_a_cable_copper
    f_sc_cable = (1.0 - f_a_cable_space_helium) * (1.0 - f_a_cable_copper)

    factor = 1.0 / (0.5 * tau_discharge + t_quench_detection)
    total_integral = (
        f_a_cable_space_helium * i_he + f_cu_cable * i_cu + f_sc_cable * i_sc
    )

    return np.sqrt(factor * f_cu_cable * total_integral)
