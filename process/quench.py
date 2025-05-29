from typing import Final

import numpy as np

# Material property parameterisations


def copper_density(temperature: float) -> float:  # noqa: ARG001
    """
    Calculate the density of cryogenic copper [kg/m^3].
    """
    return 8960.0  # No plans to include T-dependence


def copper_specific_heat_capacity(temperature: float) -> float:
    """
    Calculates the specific heat capacity of cryogenic copper at a given temperature [J/(kg·K)].

    :author M. Coleman, UKAEA

    :param float temperature: Temperature [K].
    :returns: Specific heat capacity of copper at the given temperature [J/(kg·K)].
    :rtype: float

    :notes:
        - Assumes high-purity copper with negligible impurity effects.

    :references:
        - L. Dresner, *Stability of Superconductors*, Plenum Press, New York, 1995.
    """

    gamma: Final[float] = 0.011  # [J/K²/kg] (Grueneisen)
    beta: Final[float] = 0.0011  # [J/K⁴/kg] (Debye)
    cp_300: Final[float] = 385.491  # [J/K/kg] Room-temperature specific heat

    cp_low = beta * temperature**3 + gamma * temperature
    return 1.0 / (1.0 / cp_300 + 1.0 / cp_low)


def copper_electrical_resistivity(
    temperature: float, field: float, rrr: float
) -> float:
    """
    Calculates the electrical resistivity of cryogenic copper with temperature and magnetic
    field dependence  [Ω·m].

    :author M. Coleman, UKAEA

    :param float temperature: Operating temperature [K].
    :param float field: Operating magnetic field [T].
    :param float rrr: Residual resistivity ratio (dimensionless).
    :returns: Electrical resistivity of copper [Ω·m].
    :rtype: float
    :raises ValueError: If any input is unphysical or outside the supported model range.

    :notes:
        - Resistivity increases with magnetic field due to magnetoresistance effects.

    :references:
        - J. Simon, E. S. Drexler, and R. P. Reed, *NIST Monograph 177*, "Properties of Copper and Copper Alloys
        at Cryogenic Temperatures", U.S. Government Printing Office, February 1992.
        - EFDA Material Data Compilation for Superconductor Simulation, P. Bauer, H. Rajainmaki,
        E. Salpietro, EFDA CSU, Garching, 04/18/07.
    """
    # Constants from EFDA documentation (Page 2-17)
    p1: Final[float] = 1.171e-17
    p2: Final[float] = 4.49
    p3: Final[float] = 4.5e-7
    p4: Final[float] = 3.35
    p5: Final[float] = 50.0
    p6: Final[float] = 6.428
    p7: Final[float] = 1.69e-8
    p8: Final[float] = 0.4531
    p9: Final[float] = 1.553e-8
    poly_coeffs: Final[list[float]] = [-2.662, 0.3168, 0.6229, -0.1839, 0.01827]

    # TODO: Implement PROCESS standard for kludging / warning
    t = np.clip(temperature, 4.0)
    # TODO: this function at present has the potential to return poor values at low
    # field...
    rrr = np.clip(rrr, 1.0)

    t = temperature

    # Compute rho1
    numerator = p1 * t**p2
    denominator = 1.0 + p3 * t**p4 * np.exp(-((p5 / t) ** p6))
    rho1 = numerator / denominator

    # Compute rho2
    rho2 = p7 / rrr + rho1 + p8 * (p7 * rho1) / (rrr * rho1 + p7)

    # Compute magnetic field correction factor
    x = np.log10(p9 * field / rho2)
    a = sum(c * x**i for i, c in enumerate(poly_coeffs))

    # Final resistivity
    return rho2 * (1.0 + 10**a)


def copper_irradiation_resistivity(fluence: float) -> float:
    """
    Calculates the radiation-induced electrical resistivity of copper.

    Estimates the increase in copper's electrical resistivity [Ω·m] due to neutron irradiation,
    as a function of total neutron fluence.

    :author M. Coleman, UKAEA

    :param float fluence: Total neutron fluence [n/m²].
    :returns: Radiation-induced resistivity of copper [Ω·m].
    :rtype: float
    :raises ValueError: If fluence is negative or exceeds empirical model limits.

    :notes:
        - Assumes low-temperature conditions (around 4.6 K).
        - This does not account for transmutation effects.
        - This is an additive contribution to the base residual resistivity of copper.

    :references:
        - M. Kovari, 09/11/2012, internal notes (Excel / Mathcad), Technology Program, WP12, PEX, Super-X Divertor for DEMO.
        - M. Nakagawa et al., "High-dose neutron-irradiation effects in fcc metals at 4.6 K", *Phys. Rev. B*, 16, 5285 (1977). https://doi.org/10.1103/PhysRevB.16.5285
    """
    # TODO: Check with M. Kovari and document
    c1: Final[float] = 0.00283
    c2: Final[float] = -0.0711
    c3: Final[float] = 0.77982
    res_scale: Final[float] = 1e-9
    flu_scale: Final[float] = 1e-22

    fluence_norm = flu_scale * fluence

    return res_scale * (c1 * fluence_norm**3 + c2 * fluence_norm**2 + c3 * fluence_norm)


def nb3sn_density(temperature: float) -> float:  # noqa: ARG001
    """
    Calculate the density of Nb3Sn [kg/m^3].
    """
    return 8040.0  # No plans to include T-dependence


def nb3sn_specific_heat_capacity(temperature: float) -> float:
    """
    Calculates the specific heat capacity of Nb₃Sn as a function of temperature.

    Provides the temperature-dependent specific heat capacity [J/(kg·K)] of the A15
    superconductor Nb₃Sn.

    :author M. Coleman, UKAEA

    :param float temperature: Temperature [K].
    :returns: Specific heat capacity of Nb₃Sn at the given temperature [J/(kg·K)].
    :rtype: float
    :raises ValueError: If the temperature is outside the supported range for this model.

    :notes:
        - The superconducting part is ignored, which is typical in thermal quench calculations.
        - Assumes polycrystalline, isotropic Nb₃Sn.

    :references:
        - EFDA Material Data Compilation for Superconductor Simulation, P. Bauer, H. Rajainmaki, E. Salpietro, EFDA CSU, Garching, 04/18/07.
        - ITER DRG1 Annex, Superconducting Material Database, Article 5, N 11 FDR 42 01-07-05 R 0.1.
        - V.D. Arp, Stability and Thermal Quenches in Force-Cooled Superconducting Cables, Superconducting MHD Magnet Design Conf., MIT, pp 142-157, 1980.
        - G.S. Knapp, S.D. Bader, Z. Fisk, Phonon properties of A-15 superconductors obtained from heat capacity measurements, Phys. Rev. B, 13(9), pp 3783-3789, 1976.
    """
    gamma: Final[float] = 0.1  # [J/K²/kg] (Grueneisen)
    beta: Final[float] = 0.001  # [J/K⁴/kg] (Debye)
    cp_300: Final[float] = 210.0  # [J/K/kg] Room-temperature specific heat

    # TODO: Apply PROCESS-style kludging
    temperature = np.clip(temperature, 2.0, 300.0)

    # Normally conducting specific heat capacity (i.e. ignoring transition from
    # superconducting state)
    cp_low = beta * temperature**3 + gamma * temperature
    return 1.0 / (1.0 / cp_300 + 1.0 / cp_low)


def helium_density(temperature: float) -> float:  # noqa: ARG001
    # TODO: Replace with T-dependent formulation
    return 150.0


def helium_specific_heat_capacity(temperature: float) -> float:  # noqa: ARG001
    # TODO: Replace with T-dependent formulation
    return 4750.0


# Quench model


def _quench_integrals(
    t_he_peak: float, t_max: float, field: float, rrr: float, fluence: float
) -> tuple[float, float, float]:
    """
    Calculates the material property integrals for quench protection.

    Evaluates the integrals over temperature for helium, copper, and superconductor contributions.
    These integrals are used in determining current density limits during a quench.

    :author M. Coleman, UKAEA

    :param float t_he_peak: Lower temperature bound of integration [K].
    :param float t_max: Upper temperature bound of integration [K].
    :param float field: Magnetic field [T].
    :param float rrr: Residual resistivity ratio of copper.
    :param float fluence: Neutron fluence [n/cm²] (for irradiation effects).
    :returns: Tuple of integrals for helium, copper, and superconductor contributions (I_He, I_Cu, I_sc).
    :rtype: Tuple[float, float, float]
    :raises ValueError: If temperature bounds are invalid or parameters are unphysical.

    :notes:
        - Integrals assume temperature-dependent material models are defined for the entire range [t_he_peak, t_max].
    """
    nu_irr_cu = copper_irradiation_resistivity(fluence)

    n_quad: Final[int] = 30
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)

    ihe = 0.0
    icu = 0.0
    isc = 0.0

    for xi, wi in zip(nodes, weights, strict=False):
        ti = 0.5 * (xi + 1.0) * (t_max - t_he_peak) + t_he_peak
        dti = 0.5 * wi * (t_max - t_he_peak)

        nu_cu = copper_electrical_resistivity(ti, field, rrr) + nu_irr_cu

        ihe += dti * helium_specific_heat_capacity(ti) * helium_density(ti) / nu_cu
        icu += dti * copper_specific_heat_capacity(ti) * copper_density(ti) / nu_cu
        isc += dti * nb3sn_specific_heat_capacity(ti) * nb3sn_density(ti) / nu_cu

    return ihe, icu, isc


def calculate_quench_protection_current_density(
    tau_discharge: float,
    peak_field: float,
    f_cu: float,
    f_he: float,
    t_he_peak: float,
    t_max: float,
    cu_rrr: float,
    detection_time: float,
    fluence: float,
) -> float:
    """
    Calculates the current density limited by the protection limit.

    Simplified 0-D adiabatic heat balance "hotspot criterion" model.

    Uses temperature- and magnetic field-dependent material properties to determine the
    maximum allowable winding pack current density during a quench. Assumes that the magnetic
    field does not decay over time and accounts for contributions from copper, helium, and
    superconductor materials using temperature integrals weighted by material properties.

    :author M. Coleman, UKAEA

    :param float tau_discharge: Quench discharge time constant [s].
    :param float peak_field: Magnetic field at the peak point [T].
    :param float f_cu: Fraction of conductor cross-section that is copper.
    :param float f_he: Fraction of cable occupied by helium.
    :param float t_he_peak: Peak helium temperature at quench initiation [K].
    :param float t_max: Maximum allowed conductor temperature during quench [K].
    :param float cu_rrr: Residual resistivity ratio of copper.
    :param float detection_time: Detection time delay [s].
    :param float fluence: Neutron fluence [n/m²].
    :returns: Maximum allowable winding pack current density [A/m²].

    :notes:
        - Assumes constant magnetic field over the duration of the quench.
        - Assumes the dump resistor has a constant resistance much higher
        than that of the TF coil.
        - Operates on the current-carring cross-section of a conductor. The
        jacket and insulation are ignored. The actual allowable WP current
        density must be weighted with the ratio of current-carrying cross-section
        vs. total WP cross-section (including jacket and insulation).
    """

    i_he, i_cu, i_sc = _quench_integrals(t_he_peak, t_max, peak_field, cu_rrr, fluence)

    f_cond = 1.0 - f_he  # Fraction of the cable XS area that is not helium
    f_cu_cable = f_cond * f_cu  # Fraction of the cable XS that is copper
    f_sc_cable = f_cond * (
        1.0 - f_cu
    )  # Fraction of the cable XS that is superconductor

    factor = 1.0 / (0.5 * tau_discharge + detection_time)
    total_integral = f_he * i_he + f_cu_cable * i_cu + f_sc_cable * i_sc

    return np.sqrt(factor * f_cu_cable * total_integral)
