import numpy as np
from typing import Final, Tuple


# Material property parameterisations

def copper_density(temperature: float) -> float:
    """
    Calculate the density of cryogenic copper [kg/m^3].
    """
    return 8960.0

def copper_specific_heat_capacity(temperature: float) -> float:
    """
    Calculate the specific heat capacity of cryogenic copper [J/(kg·K)] at a given temperature.

    References:
        - L. Dresner, “Stability of Superconductors”, Plenum Press, NY, 1995

    Parameters:
        temperature (float): Temperature [K].

    Returns:
        float: Specific heat capacity [J/(kg·K)].
    """
    GAMMA: Final[float] = 0.011       # [J/K²/kg]
    BETA: Final[float] = 0.0011       # [J/K⁴/kg]
    CP_300: Final[float] = 385.491    # Room-temperature specific heat [J/K/kg]

    cp_low = BETA * temperature**3 + GAMMA * temperature
    cp = 1.0 / (1.0 / CP_300 + 1.0 / cp_low)

    return cp

def copper_electrical_resistivity(temperature: float, field: float, rrr: float) -> float:
    """
    Calculate the electrical resistivity of cryogenic copper with temperature and field dependence [Ω·m].

    References:
        - NIST MONOGRAPH 177: J. Simon, E. S. Drexler, and R. P. Reed,
          "Properties of Copper and Copper Alloys at Cryogenic Temperatures",
          U.S. Government Printing Office, February 1992.
          draft version from 1987, implemented as found in:
        - EFDA Material Data Compilation for Superconductor Simulation:
          P. Bauer, H. Rajainmaki, E. Salpietro, EFDA CSU, Garching, 04/18/07.

    Parameters:
        temperature (float): Operating temperature [K].
        field (float): Operating magnetic field [T].
        rrr (float): Residual resistivity ratio (dimensionless).

    Returns:
        float: The electrical resistivity [Ω·m].
    """
    # Constants from EFDA documentation (Page 2-17)
    P1: Final[float] = 1.171e-17
    P2: Final[float] = 4.49
    P3: Final[float] = 4.5e-7
    P4: Final[float] = 3.35
    P5: Final[float] = 50.0
    P6: Final[float] = 6.428
    P7: Final[float] = 1.69e-8
    P8: Final[float] = 0.4531
    P9: Final[float] = 1.553e-8
    POLY_COEFFS: Final[list[float]] = [-2.662, 0.3168, 0.6229, -0.1839, 0.01827]

    # TODO: Implement PROCESS standard for kludging / warning
    t = np.clip(temperature, 4.0)
    # TODO: this function at present has the potential to return poor values at low 
    # field... 
    rrr = np.clip(rrr, 1.0)

    t = temperature

    # Compute rho1
    numerator = P1 * t**P2
    denominator = 1.0 + P3 * t**P4 * np.exp(-(P5 / t)**P6)
    rho1 = numerator / denominator

    # Compute rho2
    rho2 = P7 / rrr + rho1 + P8 * (P7 * rho1) / (rrr * rho1 + P7)

    # Compute magnetic field correction factor
    x = np.log10(P9 * field / rho2)
    a = sum(c * x**i for i, c in enumerate(POLY_COEFFS))

    # Final resistivity
    return rho2 * (1.0 + 10**a)


def copper_irradiation_resistivity(fluence: float) -> float:
    """
    Calculate the radiation-induced electrical resistivity of copper [Ω·m].

    References:
        - M. Kovari, 09/11/2012, internal notes (Excel / Mathcad)
          K:\Technology Program\3PT\WP12\PEX\Super-Xdivertor for DEMO\Coil calcs 7
        - M. Nakagawa et al., "High-dose neutron-irradiation effects in fcc metals at 4.6 K",
          Phys. Rev. B 16, 5285 (1977)
          https://doi.org/10.1103/PhysRevB.16.5285

    Parameters:
        fluence (float): Total neutron fluence [n/m²].

    Returns:
        float: Radiation-induced resistivity [Ω·m].
    """
    C1: Final[float] = 0.00283
    C2: Final[float] = -0.0711
    C3: Final[float] = 0.77982
    RES_SCALE: Final[float] = 1e-9
    FLU_SCALE: Final[float] = 1e-22

    fluence_norm = FLU_SCALE * fluence

    return RES_SCALE * (C1 * fluence_norm**3 + C2 * fluence_norm**2 + C3 * fluence_norm)


def nb3sn_density(temperature: float) -> float:
    """
    Calculate the density of Nb3Sn [kg/m^3].
    """
    return 8040.0


def nb3sn_specific_heat_capacity(temperature: float) -> float:
    """
    Calculate the specific heat capacity of Nb3Sn [J/(kg·K)] 
    as a function of temperature.

    References:
        EFDA Material Data Compilation for Superconductor Simulation,
        P. Bauer, H. Rajainmaki, E. Salpietro,
        EFDA CSU, Garching, 04/18/07.
    
    Notes:
        Ignoring the superconducting part for quench

    Parameters:
        temperature (float): Temperature [K].

    Returns:
        float: Specific heat capacity [J/(kg·K)].
    """
    GAMMA: Final[float] = 0.1       # [J/K²/kg] - Electronic contribution
    BETA: Final[float] = 0.001      # [J/K⁴/kg] - Lattice (Debye) contribution
    CP_HIGH: Final[float] = 210.0   # [J/K/kg] - High temperature limit

    cp_low = BETA * temperature**3 + GAMMA * temperature
    cp = 1.0 / (1.0 / CP_HIGH + 1.0 / cp_low)

    return cp

def helium_density(temperature: float) -> float:
    return 150.0


def helium_specific_heat_capacity(temperature: float) -> float:
    return 4750.0


# Quench model 

def quench_integrals(
    t_he_peak: float,
    t_max: float,
    field: float,
    rrr: float,
    fluence: float
) -> Tuple[float, float, float]:
    """
    Calculates the material property integrals for quench protection.

    Equation:
        I_He = ∫ [ρ_He(T) · c_He(T)] / η_Cu(T, B, RRR) dT
        I_Cu = ∫ [ρ_Cu(T) · c_Cu(T)] / η_Cu(T, B, RRR) dT
        I_sc = ∫ [ρ_sc(T) · c_sc(T)] / η_Cu(T, B, RRR) dT

    Parameters:
        t_he_peak (float): Lower temperature bound of integration [K].
        t_max (float): Upper temperature bound of integration [K].
        field (float): Magnetic field [T].
        rrr (float): Residual resistivity ratio of copper.
        fluence (float): Neutron fluence [n/cm²] (for irradiation effects).

    Returns:
        Tuple[float, float, float]: (I_He, I_Cu, I_sc)
    """
    nu_irr_cu = copper_irradiation_resistivity(fluence)

    N_QUAD: Final[int] = 30
    nodes, weights = np.polynomial.legendre.leggauss(N_QUAD)

    ihe = 0.0
    icu = 0.0
    isc = 0.0

    for xi, wi in zip(nodes, weights):
        ti = 0.5 * (xi + 1.0) * (t_max - t_he_peak) + t_he_peak
        dti = 0.5 * wi * (t_max - t_he_peak)

        nu_cu = copper_electrical_resistivity(ti, field, rrr) + nu_irr_cu

        ihe += dti * helium_specific_heat_capacity(ti) * helium_density(ti) / nu_cu
        icu += dti * copper_specific_heat_capacity(ti) * copper_density(ti) / nu_cu
        isc += dti * nb3sn_specific_heat_capacity(ti) * nb3sn_density(ti) / nu_cu

    return ihe, icu, isc


def calculate_quench_protection_current_density_magnetoresistive(
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
    Calculates the current density limited by the protection limit, using temperature
    and field-dependent material properties. The field is assumed not to decay with time.

    Equation:
        J = sqrt[
            (1 / (τ_discharge / 2 + t_detection)) * f_Cu *
            (f_He * I_He + f_Cu * I_Cu + f_sc * I_sc)
        ]

        where:

        I_He = ∫[T₀ to T_max] (ρ_He(T) * c_P,He(T)) / η_Cu(T, B, RRR) dT  
        I_Cu = ∫[T₀ to T_max] (ρ_Cu(T) * c_P,Cu(T)) / η_Cu(T, B, RRR) dT  
        I_sc = ∫[T₀ to T_max] (ρ_sc(T) * c_P,sc(T)) / η_Cu(T, B, RRR) dT

    Parameters:
        tau_discharge: Quench discharge time constant [s]
        peak_field: Magnetic field at the peak point [T]
        f_cu: Fraction of conductor cross-section that is copper
        f_he: Fraction of cable occupied by helium
        t_he_peak: Peak helium temperature at quench initiation [K]
        t_max: Maximum allowed conductor temperature during quench [K]
        cu_rrr: Residual resistivity ratio of copper
        detection_time: Detection time delay [s]
        fluence: Neutron fluence [n/m^2]

    Returns:
        float: Maximum allowable winding pack current density [A/m²]
    """
    i_he, i_cu, i_sc = quench_integrals(t_he_peak, t_max, peak_field, cu_rrr, fluence)

    f_cond = 1.0 - f_he                 # Fraction of the cable XS area that is not helium
    f_cu_cable = f_cond * f_cu          # Fraction of the cable XS that is copper
    f_sc_cable = f_cond * (1.0 - f_cu)  # Fraction of the cable XS that is superconductor

    factor = 1.0 / (0.5 * tau_discharge + detection_time)
    total_integral = f_he * i_he + f_cu_cable * i_cu + f_sc_cable * i_sc

    return np.sqrt(factor * f_cu_cable * total_integral)
