"""Module containing neoclassical computations
author: J Lion, IPP Greifswald
Formulas used are described in:
Beidler (2013), https://doi.org/10.1088/0029-5515/51/7/076001
"""

import numpy as np

NO_ROOTS = 30
"""Number of Gauss laguerre roots"""

densities: list[float] = None
"""Densities of the species that are considered [/m3]"""

temperatures: list[float] = None
"""Temperature of the species that are considered [J]"""

dr_densities: list[float] = None
"""Radial derivative of the density of the species [/m3]"""

dr_temperatures: list[float] = None
"""Radial derivative of the temperature of the species [J]"""

roots: list[float] = None
"""Gauss Laguerre Roots"""

weights: list[float] = None
"""Gauss Laguerre Weights"""

nu: list[float] = None
"""90-degree deflection frequency on GL roots"""

nu_star: list[float] = None
"""Dimensionless deflection frequency"""

nu_star_averaged: list[float] = None
"""Maxwellian averaged dimensionless 90-degree deflection frequency for electrons (index 1) and ions (index 2)"""

vd: list[float] = None
"""Drift velocity on GL roots"""

kt: list[float] = None
"""Thermal energy on GL roots"""

er: float = None
"""Radial electrical field [V/m]"""

iota: float = None
"""Iota (1/safety factor)"""

d11_mono: list[float] = None
"""Radial monoenergetic transport coefficient on GL roots (species dependent)"""

d11_plateau: list[float] = None
"""Toroidal monoenergetic transport coefficient as given by the stellarator
input json file as function of nu_star, normalized by the banana value.
"""

d111: list[float] = None
"""Radial integrated transport coefficient (n=1) (species dependent)"""

d112: list[float] = None
"""Radial integrated transport coefficient (n=2) (species dependent)"""

d113: list[float] = None
"""Radial integrated transport coefficient (n=3) (species dependent)"""

q_flux: list[float] = None
"""energy transport flux (J/m2)"""

gamma_flux: list[float] = None
"""energy flux from particle transport"""

d31_mono: list[float] = None
"""Toroidal monoenergetic transport coefficient"""

eps_eff: float = None
"""Epsilon effective (used in neoclassics_calc_D11_mono)"""

r_eff: float = None


def init_neoclassics_variables():
    global \
        densities, \
        temperatures, \
        dr_densities, \
        dr_temperatures, \
        roots, \
        weights, \
        nu, \
        nu_star, \
        nu_star_averaged, \
        vd, \
        kt, \
        er, \
        iota, \
        d11_mono, \
        d11_plateau, \
        d111, \
        d112, \
        d113, \
        q_flux, \
        gamma_flux, \
        d31_mono, \
        eps_eff, \
        r_eff

    densities = np.zeros(4)
    temperatures = np.zeros(4)
    dr_densities = np.zeros(4)
    dr_temperatures = np.zeros(4)
    roots = np.zeros(NO_ROOTS)
    weights = np.zeros(NO_ROOTS)
    nu = np.zeros((4, NO_ROOTS))
    nu_star = np.zeros((4, NO_ROOTS))
    nu_star_averaged = np.zeros(4)
    vd = np.zeros((4, NO_ROOTS))
    kt = np.zeros((4, NO_ROOTS))
    iota = 1.0
    d11_mono = np.zeros((4, NO_ROOTS))
    d11_plateau = np.zeros((4, NO_ROOTS))
    d111 = np.zeros(4)
    d112 = np.zeros(4)
    d113 = np.zeros(4)
    q_flux = np.zeros(4)
    gamma_flux = np.zeros(4)
    d31_mono = np.zeros(NO_ROOTS)
    eps_eff = 1e-5
    r_eff = 0.0
    er = 0.0
