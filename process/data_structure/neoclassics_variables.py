"""Module containing variables for the neoclassical models

Formulas used are described in:
Beidler (2013), https://doi.org/10.1088/0029-5515/51/7/076001
"""

from dataclasses import dataclass, field

import numpy as np

from process.core.data_structure.parameter import Parameter, PROCESSModelData

NO_ROOTS = 30
"""Number of Gauss laguerre roots"""


@dataclass(slots=True)
class NeoclassicsData(PROCESSModelData):
    """Dataclass holding neoclassics variables"""

    densities: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Densities of the species that are considered [/m3]"""

    temperatures: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Temperature of the species that are considered [J]"""

    dr_densities: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Radial derivative of the density of the species [/m3]"""

    dr_temperatures: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Radial derivative of the temperature of the species [J]"""

    roots: Parameter[list[float]] = field(default_factory=lambda: np.zeros(NO_ROOTS))
    """Gauss Laguerre Roots"""

    weights: Parameter[list[float]] = field(default_factory=lambda: np.zeros(NO_ROOTS))
    """Gauss Laguerre Weights"""

    nu: Parameter[list[float]] = field(default_factory=lambda: np.zeros((4, NO_ROOTS)))
    """90-degree deflection frequency on GL roots"""

    nu_star: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((4, NO_ROOTS))
    )
    """Dimensionless deflection frequency"""

    nu_star_averaged: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Maxwellian averaged dimensionless 90-degree deflection frequency for electrons (index 1) and ions (index 2)"""

    vd: Parameter[list[float]] = field(default_factory=lambda: np.zeros((4, NO_ROOTS)))
    """Drift velocity on GL roots"""

    kt: Parameter[list[float]] = field(default_factory=lambda: np.zeros((4, NO_ROOTS)))
    """Thermal energy on GL roots"""

    er: Parameter[float] = 0.0
    """Radial electrical field [V/m]"""

    iota: Parameter[float] = 1.0
    """Iota (1/safety factor)"""

    d11_mono: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((4, NO_ROOTS))
    )
    """Radial monoenergetic transport coefficient on GL roots (species dependent)"""

    d11_plateau: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((4, NO_ROOTS))
    )
    """Toroidal monoenergetic transport coefficient as given by the stellarator
    input json file as function of nu_star, normalised by the banana value.
    """

    d111: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Radial integrated transport coefficient (n=1) (species dependent)"""

    d112: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Radial integrated transport coefficient (n=2) (species dependent)"""

    d113: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """Radial integrated transport coefficient (n=3) (species dependent)"""

    q_flux: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """energy transport flux (J/m2)"""

    gamma_flux: Parameter[list[float]] = field(default_factory=lambda: np.zeros(4))
    """energy flux from particle transport"""

    d31_mono: Parameter[list[float]] = field(default_factory=lambda: np.zeros(NO_ROOTS))
    """Toroidal monoenergetic transport coefficient"""

    eps_eff: Parameter[float] = 1e-5
    """Epsilon effective (used in neoclassics_calc_D11_mono)"""

    r_eff: Parameter[float] = 0.0


CREATE_DICTS_FROM_DATACLASS = NeoclassicsData
