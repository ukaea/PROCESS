"""Module containing variables for the impurity radiation models"""

from dataclasses import dataclass, field
from enum import IntEnum, unique
from types import DynamicClassAttribute

import numpy as np

from process.core import constants

N_IMPURITIY_LOSS_FUNCTION_POINTS = 600
"""Number of points in the impurity loss function (L_z) tables"""


@unique
class IonSpecies(IntEnum):
    """Enumeration of ion species for impurity radiation model"""

    HYDROGEN = (1, 1, "H_", constants.M_PROTIUM_AMU)
    HELIUM = (2, 2, "He", constants.M_HELIUM_AMU)
    BERYLLIUM = (3, 4, "Be", constants.M_BERYLLIUM_AMU)
    CARBON = (4, 6, "C_", constants.M_CARBON_AMU)
    NITROGEN = (5, 7, "N_", constants.M_NITROGEN_AMU)
    OXYGEN = (6, 8, "O_", constants.M_OXYGEN_AMU)
    NEON = (7, 10, "Ne", constants.M_NEON_AMU)
    SILICON = (8, 14, "Si", constants.M_SILICON_AMU)
    ARGON = (9, 18, "Ar", constants.M_ARGON_AMU)
    IRON = (10, 26, "Fe", constants.M_IRON_AMU)
    NICKEL = (11, 28, "Ni", constants.M_NICKEL_AMU)
    KRYPTON = (12, 36, "Kr", constants.M_KRYPTON_AMU)
    XENON = (13, 54, "Xe", constants.M_XENON_AMU)
    TUNGSTEN = (14, 74, "W_", constants.M_TUNGSTEN_AMU)

    def __new__(cls, index: int, full_ion_charge: int, label: str, mass_amu: float):
        """Create a new instance of IonSpecies"""
        obj = int.__new__(cls, index)
        obj._value_ = index
        obj._full_ion_charge_ = full_ion_charge
        obj._label_ = label
        obj._mass_amu_ = mass_amu
        return obj

    @DynamicClassAttribute
    def full_ion_charge(self):
        """The full ion charge for this ion species."""
        return self._full_ion_charge_

    @DynamicClassAttribute
    def label(self):
        """The label for this ion species."""
        return self._label_

    @DynamicClassAttribute
    def mass_amu(self):
        """The atomic mass of this ion species in Atomic Mass Units (amu)."""
        return self._mass_amu_


N_IMPURITIES = len(IonSpecies)
"""Number of ion species in impurity radiation model
"""


@dataclass(slots=True)
class ImpurityRadiationData:
    """Dataclass holding impurity radiation variables"""

    radius_plasma_core_norm: float = 0.6
    """Normalised radius defining the 'core' region"""

    f_p_plasma_core_rad_reduction: float = 1.0
    """Fraction of radiation from 'core' region"""

    f_nd_impurity_electrons: list[float] = field(
        default_factory=lambda: np.array([
            1.0,
            0.1,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ])
    )

    n_charge_impurity_profile: list[float] = field(
        default_factory=lambda: np.zeros((
            N_IMPURITIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )
    """charge profile of impurities"""

    imp_label: list[str] = field(
        default_factory=lambda: np.array([
            "H_",
            "He",
            "Be",
            "C_",
            "N_",
            "O_",
            "Ne",
            "Si",
            "Ar",
            "Fe",
            "Ni",
            "Kr",
            "Xe",
            "W_",
        ])
    )
    """Labels for impurity species"""

    imp_full_ion_charge: list[int] = field(
        default_factory=lambda: np.array([
            1,
            2,
            4,
            6,
            7,
            8,
            10,
            14,
            18,
            26,
            28,
            36,
            54,
            74,
        ])
    )
    """Full ion charge of impurity species"""

    impurity_arr_label: list[str] = field(
        default_factory=lambda: np.full(N_IMPURITIES, "  ")
    )
    """Array of impurity labels for internal use"""

    impurity_arr_z: list[float] = field(default_factory=lambda: np.zeros(N_IMPURITIES))

    m_impurity_amu_array: list[float] = field(
        default_factory=lambda: np.zeros(N_IMPURITIES)
    )
    """2D array of impurity atomic masses in Atomic Mass Units (amu)"""

    f_nd_impurity_electron_array: list[float] = field(
        default_factory=lambda: np.zeros(N_IMPURITIES)
    )
    """2D array of impurity relative densities (n_imp/n_e)"""

    impurity_arr_len_tab: list[int] = field(
        default_factory=lambda: np.full(N_IMPURITIES, 0)
    )

    temp_impurity_keV_array: list[float] = field(
        default_factory=lambda: np.zeros((
            N_IMPURITIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )
    """2D array of impurity temperatures in kilo-electronvolts (keV)"""

    pden_impurity_lz_nd_temp_array: list[float] = field(
        default_factory=lambda: np.zeros((
            N_IMPURITIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )

    impurity_arr_zav: list[float] = field(
        default_factory=lambda: np.zeros((
            N_IMPURITIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )


CREATE_DICTS_FROM_DATACLASS = ImpurityRadiationData
