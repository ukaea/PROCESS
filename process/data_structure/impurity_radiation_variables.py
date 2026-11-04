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

    HYDROGEN = (1, 1, "H_", constants.M_PROTIUM_AMU, "hydrogen")
    HELIUM = (2, 2, "He", constants.M_HELIUM_AMU, "helium")
    BERYLLIUM = (3, 4, "Be", constants.M_BERYLLIUM_AMU, "beryllium")
    CARBON = (4, 6, "C_", constants.M_CARBON_AMU, "carbon")
    NITROGEN = (5, 7, "N_", constants.M_NITROGEN_AMU, "nitrogen")
    OXYGEN = (6, 8, "O_", constants.M_OXYGEN_AMU, "oxygen")
    NEON = (7, 10, "Ne", constants.M_NEON_AMU, "neon")
    SILICON = (8, 14, "Si", constants.M_SILICON_AMU, "silicon")
    ARGON = (9, 18, "Ar", constants.M_ARGON_AMU, "argon")
    IRON = (10, 26, "Fe", constants.M_IRON_AMU, "iron")
    NICKEL = (11, 28, "Ni", constants.M_NICKEL_AMU, "nickel")
    KRYPTON = (12, 36, "Kr", constants.M_KRYPTON_AMU, "krypton")
    XENON = (13, 54, "Xe", constants.M_XENON_AMU, "xenon")
    TUNGSTEN = (14, 74, "W_", constants.M_TUNGSTEN_AMU, "tungsten")

    def __new__(
        cls,
        value: int,
        full_ion_charge: int,
        label: str,
        mass_amu: float,
        full_name: str,
    ):
        """Create a new instance of IonSpecies"""
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_ion_charge_ = full_ion_charge
        obj._label_ = label
        obj._mass_amu_ = mass_amu
        obj._full_name_ = full_name
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

    @DynamicClassAttribute
    def full_name(self):
        """The full name of this ion species."""
        return self._full_name_

    @DynamicClassAttribute
    def is_impurity(self) -> bool:
        """Whether this ion species is an impurity."""
        return self not in {IonSpecies.HYDROGEN, IonSpecies.HELIUM}


N_ION_SPECIES = len(IonSpecies)
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
            N_ION_SPECIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )
    """charge profile of impurities"""

    impurity_arr_label: list[str] = field(
        default_factory=lambda: np.full(N_ION_SPECIES, "  ")
    )
    """Array of impurity labels for internal use"""

    impurity_arr_z: list[float] = field(default_factory=lambda: np.zeros(N_ION_SPECIES))

    m_impurity_amu_array: list[float] = field(
        default_factory=lambda: np.zeros(N_ION_SPECIES)
    )
    """2D array of impurity atomic masses in Atomic Mass Units (amu)"""

    f_nd_impurity_electron_array: list[float] = field(
        default_factory=lambda: np.zeros(N_ION_SPECIES)
    )
    """2D array of impurity relative densities (n_imp/n_e)"""

    impurity_arr_len_tab: list[int] = field(
        default_factory=lambda: np.full(N_ION_SPECIES, 0)
    )

    temp_impurity_keV_array: list[float] = field(
        default_factory=lambda: np.zeros((
            N_ION_SPECIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )
    """2D array of impurity temperatures in kilo-electronvolts (keV)"""

    pden_impurity_lz_nd_temp_array: list[float] = field(
        default_factory=lambda: np.zeros((
            N_ION_SPECIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )

    impurity_arr_zav: list[float] = field(
        default_factory=lambda: np.zeros((
            N_ION_SPECIES,
            N_IMPURITIY_LOSS_FUNCTION_POINTS,
        ))
    )


CREATE_DICTS_FROM_DATACLASS = ImpurityRadiationData
