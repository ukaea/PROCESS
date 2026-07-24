from dataclasses import dataclass, field

import numpy as np

N_IMPURITIES = 14
"""Number of ion species in impurity radiation model

( 1)  Hydrogen  (fraction calculated by code)
( 2)  Helium
( 3)  Beryllium
( 4)  Carbon
( 5)  Nitrogen
( 6)  Oxygen
( 7)  Neon
( 8)  Silicon
( 9)  Argon
(10)  Iron
(11)  Nickel
(12)  Krypton
(13)  Xenon
(14)  Tungsten
"""


@dataclass(slots=True)
class ImpurityRadiationData:
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
        default_factory=lambda: np.zeros((N_IMPURITIES, 200))
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

    impurity_arr_label: list[str] = field(
        default_factory=lambda: np.full(N_IMPURITIES, "  ")
    )

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
        default_factory=lambda: np.zeros((N_IMPURITIES, 200))
    )
    """2D array of impurity temperatures in kilo-electronvolts (keV)"""

    pden_impurity_lz_nd_temp_array: list[float] = field(
        default_factory=lambda: np.zeros((N_IMPURITIES, 200))
    )

    impurity_arr_zav: list[float] = field(
        default_factory=lambda: np.zeros((N_IMPURITIES, 200))
    )


CREATE_DICTS_FROM_DATACLASS = ImpurityRadiationData
