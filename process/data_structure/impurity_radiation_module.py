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

radius_plasma_core_norm: float = None
"""Normalised radius defining the 'core' region"""

coreradiationfraction: float = None
"""Fraction of radiation from 'core' region that is subtracted from the loss power"""

f_nd_impurity_electrons: list[float] = None

imp_label: list[str] = None

impurity_arr_label: list[str] = None

impurity_arr_z: list[float] = None

m_impurity_amu_array: list[float] = None
"""2D array of impurity atomic masses in Atomic Mass Units (amu)"""

f_nd_impurity_electron_array: list[float] = None
"""2D array of impurity relative densities (n_imp/n_e)"""

impurity_arr_len_tab: list[int] = None

temp_impurity_keV_array: list[float] = None
"""2D array of impurity temperatures in kilo-electronvolts (keV)"""

impurity_arr_lz_wm3: list[float] = None

impurity_arr_zav: list[float] = None


def init_impurity_radiation_module():
    global radius_plasma_core_norm
    global coreradiationfraction
    global f_nd_impurity_electrons
    global imp_label
    global impurity_arr_label
    global impurity_arr_z
    global m_impurity_amu_array
    global f_nd_impurity_electron_array
    global impurity_arr_len_tab
    global temp_impurity_keV_array
    global impurity_arr_lz_wm3
    global impurity_arr_zav

    radius_plasma_core_norm = 0.6
    coreradiationfraction = 1.0
    f_nd_impurity_electrons = np.array([
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
    imp_label = np.array([
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
    impurity_arr_label = np.full(N_IMPURITIES, "  ")
    impurity_arr_z = np.zeros(N_IMPURITIES)
    m_impurity_amu_array = np.zeros(N_IMPURITIES)
    f_nd_impurity_electron_array = np.zeros(N_IMPURITIES)
    impurity_arr_len_tab = np.full(N_IMPURITIES, 0)
    temp_impurity_keV_array = np.zeros((N_IMPURITIES, 200))
    impurity_arr_lz_wm3 = np.zeros((N_IMPURITIES, 200))
    impurity_arr_zav = np.zeros((N_IMPURITIES, 200))
