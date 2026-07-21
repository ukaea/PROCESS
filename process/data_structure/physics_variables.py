"""Module containing tokamak plasma physics routines

This module contains all the primary plasma physics routines
for a tokamak device.


Module containing global variables relating to the plasma physics
"""

from dataclasses import dataclass, field
from enum import IntEnum

import numpy as np


class PlasmaIgnitionModel(IntEnum):
    """Enum for plasma ignition models."""

    NON_IGNITED = 0
    IGNITED = 1


class DivertorNumberModels(IntEnum):
    """Enum for divertor number models. `i_single_null` is the index for this enum."""

    DOUBLE_NULL = 0
    SINGLE_NULL = 1


class ConfinementMode(IntEnum):
    """Enum for plasma confinement mode"""

    L_MODE = (0, "L")
    H_MODE = (1, "H")
    I_MODE = (2, "I")
    STELLARATOR = (3, "Stell")
    OHMIC = (4, "Ohmic")

    def __new__(cls, value: int, abbreviation: str):
        """Create a new instance of ConfinementMode.

        Parameters
        ----------
        value : int
            The enum value
        abbreviation : str
            The abbreviation of the confinement mode

        Returns
        -------
        ConfinementMode
            A new enum instance with the given value and abbreviation
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.abbreviation = abbreviation

        return obj


class ConfinementTimeModel(IntEnum):
    """Confinement time (τ_E) model types"""

    USER_INPUT = (0, "User input electron confinement   ", None)
    NEO_ALCATOR = (
        1,
        f"Neo-Alcator                ({ConfinementMode.OHMIC.abbreviation})",
        ConfinementMode.OHMIC,
    )
    MIRNOV = (
        2,
        f"Mirnov                         ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    MEREZHKIN_MUHKOVATOV = (
        3,
        f"Merezkhin-Muhkovatov    ({ConfinementMode.OHMIC.abbreviation})({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.OHMIC | ConfinementMode.L_MODE,
    )
    SHIMOMURA = (
        4,
        f"Shimomura                      ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    KAYE_GOLDSTON = (
        5,
        f"Kaye-Goldston                  ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    ITER_89P = (
        6,
        f"ITER 89-P                      ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    ITER_89_0 = (
        7,
        f"ITER 89-O                      ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    REBUT_LALLIA = (
        8,
        f"Rebut-Lallia                   ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    GOLDSTON = (
        9,
        f"Goldston                       ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    T_10 = (
        10,
        f"T10                            ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    JAERI = (
        11,
        f"JAERI / Odajima-Shimomura      ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    KAYE_BIG = (
        12,
        f"Kaye-Big Complex               ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    ITER_H90_P = (
        13,
        f"ITER H90-P                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    MINIMUM_OF_ITER_89P_AND_ITER_89_0 = (
        14,
        f"ITER 89-P & 89-O min           ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    RIEDEL_L = (
        15,
        f"Riedel                         ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    CHRISTIANSEN = (
        16,
        f"Christiansen                   ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    LACKNER_GOTTARDI = (
        17,
        f"Lackner-Gottardi               ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    NEO_KAYE = (
        18,
        f"Neo-Kaye                       ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    RIEDEL_H = (
        19,
        f"Riedel                         ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_H90_P_AMENDED = (
        20,
        f"ITER H90-P amended             ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    SUDO_ET_AL = (
        21,
        f"LHD                        ({ConfinementMode.STELLARATOR.abbreviation})",
        ConfinementMode.STELLARATOR,
    )
    GYRO_REDUCED_BOHM = (
        22,
        f"Gyro-reduced Bohm          ({ConfinementMode.STELLARATOR.abbreviation})",
        ConfinementMode.STELLARATOR,
    )
    LACKNER_GOTTARDI_STELLARATOR = (
        23,
        f"Lackner-Gottardi           ({ConfinementMode.STELLARATOR.abbreviation})",
        ConfinementMode.STELLARATOR,
    )
    ITER_93H = (
        24,
        f"ITER-93H  ELM-free             ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    TITAN_REMOVED = (
        25,
        "TITAN RFP OBSOLETE                (N/A)",
        None,
    )
    ITER_H97P = (
        26,
        f"ITER H-97P ELM-free            ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_H97P_ELMY = (
        27,
        f"ITER H-97P ELMy                ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_96P = (
        28,
        f"ITER-96P (ITER-97L)            ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    VALOVIC_ELMY = (
        29,
        f"Valovic modified ELMy          ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    KAYE = (
        30,
        f"Kaye 98 modified               ({ConfinementMode.L_MODE.abbreviation})",
        ConfinementMode.L_MODE,
    )
    ITER_PB98P_Y = (
        31,
        f"ITERH-PB98P(y)                 ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    IPB98_Y = (
        32,
        f"IPB98(y)                       ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_IPB98Y1 = (
        33,
        f"IPB98(y,1)                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_IPB98Y2 = (
        34,
        f"IPB98(y,2)                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_IPB98Y3 = (
        35,
        f"IPB98(y,3)                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITER_IPB98Y4 = (
        36,
        f"IPB98(y,4)                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ISS95_STELLARATOR = (
        37,
        f"ISS95                      ({ConfinementMode.STELLARATOR.abbreviation})",
        ConfinementMode.STELLARATOR,
    )
    ISS04_STELLARATOR = (
        38,
        f"ISS04                      ({ConfinementMode.STELLARATOR.abbreviation})",
        ConfinementMode.STELLARATOR,
    )
    DS03 = (
        39,
        f"DS03 beta-independent          ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    MURARI = (
        40,
        f'Murari "Non-power law"         ({ConfinementMode.H_MODE.abbreviation})',
        ConfinementMode.H_MODE,
    )
    PETTY08 = (
        41,
        f"Petty 2008                     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    LANG_HIGH_DENSITY = (
        42,
        f"Lang high density              ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    HUBBARD_NOMINAL = (
        43,
        f"Hubbard 2017 - nominal         ({ConfinementMode.I_MODE.abbreviation})",
        ConfinementMode.I_MODE,
    )
    HUBBARD_LOWER = (
        44,
        f"Hubbard 2017 - lower           ({ConfinementMode.I_MODE.abbreviation})",
        ConfinementMode.I_MODE,
    )
    HUBBARD_UPPER = (
        45,
        f"Hubbard 2017 - upper           ({ConfinementMode.I_MODE.abbreviation})",
        ConfinementMode.I_MODE,
    )
    MENARD_NSTX = (
        46,
        f"Menard NSTX                    ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    MENARD_NSTX_PETTY08_HYBRID = (
        47,
        f"Menard NSTX-Petty08 hybrid     ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    NSTX_GYRO_BOHM = (
        48,
        f"Buxton NSTX gyro-Bohm          ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITPA20 = (
        49,
        f"ITPA20                         ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )
    ITPA20_IL = (
        50,
        f"ITPA20-IL                      ({ConfinementMode.H_MODE.abbreviation})",
        ConfinementMode.H_MODE,
    )

    def __new__(cls, value: int, full_name: str, mode: ConfinementMode = None):
        """Create a new instance of ConfinementTimeModel.

        Parameters
        ----------
        value : int
            The enum value
        full_name : str
            The full name of the confinement time model
        mode : ConfinementMode
            The confinement mode associated with the model

        Returns
        -------
        ConfinementTimeModel
            A new enum instance with the given value and full name
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.full_name = full_name
        obj.mode = mode
        return obj


class ConfinementRadiationLossModel(IntEnum):
    """Confinement radiation loss model types"""

    FULL_RADIATION = (0, "All radiation included in loss power term")
    CORE_ONLY = (1, "Only core radiation included in loss power term")
    NO_RADIATION = (2, "No radiation included in loss power term")

    def __new__(cls, value: int, description: str):
        """Create a new instance of ConfinementRadiationLossModel.

        Parameters
        ----------
        value : int
            The enum value
        description : str
            The description of the radiation loss model

        Returns
        -------
        ConfinementRadiationLossModel
            A new enum instance with the given value and description
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        return obj


N_CONFINEMENT_SCALINGS = len(ConfinementTimeModel)
"""number of energy confinement time scaling laws"""


@dataclass(slots=True)
class PhysicsData:
    iscz: int = 0

    err242: int = 0

    err243: int = 0

    f_p_plasma_separatrix_rad: float = 0.0
    """Separatrix radiation fraction"""

    e_plasma_beta: float = 0.0
    """[J]"""

    p_plasma_heating_total_mw: float = 0.0
    """Total heating power given to the plasma (Pₕₑₐₜ) [MW]"""

    t_energy_confinement_beta: float = 0.0
    """[s]"""

    ptarmw: float = 0.0

    lambdaio: float = 0.0

    drsep: float = 0.0

    fio: float = 0.0

    fli: float = 0.0

    flo: float = 0.0

    fui: float = 0.0

    fuo: float = 0.0

    plimw: float = 0.0

    plomw: float = 0.0

    puimw: float = 0.0

    puomw: float = 0.0

    rho_star: float = 0.0

    nu_star: float = 0.0

    beta_mcdonald: float = 0.0

    itart_r: float = 0.0

    # Var in subroutine plasma_composition which requires re-initialisation on
    # each new run:
    first_call: int = 1

    m_beam_amu: float = 0.0
    """beam ion mass (amu)"""

    m_fuel_amu: float = 0.0
    """average mass of fuel portion of ions (amu)"""

    m_ions_total_amu: float = 0.0
    """average mass of all ions (amu)"""

    m_plasma_fuel_ions: float = 0.0
    """Mass of the plasma fuel ions (kg)"""

    m_plasma_ions_total: float = 0.0
    """Mass of all ions in plasma (kg)"""

    m_plasma_alpha: float = 0.0
    """Mass of the alpha particles in the plasma (kg)"""

    m_plasma_electron: float = 0.0
    """Mass of the electrons in the plasma (kg)"""

    m_plasma: float = 0.0
    """Total mass of the plasma (kg)"""

    alphaj: float = 1.0
    """Plasma current profile index (⍺ⱼ)"""  # noqa: RUF001

    alphaj_wesson: float = None
    """Wesson-like current profile index"""

    alphan: float = 0.25
    """Plasma density profile index (⍺ₙ)"""  # noqa: RUF001

    alphap: float = 0.0
    """Plasma pressure profile index (⍺ₚ)"""  # noqa: RUF001

    fusden_alpha_total: float = 0.0
    """Alpha particle production rate per unit volume, from plasma and beams [particles/m³/sec]"""

    fusden_plasma_alpha: float = 0.0
    """Alpha particle production rate per unit volume, just from plasma [particles/m³/sec]"""

    alphat: float = 0.5
    """Plasma temperature profile index (⍺ₜ)"""  # noqa: RUF001

    aspect: float = 2.907
    """Plasma aspect ratio (A) (`iteration variable 1`)"""

    beamfus0: float = 1.0
    """multiplier for beam-background fusion calculation"""

    beta_total_vol_avg: float = 0.042
    """Volume averaged total plasma beta (⟨β⟩) (`iteration variable 5`) (calculated if stellarator)"""

    beta_fast_alpha: float = 0.0
    """Fast alpha beta component (β_alpha)"""

    beta_vol_avg_max: float = 0.0
    """Max allowable volume averaged beta (⟨β⟩<)"""

    beta_vol_avg_min: float = 0.0
    """Minimum allowable volume averaged beta (⟨β⟩>)"""

    beta_beam: float = 0.0
    """Neutral beam beta component (β_beam)"""

    beta_poloidal_vol_avg: float = 0.0
    """Volume averaged poloidal beta (⟨βₚ⟩)"""

    beta_poloidal_eps: float = 0.0
    """Poloidal beta and inverse aspcet ratio product (⟨βₚ⟩*ε)"""

    beta_toroidal_vol_avg: float = 0.0
    """Volume averaged toroidal beta (⟨βₜ⟩)"""

    beta_thermal_toroidal_profile: list[float] = field(default_factory=list)
    """Toroidal beta profile"""

    beta_thermal_vol_avg: float = 0.0
    """Volume averaged thermal beta (⟨βₜₕ⟩)"""

    beta_thermal_poloidal_vol_avg: float = 0.0
    """Volume averaged poloidal thermal beta (⟨βₚₜₕ⟩)"""

    beta_thermal_toroidal_vol_avg: float = 0.0
    """Volume averaged toroidal thermal beta (⟨βₜₕ⟩)"""

    beta_norm_total: float = 0.0
    """Normalised total beta (βₙ)"""

    beta_norm_thermal: float = 0.0
    """Normalised thermal beta (βₙₜₕ)"""

    beta_norm_toroidal: float = 0.0
    """Normalised toroidal beta (βₙₜ)"""

    beta_norm_poloidal: float = 0.0
    """Normalised poloidal beta (βₙₚ)"""

    e_plasma_beta_thermal: float = 0.0
    """Plasma thermal energy derived from thermal beta"""

    e_plasma_thermal_total: float = 0.0
    """Plasma total stored thermal energy (J)"""

    eden_plasma_thermal_vol_avg: float = 0.0
    """Plasma volume averaged thermal energy density (J/m³)"""

    e_plasma_electrons_thermal: float = 0.0
    """Plasma thermal energy in electrons (J)"""

    eden_plasma_electrons_thermal_vol_avg: float = 0.0
    """Plasma volume averaged thermal energy density in electrons (J/m³)"""

    e_plasma_ions_thermal: float = 0.0
    """Plasma thermal energy in ions (J)"""

    eden_plasma_ions_thermal_vol_avg: float = 0.0
    """Plasma volume averaged thermal energy density in ions (J/m³)"""

    betbm0: float = 1.5
    """leading coefficient for NB beta fraction"""

    b_plasma_surface_poloidal_average: float = 0.0
    """Plasma surface average poloidal field (T)"""

    b_plasma_toroidal_on_axis: float = 5.68
    """Plasma toroidal field on axis (Bᴛ(R₀)) [T] (`iteration variable 2`)"""

    b_plasma_inboard_toroidal: float = 0.0
    """Plasma inboard toroidal field (Bᴛ(R₀-a)) [T]"""

    b_plasma_outboard_toroidal: float = 0.0
    """Plasma outboard toroidal field (Bᴛ(R₀+a)) [T]"""

    b_plasma_toroidal_profile: list[float] = field(default_factory=list)
    """Plasma toroidal field profile (Bᴛ(r)) [T]"""

    b_plasma_total: float = 0.0
    """Sum of plasma total toroidal + poloidal field (Bₜₒₜ) [T]"""

    e_plasma_magnetic_stored: float = 0.0
    """Plasma stored magnetic energy [J]"""

    burnup: float = 0.0
    """fractional plasma burnup"""

    burnup_in: float = 0.0
    """fractional plasma burnup user input"""

    b_plasma_vertical_required: float = 0.0
    """Vertical field needed for plasma equilibrium (Bᵥ) [T]"""

    c_beta: float = 0.5
    """Destabalisation parameter for i_beta_norm_max=4 beta limit"""

    csawth: float = 1.0
    """coeff. for sawteeth effects on burn V-s requirement"""

    f_vol_plasma: float = 1.0
    """multiplying factor for the plasma volume (normally=1)"""

    f_r_conducting_wall: float = 1.35
    """maximum ratio of conducting wall distance to plasma minor radius for
    vertical stability (`constraint equation 23`)
    """

    nd_plasma_electrons_vol_avg: float = 9.8e19
    """Plasma volume averaged electron density (⟨nₑ⟩) [/m³] (`iteration variable 6`)"""

    nd_plasma_fuel_ions_vol_avg: float = 0.0
    """Plasma volume averaged fuel ion density (⟨n_fuel⟩) [/m³]"""

    dlamee: float = 0.0
    """electron-electron coulomb logarithm"""

    dlamie: float = 0.0
    """ion-electron coulomb logarithm"""

    nd_plasma_electron_max_array: list[float] = field(
        default_factory=lambda: np.zeros(8, dtype=np.float64)
    )
    """Array of plasma electron density upper limits values (nₑ,max) [/m³]"""

    nd_plasma_alphas_thermal_vol_avg: float = 0.0
    """Plasma volume averaged thermal alpha density (⟨n_αₜₕ⟩) [/m³]"""

    nd_beam_ions: float = 0.0
    """Hot beam ion density, variable (⟨n_beam⟩) [/m³]"""

    nd_beam_ions_out: float = 0.0
    """Hot beam ion density from calculation [/m³]"""

    beta_norm_max: float = 3.5
    """Troyon-like coefficient for beta scaling"""

    beta_norm_max_wesson: float = 0.0
    """Wesson-like coefficient for beta scaling"""

    beta_norm_max_menard: float = 0.0
    """Menard-like coefficient for beta scaling"""

    beta_norm_max_original_scaling: float = 0.0
    """Original scaling coefficient for beta scaling"""

    beta_norm_max_tholerus: float = 0.0
    """Tholerus-like coefficient for beta scaling"""

    beta_norm_max_stambaugh: float = 0.0
    """Stambaugh-like coefficient for beta scaling"""

    nd_plasma_electrons_max: float = 0.0
    """Plasma electron max density limit (nₑ,max) [/m³]"""

    nd_plasma_ions_total_vol_avg: float = 0.0
    """Plasma volume averaged total ion density (⟨n_i⟩) [/m³]"""

    nd_plasma_electron_line: float = 0.0
    """Plasma line averaged electron density (⟨nₑ⟩_line) [/m³]"""

    nd_plasma_protons_vol_avg: float = 0.0
    """Plasma volume averaged proton ash density (⟨n_p⟩) [/m³]"""

    ntau: float = 0.0
    """Fusion double product [s/m³]"""

    nTtau: float = 0.0
    """Lawson triple product [keV s / m³]"""

    nd_plasma_impurities_vol_avg: float = 0.0
    """Plasma volume averaged impurity (Z > 2) ion density (⟨n_imp⟩) [/m³]"""

    gradient_length_ne: float = None
    """Max. normalised gradient length in el. density (i_plasma_pedestal==0 only)"""

    gradient_length_te: float = None
    """Max. normalised gradient length in el. temperature (i_plasma_pedestal==0 only)"""

    beta_poloidal_eps_max: float = 1.38
    """maximum (eps*beta_poloidal) (`constraint equation 6`). Note: revised issue #346
    "Operation at the tokamak equilibrium poloidal beta-limit in TFTR", 1992 Nucl. Fusion 32 1468
    """

    eps: float = 0.34399724802
    """Plasma inverse aspect ratio (ε)"""

    f_c_plasma_auxiliary: float = 0.0
    """fraction of plasma current produced by auxiliary current drive"""

    f_c_plasma_inductive: float = 0.0
    """fraction of plasma current produced inductively"""

    f_alpha_electron: float = 0.0
    """fraction of alpha energy to electrons"""

    f_p_alpha_plasma_deposited: float = 0.95
    """Fraction of alpha power deposited in plasma. Default of 0.95 taken from https://doi.org/10.1088/0029-5515/39/12/305."""

    f_alpha_ion: float = 0.0
    """fraction of alpha power to ions"""

    f_plasma_fuel_deuterium: float = 0.5
    """Plasma deuterium fuel fraction"""

    f_p_div_lower: float = 1.0
    """fraction of power to the lower divertor in double null configuration
    (`i_single_null = 0` only) (default assumes SN)
    """

    ffwal: float = 0.92
    """factor to convert plasma surface area to first wall area in neutron wall
    load calculation (`i_pflux_fw_neutron=1`)
    """

    f_nd_plasma_greenwald: float = None
    """Greenwald fraction of the line averaged electron density. The classic Greenwald
    limit value"""

    f_nd_plasma_pedestal_greenwald: float = 0.85
    """Greenwald fraction of the pedestal density
    """

    f_nd_plasma_separatrix_greenwald: float = 0.5
    """Greenwald fraction of the separatrix density
    """

    f_plasma_fuel_helium3: float = 0.0
    """Plasma Helium-3 fuel fraction"""

    figmer: float = 0.0
    """physics figure of merit (= plasma_current*aspect**sbar, where `sbar=1`)"""

    fkzohm: float = 1.0
    """Zohm elongation scaling adjustment factor (`i_plasma_geometry=2, 3`)"""

    f_plasma_fuel_tritium: float = 0.5
    """Plasma tritium fuel fraction"""

    fusden_total: float = 0.0
    """fusion reaction rate density, from beams and plasma (reactions/m3/sec)"""

    fusrat_total: float = 0.0
    """fusion reaction rate, from beams and plasma (reactions/sec)"""

    fusrat_plasma_dt_profile: list[float] = field(default_factory=list)
    """Profile of D-T fusion reaction rate in plasma, (reactions/sec)"""

    fusrat_plasma_dd_triton_profile: list[float] = field(default_factory=list)
    """Profile of D-D fusion reaction rate (tritium branch) in plasma, (reactions/sec)"""

    fusrat_plasma_dd_helion_profile: list[float] = field(default_factory=list)
    """Profile of D-D fusion reaction rate (helium branch) in plasma, (reactions/sec)"""

    fusrat_plasma_dhe3_profile: list[float] = field(default_factory=list)
    """Profile of D-3He fusion reaction rate in plasma, (reactions/sec)"""

    fusden_plasma: float = 0.0
    """fusion reaction rate, just from plasma (reactions/m3/sec)"""

    f_c_plasma_non_inductive: float = 1.0
    """fraction of the plasma current produced by non-inductive means (`iteration variable 44`)"""

    ejima_coeff: float = 0.4
    """Ejima coefficient for resistive startup V-s formula"""

    f_beta_alpha_beam_thermal: float = 0.0
    """ratio of (fast alpha + neutral beam beta) to thermal beta"""

    hfac: list[float] = field(
        default_factory=lambda: np.zeros(N_CONFINEMENT_SCALINGS, dtype=np.float64)
    )
    """H factors for an ignited plasma for each energy confinement time scaling law"""

    hfact: float = 1.0
    """H factor on energy confinement times, radiation corrected (`iteration variable 10`)."""

    hstar: float = 1.0
    """H* non-radiation corrected H factor on energy confinement times"""

    t_plasma_energy_confinement_max: float = 10.0
    """Maximum allowed energy confinement time (s)"""

    i_bootstrap_current: int = 3
    """switch for bootstrap current scaling
    - =1 ITER 1989 bootstrap scaling (high R/a only)
    - =2 for Nevins et al general scaling
    - =3 for Wilson et al numerical scaling
    - =4 for Sauter et al scaling
    - =5 for Sakai et al scaling
    - =6 for ARIES scaling
    - =7 for Andrade et al scaling
    - =8 for Hoang et al scaling
    - =9 for Wong et al scaling
    - =10 for Gi-I et al scaling
    - =11 for Gi-II et al scaling
    - =12 for Sugiyama (L-mode) et al scaling
    - =13 for Sugiyama (H-mode) et al scaling
    """

    i_beta_component: int = 0
    """switch for beta limit scaling (`constraint equation 24`)
    - =0 apply limit to total beta
    - =1 apply limit to thermal beta
    - =2 apply limit to thermal + neutral beam beta
    - =3 apply limit to toroidal beta
    """

    i_plasma_current: int = 4
    """switch for plasma current scaling to use
    - =1 Peng analytic fit
    - =2 Peng double null divertor scaling (ST)
    - =3 simple ITER scaling (k = 2.2, d = 0.6)
    - =4 later ITER scaling, a la Uckan
    - =5 Todd empirical scaling I
    - =6 Todd empirical scaling II
    - =7 Connor-Hastie model
    - =8 Sauter scaling allowing negative triangularity
    - =9 FIESTA ST fit
    """

    i_diamagnetic_current: int = 0
    """switch for diamagnetic current scaling
    - =0 Do not calculate
    - =1 Use original TART scaling
    - =2 Use SCENE scaling
    """

    i_density_limit: int = 8
    """switch for density limit to enforce (`constraint equation 5`)
    - =1 old ASDEX
    - =2 Borrass model for ITER (I)
    - =3 Borrass model for ITER (II)
    - =4 JET edge radiation
    - =5 JET simplified
    - =6 Hugill-Murakami Mq limit
    - =7 Greenwald limit
    - =8 ASDEX New
    """

    i_beta_fast_alpha: int = 1
    """switch for fast alpha pressure calculation
    - =0 ITER physics rules (Uckan) fit
    - =1 Modified fit (D. Ward) - better at high temperature
    """

    i_plasma_ignited: int = 0
    """switch for ignition assumption. Obviously, i_plasma_ignited must be zero if current drive
    is required. If i_plasma_ignited is 1, any auxiliary power is assumed to be used only during
    plasma start-up, and is excluded from all steady-state power balance calculations.
    - =0 do not assume plasma ignition
    - =1 assume ignited (but include auxiliary power in costs)</UL
    """

    i_plasma_pedestal: int = 1
    """switch for pedestal profiles:
    - =0 use original parabolic profiles
    - =1 use pedestal profile
    """

    i_pfirsch_schluter_current: int = 0
    """switch for Pfirsch-Schlüter current scaling (issue #413):
    - =0 Do not calculate
    - =1 Use SCENE scaling
    """

    nd_plasma_pedestal_electron: float = 4.0e19
    """Plasma electron density at pedestal (nₑ,pedestal) [/m³] (`i_plasma_pedestal==1)"""

    nd_plasma_separatrix_electron: float = 3.0e19
    """Plasma electron density at separatrix (nₑ,sep) [/m³] (`i_plasma_pedestal==1)"""

    i_nd_plasma_pedestal_separatrix: int = 1
    """switch for pedestal and separatrix density calculation:
    - =0 User input pedestal and separatrix density
    - =1 Calculate pedestal and separatrix density as fraction of Greenwald limit (see `f_nd_plasma_pedestal_greenwald` and `f_nd_plasma_separatrix_greenwald`)
    """

    alpha_crit: float = 0.0
    """critical ballooning parameter value"""

    nd_plasma_separatrix_electron_eich_max: float = 0.0
    """Eich critical electron density at separatrix [/m³]"""

    plasma_res_factor: float = 1.0
    """plasma resistivity pre-factor"""

    radius_plasma_pedestal_density_norm: float = 1.0
    """PPlasma normalised radius of density pedestal (ρₙ,pedestal)  (`i_plasma_pedestal==1`)"""

    radius_plasma_pedestal_temp_norm: float = 1.0
    """Plasma normalised radius of temperature pedestal (ρₜ,pedestal) (`i_plasma_pedestal==1`)"""

    rho_te_max: float = 0.0
    """r/a where the temperature gradient is largest (`i_plasma_pedestal==0`)"""

    rho_ne_max: float = 0.0
    """r/a where the density gradient is largest (`i_plasma_pedestal==0`)"""

    tbeta: float = 2.0
    """Plasma temperature profile index beta (βₜ)  (`i_plasma_pedestal==1)"""

    temp_plasma_pedestal_kev: float = 1.0
    """Plasma electron temperature of pedestal (Tₑ,pedestal) [keV] (`i_plasma_pedestal==1`)"""

    temp_plasma_separatrix_kev: float = 0.1
    """Plasma electron temperature at separatrix (Tₑ,ₛₑₚ) [keV] (`i_plasma_pedestal==1`) calculated if reinke
    criterion is used (`icc=78`)
    """

    i_beta_norm_max: int = 1
    """Switch for maximum normalised beta scaling (βₙ)"""

    i_ind_plasma_internal_norm: int = 0
    """Switch for plasma normalised internal inductance scaling (lᵢ)"""

    i_alphaj: int = 0
    """Switch for plasma current profile index scaling (αⱼ) """

    i_rad_loss: int = 1
    """switch for radiation loss term usage in power balance (see User Guide):
    - =0 total power lost is scaling power plus radiation
    - =1 total power lost is scaling power plus core radiation only
    - =2 total power lost is scaling power only, with no additional
    allowance for radiation. This is not recommended for power plant models.
    """

    i_confinement_time: int = 34
    """Switch for plasma energy confinement time scaling law (τₑ)"""

    i_plasma_wall_gap: int = 1
    """Switch for plasma-first wall clearances at the mid-plane:
    - =0 use 10% of plasma minor radius
    - =1 use input (`dr_fw_plasma_gap_inboard` and `dr_fw_plasma_gap_outboard`)
    """

    i_plasma_geometry: int = 0
    """switch for plasma elongation and triangularity calculations:
    - =0 use input kappa, triang to calculate 95% values
    - =1 scale q95_min, kappa, triang with aspect ratio (ST)
    - =2 set kappa to the natural elongation value (Zohm ITER scaling), triang input
    - =3 set kappa to the natural elongation value (Zohm ITER scaling), triang95 input
    - =4 use input kappa95, triang95 to calculate separatrix values
    - =5 use input kappa95, triang95 to calculate separatrix values based on MAST scaling (ST)
    - =6 use input kappa, triang to calculate 95% values based on MAST scaling (ST)
    - =7 use input kappa95, triang95 to calculate separatrix values based on fit to FIESTA (ST)
    - =8 use input kappa, triang to calculate 95% values based on fit to FIESTA (ST)
    - =9 set kappa to the natural elongation value, triang input
    - =10 set kappa to maximum stable value at a given aspect ratio (2.6<A<3.6)), triang input (#1399)
    - =11 set kappa Menard 2016 aspect-ratio-dependent scaling, triang input (#1439)
    - =12 set kappa Menard 1997 aspect-ratio-dependent scaling, triang input
    """

    i_plasma_shape: int = 0
    """switch for plasma boundary shape:
    - =0 use original PROCESS 2-arcs model
    - =1 use the Sauter model
    """

    itart: int = 0
    """switch for spherical tokamak (ST) models:
    - =0 use conventional aspect ratio models
    - =1 use spherical tokamak models
    """

    itartpf: int = 0
    """switch for Spherical Tokamak PF models:
    - =0 use Peng and Strickler (1986) model
    - =1 use conventional aspect ratio model
    """

    i_pflux_fw_neutron: int = 1
    """switch for neutron wall load calculation:
    - =1 use scaled plasma surface area
    - =2 use first wall area directly
    """

    plasma_square: float = 0.0
    """Plasma squareness (ζ)"""

    kappa: float = 1.792
    """Plasma separatrix elongation (κₐ)  (calculated if `i_plasma_geometry = 1-5, 7 or 9-10`)"""

    kappa95: float = 1.6
    """Plasma elongation at 95% surface (κ₉₅) (calculated if `i_plasma_geometry = 0-3, 6, or 8-10`)"""

    kappa_ipb: float = 0.0
    """Separatrix elongation calculated for IPB scalings"""

    nd_plasma_electron_on_axis: float = 0.0
    """central electron density (/m3)"""

    nd_plasma_ions_on_axis: float = 0.0
    """central ion density (/m3)"""

    m_s_limit: float = 0.3
    """margin to vertical stability"""

    pres_plasma_thermal_on_axis: float = 0.0
    """Plasma central thermal pressure (p₀) (no fast ions or beam pressure) [Pa]"""

    pres_plasma_thermal_total_profile: list[float] = field(default_factory=list)
    """Profile of total pressure in plasma [Pa]"""

    pres_plasma_electron_profile: list[float] = field(default_factory=list)
    """Profile of electron pressure in plasma [Pa]"""

    pres_plasma_ion_total_profile: list[float] = field(default_factory=list)
    """Profile of ion pressure in plasma [Pa]"""

    pres_plasma_fuel_profile: list[float] = field(default_factory=list)
    """Profile of fuel pressure in plasma [Pa]"""

    j_plasma_on_axis: float = 0.0
    """Central plasma current density (j₀) [A/m²]"""

    j_plasma_bootstrap_sauter_profile: list[float] = field(default_factory=list)
    """Profile of bootstrap current density in plasma using Sauter et al scaling [A/m²]"""

    n_plasma_profile_elements: int = 501
    """Number of elements in plasma profile"""

    pres_plasma_thermal_vol_avg: float = None
    """Volume averaged thermal plasma pressure (⟨p⟩)  (no fast ions or beam pressure) [Pa]"""

    f_dd_branching_trit: float = 0.0
    """branching ratio for DD -> T"""

    pden_plasma_alpha_mw: float = 0.0
    """Alpha power per volume just from plasma [MW/m3]"""

    pden_alpha_total_mw: float = 0.0
    """Alpha power per volume from plasma and beams [MW/m3]"""

    f_pden_alpha_electron_mw: float = 0.0
    """Alpha power per volume to electrons [MW/m3]"""

    p_fw_alpha_mw: float = 0.0
    """alpha power escaping plasma and reaching first wall (MW)"""

    f_pden_alpha_ions_mw: float = 0.0
    """alpha power per volume to ions (MW/m3)"""

    p_plasma_alpha_mw: float = 0.0
    """Alpha power from only the plasma (MW)"""

    p_alpha_total_mw: float = 0.0
    """Total alpha power from plasma and beams (MW)"""

    p_beam_alpha_mw: float = 0.0
    """alpha power from hot neutral beam ions (MW)"""

    p_beam_neutron_mw: float = 0.0
    """neutron power from hot neutral beam ions (MW)"""

    p_beam_dt_mw: float = 0.0
    """D-T fusion power from hot neutral beam ions (MW)"""

    p_non_alpha_charged_mw: float = 0.0
    """non-alpha charged particle fusion power (MW)"""

    p_charged_particle_mw: float = 0.0
    """Total charged particle fusion power [MW]"""

    pden_non_alpha_charged_mw: float = 0.0
    """Non-alpha charged particle fusion power per volume [MW/m3]"""

    f_temp_plasma_electron_density_vol_avg: float = 0.0
    """Ratio of density weighted plasma electron tempertaurature to volume averaged (Profile Factor)"""

    p_plasma_inner_rad_mw: float = 0.0
    """radiation power from inner zone (MW)"""

    pden_plasma_core_rad_mw: float = 0.0
    """total core radiation power per volume (MW/m3)"""

    p_dd_total_mw: float = 0.0
    """deuterium-deuterium fusion power (MW)"""

    p_dhe3_total_mw: float = 0.0
    """deuterium-helium3 fusion power (MW)"""

    p_plasma_separatrix_mw: float = 0.0
    """power to conducted to the divertor region (MW)"""

    p_plasma_separatrix_rmajor_mw: float = 0.0
    """Power to conducted to the divertor region per major radius (MW/m)"""

    p_div_bt_q_aspect_rmajor_mw: float = 0.0
    """EU DEMO divertor protection parameter (PₛₑₚBₜ / q₉₅AR₀)  [MWT/m]"""

    p_div_lower_separatrix_mw: float = 0.0
    """Separatrix power conducted to the lower divertor region (calculated if `i_single_null = 0`) (MW)"""

    p_div_upper_separatrix_mw: float = 0.0
    """Separatrix power conducted to the upper divertor region (calculated if `i_single_null = 0`) (MW)"""

    p_div_separatrix_max_mw: float = 0.0
    """Separatrix power conducted to the divertor with most load (calculated if `i_single_null = 0`) (MW)"""

    p_dt_total_mw: float = 0.0
    """Total deuterium-tritium fusion power, from plasma and beams [MW]"""

    p_plasma_dt_mw: float = 0.0
    """Deuterium-tritium fusion power, just from plasma [MW]"""

    p_plasma_outer_rad_mw: float = 0.0
    """radiation power from outer zone (MW)"""

    pden_plasma_outer_rad_mw: float = 0.0
    """edge radiation power per volume (MW/m3)"""

    vs_plasma_internal: float = 0.0
    """internal plasma V-s"""

    pflux_fw_rad_mw: float = 0.0
    """Nominal mean radiation load on inside surface of reactor (MW/m2)"""

    pden_ion_electron_equilibration_mw: float = 0.0
    """ion/electron equilibration power per volume (MW/m3)"""

    plasma_current: float = 0.0
    """Plasma current (Iₚ) [A]"""

    c_plasma_peng_analytic: float = 0.0
    """Peng analytic plasma current (A)"""

    c_plasma_peng_double_null: float = 0.0
    """Peng double null divertor plasma current (A)"""

    c_plasma_cyclindrical: float = 0.0
    """Cylindrical plasma current (A)"""

    c_plasma_ipdg89: float = 0.0
    """ITER IPDG89 plasma current (A)"""

    c_plasma_todd_empirical_i: float = 0.0
    """Todd empirical plasma current I (A)"""

    c_plasma_todd_empirical_ii: float = 0.0
    """Todd empirical plasma current II (A)"""

    c_plasma_connor_hastie: float = 0.0
    """Connor-Hastie plasma current (A)"""

    c_plasma_sauter: float = 0.0
    """Sauter plasma current (A)"""

    c_plasma_fiesta_st: float = 0.0
    """FIESTA ST plasma current (A)"""

    p_plasma_neutron_mw: float = 0.0
    """Neutron fusion power from just the plasma [MW]"""

    p_neutron_total_mw: float = 0.0
    """Total neutron fusion power from plasma and beams [MW]"""

    pden_neutron_total_mw: float = 0.0
    """neutron fusion power per volume from beams and plasma (MW/m3)"""

    pden_plasma_neutron_mw: float = 0.0
    """neutron fusion power per volume just from plasma (MW/m3)"""

    p_plasma_ohmic_mw: float = 0.0
    """ohmic heating power (MW)"""

    pden_plasma_ohmic_mw: float = 0.0
    """ohmic heating power per volume (MW/m3)"""

    p_plasma_loss_mw: float = 0.0
    """heating power (= transport loss power) (MW) used in confinement time calculation"""

    p_fusion_total_mw: float = 0.0
    """fusion power (MW)"""

    len_plasma_poloidal: float = 0.0
    """plasma poloidal perimeter (m)"""

    p_plasma_rad_mw: float = 0.0
    """total radiation power from inside LCFS (MW)"""

    pden_plasma_rad_mw: float = 0.0
    """total radiation power per volume (MW/m3)"""

    pradsolmw: float = 0.0
    """radiation power from SoL (MW)"""

    proton_rate_density: float = 0.0
    """Proton production rate [particles/m3/sec]"""

    psolradmw: float = 0.0
    """SOL radiation power (MW) (`stellarator only`)"""

    pden_plasma_sync_mw: float = 0.0
    """Plasma synchrotron radiation power per unit volume [MW/m³]"""

    p_plasma_sync_mw: float = 0.0
    """Total synchrotron radiation power from plasma (Pₛₙ) [MW]"""

    i_l_h_threshold: int = 19
    """switch for L-H mode power threshold scaling to use (see l_h_threshold_powers for list)"""

    p_l_h_threshold_mw: float = 0.0
    """L-H mode power threshold (MW) (chosen via i_l_h_threshold, and enforced if
    constraint equation 15 is on)
    """

    l_h_threshold_powers: list[float] = field(
        default_factory=lambda: np.zeros(21, dtype=np.float64)
    )
    """L-H power threshold for various scalings (MW)
    - =1 ITER 1996 scaling: nominal
    - =2 ITER 1996 scaling: upper bound
    - =3 ITER 1996 scaling: lower bound
    - =4 ITER 1997 scaling: excluding elongation
    - =5 ITER 1997 scaling: including elongation
    - =6 Martin 2008 scaling: nominal
    - =7 Martin 2008 scaling: 95% upper bound
    - =8 Martin 2008 scaling: 95% lower bound
    - =9 Snipes 2000 scaling: nominal
    - =10 Snipes 2000 scaling: upper bound
    - =11 Snipes 2000 scaling: lower bound
    - =12 Snipes 2000 scaling (closed divertor): nominal
    - =13 Snipes 2000 scaling (closed divertor): upper bound
    - =14 Snipes 2000 scaling (closed divertor): lower bound
    - =15 Hubbard et al. 2012 L-I threshold scaling: nominal
    - =16 Hubbard et al. 2012 L-I threshold scaling: lower bound
    - =17 Hubbard et al. 2012 L-I threshold scaling: upper bound
    - =18 Hubbard et al. 2017 L-I threshold scaling
    - =19 Martin 2008 aspect ratio corrected scaling: nominal
    - =20 Martin 2008 aspect ratio corrected scaling: 95% upper bound
    - =21 Martin 2008 aspect ratio corrected scaling: 95% lower bound
    """

    p_electron_transport_loss_mw: float = 0.0
    """electron transport power (MW)"""

    pden_electron_transport_loss_mw: float = 0.0
    """electron transport power per volume (MW/m3)"""

    p_ion_transport_loss_mw: float = 0.0
    """ion transport power (MW)"""

    pscalingmw: float = 0.0
    """Total transport power from scaling law (MW)"""

    pden_ion_transport_loss_mw: float = 0.0
    """ion transport power per volume (MW/m3)"""

    q0: float = 1.0
    """Plasma safety factor on axis (q₀)"""

    q95: float = 0.0
    """Plasma safety factor at 95% flux surface (q₉₅) (`iteration variable 18`)
    """

    molflow_plasma_fuelling_required: float = 0.0
    """plasma fuelling rate (nucleus-pairs/s)"""

    tauratio: float = 1.0
    """tauratio /1.0/ : ratio of He and pellet particle confinement times"""

    q95_min: float = 0.0
    """Plasmalower limit for edge safety factor"""

    qstar: float = 0.0
    """Plasma cylindrical safety factor (qcyl)"""

    rad_fraction_sol: float = 0.8
    """SoL radiation fraction"""

    rad_fraction_total: float = 0.0
    """Radiation fraction total = SoL + LCFS radiation / total power deposited in plasma"""

    f_nd_alpha_thermal_electron: float = 0.1
    """Thermal alpha density/electron density (⟨n_αₜₕ⟩/⟨nₑ⟩)"""

    f_nd_protium_electrons: float = 0.0
    """Seeded f_nd_protium_electrons density / electron density."""

    ind_plasma_internal_norm: float = 0.9
    """Plasma normalised internal inductance"""

    ind_plasma_internal_norm_iter_3: float = 0.0
    """Plasma normalised internal inductance (ITER type 3)"""

    ind_plasma_internal_norm_wesson: float = 0.0
    """Wesson-like plasma normalised internal inductance"""

    ind_plasma_internal_norm_menard: float = 0.0
    """Menard-like plasma normalised internal inductance"""

    ind_plasma: float = 0.0
    """plasma inductance (H)"""

    rmajor: float = 8.14
    """Plasma major radius (R₀) [m] (`iteration variable 3`)"""

    rminor: float = 0.0
    """Plasma minor radius (a) [m]"""

    f_nd_beam_electron: float = 0.005
    """hot beam density / n_e (`iteration variable 7`)"""

    f_nd_plasma_carbon_electron: float = 0.0
    """n_carbon / n_e"""

    rndfuel: float = 0.0
    """fuel burnup rate (reactions/second)"""

    f_nd_plasma_iron_argon_electron: float = 0.0
    """n_highZ / n_e"""

    f_nd_plasma_oxygen_electron: float = 0.0
    """n_oxygen / n_e"""

    f_res_plasma_neo: float = 0.0
    """neo-classical correction factor to res_plasma"""

    res_plasma: float = 0.0
    """plasma resistance (ohm)"""

    t_plasma_res_diffusion: float = 0.0
    """plasma current resistive diffusion time (s)"""

    a_plasma_surface: float = 0.0
    """plasma surface area"""

    a_plasma_surface_outboard: float = 0.0
    """outboard plasma surface area"""

    i_single_null: int = 1
    """switch for single null / double null plasma:
    - =0 for double null
    - =1 for single null (diverted side down)
    """

    f_sync_reflect: float = 0.6
    """synchrotron wall reflectivity factor"""

    t_electron_energy_confinement: float = 0.0
    """electron energy confinement time (sec)"""

    tauee_in: float = 0.0
    """Input electron energy confinement time (sec) (`i_confinement_time=48 only`)"""

    t_energy_confinement: float = 0.0
    """global thermal energy confinement time (sec)"""

    t_ion_energy_confinement: float = 0.0
    """ion energy confinement time (sec)"""

    t_alpha_confinement: float = 0.0
    """alpha particle confinement time (sec)"""

    f_t_alpha_energy_confinement: float = 0.0
    """Alpha particle to energy confinement time ratio"""

    temp_plasma_electron_vol_avg_kev: float = 12.9
    """Plasma volume averaged electron temperature (⟨Tₑ⟩) [keV] (`iteration variable 4`)"""

    temp_plasma_electron_on_axis_kev: float = 0.0
    """Plasma central electron temperature (Tₑ₀) [keV]"""

    temp_plasma_electron_density_weighted_kev: float = 0.0
    """Density weighted average electron temperature (⟨Tₑ⟩_n) [keV]"""

    temp_plasma_electron_line_avg_kev: float = None
    """Line averaged electron temperature (keV)"""

    temp_plasma_ion_vol_avg_kev: float = 12.9
    """Volume averaged ion temperature (⟨Tᵢ⟩) [keV]. N.B. calculated from temp_plasma_electron_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`"""

    temp_plasma_ion_on_axis_kev: float = 0.0
    """Plasma central ion temperature (Tᵢ₀) [keV]"""

    temp_plasma_ion_density_weighted_kev: float = 0.0
    """Plasma density weighted average ion temperature (⟨Tᵢ⟩ₙ) [keV]"""

    f_temp_plasma_ion_electron: float = 1.0
    """Plasma ratio of ion temperature to electron temperature (used to calculate temp_plasma_ion_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`)"""

    triang: float = 0.36
    """Plasma separatrix triangularity (δₐ) (calculated if `i_plasma_geometry = 1, 3-5 or 7`)"""

    triang95: float = 0.24
    """Plasma triangularity at 95% surface (δ₉₅) (calculated if `i_plasma_geometry = 0-2, 6, 8 or 9`)"""

    vol_plasma: float = 0.0
    """Plasma volume [m³]"""

    vs_plasma_burn_required: float = 0.0
    """V-s needed during flat-top (heat + burn times) (Wb)"""

    vs_plasma_ramp_required: float = 0.0
    """V-s needed during ramp-up (Wb)"""

    v_plasma_loop_burn: float = 0.0
    """Plasma loop voltage during flat-top (V)"""

    vs_plasma_ind_ramp: float = 0.0
    """Total plasma inductive flux consumption for plasma current ramp-up (Vs)(Wb)"""

    vs_plasma_res_ramp: float = 0.0
    """Plasma resistive flux consumption for plasma current ramp-up (Vs)(Wb)"""

    vs_plasma_total_required: float = 0.0
    """Total V-s needed for full plasma pulse [Wb]"""

    pflux_fw_neutron_mw: float = 0.0
    """Average FW neutron wall load [MW/m²]"""

    pflux_plasma_surface_neutron_avg_mw: float = 0.0
    """Average neutron flux at plasma surface [MW/m²]"""

    wtgpd: float = 0.0
    """Mass of fuel used per day [g]"""

    a_plasma_poloidal: float = 0.0
    """Plasma poloidal cross-sectional area [m²]"""

    n_charge_plasma_effective_vol_avg: float = 0.0
    """Volume averaged plasma effective charge (⟨Zₑ⟩)"""

    n_charge_plasma_effective_profile: list[float] = field(default_factory=list)
    """Profile of plasma effective charge"""

    n_charge_plasma_effective_mass_weighted_vol_avg: float = 0.0
    """Plasma mass-weighted volume averaged plasma effective charge (⟨Zₑ⟩ₘ)"""

    len_plasma_debye_electron_profile: list[float] = field(default_factory=list)
    """Profile of electron Debye length in plasma (m)"""

    radius_plasma_deuteron_toroidal_larmor_isotropic_profile: list[float] = field(
        default_factory=list
    )
    """Profile of deuteron toroidal Larmor radius in plasma, assuming equal speeds in all directions (m)"""

    radius_plasma_deuteron_toroidal_larmor_isotropic_vol_avg: float = 0.0
    """Volume averaged deuteron toroidal Larmor radius in plasma, assuming equal speeds in all directions (m)"""

    radius_plasma_triton_toroidal_larmor_isotropic_profile: list[float] = field(
        default_factory=list
    )
    """Profile of triton toroidal Larmor radius in plasma, assuming equal speeds in all directions (m)"""

    radius_plasma_triton_toroidal_larmor_isotropic_vol_avg: float = 0.0
    """Volume averaged triton toroidal Larmor radius in plasma, assuming equal speeds in all directions (m)"""

    len_plasma_debye_electron_vol_avg: float = 0.0
    """Volume averaged electron Debye length in plasma (m)"""

    vel_plasma_electron_profile: list[float] = field(default_factory=list)
    """Profile of electron thermal velocity in plasma (m/s)"""

    vel_plasma_deuteron_vol_avg: float = 0.0
    """Volume averaged deuteron thermal velocity in plasma (m/s)"""

    vel_plasma_electron_vol_avg: float = 0.0
    """Volume averaged electron thermal velocity in plasma (m/s)"""

    vel_plasma_deuteron_profile: list[float] = field(default_factory=list)
    """Profile of deuteron thermal velocity in plasma (m/s)"""

    vel_plasma_triton_profile: list[float] = field(default_factory=list)
    """Profile of triton thermal velocity in plasma (m/s)"""

    vel_plasma_triton_vol_avg: float = 0.0
    """Volume averaged triton thermal velocity in plasma (m/s)"""

    vel_plasma_alpha_thermal_profile: list[float] = field(default_factory=list)
    """Profile of thermal alpha particle velocity in plasma (m/s)"""

    vel_plasma_alpha_thermal_vol_avg: float = 0.0
    """Volume averaged thermal alpha particle velocity in plasma (m/s)"""

    vel_plasma_alpha_birth: float = 0.0
    """Birth velocity of alpha particles in plasma (m/s)"""

    plasma_coulomb_log_electron_electron_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-electron Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_electron_vol_avg: float = 0.0
    """Volume averaged electron-electron Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_deuteron_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-deuteron Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_deuteron_vol_avg: float = 0.0
    """Volume averaged electron-deuteron Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_triton_profile: list[float] = field(default_factory=list)
    """Profile of electron-triton Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_triton_vol_avg: float = 0.0
    """Volume averaged electron-triton Coulomb logarithm in plasma"""

    plasma_coulomb_log_deuteron_triton_profile: list[float] = field(default_factory=list)
    """Profile of deuteron-triton Coulomb logarithm in plasma"""

    plasma_coulomb_log_deuteron_triton_vol_avg: float = 0.0
    """Volume averaged deuteron-triton Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_alpha_thermal_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-alpha Coulomb logarithm in plasma"""

    plasma_coulomb_log_electron_alpha_thermal_vol_avg: float = 0.0
    """Volume averaged electron-alpha Coulomb logarithm in plasma"""

    t_plasma_electron_alpha_spitzer_slow_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-alpha Spitzer slowing down time in plasma (s)"""

    t_plasma_electron_alpha_spitzer_slow_vol_avg: float = 0.0
    """Volume averaged electron-alpha Spitzer slowing down time in plasma (s)"""

    freq_plasma_electron_profile: list[float] = field(default_factory=list)
    """Electron plasma frequency profile (Hz)"""

    freq_plasma_electron_vol_avg: float = 0.0
    """Volume averaged electron plasma frequency (Hz)"""

    freq_plasma_deuteron_profile: list[float] = field(default_factory=list)
    """Deuteron plasma frequency profile (Hz)"""

    freq_plasma_larmor_toroidal_electron_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_larmor_toroidal_electron_vol_avg: float = None
    """Volume averaged electron Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_larmor_toroidal_deuteron_profile: list[float] = field(
        default_factory=list
    )
    """Profile of deuteron Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_larmor_toroidal_deuteron_vol_avg: float = None
    """Volume averaged deuteron Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_larmor_toroidal_triton_profile: list[float] = field(default_factory=list)
    """Profile of triton Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_larmor_toroidal_triton_vol_avg: float = None
    """Volume averaged triton Larmor frequency in plasma due to toroidal magnetic field (Hz)"""

    freq_plasma_upper_hybrid_profile: list[float] = field(default_factory=list)
    """Profile of upper hybrid frequency in plasma (Hz)"""

    freq_plasma_upper_hybrid_vol_avg: float = 0.0
    """Volume averaged upper hybrid frequency in plasma (Hz)"""

    t_plasma_electron_electron_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-electron collision time in plasma (s)"""

    t_plasma_electron_electron_collision_vol_avg: float = 0.0
    """Volume averaged electron-electron collision time in plasma (s)"""

    t_plasma_electron_deuteron_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-deuteron collision time in plasma (s)"""

    t_plasma_electron_deuteron_collision_vol_avg: float = 0.0
    """Volume averaged electron-deuteron collision time in plasma (s)"""

    t_plasma_electron_triton_collision_profile: list[float] = field(default_factory=list)
    """Profile of electron-triton collision time in plasma (s)"""

    t_plasma_electron_triton_collision_vol_avg: float = 0.0
    """Volume averaged electron-triton collision time in plasma (s)"""

    t_plasma_electron_alpha_thermal_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-alpha collision time in plasma (s)"""

    t_plasma_electron_alpha_thermal_collision_vol_avg: float = 0.0
    """Volume averaged electron-alpha collision time in plasma (s)"""

    freq_plasma_electron_electron_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-electron collision frequency in plasma (Hz)"""

    freq_plasma_electron_electron_collision_vol_avg: float = 0.0
    """Volume averaged electron-electron collision frequency in plasma (Hz)"""

    freq_plasma_electron_deuteron_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-deuteron collision frequency in plasma (Hz)"""

    freq_plasma_electron_deuteron_collision_vol_avg: float = 0.0
    """Volume averaged electron-deuteron collision frequency in plasma (Hz)"""

    freq_plasma_electron_triton_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-triton collision frequency in plasma (Hz)"""

    freq_plasma_electron_triton_collision_vol_avg: float = 0.0
    """Volume averaged electron-triton collision frequency in plasma (Hz)"""

    freq_plasma_electron_alpha_thermal_collision_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-alpha collision frequency in plasma (Hz)"""

    freq_plasma_electron_alpha_thermal_collision_vol_avg: float = 0.0
    """Volume averaged electron-alpha collision frequency in plasma (Hz)"""

    len_plasma_electron_electron_mean_free_path_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-electron mean free path in plasma (m)"""

    len_plasma_electron_electron_mean_free_path_vol_avg: float = 0.0
    """Volume averaged electron-electron mean free path in plasma (m)"""

    len_plasma_electron_deuteron_mean_free_path_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-deuteron mean free path in plasma (m)"""

    len_plasma_electron_deuteron_mean_free_path_vol_avg: float = 0.0
    """Volume averaged electron-deuteron mean free path in plasma (m)"""

    len_plasma_electron_triton_mean_free_path_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-triton mean free path in plasma (m)"""

    len_plasma_electron_triton_mean_free_path_vol_avg: float = 0.0
    """Volume averaged electron-triton mean free path in plasma (m)"""

    len_plasma_electron_alpha_thermal_mean_free_path_profile: list[float] = field(
        default_factory=list
    )
    """Profile of electron-alpha mean free path in plasma (m)"""

    len_plasma_electron_alpha_thermal_mean_free_path_vol_avg: float = 0.0
    """Volume averaged electron-alpha mean free path in plasma (m)"""

    res_plasma_fuel_spitzer_profile: list[float] = field(default_factory=list)
    """Profile of plasma Spitzer resistivity due to fuel ions (ohm m)"""

    res_plasma_fuel_spitzer_vol_avg: float = 0.0
    """Volume averaged plasma Spitzer resistivity due to fuel ions (ohm m)"""

    e_plasma_alpha_fast_critical_profile: list[float] = field(default_factory=list)
    """Profile of critical energy for fast alpha particles slowing down in plasma [J]"""

    t_plasma_fast_alpha_thermalisation_profile: list[float] = field(default_factory=list)
    """Profile of fast alpha particle average thermalisation time in plasma [s]"""

    nd_plasma_alphas_fast_profile: list[float] = field(default_factory=list)
    """Profile of fast alpha particle density in plasma [/m³]"""

    f_p_plasma_alpha_fast_ions_profile: list[float] = field(default_factory=list)
    """Profile of fast alpha particle energy fraction transferred to ions in plasma (0..1)"""

    f_p_plasma_alpha_fast_electrons_profile: list[float] = field(default_factory=list)
    """Profile of fast alpha particle energy fraction transferred to electrons in plasma (0..1)"""

    dt_power_density_plasma: float = 0.0
    sigmav_dt_average: float = 0.0
    dhe3_power_density: float = 0.0
    dd_power_density: float = 0.0
    fusrat: float = 0.0


CREATE_DICTS_FROM_DATACLASS = PhysicsData
