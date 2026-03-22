"""Module containing tokamak plasma physics routines

This module contains all the primary plasma physics routines
for a tokamak device.


Module containing global variables relating to the plasma physics
"""

from dataclasses import dataclass, field
from enum import IntEnum

import numpy as np


class DivertorNumberModels(IntEnum):
    """Enum for divertor number models. `i_single_null` is the index for this enum."""

    DOUBLE_NULL = 0
    SINGLE_NULL = 1


N_CONFINEMENT_SCALINGS = 51
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
    """[W]"""

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
    """current profile index"""

    alphaj_wesson: float = None
    """Wesson-like current profile index"""

    alphan: float = 0.25
    """density profile index"""

    alphap: float = 0.0
    """pressure profile index"""

    fusden_alpha_total: float = 0.0
    """Alpha particle production rate per unit volume, from plasma and beams [particles/m3/sec]"""

    fusden_plasma_alpha: float = 0.0
    """Alpha particle production rate per unit volume, just from plasma [particles/m3/sec]"""

    alphat: float = 0.5
    """temperature profile index"""

    aspect: float = 2.907
    """aspect ratio (`iteration variable 1`)"""

    beamfus0: float = 1.0
    """multiplier for beam-background fusion calculation"""

    beta_total_vol_avg: float = 0.042
    """Volume averaged total plasma beta (`iteration variable 5`) (calculated if stellarator)"""

    beta_fast_alpha: float = 0.0
    """fast alpha beta component"""

    beta_vol_avg_max: float = 0.0
    """Max allowable volume averaged beta"""

    beta_vol_avg_min: float = 0.0
    """Minimum allowable volume averaged beta"""

    beta_beam: float = 0.0
    """neutral beam beta component"""

    beta_poloidal_vol_avg: float = 0.0
    """poloidal beta"""

    beta_poloidal_eps: float = 0.0
    """Poloidal beta and inverse aspcet ratio product"""

    beta_toroidal_vol_avg: float = 0.0
    """Plasma volume averaged toroidal beta"""

    beta_thermal_toroidal_profile: list[float] = field(default_factory=list)
    """toroidal beta profile"""

    beta_thermal_vol_avg: float = 0.0
    """Plasma volume averaged thermal beta"""

    beta_thermal_poloidal_vol_avg: float = 0.0
    """Plasma volume averaged poloidal thermal beta"""

    beta_thermal_toroidal_vol_avg: float = 0.0
    """Plasma volume averaged toloidal thermal beta"""

    beta_norm_total: float = 0.0
    """normaised total beta"""

    beta_norm_thermal: float = 0.0
    """normaised thermal beta"""

    beta_norm_toroidal: float = 0.0
    """normaised toroidal beta"""

    beta_norm_poloidal: float = 0.0
    """normaised poloidal beta"""

    e_plasma_beta_thermal: float = 0.0
    """Plasma thermal energy derived from thermal beta"""

    betbm0: float = 1.5
    """leading coefficient for NB beta fraction"""

    b_plasma_surface_poloidal_average: float = 0.0
    """Plasma surface average poloidal field (T)"""

    b_plasma_toroidal_on_axis: float = 5.68
    """Plasma toroidal field on axis (T) (`iteration variable 2`)"""

    b_plasma_inboard_toroidal: float = 0.0
    """Plasma inboard toroidal field (T)"""

    b_plasma_outboard_toroidal: float = 0.0
    """Plasma outboard toroidal field (T)"""

    b_plasma_toroidal_profile: list[float] = field(default_factory=list)
    """toroidal field profile in plasma (T)"""

    b_plasma_total: float = 0.0
    """Sum of plasma total toroidal + poloidal field (T)"""

    e_plasma_magnetic_stored: float = 0.0
    """Plasma stored magnetic energy (J)"""

    f_plasma_fuel_burnup: float = 0.0
    """fractional plasma burnup"""

    burnup_in: float = 0.0
    """fractional plasma burnup user input"""

    b_plasma_vertical_required: float = 0.0
    """Vertical field needed for plasma equilibrium (T)"""

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
    """Plasma volume averaged electron density (/m3) (`iteration variable 6`)"""

    nd_plasma_fuel_ions_vol_avg: float = 0.0
    """Plasma volume averaged fuel ion density (/m3)"""

    dlamee: float = 0.0
    """electron-electron coulomb logarithm"""

    dlamie: float = 0.0
    """ion-electron coulomb logarithm"""

    nd_plasma_electron_max_array: list[float] = field(
        default_factory=lambda: np.zeros(8, dtype=np.float64)
    )
    """Array of plasma electron density upper limits values (/m3)"""

    nd_plasma_alphas_vol_avg: float = 0.0
    """Plasma volume averaged thermal alpha density (/m3)"""

    nd_beam_ions: float = 0.0
    """hot beam ion density, variable (/m3)"""

    nd_beam_ions_out: float = 0.0
    """hot beam ion density from calculation (/m3)"""

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
    """Plasma electron max density limit (/m3)"""

    nd_plasma_ions_total_vol_avg: float = 0.0
    """Plasma volume averaged total ion density (/m3)"""

    nd_plasma_electron_line: float = 0.0
    """Plasma line averaged electron density (/m3)"""

    nd_plasma_protons_vol_avg: float = 0.0
    """Plasma volume averaged proton ash density (/m3)"""

    ntau: float = 0.0
    """Fusion double product (s/m3)"""

    nTtau: float = 0.0
    """Lawson triple product [keV s / m3]"""

    nd_plasma_impurities_vol_avg: float = 0.0
    """Plasma volume averaged impurity (Z > 2) ion density (/m3)"""

    gradient_length_ne: float = None
    """Max. normalised gradient length in el. density (i_plasma_pedestal==0 only)"""

    gradient_length_te: float = None
    """Max. normalised gradient length in el. temperature (i_plasma_pedestal==0 only)"""

    beta_poloidal_eps_max: float = 1.38
    """maximum (eps*beta_poloidal) (`constraint equation 6`). Note: revised issue #346
    "Operation at the tokamak equilibrium poloidal beta-limit in TFTR", 1992 Nucl. Fusion 32 1468
    """

    eps: float = 0.34399724802
    """inverse aspect ratio"""

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

    fusrat_plasma_dt: float = None
    """ D-T fusion reaction rate in plasma, (reactions/sec)"""

    fusrat_dt_total: float = None
    """ Total D-T fusion reaction rate from beams and plasma, (reactions/sec)"""

    fusrat_plasma_dd_helion: float = None
    """D-D fusion reaction rate (helium branch) in plasma, (reactions/sec)"""

    fusrat_plasma_dd_triton: float = None
    """D-D fusion reaction rate (tritium branch) in plasma, (reactions/sec)"""

fusrat_plasma_dd_total: float = None
"""Total D-D fusion reaction rate in plasma, (reactions/sec)"""

fusrat_plasma_dhe3: float = None
"""D-3He fusion reaction rate in plasma, (reactions/sec)"""

fusrat_neutron_production_total: float = None
"""Total neutron production rate from plasma and beams (neutrons/sec)"""


fusrat_plasma_dhe3: float = None
"""D-3He fusion reaction rate in plasma, (reactions/sec)"""

fusrat_neutron_production_total: float = None
"""Total neutron production rate from plasma and beams (neutrons/sec)"""

    fusrat_plasma_dd_helion: float = None
    """D-D fusion reaction rate (helium branch) in plasma, (reactions/sec)"""

    fusrat_plasma_dd_triton: float = None
    """D-D fusion reaction rate (tritium branch) in plasma, (reactions/sec)"""

    
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
    """electron density of pedestal [m⁻³] (`i_plasma_pedestal==1)"""

    nd_plasma_separatrix_electron: float = 3.0e19
    """electron density at separatrix [m⁻³] (`i_plasma_pedestal==1)"""

    i_nd_plasma_pedestal_separatrix: int = 1
    """switch for pedestal and separatrix density calculation:
    - =0 User input pedestal and separatrix density
    - =1 Calculate pedestal and separatrix density as fraction of Greenwald limit (see `f_nd_plasma_pedestal_greenwald` and `f_nd_plasma_separatrix_greenwald`)
    """

    alpha_crit: float = 0.0
    """critical ballooning parameter value"""

    nd_plasma_separatrix_electron_eich_max: float = 0.0
    """Eich critical electron density at separatrix [m⁻³]"""

    plasma_res_factor: float = 1.0
    """plasma resistivity pre-factor"""

    radius_plasma_pedestal_density_norm: float = 1.0
    """Normalised radius of density pedestal (`i_plasma_pedestal==1`)"""

    radius_plasma_pedestal_temp_norm: float = 1.0
    """Normalised radius of temperature pedestal (`i_plasma_pedestal==1`)"""

    rho_te_max: float = 0.0
    """r/a where the temperature gradient is largest (`i_plasma_pedestal==0`)"""

    rho_ne_max: float = 0.0
    """r/a where the density gradient is largest (`i_plasma_pedestal==0`)"""

    tbeta: float = 2.0
    """temperature profile index beta  (`i_plasma_pedestal==1)"""

    temp_plasma_pedestal_kev: float = 1.0
    """Plasma electron temperature of pedestal (keV) (`i_plasma_pedestal==1`)"""

    temp_plasma_separatrix_kev: float = 0.1
    """Plasma electron temperature at separatrix (keV) (`i_plasma_pedestal==1`) calculated if reinke
    criterion is used (`icc=78`)
    """

    i_beta_norm_max: int = 1
    """Switch for maximum normalised beta scaling:"""

    i_ind_plasma_internal_norm: int = 0
    """Switch for plasma normalised internal inductance scaling:"""

    i_alphaj: int = 0
    """Switch for current profile index scaling:"""

    i_rad_loss: int = 1
    """switch for radiation loss term usage in power balance (see User Guide):
    - =0 total power lost is scaling power plus radiation
    - =1 total power lost is scaling power plus core radiation only
    - =2 total power lost is scaling power only, with no additional
    allowance for radiation. This is not recommended for power plant models.
    """

    i_confinement_time: int = 34
    """switch for energy confinement time scaling law"""

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
    """plasma squareness used by Sauter plasma shape"""

    kappa: float = 1.792
    """plasma separatrix elongation (calculated if `i_plasma_geometry = 1-5, 7 or 9-10`)"""

    kappa95: float = 1.6
    """plasma elongation at 95% surface (calculated if `i_plasma_geometry = 0-3, 6, or 8-10`)"""

    kappa_ipb: float = 0.0
    """Separatrix elongation calculated for IPB scalings"""

    nd_plasma_electron_on_axis: float = 0.0
    """central electron density (/m3)"""

    nd_plasma_ions_on_axis: float = 0.0
    """central ion density (/m3)"""

    m_s_limit: float = 0.3
    """margin to vertical stability"""

    pres_plasma_thermal_on_axis: float = 0.0
    """Plasma central thermal pressure (no fast ions or beam pressure) (Pa)"""

    pres_plasma_thermal_total_profile: list[float] = field(default_factory=list)
    """Profile of total pressure in plasma (Pa)"""

    pres_plasma_electron_profile: list[float] = field(default_factory=list)
    """Profile of electron pressure in plasma (Pa)"""

    pres_plasma_ion_total_profile: list[float] = field(default_factory=list)
    """Profile of ion pressure in plasma (Pa)"""

    pres_plasma_fuel_profile: list[float] = field(default_factory=list)
    """Profile of fuel pressure in plasma (Pa)"""

    j_plasma_on_axis: float = 0.0
    """Central plasma current density (A/m2)"""

    j_plasma_bootstrap_sauter_profile: list[float] = field(default_factory=list)
    """Profile of bootstrap current density in plasma using Sauter et al scaling (A/m2)"""

    n_plasma_profile_elements: int = 501
    """Number of elements in plasma profile"""

    pres_plasma_thermal_vol_avg: float = None
    """Volume averaged thermal plasma pressure (no fast ions or beam pressure) (Pa)"""

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
    """EU DEMO divertor protection parameter (MW/T/m)"""

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
    """plasma current (A)"""

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
    """synchrotron radiation power per volume (MW/m3)"""

    p_plasma_sync_mw: float = 0.0
    """Total synchrotron radiation power from plasma (MW)"""

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
    """Safety factor on axis"""

    q95: float = 0.0
    """Safety factor at 95% flux surface (iteration variable 18) (unless icurr=2 (ST current scaling),
    in which case q95 = mean edge safety factor qbar)
    """

    molflow_plasma_fuelling_required: float = 0.0
    """plasma fuelling rate (nucleus-pairs/s)"""

f_plasma_particles_lcfs_recycled: float = None
"""fraction of plasma particles that are recycled at the LCFS. Recycling coefficent (R)"""

eta_plasma_fuelling: float = None
"""fuelling efficiency (fraction of fuel particles injected that become confined in the plasma)"""

molflow_plasma_fuelling_vv_injected: float = None
"""plasma fuelling rate into the vacuum vessel (nucleus-pairs/s)"""

f_molflow_plasma_fuelling_deuterium: float = None
"""fraction of plasma fuelling that is deuterium"""

f_molflow_plasma_fuelling_tritium: float = None
"""fraction of plasma fuelling that is tritium"""

f_molflow_plasma_fuelling_helium3: float = None
"""fraction of plasma fuelling that is helium-3"""
    tauratio: float = 1.0
    """tauratio /1.0/ : ratio of He and pellet particle confinement times"""

    q95_min: float = 0.0
    """lower limit for edge safety factor"""

    qstar: float = 0.0
    """cylindrical safety factor"""

    rad_fraction_sol: float = 0.8
    """SoL radiation fraction"""

    rad_fraction_total: float = 0.0
    """Radiation fraction total = SoL + LCFS radiation / total power deposited in plasma"""

    f_nd_alpha_electron: float = 0.1
    """thermal alpha density/electron density (`iteration variable 109`)"""

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
    """plasma major radius (m) (`iteration variable 3`)"""

    rminor: float = 0.0
    """plasma minor radius (m)"""

    f_nd_beam_electron: float = 0.005
    """hot beam density / n_e (`iteration variable 7`)"""


f_nd_plasma_carbon_electron: float = None
"""n_carbon / n_e"""

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

    f_alpha_energy_confinement: float = 0.0
    """alpha particle to energy confinement time ratio"""

    temp_plasma_electron_vol_avg_kev: float = 12.9
    """volume averaged electron temperature (keV) (`iteration variable 4`)"""

    temp_plasma_electron_on_axis_kev: float = 0.0
    """central electron temperature (keV)"""

    temp_plasma_electron_density_weighted_kev: float = 0.0
    """density weighted average electron temperature (keV)"""

    temp_plasma_electron_line_avg_kev: float = None
    """Line averaged electron temperature (keV)"""

    temp_plasma_ion_vol_avg_kev: float = 12.9
    """volume averaged ion temperature (keV). N.B. calculated from temp_plasma_electron_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`"""

    temp_plasma_ion_on_axis_kev: float = 0.0
    """central ion temperature (keV)"""

    temp_plasma_ion_density_weighted_kev: float = 0.0
    """density weighted average ion temperature (keV)"""

    f_temp_plasma_ion_electron: float = 1.0
    """ion temperature / electron temperature(used to calculate temp_plasma_ion_vol_avg_kev if `f_temp_plasma_ion_electron > 0.0`"""

    triang: float = 0.36
    """plasma separatrix triangularity (calculated if `i_plasma_geometry = 1, 3-5 or 7`)"""

    triang95: float = 0.24
    """plasma triangularity at 95% surface (calculated if `i_plasma_geometry = 0-2, 6, 8 or 9`)"""

    vol_plasma: float = 0.0
    """plasma volume (m3)"""

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
    """total V-s needed (Wb)"""

    pflux_fw_neutron_mw: float = 0.0
    """average neutron wall load (MW/m2)"""

    pflux_plasma_surface_neutron_avg_mw: float = 0.0
    """Average neutron flux at plasma surface (MW/m2)"""

    wtgpd: float = 0.0
    """mass of fuel used per day (g)"""

    a_plasma_poloidal: float = 0.0
    """plasma poloidal cross-sectional area [m^2]"""

    n_charge_plasma_effective_vol_avg: float = 0.0
    """Volume averaged plasma effective charge"""

    n_charge_plasma_effective_profile: list[float] = field(default_factory=list)
    """Profile of plasma effective charge"""

    n_charge_plasma_effective_mass_weighted_vol_avg: float = 0.0
    """Plasma mass-weighted volume averaged plasma effective charge"""

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

    dt_power_density_plasma: float = 0.0
    sigmav_dt_average: float = 0.0
    dhe3_power_density: float = 0.0
    dd_power_density: float = 0.0
    fusrat: float = 0.0


CREATE_DICTS_FROM_DATACLASS = PhysicsData
