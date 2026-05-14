"""
Module containing global variables relating to the current drive system
"""

from dataclasses import dataclass


@dataclass
class CurrentDriveData:
    dx_beam_duct: float = 0.58
    """width of neutral beam duct where it passes between the TF coils (m)
    T Inoue et al, Design of neutral beam system for ITER-FEAT,
    <A HREF=http://dx.doi.org/10.1016/S0920-3796(01)00339-8>
    Fusion Engineering and Design, Volumes 56-57, October 2001, Pages 517-521</A>)
    """

    big_q_plasma: float = 0.0
    """Fusion gain; P_fusion / (P_injection + P_ohmic)"""

    f_c_plasma_bootstrap: float = 0.0
    """bootstrap current fraction (enforced; see i_bootstrap_current)"""

    f_c_plasma_bootstrap_max: float = 0.9
    """maximum fraction of plasma current from bootstrap; if `f_c_plasma_bootstrap_max < 0`,
    bootstrap fraction = abs(f_c_plasma_bootstrap_max)
    """

    f_c_plasma_bootstrap_iter89: float = 0.0
    """bootstrap current fraction, ITER 1989 model"""

    f_c_plasma_bootstrap_nevins: float = 0.0
    """bootstrap current fraction, Nevins et al model"""

    f_c_plasma_bootstrap_sauter: float = 0.0
    """bootstrap current fraction, Sauter et al model"""

    f_c_plasma_bootstrap_wilson: float = 0.0
    """bootstrap current fraction, Wilson et al model"""

    f_c_plasma_bootstrap_sakai: float = 0.0
    """Bootstrap current fraction, Sakai et al model"""

    f_c_plasma_bootstrap_aries: float = 0.0
    """Bootstrap current fraction, ARIES model"""

    f_c_plasma_bootstrap_andrade: float = 0.0
    """Bootstrap current fraction, Andrade et al model"""

    f_c_plasma_bootstrap_hoang: float = 0.0
    """Bootstrap current fraction, Hoang et al model"""

    f_c_plasma_bootstrap_wong: float = 0.0
    """Bootstrap current fraction, Wong et al model"""

    bscf_gi_i: float = 0.0
    """Bootstrap current fraction, first Gi et al model"""

    bscf_gi_ii: float = 0.0
    """Bootstrap current fraction, second Gi et al model"""

    f_c_plasma_bootstrap_sugiyama_l: float = 0.0
    """Bootstrap current fraction, L-mode Sugiyama et al model"""

    f_c_plasma_bootstrap_sugiyama_h: float = 0.0
    """Bootstrap current fraction, H-mode Sugiyama et al model"""

    cboot: float = 1.0
    """bootstrap current fraction multiplier"""

    c_beam_total: float = 0.0
    """neutral beam current (A)"""

    f_c_plasma_diamagnetic_hender: float = 0.0
    """diamagnetic current fraction, Hender fit"""

    f_c_plasma_diamagnetic_scene: float = 0.0
    """diamagnetic current fraction, SCENE fit"""

    f_c_plasma_diamagnetic: float = 0.0
    """diamagnetic current fraction"""

    p_hcd_ecrh_injected_total_mw: float = 0.0
    """ECH power (MW)"""

    p_ebw_injected_mw: float = 0.0
    """Electron bernstein power (MW)"""

    p_hcd_ecrh_electric_mw: float = 0.0
    """ECH wall plug power (MW)"""

    p_hcd_ebw_electric_mw: float = 0.0
    """Electron bernstein wall plug power (MW)"""

    p_hcd_icrh_electric_mw: float = 0.0
    """Ion cyclotron wall plug power (MW)"""

    eta_cd_hcd_primary: float = 0.0
    """Current drive efficiency of primary HCD system (A/W)"""

    eta_cd_hcd_secondary: float = 0.0
    """Current drive efficiency of secondary HCD system (A/W)"""

    c_hcd_primary_driven: float = 0.0
    """Current in plasma driven by primary HCD system (A)"""

    c_hcd_secondary_driven: float = 0.0
    """Current in plasma driven by secondary HCD system (A)"""

    f_c_plasma_hcd_primary: float = 0.0
    """Fraction of plasma current driven by primary HCD system"""

    f_c_plasma_hcd_secondary: float = 0.0
    """Fraction of plasma current driven by secondary HCD system"""

    n_ecrh_harmonic: float = 2.0
    """cyclotron harmonic frequency number, used in cut-off function"""

    i_ecrh_wave_mode: int = 0
    """Switch for ECRH wave mode :
    - =0 O-mode
    - =1 X-mode
    """

    e_beam_kev: float = 1.0e3
    """neutral beam energy (keV) (`iteration variable 19`)"""

    eta_hcd_primary_injector_wall_plug: float = 0.3
    """auxiliary power wall plug to injector efficiency"""

    eta_hcd_secondary_injector_wall_plug: float = 0.3
    """secondary auxiliary power wall plug to injector efficiency"""

    eta_ecrh_injector_wall_plug: float = 0.3
    """ECH wall plug to injector efficiency"""

    eta_lowhyb_injector_wall_plug: float = 0.3
    """lower hybrid wall plug to injector efficiency"""

    eta_icrh_injector_wall_plug: float = 0.3
    """Ion cyclotron wall plug to injector efficiency"""

    eta_ebw_injector_wall_plug: float = 0.3
    """Electron bernstein wave wall plug to injector efficiency"""

    eta_beam_injector_wall_plug: float = 0.3
    """neutral beam wall plug to injector efficiency"""

    f_p_beam_injected_ions: float = 0.5
    """fraction of beam energy to ions"""

    p_beam_injected_mw: float = 0.0
    """neutral beam power entering vacuum vessel"""

    f_c_plasma_pfirsch_schluter_scene: float = 0.0
    """Pfirsch-Schlüter current fraction, SCENE fit"""

    p_beam_shine_through_mw: float = 0.0
    """neutral beam shine-through power"""

    feffcd: float = 1.0
    """current drive efficiency fudge factor (`iteration variable 47`)"""

    f_p_beam_orbit_loss: float = 0.0
    """fraction of neutral beam power lost after ionisation but before
    thermalisation (orbit loss fraction)
    """

    f_radius_beam_tangency_rmajor: float = 1.05
    """R_tangential / R_major for neutral beam injection"""

    f_beam_tritium: float = 1e-6
    """fraction of beam that is tritium"""

    eta_cd_norm_hcd_primary: float = 0.0
    """Normalised current drive efficiency for primary HCD system [(1.0e20 A)/(W m^2)]"""

    eta_cd_dimensionless_hcd_primary: float = 0.0
    """Dimensionless current drive efficiency for primary HCD system (ζ)"""

    eta_cd_norm_hcd_secondary: float = 0.0
    """Normalised current drive efficiency for secondary HCD system [(1.0e20 A)/(W m^2)]"""

    eta_cd_dimensionless_hcd_secondary: float = 0.0
    """Dimensionless current drive efficiency for secondary HCD system (ζ)"""

    eta_cd_norm_ecrh: float = 0.35
    """User input ECRH gamma (1.0e20 A/(W m^2))"""

    xi_ebw: float = 0.8
    """User scaling input for EBW plasma heating. Default 0.43"""

    i_hcd_primary: int = 5
    """Switch for current drive efficiency model:
    - =1 Fenstermacher Lower Hybrid
    - =2 Ion Cyclotron current drive
    - =3 Fenstermacher ECH
    - =4 Ehst Lower Hybrid
    - =5 ITER Neutral Beam
    - =6 new Culham Lower Hybrid model
    - =7 new Culham ECCD model
    - =8 new Culham Neutral Beam model
    - =9 RFP option removed in PROCESS (issue #508)
    - =10 ECRH user input gamma
    - =11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.
    - =12 EBW user scaling input. Scaling (S. Freethy)
    """

    i_hcd_secondary: int = 0
    """Switch for 2nd current drive efficiency model:
    - =0 No fixed current drive
    - =1 Fenstermacher Lower Hybrid
    - =2 Ion Cyclotron current drive
    - =3 Fenstermacher ECH
    - =4 Ehst Lower Hybrid
    - =5 ITER Neutral Beam
    - =6 new Culham Lower Hybrid model
    - =7 new Culham ECCD model
    - =8 new Culham Neutral Beam model
    - =9 RFP option removed in PROCESS (issue #508)
    - =10 ECRH user input gamma
    - =11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.
    - =12 EBW user scaling input. Scaling (S. Freethy)
    """

    i_hcd_calculations: int = 1
    """Switch for current drive calculation:
    - =0 turned off
    - =1 turned on
    """

    f_p_beam_shine_through: float = 0.0
    """neutral beam shine-through fraction"""

    dx_beam_shield: float = 0.5
    """neutral beam duct shielding thickness (m)"""

    p_hcd_primary_extra_heat_mw: float = 0.0
    """heating power not used for current drive (MW) (`iteration variable 11`)"""

    p_hcd_secondary_extra_heat_mw: float = 0.0
    """secondary fixed heating power not used for current drive (MW)"""

    p_hcd_injected_max: float = 150.0
    """maximum allowable value for injected power (MW) (`constraint equation 30`)"""

    p_hcd_injected_electrons_mw: float = 0.0
    """auxiliary injected power to electrons (MW)"""

    p_hcd_injected_ions_mw: float = 0.0
    """auxiliary injected power to ions (MW)"""

    p_hcd_injected_total_mw: float = 0.0
    """total auxiliary injected power (MW)"""

    p_hcd_injected_current_total_mw: float = 0.0
    """total auxiliary injected power (MW)"""

    p_hcd_secondary_injected_mw: float = 0.0
    """secondary total fixed auxiliary injected power (MW)"""

    p_hcd_primary_injected_mw: float = 0.0
    """primary auxiliary injected power (MW)"""

    f_c_plasma_internal: float = 0.0
    """plasma current fraction driven internally (Bootstrap + Diamagnetic + PS)"""

    p_hcd_lowhyb_injected_total_mw: float = 0.0
    """Total lower hybrid injection power (MW)"""

    p_hcd_icrh_injected_total_mw: float = 0.0
    """Total ion cyclotron injection power (MW)"""

    p_hcd_ebw_injected_total_mw: float = 0.0
    """Total electron bernstein wave injection power (MW)"""

    p_beam_plasma_coupled_mw: float = None
    """Total neutral beam power that is coupled to plasma after losses (MW)"""

    p_hcd_beam_injected_total_mw: float = 0.0
    """neutral beam injection power (MW)"""

    p_beam_orbit_loss_mw: float = 0.0
    """neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)"""

    f_c_plasma_pfirsch_schluter: float = 0.0
    """Pfirsch-Schlüter current fraction"""

    p_hcd_lowhyb_electric_mw: float = 0.0
    """lower hybrid wall plug power (MW)"""

    pwpnb: float = 0.0
    """neutral beam wall plug power (MW)"""

    radius_beam_tangency: float = 0.0
    """neutral beam centreline tangency radius (m)"""

    radius_beam_tangency_max: float = 0.0
    """maximum tangency radius for centreline of beam (m)"""

    n_beam_decay_lengths_core: float = 0.0
    """neutral beam e-decay lengths to plasma centre"""

    n_beam_decay_lengths_core_required: float = 3.0
    """permitted neutral beam e-decay lengths to plasma centre"""


CREATE_DICTS_FROM_DATACLASS = CurrentDriveData
