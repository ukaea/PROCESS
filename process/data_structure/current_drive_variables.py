"""
Module containing global variables relating to the current drive system
"""

dx_beam_duct: float = None
"""width of neutral beam duct where it passes between the TF coils (m)
T Inoue et al, Design of neutral beam system for ITER-FEAT,
<A HREF=http://dx.doi.org/10.1016/S0920-3796(01)00339-8>
Fusion Engineering and Design, Volumes 56-57, October 2001, Pages 517-521</A>)
"""


big_q_plasma: float = None
"""Fusion gain; P_fusion / (P_injection + P_ohmic)"""


f_c_plasma_bootstrap: float = None
"""bootstrap current fraction (enforced; see i_bootstrap_current)"""


f_c_plasma_bootstrap_max: float = None
"""maximum fraction of plasma current from bootstrap; if `f_c_plasma_bootstrap_max < 0`,
bootstrap fraction = abs(f_c_plasma_bootstrap_max)
"""


f_c_plasma_bootstrap_iter89: float = None
"""bootstrap current fraction, ITER 1989 model"""


f_c_plasma_bootstrap_nevins: float = None
"""bootstrap current fraction, Nevins et al model"""


f_c_plasma_bootstrap_sauter: float = None
"""bootstrap current fraction, Sauter et al model"""


f_c_plasma_bootstrap_wilson: float = None
"""bootstrap current fraction, Wilson et al model"""


f_c_plasma_bootstrap_sakai: float = None
"""Bootstrap current fraction, Sakai et al model"""


f_c_plasma_bootstrap_aries: float = None
"""Bootstrap current fraction, ARIES model"""


f_c_plasma_bootstrap_andrade: float = None
"""Bootstrap current fraction, Andrade et al model"""


f_c_plasma_bootstrap_hoang: float = None
"""Bootstrap current fraction, Hoang et al model"""


f_c_plasma_bootstrap_wong: float = None
"""Bootstrap current fraction, Wong et al model"""


bscf_gi_i: float = None
"""Bootstrap current fraction, first Gi et al model"""


bscf_gi_ii: float = None
"""Bootstrap current fraction, second Gi et al model"""


f_c_plasma_bootstrap_sugiyama_l: float = None
"""Bootstrap current fraction, L-mode Sugiyama et al model"""


f_c_plasma_bootstrap_sugiyama_h: float = None
"""Bootstrap current fraction, H-mode Sugiyama et al model"""


cboot: float = None
"""bootstrap current fraction multiplier"""


c_beam_total: float = None
"""neutral beam current (A)"""


f_c_plasma_diamagnetic_hender: float = None
"""diamagnetic current fraction, Hender fit"""


f_c_plasma_diamagnetic_scene: float = None
"""diamagnetic current fraction, SCENE fit"""


f_c_plasma_diamagnetic: float = None
"""diamagnetic current fraction"""


p_hcd_ecrh_injected_total_mw: float = None
"""ECH power (MW)"""


p_ebw_injected_mw: float = None
"""Electron bernstein power (MW)"""


p_hcd_ecrh_electric_mw: float = None
"""ECH wall plug power (MW)"""


p_hcd_ebw_electric_mw: float = None
"""Electron bernstein wall plug power (MW)"""


eta_cd_hcd_primary: float = None
"""Current drive efficiency of primary HCD system (A/W)"""


eta_cd_hcd_secondary: float = None
"""Current drive efficiency of secondary HCD system (A/W)"""


c_hcd_primary_driven: float = None
"""Current in plasma driven by primary HCD system (A)"""


c_hcd_secondary_driven: float = None
"""Current in plasma driven by secondary HCD system (A)"""


f_c_plasma_hcd_primary: float = None
"""Fraction of plasma current driven by primary HCD system"""


f_c_plasma_hcd_secondary: float = None
"""Fraction of plasma current driven by secondary HCD system"""


n_ecrh_harmonic: float = None
"""cyclotron harmonic frequency number, used in cut-off function"""


i_ecrh_wave_mode: int = None
"""Switch for ECRH wave mode :
- =0 O-mode
- =1 X-mode
"""


e_beam_kev: float = None
"""neutral beam energy (keV) (`iteration variable 19`)"""


eta_hcd_primary_injector_wall_plug: float = None
"""auxiliary power wall plug to injector efficiency"""


eta_hcd_secondary_injector_wall_plug: float = None
"""secondary auxiliary power wall plug to injector efficiency"""


eta_ecrh_injector_wall_plug: float = None
"""ECH wall plug to injector efficiency"""


eta_lowhyb_injector_wall_plug: float = None
"""lower hybrid wall plug to injector efficiency"""


eta_icrh_injector_wall_plug: float = None
"""Ion cyclotron wall plug to injector efficiency"""


eta_ebw_injector_wall_plug: float = None
"""Electron bernstein wave wall plug to injector efficiency"""


eta_beam_injector_wall_plug: float = None
"""neutral beam wall plug to injector efficiency"""


f_p_beam_injected_ions: float = None
"""fraction of beam energy to ions"""


p_beam_injected_mw: float = None
"""neutral beam power entering vacuum vessel"""


f_c_plasma_pfirsch_schluter_scene: float = None
"""Pfirsch-Schlüter current fraction, SCENE fit"""


p_beam_shine_through_mw: float = None
"""neutral beam shine-through power"""


feffcd: float = None
"""current drive efficiency fudge factor (`iteration variable 47`)"""


f_p_beam_orbit_loss: float = None
"""fraction of neutral beam power lost after ionisation but before
thermalisation (orbit loss fraction)
"""


f_radius_beam_tangency_rmajor: float = None
"""R_tangential / R_major for neutral beam injection"""


f_beam_tritium: float = None
"""fraction of beam that is tritium"""


eta_cd_norm_hcd_primary: float = None
"""Normalised current drive efficiency for primary HCD system [(1.0e20 A)/(W m^2)]"""

eta_cd_dimensionless_hcd_primary: float = None
"""Dimensionless current drive efficiency for primary HCD system (ζ)"""


eta_cd_norm_hcd_secondary: float = None
"""Normalised current drive efficiency for secondary HCD system [(1.0e20 A)/(W m^2)]"""

eta_cd_dimensionless_hcd_secondary: float = None
"""Dimensionless current drive efficiency for secondary HCD system (ζ)"""


eta_cd_norm_ecrh: float = None
"""User input ECRH gamma (1.0e20 A/(W m^2))"""


xi_ebw: float = None
"""User scaling input for EBW plasma heating. Default 0.43"""


i_hcd_primary: int = None
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


i_hcd_secondary: int = None
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


i_hcd_calculations: int = None
"""Switch for current drive calculation:
- =0 turned off
- =1 turned on
"""


f_p_beam_shine_through: float = None
"""neutral beam shine-through fraction"""


dx_beam_shield: float = None
"""neutral beam duct shielding thickness (m)"""


p_hcd_primary_extra_heat_mw: float = None
"""heating power not used for current drive (MW) (`iteration variable 11`)"""


p_hcd_secondary_extra_heat_mw: float = None
"""secondary fixed heating power not used for current drive (MW)"""


p_hcd_injected_max: float = None
"""maximum allowable value for injected power (MW) (`constraint equation 30`)"""


p_hcd_injected_electrons_mw: float = None
"""auxiliary injected power to electrons (MW)"""


p_hcd_injected_ions_mw: float = None
"""auxiliary injected power to ions (MW)"""


p_hcd_injected_total_mw: float = None
"""total auxiliary injected power (MW)"""


p_hcd_injected_current_total_mw: float = None
"""total auxiliary injected power (MW)"""


p_hcd_secondary_injected_mw: float = None
"""secondary total fixed auxiliary injected power (MW)"""


p_hcd_primary_injected_mw: float = None
"""primary auxiliary injected power (MW)"""


f_c_plasma_internal: float = None
"""plasma current fraction driven internally (Bootstrap + Diamagnetic + PS)"""


p_hcd_lowhyb_injected_total_mw: float = None
"""Total lower hybrid injection power (MW)"""


p_hcd_icrh_injected_total_mw: float = None
"""Total ion cyclotron injection power (MW)"""


p_hcd_ebw_injected_total_mw: float = None
"""Total electron bernstein wave injection power (MW)"""


p_beam_plasma_coupled_mw: float = None
"""Total neutral beam power that is coupled to plasma after losses (MW)"""


p_hcd_beam_injected_total_mw: float = None
"""neutral beam injection power (MW)"""


p_beam_orbit_loss_mw: float = None
"""neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)"""


f_c_plasma_pfirsch_schluter: float = None
"""Pfirsch-Schlüter current fraction"""


p_hcd_lowhyb_electric_mw: float = None
"""lower hybrid wall plug power (MW)"""


pwpnb: float = None
"""neutral beam wall plug power (MW)"""


radius_beam_tangency: float = None
"""neutral beam centreline tangency radius (m)"""


radius_beam_tangency_max: float = None
"""maximum tangency radius for centreline of beam (m)"""


n_beam_decay_lengths_core: float = None
"""neutral beam e-decay lengths to plasma centre"""


n_beam_decay_lengths_core_required: float = None
"""permitted neutral beam e-decay lengths to plasma centre"""


def init_current_drive_variables():
    """Initialise current drive variables"""
    global \
        dx_beam_duct, \
        big_q_plasma, \
        f_c_plasma_bootstrap, \
        f_c_plasma_bootstrap_max, \
        f_c_plasma_bootstrap_iter89, \
        f_c_plasma_bootstrap_nevins, \
        f_c_plasma_bootstrap_sauter, \
        f_c_plasma_bootstrap_wilson, \
        f_c_plasma_bootstrap_sakai, \
        f_c_plasma_bootstrap_aries, \
        f_c_plasma_bootstrap_andrade, \
        f_c_plasma_bootstrap_hoang, \
        f_c_plasma_bootstrap_wong, \
        bscf_gi_i, \
        bscf_gi_ii, \
        f_c_plasma_bootstrap_sugiyama_l, \
        f_c_plasma_bootstrap_sugiyama_h, \
        cboot, \
        c_beam_total, \
        f_c_plasma_diamagnetic_hender, \
        f_c_plasma_diamagnetic_scene, \
        f_c_plasma_diamagnetic, \
        p_hcd_ecrh_injected_total_mw, \
        p_ebw_injected_mw, \
        p_hcd_ecrh_electric_mw, \
        p_hcd_ebw_electric_mw, \
        eta_cd_hcd_primary, \
        eta_cd_hcd_secondary, \
        c_hcd_primary_driven, \
        c_hcd_secondary_driven, \
        f_c_plasma_hcd_primary, \
        f_c_plasma_hcd_secondary, \
        n_ecrh_harmonic, \
        i_ecrh_wave_mode, \
        e_beam_kev, \
        eta_hcd_primary_injector_wall_plug, \
        eta_hcd_secondary_injector_wall_plug, \
        eta_ecrh_injector_wall_plug, \
        eta_lowhyb_injector_wall_plug, \
        eta_icrh_injector_wall_plug, \
        eta_ebw_injector_wall_plug, \
        eta_beam_injector_wall_plug, \
        f_p_beam_injected_ions, \
        p_beam_injected_mw, \
        f_c_plasma_pfirsch_schluter_scene, \
        p_beam_shine_through_mw, \
        feffcd, \
        f_p_beam_orbit_loss, \
        f_radius_beam_tangency_rmajor, \
        f_beam_tritium, \
        eta_cd_norm_hcd_primary, \
        eta_cd_dimensionless_hcd_primary, \
        eta_cd_norm_hcd_secondary, \
        eta_cd_dimensionless_hcd_secondary, \
        eta_cd_norm_ecrh, \
        xi_ebw, \
        i_hcd_primary, \
        i_hcd_secondary, \
        i_hcd_calculations, \
        f_p_beam_shine_through, \
        dx_beam_shield, \
        p_hcd_primary_extra_heat_mw, \
        p_hcd_secondary_extra_heat_mw, \
        p_hcd_injected_max, \
        p_hcd_injected_electrons_mw, \
        p_hcd_injected_ions_mw, \
        p_hcd_injected_total_mw, \
        p_hcd_injected_current_total_mw, \
        p_hcd_secondary_injected_mw, \
        p_hcd_primary_injected_mw, \
        f_c_plasma_internal, \
        p_hcd_lowhyb_injected_total_mw, \
        p_hcd_icrh_injected_total_mw, \
        p_hcd_ebw_injected_total_mw, \
        p_beam_plasma_coupled_mw, \
        p_hcd_beam_injected_total_mw, \
        p_beam_orbit_loss_mw, \
        f_c_plasma_pfirsch_schluter, \
        p_hcd_lowhyb_electric_mw, \
        pwpnb, \
        radius_beam_tangency, \
        radius_beam_tangency_max, \
        n_beam_decay_lengths_core, \
        n_beam_decay_lengths_core_required

    dx_beam_duct = 0.58
    big_q_plasma = 0.0
    f_c_plasma_bootstrap = 0.0
    f_c_plasma_bootstrap_max = 0.9
    f_c_plasma_bootstrap_iter89 = 0.0
    f_c_plasma_bootstrap_nevins = 0.0
    f_c_plasma_bootstrap_sauter = 0.0
    f_c_plasma_bootstrap_wilson = 0.0
    f_c_plasma_bootstrap_sakai = 0.0
    f_c_plasma_bootstrap_aries = 0.0
    f_c_plasma_bootstrap_andrade = 0.0
    f_c_plasma_bootstrap_hoang = 0.0
    f_c_plasma_bootstrap_wong = 0.0
    bscf_gi_i = 0.0
    bscf_gi_ii = 0.0
    f_c_plasma_bootstrap_sugiyama_l = 0.0
    f_c_plasma_bootstrap_sugiyama_h = 0.0
    cboot = 1.0
    c_beam_total = 0.0
    f_c_plasma_diamagnetic_hender = 0.0
    f_c_plasma_diamagnetic_scene = 0.0
    f_c_plasma_diamagnetic = 0.0
    p_hcd_ecrh_injected_total_mw = 0.0
    eta_cd_hcd_primary = 0.0
    n_ecrh_harmonic = 2.0
    i_ecrh_wave_mode = 0
    e_beam_kev = 1.0e3
    eta_hcd_primary_injector_wall_plug = 0.3
    eta_hcd_secondary_injector_wall_plug = 0.3
    eta_ecrh_injector_wall_plug = 0.3
    eta_lowhyb_injector_wall_plug = 0.3
    eta_beam_injector_wall_plug = 0.3
    f_p_beam_injected_ions = 0.5
    p_beam_injected_mw = 0.0
    f_c_plasma_pfirsch_schluter_scene = 0.0
    p_beam_shine_through_mw = 0.0
    feffcd = 1.0
    f_p_beam_orbit_loss = 0.0
    f_radius_beam_tangency_rmajor = 1.05
    f_beam_tritium = 1e-6
    eta_cd_norm_hcd_primary = 0.0
    eta_cd_dimensionless_hcd_primary = 0.0
    eta_cd_norm_ecrh = 0.35
    xi_ebw = 0.8
    i_hcd_primary = 5
    i_hcd_secondary = 0
    i_hcd_calculations = 1
    f_p_beam_shine_through = 0.0
    dx_beam_shield = 0.5
    p_hcd_primary_extra_heat_mw = 0.0
    p_hcd_secondary_extra_heat_mw = 0.0
    p_hcd_injected_max = 150.0
    p_hcd_injected_electrons_mw = 0.0
    p_hcd_injected_ions_mw = 0.0
    p_hcd_injected_total_mw = 0.0
    p_hcd_injected_current_total_mw = 0.0
    p_hcd_secondary_injected_mw = 0.0
    f_c_plasma_internal = 0.0
    p_beam_orbit_loss_mw = 0.0
    f_c_plasma_pfirsch_schluter = 0.0
    pwpnb = 0.0
    radius_beam_tangency = 0.0
    radius_beam_tangency_max = 0.0
    n_beam_decay_lengths_core = 0.0
    n_beam_decay_lengths_core_required = 3.0
    eta_cd_norm_hcd_secondary = 0.0
    eta_cd_dimensionless_hcd_secondary = 0.0
    eta_cd_hcd_secondary = 0.0
    p_ebw_injected_mw = 0.0
    p_hcd_ecrh_electric_mw = 0.0
    p_hcd_ebw_electric_mw = 0.0
    c_hcd_primary_driven = 0.0
    c_hcd_secondary_driven = 0.0
    f_c_plasma_hcd_primary = 0.0
    f_c_plasma_hcd_secondary = 0.0
    eta_icrh_injector_wall_plug = 0.3
    eta_ebw_injector_wall_plug = 0.3
    p_hcd_primary_injected_mw = 0.0
    p_hcd_icrh_injected_total_mw = 0.0
    p_hcd_ebw_injected_total_mw = 0.0
    p_hcd_beam_injected_total_mw = 0.0
    p_hcd_icrh_injected_total_mw = 0.0
    p_hcd_lowhyb_electric_mw = 0.0
    p_hcd_lowhyb_injected_total_mw = 0.0
