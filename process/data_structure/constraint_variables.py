p_hcd_injected_min_mw: float = None
"""minimum auxiliary power (MW) (`constraint equation 40`)"""


beta_poloidal_max: float = None
"""maximum poloidal beta (`constraint equation 48`)"""


big_q_plasma_min: float = None
"""minimum fusion gain Q (`constraint equation 28`)"""


b_tf_inboard_max: float = None
"""maximum peak toroidal field (T) (`constraint equation 25`)"""


fdene: float = None
"""Scaling value for density limit (constraint equation 5)"""

fradpwr: float = None
"""Scaling value for radiation power upper limit (constraint equation 17)"""


fiooic: float = None
"""Constraint margin for TF coil operating current / critical current ratio
(`constraint equation 33`)
"""

q95_fixed: float = None
"""fixed safety factor q at 95% flux surface
(`constraint equation 68`)
"""

fjohc: float = None
"""Constraint margin for central solenoid current at end-of-flattop
(`constraint equation 26`)
"""


fjohc0: float = None
"""Constraint margin for central solenoid current at beginning of pulse
(`constraint equation 27`)
"""

eta_cd_norm_hcd_primary_max: float = None
"""maximum current drive gamma (`constraint equation 37`)"""


i_q95_fixed: int = None
"""Switch that allows for fixing q95 only in this constraint equation 68.
(`constraint equation 68`)
"""

pflux_fw_rad_max: float = None
"""Maximum permitted radiation wall load (MW/m^2) (`constraint equation 67`)"""

mvalim: float = None
"""maximum MVA limit (`constraint equation 19`)"""

f_p_beam_shine_through_max: float = None
"""maximum neutral beam shine-through fraction (`constraint equation 59`)"""

nflutfmax: float = None
"""max fast neutron fluence on TF coil (n/m2) (`blktmodel>0`) (`constraint equation 53`)
Also used for demontable magnets (itart = 1) and superconducting coils (i_tf_sup = 1)
and quench protection
To set the CP lifetime (`constraint equation 85`)
"""

p_plasma_separatrix_min_mw: float = None
"""Minimum p_plasma_separatrix_mw [MW] (`constraint equation 80`)"""

f_fw_rad_max: float = None
"""peaking factor for radiation wall load (`constraint equation 67`)"""

pflux_fw_rad_max_mw: float = None
"""Peak radiation wall load (MW/m^2) (`constraint equation 67`)"""

p_plant_electric_net_required_mw: float = None
"""required net electric power (MW) (`constraint equation 16`)"""

p_fusion_total_max_mw: float = None
"""maximum fusion power (MW) (`constraint equation 9`)"""

psepbqarmax: float = None
"""maximum ratio of Psep*Bt/qAR (MWT/m) (`constraint equation 68`)"""

pseprmax: float = None
"""maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
(`constraint equation 56`)
"""

ptfnucmax: float = None
"""maximum nuclear heating in TF coil (MW/m3) (`constraint equation 54`)"""

tbrmin: float = None
"""minimum tritium breeding ratio (`constraint equation 52`)"""

t_burn_min: float = None
"""minimum burn time (s) (KE - no longer itv., see issue #706)"""

t_cycle_min: float = None
"""minimum cycle time (s) (`constraint equation 42`)"""

t_current_ramp_up_min: float = None
"""minimum plasma current ramp-up time (s) (`constraint equation 41`)"""

pflux_fw_neutron_max_mw: float = None
"""allowable neutron wall-load (MW/m2) (`constraint equation 8`)"""

f_alpha_energy_confinement_min: float = None
"""Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement
times (`constraint equation 62`)
"""

zeff_max: float = None
"""maximum value for Zeff (`constraint equation 64`)"""

h_mode_threshold_margin: float = None
"""Sets the constraint bound of the L-H power threshold limit for H-mode

I.e. p_plasma_separatrix_mw / p_l_h_threshold_mw >= h_mode_threshold_margin
"""

l_mode_threshold_margin: float = None
"""Sets the constraint bound of the L-H power threshold limit.

I.e. p_l_h_threshold_mw / p_plasma_separatrix_mw >= l_mode_threshold_margin
"""


def init_constraint_variables():
    """Initialise the constraint variables"""
    global \
        p_hcd_injected_min_mw, \
        beta_poloidal_max, \
        big_q_plasma_min, \
        b_tf_inboard_max, \
        fdene, \
        fradpwr, \
        fiooic, \
        q95_fixed, \
        fjohc, \
        fjohc0, \
        eta_cd_norm_hcd_primary_max, \
        i_q95_fixed, \
        pflux_fw_rad_max, \
        mvalim, \
        f_p_beam_shine_through_max, \
        nflutfmax, \
        p_plasma_separatrix_min_mw, \
        f_fw_rad_max, \
        pflux_fw_rad_max_mw, \
        p_plant_electric_net_required_mw, \
        p_fusion_total_max_mw, \
        psepbqarmax, \
        pseprmax, \
        ptfnucmax, \
        tbrmin, \
        t_burn_min, \
        t_cycle_min, \
        t_current_ramp_up_min, \
        pflux_fw_neutron_max_mw, \
        f_alpha_energy_confinement_min, \
        zeff_max, \
        l_mode_threshold_margin, \
        h_mode_threshold_margin

    p_hcd_injected_min_mw = 0.1
    beta_poloidal_max = 0.19
    big_q_plasma_min = 10.0
    b_tf_inboard_max = 12.0
    fdene = 1.0
    fradpwr = 1.0
    fiooic = 0.7
    q95_fixed = 3.0
    fjohc = 0.7
    fjohc0 = 0.7
    eta_cd_norm_hcd_primary_max = 2.0
    i_q95_fixed = 0
    pflux_fw_rad_max = 1.0
    mvalim = 40.0
    f_p_beam_shine_through_max = 1e-3
    nflutfmax = 1.0e23
    p_plasma_separatrix_min_mw = 150.0
    f_fw_rad_max = 3.33
    pflux_fw_rad_max_mw = 0.0
    p_plant_electric_net_required_mw = 1.0e3
    p_fusion_total_max_mw = 1.5e3
    psepbqarmax = 9.5
    pseprmax = 25.0
    ptfnucmax = 1e-3
    tbrmin = 1.1
    t_burn_min = 1.0
    t_cycle_min = 0.0
    t_current_ramp_up_min = 1.0
    pflux_fw_neutron_max_mw = 1.0
    f_alpha_energy_confinement_min = 5.0
    zeff_max = 3.6
    l_mode_threshold_margin = 1.0
    h_mode_threshold_margin = 1.0
