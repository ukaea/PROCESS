from dataclasses import dataclass


@dataclass
class ConstraintData:
    p_hcd_injected_min_mw: float = 0.1
    """minimum auxiliary power (MW) (`constraint equation 40`)"""

    beta_poloidal_max: float = 0.19
    """maximum poloidal beta (`constraint equation 48`)"""

    big_q_plasma_min: float = 10.0
    """minimum fusion gain Q (`constraint equation 28`)"""

    b_tf_inboard_max: float = 12.0
    """maximum peak toroidal field (T) (`constraint equation 25`)"""

    fdene: float = 1.0
    """Scaling value for density limit (constraint equation 5)"""

    fradpwr: float = 1.0
    """Scaling value for radiation power upper limit (constraint equation 17)"""

    fiooic: float = 0.7
    """Constraint margin for TF coil operating current / critical current ratio
    (`constraint equation 33`)
    """

    q95_fixed: float = 3.0
    """fixed safety factor q at 95% flux surface
    (`constraint equation 68`)
    """

    fjohc: float = 0.7
    """Constraint margin for central solenoid current at end-of-flattop
    (`constraint equation 26`)
    """

    fjohc0: float = 0.7
    """Constraint margin for central solenoid current at beginning of pulse
    (`constraint equation 27`)
    """

    eta_cd_norm_hcd_primary_max: float = 2.0
    """maximum current drive gamma (`constraint equation 37`)"""

    i_q95_fixed: int = 0
    """Switch that allows for fixing q95 only in this constraint equation 68.
    (`constraint equation 68`)
    """

    pflux_fw_rad_max: float = 1.0
    """Maximum permitted radiation wall load (MW/m^2) (`constraint equation 67`)"""

    mvalim: float = 40.0
    """maximum MVA limit (`constraint equation 19`)"""

    f_p_beam_shine_through_max: float = 1e-3
    """maximum neutral beam shine-through fraction (`constraint equation 59`)"""

    nflutfmax: float = 1.0e23
    """max fast neutron fluence on TF coil (n/m2) (`blktmodel>0`) (`constraint equation 53`)
    Also used for demontable magnets (itart = 1) and superconducting coils (i_tf_sup = 1)
    and quench protection
    To set the CP lifetime (`constraint equation 85`)
    """

    p_plasma_separatrix_min_mw: float = 150.0
    """Minimum p_plasma_separatrix_mw [MW] (`constraint equation 80`)"""

    f_fw_rad_max: float = 3.33
    """peaking factor for radiation wall load (`constraint equation 67`)"""

    pflux_fw_rad_max_mw: float = 0.0
    """Peak radiation wall load (MW/m^2) (`constraint equation 67`)"""

    p_plant_electric_net_required_mw: float = 1.0e3
    """required net electric power (MW) (`constraint equation 16`)"""

    p_fusion_total_max_mw: float = 1.5e3
    """maximum fusion power (MW) (`constraint equation 9`)"""

    psepbqarmax: float = 9.5
    """maximum ratio of Psep*Bt/qAR (MWT/m) (`constraint equation 68`)"""

    pseprmax: float = 25.0
    """maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
    (`constraint equation 56`)
    """

    ptfnucmax: float = 1e-3
    """maximum nuclear heating in TF coil (MW/m3) (`constraint equation 54`)"""

    tbrmin: float = 1.1
    """minimum tritium breeding ratio (`constraint equation 52`)"""

    t_burn_min: float = 1.0
    """minimum burn time (s) (KE - no longer itv., see issue #706)"""

    t_cycle_min: float = 0.0
    """minimum cycle time (s) (`constraint equation 42`)"""

    t_current_ramp_up_min: float = 1.0
    """minimum plasma current ramp-up time (s) (`constraint equation 41`)"""

    pflux_fw_neutron_max_mw: float = 1.0
    """allowable neutron wall-load (MW/m2) (`constraint equation 8`)"""

    f_alpha_energy_confinement_min: float = 5.0
    """Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement
    times (`constraint equation 62`)
    """

    zeff_max: float = 3.6
    """maximum value for Zeff (`constraint equation 64`)"""

    f_h_mode_margin: float = 1.0
    """Sets the constraint bound of the L-H power threshold limit for H-mode

    I.e. p_plasma_separatrix_mw / p_l_h_threshold_mw >= f_h_mode_margin
    """

    f_l_mode_margin: float = 1.0
    """Sets the constraint bound of the L-H power threshold limit.

    I.e. p_l_h_threshold_mw / p_plasma_separatrix_mw >= f_l_mode_margin
    """


CREATE_DICTS_FROM_DATACLASS = ConstraintData
