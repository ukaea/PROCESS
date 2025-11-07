p_hcd_injected_min_mw: float = None
"""minimum auxiliary power (MW) (`constraint equation 40`)"""


beta_poloidal_max: float = None
"""maximum poloidal beta (`constraint equation 48`)"""


big_q_plasma_min: float = None
"""minimum fusion gain Q (`constraint equation 28`)"""


b_tf_inboard_max: float = None
"""maximum peak toroidal field (T) (`constraint equation 25`)"""


fp_hcd_injected_min_mw: float = None
"""f-value for minimum auxiliary power (`constraint equation 40`, `iteration variable 64`)"""


fbeta_poloidal_eps: float = None
"""f-value for epsilon beta-poloidal (`constraint equation 6`, `iteration variable 8`)"""


fbeta_poloidal: float = None
"""f-value for poloidal beta (`constraint equation 48`, `iteration variable 79`)"""


fbeta_max: float = None
"""f-value for beta limit (`constraint equation 24`, `iteration variable 36`)"""


fbeta_min: float = None
"""f-value for (lower) beta limit (`constraint equation 84`, `iteration variable 173`)"""


fc_tf_turn_max: float = None
"""f-value for TF coil current per turn upper limit
(`constraint equation 77`, `iteration variable 146`)
"""


fr_conducting_wall: float = None
"""f-value for conducting wall radius / rminor limit
(`constraint equation 23`, `iteration variable 104`)
"""


fdene: float = None
"""f-value for density limit (`constraint equation 5`, `iteration variable 9`)
(invalid if `i_plasma_pedestal=3`)
"""


fdtmp: float = None
"""f-value for first wall coolant temperature rise
(`constraint equation 38`, `iteration variable 62`)
"""


fecrh_ignition: float = None
"""f-value for ecrh ignition constraint
(`constraint equation 91`, `iteration variable 168`)
"""


fflutf: float = None
"""f-value for neutron fluence on TF coil (`constraint equation 53`, `iteration variable 92`)"""


fp_fusion_total_max_mw: float = None
"""f-value for maximum fusion power (`constraint equation 9`, `iteration variable 26`)"""


feta_cd_norm_hcd_primary_max: float = None
"""f-value for current drive gamma (`constraint equation 37`, `iteration variable 40`)"""


fpflux_div_heat_load_mw: float = None
"""f-value for divertor heat load (`constraint equation 18`, `iteration variable 27`)"""


fiooic: float = None
"""f-value for TF coil operating current / critical current ratio
(`constraint equation 33`, `iteration variable 50`)
"""


fipir: float = None
"""f-value for Ip/Irod upper limit
constraint equation icc = 46
iteration variable ixc = 72
"""


q95_fixed: float = None
"""fixed safety factor q at 95% flux surface
(`constraint equation 68`)
"""


fjohc: float = None
"""f-value for central solenoid current at end-of-flattop
(`constraint equation 26`, `iteration variable 38`)
"""


fjohc0: float = None
"""f-value for central solenoid current at beginning of pulse
(`constraint equation 27`, `iteration variable 39`)
"""


fjprot: float = None
"""f-value for TF coil winding pack current density
(`constraint equation 35`, `iteration variable 53`)
"""


fl_h_threshold: float = None
"""f-value for L-H power threshold (`constraint equation 15`, `iteration variable 103`)"""


fmva: float = None
"""f-value for maximum MVA (`constraint equation 19`, `iteration variable 30`)"""


fnbshinef: float = None
"""f-value for maximum neutral beam shine-through fraction
(`constraint equation 59`, `iteration variable 105`)
"""


fncycle: float = None
"""f-value for minimum CS coil stress load cycles
(`constraint equation 90`, `iteration variable 167`)
"""


fnesep: float = None
"""f-value for Eich critical separatrix density
(`constraint equation 76`, `iteration variable 144`)
"""


foh_stress: float = None
"""f-value for Tresca yield criterion in Central Solenoid
(`constraint equation 72`, `iteration variable 123`)
"""


fb_tf_inboard_max: float = None
"""f-value for maximum toroidal field (`constraint equation 25`, `iteration variable 35`)"""


fp_hcd_injected_max: float = None
"""f-value for injection power (`constraint equation 30`, `iteration variable 46`)"""


fp_plant_electric_net_required_mw: float = None
"""f-value for net electric power (`constraint equation 16`, `iteration variable 25`)"""


fradius_beam_tangency: float = None
"""f-value for neutral beam tangency radius limit
(`constraint equation 20`, `iteration variable 33`)
"""


fpsepbqar: float = None
"""f-value for maximum Psep*Bt/qAR limit (`constraint equation 68`, `iteration variable 117`)"""


fpsepr: float = None
"""f-value for maximum Psep/R limit (`constraint equation 56`, `iteration variable 97`)"""


fptemp: float = None
"""f-value for peak centrepost temperature (`constraint equation 44`, `iteration variable 68`)"""


fptfnuc: float = None
"""f-value for maximum TF coil nuclear heating (`constraint equation 54`, `iteration variable 95`)"""


fq95_min: float = None
"""f-value for edge safety factor (`constraint equation 45`, `iteration variable 71`)"""


fbig_q_plasma_min: float = None
"""f-value for Q (`constraint equation 28`, `iteration variable 45`)"""


fradpwr: float = None
"""f-value for core radiation power limit (`constraint equation 17`, `iteration variable 28`)"""


fpflux_fw_rad_max: float = None
"""f-value for upper limit on radiation wall load (`constr. equ. 67`, `iteration variable 116`)"""


freinke: float = None
"""f-value for Reinke detachment criterion (`constr. equ. 78`, `iteration variable 147`)"""


frminor: float = None
"""f-value for minor radius limit (`constraint equation 21`, `iteration variable 32`)"""


fstrcase: float = None
"""f-value for maximum TF coil case Tresca yield criterion
(`constraint equation 31`, `iteration variable 48`)
"""


fstrcond: float = None
"""f-value for maxiumum TF coil conduit Tresca yield criterion
(`constraint equation 32`, `iteration variable 49`)
"""


fstr_wp: float = None
"""f-value for maxiumum TF coil strain absolute value
(`constraint equation 88`, `iteration variable 165`)
"""


fmaxvvstress: float = None
"""f-value for maximum permitted stress of the VV
(`constraint equation 65`, `iteration variable 113`)
"""


ftbr: float = None
"""f-value for minimum tritium breeding ratio (`constraint equation 52`, `iteration variable 89`)"""


ft_burn_min: float = None
"""f-value for minimum burn time (`constraint equation 13`, `iteration variable 21`)"""


ft_cycle_min: float = None
"""f-value for cycle time (`constraint equation 42`, `iteration variable 67`)"""


ftmargoh: float = None
"""f-value for central solenoid temperature margin
(`constraint equation 60`, `iteration variable 106`)
"""


ftmargtf: float = None
"""f-value for TF coil temperature margin (`constraint equation 36`, `iteration variable 54`)"""


ft_current_ramp_up: float = None
"""f-value for plasma current ramp-up time (`constraint equation 41`, `iteration variable 66`)"""


ftemp_fw_max: float = None
"""f-value for first wall peak temperature (`constraint equation 39`, `iteration variable 63`)"""


fvdump: float = None
"""f-value for dump voltage (`constraint equation 34`, `iteration variable 51`)"""


fvs_plasma_total_required: float = None
"""f-value for flux-swing (V-s) requirement (STEADY STATE)
(`constraint equation 12`, `iteration variable 15`)
"""


fvvhe: float = None
"""f-value for vacuum vessel He concentration limit (`i_blanket_type = 2`)
(`constraint equation 55`, `iteration variable 96`)
"""


fpflux_fw_neutron_max_mw: float = None
"""f-value for maximum wall load (`constraint equation 8`, `iteration variable 14`)"""


fzeff_max: float = None
"""f-value for maximum n_charge_plasma_effective_vol_avg (`constraint equation 64`, `iteration variable 112`)"""


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


falpha_energy_confinement: float = None
"""f-value for lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy
confinement times (`constraint equation 62`, `iteration variable 110`)
"""


fniterpump: float = None
"""f-value for constraint that number of pumps < tfno
(`constraint equation 63`, `iteration variable 111`)
"""


zeff_max: float = None
"""maximum value for Zeff (`constraint equation 64`)"""


fpoloidalpower: float = None
"""f-value for constraint on rate of change of energy in poloidal field
(`constraint equation 66`, `iteration variable 115`)
"""


ftemp_croco_quench_max: float = None
"""TF coil quench temparature remains below temp_croco_quench_max
(`constraint equation 74`, `iteration variable 141`)
"""


def init_constraint_variables():
    """Initialise the constraint variables"""
    global p_hcd_injected_min_mw
    global beta_poloidal_max
    global big_q_plasma_min
    global b_tf_inboard_max
    global fp_hcd_injected_min_mw
    global fbeta_poloidal_eps
    global fbeta_poloidal
    global fbeta_max
    global fbeta_min
    global fc_tf_turn_max
    global fr_conducting_wall
    global fdene
    global fdtmp
    global fecrh_ignition
    global fflutf
    global fp_fusion_total_max_mw
    global feta_cd_norm_hcd_primary_max
    global fpflux_div_heat_load_mw
    global fiooic
    global fipir
    global q95_fixed
    global fjohc
    global fjohc0
    global fjprot
    global fl_h_threshold
    global fmva
    global fnbshinef
    global fncycle
    global fnesep
    global foh_stress
    global fb_tf_inboard_max
    global fp_hcd_injected_max
    global fp_plant_electric_net_required_mw
    global fradius_beam_tangency
    global fpsepbqar
    global fpsepr
    global fptemp
    global fptfnuc
    global fq95_min
    global fbig_q_plasma_min
    global fradpwr
    global fpflux_fw_rad_max
    global freinke
    global frminor
    global fstrcase
    global fstrcond
    global fstr_wp
    global fmaxvvstress
    global ftbr
    global ft_burn_min
    global ft_cycle_min
    global ftmargoh
    global ftmargtf
    global ft_current_ramp_up
    global ftemp_fw_max
    global fvdump
    global fvs_plasma_total_required
    global fvvhe
    global fpflux_fw_neutron_max_mw
    global fzeff_max
    global eta_cd_norm_hcd_primary_max
    global i_q95_fixed
    global pflux_fw_rad_max
    global mvalim
    global f_p_beam_shine_through_max
    global nflutfmax
    global p_plasma_separatrix_min_mw
    global f_fw_rad_max
    global pflux_fw_rad_max_mw
    global p_plant_electric_net_required_mw
    global p_fusion_total_max_mw
    global psepbqarmax
    global pseprmax
    global ptfnucmax
    global tbrmin
    global t_burn_min
    global t_cycle_min
    global t_current_ramp_up_min
    global pflux_fw_neutron_max_mw
    global f_alpha_energy_confinement_min
    global falpha_energy_confinement
    global fniterpump
    global zeff_max
    global fpoloidalpower
    global ftemp_croco_quench_max

    p_hcd_injected_min_mw = 0.1
    beta_poloidal_max = 0.19
    big_q_plasma_min = 10.0
    b_tf_inboard_max = 12.0
    fp_hcd_injected_min_mw = 1.0
    fbeta_poloidal_eps = 1.0
    fbeta_poloidal = 1.0
    fbeta_max = 1.0
    fbeta_min = 1.0
    fc_tf_turn_max = 1.0
    fr_conducting_wall = 1.0
    fdene = 1.0
    fdtmp = 1.0
    fflutf = 1.0
    fp_fusion_total_max_mw = 1.0
    feta_cd_norm_hcd_primary_max = 1.0
    fpflux_div_heat_load_mw = 1.0
    fiooic = 0.5
    fipir = 1.0
    q95_fixed = 3.0
    fjohc = 1.0
    fjohc0 = 1.0
    fjprot = 1.0
    fl_h_threshold = 1.0
    fmva = 1.0
    fnbshinef = 1.0
    fncycle = 1.0
    fnesep = 1.0
    foh_stress = 1.0
    fb_tf_inboard_max = 1.0
    fp_hcd_injected_max = 1.0
    fp_plant_electric_net_required_mw = 1.0
    fradius_beam_tangency = 1.0
    fpsepbqar = 1.0
    fpsepr = 1.0
    fptemp = 1.0
    fptfnuc = 1.0
    fq95_min = 1.0
    fbig_q_plasma_min = 1.0
    fradpwr = 0.99
    fpflux_fw_rad_max = 1.0
    freinke = 1.0
    frminor = 1.0
    fstrcase = 1.0
    fstrcond = 1.0
    fstr_wp = 1.0
    fmaxvvstress = 1.0
    ftbr = 1.0
    ft_burn_min = 1.0
    ft_cycle_min = 1.0
    ftmargoh = 1.0
    ftmargtf = 1.0
    ft_current_ramp_up = 1.0
    ftemp_fw_max = 1.0
    fvdump = 1.0
    fvs_plasma_total_required = 1.0
    fvvhe = 1.0
    fpflux_fw_neutron_max_mw = 1.0
    fzeff_max = 1.0
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
    falpha_energy_confinement = 1.0
    fniterpump = 1.0
    zeff_max = 3.6
    fpoloidalpower = 1.0
    ftemp_croco_quench_max = 1.0
    fecrh_ignition = 1.0
