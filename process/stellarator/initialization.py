from process.data_structure import times_variables
from process.fortran import (
    build_variables,
    current_drive_variables,
    numerics,
    pfcoil_variables,
    physics_variables,
    stellarator_variables,
    stellarator_module as st,
)

def init_stellarator_variables():
    stellarator_variables.istell = 0
    stellarator_variables.bmn = 1e-3
    stellarator_variables.f_asym = 1.0
    stellarator_variables.f_rad = 0.85
    stellarator_variables.f_w = 0.5
    stellarator_variables.fdivwet = 0.333333333333333
    stellarator_variables.flpitch = 1e-3
    stellarator_variables.hportamax = 0.0
    stellarator_variables.hportpmax = 0.0
    stellarator_variables.hporttmax = 0.0
    stellarator_variables.iotabar = 1.0
    stellarator_variables.isthtr = 3
    stellarator_variables.m_res = 5
    stellarator_variables.n_res = 5
    stellarator_variables.shear = 0.5
    stellarator_variables.vportamax = 0.0
    stellarator_variables.vportpmax = 0.0
    stellarator_variables.vporttmax = 0.0
    stellarator_variables.max_gyrotron_frequency = 1.0e9
    stellarator_variables.te0_ecrh_achievable = 1.0e2
    stellarator_variables.f_st_coil_aspect = 1.0


def init_stellarator_module():
    st.first_call = True
    st.first_call_stfwbs = True
    st.f_n = 0.0
    st.f_r = 0.0
    st.f_a = 0.0
    st.f_b = 0.0
    st.f_i = 0.0
    st.f_coil_aspect = 0.0
    st.r_coil_major = 0.0
    st.r_coil_minor = 0.0
    st.f_coil_shape = 0.0


def st_init():
    """Routine to initialise the variables relevant to stellarators
    author: P J Knight, CCFE, Culham Science Centre
    author: F Warmer, IPP Greifswald

    This routine initialises the variables relevant to stellarators.
    Many of these may override the values set in routine
    """
    if stellarator_variables.istell == 0:
        return

    numerics.boundu[0] = 40.0  # allow higher aspect ratio

    # These lines switch off tokamak specifics (solenoid, pf coils, pulses etc.).
    # Are they still up to date? (26/07/22 JL)

    # Build quantities

    build_variables.dr_cs = 0.0
    build_variables.iohcl = 0
    pfcoil_variables.f_z_cs_tf_internal = 0.0
    build_variables.dr_cs_tf_gap = 0.0
    build_variables.tfootfi = 1.0

    #  Physics quantities

    physics_variables.beta_norm_max = 0.0
    physics_variables.kappa95 = 1.0
    physics_variables.triang = 0.0
    physics_variables.q95 = 1.03

    #  Turn off current drive

    current_drive_variables.i_hcd_calculations = 0

    #  Times for different phases

    times_variables.t_precharge = 0.0
    times_variables.t_current_ramp_up = 0.0
    times_variables.t_burn = 3.15576e7  # one year
    times_variables.t_ramp_down = 0.0
    times_variables.t_pulse_repetition = (
        times_variables.t_current_ramp_up
        + times_variables.t_fusion_ramp
        + times_variables.t_burn
        + times_variables.t_ramp_down
    )
    times_variables.tdown = (
        times_variables.t_precharge
        + times_variables.t_current_ramp_up
        + times_variables.t_ramp_down
        + times_variables.t_between_pulse
    )
    times_variables.t_cycle = (
        times_variables.t_precharge
        + times_variables.t_current_ramp_up
        + times_variables.t_fusion_ramp
        + times_variables.t_burn
        + times_variables.t_ramp_down
        + times_variables.t_between_pulse
    )