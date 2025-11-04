from process.data_structure import times_variables
from process.data_structure import (
    build_variables,
    current_drive_variables,
    numerics,
    pfcoil_variables,
    physics_variables,
    stellarator_variables,
)

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
    build_variables.f_dr_tf_outboard_inboard = 1.0

    #  Physics quantities

    physics_variables.beta_norm_max = 0.0
    physics_variables.kappa95 = 1.0
    physics_variables.triang = 0.0
    physics_variables.q95 = 1.03

    #  Turn off current drive

    current_drive_variables.i_hcd_calculations = 0

    #  Times for different phases

    times_variables.t_plant_pulse_coil_precharge = 0.0
    times_variables.t_plant_pulse_plasma_current_ramp_up = 0.0
    times_variables.t_plant_pulse_burn = 3.15576e7  # one year
    times_variables.t_plant_pulse_plasma_current_ramp_down = 0.0
    times_variables.t_plant_pulse_plasma_present = (
        times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_fusion_ramp
        + times_variables.t_plant_pulse_burn
        + times_variables.t_plant_pulse_plasma_current_ramp_down
    )
    times_variables.t_plant_pulse_no_burn = (
        times_variables.t_plant_pulse_coil_precharge
        + times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_plasma_current_ramp_down
        + times_variables.t_plant_pulse_dwell
        + times_variables.t_plant_pulse_fusion_ramp
    )
    times_variables.t_plant_pulse_total = (
        times_variables.t_plant_pulse_coil_precharge
        + times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_fusion_ramp
        + times_variables.t_plant_pulse_burn
        + times_variables.t_plant_pulse_plasma_current_ramp_down
        + times_variables.t_plant_pulse_dwell
    )