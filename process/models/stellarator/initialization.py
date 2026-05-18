from process.core.model import DataStructure
from process.data_structure import (
    numerics,
    physics_variables,
    stellarator_variables,
)


def st_init(data: DataStructure):
    """Routine to initialise the variables relevant to stellarators

    This routine initialises the variables relevant to stellarators.
    Many of these may override the values set in routine

    """
    if stellarator_variables.istell == 0:
        return

    numerics.boundu[0] = 40.0  # allow higher aspect ratio

    # These lines switch off tokamak specifics (solenoid, pf coils, pulses etc.).
    # Are they still up to date? (26/07/22 JL)

    # Build quantities

    data.build.dr_cs = 0.0
    data.build.iohcl = 0
    data.pf_coil.f_z_cs_tf_internal = 0.0
    data.build.dr_cs_tf_gap = 0.0
    data.build.f_dr_tf_outboard_inboard = 1.0

    #  Physics quantities

    physics_variables.i_plasma_pedestal = 0
    physics_variables.beta_norm_max = 0.0
    physics_variables.kappa95 = 1.0
    physics_variables.triang = 0.0
    physics_variables.q95 = 1.03

    #  Turn off current drive

    data.current_drive.i_hcd_calculations = 0

    #  Times for different phases

    data.times.t_plant_pulse_coil_precharge = 0.0
    data.times.t_plant_pulse_plasma_current_ramp_up = 0.0
    data.times.t_plant_pulse_burn = 3.15576e7  # one year
    data.times.t_plant_pulse_plasma_current_ramp_down = 0.0
    data.times.t_plant_pulse_plasma_present = (
        data.times.t_plant_pulse_plasma_current_ramp_up
        + data.times.t_plant_pulse_fusion_ramp
        + data.times.t_plant_pulse_burn
        + data.times.t_plant_pulse_plasma_current_ramp_down
    )
    data.times.t_plant_pulse_no_burn = (
        data.times.t_plant_pulse_coil_precharge
        + data.times.t_plant_pulse_plasma_current_ramp_up
        + data.times.t_plant_pulse_plasma_current_ramp_down
        + data.times.t_plant_pulse_dwell
        + data.times.t_plant_pulse_fusion_ramp
    )
    data.times.t_plant_pulse_total = (
        data.times.t_plant_pulse_coil_precharge
        + data.times.t_plant_pulse_plasma_current_ramp_up
        + data.times.t_plant_pulse_fusion_ramp
        + data.times.t_plant_pulse_burn
        + data.times.t_plant_pulse_plasma_current_ramp_down
        + data.times.t_plant_pulse_dwell
    )
