"""Module to calculate quench protection limits for stellarator coils"""

import numpy as np

from process.data_structure import (
    build_variables,
    physics_variables,
    rebco_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
)


def calculate_quench_protection(coilcurrent):
    """
    Calculate quench protecion limits for stellarator coils
    Includes calculation of the vacuum vessel force density, quench protection current density
    and max dump voltage during quench
    coilcurrent : Total coils current (MA)
    """
    #
    # This copied from the tokamak module:
    # Radial position of vacuum vessel [m]
    rad_vv_in = (
        physics_variables.rmajor
        - physics_variables.rminor
        - build_variables.dr_fw_plasma_gap_inboard
        - build_variables.dr_fw_inboard
        - build_variables.dr_blkt_inboard
        - build_variables.dr_shld_blkt_gap
        - build_variables.dr_shld_inboard
    )
    rad_vv_out = (
        physics_variables.rmajor
        + physics_variables.rminor
        + build_variables.dr_fw_plasma_gap_outboard
        + build_variables.dr_fw_outboard
        + build_variables.dr_blkt_outboard
        + build_variables.dr_shld_blkt_gap
        + build_variables.dr_shld_outboard
    )

    # Stellarator version is working on the W7-X scaling, so we should actual use vv r_major
    # plasma r_major is just an approximation, but exact calculations require 3D geometry
    # Maybe it can be added to the stella_config file in the future
    rad_vv = physics_variables.rmajor

    # MN/m^3
    f_vv_actual = calculate_vv_max_force_density_from_W7X_scaling(rad_vv)

    # This approach merge stress model from tokamaks with induced force calculated from W7-X scaling
    a_vv = (rad_vv_out + rad_vv_in) / (rad_vv_out - rad_vv_in)
    zeta = 1 + ((a_vv - 1) * np.log((a_vv + 1) / (a_vv - 1)) / (2 * a_vv))

    superconducting_tf_coil_variables.vv_stress_quench = (
        zeta * f_vv_actual * 1e6 * rad_vv_in
    )

    # comparison
    # the new quench protection routine, see #1047
    tfcoil_variables.j_tf_wp_quench_heat_max = (
        calculate_quench_protection_current_density(
            tau_quench=tfcoil_variables.t_tf_superconductor_quench,
            t_detect=tfcoil_variables.t_tf_quench_detection,
            f_cu=tfcoil_variables.f_a_tf_turn_cable_copper,
            f_cond=1 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
            temp=tfcoil_variables.tftmp,
            a_cable=tfcoil_variables.a_tf_turn_cable_space_no_void,
            a_turn=tfcoil_variables.dx_tf_turn_general**2,
        )
    )

    # Also give the copper current density (copper A/m2) for REBCO quench calculations:
    rebco_variables.coppera_m2 = (
        coilcurrent
        * 1.0e6
        / (
            tfcoil_variables.a_tf_wp_conductor
            * tfcoil_variables.f_a_tf_turn_cable_copper
        )
    )

    # Max volatage during fast discharge of TF coil (V)
    # (note that tf_coil_variable is in kV, while calculation is in V)
    tfcoil_variables.v_tf_coil_dump_quench_kv = (
        max_dump_voltage(
            tf_energy_stored=(
                tfcoil_variables.e_tf_magnetic_stored_total_gj
                / tfcoil_variables.n_tf_coils
                * 1.0e9
            ),
            t_dump=tfcoil_variables.t_tf_superconductor_quench,
            current=tfcoil_variables.c_tf_turn,
        )
        / 1.0e3
    )  # turn into kV

    return f_vv_actual


def calculate_vv_max_force_density_from_W7X_scaling(rad_vv: float) -> float:
    """Actual VV force density from scaling [MN/m^3]
    Based on reference values from W-7X.
    """
    force_density_ref = 2.54  # MN/m^3
    b_ref = 3.0  # T
    i_total_ref = 1.3e6 * 50
    rminor_ref = 0.92  # m
    tau_ref = 3.0  # s
    rmajor_ref = 5.2  # m
    dr_vv_ref = 14e-3  # m, thickness of VV

    return (
        force_density_ref
        * (
            b_ref
            / physics_variables.b_plasma_toroidal_on_axis
            * i_total_ref
            / tfcoil_variables.c_tf_total
            * rminor_ref**2
            / physics_variables.rminor**2
        )
        ** (-1)
        * (
            tau_ref
            / tfcoil_variables.t_tf_superconductor_quench
            * rmajor_ref
            / rad_vv
            * dr_vv_ref
            / ((build_variables.dr_vv_inboard + build_variables.dr_vv_outboard) / 2)
        )
    )


def max_dump_voltage(tf_energy_stored: float, t_dump: float, current: float) -> float:
    """
    Max volatage during fast discharge of TF coil (V)
    tf_energy_stored : Energy stored in one TF coil (J)
    t_dump : Dump time (sec)
    current : Operating current (A)
    """
    return 2 * (tf_energy_stored / t_dump) / current


def calculate_quench_protection_current_density(
    tau_quench, t_detect, f_cu, f_cond, temp, a_cable, a_turn
):
    """
    Calculates the current density limited by the protection limit.

    Simplified 0-D adiabatic heat balance "hotspot criterion" model.

    This is slightly diffrent that tokamak version (also diffrent from the stellarator paper).
    We skip the superconduc6tor contribution (this should be more conservative in theory).
    tau_quench : Quench time (s)
    t_detect : Detection time (s)
    f_cu : Copper fraction
    f_cond : Conductor fraction
    temp : peak helium coolant temperature in TF coils and PF coils (K)
    a_cable : Cable space area (per turn)  [m2] (Includes the area of voids and central helium channel)
    a_turn : TF coil turn cross section area [m2]
    """
    temp_k = [4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 104, 114, 124]
    q_cu_array_sa2m4 = [
        1.08514e17,
        1.12043e17,
        1.12406e17,
        1.05940e17,
        9.49741e16,
        8.43757e16,
        7.56346e16,
        6.85924e16,
        6.28575e16,
        5.81004e16,
        5.40838e16,
        5.06414e16,
        4.76531e16,
    ]
    q_he_array_sa2m4 = [
        3.44562e16,
        9.92398e15,
        4.90462e15,
        2.41524e15,
        1.26368e15,
        7.51617e14,
        5.01632e14,
        3.63641e14,
        2.79164e14,
        2.23193e14,
        1.83832e14,
        1.54863e14,
        1.32773e14,
    ]

    q_he = np.interp(temp, temp_k, q_he_array_sa2m4)
    q_cu = np.interp(temp, temp_k, q_cu_array_sa2m4)

    # This leaves out the contribution from the superconductor fraction for now
    return (a_cable / a_turn) * np.sqrt(
        1
        / (0.5 * tau_quench + t_detect)
        * (f_cu**2 * f_cond**2 * q_cu + f_cu * f_cond * (1 - f_cond) * q_he)
    )
