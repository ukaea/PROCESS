"""Module containing variables for the resistive TF coil models"""

from dataclasses import dataclass


@dataclass(slots=True)
class ResistiveTFData:
    """Dataclass holding resistive TF coil variables"""

    a_res_tf_coil_conductor: float = 0.0
    """Area of resistive conductor in resistive TF coil [m2]"""

    cdtfleg: float = 0.0  # questioning this
    """TF outboard leg current density (A/m2)"""

    res_tf_leg: float = 0.0
    """TF coil leg resistance (ohm)"""

    p_cp_resistive_mw: float = 0.0
    """Peak resistive TF coil inboard leg power (MW)"""

    p_tf_joints_resistive_mw: float = 0.0
    """TF joints resistive power losses (MW)"""

    p_tf_leg_resistive_mw: float = 0.0
    """TF coil outboard leg resistive power (MW)"""

    rho_cp: float = 0.0
    """TF coil inboard leg resistivity [Ohm-m]. If `itart=0`, this variable is the
    average resistivity over the whole magnet
    """

    rho_tf_leg: float = 0.0
    """Resistivity of a TF coil leg (Ohm-m)"""

    rho_tf_bus: float = 1.86e-8
    """Resistivity of a TF coil bus (Ohm-m). Default values is for that of GLIDCOP AL-15 (C15715) at 293K"""

    frhocp: float = 1.0
    """Centrepost resistivity enhancement factor. For `itart=0`, this factor
    is used for the whole magnet
    """

    frholeg: float = 1.0
    """Outboard legs resistivity enhancement factor. Only used for `itart=1`."""

    p_tf_joints_resistive: float = 0.0
    """Calculated TF joints resistive power losses [W]"""

    vtfkv: float = 0.0
    """TF coil voltage for resistive coil including bus (kV)"""
