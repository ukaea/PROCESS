from dataclasses import dataclass, field

import numpy as np


@dataclass
class TimesData:
    pulsetimings: int = 1
    """Switch for pulse timings (if i_pulsed_plant=1):

      - =0, t_plant_pulse_plasma_current_ramp_up = Ip(MA)/0.1 t_plant_pulse_coil_precharge, t_plant_pulse_plasma_current_ramp_down = input
      - =1, t_plant_pulse_plasma_current_ramp_up = iteration var or input. t_plant_pulse_coil_precharge/t_plant_pulse_plasma_current_ramp_down max of input or t_plant_pulse_plasma_current_ramp_up"""

    t_plant_pulse_burn: list[float] = field(
        default_factory=lambda: np.array(1000.0, dtype=np.float64)
    )
    """flat-top duration (s) (calculated if `i_pulsed_plant=1`)"""

    t_burn_0: float = 0.0
    """burn time (s) - used for internal consistency"""

    t_plant_pulse_total: list[float] = field(
        default_factory=lambda: np.array(0.0, dtype=np.float64)
    )
    """Total plant pulse cycle time (s)"""

    t_plant_pulse_no_burn: float = 0.0
    """Plant pulse time spent not a flat-top burn phase (s)"""

    t_plant_pulse_dwell: float = 1800.0
    """Plant pulse dwell time before start of next pulse (s) (`iteration variable 17`)"""

    t_plant_pulse_fusion_ramp: float = 10.0
    """time for plasma temperature and density rise to full values (s)"""

    t_pulse_cumulative: list[float] = field(
        default_factory=lambda: np.zeros(6, dtype=np.float64)
    )
    """array of time points during plasma pulse (s)"""

    timelabel: list[str] = field(
        default_factory=lambda: ["Start", "BOP  ", "EOR  ", "BOF  ", "EOF  ", "EOP  "]
    )
    """array of time labels during plasma pulse (s)"""

    intervallabel: list[str] = field(
        default_factory=lambda: [
            "t_plant_pulse_coil_precharge        ",
            "t_plant_pulse_plasma_current_ramp_up  ",
            "t_plant_pulse_fusion_ramp      ",
            "t_plant_pulse_burn             ",
            "t_plant_pulse_plasma_current_ramp_down        ",
        ]
    )

    """time intervals - as strings (s)"""

    t_plant_pulse_plasma_current_ramp_up: float = 30.0
    """Plant pulse time for plasma current to ramp up to approx. full value (s) (calculated if `i_pulsed_plant=0`)
      (`iteration variable 65`)"""

    i_t_current_ramp_up: int = 0
    """Switch for plasma current ramp-up time (if i_pulsed_plant=0):
      - = 0, t_plant_pulse_plasma_current_ramp_up = t_plant_pulse_coil_precharge = t_plant_pulse_plasma_current_ramp_down = Ip(MA)/0.5
      - = 1, t_plant_pulse_plasma_current_ramp_up, t_plant_pulse_coil_precharge, t_plant_pulse_plasma_current_ramp_down are input"""

    t_plant_pulse_plasma_present: float = 0.0
    """Plant pulse time in which a plasma is present (s)"""

    t_plant_pulse_plasma_current_ramp_down: float = 15.0
    """Plant pulse time for plasma current, density, and temperature to ramp down to zero, simultaneously (s); if pulsed, = t_plant_pulse_plasma_current_ramp_up
      the CS and PF coil currents also ramp to zero at the same time"""

    t_plant_pulse_coil_precharge: float = 15.0
    """the time for the central solenoid and PF coils to ramp from zero to max current (s); if pulsed, = t_plant_pulse_plasma_current_ramp_up"""


CREATE_DICTS_FROM_DATACLASS = TimesData
