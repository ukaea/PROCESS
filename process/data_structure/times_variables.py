import numpy as np

pulsetimings: float = None
"""Switch for pulse timings (if i_pulsed_plant=1):

   - =0, t_plant_pulse_plasma_current_ramp_up = Ip(MA)/0.1 t_plant_pulse_coil_precharge, t_plant_pulse_plasma_current_ramp_down = input
   - =1, t_plant_pulse_plasma_current_ramp_up = iteration var or input. t_plant_pulse_coil_precharge/t_plant_pulse_plasma_current_ramp_down max of input or t_plant_pulse_plasma_current_ramp_up"""

t_plant_pulse_burn: float = None
"""flat-top duration (s) (calculated if `i_pulsed_plant=1`)"""

t_burn_0: float = None
"""burn time (s) - used for internal consistency"""

t_plant_pulse_total: float = None
"""Total plant pulse cycle time (s)"""

t_plant_pulse_no_burn: float = None
"""Plant pulse time spent not a flat-top burn phase (s)"""

t_plant_pulse_dwell: float = None
"""Plant pulse dwell time before start of next pulse (s) (`iteration variable 17`)"""

t_plant_pulse_fusion_ramp: float = None
"""time for plasma temperature and density rise to full values (s)"""

t_pulse_cumulative: list[float] = None
"""array of time points during plasma pulse (s)"""

timelabel: list[str] = None
"""array of time labels during plasma pulse (s)"""

intervallabel: list[str] = None
"""time intervals - as strings (s)"""

t_plant_pulse_plasma_current_ramp_up: float = None
"""Plant pulse time for plasma current to ramp up to approx. full value (s) (calculated if `i_pulsed_plant=0`)
   (`iteration variable 65`)"""

i_t_current_ramp_up: int = None
"""Switch for plasma current ramp-up time (if i_pulsed_plant=0):
   - = 0, t_plant_pulse_plasma_current_ramp_up = t_plant_pulse_coil_precharge = t_plant_pulse_plasma_current_ramp_down = Ip(MA)/0.5
   - = 1, t_plant_pulse_plasma_current_ramp_up, t_plant_pulse_coil_precharge, t_plant_pulse_plasma_current_ramp_down are input"""

t_plant_pulse_plasma_present: float = None
"""Plant pulse time in which a plasma is present (s)"""

t_plant_pulse_plasma_current_ramp_down: float = None
"""Plant pulse time for plasma current, density, and temperature to ramp down to zero, simultaneously (s); if pulsed, = t_plant_pulse_plasma_current_ramp_up
   the CS and PF coil currents also ramp to zero at the same time"""

t_plant_pulse_coil_precharge: float = None
"""the time for the central solenoid and PF coils to ramp from zero to max current (s); if pulsed, = t_plant_pulse_plasma_current_ramp_up"""


def init_times_variables():
    """Initialise plasma pulse timing variables"""
    global \
        pulsetimings, \
        t_plant_pulse_burn, \
        t_burn_0, \
        t_plant_pulse_total, \
        t_plant_pulse_no_burn, \
        t_plant_pulse_dwell, \
        t_plant_pulse_fusion_ramp, \
        t_pulse_cumulative, \
        timelabel, \
        intervallabel, \
        t_plant_pulse_plasma_current_ramp_up, \
        i_t_current_ramp_up, \
        t_plant_pulse_plasma_present, \
        t_plant_pulse_plasma_current_ramp_down, \
        t_plant_pulse_coil_precharge

    pulsetimings = 1.0
    t_plant_pulse_burn = np.array(1000.0, dtype=np.float64)
    t_burn_0 = 0.0
    t_plant_pulse_total = np.array(0.0, dtype=np.float64)
    t_plant_pulse_no_burn = 0.0
    t_plant_pulse_dwell = 1800.0
    t_plant_pulse_fusion_ramp = 10.0
    t_pulse_cumulative = np.zeros(6, dtype=np.float64)
    timelabel = ["Start", "BOP  ", "EOR  ", "BOF  ", "EOF  ", "EOP  "]
    intervallabel = [
        "t_plant_pulse_coil_precharge        ",
        "t_plant_pulse_plasma_current_ramp_up  ",
        "t_plant_pulse_fusion_ramp      ",
        "t_plant_pulse_burn             ",
        "t_plant_pulse_plasma_current_ramp_down        ",
    ]
    t_plant_pulse_plasma_current_ramp_up = 30.0
    i_t_current_ramp_up = 0
    t_plant_pulse_plasma_present = 0.0
    t_plant_pulse_plasma_current_ramp_down = 15.0
    t_plant_pulse_coil_precharge = 15.0
