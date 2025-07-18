import numpy as np

pulsetimings: float = None
"""Switch for pulse timings (if i_pulsed_plant=1):

   - =0, t_current_ramp_up = Ip(MA)/0.1 t_precharge, t_ramp_down = input
   - =1, t_current_ramp_up = iteration var or input. t_precharge/t_ramp_down max of input or t_current_ramp_up"""

t_burn: float = None
"""flat-top duration (s) (calculated if `i_pulsed_plant=1`)"""

t_burn_0: float = None
"""burn time (s) - used for internal consistency"""

t_cycle: float = None
"""full cycle time (s)"""

tdown: float = None
"""down time (s)"""

t_between_pulse: float = None
"""time between pulses in a pulsed reactor (s) (`iteration variable 17`)"""

t_fusion_ramp: float = None
"""time for plasma temperature and density rise to full values (s)"""

tim: list[float] = None
"""array of time points during plasma pulse (s)"""

timelabel: list[str] = None
"""array of time labels during plasma pulse (s)"""

intervallabel: list[str] = None
"""time intervals - as strings (s)"""

t_current_ramp_up: float = None
"""time for plasma current to ramp up to approx. full value (s) (calculated if `i_pulsed_plant=0`)
   (`iteration variable 65`)"""

i_t_current_ramp_up: int = None
"""Switch for plasma current ramp-up time (if i_pulsed_plant=0):
   - = 0, t_current_ramp_up = t_precharge = t_ramp_down = Ip(MA)/0.5
   - = 1, t_current_ramp_up, t_precharge, t_ramp_down are input"""

t_pulse_repetition: float = None
"""pulse length = t_current_ramp_up + t_fusion_ramp + t_burn + t_ramp_down"""

t_ramp_down: float = None
"""time for plasma current, density, and temperature to ramp down to zero, simultaneously (s); if pulsed, = t_current_ramp_up
   the CS and PF coil currents also ramp to zero at the same time"""

t_precharge: float = None
"""the time for the central solenoid and PF coils to ramp from zero to max current (s); if pulsed, = t_current_ramp_up"""


def init_times_variables():
    """Initialise plasma pulse timing variables"""
    global pulsetimings
    global t_burn
    global t_burn_0
    global t_cycle
    global tdown
    global t_between_pulse
    global t_fusion_ramp
    global tim
    global timelabel
    global intervallabel
    global t_current_ramp_up
    global i_t_current_ramp_up
    global t_pulse_repetition
    global t_ramp_down
    global t_precharge

    pulsetimings = 1.0
    t_burn = 1000.0
    t_burn_0 = 0.0
    t_cycle = 0.0
    tdown = 0.0
    t_between_pulse = 1800.0
    t_fusion_ramp = 10.0
    tim = np.zeros(6, dtype=np.float64)
    timelabel = ["Start", "BOP  ", "EOR  ", "BOF  ", "EOF  ", "EOP  "]
    intervallabel = [
        "t_precharge        ",
        "t_current_ramp_up  ",
        "t_fusion_ramp      ",
        "t_burn             ",
        "t_ramp_down        ",
    ]
    t_current_ramp_up = 30.0
    i_t_current_ramp_up = 0
    t_pulse_repetition = 0.0
    t_ramp_down = 15.0
    t_precharge = 15.0
