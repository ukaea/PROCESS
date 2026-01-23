gamma_he: float = None
"""ratio of specific heats for helium (`i_p_coolant_pumping=3`)"""

t_in_bb: float = None
"""temperature in FW and blanket coolant at blanket entrance (`i_p_coolant_pumping=3`) [K]"""

t_out_bb: float = None
"""temperature in FW and blanket coolant at blanket exit (`i_p_coolant_pumping=3`) [K]"""

p_he: float = None
"""pressure in FW and blanket coolant at pump exit (`i_p_coolant_pumping=3`) [Pa]"""

dp_he: float = None
"""pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

dp_fw_blkt: float = None
"""pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

dp_fw: float = None
"""pressure drop in FW coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

dp_blkt: float = None
"""pressure drop in blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

dp_liq: float = None
"""pressure drop in liquid metal blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

p_fw_blkt_coolant_pump_mw: float = None
"""mechanical pumping power for FW and blanket including heat exchanger and
pipes (`i_p_coolant_pumping=3`) [MW]
"""

f_p_fw_blkt_pump: float = None
"""Pumping power for FW and Blanket multiplier factor"""


def init_primary_pumping_variables():
    """Initialise primary pumping variables"""
    global \
        gamma_he, \
        t_in_bb, \
        t_out_bb, \
        p_he, \
        dp_he, \
        dp_fw_blkt, \
        dp_fw, \
        dp_blkt, \
        dp_liq, \
        p_fw_blkt_coolant_pump_mw, \
        f_p_fw_blkt_pump

    gamma_he = 1.667  # Ratio of specific heats  Helium
    t_in_bb = 573.13  # K
    t_out_bb = 773.13  # K
    p_he = 8.0e6  # Pa
    dp_he = 5.5e5  # Pa
    dp_fw_blkt = 1.5e5  # Pa
    dp_fw = 1.5e5  # Pa
    dp_blkt = 3.5e3  # Pa
    dp_liq = 1.0e7  # Pa
    p_fw_blkt_coolant_pump_mw = 0.0
    f_p_fw_blkt_pump = 1.0
