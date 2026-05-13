from dataclasses import dataclass


@dataclass
class PrimaryPumpingData:
    gamma_he: float = 1.667
    """ratio of specific heats for helium (`i_p_coolant_pumping=3`)"""

    t_in_bb: float = 573.13
    """temperature in FW and blanket coolant at blanket entrance (`i_p_coolant_pumping=3`) [K]"""

    t_out_bb: float = 773.13
    """temperature in FW and blanket coolant at blanket exit (`i_p_coolant_pumping=3`) [K]"""

    p_he: float = 8.0e6
    """pressure in FW and blanket coolant at pump exit (`i_p_coolant_pumping=3`) [Pa]"""

    dp_he: float = 5.5e5
    """pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_fw_blkt: float = 1.5e5
    """pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_fw: float = 1.5e5
    """pressure drop in FW coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_blkt: float = 3.5e3
    """pressure drop in blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_liq: float = 1.0e7
    """pressure drop in liquid metal blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    p_fw_blkt_coolant_pump_mw: float = 0.0
    """mechanical pumping power for FW and blanket including heat exchanger and
    pipes (`i_p_coolant_pumping=3`) [MW]
    """

    f_p_fw_blkt_pump: float = 1.0
    """Pumping power for FW and Blanket multiplier factor"""


CREATE_DICTS_FROM_DATACLASS = PrimaryPumpingData
