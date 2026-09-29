"""Module containing variables for the primary pumping models"""

from dataclasses import dataclass

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class PrimaryPumpingData(PROCESSModelData):
    """Dataclass holding primary pumping variables"""

    gamma_he: Parameter[float] = 1.667
    """ratio of specific heats for helium (`i_p_coolant_pumping=3`)"""

    t_in_bb: Parameter[float] = 573.13
    """temperature in FW and blanket coolant at blanket entrance (`i_p_coolant_pumping=3`) [K]"""

    t_out_bb: Parameter[float] = 773.13
    """temperature in FW and blanket coolant at blanket exit (`i_p_coolant_pumping=3`) [K]"""

    p_he: Parameter[float] = 8.0e6
    """pressure in FW and blanket coolant at pump exit (`i_p_coolant_pumping=3`) [Pa]"""

    dp_he: Parameter[float] = 5.5e5
    """pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_fw_blkt: Parameter[float] = 1.5e5
    """pressure drop in FW and blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_fw: Parameter[float] = 1.5e5
    """pressure drop in FW coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_blkt: Parameter[float] = 3.5e3
    """pressure drop in blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    dp_liq: Parameter[float] = 1.0e7
    """pressure drop in liquid metal blanket coolant including heat exchanger and pipes (`i_p_coolant_pumping=3`) [Pa]"""

    p_fw_blkt_coolant_pump_mw: Parameter[float] = 0.0
    """mechanical pumping power for FW and blanket including heat exchanger and
    pipes (`i_p_coolant_pumping=3`) [MW]
    """

    f_p_fw_blkt_pump: Parameter[float] = 1.0
    """Pumping power for FW and Blanket multiplier factor"""


CREATE_DICTS_FROM_DATACLASS = PrimaryPumpingData
