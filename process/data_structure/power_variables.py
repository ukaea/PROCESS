from dataclasses import dataclass, field


@dataclass
class PowerData:
    qmisc: float = 0.0

    qac: float = 0.0

    qcl: float = 0.0

    qss: float = 0.0

    p_shld_coolant_pump_elec_mw: float = 0.0

    p_div_coolant_pump_elec_mw: float = 0.0

    p_coolant_pump_total_mw: float = 0.0

    p_fw_blkt_heat_deposited_mw: float = 0.0

    p_fw_blkt_coolant_pump_elec_mw: float = 0.0

    p_blkt_breeder_pump_elec_mw: float = 0.0

    p_div_heat_deposited_mw: float = 0.0

    p_fw_heat_deposited_mw: float = 0.0

    p_blkt_heat_deposited_mw: float = 0.0

    p_blkt_liquid_breeder_heat_deposited_mw: float = 0.0

    p_shld_heat_deposited_mw: float = 0.0

    p_cp_coolant_pump_elec_mw: float = 0.0

    p_plant_core_systems_elec_mw: float = 0.0

    e_plant_net_electric_pulse_mj: float = 0.0
    """Net electric energy output per pulse (MJ)"""

    e_plant_net_electric_pulse_kwh: float = 0.0
    """Net electric energy output per pulse (kWh)"""

    f_p_div_primary_heat: float = 0.0

    delta_eta: float = 0.0

    i_div_primary_heat: float = 0.0

    p_turbine_loss_mw: float = 0.0

    p_hcd_electric_total_profile_mw: list[float] = field(default_factory=list)
    """Profile of total HCD electric power (MW) over pulse"""

    p_tf_electric_supplies_profile_mw: list[float] = None
    """Profile of total TF coil electric power (MW) over pulse"""

    p_pf_electric_supplies_profile_mw: list[float] = None
    """Profile of total PF coil electric power (MW) over pulse"""

    p_coolant_pump_elec_total_profile_mw: list[float] = None
    """Profile of total coolant pump electric power (MW) over pulse"""

    vachtmw_profile_mw: list[float] = None
    """Profile of total active vacuum pump power (MW) over pulse"""

    p_tritium_plant_electric_profile_mw: list[float] = None
    """Profile of total tritium plant electric power (MW) over pulse"""

    p_cryo_plant_electric_profile_mw: list[float] = None
    """Profile of total cryo plant electric power (MW) over pulse"""

    p_plant_electric_base_total__profile_mw: list[float] = None
    """Profile of total plant electric base power (MW) over pulse"""

    p_plant_electric_gross_profile_mw: list[float] = None
    """Profile of total plant electric gross power (MW) over pulse"""

    p_plant_electric_net_profile_mw: list[float] = None
    """Profile of total plant electric net power (MW) over pulse"""

    p_fusion_total_profile_mw: list[float] = None
    """Profile of total fusion power (MW) over pulse"""


CREATE_DICTS_FROM_DATACLASS = PowerData
