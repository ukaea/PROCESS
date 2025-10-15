qmisc: float = None

qac: float = None

qcl: float = None

qss: float = None

p_shld_coolant_pump_elec_mw: float = None

p_div_coolant_pump_elec_mw: float = None

p_coolant_pump_total_mw: float = None

p_fw_blkt_heat_deposited_mw: float = None

p_fw_blkt_coolant_pump_elec_mw: float = None

p_blkt_breeder_pump_elec_mw: float = None

p_div_heat_deposited_mw: float = None

p_fw_heat_deposited_mw: float = None

p_blkt_heat_deposited_mw: float = None

p_blkt_liquid_breeder_heat_deposited_mw: float = None

p_shld_heat_deposited_mw: float = None

p_cp_coolant_pump_elec_mw: float = None

p_plant_core_systems_elec_mw: float = None

e_plant_net_electric_pulse_mj: float = None
"""Net electric energy output per pulse (MJ)"""

e_plant_net_electric_pulse_kwh: float = None
"""Net electric energy output per pulse (kWh)"""

f_p_div_primary_heat: float = None

delta_eta: float = None

i_div_primary_heat: float = None

p_turbine_loss_mw: float = None


def init_power_variables():
    global qmisc
    qmisc = 0.0

    global qac
    qac = 0.0

    global qcl
    qcl = 0.0

    global qss
    qss = 0.0

    global p_shld_coolant_pump_elec_mw
    p_shld_coolant_pump_elec_mw = 0.0

    global p_div_coolant_pump_elec_mw
    p_div_coolant_pump_elec_mw = 0.0

    global p_coolant_pump_total_mw
    p_coolant_pump_total_mw = 0.0

    global p_fw_blkt_heat_deposited_mw
    p_fw_blkt_heat_deposited_mw = 0.0

    global p_fw_blkt_coolant_pump_elec_mw
    p_fw_blkt_coolant_pump_elec_mw = 0.0

    global p_blkt_breeder_pump_elec_mw
    p_blkt_breeder_pump_elec_mw = 0.0

    global p_div_heat_deposited_mw
    p_div_heat_deposited_mw = 0.0

    global p_fw_heat_deposited_mw
    p_fw_heat_deposited_mw = 0.0

    global p_blkt_heat_deposited_mw
    p_blkt_heat_deposited_mw = 0.0

    global p_blkt_liquid_breeder_heat_deposited_mw
    p_blkt_liquid_breeder_heat_deposited_mw = 0.0

    global p_shld_heat_deposited_mw
    p_shld_heat_deposited_mw = 0.0

    global p_cp_coolant_pump_elec_mw
    p_cp_coolant_pump_elec_mw = 0.0

    global p_plant_core_systems_elec_mw
    p_plant_core_systems_elec_mw = 0.0

    global e_plant_net_electric_pulse_mj
    e_plant_net_electric_pulse_mj = 0.0

    global e_plant_net_electric_pulse_kwh
    e_plant_net_electric_pulse_kwh = 0.0

    global f_p_div_primary_heat
    f_p_div_primary_heat = 0.0

    global delta_eta
    delta_eta = 0.0

    global i_div_primary_heat
    i_div_primary_heat = 0.0

    global p_turbine_loss_mw
    p_turbine_loss_mw = 0.0
