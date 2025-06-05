qmisc: float

qac: float

qcl: float

qss: float

p_shld_coolant_pump_elec_mw: float

p_div_coolant_pump_elec_mw: float

p_coolant_pump_total_mw: float

p_fw_blkt_heat_deposited_mw: float

p_fw_blkt_coolant_pump_elec_mw: float

p_blkt_breeder_pump_elec_mw: float

p_div_heat_deposited_mw: float

p_fw_heat_deposited_mw: float

p_blkt_heat_deposited_mw: float

p_blkt_liquid_breeder_heat_deposited_mw: float

p_shld_heat_deposited_mw: float

p_cp_coolant_pump_elec_mw: float

p_plant_core_systems_elec_mw: float

f_p_div_primary_heat: float

delta_eta: float

i_div_primary_heat: float

p_turbine_loss_mw: float


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

    global f_p_div_primary_heat
    f_p_div_primary_heat = 0.0

    global delta_eta
    delta_eta = 0.0

    global i_div_primary_heat
    i_div_primary_heat = 0.0

    global p_turbine_loss_mw
    p_turbine_loss_mw = 0.0
