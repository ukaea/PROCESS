"""
This module contains global variables relating to the heat transport system
of a fusion power plant, and also those for a hydrogen production plant.
### References
-
"""

p_plant_electric_base: float = None
"""base plant electric load (W)"""


p_cryo_plant_electric_mw: float = None
"""cryogenic plant power (MW)"""


p_cryo_plant_electric_max_mw: float = None
"""Maximum cryogenic plant power (MW)
Constraint equation icc = 87
Scan variable nwseep = 56
"""


etatf: float = None
"""AC to resistive power conversion for TF coils"""


eta_turbine: float = None
"""thermal to electric conversion efficiency if `i_thermal_electric_conversion=2`; otherwise calculated."""


etath_liq: float = None


fachtmw: float = None
"""facility heat removal (MW)"""


p_plant_electric_base_total_mw: float = None
"""total baseline power required at all times (MW)"""


fgrosbop: float = None
"""scaled fraction of gross power to balance-of-plant"""


fmgdmw: float = None
"""power to mgf (motor-generator flywheel) units (MW) (ignored if `i_pf_energy_storage_source=2`)"""


f_p_blkt_coolant_pump_total_heat: float = None
"""fraction of total blanket thermal power required to drive the blanket
coolant pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)
"""


f_p_div_coolant_pump_total_heat: float = None
"""fraction of total divertor thermal power required to drive the divertor
coolant pumps (default assumes water coolant)
"""


f_p_fw_coolant_pump_total_heat: float = None
"""fraction of total first wall thermal power required to drive the FW coolant
pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)
"""


f_p_shld_coolant_pump_total_heat: float = None
"""fraction of total shield thermal power required to drive the shield coolant
pumps (default assumes water coolant)
"""


helpow: float = None
"""Heat removal at cryogenic temperature temp_tf_cryo (W)"""


helpow_cryal: float = None
"""Heat removal at cryogenic temperature temp_cp_coolant_inlet (W)"""


p_coolant_pump_elec_total_mw: float = None
"""heat transport system electrical pump power (MW)"""


p_blkt_coolant_pump_mw: float = None
"""blanket primary coolant mechanical pumping power (MW)"""


p_blkt_breeder_pump_mw: float = None
"""blanket secondary coolant mechanical pumping power (MW)"""


htpmw_blkt_tot: float = None
"""blanket primary + secondary coolant mechanical pumping power (MW)"""


p_div_coolant_pump_mw: float = None
"""divertor coolant mechanical pumping power (MW)"""


p_fw_coolant_pump_mw: float = None
"""first wall coolant mechanical pumping power (MW)"""


p_shld_coolant_pump_mw: float = None
"""shield and vacuum vessel coolant mechanical pumping power (MW)"""


p_coolant_pump_loss_total_mw: float = None
"""Waste power lost from primary coolant pumps (MW)"""


ipowerflow: int = None
"""switch for power flow model:
- =0 pre-2014 version
- =1 comprehensive 2014 model
"""


i_shld_primary_heat: int = None
"""Switch for shield thermal power destiny:
- =0 does not contribute to energy generation cycle
- =1 contributes to energy generation cycle
"""


n_primary_heat_exchangers: int = None
"""number of primary heat exchangers"""


pacpmw: float = None
"""total pulsed power system load (MW)"""


peakmva: float = None
"""peak MVA requirement"""


p_fw_div_heat_deposited_mw: float = None
"""heat removal from first wall/divertor (MW)"""


p_plant_electric_gross_mw: float = None
"""gross electric power (MW)"""


p_hcd_electric_loss_mw: float = None
"""power dissipated in heating and current drive system (MW)"""


p_hcd_electric_total_mw: float = None
"""injector wall plug power (MW)"""


p_hcd_secondary_electric_mw: float = None
"""Secondary HCD system injector wall plug power (MW)"""


p_hcd_primary_electric_mw: float = None
"""Primary HCD system injector wall plug power (MW)"""


p_plant_electric_net_mw: float = None
"""net electric power (MW)"""


p_plant_electric_recirc_mw: float = None
"""recirculating electric power (MW)"""


priheat: float = None
"""total thermal power removed from fusion core (MW)"""


p_div_secondary_heat_mw: float = None
"""Low-grade heat lost in divertor (MW)"""


p_hcd_secondary_heat_mw: float = None
"""Low-grade heat lost into HCD apparatus (MW)"""


p_plant_secondary_heat_mw: float = None
"""Low-grade heat (MW)"""


p_shld_secondary_heat_mw: float = None
"""Low-grade heat deposited in shield (MW)"""


p_plant_primary_heat_mw: float = None
"""High-grade heat useful for electric production (MW)"""


pflux_plant_floor_electric: float = None
"""base AC power requirement per unit floor area (W/m2)"""


p_tf_electric_supplies_mw: float = None
"""total steady state TF coil AC power demand (MW)"""


tlvpmw: float = None
"""estimate of total low voltage power (MW)"""


p_tritium_plant_electric_mw: float = None
"""power required for tritium processing (MW)"""


temp_turbine_coolant_in: float = None
"""coolant temperature at turbine inlet (K) (`i_thermal_electric_conversion = 3,4`)"""


vachtmw: float = None
"""vacuum pump power (MW)"""


f_p_plant_electric_recirc: float = None
"""fraction of recirculating electric power to total electric power"""


def init_heat_transport_variables():
    """Initialise heat transport variables"""
    global \
        p_plant_electric_base, \
        p_cryo_plant_electric_mw, \
        p_cryo_plant_electric_max_mw, \
        etatf, \
        eta_turbine, \
        etath_liq, \
        fachtmw, \
        p_plant_electric_base_total_mw, \
        fgrosbop, \
        fmgdmw, \
        f_p_blkt_coolant_pump_total_heat, \
        f_p_div_coolant_pump_total_heat, \
        f_p_fw_coolant_pump_total_heat, \
        f_p_shld_coolant_pump_total_heat, \
        helpow, \
        helpow_cryal, \
        p_coolant_pump_elec_total_mw, \
        p_blkt_coolant_pump_mw, \
        p_blkt_breeder_pump_mw, \
        htpmw_blkt_tot, \
        p_div_coolant_pump_mw, \
        p_fw_coolant_pump_mw, \
        p_shld_coolant_pump_mw, \
        p_coolant_pump_loss_total_mw, \
        ipowerflow, \
        i_shld_primary_heat, \
        n_primary_heat_exchangers, \
        pacpmw, \
        peakmva, \
        p_fw_div_heat_deposited_mw, \
        p_plant_electric_gross_mw, \
        p_hcd_electric_loss_mw, \
        p_hcd_electric_total_mw, \
        p_hcd_secondary_electric_mw, \
        p_plant_electric_net_mw, \
        p_plant_electric_recirc_mw, \
        priheat, \
        p_div_secondary_heat_mw, \
        p_hcd_secondary_heat_mw, \
        p_plant_secondary_heat_mw, \
        p_shld_secondary_heat_mw, \
        p_plant_primary_heat_mw, \
        pflux_plant_floor_electric, \
        p_tf_electric_supplies_mw, \
        tlvpmw, \
        p_tritium_plant_electric_mw, \
        temp_turbine_coolant_in, \
        vachtmw, \
        f_p_plant_electric_recirc

    p_plant_electric_base = 5.0e6
    p_cryo_plant_electric_mw = 0.0
    p_cryo_plant_electric_max_mw = 50.0
    etatf = 0.9
    eta_turbine = 0.35
    etath_liq = 0.35
    fachtmw = 0.0
    p_plant_electric_base_total_mw = 0.0
    fgrosbop = 0.0
    fmgdmw = 0.0
    f_p_blkt_coolant_pump_total_heat = 0.005
    f_p_div_coolant_pump_total_heat = 0.005
    f_p_fw_coolant_pump_total_heat = 0.005
    f_p_shld_coolant_pump_total_heat = 0.005
    helpow = 0.0
    helpow_cryal = 0.0
    p_coolant_pump_elec_total_mw = 0.0
    p_blkt_coolant_pump_mw = 0.0
    p_blkt_breeder_pump_mw = 0.0
    htpmw_blkt_tot = 0.0
    p_div_coolant_pump_mw = 0.0
    p_fw_coolant_pump_mw = 0.0
    p_shld_coolant_pump_mw = 0.0
    p_coolant_pump_loss_total_mw = 0.0
    ipowerflow = 1
    i_shld_primary_heat = 1
    n_primary_heat_exchangers = 0
    pacpmw = 0.0
    peakmva = 0.0
    p_fw_div_heat_deposited_mw = 0.0
    p_plant_electric_gross_mw = 0.0
    p_hcd_electric_loss_mw = 0.0
    p_hcd_electric_total_mw = 0.0
    p_hcd_secondary_electric_mw = 0.0
    p_plant_electric_net_mw = 0.0
    p_plant_electric_recirc_mw = 0.0
    priheat = 0.0
    p_div_secondary_heat_mw = 0.0
    p_hcd_secondary_heat_mw = 0.0
    p_plant_secondary_heat_mw = 0.0
    p_shld_secondary_heat_mw = 0.0
    p_plant_primary_heat_mw = 0.0
    pflux_plant_floor_electric = 150.0
    p_tf_electric_supplies_mw = 0.0
    tlvpmw = 0.0
    p_tritium_plant_electric_mw = 15.0
    temp_turbine_coolant_in = 0.0
    vachtmw = 0.5
    f_p_plant_electric_recirc = 0.0
