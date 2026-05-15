"""
This module contains global variables relating to the heat transport system
of a fusion power plant, and also those for a hydrogen production plant.
### References
-
"""

from dataclasses import dataclass


@dataclass
class HeatTransportData:
    p_plant_electric_base: float = 5.0e6
    """base plant electric load (W)"""

    p_cryo_plant_electric_mw: float = 0.0
    """cryogenic plant power (MW)"""

    p_cryo_plant_electric_max_mw: float = 50.0
    """Maximum cryogenic plant power (MW)
    Constraint equation icc = 87
    Scan variable nwseep = 56
    """

    etatf: float = 0.9
    """AC to resistive power conversion for TF coils"""

    eta_turbine: float = 0.35
    """thermal to electric conversion efficiency if `i_thermal_electric_conversion=2`; otherwise calculated."""

    etath_liq: float = 0.35

    fachtmw: float = 0.0
    """facility heat removal (MW)"""

    p_plant_electric_base_total_mw: float = 0.0
    """total baseline power required at all times (MW)"""

    fgrosbop: float = 0.0
    """scaled fraction of gross power to balance-of-plant"""

    fmgdmw: float = 0.0
    """power to mgf (motor-generator flywheel) units (MW) (ignored if `i_pf_energy_storage_source=2`)"""

    f_p_blkt_coolant_pump_total_heat: float = 0.005
    """fraction of total blanket thermal power required to drive the blanket
    coolant pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)
    """

    f_p_div_coolant_pump_total_heat: float = 0.005
    """fraction of total divertor thermal power required to drive the divertor
    coolant pumps (default assumes water coolant)
    """

    f_p_fw_coolant_pump_total_heat: float = 0.005
    """fraction of total first wall thermal power required to drive the FW coolant
    pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)
    """

    f_p_shld_coolant_pump_total_heat: float = 0.005
    """fraction of total shield thermal power required to drive the shield coolant
    pumps (default assumes water coolant)
    """

    helpow: float = 0.0
    """Heat removal at cryogenic temperature temp_tf_cryo (W)"""

    helpow_cryal: float = 0.0
    """Heat removal at cryogenic temperature temp_cp_coolant_inlet (W)"""

    p_coolant_pump_elec_total_mw: float = 0.0
    """heat transport system electrical pump power (MW)"""

    p_blkt_coolant_pump_mw: float = 0.0
    """blanket primary coolant mechanical pumping power (MW)"""

    p_blkt_breeder_pump_mw: float = 0.0
    """blanket secondary coolant mechanical pumping power (MW)"""

    htpmw_blkt_tot: float = 0.0
    """blanket primary + secondary coolant mechanical pumping power (MW)"""

    p_div_coolant_pump_mw: float = 0.0
    """divertor coolant mechanical pumping power (MW)"""

    p_fw_coolant_pump_mw: float = 0.0
    """first wall coolant mechanical pumping power (MW)"""

    p_shld_coolant_pump_mw: float = 0.0
    """shield and vacuum vessel coolant mechanical pumping power (MW)"""

    p_coolant_pump_loss_total_mw: float = 0.0
    """Waste power lost from primary coolant pumps (MW)"""

    ipowerflow: int = 1
    """switch for power flow model:
    - =0 pre-2014 version
    - =1 comprehensive 2014 model
    """

    i_shld_primary_heat: int = 1
    """Switch for shield thermal power destiny:
    - =0 does not contribute to energy generation cycle
    - =1 contributes to energy generation cycle
    """

    n_primary_heat_exchangers: int = 0
    """number of primary heat exchangers"""

    pacpmw: float = 0.0
    """total pulsed power system load (MW)"""

    peakmva: float = 0.0
    """peak MVA requirement"""

    p_fw_div_heat_deposited_mw: float = 0.0
    """heat removal from first wall/divertor (MW)"""

    p_plant_electric_gross_mw: float = 0.0
    """gross electric power (MW)"""

    p_hcd_electric_loss_mw: float = 0.0
    """power dissipated in heating and current drive system (MW)"""

    p_hcd_electric_total_mw: float = 0.0
    """injector wall plug power (MW)"""

    p_hcd_secondary_electric_mw: float = 0.0
    """Secondary HCD system injector wall plug power (MW)"""

    p_hcd_primary_electric_mw: float = None
    """Primary HCD system injector wall plug power (MW)"""

    p_plant_electric_net_mw: float = 0.0
    """net electric power (MW)"""

    p_plant_electric_recirc_mw: float = 0.0
    """recirculating electric power (MW)"""

    priheat: float = 0.0
    """total thermal power removed from fusion core (MW)"""

    p_div_secondary_heat_mw: float = 0.0
    """Low-grade heat lost in divertor (MW)"""

    p_hcd_secondary_heat_mw: float = 0.0
    """Low-grade heat lost into HCD apparatus (MW)"""

    p_plant_secondary_heat_mw: float = 0.0
    """Low-grade heat (MW)"""

    p_shld_secondary_heat_mw: float = 0.0
    """Low-grade heat deposited in shield (MW)"""

    p_plant_primary_heat_mw: float = 0.0
    """High-grade heat useful for electric production (MW)"""

    pflux_plant_floor_electric: float = 150.0
    """base AC power requirement per unit floor area (W/m2)"""

    p_tf_electric_supplies_mw: float = 0.0
    """total steady state TF coil AC power demand (MW)"""

    tlvpmw: float = 0.0
    """estimate of total low voltage power (MW)"""

    p_tritium_plant_electric_mw: float = 15.0
    """power required for tritium processing (MW)"""

    temp_turbine_coolant_in: float = 0.0
    """coolant temperature at turbine inlet (K) (`i_thermal_electric_conversion = 3,4`)"""

    vachtmw: float = 0.5
    """vacuum pump power (MW)"""

    f_p_plant_electric_recirc: float = 0.0
    """fraction of recirculating electric power to total electric power"""


CREATE_DICTS_FROM_DATACLASS = HeatTransportData
