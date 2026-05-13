from dataclasses import dataclass


@dataclass
class BuildingsData:
    admv: float = 1.0e5
    """administration building volume (m3)"""

    admvol: float = 0.0
    """volume of administration buildings (m3)"""

    aux_build_l: float = 60.0
    """aux building supporting tokamak processes length (m)"""

    aux_build_w: float = 30.0
    """aux building supporting tokamak processes width (m)"""

    aux_build_h: float = 5.0
    """aux building supporting tokamak processes height (m)"""

    auxcool_l: float = 20.0
    """Site-Wide Auxiliary Cooling Water facility length (m)"""

    auxcool_w: float = 20.0
    """Site-Wide Auxiliary Cooling Water facility width (m)"""

    auxcool_h: float = 5.0
    """Site-Wide Auxiliary Cooling Water facility height (m)"""

    bioshld_thk: float = 2.50
    """Radial thickness of bio-shield around reactor (m)"""

    chemlab_l: float = 50.0
    """Chemistry labs and treatment buldings length (m)"""

    chemlab_w: float = 30.0
    """Chemistry labs and treatment buldings width (m)"""

    chemlab_h: float = 6.0
    """Chemistry labs and treatment buldings height (m)"""

    dz_tf_cryostat: float = 2.5
    """vertical clearance from TF coil to cryostat (m) (calculated for tokamaks)"""

    clh2: float = 15.0
    """clearance beneath TF coil to foundation (including basement) (m)"""

    control_buildings_l: float = 80.0
    """control building length (m)"""

    control_buildings_w: float = 60.0
    """control building width (m)"""

    control_buildings_h: float = 6.0
    """control building height (m)"""

    conv: float = 6.0e4
    """control building volume (m3)"""

    convol: float = 0.0
    """volume of control, protection and i&c building (m3)"""

    crane_arm_h: float = 10.0
    """vertical dimension of crane arm, operating over reactor (m)"""

    crane_clrnc_h: float = 4.0
    """horizontal clearance to building wall for crane operation (m)"""

    crane_clrnc_v: float = 3.0
    """vertical clearance for crane operation (m)"""

    cryomag_l: float = 120.0
    """Cryogenic Buildings for Magnet and Fuel Cycle length (m)"""

    cryomag_w: float = 90.0
    """Cryogenic Buildings for Magnet and Fuel Cycle width (m)"""

    cryomag_h: float = 5.0
    """Cryogenic Buildings for Magnet and Fuel Cycle height (m)"""

    cryostore_l: float = 160.0
    """Magnet Cryo Storage Tanks length, width, height (m)"""

    cryostore_w: float = 30.0
    """Magnet Cryo Storage Tanks length, width, height (m)"""

    cryostore_h: float = 20.0
    """Magnet Cryo Storage Tanks length, width, height (m)"""

    cryostat_clrnc: float = 2.5
    """vertical clearance from TF coil to cryostat (m)"""

    cryvol: float = 0.0
    """volume of cryoplant building (m3)"""

    a_plant_floor_effective: float = 0.0
    """effective total floor space (m2)"""

    elecdist_l: float = 380.0
    """Transformers and electrical distribution facilities length (m)"""

    elecdist_w: float = 350.0
    """Transformers and electrical distribution facilities width (m)"""

    elecdist_h: float = 5.0
    """Transformers and electrical distribution facilities height (m)"""

    elecload_l: float = 100.0
    """Electric (eesential and non-essential) load centres length (m)"""

    elecload_w: float = 90.0
    """Electric (eesential and non-essential) load centres width (m)"""

    elecload_h: float = 3.0
    """Electric (eesential and non-essential) load centres height (m)"""

    elecstore_l: float = 100.0
    """Energy Storage facilities length (m)"""

    elecstore_w: float = 60.0
    """Energy Storage facilities width (m)"""

    elecstore_h: float = 12.0
    """Energy Storage facilities height (m)"""

    elevol: float = 0.0
    """volume of electrical equipment building (m3)"""

    esbldgm3: float = 1.0e3
    """volume of energy storage equipment building (m3) (not used if `i_pulsed_plant=0`)"""

    fc_building_l: float = 60.0
    """Fuel Cycle facilities length (m)"""

    fc_building_w: float = 60.0
    """Fuel Cycle facilities width (m)"""

    fndt: float = 2.0
    """foundation thickness (m)"""

    gas_buildings_l: float = 25.0
    """air & gas supply (amalgamated) buildings length (m)"""

    gas_buildings_w: float = 15.0
    """air & gas supply (amalgamated) buildings width (m)"""

    gas_buildings_h: float = 5.0
    """air & gas supply (amalgamated) buildings height (m)"""

    ground_clrnc: float = 5.0
    """clearance beneath TF coil (m)"""

    hcd_building_l: float = 70.0
    """HCD building length (m)"""

    hcd_building_w: float = 40.0
    """HCD building width (m)"""

    hcd_building_h: float = 25.0
    """HCD building height (m)"""

    hccl: float = 5.0
    """clearance around components in hot cell (m)"""

    hcwt: float = 1.5
    """hot cell wall thickness (m)"""

    heat_sink_l: float = 160.0
    """heat sinks length (m)"""

    heat_sink_w: float = 80.0
    """heat sinks width (m)"""

    heat_sink_h: float = 12.0
    """heat sinks height (m)"""

    hot_sepdist: float = 2.0
    """hot cell storage component separation distance (m)"""

    hotcell_h: float = 12.0
    """hot cell storage and maintenance facility height (m)"""

    hw_storage_l: float = 20.0
    """hazardous waste storage building length, width, height (m)"""

    hw_storage_w: float = 10.0
    """hazardous waste storage building length, width, height (m)"""

    hw_storage_h: float = 5.0
    """hazardous waste storage building length, width, height (m)"""

    i_bldgs_size: int = 0
    """switch between routines estimating building sizes (0 = default; 1 = updated)"""

    i_bldgs_v: int = 0
    """switch to select verbose output for buildings (1 = verbose)"""

    ilw_smelter_l: float = 50.0
    """radioactive waste smelting facility length (m)"""

    ilw_smelter_w: float = 30.0
    """radioactive waste smelting facility width (m)"""

    ilw_smelter_h: float = 30.0
    """radioactive waste smelting facility height (m)"""

    ilw_storage_l: float = 120.0
    """ILW waste storage building length (m)"""

    ilw_storage_w: float = 100.0
    """ILW waste storage building width (m)"""

    ilw_storage_h: float = 8.0
    """ILW waste storage building height (m)"""

    llw_storage_l: float = 45.0
    """LLW waste storage building length (m)"""

    llw_storage_w: float = 20.0
    """LLW waste storage building width (m)"""

    llw_storage_h: float = 5.0
    """LLW waste storage building height (m)"""

    magnet_pulse_l: float = 105.0
    """pulsed magnet power building length (m)"""

    magnet_pulse_w: float = 40.0
    """pulsed magnet power building width (m)"""

    magnet_pulse_h: float = 5.0
    """pulsed magnet power building height (m)"""

    magnet_trains_l: float = 120.0
    """steady state magnet power trains building length (m)"""

    magnet_trains_w: float = 90.0
    """steady state magnet power trains building width (m)"""

    magnet_trains_h: float = 5.0
    """steady state magnet power trains building height (m)"""

    maint_cont_l: float = 125.0
    """maintenance control building length (m)"""

    maint_cont_w: float = 100.0
    """maintenance control building width (m)"""

    maint_cont_h: float = 6.0
    """maintenance control building height (m)"""

    mbvfac: float = 2.8
    """maintenance building volume multiplication factor"""

    nbi_sys_l: float = 225.0
    """NBI system length, width (m)"""

    nbi_sys_w: float = 185.0
    """NBI system width (m)"""

    pfbldgm3: float = 2.0e4
    """volume of PF coil power supply building (m3)"""

    pibv: float = 2.0e4
    """power injection building volume (m3)"""

    qnty_sfty_fac: float = 2.0
    """quantity safety factor for component use during plant lifetime"""

    rbvfac: float = 1.6
    """reactor building volume multiplication factor"""

    rbrt: float = 1.0
    """reactor building roof thickness (m)"""

    rbvol: float = 0.0
    """reactor building volume (m3)"""

    rbwt: float = 2.0
    """reactor building wall thickness (m)"""

    reactor_clrnc: float = 4.0
    """clearance around reactor (m)"""

    reactor_fndtn_thk: float = 2.0
    """reactor building foundation thickness (m)"""

    reactor_hall_l: float = 0.0
    """reactor building length (m)"""

    reactor_hall_w: float = 0.0
    """reactor building width (m)"""

    reactor_hall_h: float = 0.0
    """reactor building height (m)"""

    reactor_roof_thk: float = 1.0
    """reactor building roof thickness (m)"""

    reactor_wall_thk: float = 2.0
    """reactor building wall thickness (m)"""

    rmbvol: float = 0.0
    """volume of maintenance and assembly building (m3)"""

    robotics_l: float = 50.0
    """robotics buildings length (m)"""

    robotics_w: float = 30.0
    """robotics buildings width (m)"""

    robotics_h: float = 30.0
    """robotics buildings height (m)"""

    row: float = 4.0
    """clearance to building wall for crane operation (m)"""

    rxcl: float = 4.0
    """clearance around reactor (m)"""

    sec_buildings_l: float = 30.0
    """security & safety buildings length (m)"""

    sec_buildings_w: float = 25.0
    """security & safety buildings width (m)"""

    sec_buildings_h: float = 6.0
    """security & safety buildings height (m)"""

    shmf: float = 0.5
    """fraction of shield mass per TF coil to be moved in the maximum shield lift"""

    shov: float = 1.0e5
    """shops and warehouse volume (m3)"""

    shovol: float = 0.0
    """volume of shops and buildings for plant auxiliaries (m3)"""

    staff_buildings_area: float = 4.8e5
    """footprint of staff buildings (m2)"""

    staff_buildings_h: float = 5.0
    """staff buildings height (m)"""

    stcl: float = 3.0
    """clearance above crane to roof (m)"""

    tfcbv: float = 2.0e4
    """volume of TF coil power supply building (m3) (calculated if TF coils are superconducting)"""

    transp_clrnc: float = 1.0
    """transportation clearance between components (m)"""

    trcl: float = 1.0
    """transportation clearance between components (m)"""

    triv: float = 4.0e4
    """volume of tritium, fuel handling and health physics buildings (m3)"""

    turbine_hall_l: float = 109.0
    """turbine hall length (m)"""

    turbine_hall_w: float = 62.0
    """turbine hall width (m)"""

    turbine_hall_h: float = 15.0
    """turbine hall height (m)"""

    tw_storage_l: float = 90.0
    """tritiated waste storage building length (m)"""

    tw_storage_w: float = 30.0
    """tritiated waste storage building width (m)"""

    tw_storage_h: float = 5.0
    """tritiated waste storage building height (m)"""

    volrci: float = 0.0
    """internal volume of reactor building (m3)"""

    volnucb: float = 0.0
    """sum of nuclear buildings volumes (m3)"""

    warm_shop_l: float = 100.0
    """warm shop length (m)"""

    warm_shop_w: float = 50.0
    """warm shop width (m)"""

    warm_shop_h: float = 10.0
    """warm shop height (m)"""

    water_buildings_l: float = 110.0
    """water, laundry & drainage buildings length (m)"""

    water_buildings_w: float = 10.0
    """water, laundry & drainage buildings width (m)"""

    water_buildings_h: float = 5.0
    """water, laundry & drainage buildings height (m)"""

    wgt: float = 5.0e5
    """reactor building crane capacity (kg) (calculated if 0 is input)"""

    wgt2: float = 1.0e5
    """hot cell crane capacity (kg) (calculated if 0 is input)"""

    workshop_l: float = 150.0
    """[cold] workshop buildings length (m)"""

    workshop_w: float = 125.0
    """[cold] workshop buildings width (m)"""

    workshop_h: float = 10.0
    """[cold] workshop buildings height (m)"""

    wrbi: float = 0.0
    """distance from centre of machine to building wall (m)"""

    wsvol: float = 0.0
    """volume of warm shop building (m3)"""

    wsvfac: float = 1.9
    """warm shop building volume multiplication factor"""

    a_reactor_bldg: float = 8.32e3
    """Floor area of reactor building in m^2"""

    a_ee_ps_bldg: float = 2.133e4
    """Floor area of electrical equipment and power supply building in m^2"""

    a_aux_services_bldg: float = 1.0e3
    """Floor area of auxiliary services building in m^2"""

    a_hot_cell_bldg: float = 8.43e3
    """Floor area of hot cell building in m^2"""

    a_reactor_service_bldg: float = 2.44e3
    """Floor area of reactor service building in m^2"""

    a_service_water_bldg: float = 1.567e3
    """Floor area of service water building in m^2"""

    a_fuel_handling_bldg: float = 1.67e3
    """Floor area of fuel handling and storage building in m^2"""

    a_control_room_bldg: float = 2.88e3
    """Floor area of controlroom building in m^2"""

    a_ac_ps_bldg: float = 6.423e3
    """Floor area of AC power supply building in m^2"""

    a_admin_bldg: float = 2.5674e4
    """Floor area of admin building in m^2"""

    a_site_service_bldg: float = 8.3e3
    """Floor area of site service building in m^2"""

    a_cryo_inert_gas_bldg: float = 1.838e4
    """Floor area of cryogenics and inert gas storage building in m^2"""

    a_security_bldg: float = 4.552e3
    """Floor area of security building in m^2"""


CREATE_DICTS_FROM_DATACLASS = BuildingsData
