admv: float = None
"""administration building volume (m3)"""


admvol: float = None
"""volume of administration buildings (m3)"""


aux_build_l: float = None
"""aux building supporting tokamak processes length (m)"""


aux_build_w: float = None
"""aux building supporting tokamak processes width (m)"""


aux_build_h: float = None
"""aux building supporting tokamak processes height (m)"""


auxcool_l: float = None
"""Site-Wide Auxiliary Cooling Water facility length (m)"""


auxcool_w: float = None
"""Site-Wide Auxiliary Cooling Water facility width (m)"""


auxcool_h: float = None
"""Site-Wide Auxiliary Cooling Water facility height (m)"""


bioshld_thk: float = None
"""Radial thickness of bio-shield around reactor (m)"""


chemlab_l: float = None
"""Chemistry labs and treatment buldings length (m)"""


chemlab_w: float = None
"""Chemistry labs and treatment buldings width (m)"""


chemlab_h: float = None
"""Chemistry labs and treatment buldings height (m)"""


dz_tf_cryostat: float = None
"""vertical clearance from TF coil to cryostat (m) (calculated for tokamaks)"""


clh2: float = None
"""clearance beneath TF coil to foundation (including basement) (m)"""


control_buildings_l: float = None
"""control building length (m)"""


control_buildings_w: float = None
"""control building width (m)"""


control_buildings_h: float = None
"""control building height (m)"""


conv: float = None
"""control building volume (m3)"""


convol: float = None
"""volume of control, protection and i&c building (m3)"""


crane_arm_h: float = None
"""vertical dimension of crane arm, operating over reactor (m)"""


crane_clrnc_h: float = None
"""horizontal clearance to building wall for crane operation (m)"""


crane_clrnc_v: float = None
"""vertical clearance for crane operation (m)"""


cryomag_l: float = None
"""Cryogenic Buildings for Magnet and Fuel Cycle length (m)"""


cryomag_w: float = None
"""Cryogenic Buildings for Magnet and Fuel Cycle width (m)"""


cryomag_h: float = None
"""Cryogenic Buildings for Magnet and Fuel Cycle height (m)"""


cryostore_l: float = None
"""Magnet Cryo Storage Tanks length, width, height (m)"""


cryostore_w: float = None
"""Magnet Cryo Storage Tanks length, width, height (m)"""


cryostore_h: float = None
"""Magnet Cryo Storage Tanks length, width, height (m)"""


cryostat_clrnc: float = None
"""vertical clearance from TF coil to cryostat (m)"""


cryvol: float = None
"""volume of cryoplant building (m3)"""


a_plant_floor_effective: float = None
"""effective total floor space (m2)"""


elecdist_l: float = None
"""Transformers and electrical distribution facilities length (m)"""


elecdist_w: float = None
"""Transformers and electrical distribution facilities width (m)"""


elecdist_h: float = None
"""Transformers and electrical distribution facilities height (m)"""


elecload_l: float = None
"""Electric (eesential and non-essential) load centres length (m)"""


elecload_w: float = None
"""Electric (eesential and non-essential) load centres width (m)"""


elecload_h: float = None
"""Electric (eesential and non-essential) load centres height (m)"""


elecstore_l: float = None
"""Energy Storage facilities length (m)"""


elecstore_w: float = None
"""Energy Storage facilities width (m)"""


elecstore_h: float = None
"""Energy Storage facilities height (m)"""


elevol: float = None
"""volume of electrical equipment building (m3)"""


esbldgm3: float = None
"""volume of energy storage equipment building (m3) (not used if `i_pulsed_plant=0`)"""


fc_building_l: float = None
"""Fuel Cycle facilities length (m)"""


fc_building_w: float = None
"""Fuel Cycle facilities width (m)"""


fndt: float = None
"""foundation thickness (m)"""


gas_buildings_l: float = None
"""air & gas supply (amalgamated) buildings length (m)"""


gas_buildings_w: float = None
"""air & gas supply (amalgamated) buildings width (m)"""


gas_buildings_h: float = None
"""air & gas supply (amalgamated) buildings height (m)"""


ground_clrnc: float = None
"""clearance beneath TF coil (m)"""


hcd_building_l: float = None
"""HCD building length (m)"""


hcd_building_w: float = None
"""HCD building width (m)"""


hcd_building_h: float = None
"""HCD building height (m)"""


hccl: float = None
"""clearance around components in hot cell (m)"""


hcwt: float = None
"""hot cell wall thickness (m)"""


heat_sink_l: float = None
"""heat sinks length (m)"""


heat_sink_w: float = None
"""heat sinks width (m)"""


heat_sink_h: float = None
"""heat sinks height (m)"""


hot_sepdist: float = None
"""hot cell storage component separation distance (m)"""


hotcell_h: float = None
"""hot cell storage and maintenance facility height (m)"""


hw_storage_l: float = None
"""hazardous waste storage building length, width, height (m)"""


hw_storage_w: float = None
"""hazardous waste storage building length, width, height (m)"""


hw_storage_h: float = None
"""hazardous waste storage building length, width, height (m)"""


i_bldgs_size: int = None
"""switch between routines estimating building sizes (0 = default; 1 = updated)"""


i_bldgs_v: int = None
"""switch to select verbose output for buildings (1 = verbose)"""


ilw_smelter_l: float = None
"""radioactive waste smelting facility length (m)"""


ilw_smelter_w: float = None
"""radioactive waste smelting facility width (m)"""


ilw_smelter_h: float = None
"""radioactive waste smelting facility height (m)"""


ilw_storage_l: float = None
"""ILW waste storage building length (m)"""


ilw_storage_w: float = None
"""ILW waste storage building width (m)"""


ilw_storage_h: float = None
"""ILW waste storage building height (m)"""


llw_storage_l: float = None
"""LLW waste storage building length (m)"""


llw_storage_w: float = None
"""LLW waste storage building width (m)"""


llw_storage_h: float = None
"""LLW waste storage building height (m)"""


magnet_pulse_l: float = None
"""pulsed magnet power building length (m)"""


magnet_pulse_w: float = None
"""pulsed magnet power building width (m)"""


magnet_pulse_h: float = None
"""pulsed magnet power building height (m)"""


magnet_trains_l: float = None
"""steady state magnet power trains building length (m)"""


magnet_trains_w: float = None
"""steady state magnet power trains building width (m)"""


magnet_trains_h: float = None
"""steady state magnet power trains building height (m)"""


maint_cont_l: float = None
"""maintenance control building length (m)"""


maint_cont_w: float = None
"""maintenance control building width (m)"""


maint_cont_h: float = None
"""maintenance control building height (m)"""


mbvfac: float = None
"""maintenance building volume multiplication factor"""


nbi_sys_l: float = None
"""NBI system length, width (m)"""


nbi_sys_w: float = None
"""NBI system width (m)"""


pfbldgm3: float = None
"""volume of PF coil power supply building (m3)"""


pibv: float = None
"""power injection building volume (m3)"""


qnty_sfty_fac: float = None
"""quantity safety factor for component use during plant lifetime"""


rbvfac: float = None
"""reactor building volume multiplication factor"""


rbrt: float = None
"""reactor building roof thickness (m)"""


rbvol: float = None
"""reactor building volume (m3)"""


rbwt: float = None
"""reactor building wall thickness (m)"""


reactor_clrnc: float = None
"""clearance around reactor (m)"""


reactor_fndtn_thk: float = None
"""reactor building foundation thickness (m)"""


reactor_hall_l: float = None
"""reactor building length (m)"""


reactor_hall_w: float = None
"""reactor building width (m)"""


reactor_hall_h: float = None
"""reactor building height (m)"""


reactor_roof_thk: float = None
"""reactor building roof thickness (m)"""


reactor_wall_thk: float = None
"""reactor building wall thickness (m)"""


rmbvol: float = None
"""volume of maintenance and assembly building (m3)"""


robotics_l: float = None
"""robotics buildings length (m)"""


robotics_w: float = None
"""robotics buildings width (m)"""


robotics_h: float = None
"""robotics buildings height (m)"""


row: float = None
"""clearance to building wall for crane operation (m)"""


rxcl: float = None
"""clearance around reactor (m)"""


sec_buildings_l: float = None
"""security & safety buildings length (m)"""


sec_buildings_w: float = None
"""security & safety buildings width (m)"""


sec_buildings_h: float = None
"""security & safety buildings height (m)"""


shmf: float = None
"""fraction of shield mass per TF coil to be moved in the maximum shield lift"""


shov: float = None
"""shops and warehouse volume (m3)"""


shovol: float = None
"""volume of shops and buildings for plant auxiliaries (m3)"""


staff_buildings_area: float = None
"""footprint of staff buildings (m2)"""


staff_buildings_h: float = None
"""staff buildings height (m)"""


stcl: float = None
"""clearance above crane to roof (m)"""


tfcbv: float = None
"""volume of TF coil power supply building (m3) (calculated if TF coils are superconducting)"""


transp_clrnc: float = None
"""transportation clearance between components (m)"""


trcl: float = None
"""transportation clearance between components (m)"""


triv: float = None
"""volume of tritium, fuel handling and health physics buildings (m3)"""


turbine_hall_l: float = None
"""turbine hall length (m)"""


turbine_hall_w: float = None
"""turbine hall width (m)"""


turbine_hall_h: float = None
"""turbine hall height (m)"""


tw_storage_l: float = None
"""tritiated waste storage building length (m)"""


tw_storage_w: float = None
"""tritiated waste storage building width (m)"""


tw_storage_h: float = None
"""tritiated waste storage building height (m)"""


volrci: float = None
"""internal volume of reactor building (m3)"""


volnucb: float = None
"""sum of nuclear buildings volumes (m3)"""


warm_shop_l: float = None
"""warm shop length (m)"""


warm_shop_w: float = None
"""warm shop width (m)"""


warm_shop_h: float = None
"""warm shop height (m)"""


water_buildings_l: float = None
"""water, laundry & drainage buildings length (m)"""


water_buildings_w: float = None
"""water, laundry & drainage buildings width (m)"""


water_buildings_h: float = None
"""water, laundry & drainage buildings height (m)"""


wgt: float = None
"""reactor building crane capacity (kg) (calculated if 0 is input)"""


wgt2: float = None
"""hot cell crane capacity (kg) (calculated if 0 is input)"""


workshop_l: float = None
"""[cold] workshop buildings length (m)"""


workshop_w: float = None
"""[cold] workshop buildings width (m)"""


workshop_h: float = None
"""[cold] workshop buildings height (m)"""


wrbi: float = None
"""distance from centre of machine to building wall (m)"""


wsvol: float = None
"""volume of warm shop building (m3)"""


wsvfac: float = None
"""warm shop building volume multiplication factor"""


a_reactor_bldg: float = None
"""Floor area of reactor building in m^2"""


a_ee_ps_bldg: float = None
"""Floor area of electrical equipment and power supply building in m^2"""


a_aux_services_bldg: float = None
"""Floor area of auxiliary services building in m^2"""


a_hot_cell_bldg: float = None
"""Floor area of hot cell building in m^2"""


a_reactor_service_bldg: float = None
"""Floor area of reactor service building in m^2"""


a_service_water_bldg: float = None
"""Floor area of service water building in m^2"""


a_fuel_handling_bldg: float = None
"""Floor area of fuel handling and storage building in m^2"""


a_control_room_bldg: float = None
"""Floor area of controlroom building in m^2"""


a_ac_ps_bldg: float = None
"""Floor area of AC power supply building in m^2"""


a_admin_bldg: float = None
"""Floor area of admin building in m^2"""


a_site_service_bldg: float = None
"""Floor area of site service building in m^2"""


a_cryo_inert_gas_bldg: float = None
"""Floor area of cryogenics and inert gas storage building in m^2"""


a_security_bldg: float = None
"""Floor area of security building in m^2"""


def init_buildings_variables():
    global admv
    global admvol
    global aux_build_l
    global aux_build_w
    global aux_build_h
    global auxcool_l
    global auxcool_w
    global auxcool_h
    global bioshld_thk
    global chemlab_l
    global chemlab_w
    global chemlab_h
    global dz_tf_cryostat
    global clh2
    global control_buildings_l
    global control_buildings_w
    global control_buildings_h
    global conv
    global convol
    global crane_arm_h
    global crane_clrnc_h
    global crane_clrnc_v
    global cryomag_l
    global cryomag_w
    global cryomag_h
    global cryostore_l
    global cryostore_w
    global cryostore_h
    global cryostat_clrnc
    global cryvol
    global a_plant_floor_effective
    global elecdist_l
    global elecdist_w
    global elecdist_h
    global elecload_l
    global elecload_w
    global elecload_h
    global elecstore_l
    global elecstore_w
    global elecstore_h
    global elevol
    global esbldgm3
    global fc_building_l
    global fc_building_w
    global fndt
    global gas_buildings_l
    global gas_buildings_w
    global gas_buildings_h
    global ground_clrnc
    global hcd_building_l
    global hcd_building_w
    global hcd_building_h
    global hccl
    global hcwt
    global heat_sink_l
    global heat_sink_w
    global heat_sink_h
    global hot_sepdist
    global hotcell_h
    global hw_storage_l
    global hw_storage_w
    global hw_storage_h
    global i_bldgs_size
    global i_bldgs_v
    global ilw_smelter_l
    global ilw_smelter_w
    global ilw_smelter_h
    global ilw_storage_l
    global ilw_storage_w
    global ilw_storage_h
    global llw_storage_l
    global llw_storage_w
    global llw_storage_h
    global magnet_pulse_l
    global magnet_pulse_w
    global magnet_pulse_h
    global magnet_trains_l
    global magnet_trains_w
    global magnet_trains_h
    global maint_cont_l
    global maint_cont_w
    global maint_cont_h
    global mbvfac
    global nbi_sys_l
    global nbi_sys_w
    global pfbldgm3
    global pibv
    global qnty_sfty_fac
    global rbvfac
    global rbrt
    global rbvol
    global rbwt
    global reactor_clrnc
    global reactor_fndtn_thk
    global reactor_hall_l
    global reactor_hall_w
    global reactor_hall_h
    global reactor_roof_thk
    global reactor_wall_thk
    global rmbvol
    global robotics_l
    global robotics_w
    global robotics_h
    global row
    global rxcl
    global sec_buildings_l
    global sec_buildings_w
    global sec_buildings_h
    global shmf
    global shov
    global shovol
    global staff_buildings_area
    global staff_buildings_h
    global stcl
    global tfcbv
    global transp_clrnc
    global trcl
    global triv
    global turbine_hall_l
    global turbine_hall_w
    global turbine_hall_h
    global tw_storage_l
    global tw_storage_w
    global tw_storage_h
    global volrci
    global volnucb
    global warm_shop_l
    global warm_shop_w
    global warm_shop_h
    global water_buildings_l
    global water_buildings_w
    global water_buildings_h
    global wgt
    global wgt2
    global workshop_l
    global workshop_w
    global workshop_h
    global wrbi
    global wsvol
    global wsvfac
    global a_reactor_bldg
    global a_ee_ps_bldg
    global a_aux_services_bldg
    global a_hot_cell_bldg
    global a_reactor_service_bldg
    global a_service_water_bldg
    global a_fuel_handling_bldg
    global a_control_room_bldg
    global a_ac_ps_bldg
    global a_admin_bldg
    global a_site_service_bldg
    global a_cryo_inert_gas_bldg
    global a_security_bldg

    admv = 1.0e5
    admvol = 0.0
    aux_build_l = 60.0
    aux_build_w = 30.0
    aux_build_h = 5.0
    auxcool_l = 20.0
    auxcool_w = 20.0
    auxcool_h = 5.0
    bioshld_thk = 2.50
    chemlab_l = 50.0
    chemlab_w = 30.0
    chemlab_h = 6.0
    dz_tf_cryostat = 2.5
    clh2 = 15.0
    control_buildings_l = 80.0
    control_buildings_w = 60.0
    control_buildings_h = 6.0
    conv = 6.0e4
    convol = 0.0
    crane_arm_h = 10.0
    crane_clrnc_h = 4.0
    crane_clrnc_v = 3.0
    cryomag_l = 120.0
    cryomag_w = 90.0
    cryomag_h = 5.0
    cryostore_l = 160.0
    cryostore_w = 30.0
    cryostore_h = 20.0
    cryostat_clrnc = 2.5
    cryvol = 0.0
    a_plant_floor_effective = 0.0
    elecdist_l = 380.0
    elecdist_w = 350.0
    elecdist_h = 5.0
    elecload_l = 100.0
    elecload_w = 90.0
    elecload_h = 3.0
    elecstore_l = 100.0
    elecstore_w = 60.0
    elecstore_h = 12.0
    elevol = 0.0
    esbldgm3 = 1.0e3
    fc_building_l = 60.0
    fc_building_w = 60.0
    fndt = 2.0
    gas_buildings_l = 25.0
    gas_buildings_w = 15.0
    gas_buildings_h = 5.0
    ground_clrnc = 5.0
    hcd_building_l = 70.0
    hcd_building_w = 40.0
    hcd_building_h = 25.0
    hw_storage_l = 20.0
    hw_storage_w = 10.0
    hw_storage_h = 5.0
    heat_sink_l = 160.0
    heat_sink_w = 80.0
    heat_sink_h = 12.0
    hccl = 5.0
    hcwt = 1.5
    hot_sepdist = 2.0
    hotcell_h = 12.0
    i_bldgs_size = 0
    i_bldgs_v = 0
    ilw_smelter_l = 50.0
    ilw_smelter_w = 30.0
    ilw_smelter_h = 30.0
    ilw_storage_l = 120.0
    ilw_storage_w = 100.0
    ilw_storage_h = 8.0
    llw_storage_l = 45.0
    llw_storage_w = 20.0
    llw_storage_h = 5.0
    magnet_pulse_l = 105.0
    magnet_pulse_w = 40.0
    magnet_pulse_h = 5.0
    magnet_trains_l = 120.0
    magnet_trains_w = 90.0
    magnet_trains_h = 5.0
    maint_cont_l = 125.0
    maint_cont_w = 100.0
    maint_cont_h = 6.0
    mbvfac = 2.8
    nbi_sys_l = 225.0
    nbi_sys_w = 185.0
    pfbldgm3 = 2.0e4
    pibv = 2.0e4
    qnty_sfty_fac = 2.0
    rbvfac = 1.6
    rbrt = 1.0
    rbvol = 0.0
    rbwt = 2.0
    reactor_clrnc = 4.0
    reactor_fndtn_thk = 2.0
    reactor_hall_l = 0.0
    reactor_hall_w = 0.0
    reactor_hall_h = 0.0
    reactor_roof_thk = 1.0
    reactor_wall_thk = 2.0
    rmbvol = 0.0
    robotics_l = 50.0
    robotics_w = 30.0
    robotics_h = 30.0
    row = 4.0
    rxcl = 4.0
    sec_buildings_l = 30.0
    sec_buildings_w = 25.0
    sec_buildings_h = 6.0
    shmf = 0.5
    shov = 1.0e5
    shovol = 0.0
    staff_buildings_h = 5.0
    staff_buildings_area = 4.8e5
    stcl = 3.0
    tfcbv = 2.0e4
    transp_clrnc = 1.0
    trcl = 1.0
    triv = 4.0e4
    turbine_hall_l = 109.0
    turbine_hall_w = 62.0
    turbine_hall_h = 15.0
    tw_storage_l = 90.0
    tw_storage_w = 30.0
    tw_storage_h = 5.0
    volnucb = 0.0
    volrci = 0.0
    warm_shop_l = 100.0
    warm_shop_w = 50.0
    warm_shop_h = 10.0
    water_buildings_l = 110.0
    water_buildings_w = 10.0
    water_buildings_h = 5.0
    workshop_l = 150.0
    workshop_w = 125.0
    workshop_h = 10.0
    wgt = 5.0e5
    wgt2 = 1.0e5
    wrbi = 0.0
    wsvfac = 1.9
    wsvol = 0.0
    a_reactor_bldg = 8.32e3
    a_ee_ps_bldg = 2.133e4
    a_aux_services_bldg = 1.0e3
    a_hot_cell_bldg = 8.43e3
    a_reactor_service_bldg = 2.44e3
    a_service_water_bldg = 1.567e3
    a_fuel_handling_bldg = 1.67e3
    a_control_room_bldg = 2.88e3
    a_ac_ps_bldg = 6.423e3
    a_admin_bldg = 2.5674e4
    a_site_service_bldg = 8.3e3
    a_cryo_inert_gas_bldg = 1.838e4
    a_security_bldg = 4.552e3
