module buildings_variables
    !! author: J. Morris (UKAEA)
    !!
    !! Module containing global variables relating to the plant buildings
    !!
    !! GIFA = Gross Internal Floor Area
    !!
    !!### References
    !!
    !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    public

    real(dp) :: admv
    !! administration building volume (m3)

    real(dp) :: admvol
    !! volume of administration buildings (m3)

    real(dp) :: aux_build_l, aux_build_w, aux_build_h
    !! aux building supporting tokamak processes length, width, height (m)

    real(dp) :: auxcool_l, auxcool_w, auxcool_h
    !! Site-Wide Auxiliary Cooling Water facility length, width, height (m)

    real(dp) :: bioshld_thk
    !! Radial thickness of bio-shield around reactor (m)

    real(dp) :: chemlab_l, chemlab_w, chemlab_h
    !! Chemistry labs and treatment buldings length, width, height (m)

    real(dp) :: dz_tf_cryostat
    !! vertical clearance from TF coil to cryostat (m) (calculated for tokamaks)

    real(dp) :: clh2
    !! clearance beneath TF coil to foundation (including basement) (m)

    real(dp) :: control_buildings_l, control_buildings_w, control_buildings_h
    !! control building length, width, height (m)

    real(dp) :: conv
    !! control building volume (m3)

    real(dp) :: convol
    !! volume of control, protection and i&c building (m3)

    real(dp) :: crane_arm_h
    !! vertical dimension of crane arm, operating over reactor (m)

    real(dp) :: crane_clrnc_h
    !! horizontal clearance to building wall for crane operation (m)

    real(dp) :: crane_clrnc_v
    !! vertical clearance for crane operation (m)

    real(dp) :: cryomag_l, cryomag_w, cryomag_h
    !! Cryogenic Buildings for Magnet and Fuel Cycle length, width, height (m)

    real(dp) :: cryostore_l, cryostore_w, cryostore_h
    !! Magnet Cryo Storage Tanks length, width, height (m)

    real(dp) :: cryostat_clrnc
    !! vertical clearance from TF coil to cryostat (m)

    real(dp) :: cryvol
    !! volume of cryoplant building (m3)

    real(dp) :: a_plant_floor_effective
    !! effective total floor space (m2)

    real(dp) :: elecdist_l, elecdist_w, elecdist_h
    !! Transformers and electrical distribution facilities length, width, height (m)

    real(dp) :: elecload_l, elecload_w, elecload_h
    !! Electric (eesential and non-essential) load centres length, width, height (m)

    real(dp) :: elecstore_l, elecstore_w, elecstore_h
    !! Energy Storage facilities length, width, height (m)

    real(dp) :: elevol
    !! volume of electrical equipment building (m3)

    real(dp) :: esbldgm3
    !! volume of energy storage equipment building (m3) (not used if `i_pulsed_plant=0`)

    real(dp) :: fc_building_l, fc_building_w
    !! Fuel Cycle facilities length, width (m)

    real(dp) :: fndt
    !! foundation thickness (m)

    real(dp) :: gas_buildings_l, gas_buildings_w, gas_buildings_h
    !! air & gas supply (amalgamated) buildings length, width, height (m)

    real(dp) :: ground_clrnc
    !! clearance beneath TF coil (m)

    real(dp) :: hcd_building_l, hcd_building_w, hcd_building_h
    !! HCD building length, width, height (m)

    real(dp) :: hccl
    !! clearance around components in hot cell (m)

    real(dp) :: hcwt
    !! hot cell wall thickness (m)

    real(dp) :: heat_sink_l, heat_sink_w, heat_sink_h
    !! heat sinks length, width, height (m)

    real(dp) :: hot_sepdist
    !! hot cell storage component separation distance (m)

    real(dp) :: hotcell_h
    !! hot cell storage and maintenance facility height (m)

    real(dp) :: hw_storage_l, hw_storage_w, hw_storage_h
    !! hazardous waste storage building length, width, height (m)

    integer :: i_bldgs_size
    !! switch between routines estimating building sizes (0 = default; 1 = updated)

    integer :: i_bldgs_v
    !! switch to select verbose output for buildings (1 = verbose)

    real(dp) :: ilw_smelter_l, ilw_smelter_w, ilw_smelter_h
    !! radioactive waste smelting facility length, width, height (m)

    real(dp) :: ilw_storage_l, ilw_storage_w, ilw_storage_h
    !! ILW waste storage building length, width, height (m)

    real(dp) :: llw_storage_l, llw_storage_w, llw_storage_h
    !! LLW waste storage building length, width, height (m)

    real(dp) :: magnet_pulse_l, magnet_pulse_w, magnet_pulse_h
    !! pulsed magnet power building length, width, height (m)

    real(dp) :: magnet_trains_l, magnet_trains_w, magnet_trains_h
    !! steady state magnet power trains building length, width, height (m)

    real(dp) :: maint_cont_l, maint_cont_w, maint_cont_h
    !! maintenance control building length, width, height (m)

    real(dp) :: mbvfac
    !! maintenance building volume multiplication factor

    real(dp) :: nbi_sys_l, nbi_sys_w
    !! NBI system length, width (m)

    real(dp) :: pfbldgm3
    !! volume of PF coil power supply building (m3)

    real(dp) :: pibv
    !! power injection building volume (m3)

    real(dp) :: qnty_sfty_fac
    !! quantity safety factor for component use during plant lifetime

    real(dp) :: rbvfac
    !! reactor building volume multiplication factor

    real(dp) :: rbrt
    !! reactor building roof thickness (m)

    real(dp) :: rbvol
    !! reactor building volume (m3)

    real(dp) :: rbwt
    !! reactor building wall thickness (m)

    real(dp) :: reactor_clrnc
    !! clearance around reactor (m)

    real(dp) :: reactor_fndtn_thk
    !! reactor building foundation thickness (m)

    real(dp) :: reactor_hall_l, reactor_hall_w, reactor_hall_h
    !! reactor building length, width, height (m)

    real(dp) :: reactor_roof_thk
    !! reactor building roof thickness (m)

    real(dp) :: reactor_wall_thk
    !! reactor building wall thickness (m)

    real(dp) :: rmbvol
    !! volume of maintenance and assembly building (m3)

    real(dp) :: robotics_l, robotics_w, robotics_h
    !! robotics buildings length, width, height (m)

    real(dp) :: row
    !! clearance to building wall for crane operation (m)

    real(dp) :: rxcl
    !! clearance around reactor (m)

    real(dp) :: sec_buildings_l, sec_buildings_w, sec_buildings_h
    !! security & safety buildings length, width, height (m)

    real(dp) :: shmf
    !! fraction of shield mass per TF coil to be moved in the maximum shield lift

    real(dp) :: shov
    !! shops and warehouse volume (m3)

    real(dp) :: shovol
    !! volume of shops and buildings for plant auxiliaries (m3)

    real(dp) :: staff_buildings_area
    !! footprint of staff buildings (m2)

    real(dp) :: staff_buildings_h
    !! staff buildings height (m)

    real(dp) :: stcl
    !! clearance above crane to roof (m)

    real(dp) :: tfcbv
    !! volume of TF coil power supply building (m3) (calculated if TF coils are superconducting)

    real(dp) :: transp_clrnc
    !! transportation clearance between components (m)

    real(dp) :: trcl
    !! transportation clearance between components (m)

    real(dp) :: triv
    !! volume of tritium, fuel handling and health physics buildings (m3)

    real(dp) :: turbine_hall_l, turbine_hall_w, turbine_hall_h
    !! turbine hall length, width, height (m)

    real(dp) :: tw_storage_l, tw_storage_w, tw_storage_h
    !! tritiated waste storage building length, width, height (m)

    real(dp) :: volrci
    !! internal volume of reactor building (m3)

    real(dp) :: volnucb
    !! sum of nuclear buildings volumes (m3)

    real(dp) :: warm_shop_l, warm_shop_w, warm_shop_h
    !! warm shop length, width, height (m)

    real(dp) :: water_buildings_l, water_buildings_w, water_buildings_h
    !! water, laundry & drainage buildings length, width, height (m)

    real(dp) :: wgt
    !! reactor building crane capacity (kg) (calculated if 0 is input)

    real(dp) :: wgt2
    !! hot cell crane capacity (kg) (calculated if 0 is input)

    real(dp) :: workshop_l, workshop_w, workshop_h
    !! [cold] workshop buildings length, width, height (m)

    real(dp) :: wrbi
    !! distance from centre of machine to building wall (m)

    real(dp) :: wsvol
    !! volume of warm shop building (m3)

    real(dp) :: wsvfac
    !! warm shop building volume multiplication factor



    real(dp) :: a_reactor_bldg
    !! Floor area of reactor building in m^2

    real(dp) :: a_ee_ps_bldg
    !! Floor area of electrical equipment and power supply building in m^2

    real(dp) :: a_aux_services_bldg
    !! Floor area of auxiliary services building in m^2

    real(dp) :: a_hot_cell_bldg
    !! Floor area of hot cell building in m^2

    real(dp) :: a_reactor_service_bldg
    !! Floor area of reactor service building in m^2

    real(dp) :: a_service_water_bldg
    !! Floor area of service water building in m^2

    real(dp) :: a_fuel_handling_bldg
    !! Floor area of fuel handling and storage building in m^2

    real(dp) :: a_control_room_bldg
    !! Floor area of controlroom building in m^2

    real(dp) :: a_ac_ps_bldg
    !! Floor area of AC power supply building in m^2

    real(dp) :: a_admin_bldg
    !! Floor area of admin building in m^2

    real(dp) :: a_site_service_bldg
    !! Floor area of site service building in m^2

    real(dp) :: a_cryo_inert_gas_bldg
    !! Floor area of cryogenics and inert gas storage building in m^2

    real(dp) :: a_security_bldg
    !! Floor area of security building in m^2
  end module buildings_variables
