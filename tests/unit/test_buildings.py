from typing import Any, NamedTuple

import pytest

from process.buildings import Buildings
from process.data_structure import (
    build_variables,
    buildings_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    pfcoil_variables,
    physics_variables,
    tfcoil_variables,
)


@pytest.fixture
def buildings():
    """Provides Buildings object for testing.

    :returns buildings: initialised Buildings object
    :rtype: process.buildings.Buildings
    """
    return Buildings()


class BldgsSizesParam(NamedTuple):
    i_bldgs_v: Any
    a_plant_floor_effective: Any
    volnucb: Any
    bioshld_thk: Any
    reactor_wall_thk: Any
    reactor_roof_thk: Any
    reactor_fndtn_thk: Any
    reactor_clrnc: Any
    transp_clrnc: Any
    cryostat_clrnc: Any
    ground_clrnc: Any
    crane_clrnc_h: Any
    crane_arm_h: Any
    reactor_hall_l: Any
    reactor_hall_w: Any
    reactor_hall_h: Any
    nbi_sys_l: Any
    nbi_sys_w: Any
    fc_building_l: Any
    fc_building_w: Any
    warm_shop_l: Any
    warm_shop_w: Any
    warm_shop_h: Any
    workshop_l: Any
    workshop_w: Any
    workshop_h: Any
    robotics_l: Any
    robotics_w: Any
    robotics_h: Any
    maint_cont_l: Any
    maint_cont_w: Any
    maint_cont_h: Any
    turbine_hall_l: Any
    turbine_hall_w: Any
    turbine_hall_h: Any
    gas_buildings_l: Any
    gas_buildings_w: Any
    gas_buildings_h: Any
    water_buildings_l: Any
    water_buildings_w: Any
    water_buildings_h: Any
    sec_buildings_l: Any
    sec_buildings_w: Any
    sec_buildings_h: Any
    staff_buildings_area: Any
    staff_buildings_h: Any
    hcd_building_l: Any
    hcd_building_w: Any
    hcd_building_h: Any
    magnet_pulse_l: Any
    magnet_pulse_w: Any
    magnet_pulse_h: Any
    magnet_trains_l: Any
    magnet_trains_w: Any
    magnet_trains_h: Any
    control_buildings_l: Any
    control_buildings_w: Any
    control_buildings_h: Any
    ilw_smelter_l: Any
    ilw_smelter_w: Any
    ilw_smelter_h: Any
    ilw_storage_l: Any
    ilw_storage_w: Any
    ilw_storage_h: Any
    llw_storage_l: Any
    llw_storage_w: Any
    llw_storage_h: Any
    hw_storage_l: Any
    hw_storage_w: Any
    hw_storage_h: Any
    tw_storage_l: Any
    tw_storage_w: Any
    tw_storage_h: Any
    auxcool_l: Any
    auxcool_w: Any
    auxcool_h: Any
    cryomag_l: Any
    cryomag_w: Any
    cryomag_h: Any
    cryostore_l: Any
    cryostore_w: Any
    cryostore_h: Any
    elecdist_l: Any
    elecdist_w: Any
    elecdist_h: Any
    elecstore_l: Any
    elecstore_w: Any
    elecstore_h: Any
    elecload_l: Any
    elecload_w: Any
    elecload_h: Any
    chemlab_l: Any
    chemlab_w: Any
    chemlab_h: Any
    heat_sink_l: Any
    heat_sink_w: Any
    heat_sink_h: Any
    aux_build_l: Any
    aux_build_w: Any
    aux_build_h: Any
    qnty_sfty_fac: Any
    hotcell_h: Any
    hot_sepdist: Any
    i_hcd_primary: Any
    n_tf_coils: Any
    i_tf_sup: Any
    r_pf_coil_outer_max: Any
    life_plant: Any
    cplife: Any
    life_div_fpy: Any
    r_cryostat_inboard: Any
    life_blkt_fpy: Any
    z_tf_inside_half: Any
    dr_tf_inboard: Any
    dr_tf_shld_gap: Any
    dr_shld_thermal_inboard: Any
    dr_shld_thermal_outboard: Any
    dz_shld_thermal: Any
    dr_shld_inboard: Any
    dr_shld_outboard: Any
    dr_fw_plasma_gap_inboard: Any
    dr_fw_plasma_gap_outboard: Any
    dr_fw_inboard: Any
    dr_fw_outboard: Any
    dr_blkt_inboard: Any
    dr_blkt_outboard: Any
    r_cp_top: Any
    dz_divertor: Any
    rmajor: Any
    rminor: Any
    tf_radial_dim: Any
    tf_vertical_dim: Any
    outfile: Any
    iprint: Any
    expected_reactor_hall_l: Any
    expected_reactor_hall_w: Any
    expected_reactor_hall_h: Any


@pytest.mark.parametrize(
    "bldgssizesparam",
    (
        BldgsSizesParam(
            i_bldgs_v=0,
            a_plant_floor_effective=0,
            volnucb=0,
            bioshld_thk=2.5,
            reactor_wall_thk=2,
            reactor_roof_thk=1,
            reactor_fndtn_thk=2,
            reactor_clrnc=4,
            transp_clrnc=1,
            cryostat_clrnc=2.5,
            ground_clrnc=5,
            crane_clrnc_h=4,
            crane_arm_h=10,
            reactor_hall_l=0,
            reactor_hall_w=0,
            reactor_hall_h=0,
            nbi_sys_l=225,
            nbi_sys_w=185,
            fc_building_l=60,
            fc_building_w=60,
            warm_shop_l=100,
            warm_shop_w=50,
            warm_shop_h=10,
            workshop_l=150,
            workshop_w=125,
            workshop_h=10,
            robotics_l=50,
            robotics_w=30,
            robotics_h=30,
            maint_cont_l=125,
            maint_cont_w=100,
            maint_cont_h=6,
            turbine_hall_l=109,
            turbine_hall_w=62,
            turbine_hall_h=15,
            gas_buildings_l=25,
            gas_buildings_w=15,
            gas_buildings_h=5,
            water_buildings_l=110,
            water_buildings_w=10,
            water_buildings_h=5,
            sec_buildings_l=30,
            sec_buildings_w=25,
            sec_buildings_h=6,
            staff_buildings_area=480000,
            staff_buildings_h=5,
            hcd_building_l=70,
            hcd_building_w=40,
            hcd_building_h=25,
            magnet_pulse_l=105,
            magnet_pulse_w=40,
            magnet_pulse_h=5,
            magnet_trains_l=120,
            magnet_trains_w=90,
            magnet_trains_h=5,
            control_buildings_l=80,
            control_buildings_w=60,
            control_buildings_h=6,
            ilw_smelter_l=50,
            ilw_smelter_w=30,
            ilw_smelter_h=30,
            ilw_storage_l=120,
            ilw_storage_w=100,
            ilw_storage_h=8,
            llw_storage_l=45,
            llw_storage_w=20,
            llw_storage_h=5,
            hw_storage_l=20,
            hw_storage_w=10,
            hw_storage_h=5,
            tw_storage_l=90,
            tw_storage_w=30,
            tw_storage_h=5,
            auxcool_l=20,
            auxcool_w=20,
            auxcool_h=5,
            cryomag_l=120,
            cryomag_w=90,
            cryomag_h=5,
            cryostore_l=160,
            cryostore_w=30,
            cryostore_h=20,
            elecdist_l=380,
            elecdist_w=350,
            elecdist_h=5,
            elecstore_l=100,
            elecstore_w=60,
            elecstore_h=12,
            elecload_l=100,
            elecload_w=90,
            elecload_h=3,
            chemlab_l=50,
            chemlab_w=30,
            chemlab_h=6,
            heat_sink_l=160,
            heat_sink_w=80,
            heat_sink_h=12,
            aux_build_l=60,
            aux_build_w=30,
            aux_build_h=5,
            qnty_sfty_fac=2,
            hotcell_h=12,
            hot_sepdist=2,
            i_hcd_primary=10,
            n_tf_coils=16,
            i_tf_sup=1,
            r_pf_coil_outer_max=18.98258241468535,
            life_plant=40,
            cplife=0,
            life_div_fpy=0,
            r_cryostat_inboard=19.48258241468535,
            life_blkt_fpy=0,
            z_tf_inside_half=9.0730900215620327,
            dr_tf_inboard=1.208,
            dr_tf_shld_gap=0.05000000000000001,
            dr_shld_thermal_inboard=0.050000000000000003,
            dr_shld_thermal_outboard=0.050000000000000003,
            dz_shld_thermal=0.050000000000000003,
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            dr_fw_plasma_gap_inboard=0.22500000000000003,
            dr_fw_plasma_gap_outboard=0.22500000000000003,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            dr_blkt_outboard=0.98199999999999998,
            r_cp_top=4.20194118510911,
            dz_divertor=0.62100000000000011,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            tf_radial_dim=14.129464674334221,
            tf_vertical_dim=20.562180043124066,
            outfile=11,
            iprint=0,
            expected_reactor_hall_l=218.89549448811209,
            expected_reactor_hall_w=218.89549448811209,
            expected_reactor_hall_h=67.624360086248132,
        ),
        BldgsSizesParam(
            i_bldgs_v=0,
            a_plant_floor_effective=1539392.0963074313,
            volnucb=5212998.1139194397,
            bioshld_thk=2.5,
            reactor_wall_thk=2,
            reactor_roof_thk=1,
            reactor_fndtn_thk=2,
            reactor_clrnc=4,
            transp_clrnc=1,
            cryostat_clrnc=2.5,
            ground_clrnc=5,
            crane_clrnc_h=4,
            crane_arm_h=10,
            reactor_hall_l=218.89549448811209,
            reactor_hall_w=218.89549448811209,
            reactor_hall_h=67.624360086248132,
            nbi_sys_l=225,
            nbi_sys_w=185,
            fc_building_l=60,
            fc_building_w=60,
            warm_shop_l=100,
            warm_shop_w=50,
            warm_shop_h=10,
            workshop_l=150,
            workshop_w=125,
            workshop_h=10,
            robotics_l=50,
            robotics_w=30,
            robotics_h=30,
            maint_cont_l=125,
            maint_cont_w=100,
            maint_cont_h=6,
            turbine_hall_l=109,
            turbine_hall_w=62,
            turbine_hall_h=15,
            gas_buildings_l=25,
            gas_buildings_w=15,
            gas_buildings_h=5,
            water_buildings_l=110,
            water_buildings_w=10,
            water_buildings_h=5,
            sec_buildings_l=30,
            sec_buildings_w=25,
            sec_buildings_h=6,
            staff_buildings_area=480000,
            staff_buildings_h=5,
            hcd_building_l=70,
            hcd_building_w=40,
            hcd_building_h=25,
            magnet_pulse_l=105,
            magnet_pulse_w=40,
            magnet_pulse_h=5,
            magnet_trains_l=120,
            magnet_trains_w=90,
            magnet_trains_h=5,
            control_buildings_l=80,
            control_buildings_w=60,
            control_buildings_h=6,
            ilw_smelter_l=50,
            ilw_smelter_w=30,
            ilw_smelter_h=30,
            ilw_storage_l=120,
            ilw_storage_w=100,
            ilw_storage_h=8,
            llw_storage_l=45,
            llw_storage_w=20,
            llw_storage_h=5,
            hw_storage_l=20,
            hw_storage_w=10,
            hw_storage_h=5,
            tw_storage_l=90,
            tw_storage_w=30,
            tw_storage_h=5,
            auxcool_l=20,
            auxcool_w=20,
            auxcool_h=5,
            cryomag_l=120,
            cryomag_w=90,
            cryomag_h=5,
            cryostore_l=160,
            cryostore_w=30,
            cryostore_h=20,
            elecdist_l=380,
            elecdist_w=350,
            elecdist_h=5,
            elecstore_l=100,
            elecstore_w=60,
            elecstore_h=12,
            elecload_l=100,
            elecload_w=90,
            elecload_h=3,
            chemlab_l=50,
            chemlab_w=30,
            chemlab_h=6,
            heat_sink_l=160,
            heat_sink_w=80,
            heat_sink_h=12,
            aux_build_l=60,
            aux_build_w=30,
            aux_build_h=5,
            qnty_sfty_fac=2,
            hotcell_h=12,
            hot_sepdist=2,
            i_hcd_primary=10,
            n_tf_coils=16,
            i_tf_sup=1,
            r_pf_coil_outer_max=18.982980877139834,
            life_plant=40,
            cplife=0,
            life_div_fpy=6.1337250397740126,
            r_cryostat_inboard=19.482980877139834,
            life_blkt_fpy=19.216116010620578,
            z_tf_inside_half=9.0730900215620327,
            dr_tf_inboard=1.208,
            dr_tf_shld_gap=0.05000000000000001,
            dr_shld_thermal_inboard=0.050000000000000003,
            dr_shld_thermal_outboard=0.050000000000000003,
            dz_shld_thermal=0.050000000000000003,
            dr_shld_inboard=0.30000000000000004,
            dr_shld_outboard=0.80000000000000004,
            dr_fw_plasma_gap_inboard=0.22500000000000003,
            dr_fw_plasma_gap_outboard=0.22500000000000003,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_outboard=0.018000000000000002,
            dr_blkt_inboard=0.75500000000000012,
            dr_blkt_outboard=0.98199999999999998,
            r_cp_top=4.20194118510911,
            dz_divertor=0.62100000000000011,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            tf_radial_dim=14.129464674334221,
            tf_vertical_dim=20.562180043124066,
            outfile=11,
            iprint=0,
            expected_reactor_hall_l=218.897885262839,
            expected_reactor_hall_w=218.897885262839,
            expected_reactor_hall_h=67.624360086248132,
        ),
    ),
)
def test_bldgs_sizes(buildings, bldgssizesparam, monkeypatch):
    monkeypatch.setattr(buildings_variables, "i_bldgs_v", bldgssizesparam.i_bldgs_v)
    monkeypatch.setattr(
        buildings_variables,
        "a_plant_floor_effective",
        bldgssizesparam.a_plant_floor_effective,
    )
    monkeypatch.setattr(buildings_variables, "volnucb", bldgssizesparam.volnucb)
    monkeypatch.setattr(buildings_variables, "bioshld_thk", bldgssizesparam.bioshld_thk)
    monkeypatch.setattr(
        buildings_variables, "reactor_wall_thk", bldgssizesparam.reactor_wall_thk
    )
    monkeypatch.setattr(
        buildings_variables, "reactor_roof_thk", bldgssizesparam.reactor_roof_thk
    )
    monkeypatch.setattr(
        buildings_variables, "reactor_fndtn_thk", bldgssizesparam.reactor_fndtn_thk
    )
    monkeypatch.setattr(
        buildings_variables, "reactor_clrnc", bldgssizesparam.reactor_clrnc
    )
    monkeypatch.setattr(
        buildings_variables, "transp_clrnc", bldgssizesparam.transp_clrnc
    )
    monkeypatch.setattr(
        buildings_variables, "cryostat_clrnc", bldgssizesparam.cryostat_clrnc
    )
    monkeypatch.setattr(
        buildings_variables, "ground_clrnc", bldgssizesparam.ground_clrnc
    )
    monkeypatch.setattr(
        buildings_variables, "crane_clrnc_h", bldgssizesparam.crane_clrnc_h
    )
    monkeypatch.setattr(buildings_variables, "crane_arm_h", bldgssizesparam.crane_arm_h)
    monkeypatch.setattr(
        buildings_variables, "reactor_hall_l", bldgssizesparam.reactor_hall_l
    )
    monkeypatch.setattr(
        buildings_variables, "reactor_hall_w", bldgssizesparam.reactor_hall_w
    )
    monkeypatch.setattr(
        buildings_variables, "reactor_hall_h", bldgssizesparam.reactor_hall_h
    )
    monkeypatch.setattr(buildings_variables, "nbi_sys_l", bldgssizesparam.nbi_sys_l)
    monkeypatch.setattr(buildings_variables, "nbi_sys_w", bldgssizesparam.nbi_sys_w)
    monkeypatch.setattr(
        buildings_variables, "fc_building_l", bldgssizesparam.fc_building_l
    )
    monkeypatch.setattr(
        buildings_variables, "fc_building_w", bldgssizesparam.fc_building_w
    )
    monkeypatch.setattr(buildings_variables, "warm_shop_l", bldgssizesparam.warm_shop_l)
    monkeypatch.setattr(buildings_variables, "warm_shop_w", bldgssizesparam.warm_shop_w)
    monkeypatch.setattr(buildings_variables, "warm_shop_h", bldgssizesparam.warm_shop_h)
    monkeypatch.setattr(buildings_variables, "workshop_l", bldgssizesparam.workshop_l)
    monkeypatch.setattr(buildings_variables, "workshop_w", bldgssizesparam.workshop_w)
    monkeypatch.setattr(buildings_variables, "workshop_h", bldgssizesparam.workshop_h)
    monkeypatch.setattr(buildings_variables, "robotics_l", bldgssizesparam.robotics_l)
    monkeypatch.setattr(buildings_variables, "robotics_w", bldgssizesparam.robotics_w)
    monkeypatch.setattr(buildings_variables, "robotics_h", bldgssizesparam.robotics_h)
    monkeypatch.setattr(
        buildings_variables, "maint_cont_l", bldgssizesparam.maint_cont_l
    )
    monkeypatch.setattr(
        buildings_variables, "maint_cont_w", bldgssizesparam.maint_cont_w
    )
    monkeypatch.setattr(
        buildings_variables, "maint_cont_h", bldgssizesparam.maint_cont_h
    )
    monkeypatch.setattr(
        buildings_variables, "turbine_hall_l", bldgssizesparam.turbine_hall_l
    )
    monkeypatch.setattr(
        buildings_variables, "turbine_hall_w", bldgssizesparam.turbine_hall_w
    )
    monkeypatch.setattr(
        buildings_variables, "turbine_hall_h", bldgssizesparam.turbine_hall_h
    )
    monkeypatch.setattr(
        buildings_variables, "gas_buildings_l", bldgssizesparam.gas_buildings_l
    )
    monkeypatch.setattr(
        buildings_variables, "gas_buildings_w", bldgssizesparam.gas_buildings_w
    )
    monkeypatch.setattr(
        buildings_variables, "gas_buildings_h", bldgssizesparam.gas_buildings_h
    )
    monkeypatch.setattr(
        buildings_variables, "water_buildings_l", bldgssizesparam.water_buildings_l
    )
    monkeypatch.setattr(
        buildings_variables, "water_buildings_w", bldgssizesparam.water_buildings_w
    )
    monkeypatch.setattr(
        buildings_variables, "water_buildings_h", bldgssizesparam.water_buildings_h
    )
    monkeypatch.setattr(
        buildings_variables, "sec_buildings_l", bldgssizesparam.sec_buildings_l
    )
    monkeypatch.setattr(
        buildings_variables, "sec_buildings_w", bldgssizesparam.sec_buildings_w
    )
    monkeypatch.setattr(
        buildings_variables, "sec_buildings_h", bldgssizesparam.sec_buildings_h
    )
    monkeypatch.setattr(
        buildings_variables,
        "staff_buildings_area",
        bldgssizesparam.staff_buildings_area,
    )
    monkeypatch.setattr(
        buildings_variables, "staff_buildings_h", bldgssizesparam.staff_buildings_h
    )
    monkeypatch.setattr(
        buildings_variables, "hcd_building_l", bldgssizesparam.hcd_building_l
    )
    monkeypatch.setattr(
        buildings_variables, "hcd_building_w", bldgssizesparam.hcd_building_w
    )
    monkeypatch.setattr(
        buildings_variables, "hcd_building_h", bldgssizesparam.hcd_building_h
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_pulse_l", bldgssizesparam.magnet_pulse_l
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_pulse_w", bldgssizesparam.magnet_pulse_w
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_pulse_h", bldgssizesparam.magnet_pulse_h
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_trains_l", bldgssizesparam.magnet_trains_l
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_trains_w", bldgssizesparam.magnet_trains_w
    )
    monkeypatch.setattr(
        buildings_variables, "magnet_trains_h", bldgssizesparam.magnet_trains_h
    )
    monkeypatch.setattr(
        buildings_variables, "control_buildings_l", bldgssizesparam.control_buildings_l
    )
    monkeypatch.setattr(
        buildings_variables, "control_buildings_w", bldgssizesparam.control_buildings_w
    )
    monkeypatch.setattr(
        buildings_variables, "control_buildings_h", bldgssizesparam.control_buildings_h
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_smelter_l", bldgssizesparam.ilw_smelter_l
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_smelter_w", bldgssizesparam.ilw_smelter_w
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_smelter_h", bldgssizesparam.ilw_smelter_h
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_storage_l", bldgssizesparam.ilw_storage_l
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_storage_w", bldgssizesparam.ilw_storage_w
    )
    monkeypatch.setattr(
        buildings_variables, "ilw_storage_h", bldgssizesparam.ilw_storage_h
    )
    monkeypatch.setattr(
        buildings_variables, "llw_storage_l", bldgssizesparam.llw_storage_l
    )
    monkeypatch.setattr(
        buildings_variables, "llw_storage_w", bldgssizesparam.llw_storage_w
    )
    monkeypatch.setattr(
        buildings_variables, "llw_storage_h", bldgssizesparam.llw_storage_h
    )
    monkeypatch.setattr(
        buildings_variables, "hw_storage_l", bldgssizesparam.hw_storage_l
    )
    monkeypatch.setattr(
        buildings_variables, "hw_storage_w", bldgssizesparam.hw_storage_w
    )
    monkeypatch.setattr(
        buildings_variables, "hw_storage_h", bldgssizesparam.hw_storage_h
    )
    monkeypatch.setattr(
        buildings_variables, "tw_storage_l", bldgssizesparam.tw_storage_l
    )
    monkeypatch.setattr(
        buildings_variables, "tw_storage_w", bldgssizesparam.tw_storage_w
    )
    monkeypatch.setattr(
        buildings_variables, "tw_storage_h", bldgssizesparam.tw_storage_h
    )
    monkeypatch.setattr(buildings_variables, "auxcool_l", bldgssizesparam.auxcool_l)
    monkeypatch.setattr(buildings_variables, "auxcool_w", bldgssizesparam.auxcool_w)
    monkeypatch.setattr(buildings_variables, "auxcool_h", bldgssizesparam.auxcool_h)
    monkeypatch.setattr(buildings_variables, "cryomag_l", bldgssizesparam.cryomag_l)
    monkeypatch.setattr(buildings_variables, "cryomag_w", bldgssizesparam.cryomag_w)
    monkeypatch.setattr(buildings_variables, "cryomag_h", bldgssizesparam.cryomag_h)
    monkeypatch.setattr(buildings_variables, "cryostore_l", bldgssizesparam.cryostore_l)
    monkeypatch.setattr(buildings_variables, "cryostore_w", bldgssizesparam.cryostore_w)
    monkeypatch.setattr(buildings_variables, "cryostore_h", bldgssizesparam.cryostore_h)
    monkeypatch.setattr(buildings_variables, "elecdist_l", bldgssizesparam.elecdist_l)
    monkeypatch.setattr(buildings_variables, "elecdist_w", bldgssizesparam.elecdist_w)
    monkeypatch.setattr(buildings_variables, "elecdist_h", bldgssizesparam.elecdist_h)
    monkeypatch.setattr(buildings_variables, "elecstore_l", bldgssizesparam.elecstore_l)
    monkeypatch.setattr(buildings_variables, "elecstore_w", bldgssizesparam.elecstore_w)
    monkeypatch.setattr(buildings_variables, "elecstore_h", bldgssizesparam.elecstore_h)
    monkeypatch.setattr(buildings_variables, "elecload_l", bldgssizesparam.elecload_l)
    monkeypatch.setattr(buildings_variables, "elecload_w", bldgssizesparam.elecload_w)
    monkeypatch.setattr(buildings_variables, "elecload_h", bldgssizesparam.elecload_h)
    monkeypatch.setattr(buildings_variables, "chemlab_l", bldgssizesparam.chemlab_l)
    monkeypatch.setattr(buildings_variables, "chemlab_w", bldgssizesparam.chemlab_w)
    monkeypatch.setattr(buildings_variables, "chemlab_h", bldgssizesparam.chemlab_h)
    monkeypatch.setattr(buildings_variables, "heat_sink_l", bldgssizesparam.heat_sink_l)
    monkeypatch.setattr(buildings_variables, "heat_sink_w", bldgssizesparam.heat_sink_w)
    monkeypatch.setattr(buildings_variables, "heat_sink_h", bldgssizesparam.heat_sink_h)
    monkeypatch.setattr(buildings_variables, "aux_build_l", bldgssizesparam.aux_build_l)
    monkeypatch.setattr(buildings_variables, "aux_build_w", bldgssizesparam.aux_build_w)
    monkeypatch.setattr(buildings_variables, "aux_build_h", bldgssizesparam.aux_build_h)
    monkeypatch.setattr(
        buildings_variables, "qnty_sfty_fac", bldgssizesparam.qnty_sfty_fac
    )
    monkeypatch.setattr(buildings_variables, "hotcell_h", bldgssizesparam.hotcell_h)
    monkeypatch.setattr(buildings_variables, "hot_sepdist", bldgssizesparam.hot_sepdist)
    monkeypatch.setattr(
        current_drive_variables, "i_hcd_primary", bldgssizesparam.i_hcd_primary
    )
    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", bldgssizesparam.n_tf_coils)
    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", bldgssizesparam.i_tf_sup)
    monkeypatch.setattr(
        pfcoil_variables, "r_pf_coil_outer_max", bldgssizesparam.r_pf_coil_outer_max
    )
    monkeypatch.setattr(cost_variables, "life_plant", bldgssizesparam.life_plant)
    monkeypatch.setattr(cost_variables, "cplife", bldgssizesparam.cplife)
    monkeypatch.setattr(cost_variables, "life_div_fpy", bldgssizesparam.life_div_fpy)
    monkeypatch.setattr(
        fwbs_variables, "r_cryostat_inboard", bldgssizesparam.r_cryostat_inboard
    )
    monkeypatch.setattr(fwbs_variables, "life_blkt_fpy", bldgssizesparam.life_blkt_fpy)
    monkeypatch.setattr(
        build_variables, "z_tf_inside_half", bldgssizesparam.z_tf_inside_half
    )
    monkeypatch.setattr(build_variables, "dr_tf_inboard", bldgssizesparam.dr_tf_inboard)
    monkeypatch.setattr(
        build_variables, "dr_tf_shld_gap", bldgssizesparam.dr_tf_shld_gap
    )
    monkeypatch.setattr(
        build_variables,
        "dr_shld_thermal_inboard",
        bldgssizesparam.dr_shld_thermal_inboard,
    )
    monkeypatch.setattr(
        build_variables,
        "dr_shld_thermal_outboard",
        bldgssizesparam.dr_shld_thermal_outboard,
    )
    monkeypatch.setattr(
        build_variables, "dz_shld_thermal", bldgssizesparam.dz_shld_thermal
    )
    monkeypatch.setattr(
        build_variables, "dr_shld_inboard", bldgssizesparam.dr_shld_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_shld_outboard", bldgssizesparam.dr_shld_outboard
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_inboard",
        bldgssizesparam.dr_fw_plasma_gap_inboard,
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_outboard",
        bldgssizesparam.dr_fw_plasma_gap_outboard,
    )
    monkeypatch.setattr(build_variables, "dr_fw_inboard", bldgssizesparam.dr_fw_inboard)
    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", bldgssizesparam.dr_fw_outboard
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", bldgssizesparam.dr_blkt_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", bldgssizesparam.dr_blkt_outboard
    )
    monkeypatch.setattr(build_variables, "r_cp_top", bldgssizesparam.r_cp_top)
    monkeypatch.setattr(divertor_variables, "dz_divertor", bldgssizesparam.dz_divertor)
    monkeypatch.setattr(physics_variables, "rmajor", bldgssizesparam.rmajor)
    monkeypatch.setattr(physics_variables, "rminor", bldgssizesparam.rminor)

    buildings.bldgs_sizes(
        tf_radial_dim=bldgssizesparam.tf_radial_dim,
        tf_vertical_dim=bldgssizesparam.tf_vertical_dim,
        output=False,
    )

    assert buildings_variables.reactor_hall_l == pytest.approx(
        bldgssizesparam.expected_reactor_hall_l
    )
    assert buildings_variables.reactor_hall_w == pytest.approx(
        bldgssizesparam.expected_reactor_hall_w
    )
    assert buildings_variables.reactor_hall_h == pytest.approx(
        bldgssizesparam.expected_reactor_hall_h
    )


class BldgsParam(NamedTuple):
    wrbi: Any
    rxcl: Any
    trcl: Any
    row: Any
    wgt: Any
    shmf: Any
    clh2: Any
    dz_tf_cryostat: Any
    stcl: Any
    rbvfac: Any
    rbwt: Any
    rbrt: Any
    fndt: Any
    hcwt: Any
    hccl: Any
    wgt2: Any
    mbvfac: Any
    wsvfac: Any
    tfcbv: Any
    pfbldgm3: Any
    esbldgm3: Any
    pibv: Any
    a_plant_floor_effective: Any
    admvol: Any
    triv: Any
    conv: Any
    admv: Any
    shov: Any
    shovol: Any
    convol: Any
    volnucb: Any
    iprint: Any
    outfile: Any
    pfr: Any
    pfm: Any
    tfro: Any
    tfri: Any
    tfh: Any
    tfm: Any
    n_tf_coils: Any
    shro: Any
    shri: Any
    shh: Any
    shm: Any
    crr: Any
    helpow: Any
    expected_wrbi: Any
    expected_a_plant_floor_effective: Any
    expected_admvol: Any
    expected_shovol: Any
    expected_convol: Any
    expected_volnucb: Any
    expected_cryv: Any
    expected_vrci: Any
    expected_rbv: Any
    expected_rmbv: Any
    expected_wsv: Any
    expected_elev: Any


@pytest.mark.parametrize(
    "bldgsparam",
    (
        BldgsParam(
            wrbi=0,
            rxcl=4,
            trcl=1,
            row=4,
            wgt=500000,
            shmf=0.5,
            clh2=15,
            dz_tf_cryostat=5.7514039424138126,
            stcl=3,
            rbvfac=1.6000000000000001,
            rbwt=2,
            rbrt=1,
            fndt=2,
            hcwt=1.5,
            hccl=5,
            wgt2=100000,
            mbvfac=2.7999999999999998,
            wsvfac=1.8999999999999999,
            tfcbv=10601.097615432001,
            pfbldgm3=20000,
            esbldgm3=1000,
            pibv=20000,
            a_plant_floor_effective=0,
            admvol=0,
            triv=40000,
            conv=60000,
            admv=100000,
            shov=100000,
            shovol=0,
            convol=0,
            volnucb=0,
            iprint=0,
            outfile=11,
            pfr=18.98258241468535,
            pfm=1071.5897090529959,
            tfro=17.123405859443331,
            tfri=2.9939411851091102,
            tfh=20.562180043124066,
            tfm=1327.1818597762153,
            n_tf_coils=16,
            shro=13.764874193548387,
            shri=4.7423258064516141,
            shh=17.446180043124063,
            shm=2294873.8131476026,
            crr=19.48258241468535,
            helpow=77840.021662652987,
            expected_wrbi=42.612047089019569,
            expected_a_plant_floor_effective=379235.17804514873,
            expected_admvol=100000,
            expected_shovol=100000,
            expected_convol=60000,
            expected_volnucb=1812276.5359386117,
            expected_cryv=15344.903568596488,
            expected_vrci=1205439.8543893537,
            expected_rbv=1356973.2891062023,
            expected_rmbv=421473.52130148414,
            expected_wsv=130018.25667917728,
            expected_elev=51601.097615432001,
        ),
        BldgsParam(
            wrbi=42.612047089019569,
            rxcl=4,
            trcl=1,
            row=4,
            wgt=500000,
            shmf=0.5,
            clh2=15,
            dz_tf_cryostat=5.8405005070918357,
            stcl=3,
            rbvfac=1.6000000000000001,
            rbwt=2,
            rbrt=1,
            fndt=2,
            hcwt=1.5,
            hccl=5,
            wgt2=100000,
            mbvfac=2.7999999999999998,
            wsvfac=1.8999999999999999,
            tfcbv=10609.268177478583,
            pfbldgm3=20000,
            esbldgm3=1000,
            pibv=20000,
            a_plant_floor_effective=379235.17804514873,
            admvol=100000,
            triv=40000,
            conv=60000,
            admv=100000,
            shov=100000,
            shovol=100000,
            convol=60000,
            volnucb=1812276.5359386117,
            iprint=0,
            outfile=11,
            pfr=18.982980877139834,
            pfm=1073.3372194668184,
            tfro=17.123405859443331,
            tfri=2.9939411851091102,
            tfh=20.562180043124066,
            tfm=1327.9750836697808,
            n_tf_coils=16,
            shro=13.782874193548388,
            shri=4.7243258064516143,
            shh=17.446180043124063,
            shm=2297808.3935174868,
            crr=19.482980877139834,
            helpow=221493.99746816326,
            expected_wrbi=42.612445551474053,
            expected_a_plant_floor_effective=381590.59475257091,
            expected_admvol=100000,
            expected_shovol=100000,
            expected_convol=60000,
            expected_volnucb=1826281.0182016799,
            expected_cryv=25884.731838309508,
            expected_vrci=1206887.4047542624,
            expected_rbv=1358540.6868905292,
            expected_rmbv=423252.94369581528,
            expected_wsv=130255.93791329287,
            expected_elev=51609.268177478581,
        ),
    ),
)
def test_bldgs(buildings, bldgsparam, monkeypatch):
    monkeypatch.setattr(buildings_variables, "wrbi", bldgsparam.wrbi)
    monkeypatch.setattr(buildings_variables, "rxcl", bldgsparam.rxcl)
    monkeypatch.setattr(buildings_variables, "trcl", bldgsparam.trcl)
    monkeypatch.setattr(buildings_variables, "row", bldgsparam.row)
    monkeypatch.setattr(buildings_variables, "wgt", bldgsparam.wgt)
    monkeypatch.setattr(buildings_variables, "shmf", bldgsparam.shmf)
    monkeypatch.setattr(buildings_variables, "clh2", bldgsparam.clh2)
    monkeypatch.setattr(buildings_variables, "dz_tf_cryostat", bldgsparam.dz_tf_cryostat)
    monkeypatch.setattr(buildings_variables, "stcl", bldgsparam.stcl)
    monkeypatch.setattr(buildings_variables, "rbvfac", bldgsparam.rbvfac)
    monkeypatch.setattr(buildings_variables, "rbwt", bldgsparam.rbwt)
    monkeypatch.setattr(buildings_variables, "rbrt", bldgsparam.rbrt)
    monkeypatch.setattr(buildings_variables, "fndt", bldgsparam.fndt)
    monkeypatch.setattr(buildings_variables, "hcwt", bldgsparam.hcwt)
    monkeypatch.setattr(buildings_variables, "hccl", bldgsparam.hccl)
    monkeypatch.setattr(buildings_variables, "wgt2", bldgsparam.wgt2)
    monkeypatch.setattr(buildings_variables, "mbvfac", bldgsparam.mbvfac)
    monkeypatch.setattr(buildings_variables, "wsvfac", bldgsparam.wsvfac)
    monkeypatch.setattr(buildings_variables, "tfcbv", bldgsparam.tfcbv)
    monkeypatch.setattr(buildings_variables, "pfbldgm3", bldgsparam.pfbldgm3)
    monkeypatch.setattr(buildings_variables, "esbldgm3", bldgsparam.esbldgm3)
    monkeypatch.setattr(buildings_variables, "pibv", bldgsparam.pibv)
    monkeypatch.setattr(
        buildings_variables,
        "a_plant_floor_effective",
        bldgsparam.a_plant_floor_effective,
    )
    monkeypatch.setattr(buildings_variables, "admvol", bldgsparam.admvol)
    monkeypatch.setattr(buildings_variables, "triv", bldgsparam.triv)
    monkeypatch.setattr(buildings_variables, "conv", bldgsparam.conv)
    monkeypatch.setattr(buildings_variables, "admv", bldgsparam.admv)
    monkeypatch.setattr(buildings_variables, "shov", bldgsparam.shov)
    monkeypatch.setattr(buildings_variables, "shovol", bldgsparam.shovol)
    monkeypatch.setattr(buildings_variables, "convol", bldgsparam.convol)
    monkeypatch.setattr(buildings_variables, "volnucb", bldgsparam.volnucb)

    cryv, vrci, rbv, rmbv, wsv, elev = buildings.bldgs(
        output=False,
        pfr=bldgsparam.pfr,
        pfm=bldgsparam.pfm,
        tfro=bldgsparam.tfro,
        tfri=bldgsparam.tfri,
        tfh=bldgsparam.tfh,
        tfm=bldgsparam.tfm,
        n_tf_coils=bldgsparam.n_tf_coils,
        shro=bldgsparam.shro,
        shri=bldgsparam.shri,
        shh=bldgsparam.shh,
        shm=bldgsparam.shm,
        crr=bldgsparam.crr,
        helpow=bldgsparam.helpow,
    )

    assert buildings_variables.wrbi == pytest.approx(bldgsparam.expected_wrbi)
    assert buildings_variables.a_plant_floor_effective == pytest.approx(
        bldgsparam.expected_a_plant_floor_effective
    )
    assert buildings_variables.admvol == pytest.approx(bldgsparam.expected_admvol)
    assert buildings_variables.shovol == pytest.approx(bldgsparam.expected_shovol)
    assert buildings_variables.convol == pytest.approx(bldgsparam.expected_convol)
    assert buildings_variables.volnucb == pytest.approx(bldgsparam.expected_volnucb)

    assert cryv == pytest.approx(bldgsparam.expected_cryv)
    assert vrci == pytest.approx(bldgsparam.expected_vrci)
    assert rbv == pytest.approx(bldgsparam.expected_rbv)
    assert rmbv == pytest.approx(bldgsparam.expected_rmbv)
    assert wsv == pytest.approx(bldgsparam.expected_wsv)
    assert elev == pytest.approx(bldgsparam.expected_elev)
