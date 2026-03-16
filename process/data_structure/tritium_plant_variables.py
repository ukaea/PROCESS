t_div_tritium_residence: float = None
"""Tritium residence time in the divertor (s)"""

t_fw_tritium_residence: float = None
"""Tritium residence time in the first wall (s)"""

t_blkt_tritium_residence: float = None
"""Tritium residence time in the blanket (s)"""

t_heat_exchanger_tritium_residence: float = None
"""Tritium residence time in the heat exchanger (s)"""

t_tritium_extraction_system_tritium_residence: float = None
"""Tritium residence time in the tritium extraction system (s)"""

t_tritium_separation_membrane_tritium_residence: float = None
"""Tritium residence time in the tritium separation membrane (s)"""

t_plasma_fuelling_system_tritium_residence: float = None
"""Tritium residence time in the plasma fuelling system (s)"""

t_vacuum_pump_tritium_residence: float = None
"""Tritium residence time in the vacuum pump (s)"""

t_tritium_storage_tritium_residence: float = None
"""Tritium residence time in the tritium storage system (s)"""

t_isotope_separation_tritium_residence: float = None
"""Tritium residence time in the isotope separation system (s)"""

t_fuel_cleanup_tritium_residence: float = None
"""Tritium residence time in the fuel cleanup system (s)"""

t_detritiation_tritium_residence: float = None
"""Tritium residence time in the detritiation system (s)"""

m_plant_tritium_start_up: float = None
"""Mass of tritium at plant start-up (kg)"""

m_plant_tritium_start_up_minimum_required: float = None
"""Minimum mass of tritium required for plant start-up (kg)"""


def init_tritium_plant_variables():
    global \
        t_div_tritium_residence, \
        t_fw_tritium_residence, \
        t_blkt_tritium_residence, \
        t_heat_exchanger_tritium_residence, \
        t_tritium_extraction_system_tritium_residence, \
        t_tritium_separation_membrane_tritium_residence, \
        t_plasma_fuelling_system_tritium_residence, \
        t_vacuum_pump_tritium_residence, \
        t_tritium_storage_tritium_residence, \
        t_isotope_separation_tritium_residence, \
        t_fuel_cleanup_tritium_residence, \
        t_detritiation_tritium_residence, \
        m_plant_tritium_start_up, \
        m_plant_tritium_start_up_minimum_required

    t_div_tritium_residence = 100.0
    t_fw_tritium_residence = 50.0
    t_blkt_tritium_residence = 200.0
    t_heat_exchanger_tritium_residence = 20.0
    t_tritium_extraction_system_tritium_residence = 10.0
    t_tritium_separation_membrane_tritium_residence = 5.0
    t_plasma_fuelling_system_tritium_residence = 1.0
    t_vacuum_pump_tritium_residence = 0.5
    t_tritium_storage_tritium_residence = 500.0
    t_isotope_separation_tritium_residence = 1000.0
    t_fuel_cleanup_tritium_residence = 300.0
    t_detritiation_tritium_residence = 150.0
    m_plant_tritium_start_up = 10.0
    m_plant_tritium_start_up_minimum_required = 0.0
