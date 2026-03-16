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

f_blkt_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the blanket due to non-radiative processes"""

f_tritium_extraction_system_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the tritium extraction system due to non-radiative processes"""

f_fw_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the first wall due to non-radiative processes"""

f_div_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the divertor due to non-radiative processes"""

f_heat_exchanger_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the heat exchanger due to non-radiative processes"""

f_detritiation_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the detritiation system due to non-radiative processes"""

f_vacuum_pump_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the vacuum pump due to non-radiative processes"""

f_fuel_cleanup_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the fuel cleanup system due to non-radiative processes"""

f_isotope_separation_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the isotope separation system due to non-radiative processes"""

f_tritium_storage_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the tritium storage system due to non-radiative processes"""

f_fuelling_system_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the plasma fuelling system due to non-radiative processes"""

f_tritium_separation_membrane_non_rad_tritium_loss: float = None
"""Fraction of tritium lost from the tritium separation membrane due to non-radiative processes"""

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
        f_blkt_non_rad_tritium_loss, \
        f_tritium_extraction_system_non_rad_tritium_loss, \
        f_fw_non_rad_tritium_loss, \
        f_div_non_rad_tritium_loss, \
        f_heat_exchanger_non_rad_tritium_loss, \
        f_detritiation_non_rad_tritium_loss, \
        f_vacuum_pump_non_rad_tritium_loss, \
        f_fuel_cleanup_non_rad_tritium_loss, \
        f_isotope_separation_non_rad_tritium_loss, \
        f_tritium_storage_non_rad_tritium_loss, \
        f_fuelling_system_non_rad_tritium_loss, \
        f_tritium_separation_membrane_non_rad_tritium_loss, \
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
    f_blkt_non_rad_tritium_loss = 1e-4
    f_tritium_extraction_system_non_rad_tritium_loss = 1e-4
    f_fw_non_rad_tritium_loss = 0.0
    f_div_non_rad_tritium_loss = 0.0
    f_heat_exchanger_non_rad_tritium_loss = 1e-4
    f_detritiation_non_rad_tritium_loss = 1e-4
    f_vacuum_pump_non_rad_tritium_loss = 1e-4
    f_fuel_cleanup_non_rad_tritium_loss = 1e-4
    f_isotope_separation_non_rad_tritium_loss = 1e-4
    f_tritium_storage_non_rad_tritium_loss = 0.0
    f_fuelling_system_non_rad_tritium_loss = 1e-4
    f_tritium_separation_membrane_non_rad_tritium_loss = 1e-4
    m_plant_tritium_start_up = 10.0
    m_plant_tritium_start_up_minimum_required = 0.0
