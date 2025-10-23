i_vacuum_pumping: str = None
"""switch for vacuum pumping model:

 - ='old' for old detailed ETR model
 - ='simple' for simple steady-state model with comparison to ITER cryopumps
  !#TODO: old and simple not suitable names.
"""

n_iter_vacuum_pumps: float = None
"""number of high vacuum pumps (real number), each with the throughput of one
ITER cryopump (50 Pa m3 s-1), all operating at the same time (`i_vacuum_pumping='simple'`)
"""

i_vacuum_pump_type: int = None
"""switch for vacuum pump type:

   - =0 - for turbomolecular pump (magnetic bearing) with speed of 2.0 m3/s
     (1.95 for N2, 1.8 for He, 1.8 for DT)
   - =1 - for compound cryopump with nominal speed of 10.0 m3/s
     (9.0 for N2, 5.0 for He and 25.0 for DT)
"""

n_vv_vacuum_ducts: int = None
"""number of ducts (torus to pumps)"""

dlscal: float = None
"""vacuum system duct length scaling"""

pres_vv_chamber_base: float = None
"""base pressure during dwell before gas pre-fill(Pa)"""

pres_div_chamber_burn: float = None
"""divertor chamber pressure during burn (Pa)"""

molflow_vac_pumps: float = None
"""Pump throughput (molecules/s) (default is ITER value)"""

outgrat_fw: float = None
"""plasma chamber wall outgassing rate (Pa-m/s)"""

temp_vv_chamber_gas_burn_end: float = None
"""neutral gas temperature in chamber (K)"""

m_vv_vacuum_duct_shield: float = None
"""mass of vacuum duct shield (kg)"""

dia_vv_vacuum_ducts: float = None
"""diameter of duct passage (m)"""

n_vac_pumps_high: int = None
"""number of high vacuum pumps"""

i_vac_pump_dwell: int = None
"""switch for dwell pumping options:

   - =0 pumping only during t_plant_pulse_dwell
   - =1 pumping only during t_plant_pulse_coil_precharge
   - =2 pumping during t_plant_pulse_dwell + t_plant_pulse_coil_precharge

   The following are used in the Battes, Day and Rohde pump-down model
   See "Basic considerations on the pump-down time in the dwell phase of a pulsed fusion DEMO"
   http://dx.doi.org/10.1016/j.fusengdes.2015.07.011)(i_vacuum_pumping=simple')
"""

f_a_vac_pump_port_plasma_surface: float = None
"""area of one pumping port as a fraction of plasma surface area"""

volflow_vac_pumps_max: float = None
"""maximum pumping speed per unit area for deuterium & tritium, molecular flow"""

f_volflow_vac_pumps_impedance: float = None
"""effective pumping speed reduction factor due to duct impedance"""

pres_vv_chamber_dwell_start: float = None
"""initial neutral pressure at the beginning of the dwell phase (Pa)"""

outgasindex: float = None
"""outgassing decay index"""

outgasfactor: float = None
"""outgassing prefactor kw: outgassing rate at 1 s per unit area (Pa m s-1)"""


def init_vacuum_variables():
    """Initialise Vacuum variables"""
    global i_vacuum_pumping
    global n_iter_vacuum_pumps
    global i_vacuum_pump_type
    global n_vv_vacuum_ducts
    global dlscal
    global pres_vv_chamber_base
    global pres_div_chamber_burn
    global molflow_vac_pumps
    global outgrat_fw
    global temp_vv_chamber_gas_burn_end
    global m_vv_vacuum_duct_shield
    global dia_vv_vacuum_ducts
    global n_vac_pumps_high
    global i_vac_pump_dwell
    global f_a_vac_pump_port_plasma_surface
    global volflow_vac_pumps_max
    global f_volflow_vac_pumps_impedance
    global pres_vv_chamber_dwell_start
    global outgasindex
    global outgasfactor

    i_vacuum_pumping = "old"
    n_iter_vacuum_pumps = 0.0
    i_vacuum_pump_type = 1
    n_vv_vacuum_ducts = 0
    dlscal = 0.0
    pres_vv_chamber_base = 5.0e-4
    pres_div_chamber_burn = 0.36
    molflow_vac_pumps = 1.2155e22
    outgrat_fw = 1.3e-8
    temp_vv_chamber_gas_burn_end = 300.0
    m_vv_vacuum_duct_shield = 0.0
    dia_vv_vacuum_ducts = 0.0
    n_vac_pumps_high = 0
    i_vac_pump_dwell = 0
    f_a_vac_pump_port_plasma_surface = 0.0203
    volflow_vac_pumps_max = 27.3
    f_volflow_vac_pumps_impedance = 0.167
    pres_vv_chamber_dwell_start = 1.0
    outgasindex = 1.0
    outgasfactor = 0.0235
