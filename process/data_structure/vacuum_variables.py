from dataclasses import dataclass


@dataclass
class VacuumData:
    i_vacuum_pumping: str = "old"
    """switch for vacuum pumping model:

    - ='old' for old detailed ETR model
    - ='simple' for simple steady-state model with comparison to ITER cryopumps
      !#TODO: old and simple not suitable names.
    """

    n_iter_vacuum_pumps: float = 0.0
    """number of high vacuum pumps (real number), each with the throughput of one
    ITER cryopump (50 Pa m3 s-1), all operating at the same time (`i_vacuum_pumping='simple'`)
    """

    i_vacuum_pump_type: int = 1
    """switch for vacuum pump type:

      - =0 - for turbomolecular pump (magnetic bearing) with speed of 2.0 m3/s
        (1.95 for N2, 1.8 for He, 1.8 for DT)
      - =1 - for compound cryopump with nominal speed of 10.0 m3/s
        (9.0 for N2, 5.0 for He and 25.0 for DT)
    """

    n_vv_vacuum_ducts: int = 0
    """number of ducts (torus to pumps)"""

    dlscal: float = 0.0
    """vacuum system duct length scaling"""

    pres_vv_chamber_base: float = 5.0e-4
    """base pressure during dwell before gas pre-fill(Pa)"""

    pres_div_chamber_burn: float = 0.36
    """divertor chamber pressure during burn (Pa)"""

    molflow_vac_pumps: float = 1.2155e22
    """Pump throughput (molecules/s) (default is ITER value)"""

    outgrat_fw: float = 1.3e-8
    """plasma chamber wall outgassing rate (Pa-m/s)"""

    temp_vv_chamber_gas_burn_end: float = 300.0
    """neutral gas temperature in chamber (K)"""

    m_vv_vacuum_duct_shield: float = 0.0
    """mass of vacuum duct shield (kg)"""

    dia_vv_vacuum_ducts: float = 0.0
    """diameter of duct passage (m)"""

    n_vac_pumps_high: int = 0
    """number of high vacuum pumps"""

    i_vac_pump_dwell: int = 0
    """switch for dwell pumping options:

      - =0 pumping only during t_plant_pulse_dwell
      - =1 pumping only during t_plant_pulse_coil_precharge
      - =2 pumping during t_plant_pulse_dwell + t_plant_pulse_coil_precharge

      The following are used in the Battes, Day and Rohde pump-down model
      See "Basic considerations on the pump-down time in the dwell phase of a pulsed fusion DEMO"
      http://dx.doi.org/10.1016/j.fusengdes.2015.07.011)(i_vacuum_pumping=simple')
    """

    f_a_vac_pump_port_plasma_surface: float = 0.0203
    """area of one pumping port as a fraction of plasma surface area"""

    volflow_vac_pumps_max: float = 27.3
    """maximum pumping speed per unit area for deuterium & tritium, molecular flow"""

    f_volflow_vac_pumps_impedance: float = 0.167
    """effective pumping speed reduction factor due to duct impedance"""

    pres_vv_chamber_dwell_start: float = 1.0
    """initial neutral pressure at the beginning of the dwell phase (Pa)"""

    outgasindex: float = 1.0
    """outgassing decay index"""

    outgasfactor: float = 0.0235
    """outgassing prefactor kw: outgassing rate at 1 s per unit area (Pa m s-1)"""


CREATE_DICTS_FROM_DATACLASS = VacuumData
