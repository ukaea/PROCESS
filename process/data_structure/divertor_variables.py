from dataclasses import dataclass


@dataclass
class DivertorData:
    """Divertor model dataclass"""

    anginc: float = 0.262
    """angle of incidence of field line on plate (rad)"""

    deg_div_field_plate: float = 1.0
    """field line angle wrt divertor target plate (degrees)"""

    betai: float = 1.0
    """poloidal plane angle between divertor plate and leg, inboard (rad)"""

    betao: float = 1.0
    """poloidal plane angle between divertor plate and leg, outboard (rad)"""

    f_vol_div_coolant: float = 0.3
    """divertor coolant fraction"""

    den_div_structure: float = 1.0e4
    """divertor structure density (kg/m3)"""

    dz_divertor: float = 0.2
    """divertor structure vertical thickness (m)"""

    m_div_plate: float = 0.0
    """divertor plate mass (kg)"""

    dx_div_plate: float = 0.035
    """divertor plate thickness (m) (from Spears, Sept 1990)"""

    a_div_surface_total: float = 0.0
    """divertor surface area (m2)"""

    fdiva: float = 1.11
    """divertor area fudge factor (for ITER, Sept 1990)"""

    f_div_flux_expansion: float = 2.0
    """The plasma flux expansion in the divertor (default 2; Wade 2020)"""

    pflux_div_heat_load_mw: float = 0.0
    """divertor heat load (MW/m2)"""

    i_div_heat_load: int = 2
    """switch for user input pflux_div_heat_load_mw:

    - = 0: divtart model turned off and user inputs pflux_div_heat_load_mw
    - = 1: divtart model calculates pflux_div_heat_load_mw
    - = 2: divwade model calculates pflux_div_heat_load_mw"""

    pflux_div_heat_load_max_mw: float = 5.0
    """heat load limit (MW/m2)"""

    prn1: float = 0.285
    """n-scrape-off / n-average plasma; (input for `i_plasma_pedestal=0`, = nd_plasma_separatrix_electron/nd_plasma_electrons_vol_avg if `i_plasma_pedestal>=1`)"""

    tdiv: float = 2.0
    """temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)"""

    xpertin: float = 2.0
    """perpendicular heat transport coefficient (m2/s)"""

    p_div_lower_nuclear_heat_mw: float = 0.0
    """Lower divertor neutron nuclear heat load on (MW)"""

    p_div_upper_nuclear_heat_mw: float = 0.0
    """Upper divertor neutron nuclear heat load on (MW)"""

    p_div_upper_rad_mw: float = 0.0
    """Upper divertor incident radiation power radiation power (MW)"""

    p_div_lower_rad_mw: float = 0.0
    """Lower divertor incident radiation power radiation power (MW)"""

    n_divertors: int = 2
    """Number of divertors (calculated from `i_single_null`)"""


CREATE_DICTS_FROM_DATACLASS = DivertorData
