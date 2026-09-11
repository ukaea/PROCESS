"""Module containing variables for the divertor models"""

from dataclasses import dataclass
from enum import IntEnum, unique

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@unique
class DivertorHeatLoadModel(IntEnum):
    """Divertor heat load model enumeration, controlled' by `i_div_heat_load`"""

    USER_INPUT = 0
    """User input for divertor heat load"""

    PENG_CHAMBER = 1
    """Divertor heat load model based on Peng chamber"""

    WADE = 2
    """Divertor heat load model based on Wade (Wade 2020)"""


@dataclass(slots=True)
class DivertorData(PROCESSModelData):
    """Dataclass holding divertor variables"""

    anginc: Parameter[float] = 0.262
    """angle of incidence of field line on plate (rad)"""

    deg_div_field_plate: Parameter[float] = 1.0
    """field line angle wrt divertor target plate (degrees)"""

    betai: Parameter[float] = 1.0
    """poloidal plane angle between divertor plate and leg, inboard (rad)"""

    betao: Parameter[float] = 1.0
    """poloidal plane angle between divertor plate and leg, outboard (rad)"""

    f_vol_div_coolant: Parameter[float] = 0.3
    """divertor coolant fraction"""

    den_div_structure: Parameter[float] = 1.0e4
    """divertor structure density (kg/m3)"""

    dz_divertor: Parameter[float] = 0.2
    """divertor structure vertical thickness (m)"""

    m_div_plate: Parameter[float] = 0.0
    """divertor plate mass (kg)"""

    dx_div_plate: Parameter[float] = 0.035
    """divertor plate thickness (m) (from Spears, Sept 1990)"""

    a_div_surface_total: Parameter[float] = 0.0
    """divertor surface area (m2)"""

    fdiva: Parameter[float] = 1.11
    """divertor area fudge factor (for ITER, Sept 1990)"""

    f_div_flux_expansion: Parameter[float] = 2.0
    """The plasma flux expansion in the divertor (default 2; Wade 2020)"""

    pflux_div_heat_load_mw: Parameter[float] = 0.0
    """divertor heat load (MW/m2)"""

    i_div_heat_load: int = 2
    """switch for user input pflux_div_heat_load_mw:

    - = 0: User input
    - = 1: Peng chamber model
    - = 2: Wade model"""

    pflux_div_heat_load_max_mw: Parameter[float] = 5.0
    """heat load limit (MW/m2)"""

    prn1: Parameter[float] = 0.285
    """n-scrape-off / n-average plasma; (input for `i_plasma_pedestal=0`, = nd_plasma_separatrix_electron/nd_plasma_electrons_vol_avg if `i_plasma_pedestal>=1`)"""

    tdiv: Parameter[float] = 2.0
    """temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)"""

    xpertin: Parameter[float] = 2.0
    """perpendicular heat transport coefficient (m2/s)"""

    p_div_lower_nuclear_heat_mw: Parameter[float] = 0.0
    """Lower divertor neutron nuclear heat load on (MW)"""

    p_div_upper_nuclear_heat_mw: Parameter[float] = 0.0
    """Upper divertor neutron nuclear heat load on (MW)"""

    p_div_upper_rad_mw: Parameter[float] = 0.0
    """Upper divertor incident radiation power radiation power (MW)"""

    p_div_lower_rad_mw: Parameter[float] = 0.0
    """Lower divertor incident radiation power radiation power (MW)"""

    n_divertors: int = 2
    """Number of divertors (calculated from `i_single_null`)"""

    deg_div_poloidal_plasma: Parameter[float] = 0.0
    """Divertor poloidal angle subtended by plasma (degrees)"""


CREATE_DICTS_FROM_DATACLASS = DivertorData
