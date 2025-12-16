anginc: float = None
"""angle of incidence of field line on plate (rad)"""

deg_div_field_plate: float = None
"""field line angle wrt divertor target plate (degrees)"""

betai: float = None
"""poloidal plane angle between divertor plate and leg, inboard (rad)"""

betao: float = None
"""poloidal plane angle between divertor plate and leg, outboard (rad)"""

f_vol_div_coolant: float = None
"""divertor coolant fraction"""

den_div_structure: float = None
"""divertor structure density (kg/m3)"""

dz_divertor: float = None
"""divertor structure vertical thickness (m)"""

m_div_plate: float = None
"""divertor plate mass (kg)"""

dx_div_plate: float = None
"""divertor plate thickness (m) (from Spears, Sept 1990)"""

a_div_surface_total: float = None
"""divertor surface area (m2)"""

fdiva: float = None
"""divertor area fudge factor (for ITER, Sept 1990)"""

f_div_flux_expansion: float = None
"""The plasma flux expansion in the divertor (default 2; Wade 2020)"""

pflux_div_heat_load_mw: float = None
"""divertor heat load (MW/m2)"""

i_div_heat_load: int = None
"""switch for user input pflux_div_heat_load_mw:

 - = 0: divtart model turned off and user inputs pflux_div_heat_load_mw
 - = 1: divtart model calculates pflux_div_heat_load_mw
 - = 2: divwade model calculates pflux_div_heat_load_mw"""

pflux_div_heat_load_max_mw: float = None
"""heat load limit (MW/m2)"""

prn1: float = None
"""n-scrape-off / n-average plasma; (input for `i_plasma_pedestal=0`, = nd_plasma_separatrix_electron/nd_plasma_electrons_vol_avg if `i_plasma_pedestal>=1`)"""

tdiv: float = None
"""temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)"""

xpertin: float = None
"""perpendicular heat transport coefficient (m2/s)"""


def init_divertor_variables():
    global \
        anginc, \
        deg_div_field_plate, \
        betai, \
        betao, \
        f_vol_div_coolant, \
        den_div_structure, \
        dz_divertor, \
        m_div_plate, \
        dx_div_plate, \
        a_div_surface_total, \
        fdiva, \
        f_div_flux_expansion, \
        pflux_div_heat_load_mw, \
        i_div_heat_load, \
        pflux_div_heat_load_max_mw, \
        prn1, \
        tdiv, \
        xpertin

    anginc = 0.262
    deg_div_field_plate = 1.0
    betai = 1.0
    betao = 1.0
    f_vol_div_coolant = 0.3
    den_div_structure = 1.0e4
    dz_divertor = 0.2
    m_div_plate = 0.0
    dx_div_plate = 0.035
    a_div_surface_total = 0.0
    fdiva = 1.11
    f_div_flux_expansion = 2.0
    pflux_div_heat_load_mw = 0.0
    i_div_heat_load = 2
    pflux_div_heat_load_max_mw = 5.0
    prn1 = 0.285
    tdiv = 2.0
    xpertin = 2.0
