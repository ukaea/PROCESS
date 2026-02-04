# These variables were from stellarator.f90
f_st_n_coils: float = None
"""Actual number of coils to reference value from stella_config file"""

f_st_rmajor: float = None
"""Actual major radius to reference value from stella_config file"""

f_st_aspect: float = None
"""Actual aspect ratio to reference value from stella_config file"""

f_st_coil_aspect:float = None
"""Scaling factor for (stellarator major radius / coil radius ratio)"""

f_st_b: float = None
"""Actual b_plasma_toroidal_on_axis to reference value from stella_config file """

f_st_i_total: float = None
"""Actual totail coil current to reference value from stella_config file"""

f_st_rminor: float = None
"""Actual minor radius to reference value from stella_config file"""

first_call: bool = None

first_call_stfwbs: bool = None

# These variables were from stellarator_variables.f90
istell: int = None
"""Switch for stellarator option (set via `device.dat`):
   - =0 use tokamak model
   - =1 use stellarator model: Helias5
   - =2 use stellarator model: Helias4
   - =3 use stellarator model: Helias3
   - =4 use stellarator model: Wendelstein 7-X with 50 Coils
   - =5 use stellarator model: Wendelstein 7-X with 30 Coils
   - =6 use stellarator model: Use stella_conf.json file (any modulear stellarator, see documentation)
"""

bmn: float = None
"""relative radial field perturbation"""

f_asym: float = None
"""divertor heat load peaking factor"""

f_rad: float = None
"""radiated power fraction in SOL"""

f_w: float = None
"""island size fraction factor"""

fdivwet: float = None
"""wetted fraction of the divertor area"""

flpitch: float = None
"""field line pitch (rad)"""

hportamax: float = None
"""maximum available area for horizontal ports (m2)"""

hportpmax: float = None
"""maximum available poloidal extent for horizontal ports (m)"""

hporttmax: float = None
"""maximum available toroidal extent for horizontal ports (m)"""

iotabar: float = None
"""rotational transform (reciprocal of tokamak q) for stellarator confinement time scaling laws"""

isthtr: int = None
"""Switch for stellarator auxiliary heating method:
   - = 1electron cyclotron resonance heating
   - = 2lower hybrid heating
   - = 3neutral beam injection
"""

m_res: int = None
"""poloidal resonance number (1)"""

max_gyrotron_frequency: float = None
"""Maximal available gyrotron frequency (input parameter) (Hz)"""

n_res: int = None
"""toroidal resonance number (1)"""

shear: float = None
"""magnetic shear, derivative of iotabar (1)"""

te0_ecrh_achievable: float = None
"""maximal central electron temperature as achievable by the ECRH, input. (keV)"""

vportamax: float = None
"""maximum available area for vertical ports (m2)"""

vportpmax: float = None
"""maximum available poloidal extent for vertical ports (m)"""

vporttmax: float = None
"""maximum available toroidal extent for vertical ports (m)"""

powerht_constraint: float = None

powerscaling_constraint: float = None


def init_stellarator_variables():
    global \
        first_call, \
        first_call_stfwbs, \
        f_st_n_coils, \
        f_st_rmajor, \
        f_st_rminor, \
        f_st_b, \
        f_st_i_total, \
        istell, \
        bmn, \
        f_asym, \
        f_rad, \
        f_w, \
        f_st_coil_aspect, \
        fdivwet, \
        flpitch, \
        hportamax, \
        hportpmax, \
        hporttmax, \
        iotabar, \
        isthtr, \
        m_res, \
        n_res, \
        shear, \
        vportamax, \
        vportpmax, \
        vporttmax, \
        max_gyrotron_frequency, \
        te0_ecrh_achievable, \
        powerht_constraint, \
        powerscaling_constraint

    first_call = True
    first_call_stfwbs = True
    f_st_n_coils = 0.0
    f_st_rmajor = 0.0
    f_st_rminor = 0.0
    f_st_b = 0.0
    f_st_i_total = 0.0
    istell = 0
    bmn = 1e-3
    f_asym = 1.0
    f_rad = 0.85
    f_w = 0.5
    fdivwet = 0.333333333333333
    flpitch = 1e-3
    hportamax = 0.0
    hportpmax = 0.0
    hporttmax = 0.0
    iotabar = 1.0
    isthtr = 1
    m_res = 5
    n_res = 5
    shear = 0.5
    vportamax = 0.0
    vportpmax = 0.0
    vporttmax = 0.0
    max_gyrotron_frequency = 1.0e9
    te0_ecrh_achievable = 1.0e2
    powerht_constraint = 0.0
    powerscaling_constraint = 0.0
