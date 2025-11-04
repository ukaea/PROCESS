# These variables were from stellarator.f90
f_n: float = None

f_r: float = None

f_aspect: float = None

f_st_coil_aspect:float = None

f_b: float = None

f_i: float = None

f_a: float = None

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
    global first_call
    global first_call_stfwbs
    global f_n
    global f_r
    global f_a
    global f_b
    global f_i
    global istell
    global bmn
    global f_asym
    global f_rad
    global f_w
    global f_st_coil_aspect
    global fdivwet
    global flpitch
    global hportamax
    global hportpmax
    global hporttmax
    global iotabar
    global isthtr
    global m_res
    global n_res
    global shear
    global vportamax
    global vportpmax
    global vporttmax
    global max_gyrotron_frequency
    global te0_ecrh_achievable

    first_call = True
    first_call_stfwbs = True
    f_n = 0.0
    f_r = 0.0
    f_a = 0.0
    f_b = 0.0
    f_i = 0.0
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
