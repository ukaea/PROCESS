residual_sig_hoop: float = None
"""residual hoop stress in strucutal material (Pa)"""

n_cycle: float = None
"""Allowable number of cycles for CS stress model"""

n_cycle_min: float = None
"""Minimum allowable number of cycles for CS stress model"""

t_crack_radial: float = None
"""Initial depth of crack in thickness of conduit (m)"""

t_crack_vertical: float = None
"""Inital vertical crack size (m)"""

dr_cs_turn_conduit: float = None
"""Thickness of CS conductor conduit (m)"""

dz_cs_turn_conduit: float = None
"""Vertical thickness of CS conductor conduit (m)"""

bkt_life_csf: float = None
"""Switch to pass bkt_life cycles to n_cycle_min"""

sf_vertical_crack: float = None
"""Safety factor for vertical crack size (-)"""

sf_radial_crack: float = None
"""Safety factor for radial crack size (-)"""

sf_fast_fracture: float = None
"""safety factor for stress intensity factor (-)"""

paris_coefficient: float = None
"""Paris equation material coefficient (-)"""

paris_power_law: float = None
"""Paris equation material power law (-)"""

walker_coefficient: float = None
"""walker coefficent (-)"""

fracture_toughness: float = None
"""fracture toughness (MPa m^1/2)"""


def init_cs_fatigue_variables():
    global residual_sig_hoop
    residual_sig_hoop = 2.4e8

    global t_crack_radial
    t_crack_radial = 6.0e-3

    global t_crack_vertical
    t_crack_vertical = 0.89e-3

    global n_cycle
    n_cycle = 0.0

    global n_cycle_min
    n_cycle_min = 2.0e4

    global dz_cs_turn_conduit
    dz_cs_turn_conduit = 0.022

    global dr_cs_turn_conduit
    dr_cs_turn_conduit = 0.07

    global bkt_life_csf
    bkt_life_csf = 0.0

    global sf_vertical_crack
    sf_vertical_crack = 2.0

    global sf_radial_crack
    sf_radial_crack = 2.0

    global sf_fast_fracture
    sf_fast_fracture = 1.5

    global paris_coefficient
    paris_coefficient = 65.0e-14

    global paris_power_law
    paris_power_law = 3.5

    global walker_coefficient
    walker_coefficient = 0.436

    global fracture_toughness
    fracture_toughness = 2.0e2
