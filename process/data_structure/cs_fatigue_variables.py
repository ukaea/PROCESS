from dataclasses import dataclass


@dataclass
class CSFatigueData:
    residual_sig_hoop: float = 2.4e8
    """residual hoop stress in strucutal material (Pa)"""

    n_cycle: float = 0.0
    """Allowable number of cycles for CS stress model"""

    n_cycle_min: float = 2.0e4
    """Minimum allowable number of cycles for CS stress model"""

    t_crack_radial: float = 6.0e-3
    """Initial depth of crack in thickness of conduit (m)"""

    t_crack_vertical: float = 0.89e-3
    """Inital vertical crack size (m)"""

    dr_cs_turn_conduit: float = 0.07
    """Thickness of CS conductor conduit (m)"""

    dz_cs_turn_conduit: float = 0.022
    """Vertical thickness of CS conductor conduit (m)"""

    bkt_life_csf: float = 0.0
    """Switch to pass bkt_life cycles to n_cycle_min"""

    sf_vertical_crack: float = 2.0
    """Safety factor for vertical crack size (-)"""

    sf_radial_crack: float = 2.0
    """Safety factor for radial crack size (-)"""

    sf_fast_fracture: float = 1.5
    """safety factor for stress intensity factor (-)"""

    paris_coefficient: float = 65.0e-14
    """Paris equation material coefficient (-)"""

    paris_power_law: float = 3.5
    """Paris equation material power law (-)"""

    walker_coefficient: float = 0.436
    """walker coefficent (-)"""

    fracture_toughness: float = 2.0e2
    """fracture toughness (MPa m^1/2)"""


CREATE_DICTS_FROM_DATACLASS = CSFatigueData
