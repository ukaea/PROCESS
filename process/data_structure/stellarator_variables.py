"""Module containing variables for the stellarator models"""

from dataclasses import dataclass
from enum import IntEnum, unique


@unique
class StellaratorModel(IntEnum):
    """Enum for stellarator models"""

    DISABLED = 0
    """Disabled stellarator model (use tokamak model)"""

    HELIAS_5 = 1
    """Helias5 stellarator model"""

    HELIAS_4 = 2
    """Helias4 stellarator model"""

    HELIAS_3 = 3
    """Helias3 stellarator model"""

    W7X_50_COILS = 4
    """Wendelstein 7-X stellarator model with 50 coils"""

    W7X_30_COILS = 5
    """Wendelstein 7-X stellarator model with 30 coils"""

    CUSTOM = 6
    """Custom stellarator model using stella_conf.json file"""


@dataclass(slots=True)
class StellaratorData:
    """Dataclass holding stellarator variables"""

    f_st_n_coils: float = 0.0
    """Actual number of coils to reference value from stella_config file"""

    f_st_rmajor: float = 0.0
    """Actual major radius to reference value from stella_config file"""

    f_st_aspect: float = 0.0
    """Actual aspect ratio to reference value from stella_config file"""

    f_st_coil_aspect: float = 1.0
    """Scaling factor for (stellarator major radius / coil radius ratio)"""

    f_st_b: float = 0.0
    """Actual b_plasma_toroidal_on_axis to reference value from stella_config file """

    f_st_i_total: float = 0.0
    """Actual total coil current to reference value from stella_config file"""

    f_st_rminor: float = 0.0
    """Actual minor radius to reference value from stella_config file"""

    f_coil_shape: float = 0.0
    """Parameter required for coil scaling
    (min_plasma_coil_distance + stella_config_rminor_ref) / stella_config_coil_rminor
    """

    first_call: bool = True

    first_call_stfwbs: bool = True

    r_coil_minor: float = 0.0
    """Coil minor radius (m)"""

    r_coil_major: float = 0.0
    """Coil major radius (m)"""

    istell: int = 0
    """Switch for stellarator option `StellaratorModel`:
        - `StellaratorModel.DISABLED` use tokamak model
        - `StellaratorModel.HELIAS_5` use stellarator model: Helias5
        - `StellaratorModel.HELIAS_4` use stellarator model: Helias4
        - `StellaratorModel.HELIAS_3` use stellarator model: Helias3
        - `StellaratorModel.W7X_50_COILS` use stellarator model: Wendelstein 7-X with 50 Coils
        - `StellaratorModel.W7X_30_COILS` use stellarator model: Wendelstein 7-X with 30 Coils
        - `StellaratorModel.CUSTOM` use stellarator model: Use stella_conf.json file (any modular stellarator, see documentation)
    """

    bmn: float = 1e-3
    """relative radial field perturbation"""

    f_asym: float = 1.0
    """divertor heat load peaking factor"""

    f_rad: float = 0.85
    """radiated power fraction in SOL"""

    f_w: float = 0.5
    """island size fraction factor"""

    fdivwet: float = 0.333333333333333
    """wetted fraction of the divertor area"""

    flpitch: float = 1e-3
    """field line pitch (rad)"""

    hportamax: float = 0.0
    """maximum available area for horizontal ports (m2)"""

    hportpmax: float = 0.0
    """maximum available poloidal extent for horizontal ports (m)"""

    hporttmax: float = 0.0
    """maximum available toroidal extent for horizontal ports (m)"""

    iotabar: float = 1.0
    """rotational transform (reciprocal of tokamak q) for stellarator confinement time scaling laws"""

    isthtr: int = 1
    """Switch for stellarator auxiliary heating method:
        - = 1electron cyclotron resonance heating
        - = 2lower hybrid heating
        - = 3neutral beam injection
    """

    m_res: int = 5
    """poloidal resonance number (1)"""

    max_gyrotron_frequency: float = 1.0e9
    """Maximal available gyrotron frequency (input parameter) (Hz)"""

    n_res: int = 5
    """toroidal resonance number (1)"""

    shear: float = 0.5
    """magnetic shear, derivative of iotabar (1)"""

    te0_ecrh_achievable: float = 1.0e2
    """maximal central electron temperature as achievable by the ECRH, input. (keV)"""

    vportamax: float = 0.0
    """maximum available area for vertical ports (m2)"""

    vportpmax: float = 0.0
    """maximum available poloidal extent for vertical ports (m)"""

    vporttmax: float = 0.0
    """maximum available toroidal extent for vertical ports (m)"""

    powerht_constraint: float = 0.0

    powerscaling_constraint: float = 0.0


CREATE_DICTS_FROM_DATACLASS = StellaratorData
