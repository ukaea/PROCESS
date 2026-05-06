from dataclasses import dataclass


@dataclass
class FirstWallData:
    a_fw_total_full_coverage: float = 0.0
    """First wall total surface area with no holes or ports [m^2]"""

    a_fw_inboard_full_coverage: float = 0.0
    """Inboard first wall surface area with no holes or ports [m^2]"""

    a_fw_outboard_full_coverage: float = 0.0
    """Outboard first wall surface area with no holes or ports [m^2]"""

    a_fw_total: float = 0.0
    """First wall total surface area [m^2]"""

    a_fw_inboard: float = 0.0
    """Inboard first wall surface area [m^2]"""

    a_fw_outboard: float = 0.0
    """Outboard first wall surface area [m^2]"""


CREATE_DICTS_FROM_DATACLASS = FirstWallData
