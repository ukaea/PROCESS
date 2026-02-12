a_fw_total_full_coverage: float = None
"""First wall total surface area with no holes or ports [m^2]"""


a_fw_inboard_full_coverage: float = None
"""Inboard first wall surface area with no holes or ports [m^2]"""


a_fw_outboard_full_coverage: float = None
"""Outboard first wall surface area with no holes or ports [m^2]"""


a_fw_total: float = None
"""First wall total surface area [m^2]"""

a_fw_inboard: float = None
"""Inboard first wall surface area [m^2]"""

a_fw_outboard: float = None
"""Outboard first wall surface area [m^2]"""


def init_first_wall_variables():
    """Initializes first wall variables to None"""
    global \
        a_fw_total_full_coverage, \
        a_fw_inboard_full_coverage, \
        a_fw_outboard_full_coverage, \
        a_fw_total, \
        a_fw_inboard, \
        a_fw_outboard

    a_fw_total_full_coverage = 0.0
    a_fw_inboard_full_coverage = 0.0
    a_fw_outboard_full_coverage = 0.0
    a_fw_total = 0.0
    a_fw_inboard = 0.0
    a_fw_outboard = 0.0
