"""Module containing variables for the first wall models"""

from dataclasses import dataclass

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class FirstWallData(PROCESSModelData):
    """Dataclass holding first wall variables"""

    a_fw_total_full_coverage: Parameter[float] = 0.0
    """First wall total surface area with no holes or ports [m^2]"""

    a_fw_inboard_full_coverage: Parameter[float] = 0.0
    """Inboard first wall surface area with no holes or ports [m^2]"""

    a_fw_outboard_full_coverage: Parameter[float] = 0.0
    """Outboard first wall surface area with no holes or ports [m^2]"""

    a_fw_total: Parameter[float] = 0.0
    """First wall total surface area [m^2]"""

    a_fw_inboard: Parameter[float] = 0.0
    """Inboard first wall surface area [m^2]"""

    a_fw_outboard: Parameter[float] = 0.0
    """Outboard first wall surface area [m^2]"""


CREATE_DICTS_FROM_DATACLASS = FirstWallData
