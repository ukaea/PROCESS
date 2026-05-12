from dataclasses import dataclass


@dataclass
class StructureData:
    """Data describing the structure."""

    aintmass: float = 0.0
    """intercoil structure mass (kg)"""

    clgsmass: float = 0.0
    """gravity support structure for TF coil, PF coil and intercoil support systems (kg)"""

    coldmass: float = 0.0
    """total mass of components at cryogenic temperatures (kg)"""

    fncmass: float = 0.0
    """PF coil outer support fence mass (kg)"""

    gsmass: float = 0.0
    """reactor core gravity support mass (kg)"""


CREATE_DICTS_FROM_DATACLASS = StructureData
