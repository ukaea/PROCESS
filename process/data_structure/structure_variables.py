"""Module containing variables for the structure models"""

from dataclasses import dataclass

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class StructureData(PROCESSModelData):
    """Dataclass holding structure variables"""

    aintmass: Parameter[float] = 0.0
    """intercoil structure mass (kg)"""

    clgsmass: Parameter[float] = 0.0
    """gravity support structure for TF coil, PF coil and intercoil support systems (kg)"""

    coldmass: Parameter[float] = 0.0
    """total mass of components at cryogenic temperatures (kg)"""

    fncmass: Parameter[float] = 0.0
    """PF coil outer support fence mass (kg)"""

    gsmass: Parameter[float] = 0.0
    """reactor core gravity support mass (kg)"""


CREATE_DICTS_FROM_DATACLASS = StructureData
