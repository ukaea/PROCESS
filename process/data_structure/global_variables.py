from dataclasses import dataclass


@dataclass(slots=True)
class GlobalData:
    icase: str = "Steady-state tokamak model"
    """Power plant type"""

    runtitle: str = "Run Title (change this line using input variable 'runtitle')"
    """A short descriptive title for the run"""

    maxcal: int = 200
    """Maximum number of solver iterations"""

    fileprefix: str = ""
    """Path to input file"""

    output_prefix: str = ""
    """Output file path prefix"""

    xlabel: str = ""
    """Scan parameter description label"""

    vlabel: str = ""
    """Scan value name label"""

    xlabel_2: str = ""
    """Scan parameter description label (2nd dimension)"""

    vlabel_2: str = ""
    """Scan value name label (2nd dimension)"""

    iscan_global: int = 0

    convergence_parameter: float = 0.0
    """VMCON convergence parameter 'sum'"""


CREATE_DICTS_FROM_DATACLASS = GlobalData
