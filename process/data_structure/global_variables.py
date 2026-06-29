from dataclasses import dataclass, field


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

    xlabel: list[str] = field(default_factory=lambda: [""])
    """Scan parameters description label"""

    vlabel: list[str] = field(default_factory=lambda: [""])
    """Scan values name label"""

    convergence_parameter: float = 0.0
    """VMCON convergence parameter 'sum'"""


CREATE_DICTS_FROM_DATACLASS = GlobalData
