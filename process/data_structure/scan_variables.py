"""Module containing variables required to perform a parameter scan"""

from dataclasses import dataclass, field

import numpy as np

IPNSCNS = 1000
"""Maximum number of scan points"""


IPNSCNV = 82
"""Number of available scan variables"""


NOUTVARS = 84


@dataclass(slots=True)
class ScanData:
    """Dataclass holding scan variables"""

    scan_dim: int = 1
    """1-D or 2-D scan switch (1=1D, 2=2D)"""

    isweep: int = 0
    """Number of scan points to calculate"""

    isweep_2: int = 0
    """Number of 2D scan points to calculate"""

    nsweep: int = 1
    """Switch denoting quantity to scan
    """

    nsweep_2: int = 3
    """nsweep_2 /3/ : switch denoting quantity to scan for 2D scan:"""

    sweep: list[float] = field(
        default_factory=lambda: np.zeros(IPNSCNS, dtype=np.float64)
    )
    """sweep(IPNSCNS) /../: actual values to use in scan"""

    sweep_2: list[float] = field(
        default_factory=lambda: np.zeros(IPNSCNS, dtype=np.float64)
    )
    """sweep_2(IPNSCNS) /../: actual values to use in 2D scan"""

    # Vars in subroutines scan_1d and scan_2d requiring re-initialising before
    # each new run

    first_call_1d: bool = True

    first_call_2d: bool = True


CREATE_DICTS_FROM_DATACLASS = ScanData
