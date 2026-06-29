"""Module containing variables required to perform a parameter scan"""

from dataclasses import dataclass, field

import numpy as np

from process.core.exceptions import ProcessValueError

IPNSCNS = 1000
"""Maximum number of scan points"""


@dataclass(slots=True)
class ScanData:
    """Dataclass holding scan variables"""

    scan_dim: int = 1
    """1-D or 2-D scan switch (1=1D, 2=2D)"""

    isweep: list[int] | int = 1
    """Number of scan points to calculate"""

    nsweep: list[int] | int | None = None
    """Switch denoting quantity to scan
    see `process.core.scan.ScanVariables` for available options
    """

    sweep: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=np.float64))
    """Actual values to use in scan"""

    def __post_init__(self):
        if isinstance(self.isweep, int):
            # avoid old 0 default
            self.isweep = [self.isweep or 1]

        if len(self.isweep) > 2 or len(self.sweep.shape) > 2:
            raise NotImplementedError("N-D Scans not currently supported")

        if max(self.isweep) > IPNSCNS:
            raise ProcessValueError(
                "Illegal value of isweep",
                isweep=self.isweep,
                IPNSCNS=IPNSCNS,
            )
        if self.nsweep != len(self.isweep):
            raise ValueError(
                "Number of sweep variables not equal to scan point dimensions"
            )
        if self.sweep.shape != self.isweep:
            if self.isweep != 1:
                self.isweep = list(self.sweep.shape)
            else:
                print("Unset sweep values set to zero")
                # TODO append to size instead of resetting
                self.sweep = np.zeros(self.isweep, dtype=np.float64)

        self.nsweep = np.asarray(self.nsweep, dtype=int)
        self.isweep = np.asarray(self.isweep, dtype=int)


CREATE_DICTS_FROM_DATACLASS = ScanData
