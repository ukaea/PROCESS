from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class Cost2015Data:

    mean_electric_output: float = 0.0
    """Average plant electrical output factoring in availability (MW)"""

    annual_electric_output: float = 0.0
    """Plant annual electrical output (MW.hrs)"""

    maintenance: float = 0.0
    """Annual maintenance cost (M$)"""

    total_costs: float = 0.0
    """Total capital cost (M$)"""

    s_label: list[str] = field(
        default_factory=lambda: np.array(["not used"] * 100, dtype=object)
    )

    s_kref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))
    """Array for cost scaling value reference"""

    s_k: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))
    """Array for cost scaling value calculated"""

    s_cref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))
    """Array for cost reference values"""

    s_cost: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))
    """Array for cost values"""

    s_cost_factor: list[float] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )
    """Array for cost scaling factors"""


CREATE_DICTS_FROM_DATACLASS = Cost2015Data
