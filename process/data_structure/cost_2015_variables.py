from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class Cost2015Data:

    """Average plant electrical output factoring in availability (MW)"""
    mean_electric_output: float = 0.0

    """Plant annual electrical output (MW.hrs)"""
    annual_electric_output: float = 0.0

    """Annual maintenance cost (M$)"""
    maintenance: float = 0.0

    """Total capital cost (M$)"""
    total_costs: float = 0.0

    s_label: list[str] = field(
        default_factory=lambda: np.array(["not used"] * 100, dtype=object)
    )

    """arrays for cost scaling values: k_ref, k, c_ref, cost, cost_factor"""
    s_kref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_k: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cost: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cost_factor: list[float] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )


CREATE_DICTS_FROM_DATACLASS = Cost2015Data
