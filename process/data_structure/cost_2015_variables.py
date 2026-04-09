from dataclasses import dataclass, field

import numpy as np


@dataclass
class Cost2015Data:
    mean_electric_output: float = 0.0

    annual_electric_output: float = 0.0

    maintenance: float = 0.0

    total_costs: float = 0.0

    s_label: list[str] = field(
        default_factory=lambda: np.array(["not used"] * 100, dtype=object)
    )

    s_kref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_k: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cref: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cost: list[float] = field(default_factory=lambda: np.zeros(100, dtype=np.float64))

    s_cost_factor: list[float] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )


CREATE_DICTS_FROM_DATACLASS = Cost2015Data
