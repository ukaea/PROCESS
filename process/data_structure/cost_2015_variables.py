"""Module containing variables for the costs 2015 models"""

from dataclasses import dataclass, field

import numpy as np

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class Cost2015Data(PROCESSModelData):
    """Dataclass holding cost 2015 variables"""

    mean_electric_output: Parameter[float] = 0.0

    annual_electric_output: Parameter[float] = 0.0

    maintenance: Parameter[float] = 0.0

    total_costs: Parameter[float] = 0.0

    s_label: list[str] = field(
        default_factory=lambda: np.array(["not used"] * 100, dtype=object)
    )

    s_kref: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )

    s_k: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )

    s_cref: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )

    s_cost: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )

    s_cost_factor: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros(100, dtype=np.float64)
    )


CREATE_DICTS_FROM_DATACLASS = Cost2015Data
