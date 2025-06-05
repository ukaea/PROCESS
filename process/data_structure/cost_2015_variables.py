import numpy as np

mean_electric_output: float

annual_electric_output: float

maintenance: float

total_costs: float

s_label: list[str]

s_kref: list[float]

s_k: list[float]

s_cref: list[float]

s_cost: list[float]

s_cost_factor: list[float]


def init_cost_2015_variables():
    global mean_electric_output
    mean_electric_output = 0.0

    global annual_electric_output
    annual_electric_output = 0.0

    global maintenance
    maintenance = 0.0

    global total_costs
    total_costs = 0.0

    global s_label
    s_label = np.array(["not used"] * 100, dtype=object)

    global s_kref
    s_kref = np.zeros(100, dtype=np.float64)

    global s_k
    s_k = np.zeros(100, dtype=np.float64)

    global s_cref
    s_cref = np.zeros(100, dtype=np.float64)

    global s_cost
    s_cost = np.zeros(100, dtype=np.float64)

    global s_cost_factor
    s_cost_factor = np.zeros(100, dtype=np.float64)
