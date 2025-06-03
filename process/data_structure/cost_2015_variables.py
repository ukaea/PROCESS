import numpy as np

mean_electric_output: float = 0.0

annual_electric_output: float = 0.0

maintenance: float = 0.0

total_costs: float = 0.0

s_label: np.ndarray = np.array(["not used"] * 100, dtype=object)

s_kref: np.ndarray = np.zeros(100, dtype=np.float64)

s_k: np.ndarray = np.zeros(100, dtype=np.float64)

s_cref: np.ndarray = np.zeros(100, dtype=np.float64)

s_cost: np.ndarray = np.zeros(100, dtype=np.float64)

s_cost_factor: np.ndarray = np.zeros(100, dtype=np.float64)
