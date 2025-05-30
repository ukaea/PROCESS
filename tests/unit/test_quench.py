import pytest
import numpy as np
from process.quench import _copper_irradiation_resistivity


def test_copper_irradiation(fluence, nu_data):
    fluence = [
        0.24, 0.45, 0.60, 0.79, 1.06, 1.26, 1.56, 1.89, 2.23, 2.68,
        2.92, 3.13, 3.30, 3.59, 3.90, 4.20, 4.71, 5.04, 5.34, 5.56,
        5.82, 6.03, 6.23, 6.58, 6.86, 7.38, 7.75, 8.11, 8.40, 8.51
    ]

    nu_irr_data = [
        16.32, 33.80, 44.12, 58.44, 76.54, 89.38, 106.13, 124.96, 142.84, 163.72,
        173.59, 181.75, 188.98, 200.51, 211.92, 221.67, 238.77, 249.15, 257.36, 262.83,
        268.51, 273.56, 278.82, 286.39, 293.10, 303.27, 309.51, 316.05, 320.58, 322.19
    ]

    nu_irr = _copper_irradiation_resistivity(fluence)
    assert(np.allclose(nu_irr, nu_irr_data))