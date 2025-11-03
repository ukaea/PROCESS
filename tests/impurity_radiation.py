'''
This is not a real test, it was made to manually compare some selected values'''

import process.impurity_radiation as impurity_radiation
from process.fortran import impurity_radiation_module
import numpy as np
from typing import NamedTuple

class PimpdenParam(NamedTuple):
    imp_element_index: int = 0
    ne: np.array = np.array
    te: np.array = np.array
    expected_pimpden: np.array = np.array

pimden_parameters = PimpdenParam(
    imp_element_index=0,
    ne=np.array([
        9.42593370e19,
        9.37237672e19,
        9.21170577e19,
        8.94392086e19,
        8.56902197e19,
        8.08700913e19,
        7.49788231e19,
        6.80164153e19,
        5.99828678e19,
        3.28986749e19,
    ]),
    te=np.array([
        27.73451868,
        27.25167194,
        25.82164396,
        23.50149071,
        20.39190536,
        16.64794796,
        12.50116941,
        8.31182764,
        4.74643357,
        0.1,
    ]),
    expected_pimpden=np.array([
        25483.040634309407,
        24983.364799017138,
        23519.36229676814,
        21187.36013272842,
        18173.71029818293,
        14685.542994819023,
        11005.497709894435,
        7448.7783515380615,
        4440.090318064716,
        294.54192663787137,
    ]),
)

impurity_radiation_module.init_impurity_radiation_module()
impurity_radiation.initialise_imprad()

# test = impurity_radiation.pimpden(0, pimden_parameters.ne, pimden_parameters.te)
# print(test)

temp = 1.4

test = impurity_radiation.pimpden(0, np.array([1e20]), np.array([temp]))
print(test)

no = 13
impurity_radiation_module.impurity_arr_frac[no - 1] = 1e-4
test = impurity_radiation.pimpden(no-1, np.array([1e20]), np.array([temp]))
print(test)

no = 14
impurity_radiation_module.impurity_arr_frac[no - 1] = 1e-5
test = impurity_radiation.pimpden(no-1, np.array([1e20]), np.array([temp]))
print(test)
