"""Unit tests for the impurity_radiation.f90.py module."""

from typing import NamedTuple

import numpy as np
import pytest

import process.models.physics.impurity_radiation as impurity_radiation
from process.data_structure import impurity_radiation_module


@pytest.fixture(autouse=True)
def initialise_impurity_radiation():
    impurity_radiation_module.init_impurity_radiation_module()
    impurity_radiation.initialise_imprad()


class PimpdenParam(NamedTuple):
    imp_element_index: int = 0
    ne: np.array = np.array
    te: np.array = np.array
    expected_pimpden: np.array = np.array


def test_pimpden():
    """Tests `pimpden` function.

    :param imp_element_index: impurity element
    :type imp_element_index: float

    :param ne:  electron density (/m3).
    :type ne: float

    :param te: electron temperature (keV)
    :type te: float

    :param expected_pimpden: Total impurity radiation density (W/m3)
    :type expected_pimpden: float
    """
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

    pimpden = impurity_radiation.pimpden(
        pimden_parameters.imp_element_index, pimden_parameters.ne, pimden_parameters.te
    )
    assert pytest.approx(pimpden) == pimden_parameters.expected_pimpden


class FradcoreParam(NamedTuple):
    rho: np.array = np.array
    radius_plasma_core_norm: float = 0.0
    f_p_plasma_core_rad_reduction: float = 0.0
    expected_fradcore: np.array = np.array


def test_fradcore():
    """Tests `fradcore` function.

    :param rho: normalised minor radius
    :type rho: np.array

    :param radius_plasma_core_norm:  normalised core radius
    :type radius_plasma_core_norm: float

    :param f_p_plasma_core_rad_reduction: fraction of core radiation
    :type f_p_plasma_core_rad_reduction: float
    :param expected_fradcore: Function to calculate core radiation fraction
    """
    fradcoreparam = FradcoreParam(
        rho=np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
        radius_plasma_core_norm=0.75000000000000011,
        f_p_plasma_core_rad_reduction=0.60000000000000009,
        expected_fradcore=np.array([
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.6,
            0.0,
            0.0,
        ]),
    )
    fradcore = impurity_radiation.fradcore(
        fradcoreparam.rho,
        fradcoreparam.radius_plasma_core_norm,
        fradcoreparam.f_p_plasma_core_rad_reduction,
    )
    assert pytest.approx(fradcore) == fradcoreparam.expected_fradcore


class ZavofteParam(NamedTuple):
    imp_element_index: int = 0
    te: np.array = np.array
    expected_zav_of_te: np.array = np.array


def test_zav_of_te():
    """Tests `Zav_of_te` function.

    :param imp_element_index: impurity element
    :type imp_element_index: float

    :param te:  electron temperature (keV)
    :type te: np.array

    :param expected_zav_of_te: Electron temperature dependent average atomic number
    :type expected_zav_of_te: np.array
    """
    zavofteparam = ZavofteParam(
        imp_element_index=0,
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
        expected_zav_of_te=np.array([
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
            1.00000000000001,
        ]),
    )
    zav_of_te = impurity_radiation.zav_of_te(
        zavofteparam.imp_element_index, zavofteparam.te
    )

    assert pytest.approx(zav_of_te) == zavofteparam.expected_zav_of_te
