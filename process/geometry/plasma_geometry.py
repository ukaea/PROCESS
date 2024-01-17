"""
Calculate plasma elongation and radial and vertical coordinates for the geometry of the plasma
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class PlasmaGeometry:
    """Holds radial and vertical coordinates for the geometry of a plasma"""

    rs: np.ndarray
    """radial plasma coordinates"""
    zs: np.ndarray
    """vertical plasma coordinates"""
    kappa: float
    """plasma separatrix elongation"""


def plasma_geometry(
    r_0: float, a: float, triang_95: float, kappa_95: float, i_single_null: bool
) -> PlasmaGeometry:
    """Calculates radial and vertical distances and plasma elongation for the geometry of the plasma

    :param r_0: plasma major radius
    :type r_0: float
    :param a: plasma minor radius
    :type a: float
    :param triang_95: plasma triangularity at 95% surface
    :type triang_95: float
    :param kappa_95: plasma elongation at 95% surface
    :type kappa_95: float
    :param i_single_null: switch for single null / double null plasma
    :type i_single_null: bool
    :return: dataclass returning plasma elongation and radial and vertical coordinates of the plasma
    :rtype: PlasmaGeometry
    """
    delta = 1.5 * triang_95
    kappa = (1.1 * kappa_95) + 0.04

    x1 = (2.0 * r_0 * (1.0 + delta) - a * (delta**2 + kappa**2 - 1.0)) / (
        2.0 * (1.0 + delta)
    )
    x2 = (2.0 * r_0 * (delta - 1.0) - a * (delta**2 + kappa**2 - 1.0)) / (
        2.0 * (delta - 1.0)
    )
    r1 = 0.5 * math.sqrt(
        (a**2 * ((delta + 1.0) ** 2 + kappa**2) ** 2) / ((delta + 1.0) ** 2)
    )
    r2 = 0.5 * math.sqrt(
        (a**2 * ((delta - 1.0) ** 2 + kappa**2) ** 2) / ((delta - 1.0) ** 2)
    )
    theta1 = np.arcsin((kappa * a) / r1)
    theta2 = np.arcsin((kappa * a) / r2)
    inang = 1.0 / r1
    outang = 1.5 / r2
    if i_single_null == 0:
        angs1 = np.linspace(
            -(inang + theta1) + np.pi, (inang + theta1) + np.pi, 256, endpoint=True
        )
        angs2 = np.linspace(-(outang + theta2), (outang + theta2), 256, endpoint=True)
    else:
        angs1 = np.linspace(
            -theta1 + np.pi, (inang + theta1) + np.pi, 256, endpoint=True
        )
        angs2 = np.linspace(-(outang + theta2), theta2, 256, endpoint=True)

    xs1 = -(r1 * np.cos(angs1) - x1)
    ys1 = r1 * np.sin(angs1)
    xs2 = -(r2 * np.cos(angs2) - x2)
    ys2 = r2 * np.sin(angs2)

    rs = [xs1, xs2]
    zs = [ys1, ys2]

    return PlasmaGeometry(rs=rs, zs=zs, kappa=kappa)
