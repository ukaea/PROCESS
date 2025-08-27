"""
Calculate plasma elongation and radial and vertical coordinates for the geometry of the plasma
"""

import math
from dataclasses import dataclass

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
    rmajor: float,
    rminor: float,
    triang: float,
    kappa: float,
    i_single_null: int,
    i_plasma_shape: int,
    square: float,
) -> PlasmaGeometry:
    """
    Calculates radial and vertical distances and plasma elongation for the geometry of the plasma.

    This function computes the radial and vertical coordinates of the plasma boundary, as well as the plasma elongation,
    based on the given major radius, minor radius, triangularity, and elongation at 95% of the plasma surface. It also
    considers whether the plasma configuration is single null or double null.

    :param rmajor: Plasma major radius.
    :type rmajor: float
    :param rminor: Plasma minor radius.
    :type rminor: float
    :param triang: Plasma triangularity at separatrix.
    :type triang: float
    :param kappa: Plasma elongation at separatrix.
    :type kappa: float
    :param i_single_null: Switch for single null (1) or double null (0) plasma configuration.
    :type i_single_null: int
    :param i_plasma_shape: Switch for plasma shape (0 for double arc, 1 for Sauter).
    :type i_plasma_shape: int
    :param square: Square term for Sauter plasma shape.
    :type square: float
    :returns: A dataclass containing the plasma elongation and the radial and vertical coordinates of the plasma.
    :rtype: PlasmaGeometry

    """

    # Original PROCESS double arc plasma shape
    if i_plasma_shape == 0:
        x1 = (2.0 * rmajor * (1.0 + triang) - rminor * (triang**2 + kappa**2 - 1.0)) / (
            2.0 * (1.0 + triang)
        )
        x2 = (2.0 * rmajor * (triang - 1.0) - rminor * (triang**2 + kappa**2 - 1.0)) / (
            2.0 * (triang - 1.0)
        )
        r1 = 0.5 * math.sqrt(
            (rminor**2 * ((triang + 1.0) ** 2 + kappa**2) ** 2) / ((triang + 1.0) ** 2)
        )
        r2 = 0.5 * math.sqrt(
            (rminor**2 * ((triang - 1.0) ** 2 + kappa**2) ** 2) / ((triang - 1.0) ** 2)
        )
        theta1 = np.arcsin((kappa * rminor) / r1)
        theta2 = np.arcsin((kappa * rminor) / r2)
        inang = 1.0 / r1
        outang = 1.5 / r2
        if i_single_null == 0:
            angs1 = np.linspace(
                -(inang + theta1) + np.pi, (inang + theta1) + np.pi, 500, endpoint=True
            )
            angs2 = np.linspace(
                -(outang + theta2), (outang + theta2), 500, endpoint=True
            )
        else:
            angs1 = np.linspace(
                -theta1 + np.pi, (inang + theta1) + np.pi, 500, endpoint=True
            )
            angs2 = np.linspace(-(outang + theta2), theta2, 500, endpoint=True)

        xs1 = -(r1 * np.cos(angs1) - x1)
        ys1 = r1 * np.sin(angs1)
        xs2 = -(r2 * np.cos(angs2) - x2)
        ys2 = r2 * np.sin(angs2)

        rs = [xs1, xs2]
        zs = [ys1, ys2]

        return PlasmaGeometry(rs=rs, zs=zs, kappa=kappa)

    # Sauter plasma shape
    if i_plasma_shape == 1:
        x = np.linspace(-np.pi, np.pi, 256)

        # Sauter
        return PlasmaGeometry(
            rs=rmajor
            + rminor * np.cos(x + triang * np.sin(x) - square * np.sin(2 * x)),
            zs=kappa * rminor * np.sin(x + square * np.sin(2 * x)),
            kappa=kappa,
        )

    # Explicit return statement at the end of the function
    return None
