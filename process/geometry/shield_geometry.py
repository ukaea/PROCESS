"""
Calculate radial and vertical coordinates for the geometry of the shield
"""
from dataclasses import dataclass
from typing import List
import numpy as np
from process.geometry.utils import plotdh


@dataclass
class ShieldGeometry:
    """Holds radial and vertical coordinates for the geometry of a shield"""

    rs: List[List[float]]
    """outer and inner radial coordinates of shield"""
    zs: List[List[float]]
    """outer and inner vertical coordinates of shield"""


def shield_geometry_single_null(
    cumulative_upper: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
    cumulative_lower: dict,
) -> ShieldGeometry:
    """Calculates radial and vertical distances for the geometry of shield in a single null configuration

    :param cumulative_upper: cumulative vertical thicknesses of components above the midplane
    :type cumulative_upper: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :return: ShieldGeometry - dataclass returning shield radial and vertical coordinates
    :rtype: DataClass
    """
    # Upper shield
    # Side furthest from plasma
    kapx = cumulative_upper["shldtth"] / rminx_far
    rs_upper_1, zs_upper_1 = plotdh(radx_far, rminx_far, triang, kapx)

    # Side nearest to plasma
    kapx = (cumulative_upper["vvblgap"]) / rminx_near
    rs_upper_2, zs_upper_2 = plotdh(radx_near, rminx_near, triang, kapx)

    # Lower shield
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = shield_geometry_lower(
        cumulative_lower=cumulative_lower,
        radx_far=radx_far,
        radx_near=radx_near,
        rminx_far=rminx_far,
        rminx_near=rminx_near,
        triang=triang,
    )

    rs = np.concatenate([rs_lower_2, rs_lower_1[::-1], rs_upper_1, rs_upper_2[::-1]])
    zs = np.concatenate([zs_lower_2, zs_lower_1[::-1], zs_upper_1, zs_upper_2[::-1]])
    return ShieldGeometry(
        rs=rs,
        zs=zs,
    )


def shield_geometry_lower(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> tuple:
    """Calculates radial and vertical distances for the geometry of section of shield below the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :return: radial and vertical coordinates for shield geometry below the midplane
    :rtype: tuple
    """
    # Side furthest from plasma
    kapx = cumulative_lower["shldlth"] / rminx_far
    rs_lower_1, zs_lower_1 = plotdh(radx_far, rminx_far, triang, kapx)

    # Side nearest to plasma
    kapx = (cumulative_lower["divfix"]) / rminx_near
    rs_lower_2, zs_lower_2 = plotdh(radx_near, rminx_near, triang, kapx)

    return rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2


def shield_geometry_double_null(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> ShieldGeometry:
    """Calculates radial and vertical distances for the geometry of shield in a double null configuration
    In a double null configuration, the geometry of the lower shield is reflected across the midplane to create the section of shield above the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :return: ShieldGeometry - dataclass returning shield radial and vertical coordinates
    :rtype: DataClass
    """
    # Lower shield
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = shield_geometry_lower(
        cumulative_lower=cumulative_lower,
        radx_far=radx_far,
        radx_near=radx_near,
        rminx_far=rminx_far,
        rminx_near=rminx_near,
        triang=triang,
    )

    rs_lower = np.concatenate([rs_lower_1, rs_lower_2[::-1]])
    zs_lower = np.concatenate([zs_lower_1, zs_lower_2[::-1]])

    # Upper shield
    rs_upper = rs_lower
    zs_upper = -1 * zs_lower

    rs = np.concatenate([rs_lower[::-1], rs_upper])
    zs = np.concatenate([zs_lower[::-1], zs_upper])

    return ShieldGeometry(rs=rs, zs=zs)
