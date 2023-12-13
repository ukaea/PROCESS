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

    rs_1: List[float]
    """outer radial coordinates of shield"""
    zs_1: List[float]
    """outer vertical coordinates of shield"""
    rs_2: List[float]
    """inner radial coordinates of shield"""
    zs_2: List[float]
    """inner vertical coordinates of shield"""
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
) -> ShieldGeometry:
    """Calculates radial and vertical distances for the geometry of section of blanket above the midplane in a single null configuration

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
    :return: ShieldGeometry - dataclass returning shield radial and vertical coordinates
    :rtype: DataClass
    """

    # Side furthest from plasma
    kapx = cumulative_upper["shldtth"] / rminx_far
    rs_1, zs_1 = plotdh(radx_far, rminx_far, triang, kapx)

    # Side nearest to plasma
    kapx = (cumulative_upper["vvblgap"]) / rminx_near
    rs_2, zs_2 = plotdh(radx_near, rminx_near, triang, kapx)

    rs = np.concatenate([rs_1, rs_2[::-1]])
    zs = np.concatenate([zs_1, zs_2[::-1]])

    return ShieldGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )


def shield_geometry(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> ShieldGeometry:
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
    :return: ShieldGeometry - dataclass returning shield radial and vertical coordinates
    :rtype: DataClass
    """

    # Side furthest from plasma
    kapx = cumulative_lower["shldlth"] / rminx_far
    rs_1, zs_1 = plotdh(radx_far, rminx_far, triang, kapx)

    # Side nearest to plasma
    kapx = (cumulative_lower["divfix"]) / rminx_near
    rs_2, zs_2 = plotdh(radx_near, rminx_near, triang, kapx)

    rs = np.concatenate([rs_1, rs_2[::-1]])
    zs = np.concatenate([zs_1, zs_2[::-1]])

    return ShieldGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )
