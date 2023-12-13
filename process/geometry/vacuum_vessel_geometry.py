"""
Calculate radial and vertical coordinates for the geometry of the vacuum vessel
"""
from typing import List
from dataclasses import dataclass
import numpy as np
from process.geometry.utils import plotdh


@dataclass
class VacuumVesselGeometry:
    """Holds radial and vertical coordinates for the geometry of a vacuum vessel"""

    rs_1: List[float]
    """outer radial coordinates of vacuum vessel"""
    zs_1: List[float]
    """outer vertical coordinates of vacuum vessel"""
    rs_2: List[float]
    """inner radial coordinates of vacuum vessel"""
    zs_2: List[float]
    """inner vertical coordinates of vacuum vessel"""
    rs: List[List[float]]
    """outer and inner radial coordinates of vacuum vessel"""
    zs: List[List[float]]
    """outer and inner vertical coordinates of vacuum vessel"""


def vacuum_vessel_geometry_single_null(
    cumulative_upper: dict,
    upper: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> VacuumVesselGeometry:
    """Calculates radial and vertical distances for the geometry of section of vacuum vessel above the midplane in a single null configuration

    :param cumulative_upper: cumulative vertical thicknesses of components above the midplane
    :type cumulative_upper: dict
    :param upper: vertical thicknesses of components above the midplane
    :type upper: dict
    :param triang: plasma triangularity
    :type triang: float
    :param radx_outer: outboard radius of outer surface of vacuum vessel
    :type radx_outer: float
    :param rminx_outer: inboard radius of outer surface of vacuum vessel
    :type rminx_outer: float
    :param radx_inner: outboard radius of inner surface of vacuum vessel
    :type radx_inner: float
    :param rminx_inner: inboard radius of inner surface of vacuum vessel
    :type rminx_inner: float
    :return: VacuumVesselGeometry - dataclass returning vacuum vessel radial and vertical coordinates
    :rtype: DataClass
    """

    kapx = cumulative_upper["d_vv_top"] / rminx_outer
    rs_1, zs_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    kapx = (
        float(cumulative_upper["d_vv_top"]) - float(upper["d_vv_top"])
    ) / rminx_inner
    rs_2, zs_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    rs = np.concatenate([rs_1, rs_2[::-1]])
    zs = np.concatenate([zs_1, zs_2[::-1]])

    return VacuumVesselGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )


def vacuum_vessel_geometry(
    cumulative_lower: dict,
    lower: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> VacuumVesselGeometry:
    """Calculates radial and vertical distances for the geometry of section of vacuum vessel below the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param lower: vertical thicknesses of components below the midplane
    :type lower: dict
    :param triang: plasma triangularity
    :type triang: float
    :param radx_outer: outboard radius of outer surface of vacuum vessel
    :type radx_outer: float
    :param rminx_outer: inboard radius of outer surface of vacuum vessel
    :type rminx_outer: float
    :param radx_inner: outboard radius of inner surface of vacuum vessel
    :type radx_inner: float
    :param rminx_inner: inboard radius of inner surface of vacuum vessel
    :type rminx_inner: float
    :return: VacuumVesselGeometry - dataclass returning vacuum vessel radial and vertical coordinates
    :rtype: DataClass
    """

    kapx = cumulative_lower["d_vv_bot"] / rminx_outer
    rs_1, zs_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    kapx = (
        float(cumulative_lower["d_vv_bot"]) + float(lower["d_vv_bot"])
    ) / rminx_inner
    rs_2, zs_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    rs = np.concatenate([rs_1, rs_2[::-1]])
    zs = np.concatenate([zs_1, zs_2[::-1]])

    return VacuumVesselGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )
