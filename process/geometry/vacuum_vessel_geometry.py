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
    cumulative_lower: dict,
    lower: dict,
) -> VacuumVesselGeometry:
    """Calculates radial and vertical distances for the geometry of vacuum vessel in a single null configuration

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
    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param lower: vertical thicknesses of components below the midplane
    :type lower: dict
    :return: VacuumVesselGeometry - dataclass returning vacuum vessel radial and vertical coordinates
    :rtype: DataClass
    """
    # Upper vacuum vessel
    kapx = cumulative_upper["d_vv_top"] / rminx_outer
    rs_upper_1, zs_upper_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    kapx = (
        float(cumulative_upper["d_vv_top"]) - float(upper["d_vv_top"])
    ) / rminx_inner
    rs_upper_2, zs_upper_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    # Lower vacuum vessel
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = vacuum_vessel_geometry_lower(
        cumulative_lower=cumulative_lower,
        lower=lower,
        radx_inner=radx_inner,
        radx_outer=radx_outer,
        rminx_inner=rminx_inner,
        rminx_outer=rminx_outer,
        triang=triang,
    )

    rs = np.concatenate([rs_lower_2, rs_lower_1[::-1], rs_upper_1, rs_upper_2[::-1]])
    zs = np.concatenate([zs_lower_2, zs_lower_1[::-1], zs_upper_1, zs_upper_2[::-1]])
    return VacuumVesselGeometry(
        rs=rs,
        zs=zs,
    )


def vacuum_vessel_geometry_lower(
    cumulative_lower: dict,
    lower: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> tuple:
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
    :return: radial and vertical coordinates for vacuum vessel geometry below the midplane
    :rtype: tuple
    """
    kapx = cumulative_lower["d_vv_bot"] / rminx_outer
    rs_lower_1, zs_lower_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    kapx = (
        float(cumulative_lower["d_vv_bot"]) + float(lower["d_vv_bot"])
    ) / rminx_inner
    rs_lower_2, zs_lower_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    return rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2


def vacuum_vessel_geometry_double_null(
    cumulative_lower: dict,
    lower: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> VacuumVesselGeometry:
    """Calculates radial and vertical distances for the geometry of vacuum vessel in a double null configuration
    In a double null configuration, the geometry of the lower vacuum vessel is reflected across the midplane to create the section of vacuum vessel above the midplane

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
    # Lower vacuum vessel
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = vacuum_vessel_geometry_lower(
        cumulative_lower=cumulative_lower,
        lower=lower,
        radx_inner=radx_inner,
        radx_outer=radx_outer,
        rminx_inner=rminx_inner,
        rminx_outer=rminx_outer,
        triang=triang,
    )
    rs_lower = np.concatenate([rs_lower_1, rs_lower_2[::-1]])
    zs_lower = np.concatenate([zs_lower_1, zs_lower_2[::-1]])

    # Upper vacuum vessel
    rs_upper = rs_lower
    zs_upper = -1 * zs_lower

    rs = np.concatenate([rs_lower[::-1], rs_upper])
    zs = np.concatenate([zs_lower[::-1], zs_upper])

    return VacuumVesselGeometry(rs=rs, zs=zs)
