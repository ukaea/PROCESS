"""
Calculate radial and vertical coordinates for the geometry of the first wall
"""
from dataclasses import dataclass
from typing import List
import numpy as np
from process.geometry.utils import plotdh, plotdhgap


#
@dataclass
class FirstWallGeometry:
    """Holds radial and vertical coordinates for the geometry of a first wall"""

    rs: List[List[float]]
    """radial first wall coordinates"""
    zs: List[List[float]]
    """vertical first wall coordinates"""


def first_wall_geometry_single_null(
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
    cumulative_upper: dict,
    triang: float,
    cumulative_lower: dict,
    blnktth: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
) -> FirstWallGeometry:
    """Calculates radial and vertical distances for the geometry of first wall in a single null configuration

    :param radx_outer: outboard radius of outer surface of first wall
    :type radx_outer: float
    :param rminx_outer: inboard radius of outer surface of first wall
    :type rminx_outer: float
    :param radx_inner: outboard radius of inner surface of first wall
    :type radx_inner: float
    :param rminx_inner: inboard radius of inner surface of first wall
    :type rminx_inner: float
    :param cumulative_upper: cumulative vertical thicknesses of components above the midplane
    :type cumulative_upper: dict
    :param triang: plasma triangularity
    :type triang: float
    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param blnktth: top blanket vertical thickness
    :type blnktth: float
    :param c_blnkith: inboard blanket vertical thickness
    :type c_blnkith: float
    :param c_fwoth: outboard first wall vertical thickness
    :type c_fwoth: float
    :param fwith: inboard first wall radial thickness
    :type fwith: float
    :param fwoth: outboard first wall radial thickness
    :type fwoth: float
    :param tfwvt: top first wall vertical thickness
    :type tfwvt: float
    :return: FirstWallGeometry - dataclass returning first wall radial and vertical coordinates
    :rtype: DataClass
    """
    # Upper first wall: outer surface
    kapx = cumulative_upper["fwtth"] / rminx_outer
    rs_1, zs_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    # Upper first wall: inner surface
    rs_2, zs_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    # Lower first wall
    divgap = cumulative_lower["divfix"]
    top_point = divgap + blnktth
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        fwith=fwith,
        fwoth=fwoth,
        tfwvt=tfwvt,
        top_point=top_point,
    )

    rs = np.concatenate([rs_1, rs_lower_2, rs_2[::-1], rs_lower_1[::-1]])
    zs = np.concatenate([zs_1, zs_lower_2, zs_2[::-1], zs_lower_1[::-1]])

    return FirstWallGeometry(rs=rs, zs=zs)


def first_wall_geometry_lower(
    triang: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
    top_point: float,
) -> tuple:
    """Calculates radial and vertical distances for the geometry of section of first wall below the midplane

    :param triang: plasma triangularity
    :type triang: float
    :param c_blnkith: inboard blanket vertical thickness
    :type c_blnkith: float
    :param c_fwoth: outboard first wall vertical thickness
    :type c_fwoth: float
    :param fwith: inboard first wall radial thickness
    :type fwith: float
    :param fwoth: outboard first wall radial thickness
    :type fwoth: float
    :param tfwvt: top first wall vertical thickness
    :type tfwvt: float
    :param top_point: top point for plotdhgap, equal to
    :type top_point: float
    :return: radial and vertical coordinates for first wall geometry below the midplane
    :rtype: tuple
    """
    # Lower first wall
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_blnkith,
        c_fwoth,
        fwith,
        fwoth,
        top_point,
        -tfwvt,
        triang,
    )
    rs_1 = np.concatenate([rs1, rs2[::-1]])
    zs_1 = np.concatenate([zs1, zs2[::-1]])
    rs_2 = np.concatenate([rs3, rs4[::-1]])
    zs_2 = -np.concatenate([zs3, zs4[::-1]])

    return rs_1, zs_1, rs_2, zs_2


def first_wall_geometry_double_null(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
) -> FirstWallGeometry:
    """Calculates radial and vertical distances for the geometry of first wall in a double null configuration
    In a double null configuration, the geometry of the lower first wall is reflected across the midplane to create the section of first wall above the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param triang: plasma triangularity
    :type triang: float
    :param blnktth: top blanket vertical thickness
    :type blnktth: float
    :param c_blnkith: inboard blanket vertical thickness
    :type c_blnkith: float
    :param c_fwoth: outboard first wall vertical thickness
    :type c_fwoth: float
    :param fwith: inboard first wall radial thickness
    :type fwith: float
    :param fwoth: outboard first wall radial thickness
    :type fwoth: float
    :param tfwvt: top first wall vertical thickness
    :type tfwvt: float
    :return: FirstWallGeometry - dataclass returning first wall radial and vertical coordinates
    :rtype: DataClass
    """
    # Lower first wall
    divgap = cumulative_lower["divfix"]
    top_point = divgap + blnktth
    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        fwith=fwith,
        fwoth=fwoth,
        tfwvt=tfwvt,
        top_point=top_point,
    )

    # Upper first wall
    top_point = -1 * top_point
    rs_upper_1, zs_upper_1, rs_upper_2, zs_upper_2 = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        fwith=fwith,
        fwoth=fwoth,
        tfwvt=tfwvt,
        top_point=top_point,
    )

    rs_1 = np.concatenate([rs_upper_1, rs_lower_1[::-1]])
    rs_2 = np.concatenate([rs_upper_2, rs_lower_2[::-1]])
    zs_1 = np.concatenate([zs_upper_1, zs_lower_1[::-1]])
    zs_2 = np.concatenate([zs_upper_2, zs_lower_2[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return FirstWallGeometry(rs=rs, zs=zs)
