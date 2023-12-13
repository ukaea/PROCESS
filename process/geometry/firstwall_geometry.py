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
) -> FirstWallGeometry:
    """Calculates radial and vertical distances for the geometry of section of first wall above the midplane in a single null configuration

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
    :return: FirstWallGeometry - dataclass returning first wall radial and vertical coordinates
    :rtype: DataClass
    """
    # Upper first wall: outer surface
    kapx = cumulative_upper["fwtth"] / rminx_outer
    rs_1, zs_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    # Upper first wall: inner surface
    rs_2, zs_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    rs_3 = np.concatenate([rs_1, rs_2[::-1]])
    zs_3 = np.concatenate([zs_1, zs_2[::-1]])

    rs = [rs_1, rs_2, rs_3]
    zs = [zs_1, zs_2, zs_3]

    return FirstWallGeometry(rs=rs, zs=zs)


def first_wall_geometry(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
) -> FirstWallGeometry:
    """Calculates radial and vertical distances for the geometry of section of first wall below the midplane

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
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_blnkith,
        c_fwoth,
        fwith,
        fwoth,
        divgap + blnktth,
        -tfwvt,
        triang,
    )
    rs_1 = np.concatenate([rs1, rs2[::-1]])
    zs_1 = np.concatenate([zs1, zs2[::-1]])
    rs_2 = np.concatenate([rs3, rs4[::-1]])
    zs_2 = -np.concatenate([zs3, zs4[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return FirstWallGeometry(rs=rs, zs=zs)


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
    """Calculates radial and vertical distances for the geometry of section of first wall below the midplane

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
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_blnkith,
        c_fwoth,
        fwith,
        fwoth,
        -(divgap + blnktth),
        -tfwvt,
        triang,
    )
    rs_1 = np.concatenate([rs1, rs2[::-1]])
    zs_1 = np.concatenate([zs1, zs2[::-1]])
    rs_2 = np.concatenate([rs3, rs4[::-1]])
    zs_2 = -np.concatenate([zs3, zs4[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return FirstWallGeometry(rs=rs, zs=zs)
