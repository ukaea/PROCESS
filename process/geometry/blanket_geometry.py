"""
Calculate radial and vertical coordinates for the geometry of the blanket
"""
from dataclasses import dataclass
from typing import List
import numpy as np
from process.geometry.utils import plotdh, plotdhgap


@dataclass
class BlanketGeometry:
    """Holds radial and vertical coordinates for the geometry of a blanket"""

    rs: List[List[float]]
    """radial blanket coordinates"""
    zs: List[List[float]]
    """vertical blanket coordinates"""


def blanket_geometry_single_null(
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
    cumulative_upper: dict,
    triang: float,
    cumulative_lower: dict,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
) -> BlanketGeometry:
    """Calculates radial and vertical distances for the geometry of blanket in a single null configuration

    :param radx_outer: outboard radius of outer surface of blanket
    :type radx_outer: float
    :param rminx_outer: inboard radius of outer surface of blanket
    :type rminx_outer: float
    :param radx_inner: outboard radius of inner surface of blanket
    :type radx_inner: float
    :param rminx_inner: inboard radius of inner surface of blanket
    :type rminx_inner: float
    :param cumulative_upper: cumulative vertical thicknesses of components above the midplane
    :type cumulative_upper: dict
    :param triang: plasma triangularity
    :type triang: float
    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param blnktth: top blanket vertical thickness
    :type blnktth: float
    :param c_shldith: inboard shield thickness
    :type c_shldith: float
    :param c_blnkoth: outboard blanket radial thickness
    :type c_blnkoth: float
    :param blnkith: inboard blanket radial thickness
    :type blnkith: float
    :param blnkoth: outboard blanket radial thickness
    :type blnkoth: float
    :return: BlanketGeometry - dataclass returning blanket radial and vertical coordinates
    :rtype: DataClass
    """
    # Upper blanket outer surface
    kapx = cumulative_upper["blnktth"] / rminx_outer
    rs_upper_1, zs_upper_1 = plotdh(radx_outer, rminx_outer, triang, kapx)

    # Upper blanket inner surface
    kapx = cumulative_upper["fwtth"] / rminx_inner
    rs_upper_2, zs_upper_2 = plotdh(radx_inner, rminx_inner, triang, kapx)

    # Lower blanket
    divgap = cumulative_lower["divfix"]

    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = blanket_geometry_lower(
        triang=triang,
        blnktth=blnktth,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        blnkith=blnkith,
        blnkoth=blnkoth,
        divgap=divgap,
    )

    rs = np.concatenate([rs_upper_1, rs_lower_2, rs_upper_2[::-1], rs_lower_1[::-1]])
    zs = np.concatenate([zs_upper_1, zs_lower_2, zs_upper_2[::-1], zs_lower_1[::-1]])

    return BlanketGeometry(rs=rs, zs=zs)


def blanket_geometry_lower(
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
    divgap: float,
) -> tuple:
    """Calculates radial and vertical distances for the geometry of section of blanket below the midplane

    :param triang: plasma triangularity
    :type triang: float
    :param blnktth: top blanket vertical thickness
    :type blnktth: float
    :param c_shldith: inboard shield thickness
    :type c_shldith: float
    :param c_blnkoth: outboard blanket radial thickness
    :type c_blnkoth: float
    :param blnkith: inboard blanket radial thickness
    :type blnkith: float
    :param blnkoth: outboard blanket radial thickness
    :type blnkoth: float
    :param divgap: divertor structure vertical thickness
    :type divgap: float
    :return: radial and vertical coordinates for blanket geometry below the midplane
    :rtype: tuple
    """
    # Lower blanket
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_shldith, c_blnkoth, blnkith, blnkoth, divgap, -blnktth, triang
    )
    rs_lower_1 = np.concatenate([rs1, rs2[::-1]])
    zs_lower_1 = np.concatenate([zs1, zs2[::-1]])
    rs_lower_2 = np.concatenate([rs3, rs4[::-1]])
    zs_lower_2 = -np.concatenate([zs3, zs4[::-1]])

    return rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2


def blanket_geometry_double_null(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
) -> BlanketGeometry:
    """Calculates radial and vertical distances for the geometry of blanket in a double null configuration
    In a double null configuration, the geometry of the lower blanket is reflected across the midplane to create the section of blanket above the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param triang: plasma triangularity
    :type triang: float
    :param blnktth: top blanket vertical thickness
    :type blnktth: float
    :param c_shldith: inboard shield thickness
    :type c_shldith: float
    :param c_blnkoth: outboard blanket radial thickness
    :type c_blnkoth: float
    :param blnkith: inboard blanket radial thickness
    :type blnkith: float
    :param blnkoth: outboard blanket radial thickness
    :type blnkoth: float
    :return: BlanketGeometry - dataclass returning blanket radial and vertical coordinates
    :rtype: DataClass
    """
    # Lower blanket
    divgap = cumulative_lower["divfix"]

    rs_lower_1, zs_lower_1, rs_lower_2, zs_lower_2 = blanket_geometry_lower(
        triang=triang,
        blnktth=blnktth,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        blnkith=blnkith,
        blnkoth=blnkoth,
        divgap=divgap,
    )

    # Upper blanket
    divgap = -1 * divgap

    rs_upper_1, zs_upper_1, rs_upper_2, zs_upper_2 = blanket_geometry_lower(
        triang=triang,
        blnktth=blnktth,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        blnkith=blnkith,
        blnkoth=blnkoth,
        divgap=divgap,
    )
    rs_1 = np.concatenate([rs_upper_1, rs_lower_1[::-1]])
    rs_2 = np.concatenate([rs_upper_2, rs_lower_2[::-1]])
    zs_1 = np.concatenate([zs_upper_1, zs_lower_1[::-1]])
    zs_2 = np.concatenate([zs_upper_2, zs_lower_2[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return BlanketGeometry(rs=rs, zs=zs)
