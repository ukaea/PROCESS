"""
Calculate radial and vertical coordinates for the geometry of the blanket
"""
from typing import Tuple
import numpy as np
from process.geometry.utils import dh_vertices, dhgap_vertices
from process.geometry.geometry_parameterisations import ArbitraryGeometry


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
) -> ArbitraryGeometry:
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
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Upper blanket outer surface
    kapx = cumulative_upper["blnktth"] / rminx_outer
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    # Upper blanket inner surface
    kapx = cumulative_upper["fwtth"] / rminx_inner
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    # Lower blanket
    divgap = cumulative_lower["divfix"]

    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = blanket_geometry_lower(
        triang=triang,
        blnktth=blnktth,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        blnkith=blnkith,
        blnkoth=blnkoth,
        divgap=divgap,
    )

    rs = np.concatenate(
        [
            rs_upper_outboard,
            rs_lower_inboard,
            rs_upper_inboard[::-1],
            rs_lower_outboard[::-1],
        ]
    )
    zs = np.concatenate(
        [
            zs_upper_outboard,
            zs_lower_inboard,
            zs_upper_inboard[::-1],
            zs_lower_outboard[::-1],
        ]
    )

    return ArbitraryGeometry(rs=rs, zs=zs)


def blanket_geometry_lower(
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
    divgap: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    :return: tuple containing the R coordinates for the outboard, Z coordinates for the outboard, R coordinates for the inboard, Z coordinates for the inboard of the blanket geometry below the midplane
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    # Lower blanket
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = dhgap_vertices(
        c_shldith, c_blnkoth, blnkith, blnkoth, divgap, -blnktth, triang
    )
    # outboard radial and vertical coordinates
    rs_lower_outboard = np.concatenate([rs1, rs2[::-1]])
    zs_lower_outboard = np.concatenate([zs1, zs2[::-1]])
    # inboard radial and vertical coordinates
    rs_lower_inboard = np.concatenate([rs3, rs4[::-1]])
    zs_lower_inboard = -np.concatenate([zs3, zs4[::-1]])

    return rs_lower_outboard, zs_lower_outboard, rs_lower_inboard, zs_lower_inboard


def blanket_geometry_double_null(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
) -> ArbitraryGeometry:
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
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Lower blanket
    divgap = cumulative_lower["divfix"]

    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = blanket_geometry_lower(
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

    (
        rs_upper_outboard,
        zs_upper_outboard,
        rs_upper_inboard,
        zs_upper_inboard,
    ) = blanket_geometry_lower(
        triang=triang,
        blnktth=blnktth,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        blnkith=blnkith,
        blnkoth=blnkoth,
        divgap=divgap,
    )
    rs_1 = np.concatenate([rs_upper_outboard, rs_lower_outboard[::-1]])
    rs_2 = np.concatenate([rs_upper_inboard, rs_lower_inboard[::-1]])
    zs_1 = np.concatenate([zs_upper_outboard, zs_lower_outboard[::-1]])
    zs_2 = np.concatenate([zs_upper_inboard, zs_lower_inboard[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return ArbitraryGeometry(rs=rs, zs=zs)
