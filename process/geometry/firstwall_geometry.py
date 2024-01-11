"""
Calculate radial and vertical coordinates for the geometry of the first wall
"""
from typing import Tuple
import numpy as np
from process.geometry.utils import dh_vertices, dhgap_vertices
from process.geometry.geometry_parameterisations import ArbitraryGeometry


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
) -> ArbitraryGeometry:
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
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Upper first wall: outer surface
    kapx = cumulative_upper["fwtth"] / rminx_outer
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    # Upper first wall: inner surface
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    # Lower first wall
    divgap = cumulative_lower["divfix"]
    top_point = divgap + blnktth
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        fwith=fwith,
        fwoth=fwoth,
        tfwvt=tfwvt,
        top_point=top_point,
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


def first_wall_geometry_lower(
    triang: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
    top_point: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    :return: tuple containing the R coordinates for the outboard, Z coordinates for the outboard, R coordinates for the inboard, Z coordinates for the inboard of the first wall geometry below the midplane
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    # Lower first wall
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = dhgap_vertices(
        c_blnkith,
        c_fwoth,
        fwith,
        fwoth,
        top_point,
        -tfwvt,
        triang,
    )
    # outboard radial and vertical coordinates
    rs_lower_outboard = np.concatenate([rs1, rs2[::-1]])
    zs_lower_outboard = np.concatenate([zs1, zs2[::-1]])
    # inboard radial and vertical coordinates
    rs_lower_inboard = np.concatenate([rs3, rs4[::-1]])
    zs_lower_inboard = -np.concatenate([zs3, zs4[::-1]])

    return rs_lower_outboard, zs_lower_outboard, rs_lower_inboard, zs_lower_inboard


def first_wall_geometry_double_null(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_blnkith: float,
    c_fwoth: float,
    fwith: float,
    fwoth: float,
    tfwvt: float,
) -> ArbitraryGeometry:
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
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Lower first wall
    divgap = cumulative_lower["divfix"]
    top_point = divgap + blnktth
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = first_wall_geometry_lower(
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
    (
        rs_upper_outboard,
        zs_upper_outboard,
        rs_upper_inboard,
        zs_upper_inboard,
    ) = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        fwith=fwith,
        fwoth=fwoth,
        tfwvt=tfwvt,
        top_point=top_point,
    )

    rs_1 = np.concatenate([rs_upper_outboard, rs_lower_outboard[::-1]])
    rs_2 = np.concatenate([rs_upper_inboard, rs_lower_inboard[::-1]])
    zs_1 = np.concatenate([zs_upper_outboard, zs_lower_outboard[::-1]])
    zs_2 = np.concatenate([zs_upper_inboard, zs_lower_inboard[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return ArbitraryGeometry(rs=rs, zs=zs)
