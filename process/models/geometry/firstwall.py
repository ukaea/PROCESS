"""
Calculate radial and vertical coordinates for the geometry of the first wall
"""

import numpy as np

from process.models.geometry.parameterisations import ArbitraryGeometry
from process.models.geometry.utils import dh_vertices, dhgap_vertices


def first_wall_geometry_single_null(
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
    cumulative_upper: dict,
    triang: float,
    cumulative_lower: dict,
    dz_blkt_upper: float,
    c_blnkith: float,
    c_fwoth: float,
    dr_fw_inboard: float,
    dr_fw_outboard: float,
    tfwvt: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of first wall in a single null configuration

    Parameters
    ----------
    radx_outer:
        outboard radius of outer surface of first wall
    rminx_outer:
        inboard radius of outer surface of first wall
    radx_inner:
        outboard radius of inner surface of first wall
    rminx_inner:
        inboard radius of inner surface of first wall
    cumulative_upper:
        cumulative vertical thicknesses of components above the midplane
    triang:
        plasma triangularity
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    dz_blkt_upper:
        top blanket vertical thickness
    c_blnkith:
        inboard blanket vertical thickness
    c_fwoth:
        outboard first wall vertical thickness
    dr_fw_inboard:
        inboard first wall radial thickness
    dr_fw_outboard:
        outboard first wall radial thickness
    tfwvt:
        top first wall vertical thickness

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Upper first wall: outer surface
    kapx = cumulative_upper["dz_fw_upper"] / rminx_outer
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    # Upper first wall: inner surface
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    # Lower first wall
    divgap = cumulative_lower["dz_divertor"]
    top_point = divgap + dz_blkt_upper
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        dr_fw_inboard=dr_fw_inboard,
        dr_fw_outboard=dr_fw_outboard,
        tfwvt=tfwvt,
        top_point=top_point,
    )

    rs = np.concatenate([
        rs_upper_outboard,
        rs_lower_inboard,
        rs_upper_inboard[::-1],
        rs_lower_outboard[::-1],
    ])
    zs = np.concatenate([
        zs_upper_outboard,
        zs_lower_inboard,
        zs_upper_inboard[::-1],
        zs_lower_outboard[::-1],
    ])

    return ArbitraryGeometry(rs=rs, zs=zs)


def first_wall_geometry_lower(
    triang: float,
    c_blnkith: float,
    c_fwoth: float,
    dr_fw_inboard: float,
    dr_fw_outboard: float,
    tfwvt: float,
    top_point: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates radial and vertical distances for the geometry of section of first wall below the midplane

    Parameters
    ----------
    triang:
        plasma triangularity
    c_blnkith:
        inboard blanket vertical thickness
    c_fwoth:
        outboard first wall vertical thickness
    dr_fw_inboard:
        inboard first wall radial thickness
    dr_fw_outboard:
        outboard first wall radial thickness
    tfwvt:
        top first wall vertical thickness
    top_point:
        top point for plotdhgap, equal to

    Returns
    -------
    :
        tuple containing the R coordinates for the outboard, Z coordinates for the outboard,
        R coordinates for the inboard, Z coordinates for the inboard of the first wall geometry below the midplane
    """
    # Lower first wall
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = dhgap_vertices(
        c_blnkith,
        c_fwoth,
        dr_fw_inboard,
        dr_fw_outboard,
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
    dz_blkt_upper: float,
    c_blnkith: float,
    c_fwoth: float,
    dr_fw_inboard: float,
    dr_fw_outboard: float,
    tfwvt: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of first wall in a double null configuration
    In a double null configuration, the geometry of the lower first wall is reflected across the midplane to create the section of first wall above the midplane

    Parameters
    ----------
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    triang:
        plasma triangularity
    dz_blkt_upper:
        top blanket vertical thickness
    c_blnkith:
        inboard blanket vertical thickness
    c_fwoth:
        outboard first wall vertical thickness
    dr_fw_inboard:
        inboard first wall radial thickness
    dr_fw_outboard:
        outboard first wall radial thickness
    tfwvt:
        top first wall vertical thickness

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Lower first wall
    divgap = cumulative_lower["dz_divertor"]
    top_point = divgap + dz_blkt_upper
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = first_wall_geometry_lower(
        triang=triang,
        c_blnkith=c_blnkith,
        c_fwoth=c_fwoth,
        dr_fw_inboard=dr_fw_inboard,
        dr_fw_outboard=dr_fw_outboard,
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
        dr_fw_inboard=dr_fw_inboard,
        dr_fw_outboard=dr_fw_outboard,
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
