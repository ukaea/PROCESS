"""
Calculate radial and vertical coordinates for the geometry of the blanket
"""

import numpy as np

from process.geometry.geometry_parameterisations import ArbitraryGeometry
from process.geometry.utils import dh_vertices, dhgap_vertices


def blanket_geometry_single_null(
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
    cumulative_upper: dict,
    triang: float,
    cumulative_lower: dict,
    dz_blkt_upper: float,
    c_shldith: float,
    c_blnkoth: float,
    dr_blkt_inboard: float,
    dr_blkt_outboard: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of the blanket in a single null configuration

    Parameters
    ----------
    radx_outer:
        outboard radius of outer surface of blanket
    rminx_outer:
        inboard radius of outer surface of blanket
    radx_inner:
        outboard radius of inner surface of blanket
    rminx_inner:
        inboard radius of inner surface of blanket
    cumulative_upper:
        cumulative vertical thicknesses of components above the midplane
    triang:
        plasma triangularity
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    dz_blkt_upper:
        top blanket vertical thickness
    c_shldith:
        inboard shield thickness
    c_blnkoth:
        outboard blanket radial thickness
    dr_blkt_inboard:
        inboard blanket radial thickness
    dr_blkt_outboard:
        outboard blanket radial thickness

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Upper blanket outer surface
    kapx = cumulative_upper["dz_blkt_upper"] / rminx_outer
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    # Upper blanket inner surface
    kapx = cumulative_upper["dz_fw_upper"] / rminx_inner
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    # Lower blanket
    divgap = cumulative_lower["dz_divertor"]

    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = blanket_geometry_lower(
        triang=triang,
        dz_blkt_upper=dz_blkt_upper,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        dr_blkt_inboard=dr_blkt_inboard,
        dr_blkt_outboard=dr_blkt_outboard,
        divgap=divgap,
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


def blanket_geometry_lower(
    triang: float,
    dz_blkt_upper: float,
    c_shldith: float,
    c_blnkoth: float,
    dr_blkt_inboard: float,
    dr_blkt_outboard: float,
    divgap: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates radial and vertical distances for the geometry of section of blanket below the midplane

    Parameters
    ----------
    triang:
        plasma triangularity
    dz_blkt_upper:
        top blanket vertical thickness
    c_shldith:
        inboard shield thickness
    c_blnkoth:
        outboard blanket radial thickness
    dr_blkt_inboard:
        inboard blanket radial thickness
    dr_blkt_outboard:
        outboard blanket radial thickness
    divgap:
        divertor structure vertical thickness

    Returns
    -------
    :
        tuple containing the R coordinates for the outboard, Z coordinates for the outboard, R coordinates for the inboard, Z coordinates for the inboard of the blanket geometry below the midplane
    """
    # Lower blanket
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = dhgap_vertices(
        c_shldith,
        c_blnkoth,
        dr_blkt_inboard,
        dr_blkt_outboard,
        divgap,
        -dz_blkt_upper,
        triang,
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
    dz_blkt_upper: float,
    c_shldith: float,
    c_blnkoth: float,
    dr_blkt_inboard: float,
    dr_blkt_outboard: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of blanket in a double null configuration
    In a double null configuration, the geometry of the lower blanket is reflected across the midplane to create the section of blanket above the midplane

    Parameters
    ----------
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    triang:
        plasma triangularity
    dz_blkt_upper:
        top blanket vertical thickness
    c_shldith:
        inboard shield thickness
    c_blnkoth:
        outboard blanket radial thickness
    dr_blkt_inboard:
        inboard blanket radial thickness
    dr_blkt_outboard:
        outboard blanket radial thickness

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Lower blanket
    divgap = cumulative_lower["dz_divertor"]

    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = blanket_geometry_lower(
        triang=triang,
        dz_blkt_upper=dz_blkt_upper,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        dr_blkt_inboard=dr_blkt_inboard,
        dr_blkt_outboard=dr_blkt_outboard,
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
        dz_blkt_upper=dz_blkt_upper,
        c_shldith=c_shldith,
        c_blnkoth=c_blnkoth,
        dr_blkt_inboard=dr_blkt_inboard,
        dr_blkt_outboard=dr_blkt_outboard,
        divgap=divgap,
    )
    rs_1 = np.concatenate([rs_upper_outboard, rs_lower_outboard[::-1]])
    rs_2 = np.concatenate([rs_upper_inboard, rs_lower_inboard[::-1]])
    zs_1 = np.concatenate([zs_upper_outboard, zs_lower_outboard[::-1]])
    zs_2 = np.concatenate([zs_upper_inboard, zs_lower_inboard[::-1]])
    rs = [rs_1, rs_2]
    zs = [zs_1, zs_2]

    return ArbitraryGeometry(rs=rs, zs=zs)
