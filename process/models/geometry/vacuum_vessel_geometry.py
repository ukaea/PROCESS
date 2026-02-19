"""
Calculate radial and vertical coordinates for the geometry of the vacuum vessel
"""

import numpy as np

from process.models.geometry.geometry_parameterisations import ArbitraryGeometry
from process.models.geometry.utils import dh_vertices


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
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of the vacuum vessel in a single null configuration

    Parameters
    ----------
    cumulative_upper:
        cumulative vertical thicknesses of components above the midplane
    upper:
        vertical thicknesses of components above the midplane
    triang:
        plasma triangularity
    radx_outer:
        outboard radius of outer surface of vacuum vessel
    rminx_outer:
        inboard radius of outer surface of vacuum vessel
    radx_inner:
        outboard radius of inner surface of vacuum vessel
    rminx_inner:
        inboard radius of inner surface of vacuum vessel
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    lower:
        vertical thicknesses of components below the midplane

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Upper vacuum vessel
    kapx = cumulative_upper["dz_vv_upper"] / rminx_outer
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    kapx = (
        float(cumulative_upper["dz_vv_upper"]) - float(upper["dz_vv_upper"])
    ) / rminx_inner
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    # Lower vacuum vessel
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = vacuum_vessel_geometry_lower(
        cumulative_lower=cumulative_lower,
        lower=lower,
        radx_inner=radx_inner,
        radx_outer=radx_outer,
        rminx_inner=rminx_inner,
        rminx_outer=rminx_outer,
        triang=triang,
    )

    rs = np.concatenate([
        rs_lower_inboard,
        rs_lower_outboard[::-1],
        rs_upper_outboard,
        rs_upper_inboard[::-1],
    ])
    zs = np.concatenate([
        zs_lower_inboard,
        zs_lower_outboard[::-1],
        zs_upper_outboard,
        zs_upper_inboard[::-1],
    ])
    return ArbitraryGeometry(
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates radial and vertical distances for the geometry of section of vacuum vessel below the midplane

    Parameters
    ----------
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    lower:
        vertical thicknesses of components below the midplane
    triang:
        plasma triangularity
    radx_outer:
        outboard radius of outer surface of vacuum vessel
    rminx_outer:
        inboard radius of outer surface of vacuum vessel
    radx_inner:
        outboard radius of inner surface of vacuum vessel
    rminx_inner:
        inboard radius of inner surface of vacuum vessel

    Returns
    -------
    :
        tuple containing the R coordinates for the outboard, Z coordinates for the outboard, R coordinates for the inboard, Z coordinates for the inboard of the vacuum vessel geometry below the midplane
    """
    kapx = cumulative_lower["dz_vv_lower"] / rminx_outer
    rs_lower_outboard, zs_lower_outboard = dh_vertices(
        radx_outer, rminx_outer, triang, kapx
    )

    kapx = (
        float(cumulative_lower["dz_vv_lower"]) + float(lower["dz_vv_lower"])
    ) / rminx_inner
    rs_lower_inboard, zs_lower_inboard = dh_vertices(
        radx_inner, rminx_inner, triang, kapx
    )

    return rs_lower_outboard, zs_lower_outboard, rs_lower_inboard, zs_lower_inboard


def vacuum_vessel_geometry_double_null(
    cumulative_lower: dict,
    lower: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of vacuum vessel in a double null configuration
    In a double null configuration, the geometry of the lower vacuum vessel is reflected across the midplane to create the section of vacuum vessel above the midplane

    Parameters
    ----------
    cumulative_lower:
        cumulative vertical thicknesses of components below the midplane
    lower:
        vertical thicknesses of components below the midplane
    triang:
        plasma triangularity
    radx_outer:
        outboard radius of outer surface of vacuum vessel
    rminx_outer:
        inboard radius of outer surface of vacuum vessel
    radx_inner:
        outboard radius of inner surface of vacuum vessel
    rminx_inner:
        inboard radius of inner surface of vacuum vessel

    Returns
    -------
    ArbitraryGeometry
        dataclass returning radial and vertical coordinates
    """
    # Lower vacuum vessel
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = vacuum_vessel_geometry_lower(
        cumulative_lower=cumulative_lower,
        lower=lower,
        radx_inner=radx_inner,
        radx_outer=radx_outer,
        rminx_inner=rminx_inner,
        rminx_outer=rminx_outer,
        triang=triang,
    )
    rs_lower = np.concatenate([rs_lower_outboard, rs_lower_inboard[::-1]])
    zs_lower = np.concatenate([zs_lower_outboard, zs_lower_inboard[::-1]])

    # Upper vacuum vessel
    rs_upper = rs_lower
    zs_upper = -1 * zs_lower

    rs = np.concatenate([rs_lower[::-1], rs_upper])
    zs = np.concatenate([zs_lower[::-1], zs_upper])

    return ArbitraryGeometry(rs=rs, zs=zs)
