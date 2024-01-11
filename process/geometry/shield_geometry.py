"""
Calculate radial and vertical coordinates for the geometry of the shield
"""
from typing import Tuple
import numpy as np
from process.geometry.utils import dh_vertices
from process.geometry.geometry_parameterisations import ArbitraryGeometry


def shield_geometry_single_null(
    cumulative_upper: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
    cumulative_lower: dict,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of shield in a single null configuration

    :param cumulative_upper: cumulative vertical thicknesses of components above the midplane
    :type cumulative_upper: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Upper shield
    # Side furthest from plasma
    kapx = cumulative_upper["shldtth"] / rminx_far
    rs_upper_outboard, zs_upper_outboard = dh_vertices(
        radx_far, rminx_far, triang, kapx
    )

    # Side nearest to plasma
    kapx = (cumulative_upper["vvblgap"]) / rminx_near
    rs_upper_inboard, zs_upper_inboard = dh_vertices(
        radx_near, rminx_near, triang, kapx
    )

    # Lower shield
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = shield_geometry_lower(
        cumulative_lower=cumulative_lower,
        radx_far=radx_far,
        radx_near=radx_near,
        rminx_far=rminx_far,
        rminx_near=rminx_near,
        triang=triang,
    )

    rs = np.concatenate(
        [
            rs_lower_inboard,
            rs_lower_outboard[::-1],
            rs_upper_outboard,
            rs_upper_inboard[::-1],
        ]
    )
    zs = np.concatenate(
        [
            zs_lower_inboard,
            zs_lower_outboard[::-1],
            zs_upper_outboard,
            zs_upper_inboard[::-1],
        ]
    )
    return ArbitraryGeometry(
        rs=rs,
        zs=zs,
    )


def shield_geometry_lower(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates radial and vertical distances for the geometry of section of shield below the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :return: tuple containing the R coordinates for the outboard, Z coordinates for the outboard, R coordinates for the inboard, Z coordinates for the inboard of the shield geometry below the midplane
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    # Side furthest from plasma
    kapx = cumulative_lower["shldlth"] / rminx_far
    rs_lower_outboard, zs_lower_outboard = dh_vertices(
        radx_far, rminx_far, triang, kapx
    )

    # Side nearest to plasma
    kapx = (cumulative_lower["divfix"]) / rminx_near
    rs_lower_inboard, zs_lower_inboard = dh_vertices(
        radx_near, rminx_near, triang, kapx
    )

    return rs_lower_outboard, zs_lower_outboard, rs_lower_inboard, zs_lower_inboard


def shield_geometry_double_null(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> ArbitraryGeometry:
    """Calculates radial and vertical distances for the geometry of shield in a double null configuration
    In a double null configuration, the geometry of the lower shield is reflected across the midplane to create the section of shield above the midplane

    :param cumulative_lower: cumulative vertical thicknesses of components below the midplane
    :type cumulative_lower: dict
    :param radx_far: outboard radius of outer surface of shield
    :type radx_far: float
    :param rminx_far: inboard radius of outer surface of shield
    :type rminx_far: float
    :param radx_near: outboard radius of inner surface of shield
    :type radx_near: float
    :param rminx_near: inboard radius of inner surface of shield
    :type rminx_near: float
    :param triang: plasma triangularity
    :type triang: float
    :return: dataclass returning radial and vertical coordinates
    :rtype: ArbitraryGeometry
    """
    # Lower shield
    (
        rs_lower_outboard,
        zs_lower_outboard,
        rs_lower_inboard,
        zs_lower_inboard,
    ) = shield_geometry_lower(
        cumulative_lower=cumulative_lower,
        radx_far=radx_far,
        radx_near=radx_near,
        rminx_far=rminx_far,
        rminx_near=rminx_near,
        triang=triang,
    )

    rs_lower = np.concatenate([rs_lower_outboard, rs_lower_inboard[::-1]])
    zs_lower = np.concatenate([zs_lower_outboard, zs_lower_inboard[::-1]])

    # Upper shield
    rs_upper = rs_lower
    zs_upper = -1 * zs_lower

    rs = np.concatenate([rs_lower[::-1], rs_upper])
    zs = np.concatenate([zs_lower[::-1], zs_upper])

    return ArbitraryGeometry(rs=rs, zs=zs)
