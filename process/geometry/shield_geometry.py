from dataclasses import dataclass
from typing import List
import numpy as np
from process.geometry.utils import plotdh


@dataclass
class ShieldGeometry:
    rs_1: List[float]
    """radial spherical radius 1"""
    zs_1: List[float]
    """vertical distance 1"""
    rs_2: List[float]
    """radial spherical radius 2"""
    zs_2: List[float]
    """vertical distance 2"""
    rs: List[float]
    """radial spherical radii"""
    zs: List[float]
    """vertical distances"""


def shield_geometry_single_null(
    cumulative_upper: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> ShieldGeometry:
    temp_array = []

    # Side furthese from plasma
    kapx = cumulative_upper["shldtth"] / rminx_far
    (rs_1, zs_1) = plotdh(radx_far, rminx_far, triang, kapx)
    temp_array.append(rs_1)
    temp_array.append(zs_1)

    # Side nearest to plasma
    kapx = (cumulative_upper["vvblgap"]) / rminx_near
    (rs_2, zs_2) = plotdh(radx_near, rminx_near, triang, kapx)
    temp_array.append(rs_2)
    temp_array.append(zs_2)

    # Single null: Draw top half from output
    rs = np.concatenate([temp_array[0], temp_array[2][::-1]])
    zs = np.concatenate([temp_array[1], temp_array[3][::-1]])

    return ShieldGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )


def shield_geometry(
    cumulative_lower: dict,
    radx_far: float,
    rminx_far: float,
    radx_near: float,
    rminx_near: float,
    triang: float,
) -> ShieldGeometry:
    temp_array = []
    # Side furthest from plasma
    kapx = cumulative_lower["shldlth"] / rminx_far
    (rs_1, zs_1) = plotdh(radx_far, rminx_far, triang, kapx)
    temp_array.append(rs_1)
    temp_array.append(zs_1)

    # Side nearest to plasma
    kapx = (cumulative_lower["divfix"]) / rminx_near
    (rs_2, zs_2) = plotdh(radx_near, rminx_near, triang, kapx)
    temp_array.append(rs_2)
    temp_array.append(zs_2)

    # Double null: Reflect bottom half to top
    rs = np.concatenate([temp_array[0], temp_array[2][::-1]])
    zs = np.concatenate([temp_array[1], temp_array[3][::-1]])

    return ShieldGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )
