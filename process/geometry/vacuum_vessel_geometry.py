from typing import List
from dataclasses import dataclass
import numpy as np
from process.geometry.utils import plotdh


@dataclass
class VacuumVesselGeometry:
    rs_1: List[float]
    zs_1: List[float]
    rs_2: List[float]
    zs_2: List[float]
    rs: List[float]
    zs: List[float]


def vacuum_vessel_geometry_single_null(
    cumulative_upper: dict,
    upper: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
) -> VacuumVesselGeometry:
    temp_array = []

    kapx = cumulative_upper["d_vv_top"] / rminx_outer
    (rs_1, zs_1) = plotdh(radx_outer, rminx_outer, triang, kapx)
    temp_array.append(rs_1)
    temp_array.append(zs_1)

    kapx = (
        float(cumulative_upper["d_vv_top"]) - float(upper["d_vv_top"])
    ) / rminx_inner
    (rs_2, zs_2) = plotdh(radx_inner, rminx_inner, triang, kapx)
    temp_array.append(rs_2)
    temp_array.append(zs_2)

    rs = np.concatenate([temp_array[0], temp_array[2][::-1]])
    zs = np.concatenate([temp_array[1], temp_array[3][::-1]])

    return VacuumVesselGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )


def vacuum_vessel_geometry(
    cumulative_lower: dict,
    lower: dict,
    triang: float,
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
):
    temp_array = []

    kapx = cumulative_lower["d_vv_bot"] / rminx_outer
    (rs_1, zs_1) = plotdh(radx_outer, rminx_outer, triang, kapx)
    temp_array.append(rs_1)
    temp_array.append(zs_1)

    kapx = (
        float(cumulative_lower["d_vv_bot"]) + float(lower["d_vv_bot"])
    ) / rminx_inner
    (rs_2, zs_2) = plotdh(radx_inner, rminx_inner, triang, kapx)
    temp_array.append(rs_2)
    temp_array.append(zs_2)

    rs = np.concatenate([temp_array[0], temp_array[2][::-1]])
    zs = np.concatenate([temp_array[1], temp_array[3][::-1]])

    return VacuumVesselGeometry(
        rs_1=rs_1,
        zs_1=zs_1,
        rs_2=rs_2,
        zs_2=zs_2,
        rs=rs,
        zs=zs,
    )
