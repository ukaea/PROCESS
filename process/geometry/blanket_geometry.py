from dataclasses import dataclass
from typing import List
import numpy as np
from process.geometry.utils import plotdh, plotdhgap


# TODO documentation
@dataclass
class BlanketGeometrySingleNull:
    rs: List[float]
    zs: List[float]


def blanket_geometry_single_null(
    radx_outer: float,
    rminx_outer: float,
    radx_inner: float,
    rminx_inner: float,
    cumulative_upper: dict,
    triang: float,
) -> BlanketGeometrySingleNull:
    point_array = []
    # Upper blanket
    # upper blanket outer surface
    kapx = cumulative_upper["blnktth"] / rminx_outer
    (rs_1, zs_1) = plotdh(radx_outer, rminx_outer, triang, kapx)
    # axis.plot(rs_1, zs_1, color="black", lw=thin)
    point_array.append(rs_1)
    point_array.append(zs_1)

    # upper blanket inner surface
    kapx = cumulative_upper["fwtth"] / rminx_inner
    (rs_2, zs_2) = plotdh(radx_inner, rminx_inner, triang, kapx)
    point_array.append(rs_2)
    point_array.append(zs_2)

    # Plot upper blanket
    rs_3 = np.concatenate([point_array[0], point_array[2][::-1]])
    zs_3 = np.concatenate([point_array[1], point_array[3][::-1]])

    rs = [rs_1, rs_2, rs_3]
    zs = [zs_1, zs_2, zs_3]

    return BlanketGeometrySingleNull(
        rs=rs,
        zs=zs,
    )


@dataclass
class BlanketGeometry:
    rs1: List[float]
    rs2: List[float]
    rs3: List[float]
    rs4: List[float]
    zs1: List[float]
    zs2: List[float]
    zs3: List[float]
    zs4: List[float]


def blanket_geometry(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
) -> BlanketGeometry:
    # Lower blanket
    divgap = cumulative_lower["divfix"]
    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_shldith, c_blnkoth, blnkith, blnkoth, divgap, -blnktth, triang
    )

    return BlanketGeometry(
        rs1=rs1,
        rs2=rs2,
        rs3=rs3,
        rs4=rs4,
        zs1=zs1,
        zs2=zs2,
        zs3=zs3,
        zs4=zs4,
    )


@dataclass
class BlanketGeometryDoubleNull:
    rs1: list
    rs2: list
    rs3: list
    rs4: list
    zs1: list
    zs2: list
    zs3: list
    zs4: list


# i_single_null == 0 part
def blanket_geometry_double_null(
    cumulative_lower: dict,
    triang: float,
    blnktth: float,
    c_shldith: float,
    c_blnkoth: float,
    blnkith: float,
    blnkoth: float,
) -> BlanketGeometryDoubleNull:
    divgap = cumulative_lower["divfix"]

    rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4 = plotdhgap(
        c_shldith,
        c_blnkoth,
        blnkith,
        blnkoth,
        -divgap,
        -blnktth,
        triang,
    )
    return BlanketGeometryDoubleNull(
        rs1=rs1,
        rs2=rs2,
        rs3=rs3,
        rs4=rs4,
        zs1=zs1,
        zs2=zs2,
        zs3=zs3,
        zs4=zs4,
    )
