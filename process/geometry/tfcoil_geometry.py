"""
Calculate radial and vertical coordinates for the geometry of the tf coils
"""
from typing import List, Tuple
from process.geometry.geometry_parameterisations import RectangleGeometry
from process.geometry.utils import ellips_fill_vertices


def tfcoil_geometry_rectangular_shape(
    x1: float,
    x2: float,
    x4: float,
    x5: float,
    y1: float,
    y2: float,
    y4: float,
    y5: float,
    tfcth: float,
) -> List[RectangleGeometry]:
    """Calculates rectangular geometries for tf coils in a picture frame/rectangular shape parametrization

    :param x1: radial location of arc point 1
    :type x1: float
    :param x2: radial location of arc point 2
    :type x2: float
    :param x4: radial location of arc point 4
    :type x4: float
    :param x5: radial location of arc point 5
    :type x5: float
    :param y1: vertical location of arc point 1
    :type y1: float
    :param y2: vertical location of arc point 2
    :type y2: float
    :param y4: vertical location of arc point 4
    :type y4: float
    :param y5: vertical location of arc point 5
    :type y5: float
    :param tfcth: inboard tf coil thickness
    :type tfcth: float
    :return: list of RectangleGeometry - dataclass returning rectangular geometry parameters
    :rtype: List[RectangleGeometry]
    """
    return_rects = []
    # In this geometry, the tf coil is represented by 4 rectangular sections as follows:

    # rectangle representing the inboard part of the tf coil
    return_rects.append(
        RectangleGeometry(
            anchor_x=x5 - tfcth,
            anchor_z=y5 - tfcth,
            width=tfcth,
            height=(y1 - y5 + 2.0 * tfcth),
        )
    )
    # rectangle representing the outboard part of the tf coil
    return_rects.append(
        RectangleGeometry(
            anchor_x=x4,
            anchor_z=y4 - tfcth,
            width=tfcth,
            height=(y2 - y4 + 2.0 * tfcth),
        )
    )
    # rectangle representing the lower horizontal part of the tf coil
    return_rects.append(
        RectangleGeometry(anchor_x=x5, anchor_z=y5 - tfcth, width=x4 - x5, height=tfcth)
    )
    # rectangle representing the upper horizontal part of the tf coil
    return_rects.append(
        RectangleGeometry(anchor_x=x1, anchor_z=y1, width=(x2 - x1), height=tfcth)
    )

    return return_rects


def tfcoil_geometry_d_shape(
    x1: float,
    x2: float,
    x3: float,
    x4: float,
    x5: float,
    y1: float,
    y2: float,
    y4: float,
    y5: float,
    tfcth: float,
    rtangle: float,
    rtangle2: float,
) -> Tuple[List[RectangleGeometry], List[List[Tuple[float, float]]]]:
    """Calculates radial and vertical distances for the geometry of the tf coils in a D-shape parametrization

    :param x1: radial location of arc point 1
    :type x1: float
    :param x2: radial location of arc point 2
    :type x2: float
    :param x3: radial location of arc point 3
    :type x3: float
    :param x4: radial location of arc point 4
    :type x4: float
    :param x5: radial location of arc point 5
    :type x5: float
    :param y1: vertical location of arc point 1
    :type y1: float
    :param y2: vertical location of arc point 2
    :type y2: float
    :param y4: vertical location of arc point 4
    :type y4: float
    :param y5: vertical location of arc point 5
    :type y5: float
    :param tfcth: inboard tf coil thickness
    :type tfcth: float
    :param rtangle: angle used in tf coil parametrization
    :type rtangle: float
    :param rtangle2: angle used in tf coil parametrization
    :type rtangle2: float
    :return: radial and vertical coordinates for tf coils
    :rtype: Tuple[List[RectangleGeometry], List[List[Tuple[float, float]]]]
    """
    return_rects = []
    return_verts = []
    # Inboard upper arc
    x0 = x2
    y0 = y1
    a1 = x2 - x1
    b1 = y2 - y1
    a2 = a1 + tfcth
    b2 = b1 + tfcth
    verts1 = ellips_fill_vertices(
        a1=a1,
        a2=a2,
        b1=b1,
        b2=b2,
        x0=x0,
        y0=y0,
        ang1=rtangle,
        ang2=rtangle2,
    )
    # Outboard upper arc
    x0 = x2
    y0 = 0
    a1 = x3 - x2
    b1 = y2
    a2 = a1 + tfcth
    b2 = b1 + tfcth
    verts2 = ellips_fill_vertices(
        a1=a1,
        a2=a2,
        b1=b1,
        b2=b2,
        x0=x0,
        y0=y0,
        ang1=0,
        ang2=rtangle,
    )
    # Inboard lower arc
    x0 = x4
    y0 = y5
    a1 = x4 - x5
    b1 = y5 - y4
    a2 = a1 + tfcth
    b2 = b1 + tfcth
    verts3 = ellips_fill_vertices(
        a1=a1,
        a2=a2,
        b1=b1,
        b2=b2,
        x0=x0,
        y0=y0,
        ang1=-rtangle,
        ang2=-1 * rtangle2,
    )
    # Outboard lower arc
    x0 = x4
    y0 = 0
    a1 = x3 - x2
    b1 = -y4
    a2 = a1 + tfcth
    b2 = b1 + tfcth
    verts4 = ellips_fill_vertices(
        a1=a1,
        a2=a2,
        b1=b1,
        b2=b2,
        x0=x0,
        y0=y0,
        ang1=0,
        ang2=-rtangle,
    )
    # Vertical leg
    return_rects.append(
        RectangleGeometry(
            anchor_x=x5 - tfcth, anchor_z=y5, width=tfcth, height=(y1 - y5)
        )
    )
    return_verts = [verts1, verts2, verts3, verts4]

    return return_rects, return_verts
