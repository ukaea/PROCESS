from typing import List, Tuple
from process.geometry.geometry_parameterisations import RectangleGeometry
from process.geometry.utils import ellips_fill


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
    return_rects = []
    return_rects.append(
        RectangleGeometry(
            center_x=x5 - tfcth,
            center_z=y5 - tfcth,
            width=tfcth,
            height=(y1 - y5 + 2.0 * tfcth),
        )
    )
    return_rects.append(
        RectangleGeometry(
            center_x=x4,
            center_z=y4 - tfcth,
            width=tfcth,
            height=(y2 - y4 + 2.0 * tfcth),
        )
    )
    return_rects.append(
        RectangleGeometry(center_x=x5, center_z=y5 - tfcth, width=x4 - x5, height=tfcth)
    )
    return_rects.append(
        RectangleGeometry(center_x=x1, center_z=y1, width=(x2 - x1), height=tfcth)
    )

    return return_rects


# TODO ask Jon about fn name
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
) -> Tuple[List[RectangleGeometry], list]:
    return_rects = []
    return_verts = []
    # Inboard upper arc
    x0 = x2
    y0 = y1
    a1 = x2 - x1
    b1 = y2 - y1
    a2 = a1 + tfcth
    b2 = b1 + tfcth
    verts1 = ellips_fill(
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
    verts2 = ellips_fill(
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
    verts3 = ellips_fill(
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
    verts4 = ellips_fill(
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
    # Bottom left corner
    return_rects.append(
        RectangleGeometry(
            center_x=x5 - tfcth, center_z=y5, width=tfcth, height=(y1 - y5)
        )
    )
    return_verts = [verts1, verts2, verts3, verts4]

    return return_rects, return_verts
