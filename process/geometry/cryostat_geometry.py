from typing import List

from process.geometry.geometry_parameterisations import RectangleGeometry


def cryostat_geometry(
    rdewex: float, ddwex: float, zdewex: float
) -> List[RectangleGeometry]:
    rect1 = RectangleGeometry(
        center_x=rdewex, center_z=0, width=ddwex, height=(zdewex + ddwex)
    )
    rect2 = RectangleGeometry(
        center_x=rdewex, center_z=0, width=ddwex, height=-(zdewex + ddwex)
    )

    rect3 = RectangleGeometry(center_x=0, center_z=zdewex, width=rdewex, height=ddwex)

    rect4 = RectangleGeometry(center_x=0, center_z=-zdewex, width=rdewex, height=-ddwex)
    return_rects = [rect1, rect2, rect3, rect4]
    return return_rects
