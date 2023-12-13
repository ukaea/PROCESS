"""
Calculate cryostat geometries
"""
from typing import List

from process.geometry.geometry_parameterisations import RectangleGeometry


def cryostat_geometry(
    rdewex: float, ddwex: float, zdewex: float
) -> List[RectangleGeometry]:
    """Calculates rectangular geometries of the cryostat

    :param rdewex: cryostat internal radius
    :type rdewex: float
    :param ddwex: external cryostat thickness
    :type ddwex: float
    :param zdewex: cryostat internal half-height
    :type zdewex: float
    :return: list of RectangleGeometry - dataclass returning rectangular geometry parameters
    :rtype: List[RectangleGeometry]
    """
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
