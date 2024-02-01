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
    # The cryostat is represented by 4 rectangular sections as follows:

    # rectangle representing vertical part of cryostat above the midplane
    rect1 = RectangleGeometry(
        anchor_x=rdewex, anchor_z=0, width=ddwex, height=(zdewex + ddwex)
    )

    # rectangle representing vertical part of cryostat below the midplane
    rect2 = RectangleGeometry(
        anchor_x=rdewex, anchor_z=0, width=ddwex, height=-(zdewex + ddwex)
    )

    # rectangle representing horizontal part of cryostat above the midplane
    rect3 = RectangleGeometry(anchor_x=0, anchor_z=zdewex, width=rdewex, height=ddwex)

    # rectangle representing horizontal part of cryostat below the midplane
    rect4 = RectangleGeometry(anchor_x=0, anchor_z=-zdewex, width=rdewex, height=-ddwex)
    return_rects = [rect1, rect2, rect3, rect4]
    return return_rects
