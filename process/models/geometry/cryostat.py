"""
Calculate cryostat geometries
"""

from process.models.geometry.parameterisations import RectangleGeometry


def cryostat_geometry(
    r_cryostat_inboard: float, dr_cryostat: float, z_cryostat_half_inside: float
) -> list[RectangleGeometry]:
    """Calculates rectangular geometries of the cryostat

    Parameters
    ----------
    r_cryostat_inboard:
        cryostat internal radius
    dr_cryostat:
        external cryostat thickness
    z_cryostat_half_inside:
        cryostat internal half-height

    Returns
    -------
    :
        list of RectangleGeometry - dataclass returning rectangular geometry parameters
    """
    # The cryostat is represented by 4 rectangular sections as follows:

    # rectangle representing vertical part of cryostat above the midplane
    rect1 = RectangleGeometry(
        anchor_x=r_cryostat_inboard,
        anchor_z=0,
        width=dr_cryostat,
        height=(z_cryostat_half_inside + dr_cryostat),
    )

    # rectangle representing vertical part of cryostat below the midplane
    rect2 = RectangleGeometry(
        anchor_x=r_cryostat_inboard,
        anchor_z=0,
        width=dr_cryostat,
        height=-(z_cryostat_half_inside + dr_cryostat),
    )

    # rectangle representing horizontal part of cryostat above the midplane
    rect3 = RectangleGeometry(
        anchor_x=0,
        anchor_z=z_cryostat_half_inside,
        width=r_cryostat_inboard,
        height=dr_cryostat,
    )

    # rectangle representing horizontal part of cryostat below the midplane
    rect4 = RectangleGeometry(
        anchor_x=0,
        anchor_z=-z_cryostat_half_inside,
        width=r_cryostat_inboard,
        height=-dr_cryostat,
    )
    return [rect1, rect2, rect3, rect4]
