"""
Calculate radial and vertical coordinates for the geometry of the pf coils and central coil
"""

import numpy as np

from process.models.geometry.parameterisations import RectangleGeometry


def pfcoil_geometry(
    coils_r: list[float],
    coils_z: list[float],
    coils_dr: list[float],
    coils_dz: list[float],
    dr_bore: float,
    dr_cs: float,
    ohdz: float,
) -> tuple[np.ndarray, np.ndarray, RectangleGeometry]:
    """Calculates radial and vertical distances for the geometry of the pf coils and central coil

    Parameters
    ----------
    coils_r:
        list of pf coil radii
    coils_z:
        list of pf coil vertical positions
    coils_dr:
        list of pf coil radial thicknesses
    coils_dz:
        list of pf coil vertical thicknesses
    dr_bore:
        central solenoid inboard radius
    dr_cs:
        central solenoid thickness
    ohdz:
        central solenoid vertical thickness

    Returns
    -------
    :
        tuple containing radial and vertical coordinates for pf coils, and dataclass returning coordinates representing a rectangular geometry used to plot the central coil
    """
    r_points = []
    z_points = []
    for i in range(len(coils_r)):
        r_1 = float(coils_r[i]) - 0.5 * float(coils_dr[i])
        z_1 = float(coils_z[i]) - 0.5 * float(coils_dz[i])
        z_2 = float(coils_z[i]) + 0.5 * float(coils_dz[i])
        r_2 = float(coils_r[i]) + 0.5 * float(coils_dr[i])
        r_points.append([r_1, r_1, r_2, r_2, r_1])
        z_points.append([z_1, z_2, z_2, z_1, z_1])

    central_coil = RectangleGeometry(
        anchor_x=dr_bore, anchor_z=(-ohdz / 2), width=dr_cs, height=ohdz
    )

    return r_points, z_points, central_coil
