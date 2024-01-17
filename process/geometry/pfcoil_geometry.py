"""
Calculate radial and vertical coordinates for the geometry of the pf coils and central coil
"""
from typing import List, Tuple
import numpy as np
from process.geometry.geometry_parameterisations import RectangleGeometry


def pfcoil_geometry(
    coils_r: List[float],
    coils_z: List[float],
    coils_dr: List[float],
    coils_dz: List[float],
    bore: float,
    ohcth: float,
    ohdz: float,
) -> Tuple[np.ndarray, np.ndarray, RectangleGeometry]:
    """Calculates radial and vertical distances for the geometry of the pf coils and central coil

    :param coils_r: list of pf coil radii
    :type coils_r: List[float]
    :param coils_z: list of pf coil vertical positions
    :type coils_z: List[float]
    :param coils_dr: list of pf coil radial thicknesses
    :type coils_dr: List[float]
    :param coils_dz: list of pf coil vertical thicknesses
    :type coils_dz: List[float]
    :param bore: central solenoid inboard radius
    :type bore: float
    :param ohcth: central solenoid thickness
    :type ohcth: float
    :param ohdz: central solenoid vertical thickness
    :type ohdz: float
    :return: tuple containing radial and vertical coordinates for pf coils, and dataclass returning coordinates representing a rectangular geometry used to plot the central coil
    :rtype: Tuple[np.ndarray, np.ndarray, RectangleGeometry]
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
        anchor_x=bore, anchor_z=(-ohdz / 2), width=ohcth, height=ohdz
    )

    return r_points, z_points, central_coil
