from typing import List

from process.geometry.geometry_parameterisations import RectangleGeometry


def pfcoil_geometry(
    coils_r: List[float],
    coils_z: List[float],
    coils_dr: List[float],
    coils_dz: List[float],
    bore: float,
    ohcth: float,
    ohdz: float,
) -> tuple:
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
        center_x=bore, center_z=(-ohdz / 2), width=ohcth, height=ohdz
    )

    return r_points, z_points, central_coil
