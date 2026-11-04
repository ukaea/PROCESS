"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import interp1d

from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.geometry.plasma import plasma_geometry
from process.models.physics.impurity_radiation import read_impurity_file

if TYPE_CHECKING:
    from process.core.io.mfile import MFile

__all__ = ["interp1d_profile", "profiles_with_pedestal"]


def read_imprad_data(_skiprows, data_path):
    """Function to read all data needed for creation of radiation profile

    Parameters
    ----------
    _skiprows :
        number of rows to skip when reading impurity data files
    data_path :
        path to impurity data

    """
    imp_label = ImpurityRadiationData().imp_label
    lzdata = []

    for label in imp_label:
        file_iden = data_path + label.ljust(3, "_")

        Te = None
        lz = None
        zav = None

        for header in read_impurity_file(file_iden + "lz_tau.dat"):
            if "Te[eV]" in header.content:
                Te = np.asarray(header.data, dtype=float)

            if "infinite confinement" in header.content:
                lz = np.asarray(header.data, dtype=float)
        for header in read_impurity_file(file_iden + "z_tau.dat"):
            if "infinite confinement" in header.content:
                zav = np.asarray(header.data, dtype=float)

        lzdata.append(np.column_stack([Te, lz, zav]))

    # then switch string to floats
    return np.array(lzdata, dtype=float)


def profiles_with_pedestal(mfile, scan: int):
    """Calculate profiles with pedestal"""
    alphan = mfile.get("alphan", scan=scan)
    alphat = mfile.get("alphat", scan=scan)
    nd_plasma_electron_on_axis = mfile.get("nd_plasma_electron_on_axis", scan=scan)
    temp_plasma_electron_on_axis_kev = mfile.get(
        "temp_plasma_electron_on_axis_kev", scan=scan
    )

    radius_plasma_pedestal_temp_norm = mfile.get(
        "radius_plasma_pedestal_temp_norm", scan=scan
    )

    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))
    i_plasma_pedestal = mfile.get("i_plasma_pedestal", scan=scan)
    nd_plasma_pedestal_electron = mfile.get("nd_plasma_pedestal_electron", scan=scan)
    radius_plasma_pedestal_density_norm = mfile.get(
        "radius_plasma_pedestal_density_norm", scan=scan
    )
    ne0 = mfile.get("nd_plasma_electron_on_axis", scan=scan)
    rho = np.linspace(0, 1.0, n_plasma_profile_elements)
    nd_plasma_separatrix_electron = mfile.get("nd_plasma_separatrix_electron", scan=scan)
    temp_plasma_pedestal_electron_kev = mfile.get(
        "temp_plasma_pedestal_electron_kev", scan=scan
    )
    temp_plasma_separatrix_electron_kev = mfile.get(
        "temp_plasma_separatrix_electron_kev", scan=scan
    )
    tbeta = mfile.get("tbeta", scan=scan)
    te0 = mfile.get("temp_plasma_electron_on_axis_kev", scan=scan)

    if i_plasma_pedestal == 0:
        # Initialise the radius

        # The density profile
        ne = nd_plasma_electron_on_axis * (1 - rho**2) ** alphan

        # The temperature profile
        te = temp_plasma_electron_on_axis_kev * (1 - rho**2) ** alphat

    # Profiles with pedestal
    elif i_plasma_pedestal == 1:
        # The density and temperature profile
        # Initiliase empty normalised array with zeros
        ne = np.zeros_like(rho)
        te = np.zeros_like(rho)
        # Reconstruct the temperature and density profiles with pedestal
        for q in range(rho.shape[0]):
            # Core density region
            if rho[q] <= radius_plasma_pedestal_density_norm:
                ne[q] = (
                    nd_plasma_pedestal_electron
                    + (ne0 - nd_plasma_pedestal_electron)
                    * (1 - rho[q] ** 2 / radius_plasma_pedestal_density_norm**2)
                    ** alphan
                )
            else:
                # Pedestal density region
                ne[q] = nd_plasma_separatrix_electron + (
                    nd_plasma_pedestal_electron - nd_plasma_separatrix_electron
                ) * (1 - rho[q]) / (1 - radius_plasma_pedestal_density_norm)

            # Core temperature region
            if rho[q] <= radius_plasma_pedestal_temp_norm:
                te[q] = (
                    temp_plasma_pedestal_electron_kev
                    + (te0 - temp_plasma_pedestal_electron_kev)
                    * (1 - (rho[q] / radius_plasma_pedestal_temp_norm) ** tbeta)
                    ** alphat
                )
            else:
                # Pedestal temperature region
                te[q] = temp_plasma_separatrix_electron_kev + (
                    temp_plasma_pedestal_electron_kev
                    - temp_plasma_separatrix_electron_kev
                ) * (1 - rho[q]) / (1 - radius_plasma_pedestal_temp_norm)

    return rho, ne, te


def interp1d_profile(profile, mfile: MFile, scan: int):
    """Interpolate profile over a grid"""
    # Get plasma geometry and boundary
    pg = plasma_geometry(
        rmajor=mfile.get("rmajor", scan=scan),
        rminor=mfile.get("rminor", scan=scan),
        triang=mfile.get("triang", scan=scan),
        kappa=mfile.get("kappa", scan=scan),
        i_single_null=mfile.get("i_single_null", scan=scan),
        i_plasma_shape=mfile.get("i_plasma_shape", scan=scan),
        square=mfile.get("plasma_square", scan=scan),
    )

    # Create a grid of (R, Z) points inside the plasma boundary
    rho = np.linspace(0, 1, 500)
    theta = np.linspace(0, 2 * np.pi, 720)
    rho_grid, theta_grid = np.meshgrid(rho, theta)

    # Map (rho, theta) to (R, Z) using plasma boundary shape
    # For each theta, get boundary (R, Z), then scale by rho
    bdry_r = pg.rs
    bdry_z = pg.zs
    # Interpolate boundary for all theta
    bdry_theta = np.arctan2(bdry_z - pg.zs.mean(), bdry_r - pg.rs.mean())
    # Ensure bdry_theta is monotonic and covers [0, 2pi]
    bdry_theta = np.unwrap(bdry_theta)
    # Sort bdry_theta and corresponding r/z for monotonic interpolation
    sort_idx = np.argsort(bdry_theta)
    bdry_theta = bdry_theta[sort_idx]
    bdry_r = bdry_r[sort_idx]
    bdry_z = bdry_z[sort_idx]
    # Extend boundary to cover full [0, 2pi] if needed
    if bdry_theta[0] > 0 or bdry_theta[-1] < 2 * np.pi:
        bdry_theta = np.concatenate(([0], bdry_theta, [2 * np.pi]))
        bdry_r = np.concatenate(([bdry_r[0]], bdry_r, [bdry_r[-1]]))
        bdry_z = np.concatenate(([bdry_z[0]], bdry_z, [bdry_z[-1]]))
    # Map theta to boundary r/z
    f_r = interp1d(
        bdry_theta,
        bdry_r,
        kind="linear",
        fill_value="extrapolate",
        assume_sorted=True,
    )
    # Map theta to boundary z
    f_z = interp1d(
        bdry_theta,
        bdry_z,
        kind="linear",
        fill_value="extrapolate",
        assume_sorted=True,
    )
    # For each (theta, rho), get boundary (R, Z), then scale by rho
    # Use the boundary center for scaling, not mean, to avoid vertical offset
    r_center = mfile.get("rmajor", scan=scan)
    z_center = pg.zs.mean()
    r_grid = r_center + (f_r(theta_grid) - r_center) * rho_grid
    z_grid = z_center + (f_z(theta_grid) - z_center) * rho_grid

    # Interpolate profile for each rho
    profile_grid = np.interp(rho_grid, np.linspace(0, 1, len(profile)), profile)

    return profile_grid, r_grid, z_grid
