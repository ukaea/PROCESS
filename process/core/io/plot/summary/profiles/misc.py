"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import interp1d

from process.core.io.plot.summary.common import box_style
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.geometry.plasma import plasma_geometry
from process.models.physics.impurity_radiation import read_impurity_file
from process.models.physics.profiles import calculate_profile_shell_contributions

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


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


def plot_line_brem_power_density_profile(
    axis: plt.Axes, mfile: MFile, scan: int, impp: str, demo_ranges: bool
):
    """Function to plot Line and Bremsstrahlung radiation power density [MW/m³] profile.

    Parameters
    ----------
    axis : plt.Axes
        axis object to add plot to
    mfile : MFile
        MFile object containing plasma and impurity profile information.
    scan : int
        scan number to use
    impp : str
        impurity path
    demo_ranges : bool
        whether to use fixed demo ranges for the plot

    """
    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_ylabel(r"$P_{\mathrm{rad}}$ $[\mathrm{MW.m}^{-3}]$")
    axis.set_title("Raw Data: Line & Bremsstrahlung Radiation Density Profile")

    # read in the impurity data
    imp_data = read_imprad_data(_skiprows=2, data_path=impp)

    # find impurity densities
    imp_frac = np.array([
        mfile.get(f"f_nd_impurity_electrons({i:02d})", scan=scan) for i in range(1, 15)
    ])

    rho, nd_electron, temp_electron_kev = profiles_with_pedestal(mfile, scan)
    # imp_data Te values are in eV, but te from profiles_with_pedestal is in keV
    temp_electron_ev = temp_electron_kev * 1.0e3

    # Intailise the radiation profile arrays
    pden_rad_array = np.zeros([imp_data.shape[0], temp_electron_kev.shape[0]])
    lz = np.zeros([imp_data.shape[0], temp_electron_kev.shape[0]])
    pden_total_profile = np.zeros(temp_electron_kev.shape[0])

    # Intailise the impurity radiation profile
    for temp_point in range(temp_electron_kev.shape[0]):
        for impurity in range(imp_data.shape[0]):
            if temp_electron_ev[temp_point] <= imp_data[impurity][0][0]:
                lz[impurity][temp_point] = imp_data[impurity][0][1]
            elif (
                temp_electron_ev[temp_point]
                >= imp_data[impurity][imp_data.shape[1] - 1][0]
            ):
                lz[impurity][temp_point] = imp_data[impurity][imp_data.shape[1] - 1][1]
            else:
                # Use np.interp for log-log interpolation
                log_te_data = np.log([row[0] for row in imp_data[impurity]])
                log_lz_data = np.log([row[1] for row in imp_data[impurity]])
                lz[impurity][temp_point] = np.exp(
                    np.interp(
                        np.log(temp_electron_ev[temp_point]), log_te_data, log_lz_data
                    )
                )
                pden_rad_array[impurity][temp_point] = (
                    imp_frac[impurity]
                    * nd_electron[temp_point]
                    * nd_electron[temp_point]
                    * lz[impurity][temp_point]
                )

        for l_ in range(imp_data.shape[0]):
            pden_total_profile[temp_point] += pden_rad_array[l_][temp_point] * 1.0e-6

        # Plot the total radiation profile and individual impurity contributions
        axis.plot(rho, pden_total_profile, label="Total", linestyle="dotted")
        axis.plot(rho, pden_rad_array[0] * 1.0e-6, label="H")
        axis.plot(rho, pden_rad_array[1] * 1.0e-6, label="He")

        # Plot the remaining impurity contributions if their fraction is significant
        imp_labels = ImpurityRadiationData().imp_label
        for ind in range(2, imp_data.shape[0]):
            if imp_frac[ind] > 1.0e-30:
                axis.plot(
                    rho,
                    pden_rad_array[ind] * 1.0e-6,
                    label=imp_labels[ind].replace("_", ""),
                )

    axis.minorticks_on()
    # Plot a vertical line at the core region radius
    core_radius = mfile.get("radius_plasma_core_norm", scan=scan)

    # Plot a vertical line at the core region radius
    axis.axvline(x=core_radius, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    # Plot a box in the bottom left with f_{core,reduce} and \rho_{core}
    axis.text(
        0.02,
        0.02,
        rf"$f_{{\text{{core,reduce}}}}$ =  {1.0}\n"
        rf"$\rho_{{\text{{core}}}}$ =  {core_radius:.3f}",
        transform=axis.transAxes,
        fontsize=8,
        verticalalignment="bottom",
        bbox=box_style("khaki", alpha=0.8, linewidth=None),
    )

    # Ranges
    # ---
    axis.legend(loc="upper left", bbox_to_anchor=(-0.1, -0.1), ncol=4)
    axis.set_xlim(0, 1.0)
    axis.set_yscale("log")
    axis.yaxis.grid(True, which="both", alpha=0.2)
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        axis.set_ylim(1e-4, 0.5)

    # Adaptive ranges
    else:
        axis.set_ylim(1e-4, axis.get_ylim()[1])


def plot_line_brem_power_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    impp: str,
):
    """Function to plot Line and Bremsstrahlung radiation power [MW] profile.

    Parameters
    ----------
    axis : plt.Axes
        axis object to add plot to
    mfile : MFile
        MFile object containing plasma and impurity profile information.
    scan : int
        scan number to use
    impp : str
        impurity path

    """
    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_ylabel(r"$P_{\mathrm{rad}}$ $[\mathrm{MW}]$")
    axis.set_title("Line & Bremsstrahlung Radiation Power Profile")

    # read in the impurity data
    imp_data = read_imprad_data(_skiprows=2, data_path=impp)

    # find impurity densities
    imp_frac = np.array([
        mfile.get(f"f_nd_impurity_electrons({i:02d})", scan=scan) for i in range(1, 15)
    ])
    vol_plasma = mfile.get("vol_plasma", scan=scan)
    p_plasma_rad_imps_mw = mfile.get("p_plasma_rad_mw", scan=scan) - mfile.get(
        "p_plasma_sync_mw", scan=scan
    )

    rho, nd_electron, temp_electron_kev = profiles_with_pedestal(mfile=mfile, scan=scan)
    # imp_data Te values are in eV, but te from profiles_with_pedestal is in keV
    temp_electron_ev = temp_electron_kev * 1.0e3

    # Intailise the radiation profile arrays
    p_rad_array = np.zeros([imp_data.shape[0], temp_electron_kev.shape[0]])
    lz = np.zeros([imp_data.shape[0], temp_electron_kev.shape[0]])
    p_total_profile = np.zeros(temp_electron_kev.shape[0])

    # Intailise the impurity radiation profile
    for temp_point in range(temp_electron_kev.shape[0]):
        for impurity in range(imp_data.shape[0]):
            if temp_electron_ev[temp_point] <= imp_data[impurity][0][0]:
                lz[impurity][temp_point] = imp_data[impurity][0][1]
            elif (
                temp_electron_ev[temp_point]
                >= imp_data[impurity][imp_data.shape[1] - 1][0]
            ):
                lz[impurity][temp_point] = imp_data[impurity][imp_data.shape[1] - 1][1]
            else:
                # Use np.interp for log-log interpolation
                log_te_data = np.log([row[0] for row in imp_data[impurity]])
                log_lz_data = np.log([row[1] for row in imp_data[impurity]])
                lz[impurity][temp_point] = np.exp(
                    np.interp(
                        np.log(temp_electron_ev[temp_point]), log_te_data, log_lz_data
                    )
                )

            # Calculate the absolue radiation power for each volumetric shell for each
            # impurity
            p_rad_array[impurity] = calculate_profile_shell_contributions(
                profile_x=rho,
                profile_y=(
                    imp_frac[impurity] * nd_electron * nd_electron * lz[impurity]
                ),
                vol_plasma=vol_plasma,
                profile_dx=rho[1] - rho[0],
            )

        for l_ in range(imp_data.shape[0]):
            p_total_profile[temp_point] += p_rad_array[l_][temp_point] * 1.0e-6

    axis.plot(
        rho, p_total_profile, label="$P_{\\Sigma\\text{Impurities}}$", linestyle="dotted"
    )

    # Plot individual impurity radiation profiles
    axis.plot(rho, p_rad_array[0] * 1.0e-6, label="H")
    axis.plot(rho, p_rad_array[1] * 1.0e-6, label="He")
    # Only plot impurities with a significant fraction
    for ind in range(2, imp_data.shape[0]):
        if imp_frac[ind] > 1.0e-30:
            axis.plot(
                rho,
                p_rad_array[ind] * 1.0e-6,
                label=ImpurityRadiationData().imp_label[ind].replace("_", ""),
            )

    axis.plot(
        rho,
        np.cumsum(p_total_profile),
        label="$\\Sigma P_{\\Sigma\\text{Impurities}}$",
        color="black",
    )
    axis.axhline(
        y=p_plasma_rad_imps_mw,
        color="black",
        linestyle="--",
        label="$P_{\\text{total}}$",
    )

    axis.minorticks_on()
    # Plot a vertical line at the core region radius
    core_radius = mfile.get("radius_plasma_core_norm", scan=scan)

    # Plot a vertical line at the core region radius
    axis.axvline(x=core_radius, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    # Plot a box in the bottom left with f_{core,reduce} and \rho_{core}
    axis.text(
        0.05,
        0.02,
        rf"$f_{{\text{{core,reduce}}}}$ =  {1.0}"
        "\n"
        rf"$\rho_{{\text{{core}}}}$ =  {core_radius:.3f}",
        transform=axis.transAxes,
        fontsize=8,
        verticalalignment="bottom",
        bbox=box_style("khaki", alpha=0.8, linewidth=None),
    )

    cumulative_thermal_energy_mj = np.cumsum(p_total_profile)
    half_thermal_energy_mj = 0.5 * p_plasma_rad_imps_mw
    half_thermal_energy_position = np.interp(
        half_thermal_energy_mj,
        cumulative_thermal_energy_mj,
        np.linspace(0, 1, temp_electron_kev.shape[0]),
    )
    axis.axhline(
        y=half_thermal_energy_mj,
        label="$50\\%\\ P_{\\text{total}}$",
        color="tab:green",
        linestyle=":",
    )
    axis.axvline(
        x=half_thermal_energy_position,
        color="tab:green",
        linestyle=":",
    )

    # Ranges
    # ---
    axis.set_xlim(0, 1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), ncol=1)
    axis.set_yscale("log")
    axis.yaxis.grid(True, which="both", alpha=0.2)


def plot_line_brem_loss_function_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    impp: str,
):
    """Function to plot Line and Bremsstrahlung loss function (L_z) profile.

    Parameters
    ----------
    axis : plt.Axes
        axis object to add plot to
    mfile : MFile
        MFile object containing plasma and impurity profile information.
    scan : int
        scan number to use
    impp : str
        impurity path

    """
    # read in the impurity data
    imp_data = read_imprad_data(_skiprows=2, data_path=impp)

    # find impurity densities
    imp_frac = np.array([
        mfile.get(f"f_nd_impurity_electrons({i:02d})", scan=scan) for i in range(1, 15)
    ])

    rho, _, temp_electron_kev = profiles_with_pedestal(mfile, scan)
    # imp_data Te values are in eV, but te from profiles_with_pedestal is in keV
    temp_electron_ev = temp_electron_kev * 1.0e3

    # Intailise the radiation profile arrays
    lz = np.zeros([imp_data.shape[0], temp_electron_kev.shape[0]])

    # Intailise the impurity radiation profile
    for temp_point in range(temp_electron_kev.shape[0]):
        for impurity in range(imp_data.shape[0]):
            if temp_electron_ev[temp_point] <= imp_data[impurity][0][0]:
                lz[impurity][temp_point] = imp_data[impurity][0][1]
            elif (
                temp_electron_ev[temp_point]
                >= imp_data[impurity][imp_data.shape[1] - 1][0]
            ):
                lz[impurity][temp_point] = imp_data[impurity][imp_data.shape[1] - 1][1]
            else:
                # Use np.interp for log-log interpolation
                log_te_data = np.log([row[0] for row in imp_data[impurity]])
                log_lz_data = np.log([row[1] for row in imp_data[impurity]])
                lz[impurity][temp_point] = np.exp(
                    np.interp(
                        np.log(temp_electron_ev[temp_point]), log_te_data, log_lz_data
                    )
                )

    # Plot the radiation loss function profiles
    axis.plot(rho, lz[0], label="H")
    axis.plot(rho, lz[1], label="He")
    # Plot the remaining impurities if their fraction is significant
    impurity_data = ImpurityRadiationData()
    for ind in range(2, imp_data.shape[0]):
        if imp_frac[ind] > 1.0e-30:
            axis.plot(
                rho,
                lz[ind],
                label=impurity_data.imp_label[ind].replace("_", ""),
            )
    axis.legend(loc="best", ncol=4)
    axis.minorticks_on()

    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_ylabel(r"$L_z$ $[\mathrm{W}\mathrm{m}^3]$")
    axis.set_title("Line & Bremsstrahlung Loss Function ($L_z$) Profiles")
    axis.set_xlim(0, 1.0)
    axis.set_yscale("log")
    axis.yaxis.grid(True, which="both", alpha=0.2)


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


__all__ = [
    "interp1d_profile",
    "plot_line_brem_loss_function_profile",
    "plot_line_brem_power_density_profile",
    "profiles_with_pedestal",
]
