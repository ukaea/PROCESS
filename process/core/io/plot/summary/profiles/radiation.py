"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.common import box_style
from process.core.io.plot.summary.profiles.misc import (
    interp1d_profile,
    profiles_with_pedestal,
    read_imprad_data,
)
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.physics.profiles import calculate_profile_shell_contributions

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes

    from process.core.io.mfile import MFile


__all__ = [
    "plot_line_brem_loss_function_profile",
    "plot_line_brem_power_density_profile",
    "plot_line_brem_power_profile",
    "plot_plasma_pressure_gradient_profiles",
    "plot_rad_density_contour",
]


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

    rho, nd_electron, temp_electron_kev = profiles_with_pedestal(mfile=mfile, scan=scan)
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
        bbox={"boxstyle": "round", "facecolor": "khaki", "alpha": 0.8},
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


def plot_rad_density_contour(axis: Axes, mfile: MFile, scan: int, impp: str):
    """Plots the contour of line and bremsstrahlung radiation density [MW/m³]
    for a plasma cross-section.

    This function reads impurity and plasma profile data, computes the radiation
    density profile, interpolates it onto a 2D grid, and plots the upper and lower
    half contours on the provided axis.

    Parameters
    ----------
    axis:
        The matplotlib axis object to plot the contours on.
    mfile:
        Data object containing plasma and impurity profile information.
    scan:
        The scan index to extract profile data for plotting.
    impp:
        The impurity data path

    Notes
    -----
    - The function assumes the existence of several global or previously defined
       variables and functions, such as `read_imprad_data`, `interp1d_profile`,
       and plasma pedestal parameters.
    - The plotted contours represent the radiation density in units of [MW/m³]
    - The function adds colorbar, axis labels, title, and core reduction annotation to
       the plot.
    """
    rminor = mfile.get("rminor", scan=scan)
    rmajor = mfile.get("rmajor", scan=scan)
    # Read in the impurity data
    imp_data = read_imprad_data(2, impp)
    # imp data is a 3D array with shape (num_impurities, num_temp_points, (temp, lz,
    # zav))

    # Find the relative number density of each impurity
    imp_frac = np.array([
        mfile.get(f"f_nd_impurity_electrons({i:02d})", scan=scan) for i in range(1, 15)
    ])

    # Initialize the radius
    rho, ne, te = profiles_with_pedestal(mfile, scan)
    # imp_data Te values are in eV, but te from profiles_with_pedestal is in keV
    te_ev = te * 1.0e3

    # Intailise the radiation profile arrays
    pimpden = np.zeros([imp_data.shape[0], te.shape[0]])
    lz = np.zeros([imp_data.shape[0], te.shape[0]])
    prad = np.zeros(te.shape[0])

    # Intailise the impurity radiation profile
    for rho in range(te.shape[0]):
        # imp data is a 3D array with shape (num_impurities, num_temp_points, (temp, lz,
        # zav))
        for impurity in range(imp_data.shape[0]):
            # Check if profile temperature is lower than dataset minimum.
            # If so, use the minimum loss function value
            if te_ev[rho] <= imp_data[impurity][0][0]:
                lz[impurity][rho] = imp_data[impurity][0][1]

            # Check if profile temperature is higher than dataset maximum.
            # If so, use the maximum loss function value
            elif te_ev[rho] >= imp_data[impurity][imp_data.shape[1] - 1][0]:
                lz[impurity][rho] = imp_data[impurity][imp_data.shape[1] - 1][1]
            else:
                # If profile valie is within dataset range, use log-log interpolation to
                # find value for loss function
                log_te_data = np.log([row[0] for row in imp_data[impurity]])
                log_lz_data = np.log([row[1] for row in imp_data[impurity]])
                lz[impurity][rho] = np.exp(
                    np.interp(np.log(te_ev[rho]), log_te_data, log_lz_data)
                )
            # Find the power density for each impurity at each rho
            pimpden[impurity][rho] = (
                imp_frac[impurity] * ne[rho] * ne[rho] * lz[impurity][rho]
            )

        for impurity in range(imp_data.shape[0]):
            prad[rho] += pimpden[impurity][rho] * 1.0e-6

    p_rad_grid, r_grid, z_grid = interp1d_profile(prad, mfile, scan)

    # Plot the upper half contour
    p_rad_upper = axis.contourf(
        r_grid, z_grid, p_rad_grid, levels=25, cmap="turbo", zorder=2
    )
    # Plot the lower half contour (mirror)
    axis.contourf(r_grid, -z_grid, p_rad_grid, levels=25, cmap="turbo", zorder=2)

    axis.figure.colorbar(
        p_rad_upper,
        ax=axis,
        label=r"$P_{\mathrm{rad}}$ $[\mathrm{MW.m}^{-3}]$",
        location="left",
        anchor=(-0.25, 0.5),
    )

    axis.set_xlabel("R [m]")
    axis.set_xlim(rmajor - 1.2 * rminor, rmajor + 1.2 * rminor)
    axis.set_ylim(
        -1.2 * rminor * mfile.get("kappa", scan=scan),
        1.2 * mfile.get("kappa", scan=scan) * rminor,
    )
    axis.set_ylabel("Z [m]")
    axis.set_title("Line & Bremsstrahlung Radiation Density Contours")
    axis.plot(
        rmajor,
        0,
        marker="o",
        color="red",
        markersize=6,
        markeredgecolor="black",
        zorder=100,
    )
    # enable minor ticks and grid for clearer reading
    axis.minorticks_on()
    axis.grid(True, which="major", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)

    axis.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.5, zorder=1)
    props_core_reduce = {
        "boxstyle": "round",
        "facecolor": "khaki",
        "alpha": 0.8,
    }
    axis.text(
        0.02,
        0.02,
        rf"$f_{{\text{{core,reduce}}}}$ =  {1.0}",
        transform=axis.transAxes,
        fontsize=8,
        verticalalignment="bottom",
        bbox=props_core_reduce,
    )
    # make minor ticks visible on all sides and draw ticks inward for compact look
    axis.tick_params(which="both", direction="in", top=True, right=True)


def plot_plasma_pressure_gradient_profiles(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot plasma pressure gradient profiles"""
    # Get the plasma pressure profiles
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    pres_plasma_profile = [
        mfile.get(f"pres_plasma_electron_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_ion = [
        mfile.get(f"pres_plasma_ion_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_total = [
        mfile.get(f"pres_plasma_thermal_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_fuel = [
        mfile.get(f"pres_plasma_fuel_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_kpa = np.array(pres_plasma_profile) / 1000.0
    pres_plasma_profile_ion_kpa = np.array(pres_plasma_profile_ion) / 1000.0
    pres_plasma_profile_fuel_kpa = np.array(pres_plasma_profile_fuel) / 1000.0
    pres_plasma_profile_total_kpa = np.array(pres_plasma_profile_total) / 1000.0

    # Calculate the normalised radius
    rho = np.linspace(0, 1, len(pres_plasma_profile_kpa))

    # Compute gradients using numpy.gradient
    grad_electron = np.gradient(pres_plasma_profile_kpa, rho)
    grad_ion = np.gradient(pres_plasma_profile_ion_kpa, rho)
    grad_total = np.gradient(pres_plasma_profile_total_kpa, rho)
    grad_fuel = np.gradient(pres_plasma_profile_fuel_kpa, rho)

    axis.plot(rho, grad_electron, color="blue", label="Electron")
    axis.plot(rho, grad_ion, color="red", label="Ion")
    axis.plot(rho, grad_total, color="green", label="Total")
    axis.plot(rho, grad_fuel, color="orange", label="Fuel")
    axis.set_xlabel("$\\rho$ [r/a]")
    axis.set_ylabel("$dP/dr$ [kPa / m]")
    axis.minorticks_on()
    axis.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
    axis.set_title("Plasma Thermal Pressure Gradient Profiles")
    axis.grid(True, linestyle="--", alpha=0.5)
    axis.set_xlim(0, 1.025)
    axis.legend()
