"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.common import (
    add_colourbar,
)
from process.core.io.plot.summary.profiles.misc import (
    interp1d_profile,
    profiles_with_pedestal,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)
from process.models.engineering.materials import (
    poisson_steel,
)
from process.models.pfcoil import N_CS_STRESS_PROFILE_POINTS, CSCoil
from process.models.physics.impurity_radiation import read_impurity_file

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes

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
    label = [
        "H_",
        "He",
        "Be",
        "C_",
        "N_",
        "O_",
        "Ne",
        "Si",
        "Ar",
        "Fe",
        "Ni",
        "Kr",
        "Xe",
        "W_",
    ]
    lzdata = [0.0 for x in range(len(label))]

    for i in range(len(label)):
        file_iden = data_path + label[i].ljust(3, "_")

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

        lzdata[i] = np.column_stack([Te, lz, zav])

    # then switch string to floats
    return np.array(lzdata, dtype=float)


def plot_rad_contour(axis: Axes, mfile: MFile, scan: int, impp: str):
    """Plots the contour of line and bremsstrahlung radiation density for a plasma
    cross-section.

    This function reads impurity and plasma profile data, computes the radiation density
    profile,
    interpolates it onto a 2D grid, and plots the upper and lower half contours on the
    provided axis.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The matplotlib axis object to plot the contours on.
    mfile : Any
        Data object containing plasma and impurity profile information.
    scan : int
        The scan index to extract profile data for plotting.
    impp : str
        The impurity data path

    Notes
    -----
    - The function assumes the existence of several global or previously defined
    variables and functions,
        such as `read_imprad_data`, `interp1d_profile`, and plasma pedestal parameters.
    - The plotted contours represent the radiation density in units of MW.m^-3.
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
            if te[rho] <= imp_data[impurity][0][0]:
                lz[impurity][rho] = imp_data[impurity][0][1]

            # Check if profile temperature is higher than dataset maximum.
            # If so, use the maximum loss function value
            elif te[rho] >= imp_data[impurity][imp_data.shape[1] - 1][0]:
                lz[impurity][rho] = imp_data[impurity][imp_data.shape[1] - 1][1]
            else:
                # If profile valie is within dataset range, use log-log interpolation to
                # find value for loss function
                log_te_data = np.log([row[0] for row in imp_data[impurity]])
                log_lz_data = np.log([row[1] for row in imp_data[impurity]])
                lz[impurity][rho] = np.exp(
                    np.interp(np.log(te[rho]), log_te_data, log_lz_data)
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
        r_grid, z_grid, p_rad_grid, levels=50, cmap="plasma", zorder=2
    )
    # Plot the lower half contour (mirror)
    axis.contourf(r_grid, -z_grid, p_rad_grid, levels=50, cmap="plasma", zorder=2)

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
    draw_text(
        axis,
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


def plot_larmor_radius_profile(axis: plt.Axes, mfile_data: MFile, scan: int):
    """Plot the Larmor radius profile on the given axis."""
    radius_plasma_deuteron_larmor_profile = [
        mfile_data.data[
            f"radius_plasma_deuteron_toroidal_larmor_isotropic_profile{i}"
        ].get_scan(scan)
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    radius_plasma_triton_larmor_profile = [
        mfile_data.data[
            f"radius_plasma_triton_toroidal_larmor_isotropic_profile{i}"
        ].get_scan(scan)
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    radius_plasma_deuteron_larmor_profile_mm = [
        radius * 1e3 for radius in radius_plasma_deuteron_larmor_profile
    ]

    radius_plasma_triton_larmor_profile_mm = [
        radius * 1e3 for radius in radius_plasma_triton_larmor_profile
    ]

    axis.plot(
        np.linspace(-1, 1, len(radius_plasma_deuteron_larmor_profile_mm)),
        radius_plasma_deuteron_larmor_profile_mm,
        color="red",
        linestyle="-",
        label=r"$\rho_{Larmor,toroidal,D}$",
    )

    axis.plot(
        np.linspace(-1, 1, len(radius_plasma_triton_larmor_profile_mm)),
        radius_plasma_triton_larmor_profile_mm,
        color="green",
        linestyle="-",
        label=r"$\rho_{Larmor,toroidal,T}$",
    )

    axis.set_ylabel(r"Larmor Radii [mm]")
    axis.set_title(r" Toroidal Larmor Radii ($v_{\perp}^2 = 2v_{th}^2$)")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_cs_radial_stress_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    j_cs: float,
    b_cs_inner: float,
):
    """Plot CS radial stress profile"""
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)

    radii = np.linspace(r_cs_inner, r_cs_outer, num=25)
    stress_values = np.array([
        CSCoil.calculate_cs_radial_stress(
            r_stress_point=radius,
            r_cs_inner=r_cs_inner,
            r_cs_outer=r_cs_outer,
            j_cs=j_cs,
            b_cs_inner=b_cs_inner,
            f_poisson_cs_structure=poisson_steel,
        )
        for radius in radii
    ])

    axis.plot(
        radii,
        stress_values / 1e6,
        linewidth=2,
        label="$\\sigma_{r}$,Radial Stress",
    )
    max_idx = np.argmax(np.abs(stress_values))
    max_radius = radii[max_idx]
    max_stress = stress_values[max_idx] / 1e6
    axis.axvline(max_radius, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.axhline(max_stress, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.set_xlabel("Radial Position (m)")
    axis.set_ylabel("Radial Stress (MPa)")
    axis.minorticks_on()
    axis.grid(True, alpha=0.3)
    axis.set_title("CS Radial Stress at BOP")


def plot_cs_radial_stress_contour_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    j_cs: float,
    b_cs_inner: float,
    colorbar_axis: plt.Axes | None = None,
):
    """Plot CS radial stress contour profile"""
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)

    # Create 2D grid for contour plot: radial and vertical dimensions
    n_radial = 50
    radial_grid = np.linspace(r_cs_inner, r_cs_outer, n_radial)
    height_grid = np.linspace(
        -dz_cs_full / 2, dz_cs_full / 2, N_CS_STRESS_PROFILE_POINTS
    )

    # Create meshgrid for filled contour
    r, z = np.meshgrid(radial_grid, height_grid)

    # Calculate radial stress across the 2D grid
    stress_data = np.zeros((len(height_grid), n_radial))
    for i in range(len(height_grid)):
        for j in range(n_radial):
            stress_data[i, j] = (
                CSCoil.calculate_cs_radial_stress(
                    r_stress_point=radial_grid[j],
                    r_cs_inner=r_cs_inner,
                    r_cs_outer=r_cs_outer,
                    j_cs=j_cs,
                    b_cs_inner=b_cs_inner,
                    f_poisson_cs_structure=poisson_steel,
                )
                / 1e6
            )

    # Plot filled contour of stress distribution
    contour_fill = axis.contourf(r, z, stress_data, levels=15, cmap="RdYlBu")
    contour_lines = axis.contour(
        r,
        z,
        stress_data,
        levels=[stress_data.max()],
        colors="black",
        linewidths=0.5,
        alpha=0.4,
    )
    axis.clabel(contour_lines, inline=True, fontsize=8)

    # Plot CS outline
    axis.plot(
        [r_cs_inner, r_cs_inner],
        [-dz_cs_full / 2, dz_cs_full / 2],
        "k-",
        linewidth=2,
        label="CS Inner",
    )
    axis.plot(
        [r_cs_outer, r_cs_outer],
        [-dz_cs_full / 2, dz_cs_full / 2],
        "k-",
        linewidth=2,
        label="CS Outer",
    )
    axis.plot(
        [r_cs_inner, r_cs_outer],
        [dz_cs_full / 2, dz_cs_full / 2],
        "k-",
        linewidth=2,
    )
    axis.plot(
        [r_cs_inner, r_cs_outer],
        [-dz_cs_full / 2, -dz_cs_full / 2],
        "k-",
        linewidth=2,
    )

    cbar = add_colourbar(contour_fill, axis, colorbar_axis)
    cbar.set_label("Radial Stress (MPa)")

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.minorticks_on()
    axis.set_xlim(r_cs_inner * 0.9, r_cs_outer * 1.1)
    axis.set_ylim((-dz_cs_full / 2) * 1.1, (dz_cs_full / 2) * 1.1)
    axis.grid(True, alpha=0.3)


__all__ = [
    "plot_cs_radial_stress_contour_profile",
    "plot_cs_radial_stress_profile",
    "plot_larmor_radius_profile",
    "plot_plasma_pressure_gradient_profiles",
    "plot_rad_contour",
    "read_imprad_data",
]
