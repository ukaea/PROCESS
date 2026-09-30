"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.common import (
    add_colourbar,
    get_pulse_timings,
)
from process.models.engineering.materials import (
    calculate_tresca_stress,
    calculate_von_mises_stress,
    poisson_steel,
)
from process.models.pfcoil import N_CS_STRESS_PROFILE_POINTS, CSCoil

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_cs_stress_time_profile(axis: plt.Axes, mfile: MFile, scan: int) -> None:
    """Function to plot the time profile of the CS stress during the pulse."""
    pulse_timings = get_pulse_timings(mfile, scan)

    stress_z_cs_self_midplane_profile = np.zeros(pulse_timings.n_pf_active_points_total)
    for i in range(pulse_timings.n_pf_active_points_total):
        stress_z_cs_self_midplane_profile[i] = mfile.get(
            f"stress_z_cs_self_midplane_profile[{i}]", scan=scan
        )

    # Plot stress vs time
    axis.plot(
        pulse_timings.pf_active_cumulative,
        stress_z_cs_self_midplane_profile / 1e6,
        "o-",
        linewidth=2,
        markersize=4,
        label="$\\sigma_{z}$,Axial Stress",
    )
    axis.set_xlabel("Pulse Time (s)")
    axis.set_ylabel("Midplane Axial Stress (MPa)")
    axis.minorticks_on()
    axis.legend(loc="best")
    axis.set_title("CS Midplane Axial Stress Time Profile")
    axis.grid(True, alpha=0.3)


def plot_cs_hoop_stress_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    j_cs: float,
    b_cs_inner: float,
):
    """Plot CS hoop stress profile"""
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)

    radii = np.linspace(r_cs_inner, r_cs_outer, num=10)
    stress_values = np.array([
        CSCoil.calculate_cs_hoop_stress(
            r_stress_point=radius,
            r_cs_inner=r_cs_inner,
            r_cs_outer=r_cs_outer,
            j_cs=j_cs,
            b_cs_inner=b_cs_inner,
            f_poisson_cs_structure=poisson_steel,
            f_a_cs_turn_steel=mfile.get("f_a_cs_turn_steel", scan=scan),
        )
        for radius in radii
    ])

    axis.plot(
        radii,
        stress_values / 1e6,
        linewidth=2,
        label="$\\sigma_{\\theta}$,Hoop Stress",
    )
    max_idx = np.argmax(np.abs(stress_values))
    max_radius = radii[max_idx]
    max_stress = stress_values[max_idx] / 1e6
    axis.axvline(max_radius, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.axhline(max_stress, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.set_xlabel("Radial Position (m)")
    axis.set_ylabel("Hoop Stress (MPa)")
    axis.minorticks_on()
    axis.set_title("CS Hoop Stress at BOP")
    axis.grid(True, alpha=0.3)


def plot_cs_hoop_stress_contour_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    j_cs: float,
    b_cs_inner: float,
    colorbar_axis: plt.Axes | None = None,
):
    """Plot CS hoop stress contour profile"""
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)
    f_a_cs_turn_steel = mfile.get("f_a_cs_turn_steel", scan=scan)

    # Create 2D grid for contour plot: radial and vertical dimensions
    n_radial = 50
    radial_grid = np.linspace(r_cs_inner, r_cs_outer, n_radial)
    height_grid = np.linspace(
        -dz_cs_full / 2, dz_cs_full / 2, N_CS_STRESS_PROFILE_POINTS
    )

    # Create meshgrid for filled contour
    r, z = np.meshgrid(radial_grid, height_grid)

    # Calculate hoop stress across the 2D grid
    stress_data = np.zeros((len(height_grid), n_radial))
    for i in range(len(height_grid)):
        for j in range(n_radial):
            stress_data[i, j] = (
                CSCoil.calculate_cs_hoop_stress(
                    r_stress_point=radial_grid[j],
                    r_cs_inner=r_cs_inner,
                    r_cs_outer=r_cs_outer,
                    j_cs=j_cs,
                    b_cs_inner=b_cs_inner,
                    f_poisson_cs_structure=poisson_steel,
                    f_a_cs_turn_steel=f_a_cs_turn_steel,
                )
                / 1e6
            )

    # Plot filled contour of stress distribution
    contour_fill = axis.contourf(r, z, stress_data, levels=15, cmap="RdYlBu_r")
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
    cbar.set_label("Hoop Stress (MPa)")

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.minorticks_on()
    axis.set_xlim(r_cs_inner * 0.9, r_cs_outer * 1.1)
    axis.set_ylim((-dz_cs_full / 2) * 1.1, (dz_cs_full / 2) * 1.1)
    axis.grid(True, alpha=0.3)


def plot_cs_vertical_stress_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
):
    """Plot CS vertical stress profile"""
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)

    stress_z_profile = np.array([
        float(mfile.data[f"stress_z_cs_self_profile_{i}"].get_scan(scan)) / 1e6
        for i in range(N_CS_STRESS_PROFILE_POINTS)
    ])
    z_positions = np.linspace(-dz_cs_full / 2, dz_cs_full / 2, len(stress_z_profile))

    axis.plot(
        stress_z_profile,
        z_positions,
        linewidth=2,
        label="$\\sigma_{z}$,Vertical Stress",
    )
    max_idx = np.argmax(np.abs(stress_z_profile))
    max_stress = stress_z_profile[max_idx]
    max_z = z_positions[max_idx]
    axis.axvline(max_stress, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.axhline(max_z, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
    axis.set_xlabel("Vertical Stress (MPa)")
    axis.set_ylabel("Z [m]")
    axis.minorticks_on()
    axis.grid(True, alpha=0.3)
    axis.set_title("CS Vertical Stress at BOP")


def plot_vertical_stress_contour_profile(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colorbar_axis: plt.Axes | None = None,
):
    """Vertical stress contour plot"""
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)

    stress_z_profile = [
        float(mfile.data[f"stress_z_cs_self_profile_{i}"].get_scan(scan)) / 1e6
        for i in range(N_CS_STRESS_PROFILE_POINTS)
    ]

    # Create 2D grid for contour plot: radial and vertical dimensions
    n_radial = 50
    radial_grid = np.linspace(r_cs_inner, r_cs_outer, n_radial)
    height_grid = np.linspace(-dz_cs_full / 2, dz_cs_full / 2, len(stress_z_profile))

    # Create meshgrid for filled contour
    r, z = np.meshgrid(radial_grid, height_grid)

    # Interpolate stress values across radial direction (assume linear variation)
    stress_data = np.zeros((len(stress_z_profile), n_radial))
    for i, stress_val in enumerate(stress_z_profile):
        stress_data[i, :] = stress_val

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
    cbar.set_label("Vertical Stress (MPa)")

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.minorticks_on()
    axis.set_xlim(r_cs_inner * 0.9, r_cs_outer * 1.1)
    axis.set_ylim((-dz_cs_full / 2) * 1.1, (dz_cs_full / 2) * 1.1)
    axis.grid(True, alpha=0.3)


def plot_cs_tresca_2d_contour(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colorbar_axis: plt.Axes | None = None,
):
    """CS Tresca stress contour plot"""
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)
    j_cs = mfile.get("j_cs_pulse_start", scan=scan)
    b_cs_inner = mfile.get("b_cs_peak_pulse_start", scan=scan)
    f_a_cs_turn_steel = mfile.get("f_a_cs_turn_steel", scan=scan)

    stress_z_profile = np.array([
        float(mfile.data[f"stress_z_cs_self_profile_{i}"].get_scan(scan))
        for i in range(N_CS_STRESS_PROFILE_POINTS)
    ])

    # Create 2D grid for contour plot: radial and vertical dimensions
    n_radial = 50
    radial_grid = np.linspace(r_cs_inner, r_cs_outer, n_radial)
    height_grid = np.linspace(-dz_cs_full / 2, dz_cs_full / 2, len(stress_z_profile))

    # Create meshgrid for filled contour
    r, z = np.meshgrid(radial_grid, height_grid)

    # Calculate Tresca stress across the coil cross-section.
    tresca_data = np.zeros((len(height_grid), n_radial))
    for i, stress_z in enumerate(stress_z_profile):
        for j, radius in enumerate(radial_grid):
            stress_hoop = CSCoil.calculate_cs_hoop_stress(
                r_stress_point=radius,
                r_cs_inner=r_cs_inner,
                r_cs_outer=r_cs_outer,
                j_cs=j_cs,
                b_cs_inner=b_cs_inner,
                f_poisson_cs_structure=poisson_steel,
                f_a_cs_turn_steel=f_a_cs_turn_steel,
            )
            stress_radial = CSCoil.calculate_cs_radial_stress(
                r_stress_point=radius,
                r_cs_inner=r_cs_inner,
                r_cs_outer=r_cs_outer,
                j_cs=j_cs,
                b_cs_inner=b_cs_inner,
                f_poisson_cs_structure=poisson_steel,
            )
            tresca_data[i, j] = (
                calculate_tresca_stress(
                    stress_x=stress_hoop,
                    stress_y=stress_z,
                    stress_z=stress_radial,
                )
                / 1e6
            )

    # Plot filled contour of Tresca stress distribution
    contour_lines = axis.contour(
        r,
        z,
        tresca_data,
        levels=[tresca_data.max()],
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

    contour_fill = axis.contourf(r, z, tresca_data, levels=15, cmap="RdYlBu_r")
    cbar = add_colourbar(contour_fill, axis, colorbar_axis)
    cbar.set_label("Tresca Stress (MPa)")

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.minorticks_on()
    axis.set_xlim(r_cs_inner * 0.9, r_cs_outer * 1.1)
    axis.set_ylim((-dz_cs_full / 2) * 1.1, (dz_cs_full / 2) * 1.1)
    axis.grid(True, alpha=0.3)
    axis.set_title("CS Tresca Stress Contour at BOP")


def plot_cs_von_mises_2d_contour(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colorbar_axis: plt.Axes | None = None,
):
    """CS Von Mises stress contour plot"""
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)
    r_cs_inner = mfile.get("r_cs_inner", scan=scan)
    r_cs_outer = mfile.get("r_cs_outer", scan=scan)
    j_cs = mfile.get("j_cs_pulse_start", scan=scan)
    b_cs_inner = mfile.get("b_cs_peak_pulse_start", scan=scan)
    f_a_cs_turn_steel = mfile.get("f_a_cs_turn_steel", scan=scan)

    stress_z_profile = np.array([
        float(mfile.data[f"stress_z_cs_self_profile_{i}"].get_scan(scan))
        for i in range(N_CS_STRESS_PROFILE_POINTS)
    ])

    # Create 2D grid for contour plot: radial and vertical dimensions
    n_radial = 50
    radial_grid = np.linspace(r_cs_inner, r_cs_outer, n_radial)
    height_grid = np.linspace(-dz_cs_full / 2, dz_cs_full / 2, len(stress_z_profile))

    # Create meshgrid for filled contour
    r, z = np.meshgrid(radial_grid, height_grid)

    # Calculate Von Mises stress across the coil cross-section.
    von_mises_data = np.zeros((len(height_grid), n_radial))
    for i, stress_z in enumerate(stress_z_profile):
        for j, radius in enumerate(radial_grid):
            stress_hoop = CSCoil.calculate_cs_hoop_stress(
                r_stress_point=radius,
                r_cs_inner=r_cs_inner,
                r_cs_outer=r_cs_outer,
                j_cs=j_cs,
                b_cs_inner=b_cs_inner,
                f_poisson_cs_structure=poisson_steel,
                f_a_cs_turn_steel=f_a_cs_turn_steel,
            )
            stress_radial = CSCoil.calculate_cs_radial_stress(
                r_stress_point=radius,
                r_cs_inner=r_cs_inner,
                r_cs_outer=r_cs_outer,
                j_cs=j_cs,
                b_cs_inner=b_cs_inner,
                f_poisson_cs_structure=poisson_steel,
            )
            von_mises_data[i, j] = (
                calculate_von_mises_stress(
                    stress_x=stress_hoop,
                    stress_y=stress_z,
                    stress_z=stress_radial,
                    stress_shear_xy=0.0,
                    stress_shear_yz=0.0,
                    stress_shear_zx=0.0,
                )
                / 1e6
            )

    # Plot filled contour of Von Mises stress distribution
    contour_lines = axis.contour(
        r,
        z,
        von_mises_data,
        levels=[von_mises_data.max()],
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

    contour_fill = axis.contourf(r, z, von_mises_data, levels=15, cmap="RdYlBu_r")
    cbar = add_colourbar(contour_fill, axis, colorbar_axis)
    cbar.set_label("Von Mises Stress (MPa)")

    axis.set_xlabel("R [m]")
    axis.minorticks_on()
    axis.set_xlim(r_cs_inner * 0.9, r_cs_outer * 1.1)
    axis.set_ylim((-dz_cs_full / 2) * 1.1, (dz_cs_full / 2) * 1.1)
    axis.grid(True, alpha=0.3)
    axis.set_title("CS Von Mises Stress Contour at BOP")


__all__ = [
    "plot_cs_hoop_stress_contour_profile",
    "plot_cs_hoop_stress_profile",
    "plot_cs_stress_time_profile",
    "plot_cs_tresca_2d_contour",
    "plot_cs_vertical_stress_profile",
    "plot_cs_von_mises_2d_contour",
    "plot_vertical_stress_contour_profile",
]
