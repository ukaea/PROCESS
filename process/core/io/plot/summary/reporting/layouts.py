"""Reporting functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

from process.core.io.plot.summary.constants import (
    BLANKET_COLOUR,
    FIRSTWALL_COLOUR,
    PLASMA_COLOUR,
    SHIELD_COLOUR,
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    VESSEL_COLOUR,
)

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes

    from process.core.io.mfile import MFile


def plot_upper_vertical_build(
    axis: plt.Axes, mfile: MFile, colour_scheme: Literal[1, 2]
):
    """Plots the upper vertical build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's vertical
    build
    (such as plasma, first wall, divertor, shield, vacuum vessel, thermal shield, TF
    coil, etc.)
    as a vertical stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile`, and each segment is color-coded and labeled accordingly.

    Parameters
    ----------
    axis:
        The matplotlib axis on which to plot the vertical build.
    mfile:
        An object containing the machine build data, with required fields for each
        vertical component.
    colour_scheme:
        Colour scheme index to use for component colors.

    Notes
    -----
    This function modifies the provided axis in-place and does not return a value.
    - Components with zero thickness are omitted from the plot.
    - The legend displays the name and thickness (in meters) of each component.
    """
    if mfile.get("i_single_null", scan=-1) == 1:
        upper_vertical_variables = [
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_fw_upper",
            "dz_blkt_upper",
            "dr_shld_blkt_gap",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
            "dz_tf_cryostat",
        ]
        upper_vertical_labels = [
            "Plasma Height",
            "First Wall - Plasma Gap",
            "First Wall Upper",
            "Blanket Upper",
            "Shield-Blanket Gap",
            "Shield Upper",
            "Vacuum Vessel Upper",
            "Shield-VV Gap",
            "Thermal Shield",
            "TF Coil - Shield Gap",
            "TF Coil",
            "TF Coil - Cryostat gap",
        ]
        upper_vertical_colours = [
            PLASMA_COLOUR[colour_scheme - 1],
            "white",
            FIRSTWALL_COLOUR[colour_scheme - 1],
            BLANKET_COLOUR[colour_scheme - 1],
            "white",
            SHIELD_COLOUR[colour_scheme - 1],
            VESSEL_COLOUR[colour_scheme - 1],
            "white",
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
            "white",
            (
                TFC_COLOUR[colour_scheme - 1]
                if mfile.get("i_tf_sup", scan=-1) != 0
                else "#b87333"
            ),
            "white",
        ]
    # Double null case
    else:
        upper_vertical_variables = [
            "z_plasma_xpoint_upper",
            "dz_xpoint_divertor",
            "dz_divertor",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
            "dz_tf_cryostat",
        ]
        upper_vertical_labels = [
            "Plasma Height",
            "Plasma - Divertor Gap",
            "Divertor Upper",
            "Shield Upper",
            "Vacuum Vessel Upper",
            "Shield-VV Gap",
            "Thermal Shield",
            "TF Coil - Shield Gap",
            "TF Coil",
            "TF Coil - Cryostat gap",
        ]
        upper_vertical_colours = [
            PLASMA_COLOUR[colour_scheme - 1],
            "white",
            "black",
            SHIELD_COLOUR[colour_scheme - 1],
            VESSEL_COLOUR[colour_scheme - 1],
            "white",
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
            "white",
            (
                TFC_COLOUR[colour_scheme - 1]
                if mfile.get("i_tf_sup", scan=-1) != 0
                else "#b87333"
            ),
            "white",
        ]

    # Get thicknesses for each layer
    upper_vertical_build = np.array([
        mfile.get(rl, scan=-1) for rl in upper_vertical_variables
    ])

    # Remove build parts equal to zero
    mask = ~(upper_vertical_build == 0.0)  # noqa: RUF069
    filtered_build = upper_vertical_build[mask]
    filtered_labels = [lbl for i, lbl in enumerate(upper_vertical_labels) if mask[i]]
    filtered_colors = [col for i, col in enumerate(upper_vertical_colours) if mask[i]]
    filtered_vars = [v for i, v in enumerate(upper_vertical_variables) if mask[i]]

    # Compute cumulative positions (bottoms) for stacking
    bottoms = np.zeros_like(filtered_build)
    for i in range(1, len(filtered_build)):
        bottoms[i] = bottoms[i - 1] + filtered_build[i - 1]

    # Plot each layer as a bar, stacking upwards from zero
    for kk in range(len(filtered_build)):
        axis.bar(
            0,
            filtered_build[kk],
            bottom=bottoms[kk],
            width=0.8,
            label=(
                f"{filtered_labels[kk]}\n[{filtered_vars[kk]}]\n{filtered_build[kk]:.3f} m"  # noqa: E501
            ),
            color=filtered_colors[kk],
            edgecolor="black",
            linewidth=0.05,
        )

    axis.set_xticks([])
    axis.legend(
        bbox_to_anchor=(0, 0),
        loc="upper left",
        ncol=6,
    )
    axis.minorticks_on()
    axis.set_ylabel("Height [m]")
    axis.title.set_text("Upper Vertical Build")


def draw_bend(
    ax: Axes,
    elbow_radius: float,
    theta_span: float,
    radius_pipe: float,
    title: str = "Bend",
    alpha: float = 0.8,
):
    """
    Draws a circular pipe bend with centerline and inner/outer boundaries.

    Parameters
    ----------
    ax:
        Target axes for plotting.
    elbow_radius:
        Radius of the elbow [m].
    theta_span:
        Array of angles [0, θ] where θ is pi/2 or pi [rad]
    radius_pipe:
        Pipe radius (fallback to 0.1m if not provided) [m]
    title:
        Plot title string.
    alpha:
        fill opacity
    """
    # Convert all inputs to mm
    elbow_radius_mm = elbow_radius * 1000
    pipe_radius_mm = radius_pipe * 1000

    theta = np.linspace(0, theta_span, 100)
    x_center = elbow_radius_mm * np.cos(theta)
    y_center = elbow_radius_mm * np.sin(theta)

    # Outer and inner walls (offset by ± pipe radius in mm)
    x_outer = (elbow_radius_mm + pipe_radius_mm) * np.cos(theta)
    y_outer = (elbow_radius_mm + pipe_radius_mm) * np.sin(theta)
    x_inner = (elbow_radius_mm - pipe_radius_mm) * np.cos(theta)
    y_inner = (elbow_radius_mm - pipe_radius_mm) * np.sin(theta)

    # Plot
    ax.plot(x_center, y_center, color="black", linestyle="--", label="Centerline")
    ax.plot(x_outer, y_outer, color="black")
    ax.plot(x_inner, y_inner, color="black")
    ax.fill(
        np.concatenate([x_outer, x_inner[::-1]]),
        np.concatenate([y_outer, y_inner[::-1]]),
        color="lightgrey",
        alpha=alpha,
    )

    ax.set_aspect("equal")
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.3)

    # Legend: Centerline + pipe radius info
    legend_text = (
        f"Centerline\nPipe radius: {pipe_radius_mm:.2f} mm\nElbow radius:"
        f" {elbow_radius_mm:.2f} mm"
    )
    ax.legend([legend_text], loc="upper right")


__all__ = ["draw_bend", "plot_upper_vertical_build"]
