"""Geometry functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

from process.core.io.plot.summary.common import (
    box_style,
    text_layout,
)
from process.core.io.plot.summary.geometry.poloidal import (
    plot_blanket,
    plot_firstwall,
)
from process.core.io.plot.summary.plasma.physics import (
    plot_plasma,
)
from process.core.io.plot.summary.reporting.layouts import (
    draw_bend,
)
from process.data_structure.physics_variables import DivertorNumberModels

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_blkt_pipe_bends(fig, m_file, scan: int):
    """Plot the blanket pipe bends on the given axis, with axes in mm.

    Parameters
    ----------
    fig :

    m_file :

    scan: int :

    """
    ax_90 = fig.add_subplot(341)
    ax_180 = fig.add_subplot(342)

    r = m_file.get("radius_blkt_channel", scan=scan)
    fallback_radius = 0.1  # meters

    elbow_radius_90 = (
        m_file.get("radius_blkt_channel_90_bend", scan=scan) or fallback_radius
    )
    elbow_radius_180 = (
        m_file.get("radius_blkt_channel_180_bend", scan=scan) or fallback_radius
    )

    draw_bend(ax_90, elbow_radius_90, np.pi / 2, r, title="Blanket Pipe 90° Bend")
    draw_bend(ax_180, elbow_radius_180, np.pi, r, title="Blanket Pipe 180° Bend")


def plot_blkt_structure(
    ax: plt.Axes,
    fig: plt.Figure,
    m_file: MFile,
    scan: int,
    radial_build: dict[str, float],
    colour_scheme: Literal[1, 2],
):
    """Plot the blkt structure and relevant angles"""
    # MFILE variables needed to plot the blkt structure and angles
    rmajor = m_file.get("rmajor", scan=scan)
    rminor = m_file.get("rminor", scan=scan)
    dr_fw_plasma_gap_outboard = m_file.get("dr_fw_plasma_gap_outboard", scan=scan)
    dr_fw_plasma_gap_inboard = m_file.get("dr_fw_plasma_gap_inboard", scan=scan)
    dr_fw_inboard = m_file.get("dr_fw_inboard", scan=scan)
    dr_fw_outboard = m_file.get("dr_fw_outboard", scan=scan)
    dr_blkt_outboard = m_file.get("dr_blkt_outboard", scan=scan)
    dr_blkt_inboard = m_file.get("dr_blkt_inboard", scan=scan)
    dz_blkt_half = m_file.get("dz_blkt_half", scan=scan)
    deg_blkt_outboard_poloidal_plasma = m_file.get(
        "deg_blkt_outboard_poloidal_plasma", scan=scan
    )
    deg_blkt_inboard_poloidal_plasma = m_file.get(
        "deg_blkt_inboard_poloidal_plasma", scan=scan
    )
    f_deg_blkt_outboard_poloidal_plasma = m_file.get(
        "f_deg_blkt_outboard_poloidal_plasma", scan=scan
    )
    f_deg_blkt_inboard_poloidal_plasma = m_file.get(
        "f_deg_blkt_inboard_poloidal_plasma", scan=scan
    )
    deg_div_poloidal_plasma = m_file.get("deg_div_poloidal_plasma", scan=scan)
    f_ster_div_single = m_file.get("f_ster_div_single", scan=scan)
    i_single_null = m_file.get("i_single_null", scan=scan)

    # ======================

    plot_blanket(ax, m_file, scan, radial_build, colour_scheme)
    plot_plasma(ax, m_file, scan, colour_scheme)
    plot_firstwall(ax, m_file, scan, radial_build, colour_scheme)

    ax.set_xlabel("Radial position [m]")
    ax.set_ylabel("Vertical position [m]")
    ax.set_title("Blanket and First Wall Poloidal Cross-Section")
    ax.minorticks_on()
    ax.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)

    r_blkt_outboard_out = (
        rmajor + rminor + dr_fw_outboard + dr_fw_plasma_gap_outboard + dr_blkt_outboard
    )
    r_blkt_inboard_in = (
        rmajor - rminor - dr_fw_plasma_gap_inboard - dr_fw_inboard - dr_blkt_inboard
    )
    r_fw_outboard_in = r_blkt_outboard_out - dr_blkt_outboard - dr_fw_outboard
    r_fw_inboard_out = r_blkt_inboard_in + dr_blkt_inboard + dr_fw_inboard

    # Plot a horizontal line at dz_blkt_half (blanket half height)
    for dz_blkt in (dz_blkt_half, -dz_blkt_half):
        ax.axhline(
            dz_blkt,
            color="purple",
            linestyle="--",
            linewidth=1.5,
            label="Blanket Half Height",
        )

    if DivertorNumberModels(i_single_null) == DivertorNumberModels.DOUBLE_NULL:
        # Plot arrows for the outboard blanket angles
        ax.annotate(
            "",
            xy=(rmajor, 0),
            xytext=(rmajor, dz_blkt_half),
            arrowprops={"arrowstyle": "<-", "color": "purple"},
            zorder=5,
        )
    # If single null then only plot the lower arrow for the outboard blanket angle
    ax.annotate(
        "",
        xy=(rmajor, 0),
        xytext=(rmajor, -dz_blkt_half),
        arrowprops={"arrowstyle": "<-", "color": "purple"},
        zorder=5,
    )

    # Plot arc showing the angle between the two outboard blanket arrows
    arc_radius = 1.0

    # 3 to 6 o'clock position is -90 degrees,
    angle_start = -90.0
    match DivertorNumberModels(i_single_null):
        case DivertorNumberModels.SINGLE_NULL:
            angle_end = 90.0 + deg_div_poloidal_plasma
        case DivertorNumberModels.DOUBLE_NULL:
            # 3 to 12 o'clock position is +90 degrees
            angle_end = 90.0

    theta = np.linspace(np.deg2rad(angle_start), np.deg2rad(angle_end), 50)
    arc_x = rmajor + arc_radius * np.cos(theta)
    arc_y = arc_radius * np.sin(theta)

    ax.plot(arc_x, arc_y, color="purple", linewidth=2)

    # Add angle label at the arc
    mid_angle = np.deg2rad((angle_start + angle_end) / 2)
    label_radius = arc_radius * 1.8
    label_x = rmajor + label_radius * np.cos(mid_angle)
    label_y = label_radius * np.sin(mid_angle)

    # Plot the info box for the outboard blanket
    ax.text(
        label_x,
        label_y,
        f"{deg_blkt_outboard_poloidal_plasma:.1f}°\n({f_deg_blkt_outboard_poloidal_plasma * 100:.1f}%)",  # noqa: E501
        fontsize=7,
        color="purple",
        ha="center",
        va="center",
        weight="bold",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "purple",
            "linewidth": 1.5,
        },
    )

    # Plot arrows for the inboard blanket angles
    for dz_blkt in (dz_blkt_half, -dz_blkt_half):
        ax.annotate(
            "",
            xy=(rmajor, 0),
            xytext=(r_fw_inboard_out, dz_blkt),
            arrowprops={"arrowstyle": "<-", "color": "green"},
            zorder=5,
        )

    # Plot arc showing the angle between the two inboard blanket arrows
    arc_radius = 1.0
    angle_start = -deg_blkt_inboard_poloidal_plasma / 2
    angle_end = deg_blkt_inboard_poloidal_plasma / 2

    theta = np.linspace(np.deg2rad(angle_start), np.deg2rad(angle_end), 50)
    arc_x = rmajor - arc_radius * np.cos(theta)
    arc_y = arc_radius * np.sin(theta)

    ax.plot(arc_x, arc_y, color="green", linewidth=2)

    # Add angle label at the arc
    mid_angle = np.deg2rad((angle_start + angle_end) / 2)
    label_radius = arc_radius * 1.8
    label_x = rmajor - label_radius * np.cos(mid_angle)
    label_y = label_radius * np.sin(mid_angle)

    # Plot the info box for the inboard blanket
    ax.text(
        label_x,
        label_y,
        f"{deg_blkt_inboard_poloidal_plasma:.1f}°\n({f_deg_blkt_inboard_poloidal_plasma * 100:.1f}%)",  # noqa: E501
        fontsize=7,
        color="green",
        ha="center",
        va="center",
        weight="bold",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "green",
            "linewidth": 1.5,
        },
        zorder=5,
    )

    # Plot arrows for the divertor angles
    # If double null then plot the upper also
    if DivertorNumberModels(i_single_null) == DivertorNumberModels.DOUBLE_NULL:
        # Plot arc showing the angle between the two arrows (divertor angle)
        arc_radius = 1.5
        # 3 to 12 o'clock position is +90 degrees,
        angle_start = 90.0
        angle_end = 90.0 + deg_div_poloidal_plasma

        theta = np.linspace(np.deg2rad(angle_start), np.deg2rad(angle_end), 50)
        arc_x = rmajor + arc_radius * np.cos(theta)
        arc_y = arc_radius * np.sin(theta)

        ax.plot(arc_x, arc_y, color="black", linewidth=2)

        # Add angle label at the arc
        mid_angle = np.deg2rad((angle_start + angle_end) / 2)
        label_radius = arc_radius * 1.8
        label_x = rmajor + label_radius * np.cos(mid_angle)
        label_y = label_radius * np.sin(mid_angle)

        ax.text(
            label_x,
            label_y,
            f"{deg_div_poloidal_plasma:.1f}°\n({f_ster_div_single * 100:.1f}%)",
            fontsize=7,
            color="black",
            ha="center",
            va="center",
            weight="bold",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.8,
                "edgecolor": "black",
                "linewidth": 1.5,
            },
            zorder=5,
        )

    # Plot arc showing the angle between the two arrows for the lower divertor (divertor
    # angle)
    arc_radius = 1.5
    # 3 to 6 o'clock is -90 degrees
    angle_start = -90.0
    angle_end = angle_start - deg_div_poloidal_plasma

    theta = np.linspace(np.deg2rad(angle_start), np.deg2rad(angle_end), 50)
    arc_x = rmajor + arc_radius * np.cos(theta)
    arc_y = arc_radius * np.sin(theta)

    ax.plot(arc_x, arc_y, color="black", linewidth=2)

    # Add angle label at the arc
    mid_angle = np.deg2rad((angle_start + angle_end) / 2)
    label_radius = arc_radius * 1.8
    label_x = rmajor + label_radius * np.cos(mid_angle)
    label_y = label_radius * np.sin(mid_angle)

    # Plot the info box for the lower divertor angle
    ax.text(
        label_x,
        label_y,
        f"{deg_div_poloidal_plasma:.1f}°\n({f_ster_div_single * 100:.1f}%)",
        fontsize=7,
        color="black",
        ha="center",
        va="center",
        weight="bold",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "black",
            "linewidth": 1.5,
        },
        zorder=5,
    )

    # Plot vertical lines at the inner and outer radial boundaries of the blanket
    linestyle = {
        "color": "black",
        "linestyle": "--",
        "linewidth": 1.5,
        "zorder": 10,
    }
    ax.axvline(r_blkt_inboard_in, **linestyle)
    ax.axvline(r_blkt_outboard_out, **linestyle)
    ax.axvline(r_fw_inboard_out, **linestyle)
    ax.axvline(r_fw_outboard_in, **linestyle)

    ax.axvline(
        rmajor,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="Major Radius $R_0$",
    )

    # Plot midplane line (horizontal dashed line at Z=0)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1.5, label="Midplane")

    textstr_blkt_areas = (
        "$\\mathbf{Blanket \\ Areas:}$\n\nInboard blanket, with holes and"
        f" gaps: {m_file.get('a_blkt_inboard_surface', scan=scan):,.3f}"
        " $\\text{m}^2$\nOutboard blanket, with holes and gaps:"
        f" {m_file.get('a_blkt_outboard_surface', scan=scan):,.3f}"
        " $\\text{m}^2$\nTotal blanket, with holes and gaps:"
        f" {m_file.get('a_blkt_total_surface', scan=scan):,.3f}"
        " $\\text{m}^2$\n\nInboard blanket, full coverage:"
        f" {m_file.get('a_blkt_inboard_surface_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^2$\nOutboard blanket, full coverage:"
        f" {m_file.get('a_blkt_outboard_surface_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^2$\nTotal blanket, full coverage:"
        f" {m_file.get('a_blkt_total_surface_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^2$ "
    )

    ax.text(
        0.05,
        0.3,
        textstr_blkt_areas,
        **text_layout(fig),
        bbox=box_style("wheat"),
    )

    textstr_blkt_volumes = (
        "$\\mathbf{Blanket \\ Volumes:}$\n\nInboard blanket, with holes and"
        f" gaps: {m_file.get('vol_blkt_inboard', scan=scan):,.3f}"
        " $\\text{m}^3$\nOutboard blanket, with holes and gaps:"
        f" {m_file.get('vol_blkt_outboard', scan=scan):,.3f}"
        " $\\text{m}^3$\nTotal blanket, with holes and gaps:"
        f" {m_file.get('vol_blkt_total', scan=scan):,.3f}"
        " $\\text{m}^3$\n\nInboard blanket, full coverage:"
        f" {m_file.get('vol_blkt_inboard_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^3$\nOutboard blanket, full coverage:"
        f" {m_file.get('vol_blkt_outboard_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^3$\nTotal blanket, full coverage:"
        f" {m_file.get('vol_blkt_total_full_coverage', scan=scan):,.3f}"
        " $\\text{m}^3$ "
    )

    ax.text(
        0.05,
        0.05,
        textstr_blkt_volumes,
        **text_layout(fig),
        bbox=box_style("wheat"),
    )


__all__ = ["plot_blkt_pipe_bends", "plot_blkt_structure"]
