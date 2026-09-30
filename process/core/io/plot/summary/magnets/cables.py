"""Magnets functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.patches import Rectangle

from process.core.io.plot.summary.common import (
    box_style,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_cable_in_conduit_cable(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plots TF coil CICC cable cross-section.

    Parameters
    ----------
    axis: plt.Axes :

    fig :

    mfile: MFile :

    scan: int :

    """
    dia_tf_turn_superconducting_cable = mfile.get(
        "dia_tf_turn_superconducting_cable", scan=scan
    )
    f_a_tf_turn_cable_copper = mfile.get("f_a_tf_turn_cable_copper", scan=scan)

    # Convert to mm
    dia_mm = dia_tf_turn_superconducting_cable * 1000
    radius_superconductor_mm = np.sqrt(1 - f_a_tf_turn_cable_copper) * (dia_mm / 2)

    # Draw the outer copper circle
    circle_copper_surrounding = patches.Circle(
        (0, 0),
        dia_mm / 2,
        facecolor="#b87333",  # copper color
        edgecolor="#8B4000",  # darker copper edge
        linewidth=0.1,
        alpha=0.8,
        label="Copper",
        zorder=1,
    )
    axis.add_patch(circle_copper_surrounding)

    # Draw the inner superconductor circle
    circle_central_conductor = patches.Circle(
        (0, 0),
        radius_superconductor_mm,
        facecolor="black",
        linewidth=0.3,
        alpha=0.7,
        label="Superconductor",
        zorder=2,
    )
    axis.add_patch(circle_central_conductor)

    # Convert cable diameter to mm
    cable_diameter_mm = mfile.get("dia_tf_turn_superconducting_cable", scan=scan) * 1000
    # Convert lengths from meters to kilometers for display
    len_tf_coil_superconductor_km = (
        mfile.get("len_tf_coil_superconductor", scan=scan) / 1000.0
    )
    len_tf_superconductor_total_km = (
        mfile.get("len_tf_superconductor_total", scan=scan) / 1000.0
    )

    textstr_cable = (
        f"$\\mathbf{{Cable:}}$\n\nCable diameter: {cable_diameter_mm:,.4f}"
        " mm\nCopper area fraction:"
        f" {mfile.get('f_a_tf_turn_cable_copper', scan=scan):.4f}\nNumber of"
        " strands per turn:"
        f" {int(mfile.get('n_tf_turn_superconducting_cables', scan=scan)):,}\nLength"
        f" of superconductor per coil: {len_tf_coil_superconductor_km:,.2f}"
        " km\nTotal length of superconductor in all coils:"
        f" {len_tf_superconductor_total_km:,.2f} km\n"
    )
    draw_text(
        axis,
        0.4,
        0.3,
        textstr_cable,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("#cccccc"),
    )

    axis.set_aspect("equal")
    axis.set_xlim(-dia_mm / 1.5, dia_mm / 1.5)
    axis.set_ylim(-dia_mm / 1.5, dia_mm / 1.5)
    axis.set_title("TF CICC Cable Cross-Section")
    axis.minorticks_on()
    axis.legend(loc="upper right")
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    axis.set_xlabel("X [mm]")
    axis.set_ylabel("Y [mm]")


def plot_hts_tape_geometry(
    axis,
    r_left: float,
    z_bottom: float,
    dr_hts_tape: float,
    dx_hts_tape_rebco: float,
    dx_hts_tape_copper: float,
    dx_hts_tape_hastelloy: float,
    show_legend: bool = True,
):
    """Plot HTS tape geometry"""
    legend_label = None if show_legend else "_nolegend_"
    # Plot a rectangular tape stack in the middle
    rect = Rectangle(
        (r_left, z_bottom),
        width=dr_hts_tape,
        height=dx_hts_tape_copper / 2,
        edgecolor=None,
        facecolor="#B87333",
        linewidth=2,
        label="Copper" if show_legend else legend_label,
    )
    axis.add_patch(rect)
    rect = Rectangle(
        (r_left, z_bottom + dx_hts_tape_copper / 2),
        width=dr_hts_tape,
        height=dx_hts_tape_hastelloy / 2,
        edgecolor=None,
        facecolor="grey",
        linewidth=2,
        label="Hastelloy" if show_legend else legend_label,
    )
    axis.add_patch(rect)
    rect = Rectangle(
        (
            r_left,
            z_bottom + dx_hts_tape_copper / 2 + dx_hts_tape_hastelloy / 2,
        ),
        width=dr_hts_tape,
        height=dx_hts_tape_rebco,
        edgecolor=None,
        facecolor="blue",
        linewidth=2,
        label="REBCO" if show_legend else legend_label,
    )
    axis.add_patch(rect)
    rect = Rectangle(
        (
            r_left,
            z_bottom
            + dx_hts_tape_copper / 2
            + dx_hts_tape_hastelloy / 2
            + dx_hts_tape_rebco,
        ),
        width=dr_hts_tape,
        height=dx_hts_tape_hastelloy / 2,
        edgecolor=None,
        facecolor="grey",
        linewidth=2,
        label="Hastelloy" if show_legend else legend_label,
    )
    axis.add_patch(rect)
    rect = Rectangle(
        (
            r_left,
            z_bottom
            + dx_hts_tape_copper / 2
            + dx_hts_tape_hastelloy / 2
            + dx_hts_tape_rebco
            + dx_hts_tape_hastelloy / 2,
        ),
        width=dr_hts_tape,
        height=dx_hts_tape_copper / 2,
        edgecolor=None,
        facecolor="#B87333",
        linewidth=2,
        label="Copper" if show_legend else legend_label,
    )
    axis.add_patch(rect)

    axis.set_title("HTS Tape Geometry")
    axis.grid(True)
    axis.set_xlabel("X-axis (m)")
    axis.set_ylabel("Y-axis (m)")
    axis.set_xlim(r_left * 0.9, dr_hts_tape * 1.1)
    axis.set_ylim(
        z_bottom * 0.9,
        (dx_hts_tape_copper + dx_hts_tape_hastelloy + dx_hts_tape_rebco) * 1.1,
    )
    axis.minorticks_on()
    axis.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
    if show_legend:
        axis.legend(loc="upper right")


__all__ = ["plot_cable_in_conduit_cable", "plot_hts_tape_geometry"]
