"""Common functions for PROCESS summary plots."""

from __future__ import annotations

from importlib import resources
from typing import TYPE_CHECKING, Literal

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.patches import FancyBboxPatch

from process.core.io.plot.summary.constants import (
    BLANKET_COLOUR,
    CRYOSTAT_COLOUR,
    CSCOMPRESSION_COLOUR,
    FIRSTWALL_COLOUR,
    NBSHIELD_COLOUR,
    PLASMA_COLOUR,
    SOLENOID_COLOUR,
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    VESSEL_COLOUR,
)
from process.models.pulse import PulseTimings

if TYPE_CHECKING:
    from process.core.io.mfile import MFile

__all__ = [
    "add_colourbar",
    "box_style",
    "color_key",
    "get_pulse_timings",
    "load_plot_image",
    "place_plot_image",
    "setup_axis",
    "text_layout",
]


def load_plot_image(name: str):
    """Load an image bundled with the PROCESS plotting resources."""
    image = resources.files("process.core.io.plot.images").joinpath(name)
    with image.open("rb") as image_file:
        return mpimg.imread(image_file)


def place_plot_image(axis, name: str, bounds, *, zorder: int = 10):
    """Place a bundled image in an inset axes and hide its frame."""
    image_axis = axis.inset_axes(bounds, transform=axis.transAxes, zorder=zorder)
    image_axis.imshow(load_plot_image(name))
    image_axis.axis("off")
    return image_axis


def get_pulse_timings(mfile, scan: int) -> PulseTimings:
    """Construct pulse timings from an MFILE scan."""
    keys = (
        "t_plant_pulse_coil_precharge",
        "t_plant_pulse_plasma_current_ramp_up",
        "t_plant_pulse_fusion_ramp",
        "t_plant_pulse_burn",
        "t_plant_pulse_plasma_current_ramp_down",
        "t_plant_pulse_dwell",
    )
    return PulseTimings(**{key: mfile.get(key, scan=scan) for key in keys})


def box_style(colour: str, alpha: float = 1.0, linewidth: float | None = 2, **kwargs):
    """Matplotlib boxstyle"""
    return {
        "boxstyle": "round",
        "facecolor": colour,
        "alpha": alpha,
        "linewidth": linewidth,
        **kwargs,
    }


def text_layout(fig, h_align="left", v_align="bottom"):
    """Matplotlib text layout"""
    return {
        "fontsize": 9,
        "verticalalignment": v_align,
        "horizontalalignment": h_align,
        "transform": fig.transFigure,
    }


def setup_axis(axis, xmin, xmax, ymin, ymax):
    """Matplotlib axis setup"""
    axis.set_ylim(ymin, ymax)
    axis.set_xlim(xmin, xmax)
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)


def add_colourbar(contour_fill, axis, colourbar_axis):
    """Add colourbar to axis"""
    # Use a dedicated colorbar axes when provided so the main axes width is unchanged.
    if colourbar_axis is None:
        return axis.figure.colorbar(contour_fill, ax=axis, pad=0.02)
    return axis.figure.colorbar(contour_fill, cax=colourbar_axis)


def add_panel(
    axis,
    x: float,
    top: float,
    width: float,
    height: float,
    text: str,
    style: dict,
    text_x_offset: float = 0.004,
):
    """Add one fixed-size box and its unchanged original text."""
    box_expansion = 0.002
    left = x - box_expansion
    bottom = top - height - box_expansion
    expanded_width = width + 2 * box_expansion
    expanded_height = height + 2 * box_expansion
    fig = axis.get_figure()
    patch = FancyBboxPatch(
        (left, bottom),
        expanded_width,
        expanded_height,
        boxstyle="round,pad=0,rounding_size=0.003",
        transform=fig.transFigure,
        facecolor=style.get("facecolor", "white"),
        edgecolor=style.get("edgecolor", "black"),
        linewidth=style.get("linewidth", 1.0),
        alpha=style.get("alpha", 1.0),
        clip_on=False,
        zorder=20,
    )
    axis.add_patch(patch)

    return axis.text(
        x + text_x_offset,
        top - 0.004,
        text,
        color="black",
        clip_on=False,
        zorder=21,
        linespacing=1,
        **text_layout(fig, v_align="top"),
    )


def add_symbol(axis, panel_name: str, x: float, top: float):
    """Add a large symbol at its original figure-relative position."""
    return axis.text(
        x,
        top,
        panel_name,
        fontsize=23,
        color="black",
        verticalalignment="top",
        transform=axis.get_figure().transFigure,
        clip_on=False,
        zorder=22,
    )


def color_key(axis: plt.Axes, mfile: MFile, scan: int, colour_scheme: Literal[1, 2]):
    """Function to plot the colour key"""
    axis.set_ylim(0, 10)
    axis.set_xlim(0, 10)
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    labels = [
        ("CS coil", SOLENOID_COLOUR[colour_scheme - 1]),
        ("CS comp", CSCOMPRESSION_COLOUR[colour_scheme - 1]),
        (
            "TF coil",
            (
                TFC_COLOUR[colour_scheme - 1]
                if mfile.get("i_tf_sup", scan=scan) != 0
                else "#b87333"
            ),
        ),
        ("Thermal shield", THERMAL_SHIELD_COLOUR[colour_scheme - 1]),
        ("VV & shield", VESSEL_COLOUR[colour_scheme - 1]),
        ("Blanket", BLANKET_COLOUR[colour_scheme - 1]),
        ("First wall", FIRSTWALL_COLOUR[colour_scheme - 1]),
        ("Plasma", PLASMA_COLOUR[colour_scheme - 1]),
        ("PF coils", "none"),
        ("Divertor", "black"),
    ]

    if (mfile.get("i_hcd_primary", scan=scan) in {5, 8}) or (
        mfile.get("i_hcd_secondary", scan=scan) in {5, 8}
    ):
        labels.extend((
            ("NB duct shield", NBSHIELD_COLOUR[colour_scheme - 1]),
            ("Cryostat", CRYOSTAT_COLOUR[colour_scheme - 1]),
        ))
    else:
        labels.append(("Cryostat", CRYOSTAT_COLOUR[colour_scheme - 1]))

    for i, (text, color) in enumerate(labels):
        row = i // 4
        col = i % 4
        y_pos = 9 - row * 1.5
        x_pos = col * 2.5

        axis.text(x_pos, y_pos, text, ha="left", va="top", size="small")
        axis.add_patch(
            patches.Rectangle(
                (x_pos + 1.5, y_pos - 0.35),
                0.5,
                0.4,
                lw=0 if color != "none" else 1,
                facecolor=color if color != "none" else "none",
                edgecolor="black" if color == "none" else "none",
            )
        )
