"""Geometry functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.path import Path as mplPath

from process.core.io.plot.summary.constants import (
    BLANKET_COLOUR,
    CRYOSTAT_COLOUR,
    CSCOMPRESSION_COLOUR,
    FIRSTWALL_COLOUR,
    NBSHIELD_COLOUR,
    PLASMA_COLOUR,
    SHIELD_COLOUR,
    SOLENOID_COLOUR,
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    VESSEL_COLOUR,
    rtangle,
)
from process.core.io.plot.summary.geometry.build import (
    cumulative_radial_build2,
)
from process.core.io.plot.summary.magnets import (
    TF_outboard,
)
from process.models.physics.current_drive import (
    CurrentDriveMethodType,
    CurrentDriveModel,
)
from process.models.tfcoil.base import TFConductorModel

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def toroidal_cross_section(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    demo_ranges: bool,
    colour_scheme: Literal[1, 2],
):
    """Function to plot toroidal cross-section"""
    axis.set_xlabel("R [m]")
    axis.set_ylabel("X [m]")
    axis.set_title("Toroidal Cross-Section")
    axis.minorticks_on()
    axis.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.2)

    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    r_cryostat_inboard = mfile.get("r_cryostat_inboard", scan=scan)
    dr_cryostat = mfile.get("dr_cryostat", scan=scan)
    n_tf_coils = mfile.get("n_tf_coils", scan=scan)
    if (
        CurrentDriveModel(mfile.get("i_hcd_primary", scan=scan)).method
        == CurrentDriveMethodType.NEUTRAL_BEAM
        or CurrentDriveModel(mfile.get("i_hcd_secondary", scan=scan)).method
        == CurrentDriveMethodType.NEUTRAL_BEAM
    ):
        dx_beam_shield = mfile.get("dx_beam_shield", scan=scan)
        dx_beam_duct = mfile.get("dx_beam_duct", scan=scan)
        radius_beam_tangency = mfile.get("radius_beam_tangency", scan=scan)
    else:
        dx_beam_shield = 0
        dx_beam_duct = 0
        radius_beam_tangency = 0

    dr_tf_outboard = mfile.get("dr_tf_outboard", scan=scan)
    full_angle = 2 * np.pi
    arc(axis, rmajor, theta2=full_angle, style="dashed")

    # Colour in the main components
    for v, colours in [
        ("dr_cs", SOLENOID_COLOUR[colour_scheme - 1]),
        ("dr_cs_precomp", CSCOMPRESSION_COLOUR[colour_scheme - 1]),
        (
            "dr_tf_inboard",
            (
                TFC_COLOUR[colour_scheme - 1]
                if TFConductorModel(mfile.get("i_tf_sup", scan=scan))
                != TFConductorModel.WATER_COOLED_COPPER
                else "#b87333"
            ),
        ),
        ("dr_shld_thermal_inboard", THERMAL_SHIELD_COLOUR[colour_scheme - 1]),
        ("dr_vv_inboard", VESSEL_COLOUR[colour_scheme - 1]),
        ("dr_shld_inboard", VESSEL_COLOUR[colour_scheme - 1]),
        ("dr_blkt_inboard", BLANKET_COLOUR[colour_scheme - 1]),
        ("dr_fw_inboard", FIRSTWALL_COLOUR[colour_scheme - 1]),
        ("dr_fw_outboard", FIRSTWALL_COLOUR[colour_scheme - 1]),
        ("dr_blkt_outboard", BLANKET_COLOUR[colour_scheme - 1]),
        ("dr_shld_outboard", SHIELD_COLOUR[colour_scheme - 1]),
        ("dr_vv_outboard", VESSEL_COLOUR[colour_scheme - 1]),
        ("dr_shld_thermal_outboard", THERMAL_SHIELD_COLOUR[colour_scheme - 1]),
    ]:
        r2, r1 = cumulative_radial_build2(v, mfile, scan)
        arc_fill(axis, r1, r2, color=colours, theta2=full_angle + 1)

    arc_fill(
        axis,
        rmajor - rminor,
        rmajor + rminor,
        color=PLASMA_COLOUR[colour_scheme - 1],
        theta2=full_angle + 1,
    )

    arc_fill(
        axis,
        r_cryostat_inboard,
        r_cryostat_inboard + dr_cryostat,
        color=CRYOSTAT_COLOUR[colour_scheme - 1],
        theta2=full_angle + 1,
    )

    # Segment the TF coil inboard
    # Calculate centrelines
    spacing = 2 * np.pi / n_tf_coils
    coil_indices = np.arange(int(n_tf_coils))

    r1, _ = cumulative_radial_build2("dr_cs_tf_gap", mfile, scan)
    r2, _ = cumulative_radial_build2("dr_tf_inboard", mfile, scan)
    r4, r3 = cumulative_radial_build2("dr_tf_outboard", mfile, scan)

    # Coil width
    w = r2 * np.tan(spacing / 2)
    for ang in (coil_indices * spacing) - spacing / 2:
        axis.plot(
            [r1 * np.cos(ang), r2 * np.cos(ang)],
            [r1 * np.sin(ang), r2 * np.sin(ang)],
            color="black",
        )

    for item in coil_indices:
        # Neutral beam shielding
        TF_outboard(
            axis,
            item,
            n_tf_coils=n_tf_coils,
            r3=r3,
            r4=r4,
            w=w + dx_beam_shield,
            facecolor=NBSHIELD_COLOUR[colour_scheme - 1],
        )
        # Overlay TF coil segments
        TF_outboard(
            axis,
            item,
            n_tf_coils=n_tf_coils,
            r3=r3,
            r4=r4,
            w=w,
            facecolor=(
                TFC_COLOUR[colour_scheme - 1]
                if TFConductorModel(mfile.get("i_tf_sup", scan=scan))
                != TFConductorModel.WATER_COOLED_COPPER
                else "#b87333"
            ),
        )

    i_hcd_primary = mfile.get("i_hcd_primary", scan=scan)
    if CurrentDriveModel(i_hcd_primary).method == CurrentDriveMethodType.NEUTRAL_BEAM:
        # Neutral beam geometry. See docs for diagram.
        a = w + dx_beam_shield
        b = dr_tf_outboard
        d = r3
        e = np.sqrt(a**2 + (d + b) ** 2)

        # Beam edges from centreline
        half_duct = 0.5 * dx_beam_duct
        r_beam_inner = radius_beam_tangency - half_duct
        r_beam_outer = radius_beam_tangency + half_duct

        def calc_xy(rt, e=e):
            arg = np.clip(rt / e, -1.0, 1.0)
            beta = np.arccos(arg)
            x = rt * np.cos(beta)
            y = rt * np.sin(beta)
            return x, y

        # Tangency points
        x_beam_inner, y_beam_inner = calc_xy(r_beam_inner)
        x_beam_outer, y_beam_outer = calc_xy(r_beam_outer)

        # TF-side positions (beam sits inside shield)
        x0 = r4
        y0_beam_inner = w + dx_beam_shield
        y0_beam_outer = y0_beam_inner + dx_beam_duct

        # Centreline tangency point
        x_beam_centre, y_beam_centre = calc_xy(radius_beam_tangency)
        y0_beam_centre = y0_beam_inner + 0.5 * dx_beam_duct

        # Draw beam duct boundaries
        axis.plot(
            [x_beam_inner, x0],
            [y_beam_inner, y0_beam_inner],
            linestyle="dotted",
            color="black",
        )
        axis.plot(
            [x_beam_outer, x0],
            [y_beam_outer, y0_beam_outer],
            linestyle="dotted",
            color="black",
        )
        # Draw beam centreline
        axis.plot(
            [x_beam_centre, x0],
            [y_beam_centre, y0_beam_centre],
            linestyle="--",
            color="black",
            linewidth=1.5,
        )

    # Draw dividing lines in the blanket (inboard modules, toroidal direction)
    n_blkt_inboard_modules_toroidal = mfile.get(
        "n_blkt_inboard_modules_toroidal", scan=scan
    )
    if n_blkt_inboard_modules_toroidal > 1:
        # Calculate the angular spacing for each module
        spacing = full_angle / (n_blkt_inboard_modules_toroidal)
        r1, _ = cumulative_radial_build2("dr_shld_inboard", mfile, scan)
        r2, _ = cumulative_radial_build2("dr_blkt_inboard", mfile, scan)
        for i in range(int(n_blkt_inboard_modules_toroidal)):
            ang = i * spacing
            # Draw a line from r1 to r2 at angle ang
            axis.plot(
                [r1 * np.cos(ang), r2 * np.cos(ang)],
                [r1 * np.sin(ang), r2 * np.sin(ang)],
                color="black",
                linestyle="-",
                linewidth=1.5,
                zorder=100,
            )

    # Draw dividing lines in the blanket (outboard modules, toroidal direction)
    n_blkt_outboard_modules_toroidal = mfile.get(
        "n_blkt_outboard_modules_toroidal", scan=scan
    )
    if n_blkt_outboard_modules_toroidal > 1:
        # Calculate the angular spacing for each module
        spacing = full_angle / (n_blkt_outboard_modules_toroidal)
        r1, _ = cumulative_radial_build2("dr_fw_outboard", mfile, scan)
        r2, _ = cumulative_radial_build2("dr_blkt_outboard", mfile, scan)
        for i in range(int(n_blkt_outboard_modules_toroidal)):
            ang = i * spacing
            # Draw a line from r1 to r2 at angle ang
            axis.plot(
                [r1 * np.cos(ang), r2 * np.cos(ang)],
                [r1 * np.sin(ang), r2 * np.sin(ang)],
                color="black",
                linestyle="-",
                linewidth=1.5,
                zorder=100,
            )

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        axis.set_ylim(0, 20)
        axis.set_xlim(0, 20)

    # Adaptive ranges
    else:
        axis.set_ylim(0.0, axis.get_ylim()[1])
        axis.set_xlim(0.0, axis.get_xlim()[1])


def arc(axis: plt.Axes, r, theta1=0, theta2=rtangle, style="solid"):
    """Plots an arc.

    Parameters
    ----------
    axis :
        plot object
    r :
        radius
    theta1 :
        starting polar angle (Default value = 0)
    theta2 :
        finishing polar angle (Default value = rtangle)
    axis: plt.Axes :

    style :
         (Default value = "solid")
    """
    angs = np.linspace(theta1, theta2)
    xs = r * np.cos(angs)
    ys = r * np.sin(angs)
    axis.plot(xs, ys, linestyle=style, color="black", lw=0.2)


def arc_fill(axis: plt.Axes, r1, r2, color="pink", theta1=0, theta2=rtangle):
    """Fills the space between two quarter circles.

    Parameters
    ----------
    axis :
        plot object
    r1 :
        r2 radii to be filled
    axis: plt.Axes :

    r2 :

    color :
         (Default value = "pink")
    """
    angs = np.linspace(theta1, theta2, endpoint=True)
    xs1 = r1 * np.cos(angs)
    ys1 = r1 * np.sin(angs)
    angs = np.linspace(theta2, theta1, endpoint=True)
    xs2 = r2 * np.cos(angs)
    ys2 = r2 * np.sin(angs)
    verts = list(zip(xs1, ys1, strict=False))
    verts.extend(list(zip(xs2, ys2, strict=False)))
    path = mplPath(verts, closed=True)
    patch = patches.PathPatch(path, facecolor=color, lw=0)
    axis.add_patch(patch)


__all__ = ["arc", "arc_fill", "toroidal_cross_section"]
