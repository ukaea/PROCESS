"""Magnets functions for PROCESS summary plots."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.patches import Circle, Rectangle
from matplotlib.path import Path as mplPath

from process.core.io.plot.summary.common import (
    box_style,
)
from process.core.io.plot.summary.constants import (
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    rtangle,
    rtangle2,
)
from process.core.io.plot.summary.magnets.cables import (
    plot_hts_tape_geometry,
)
from process.core.io.plot.summary.rendering import (
    draw_annotation,
    draw_text,
)
from process.data_structure.superconducting_tf_coil_variables import (
    TFWPIntegerTurnType,
)
from process.models.geometry.tfcoil import (
    tfcoil_geometry_d_shape,
    tfcoil_geometry_rectangular_shape,
)
from process.models.superconductors import SuperconductorModel
from process.models.tfcoil import quench
from process.models.tfcoil.base import TFCoilShapeModel, TFPlasmaCaseType

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def TF_outboard(axis: plt.Axes, item, n_tf_coils, r3, r4, w, facecolor):
    """Plot outboard TF coils"""
    spacing = 2 * np.pi / n_tf_coils
    ang = item * spacing
    dx = w * np.sin(ang)
    dy = w * np.cos(ang)
    x1 = r3 * np.cos(ang) + dx
    y1 = r3 * np.sin(ang) - dy
    x2 = r4 * np.cos(ang) + dx
    y2 = r4 * np.sin(ang) - dy
    x3 = r4 * np.cos(ang) - dx
    y3 = r4 * np.sin(ang) + dy
    x4 = r3 * np.cos(ang) - dx
    y4 = r3 * np.sin(ang) + dy
    verts = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]
    path = mplPath(verts, closed=True)
    patch = patches.PathPatch(path, facecolor=facecolor, lw=0)
    axis.add_patch(patch)


def plot_tf_coils(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    mirror_negative_x: bool = False,
):
    """Function to plot TF coils

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    # Arc points
    # MDK Only 4 points now required for elliptical arcs
    x1 = mfile.get("r_tf_arc(1)", scan=scan)
    y1 = mfile.get("z_tf_arc(1)", scan=scan)
    x2 = mfile.get("r_tf_arc(2)", scan=scan)
    y2 = mfile.get("z_tf_arc(2)", scan=scan)
    x3 = mfile.get("r_tf_arc(3)", scan=scan)
    y3 = mfile.get("z_tf_arc(3)", scan=scan)
    x4 = mfile.get("r_tf_arc(4)", scan=scan)
    y4 = mfile.get("z_tf_arc(4)", scan=scan)
    x5 = mfile.get("r_tf_arc(5)", scan=scan)
    y5 = mfile.get("z_tf_arc(5)", scan=scan)

    dr_tf_inboard = mfile.get("dr_tf_inboard", scan=scan)
    dr_tf_outboard = mfile.get("dr_tf_outboard", scan=scan)
    dr_shld_thermal_inboard = mfile.get("dr_shld_thermal_inboard", scan=scan)
    dr_shld_thermal_outboard = mfile.get("dr_shld_thermal_outboard", scan=scan)
    dr_tf_shld_gap = mfile.get("dr_tf_shld_gap", scan=scan)
    if y3 != 0:
        print("TF coil geometry: The value of z_tf_arc(3) is not zero, but should be.")

    if dr_shld_thermal_inboard != dr_shld_thermal_outboard:
        print(
            "dr_shld_thermal_inboard and dr_shld_thermal_outboard are"
            " different. Using dr_shld_thermal_inboardfor the poloidal plot of"
            " the thermal shield."
        )

    for offset, colour in (
        (
            dr_shld_thermal_inboard + dr_tf_shld_gap,
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
        ),
        (dr_tf_shld_gap, "white"),
        (
            0.0,
            (
                TFC_COLOUR[colour_scheme - 1]
                if mfile.get("i_tf_sup", scan=scan) != 0
                else "#b87333"
            ),
        ),
    ):
        # Check for TF coil shape
        if "i_tf_shape" in mfile.data:
            i_tf_shape = int(mfile.get("i_tf_shape", scan=scan))
        else:
            i_tf_shape = 1

        if i_tf_shape == TFCoilShapeModel.PICTURE_FRAME:
            rects = tfcoil_geometry_rectangular_shape(
                x1=x1,
                x2=x2,
                x4=x4,
                x5=x5,
                y1=y1,
                y2=y2,
                y4=y4,
                y5=y5,
                dr_tf_inboard=dr_tf_inboard,
                dr_tf_outboard=dr_tf_outboard,
                offset_in=offset,
            )

        else:
            rects, verts = tfcoil_geometry_d_shape(
                x1=x1,
                x2=x2,
                x3=x3,
                x4=x4,
                x5=x5,
                y1=y1,
                y2=y2,
                y4=y4,
                y5=y5,
                dr_tf_inboard=dr_tf_inboard,
                rtangle=rtangle,
                rtangle2=rtangle2,
                offset_in=offset,
            )

            for vert in verts:
                # Mirror vertices if needed
                mirrored_vert = [[x_scale * point[0], point[1]] for point in vert]
                path = mplPath(mirrored_vert, closed=True)
                patch = patches.PathPatch(path, facecolor=colour, lw=0)
                axis.add_patch(patch)

        for rec in rects:
            axis.add_patch(
                patches.Rectangle(
                    xy=(x_scale * rec.anchor_x, rec.anchor_z),
                    width=x_scale * rec.width,
                    height=rec.height,
                    facecolor=colour,
                )
            )


def plot_superconducting_tf_wp(axis: plt.Axes, mfile: MFile, scan: int, fig):
    """Plots inboard TF coil and winding pack.

    Parameters
    ----------
    axis : matplotlib.axes object
        Axis object to plot to.
    mfile : MFILE data object
        Object containing data for the plot.
    scan : int
        Scan number to use.
    """
    # Import the TF variables
    r_tf_inboard_in = mfile.get("r_tf_inboard_in", scan=scan)
    r_tf_inboard_out = mfile.get("r_tf_inboard_out", scan=scan)
    dx_tf_wp_primary_toroidal = mfile.get("dx_tf_wp_primary_toroidal", scan=scan)
    dx_tf_side_case_peak = mfile.get("dx_tf_side_case_peak", scan=scan)
    dx_tf_wp_secondary_toroidal = mfile.get("dx_tf_wp_secondary_toroidal", scan=scan)
    dr_tf_wp_with_insulation = mfile.get("dr_tf_wp_with_insulation", scan=scan)
    r_tf_wp_inboard_inner = mfile.get("r_tf_wp_inboard_inner", scan=scan)
    dx_tf_wp_insulation = mfile.get("dx_tf_wp_insulation", scan=scan)
    n_tf_coil_turns = round(mfile.get("n_tf_coil_turns", scan=scan))
    i_tf_wp_geom = round(mfile.get("i_tf_wp_geom", scan=scan))
    i_tf_sup = round(mfile.get("i_tf_sup", scan=scan))
    i_tf_case_geom = mfile.get("i_tf_case_geom", scan=scan)
    i_tf_turns_integer = mfile.get("i_tf_turns_integer", scan=scan)
    b_tf_inboard_peak_symmetric = mfile.get("b_tf_inboard_peak_symmetric", scan=scan)
    b_tf_inboard_peak_with_ripple = mfile.get("b_tf_inboard_peak_with_ripple", scan=scan)
    f_b_tf_inboard_peak_ripple_symmetric = mfile.get(
        "f_b_tf_inboard_peak_ripple_symmetric", scan=scan
    )
    r_b_tf_inboard_peak = mfile.get("r_b_tf_inboard_peak", scan=scan)
    dx_tf_wp_insertion_gap = mfile.get("dx_tf_wp_insertion_gap", scan=scan)
    r_tf_wp_inboard_outer = mfile.get("r_tf_wp_inboard_outer", scan=scan)
    r_tf_wp_inboard_centre = mfile.get("r_tf_wp_inboard_centre", scan=scan)

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        turn_layers = mfile.get("n_tf_wp_layers", scan=scan)
        turn_pancakes = mfile.get("n_tf_wp_pancakes", scan=scan)

    # Superconducting coil check
    if i_tf_sup == 1:
        axis.add_patch(
            Circle(
                (0, 0),
                r_tf_inboard_in,
                facecolor="none",
                edgecolor="black",
                linestyle="--",
            ),
        )

        if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
            axis.add_patch(
                Circle(
                    (0, 0),
                    r_tf_inboard_out,
                    facecolor="none",
                    edgecolor="black",
                    linestyle="--",
                ),
            )

        # Equations for plotting the TF case
        rad_tf_coil_inboard_toroidal_half = mfile.get(
            "rad_tf_coil_inboard_toroidal_half", scan=scan
        )

        # X points for inboard case curve
        x11 = r_tf_inboard_in * np.cos(
            np.linspace(
                rad_tf_coil_inboard_toroidal_half,
                -rad_tf_coil_inboard_toroidal_half,
                256,
                endpoint=True,
            )
        )
        # Y points for inboard case curve
        y11 = r_tf_inboard_in * np.sin(
            np.linspace(
                rad_tf_coil_inboard_toroidal_half,
                -rad_tf_coil_inboard_toroidal_half,
                256,
                endpoint=True,
            )
        )
        # Check for plasma side case type
        if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
            # Rounded case

            # X points for outboard case curve
            x12 = r_tf_inboard_out * np.cos(
                np.linspace(
                    rad_tf_coil_inboard_toroidal_half,
                    -rad_tf_coil_inboard_toroidal_half,
                    256,
                    endpoint=True,
                )
            )

        elif i_tf_case_geom == TFPlasmaCaseType.STRAIGHT:
            # Flat case

            # X points for outboard case
            x12 = np.full(256, r_tf_inboard_out)
        else:
            raise NotImplementedError("i_tf_case_geom must be 0 or 1")

        # Y points for outboard case
        y12 = r_tf_inboard_out * np.sin(
            np.linspace(
                rad_tf_coil_inboard_toroidal_half,
                -rad_tf_coil_inboard_toroidal_half,
                256,
                endpoint=True,
            )
        )

        # Cordinates of the top and bottom of case curves,
        # used to plot the lines connecting the inside and outside of the case
        y13 = [y11[0], y12[0]]
        x13 = [x11[0], x12[0]]
        y14 = [y11[-1], y12[-1]]
        x14 = [x11[-1], x12[-1]]

        # Plot the case outline
        axis.plot(x11, y11, color="black")
        axis.plot(x12, y12, color="black")
        axis.plot(x13, y13, color="black")
        axis.plot(x14, y14, color="black")

        # Fill in the case segemnts

        # Upper main
        if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                    (r_tf_inboard_out * np.cos(rad_tf_coil_inboard_toroidal_half)),
                ],
                y13,
                color="grey",
                alpha=0.25,
            )
            # Lower main
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                    (r_tf_inboard_out * np.cos(rad_tf_coil_inboard_toroidal_half)),
                ],
                y14,
                color="grey",
                alpha=0.25,
            )
            axis.fill_between(
                x12,
                y12,
                color="grey",
                alpha=0.25,
            )
        elif i_tf_case_geom == TFPlasmaCaseType.STRAIGHT:
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                    (r_tf_inboard_out),
                ],
                y13,
                color="grey",
                alpha=0.25,
            )
            # Lower main
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                    (r_tf_inboard_out),
                ],
                y14,
                color="grey",
                alpha=0.25,
            )

        # Removes ovelapping colours on inner nose case
        axis.fill_between(
            x11,
            y11,
            color="white",
            alpha=1.0,
        )

        # Centre line for relative reference
        axis.axhline(y=0.0, color="r", linestyle="--", linewidth=0.25)

        # ================================================================

        # Plot the rectangular WP
        if i_tf_wp_geom == 0:
            if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
                long_turns = round(turn_layers)
                short_turns = round(turn_pancakes)
            else:
                wp_side_ratio = (
                    dr_tf_wp_with_insulation
                    - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap))
                ) / (
                    dx_tf_wp_primary_toroidal
                    - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap))
                )  # row to height
                side_unit = n_tf_coil_turns / wp_side_ratio
                root_turns = round(np.sqrt(side_unit), 1)
                long_turns = round(root_turns * wp_side_ratio)
                short_turns = round(root_turns)

            # Plots the surrounding insualtion
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_inner,
                        -(0.5 * dx_tf_wp_primary_toroidal),
                    ),
                    dr_tf_wp_with_insulation,
                    dx_tf_wp_primary_toroidal,
                    color="darkgreen",
                ),
            )
            # Plots the WP inside the insulation
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_inner
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                        -(0.5 * dx_tf_wp_primary_toroidal)
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                    ),
                    (
                        dr_tf_wp_with_insulation
                        - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap))
                    ),
                    (
                        dx_tf_wp_primary_toroidal
                        - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap))
                    ),
                    color="blue",
                )
            )
            # Dvides the WP up into the turn segments
            for i in range(1, long_turns):
                axis.plot(
                    [
                        (
                            r_tf_wp_inboard_inner
                            + dx_tf_wp_insulation
                            + dx_tf_wp_insertion_gap
                        )
                        + i
                        * (
                            (
                                dr_tf_wp_with_insulation
                                - dx_tf_wp_insulation
                                - dx_tf_wp_insertion_gap
                            )
                            / long_turns
                        ),
                        (
                            r_tf_wp_inboard_inner
                            + dx_tf_wp_insulation
                            + dx_tf_wp_insertion_gap
                        )
                        + i
                        * (
                            (
                                dr_tf_wp_with_insulation
                                - dx_tf_wp_insulation
                                - dx_tf_wp_insertion_gap
                            )
                            / long_turns
                        ),
                    ],
                    [
                        -0.5 * dx_tf_wp_primary_toroidal
                        + (dx_tf_wp_insulation + dx_tf_wp_insertion_gap),
                        0.5 * dx_tf_wp_primary_toroidal
                        - (dx_tf_wp_insulation + dx_tf_wp_insertion_gap),
                    ],
                    color="white",
                    linewidth="0.25",
                    linestyle="dashed",
                )

            for i in range(1, short_turns):
                axis.plot(
                    [
                        (
                            r_tf_wp_inboard_inner
                            + dx_tf_wp_insulation
                            + dx_tf_wp_insertion_gap
                        ),
                        (
                            r_tf_wp_inboard_outer
                            - dx_tf_wp_insulation
                            - dx_tf_wp_insertion_gap
                        ),
                    ],
                    [
                        (
                            -0.5 * dx_tf_wp_primary_toroidal
                            + dx_tf_wp_insulation
                            + dx_tf_wp_insertion_gap
                        )
                        + (
                            i
                            * (
                                dx_tf_wp_primary_toroidal
                                - dx_tf_wp_insulation
                                - dx_tf_wp_insertion_gap
                            )
                            / short_turns
                        ),
                        (
                            -0.5 * dx_tf_wp_primary_toroidal
                            + dx_tf_wp_insulation
                            + dx_tf_wp_insertion_gap
                        )
                        + (
                            i
                            * (
                                dx_tf_wp_primary_toroidal
                                - dx_tf_wp_insulation
                                - dx_tf_wp_insertion_gap
                            )
                            / short_turns
                        ),
                    ],
                    color="white",
                    linewidth="0.25",
                    linestyle="dashed",
                )

        # ================================================================

        # Plot the double rectangle winding pack
        if i_tf_wp_geom == 1:
            # Inner WP insulation
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_inner,
                        -(0.5 * dx_tf_wp_secondary_toroidal),
                    ),
                    (dr_tf_wp_with_insulation / 2) + (dx_tf_wp_insulation),
                    dx_tf_wp_secondary_toroidal,
                    color="darkgreen",
                ),
            )

            # Outer WP insulation
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_centre,
                        -(0.5 * dx_tf_wp_primary_toroidal),
                    ),
                    (dr_tf_wp_with_insulation / 2),
                    dx_tf_wp_primary_toroidal,
                    color="darkgreen",
                ),
            )

            # Outer WP
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_centre
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                        -(0.5 * dx_tf_wp_primary_toroidal)
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                    ),
                    (dr_tf_wp_with_insulation / 2)
                    - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)),
                    dx_tf_wp_primary_toroidal
                    - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)),
                    color="blue",
                ),
            )
            # Inner WP
            axis.add_patch(
                Rectangle(
                    (
                        r_tf_wp_inboard_inner
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                        -(0.5 * dx_tf_wp_secondary_toroidal)
                        + dx_tf_wp_insulation
                        + dx_tf_wp_insertion_gap,
                    ),
                    (dr_tf_wp_with_insulation / 2),
                    dx_tf_wp_secondary_toroidal
                    - (2 * (dx_tf_wp_insulation + dx_tf_wp_insertion_gap)),
                    color="blue",
                ),
            )

        # ================================================================

        # Trapezium WP
        if i_tf_wp_geom == 2:
            # WP insulation
            x = [
                r_tf_wp_inboard_inner,
                r_tf_wp_inboard_inner,
                r_tf_wp_inboard_outer,
                r_tf_wp_inboard_outer,
            ]
            y = [
                (-0.5 * dx_tf_wp_secondary_toroidal),
                (0.5 * dx_tf_wp_secondary_toroidal),
                (0.5 * dx_tf_wp_primary_toroidal),
                (-0.5 * dx_tf_wp_primary_toroidal),
            ]
            axis.add_patch(
                patches.Polygon(
                    xy=list(zip(x, y, strict=False)),
                    color="darkgreen",
                )
            )

            # WP
            x = [
                r_tf_wp_inboard_inner + dx_tf_wp_insulation + dx_tf_wp_insertion_gap,
                r_tf_wp_inboard_inner + dx_tf_wp_insulation + dx_tf_wp_insertion_gap,
                (r_tf_wp_inboard_outer - dx_tf_wp_insulation - dx_tf_wp_insertion_gap),
                (r_tf_wp_inboard_outer - dx_tf_wp_insulation - dx_tf_wp_insertion_gap),
            ]
            y = [
                (
                    -0.5 * dx_tf_wp_secondary_toroidal
                    + dx_tf_wp_insulation
                    + dx_tf_wp_insertion_gap
                ),
                (
                    0.5 * dx_tf_wp_secondary_toroidal
                    - dx_tf_wp_insulation
                    - dx_tf_wp_insertion_gap
                ),
                (
                    0.5 * dx_tf_wp_primary_toroidal
                    - dx_tf_wp_insulation
                    - dx_tf_wp_insertion_gap
                ),
                (
                    -0.5 * dx_tf_wp_primary_toroidal
                    + dx_tf_wp_insulation
                    + dx_tf_wp_insertion_gap
                ),
            ]
            axis.add_patch(
                patches.Polygon(
                    xy=list(zip(x, y, strict=False)),
                    color="blue",
                )
            )

        # Plot a dot for the location of the peak field
        axis.plot(
            r_b_tf_inboard_peak,
            0,
            marker="o",
            color="red",
            label=(
                "Peak axisymmetric field:"
                f" {b_tf_inboard_peak_symmetric:.3f} T\n"
                "Peak non-axisymmetric field with ripple: "
                f"{b_tf_inboard_peak_with_ripple:.3f} T\n"
                "$\\frac{B_{\\text{axisymmetric}}}{B_{\\text{non-axisymmetric}}}$: "
                f"{f_b_tf_inboard_peak_ripple_symmetric:.3f}\n"
                f"$r_{{\\text{{peak}}}}$={r_b_tf_inboard_peak:.3f} m"
            ),
        )

        # Plot a horizontal line at y = dx_tf_wp_inner_toroidal
        axis.axhline(
            y=dx_tf_wp_secondary_toroidal / 2,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        # Plot a horizontal line at y = dx_tf_wp_inner_toroidal
        axis.axhline(
            y=-dx_tf_wp_secondary_toroidal / 2,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axhline(
            y=dx_tf_wp_primary_toroidal / 2,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axhline(
            y=-dx_tf_wp_primary_toroidal / 2,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        # Max toroidal width including side case
        axis.axhline(
            y=(dx_tf_wp_primary_toroidal / 2) + dx_tf_side_case_peak,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )

        axis.axhline(
            y=-(dx_tf_wp_primary_toroidal / 2) - dx_tf_side_case_peak,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )

        axis.axvline(
            x=r_tf_inboard_in,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axvline(
            x=r_tf_wp_inboard_inner,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axvline(
            x=r_tf_wp_inboard_outer,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axvline(
            x=r_tf_wp_inboard_centre,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )
        axis.axvline(
            x=r_tf_inboard_out,
            color="black",
            linestyle="--",
            linewidth=0.6,
            alpha=0.5,
        )

        # Add info about the steel casing surrounding the WP
        textstr_casing = (
            "$\\mathbf{Casing:}$\n\nCoil half angle:"
            f" {mfile.get('rad_tf_coil_inboard_toroidal_half', scan=scan):.3f}"
            " radians\n\n$\\text{Full Coil Case:}$\n$r_{start}"
            " \\rightarrow r_{end}$:"
            f" {mfile.get('r_tf_inboard_in', scan=scan):.3f} $\\rightarrow$"
            f" {mfile.get('r_tf_inboard_out', scan=scan):.3f} m\n$\\Delta r$:"
            f" {mfile.get('dr_tf_inboard', scan=scan):.3f} m\nArea of casing"
            f" around WP: {mfile.get('a_tf_coil_inboard_case', scan=scan):.3f}"
            " $\\mathrm{m}^2$\n\n$\\text{Nose Case:}$\n$r_{start}"
            " \\rightarrow r_{end}$:"
            f" {mfile.get('r_tf_inboard_in', scan=scan):.3f} $\\rightarrow$"
            f" {mfile.get('r_tf_wp_inboard_inner', scan=scan):.3f} m\n$\\Delta"
            f" r$: {mfile.get('dr_tf_nose_case', scan=scan):.3f} m\n$A$:"
            f" {mfile.get('a_tf_coil_nose_case', scan=scan):.3f}"
            " $\\mathrm{m}^2$\n\n$\\text{Plasma Case:}$\n$r_{start}"
            " \\rightarrow r_{end}$:"
            f" {mfile.get('r_tf_wp_inboard_outer', scan=scan):.3f}"
            f" $\\rightarrow$ {mfile.get('r_tf_inboard_out', scan=scan):.3f}"
            f" m\n$\\Delta r$: {mfile.get('dr_tf_plasma_case', scan=scan):.3f}"
            f" m\n$A$: {mfile.get('a_tf_plasma_case', scan=scan):.3f}"
            " $\\mathrm{m}^2$\n\n$\\text{Side Case:}$\nMinimum $\\Delta"
            f" r$: {mfile.get('dx_tf_side_case_min', scan=scan):.3f}"
            " m\nAverage $\\Delta r$:"
            f" {mfile.get('dx_tf_side_case_average', scan=scan):.3f} m\nMax"
            " $\\Delta r$:"
            f" {mfile.get('dx_tf_side_case_peak', scan=scan):.3f} m"
        )
        draw_text(
            axis,
            0.55,
            0.975,
            textstr_casing,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="left",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "grey",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add info about the steel casing surrounding the WP
        textstr_wp_insulation = (
            "$\\mathbf{Ground \\ Insulation:}$\n\nArea of insulation around"
            f" WP: {mfile.get('a_tf_wp_ground_insulation', scan=scan):.3f}"
            " $\\mathrm{m}^2$\n$\\Delta r$:"
            f" {mfile.get('dx_tf_wp_insulation', scan=scan):.4f} m\n\nWP"
            " Insertion Gap:\n$\\Delta r$:"
            f" {mfile.get('dx_tf_wp_insertion_gap', scan=scan):.4f} m"
        )
        draw_text(
            axis,
            0.55,
            0.575,
            textstr_wp_insulation,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="left",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "green",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add info about the Winding Pack
        textstr_wp = (
            "$\\mathbf{Winding \\  Pack:}$\n\n$N_{\\text{turns}}$:"
            f" {int(mfile.get('n_tf_coil_turns', scan=scan))}"
            " turns\n$r_{start} \\rightarrow r_{end}$:"
            f" {mfile.get('r_tf_wp_inboard_inner', scan=scan):.3f}"
            " $\\rightarrow$"
            f" {mfile.get('r_tf_wp_inboard_outer', scan=scan):.3f} m\n$\\Delta"
            f" r$: {mfile.get('dr_tf_wp_with_insulation', scan=scan):.3f}"
            " m\n\n$A$, with insulation:"
            f" {mfile.get('a_tf_wp_with_insulation', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n$A$, no insulation:"
            f" {mfile.get('a_tf_wp_no_insulation', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n$A$, total turn insulation:"
            f" {mfile.get('a_tf_coil_wp_turn_insulation', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n$A$, total turn steel:"
            f" {mfile.get('a_tf_wp_steel', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n$A$, total conductor:"
            f" {mfile.get('a_tf_wp_conductor', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n$A$, total non-cooling void:"
            f" {mfile.get('a_tf_wp_extra_void', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n\nPrimary WP:\n$\\Delta x$:"
            f" {mfile.get('dx_tf_wp_primary_toroidal', scan=scan):.4f}"
            " m\n\nSecondary WP:\n$\\Delta x$:"
            f" {mfile.get('dx_tf_wp_secondary_toroidal', scan=scan):.4f}"
            " m\n\n$J$ no insulation:"
            f" {mfile.get('j_tf_wp', scan=scan) / 1e6:.4f} MA/m$^2$"
        )

        draw_text(
            axis,
            0.775,
            0.95,
            textstr_wp,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="left",
            color="white",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "blue",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add info about the Winding Pack
        textstr_general_info = (
            "$\\mathbf{General \\ info:}$\n\n$N_{\\text{TF,coil}}$:"
            f" {mfile.get('n_tf_coils', scan=scan)}\nSelf inductance of single"
            f" coil: {mfile.get('ind_tf_coil', scan=scan) * 1e6:.4f}"
            " $\\mu$H\nStored energy of all coils:"
            f" {mfile.get('e_tf_magnetic_stored_total_gj', scan=scan):.4f}"
            " GJ\nStored energy of a single coil:"
            f" {mfile.get('e_tf_coil_magnetic_stored', scan=scan) / 1e9:.2f}"
            " GJ\nTotal area of steel in coil:"
            f" {mfile.get('a_tf_coil_inboard_steel', scan=scan):.4f}"
            " $\\mathrm{m}^2$\nTotal area fraction of steel:"
            f" {mfile.get('f_a_tf_coil_inboard_steel', scan=scan):.4f}\nTotal"
            " area fraction of insulation:"
            f" {mfile.get('f_a_tf_coil_inboard_insulation', scan=scan):.4f}\n$A$,"
            " all insulation in coil:"
            f" {mfile.get('a_tf_coil_inboard_insulation', scan=scan):.4f}"
            " $\\mathrm{m}^2$\n"
        )
        draw_text(
            axis,
            0.775,
            0.58,
            textstr_general_info,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="left",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "wheat",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        axis.minorticks_on()
        axis.set_xlim(r_tf_inboard_in * 0.8, r_tf_inboard_out * 1.1)
        axis.set_ylim((y14[-1] * 1.25), (-y14[-1] * 1.25))

        axis.set_title("Top-down view of inboard TF coil at midplane")
        axis.set_xlabel("Radial distance [m]")
        axis.set_ylabel("Toroidal distance [m]")
        axis.legend(loc="upper left")


def plot_resistive_tf_wp(axis: plt.Axes, mfile: MFile, scan: int, fig):
    """Plots inboard TF coil and winding pack.

    Parameters
    ----------
    axis : matplotlib.axes object
        Axis object to plot to.
    mfile : MFILE data object
        Object containing data for the plot.
    scan : int
        Scan number to use.
    """
    # Import the TF variables
    r_tf_inboard_in = mfile.get("r_tf_inboard_in", scan=scan)
    r_tf_inboard_out = mfile.get("r_tf_inboard_out", scan=scan)

    r_tf_wp_inboard_inner = mfile.get("r_tf_wp_inboard_inner", scan=scan)
    i_tf_case_geom = mfile.get("i_tf_case_geom", scan=scan)
    b_tf_inboard_peak_symmetric = mfile.get("b_tf_inboard_peak_symmetric", scan=scan)
    r_b_tf_inboard_peak = mfile.get("r_b_tf_inboard_peak", scan=scan)
    r_tf_wp_inboard_outer = mfile.get("r_tf_wp_inboard_outer", scan=scan)
    r_tf_wp_inboard_centre = mfile.get("r_tf_wp_inboard_centre", scan=scan)
    dx_tf_wp_insulation = mfile.get("dx_tf_wp_insulation", scan=scan)

    axis.add_patch(
        Circle(
            (0, 0),
            r_tf_inboard_in,
            facecolor="none",
            edgecolor="black",
            linestyle="--",
        ),
    )

    if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
        axis.add_patch(
            Circle(
                (0, 0),
                r_tf_inboard_out,
                facecolor="none",
                edgecolor="black",
                linestyle="--",
            ),
        )

    # Equations for plotting the TF case
    rad_tf_coil_inboard_toroidal_half = mfile.get(
        "rad_tf_coil_inboard_toroidal_half", scan=scan
    )

    # X points for inboard case curve
    x11 = r_tf_inboard_in * np.cos(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            256,
            endpoint=True,
        )
    )
    # Y points for inboard case curve
    y11 = r_tf_inboard_in * np.sin(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            256,
            endpoint=True,
        )
    )
    # Check for plasma side case type
    if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
        # Rounded case

        # X points for outboard case curve
        x12 = r_tf_inboard_out * np.cos(
            np.linspace(
                rad_tf_coil_inboard_toroidal_half,
                -rad_tf_coil_inboard_toroidal_half,
                256,
                endpoint=True,
            )
        )

    elif i_tf_case_geom == TFPlasmaCaseType.STRAIGHT:
        # Flat case

        # X points for outboard case
        x12 = np.full(256, r_tf_inboard_out)
    else:
        raise NotImplementedError("i_tf_case_geom must be 0 or 1")

    # Y points for outboard case
    y12 = r_tf_inboard_out * np.sin(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            256,
            endpoint=True,
        )
    )

    # Cordinates of the top and bottom of case curves,
    # used to plot the lines connecting the inside and outside of the case
    y13 = [y11[0], y12[0]]
    x13 = [x11[0], x12[0]]
    y14 = [y11[-1], y12[-1]]
    x14 = [x11[-1], x12[-1]]

    # Plot the case outline
    axis.plot(x11, y11, color="black")
    axis.plot(x12, y12, color="black")
    axis.plot(x13, y13, color="black")
    axis.plot(x14, y14, color="black")

    # Fill in the case segemnts

    # Upper main
    if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
        axis.fill_between(
            [
                (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_inboard_out * np.cos(rad_tf_coil_inboard_toroidal_half)),
            ],
            y13,
            color="grey",
            alpha=0.25,
        )
        # Lower main
        axis.fill_between(
            [
                (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_inboard_out * np.cos(rad_tf_coil_inboard_toroidal_half)),
            ],
            y14,
            color="grey",
            alpha=0.25,
        )
        axis.fill_between(
            x12,
            y12,
            color="grey",
            alpha=0.25,
        )
    elif i_tf_case_geom == TFPlasmaCaseType.STRAIGHT:
        axis.fill_between(
            [
                (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_inboard_out),
            ],
            y13,
            color="grey",
            alpha=0.25,
        )
        # Lower main
        axis.fill_between(
            [
                (r_tf_inboard_in * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_inboard_out),
            ],
            y14,
            color="grey",
            alpha=0.25,
        )

    # Removes ovelapping colours on inner nose case
    axis.fill_between(
        x11,
        y11,
        color="white",
        alpha=1.0,
    )

    # Centre line for relative reference
    axis.axhline(y=0.0, color="r", linestyle="--", linewidth=0.25)

    # ================================================================

    # Plot the WP insulation

    # X points for inboard insulation curve
    x11 = r_tf_wp_inboard_inner * np.cos(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            500,
            endpoint=True,
        )
    )
    # Y points for inboard insulation curve
    y11 = r_tf_wp_inboard_inner * np.sin(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            500,
            endpoint=True,
        )
    )

    # X points for outboard insulation curve
    x12 = r_tf_wp_inboard_outer * np.cos(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            500,
            endpoint=True,
        )
    )

    # Y points for outboard insulation curve
    y12 = r_tf_wp_inboard_outer * np.sin(
        np.linspace(
            rad_tf_coil_inboard_toroidal_half,
            -rad_tf_coil_inboard_toroidal_half,
            500,
            endpoint=True,
        )
    )

    # Cordinates of the top and bottom of WP insulation curves,
    y13 = [y11[0], y12[0]]
    x13 = [x11[0], x12[0]]
    y14 = [y11[-1], y12[-1]]
    x14 = [x11[-1], x12[-1]]

    # Plot the insualtion outline
    axis.plot(x11, y11, color="black")
    axis.plot(x12, y12, color="black")
    axis.plot(x13, y13, color="black")
    axis.plot(x14, y14, color="black")

    # Upper main
    if i_tf_case_geom == TFPlasmaCaseType.CIRCULAR:
        axis.fill_between(
            [
                (r_tf_wp_inboard_inner * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_wp_inboard_outer * np.cos(rad_tf_coil_inboard_toroidal_half)),
            ],
            y13,
            color="green",
        )
        # Lower main
        axis.fill_between(
            [
                (r_tf_wp_inboard_inner * np.cos(rad_tf_coil_inboard_toroidal_half)),
                (r_tf_wp_inboard_outer * np.cos(rad_tf_coil_inboard_toroidal_half)),
            ],
            y14,
            color="green",
        )
        axis.fill_between(
            x12,
            y12,
            color="green",
        )

    # ================================================================

    # Plot the WP

    # The winding pack should be inside the insulation, so subtract dx_tf_wp_insulation
    # from both the inner and outer radii.
    # The angular extent should also be reduced by the insulation thickness, i.e., the
    # winding pack does not extend all the way to the top/bottom.

    # Calculate the reduced angle for the winding pack (subtract insulation thickness in
    # arc length, convert to angle)
    # arc_length = r * angle => angle = arc_length / r
    # So, for both inner and outer radii, compute the angle offset due to insulation
    # thickness
    angle_offset_inner = (
        dx_tf_wp_insulation / r_tf_wp_inboard_inner if r_tf_wp_inboard_inner > 0 else 0
    )
    angle_offset_outer = (
        dx_tf_wp_insulation / r_tf_wp_inboard_outer if r_tf_wp_inboard_outer > 0 else 0
    )

    # Use the maximum angle offset to ensure the winding pack stays within the insulation
    angle_offset = max(angle_offset_inner, angle_offset_outer)

    # Define the angular range for the winding pack
    theta_start = rad_tf_coil_inboard_toroidal_half - angle_offset
    theta_end = -rad_tf_coil_inboard_toroidal_half + angle_offset
    theta_vals = np.linspace(theta_start, theta_end, 256, endpoint=True)

    # X and Y points for inboard and outboard winding pack curves
    x11 = (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.cos(theta_vals)
    y11 = (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.sin(theta_vals)
    x12 = (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.cos(theta_vals)
    y12 = (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.sin(theta_vals)

    # Cordinates of the top and bottom of WP curves,
    y13 = [y11[0], y12[0]]
    x13 = [x11[0], x12[0]]
    y14 = [y11[-1], y12[-1]]
    x14 = [x11[-1], x12[-1]]

    # Plot the winding pack outline
    axis.plot(x11, y11, color="black")
    axis.plot(x12, y12, color="black")
    axis.plot(x13, y13, color="black")
    axis.plot(x14, y14, color="black")

    # Choose color based on i_tf_sup: copper for resistive, aluminium for cryo
    # light steel blue (cryo aluminium) or copper color
    wp_color = "#b0c4de" if mfile.get("i_tf_sup", scan=scan) == 2 else "#b87333"

    axis.fill_between(
        [
            (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.cos(theta_vals[0]),
            (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.cos(theta_vals[0]),
        ],
        y13,
        color=wp_color,
    )
    # Lower main
    axis.fill_between(
        [
            (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.cos(theta_vals[-1]),
            (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.cos(theta_vals[-1]),
        ],
        y14,
        color=wp_color,
    )
    axis.fill_between(x12, y12, color=wp_color)

    # ================================================================

    # Divide the winding pack into toroidal segments based on n_tf_coil_turns
    n_turns = int(mfile.get("n_tf_coil_turns", scan=scan))
    if n_turns > 0:
        # Calculate the angular extent for each turn
        theta_start = rad_tf_coil_inboard_toroidal_half - angle_offset
        theta_end = -rad_tf_coil_inboard_toroidal_half + angle_offset
        theta_vals = np.linspace(theta_start, theta_end, 256, endpoint=True)

        # For each turn, plot a radial line at the corresponding angle
        turn_angles = np.linspace(theta_start, theta_end, n_turns + 1)
        for t in range(1, n_turns):
            angle = turn_angles[t]
            # Inner and outer points for this turn
            x_in = (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.cos(angle)
            y_in = (r_tf_wp_inboard_inner + dx_tf_wp_insulation) * np.sin(angle)
            x_out = (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.cos(angle)
            y_out = (r_tf_wp_inboard_outer - dx_tf_wp_insulation) * np.sin(angle)
            axis.plot(
                [x_in, x_out],
                [y_in, y_out],
                color="white",
                linewidth=0.5,
                linestyle="--",
            )

    # ================================================================

    # Plot a dot for the location of the peak field
    axis.plot(
        r_b_tf_inboard_peak,
        0,
        marker="o",
        color="red",
        label=(
            f"Peak Field: {b_tf_inboard_peak_symmetric:.2f}"
            f" T\nr={r_b_tf_inboard_peak:.3f} m"
        ),
    )

    x_kwargs = {
        "color": "black",
        "linestyle": "--",
        "linewidth": 0.6,
        "alpha": 0.5,
    }
    axis.axvline(x=r_tf_inboard_in, **x_kwargs)
    axis.axvline(x=r_tf_wp_inboard_inner, **x_kwargs)
    axis.axvline(x=r_tf_wp_inboard_outer, **x_kwargs)
    axis.axvline(x=r_tf_wp_inboard_centre, **x_kwargs)
    axis.axvline(x=r_tf_inboard_out, **x_kwargs)

    axis.minorticks_on()
    axis.set_xlim(0.0, r_tf_inboard_out * 1.1)
    axis.set_ylim((y14[-1] * 1.65), (-y14[-1] * 1.65))

    axis.set_title("Top-down view of inboard TF coil at midplane")
    axis.set_xlabel("Radial distance [m]")
    axis.set_ylabel("Toroidal distance [m]")
    axis.legend(loc="upper left")

    draw_text(
        axis,
        0.05,
        0.975,
        "*Turn insulation and cooling pipes not shown",
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
        transform=fig.transFigure,
    )


def plot_resistive_tf_info(axis: plt.Axes, mfile: MFile, scan: int, fig):
    """Plot info about the resistive TF coils"""
    # Add info about the steel casing surrounding the WP
    textstr_casing = (
        "$\\mathbf{Casing:}$\n\nCoil half angle:"
        f" {mfile.get('rad_tf_coil_inboard_toroidal_half', scan=scan):.3f}"
        " radians\n\n$\\text{Full Coil Case:}$\n$r_{start} \\rightarrow"
        f" r_{{end}}$: {mfile.get('r_tf_inboard_in', scan=scan):.3f}"
        f" $\\rightarrow$ {mfile.get('r_tf_inboard_out', scan=scan):.3f}"
        f" m\n$\\Delta r$: {mfile.get('dr_tf_inboard', scan=scan):.3f} m\nArea"
        " of casing around WP:"
        f" {mfile.get('a_tf_coil_inboard_case', scan=scan):.3f}"
        " $\\mathrm{m}^2$\n\n$\\text{Nose Case:}$\n$r_{start}"
        " \\rightarrow r_{end}$:"
        f" {mfile.get('r_tf_inboard_in', scan=scan):.3f} $\\rightarrow$"
        f" {mfile.get('r_tf_wp_inboard_inner', scan=scan):.3f} m\n$\\Delta r$:"
        f" {mfile.get('dr_tf_nose_case', scan=scan):.3f} m\n$A$:"
        f" {mfile.get('a_tf_coil_nose_case', scan=scan):.4f}"
        " $\\mathrm{m}^2$\n\n$\\text{Plasma Case:}$\n$r_{start}"
        " \\rightarrow r_{end}$:"
        f" {mfile.get('r_tf_wp_inboard_outer', scan=scan):.3f} $\\rightarrow$"
        f" {mfile.get('r_tf_inboard_out', scan=scan):.3f} m\n$\\Delta r$:"
        f" {mfile.get('dr_tf_plasma_case', scan=scan):.3f} m\n$A$:"
        f" {mfile.get('a_tf_plasma_case', scan=scan):.3f} $\\mathrm{{m}}^2$"
    )
    draw_text(
        axis,
        0.775,
        0.925,
        textstr_casing,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("grey"),
    )

    # Add info about the steel casing surrounding the WP
    textstr_wp_insulation = (
        "$\\mathbf{Insulation:}$\n\nArea of insulation around WP:"
        f" {mfile.get('a_tf_wp_ground_insulation', scan=scan):.3f}"
        " $\\mathrm{m}^2$\n$\\Delta r$:"
        f" {mfile.get('dx_tf_wp_insulation', scan=scan):.4f}"
        " m\n\n$\\text{Turn Insulation:}$\n$\\Delta r$:"
        f" {mfile.get('dx_tf_turn_insulation', scan=scan):.4f} m"
    )
    draw_text(
        axis,
        0.775,
        0.62,
        textstr_wp_insulation,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "green",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add info about the Winding Pack
    textstr_wp = (
        "$\\mathbf{Winding Pack:}$\n\n$N_{\\text{turns}}$:"
        f" {int(mfile.get('n_tf_coil_turns', scan=scan))} turns\n$r_{{start}}"
        " \\rightarrow r_{end}$:"
        f" {mfile.get('r_tf_wp_inboard_inner', scan=scan):.3f} $\\rightarrow$"
        f" {mfile.get('r_tf_wp_inboard_outer', scan=scan):.3f} m\n$\\Delta r$:"
        f" {mfile.get('dr_tf_wp_with_insulation', scan=scan):.3f} m\n$A$, with"
        f" insulation: {mfile.get('a_tf_wp_with_insulation', scan=scan):.3f}"
        " $\\mathrm{m}^2$\n$A$, no insulation:"
        f" {mfile.get('a_tf_wp_no_insulation', scan=scan):.3f}"
        " $\\mathrm{m}^2$\n\nCurrent per turn:"
        f" {mfile.get('c_tf_turn', scan=scan) / 1e3:.3f}"
        " $\\mathrm{kA}$\nResistive conductor per coil:"
        f" {mfile.get('a_res_tf_coil_conductor', scan=scan):.3f}"
        " $\\mathrm{m}^2$\nCoolant area void fraction per turn:"
        f" {mfile.get('fcoolcp', scan=scan):.3f}"
    )
    draw_text(
        axis,
        0.77,
        0.475,
        textstr_wp,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        color="white",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "blue",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add info about the Winding Pack
    textstr_general_info = (
        "$\\mathbf{General \\ info:}$\n\nSelf inductance:"
        f" {mfile.get('ind_tf_coil', scan=scan) * 1e6:.4f} $\\mu$H\nStored"
        " energy of all coils:"
        f" {mfile.get('e_tf_magnetic_stored_total_gj', scan=scan):.4f} GJ\n"
    )
    draw_text(
        axis,
        0.55,
        0.475,
        textstr_general_info,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add info about the Winding Pack
    textstr_cooling = (
        "$\\mathbf{Cooling \\ info:}$\n\nCoolant inlet temperature:"
        f" {mfile.get('temp_cp_coolant_inlet', scan=scan):.2f} K\nCoolant"
        f" temperature rise: {mfile.get('dtemp_cp_coolant', scan=scan):.2f}"
        " K\nCoolant velocity:"
        f" {mfile.get('vel_cp_coolant_midplane', scan=scan):.2f}"
        " $\\mathrm{ms^{-1}}$\n\nAverage CP temperature:"
        f" {mfile.get('temp_cp_average', scan=scan):.2f} K\nCP resistivity:"
        f" {mfile.get('rho_cp', scan=scan):.2e} $\\Omega \\mathrm{{m}}$\nLeg"
        f" resistivity: {mfile.get('rho_tf_leg', scan=scan):.2e} $\\Omega"
        " \\mathrm{m}$\nLeg resistance:"
        f" {mfile.get('res_tf_leg', scan=scan):.2e} $\\Omega$\nCP resistive"
        f" losses: {mfile.get('p_cp_resistive', scan=scan):,.2f}"
        " $\\mathrm{W}$\nLeg resistive losses:"
        f" {mfile.get('p_tf_leg_resistive', scan=scan):,.2f}"
        " $\\mathrm{W}$\nJoints resistive losses:"
        f" {mfile.get('p_tf_joints_resistive', scan=scan):,.2f}"
        " $\\mathrm{W}$\n"
    )
    draw_text(
        axis,
        0.55,
        0.35,
        textstr_cooling,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("wheat"),
    )


def plot_tf_cable_in_conduit_turn(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plots inboard TF coil CICC individual turn structure.

    Parameters
    ----------
    axis : matplotlib.axes object
        Axis object to plot to.
    mfile : MFILE data object
        Object containing data for the plot.
    scan : int
        Scan number to use.
    """

    def _pack_strands_rectangular_with_obstacles(
        cable_space_bounds,
        pipe_center,
        pipe_radius,
        strand_diameter,
        void_fraction,
        n_strands,
        axis,
        corner_radius,
        f_a_tf_turn_cable_copper,
    ):
        """Pack circular strands in rectangular space with cooling pipe obstacle

        Parameters
        ----------
        cable_space_bounds :

        pipe_center :

        pipe_radius :

        strand_diameter :

        void_fraction :

        n_strands :

        axis :

        corner_radius :

        f_a_tf_turn_cable_copper :

        """
        x, y, width, height = cable_space_bounds

        radius = strand_diameter / 2
        placed_strands = []
        attempts = 0

        pipe_x, pipe_y = pipe_center

        # Hexagonal packing parameters
        # Calculate the spacing between strand centers for the desired void fraction
        # For hexagonal packing, packing fraction = pi/(2*sqrt(3)) ~ 0.9069
        # To achieve a lower packing fraction (higher void fraction), increase spacing
        ideal_packing_fraction = np.pi / (2 * np.sqrt(3))
        target_packing_fraction = 1 - void_fraction
        spacing_factor = np.sqrt(ideal_packing_fraction / target_packing_fraction)
        strand_spacing = strand_diameter * spacing_factor

        # Number of rows and columns that fit in the cable space
        n_rows = int((height - 2 * radius) // (strand_spacing * np.sqrt(3) / 2))
        n_cols = int((width - 2 * radius) // strand_spacing)

        # Calculate the radius of the inner superconductor circle based on the copper
        # area fraction
        # Area_superconductor = (1 - f_a_tf_turn_cable_copper) * Area_strand
        # Area_strand = pi * radius^2
        # So, radius_superconductor = sqrt(1 - f_a_tf_turn_cable_copper) * radius
        radius_superconductor = np.sqrt(1 - f_a_tf_turn_cable_copper) * radius

        # Generate hexagonal grid positions
        for row in range(n_rows):
            y_pos = (y + radius + row * strand_spacing * np.sqrt(3) / 2) * 1.07
            x_offset = strand_spacing / 2 if row % 2 else 0
            for col in range(n_cols):
                candidate_x = (x + radius + col * strand_spacing + x_offset) * 1.05
                candidate_y = y_pos

                # Check if within bounds
                if candidate_x > x + width - radius or candidate_y > y + height - radius:
                    continue

                # Check collision with cooling pipe
                pipe_distance = np.sqrt(
                    (candidate_x - pipe_x) ** 2 + (candidate_y - pipe_y) ** 2
                )
                if pipe_distance < (pipe_radius + radius):
                    continue

                # Check collision with corners if rounded
                if corner_radius > 0:
                    corners = [
                        (x + corner_radius, y + corner_radius),  # bottom-left
                        (
                            x + width - corner_radius,
                            y + corner_radius,
                        ),  # bottom-right
                        (
                            x + width - corner_radius,
                            y + height - corner_radius,
                        ),  # top-right
                        (
                            x + corner_radius,
                            y + height - corner_radius,
                        ),  # top-left
                    ]
                    if (
                        (
                            candidate_x < corners[0][0]
                            and candidate_y < corners[0][1]
                            and np.sqrt(
                                (candidate_x - corners[0][0]) ** 2
                                + (candidate_y - corners[0][1]) ** 2
                            )
                            > corner_radius - radius
                        )
                        or (
                            candidate_x > corners[1][0]
                            and candidate_y < corners[1][1]
                            and np.sqrt(
                                (candidate_x - corners[1][0]) ** 2
                                + (candidate_y - corners[1][1]) ** 2
                            )
                            > corner_radius - radius
                        )
                        or (
                            candidate_x > corners[2][0]
                            and candidate_y > corners[2][1]
                            and np.sqrt(
                                (candidate_x - corners[2][0]) ** 2
                                + (candidate_y - corners[2][1]) ** 2
                            )
                            > corner_radius - radius
                        )
                        or (
                            candidate_x < corners[3][0]
                            and candidate_y > corners[3][1]
                            and np.sqrt(
                                (candidate_x - corners[3][0]) ** 2
                                + (candidate_y - corners[3][1]) ** 2
                            )
                            > corner_radius - radius
                        )
                    ):
                        continue

                # Check collision with existing strands
                collision = False
                for existing_x, existing_y in placed_strands:
                    distance = np.sqrt(
                        (candidate_x - existing_x) ** 2 + (candidate_y - existing_y) ** 2
                    )
                    if distance < strand_diameter:
                        collision = True
                        break

                if not collision:
                    placed_strands.append((candidate_x, candidate_y))
                    # Plot the strand
                    circle_copper_surrounding = Circle(
                        (candidate_x, candidate_y),
                        radius,
                        facecolor="#b87333",  # copper color
                        edgecolor="#8B4000",  # darker copper edge
                        linewidth=0.1,
                        alpha=0.8,
                    )
                    axis.add_patch(circle_copper_surrounding)

                    circle_central_conductor = Circle(
                        (candidate_x, candidate_y),
                        radius_superconductor,
                        facecolor="black",
                        linewidth=0.3,
                        alpha=0.5,
                    )
                    axis.add_patch(circle_central_conductor)

                if len(placed_strands) >= n_strands:
                    break
            if len(placed_strands) >= n_strands:
                break

        attempts = n_rows * n_cols

        return len(placed_strands), attempts

    # Import the TF turn variables then multiply into mm
    i_tf_turns_integer = mfile.get("i_tf_turns_integer", scan=scan)
    # If integer turns switch is on then the turns can have non square dimensions
    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        turn_width = mfile.get("dr_tf_turn", scan=scan)
        turn_height = mfile.get("dx_tf_turn", scan=scan)
        cable_space_width_radial = mfile.get("dr_tf_turn_cable_space", scan=scan)
        cable_space_width_toroidal = mfile.get("dx_tf_turn_cable_space", scan=scan)

    elif TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        turn_width = mfile.get("dx_tf_turn_general", scan=scan)
        cable_space_width = mfile.get("dx_tf_turn_cable_space_average", scan=scan)

    he_pipe_diameter = mfile.get("dia_tf_turn_coolant_channel", scan=scan)
    steel_thickness = mfile.get("dx_tf_turn_steel", scan=scan)
    insulation_thickness = mfile.get("dx_tf_turn_insulation", scan=scan)

    a_tf_turn_cable_space_no_void = mfile.get("a_tf_turn_cable_space_no_void", scan=scan)
    radius_tf_turn_cable_space_corners = mfile.get(
        "radius_tf_turn_cable_space_corners", scan=scan
    )

    a_tf_wp_coolant_channels = mfile.get("a_tf_wp_coolant_channels", scan=scan)

    f_a_tf_turn_cable_space_extra_void = mfile.get(
        "f_a_tf_turn_cable_space_extra_void", scan=scan
    )
    a_tf_turn_steel = mfile.get("a_tf_turn_steel", scan=scan)
    a_tf_turn_cable_space_effective = mfile.get(
        "a_tf_turn_cable_space_effective", scan=scan
    )

    # Plot the total turn shape
    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        axis.add_patch(
            Rectangle(
                (0, 0),
                turn_width,
                turn_width,
                facecolor="red",
                edgecolor="black",
            ),
        )
        # Plot the steel conduit
        axis.add_patch(
            Rectangle(
                (insulation_thickness, insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                facecolor="grey",
                edgecolor="black",
            ),
        )

        # Plot the cable space with rounded corners
        axis.add_patch(
            patches.FancyBboxPatch(
                (
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                boxstyle=patches.BoxStyle(
                    "Round",
                    pad=0,
                    rounding_size=radius_tf_turn_cable_space_corners,
                ),
                facecolor="royalblue",
                edgecolor="black",
            ),
        )

        # Plot dashed line around the cable space
        axis.add_patch(
            Rectangle(
                (
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                facecolor="none",
                edgecolor="black",
                linestyle="--",
                linewidth=1.2,
                alpha=0.5,
            ),
        )
        # Plot the coolant channel
        axis.add_patch(
            Circle(
                ((turn_width / 2), (turn_width / 2)),
                he_pipe_diameter / 2,
                facecolor="white",
                edgecolor="black",
            ),
        )

        # Cable strand packing parameters
        strand_diameter = mfile.get("dia_tf_turn_superconducting_cable", scan=scan)
        void_fraction = mfile.get("f_a_tf_turn_cable_space_extra_void", scan=scan)

        # Cable space bounds
        cable_bounds = [
            insulation_thickness + steel_thickness,
            insulation_thickness + steel_thickness,
            turn_width - 2 * (insulation_thickness + steel_thickness),
            turn_width - 2 * (insulation_thickness + steel_thickness),
        ]

        # Pack strands if significant void fraction
        if void_fraction > 0.001:
            _n_strands, _attempts = _pack_strands_rectangular_with_obstacles(
                cable_space_bounds=cable_bounds,
                pipe_center=(
                    turn_width / 2,
                    (
                        turn_width
                        if TFWPIntegerTurnType(i_tf_turns_integer)
                        == TFWPIntegerTurnType.NON_INTEGER
                        else turn_height
                    )
                    / 2,
                ),
                pipe_radius=he_pipe_diameter / 2,
                strand_diameter=strand_diameter,
                void_fraction=void_fraction,
                axis=axis,
                corner_radius=radius_tf_turn_cable_space_corners,
                n_strands=mfile.get("n_tf_turn_superconducting_cables", scan=scan),
                f_a_tf_turn_cable_copper=mfile.get(
                    "f_a_tf_turn_cable_copper", scan=scan
                ),
            )

        axis.set_xlim(-turn_width * 0.05, turn_width * 1.05)
        axis.set_ylim(-turn_width * 0.05, turn_width * 1.05)

    # Non square turns
    elif TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        axis.add_patch(
            Rectangle(
                (0, 0),
                turn_width,
                turn_height,
                facecolor="red",
                edgecolor="black",
            ),
        )

        # Plot the steel conduit
        axis.add_patch(
            Rectangle(
                (insulation_thickness, insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                (turn_height - 2 * insulation_thickness),
                facecolor="grey",
                edgecolor="black",
            ),
        )

        # Plot the cable space with rounded corners
        axis.add_patch(
            patches.FancyBboxPatch(
                (
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_height - 2 * (insulation_thickness + steel_thickness)),
                boxstyle=patches.BoxStyle(
                    "Round",
                    pad=0,
                    rounding_size=radius_tf_turn_cable_space_corners,
                ),
                facecolor="royalblue",
                edgecolor="black",
            ),
        )
        # Plot dashed line around the cable space
        axis.add_patch(
            Rectangle(
                (
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_height - 2 * (insulation_thickness + steel_thickness)),
                facecolor="none",
                edgecolor="black",
                linestyle="--",
                linewidth=1.0,
                alpha=0.5,
            ),
        )
        axis.add_patch(
            Circle(
                ((turn_width / 2), (turn_height / 2)),
                he_pipe_diameter / 2,
                facecolor="white",
                edgecolor="black",
            ),
        )

        # Cable space bounds
        cable_bounds = [
            insulation_thickness + steel_thickness,
            insulation_thickness + steel_thickness,
            turn_width - 2 * (insulation_thickness + steel_thickness),
            turn_height - 2 * (insulation_thickness + steel_thickness),
        ]

        # Cable strand packing parameters
        strand_diameter = mfile.get("dia_tf_turn_superconducting_cable", scan=scan)
        void_fraction = mfile.get("f_a_tf_turn_cable_space_extra_void", scan=scan)

        # Pack strands if significant void fraction
        if void_fraction > 0.001:
            _, _ = _pack_strands_rectangular_with_obstacles(
                cable_space_bounds=cable_bounds,
                pipe_center=(
                    turn_width / 2,
                    turn_height / 2,
                ),
                pipe_radius=he_pipe_diameter / 2,
                strand_diameter=strand_diameter,
                void_fraction=void_fraction,
                axis=axis,
                corner_radius=radius_tf_turn_cable_space_corners,
                n_strands=mfile.get("n_tf_turn_superconducting_cables", scan=scan),
                f_a_tf_turn_cable_copper=mfile.get(
                    "f_a_tf_turn_cable_copper", scan=scan
                ),
            )

        axis.set_xlim(-turn_width * 0.05, turn_width * 1.05)
        axis.set_ylim(-turn_height * 0.05, turn_height * 1.05)

    axis.minorticks_on()
    axis.set_title("WP Turn Structure")
    axis.set_xlabel("r [m]")
    axis.set_ylabel("x [m]")

    # Add info about the steel casing surrounding the WP
    textstr_turn_insulation = (
        f"$\\mathbf{{Turn \\ Insulation:}}$\n\n$\\Delta r:${insulation_thickness:.3e} m"
    )

    draw_text(
        axis,
        0.4,
        0.9,
        textstr_turn_insulation,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("red"),
    )

    # Add info about the steel casing surrounding the WP
    textstr_turn_steel = (
        f"$\\mathbf{{Steel \\ Conduit:}}$\n\n$\\Delta r:${steel_thickness:.3e}"
        f" m\n$A$: {a_tf_turn_steel:.3e} m$^2$"
    )

    draw_text(
        axis,
        0.65,
        0.9,
        textstr_turn_steel,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("grey"),
    )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        # Add info about the steel casing surrounding the WP
        textstr_turn_cable_space = (
            "$\\mathbf{Cable \\ Space:}$\n\n$\\Delta r:$"
            f" {cable_space_width:.3e} m\nCorner radius, $r$:"
            f" {radius_tf_turn_cable_space_corners:.3e} m\nCable area with no"
            f" cooling\nchannel or gaps: {a_tf_turn_cable_space_no_void:.3e}"
            " m$^2$\nExtra cable space area void fraction:"
            f" {f_a_tf_turn_cable_space_extra_void}\nTrue cable space area:"
            f" {a_tf_turn_cable_space_effective:.3e} m$^2$"
        )
    elif TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        textstr_turn_cable_space = (
            "$\\mathbf{Cable \\ Space:}$\n\nCable space:\n$\\Delta r$:"
            f" {cable_space_width_radial:.3e} m\n$\\Delta x$:"
            f" {cable_space_width_toroidal:.3e} m\nCorner radius, $r$:"
            f" {radius_tf_turn_cable_space_corners:.3e} m\nCable area with no"
            f" cooling channel or gaps: {a_tf_turn_cable_space_no_void:.3e}"
            " m$^2$\nExtra cable space area void fraction:"
            f" {f_a_tf_turn_cable_space_extra_void}\nTrue cable space area:"
            f" {a_tf_turn_cable_space_effective:.3e} m$^2$"
        )

    draw_text(
        axis,
        0.40,
        0.7,
        textstr_turn_cable_space,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("royalblue"),
    )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        textstr_turn = (
            "$\\mathbf{Turn:}$\n\n"
            f"$\\Delta r$: {turn_width:.3e} m\n"
            f"$\\Delta x$: {turn_width:.3e} m"
        )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        textstr_turn = (
            "$\\mathbf{Turn:}$\n\n"
            f"$\\Delta r$: {turn_width:.3e} m\n"
            f"$\\Delta x$: {turn_height:.3e} m"
        )

    draw_text(
        axis,
        0.525,
        0.9,
        textstr_turn,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("wheat"),
    )

    # Add info about the steel casing surrounding the WP
    textstr_turn_cooling = (
        f"$\\mathbf{{Cooling:}}$\n\n$\\varnothing$: {he_pipe_diameter:.3e}"
        " m\nTotal area of all coolant channels:"
        f" {a_tf_wp_coolant_channels:.4f} m$^2$"
    )

    draw_text(
        axis,
        0.45,
        0.8,
        textstr_turn_cooling,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("white"),
    )

    textstr_superconductor = (
        "$\\mathbf{Superconductor:}$\n\nSuperconductor"
        f" used:\n{SuperconductorModel(mfile.get('i_tf_sc_mat', scan=scan)).full_name}\nCritical"  # noqa: E501
        " field at zero\ntemperature and strain:"
        f" {mfile.get('b_tf_superconductor_critical_zero_temp_strain', scan=scan):.4f}"
        " T\nCritical temperature at\nzero field and strain:"
        f" {mfile.get('temp_tf_superconductor_critical_zero_field_strain', scan=scan):.4f}"  # noqa: E501
        f" K\nTemperature at conductor: {mfile.get('tftmp', scan=scan):.4f}"
        " K\nField at conductor:"
        f" {mfile.get('b_tf_inboard_peak_with_ripple', scan=scan):.4f}"
        " T\nSuperconductor critical current density at\noperating"
        " conditions:"
        f" {mfile.get('j_tf_superconductor_critical', scan=scan):.2e}"
        " A/m$^2$\n$I_{\\text{TF,turn critical}}$:"
        f" {mfile.get('c_turn_cables_critical', scan=scan):,.2f}"
        " A\n$I_{\\text{TF,turn}}$:"
        f" {mfile.get('c_tf_turn', scan=scan):,.2f} A\nCritcal current ratio:"
        f" {mfile.get('f_c_tf_turn_operating_critical', scan=scan):,.4f}\nSuperconductor"
        " temperature\nmargin:"
        f" {mfile.get('temp_tf_superconductor_margin', scan=scan):,.4f}"
        " K\n\n$\\mathbf{Quench:}$\n\nQuench dump time:"
        f" {mfile.get('t_tf_superconductor_quench', scan=scan):.4f} s\nQuench"
        f" detection time: {mfile.get('t_tf_quench_detection', scan=scan):.4f}"
        " s\nUser input max temperature\nduring quench:"
        f" {mfile.get('temp_tf_conductor_quench_max', scan=scan):.2f}"
        " K\nRequired maxium WP current\ndensity for heat"
        f" protection:\n{mfile.get('j_tf_wp_quench_heat_max', scan=scan):.2e}"
        " A/m$^2$\n"
    )
    draw_text(
        axis,
        0.75,
        0.9,
        textstr_superconductor,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("#6dd3f7"),  # light blue for superconductors
    )


def plot_tf_croco_turn(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plots inboard TF coil CICC individual turn structure with croco cable layout."""
    # Import the TF turn variables then multiply into mm
    i_tf_turns_integer = mfile.get("i_tf_turns_integer", scan=scan)
    # If integer turns switch is on then the turns can have non square dimensions
    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        turn_width = mfile.get("dr_tf_turn", scan=scan)
        turn_height = mfile.get("dx_tf_turn", scan=scan)
        cable_space_width_radial = mfile.get("dr_tf_turn_cable_space", scan=scan)
        cable_space_width_toroidal = mfile.get("dx_tf_turn_cable_space", scan=scan)

    elif TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        turn_width = mfile.get("dx_tf_turn_general", scan=scan)
        cable_space_width = mfile.get("dx_tf_turn_cable_space_average", scan=scan)

    steel_thickness = mfile.get("dx_tf_turn_steel", scan=scan)
    insulation_thickness = mfile.get("dx_tf_turn_insulation", scan=scan)

    a_tf_turn_cable_space_no_void = mfile.get("a_tf_turn_cable_space_no_void", scan=scan)
    radius_tf_turn_cable_space_corners = mfile.get(
        "radius_tf_turn_cable_space_corners", scan=scan
    )

    a_tf_wp_coolant_channels = mfile.get("a_tf_wp_coolant_channels", scan=scan)

    f_a_tf_turn_cable_space_extra_void = mfile.get(
        "f_a_tf_turn_cable_space_extra_void", scan=scan
    )
    a_tf_turn_steel = mfile.get("a_tf_turn_steel", scan=scan)
    a_tf_turn_cable_space_effective = mfile.get(
        "a_tf_turn_cable_space_effective", scan=scan
    )

    he_pipe_diameter = mfile.get("dia_tf_turn_coolant_channel", scan=scan)
    dia_tf_turn_croco_cable = mfile.get("dia_tf_turn_croco_cable", scan=scan)

    # Plot the total turn shape
    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        axis.add_patch(
            Rectangle(
                (0, 0),
                turn_width,
                turn_width,
                facecolor="red",
                edgecolor="black",
            ),
        )
        # Plot the steel conduit
        axis.add_patch(
            Rectangle(
                (insulation_thickness, insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                facecolor="grey",
                edgecolor="black",
            ),
        )

        # Plot the central cable space and copper cylinder
        for rad, col in [
            (1.5 * dia_tf_turn_croco_cable, "white"),
            (dia_tf_turn_croco_cable / 2, "#B87333"),
        ]:
            axis.add_patch(
                Circle(
                    ((turn_width / 2), (turn_width / 2)),
                    rad,
                    facecolor=col,
                    edgecolor="black",
                    linewidth=1.2,
                ),
            )

        # Plot six surrounding Croco cables in a hexagonal layout.
        center_x = turn_width / 2
        center_y = turn_width / 2
        ring_radius = dia_tf_turn_croco_cable
        for angle in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            plot_corc_cable_geometry(
                axis=axis,
                r_centre=center_x + ring_radius * np.cos(angle),
                z_centre=center_y + ring_radius * np.sin(angle),
                dia_croco_strand=mfile.get("dia_tf_turn_croco_cable", scan=scan),
                dx_croco_strand_copper=mfile.get("dx_tf_croco_strand_copper", scan=scan),
                dr_hts_tape=mfile.get("dr_tf_hts_tape", scan=scan),
                dx_croco_strand_tape_stack=mfile.get(
                    "dx_tf_croco_strand_tape_stack", scan=scan
                ),
                n_croco_strand_hts_tapes=mfile.get(
                    "n_tf_croco_strand_hts_tapes", scan=scan
                ),
                dx_hts_tape_rebco=mfile.get("dx_tf_hts_tape_rebco", scan=scan),
                dx_hts_tape_copper=mfile.get("dx_tf_hts_tape_copper", scan=scan),
                dx_hts_tape_hastelloy=mfile.get("dx_tf_hts_tape_hastelloy", scan=scan),
                show_legend=False,
            )

    axis.minorticks_on()
    axis.set_title("WP Turn Structure")
    axis.set_xlim(-turn_width * 0.025, turn_width * 1.025)
    axis.set_ylim(-turn_width * 0.025, turn_width * 1.025)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("r [m]")
    axis.set_ylabel("x [m]")

    # Add info about the steel casing surrounding the WP
    textstr_turn_insulation = (
        f"$\\mathbf{{Turn \\ Insulation:}}$\n\n$\\Delta r:${insulation_thickness:.3e} m"
    )

    draw_text(
        axis,
        0.4,
        0.9,
        textstr_turn_insulation,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("red"),
    )

    # Add info about the steel casing surrounding the WP
    textstr_turn_steel = (
        f"$\\mathbf{{Steel \\ Conduit:}}$\n\n$\\Delta r:${steel_thickness:.3e}"
        f" m\n$A$: {a_tf_turn_steel:.3e} m$^2$"
    )

    draw_text(
        axis,
        0.65,
        0.9,
        textstr_turn_steel,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("grey"),
    )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        # Add info about the steel casing surrounding the WP
        textstr_turn_cable_space = (
            "$\\mathbf{Cable \\ Space:}$\n\n$\\Delta r:$"
            f" {cable_space_width:.3e} m\nCorner radius, $r$:"
            f" {radius_tf_turn_cable_space_corners:.3e} m\nCable area with no"
            f" cooling\nchannel or gaps: {a_tf_turn_cable_space_no_void:.3e}"
            " m$^2$\nExtra cable space area void fraction:"
            f" {f_a_tf_turn_cable_space_extra_void}\nTrue cable space area:"
            f" {a_tf_turn_cable_space_effective:.3e} m$^2$"
        )
    elif TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        textstr_turn_cable_space = (
            "$\\mathbf{Cable \\ Space:}$\n\nCable space:\n$\\Delta r$:"
            f" {cable_space_width_radial:.3e} m\n$\\Delta x$:"
            f" {cable_space_width_toroidal:.3e} m\nCorner radius, $r$:"
            f" {radius_tf_turn_cable_space_corners:.3e} m\nCable area with no"
            f" cooling channel or gaps: {a_tf_turn_cable_space_no_void:.3e}"
            " m$^2$\nExtra cable space area void fraction:"
            f" {f_a_tf_turn_cable_space_extra_void}\nTrue cable space area:"
            f" {a_tf_turn_cable_space_effective:.3e} m$^2$"
        )

    draw_text(
        axis,
        0.40,
        0.7,
        textstr_turn_cable_space,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("royalblue"),
    )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.NON_INTEGER:
        textstr_turn = (
            "$\\mathbf{Turn:}$\n\n"
            f"$\\Delta r$: {turn_width:.3e} m\n"
            f"$\\Delta x$: {turn_width:.3e} m"
        )

    if TFWPIntegerTurnType(i_tf_turns_integer) == TFWPIntegerTurnType.INTEGER:
        textstr_turn = (
            "$\\mathbf{Turn:}$\n\n"
            f"$\\Delta r$: {turn_width:.3e} m\n"
            f"$\\Delta x$: {turn_height:.3e} m"
        )

    draw_text(
        axis,
        0.525,
        0.9,
        textstr_turn,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("wheat"),
    )

    # Add info about the steel casing surrounding the WP
    textstr_turn_cooling = (
        f"$\\mathbf{{Cooling:}}$\n\n$\\varnothing$: {he_pipe_diameter:.3e}"
        " m\nTotal area of all coolant channels:"
        f" {a_tf_wp_coolant_channels:.4f} m$^2$"
    )

    draw_text(
        axis,
        0.45,
        0.8,
        textstr_turn_cooling,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("white"),
    )

    textstr_superconductor = (
        "$\\mathbf{Superconductor:}$\n\nSuperconductor used:"
        f" {SuperconductorModel(mfile.get('i_tf_sc_mat', scan=scan)).full_name}\nCritical"  # noqa: E501
        " field at zero\ntemperature and strain:"
        f" {mfile.get('b_tf_superconductor_critical_zero_temp_strain', scan=scan):.4f}"
        " T\nCritical temperature at\nzero field and strain:"
        f" {mfile.get('temp_tf_superconductor_critical_zero_field_strain', scan=scan):.4f}"  # noqa: E501
        f" K\nTemperature at conductor: {mfile.get('tftmp', scan=scan):.4f}"
        " K\nField at conductor:"
        f" {mfile.get('b_tf_inboard_peak_with_ripple', scan=scan):.4f}"
        " T\nSuperconductor critical current density at\noperating"
        " conditions:"
        f" {mfile.get('j_tf_superconductor_critical', scan=scan):.2e}"
        " A/m$^2$\n$I_{\\text{TF,turn critical}}$:"
        f" {mfile.get('c_turn_cables_critical', scan=scan):,.2f}"
        " A\n$I_{\\text{TF,turn}}$:"
        f" {mfile.get('c_tf_turn', scan=scan):,.2f} A\nCritcal current ratio:"
        f" {mfile.get('f_c_tf_turn_operating_critical', scan=scan):,.4f}\nSuperconductor"
        " temperature\nmargin:"
        f" {mfile.get('temp_tf_superconductor_margin', scan=scan):,.4f}"
        " K\n\n$\\mathbf{Quench:}$\n\nQuench dump time:"
        f" {mfile.get('t_tf_superconductor_quench', scan=scan):.4e} s\nQuench"
        f" detection time: {mfile.get('t_tf_quench_detection', scan=scan):.4e}"
        " s\nUser input max temperature\nduring quench:"
        f" {mfile.get('temp_tf_conductor_quench_max', scan=scan):.2f}"
        " K\nRequired maxium WP current\ndensity for heat"
        f" protection:\n{mfile.get('j_tf_wp_quench_heat_max', scan=scan):.2e}"
        " A/m$^2$\n"
    )
    draw_text(
        axis,
        0.75,
        0.9,
        textstr_superconductor,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("#6dd3f7"),
    )


def plot_tf_coil_structure(axis: plt.Axes, mfile: MFile, scan: int, colour_scheme=1):
    """Plot the TF coil poloidal cross-section"""
    plot_tf_coils(axis, mfile, scan, colour_scheme)

    x1 = mfile.get("r_tf_arc(1)", scan=scan)
    y1 = mfile.get("z_tf_arc(1)", scan=scan)
    x2 = mfile.get("r_tf_arc(2)", scan=scan)
    y2 = mfile.get("z_tf_arc(2)", scan=scan)
    x3 = mfile.get("r_tf_arc(3)", scan=scan)
    y3 = mfile.get("z_tf_arc(3)", scan=scan)
    x4 = mfile.get("r_tf_arc(4)", scan=scan)
    y4 = mfile.get("z_tf_arc(4)", scan=scan)
    x5 = mfile.get("r_tf_arc(5)", scan=scan)
    y5 = mfile.get("z_tf_arc(5)", scan=scan)

    z_tf_inside_half = mfile.get("z_tf_inside_half", scan=scan)
    z_tf_top = mfile.get("z_tf_top", scan=scan)
    dr_tf_inboard = mfile.get("dr_tf_inboard", scan=scan)
    r_tf_inboard_out = mfile.get("r_tf_inboard_out", scan=scan)
    r_tf_outboard_in = mfile.get("r_tf_outboard_in", scan=scan)
    r_tf_inboard_in = mfile.get("r_tf_inboard_in", scan=scan)
    dr_tf_outboard = mfile.get("dr_tf_outboard", scan=scan)
    len_tf_coil = mfile.get("len_tf_coil", scan=scan)
    dz_tf_upper_lower_midplane = mfile.get("dz_tf_upper_lower_midplane", scan=scan)

    # Plot the points as black dots, number them, and connect them with lines
    xs = [x1, x2, x3, x4, x5]
    ys = [y1, y2, y3, y4, y5]
    labels = []
    for i, (x, y) in enumerate(zip(xs, ys, strict=False), 1):
        axis.plot(x, y, "ko", markersize=8)
        draw_text(
            axis,
            x,
            y,
            str(i),
            color="red",
            fontsize=5,
            ha="center",
            va="center",
            fontweight="bold",
        )
        labels.append(f"TF Arc Point {i}: ({x:.2f}, {y:.2f})")

    # =========================================================

    # If D-shaped coil, plot the full internal height arrow
    if mfile.get("i_tf_shape", scan=scan) == 1:
        # Arrow for internal coil width
        draw_annotation(
            axis,
            "",
            xy=(x2, y2),
            xytext=(x4, y4),
            arrowprops={"arrowstyle": "<->", "color": "black"},
        )

        # Add a label for the internal coil width
        draw_text(
            axis,
            x2,
            0.0,
            f"{y2 - y4:.3f} m",
            fontsize=7,
            color="black",
            rotation=270,
            verticalalignment="center",
            horizontalalignment="center",
            bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
            zorder=100,  # Ensure label is on top of all plots
        )

    # ==========================================================

    # Arrow for the full TF coil height
    if mfile.get("i_tf_shape", scan=scan) == 1:
        x = x2 * 0.9
    elif mfile.get("i_tf_shape", scan=scan) == 2:
        x = (x2 - x1) / 2

    draw_annotation(
        axis,
        "",
        xy=(x, y2 + dr_tf_inboard),
        xytext=(x, y4 - dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the full TF coil height
    draw_text(
        axis,
        x,
        0.0,
        f"{((y2 + 2 * dr_tf_inboard) - y4):.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=101,  # Ensure label is on top of all plots
    )

    # ==========================================================

    # Arrow for top half height of TF coil
    draw_annotation(
        axis,
        "",
        xy=(-2.0, 0),
        xytext=(-2.0, y2 + dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y2 + dr_tf_inboard, color="black", linestyle="--", linewidth=1)

    # Add a label for top of TF coil
    draw_text(
        axis,
        -2.0,
        (y2 + dr_tf_inboard) / 2,
        f"{y2 + dr_tf_inboard:.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # ==========================================================

    # Arrow for bottom half height of TF coil
    draw_annotation(
        axis,
        "",
        xy=(-2.0, 0),
        xytext=(-2.0, y4 - dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y4 - dr_tf_inboard, color="black", linestyle="--", linewidth=1)

    # Add a label for top of TF coil
    draw_text(
        axis,
        -2.0,
        -z_tf_top / 2,
        f"{y4 - dr_tf_inboard:.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # Arrow for top inside internal height
    draw_annotation(
        axis,
        "",
        xy=(-1.0, 0),
        xytext=(-1.0, y2),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y2, color="black", linestyle="--", linewidth=1)

    # Add a label for height of top internal height
    draw_text(
        axis,
        -1.0,
        y2 / 2,
        f"{y2:.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # =========================================================

    # Arrow for coil internal height
    draw_annotation(
        axis,
        "",
        xy=(-1.0, 0),  # Inner plasma edge
        xytext=(-1.0, -z_tf_inside_half),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=-z_tf_inside_half, color="black", linestyle="--", linewidth=1)

    # Add a label for coil internal height
    draw_text(
        axis,
        -1.0,
        -z_tf_inside_half / 2,
        f"{z_tf_inside_half:.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )
    # =========================================================

    # Arrow for internal coil width
    draw_annotation(
        axis,
        "",
        xy=(r_tf_inboard_out, -z_tf_inside_half / 12),
        xytext=(r_tf_outboard_in, -z_tf_inside_half / 12),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the internal coil width
    draw_text(
        axis,
        (r_tf_inboard_out + r_tf_outboard_in) / 1.5,
        -z_tf_inside_half / 12,
        f"{mfile.get('dr_tf_internal_midplane', scan=scan):.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=100,  # Ensure label is on top of all plots
    )

    # =============================================================

    # Arrow for full coil width
    draw_annotation(
        axis,
        "",
        xy=(r_tf_inboard_in, 0.0),
        xytext=(r_tf_outboard_in + dr_tf_outboard, 0.0),
        arrowprops={"arrowstyle": "<|-|>", "color": "black"},
        zorder=100,  # Ensure label is on top of all plots
    )

    # Add a label for the full coil width
    draw_text(
        axis,
        (r_tf_inboard_out + r_tf_outboard_in) / 1.5,
        0.0,
        f"{mfile.get('dr_tf_full_midplane', scan=scan):.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=100,  # Ensure label is on top of all plots
    )

    # =============================================================

    # Plot vertical lines for the inboard TF coil start and end
    axis.axvline(
        r_tf_inboard_in,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
        label="TF Inboard Start",
    )
    axis.axvline(
        r_tf_inboard_out,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
        label="TF Inboard End",
    )
    # Plot vertical lines for the outboard TF coil start and end
    axis.axvline(
        r_tf_outboard_in,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
        label="TF Outboard Start",
    )
    axis.axvline(
        r_tf_outboard_in + dr_tf_outboard,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
        label="TF Outboard End",
    )

    # Add a label for the inboard thickness
    draw_text(
        axis,
        r_tf_inboard_in,
        (y4 - dr_tf_inboard) * 1.1,
        rf"$\Delta r = ${dr_tf_inboard:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # Add a label for the outboard thickness
    draw_text(
        axis,
        r_tf_outboard_in,
        (y4 - dr_tf_inboard) * 1.1,
        rf"$\Delta r = ${dr_tf_outboard:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # ==============================================================

    # Add a label for the length of the coil
    draw_text(
        axis,
        (r_tf_outboard_in + 2 * dr_tf_outboard),
        0.0,
        rf"Length of coil = {len_tf_coil:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=100,  # Ensure label is on top of all plots
    )

    # ==============================================================

    # Add a label for the length of the coil
    draw_text(
        axis,
        (r_tf_outboard_in + 2 * dr_tf_outboard),
        -1.0,
        f"$\\Delta Z$ upper and lower to midplane = {dz_tf_upper_lower_midplane:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # ==============================================================

    # Add arow for inboard coil radius
    draw_annotation(
        axis,
        "",
        xy=(r_tf_inboard_in, 0),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "->", "color": "black"},
    )

    # Add label for inboard coil radius
    draw_text(
        axis,
        r_tf_inboard_in / 2,
        0.0,
        f"{r_tf_inboard_in:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=101,  # Ensure label is on top of all plots
    )

    # =============================================================

    # ==============================================================

    if mfile.get("i_tf_shape", scan=scan) == 1:
        # Add arow for inboard coil radius
        draw_annotation(
            axis,
            "",
            xy=(r_tf_outboard_in + dr_tf_outboard, y2 + dr_tf_inboard),
            xytext=(0, y2 + dr_tf_inboard),
            arrowprops={"arrowstyle": "->", "color": "black"},
        )

        # Add label for inboard coil radius
        draw_text(
            axis,
            r_tf_inboard_in / 2,
            y2 + dr_tf_inboard,
            f"{r_tf_outboard_in + dr_tf_outboard:.3f} m",
            fontsize=7,
            color="black",
            verticalalignment="center",
            horizontalalignment="center",
            bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        )

        axis.plot(
            0,
            y2 + dr_tf_inboard,
            marker="o",
            color="black",
            markersize=7,
            zorder=100,
        )

    # ==============================================================

    y_center = y2 - ((y2 - y4) / 2)
    # also draw a red horizontal line at the same vertical centre
    axis.axhline(y=y_center, color="red", linestyle="--", linewidth=1.0, zorder=5)

    # Add a label the plasma and TF vertical centre distance offset
    draw_text(
        axis,
        (r_tf_outboard_in + 2 * dr_tf_outboard),
        -2.0,
        "$\\Delta Z$ coil centre to plasma centre ="
        f" {mfile.get('dz_tf_plasma_centre_offset', scan=scan):.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # =============================================================

    # Plot a red dot at (0,0)
    axis.plot(0, 0, marker="o", color="red", markersize=7)

    # Plot a red dashed vertical line at R=0
    axis.axvline(0, color="red", linestyle="--", linewidth=1)

    # Add centre line at
    axis.axhline(y=0, color="red", linestyle="--", linewidth=1)
    axis.set_xlim(-3.0, (r_tf_outboard_in + dr_tf_outboard) * 1.4)
    axis.set_ylim((y4 - dr_tf_inboard) * 1.2, (y2 + dr_tf_inboard) * 1.2)
    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("TF Coil Poloidal Cross-Section")
    axis.minorticks_on()
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    # Move the legend to above the plot
    axis.legend(labels, loc="upper center", bbox_to_anchor=(1.01, 0.85), ncol=1)


def plot_tf_stress(axis: plt.Axes, mfile: MFile):
    """Function to plot the TF coil stress from the SIG_TF.json file.

    Input file:
    SIG_TF.json

    Parameters
    ----------
    axis: plt.Axes :

    mfile: MFile :

    """
    # Step 1 : Data extraction
    # ----------------------------------------------------------------------------------------------  # noqa: E501
    # Number of physical quantity value per coil layer
    n_radial_array_layer = 0

    # Physical quantities : full vectors
    radius = []
    radial_smeared_stress = []
    toroidal_smeared_stress = []
    vertical_smeared_stress = []
    tresca_smeared_stress = []
    radial_stress = []
    toroidal_stress = []
    vertical_stress = []
    vm_stress = []
    tresca_stress = []
    cea_tresca_stress = []
    radial_strain = []
    toroidal_strain = []
    vertical_strain = []
    radial_displacement = []

    # Physical quantity : WP stress
    wp_vertical_stress = []

    # Physical quantity : values at layer border
    bound_radius = []
    bound_radial_smeared_stress = []
    bound_toroidal_smeared_stress = []
    bound_vertical_smeared_stress = []
    bound_tresca_smeared_stress = []
    bound_radial_stress = []
    bound_toroidal_stress = []
    bound_vertical_stress = []
    bound_vm_stress = []
    bound_tresca_stress = []
    bound_cea_tresca_stress = []
    bound_radial_strain = []
    bound_toroidal_strain = []
    bound_vertical_strain = []
    bound_radial_displacement = []

    with open(
        mfile.filename.with_name(mfile.filename.name.replace("MFILE.DAT", "SIG_TF.json"))
    ) as f:
        sig_data = json.load(f)

    # Getting the data to be plotted
    n_radial_array_layer = sig_data["Points per layers"]
    n_points = len(sig_data["Radius (m)"])
    n_layers = int(n_points / n_radial_array_layer)
    for ii in range(n_layers):
        # Full vector
        radius.append([])
        radial_stress.append([])
        toroidal_stress.append([])
        vertical_stress.append([])
        radial_smeared_stress.append([])
        toroidal_smeared_stress.append([])
        vertical_smeared_stress.append([])
        vm_stress.append([])
        tresca_stress.append([])
        cea_tresca_stress.append([])
        radial_displacement.append([])

        for jj in range(n_radial_array_layer):
            radius[ii].append(sig_data["Radius (m)"][ii * n_radial_array_layer + jj])
            radial_stress[ii].append(
                sig_data["Radial stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            toroidal_stress[ii].append(
                sig_data["Toroidal stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            if len(sig_data["Vertical stress (MPa)"]) == 1:
                vertical_stress[ii].append(sig_data["Vertical stress (MPa)"][0])
            else:
                vertical_stress[ii].append(
                    sig_data["Vertical stress (MPa)"][ii * n_radial_array_layer + jj]
                )
            radial_smeared_stress[ii].append(
                sig_data["Radial smear stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            toroidal_smeared_stress[ii].append(
                sig_data["Toroidal smear stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            vertical_smeared_stress[ii].append(
                sig_data["Vertical smear stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            vm_stress[ii].append(
                sig_data["Von-Mises stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            tresca_stress[ii].append(
                sig_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            cea_tresca_stress[ii].append(
                sig_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            radial_displacement[ii].append(
                sig_data["rad. displacement (mm)"][ii * n_radial_array_layer + jj]
            )

        # Layer lower boundaries values
        bound_radius.append(sig_data["Radius (m)"][ii * n_radial_array_layer])
        bound_radial_stress.append(
            sig_data["Radial stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_toroidal_stress.append(
            sig_data["Toroidal stress (MPa)"][ii * n_radial_array_layer]
        )
        if len(sig_data["Vertical stress (MPa)"]) == 1:
            bound_vertical_stress.append(sig_data["Vertical stress (MPa)"][0])
        else:
            bound_vertical_stress.append(
                sig_data["Vertical stress (MPa)"][ii * n_radial_array_layer]
            )
        bound_radial_smeared_stress.append(
            sig_data["Radial smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_toroidal_smeared_stress.append(
            sig_data["Toroidal smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_vertical_smeared_stress.append(
            sig_data["Vertical smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_vm_stress.append(
            sig_data["Von-Mises stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_tresca_stress.append(
            sig_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_cea_tresca_stress.append(
            sig_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_radial_displacement.append(
            sig_data["rad. displacement (mm)"][ii * n_radial_array_layer]
        )

        # Layer upper boundaries values
        bound_radius.append(sig_data["Radius (m)"][(ii + 1) * n_radial_array_layer - 1])
        bound_radial_stress.append(
            sig_data["Radial stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_toroidal_stress.append(
            sig_data["Toroidal stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        if len(sig_data["Vertical stress (MPa)"]) == 1:
            bound_vertical_stress.append(sig_data["Vertical stress (MPa)"][0])
        else:
            bound_vertical_stress.append(
                sig_data["Vertical stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
            )
        bound_radial_smeared_stress.append(
            sig_data["Radial smear stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_toroidal_smeared_stress.append(
            sig_data["Toroidal smear stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_vertical_smeared_stress.append(
            sig_data["Vertical smear stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_vm_stress.append(
            sig_data["Von-Mises stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_tresca_stress.append(
            sig_data["CEA Tresca stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_cea_tresca_stress.append(
            sig_data["CEA Tresca stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_radial_displacement.append(
            sig_data["rad. displacement (mm)"][(ii + 1) * n_radial_array_layer - 1]
        )

    # TRESCA smeared stress [MPa]
    for ii in range(n_layers):
        tresca_smeared_stress.append([])

        bound_tresca_smeared_stress.extend([
            max(
                abs(radial_smeared_stress[ii][0]),
                abs(toroidal_smeared_stress[ii][0]),
            )
            + vertical_smeared_stress[ii][0],
            max(
                abs(radial_smeared_stress[ii][n_radial_array_layer - 1]),
                abs(toroidal_smeared_stress[ii][n_radial_array_layer - 1]),
            )
            + vertical_smeared_stress[ii][n_radial_array_layer - 1],
        ])
        for jj in range(n_radial_array_layer):
            tresca_smeared_stress[ii].append(
                max(
                    abs(radial_smeared_stress[ii][jj]),
                    abs(toroidal_smeared_stress[ii][jj]),
                )
                + vertical_smeared_stress[ii][jj]
            )

    # Strains
    if len(sig_data) > 16:
        for ii in range(n_layers):
            radial_strain.append([])
            toroidal_strain.append([])
            vertical_strain.append([])

            bound_radial_strain.extend([
                sig_data["Radial strain"][ii * n_radial_array_layer],
                sig_data["Radial strain"][(ii + 1) * n_radial_array_layer - 1],
            ])
            bound_toroidal_strain.extend([
                sig_data["Toroidal strain"][ii * n_radial_array_layer],
                sig_data["Toroidal strain"][(ii + 1) * n_radial_array_layer - 1],
            ])
            bound_vertical_strain.extend([
                sig_data["Vertical strain"][ii * n_radial_array_layer],
                sig_data["Vertical strain"][(ii + 1) * n_radial_array_layer - 1],
            ])
            for jj in range(n_radial_array_layer):
                radial_strain[ii].append(
                    sig_data["Radial strain"][ii * n_radial_array_layer + jj]
                )
                toroidal_strain[ii].append(
                    sig_data["Toroidal strain"][ii * n_radial_array_layer + jj]
                )
                vertical_strain[ii].append(
                    sig_data["Vertical strain"][ii * n_radial_array_layer + jj]
                )

                if "WP smeared stress (MPa)" in sig_data:
                    wp_vertical_stress.append(sig_data["WP smeared stress (MPa)"][jj])

    axis_tick_size = 12
    legend_size = 10
    mark_size = 10
    line_width = 3.5

    # PLOT 1 : Stress summary
    # ------------------------

    ax = axis[0]
    for ii in range(n_layers):
        ax.plot(
            radius[ii],
            radial_stress[ii],
            "-",
            linewidth=line_width,
            color="lightblue",
        )
        ax.plot(
            radius[ii],
            toroidal_stress[ii],
            "-",
            linewidth=line_width,
            color="wheat",
        )
        ax.plot(
            radius[ii],
            vertical_stress[ii],
            "-",
            linewidth=line_width,
            color="lightgrey",
        )
        ax.plot(
            radius[ii],
            tresca_stress[ii],
            "-",
            linewidth=line_width,
            color="pink",
        )
        ax.plot(
            radius[ii],
            vm_stress[ii],
            "-",
            linewidth=line_width,
            color="violet",
        )
    ax.plot(
        radius[0],
        radial_stress[0],
        "--",
        color="dodgerblue",
        label=r"$\sigma_{rr}$",
    )
    ax.plot(
        radius[0],
        toroidal_stress[0],
        "--",
        color="orange",
        label=r"$\sigma_{\theta\theta}$",
    )
    ax.plot(
        radius[0],
        vertical_stress[0],
        "--",
        color="mediumseagreen",
        label=r"$\sigma_{zz}$",
    )
    ax.plot(
        radius[0],
        tresca_stress[0],
        "-",
        color="crimson",
        label=r"$\sigma_{TRESCA}$",
    )
    ax.plot(
        radius[0],
        vm_stress[0],
        "-",
        color="darkviolet",
        label=r"$\sigma_{Von\ mises}$",
    )
    for ii in range(1, n_layers):
        ax.plot(radius[ii], radial_stress[ii], "--", color="dodgerblue")
        ax.plot(radius[ii], toroidal_stress[ii], "--", color="orange")
        ax.plot(radius[ii], vertical_stress[ii], "--", color="mediumseagreen")
        ax.plot(radius[ii], tresca_stress[ii], "-", color="crimson")
        ax.plot(radius[ii], vm_stress[ii], "-", color="darkviolet")
    ax.plot(
        bound_radius,
        bound_radial_stress,
        "|",
        markersize=mark_size,
        color="dodgerblue",
    )
    ax.plot(
        bound_radius,
        bound_toroidal_stress,
        "|",
        markersize=mark_size,
        color="orange",
    )
    ax.plot(
        bound_radius,
        bound_vertical_stress,
        "|",
        markersize=mark_size,
        color="mediumseagreen",
    )
    ax.plot(
        bound_radius,
        bound_tresca_stress,
        "|",
        markersize=mark_size,
        color="crimson",
    )
    ax.plot(
        bound_radius,
        bound_vm_stress,
        "|",
        markersize=mark_size,
        color="darkviolet",
    )
    ax.grid(True)
    ax.set_ylabel(r"$\sigma$ [$MPa$]", fontsize=axis_tick_size)
    ax.set_title("Structure Stress Summary")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=legend_size)

    # PLOT 2 : Smeared stress summary
    # ------------------------
    ax = axis[1]
    for ii in range(n_layers):
        ax.plot(
            radius[ii],
            radial_smeared_stress[ii],
            "-",
            linewidth=line_width,
            color="lightblue",
        )
        ax.plot(
            radius[ii],
            toroidal_smeared_stress[ii],
            "-",
            linewidth=line_width,
            color="wheat",
        )
        ax.plot(
            radius[ii],
            vertical_smeared_stress[ii],
            "-",
            linewidth=line_width,
            color="lightgrey",
        )
        ax.plot(
            radius[ii],
            tresca_smeared_stress[ii],
            "-",
            linewidth=line_width,
            color="pink",
        )
    ax.plot(
        radius[0],
        radial_smeared_stress[0],
        "--",
        color="dodgerblue",
        label=r"$\sigma_{rr}^\mathrm{smeared}$",
    )
    ax.plot(
        radius[0],
        toroidal_smeared_stress[0],
        "--",
        color="orange",
        label=r"$\sigma_{\theta\theta}^\mathrm{smeared}$",
    )
    ax.plot(
        radius[0],
        vertical_smeared_stress[0],
        "--",
        color="mediumseagreen",
        label=r"$\sigma_{zz}^\mathrm{smeared}$",
    )
    ax.plot(
        radius[0],
        tresca_smeared_stress[0],
        "-",
        color="crimson",
        label=r"$\sigma_{TRESCA}^\mathrm{smeared}$",
    )
    for ii in range(1, n_layers):
        ax.plot(radius[ii], radial_smeared_stress[ii], "--", color="dodgerblue")
        ax.plot(radius[ii], toroidal_smeared_stress[ii], "--", color="orange")
        ax.plot(
            radius[ii],
            vertical_smeared_stress[ii],
            "--",
            color="mediumseagreen",
        )
        ax.plot(radius[ii], tresca_smeared_stress[ii], "-", color="crimson")
    ax.plot(
        bound_radius,
        bound_radial_smeared_stress,
        "|",
        markersize=mark_size,
        color="dodgerblue",
    )
    ax.plot(
        bound_radius,
        bound_toroidal_smeared_stress,
        "|",
        markersize=mark_size,
        color="orange",
    )
    ax.plot(
        bound_radius,
        bound_vertical_smeared_stress,
        "|",
        markersize=mark_size,
        color="mediumseagreen",
    )
    ax.plot(
        bound_radius,
        bound_tresca_smeared_stress,
        "|",
        markersize=mark_size,
        color="crimson",
    )
    ax.grid(True)
    ax.set_ylabel(r"$\sigma$ [$MPa$]", fontsize=axis_tick_size)
    ax.set_title("Smeared Stress Summary")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=legend_size)

    # PLOT 4 : Displacement
    # ----------------------
    ax = axis[2]
    ax.plot(radius[0], radial_displacement[0], color="dodgerblue")
    for ii in range(1, n_layers):
        ax.plot(radius[ii], radial_displacement[ii], color="dodgerblue")
    ax.grid(True)
    ax.set_ylabel(r"$u_{r}$ [mm]", fontsize=axis_tick_size)
    ax.set_xlabel(r"$R$ [$m$]", fontsize=axis_tick_size)
    ax.set_title("Radial Displacement")
    # Only set legend for the last plot if needed

    # Set x-label only on the last axis
    axis[2].set_xlabel(r"$R$ [$m$]", fontsize=axis_tick_size)

    # Set minor ticks on for all axes
    for ax in axis:
        ax.minorticks_on()
    # Set x-ticks and y-ticks font size for all axes
    for ax in axis:
        ax.tick_params(axis="x", labelsize=axis_tick_size)
        ax.tick_params(axis="y", labelsize=axis_tick_size)
    plt.tight_layout()


def plot_corc_cable_geometry(
    axis,
    r_centre: float,
    z_centre: float,
    dia_croco_strand: float,
    dx_croco_strand_copper: float,
    dr_hts_tape: float,
    dx_croco_strand_tape_stack: float,
    n_croco_strand_hts_tapes: int,
    dx_hts_tape_rebco: float,
    dx_hts_tape_copper: float,
    dx_hts_tape_hastelloy: float,
    show_legend: bool = True,
):
    """Plot the geometry of a CroCo strand cable.

    Parameters
    ----------
    axis : matplotlib.axes._axes.Axes
        The matplotlib axis to plot on.
    r_centre : float
        Radial position of the strand centre (in meters).
    z_centre : float
        Vertical position of the strand centre (in meters).
    dia_croco_strand : float
        Diameter of the CroCo strand (in meters).
    dx_croco_strand_copper : float
        Thickness of the copper layer (in meters).
    dr_hts_tape : float
        Radius of the HTS tape stack (in meters).
    dx_croco_strand_tape_stack : float
        Height of the HTS tape stack (in meters).
    n_croco_strand_hts_tapes : int
        Number of HTS tape layers in the stack.
    """
    legend_label = None if show_legend else "_nolegend_"

    # Plot a circle with the given diameter and copper edges
    circle = Circle(
        (r_centre, z_centre),
        radius=(dia_croco_strand / 2),
        edgecolor="black",
        facecolor="#B87333",
        linewidth=0.5,
        label="Copper jacket" if show_legend else legend_label,
    )
    axis.add_patch(circle)

    # Plot an inner circle with copper edges
    circle = Circle(
        (r_centre, z_centre),
        radius=((dia_croco_strand / 2) - dx_croco_strand_copper),
        edgecolor="grey",
        facecolor="grey",
        linewidth=2,
        label="Solder" if show_legend else legend_label,
    )
    axis.add_patch(circle)

    # Plot a rectangular tape stack in the middle
    rect = Rectangle(
        (
            r_centre - dr_hts_tape / 2,
            z_centre - dx_croco_strand_tape_stack / 2,
        ),
        width=dr_hts_tape,
        height=dx_croco_strand_tape_stack,
        edgecolor="black",
        facecolor=None,
        linewidth=0.1,
        alpha=0.5,
        linestyle="--",
        label="HTS Tape Stack" if show_legend else legend_label,
    )
    axis.add_patch(rect)

    # Slice the tape stack into n_croco_strand_hts_tapes layers
    for i in range(int(n_croco_strand_hts_tapes)):
        y_start = (
            z_centre
            - (dx_croco_strand_tape_stack / 2)
            + i * (dx_croco_strand_tape_stack / n_croco_strand_hts_tapes)
        )
        plot_hts_tape_geometry(
            axis=axis,
            r_left=r_centre - (dr_hts_tape / 2),
            z_bottom=y_start,
            dr_hts_tape=dr_hts_tape,
            dx_hts_tape_rebco=dx_hts_tape_rebco,
            dx_hts_tape_copper=dx_hts_tape_copper,
            dx_hts_tape_hastelloy=dx_hts_tape_hastelloy,
            show_legend=False,
        )

    axis.set_xlim(-dia_croco_strand * 0.75, dia_croco_strand * 0.75)
    axis.set_ylim(-dia_croco_strand * 0.75, dia_croco_strand * 0.75)
    axis.set_aspect("equal", adjustable="datalim")
    axis.set_title("CroCo Strand Geometry")
    axis.grid(True)
    axis.set_xlabel("X-axis (m)")
    axis.set_ylabel("Y-axis (m)")
    axis.minorticks_on()
    if show_legend:
        axis.legend(loc="upper right")


def plot_tf_corc_cable_summary_box(axis, fig, mfile: MFile, scan: int):
    """Plot TF CORC cable summary box"""
    textstr_cable = (
        "$\\mathbf{CroCo \\ Cable:}$\n\nCable diameter:"
        f" {mfile.get('dia_tf_turn_croco_cable', scan=scan) * 1e3:,.4f}"
        " mm\nCopper width:"
        f" {mfile.get('dx_tf_croco_strand_copper', scan=scan) * 1e3:,.4f}"
        " mm\nDiameter of solder tape region:"
        f" {mfile.get('dia_tf_croco_strand_tape_region', scan=scan) * 1e3:,.4f}"
        " mm\nHeight of tape stack:"
        f" {mfile.get('dx_tf_croco_strand_tape_stack', scan=scan) * 1e3:,.4f}"
        " mm\nWidth of HTS tape / tape stack:"
        f" {mfile.get('dr_tf_hts_tape', scan=scan) * 1e3:,.4f} mm\nNumber of"
        " HTS tape layers:"
        f" {int(mfile.get('n_tf_croco_strand_hts_tapes', scan=scan))}\n\nTotal"
        " copper area:"
        f" {mfile.get('a_tf_croco_strand_copper_total', scan=scan) * 1e6:,.4f}"
        " mm²\nTotal hastelloy area:"
        f" {mfile.get('a_tf_croco_strand_hastelloy', scan=scan) * 1e6:,.4f}"
        " mm²\nTotal solder area:"
        f" {mfile.get('a_tf_croco_strand_solder', scan=scan) * 1e6:,.4f}"
        " mm²\nTotal superconductor area:"
        f" {mfile.get('a_tf_croco_strand_rebco', scan=scan) * 1e6:,.4f}"
        " mm²\nTotal strand area:"
        f" {mfile.get('a_tf_croco_strand', scan=scan) * 1e6:,.4f} mm²\n"
    )

    draw_text(
        axis,
        0.4,
        0.4,
        textstr_cable,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("#cccccc"),  # grayish color
    )


def plot_quench_time_evolution(
    tau_discharge: float,
    b_peak: float,
    f_a_cable_copper: float,
    f_a_cable_space_helium: float,
    temp_he_peak: float,
    temp_quench_max: float,
    cu_rrr: float,
    t_quench_detection: float,
    fluence: float,
    j_operating: float,
    a_tf_turn_cable_space: float,
    a_tf_turn: float,
    n_points: int = 500,
    axes_1: plt.Axes | None = None,
    axes_2: plt.Axes | None = None,
    show: bool = False,
) -> None:
    """Plots the time evolution of the quench model hotspot temperature and current.

    Visualises the adiabatic hotspot temperature rise and exponentially decaying
    current during a quench, highlighting the quench detection time.

    Parameters
    ----------
    tau_discharge:
        Quench discharge time constant [s].
    b_peak:
        Magnetic field at the peak point [T].
    f_a_cable_copper:
        Fraction of cable cross-section that is copper.
    f_a_cable_space_helium:
        Fraction of cable space occupied by helium.
    temp_he_peak:
        Peak helium temperature at quench initiation [K].
    temp_quench_max:
        Maximum allowed conductor temperature during quench [K].
    cu_rrr:
        Residual resistivity ratio of copper.
    t_quench_detection:
        Detection time delay [s].
    fluence:
        Neutron fluence [n/m²].
    j_operating:
        Operating current density [A/m²] to compare against the quench protection limit.
    a_tf_turn_cable_space:
        Area of the TF turn cable space [m²].
    a_tf_turn:
        Area of the TF turn [m²].
    n_points:
        Number of time points for the plot.
    axes_1:
        Optional axis for the current density panel.
    axes_2:
        Optional axis for the hotspot temperature panel.
    show:
        Whether to display the plot with Matplotlib. Defaults to False to avoid
        GUI backend warnings in non-interactive environments.

    Raises
    ------
    ValueError
        If only one set of axes is provided, instead of both or neither
    """
    figure = None
    if axes_1 is None and axes_2 is None:
        figure, (axes_1, axes_2) = plt.subplots(2, 1, sharex=True)
    elif axes_1 is None or axes_2 is None:
        msg = "Both axes_1 and axes_2 must be provided together, or neither."
        raise ValueError(msg)

    fluence = np.clip(fluence, 0.0, 1.5e23)

    j_max = (
        a_tf_turn_cable_space / a_tf_turn
    ) * quench.calculate_quench_protection_current_density(
        tau_discharge=tau_discharge,
        b_peak=b_peak,
        f_a_cable_copper=f_a_cable_copper,
        f_a_cable_space_helium=f_a_cable_space_helium,
        temp_he_peak=temp_he_peak,
        temp_quench_max=temp_quench_max,
        cu_rrr=cu_rrr,
        t_quench_detection=t_quench_detection,
        fluence=fluence,
    )

    fluence_1e23 = 1e23
    j_max_1e23 = (
        a_tf_turn_cable_space / a_tf_turn
    ) * quench.calculate_quench_protection_current_density(
        tau_discharge=tau_discharge,
        b_peak=b_peak,
        f_a_cable_copper=f_a_cable_copper,
        f_a_cable_space_helium=f_a_cable_space_helium,
        temp_he_peak=temp_he_peak,
        temp_quench_max=temp_quench_max,
        cu_rrr=cu_rrr,
        t_quench_detection=t_quench_detection,
        fluence=fluence_1e23,
    )

    # Time axis: from 0 to ~4 time constants after discharge begins at detection.
    # This ensures later annotations/interpolations at t_quench_detection + tau_discharge
    # and beyond remain within the sampled domain.
    t_end = max(4.0 * tau_discharge, t_quench_detection + 4.0 * tau_discharge)
    times = np.linspace(0.0, t_end, n_points)

    # Current density decays exponentially after detection
    decay = np.exp(-(times - t_quench_detection) / tau_discharge)

    j_profile_required, j_profile_required_1e23, j_profile_real = [
        np.where(times < t_quench_detection, j0, j0 * decay)
        for j0 in (j_max, j_max_1e23, j_operating)
    ]

    # Adiabatic hotspot temperature: integrate heat balance over time
    # T(t) is found by inverting: integral_{T0}^{T(t)} [sum(rho*cp)] / rho_cu dT =
    # integral_0^t J² dt
    # We accumulate the (∫J² dt) and map it to temperature via the precomputed integral.
    f_cu_cable = (1.0 - f_a_cable_space_helium) * f_a_cable_copper
    f_sc_cable = (1.0 - f_a_cable_space_helium) * (1.0 - f_a_cable_copper)

    # Build a temperature lookup: cumulative integral from t_he_peak to T
    temp_array, cum_integral = quench._build_cumulative_quench_integral(
        temp_he_peak=temp_he_peak,
        temp_quench_max=temp_quench_max,
        field=b_peak,
        rrr=cu_rrr,
        fluence=fluence,
        f_a_cable_space_helium=f_a_cable_space_helium,
        f_cu_cable=f_cu_cable,
        f_sc_cable=f_sc_cable,
    )
    temp_array_1e23, cum_integral_1e23 = quench._build_cumulative_quench_integral(
        temp_he_peak=temp_he_peak,
        temp_quench_max=temp_quench_max,
        field=b_peak,
        rrr=cu_rrr,
        fluence=fluence_1e23,
        f_a_cable_space_helium=f_a_cable_space_helium,
        f_cu_cable=f_cu_cable,
        f_sc_cable=f_sc_cable,
    )

    # Numerically integrate J² dt over time to get MIIT (Mega-Ampere²-seconds) at
    # each time step
    dt = times[1] - times[0]
    miit_required = np.cumsum(j_profile_required**2) * dt
    miit_required_1e23 = np.cumsum(j_profile_required_1e23**2) * dt
    miit_real = np.cumsum(j_profile_real**2) * dt

    # Convert the cable-space thermal integral to winding-pack basis to match
    # j_profile_*.
    area_ratio = a_tf_turn_cable_space / a_tf_turn
    scaled_integral = (area_ratio**2) * f_cu_cable * cum_integral
    scaled_integral_1e23 = (area_ratio**2) * f_cu_cable * cum_integral_1e23
    hotspot_temp_required = np.interp(miit_required, scaled_integral, temp_array)
    hotspot_temp_required_1e23 = np.interp(
        miit_required_1e23, scaled_integral_1e23, temp_array_1e23
    )
    hotspot_temp_real = np.interp(miit_real, scaled_integral, temp_array)

    # --- Current density panel ---
    axes_1.plot(
        times,
        j_profile_required,
        color="darkorange",
        linewidth=2,
        label=(
            f"Max allowed current density for protection (fluence = {fluence:.2e} n/m²)"
        ),
    )
    axes_1.plot(
        times,
        j_profile_required_1e23,
        color="darkorange",
        linewidth=2,
        linestyle="--",
        label=("Max allowed current density for protection (fluence = 1e23 n/m²)"),
    )
    axes_1.plot(
        times,
        j_profile_real,
        color="blue",
        linewidth=2,
        label="Operating current density",
    )
    axes_1.axvline(
        t_quench_detection,
        color="crimson",
        linestyle="--",
        linewidth=1.5,
        label=f"Detection time ({t_quench_detection:.1f} s)",
    )
    axes_1.axvspan(
        0,
        t_quench_detection,
        alpha=0.08,
        color="crimson",
        label="Pre-detection phase",
    )
    axes_1.set_ylabel("Current density [A/m²]")
    axes_1.legend(fontsize=9)
    axes_1.grid(True, alpha=0.3)
    axes_1.set_title(
        "TF Coil Quench Protection: Current Density and Hotspot Temperature Evolution"
    )

    # --- Temperature panel ---
    axes_2.plot(
        times,
        hotspot_temp_required,
        color="darkorange",
        linewidth=2,
        label=(
            f"Hotspot temperature at protection limit (fluence = {fluence:.2e} n/m²)"
        ),
    )
    axes_2.plot(
        times,
        hotspot_temp_required_1e23,
        color="darkorange",
        linewidth=2,
        linestyle="--",
        label="Hotspot temperature at protection limit (fluence = 1e23 n/m²)",
    )
    axes_2.plot(
        times,
        hotspot_temp_real,
        color="blue",
        linewidth=2,
        label="Operating hotspot temperature",
    )

    axes_2.axvline(
        t_quench_detection,
        color="crimson",
        linestyle="--",
        linewidth=1.5,
        label=f"$t_{{\\text{{detect}}}}$ ({t_quench_detection:.2f} s)",
    )
    axes_2.axvspan(0, t_quench_detection, alpha=0.08, color="crimson")
    axes_2.axhline(
        temp_quench_max,
        color="grey",
        linestyle=":",
        linewidth=1.5,
        label=f"$T_{{\\text{{max}}}}$ = {temp_quench_max} K",
    )
    axes_2.set_xlabel("Time [s]")
    axes_2.set_ylabel("Temperature [K]")
    axes_2.legend(fontsize=9)
    axes_2.grid(True, alpha=0.3)

    # Mark tau_discharge after detection time with vertical and horizontal lines
    tau_time = t_quench_detection + tau_discharge
    tau_j = j_max * np.exp(
        -1
    )  # current density at t = t_quench_detection + tau_discharge
    tau_temp = float(np.interp(tau_time, times, hotspot_temp_required))

    for ax, val, label in [
        (
            axes_1,
            tau_j,
            f"$J$ at $\\tau_{{\\text{{discharge}}}}$ ({tau_j:.2e} A/m²)",
        ),
        (
            axes_2,
            tau_temp,
            f"$T$ at $\\tau_{{\\text{{discharge}}}}$ ({tau_temp:.1f} K)",
        ),
    ]:
        ax.axvline(
            tau_time,
            color="forestgreen",
            linestyle="--",
            linewidth=1.5,
            label=(
                "$t_{\\text{detect}} + \\tau_{\\text{discharge}}$"
                f" ({tau_time:.2f} s)"
            ),
        )
        ax.axhline(
            val,
            color="forestgreen",
            linestyle=":",
            linewidth=1.5,
            label=label,
        )
    axes_1.legend(fontsize=9)
    axes_1.minorticks_on()
    axes_2.legend(fontsize=9)
    axes_2.minorticks_on()

    if figure is not None:
        figure.tight_layout()
    else:
        plt.tight_layout()

    if show:
        plt.show()


__all__ = [
    "TF_outboard",
    "plot_corc_cable_geometry",
    "plot_quench_time_evolution",
    "plot_resistive_tf_info",
    "plot_resistive_tf_wp",
    "plot_superconducting_tf_wp",
    "plot_tf_cable_in_conduit_turn",
    "plot_tf_coil_structure",
    "plot_tf_coils",
    "plot_tf_corc_cable_summary_box",
    "plot_tf_croco_turn",
    "plot_tf_stress",
]
