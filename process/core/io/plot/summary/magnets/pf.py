"""Magnets functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
from matplotlib import patches

from process.core.io.plot.summary.constants import CSCOMPRESSION_COLOUR, SOLENOID_COLOUR
from process.core.io.plot.summary.radial_build import cumulative_radial_build2
from process.models.geometry.pfcoil import pfcoil_geometry

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_pf_coils(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    mirror_negative_x: bool = False,
):
    """Function to plot PF coils

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

    coils_r = []
    coils_z = []
    coils_dr = []
    coils_dz = []
    coil_text = []

    dr_cs_bore = mfile.get("dr_cs_bore", scan=scan)
    dr_cs = mfile.get("dr_cs", scan=scan)
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)

    # Number of coils, both PF and CS
    number_of_coils = 0
    for item in mfile.data:
        if "r_pf_coil_middle[" in item:
            number_of_coils += 1

    # Check for Central Solenoid
    iohcl = mfile.get("iohcl", scan=scan) if "iohcl" in mfile.data else 1

    # If Central Solenoid present, ignore last entry in for loop
    # The last entry will be the OH coil in this case
    noc = number_of_coils - 1 if iohcl == 1 else number_of_coils

    for coil in range(noc):
        coils_r.append(mfile.get(f"r_pf_coil_middle[{coil + 1:01}]", scan=scan))
        coils_z.append(mfile.get(f"z_pf_coil_middle[{coil + 1:01}]", scan=scan))
        coils_dr.append(mfile.get(f"pfdr({coil + 1:01})", scan=scan))
        coils_dz.append(mfile.get(f"pfdz({coil + 1:01})", scan=scan))
        coil_text.append(str(coil + 1))

    r_points, z_points, central_coil = pfcoil_geometry(
        coils_r=coils_r,
        coils_z=coils_z,
        coils_dr=coils_dr,
        coils_dz=coils_dz,
        dr_cs_bore=dr_cs_bore,
        dr_cs=dr_cs,
        ohdz=dz_cs_full,
    )

    # Plot CS compression structure
    r_precomp_outer, r_precomp_inner = cumulative_radial_build2(
        "dr_cs_precomp", mfile, scan
    )
    axis.add_patch(
        patches.Rectangle(
            xy=(x_scale * r_precomp_inner, central_coil.anchor_z),
            width=(x_scale * (r_precomp_outer - r_precomp_inner)),
            height=central_coil.height,
            facecolor=CSCOMPRESSION_COLOUR[colour_scheme - 1],
        )
    )

    # Get axis height for fontsize scaling
    axis_height = (
        axis
        .get_window_extent()
        .transformed(axis.figure.dpi_scale_trans.inverted())
        .height
    )

    for i in range(len(coils_r)):
        mirrored_r_points = [x_scale * r for r in r_points[i]]
        axis.plot(mirrored_r_points, z_points[i], color="black")
        # Scale fontsize relative to axis height and coil size
        fontsize = max(6, axis_height * abs(coils_dr[i] * coils_dz[i]) * 1.5)
        axis.text(
            x_scale * coils_r[i],
            coils_z[i] - 0.05,
            coil_text[i],
            ha="center",
            va="center",
            fontsize=fontsize,
        )
    axis.add_patch(
        patches.Rectangle(
            xy=(x_scale * central_coil.anchor_x, central_coil.anchor_z),
            width=x_scale * central_coil.width,
            height=central_coil.height,
            facecolor=SOLENOID_COLOUR[colour_scheme - 1],
            edgecolor="black",
            linewidth=1,
        )
    )
    axis.add_patch(
        patches.Rectangle(
            xy=(0.0, central_coil.anchor_z),
            width=x_scale * central_coil.anchor_x,
            height=central_coil.height,
            facecolor="grey",
            alpha=0.5,
        )
    )


def plot_pf_dimensions(
    axis: plt.Axes, mfile: MFile, scan: int, colour_scheme: Literal[1, 2] = 1
) -> None:
    """Plot the PF coil dimensions on the given axis."""
    r_pf_coil_middle = []
    z_pf_coil_middle = []
    radial_thicknesses = []
    vertical_thicknesses = []
    iohcl = mfile.get("iohcl", scan=scan) if "iohcl" in mfile.data else 1
    x = 1 if iohcl == 0 else 2
    for coil in range(int(mfile.get("n_pf_cs_plasma_circuits", scan=scan) - x)):
        r_pf_coil_middle.append(mfile.get(f"r_pf_coil_middle[{coil + 1}]", scan=scan))
        z_pf_coil_middle.append(mfile.get(f"z_pf_coil_middle[{coil + 1}]", scan=scan))
        radial_thicknesses.append(mfile.get(f"pfdr({coil + 1})", scan=scan))
        vertical_thicknesses.append(mfile.get(f"pfdz({coil + 1})", scan=scan))

    plot_pf_coils(axis=axis, mfile=mfile, scan=scan, colour_scheme=colour_scheme)

    if r_pf_coil_middle:
        for r_middle, z_middle, dr_coil, dz_coil in zip(
            r_pf_coil_middle,
            z_pf_coil_middle,
            radial_thicknesses,
            vertical_thicknesses,
            strict=False,
        ):
            half_radial_thickness = dr_coil / 2
            half_vertical_thickness = dz_coil / 2
            coil_left = r_middle - half_radial_thickness
            coil_right = r_middle + half_radial_thickness
            coil_bottom = z_middle - half_vertical_thickness
            coil_top = z_middle + half_vertical_thickness

            for x_position in (coil_left, r_middle, coil_right):
                axis.axvline(
                    x=x_position,
                    color="r",
                    linewidth=0.8,
                    linestyle="--" if x_position == r_middle else "-",
                    alpha=0.3,
                    zorder=4,
                )

            for y_position in (coil_bottom, z_middle, coil_top):
                axis.axhline(
                    y=y_position,
                    xmax=coil_left,
                    color="r",
                    linewidth=0.8,
                    linestyle="--" if y_position == z_middle else "-",
                    alpha=0.3,
                    zorder=4,
                )

            axis.annotate(
                f"({r_middle:.3f}, {z_middle:.3f})",
                xy=(coil_left * 0.925, z_middle),
                ha="right",
                va="center",
                fontsize=8,
                zorder=6,
                bbox={
                    "boxstyle": "round,pad=0.2",
                    "fc": "white",
                    "alpha": 1.0,
                    "ec": "none",
                },
            )

            radial_arrow_y = coil_bottom if z_middle < 0 else coil_top
            radial_label_offset = (0, -24) if z_middle < 0 else (0, 4)
            radial_label_va = "top" if z_middle < 0 else "bottom"
            axis.annotate(
                "",
                xy=(coil_left, radial_arrow_y),
                xytext=(coil_right, radial_arrow_y),
                arrowprops={
                    "arrowstyle": "<->",
                    "linewidth": 0.8,
                    "color": "red",
                    "shrinkA": 0,
                    "shrinkB": 0,
                },
                zorder=5,
            )
            axis.annotate(
                f"ΔR={abs(dr_coil):.3f}",
                xy=(r_middle, radial_arrow_y),
                xytext=radial_label_offset,
                textcoords="offset points",
                ha="center",
                va=radial_label_va,
                fontsize=8,
                zorder=6,
                bbox={
                    "boxstyle": "round,pad=0.2",
                    "fc": "white",
                    "alpha": 1.0,
                    "ec": "none",
                },
            )

            vertical_arrow_x = coil_right
            axis.annotate(
                "",
                xy=(vertical_arrow_x, coil_bottom),
                xytext=(vertical_arrow_x, coil_top),
                arrowprops={
                    "arrowstyle": "<->",
                    "linewidth": 0.8,
                    "color": "red",
                    "shrinkA": 0,
                    "shrinkB": 0,
                },
            )
            axis.annotate(
                f"ΔZ={abs(dz_coil):.3f}",
                xy=(vertical_arrow_x, z_middle),
                xytext=(4, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=8,
                zorder=6,
                bbox={
                    "boxstyle": "round,pad=0.2",
                    "fc": "white",
                    "alpha": 1.0,
                    "ec": "none",
                },
            )

    axis.set_title("PF Coil Dimensions")
    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_xlim(left=0.0)
    axis.minorticks_on()
    axis.grid(True, alpha=0.3)
    axis.set_aspect("equal", adjustable="box")


__all__ = ["plot_pf_coils", "plot_pf_dimensions"]
