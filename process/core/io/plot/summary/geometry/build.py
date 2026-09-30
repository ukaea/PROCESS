"""Geometry functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

from process.core.io.plot.summary.common import (
    setup_axis,
)
from process.core.io.plot.summary.constants import (
    BLANKET_COLOUR,
    CSCOMPRESSION_COLOUR,
    FIRSTWALL_COLOUR,
    PLASMA_COLOUR,
    RADIAL_BUILD,
    SHIELD_COLOUR,
    SOLENOID_COLOUR,
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    VESSEL_COLOUR,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)
from process.core.io.plot.summary.reporting import (
    plot_info,
)
from process.data_structure.build_variables import TFCSRadialConfiguration

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def cumulative_radial_build(section, mfile: MFile, scan: int):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Parameters
    ----------
    section :
        section of the radial build to go up to
    mfile :
        MFILE data object
    scan :
        scan number to use

    Returns
    -------
    :
        cumulative_build:cumulative radial build up to section given
    """
    complete = False
    cumulative_build = 0
    for item in RADIAL_BUILD:
        if item in {"rminori", "rminoro"}:
            cumulative_build += mfile.get("rminor", scan=scan)
        elif item in {"vvblgapi", "vvblgapo"}:
            cumulative_build += mfile.get("dr_shld_blkt_gap", scan=scan)
        elif "dr_vv_inboard" in item:
            cumulative_build += mfile.get("dr_vv_inboard", scan=scan)
        elif "dr_vv_outboard" in item:
            cumulative_build += mfile.get("dr_vv_outboard", scan=scan)
        else:
            cumulative_build += mfile.get(item, scan=scan)
        if item == section:
            complete = True
            break

    if complete is False:
        print("radial build parameter ", section, " not found")
    return cumulative_build


def cumulative_radial_build2(section, mfile: MFile, scan: int):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Parameters
    ----------
    section :
        section of the radial build to go up to
    mfile :
        MFILE data object
    scan :
        scan number to use

    Returns
    -------
    :
        cumulative_build --> cumulative radial build up to and including
        section given
        previous         --> cumulative radial build up to section given
    """
    cumulative_build = 0
    build = 0
    for item in RADIAL_BUILD:
        if item in {"rminori", "rminoro"}:
            build = mfile.get("rminor", scan=scan)
        elif item in {"vvblgapi", "vvblgapo"}:
            build = mfile.get("dr_shld_blkt_gap", scan=scan)
        elif "dr_vv_inboard" in item:
            build = mfile.get("dr_vv_inboard", scan=scan)
        elif "dr_vv_outboard" in item:
            build = mfile.get("dr_vv_outboard", scan=scan)
        else:
            build = mfile.get(item, scan=scan)
        cumulative_build += build
        if item == section:
            break
    previous = cumulative_build - build
    return (cumulative_build, previous)


def plot_geometry_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot geometry info

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    draw_text(axis, -0.05, 1, "Geometry:", ha="left", va="center")
    setup_axis(axis, xmin, xmax, ymin, ymax)

    in_blanket_thk = mfile.get("dr_shld_inboard", scan=scan) + mfile.get(
        "dr_blkt_inboard", scan=scan
    )
    out_blanket_thk = mfile.get("dr_shld_outboard", scan=scan) + mfile.get(
        "dr_blkt_outboard", scan=scan
    )

    data = [
        ("rmajor", "$R_0$", "m"),
        ("rminor", "a", "m"),
        ("aspect", "A", ""),
        ("kappa95", r"$\kappa_{95}$", ""),
        ("triang95", r"$\delta_{95}$", ""),
        ("a_plasma_surface", "Plasma surface area", "m$^2$"),
        ("a_plasma_poloidal", "Plasma cross-sectional area", "m$^2$"),
        ("vol_plasma", "Plasma volume", "m$^3$"),
        ("n_tf_coils", "No. of TF coils", ""),
        (in_blanket_thk, "Inboard blanket+shield", "m"),
        ("dr_inboard_build", "Inboard build thickness", "m"),
        (out_blanket_thk, "Outboard blanket+shield", "m"),
    ]

    plot_info(axis, data, mfile, scan)


def plot_radial_build(axis: plt.Axes, mfile: MFile, colour_scheme: Literal[1, 2]):
    """Plots the radial build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's radial
    build
    (such as central solenoid, toroidal field coils, vacuum vessel, shields, blankets,
    etc.)
    as a horizontal stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile`, and each segment is color-coded and labeled accordingly.

    If the toroidal field coil is inside the central solenoid (as indicated by the
    "i_tf_inside_cs" flag in `mfile`), the order and labels of the components are
    adjusted accordingly.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The matplotlib axis on which to plot the radial build.
    mfile : MFile
        An object containing the machine build data, with required fields for each
        radial component and the "i_tf_inside_cs" flag.
    colour_scheme:

    Notes
    -----
    This function modifies the provided axis in-place and does not return a value.
    - Components with zero thickness are omitted from the plot.
    - The legend displays the name and thickness (in meters) of each component.
    """
    radial_variables = [
        "dr_bore",
        "dr_cs",
        "dr_cs_precomp",
        "dr_cs_tf_gap",
        "dr_tf_inboard",
        "dr_tf_shld_gap",
        "dr_shld_thermal_inboard",
        "dr_shld_vv_gap_inboard",
        "dr_vv_inboard",
        "dr_shld_inboard",
        "dr_shld_blkt_gap",
        "dr_blkt_inboard",
        "dr_fw_inboard",
        "dr_fw_plasma_gap_inboard",
        "rminor",
        "dr_fw_plasma_gap_outboard",
        "dr_fw_outboard",
        "dr_blkt_outboard",
        "dr_shld_blkt_gap",
        "dr_vv_outboard",
        "dr_shld_outboard",
        "dr_shld_vv_gap_outboard",
        "dr_shld_thermal_outboard",
        "dr_tf_shld_gap",
        "dr_tf_outboard",
    ]
    if int(mfile.get("i_tf_inside_cs", scan=-1)) == TFCSRadialConfiguration.TF_INSIDE_CS:
        radial_variables[1] = "dr_tf_inboard"
        radial_variables[2] = "dr_cs_tf_gap"
        radial_variables[3] = "dr_cs"
        radial_variables[4] = "dr_cs_precomp"
        radial_variables[5] = "dr_tf_shld_gap"

    radial_build = [[mfile.get(rl, scan=-1) for rl in radial_variables]]

    radial_build = np.array(radial_build)

    for kk in range(radial_build.shape[0]):
        radial_build[kk, 14] *= 2.0

    radial_build = np.transpose(radial_build)
    # ====================

    radial_labels = [
        "Machine Bore",
        "Central Solenoid",
        "CS precompression",
        "CS Coil gap",
        "TF Coil Inboard Leg",
        "TF Coil gap",
        "Inboard Thermal Shield",
        "Gap",
        "Inboard VV",
        "Inboard Shield",
        "Gap",
        "Inboard Blanket",
        "Inboard First Wall",
        "Inboard SOL",
        "Plasma",
        "Outboard SOL",
        "Outboard First Wall",
        "Outboard Blanket",
        "Gap",
        "Outboard VV",
        "Outboard Shield",
        "Gap",
        "Outboard Thermal Shield",
        "Gap",
        "TF Coil Outboard Leg",
    ]
    if int(mfile.get("i_tf_inside_cs", scan=-1)) == TFCSRadialConfiguration.TF_INSIDE_CS:
        radial_labels[1] = "TF Coil Inboard Leg"
        radial_labels[2] = "CS Coil gap"
        radial_labels[3] = "Central Solenoid"
        radial_labels[4] = "CS precompression"
        radial_labels[5] = "TF Coil gap"

    radial_color = [
        "white",
        SOLENOID_COLOUR[colour_scheme - 1],
        CSCOMPRESSION_COLOUR[colour_scheme - 1],
        "white",
        (
            TFC_COLOUR[colour_scheme - 1]
            if mfile.get("i_tf_sup", scan=-1) != 0
            else "#b87333"
        ),
        "white",
        THERMAL_SHIELD_COLOUR[colour_scheme - 1],
        "white",
        VESSEL_COLOUR[colour_scheme - 1],
        SHIELD_COLOUR[colour_scheme - 1],
        "white",
        BLANKET_COLOUR[colour_scheme - 1],
        FIRSTWALL_COLOUR[colour_scheme - 1],
        "white",
        PLASMA_COLOUR[colour_scheme - 1],
        "white",
        FIRSTWALL_COLOUR[colour_scheme - 1],
        BLANKET_COLOUR[colour_scheme - 1],
        "white",
        VESSEL_COLOUR[colour_scheme - 1],
        SHIELD_COLOUR[colour_scheme - 1],
        "white",
        THERMAL_SHIELD_COLOUR[colour_scheme - 1],
        "white",
        (
            TFC_COLOUR[colour_scheme - 1]
            if mfile.get("i_tf_sup", scan=-1) != 0
            else "#b87333"
        ),
    ]
    if int(mfile.get("i_tf_inside_cs", scan=-1)) == TFCSRadialConfiguration.TF_INSIDE_CS:
        radial_color[1] = (
            TFC_COLOUR[colour_scheme - 1]
            if mfile.get("i_tf_sup", scan=-1) != 0
            else "#b87333"
        )
        radial_color[2] = "white"
        radial_color[3] = SOLENOID_COLOUR[colour_scheme - 1]
        radial_color[4] = CSCOMPRESSION_COLOUR[colour_scheme - 1]
        radial_color[5] = "white"

    lower = np.zeros(radial_build.shape[1])
    for kk in range(radial_build.shape[0]):
        axis.barh(
            0,
            radial_build[kk, :],
            left=lower,
            height=0.8,
            label=(
                f"{radial_labels[kk]}\n[{radial_variables[kk]}]\n{radial_build[kk][0]:.3f} m"  # noqa: E501
            ),
            color=radial_color[kk],
            edgecolor="black",
            linewidth=0.05,
        )
        lower += radial_build[kk, :]

    axis.set_yticks([])

    axis.legend(
        bbox_to_anchor=(0.5, -0.1),
        loc="upper center",
        ncol=5,
    )
    # Plot a vertical dashed line at rmajor
    axis.axvline(
        mfile.get("rmajor", scan=-1),
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Major Radius $R_0$",
    )
    axis.minorticks_on()
    axis.set_xlabel("Radius [m]")


__all__ = [
    "cumulative_radial_build",
    "cumulative_radial_build2",
    "plot_geometry_info",
    "plot_radial_build",
]
