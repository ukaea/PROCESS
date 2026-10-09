"""Reporting functions for PROCESS summary plots."""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING, Any, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from process.core.io.mfile import MFile, MFileErrorClass
from process.core.io.plot.summary.common import (
    box_style,
    setup_axis,
)
from process.core.io.plot.summary.geometry.poloidal import (
    poloidal_cross_section,
)
from process.core.io.plot.summary.plasma.physics import (
    plot_plasma,
)
from process.core.io.plot.summary.reporting.text import plot_info
from process.data_structure.numerics import FiguresOfMerit, PROCESSRunMode
from process.data_structure.physics_variables import DivertorNumberModels

if TYPE_CHECKING:
    from process.core.io.plot.summary.reporting.misc import (
        RadialBuild,
    )


def plot_header(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot header info: date, rutitle etc

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    setup_axis(axis, xmin=0, xmax=1, ymin=-16, ymax=1)

    data2 = [
        (f"!{mfile.get('runtitle', scan=-1)}", "Run title", ""),
        (f"!{mfile.get('procver', scan=-1)}", "PROCESS Version", ""),
        (f"!{mfile.get('date', scan=-1)}", "Date:", ""),
        (f"!{mfile.get('time', scan=-1)}", "Time:", ""),
        (f"!{mfile.get('username', scan=-1)}", "User:", ""),
        (
            ("!Evaluation", "Run type", "")
            if isinstance(mfile.data["i_figure_merit"], MFileErrorClass)
            else (
                (
                    f"!{FiguresOfMerit(abs(int(mfile.get('i_figure_merit', scan=-1)))).description}"  # noqa: E501
                ),
                "Optimising:",
                "",
            )
        ),
    ]

    axis.text(-0.05, 4.0, "Colour Legend:", ha="left", va="center")
    axis.text(
        0.0,
        3.0,
        "ITR --> Iteration variable",
        color="red",
        ha="left",
        va="center",
    )
    axis.text(
        0.0,
        2.0,
        "OP  --> Output variable",
        color="blue",
        ha="left",
        va="center",
    )

    H = mfile.get("f_nd_impurity_electrons(01)", scan=scan)
    He = mfile.get("f_nd_impurity_electrons(02)", scan=scan)
    Be = mfile.get("f_nd_impurity_electrons(03)", scan=scan)
    C = mfile.get("f_nd_impurity_electrons(04)", scan=scan)
    N = mfile.get("f_nd_impurity_electrons(05)", scan=scan)
    O = mfile.get("f_nd_impurity_electrons(06)", scan=scan)  # noqa: E741
    Ne = mfile.get("f_nd_impurity_electrons(07)", scan=scan)
    Si = mfile.get("f_nd_impurity_electrons(08)", scan=scan)
    Ar = mfile.get("f_nd_impurity_electrons(09)", scan=scan)
    Fe = mfile.get("f_nd_impurity_electrons(10)", scan=scan)
    Ni = mfile.get("f_nd_impurity_electrons(11)", scan=scan)
    Kr = mfile.get("f_nd_impurity_electrons(12)", scan=scan)
    Xe = mfile.get("f_nd_impurity_electrons(13)", scan=scan)
    W = mfile.get("f_nd_impurity_electrons(14)", scan=scan)

    data = [("", "", ""), ("", "", "")]
    count = 0

    data = [*data, (H, "D + T", "")]
    count += 1

    data = [*data, (He, "He", "")]
    count += 1
    if Be > 1e-10:
        data = [*data, (Be, "Be", "")]
        count += +1
    if C > 1e-10:
        data = [*data, (C, "C", "")]
        count += 1
    if N > 1e-10:
        data = [*data, (N, "N", "")]
        count += 1
    if O > 1e-10:
        data = [*data, (O, "O", "")]
        count += 1
    if Ne > 1e-10:
        data = [*data, (Ne, "Ne", "")]
        count += 1
    if Si > 1e-10:
        data = [*data, (Si, "Si", "")]
        count += 1
    if Ar > 1e-10:
        data = [*data, (Ar, "Ar", "")]
        count += 1
    if Fe > 1e-10:
        data = [*data, (Fe, "Fe", "")]
        count += 1
    if Ni > 1e-10:
        data = [*data, (Ni, "Ni", "")]
        count += 1
    if Kr > 1e-10:
        data = [*data, (Kr, "Kr", "")]
        count += 1
    if Xe > 1e-10:
        data = [*data, (Xe, "Xe", "")]
        count += 1
    if W > 1e-10:
        data = [*data, (W, "W", "")]
        count += 1

    if count > 11:
        data = [
            ("", "", ""),
            ("", "", ""),
            ("", "More than 11 impurities", ""),
        ]
    else:
        axis.text(-0.05, -6.4, "Plasma composition:", ha="left", va="center")
        axis.text(
            -0.05,
            -7.2,
            "Number densities relative to electron density:",
            ha="left",
            va="center",
        )
    data2 += data

    plot_info(axis, data2, mfile, scan)


def plot_separatrix_power_split(axis: plt.Axes, mfile: MFile, scan: int, colour_scheme):
    """Plot separatrix power split fractions as a bar chart."""
    plot_plasma(axis=axis, mfile=mfile, scan=scan, colour_scheme=colour_scheme)
    rmajor, rminor, kappa, dr_sep = mfile.get_variables(
        "rmajor",
        "rminor",
        "kappa",
        "dr_plasma_outboard_midplane_separatrix_separation",
        scan=scan,
    )

    plasma_scale = max(rminor, abs(kappa * rminor), 1e-6)
    scale_factor = min(max(plasma_scale / 2.0, 0.7), 1.0)
    text_fontsize = 9 * scale_factor

    is_double_null = (
        DivertorNumberModels(mfile.get("i_single_null", scan=scan))
        == DivertorNumberModels.DOUBLE_NULL
    )
    p_sep = mfile.get("p_plasma_separatrix_mw", scan=scan)
    f_outboard = mfile.get("f_p_div_outboard_separatrix", scan=scan)
    f_inboard = mfile.get("f_p_div_inboard_separatrix", scan=scan)
    p_outboard = p_sep * f_outboard
    p_inboard = p_sep * f_inboard
    p_lower_inboard = mfile.get("p_div_lower_inboard_separatrix_mw", scan=scan)
    p_lower_outboard = mfile.get("p_div_lower_outboard_separatrix_mw", scan=scan)

    power_values = [
        p_sep,
        p_outboard,
        p_inboard,
        p_lower_inboard,
        p_lower_outboard,
    ]

    p_upper_inboard = None
    p_upper_outboard = None
    if is_double_null:
        p_upper_inboard = mfile.get("p_div_upper_inboard_separatrix_mw", scan=scan)
        p_upper_outboard = mfile.get("p_div_upper_outboard_separatrix_mw", scan=scan)
        power_values.extend([p_upper_inboard, p_upper_outboard])

    power_min = min(power_values)
    power_max = max(power_values)
    colour_map = mpl.colormaps["coolwarm"]

    def make_bbox_props(power: float) -> dict[str, Any]:
        norm_power = (
            1.0
            if np.isclose(power_max, power_min)
            else (power - power_min) / (power_max - power_min)
        )
        return {
            "boxstyle": f"round,pad={0.3 * scale_factor:.3f}",
            "facecolor": colour_map(norm_power),
            "alpha": 1.0,
            "linewidth": 2 * scale_factor,
            "edgecolor": "black",
        }

    centre_pos = (rmajor, 0.0)
    outboard_pos = (rmajor + rminor, 0.0)
    inboard_pos = (rmajor - rminor, 0.0)
    lower_inboard_pos = (rmajor - rminor, -kappa * rminor)
    lower_outboard_pos = (rmajor + rminor, -kappa * rminor)
    upper_inboard_pos = (rmajor - rminor, kappa * rminor)
    upper_outboard_pos = (rmajor + rminor, kappa * rminor)

    axis.text(
        *centre_pos,
        f"$P_{{\\mathrm{{sep}}}} = {p_sep:.3f}$ MW",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=make_bbox_props(p_sep),
        zorder=101,
    )
    axis.text(
        *outboard_pos,
        f"$f_{{\\mathrm{{outboard}}}} = {f_outboard:.3f}$\n"
        f"$\\Delta r_{{\\mathrm{{sep}}}} = {dr_sep:.3f}$ m",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=make_bbox_props(p_outboard),
        zorder=101,
    )
    axis.text(
        *inboard_pos,
        f"$f_{{\\mathrm{{inboard}}}} = {f_inboard:.3f}$",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=make_bbox_props(p_inboard),
        zorder=101,
    )
    axis.text(
        *lower_inboard_pos,
        "$f_{\\mathrm{lower\\ inboard}} ="
        f" {mfile.get('f_p_div_lower_inboard_separatrix', scan=scan):.3f}$\n$P_{{\\mathrm{{lower\\"  # noqa: E501
        f" inboard}}}} = {p_lower_inboard:.3f}$ MW",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=make_bbox_props(p_lower_inboard),
        zorder=101,
    )
    axis.text(
        *lower_outboard_pos,
        "$f_{\\mathrm{lower\\ outboard}} ="
        f" {mfile.get('f_p_div_lower_outboard_separatrix', scan=scan):.3f}$\n$P_{{\\mathrm{{lower\\"  # noqa: E501
        f" outboard}}}} = {p_lower_outboard:.3f}$ MW",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=make_bbox_props(p_lower_outboard),
        zorder=101,
    )
    if is_double_null:
        axis.text(
            *upper_inboard_pos,
            "$f_{\\mathrm{upper\\ inboard}} ="
            f" {mfile.get('f_p_div_upper_inboard_separatrix', scan=scan):.3f}$\n$P_{{\\mathrm{{upper\\"  # noqa: E501
            f" inboard}}}} = {p_upper_inboard:.3f}$ MW",
            fontsize=text_fontsize,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=make_bbox_props(p_upper_inboard),
            zorder=101,
        )
        axis.text(
            *upper_outboard_pos,
            "$f_{\\mathrm{upper\\ outboard}} ="
            f" {mfile.get('f_p_div_upper_outboard_separatrix', scan=scan):.3f}$\n$P_{{\\mathrm{{upper\\"  # noqa: E501
            f" outboard}}}} = {p_upper_outboard:.3f}$ MW",
            fontsize=text_fontsize,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=make_bbox_props(p_upper_outboard),
            zorder=101,
        )

    arrow_props = {
        "arrowstyle": "->",
        "color": "red",
        "linewidth": 3 * scale_factor,
        "shrinkA": 14 * scale_factor,
        "shrinkB": 14 * scale_factor,
        "mutation_scale": 12 * scale_factor,
    }
    axis.annotate(
        "",
        xy=outboard_pos,
        xytext=centre_pos,
        arrowprops=arrow_props,
        zorder=1,
    )
    axis.annotate(
        "",
        xy=inboard_pos,
        xytext=centre_pos,
        arrowprops=arrow_props,
        zorder=1,
    )
    axis.annotate(
        "",
        xy=lower_outboard_pos,
        xytext=outboard_pos,
        arrowprops={
            **arrow_props,
            "connectionstyle": "angle3,angleA=0,angleB=-90",
        },
        zorder=102,
    )
    axis.annotate(
        "",
        xy=lower_inboard_pos,
        xytext=inboard_pos,
        arrowprops={
            **arrow_props,
            "connectionstyle": "angle3,angleA=180,angleB=-90",
        },
        zorder=102,
    )
    if is_double_null:
        axis.annotate(
            "",
            xy=upper_outboard_pos,
            xytext=outboard_pos,
            arrowprops={
                **arrow_props,
                "connectionstyle": "angle3,angleA=0,angleB=90",
            },
            zorder=102,
        )
        axis.annotate(
            "",
            xy=upper_inboard_pos,
            xytext=inboard_pos,
            arrowprops={
                **arrow_props,
                "connectionstyle": "angle3,angleA=180,angleB=90",
            },
            zorder=102,
        )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.get_xaxis().set_ticks([])
    axis.get_yaxis().set_ticks([])


def plot_cover_page(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    fig,
    radial_build: RadialBuild,
    colour_scheme: Literal[1, 2],
):
    """Plots a cover page for the PROCESS run, including run title, date, user, and
    summary info.

    Parameters
    ----------
    axis : plt.Axes
        The matplotlib axis object to plot on.
    mfile : MFile
        The MFILE data object containing run info.
    scan : int
        The scan number to use for extracting data.
    fig : plt.Figure
        The matplotlib figure object for additional annotations.
    radial_build:

    colour_scheme:

    """
    axis.axis("off")
    title = mfile.get("runtitle", scan=-1)
    date = mfile.get("date", scan=-1)
    time = mfile.get("time", scan=-1)
    user = mfile.get("username", scan=-1)
    procver = mfile.get("procver", scan=-1)
    tagno = mfile.get("tagno", scan=-1)
    branch_name = mfile.get("branch_name", scan=-1)
    fileprefix = mfile.get("fileprefix", scan=-1)
    optmisation_switch = int(mfile.get("i_process_run_mode", scan=-1))
    figure_merit_switch = mfile.get("i_figure_merit", scan=-1) or "N/A"
    ifail = mfile.get("ifail", scan=-1)
    nvars = mfile.get("n_iteration_variables", scan=-1)
    # Objective_function_name
    objf_name = mfile.get("objf_name", scan=-1)
    # Square_root_of_the_sum_of_squares_of_the_constraint_residuals
    sqsumsq = mfile.get("sqsumsq", scan=-1)
    # VMCON_convergence_parameter
    convergence_parameter = mfile.get("convergence_parameter", scan=-1) or "N/A"
    # Number_of_optimising_solver_iterations
    n_solver_iterations = int(mfile.get("n_solver_iterations", scan=-1)) or "N/A"

    # Objective name with minimising/maximising
    if isinstance(figure_merit_switch, str):
        objective_text = ""
    elif figure_merit_switch >= 0:
        figure_merit_switch = int(figure_merit_switch)
        objective_text = f"  -> Minimising: {objf_name}"
    else:
        figure_merit_switch = int(figure_merit_switch)
        objective_text = f"  -> Maximising: {objf_name}"

    axis.text(
        0.1,
        0.85,
        "PROCESS Run Summary",
        fontsize=28,
        ha="left",
        va="center",
        transform=fig.transFigure,
    )

    # Box 1: Run Info
    run_info = (
        f"• Run Title: {title}\n"
        f"• Date: {date}   Time: {time}\n"
        f"• User: {user}\n"
        f"• PROCESS Version: {procver}"
    )
    axis.text(
        0.1,
        0.72,
        run_info,
        fontsize=16,
        ha="left",
        va="top",
        transform=fig.transFigure,
        bbox=box_style("#e0f7fa"),
    )

    # Box 2: File/Branch Info
    # Wrap the whole "Branch Name: ..." line if too long
    max_line_len = 60
    branch_line = textwrap.fill(f"• Branch Name: {branch_name}", max_line_len)
    fileprefix = textwrap.fill(f"File Prefix: {fileprefix}", max_line_len)

    file_info = f"• Tag Number: {tagno}\n{branch_line}\n• {fileprefix}"
    axis.text(
        0.1,
        0.57,
        file_info,
        fontsize=14,
        ha="left",
        va="top",
        transform=fig.transFigure,
        bbox=box_style("#fffde7"),
    )

    # Box 3: Run Settings
    settings_info = (
        f"• Optimisation Switch: {int(optmisation_switch)}\n"
        f"     {PROCESSRunMode(int(optmisation_switch)).description}\n"
        f"• Figure of Merit Switch (i_figure_merit): {figure_merit_switch}\n"
        f"     {objective_text}\n"
        f"• Fail Status (ifail): {int(ifail)}\n"
        f"• Number of Iteration Variables: {int(nvars)}\n"
        f"• Constraint Residuals (sqrt sum sq): {sqsumsq}\n"
        f"• Convergence Parameter: {convergence_parameter}\n"
        f"• Solver Iterations: {n_solver_iterations}\n"
        f"• Runtime: {mfile.get('process_runtime', scan=-1):.6f} seconds"
    )
    axis.text(
        0.1,
        0.46,
        settings_info,
        fontsize=14,
        ha="left",
        va="top",
        transform=fig.transFigure,
        bbox=box_style("#f3e5f5"),
    )

    axis.text(
        0.1,
        0.15,
        "For more information, see the following pages.",
        fontsize=12,
        ha="left",
        va="center",
        transform=fig.transFigure,
        color="gray",
    )

    # Add a small poloidal cross-section inset on the cover page
    inset_ax = fig.add_axes([0.55, 0.2, 0.55, 0.55], aspect="equal")
    poloidal_cross_section(
        inset_ax,
        mfile,
        scan,
        demo_ranges=False,
        radial_build=radial_build,
        colour_scheme=colour_scheme,
    )
    inset_ax.set_title("")  # Remove the plot title
    inset_ax.axis("off")


__all__ = [
    "plot_cover_page",
    "plot_header",
    "plot_separatrix_power_split",
]
