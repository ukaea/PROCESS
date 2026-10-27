"""Magnets functions for PROCESS summary plots."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.patches import Circle, Rectangle
from matplotlib.path import Path as mplPath

from process.core.io.plot.summary.common import add_panel, box_style, text_layout
from process.core.io.plot.summary.constants import (
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    rtangle,
    rtangle2,
)
from process.core.io.plot.summary.magnets.cables import plot_hts_tape_geometry
from process.data_structure.superconducting_tf_coil_variables import TFWPIntegerTurnType
from process.models.geometry.tfcoil import (
    tfcoil_geometry_d_shape,
    tfcoil_geometry_rectangular_shape,
)
from process.models.superconductors import SuperconductorModel
from process.models.tfcoil import quench
from process.models.tfcoil.base import TFCoilShapeModel, TFPlasmaCaseType

if TYPE_CHECKING:
    from process.core.io.mfile import MFile

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


def _reader(mfile: MFile, scan: int):
    """Bind the scan without caching or changing MFile lookup semantics."""

    def get(name):
        return mfile.get(name, scan=scan)

    return get


def _arc_points(mfile: MFile, scan: int):
    get = _reader(mfile, scan)
    return tuple((get(f"r_tf_arc({i})"), get(f"z_tf_arc({i})")) for i in range(1, 6))


def _reference_lines(axis, *, xs=(), ys=(), **overrides):
    style = {"color": "black", "linestyle": "--", "linewidth": 0.6, "alpha": 0.5}
    style.update(overrides)
    for x in xs:
        axis.axvline(x=x, **style)
    for y in ys:
        axis.axhline(y=y, **style)


def _wp_reference_radii(get):
    return tuple(
        get(key)
        for key in (
            "r_tf_inboard_in",
            "r_tf_wp_inboard_inner",
            "r_tf_wp_inboard_outer",
            "r_tf_wp_inboard_centre",
            "r_tf_inboard_out",
        )
    )


def _finish_wp_axis(axis, xlim, ylim):
    axis.minorticks_on()
    axis.set_xlim(*xlim)
    axis.set_ylim(*ylim)
    axis.set_title("Top-down view of inboard TF coil at midplane")
    axis.set_xlabel("Radial distance [m]")
    axis.set_ylabel("Toroidal distance [m]")
    axis.legend(loc="upper left")


def _curves(inner_radius, outer_radius, half_angle, *, points=256, straight=False):
    theta = np.linspace(half_angle, -half_angle, points, endpoint=True)
    cosine, sine = np.cos(theta), np.sin(theta)
    return (
        inner_radius * cosine,
        inner_radius * sine,
        np.full(points, outer_radius) if straight else outer_radius * cosine,
        outer_radius * sine,
    )


def _outline(axis, curves):
    inner_x, inner_y, outer_x, outer_y = curves
    for x, y in (
        (inner_x, inner_y),
        (outer_x, outer_y),
        ([inner_x[0], outer_x[0]], [inner_y[0], outer_y[0]]),
        ([inner_x[-1], outer_x[-1]], [inner_y[-1], outer_y[-1]]),
    ):
        axis.plot(x, y, color="black")


def _fill_curves(axis, curves, colour, *, alpha=1.0, fill_outer=True):
    inner_x, inner_y, outer_x, outer_y = curves
    for endpoint in (0, -1):
        axis.fill_between(
            [inner_x[endpoint], outer_x[endpoint]],
            [inner_y[endpoint], outer_y[endpoint]],
            color=colour,
            alpha=alpha,
        )
    if fill_outer:
        axis.fill_between(outer_x, outer_y, color=colour, alpha=alpha)


def _plot_tf_case(axis, inner_radius, outer_radius, half_angle, case_type):
    """Draw the shared case, preserving the original fill/masking sequence."""
    circular = case_type == TFPlasmaCaseType.CIRCULAR
    if not circular and case_type != TFPlasmaCaseType.STRAIGHT:
        raise NotImplementedError("i_tf_case_geom must be 0 or 1")
    for radius in (inner_radius, outer_radius) if circular else (inner_radius,):
        axis.add_patch(
            Circle(
                (0, 0),
                radius,
                facecolor="none",
                edgecolor="black",
                linestyle="--",
            )
        )
    curves = _curves(
        inner_radius,
        outer_radius,
        half_angle,
        straight=not circular,
    )
    _outline(axis, curves)
    _fill_curves(axis, curves, "grey", alpha=0.25, fill_outer=circular)
    axis.fill_between(curves[0], curves[1], color="white", alpha=1.0)
    axis.axhline(y=0.0, color="r", linestyle="--", linewidth=0.25)
    return float(curves[3][-1])


def TF_outboard(axis: plt.Axes, item, n_tf_coils, r3, r4, w, facecolor):
    """Plot an outboard TF coil."""
    angle = item * 2 * np.pi / n_tf_coils
    cosine, sine = np.cos(angle), np.sin(angle)
    dx, dy = w * sine, w * cosine
    vertices = [
        (r3 * cosine + dx, r3 * sine - dy),
        (r4 * cosine + dx, r4 * sine - dy),
        (r4 * cosine - dx, r4 * sine + dy),
        (r3 * cosine - dx, r3 * sine + dy),
    ]
    vertices.append(vertices[0])
    axis.add_patch(
        patches.PathPatch(
            mplPath(vertices, closed=True),
            facecolor=facecolor,
            lw=0,
        )
    )


def plot_tf_coils(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    mirror_negative_x: bool = False,
):
    """Plot the TF coils and thermal shield, optionally mirrored in x."""
    get = _reader(mfile, scan)
    (x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5) = _arc_points(mfile, scan)
    x_scale = -1 if mirror_negative_x else 1
    inboard = get("dr_tf_inboard")
    outboard = get("dr_tf_outboard")
    shield = get("dr_shld_thermal_inboard")
    gap = get("dr_tf_shld_gap")
    if y3 != 0:
        print("TF coil geometry: The value of z_tf_arc(3) is not zero, but should be.")
    if shield != get("dr_shld_thermal_outboard"):
        print(
            "dr_shld_thermal_inboard and dr_shld_thermal_outboard are different. "
            "Using dr_shld_thermal_inboard for the poloidal plot of the thermal shield."
        )
    shape = int(get("i_tf_shape")) if "i_tf_shape" in mfile.data else 1
    coil_colour = TFC_COLOUR[colour_scheme - 1] if get("i_tf_sup") != 0 else "#b87333"
    geometry = {
        "x1": x1,
        "x2": x2,
        "x4": x4,
        "x5": x5,
        "y1": y1,
        "y2": y2,
        "y4": y4,
        "y5": y5,
        "dr_tf_inboard": inboard,
    }
    for offset, colour in (
        (shield + gap, THERMAL_SHIELD_COLOUR[colour_scheme - 1]),
        (gap, "white"),
        (0.0, coil_colour),
    ):
        if shape == TFCoilShapeModel.PICTURE_FRAME:
            rects = tfcoil_geometry_rectangular_shape(
                **geometry,
                dr_tf_outboard=outboard,
                offset_in=offset,
            )
        else:
            rects, vertices = tfcoil_geometry_d_shape(
                **geometry,
                x3=x3,
                rtangle=rtangle,
                rtangle2=rtangle2,
                offset_in=offset,
            )
            for vertex in vertices:
                path = mplPath([(x_scale * x, y) for x, y in vertex], closed=True)
                axis.add_patch(patches.PathPatch(path, facecolor=colour, lw=0))
        for rectangle in rects:
            axis.add_patch(
                Rectangle(
                    (x_scale * rectangle.anchor_x, rectangle.anchor_z),
                    x_scale * rectangle.width,
                    rectangle.height,
                    facecolor=colour,
                )
            )


def _casing_text(get, *, superconducting=False):
    nose_precision = 3 if superconducting else 4
    text = (
        "$\\mathbf{Casing:}$\n\nCoil half angle:"
        f" {get('rad_tf_coil_inboard_toroidal_half'):.3f} radians\n\n"
    )
    for heading, start, end, thickness, area, precision in (
        (
            "Full Coil Case",
            "r_tf_inboard_in",
            "r_tf_inboard_out",
            "dr_tf_inboard",
            "a_tf_coil_inboard_case",
            3,
        ),
        (
            "Nose Case",
            "r_tf_inboard_in",
            "r_tf_wp_inboard_inner",
            "dr_tf_nose_case",
            "a_tf_coil_nose_case",
            nose_precision,
        ),
        (
            "Plasma Case",
            "r_tf_wp_inboard_outer",
            "r_tf_inboard_out",
            "dr_tf_plasma_case",
            "a_tf_plasma_case",
            3,
        ),
    ):
        area_label = "Area of casing around WP" if heading == "Full Coil Case" else "$A$"
        text += (
            f"$\\text{{{heading}:}}$\n$r_{{start}} \\rightarrow r_{{end}}$: "
            f"{get(start):.3f} $\\rightarrow$ {get(end):.3f} m\n"
            f"$\\Delta r$: {get(thickness):.3f} m\n"
            f"{area_label}: {get(area):.{precision}f} $\\mathrm{{m}}^2$"
        )
        if heading != "Plasma Case":
            text += "\n\n"
    if superconducting:
        text += (
            "\n\n$\\text{Side Case:}$\nMinimum $\\Delta r$:"
            f" {get('dx_tf_side_case_min'):.3f} m\nAverage $\\Delta r$:"
            f" {get('dx_tf_side_case_average'):.3f} m\nMax $\\Delta r$:"
            f" {get('dx_tf_side_case_peak'):.3f} m"
        )
    return text


def _wp_header(get):
    return (
        "$\\mathbf{Winding Pack:}$\n\n$N_{\\text{turns}}$:"
        f" {int(get('n_tf_coil_turns'))} turns\n"
        "$r_{start} \\rightarrow r_{end}$:"
        f" {get('r_tf_wp_inboard_inner'):.3f} $\\rightarrow$"
        f" {get('r_tf_wp_inboard_outer'):.3f} m\n$\\Delta r$:"
        f" {get('dr_tf_wp_with_insulation'):.3f} m\n"
    )


def _superconducting_wp_panels(axis, get):
    casing = _casing_text(get, superconducting=True)
    insulation = (
        "$\\mathbf{Ground \\ Insulation:}$\n\nArea of insulation around WP:"
        f" {get('a_tf_wp_ground_insulation'):.3f} $\\mathrm{{m}}^2$\n$\\Delta r$:"
        f" {get('dx_tf_wp_insulation'):.4f} m\n\nWP Insertion Gap:\n$\\Delta r$:"
        f" {get('dx_tf_wp_insertion_gap'):.4f} m"
    )
    wp = _wp_header(get) + "\n"
    for label, key in (
        ("with insulation", "a_tf_wp_with_insulation"),
        ("no insulation", "a_tf_wp_no_insulation"),
        ("total turn insulation", "a_tf_coil_wp_turn_insulation"),
        ("total turn steel", "a_tf_wp_steel"),
        ("total conductor", "a_tf_wp_conductor"),
        ("total non-cooling void", "a_tf_wp_extra_void"),
    ):
        wp += f"$A$, {label}: {get(key):.4f} $\\mathrm{{m}}^2$\n"
    wp += (
        f"\nPrimary WP:\n$\\Delta x$: {get('dx_tf_wp_primary_toroidal'):.4f} m\n"
        f"\nSecondary WP:\n$\\Delta x$: {get('dx_tf_wp_secondary_toroidal'):.4f} m\n"
        f"\n$J$ no insulation: {get('j_tf_wp') / 1e6:.4f} MA/m$^2$"
    )
    general = (
        "$\\mathbf{General \\ info:}$\n\n$N_{\\text{TF,coil}}$:"
        f" {get('n_tf_coils')}\nSelf inductance of single coil:"
        f" {get('ind_tf_coil') * 1e6:.4f} $\\mu$H\n"
        "Stored energy of all coils:"
        f" {get('e_tf_magnetic_stored_total_gj'):.4f} GJ\n"
        "Stored energy of a single coil:"
        f" {get('e_tf_coil_magnetic_stored') / 1e9:.2f} GJ\n"
        "Total area of steel in coil:"
        f" {get('a_tf_coil_inboard_steel'):.4f} $\\mathrm{{m}}^2$\n"
        "Total area fraction of steel:"
        f" {get('f_a_tf_coil_inboard_steel'):.4f}\nTotal area fraction of insulation:"
        f" {get('f_a_tf_coil_inboard_insulation'):.4f}\n$A$, all insulation in coil:"
        f" {get('a_tf_coil_inboard_insulation'):.4f} $\\mathrm{{m}}^2$\n"
    )
    for x, y, width, height, text, colour, foreground in (
        (0.5, 0.95, 0.215, 0.38, casing, "grey", "black"),
        (0.5, 0.55, 0.215, 0.12, insulation, "green", "black"),
        (0.725, 0.95, 0.2, 0.33, wp, "blue", "white"),
        (0.725, 0.55, 0.22, 0.18, general, "wheat", "black"),
    ):
        add_panel(
            axis,
            x,
            y,
            width,
            height,
            text,
            box_style(colour),
            text_kwargs={"color": foreground},
        )


def _plot_superconducting_pack(axis, get):
    inner = get("r_tf_wp_inboard_inner")
    outer = get("r_tf_wp_inboard_outer")
    centre = get("r_tf_wp_inboard_centre")
    radial = get("dr_tf_wp_with_insulation")
    primary = get("dx_tf_wp_primary_toroidal")
    secondary = get("dx_tf_wp_secondary_toroidal")
    insulation = get("dx_tf_wp_insulation")
    inset = insulation + get("dx_tf_wp_insertion_gap")
    geometry = round(get("i_tf_wp_geom"))

    if geometry == 0:
        if TFWPIntegerTurnType(get("i_tf_turns_integer")) == TFWPIntegerTurnType.INTEGER:
            long_turns = round(get("n_tf_wp_layers"))
            short_turns = round(get("n_tf_wp_pancakes"))
        else:
            ratio = (radial - 2 * inset) / (primary - 2 * inset)
            root_turns = round(np.sqrt(round(get("n_tf_coil_turns")) / ratio), 1)
            long_turns, short_turns = round(root_turns * ratio), round(root_turns)
        rectangles = (
            (inner, -primary / 2, radial, primary, "darkgreen"),
            (
                inner + inset,
                -primary / 2 + inset,
                radial - 2 * inset,
                primary - 2 * inset,
                "blue",
            ),
        )
    elif geometry == 1:
        # Keep the asymmetric inner/outer widths from the original geometry.
        rectangles = (
            (inner, -secondary / 2, radial / 2 + insulation, secondary, "darkgreen"),
            (centre, -primary / 2, radial / 2, primary, "darkgreen"),
            (
                centre + inset,
                -primary / 2 + inset,
                radial / 2 - 2 * inset,
                primary - 2 * inset,
                "blue",
            ),
            (
                inner + inset,
                -secondary / 2 + inset,
                radial / 2,
                secondary - 2 * inset,
                "blue",
            ),
        )
    elif geometry == 2:
        for offset, colour in ((0.0, "darkgreen"), (inset, "blue")):
            axis.add_patch(
                patches.Polygon(
                    [
                        (inner + offset, -secondary / 2 + offset),
                        (inner + offset, secondary / 2 - offset),
                        (outer - offset, primary / 2 - offset),
                        (outer - offset, -primary / 2 + offset),
                    ],
                    color=colour,
                )
            )
        return
    else:
        return

    for x, y, width, height, colour in rectangles:
        axis.add_patch(Rectangle((x, y), width, height, color=colour))
    if geometry == 0:
        # Preserve the original spacing formula, including its single-inset span.
        style = {"color": "white", "linewidth": 0.25, "linestyle": "dashed"}
        for i in range(1, long_turns):
            x = inner + inset + i * (radial - inset) / long_turns
            axis.plot([x, x], [-primary / 2 + inset, primary / 2 - inset], **style)
        for i in range(1, short_turns):
            y = -primary / 2 + inset + i * (primary - inset) / short_turns
            axis.plot([inner + inset, outer - inset], [y, y], **style)


def plot_superconducting_tf_wp(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot the superconducting inboard TF case and winding pack."""
    get = _reader(mfile, scan)
    if round(get("i_tf_sup")) != 1:
        return
    inner, outer = get("r_tf_inboard_in"), get("r_tf_inboard_out")
    lower_y = _plot_tf_case(
        axis,
        inner,
        outer,
        get("rad_tf_coil_inboard_toroidal_half"),
        get("i_tf_case_geom"),
    )
    _plot_superconducting_pack(axis, get)
    peak_radius = get("r_b_tf_inboard_peak")
    axis.plot(
        peak_radius,
        0,
        marker="o",
        color="red",
        label=(
            f"Peak axisymmetric field: {get('b_tf_inboard_peak_symmetric'):.3f} T\n"
            "Peak non-axisymmetric field with ripple: "
            f"{get('b_tf_inboard_peak_with_ripple'):.3f} T\n"
            "$\\frac{B_{\\text{axisymmetric}}}{B_{\\text{non-axisymmetric}}}$: "
            f"{get('f_b_tf_inboard_peak_ripple_symmetric'):.3f}\n"
            f"$r_{{\\text{{peak}}}}$={peak_radius:.3f} m"
        ),
    )
    half_widths = (
        get("dx_tf_wp_secondary_toroidal") / 2,
        get("dx_tf_wp_primary_toroidal") / 2,
        get("dx_tf_wp_primary_toroidal") / 2 + get("dx_tf_side_case_peak"),
    )
    _reference_lines(
        axis,
        xs=_wp_reference_radii(get),
        ys=tuple(sign * width for width in half_widths for sign in (1, -1)),
    )
    _superconducting_wp_panels(axis, get)
    _finish_wp_axis(axis, (inner * 0.8, outer * 1.1), (lower_y * 1.25, -lower_y * 1.25))


def plot_resistive_tf_wp(axis: plt.Axes, mfile: MFile, scan: int, fig):
    """Plot the resistive inboard TF case and winding pack."""
    get = _reader(mfile, scan)
    inner, outer = get("r_tf_wp_inboard_inner"), get("r_tf_wp_inboard_outer")
    half_angle = get("rad_tf_coil_inboard_toroidal_half")
    case_type = get("i_tf_case_geom")
    _plot_tf_case(
        axis,
        get("r_tf_inboard_in"),
        get("r_tf_inboard_out"),
        half_angle,
        case_type,
    )
    insulation_curves = _curves(inner, outer, half_angle, points=500)
    _outline(axis, insulation_curves)
    if case_type == TFPlasmaCaseType.CIRCULAR:
        _fill_curves(axis, insulation_curves, "green")

    insulation = get("dx_tf_wp_insulation")
    angle_offset = max(
        insulation / radius if radius > 0 else 0 for radius in (inner, outer)
    )
    wp_half_angle = half_angle - angle_offset
    wp_inner, wp_outer = inner + insulation, outer - insulation
    wp_curves = _curves(wp_inner, wp_outer, wp_half_angle)
    _outline(axis, wp_curves)
    colour = "#b0c4de" if get("i_tf_sup") == 2 else "#b87333"
    _fill_curves(axis, wp_curves, colour)

    n_turns = int(get("n_tf_coil_turns"))
    if n_turns > 0:
        for angle in np.linspace(wp_half_angle, -wp_half_angle, n_turns + 1)[1:-1]:
            axis.plot(
                [wp_inner * np.cos(angle), wp_outer * np.cos(angle)],
                [wp_inner * np.sin(angle), wp_outer * np.sin(angle)],
                color="white",
                linewidth=0.5,
                linestyle="--",
            )
    peak_radius = get("r_b_tf_inboard_peak")
    axis.plot(
        peak_radius,
        0,
        marker="o",
        color="red",
        label=f"Peak Field: {get('b_tf_inboard_peak_symmetric'):.2f} T\n"
        f"r={peak_radius:.3f} m",
    )
    _reference_lines(axis, xs=_wp_reference_radii(get))
    lower_y = wp_curves[3][-1]
    _finish_wp_axis(
        axis,
        (0.0, get("r_tf_inboard_out") * 1.1),
        (lower_y * 1.65, -lower_y * 1.65),
    )
    axis.text(
        0.05,
        0.975,
        "*Turn insulation and cooling pipes not shown",
        fontsize=9,
        va="top",
        ha="left",
        color="black",
        transform=fig.transFigure,
    )


def plot_resistive_tf_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot information about the resistive TF coils."""
    get = _reader(mfile, scan)
    insulation = (
        "$\\mathbf{Insulation:}$\n\nArea of insulation around WP:"
        f" {get('a_tf_wp_ground_insulation'):.3f} $\\mathrm{{m}}^2$\n$\\Delta r$:"
        f" {get('dx_tf_wp_insulation'):.4f} m\n\n$\\text{{Turn Insulation:}}$\n"
        "$\\Delta r$:"
        f" {get('dx_tf_turn_insulation'):.4f} m"
    )
    wp = _wp_header(get) + (
        f"$A$, with insulation: {get('a_tf_wp_with_insulation'):.3f} $\\mathrm{{m}}^2$\n"
        f"$A$, no insulation: {get('a_tf_wp_no_insulation'):.3f} $\\mathrm{{m}}^2$\n\n"
        f"Current per turn: {get('c_tf_turn') / 1e3:.3f} $\\mathrm{{kA}}$\n"
        "Resistive conductor per coil: "
        f"{get('a_res_tf_coil_conductor'):.3f} $\\mathrm{{m}}^2$\n"
        f"Coolant area void fraction per turn: {get('fcoolcp'):.3f}"
    )
    general = (
        "$\\mathbf{General \\ info:}$\n\nSelf inductance:"
        f" {get('ind_tf_coil') * 1e6:.4f} $\\mu$H\nStored energy of all coils:"
        f" {get('e_tf_magnetic_stored_total_gj'):.4f} GJ\n"
    )
    cooling = "$\\mathbf{Cooling \\ info:}$\n\n"
    for label, key, fmt, unit in (
        ("Coolant inlet temperature", "temp_cp_coolant_inlet", ".2f", "K"),
        ("Coolant temperature rise", "dtemp_cp_coolant", ".2f", "K"),
        ("Coolant velocity", "vel_cp_coolant_midplane", ".2f", "$\\mathrm{ms^{-1}}$\n"),
        ("Average CP temperature", "temp_cp_average", ".2f", "K"),
        ("CP resistivity", "rho_cp", ".2e", "$\\Omega \\mathrm{m}$"),
        ("Leg resistivity", "rho_tf_leg", ".2e", "$\\Omega \\mathrm{m}$"),
        ("Leg resistance", "res_tf_leg", ".2e", "$\\Omega$"),
        ("CP resistive losses", "p_cp_resistive", ",.2f", "$\\mathrm{W}$"),
        ("Leg resistive losses", "p_tf_leg_resistive", ",.2f", "$\\mathrm{W}$"),
        ("Joints resistive losses", "p_tf_joints_resistive", ",.2f", "$\\mathrm{W}$"),
    ):
        cooling += f"{label}: {get(key):{fmt}} {unit}\n"
    for x, y, width, height, text, colour in (
        (0.775, 0.925, 0.19, 0.31, _casing_text(get), "grey"),
        (0.775, 0.62, 0.19, 0.09, insulation, "green"),
        (0.77, 0.475, 0.19, 0.3, wp, "blue"),
        (0.55, 0.475, 0.19, 0.1, general, "wheat"),
        (0.55, 0.35, 0.19, 0.22, cooling, "wheat"),
    ):
        add_panel(axis, x, y, width, height, text, box_style(colour))


@dataclass(frozen=True)
class _TurnGeometry:
    integer: bool
    width: float
    height: float
    cable_radial: float
    cable_toroidal: float
    steel: float
    insulation: float
    corner_radius: float
    pipe_diameter: float

    @classmethod
    def from_mfile(cls, mfile, scan):
        get = _reader(mfile, scan)
        integer = (
            TFWPIntegerTurnType(get("i_tf_turns_integer")) == TFWPIntegerTurnType.INTEGER
        )
        if integer:
            width, height = get("dr_tf_turn"), get("dx_tf_turn")
            radial, toroidal = (
                get("dr_tf_turn_cable_space"),
                get("dx_tf_turn_cable_space"),
            )
        else:
            width = height = get("dx_tf_turn_general")
            radial = toroidal = get("dx_tf_turn_cable_space_average")
        return cls(
            integer,
            width,
            height,
            radial,
            toroidal,
            get("dx_tf_turn_steel"),
            get("dx_tf_turn_insulation"),
            get("radius_tf_turn_cable_space_corners"),
            get("dia_tf_turn_coolant_channel"),
        )

    @property
    def center(self):
        return self.width / 2, self.height / 2

    @property
    def cable_bounds(self):
        inset = self.insulation + self.steel
        return inset, inset, self.width - 2 * inset, self.height - 2 * inset


def _plot_turn_envelope(axis, turn):
    for inset, colour in ((0.0, "red"), (turn.insulation, "grey")):
        axis.add_patch(
            Rectangle(
                (inset, inset),
                turn.width - 2 * inset,
                turn.height - 2 * inset,
                facecolor=colour,
                edgecolor="black",
            )
        )


def _finish_turn_axis(axis, turn, *, margin=0.05, equal=False):
    axis.minorticks_on()
    axis.set_title("WP Turn Structure")
    axis.set_xlim(-turn.width * margin, turn.width * (1 + margin))
    axis.set_ylim(-turn.height * margin, turn.height * (1 + margin))
    if equal:
        axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("r [m]")
    axis.set_ylabel("x [m]")


def _superconductor_text(get, *, croco=False):
    time_format = ".4e" if croco else ".4f"
    separator = " " if croco else "\n"
    text = (
        "$\\mathbf{Superconductor:}$\n\nSuperconductor used:"
        f"{separator}{SuperconductorModel(get('i_tf_sc_mat')).full_name}\n"
    )
    for label, key, fmt, unit in (
        (
            "Critical field at zero\ntemperature and strain",
            "b_tf_superconductor_critical_zero_temp_strain",
            ".4f",
            "T",
        ),
        (
            "Critical temperature at\nzero field and strain",
            "temp_tf_superconductor_critical_zero_field_strain",
            ".4f",
            "K",
        ),
        ("Temperature at conductor", "tftmp", ".4f", "K"),
        ("Field at conductor", "b_tf_inboard_peak_with_ripple", ".4f", "T"),
        (
            "Superconductor critical current density at\noperating conditions",
            "j_tf_superconductor_critical",
            ".2e",
            "A/m$^2$",
        ),
        ("$I_{\\text{TF,turn critical}}$", "c_turn_cables_critical", ",.2f", "A"),
        ("$I_{\\text{TF,turn}}$", "c_tf_turn", ",.2f", "A"),
        ("Critical current ratio", "f_c_tf_turn_operating_critical", ",.4f", ""),
        (
            "Superconductor temperature\nmargin",
            "temp_tf_superconductor_margin",
            ",.4f",
            "K",
        ),
    ):
        text += f"{label}: {get(key):{fmt}}{(' ' + unit) if unit else ''}\n"
    text += "\n$\\mathbf{Quench:}$\n\n"
    for label, key, fmt, unit in (
        ("Quench dump time", "t_tf_superconductor_quench", time_format, "s"),
        ("Quench detection time", "t_tf_quench_detection", time_format, "s"),
        (
            "User input max temperature\nduring quench",
            "temp_tf_conductor_quench_max",
            ".2f",
            "K",
        ),
        (
            "Required maximum WP current\ndensity for heat protection",
            "j_tf_wp_quench_heat_max",
            ".2e",
            "A/m$^2$",
        ),
    ):
        text += f"{label}: {get(key):{fmt}} {unit}\n"
    return text


def _turn_panel_texts(get, turn, *, croco=False):
    cable_dimensions = (
        f"Cable space:\n$\\Delta r$: {turn.cable_radial:.3e} m\n"
        f"$\\Delta x$: {turn.cable_toroidal:.3e} m\n"
        if turn.integer
        else f"$\\Delta r:$ {turn.cable_radial:.3e} m\n"
    )
    return {
        "insulation": (
            f"$\\mathbf{{Turn \\ Insulation:}}$\n\n$\\Delta r:${turn.insulation:.3e} m"
        ),
        "steel": (
            "$\\mathbf{Steel \\ Conduit:}$\n\n"
            f"$\\Delta r:${turn.steel:.3e} m\n$A$: {get('a_tf_turn_steel'):.3e} m$^2$"
        ),
        "cable": (
            "$\\mathbf{Cable \\ Space:}$\n\n"
            + cable_dimensions
            + f"Corner radius, $r$: {turn.corner_radius:.3e} m\n"
            + "Cable area with no cooling channel or gaps: "
            + f"{get('a_tf_turn_cable_space_no_void'):.3e} m$^2$\n"
            + "Extra cable space area void fraction: "
            + f"{get('f_a_tf_turn_cable_space_extra_void')}\nTrue cable space area: "
            + f"{get('a_tf_turn_cable_space_effective'):.3e} m$^2$"
        ),
        "turn": (
            "$\\mathbf{Turn:}$\n\n"
            f"$\\Delta r$: {turn.width:.3e} m\n$\\Delta x$: {turn.height:.3e} m"
        ),
        "cooling": (
            "$\\mathbf{Cooling:}$\n\n"
            f"$\\varnothing$: {turn.pipe_diameter:.3e} m\n"
            "Total area of all coolant channels: "
            f"{get('a_tf_wp_coolant_channels'):.4f} m$^2$"
        ),
        "superconductor": _superconductor_text(get, croco=croco),
    }


def _plot_turn_panels(axis, get, turn, *, fig=None):
    croco = fig is not None
    texts = _turn_panel_texts(get, turn, croco=croco)
    for key, x, y, width, height, colour in (
        ("insulation", 0.39, 0.90, 0.115, 0.075, "red"),
        ("cable", 0.39, 0.70 if croco else 0.67, 0.31, 0.195, "royalblue"),
        ("cooling", 0.39, 0.80 if croco else 0.79, 0.25, 0.070, "white"),
        ("steel", 0.632, 0.90, 0.115, 0.075, "grey"),
        ("turn", 0.51, 0.90, 0.115, 0.075, "wheat"),
        ("superconductor", 0.76, 0.90, 0.235, 0.47, "#6dd3f7"),
    ):
        add_panel(axis, x, y, width, height, texts[key], box_style(colour))


def _outside_rounded_corner(cx, cy, x, y, width, height, corner_radius, radius):
    """Apply the original four corner checks without repeating each expression."""
    for corner_x, corner_y, sx, sy in (
        (x + corner_radius, y + corner_radius, -1, -1),
        (x + width - corner_radius, y + corner_radius, 1, -1),
        (x + width - corner_radius, y + height - corner_radius, 1, 1),
        (x + corner_radius, y + height - corner_radius, -1, 1),
    ):
        if (
            sx * (cx - corner_x) > 0
            and sy * (cy - corner_y) > 0
            and np.hypot(cx - corner_x, cy - corner_y) > corner_radius - radius
        ):
            return True
    return False


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
    """Pack strands using the original grid, spacing and collision criteria."""
    x, y, width, height = cable_space_bounds
    radius = strand_diameter / 2
    ideal_fraction = np.pi / (2 * np.sqrt(3))
    spacing = strand_diameter * np.sqrt(ideal_fraction / (1 - void_fraction))
    n_rows = int((height - 2 * radius) // (spacing * np.sqrt(3) / 2))
    n_cols = int((width - 2 * radius) // spacing)
    sc_radius = np.sqrt(1 - f_a_tf_turn_cable_copper) * radius
    placed = []
    for row in range(n_rows):
        cy = (y + radius + row * spacing * np.sqrt(3) / 2) * 1.07
        offset = spacing / 2 if row % 2 else 0
        for col in range(n_cols):
            cx = (x + radius + col * spacing + offset) * 1.05
            if cx > x + width - radius or cy > y + height - radius:
                continue
            if np.hypot(cx - pipe_center[0], cy - pipe_center[1]) < pipe_radius + radius:
                continue
            if corner_radius > 0 and _outside_rounded_corner(
                cx,
                cy,
                x,
                y,
                width,
                height,
                corner_radius,
                radius,
            ):
                continue
            if not any(
                np.hypot(cx - px, cy - py) < strand_diameter for px, py in placed
            ):
                placed.append((cx, cy))
                axis.add_patch(
                    Circle(
                        (cx, cy),
                        radius,
                        facecolor="#b87333",
                        edgecolor="#8B4000",
                        linewidth=0.1,
                        alpha=0.8,
                    )
                )
                axis.add_patch(
                    Circle(
                        (cx, cy),
                        sc_radius,
                        facecolor="black",
                        linewidth=0.3,
                        alpha=0.5,
                    )
                )
            if len(placed) >= n_strands:
                break
        if len(placed) >= n_strands:
            break
    return len(placed), n_rows * n_cols


def plot_tf_cable_in_conduit_turn(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot the CICC turn envelope, coolant channel and packed strands."""
    get = _reader(mfile, scan)
    turn = _TurnGeometry.from_mfile(mfile, scan)
    _plot_turn_envelope(axis, turn)
    x, y, width, height = turn.cable_bounds
    axis.add_patch(
        Rectangle(
            (x, y),
            width,
            height,
            facecolor="none",
            edgecolor="black",
            linestyle="--",
            linewidth=1.0 if turn.integer else 1.2,
            alpha=0.5,
        )
    )

    strand_diameter = get("dia_tf_turn_superconducting_cable")
    void_fraction = get("f_a_tf_turn_cable_space_extra_void")
    if void_fraction > 0.001:
        _pack_strands_rectangular_with_obstacles(
            cable_space_bounds=turn.cable_bounds,
            pipe_center=turn.center,
            pipe_radius=turn.pipe_diameter / 2,
            strand_diameter=strand_diameter,
            void_fraction=void_fraction,
            axis=axis,
            corner_radius=turn.corner_radius,
            n_strands=get("n_tf_turn_superconducting_cables"),
            f_a_tf_turn_cable_copper=get("f_a_tf_turn_cable_copper"),
        )
    axis.add_patch(
        patches.FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle=patches.BoxStyle("Round", pad=0, rounding_size=turn.corner_radius),
            facecolor="royalblue",
            edgecolor="black",
            alpha=0.4,
        )
    )
    axis.add_patch(
        Circle(
            turn.center,
            turn.pipe_diameter / 2,
            facecolor="white",
            edgecolor="black",
        )
    )
    _finish_turn_axis(axis, turn)
    _plot_turn_panels(axis, get, turn)


def plot_tf_croco_turn(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plot a CroCo turn; retain the original non-integer-only cable layout."""
    get = _reader(mfile, scan)
    turn = _TurnGeometry.from_mfile(mfile, scan)
    if not turn.integer:
        _plot_turn_envelope(axis, turn)
        diameter = get("dia_tf_turn_croco_cable")
        for radius, colour in ((1.5 * diameter, "white"), (diameter / 2, "#B87333")):
            axis.add_patch(
                Circle(
                    turn.center,
                    radius,
                    facecolor=colour,
                    edgecolor="black",
                    linewidth=1.2,
                )
            )
        cable_arguments = {
            "dia_croco_strand": diameter,
            "dx_croco_strand_copper": get("dx_tf_croco_strand_copper"),
            "dr_hts_tape": get("dr_tf_hts_tape"),
            "dx_croco_strand_tape_stack": get("dx_tf_croco_strand_tape_stack"),
            "n_croco_strand_hts_tapes": get("n_tf_croco_strand_hts_tapes"),
            "dx_hts_tape_rebco": get("dx_tf_hts_tape_rebco"),
            "dx_hts_tape_copper": get("dx_tf_hts_tape_copper"),
            "dx_hts_tape_hastelloy": get("dx_tf_hts_tape_hastelloy"),
        }
        for angle in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            plot_corc_cable_geometry(
                axis=axis,
                r_centre=turn.center[0] + diameter * np.cos(angle),
                z_centre=turn.center[1] + diameter * np.sin(angle),
                **cable_arguments,
                show_legend=False,
            )
    # The original CroCo plot uses square limits even for integer turns.
    _finish_turn_axis(axis, turn, margin=0.025, equal=True)
    axis.set_ylim(-turn.width * 0.025, turn.width * 1.025)
    _plot_turn_panels(axis, get, turn, fig=fig)


def _dimension(
    axis,
    start,
    end,
    label_position,
    label,
    *,
    rotation=0,
    arrowstyle="<->",
    text_zorder=None,
    arrow_zorder=None,
):
    arrow_options = {} if arrow_zorder is None else {"zorder": arrow_zorder}
    text_options = {} if text_zorder is None else {"zorder": text_zorder}
    axis.annotate(
        "",
        xy=start,
        xytext=end,
        arrowprops={"arrowstyle": arrowstyle, "color": "black"},
        **arrow_options,
    )
    axis.text(
        *label_position,
        label,
        fontsize=7,
        color="black",
        rotation=rotation,
        va="center",
        ha="center",
        bbox=box_style("pink", linewidth=0),
        **text_options,
    )


def plot_tf_coil_structure(axis: plt.Axes, mfile: MFile, scan: int, colour_scheme=1):
    """Plot the TF coil poloidal cross-section and its dimensions."""
    plot_tf_coils(axis, mfile, scan, colour_scheme)
    get = _reader(mfile, scan)
    points = _arc_points(mfile, scan)
    (x1, _), (x2, y2), _, (x4, y4), _ = points
    labels = []
    for i, (x, y) in enumerate(points, 1):
        axis.plot(x, y, "ko", markersize=8)
        axis.text(
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
    shape = get("i_tf_shape")
    thickness = get("dr_tf_inboard")
    out_thickness = get("dr_tf_outboard")
    inner = get("r_tf_inboard_in")
    inner_out = get("r_tf_inboard_out")
    outer_in = get("r_tf_outboard_in")
    outer_out = outer_in + out_thickness
    inside_half = get("z_tf_inside_half")
    top, bottom = y2 + thickness, y4 - thickness
    if shape == 1:
        _dimension(
            axis,
            (x2, y2),
            (x4, y4),
            (x2, 0),
            f"{y2 - y4:.3f} m",
            rotation=270,
            text_zorder=100,
        )
        height_x = x2 * 0.9
    elif shape == 2:
        height_x = (x2 - x1) / 2
        # Preserve the original picture-frame-only full-height arrow.
        axis.annotate(
            "",
            xy=(height_x, top),
            xytext=(height_x, bottom),
            arrowprops={"arrowstyle": "<->", "color": "black"},
        )
    else:
        raise NotImplementedError("i_tf_shape must be 1 or 2")
    axis.text(
        height_x,
        0,
        f"{y2 + 2 * thickness - y4:.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        va="center",
        ha="center",
        bbox=box_style("pink", linewidth=0),
        zorder=101,
    )
    for x, endpoint, label_y, value in (
        (-2.0, top, top / 2, top),
        (-2.0, bottom, -get("z_tf_top") / 2, bottom),
        (-1.0, y2, y2 / 2, y2),
        (-1.0, -inside_half, -inside_half / 2, inside_half),
    ):
        _dimension(
            axis, (x, 0), (x, endpoint), (x, label_y), f"{value:.3f} m", rotation=270
        )
        axis.axhline(y=endpoint, color="black", linestyle="--", linewidth=1)
    label_x = (inner_out + outer_in) / 1.5
    _dimension(
        axis,
        (inner_out, -inside_half / 12),
        (outer_in, -inside_half / 12),
        (label_x, -inside_half / 12),
        f"{get('dr_tf_internal_midplane'):.3f} m",
        text_zorder=100,
    )
    _dimension(
        axis,
        (inner, 0),
        (outer_out, 0),
        (label_x, 0),
        f"{get('dr_tf_full_midplane'):.3f} m",
        arrowstyle="<|-|>",
        text_zorder=100,
        arrow_zorder=100,
    )
    for x, label in (
        (inner, "TF Inboard Start"),
        (inner_out, "TF Inboard End"),
        (outer_in, "TF Outboard Start"),
        (outer_out, "TF Outboard End"),
    ):
        axis.axvline(
            x, color="black", linestyle="--", linewidth=1, alpha=0.5, label=label
        )
    text_style = {"fontsize": 7, "color": "black", "va": "center"}
    for x, value in ((inner, thickness), (outer_in, out_thickness)):
        axis.text(
            x,
            bottom * 1.1,
            rf"$\Delta r = ${value:.3f} m",
            **text_style,
            bbox=box_style("pink", linewidth=0),
        )
    info_x = outer_in + 2 * out_thickness
    for y, text, zorder in (
        (0.0, f"Length of coil = {get('len_tf_coil'):.3f} m", 100),
        (
            -1.0,
            (
                "$\\Delta Z$ upper and lower to midplane ="
                f" {get('dz_tf_upper_lower_midplane'):.3f} m"
            ),
            None,
        ),
        (
            -2.0,
            (
                "$\\Delta Z$ coil centre to plasma centre ="
                f" {get('dz_tf_plasma_centre_offset'):.3f} m"
            ),
            None,
        ),
    ):
        options = {} if zorder is None else {"zorder": zorder}
        axis.text(
            info_x, y, text, **text_style, bbox=box_style("pink", linewidth=0), **options
        )
    _dimension(
        axis,
        (inner, 0),
        (0, 0),
        (inner / 2, 0),
        f"{inner:.3f} m",
        arrowstyle="->",
        text_zorder=101,
    )
    if shape == 1:
        _dimension(
            axis,
            (outer_out, top),
            (0, top),
            (inner / 2, top),
            f"{outer_out:.3f} m",
            arrowstyle="->",
        )
        axis.plot(0, top, marker="o", color="black", markersize=7, zorder=100)
    axis.axhline(y=(y2 + y4) / 2, color="red", linestyle="--", linewidth=1.0, zorder=5)
    axis.plot(0, 0, marker="o", color="red", markersize=7)
    axis.axvline(0, color="red", linestyle="--", linewidth=1)
    axis.axhline(y=0, color="red", linestyle="--", linewidth=1)
    axis.set_xlim(-3.0, outer_out * 1.4)
    axis.set_ylim(bottom * 1.2, top * 1.2)
    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("TF Coil Poloidal Cross-Section")
    axis.minorticks_on()
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    axis.legend(labels, loc="upper center", bbox_to_anchor=(1.01, 0.85), ncol=1)


def _plot_stress_summary(axis, radius, series, title):
    """Draw all background curves before foreground curves and boundary markers."""
    for layer in range(radius.shape[0]):
        for values, background, _, _, _ in series:
            axis.plot(radius[layer], values[layer], "-", linewidth=3.5, color=background)
    for layer in range(radius.shape[0]):
        for values, _, colour, linestyle, label in series:
            axis.plot(
                radius[layer],
                values[layer],
                linestyle,
                color=colour,
                label=label if layer == 0 else None,
            )
    boundary_radius = radius[:, [0, -1]].ravel()
    for values, _, colour, _, _ in series:
        axis.plot(
            boundary_radius, values[:, [0, -1]].ravel(), "|", markersize=10, color=colour
        )
    axis.grid(True)
    axis.set_ylabel(r"$\sigma$ [$MPa$]", fontsize=12)
    axis.set_title(title)
    axis.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=10)


def plot_tf_stress(axis, mfile: MFile):
    """Plot layered stresses and displacement from SIG_TF.json.

    Incomplete layers are rejected rather than silently truncated. Unused strain
    and duplicate boundary arrays from the original implementation are omitted.

    Raises
    ------
    ValueError
        layers not 1 or the length of the number of points
    """
    filename = mfile.filename.with_name(
        mfile.filename.name.replace("MFILE.DAT", "SIG_TF.json")
    )
    with filename.open() as stream:
        data = json.load(stream)
    points_per_layer = int(data["Points per layers"])
    n_points = len(data["Radius (m)"])
    if points_per_layer <= 0 or n_points == 0:
        raise ValueError("Stress data must contain non-empty layers.")
    if n_points % points_per_layer:
        raise ValueError("Stress data contains an incomplete layer.")

    def layered(key, *, allow_scalar=False):
        values = np.asarray(data[key], dtype=float)
        if allow_scalar and values.size == 1:
            values = np.full(n_points, values.item())
        if values.ndim != 1 or values.size != n_points:
            raise ValueError(
                f"{key!r} must contain {n_points} values"
                + (" or one scalar value." if allow_scalar else ".")
            )
        return values.reshape(-1, points_per_layer)

    radius = layered("Radius (m)")
    radial = layered("Radial smear stress (MPa)")
    toroidal = layered("Toroidal smear stress (MPa)")
    vertical = layered("Vertical smear stress (MPa)")
    tresca = np.maximum(np.abs(radial), np.abs(toroidal)) + vertical
    structure = (
        (
            layered("Radial stress (MPa)"),
            "lightblue",
            "dodgerblue",
            "--",
            r"$\sigma_{rr}$",
        ),
        (
            layered("Toroidal stress (MPa)"),
            "wheat",
            "orange",
            "--",
            r"$\sigma_{\theta\theta}$",
        ),
        (
            layered("Vertical stress (MPa)", allow_scalar=True),
            "lightgrey",
            "mediumseagreen",
            "--",
            r"$\sigma_{zz}$",
        ),
        (
            layered("CEA Tresca stress (MPa)"),
            "pink",
            "crimson",
            "-",
            r"$\sigma_{TRESCA}$",
        ),
        (
            layered("Von-Mises stress (MPa)"),
            "violet",
            "darkviolet",
            "-",
            r"$\sigma_{Von\ mises}$",
        ),
    )
    smeared = (
        (radial, "lightblue", "dodgerblue", "--", r"$\sigma_{rr}^\mathrm{smeared}$"),
        (toroidal, "wheat", "orange", "--", r"$\sigma_{\theta\theta}^\mathrm{smeared}$"),
        (
            vertical,
            "lightgrey",
            "mediumseagreen",
            "--",
            r"$\sigma_{zz}^\mathrm{smeared}$",
        ),
        (tresca, "pink", "crimson", "-", r"$\sigma_{TRESCA}^\mathrm{smeared}$"),
    )
    _plot_stress_summary(axis[0], radius, structure, "Structure Stress Summary")
    _plot_stress_summary(axis[1], radius, smeared, "Smeared Stress Summary")
    displacement = layered("rad. displacement (mm)")
    for layer_radius, layer_displacement in zip(radius, displacement, strict=True):
        axis[2].plot(layer_radius, layer_displacement, color="dodgerblue")
    axis[2].grid(True)
    axis[2].set_ylabel(r"$u_{r}$ [mm]", fontsize=12)
    axis[2].set_xlabel(r"$R$ [$m$]", fontsize=12)
    axis[2].set_title("Radial Displacement")
    for ax in axis:
        ax.minorticks_on()
        ax.tick_params(axis="both", labelsize=12)
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
    """Plot a CroCo strand cable, its jacket, solder and HTS tape stack."""
    for radius, edge, face, linewidth, label in (
        (dia_croco_strand / 2, "black", "#B87333", 0.5, "Copper jacket"),
        (dia_croco_strand / 2 - dx_croco_strand_copper, "grey", "grey", 2, "Solder"),
    ):
        axis.add_patch(
            Circle(
                (r_centre, z_centre),
                radius,
                edgecolor=edge,
                facecolor=face,
                linewidth=linewidth,
                label=label if show_legend else "_nolegend_",
            )
        )
    axis.add_patch(
        Rectangle(
            (r_centre - dr_hts_tape / 2, z_centre - dx_croco_strand_tape_stack / 2),
            dr_hts_tape,
            dx_croco_strand_tape_stack,
            edgecolor="black",
            facecolor=None,
            linewidth=0.1,
            alpha=0.5,
            linestyle="--",
            label="HTS Tape Stack" if show_legend else "_nolegend_",
        )
    )
    for i in range(int(n_croco_strand_hts_tapes)):
        y_start = (
            z_centre
            - dx_croco_strand_tape_stack / 2
            + i * dx_croco_strand_tape_stack / n_croco_strand_hts_tapes
        )
        plot_hts_tape_geometry(
            axis=axis,
            r_left=r_centre - dr_hts_tape / 2,
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
    """Plot the TF CroCo cable summary box."""
    get = _reader(mfile, scan)
    text = "$\\mathbf{CroCo \\ Cable:}$\n\n"
    for label, key in (
        ("Cable diameter", "dia_tf_turn_croco_cable"),
        ("Copper width", "dx_tf_croco_strand_copper"),
        ("Diameter of solder tape region", "dia_tf_croco_strand_tape_region"),
        ("Height of tape stack", "dx_tf_croco_strand_tape_stack"),
        ("Width of HTS tape / tape stack", "dr_tf_hts_tape"),
    ):
        text += f"{label}: {get(key) * 1e3:,.4f} mm\n"
    text += f"Number of HTS tape layers: {int(get('n_tf_croco_strand_hts_tapes'))}\n\n"
    for label, key in (
        ("Total copper area", "a_tf_croco_strand_copper_total"),
        ("Total hastelloy area", "a_tf_croco_strand_hastelloy"),
        ("Total solder area", "a_tf_croco_strand_solder"),
        ("Total superconductor area", "a_tf_croco_strand_rebco"),
        ("Total strand area", "a_tf_croco_strand"),
    ):
        text += f"{label}: {get(key) * 1e6:,.4f} mm²\n"
    axis.text(
        0.4, 0.4, text, **text_layout(fig, v_align="top"), bbox=box_style("#cccccc")
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
    """Plot quench current density and adiabatic hotspot temperature evolution.

    Both axes must be supplied together, or neither. The thermal integration,
    area scaling, fluence clipping and sampling follow the original model.

    Raises
    ------
    ValueError
        Axes not all provided
    """
    figure = None
    if axes_1 is None and axes_2 is None:
        figure, (axes_1, axes_2) = plt.subplots(2, 1, sharex=True)
    elif axes_1 is None or axes_2 is None:
        raise ValueError("Both axes_1 and axes_2 must be provided together, or neither.")

    fluence = np.clip(fluence, 0.0, 1.5e23)
    area_ratio = a_tf_turn_cable_space / a_tf_turn
    protection_arguments = {
        "tau_discharge": tau_discharge,
        "b_peak": b_peak,
        "f_a_cable_copper": f_a_cable_copper,
        "f_a_cable_space_helium": f_a_cable_space_helium,
        "temp_he_peak": temp_he_peak,
        "temp_quench_max": temp_quench_max,
        "cu_rrr": cu_rrr,
        "t_quench_detection": t_quench_detection,
    }
    j_max, j_max_1e23 = (
        area_ratio
        * quench.calculate_quench_protection_current_density(
            **protection_arguments,
            fluence=value,
        )
        for value in (fluence, 1e23)
    )
    t_end = max(4.0 * tau_discharge, t_quench_detection + 4.0 * tau_discharge)
    times = np.linspace(0.0, t_end, n_points)
    decay = np.exp(-(times - t_quench_detection) / tau_discharge)
    profiles = [
        np.where(times < t_quench_detection, initial, initial * decay)
        for initial in (j_max, j_max_1e23, j_operating)
    ]
    f_cu_cable = (1.0 - f_a_cable_space_helium) * f_a_cable_copper
    f_sc_cable = (1.0 - f_a_cable_space_helium) * (1.0 - f_a_cable_copper)
    integral_arguments = {
        "temp_he_peak": temp_he_peak,
        "temp_quench_max": temp_quench_max,
        "field": b_peak,
        "rrr": cu_rrr,
        "f_a_cable_space_helium": f_a_cable_space_helium,
        "f_cu_cable": f_cu_cable,
        "f_sc_cable": f_sc_cable,
    }
    lookups = [
        quench._build_cumulative_quench_integral(**integral_arguments, fluence=value)
        for value in (fluence, 1e23)
    ]
    dt = times[1] - times[0]
    temperatures = []
    for profile, lookup_index in zip(profiles, (0, 1, 0), strict=True):
        temp_array, cumulative_integral = lookups[lookup_index]
        scaled_integral = area_ratio**2 * f_cu_cable * cumulative_integral
        miit = np.cumsum(profile**2) * dt
        temperatures.append(np.interp(miit, scaled_integral, temp_array))

    current_labels = (
        f"Max allowed current density for protection (fluence = {fluence:.2e} n/m²)",
        "Max allowed current density for protection (fluence = 1e23 n/m²)",
        "Operating current density",
    )
    temperature_labels = (
        f"Hotspot temperature at protection limit (fluence = {fluence:.2e} n/m²)",
        "Hotspot temperature at protection limit (fluence = 1e23 n/m²)",
        "Operating hotspot temperature",
    )
    styles = (("darkorange", "-"), ("darkorange", "--"), ("blue", "-"))
    for ax, values, labels in (
        (axes_1, profiles, current_labels),
        (axes_2, temperatures, temperature_labels),
    ):
        for values_i, label, (colour, linestyle) in zip(
            values, labels, styles, strict=True
        ):
            ax.plot(
                times,
                values_i,
                color=colour,
                linewidth=2,
                linestyle=linestyle,
                label=label,
            )
        ax.grid(True, alpha=0.3)

    axes_1.axvline(
        t_quench_detection,
        color="crimson",
        linestyle="--",
        linewidth=1.5,
        label=f"Detection time ({t_quench_detection:.1f} s)",
    )
    axes_1.axvspan(
        0, t_quench_detection, alpha=0.08, color="crimson", label="Pre-detection phase"
    )
    axes_1.set_ylabel("Current density [A/m²]")
    axes_1.set_title(
        "TF Coil Quench Protection: Current Density and Hotspot Temperature Evolution"
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

    tau_time = t_quench_detection + tau_discharge
    tau_j = j_max * np.exp(-1)
    tau_temp = float(np.interp(tau_time, times, temperatures[0]))
    for ax, value, label in (
        (axes_1, tau_j, f"$J$ at $\\tau_{{\\text{{discharge}}}}$ ({tau_j:.2e} A/m²)"),
        (axes_2, tau_temp, f"$T$ at $\\tau_{{\\text{{discharge}}}}$ ({tau_temp:.1f} K)"),
    ):
        ax.axvline(
            tau_time,
            color="forestgreen",
            linestyle="--",
            linewidth=1.5,
            label="$t_{\\text{detect}} + \\tau_{\\text{discharge}}$"
            f" ({tau_time:.2f} s)",
        )
        ax.axhline(value, color="forestgreen", linestyle=":", linewidth=1.5, label=label)
        ax.legend(fontsize=9)
        ax.minorticks_on()
    if figure is not None:
        figure.tight_layout()
    else:
        plt.tight_layout()
    if show:
        plt.show()
