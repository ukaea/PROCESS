#!/usr/bin/env python
"""

  PROCESS plot_proc using process_io_lib functions and MFILE.DAT

  James Morris
  13/04/2014
  CCFE
  Revised by Michael Kovari, 7/1/2016

  24/11/2021: Global dictionary variables moved within the functions
              to avoid cyclic dependencies. This is because the dicts
              generation script imports, and inspects, process.

"""

import os
import sys
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np

import process.io.mfile as mf

from process.geometry.shield_geometry import (
    shield_geometry_single_null,
    shield_geometry_double_null,
)
from process.geometry.plasma_geometry import plasma_geometry
from process.geometry.vacuum_vessel_geometry import (
    vacuum_vessel_geometry_single_null,
    vacuum_vessel_geometry_double_null,
)
from process.geometry.blanket_geometry import (
    blanket_geometry_single_null,
    blanket_geometry_double_null,
)
from process.geometry.cryostat_geometry import cryostat_geometry
from process.geometry.tfcoil_geometry import (
    tfcoil_geometry_rectangular_shape,
    tfcoil_geometry_d_shape,
)
from process.geometry.pfcoil_geometry import pfcoil_geometry
from process.geometry.firstwall_geometry import (
    first_wall_geometry_single_null,
    first_wall_geometry_double_null,
)
from process.impurity_radiation import read_impurity_file
from process.io.python_fortran_dicts import get_dicts

if os.name == "posix" and "DISPLAY" not in os.environ:
    matplotlib.use("Agg")
matplotlib.rcParams["figure.max_open_warning"] = 40

if sys.version_info >= (3, 7):
    from importlib import resources
else:
    import importlib_resources as resources


solenoid = "pink"
cscompression = "red"
tfc = "cyan"
thermal_shield = "gray"
vessel = "green"
shield = "green"
blanket = "magenta"
plasma = "khaki"
cryostat = "red"
firstwall = "darkblue"
winding = "blue"
nbshield_colour = "gray"

thin = 0.0

RADIAL_BUILD = [
    "bore",
    "ohcth",
    "precomp",
    "gapoh",
    "tfcth",
    "tftsgap",
    "thshield_ib",
    "gapds",
    "d_vv_in",
    "shldith",
    "vvblgapi",
    "blnkith",
    "fwith",
    "scrapli",
    "rminori",
    "rminoro",
    "scraplo",
    "fwoth",
    "blnkoth",
    "vvblgapo",
    "shldoth",
    "d_vv_out",
    "gapsto",
    "thshield_ob",
    "tftsgap",
    "tfthko",
]

vertical_upper = [
    "rminor*kappa",
    "vgaptop",
    "fwtth",
    "blnktth",
    "vvblgap",
    "shldtth",
    "d_vv_top",
    "vgap2",
    "thshield_vb",
    "tftsgap",
    "tfcth",
]

vertical_lower = [
    "rminor*kappa",
    "vgap",
    "divfix",
    "shldlth",
    "d_vv_bot",
    "vgap2",
    "thshield_vb",
    "tftsgap",
    "tfcth",
]

ANIMATION_INFO = [
    ("rmajor", "Major radius", "m"),
    ("rminor", "Minor radius", "m"),
    ("aspect", "Aspect ratio", ""),
]

rtangle = np.pi / 2
rtangle2 = 2 * rtangle


def plot_plasma(axis, mfile_data, scan):
    """Plots the plasma boundary arcs.

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use

    """

    r_0 = mfile_data.data["rmajor"].get_scan(scan)
    a = mfile_data.data["rminor"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)
    kappa_95 = mfile_data.data["kappa95"].get_scan(scan)
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)

    pg = plasma_geometry(
        r_0=r_0,
        a=a,
        triang_95=triang_95,
        kappa_95=kappa_95,
        i_single_null=i_single_null,
    )

    axis.plot(pg.rs[0], pg.zs[0], color="black")
    axis.plot(pg.rs[1], pg.zs[1], color="black")
    axis.fill_betweenx(
        pg.zs[0],
        pg.rs[0],
        pg.rs[1],
        where=(pg.rs[1] < pg.rs[0])
        & (pg.zs[0] > (-a * pg.kappa))
        & (pg.zs[0] < (a * pg.kappa)),
        color=plasma,
    )
    axis.fill_betweenx(
        pg.zs[0], pg.rs[0], pg.rs[1], where=(pg.rs[1] > pg.rs[0]), color="none"
    )


def plot_centre_cross(axis, mfile_data, scan):
    """Function to plot centre cross on plot

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
    """
    rmajor = mfile_data.data["rmajor"].get_scan(scan)
    axis.plot(
        [rmajor - 0.25, rmajor + 0.25, rmajor, rmajor, rmajor],
        [0, 0, 0, 0.25, -0.25],
        color="black",
    )


def cumulative_radial_build(section, mfile_data, scan):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Arguments:
        section --> section of the radial build to go up to
        mfile_data --> MFILE data object
        scan --> scan number to use

    Returns:
        cumulative_build --> cumulative radial build up to section given

    """
    complete = False
    cumulative_build = 0
    for item in RADIAL_BUILD:
        if item == "rminori" or item == "rminoro":
            cumulative_build += mfile_data.data["rminor"].get_scan(scan)
        elif item == "vvblgapi" or item == "vvblgapo":
            cumulative_build += mfile_data.data["vvblgap"].get_scan(scan)
        elif "d_vv_in" in item:
            cumulative_build += mfile_data.data["d_vv_in"].get_scan(scan)
        elif "d_vv_out" in item:
            cumulative_build += mfile_data.data["d_vv_out"].get_scan(scan)
        else:
            cumulative_build += mfile_data.data[item].get_scan(scan)
        if item == section:
            complete = True
            break

    if complete is False:
        print("radial build parameter ", section, " not found")
    return cumulative_build


def cumulative_radial_build2(section, mfile_data, scan):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Arguments:
        section --> section of the radial build to go up to
        mfile_data --> MFILE data object
        scan --> scan number to use

    Returns:
        cumulative_build --> cumulative radial build up to and including
                             section given
        previous         --> cumulative radial build up to section given

    """
    cumulative_build = 0
    build = 0
    for item in RADIAL_BUILD:
        if item == "rminori" or item == "rminoro":
            build = mfile_data.data["rminor"].get_scan(scan)
        elif item == "vvblgapi" or item == "vvblgapo":
            build = mfile_data.data["vvblgap"].get_scan(scan)
        elif "d_vv_in" in item:
            build = mfile_data.data["d_vv_in"].get_scan(scan)
        elif "d_vv_out" in item:
            build = mfile_data.data["d_vv_out"].get_scan(scan)
        else:
            build = mfile_data.data[item].get_scan(scan)
        cumulative_build += build
        if item == section:
            break
    previous = cumulative_build - build
    return (cumulative_build, previous)


def poloidal_cross_section(axis, mfile_data, scan, demo_ranges):
    """Function to plot poloidal cross-section

    Arguments:
      axis --> axis object to add plot to
      mfile_data --> MFILE data object
      scan --> scan number to use

    """

    axis.set_xlabel("R / m")
    axis.set_ylabel("Z / m")
    axis.set_title("Poloidal cross-section")

    plot_vacuum_vessel(axis, mfile_data, scan)
    plot_shield(axis, mfile_data, scan)
    plot_blanket(axis, mfile_data, scan)
    plot_firstwall(axis, mfile_data, scan)

    plot_plasma(axis, mfile_data, scan)
    plot_centre_cross(axis, mfile_data, scan)
    plot_cryostat(axis, mfile_data, scan)

    plot_tf_coils(axis, mfile_data, scan)
    plot_pf_coils(axis, mfile_data, scan)

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        axis.set_ylim([-15, 15])
        axis.set_xlim([0, 20])

    # Adapatative ranges
    else:
        axis.set_xlim([0, axis.get_xlim()[1]])
    # ---


def plot_cryostat(axis, mfile_data, scan):
    """Function to plot cryostat in poloidal cross-section"""

    rects = cryostat_geometry(rdewex=rdewex, ddwex=ddwex, zdewex=zdewex)

    for rec in rects:
        axis.add_patch(
            patches.Rectangle(
                xy=(rec.anchor_x, rec.anchor_z),
                width=rec.width,
                height=rec.height,
                facecolor=cryostat,
            )
        )


def color_key(axis):
    """Function to plot the colour key
    Arguments:
      axis --> object to add plot to
    """
    axis.set_ylim([0, 10])
    axis.set_xlim([0, 10])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    axis.text(-5, 10, "CS coil", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 9.7], 1, 0.4, lw=0, facecolor=solenoid))

    axis.text(-5, 9, "CS comp", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 8.7], 1, 0.4, lw=0, facecolor=cscompression))

    axis.text(-5, 8, "TF coil", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 7.7], 1, 0.4, lw=0, facecolor=tfc))

    axis.text(-5, 7, "Th shield", ha="left", va="top", size="medium")
    axis.add_patch(
        patches.Rectangle([0.2, 6.7], 1, 0.4, lw=0, facecolor=thermal_shield)
    )

    axis.text(-5, 6, "VV & shield", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 5.7], 1, 0.4, lw=0, facecolor=vessel))

    axis.text(-5, 5, "Blanket", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 4.7], 1, 0.4, lw=0, facecolor=blanket))

    axis.text(-5, 4, "First wall", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 3.7], 1, 0.4, lw=0, facecolor=firstwall))

    axis.text(-5, 3, "Plasma", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, 2.7], 1, 0.4, lw=0, facecolor=plasma))

    axis.text(-5, 2, "PF coils", ha="left", va="top", size="medium")
    axis.add_patch(
        patches.Rectangle([0.2, 1.7], 1, 0.4, lw=1, facecolor="none", edgecolor="black")
    )

    axis.text(-5, 1, "NB duct shield", ha="left", va="top", size="medium")
    axis.add_patch(
        patches.Rectangle([0.2, 0.7], 1, 0.4, lw=0, facecolor=nbshield_colour)
    )

    axis.text(-5, 0.1, "cryostat", ha="left", va="top", size="medium")
    axis.add_patch(patches.Rectangle([0.2, -0.3], 1, 0.4, lw=0, facecolor=cryostat))


def toroidal_cross_section(axis, mfile_data, scan, demo_ranges):
    """Function to plot toroidal cross-section
    Arguments:
      axis --> axis object to add plot to
      mfile_data --> MFILE data object
      scan --> scan number to use
    """

    # Check for Copper magnets
    if "i_tf_sup" in mfile_data.data.keys():
        i_tf_sup = int(mfile_data.data["i_tf_sup"].get_scan(scan))
    else:
        i_tf_sup = int(1)

    if "i_tf_turns_integer" in mfile_data.data.keys():
        i_tf_turns_integer = int(mfile_data.data["i_tf_turns_integer"].get_scan(scan))
    else:
        i_tf_turns_integer = int(0)

    axis.set_xlabel("x / m")
    axis.set_ylabel("y / m")
    axis.set_title("Toroidal cross-section")

    arc(axis, rmajor, style="dashed")

    # Colour in the main components
    r2, r1 = cumulative_radial_build2("ohcth", mfile_data, scan)
    arc_fill(axis, r1, r2, color=solenoid)

    r2, r1 = cumulative_radial_build2("precomp", mfile_data, scan)
    arc_fill(axis, r1, r2, color=cscompression)

    r2, r1 = cumulative_radial_build2("tfcth", mfile_data, scan)
    arc_fill(axis, r1, r2, color=tfc)

    r2, r1 = cumulative_radial_build2("thshield_ib", mfile_data, scan)
    arc_fill(axis, r1, r2, color=thermal_shield)

    r2, r1 = cumulative_radial_build2("d_vv_in", mfile_data, scan)
    arc_fill(axis, r1, r2, color=vessel)

    r2, r1 = cumulative_radial_build2("shldith", mfile_data, scan)
    arc_fill(axis, r1, r2, color=vessel)

    r2, r1 = cumulative_radial_build2("blnkith", mfile_data, scan)
    arc_fill(axis, r1, r2, color=blanket)

    r2, r1 = cumulative_radial_build2("fwith", mfile_data, scan)
    arc_fill(axis, r1, r2, color=firstwall)

    arc_fill(axis, rmajor - rminor, rmajor + rminor, color=plasma)

    r2, r1 = cumulative_radial_build2("fwoth", mfile_data, scan)
    arc_fill(axis, r1, r2, color=firstwall)

    r2, r1 = cumulative_radial_build2("blnkoth", mfile_data, scan)
    arc_fill(axis, r1, r2, color=blanket)

    r2, r1 = cumulative_radial_build2("shldoth", mfile_data, scan)
    arc_fill(axis, r1, r2, color=shield)

    r2, r1 = cumulative_radial_build2("d_vv_out", mfile_data, scan)
    arc_fill(axis, r1, r2, color=vessel)

    r2, r1 = cumulative_radial_build2("thshield_ob", mfile_data, scan)
    arc_fill(axis, r1, r2, color=thermal_shield)

    arc_fill(axis, rdewex, rdewex + ddwex, color=cryostat)

    # Segment the TF coil inboard
    # Calculate centrelines
    n = int(n_tf / 4) + 1
    spacing = 2 * np.pi / n_tf
    i = np.arange(0, n)

    ang = i * spacing
    angl = ang - spacing / 2
    angu = ang + spacing / 2
    r1, null = cumulative_radial_build2("gapoh", mfile_data, scan)
    r2, null = cumulative_radial_build2("tfcth", mfile_data, scan)
    r4, r3 = cumulative_radial_build2("tfthko", mfile_data, scan)

    # Coil width
    w = r2 * np.tan(spacing / 2)
    xi = r1 * np.cos(angl)
    yi = r1 * np.sin(angl)
    xo = r2 * np.cos(angl)
    yo = r2 * np.sin(angl)
    axis.plot((xi, xo), (yi, yo), color="black")
    xi = r1 * np.cos(angu)
    yi = r1 * np.sin(angu)
    xo = r2 * np.cos(angu)
    yo = r2 * np.sin(angu)
    axis.plot((xi, xo), (yi, yo), color="black")
    # Annotate plot.
    axis.text(
        rmajor * np.cos(0.3),
        rmajor * np.sin(0.3),
        "plasma",
        fontsize=(12),
        ha="center",
        va="center",
    )
    axis.text(
        (rdewex + ddwex) / 1.41,
        (rdewex + ddwex) / 1.41,
        "cryostat",
        fontsize=(10),
        ha="left",
        va="bottom",
    )

    for item in i:
        # Neutral beam shielding
        TF_outboard(
            axis,
            item,
            n_tf=n_tf,
            r3=r3,
            r4=r4,
            w=w + nbshield,
            facecolor=nbshield_colour,
        )
        # Overlay TF coil segments
        TF_outboard(axis, item, n_tf=n_tf, r3=r3, r4=r4, w=w, facecolor="cyan")

    # Winding pack : inboard (superconducor only)
    if i_tf_sup == 1:
        # Inboard
        if i_tf_turns_integer == 1:
            rect = patches.Rectangle(
                [r1 + thkcas + tinstf, 0], dr_tf_wp, wwp1 / 2, lw=0, facecolor=winding
            )
            axis.add_patch(rect)
        else:
            rect = patches.Rectangle(
                [r1 + thkcas + tinstf, 0],
                dr_tf_wp / 2,
                wwp2 / 2,
                lw=0,
                facecolor=winding,
            )
            axis.add_patch(rect)

            rect = patches.Rectangle(
                [r1 + thkcas + tinstf + dr_tf_wp / 2, 0],
                dr_tf_wp / 2,
                wwp1 / 2,
                lw=0,
                facecolor=winding,
            )
            axis.add_patch(rect)

        # Outboard
        if i_tf_turns_integer == 1:
            rect = patches.Rectangle(
                [r3 + casthi + tinstf, 0], dr_tf_wp, wwp1 / 2, lw=0, facecolor=winding
            )
            axis.add_patch(rect)
        else:
            rect = patches.Rectangle(
                [r3 + casthi + tinstf, 0],
                dr_tf_wp / 2,
                wwp1 / 2,
                lw=0,
                facecolor=winding,
            )
            axis.add_patch(rect)
            rect = patches.Rectangle(
                [r3 + casthi + tinstf + dr_tf_wp / 2, 0],
                dr_tf_wp / 2,
                wwp2 / 2,
                lw=0,
                facecolor=winding,
            )
            axis.add_patch(rect)

    iefrf = mfile_data.data["iefrf"].get_scan(scan)
    if (iefrf == 5) or (iefrf == 8):
        # Neutral beam geometry
        a = w
        b = tfthko
        c = beamwd + 2 * nbshield
        d = r3
        e = np.sqrt(a**2 + (d + b) ** 2)
        # Coordinates of the inner and outer edges of the beam at its tangency point
        rinner = rtanbeam - beamwd
        router = rtanbeam + beamwd
        beta = np.arccos(rinner / e)
        xinner = rinner * np.cos(beta)
        yinner = rinner * np.sin(beta)
        xouter = router * np.cos(beta)
        youter = router * np.sin(beta)
        # Corner of TF coils
        xcorner = r4
        ycorner = w + nbshield
        axis.plot(
            [xinner, xcorner], [yinner, ycorner], linestyle="dotted", color="black"
        )
        x = xcorner + c * np.cos(beta) - nbshield * np.cos(beta)
        y = ycorner + c * np.sin(beta) - nbshield * np.sin(beta)
        axis.plot([xouter, x], [youter, y], linestyle="dotted", color="black")

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        axis.set_ylim([0, 20])
        axis.set_xlim([0, 20])

    # Adapatative ranges
    else:
        axis.set_ylim([0.0, axis.get_ylim()[1]])
        axis.set_xlim([0.0, axis.get_xlim()[1]])
    # ---


def TF_outboard(axis, item, n_tf, r3, r4, w, facecolor):
    spacing = 2 * np.pi / n_tf
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
    path = Path(verts, closed=True)
    patch = patches.PathPatch(path, facecolor=facecolor, lw=0)
    axis.add_patch(patch)


def arc(axis, r, theta1=0, theta2=rtangle, style="solid"):
    """Plots an arc.

    Arguments

    axis: plot object
    r: radius
    theta1: starting polar angle
    theta2: finishing polar angle

    """
    angs = np.linspace(theta1, theta2)
    xs = r * np.cos(angs)
    ys = r * np.sin(angs)
    axis.plot(xs, ys, linestyle=style, color="black", lw=0.2)


def arc_fill(axis, r1, r2, color="pink"):
    """Fills the space between two quarter circles.

    Arguments

    axis: plot object
    r1, r2 radii to be filled

    """
    angs = np.linspace(0, rtangle, endpoint=True)
    xs1 = r1 * np.cos(angs)
    ys1 = r1 * np.sin(angs)
    angs = np.linspace(rtangle, 0, endpoint=True)
    xs2 = r2 * np.cos(angs)
    ys2 = r2 * np.sin(angs)
    verts = list(zip(xs1, ys1))
    verts.extend(list(zip(xs2, ys2)))
    endpoint = [(r2, 0)]
    verts.extend(endpoint)
    path = Path(verts, closed=True)
    patch = patches.PathPatch(path, facecolor=color, lw=0)
    axis.add_patch(patch)


def plot_nprofile(prof, demo_ranges):
    """Function to plot density profile
    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel("r/a")
    prof.set_ylabel(r"$n_{e}\cdot 10^{19}$ $[\mathrm{m}^{-3}]$")
    prof.set_title("Density profile")

    if ipedestal == 1:
        rhocore1 = np.linspace(0, 0.95 * rhopedn)
        rhocore2 = np.linspace(0.95 * rhopedn, rhopedn)
        rhocore = np.append(rhocore1, rhocore2)
        ncore = neped + (ne0 - neped) * (1 - rhocore**2 / rhopedn**2) ** alphan

        rhosep = np.linspace(rhopedn, 1)
        nsep = nesep + (neped - nesep) * (1 - rhosep) / (1 - min(0.9999, rhopedn))

        rho = np.append(rhocore, rhosep)
        ne = np.append(ncore, nsep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        ne = ne0 * (1 - rho**2) ** alphan
    ne = ne / 1e19
    prof.plot(rho, ne)

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    prof.set_xlim([0, 1])
    if demo_ranges:
        prof.set_ylim([0, 20])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])
    # ---


def plot_tprofile(prof, demo_ranges):
    """Function to plot temperature profile
    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel("r/a")
    prof.set_ylabel("$T_{e}$ [keV]")
    prof.set_title("Temperature profile")

    if ipedestal == 1:
        rhocore1 = np.linspace(0, 0.9 * rhopedt)
        rhocore2 = np.linspace(0.9 * rhopedt, rhopedt)
        rhocore = np.append(rhocore1, rhocore2)
        tcore = teped + (te0 - teped) * (1 - (rhocore / rhopedt) ** tbeta) ** alphat

        rhosep = np.linspace(rhopedt, 1)
        tsep = tesep + (teped - tesep) * (1 - rhosep) / (1 - min(0.9999, rhopedt))

        rho = np.append(rhocore, rhosep)
        te = np.append(tcore, tsep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        te = te0 * (1 - rho**2) ** alphat
    prof.plot(rho, te)

    # Ranges
    # ---
    prof.set_xlim([0, 1])
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim([0, 50])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])
    # ---


def plot_qprofile(prof, demo_ranges):
    """Function to plot q profile, formula taken from Nevins bootstrap model.

    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel("r/a")
    prof.set_ylabel("q(r)")
    prof.set_title("q profile")

    rho = np.linspace(0, 1)
    q_r_nevin = q0 + (q95 - q0) * (rho + rho * rho + rho**3) / (3.0)
    q_r_sauter = q0 + (q95 - q0) * (rho * rho)

    prof.plot(rho, q_r_nevin, label="Nevins")
    prof.plot(rho, q_r_sauter, label="Sauter")
    prof.legend()

    # Ranges
    # ---
    prof.set_xlim([0, 1])
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim([0, 10])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])
    # ---


def read_imprad_data(skiprows, data_path):
    """Function to read all data needed for creation of radiation profile

    Arguments:
        skiprows --> number of rows to skip when reading impurity data files
        data_path --> path to impurity data
    """
    label = [
        "H_",
        "He",
        "Be",
        "C_",
        "N_",
        "O_",
        "Ne",
        "Si",
        "Ar",
        "Fe",
        "Ni",
        "Kr",
        "Xe",
        "W_",
    ]
    lzdata = [0.0 for x in range(len(label))]
    # DATAFILENAME = p DATAPATH +

    for i in range(len(label)):
        file_iden = data_path + label[i].ljust(3, "_")

        Te = None
        lz = None
        zav = None

        for header in read_impurity_file(file_iden + "lz_tau.dat"):
            if "Te[eV]" in header.content:
                Te = np.asarray(header.data, dtype=float)

            if "infinite confinement" in header.content:
                lz = np.asarray(header.data, dtype=float)
        for header in read_impurity_file(file_iden + "z_tau.dat"):
            if "infinite confinement" in header.content:
                zav = np.asarray(header.data, dtype=float)

        lzdata[i] = np.column_stack((Te, lz, zav))

    # then switch string to floats
    impdata = np.array(lzdata, dtype=float)
    return impdata


def synchrotron_rad():
    """Function for Synchrotron radiation power calculation from Albajar, Nuclear Fusion 41 (2001) 665
      Fidone, Giruzzi, Granata, Nuclear Fusion 41 (2001) 1755

    Arguments:
    """
    # tbet is betaT in Albajar, not to be confused with plasma beta

    tbet = 2.0
    # rpow is the(1-Rsyn) power dependence based on plasma shape
    # (see Fidone)
    rpow = 0.62
    kap = vol / (2.0 * 3.1415**2 * rmajor * rminor**2)

    # No account is taken of pedestal profiles here, other than use of
    # the correct ne0 and te0...
    de2o = 1.0e-20 * ne0
    pao = 6.04e3 * (rminor * de2o) / bt
    gfun = 0.93 * (1.0 + 0.85 * np.exp(-0.82 * rmajor / rminor))
    kfun = (alphan + 3.87e0 * alphat + 1.46) ** (-0.79)
    kfun = kfun * (1.98 + alphat) ** 1.36 * tbet**2.14
    kfun = kfun * (tbet**1.53 + 1.87 * alphat - 0.16) ** (-1.33)
    dum = 1.0 + 0.12 * (te0 / (pao**0.41)) * (1.0 - ssync) ** 0.41
    # Very high T modification, from Fidone
    dum = dum ** (-1.51)

    psync = 3.84e-8 * (1.0e0 - ssync) ** rpow * rmajor * rminor**1.38
    psync = psync * kap**0.79 * bt**2.62 * de2o**0.38
    psync = psync * te0 * (16.0 + te0) ** 2.61 * dum * gfun * kfun

    # psyncpv should be per unit volume
    # Albajar gives it as total
    psyncpv = psync / vol
    print("psyncpv = ", psyncpv * vol)  # matches the out.dat file

    return psyncpv


def plot_radprofile(prof, mfile_data, scan, impp, demo_ranges) -> float:
    """Function to plot radiation profile, formula taken from ???.

    Arguments:
      prof --> axis object to add plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use
      impp --> impurity path
    """

    prof.set_xlabel("r/a")
    prof.set_ylabel(r"$P_{\mathrm{rad}}$ $[\mathrm{MW.m}^{-3}]$")
    prof.set_title("Radiation profile")

    # read in the impurity data
    imp_data = read_imprad_data(2, impp)

    # find impurity densities
    imp_frac = np.array(
        [
            mfile_data.data["fimp(01)"].get_scan(scan),
            mfile_data.data["fimp(02)"].get_scan(scan),
            mfile_data.data["fimp(03)"].get_scan(scan),
            mfile_data.data["fimp(04)"].get_scan(scan),
            mfile_data.data["fimp(05)"].get_scan(scan),
            mfile_data.data["fimp(06)"].get_scan(scan),
            mfile_data.data["fimp(07)"].get_scan(scan),
            mfile_data.data["fimp(08)"].get_scan(scan),
            mfile_data.data["fimp(09)"].get_scan(scan),
            mfile_data.data["fimp(10)"].get_scan(scan),
            mfile_data.data["fimp(11)"].get_scan(scan),
            mfile_data.data["fimp(12)"].get_scan(scan),
            mfile_data.data["fimp(13)"].get_scan(scan),
            mfile_data.data["fimp(14)"].get_scan(scan),
        ]
    )

    if ipedestal == 0:
        # Intialise the radius
        rho = np.linspace(0, 1.0)

        # The density profile
        ne = ne0 * (1 - rho**2) ** alphan

        # The temperature profile
        te = te0 * (1 - rho**2) ** alphat

    if ipedestal == 1:
        # Intialise the normalised radius
        rhoped = (rhopedn + rhopedt) / 2.0
        rhocore1 = np.linspace(0, 0.95 * rhoped)
        rhocore2 = np.linspace(0.95 * rhoped, rhoped)
        rhocore = np.append(rhocore1, rhocore2)
        rhosep = np.linspace(rhoped, 1)
        rho = np.append(rhocore, rhosep)

        # The density and temperature profile
        # done in such away as to allow for plotting pedestals
        # with different rhopedn and rhopedt
        ne = np.zeros(rho.shape[0])
        te = np.zeros(rho.shape[0])
        for q in range(rho.shape[0]):
            if rho[q] <= rhopedn:
                ne[q] = (
                    neped + (ne0 - neped) * (1 - rho[q] ** 2 / rhopedn**2) ** alphan
                )
            else:
                ne[q] = nesep + (neped - nesep) * (1 - rho[q]) / (
                    1 - min(0.9999, rhopedn)
                )

            if rho[q] <= rhopedt:
                te[q] = (
                    teped + (te0 - teped) * (1 - (rho[q] / rhopedt) ** tbeta) ** alphat
                )
            else:
                te[q] = tesep + (teped - tesep) * (1 - rho[q]) / (
                    1 - min(0.9999, rhopedt)
                )

        # ncore = neped + (ne0-neped) * (1-rhocore**2/rhopedn**2)**alphan
        # nsep = nesep + (neped-nesep) * (1-rhosep)/(1-min(0.9999, rhopedn))
        # ne = np.append(ncore, nsep)

        # The temperatue profile
        # tcore = teped + (te0-teped) * (1-(rhocore/rhopedt)**tbeta)**alphat
        # tsep = tesep + (teped-tesep)* (1-rhosep)/(1-min(0.9999,rhopedt))
        # te = np.append(tcore,tsep)

    # Intailise the radiation profile arrays
    pimpden = np.zeros([imp_data.shape[0], te.shape[0]])
    lz = np.zeros([imp_data.shape[0], te.shape[0]])
    prad = np.zeros(te.shape[0])

    # psyncpv = synchrotron_rad()

    # Intailise the impurity radiation profile
    for k in range(te.shape[0]):
        for i in range(imp_data.shape[0]):
            if te[k] <= imp_data[i][0][0]:
                lz[i][k] = imp_data[i][0][1]
            elif te[k] >= imp_data[i][imp_data.shape[1] - 1][0]:
                lz[i][k] = imp_data[i][imp_data.shape[1] - 1][1]
            else:
                for j in range(imp_data.shape[1] - 1):
                    # Linear interpolation in log-log space
                    if (te[k] > imp_data[i][j][0]) and (te[k] <= imp_data[i][j + 1][0]):
                        yi = np.log(imp_data[i][j][1])
                        xi = np.log(imp_data[i][j][0])
                        c = (np.log(imp_data[i][j + 1][1]) - yi) / (
                            np.log(imp_data[i][j + 1][0]) - xi
                        )
                        lz[i][k] = np.exp(yi + c * (np.log(te[k]) - xi))
                        # Zav[i][k] = imp_data[i][j][2]
            # The impurity radiation
            pimpden[i][k] = imp_frac[i] * ne[k] * ne[k] * lz[i][k]
            # The Bremsstrahlung
            # pbremden[i][k] = imp_frac[i] * ne[k] * ne[k] * Zav[i][k] * Zav[i][k] * 5.355e-37 * np.sqrt(te[k])

        for l in range(imp_data.shape[0]):  # noqa: E741
            prad[k] = prad[k] + pimpden[l][k] * 2.0e-6
            # pbrem[k] = pbrem[k] + pbremden[l][k] * 2.0e-6

    # benchmark prad again outfile so mod prad
    # pbremint = (rho[1:] * pbrem[1:]) @ drho
    # pradint = prad[1:] @ drho * 2.0e-5
    # pbremint = pbrem[1:] @ drho * 2.0e-5

    # print('prad = ',prad)
    # print('pbrem = ',pbrem)
    # print(1.0e32*lz[12])
    # print('pradpv = ',pradint)
    # print('pbrempv = ',pbremint)
    # print('pbremmw = ',pbremint*vol)
    # print('pradmw = ', pradint*vol, 'MW') # pimp = pline + pbrem

    prof.plot(rho, prad, label="Total")
    prof.plot(rho, pimpden[0] * 2.0e-6, label="H")
    prof.plot(rho, pimpden[1] * 2.0e-6, label="He")
    if imp_frac[2] > 1.0e-30:
        prof.plot(rho, pimpden[2] * 2.0e-6, label="Be")
    if imp_frac[3] > 1.0e-30:
        prof.plot(rho, pimpden[3] * 2.0e-6, label="C")
    if imp_frac[4] > 1.0e-30:
        prof.plot(rho, pimpden[4] * 2.0e-6, label="N")
    if imp_frac[5] > 1.0e-30:
        prof.plot(rho, pimpden[5] * 2.0e-6, label="O")
    if imp_frac[6] > 1.0e-30:
        prof.plot(rho, pimpden[6] * 2.0e-6, label="Ne")
    if imp_frac[7] > 1.0e-30:
        prof.plot(rho, pimpden[7] * 2.0e-6, label="Si")
    if imp_frac[8] > 1.0e-30:
        prof.plot(rho, pimpden[8] * 2.0e-6, label="Ar")
    if imp_frac[9] > 1.0e-30:
        prof.plot(rho, pimpden[9] * 2.0e-6, label="Fe")
    if imp_frac[10] > 1.0e-30:
        prof.plot(rho, pimpden[10] * 2.0e-6, label="Ni")
    if imp_frac[11] > 1.0e-30:
        prof.plot(rho, pimpden[11] * 2.0e-6, label="Kr")
    if imp_frac[12] > 1.0e-30:
        prof.plot(rho, pimpden[12] * 2.0e-6, label="Xe")
    if imp_frac[13] > 1.0e-30:
        prof.plot(rho, pimpden[13] * 2.0e-6, label="W")
    prof.legend()

    # Ranges
    # ---
    prof.set_xlim([0, 1])
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim([0, 0.5])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])
    # ---


def plot_vacuum_vessel(axis, mfile_data, scan):
    """Function to plot vacuum vessel

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
    """
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)

    # Outer side (furthest from plasma)
    radx_outer = (
        cumulative_radial_build("d_vv_out", mfile_data, scan)
        + cumulative_radial_build("gapds", mfile_data, scan)
    ) / 2.0
    rminx_outer = (
        cumulative_radial_build("d_vv_out", mfile_data, scan)
        - cumulative_radial_build("gapds", mfile_data, scan)
    ) / 2.0

    # Inner side (nearest to the plasma)
    radx_inner = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        + cumulative_radial_build("d_vv_in", mfile_data, scan)
    ) / 2.0
    rminx_inner = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        - cumulative_radial_build("d_vv_in", mfile_data, scan)
    ) / 2.0

    if i_single_null == 1:
        vvg_single_null = vacuum_vessel_geometry_single_null(
            cumulative_upper=cumulative_upper,
            upper=upper,
            triang=triang_95,
            radx_outer=radx_outer,
            rminx_outer=rminx_outer,
            radx_inner=radx_inner,
            rminx_inner=rminx_inner,
            cumulative_lower=cumulative_lower,
            lower=lower,
        )

        axis.plot(
            vvg_single_null.rs,
            vvg_single_null.zs,
            color="black",
            lw=thin,
        )

        axis.fill(
            vvg_single_null.rs,
            vvg_single_null.zs,
            color=vessel,
        )

    if i_single_null == 0:
        vvg_double_null = vacuum_vessel_geometry_double_null(
            cumulative_lower=cumulative_lower,
            lower=lower,
            radx_inner=radx_inner,
            radx_outer=radx_outer,
            rminx_inner=rminx_inner,
            rminx_outer=rminx_outer,
            triang=triang_95,
        )
        axis.plot(vvg_double_null.rs, vvg_double_null.zs, color="black", lw=thin)

        axis.fill(
            vvg_double_null.rs,
            vvg_double_null.zs,
            color=vessel,
        )


def plot_shield(axis, mfile_data, scan):
    """Function to plot shield

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
    """
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)

    # Side furthest from plasma
    radx_far = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        + cumulative_radial_build("d_vv_in", mfile_data, scan)
    ) / 2.0
    rminx_far = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        - cumulative_radial_build("d_vv_in", mfile_data, scan)
    ) / 2.0

    # Side nearest to the plasma
    radx_near = (
        cumulative_radial_build("vvblgapo", mfile_data, scan)
        + cumulative_radial_build("shldith", mfile_data, scan)
    ) / 2.0
    rminx_near = (
        cumulative_radial_build("vvblgapo", mfile_data, scan)
        - cumulative_radial_build("shldith", mfile_data, scan)
    ) / 2.0

    if i_single_null == 1:
        sg_single_null = shield_geometry_single_null(
            cumulative_upper=cumulative_upper,
            radx_far=radx_far,
            rminx_far=rminx_far,
            radx_near=radx_near,
            rminx_near=rminx_near,
            triang=triang_95,
            cumulative_lower=cumulative_lower,
        )
        axis.plot(sg_single_null.rs, sg_single_null.zs, color="black", lw=thin)
        axis.fill(sg_single_null.rs, sg_single_null.zs, color=shield)

    if i_single_null == 0:
        sg_double_null = shield_geometry_double_null(
            cumulative_lower=cumulative_lower,
            radx_far=radx_far,
            radx_near=radx_near,
            rminx_far=rminx_far,
            rminx_near=rminx_near,
            triang=triang_95,
        )
        axis.plot(sg_double_null.rs, sg_double_null.zs, color="black", lw=thin)
        axis.fill(sg_double_null.rs, sg_double_null.zs, color=shield)


def plot_blanket(axis, mfile_data, scan) -> None:
    """Function to plot blanket

    Arguments:
      axis --> axis object to plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use

    """
    # Single null: Draw top half from output
    # Double null: Reflect bottom half to top
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)
    blnktth = mfile_data.data["blnktth"].get_scan(scan)
    c_shldith = cumulative_radial_build("shldith", mfile_data, scan)
    c_blnkoth = cumulative_radial_build("blnkoth", mfile_data, scan)

    if i_single_null == 1:
        # Upper blanket: outer surface
        radx_outer = (
            cumulative_radial_build("blnkoth", mfile_data, scan)
            + cumulative_radial_build("vvblgapi", mfile_data, scan)
        ) / 2.0
        rminx_outer = (
            cumulative_radial_build("blnkoth", mfile_data, scan)
            - cumulative_radial_build("vvblgapi", mfile_data, scan)
        ) / 2.0

        # Upper blanket: inner surface
        radx_inner = (
            cumulative_radial_build("fwoth", mfile_data, scan)
            + cumulative_radial_build("blnkith", mfile_data, scan)
        ) / 2.0
        rminx_inner = (
            cumulative_radial_build("fwoth", mfile_data, scan)
            - cumulative_radial_build("blnkith", mfile_data, scan)
        ) / 2.0
        bg_single_null = blanket_geometry_single_null(
            radx_outer=radx_outer,
            rminx_outer=rminx_outer,
            radx_inner=radx_inner,
            rminx_inner=rminx_inner,
            cumulative_upper=cumulative_upper,
            triang=triang_95,
            cumulative_lower=cumulative_lower,
            blnktth=blnktth,
            c_shldith=c_shldith,
            c_blnkoth=c_blnkoth,
            blnkith=blnkith,
            blnkoth=blnkoth,
        )

        # Plot blanket
        axis.plot(
            bg_single_null.rs,
            bg_single_null.zs,
            color="black",
            lw=thin,
        )

        axis.fill(bg_single_null.rs, bg_single_null.zs, color=blanket)

    if i_single_null == 0:
        bg_double_null = blanket_geometry_double_null(
            cumulative_lower=cumulative_lower,
            triang=triang_95,
            blnktth=blnktth,
            c_shldith=c_shldith,
            c_blnkoth=c_blnkoth,
            blnkith=blnkith,
            blnkoth=blnkoth,
        )
        # Plot blanket
        axis.plot(bg_double_null.rs[0], bg_double_null.zs[0], color="black", lw=thin)
        axis.plot(bg_double_null.rs[1], bg_double_null.zs[1], color="black", lw=thin)
        axis.fill(bg_double_null.rs[0], bg_double_null.zs[0], color=blanket)
        axis.fill(bg_double_null.rs[1], bg_double_null.zs[1], color=blanket)


def plot_firstwall(axis, mfile_data, scan):
    """Function to plot first wall

    Arguments:
      axis --> axis object to plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use

    """
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)
    blnktth = mfile_data.data["blnktth"].get_scan(scan)
    tfwvt = mfile_data.data["fwtth"].get_scan(scan)
    c_blnkith = cumulative_radial_build("blnkith", mfile_data, scan)
    c_fwoth = cumulative_radial_build("fwoth", mfile_data, scan)

    if i_single_null == 1:
        # Upper first wall: outer surface
        radx_outer = (
            cumulative_radial_build("fwoth", mfile_data, scan)
            + cumulative_radial_build("blnkith", mfile_data, scan)
        ) / 2.0
        rminx_outer = (
            cumulative_radial_build("fwoth", mfile_data, scan)
            - cumulative_radial_build("blnkith", mfile_data, scan)
        ) / 2.0

        # Upper first wall: inner surface
        radx_inner = (
            cumulative_radial_build("scraplo", mfile_data, scan)
            + cumulative_radial_build("fwith", mfile_data, scan)
        ) / 2.0
        rminx_inner = (
            cumulative_radial_build("scraplo", mfile_data, scan)
            - cumulative_radial_build("fwith", mfile_data, scan)
        ) / 2.0

        fwg_single_null = first_wall_geometry_single_null(
            radx_outer=radx_outer,
            rminx_outer=rminx_outer,
            radx_inner=radx_inner,
            rminx_inner=rminx_inner,
            cumulative_upper=cumulative_upper,
            triang=triang_95,
            cumulative_lower=cumulative_lower,
            blnktth=blnktth,
            c_blnkith=c_blnkith,
            c_fwoth=c_fwoth,
            fwith=fwith,
            fwoth=fwoth,
            tfwvt=tfwvt,
        )

        # Plot first wall
        axis.plot(fwg_single_null.rs, fwg_single_null.zs, color="black", lw=thin)
        axis.fill(fwg_single_null.rs, fwg_single_null.zs, color=firstwall)

    if i_single_null == 0:
        fwg_double_null = first_wall_geometry_double_null(
            cumulative_lower=cumulative_lower,
            triang=triang_95,
            blnktth=blnktth,
            c_blnkith=c_blnkith,
            c_fwoth=c_fwoth,
            fwith=fwith,
            fwoth=fwoth,
            tfwvt=tfwvt,
        )
        # Plot blanket
        axis.plot(fwg_double_null.rs[0], fwg_double_null.zs[0], color="black", lw=thin)
        axis.plot(fwg_double_null.rs[1], fwg_double_null.zs[1], color="black", lw=thin)
        axis.fill(fwg_double_null.rs[0], fwg_double_null.zs[0], color=blanket)
        axis.fill(fwg_double_null.rs[1], fwg_double_null.zs[1], color=blanket)


def angle_check(angle1, angle2):
    """Function to perform TF coil angle check"""
    if angle1 > 1:
        angle1 = 1
    if angle1 < -1:
        angle1 = -1
    if angle2 > 1:
        angle2 = 1
    if angle2 < -1:
        angle2 = -1
    return angle1, angle2


def plot_tf_coils(axis, mfile_data, scan):
    """Function to plot TF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    # Arc points
    # MDK Only 4 points now required for elliptical arcs
    x1 = mfile_data.data["xarc(1)"].get_scan(scan)
    y1 = mfile_data.data["yarc(1)"].get_scan(scan)
    x2 = mfile_data.data["xarc(2)"].get_scan(scan)
    y2 = mfile_data.data["yarc(2)"].get_scan(scan)
    x3 = mfile_data.data["xarc(3)"].get_scan(scan)
    y3 = mfile_data.data["yarc(3)"].get_scan(scan)
    x4 = mfile_data.data["xarc(4)"].get_scan(scan)
    y4 = mfile_data.data["yarc(4)"].get_scan(scan)
    x5 = mfile_data.data["xarc(5)"].get_scan(scan)
    y5 = mfile_data.data["yarc(5)"].get_scan(scan)
    if y3 != 0:
        print("TF coil geometry: The value of yarc(3) is not zero, but should be.")

    # Check for TF coil shape
    if "i_tf_shape" in mfile_data.data.keys():
        i_tf_shape = int(mfile_data.data["i_tf_shape"].get_scan(scan))
    else:
        i_tf_shape = int(1)

    if i_tf_shape == 2:
        rects = tfcoil_geometry_rectangular_shape(
            x1=x1,
            x2=x2,
            x4=x4,
            x5=x5,
            y1=y1,
            y2=y2,
            y4=y4,
            y5=y5,
            tfcth=tfcth,
        )

        for rec in rects:
            axis.add_patch(
                patches.Rectangle(
                    xy=(rec.anchor_x, rec.anchor_z),
                    width=rec.width,
                    height=rec.height,
                    facecolor=tfc,
                )
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
            tfcth=tfcth,
            rtangle=rtangle,
            rtangle2=rtangle2,
        )

        for vert in verts:
            path = Path(vert, closed=True)
            patch = patches.PathPatch(path, facecolor=tfc, lw=0)
            axis.add_patch(patch)

        for rec in rects:
            axis.add_patch(
                patches.Rectangle(
                    xy=(rec.anchor_x, rec.anchor_z),
                    width=rec.width,
                    height=rec.height,
                    facecolor=tfc,
                )
            )


def plot_pf_coils(axis, mfile_data, scan):
    """Function to plot PF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use
    """

    coils_r = []
    coils_z = []
    coils_dr = []
    coils_dz = []
    coil_text = []

    bore = mfile_data.data["bore"].get_scan(scan)
    ohcth = mfile_data.data["ohcth"].get_scan(scan)
    ohdz = mfile_data.data["ohdz"].get_scan(scan)

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile_data.data.keys():
        if "rpf[" in item:
            number_of_coils += 1

    # Check for Central Solenoid
    if "iohcl" in mfile_data.data.keys():
        iohcl = mfile_data.data["iohcl"].get_scan(scan)
    else:
        iohcl = 1

    # If Central Solenoid present, ignore last entry in for loop
    # The last entry will be the OH coil in this case
    if iohcl == 0:
        noc = number_of_coils + 1
    else:
        noc = number_of_coils

    for coil in range(0, noc):
        coils_r.append(mfile_data.data["rpf[{:01}]".format(coil)].get_scan(scan))
        coils_z.append(mfile_data.data["zpf[{:01}]".format(coil)].get_scan(scan))
        coils_dr.append(mfile_data.data["pfdr({:01})".format(coil)].get_scan(scan))
        coils_dz.append(mfile_data.data["pfdz({:01})".format(coil)].get_scan(scan))
        coil_text.append(str(coil + 1))

    r_points, z_points, central_coil = pfcoil_geometry(
        coils_r=coils_r,
        coils_z=coils_z,
        coils_dr=coils_dr,
        coils_dz=coils_dz,
        bore=bore,
        ohcth=ohcth,
        ohdz=ohdz,
    )

    for i in range(len(coils_r)):
        axis.plot(r_points[i], z_points[i], color="black")
        axis.text(
            coils_r[i],
            coils_z[i],
            coil_text[i],
            ha="center",
            va="center",
            fontsize=5 * abs((coils_dr[i] * coils_dz[i])),
        )
    axis.add_patch(
        patches.Rectangle(
            xy=(central_coil.anchor_x, central_coil.anchor_z),
            width=central_coil.width,
            height=central_coil.height,
            facecolor="pink",
        )
    )


def plot_info(axis, data, mfile_data, scan):
    """Function to plot data in written form on a matplotlib plot.

    Arguments:
        axis --> axis object to plot to
        data --> plot information
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    eqpos = 0.7
    for i in range(len(data)):
        colorflag = "black"
        if mfile_data.data[data[i][0]].exists:
            if mfile_data.data[data[i][0]].var_flag == "ITV":
                colorflag = "red"
            elif mfile_data.data[data[i][0]].var_flag == "OP":
                colorflag = "blue"
        axis.text(0, -i, data[i][1], color=colorflag, ha="left", va="center")
        if isinstance(data[i][0], str):
            if data[i][0] == "":
                axis.text(eqpos, -i, "\n", ha="left", va="center")
            elif data[i][0][0] == "#":
                axis.text(
                    -0.05, -i, "{}\n".format(data[i][0][1:]), ha="left", va="center"
                )
            elif data[i][0][0] == "!":
                value = data[i][0][1:]
                axis.text(
                    0.4,
                    -i,
                    "-->  " + str(value.replace('"', "")) + " " + data[i][2],
                    ha="left",
                    va="center",
                )
            else:
                if mfile_data.data[data[i][0]].exists:
                    dat = mfile_data.data[data[i][0]].get_scan(scan)
                    if isinstance(dat, str):
                        value = dat
                    else:
                        value = "{:.4g}".format(
                            mfile_data.data[data[i][0]].get_scan(scan)
                        )
                    if "alpha" in data[i][0]:
                        value = str(float(value) + 1.0)
                    axis.text(
                        eqpos,
                        -i,
                        "= " + value + " " + data[i][2],
                        color=colorflag,
                        ha="left",
                        va="center",
                    )
                else:
                    mfile_data.data[data[i][0]].get_scan(-1)
                    axis.text(
                        eqpos,
                        -i,
                        "=" + "ERROR! Var missing",
                        color=colorflag,
                        ha="left",
                        va="center",
                    )
        else:
            dat = data[i][0]
            if isinstance(dat, str):
                value = dat
            else:
                value = "{:.4g}".format(data[i][0])
            axis.text(
                eqpos,
                -i,
                "= " + value + " " + data[i][2],
                color=colorflag,
                ha="left",
                va="center",
            )


def plot_header(axis, mfile_data, scan):
    """Function to plot header info: date, rutitle etc

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    data2 = [
        ("!" + str(mfile_data.data["runtitle"].get_scan(-1)), "Run title", ""),
        ("!" + str(mfile_data.data["procver"].get_scan(-1)), "PROCESS Version", ""),
        ("!" + mfile_data.data["date"].get_scan(-1), "Date:", ""),
        ("!" + mfile_data.data["time"].get_scan(-1), "Time:", ""),
        ("!" + mfile_data.data["username"].get_scan(-1), "User:", ""),
        (
            "!"
            + dicts["DICT_OPTIMISATION_VARS"][
                str(abs(int(mfile_data.data["minmax"].get_scan(-1))))
            ],
            "Optimising:",
            "",
        ),
    ]

    H = mfile_data.data["fimp(01)"].get_scan(scan)
    He = mfile_data.data["fimp(02)"].get_scan(scan)
    Be = mfile_data.data["fimp(03)"].get_scan(scan)
    C = mfile_data.data["fimp(04)"].get_scan(scan)
    N = mfile_data.data["fimp(05)"].get_scan(scan)
    O = mfile_data.data["fimp(06)"].get_scan(scan)  # noqa: E741
    Ne = mfile_data.data["fimp(07)"].get_scan(scan)
    Si = mfile_data.data["fimp(08)"].get_scan(scan)
    Ar = mfile_data.data["fimp(09)"].get_scan(scan)
    Fe = mfile_data.data["fimp(10)"].get_scan(scan)
    Ni = mfile_data.data["fimp(11)"].get_scan(scan)
    Kr = mfile_data.data["fimp(12)"].get_scan(scan)
    Xe = mfile_data.data["fimp(13)"].get_scan(scan)
    W = mfile_data.data["fimp(14)"].get_scan(scan)

    data = [("", "", ""), ("", "", "")]
    count = 0

    data = data + [(H, "D + T", "")]
    count += 1

    data = data + [(He, "He", "")]
    count += 1
    if Be > 1e-10:
        data = data + [(Be, "Be", "")]
        count += +1
    if C > 1e-10:
        data = data + [(C, "C", "")]
        count += 1
    if N > 1e-10:
        data = data + [(N, "N", "")]
        count += 1
    if O > 1e-10:
        data = data + [(O, "O", "")]
        count += 1
    if Ne > 1e-10:
        data = data + [(Ne, "Ne", "")]
        count += 1
    if Si > 1e-10:
        data = data + [(Si, "Si", "")]
        count += 1
    if Ar > 1e-10:
        data = data + [(Ar, "Ar", "")]
        count += 1
    if Fe > 1e-10:
        data = data + [(Fe, "Fe", "")]
        count += 1
    if Ni > 1e-10:
        data = data + [(Ni, "Ni", "")]
        count += 1
    if Kr > 1e-10:
        data = data + [(Kr, "Kr", "")]
        count += 1
    if Xe > 1e-10:
        data = data + [(Xe, "Xe", "")]
        count += 1
    if W > 1e-10:
        data = data + [(W, "W", "")]
        count += 1

    if count > 11:
        data = [("", "", ""), ("", "", ""), ("", "More than 11 impurities", "")]
    else:
        axis.text(-0.05, -6.4, "Plasma composition:", ha="left", va="center")
        axis.text(
            -0.05,
            -7.2,
            "Number densities relative to electron density:",
            ha="left",
            va="center",
        )
    data2 = data2 + data

    axis.text(-0.05, -12.6, "Colour Legend:", ha="left", va="center")
    axis.text(0.0, -13.4, "ITR", color="red", ha="left", va="center")
    axis.text(0.0, -14.2, "OP", color="blue", ha="left", va="center")

    plot_info(axis, data2, mfile_data, scan)


def plot_geometry_info(axis, mfile_data, scan):
    """Function to plot geometry info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Geometry:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    in_blanket_thk = mfile_data.data["shldith"].get_scan(scan) + mfile_data.data[
        "blnkith"
    ].get_scan(scan)
    out_blanket_thk = mfile_data.data["shldoth"].get_scan(scan) + mfile_data.data[
        "blnkoth"
    ].get_scan(scan)

    data = [
        ("rmajor", "$R_0$", "m"),
        ("rminor", "a", "m"),
        ("aspect", "A", ""),
        ("kappa95", r"$\kappa_{95}$", ""),
        ("triang95", r"$\delta_{95}$", ""),
        ("sarea", "Surface area", "m$^2$"),
        ("vol", "Plasma volume", "m$^3$"),
        ("n_tf", "No. of TF coils", ""),
        (in_blanket_thk, "inboard blanket+shield", "m"),
        (out_blanket_thk, "ouboard blanket+shield", "m"),
        ("powfmw", "Fusion power", "MW"),
        ("bigq", "$Q$", ""),
        ("", "", ""),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_physics_info(axis, mfile_data, scan):
    """Function to plot geometry info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Physics:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    nong = mfile_data.data["dnla"].get_scan(scan) / mfile_data.data[
        "dlimit(7)"
    ].get_scan(scan)

    dnz = mfile_data.data["dnz"].get_scan(scan) / mfile_data.data["dene"].get_scan(scan)

    tepeak = mfile_data.data["te0"].get_scan(scan) / mfile_data.data["te"].get_scan(
        scan
    )

    nepeak = mfile_data.data["ne0"].get_scan(scan) / mfile_data.data["dene"].get_scan(
        scan
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    if "plhthresh" in mfile_data.data.keys():
        pthresh = mfile_data.data["plhthresh"].get_scan(scan)
    else:
        pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)

    data = [
        ("plascur/1d6", "$I_p$", "MA"),
        ("bt", "Vacuum $B_T$ at $R_0$", "T"),
        ("q95", r"$q_{\mathrm{95}}$", ""),
        ("normalised_thermal_beta", r"$\beta_N$, thermal", "% m T MA$^{-1}$"),
        ("normalised_toroidal_beta", r"$\beta_N$, toroidal", "% m T MA$^{-1}$"),
        ("thermal_poloidal_beta", r"$\beta_P$, thermal", ""),
        ("betap", r"$\beta_P$, total", ""),
        ("te", r"$< t_e >$", "keV"),
        ("dene", r"$< n_e >$", "m$^{-3}$"),
        (nong, r"$< n_{\mathrm{e,line}} >/n_G$", ""),
        (tepeak, r"$T_{e0}/ < T_e >$", ""),
        (nepeak, r"$n_{e0}/ < n_{\mathrm{e, vol}} >$", ""),
        ("zeff", r"$Z_{\mathrm{eff}}$", ""),
        (dnz, r"$n_Z/ < n_{\mathrm{e, vol}} >$", ""),
        ("taueff", r"$\tau_e$", "s"),
        ("hfact", "H-factor", ""),
        (pthresh, "H-mode threshold", "MW"),
        ("tauelaw", "Scaling law", ""),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_magnetics_info(axis, mfile_data, scan):
    """Function to plot magnet info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Check for Copper magnets
    if "i_tf_sup" in mfile_data.data.keys():
        i_tf_sup = int(mfile_data.data["i_tf_sup"].get_scan(scan))
    else:
        i_tf_sup = int(1)

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Coil currents etc:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile_data.data.keys():
        if "rpf[" in item:
            number_of_coils += 1

    pf_info = []
    for i in range(1, number_of_coils):
        if i % 2 != 0:
            pf_info.append(
                (
                    mfile_data.data["ric[{:01}]".format(i)].get_scan(scan),
                    "PF {}".format(i),
                )
            )

    if len(pf_info) > 2:
        pf_info_3_a = pf_info[2][0]
        pf_info_3_b = pf_info[2][1]
    else:
        pf_info_3_a = ""
        pf_info_3_b = ""

    tburn = mfile_data.data["tburn"].get_scan(scan) / 3600.0

    if "i_tf_bucking" in mfile_data.data.keys():
        i_tf_bucking = int(mfile_data.data["i_tf_bucking"].get_scan(scan))
    else:
        i_tf_bucking = int(1)

    # Get superconductor material (i_tf_sc_mat)
    # If i_tf_sc_mat not present, assume resistive
    if "i_tf_sc_mat" in mfile_data.data.keys():
        i_tf_sc_mat = int(mfile_data.data["i_tf_sc_mat"].get_scan(scan))
    else:
        i_tf_sc_mat = int(0)

    if i_tf_sc_mat > 0:
        tftype = dicts["DICT_TF_TYPE"][
            str(int(mfile_data.data["i_tf_sc_mat"].get_scan(scan)))
        ]
    else:
        tftype = "Resistive"

    vssoft = mfile_data.data["vsres"].get_scan(scan) + mfile_data.data[
        "vsind"
    ].get_scan(scan)

    sig_case = 1.0e-6 * mfile_data.data[
        "sig_tf_tresca_max({})".format(i_tf_bucking)
    ].get_scan(scan)
    sig_cond = 1.0e-6 * mfile_data.data[
        "sig_tf_tresca_max({})".format(i_tf_bucking + 1)
    ].get_scan(scan)
    alstrtf = 1.0e-6 * mfile_data.data["alstrtf"].get_scan(scan)

    if i_tf_sup == 1:
        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vstot", "Available flux swing", "Wb"),
            (tburn, "Burn time", "hrs"),
            ("", "", ""),
            ("#TF coil type is {}".format(tftype), "", ""),
            ("bmaxtfrp", "Peak field at conductor (w. rip.)", "T"),
            ("iooic", r"I/I$_{\mathrm{crit}}$", ""),
            ("tmargtf", "TF Temperature margin", "K"),
            ("tmargoh", "CS Temperature margin", "K"),
            (sig_cond, "TF Cond max TRESCA stress", "MPa"),
            (sig_case, "TF Case max TRESCA stress", "MPa"),
            (alstrtf, "Allowable stress", "Pa"),
            ("whttf/n_tf", "Mass per TF coil", "kg"),
        ]

    else:
        n_tf = mfile_data.data["n_tf"].get_scan(scan)
        prescp = 1.0e-6 * mfile_data.data["prescp"].get_scan(scan)
        presleg = 1.0e-6 * mfile_data.data["presleg"].get_scan(scan)
        pres_joints = 1.0e-6 * mfile_data.data["pres_joints"].get_scan(scan)
        fcoolcp = 100.0 * mfile_data.data["fcoolcp"].get_scan(scan)

        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vstot", "Available flux swing", "Wb"),
            (tburn, "Burn time", "hrs"),
            ("", "", ""),
            ("#TF coil type is {}".format(tftype), "", ""),
            ("bmaxtf", "Peak field at conductor (w. rip.)", "T"),
            ("ritfc", "TF coil currents sum", "A"),
            ("", "", ""),
            ("#TF coil forces/stresses", "", ""),
            (sig_cond, "TF conductor max TRESCA stress", "MPa"),
            (sig_case, "TF bucking max TRESCA stress", "MPa"),
            (alstrtf, "conductor Allowable stress", "Pa"),
            (fcoolcp, "CP cooling fraction", "%"),
            ("vcool", "Maximum coolant flow speed", "m.s$^{-1}$"),
            (prescp, "CP Resisitive heating", "MW"),
            (presleg * n_tf, "legs Resisitive heating (all legs)", "MW"),
            (pres_joints, "TF joints resisitive heating ", "MW"),
        ]

    plot_info(axis, data, mfile_data, scan)


def plot_power_info(axis, mfile_data, scan):
    """Function to plot power info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Power flows:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    gross_eff = 100.0 * (
        mfile_data.data["pgrossmw"].get_scan(scan)
        / mfile_data.data["pthermmw"].get_scan(scan)
    )

    net_eff = 100.0 * (
        (
            mfile_data.data["pgrossmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
        / (
            mfile_data.data["pthermmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
    )

    plant_eff = 100.0 * (
        mfile_data.data["pnetelmw"].get_scan(scan)
        / mfile_data.data["powfmw"].get_scan(scan)
    )

    # Define appropriate pedestal and impurity parameters
    coredescription = ("coreradius", "Normalised radius of 'core' region", "")
    if ipedestal == 1:
        ped_height = ("neped", "Electron density at pedestal", "m$^{-3}$")
        ped_pos = ("rhopedn", "r/a at density pedestal", "")
    else:
        ped_height = ("", "No pedestal model used", "")
        ped_pos = ("", "", "")

    crypmw = mfile_data.data["crypmw"].get_scan(scan)

    data = [
        ("wallmw", "Nominal neutron wall load", "MW m$^{-2}$"),
        coredescription,
        ped_height,
        ped_pos,
        ("ralpne", "Helium fraction", ""),
        ("pinnerzoneradmw", "inner zone radiation", "MW"),
        ("pradmw", "Total radiation in LCFS", "MW"),
        ("pnucblkt", "Nuclear heating in blanket", "MW"),
        ("pnucshld", "Nuclear heating in shield", "MW"),
        (crypmw, "TF cryogenic power", "MW"),
        ("pdivt", "Power to divertor", "MW"),
        ("divlife", "Divertor life", "years"),
        ("pthermmw", "Primary (high grade) heat", "MW"),
        (gross_eff, "Gross cycle efficiency", "%"),
        (net_eff, "Net cycle efficiency", "%"),
        ("pgrossmw", "Gross electric power", "MW"),
        ("pnetelmw", "Net electric power", "MW"),
        (
            plant_eff,
            "Fusion-to-electric efficiency "
            + r"$\frac{P_{\mathrm{e,net}}}{P_{\mathrm{fus}}}$",
            "%",
        ),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_current_drive_info(axis, mfile_data, scan):
    """Function to plot current drive info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1
    iefrf = mfile_data.data["iefrf"].get_scan(scan)
    nbi = False
    ecrh = False
    ebw = False
    if (iefrf == 5) or (iefrf == 8):
        nbi = True
        axis.text(-0.05, 1, "Neutral Beam Current Drive:", ha="left", va="center")
    if (iefrf == 3) or (iefrf == 7) or (iefrf == 10) or (iefrf == 11) or (iefrf == 13):
        ecrh = True
        axis.text(-0.05, 1, "Electron Cyclotron Current Drive:", ha="left", va="center")
    if iefrf == 12:
        ebw = True
        axis.text(-0.05, 1, "Electron Bernstein Wave Drive:", ha="left", va="center")
    if (iefrf == 1) or (iefrf == 2) or (iefrf == 4) or (iefrf == 6) or (iefrf == 9):
        print(
            "Options 1, 2, 4, 6 and 9 not implemented yet in this python script plot_proc.py\n"
        )
        print("NEEDS TO BE IMPLEMENTED in plot_current_drive_info subroutine!!\n")

    if "iefrffix" in mfile_data.data.keys():
        secondary_heating = ""
        iefrffix = mfile_data.data["iefrffix"].get_scan(scan)

        if (iefrffix == 5) or (iefrffix == 8):
            secondary_heating = "NBI"
        if (
            (iefrffix == 3)
            or (iefrffix == 7)
            or (iefrffix == 10)
            or (iefrffix == 11)
            or (iefrffix == 13)
        ):
            secondary_heating = "ECH"
        if iefrffix == 12:
            secondary_heating = "EBW"
        if (
            (iefrffix == 1)
            or (iefrffix == 2)
            or (iefrffix == 4)
            or (iefrffix == 6)
            or (iefrffix == 9)
        ):
            print(
                "Options 1, 2, 4, 6 and 9 not implemented yet in this python script plot_proc.py\n"
            )
            print("NEEDS TO BE IMPLEMENTED in plot_current_drive_info subroutine!!\n")

    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    pinjie = mfile_data.data["pinjmw"].get_scan(scan)
    pinjmwfix = mfile_data.data["pinjmwfix"].get_scan(scan)
    pdivt = mfile_data.data["pdivt"].get_scan(scan)
    pdivr = pdivt / mfile_data.data["rmajor"].get_scan(scan)

    pdivnr = (
        1.0e20
        * mfile_data.data["pdivt"].get_scan(scan)
        / (
            mfile_data.data["rmajor"].get_scan(scan)
            * mfile_data.data["dene"].get_scan(scan)
        )
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    if "plhthresh" in mfile_data.data.keys():
        pthresh = mfile_data.data["plhthresh"].get_scan(scan)
    else:
        pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)
    flh = pdivt / pthresh

    powerht = mfile_data.data["powerht"].get_scan(scan)
    psync = mfile_data.data["psyncpv*vol"].get_scan(scan)
    pbrem = mfile_data.data["pinnerzoneradmw"].get_scan(scan)
    hfact = mfile_data.data["hfact"].get_scan(scan)
    hstar = hfact * (powerht / (powerht + psync + pbrem)) ** 0.31

    if ecrh:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootipf", "Bootstrap fraction", ""),
            ("faccd", "Auxiliary fraction", ""),
            ("facoh", "Inductive fraction", ""),
            ("powerht", "Plasma heating used for H factor", "MW"),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{<n> R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data.keys():
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if nbi:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootipf", "Bootstrap fraction", ""),
            ("faccd", "Auxiliary fraction", ""),
            ("facoh", "Inductive fraction", ""),
            ("gamnb", "NB gamma", "$10^{20}$ A W$^{-1}$ m$^{-2}$"),
            ("enbeam", "NB energy", "keV"),
            ("powerht", "Plasma heating used for H factor", "MW"),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{<n> R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data.keys():
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if ebw:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootipf", "Bootstrap fraction", ""),
            ("faccd", "Auxiliary fraction", ""),
            ("facoh", "Inductive fraction", ""),
            ("powerht", "Plasma heating used for H factor", "MW"),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{<n> R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data.keys():
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    coe = mfile_data.data["coe"].get_scan(scan)
    if coe == 0.0:
        data.append(("", "", ""))
        data.append(("#Costs", "", ""))
        data.append(("", "Cost output not selected", ""))
    else:
        data.append(("", "", ""))
        data.append(("#Costs", "", ""))
        data.append((coe, "Cost of electricity", r"\$/MWh"))

    plot_info(axis, data, mfile_data, scan)


def main_plot(
    fig1,
    fig2,
    m_file_data,
    scan,
    imp="../data/lz_non_corona_14_elements/",
    demo_ranges=False,
):
    """Function to create radial and vertical build plot on given figure.

    Arguments:
      fig1 --> figure object to add plot to.
      fig2 --> figure object to add plot to.
      m_file_data --> MFILE.DAT data to read
      scan --> scan to read from MFILE.DAT
      imp --> path to impurity data
    """

    # Checking the impurity data folder
    # Get path to impurity data dir
    # TODO use Path objects throughout module, not strings
    with resources.path(
        "process.data.lz_non_corona_14_elements", "Ar_lz_tau.dat"
    ) as imp_path:
        data_folder = str(imp_path.parent) + "/"

    if os.path.isdir(data_folder):
        imp = data_folder
    else:
        print(
            "\033[91m Warning : Impossible to recover impurity data, try running the macro in the main/utility folder"
        )
        print("          -> No impurity plot done\033[0m")

    # Plot poloidal cross-section
    plot_1 = fig2.add_subplot(221, aspect="equal")
    poloidal_cross_section(plot_1, m_file_data, scan, demo_ranges)

    # Plot toroidal cross-section
    plot_2 = fig2.add_subplot(222, aspect="equal")
    toroidal_cross_section(plot_2, m_file_data, scan, demo_ranges)

    # Plot color key
    plot_3 = fig2.add_subplot(241)
    color_key(plot_3)

    # Plot density profiles
    plot_4 = fig2.add_subplot(234)  # , aspect= 0.05)
    plot_nprofile(plot_4, demo_ranges)

    # Plot temperature profiles
    plot_5 = fig2.add_subplot(235)  # , aspect= 1/35)
    plot_tprofile(plot_5, demo_ranges)

    # plot_qprofile(plot_6)
    plot_6 = fig2.add_subplot(236)  # , aspect=2)
    if os.path.isdir(imp):
        plot_radprofile(plot_6, m_file_data, scan, imp, demo_ranges)

    # plot_7 =
    # plot_radprofile(plot_7)

    # Setup params for text plots
    plt.rcParams.update({"font.size": 8})

    # Plot header info
    plot_1 = fig1.add_subplot(231)
    plot_header(plot_1, m_file_data, scan)

    # Geometry
    plot_2 = fig1.add_subplot(232)
    plot_geometry_info(plot_2, m_file_data, scan)

    # Physics
    plot_3 = fig1.add_subplot(233)
    plot_physics_info(plot_3, m_file_data, scan)

    # Magnetics
    plot_4 = fig1.add_subplot(234)
    plot_magnetics_info(plot_4, m_file_data, scan)

    # power/flow economics
    plot_5 = fig1.add_subplot(235)
    plot_power_info(plot_5, m_file_data, scan)

    # Current drive
    plot_6 = fig1.add_subplot(236)
    plot_current_drive_info(plot_6, m_file_data, scan)
    fig1.subplots_adjust(wspace=0.25)


def save_plots(m_file_data, scan):
    """Function to recreate and save individual plots."""

    fig = plt.figure(figsize=(12, 9), dpi=80)

    # Plot poloidal cross-section
    pol = fig.add_subplot(111, aspect="equal")
    poloidal_cross_section(pol, m_file_data, scan)

    # Plot TF coils
    plot_tf_coils(pol, m_file_data, scan)

    # Plot PF coils
    plot_pf_coils(pol, m_file_data, scan)

    fig.savefig("psection.svg", format="svg", dpi=1200)

    # Plot toroidal cross-section
    fig = plt.figure(figsize=(12, 9), dpi=80)
    tor = fig.add_subplot(222, aspect="equal")
    toroidal_cross_section(tor)
    fig.savefig("tsection.svg", format="svg", dpi=1200)

    # Plot color key
    fig = plt.figure(figsize=(12, 9), dpi=80)
    plot = fig.add_subplot(241)
    color_key(plot)
    fig.savefig("color_key.svg", format="svg", dpi=1200)

    # Plot profiles
    fig = plt.figure(figsize=(12, 9), dpi=80)
    plot = fig.add_subplot(223, aspect=0.05)
    plot_nprofile(plot, False)
    fig.savefig("nprofile.svg", format="svg", dpi=1200)

    fig = plt.figure(figsize=(12, 9), dpi=80)
    plot = fig.add_subplot(224, aspect=1 / 35)
    plot_tprofile(plot, False)
    fig.savefig("tprofile.svg", format="svg", dpi=1200)


def test(f):
    """Test Function

    :param f: filename to test
    """

    try:
        # read MFILE
        m_file = mf.MFile(filename=f)
        scan = -1

        # Check for Copper magnets
        if "i_tf_sup" in m_file.data.keys():
            i_tf_sup = int(m_file.data["i_tf_sup"].get_scan(scan))
        else:
            i_tf_sup = int(1)

        # Check integer turns
        if "i_tf_turns_integer" in m_file.data.keys():
            i_tf_turns_integer = int(m_file.data["i_tf_turns_integer"].get_scan(scan))
        else:
            i_tf_turns_integer = int(0)

        global bore
        bore = m_file.data["bore"].get_scan(scan)
        global ohcth
        ohcth = m_file.data["ohcth"].get_scan(scan)
        global gapoh
        gapoh = m_file.data["gapoh"].get_scan(scan)
        global tfcth
        tfcth = m_file.data["tfcth"].get_scan(scan)
        global gapds
        gapds = m_file.data["gapds"].get_scan(scan)
        global d_vv_in
        d_vv_in = m_file.data["d_vv_in"].get_scan(scan)
        global shldith
        shldith = m_file.data["shldith"].get_scan(scan)
        global blnkith
        blnkith = m_file.data["blnkith"].get_scan(scan)
        global fwith
        fwith = m_file.data["fwith"].get_scan(scan)
        global scrapli
        scrapli = m_file.data["scrapli"].get_scan(scan)
        global rmajor
        rmajor = m_file.data["rmajor"].get_scan(scan)
        global rminor
        rminor = m_file.data["rminor"].get_scan(scan)
        global scraplo
        scraplo = m_file.data["scraplo"].get_scan(scan)
        global fwoth
        fwoth = m_file.data["fwoth"].get_scan(scan)
        global blnkoth
        blnkoth = m_file.data["blnkoth"].get_scan(scan)
        global shldoth
        shldoth = m_file.data["shldoth"].get_scan(scan)
        global d_vv_out
        d_vv_out = m_file.data["d_vv_out"].get_scan(scan)
        global gapsto
        gapsto = m_file.data["gapsto"].get_scan(scan)
        global tfthko
        tfthko = m_file.data["tfthko"].get_scan(scan)
        global rdewex
        rdewex = m_file.data["rdewex"].get_scan(scan)
        global ddwex
        ddwex = m_file.data["ddwex"].get_scan(scan)
        global zdewex
        zdewex = m_file.data["zdewex"].get_scan(scan)
        global n_tf
        n_tf = m_file.data["n_tf"].get_scan(scan)

        if i_tf_sup == 1:
            global wwp1
            wwp1 = m_file.data["wwp1"].get_scan(scan)
            if i_tf_turns_integer == 0:
                global wwp2
                wwp2 = m_file.data["wwp2"].get_scan(scan)
            global dr_tf_wp
            dr_tf_wp = m_file.data["dr_tf_wp"].get_scan(scan)
            global tinstf
            tinstf = m_file.data["tinstf"].get_scan(scan)
            global thkcas
            thkcas = m_file.data["thkcas"].get_scan(scan)

            # To be re-inergrated to resistives when in-plane stresses is integrated
            global casthi
            casthi = m_file.data["casthi"].get_scan(scan)

        global nbshield
        nbshield = m_file.data["nbshield"].get_scan(scan)
        global rtanbeam
        rtanbeam = m_file.data["rtanbeam"].get_scan(scan)
        global rtanmax
        rtanmax = m_file.data["rtanmax"].get_scan(scan)
        global beamwd
        beamwd = m_file.data["beamwd"].get_scan(scan)
        # # Pedestal profile parameters
        global ipedestal
        ipedestal = m_file.data["ipedestal"].get_scan(scan)
        global neped
        neped = m_file.data["neped"].get_scan(scan)
        global nesep
        nesep = m_file.data["nesep"].get_scan(scan)
        global rhopedn
        rhopedn = m_file.data["rhopedn"].get_scan(scan)
        global rhopedt
        rhopedt = m_file.data["rhopedt"].get_scan(scan)
        global tbeta
        tbeta = m_file.data["tbeta"].get_scan(scan)
        global teped
        teped = m_file.data["teped"].get_scan(scan)
        global tesep
        tesep = m_file.data["tesep"].get_scan(scan)
        global alphan
        alphan = m_file.data["alphan"].get_scan(scan)
        global alphat
        alphat = m_file.data["alphat"].get_scan(scan)
        global ne0
        ne0 = m_file.data["ne0"].get_scan(scan)
        global te0
        te0 = m_file.data["te0"].get_scan(scan)
        # # Plasma
        global triang
        triang = m_file.data["triang95"].get_scan(scan)
        global alphaj
        alphaj = m_file.data["alphaj"].get_scan(scan)
        global q0
        q0 = m_file.data["q0"].get_scan(scan)
        global q95
        q95 = m_file.data["q95"].get_scan(scan)

        # Build the dictionaries of radial and vertical build values and cumulative
        # values
        global radial
        radial = dict()
        global cumulative_radial
        cumulative_radial = dict()
        subtotal = 0
        for item in RADIAL_BUILD:
            if item == "rminori" or item == "rminoro":
                build = m_file.data["rminor"].get_scan(scan)
            elif item == "vvblgapi" or item == "vvblgapo":
                build = m_file.data["vvblgap"].get_scan(scan)
            elif "d_vv_in" in item:
                build = m_file.data["d_vv_in"].get_scan(scan)
            elif "d_vv_out" in item:
                build = m_file.data["d_vv_out"].get_scan(scan)
            else:
                build = m_file.data[item].get_scan(scan)

        radial[item] = build
        subtotal += build
        cumulative_radial[item] = subtotal

        global upper
        upper = dict()
        global cumulative_upper
        cumulative_upper = dict()
        subtotal = 0
        for item in vertical_upper:
            upper[item] = m_file.data[item].get_scan(scan)
            subtotal += upper[item]
            cumulative_upper[item] = subtotal

        global lower
        lower = dict()
        global cumulative_lower
        cumulative_lower = dict()
        subtotal = 0
        for item in vertical_lower:
            lower[item] = m_file.data[item].get_scan(scan)
            subtotal -= lower[item]
            cumulative_lower[item] = subtotal

        colour_dict = {}
        colour_dict["ohcth"] = solenoid
        colour_dict["tfcth"] = tfc
        colour_dict["thshield_ib"] = thermal_shield
        colour_dict["thshield_ob"] = thermal_shield
        colour_dict["thshield_vb"] = thermal_shield
        colour_dict["d_vv_in"] = vessel
        colour_dict["d_vv_out"] = vessel
        colour_dict["d_vv_top"] = vessel
        colour_dict["d_vv_bot"] = vessel
        colour_dict["shldith"] = shield
        colour_dict["blnkith"] = blanket
        colour_dict["rminor"] = plasma
        colour_dict["fwith"] = firstwall
        colour_dict["fwoth"] = firstwall

        # create main plot
        page1 = plt.figure(figsize=(12, 9), dpi=80)
        page2 = plt.figure(figsize=(12, 9), dpi=80)

        # run main_plot
        main_plot(page1, page2, m_file, scan=scan)

        # with bpdf.PdfPages(args.o) as pdf:
        # with bpdf.PdfPages("ref.SUMMARY.pdf") as pdf:
        #    pdf.savefig(page1)
        #    pdf.savefig(page2)
        # plt.show()

        # # tidy up to avoid memory issues
        # del page1
        # del page2
        plt.close(page1)
        plt.close(page2)

        return True
    except Exception:
        print("FTest failure for file : {}".format(f))
        return False


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Produces a two page summary of the PROCESS MFILE output, using the MFILE.  "
        "For info please see https://github.com/ukaea/PROCESS?tab=readme-ov-file#contacts "
    )

    parser.add_argument(
        "-f",
        metavar="FILENAME",
        type=str,
        default="",
        help="specify input/output file path",
    )
    parser.add_argument("-s", "--show", help="show plot", action="store_true")

    parser.add_argument("-n", type=int, help="Which scan to plot?")

    parser.add_argument(
        "-d",
        "--DEMO_ranges",
        help="Uses the DEMO dimensions as ranges for all graphics",
        action="store_true",
    )

    return parser.parse_args(args)


def main(args=None):
    # TODO The use of globals here isn't ideal, but is required to get main()
    # working with minimal changes. Should be converted to class structure
    args = parse_args(args)

    # read MFILE
    if args.f != "":
        m_file = mf.MFile(args.f)
    else:
        m_file = mf.MFile("MFILE.DAT")

    if args.n:
        scan = args.n
    else:
        scan = -1

    if args.DEMO_ranges:
        demo_ranges = True
    else:
        demo_ranges = False

    # Check for Copper magnets
    if "i_tf_sup" in m_file.data.keys():
        i_tf_sup = int(m_file.data["i_tf_sup"].get_scan(scan))
    else:
        i_tf_sup = int(1)

    # Check integer turns
    if "i_tf_turns_integer" in m_file.data.keys():
        i_tf_turns_integer = int(m_file.data["i_tf_turns_integer"].get_scan(scan))
    else:
        i_tf_turns_integer = int(0)

    global bore
    global ohcth
    global gapoh
    global tfcth
    global gapds
    global ddwi
    global shldith
    global blnkith
    global fwith
    global scrapli
    global rmajor
    global rminor
    global scraplo
    global fwoth
    global blnkoth
    global shldoth
    global ddwi
    global gapsto
    global tfthko
    global rdewex
    global zdewex
    global ddwex

    bore = m_file.data["bore"].get_scan(scan)
    ohcth = m_file.data["ohcth"].get_scan(scan)
    gapoh = m_file.data["gapoh"].get_scan(scan)
    tfcth = m_file.data["tfcth"].get_scan(scan)
    gapds = m_file.data["gapds"].get_scan(scan)
    shldith = m_file.data["shldith"].get_scan(scan)
    blnkith = m_file.data["blnkith"].get_scan(scan)
    fwith = m_file.data["fwith"].get_scan(scan)
    scrapli = m_file.data["scrapli"].get_scan(scan)
    rmajor = m_file.data["rmajor"].get_scan(scan)
    rminor = m_file.data["rminor"].get_scan(scan)
    scraplo = m_file.data["scraplo"].get_scan(scan)
    fwoth = m_file.data["fwoth"].get_scan(scan)
    blnkoth = m_file.data["blnkoth"].get_scan(scan)
    shldoth = m_file.data["shldoth"].get_scan(scan)
    gapsto = m_file.data["gapsto"].get_scan(scan)
    tfthko = m_file.data["tfthko"].get_scan(scan)
    rdewex = m_file.data["rdewex"].get_scan(scan)
    zdewex = m_file.data["zdewex"].get_scan(scan)
    ddwex = m_file.data["ddwex"].get_scan(scan)

    # Magnets related
    global n_tf
    global wwp1
    global wwp2
    global dr_tf_wp
    global tinstf
    global thkcas
    global casthi

    n_tf = m_file.data["n_tf"].get_scan(scan)
    if i_tf_sup == 1:  # If superconducting magnets
        wwp1 = m_file.data["wwp1"].get_scan(scan)
        if i_tf_turns_integer == 0:
            wwp2 = m_file.data["wwp2"].get_scan(scan)
        dr_tf_wp = m_file.data["dr_tf_wp"].get_scan(scan)
        tinstf = m_file.data["tinstf"].get_scan(scan)
        thkcas = m_file.data["thkcas"].get_scan(scan)

        # To be re-inergrated to resistives when in-plane stresses is integrated
        casthi = m_file.data["casthi"].get_scan(scan)

    global nbshield
    global rtanbeam
    global rtanmax
    global beamwd

    nbshield = m_file.data["nbshield"].get_scan(scan)
    rtanbeam = m_file.data["rtanbeam"].get_scan(scan)
    rtanmax = m_file.data["rtanmax"].get_scan(scan)
    beamwd = m_file.data["beamwd"].get_scan(scan)

    # Pedestal profile parameters
    global ipedestal
    global neped
    global nesep
    global rhopedn
    global rhopedt
    global tbeta
    global teped
    global tesep
    global alphan
    global alphat
    global ne0
    global te0

    ipedestal = m_file.data["ipedestal"].get_scan(scan)
    neped = m_file.data["neped"].get_scan(scan)
    nesep = m_file.data["nesep"].get_scan(scan)
    rhopedn = m_file.data["rhopedn"].get_scan(scan)
    rhopedt = m_file.data["rhopedt"].get_scan(scan)
    tbeta = m_file.data["tbeta"].get_scan(scan)
    teped = m_file.data["teped"].get_scan(scan)
    tesep = m_file.data["tesep"].get_scan(scan)
    alphan = m_file.data["alphan"].get_scan(scan)
    alphat = m_file.data["alphat"].get_scan(scan)
    ne0 = m_file.data["ne0"].get_scan(scan)
    te0 = m_file.data["te0"].get_scan(scan)

    # Plasma
    global triang
    global alphaj
    global q0
    global q95
    global kallenbach_switch

    triang = m_file.data["triang95"].get_scan(scan)
    alphaj = m_file.data["alphaj"].get_scan(scan)
    q0 = m_file.data["q0"].get_scan(scan)
    q95 = m_file.data["q95"].get_scan(scan)
    kallenbach_switch = m_file.data["kallenbach_switch"].get_scan(scan)

    # Radial position  -- 0
    # Electron density -- 1
    # Electron temperature -- 2
    # Ion temperature -- 3
    # Deuterium density -- 4
    # Tritium density -- 5
    # BS current density(MA/m^2) -- 6
    # CD current dens(MA/m^2) -- 7
    # Total current dens(MA/m^2) -- 8
    # Poloidal current(R*Bp)(T.m) -- 9
    # Safety factor q -- 10
    # Volume (m^3) -- 11
    # dVolume/dr (m^2) -- 12
    # Plasma conductivity(MA/(V.m) -- 13
    # Alpha press(keV*10^10 m^-3) -- 14
    # Ion dens(10^19 m^-3) -- 15
    # Poloidal flux (Wb) -- 16
    # rad profile
    global ssync
    global bt
    global vol
    ssync = m_file.data["ssync"].get_scan(scan)
    bt = m_file.data["bt"].get_scan(scan)
    vol = m_file.data["vol"].get_scan(scan)

    # Build the dictionaries of radial and vertical build values and cumulative values
    radial = {}
    cumulative_radial = {}
    subtotal = 0
    for item in RADIAL_BUILD:
        if item == "rminori" or item == "rminoro":
            build = m_file.data["rminor"].get_scan(scan)
        elif item == "vvblgapi" or item == "vvblgapo":
            build = m_file.data["vvblgap"].get_scan(scan)
        elif "d_vv_in" in item:
            build = m_file.data["d_vv_in"].get_scan(scan)
        elif "d_vv_out" in item:
            build = m_file.data["d_vv_out"].get_scan(scan)
        else:
            build = m_file.data[item].get_scan(scan)

    radial[item] = build
    subtotal += build
    cumulative_radial[item] = subtotal

    global upper
    global cumulative_upper
    upper = {}
    cumulative_upper = {}
    subtotal = 0
    for item in vertical_upper:
        upper[item] = m_file.data[item].get_scan(scan)
        subtotal += upper[item]
        cumulative_upper[item] = subtotal

    global lower
    global cumulative_lower
    lower = {}
    cumulative_lower = {}
    subtotal = 0
    for item in vertical_lower:
        lower[item] = m_file.data[item].get_scan(scan)
        subtotal -= lower[item]
        cumulative_lower[item] = subtotal

    colour_dict = {}
    colour_dict["ohcth"] = solenoid
    colour_dict["tfcth"] = tfc
    colour_dict["thshield_ib"] = thermal_shield
    colour_dict["thshield_ob"] = thermal_shield
    colour_dict["thshield_vb"] = thermal_shield
    colour_dict["d_vv_in"] = vessel
    colour_dict["d_vv_out"] = vessel
    colour_dict["d_vv_top"] = vessel
    colour_dict["d_vv_bot"] = vessel
    colour_dict["shldith"] = shield
    colour_dict["blnkith"] = blanket
    colour_dict["rminor"] = plasma
    colour_dict["fwith"] = firstwall
    colour_dict["fwoth"] = firstwall

    # read MFILE
    # m_file = mf.MFile(args.f)
    # scan = scan

    # create main plot
    page1 = plt.figure(figsize=(12, 9), dpi=80)
    page2 = plt.figure(figsize=(12, 9), dpi=80)

    # run main_plot
    main_plot(page1, page2, m_file, scan=scan, demo_ranges=demo_ranges)

    # with bpdf.PdfPages(args.o) as pdf:
    with bpdf.PdfPages(args.f + "SUMMARY.pdf") as pdf:
        pdf.savefig(page1)
        pdf.savefig(page2)

    # show fig if option used
    if args.show:
        plt.show(block=True)

    # This bit doesn't work - the argument is not recognised for some reason.:
    # if args.svg:
    #    save_plots(m_file)
    plt.close(page1)
    plt.close(page2)


if __name__ == "__main__":
    main()
