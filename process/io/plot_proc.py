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

import argparse
import os
from argparse import RawTextHelpFormatter
from importlib import resources

import matplotlib as mpl
import matplotlib.backends.backend_pdf as bpdf
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle
from matplotlib.path import Path

import process.io.mfile as mf
from process.geometry.blanket_geometry import (
    blanket_geometry_double_null,
    blanket_geometry_single_null,
)
from process.geometry.cryostat_geometry import cryostat_geometry
from process.geometry.firstwall_geometry import (
    first_wall_geometry_double_null,
    first_wall_geometry_single_null,
)
from process.geometry.pfcoil_geometry import pfcoil_geometry
from process.geometry.plasma_geometry import plasma_geometry
from process.geometry.shield_geometry import (
    shield_geometry_double_null,
    shield_geometry_single_null,
)
from process.geometry.tfcoil_geometry import (
    tfcoil_geometry_d_shape,
    tfcoil_geometry_rectangular_shape,
)
from process.geometry.vacuum_vessel_geometry import (
    vacuum_vessel_geometry_double_null,
    vacuum_vessel_geometry_single_null,
)
from process.impurity_radiation import read_impurity_file
from process.io.python_fortran_dicts import get_dicts

if os.name == "posix" and "DISPLAY" not in os.environ:
    mpl.use("Agg")
mpl.rcParams["figure.max_open_warning"] = 40


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Produces a three page summary of the PROCESS MFILE output, using the MFILE.  "
        "For info please see https://github.com/ukaea/PROCESS?tab=readme-ov-file#contacts ",
        formatter_class=RawTextHelpFormatter,
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

    parser.add_argument(
        "-c",
        "--colour",
        type=int,
        help=(
            "Which colour scheme to use for cross-section plots\n"
            "1: Original PROCESS (default)\n"
            "2: BLUEMIRA"
        ),
        default=1,
    )

    return parser.parse_args(args)


# Colours are PROCESS defualt, BLUEMIRA
SOLENOID_COLOUR = ["pink", "#1764ab"]
CSCOMPRESSION_COLOUR = ["maroon", "#33CCCC"]
TFC_COLOUR = ["cyan", "#084a91"]
THERMAL_SHIELD_COLOUR = ["gray", "#e3eef9"]
VESSEL_COLOUR = ["green", "#b7d4ea"]
SHIELD_COLOUR = ["green", "#94c4df"]
BLANKET_COLOUR = ["magenta", "#4a98c9"]
PLASMA_COLOUR = ["khaki", "#cc8acc"]
CRYOSTAT_COLOUR = ["red", "#2e7ebc"]
FIRSTWALL_COLOUR = ["darkblue", "darkblue"]
NBSHIELD_COLOUR = ["black", "black"]

thin = 0.0

RADIAL_BUILD = [
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
    "vvblgapi",
    "dr_blkt_inboard",
    "dr_fw_inboard",
    "dr_fw_plasma_gap_inboard",
    "rminori",
    "rminoro",
    "dr_fw_plasma_gap_outboard",
    "dr_fw_outboard",
    "dr_blkt_outboard",
    "vvblgapo",
    "dr_shld_outboard",
    "dr_vv_outboard",
    "dr_shld_vv_gap_outboard",
    "dr_shld_thermal_outboard",
    "dr_tf_shld_gap",
    "dr_tf_outboard",
]

vertical_lower = [
    "rminor*kappa",
    "vgap_xpoint_divertor",
    "divfix",
    "shldlth",
    "d_vv_bot",
    "vgap_vv_thermalshield",
    "thshield_vb",
    "dr_tf_shld_gap",
    "dr_tf_inboard",
]

ANIMATION_INFO = [
    ("rmajor", "Major radius", "m"),
    ("rminor", "Minor radius", "m"),
    ("aspect", "Aspect ratio", ""),
]

rtangle = np.pi / 2
rtangle2 = 2 * rtangle


def plot_plasma(axis, mfile_data, scan, colour_scheme):
    """Plots the plasma boundary arcs.

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
        colour_scheme --> colour scheme to use for plots

    """

    r_0 = mfile_data.data["rmajor"].get_scan(scan)
    a = mfile_data.data["rminor"].get_scan(scan)
    triang = mfile_data.data["triang"].get_scan(scan)
    kappa = mfile_data.data["kappa"].get_scan(scan)
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    i_plasma_shape = mfile_data.data["i_plasma_shape"].get_scan(scan)
    plasma_square = mfile_data.data["plasma_square"].get_scan(scan)

    pg = plasma_geometry(
        rmajor=r_0,
        rminor=a,
        triang=triang,
        kappa=kappa,
        i_single_null=i_single_null,
        i_plasma_shape=i_plasma_shape,
        square=plasma_square,
    )
    if i_plasma_shape == 0:
        # Plot the 2 plasma outline arcs.
        axis.plot(pg.rs[0], pg.zs[0], color="black")
        axis.plot(pg.rs[1], pg.zs[1], color="black")

        # Set triang_95 to stop plotting plasma past boundary
        # Assume IPDG scaling
        triang_95 = triang / 1.5

        # Colour in right side of plasma
        axis.fill_between(
            x=pg.rs[0],
            y1=pg.zs[0],
            where=(pg.rs[0] > r_0 - (triang_95 * a * 1.5)),
            color=PLASMA_COLOUR[colour_scheme - 1],
        )
        # Colour in left side of plasma
        axis.fill_between(
            x=pg.rs[1],
            y1=pg.zs[1],
            where=(pg.rs[1] < r_0 - (triang_95 * a * 1.5)),
            color=PLASMA_COLOUR[colour_scheme - 1],
        )

    elif i_plasma_shape == 1:
        axis.plot(pg.rs, pg.zs, color="black")
        axis.fill(pg.rs, pg.zs, color=PLASMA_COLOUR[colour_scheme - 1])


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
            cumulative_build += mfile_data.data["dr_shld_blkt_gap"].get_scan(scan)
        elif "dr_vv_inboard" in item:
            cumulative_build += mfile_data.data["dr_vv_inboard"].get_scan(scan)
        elif "dr_vv_outboard" in item:
            cumulative_build += mfile_data.data["dr_vv_outboard"].get_scan(scan)
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
            build = mfile_data.data["dr_shld_blkt_gap"].get_scan(scan)
        elif "dr_vv_inboard" in item:
            build = mfile_data.data["dr_vv_inboard"].get_scan(scan)
        elif "dr_vv_outboard" in item:
            build = mfile_data.data["dr_vv_outboard"].get_scan(scan)
        else:
            build = mfile_data.data[item].get_scan(scan)
        cumulative_build += build
        if item == section:
            break
    previous = cumulative_build - build
    return (cumulative_build, previous)


def poloidal_cross_section(axis, mfile_data, scan, demo_ranges, colour_scheme):
    """Function to plot poloidal cross-section

    Arguments:
      axis --> axis object to add plot to
      mfile_data --> MFILE data object
      scan --> scan number to use
      colour_scheme --> colour scheme to use for plots

    """

    axis.set_xlabel("R / m")
    axis.set_ylabel("Z / m")
    axis.set_title("Poloidal cross-section")
    axis.minorticks_on()

    plot_vacuum_vessel(axis, mfile_data, scan, colour_scheme)
    plot_shield(axis, mfile_data, scan, colour_scheme)
    plot_blanket(axis, mfile_data, scan, colour_scheme)
    plot_firstwall(axis, mfile_data, scan, colour_scheme)

    plot_plasma(axis, mfile_data, scan, colour_scheme)
    plot_centre_cross(axis, mfile_data, scan)
    plot_cryostat(axis, mfile_data, scan, colour_scheme)

    plot_tf_coils(axis, mfile_data, scan, colour_scheme)
    plot_pf_coils(axis, mfile_data, scan, colour_scheme)

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


def plot_current_profiles_over_time(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int
) -> None:
    """
    Plots the current profiles over time for PF circuits, CS coil, and plasma.

    Arguments:
        axis (plt.Axes): Axis object to plot to.
        mfile_data (mf.MFile): MFILE data object.
        scan (int): Scan number to use.
    """
    t_precharge = mfile_data.data["t_precharge"].get_scan(scan)
    t_current_ramp_up = mfile_data.data["t_current_ramp_up"].get_scan(scan)
    t_fusion_ramp = mfile_data.data["t_fusion_ramp"].get_scan(scan)
    t_burn = mfile_data.data["t_burn"].get_scan(scan)
    t_ramp_down = mfile_data.data["t_ramp_down"].get_scan(scan)

    # Define a cumulative sum list for each point in the pulse
    t_steps = np.cumsum([
        0,
        t_precharge,
        t_current_ramp_up,
        t_fusion_ramp,
        t_burn,
        t_ramp_down,
    ])

    # Find the number of PF circuits, ncirt includes the CS and plasma circuits
    ncirt = mfile_data.data["ncirt"].get_scan(scan)

    # Extract PF circuit times
    # ncirt contains the CS and plasma at the end so we subtract 2
    pf_circuits = {}
    for i in range(int(ncirt - 2)):
        pf_circuits[f"PF Circuit {i}"] = [
            mfile_data.data[f"pfc{i}t{j}"].get_scan(scan) for j in range(6)
        ]
        # Change from 0 to 1 index to align with poloidal cross-section plot numbering
        axis.plot(
            t_steps,
            pf_circuits[f"PF Circuit {i}"],
            label=f"PF Coil {i + 1}",
            linestyle="--",
        )

    # Since CS may not always be present try to retireve values
    try:
        cs_circuit = [mfile_data.data[f"cs_t{i}"].get_scan(scan) for i in range(6)]
        axis.plot(t_steps, cs_circuit, label="CS Coil", linestyle="--")
    except KeyError:
        pass

    # Plasma current values
    plasmat1 = mfile_data.data["plasmat1"].get_scan(scan)
    plasmat2 = mfile_data.data["plasmat2"].get_scan(scan)
    plasmat3 = mfile_data.data["plasmat3"].get_scan(scan)
    plasmat4 = mfile_data.data["plasmat4"].get_scan(scan)
    plasmat5 = mfile_data.data["plasmat5"].get_scan(scan)

    # x-coirdinates for the plasma current
    x_plasma = t_steps[1:]
    # x-coirdinates for the plasma current
    y_plasma = [plasmat1, plasmat2, plasmat3, plasmat4, plasmat5]

    # Plot the plasma current
    axis.plot(x_plasma, y_plasma, "black", linewidth=2, label="Plasma")

    # Move the x-axis to 0 on the y-axis
    axis.spines["bottom"].set_position("zero")

    # Annotate key points
    # Create a secondary x-axis for annotations
    secax = axis.secondary_xaxis("bottom")
    secax.set_xticks(t_steps)
    secax.set_xticklabels(
        [
            "Precharge",
            r"$I_{\text{P}}$ Ramp-Up",
            "Fusion Ramp",
            "Burn",
            "Ramp Down",
            "Between Pulse",
        ],
        rotation=60,
    )
    secax.tick_params(axis="x", which="major")

    # Add axis labels
    axis.set_xlabel("Time [s]", fontsize=12)
    axis.xaxis.set_label_coords(1.05, 0.5)
    axis.set_ylabel("Current [A]", fontsize=12)

    # Add a title
    axis.set_title("Current Profiles Over Time", fontsize=14)

    # Add a legend
    axis.legend()

    axis.set_yscale("symlog")

    # Add a grid for better readability
    axis.grid(True, linestyle="--", alpha=0.6)


def plot_cryostat(axis, _mfile_data, _scan, colour_scheme):
    """Function to plot cryostat in poloidal cross-section"""

    rects = cryostat_geometry(
        r_cryostat_inboard=r_cryostat_inboard,
        dr_cryostat=dr_cryostat,
        z_cryostat_half_inside=z_cryostat_half_inside,
    )

    for rec in rects:
        axis.add_patch(
            patches.Rectangle(
                xy=(rec.anchor_x, rec.anchor_z),
                width=rec.width,
                height=rec.height,
                facecolor=CRYOSTAT_COLOUR[colour_scheme - 1],
            )
        )


def color_key(axis, mfile_data, scan, colour_scheme):
    """Function to plot the colour key
    Arguments:
      axis --> object to add plot to
      colour_scheme --> colour scheme to use for plots
    """

    axis.set_ylim([0, 10])
    axis.set_xlim([0, 10])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    labels = [
        ("CS coil", SOLENOID_COLOUR[colour_scheme - 1]),
        ("CS comp", CSCOMPRESSION_COLOUR[colour_scheme - 1]),
        ("TF coil", TFC_COLOUR[colour_scheme - 1]),
        ("Thermal shield", THERMAL_SHIELD_COLOUR[colour_scheme - 1]),
        ("VV & shield", VESSEL_COLOUR[colour_scheme - 1]),
        ("Blanket", BLANKET_COLOUR[colour_scheme - 1]),
        ("First wall", FIRSTWALL_COLOUR[colour_scheme - 1]),
        ("Plasma", PLASMA_COLOUR[colour_scheme - 1]),
        ("PF coils", "none"),
    ]

    if (mfile_data.data["iefrf"].get_scan(scan) in [5, 8]) or (
        mfile_data.data["iefrffix"].get_scan(scan) in [5, 8]
    ):
        labels.append(("NB duct shield", NBSHIELD_COLOUR[colour_scheme - 1]))
        labels.append(("Cryostat", CRYOSTAT_COLOUR[colour_scheme - 1]))
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
                [x_pos + 1.5, y_pos - 0.35],
                0.5,
                0.4,
                lw=0 if color != "none" else 1,
                facecolor=color if color != "none" else "none",
                edgecolor="black" if color == "none" else "none",
            )
        )


def toroidal_cross_section(axis, mfile_data, scan, demo_ranges, colour_scheme):
    """Function to plot toroidal cross-section
    Arguments:
      axis --> axis object to add plot to
      mfile_data --> MFILE data object
      scan --> scan number to use
      colour_scheme --> colour scheme to use for plots
    """

    axis.set_xlabel("x / m")
    axis.set_ylabel("y / m")
    axis.set_title("Toroidal cross-section")
    axis.minorticks_on()

    arc(axis, rmajor, style="dashed")

    # Colour in the main components
    r2, r1 = cumulative_radial_build2("dr_cs", mfile_data, scan)
    arc_fill(axis, r1, r2, color=SOLENOID_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_cs_precomp", mfile_data, scan)
    arc_fill(axis, r1, r2, color=CSCOMPRESSION_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_tf_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=TFC_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_shld_thermal_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=THERMAL_SHIELD_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_vv_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=VESSEL_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_shld_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=VESSEL_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_blkt_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=BLANKET_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_fw_inboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=FIRSTWALL_COLOUR[colour_scheme - 1])

    arc_fill(
        axis, rmajor - rminor, rmajor + rminor, color=PLASMA_COLOUR[colour_scheme - 1]
    )

    r2, r1 = cumulative_radial_build2("dr_fw_outboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=FIRSTWALL_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_blkt_outboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=BLANKET_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_shld_outboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=SHIELD_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_vv_outboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=VESSEL_COLOUR[colour_scheme - 1])

    r2, r1 = cumulative_radial_build2("dr_shld_thermal_outboard", mfile_data, scan)
    arc_fill(axis, r1, r2, color=THERMAL_SHIELD_COLOUR[colour_scheme - 1])

    arc_fill(
        axis,
        r_cryostat_inboard,
        r_cryostat_inboard + dr_cryostat,
        color=CRYOSTAT_COLOUR[colour_scheme - 1],
    )

    # Segment the TF coil inboard
    # Calculate centrelines
    n = int(n_tf_coils / 4) + 1
    spacing = 2 * np.pi / n_tf_coils
    i = np.arange(0, n)

    ang = i * spacing
    angl = ang - spacing / 2
    angu = ang + spacing / 2
    r1, null = cumulative_radial_build2("dr_cs_tf_gap", mfile_data, scan)
    r2, null = cumulative_radial_build2("dr_tf_inboard", mfile_data, scan)
    r4, r3 = cumulative_radial_build2("dr_tf_outboard", mfile_data, scan)

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

    for item in i:
        # Neutral beam shielding
        TF_outboard(
            axis,
            item,
            n_tf_coils=n_tf_coils,
            r3=r3,
            r4=r4,
            w=w + nbshield,
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
            facecolor=TFC_COLOUR[colour_scheme - 1],
        )

    iefrf = mfile_data.data["iefrf"].get_scan(scan)
    if (iefrf == 5) or (iefrf == 8):
        # Neutral beam geometry
        a = w
        b = dr_tf_outboard
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


def TF_outboard(axis, item, n_tf_coils, r3, r4, w, facecolor):
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
    verts = list(zip(xs1, ys1, strict=False))
    verts.extend(list(zip(xs2, ys2, strict=False)))
    endpoint = [(r2, 0)]
    verts.extend(endpoint)
    path = Path(verts, closed=True)
    patch = patches.PathPatch(path, facecolor=color, lw=0)
    axis.add_patch(patch)


def plot_neprofile(prof, demo_ranges):
    """Function to plot density profile
    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel(r"$n $ $[10^{19} \mathrm{m}^{-3}]$")
    prof.set_title("Density profile")

    if ipedestal == 1:
        rhocore = np.linspace(0, rhopedn)
        necore = neped + (ne0 - neped) * (1 - rhocore**2 / rhopedn**2) ** alphan
        nicore = necore * (nd_fuel_ions / dene)

        rhosep = np.linspace(rhopedn, 1)
        neesep = nesep + (neped - nesep) * (1 - rhosep) / (1 - min(0.9999, rhopedn))
        nisep = neesep * (nd_fuel_ions / dene)

        rho = np.append(rhocore, rhosep)
        ne = np.append(necore, neesep)
        ni = np.append(nicore, nisep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        ne = ne0 * (1 - rho**2) ** alphan
    ne = ne / 1e19
    ni = ni / 1e19
    prof.plot(rho, ni, label=r"$n_{\text{i,fuel}}$", color="red")
    prof.plot(rho, ne, label="$n_{e}$", color="blue")
    prof.legend()

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    prof.set_xlim([0, 1])
    if demo_ranges:
        prof.set_ylim([0, 20])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])

    if ipedestal != 0:
        # Print pedestal lines
        prof.axhline(
            y=neped / 1e19,
            xmax=rhopedn,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        prof.vlines(
            x=rhopedn,
            ymin=0.0,
            ymax=neped / 1e19,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        prof.minorticks_on()

        # Add text box with density profile parameters
        textstr_density = "\n".join((
            rf"$n_{{\text{{e,0}}}}$: {ne0:.3e} m$^{{-3}}$"
            rf"$\hspace{{4}} \alpha_{{\text{{n}}}}$: {alphan:.3f}",
            rf"$n_{{\text{{e,ped}}}}$: {neped:.3e} m$^{{-3}}$"
            r"$ \hspace{3} \frac{\langle n_i \rangle}{\langle n_e \rangle}$: "
            f"{nd_fuel_ions / dene:.3f}",
            rf"$f_{{\text{{GW e,ped}}}}$: {fgwped_out:.3f}",
            rf"$\rho_{{\text{{ped,n}}}}$: {rhopedn:.3f}",
            rf"$n_{{\text{{e,sep}}}}$: {nesep:.3e} m$^{{-3}}$",
            rf"$f_{{\text{{GW e,sep}}}}$: {fgwsep_out:.3f}",
        ))

        props_density = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
        prof.text(
            0.0,
            -0.16,
            textstr_density,
            transform=prof.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=props_density,
        )

    # ---


def plot_jprofile(prof):
    """Function to plot density profile
    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel(r"Current density $[kA/m^2]$")
    prof.set_title("$J$ profile")
    prof.minorticks_on()
    prof.set_xlim([0, 1.0])

    rho = np.linspace(0, 1)
    y2 = (j_plasma_0 * (1 - rho**2) ** alphaj) / 1e3

    prof.plot(rho, y2, label="$n_{i}$", color="red")

    textstr_j = "\n".join((
        r"$j_0$: " + f"{y2[0]:.3f} kA m$^{{-2}}$\n",
        r"$\alpha_J$: " + f"{alphaj:.3f}",
    ))

    props_j = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    prof.text(
        1.1,
        0.75,
        textstr_j,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_j,
    )

    prof.text(
        0.05,
        0.04,
        "*Current profile is assumed to be parabolic",
        fontsize=10,
        ha="left",
        transform=plt.gcf().transFigure,
    )


def plot_teprofile(prof, demo_ranges):
    """Function to plot temperature profile
    Arguments:
      prof --> axis object to add plot to
    """

    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel("$T$ [keV]")
    prof.set_title("Temperature profile")

    if ipedestal == 1:
        rhocore = np.linspace(0.0, rhopedt)
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
    prof.plot(rho, te, color="blue", label="$T_{e}$")
    prof.plot(rho, te[:] * tratio, color="red", label="$T_{i}$")
    prof.legend()

    # Ranges
    # ---
    prof.set_xlim([0, 1])
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim([0, 50])

    # Adapatative ranges
    else:
        prof.set_ylim([0, prof.get_ylim()[1]])

    if ipedestal != 0:
        # Plot pedestal lines
        prof.axhline(
            y=teped, xmax=rhopedt, color="r", linestyle="-", linewidth=0.4, alpha=0.4
        )
        prof.vlines(
            x=rhopedt,
            ymin=0.0,
            ymax=teped,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        prof.minorticks_on()

    # Add text box with temperature profile parameters
    textstr_temperature = "\n".join((
        rf"$T_{{\text{{e,0}}}}$: {te0:.3f} keV"
        rf"$\hspace{{4}} \alpha_{{\text{{T}}}}$: {alphat:.3f}",
        rf"$T_{{\text{{e,ped}}}}$: {teped:.3f} keV"
        r"$ \hspace{4} \frac{\langle T_i \rangle}{\langle T_e \rangle}$: "
        f"{tratio:.3f}",
        rf"$\rho_{{\text{{ped,T}}}}$: {rhopedt:.3f}",
        rf"$T_{{\text{{e,sep}}}}$: {tesep:.3f} keV",
    ))

    props_temperature = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    prof.text(
        0.0,
        -0.16,
        textstr_temperature,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_temperature,
    )
    # ---


def plot_qprofile(prof, demo_ranges):
    """Function to plot q profile, formula taken from Nevins bootstrap model.

    Arguments:
      prof --> axis object to add plot to
    """
    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel("$q$")
    prof.set_title("$q$ profile")
    prof.minorticks_on()

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
        prof.set_ylim([0, q95 * 1.2])

    prof.text(
        0.6,
        0.04,
        "*Profile is not calculated, only $q_0$ and $q_{95}$ are known.",
        fontsize=10,
        ha="left",
        transform=plt.gcf().transFigure,
    )
    # ---

    textstr_q = "\n".join((
        r"$q_0$: " + f"{q0:.3f}\n",
        r"$q_{95}$: " + f"{q95:.3f}",
    ))

    props_q = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    prof.text(
        -0.4,
        0.75,
        textstr_q,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_q,
    )


def read_imprad_data(_skiprows, data_path):
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
    return np.array(lzdata, dtype=float)


def plot_radprofile(prof, mfile_data, scan, impp, demo_ranges) -> float:
    """Function to plot radiation profile, formula taken from ???.

    Arguments:
      prof --> axis object to add plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use
      impp --> impurity path
    """

    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel(r"$P_{\mathrm{rad}}$ $[\mathrm{MW.m}^{-3}]$")
    prof.set_title("Line & Bremsstrahlung radiation profile")

    # read in the impurity data
    imp_data = read_imprad_data(2, impp)

    # find impurity densities
    imp_frac = np.array([
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
    ])

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
                ne[q] = neped + (ne0 - neped) * (1 - rho[q] ** 2 / rhopedn**2) ** alphan
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

    # Intailise the radiation profile arrays
    pimpden = np.zeros([imp_data.shape[0], te.shape[0]])
    lz = np.zeros([imp_data.shape[0], te.shape[0]])
    prad = np.zeros(te.shape[0])

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

        for l in range(imp_data.shape[0]):  # noqa: E741
            prad[k] = prad[k] + pimpden[l][k] * 2.0e-6

    prof.plot(rho, prad, label="Total", linestyle="dotted")
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
    prof.legend(loc="upper left", bbox_to_anchor=(-0.1, -0.1), ncol=4)
    prof.minorticks_on()

    # Ranges
    # ---
    prof.set_xlim([0, 1.0])
    prof.set_yscale("log")
    prof.yaxis.grid(True, which="both", alpha=0.2)
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim([1e-6, 0.5])

    # Adapatative ranges
    else:
        prof.set_ylim([1e-6, prof.get_ylim()[1]])
    # ---


def plot_vacuum_vessel(axis, mfile_data, scan, colour_scheme):
    """Function to plot vacuum vessel

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
        colour_scheme --> colour scheme to use for plots
    """

    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)

    # Outer side (furthest from plasma)
    radx_outer = (
        cumulative_radial_build("dr_vv_outboard", mfile_data, scan)
        + cumulative_radial_build("dr_shld_vv_gap_inboard", mfile_data, scan)
    ) / 2.0
    rminx_outer = (
        cumulative_radial_build("dr_vv_outboard", mfile_data, scan)
        - cumulative_radial_build("dr_shld_vv_gap_inboard", mfile_data, scan)
    ) / 2.0

    # Inner side (nearest to the plasma)
    radx_inner = (
        cumulative_radial_build("dr_shld_outboard", mfile_data, scan)
        + cumulative_radial_build("dr_vv_inboard", mfile_data, scan)
    ) / 2.0
    rminx_inner = (
        cumulative_radial_build("dr_shld_outboard", mfile_data, scan)
        - cumulative_radial_build("dr_vv_inboard", mfile_data, scan)
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
            color=VESSEL_COLOUR[colour_scheme - 1],
            lw=0.01,
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
            color=VESSEL_COLOUR[colour_scheme - 1],
            lw=0.01,
        )


def plot_shield(axis, mfile_data, scan, colour_scheme):
    """Function to plot shield

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
        colour_scheme --> colour scheme to use for plots
    """

    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)

    # Side furthest from plasma
    radx_far = (
        cumulative_radial_build("dr_shld_outboard", mfile_data, scan)
        + cumulative_radial_build("dr_vv_inboard", mfile_data, scan)
    ) / 2.0
    rminx_far = (
        cumulative_radial_build("dr_shld_outboard", mfile_data, scan)
        - cumulative_radial_build("dr_vv_inboard", mfile_data, scan)
    ) / 2.0

    # Side nearest to the plasma
    radx_near = (
        cumulative_radial_build("vvblgapo", mfile_data, scan)
        + cumulative_radial_build("dr_shld_inboard", mfile_data, scan)
    ) / 2.0
    rminx_near = (
        cumulative_radial_build("vvblgapo", mfile_data, scan)
        - cumulative_radial_build("dr_shld_inboard", mfile_data, scan)
    ) / 2.0

    if i_single_null == 1:
        shield_geometry = shield_geometry_single_null(
            cumulative_upper=cumulative_upper,
            radx_far=radx_far,
            rminx_far=rminx_far,
            radx_near=radx_near,
            rminx_near=rminx_near,
            triang=triang_95,
            cumulative_lower=cumulative_lower,
        )
    else:
        shield_geometry = shield_geometry_double_null(
            cumulative_lower=cumulative_lower,
            radx_far=radx_far,
            radx_near=radx_near,
            rminx_far=rminx_far,
            rminx_near=rminx_near,
            triang=triang_95,
        )

    axis.plot(shield_geometry.rs, shield_geometry.zs, color="black", lw=thin)
    axis.fill(
        shield_geometry.rs,
        shield_geometry.zs,
        color=SHIELD_COLOUR[colour_scheme - 1],
        lw=0.01,
    )


def plot_blanket(axis, mfile_data, scan, colour_scheme) -> None:
    """Function to plot blanket

    Arguments:
      axis --> axis object to plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use
      colour_scheme --> colour scheme to use for plots

    """

    # Single null: Draw top half from output
    # Double null: Reflect bottom half to top
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)
    if int(i_single_null) == 1:
        blnktth = mfile_data.data["blnktth"].get_scan(scan)
    else:
        blnktth = 0.0

    c_shldith = cumulative_radial_build("dr_shld_inboard", mfile_data, scan)
    c_blnkoth = cumulative_radial_build("dr_blkt_outboard", mfile_data, scan)

    if i_single_null == 1:
        # Upper blanket: outer surface
        radx_outer = (
            cumulative_radial_build("dr_blkt_outboard", mfile_data, scan)
            + cumulative_radial_build("vvblgapi", mfile_data, scan)
        ) / 2.0
        rminx_outer = (
            cumulative_radial_build("dr_blkt_outboard", mfile_data, scan)
            - cumulative_radial_build("vvblgapi", mfile_data, scan)
        ) / 2.0

        # Upper blanket: inner surface
        radx_inner = (
            cumulative_radial_build("dr_fw_outboard", mfile_data, scan)
            + cumulative_radial_build("dr_blkt_inboard", mfile_data, scan)
        ) / 2.0
        rminx_inner = (
            cumulative_radial_build("dr_fw_outboard", mfile_data, scan)
            - cumulative_radial_build("dr_blkt_inboard", mfile_data, scan)
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
            dr_blkt_inboard=dr_blkt_inboard,
            dr_blkt_outboard=dr_blkt_outboard,
        )

        # Plot blanket
        axis.plot(
            bg_single_null.rs,
            bg_single_null.zs,
            color="black",
            lw=thin,
        )

        axis.fill(
            bg_single_null.rs,
            bg_single_null.zs,
            color=BLANKET_COLOUR[colour_scheme - 1],
            lw=0.01,
        )

    if i_single_null == 0:
        bg_double_null = blanket_geometry_double_null(
            cumulative_lower=cumulative_lower,
            triang=triang_95,
            blnktth=blnktth,
            c_shldith=c_shldith,
            c_blnkoth=c_blnkoth,
            dr_blkt_inboard=dr_blkt_inboard,
            dr_blkt_outboard=dr_blkt_outboard,
        )
        # Plot blanket
        axis.plot(bg_double_null.rs[0], bg_double_null.zs[0], color="black", lw=thin)
        axis.fill(
            bg_double_null.rs[0],
            bg_double_null.zs[0],
            color=BLANKET_COLOUR[colour_scheme - 1],
            lw=0.01,
        )
        if dr_blkt_inboard > 0.0:
            # only plot inboard blanket if inboard blanket thickness > 0
            axis.plot(
                bg_double_null.rs[1], bg_double_null.zs[1], color="black", lw=thin
            )
            axis.fill(
                bg_double_null.rs[1],
                bg_double_null.zs[1],
                color=BLANKET_COLOUR[colour_scheme - 1],
                lw=0.01,
            )


def plot_firstwall(axis, mfile_data, scan, colour_scheme):
    """Function to plot first wall

    Arguments:
      axis --> axis object to plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use
      colour_scheme --> colour scheme to use for plots

    """

    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)
    triang_95 = mfile_data.data["triang95"].get_scan(scan)
    if int(i_single_null) == 1:
        blnktth = mfile_data.data["blnktth"].get_scan(scan)
        tfwvt = mfile_data.data["fwtth"].get_scan(scan)
    else:
        blnktth = tfwvt = 0.0

    c_blnkith = cumulative_radial_build("dr_blkt_inboard", mfile_data, scan)
    c_fwoth = cumulative_radial_build("dr_fw_outboard", mfile_data, scan)

    if i_single_null == 1:
        # Upper first wall: outer surface
        radx_outer = (
            cumulative_radial_build("dr_fw_outboard", mfile_data, scan)
            + cumulative_radial_build("dr_blkt_inboard", mfile_data, scan)
        ) / 2.0
        rminx_outer = (
            cumulative_radial_build("dr_fw_outboard", mfile_data, scan)
            - cumulative_radial_build("dr_blkt_inboard", mfile_data, scan)
        ) / 2.0

        # Upper first wall: inner surface
        radx_inner = (
            cumulative_radial_build("dr_fw_plasma_gap_outboard", mfile_data, scan)
            + cumulative_radial_build("dr_fw_inboard", mfile_data, scan)
        ) / 2.0
        rminx_inner = (
            cumulative_radial_build("dr_fw_plasma_gap_outboard", mfile_data, scan)
            - cumulative_radial_build("dr_fw_inboard", mfile_data, scan)
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
            dr_fw_inboard=dr_fw_inboard,
            dr_fw_outboard=dr_fw_outboard,
            tfwvt=tfwvt,
        )

        # Plot first wall
        axis.plot(fwg_single_null.rs, fwg_single_null.zs, color="black", lw=thin)
        axis.fill(
            fwg_single_null.rs,
            fwg_single_null.zs,
            color=FIRSTWALL_COLOUR[colour_scheme - 1],
            lw=0.01,
        )

    if i_single_null == 0:
        fwg_double_null = first_wall_geometry_double_null(
            cumulative_lower=cumulative_lower,
            triang=triang_95,
            blnktth=blnktth,
            c_blnkith=c_blnkith,
            c_fwoth=c_fwoth,
            dr_fw_inboard=dr_fw_inboard,
            dr_fw_outboard=dr_fw_outboard,
            tfwvt=tfwvt,
        )
        # Plot first wall
        axis.plot(fwg_double_null.rs[0], fwg_double_null.zs[0], color="black", lw=thin)
        axis.plot(fwg_double_null.rs[1], fwg_double_null.zs[1], color="black", lw=thin)
        axis.fill(
            fwg_double_null.rs[0],
            fwg_double_null.zs[0],
            color=FIRSTWALL_COLOUR[colour_scheme - 1],
            lw=0.01,
        )
        axis.fill(
            fwg_double_null.rs[1],
            fwg_double_null.zs[1],
            color=FIRSTWALL_COLOUR[colour_scheme - 1],
            lw=0.01,
        )


def plot_tf_coils(axis, mfile_data, scan, colour_scheme):
    """Function to plot TF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use
        colour_scheme --> colour scheme to use for plots

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
    dr_shld_thermal_inboard = mfile_data.data["dr_shld_thermal_inboard"].get_scan(scan)
    dr_shld_thermal_outboard = mfile_data.data["dr_shld_thermal_outboard"].get_scan(
        scan
    )
    dr_tf_shld_gap = mfile_data.data["dr_tf_shld_gap"].get_scan(scan)
    if y3 != 0:
        print("TF coil geometry: The value of yarc(3) is not zero, but should be.")

    if dr_shld_thermal_inboard != dr_shld_thermal_outboard:
        print(
            "dr_shld_thermal_inboard and dr_shld_thermal_outboard are different. Using dr_shld_thermal_inboard"
            "for the poloidal plot of the thermal shield."
        )

    for offset, colour in (
        (
            dr_shld_thermal_inboard + dr_tf_shld_gap,
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
        ),
        (dr_tf_shld_gap, "white"),
        (0.0, TFC_COLOUR[colour_scheme - 1]),
    ):
        # Check for TF coil shape
        if "i_tf_shape" in mfile_data.data:
            i_tf_shape = int(mfile_data.data["i_tf_shape"].get_scan(scan))
        else:
            i_tf_shape = 1

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
                dr_tf_inboard=dr_tf_inboard,
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
                path = Path(vert, closed=True)
                patch = patches.PathPatch(path, facecolor=colour, lw=0)
                axis.add_patch(patch)

        for rec in rects:
            axis.add_patch(
                patches.Rectangle(
                    xy=(rec.anchor_x, rec.anchor_z),
                    width=rec.width,
                    height=rec.height,
                    facecolor=colour,
                )
            )


def plot_tf_wp(axis, mfile_data, scan: int) -> None:
    """
    Plots inboard TF coil and winding pack.
    Author: C. Ashe

    Parameters
    ----------
    axis : matplotlib.axes object
        Axis object to plot to.
    mfile_data : MFILE data object
        Object containing data for the plot.
    scan : int
        Scan number to use.

    Returns
    -------
    None
    """

    # Import the TF variables
    r_tf_inboard_in = mfile_data.data["r_tf_inboard_in"].get_scan(scan)
    r_tf_inboard_out = mfile_data.data["r_tf_inboard_out"].get_scan(scan)
    wp_toridal_dxbig = mfile_data.data["wwp1"].get_scan(scan)
    wp_toridal_dxsmall = mfile_data.data["wwp2"].get_scan(scan)
    dr_tf_wp = mfile_data.data["dr_tf_wp"].get_scan(scan)
    side_case_dx = mfile_data.data["casths"].get_scan(scan)
    wp_inner = mfile_data.data["r_wp_inner"].get_scan(scan)
    tinstf = mfile_data.data["tinstf"].get_scan(scan)
    turns = round(mfile_data.data["n_tf_turn"].get_scan(scan))
    wp_shape = round(mfile_data.data["i_tf_wp_geom"].get_scan(scan))
    cond_type = round(mfile_data.data["i_tf_sup"].get_scan(scan))
    nose_thickness = mfile_data.data["thkcas"].get_scan(scan)
    side_thickness = mfile_data.data["casths"].get_scan(scan)
    case_plasma = mfile_data.data["i_tf_case_geom"].get_scan(scan)
    jwptf = round(mfile_data.data["jwptf"].get_scan(scan)) / 1e6
    tf_thickness = mfile_data.data["dr_tf_inboard"].get_scan(scan)
    integer_turns = mfile_data.data["i_tf_turns_integer"].get_scan(scan)

    if integer_turns == 1:
        turn_layers = mfile_data.data["n_layer"].get_scan(scan)
        turn_pancakes = mfile_data.data["n_pancake"].get_scan(scan)

    # Superconducting coil check
    if cond_type == 1:
        axis.add_patch(
            Circle(
                [0, 0],
                r_tf_inboard_in,
                facecolor="none",
                edgecolor="black",
                linestyle="--",
            ),
        )

        if case_plasma == 0:
            axis.add_patch(
                Circle(
                    [0, 0],
                    r_tf_inboard_out,
                    facecolor="none",
                    edgecolor="black",
                    linestyle="--",
                ),
            )

        # Equations for plotting the TF case
        half_case_angle = np.arctan(
            (side_case_dx + (0.5 * wp_toridal_dxbig)) / wp_inner
        )

        # X points for inboard case curve
        x11 = r_tf_inboard_in * np.cos(
            np.linspace(half_case_angle, -half_case_angle, 256, endpoint=True)
        )
        # Y points for inboard case curve
        y11 = r_tf_inboard_in * np.sin(
            np.linspace(half_case_angle, -half_case_angle, 256, endpoint=True)
        )
        # Check for plasma side case type
        if case_plasma == 0:
            # Rounded case

            # X points for outboard case curve
            x12 = r_tf_inboard_out * np.cos(
                np.linspace(half_case_angle, -half_case_angle, 256, endpoint=True)
            )

        elif case_plasma == 1:
            # Flat case

            # X points for outboard case
            x12 = np.full(256, r_tf_inboard_out)
        else:
            raise NotImplementedError("case_plasma must be 0 or 1")

        # Y points for outboard case
        y12 = r_tf_inboard_out * np.sin(
            np.linspace(half_case_angle, -half_case_angle, 256, endpoint=True)
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
        case_plasma_label = (
            f"Case: \n{nose_thickness:.4f} m nose thickness \n"
            f"{side_thickness:.4f} m sidewall thickness \n$"
            f"\\Delta$R = {tf_thickness:.4f} m \n "
        )
        if case_plasma == 0:
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(half_case_angle)),
                    (r_tf_inboard_out * np.cos(half_case_angle)),
                ],
                y13,
                color="grey",
                alpha=0.25,
                label=case_plasma_label,
            )
            # Lower main
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(half_case_angle)),
                    (r_tf_inboard_out * np.cos(half_case_angle)),
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
        elif case_plasma == 1:
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(half_case_angle)),
                    (r_tf_inboard_out),
                ],
                y13,
                color="grey",
                alpha=0.25,
                label=case_plasma_label,
            )
            # Lower main
            axis.fill_between(
                [
                    (r_tf_inboard_in * np.cos(half_case_angle)),
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

        # Plot the rectangular WP
        if wp_shape == 0:
            if integer_turns == 1:
                long_turns = round(turn_layers)
                short_turns = round(turn_pancakes)
            else:
                wp_side_ratio = (dr_tf_wp - (2 * tinstf)) / (
                    wwp1 - (2 * tinstf)
                )  # row to height
                side_unit = turns / wp_side_ratio
                root_turns = round(np.sqrt(side_unit), 1)
                long_turns = round(root_turns * wp_side_ratio)
                short_turns = round(root_turns)

            # Plots the surrounding insualtion
            axis.add_patch(
                Rectangle(
                    (wp_inner, -(0.5 * wp_toridal_dxbig)),
                    dr_tf_wp,
                    wp_toridal_dxbig,
                    color="darkgreen",
                    label=f"Insulation: \n{tinstf * 1000} mm thickness \n",
                ),
            )
            # Plots the WP inside the insulation
            axis.add_patch(
                Rectangle(
                    (wp_inner + tinstf, -(0.5 * wp_toridal_dxbig) + tinstf),
                    (dr_tf_wp - (2 * tinstf)),
                    (wp_toridal_dxbig - (2 * tinstf)),
                    color="blue",
                    label=(
                        f"Winding pack:  \n{turns} turns \n{jwptf:.4f} MA/m$^2$ \n$"
                        f"\\Delta$R= {dr_tf_wp:.4f} m \n  "
                    ),
                )
            )
            # Dvides the WP up into the turn segments
            for i in range(1, long_turns):
                axis.plot(
                    [
                        (wp_inner + tinstf) + i * (dr_tf_wp / long_turns),
                        (wp_inner + tinstf) + i * (dr_tf_wp / long_turns),
                    ],
                    [
                        -0.5 * (wp_toridal_dxbig - 2 * tinstf),
                        0.5 * (wp_toridal_dxbig - 2 * tinstf),
                    ],
                    color="white",
                    linewidth="0.25",
                    linestyle="dashed",
                )

            for i in range(1, short_turns):
                axis.plot(
                    [(wp_inner + tinstf), (wp_inner - tinstf + dr_tf_wp)],
                    [
                        (-0.5 * wp_toridal_dxbig)
                        + (i * wp_toridal_dxbig / short_turns),
                        (-0.5 * wp_toridal_dxbig)
                        + (i * wp_toridal_dxbig / short_turns),
                    ],
                    color="white",
                    linewidth="0.25",
                    linestyle="dashed",
                )

        # Plot the double rectangle winding pack
        if wp_shape == 1:
            # Inner WP insulation
            axis.add_patch(
                Rectangle(
                    (
                        wp_inner - tinstf,
                        -(0.5 * wp_toridal_dxsmall) - 2 * tinstf,
                    ),
                    (dr_tf_wp / 2) + (tinstf),
                    wp_toridal_dxsmall + (tinstf),
                    color="darkgreen",
                    label=f"Insulation: \n{tinstf * 1000} mm thickness \n",
                ),
            )

            # Outer WP insulation
            axis.add_patch(
                Rectangle(
                    (
                        wp_inner + (0.5 * dr_tf_wp) - tinstf,
                        -(0.5 * wp_toridal_dxbig) - tinstf,
                    ),
                    (dr_tf_wp / 2) + (tinstf),
                    wp_toridal_dxbig + (tinstf),
                    color="darkgreen",
                ),
            )

            # Outer WP
            axis.add_patch(
                Rectangle(
                    (wp_inner + (0.5 * dr_tf_wp), -(0.5 * wp_toridal_dxbig)),
                    (dr_tf_wp / 2) - (2 * tinstf),
                    wp_toridal_dxbig - (2 * tinstf),
                    color="blue",
                    label=(
                        f"Winding pack: \n{turns} turns \n{jwptf:.4f} MA/m$^2$ \n$"
                        f"\\Delta$R= {dr_tf_wp:.4f} m \n  "
                    ),
                ),
            )
            # Inner WP
            axis.add_patch(
                Rectangle(
                    (
                        wp_inner + tinstf,
                        -(0.5 * wp_toridal_dxsmall) - tinstf,
                    ),
                    (dr_tf_wp / 2) - (2 * tinstf),
                    wp_toridal_dxsmall - (2 * tinstf),
                    color="blue",
                ),
            )

        # Trapezium WP
        if wp_shape == 2:
            # WP insulation
            x = [wp_inner, wp_inner, (wp_inner + dr_tf_wp), (wp_inner + dr_tf_wp)]
            y = [
                (-0.5 * wp_toridal_dxsmall),
                (0.5 * wp_toridal_dxsmall),
                (0.5 * wp_toridal_dxbig),
                (-0.5 * wp_toridal_dxbig),
            ]
            axis.add_patch(
                patches.Polygon(
                    xy=list(zip(x, y, strict=False)),
                    color="darkgreen",
                    label=f"Insulation: \n{tinstf * 1000} mm thickness \n",
                )
            )

            # WP
            x = [
                wp_inner + tinstf,
                wp_inner + tinstf,
                (wp_inner + dr_tf_wp - tinstf),
                (wp_inner + dr_tf_wp - tinstf),
            ]
            y = [
                (-0.5 * wp_toridal_dxsmall + tinstf),
                (0.5 * wp_toridal_dxsmall - tinstf),
                (0.5 * wp_toridal_dxbig - tinstf),
                (-0.5 * wp_toridal_dxbig + tinstf),
            ]
            axis.add_patch(
                patches.Polygon(
                    xy=list(zip(x, y, strict=False)),
                    color="blue",
                    label=(
                        f"Winding pack: \n{turns} turns \n{jwptf:.4f} MA/m$^2$ \n"
                        f"$\\Delta$R= {dr_tf_wp:.4f} m \n  "
                    ),
                )
            )

        axis.minorticks_on()
        axis.set_xlim(0.0, r_tf_inboard_out * 1.1)
        axis.set_ylim((y14[-1] * 1.25), (-y14[-1] * 1.25))

        axis.set_title("Top-down view of inboard TF coil at midplane")
        axis.set_xlabel("Radial distance [m]")
        axis.set_ylabel("Toroidal distance [m]")
        axis.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")


def plot_tf_turn(axis, mfile_data, scan: int) -> None:
    """
    Plots inboard TF coil individual turn structure.
    Author: C. Ashe

    Parameters
    ----------
    axis : matplotlib.axes object
        Axis object to plot to.
    mfile_data : MFILE data object
        Object containing data for the plot.
    scan : int
        Scan number to use.

    Returns
    -------
    None
    """

    # Import the TF turn variables then multiply into mm
    integer_turns = mfile_data.data["i_tf_turns_integer"].get_scan(scan)
    # If integer turns switch is on then the turns can have non square dimensions
    if integer_turns == 1:
        turn_width = round(mfile_data.data["t_turn_radial"].get_scan(scan) * 1e3, 5)
        turn_height = round(mfile_data.data["t_turn_toroidal"].get_scan(scan) * 1e3, 5)
        cable_space_width_radial = round(
            mfile_data.data["t_cable_radial"].get_scan(scan) * 1e3, 5
        )
        cable_space_width_toroidal = round(
            mfile_data.data["t_cable_toroidal"].get_scan(scan) * 1e3, 5
        )
    elif integer_turns == 0:
        turn_width = round(mfile_data.data["t_turn_tf"].get_scan(scan) * 1e3, 5)
        cable_space_width = round(mfile_data.data["t_cable"].get_scan(scan) * 1e3, 5)

    he_pipe_diameter = round(mfile_data.data["dhecoil"].get_scan(scan) * 1e3, 5)
    steel_thickness = round(mfile_data.data["thwcndut"].get_scan(scan) * 1e3, 5)
    insulation_thickness = round(mfile_data.data["thicndut"].get_scan(scan) * 1e3, 5)
    internal_cable_space = round(mfile_data.data["acstf"].get_scan(scan) * 1e6, 5)
    cpttf = mfile_data.data["cpttf"].get_scan(scan)

    # Plot the total turn shape
    if integer_turns == 0:
        axis.add_patch(
            Rectangle(
                [0, 0],
                turn_width,
                turn_width,
                facecolor="red",
                label=f"Inter-turn insulation: \n{insulation_thickness} mm thickness",
                edgecolor="black",
            ),
        )
        # Plot the steel conduit
        axis.add_patch(
            Rectangle(
                [insulation_thickness, insulation_thickness],
                (turn_width - 2 * insulation_thickness),
                (turn_width - 2 * insulation_thickness),
                facecolor="grey",
                label=f"Steel Conduit: \n{steel_thickness} mm thickness",
                edgecolor="black",
            ),
        )

        # Plot the cable space
        axis.add_patch(
            Rectangle(
                [
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ],
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                facecolor="royalblue",
                label=f"Cable space: \n{cable_space_width} mm width \n{internal_cable_space} mm$^2$",
                edgecolor="black",
            ),
        )
        axis.add_patch(
            Circle(
                [(turn_width / 2), (turn_width / 2)],
                he_pipe_diameter / 2,
                facecolor="white",
                label=f"Cooling pipe: \n{he_pipe_diameter} mm diameter \n \n Current per turn: {cpttf:.2f} A",
                edgecolor="black",
            ),
        )
        axis.set_xlim(-turn_width * 0.05, turn_width * 1.05)
        axis.set_ylim(-turn_width * 0.05, turn_width * 1.05)

    # Non square turns
    elif integer_turns == 1:
        axis.add_patch(
            Rectangle(
                [0, 0],
                turn_width,
                turn_height,
                facecolor="red",
                label=f"Inter-turn insulation: \n{insulation_thickness} mm thickness",
                edgecolor="black",
            ),
        )

        # Plot the steel conduit
        axis.add_patch(
            Rectangle(
                [insulation_thickness, insulation_thickness],
                (turn_width - 2 * insulation_thickness),
                (turn_height - 2 * insulation_thickness),
                facecolor="grey",
                label=f"Steel Conduit: \n{steel_thickness} mm thickness",
                edgecolor="black",
            ),
        )

        # Plot the cable space
        axis.add_patch(
            Rectangle(
                [
                    insulation_thickness + steel_thickness,
                    insulation_thickness + steel_thickness,
                ],
                (turn_width - 2 * (insulation_thickness + steel_thickness)),
                (turn_height - 2 * (insulation_thickness + steel_thickness)),
                facecolor="royalblue",
                label=(
                    f"Cable space: \n{cable_space_width_radial} mm radial width \n"
                    f"{cable_space_width_toroidal} mm toroidal width \n{internal_cable_space} mm$^2$"
                ),
                edgecolor="black",
            ),
        )
        axis.add_patch(
            Circle(
                [(turn_width / 2), (turn_height / 2)],
                he_pipe_diameter / 2,
                facecolor="white",
                label=f"Cooling pipe: \n{he_pipe_diameter} mm diameter \n \n Current per turn: {cpttf:.2f} A",
                edgecolor="black",
            ),
        )

        axis.set_xlim(-turn_width * 0.05, turn_width * 1.05)
        axis.set_ylim(-turn_height * 0.05, turn_height * 1.05)

    axis.minorticks_on()
    axis.set_title("WP Turn Structure")
    axis.set_xlabel("X [mm]")
    axis.set_ylabel("Y [mm]")
    axis.legend(loc="upper right", bbox_to_anchor=(2.0, 1.0))


def plot_pf_coils(axis, mfile_data, scan, colour_scheme):
    """Function to plot PF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use
        colour_scheme --> colour scheme to use for plots
    """

    coils_r = []
    coils_z = []
    coils_dr = []
    coils_dz = []
    coil_text = []

    dr_bore = mfile_data.data["dr_bore"].get_scan(scan)
    dr_cs = mfile_data.data["dr_cs"].get_scan(scan)
    ohdz = mfile_data.data["ohdz"].get_scan(scan)

    # Number of coils, both PF and CS
    number_of_coils = 0
    for item in mfile_data.data:
        if "rpf[" in item:
            number_of_coils += 1

    # Check for Central Solenoid
    iohcl = mfile_data.data["iohcl"].get_scan(scan) if "iohcl" in mfile_data.data else 1

    # If Central Solenoid present, ignore last entry in for loop
    # The last entry will be the OH coil in this case
    noc = number_of_coils - 1 if iohcl == 1 else number_of_coils

    for coil in range(noc):
        coils_r.append(mfile_data.data[f"rpf[{coil:01}]"].get_scan(scan))
        coils_z.append(mfile_data.data[f"zpf[{coil:01}]"].get_scan(scan))
        coils_dr.append(mfile_data.data[f"pfdr({coil:01})"].get_scan(scan))
        coils_dz.append(mfile_data.data[f"pfdz({coil:01})"].get_scan(scan))
        coil_text.append(str(coil + 1))

    r_points, z_points, central_coil = pfcoil_geometry(
        coils_r=coils_r,
        coils_z=coils_z,
        coils_dr=coils_dr,
        coils_dz=coils_dz,
        dr_bore=dr_bore,
        dr_cs=dr_cs,
        ohdz=ohdz,
    )

    # Plot CS compression structure
    r_precomp_outer, r_precomp_inner = cumulative_radial_build2(
        "dr_cs_precomp", mfile_data, scan
    )
    axis.add_patch(
        patches.Rectangle(
            xy=(r_precomp_inner, central_coil.anchor_z),
            width=r_precomp_outer - r_precomp_inner,
            height=central_coil.height,
            facecolor=CSCOMPRESSION_COLOUR[colour_scheme - 1],
        )
    )

    for i in range(len(coils_r)):
        axis.plot(r_points[i], z_points[i], color="black")
        axis.text(
            coils_r[i],
            coils_z[i] - 0.05,
            coil_text[i],
            ha="center",
            va="center",
            fontsize=8.5 * abs(coils_dr[i] * coils_dz[i]),
        )
    axis.add_patch(
        patches.Rectangle(
            xy=(central_coil.anchor_x, central_coil.anchor_z),
            width=central_coil.width,
            height=central_coil.height,
            facecolor=SOLENOID_COLOUR[colour_scheme - 1],
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
    eqpos = 0.75
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
                axis.text(-0.05, -i, f"{data[i][0][1:]}\n", ha="left", va="center")
            elif data[i][0][0] == "!":
                value = data[i][0][1:].replace('"', "")
                axis.text(
                    0.4,
                    -i,
                    f"-->  {value} {data[i][2]}",
                    ha="left",
                    va="center",
                )
            else:
                if mfile_data.data[data[i][0]].exists:
                    dat = mfile_data.data[data[i][0]].get_scan(scan)
                    if isinstance(dat, str):
                        value = dat
                    else:
                        value = f"{mfile_data.data[data[i][0]].get_scan(scan):.4g}"
                    if "alpha" in data[i][0]:
                        value = str(float(value) + 1.0)
                    axis.text(
                        eqpos,
                        -i,
                        f"= {value} {data[i][2]}",
                        color=colorflag,
                        ha="left",
                        va="center",
                    )
                else:
                    mfile_data.data[data[i][0]].get_scan(-1)
                    axis.text(
                        eqpos,
                        -i,
                        "= ERROR! Var missing",
                        color=colorflag,
                        ha="left",
                        va="center",
                    )
        else:
            dat = data[i][0]
            value = dat if isinstance(dat, str) else f"{data[i][0]:.4g}"
            axis.text(
                eqpos,
                -i,
                f"= {value} {data[i][2]}",
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
        (f"!{mfile_data.data['runtitle'].get_scan(-1)}", "Run title", ""),
        (f"!{mfile_data.data['procver'].get_scan(-1)}", "PROCESS Version", ""),
        (f"!{mfile_data.data['date'].get_scan(-1)}", "Date:", ""),
        (f"!{mfile_data.data['time'].get_scan(-1)}", "Time:", ""),
        (f"!{mfile_data.data['username'].get_scan(-1)}", "User:", ""),
        (
            f"!{dicts['DICT_OPTIMISATION_VARS'][str(abs(int(mfile_data.data['minmax'].get_scan(-1))))]}",
            "Optimising:",
            "",
        ),
    ]

    axis.text(-0.05, 4.0, "Colour Legend:", ha="left", va="center")
    axis.text(
        0.0, 3.0, "ITR --> Iteration variable", color="red", ha="left", va="center"
    )
    axis.text(0.0, 2.0, "OP  --> Output variable", color="blue", ha="left", va="center")

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

    in_blanket_thk = mfile_data.data["dr_shld_inboard"].get_scan(
        scan
    ) + mfile_data.data["dr_blkt_inboard"].get_scan(scan)
    out_blanket_thk = mfile_data.data["dr_shld_outboard"].get_scan(
        scan
    ) + mfile_data.data["dr_blkt_outboard"].get_scan(scan)

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

    nd_impurities = mfile_data.data["nd_impurities"].get_scan(scan) / mfile_data.data[
        "dene"
    ].get_scan(scan)

    tepeak = mfile_data.data["te0"].get_scan(scan) / mfile_data.data["te"].get_scan(
        scan
    )

    nepeak = mfile_data.data["ne0"].get_scan(scan) / mfile_data.data["dene"].get_scan(
        scan
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    if "plhthresh" in mfile_data.data:
        pthresh = mfile_data.data["plhthresh"].get_scan(scan)
    else:
        pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)

    data = [
        ("fusion_power", "Fusion power", "MW"),
        ("bigq", "$Q_{p}$", ""),
        ("plasma_current_ma", "$I_p$", "MA"),
        ("bt", "Vacuum $B_T$ at $R_0$", "T"),
        ("q95", r"$q_{\mathrm{95}}$", ""),
        ("beta_norm_thermal", r"$\beta_N$, thermal", "% m T MA$^{-1}$"),
        ("beta_norm_toroidal", r"$\beta_N$, toroidal", "% m T MA$^{-1}$"),
        ("beta_thermal_poloidal", r"$\beta_P$, thermal", ""),
        ("beta_poloidal", r"$\beta_P$, total", ""),
        ("te", r"$\langle T_e \rangle$", "keV"),
        ("dene", r"$\langle n_e \rangle$", "m$^{-3}$"),
        (nong, r"$\langle n_{\mathrm{e,line}} \rangle \ / \ n_G$", ""),
        (tepeak, r"$T_{e0} \ / \ \langle T_e \rangle$", ""),
        (nepeak, r"$n_{e0} \ / \ \langle n_{\mathrm{e, vol}} \rangle$", ""),
        ("zeff", r"$Z_{\mathrm{eff}}$", ""),
        (nd_impurities, r"$n_Z \ / \  \langle n_{\mathrm{e, vol}} \rangle$", ""),
        ("t_energy_confinement", r"$\tau_e$", "s"),
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
    if "i_tf_sup" in mfile_data.data:
        i_tf_sup = int(mfile_data.data["i_tf_sup"].get_scan(scan))
    else:
        i_tf_sup = 1

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
    for item in mfile_data.data:
        if "rpf[" in item:
            number_of_coils += 1

    pf_info = [
        (
            mfile_data.data[f"ric[{i:01}]"].get_scan(scan),
            f"PF {i}",
        )
        for i in range(1, number_of_coils)
        if i % 2 != 0
    ]

    if len(pf_info) > 2:
        pf_info_3_a = pf_info[2][0]
        pf_info_3_b = pf_info[2][1]
    else:
        pf_info_3_a = ""
        pf_info_3_b = ""

    t_burn = mfile_data.data["t_burn"].get_scan(scan) / 3600.0

    if "i_tf_bucking" in mfile_data.data:
        i_tf_bucking = int(mfile_data.data["i_tf_bucking"].get_scan(scan))
    else:
        i_tf_bucking = 1

    # Get superconductor material (i_tf_sc_mat)
    # If i_tf_sc_mat not present, assume resistive
    if "i_tf_sc_mat" in mfile_data.data:
        i_tf_sc_mat = int(mfile_data.data["i_tf_sc_mat"].get_scan(scan))
    else:
        i_tf_sc_mat = 0

    if i_tf_sc_mat > 0:
        tftype = dicts["DICT_TF_TYPE"][
            str(int(mfile_data.data["i_tf_sc_mat"].get_scan(scan)))
        ]
    else:
        tftype = "Resistive Copper"

    vssoft = mfile_data.data["vsres"].get_scan(scan) + mfile_data.data[
        "vsind"
    ].get_scan(scan)

    sig_case = 1.0e-6 * mfile_data.data[f"sig_tf_tresca_max({i_tf_bucking})"].get_scan(
        scan
    )
    sig_cond = 1.0e-6 * mfile_data.data[
        f"sig_tf_tresca_max({i_tf_bucking + 1})"
    ].get_scan(scan)

    if i_tf_sup == 1:
        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vstot", "Available flux swing", "Wb"),
            (t_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            ("bmaxtfrp", "Peak field at conductor (w. rip.)", "T"),
            ("iooic", r"I/I$_{\mathrm{crit}}$", ""),
            ("tmargtf", "TF Temperature margin", "K"),
            ("tmargoh", "CS Temperature margin", "K"),
            (sig_cond, "TF Cond max TRESCA stress", "MPa"),
            (sig_case, "TF Case max TRESCA stress", "MPa"),
            ("whttf/n_tf_coils", "Mass per TF coil", "kg"),
        ]

    else:
        p_cp_resistive = 1.0e-6 * mfile_data.data["p_cp_resistive"].get_scan(scan)
        p_tf_leg_resistive = 1.0e-6 * mfile_data.data["p_tf_leg_resistive"].get_scan(
            scan
        )
        pres_joints = 1.0e-6 * mfile_data.data["pres_joints"].get_scan(scan)
        fcoolcp = 100.0 * mfile_data.data["fcoolcp"].get_scan(scan)

        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vstot", "Available flux swing", "Wb"),
            (t_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            ("bmaxtf", "Peak field at conductor (w. rip.)", "T"),
            ("c_tf_total", "TF coil currents sum", "A"),
            ("", "", ""),
            ("#TF coil forces/stresses", "", ""),
            (sig_cond, "TF conductor max TRESCA stress", "MPa"),
            (sig_case, "TF bucking max TRESCA stress", "MPa"),
            (fcoolcp, "CP cooling fraction", "%"),
            ("vcool", "Maximum coolant flow speed", "ms$^{-1}$"),
            (p_cp_resistive, "CP Resisitive heating", "MW"),
            (
                p_tf_leg_resistive,
                "legs Resisitive heating (all legs)",
                "MW",
            ),
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
        / mfile_data.data["fusion_power"].get_scan(scan)
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
        ("p_plasma_inner_rad_mw", "Inner zone radiation", "MW"),
        ("p_plasma_rad_mw", "Total radiation in LCFS", "MW"),
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
            r"Fusion-to-electric efficiency $\frac{P_{\mathrm{e,net}}}{P_{\mathrm{fus}}}$",
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
    lhcd = False
    iccd = False

    if (iefrf == 5) or (iefrf == 8):
        nbi = True
        axis.text(-0.05, 1, "Neutral Beam Current Drive:", ha="left", va="center")
    if (iefrf == 3) or (iefrf == 7) or (iefrf == 10) or (iefrf == 11) or (iefrf == 13):
        ecrh = True
        axis.text(-0.05, 1, "Electron Cyclotron Current Drive:", ha="left", va="center")
    if iefrf == 12:
        ebw = True
        axis.text(-0.05, 1, "Electron Bernstein Wave Drive:", ha="left", va="center")
    if iefrf in [1, 4, 6]:
        lhcd = True
        axis.text(
            -0.05,
            1,
            "Lower Hybrid Current Drive:",
            ha="left",
            va="center",
        )
    if iefrf == 2:
        iccd = True
        axis.text(-0.05, 1, "Ion Cyclotron Current Drive:", ha="left", va="center")

    if "iefrffix" in mfile_data.data:
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
        if iefrffix in [1, 4, 6]:
            secondary_heating = "LHCD"
        if iefrffix == 2:
            secondary_heating = "ICCD"

    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    pinjie = mfile_data.data["pinjmw"].get_scan(scan)
    pdivt = mfile_data.data["pdivt"].get_scan(scan)
    pdivr = pdivt / mfile_data.data["rmajor"].get_scan(scan)

    if mfile_data.data["iefrffix"].get_scan(scan) != 0:
        pinjmwfix = mfile_data.data["pinjmwfix"].get_scan(scan)

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
    if "plhthresh" in mfile_data.data:
        pthresh = mfile_data.data["plhthresh"].get_scan(scan)
    else:
        pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)
    flh = pdivt / pthresh

    hstar = mfile_data.data["hstar"].get_scan(scan)

    if ecrh:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootstrap_current_fraction", "Bootstrap fraction", ""),
            ("aux_current_fraction", "Auxiliary fraction", ""),
            ("inductive_current_fraction", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "effcd",
                "Current drive efficiency",
                "A W$^{-1}$",
            ),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{\langle n \rangle R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        # iefrffix is now always in the MFILE with = 0 meaning no fixed heating
        if mfile_data.data["iefrffix"].get_scan(scan) != 0:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if nbi:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootstrap_current_fraction", "Bootstrap fraction", ""),
            ("aux_current_fraction", "Auxiliary fraction", ""),
            ("inductive_current_fraction", "Inductive fraction", ""),
            ("gamnb", "NB gamma", "$10^{20}$ A W$^{-1}$ m$^{-2}$"),
            ("beam_energy", "NB energy", "keV"),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{\langle n \rangle R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if mfile_data.data["iefrffix"].get_scan(scan) != 0:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if ebw:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootstrap_current_fraction", "Bootstrap fraction", ""),
            ("aux_current_fraction", "Auxiliary fraction", ""),
            ("inductive_current_fraction", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "gamcd",
                "Normalised current drive efficiency",
                "(10$^{20}$ A/(Wm$^{2}$))",
            ),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{\langle n \rangle R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if lhcd:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootstrap_current_fraction", "Bootstrap fraction", ""),
            ("aux_current_fraction", "Auxiliary fraction", ""),
            ("inductive_current_fraction", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "gamcd",
                "Normalised current drive efficiency",
                "(10$^{20}$ A/(Wm$^{2}$))",
            ),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{\langle n \rangle R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if iccd:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("pheat", "Power for heating only", "MW"),
            ("bootstrap_current_fraction", "Bootstrap fraction", ""),
            ("aux_current_fraction", "Auxiliary fraction", ""),
            ("inductive_current_fraction", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "gamcd",
                "Normalised current drive efficiency",
                "(10$^{20}$ A/(Wm$^{2}$))",
            ),
            (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
            (
                pdivnr,
                r"$\frac{P_{\mathrm{div}}}{\langle n \rangle R_{0}}$",
                r"$\times 10^{-20}$ MW m$^{2}$",
            ),
            (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
            (hstar, "H* (non-rad. corr.)", ""),
        ]
        if "iefrffix" in mfile_data.data:
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


def plot_bootstrap_comparison(axis, mfile_data, scan):
    """Function to plot a scatter box plot of bootstrap current fractions.

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use
    """

    boot_ipdg = mfile_data.data["bscf_iter89"].get_scan(scan)
    boot_sauter = mfile_data.data["bscf_sauter"].get_scan(scan)
    boot_nenins = mfile_data.data["bscf_nevins"].get_scan(scan)
    boot_wilson = mfile_data.data["bscf_wilson"].get_scan(scan)
    boot_sakai = mfile_data.data["bscf_sakai"].get_scan(scan)
    boot_aries = mfile_data.data["bscf_aries"].get_scan(scan)
    boot_andrade = mfile_data.data["bscf_andrade"].get_scan(scan)
    boot_hoang = mfile_data.data["bscf_hoang"].get_scan(scan)
    boot_wong = mfile_data.data["bscf_wong"].get_scan(scan)
    boot_gi_I = mfile_data.data["bscf_gi_i"].get_scan(scan)  # noqa: N806
    boot_gi_II = mfile_data.data["bscf_gi_ii"].get_scan(scan)  # noqa: N806

    # Data for the box plot
    data = {
        "IPDG": boot_ipdg,
        "Sauter": boot_sauter,
        "Nevins": boot_nenins,
        "Wilson": boot_wilson,
        "Sakai": boot_sakai,
        "ARIES": boot_aries,
        "Andrade": boot_andrade,
        "Hoang": boot_hoang,
        "Wong": boot_wong,
        "Gi-I": boot_gi_I,
        "Gi-II": boot_gi_II,
    }
    # Create the violin plot
    axis.violinplot(data.values(), showextrema=False)

    # Create the box plot
    axis.boxplot(
        data.values(), showfliers=True, showmeans=True, meanline=True, widths=0.3
    )

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_bootstrap = np.mean(data_values)
    std_bootstrap = np.std(data_values)
    median_bootstrap = np.median(data_values)

    # Plot average, standard deviation, and median as text
    axis.text(
        1.02, 0.2, f"Average: {avg_bootstrap:.4f}", transform=axis.transAxes, fontsize=9
    )
    axis.text(
        1.02,
        0.15,
        f"Standard Dev: {std_bootstrap:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        1.02,
        0.1,
        f"Median: {median_bootstrap:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Bootstrap Current Fraction Comparison")
    axis.set_ylabel("Bootstrap Current Fraction")
    axis.set_xlim([0.5, 1.5])
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def plot_h_threshold_comparison(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int, u_seed=None
) -> None:
    """
    Function to plot a scatter box plot of L-H threshold power comparisons.

    Arguments:
        axis (plt.Axes): Axis object to plot to.
        mfile_data (mf.MFile): MFILE data object.
        scan (int): Scan number to use.
    """
    iter_nominal = mfile_data.data["pthrmw(1)"].get_scan(scan)
    iter_upper = mfile_data.data["pthrmw(2)"].get_scan(scan)
    iter_lower = mfile_data.data["pthrmw(3)"].get_scan(scan)
    iter_1997_1 = mfile_data.data["pthrmw(4)"].get_scan(scan)
    iter_1997_2 = mfile_data.data["pthrmw(5)"].get_scan(scan)
    martin_nominal = mfile_data.data["pthrmw(6)"].get_scan(scan)
    martin_upper = mfile_data.data["pthrmw(7)"].get_scan(scan)
    martin_lower = mfile_data.data["pthrmw(8)"].get_scan(scan)
    snipes_nominal = mfile_data.data["pthrmw(9)"].get_scan(scan)
    snipes_upper = mfile_data.data["pthrmw(10)"].get_scan(scan)
    snipes_lower = mfile_data.data["pthrmw(11)"].get_scan(scan)
    snipes_closed_nominal = mfile_data.data["pthrmw(12)"].get_scan(scan)
    snipes_closed_upper = mfile_data.data["pthrmw(13)"].get_scan(scan)
    snipes_closed_lower = mfile_data.data["pthrmw(14)"].get_scan(scan)
    hubbard_nominal = mfile_data.data["pthrmw(15)"].get_scan(scan)
    hubbard_lower = mfile_data.data["pthrmw(16)"].get_scan(scan)
    hubbard_upper = mfile_data.data["pthrmw(17)"].get_scan(scan)
    hubbard_2017 = mfile_data.data["pthrmw(18)"].get_scan(scan)
    martin_aspect_nominal = mfile_data.data["pthrmw(19)"].get_scan(scan)
    martin_aspect_upper = mfile_data.data["pthrmw(20)"].get_scan(scan)
    martin_aspect_lower = mfile_data.data["pthrmw(21)"].get_scan(scan)

    # Data for the box plot
    data = {
        "ITER 1996 Nominal": iter_nominal,
        "ITER 1996 Upper": iter_upper,
        "ITER 1996 Lower": iter_lower,
        "ITER 1997 (1)": iter_1997_1,
        "ITER 1997 (2)": iter_1997_2,
        "Martin Nominal": martin_nominal,
        "Martin Upper": martin_upper,
        "Martin Lower": martin_lower,
        "Snipes Nominal": snipes_nominal,
        "Snipes Upper": snipes_upper,
        "Snipes Lower": snipes_lower,
        "Snipes Closed Divertor Nominal": snipes_closed_nominal,
        "Snipes Closed Divertor Upper": snipes_closed_upper,
        "Snipes Closed Divertor Lower": snipes_closed_lower,
        "Hubbard Nominal (I-mode)": hubbard_nominal,
        "Hubbard Lower (I-mode)": hubbard_lower,
        "Hubbard Upper (I-mode)": hubbard_upper,
        "Hubbard 2017 (I-mode)": hubbard_2017,
        "Martin Aspect Corrected Nominal": martin_aspect_nominal,
        "Martin Aspect Corrected Upper": martin_aspect_upper,
        "Martin Aspect Corrected Lower": martin_aspect_lower,
    }

    # Create the violin plot
    axis.violinplot(data.values(), showextrema=False)

    # Create the box plot
    axis.boxplot(
        data.values(), showfliers=True, showmeans=True, meanline=True, widths=0.3
    )

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    generator = np.random.default_rng(seed=u_seed)
    x_values = generator.normal(loc=1, scale=0.01, size=len(data.values()))
    for index, (key, value) in enumerate(data.items()):
        if "ITER 1996" in key:
            color = "blue"
        elif "ITER 1997" in key:
            color = "cyan"
        elif "Martin" in key and "Aspect" not in key:
            color = "green"
        elif "Snipes" in key and "Closed" not in key:
            color = "red"
        elif "Snipes Closed" in key:
            color = "orange"
        elif "Martin Aspect" in key:
            color = "yellow"
        elif "Hubbard" in key and "2017" not in key:
            color = "purple"
        elif "Hubbard 2017" in key:
            color = "magenta"
        else:
            color = colors[index]
        axis.scatter(x_values[index], value, color=color, label=key, alpha=1.0)
        axis.legend(loc="upper left", bbox_to_anchor=(-1.1, 1), ncol=2)

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_threshold = np.mean(data_values)
    std_threshold = np.std(data_values)
    median_threshold = np.median(data_values)

    # Plot average, standard deviation, and median as text
    axis.text(
        -0.45,
        0.15,
        f"Average: {avg_threshold:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        -0.45,
        0.1,
        f"Standard Dev: {std_threshold:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        -0.45,
        0.05,
        f"Median: {median_threshold:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("L-H Threshold Comparison")
    axis.set_ylabel("L-H threshold power [MW]")
    axis.set_xlim([0.5, 1.5])
    axis.set_xticks([])
    axis.set_xticklabels([])

    # Add background color
    axis.set_facecolor("#f0f0f0")


def plot_density_limit_comparison(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int
) -> None:
    """
    Function to plot a scatter box plot of different density limit comparisons.

    Arguments:
        axis (plt.Axes): Axis object to plot to.
        mfile_data (mf.MFile): MFILE data object.
        scan (int): Scan number to use.
    """
    old_asdex = mfile_data.data["dlimit(1)"].get_scan(scan)
    borrass_iter_i = mfile_data.data["dlimit(2)"].get_scan(scan)
    borrass_iter_ii = mfile_data.data["dlimit(3)"].get_scan(scan)
    jet_edge_radiation = mfile_data.data["dlimit(4)"].get_scan(scan)
    jet_simplified = mfile_data.data["dlimit(5)"].get_scan(scan)
    hugill_murakami = mfile_data.data["dlimit(6)"].get_scan(scan)
    greenwald = mfile_data.data["dlimit(7)"].get_scan(scan)
    asdex_new = mfile_data.data["dlimit(8)"].get_scan(scan)

    # Data for the box plot
    data = {
        "Old ASDEX": old_asdex,
        "Borrass ITER I": borrass_iter_i,
        "Borrass ITER II": borrass_iter_ii,
        "JET Edge Radiation": jet_edge_radiation,
        "JET Simplified": jet_simplified,
        "Hugill-Murakami": hugill_murakami,
        "Greenwald": greenwald,
        "ASDEX New": asdex_new,
    }

    # Create the violin plot
    axis.violinplot(data.values(), showextrema=False)

    # Create the box plot
    axis.boxplot(
        data.values(), showfliers=True, showmeans=True, meanline=True, widths=0.3
    )

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_density_limit = np.mean(data_values)
    std_density_limit = np.std(data_values)
    median_density_limit = np.median(data_values)

    # Plot average, standard deviation, and median as text
    axis.text(
        1.02,
        0.2,
        rf"Average: {avg_density_limit * 1e-20:.4f} $\times 10^{{20}}$",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        1.02,
        0.15,
        rf"Standard Dev: {std_density_limit * 1e-20:.4f} $\times 10^{{20}}$",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        1.02,
        0.1,
        rf"Median: {median_density_limit * 1e-20:.4f} $\times 10^{{20}}$",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_yscale("log")
    axis.set_title("Density Limit Comparison")
    axis.set_ylabel(r"Density Limit [$10^{20}$ m$^{-3}$]")
    axis.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x * 1e-20:.1f}"))
    axis.set_xlim([0.5, 1.5])
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def main_plot(
    fig1,
    fig2,
    fig3,
    fig4,
    fig5,
    fig6,
    fig7,
    m_file_data,
    scan,
    imp="../data/lz_non_corona_14_elements/",
    demo_ranges=False,
    colour_scheme=1,
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
    plot_1 = fig3.add_subplot(121, aspect="equal")
    poloidal_cross_section(plot_1, m_file_data, scan, demo_ranges, colour_scheme)

    # Plot toroidal cross-section
    plot_2 = fig3.add_subplot(122, aspect="equal")
    toroidal_cross_section(plot_2, m_file_data, scan, demo_ranges, colour_scheme)
    # fig3.subplots_adjust(bottom=-0.2, top = 0.9, left = 0.1, right = 0.9)

    # Plot color key
    plot_3 = fig3.add_subplot(222)
    plot_3.set_position([0.5, 0.5, 0.5, 0.5])
    color_key(plot_3, m_file_data, scan, colour_scheme)

    # Plot density profiles
    plot_4 = fig2.add_subplot(231)  # , aspect= 0.05)
    plot_4.set_position([0.075, 0.55, 0.25, 0.4])
    plot_neprofile(plot_4, demo_ranges)

    # Plot temperature profiles
    plot_5 = fig2.add_subplot(232)
    plot_5.set_position([0.375, 0.55, 0.25, 0.4])
    plot_teprofile(plot_5, demo_ranges)

    # Plot impurity profiles
    plot_8 = fig2.add_subplot(233)
    plot_8.set_position([0.7, 0.45, 0.25, 0.5])
    plot_radprofile(plot_8, m_file_data, scan, imp, demo_ranges)

    # Plot current density profile
    plot_7 = fig2.add_subplot(4, 3, 10)
    plot_7.set_position([0.075, 0.125, 0.25, 0.15])
    plot_jprofile(plot_7)

    # Plot q profile
    plot_6 = fig2.add_subplot(4, 3, 12)
    plot_6.set_position([0.7, 0.125, 0.25, 0.15])
    plot_qprofile(plot_6, demo_ranges)

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
    fig1.subplots_adjust(wspace=0.25, hspace=0.25)

    # Can only plot WP and turn structure if superconducting coil at the moment
    if m_file_data.data["i_tf_sup"].get_scan(scan) == 1:
        # TF coil with WP
        plot_7 = fig4.add_subplot(211, aspect="equal")
        plot_7.set_position([0.05, 0.5, 0.8, 0.4])
        plot_tf_wp(plot_7, m_file_data, scan)

        # TF coil turn structure
        plot_8 = fig4.add_subplot(325, aspect="equal")
        plot_8.set_position([0.1, 0.1, 0.3, 0.3])
        plot_tf_turn(plot_8, m_file_data, scan)

    plot_9 = fig5.add_subplot(221)
    plot_bootstrap_comparison(plot_9, m_file_data, scan)

    plot_10 = fig5.add_subplot(224)
    plot_h_threshold_comparison(plot_10, m_file_data, scan)

    plot_11 = fig6.add_subplot(221)
    plot_density_limit_comparison(plot_11, m_file_data, scan)

    plot_12 = fig7.add_subplot(111)
    plot_current_profiles_over_time(plot_12, m_file_data, scan)


def main(args=None):
    # TODO The use of globals here isn't ideal, but is required to get main()
    # working with minimal changes. Should be converted to class structure
    args = parse_args(args)
    colour_scheme = int(args.colour)
    # read MFILE
    m_file = mf.MFile(args.f) if args.f != "" else mf.MFile("MFILE.DAT")

    scan = args.n if args.n else -1

    demo_ranges = bool(args.DEMO_ranges)

    # Check for Copper magnets
    if "i_tf_sup" in m_file.data:
        i_tf_sup = int(m_file.data["i_tf_sup"].get_scan(scan))
    else:
        i_tf_sup = 1

    # Check WP configuration
    if "i_tf_wp_geom" in m_file.data:
        i_tf_wp_geom = int(m_file.data["i_tf_wp_geom"].get_scan(scan))
    else:
        i_tf_wp_geom = 0

    global dr_bore
    global dr_cs
    global dr_cs_tf_gap
    global dr_tf_inboard
    global dr_shld_vv_gap_inboard
    global ddwi
    global dr_shld_inboard
    global dr_blkt_inboard
    global dr_fw_inboard
    global dr_fw_plasma_gap_inboard
    global rmajor
    global rminor
    global dr_fw_plasma_gap_outboard
    global dr_fw_outboard
    global dr_blkt_outboard
    global dr_shld_outboard
    global ddwi
    global dr_shld_vv_gap_outboard
    global dr_tf_outboard
    global r_cryostat_inboard
    global z_cryostat_half_inside
    global dr_cryostat
    global j_plasma_0

    dr_bore = m_file.data["dr_bore"].get_scan(scan)
    dr_cs = m_file.data["dr_cs"].get_scan(scan)
    dr_cs_tf_gap = m_file.data["dr_cs_tf_gap"].get_scan(scan)
    dr_tf_inboard = m_file.data["dr_tf_inboard"].get_scan(scan)
    dr_shld_vv_gap_inboard = m_file.data["dr_shld_vv_gap_inboard"].get_scan(scan)
    dr_shld_inboard = m_file.data["dr_shld_inboard"].get_scan(scan)
    dr_blkt_inboard = m_file.data["dr_blkt_inboard"].get_scan(scan)
    dr_fw_inboard = m_file.data["dr_fw_inboard"].get_scan(scan)
    dr_fw_plasma_gap_inboard = m_file.data["dr_fw_plasma_gap_inboard"].get_scan(scan)
    rmajor = m_file.data["rmajor"].get_scan(scan)
    rminor = m_file.data["rminor"].get_scan(scan)
    dr_fw_plasma_gap_outboard = m_file.data["dr_fw_plasma_gap_outboard"].get_scan(scan)
    dr_fw_outboard = m_file.data["dr_fw_outboard"].get_scan(scan)
    dr_blkt_outboard = m_file.data["dr_blkt_outboard"].get_scan(scan)
    dr_shld_outboard = m_file.data["dr_shld_outboard"].get_scan(scan)
    dr_shld_vv_gap_outboard = m_file.data["dr_shld_vv_gap_outboard"].get_scan(scan)
    dr_tf_outboard = m_file.data["dr_tf_outboard"].get_scan(scan)
    r_cryostat_inboard = m_file.data["r_cryostat_inboard"].get_scan(scan)
    z_cryostat_half_inside = m_file.data["z_cryostat_half_inside"].get_scan(scan)
    dr_cryostat = m_file.data["dr_cryostat"].get_scan(scan)
    j_plasma_0 = m_file.data["j_plasma_0"].get_scan(scan)

    # Magnets related
    global n_tf_coils
    global wwp1
    global wwp2
    global dr_tf_wp
    global tinstf
    global thkcas
    global casthi

    n_tf_coils = m_file.data["n_tf_coils"].get_scan(scan)
    if i_tf_sup == 1:  # If superconducting magnets
        wwp1 = m_file.data["wwp1"].get_scan(scan)
        if i_tf_wp_geom == 1:
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

    iefrf = int(m_file.data["iefrf"].get_scan(scan))
    iefrffix = int(m_file.data["iefrffix"].get_scan(scan))

    if (iefrf in [5, 8]) or (iefrffix in [5, 8]):
        nbshield = m_file.data["nbshield"].get_scan(scan)
        rtanbeam = m_file.data["rtanbeam"].get_scan(scan)
        rtanmax = m_file.data["rtanmax"].get_scan(scan)
        beamwd = m_file.data["beamwd"].get_scan(scan)
    else:
        nbshield = rtanbeam = rtanmax = beamwd = 0.0

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
    global nd_fuel_ions
    global dene
    global te0
    global ti
    global te
    global fgwped_out
    global fgwsep_out
    global tratio

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
    nd_fuel_ions = m_file.data["nd_fuel_ions"].get_scan(scan)
    dene = m_file.data["dene"].get_scan(scan)
    te0 = m_file.data["te0"].get_scan(scan)
    ti = m_file.data["ti"].get_scan(scan)
    te = m_file.data["te"].get_scan(scan)
    fgwped_out = m_file.data["fgwped_out"].get_scan(scan)
    fgwsep_out = m_file.data["fgwsep_out"].get_scan(scan)
    tratio = m_file.data["tratio"].get_scan(scan)

    # Plasma
    global triang
    global alphaj
    global q0
    global q95
    global plasma_current_MA
    global a_plasma_poloidal

    triang = m_file.data["triang95"].get_scan(scan)
    alphaj = m_file.data["alphaj"].get_scan(scan)
    q0 = m_file.data["q0"].get_scan(scan)
    q95 = m_file.data["q95"].get_scan(scan)
    plasma_current_MA = m_file.data["plasma_current_ma"].get_scan(scan)
    a_plasma_poloidal = m_file.data["a_plasma_poloidal"].get_scan(scan)

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
    global f_sync_reflect
    global bt
    global vol_plasma
    f_sync_reflect = m_file.data["f_sync_reflect"].get_scan(scan)
    bt = m_file.data["bt"].get_scan(scan)
    vol_plasma = m_file.data["vol_plasma"].get_scan(scan)

    # Build the dictionaries of radial and vertical build values and cumulative values
    global vertical_upper
    if int(m_file.data["i_single_null"].get_scan(scan)) == 0:
        vertical_upper = [
            "rminor*kappa",
            "vgaptop",
            "divfix",
            "shldtth",
            "d_vv_top",
            "vgap_vv_thermalshield",
            "thshield_vb",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
        ]
    else:
        vertical_upper = [
            "rminor*kappa",
            "vgaptop",
            "fwtth",
            "blnktth",
            "dr_shld_blkt_gap",
            "shldtth",
            "d_vv_top",
            "vgap_vv_thermalshield",
            "thshield_vb",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
        ]

    radial = {}
    cumulative_radial = {}
    subtotal = 0
    for item in RADIAL_BUILD:
        if item == "rminori" or item == "rminoro":
            build = m_file.data["rminor"].get_scan(scan)
        elif item == "vvblgapi" or item == "vvblgapo":
            build = m_file.data["dr_shld_blkt_gap"].get_scan(scan)
        elif "dr_vv_inboard" in item:
            build = m_file.data["dr_vv_inboard"].get_scan(scan)
        elif "dr_vv_outboard" in item:
            build = m_file.data["dr_vv_outboard"].get_scan(scan)
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

    # read MFILE
    # m_file = mf.MFile(args.f)
    # scan = scan

    # create main plot
    page1 = plt.figure(figsize=(12, 9), dpi=80)
    page2 = plt.figure(figsize=(12, 9), dpi=80)
    page3 = plt.figure(figsize=(12, 9), dpi=80)
    page4 = plt.figure(figsize=(12, 9), dpi=80)
    page5 = plt.figure(figsize=(12, 9), dpi=80)
    page6 = plt.figure(figsize=(12, 9), dpi=80)
    page7 = plt.figure(figsize=(12, 9), dpi=80)

    # run main_plot
    main_plot(
        page1,
        page2,
        page3,
        page4,
        page5,
        page6,
        page7,
        m_file,
        scan=scan,
        demo_ranges=demo_ranges,
        colour_scheme=colour_scheme,
    )

    # with bpdf.PdfPages(args.o) as pdf:
    with bpdf.PdfPages(args.f + "SUMMARY.pdf") as pdf:
        pdf.savefig(page1)
        pdf.savefig(page2)
        pdf.savefig(page3)
        pdf.savefig(page4)
        pdf.savefig(page5)
        pdf.savefig(page6)
        pdf.savefig(page7)

    # show fig if option used
    if args.show:
        plt.show(block=True)

    plt.close(page1)
    plt.close(page2)
    plt.close(page3)
    plt.close(page4)
    plt.close(page5)
    plt.close(page6)
    plt.close(page7)


if __name__ == "__main__":
    main()
