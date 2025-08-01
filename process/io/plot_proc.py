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
import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle
from matplotlib.path import Path

import process.confinement_time as confine
import process.io.mfile as mf
from process.fortran import physics_variables
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
from process.objectives import OBJECTIVE_NAMES
from process.superconducting_tf_coil import SUPERCONDUCTING_TF_TYPES
from process.utilities.f2py_string_patch import f2py_compatible_to_string

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
    "z_plasma_xpoint_lower",
    "dz_xpoint_divertor",
    "dz_divertor",
    "dz_shld_lower",
    "dz_vv_lower",
    "dz_shld_vv_gap",
    "dz_shld_thermal",
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


def plot_main_power_flow(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int, fig: plt.Figure
) -> None:
    """
    Plots the main power flow diagram for the fusion reactor, including plasma, heating and current drive,
    first wall, blanket, vacuum vessel, divertor, coolant pumps, turbine, generator, and auxiliary systems.
    Annotates the diagram with power values and draws arrows to indicate power flows.

    Args:
        axis (plt.Axes): The matplotlib axis object to plot on.
        mfile_data (mf.MFile): The MFILE data object containing power flow parameters.
        scan (int): The scan number to use for extracting data.
        fig (plt.Figure): The matplotlib figure object for additional annotations.
    """

    axis.text(
        0.05,
        0.95,
        "* Components do not represent the design",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )

    # ==========================================
    # Plasma
    # ===========================================

    # Load the plasma image
    plasma = mpimg.imread(
        resources.path("process.io", "plasma.png")
    )  # Use importlib.resources to locate the image

    # Display the plasma image over the figure, not the axes
    new_ax = axis.inset_axes(
        [-0.15, 0.6, 0.45, 0.45], transform=axis.transAxes, zorder=1
    )
    new_ax.imshow(plasma)
    new_ax.axis("off")

    # Add fusion power to plasma
    axis.text(
        0.22,
        0.75,
        "$P_{{fus}}$\n"
        f"{mfile_data.data['p_fusion_total_mw'].get_scan(scan):.2f} MW",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )
    # Load the neutron image
    neutron = mpimg.imread(resources.path("process.io", "neutron.png"))

    new_ax = axis.inset_axes(
        [0.2, 0.85, 0.03, 0.03], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(neutron)
    new_ax.axis("off")

    # Add lost alpha power
    axis.text(
        0.22,
        0.81,
        "$P_{\\alpha,{loss}}$\n"
        f"{mfile_data.data['p_fw_alpha_mw'].get_scan(scan):.2f} MW",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )

    # Add radiation power to plasma
    axis.text(
        0.22,
        0.69,
        f"$P_{{{{rad}}}}$\n{mfile_data.data['p_plasma_rad_mw'].get_scan(scan):.2f} MW",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )

    # Add photon image to plasma
    axis.text(
        0.34,
        0.71,
        "$\\gamma$",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=12,
    )

    # Draw from gamma arrow bend towards divertor
    axis.annotate(
        "",
        xy=(0.35, 0.55),
        xytext=(0.35, 0.695),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "blue",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Add separatrix power to plasma
    axis.text(
        0.22,
        0.63,
        f"$P_{{{{sep}}}}$\n{mfile_data.data['p_plasma_separatrix_mw'].get_scan(scan):.2f} MW",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )

    # Draw from separatrix power to arrow bend
    axis.annotate(
        "",
        xy=(0.3725, 0.65),
        xytext=(0.3, 0.65),
        xycoords=fig.transFigure,
        arrowprops={
            "color": "pink",
            "arrowstyle": "-",  # No arrow head
            "linewidth": 2.0,
        },
    )

    # Draw from separatrix arrow bend to the divertor
    axis.annotate(
        "",
        xy=(0.37, 0.55),
        xytext=(0.37, 0.65),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",  # solid filled head
            "color": "pink",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw neutron arrow from plasma
    axis.annotate(
        "",
        xy=(0.95, 0.76),
        xytext=(0.31, 0.76),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "grey",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw arrow from main neutron arrow down to divertor
    axis.annotate(
        "",
        xy=(0.39, 0.55),
        xytext=(0.39, 0.76),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "grey",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw radiation arrow from plasma
    axis.annotate(
        "",
        xy=(0.56, 0.695),
        xytext=(0.3, 0.695),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "blue",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Load the alpha particle image
    alpha = mpimg.imread(
        resources.path("process.io", "alpha_particle.PNG")
    )  # Use importlib.resources to locate the image

    # Display the alpha particle image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.16, 0.95, 0.025, 0.025], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(alpha)
    new_ax.axis("off")

    # Hide the axes for a cleaner look
    axis.axis("off")

    # Draw alpha particle arrow from plasma
    axis.annotate(
        "",
        xy=(0.56, 0.83),
        xytext=(0.3, 0.83),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "red",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plot neutron power from plasma to box
    axis.text(
        0.37,
        0.775,
        f"$P_{{\\text{{neutron}}}}$:\n{mfile_data.data['p_neutron_total_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "grey",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # ===========================================

    # =========================================
    # Heating and current drive systems
    # =========================================

    # Add HCD primary injected power
    axis.text(
        0.0725,
        0.83,
        f"$P_{{\\text{{HCD,primary}}}}$: {mfile_data.data['p_hcd_primary_injected_mw'].get_scan(scan) + mfile_data.data['p_hcd_primary_extra_heat_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add HCD secondary injected power
    axis.text(
        0.0725,
        0.725,
        f"$P_{{\\text{{HCD,secondary}}}}$: {mfile_data.data['p_hcd_secondary_injected_mw'].get_scan(scan) + mfile_data.data['p_hcd_secondary_extra_heat_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Load the HCD injector image
    hcd_injector_1 = hcd_injector_2 = mpimg.imread(
        resources.path("process.io", "hcd_injector.png")
    )  # Use importlib.resources to locate the image

    # Display the injector image over the figure, not the axes
    new_ax = axis.inset_axes(
        [-0.2, 0.8, 0.15, 0.15], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(hcd_injector_1)
    new_ax.axis("off")
    new_ax = axis.inset_axes(
        [-0.2, 0.5, 0.15, 0.5], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(hcd_injector_2)
    new_ax.axis("off")

    # Draw a dashed line with an arrow tip coming from the left of each injector
    for y in [0.875, 0.75]:
        axis.annotate(
            "",
            xy=(-0.2, y),
            xytext=(-0.28, y),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

    # Plot line from HCD power supply to bend for injected
    axis.plot(
        [-0.28, -0.28],
        [0.875, 0.5],
        transform=axis.transAxes,
        color="black",
        linewidth=1.5,
        zorder=3,
        clip_on=False,
    )

    # Plot the HCD power supply box
    axis.text(
        0.04,
        0.45,
        "\n\nH&CD Power Supply\n\n",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
        zorder=4,
    )

    # Draw arrow from HCD box going to primary HCD losses
    axis.annotate(
        "",
        xy=(0.2, 0.5),
        xytext=(0.1, 0.5),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.2",
            "color": "black",
            "linestyle": "--",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plot electric power losses for secondary HCD
    axis.text(
        0.2,
        0.435,
        f"$P_{{\\text{{secondary,loss}}}}$:\n{mfile_data.data['p_hcd_secondary_electric_mw'].get_scan(scan) * (1.0 - mfile_data.data['eta_hcd_secondary_injector_wall_plug'].get_scan(scan)):.2f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Draw an arrow from HCD secondary losses to the total secondary heat power
    axis.annotate(
        "",
        xy=(0.25, 0.3),
        xytext=(0.25, 0.43),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Draw an arrow from HCD primary losses bend to the total secondary heat power
    axis.annotate(
        "",
        xy=(0.28, 0.3),
        xytext=(0.28, 0.5),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Draw line from HCD primary losses to the arrow bend
    axis.annotate(
        "",
        xy=(0.26, 0.5),
        xytext=(0.28, 0.5),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "black",
            "linestyle": "--",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw arrow frim HCD power supply to secondary HCD losses
    axis.annotate(
        "",
        xy=(0.2, 0.46),
        xytext=(0.1, 0.46),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.2",
            "color": "black",
            "linestyle": "--",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plot electric power losses for primary HCD
    axis.text(
        0.2,
        0.485,
        f"$P_{{\\text{{primary,loss}}}}$:\n{mfile_data.data['p_hcd_primary_electric_mw'].get_scan(scan) * (1.0 - mfile_data.data['eta_hcd_primary_injector_wall_plug'].get_scan(scan)):.2f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Draw arrow from HCD primary electric box to HCD power supply box
    axis.annotate(
        "",
        xy=(0.06, 0.45),
        xytext=(0.06, 0.38),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "->",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
        },
    )

    # Draw arrow from HCD secondary electric box to HCD power supply box
    axis.annotate(
        "",
        xy=(0.12, 0.45),
        xytext=(0.12, 0.38),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "->",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
        },
    )

    # Plot HCD secondary losses box
    axis.text(
        0.12,
        0.35,
        f"$P_{{\\text{{secondary}}}}$:\n{mfile_data.data['p_hcd_secondary_electric_mw'].get_scan(scan):.2f} MWe \n$\\eta$: {mfile_data.data['eta_hcd_secondary_injector_wall_plug'].get_scan(scan):.2f}",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Plot HCD primary electric box
    axis.text(
        0.025,
        0.35,
        f"$P_{{\\text{{primary}}}}$:\n{mfile_data.data['p_hcd_primary_electric_mw'].get_scan(scan):.2f} MWe\n$\\eta$: {mfile_data.data['eta_hcd_primary_injector_wall_plug'].get_scan(scan):.2f}",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # =============================================

    # =============================================
    # Low grade heat total
    # =============================================

    # Plot box of total low grade secondary heat
    axis.text(
        0.325,
        0.225,
        f"\n\nTotal Low Grade Secondary Heat\n\n {mfile_data.data['p_plant_secondary_heat_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="center",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
        zorder=4,
    )

    # =============================================

    # ==========================================
    # Power conversion systems
    # ===========================================

    # Load the turbine image
    turbine = mpimg.imread(
        resources.path("process.io", "turbine.png")
    )  # Use importlib.resources to locate the image

    # Display the turbine image over the figure, not the axes
    new_ax = axis.inset_axes(
        [1.1, 0.0, 0.15, 0.15], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(turbine)
    new_ax.axis("off")

    # Plot the total primary thermal power box
    axis.text(
        0.9,
        0.25,
        f"$P_{{\\text{{primary,thermal}}}}$:\n{mfile_data.data['p_plant_primary_heat_mw'].get_scan(scan):.2f} MW \n$\\eta_{{\\text{{turbine}}}}$: {mfile_data.data['eta_turbine'].get_scan(scan):.3f}",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Draw arrow from bend to turbine inlet
    axis.annotate(
        "",
        xy=(0.925, 0.165),
        xytext=(0.96, 0.165),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Total primary thermal to turbine inlet line bend
    axis.annotate(
        "",
        xy=(0.96, 0.245),
        xytext=(0.96, 0.1625),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Load the generator image
    generator = mpimg.imread(
        resources.path("process.io", "generator.png")
    )  # Use importlib.resources to locate the image

    # Display the generator image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.96, 0.0, 0.15, 0.15], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(generator)
    new_ax.axis("off")

    # Generator to gross electric power
    axis.annotate(
        "",
        xy=(0.75, 0.17),
        xytext=(0.79, 0.17),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Generator labels
    axis.text(
        0.79,
        0.16,
        "Generator",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        zorder=20,
    )

    # Connector from turbine to generator
    axis.annotate(
        "",
        xy=(0.85, 0.17),
        xytext=(0.925, 0.17),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "black",
            "linewidth": 7.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Turbine to loss power
    axis.annotate(
        "",
        xy=(0.91, 0.08),
        xytext=(0.91, 0.13),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "dashed",
        },
    )

    # Load the pylon image
    pylon = mpimg.imread(
        resources.path("process.io", "pylon.png")
    )  # Use importlib.resources to locate the image

    # Display the pylon image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.925, -0.1, 0.1, 0.1], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(pylon)
    new_ax.axis("off")

    # Plot the gross electric power box
    axis.text(
        0.68,
        0.15,
        f"$P_{{\\text{{gross}}}}$:\n{mfile_data.data['p_plant_electric_gross_mw'].get_scan(scan):.2f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lime",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Gross to net electric power
    axis.annotate(
        "",
        xy=(0.72, 0.08),
        xytext=(0.72, 0.15),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plot the turbine loss box
    axis.text(
        0.875,
        0.05,
        f"$P_{{\\text{{loss}}}}$:\n{mfile_data.data['p_turbine_loss_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Shield primary thermal to plant total primary thermal arrow
    axis.annotate(
        "",
        xy=(0.95, 0.3),
        xytext=(0.95, 0.55),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plot the net electric power box
    axis.text(
        0.68,
        0.05,
        f"$P_{{\\text{{net,electric}}}}$:\n{mfile_data.data['p_plant_electric_net_mw'].get_scan(scan):.2f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lime",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Plot the recirculated electric power box
    axis.text(
        0.575,
        0.14,
        (
            f"$P_{{\\text{{recirc,electric}}}}$:\n{mfile_data.data['p_plant_electric_recirc_mw'].get_scan(scan):.2f} MWe\n"
            f"$f_{{\\text{{recirc}}}}$:\n{mfile_data.data['f_p_plant_electric_recirc'].get_scan(scan):.2f}"
        ),
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lime",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Gross to recirculated power arrow
    axis.annotate(
        "",
        xy=(0.635, 0.17),
        xytext=(0.675, 0.17),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Recirculated to pumps electric
    axis.annotate(
        "",
        xy=(0.7, 0.225),
        xytext=(0.635, 0.185),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Recirculated power to HCD secondary electric arrow bend
    axis.annotate(
        "",
        xy=(0.14, 0.2),
        xytext=(0.57, 0.2),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Recirculated power to HCD primary electric arrow bend
    axis.annotate(
        "",
        xy=(0.08, 0.18),
        xytext=(0.57, 0.18),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Arrow to primary HCD electric from bend
    axis.annotate(
        "",
        xy=(0.08, 0.35),
        xytext=(0.08, 0.1775),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Arrow to secondary HCD electric from bend
    axis.annotate(
        "",
        xy=(0.14, 0.35),
        xytext=(0.14, 0.2),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # ==========================================

    # ================================
    # First wall, blanket and shield
    # ================================

    # Load the first wall image
    fw = mpimg.imread(
        resources.path("process.io", "fw.png")
    )  # Use importlib.resources to locate the image

    # Display the first wall image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.4, 0.625, 0.4, 0.4], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(fw)
    new_ax.axis("off")

    # Add first wall label above image
    axis.text(
        0.5,
        0.9,
        "First Wall",
        fontsize=11,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
    )

    # Alpha power incident on first wall box
    axis.text(
        0.46,
        0.85,
        f"$P_{{\\text{{FW, }}\\alpha}}$:\n{mfile_data.data['p_fw_alpha_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "red",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Neutron power incident on first wall box
    axis.text(
        0.46,
        0.775,
        f"$P_{{\\text{{FW,nuclear}}}}$:\n{mfile_data.data['p_fw_nuclear_heat_total_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "grey",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Plot radiation power incident on first wall box
    axis.text(
        0.46,
        0.71,
        f"$P_{{\\text{{FW,rad}}}}$:\n{mfile_data.data['p_fw_rad_total_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "dodgerblue",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Draw arrow from FW to heat depsoited box
    axis.annotate(
        "",
        xy=(0.61, 0.585),
        xytext=(0.61, 0.65),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw arrow from Blanket to heat deposited box
    axis.annotate(
        "",
        xy=(0.81, 0.585),
        xytext=(0.81, 0.63),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw arrow from shield to heat deposited box
    axis.annotate(
        "",
        xy=(0.92, 0.59),
        xytext=(0.92, 0.62),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # First wall heat deposited box
    axis.text(
        0.5,
        0.555,
        f"Primary thermal\n(inc pump): {mfile_data.data['p_fw_heat_deposited_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "linewidth": 2,
        },
    )

    # Blanket heat deposited box
    axis.text(
        0.7,
        0.555,
        f"Primary thermal\n(inc pump): {mfile_data.data['p_blkt_heat_deposited_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "linewidth": 2,
        },
    )

    # Shield heat deposited box
    axis.text(
        0.875,
        0.555,
        f"Primary thermal:\n{mfile_data.data['p_shld_heat_deposited_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "linewidth": 2,
        },
    )

    # Draw arrow from FW primary heat box to blanket and FW primary heat deposited box
    axis.annotate(
        "",
        xy=(0.65, 0.52),
        xytext=(0.62, 0.55),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw arrow from blanket primary heat box to blanket and FW primary heat deposited box
    axis.annotate(
        "",
        xy=(0.68, 0.52),
        xytext=(0.7, 0.55),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Draw a downward arrow from the primary thermal box to the right side of the generator
    axis.annotate(
        "",
        xy=(0.825, 0.57),
        xytext=(0.87, 0.57),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Connect blanket thermal heat deposited to the shield heat deposited
    axis.annotate(
        "",
        xy=(0.625, 0.57),
        xytext=(0.695, 0.57),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Connect first wall thermal heat deposited to the blanket heat deposited
    axis.annotate(
        "",
        xy=(0.56, 0.52),
        xytext=(0.56, 0.55),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "orange",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # FW and blanket heat deposited box
    axis.text(
        0.6,
        0.49,
        f"Primary thermal (inc pump): {mfile_data.data['p_fw_blkt_heat_deposited_mw'].get_scan(scan):.2f} MWth\n",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Load the blanket image
    blanket = mpimg.imread(
        resources.path("process.io", "blanket_with_coolant.png")
    )  # Use importlib.resources to locate the image

    # Display the blanket image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.75, 0.625, 0.4, 0.4], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(blanket)
    new_ax.axis("off")

    # Add blanket label above image
    axis.text(
        0.7,
        0.9,
        "Blanket",
        fontsize=11,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
    )

    # Plot the nuclear heat total from blanket
    axis.text(
        0.625,
        0.775,
        (
            f"$P_{{\\text{{Blkt,nuclear}}}}$:\n{mfile_data.data['p_blkt_nuclear_heat_total_mw'].get_scan(scan):.2f} MW \n"
            f"$P_{{\\text{{Blkt,multiplication}}}}$:\n{mfile_data.data['p_blkt_multiplication_mw'].get_scan(scan):.2f} MW\n"
            f"$f_{{\\text{{multiplication}}}}$:\n{mfile_data.data['f_p_blkt_multiplication'].get_scan(scan):.2f}"
        ),
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "grey",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Load the vacuum vessel image
    vv = mpimg.imread(
        resources.path("process.io", "vv.png")
    )  # Use importlib.resources to locate the image

    # Display the vacuum vessel image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.975, 0.625, 0.4, 0.4], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(vv)
    new_ax.axis("off")

    # Add vacuum vessel label above image
    axis.text(
        0.85,
        0.9,
        "Vacuum Vessel",
        fontsize=11,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
    )

    # Plot the secondary heat from the shield
    axis.text(
        0.38,
        0.375,
        f"$P_{{\\text{{shld,secondary}}}}$:\n{mfile_data.data['p_shld_secondary_heat_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Shield secondary power box to secondary heat total
    axis.annotate(
        "",
        xy=(0.4, 0.3),
        xytext=(0.4, 0.37),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Arrow from shield bend to sheidl secondary heat
    axis.annotate(
        "",
        xy=(0.445, 0.39),
        xytext=(0.85, 0.39),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Line from shield to arrow bend for secondary heat
    axis.annotate(
        "",
        xy=(0.85, 0.385),
        xytext=(0.85, 0.625),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # ============================================
    # Divertor
    # ============================================

    axis.text(
        0.325,
        0.48,
        "Divertor",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=1000,  # bring to front
        fontsize=11,
        color="white",  # make text white
    )

    # Load the divertor image
    divertor = mpimg.imread(
        resources.path("process.io", "divertor.png")
    )  # Use importlib.resources to locate the image

    # Display the divertor image over the figure, not the axes
    new_ax = axis.inset_axes([0.1, 0.4, 0.3, 0.25], transform=axis.transAxes, zorder=10)
    new_ax.imshow(divertor)
    new_ax.axis("off")

    # Total divertor radiation power box
    axis.text(
        0.29,
        0.57,
        f"$P_{{\\text{{div,rad}}}}$:\n{mfile_data.data['p_div_rad_total_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "dodgerblue",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Divertor nuclear heat total box
    axis.text(
        0.4,
        0.58,
        f"$P_{{\\text{{div,nuclear}}}}$:\n{mfile_data.data['p_div_nuclear_heat_total_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "grey",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Divertor primary thermal heat deposited box
    axis.text(
        0.44,
        0.46,
        (
            f"Primary thermal (inc pump):\n{mfile_data.data['p_div_heat_deposited_mw'].get_scan(scan):.2f} MWth\n"
            f"Solid angle fraction: {mfile_data.data['f_ster_div_single'].get_scan(scan):.3f}\n"
            f"Primary heat fraction: {mfile_data.data['f_p_div_primary_heat'].get_scan(scan):.3f}"
        ),
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "linewidth": 2,
        },
        zorder=100,
    )

    # Divertor secondary heat box
    axis.text(
        0.3,
        0.375,
        f"$P_{{\\text{{div,secondary}}}}$:\n{mfile_data.data['p_div_secondary_heat_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Divertor to divertor secondary heat arrow
    axis.annotate(
        "",
        xy=(0.33, 0.405),
        xytext=(0.33, 0.5),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Divertor to divertor primary thermal heat arrow
    axis.annotate(
        "",
        xy=(0.445, 0.5),
        xytext=(0.4, 0.5),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "orange",
            "linewidth": 2.0,
            "zorder": 50,
            "fill": True,
        },
    )

    # Divertor secondary heat to total secondary heat arrow
    axis.annotate(
        "",
        xy=(0.33, 0.3),
        xytext=(0.33, 0.375),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # ===========================================

    # ===========================================
    # Coolant pumps
    # ===========================================

    # Divertor coolant pump box
    axis.text(
        0.55,
        0.33,
        f"$P_{{\\text{{div,pump}}}}$: {mfile_data.data['p_div_coolant_pump_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Divertor pump box to divertor primary heat deposited box
    axis.annotate(
        "",
        xy=(0.57, 0.46),
        xytext=(0.57, 0.35),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Coolant pumps total to divertor pump box
    axis.annotate(
        "",
        xy=(0.64, 0.34),
        xytext=(0.7, 0.34),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Pumps total to shield bump box arrow
    axis.annotate(
        "",
        xy=(0.875, 0.34),
        xytext=(0.81, 0.34),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Shield coolant pump box
    axis.text(
        0.875,
        0.325,
        f"$P_{{\\text{{shld,pump}}}}$:\n{mfile_data.data['p_shld_coolant_pump_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # FW and Blanket coolant pumps total
    axis.text(
        0.725,
        0.4,
        f"$P_{{\\text{{FW + Blkt}}}}$:\n{mfile_data.data['p_fw_blkt_coolant_pump_mw'].get_scan(scan):.2f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # FW and Blanket coolant pumps total to FW and Blanket heat deposited box
    axis.annotate(
        "",
        xy=(0.75, 0.49),
        xytext=(0.75, 0.44),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 3.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Coolant pumps total to blanket and FW pump
    axis.annotate(
        "",
        xy=(0.75, 0.4),
        xytext=(0.75, 0.36),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Shield pump to sheild primary thermal
    axis.annotate(
        "",
        xy=(0.9, 0.54),
        xytext=(0.9, 0.36),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Coolant pumps total electric box
    axis.text(
        0.7,
        0.225,
        (
            f"Coolant pumps electric:\n{mfile_data.data['p_coolant_pump_elec_total_mw'].get_scan(scan):.3f} MWe\n"
            f"$\\eta$: {mfile_data.data['eta_coolant_pump_electric'].get_scan(scan):.3f}"
        ),
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lime",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Coolant pumps total
    axis.text(
        0.7,
        0.325,
        f"Coolant pumps total:\n{mfile_data.data['p_coolant_pump_total_mw'].get_scan(scan):.3f} MW",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Electric recirculated to pumps total arrow
    axis.annotate(
        "",
        xy=(0.75, 0.325),
        xytext=(0.75, 0.275),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Coolant pumps losses total box
    axis.text(
        0.5,
        0.235,
        f"Coolant pumps losses total:\n{mfile_data.data['p_coolant_pump_loss_total_mw'].get_scan(scan):.3f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 0.8,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # Coolant electric to pump losses arrow
    axis.annotate(
        "",
        xy=(0.645, 0.25),
        xytext=(0.695, 0.25),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # Coolant losses to secondary heat total arrow
    axis.annotate(
        "",
        xy=(0.405, 0.25),
        xytext=(0.4975, 0.25),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # ============================================

    # ===========================================
    # Plant core systems
    # ===========================================

    # Cryo Plant box
    axis.text(
        0.45,
        0.05,
        f"Cryo Plant:\n{mfile_data.data['p_cryo_plant_electric_mw'].get_scan(scan):.3f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "burlywood",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Recirculated power to cryo plant arrow
    axis.annotate(
        "",
        xy=(0.49, 0.09),
        xytext=(0.575, 0.14),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Tritium Plant box
    axis.text(
        0.35,
        0.05,
        f"Tritium Plant:\n{mfile_data.data['p_tritium_plant_electric_mw'].get_scan(scan):.3f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "burlywood",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # # Recirculated power to tritium plant arrow
    axis.annotate(
        "",
        xy=(0.4, 0.09),
        xytext=(0.575, 0.15),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Vacuum Pumps box
    axis.text(
        0.55,
        0.05,
        f"Vacuum pumps:\n{mfile_data.data['vachtmw'].get_scan(scan):.3f} MWe",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "burlywood",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Recirculated power to vacuum pumps arrow
    axis.annotate(
        "",
        xy=(0.6, 0.08),
        xytext=(0.6, 0.14),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Plant base load box
    axis.text(
        0.15,
        0.05,
        (
            f"Plant base load:\n{mfile_data.data['p_plant_electric_base_total_mw'].get_scan(scan):.3f} MWe\n"
            f"Minimum base load:\n{mfile_data.data['p_plant_electric_base'].get_scan(scan) * 1.0e-6:.3f} MWe\n"
            f"Plant floor power density:\n{mfile_data.data['pflux_plant_floor_electric'].get_scan(scan) * 1.0e-3:.3f} kW/m^2"
        ),
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "burlywood",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # Recirculated power to plant base load arrow
    axis.annotate(
        "",
        xy=(0.28, 0.09),
        xytext=(0.575, 0.17),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # HCD secondary heat box
    axis.text(
        0.46,
        0.285,
        f"$P_{{\\text{{HCD,loss}}}}$:\n{mfile_data.data['p_hcd_secondary_heat_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 0.8,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # FW to HCD secondary heat arrow
    axis.annotate(
        "",
        xy=(0.47, 0.32),
        xytext=(0.47, 0.65),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # HCD loss to total secondary heat
    axis.annotate(
        "",
        xy=(0.41, 0.295),
        xytext=(0.455, 0.295),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # TF nuclear heat box
    axis.text(
        0.155,
        0.25,
        f"$P_{{\\text{{TF,nuclear}}}}$:\n{mfile_data.data['p_tf_nuclear_heat_mw'].get_scan(scan):.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
            "linestyle": "dashed",
        },
    )

    # TF nuclear heat to secondary heat total box arrow
    axis.annotate(
        "",
        xy=(0.245, 0.265),
        xytext=(0.215, 0.265),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
            "linestyle": "--",
        },
    )

    # ===========================================


def plot_main_plasma_information(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int, colour_scheme: int, fig: plt.Figure
) -> None:
    """
    Plots the main plasma information including plasma shape, geometry, currents, heating,
    confinement, and other relevant plasma parameters.

    Args:
        axis (plt.Axes): The matplotlib axis object to plot on.
        mfile_data (mf.MFile): The MFILE data object containing plasma parameters.
        scan (int): The scan number to use for extracting data.
        colour_scheme (int): The colour scheme to use for plots.
        fig (plt.Figure): The matplotlib figure object for additional annotations.
    """
    # Import key variables
    triang = mfile_data.data["triang"].get_scan(scan)
    kappa = mfile_data.data["kappa"].get_scan(scan)

    # Remove the axes
    axis.axis("off")

    # Plot the main plasma shape
    plot_plasma(axis, mfile_data, scan, colour_scheme)

    # Get the plasma permieter points for the core plasma region
    pg = plasma_geometry(
        rmajor=mfile_data.data["rmajor"].get_scan(scan),
        rminor=mfile_data.data["rminor"].get_scan(scan)
        * mfile_data.data["radius_plasma_core_norm"].get_scan(scan),
        triang=mfile_data.data["triang"].get_scan(scan),
        kappa=mfile_data.data["kappa"].get_scan(scan),
        i_single_null=mfile_data.data["i_single_null"].get_scan(scan),
        i_plasma_shape=1,
        square=mfile_data.data["plasma_square"].get_scan(scan),
    )
    # Plot the core plasma boundary line
    axis.plot(pg.rs, pg.zs, color="black", linestyle="--")

    # Plot the centre of the plasma
    axis.plot(rmajor, 0, "r+", markersize=20, markeredgewidth=2)

    # =========================================

    # Draw a double-ended arrow from the inner plasma edge to the center
    axis.annotate(
        "",
        xy=(rmajor - rminor, 0),  # Inner plasma edge
        xytext=(rmajor, 0),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the minor radius
    axis.text(
        rmajor - rminor / 2,
        -rminor * kappa * 0.08,
        f"$a$: {rminor:.2f} m",
        fontsize=9,
        color="black",
        ha="center",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )

    # ============================================

    # Draw a single-ended arrow from the machien centre to the plasma center
    axis.annotate(
        "",
        xy=(axis.get_xlim()[0], -rminor * 0.3 * kappa),  # Inner plasma edge
        xytext=(rmajor, -rminor * 0.3 * kappa),  # Center
        arrowprops={"arrowstyle": "<-", "color": "black"},
    )

    # Add a label for the major radius
    axis.text(
        rmajor - rminor / 2,
        -rminor * kappa * 0.25,
        f"$R_0$: {rmajor:.2f} m",
        fontsize=9,
        color="black",
        ha="center",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )

    # ============================================

    # Draw a double-ended arrow from the xpoint to the center to show elongation
    axis.annotate(
        "",
        xy=(rmajor - rminor * triang, kappa * rminor),  # Inner plasma edge
        xytext=(rmajor - rminor * triang, 0),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Write the elongation beside the vertical line, position relative to figure axes
    axis.text(
        0.3,
        0.75,
        f"$\\kappa$: {mfile_data.data['kappa'].get_scan(scan):.2f}",
        fontsize=9,
        color="black",
        rotation=270,
        verticalalignment="center",
        transform=axis.transAxes,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )

    # =============================================

    # Draw a double-ended arrow from the inner plasma edge to the center
    axis.annotate(
        "",
        xy=(rmajor - rminor * triang, kappa * rminor * 0.25),  # Inner plasma edge
        xytext=(rmajor, kappa * rminor * 0.25),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Write the triangularity to the left of the cross, position relative to figure axes
    axis.text(
        rmajor - (rminor * triang * 0.75),
        kappa * rminor * 0.3,
        f"$\\delta$: {mfile_data.data['triang'].get_scan(scan):.2f}",
        fontsize=9,
        color="black",
        rotation=0,
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )

    # =============================================

    radius_plasma_core_norm = mfile_data.data["radius_plasma_core_norm"].get_scan(scan)

    # Draw a double-ended arrow for the plasma core region
    axis.annotate(
        "",
        xy=(rmajor, -rminor * 0.1 * kappa),  # Inner plasma edge
        xytext=(rmajor + (rminor * radius_plasma_core_norm), -rminor * 0.1 * kappa),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    # Add a label for core region
    axis.text(
        rmajor + (rminor * radius_plasma_core_norm / 4),
        -rminor * kappa * 0.15,
        f"$\\rho_{{\\text{{core}}}}$: {radius_plasma_core_norm:.2f}",
        fontsize=9,
        color="black",
        rotation=0,
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )

    # ================================================

    # Add plasma volume, areas and shaping information
    textstr_plasma = (
        f"$\\mathbf{{Shaping:}}$\n \n"
        f"$\\kappa_{{95}}$: {mfile_data.data['kappa95'].get_scan(scan):.2f} | $\\delta_{{95}}$: {mfile_data.data['triang95'].get_scan(scan):.2f} | $\\zeta$: {mfile_data.data['plasma_square'].get_scan(scan):.2f}\n"
        f"A: {mfile_data.data['aspect'].get_scan(scan):.2f}\n"
        f"$ V_{{\\text{{p}}}}:$ {mfile_data.data['vol_plasma'].get_scan(scan):.2f}$ \\ m^3$\n"
        f"$ A_{{\\text{{p,surface}}}}:$ {mfile_data.data['a_plasma_surface'].get_scan(scan):.2f}$ \\ m^2$\n"
        f"$ A_{{\\text{{p_poloidal}}}}:$ {mfile_data.data['a_plasma_poloidal'].get_scan(scan):.2f}$ \\ m^2$\n"
        f"$ L_{{\\text{{p_poloidal}}}}:$ {mfile_data.data['len_plasma_poloidal'].get_scan(scan):.2f}$ \\ m$"
    )

    axis.text(
        0.6,
        0.825,
        textstr_plasma,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # ============================================

    # Draw a red arrow coming from the right and pointing at the plasma
    axis.annotate(
        "",
        xy=(rmajor + (rminor * 0.8), kappa * rminor * 0.2),  # Pointing at the plasma
        xytext=(
            rmajor + (rminor * 1.4),
            kappa * rminor * 0.2,
        ),  # Starting point of the arrow
        arrowprops={"facecolor": "red", "edgecolor": "red", "lw": 2},
    )

    # Draw a red arrow coming from the right and pointing at the plasma
    axis.annotate(
        "",
        xy=(rmajor + (rminor * 0.8), -kappa * rminor * 0.2),  # Pointing at the plasma
        xytext=(
            rmajor + (rminor * 1.4),
            -kappa * rminor * 0.2,
        ),  # Starting point of the arrow
        arrowprops={"facecolor": "red", "edgecolor": "red", "lw": 2},
    )

    i_hcd_primary = mfile_data.data["i_hcd_primary"].get_scan(scan)
    i_hcd_secondary = mfile_data.data["i_hcd_secondary"].get_scan(scan)

    # Determine heating type for primary and secondary systems
    if i_hcd_primary in [5, 8]:
        primary_heating = "NBI"
    elif i_hcd_primary in [3, 7, 10, 11, 13]:
        primary_heating = "ECRH"
    elif i_hcd_primary == 12:
        primary_heating = "EBW"
    elif i_hcd_primary in [1, 4, 6]:
        primary_heating = "LHCD"
    elif i_hcd_primary == 2:
        primary_heating = "ICCD"
    else:
        primary_heating = ""

    if i_hcd_secondary in [5, 8]:
        secondary_heating = "NBI"
    elif i_hcd_secondary in [3, 7, 10, 11, 13]:
        secondary_heating = "ECRH"
    elif i_hcd_secondary == 12:
        secondary_heating = "EBW"
    elif i_hcd_secondary in [1, 4, 6]:
        secondary_heating = "LHCD"
    elif i_hcd_secondary == 2:
        secondary_heating = "ICCD"
    else:
        secondary_heating = ""

    # Add heating and current drive information
    textstr_hcd = (
        f"$\\mathbf{{Heating \\ & \\ current \\ drive:}}$\n \n"
        f"Total injected heat: {mfile_data.data['p_hcd_injected_total_mw'].get_scan(scan):.3f} MW                       \n"
        f"Ohmic heating power: {mfile_data.data['p_plasma_ohmic_mw'].get_scan(scan):.3f} MW         \n\n"
        f"$\\mathbf{{Primary \\ system: {primary_heating}}}$ \n"
        f"Current driving power {mfile_data.data['p_hcd_primary_injected_mw'].get_scan(scan):.4f} MW\n"
        f"Extra heat power: {mfile_data.data['p_hcd_primary_extra_heat_mw'].get_scan(scan):.4f} MW\n"
        f"$\\gamma_{{\\text{{CD,prim}}}}$: {mfile_data.data['eta_cd_hcd_primary'].get_scan(scan):.4f} A/W\n"
        f"$\\eta_{{\\text{{CD,prim}}}}$: {mfile_data.data['eta_cd_norm_hcd_primary'].get_scan(scan):.2f} $\\times 10^{{20}}  \\mathrm{{A}} / \\mathrm{{Wm}}^2$\n"
        f"Current driven by primary: {mfile_data.data['c_hcd_primary_driven'].get_scan(scan) / 1e6:.3f} MA\n\n"
        f"$\\mathbf{{Secondary \\ system: {secondary_heating}}}$ \n"
        f"Current driving power {mfile_data.data['p_hcd_secondary_injected_mw'].get_scan(scan):.4f} MW\n"
        f"Extra heat power: {mfile_data.data['p_hcd_secondary_extra_heat_mw'].get_scan(scan):.4f} MW\n"
        f"$\\gamma_{{\\text{{CD,sec}}}}$: {mfile_data.data['eta_cd_hcd_secondary'].get_scan(scan):.4f} A/W\n"
        f"$\\eta_{{\\text{{CD,sec}}}}$: {mfile_data.data['eta_cd_norm_hcd_secondary'].get_scan(scan):.2f} $\\times 10^{{20}}  \\mathrm{{A}} / \\mathrm{{Wm}}^2$\n"
        f"Current driven by secondary: {mfile_data.data['c_hcd_secondary_driven'].get_scan(scan) / 1e6:.3f} MA\n"
    )

    axis.text(
        0.66,
        0.6,
        textstr_hcd,
        fontsize=9,
        verticalalignment="top",
        transform=plt.gcf().transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "paleturquoise",
            "alpha": 1.0,
            "linewidth": 2,
            "edgecolor": "black",  # Set box outline to black
        },
    )

    # Add injected power label
    axis.text(
        0.92,
        0.625,
        "$P_{\\text{inj}}$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ================================================

    # Add beta information
    textstr_beta = (
        f"$\\mathbf{{Beta \\ Information:}}$\n \n"
        f"Total beta,$ \\ \\beta$: {mfile_data.data['beta'].get_scan(scan):.4f}\n"
        f"Thermal beta,$ \\ \\beta_{{\\text{{thermal}}}}$: {mfile_data.data['beta_thermal'].get_scan(scan):.4f}\n"
        f"Toroidal beta,$ \\ \\beta_{{\\text{{t}}}}$: {mfile_data.data['beta_toroidal'].get_scan(scan):.4f}\n"
        f"Poloidal beta,$ \\ \\beta_{{\\text{{p}}}}$: {mfile_data.data['beta_poloidal'].get_scan(scan):.4f}\n"
        f"Fast-alpha beta,$ \\ \\beta_{{\\alpha}}$: {mfile_data.data['beta_fast_alpha'].get_scan(scan):.4f}\n"
        f"Normalised total beta,$ \\ \\beta_{{\\text{{N}}}}$: {mfile_data.data['beta_norm_total'].get_scan(scan):.4f}\n"
        f"Normalised thermal beta,$ \\ \\beta_{{\\text{{N,thermal}}}}$: {mfile_data.data['beta_norm_thermal'].get_scan(scan):.4f}\n"
    )

    axis.text(
        0.025,
        0.95,
        textstr_beta,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightblue",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add beta label
    axis.text(
        0.215,
        0.925,
        "$\\beta$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ================================================

    # Add volt-second information
    textstr_volt_second = (
        f"$\\mathbf{{Volt-second \\ requirements:}}$\n \n"
        f"Total volt-second consumption: {mfile_data.data['vs_plasma_total_required'].get_scan(scan):.4f} Vs                \n"
        f"  - Internal volt-seconds: {mfile_data.data['vs_plasma_internal'].get_scan(scan):.4f} Vs\n"
        f"  - Volt-seconds needed for burn: {mfile_data.data['vs_plasma_burn_required'].get_scan(scan):.4f} Vs\n"
        f"$V_{{\\text{{loop}}}}$: {mfile_data.data['v_plasma_loop_burn'].get_scan(scan):.4f} V\n"
        f"$\\Omega_{{\\text{{p}}}}$: {mfile_data.data['res_plasma'].get_scan(scan):.4e} $\\Omega$\n"
        f"Plasma resistive diffusion time: {mfile_data.data['t_plasma_res_diffusion'].get_scan(scan):.4f} s\n"
        f"Plasma inductance: {mfile_data.data['ind_plasma'].get_scan(scan):.4e} H\n"
        f"Plasma normalised internal inductance: {mfile_data.data['ind_plasma_internal_norm'].get_scan(scan):.4f}\n"
    )

    axis.text(
        0.025,
        0.77,
        textstr_volt_second,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightgreen",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add volt second label
    axis.text(
        0.27,
        0.75,
        "Vs",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =========================================

    # Add divertor information
    textstr_div = (
        f"\n$P_{{\\text{{sep}}}}$: {mfile_data.data['p_plasma_separatrix_mw'].get_scan(scan):.2f} MW           \n"
        f"$\\frac{{P_{{\\text{{sep}}}}}}{{R}}$: {mfile_data.data['p_plasma_separatrix_mw/rmajor'].get_scan(scan):.2f} MW/m               \n"
        f"$\\frac{{P_{{\\text{{sep}}}}}}{{B_T  q_a  R}}$: {mfile_data.data['pdivtbt_over_qar'].get_scan(scan):.2f} MW T/m               "
    )

    axis.text(
        0.35,
        0.12,
        textstr_div,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "orange",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add divertor label
    axis.text(
        0.45,
        0.1,
        "$P_{\\text{div}}$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ================================================

    # Add confinement information
    textstr_confinement = (
        f"$\\mathbf{{Confinement:}}$\n \n"
        f"Confinement scaling law: {mfile_data.data['tauelaw'].get_scan(scan)}\n"
        f"Confinement $H$ factor: {mfile_data.data['hfact'].get_scan(scan):.4f}\n"
        f"Energy confinement time from scaling: {mfile_data.data['t_energy_confinement'].get_scan(scan):.4f} s\n"
        f"Fusion double product: {mfile_data.data['ntau'].get_scan(scan):.4e} s/m³\n"
        f"Lawson Triple product: {mfile_data.data['nttau'].get_scan(scan):.4e} keV·s/m³\n"
        f"Transport loss power assumed in scaling law: {mfile_data.data['p_plasma_loss_mw'].get_scan(scan):.4f} MW\n"
        f"Plasma thermal energy (inc. $\\alpha$), $W$: {mfile_data.data['e_plasma_beta'].get_scan(scan) / 1e9:.4f} GJ\n"
        f"Alpha particle confinement time: {mfile_data.data['t_alpha_confinement'].get_scan(scan):.4f} s"
    )

    axis.text(
        0.025,
        0.575,
        textstr_confinement,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "gainsboro",  # Changed to a not normal color (Aquamarine)
            "alpha": 1.0,
            "linewidth": 2,
            "edgecolor": "black",  # Set box outline to black
        },
    )

    # Add tau label
    axis.text(
        0.3,
        0.55,
        "$\\tau_{\\text{e}} $",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =========================================

    # Load the neutron image
    with resources.path(
        "process.io", "alpha_particle.PNG"
    ) as alpha_particle_image_path:
        # Use importlib.resources to locate the image
        alpha_particle = mpimg.imread(alpha_particle_image_path.open("rb"))

    # Display the neutron image over the figure, not the axes
    new_ax = axis.inset_axes(
        [0.975, 0.275, 0.075, 0.075], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(alpha_particle)
    new_ax.axis("off")

    axis.annotate(
        "",
        xy=(rmajor + rminor, -rminor * kappa * 0.55),  # Pointing at the plasma
        xytext=(
            rmajor + 0.2 * rminor,
            -rminor * kappa * 0.25,
        ),
        arrowprops={"facecolor": "red", "edgecolor": "grey", "lw": 1},
    )

    textstr_alpha = (
        f"$P_{{\\alpha,\\text{{loss}}}}$ {mfile_data.data['p_fw_alpha_mw'].get_scan(scan):.2f} MW \n"
        f"$f_{{\\alpha,\\text{{coupled}}}}$ {mfile_data.data['f_p_alpha_plasma_deposited'].get_scan(scan):.2f}"
    )

    axis.text(
        1.0,
        0.275,
        textstr_alpha,
        fontsize=9,
        verticalalignment="top",
        transform=axis.transAxes,
        bbox={
            "boxstyle": "round",
            "facecolor": "red",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # =========================================
    with resources.path("process.io", "neutron.png") as neutron_image_path:
        neutron = mpimg.imread(neutron_image_path.open("rb"))
    new_ax = axis.inset_axes(
        [0.975, 0.75, 0.075, 0.075], transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(neutron)
    new_ax.axis("off")

    # Draw a red arrow coming from the right and pointing at the plasma
    axis.annotate(
        "",
        xy=(rmajor + rminor, rminor * kappa * 0.65),  # Pointing at the plasma
        xytext=(
            rmajor,
            rminor * kappa * 0.5,
        ),
        arrowprops={"facecolor": "grey", "edgecolor": "grey", "lw": 1},
    )

    textstr_neutron = f"$P_{{\\text{{n,total}}}}$ {mfile_data.data['p_neutron_total_mw'].get_scan(scan):.2f} MW"

    axis.text(
        0.75,
        0.85,
        textstr_neutron,
        fontsize=9,
        verticalalignment="top",
        transform=axis.transAxes,
        bbox={
            "boxstyle": "round",
            "facecolor": "grey",
            "alpha": 0.8,
            "linewidth": 2,
        },
    )

    # ===============================================

    # Add fusion reaction information
    textstr_reactions = (
        f"$\\mathbf{{Fusion \\ Reactions:}}$\n \n"
        f"Fuel mixture: \n"
        f"|  D: {mfile_data.data['f_deuterium'].get_scan(scan):.2f}  |  T: {mfile_data.data['f_tritium'].get_scan(scan):.2f}  |  3He: {mfile_data.data['f_helium3'].get_scan(scan):.2f}  |\n\n"
        f"Fusion Power, $P_{{\\text{{fus}}}}:$ {mfile_data.data['p_fusion_total_mw'].get_scan(scan):.2f} MW\n"
        f"D-T Power, $P_{{\\text{{fus,DT}}}}:$ {mfile_data.data['p_dt_total_mw'].get_scan(scan):.2f} MW\n"
        f"D-D Power, $P_{{\\text{{fus,DD}}}}:$ {mfile_data.data['p_dd_total_mw'].get_scan(scan):.2f} MW\n"
        f"D-3He Power, $P_{{\\text{{fus,D3He}}}}:$ {mfile_data.data['p_dhe3_total_mw'].get_scan(scan):.2f} MW\n"
        f"Alpha Power, $P_{{\\alpha}}:$ {mfile_data.data['p_alpha_total_mw'].get_scan(scan):.2f} MW"
    )

    axis.text(
        0.025,
        0.4,
        textstr_reactions,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "red",
            "alpha": 0.6,
            "linewidth": 2,
        },
    )

    # ================================================

    # Add fuelling information
    textstr_fuelling = (
        f"$\\mathbf{{Fuelling:}}$\n \n"
        f"Plasma mass: {mfile_data.data['m_plasma'].get_scan(scan) * 1000:.4f} g\n"
        f"   - Average mass of all plasma ions: {mfile_data.data['m_ions_total_amu'].get_scan(scan):.3f} amu\n"
        f"Fuel mass: {mfile_data.data['m_plasma_fuel_ions'].get_scan(scan) * 1000:.4f} g\n"
        f"   - Average mass of all fuel ions: {mfile_data.data['m_fuel_amu'].get_scan(scan):.3f} amu\n\n"
        f"Fueling rate: {mfile_data.data['qfuel'].get_scan(scan):.3e} nucleus-pairs/s\n"
        f"Fuel burn-up rate: {mfile_data.data['rndfuel'].get_scan(scan):.3e} reactions/s \n"
        f"Burn-up fraction: {mfile_data.data['burnup'].get_scan(scan):.4f} \n"
    )

    axis.text(
        0.025,
        0.22,
        textstr_fuelling,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "khaki",
            "alpha": 1.0,
            "linewidth": 2,
            "edgecolor": "black",  # Set box outline to black
        },
    )

    # ================================================

    # Add ion density information
    textstr_ions = (
        f"             $\\mathbf{{Ion \\ to \\ electron}}$\n"
        f"             $\\mathbf{{relative \\ number}}$\n"
        f"             $\\mathbf{{densities:}}$\n \n"
        f"             Effective charge: {mfile_data.data['zeff'].get_scan(scan):.3f}\n\n"
        f"             H:    {mfile_data.data['fimp(01)'].get_scan(scan):.4e}\n"
        f"             He:   {mfile_data.data['fimp(02)'].get_scan(scan):.4e}\n"
        f"             Be:   {mfile_data.data['fimp(03)'].get_scan(scan):.4e}\n"
        f"             C:    {mfile_data.data['fimp(04)'].get_scan(scan):.4e}\n"
        f"             N:    {mfile_data.data['fimp(05)'].get_scan(scan):.4e}\n"
        f"             O:    {mfile_data.data['fimp(06)'].get_scan(scan):.4e}\n"
        f"             Ne:   {mfile_data.data['fimp(07)'].get_scan(scan):.4e}\n"
        f"             Si:   {mfile_data.data['fimp(08)'].get_scan(scan):.4e}\n"
        f"             Ar:   {mfile_data.data['fimp(09)'].get_scan(scan):.4e}\n"
        f"             Fe:   {mfile_data.data['fimp(10)'].get_scan(scan):.4e}\n"
        f"             Ni:   {mfile_data.data['fimp(11)'].get_scan(scan):.4e}\n"
        f"             Kr:   {mfile_data.data['fimp(12)'].get_scan(scan):.4e}\n"
        f"             Xe:   {mfile_data.data['fimp(13)'].get_scan(scan):.4e}\n"
        f"              W:    {mfile_data.data['fimp(14)'].get_scan(scan):.4e}"
    )

    axis.text(
        0.805,
        0.335,
        textstr_ions,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "olivedrab",
            "alpha": 0.7,
            "linewidth": 2,
        },
    )

    # Add ion charge label
    axis.text(
        0.815,
        0.29,
        "$Z$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ================================================

    # Add plasma current information
    textstr_currents = (
        f"          $\\mathbf{{Plasma\\ currents:}}$\n\n"
        f"          Plasma current {mfile_data.data['plasma_current_ma'].get_scan(scan):.4f} MA\n"
        f"            - Bootstrap fraction {mfile_data.data['f_c_plasma_bootstrap'].get_scan(scan):.4f}\n"
        f"            - Diamagnetic fraction {mfile_data.data['f_c_plasma_diamagnetic'].get_scan(scan):.4f}\n"
        f"            - Pfirsch-Schlüter fraction {mfile_data.data['f_c_plasma_pfirsch_schluter'].get_scan(scan):.4f}\n"
        f"            - Auxiliary fraction {mfile_data.data['f_c_plasma_auxiliary'].get_scan(scan):.4f}\n"
        f"            - Inductive fraction {mfile_data.data['f_c_plasma_inductive'].get_scan(scan):.4f}"
    )

    axis.text(
        0.765,
        0.95,
        textstr_currents,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "#C8A2C8",  # Hex code for lilac color
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add plasma current label
    axis.text(
        0.78,
        0.9,
        "$I_{\\text{p}} $",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =================================================

    # Add magnetic field information
    textstr_fields = (
        f"$\\mathbf{{Magnetic\\ fields:}}$\n\n"
        f"Toroidal field at $R_0$, $B_{{T}}$: {mfile_data.data['bt'].get_scan(scan):.4f} T                  \n"
        f"Average poloidal field, $B_{{p}}$: {mfile_data.data['bp'].get_scan(scan):.4f} T              \n"
        f"Total field, $B_{{tot}}$: {mfile_data.data['btot'].get_scan(scan):.4f} T                \n"
        f"Vertical field, $B_{{vert}}$: {mfile_data.data['bvert'].get_scan(scan):.4f} T"
    )

    axis.text(
        0.55,
        0.13,
        textstr_fields,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "royalblue",  # Hex code for lilac color
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # Add magnetic field label
    axis.text(
        0.75,
        0.1,
        "$B$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ===========================================

    # Add radiation information
    textstr_radiation = (
        f"           $\\mathbf{{Radiation:}}$\n\n"
        f"           Total radiation power {mfile_data.data['p_plasma_rad_mw'].get_scan(scan):.4f} MW\n"
        f"           Core radiation power {mfile_data.data['p_plasma_inner_rad_mw'].get_scan(scan):.4f} MW\n"
        f"           Edge radiation power {mfile_data.data['p_plasma_outer_rad_mw'].get_scan(scan):.4f} MW\n"
        f"           Synchrotron radiation power {mfile_data.data['p_plasma_sync_mw'].get_scan(scan):.4f} MW\n"
        f"           Synchrotron wall reflectivity {mfile_data.data['f_sync_reflect'].get_scan(scan):.4f}\n"
    )

    axis.text(
        0.72,
        0.81,
        textstr_radiation,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lavender",
            "alpha": 1.0,
            "linewidth": 2,
            "edgecolor": "black",  # Set box outline to black
        },
    )

    # Add radiation label
    axis.text(
        0.725,
        0.76,
        "$\\gamma$",
        fontsize=23,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # ============================================

    # Add L-H threshold information
    textstr_lh = (
        f"$\\mathbf{{L-H \\ threshold:}}$\n\n"
        f"$P_{{\\text{{L-H}}}}:$ {mfile_data.data['p_l_h_threshold_mw'].get_scan(scan):.4f} MW\n"
    )

    axis.text(
        0.22,
        0.4,
        textstr_lh,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "peachpuff",  # Changed color to navy
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    # ======================================

    # Add density limit information
    textstr_density_limit = (
        f"$\\mathbf{{Density \\ limit:}}$\n\n"
        f"$n_{{\\text{{e,limit}}}}: {mfile_data.data['dnelimt'].get_scan(scan):.3e} \\ m^{{-3}}$\n"
    )

    axis.text(
        0.22,
        0.32,
        textstr_density_limit,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "pink",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )


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

    # Find the number of PF circuits, n_pf_cs_plasma_circuits includes the CS and plasma circuits
    n_pf_cs_plasma_circuits = mfile_data.data["n_pf_cs_plasma_circuits"].get_scan(scan)

    # Extract PF circuit times
    # n_pf_cs_plasma_circuits contains the CS and plasma at the end so we subtract 2
    pf_circuits = {}
    for i in range(int(n_pf_cs_plasma_circuits - 2)):
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

    if (mfile_data.data["i_hcd_primary"].get_scan(scan) in [5, 8]) or (
        mfile_data.data["i_hcd_secondary"].get_scan(scan) in [5, 8]
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
            facecolor=TFC_COLOUR[colour_scheme - 1],
        )

    i_hcd_primary = mfile_data.data["i_hcd_primary"].get_scan(scan)
    if (i_hcd_primary == 5) or (i_hcd_primary == 8):
        # Neutral beam geometry
        a = w
        b = dr_tf_outboard
        c = dx_beam_duct + 2 * dx_beam_shield
        d = r3
        e = np.sqrt(a**2 + (d + b) ** 2)
        # Coordinates of the inner and outer edges of the beam at its tangency point
        rinner = radius_beam_tangency - dx_beam_duct
        router = radius_beam_tangency + dx_beam_duct
        beta = np.arccos(rinner / e)
        xinner = rinner * np.cos(beta)
        yinner = rinner * np.sin(beta)
        xouter = router * np.cos(beta)
        youter = router * np.sin(beta)
        # Corner of TF coils
        xcorner = r4
        ycorner = w + dx_beam_shield
        axis.plot(
            [xinner, xcorner], [yinner, ycorner], linestyle="dotted", color="black"
        )
        x = xcorner + c * np.cos(beta) - dx_beam_shield * np.cos(beta)
        y = ycorner + c * np.sin(beta) - dx_beam_shield * np.sin(beta)
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


def plot_n_profiles(prof, demo_ranges, mfile_data, scan):
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
        ni = (ne0 * (nd_fuel_ions / dene)) * (1 - rho**2) ** alphan
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
        rf"$\langle n_{{\text{{e}}}} \rangle$: {mfile_data.data['dene'].get_scan(scan):.3e} m$^{{-3}}$",
        rf"$n_{{\text{{e,0}}}}$: {ne0:.3e} m$^{{-3}}$"
        rf"$\hspace{{4}} \alpha_{{\text{{n}}}}$: {alphan:.3f}",
        rf"$n_{{\text{{e,ped}}}}$: {neped:.3e} m$^{{-3}}$"
        r"$ \hspace{3} \frac{\langle n_i \rangle}{\langle n_e \rangle}$: "
        f"{nd_fuel_ions / dene:.3f}",
        rf"$f_{{\text{{GW e,ped}}}}$: {fgwped_out:.3f}"
        r"$ \hspace{7} \frac{n_{e,0}}{\langle n_e \rangle}$: "
        f"{ne0 / dene:.3f}",
        rf"$\rho_{{\text{{ped,n}}}}$: {rhopedn:.3f}"
        r"$ \hspace{8} \frac{\overline{n_{e}}}{n_{\text{GW}}}$: "
        f"{mfile_data.data['nd_electron_line'].get_scan(scan) / mfile_data.data['dlimit(7)'].get_scan(scan):.3f}",
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
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)

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
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)


def plot_t_profiles(prof, demo_ranges, mfile_data, scan):
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

    te = mfile_data.data["te"].get_scan(scan)
    # Add text box with temperature profile parameters
    textstr_temperature = "\n".join((
        rf"$\langle T_{{\text{{e}}}} \rangle$: {mfile_data.data['te'].get_scan(scan):.3f} keV",
        rf"$T_{{\text{{e,0}}}}$: {te0:.3f} keV"
        rf"$\hspace{{4}} \alpha_{{\text{{T}}}}$: {alphat:.3f}",
        rf"$T_{{\text{{e,ped}}}}$: {teped:.3f} keV"
        r"$ \hspace{4} \frac{\langle T_i \rangle}{\langle T_e \rangle}$: "
        f"{tratio:.3f}",
        rf"$\rho_{{\text{{ped,T}}}}$: {rhopedt:.3f}"
        r"$ \hspace{6} \frac{T_{e,0}}{\langle T_e \rangle}$: "
        f"{te0 / te:.3f}",
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
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    # ---


def plot_qprofile(prof, demo_ranges, mfile_data, scan):
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
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    # ---

    textstr_q = "\n".join((
        r"$q_0$: " + f"{q0:.3f}\n",
        r"$q_{95}$: " + f"{q95:.3f}\n",
        r"$q_{\text{cyl}}$: " + f"{mfile_data.data['qstar'].get_scan(scan):.3f}",
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
        dz_blkt_upper = mfile_data.data["dz_blkt_upper"].get_scan(scan)
    else:
        dz_blkt_upper = 0.0

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
            dz_blkt_upper=dz_blkt_upper,
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
            dz_blkt_upper=dz_blkt_upper,
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


def plot_first_wall_top_down_cross_section(axis, mfile_data, scan):
    # Import required variables
    radius_fw_channel = mfile_data.data["radius_fw_channel"].get_scan(scan) * 100
    dr_fw_wall = mfile_data.data["dr_fw_wall"].get_scan(scan) * 100
    dx_fw_module = mfile_data.data["dx_fw_module"].get_scan(scan) * 100

    # Flot first module
    axis.add_patch(
        patches.Rectangle(
            xy=(0, 0),
            width=dx_fw_module,
            height=2 * (dr_fw_wall + radius_fw_channel),
            edgecolor="black",
            facecolor="gray",
        )
    )

    # Plot cooling channel in first module
    axis.add_patch(
        patches.Circle(
            xy=(dx_fw_module / 2, dr_fw_wall + radius_fw_channel),
            radius=radius_fw_channel,
            edgecolor="black",
            facecolor="#b87333",
        )
    )

    # Plot second module
    axis.add_patch(
        patches.Rectangle(
            xy=(dx_fw_module, 0),
            width=dx_fw_module,
            height=2 * (dr_fw_wall + radius_fw_channel),
            edgecolor="black",
            facecolor="gray",
        )
    )

    # Plot cooling channel in second module
    axis.add_patch(
        patches.Circle(
            xy=(dx_fw_module + dx_fw_module / 2, dr_fw_wall + radius_fw_channel),
            radius=radius_fw_channel,
            edgecolor="black",
            facecolor="#b87333",
        )
    )

    # Draw radius line in the second circle
    axis.plot(
        [
            dx_fw_module + dx_fw_module / 2,
            dx_fw_module + dx_fw_module / 2 + radius_fw_channel * np.cos(np.pi / 4),
        ],
        [
            dr_fw_wall + radius_fw_channel,
            dr_fw_wall + radius_fw_channel + radius_fw_channel * np.sin(np.pi / 4),
        ],
        color="black",
        linestyle="--",
        label=f"$r_{{channel}}$ = {radius_fw_channel:.3f} cm",
    )

    # Draw width line below the second module
    axis.plot(
        [0, 0],
        [0, 0],
        color="black",
        label=f"$w_{{module}}$ = {dx_fw_module:.3f} cm",
    )
    axis.annotate(
        "",
        xy=(dx_fw_module, -0.2),
        xytext=(2 * dx_fw_module, -0.2),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Draw dotted line above the channel
    axis.plot(
        [dx_fw_module * 1.5, dx_fw_module * 1.5],
        [2 * radius_fw_channel + dr_fw_wall, 2 * (radius_fw_channel + dr_fw_wall)],
        color="black",
        linestyle="dotted",
        label=rf"$\Delta r_{{wall}}$ = {dr_fw_wall:.3f} cm",
    )

    # Draw dotted line below the channel
    axis.plot(
        [dx_fw_module * 1.5, dx_fw_module * 1.5],
        [0, dr_fw_wall],
        color="black",
        linestyle="dotted",
    )
    # Plot a dot in the center of the second channel
    axis.plot(
        dx_fw_module + dx_fw_module / 2,
        dr_fw_wall + radius_fw_channel,
        marker="o",
        color="black",
    )

    # Add the legend to the plot
    axis.legend()
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    axis.set_xlabel("X [cm]")
    axis.set_ylabel("R [cm]")
    axis.set_title("First Wall Top-Down Cross Section")
    axis.set_xlim([-1, 2 * dx_fw_module + 1])
    axis.set_ylim([-1, 2 * (dr_fw_wall + radius_fw_channel) + 1])


def plot_first_wall_poloidal_cross_section(axis, mfile_data, scan):
    # Import required variables
    radius_fw_channel = mfile_data.data["radius_fw_channel"].get_scan(scan)
    dr_fw_wall = mfile_data.data["dr_fw_wall"].get_scan(scan)
    dx_fw_module = mfile_data.data["dx_fw_module"].get_scan(scan)
    len_fw_channel = mfile_data.data["len_fw_channel"].get_scan(scan)
    temp_fw_coolant_in = mfile_data.data["temp_fw_coolant_in"].get_scan(scan)
    temp_fw_coolant_out = mfile_data.data["temp_fw_coolant_out"].get_scan(scan)
    i_fw_coolant_type = mfile_data.data["i_fw_coolant_type"].get_scan(scan)
    temp_fw_peak = mfile_data.data["temp_fw_peak"].get_scan(scan)
    pres_fw_coolant = mfile_data.data["pres_fw_coolant"].get_scan(scan)
    n_fw_outboard_channels = mfile_data.data["n_fw_outboard_channels"].get_scan(scan)
    n_fw_inboard_channels = mfile_data.data["n_fw_inboard_channels"].get_scan(scan)

    # Plot first wall structure facing plasma
    axis.add_patch(
        patches.Rectangle(
            xy=(0, 0),
            width=dr_fw_wall,
            height=len_fw_channel,
            edgecolor="black",
            facecolor="gray",
        )
    )

    # Plot the cooling channel
    axis.add_patch(
        patches.Rectangle(
            xy=(dr_fw_wall, 0),
            width=2 * radius_fw_channel,
            height=len_fw_channel,
            edgecolor="black",
            facecolor="#b87333",  # Copper color
        )
    )

    # Plot the back wall of the first wall
    axis.add_patch(
        patches.Rectangle(
            xy=(dr_fw_wall + 2 * radius_fw_channel, 0),
            width=dr_fw_wall,
            height=len_fw_channel,
            edgecolor="black",
            facecolor="grey",
        )
    )

    # Draw an upward pointing arrow
    axis.arrow(
        dx_fw_module + 0.5 * dr_fw_wall,
        dr_fw_wall + radius_fw_channel,
        0,
        len_fw_channel / 6,
        head_width=dr_fw_wall,
        head_length=len_fw_channel / 20,
        fc="black",
        ec="black",
    )

    # Add the inlet temperature beside the arrow
    axis.text(
        dx_fw_module + 2 * dr_fw_wall,
        dr_fw_wall + radius_fw_channel + len_fw_channel / 6,
        f"$T_{{inlet}} = ${temp_fw_coolant_in:.2f} K",
        ha="left",
        va="bottom",
        fontsize=10,
        color="black",
    )

    # Draw a right pointing arrow
    axis.arrow(
        dx_fw_module + 0.5 * dr_fw_wall,
        len_fw_channel,
        2 * dr_fw_wall,
        0,
        head_width=len_fw_channel / 30,
        head_length=dr_fw_wall,
        fc="black",
        ec="black",
        linewidth=5,  # Thicker stem
    )

    # Add the outlet temperature beside the arrow
    axis.text(
        dx_fw_module + 0.5 * dr_fw_wall,
        len_fw_channel * 0.9,
        f"$T_{{outlet}} = ${temp_fw_coolant_out:.2f} K",
        ha="left",
        va="bottom",
        fontsize=10,
        color="black",
    )

    textstr_fw = "\n".join((
        rf"Coolant type: {i_fw_coolant_type}",
        rf"$T_{{FW,peak}}$: {temp_fw_peak:.3f} K",
        rf"$P_{{FW}}$: {pres_fw_coolant / 1e3:.3f} kPa",
        rf"$P_{{FW}}$: {pres_fw_coolant / 1e5:.3f} bar",
        rf"$N_{{outboard}}$: {n_fw_outboard_channels}",
        rf"$N_{{inboard}}$: {n_fw_inboard_channels}",
    ))

    props_fw = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    axis.text(
        -0.7,
        0.05,
        textstr_fw,
        transform=axis.transAxes,
        fontsize=11,
        verticalalignment="bottom",
        bbox=props_fw,
    )

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("First Wall Poloidal Cross Section")
    axis.set_xlim([-0.01, (dx_fw_module + radius_fw_channel * 2) + 0.01])
    axis.set_ylim([-0.2, len_fw_channel + 0.2])


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
        dz_blkt_upper = mfile_data.data["dz_blkt_upper"].get_scan(scan)
        tfwvt = mfile_data.data["dz_fw_upper"].get_scan(scan)
    else:
        dz_blkt_upper = tfwvt = 0.0

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
            dz_blkt_upper=dz_blkt_upper,
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
            dz_blkt_upper=dz_blkt_upper,
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
    x1 = mfile_data.data["r_tf_arc(1)"].get_scan(scan)
    y1 = mfile_data.data["z_tf_arc(1)"].get_scan(scan)
    x2 = mfile_data.data["r_tf_arc(2)"].get_scan(scan)
    y2 = mfile_data.data["z_tf_arc(2)"].get_scan(scan)
    x3 = mfile_data.data["r_tf_arc(3)"].get_scan(scan)
    y3 = mfile_data.data["z_tf_arc(3)"].get_scan(scan)
    x4 = mfile_data.data["r_tf_arc(4)"].get_scan(scan)
    y4 = mfile_data.data["z_tf_arc(4)"].get_scan(scan)
    x5 = mfile_data.data["r_tf_arc(5)"].get_scan(scan)
    y5 = mfile_data.data["z_tf_arc(5)"].get_scan(scan)
    dr_shld_thermal_inboard = mfile_data.data["dr_shld_thermal_inboard"].get_scan(scan)
    dr_shld_thermal_outboard = mfile_data.data["dr_shld_thermal_outboard"].get_scan(
        scan
    )
    dr_tf_shld_gap = mfile_data.data["dr_tf_shld_gap"].get_scan(scan)
    if y3 != 0:
        print("TF coil geometry: The value of z_tf_arc(3) is not zero, but should be.")

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


def plot_tf_wp(axis, mfile_data, scan: int, fig) -> None:
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
    dx_tf_wp_primary_toroidal = mfile_data.data["dx_tf_wp_primary_toroidal"].get_scan(
        scan
    )
    dx_tf_wp_secondary_toroidal = mfile_data.data[
        "dx_tf_wp_secondary_toroidal"
    ].get_scan(scan)
    dr_tf_wp_with_insulation = mfile_data.data["dr_tf_wp_with_insulation"].get_scan(
        scan
    )
    r_tf_wp_inboard_inner = mfile_data.data["r_tf_wp_inboard_inner"].get_scan(scan)
    dx_tf_wp_insulation = mfile_data.data["dx_tf_wp_insulation"].get_scan(scan)
    n_tf_coil_turns = round(mfile_data.data["n_tf_coil_turns"].get_scan(scan))
    i_tf_wp_geom = round(mfile_data.data["i_tf_wp_geom"].get_scan(scan))
    i_tf_sup = round(mfile_data.data["i_tf_sup"].get_scan(scan))
    i_tf_case_geom = mfile_data.data["i_tf_case_geom"].get_scan(scan)
    i_tf_turns_integer = mfile_data.data["i_tf_turns_integer"].get_scan(scan)
    b_tf_inboard_peak = mfile_data.data["b_tf_inboard_peak"].get_scan(scan)
    r_b_tf_inboard_peak = mfile_data.data["r_b_tf_inboard_peak"].get_scan(scan)
    dx_tf_wp_insertion_gap = mfile_data.data["dx_tf_wp_insertion_gap"].get_scan(scan)
    r_tf_wp_inboard_outer = mfile_data.data["r_tf_wp_inboard_outer"].get_scan(scan)
    r_tf_wp_inboard_centre = mfile_data.data["r_tf_wp_inboard_centre"].get_scan(scan)

    if i_tf_turns_integer == 1:
        turn_layers = mfile_data.data["n_layer"].get_scan(scan)
        turn_pancakes = mfile_data.data["n_pancake"].get_scan(scan)

    # Superconducting coil check
    if i_tf_sup == 1:
        axis.add_patch(
            Circle(
                [0, 0],
                r_tf_inboard_in,
                facecolor="none",
                edgecolor="black",
                linestyle="--",
            ),
        )

        if i_tf_case_geom == 0:
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
        rad_tf_coil_inboard_toroidal_half = mfile_data.data[
            "rad_tf_coil_inboard_toroidal_half"
        ].get_scan(scan)

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
        if i_tf_case_geom == 0:
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

        elif i_tf_case_geom == 1:
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
        if i_tf_case_geom == 0:
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
        elif i_tf_case_geom == 1:
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
            if i_tf_turns_integer == 1:
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
                    (r_tf_wp_inboard_inner, -(0.5 * dx_tf_wp_primary_toroidal)),
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
            label=f"Peak Field: {b_tf_inboard_peak:.2f} T\n"
            f"r={r_b_tf_inboard_peak:.3f} m",
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
            f"$\\mathbf{{Casing:}}$\n \n"
            f"Coil half angle: {mfile_data.data['rad_tf_coil_inboard_toroidal_half'].get_scan(scan):.3f} radians\n\n"
            f"$\\text{{Full Coil Case:}}$\n"
            f"$r_{{start}} \\rightarrow r_{{end}}$: {mfile_data.data['r_tf_inboard_in'].get_scan(scan):.3f} $\\rightarrow$ {mfile_data.data['r_tf_inboard_out'].get_scan(scan):.3f} m\n"
            f"$\\Delta r$: {mfile_data.data['dr_tf_inboard'].get_scan(scan):.3f} m\n"
            f"Area of casing around WP: {mfile_data.data['a_tf_coil_inboard_case'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n\n"
            f"$\\text{{Nose Case:}}$\n"
            f"$r_{{start}} \\rightarrow r_{{end}}$: {mfile_data.data['r_tf_inboard_in'].get_scan(scan):.3f} $\\rightarrow$ {mfile_data.data['r_tf_wp_inboard_inner'].get_scan(scan):.3f} m\n"
            f"$\\Delta r$: {mfile_data.data['dr_tf_nose_case'].get_scan(scan):.3f} m\n"
            f"$A$: {mfile_data.data['a_tf_coil_nose_case'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n\n"
            f"$\\text{{Plasma Case:}}$\n"
            f"$r_{{start}} \\rightarrow r_{{end}}$: {mfile_data.data['r_tf_wp_inboard_outer'].get_scan(scan):.3f} $\\rightarrow$ {mfile_data.data['r_tf_inboard_out'].get_scan(scan):.3f} m\n"
            f"$\\Delta r$: {mfile_data.data['dr_tf_plasma_case'].get_scan(scan):.3f} m\n"
            f"$A$: {mfile_data.data['a_tf_plasma_case'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n\n"
            f"$\\text{{Side Case:}}$\n"
            f"Minimum $\\Delta r$: {mfile_data.data['dx_tf_side_case_min'].get_scan(scan):.3f} m\n"
            f"Average $\\Delta r$: {mfile_data.data['dx_tf_side_case_average'].get_scan(scan):.3f} m"
        )
        axis.text(
            0.775,
            0.925,
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
            f"$\\mathbf{{Insulation:}}$\n \n"
            f"Area of insulation around WP: {mfile_data.data['a_tf_wp_ground_insulation'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n"
            f"$\\Delta r$: {mfile_data.data['dx_tf_wp_insulation'].get_scan(scan):.4f} m\n\n"
            f"WP Insertion Gap:\n"
            f"$\\Delta r$: {mfile_data.data['dx_tf_wp_insertion_gap'].get_scan(scan):.4f} m"
        )
        axis.text(
            0.775,
            0.55,
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
            f"$\\mathbf{{Winding Pack:}}$\n \n"
            f"$N_{{\\text{{turns}}}}$: "
            f"{int(mfile_data.data['n_tf_coil_turns'].get_scan(scan))} turns\n"
            f"$r_{{start}} \\rightarrow r_{{end}}$: {mfile_data.data['r_tf_wp_inboard_inner'].get_scan(scan):.3f} $\\rightarrow$ {mfile_data.data['r_tf_wp_inboard_outer'].get_scan(scan):.3f} m\n"
            f"$\\Delta r$: {mfile_data.data['dr_tf_wp_with_insulation'].get_scan(scan):.3f} m\n"
            f"$A$, with insulation: {mfile_data.data['a_tf_wp_with_insulation'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n"
            f"$A$, no insulation: {mfile_data.data['a_tf_wp_no_insulation'].get_scan(scan):.3f} $\\mathrm{{m}}^2$\n\n"
            f"Primary WP:\n"
            f"$\\Delta x$: {mfile_data.data['dx_tf_wp_primary_toroidal'].get_scan(scan):.4f} m\n\n"
            f"Secondary WP:\n"
            f"$\\Delta x$: {mfile_data.data['dx_tf_wp_secondary_toroidal'].get_scan(scan):.4f} m\n\n"
            f"$J$ no insulation: {mfile_data.data['j_tf_wp'].get_scan(scan) / 1e6:.4f} MA/m$^2$"
        )
        axis.text(
            0.775,
            0.425,
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
            f"$\\mathbf{{General \\ info:}}$\n \n"
            f"Self inductance: {mfile_data.data['ind_tf_coil'].get_scan(scan) * 1e6:.4f} $\\mu$H\n"
            f"Stored energy of all coils: {mfile_data.data['e_tf_magnetic_stored_total_gj'].get_scan(scan):.4f} GJ\n"
        )
        axis.text(
            0.55,
            0.425,
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
        axis.set_xlim(0.0, r_tf_inboard_out * 1.1)
        axis.set_ylim((y14[-1] * 1.25), (-y14[-1] * 1.25))

        axis.set_title("Top-down view of inboard TF coil at midplane")
        axis.set_xlabel("Radial distance [m]")
        axis.set_ylabel("Toroidal distance [m]")
        axis.legend(loc="upper left")


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
    i_tf_turns_integer = mfile_data.data["i_tf_turns_integer"].get_scan(scan)
    # If integer turns switch is on then the turns can have non square dimensions
    if i_tf_turns_integer == 1:
        turn_width = round(mfile_data.data["dr_tf_turn"].get_scan(scan) * 1e3, 5)
        turn_height = round(mfile_data.data["dx_tf_turn"].get_scan(scan) * 1e3, 5)
        cable_space_width_radial = round(
            mfile_data.data["dr_tf_turn_cable_space"].get_scan(scan) * 1e3, 5
        )
        cable_space_width_toroidal = round(
            mfile_data.data["dx_tf_turn_cable_space"].get_scan(scan) * 1e3, 5
        )
    elif i_tf_turns_integer == 0:
        turn_width = round(mfile_data.data["t_turn_tf"].get_scan(scan) * 1e3, 5)
        cable_space_width = round(
            mfile_data.data["dx_tf_turn_cable_space_average"].get_scan(scan) * 1e3, 5
        )

    he_pipe_diameter = round(
        mfile_data.data["dia_tf_turn_coolant_channel"].get_scan(scan) * 1e3, 5
    )
    steel_thickness = round(mfile_data.data["dx_tf_turn_steel"].get_scan(scan) * 1e3, 5)
    insulation_thickness = round(
        mfile_data.data["dx_tf_turn_insulation"].get_scan(scan) * 1e3, 5
    )
    internal_cable_space = round(
        mfile_data.data["a_tf_turn_cable_space_no_void"].get_scan(scan) * 1e6, 5
    )
    c_tf_turn = mfile_data.data["c_tf_turn"].get_scan(scan)

    # Plot the total turn shape
    if i_tf_turns_integer == 0:
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
                label=f"Cooling pipe: \n{he_pipe_diameter} mm diameter \n \n Current per turn: {c_tf_turn:.2f} A",
                edgecolor="black",
            ),
        )
        axis.set_xlim(-turn_width * 0.05, turn_width * 1.05)
        axis.set_ylim(-turn_width * 0.05, turn_width * 1.05)

    # Non square turns
    elif i_tf_turns_integer == 1:
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
                label=f"Cooling pipe: \n{he_pipe_diameter} mm diameter \n \n Current per turn: {c_tf_turn:.2f} A",
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
        if "r_pf_coil_middle[" in item:
            number_of_coils += 1

    # Check for Central Solenoid
    iohcl = mfile_data.data["iohcl"].get_scan(scan) if "iohcl" in mfile_data.data else 1

    # If Central Solenoid present, ignore last entry in for loop
    # The last entry will be the OH coil in this case
    noc = number_of_coils - 1 if iohcl == 1 else number_of_coils

    for coil in range(noc):
        coils_r.append(mfile_data.data[f"r_pf_coil_middle[{coil:01}]"].get_scan(scan))
        coils_z.append(mfile_data.data[f"z_pf_coil_middle[{coil:01}]"].get_scan(scan))
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
            f"!{OBJECTIVE_NAMES[abs(int(mfile_data.data['minmax'].get_scan(-1)))]}",
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

    nong = mfile_data.data["nd_electron_line"].get_scan(scan) / mfile_data.data[
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
    if "p_l_h_threshold_mw" in mfile_data.data:
        pthresh = mfile_data.data["p_l_h_threshold_mw"].get_scan(scan)
    else:
        pthresh = mfile_data.data["l_h_threshold_powers(6)"].get_scan(scan)

    data = [
        ("p_fusion_total_mw", "Fusion power", "MW"),
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
        if "r_pf_coil_middle[" in item:
            number_of_coils += 1

    pf_info = [
        (
            mfile_data.data[f"c_pf_cs_coils_peak_ma[{i:01}]"].get_scan(scan),
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
        tftype = SUPERCONDUCTING_TF_TYPES[
            int(mfile_data.data["i_tf_sc_mat"].get_scan(scan))
        ]
    else:
        tftype = "Resistive Copper"

    vssoft = mfile_data.data["vs_plasma_res_ramp"].get_scan(scan) + mfile_data.data[
        "vs_plasma_ind_ramp"
    ].get_scan(scan)

    sig_case = 1.0e-6 * mfile_data.data[f"s_shear_tf_peak({i_tf_bucking})"].get_scan(
        scan
    )
    sig_cond = 1.0e-6 * mfile_data.data[
        f"s_shear_tf_peak({i_tf_bucking + 1})"
    ].get_scan(scan)

    if i_tf_sup == 1:
        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vs_cs_pf_total_pulse", "Available flux swing", "Wb"),
            (t_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            ("bmaxtfrp", "Peak field at conductor (w. rip.)", "T"),
            ("iooic", r"I/I$_{\mathrm{crit}}$", ""),
            ("tmargtf", "TF Temperature margin", "K"),
            ("temp_cs_margin", "CS Temperature margin", "K"),
            (sig_cond, "TF Cond max TRESCA stress", "MPa"),
            (sig_case, "TF Case max TRESCA stress", "MPa"),
            ("m_tf_coils_total/n_tf_coils", "Mass per TF coil", "kg"),
        ]

    else:
        p_cp_resistive = 1.0e-6 * mfile_data.data["p_cp_resistive"].get_scan(scan)
        p_tf_leg_resistive = 1.0e-6 * mfile_data.data["p_tf_leg_resistive"].get_scan(
            scan
        )
        p_tf_joints_resistive = 1.0e-6 * mfile_data.data[
            "p_tf_joints_resistive"
        ].get_scan(scan)
        fcoolcp = 100.0 * mfile_data.data["fcoolcp"].get_scan(scan)

        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vs_cs_pf_total_pulse", "Available flux swing", "Wb"),
            (t_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            ("b_tf_inboard_peak", "Peak field at conductor (w. rip.)", "T"),
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
            (p_tf_joints_resistive, "TF joints resisitive heating ", "MW"),
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
        mfile_data.data["p_plant_electric_gross_mw"].get_scan(scan)
        / mfile_data.data["p_plant_primary_heat_mw"].get_scan(scan)
    )

    net_eff = 100.0 * (
        (
            mfile_data.data["p_plant_electric_gross_mw"].get_scan(scan)
            - mfile_data.data["p_coolant_pump_elec_total_mw"].get_scan(scan)
        )
        / (
            mfile_data.data["p_plant_primary_heat_mw"].get_scan(scan)
            - mfile_data.data["p_coolant_pump_elec_total_mw"].get_scan(scan)
        )
    )

    plant_eff = 100.0 * (
        mfile_data.data["p_plant_electric_net_mw"].get_scan(scan)
        / mfile_data.data["p_fusion_total_mw"].get_scan(scan)
    )

    # Define appropriate pedestal and impurity parameters
    coredescription = (
        "radius_plasma_core_norm",
        "Normalised radius of 'core' region",
        "",
    )
    if ipedestal == 1:
        ped_height = ("neped", "Electron density at pedestal", "m$^{-3}$")
        ped_pos = ("rhopedn", "r/a at density pedestal", "")
    else:
        ped_height = ("", "No pedestal model used", "")
        ped_pos = ("", "", "")

    p_cryo_plant_electric_mw = mfile_data.data["p_cryo_plant_electric_mw"].get_scan(
        scan
    )

    data = [
        ("pflux_fw_neutron_mw", "Nominal neutron wall load", "MW m$^{-2}$"),
        coredescription,
        ped_height,
        ped_pos,
        ("p_plasma_inner_rad_mw", "Inner zone radiation", "MW"),
        ("p_plasma_rad_mw", "Total radiation in LCFS", "MW"),
        ("p_blkt_nuclear_heat_total_mw", "Nuclear heating in blanket", "MW"),
        ("p_shld_nuclear_heat_mw", "Nuclear heating in shield", "MW"),
        (p_cryo_plant_electric_mw, "TF cryogenic power", "MW"),
        ("p_plasma_separatrix_mw", "Power to divertor", "MW"),
        ("divlife", "Divertor life", "years"),
        ("p_plant_primary_heat_mw", "Primary (high grade) heat", "MW"),
        (gross_eff, "Gross cycle efficiency", "%"),
        (net_eff, "Net cycle efficiency", "%"),
        ("p_plant_electric_gross_mw", "Gross electric power", "MW"),
        ("p_plant_electric_net_mw", "Net electric power", "MW"),
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
    i_hcd_primary = mfile_data.data["i_hcd_primary"].get_scan(scan)
    nbi = False
    ecrh = False
    ebw = False
    lhcd = False
    iccd = False

    if (i_hcd_primary == 5) or (i_hcd_primary == 8):
        nbi = True
        axis.text(-0.05, 1, "Neutral Beam Current Drive:", ha="left", va="center")
    if (
        (i_hcd_primary == 3)
        or (i_hcd_primary == 7)
        or (i_hcd_primary == 10)
        or (i_hcd_primary == 11)
        or (i_hcd_primary == 13)
    ):
        ecrh = True
        axis.text(-0.05, 1, "Electron Cyclotron Current Drive:", ha="left", va="center")
    if i_hcd_primary == 12:
        ebw = True
        axis.text(-0.05, 1, "Electron Bernstein Wave Drive:", ha="left", va="center")
    if i_hcd_primary in [1, 4, 6]:
        lhcd = True
        axis.text(
            -0.05,
            1,
            "Lower Hybrid Current Drive:",
            ha="left",
            va="center",
        )
    if i_hcd_primary == 2:
        iccd = True
        axis.text(-0.05, 1, "Ion Cyclotron Current Drive:", ha="left", va="center")

    if "i_hcd_secondary" in mfile_data.data:
        secondary_heating = ""
        i_hcd_secondary = mfile_data.data["i_hcd_secondary"].get_scan(scan)

        if (i_hcd_secondary == 5) or (i_hcd_secondary == 8):
            secondary_heating = "NBI"
        if (
            (i_hcd_secondary == 3)
            or (i_hcd_secondary == 7)
            or (i_hcd_secondary == 10)
            or (i_hcd_secondary == 11)
            or (i_hcd_secondary == 13)
        ):
            secondary_heating = "ECH"
        if i_hcd_secondary == 12:
            secondary_heating = "EBW"
        if i_hcd_secondary in [1, 4, 6]:
            secondary_heating = "LHCD"
        if i_hcd_secondary == 2:
            secondary_heating = "ICCD"

    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    pinjie = mfile_data.data["p_hcd_injected_total_mw"].get_scan(scan)
    p_plasma_separatrix_mw = mfile_data.data["p_plasma_separatrix_mw"].get_scan(scan)
    pdivr = p_plasma_separatrix_mw / mfile_data.data["rmajor"].get_scan(scan)

    if mfile_data.data["i_hcd_secondary"].get_scan(scan) != 0:
        pinjmwfix = mfile_data.data["pinjmwfix"].get_scan(scan)

    pdivnr = (
        1.0e20
        * mfile_data.data["p_plasma_separatrix_mw"].get_scan(scan)
        / (
            mfile_data.data["rmajor"].get_scan(scan)
            * mfile_data.data["dene"].get_scan(scan)
        )
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    if "p_l_h_threshold_mw" in mfile_data.data:
        pthresh = mfile_data.data["p_l_h_threshold_mw"].get_scan(scan)
    else:
        pthresh = mfile_data.data["l_h_threshold_powers(6)"].get_scan(scan)
    flh = p_plasma_separatrix_mw / pthresh

    hstar = mfile_data.data["hstar"].get_scan(scan)

    if ecrh:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
            ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
            ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
            ("f_c_plasma_inductive", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "eta_cd_hcd_primary",
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
        # i_hcd_secondary is now always in the MFILE with = 0 meaning no fixed heating
        if mfile_data.data["i_hcd_secondary"].get_scan(scan) != 0:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if nbi:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
            ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
            ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
            ("f_c_plasma_inductive", "Inductive fraction", ""),
            ("gamnb", "NB gamma", "$10^{20}$ A W$^{-1}$ m$^{-2}$"),
            ("e_beam_kev", "NB energy", "keV"),
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
        if mfile_data.data["i_hcd_secondary"].get_scan(scan) != 0:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if ebw:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
            ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
            ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
            ("f_c_plasma_inductive", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "eta_cd_norm_hcd_primary",
                "Normalised current drive efficiency of primary HCD system",
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
        if "i_hcd_secondary" in mfile_data.data:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if lhcd:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
            ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
            ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
            ("f_c_plasma_inductive", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "eta_cd_norm_hcd_primary",
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
        if "i_hcd_secondary" in mfile_data.data:
            data.insert(
                1, ("pinjmwfix", f"{secondary_heating} secondary auxiliary power", "MW")
            )
            data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
            data.insert(2, (pinjie, "Total auxillary power", "MW"))

    if iccd:
        data = [
            (pinjie, "Steady state auxiliary power", "MW"),
            ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
            ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
            ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
            ("f_c_plasma_inductive", "Inductive fraction", ""),
            ("p_plasma_loss_mw", "Plasma heating used for H factor", "MW"),
            (
                "eta_cd_norm_hcd_primary",
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
        if "i_hcd_secondary" in mfile_data.data:
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

    boot_ipdg = mfile_data.data["f_c_plasma_bootstrap_iter89"].get_scan(scan)
    boot_sauter = mfile_data.data["f_c_plasma_bootstrap_sauter"].get_scan(scan)
    boot_nenins = mfile_data.data["f_c_plasma_bootstrap_nevins"].get_scan(scan)
    boot_wilson = mfile_data.data["f_c_plasma_bootstrap_wilson"].get_scan(scan)
    boot_sakai = mfile_data.data["f_c_plasma_bootstrap_sakai"].get_scan(scan)
    boot_aries = mfile_data.data["f_c_plasma_bootstrap_aries"].get_scan(scan)
    boot_andrade = mfile_data.data["f_c_plasma_bootstrap_andrade"].get_scan(scan)
    boot_hoang = mfile_data.data["f_c_plasma_bootstrap_hoang"].get_scan(scan)
    boot_wong = mfile_data.data["f_c_plasma_bootstrap_wong"].get_scan(scan)
    boot_gi_I = mfile_data.data["bscf_gi_i"].get_scan(scan)  # noqa: N806
    boot_gi_II = mfile_data.data["bscf_gi_ii"].get_scan(scan)  # noqa: N806
    boot_sugiyama_l = mfile_data.data["f_c_plasma_bootstrap_sugiyama_l"].get_scan(scan)
    boot_sugiyama_h = mfile_data.data["f_c_plasma_bootstrap_sugiyama_h"].get_scan(scan)

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
        "Sugiyama (L-mode)": boot_sugiyama_l,
        "Sugiyama (H-mode)": boot_sugiyama_h,
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
    iter_nominal = mfile_data.data["l_h_threshold_powers(1)"].get_scan(scan)
    iter_upper = mfile_data.data["l_h_threshold_powers(2)"].get_scan(scan)
    iter_lower = mfile_data.data["l_h_threshold_powers(3)"].get_scan(scan)
    iter_1997_1 = mfile_data.data["l_h_threshold_powers(4)"].get_scan(scan)
    iter_1997_2 = mfile_data.data["l_h_threshold_powers(5)"].get_scan(scan)
    martin_nominal = mfile_data.data["l_h_threshold_powers(6)"].get_scan(scan)
    martin_upper = mfile_data.data["l_h_threshold_powers(7)"].get_scan(scan)
    martin_lower = mfile_data.data["l_h_threshold_powers(8)"].get_scan(scan)
    snipes_nominal = mfile_data.data["l_h_threshold_powers(9)"].get_scan(scan)
    snipes_upper = mfile_data.data["l_h_threshold_powers(10)"].get_scan(scan)
    snipes_lower = mfile_data.data["l_h_threshold_powers(11)"].get_scan(scan)
    snipes_closed_nominal = mfile_data.data["l_h_threshold_powers(12)"].get_scan(scan)
    snipes_closed_upper = mfile_data.data["l_h_threshold_powers(13)"].get_scan(scan)
    snipes_closed_lower = mfile_data.data["l_h_threshold_powers(14)"].get_scan(scan)
    hubbard_nominal = mfile_data.data["l_h_threshold_powers(15)"].get_scan(scan)
    hubbard_lower = mfile_data.data["l_h_threshold_powers(16)"].get_scan(scan)
    hubbard_upper = mfile_data.data["l_h_threshold_powers(17)"].get_scan(scan)
    hubbard_2017 = mfile_data.data["l_h_threshold_powers(18)"].get_scan(scan)
    martin_aspect_nominal = mfile_data.data["l_h_threshold_powers(19)"].get_scan(scan)
    martin_aspect_upper = mfile_data.data["l_h_threshold_powers(20)"].get_scan(scan)
    martin_aspect_lower = mfile_data.data["l_h_threshold_powers(21)"].get_scan(scan)

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


def plot_confinement_time_comparison(
    axis: plt.Axes, mfile_data: mf.MFile, scan: int, u_seed=None
) -> None:
    """
    Function to plot a scatter box plot of confinement time comparisons.

    Arguments:
        axis (plt.Axes): Axis object to plot to.
        mfile_data (mf.MFile): MFILE data object.
        scan (int): Scan number to use.
    """
    rminor = mfile_data.data["rminor"].get_scan(scan)
    rmajor = mfile_data.data["rmajor"].get_scan(scan)
    c_plasma_ma = mfile_data.data["plasma_current_ma"].get_scan(scan)
    kappa95 = mfile_data.data["kappa95"].get_scan(scan)
    dnla20 = mfile_data.data["nd_electron_line"].get_scan(scan) / 1e20
    afuel = mfile_data.data["m_fuel_amu"].get_scan(scan)
    bt = mfile_data.data["bt"].get_scan(scan)
    p_plasma_separatrix_mw = mfile_data.data["p_plasma_separatrix_mw"].get_scan(scan)
    kappa = mfile_data.data["kappa"].get_scan(scan)
    aspect = mfile_data.data["aspect"].get_scan(scan)
    dnla19 = mfile_data.data["nd_electron_line"].get_scan(scan) / 1e19
    kappa_ipb = mfile_data.data["kappa_ipb"].get_scan(scan)
    triang = mfile_data.data["triang"].get_scan(scan)
    m_ions_total_amu = mfile_data.data["m_ions_total_amu"].get_scan(scan)

    # Calculate confinement times using the scan data
    iter_89p = confine.iter_89p_confinement_time(
        pcur=c_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        dnla20=dnla20,
        bt=bt,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_89_0 = confine.iter_89_0_confinement_time(
        pcur=c_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        dnla20=dnla20,
        bt=bt,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_h90_p = confine.iter_h90_p_confinement_time(
        pcur=c_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        dnla20=dnla20,
        bt=bt,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_h90_p_amended = confine.iter_h90_p_amended_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        afuel=afuel,
        rmajor=rmajor,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        kappa=kappa,
    )
    iter_93h = confine.iter_93h_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        afuel=afuel,
        rmajor=rmajor,
        dnla20=dnla20,
        aspect=aspect,
        kappa=kappa,
    )
    iter_h97p = confine.iter_h97p_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        dnla19=dnla19,
        rmajor=rmajor,
        aspect=aspect,
        kappa=kappa,
        afuel=afuel,
    )
    iter_h97p_elmy = confine.iter_h97p_elmy_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        dnla19=dnla19,
        rmajor=rmajor,
        aspect=aspect,
        kappa=kappa,
        afuel=afuel,
    )
    iter_96p = confine.iter_96p_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        kappa95=kappa95,
        rmajor=rmajor,
        aspect=aspect,
        dnla19=dnla19,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_pb98py = confine.iter_pb98py_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa=kappa,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y = confine.iter_ipb98y_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa=kappa,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y1 = confine.iter_ipb98y1_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y2 = confine.iter_ipb98y2_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y3 = confine.iter_ipb98y3_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y4 = confine.iter_ipb98y4_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    petty08 = confine.petty08_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
    )
    menard_nstx = confine.menard_nstx_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    menard_nstx_petty08 = confine.menard_nstx_petty08_hybrid_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    itpa20 = confine.itpa20_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        dnla19=dnla19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        triang=triang,
        kappa_ipb=kappa_ipb,
        eps=(1 / aspect),
        aion=m_ions_total_amu,
    )
    itpa20_ilc = confine.itpa20_il_confinement_time(
        pcur=c_plasma_ma,
        bt=bt,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        dnla19=dnla19,
        aion=m_ions_total_amu,
        rmajor=rmajor,
        triang=triang,
        kappa_ipb=kappa_ipb,
    )

    # Data for the box plot
    data = {
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[6])}": iter_89p,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[7])}": iter_89_0,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[13])}": iter_h90_p,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[20])}": iter_h90_p_amended,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[24])}": iter_93h,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[26])}": iter_h97p,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[27])}": iter_h97p_elmy,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[28])}": iter_96p,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[31])}": iter_pb98py,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[32])}": iter_ipb98y,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[33])}": iter_ipb98y1,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[34])}": iter_ipb98y2,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[35])}": iter_ipb98y3,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[36])}": iter_ipb98y4,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[41])}": petty08,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[46])}": menard_nstx,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[47])}": menard_nstx_petty08,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[49])}": itpa20,
        rf"{f2py_compatible_to_string(physics_variables.labels_confinement_scalings[50])}": itpa20_ilc,
    }

    # Create the violin plot
    axis.violinplot(data.values(), showextrema=False)

    # Create the box plot
    axis.boxplot(
        data.values(), showfliers=True, showmeans=True, meanline=True, widths=0.3
    )

    # Scatter plot for each data point
    # Use a set of distinct colors for better differentiation
    distinct_colors = [
        "#1f77b4",  # blue
        "#ff7f0e",  # orange
        "#2ca02c",  # green
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#bcbd22",  # olive
        "#17becf",  # cyan
        "#aec7e8",  # light blue
        "#ffbb78",  # light orange
        "#98df8a",  # light green
        "#ff9896",  # light red
        "#c5b0d5",  # light purple
        "#c49c94",  # light brown
        "#f7b6d2",  # light pink
        "#c7c7c7",  # light gray
        "#dbdb8d",  # light olive
        "#9edae5",  # light cyan
    ]
    generator = np.random.default_rng(seed=u_seed)
    x_values = generator.normal(loc=1, scale=0.035, size=len(data.values()))
    for index, (key, value) in enumerate(data.items()):
        if "Hubbard" in key and "2017" not in key:
            color = "#800080"  # strong purple
        else:
            color = distinct_colors[index % len(distinct_colors)]
        axis.scatter(
            x_values[index],
            value,
            color=color,
            label=key,
            alpha=1.0,
            edgecolor="black",
            linewidth=0.7,
        )
    axis.legend(loc="upper left", bbox_to_anchor=(-1.3, 0.75), ncol=2)

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_threshold = np.mean(data_values)
    std_threshold = np.std(data_values)
    median_threshold = np.median(data_values)

    # Plot average, standard deviation, and median as text
    axis.text(
        0.7,
        1.25,
        f"Average: {avg_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        0.7,
        1.2,
        f"Standard Dev: {std_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        0.7,
        1.15,
        f"Median: {median_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        0.75,
        -0.05,
        r"$H \ factor = 1.0$",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Confinement time ($\\tau_{\\text{E}}$) Comparison")
    axis.set_ylabel("Confinement time, $\\tau_{\\text{E}}$ [s]")
    axis.set_xlim([0.5, 1.5])
    axis.set_xticks([])
    axis.set_xticklabels([])

    # Add background color
    axis.set_facecolor("#f0f0f0")


def plot_radial_build(axis: plt.Axes, mfile_data: mf.MFile, colour_scheme) -> None:
    """
    Plots the radial build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's radial build
    (such as central solenoid, toroidal field coils, vacuum vessel, shields, blankets, etc.)
    as a horizontal stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile_data`, and each segment is color-coded and labeled accordingly.

    If the toroidal field coil is inside the central solenoid (as indicated by the
    "i_tf_inside_cs" flag in `mfile_data`), the order and labels of the components are
    adjusted accordingly.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The matplotlib axis on which to plot the radial build.
    mfile_data : mf.MFile
        An object containing the machine build data, with required fields for each
        radial component and the "i_tf_inside_cs" flag.

    Returns
    -------
    None
        This function modifies the provided axis in-place and does not return a value.

    Notes
    -----
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
    if int(mfile_data.data["i_tf_inside_cs"].get_scan(-1)) == 1:
        radial_variables[1] = "dr_tf_inboard"
        radial_variables[2] = "dr_cs_tf_gap"
        radial_variables[3] = "dr_cs"
        radial_variables[4] = "dr_cs_precomp"
        radial_variables[5] = "dr_tf_shld_gap"

    radial_build = [[mfile_data.data[rl].get_scan(-1) for rl in radial_variables]]

    radial_build = np.array(radial_build)

    for kk in range(radial_build.shape[0]):
        radial_build[kk, 14] = 2.0 * radial_build[kk, 14]

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
    if int(mfile_data.data["i_tf_inside_cs"].get_scan(-1)) == 1:
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
        TFC_COLOUR[colour_scheme - 1],
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
        TFC_COLOUR[colour_scheme - 1],
    ]
    if int(mfile_data.data["i_tf_inside_cs"].get_scan(-1)) == 1:
        radial_color[1] = TFC_COLOUR[colour_scheme - 1]
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
            label=f"{radial_labels[kk]}\n[{radial_variables[kk]}]\n{radial_build[kk][0]:.3f} m",
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
        mfile_data.data["rmajor"].get_scan(-1),
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Major Radius $R_0$",
    )
    axis.minorticks_on()
    axis.set_xlabel("Radius [m]")


def plot_lower_vertical_build(
    axis: plt.Axes, mfile_data: mf.MFile, colour_scheme
) -> None:
    """
    Plots the lower vertical build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's vertical build
    (such as plasma, first wall, divertor, shield, vacuum vessel, thermal shield, TF coil, etc.)
    as a vertical stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile_data`, and each segment is color-coded and labeled accordingly.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The matplotlib axis on which to plot the vertical build.
    mfile_data : mf.MFile
        An object containing the machine build data, with required fields for each
        vertical component.
    colour_scheme : int
        Colour scheme index to use for component colors.

    Returns
    -------
    None
        This function modifies the provided axis in-place and does not return a value.

    Notes
    -----
    - Components with zero thickness are omitted from the plot.
    - The legend displays the name and thickness (in meters) of each component.
    """

    lower_vertical_variables = [
        "z_plasma_xpoint_upper",
        "dz_xpoint_divertor",
        "dz_divertor",
        "dz_shld_upper",
        "dz_vv_upper",
        "dz_shld_vv_gap",
        "dz_shld_thermal",
        "dr_tf_shld_gap",
        "dr_tf_inboard",
        "dz_tf_cryostat",
    ]

    lower_vertical_build = [
        [mfile_data.data[rl].get_scan(-1) for rl in lower_vertical_variables]
    ]

    lower_vertical_build = np.array(lower_vertical_build)

    lower_vertical_build = np.transpose(lower_vertical_build)

    lower_vertical_labels = [
        "Plasma Height",
        "Plasma - Divertor Gap",
        "Divertor",
        "Shield",
        "Vacuum Vessel",
        "Shield - VV Gap",
        "Thermal shield",
        "TF Coil - Shield Gap",
        "TF Coil",
        "TF Coil - Cryostat gap",
    ]

    lower_vertical_color = [
        PLASMA_COLOUR[colour_scheme - 1],
        "white",
        "black",
        SHIELD_COLOUR[colour_scheme - 1],
        VESSEL_COLOUR[colour_scheme - 1],
        "white",
        THERMAL_SHIELD_COLOUR[colour_scheme - 1],
        "white",
        TFC_COLOUR[colour_scheme - 1],
        "white",
    ]

    # Remove build parts equal to zero
    mask = ~(lower_vertical_build[:, 0] == 0.0)
    filtered_vertical_build = lower_vertical_build[mask]
    filtered_labels = [lbl for i, lbl in enumerate(lower_vertical_labels) if mask[i]]
    filtered_colors = [col for i, col in enumerate(lower_vertical_color) if mask[i]]

    bottom = np.zeros(filtered_vertical_build.shape[1])
    for kk in range(filtered_vertical_build.shape[0]):
        axis.bar(
            np.arange(filtered_vertical_build.shape[1]),
            -filtered_vertical_build[kk, :],
            bottom=bottom,
            width=0.8,
            label=f"{filtered_labels[kk]}\n[{lower_vertical_variables[kk]}]\n{filtered_vertical_build[kk][0]:.3f} m",
            color=filtered_colors[kk],
            edgecolor="black",
            linewidth=0.05,
        )
        bottom -= filtered_vertical_build[kk, :]

    axis.set_xticks([])
    axis.legend(
        bbox_to_anchor=(0, 0),
        loc="upper left",
        ncol=5,
    )
    axis.minorticks_on()
    axis.set_ylabel("Height [m]")
    axis.title.set_text("Lower Vertical Build")


def plot_upper_vertical_build(
    axis: plt.Axes, mfile_data: mf.MFile, colour_scheme
) -> None:
    """
    Plots the upper vertical build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's vertical build
    (such as plasma, first wall, divertor, shield, vacuum vessel, thermal shield, TF coil, etc.)
    as a vertical stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile_data`, and each segment is color-coded and labeled accordingly.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        The matplotlib axis on which to plot the vertical build.
    mfile_data : mf.MFile
        An object containing the machine build data, with required fields for each
        vertical component.
    colour_scheme : int
        Colour scheme index to use for component colors.

    Returns
    -------
    None
        This function modifies the provided axis in-place and does not return a value.

    Notes
    -----
    - Components with zero thickness are omitted from the plot.
    - The legend displays the name and thickness (in meters) of each component.
    """
    if mfile_data.data["i_single_null"].get_scan(-1) == 1:
        upper_vertical_variables = [
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_fw_upper",
            "dz_blkt_upper",
            "dr_shld_blkt_gap",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
            "dz_tf_cryostat",
        ]
        upper_vertical_labels = [
            "Plasma Height",
            "First Wall - Plasma Gap",
            "First Wall Upper",
            "Blanket Upper",
            "Shield-Blanket Gap",
            "Shield Upper",
            "Vacuum Vessel Upper",
            "Shield-VV Gap",
            "Thermal Shield",
            "TF Coil - Shield Gap",
            "TF Coil",
            "TF Coil - Cryostat gap",
        ]
        upper_vertical_colours = [
            PLASMA_COLOUR[colour_scheme - 1],
            "white",
            FIRSTWALL_COLOUR[colour_scheme - 1],
            BLANKET_COLOUR[colour_scheme - 1],
            "white",
            SHIELD_COLOUR[colour_scheme - 1],
            VESSEL_COLOUR[colour_scheme - 1],
            "white",
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
            "white",
            TFC_COLOUR[colour_scheme - 1],
            "white",
        ]
    # Double null case
    else:
        upper_vertical_variables = [
            "z_plasma_xpoint_upper",
            "dz_xpoint_divertor",
            "dz_divertor",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
            "dz_tf_cryostat",
        ]
        upper_vertical_labels = [
            "Plasma Height",
            "Plasma - Divertor Gap",
            "Divertor Upper",
            "Shield Upper",
            "Vacuum Vessel Upper",
            "Shield-VV Gap",
            "Thermal Shield",
            "TF Coil - Shield Gap",
            "TF Coil",
            "TF Coil - Cryostat gap",
        ]
        upper_vertical_colours = [
            PLASMA_COLOUR[colour_scheme - 1],
            "white",
            "black",
            SHIELD_COLOUR[colour_scheme - 1],
            VESSEL_COLOUR[colour_scheme - 1],
            "white",
            THERMAL_SHIELD_COLOUR[colour_scheme - 1],
            "white",
            TFC_COLOUR[colour_scheme - 1],
            "white",
        ]

    # Get thicknesses for each layer
    upper_vertical_build = np.array([
        mfile_data.data[rl].get_scan(-1) for rl in upper_vertical_variables
    ])

    # Remove build parts equal to zero
    mask = ~(upper_vertical_build == 0.0)
    filtered_build = upper_vertical_build[mask]
    filtered_labels = [lbl for i, lbl in enumerate(upper_vertical_labels) if mask[i]]
    filtered_colors = [col for i, col in enumerate(upper_vertical_colours) if mask[i]]
    filtered_vars = [v for i, v in enumerate(upper_vertical_variables) if mask[i]]

    # Compute cumulative positions (bottoms) for stacking
    bottoms = np.zeros_like(filtered_build)
    for i in range(1, len(filtered_build)):
        bottoms[i] = bottoms[i - 1] + filtered_build[i - 1]

    # Plot each layer as a bar, stacking upwards from zero
    for kk in range(len(filtered_build)):
        axis.bar(
            0,
            filtered_build[kk],
            bottom=bottoms[kk],
            width=0.8,
            label=f"{filtered_labels[kk]}\n[{filtered_vars[kk]}]\n{filtered_build[kk]:.3f} m",
            color=filtered_colors[kk],
            edgecolor="black",
            linewidth=0.05,
        )

    axis.set_xticks([])
    axis.legend(
        bbox_to_anchor=(0, 0),
        loc="upper left",
        ncol=6,
    )
    axis.minorticks_on()
    axis.set_ylabel("Height [m]")
    axis.title.set_text("Upper Vertical Build")


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


def plot_cs_coil_structure(axis, fig, mfile_data, scan, colour_scheme=1):
    """Function to plot the coil structure of the CS.

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use
        demo_ranges --> whether to use demo ranges for the plot
        colour_scheme --> colour scheme to use for the plot
    """
    # Get CS coil parameters
    dr_cs = mfile_data.data["dr_cs"].get_scan(scan)
    dz_cs = mfile_data.data["ohdz"].get_scan(scan)
    dr_bore = mfile_data.data["dr_bore"].get_scan(scan)

    # Plot the right side of the CS
    right_cs = patches.Rectangle(
        (dr_bore, -dz_cs / 2),
        dr_cs,
        dz_cs,
        edgecolor="black",
        facecolor=SOLENOID_COLOUR[colour_scheme - 1],
        lw=1.5,
        label="CS Coil",
    )
    axis.add_patch(right_cs)

    # Plot the bore of the machine
    bore_rect = patches.Rectangle(
        (-dr_bore, -dz_cs / 2),
        dr_bore * 2,
        dz_cs,
        edgecolor="black",
        facecolor="lightgrey",
        lw=1.0,
        label="Bore",
    )
    axis.add_patch(bore_rect)

    left_cs = patches.Rectangle(
        (-dr_bore - dr_cs, -dz_cs / 2),
        dr_cs,
        dz_cs,
        edgecolor="black",
        facecolor=SOLENOID_COLOUR[colour_scheme - 1],
        lw=1.5,
        label="CS Coil",
    )
    axis.add_patch(left_cs)

    # Draw vertical lines to represent CS turns
    # Get the turn width (radial thickness of each turn)
    dr_cs_turn = mfile_data.data["dr_cs_turn"].get_scan(scan)
    dz_cs_turn = mfile_data.data["dz_cs_turn"].get_scan(scan)
    # Number of vertical lines (number of turns)
    if dr_cs_turn > 0:
        n_lines = int(dr_cs / dr_cs_turn)
        for i in range(1, n_lines):
            x = dr_bore + i * dr_cs_turn
            axis.plot(
                [x, x],
                [-dz_cs / 2, dz_cs / 2],
                color="black",
                linestyle="--",
                linewidth=0.2,
            )
            x_left = -dr_bore - dr_cs + i * dr_cs_turn
            axis.plot(
                [x_left, x_left],
                [-dz_cs / 2, dz_cs / 2],
                color="black",
                linestyle="--",
                linewidth=0.2,
            )
    # Plot horizontal lines (along Z) for each turn
    if dz_cs_turn > 0:
        n_hlines = int(dz_cs / dz_cs_turn)
        for j in range(1, n_hlines):
            y = -dz_cs / 2 + j * dz_cs_turn
            # Right CS
            axis.plot(
                [dr_bore, dr_bore + dr_cs],
                [y, y],
                color="black",
                linestyle="--",
                linewidth=0.2,
            )
            # Left CS
            axis.plot(
                [-dr_bore - dr_cs, -dr_bore],
                [y, y],
                color="black",
                linestyle="--",
                linewidth=0.2,
            )

    textstr_cs = (
        f"$\\mathbf{{Coil \\ parameters:}}$\n \n"
        f"CS height vs TF internal height: {mfile_data.data['f_z_cs_tf_internal'].get_scan(scan):.2f}\n"
        f"$N_{{\\text{{turns}}}}:$ {mfile_data.data['n_pf_coil_turns[n_cs_pf_coils-1]'].get_scan(scan):.2f}\n"
        f"$I_{{\\text{{peak}}}}:$ {mfile_data.data['c_pf_cs_coils_peak_ma[n_cs_pf_coils-1]'].get_scan(scan):.3f}$ \\ MA$\n"
        f"$B_{{\\text{{peak}}}}:$ {mfile_data.data['b_pf_coil_peak[n_cs_pf_coils-1]'].get_scan(scan):.3f}$ \\ T$\n"
    )

    axis.text(
        0.6,
        0.825,
        textstr_cs,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("Central Solenoid Cross-Section")
    axis.grid(True, linestyle="--", alpha=0.3)
    axis.set_xlim(-(dr_bore + dr_cs + 0.1), (dr_bore + dr_cs + 0.1))
    axis.set_ylim(-dz_cs / 2 - 0.1, dz_cs / 2 + 0.1)


def plot_cs_turn_structure(axis, fig, mfile_data, scan):
    a_cs_turn = mfile_data.data["a_cs_turn"].get_scan(scan)
    dz_cs_turn = mfile_data.data["dz_cs_turn"].get_scan(scan)
    dr_cs_turn = mfile_data.data["dr_cs_turn"].get_scan(scan)

    ld_ratio_cst = mfile_data.data["ld_ratio_cst"].get_scan(scan)
    radius_cs_turn_cable_space = mfile_data.data["radius_cs_turn_cable_space"].get_scan(
        scan
    )
    t_structural_vertical = mfile_data.data["t_structural_vertical"].get_scan(scan)
    t_structural_radial = mfile_data.data["t_structural_radial"].get_scan(scan)
    r_out_cst = mfile_data.data["r_out_cst"].get_scan(scan)

    # Plot the CS turn as a rectangle representing the conductor cross-section
    # Assume dz_cs_turn is the diameter and dr_cs_turn is the length of the conductor cross-section

    # Draw the conductor cross-section as a rectangle
    axis.add_patch(
        patches.FancyBboxPatch(
            (0, 0),
            dr_cs_turn,
            dz_cs_turn,
            boxstyle=patches.BoxStyle("Round", pad=0, rounding_size=r_out_cst),
            edgecolor="black",
            facecolor="grey",
            lw=1.5,
            label="CS Turn Steel conduit",
        )
    )

    # Draw the conductor cross-section as a rectangle
    axis.add_patch(
        patches.Rectangle(
            (t_structural_radial + radius_cs_turn_cable_space, t_structural_vertical),
            dr_cs_turn - ((2 * t_structural_radial) + (2 * radius_cs_turn_cable_space)),
            2 * radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            label="CS Turn Cable Space",
            zorder=2,
        )
    )
    # Plot the right hand circle for the CS turn cable space
    axis.add_patch(
        patches.Circle(
            (
                (dr_cs_turn - t_structural_radial - radius_cs_turn_cable_space),
                dz_cs_turn / 2,
            ),
            radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            zorder=3,
        )
    )
    # Plot the left hand circle for the CS turn cable space
    axis.add_patch(
        patches.Circle(
            ((t_structural_radial + radius_cs_turn_cable_space), dz_cs_turn / 2),
            radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            zorder=3,
        )
    )

    # Add plasma volume, areas and shaping information
    textstr_turn = (
        f"$\\mathbf{{Turn \\ structure:}}$\n\n$A:$ {a_cs_turn:.4e}$ \\ \\text{{m}}^2$\n"
        f"Turn width: {dr_cs_turn:.4e}$ \\ \\text{{m}}$\n"
        f"Turn height: {dz_cs_turn:.4e}$ \\ \\text{{m}}$\n"
        f"Turn width to height ratio: {ld_ratio_cst:.3f}\n"
        f"Steel conduit width: {t_structural_radial:.4e}$ \\ \\text{{m}}$\n"
        f"Radius of turn cable space: {radius_cs_turn_cable_space:.4e}$ \\ \\text{{m}}$\n"
        f"Radius of turn corner: {r_out_cst:.4e}$ \\ \\text{{m}}$\n"
    )

    axis.text(
        0.6,
        0.5,
        textstr_turn,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox={
            "boxstyle": "round",
            "facecolor": "lightyellow",
            "alpha": 1.0,
            "linewidth": 2,
        },
    )

    axis.set_xlim(-dr_cs_turn * 0.2, dr_cs_turn * 1.2)
    axis.set_ylim(-dz_cs_turn * 0.3, dz_cs_turn * 1.3)
    axis.set_aspect("equal")
    axis.set_xlabel("Length [m]")
    axis.set_ylabel("Height [m]")
    axis.set_title("CS Turn Conductor Cross-Section")
    axis.legend()
    axis.grid(True, linestyle="--", alpha=0.3)


def plot_tf_coil_structure(axis, mfile_data, scan, colour_scheme=1):
    # Plot the TF coil poloidal cross-section
    plot_tf_coils(axis, mfile_data, scan, colour_scheme)

    x1 = mfile_data.data["r_tf_arc(1)"].get_scan(scan)
    y1 = mfile_data.data["z_tf_arc(1)"].get_scan(scan)
    x2 = mfile_data.data["r_tf_arc(2)"].get_scan(scan)
    y2 = mfile_data.data["z_tf_arc(2)"].get_scan(scan)
    x3 = mfile_data.data["r_tf_arc(3)"].get_scan(scan)
    y3 = mfile_data.data["z_tf_arc(3)"].get_scan(scan)
    x4 = mfile_data.data["r_tf_arc(4)"].get_scan(scan)
    y4 = mfile_data.data["z_tf_arc(4)"].get_scan(scan)
    x5 = mfile_data.data["r_tf_arc(5)"].get_scan(scan)
    y5 = mfile_data.data["z_tf_arc(5)"].get_scan(scan)

    z_tf_inside_half = mfile_data.data["z_tf_inside_half"].get_scan(scan)
    z_tf_top = mfile_data.data["z_tf_top"].get_scan(scan)
    dr_tf_inboard = mfile_data.data["dr_tf_inboard"].get_scan(scan)
    r_tf_inboard_out = mfile_data.data["r_tf_inboard_out"].get_scan(scan)
    r_tf_outboard_in = mfile_data.data["r_tf_outboard_in"].get_scan(scan)
    r_tf_inboard_in = mfile_data.data["r_tf_inboard_in"].get_scan(scan)
    dr_tf_outboard = mfile_data.data["dr_tf_outboard"].get_scan(scan)
    len_tf_coil = mfile_data.data["len_tf_coil"].get_scan(scan)

    # Plot the points as black dots, number them, and connect them with lines
    xs = [x1, x2, x3, x4, x5]
    ys = [y1, y2, y3, y4, y5]
    labels = []
    for i, (x, y) in enumerate(zip(xs, ys, strict=False), 1):
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

    # =========================================================

    # If D-shaped coil, plot the full internal height arrow
    if mfile_data.data["i_tf_shape"].get_scan(scan) == 1:
        # Arrow for internal coil width
        axis.annotate(
            "",
            xy=(x2, y2),
            xytext=(x4, y4),
            arrowprops={"arrowstyle": "<->", "color": "black"},
        )

        # Add a label for the internal coil width
        axis.text(
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
    if mfile_data.data["i_tf_shape"].get_scan(scan) == 1:
        x = x2 * 0.9
    elif mfile_data.data["i_tf_shape"].get_scan(scan) == 2:
        x = (x2 - x1) / 2

    axis.annotate(
        "",
        xy=(x, y2 + dr_tf_inboard),
        xytext=(x, y4 - dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the full TF coil height
    axis.text(
        x,
        0.0,
        f"{((y2 + 2 * dr_tf_inboard) - y4):.3f} m",
        fontsize=7,
        color="black",
        rotation=270,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
        zorder=100,  # Ensure label is on top of all plots
    )

    # ==========================================================

    # Arrow for top half height of TF coil
    axis.annotate(
        "",
        xy=(-2.0, 0),
        xytext=(-2.0, y2 + dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y2 + dr_tf_inboard, color="black", linestyle="--", linewidth=1)

    # Add a label for top of TF coil
    axis.text(
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
    axis.annotate(
        "",
        xy=(-2.0, 0),
        xytext=(-2.0, y4 - dr_tf_inboard),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y4 - dr_tf_inboard, color="black", linestyle="--", linewidth=1)

    # Add a label for top of TF coil
    axis.text(
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
    axis.annotate(
        "",
        xy=(-1.0, 0),
        xytext=(-1.0, y2),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=y2, color="black", linestyle="--", linewidth=1)

    # Add a label for height of top internal height
    axis.text(
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
    axis.annotate(
        "",
        xy=(-1.0, 0),  # Inner plasma edge
        xytext=(-1.0, -z_tf_inside_half),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    axis.axhline(y=-z_tf_inside_half, color="black", linestyle="--", linewidth=1)

    # Add a label for coil internal height
    axis.text(
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
    axis.annotate(
        "",
        xy=(r_tf_inboard_out, -z_tf_inside_half / 12),
        xytext=(r_tf_outboard_in, -z_tf_inside_half / 12),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the internal coil width
    axis.text(
        (r_tf_inboard_out + r_tf_outboard_in) / 1.5,
        -z_tf_inside_half / 12,
        f"{r_tf_outboard_in - r_tf_inboard_in:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # =============================================================

    # Arrow for full coil width
    axis.annotate(
        "",
        xy=(r_tf_inboard_in, 0.0),
        xytext=(r_tf_outboard_in + dr_tf_outboard, 0.0),
        arrowprops={"arrowstyle": "<|-|>", "color": "black"},
    )

    # Add a label for the full coil width
    axis.text(
        (r_tf_inboard_out + r_tf_outboard_in) / 1.5,
        0.0,
        f"{(r_tf_outboard_in + dr_tf_outboard) - r_tf_inboard_in:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
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
    axis.text(
        r_tf_inboard_in,
        (y4 - dr_tf_inboard) * 1.1,
        rf"$\Delta r = ${dr_tf_inboard:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # Add a label for the outboard thickness
    axis.text(
        r_tf_outboard_in,
        (y4 - dr_tf_inboard) * 1.1,
        rf"$\Delta r = ${dr_tf_outboard:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # ==============================================================

    # Add a label for the inboard thickness
    axis.text(
        (r_tf_outboard_in + 2 * dr_tf_outboard),
        0.0,
        rf"Length of coil = {len_tf_coil:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # ==============================================================

    # Add arow for inboard coil radius
    axis.annotate(
        "",
        xy=(r_tf_inboard_in, 0),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "->", "color": "black"},
    )

    # Add label for inboard coil radius
    axis.text(
        r_tf_inboard_in / 2,
        0.0,
        f"{r_tf_inboard_in:.3f} m",
        fontsize=7,
        color="black",
        verticalalignment="center",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "pink", "alpha": 1.0},
    )

    # =============================================================

    # ==============================================================

    if mfile_data.data["i_tf_shape"].get_scan(scan) == 1:
        # Add arow for inboard coil radius
        axis.annotate(
            "",
            xy=(r_tf_outboard_in + dr_tf_outboard, y2 + dr_tf_inboard),
            xytext=(0, y2 + dr_tf_inboard),
            arrowprops={"arrowstyle": "->", "color": "black"},
        )

        # Add label for inboard coil radius
        axis.text(
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
            0, y2 + dr_tf_inboard, marker="o", color="black", markersize=7, zorder=100
        )

    # ==============================================================

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


def plot_iteration_variables(axis, m_file_data, scan):
    """Plot the iteration variables for a given figure."""
    n_itvars = int(m_file_data.data["nvar"].get_scan(scan))

    y_labels = []
    y_pos = []
    n_plot = 0

    # Build a mapping from itvar index to its name (description)
    itvar_names = {}
    for var in m_file_data.data:
        if var.startswith("itvar"):
            idx = int(var[5:])  # e.g. "itvar001" -> 1
            itvar_names[idx] = m_file_data.data[var].var_description

    for n_plot, n in enumerate(range(1, n_itvars + 1)):
        itvar_final = m_file_data.data[f"itvar{n:03d}"].get_scan(scan)
        itvar_upper = m_file_data.data[f"boundu{n:03d}"].get_scan(scan)
        itvar_lower = m_file_data.data[f"boundl{n:03d}"].get_scan(scan)

        # Use the variable name if available, else fallback to "itvarXXX"
        var_label = itvar_names.get(n, f"itvar{n:03d}")

        # Plot the final value as a vertical marker, normalised to upper bound
        norm_final = (
            (itvar_final - itvar_lower) / (itvar_upper - itvar_lower)
            if itvar_final != itvar_lower
            else 0
        )
        if np.isclose(norm_final, 1.0, atol=1e-6):
            # Plot the lower bound at x=0
            axis.plot(
                1,
                n_plot,
                "s",
                color="black",
                markersize=8,
                label="Lower Bound" if n_plot == 0 else "",
            )
        elif np.isclose(norm_final, 0.0, atol=1e-6):
            # Plot the upper bound at x=1
            axis.plot(
                0,
                n_plot,
                "s",
                color="black",
                markersize=8,
                label="Upper Bound" if n_plot == 0 else "",
            )
        # Draw a horizontal bar from 0 to norm_final at y=n_plot
        else:
            axis.barh(
                n_plot,
                norm_final,
                left=0,
                height=0.8,
                color="blue",
                alpha=0.7,
                label="Final Value" if n_plot == 0 else "",
            )
        # Plot the value as a number at x = 0.5
        axis.text(
            0.5,
            n_plot,
            f"{itvar_final:.8g}",
            va="center",
            ha="center",
            fontsize=8,
            color="black",
        )
        # Plot the value of the upper bound to the right of x=1
        axis.text(
            1.05,
            n_plot,
            f"{itvar_upper:.3g}",
            va="center",
            ha="left",
            fontsize=8,
            color="gray",
        )
        # Plot the value of the lower bound to the left of x=0
        axis.text(
            -0.05,
            n_plot,
            f"{itvar_lower:.3g}",
            va="center",
            ha="right",
            fontsize=8,
            color="gray",
        )
        y_labels.append(var_label)
        y_pos.append(n_plot)

    # Plot vertical lines at x=0 and x=1 to indicate bounds
    axis.axvline(0, color="darkgreen", linewidth=2, zorder=0)
    axis.axvline(1, color="red", linewidth=2, zorder=0)
    axis.set_yticks(y_pos)
    axis.set_yticklabels(y_labels)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_xlim(-0.2, 1.2)  # Normalised bounds
    axis.set_title("Iteration Variables Final Values and Bounds")


def main_plot(
    fig1,
    fig2,
    fig3,
    fig4,
    fig5,
    fig6,
    fig7,
    fig8,
    fig9,
    fig10,
    fig11,
    fig12,
    fig13,
    fig14,
    fig15,
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

    plot_7 = fig2.add_subplot(121)
    plot_7.set_position([0.175, 0.1, 0.35, 0.8])  # Move plot slightly to the right
    plot_iteration_variables(plot_7, m_file_data, scan)

    # Plot main plasma information
    plot_8 = fig3.add_subplot(111, aspect="equal")
    plot_main_plasma_information(plot_8, m_file_data, scan, colour_scheme, fig3)

    # Plot density profiles
    plot_9 = fig4.add_subplot(231)  # , aspect= 0.05)
    plot_9.set_position([0.075, 0.55, 0.25, 0.4])
    plot_n_profiles(plot_9, demo_ranges, m_file_data, scan)

    # Plot temperature profiles
    plot_10 = fig4.add_subplot(232)
    plot_10.set_position([0.375, 0.55, 0.25, 0.4])
    plot_t_profiles(plot_10, demo_ranges, m_file_data, scan)

    # Plot impurity profiles
    plot_11 = fig4.add_subplot(233)
    plot_11.set_position([0.7, 0.45, 0.25, 0.5])
    plot_radprofile(plot_11, m_file_data, scan, imp, demo_ranges)

    # Plot current density profile
    plot_12 = fig4.add_subplot(4, 3, 10)
    plot_12.set_position([0.075, 0.125, 0.25, 0.15])
    plot_jprofile(plot_12)

    # Plot q profile
    plot_13 = fig4.add_subplot(4, 3, 12)
    plot_13.set_position([0.7, 0.125, 0.25, 0.15])
    plot_qprofile(plot_13, demo_ranges, m_file_data, scan)

    # Plot poloidal cross-section
    plot_14 = fig5.add_subplot(121, aspect="equal")
    poloidal_cross_section(plot_14, m_file_data, scan, demo_ranges, colour_scheme)

    # Plot toroidal cross-section
    plot_15 = fig5.add_subplot(122, aspect="equal")
    toroidal_cross_section(plot_15, m_file_data, scan, demo_ranges, colour_scheme)
    # fig4.subplots_adjust(bottom=-0.2, top = 0.9, left = 0.1, right = 0.9)

    # Plot color key
    plot_16 = fig5.add_subplot(222)
    plot_16.set_position([0.5, 0.5, 0.5, 0.5])
    color_key(plot_16, m_file_data, scan, colour_scheme)

    plot_17 = fig6.add_subplot(211)
    plot_17.set_position([0.1, 0.33, 0.8, 0.6])  # x0, y0, width, height (2/3 vertical)
    plot_radial_build(plot_17, m_file_data, colour_scheme)

    # Make each axes smaller vertically to leave room for the legend
    plot_175 = fig7.add_subplot(211)
    plot_175.set_position([0.1, 0.61, 0.8, 0.32])  # x0, y0, width, height

    plot_17 = fig7.add_subplot(212)
    plot_17.set_position([0.1, 0.13, 0.8, 0.32])  # x0, y0, width, height
    plot_upper_vertical_build(plot_175, m_file_data, colour_scheme)
    plot_lower_vertical_build(plot_17, m_file_data, colour_scheme)

    # Can only plot WP and turn structure if superconducting coil at the moment
    if m_file_data.data["i_tf_sup"].get_scan(scan) == 1:
        # TF coil with WP
        plot_19 = fig8.add_subplot(211, aspect="equal")
        plot_19.set_position([0.05, 0.5, 0.7, 0.4])
        plot_tf_wp(plot_19, m_file_data, scan, fig8)

        # TF coil turn structure
        plot_20 = fig8.add_subplot(325, aspect="equal")
        plot_20.set_position([0.025, 0.1, 0.3, 0.3])
        plot_tf_turn(plot_20, m_file_data, scan)

    plot_21 = fig9.add_subplot(111, aspect="equal")
    plot_tf_coil_structure(plot_21, m_file_data, scan, colour_scheme)

    plot_22 = fig10.add_subplot(221)
    plot_bootstrap_comparison(plot_22, m_file_data, scan)

    plot_23 = fig10.add_subplot(224)
    plot_h_threshold_comparison(plot_23, m_file_data, scan)

    plot_24 = fig11.add_subplot(221)
    plot_density_limit_comparison(plot_24, m_file_data, scan)

    plot_25 = fig11.add_subplot(224)
    plot_confinement_time_comparison(plot_25, m_file_data, scan)

    plot_26 = fig12.add_subplot(111)
    plot_current_profiles_over_time(plot_26, m_file_data, scan)

    plot_27 = fig13.add_subplot(121, aspect="equal")
    plot_cs_coil_structure(plot_27, fig13, m_file_data, scan)

    plot_28 = fig13.add_subplot(224, aspect="equal")
    plot_cs_turn_structure(plot_28, fig13, m_file_data, scan)

    plot_29 = fig14.add_subplot(221, aspect="equal")
    plot_first_wall_top_down_cross_section(plot_29, m_file_data, scan)

    plot_30 = fig14.add_subplot(122)
    plot_first_wall_poloidal_cross_section(plot_30, m_file_data, scan)

    plot_31 = fig15.add_subplot(111, aspect="equal")
    plot_main_power_flow(plot_31, m_file_data, scan, fig15)


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
    global dx_tf_wp_primary_toroidal
    global dx_tf_wp_secondary_toroidal
    global dr_tf_wp_with_insulation
    global dx_tf_wp_insulation
    global dr_tf_nose_case
    global dr_tf_plasma_case

    n_tf_coils = m_file.data["n_tf_coils"].get_scan(scan)
    if i_tf_sup == 1:  # If superconducting magnets
        dx_tf_wp_primary_toroidal = m_file.data["dx_tf_wp_primary_toroidal"].get_scan(
            scan
        )
        if i_tf_wp_geom == 1:
            dx_tf_wp_secondary_toroidal = m_file.data[
                "dx_tf_wp_secondary_toroidal"
            ].get_scan(scan)
        dr_tf_wp_with_insulation = m_file.data["dr_tf_wp_with_insulation"].get_scan(
            scan
        )
        dx_tf_wp_insulation = m_file.data["dx_tf_wp_insulation"].get_scan(scan)
        dr_tf_nose_case = m_file.data["dr_tf_nose_case"].get_scan(scan)

        # To be re-inergrated to resistives when in-plane stresses is integrated
        dr_tf_plasma_case = m_file.data["dr_tf_plasma_case"].get_scan(scan)

    global dx_beam_shield
    global radius_beam_tangency
    global radius_beam_tangency_max
    global dx_beam_duct

    i_hcd_primary = int(m_file.data["i_hcd_primary"].get_scan(scan))
    i_hcd_secondary = int(m_file.data["i_hcd_secondary"].get_scan(scan))

    if (i_hcd_primary in [5, 8]) or (i_hcd_secondary in [5, 8]):
        dx_beam_shield = m_file.data["dx_beam_shield"].get_scan(scan)
        radius_beam_tangency = m_file.data["radius_beam_tangency"].get_scan(scan)
        radius_beam_tangency_max = m_file.data["radius_beam_tangency_max"].get_scan(
            scan
        )
        dx_beam_duct = m_file.data["dx_beam_duct"].get_scan(scan)
    else:
        dx_beam_shield = radius_beam_tangency = radius_beam_tangency_max = (
            dx_beam_duct
        ) = 0.0

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
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_divertor",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
            "dr_tf_shld_gap",
            "dr_tf_inboard",
        ]
    else:
        vertical_upper = [
            "z_plasma_xpoint_upper",
            "dz_fw_plasma_gap",
            "dz_fw_upper",
            "dz_blkt_upper",
            "dr_shld_blkt_gap",
            "dz_shld_upper",
            "dz_vv_upper",
            "dz_shld_vv_gap",
            "dz_shld_thermal",
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
    page8 = plt.figure(figsize=(12, 9), dpi=80)
    page9 = plt.figure(figsize=(12, 9), dpi=80)
    page10 = plt.figure(figsize=(12, 9), dpi=80)
    page11 = plt.figure(figsize=(12, 9), dpi=80)
    page12 = plt.figure(figsize=(12, 9), dpi=80)
    page13 = plt.figure(figsize=(12, 9), dpi=80)
    page14 = plt.figure(figsize=(12, 9), dpi=80)
    page15 = plt.figure(figsize=(12, 9), dpi=80)

    # run main_plot
    main_plot(
        page1,
        page2,
        page3,
        page4,
        page5,
        page6,
        page7,
        page8,
        page9,
        page10,
        page11,
        page12,
        page13,
        page14,
        page15,
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
        pdf.savefig(page8)
        pdf.savefig(page9)
        pdf.savefig(page10)
        pdf.savefig(page11)
        pdf.savefig(page12)
        pdf.savefig(page13)
        pdf.savefig(page14)
        pdf.savefig(page15)

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
    plt.close(page8)
    plt.close(page9)
    plt.close(page10)
    plt.close(page11)
    plt.close(page12)
    plt.close(page13)
    plt.close(page14)
    plt.close(page15)


if __name__ == "__main__":
    main()
