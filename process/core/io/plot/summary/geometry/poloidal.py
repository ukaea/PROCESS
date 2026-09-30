"""Geometry functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

from process.core.io.plot.summary.constants import (
    BLANKET_COLOUR,
    CRYOSTAT_COLOUR,
    FIRSTWALL_COLOUR,
    SHIELD_COLOUR,
    VESSEL_COLOUR,
    thin,
)
from process.core.io.plot.summary.magnets.pf import (
    plot_pf_coils,
)
from process.core.io.plot.summary.magnets.tf import (
    plot_tf_coils,
)
from process.core.io.plot.summary.plasma.physics import (
    plot_plasma,
)
from process.core.io.plot.summary.radial_build import (
    cumulative_radial_build,
)
from process.core.io.plot.summary.rendering import (
    draw_annotation,
    draw_text,
)
from process.core.io.plot.summary.reporting.misc import (
    plot_centre_cross,
)
from process.data_structure.physics_variables import DivertorNumberModels
from process.models.geometry.blanket import (
    blanket_geometry_double_null,
    blanket_geometry_single_null,
)
from process.models.geometry.cryostat import cryostat_geometry
from process.models.geometry.firstwall import (
    first_wall_geometry_double_null,
    first_wall_geometry_single_null,
)
from process.models.geometry.shield import (
    shield_geometry_double_null,
    shield_geometry_single_null,
)
from process.models.geometry.vacuum_vessel import (
    vacuum_vessel_geometry_double_null,
    vacuum_vessel_geometry_single_null,
)

if TYPE_CHECKING:
    from process.core.io.mfile import MFile
    from process.core.io.plot.summary.reporting.misc import (
        RadialBuild,
    )


def poloidal_cross_section(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    demo_ranges: bool,
    radial_build: RadialBuild,
    colour_scheme: Literal[1, 2],
):
    """Function to plot poloidal cross-section

    Parameters
    ----------
    axis :
        axis object to add plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    demo_ranges:

    colour_scheme :
        colour scheme to use for plots
    """
    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("Poloidal Cross-Section")
    axis.minorticks_on()
    axis.grid(which="both", linestyle="--", linewidth=0.5, alpha=0.2)

    plot_vacuum_vessel_and_divertor(axis, mfile, scan, radial_build, colour_scheme)
    plot_shield(axis, mfile, scan, radial_build, colour_scheme)
    plot_blanket(axis, mfile, scan, radial_build, colour_scheme)
    plot_firstwall(axis, mfile, scan, radial_build, colour_scheme)

    plot_plasma(axis, mfile, scan, colour_scheme)
    plot_centre_cross(axis, mfile, scan)
    plot_cryostat(axis, mfile, scan, colour_scheme)

    plot_tf_coils(axis, mfile, scan, colour_scheme)
    plot_pf_coils(axis, mfile, scan, colour_scheme)

    if demo_ranges:
        axis.set_ylim(-15, 15)
        axis.set_xlim(0, 20)

    else:
        axis.set_xlim(0, axis.get_xlim()[1])


def plot_full_machine_poloidal_cross_section(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    radial_build: RadialBuild,
    colour_scheme: Literal[1, 2],
):
    """Function to plot full machine poloidal cross-section, including mirrored negative
    x-axis

    Parameters
    ----------
    axis :
        axis object to add plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    radial_build :
        radial build data
    colour_scheme :
        colour scheme to use for plots
    """
    plot_vacuum_vessel_and_divertor(axis, mfile, scan, radial_build, colour_scheme)
    plot_vacuum_vessel_and_divertor(
        axis, mfile, scan, radial_build, colour_scheme, mirror_negative_x=True
    )
    plot_shield(axis, mfile, scan, radial_build, colour_scheme)
    plot_shield(axis, mfile, scan, radial_build, colour_scheme, mirror_negative_x=True)

    plot_blanket(axis, mfile, scan, radial_build, colour_scheme)
    plot_blanket(axis, mfile, scan, radial_build, colour_scheme, mirror_negative_x=True)
    plot_firstwall(axis, mfile, scan, radial_build, colour_scheme)
    plot_firstwall(
        axis, mfile, scan, radial_build, colour_scheme, mirror_negative_x=True
    )
    plot_plasma(axis, mfile, scan, colour_scheme)
    plot_plasma(axis, mfile, scan, colour_scheme, mirror_negative_x=True)
    plot_centre_cross(axis, mfile, scan)
    plot_centre_cross(axis, mfile, scan, mirror_negative_x=True)
    plot_cryostat(axis, mfile, scan, colour_scheme)
    plot_cryostat(axis, mfile, scan, colour_scheme, mirror_negative_x=True)
    plot_tf_coils(axis, mfile, scan, colour_scheme)
    plot_tf_coils(axis, mfile, scan, colour_scheme, mirror_negative_x=True)
    plot_pf_coils(axis, mfile, scan, colour_scheme)
    plot_pf_coils(axis, mfile, scan, colour_scheme, mirror_negative_x=True)

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_aspect("equal")
    axis.minorticks_on()
    axis.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)


def plot_cryostat(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    mirror_negative_x: bool = False,
):
    """Function to plot cryostat in poloidal cross-section

    Parameters
    ----------
    axis : plt.Axes
        axis object to plot to
    mfile : MFile
        MFILE data object
    scan : int
        scan number to use
    colour_scheme : Literal[1, 2]
        colour scheme to use for plots
    mirror_negative_x : bool
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    rects = cryostat_geometry(
        r_cryostat_inboard=mfile.get("r_cryostat_inboard", scan=scan),
        dr_cryostat=mfile.get("dr_cryostat", scan=scan),
        z_cryostat_half_inside=mfile.get("z_cryostat_half_inside", scan=scan),
    )

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    for rec in rects:
        axis.add_patch(
            patches.Rectangle(
                xy=(x_scale * rec.anchor_x, rec.anchor_z),
                width=x_scale * rec.width,
                height=rec.height,
                facecolor=CRYOSTAT_COLOUR[colour_scheme - 1],
            )
        )


def plot_vacuum_vessel_and_divertor(
    axis,
    mfile: MFile,
    scan,
    radial_build,
    colour_scheme,
    mirror_negative_x: bool = False,
):
    """Function to plot vacuum vessel and divertor boxes

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    radial_build :

    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    cumulative_upper = radial_build.cumulative_upper
    cumulative_lower = radial_build.cumulative_lower
    upper = radial_build.upper
    lower = radial_build.lower

    i_single_null = int(mfile.get("i_single_null", scan=scan))
    triang_95 = mfile.get("triang95", scan=scan)
    dz_divertor = mfile.get("dz_divertor", scan=scan)
    dz_xpoint_divertor = mfile.get("dz_xpoint_divertor", scan=scan)
    kappa = mfile.get("kappa", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    dr_vv_inboard = mfile.get("dr_vv_inboard", scan=scan)
    dr_vv_outboard = mfile.get("dr_vv_outboard", scan=scan)
    dr_shld_inboard = mfile.get("dr_shld_inboard", scan=scan)
    dr_shld_outboard = mfile.get("dr_shld_outboard", scan=scan)
    dr_blkt_inboard = mfile.get("dr_blkt_inboard", scan=scan)
    dr_blkt_outboard = mfile.get("dr_blkt_outboard", scan=scan)

    # Outer side (furthest from plasma)
    radx_outer = (
        cumulative_radial_build("dr_vv_outboard", mfile, scan)
        + cumulative_radial_build("dr_shld_vv_gap_inboard", mfile, scan)
    ) / 2.0
    rminx_outer = (
        cumulative_radial_build("dr_vv_outboard", mfile, scan)
        - cumulative_radial_build("dr_shld_vv_gap_inboard", mfile, scan)
    ) / 2.0

    # Inner side (nearest to the plasma)
    radx_inner = (
        cumulative_radial_build("dr_shld_outboard", mfile, scan)
        + cumulative_radial_build("dr_vv_inboard", mfile, scan)
    ) / 2.0
    rminx_inner = (
        cumulative_radial_build("dr_shld_outboard", mfile, scan)
        - cumulative_radial_build("dr_vv_inboard", mfile, scan)
    ) / 2.0

    z_divertor_lower_top = (-kappa * rminor) - dz_xpoint_divertor
    z_divertor_lower_bottom = z_divertor_lower_top - dz_divertor

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    match DivertorNumberModels(i_single_null):
        case DivertorNumberModels.SINGLE_NULL:
            z_divertor_upper_bottom = None
            z_divertor_upper_top = None
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
                x_scale * np.array(vvg_single_null.rs),
                vvg_single_null.zs,
                color="black",
                lw=thin,
                zorder=5,
            )

            axis.fill(
                x_scale * np.array(vvg_single_null.rs),
                vvg_single_null.zs,
                color=VESSEL_COLOUR[colour_scheme - 1],
                lw=0.01,
                zorder=5,
            )

            # Find indices where vessel boundary is between z_divertor_bottom and
            # z_divertor_top
            # Find the min and max R values of the vessel boundary between the divertor
            # lines
            mask = (vvg_single_null.zs >= z_divertor_lower_bottom) & (
                vvg_single_null.zs <= z_divertor_lower_top
            )
            # Get the min/max R for the region between the divertor lines
            r_min = (
                np.min(vvg_single_null.rs[mask])
                + dr_vv_inboard
                + dr_shld_inboard
                + (dr_blkt_inboard * 0.5)
            )
            r_max = (
                np.max(vvg_single_null.rs[mask])
                - dr_vv_outboard
                - dr_shld_outboard
                - (dr_blkt_outboard * 0.5)
            )
            # Draw a rectangle (box) between the two lines and inside the vessel
            axis.add_patch(
                patches.Rectangle(
                    (
                        x_scale * r_min,
                        z_divertor_lower_bottom,
                    ),
                    x_scale * (r_max - r_min),
                    z_divertor_lower_top - z_divertor_lower_bottom,
                    facecolor="black",
                    alpha=0.8,
                    zorder=1,
                )
            )

        case DivertorNumberModels.DOUBLE_NULL:
            z_divertor_upper_bottom = (kappa * rminor) + dz_xpoint_divertor
            z_divertor_upper_top = z_divertor_upper_bottom + dz_divertor
            vvg_double_null = vacuum_vessel_geometry_double_null(
                cumulative_lower=cumulative_lower,
                lower=lower,
                radx_inner=radx_inner,
                radx_outer=radx_outer,
                rminx_inner=rminx_inner,
                rminx_outer=rminx_outer,
                triang=triang_95,
            )
            axis.plot(
                x_scale * np.array(vvg_double_null.rs),
                vvg_double_null.zs,
                color="black",
                lw=thin,
                zorder=5,
            )

            axis.fill(
                x_scale * np.array(vvg_double_null.rs),
                vvg_double_null.zs,
                color=VESSEL_COLOUR[colour_scheme - 1],
                lw=0.01,
                zorder=5,
            )

            # Plot lower divertor
            # Find indices where vessel boundary is between z_divertor_bottom and
            # z_divertor_top
            # Find the min and max R values of the vessel boundary between the divertor
            # lines
            mask = (vvg_double_null.zs >= z_divertor_lower_bottom) & (
                vvg_double_null.zs <= z_divertor_lower_top
            )
            # Get the min/max R for the region between the divertor lines
            r_min = (
                np.min(vvg_double_null.rs[mask])
                + dr_vv_inboard
                + dr_shld_inboard
                + (dr_blkt_inboard * 0.5)
            )
            r_max = (
                np.max(vvg_double_null.rs[mask])
                - dr_vv_outboard
                - dr_shld_outboard
                - (dr_blkt_outboard * 0.5)
            )
            # Draw a rectangle (box) between the two lines and inside the vessel
            axis.add_patch(
                patches.Rectangle(
                    (
                        x_scale * r_min,
                        z_divertor_lower_bottom,
                    ),
                    x_scale * (r_max - r_min),
                    z_divertor_lower_top - z_divertor_lower_bottom,
                    facecolor="black",
                    alpha=0.8,
                    zorder=1,
                )
            )
            # Plot upper divertor
            # Find indices where vessel boundary is between z_divertor_bottom and
            # z_divertor_top
            # Find the min and max R values of the vessel boundary between the divertor
            # lines
            mask = (vvg_double_null.zs >= z_divertor_upper_bottom) & (
                vvg_double_null.zs <= z_divertor_upper_top
            )
            # Get the min/max R for the region between the divertor lines
            r_min = (
                np.min(vvg_double_null.rs[mask])
                + dr_vv_inboard
                + dr_shld_inboard
                + (dr_blkt_inboard * 0.5)
            )
            r_max = (
                np.max(vvg_double_null.rs[mask])
                - dr_vv_outboard
                - dr_shld_outboard
                - (dr_blkt_outboard * 0.5)
            )
            # Draw a rectangle (box) between the two lines and inside the vessel
            axis.add_patch(
                patches.Rectangle(
                    (
                        x_scale * r_min,
                        z_divertor_upper_bottom,
                    ),
                    x_scale * (r_max - r_min),
                    z_divertor_upper_top - z_divertor_upper_bottom,
                    facecolor="black",
                    alpha=0.8,
                    zorder=1,
                )
            )


def plot_shield(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    radial_build,
    colour_scheme,
    mirror_negative_x: bool = False,
):
    """Function to plot shield

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    radial_build :

    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)

    """
    cumulative_upper = radial_build.cumulative_upper
    cumulative_lower = radial_build.cumulative_lower

    i_single_null = mfile.get("i_single_null", scan=scan)
    triang_95 = mfile.get("triang95", scan=scan)

    # Side furthest from plasma
    radx_far = (
        cumulative_radial_build("dr_shld_outboard", mfile, scan)
        + cumulative_radial_build("dr_vv_inboard", mfile, scan)
    ) / 2.0
    rminx_far = (
        cumulative_radial_build("dr_shld_outboard", mfile, scan)
        - cumulative_radial_build("dr_vv_inboard", mfile, scan)
    ) / 2.0

    # Side nearest to the plasma
    radx_near = (
        cumulative_radial_build("vvblgapo", mfile, scan)
        + cumulative_radial_build("dr_shld_inboard", mfile, scan)
    ) / 2.0
    rminx_near = (
        cumulative_radial_build("vvblgapo", mfile, scan)
        - cumulative_radial_build("dr_shld_inboard", mfile, scan)
    ) / 2.0

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    match DivertorNumberModels(i_single_null):
        case DivertorNumberModels.SINGLE_NULL:
            shield_geometry = shield_geometry_single_null(
                cumulative_upper=cumulative_upper,
                radx_far=radx_far,
                rminx_far=rminx_far,
                radx_near=radx_near,
                rminx_near=rminx_near,
                triang=triang_95,
                cumulative_lower=cumulative_lower,
            )
        case DivertorNumberModels.DOUBLE_NULL:
            shield_geometry = shield_geometry_double_null(
                cumulative_lower=cumulative_lower,
                radx_far=radx_far,
                radx_near=radx_near,
                rminx_far=rminx_far,
                rminx_near=rminx_near,
                triang=triang_95,
            )

    axis.plot(
        x_scale * np.array(shield_geometry.rs),
        shield_geometry.zs,
        color="black",
        lw=thin,
    )
    axis.fill(
        x_scale * np.array(shield_geometry.rs),
        shield_geometry.zs,
        color=SHIELD_COLOUR[colour_scheme - 1],
        lw=0.01,
    )


def plot_blanket(
    axis: plt.Axes,
    mfile: MFile,
    scan,
    radial_build,
    colour_scheme,
    mirror_negative_x: bool = False,
):
    """Function to plot blanket

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    radial_build :

    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    cumulative_upper = radial_build.cumulative_upper
    cumulative_lower = radial_build.cumulative_lower

    dr_blkt_inboard = mfile.get("dr_blkt_inboard", scan=scan)
    dr_blkt_outboard = mfile.get("dr_blkt_outboard", scan=scan)
    # Single null: Draw top half from output
    # Double null: Reflect bottom half to top
    i_single_null = mfile.get("i_single_null", scan=scan)
    triang_95 = mfile.get("triang95", scan=scan)
    if int(i_single_null) == 1:
        dz_blkt_upper = mfile.get("dz_blkt_upper", scan=scan)
    else:
        dz_blkt_upper = 0.0

    c_shldith = cumulative_radial_build("dr_shld_inboard", mfile, scan)
    c_blnkoth = cumulative_radial_build("dr_blkt_outboard", mfile, scan)

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    match DivertorNumberModels(i_single_null):
        case DivertorNumberModels.SINGLE_NULL:
            # Upper blanket: outer surface
            radx_outer = (
                cumulative_radial_build("dr_blkt_outboard", mfile, scan)
                + cumulative_radial_build("vvblgapi", mfile, scan)
            ) / 2.0
            rminx_outer = (
                cumulative_radial_build("dr_blkt_outboard", mfile, scan)
                - cumulative_radial_build("vvblgapi", mfile, scan)
            ) / 2.0

            # Upper blanket: inner surface
            radx_inner = (
                cumulative_radial_build("dr_fw_outboard", mfile, scan)
                + cumulative_radial_build("dr_blkt_inboard", mfile, scan)
            ) / 2.0
            rminx_inner = (
                cumulative_radial_build("dr_fw_outboard", mfile, scan)
                - cumulative_radial_build("dr_blkt_inboard", mfile, scan)
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
                x_scale * np.array(bg_single_null.rs),
                bg_single_null.zs,
                color="black",
                lw=thin,
                zorder=5,
            )

            axis.fill(
                x_scale * np.array(bg_single_null.rs),
                bg_single_null.zs,
                color=BLANKET_COLOUR[colour_scheme - 1],
                lw=0.01,
                zorder=5,
            )

        case DivertorNumberModels.DOUBLE_NULL:
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
            axis.plot(
                x_scale * np.array(bg_double_null.rs[0]),
                bg_double_null.zs[0],
                color="black",
                lw=thin,
            )
            axis.fill(
                x_scale * np.array(bg_double_null.rs[0]),
                bg_double_null.zs[0],
                color=BLANKET_COLOUR[colour_scheme - 1],
                lw=0.01,
                zorder=5,
            )
            if dr_blkt_inboard > 0.0:
                # only plot inboard blanket if inboard blanket thickness > 0
                axis.plot(
                    x_scale * np.array(bg_double_null.rs[1]),
                    bg_double_null.zs[1],
                    color="black",
                    lw=thin,
                    zorder=5,
                )
                axis.fill(
                    x_scale * np.array(bg_double_null.rs[1]),
                    bg_double_null.zs[1],
                    color=BLANKET_COLOUR[colour_scheme - 1],
                    lw=0.01,
                    zorder=5,
                )


def plot_first_wall_top_down_cross_section(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot first wall top down cross-section"""
    # Import required variables
    radius_fw_channel = mfile.get("radius_fw_channel", scan=scan) * 100
    dr_fw_wall = mfile.get("dr_fw_wall", scan=scan) * 100
    dx_fw_module = mfile.get("dx_fw_module", scan=scan) * 100

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
            xy=(
                dx_fw_module + dx_fw_module / 2,
                dr_fw_wall + radius_fw_channel,
            ),
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
    draw_annotation(
        axis,
        "",
        xy=(dx_fw_module, -0.2),
        xytext=(2 * dx_fw_module, -0.2),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Draw dotted line above the channel
    axis.plot(
        [dx_fw_module * 1.5, dx_fw_module * 1.5],
        [
            2 * radius_fw_channel + dr_fw_wall,
            2 * (radius_fw_channel + dr_fw_wall),
        ],
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
    axis.set_xlim(-1, 2 * dx_fw_module + 1)
    axis.set_ylim(-1, 2 * (dr_fw_wall + radius_fw_channel) + 1)


def plot_first_wall_poloidal_cross_section(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot first wall poloidal cross-section"""
    # Import required variables
    radius_fw_channel = mfile.get("radius_fw_channel", scan=scan)
    dr_fw_wall = mfile.get("dr_fw_wall", scan=scan)
    dx_fw_module = mfile.get("dx_fw_module", scan=scan)
    len_fw_channel = mfile.get("len_fw_channel", scan=scan)
    temp_fw_coolant_in = mfile.get("temp_fw_coolant_in", scan=scan)
    temp_fw_coolant_out = mfile.get("temp_fw_coolant_out", scan=scan)
    i_fw_coolant_type = mfile.get("i_fw_coolant_type", scan=scan).strip("'\"")
    temp_fw_peak = mfile.get("temp_fw_peak", scan=scan)
    pres_fw_coolant = mfile.get("pres_fw_coolant", scan=scan)
    n_fw_outboard_channels = mfile.get("n_fw_outboard_channels", scan=scan)
    n_fw_inboard_channels = mfile.get("n_fw_inboard_channels", scan=scan)

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
    draw_text(
        axis,
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
    draw_text(
        axis,
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
        rf"$T_{{FW,peak}}$: {temp_fw_peak:,.3f} K",
        rf"$P_{{FW}}$: {pres_fw_coolant / 1e3:,.3f} kPa",
        rf"$P_{{FW}}$: {pres_fw_coolant / 1e5:,.3f} bar",
        rf"$N_{{outboard}}$: {n_fw_outboard_channels}",
        rf"$N_{{inboard}}$: {n_fw_inboard_channels}",
    ))

    props_fw = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    draw_text(
        axis,
        -0.5,
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
    axis.set_xlim(-0.01, (dx_fw_module + radius_fw_channel * 2) + 0.01)
    axis.set_ylim(-0.2, len_fw_channel + 0.2)


def plot_firstwall(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    radial_build,
    colour_scheme,
    mirror_negative_x: bool = False,
):
    """Function to plot first wall

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    radial_build :

    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    cumulative_upper = radial_build.cumulative_upper
    cumulative_lower = radial_build.cumulative_lower

    i_single_null = mfile.get("i_single_null", scan=scan)
    triang_95 = mfile.get("triang95", scan=scan)
    if int(i_single_null) == 1:
        dz_blkt_upper = mfile.get("dz_blkt_upper", scan=scan)
        tfwvt = mfile.get("dz_fw_upper", scan=scan)
    else:
        dz_blkt_upper = tfwvt = 0.0

    c_blnkith = cumulative_radial_build("dr_blkt_inboard", mfile, scan)
    c_fwoth = cumulative_radial_build("dr_fw_outboard", mfile, scan)

    dr_fw_inboard = mfile.get("dr_fw_inboard", scan=scan)
    dr_fw_outboard = mfile.get("dr_fw_outboard", scan=scan)

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    match DivertorNumberModels(i_single_null):
        case DivertorNumberModels.SINGLE_NULL:
            # Upper first wall: outer surface
            radx_outer = (
                cumulative_radial_build("dr_fw_outboard", mfile, scan)
                + cumulative_radial_build("dr_blkt_inboard", mfile, scan)
            ) / 2.0
            rminx_outer = (
                cumulative_radial_build("dr_fw_outboard", mfile, scan)
                - cumulative_radial_build("dr_blkt_inboard", mfile, scan)
            ) / 2.0

            # Upper first wall: inner surface
            radx_inner = (
                cumulative_radial_build("dr_fw_plasma_gap_outboard", mfile, scan)
                + cumulative_radial_build("dr_fw_inboard", mfile, scan)
            ) / 2.0
            rminx_inner = (
                cumulative_radial_build("dr_fw_plasma_gap_outboard", mfile, scan)
                - cumulative_radial_build("dr_fw_inboard", mfile, scan)
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
            axis.plot(
                x_scale * np.array(fwg_single_null.rs),
                fwg_single_null.zs,
                color="black",
                lw=thin,
            )
            axis.fill(
                x_scale * np.array(fwg_single_null.rs),
                fwg_single_null.zs,
                color=FIRSTWALL_COLOUR[colour_scheme - 1],
                lw=0.01,
            )

        case DivertorNumberModels.DOUBLE_NULL:
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
            axis.plot(
                x_scale * np.array(fwg_double_null.rs[0]),
                fwg_double_null.zs[0],
                color="black",
                lw=thin,
            )
            axis.plot(
                x_scale * np.array(fwg_double_null.rs[1]),
                fwg_double_null.zs[1],
                color="black",
                lw=thin,
            )
            axis.fill(
                x_scale * np.array(fwg_double_null.rs[0]),
                fwg_double_null.zs[0],
                color=FIRSTWALL_COLOUR[colour_scheme - 1],
                lw=0.01,
            )
            axis.fill(
                x_scale * np.array(fwg_double_null.rs[1]),
                fwg_double_null.zs[1],
                color=FIRSTWALL_COLOUR[colour_scheme - 1],
                lw=0.01,
            )


__all__ = [
    "plot_blanket",
    "plot_cryostat",
    "plot_first_wall_poloidal_cross_section",
    "plot_first_wall_top_down_cross_section",
    "plot_firstwall",
    "plot_full_machine_poloidal_cross_section",
    "plot_shield",
    "plot_vacuum_vessel_and_divertor",
    "poloidal_cross_section",
]
