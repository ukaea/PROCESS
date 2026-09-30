"""Plasma functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary.common import (
    box_style,
)
from process.core.io.plot.summary.constants import (
    PLASMA_COLOUR,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)
from process.models.build import Build
from process.models.geometry.plasma import plasma_geometry
from process.models.physics.physics import (
    BetaNormMaxModel,
)
from process.models.physics.plasma_current import (
    PlasmaCurrentModel,
)
from process.models.physics.plasma_geometry import (
    PlasmaShapeModelType,
)

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_plasma(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    mirror_negative_x: bool = False,
):
    """Plots the plasma boundary arcs.

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    colour_scheme :
        colour scheme to use for plots
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)


    Raises
    ------
    ValueError
        If an unsupported plasma shape model type is encountered.
    """
    r_0, a, triang, kappa, i_single_null, i_plasma_shape, plasma_square = (
        mfile.get_variables(
            "rmajor",
            "rminor",
            "triang",
            "kappa",
            "i_single_null",
            "i_plasma_shape",
            "plasma_square",
            scan=scan,
        )
    )

    pg = plasma_geometry(
        rmajor=r_0,
        rminor=a,
        triang=triang,
        kappa=kappa,
        i_single_null=i_single_null,
        i_plasma_shape=i_plasma_shape,
        square=plasma_square,
    )

    # Apply mirror transformation if requested
    x_scale = -1 if mirror_negative_x else 1

    match PlasmaShapeModelType(i_plasma_shape):
        case PlasmaShapeModelType.PROCESS_ORIGINAL:
            # Plot the 2 plasma outline arcs.
            axis.plot(x_scale * np.array(pg.rs[0]), pg.zs[0], color="black")
            axis.plot(x_scale * np.array(pg.rs[1]), pg.zs[1], color="black")

            # Set triang_95 to stop plotting plasma past boundary
            # Assume IPDG scaling
            triang_95 = triang / 1.5

            # Colour in right side of plasma
            axis.fill_between(
                x=x_scale * np.array(pg.rs[0]),
                y1=pg.zs[0],
                where=(pg.rs[0] > r_0 - (triang_95 * a * 1.5)),
                color=PLASMA_COLOUR[colour_scheme - 1],
            )
            # Colour in left side of plasma
            axis.fill_between(
                x=x_scale * np.array(pg.rs[1]),
                y1=pg.zs[1],
                where=(pg.rs[1] < r_0 - (triang_95 * a * 1.5)),
                color=PLASMA_COLOUR[colour_scheme - 1],
            )

        case PlasmaShapeModelType.SAUTER:
            axis.plot(x_scale * np.array(pg.rs), pg.zs, color="black")
            axis.fill(
                x_scale * np.array(pg.rs),
                pg.zs,
                color=PLASMA_COLOUR[colour_scheme - 1],
            )
        case _:
            raise ValueError(f"Unsupported plasma shape model type: {i_plasma_shape}")


def plot_plasma_current_comparison(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot a scatter box plot of different plasma current comparisons.

    Parameters
    ----------
    axis :
        Axis object to plot to.
    mfile :
        MFILE data object.
    scan :
        Scan number to use.
    """
    c_plasma_peng_analytic = mfile.get("c_plasma_peng_analytic", scan=scan)
    c_plasma_peng_double_null = mfile.get("c_plasma_peng_double_null", scan=scan)
    c_plasma_cyclindrical = mfile.get("c_plasma_cyclindrical", scan=scan)
    c_plasma_ipdg89 = mfile.get("c_plasma_ipdg89", scan=scan)
    c_plasma_todd_empirical_i = mfile.get("c_plasma_todd_empirical_i", scan=scan)
    c_plasma_todd_empirical_ii = mfile.get("c_plasma_todd_empirical_ii", scan=scan)
    c_plasma_connor_hastie = mfile.get("c_plasma_connor_hastie", scan=scan)
    c_plasma_sauter = mfile.get("c_plasma_sauter", scan=scan)
    c_plasma_fiesta_st = mfile.get("c_plasma_fiesta_st", scan=scan)

    # Data for the box plot
    data = {
        f"{PlasmaCurrentModel.PENG_ANALYTIC_FIT.full_name}": (c_plasma_peng_analytic),
        f"{PlasmaCurrentModel.PENG_DIVERTOR_SCALING.full_name}": (
            c_plasma_peng_double_null
        ),
        f"{PlasmaCurrentModel.ITER_SCALING.full_name}": c_plasma_cyclindrical,
        f"{PlasmaCurrentModel.IPDG89_SCALING.full_name}": c_plasma_ipdg89,
        f"{PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_I.full_name}": (
            c_plasma_todd_empirical_i
        ),
        f"{PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_II.full_name}": (
            c_plasma_todd_empirical_ii
        ),
        f"{PlasmaCurrentModel.CONNOR_HASTIE_MODEL.full_name}": (c_plasma_connor_hastie),
        f"{PlasmaCurrentModel.SAUTER_SCALING.full_name}": c_plasma_sauter,
        f"{PlasmaCurrentModel.FIESTA_ST_SCALING.full_name}": (c_plasma_fiesta_st),
    }

    # Create the violin plot
    data_values = list(data.values())
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(-0.9, 1))

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_density_limit = np.mean(data_values)
    std_density_limit = np.std(data_values)
    median_density_limit = np.median(data_values)

    # Plot average, standard deviation, and median as text
    draw_text(
        axis,
        -0.45,
        0.15,
        rf"Average: {avg_density_limit * 1e-6:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        -0.45,
        0.1,
        rf"Standard Dev: {std_density_limit * 1e-6:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        -0.45,
        0.05,
        rf"Median: {median_density_limit * 1e-6:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Plasma Current ($I_p$) Comparison")
    axis.set_ylabel(r"Plasma Current [MA]")
    axis.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x * 1e-6:.1f}"))
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def plot_max_normalised_beta_comparison(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot a scatter box plot of different max normalised beta comparisons.

    Parameters
    ----------
    axis :
        Axis object to plot to.
    mfile :
        MFILE data object.
    scan :
        Scan number to use.
    """
    beta_norm_max_wesson = mfile.get("beta_norm_max_wesson", scan=scan)
    beta_norm_max_original_scaling = mfile.get(
        "beta_norm_max_original_scaling", scan=scan
    )
    beta_norm_max_menard = mfile.get("beta_norm_max_menard", scan=scan)
    beta_norm_max_tholerus = mfile.get("beta_norm_max_tholerus", scan=scan)
    beta_norm_max_stambaugh = mfile.get("beta_norm_max_stambaugh", scan=scan)

    # Data for the box plot
    data = {
        f"{BetaNormMaxModel.WESSON.full_name}": beta_norm_max_wesson,
        f"{BetaNormMaxModel.ORIGINAL_SCALING.full_name}": (
            beta_norm_max_original_scaling
        ),
        f"{BetaNormMaxModel.MENARD.full_name}": beta_norm_max_menard,
        f"{BetaNormMaxModel.THOLERUS.full_name}": beta_norm_max_tholerus,
        f"{BetaNormMaxModel.STAMBAUGH.full_name}": beta_norm_max_stambaugh,
    }
    data_values = list(data.values())
    # Create the violin plot
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1.1, 1))

    # Calculate average, standard deviation, and median
    data_values = list(data.values())
    avg_beta_norm_max = np.mean(data_values)
    std_beta_norm_max = np.std(data_values)
    median_beta_norm_max = np.median(data_values)

    # Plot average, standard deviation, and median as text
    draw_text(
        axis,
        1.1,
        0.15,
        rf"Average: {avg_beta_norm_max:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.1,
        0.1,
        rf"Standard Dev: {std_beta_norm_max:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.1,
        0.05,
        rf"Median: {median_beta_norm_max:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Max Normalised Beta ($\\beta_N$) Comparison")
    axis.set_ylabel("Max Normalised Beta $\\beta_N$ [unitless]")
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def reaction_plot_grid(
    rminor,
    rmajor,
    kappa,
    r_grid,
    z_grid,
    grid,
    ax,
    fractions=(0.25, 0.5, 0.75),
    colours=("blue", "yellow", "red"),
):
    """Plot fusion reaction rate density"""
    # Mask points outside the plasma boundary (optional, but grid is inside by
    # construction)
    # Plot filled contour

    upper = ax.contourf(r_grid, z_grid, grid, levels=50, cmap="plasma", zorder=2)
    ax.contourf(r_grid, -z_grid, grid, levels=50, cmap="plasma", zorder=2)

    ax.figure.colorbar(
        upper,
        ax=ax,
        label="Fusion Rate Density [reactions/m³/sec]",
        location="left",
        anchor=(-0.25, 0.5),
    )

    ax.set_xlabel("R [m]")
    ax.set_xlim(rmajor - 1.2 * rminor, rmajor + 1.2 * rminor)
    ax.set_ylim(-1.2 * rminor * kappa, 1.2 * kappa * rminor)
    ax.set_ylabel("Z [m]")
    ax.plot(
        rmajor,
        0,
        marker="o",
        color="red",
        markersize=6,
        markeredgecolor="black",
        zorder=100,
    )
    # enable minor ticks and grid for clearer reading
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.5, zorder=1)
    # make minor ticks visible on all sides and draw ticks inward for compact look
    ax.tick_params(which="both", direction="in", top=True, right=True)

    # draw contours at % of the DT peak value (both top and mirrored bottom)
    peak = np.nanmax(grid)
    if peak > 0:
        c_kwargs = {
            "levels": [f * peak for f in fractions],
            "colors": colours,
            "linewidths": 1.5,
        }
        # distinct colours for each level

        # top and mirrored bottom contours (no clabel calls — keep only legend)
        ax.contour(r_grid, z_grid, grid, **c_kwargs)
        ax.contour(r_grid, -z_grid, grid, **c_kwargs)

        # create legend entries (use Line2D proxies so we get one entry per requested
        # level)
        legend_handles = [mpl.lines.Line2D([0], [0], color=c, lw=2) for c in colours]
        legend_labels = ["25% peak", "50% peak", "75% peak"]
        ax.legend(legend_handles, legend_labels, loc="upper right", fontsize=8)


def plot_magnetic_fields_in_plasma(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot magnetic field profiles inside the plasma boundary"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    # Get toroidal magnetic field profile (in Tesla)
    b_plasma_toroidal_profile = [
        mfile.get(f"b_plasma_toroidal_profile{i}", scan=scan)
        for i in range(2 * n_plasma_profile_elements)
    ]

    # Get major and minor radius for x-axis in metres
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)

    # Plot magnetic field first (background)
    axis.plot(
        np.linspace(rmajor - rminor, rmajor + rminor, len(b_plasma_toroidal_profile)),
        b_plasma_toroidal_profile,
        color="blue",
        label="Toroidal B-field [T]",
        linewidth=2,
    )

    # Plot plasma on top of magnetic field, displaced vertically by bt
    plot_plasma(axis, mfile, scan, colour_scheme=1)

    # Plot plasma centre dot
    axis.plot(rmajor, 0, marker="o", color="red", markersize=8, label="Plasma Centre")

    v_kwargs = {"color": "green", "linestyle": "--", "linewidth": 1.0}

    # Plot vertical lines at plasma edge
    axis.axvline(rmajor - rminor, **v_kwargs)
    axis.axvline(rmajor + rminor, **v_kwargs)

    h_kwargs = {"color": "blue", "linestyle": "--", "linewidth": 1.0}

    # Plot horizontal line for toroidal magnetic field at plasma inboard
    axis.axhline(mfile.get(f"b_plasma_toroidal_profile{0}", scan=scan), **h_kwargs)

    # Plot horizontal line for toroidal magnetic field at plasma centre
    axis.axhline(mfile.get("b_plasma_toroidal_on_axis", scan=scan), **h_kwargs)

    # Plot horizontal line for toroidal magnetic field at plasma outboard
    axis.axhline(b_plasma_toroidal_profile[-1], **h_kwargs)

    # Text box for inboard toroidal field
    draw_text(
        axis,
        0.1,
        0.025,
        f"$B_{{\\text{{T,inboard}}}}={mfile.get('b_plasma_inboard_toroidal', scan=scan):.2f}$ T\n"  # noqa: E501
        f"$B_{{\\text{{total,inboard}}}}={mfile.get('b_plasma_inboard_total', scan=scan):.2f}$ T",  # noqa: E501
        verticalalignment="center",
        horizontalalignment="center",
        transform=axis.transAxes,
        bbox=box_style("wheat"),
    )

    # Text box for outboard toroidal field
    draw_text(
        axis,
        0.9,
        0.1,
        f"$B_{{\\text{{T,outboard}}}}={mfile.get('b_plasma_outboard_toroidal', scan=scan):.2f}$ T\n"  # noqa: E501
        f"$B_{{\\text{{total,outboard}}}}={mfile.get('b_plasma_outboard_total', scan=scan):.2f}$ T",  # noqa: E501
        verticalalignment="center",
        horizontalalignment="center",
        transform=axis.transAxes,
        bbox=box_style("wheat"),
    )

    axis.set_xlabel("Radial Position [m]")
    axis.set_ylabel("Toroidal Magnetic Field [T]")
    axis.set_title("Toroidal Magnetic Field Profile in Plasma")
    axis.minorticks_on()
    # Enable grid for both major and minor ticks
    axis.grid(which="both", linestyle="--", alpha=0.5)
    axis.grid(which="minor", linestyle=":", alpha=0.3)
    axis.legend(loc="lower right")
    axis.set_xlim(rmajor - 1.25 * rminor, rmajor + 1.25 * rminor)


def plot_plasma_outboard_toroidal_ripple_map(fig, mfile: MFile, scan: int):
    """Plot plasma outboard toroidal ripple map"""
    r_tf_outboard_mid = mfile.get("r_tf_outboard_mid", scan=scan)
    n_tf_coils = mfile.get("n_tf_coils", scan=scan)
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    r_tf_wp_inboard_inner = mfile.get("r_tf_wp_inboard_inner", scan=scan)
    r_tf_wp_inboard_centre = mfile.get("r_tf_wp_inboard_centre", scan=scan)
    r_tf_wp_inboard_outer = mfile.get("r_tf_wp_inboard_outer", scan=scan)
    dx_tf_wp_primary_toroidal = mfile.get("dx_tf_wp_primary_toroidal", scan=scan)
    i_tf_shape = mfile.get("i_tf_shape", scan=scan)
    i_tf_sup = mfile.get("i_tf_sup", scan=scan)
    dx_tf_wp_insulation = mfile.get("dx_tf_wp_insulation", scan=scan)
    dx_tf_wp_insertion_gap = mfile.get("dx_tf_wp_insertion_gap", scan=scan)
    ripple_b_tf_plasma_edge_max = mfile.get("ripple_b_tf_plasma_edge_max", scan=scan)
    i_tf_wp_geom = round(mfile.get("i_tf_wp_geom", scan=scan))

    build = Build()

    r_nom = r_tf_outboard_mid
    dx_nom = dx_tf_wp_primary_toroidal if dx_tf_wp_primary_toroidal is not None else 0.0

    # Simple ±20% scan around nominal values for r and dx
    r_min = r_nom * 0.9
    r_max = r_nom * 1.1

    if dx_nom > 0:
        dx_min = dx_nom * 0.8
        dx_max = dx_nom * 1.2
    else:
        # fallback sensible small range if nominal is zero
        dx_min = 1e-3
        dx_max = 1e-2

    n_r = 50
    n_dx = 50
    r_vals = np.linspace(r_min, r_max, n_r)
    dx_vals = np.linspace(dx_min, dx_max, n_dx)

    rg, dxg = np.meshgrid(r_vals, dx_vals)

    # prepare metric array to hold ripple metric for each (r, dx) pair
    metric = np.full(rg.shape, np.nan, dtype=float)

    for ii in range(rg.shape[0]):
        for jj in range(rg.shape[1]):
            r_test = float(rg[ii, jj])
            dx_test = float(dxg[ii, jj])

            try:
                rip, _, _ = build.plasma_outboard_edge_toroidal_ripple(
                    ripple_b_tf_plasma_edge_max=0.05,
                    r_tf_outboard_mid=r_test,
                    n_tf_coils=int(n_tf_coils),
                    rmajor=rmajor,
                    rminor=rminor,
                    r_tf_wp_inboard_inner=r_tf_wp_inboard_inner,
                    r_tf_wp_inboard_centre=r_tf_wp_inboard_centre,
                    r_tf_wp_inboard_outer=r_tf_wp_inboard_outer,
                    dx_tf_wp_primary_toroidal=dx_test,
                    i_tf_shape=i_tf_shape,
                    i_tf_sup=i_tf_sup,
                    dx_tf_wp_insulation=dx_tf_wp_insulation,
                    dx_tf_wp_insertion_gap=dx_tf_wp_insertion_gap,
                    i_tf_wp_geom=i_tf_wp_geom,
                )
            except (ValueError, ZeroDivisionError, OverflowError, TypeError):
                # Only catch expected numeric/validation errors from the ripple
                # calculation;
                # let other exceptions propagate so they can be diagnosed.
                rip = np.nan
            metric[ii, jj] = rip

    # Create two subplots that share the same x axis
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

    # Make contour plot of the ripple metric (r vs dx) on ax1
    if np.all(np.isnan(metric)):
        ax1.text(
            0.5,
            0.5,
            "No valid ripple data (r vs dx)",
            ha="center",
            va="center",
        )
    else:
        vmin = np.nanmin(metric)
        vmax = np.nanmax(metric)

        # Guard against degenerate range
        if np.isclose(vmin, vmax, atol=1e-12) or np.isnan(vmin) or np.isnan(vmax):
            vmin -= 0.25
            vmax += 0.25

        # Smooth filled contour levels
        levels = np.linspace(vmin, vmax, 50)
        cf = ax1.contourf(rg, dxg, metric, levels=levels, cmap="plasma", extend="both")

        # Contour lines only at 0.5 increments
        step = 0.5
        start = np.floor(vmin / step) * step
        end = np.ceil(vmax / step) * step
        contour_levels = np.arange(start, end + 1e-12, step)

        # Fallback if contour_levels is empty for some reason
        if contour_levels.size < 2:
            contour_levels = np.array([vmin, vmax])

        contours = ax1.contour(
            rg,
            dxg,
            metric,
            levels=contour_levels,
            colors="k",
            linewidths=0.5,
            alpha=0.7,
        )
        ax1.clabel(contours, inline=True, fontsize=8, fmt="%.2f%%", colors="white")
        # Overlay contour line at the specified target ripple value

        target = float(ripple_b_tf_plasma_edge_max)

        if target is not None and not np.isnan(target):
            # Check if target lies within computed metric range
            if (target >= vmin) and (target <= vmax):
                c_target = ax1.contour(
                    rg,
                    dxg,
                    metric,
                    levels=[target],
                    colors="white",
                    linewidths=2.0,
                    linestyles="--",
                    zorder=20,
                )
                ax1.clabel(
                    c_target,
                    inline=True,
                    fmt={target: f"Input Max {target:.2f}%"},
                    fontsize=8,
                    colors="white",
                )
            else:
                # annotate that target is outside plotted range
                ax1.text(
                    0.02,
                    0.98,
                    f"Target ripple {target:.2f}% outside plot range"
                    f" [{vmin:.2f},{vmax:.2f}]",
                    transform=ax1.transAxes,
                    color="white",
                    fontsize=8,
                    va="top",
                    bbox={"facecolor": "black", "alpha": 0.6, "pad": 2},
                )

        # Colourbar with 0.5 increments (use the same contour_levels as for the contour
        # lines)
        ticks = contour_levels
        # Fallback to sensible ticks if contour_levels is not appropriate
        if ticks.size == 0 or np.isnan(ticks).all():
            ticks = np.linspace(vmin, vmax, 5)
        cb = ax1.figure.colorbar(
            cf, ax=ax1, label="Plasma Outboard Toroidal Ripple", ticks=ticks
        )
        cb.ax.set_yticklabels([f"{t:.2f}%" for t in ticks])

        # mark nominal point
        ax1.scatter(
            [r_nom],
            [dx_nom],
            color="white",
            edgecolor="black",
            s=200,
            linewidths=1.5,
            marker="o",
            zorder=10,
            label="Design Point",
        )
        ax1.set_xlabel("Outboard TF leg centre [m]")
        ax1.set_ylabel("WP Toroidal Width [m]")
        ax1.legend(loc="upper right")

    # ---------------------------------------------------------------------
    # Second plot: scan number of TF coils vs r_tf_outboard_mid (keep dx at nominal)
    # ---------------------------------------------------------------------
    # Determine a sensible integer range of TF coils to scan around nominal
    n_nom = int(n_tf_coils)
    span = max(2, int(min(12, n_nom // 2)))  # choose a span based on nominal
    n_min = max(10, n_nom - span)
    n_max = n_nom + span
    n_vals = np.arange(n_min, n_max + 1, dtype=int)

    n_r2 = 60
    r_vals2 = np.linspace(r_min, r_max, n_r2)
    rg2, ng2 = np.meshgrid(r_vals2, n_vals)

    metric2 = np.full(rg2.shape, np.nan, dtype=float)

    for ii in range(rg2.shape[0]):
        for jj in range(rg2.shape[1]):
            r_test = float(rg2[ii, jj])
            n_test = int(ng2[ii, jj])
            try:
                rip, *_ = build.plasma_outboard_edge_toroidal_ripple(
                    ripple_b_tf_plasma_edge_max=0.05,
                    r_tf_outboard_mid=r_test,
                    n_tf_coils=n_test,
                    rmajor=rmajor,
                    rminor=rminor,
                    r_tf_wp_inboard_inner=r_tf_wp_inboard_inner,
                    r_tf_wp_inboard_centre=r_tf_wp_inboard_centre,
                    r_tf_wp_inboard_outer=r_tf_wp_inboard_outer,
                    dx_tf_wp_primary_toroidal=dx_nom,
                    i_tf_shape=i_tf_shape,
                    i_tf_sup=i_tf_sup,
                    dx_tf_wp_insulation=dx_tf_wp_insulation,
                    dx_tf_wp_insertion_gap=dx_tf_wp_insertion_gap,
                    i_tf_wp_geom=i_tf_wp_geom,
                )
            except (ValueError, ZeroDivisionError, OverflowError, TypeError):
                # Only catch expected numeric/validation errors from the ripple
                # calculation;
                # let other exceptions propagate so they can be diagnosed.
                rip = np.nan
            metric2[ii, jj] = rip

    # Plot the second metric on the bottom axes (ax2) so it shares x-axis with ax1
    if np.all(np.isnan(metric2)):
        ax2.text(
            0.5,
            0.5,
            "No valid ripple data (r vs n_tf_coils)",
            ha="center",
            va="center",
        )
    else:
        vmin2 = np.nanmin(metric2)
        vmax2 = np.nanmax(metric2)

        # filled contour levels (smooth shading)
        levels2 = np.linspace(vmin2, vmax2, 40)
        cf2 = ax2.contourf(
            rg2, ng2, metric2, levels=levels2, cmap="viridis", extend="both"
        )

        # contour lines only at 0.5 steps
        step = 0.5
        start = np.floor(vmin2 / step) * step
        end = np.ceil(vmax2 / step) * step
        contour_levels = np.arange(start, end + 1e-12, step)

        # fallback if arange returned empty (very small range)
        if contour_levels.size == 0:
            contour_levels = np.array([vmin2, vmax2])

        contours2 = ax2.contour(
            rg2,
            ng2,
            metric2,
            levels=contour_levels,
            colors="k",
            linewidths=0.5,
            alpha=0.7,
        )
        ax2.clabel(contours2, inline=True, fontsize=8, fmt="%.2f%%", colors="white")

        target2 = float(ripple_b_tf_plasma_edge_max)

        if target2 is not None and not np.isnan(target2):
            if (target2 >= vmin2) and (target2 <= vmax2):
                c_target2 = ax2.contour(
                    rg2,
                    ng2,
                    metric2,
                    levels=[target2],
                    colors="white",
                    linewidths=2.0,
                    linestyles="--",
                    zorder=20,
                )
            ax2.clabel(
                c_target2,
                inline=True,
                fmt={target2: f"Input Max {target2:.2f}%"},
                fontsize=8,
                colors="white",
            )
        else:
            ax2.text(
                0.02,
                0.98,
                f"Target ripple {target2:.2f}% outside plot range"
                f" [{vmin2:.2f},{vmax2:.2f}]",
                transform=ax2.transAxes,
                color="white",
                fontsize=8,
                va="top",
                bbox={"facecolor": "black", "alpha": 0.6, "pad": 2},
            )
        # colorbar with 0.5 increments
        # ensure contour_levels exists and is in 0.5 steps (constructed above)
        ticks = contour_levels
        cb2 = ax2.figure.colorbar(
            cf2, ax=ax2, label="Plasma Outboard Toroidal Ripple", ticks=ticks
        )
        cb2.ax.set_yticklabels([f"{t:.2f}%" for t in ticks])

        # nominal markers
        ax2.scatter(
            [r_nom],
            [n_nom],
            color="white",
            edgecolor="black",
            s=300,
            linewidths=1.5,
            marker="o",
            zorder=10,
            label="Design Point",
        )
        ax2.set_xlabel("Outboard TF leg centre [m]")
        ax2.set_ylabel("Number of TF coils")
        ax2.set_yticks(n_vals)
        ax2.legend(loc="upper right")

    # Improve layout
    fig.tight_layout()


def plot_plasma_coloumb_logarithms(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the plasma coloumb logarithms on the given axis."""
    plasma_coulomb_log_electron_electron_profile = [
        mfile_data.data[f"plasma_coulomb_log_electron_electron_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    plasma_coulomb_log_electron_deuteron_profile = [
        mfile_data.data[f"plasma_coulomb_log_electron_deuteron_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    plasma_coulomb_log_electron_triton_profile = [
        mfile_data.data[f"plasma_coulomb_log_electron_triton_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    plasma_coulomb_log_deuteron_triton_profile = [
        mfile_data.data[f"plasma_coulomb_log_deuteron_triton_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    plasma_coulomb_log_electron_alpha_thermal_profile = [
        mfile_data.data[
            f"plasma_coulomb_log_electron_alpha_thermal_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(plasma_coulomb_log_electron_electron_profile)),
        plasma_coulomb_log_electron_electron_profile,
        color="blue",
        linestyle="-",
        label=r"$ln \Lambda_{e-e}$",
    )

    axis.plot(
        np.linspace(0, 1, len(plasma_coulomb_log_electron_deuteron_profile)),
        plasma_coulomb_log_electron_deuteron_profile,
        color="pink",
        linestyle="-",
        label=r"$ln \Lambda_{e-D}$",
    )

    axis.plot(
        np.linspace(0, 1, len(plasma_coulomb_log_electron_triton_profile)),
        plasma_coulomb_log_electron_triton_profile,
        color="green",
        linestyle="-",
        label=r"$ln \Lambda_{e-T}$",
    )

    axis.plot(
        np.linspace(0, 1, len(plasma_coulomb_log_deuteron_triton_profile)),
        plasma_coulomb_log_deuteron_triton_profile,
        color="orange",
        linestyle="-",
        label=r"$ln \Lambda_{D-T}$",
    )

    axis.plot(
        np.linspace(0, 1, len(plasma_coulomb_log_electron_alpha_thermal_profile)),
        plasma_coulomb_log_electron_alpha_thermal_profile,
        color="red",
        linestyle="-",
        label=r"$ln \Lambda_{e-\alpha,thermal}$",
    )

    axis.set_ylabel("Coulomb Logarithm")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


__all__ = [
    "plot_magnetic_fields_in_plasma",
    "plot_max_normalised_beta_comparison",
    "plot_plasma",
    "plot_plasma_coloumb_logarithms",
    "plot_plasma_current_comparison",
    "plot_plasma_outboard_toroidal_ripple_map",
    "reaction_plot_grid",
]
