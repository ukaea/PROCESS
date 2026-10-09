"""Power Flow functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

from process.core.io.plot.summary.common import (
    box_style,
    load_plot_image,
    setup_axis,
    text_layout,
)
from process.core.io.plot.summary.reporting.text import plot_info

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_main_power_flow(axis: plt.Axes, mfile: MFile, scan: int, fig: plt.Figure):
    """Plots the main power flow diagram for the fusion reactor, including plasma,
    heating and current drive,
    first wall, blanket, vacuum vessel, divertor, coolant pumps, turbine, generator, and
    auxiliary systems.
    Annotates the diagram with power values and draws arrows to indicate power flows.

    Parameters
    ----------
    axis:
        The matplotlib axis object to plot on.
    mfile:
        The MFILE data object containing power flow parameters.
    scan:
        The scan number to use for extracting data.
    fig:
        The matplotlib figure object for additional annotations.
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
    plasma = load_plot_image("plasma.png")

    # Display the plasma image over the figure, not the axes
    new_ax = axis.inset_axes(
        (-0.15, 0.6, 0.45, 0.45), transform=axis.transAxes, zorder=1
    )
    new_ax.imshow(plasma)
    new_ax.axis("off")

    # Add fusion power to plasma
    axis.text(
        0.22,
        0.75,
        f"$P_{{{{fus}}}}$\n{mfile.get('p_fusion_total_mw', scan=scan):.2f} MW",
        transform=fig.transFigure,
        horizontalalignment="left",
        verticalalignment="bottom",
        zorder=2,
        fontsize=11,
    )
    # Load the neutron image
    neutron = load_plot_image("neutron.png")

    new_ax = axis.inset_axes(
        (0.2, 0.85, 0.03, 0.03), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(neutron)
    new_ax.axis("off")

    # Add lost alpha power
    axis.text(
        0.22,
        0.81,
        f"$P_{{\\alpha,{{loss}}}}$\n{mfile.get('p_fw_alpha_mw', scan=scan):,.2f} MW",
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
        f"$P_{{{{rad}}}}$\n{mfile.get('p_plasma_rad_mw', scan=scan):,.2f} MW",
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
        f"$P_{{{{sep}}}}$\n{mfile.get('p_plasma_separatrix_mw', scan=scan):,.2f} MW",
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
            "arrowstyle": ("-|>,head_length=1,head_width=0.3"),  # solid filled head
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
    alpha = load_plot_image("alpha_particle.png")

    # Display the alpha particle image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.16, 0.95, 0.025, 0.025), transform=axis.transAxes, zorder=10
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
        f"$P_{{\\text{{neutron}}}}$:\n{mfile.get('p_neutron_total_mw', scan=scan):,.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("grey", alpha=0.8),
    )

    # ===========================================

    # =========================================
    # Heating and current drive systems
    # =========================================

    # Add HCD primary injected power
    axis.text(
        0.0725,
        0.83,
        "$P_{\\text{HCD,primary}}$:"
        f" {mfile.get('p_hcd_primary_injected_mw', scan=scan) + mfile.get('p_hcd_primary_extra_heat_mw', scan=scan):.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # Add HCD secondary injected power
    axis.text(
        0.0725,
        0.725,
        "$P_{\\text{HCD,secondary}}$:"
        f" {mfile.get('p_hcd_secondary_injected_mw', scan=scan) + mfile.get('p_hcd_secondary_extra_heat_mw', scan=scan):.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # Load the HCD injector image
    hcd_injector_1 = hcd_injector_2 = load_plot_image("hcd_injector.png")

    # Display the injector image over the figure, not the axes
    new_ax = axis.inset_axes(
        (-0.2, 0.8, 0.15, 0.15), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(hcd_injector_1)
    new_ax.axis("off")
    new_ax = axis.inset_axes((-0.2, 0.5, 0.15, 0.5), transform=axis.transAxes, zorder=10)
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
        **text_layout(fig),
        bbox=box_style("lightyellow"),
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
        f"$P_{{\\text{{secondary,loss}}}}$:\n{mfile.get('p_hcd_secondary_electric_mw', scan=scan) * (1.0 - mfile.get('eta_hcd_secondary_injector_wall_plug', scan=scan)):.2f} MWe",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lightblue", linestyle="dashed"),
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
        f"$P_{{\\text{{primary,loss}}}}$:\n{mfile.get('p_hcd_primary_electric_mw', scan=scan) * (1.0 - mfile.get('eta_hcd_primary_injector_wall_plug', scan=scan)):.2f} MWe",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lightblue", linestyle="dashed"),
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
        f"$P_{{\\text{{secondary}}}}$:\n{mfile.get('p_hcd_secondary_electric_mw', scan=scan):.2f}"  # noqa: E501
        " MWe\n$\\eta$:"
        f" {mfile.get('eta_hcd_secondary_injector_wall_plug', scan=scan):.2f}",
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # Plot HCD primary electric box
    axis.text(
        0.025,
        0.35,
        f"$P_{{\\text{{primary}}}}$:\n{mfile.get('p_hcd_primary_electric_mw', scan=scan):.2f}"  # noqa: E501
        " MWe\n$\\eta$:"
        f" {mfile.get('eta_hcd_primary_injector_wall_plug', scan=scan):.2f}",
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # =============================================

    # =============================================
    # Low grade heat total
    # =============================================

    # Plot box of total low grade secondary heat
    axis.text(
        0.325,
        0.225,
        "\n\nTotal Low Grade Secondary Heat\n\n"
        f" {mfile.get('p_plant_secondary_heat_mw', scan=scan):,.2f} MWth",
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="center",
        transform=fig.transFigure,
        bbox=box_style("lightblue", linestyle="dashed"),
        zorder=4,
    )

    # =============================================

    # ==========================================
    # Power conversion systems
    # ===========================================

    # Load the turbine image
    turbine = load_plot_image("turbine.png")

    # Display the turbine image over the figure, not the axes
    new_ax = axis.inset_axes((1.1, 0.0, 0.15, 0.15), transform=axis.transAxes, zorder=10)
    new_ax.imshow(turbine)
    new_ax.axis("off")

    # Plot the total primary thermal power box
    axis.text(
        0.9,
        0.25,
        f"$P_{{\\text{{primary,thermal}}}}$:\n{mfile.get('p_plant_primary_heat_mw', scan=scan):,.2f}"  # noqa: E501
        " MW\n$\\eta_{\\text{turbine}}$:"
        f" {mfile.get('eta_turbine', scan=scan):.3f}",
        **text_layout(fig),
        bbox=box_style("orange"),
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
    generator = load_plot_image("generator.png")

    # Display the generator image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.96, 0.0, 0.15, 0.15), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(generator)
    new_ax.axis("off")

    # Generator to gross electric power
    axis.annotate(
        "",
        xy=(0.745, 0.17),
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
        **text_layout(fig),
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
    pylon = load_plot_image("pylon.png")

    # Display the pylon image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.925, -0.1, 0.1, 0.1), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(pylon)
    new_ax.axis("off")

    # Plot the gross electric power box
    axis.text(
        0.68,
        0.15,
        f"$P_{{\\text{{gross}}}}$:\n{mfile.get('p_plant_electric_gross_mw', scan=scan):,.2f} MWe",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lime"),
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
        f"$P_{{\\text{{loss}}}}$:\n{mfile.get('p_turbine_loss_mw', scan=scan):,.2f}"
        " MWth",
        **text_layout(fig),
        bbox=box_style("orange", linestyle="dashed"),
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
        f"$P_{{\\text{{net,electric}}}}$:\n{mfile.get('p_plant_electric_net_mw', scan=scan):,.2f} MWe",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lime"),
    )

    # Plot the recirculated electric power box
    axis.text(
        0.575,
        0.14,
        f"$P_{{\\text{{recirc,electric}}}}$:\n{mfile.get('p_plant_electric_recirc_mw', scan=scan):,.2f}"  # noqa: E501
        " MWe\n"
        f"$f_{{\\text{{recirc}}}}$:\n{mfile.get('f_p_plant_electric_recirc', scan=scan):,.2f}",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("lime"),
    )

    # Gross to recirculated power arrow
    axis.annotate(
        "",
        xy=(0.64, 0.17),
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
        xytext=(0.645, 0.185),
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
    fw = load_plot_image("fw.png")

    # Display the first wall image over the figure, not the axes
    new_ax = axis.inset_axes((0.4, 0.625, 0.4, 0.4), transform=axis.transAxes, zorder=10)
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
        "$P_{\\text{FW,"
        f" }}\\alpha}}$:\n{mfile.get('p_fw_alpha_mw', scan=scan):.2f} MW",
        **text_layout(fig),
        bbox=box_style("red"),
    )

    # Neutron power incident on first wall box
    axis.text(
        0.46,
        0.775,
        f"$P_{{\\text{{FW,nuclear}}}}$:\n{mfile.get('p_fw_nuclear_heat_total_mw', scan=scan):,.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("grey", alpha=0.8),
    )

    # Plot radiation power incident on first wall box
    axis.text(
        0.46,
        0.71,
        f"$P_{{\\text{{FW,rad}}}}$:\n{mfile.get('p_fw_rad_total_mw', scan=scan):,.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("dodgerblue", alpha=0.8),
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
        "Primary thermal\n(inc pump):"
        f" {mfile.get('p_fw_heat_deposited_mw', scan=scan):,.2f} MWth",
        **text_layout(fig),
        bbox=box_style("orange"),
    )

    # Blanket heat deposited box
    axis.text(
        0.7,
        0.555,
        "Primary thermal\n(inc pump):"
        f" {mfile.get('p_blkt_heat_deposited_mw', scan=scan):,.2f} MWth",
        **text_layout(fig),
        bbox=box_style("orange"),
    )

    # Shield heat deposited box
    axis.text(
        0.875,
        0.555,
        f"Primary thermal:\n{mfile.get('p_shld_heat_deposited_mw', scan=scan):.2f} MWth",
        **text_layout(fig),
        bbox=box_style("orange"),
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

    # Draw arrow from blanket primary heat box to blanket and FW primary heat deposited
    # box
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

    # Draw a downward arrow from the primary thermal box to the right side of the
    # generator
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
        "Primary thermal (inc pump):"
        f" {mfile.get('p_fw_blkt_heat_deposited_mw', scan=scan):,.2f} MWth\n",
        **text_layout(fig),
        bbox=box_style("orange"),
    )

    # Load the blanket image
    blanket = load_plot_image("blanket_with_coolant.png")

    # Display the blanket image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.75, 0.625, 0.4, 0.4), transform=axis.transAxes, zorder=10
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
        f"$P_{{\\text{{Blkt,nuclear}}}}$:\n{mfile.get('p_blkt_nuclear_heat_total_mw', scan=scan):,.2f}"  # noqa: E501
        " MW\n"
        f"$P_{{\\text{{Blkt,multiplication}}}}$:\n{mfile.get('p_blkt_multiplication_mw', scan=scan):,.2f}"  # noqa: E501
        " MW\n"
        f"$f_{{\\text{{multiplication}}}}$:\n{mfile.get('f_p_blkt_multiplication', scan=scan):,.2f}",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("grey", alpha=0.8),
    )

    # Load the vacuum vessel image
    vv = load_plot_image("vv.png")

    # Display the vacuum vessel image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.975, 0.625, 0.4, 0.4), transform=axis.transAxes, zorder=10
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
        f"$P_{{\\text{{shld,secondary}}}}$:\n{mfile.get('p_shld_secondary_heat_mw', scan=scan):,.2f}"  # noqa: E501
        " MWth",
        **text_layout(fig),
        bbox=box_style("lightblue", linestyle="dashed"),
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
    divertor = load_plot_image("divertor.png")

    # Display the divertor image over the figure, not the axes
    new_ax = axis.inset_axes((0.1, 0.4, 0.3, 0.25), transform=axis.transAxes, zorder=10)
    new_ax.imshow(divertor)
    new_ax.axis("off")

    # Total divertor radiation power box
    axis.text(
        0.29,
        0.57,
        f"$P_{{\\text{{div,rad}}}}$:\n{mfile.get('p_div_rad_total_mw', scan=scan):,.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("dodgerblue", alpha=0.8),
    )

    # Divertor nuclear heat total box
    axis.text(
        0.4,
        0.58,
        f"$P_{{\\text{{div,nuclear}}}}$:\n{mfile.get('p_div_nuclear_heat_total_mw', scan=scan):,.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("grey", alpha=0.8),
    )

    # Divertor primary thermal heat deposited box
    axis.text(
        0.44,
        0.46,
        "Primary thermal (inc"
        f" pump):\n{mfile.get('p_div_heat_deposited_mw', scan=scan):.2f}"
        " MWth\nSolid angle fraction:"
        f" {mfile.get('f_ster_div_single', scan=scan):.3f}\nPrimary heat"
        f" fraction: {mfile.get('f_p_div_primary_heat', scan=scan):.3f}",
        **text_layout(fig),
        bbox=box_style("orange"),
        zorder=100,
    )

    # Divertor secondary heat box
    axis.text(
        0.3,
        0.375,
        f"$P_{{\\text{{div,secondary}}}}$:\n{mfile.get('p_div_secondary_heat_mw', scan=scan):.2f}"  # noqa: E501
        " MWth",
        **text_layout(fig),
        bbox=box_style("lightblue", linestyle="dashed"),
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
        "$P_{\\text{div,pump}}$:"
        f" {mfile.get('p_div_coolant_pump_mw', scan=scan):.2f} MW",
        **text_layout(fig),
        bbox=box_style("wheat", alpha=0.8),
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
        f"$P_{{\\text{{shld,pump}}}}$:\n{mfile.get('p_shld_coolant_pump_mw', scan=scan):.2f} MW",  # noqa: E501
        **text_layout(fig),
        bbox=box_style("wheat", alpha=0.8),
    )

    # FW and Blanket coolant pumps total
    axis.text(
        0.725,
        0.4,
        "$P_{\\text{FW +"
        f" Blkt}}}}$:\n{mfile.get('p_fw_blkt_coolant_pump_mw', scan=scan):.2f} MW",
        **text_layout(fig),
        bbox=box_style("wheat", alpha=0.8),
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
        "Coolant pumps"
        f" electric:\n{mfile.get('p_coolant_pump_elec_total_mw', scan=scan):.3f}"
        " MWe\n$\\eta$:"
        f" {mfile.get('eta_coolant_pump_electric', scan=scan):.3f}",
        **text_layout(fig),
        bbox=box_style("lime", alpha=0.8),
    )

    # Coolant pumps total
    axis.text(
        0.7,
        0.325,
        "Coolant pumps"
        f" total:\n{mfile.get('p_coolant_pump_total_mw', scan=scan):.3f} MW",
        **text_layout(fig),
        bbox=box_style("wheat", alpha=0.8),
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
        "Coolant pumps losses"
        f" total:\n{mfile.get('p_coolant_pump_loss_total_mw', scan=scan):.3f}"
        " MWth",
        **text_layout(fig),
        bbox=box_style("lightblue", alpha=0.8, linestyle="dashed"),
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
        0.49,
        0.05,
        f"Cryo Plant:\n{mfile.get('p_cryo_plant_electric_mw', scan=scan):.3f} MWe",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # Recirculated power to cryo plant arrow
    axis.annotate(
        "",
        xy=(0.525, 0.075),
        xytext=(0.525, 0.1625),
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
        0.4,
        0.05,
        f"Tritium Plant:\n{mfile.get('p_tritium_plant_electric_mw', scan=scan):.3f} MWe",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # # Recirculated power to tritium plant arrow
    axis.annotate(
        "",
        xy=(0.44, 0.075),
        xytext=(0.44, 0.1625),
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
        0.575,
        0.05,
        f"Vacuum pumps:\n{mfile.get('vachtmw', scan=scan):.3f} MWe",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # Recirculated power to vacuum pumps arrow
    axis.annotate(
        "",
        xy=(0.62, 0.08),
        xytext=(0.62, 0.1375),
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
        0.085,
        0.075,
        "Plant base"
        f" load:\n{mfile.get('p_plant_electric_base_total_mw', scan=scan):.3f}"
        " MWe\nMinimum base"
        f" load:\n{mfile.get('p_plant_electric_base', scan=scan) * 1.0e-6:.3f}"
        " MWe\nPlant floor power"
        f" density:\n{mfile.get('pflux_plant_floor_electric', scan=scan) * 1.0e-3:.3f}"
        " kW$\\text{m}^{-2}$",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # TF coil power box
    axis.text(
        0.325,
        0.075,
        f"TF coils:\n{mfile.get('p_tf_electric_supplies_mw', scan=scan):.3f} MWe",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # PF coil power box
    axis.text(
        0.25,
        0.05,
        f"PF coils:\n{mfile.get('p_pf_electric_supplies_mw', scan=scan):.3f} MWe",
        **text_layout(fig),
        bbox=box_style("burlywood", alpha=0.8),
    )

    # Recirculated power to TF,PF and plant base arrow bend
    axis.annotate(
        "",
        xy=(0.22, 0.16),
        xytext=(0.574, 0.16),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 1.5,
            "zorder": 5,
            "fill": True,
        },
    )

    # Recirculated power to  PF
    axis.annotate(
        "",
        xy=(0.28, 0.075),
        xytext=(0.28, 0.1625),
        xycoords=fig.transFigure,
        arrowprops={
            "arrowstyle": "-|>,head_length=1,head_width=0.3",
            "color": "black",
            "linewidth": 2.0,
            "zorder": 5,
            "fill": True,
        },
    )

    # Recirculated power to TF
    axis.annotate(
        "",
        xy=(0.35, 0.1),
        xytext=(0.35, 0.1625),
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
        f"$P_{{\\text{{HCD,loss}}}}$:\n{mfile.get('p_hcd_secondary_heat_mw', scan=scan):.2f}"  # noqa: E501
        " MWth",
        **text_layout(fig),
        bbox=box_style("lightblue", alpha=0.8, linestyle="dashed"),
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
        f"$P_{{\\text{{TF,nuclear}}}}$:\n{mfile.get('p_tf_nuclear_heat_mw', scan=scan):.2f}"  # noqa: E501
        " MWth",
        **text_layout(fig),
        bbox=box_style("lightblue", linestyle="dashed"),
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


def plot_power_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot power info

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    axis.text(-0.05, 1, "Power flows:", ha="left", va="center")
    setup_axis(axis, xmin=0, xmax=1, ymin=-16, ymax=1)

    gross_eff = 100.0 * (
        mfile.get("p_plant_electric_gross_mw", scan=scan)
        / mfile.get("p_plant_primary_heat_mw", scan=scan)
    )

    net_eff = 100.0 * (
        (
            mfile.get("p_plant_electric_gross_mw", scan=scan)
            - mfile.get("p_coolant_pump_elec_total_mw", scan=scan)
        )
        / (
            mfile.get("p_plant_primary_heat_mw", scan=scan)
            - mfile.get("p_coolant_pump_elec_total_mw", scan=scan)
        )
    )

    plant_eff = 100.0 * (
        mfile.get("p_plant_electric_net_mw", scan=scan)
        / mfile.get("p_fusion_total_mw", scan=scan)
    )

    # Define appropriate pedestal and impurity parameters
    coredescription = (
        "radius_plasma_core_norm",
        "Normalised radius of 'core' region",
        "",
    )
    if mfile.get("i_plasma_pedestal", scan=scan) == 1:
        ped_height = (
            "nd_plasma_pedestal_electron",
            "Electron density at pedestal",
            "m$^{-3}$",
        )
        ped_pos = (
            "radius_plasma_pedestal_density_norm",
            "r/a at density pedestal",
            "",
        )
    else:
        ped_height = ("", "No pedestal model used", "")
        ped_pos = ("", "", "")

    p_cryo_plant_electric_mw = mfile.get("p_cryo_plant_electric_mw", scan=scan)

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
        ("life_div_fpy", "Divertor life", "years"),
        ("p_plant_primary_heat_mw", "Primary (high grade) heat", "MW"),
        (gross_eff, "Gross cycle efficiency", "%"),
        (net_eff, "Net cycle efficiency", "%"),
        ("p_plant_electric_gross_mw", "Gross electric power", "MW"),
        ("p_plant_electric_net_mw", "Net electric power", "MW"),
        (
            plant_eff,
            (
                r"Fusion-to-electric efficiency"
                r" $\frac{P_{\mathrm{e,net}}}{P_{\mathrm{fus}}}$"
            ),
            "%",
        ),
    ]

    plot_info(axis, data, mfile, scan)


def plot_blanket_coolant_properties(fig: plt.Figure, m_file: MFile, scan: int):
    """Combined plot of blanket coolant channel structure and properties."""
    for side, x_position in (("inboard", 0.1), ("outboard", 0.5)):

        def get(variable: str):
            return m_file.get(variable, scan=scan)

        text = (
            f"$\\mathbf{{{side.capitalize()} \\ blanket:}}$\n \n"
            "Radius of blanket channel: "
            f"{m_file.get('radius_blkt_channel', scan=scan):.4f} m\n"
            "Channel roughness ($\\epsilon$): "
            f"{m_file.get('roughness_fw_channel', scan=scan):.4e} m\n\n"
            "Radial coolant channel length: "
            f"{get(f'len_blkt_{side}_coolant_channel_radial'):.4f} m\n"
            "Poloidal coolant channel length: "
            f"{get(f'len_blkt_{side}_segment_poloidal'):.4f} m\n"
            "Number of radial channels: "
            f"{get(f'n_blkt_{side}_module_coolant_sections_radial')}\n"
            "Number of poloidal channels: "
            f"{get(f'n_blkt_{side}_module_coolant_sections_poloidal')}\n"
            "Total length of coolant channel straight sections: "
            f"{get(f'len_blkt_{side}_channel_total'):.4f} m\n\n"
            "Pressure drop for straight sections: "
            f"{get(f'dpres_blkt_{side}_coolant_channel_straight_total'):,.2f} Pa\n"
            "Pressure drop for 90° bends: "
            f"{get(f'dpres_blkt_{side}_coolant_channel_90_bend'):,.2f} Pa\n"
            "Total pressure drop for 90° bends: "
            f"{get(f'dpres_blkt_{side}_coolant_channel_90_bends_total'):,.2f} Pa\n"
            "Pressure drop for 180° bends: "
            f"{get(f'dpres_blkt_{side}_coolant_channel_180_bend'):,.2f} Pa\n"
            "Total pressure drop for 180° bends: "
            f"{get(f'dpres_blkt_{side}_coolant_channel_180_bends_total'):,.2f} Pa\n"
            "Total pressure drop for all bends: "
            f"{get(f'dpres_blkt_{side}_bends_total'):,.2f} Pa\n\n"
            "Reynolds number ($Re$): "
            f"{get(f'reynolds_blkt_{side}_coolant'):,.4f}\n"
            "Darcy Friction factor ($f$): "
            f"{get(f'darcy_frict_blkt_{side}_coolant'):.4f}\n\n"
            "Friction drop coefficient for straight sections: "
            f"{get(f'f_straight_blkt_{side}_coolant'):.4f}\n"
            "Friction drop coefficient for 90° bends: "
            f"{get(f'f_elbow_blkt_{side}_90_bend'):.4f}\n"
            "Friction drop coefficient for 180° bends: "
            f"{get(f'f_elbow_blkt_{side}_180_bend'):.4f}\n\n"
            "Total coolant mass flow rate: "
            f"{get(f'mflow_blkt_{side}_coolant'):.4f} kg/s\n"
            "Coolant mass flow rate in single channel: "
            f"{get(f'mflow_blkt_{side}_coolant_channel'):.4f} kg/s\n"
            "Coolant velocity in single channel: "
            f"{get(f'vel_blkt_{side}_coolant'):.4f} m/s"
        )

        fig.text(
            x_position,
            0.5,
            text,
            **text_layout(fig, v_align="top"),
            bbox=box_style("wheat"),
        )


__all__ = ["plot_blanket_coolant_properties", "plot_main_power_flow", "plot_power_info"]
