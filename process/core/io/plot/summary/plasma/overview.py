"""Plasma functions for PROCESS summary plots."""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING, Literal, TypedDict

import matplotlib.pyplot as plt

from process.core.io.plot.summary.common import (
    box_style,
    load_plot_image,
)
from process.core.io.plot.summary.constants import (
    white_box,
)
from process.core.io.plot.summary.plasma.physics import (
    plot_plasma,
)
from process.core.io.plot.summary.rendering import (
    draw_annotation,
    draw_text,
)
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.geometry.plasma import plasma_geometry
from process.models.physics.bootstrap_current import (
    BootstrapCurrentFractionModel,
)
from process.models.physics.current_drive import CurrentDriveModel
from process.models.physics.density_limit import DensityLimitModel
from process.models.physics.l_h_transition import (
    PlasmaConfinementTransitionModel,
)
from process.models.physics.physics import (
    BetaComponentLimits,
    BetaNormMaxModel,
    IndInternalNormModel,
)
from process.models.physics.plasma_current import (
    PlasmaCurrentModel,
    PlasmaDiamagneticCurrentModel,
)
from process.models.physics.plasma_geometry import (
    PlasmaGeometryModelType,
)

if TYPE_CHECKING:
    from matplotlib.transforms import Transform

    from process.core.io.mfile import MFile


def plot_main_plasma_information(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    fig: plt.Figure,
):
    """Plots the main plasma information including plasma shape, geometry, currents,
    heating,
    confinement, and other relevant plasma parameters.

    Parameters
    ----------
    axis : plt.Axes
        The matplotlib axis object to plot on.
    mfile : MFile
        The MFILE data object containing plasma parameters.
    scan : int
        The scan number to use for extracting data.
    colour_scheme : int
        The colour scheme to use for plots.
    fig : plt.Figure
        The matplotlib figure object for additional annotations.
    """
    # Import key variables
    triang = mfile.get("triang", scan=scan)
    kappa = mfile.get("kappa", scan=scan)

    # Remove the axes
    axis.axis("off")

    # Plot the main plasma shape
    plot_plasma(axis, mfile, scan, colour_scheme)

    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    # Get the plasma permieter points for the core plasma region
    pg = plasma_geometry(
        rmajor=rmajor,
        rminor=mfile.get("rminor", scan=scan)
        * mfile.get("radius_plasma_core_norm", scan=scan),
        triang=mfile.get("triang", scan=scan),
        kappa=mfile.get("kappa", scan=scan),
        i_single_null=mfile.get("i_single_null", scan=scan),
        i_plasma_shape=1,
        square=mfile.get("plasma_square", scan=scan),
    )
    # Plot the core plasma boundary line
    axis.plot(pg.rs, pg.zs, color="black", linestyle="--")

    # Plot the centre of the plasma
    axis.plot(rmajor, 0, "r+", markersize=20, markeredgewidth=2)

    # Add Q plasma information box
    draw_text(
        axis,
        0.725,
        0.175,
        f"$Q_{{\\text{{plasma}}}}$: {mfile.get('big_q_plasma', scan=scan):.2f}",
        fontsize=15,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=white_box,
        transform=fig.transFigure,
    )

    # =========================================

    # Draw a double-ended arrow from the inner plasma edge to the center
    draw_annotation(
        axis,
        "",
        xy=(rmajor - rminor, 0),  # Inner plasma edge
        xytext=(rmajor, 0),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Add a label for the minor radius
    draw_text(
        axis,
        rmajor - rminor / 2,
        -rminor * kappa * 0.08,
        f"$a$: {rminor:.2f} m",
        fontsize=9,
        color="black",
        ha="center",
        bbox=white_box,
    )

    # ============================================

    # Draw a single-ended arrow from the machien centre to the plasma center
    draw_annotation(
        axis,
        "",
        xy=(axis.get_xlim()[0], -rminor * 0.3 * kappa),  # Inner plasma edge
        xytext=(rmajor, -rminor * 0.3 * kappa),  # Center
        arrowprops={"arrowstyle": "<-", "color": "black"},
    )

    # Add a label for the major radius
    draw_text(
        axis,
        rmajor - rminor / 2,
        -rminor * kappa * 0.25,
        f"$R_0$: {rmajor:.2f} m",
        fontsize=9,
        color="black",
        ha="center",
        bbox=white_box,
    )

    # ============================================

    # Draw a double-ended arrow from the xpoint to the center to show elongation
    draw_annotation(
        axis,
        "",
        xy=(rmajor - rminor * triang, kappa * rminor),  # Inner plasma edge
        xytext=(rmajor - rminor * triang, 0),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Write the elongation beside the vertical line, position relative to figure axes
    draw_text(
        axis,
        0.3,
        0.75,
        f"$\\kappa$: {mfile.get('kappa', scan=scan):.2f}",
        fontsize=9,
        color="black",
        rotation=270,
        verticalalignment="center",
        transform=axis.transAxes,
        bbox=white_box,
    )

    # =============================================

    # Draw a double-ended arrow from the inner plasma edge to the center
    draw_annotation(
        axis,
        "",
        xy=(
            rmajor - rminor * triang,
            kappa * rminor * 0.25,
        ),  # Inner plasma edge
        xytext=(rmajor, kappa * rminor * 0.25),  # Center
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )

    # Write the triangularity to the left of the cross, position relative to figure axes
    draw_text(
        axis,
        rmajor - (rminor * triang * 0.75),
        kappa * rminor * 0.3,
        f"$\\delta$: {mfile.get('triang', scan=scan):.2f}",
        fontsize=9,
        color="black",
        rotation=0,
        verticalalignment="center",
        bbox=white_box,
    )

    # =============================================

    radius_plasma_core_norm = mfile.get("radius_plasma_core_norm", scan=scan)

    # Draw a double-ended arrow for the plasma core region
    draw_annotation(
        axis,
        "",
        xy=(rmajor, -rminor * 0.1 * kappa),  # Inner plasma edge
        xytext=(
            rmajor + (rminor * radius_plasma_core_norm),
            -rminor * 0.1 * kappa,
        ),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    # Add a label for core region
    draw_text(
        axis,
        rmajor + (rminor * radius_plasma_core_norm / 4),
        -rminor * kappa * 0.15,
        f"$\\rho_{{\\text{{core}}}}$: {radius_plasma_core_norm:.2f}",
        fontsize=9,
        color="black",
        rotation=0,
        verticalalignment="center",
        bbox=white_box,
    )

    # ================================================

    # Add plasma volume, areas and shaping information

    geom_type = PlasmaGeometryModelType(mfile.get("i_plasma_geometry", scan=scan))

    textstr_plasma = (
        "$\\mathbf{Shaping:}$\n\n$\\kappa_{95}$:"
        f" {mfile.get('kappa95', scan=scan):.2f}"
        f" ({geom_type.kappa95_model.description}) | $\\delta_{{95}}$:"
        f" {mfile.get('triang95', scan=scan):.2f}"
        f" ({geom_type.triang95_model.description}) | $\\zeta$:"
        f" {mfile.get('plasma_square', scan=scan):.2f}\n$\\kappa$:"
        f" {mfile.get('kappa', scan=scan):.2f}"
        f" ({geom_type.kappa_model.description}) | $\\delta$:"
        f" {mfile.get('triang', scan=scan):.2f}"
        f" ({geom_type.triang_model.description}) | A:"
        f" {mfile.get('aspect', scan=scan):.2f}\n$ V_{{\\text{{p}}}}:$"
        f" {mfile.get('vol_plasma', scan=scan):,.2f}$ \\ \\text{{m}}^3$ | $"
        " A_{\\text{p,surface}}:$"
        f" {mfile.get('a_plasma_surface', scan=scan):,.2f}$ \\ \\text{{m}}^2$"
        " | $ A_{\\text{p,poloidal}}:$"
        f" {mfile.get('a_plasma_poloidal', scan=scan):,.3f}$ \\"
        " \\text{m}^2$\n$ L_{\\text{p,poloidal}}:$"
        f" {mfile.get('len_plasma_poloidal', scan=scan):,.3f}$ \\ \\text{{m}}$"
    )

    draw_text(
        axis,
        0.365,
        0.975,
        textstr_plasma,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=box_style("lightyellow"),
    )

    # ============================================

    # Draw a red arrow coming from the right and pointing at the plasma
    for kap in (-kappa, kappa):
        draw_annotation(
            axis,
            "",
            # Pointing at plasma
            xy=(rmajor + (rminor * 0.8), kap * rminor * 0.2),
            # Starting point of arrow
            xytext=(rmajor + (rminor * 1.4), kap * rminor * 0.2),
            arrowprops={"facecolor": "red", "edgecolor": "red", "lw": 2},
        )

    i_hcd_primary = mfile.get("i_hcd_primary", scan=scan)
    i_hcd_secondary = mfile.get("i_hcd_secondary", scan=scan)

    # Add heating and current drive information
    textstr_hcd = (
        "$\\mathbf{Heating \\ & \\ current \\ drive:}$\n\nTotal injected"
        f" heat: {mfile.get('p_hcd_injected_total_mw', scan=scan):.3f}"
        " MW\nOhmic heating power:"
        f" {mfile.get('p_plasma_ohmic_mw', scan=scan):.3f}"
        " MW\n\n$\\mathbf{Primary \\ system:"
        f" {CurrentDriveModel(i_hcd_primary).abbreviation}}}$\nCurrent driving"
        f" power {mfile.get('p_hcd_primary_injected_mw', scan=scan):.4f}"
        " MW\nExtra heat power:"
        f" {mfile.get('p_hcd_primary_extra_heat_mw', scan=scan):.4f}"
        " MW\n$\\eta_{\\text{CD,prim}}$:"
        f" {mfile.get('eta_cd_hcd_primary', scan=scan):.4f} A/W  |  "
        " $\\langle\\zeta_{\\text{CD,prim}}\\rangle$:"
        f" {mfile.get('eta_cd_dimensionless_hcd_primary', scan=scan):.4f}\n$\\gamma_{{\\text{{CD,prim}}}}$:"  # noqa: E501
        f" {mfile.get('eta_cd_norm_hcd_primary', scan=scan):.4f} $\\times"
        " 10^{20}  \\mathrm{A} / \\mathrm{Wm}^2$\nCurrent driven by"
        f" primary: {mfile.get('c_hcd_primary_driven', scan=scan) / 1e6:.3f}"
        " MA\n\n$\\mathbf{Secondary \\ system:"
        f" {CurrentDriveModel(i_hcd_secondary).abbreviation}}}$\nCurrent"
        " driving power"
        f" {mfile.get('p_hcd_secondary_injected_mw', scan=scan):.4f} MW\nExtra"
        " heat power:"
        f" {mfile.get('p_hcd_secondary_extra_heat_mw', scan=scan):.4f}"
        " MW\n$\\eta_{\\text{CD,sec}}$:"
        f" {mfile.get('eta_cd_hcd_secondary', scan=scan):.4f} A/W  |  "
        " $\\langle\\zeta_{\\text{CD,sec}}\\rangle$:"
        f" {mfile.get('eta_cd_dimensionless_hcd_secondary', scan=scan):.4f}\n$\\gamma_{{\\text{{CD,sec}}}}$:"  # noqa: E501
        f" {mfile.get('eta_cd_norm_hcd_secondary', scan=scan):.4f} $\\times"
        " 10^{20}  \\mathrm{A} / \\mathrm{Wm}^2$\nCurrent driven by"
        " secondary:"
        f" {mfile.get('c_hcd_secondary_driven', scan=scan) / 1e6:.3f} MA"
    )

    draw_text(
        axis,
        0.73,
        0.675,
        textstr_hcd,
        fontsize=9,
        verticalalignment="top",
        transform=plt.gcf().transFigure,
        bbox=box_style("paleturquoise") | {"edgecolor": "black"},
    )

    class TextArgs(TypedDict):
        fontsize: int
        verticalalignment: str
        transform: Transform

    text_args = TextArgs({
        "fontsize": 23,
        "verticalalignment": "top",
        "transform": fig.transFigure,
    })

    # Add injected power label
    draw_text(axis, 0.92, 0.625, "$P_{\\text{inj}}$", **text_args)

    # ================================================

    # Add beta information
    textstr_beta = (
        "$\\mathbf{Beta \\ Information:}$\n\nTotal beta,$ \\ \\langle"
        " \\beta \\rangle$:"
        f" {mfile.get('beta_total_vol_avg', scan=scan):.4f}\nThermal beta,$ \\"
        " \\langle \\beta_{\\text{thermal}} \\rangle$:"
        f" {mfile.get('beta_thermal_vol_avg', scan=scan):.4f}\nToroidal beta,$"
        " \\ \\langle \\beta_{\\text{t}} \\rangle$:"
        f" {mfile.get('beta_toroidal_vol_avg', scan=scan):.4f}\nPoloidal"
        " beta,$ \\ \\langle \\beta_{\\text{p}} \\rangle$:"
        f" {mfile.get('beta_poloidal_vol_avg', scan=scan):.4f}\nFast-alpha"
        " beta,$ \\ \\langle \\beta_{\\alpha} \\rangle$:"
        f" {mfile.get('beta_fast_alpha', scan=scan):.4f}\nUpper limit on"
        f" {BetaComponentLimits(int(mfile.get('i_beta_component', scan=scan))).full_name}:"  # noqa: E501
        " $ \\langle \\beta \\rangle$:"
        f" {mfile.get('beta_vol_avg_max', scan=scan):.4f}\nNormalised total"
        " beta,$ \\ \\beta_{\\text{N}}$:"
        f" {mfile.get('beta_norm_total', scan=scan):.4f}\nNormalised thermal"
        " beta,$ \\ \\beta_{\\text{N,thermal}}$:"
        f" {mfile.get('beta_norm_thermal', scan=scan):.4f}\nMaximum normalised"
        " beta"
        f" ({BetaNormMaxModel(int(mfile.get('i_beta_norm_max', scan=scan))).full_name}),$"  # noqa: E501
        " \\ \\beta_{\\text{N,max}}$:"
        f" {mfile.get('beta_norm_max', scan=scan):.4f}"
    )

    draw_text(
        axis,
        0.025,
        0.975,
        textstr_beta,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("lightblue"),
    )

    # Add beta label
    draw_text(axis, 0.27, 0.94, "$\\beta$", **text_args)

    # ================================================

    # Add volt-second information
    textstr_volt_second = (
        "$\\mathbf{Volt-second \\ requirements:}$\n\nTotal volt-second"
        f" consumption: {mfile.get('vs_plasma_total_required', scan=scan):.4f}"
        " Vs\n  - Internal volt-seconds:"
        f" {mfile.get('vs_plasma_internal', scan=scan):.4f} Vs\n  -"
        " Volt-seconds needed for burn:"
        f" {mfile.get('vs_plasma_burn_required', scan=scan):.4f} Vs\n  -"
        " Volt-seconds needed for ramp:"
        f" {mfile.get('vs_plasma_ramp_required', scan=scan):.4f} Vs |"
        " $C_{\\text{ejima}}$:"
        f" {mfile.get('ejima_coeff', scan=scan):.4f}\n$V_{{\\text{{loop}}}}$:"
        f" {mfile.get('v_plasma_loop_burn', scan=scan):.4f}"
        " V\n$\\Omega_{\\text{p}}$:"
        f" {mfile.get('res_plasma', scan=scan):.4e} $\\Omega$\nPlasma"
        " resistive diffusion time:"
        f" {mfile.get('t_plasma_res_diffusion', scan=scan):,.4f} s\nPlasma"
        f" inductance: {mfile.get('ind_plasma', scan=scan):.4e} H | ITER"
        " $l_i(3)$:"
        f" {mfile.get('ind_plasma_internal_norm_iter_3', scan=scan):.4f}\nPlasma"
        " stored magnetic energy:"
        f" {mfile.get('e_plasma_magnetic_stored', scan=scan) / 1e9:.4f}"
        " GJ\nPlasma normalised internal inductance, $l_i$"
        f" ({IndInternalNormModel(int(mfile.get('i_ind_plasma_internal_norm', scan=scan))).full_name})"  # noqa: E501
        f" :{mfile.get('ind_plasma_internal_norm', scan=scan):.3f}"
    )

    draw_text(
        axis,
        0.025,
        0.78,
        textstr_volt_second,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("lightgreen"),
    )

    # Add volt second label
    draw_text(axis, 0.30, 0.77, "Vs", **text_args)

    # =========================================

    # Add divertor information
    textstr_div = (
        "\n$P_{\\text{sep}}$:"
        f" {mfile.get('p_plasma_separatrix_mw', scan=scan):.2f}"
        " MW\n$\\frac{P_{\\text{sep}}}{R}$:"
        f" {mfile.get('p_plasma_separatrix_rmajor_mw', scan=scan):.2f}"
        " MW/m\n$\\frac{P_{\\text{sep}}B_T}{q_{95} A  R}$:"
        f" {mfile.get('p_div_bt_q_aspect_rmajor_mw', scan=scan):.2f} MW T/m   "
        "            "
    )

    draw_text(
        axis,
        0.35,
        0.12,
        textstr_div,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("orange"),
    )

    # Add divertor label
    draw_text(axis, 0.45, 0.1, "$P_{\\text{div}}$", **text_args)

    # ================================================

    # Add confinement information
    textstr_confinement = (
        "$\\mathbf{Confinement:}$\n\nConfinement scaling law:"
        f" {mfile.get('tauelaw', scan=scan)}\nConfinement $H$ factor:"
        f" {mfile.get('hfact', scan=scan):.4f}\nEnergy confinement time from"
        f" scaling: {mfile.get('t_energy_confinement', scan=scan):.4f}"
        f" s\nFusion double product: {mfile.get('ntau', scan=scan):.4e}"
        f" s/m³\nLawson Triple product: {mfile.get('nttau', scan=scan):.4e}"
        " keV·s/m³\nTransport loss power assumed in scaling law:"
        f" {mfile.get('p_plasma_loss_mw', scan=scan):.4f} MW\nPlasma thermal"
        " energy (inc. $\\alpha$), $W$:"
        f" {mfile.get('e_plasma_beta', scan=scan) / 1e9:.4f} GJ\nAlpha"
        " particle confinement time:"
        f" {mfile.get('t_alpha_confinement', scan=scan):.4f} s |"
        " $\\tau_{\\alpha}/\\tau_{e}$:"
        f" {mfile.get('f_t_alpha_energy_confinement', scan=scan):.4f}"
    )

    draw_text(
        axis,
        0.025,
        0.57,
        textstr_confinement,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        # Changed to a not normal color (Aquamarine)
        bbox=box_style("gainsboro") | {"edgecolor": "black"},
    )

    # Add tau label
    draw_text(axis, 0.3, 0.55, "$\\tau_{\\text{e}} $", **text_args)

    # =========================================

    # Load the neutron image
    alpha_particle = load_plot_image("alpha_particle.png")

    # Display the neutron image over the figure, not the axes
    new_ax = axis.inset_axes(
        (0.975, 0.275, 0.075, 0.075), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(alpha_particle)
    new_ax.axis("off")

    draw_annotation(
        axis,
        "",
        xy=(rmajor + rminor, -rminor * kappa * 0.55),  # Pointing at the plasma
        xytext=(rmajor + 0.2 * rminor, -rminor * kappa * 0.25),
        arrowprops={"facecolor": "red", "edgecolor": "grey", "lw": 1},
    )

    textstr_alpha = (
        "$P_{\\alpha,\\text{loss}}$"
        f" {mfile.get('p_fw_alpha_mw', scan=scan):.2f}"
        " MW\n$f_{\\alpha,\\text{coupled}}$"
        f" {mfile.get('f_p_alpha_plasma_deposited', scan=scan):.2f}"
    )

    draw_text(
        axis,
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
    neutron = load_plot_image("neutron.png")
    new_ax = axis.inset_axes(
        (0.975, 0.75, 0.075, 0.075), transform=axis.transAxes, zorder=10
    )
    new_ax.imshow(neutron)
    new_ax.axis("off")

    # Draw a red arrow coming from the right and pointing at the plasma
    draw_annotation(
        axis,
        "",
        xy=(rmajor + rminor, rminor * kappa * 0.65),  # Pointing at the plasma
        xytext=(rmajor, rminor * kappa * 0.5),
        arrowprops={"facecolor": "grey", "edgecolor": "grey", "lw": 1},
    )

    textstr_neutron = (
        "$P_{\\text{n,total}}$"
        f" {mfile.get('p_neutron_total_mw', scan=scan):.2f}"
        " MW\n$\\phi_{\\text{n,avg}}$"
        f" {mfile.get('pflux_plasma_surface_neutron_avg_mw', scan=scan):.3f}"
        " MW/m²"
    )

    draw_text(
        axis,
        0.775,
        0.875,
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
        "$\\mathbf{Fusion \\ Reactions:}$\n\nFuel mixture:\n|  D:"
        f" {mfile.get('f_plasma_fuel_deuterium', scan=scan):.2f}  |  T:"
        f" {mfile.get('f_plasma_fuel_tritium', scan=scan):.2f}  |  3He:"
        f" {mfile.get('f_plasma_fuel_helium3', scan=scan):.2f}  |\n\nFusion"
        " Power, $P_{\\text{fus}}:$"
        f" {mfile.get('p_fusion_total_mw', scan=scan):,.2f} MW\nD-T Power,"
        " $P_{\\text{fus,DT}}:$"
        f" {mfile.get('p_dt_total_mw', scan=scan):,.2f} MW\nD-D Power,"
        " $P_{\\text{fus,DD}}:$"
        f" {mfile.get('p_dd_total_mw', scan=scan):,.2f} MW\nD-3He Power,"
        " $P_{\\text{fus,D3He}}:$"
        f" {mfile.get('p_dhe3_total_mw', scan=scan):,.2f} MW\nAlpha Power,"
        f" $P_{{\\alpha}}:$ {mfile.get('p_alpha_total_mw', scan=scan):,.2f} MW"
    )

    draw_text(
        axis,
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
        "$\\mathbf{Fuelling:}$\n\nPlasma mass:"
        f" {mfile.get('m_plasma', scan=scan) * 1000:.4f} g\n   - Average mass"
        f" of all plasma ions: {mfile.get('m_ions_total_amu', scan=scan):.3f}"
        " amu\nFuel mass:"
        f" {mfile.get('m_plasma_fuel_ions', scan=scan) * 1000:.4f} g\n   -"
        " Average mass of all fuel ions:"
        f" {mfile.get('m_fuel_amu', scan=scan):.3f} amu\n\nFueling rate:"
        f" {mfile.get('molflow_plasma_fuelling_required', scan=scan):.3e}"
        " nucleus-pairs/s\nFuel burn-up rate:"
        f" {mfile.get('rndfuel', scan=scan):.3e} reactions/s\nBurn-up"
        f" fraction: {mfile.get('burnup', scan=scan):.4f}\n"
    )

    draw_text(
        axis,
        0.025,
        0.22,
        textstr_fuelling,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("khaki") | {"edgecolor": "black"},
    )

    # ================================================

    # Add ion density information
    impurity_data = ImpurityRadiationData()
    textstr_ions = (
        f"             $\\mathbf{{Ion \\ to \\ electron}}$\n"
        f"             $\\mathbf{{relative \\ number}}$\n"
        f"             $\\mathbf{{densities:}}$\n\n"
        "             Effective charge: "
        f"{mfile.get('n_charge_plasma_effective_vol_avg', scan=scan):.3f}\n\n"
        + "\n".join(
            f"             {label.replace('_', '') + ':':<6}"
            f"{mfile.get(f'f_nd_impurity_electrons({index:02d})', scan=scan):.4e}"
            for index, label in enumerate(impurity_data.imp_label[:14], start=1)
        )
    )

    draw_text(
        axis,
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
    draw_text(axis, 0.815, 0.29, "$Z$", **text_args)

    # ================================================

    # Add plasma current information
    textstr_currents = (
        "$\\mathbf{Plasma\\ currents:}$\n\nPlasma current"
        f" ({PlasmaCurrentModel(int(mfile.get('i_plasma_current', scan=scan))).full_name}):"  # noqa: E501
        f" {mfile.get('plasma_current_ma', scan=scan):.4f} MA\n  - Bootstrap"
        " fraction"
        f" ({BootstrapCurrentFractionModel(int(mfile.get('i_bootstrap_current', scan=scan))).full_name}):"  # noqa: E501
        f" {mfile.get('f_c_plasma_bootstrap', scan=scan):.4f}\n  - Diamagnetic"
        " fraction"
        f" ({PlasmaDiamagneticCurrentModel(int(mfile.get('i_diamagnetic_current', scan=scan))).full_name}):"  # noqa: E501
        f" {mfile.get('f_c_plasma_diamagnetic', scan=scan):.4f}\n  -"
        " Pfirsch-Schlüter fraction"
        f" {mfile.get('f_c_plasma_pfirsch_schluter', scan=scan):.4f}\n  -"
        " Auxiliary fraction"
        f" {mfile.get('f_c_plasma_auxiliary', scan=scan):.4f}\n  - Inductive"
        f" fraction {mfile.get('f_c_plasma_inductive', scan=scan):.4f}"
    )

    draw_text(
        axis,
        0.72,
        0.975,
        textstr_currents,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("#C8A2C8"),  # Hex code for lilac color
    )

    # Add plasma current label
    draw_text(axis, 0.93, 0.9, "$I_{\\text{p}} $", **text_args)

    # Add magnetic field information
    textstr_fields = (
        "$\\mathbf{Magnetic\\ fields:}$\n\nToroidal field at $R_0$,"
        f" $B_{{T}}$: {mfile.get('b_plasma_toroidal_on_axis', scan=scan):.4f}"
        " T\n  Ripple at outboard , $\\delta$:"
        f" {mfile.get('ripple_b_tf_plasma_edge', scan=scan):.2f}%\nSurface"
        " average poloidal field, $\\langle B_{p}(a) \\rangle$:"
        f" {mfile.get('b_plasma_surface_poloidal_average', scan=scan):.4f}"
        " T\nTotal field, $B_{tot}$:"
        f" {mfile.get('b_plasma_total', scan=scan):.4f} T\nVertical field,"
        " $B_{vert}$:"
        f" {mfile.get('b_plasma_vertical_required', scan=scan):.4f} T"
    )

    draw_text(
        axis,
        0.5325,
        0.14,
        textstr_fields,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("royalblue"),
    )

    # Add magnetic field label
    draw_text(axis, 0.75, 0.12, "$B$", **text_args)

    # Add radiation information
    textstr_radiation = (
        "           $\\mathbf{Radiation:}$\n\n           Total radiation"
        f" power {mfile.get('p_plasma_rad_mw', scan=scan):.4f} MW\n          "
        " Separatrix radiation fraction"
        f" {mfile.get('f_p_plasma_separatrix_rad', scan=scan):.4f}\n          "
        " Core radiation power"
        f" {mfile.get('p_plasma_inner_rad_mw', scan=scan):.4f} MW\n           "
        "   - $f_{\\text{core,reduce}}$"
        f" {mfile.get('f_p_plasma_core_rad_reduction', scan=scan):.4f}\n      "
        "     Edge radiation power"
        f" {mfile.get('p_plasma_outer_rad_mw', scan=scan):.4f} MW\n          "
        " Synchrotron radiation power"
        f" {mfile.get('p_plasma_sync_mw', scan=scan):.4f} MW\n          "
        " Synchrotron wall reflectivity"
        f" {mfile.get('f_sync_reflect', scan=scan):.4f}"
    )

    draw_text(
        axis,
        0.72,
        0.83,
        textstr_radiation,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("lavender") | {"edgecolor": "black"},
    )

    # Add radiation label
    draw_text(axis, 0.725, 0.78, "$\\gamma$", **text_args)

    # Add L-H threshold information
    model_name = PlasmaConfinementTransitionModel(
        int(mfile.get("i_l_h_threshold", scan=scan))
    ).full_name

    # Wrap long model names to new line
    if len(model_name) > 20:
        model_name = "\n".join(textwrap.wrap(model_name, width=20))

    textstr_lh = (
        "$\\mathbf{L-H \\"
        f" threshold:}}$\n{model_name}\n\n$P_{{\\text{{L-H}}}}:$"
        f" {mfile.get('p_l_h_threshold_mw', scan=scan):.4f} MW"
    )

    draw_text(
        axis,
        0.22,
        0.4,
        textstr_lh,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("peachpuff"),
    )

    # Add density limit information
    textstr_density_limit = (
        "$\\mathbf{Density \\"
        f" limit:}}$\n({DensityLimitModel(int(mfile.get('i_density_limit', scan=scan))).full_name})\n$n_{{\\text{{e,limit}}}}:"  # noqa: E501
        f" {mfile.get('nd_plasma_electrons_max', scan=scan):.3e} \\"
        " m^{-3}$\n$f_{\\text{GW}}$:"
        f" {mfile.get('f_nd_plasma_greenwald', scan=scan):.4f}"
    )

    draw_text(
        axis,
        0.22,
        0.31,
        textstr_density_limit,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("pink"),
    )


def plot_detailed_plasma_parameters(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Function to plot detailed plasma parameters from physics data.

    Parameters
    ----------
    axis : plt.Axes
        Axis object to plot to
    fig : plt.Figure
        Figure object for text placement
    mfile : MFile
        MFILE data object
    scan : int
        Scan number to use
    """
    textstr_debye = (
        "$\\mathbf{Debye \\"
        " Lengths:}$\n\n$\\langle\\lambda_{Debye,e}\\rangle$:"
        f" {mfile.get('len_plasma_debye_electron_vol_avg', scan=scan):.4e} m"
    )

    textstr_larmor = (
        "$\\mathbf{Larmor \\ Radii:}$\n\n"
        "$\\langle\\rho_{Larmor,toroidal,D}\\rangle$:"
        f" {mfile.get('radius_plasma_deuteron_toroidal_larmor_isotropic_vol_avg', scan=scan):.4e} m\n"  # noqa: E501
        "$\\langle\\rho_{Larmor,toroidal,T}\\rangle$:"
        f" {mfile.get('radius_plasma_triton_toroidal_larmor_isotropic_vol_avg', scan=scan):.4e} m"  # noqa: E501
    )

    textstr_velocities = (
        "$\\mathbf{Velocities:}$\n\n$\\langle v_{e}\\rangle$:"
        f" {mfile.get('vel_plasma_electron_vol_avg', scan=scan):.4e}"
        " m/s\n$\\langle v_{D}\\rangle$:"
        f" {mfile.get('vel_plasma_deuteron_vol_avg', scan=scan):.4e}"
        " m/s\n$\\langle v_{T}\\rangle$:"
        f" {mfile.get('vel_plasma_triton_vol_avg', scan=scan):.4e}"
        " m/s\n$\\langle v_{\\alpha,thermal}\\rangle$:"
        f" {mfile.get('vel_plasma_alpha_thermal_vol_avg', scan=scan):.4e}"
        " m/s\n$v_{\\alpha,birth}$:"
        f" {mfile.get('vel_plasma_alpha_birth', scan=scan):.4e} m/s"
    )

    textstr_frequencies = (
        "$\\mathbf{Frequencies:}$\n\n$\\langle\\omega_{p,e}\\rangle$:"
        f" {mfile.get('freq_plasma_electron_vol_avg', scan=scan):.4e}"
        " Hz\n$\\langle f_{Larmor,toroidal,e}\\rangle$:"
        f" {mfile.get('freq_plasma_larmor_toroidal_electron_vol_avg', scan=scan):.4e}"
        " Hz\n$\\langle f_{Larmor,toroidal,D}\\rangle$:"
        f" {mfile.get('freq_plasma_larmor_toroidal_deuteron_vol_avg', scan=scan):.4e}"
        " Hz\n$\\langle f_{Larmor,toroidal,T}\\rangle$:"
        f" {mfile.get('freq_plasma_larmor_toroidal_triton_vol_avg', scan=scan):.4e}"
        " Hz\n$\\langle\\omega_{UH,e}\\rangle$:"
        f" {mfile.get('freq_plasma_upper_hybrid_vol_avg', scan=scan):.4e} Hz"
    )

    textstr_coulomb = (
        "$\\mathbf{Coulomb \\ Logarithms:}$\n\n"
        "$\\langle\\ln \\Lambda_{e-e}\\rangle$:"
        f" {mfile.get('plasma_coulomb_log_electron_electron_vol_avg', scan=scan):.4f}\n"
        "$\\langle\\ln \\Lambda_{e-D}\\rangle$:"
        f" {mfile.get('plasma_coulomb_log_electron_deuteron_vol_avg', scan=scan):.4f}\n"
        "$\\langle\\ln \\Lambda_{e-T}\\rangle$:"
        f" {mfile.get('plasma_coulomb_log_electron_triton_vol_avg', scan=scan):.4f}\n"
        "$\\langle\\ln \\Lambda_{D-T}\\rangle$:"
        f" {mfile.get('plasma_coulomb_log_deuteron_triton_vol_avg', scan=scan):.4f}\n"
        "$\\langle\\ln \\Lambda_{e-\\alpha}\\rangle$:"
        f" {mfile.get('plasma_coulomb_log_electron_alpha_thermal_vol_avg', scan=scan):.4f}"  # noqa: E501
    )

    textstr_collision_times = (
        "$\\mathbf{Collision \\ Times:}$\n\n"
        "$\\langle\\tau_{e-e}\\rangle$:"
        f" {mfile.get('t_plasma_electron_electron_collision_vol_avg', scan=scan):.4e} s\n"  # noqa: E501
        "$\\langle\\tau_{e-D}\\rangle$:"
        f" {mfile.get('t_plasma_electron_deuteron_collision_vol_avg', scan=scan):.4e} s\n"  # noqa: E501
        "$\\langle\\tau_{e-T}\\rangle$:"
        f" {mfile.get('t_plasma_electron_triton_collision_vol_avg', scan=scan):.4e} s\n"
        "$\\langle\\tau_{e-\\alpha}\\rangle$:"
        f" {mfile.get('t_plasma_electron_alpha_thermal_collision_vol_avg', scan=scan):.4e} s"  # noqa: E501
    )

    textstr_collision_freq = (
        "$\\mathbf{Collision \\ Frequencies:}$\n\n"
        "$\\langle\\nu_{e-e}\\rangle$:"
        f" {mfile.get('freq_plasma_electron_electron_collision_vol_avg', scan=scan):.4e}"
        " Hz\n"
        "$\\langle\\nu_{e-D}\\rangle$:"
        f" {mfile.get('freq_plasma_electron_deuteron_collision_vol_avg', scan=scan):.4e}"
        " Hz\n"
        "$\\langle\\nu_{e-T}\\rangle$:"
        f" {mfile.get('freq_plasma_electron_triton_collision_vol_avg', scan=scan):.4e}"
        " Hz\n"
        "$\\langle\\nu_{e-\\alpha}\\rangle$:"
        f" {mfile.get('freq_plasma_electron_alpha_thermal_collision_vol_avg', scan=scan):.4e} Hz"  # noqa: E501
    )

    textstr_mfp = (
        "$\\mathbf{Mean \\ Free \\ Paths:}$\n\n"
        "$\\langle\\lambda_{mfp,e-e}\\rangle$:"
        f" {mfile.get('len_plasma_electron_electron_mean_free_path_vol_avg', scan=scan):.4e} m\n"  # noqa: E501
        "$\\langle\\lambda_{mfp,e-D}\\rangle$:"
        f" {mfile.get('len_plasma_electron_deuteron_mean_free_path_vol_avg', scan=scan):.4e} m\n"  # noqa: E501
        "$\\langle\\lambda_{mfp,e-T}\\rangle$:"
        f" {mfile.get('len_plasma_electron_triton_mean_free_path_vol_avg', scan=scan):.4e} m\n"  # noqa: E501
        "$\\langle\\lambda_{mfp,e-\\alpha}\\rangle$:"
        f" {mfile.get('len_plasma_electron_alpha_thermal_mean_free_path_vol_avg', scan=scan):.4e} m"  # noqa: E501
    )

    textstr_spitzer = (
        "$\\mathbf{Spitzer \\ Slowing \\"
        " Down:}$\n\n$\\langle\\tau_{e-\\alpha,Spitzer}\\rangle$:"
        f" {mfile.get('t_plasma_electron_alpha_spitzer_slow_vol_avg', scan=scan):.4e} s"
    )

    textstr_resistivity = (
        "$\\mathbf{Resistivities:}$\n\n$\\langle\\eta_{Spitzer}\\rangle$:"
        f" {mfile.get('res_plasma_fuel_spitzer_vol_avg', scan=scan):.4e}"
        " $\\Omega\\mathrm{m}$"
    )

    light_yellow_box = {
        "boxstyle": "round",
        "facecolor": "lightyellow",
        "alpha": 1.0,
        "linewidth": 2,
    }

    draw_text(
        axis,
        0.05,
        0.45,
        textstr_debye,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_yellow_box,
    )

    draw_text(
        axis,
        0.25,
        0.45,
        textstr_larmor,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_yellow_box,
    )

    draw_text(
        axis,
        0.45,
        0.45,
        textstr_velocities,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_yellow_box,
    )

    light_cyan_box = {
        "boxstyle": "round",
        "facecolor": "lightcyan",
        "alpha": 1.0,
        "linewidth": 2,
    }

    draw_text(
        axis,
        0.05,
        0.31,
        textstr_frequencies,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_cyan_box,
    )

    draw_text(
        axis,
        0.25,
        0.31,
        textstr_coulomb,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_cyan_box,
    )

    draw_text(
        axis,
        0.45,
        0.31,
        textstr_collision_times,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_cyan_box,
    )

    light_green_box = {
        "boxstyle": "round",
        "facecolor": "lightgreen",
        "alpha": 1.0,
        "linewidth": 2,
    }

    draw_text(
        axis,
        0.05,
        0.17,
        textstr_collision_freq,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_green_box,
    )

    draw_text(
        axis,
        0.25,
        0.17,
        textstr_mfp,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_green_box,
    )

    draw_text(
        axis,
        0.45,
        0.17,
        textstr_spitzer + "\n" + textstr_resistivity,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        transform=fig.transFigure,
        bbox=light_green_box,
    )

    axis.axis("off")


__all__ = ["plot_detailed_plasma_parameters", "plot_main_plasma_information"]
