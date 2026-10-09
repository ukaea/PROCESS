"""Plasma functions for PROCESS summary plots."""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING, Literal

from matplotlib.patches import FancyBboxPatch

from process.core.io.plot.summary.common import box_style, load_plot_image, text_layout
from process.core.io.plot.summary.plasma.physics import plot_plasma
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.geometry.plasma import plasma_geometry
from process.models.physics.bootstrap_current import BootstrapCurrentFractionModel
from process.models.physics.current_drive import CurrentDriveModel
from process.models.physics.density_limit import DensityLimitModel
from process.models.physics.l_h_transition import PlasmaConfinementTransitionModel
from process.models.physics.physics import (
    BetaComponentLimits,
    BetaNormMaxModel,
    IndInternalNormModel,
)
from process.models.physics.plasma_current import (
    PlasmaCurrentModel,
    PlasmaDiamagneticCurrentModel,
)
from process.models.physics.plasma_geometry import PlasmaGeometryModelType

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plasma_main_page(
    axis,
    fig,
    rmajor,
    rminor,
    triang,
    radius_plasma_core_norm,
    kappa,
    i_single_null,
    plasma_square,
    big_q_plasma,
    p_fw_alpha_mw,
    f_p_alpha_plasma_deposited,
    p_neutron_total_mw,
    pflux_plasma_surface_neutron_avg_mw,
):
    """Plot plasma from main page"""
    white_box = box_style("white", linewidth=None)
    pg = plasma_geometry(
        rmajor=rmajor,
        rminor=rminor * radius_plasma_core_norm,
        triang=triang,
        kappa=kappa,
        i_single_null=i_single_null,
        i_plasma_shape=1,
        square=plasma_square,
    )
    axis.plot(pg.rs, pg.zs, color="black", linestyle="--")
    axis.plot(rmajor, 0, "r+", markersize=20, markeredgewidth=2)

    axis.text(
        0.725,
        0.175,
        f"$Q_{{\\text{{plasma}}}}$: {big_q_plasma:.2f}",
        fontsize=15,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=white_box,
        transform=fig.transFigure,
    )

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
        bbox=white_box,
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
        bbox=white_box,
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
        f"$\\kappa$: {kappa:.2f}",
        fontsize=9,
        color="black",
        rotation=270,
        verticalalignment="center",
        transform=axis.transAxes,
        bbox=white_box,
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
        f"$\\delta$: {triang:.2f}",
        fontsize=9,
        color="black",
        verticalalignment="center",
        bbox=white_box,
    )

    # Draw a double-ended arrow for the plasma core region
    axis.annotate(
        "",
        xy=(rmajor, -rminor * 0.1 * kappa),
        xytext=(rmajor + rminor * radius_plasma_core_norm, -rminor * 0.1 * kappa),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    # Add a label for core region
    axis.text(
        rmajor + (rminor * radius_plasma_core_norm / 4),
        -rminor * kappa * 0.15,
        f"$\\rho_{{\\text{{core}}}}$: {radius_plasma_core_norm:.2f}",
        fontsize=9,
        color="black",
        verticalalignment="center",
        bbox=white_box,
    )

    alpha_particle = load_plot_image("alpha_particle.png")
    image_axis = axis.inset_axes(
        (0.975, 0.275, 0.075, 0.075),
        transform=axis.transAxes,
        zorder=10,
    )
    image_axis.imshow(alpha_particle)
    image_axis.axis("off")

    axis.annotate(
        "",
        xy=(rmajor + rminor, -rminor * kappa * 0.55),
        xytext=(rmajor + 0.2 * rminor, -rminor * kappa * 0.25),
        arrowprops={"facecolor": "red", "edgecolor": "grey", "lw": 1},
    )
    axis.text(
        1.0,
        0.275,
        f"$P_{{\\alpha,\\text{{loss}}}}$ {p_fw_alpha_mw:.2f} MW\n"
        f"$f_{{\\alpha,\\text{{coupled}}}}$ "
        f"{f_p_alpha_plasma_deposited:.2f}",
        fontsize=9,
        color="black",
        verticalalignment="top",
        transform=axis.transAxes,
        bbox=box_style("red"),
    )

    neutron = load_plot_image("neutron.png")
    image_axis = axis.inset_axes(
        (0.975, 0.75, 0.075, 0.075),
        transform=axis.transAxes,
        zorder=10,
    )
    image_axis.imshow(neutron)
    image_axis.axis("off")

    axis.annotate(
        "",
        xy=(rmajor + rminor, rminor * kappa * 0.65),
        xytext=(rmajor, rminor * kappa * 0.5),
        arrowprops={"facecolor": "grey", "edgecolor": "grey", "lw": 1},
    )
    axis.text(
        0.775,
        0.875,
        f"$P_{{\\text{{n,total}}}}$ {p_neutron_total_mw:.2f} MW\n"
        f"$\\phi_{{\\text{{n,avg}}}}$ "
        f"{pflux_plasma_surface_neutron_avg_mw:.3f} MW/m²",
        fontsize=9,
        color="black",
        verticalalignment="top",
        transform=axis.transAxes,
        bbox=box_style("grey", alpha=0.8),
    )
    for kap in (-kappa, kappa):
        axis.annotate(
            "",
            xy=(rmajor + rminor * 0.8, kap * rminor * 0.2),
            xytext=(rmajor + rminor * 1.4, kap * rminor * 0.2),
            arrowprops={"facecolor": "red", "edgecolor": "red", "lw": 2},
        )


def plot_main_plasma_information(
    axis: plt.Axes,
    mfile: MFile,
    scan: int,
    colour_scheme: Literal[1, 2],
    fig: plt.Figure,
):
    """Plots the main plasma information including plasma shape, geometry, currents,
    heating, confinement, and other relevant plasma parameters.

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

    def value(name: str):
        return mfile.get(name, scan=scan)

    def add_panel(
        x: float,
        top: float,
        width: float,
        height: float,
        text: str,
        style: dict,
        text_x_offset: float = 0.004,
    ):
        """Add one fixed-size box and its unchanged original text."""
        box_expansion = 0.002
        left = x - box_expansion
        bottom = top - height - box_expansion
        expanded_width = width + 2 * box_expansion
        expanded_height = height + 2 * box_expansion

        patch = FancyBboxPatch(
            (left, bottom),
            expanded_width,
            expanded_height,
            boxstyle="round,pad=0,rounding_size=0.003",
            transform=fig.transFigure,
            facecolor=style.get("facecolor", "white"),
            edgecolor=style.get("edgecolor", "black"),
            linewidth=style.get("linewidth", 1.0),
            alpha=style.get("alpha", 1.0),
            clip_on=False,
            zorder=20,
        )
        axis.add_patch(patch)

        return axis.text(
            x + text_x_offset,
            top - 0.004,
            text,
            fontsize=9,
            color="black",
            verticalalignment="top",
            horizontalalignment="left",
            transform=fig.transFigure,
            clip_on=False,
            zorder=21,
        )

    def add_symbol(panel_name: str, x: float, top: float):
        """Add a large symbol at its original figure-relative position."""
        return axis.text(
            x,
            top,
            panel_name,
            fontsize=23,
            color="black",
            verticalalignment="top",
            transform=fig.transFigure,
            clip_on=False,
            zorder=22,
        )

    triang = value("triang")
    kappa = value("kappa")
    rmajor = value("rmajor")
    rminor = value("rminor")
    radius_plasma_core_norm = value("radius_plasma_core_norm")

    axis.axis("off")
    plot_plasma(axis, mfile, scan, colour_scheme)
    plasma_main_page(
        axis,
        fig,
        rmajor,
        rminor,
        triang,
        radius_plasma_core_norm,
        kappa,
        value("i_single_null"),
        value("plasma_square"),
        value("big_q_plasma"),
        value("p_fw_alpha_mw"),
        value("f_p_alpha_plasma_deposited"),
        value("p_neutron_total_mw"),
        value("pflux_plasma_surface_neutron_avg_mw"),
    )

    geom_type = PlasmaGeometryModelType(value("i_plasma_geometry"))
    textstr_plasma = (
        "$\\mathbf{Shaping:}$\n\n"
        f"$\\kappa_{{95}}$: {value('kappa95'):.2f} "
        f"({geom_type.kappa95_model.description}) | "
        f"$\\delta_{{95}}$: {value('triang95'):.2f} "
        f"({geom_type.triang95_model.description}) | "
        f"$\\zeta$: {value('plasma_square'):.2f}\n"
        f"$\\kappa$: {kappa:.2f} ({geom_type.kappa_model.description}) | "
        f"$\\delta$: {triang:.2f} ({geom_type.triang_model.description}) | "
        f"A: {value('aspect'):.2f}\n"
        f"$V_{{\\text{{p}}}}$: {value('vol_plasma'):,.2f} $\\mathrm{{m}}^3$ | "
        "$A_{\\text{p,surface}}$: "
        f"{value('a_plasma_surface'):,.2f} $\\mathrm{{m}}^2$ | "
        "$A_{\\text{p,poloidal}}$: "
        f"{value('a_plasma_poloidal'):,.3f} $\\mathrm{{m}}^2$\n"
        "$L_{\\text{p,poloidal}}$: "
        f"{value('len_plasma_poloidal'):,.3f} $\\mathrm{{m}}$"
    )
    add_panel(0.365, 0.975, 0.340, 0.110, textstr_plasma, box_style("lightyellow"))

    i_hcd_primary = value("i_hcd_primary")
    i_hcd_secondary = value("i_hcd_secondary")
    cd_1 = CurrentDriveModel(i_hcd_primary)
    cd_2 = CurrentDriveModel(i_hcd_secondary)
    textstr_hcd = (
        "$\\mathbf{Heating\\ &\\ current\\ drive:}$\n\n"
        f"Total injected heat: {value('p_hcd_injected_total_mw'):.3f} MW\n"
        f"Ohmic heating power: {value('p_plasma_ohmic_mw'):.3f} MW\n\n"
        f"$\\mathbf{{Primary\\ system: {cd_1.abbreviation}}}$\n"
        f"Current driving power: {value('p_hcd_primary_injected_mw'):.4f} MW\n"
        f"Extra heat power: {value('p_hcd_primary_extra_heat_mw'):.4f} MW\n"
        f"$\\eta_{{\\text{{CD,prim}}}}$: {value('eta_cd_hcd_primary'):.4f} A/W | "
        f"$\\langle\\zeta_{{\\text{{CD,prim}}}}\\rangle$: "
        f"{value('eta_cd_dimensionless_hcd_primary'):.4f}\n"
        f"$\\gamma_{{\\text{{CD,prim}}}}$: {value('eta_cd_norm_hcd_primary'):.4f} "
        "$\\times 10^{20}\\ \\mathrm{A}/\\mathrm{Wm}^2$\n"
        f"Current driven by primary: {value('c_hcd_primary_driven') / 1e6:.3f} MA\n\n"
        f"$\\mathbf{{Secondary\\ system: {cd_2.abbreviation}}}$\n"
        f"Current driving power: {value('p_hcd_secondary_injected_mw'):.4f} MW\n"
        f"Extra heat power: {value('p_hcd_secondary_extra_heat_mw'):.4f} MW\n"
        f"$\\eta_{{\\text{{CD,sec}}}}$: {value('eta_cd_hcd_secondary'):.4f} A/W | "
        f"$\\langle\\zeta_{{\\text{{CD,sec}}}}\\rangle$: "
        f"{value('eta_cd_dimensionless_hcd_secondary'):.4f}\n"
        f"$\\gamma_{{\\text{{CD,sec}}}}$: {value('eta_cd_norm_hcd_secondary'):.4f} "
        "$\\times 10^{20}\\ \\mathrm{A}/\\mathrm{Wm}^2$\n"
        f"Current driven by secondary: {value('c_hcd_secondary_driven') / 1e6:.3f} MA"
    )
    add_panel(
        0.730,
        0.675,
        0.245,
        0.300,
        textstr_hcd,
        box_style("paleturquoise") | {"edgecolor": "black"},
    )
    add_symbol("$P_{\\text{inj}}$", 0.92, 0.625)

    beta_component = BetaComponentLimits(value("i_beta_component")).full_name
    beta_norm_model = BetaNormMaxModel(value("i_beta_norm_max")).full_name
    textstr_beta = (
        "$\\mathbf{Beta\\ Information:}$\n\n"
        f"Total beta, $\\langle\\beta\\rangle$: {value('beta_total_vol_avg'):.4f}\n"
        f"Thermal beta, $\\langle\\beta_{{\\text{{thermal}}}}\\rangle$: "
        f"{value('beta_thermal_vol_avg'):.4f}\n"
        f"Toroidal beta, $\\langle\\beta_{{\\text{{t}}}}\\rangle$: "
        f"{value('beta_toroidal_vol_avg'):.4f}\n"
        f"Poloidal beta, $\\langle\\beta_{{\\text{{p}}}}\\rangle$: "
        f"{value('beta_poloidal_vol_avg'):.4f}\n"
        f"Fast-alpha beta, $\\langle\\beta_{{\\alpha}}\\rangle$: "
        f"{value('beta_fast_alpha'):.4f}\n"
        f"Upper limit on {beta_component}, $\\langle\\beta\\rangle$: "
        f"{value('beta_vol_avg_max'):.4f}\n"
        f"Normalised total beta, $\\beta_{{\\text{{N}}}}$: "
        f"{value('beta_norm_total'):.4f}\n"
        f"Normalised thermal beta, $\\beta_{{\\text{{N,thermal}}}}$: "
        f"{value('beta_norm_thermal'):.4f}\n"
        f"Maximum normalised beta ({beta_norm_model}), $\\beta_{{\\text{{N,max}}}}$: "
        f"{value('beta_norm_max'):.4f}"
    )
    add_panel(0.025, 0.975, 0.325, 0.180, textstr_beta, box_style("lightblue"))
    add_symbol("$\\beta$", 0.27, 0.94)

    inductance_model = IndInternalNormModel(
        value("i_ind_plasma_internal_norm")
    ).full_name
    textstr_volt_second = (
        "$\\mathbf{Volt-second\\ requirements:}$\n\n"
        f"Total volt-second consumption: {value('vs_plasma_total_required'):.4f} Vs\n"
        f"  - Internal volt-seconds: {value('vs_plasma_internal'):.4f} Vs\n"
        f"  - Volt-seconds needed for burn: {value('vs_plasma_burn_required'):.4f} Vs\n"
        f"  - Volt-seconds needed for ramp: {value('vs_plasma_ramp_required'):.4f} Vs | "
        f"$C_{{\\text{{ejima}}}}$: {value('ejima_coeff'):.4f}\n"
        f"$V_{{\\text{{loop}}}}$: {value('v_plasma_loop_burn'):.4f} V\n"
        f"$\\Omega_{{\\text{{p}}}}$: {value('res_plasma'):.4e} $\\Omega$\n"
        f"Plasma resistive diffusion time: {value('t_plasma_res_diffusion'):,.4f} s\n"
        f"Plasma inductance: {value('ind_plasma'):.4e} H | "
        f"ITER $l_i(3)$: {value('ind_plasma_internal_norm_iter_3'):.4f}\n"
        "Plasma stored magnetic energy: "
        f"{value('e_plasma_magnetic_stored') / 1e9:.4f} GJ\n"
        f"Plasma normalised internal inductance, $l_i$ ({inductance_model}): "
        f"{value('ind_plasma_internal_norm'):.3f}"
    )
    add_panel(0.025, 0.780, 0.335, 0.195, textstr_volt_second, box_style("lightgreen"))
    add_symbol("Vs", 0.30, 0.77)

    textstr_div = (
        f"$P_{{\\text{{sep}}}}$: {value('p_plasma_separatrix_mw'):.2f} MW\n"
        f"$\\frac{{P_{{\\text{{sep}}}}}}{{R}}$: "
        f"{value('p_plasma_separatrix_rmajor_mw'):.2f} MW/m\n"
        f"$\\frac{{P_{{\\text{{sep}}}}B_T}}{{q_{{95}} A R}}$: "
        f"{value('p_div_bt_q_aspect_rmajor_mw'):.2f} MW T/m"
    )
    add_panel(0.350, 0.120, 0.165, 0.095, textstr_div, box_style("orange"))
    add_symbol("$P_{\\text{div}}$", 0.45, 0.10)

    textstr_confinement = (
        "$\\mathbf{Confinement:}$\n\n"
        f"Confinement scaling law: {value('tauelaw')}\n"
        f"Confinement $H$ factor: {value('hfact'):.4f}\n"
        f"Energy confinement time from scaling: {value('t_energy_confinement'):.4f} s\n"
        f"Fusion double product: {value('ntau'):.4e} s/m³\n"
        f"Lawson Triple product: {value('nttau'):.4e} keV·s/m³\n"
        "Transport loss power assumed in scaling law: "
        f"{value('p_plasma_loss_mw'):.4f} MW\n"
        f"Plasma thermal energy (inc. $\\alpha$), $W$: "
        f"{value('e_plasma_beta') / 1e9:.4f} GJ\n"
        f"Alpha particle confinement time: {value('t_alpha_confinement'):.4f} s | "
        f"$\\tau_{{\\alpha}}/\\tau_{{e}}$: {value('f_t_alpha_energy_confinement'):.4f}"
    )
    add_panel(
        0.025,
        0.570,
        0.325,
        0.155,
        textstr_confinement,
        box_style("gainsboro") | {"edgecolor": "black"},
    )
    add_symbol("$\\tau_{\\text{e}}$", 0.30, 0.55)

    textstr_reactions = (
        "$\\mathbf{Fusion\\ Reactions:}$\n\n"
        "Fuel mixture:\n"
        f"| D: {value('f_plasma_fuel_deuterium'):.2f} | "
        f"T: {value('f_plasma_fuel_tritium'):.2f} | "
        f"3He: {value('f_plasma_fuel_helium3'):.2f} |\n\n"
        f"Fusion Power, $P_{{\\text{{fus}}}}$: {value('p_fusion_total_mw'):,.2f} MW\n"
        f"D-T Power, $P_{{\\text{{fus,DT}}}}$: {value('p_dt_total_mw'):,.2f} MW\n"
        f"D-D Power, $P_{{\\text{{fus,DD}}}}$: {value('p_dd_total_mw'):,.2f} MW\n"
        f"D-3He Power, $P_{{\\text{{fus,D3He}}}}$: {value('p_dhe3_total_mw'):,.2f} MW\n"
        f"Alpha Power, $P_{{\\alpha}}$: {value('p_alpha_total_mw'):,.2f} MW"
    )
    add_panel(
        0.025,
        0.400,
        0.185,
        0.165,
        textstr_reactions,
        {"boxstyle": "round", "facecolor": "red", "alpha": 0.6, "linewidth": 2},
    )

    textstr_fuelling = (
        "$\\mathbf{Fuelling:}$\n\n"
        f"Plasma mass: {value('m_plasma') * 1000:.4f} g\n"
        f"   - Average mass of all plasma ions: {value('m_ions_total_amu'):.3f} amu\n"
        f"Fuel mass: {value('m_plasma_fuel_ions') * 1000:.4f} g\n"
        f"   - Average mass of all fuel ions: {value('m_fuel_amu'):.3f} amu\n\n"
        "Fueling rate: "
        f"{value('molflow_plasma_fuelling_required'):.3e} nucleus-pairs/s\n"
        f"Fuel burn-up rate: {value('rndfuel'):.3e} reactions/s\n"
        f"Burn-up fraction: {value('burnup'):.4f}"
    )
    add_panel(
        0.025,
        0.220,
        0.255,
        0.175,
        textstr_fuelling,
        box_style("khaki") | {"edgecolor": "black"},
    )

    impurity_data = ImpurityRadiationData()
    textstr_ions = (
        "$\\mathbf{Ion\\ to\\ electron}$\n"
        "$\\mathbf{relative\\ number}$\n"
        "$\\mathbf{densities:}$\n\n"
        f"Effective charge: {value('n_charge_plasma_effective_vol_avg'):.3f}\n\n"
        + "\n".join(
            f"{label.replace('_', '') + ':':<6}"
            f"{value(f'f_nd_impurity_electrons({index:02d})'):.4e}"
            for index, label in enumerate(impurity_data.imp_label[:14], start=1)
        )
    )
    add_panel(
        0.805,
        0.335,
        0.170,
        0.310,
        textstr_ions,
        {"boxstyle": "round", "facecolor": "olivedrab", "alpha": 0.7, "linewidth": 2},
        text_x_offset=0.045,
    )
    add_symbol("$Z$", 0.815, 0.29)

    current_model = PlasmaCurrentModel(value("i_plasma_current")).full_name
    bootstrap_model = BootstrapCurrentFractionModel(
        value("i_bootstrap_current")
    ).full_name
    diamagnetic_model = PlasmaDiamagneticCurrentModel(
        value("i_diamagnetic_current")
    ).full_name
    textstr_currents = (
        "$\\mathbf{Plasma\\ currents:}$\n\n"
        f"Plasma current ({current_model}): {value('plasma_current_ma'):.4f} MA\n"
        f"  - Bootstrap fraction ({bootstrap_model}): "
        f"{value('f_c_plasma_bootstrap'):.4f}\n"
        f"  - Diamagnetic fraction ({diamagnetic_model}): "
        f"{value('f_c_plasma_diamagnetic'):.4f}\n"
        f"  - Pfirsch-Schlüter fraction: {value('f_c_plasma_pfirsch_schluter'):.4f}\n"
        f"  - Auxiliary fraction: {value('f_c_plasma_auxiliary'):.4f}\n"
        f"  - Inductive fraction: {value('f_c_plasma_inductive'):.4f}"
    )
    add_panel(0.720, 0.975, 0.255, 0.130, textstr_currents, box_style("#C8A2C8"))
    add_symbol("$I_{\\text{p}}$", 0.93, 0.90)

    textstr_fields = (
        "$\\mathbf{Magnetic\\ fields:}$\n\n"
        f"Toroidal field at $R_0$, $B_T$: {value('b_plasma_toroidal_on_axis'):.4f} T\n"
        f"  Ripple at outboard, $\\delta$: {value('ripple_b_tf_plasma_edge'):.2f}%\n"
        f"Surface average poloidal field, $\\langle B_p(a)\\rangle$: "
        f"{value('b_plasma_surface_poloidal_average'):.4f} T\n"
        f"Total field, $B_{{\\text{{tot}}}}$: {value('b_plasma_total'):.4f} T\n"
        f"Vertical field, $B_{{\\text{{vert}}}}$: "
        f"{value('b_plasma_vertical_required'):.4f} T"
    )
    add_panel(0.5325, 0.140, 0.2575, 0.115, textstr_fields, box_style("royalblue"))
    add_symbol("$B$", 0.75, 0.12)

    textstr_radiation = (
        "$\\mathbf{Radiation:}$\n\n"
        f"Total radiation power: {value('p_plasma_rad_mw'):.4f} MW\n"
        f"Separatrix radiation fraction: {value('f_p_plasma_separatrix_rad'):.4f}\n"
        f"Core radiation power: {value('p_plasma_inner_rad_mw'):.4f} MW\n"
        "   - $f_{\\text{core,reduce}}$: "
        f"{value('f_p_plasma_core_rad_reduction'):.4f}\n"
        f"Edge radiation power: {value('p_plasma_outer_rad_mw'):.4f} MW\n"
        f"Synchrotron radiation power: {value('p_plasma_sync_mw'):.4f} MW\n"
        f"Synchrotron wall reflectivity: {value('f_sync_reflect'):.4f}"
    )
    add_panel(
        0.720,
        0.830,
        0.255,
        0.145,
        textstr_radiation,
        box_style("lavender") | {"edgecolor": "black"},
        text_x_offset=0.032,
    )
    add_symbol("$\\gamma$", 0.725, 0.78)

    model_name = PlasmaConfinementTransitionModel(value("i_l_h_threshold")).full_name
    if len(model_name) > 20:
        model_name = "\n".join(textwrap.wrap(model_name, width=20))
    textstr_lh = (
        "$\\mathbf{L-H\\ threshold:}$\n"
        f"{model_name}\n\n"
        f"$P_{{\\text{{L-H}}}}$: {value('p_l_h_threshold_mw'):.4f} MW"
    )
    add_panel(0.220, 0.400, 0.110, 0.080, textstr_lh, box_style("peachpuff"))

    density_model = DensityLimitModel(value("i_density_limit")).full_name
    textstr_density_limit = (
        "$\\mathbf{Density\\ limit:}$\n"
        f"({density_model})\n"
        f"$n_{{\\text{{e,limit}}}}$: {value('nd_plasma_electrons_max'):.3e} "
        "$\\mathrm{m}^{-3}$\n"
        f"$f_{{\\text{{GW}}}}$: {value('f_nd_plasma_greenwald'):.4f}"
    )
    add_panel(0.220, 0.310, 0.130, 0.075, textstr_density_limit, box_style("pink"))


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

    light_yellow_box = box_style("lightyellow")

    axis.text(
        0.05,
        0.45,
        textstr_debye,
        **text_layout(fig, v_align="top"),
        bbox=light_yellow_box,
    )

    axis.text(
        0.25,
        0.45,
        textstr_larmor,
        **text_layout(fig, v_align="top"),
        bbox=light_yellow_box,
    )

    axis.text(
        0.45,
        0.45,
        textstr_velocities,
        **text_layout(fig, v_align="top"),
        bbox=light_yellow_box,
    )

    light_cyan_box = box_style("lightcyan")
    axis.text(
        0.05,
        0.31,
        textstr_frequencies,
        **text_layout(fig, v_align="top"),
        bbox=light_cyan_box,
    )

    axis.text(
        0.25,
        0.31,
        textstr_coulomb,
        **text_layout(fig, v_align="top"),
        bbox=light_cyan_box,
    )

    axis.text(
        0.45,
        0.31,
        textstr_collision_times,
        **text_layout(fig, v_align="top"),
        bbox=light_cyan_box,
    )

    light_green_box = box_style("lightgreen")

    axis.text(
        0.05,
        0.17,
        textstr_collision_freq,
        **text_layout(fig, v_align="top"),
        bbox=light_green_box,
    )

    axis.text(
        0.25,
        0.17,
        textstr_mfp,
        **text_layout(fig, v_align="top"),
        bbox=light_green_box,
    )

    axis.text(
        0.45,
        0.17,
        textstr_spitzer + "\n" + textstr_resistivity,
        **text_layout(fig, v_align="top"),
        bbox=light_green_box,
    )

    axis.axis("off")


__all__ = ["plot_detailed_plasma_parameters", "plot_main_plasma_information"]
