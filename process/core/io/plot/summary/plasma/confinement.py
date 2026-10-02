"""Plasma functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary.rendering import draw_text
from process.data_structure.physics_variables import (
    ConfinementTimeModel,
    OutbordSOLPowerDecayLengthModel,
)
from process.models.physics.confinement_time import PlasmaConfinementTime
from process.models.physics.exhaust import calculate_brunner_divertor_power_splits

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_sol_power_decay_length_comparison(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot a scatter box plot of SOL power decay lengths (λ_q).

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    """
    len_plasma_sol_eich13_power_decay_mm = (
        mfile.get("len_plasma_sol_eich13_power_decay", scan=scan) * 1e3
    )
    len_plasma_sol_mast14_power_decay_1_mm = (
        mfile.get("len_plasma_sol_mast14_power_decay_1", scan=scan) * 1e3
    )
    len_plasma_sol_mast14_power_decay_2_mm = (
        mfile.get("len_plasma_sol_mast14_power_decay_2", scan=scan) * 1e3
    )
    len_plasma_sol_eich11_jet_power_decay_mm = (
        mfile.get("len_plasma_sol_eich11_jet_power_decay", scan=scan) * 1e3
    )
    len_plasma_sol_eich11_jet_asdex_power_decay_mm = (
        mfile.get("len_plasma_sol_eich11_jet_asdex_power_decay", scan=scan) * 1e3
    )
    # Data for the box plot
    data = {
        f"{OutbordSOLPowerDecayLengthModel.EICH_2013.description}": (
            len_plasma_sol_eich13_power_decay_mm
        ),
        f"{OutbordSOLPowerDecayLengthModel.MAST_2014_1.description}": (
            len_plasma_sol_mast14_power_decay_1_mm
        ),
        f"{OutbordSOLPowerDecayLengthModel.MAST_2014_2.description}": (
            len_plasma_sol_mast14_power_decay_2_mm
        ),
        f"{OutbordSOLPowerDecayLengthModel.EICH_2011_JET.description}": (
            len_plasma_sol_eich11_jet_power_decay_mm
        ),
        f"{OutbordSOLPowerDecayLengthModel.EICH_2011_JET_ASDEX.description}": (
            len_plasma_sol_eich11_jet_asdex_power_decay_mm
        ),
    }
    data_values = list(data.values())

    # Create the violin plot
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data_values)))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Calculate average, standard deviation, and median
    avg_decay_length = np.mean(data_values)
    std_decay_length = np.std(data_values)
    median_decay_length = np.median(data_values)

    # Plot average, standard deviation, and median as text
    draw_text(
        axis,
        1.02,
        0.2,
        f"Average: {avg_decay_length:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.02,
        0.15,
        f"Standard Dev: {std_decay_length:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.02,
        0.1,
        f"Median: {median_decay_length:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("SOL Power Decay Length ($\\lambda_q$) Comparison")
    axis.set_ylabel("Power Decay Length [mm]")
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def plot_brunner_divertor_power_split_comparison_stackplot(
    axis: plt.Axes, mfile: MFile, scan: int
):
    """Plot Brunner divertor power split fractions as a stack plot over dr_sep."""
    # Use the case decay length when available; fall back to 1 mm if absent.

    len_plasma_sol_outboard_pd = mfile.get("len_sol_outboard_power_decay", scan=scan)
    len_plasma_sol_inboard_pd = mfile.get("len_sol_inboard_power_decay", scan=scan)
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, 4))

    dr_sep_values = np.linspace(
        -5 * len_plasma_sol_outboard_pd,
        5 * len_plasma_sol_outboard_pd,
        200,
    )
    f_p_inboard_lower = np.zeros_like(dr_sep_values)
    f_p_inboard_upper = np.zeros_like(dr_sep_values)
    f_p_outboard_lower = np.zeros_like(dr_sep_values)
    f_p_outboard_upper = np.zeros_like(dr_sep_values)

    for idx, dr_sep in enumerate(dr_sep_values):
        div_power_splits = calculate_brunner_divertor_power_splits(
            dr_outboard_midplane_sep=dr_sep,
            len_plasma_sol_outboard_power_decay=len_plasma_sol_outboard_pd,
            len_plasma_sol_inboard_power_decay=len_plasma_sol_inboard_pd,
        )
        f_p_inboard_lower[idx] = div_power_splits.f_p_div_inboard_lower_separatrix
        f_p_inboard_upper[idx] = div_power_splits.f_p_div_inboard_upper_separatrix
        f_p_outboard_lower[idx] = div_power_splits.f_p_div_outboard_lower_separatrix
        f_p_outboard_upper[idx] = div_power_splits.f_p_div_outboard_upper_separatrix

    axis.stackplot(
        dr_sep_values,
        f_p_inboard_lower,
        f_p_inboard_upper,
        f_p_outboard_lower,
        f_p_outboard_upper,
        labels=[
            "$f_{P,\\mathrm{in,lower}}$",
            "$f_{P,\\mathrm{in,upper}}$",
            "$f_{P,\\mathrm{out,lower}}$",
            "$f_{P,\\mathrm{out,upper}}$",
        ],
        colors=colors,
        alpha=0.9,
    )

    axis.axvline(
        mfile.get("dr_plasma_outboard_midplane_separatrix_separation", scan=scan),
        color="k",
        linestyle="--",
        linewidth=1.0,
        alpha=0.5,
        label="$\u0394 r_{\\mathrm{sep}}$",
    )
    axis.set_ylim(0.0, 1.0)
    axis.set_xlim(-5 * len_plasma_sol_outboard_pd, 5 * len_plasma_sol_outboard_pd)
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.35)
    axis.set_title("Brunner Divertor Power Split Fractions")
    axis.set_xlabel("$\\Delta r_{\\mathrm{sep}}$ [m]")
    axis.set_ylabel("Power split fraction, $f_P$")
    axis.legend(loc="upper left", fontsize=8)


def plot_confinement_time_comparison(
    axis: plt.Axes, mfile: MFile, scan: int, u_seed=None
):
    """Function to plot a scatter box plot of confinement time comparisons.

    Parameters
    ----------
    axis :
        Axis object to plot to.
    mfile :
         MFILE data object.
    scan :
        Scan number to use.
    u_seed :
         (Default value = None)
    """
    rminor = mfile.get("rminor", scan=scan)
    rmajor = mfile.get("rmajor", scan=scan)
    cur_plasma_ma = mfile.get("plasma_current_ma", scan=scan)
    kappa95 = mfile.get("kappa95", scan=scan)
    nd_plasma_electron_line_20 = mfile.get("nd_plasma_electron_line", scan=scan) / 1e20
    afuel = mfile.get("m_fuel_amu", scan=scan)
    b_plasma_toroidal_on_axis = mfile.get("b_plasma_toroidal_on_axis", scan=scan)
    p_plasma_separatrix_mw = mfile.get("p_plasma_separatrix_mw", scan=scan)
    kappa = mfile.get("kappa", scan=scan)
    aspect = mfile.get("aspect", scan=scan)
    nd_plasma_electron_line_19 = mfile.get("nd_plasma_electron_line", scan=scan) / 1e19
    kappa_ipb = mfile.get("kappa_ipb", scan=scan)
    triang = mfile.get("triang", scan=scan)
    m_ions_total_amu = mfile.get("m_ions_total_amu", scan=scan)

    confine = PlasmaConfinementTime()

    # Calculate confinement times using the scan data
    iter_89p = confine.iter_89p_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        nd_plasma_electron_line_20=nd_plasma_electron_line_20,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_89_0 = confine.iter_89_0_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        nd_plasma_electron_line_20=nd_plasma_electron_line_20,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_h90_p = confine.iter_h90_p_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        rmajor=rmajor,
        rminor=rminor,
        kappa=kappa,
        nd_plasma_electron_line_20=nd_plasma_electron_line_20,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_h90_p_amended = confine.iter_h90_p_amended_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        afuel=afuel,
        rmajor=rmajor,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        kappa=kappa,
    )
    iter_93h = confine.iter_93h_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        afuel=afuel,
        rmajor=rmajor,
        nd_plasma_electron_line_20=nd_plasma_electron_line_20,
        aspect=aspect,
        kappa=kappa,
    )
    iter_h97p = confine.iter_h97p_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        rmajor=rmajor,
        aspect=aspect,
        kappa=kappa,
        afuel=afuel,
    )
    iter_h97p_elmy = confine.iter_h97p_elmy_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        rmajor=rmajor,
        aspect=aspect,
        kappa=kappa,
        afuel=afuel,
    )
    iter_96p = confine.iter_96p_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        kappa95=kappa95,
        rmajor=rmajor,
        aspect=aspect,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        afuel=afuel,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
    )
    iter_pb98py = confine.iter_pb98py_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa=kappa,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y = confine.iter_ipb98y_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa=kappa,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y1 = confine.iter_ipb98y1_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y2 = confine.iter_ipb98y2_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y3 = confine.iter_ipb98y3_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    iter_ipb98y4 = confine.iter_ipb98y4_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    petty08 = confine.petty08_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
    )
    menard_nstx = confine.menard_nstx_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    menard_nstx_petty08 = confine.menard_nstx_petty08_hybrid_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        kappa_ipb=kappa_ipb,
        aspect=aspect,
        afuel=afuel,
    )
    itpa20 = confine.itpa20_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        triang=triang,
        kappa_ipb=kappa_ipb,
        eps=(1 / aspect),
        aion=m_ions_total_amu,
    )
    itpa20_ilc = confine.itpa20_il_confinement_time(
        cur_plasma_ma=cur_plasma_ma,
        b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
        p_plasma_loss_mw=p_plasma_separatrix_mw,
        nd_plasma_electron_line_19=nd_plasma_electron_line_19,
        aion=m_ions_total_amu,
        rmajor=rmajor,
        triang=triang,
        kappa_ipb=kappa_ipb,
    )

    # Data for the box plot
    data = {
        rf"{ConfinementTimeModel.ITER_89P.full_name}": iter_89p,
        rf"{ConfinementTimeModel.ITER_89_0.full_name}": iter_89_0,
        rf"{ConfinementTimeModel.ITER_H90_P.full_name}": iter_h90_p,
        rf"{ConfinementTimeModel.ITER_H90_P_AMENDED.full_name}": (iter_h90_p_amended),
        rf"{ConfinementTimeModel.ITER_93H.full_name}": iter_93h,
        rf"{ConfinementTimeModel.ITER_H97P.full_name}": iter_h97p,
        rf"{ConfinementTimeModel.ITER_H97P_ELMY.full_name}": iter_h97p_elmy,
        rf"{ConfinementTimeModel.ITER_96P.full_name}": iter_96p,
        rf"{ConfinementTimeModel.ITER_PB98P_Y.full_name}": iter_pb98py,
        rf"{ConfinementTimeModel.IPB98_Y.full_name}": iter_ipb98y,
        rf"{ConfinementTimeModel.ITER_IPB98Y1.full_name}": iter_ipb98y1,
        rf"{ConfinementTimeModel.ITER_IPB98Y2.full_name}": iter_ipb98y2,
        rf"{ConfinementTimeModel.ITER_IPB98Y3.full_name}": iter_ipb98y3,
        rf"{ConfinementTimeModel.ITER_IPB98Y4.full_name}": iter_ipb98y4,
        rf"{ConfinementTimeModel.PETTY08.full_name}": petty08,
        rf"{ConfinementTimeModel.MENARD_NSTX.full_name}": menard_nstx,
        rf"{ConfinementTimeModel.MENARD_NSTX_PETTY08_HYBRID.full_name}": (
            menard_nstx_petty08
        ),
        rf"{ConfinementTimeModel.ITPA20.full_name}": itpa20,
        rf"{ConfinementTimeModel.ITPA20_IL.full_name}": itpa20_ilc,
    }
    data_values = list(data.values())

    # Create the violin plot
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

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
    avg_threshold = np.mean(data_values)
    std_threshold = np.std(data_values)
    median_threshold = np.median(data_values)

    # Plot average, standard deviation, and median as text
    draw_text(
        axis,
        0.7,
        1.25,
        f"Average: {avg_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        0.7,
        1.2,
        f"Standard Dev: {std_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        0.7,
        1.15,
        f"Median: {median_threshold:.4f} s",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        0.75,
        -0.05,
        r"$H \ factor = 1.0$",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Confinement time ($\\tau_{\\text{E}}$) Comparison")
    axis.set_ylabel("Confinement time, $\\tau_{\\text{E}}$ [s]")
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])

    # Add background color
    axis.set_facecolor("#f0f0f0")


__all__ = [
    "plot_brunner_divertor_power_split_comparison_stackplot",
    "plot_confinement_time_comparison",
    "plot_sol_power_decay_length_comparison",
]
