"""Scrape off layer functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary.constants import PLASMA_COLOUR
from process.core.io.plot.summary.plasma.physics import plot_plasma
from process.data_structure.physics_variables import OutbordSOLPowerDecayLengthModel
from process.models.physics.exhaust import calculate_brunner_divertor_power_splits
from process.models.physics.scrape_off_layer import ScrapeOffLayer

if TYPE_CHECKING:
    from process.core.io.mfile import MFile

__all__ = [
    "plot_brunner_divertor_power_split_comparison_stackplot",
    "plot_div_lower_outboard_eich_target_profile",
    "plot_midplane_near_sol_radial_profile",
    "plot_sol_power_decay_length_comparison",
    "plot_sol_power_flux_profiles",
]


def plot_midplane_near_sol_radial_profile(
    axis: plt.Axes, mfile: MFile, scan: int, colour_scheme: int
):
    """Function to plot the radial profile of the near SOL at the midplane."""
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    len_sol_outboard_power_decay = mfile.get("len_sol_outboard_power_decay", scan=scan)
    r = np.linspace(
        (rmajor + rminor), (rmajor + rminor) + (7 * len_sol_outboard_power_decay), 100
    )

    radial_profile = (
        ScrapeOffLayer().calculate_outboard_midplane_near_sol_radial_profile(
            rmajor=rmajor,
            rminor=rminor,
            len_plasma_sol_power_decay=len_sol_outboard_power_decay,
            pflux_plasma_outboard_sol_parallel_mw=mfile.get(
                "pflux_plasma_outboard_sol_parallel_mw", scan=scan
            ),
            r=r,
        )
    )

    x_ref = rmajor + rminor + len_sol_outboard_power_decay
    plasma_boundary = rmajor + rminor

    axis.axvspan(
        0,
        plasma_boundary,
        color=PLASMA_COLOUR[colour_scheme - 1],
        alpha=0.35,
        label="Plasma",
    )
    axis.plot(r, radial_profile, label=r"$q_{||}$ profile")
    axis.axvline(
        x_ref,
        color="black",
        linestyle="--",
        linewidth=1,
        label=r"$\lambda_{q,\mathrm{out}}$",
    )
    axis.grid()
    axis.legend()
    axis.set_xlim(
        (rmajor + rminor - (3 * len_sol_outboard_power_decay)),
        (rmajor + rminor) + (7 * len_sol_outboard_power_decay),
    )
    axis.set_title(r"Midplane Near SOL Radial Profile")
    axis.minorticks_on()
    axis.tick_params(axis="x", labelbottom=False)
    axis.set_ylabel(r"$q_{||}$ [MW/m$^2$]")


def plot_div_lower_outboard_eich_target_profile(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot the Eich target profile at the lower outboard divertor."""
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    len_plasma_sol_power_decay = mfile.get(
        "len_plasma_sol_eich13_power_decay", scan=scan
    )
    len_div_outboard_lower_power_spreading = mfile.get(
        "len_div_outboard_lower_power_spreading", scan=scan
    )
    f_b_flux_expansion = mfile.get("f_b_div_outboard_lower_flux_expansion", scan=scan)
    r = np.linspace(
        (rmajor + rminor)
        - ((f_b_flux_expansion / 2) * len_plasma_sol_power_decay) * f_b_flux_expansion,
        (rmajor + rminor) + (3 * len_plasma_sol_power_decay) * f_b_flux_expansion,
        200,
    )

    pflux_target_profile = ScrapeOffLayer().calculate_eich_target_heat_flux_profile(
        rmajor=rmajor,
        rminor=rminor,
        pflux_plasma_sol_parallel_mw=mfile.get(
            "pflux_plasma_outboard_sol_parallel_mw", scan=scan
        ),
        len_plasma_sol_power_decay=mfile.get("len_sol_outboard_power_decay", scan=scan),
        f_b_div_flux_expansion=f_b_flux_expansion,
        len_plasma_sol_power_spreading=len_div_outboard_lower_power_spreading,
        pflux_target_background_heat_flux_mw=0.0,
        r=r,
    )
    peak_idx = np.argmax(pflux_target_profile)
    peak_r = r[peak_idx]
    peak_q = pflux_target_profile[peak_idx]

    axis.plot(r, pflux_target_profile)
    axis.axvline(peak_r, color="black", linestyle="--", linewidth=1)
    axis.axhline(peak_q, color="black", linestyle="--", linewidth=1)
    axis.text(
        0.8,
        0.9,
        f"$f_x$ = {f_b_flux_expansion:.2f}\n$S$ = "
        f"{len_div_outboard_lower_power_spreading * 1e3:.3f} mm",
        transform=axis.transAxes,
        ha="left",
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1.0},
    )
    axis.grid()
    axis.minorticks_on()
    axis.set_title(r"Lower Outboard Eich Target Unmitigated Parallel Heat Flux Profile")
    axis.set_xlabel("Radial Position [m]")
    axis.set_xlim(r[0], r[-1])
    axis.set_ylabel(r"$q_{||,t,\mathrm{unmitigated}}$ [MW/m$^2$]")


def plot_sol_power_flux_profiles(axis: plt.Axes, mfile: MFile, scan: int, colour_scheme):
    """Plot separatrix power split fractions as a bar chart."""
    plot_plasma(axis=axis, mfile=mfile, scan=scan, colour_scheme=colour_scheme)
    rmajor, rminor, kappa = mfile.get_variables(
        "rmajor",
        "rminor",
        "kappa",
        scan=scan,
    )
    len_sol_outboard_power_decay = mfile.get("len_sol_outboard_power_decay", scan=scan)
    a_plasma_outboard_sol_parallel = mfile.get(
        "a_plasma_outboard_sol_parallel", scan=scan
    )
    pflux_plasma_outboard_sol_parallel_mw = mfile.get(
        "pflux_plasma_outboard_sol_parallel_mw", scan=scan
    )
    plasma_scale = max(rminor, abs(kappa * rminor), 1e-6)
    scale_factor = min(max(plasma_scale / 2.0, 0.7), 1.0)
    text_fontsize = 9 * scale_factor
    sol_power_decay_model_description = OutbordSOLPowerDecayLengthModel(
        mfile.get("i_len_sol_outboard_power_decay", scan=scan)
    ).description

    outboard_pos = (rmajor + rminor, 0.0)

    axis.text(
        *outboard_pos,
        f"$\\lambda_q = {len_sol_outboard_power_decay * 1e3:.3f}$ mm "
        f"({sol_power_decay_model_description})\n"
        f"$A_{{||}} = {a_plasma_outboard_sol_parallel:.4f}$ m$^2$\n"
        f"$q_{{||}} = {pflux_plasma_outboard_sol_parallel_mw:,.2f}$ MW/m$^2$",
        fontsize=text_fontsize,
        verticalalignment="center",
        horizontalalignment="center",
        bbox={
            "boxstyle": f"round,pad={0.3 * scale_factor:.3f}",
            "alpha": 1.0,
            "linewidth": 2 * scale_factor,
            "edgecolor": "black",
        },
        zorder=101,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.get_xaxis().set_ticks([])
    axis.get_yaxis().set_ticks([])


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
    axis.text(
        1.02,
        0.2,
        f"Average: {avg_decay_length:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
        1.02,
        0.15,
        f"Standard Dev: {std_decay_length:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    axis.text(
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
    len_plasma_sol_outboard_pd = mfile.get("len_sol_outboard_power_decay", scan=scan)
    len_plasma_sol_inboard_pd = mfile.get("len_sol_inboard_power_decay", scan=scan)
    # Use the case decay length when available; fall back to 1 mm if absent.
    if not np.isfinite(len_plasma_sol_outboard_pd) or len_plasma_sol_outboard_pd == 0.0:  # noqa: RUF069
        len_plasma_sol_outboard_pd = 0.001

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
