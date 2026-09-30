"""Plasma functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary.common import (
    setup_axis,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)
from process.core.io.plot.summary.reporting import (
    plot_info,
)

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_current_drive_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot current drive info

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    setup_axis(axis, xmin=0, xmax=1, ymin=-16, ymax=1)

    i_hcd_primary = mfile.get("i_hcd_primary", scan=scan)

    if nbi := (i_hcd_primary in {5, 8}):
        draw_text(
            axis,
            -0.05,
            1,
            "Neutral Beam Current Drive:",
            ha="left",
            va="center",
        )
    if ecrh := (i_hcd_primary in {3, 7, 10, 11, 13}):
        draw_text(
            axis,
            -0.05,
            1,
            "Electron Cyclotron Current Drive:",
            ha="left",
            va="center",
        )
    if ebw := (i_hcd_primary == 12):
        draw_text(
            axis,
            -0.05,
            1,
            "Electron Bernstein Wave Drive:",
            ha="left",
            va="center",
        )
    if lhcd := (i_hcd_primary in {1, 4, 6}):
        draw_text(
            axis,
            -0.05,
            1,
            "Lower Hybrid Current Drive:",
            ha="left",
            va="center",
        )
    if iccd := (i_hcd_primary == 2):
        draw_text(
            axis,
            -0.05,
            1,
            "Ion Cyclotron Current Drive:",
            ha="left",
            va="center",
        )

    i_hcd_secondary = mfile.get("i_hcd_secondary", scan=scan) or 0
    if i_hcd_secondary in {5, 8}:
        secondary_heating = "NBI"
    elif i_hcd_secondary in {3, 7, 10, 11, 13}:
        secondary_heating = "ECH"
    elif i_hcd_secondary == 12:
        secondary_heating = "EBW"
    elif i_hcd_secondary in {1, 4, 6}:
        secondary_heating = "LHCD"
    elif i_hcd_secondary == 2:
        secondary_heating = "ICCD"
    else:
        secondary_heating = ""

    pinjie = mfile.get("p_hcd_injected_total_mw", scan=scan)
    p_plasma_separatrix_mw = mfile.get("p_plasma_separatrix_mw", scan=scan)
    pdivr = p_plasma_separatrix_mw / mfile.get("rmajor", scan=scan)

    if mfile.get("i_hcd_secondary", scan=scan) != 0:
        pinjmwfix = mfile.get("pinjmwfix", scan=scan)

    pdivnr = (
        1.0e20
        * mfile.get("p_plasma_separatrix_mw", scan=scan)
        / (
            mfile.get("rmajor", scan=scan)
            * mfile.get("nd_plasma_electrons_vol_avg", scan=scan)
        )
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    pthresh_name = (
        "p_l_h_threshold_mw"
        if "p_l_h_threshold_mw" in mfile.data
        else "l_h_threshold_powers(6)"
    )
    pthresh = mfile.get(pthresh_name, scan=scan)
    flh = p_plasma_separatrix_mw / pthresh

    hstar = mfile.get("hstar", scan=scan)

    data = [
        (pinjie, "Steady state auxiliary power", "MW"),
        ("p_hcd_primary_extra_heat_mw", "Power for heating only", "MW"),
        ("f_c_plasma_bootstrap", "Bootstrap fraction", ""),
        ("f_c_plasma_auxiliary", "Auxiliary fraction", ""),
        ("f_c_plasma_inductive", "Inductive fraction", ""),
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
    # Optional override based on condition
    field_overrides = {
        "ecrh": (
            "eta_cd_hcd_primary",
            r"$\frac{P_{\mathrm{div}}}{R_{0}}$",
            "A W$^{-1}$",
        ),
        "nbi": (
            ("gamnb", "NB gamma", "$10^{20}$ A W$^{-1}$ m$^{-2}$"),
            ("e_beam_kev", "NB energy", "keV"),
        ),
        "ebw": (
            "eta_cd_norm_hcd_primary",
            "Normalised current drive efficiency of primary HCD system",
            "(10$^{20}$ A/(Wm$^{2}$))",
        ),
        "lhcd": (
            "eta_cd_norm_hcd_primary",
            "Normalised current drive efficiency",
            "(10$^{20}$ A/(Wm$^{2}$))",
        ),
        "iccd": (
            "eta_cd_norm_hcd_primary",
            "Normalised current drive efficiency",
            "(10$^{20}$ A/(Wm$^{2}$))",
        ),
    }

    if ecrh:
        data.insert(6, field_overrides["ecrh"])
    elif nbi:
        data.insert(6, field_overrides["nbi"][0])
        data.insert(7, field_overrides["nbi"][1])
    elif ebw:
        data.insert(6, field_overrides["ebw"])
    elif lhcd:
        data.insert(6, field_overrides["lhcd"])
    elif iccd:
        data.insert(6, field_overrides["iccd"])

    # Secondary heating logic — common across all cases
    if mfile.get("i_hcd_secondary", scan=scan) != 0:
        data.insert(
            1,
            (
                "pinjmwfix",
                f"{secondary_heating} secondary auxiliary power",
                "MW",
            ),
        )
        data[0] = ((pinjie - pinjmwfix), "Primary auxiliary power", "MW")
        data.insert(2, (pinjie, "Total auxillary power", "MW"))

    coe = mfile.get("coe", scan=scan)
    data.extend((
        ("", "", ""),
        ("#Costs", "", ""),
        (
            ("", "Cost output not selected", "")
            if coe == 0.0  # noqa: RUF069
            else (coe, "Cost of electricity", r"\$/MWh")
        ),
    ))

    plot_info(axis, data, mfile, scan)


def plot_bootstrap_comparison(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot a scatter box plot of bootstrap current fractions.

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    """
    # Data for the box plot
    data = {
        "IPDG": mfile.get("f_c_plasma_bootstrap_iter89", scan=scan),
        "Sauter": mfile.get("f_c_plasma_bootstrap_sauter", scan=scan),
        "Nevins": mfile.get("f_c_plasma_bootstrap_nevins", scan=scan),
        "Wilson": mfile.get("f_c_plasma_bootstrap_wilson", scan=scan),
        "Sakai": mfile.get("f_c_plasma_bootstrap_sakai", scan=scan),
        "ARIES": mfile.get("f_c_plasma_bootstrap_aries", scan=scan),
        "Andrade": mfile.get("f_c_plasma_bootstrap_andrade", scan=scan),
        "Hoang": mfile.get("f_c_plasma_bootstrap_hoang", scan=scan),
        "Wong": mfile.get("f_c_plasma_bootstrap_wong", scan=scan),
        "Gi-I": mfile.get("bscf_gi_i", scan=scan),
        "Gi-II": mfile.get("bscf_gi_ii", scan=scan),
        "Sugiyama (L-mode)": mfile.get("f_c_plasma_bootstrap_sugiyama_l", scan=scan),
        "Sugiyama (H-mode)": mfile.get("f_c_plasma_bootstrap_sugiyama_h", scan=scan),
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
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Calculate average, standard deviation, and median
    avg_bootstrap = np.mean(data_values)
    std_bootstrap = np.std(data_values)
    median_bootstrap = np.median(data_values)

    # Plot average, standard deviation, and median as text
    draw_text(
        axis,
        1.02,
        0.2,
        f"Average: {avg_bootstrap:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.02,
        0.15,
        f"Standard Dev: {std_bootstrap:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )
    draw_text(
        axis,
        1.02,
        0.1,
        f"Median: {median_bootstrap:.4f}",
        transform=axis.transAxes,
        fontsize=9,
    )

    axis.set_title("Bootstrap Current Fraction ($f_\\text{BS}$) Comparison")
    axis.set_ylabel("Bootstrap Current Fraction")
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


__all__ = ["plot_bootstrap_comparison", "plot_current_drive_info"]
