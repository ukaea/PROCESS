"""Reporting functions for PROCESS summary plots."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary.constants import (
    PLASMA_COLOUR,
    SHIELD_COLOUR,
    TFC_COLOUR,
    THERMAL_SHIELD_COLOUR,
    VESSEL_COLOUR,
)
from process.core.io.plot.summary.reporting.layouts import (
    draw_bend,
)
from process.models.physics.current_drive import (
    ElectronBernstein,
    ElectronCyclotron,
)

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


@dataclass
class RadialBuild:
    """Dataclass containing radial build dictionaries"""

    upper: dict[str, float]
    lower: dict[str, float]
    radial: dict[str, float]

    cumulative_upper: dict[str, float]
    cumulative_lower: dict[str, float]
    cumulative_radial: dict[str, float]


def plot_centre_cross(
    axis: plt.Axes, mfile: MFile, scan: int, mirror_negative_x: bool = False
):
    """Function to plot centre cross on plot

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE data object
    scan :
        scan number to use
    mirror_negative_x :
        if True, mirror the plot to the negative x-axis (Default value = False)
    """
    rmajor = mfile.get("rmajor", scan=scan)
    x_scale = -1 if mirror_negative_x else 1
    axis.plot(
        x_scale * np.array([rmajor - 0.25, rmajor + 0.25, rmajor, rmajor, rmajor]),
        [0, 0, 0, 0.25, -0.25],
        color="black",
    )


def plot_h_threshold_comparison(axis: plt.Axes, mfile: MFile, scan: int, u_seed=None):
    """Function to plot a scatter box plot of L-H threshold power comparisons.

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
    # Data for the box plot
    data = {
        "ITER 1996 Nominal": mfile.get("l_h_threshold_powers(1)", scan=scan),
        "ITER 1996 Upper": mfile.get("l_h_threshold_powers(2)", scan=scan),
        "ITER 1996 Lower": mfile.get("l_h_threshold_powers(3)", scan=scan),
        "ITER 1997 (1)": mfile.get("l_h_threshold_powers(4)", scan=scan),
        "ITER 1997 (2)": mfile.get("l_h_threshold_powers(5)", scan=scan),
        "Martin Nominal": mfile.get("l_h_threshold_powers(6)", scan=scan),
        "Martin Upper": mfile.get("l_h_threshold_powers(7)", scan=scan),
        "Martin Lower": mfile.get("l_h_threshold_powers(8)", scan=scan),
        "Snipes Nominal": mfile.get("l_h_threshold_powers(9)", scan=scan),
        "Snipes Upper": mfile.get("l_h_threshold_powers(10)", scan=scan),
        "Snipes Lower": mfile.get("l_h_threshold_powers(11)", scan=scan),
        "Snipes Closed Divertor Nominal": mfile.get(
            "l_h_threshold_powers(12)", scan=scan
        ),
        "Snipes Closed Divertor Upper": mfile.get("l_h_threshold_powers(13)", scan=scan),
        "Snipes Closed Divertor Lower": mfile.get("l_h_threshold_powers(14)", scan=scan),
        "Hubbard Nominal (I-mode)": mfile.get("l_h_threshold_powers(15)", scan=scan),
        "Hubbard Lower (I-mode)": mfile.get("l_h_threshold_powers(16)", scan=scan),
        "Hubbard Upper (I-mode)": mfile.get("l_h_threshold_powers(17)", scan=scan),
        "Hubbard 2017 (I-mode)": mfile.get("l_h_threshold_powers(18)", scan=scan),
        "Martin Aspect Corrected Nominal": mfile.get(
            "l_h_threshold_powers(19)", scan=scan
        ),
        "Martin Aspect Corrected Upper": mfile.get(
            "l_h_threshold_powers(20)", scan=scan
        ),
        "Martin Aspect Corrected Lower": mfile.get(
            "l_h_threshold_powers(21)", scan=scan
        ),
    }
    data_values = list(data.values())
    # Create the violin plot
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data_values)))
    generator = np.random.default_rng(seed=u_seed)
    x_values = generator.normal(loc=1, scale=0.01, size=len(data_values))
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

    axis.set_title("L-H Threshold ($P_\\text{LH}$) Comparison")
    axis.set_ylabel("L-H threshold power [MW]")
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])

    # Add background color
    axis.set_facecolor("#f0f0f0")


def plot_lower_vertical_build(
    axis: plt.Axes, mfile: MFile, colour_scheme: Literal[1, 2]
):
    """Plots the lower vertical build of a fusion device on the given matplotlib axis.

    This function visualizes the different layers/components of the machine's vertical
    build
    (such as plasma, first wall, divertor, shield, vacuum vessel, thermal shield, TF
    coil, etc.)
    as a vertical stacked bar chart. The thickness of each layer is extracted from the
    provided `mfile`, and each segment is color-coded and labeled accordingly.

    Parameters
    ----------
    axis :
        The matplotlib axis on which to plot the vertical build.
    mfile :
        An object containing the machine build data, with required fields for each
        vertical component.
    colour_scheme :
        Colour scheme index to use for component colors.


    Notes
    -----
    This function modifies the provided axis in-place and does not return a value.
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

    lower_vertical_build = [[mfile.get(rl, scan=-1) for rl in lower_vertical_variables]]

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
        (
            TFC_COLOUR[colour_scheme - 1]
            if mfile.get("i_tf_sup", scan=-1) != 0
            else "#b87333"
        ),
        "white",
    ]

    # Remove build parts equal to zero
    mask = ~(lower_vertical_build[:, 0] == 0.0)  # noqa: RUF069
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
            label=(
                f"{filtered_labels[kk]}\n[{lower_vertical_variables[kk]}]\n{filtered_vertical_build[kk][0]:.3f} m"  # noqa: E501
            ),
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


def plot_density_limit_comparison(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot a scatter box plot of different density limit comparisons.

    Parameters
    ----------
    axis :
        Axis object to plot to.
    mfile :
         MFILE data object.
    scan :
        Scan number to use.
    """
    old_asdex = mfile.get("nd_plasma_electron_max_array(1)", scan=scan)
    borrass_iter_i = mfile.get("nd_plasma_electron_max_array(2)", scan=scan)
    borrass_iter_ii = mfile.get("nd_plasma_electron_max_array(3)", scan=scan)
    jet_edge_radiation = mfile.get("nd_plasma_electron_max_array(4)", scan=scan)
    jet_simplified = mfile.get("nd_plasma_electron_max_array(5)", scan=scan)
    hugill_murakami = mfile.get("nd_plasma_electron_max_array(6)", scan=scan)
    greenwald = mfile.get("nd_plasma_electron_max_array(7)", scan=scan)
    asdex_new = mfile.get("nd_plasma_electron_max_array(8)", scan=scan)

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
    data_values = list(data.values())

    # Create the violin plot
    axis.violinplot(data_values, showextrema=False)

    # Create the box plot
    axis.boxplot(data_values, showfliers=True, showmeans=True, meanline=True, widths=0.3)

    # Scatter plot for each data point
    colors = plt.cm.plasma(np.linspace(0, 1, len(data.values())))
    for index, (key, value) in enumerate(data.items()):
        axis.scatter(1, value, color=colors[index], label=key, alpha=1.0)
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Calculate average, standard deviation, and median
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
    axis.set_xlim(0.5, 1.5)
    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.set_facecolor("#f0f0f0")


def plot_iteration_variables(axis: plt.Axes, m_file: MFile, scan: int):
    """Plot the iteration variables and where they lay in their bounds on a given axes

    Parameters
    ----------
    axis: plt.Axes :

    m_file: MFile :

    scan: int :

    """
    # Get total number of iteration variables
    n_itvars = int(m_file.get("n_iteration_variables", scan=scan))

    y_labels = []
    y_pos = []
    n_plot = 0

    # Build a mapping from itvar index to its name (description)
    itvar_names = {}
    for var in m_file.data:
        if var.startswith("itvar"):
            idx = int(var[5:])  # e.g. "itvar001" -> 1
            itvar_names[idx] = m_file.data[var].var_description

    for n_plot, n in enumerate(range(1, n_itvars + 1)):
        # Get the final value of the iteration variable, its bounds, and relative change
        itvar_final = m_file.get(f"itvar{n:03d}", scan=scan)
        itvar_upper = m_file.get(f"boundu{n:03d}", scan=scan)
        itvar_lower = m_file.get(f"boundl{n:03d}", scan=scan)
        itvar_relative_change = m_file.get(f"xcm{n:03d}", scan=scan)
        final_value_normalised = m_file.get(f"nitvar{n:03d}", scan=scan)

        # Use the variable name if available, else fallback to "itvarXXX"
        var_label = itvar_names.get(n, f"itvar{n:03d}")

        norm_relative_change = (
            ((itvar_final / itvar_relative_change) - itvar_lower)
            / (itvar_upper - itvar_lower)
            if itvar_final != itvar_lower
            else 0
        )

        # Plot square marker at the final value if at bounds
        if np.isclose(final_value_normalised, 1.0, atol=1e-3):
            axis.plot(
                1,
                n_plot,
                "s",
                color="black",
                markersize=8,
                label="Lower Bound" if n_plot == 0 else "",
            )
        elif np.isclose(final_value_normalised, 0.0, atol=1e-3):
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
                final_value_normalised,
                left=0,
                height=1.0,
                color="blue",
                edgecolor="black",
                linewidth=1.5,
                alpha=0.7,
                label="Final Value" if n_plot == 0 else "",
            )

        # Plot scatter point for normalised relative change
        axis.scatter(
            norm_relative_change,
            n_plot,
            color="black",
            marker="o",
            linewidths=2,
            alpha=1.0,
            label="Initial Value" if n_plot == 0 else "",
        )

        # Draw an arrow from the initial value to the final value
        axis.annotate(
            "",
            xy=(final_value_normalised, n_plot),
            xytext=(norm_relative_change, n_plot),
            arrowprops={
                "arrowstyle": "->",
                "color": "black",
                "linestyle": "--",
                "linewidth": 1.0,
                "alpha": 0.9,
            },
        )
        # Plot the value as a number at x = 0.5
        axis.text(
            0.5,
            n_plot,
            f"{itvar_final:,.8g}",
            va="center",
            ha="center",
            fontsize=10,
            color=(
                "orange"
                if np.isclose(final_value_normalised, 1.0, atol=1e-3)
                or np.isclose(final_value_normalised, 0.0, atol=1e-3)
                else "green"
            ),
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.8,
                "edgecolor": "white",
                "linewidth": 1,
            },
        )

        # Plot the value of the upper bound to the right of x=1
        axis.text(
            1.05,
            n_plot,
            f"{itvar_upper:,.3g}",
            va="center",
            ha="left",
            fontsize=10,
            color="gray",
        )
        # Plot the value of the lower bound to the left of x=0
        axis.text(
            -0.05,
            n_plot,
            f"{itvar_lower:,.3g}",
            va="center",
            ha="right",
            fontsize=10,
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
    axis.set_facecolor("#f5f5f5")
    axis.set_xlim(-0.2, 1.2)  # Normalised bounds
    axis.set_title("Iteration Variables Bounds")
    axis.set_xticks(np.arange(0, 1.0, 0.1))
    axis.grid(True, axis="x", linestyle="--", alpha=0.3)
    axis.legend(loc="upper left", bbox_to_anchor=(-0.15, 1.05), ncol=1)


def plot_fw_90_deg_pipe_bend(ax, m_file, scan: int):
    """Plot the first wall pipe 90 degree bend on the given axis, with axes in mm.

    Parameters
    ----------
    ax :

    m_file :

    scan: int :

    """
    # Get pipe radius from m_file, fallback to 0.1 m
    r = m_file.get("radius_fw_channel", scan=scan)
    elbow_radius = m_file.get("radius_fw_channel_90_bend", scan=scan)

    draw_bend(
        ax,
        elbow_radius,
        np.pi / 2,
        r,
        title="First Wall Pipe 90° Bend",
        alpha=1.0,
    )


def plot_ebw_ecrh_coupling_graph(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot EBW and ECRH coupling efficiency graph"""
    ebw = ElectronBernstein(plasma_profile=0)
    ecrg = ElectronCyclotron(plasma_profile=0)
    b_on_axis = mfile.get("b_plasma_toroidal_on_axis", scan=scan)
    bs = np.linspace(0.0, b_on_axis + 2.0, 500)
    # Use a color map for harmonics
    colors = ["red", "green", "blue"]
    linestyles = ["-", "--"]  # EBW: solid, ECRH: dashed

    for idx, n_harmonic in enumerate(range(1, 4)):
        eta_ebw_vals = []
        # For ECRH, store results for both wave modes (0: O-mode, 1: X-mode)
        eta_ecrh_vals_omode = []
        eta_ecrh_vals_xmode = []
        for b in bs:
            eta_ebw = ebw.electron_bernstein_freethy(
                te=mfile.get("temp_plasma_electron_vol_avg_kev", scan=scan),
                rmajor=mfile.get("rmajor", scan=scan),
                dene20=mfile.get("nd_plasma_electrons_vol_avg", scan=scan) / 1e20,
                b_plasma_toroidal_on_axis=b,
                n_ecrh_harmonic=n_harmonic,
                xi_ebw=mfile.get("xi_ebw", scan=scan),
            )
            eta_ecrh_omode = ecrg.electron_cyclotron_freethy(
                te=mfile.get("temp_plasma_electron_vol_avg_kev", scan=scan),
                zeff=mfile.get("n_charge_plasma_effective_vol_avg", scan=scan),
                rmajor=mfile.get("rmajor", scan=scan),
                nd_plasma_electrons_vol_avg=mfile.get(
                    "nd_plasma_electrons_vol_avg", scan=scan
                ),
                b_plasma_toroidal_on_axis=b,
                n_ecrh_harmonic=n_harmonic,
                i_ecrh_wave_mode=0,  # O-mode
            )
            eta_ecrh_xmode = ecrg.electron_cyclotron_freethy(
                te=mfile.get("temp_plasma_electron_vol_avg_kev", scan=scan),
                zeff=mfile.get("n_charge_plasma_effective_vol_avg", scan=scan),
                rmajor=mfile.get("rmajor", scan=scan),
                nd_plasma_electrons_vol_avg=mfile.get(
                    "nd_plasma_electrons_vol_avg", scan=scan
                ),
                b_plasma_toroidal_on_axis=b,
                n_ecrh_harmonic=n_harmonic,
                i_ecrh_wave_mode=1,  # X-mode
            )
            eta_ebw_vals.append(eta_ebw)
            eta_ecrh_vals_omode.append(eta_ecrh_omode)
            eta_ecrh_vals_xmode.append(eta_ecrh_xmode)
        # EBW: solid, ECRH O-mode: dashed, ECRH X-mode: dotted, same color for same
        # harmonic
        axis.plot(
            bs,
            eta_ebw_vals,
            label=f"EBW (harmonic {n_harmonic})",
            color=colors[idx],
            linestyle=linestyles[0],
        )
        axis.plot(
            bs,
            eta_ecrh_vals_omode,
            label=f"ECRH O-mode (harmonic {n_harmonic})",
            color=colors[idx],
            linestyle="--",
        )
        axis.plot(
            bs,
            eta_ecrh_vals_xmode,
            label=f"ECRH X-mode (harmonic {n_harmonic})",
            color=colors[idx],
            linestyle=":",
        )
    axis.set_xlabel("On axis toroidal B-field [T]")
    axis.set_ylabel("Current drive efficiency [A/W]")
    axis.set_title("EBW/ECRH Coupling Efficiency vs Toroidal B-field")
    axis.legend()
    axis.grid(True)
    # Plot a vertical line at the on-axis value of the toroidal B-field
    b_on_axis = mfile.get("b_plasma_toroidal_on_axis", scan=scan)
    axis.axvline(
        b_on_axis,
        color="black",
        linestyle="-",
        linewidth=2.5,
        label="On-axis $B_T$",
    )
    axis.minorticks_on()


__all__ = [
    "RadialBuild",
    "plot_centre_cross",
    "plot_density_limit_comparison",
    "plot_ebw_ecrh_coupling_graph",
    "plot_fw_90_deg_pipe_bend",
    "plot_h_threshold_comparison",
    "plot_iteration_variables",
    "plot_lower_vertical_build",
]
