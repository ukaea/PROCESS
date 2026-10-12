"""Reporting functions for PROCESS summary plots."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile

__all__ = ["RadialBuild", "plot_centre_cross", "plot_iteration_variables"]


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
