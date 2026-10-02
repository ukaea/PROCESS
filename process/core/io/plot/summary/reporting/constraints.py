"""Reporting functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.rendering import draw_text

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_equality_constraint_equations(axis: plt.Axes, m_file_data: MFile, scan: int):
    """Plot the equality constraints for a solution and their normalised residuals

    Parameters
    ----------
    axis: plt.Axes :

    m_file_data: MFile :

    scan: int :

    """
    y_labels = []
    y_pos = []

    # Build a mapping from itvar index to its name (description)
    con_names = {}
    con_numbers = {}
    for var in m_file_data.data:
        if var.startswith("eq_con"):
            idx = int(var[6:])  # e.g. "itvar001" -> 1
            con_names[idx] = m_file_data.data[var].var_description
            con_numbers[idx] = idx

    for n_plot, n in enumerate(con_numbers.values()):
        # Constraint value needed
        con_value = m_file_data.data[f"val_eq_con{n:03d}"].get_scan(scan)

        # Use the variable name if available, else fallback to "eq_conXXX"
        var_label = con_names.get(n, f"eq_con{n:03d}")

        # Normalized residual of the constraint
        con_norm_residual = m_file_data.data[f"eq_con{n:03d}"].get_scan(scan)

        # Unit type of the constraint
        con_units_raw = m_file_data.data[f"eq_units_con{n:03d}"].get_scan(scan)
        con_units = str(con_units_raw).strip("'`")

        # Remove '_normalised_residue' from the label if present
        if isinstance(var_label, str) and var_label.endswith("_normalised_residue"):
            var_label = var_label.replace("_normalised_residue", "")

            # Remove trailing underscores and replace underscores between words with
            # spaces
            var_label = var_label.rstrip("_").replace("_", " ")

        # Plot the normalised residual as a bar
        axis.barh(
            n_plot,
            con_norm_residual,
            height=0.6,
            color="blue",
            label="Normalized Residual" if n_plot == 0 else "",
            align="center",
        )

        # Add the value as a number to the right of the bar
        draw_text(
            axis,
            con_norm_residual + 0.52,
            n_plot,
            f"{con_norm_residual:.8g}",
            va="center",
            ha="left",
            fontsize=8,
            color="blue",
        )

        # Add the constraint value as text to the left of the y-axis
        draw_text(
            axis,
            0.45,
            n_plot,
            f"{con_value:.8g} {con_units}",
            va="center",
            ha="right",
            fontsize=8,
            color="black",
        )

        y_labels.append(var_label)
        y_pos.append(n_plot)

    axis.axvline(0.5, color="red", linewidth=2, zorder=0)
    axis.set_yticks(y_pos)
    axis.set_yticklabels(y_labels)
    axis.set_facecolor("#f5f5f5")
    axis.set_xlim(-0.4, 1.2)  # Normalised bounds
    axis.set_title("Equality Constraint Equations")
    axis.set_xticks([])
    axis.legend()


def plot_inequality_constraint_equations(axis: plt.Axes, m_file: MFile, scan: int):
    """Plot the inequality constraints for a solution and where they lay within their
    bounds

    Parameters
    ----------
    axis: plt.Axes :

    m_file: MFile :

    scan: int :

    """
    y_labels = []
    y_pos = []

    # Build a mapping from itvar index to its name (description)
    con_names = {}
    con_numbers = {}
    for var in m_file.data:
        if var.startswith("ineq_con"):
            idx = int(var[8:])  # e.g. "ineq_con001" -> 1
            con_names[idx] = m_file.data[var].var_description
            con_numbers[idx] = idx

    for n_plot, n in enumerate(con_numbers.values()):
        # Constraint value/bound
        con_bound = m_file.data[f"ineq_bound_con{n:03d}"].get_scan(scan)

        # Value of constraint variable
        con_value = m_file.data[f"ineq_value_con{n:03d}"].get_scan(scan)

        # Constraint symbol can be `<=` for an upper limit or `>=` for a lower limit
        con_symbol = m_file.data[f"ineq_symbol_con{n:03d}"].get_scan(scan)

        # Use the variable name if available, else fallback to "ineq_conXXX"
        var_label = con_names.get(n, f"ineq_con{n:03d}")

        # Normalized residual of the constraint
        con_residual_norm = m_file.data[f"ineq_con{n:03d}"].get_scan(scan)

        # Unit type of the constraint
        con_units = m_file.data[f"ineq_units_con{n:03d}"].get_scan(scan).strip("'`")

        # Add a vertical line at the normalised constraint bounds of 0 and 1
        axis.axvline(
            0.0,
            color="red",
            linestyle="--",
            linewidth=1.5,
            zorder=0,
        )

        axis.axvline(
            1.0,
            color="red",
            linestyle="--",
            linewidth=1.5,
            zorder=0,
        )

        # Remove '_normalised_residue' from the label if present
        if isinstance(var_label, str) and var_label.endswith("_normalised_residue"):
            var_label = var_label.replace("_normalised_residue", "")
            var_label = var_label.rstrip("_").replace("_", " ")

        # Calculate the normalised constraint threshold depending if the constraint is an
        # upper
        # or lower limit
        if con_symbol == "'<='":
            normalised_value = 1 - con_residual_norm
            bar_left = normalised_value
            bar_width = 1 - normalised_value
        else:
            # For a lower limit, the normalised value is the residual itself
            normalised_value = con_residual_norm
            bar_left = 0
            # Set the bar width to be 1/10 times the normalised value,
            # but cap it at 1.0 to avoid overly long bars
            bar_width = min(normalised_value * 0.1, 1.0)

        # If the constraint value is very close to the bound then plot a square marker at
        # the bound
        if np.isclose(normalised_value, 1.0, atol=1e-3):
            axis.plot(
                1,
                n_plot,
                "s",
                color="black",
                markersize=8,
                zorder=5,
            )
        elif np.isclose(normalised_value, 0.0, atol=1e-3):
            axis.plot(
                0,
                n_plot,
                "s",
                color="black",
                markersize=8,
                zorder=5,
            )

        else:
            # If constraint value is not very close to bound then plot bar as normal
            axis.barh(
                n_plot,
                bar_width,
                left=bar_left,
                color="blue",
                edgecolor="black",
                linewidth=1.5,
                height=1.0,
                alpha=0.7,
                label="Constraint Value" if n_plot == 0 else "",
            )

        # Plot the value as a number at x = 0.5
        draw_text(
            axis,
            0.5,
            n_plot,
            f"{con_value:,.8g} {con_units}",
            va="center",
            ha="center",
            fontsize=8,
            color=(
                "orange"
                if np.isclose(normalised_value, 1.0, atol=1e-3)
                or np.isclose(normalised_value, 0.0, atol=1e-3)
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
        # Annoate the bound value depending if it is an upper or lower limit
        if con_symbol == "'<='":
            # Add the constraint symbol and bound as text
            draw_text(
                axis,
                1.02,  # Position text slightly to the right of the normalised bound
                n_plot,
                f"$\\leq$ {con_bound:,.8g} {con_units}",
                va="center",
                ha="left",
                fontsize=8,
                color="black",
            )
        else:  # con_symbol == ">="
            draw_text(
                axis,
                -0.025,  # Position text slightly to the left of the normalised bound
                n_plot,
                f"$\\geq$ {con_bound:,.8g} {con_units}",
                va="center",
                ha="right",
                fontsize=8,
                color="black",
            )

        y_labels.append(var_label)
        y_pos.append(n_plot)

    axis.set_yticks(y_pos)
    axis.set_yticklabels(y_labels)
    axis.set_title("Inequality Constraint Equations")
    axis.set_xlim(-0.3, 1.275)
    axis.set_xticks([])
    axis.set_facecolor("#f5f5f5")
    axis.set_xticks(np.arange(0, 1.0, 0.1))
    axis.grid(True, axis="x", linestyle="--", alpha=0.3)
    axis.set_xticklabels([])


__all__ = [
    "plot_equality_constraint_equations",
    "plot_inequality_constraint_equations",
]
