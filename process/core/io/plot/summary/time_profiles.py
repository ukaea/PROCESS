"""Time Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.common import (
    box_style,
    get_pulse_timings,
)
from process.core.io.plot.summary.magnets import (
    secs_to_hms,
)
from process.core.io.plot.summary.rendering import (
    draw_text,
)

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_current_profiles_over_time(axis: plt.Axes, mfile: MFile, scan: int):
    """Plots the current profiles over time for PF circuits, CS coil, and plasma."""
    pulse_timings = get_pulse_timings(mfile, scan)

    # Find the number of PF circuits, n_pf_cs_plasma_circuits includes the CS and plasma
    # circuits
    n_pf_cs_plasma_circuits = mfile.get("n_pf_cs_plasma_circuits", scan=scan)

    # Extract PF circuit times
    # n_pf_cs_plasma_circuits contains the CS and plasma at the end so we subtract 2
    for i in range(int(n_pf_cs_plasma_circuits - 2)):
        circuit_current = [
            mfile.get(f"pfc{i}t{j}", scan=scan)
            for j in range(pulse_timings.n_pf_active_points_total)
        ]
        # Change from 0 to 1 index to align with poloidal cross-section plot numbering
        axis.plot(
            pulse_timings.pf_active_cumulative,
            circuit_current,
            label=f"PF Coil {i + 1}",
            linestyle="--",
        )

    # Since CS may not always be present try to retrieve values
    try:
        cs_circuit = [
            mfile.get(f"cs_t{i}", scan=scan)
            for i in range(pulse_timings.n_pf_active_points_total)
        ]
        axis.plot(
            pulse_timings.pf_active_cumulative,
            cs_circuit,
            label="CS Coil",
            linestyle="--",
        )
    except KeyError:
        pass

    # Plasma current values
    plasmat1 = mfile.get("plasmat1", scan=scan)
    plasmat2 = mfile.get("plasmat2", scan=scan)
    plasmat3 = mfile.get("plasmat3", scan=scan)
    plasmat4 = mfile.get("plasmat4", scan=scan)
    plasmat5 = mfile.get("plasmat5", scan=scan)

    # x-coordinates for the plasma current
    x_plasma = pulse_timings.pf_active_cumulative[1:]
    # x-coordinates for the plasma current
    y_plasma = [plasmat1, plasmat2, plasmat3, plasmat4, plasmat5]

    # Plot the plasma current
    axis.plot(x_plasma, y_plasma, "black", linewidth=2, label="Plasma")

    # Move the x-axis to 0 on the y-axis
    axis.spines["bottom"].set_position("zero")

    # Annotate key points
    # Create a secondary x-axis for annotations
    secax = axis.secondary_xaxis("bottom")
    # Exclude the dwell point so tick positions and labels remain aligned.
    secax.set_xticks(pulse_timings.pf_active_cumulative[:-1])
    secax.set_xticklabels(
        pulse_timings.POINT_LABELS[
            :-1
        ],  # Exclude the last label as it corresponds to the dwell period
        rotation=60,
    )
    secax.tick_params(axis="x", which="major")

    # Add axis labels
    axis.set_xlabel("Time [s]", fontsize=12)
    axis.xaxis.set_label_coords(1.05, 0.5)
    axis.set_ylabel("Current [A]", fontsize=12)

    # Add a title
    axis.set_title("Current Profiles Over Time", fontsize=14)

    # Add a legend
    axis.legend()

    axis.set_yscale("symlog")

    # Add a grid for better readability
    axis.grid(True, linestyle="--", alpha=0.6)


def plot_system_power_profiles_over_time(axis: plt.Axes, mfile: MFile, scan: int, fig):
    """Plots the power profiles over time for various systems."""
    pulse_timings = get_pulse_timings(mfile, scan)

    # Create empty arrays for the power at each time step for each system
    power_profiles = {
        "Fusion Power": np.zeros(pulse_timings.n_pulse_points_total),
        "Plant Base Load": np.zeros(pulse_timings.n_pulse_points_total),
        "Cryo Plant": np.zeros(pulse_timings.n_pulse_points_total),
        "Tritium Plant": np.zeros(pulse_timings.n_pulse_points_total),
        "Vacuum Pumps": np.zeros(pulse_timings.n_pulse_points_total),
        "TF Coil Supplies": np.zeros(pulse_timings.n_pulse_points_total),
        "PF Coil Supplies": np.zeros(pulse_timings.n_pulse_points_total),
        "Coolant Pump Elec Total": np.zeros(pulse_timings.n_pulse_points_total),
        "HCD Electric Total": np.zeros(pulse_timings.n_pulse_points_total),
        "Gross Electric Power": np.zeros(pulse_timings.n_pulse_points_total),
        "Net Electric Power": np.zeros(pulse_timings.n_pulse_points_total),
    }

    # Fill power_profiles arrays using vectorized assignment
    for label, key in [
        ("Fusion Power", "p_fusion_total_profile_mw"),
        ("Gross Electric Power", "p_plant_electric_gross_profile_mw"),
        ("Net Electric Power", "p_plant_electric_net_profile_mw"),
        ("Plant Base Load", "p_plant_electric_base_total_profile_mw"),
        ("Cryo Plant", "p_cryo_plant_electric_profile_mw"),
        ("Tritium Plant", "p_tritium_plant_electric_profile_mw"),
        ("Vacuum Pumps", "vachtmw_profile_mw"),
        ("TF Coil Supplies", "p_tf_electric_supplies_profile_mw"),
        ("PF Coil Supplies", "p_pf_electric_supplies_profile_mw"),
        ("Coolant Pump Elec Total", "p_coolant_pump_elec_total_profile_mw"),
        ("HCD Electric Total", "p_hcd_electric_total_profile_mw"),
    ]:
        for time in range(pulse_timings.n_pulse_points_total):
            power_profiles[label][time] = mfile.get(f"{key}{time}", scan=scan)

    # Define line styles for each system
    # All net drains (negative power flows) use the same line style: dashed
    line_styles = {
        "Fusion Power": ":",
        "Plant Base Load": "--",
        "Cryo Plant": "--",
        "Tritium Plant": "--",
        "Vacuum Pumps": "--",
        "TF Coil Supplies": "--",
        "PF Coil Supplies": "--",
        "Coolant Pump Elec Total": "--",
        "HCD Electric Total": "--",
        "Gross Electric Power": "-",
        "Net Electric Power": "-",
    }

    # Plot each system's power profile over time with different line styles
    for label, powers in power_profiles.items():
        style = line_styles.get(label, "-")
        axis.plot(
            pulse_timings.total_pulse_cumulative,
            powers,
            label=label,
            linestyle=style,
        )

    # Move the x-axis to 0 on the y-axis
    axis.spines["bottom"].set_position("zero")

    # Annotate key points
    # Create a secondary x-axis for annotations
    secax = axis.secondary_xaxis("bottom")
    # Label phase starts only (exclude final end-of-dwell point).
    secax.set_xticks(pulse_timings.total_pulse_cumulative[:-1])
    secax.set_xticklabels(
        pulse_timings.POINT_LABELS,
        rotation=60,
    )
    secax.tick_params(axis="x", which="major")

    # Add axis labels
    axis.set_xlabel("Time [s]", fontsize=12)
    axis.xaxis.set_label_coords(1.05, 0.5)
    axis.set_ylabel("Power [MW]", fontsize=12)

    # Add a title
    axis.set_title("System Power Over Time", fontsize=14)

    # Add a legend
    axis.legend()

    axis.set_yscale("symlog")
    axis.minorticks_on()
    axis.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)

    # Add a grid for better readability
    axis.grid(True, linestyle="--", alpha=0.6)

    # Add energy produced info
    textstr_energy = (
        "$\\mathbf{Energy \\ Production:}$\n\nEnergy produced over whole"
        f" pulse: {mfile.get('e_plant_net_electric_pulse_mj', scan=scan):,.4f}"
        " MJ\nEnergy produced over whole pulse:"
        f" {mfile.get('e_plant_net_electric_pulse_kwh', scan=scan):,.4f} kWh\n"
    )

    draw_text(
        axis,
        0.075,
        0.2,
        textstr_energy,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("grey"),
    )

    # Add energy produced info

    textstr_times = (
        "$\\mathbf{Pulse \\ Timings:}$\n\nCoil precharge,"
        " $t_{\\text{precharge}}$:       "
        f" {mfile.get('t_plant_pulse_coil_precharge', scan=scan):,.1f} s "
        f" ({secs_to_hms(mfile.get('t_plant_pulse_coil_precharge', scan=scan))})\nCurrent"  # noqa: E501
        " ramp up, $t_{\\text{current ramp}}$: "
        f" {mfile.get('t_plant_pulse_plasma_current_ramp_up', scan=scan):,.1f}"
        f" s  ({secs_to_hms(mfile.get('t_plant_pulse_plasma_current_ramp_up', scan=scan))})\nFusion"  # noqa: E501
        " ramp, $t_{\\text{fusion ramp}}$:         "
        f" {mfile.get('t_plant_pulse_fusion_ramp', scan=scan):,.1f} s "
        f" ({secs_to_hms(mfile.get('t_plant_pulse_fusion_ramp', scan=scan))})\nBurn,"
        " $t_{\\text{burn}}$:                             "
        f" {mfile.get('t_plant_pulse_burn', scan=scan):,.1f} s "
        f" ({secs_to_hms(mfile.get('t_plant_pulse_burn', scan=scan))})\nRamp"
        " down, $t_{\\text{ramp down}}$:          "
        f" {mfile.get('t_plant_pulse_plasma_current_ramp_down', scan=scan):,.1f}"
        f" s  ({secs_to_hms(mfile.get('t_plant_pulse_plasma_current_ramp_down', scan=scan))})\nBetween"  # noqa: E501
        " pulse, $t_{\\text{between pulse}}$:  "
        f" {mfile.get('t_plant_pulse_dwell', scan=scan):,.1f} s "
        f" ({secs_to_hms(mfile.get('t_plant_pulse_dwell', scan=scan))})\n\nTotal"
        " pulse length, $t_{\\text{cycle}}$:       "
        f" {mfile.get('t_plant_pulse_total', scan=scan):,.1f} s "
        f" ({secs_to_hms(mfile.get('t_plant_pulse_total', scan=scan))})\n"
    )

    draw_text(
        axis,
        0.6,
        0.225,
        textstr_times,
        fontsize=9,
        verticalalignment="top",
        transform=fig.transFigure,
        bbox=box_style("grey"),
    )


__all__ = [
    "plot_current_profiles_over_time",
    "plot_system_power_profiles_over_time",
]
