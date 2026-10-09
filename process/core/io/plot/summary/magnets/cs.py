"""Magnets functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

from process.core.io.plot.summary.common import box_style, setup_axis, text_layout
from process.core.io.plot.summary.constants import SOLENOID_COLOUR
from process.core.io.plot.summary.reporting.text import plot_info
from process.data_structure.pfcoil_variables import NFIXMX
from process.models.superconductors import SuperconductorModel

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def secs_to_hms(s):
    """Convert seconds to 'Hh Mm Ss' string."""
    s = float(s)
    return f"{int(s // 3600)}h {int((s % 3600) // 60)}m {int(s % 60)}s"


def plot_physics_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot geometry info

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    axis.text(-0.05, 1, "Physics:", ha="left", va="center")
    setup_axis(axis, xmin=0, xmax=1, ymin=-16, ymax=1)

    nong = mfile.get("nd_plasma_electron_line", scan=scan) / mfile.get(
        "nd_plasma_electron_max_array(7)", scan=scan
    )

    nd_plasma_impurities_vol_avg = mfile.get(
        "nd_plasma_impurities_vol_avg", scan=scan
    ) / mfile.get("nd_plasma_electrons_vol_avg", scan=scan)

    tepeak = mfile.get("temp_plasma_electron_on_axis_kev", scan=scan) / mfile.get(
        "temp_plasma_electron_vol_avg_kev", scan=scan
    )

    nepeak = mfile.get("nd_plasma_electron_on_axis", scan=scan) / mfile.get(
        "nd_plasma_electrons_vol_avg", scan=scan
    )

    # Assume Martin scaling if pthresh is not printed
    # Accounts for pthresh not being written prior to issue #679 and #680
    if "p_l_h_threshold_mw" in mfile.data:
        pthresh = mfile.get("p_l_h_threshold_mw", scan=scan)
    else:
        pthresh = mfile.get("l_h_threshold_powers(6)", scan=scan)

    data = [
        ("p_fusion_total_mw", "Fusion power", "MW"),
        ("big_q_plasma", "$Q_{p}$", ""),
        ("plasma_current_ma", "$I_p$", "MA"),
        ("b_plasma_toroidal_on_axis", "Vacuum $B_T$ at $R_0$", "T"),
        ("q95", r"$q_{\mathrm{95}}$", ""),
        ("beta_norm_thermal", r"$\beta_N$, thermal", "% m T MA$^{-1}$"),
        ("beta_norm_toroidal", r"$\beta_N$, toroidal", "% m T MA$^{-1}$"),
        ("beta_thermal_poloidal_vol_avg", r"$\beta_P$, thermal", ""),
        ("beta_poloidal_vol_avg", r"$\beta_P$, total", ""),
        ("temp_plasma_electron_vol_avg_kev", r"$\langle T_e \rangle$", "keV"),
        ("nd_plasma_electrons_vol_avg", r"$\langle n_e \rangle$", "m$^{-3}$"),
        (nong, r"$\langle n_{\mathrm{e,line}} \rangle \ / \ n_G$", ""),
        (tepeak, r"$T_{e0} \ / \ \langle T_e \rangle$", ""),
        (nepeak, r"$n_{e0} \ / \ \langle n_{\mathrm{e, vol}} \rangle$", ""),
        ("n_charge_plasma_effective_vol_avg", r"$Z_{\mathrm{eff}}$", ""),
        (
            nd_plasma_impurities_vol_avg,
            r"$n_Z \ / \  \langle n_{\mathrm{e, vol}} \rangle$",
            "",
        ),
        ("t_energy_confinement", r"$\tau_e$", "s"),
        ("hfact", "H-factor", ""),
        (pthresh, "H-mode threshold", "MW"),
        ("tauelaw", "Scaling law", ""),
    ]

    plot_info(axis, data, mfile, scan)


def plot_magnetics_info(axis: plt.Axes, mfile: MFile, scan: int):
    """Function to plot magnet info

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    """
    # Check for Copper magnets
    i_tf_sup = int(mfile.get("i_tf_sup", scan=scan)) if "i_tf_sup" in mfile.data else 1

    axis.text(-0.05, 1, "Coil currents etc:", ha="left", va="center")
    setup_axis(axis, xmin=0, xmax=1, ymin=-16, ymax=1)

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile.data:
        if "r_pf_coil_middle[" in item:
            number_of_coils += 1

    pf_info = [
        (
            mfile.get(f"c_pf_cs_coils_peak_ma[{i:01}]", scan=scan),
            f"PF {i}",
        )
        for i in range(1, number_of_coils)
        if i % 2 != 0
    ]

    if len(pf_info) > 2:
        pf_info_3_a = pf_info[2][0]
        pf_info_3_b = pf_info[2][1]
    else:
        pf_info_3_a = ""
        pf_info_3_b = ""

    t_plant_pulse_burn = mfile.get("t_plant_pulse_burn", scan=scan) / 3600.0

    i_tf_bucking = (
        int(mfile.get("i_tf_bucking", scan=scan)) if "i_tf_bucking" in mfile.data else 1
    )

    # Get superconductor material (i_tf_sc_mat)
    # If i_tf_sc_mat not present, assume resistive
    i_tf_sc_mat = (
        int(mfile.get("i_tf_sc_mat", scan=scan)) if "i_tf_sc_mat" in mfile.data else 0
    )

    tftype = (
        SuperconductorModel(int(mfile.get("i_tf_sc_mat", scan=scan))).full_name
        if i_tf_sc_mat > 0
        else "Resistive Copper"
    )

    vssoft = mfile.get("vs_plasma_res_ramp", scan=scan) + mfile.get(
        "vs_plasma_ind_ramp", scan=scan
    )

    sig_case = 1.0e-6 * mfile.get(f"s_shear_tf_peak({i_tf_bucking})", scan=scan)
    sig_cond = 1.0e-6 * mfile.get(f"s_shear_tf_peak({i_tf_bucking + 1})", scan=scan)

    if i_tf_sup == 1:
        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vs_cs_pf_total_pulse", "Available flux swing", "Wb"),
            (t_plant_pulse_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            (
                "b_tf_inboard_peak_with_ripple",
                "Peak field at conductor (w. rip.)",
                "T",
            ),
            ("f_c_tf_turn_operating_critical", r"I/I$_{\mathrm{crit}}$", ""),
            ("temp_tf_superconductor_margin", "TF Temperature margin", "K"),
            ("temp_cs_superconductor_margin", "CS Temperature margin", "K"),
            (sig_cond, "TF Cond max TRESCA stress", "MPa"),
            (sig_case, "TF Case max TRESCA stress", "MPa"),
            ("m_tf_coils_total/n_tf_coils", "Mass per TF coil", "kg"),
        ]

    else:
        p_cp_resistive = 1.0e-6 * mfile.get("p_cp_resistive", scan=scan)
        p_tf_leg_resistive = 1.0e-6 * mfile.get("p_tf_leg_resistive", scan=scan)
        p_tf_joints_resistive = 1.0e-6 * mfile.get("p_tf_joints_resistive", scan=scan)
        fcoolcp = 100.0 * mfile.get("fcoolcp", scan=scan)

        data = [
            (pf_info[0][0], pf_info[0][1], "MA"),
            (pf_info[1][0], pf_info[1][1], "MA"),
            (pf_info_3_a, pf_info_3_b, "MA"),
            (vssoft, "Startup flux swing", "Wb"),
            ("vs_cs_pf_total_pulse", "Available flux swing", "Wb"),
            (t_plant_pulse_burn, "Burn time", "hrs"),
            ("", "", ""),
            (f"#TF coil type is {tftype}", "", ""),
            (
                "b_tf_inboard_peak_symmetric",
                "Peak field at conductor (w. rip.)",
                "T",
            ),
            ("c_tf_total", "TF coil currents sum", "A"),
            ("", "", ""),
            ("#TF coil forces/stresses", "", ""),
            (sig_cond, "TF conductor max TRESCA stress", "MPa"),
            (sig_case, "TF bucking max TRESCA stress", "MPa"),
            (fcoolcp, "CP cooling fraction", "%"),
            (
                "vel_cp_coolant_midplane",
                "Maximum coolant flow speed",
                "ms$^{-1}$",
            ),
            (p_cp_resistive, "CP resistive heating", "MW"),
            (
                p_tf_leg_resistive,
                "legs resistive heating (all legs)",
                "MW",
            ),
            (p_tf_joints_resistive, "TF joints resistive heating ", "MW"),
        ]

    plot_info(axis, data, mfile, scan)


def plot_cs_coil_structure(
    axis: plt.Axes, fig, mfile: MFile, scan: int, colour_scheme=1
):
    """Function to plot the coil structure of the CS.

    Parameters
    ----------
    axis :
        axis object to plot to
    mfile :
        MFILE
    scan :
        scan number to use
    colour_scheme :
        colour scheme to use for the plot (Default value = 1)

    """
    # Get CS coil parameters
    dr_cs = mfile.get("dr_cs", scan=scan)
    dr_cs_full = mfile.get("dr_cs_full", scan=scan)
    dz_cs_full = mfile.get("dz_cs_full", scan=scan)
    dz_cs = mfile.get("dz_cs_full", scan=scan)
    dr_cs_bore = mfile.get("dr_cs_bore", scan=scan)
    r_cs_current_filaments_array = [
        mfile.get(f"r_pf_cs_current_filaments{i}", scan=scan) for i in range(NFIXMX)
    ]
    z_cs_current_filaments_array = [
        mfile.get(f"z_pf_cs_current_filaments{i}", scan=scan) for i in range(NFIXMX)
    ]

    # Plot the right side of the CS
    right_cs = patches.Rectangle(
        (dr_cs_bore, -dz_cs / 2),
        dr_cs,
        dz_cs,
        edgecolor="black",
        facecolor=SOLENOID_COLOUR[colour_scheme - 1],
        lw=1.5,
    )
    axis.add_patch(right_cs)

    # Plot the bore of the machine
    bore_rect = patches.Rectangle(
        (-dr_cs_bore, -dz_cs / 2),
        dr_cs_bore * 2,
        dz_cs,
        edgecolor="black",
        facecolor="lightgrey",
        lw=1.0,
    )
    axis.add_patch(bore_rect)

    left_cs = patches.Rectangle(
        (-dr_cs_bore - dr_cs, -dz_cs / 2),
        dr_cs,
        dz_cs,
        edgecolor="black",
        facecolor=SOLENOID_COLOUR[colour_scheme - 1],
        lw=1.5,
    )
    axis.add_patch(left_cs)

    # Draw vertical lines to represent CS turns
    # Get the turn width (radial thickness of each turn)
    dr_cs_turn = mfile.get("dr_cs_turn", scan=scan)
    dz_cs_turn = mfile.get("dz_cs_turn", scan=scan)
    # Number of vertical lines (number of turns)
    t_kwargs = {"color": "black", "linestyle": "--", "linewidth": 0.2}
    if dr_cs_turn > 0:
        n_lines = int(dr_cs / dr_cs_turn)
        for i in range(1, n_lines):
            x = dr_cs_bore + i * dr_cs_turn
            axis.plot([x, x], [-dz_cs / 2, dz_cs / 2], **t_kwargs)
            x_left = -dr_cs_bore - dr_cs + i * dr_cs_turn
            axis.plot([x_left, x_left], [-dz_cs / 2, dz_cs / 2], **t_kwargs)
    # Plot horizontal lines (along Z) for each turn
    if dz_cs_turn > 0:
        n_hlines = int(dz_cs / dz_cs_turn)
        for j in range(1, n_hlines):
            y = -dz_cs / 2 + j * dz_cs_turn
            # Right CS
            axis.plot([dr_cs_bore, dr_cs_bore + dr_cs], [y, y], **t_kwargs)
            # Left CS
            axis.plot([-dr_cs_bore - dr_cs, -dr_cs_bore], [y, y], **t_kwargs)

        l_kwargs = {
            "color": "black",
            "linestyle": "--",
            "linewidth": 0.6,
            "alpha": 0.5,
        }

        # Plot a horizontal line at y = 0.0
        axis.axhline(y=0.0, **l_kwargs)
        # Plot a vertical line at x = 0.0
        axis.axvline(x=0.0, **l_kwargs)
        # Plot a vertical line at x = dr_cs_bore
        axis.axvline(x=dr_cs_bore, **l_kwargs)
        # Plot a vertical line at x = -dr_cs_bore
        axis.axvline(x=-dr_cs_bore, **l_kwargs)
        # Plot a vertical line at x = dr_cs_bore + dr_cs
        axis.axvline(x=(dr_cs_bore + dr_cs), **l_kwargs)
        # Plot a vertical line at x = -dr_cs_bore - dr_cs
        axis.axvline(x=-(dr_cs_bore + dr_cs), **l_kwargs)
        # Plot a vertical line at y= dz_cs / 2
        axis.axhline(y=(dz_cs / 2), **l_kwargs)
        # Plot a vertical line at y= -dz_cs / 2
        axis.axhline(y=-(dz_cs / 2), **l_kwargs)

        # Plot a vertical line at x = r_cs_middle
        axis.axvline(x=mfile.get("r_cs_middle", scan=scan), **l_kwargs)
        # Plot a vertical line at x= -r_cs_middle
        axis.axvline(x=-mfile.get("r_cs_middle", scan=scan), **l_kwargs)

        # Arrow for coil width
        axis.annotate(
            "",
            xy=(0, (dz_cs_full / 2)),
            xytext=(0, -(dz_cs_full / 2)),
            arrowprops={"arrowstyle": "<->", "color": "black"},
        )

        # Add a label for full coil width
        axis.text(
            0.0,
            -(dz_cs_full / 4),
            f"{dz_cs_full:.3f} m",
            fontsize=7,
            color="black",
            rotation=270,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=box_style("pink", linewidth=0),
        )

        # Arrow for coil width
        axis.annotate(
            "",
            xy=(-(dr_cs_full / 2), (dz_cs_full / 4)),
            xytext=((dr_cs_full / 2), (dz_cs_full / 4)),
            arrowprops={"arrowstyle": "<->", "color": "black"},
        )

        # Add a label for full coil width
        axis.text(
            0.0,
            (dz_cs_full / 4),
            f"{dr_cs_full:.3f} m",
            fontsize=7,
            color="black",
            rotation=0,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=box_style("pink", linewidth=0),
        )

    textstr_cs = (
        "$\\mathbf{Coil \\ parameters:}$\n\nCS height vs TF internal"
        f" height: {mfile.get('f_z_cs_tf_internal', scan=scan):.2f}\nCS"
        f" thickness: {mfile.get('dr_cs', scan=scan):.4f} m\nCS radial middle:"
        f" {mfile.get('r_cs_middle', scan=scan):.4f} m\nCS full height:"
        f" {mfile.get('dz_cs_full', scan=scan):.4f} m\nCS full width:"
        f" {mfile.get('dr_cs_full', scan=scan):.4f} m\nCS poloidal area:"
        f" {mfile.get('a_cs_poloidal', scan=scan):.4f} m$^2$\nCS top-down"
        f" toroidal area: {mfile.get('a_cs_toroidal', scan=scan):.4f}"
        " m$^2$\n$N_{\\text{turns}}:$"
        f" {mfile.get('n_pf_coil_turns[n_cs_pf_coils-1]', scan=scan):,.2f}\n$I_{{\\text{{peak}}}}:$"  # noqa: E501
        f" {mfile.get('c_pf_cs_coils_peak_ma[n_cs_pf_coils-1]', scan=scan):.3f}"
        " MA\n$B_{\\text{peak}}:$"
        f" {mfile.get('b_pf_coil_peak[n_cs_pf_coils-1]', scan=scan):.3f}"
        " T\n$F_{\\text{z,self,peak}}:$"
        f" {mfile.get('forc_z_cs_self_peak_midplane', scan=scan) / 1e6:.3f}"
        " MN\n$\\sigma_{\\text{z,self,peak}}:$"
        f" {mfile.get('stress_z_cs_self_peak_midplane', scan=scan) / 1e6:.3f}"
        " MPa\n$\\sigma_{\\text{mises,peak}}:$"
        f" {mfile.get('stress_mises_cs_peak', scan=scan) / 1e6:.3f}"
        " MPa\n$\\tau_{\\text{shear,peak}}:$"
        f" {mfile.get('stress_shear_cs_peak', scan=scan) / 1e6:.3f} MPa "
    )

    axis.text(
        0.5,
        0.6,
        textstr_cs,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # Plot the current filament points as blue dots and label them

    axis.plot(
        r_cs_current_filaments_array,
        z_cs_current_filaments_array,
        "bo",
        markersize=2,
        label="CS, PF and Plasma Current Filaments",
    )

    axis.set_xlabel("R [m]")
    axis.set_ylabel("Z [m]")
    axis.set_title("Central Solenoid Poloidal Cross-Section")
    axis.grid(True, linestyle="--", alpha=0.3)
    axis.minorticks_on()
    axis.legend()


def plot_cs_turn_structure(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plot the CS turn structure"""
    a_cs_turn = mfile.get("a_cs_turn", scan=scan)
    dz_cs_turn = mfile.get("dz_cs_turn", scan=scan)
    dr_cs_turn = mfile.get("dr_cs_turn", scan=scan)

    f_dr_dz_cs_turn = mfile.get("f_dr_dz_cs_turn", scan=scan)
    radius_cs_turn_cable_space = mfile.get("radius_cs_turn_cable_space", scan=scan)
    dz_cs_turn_conduit = mfile.get("dz_cs_turn_conduit", scan=scan)
    dr_cs_turn_conduit = mfile.get("dr_cs_turn_conduit", scan=scan)
    radius_cs_turn_corners = mfile.get("radius_cs_turn_corners", scan=scan)
    f_a_cs_turn_steel = mfile.get("f_a_cs_turn_steel", scan=scan)

    # Plot the CS turn as a rectangle representing the conductor cross-section
    # Assume dz_cs_turn is the diameter and dr_cs_turn is the length of the conductor
    # cross-section

    # Draw the conductor cross-section as a rectangle
    axis.add_patch(
        patches.FancyBboxPatch(
            (0, 0),
            dr_cs_turn,
            dz_cs_turn,
            boxstyle=patches.BoxStyle(
                "Round", pad=0, rounding_size=radius_cs_turn_corners
            ),
            edgecolor="black",
            facecolor="grey",
            lw=1.5,
            label="CS Turn Steel Conduit",
        )
    )

    # Draw the conductor cross-section as a rectangle
    axis.add_patch(
        patches.Rectangle(
            (
                dr_cs_turn_conduit + radius_cs_turn_cable_space,
                dz_cs_turn_conduit,
            ),
            dr_cs_turn - ((2 * dr_cs_turn_conduit) + (2 * radius_cs_turn_cable_space)),
            2 * radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            label="CS Turn Cable Space",
            zorder=2,
        )
    )
    # Plot the right hand circle for the CS turn cable space
    axis.add_patch(
        patches.Circle(
            (
                (dr_cs_turn - dr_cs_turn_conduit - radius_cs_turn_cable_space),
                dz_cs_turn / 2,
            ),
            radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            zorder=3,
        )
    )
    # Plot the left hand circle for the CS turn cable space
    axis.add_patch(
        patches.Circle(
            (
                (dr_cs_turn_conduit + radius_cs_turn_cable_space),
                dz_cs_turn / 2,
            ),
            radius_cs_turn_cable_space,
            facecolor="white",
            lw=1.5,
            zorder=3,
        )
    )

    # Add plasma volume, areas and shaping information
    textstr_turn = (
        f"$\\mathbf{{Turn \\ structure:}}$\n\n$A:$ {a_cs_turn:.4e}$ \\"
        f" \\text{{m}}^2$\nTurn width: {dr_cs_turn:.4e}$ \\ \\text{{m}}$\nTurn"
        f" height: {dz_cs_turn:.4e}$ \\ \\text{{m}}$\nTurn width to height"
        f" ratio: {f_dr_dz_cs_turn:.3f}\nSteel conduit width:"
        f" {dr_cs_turn_conduit:.4e}$ \\ \\text{{m}}$\nRadius of turn cable"
        f" space: {radius_cs_turn_cable_space:.4e}$ \\ \\text{{m}}$\nRadius of"
        f" turn corner: {radius_cs_turn_corners:.4e}$ \\"
        " \\text{m}$\nFraction of turn area that is steel:"
        f" {f_a_cs_turn_steel:.4f}\n"
    )

    axis.text(
        0.7,
        0.375,
        textstr_turn,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    axis.set_xlim(-dr_cs_turn * 0.2, dr_cs_turn * 1.2)
    axis.set_ylim(-dz_cs_turn * 0.3, dz_cs_turn * 1.3)
    axis.set_aspect("equal")
    axis.set_xlabel("Length [m]")
    axis.set_ylabel("Height [m]")
    axis.set_title("CS Turn Conductor Cross-Section")
    cs_legend = axis.legend(loc="upper right", bbox_to_anchor=(0.7, -0.25))
    cs_legend.get_frame().set_edgecolor("black")
    axis.grid(True, linestyle="--", alpha=0.3)


def plot_pf_cs_plasma_mutual_inductance(
    axis: plt.Axes, m_file: MFile, scan: int
) -> None:
    """Plot the mutual inductance between the plasma and PF/CS coils.

    Parameters
    ----------
    axis : plt.Axes
        Axis to plot on
    m_file : MFile
        MFILE data object
    scan : int
        Scan number to read from MFILE

    """
    n_pf_cs_plasma_circuits = int(m_file.get("n_pf_cs_plasma_circuits", scan=scan))
    mutual_inductance = np.zeros((n_pf_cs_plasma_circuits, n_pf_cs_plasma_circuits))
    iohcl = int(m_file.get("iohcl", scan=scan))

    for coil in range(n_pf_cs_plasma_circuits):
        for circuit in range(n_pf_cs_plasma_circuits):
            mutual_inductance[coil, circuit] = m_file.get(
                f"ind_pf_cs_plasma_mutual[{coil},_{circuit}]",
                scan=scan,
            )

    # Create lower triangular matrix
    mutual_inductance = np.tril(mutual_inductance)
    im = axis.imshow(mutual_inductance, cmap="RdBu_r", aspect="auto")
    axis.set_xlabel("Circuit")
    axis.set_ylabel("Circuit")
    axis.set_title("PF/CS Plasma Mutual Inductance")
    axis.set_xticks(range(n_pf_cs_plasma_circuits))
    axis.set_yticks(range(n_pf_cs_plasma_circuits))
    labels = list(range(1, n_pf_cs_plasma_circuits + 1))

    if iohcl == 1:
        labels[-2] = "CS"
    labels[-1] = "Plasma"
    axis.set_xticklabels(labels)
    axis.set_yticklabels(labels)

    # Add boxes around each cell
    for i in range(n_pf_cs_plasma_circuits):
        for j in range(n_pf_cs_plasma_circuits):
            if mutual_inductance[i, j] != 0:
                axis.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor="black",
                        linewidth=0.5,
                    )
                )
                # Add text annotation with values
                axis.text(
                    j,
                    i,
                    f"{mutual_inductance[i, j]:.3e}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8,
                )

    axis.get_figure().colorbar(im, ax=axis, label="Mutual Inductance (H)")


__all__ = [
    "plot_cs_coil_structure",
    "plot_cs_turn_structure",
    "plot_magnetics_info",
    "plot_pf_cs_plasma_mutual_inductance",
    "plot_physics_info",
    "secs_to_hms",
]
