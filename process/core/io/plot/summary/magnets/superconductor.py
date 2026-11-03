"""Superconductor critical surface state plots for PROCESS summary plots."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.common import box_style
from process.models.superconductors import SuperconductorModel

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile

__all__ = ["plot_tf_superconductor_state"]

N_GRID = 60
OPERATING_COLOUR = "tab:red"


def _evaluate(model: SuperconductorModel, b, temp, strain: float, **kwargs):
    """Evaluate the critical current density [A/m²] over broadcast `b` and `temp`.

    Points where a model is outside its range of validity are returned as NaN.
    """
    b, temp = np.broadcast_arrays(np.asarray(b, float), np.asarray(temp, float))
    out = np.full(b.shape, np.nan)
    # The models log out-of-range warnings, which are expected when scanning a grid
    previous = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    try:
        with np.errstate(all="ignore"):
            for idx in np.ndindex(b.shape):
                try:
                    out[idx] = model.critical_current_density(
                        b[idx], temp[idx], strain, **kwargs
                    )
                except Exception:  # noqa: BLE001
                    out[idx] = np.nan
    finally:
        logging.disable(previous)
    return out


def _crossing(x: np.ndarray, y: np.ndarray, level: float) -> float | None:
    """First `x` (linearly interpolated) where `y` falls through `level`."""
    valid = np.isfinite(y)
    x, y = x[valid], y[valid]
    below = np.nonzero((y[:-1] >= level) & (y[1:] < level))[0]
    if below.size == 0:
        return None
    i = below[0]
    return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))


def plot_tf_superconductor_state(fig: plt.Figure, mfile: MFile, scan: int):
    """Plot the TF superconductor operating state on 2D slices of its critical surface.

    Three panels are drawn: the B-T plane with the J_c contours and the critical
    boundary at the operating current density, and the J_c-B and J_c-T load lines.
    The operating point is marked on each.

    Parameters
    ----------
    fig:
        Figure to draw on.
    mfile:
        MFile containing the PROCESS run.
    scan:
        Scan number to plot.
    """
    model = SuperconductorModel(mfile.get("i_tf_sc_mat", scan=scan))
    b_op = mfile.get("b_tf_inboard_peak_with_ripple", scan=scan)
    t_op = mfile.get("tftmp", scan=scan)
    j_crit_op = mfile.get("j_tf_superconductor_critical", scan=scan)
    j_op = j_crit_op * mfile.get("f_c_tf_turn_operating_critical", scan=scan)

    # str_wp is only populated when the self-consistent strain is in use
    strain = mfile.get("str_wp", scan=scan) or mfile.get("str_tf_con_res", scan=scan)

    kwargs = {}
    b_c20 = mfile.get("b_tf_superconductor_critical_zero_temp_strain", scan=scan)
    t_c0 = mfile.get("temp_tf_superconductor_critical_zero_field_strain", scan=scan)
    if model.b_crit_zero_field_strain is not None and b_c20 and t_c0:
        kwargs.update(b_c20max=b_c20, temp_c0max=t_c0)
    if model == SuperconductorModel.HAZELTON_ZHAI_REBCO:
        kwargs.update(
            dr_hts_tape=mfile.get("dr_tf_hts_tape", scan=scan),
            dx_hts_tape_rebco=mfile.get("dx_tf_hts_tape_rebco", scan=scan),
            dx_hts_tape_total=mfile.get("dx_tf_hts_tape_total", scan=scan),
        )

    # Models fitted to data (Bi-2212) have no defined critical field or temperature
    b_max = 1.6 * b_op
    t_max = max(20.0 if not t_c0 else min(t_c0, 0.4 * t_c0 + 12.0), 2.0 * t_op)

    bs = np.linspace(0.3 * b_op, b_max, N_GRID)
    ts = np.linspace(2.0, t_max, N_GRID)
    b_grid, t_grid = np.meshgrid(bs, ts)
    j_grid = _evaluate(model, b_grid, t_grid, strain, **kwargs)

    # Critical-current margin: J_c along B at fixed T and along T at fixed B
    j_vs_b = _evaluate(model, bs, t_op, strain, **kwargs)
    j_vs_t = _evaluate(model, b_op, ts, strain, **kwargs)
    b_cs = _crossing(bs, j_vs_b, j_op)
    t_cs = _crossing(ts, j_vs_t, j_op)

    ax_bt = fig.add_subplot(221)
    ax_jb = fig.add_subplot(223)
    ax_jt = fig.add_subplot(224)
    ax_bt.set_position([0.08, 0.5, 0.55, 0.43])
    ax_jb.set_position([0.08, 0.08, 0.38, 0.32])
    ax_jt.set_position([0.57, 0.08, 0.38, 0.32])

    # B-T plane
    j_plot_max = 3.0 * j_op  # Keeps the colour scale useful near the operating point
    levels = np.linspace(0.0, j_plot_max / 1.0e6, 16)
    contour = ax_bt.contourf(
        b_grid, t_grid, j_grid / 1.0e6, levels=levels, cmap="viridis", extend="max"
    )
    fig.colorbar(contour, ax=ax_bt, label="$J_c$ [MA/m$^2$]")
    if np.nanmin(j_grid) < j_op < np.nanmax(j_grid):
        ax_bt.contour(
            b_grid,
            t_grid,
            j_grid,
            levels=[j_op],
            colors="white",
            linestyles="--",
            linewidths=2,
        )
    ax_bt.plot(
        [], [], "--", c="white", label=f"$J_c = J_{{op}}$ ({j_op / 1e6:.0f} MA/m²)"
    )
    ax_bt.plot(
        b_op,
        t_op,
        "*",
        c=OPERATING_COLOUR,
        ms=16,
        mec="white",
        label="Operating point",
    )
    ax_bt.set_xlabel("Magnetic field $B$ [T]")
    ax_bt.set_ylabel("Temperature $T$ [K]")
    ax_bt.set_title(
        f"{model.full_name}: critical surface at $\\varepsilon$ = {strain:.2e}"
    )
    ax_bt.legend(loc="upper right", facecolor="lightgrey", framealpha=0.85)

    # J_c(B) load line
    ax_jb.plot(bs, j_vs_b / 1.0e6, c="tab:blue", label=f"$J_c(B)$, T = {t_op:.2f} K")
    ax_jb.plot(b_op, j_op / 1.0e6, "*", c=OPERATING_COLOUR, ms=14, label="Operating")
    if b_cs is not None:
        ax_jb.axvline(b_cs, c="grey", ls=":")
        ax_jb.annotate(
            f"$B_{{margin}}$ = {b_cs - b_op:.2f} T",
            (b_cs, j_op / 1.0e6),
            xytext=(0.55, 0.8),
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "->"},
            bbox=box_style("white"),
        )
    ax_jb.set_xlabel("Magnetic field $B$ [T]")
    ax_jb.set_ylabel("$J_c$ [MA/m$^2$]")
    ax_jb.set_ylim(0.0, j_plot_max / 1.0e6)
    ax_jb.legend()
    ax_jb.grid(alpha=0.3)

    # J_c(T) load line
    ax_jt.plot(ts, j_vs_t / 1.0e6, c="tab:blue", label=f"$J_c(T)$, B = {b_op:.2f} T")
    ax_jt.plot(t_op, j_op / 1.0e6, "*", c=OPERATING_COLOUR, ms=14, label="Operating")
    if t_cs is not None:
        ax_jt.axvline(t_cs, c="grey", ls=":")
        ax_jt.annotate(
            f"$T_{{cs}}$ = {t_cs:.2f} K\n$\\Delta T$ = {t_cs - t_op:.2f} K",
            (t_cs, j_op / 1.0e6),
            xytext=(0.55, 0.65),
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "->"},
            bbox=box_style("white"),
        )
    ax_jt.set_xlabel("Temperature $T$ [K]")
    ax_jt.set_ylabel("$J_c$ [MA/m$^2$]")
    ax_jt.set_ylim(0.0, j_plot_max / 1.0e6)
    ax_jt.legend()
    ax_jt.grid(alpha=0.3)
