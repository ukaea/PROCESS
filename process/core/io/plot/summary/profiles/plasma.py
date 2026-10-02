"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from process.core import constants
from process.core.io.plot.summary.common import box_style, text_layout
from process.core.io.plot.summary.plasma.physics import reaction_plot_grid
from process.core.io.plot.summary.profiles.misc import interp1d_profile
from process.core.io.plot.summary.rendering import draw_text
from process.data_structure.impurity_radiation_variables import ImpurityRadiationData
from process.models.physics.profiles import PlasmaProfileShapeType

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def plot_n_profiles(prof, demo_ranges: bool, mfile: MFile, scan: int):
    """Function to plot density profile

    Parameters
    ----------
    prof :
        axis object to add plot to
    demo_ranges: bool :

    mfile: MFile :

    scan: int :

    """
    nd_alphas = mfile.get("nd_plasma_alphas_thermal_vol_avg", scan=scan)
    nd_protons = mfile.get("nd_plasma_protons_vol_avg", scan=scan)
    nd_impurities = mfile.get("nd_plasma_impurities_vol_avg", scan=scan)
    nd_ions_total = mfile.get("nd_plasma_ions_total_vol_avg", scan=scan)
    nd_fuel_ions = mfile.get("nd_plasma_fuel_ions_vol_avg", scan=scan)
    alphan = mfile.get("alphan", scan=scan)
    f_nd_plasma_pedestal_greenwald = mfile.get(
        "f_nd_plasma_pedestal_greenwald", scan=scan
    )
    f_nd_plasma_separatrix_greenwald = mfile.get(
        "f_nd_plasma_separatrix_greenwald", scan=scan
    )
    nd_plasma_electrons_vol_avg = mfile.get("nd_plasma_electrons_vol_avg", scan=scan)
    # find impurity densities
    imp_frac = np.array([
        mfile.get(f"f_nd_impurity_electrons({i:02d})", scan=scan) for i in range(1, 15)
    ])

    nd_plasma_separatrix_electron = mfile.get("nd_plasma_separatrix_electron", scan=scan)

    ax_main = prof.add_subplot(631)
    ax_main.set_position([0.075, 0.625, 0.25, 0.325])
    ax_impurity = prof.add_subplot(634, sharex=ax_main)
    ax_impurity.set_position([0.075, 0.275, 0.25, 0.325])
    ax_main.tick_params(labelbottom=False)

    ax_impurity.set_xlabel(r"$\rho \quad [r/a]$")
    ax_main.set_ylabel(r"$n \ [10^{19}\ \mathrm{m}^{-3}]$")
    ax_impurity.set_ylabel(r"$n \ [10^{16}\ \mathrm{m}^{-3}]$")
    ax_main.set_title("Density profile")

    i_plasma_pedestal = mfile.get("i_plasma_pedestal", scan=scan)
    nd_plasma_pedestal_electron = mfile.get("nd_plasma_pedestal_electron", scan=scan)
    ne0 = mfile.get("nd_plasma_electron_on_axis", scan=scan)
    nd_plasma_electrons_vol_avg = mfile.get("nd_plasma_electrons_vol_avg", scan=scan)
    radius_plasma_pedestal_density_norm = mfile.get(
        "radius_plasma_pedestal_density_norm", scan=scan
    )
    ne0 = mfile.get("nd_plasma_electron_on_axis", scan=scan)
    n_plasma_profile_elements = mfile.get("n_plasma_profile_elements", scan=scan)

    # build electron profile and species profiles (scale with electron profile shape)
    if i_plasma_pedestal == 1:
        rho = np.linspace(0, 1.0, int(n_plasma_profile_elements))
        ne = np.zeros_like(rho)

        for i in range(len(rho)):
            if rho[i] <= radius_plasma_pedestal_density_norm:
                ne[i] = (
                    nd_plasma_pedestal_electron
                    + (ne0 - nd_plasma_pedestal_electron)
                    * (1 - rho[i] ** 2 / radius_plasma_pedestal_density_norm**2)
                    ** alphan
                )
            else:
                ne[i] = nd_plasma_separatrix_electron + (
                    nd_plasma_pedestal_electron - nd_plasma_separatrix_electron
                ) * (1 - rho[i]) / (1 - min(0.9999, radius_plasma_pedestal_density_norm))
    else:
        rho = np.linspace(0, 1.0, n_plasma_profile_elements)
        ne = ne0 * (1 - rho**2) ** alphan

    # species profiles scaled by their average fraction relative to electrons

    if nd_plasma_electrons_vol_avg != 0:
        fracs = (
            np.array([
                nd_fuel_ions,
                nd_alphas,
                nd_protons,
                nd_impurities,
                nd_ions_total,
                nd_plasma_electrons_vol_avg,
            ])
            / nd_plasma_electrons_vol_avg
        )
    else:
        fracs = np.zeros(5)

    # build species density profiles from electron profile and fractions
    # fracs = [fuel, alpha, protons, impurities, ions_total]
    # Create a density profile for each species by multiplying ne by each fraction in
    # fracs
    density_profiles = np.array([ne * frac for frac in fracs])

    # convert to 1e19 m^-3 units for plotting (vectorised)
    density_profiles_plotting = density_profiles / 1e19

    ax_main.plot(
        rho,
        density_profiles_plotting[0],
        label=r"$n_{\text{fuel}}$",
        color="#2ca02c",
        linewidth=1.5,
    )
    ax_main.plot(
        rho,
        density_profiles_plotting[1],
        label=r"$n_{\alpha,\text{thermal}}$",
        color="#d62728",
        linewidth=1.5,
    )
    ax_impurity.plot(
        rho,
        density_profiles_plotting[2] * 1e3,
        label=r"$n_{p}$",
        color="#17becf",
        linewidth=1.5,
    )
    ax_impurity.plot(
        rho,
        density_profiles_plotting[3] * 1e3,
        label=r"$n_{imp,total}$",
        color="#9467bd",
        linewidth=2.5,
        linestyle="dotted",
    )
    ax_main.plot(
        rho,
        density_profiles_plotting[4],
        label=r"$n_{i,total}$",
        color="#ff7f0e",
        linewidth=1.5,
    )
    ax_main.plot(
        rho,
        density_profiles_plotting[5],
        label=r"$n_{e}$",
        color="blue",
        linewidth=1.5,
    )

    imp_labels = ImpurityRadiationData().imp_label
    for ind in range(2, imp_frac.shape[0]):
        lbl = imp_labels[ind].replace("_", "")
        if imp_frac[ind] > 1.0e-30:
            ax_impurity.plot(
                rho, imp_frac[ind] * ne / 1e16, label=rf"$n_{{\text{{{lbl}}}}}$"
            )

    ax_main.legend(loc="best")
    ax_impurity.legend(loc="best")

    # Ranges
    # ---
    # DEMO : Fixed ranges for comparison
    ax_main.set_xlim(0, 1)
    ax_impurity.set_xlim(0, 1)
    if demo_ranges:
        ax_main.set_ylim(0, 20)

    # Adaptive ranges
    else:
        ax_main.set_ylim(0, ax_main.get_ylim()[1])
        # Use logarithmic scale for impurity axis if any impurity values are very small
        impurity_data = [
            imp_frac[i] * ne / 1e16
            for i in range(len(imp_frac))
            if imp_frac[i] > 1.0e-30
        ]
        if impurity_data and np.min(impurity_data) / np.max(impurity_data) < 0.01:
            # If range spans more than 100x, use log scale
            ax_impurity.set_yscale("log")
        ax_impurity.set_ylim(1e-3, ax_impurity.get_ylim()[1])

    if i_plasma_pedestal != 0:
        # Print pedestal lines
        ax_main.axhline(
            y=nd_plasma_pedestal_electron / 1e19,
            xmax=radius_plasma_pedestal_density_norm,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        ax_main.vlines(
            x=radius_plasma_pedestal_density_norm,
            ymin=0.0,
            ymax=nd_plasma_pedestal_electron / 1e19,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
    ax_main.minorticks_on()
    ax_impurity.minorticks_on()

    # Add text box with density profile parameters
    textstr_density = "\n".join((
        (
            r"$\langle n_{\text{e}} \rangle$:"
            rf" {nd_plasma_electrons_vol_avg:.3e}"
            r" m$^{-3}$"
            r"$\hspace{4} \overline{n_{e}}$:"
            rf" {mfile.get('nd_plasma_electron_line', scan=scan):.3e}"
            r" m$^{-3}$"
        ),
        (
            rf"$n_{{\text{{e,0}}}}$: {ne0:.3e} m$^{{-3}}$"
            rf"$\hspace{{4}} \alpha_{{\text{{n}}}}$: {alphan:.3f}"
        ),
        (
            rf"$n_{{\text{{e,ped}}}}$: {nd_plasma_pedestal_electron:.3e}"
            r" m$^{-3}$"
            r"$ \hspace{3} \frac{\langle n_i \rangle}{\langle n_e"
            r" \rangle}$: "
            f"{nd_fuel_ions / nd_plasma_electrons_vol_avg:.3f}"
        ),
        (
            r"$f_{\text{GW e,ped}}$:"
            rf" {f_nd_plasma_pedestal_greenwald:.3f}"
            r"$ \hspace{7} \frac{n_{e,0}}{\langle n_e \rangle}$: "
            f"{ne0 / nd_plasma_electrons_vol_avg:.3f}"
        ),
        (
            r"$\rho_{\text{ped,n}}$:"
            rf" {radius_plasma_pedestal_density_norm:.3f}"
            r"$ \hspace{8} \frac{\overline{n_{e}}}{n_{\text{GW}}}$: "
            f"{mfile.get('nd_plasma_electron_line', scan=scan) / mfile.get('nd_plasma_electron_max_array(7)', scan=scan):.3f}"  # noqa: E501
        ),
        (
            rf"$n_{{\text{{e,sep}}}}$: {nd_plasma_separatrix_electron:.3e}"
            r" m$^{-3}$"
        ),
        (
            r"$f_{\text{GW e,sep}}$:"
            rf" {f_nd_plasma_separatrix_greenwald:.3f}"
        ),
    ))

    props_density = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    ax_main.text(
        -0.05,
        -0.175,
        textstr_density,
        transform=ax_impurity.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_density,
    )

    textstr_ions = "\n".join((
        (
            r"$\langle n_{\text{ions-total}} \rangle $: "
            f"{mfile.get('nd_plasma_ions_total_vol_avg', scan=scan):.3e}"
            " m$^{-3}$"
        ),
        (
            r"$\langle n_{\text{fuel}} \rangle $: "
            f"{mfile.get('nd_plasma_fuel_ions_vol_avg', scan=scan):.3e}"
            " m$^{-3}$"
        ),
        (
            r"$\langle n_{\alpha,\text{thermal}} \rangle $: "
            f"{mfile.get('nd_plasma_alphas_thermal_vol_avg', scan=scan):.3e}"
            " m$^{-3}$"
        ),
        (
            r"$\langle n_{\text{impurities}} \rangle $: "
            f"{mfile.get('nd_plasma_impurities_vol_avg', scan=scan):.3e}"
            " m$^{-3}$"
        ),
        (
            r"$\langle n_{\text{protons}} \rangle $:"
            f"{mfile.get('nd_plasma_protons_vol_avg', scan=scan):.3e}"
            " m$^{-3}$"
        ),
    ))

    ax_impurity.text(
        1.2,
        0.05,
        textstr_ions,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=ax_impurity.transAxes,
        bbox={
            "boxstyle": "round",
            "facecolor": "wheat",
            "alpha": 0.5,
        },
    )

    ax_main.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    ax_impurity.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)


def plot_jprofile(prof, mfile: MFile, scan: int):
    """Function to plot density profile

    Parameters
    ----------
    prof :
        axis object to add plot to
    mfile: MFile :

    scan: int :

    """
    alphaj = mfile.get("alphaj", scan=scan)
    j_plasma_0 = mfile.get("j_plasma_on_axis", scan=scan)
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    j_plasma_bootstrap_sauter_profile = [
        mfile.get(f"j_plasma_bootstrap_sauter_profile{i}", scan=scan) / 1000.0
        for i in range(n_plasma_profile_elements - 3)
    ]

    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel(r"Current density $[kA/m^2]$")
    prof.set_title("$J$ profile")
    prof.minorticks_on()
    prof.set_xlim(0, 1.0)

    rho = np.linspace(0, 1)
    y2 = (j_plasma_0 * (1 - rho**2) ** alphaj) / 1e3

    prof.plot(rho, y2, color="red")

    prof.plot(
        np.linspace(0, 1, n_plasma_profile_elements - 3),
        j_plasma_bootstrap_sauter_profile,
        label="Sauter Bootstrap",
        color="green",
        linestyle="--",
    )
    prof.legend()

    textstr_j = "\n".join((
        r"$j_0$: " + f"{y2[0]:.3f} kA m$^{{-2}}$\n",
        r"$\alpha_J$: " + f"{alphaj:.3f}",
    ))

    props_j = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    prof.text(
        0.65,
        1.6,
        textstr_j,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_j,
    )

    prof.text(
        0.35,
        0.04,
        "*Current profile is assumed to be parabolic",
        fontsize=10,
        ha="left",
        transform=plt.gcf().transFigure,
    )
    prof.text(
        0.35,
        0.02,
        "*Bootstrap profile is for representation only",
        fontsize=10,
        ha="left",
        transform=plt.gcf().transFigure,
    )
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)


def plot_t_profiles(prof, demo_ranges: bool, mfile: MFile, scan: int):
    """Function to plot temperature profile

    Parameters
    ----------
    prof :
        axis object to add plot to
    demo_ranges: bool :

    mfile: MFile :

    scan: int :

    """
    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel("$T$ [keV]")
    prof.set_title("Temperature profile")

    alphat = mfile.get("alphat", scan=scan)
    radius_plasma_pedestal_temp_norm = mfile.get(
        "radius_plasma_pedestal_temp_norm", scan=scan
    )

    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))
    i_plasma_pedestal = mfile.get("i_plasma_pedestal", scan=scan)
    rho = np.linspace(0, 1.0, n_plasma_profile_elements)
    temp_plasma_pedestal_electron_kev = mfile.get(
        "temp_plasma_pedestal_electron_kev", scan=scan
    )
    temp_plasma_separatrix_electron_kev = mfile.get(
        "temp_plasma_separatrix_electron_kev", scan=scan
    )
    f_temp_plasma_ion_electron = mfile.get("f_temp_plasma_ion_electron", scan=scan)
    tbeta = mfile.get("tbeta", scan=scan)
    te0 = mfile.get("temp_plasma_electron_on_axis_kev", scan=scan)

    if i_plasma_pedestal == 1:
        rhocore = np.linspace(0.0, radius_plasma_pedestal_temp_norm)
        tcore = (
            temp_plasma_pedestal_electron_kev
            + (te0 - temp_plasma_pedestal_electron_kev)
            * (1 - (rhocore / radius_plasma_pedestal_temp_norm) ** tbeta) ** alphat
        )

        rhosep = np.linspace(radius_plasma_pedestal_temp_norm, 1)
        tsep = temp_plasma_separatrix_electron_kev + (
            temp_plasma_pedestal_electron_kev - temp_plasma_separatrix_electron_kev
        ) * (1 - rhosep) / (1 - min(0.9999, radius_plasma_pedestal_temp_norm))

        rho = np.append(rhocore, rhosep)
        te = np.append(tcore, tsep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        te = te0 * (1 - rho**2) ** alphat
    prof.plot(rho, te, color="blue", label="$T_{e}$")
    prof.plot(rho, te[:] * f_temp_plasma_ion_electron, color="red", label="$T_{i}$")
    prof.legend()

    # Ranges
    # ---
    prof.set_xlim(0, 1)
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim(0, 50)

    # Adaptive ranges
    else:
        prof.set_ylim(0, prof.get_ylim()[1])

    if i_plasma_pedestal != 0:
        # Plot pedestal lines
        prof.axhline(
            y=temp_plasma_pedestal_electron_kev,
            xmax=radius_plasma_pedestal_temp_norm,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        prof.vlines(
            x=radius_plasma_pedestal_temp_norm,
            ymin=0.0,
            ymax=temp_plasma_pedestal_electron_kev,
            color="r",
            linestyle="-",
            linewidth=0.4,
            alpha=0.4,
        )
        prof.minorticks_on()

    te = mfile.get("temp_plasma_electron_vol_avg_kev", scan=scan)
    # Add text box with temperature profile parameters
    textstr_temperature = "\n".join((
        (
            r"$\langle T_{\text{e}} \rangle_\text{V}$: "
            rf" {mfile.get('temp_plasma_electron_vol_avg_kev', scan=scan):.3f} keV"
            r"$\hspace{2} \langle T_{\text{e}} \rangle_\text{n}$:"
            rf" {mfile.get('temp_plasma_electron_density_weighted_kev', scan=scan):.3f} keV"  # noqa: E501
            r"$\hspace{2} \overline{T_{e}}$:"
            rf" {mfile.get('temp_plasma_electron_line_avg_kev', scan=scan):.3f} keV"
        ),
        (
            rf"$T_{{\text{{e,0}}}}$:    {te0:.3f} keV"
            rf"$\hspace{{3}} \alpha_{{\text{{T}}}}$: {alphat:.3f}    "
            r"$\hspace{3} \langle T_{\text{i}} \rangle_\text{V}$:"
            rf" {mfile.get('temp_plasma_ion_vol_avg_kev', scan=scan):.3f} keV"
        ),
        (
            r"$T_{\text{e,ped}}$:"
            rf" {temp_plasma_pedestal_electron_kev:.3f} keV"
            r"$  \hspace{3} \frac{\langle T_i \rangle}{\langle T_e"
            r" \rangle}$: "
            f"{f_temp_plasma_ion_electron:.3f}   "
            "$\\hspace{4} T_{\\text{i,0}}$:"
            f" {mfile.get('temp_plasma_ion_on_axis_kev', scan=scan):.3f} keV"
        ),
        (
            r"$\rho_{\text{ped,T}}$:"
            rf" {radius_plasma_pedestal_temp_norm:.3f}"
            r"$ \hspace{5} \frac{T_{e,0}}{\langle T_e \rangle}$: "
            f"{mfile.get('f_temp_plasma_electron_on_axis_vol_avg', scan=scan):.3f}  "
            "$\\hspace{4} T_{\\text{i,ped}}$:"
            f" {mfile.get('temp_plasma_pedestal_ion_kev', scan=scan):.3f} keV"
        ),
        (
            r"$T_{\text{e,sep}}$:"
            rf" {temp_plasma_separatrix_electron_kev:.3f} keV"
            r"$\hspace{3} \frac{{{\langle T_e \rangle_n}}}{{{\langle T_e"
            r" \rangle_V}}}$: "
            f"{mfile.get('f_temp_plasma_electron_density_vol_avg', scan=scan):.3f}"
            "$\\hspace{4} T_{\\text{i,sep}}$:"
            f" {mfile.get('temp_plasma_separatrix_ion_kev', scan=scan):.3f} keV"
        ),
    ))

    props_temperature = {
        "boxstyle": "round",
        "facecolor": "wheat",
        "alpha": 0.5,
    }
    prof.text(
        -0.1,
        -0.125,
        textstr_temperature,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_temperature,
    )
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)


def plot_qprofile(prof, demo_ranges: bool, mfile: MFile, scan: int):
    """Function to plot q profile, formula taken from Nevins bootstrap model.

    Parameters
    ----------
    prof :
        axis object to add plot to
    demo_ranges: bool :

    mfile: MFile :

    scan: int :

    """
    prof.set_xlabel(r"$\rho \quad [r/a]$")
    prof.set_ylabel("$q$")
    prof.set_title("$q$ profile")
    prof.minorticks_on()

    rho = np.linspace(0, 1)
    q0 = mfile.get("q0", scan=scan)
    q95 = mfile.get("q95", scan=scan)

    q_r_nevin = q0 + (q95 - q0) * (rho + rho * rho + rho**3) / (3.0)
    q_r_sauter = q0 + (q95 - q0) * (rho * rho)

    prof.plot(rho, q_r_nevin, label="Nevins")
    prof.plot(rho, q_r_sauter, label="Sauter")
    prof.legend()

    # Ranges
    # ---
    prof.set_xlim(0, 1)
    # DEMO : Fixed ranges for comparison
    if demo_ranges:
        prof.set_ylim(0, 10)

    # Adaptive ranges
    else:
        prof.set_ylim(0, q95 * 1.2)

    prof.text(
        0.6,
        0.04,
        "*Profile is not calculated, only $q_0$ and $q_{95}$ are known.",
        fontsize=10,
        ha="left",
        transform=plt.gcf().transFigure,
    )
    prof.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.2)
    # ---

    textstr_q = " | ".join((
        r"$q_0$: " + f"{q0:.3f}",
        r"$q_{95}$: " + f"{q95:.3f}",
        r"$q_{\text{cyl}}$: " + f"{mfile.get('qstar', scan=scan):.3f}",
    ))

    props_q = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}
    prof.text(
        0.0,
        1.4,
        textstr_q,
        transform=prof.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=props_q,
    )


def plot_fusion_rate_profiles(axis: plt.Axes, fig, mfile: MFile, scan: int):
    """Plot the fusion rate density profiles on the given axis"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    fusden_plasma_dt_profile = [
        mfile.get(f"fusden_plasma_dt_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]

    fusden_plasma_dd_triton_profile = [
        mfile.get(f"fusden_plasma_dd_triton_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]

    fusden_plasma_dd_helion_profile = [
        mfile.get(f"fusden_plasma_dd_helion_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    fusden_plasma_dhe3_profile = [
        mfile.get(f"fusden_plasma_dhe3_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]

    fusrat_plasma_total_profile = [
        fusden_plasma_dt_profile[i]
        + fusden_plasma_dd_triton_profile[i]
        + fusden_plasma_dd_helion_profile[i]
        + fusden_plasma_dhe3_profile[i]
        for i in range(len(fusden_plasma_dt_profile))
    ]

    axis.spines["left"].set_color("red")
    axis.yaxis.label.set_color("black")
    axis.tick_params(axis="y", colors="red")

    # Plot fusion rates (dashed lines, left axis) with axis color and different
    # linestyles
    axis.plot(
        np.linspace(0, 1, len(fusden_plasma_dt_profile)),
        fusden_plasma_dt_profile,
        color=axis.spines["left"].get_edgecolor(),
        linestyle="-",
        label=r"$\mathrm{D-T}$",
    )
    axis.plot(
        np.linspace(0, 1, len(fusden_plasma_dd_triton_profile)),
        fusden_plasma_dd_triton_profile,
        color=axis.spines["left"].get_edgecolor(),
        linestyle=":",
        label=r"$\mathrm{D-D \ Triton}$",
    )
    axis.plot(
        np.linspace(0, 1, len(fusden_plasma_dd_helion_profile)),
        fusden_plasma_dd_helion_profile,
        color=axis.spines["left"].get_edgecolor(),
        linestyle="-.",
        label=r"$\mathrm{D-D \ Helion}$",
    )
    axis.plot(
        np.linspace(0, 1, len(fusden_plasma_dhe3_profile)),
        fusden_plasma_dhe3_profile,
        color=axis.spines["left"].get_edgecolor(),
        linestyle="--",
        label=r"$\mathrm{D-3He}$",
    )
    axis.plot(
        np.linspace(0, 1, len(fusrat_plasma_total_profile)),
        fusrat_plasma_total_profile,
        color=axis.spines["left"].get_edgecolor(),
        linestyle="None",
        marker="d",
        markersize=1,
        label=r"Total",
    )

    # Show the plasma volume-averaged rate density and its position on the
    # profile.
    profile_positions = np.linspace(0, 1, len(fusrat_plasma_total_profile))
    profile_rates = np.asarray(fusrat_plasma_total_profile)
    average_rate = mfile.get("fusden_plasma_vol_avg", scan=scan)
    axis.axhline(
        average_rate,
        color="black",
        linestyle="--",
        linewidth=0.9,
        label="Plasma volume average",
    )

    average_position = profile_positions[
        np.nanargmin(np.abs(profile_rates - average_rate))
    ]
    axis.axvline(
        average_position,
        color="black",
        linestyle="--",
        linewidth=0.9,
    )

    # Plot fusion power (solid lines, right axis) with axis color and different
    # linestyles
    ax2 = axis.twinx()
    ax2.spines["right"].set_color("blue")
    ax2.yaxis.label.set_color("black")
    ax2.tick_params(axis="y", colors="blue")
    ax2.plot(
        np.linspace(0, 1, len(fusden_plasma_dt_profile)),
        np.array(fusden_plasma_dt_profile) * constants.D_T_ENERGY,
        color=ax2.spines["right"].get_edgecolor(),
        linestyle="-",
    )

    ax2.plot(
        np.linspace(0, 1, len(fusden_plasma_dd_triton_profile)),
        np.array(fusden_plasma_dd_triton_profile) * constants.DD_TRITON_ENERGY,
        color=ax2.spines["right"].get_edgecolor(),
        linestyle=":",
    )
    ax2.plot(
        np.linspace(0, 1, len(fusden_plasma_dd_helion_profile)),
        np.array(fusden_plasma_dd_helion_profile) * constants.DD_HELIUM_ENERGY,
        color=ax2.spines["right"].get_edgecolor(),
        linestyle="-.",
    )
    ax2.plot(
        np.linspace(0, 1, len(fusden_plasma_dhe3_profile)),
        np.array(fusden_plasma_dhe3_profile) * constants.D_HELIUM_ENERGY,
        color=ax2.spines["right"].get_edgecolor(),
        linestyle="--",
    )
    ax2.plot(
        np.linspace(0, 1, len(fusrat_plasma_total_profile)),
        (
            np.array(fusden_plasma_dhe3_profile) * constants.D_HELIUM_ENERGY
            + np.array(fusden_plasma_dd_helion_profile) * constants.DD_HELIUM_ENERGY
            + np.array(fusden_plasma_dd_triton_profile) * constants.DD_TRITON_ENERGY
            + np.array(fusden_plasma_dt_profile) * constants.D_T_ENERGY
        ),
        color=ax2.spines["right"].get_edgecolor(),
        linestyle="None",
        marker="d",
        markersize=1,
        label=r"Total",
    )

    # =================================================

    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.set_ylabel("Fusion Rate Density [reactions/m³/sec]")
    axis.legend(
        loc="lower left",
        edgecolor="black",
        facecolor="white",
        labelcolor="black",
        framealpha=1.0,
        frameon=True,
    )
    axis.set_yscale("log")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.set_xlim(0, 1.025)
    axis.minorticks_on()
    axis.set_ylim(1e10, 1e23)
    axis.yaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=10))
    axis.yaxis.set_minor_locator(
        plt.LogLocator(base=10.0, subs=np.arange(1, 10) * 0.1, numticks=100)
    )
    axis.tick_params(axis="y", which="minor", colors="red")

    ax2.set_title("Fusion Rate and Fusion Power Density Profiles")
    ax2.set_ylabel("Fusion Power Density [W/m³]")
    ax2.set_yscale("log")
    ax2.minorticks_on()
    ax2.yaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=10))
    ax2.yaxis.set_minor_locator(
        plt.LogLocator(base=10.0, subs=np.arange(1, 10) * 0.1, numticks=100)
    )
    ax2.tick_params(axis="y", which="minor", colors="blue")

    # =================================================

    # Add plasma volume, areas and shaping information
    textstr_general = (
        f"Total fusion rate: {mfile.get('fusrat_total', scan=scan):.4e}"
        " reactions/s\nTotal volume averaged fusion rate density:"
        f" {mfile.get('fusden_total_vol_avg', scan=scan):.4e}"
        " reactions/m3/s\nPlasma volume averaged fusion rate density:"
        f" {mfile.get('fusden_plasma_vol_avg', scan=scan):.4e}"
        " reactions/m3/s\n"
    )

    draw_text(
        axis,
        0.05,
        0.85,
        textstr_general,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    # ============================================================================

    textstr_dt = (
        f"Total fusion power: {mfile.get('p_dt_total_mw', scan=scan):,.2f}"
        " MW\nPlasma fusion power:"
        f" {mfile.get('p_plasma_dt_mw', scan=scan):,.2f} MW\nVolume-averaged"
        " fusion power density: plasma:"
        f" {mfile.get('pden_plasma_dt_vol_avg_mw', scan=scan):,.3f}"
        " MW/m³\nBeam fusion power:"
        f" {mfile.get('p_beam_dt_mw', scan=scan):,.2f} MW\n"
    )

    draw_text(
        axis,
        0.05,
        0.75,
        textstr_dt,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    draw_text(
        axis,
        0.24,
        0.8,
        "$\\text{D - T}$",
        fontsize=20,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =================================================

    textstr_dd = (
        f"Total fusion power: {mfile.get('p_dd_total_mw', scan=scan):,.2f}"
        " MW\nVolume-averaged total power density:"
        f" {mfile.get('pden_dd_total_vol_avg_mw', scan=scan):,.3e}"
        " MW/m³\nTritium branching ratio:"
        f" {mfile.get('f_dd_branching_trit', scan=scan):.4f}\n"
    )

    draw_text(
        axis,
        0.05,
        0.65,
        textstr_dd,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    draw_text(
        axis,
        0.22,
        0.685,
        "$\\text{D - D}$",
        fontsize=20,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =================================================

    textstr_dhe3 = (
        f"Total fusion power: {mfile.get('p_dhe3_total_mw', scan=scan):,.2f}"
        " MW\n\nVolume-averaged total power density:"
        f" {mfile.get('pden_dhe3_total_vol_avg_mw', scan=scan):,.3e} MW/m³\n\n"
    )

    draw_text(
        axis,
        0.05,
        0.55,
        textstr_dhe3,
        **text_layout(fig),
        bbox=box_style("lightyellow"),
    )

    draw_text(
        axis,
        0.21,
        0.59,
        "$\\text{D - 3He}$",
        fontsize=20,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =================================================

    textstr_alpha = (
        f"Total power: {mfile.get('p_alpha_total_mw', scan=scan):.2f}"
        f" MW\nPlasma power: {mfile.get('p_plasma_alpha_mw', scan=scan):.2f}"
        f" MW\nBeam power: {mfile.get('p_beam_alpha_mw', scan=scan):.2f}"
        " MW\n\nVolume-averaged rate density total:"
        f" {mfile.get('fusden_alpha_total_vol_avg', scan=scan):.4e}"
        " particles/m3/sec\nVolume-averaged rate density, plasma:"
        f" {mfile.get('fusden_plasma_alpha_vol_avg', scan=scan):.4e}"
        " particles/m3/sec\n\nVolume-averaged total power density:"
        f" {mfile.get('pden_alpha_total_vol_avg_mw', scan=scan):.4e}"
        " MW/m3\nVolume-averaged plasma power density:"
        f" {mfile.get('pden_plasma_alpha_vol_avg_mw', scan=scan):.4e}"
        " MW/m3\n\nPower per unit volume transferred to electrons:"
        f" {mfile.get('f_pden_alpha_electron_mw', scan=scan):.4e} MW/m3\nPower"
        " per unit volume transferred to ions:"
        f" {mfile.get('f_pden_alpha_ions_mw', scan=scan):.4e} MW/m3\n\n"
    )

    draw_text(
        axis,
        0.05,
        0.25,
        textstr_alpha,
        **text_layout(fig),
        bbox=box_style("red"),
    )

    draw_text(
        axis,
        0.35,
        0.45,
        "$\\alpha$",
        fontsize=22,
        verticalalignment="top",
        transform=fig.transFigure,
    )

    # =================================================

    textstr_neutron = (
        f"Total power: {mfile.get('p_neutron_total_mw', scan=scan):,.2f}"
        " MW\nPlasma power:"
        f" {mfile.get('p_plasma_neutron_mw', scan=scan):,.2f} MW\nBeam power:"
        f" {mfile.get('p_beam_neutron_mw', scan=scan):,.2f}"
        " MW\n\nVolume-averaged total power density:"
        f" {mfile.get('pden_neutron_total_vol_avg_mw', scan=scan):,.4e}"
        " MW/m3\nVolume-averaged plasma power density:"
        f" {mfile.get('pden_plasma_neutron_vol_avg_mw', scan=scan):,.4e}"
        " MW/m3\n"
    )

    draw_text(
        axis,
        0.05,
        0.1,
        textstr_neutron,
        **text_layout(fig),
        bbox=box_style("grey"),
    )

    draw_text(
        axis,
        0.25,
        0.2,
        "$n$",
        fontsize=20,
        verticalalignment="top",
        transform=fig.transFigure,
    )


def plot_plasma_pressure_profiles(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot the plasma pressure profiles on the given axis"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    pres_plasma_profile = [
        mfile.get(f"pres_plasma_electron_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_ion = [
        mfile.get(f"pres_plasma_ion_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_thermal_total_profile = [
        mfile.get(f"pres_plasma_thermal_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_fuel = [
        mfile.get(f"pres_plasma_fuel_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_kpa = [p / 1000.0 for p in pres_plasma_profile]
    pres_plasma_profile_ion_kpa = [p / 1000.0 for p in pres_plasma_profile_ion]
    pres_plasma_profile_fuel_kpa = [p / 1000.0 for p in pres_plasma_profile_fuel]
    pres_plasma_profile_total_kpa = [
        p / 1000.0 for p in pres_plasma_thermal_total_profile
    ]

    axis.plot(
        np.linspace(0, 1, len(pres_plasma_profile_kpa)),
        pres_plasma_profile_kpa,
        color="blue",
        label="Electron",
    )
    axis.plot(
        np.linspace(0, 1, len(pres_plasma_profile_ion_kpa)),
        pres_plasma_profile_ion_kpa,
        color="Red",
        label="Ion-total",
    )
    axis.plot(
        np.linspace(0, 1, len(pres_plasma_profile_fuel_kpa)),
        pres_plasma_profile_fuel_kpa,
        color="orange",
        label="Fuel",
    )
    axis.plot(
        np.linspace(0, 1, len(pres_plasma_profile_total_kpa)),
        pres_plasma_profile_total_kpa,
        color="green",
        label="Total",
    )

    # Plot horizontal line for volume-average thermal pressure (converted to kPa)
    p_vol_kpa = mfile.get("pres_plasma_thermal_vol_avg", scan=scan) / 1000.0
    axis.axhline(
        p_vol_kpa,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Volume avg",
        zorder=5,
    )

    axis.set_xlabel("$\\rho$ [r/a]")
    axis.set_ylabel("Thermal Pressure [kPa]")
    axis.minorticks_on()
    axis.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
    axis.set_title("Plasma Thermal Pressure Profiles")
    axis.grid(True, linestyle="--", alpha=0.5)
    axis.set_xlim(0, 1.025)
    axis.set_ylim(bottom=0)
    axis.legend()

    textstr_pressure = "\n".join((
        (
            r"$p_0$:"
            rf" {mfile.get('pres_plasma_thermal_on_axis', scan=scan) / 1000:,.3f} kPa"
            r"$\hspace{2} \frac{p_0}{\langle p_{\text{total}}"
            r" \rangle_\text{V}}$:"
            rf" {mfile.get('f_pres_plasma_thermal_on_axis_vol_avg', scan=scan):,.3f}"
        ),
        (
            r"$\langle p_{\text{total}} \rangle_\text{V}$:"
            rf" {mfile.get('pres_plasma_thermal_vol_avg', scan=scan) / 1000:,.3f} kPa"
        ),
    ))

    draw_text(
        axis,
        0.5,
        1.2,
        textstr_pressure,
        transform=axis.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="center",
        bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5},
    )

    if (
        int(mfile.get("i_plasma_pedestal", scan=scan))
        == PlasmaProfileShapeType.PEDESTAL_PROFILE
    ):
        textstr_pressure_pedestal = "\n".join((
            (
                r"$p_{\text{ped}}$:"
                rf" {mfile.get('pres_plasma_pedestal_thermal', scan=scan) / 1000:,.3f} kPa"  # noqa: E501
            ),
            (
                r"$p_{\text{sep}}$:"
                rf" {mfile.get('pres_plasma_separatrix_thermal', scan=scan) / 1000:,.3f} kPa"  # noqa: E501
            ),
        ))

        draw_text(
            axis,
            0.9,
            1.2,
            textstr_pressure_pedestal,
            transform=axis.transAxes,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="center",
            bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5},
        )


def plot_plasma_poloidal_pressure_contours(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot plasma poloidal pressure contours inside the plasma boundary.

    This function visualizes the poloidal pressure distribution inside the plasma
    boundary
    by interpolating the pressure profile onto a grid defined by the plasma geometry.
    The pressure is shown as filled contours, with the plasma boundary overlaid.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        Matplotlib axis object to plot on.
    mfile : mfile: MFile
        MFILE data object containing plasma and geometry data.
    scan : int
        Scan number to use for extracting data.
    """
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    # Get pressure profile (function of normalised radius rho, 0..1)
    pres_plasma_electron_profile = [
        mfile.get(f"pres_plasma_electron_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    pres_plasma_profile_ion = [
        mfile.get(f"pres_plasma_ion_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]

    # Convert pressure to kPa
    pres_plasma_electron_profile_kpa = [p / 1000.0 for p in pres_plasma_electron_profile]
    pres_plasma_profile_ion_kpa = [p / 1000.0 for p in pres_plasma_profile_ion]
    pres_plasma_profile = [
        e + i
        for e, i in zip(
            pres_plasma_electron_profile_kpa,
            pres_plasma_profile_ion_kpa,
            strict=False,
        )
    ]

    pressure_grid, r_grid, z_grid = interp1d_profile(pres_plasma_profile, mfile, scan)

    # Mask points outside the plasma boundary (optional, but grid is inside by
    # construction)
    # Plot filled contour
    c = axis.contourf(r_grid, -z_grid, pressure_grid, levels=50, cmap="plasma")

    # Add colorbar for pressure (now in kPa)
    # You can control the location using the 'location' argument ('left', 'right', 'top',
    # 'bottom')
    # For more control, use 'ax' or 'fraction', 'pad', etc.
    # Example: location="right", pad=0.05, fraction=0.05
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)

    axis.figure.colorbar(
        c,
        ax=axis,
        label="Pressure [kPa]",
        location="left",
        anchor=(-0.25, 0.5),
    )

    axis.set_aspect("equal")
    axis.set_xlabel("R [m]")
    axis.set_xlim(rmajor - 1.2 * rminor, rmajor + 1.2 * rminor)
    axis.set_ylim(
        -1.2 * rminor * mfile.get("kappa", scan=scan),
        1.2 * mfile.get("kappa", scan=scan) * rminor,
    )
    axis.set_ylabel("Z [m]")
    axis.set_title("Plasma Poloidal Pressure Contours")
    axis.plot(
        rmajor,
        0,
        marker="o",
        color="red",
        markersize=6,
        markeredgecolor="black",
        zorder=100,
    )


def plot_fusion_rate_contours(fig1, fig2, mfile: MFile, scan: int):
    """Plot fusion rate density contours"""
    rmajor = mfile.get("rmajor", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    kappa = mfile.get("kappa", scan=scan)
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    def fusrat(name):
        fusrat_dat = [
            mfile.get(f"fusrat_plasma_{name}_profile{i}", scan=scan)
            for i in range(n_plasma_profile_elements)
        ]
        return interp1d_profile(fusrat_dat, mfile, scan)

    dt_grid, _r_grid, _z_grid = fusrat("dt")
    dd_triton_grid, _r_grid, _z_grid = fusrat(" dd_triton ")
    dd_helion_grid, _r_grid, _z_grid = fusrat(" dd_helion ")
    dhe3_grid, r_grid, z_grid = fusrat(" dhe3")

    dt_axes = fig1.add_subplot(121, aspect="equal")
    dd_triton_axes = fig1.add_subplot(122, aspect="equal")
    dd_helion_axes = fig2.add_subplot(121, aspect="equal")
    dhe3_axes = fig2.add_subplot(122, aspect="equal")

    dt_axes.set_title("D+T -> 4He + n Fusion Rate Density Contours")
    reaction_plot_grid(rminor, rmajor, kappa, r_grid, z_grid, dt_grid, dt_axes)

    dd_triton_axes.set_title("D+D -> T + p Fusion Rate Density Contours")
    reaction_plot_grid(
        rminor, rmajor, kappa, r_grid, z_grid, dd_triton_grid, dd_triton_axes
    )
    dd_helion_axes.set_title("D+D -> 3He + n Fusion Rate Density Contours")
    reaction_plot_grid(
        rminor, rmajor, kappa, r_grid, z_grid, dd_helion_grid, dd_helion_axes
    )
    dhe3_axes.set_title("D+3He -> 4He + n Fusion Rate Density Contours")
    reaction_plot_grid(rminor, rmajor, kappa, r_grid, z_grid, dhe3_grid, dhe3_axes)


def plot_beta_profiles(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot the beta profiles on the given axis"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    beta_plasma_toroidal_profile = [
        mfile.get(f"beta_thermal_toroidal_profile{i}", scan=scan)
        for i in range(2 * n_plasma_profile_elements)
    ]

    axis.plot(
        np.linspace(-1, 1, 2 * n_plasma_profile_elements),
        beta_plasma_toroidal_profile,
        color="blue",
        label="$\\beta_t$",
    )

    axis.axhline(
        mfile.get("beta_thermal_toroidal_vol_avg", scan=scan),
        color="blue",
        linestyle="--",
        linewidth=1.0,
        label="$\\langle \\beta_t \\rangle_{\\text{V}}$",
    )

    axis.set_xlabel("$\\rho$ [r/a]")
    axis.set_ylabel("$\\beta$")
    axis.minorticks_on()
    axis.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
    axis.set_title("Thermal Beta Profiles")
    axis.legend()
    axis.axvline(x=0, color="black", linestyle="--", linewidth=1)
    axis.grid(True, linestyle="--", alpha=0.5)
    axis.set_ylim(bottom=0.0)


def plot_plasma_effective_charge_profile(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot plasma effective charge profile"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    n_charge_plasma_effective_vol_avg = mfile.get(
        "n_charge_plasma_effective_vol_avg", scan=scan
    )

    n_charge_plasma_effective_profile = [
        mfile.get(f"n_charge_plasma_effective_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]

    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        n_charge_plasma_effective_profile,
    )

    axis.hlines(
        n_charge_plasma_effective_vol_avg,
        xmin=0,
        xmax=1,
        colors="red",
        linestyles="--",
        label=(
            "Volume-Averaged $Z_{\\text{eff}}$ ="
            f" {n_charge_plasma_effective_vol_avg:.2f}"
        ),
    )

    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_ylabel("Effective Charge ($Z_{\\text{eff}}$)")
    axis.set_title("Plasma Effective Charge Profile")
    axis.minorticks_on()
    axis.set_xlim(0, 1.025)
    axis.grid(which="both", linestyle="--", alpha=0.5)
    axis.legend()


def plot_plasma_thermal_energy_profiles(axis, m_file: MFile, scan: int):
    """Function to plot plasma thermal energy profiles on the given axis.

    Parameters
    ----------
    axis :
        Matplotlib axis to plot on
    m_file :
        MFILE
    scan :
        scan to read from MFILE
    """
    n_plasma_profile_elements = int(m_file.get("n_plasma_profile_elements", scan=scan))

    eden_plasma_electrons_thermal_profile_mj = [
        m_file.get(f"eden_plasma_electrons_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]
    eden_plasma_ions_thermal_profile_mj = [
        m_file.get(f"eden_plasma_ions_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]
    eden_plasma_thermal_profile_mj = [
        m_file.get(f"eden_plasma_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]
    e_plasma_electrons_thermal_profile_mj = [
        m_file.get(f"e_plasma_electrons_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]

    e_plasma_ions_thermal_profile_mj = [
        m_file.get(f"e_plasma_ions_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]
    e_plasma_thermal_profile_mj = [
        m_file.get(f"e_plasma_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]

    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        e_plasma_electrons_thermal_profile_mj,
        label="$W_{\\text{e}}$",
        color="tab:blue",
        linestyle=":",
    )
    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        e_plasma_ions_thermal_profile_mj,
        label="$W_{\\text{i}}$",
        color="tab:blue",
        linestyle="--",
    )
    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        e_plasma_thermal_profile_mj,
        label="$W_{\\text{total}}$",
        color="tab:blue",
        linestyle="-",
    )

    total_thermal_energy_max_index = int(np.argmax(e_plasma_thermal_profile_mj))
    total_thermal_energy_max_rho = np.linspace(0, 1, n_plasma_profile_elements)[
        total_thermal_energy_max_index
    ]
    total_thermal_energy_max = e_plasma_thermal_profile_mj[
        total_thermal_energy_max_index
    ]
    axis.axvline(
        total_thermal_energy_max_rho,
        color="tab:red",
        alpha=0.7,
        label="$W_{\\text{total, peak}}$",
    )
    axis.axhline(
        total_thermal_energy_max,
        color="tab:red",
        alpha=0.7,
    )

    density_axis = axis.twinx()
    density_axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        eden_plasma_electrons_thermal_profile_mj,
        label="$W_{\\text{density, e}}$",
        color="tab:orange",
        linestyle=":",
    )
    density_axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        eden_plasma_ions_thermal_profile_mj,
        label="$W_{\\text{density, i}}$",
        color="tab:orange",
        linestyle="--",
    )
    density_axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        eden_plasma_thermal_profile_mj,
        label="$W_{\\text{density, total}}$",
        color="tab:orange",
        linestyle="-",
    )

    axis.grid(True, alpha=0.3)
    axis.minorticks_on()
    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_xlim(left=0.0, right=1.0)
    axis.set_ylabel(
        "Thermal Energy [MJ]",
        color="tab:blue",
    )
    axis.tick_params(axis="y", colors="tab:blue")
    density_axis.set_ylabel(
        "Thermal Energy Density [MJ/m$^3$]",
        color="tab:orange",
    )
    density_axis.tick_params(axis="y", colors="tab:orange")
    handles, labels = axis.get_legend_handles_labels()
    density_handles, density_labels = density_axis.get_legend_handles_labels()
    axis.legend(handles + density_handles, labels + density_labels)


def plot_cumulative_plasma_thermal_energy_profiles(axis, m_file: MFile, scan: int):
    """Function to plot the cumulative plasma thermal energy profiles on the given axis.

    Parameters
    ----------
    axis :
        Matplotlib axis to plot on
    m_file :
        MFILE
    scan :
        scan to read from MFILE
    """
    n_plasma_profile_elements = int(m_file.get("n_plasma_profile_elements", scan=scan))
    e_plasma_thermal_total_mj = m_file.get("e_plasma_thermal_total", scan=scan) / 1e6
    e_plasma_electrons_thermal_profile_mj = [
        m_file.get(f"e_plasma_electrons_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]

    e_plasma_ions_thermal_profile_mj = [
        m_file.get(f"e_plasma_ions_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]
    e_plasma_thermal_profile_mj = [
        m_file.get(f"e_plasma_thermal_profile{i}", scan=scan) / 1e6
        for i in range(n_plasma_profile_elements)
    ]

    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        np.cumsum(e_plasma_electrons_thermal_profile_mj),
        label="$\\Sigma W_{\\text{e}}$",
        color="tab:blue",
        linestyle=":",
    )
    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        np.cumsum(e_plasma_ions_thermal_profile_mj),
        label="$\\Sigma W_{\\text{i}}$",
        color="tab:blue",
        linestyle="--",
    )
    axis.plot(
        np.linspace(0, 1, n_plasma_profile_elements),
        np.cumsum(e_plasma_thermal_profile_mj),
        label="$\\Sigma W_{\\text{total}}$",
        color="tab:blue",
        linestyle="-",
    )
    axis.axhline(
        y=e_plasma_thermal_total_mj,
        label="$W_{\\text{thermal,total}}$",
        color="tab:red",
        linestyle="--",
    )
    cumulative_thermal_energy_mj = np.cumsum(e_plasma_thermal_profile_mj)
    half_thermal_energy_mj = 0.5 * e_plasma_thermal_total_mj
    half_thermal_energy_position = np.interp(
        half_thermal_energy_mj,
        cumulative_thermal_energy_mj,
        np.linspace(0, 1, n_plasma_profile_elements),
    )
    axis.axhline(
        y=half_thermal_energy_mj,
        label="$50\\%\\ W_{\\text{thermal,total}}$",
        color="tab:green",
        linestyle=":",
    )
    axis.axvline(
        x=half_thermal_energy_position,
        color="tab:green",
        linestyle=":",
    )

    axis.legend()
    axis.set_title("Thermal Energy Profiles and Cumulative Distribution")
    axis.grid(True, alpha=0.3)
    axis.minorticks_on()
    axis.tick_params(axis="x", labelbottom=False)
    axis.set_xlim(left=0.0, right=1.0)
    axis.set_ylabel(
        "Cumulative Thermal Energy [MJ]",
    )


__all__ = [
    "plot_beta_profiles",
    "plot_cumulative_plasma_thermal_energy_profiles",
    "plot_fusion_rate_contours",
    "plot_fusion_rate_profiles",
    "plot_jprofile",
    "plot_n_profiles",
    "plot_plasma_effective_charge_profile",
    "plot_plasma_poloidal_pressure_contours",
    "plot_plasma_pressure_profiles",
    "plot_plasma_thermal_energy_profiles",
    "plot_qprofile",
    "plot_t_profiles",
]
