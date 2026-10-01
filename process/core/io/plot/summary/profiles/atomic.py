"""Profiles functions for PROCESS summary plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.data_structure.impurity_radiation_variables import (
    N_IMPURITIES,
    ImpurityRadiationData,
)

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_ion_charge_profile(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot ion charge profile"""
    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    # find impurity densities
    imp_frac = np.array([
        mfile.get("f_nd_impurity_electrons(01)", scan=scan),
        mfile.get("f_nd_impurity_electrons(02)", scan=scan),
        mfile.get("f_nd_impurity_electrons(03)", scan=scan),
        mfile.get("f_nd_impurity_electrons(04)", scan=scan),
        mfile.get("f_nd_impurity_electrons(05)", scan=scan),
        mfile.get("f_nd_impurity_electrons(06)", scan=scan),
        mfile.get("f_nd_impurity_electrons(07)", scan=scan),
        mfile.get("f_nd_impurity_electrons(08)", scan=scan),
        mfile.get("f_nd_impurity_electrons(09)", scan=scan),
        mfile.get("f_nd_impurity_electrons(10)", scan=scan),
        mfile.get("f_nd_impurity_electrons(11)", scan=scan),
        mfile.get("f_nd_impurity_electrons(12)", scan=scan),
        mfile.get("f_nd_impurity_electrons(13)", scan=scan),
        mfile.get("f_nd_impurity_electrons(14)", scan=scan),
    ])

    n_charge_plasma_profile = []
    impurity_data = ImpurityRadiationData()
    for imp in range(N_IMPURITIES):
        if imp_frac[imp] > 1.0e-30:
            profile = [
                mfile.get(f"n_charge_plasma_profile{imp}_{i}", scan=scan)
                for i in range(n_plasma_profile_elements)
            ]
            n_charge_plasma_profile.append(profile)
            z_max = impurity_data.imp_full_ion_charge[imp]
            # Calculate relative ionisation state as percent of full ionisation
            rel_ion_state = [
                100.0 * (val / z_max if z_max > 0 else 0) for val in profile
            ]
            avg_ionisation = np.mean(rel_ion_state)
            axis.plot(
                np.linspace(0, 1, n_plasma_profile_elements),
                rel_ion_state,
                label=f"{impurity_data.imp_label[imp].replace('_', '')} (Z={z_max}): "
                f"avg {avg_ionisation:.1f}%",
            )
    axis.set_ylabel("Relative Ionisation State [% of $Z$]")
    axis.legend()
    axis.set_xlim(0, 1.025)
    axis.set_xlabel(r"$\rho \quad [r/a]$")
    axis.set_title("Impurity Ion Charge State Profiles")
    axis.minorticks_on()
    axis.grid(which="both", linestyle="--", alpha=0.5)


def plot_debye_length_profile(axis: plt.Axes, mfile_data: MFile, scan: int):
    """Plot the Debye length profile on the given axis."""
    len_plasma_debye_electron_profile = [
        mfile_data.data[f"len_plasma_debye_electron_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    # Convert to micrometres (1e-6 m)
    len_plasma_debye_electron_profile_um = [
        length * 1e6 for length in len_plasma_debye_electron_profile
    ]

    axis.plot(
        np.linspace(0, 1, len(len_plasma_debye_electron_profile_um)),
        len_plasma_debye_electron_profile_um,
        color="blue",
        linestyle="-",
        label=r"$\lambda_{Debye,e}$",
    )

    axis.set_ylabel(r"Debye Length [$\mu$m]")

    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.set_xlim(0, 1.025)
    axis.minorticks_on()
    axis.legend()


def plot_velocity_profile(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the electron thermal velocity profile on the given axis."""
    vel_plasma_electron_profile = [
        mfile_data.data[f"vel_plasma_electron_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]
    vel_plasma_deuteron_profile = [
        mfile_data.data[f"vel_plasma_deuteron_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]
    vel_plasma_triton_profile = [
        mfile_data.data[f"vel_plasma_triton_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]
    vel_plasma_alpha_thermal_profile = [
        mfile_data.data[f"vel_plasma_alpha_thermal_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    vel_plasma_alpha_birth = mfile_data.data["vel_plasma_alpha_birth"].get_scan(scan)

    axis.plot(
        np.linspace(0, 1, len(vel_plasma_electron_profile)),
        vel_plasma_electron_profile,
        color="blue",
        linestyle="-",
        label=r"$v_{e}$",
    )
    axis.plot(
        np.linspace(0, 1, len(vel_plasma_deuteron_profile)),
        vel_plasma_deuteron_profile,
        color="pink",
        linestyle="-",
        label=r"$v_{D}$",
    )
    axis.plot(
        np.linspace(0, 1, len(vel_plasma_triton_profile)),
        vel_plasma_triton_profile,
        color="green",
        linestyle="-",
        label=r"$v_{T}$",
    )
    axis.plot(
        np.linspace(0, 1, len(vel_plasma_alpha_thermal_profile)),
        vel_plasma_alpha_thermal_profile,
        color="red",
        linestyle="-",
        label=r"$v_{\alpha,thermal}$",
    )
    axis.axhline(
        vel_plasma_alpha_birth,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=r"$v_{\alpha,birth}$",
    )

    axis.set_yscale("log")
    axis.set_ylabel("Velocity [m/s]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.set_xlim(0, 1.025)
    axis.minorticks_on()
    axis.legend()


def plot_electron_frequency_profile(
    axis: plt.Axes, mfile_data: MFile, scan: int
) -> None:
    """Plot the electron thermal frequency profile on the given axis."""
    freq_plasma_electron_profile = [
        mfile_data.data[f"freq_plasma_electron_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]
    freq_plasma_larmor_toroidal_electron_profile = [
        mfile_data.data[f"freq_plasma_larmor_toroidal_electron_profile{i}"].get_scan(
            scan
        )
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    freq_plasma_upper_hybrid_electron_profile = [
        mfile_data.data[f"freq_plasma_upper_hybrid_profile{i}"].get_scan(scan)
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_larmor_toroidal_electron_profile)),
        np.array(freq_plasma_larmor_toroidal_electron_profile) / 1e9,
        color="red",
        linestyle="-",
        label=r"$f_{Larmor,toroidal,e}$ | Fundamental",
    )

    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_larmor_toroidal_electron_profile)),
        2 * np.array(freq_plasma_larmor_toroidal_electron_profile) / 1e9,
        color="red",
        linestyle="--",
        label=r"$f_{Larmor,toroidal,e}$ | 2nd harmonic",
    )

    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_larmor_toroidal_electron_profile)),
        3 * np.array(freq_plasma_larmor_toroidal_electron_profile) / 1e9,
        color="red",
        linestyle=":",
        label=r"$f_{Larmor,toroidal,e}$ | 3rd harmonic",
    )

    x = np.linspace(0, 1, len(freq_plasma_electron_profile))
    y = np.array(freq_plasma_electron_profile) / 1e9
    # original curve
    axis.plot(
        x,
        y,
        color="blue",
        linestyle="-",
        label=r"$\omega_{p,e}$ | Plasma Frequency",
    )
    # mirrored across the y-axis (drawn at negative rho)
    axis.plot(-x, y, color="blue", linestyle="-", label="_nolegend_")

    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_upper_hybrid_electron_profile)),
        np.array(freq_plasma_upper_hybrid_electron_profile) / 1e9,
        color="purple",
        linestyle="-",
        label=r"$\omega_{UH,e}$ | Upper Hybrid",
    )

    axis.set_xlim(-1.025, 1.025)
    axis.set_ylim(None, max(freq_plasma_larmor_toroidal_electron_profile) / 1e9 * 1.6)

    axis.set_xlabel("$\\rho$ [r/a]")
    axis.set_ylabel("Frequency [GHz]")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)

    # Add secondary x-axis showing radius in metres below the primary axis
    ax2 = axis.twiny()
    rmajor = mfile_data.get("rmajor", scan=scan)
    rminor = mfile_data.get("rminor", scan=scan)

    # Convert normalized radius to actual radius
    # rho ranges from -1 to 1, which corresponds to r = rmajor - rminor to rmajor +
    # rminor
    rho_ticks = np.array([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
    r_ticks = rmajor + rho_ticks * rminor

    ax2.set_xticks(rho_ticks)
    ax2.set_xticklabels([f"{r:.2f}" for r in r_ticks])
    ax2.set_xlabel("Radius [m]")
    ax2.minorticks_on()
    ax2.set_xlim(axis.get_xlim())

    # Move secondary axis to the bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("outward", 30))

    axis.legend()


def plot_ion_frequency_profile(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the ion thermal frequency profile on the given axis."""
    freq_plasma_larmor_toroidal_deuteron_profile = [
        mfile_data.data[f"freq_plasma_larmor_toroidal_deuteron_profile{i}"].get_scan(
            scan
        )
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    freq_plasma_larmor_toroidal_triton_profile = [
        mfile_data.data[f"freq_plasma_larmor_toroidal_triton_profile{i}"].get_scan(scan)
        for i in range(
            2 * int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan))
        )
    ]

    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_larmor_toroidal_deuteron_profile)),
        np.array(freq_plasma_larmor_toroidal_deuteron_profile) / 1e6,
        color="red",
        linestyle="-",
        label=r"$f_{Larmor,toroidal,D}$",
    )
    axis.plot(
        np.linspace(-1, 1, len(freq_plasma_larmor_toroidal_triton_profile)),
        np.array(freq_plasma_larmor_toroidal_triton_profile) / 1e6,
        color="green",
        linestyle="-",
        label=r"$f_{Larmor,toroidal,T}$",
    )

    axis.set_ylabel("Frequency [MHz]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_collision_time_profile(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the plasma collision times on the given axis."""
    t_plasma_electron_electron_collision_profile = [
        mfile_data.data[f"t_plasma_electron_electron_collision_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    t_plasma_electron_deuteron_collision_profile = [
        mfile_data.data[f"t_plasma_electron_deuteron_collision_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    t_plasma_electron_triton_collision_profile = [
        mfile_data.data[f"t_plasma_electron_triton_collision_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    t_plasma_electron_alpha_thermal_collision_profile = [
        mfile_data.data[
            f"t_plasma_electron_alpha_thermal_collision_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(t_plasma_electron_electron_collision_profile)),
        t_plasma_electron_electron_collision_profile,
        color="blue",
        linestyle="-",
        label=r"$\tau_{e-e}$",
    )

    axis.plot(
        np.linspace(0, 1, len(t_plasma_electron_deuteron_collision_profile)),
        t_plasma_electron_deuteron_collision_profile,
        color="pink",
        linestyle="-",
        label=r"$\tau_{e-D}$",
    )

    axis.plot(
        np.linspace(0, 1, len(t_plasma_electron_triton_collision_profile)),
        t_plasma_electron_triton_collision_profile,
        color="green",
        linestyle="-",
        label=r"$\tau_{e-T}$",
    )

    axis.plot(
        np.linspace(0, 1, len(t_plasma_electron_alpha_thermal_collision_profile)),
        t_plasma_electron_alpha_thermal_collision_profile,
        color="red",
        linestyle="-",
        label=r"$\tau_{e-\alpha,thermal}$",
    )

    axis.set_yscale("log")
    axis.set_ylabel("Collision Time [s]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_collision_frequency_profile(
    axis: plt.Axes, mfile_data: MFile, scan: int
) -> None:
    """Plot the plasma collision frequencies on the given axis."""
    freq_plasma_electron_electron_collision_profile = [
        mfile_data.data[f"freq_plasma_electron_electron_collision_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    freq_plasma_electron_deuteron_collision_profile = [
        mfile_data.data[f"freq_plasma_electron_deuteron_collision_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    freq_plasma_electron_triton_collision_profile = [
        mfile_data.data[f"freq_plasma_electron_triton_collision_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    freq_plasma_electron_alpha_thermal_collision_profile = [
        mfile_data.data[
            f"freq_plasma_electron_alpha_thermal_collision_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(freq_plasma_electron_electron_collision_profile)),
        freq_plasma_electron_electron_collision_profile,
        color="blue",
        linestyle="-",
        label=r"$\nu_{e-e}$",
    )

    axis.plot(
        np.linspace(0, 1, len(freq_plasma_electron_deuteron_collision_profile)),
        freq_plasma_electron_deuteron_collision_profile,
        color="pink",
        linestyle="-",
        label=r"$\nu_{e-D}$",
    )

    axis.plot(
        np.linspace(0, 1, len(freq_plasma_electron_triton_collision_profile)),
        freq_plasma_electron_triton_collision_profile,
        color="green",
        linestyle="-",
        label=r"$\nu_{e-T}$",
    )

    axis.plot(
        np.linspace(0, 1, len(freq_plasma_electron_alpha_thermal_collision_profile)),
        freq_plasma_electron_alpha_thermal_collision_profile,
        color="red",
        linestyle="-",
        label=r"$\nu_{e-\alpha,thermal}$",
    )
    axis.set_yscale("log")
    axis.set_ylabel("Collision Frequency [Hz]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_mean_free_path_profile(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the plasma mean free path on the given axis."""
    len_plasma_electron_electron_mean_free_path_profile = [
        mfile_data.data[
            f"len_plasma_electron_electron_mean_free_path_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    len_plasma_electron_deuteron_mean_free_path_profile = [
        mfile_data.data[
            f"len_plasma_electron_deuteron_mean_free_path_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    len_plasma_electron_triton_mean_free_path_profile = [
        mfile_data.data[
            f"len_plasma_electron_triton_mean_free_path_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    len_plasma_electron_alpha_thermal_mean_free_path_profile = [
        mfile_data.data[
            f"len_plasma_electron_alpha_thermal_mean_free_path_profile{i}"
        ].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(len_plasma_electron_electron_mean_free_path_profile)),
        len_plasma_electron_electron_mean_free_path_profile,
        color="blue",
        linestyle="-",
        label=r"$\lambda_{mfp,e-e}$",
    )

    axis.plot(
        np.linspace(0, 1, len(len_plasma_electron_deuteron_mean_free_path_profile)),
        len_plasma_electron_deuteron_mean_free_path_profile,
        color="pink",
        linestyle="-",
        label=r"$\lambda_{mfp,e-D}$",
    )

    axis.plot(
        np.linspace(0, 1, len(len_plasma_electron_triton_mean_free_path_profile)),
        len_plasma_electron_triton_mean_free_path_profile,
        color="green",
        linestyle="-",
        label=r"$\lambda_{mfp,e-T}$",
    )
    axis.plot(
        np.linspace(0, 1, len(len_plasma_electron_alpha_thermal_mean_free_path_profile)),
        len_plasma_electron_alpha_thermal_mean_free_path_profile,
        color="red",
        linestyle="-",
        label=r"$\lambda_{mfp,e-\alpha,thermal}$",
    )
    axis.set_yscale("log")
    axis.set_ylabel("Mean Free Path [m]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_ion_slowing_down_time_profile(
    axis: plt.Axes, mfile_data: MFile, scan: int
) -> None:
    """Plot the plasma Spitzer slowing down time on the given axis."""
    t_plasma_electron_alpha_spitzer_slow_profile = [
        mfile_data.data[f"t_plasma_electron_alpha_spitzer_slow_profile{i}"].get_scan(
            scan
        )
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(t_plasma_electron_alpha_spitzer_slow_profile)),
        t_plasma_electron_alpha_spitzer_slow_profile,
        color="red",
        linestyle="-",
        label=r"$\tau_{e-\alpha,Spitzer}$",
    )

    axis.set_yscale("log")
    axis.set_ylabel("Spitzer Slowing Down Time [s]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


def plot_resistivity_profile(axis: plt.Axes, mfile_data: MFile, scan: int) -> None:
    """Plot the plasma resistivity on the given axis."""
    res_plasma_fuel_spitzer_profile = [
        mfile_data.data[f"res_plasma_fuel_spitzer_profile{i}"].get_scan(scan)
        for i in range(int(mfile_data.data["n_plasma_profile_elements"].get_scan(scan)))
    ]

    axis.plot(
        np.linspace(0, 1, len(res_plasma_fuel_spitzer_profile)),
        res_plasma_fuel_spitzer_profile,
        color="red",
        linestyle="-",
        label=r"$\eta_{Spitzer-fuel}$",
    )

    axis.set_yscale("log")
    axis.set_ylabel("Resistivity [Ohm m]")
    axis.set_xlabel("$\\rho \\ [r/a]$")
    axis.grid(True, which="both", linestyle="--", alpha=0.5)
    axis.minorticks_on()
    axis.legend()


__all__ = [
    "plot_collision_frequency_profile",
    "plot_collision_time_profile",
    "plot_debye_length_profile",
    "plot_electron_frequency_profile",
    "plot_ion_charge_profile",
    "plot_ion_frequency_profile",
    "plot_ion_slowing_down_time_profile",
    "plot_mean_free_path_profile",
    "plot_resistivity_profile",
    "plot_velocity_profile",
]
