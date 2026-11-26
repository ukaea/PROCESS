"""Plotting utilities for bootstrap current summary plots using MFile data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.plasma.physics import plot_plasma
from process.core.io.plot.summary.profiles.misc import profiles_with_pedestal
from process.models.physics.bootstrap_current import SauterBootstrapCurrent

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_sauter_bootstrap_summary(fig: plt.Figure, mfile: MFile, scan: int):
    """Plot the Sauter bootstrap current summary for a given MFile and scan."""
    rho, nd_electron_profile, temp_electron_profile_kev = profiles_with_pedestal(
        mfile=mfile, scan=scan
    )
    sauter = SauterBootstrapCurrent()
    a_plasma_poloidal = mfile.get("a_plasma_poloidal", scan=scan)
    rminor = mfile.get("rminor", scan=scan)
    rmajor = mfile.get("rmajor", scan=scan)

    sauter_params = sauter.bootstrap_fraction_sauter(
        a_plasma_poloidal=a_plasma_poloidal,
        rminor=rminor,
        rmajor=rmajor,
        triang=mfile.get("triang", scan=scan),
        n_charge_plasma_effective_vol_avg=mfile.get(
            "n_charge_plasma_effective_vol_avg", scan=scan
        ),
        nd_plasma_electron_profile=nd_electron_profile,
        roa=rho,
        temp_plasma_electron_profile_kev=temp_electron_profile_kev,
        temp_plasma_ion_profile_kev=temp_electron_profile_kev
        * mfile.get("f_temp_plasma_ion_electron", scan=scan),
        nd_plasma_ions_total_vol_avg=mfile.get(
            "nd_plasma_ions_total_vol_avg", scan=scan
        ),
        nd_plasma_electrons_vol_avg=mfile.get("nd_plasma_electrons_vol_avg", scan=scan),
        q0=mfile.get("q0", scan=scan),
        q95=mfile.get("q95", scan=scan),
        m_ions_total_amu=mfile.get("m_ions_total_amu", scan=scan),
        f_plasma_fuel_helium3=mfile.get("f_plasma_fuel_helium3", scan=scan),
        b_plasma_toroidal_on_axis=mfile.get("b_plasma_toroidal_on_axis", scan=scan),
        plasma_current=mfile.get("plasma_current", scan=scan),
    )
    grid = fig.add_gridspec(3, 3)
    j_profile = fig.add_subplot(grid[0, 0])

    j_profile.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.j_bootstrap,
        label="Sauter Bootstrap Current Density",
    )

    axes_trapped = fig.add_subplot(grid[0, 1])

    axes_trapped.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.f_trapped_particles,
        label="$f_{trapped}$",
    )

    axes_coefficients = fig.add_subplot(grid[2, 0])
    axes_coefficients.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.l31_coefficient,
        label=r"$\mathcal{L}_{31}$ Coefficient",
    )
    axes_coefficients.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.l32_coefficient,
        label=r"$\mathcal{L}_{32}$ Coefficient",
    )
    axes_coefficients.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.l34_coefficient,
        label=r"$\mathcal{L}_{34}$ Coefficient",
    )
    axes_coefficients.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.alpha,
        label="$\\alpha$ Coefficient",
    )
    axes_coefficients.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.alpha_0,
        label="$\\alpha_0$ Coefficient",
    )

    axes_collisionality = fig.add_subplot(grid[1, 1])
    axes_collisionality.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.normalised_electron_collisionality,
        label="$\\nu_{e}^*$",
    )
    axes_collisionality.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.normalised_ion_collisionality,
        label="$\\nu_{i}^*$",
    )

    axes_aspect_ratio = fig.add_subplot(grid[1, 0])
    axes_aspect_ratio.plot(
        rho,
        sauter_params.inverse_aspect,
        label="$\\epsilon$",
    )
    axes_aspect_ratio.set_xlabel("Normalized Radius")
    axes_aspect_ratio.set_ylabel("Inverse Aspect Ratio")
    axes_aspect_ratio.grid(True)
    axes_aspect_ratio.legend()

    axes_aspect_ratio.plot(
        rho,
        sauter_params.sqrt_inverse_aspect,
        label="$\\sqrt{\\epsilon}$",
    )
    axes_aspect_ratio.legend()

    axes_collisionality.set_xlabel("Normalized Radius")
    axes_collisionality.set_ylabel("Collisionality")
    axes_collisionality.grid(True)
    axes_collisionality.legend()

    axes_coefficients.legend()
    axes_coefficients.set_xlabel("Normalized Radius")
    axes_coefficients.set_ylabel("Coefficient Value")
    axes_coefficients.grid(True)
    j_profile.set_xlabel("Normalized Radius")
    j_profile.set_ylabel("Bootstrap Current Density")
    j_profile.grid(True)
    j_profile.legend()
    axes_trapped.set_xlabel("Normalized Radius")
    axes_trapped.set_ylabel("Trapped Fraction")
    axes_trapped.grid(True)
    axes_trapped.legend()
    q_profile = fig.add_subplot(grid[2, 1])
    q_profile.plot(
        rho,
        sauter_params.q_profile,
        label="Assumed $q$ profile",
    )
    q_profile.set_xlabel("Normalized Radius")
    q_profile.set_ylabel("Safety Factor")
    q_profile.grid(True)
    q_profile.legend()

    geometry = fig.add_subplot(grid[:, 2], aspect="equal")
    plot_plasma(geometry, mfile, scan, colour_scheme=1)
    for index, boundary in enumerate(geometry.lines):
        boundary.set_label("PROCESS plasma outline" if index == 0 else "_nolegend_")
        boundary.set_linestyle("--")
        boundary.set_zorder(1)
    for fill in (*geometry.collections, *geometry.patches):
        fill.set_alpha(0.25)
        fill.set_zorder(0)
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    kappa_effective = a_plasma_poloidal / (np.pi * rminor**2)
    for radius_fraction in (0.25, 0.5, 0.75, 1.0):
        radius = radius_fraction * rminor
        geometry.plot(
            rmajor + radius * np.cos(theta),
            kappa_effective * radius * np.sin(theta),
            label=f"Assumed $r/a={radius_fraction:g}$",
        )
    geometry.plot(rmajor, 0, "k+", label="Centre")
    geometry.set_xlabel("$R$ [m]")
    geometry.set_ylabel("$Z$ [m]")
    geometry.set_title(
        "Assumed flux surfaces in plasma\n"
        rf"Effective elongation $\kappa_{{\rm eff}} = {kappa_effective:.3f}$"
    )
    geometry.grid(True)
    geometry.legend(fontsize="small", loc="best")
    fig.suptitle("Sauter Model Graphs")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for axis in (
        j_profile,
        axes_trapped,
        axes_aspect_ratio,
        axes_collisionality,
        axes_coefficients,
        q_profile,
    ):
        axis.set_position(axis.get_position().translated(0, 0.015))
    geometry.set_position(geometry.get_position().translated(-0.008, 0))
