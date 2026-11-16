"""Plotting utilities for bootstrap current summary plots using MFile data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

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

    sauter_params = sauter.bootstrap_fraction_sauter(
        a_plasma_poloidal=mfile.get("a_plasma_poloidal"),
        rminor=mfile.get("rminor"),
        rmajor=mfile.get("rmajor"),
        triang=mfile.get("triang"),
        n_charge_plasma_effective_vol_avg=mfile.get("n_charge_plasma_effective_vol_avg"),
        nd_plasma_electron_profile=nd_electron_profile,
        roa=rho,
        temp_plasma_electron_profile_kev=temp_electron_profile_kev,
        temp_plasma_ion_profile_kev=temp_electron_profile_kev,
        nd_plasma_ions_total_vol_avg=mfile.get("nd_plasma_ions_total_vol_avg"),
        nd_plasma_electrons_vol_avg=mfile.get("nd_plasma_electrons_vol_avg"),
        q0=mfile.get("q0"),
        q95=mfile.get("q95"),
        m_ions_total_amu=mfile.get("m_ions_total_amu"),
        f_plasma_fuel_helium3=mfile.get("f_plasma_fuel_helium3"),
        b_plasma_toroidal_on_axis=mfile.get("b_plasma_toroidal_on_axis"),
        plasma_current=mfile.get("plasma_current"),
    )
    j_profile = fig.add_subplot(221)

    j_profile.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.j_bootstrap,
        label="Sauter Bootstrap Current",
    )

    axes_trapped = fig.add_subplot(322)

    axes_trapped.plot(
        np.linspace(0, 1, len(rho) - 2),
        sauter_params.f_trapped_particles,
        label="$f_{trapped}$",
    )

    axes_coefficients = fig.add_subplot(223)
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

    axes_collisionality = fig.add_subplot(324)
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
    axes_collisionality.set_xlabel("Normalized Radius")
    axes_collisionality.set_ylabel("Collisionality")
    axes_collisionality.grid(True)
    axes_collisionality.legend()

    axes_coefficients.legend()
    axes_coefficients.set_xlabel("Normalized Radius")
    axes_coefficients.set_ylabel("Coefficient Value")
    axes_coefficients.grid(True)
    j_profile.set_xlabel("Normalized Radius")
    j_profile.set_ylabel("Bootstrap Current")
    j_profile.grid(True)
    j_profile.legend()
    axes_trapped.set_xlabel("Normalized Radius")
    axes_trapped.set_ylabel("Trapped Fraction")
    axes_trapped.grid(True)
    axes_trapped.legend()
    q_profile = fig.add_subplot(326)
    q_profile.plot(
        rho,
        sauter_params.q_profile,
        label="Assumed $q$ profile",
    )
    q_profile.set_xlabel("Normalized Radius")
    q_profile.set_ylabel("Safety Factor")
    q_profile.grid(True)
    q_profile.legend()
