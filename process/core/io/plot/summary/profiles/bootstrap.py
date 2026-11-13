"""Plotting utilities for bootstrap current summary plots using MFile data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from process.core.io.plot.summary.profiles.misc import profiles_with_pedestal
from process.models.physics.bootstrap_current import SauterBootstrapCurrent

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_sauter_bootstrap_summary(axis: plt.Axes, mfile: MFile, scan: int):
    """Plot the Sauter bootstrap current summary for a given MFile and scan."""
    rho, nd_electron_profile, temp_electron_profile_kev = profiles_with_pedestal(
        mfile=mfile, scan=scan
    )
    sauter = SauterBootstrapCurrent()

    _, j_boot_profile = sauter.bootstrap_fraction_sauter(
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

    axis.plot(
        np.linspace(0, 1, len(rho) - 2), j_boot_profile, label="Sauter Bootstrap Current"
    )
