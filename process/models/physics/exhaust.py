import logging

import matplotlib.pyplot as plt
import numpy as np

import process.core.io.mfile as mf
from process.core import constants

logger = logging.getLogger(__name__)


class PlasmaExhaust:
    """Class to hold plasma exhaust calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    @staticmethod
    def calculate_separatrix_power(
        f_p_alpha_plasma_deposited: float,
        p_alpha_total_mw: float,
        p_non_alpha_charged_mw: float,
        p_hcd_injected_total_mw: float,
        p_plasma_ohmic_mw: float,
        p_plasma_rad_mw: float,
    ) -> float:
        """
        Calculate the power crossing the separatrix (P_sep).

        Parameters
        ----------
        f_p_alpha_plasma_deposited : float
            Fraction of alpha power deposited in plasma.
        p_alpha_total_mw : float
            Total alpha power produced (MW).
        p_non_alpha_charged_mw : float
            Power from non-alpha charged particles (MW).
        p_hcd_injected_total_mw : float
            Total power injected by heating and current drive (MW).
        p_plasma_ohmic_mw : float
            Ohmic heating power (MW).
        p_plasma_rad_mw : float
            Radiated power from plasma (MW).

        Returns
        -------
        float
            Power crossing the separatrix (MW).
        """

        return (
            f_p_alpha_plasma_deposited * p_alpha_total_mw
            + p_non_alpha_charged_mw
            + p_hcd_injected_total_mw
            + p_plasma_ohmic_mw
            - p_plasma_rad_mw
        )

    @staticmethod
    def calculate_psep_over_r_metric(
        p_plasma_separatrix_mw: float, rmajor: float
    ) -> float:
        """
        Calculate the power crossing the separatrix per unit major radius (P_sep/R).

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (MW).
        rmajor : float
            Plasma major radius (m).

        Returns
        -------
        float
            Power crossing the separatrix per unit major radius (MW/m).
        """
        return p_plasma_separatrix_mw / rmajor

    @staticmethod
    def calculate_eu_demo_re_attachment_metric(
        p_plasma_separatrix_mw: float,
        b_plasma_toroidal_on_axis: float,
        q95: float,
        aspect: float,
        rmajor: float,
    ) -> float:
        """Calculate the EU DEMO divertor protection re-attachment metric for plasma exhaust.

        Parameters
        ----------
        p_plasma_separatrix_mw : float
            Power crossing the separatrix (MW).
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field on axis (T).
        q95 : float
            Safety factor at 95% flux surface.
        aspect : float
            Aspect ratio of the plasma.
        rmajor : float
            Plasma major radius (m).

        Returns
        -------
        float
            EU DEMO re-attachment metric (MW T /m).

        References
        ----------
        - M. Siccinio, G. Federici, R. Kembleton, H. Lux, F. Maviglia, and J. Morris,
          "Figure of merit for divertor protection in the preliminary design of the EU-DEMO reactor,"
          Nuclear Fusion, vol. 59, no. 10, pp. 106026-106026, Jul. 2019,
          doi: https://doi.org/10.1088/1741-4326/ab3153.

        - H. Zohm et al.,
          "A stepladder approach to a tokamak fusion power plant,"
          Nuclear Fusion, vol. 57, no. 8, pp. 086002-086002, May 2017,
          doi: https://doi.org/10.1088/1741-4326/aa739e.
        """

        return (p_plasma_separatrix_mw * b_plasma_toroidal_on_axis) / (
            q95 * aspect * rmajor
        )

    @staticmethod
    def calculate_plasma_tritium_flow_rate(
        eta_plasma_fuelling,
        molflow_plasma_fuelling_vv_injected,
        fusrat_dt_total,
        fusrat_plasma_dd_triton,
        t_energy_confinement,
        f_plasma_particles_lcfs_recycled,
        nd_plasma_fuel_ions_vol_avg,
        vol_plasma,
        f_plasma_fuel_tritium,
    ):
        """Calculate the tritium flow rate in the plasma exhaust."""
        # Assuming 50/50 D-T fuel mix, the tritium fuelling rate is half the total fuelling rate, and the tritium loss includes both D and T contributions from fusion reactions.
        return (
            (eta_plasma_fuelling * molflow_plasma_fuelling_vv_injected / 2)
            - fusrat_dt_total
            + fusrat_plasma_dd_triton
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_tritium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_deuterium_flow_rate(
        eta_plasma_fuelling,
        molflow_plasma_fuelling_vv_injected,
        fusrat_dt_total,
        fusrat_plasma_dd_total,
        t_energy_confinement,
        f_plasma_particles_lcfs_recycled,
        nd_plasma_fuel_ions_vol_avg,
        vol_plasma,
        f_plasma_fuel_deuterium,
    ):
        """Calculate the deuterium flow rate in the plasma exhaust."""
        # Assuming 50/50 D-T fuel mix, the deuterium fuelling rate is half the total fuelling rate, and the deuterium loss includes both D and T contributions from fusion reactions.
        return (
            (eta_plasma_fuelling * molflow_plasma_fuelling_vv_injected / 2)
            - fusrat_dt_total
            - 2 * fusrat_plasma_dd_total
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_deuterium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_alphas_flow_rate(
        fusrat_dt_total,
        t_energy_confinement,
        f_t_alpha_energy_confinement,
        f_plasma_particles_lcfs_recycled,
        nd_plasma_alphas_vol_avg,
        vol_plasma,
    ):
        """Calculate the alpha particle flow rate in the plasma exhaust."""

        # Alpha particle balance

        return fusrat_dt_total - (nd_plasma_alphas_vol_avg * vol_plasma) / (
            (t_energy_confinement * f_t_alpha_energy_confinement)
            / (1 - f_plasma_particles_lcfs_recycled)
        )

    def plot_tritium_flow_contour(self, axis: plt.Axes, mfile: mf.MFile, scan: int):
        """Plot contour of tritium flow rate vs recycling and fuelling rate."""

        recycling_range = np.linspace(0.01, 0.99, 20)
        fuelling_range = np.linspace(0.01, 1.0, 20)
        tritium_flow = np.zeros((len(recycling_range), len(fuelling_range)))

        for i, recycling in enumerate(recycling_range):
            for j, fuelling in enumerate(fuelling_range):
                tritium_flow[i, j] = self.calculate_plasma_tritium_flow_rate(
                    eta_plasma_fuelling=fuelling,
                    molflow_plasma_fuelling_vv_injected=mfile.get(
                        "molflow_plasma_fuelling_vv_injected", scan=scan
                    ),
                    fusrat_dt_total=mfile.get("fusrat_dt_total", scan=scan),
                    fusrat_plasma_dd_triton=mfile.get(
                        "fusrat_plasma_dd_triton", scan=scan
                    ),
                    t_energy_confinement=mfile.get("t_energy_confinement", scan=scan),
                    f_plasma_particles_lcfs_recycled=recycling,
                    nd_plasma_fuel_ions_vol_avg=mfile.get(
                        "nd_plasma_fuel_ions_vol_avg", scan=scan
                    ),
                    vol_plasma=mfile.get("vol_plasma", scan=scan),
                    f_plasma_fuel_tritium=mfile.get("f_plasma_fuel_tritium", scan=scan),
                )

        contour = axis.contourf(
            fuelling_range, recycling_range, tritium_flow, levels=15, cmap="RdBu_r"
        )
        axis.contour(
            fuelling_range,
            recycling_range,
            tritium_flow,
            levels=[0],
            colors="black",
            linewidths=2,
        )

        # Plot star for mfile values
        recycling_mfile = mfile.get("f_plasma_particles_lcfs_recycled", scan=scan)
        fuelling_mfile = mfile.get("eta_plasma_fuelling", scan=scan)
        axis.plot(
            fuelling_mfile,
            recycling_mfile,
            marker="*",
            markersize=15,
            color="yellow",
            markeredgecolor="black",
            markeredgewidth=1.5,
        )

        axis.set_xlabel("Fuelling Rate Efficiency ($\\eta_{\\text{fuelling}}$)")
        axis.set_ylabel("Recycling Fraction [$R$]")
        axis.set_title("Plasma Tritium Flow Rate (particles/s)")
        axis.minorticks_on()
        axis.grid(True, which="major", linestyle="-", alpha=0.7)
        axis.grid(True, which="minor", linestyle=":", alpha=0.4)
        plt.colorbar(contour, ax=axis, label="Tritium Flow Rate")

    def plot_deuterium_flow_contour(self, axis: plt.Axes, mfile: mf.MFile, scan: int):
        """Plot contour of deuterium flow rate vs recycling and fuelling rate."""

        recycling_range = np.linspace(0.01, 0.99, 20)
        fuelling_range = np.linspace(0.01, 1.0, 20)
        deuterium_flow = np.zeros((len(recycling_range), len(fuelling_range)))

        for i, recycling in enumerate(recycling_range):
            for j, fuelling in enumerate(fuelling_range):
                deuterium_flow[i, j] = self.calculate_plasma_deuterium_flow_rate(
                    eta_plasma_fuelling=fuelling,
                    molflow_plasma_fuelling_vv_injected=mfile.get(
                        "molflow_plasma_fuelling_vv_injected", scan=scan
                    ),
                    fusrat_dt_total=mfile.get("fusrat_dt_total", scan=scan),
                    fusrat_plasma_dd_total=mfile.get(
                        "fusrat_plasma_dd_total", scan=scan
                    ),
                    t_energy_confinement=mfile.get("t_energy_confinement", scan=scan),
                    f_plasma_particles_lcfs_recycled=recycling,
                    nd_plasma_fuel_ions_vol_avg=mfile.get(
                        "nd_plasma_fuel_ions_vol_avg", scan=scan
                    ),
                    vol_plasma=mfile.get("vol_plasma", scan=scan),
                    f_plasma_fuel_deuterium=mfile.get(
                        "f_plasma_fuel_deuterium", scan=scan
                    ),
                )

        contour = axis.contourf(
            fuelling_range, recycling_range, deuterium_flow, levels=15, cmap="RdBu_r"
        )
        axis.contour(
            fuelling_range,
            recycling_range,
            deuterium_flow,
            levels=[0],
            colors="black",
            linewidths=2,
        )

        # Plot star for mfile values
        recycling_mfile = mfile.get("f_plasma_particles_lcfs_recycled", scan=scan)
        fuelling_mfile = mfile.get("eta_plasma_fuelling", scan=scan)
        axis.plot(
            fuelling_mfile,
            recycling_mfile,
            marker="*",
            markersize=15,
            color="yellow",
            markeredgecolor="black",
            markeredgewidth=1.5,
        )

        axis.set_xlabel("Fuelling Rate Efficiency ($\\eta_{\\text{fuelling}}$)")
        axis.set_ylabel("Recycling Fraction [$R$]")
        axis.set_title("Plasma Deuterium Flow Rate (particles/s)")
        axis.minorticks_on()
        axis.grid(True, which="major", linestyle="-", alpha=0.7)
        axis.grid(True, which="minor", linestyle=":", alpha=0.4)
        plt.colorbar(contour, ax=axis, label="Deuterium Flow Rate")

    def plot_alpha_flow_contour(self, axis: plt.Axes, mfile: mf.MFile, scan: int):
        """Plot contour of alpha particle flow rate vs recycling and fuelling rate."""

        recycling_range = np.linspace(0.01, 0.99, 20)
        f_t_alpha_energy_confinement_range = np.linspace(4.0, 10.0, 20)
        alpha_flow = np.zeros((
            len(recycling_range),
            len(f_t_alpha_energy_confinement_range),
        ))

        for i, recycling in enumerate(recycling_range):
            for j, f_t_alpha_energy_confinement in enumerate(
                f_t_alpha_energy_confinement_range
            ):
                alpha_flow[i, j] = self.calculate_plasma_alphas_flow_rate(
                    fusrat_dt_total=mfile.get("fusrat_dt_total", scan=scan),
                    t_energy_confinement=mfile.get("t_energy_confinement", scan=scan),
                    f_plasma_particles_lcfs_recycled=recycling,
                    nd_plasma_alphas_vol_avg=mfile.get(
                        "nd_plasma_alphas_vol_avg", scan=scan
                    ),
                    vol_plasma=mfile.get("vol_plasma", scan=scan),
                    f_t_alpha_energy_confinement=f_t_alpha_energy_confinement,
                )

        contour = axis.contourf(
            f_t_alpha_energy_confinement_range,
            recycling_range,
            alpha_flow,
            levels=15,
            cmap="RdBu_r",
        )
        axis.contour(
            f_t_alpha_energy_confinement_range,
            recycling_range,
            alpha_flow,
            levels=[0],
            colors="black",
            linewidths=2,
        )

        # Plot star for mfile values
        recycling_mfile = mfile.get("f_plasma_particles_lcfs_recycled", scan=scan)
        f_t_alpha_mfile = mfile.get("f_alpha_energy_confinement", scan=scan)
        axis.plot(
            f_t_alpha_mfile,
            recycling_mfile,
            marker="*",
            markersize=15,
            color="yellow",
            markeredgecolor="black",
            markeredgewidth=1.5,
        )

        axis.set_xlabel(
            "Alpha to Energy Confinement Time Ratio ($f_{\\alpha, \\text{energy confinement}}$)"
        )
        axis.set_ylabel("Recycling Fraction [$R$]")
        axis.set_title("Plasma Alpha Particle Flow Rate (particles/s)")
        axis.minorticks_on()
        axis.grid(True, which="major", linestyle="-", alpha=0.7)
        axis.grid(True, which="minor", linestyle=":", alpha=0.4)
        plt.colorbar(contour, ax=axis, label="Alpha Particle Flow Rate")
