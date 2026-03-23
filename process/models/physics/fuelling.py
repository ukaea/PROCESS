import logging

import matplotlib.pyplot as plt
import numpy as np

import process.core.io.mfile as mf
from process.core import constants

logger = logging.getLogger(__name__)


class PlasmaFuelling:
    """Class to hold plasma fuelling calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    @staticmethod
    def calculate_plasma_tritium_flow_rate(
        f_molflow_plasma_fuelling_tritium,
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
        return (
            (
                f_molflow_plasma_fuelling_tritium
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            - fusrat_dt_total
            + fusrat_plasma_dd_triton
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_tritium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_deuterium_flow_rate(
        f_molflow_plasma_fuelling_deuterium,
        eta_plasma_fuelling,
        molflow_plasma_fuelling_vv_injected,
        fusrat_dt_total,
        fusrat_plasma_dhe3,
        fusrat_plasma_dd_total,
        t_energy_confinement,
        f_plasma_particles_lcfs_recycled,
        nd_plasma_fuel_ions_vol_avg,
        vol_plasma,
        f_plasma_fuel_deuterium,
    ):
        """Calculate the deuterium flow rate in the plasma exhaust."""
        return (
            (
                f_molflow_plasma_fuelling_deuterium
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            - fusrat_dt_total
            - 2 * fusrat_plasma_dd_total
            - fusrat_plasma_dhe3
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_deuterium)
                / (t_energy_confinement / (1 - f_plasma_particles_lcfs_recycled))
            )
        )

    @staticmethod
    def calculate_plasma_alphas_flow_rate(
        fusrat_dt_total,
        fusrat_plasma_dhe3,
        t_energy_confinement,
        f_t_alpha_energy_confinement,
        nd_plasma_alphas_vol_avg,
        vol_plasma,
    ):
        """Calculate the alpha particle flow rate in the plasma exhaust."""

        # Alpha particle balance

        return (
            fusrat_dt_total
            + fusrat_plasma_dhe3
            - (nd_plasma_alphas_vol_avg * vol_plasma)
            / (t_energy_confinement * f_t_alpha_energy_confinement)
        )

    def plot_tritium_flow_contour(self, axis: plt.Axes, mfile: mf.MFile, scan: int):
        """Plot contour of tritium flow rate vs recycling and fuelling rate."""

        recycling_range = np.linspace(0.01, 0.99, 20)
        fuelling_range = np.linspace(0.01, 1.0, 20)
        tritium_flow = np.zeros((len(recycling_range), len(fuelling_range)))

        for i, recycling in enumerate(recycling_range):
            for j, fuelling in enumerate(fuelling_range):
                tritium_flow[i, j] = self.calculate_plasma_tritium_flow_rate(
                    f_molflow_plasma_fuelling_tritium=mfile.get(
                        "f_molflow_plasma_fuelling_tritium", scan=scan
                    ),
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
                    f_molflow_plasma_fuelling_deuterium=mfile.get(
                        "f_molflow_plasma_fuelling_deuterium", scan=scan
                    ),
                    eta_plasma_fuelling=fuelling,
                    molflow_plasma_fuelling_vv_injected=mfile.get(
                        "molflow_plasma_fuelling_vv_injected", scan=scan
                    ),
                    fusrat_dt_total=mfile.get("fusrat_dt_total", scan=scan),
                    fusrat_plasma_dhe3=mfile.get("fusrat_plasma_dhe3", scan=scan),
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

        fusion_dt_range = np.linspace(1e19, 5e21, 20)
        f_t_alpha_energy_confinement_range = np.linspace(2.0, 10.0, 20)
        alpha_flow = np.zeros((
            len(fusion_dt_range),
            len(f_t_alpha_energy_confinement_range),
        ))

        for i, fusion_dt in enumerate(fusion_dt_range):
            for j, f_t_alpha_energy_confinement in enumerate(
                f_t_alpha_energy_confinement_range
            ):
                alpha_flow[i, j] = self.calculate_plasma_alphas_flow_rate(
                    fusrat_dt_total=fusion_dt,
                    fusrat_plasma_dhe3=mfile.get("fusrat_plasma_dhe3", scan=scan),
                    t_energy_confinement=mfile.get("t_energy_confinement", scan=scan),
                    nd_plasma_alphas_vol_avg=mfile.get(
                        "nd_plasma_alphas_vol_avg", scan=scan
                    ),
                    vol_plasma=mfile.get("vol_plasma", scan=scan),
                    f_t_alpha_energy_confinement=f_t_alpha_energy_confinement,
                )

        contour = axis.contourf(
            f_t_alpha_energy_confinement_range,
            fusion_dt_range,
            alpha_flow,
            levels=15,
            cmap="RdBu_r",
        )
        axis.contour(
            f_t_alpha_energy_confinement_range,
            fusion_dt_range,
            alpha_flow,
            levels=[0],
            colors="black",
            linewidths=2,
        )

        # Plot star for mfile values
        fusion_dt_mfile = mfile.get("fusrat_dt_total", scan=scan)
        f_t_alpha_mfile = mfile.get("f_alpha_energy_confinement", scan=scan)
        axis.plot(
            f_t_alpha_mfile,
            fusion_dt_mfile,
            marker="*",
            markersize=15,
            color="yellow",
            markeredgecolor="black",
            markeredgewidth=1.5,
        )

        axis.set_xlabel(
            "Alpha to Energy Confinement Time Ratio ($f_{\\alpha, \\text{energy confinement}}$)"
        )
        axis.set_ylabel("Fusion DT Rate [$\\text{particles/s}$]")
        axis.set_title("Plasma Alpha Particle Flow Rate (particles/s)")
        axis.minorticks_on()
        axis.grid(True, which="major", linestyle="-", alpha=0.7)
        axis.grid(True, which="minor", linestyle=":", alpha=0.4)
        plt.colorbar(contour, ax=axis, label="Alpha Particle Flow Rate")
