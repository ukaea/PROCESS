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
    def calculate_fuel_burnup_fraction(
        fusrat_total: float, molflow_plasma_fuelling_vv_injected: float
    ) -> float:
        """Calculate the fuel burnup fraction

        Parameters
        ----------
        fusrat_total : float
            Total fusion rate (particles/s).
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate into vacuum vessel (particles/s).

        Returns
        -------
        float            Fuel burnup fraction (dimensionless).

        """

        return fusrat_total / molflow_plasma_fuelling_vv_injected

    @staticmethod
    def calculate_plasma_tritium_flow_rate(
        f_molflow_plasma_fuelling_tritium: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_dt_total: float,
        fusrat_plasma_dd_triton: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_tritium: float,
    ) -> float:
        """Calculate the tritium flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_tritium : float
            Fraction of tritium in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        fusrat_plasma_dd_triton : float
            Tritium production rate from DD fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_tritium : float
            Fraction of tritium in the plasma fuel.

        Returns
        -------
        float
            Tritium flow rate in the plasma exhaust (particles/s).

        """

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
        f_molflow_plasma_fuelling_deuterium: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_dt_total: float,
        fusrat_plasma_dhe3: float,
        fusrat_plasma_dd_total: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_deuterium: float,
    ) -> float:
        """Calculate the deuterium flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_deuterium : float
            Fraction of deuterium in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_dt_total : float
            Total DT fusion rate (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).
        fusrat_plasma_dd_total : float
            Total deuterium consumption rate from DD fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_deuterium : float
            Fraction of deuterium in the plasma fuel.

        Returns
        -------
        float
            Deuterium flow rate in the plasma exhaust (particles/s).


        """
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
    def calculate_plasma_helium3_flow_rate(
        f_molflow_plasma_fuelling_helium3: float,
        eta_plasma_fuelling: float,
        molflow_plasma_fuelling_vv_injected: float,
        fusrat_plasma_dhe3: float,
        t_energy_confinement: float,
        f_plasma_particles_lcfs_recycled: float,
        nd_plasma_fuel_ions_vol_avg: float,
        vol_plasma: float,
        f_plasma_fuel_helium3: float,
    ) -> float:
        """Calculate the helium-3 flow rate in the plasma exhaust.

        Parameters
        ----------
        f_molflow_plasma_fuelling_helium3 : float
            Fraction of helium-3 in the plasma fuelling.
        eta_plasma_fuelling : float
            Fuelling rate efficiency.
        molflow_plasma_fuelling_vv_injected : float
            Total fuelling rate (particles/s).
        fusrat_plasma_dhe3 : float
            Deuterium consumption rate from D-He3 fusion (particles/s).
        t_energy_confinement : float
            Energy confinement time (s).
        f_plasma_particles_lcfs_recycled : float
            Fraction of plasma particles recycled at the LCFS.
        nd_plasma_fuel_ions_vol_avg : float
            Volume-averaged density of fuel ions in the plasma (particles/m^3).
        vol_plasma : float
            Plasma volume (m^3).
        f_plasma_fuel_helium3 : float
            Fraction of helium-3 in the plasma fuel.

        Returns
        -------
        float
            Helium-3 flow rate in the plasma exhaust (particles/s).

        """

        return (
            (
                f_molflow_plasma_fuelling_helium3
                * eta_plasma_fuelling
                * molflow_plasma_fuelling_vv_injected
            )
            + fusrat_plasma_dhe3
            - (
                (nd_plasma_fuel_ions_vol_avg * vol_plasma * f_plasma_fuel_helium3)
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

    def plot_fuelling_info(self, fig: plt.Figure, mfile: mf.MFile, scan: int):
        """Plot fuelling information."""
        msg = (
            f"$\\mathbf{{Plasma \\ Fuelling \\ Information:}}$\n\n"
            f"Total fuelling rate:"
            f"{mfile.get('molflow_plasma_fuelling_vv_injected', scan=scan):.4e} particles/s\n"
            f"Fuelling Rate Efficiency ($\\eta_{{\\text{{fuelling}}}}$): "
            f"{mfile.get('eta_plasma_fuelling', scan=scan):.4f}\n"
            f"Recycling Fraction ($R$): "
            f"{mfile.get('f_plasma_particles_lcfs_recycled', scan=scan):.4f}\n\n"
            f"Fraction of Tritium Fuelling: "
            f"{mfile.get('f_molflow_plasma_fuelling_tritium', scan=scan):.4f}\n"
            f"Fraction of Deuterium Fuelling: "
            f"{mfile.get('f_molflow_plasma_fuelling_deuterium', scan=scan):.4f}\n"
            f"Fraction of 3-Helium Fuelling: "
            f"{mfile.get('f_molflow_plasma_fuelling_helium3', scan=scan):.4f}"
        )
        fig.text(
            0.75,
            0.25,
            msg,
            ha="center",
            va="center",
            transform=fig.transFigure,
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 1.0},
        )
