from importlib import resources

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

import process.core.io.mfile as mf
import process.data_structure.tritium_plant_variables as tritium
from process.core import constants
from process.core import process_output as po


class TritiumPlantMeschini:
    """Class to model the tritium plant inventory and flow rates.


        :notes: The original paper labels the systems as follows:
                # 1 = Blanket
                # 2 = Tritium Extraction System
                # 3 = First Wall
                # 4 = Divertor
                # 5 = Heat Exchanger
                # 6 = Detritiation System
                # 7 = Vacuum Pump
                # 8 = Fuel cleanup
                # 9 = Isotope Separation System
                # 10 Storage and management
                # 11 = Fuelling system
                # 12 = Tritium sepration membrane

        :reference:

            - Samuele Meschini, S. Ferry, Rémi Delaporte-Mathurin, and D. G. Whyte,
            “Modeling and analysis of the tritium fuel cycle for ARC- and STEP-class D-T fusion power plants,”
            Nuclear Fusion, vol. 63, no. 12, pp. 126005-126005, Sep. 2023,
            doi: https://doi.org/10.1088/1741-4326/acf3fc.
    ‌

    """

    def __init__(self) -> None:
        self.outfile: int = constants.NOUT

    def run(self) -> None:
        """Run the tritium plant model."""

        t_end = 20000  # End time for simulation (s)
        dt = 1  # Time step (s)

        RATE_T_DECAY = 1.78e-9  # Tritium decay rate (1/s)

        TBR = 0.5
        # Tritium burn efficiency in the plasma
        TBE = 0.05

        # Direct internal recycling fraction
        f_dir = 0.1

        # Tritium residence time in the ith component
        tau_1 = 5.0
        tau_2 = 10.0
        tau_3 = 15.0
        tau_4 = 20.0
        tau_5 = 25.0
        tau_6 = 30.0
        tau_7 = 35.0
        tau_8 = 40.0
        tau_9 = 45.0
        tau_10 = 50.0
        tau_11 = 55.0
        tau_12 = 60.0

        n_t_burn = 3e-6  # Tritium burn rate in the plasma (kg/s)

        i_startup = 5.0

        ETA_2 = 0.5

        # Flow rate fractions between components
        f_5_1 = 0.33
        f_5_3 = 0.33
        f_5_6 = 1e-4
        f_9_6 = 0.1
        f_p_3 = 1e-4
        f_p_4 = 1e-4
        f_5_4 = 0.33

        # The non-radioactive loss fraction has been assumed asb in this work.
        epsilon_1 = 1e-4
        epsilon_2 = 1e-4
        epsilon_3 = 0.0
        epsilon_4 = 0.0
        epsilon_5 = 1e-4
        epsilon_6 = 1e-4
        epsilon_7 = 1e-4
        epsilon_8 = 1e-4
        epsilon_9 = 1e-4
        epsilon_10 = 0.0
        epsilon_11 = 1e-4
        epsilon_12 = 1e-4

        # Tritium inventory in the ith component (kg)
        # Initialise plant to be zero and clean with startup inventory in storage and management
        m_tritium_component_1 = 0.0
        m_tritium_component_2 = 0.0
        m_tritium_component_3 = 0.0
        m_tritium_component_4 = 0.0
        m_tritium_component_5 = 0.0
        m_tritium_component_6 = 0.0
        m_tritium_component_7 = 0.0
        m_tritium_component_8 = 0.0
        m_tritium_component_9 = 0.0
        m_tritium_component_10 = i_startup
        m_tritium_component_11 = 0.0
        m_tritium_component_12 = 0.0

        n_steps = int(t_end / dt) + 1
        times = [i * dt for i in range(n_steps)]

        # Initialize inventory arrays
        inventories = {
            "i_1": [m_tritium_component_1],
            "i_2": [m_tritium_component_2],
            "i_3": [m_tritium_component_3],
            "i_4": [m_tritium_component_4],
            "i_5": [m_tritium_component_5],
            "i_6": [m_tritium_component_6],
            "i_7": [m_tritium_component_7],
            "i_8": [m_tritium_component_8],
            "i_9": [m_tritium_component_9],
            "i_10": [m_tritium_component_10],
            "i_11": [m_tritium_component_11],
            "i_12": [m_tritium_component_12],
        }

        # Current values
        curr = [
            m_tritium_component_1,
            m_tritium_component_2,
            m_tritium_component_3,
            m_tritium_component_4,
            m_tritium_component_5,
            m_tritium_component_6,
            m_tritium_component_7,
            m_tritium_component_8,
            m_tritium_component_9,
            m_tritium_component_10,
            m_tritium_component_11,
            m_tritium_component_12,
        ]

        for _ in range(1, n_steps):
            # Calculate rates
            rate_1 = self.calculate_component_1_time_derivative(
                m_tritium_component_1=m_tritium_component_1,
                m_tritium_component_3=m_tritium_component_3,
                m_tritium_component_4=m_tritium_component_4,
                m_tritium_component_5=m_tritium_component_5,
                t_tritium_residence_component_1=tau_1,
                t_tritium_residence_component_3=tau_3,
                t_tritium_residence_component_4=tau_4,
                t_tritium_residence_component_5=tau_5,
                tbr=TBR,
                fusrat_tritium_kg=n_t_burn,
                f_5_1=f_5_1,
                epsilon_1=epsilon_1,
            )
            rate_2 = self.calculate_component_2_time_derivative(
                m_tritium_component_1=m_tritium_component_1,
                m_tritium_component_2=m_tritium_component_2,
                t_tritium_residence_component_1=tau_1,
                t_tritium_residence_component_2=tau_2,
                epsilon_2=epsilon_2,
            )
            rate_3 = self.calculate_component_3_time_derivative(
                m_tritium_component_3=m_tritium_component_3,
                m_tritium_component_5=m_tritium_component_5,
                t_tritium_residence_component_3=tau_3,
                t_tritium_residence_component_5=tau_5,
                f_p_3=f_p_3,
                f_5_3=f_5_3,
                fusrat_tritium_kg=n_t_burn,
                epsilon_3=epsilon_3,
            )
            rate_4 = self.calculate_component_4_time_derivative(
                m_tritium_component_4=m_tritium_component_4,
                m_tritium_component_5=m_tritium_component_5,
                t_tritium_residence_component_4=tau_4,
                t_tritium_residence_component_5=tau_5,
                f_p_4=f_p_4,
                f_5_4=f_5_4,
                fusrat_tritium_kg=n_t_burn,
                epsilon_4=epsilon_4,
            )
            rate_5 = self.calculate_component_5_time_derivative(
                m_tritium_component_2=m_tritium_component_2,
                m_tritium_component_5=m_tritium_component_5,
                t_tritium_residence_component_2=tau_2,
                t_tritium_residence_component_5=tau_5,
                epsilon_5=epsilon_5,
            )
            rate_6 = self.calculate_component_6_time_derivative(
                m_tritium_component_5=m_tritium_component_5,
                m_tritium_component_6=m_tritium_component_6,
                m_tritium_component_9=m_tritium_component_9,
                t_tritium_residence_component_5=tau_5,
                t_tritium_residence_component_6=tau_6,
                t_tritium_residence_component_9=tau_9,
                f_5_6=f_5_6,
                f_9_6=f_9_6,
                epsilon_6=epsilon_6,
            )
            rate_7 = self.calculate_component_7_time_derivative(
                m_tritium_component_7=m_tritium_component_7,
                t_tritium_residence_component_7=tau_7,
                f_p_3=f_p_3,
                f_p_4=f_p_4,
                fusrat_tritium_kg=n_t_burn,
                epsilon_7=epsilon_7,
            )
            rate_8 = self.calculate_component_8_time_derivative(
                m_tritium_component_7=m_tritium_component_7,
                m_tritium_component_8=m_tritium_component_8,
                t_tritium_residence_component_7=tau_7,
                t_tritium_residence_component_8=tau_8,
                f_dir=f_dir,
                epsilon_8=epsilon_8,
            )
            rate_9 = self.calculate_component_9_time_derivative(
                m_tritium_component_6=m_tritium_component_6,
                m_tritium_component_8=m_tritium_component_8,
                m_tritium_component_9=m_tritium_component_9,
                t_tritium_residence_component_6=tau_6,
                t_tritium_residence_component_8=tau_8,
                t_tritium_residence_component_9=tau_9,
                epsilon_9=epsilon_9,
            )
            rate_10 = self.calculate_component_10_time_derivative(
                m_tritium_component_7=m_tritium_component_7,
                m_tritium_component_9=m_tritium_component_9,
                m_tritium_component_10=m_tritium_component_10,
                m_tritium_component_12=m_tritium_component_12,
                t_tritium_residence_component_7=tau_7,
                t_tritium_residence_component_9=tau_9,
                t_tritium_residence_component_12=tau_12,
                f_9_6=f_9_6,
                f_dir=f_dir,
                n_t_burn=n_t_burn,
            )
            rate_12 = self.calculate_component_12_time_derivative(
                m_tritium_component_2=m_tritium_component_2,
                m_tritium_component_12=m_tritium_component_12,
                t_tritium_residence_component_2=tau_2,
                t_tritium_residence_component_12=tau_12,
                epsilon_12=epsilon_12,
            )
            # Update inventories using Euler method
            curr[0] += rate_1 * dt
            curr[1] += rate_2 * dt
            curr[2] += rate_3 * dt
            curr[3] += rate_4 * dt
            curr[4] += rate_5 * dt
            curr[5] += rate_6 * dt
            curr[6] += rate_7 * dt
            curr[7] += rate_8 * dt
            curr[8] += rate_9 * dt
            curr[9] += rate_10 * dt
            curr[11] += rate_12 * dt
            # i_11 remains 0

            # Store values
            for idx, key in enumerate([
                "i_1",
                "i_2",
                "i_3",
                "i_4",
                "i_5",
                "i_6",
                "i_7",
                "i_8",
                "i_9",
                "i_10",
                "i_11",
                "i_12",
            ]):
                inventories[key].append(curr[idx])

        return times, inventories

    def calculate_component_1_time_derivative(
        m_tritium_component_1: float,
        m_tritium_component_3: float,
        m_tritium_component_4: float,
        m_tritium_component_5: float,
        t_tritium_residence_component_1: float,
        t_tritium_residence_component_3: float,
        t_tritium_residence_component_4: float,
        t_tritium_residence_component_5: float,
        tbr: float,
        fusrat_tritium_kg: float,
        f_5_1: float,
        epsilon_1: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in blanket (component 1)."""
        return (
            tbr * fusrat_tritium_kg
            + (m_tritium_component_3 / t_tritium_residence_component_3)
            + (m_tritium_component_4 / t_tritium_residence_component_4)
            + f_5_1 * (m_tritium_component_5 / t_tritium_residence_component_5)
            - (
                m_tritium_component_1
                * (((1 + epsilon_1) / t_tritium_residence_component_1) + RATE_T_DECAY)
            )
        )

    def calculate_component_2_time_derivative(
        m_tritium_component_1: float,
        m_tritium_component_2: float,
        t_tritium_residence_component_1: float,
        t_tritium_residence_component_2: float,
        epsilon_2: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in tritium extraction system (component 2)."""
        return (m_tritium_component_1 / t_tritium_residence_component_1) - (
            m_tritium_component_2
            * (((1 + epsilon_2) / t_tritium_residence_component_2) + RATE_T_DECAY)
        )

    def calculate_component_3_time_derivative(
        m_tritium_component_3: float,
        m_tritium_component_5: float,
        t_tritium_residence_component_3: float,
        t_tritium_residence_component_5: float,
        f_p_3: float,
        f_5_3: float,
        fusrat_tritium_kg: float,
        epsilon_3: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in first wall (component 3)."""
        return (
            f_p_3 * (fusrat_tritium_kg / TBE)
            + f_5_3 * (m_tritium_component_5 / t_tritium_residence_component_5)
            - (
                m_tritium_component_3
                * (((1 + epsilon_3) / t_tritium_residence_component_3) + RATE_T_DECAY)
            )
        )

    def calculate_component_4_time_derivative(
        m_tritium_component_4: float,
        m_tritium_component_5: float,
        t_tritium_residence_component_4: float,
        t_tritium_residence_component_5: float,
        f_p_4: float,
        f_5_4: float,
        fusrat_tritium_kg: float,
        epsilon_4: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in divertor (component 4)."""
        return (
            f_p_4 * (fusrat_tritium_kg / TBE)
            + f_5_4 * (m_tritium_component_5 / t_tritium_residence_component_5)
            - (
                m_tritium_component_4
                * (((1 + epsilon_4) / t_tritium_residence_component_4) + RATE_T_DECAY)
            )
        )

    def calculate_component_5_time_derivative(
        m_tritium_component_2: float,
        m_tritium_component_5: float,
        t_tritium_residence_component_2: float,
        t_tritium_residence_component_5: float,
        epsilon_5: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in heat exchanger (component 5)."""
        return (1 - ETA_2) * (
            m_tritium_component_2 / t_tritium_residence_component_2
        ) - (
            m_tritium_component_5
            * (((1 + epsilon_5) / t_tritium_residence_component_5) + RATE_T_DECAY)
        )

    def calculate_component_6_time_derivative(
        m_tritium_component_5: float,
        m_tritium_component_6: float,
        m_tritium_component_9: float,
        t_tritium_residence_component_5: float,
        t_tritium_residence_component_6: float,
        t_tritium_residence_component_9: float,
        f_5_6: float,
        f_9_6: float,
        epsilon_6: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in detritiation system (component 6)."""
        return (
            f_5_6 * (m_tritium_component_5 / t_tritium_residence_component_5)
            + f_9_6 * (m_tritium_component_9 / t_tritium_residence_component_9)
            - (
                m_tritium_component_6
                * (((1 + epsilon_6) / t_tritium_residence_component_6) + RATE_T_DECAY)
            )
        )

    def calculate_component_7_time_derivative(
        m_tritium_component_7: float,
        t_tritium_residence_component_7: float,
        f_p_3: float,
        f_p_4: float,
        fusrat_tritium_kg: float,
        epsilon_7: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in vacuum pump (component 7)."""
        return (1 - TBE - f_p_3 - f_p_4) * (fusrat_tritium_kg / TBE) - (
            m_tritium_component_7
            * (((1 + epsilon_7) / t_tritium_residence_component_7) + RATE_T_DECAY)
        )

    def calculate_component_8_time_derivative(
        m_tritium_component_7: float,
        m_tritium_component_8: float,
        t_tritium_residence_component_7: float,
        t_tritium_residence_component_8: float,
        f_dir: float,
        epsilon_8: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in fuel cleanup (component 8)."""
        return (1 - f_dir) * (
            m_tritium_component_7 / t_tritium_residence_component_7
        ) - (
            m_tritium_component_8
            * (((1 + epsilon_8) / t_tritium_residence_component_8) + RATE_T_DECAY)
        )

    def calculate_component_9_time_derivative(
        m_tritium_component_6: float,
        m_tritium_component_8: float,
        m_tritium_component_9: float,
        t_tritium_residence_component_6: float,
        t_tritium_residence_component_8: float,
        t_tritium_residence_component_9: float,
        epsilon_9: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in isotope separation system (component 9)."""
        return (
            (m_tritium_component_6 / t_tritium_residence_component_6)
            + (m_tritium_component_8 / t_tritium_residence_component_8)
            - (
                m_tritium_component_9
                * (((1 + epsilon_9) / t_tritium_residence_component_9) + RATE_T_DECAY)
            )
        )

    def calculate_component_10_time_derivative(
        m_tritium_component_7: float,
        m_tritium_component_9: float,
        m_tritium_component_10: float,
        m_tritium_component_12: float,
        t_tritium_residence_component_7: float,
        t_tritium_residence_component_9: float,
        t_tritium_residence_component_12: float,
        f_9_6: float,
        f_dir: float,
        n_t_burn: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in storage and management (component 10)."""
        return (
            (1 - f_9_6) * (m_tritium_component_9 / t_tritium_residence_component_9)
            + f_dir * (m_tritium_component_7 / t_tritium_residence_component_7)
            + (m_tritium_component_12 / t_tritium_residence_component_12)
            - (n_t_burn / TBE)
            - (RATE_T_DECAY * m_tritium_component_10)
        )

    def calculate_component_12_time_derivative(
        m_tritium_component_2: float,
        m_tritium_component_12: float,
        t_tritium_residence_component_2: float,
        t_tritium_residence_component_12: float,
        epsilon_12: float,
    ) -> float:
        """Calculate time derivative of tritium inventory in tritium separation membrane (component 12)."""
        return ETA_2 * (m_tritium_component_2 / t_tritium_residence_component_2) - (
            m_tritium_component_12
            * (((1 + epsilon_12) / t_tritium_residence_component_12) + RATE_T_DECAY)
        )

    def plot_tritium_systems_overview(
        axis: plt.Axes, m_file: mf.MFile, scan: int, fig: plt.Figure
    ):
        """Plot an overview of the tritium systems inventory and flow rates."""

        axis.text(
            0.1,
            0.9,
            "$\\epsilon$ = Non radioactive loss fraction \n$\\tau = $ Tritium residence time (s)",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Load the plasma image
        with resources.path("process.core.io", "plasma.png") as img_path:
            plasma = mpimg.imread(img_path.open("rb"))

        # Display the plasma image over the figure, not the axes
        new_ax = axis.inset_axes(
            [0.35, 0.35, 0.30, 0.30], transform=axis.transAxes, zorder=1
        )
        new_ax.imshow(plasma)
        new_ax.axis("off")
        axis.axis("off")  # Hide the main axes

        # Add an arrow from the plasma to the first wall box
        axis.annotate(
            "",
            xy=(0.65, 0.5),
            xytext=(0.55, 0.5),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Draw line from bottom of plasma to inline with divertor box height
        axis.annotate(
            "",
            xy=(0.5, 0.25),
            xytext=(0.5, 0.35),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 5,
                "fill": True,
            },
        )

        # Draw arrow from line above towards the divertor box
        axis.annotate(
            "",
            xy=(0.65, 0.25),
            xytext=(0.5, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add an arrow from the plasma to vacuum pump box
        axis.annotate(
            "",
            xy=(0.35, 0.5),
            xytext=(0.45, 0.5),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =====================================================

        # Add First Wall box
        axis.text(
            0.65,
            0.5,
            f"First Wall\n\n $\\epsilon$ = \n $\\tau = {tritium.t_fw_tritium_residence} s $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add an arrow from the first wall box to the blanket box
        axis.annotate(
            "",
            xy=(0.675, 0.665),
            xytext=(0.675, 0.55),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add blanket box
        axis.text(
            0.65,
            0.65,
            "Blanket\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # ============================================

        # Add divertor box
        axis.text(
            0.65,
            0.3,
            "Divertor\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add a small line from the right side of the divertor box (no arrow)
        axis.annotate(
            "",
            xy=(0.71, 0.25),
            xytext=(0.74, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add a vertical line stemming from the line above
        axis.annotate(
            "",
            xy=(0.73, 0.6),
            xytext=(0.73, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add a small line from the right side of the divertor box (no arrow)
        axis.annotate(
            "",
            xy=(0.71, 0.6),
            xytext=(0.74, 0.6),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.71, 0.675),
            xytext=(0.71, 0.6),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # ============================================

        # Add tritium extraction system box
        axis.text(
            0.8,
            0.5,
            "Tritium Extraction System\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add an arrow from tritium extraction system to the heat exchanger box
        axis.annotate(
            "",
            xy=(0.875, 0.15),
            xytext=(0.875, 0.47),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add line coming from right of tritum extraction system
        axis.annotate(
            "",
            xy=(0.95, 0.5),
            xytext=(1.0, 0.5),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(1.0, 0.9),
            xytext=(1.0, 0.5),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.75, 0.9),
            xytext=(1.0, 0.9),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # ==========================================

        # Add heat exchanger box
        axis.text(
            0.8,
            0.2,
            "Heat Exchanger\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "lightyellow",
                "alpha": 1.0,
                "linewidth": 2,
            },
        )

        # Add line from bottom of heat exchanger box to
        axis.annotate(
            "",
            xy=(0.9, 0.1),
            xytext=(0.9, 0.0),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add line from bottom of heat exchanger box to
        axis.annotate(
            "",
            xy=(0.2, 0.0),
            xytext=(0.9, 0.0),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.3, 0.825),
            xytext=(0.3, 0.9),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.2, 0.2),
            xytext=(0.2, 0.0),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =====================================================

        # Add fuelling system box
        axis.text(
            0.5,
            0.7,
            f"Fuelling System\n\n $\\epsilon$ = \n $\\tau = {m_file.get('t_plasma_fuelling_system_tritium_residence', scan=scan)} s $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        # Add an arrow from the plasma to the first wall box
        axis.annotate(
            "",
            xy=(0.5, 0.6),
            xytext=(0.5, 0.65),
            xycoords=fig.transFigure,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =====================================================

        # Add tritium seperation membrane box
        axis.text(
            0.6,
            0.8,
            "Tritium Separation Membrane\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        # Add line from tritium separation membrane to just above storage and management system box
        axis.annotate(
            "",
            xy=(0.5, 0.9),
            xytext=(0.3, 0.9),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add an arrow from the line above to the storage and management system box
        axis.annotate(
            "",
            xy=(0.3, 0.825),
            xytext=(0.3, 0.9),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # ============================================

        # Add vacuum pump box
        axis.text(
            0.35,
            0.5,
            "Vacuum Pump\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        # Add an arrow from vaccum pump to the storage
        axis.annotate(
            "",
            xy=(0.3, 0.7),
            xytext=(0.3, 0.55),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "linestyle": "--",
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add an arrow from vaccum pump to fuel cleanup system
        axis.annotate(
            "",
            xy=(0.1, 0.5),
            xytext=(0.25, 0.5),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # ==============================================

        # Add fuel cleanup system box
        axis.text(
            0.15,
            0.5,
            "Fuel Cleanup System\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        # Add an arrow from fuel cleanup system to isotope separation system
        axis.annotate(
            "",
            xy=(0.05, 0.7),
            xytext=(0.05, 0.55),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # Add an arrow from fuel cleanup system to storage and management system
        axis.annotate(
            "",
            xy=(0.2, 0.75),
            xytext=(0.125, 0.75),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =================================================

        # Add isotope speration system box
        axis.text(
            0.15,
            0.7,
            "Isotope Separation System\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        axis.annotate(
            "",
            xy=(-0.1, 0.75),
            xytext=(-0.065, 0.75),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(-0.1, 0.75),
            xytext=(-0.1, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.175, 0.25),
            xytext=(-0.1, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # ===================================================

        # Add detritiation system box
        axis.text(
            0.3,
            0.3,
            "Detritiation System\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        axis.annotate(
            "",
            xy=(-0.065, 0.8),
            xytext=(-0.15, 0.8),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(-0.15, 0.8),
            xytext=(-0.15, 0.25),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        axis.annotate(
            "",
            xy=(0.175, 0.225),
            xytext=(-0.15, 0.225),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =================================================

        # Add storage and management system box
        axis.text(
            0.35,
            0.7,
            "Storage and\n Management System\n\n $\\epsilon$ = \n $\\tau = $",
            fontsize=9,
            verticalalignment="center",
            horizontalalignment="center",
            transform=fig.transFigure,
            bbox={
                "boxstyle": "round",
                "facecolor": "cyan",
                "alpha": 0.8,
                "linewidth": 2,
            },
        )

        # Add an arrow from storage and management system to fuel system
        axis.annotate(
            "",
            xy=(0.45, 0.75),
            xytext=(0.375, 0.75),
            xycoords=axis.transAxes,
            arrowprops={
                "arrowstyle": "-|>,head_length=1,head_width=0.3",
                "color": "black",
                "linewidth": 1.5,
                "zorder": 11,
            },
            annotation_clip=False,
        )

        # =================================================

    def output_tritium_plant_info(self):
        """Output tritium plant information to the output file."""
        po.oheadr(self.outfile, "Tritium Plant Information")

        po.ovarin(
            self.outfile,
            "Divertor Tritium Residence Time (s)",
            "(t_div_tritium_residence)",
            tritium.t_div_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "First Wall Tritium Residence Time (s)",
            "(t_fw_tritium_residence)",
            tritium.t_fw_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Blanket Tritium Residence Time (s)",
            "(t_blkt_tritium_residence)",
            tritium.t_blkt_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Heat Exchanger Tritium Residence Time (s)",
            "(t_heat_exchanger_tritium_residence)",
            tritium.t_heat_exchanger_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Tritium Extraction System Residence Time (s)",
            "(t_tritium_extraction_system_tritium_residence)",
            tritium.t_tritium_extraction_system_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Tritium Separation Membrane Residence Time (s)",
            "(t_tritium_separation_membrane_tritium_residence)",
            tritium.t_tritium_separation_membrane_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Plasma Fuelling System Tritium Residence Time (s)",
            "(t_plasma_fuelling_system_tritium_residence)",
            tritium.t_plasma_fuelling_system_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Vacuum Pump Tritium Residence Time (s)",
            "(t_vacuum_pump_tritium_residence)",
            tritium.t_vacuum_pump_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Tritium Storage Residence Time (s)",
            "(t_tritium_storage_tritium_residence)",
            tritium.t_tritium_storage_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Isotope Separation Tritium Residence Time (s)",
            "(t_isotope_separation_tritium_residence)",
            tritium.t_isotope_separation_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Fuel Cleanup Tritium Residence Time (s)",
            "(t_fuel_cleanup_tritium_residence)",
            tritium.t_fuel_cleanup_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Detritiation System Tritium Residence Time (s)",
            "(t_detritiation_tritium_residence)",
            tritium.t_detritiation_tritium_residence,
        )
        po.ovarin(
            self.outfile,
            "Plant Tritium Start-up Inventory (kg)",
            "(m_plant_tritium_start_up)",
            tritium.m_plant_tritium_start_up,
        )
        po.ovarre(
            self.outfile,
            "Plant Tritium Start-up Minimum Required (kg)",
            "(m_plant_tritium_start_up_minimum_required)",
            tritium.m_plant_tritium_start_up_minimum_required,
            "OP",
        )


# RATE_T_DECAY = 1.78e-9  # Tritium decay rate (1/s)

# TBR = 0.5
# # Tritium burn efficiency in the plasma
# TBE = 0.05

# # Direct internal recycling fraction
# f_dir = 0.1

# # Tritium residence time in the ith component
# tau_1 = 5.0
# tau_2 = 10.0
# tau_3 = 15.0
# tau_4 = 20.0
# tau_5 = 25.0
# tau_6 = 30.0
# tau_7 = 35.0
# tau_8 = 40.0
# tau_9 = 45.0
# tau_10 = 50.0
# tau_11 = 55.0
# tau_12 = 60.0

# n_t_burn = 3e-6  # Tritium burn rate in the plasma (kg/s)

# i_startup = 5.0

# ETA_2 = 0.5

# # Flow rate fractions between components
# f_5_1 = 0.33
# f_5_3 = 0.33
# f_5_6 = 1e-4
# f_9_6 = 0.1
# f_p_3 = 1e-4
# f_p_4 = 1e-4
# f_5_4 = 0.33

# # The non-radioactive loss fraction has been assumed asb in this work.
# epsilon_1 = 1e-4
# epsilon_2 = 1e-4
# epsilon_3 = 0.0
# epsilon_4 = 0.0
# epsilon_5 = 1e-4
# epsilon_6 = 1e-4
# epsilon_7 = 1e-4
# epsilon_8 = 1e-4
# epsilon_9 = 1e-4
# epsilon_10 = 0.0
# epsilon_11 = 1e-4
# epsilon_12 = 1e-4

# # Tritium inventory in the ith component (kg)

# m_tritium_component_1 = 0.0
# m_tritium_component_2 = 0.0
# m_tritium_component_3 = 0.0
# m_tritium_component_4 = 0.0
# m_tritium_component_5 = 0.0
# m_tritium_component_6 = 0.0
# m_tritium_component_7 = 0.0
# m_tritium_component_8 = 0.0
# m_tritium_component_9 = 0.0
# m_tritium_component_10 = i_startup
# m_tritium_component_11 = 0.0
# m_tritium_component_12 = 0.0


# def calculate_rates(
#     m_tritium_component_1,
#     m_tritium_component_2,
#     m_tritium_component_3,
#     m_tritium_component_4,
#     m_tritium_component_5,
#     m_tritium_component_6,
#     m_tritium_component_7,
#     m_tritium_component_8,
#     m_tritium_component_9,
#     m_tritium_component_10,
#     m_tritium_component_12,
# ):
#     """Calculate rates of change for all components."""
#     d_m_tritium_component_1_dt = (
#         TBR * n_t_burn
#         + (m_tritium_component_3 / tau_3)
#         + (m_tritium_component_4 / tau_4)
#         + f_5_1 * (m_tritium_component_5 / tau_5)
#         - (m_tritium_component_1 * (((1 + epsilon_1) / tau_1) + RATE_T_DECAY))
#     )

#     d_m_tritium_component_2_dt = (m_tritium_component_1 / tau_1) - (
#         m_tritium_component_2 * (((1 + epsilon_2) / tau_2) + RATE_T_DECAY)
#     )

#     d_m_tritium_component_3_dt = (
#         f_p_3 * (n_t_burn / TBE)
#         + f_5_3 * (m_tritium_component_5 / tau_5)
#         - (m_tritium_component_3 * (((1 + epsilon_3) / tau_3) + RATE_T_DECAY))
#     )

#     d_m_tritium_component_4_dt = (
#         f_p_4 * (n_t_burn / TBE)
#         + f_5_4 * (m_tritium_component_5 / tau_5)
#         - (m_tritium_component_4 * (((1 + epsilon_4) / tau_4) + RATE_T_DECAY))
#     )

#     d_m_tritium_component_5_dt = (1 - ETA_2) * (m_tritium_component_2 / tau_2) - (
#         m_tritium_component_5 * (((1 + epsilon_5) / tau_5) + RATE_T_DECAY)
#     )

#     d_m_tritium_component_6_dt = (
#         f_5_6 * (m_tritium_component_5 / tau_5)
#         + f_9_6 * (m_tritium_component_9 / tau_9)
#         - (m_tritium_component_6 * (((1 + epsilon_6) / tau_6) + RATE_T_DECAY))
#     )

#     d_m_tritium_component_7_dt = (1 - TBE - f_p_3 - f_p_4) * (n_t_burn / TBE) - (
#         m_tritium_component_7 * (((1 + epsilon_7) / tau_7) + RATE_T_DECAY)
#     )

#     d_m_tritium_component_8_dt = (1 - f_dir) * (m_tritium_component_7 / tau_7) - (
#         m_tritium_component_8 * (((1 + epsilon_8) / tau_8) + RATE_T_DECAY)
#     )

#     d_m_tritium_component_9_dt = (
#         (m_tritium_component_6 / tau_6)
#         + (m_tritium_component_8 / tau_8)
#         - (m_tritium_component_9 * (((1 + epsilon_9) / tau_9) + RATE_T_DECAY))
#     )

#     d_m_tritium_component_10_dt = (
#         (1 - f_9_6) * (m_tritium_component_9 / tau_9)
#         + f_dir * (m_tritium_component_7 / tau_7)
#         + (m_tritium_component_12 / tau_12)
#         - (n_t_burn / TBE)
#         - (RATE_T_DECAY * m_tritium_component_10)
#     )

#     d_m_tritium_component_12_dt = ETA_2 * (m_tritium_component_2 / tau_2) - (
#         m_tritium_component_12 * (((1 + epsilon_12) / tau_12) + RATE_T_DECAY)
#     )

#     return (
#         d_m_tritium_component_1_dt,
#         d_m_tritium_component_2_dt,
#         d_m_tritium_component_3_dt,
#         d_m_tritium_component_4_dt,
#         d_m_tritium_component_5_dt,
#         d_m_tritium_component_6_dt,
#         d_m_tritium_component_7_dt,
#         d_m_tritium_component_8_dt,
#         d_m_tritium_component_9_dt,
#         d_m_tritium_component_10_dt,
#         d_m_tritium_component_12_dt,
#     )


# def simulate_tritium_inventory(t_end, dt):
#     """
#     Simulate tritium inventory evolution over time using Euler method.

#     Parameters:
#     t_end: End time for simulation (s)
#     dt: Time step (s)

#     Returns:
#     times: Array of time points
#     inventories: Array of inventories for each component at each time point
#     """
#     n_steps = int(t_end / dt) + 1
#     times = [i * dt for i in range(n_steps)]

#     # Initialize inventory arrays
#     inventories = {
#         "i_1": [m_tritium_component_1],
#         "i_2": [m_tritium_component_2],
#         "i_3": [m_tritium_component_3],
#         "i_4": [m_tritium_component_4],
#         "i_5": [m_tritium_component_5],
#         "i_6": [m_tritium_component_6],
#         "i_7": [m_tritium_component_7],
#         "i_8": [m_tritium_component_8],
#         "i_9": [m_tritium_component_9],
#         "i_10": [m_tritium_component_10],
#         "i_11": [m_tritium_component_11],
#         "i_12": [m_tritium_component_12],
#     }

#     # Current values
#     curr = [
#         m_tritium_component_1,
#         m_tritium_component_2,
#         m_tritium_component_3,
#         m_tritium_component_4,
#         m_tritium_component_5,
#         m_tritium_component_6,
#         m_tritium_component_7,
#         m_tritium_component_8,
#         m_tritium_component_9,
#         m_tritium_component_10,
#         m_tritium_component_11,
#         m_tritium_component_12,
#     ]

#     for step in range(1, n_steps):
#         # Calculate rates
#         rates = calculate_rates(
#             curr[0],
#             curr[1],
#             curr[2],
#             curr[3],
#             curr[4],
#             curr[5],
#             curr[6],
#             curr[7],
#             curr[8],
#             curr[9],
#             curr[11],
#         )

#         # Update inventories using Euler method
#         curr[0] += rates[0] * dt
#         curr[1] += rates[1] * dt
#         curr[2] += rates[2] * dt
#         curr[3] += rates[3] * dt
#         curr[4] += rates[4] * dt
#         curr[5] += rates[5] * dt
#         curr[6] += rates[6] * dt
#         curr[7] += rates[7] * dt
#         curr[8] += rates[8] * dt
#         curr[9] += rates[9] * dt
#         curr[11] += rates[10] * dt
#         # i_11 remains 0

#         # Store values
#         for idx, key in enumerate([
#             "i_1",
#             "i_2",
#             "i_3",
#             "i_4",
#             "i_5",
#             "i_6",
#             "i_7",
#             "i_8",
#             "i_9",
#             "i_10",
#             "i_11",
#             "i_12",
#         ]):
#             inventories[key].append(curr[idx])

#     return times, inventories
