"""
Library of Sankey plotting routine
"""

from collections.abc import Iterable
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.sankey import Sankey
from numpy import sqrt
from scipy.optimize import minimize

from process.core.io.mfile.mfile import MFile

try:
    import plotly.graph_objects as go

    PLOT_SANKEY = True
except ImportError:
    PLOT_SANKEY = False


def plot_sankey_plotly(m_file):
    if not PLOT_SANKEY:
        print(
            "\nPlotly is not installed, unable to create sankey diagram!\n"
            "Install plotly by installing the optional 'plotly' dependency "
            "e.g. \"pip install -e '.[plotly]'\""
        )
        return None
    return plotly(power_balance_sankey(m_file), m_file)


def power_balance_sankey(m_file):
    m_file = MFile(m_file)
    p_hcd_injected_total_mw = m_file.get("p_hcd_injected_total_mw", scan=-1)
    p_plasma_ohmic_mw = m_file.get("p_plasma_ohmic_mw", scan=-1)
    p_alpha_total_mw = m_file.get("p_alpha_total_mw", scan=-1)
    p_neutron_total_mw = m_file.get("p_neutron_total_mw", scan=-1)
    p_plasma_rad_mw = m_file.get("p_plasma_rad_mw", scan=-1)
    p_fw_rad_total_mw = m_file.get("p_fw_rad_total_mw", scan=-1)
    p_fw_alpha_mw = p_alpha_total_mw * (
        1 - m_file.get("f_p_alpha_plasma_deposited", scan=-1)
    )
    p_blkt_nuclear_heat_total_mw = m_file.get("p_blkt_nuclear_heat_total_mw", scan=-1)

    # Define node labels (linearized flow)
    labels = [
        "H&CD injector",  # 0
        "Ohmic",  # 1
        "Plasma Fusion Power",  # 2
        "Alpha particles",  # 3
        "Neutrons",  # 4
        "Radiation",  # 5
        "First Wall",  # 6
        "Blanket",  # 7
        "Divertor",  # 8
        "FW+Blkt",  # 9
        "Primary Thermal",  # 10
        "Turbine",  # 11
        "Gross Electric",  # 12
        "Net Electric",  # 13
        "HCD Electric Power",  # 14
        "HCD electric losses",  # 15
        "Core systems",  # 16
        "Cryo plant",  # 17
        "Base plant load",  # 18
        "TF power supplies",  # 19
        "PF power supplies",  # 20
        "Vacuum pumps",  # 21
        "Tritium plant",  # 22
        "Coolant pumps electric",  # 23
        "Coolant pump electric losses",  # 24
        "Divertor pump",  # 25
        "FW+Blkt pumps",  # 26
        "Shield pump",  # 27
        "Shield",  # 28
        "Secondary heat",  # 29
        "TF nuclear heat",  # 30
        "H&CD & Diagnostics",  # 31
        "Total Secondary Heat",  # 32
        "Turbine Loss",  # 33
        "Blanket neutron multiplication",  # 34
    ]

    # Define links (source, target, value) for a more linear flow
    sources = [
        0,  # 0: H&CD to Fusion
        1,  # 1: Ohmic to Fusion
        2,  # 2: Fusion to Alpha
        2,  # 3: Fusion to Neutrons
        2,  # 4: Fusion to Radiation
        3,  # 5: Alpha to First Wall
        4,  # 6: Neutrons to Blanket
        5,  # 7: Radiation to First Wall
        4,  # 8: Neutrons to Divertor
        5,  # 9: Radiation to Divertor
        6,  # 10: First Wall to FW+Blkt
        7,  # 11: Blanket to FW+Blkt
        8,  # 12: Divertor to FW+Blkt
        9,  # 13: FW+Blkt to Primary Thermal
        10,  # 14: Primary Thermal to Turbine
        11,  # 15: Turbine to Gross Electric
        12,  # 16: Gross Electric to Net Electric
        12,  # 17: Gross Electric to HCD Electric Power
        14,  # 18: HCD Electric Power to HCD electric losses
        14,  # 19: HCD Electric Power to H&CD
        12,  # 20: Gross Electric to Core systems
        16,  # 21: Core systems to Cryo plant
        16,  # 22: Core systems to Base plant load
        16,  # 23: Core systems to TF coils
        16,  # 24: Core systems to PF coils
        16,  # 25: Core systems to Vacuum pumps
        16,  # 26: Core systems to Tritium plant
        12,  # 27: Gross Electric to Coolant pumps electric
        23,  # 28: Coolant pumps electric to Coolant pump electric losses
        23,  # 29: Coolant pumps electric to Divertor pump
        23,  # 30: Coolant pumps electric to FW+Blkt pumps
        26,  # 31: FW+Blkt pumps to FW+Blkt
        25,  # 32: Divertor pump to Divertor
        23,  # 33: Coolant pumps electric to Shield pump
        27,  # 34: Shield pump to Shield
        28,  # 35: Shield to primary thermal
        4,  # 36: Neutrons to shield
        17,  # 37: Cryo plant to secondary heat
        18,  # 38: Base plant load to secondary heat
        19,  # 39: TF coils to secondary heat
        20,  # 40: PF coils to secondary heat
        21,  # 41: Vacuum pumps to secondary heat
        22,  # 42: Tritium plant to secondary heat
        4,  # 43: Neutrons to tf
        30,  # 44: TF nuclear heat to secondary heat
        15,  # 45: HCD electric losses to secondary heat
        24,  # 46: Coolant pumps electric to secondary heat
        6,  # 47: FW pump to primary heat, Should only show if FW and Bkt pumps are separate
        7,  # 48: Blkt pump to primary heat, Should only show if FW and Blkt pumps are separate
        2,  # 49 Should show in beams are present
        2,  # 50:  Should show in beams are present
        4,  # 51 Neutrons to CP shield, should only show if CP shield is present
        2,  # 52 Plasma separatrix power to divertor
        8,  # 53 Divertor secondary heat,
        28,  # 54 Shield secondary heat
        4,  # 55 Neutron power to H&CD & Diagnostics
        5,  # 56: Radiation to H&CD & Diagnostics
        29,  # 57: Total Secondary Heat
        31,  # 58: H&CD & Diagnostics secondary heat
        11,  # 59: Turbine Loss
        4,  # 60: FW nuclear heat
        3,  # 61: Alpha particles back to plasma
        34,  # 62: Blanket neutron multiplication
    ]
    targets = [
        2,  # 0: H&CD to Fusion
        2,  # 1: Ohmic to Fusion
        3,  # 2: Fusion to Alpha
        4,  # 3: Fusion to Neutrons
        5,  # 4: Fusion to Radiation
        6,  # 5: Alpha to First Wall
        7,  # 6: Neutrons to Blanket
        6,  # 7: Radiation to First Wall
        8,  # 8: Neutrons to Divertor
        8,  # 9: Radiation to Divertor
        9,  # 10: First Wall to FW+Blkt
        9,  # 11: Blanket to FW+Blkt
        10,  # 12: Divertor to FW+Blkt
        10,  # 13: FW+Blkt to Primary Thermal
        11,  # 14: Primary Thermal to Turbine
        12,  # 15: Turbine to Gross Electric
        13,  # 16: Gross Electric to Net Electric
        14,  # 17: Gross Electric to HCD Electric Power
        15,  # 18: HCD Electric Power to HCD electric losses
        0,  # 19: HCD Electric Power to H&CD
        16,  # 20: Gross Electric to Core systems
        17,  # 21: Core systems to Cryo plant
        18,  # 22: Core systems to Base plant load
        19,  # 23: Core systems to TF coils
        20,  # 24: Core systems to PF coils
        21,  # 25: Core systems to Vacuum pumps
        22,  # 26: Core systems to Tritium plant
        23,  # 27: Gross Electric to Coolant pumps electric
        24,  # 28: Coolant pumps electric to Coolant pump electric losses
        25,  # 29: Coolant pumps electric to Divertor pump
        26,  # 30: Coolant pumps electric to FW+Blkt pumps
        9,  # 31: FW+Blkt pumps to FW+Blkt
        8,  # 32: Divertor pump to Divertor
        27,  # 33: Coolant pumps electric to Shield pump
        28,  # 34: Shield pump to Shield
        10,  # 35: Shield to primary thermal
        28,  # 36: Neutrons to shield
        29,  # 37: Cryo plant to secondary heat
        29,  # 38: Base plant load to secondary heat
        29,  # 39: TF coils to secondary heat
        29,  # 40: PF coils to secondary heat
        29,  # 41: Vacuum pumps to secondary heat
        29,  # 42: Tritium plant to secondary heat
        30,  # 43: Neutrons to tf
        29,  # 44: TF nuclear heat to secondary heat
        29,  # 45: HCD electric losses to secondary heat
        29,  # 46: Coolant pumps electric to secondary heat
        9,  # 47: FW pump to primary heat, Should only show if FW and Bkt pumps are separate
        9,  # 48: Blkt pump to primary heat, Should only show if FW and Blkt pumps are separate
        6,  # 49 Should show in beams are present
        6,  # 50:  Should show in beams are present
        28,  # 51 Neutrons to CP shield, should only show if CP shield is present
        8,  # 52 Plasma separatrix power to divertor
        29,  # 53 Divertor secondary heat,
        29,  # 54 Shield secondary heat
        31,  # 55 Neutron power to H&CD & Diagnostics
        31,  # 56: Radiation to H&CD & Diagnostics
        32,  # 57: Total Secondary Heat
        32,  # 58: H&CD & Diagnostics secondary heat
        33,  # 59: Turbine Loss
        6,  # 60: FW nuclear heat
        2,  # 61: Alpha particles back to plasma
        7,  # 62: Blanket neutron multiplication
    ]
    values = [
        p_hcd_injected_total_mw,  # 0
        p_plasma_ohmic_mw,  # 1
        p_alpha_total_mw,  # 2
        p_neutron_total_mw,  # 3
        p_plasma_rad_mw,  # 4
        p_fw_alpha_mw,  # 5
        p_blkt_nuclear_heat_total_mw
        - m_file.get("p_blkt_multiplication_mw", scan=-1),  # 6
        p_fw_rad_total_mw,  # 7
        m_file.get("p_div_nuclear_heat_total_mw", scan=-1),  # 8
        m_file.get("p_div_rad_total_mw", scan=-1),  # 9
        m_file.get("p_fw_heat_deposited_mw", scan=-1),  # 10
        m_file.get("p_blkt_heat_deposited_mw", scan=-1),  # 11
        m_file.get("p_div_heat_deposited_mw", scan=-1),  # 12
        m_file.get("p_fw_blkt_heat_deposited_mw", scan=-1),  # 13
        m_file.get("p_plant_primary_heat_mw", scan=-1),  # 14
        m_file.get("p_plant_electric_gross_mw", scan=-1),  # 15
        m_file.get("p_plant_electric_net_mw", scan=-1),  # 16
        m_file.get("p_hcd_electric_total_mw", scan=-1),  # 17
        m_file.get("p_hcd_electric_loss_mw", scan=-1),  # 18
        p_hcd_injected_total_mw,  # 19
        m_file.get("p_plant_core_systems_elec_mw", scan=-1),  # 20
        m_file.get("p_cryo_plant_electric_mw", scan=-1),  # 21
        m_file.get("p_plant_electric_base_total_mw", scan=-1),  # 22
        m_file.get("p_tf_electric_supplies_mw", scan=-1),  # 23
        m_file.get("p_pf_electric_supplies_mw", scan=-1),  # 24
        m_file.get("vachtmw", scan=-1),  # 25
        m_file.get("p_tritium_plant_electric_mw", scan=-1),  # 26
        m_file.get("p_coolant_pump_elec_total_mw", scan=-1),  # 27
        m_file.get("p_coolant_pump_loss_total_mw", scan=-1),  # 28
        m_file.get("p_div_coolant_pump_mw", scan=-1),  # 29
        m_file.get("p_fw_blkt_coolant_pump_mw", scan=-1),  # 30
        m_file.get("p_fw_blkt_coolant_pump_mw", scan=-1),  # 31
        m_file.get("p_div_coolant_pump_mw", scan=-1),  # 32
        m_file.get("p_shld_coolant_pump_mw", scan=-1),  # 33
        m_file.get("p_shld_coolant_pump_mw", scan=-1),  # 34
        m_file.get("p_shld_heat_deposited_mw", scan=-1),  # 35
        m_file.get("p_shld_nuclear_heat_mw", scan=-1),  # 36
        m_file.get("p_cryo_plant_electric_mw", scan=-1),  # 37
        m_file.get("p_plant_electric_base_total_mw", scan=-1),  # 38
        m_file.get("p_tf_electric_supplies_mw", scan=-1),  # 39
        m_file.get("p_pf_electric_supplies_mw", scan=-1),  # 40
        m_file.get("vachtmw", scan=-1),  # 41
        m_file.get("p_tritium_plant_electric_mw", scan=-1),  # 42
        m_file.get("p_tf_nuclear_heat_mw", scan=-1),  # 43
        m_file.get("p_tf_nuclear_heat_mw", scan=-1),  # 44
        m_file.get("p_hcd_electric_loss_mw", scan=-1),  # 45
        m_file.get("p_coolant_pump_loss_total_mw", scan=-1),  # 46
        #
        # Should only show if FW and Bkt pumps are seperate
        m_file.get("p_fw_coolant_pump_mw", scan=-1),  # 47
        m_file.get("p_blkt_coolant_pump_mw", scan=-1),  # 48
        #
        # Should show in beams are present
        m_file.get("p_beam_shine_through_mw", scan=-1),  # 49
        m_file.get("p_beam_orbit_loss_mw", scan=-1),  # 50
        #
        # Neutrons to CP shield, should only show if CP shield is present
        m_file.get("p_cp_shield_nuclear_heat_mw", scan=-1),  # 51
        #
        m_file.get("p_plasma_separatrix_mw", scan=-1),  # 52
        m_file.get("p_div_secondary_heat_mw", scan=-1),  # 53
        m_file.get("p_shld_secondary_heat_mw", scan=-1),  # 54
        m_file.get("p_fw_hcd_nuclear_heat_mw", scan=-1),  #
        m_file.get("p_fw_hcd_rad_total_mw", scan=-1),  # 56
        m_file.get("p_plant_secondary_heat_mw", scan=-1),  # 57
        m_file.get("p_hcd_secondary_heat_mw", scan=-1),  # 58
        m_file.get("p_turbine_loss_mw", scan=-1),  # 59
        m_file.get("p_fw_nuclear_heat_total_mw", scan=-1),  # 60
        #
        # Alpha particles back to plasma
        p_alpha_total_mw * m_file.get("f_p_alpha_plasma_deposited", scan=-1),  # 61
        m_file.get("p_blkt_multiplication_mw", scan=-1),
    ]

    # Define colors for each node (hex or rgba)
    node_colors = [
        "#1f77b4",  # 0: H&CD injector
        "#ff7f0e",  # 1: Ohmic
        "#2ca02c",  # 2: Plasma Fusion Power
        "#d62728",  # 3: Alpha particles
        "#9467bd",  # 4: Neutrons
        "#8c564b",  # 5: Radiation
        "#e377c2",  # 6: First Wall
        "#7f7f7f",  # 7: Blanket
        "#bcbd22",  # 8: Divertor
        "#17becf",  # 9: FW+Blkt
        "#aec7e8",  # 10: Primary Thermal
        "#ffbb78",  # 11: Turbine
        "#98df8a",  # 12: Gross Electric
        "#ff9896",  # 13: Net Electric
        "#c5b0d5",  # 14: HCD Electric Power
        "#c49c94",  # 15: HCD electric losses
        "#f7b6d2",  # 16: Core systems
        "#c7c7c7",  # 17: Cryo plant
        "#dbdb8d",  # 18: Base plant load
        "#9edae5",  # 19: TF coils
        "#393b79",  # 20: PF coils
        "#637939",  # 21: Vacuum pumps
        "#8c6d31",  # 22: Tritium plant
        "#843c39",  # 23: Coolant pumps electric
        "#7b4173",  # 24: Coolant pump electric losses
        "#5254a3",  # 25: Divertor pump
        "#6b6ecf",  # 26: FW+Blkt pumps
        "#b5cf6b",  # 27: Shield pump
        "#cedb9c",  # 28: Shield
        "#9c9ede",  # 29: Secondary heat
        "#e7ba52",  # 30: TF nuclear heat
        "#ad494a",  # 31: H&CD & Diagnostics
        "#a55194",  # 32: Total Secondary Heat
        "#393b79",  # 33: Turbine Loss
        "#637939",  # 34: Blanket neutron multiplication
    ]

    # Assign link colors to match their source node
    link_colors = [node_colors[src] for src in sources]

    # Add value labels to the links
    value_labels = [f"{v:.3f} MW" for v in values]

    return {
        "type": "sankey",
        "node": {
            "pad": 30,
            "thickness": 20,
            "line": {"color": "black", "width": 0.5},
            "label": labels,
            "color": node_colors,
        },
        "link": {
            "source": sources,
            "target": targets,
            "value": values,
            "label": value_labels,
            "color": link_colors,
        },
    }


def plotly(sankey_dict, m_file):
    fig = go.Figure(data=[sankey_dict])

    fig.update_layout({
        "title_text": "Fusion Power Balance Sankey Diagram",
        "font_size": 7,
        "autosize": True,
        "margin": {"l": 40, "r": 40, "t": 40, "b": 40},
    })
    # Strip 'MFILE' from the filename for the HTML output
    html_output_path = (
        Path(m_file)
        .with_stem(Path(m_file).stem.replace("MFILE", "plotly_sankey"))
        .with_suffix(".html")
    )
    fig.write_html(str(html_output_path))
    print(f"Interactive Sankey diagram saved to {html_output_path}")
    return fig


class SuperSankey(Sankey):
    """
    Originally from Bluemira

    A sub-class of the Sankey diagram class from matplotlib, which is capable
    of connecting two blocks, instead of just one. This is done using a cute
    sledgehammer approach, using optimisation. Basically, the Sankey object
    is quite complex, and it makes it very hard to calculate the exact lengths
    required to connect two sub-diagrams.
    """

    def add(
        self,
        patchlabel: str = "",
        flows: Iterable[float] | None = None,
        orientations: Iterable[float] | None = None,
        labels: str | list[str | None] | None = "",
        trunklength: float = 1.0,
        pathlengths: float | list[float | None] = 0.25,
        prior: int | None = None,
        future: int | None = None,
        connect: tuple[int, int] | list[tuple[int, int]] = (0, 0),
        rotation: float = 0,
        **kwargs,
    ):
        __doc__ = super().__doc__  # noqa: F841, A001
        # Here we first check if the "add" method has received arguments that
        # the Sankey class can't handle.
        if future is None:
            # There is only one connection, Sankey knows how to do this
            super().add(
                patchlabel,
                flows,
                orientations,
                labels,
                trunklength,
                pathlengths,
                prior,
                connect,
                rotation,
                **kwargs,
            )
        else:
            # There are two connections, use new method
            self._double_connect(
                patchlabel,
                flows,
                orientations,
                labels,
                trunklength,
                pathlengths,
                prior,
                future,
                connect,
                rotation,
                **kwargs,
            )

    def _double_connect(
        self,
        patchlabel: str,
        flows: Iterable[float] | None,
        orientations: Iterable[float] | None,
        labels: str | list[str | None] | None,
        trunklength: float,
        pathlengths: list[float],
        prior: int | None,
        future: int | None,
        connect: list[tuple[int, int]],
        rotation: float,
        **kwargs,
    ):
        """
        Handles two connections in a Sankey diagram.

        Parameters
        ----------
        future:
            The index of the diagram to connect to
        connect:
            The list of (int, int) connections.
            - connect[0] is a (prior, this) tuple indexing the flow of the
            prior diagram and the flow of this diagram to connect.
            - connect[1] is a (future, this) tuple indexing of the flow of the
            future diagram and the flow of this diagram to connect.

        See Also
        --------
        Sankey.add for a full description of the various args and kwargs

        """
        # Get the optimum deltas
        dx, dy = self._opt_connect(
            flows, orientations, prior, future, connect, trunklength=trunklength
        )
        # Replace
        pathlengths[0] = dx
        pathlengths[-1] = dy
        self.add(
            patchlabel=patchlabel,
            labels=labels,
            flows=flows,
            orientations=orientations,
            prior=prior,
            connect=connect[0],
            trunklength=trunklength,
            pathlengths=pathlengths,
            rotation=rotation,
            facecolor=kwargs.get("facecolor"),
        )

    def _opt_connect(
        self,
        flows: Iterable[float] | None,
        orient: Iterable[float] | None,
        prior: int | None,
        future: int | None,
        connect: list[tuple[int, int]],
        trunklength: float,
    ) -> tuple[float, float]:
        """
        Optimises the second connection between Sankey diagrams.

        Returns
        -------
        dx:
            The x pathlength to use to match the tips
        dy:
            The y pathlength to use to match the tips

        Notes
        -----
        This is because Sankey is very complicated, and makes it hard to work
        out the positions of things prior to adding them to the diagrams.
        Because we are bizarrely using a plotting function as a minimisation
        objective, we need to make sure we clean the plot on every call.
        """
        future_index, this_f_index = connect[1]
        labels = [None] * len(flows)
        pathlengths = [0.0] * len(flows)

        # Make a local copy of the Sankey.extent attribute to override any
        # modifications during optimisation
        extent = deepcopy(self.extent)

        def minimise_dxdy(x_opt):
            """
            Minimisation function for the spatial difference between the target
            tip and the actual tip.

            Parameters
            ----------
            x_opt: array_like
                The vector of d_x, d_y delta-vectors to match tip positions

            Returns
            -------
            delta: float
                The sum of the absolute differences
            """
            tip2 = self.diagrams[future].tips[future_index]
            pathlengths[0] = x_opt[0]
            pathlengths[-1] = x_opt[1]
            self.add(
                trunklength=trunklength,
                pathlengths=pathlengths,
                flows=flows,
                prior=prior,
                connect=connect[0],
                orientations=orient,
                labels=labels,
                facecolor="#00000000",
            )
            new_tip = self.diagrams[-1].tips[this_f_index].copy()
            # Clean sankey plot
            self.diagrams.pop()
            self.ax.patches[-1].remove()
            return np.sum(np.abs(tip2 - new_tip))

        x0 = np.zeros(2)
        result = minimize(minimise_dxdy, x0, method="SLSQP")
        self.extent = extent  # Finish clean-up
        return result.x


def plot_sankey(
    mfilename=Path("MFILE.DAT"), format_: str = "pdf"
):  # Plot simplified power flow Sankey Diagram
    # ------------------------------- Pulling values from the MFILE -------------------------------
    mfilename = Path(mfilename)
    m_file = MFile(mfilename)

    variables = [
        # Used in [PLASMA]
        "p_fusion_total_mw",  # Fusion Power (MW)
        "p_hcd_injected_total_mw",  # Total auxiliary injected Power (MW)
        "p_plasma_ohmic_mw",  # Ohmic heating Power (MW)
        # Used in [DEPOSITION]
        "p_plasma_rad_mw",  # Total radiation Power (MW)
        "f_ster_div_single",  # Area fraction taken up by divertor
        "2*f_ster_div_single",  # Area fraction taken up by double null divertor
        "f_a_fw_outboard_hcd",  # Area fraction covered by HCD and diagnostics
        "p_plasma_separatrix_mw",  # power to conducted to the divertor region (MW)
        "p_div_nuclear_heat_total_mw",  # nuclear heating in the divertor (MW)
        "p_fw_nuclear_heat_total_mw",  # nuclear heating in the first wall (MW)
        "p_blkt_nuclear_heat_total_mw",  # nuclear heating in the blanket (MW)
        "p_shld_nuclear_heat_mw",  # nuclear heating in the shield (MW)
        "p_cp_shield_nuclear_heat_mw",  # nuclear heating in the CP shield (MW)
        "p_blkt_multiplication_mw",  # Blanket energy multiplication (MW)
        "p_alpha_total_mw",  # Alpha power (MW)
        "f_p_alpha_plasma_deposited",  # Fraction of alpha power deposited in plasma
        "itart",  # switch for spherical tokamak (ST) models
        # Used in [BLANKETSETC]
        "p_fw_blkt_heat_deposited_mw",  # Heat for electricity (MW)
        "p_fw_blkt_coolant_pump_mw",  # 1st wall & blanket pumping (MW)
        # Used in [PRIMARY]
        "p_plant_electric_gross_mw",  # gross electric power (MW)
        # Used in [NET]
        "p_plant_electric_net_mw",  # net electric power (MW)
        # Used in [RECIRC]
        "p_cryo_plant_electric_mw",  # cryogenic plant power (MW)
        "fachtmw",  # facility heat removal (MW)
        "p_tf_electric_supplies_mw",  # total steady state TF coil AC power demand (MW)
        "p_tritium_plant_electric_mw",  # power required for tritium processing (MW)
        "vachtmw",  # vacuum pump power (MW)
        "p_pf_electric_supplies_mw",  # Total mean wall plug power for PFC & CS (MW)
        "p_hcd_electric_total_mw",  # injector wall plug power (MW)
        "p_coolant_pump_elec_total_mw",  # heat transport system electrical pump power (MW)
        "p_cp_coolant_pump_elec",  # pumping power
    ]
    (
        p_fusion_total_mw,
        p_hcd_injected_total_mw,
        p_plasma_ohmic_mw,
        p_plasma_rad_mw,
        f_ster_div_single,
        fdiv_2,
        f_a_fw_outboard_hcd,
        p_plasma_separatrix_mw,
        p_div_nuclear_heat_total_mw,
        p_fw_nuclear_heat_total_mw,
        p_blkt_nuclear_heat_total_mw,
        p_shld_nuclear_heat_mw,
        p_cp_shield_nuclear_heat_mw,
        p_blkt_multiplication_mw,
        p_alpha_total_mw,
        f_p_alpha_plasma_deposited,
        itart,
        p_fw_blkt_heat_deposited_mw,
        p_fw_blkt_coolant_pump_mw,
        p_plant_electric_gross_mw,
        p_plant_electric_net_mw,
        p_cryo_plant_electric_mw,
        fachtmw,
        p_tf_electric_supplies_mw,
        p_tritium_plant_electric_mw,
        vachtmw,
        p_pf_electric_supplies_mw,
        p_hcd_electric_total_mw,
        p_coolant_pump_elec_total_mw,
        p_cp_coolant_pump_elec,
    ) = m_file.get_variables(*variables, scan=-1)

    p_cp_coolant_pump_elec_mw = p_cp_coolant_pump_elec / 1e6

    # Total Power in plasma (MW)
    totalplasma = p_fusion_total_mw + p_hcd_injected_total_mw + p_plasma_ohmic_mw

    if fdiv_2 > 0:  # Takes into account old MFILE representation of double null divertor
        f_ster_div_single = fdiv_2

    # Radiation deposited on the divertor (MW)
    p_div_rad_total_mw = p_plasma_rad_mw * f_ster_div_single
    # Radiation deposited on HCD and diagnostics (MW)
    p_fw_hcd_rad_total_mw = p_plasma_rad_mw * f_a_fw_outboard_hcd
    # Radiation deposited in the blanket (MW)
    p_fw_rad_total_mw = p_plasma_rad_mw - p_div_rad_total_mw - p_fw_hcd_rad_total_mw
    # Alpha power hitting 1st wall (MW)
    p_fw_alpha_mw = p_alpha_total_mw * (1 - f_p_alpha_plasma_deposited)

    # Power deposited on divertor (MW)
    totaldivetc = (
        p_plasma_separatrix_mw + p_div_nuclear_heat_total_mw + p_div_rad_total_mw
    )
    # Power deposited on Blanket (MW)
    totalblktetc = (
        p_fw_nuclear_heat_total_mw
        + p_blkt_nuclear_heat_total_mw
        + p_shld_nuclear_heat_mw
        + p_fw_rad_total_mw
        + p_fw_alpha_mw
        - p_blkt_multiplication_mw
    )

    if itart == 0:
        # Power deposited in CP (MW) (None here)
        totalcpetc = 0.0
    elif itart == 1:
        # Power deposited in CP (MW)
        totalcpetc = p_cp_shield_nuclear_heat_mw

    # Heat - pumping power (MW)
    pthermmw_p = p_fw_blkt_heat_deposited_mw - p_fw_blkt_coolant_pump_mw

    # Recirculating power (MW)
    p_plant_electric_recirc_mw = p_plant_electric_gross_mw - p_plant_electric_net_mw

    # Energy required for rest of power plant (MW)
    p_plant_core_systems_elec_mw = (
        p_cryo_plant_electric_mw
        + fachtmw
        + p_tf_electric_supplies_mw
        + p_tritium_plant_electric_mw
        + vachtmw
        + p_pf_electric_supplies_mw
        + p_cp_coolant_pump_elec_mw
    )

    # -------------------------------- Visual Settings ------------------------------------

    plt.rcParams.update({"font.size": 9})  # Setting font size to 9
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], frameon=False)
    sankey = SuperSankey(
        ax=ax, unit="MW", margin=0.0, format="%1.0f", scale=1.0 / (totalplasma)
    )
    trunk = 0.7
    len1 = 0.5
    len2 = 0.8
    # --------------------------------------- PLASMA - 0 --------------------------------------

    # Fusion power, Injected power + ohmic power, - total plasma power
    plasma = [
        p_fusion_total_mw,
        p_hcd_injected_total_mw + p_plasma_ohmic_mw,
        -totalplasma,
    ]
    sankey.add(
        flows=plasma,
        orientations=[0, -1, 0],  # [right(in), down(in), right(out)]
        pathlengths=[
            len1,
            len2,
            -0.1 + len1,
        ],  # 'Plasma Heating' adjust
        trunklength=trunk,
        labels=["Fusion Power", None, "Plasma"],
    )

    # --------------------------------- ENERGY DEPOSITION - 1 ---------------------------------

    # Plasma power, - divertor deposited power, - blanket deposited power
    deposition = [totalplasma, -totalblktetc - totaldivetc - totalcpetc]
    # Check if difference >2 between plasma and divertor + blanket
    if sqrt(sum(deposition) ** 2) > 2:
        print(
            "\ncomponents power balance difference =",
            totalplasma - totaldivetc - totalblktetc - totalcpetc,
        )
    sankey.add(
        flows=deposition,
        orientations=[0, 0],  # [right(in), up(in), right(out)]
        prior=0,  # PLASMA
        connect=(2, 0),  # Plasma --> None
        pathlengths=[0.2, len2],  # 'Plasma Heating' adjust
        trunklength=trunk,
        labels=[None, "Blanket/etc."],
    )

    # -------------------------------------- BLANKET - 2 --------------------------------------

    # Blanket deposited power, blanket energy multiplication, - primary heat
    blanketsetc = [
        totalblktetc + totaldivetc + totalcpetc,
        p_blkt_multiplication_mw,
        -pthermmw_p - totaldivetc - totalcpetc - p_shld_nuclear_heat_mw,
    ]
    # Check if difference >2 between primary heat and blanket + blanket multiplication
    if sqrt(sum(blanketsetc) ** 2) > 2:
        print(
            "blankets etc. power balance",
            totalblktetc + p_blkt_multiplication_mw,
            -pthermmw_p - p_shld_nuclear_heat_mw,
        )
    sankey.add(
        flows=blanketsetc,
        orientations=[0, -1, 0],  # [right(in), down(in), right(out)]
        prior=1,  # DEPOSITION
        connect=(1, 0),  # Blanket/etc. --> None
        pathlengths=[len1, len1 / 2, 0.0],
        trunklength=trunk,
        labels=[None, "Energy Mult.", "Primary Heat"],
    )

    # ------------------------------------- HEAT LOSS - 3 -------------------------------------

    # Primary heat, -Gross electric power, -difference (loss)
    primary = [
        pthermmw_p + totaldivetc + totalcpetc + p_shld_nuclear_heat_mw,
        -p_plant_electric_gross_mw,
        -pthermmw_p
        + p_plant_electric_gross_mw
        - totaldivetc
        - totalcpetc
        - p_shld_nuclear_heat_mw,
    ]
    sankey.add(
        flows=primary,
        orientations=[0, -1, 0],  # [right(in), down(out), right(out)]
        prior=2,  # BLANKETSETC
        connect=(2, 0),  # Primary Heat --> None
        pathlengths=[len2 / 4, len2, len1 / 2],
        trunklength=trunk,
        labels=[None, "Gross electric", "Losses"],
    )

    # ------------------------------------ ELECTRICITY - 4 ------------------------------------

    # If net electric is +ve or -ve changes the flow organisation
    if p_plant_electric_net_mw >= 0:  # net electric is +ve
        # Gross electric power, -net electric power, -recirculated power
        net = [
            p_plant_electric_gross_mw,
            -p_plant_electric_net_mw,
            -p_plant_electric_recirc_mw,
        ]
        sankey.add(
            flows=net,
            orientations=[0, 0, -1],  # [down(in), down(out), left(out)]
            prior=3,  # PRIMARY
            connect=(1, 0),  # Gross electric --> None
            pathlengths=[len2 / 4, len1 / 2, 3 * len1],
            trunklength=trunk,
            labels=[None, "Net elec.", "Recirc. Power"],
        )
    elif p_plant_electric_net_mw < 0:  # net electric is -ve
        # Gross electric power, -net electric power, -recirculated power
        net = [
            -p_plant_electric_net_mw,
            p_plant_electric_gross_mw,
            -p_plant_electric_recirc_mw,
        ]
        sankey.add(
            flows=net,
            orientations=[0, -1, 0],  # [left(in), down(in), left(out)]
            prior=3,  # PRIMARY
            connect=(1, 1),  # Gross electric --> None
            pathlengths=[len1 / 2, 2 * len1, len1],
            trunklength=trunk,
            labels=["Net elec.", None, "Recirc. Power"],
        )

    # -------------------------------- RECIRCULATING POWER - 5 --------------------------------

    # Recirculated power, -Core Systems, -Heating System
    recirc = [
        p_plant_electric_recirc_mw,
        -p_plant_core_systems_elec_mw - p_coolant_pump_elec_total_mw,
        -p_hcd_electric_total_mw + p_cp_coolant_pump_elec_mw,
    ]
    # Check if difference >2 between recirculated power and the output sum
    if sum(recirc) ** 2 > 2:
        print(
            "Recirc. Power Balance",
            p_plant_electric_recirc_mw,
            -p_plant_core_systems_elec_mw
            + p_cp_coolant_pump_elec_mw
            - p_hcd_electric_total_mw
            - p_coolant_pump_elec_total_mw,
        )
    sankey.add(
        flows=recirc,
        orientations=[0, 1, 0],  # [left(in), down(out), left(out)]
        prior=4,  # NET
        connect=(2, 0),  # Recirc. Power --> None
        pathlengths=[0.1, len1 / 2, len2],
        trunklength=trunk * 1.2,
        labels=[None, "Core Systems", "Heating System"],
    )

    # --------------------------------------- LOSSES - 6 --------------------------------------

    # HCD: Heating system, -Plasma heating, -losses
    hcd = [
        p_hcd_electric_total_mw - p_cp_coolant_pump_elec_mw,
        -p_hcd_injected_total_mw,
        -p_hcd_electric_total_mw + p_hcd_injected_total_mw + p_cp_coolant_pump_elec_mw,
    ]
    sankey.add(
        flows=hcd,
        orientations=[0, 0, -1],  # [left(in), up(out), left(out)]
        prior=5,  # RECIRC
        future=0,
        connect=[(2, 0), (1, 2)],  # Heating System --> None
        pathlengths=[None, len1, None],  # 'Plasma Heating' adjust
        trunklength=trunk,
        labels=[None, "Losses", "Plasma Heating"],
    )

    # Collecting Sankey diagram and applying a condensed layout
    diagrams = sankey.finish()
    fig.tight_layout()

    # --------------------------------------- Label Positioning ---------------------------------------

    # Munipulating the positioning of the branch labels
    # -ve to left and down; +ve to right and up
    # pos[0] = x-axis; pos[1] = y-axis
    for d in diagrams:
        for y, t in enumerate(d.texts):
            pos = tuple(np.ndarray.tolist(d.tips[y]))
            t.set_position(pos)
            if t == diagrams[0].texts[0]:  # Fusion Power
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] - 0.35,
                    pos[1] + 0.5 * (p_fusion_total_mw / totalplasma) + 0.2,
                ))
            if t == diagrams[0].texts[2]:  # Plasma
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.25, pos[1]))
            if t == diagrams[1].texts[1]:  # Blanket/etc.
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.2, pos[1]))
            if t == diagrams[2].texts[1]:  # Energy Mult.
                t.set_position((pos[0], pos[1] - 0.3))
            if t == diagrams[2].texts[2]:  # Primary Heat
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.25, pos[1]))
            if t == diagrams[3].texts[1]:  # Gross Electric
                t.set_horizontalalignment("right")
                t.set_position((
                    pos[0] - 0.5 * (p_plant_electric_gross_mw / totalplasma) - 0.1,
                    pos[1] + 0.1,
                ))
            if t == diagrams[3].texts[2]:  # Losses
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.2, pos[1]))
            if p_plant_electric_net_mw >= 1:
                if t == diagrams[4].texts[1]:  # Net electric
                    t.set_horizontalalignment("center")
                    t.set_position((pos[0], pos[1] - 0.2))
            elif (
                p_plant_electric_net_mw < 1 and t == diagrams[4].texts[0]
            ):  # Net electric
                t.set_horizontalalignment("left")
                t.set_position((pos[0] + 0.2, pos[1]))
            if t == diagrams[4].texts[2]:  # Recirc. Power
                if p_plant_electric_net_mw >= 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (p_plant_electric_recirc_mw / totalplasma) + 0.2,
                    ))
                elif p_plant_electric_net_mw < 1:
                    t.set_horizontalalignment("left")
                    t.set_position((pos[0] + 0.2, pos[1]))
            if t == diagrams[5].texts[1]:  # Core Systems
                t.set_position((pos[0], pos[1] - 0.2))
            if t == diagrams[5].texts[2]:  # Heating System
                if p_plant_electric_net_mw >= 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (p_hcd_electric_total_mw / totalplasma) + 0.2,
                    ))
                if p_plant_electric_net_mw < 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (p_hcd_electric_total_mw / totalplasma) + 0.2,
                    ))
            if t == diagrams[6].texts[2]:  # Plasma Heating
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] + 0.5 * (p_hcd_injected_total_mw / totalplasma) + 0.1,
                    pos[1] - 0.05,
                ))
            if t == diagrams[6].texts[1]:  # Losses
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] + 0.15,
                    pos[1]
                    - 0.5
                    * ((p_hcd_electric_total_mw - p_hcd_injected_total_mw) / totalplasma)
                    - 0.2,
                ))

    # Get directory of mfile
    fig.savefig(mfilename.parent / f"SankeyPowerFlow.{format_}")

    plt.show()
    return fig
