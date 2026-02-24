"""
Library of Sankey plotting routine
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.sankey import Sankey
from numpy import sqrt

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
    html_output_path = m_file.with_stem(
        m_file.stem.replace("MFILE", "plotly_sankey")
    ).with_suffix(".html")
    fig.write_html(str(html_output_path))
    print(f"Interactive Sankey diagram saved to {html_output_path}")
    return fig


def plot_full_sankey(
    mfilename="MFILE.DAT",
):  # Plots the power flow from PROCESS as a Sankey Diagram
    # ------------------------------- Pulling values from the MFILE -------------------------------

    m_file = MFile(mfilename)
    variables = [
        # Used in [PLASMA]
        "p_fusion_total_mw",  # Fusion Power (MW)
        "p_hcd_injected_total_mw",  # Total auxiliary injected Power (MW)
        "p_plasma_ohmic_mw",  # Ohmic heating Power (MW)
        "p_neutron_total_mw",  # Neutron fusion power (MW)
        "p_non_alpha_charged_mw",  # Non-alpha charged particle power (MW)
        "p_alpha_total_mw",  # Alpha power (MW)
        # Used in [NEUTRONICS]
        "p_blkt_multiplication_mw",  # Energy multiplication in blanket (MW)
        "p_blkt_nuclear_heat_total_mw",  # Total Nuclear heating in the blanket (MW)
        "p_div_nuclear_heat_total_mw",  # Nuclear heating in the divertor (MW)
        "p_fw_nuclear_heat_total_mw",  # Nuclear heating in the first wall (MW)
        "p_shld_nuclear_heat_mw",  # Nuclear heating in the shield (MW)
        "p_tf_nuclear_heat_mw",  # Nuclear heating in the TF coil (MW)
        # Used in [CHARGEP]
        "p_plasma_separatrix_mw",  # Charged particle power deposited on divertor (MW)
        "f_p_alpha_plasma_deposited",  # Fraction of alpha power deposited in plasma
        "p_plasma_rad_mw",  # Total radiation Power (MW)
        # Used in [RADIATION]
        "f_ster_div_single"
        "f_a_fw_outboard_hcd"
        # Used in [DIVERTOR]
        "p_div_coolant_pump_mw",  # Divertor coolant pumping power
        "p_div_heat_deposited_mw",  # Total power extracted from divertor (MW)
        # Used in [FIRST_WALL]
        "p_fw_blkt_heat_deposited_mw",  # Power extracted blanket & FW (MW)
        "p_fw_blkt_coolant_pump_mw",  # Pump Power in FW and blanket (MW)
    ]
    (
        p_fusion_total_mw,
        p_hcd_injected_total_mw,
        p_plasma_ohmic_mw,
        p_neutron_total_mw,
        p_non_alpha_charged_mw,
        p_alpha_total_mw,
        p_blkt_multiplication_mw,
        p_blkt_nuclear_heat_total_mw,
        p_div_nuclear_heat_total_mw,
        p_fw_nuclear_heat_total_mw,
        p_shld_nuclear_heat_mw,
        p_tf_nuclear_heat_mw,
        p_plasma_separatrix_mw,
        f_p_alpha_plasma_deposited,
        p_plasma_rad_mw,
        f_ster_div_single,
        f_a_fw_outboard_hcd,
        p_div_coolant_pump_mw,
        p_div_heat_deposited_mw,
        p_fw_blkt_heat_deposited_mw,
        p_fw_blkt_coolant_pump_mw,
    ) = m_file.get_variables(*variables, scan=-1)

    # Used in [PLASMA]
    # Total Power in plasma (MW)
    totalplasma = p_fusion_total_mw + p_hcd_injected_total_mw + p_plasma_ohmic_mw
    # The ohmic and charged particle power (MW)
    pcharohmmw = p_non_alpha_charged_mw + p_plasma_ohmic_mw
    # Alpha particle and HC&D power (MW)
    palpinjmw = p_alpha_total_mw + p_hcd_injected_total_mw

    # Used in [NEUTRONICS]
    # External nuclear heating in blanket (MW)
    pnucemblkt = p_blkt_nuclear_heat_total_mw - p_blkt_multiplication_mw

    # Used in [CHARGEP]
    # Alpha particles hitting first wall (MW)
    p_fw_alpha_mw = p_alpha_total_mw * (1 - f_p_alpha_plasma_deposited)

    # Used in [RADIATION]
    # Radiation deposited on the divertor (MW)
    p_div_rad_total_mw = p_plasma_rad_mw * f_ster_div_single
    # Radiation deposited on HCD (MW)
    p_fw_hcd_rad_total_mw = p_plasma_rad_mw * f_a_fw_outboard_hcd
    # Radiation deposited in the FW (MW)
    p_fw_rad_total_mw = p_plasma_rad_mw - p_div_rad_total_mw - p_fw_hcd_rad_total_mw

    # Used in [FIRST_WALL]
    htpmwblkt = p_fw_blkt_coolant_pump_mw / 2  # Pump power in blanket (MW)
    htpmwfw = p_fw_blkt_coolant_pump_mw / 2  # Pump power in FW (MW)
    p_fw_heat_deposited_mw = (
        p_fw_blkt_heat_deposited_mw - htpmwblkt - p_blkt_nuclear_heat_total_mw
    )  # Power extracted 1st wall (MW)
    # porbitloss = m_file.data['porbitloss'].get_scan(-1) # Charged P. on FW before thermalising
    # p_beam_shine_through_mw = m_file.data['p_beam_shine_through_mw'].get_scan(-1) # Injection shine-through to 1st wall

    # Initialising x and y variables for adjusting 'Plasma Heating' branch tip location
    y_adj_1 = 0
    y_adj_2 = 0

    # Loop 1 to get 'Plasma Heating' branch tip coords; loop 2 to match 'PLASMA' branch
    for _ in range(2):
        # The visual settings of the Sankey Plot
        plt.rcParams.update({"font.size": 9})
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], frameon=False)
        sankey = Sankey(
            ax=ax, unit="MW", margin=0.5, format="%1.0f", scale=1.0 / (totalplasma)
        )

        # --------------------------------------- PLASMA - 0 --------------------------------------

        # Fusion, Injected, Ohmic, -Charged P.-Ohmic, -Alphas-Injected, -Neutrons
        plasma = [
            p_fusion_total_mw,
            p_hcd_injected_total_mw,
            p_plasma_ohmic_mw,
            -pcharohmmw,
            -palpinjmw,
            -p_neutron_total_mw,
        ]
        sankey.add(
            flows=plasma,
            # [left(in), down(in), down(in), up(out), up(out), right(out)]
            orientations=[0, -1, -1, 1, 1, 0],
            trunklength=0.5,
            pathlengths=[0.5, 0.25, 0.25, 0.75, 0.25 + 0.5 * y_adj_1, 0.0],
            # labels=["Fusion","H&CD", "Ohmic", "Charged P.", "Alphas", "Neutrons"])
            labels=[None, None, None, None, None, None],
        )

        # Check to see if the fusion components balance
        if _ == 0 and sqrt(sum(plasma) ** 2) > 0.1:
            print("FUSION power balance =", sum(plasma), "\n")
            exit()

        if _ == 1:
            print(sankey.finish()[0])
        if _ == 1:
            print(sankey.finish()[0].patch)
        if _ == 1:
            print(type(sankey.finish()[0].patch))

        # ------------------------------------- NEUTRONICS - 1 ------------------------------------

        # Neutrons, -Divertor, -1st wall, -Shield, -TF coils, -Blanket+Energy Mult.
        neutrons = [
            p_neutron_total_mw,
            -p_div_nuclear_heat_total_mw,
            -p_fw_nuclear_heat_total_mw,
            -p_shld_nuclear_heat_mw,
            -p_tf_nuclear_heat_mw,
            -pnucemblkt,
        ]
        sankey.add(
            flows=neutrons,
            # left(in), up(out), up(out), up(out), up(out), right(out)
            orientations=[0, 1, 1, 1, 1, 0],
            trunklength=0.5,
            pathlengths=[0.3, 0.25, 0.25, 0.25, 0.25, 0.15],
            prior=0,  # PLASMA
            connect=(5, 0),  # Neutrons
            # labels=["Neutrons", "Divertor", "1st Wall", "Shield", "TF coils", "Blanket"])
            labels=[None, None, None, None, None, None],
        )

        # Checking to see if the neutronics components balance
        if _ == 0 and sqrt(sum(neutrons) ** 2) > 0.1:
            print("NEUTRONS power balance =", sum(neutrons), "\n")
            exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff1_1 = check[0].flows[5] + check[1].flows[0]
            plt.close()
            if diff1_1 > 0.1:
                print("Neutrons [0][5] and [1][0] difference =", diff1_1)
                exit()

        # --------------------------------- CHARGED PARTICLES - 2 ---------------------------------

        # Charge P.+Ohmic, Alpha+Injected, -Divertor, -1st Wall, -Photons
        chargedp = [
            pcharohmmw,
            palpinjmw,
            -p_plasma_separatrix_mw,
            -p_fw_alpha_mw,
            -p_plasma_rad_mw,
        ]
        sankey.add(
            flows=chargedp,
            # down(in), down(in), up(out), up(out), right(out)
            orientations=[-1, -1, 1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.75, 0.25 + 0.5 * y_adj_1, 0.25, 0.25, 0.25],
            prior=0,  # PLASMA
            connect=(3, 0),  # Charged P.+Ohmic
            # labels=["Charged P.", "Alphas", "Divertor", "1st Wall", "Photons"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(chargedp) ** 2) > 0.1:
            print("CHARGEDP power balance =", sum(chargedp))
            exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff2_1 = check[0].flows[3] + check[2].flows[0]
            diff2_2 = check[0].flows[4] + check[2].flows[1]
            plt.close()
            if diff2_1 > 0.1:
                print("Charged P.+Ohmic [0][3] and [2][0] difference =", diff2_1)
                exit()
            if diff2_2 > 0.1:
                print("Alphas+Injected [0][4] and [2][1] difference =", diff2_2)
                exit()

        # ------------------------------------- RADIATION - 3 -------------------------------------

        # Photons, -1st Wall, -Divertor, -H&CD
        radiation = [
            p_plasma_rad_mw,
            -p_fw_rad_total_mw,
            -p_div_rad_total_mw,
            -p_fw_hcd_rad_total_mw,
        ]
        sankey.add(
            flows=radiation,
            # right(in), up(out), up(out), up(out)
            orientations=[
                0,
                -1,
                1,
                1,
            ],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25],
            prior=2,  # CHARGED PARTICLES
            connect=(4, 0),  # Charged P.
            # labels=["Photons", "1st Wall", "Divertor", "H&CD"])
            labels=[None, None, None, None],
        )

        if _ == 0 and sqrt(sum(radiation) ** 2) > 0.1:
            print("RADIATION power balance =", sum(radiation))
            exit()

        if _ == 0:
            check = sankey.finish()
            diff3_1 = check[2].flows[4] + check[3].flows[0]
            plt.close()
            if diff3_1 > 0.1:
                print("Photons [2][4] and [3][0] difference =", diff3_1)
                exit()

        # -------------------------------------- DIVERTOR - 4 -------------------------------------

        # Charged P., Neutrons, Photons, Coolant Pumping, Total Divertor
        divertor = [
            p_plasma_separatrix_mw,
            p_div_nuclear_heat_total_mw,
            p_div_rad_total_mw,
            p_div_coolant_pump_mw,
            -p_div_heat_deposited_mw,
        ]
        sankey.add(
            flows=divertor,
            # down(in), up(in), down(in), up(in), right(out)
            orientations=[-1, -1, -1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25 - 0.5 * y_adj_2, 0.25],
            prior=2,  # CHARGED PARTICLES
            connect=(2, 0),  # Charged P. --> None
            # labels=["Charged P.", "Neutrons", "Photons", "Coolant Pumping", "Divertor Power"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(divertor) ** 2) > 0.1:
            print("DIVERTOR power balance =", sum(divertor))
            exit()

        if _ == 0:
            check = sankey.finish()
            diff4_1 = check[1].flows[1] + check[4].flows[0]
            diff4_2 = check[2].flows[3] + check[4].flows[3]
            plt.close()
            if diff4_1 > 0.1:
                print("Neutrons [1][1] and [4][0] difference =", diff4_1)
                exit()
            if diff4_2 > 0.1:
                print("Charged P. [2][3] and [4][3] difference =", diff4_2)
                exit()

        # ---------------------------------------- 1ST WALL - 5 ---------------------------------------

        # Alphas, Neutrons, Photons, Coolant Pumping, Total 1st Wall
        first_wall = [
            p_fw_alpha_mw,
            p_fw_nuclear_heat_total_mw,
            p_fw_rad_total_mw,
            htpmwfw,
            -p_fw_heat_deposited_mw,
        ]
        sankey.add(
            flows=first_wall,
            orientations=[0, -1, 1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25, 0.25],
            prior=1,
            connect=(2, 1),
            # labels=["Alphas", "Neutrons", "Radiation", "Coolant Pumping", "FW Power"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(first_wall) ** 2) > 0.1:
            print("FIRST_WALL power balance =", sum(first_wall))
            exit()
        """# -------------------------------------- BLANKET - 6 --------------------------------------

        # Blanket - Energy mult., Energy Mult., pumping power, Blanket
        BLANKET = [pnucemblkt, p_blkt_multiplication_mw, htpmwblkt, -p_blkt_heat_deposited_mw]
        sankey.add(flows=BLANKET,
                   # left(in), down(in), down(in), right(out)
                   orientations=[0, -1, -1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25, 0.25, 0.25],
                   #prior=1, # NEUTRONICS
                   #connect=(1, 0), # Blanket --> None
                   labels=[None, "Energy Mult.", "Coolant Pumping", "Blanket"])

        # Checking to see if the blanket components balance
        if _ == 0:
            if sqrt(sum(BLANKET)**2) > 0.1:
                print("BLANKET power balance =", sum(BLANKET), "\n")
                exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff = check[1].flows[1]+check[3].flows[0]
            if diff > 0.1:
                print("The difference between [1][1] and [3][0] =", diff)
                exit()"""
        """# --------------------------------------- SHIELD - 7 --------------------------------------

        # Neutrons, Coolant pumping, Total power
        SHIELD = [p_shld_nuclear_heat_mw, p_shld_coolant_pump_mw, -p_shld_heat_deposited_mw]
        sankey.add(flows=SHIELD,
                   orientations=[-1, -1, 1],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["Neutrons", "Coolant Pumping", "Shield Power"])

        if _ == 0:
            if sqrt(sum(SHIELD)**2) > 0.1:
                print("SHIELD power balance =", sum(SHIELD))
                exit()"""
        """# ------------------------------------ PRIMARY HEAT - 7 -----------------------------------

        # 1st wall, Blanket, Shield, Divertor, Total thermal power
        HEAT = [p_fw_heat_deposited_mw, p_blkt_heat_deposited_mw, p_shld_heat_deposited_mw, p_div_heat_deposited_mw, -p_plant_primary_heat_mw]
        sankey.add(flows=HEAT,
                   orientations=[1, 0, -1, 1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25, 0.25, 0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["1st Wall", "Blanket", "Shield", "Divertor", "Total Power"])

        if _ == 0:
            if sqrt(sum(HEAT)**2) > 0.1:
                print("PRIMARY power balance =", sum(HEAT))
                exit()"""
        """# ------------------------------- ELECTRICITY CONVERSION - 8 ------------------------------

        # Total thermal, Elctricty conversion loss, Gross Electricity
        GROSS = [p_plant_primary_heat_mw, -pelectloss, -p_plant_electric_gross_mw]
        sankey.add(flows=GROSS,
                   orientations=[0, -1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["Thermal Power", "Conversion loss", "Gross Electricity"])

        if _ == 0:
            if sqrt(sum(GROSS)**2) > 0.1:
                print("GROSS power balance =", sum(GROSS))
                exit()"""

        # ------------------------------ RECIRCULATED ELECTRICITY - 9 -----------------------------
        """# ---------------------------------------- HCD - 11 ----------------------------------------

        # HCD loss + injected, -injected, -HCD loss
        HCD = [p_hcd_electric_loss_mw+p_hcd_injected_total_mw, -p_hcd_injected_total_mw, -p_hcd_electric_loss_mw]
        assert(sum(HCD)**2 < 0.5)
        sankey.add(flows=HCD,
                   # [down(in), up(out), down(out)]
                   orientations=[-1, 1, -1],
                   #prior=0, # PLASMA
                   #connect=(1, 1), # H&CD --> None
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25, 0.25],
                   labels=['H&CD power', None, 'H&CD loss'])"""

        fig.tight_layout()

        if _ == 0:
            plt.close()

        # Matching PLASMA and CHARGED PARTICLES 'Alphas' branches
        # x_adj_1, y_adj_1 = diagrams[2].tips[1] - diagrams[0].tips[4]
        # Matching CHARGED PARTICLES and DIVERTOR 'Charged P.' branches
        # x_adj_2, y_adj_2 = diagrams[4].tips[3] - diagrams[2].tips[3]
        # x_adj_3, y_adj_3 = diagrams[3].tips[3] - diagrams[4].tips[0]

    # --------------------------------------- Label Positioning ---------------------------------------

    # Munipulating the positioning of the branch labels
    # -ve to left and down; +ve to right and up
    # pos[0] = x-axis; pos[1] = y-axis
    """for d in diagrams:
        y = 0
        for t in d.texts:
            pos = tuple(np.ndarray.tolist(d.tips[y]))
            t.set_position(pos)
            if t == diagrams[0].texts[0]: # Fusion Power
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[0].texts[1]: # H&CD
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_hcd_injected_total_mw/totalplasma)-0.05,pos[1]))
            if t == diagrams[0].texts[2]: # Ohmic
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*(p_plasma_ohmic_mw/totalplasma)+0.05,pos[1]))
            if t == diagrams[0].texts[3]: # Neutrons
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[0].texts[4]: # Charged Particles
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_non_alpha_charged_mw/totalplasma)-0.05,pos[1]))
            if t == diagrams[0].texts[5]: # Alphas
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*(p_alpha_total_mw/totalplasma)+0.05,pos[1]-0.1))
            if t == diagrams[1].texts[0]: # H&CD power
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*((p_hcd_electric_loss_mw+p_hcd_injected_total_mw)/totalplasma)-0.05,pos[1]))
            if t == diagrams[1].texts[2]: # H&CD losses
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+(p_hcd_electric_loss_mw/totalplasma)+0.05,pos[1]))
            if t == diagrams[2].texts[1]: # Energy Multiplication
                t.set_horizontalalignment('center')
                t.set_position((pos[0],pos[1]-0.2))
            if t == diagrams[2].texts[2]: # Blanket
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[2].texts[3]: # Divertor
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_div_nuclear_heat_total_mw/totalplasma)-0.1,pos[1]))
            if t == diagrams[3].texts[2]: # Rad.FW
                t.set_horizontalalignment('right')
                t.set_position((pos[0],pos[1]+0.5*(p_fw_rad_total_mw/totalplasma)+0.15))
            if t == diagrams[3].texts[3]: # Charged P.
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*((p_plasma_separatrix_mw+p_fw_alpha_mw)/totalplasma)+0.1,pos[1]+0.05))
            if t == diagrams[3].texts[4]: # Rad. Div.
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_div_rad_total_mw/totalplasma)-0.1,pos[1]))
            y += 1"""


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
        p_plant_electric_recirc_mw,
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

    # Initialising x and y variables for adjusting 'Plasma Heating' branch tip location
    x_adj, y_adj = 0, 0

    # Loop 1 to get 'Plasma Heating' branch tip coords; loop 2 to match 'PLASMA' branch
    for _ in range(2):
        # ------------------------------------ Visual Settings ------------------------------------

        plt.rcParams.update({"font.size": 9})  # Setting font size to 9
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], frameon=False)
        sankey = Sankey(
            ax=ax, unit="MW", margin=0.0, format="%1.0f", scale=1.0 / (totalplasma)
        )

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
                0.5,
                0.8 + 0.5 * y_adj,
                -0.1 + 0.5 * x_adj,
            ],  # 'Plasma Heating' adjust
            labels=["Fusion Power", None, "Plasma"],
        )

        # --------------------------------- ENERGY DEPOSITION - 1 ---------------------------------

        # Plasma power, - divertor deposited power, - blanket deposited power
        deposition = [totalplasma, -totalblktetc - totaldivetc - totalcpetc]
        # Check if difference >2 between plasma and divertor + blanket
        if _ == 1 and sqrt(sum(deposition) ** 2) > 2:
            print(
                "\ncomponents power balance difference =",
                totalplasma - totaldivetc - totalblktetc - totalcpetc,
            )
        sankey.add(
            flows=deposition,
            orientations=[0, 0],  # [right(in), up(in), right(out)]
            prior=0,  # PLASMA
            connect=(2, 0),  # Plasma --> None
            pathlengths=[0.2, 0.2 + 0.5 * x_adj],  # 'Plasma Heating' adjust
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
        if _ == 1 and sqrt(sum(blanketsetc) ** 2) > 2:
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
            pathlengths=[0.5, 0.25, 0.0],
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
            pathlengths=[0.2, 0.7, 0.4],
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
                pathlengths=[0.1, 0.25, 1.5],
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
                pathlengths=[0.25, 1.0, 0.5],
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
            pathlengths=[0.1, 0.25, 0.8],
            labels=[None, "Core Systems", "Heating System"],
        )

        # --------------------------------------- LOSSES - 6 --------------------------------------

        # HCD: Heating system, -Plasma heating, -losses
        hcd = [
            p_hcd_electric_total_mw - p_cp_coolant_pump_elec_mw,
            -p_hcd_injected_total_mw,
            -p_hcd_electric_total_mw
            + p_hcd_injected_total_mw
            + p_cp_coolant_pump_elec_mw,
        ]
        sankey.add(
            flows=hcd,
            orientations=[0, -1, 0],  # [left(in), up(out), left(out)]
            prior=5,  # RECIRC
            connect=(2, 0),  # Heating System --> None
            pathlengths=[0.5, 0.8 + 0.5 * y_adj, 0.4],  # 'Plasma Heating' adjust
            labels=[None, "Plasma Heating", "Losses"],
        )

        # Collecting Sankey diagram and applying a condensed layout
        diagrams = sankey.finish()
        fig.tight_layout()

        # Difference in branch tip locations for 'Plasma Heating'
        x_adj, y_adj = diagrams[0].tips[1] - diagrams[6].tips[1]

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
            if t == diagrams[6].texts[1]:  # Plasma Heating
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] + 0.5 * (p_hcd_injected_total_mw / totalplasma) + 0.1,
                    pos[1] - 0.05,
                ))
            if t == diagrams[6].texts[2]:  # Losses
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
