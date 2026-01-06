import argparse
import pathlib
import re

try:
    import plotly.graph_objects as go

    PLOT_SANKEY = True
except ImportError:
    PLOT_SANKEY = False

from process.io.mfile import MFile


def main(args=None):
    ###########################################################
    # Usage

    if not PLOT_SANKEY:
        print(
            "\nPlotly is not installed, unable to create sankey diagram!\n"
            "Install plotly by installing the optional 'plotly' dependency "
            "e.g. \"pip install -e '.[plotly]'\""
        )
        return

    parser = argparse.ArgumentParser(
        description="Program to plot\
     the power flow in PROCESS using a Sankey diagram."
    )

    parser.add_argument("-e", "--end", default="pdf", help="file format, default = pdf")

    parser.add_argument(
        "-m", "--mfile", default="MFILE.DAT", help="mfile name, default = MFILE.DAT"
    )

    args = parser.parse_args(args)

    #########################################################
    # main program

    plot_power_balance_sankey(args.mfile)


def plot_power_balance_sankey(m_file):
    m_file = MFile(m_file)
    p_hcd_injected_total_mw = m_file.data["p_hcd_injected_total_mw"].get_scan(-1)
    p_plasma_ohmic_mw = m_file.data["p_plasma_ohmic_mw"].get_scan(-1)
    p_alpha_total_mw = m_file.data["p_alpha_total_mw"].get_scan(-1)
    p_neutron_total_mw = m_file.data["p_neutron_total_mw"].get_scan(-1)
    p_plasma_rad_mw = m_file.data["p_plasma_rad_mw"].get_scan(-1)
    p_fw_rad_total_mw = m_file.data["p_fw_rad_total_mw"].get_scan(-1)
    p_fw_alpha_mw = p_alpha_total_mw * (
        1 - m_file.data["f_p_alpha_plasma_deposited"].get_scan(-1)
    )
    p_blkt_nuclear_heat_total_mw = m_file.data["p_blkt_nuclear_heat_total_mw"].get_scan(
        -1
    )

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
        - m_file.data["p_blkt_multiplication_mw"].get_scan(-1),  # 6
        p_fw_rad_total_mw,  # 7
        m_file.data["p_div_nuclear_heat_total_mw"].get_scan(-1),  # 8
        m_file.data["p_div_rad_total_mw"].get_scan(-1),  # 9
        m_file.data["p_fw_heat_deposited_mw"].get_scan(-1),  # 10
        m_file.data["p_blkt_heat_deposited_mw"].get_scan(-1),  # 11
        m_file.data["p_div_heat_deposited_mw"].get_scan(-1),  # 12
        m_file.data["p_fw_blkt_heat_deposited_mw"].get_scan(-1),  # 13
        m_file.data["p_plant_primary_heat_mw"].get_scan(-1),  # 14
        m_file.data["p_plant_electric_gross_mw"].get_scan(-1),  # 15
        m_file.data["p_plant_electric_net_mw"].get_scan(-1),  # 16
        m_file.data["p_hcd_electric_total_mw"].get_scan(-1),  # 17
        m_file.data["p_hcd_electric_loss_mw"].get_scan(-1),  # 18
        p_hcd_injected_total_mw,  # 19
        m_file.data["p_plant_core_systems_elec_mw"].get_scan(-1),  # 20
        m_file.data["p_cryo_plant_electric_mw"].get_scan(-1),  # 21
        m_file.data["p_plant_electric_base_total_mw"].get_scan(-1),  # 22
        m_file.data["p_tf_electric_supplies_mw"].get_scan(-1),  # 23
        m_file.data["p_pf_electric_supplies_mw"].get_scan(-1),  # 24
        m_file.data["vachtmw"].get_scan(-1),  # 25
        m_file.data["p_tritium_plant_electric_mw"].get_scan(-1),  # 26
        m_file.data["p_coolant_pump_elec_total_mw"].get_scan(-1),  # 27
        m_file.data["p_coolant_pump_loss_total_mw"].get_scan(-1),  # 28
        m_file.data["p_div_coolant_pump_mw"].get_scan(-1),  # 29
        m_file.data["p_fw_blkt_coolant_pump_mw"].get_scan(-1),  # 30
        m_file.data["p_fw_blkt_coolant_pump_mw"].get_scan(-1),  # 31
        m_file.data["p_div_coolant_pump_mw"].get_scan(-1),  # 32
        m_file.data["p_shld_coolant_pump_mw"].get_scan(-1),  # 33
        m_file.data["p_shld_coolant_pump_mw"].get_scan(-1),  # 34
        m_file.data["p_shld_heat_deposited_mw"].get_scan(-1),  # 35
        m_file.data["p_shld_nuclear_heat_mw"].get_scan(-1),  # 36
        m_file.data["p_cryo_plant_electric_mw"].get_scan(-1),  # 37
        m_file.data["p_plant_electric_base_total_mw"].get_scan(-1),  # 38
        m_file.data["p_tf_electric_supplies_mw"].get_scan(-1),  # 39
        m_file.data["p_pf_electric_supplies_mw"].get_scan(-1),  # 40
        m_file.data["vachtmw"].get_scan(-1),  # 41
        m_file.data["p_tritium_plant_electric_mw"].get_scan(-1),  # 42
        m_file.data["p_tf_nuclear_heat_mw"].get_scan(-1),  # 43
        m_file.data["p_tf_nuclear_heat_mw"].get_scan(-1),  # 44
        m_file.data["p_hcd_electric_loss_mw"].get_scan(-1),  # 45
        m_file.data["p_coolant_pump_loss_total_mw"].get_scan(-1),  # 46
        m_file.data["p_fw_coolant_pump_mw"].get_scan(
            -1
        ),  # 47  Should only show if FW and Bkt pumps are seperate
        m_file.data["p_blkt_coolant_pump_mw"].get_scan(
            -1
        ),  # 48  Should only show if FW and Blkt pumps are seperate
        m_file.data["p_beam_shine_through_mw"].get_scan(
            -1
        ),  # 49 Should show in beams are present
        m_file.data["p_beam_orbit_loss_mw"].get_scan(
            -1
        ),  # 50 Should show in beams are present
        m_file.data["p_cp_shield_nuclear_heat_mw"].get_scan(
            -1
        ),  # 51 Neutrons to CP shield, should only show if CP shield is present
        m_file.data["p_plasma_separatrix_mw"].get_scan(
            -1
        ),  # 52 Plasma separatrix power to divertor
        m_file.data["p_div_secondary_heat_mw"].get_scan(
            -1
        ),  # 53 Divertor secondary heat,
        m_file.data["p_shld_secondary_heat_mw"].get_scan(-1),  # 54 Shield secondary heat
        m_file.data["p_fw_hcd_nuclear_heat_mw"].get_scan(
            -1
        ),  # 55 Neutron power to H&CD & Diagnostics
        m_file.data["p_fw_hcd_rad_total_mw"].get_scan(
            -1
        ),  # 56: Radiation to H&CD & Diagnostics
        m_file.data["p_plant_secondary_heat_mw"].get_scan(
            -1
        ),  # 57: Total Secondary Heat
        m_file.data["p_hcd_secondary_heat_mw"].get_scan(
            -1
        ),  # 58: H&CD & Diagnostics secondary heat
        m_file.data["p_turbine_loss_mw"].get_scan(-1),  # 59: Turbine Loss
        m_file.data["p_fw_nuclear_heat_total_mw"].get_scan(-1),  # 60: FW nuclear heat
        p_alpha_total_mw
        * m_file.data["f_p_alpha_plasma_deposited"].get_scan(
            -1
        ),  # 61: Alpha particles back to plasma
        m_file.data["p_blkt_multiplication_mw"].get_scan(-1),
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

    sankey_dict = {
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
    fig = go.Figure(data=[sankey_dict])

    fig.update_layout({
        "title_text": "Fusion Power Balance Sankey Diagram",
        "font_size": 7,
        "autosize": True,
        "margin": {"l": 40, "r": 40, "t": 40, "b": 40},
    })
    # Strip 'MFILE' from the filename for the HTML output
    # Remove the character before "MFILE" and "MFILE" itself from the filename
    html_output_path = pathlib.Path(
        re.sub(r"(.)?[ \.\_]?MFILE", r"\1_plotly_sankey", m_file.filename)
    ).with_suffix(".html")
    fig.write_html(str(html_output_path))
    print(f"Interactive Sankey diagram saved to {html_output_path}")
    return fig


if __name__ == "__main__":
    main()
