import argparse
import pathlib

import matplotlib as mpl
import plotly.graph_objects as go

from process.io.mfile import MFile

mpl.use("Agg")


def main(args=None):
    ###########################################################
    # Usage

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

    fig = plot_power_balance_sankey(args.mfile)

    # Get directory of mfile
    mfile_path = pathlib.Path(args.mfile).resolve()
    mfile_dir = mfile_path.parent
    output_path = mfile_dir / f"SankeyPowerFlow.{args.end}"

    # Save the Plotly figure as a static image
    fig.write_image(str(output_path))
    print(f"Sankey diagram saved to {output_path}")


def plot_power_balance_sankey(m_file):
    m_file = MFile(m_file)
    # Example: Extract values from m_file as in your matplotlib code
    # Replace these with actual values from your MFile object
    p_fusion_total_mw = m_file.data["p_fusion_total_mw"].get_scan(-1)
    p_hcd_injected_total_mw = m_file.data["p_hcd_injected_total_mw"].get_scan(-1)
    p_plasma_ohmic_mw = m_file.data["p_plasma_ohmic_mw"].get_scan(-1)
    p_alpha_total_mw = m_file.data["p_alpha_total_mw"].get_scan(-1)
    p_neutron_total_mw = m_file.data["p_neutron_total_mw"].get_scan(-1)
    p_plasma_rad_mw = m_file.data["p_plasma_rad_mw"].get_scan(-1)
    p_fw_nuclear_heat_total_mw = m_file.data["p_fw_nuclear_heat_total_mw"].get_scan(-1)
    p_fw_rad_total_mw = m_file.data["p_fw_rad_total_mw"].get_scan(-1)
    p_fw_alpha_mw = p_alpha_total_mw * (1 - m_file.data["f_alpha_plasma"].get_scan(-1))
    p_fw_heat_deposited_mw = m_file.data["p_fw_blkt_heat_deposited_mw"].get_scan(-1)
    p_blkt_nuclear_heat_total_mw = m_file.data["p_blkt_nuclear_heat_total_mw"].get_scan(
        -1
    )
    p_blkt_multiplication_mw = m_file.data["p_blkt_multiplication_mw"].get_scan(-1)
    p_plant_electric_gross_mw = m_file.data["p_plant_electric_gross_mw"].get_scan(-1)
    p_plant_electric_net_mw = m_file.data["p_plant_electric_net_mw"].get_scan(-1)

    # Define node labels (linearized flow)
    labels = [
        "0: H&CD injector",
        "1: Ohmic",
        "2: Plasma Fusion Power",
        "3: Alpha particles",
        "4: Neutrons",
        "5: Radiation",
        "6: First Wall",
        "7: Blanket",
        "8: Divertor",
        "9: FW+Blkt",
        "10: Primary Thermal",
        "11: Turbine",
        "12: Gross Electric",
        "13: Net Electric",
        "14: HCD Electric Power",
        "15: HCD electric losses",
        "16: Core systems",
        "17: Cryo plant",
        "18: Base plant load",
        "19: TF coils",
        "20: PF coils",
        "21: Vacuum pumps",
        "22: Tritium plant",
        "23: Coolant pumps electric",
        "24: Coolant pump electric losses",
        "25: Divertor pump",
        "26: FW+Blkt pumps",
        "27: Shield pump",
        "28: Shield",
        "29: Secondary heat",
        "30: TF nuclear heat",
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
    ]
    values = [
        p_hcd_injected_total_mw,  # 0
        p_plasma_ohmic_mw,  # 1
        p_alpha_total_mw,  # 2
        p_neutron_total_mw,  # 3
        p_plasma_rad_mw,  # 4
        p_fw_alpha_mw,  # 5
        p_blkt_nuclear_heat_total_mw,  # 6
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
    ]
    # Add value labels to the links
    value_labels = [f"{v:.3f} MW" for v in values]

    sankey_dict = {
        "type": "sankey",
        "node": {
            "pad": 15,
            "thickness": 20,
            "line": {"color": "black", "width": 0.5},
            "label": labels,
        },
        "link": {
            "source": sources,
            "target": targets,
            "value": values,
            "label": value_labels,
        },
    }
    fig = go.Figure(data=[sankey_dict])

    fig.update_layout(title_text="Fusion Power Balance Sankey Diagram", font_size=7)
    # Save as interactive HTML file
    fig.update_layout(
        title_text="Fusion Power Balance Sankey Diagram",
        font_size=7,
        autosize=True,
    )
    html_output_path = pathlib.Path(m_file.filename).with_suffix(".html")
    fig.write_html(str(html_output_path))
    print(f"Interactive Sankey diagram saved to {html_output_path}")
    return fig


if __name__ == "__main__":
    main()
