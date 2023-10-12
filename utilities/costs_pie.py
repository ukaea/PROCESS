#!/usr/bin/env python
"""
Code to display the cost breakdown as a pie chart

Stuart Muldrew (stuart.muldrew@ukaea.uk)
06/09/2018

History
04/04/2019 SIM Added step_cost_model

"""

# Imported libraries
import argparse
import process.io.mfile as mf
import matplotlib.pyplot as plt


def orig_cost_model():
    """

    Plot pie chart for the orginal 1990 cost model.
    Two plots produced: (1) Breakdown of the direct costs and (2) Direct, indirect, etc.

    """
    # Read Cost Values
    c21 = m_file.data["c21"].get_scan(-1)  # Site and Buildings
    c221 = m_file.data["c221"].get_scan(-1)  # Reactor Systems
    c222 = m_file.data["c222"].get_scan(-1)  # Magnets
    c223 = m_file.data["c223"].get_scan(-1)  # Power Injection
    c224 = m_file.data["c224"].get_scan(-1)  # Vacuum Systems
    c225 = m_file.data["c225"].get_scan(-1)  # Power Conditioning
    c226 = m_file.data["c226"].get_scan(-1)  # Heat Transport System
    c227 = m_file.data["c227"].get_scan(-1)  # Fuel Handling System
    c228 = m_file.data["c228"].get_scan(-1)  # Instrumentation and Control
    c229 = m_file.data["c229"].get_scan(-1)  # Maintenance Equipment
    c23 = m_file.data["c23"].get_scan(-1)  # Turbine Plant Equipment
    c24 = m_file.data["c24"].get_scan(-1)  # Electric Plant Equipment
    c25 = m_file.data["c25"].get_scan(-1)  # Miscellaneous Plant Equipment
    c26 = m_file.data["c26"].get_scan(-1)  # Heat Rejection System

    cdirt = m_file.data["cdirt"].get_scan(-1)  # Plant Direct Cost
    c9 = m_file.data["c9"].get_scan(-1)  # Indirect Cost
    ccont = m_file.data["ccont"].get_scan(-1)  # Total Contingency

    # Interest during construction is linked to ireactor = 1
    if "moneyint" in m_file.data.keys():
        moneyint = m_file.data["moneyint"].get_scan(-1)  # Interest during Construction
        labels2 = [
            "Plant Direct Cost",
            "Indirect Cost",
            "Total Contingency",
            "Interest during Construction",
        ]
        sizes2 = [cdirt, c9, ccont, moneyint]
    else:
        labels2 = ["Plant Direct Cost", "Indirect Cost", "Total Contingency"]
        sizes2 = [cdirt, c9, ccont]

    # No turbines if ireactor = 0
    if c23 > 1.0e-3:
        labels = [
            "Magnets and Power Conditioning",
            "Site and Buildings",
            "Maintenance Equipment",
            "Power Injection",
            "Reactor Systems",
            "Fuel Handling System",
            "Instrumentation and Control",
            "Turbine Plant Equipment",
            "Heat Transport System",
            "Other",
        ]

        sizes = [
            c222 + c225,
            c21,
            c229,
            c223,
            c221,
            c227,
            c228,
            c23,
            c226,
            c224 + c24 + c25 + c26,
        ]

    else:
        labels = [
            "Magnets and Power Conditioning",
            "Site and Buildings",
            "Maintenance Equipment",
            "Power Injection",
            "Reactor Systems",
            "Fuel Handling System",
            "Instrumentation and Control",
            "Heat Transport System",
            "Other",
        ]

        sizes = [
            c222 + c225,
            c21,
            c229,
            c223,
            c221,
            c227,
            c228,
            c226,
            c224 + c24 + c25 + c26,
        ]

    # Setup figures
    # Plot direct cost items
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.pie(sizes, labels=labels, autopct="%1.1f%%")
    ax1.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Plot overall breakdown
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.pie(sizes2, labels=labels2, autopct="%1.1f%%")
    ax2.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Save figures if option selected
    if args.save:
        fig1.savefig("direct_cost_pie.pdf")
        fig2.savefig("cost_pie.pdf")
    plt.show()


def new_cost_model():
    """

    Plot pie chart for the new 2014 cost model.

    """
    # Read Cost Values
    s09 = m_file.data["s09"].get_scan(-1)  # Buildings
    s13 = m_file.data["s13"].get_scan(-1)  # Land
    s21 = m_file.data["s21"].get_scan(-1)  # TF Coils
    s27 = m_file.data["s27"].get_scan(-1)  # First wall and blanket
    s31 = m_file.data["s31"].get_scan(-1)  # Active maintenance and remote handling
    s34 = m_file.data["s34"].get_scan(-1)  # Vacuum vessel and liquid nitrogen plant
    s35 = m_file.data["s35"].get_scan(-1)  # System for converting heat to electricity
    s36 = m_file.data["s36"].get_scan(-1)  # CS and PF coils
    s51 = m_file.data["s51"].get_scan(-1)  # Cryoplant and distribution
    s52 = m_file.data["s52"].get_scan(-1)  # Electrical power supply and distribution
    s59 = m_file.data["s59"].get_scan(-1)  # Additional project expenditure
    s61 = m_file.data["s61"].get_scan(-1)  # Remaining subsystems

    labels = [
        "Land and Buildings",
        "TF Coils",
        "First wall and blanket",
        "Active maintenance and remote handling",
        "Vacuum vessel and liquid nitrogen plant",
        "CS and PF coils",
        "Cryoplant and distribution",
        "Electrical power supply and distribution",
        "Additional project expenditure",
        "Other subsystems",
    ]

    # Split up Remaining Subsystems as it is too large
    sizes = [
        s09 + s13,
        s21,
        s27,
        s31,
        s34,
        s36,
        s51,
        s52,
        s59,
        s35 + s61 - s36 - s51 - s52 - s59,
    ]

    # Setup figure
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.pie(sizes, labels=labels, autopct="%1.1f%%")
    ax1.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Save figures if option selected
    if args.save:
        fig1.savefig("cost_pie.pdf")
    plt.show()


def step_cost_model():
    """

    Plot pie chart for the STEP cost model.
    Two plots produced: (1) Breakdown of the direct costs and (2) Direct, indirect, etc.

    """
    # Read Cost Values
    step20 = m_file.data["step20"].get_scan(-1)  # Land and Rights
    step21 = m_file.data["step21"].get_scan(
        -1
    )  # Building and Site Service Infrastructure
    step2201 = m_file.data["step2201"].get_scan(-1)  # Reactor Equipment
    step220101 = m_file.data["step220101"].get_scan(-1)  # Blanket and First Wall
    step220102 = m_file.data["step220102"].get_scan(-1)  # Shield
    step22010301 = m_file.data["step22010301"].get_scan(-1)  # TF Coils
    step22010302 = m_file.data["step22010302"].get_scan(-1)  # PF Coils
    step22010303 = m_file.data["step22010303"].get_scan(-1)  # Central Solenoid
    step22010304 = m_file.data["step22010304"].get_scan(-1)  # Control Coils
    step220103 = step22010301 + step22010302 + step22010303 + step22010304  # Magnets
    step2202 = m_file.data["step2202"].get_scan(-1)  # Heat Transfer System
    step22 = m_file.data["step22"].get_scan(-1)  # Reactor Plant Equipment
    step23 = m_file.data["step23"].get_scan(-1)  # Turbine Plant Equipment
    step24 = m_file.data["step24"].get_scan(-1)  # Electric Plant Equipment
    step25 = m_file.data["step25"].get_scan(-1)  # Miscellaneous Plant Equipment

    cdirt = m_file.data["cdirt"].get_scan(-1)  # Plant Direct Cost
    step91 = m_file.data["step91"].get_scan(
        -1
    )  # Construction Facilities, Equipment and Services
    step92 = m_file.data["step92"].get_scan(
        -1
    )  # Engineering and Costruction Management Services
    step93 = m_file.data["step93"].get_scan(-1)  # Other Costs

    # Interest during construction is linked to ireactor = 1
    if "moneyint" in m_file.data.keys():
        moneyint = m_file.data["moneyint"].get_scan(-1)  # Interest during Construction
        labels2 = [
            "Plant Direct Cost",
            "Construction Facilities,\n Equipment and\n Services",
            "Engineering and\n Costruction Management\n Services",
            "Other Costs",
            "Interest during Construction",
        ]
        sizes2 = [cdirt, step91, step92, step93, moneyint]
    else:
        labels2 = [
            "Plant Direct Cost",
            "Construction Facilities,\n Equipment and\n Services",
            "Engineering and\n Costruction Management\n Services",
            "Other Costs",
        ]
        sizes2 = [cdirt, step91, step92, step93]

    # Direct Cost Breakdown
    labels = [
        "Site and Buildings",
        "Blanket and First Wall",
        "Shield",
        "Magnets",
        "Other Reactor Equipment",
        "Heat Transfer System",
        "Other Reactor Plant Equipment",
        "Turbine Plant Equipment",
        "Electric Plant Equipment",
        "Miscellaneous Plant Equipment",
    ]

    sizes = [
        step20 + step21,
        step220101,
        step220102,
        step220103,
        step2201 - step220101 - step220102 - step220103,
        step2202,
        step22 - step2201 - step2202,
        step23,
        step24,
        step25,
    ]

    # Setup figures
    # Plot direct cost items
    fig1, ax1 = plt.subplots(figsize=(9, 5))
    ax1.pie(sizes, labels=labels, autopct="%1.1f%%")
    ax1.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Plot overall breakdown
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.pie(sizes2, labels=labels2, autopct="%1.1f%%")
    ax2.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Save figures if option selected
    if args.save:
        fig1.savefig("direct_cost_pie.pdf")
        fig2.savefig("cost_pie.pdf")
    plt.show()


# Main code
if __name__ == "__main__":
    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Displays the cost breakdown as a pie chart.  "
        "For more information contact Stuart.Muldrew@ukaea.uk"
    )

    parser.add_argument(
        "-f",
        metavar="MFILE",
        type=str,
        default="MFILE.DAT",
        help="specify the MFILE (default=MFILE.DAT)",
    )

    parser.add_argument("-s", "--save", help="save figure", action="store_true")

    args = parser.parse_args()

    m_file = mf.MFile(args.f)

    # Check which cost model is being used
    if "c21" in m_file.data.keys():
        orig_cost_model()
    elif "s01" in m_file.data.keys():
        new_cost_model()
    elif "step20" in m_file.data.keys():
        step_cost_model()
    else:
        print("ERROR: Cannot identify cost data, check MFILE!")
