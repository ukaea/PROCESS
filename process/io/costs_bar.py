#!/usr/bin/env python
"""

Code to produce costs bar chart.
Can take multiple input files.

Stuart Muldrew (stuart.muldrew@ukaea.uk)
11/09/2018

History
04/04/2019 SIM Added comp_step

"""

# Imported libraries
import argparse
import process.io.mfile as mf
import matplotlib.pyplot as plt
import numpy as np
import sys
from typing import List


def comp_orig(args, mfile_list: List[str], inflate: float) -> None:
    """

    Plot bar chart for the orginal 1990 cost model.
    Two plots produced: (1) Breakdown of the direct costs and (2) Direct, indirect, etc.

    """
    # Setup figures
    labels = [
        "Magnets and\n Power Conditioning",
        "Site and Buildings",
        "Maintenance\n Equipment",
        "Power Injection",
        "Reactor Systems",
        "Fuel Handling\n System",
        "Instrumentation\n and Control",
        "Turbine Plant\n Equipment",
        "Heat Transport\n System",
        "Other",
    ]
    labels2 = [
        "Plant Direct\n Cost",
        "Indirect\n Cost",
        "Total\n Contingency",
        "Interest during\n Construction",
    ]
    index = np.arange(len(labels))
    index2 = np.arange(len(labels2))
    bar_width = 0.7 / len(mfile_list)
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()

    # Read cost data
    for id, item in enumerate(mfile_list):
        cost = np.zeros(18)
        cost[0] = item.data["c21"].get_scan(-1)  # Site and Buildings
        cost[1] = item.data["c221"].get_scan(-1)  # Reactor Systems
        cost[2] = item.data["c222"].get_scan(-1)  # Magnets
        cost[3] = item.data["c223"].get_scan(-1)  # Power Injection
        cost[4] = item.data["c224"].get_scan(-1)  # Vacuum Systems
        cost[5] = item.data["c225"].get_scan(-1)  # Power Conditioning
        cost[6] = item.data["c226"].get_scan(-1)  # Heat Transport System
        cost[7] = item.data["c227"].get_scan(-1)  # Fuel Handling System
        cost[8] = item.data["c228"].get_scan(-1)  # Instrumentation and Control
        cost[9] = item.data["c229"].get_scan(-1)  # Maintenance Equipment
        cost[10] = item.data["c23"].get_scan(-1)  # Turbine Plant Equipment
        cost[11] = item.data["c24"].get_scan(-1)  # Electric Plant Equipment
        cost[12] = item.data["c25"].get_scan(-1)  # Miscellaneous Plant Equipment
        cost[13] = item.data["c26"].get_scan(-1)  # Heat Rejection System
        cost[14] = item.data["cdirt"].get_scan(-1)  # Plant Direct Cost
        cost[15] = item.data["c9"].get_scan(-1)  # Indirect Cost
        cost[16] = item.data["ccont"].get_scan(-1)  # Total Contingency
        cost[17] = item.data["moneyint"].get_scan(-1)  # Interest during Construction

        # Explain why moneyint is missing
        if "moneyint" not in item.data.keys():
            print(
                "Interest during construction (moneyint) is only calculated for ireactor = 1"
            )

        # Inflate costs using value parsed if specified
        if args.inf:
            cost = inflate * cost

        # Simplify grouping
        sizes = [
            cost[2] + cost[5],
            cost[0],
            cost[9],
            cost[3],
            cost[1],
            cost[7],
            cost[8],
            cost[10],
            cost[6],
            cost[4] + cost[11] + cost[12] + cost[13],
        ]

        # Direct, indirect costs etc. for second plot
        sizes2 = [cost[14], cost[15], cost[16], cost[17]]

        # Plot bar charts
        ax.bar(index + id * bar_width, sizes, bar_width, label=args.f[id])
        ax2.bar(index2 + id * bar_width, sizes2, bar_width, label=args.f[id])

    # Plot labels
    ax.set_xticks(index + (len(mfile_list) - 1) * 0.5 * bar_width)
    ax2.set_xticks(index2 + (len(mfile_list) - 1) * 0.5 * bar_width)
    ax.set_xticklabels(labels, rotation=90)
    ax2.set_xticklabels(labels2, rotation=90)
    ax.legend()
    ax2.legend()

    # Adjust axis label depending on if inflation factor is used
    if args.inf:
        ax.set_ylabel("%.2f x (1990 M$)" % inflate)
        ax2.set_ylabel("%.2f x (1990 M$)" % inflate)
    else:
        ax.set_ylabel("1990 M$")
        ax2.set_ylabel("1990 M$")

    fig.tight_layout()
    fig2.tight_layout()

    # Save plots if option selected
    if args.save:
        fig.savefig("direct_cost_bar.pdf")
        fig2.savefig("cost_bar.pdf")
    else:
        plt.show()


def comp_new(args, mfile_list: List[str], inflate: float):
    """

    Plot bar chart for the new 2014 cost model.

    """
    # Setup figures
    labels = [
        "Land and Buildings",
        "TF Coils",
        "First wall and blanket",
        "Active maintenance\n and remote handling",
        "Vacuum vessel and\n liquid nitrogen plant",
        "CS and PF coils",
        "Cryoplant and\n distribution",
        "Electrical power supply\n and distribution",
        "Additional project\n expenditure",
        "Other subsystems",
    ]
    index = np.arange(len(labels))
    bar_width = 0.7 / len(mfile_list)
    fig, ax = plt.subplots()

    # Read cost data
    for id, item in enumerate(mfile_list):
        cost = np.zeros(12)
        cost[0] = item.data["s09"].get_scan(-1)  # Buildings
        cost[1] = item.data["s13"].get_scan(-1)  # Land
        cost[2] = item.data["s21"].get_scan(-1)  # TF Coils
        cost[3] = item.data["s27"].get_scan(-1)  # First wall and blanket
        cost[4] = item.data["s31"].get_scan(
            -1
        )  # Active maintenance and remote handling
        cost[5] = item.data["s34"].get_scan(
            -1
        )  # Vacuum vessel and liquid nitrogen plant
        cost[6] = item.data["s35"].get_scan(
            -1
        )  # System for converting heat to electricity
        cost[7] = item.data["s36"].get_scan(-1)  # CS and PF coils
        cost[8] = item.data["s51"].get_scan(-1)  # Cryoplant and distribution
        cost[9] = item.data["s52"].get_scan(
            -1
        )  # Electrical power supply and distribution
        cost[10] = item.data["s59"].get_scan(-1)  # Additional project expenditure
        cost[11] = item.data["s61"].get_scan(-1)  # Remaining subsystems

        # Inflate costs using value parsed if specified
        if args.inf:
            cost = inflate * cost

        # Split up Remaining Subsystems as it is too large
        sizes = [
            cost[0] + cost[1],
            cost[2],
            cost[3],
            cost[4],
            cost[5],
            cost[7],
            cost[8],
            cost[9],
            cost[10],
            cost[6] + cost[11] - cost[7] - cost[8] - cost[9] - cost[10],
        ]

        # Plot bar chart
        ax.bar(index + id * bar_width, sizes, bar_width, label=args.f[id])

    # Plot labels
    ax.set_xticks(index + (len(mfile_list) - 1) * 0.5 * bar_width)
    ax.set_xticklabels(labels, rotation=90)
    ax.legend()

    # Adjust axis label depending on if inflation factor is used
    if args.inf:
        ax.set_ylabel("%.2f x (2014 M$)" % inflate)
    else:
        ax.set_ylabel("2014 M$")

    fig.tight_layout()

    # Save plots if option selected
    if args.save:
        fig.savefig("cost_bar.pdf")
    else:
        plt.show()


def main(args=None):
    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Displays the cost breakdown as a bar chart.  "
        "Multiple MFILEs can be given and will be plotted on the same chart.  "
        "For more information contact Stuart.Muldrew@ukaea.uk"
    )

    parser.add_argument(
        "-f", metavar="f", type=str, nargs="+", help="specify the MFILE(s) to plot"
    )

    parser.add_argument("-s", "--save", help="save figure", action="store_true")

    parser.add_argument(
        "-inf", type=float, help="Inflation Factor (multiplies costs)", default=1.0
    )

    args = parser.parse_args(args)

    # Get inflation factor if specified
    inflate = args.inf

    # Get file names
    mfile_list = list()
    for item in args.f:
        mfile_list.append(mf.MFile(filename=item))

    # Check which cost model is being used
    if "c21" in mfile_list[0].data.keys():
        # Check all MFILEs use original cost model
        for item in mfile_list:
            if "c21" not in item.data.keys():
                sys.exit("ERROR: Inconsistent cost models used between MFILEs!")

        comp_orig(args=args, mfile_list=mfile_list, inflate=inflate)

    elif "s01" in mfile_list[0].data.keys():
        # Check all MFILEs use new cost model
        for item in mfile_list:
            if "s01" not in item.data.keys():
                sys.exit("ERROR: Inconsistent cost models used between MFILEs!")

        comp_new(args=args, mfile_list=mfile_list, inflate=inflate)

    else:
        print("ERROR: Failed to identify cost data, check MFILE!")


if __name__ == "__main__":
    main()
