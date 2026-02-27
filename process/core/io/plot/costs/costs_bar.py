"""

Code to produce costs bar chart.
Can take multiple input files.

"""

from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np

import process.core.io.mfile.mfile as mf


def _format_fig(ax, fig, label, save, filename, index, inflate, ylabel_suffix, n_mfiles):
    ax.set_xticks(index + (n_mfiles - 1) * 0.5 * (0.7 / n_mfiles))

    ax.set_xticklabels(label, rotation=90)
    ax.legend()

    if inflate:
        ax.set_ylabel(f"{inflate:.2f} x {ylabel_suffix}")
    else:
        ax.set_ylabel(ylabel_suffix)

    fig.tight_layout()

    if save:
        fig.savefig(filename)
    else:
        plt.show()


def cost_comp_1990(
    mfile_list: list[mf.MFile], inflate: float = 1, save: bool = False
) -> None:
    """
    Plot bar chart for the orginal 1990 cost model.
    Two plots produced: (1) Breakdown of the direct costs and (2) Direct, indirect, etc.
    """
    fnames = ["direct_cost_bar.pdf", "cost_bar.pdf"]
    n_mfiles = len(mfile_list)
    bar_width = 0.7 / n_mfiles
    ylabel_suffix = "(1990 M$)"
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()

    variables = (
        "c21",  # Site and Buildings
        "c221",  # Reactor Systems
        "c222",  # Magnets
        "c223",  # Power Injection
        "c224",  # Vacuum Systems
        "c225",  # Power Conditioning
        "c226",  # Heat Transport System
        "c227",  # Fuel Handling System
        "c228",  # Instrumentation and Control
        "c229",  # Maintenance Equipment
        "c23",  # Turbine Plant Equipment
        "c24",  # Electric Plant Equipment
        "c25",  # Miscellaneous Plant Equipment
        "c26",  # Heat Rejection System
        "cdirt",  # Plant Direct Cost
        "c9",  # Indirect Cost
        "ccont",  # Total Contingency
        "moneyint",  # Interest during Construction
    )
    labels = [
        [
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
        ],
        [
            "Plant Direct\n Cost",
            "Indirect\n Cost",
            "Total\n Contingency",
            "Interest during\n Construction",
        ],
    ]
    index = np.arange(len(labels[0]))
    index2 = np.arange(len(labels[1]))

    # Read cost data
    for identity, item in enumerate(mfile_list):
        cost = np.array(item.get_variables(*variables, scan=-1), dtype=float)

        # Explain why moneyint is missing
        if "moneyint" not in item.data:
            print(
                "Interest during construction (moneyint) is only calculated for ireactor = 1"
            )

        if inflate:  # Inflate costs using value parsed if specified
            cost *= inflate

        # Simplify grouping
        sizes = [
            cost[2] + cost[5],
            *itemgetter(0, 9, 3, 1, 7, 8, 10, 6)(cost),
            sum(itemgetter(4, 11, 12, 13)(cost)),
        ]

        # Direct, indirect costs etc. for second plot
        sizes2 = itemgetter(14, 15, 16, 17)(cost)

        # Plot bar charts
        ax.bar(index + identity * bar_width, sizes, bar_width, label=item.filename)
        ax2.bar(index2 + identity * bar_width, sizes2, bar_width, label=item.filename)

    for ax_, fig_, lab, save_name, ind in zip(
        [ax, ax2], [fig, fig2], labels, fnames, [index, index2], strict=True
    ):
        _format_fig(
            ax_, fig_, lab, save, save_name, ind, inflate, ylabel_suffix, n_mfiles
        )


def cost_comp_2014(mfile_list: list[mf.MFile], inflate: float = 1, save: bool = False):
    """Plot bar chart for the new 2014 cost model."""
    variables = (
        "s09",  # Buildings
        "s13",  # Land
        "s21",  # TF Coils
        "s27",  # First wall and blanket
        "s31",  # Active maintenance and remote handling
        "s34",  # Vacuum vessel and liquid nitrogen plant
        "s35",  # System for converting heat to electricity
        "s36",  # CS and PF coils
        "s51",  # Cryoplant and distribution
        "s52",  # Electrical power supply and distribution
        "s59",  # Additional project expenditure
        "s61",  # Remaining subsystems
    )
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
    n_mfiles = len(mfile_list)
    bar_width = 0.7 / n_mfiles
    index = np.arange(len(labels))
    fig, ax = plt.subplots()

    # Read cost data
    for identity, item in enumerate(mfile_list):
        cost = np.array(item.get_variables(*variables, scan=-1), dtype=float)

        # Inflate costs using value parsed if specified
        if inflate:
            cost *= inflate

        # Split up Remaining Subsystems as it is too large
        sizes = [
            cost[0] + cost[1],
            *cost[2:6],
            *cost[7:11],
            cost[6] + cost[11] - cost[7] - cost[8] - cost[9] - cost[10],
        ]

        # Plot bar chart
        ax.bar(index + identity * bar_width, sizes, bar_width, label=item.filename)

    _format_fig(
        ax, fig, labels, save, "cost_bar.pdf", index, inflate, "(2014 M$)", n_mfiles
    )
