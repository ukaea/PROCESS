import sys

import click

import process.io.mfile.mfile as mf
from process.io.plot.costs.costs_bar import cost_comp_1990, cost_comp_2014
from process.io.plot.costs.costs_pie import new_cost_model, orig_cost_model
from process.io.tools import mfile_arg, mfile_opt, save

save = save("Save figure")


@click.group()
def costs():
    """Cost plotting utilities"""


@costs.command("pie", no_args_is_help=True)
@mfile_opt(exists=True)
@save
def pie_plot(mfile, save):
    """Displays the cost breakdown as a pie chart."""

    m_file = mf.MFile(mfile)

    # Check which cost model is being used
    if "c21" in m_file.data:
        orig_cost_model(m_file, save)
    elif "s01" in m_file.data:
        new_cost_model(m_file, save)
    else:
        print("ERROR: Cannot identify cost data, check MFILE!")


@costs.command("bar", no_args_is_help=True)
@mfile_arg
@save
@click.option(
    "-inf",
    "--inflate",
    type=float,
    help="Inflation Factor (multiplies costs)",
    default=1.0,
)
def bar_plot(mfile, save, inflate):
    """Displays the cost breakdown as a bar chart.

    Multiple MFILEs can be given and will be plotted on the same chart.
    """
    # Get file names
    mfile_list = [mf.MFile(filename=item) for item in mfile]

    # Check which cost model is being used
    if "c21" in mfile_list[0].data:
        # Check all MFILEs use original cost model
        for item in mfile_list:
            if "c21" not in item.data:
                sys.exit("ERROR: Inconsistent cost models used between MFILEs!")

        cost_comp_1990(mfile_list=mfile_list, inflate=inflate, save=save)

    elif "s01" in mfile_list[0].data:
        # Check all MFILEs use new cost model
        for item in mfile_list:
            if "s01" not in item.data:
                sys.exit("ERROR: Inconsistent cost models used between MFILEs!")

        cost_comp_2014(mfile_list=mfile_list, inflate=inflate, save=save)

    else:
        print("ERROR: Failed to identify cost data, check MFILE!")
