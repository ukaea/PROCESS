"""
Code to display the power flow of a PROCESS run in a Sankey diagram

Input file:
MFILE.DAT
"""

import click

from process.io.plot.sankey.sankey_funcs import plot_sankey, plot_sankey_plotly
from process.io.tools import mfile_opt


@click.command("sankey", no_args_is_help=True)
@mfile_opt
@click.option("-fmt", "--format", "format_", default="pdf", help="file format []")
def sankey(mfile, format_):
    """Plot the power flow in PROCESS using a Sankey diagram."""
    if format_ in {"html", "plotly"}:
        out = plot_sankey_plotly(mfile)
        if out is not None:
            return out

    return plot_sankey(mfile, format_)
