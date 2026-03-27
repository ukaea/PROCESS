import json

import click

from process.core.io.mfile import MFile
from process.core.io.mfile.mfile_comparison import compare_mfiles
from process.core.io.tools import mfile_arg, mfile_opt, save, scan_opt


@click.group()
def mfile():
    """MFile tools"""


@mfile.command("convert", no_args_is_help=True)
@mfile_opt(exists=True)
@click.option(
    "-v",
    "--variables",
    type=str,
    help="Optional list of variables or json file with list of variables to extract",
)
@click.option(
    "-fmt",
    "--format",
    "format_",
    type=click.Choice(["json", "csv", "toml"]),
)
@scan_opt
@click.option("--verbose", is_flag=True)
def convert(mfile, variables, format_, scan, verbose):
    """Convert MFile to other formats."""
    if variables.endswith(".json"):
        with open(variables) as f:
            variables = json.load(f)["variables"]
    else:
        variables = list(filter(None, variables.replace(" ", ":").split(":")))

    getattr(MFile(mfile), f"to_{format_}")(
        keys_to_write=variables, scan=scan, verbose=verbose
    )


@mfile.command("compare", no_args_is_help=True)
@mfile_arg
@save("Save output to file called comp.txt")
@click.option(
    "-t",
    "--comparison-type",
    "comparison",
    type=click.Choice(["defaults", "baseline", "blanket", "generic", "all"]),
    default="all",
    help="Format to save the eqdsk file in.",
)
@click.option("-v", "--verbose", default=False, is_flag=True)
@click.option("--acc", type=float, default=5.0)
def compare(mfiles, save, comparison, acc, verbose):
    """Produce a comparison between two PROCESS MFILEs.

    User Can specify level of differences to show.

    """
    compare_mfiles(mfiles, comparison, acc, save, verbose)
