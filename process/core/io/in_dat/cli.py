import click

from process.core.io.in_dat.create import write_indat
from process.core.io.tools import indat_opt, mfile_opt


@click.command("indat", no_args_is_help=True)
@mfile_opt(exists=True)
@indat_opt()
@click.option(
    "-o",
    "--indat-out",
    "indat_out",
    type=str,
    default="new_IN.DAT",
    help="IN.DAT to write out",
)
@click.option(
    "-fpi",
    "--feasuble-point-index",
    type=int,
    default=-1,
    help="Create indat from the Nth feasible point in mfile",
)
def new_indat(mfile, indat, indat_out, feasuble_point_index):
    """Creates a new IN.DAT using MFILE.DAT iteration variables."""

    write_indat(mfile, indat, indat_out, feasuble_point_index)
