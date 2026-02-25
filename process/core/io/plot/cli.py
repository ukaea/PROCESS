from pathlib import Path

import click

from process.io.plot.plot_proc import setup_plot
from process.io.plot.plot_scans import plot_scan
from process.io.plot.plot_stress_tf import plot_stress
from process.io.plot.sankey import plot_sankey, plot_sankey_plotly
from process.io.tools import LazyGroup, mfile_arg, mfile_opt, split_callback


@click.group(
    cls=LazyGroup,
    lazy_subcommands={"costs": "process.io.plot.costs.cli.costs"},
)
def plot():
    """Plotting utilities for PROCESS"""


@plot.command("sankey", no_args_is_help=True)
@mfile_opt(exists=True)
@click.option("-fmt", "--format", "format_", default="pdf", help="file format []")
def sankey(mfile, format_):
    """Plot the power flow in PROCESS using a Sankey diagram."""
    if format_ in {"html", "plotly"}:
        out = plot_sankey_plotly(mfile)
        if out is not None:
            return out

    return plot_sankey(mfile, format_)


@plot.command("scans", no_args_is_help=True)
@mfile_arg
# At least one output variable must be supplied in order to plot
@click.option(
    "-yv",
    "--y-vars",
    "output_names",
    callback=split_callback,
    required=True,
    help=(
        "Select the output variables\nMore than one output can be plotted "
        "eg: -yv 'var1:var2'\nA separate plot will be created for each "
        "inputs"
    ),
)
@click.option(
    "-yv2",
    "--y-vars2",
    "output_names2",
    callback=split_callback,
    default="",
    help=(
        "Select the 2nd axis output variable\n "
        "eg: -yv2 'var'\n 2nd variable will be plotted on shared figure "
        "inputs"
    ),
)
@click.option(
    "-o",
    "--outputdir",
    default=Path.cwd(),
    type=click.Path(),
    help="Output directory for plots, defaults to current working directory.",
)
@click.option(
    "-out",
    "--term-output",
    is_flag=True,
    help="Option to show scans values on terminal",
)
@click.option(
    "-sf",
    "--save-format",
    default="pdf",
    help="Output format (default='pdf') ",
)
@click.option(
    "-afs",
    "--axis-font-size",
    default=18,
    help="Axis label font size selection (default=18)",
    type=int,
)
@click.option(
    "-ats",
    "--axis-tick-size",
    default=16,
    help="Axis tick label font size selection (default=16)",
    type=int,
)
@click.option(
    "-x%",
    "--x-axis-percent",
    is_flag=True,
    help=("Used to set the x axis ticks to percentages in place of absolute \nvalues."),
)
@click.option(
    "-xm",
    "--x-axis-max",
    callback=split_callback,
    default="",
    help=(
        "Used to set the x value corresponding to 100 percent when "
        "converting from absolute to percent values. Multiple values separated with :"
    ),
)
@click.option(
    "-xr",
    "--x-axis-range",
    callback=split_callback,
    default="",
    help=("Used to set the range for x axis"),
)
@click.option(
    "-y%",
    "--y-axis-percent",
    is_flag=True,
    help=("Used to set the y axis ticks to percentages in place of absolute \nvalues."),
)
@click.option(
    "-y2%",
    "--y-axis2-percent",
    "y_axis_percent2",
    is_flag=True,
    help=(
        "Used to set the y axis ticks to percentages in place of absolute \nvalues. For the twinned axis if present."
    ),
)
@click.option(
    "-ym",
    "--y-axis-max",
    callback=split_callback,
    default="",
    help=(
        "Used to set the y value corresponding to 100 percent when \nconverting from absolute to percent values."
    ),
)
@click.option(
    "-ym2",
    "--y-axis2-max",
    callback=split_callback,
    default="",
    help=(
        "Used to set the y value corresponding to 100 percent when \nconverting from absolute to percent values."
        "For the twinned axis if present."
    ),
)
@click.option(
    "-yr",
    "--y-axis-range",
    callback=split_callback,
    default="",
    help=("Used to set the range for y axis"),
)
@click.option(
    "-yr2",
    "--y-axis2-range",
    "y_axis_range2",
    callback=split_callback,
    default="",
    help=("Used to set the range for y axis. For the twinned axis if present."),
)
@click.option(
    "-ln",
    "--label-name",
    default="",
    callback=split_callback,
    help=(
        "Label names for plot legend. If multiple input files used then \n"
        "list the same number of label names eg: -nl 'leg1 leg2'\n"
        "(default = MFile file name) "
    ),
)
@click.option(
    "-2DC",
    "--two-dimensional-contour",
    "twoD_contour",
    is_flag=True,
    help=(
        "Option to plot 2D scans as a coloured contour plot instead of a line plot \n  "
        "Note: Non convergent points will show up with a value of zero \n "
        "Note: The scan paramters must both be in increasing orderl \n "
    ),
)
@click.option(
    "-stc",
    "--stack-plots",
    is_flag=True,
    help=(
        "Option to plot multiple 1D plots in a column of subplots \n  "
        "Variables will be plotted in order of input"
    ),
)
def plot_scans_cli(
    mfiles,
    output_names,
    output_names2,
    outputdir,
    term_output,
    save_format,
    axis_font_size,
    axis_tick_size,
    x_axis_percent,
    x_axis_max,
    x_axis_range,
    y_axis_percent,
    y_axis_percent2,
    y_axis_max,
    y_axis2_max,
    y_axis_range,
    y_axis_range2,
    label_name,
    twod_contour,
    stack_plots,
):
    """Plot optimisation information"""
    return plot_scan(
        mfiles,
        output_names,
        output_names2,
        outputdir,
        term_output,
        save_format,
        axis_font_size,
        axis_tick_size,
        x_axis_percent,
        list(map(float, x_axis_max)),
        list(map(float, x_axis_range)),
        y_axis_percent,
        y_axis_percent2,
        list(map(float, y_axis_max)),
        list(map(float, y_axis2_max)),
        list(map(float, y_axis_range)),
        list(map(float, y_axis_range2)),
        label_name,
        twod_contour,
        stack_plots,
    )


@plot.command("tf-stress", no_args_is_help=True)
@click.option(
    "-p",
    "--plot-selec",
    multiple=True,
    default=["all"],
    type=click.Choice(["all", "sig", "disp", "strain", "sm_sig"]),
    help="""\b
Plot selection string :
- If it containts 'sig'      -> Stress radial dependency
- If it containts 'strain'   -> Strain
- If it containts 'disp'     -> Displacement
- If it containts 'all'      -> all the mentioned plots (default value)
""",
)
@click.option(
    "-sf",
    "--save-format",
    default="pdf",
    help="output format (default='pdf') ",
)
@click.option(
    "-as",
    "--axis-font-size",
    default=18,
    help="Axis label font size selection (default=18)",
    type=int,
)
@click.option(
    "-out",
    "--term-output",
    is_flag=True,
    help="Option to show stress on terminal output",
)
@click.option(
    "-f",
    "--input-file",
    default="SIG_TF.json",
    help="specify input file path (default = SIG_TF.json)",
)
def plot_tf_stress(plot_selec, save_format, axis_font_size, term_output, input_file):
    """TF coil inboard mid-plane stress/strain summary plots"""
    plot_stress(plot_selec, save_format, axis_font_size, term_output, input_file)


@plot.command("summary", no_args_is_help=True)
@mfile_opt(exists=True)
@click.option("-n", "scan", type=int, default=-1, help="Which scan to plot?")
@click.option(
    "-d",
    "--DEMO-ranges",
    "demo_ranges",
    help="Uses the DEMO dimensions as ranges for all graphics",
    is_flag=True,
)
@click.option(
    "-c",
    "--colour",
    help=(
        "Which colour scheme to use for cross-section plots\n"
        "1: Original PROCESS (default)\n"
        "2: BLUEMIRA"
    ),
    default=1,
    type=click.Choice([1, 2]),
)
@click.option(
    "-o",
    "--output-format",
    help=(
        "Output file format\npdf: pdf output (default)\npng: png output\nnone: no output file written"
    ),
    default="pdf",
    type=click.Choice(["pdf", "png", "none"]),
)
@click.option("-s", "--show", help="show plot", is_flag=True)
def plot_proc(mfile, scan, demo_ranges, colour, output_format, show):
    """Produces a summary of the PROCESS MFILE output."""
    setup_plot(mfile, scan, demo_ranges, colour, output_format, show)
