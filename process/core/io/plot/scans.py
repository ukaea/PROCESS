"""Python utility for plotting the output of a PROCESS scan.

Depending of the type of scans, different actions will be taken:
1D SCANS: a simple graph using the scanned variable for x axis
and the selected variable on the y axis.
- Any number of output variables can be selected, a plot will be
made for each
- Several inputs files can be used at the same time if the same variable
is scanned. The different runs results will be plotted in the same
graph.
- If several inputs are used, the folder name or the file is used as
a legend

- 2D SCANS: n_scan_1 graph will be plotted using the second scanned variable
as x axis and the selected output as y axis
- Only one 2D scan can be ploted at once.

Performed checks:
- Non converged points are not plotted
- Only outputs existing in the MFILE.DAT are plotted
- No plot is made if the MFILE does not exists
- If the file is a folder, the contained MFILE is used as an input.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, PercentFormatter

from process.core.io.mfile import MFile
from process.core.io.mfile.cli import mfile
from process.core.io.variable_metadata import var_dicts as meta
from process.core.scan import ScanVariables

if TYPE_CHECKING:
    from matproplib.axes import Axes


@dataclass
class AxisData:
    name: str
    percent: bool
    max_: Sequence[float]
    range_: Sequence[float]
    tick_size: float
    font_size: float
    legend_size: float = 12


class AxisChoice(Enum):
    X = auto()
    Y = auto()

    def axis(self, ax):
        return getattr(ax, f"{self.name.lower()}axis")

    def set_lim(self, ax, lower, upper):
        getattr(ax, f"set_{self.name.lower()}lim")(lower, upper)


def get_list_padded(inp, names):
    target_len = len(names)
    inp_array = np.array(inp, dtype=float)
    if (i_len := len(inp_array)) < target_len:
        if i_len == 0:
            return [None] * target_len

        return np.concatenate((
            inp_array,
            np.full(target_len - i_len, inp_array[-1], dtype=float),
        ))
    return inp_array[:target_len]


def value_checks(
    scan_var: ScanVariables,
    scan_2_var: ScanVariables | None,
    m_file: MFile,
    input_files: Sequence[Path],
):
    ve_string = (
        "`{}` does not exist in PROCESS dicts\n"
        " The scan variable is probably an upper/lower boundary\n"
        " Please modify 'nsweep_dict' dict with the constrained var"
    )
    # Check if the scan variable is present
    if scan_var.name not in m_file.data:
        raise ValueError(ve_string.format(scan_var.name))

    # Check if the second scan variable is present
    if scan_2_var is not None:
        if scan_2_var.name not in m_file.data:
            raise ValueError(ve_string.format(scan_2_var.name))

        if len(input_files) > 1:
            raise ValueError("Only one input file can be used for 2D scans")


def array_check(output_name: str, m_file: MFile) -> bool:
    # Check if the output variable exists in the MFILE
    if output_name not in m_file.data:
        print(
            f"Warning : `{output_name}` does not exist in PROCESS dicts\n"
            f"Warning : `{output_name}` will not be output"
        )
        return False
    return True


def create_o_array(
    n_scan: int, m_file: MFile, output_name: str, conv_i: list[int]
) -> np.ndarray:
    return np.array([m_file.get(output_name, scan=conv_i[ii]) for ii in range(n_scan)])


def get_label(name: str) -> str:
    return meta[name].latex if name in meta else f"{name}"


def axis_manipulation(ax: Axes, axis: AxisData, index: int, contour: np.ndarray):

    an = AxisChoice[axis.name.upper()]

    if len(axis.range_) > 0:
        divisions = (axis.range_[1] - axis.range_[0]) / 10
    if axis.percent:
        if axis.max_[index] is None:
            axis.max_[index] = max(np.abs(contour))
        ticks = PercentFormatter(axis.max_[index])
        if len(axis.range_) > 0:
            scale = axis.max_[index] / 100
            divisions = 5 * math.ceil(divisions / 5) * scale
            range_ = (axis.range_[0] * scale, axis.range_[1] * scale)
        an.axis(ax).set_major_formatter(ticks)

    if len(axis.range_) > 0:
        if axis.percent is False:
            range_ = axis.range_
        an.set_lim(ax, range_[0], range_[1])
        an.axis(ax).set_major_locator(MultipleLocator(divisions))

    ax.figure.tight_layout()
    ax.tick_params(axis=an.name.lower(), labelsize=axis.tick_size)


def plot_scan(
    mfiles: Path | Iterable[Path],
    output_names: Sequence[str] = (),
    output_names2: Sequence[str] = (),
    outputdir: Path | None = None,
    term_output: bool = False,
    save_format: str = "pdf",
    axis_font_size: float = 18,
    axis_tick_size: float = 16,
    x_axis_percent: bool = False,
    x_axis_max: Sequence[float] = (),
    x_axis_range: Sequence[float] = (),
    y_axis_percent: bool = False,
    y_axis_max: Sequence[float] = (),
    y_axis_range: Sequence[float] = (),
    y_axis_percent2: bool = False,
    y_axis2_max: Sequence[float] = (),
    y_axis_range2: Sequence[float] = (),
    label_name: Sequence[str] = (),
    twod_contour: bool = False,
    stack_plots: bool = False,
):
    """Main plot scans script."""
    outputdir = outputdir or Path.cwd()
    input_files = mfiles if isinstance(mfiles, Iterable) else [mfiles]

    # If the input file is a directory, add MFILE.DAT
    for ii, if_ in enumerate(input_files):
        if if_.is_dir():
            input_files[ii] = if_ / "MFILE.DAT"

    # Getting the scanned variable name
    m_file = MFile(filename=input_files[-1])
    nsweep_ref = int(m_file.get("nsweep", scan=-1))
    scan_var = ScanVariables(nsweep_ref)

    # Get the eventual second scan variable
    scan_2_var = (
        ScanVariables(int(m_file.get("nsweep_2", scan=-1)))
        if "nsweep_2" in m_file.data
        else None
    )

    value_checks(scan_var, scan_2_var, m_file, input_files)

    x_max = get_list_padded(x_axis_max, output_names)
    x_axis = AxisData(
        "x", x_axis_percent, x_max, x_axis_range, axis_tick_size, axis_font_size
    )

    y_max = get_list_padded(y_axis_max, output_names)
    y_axis = AxisData(
        "y", y_axis_percent, y_max, y_axis_range, axis_tick_size, axis_font_size
    )
    if scan_2_var is None:
        y_axis2 = AxisData(
            "y",
            y_axis_percent2,
            (
                get_list_padded(y_axis2_max, output_names)
                if len(output_names2) > 0
                else y_axis2_max
            ),
            y_axis_range2,
            axis_tick_size,
            axis_font_size,
        )

        scan_var_array, output_arrays, output_arrays2 = oned_scan(
            input_files,
            nsweep_ref,
            scan_var,
            output_names,
            output_names2,
            term_output=term_output,
        )
        plot_1d_scan(
            input_files,
            m_file,
            output_names,
            output_names2,
            scan_var,
            outputdir,
            save_format,
            label_name,
            scan_var_array,
            output_arrays,
            output_arrays2,
            x_axis,
            y_axis,
            y_axis2,
            stack_plots=stack_plots,
        )
    else:
        twod_scan(
            input_files,
            scan_var,
            scan_2_var,
            output_names,
            outputdir,
            save_format,
            x_axis,
            y_axis,
            twod_contour=twod_contour,
        )


def oned_scan(
    input_files: Sequence[Path],
    nsweep_ref: int,
    scan_var: ScanVariables,
    output_names: Sequence[str],
    output_names2: Sequence[str],
    *,
    term_output: bool,
) -> tuple[np.ndarray, ...]:
    # input file, output_name, scan
    output_arrays = {}
    # input file, output_name2, scan
    output_arrays2 = {}
    # input_file, scan
    scan_var_array = {}
    for input_file in input_files:
        # Opening the MFILE.DAT
        m_file = MFile(filename=input_file)
        n_scan = int(m_file.get("isweep", scan=-1))

        # Check if the the scan variable is the same for all inputs
        # Same scan var
        nsweep = int(m_file.get("nsweep", scan=-1))
        if nsweep != nsweep_ref:
            raise ValueError("You must use inputs files with the same scan variables\n")

        if "nsweep_2" in m_file.data:
            raise ValueError("You cannot mix 1D with 2D scans\nERROR : Exiting")

        # Only selecting the scans that has converged
        # Converged indexes
        conv_i = []
        for ii in range(n_scan):
            ifail = m_file.get("ifail", scan=ii + 1)
            if ifail == 1:
                conv_i.append(ii + 1)
            else:
                failed_value = scan_var.get_val(m_file, scan=ii + 1)
                print(
                    f"Warning : Non-convergent scan point : {scan_var.name} = {failed_value}\n"
                    "Warning : This point will not be shown."
                )

        # Updating the number of scans
        n_scan = len(conv_i)
        scan_var_array[input_file] = np.array([
            scan_var.get_val(mfile, scan=conv_i[ii]) for ii in range(n_scan)
        ])
        output_arrays[input_file] = {
            output_name: create_o_array(n_scan, m_file, output_name, conv_i)
            for output_name in output_names
            if array_check(output_name, m_file)
        }
        output_arrays2[input_file] = {
            output_name2: create_o_array(n_scan, m_file, output_name2, conv_i)
            for output_name2 in output_names2
            if array_check(output_name2, m_file)
        }
        # Terminal output
        if term_output:
            print(
                f"\n{input_file} scan output\n\nX-axis:\n"
                f"scan var {scan_var.name} : {scan_var_array[input_file]}\n\nY-axis:"
                + "\n".join(
                    f"{output_name} : {output_arrays[input_file][output_name]}"
                    for output_name in output_names
                    if output_name in m_file.data
                )
                + "\n"
            )
            if len(output_names2) > 0:
                last_name = output_names2[-1]
                print(
                    f"Y2-Axis\n  {last_name} : {output_arrays2[input_file][last_name]}\n"
                )
    return scan_var_array, output_arrays, output_arrays2


def plot_1d_scan(
    input_files: Sequence[Path],
    m_file: MFile,
    output_names: Sequence[str],
    output_names2: Sequence[str],
    scan_var: ScanVariables,
    outputdir: Path,
    save_format: str,
    label_name: Sequence[str],
    scan_var_array: np.ndarray,
    output_arrays: np.ndarray,
    output_arrays2: np.ndarray,
    x_axis: AxisData,
    y_axis: AxisData,
    y_axis2: AxisData,
    *,
    stack_plots: bool,
):

    if stack_plots:
        # check stack plots will work
        if len(output_names) <= 1:
            raise ValueError("stack_plots requires at least two output variables")
        # Create subplots only once for the first output
        fig, axs = plt.subplots(
            len(output_names),
            1,
            figsize=(8.0, (3.5 + (1 * len(output_names)))),
            sharex=True,
        )
        fig.subplots_adjust(hspace=0.0)

    colour = (  # be careful changing this
        ("blue" if len(output_names2) > 0 else None)
        if len(output_names2) <= 0 or stack_plots
        else ("blue" if len(input_files) == 1 else None)
    )

    for index, output_name in enumerate(output_names):
        # reset counter for label_name
        kk = 0

        if output_name not in m_file.data:
            continue

        if stack_plots:
            ax_ = axs[index]
            ax = axs[output_names.index(output_name)]
        else:
            fig, ax = plt.subplots()
            if len(output_names2) > 0:
                ax2 = ax.twinx()
            ax_ = ax

        for input_file in input_files:
            if len(label_name) == 0:
                labl = input_file.name
            else:
                labl = label_name[kk]
                kk += 1

            ax_.plot(
                scan_var_array[input_file],
                output_arrays[input_file][output_name],
                "--o",
                color=colour,
                label=labl,
            )

            axis_manipulation(
                ax, axis=x_axis, index=index, contour=scan_var_array[input_file]
            )
            axis_manipulation(
                ax,
                axis=y_axis,
                index=index,
                contour=output_arrays[input_file][output_name],
            )

            if len(output_names2) > 0:
                for output_name2 in output_names2:
                    yval = output_arrays2[input_file][output_name2]
                    colour = "red" if len(input_files) == 1 else None
                    ax2.plot(
                        scan_var_array[input_file],
                        yval,
                        "--o",
                        color=colour,
                        label=labl,
                    )
                    ax2.set_ylabel(
                        get_label(output_name2),
                        fontsize=y_axis.font_size,
                        color=colour or "black",
                    )
                    axis_manipulation(ax2, axis=y_axis2, index=index, contour=yval)
        if len(output_names2) > 0:
            ax2.yaxis.grid(True)
            ax.xaxis.grid(True)
            ax.set_ylabel(
                get_label(output_name),
                fontsize=y_axis.font_size,
                color="blue" if len(input_files) == 1 else "black",
            )
            ax.set_xlabel(
                get_label(scan_var.name),
                fontsize=x_axis.font_size,
            )
            if len(input_files) != 1:
                fig.legend(loc="best", fontsize=x_axis.legend_size)
        elif stack_plots:
            ax.minorticks_on()
            ax.grid(True)
            ax.set_ylabel(get_label(output_name))
            ax.set_xlabel(get_label(scan_var.name), fontsize=x_axis.font_size)

            ymin, ymax = ax.get_ylim()
            if ymin < 0 and ymax > 0:
                mod_min = ymin * 1.1
                mod_max = ymax * 1.1
            elif ymin >= 0:
                mod_min = ymin * 0.9
                mod_max = ymax * 1.1
            else:
                mod_min = ymin * 1.1
                mod_max = ymax * 0.9
            ax.set_ylim(mod_min, mod_max)

            if len(input_files) > 1:
                fig.legend(
                    loc="lower center",
                    fontsize=x_axis.legend_size,
                    bbox_to_anchor=(0.5, -0.5 - (0.1 * len(output_names))),
                    fancybox=True,
                    shadow=False,
                    ncol=len(input_files),
                    columnspacing=0.8,
                )

        else:
            ax.grid(True)
            ax.set_ylabel(
                get_label(output_name),
                fontsize=x_axis.font_size,
                color="red" if len(output_names2) > 0 else "black",
            )
            ax.set_xlabel(get_label(scan_var.name), fontsize=x_axis.font_size)

            fig.title(
                f"{get_label(output_name)} vs {get_label(scan_var.name)}",
                fontsize=x_axis.font_size,
            )
            if len(input_files) != 1:
                fig.legend(loc="best", fontsize=x_axis.legend_size)

        ax.tick_params(axis=x_axis.name.lower(), labelsize=x_axis.tick_size)
        ax.tick_params(axis=y_axis.name.lower(), labelsize=y_axis.tick_size)
        fig.tight_layout()

        # Output file naming
        # This uses exclusively the last output_name2 defined in an earlier loop ignoring all other output_name2s...
        if output_name == "plasma_current_MA":
            extra_str = f"plasma_current{f'_vs_{output_name2}' if len(output_names2) > 0 else ''}"
        elif stack_plots and output_names[-1] == output_name:
            extra_str = f"{output_name}{f'_vs_{output_name2}' if len(output_names2) > 0 else '_vs_'.join(output_names)}"
        else:
            extra_str = (
                f"{output_name}{f'_vs_{output_name2}' if len(output_names2) > 0 else ''}"
            )

        if (not stack_plots) or (stack_plots and output_names[-1] == output_name):
            fig.savefig(
                f"{outputdir}/scan_{scan_var.name}_vs_{extra_str}.{save_format}",
                dpi=300,
            )
            plt.show()


def twod_scan(
    input_files: Sequence[Path],
    scan_var: ScanVariables,
    scan_2_var: ScanVariables,
    output_names: Sequence[str],
    outputdir: Path,
    save_format: str,
    x_axis: AxisData,
    y_axis: AxisData,
    *,
    twod_contour: bool,
):
    m_file = MFile(filename=input_files[0])

    # Number of scan points
    n_scan_1 = int(m_file.get("isweep", scan=-1))
    n_scan_2 = int(m_file.get("isweep_2", scan=-1))
    # Selecting the converged runs only
    contour_conv_ij = []  # List of non-converged scan point numbers
    conv_ij = []  # 2D array of converged scan point numbers (sweep = rows, sweep_2 = columns)
    ii_jj = 0
    for ii in range(n_scan_1):
        conv_ij.append([])
        for _jj in range(n_scan_2):
            ii_jj += 1  # Represents the scan point number in the MFILE
            ifail = m_file.get("ifail", scan=ii_jj)
            if ifail == 1:
                conv_ij[ii].append(ii_jj)  # Only appends scan number if scan converged
                contour_conv_ij.append(ii_jj)
            else:
                failed_value_1 = scan_var.get_val(m_file, scan=ii_jj)
                failed_value_2 = scan_2_var.get_val(m_file, scan=ii_jj)
                print(
                    f"Warning : Non-convergent scan point : ({scan_var.name},{scan_2_var.name}) "
                    f"= ({failed_value_1},{failed_value_2})\n"
                    "Warning : This point will not be shown."
                )

    for index, output_name in enumerate(output_names):
        if output_name not in m_file.data:
            print(
                f"Warning : `{output_name}` does not exist in PROCESS dicts\n"
                f"Warning : `{output_name}` will not be output"
            )
            continue

        fig, ax = plt.subplots()
        x_contour = [scan_2_var.get_val(m_file, scan=i + 1) for i in range(n_scan_2)]

        if twod_contour:
            y_contour = [
                scan_var.get_val(m_file, scan=i + 1)
                for i in range(1, n_scan_1 * n_scan_2, n_scan_2)
            ]

            output_contour_z = np.zeros((n_scan_1, n_scan_2))
            for i in contour_conv_ij:
                ind1 = (i - 1) // n_scan_2
                ind2 = (i - 1) % n_scan_2
                output_contour_z[ind1][(ind2 if ind1 % 2 == 0 else (-ind2 - 1))] = (
                    m_file.get(output_name, scan=i)
                )
            flat_output_z = output_contour_z.flatten()
            flat_output_z.sort()

            levels = np.linspace(
                next(filter(lambda i: i > 0.0, flat_output_z)), flat_output_z.max(), 50
            )
            contour = ax.contourf(x_contour, y_contour, output_contour_z, levels=levels)

            fig.colorbar(contour).set_label(
                label=get_label(output_name), size=y_axis.font_size
            )
            ax.set_ylabel(get_label(scan_var.name), fontsize=y_axis.font_size)

        else:
            y_contour = [m_file.get(output_name, scan=i + 1) for i in range(n_scan_2)]

            # conv_j is an array element containing the converged scan numbers
            for conv_j in conv_ij:
                # Scanned variables
                scan_1_var_array, scan_2_var_array, output_array = np.array([
                    (
                        scan_var.get_val(m_file, scan=conv_j[jj]),
                        scan_2_var.get_val(m_file, scan=conv_j[jj]),
                        m_file.get(output_name, scan=conv_j[jj]),
                    )
                    for jj in range(len(conv_j))
                ]).T

                label = f"{get_label(scan_var.name)} = {scan_1_var_array[0]}"
                ax.plot(scan_2_var_array, output_array, "--o", label=label)

            ax.set_ylabel(get_label(output_name), fontsize=y_axis.font_size)
            fig.legend(loc="best", fontsize=x_axis.legend_size)

        ax.set_xlabel(get_label(scan_2_var.name), fontsize=x_axis.font_size)

        axis_manipulation(ax, axis=x_axis, index=index, contour=x_contour)
        axis_manipulation(ax, axis=y_axis, index=index, contour=y_contour)
        ax.grid(True)
        fname = f"scan_{output_name}_vs_{scan_var.name}_{scan_2_var.name}.{save_format}"
        fig.savefig(outputdir / fname)
        plt.show()
