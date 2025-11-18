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

import argparse
import math
import os
from argparse import RawTextHelpFormatter
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# PROCESS libraries
import process.io.mfile as mf
from process.io.variable_metadata import var_dicts as meta


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Plot optimization information",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "-f",
        "--input_files",
        default="MFILE.DAT",
        help=(
            "Specify input file(s) path(s) (default = MFILE.DAT)\n"
            "More than one input file can be used eg: -f 'A_MFILE.DAT "
            "B_MFILE.DAT'.\nYou can only specify the folder containing the "
            "MFILE.\nThe different files scan will be plotted on the same "
            "graph.\nThe scans must use the same scan variation."
        ),
    )

    # At least one output variable must be supplied in order to plot
    parser.add_argument(
        "-yv",
        "--y_vars",
        required=True,
        help=(
            "Select the output variables\nMore than one output can be plotted "
            "eg: -yv 'var1 var2'\nA separate plot will be created for each "
            "inputs"
        ),
    )

    parser.add_argument(
        "-yv2",
        "--y_vars2",
        default="",
        help=(
            "Select the 2nd axis output variable\n "
            "eg: -yv2 'var'\n 2nd variable will be plotted on shared figure "
            "inputs"
        ),
    )

    parser.add_argument(
        "-o",
        "--outputdir",
        default=Path.cwd(),
        help="Output directory for plots, defaults to current working directory.",
    )

    parser.add_argument(
        "-out",
        "--term_output",
        action="store_true",
        help="Option to show scans values on terminal",
    )

    parser.add_argument(
        "-sf",
        "--save_format",
        nargs="?",
        default="pdf",
        help="Output format (default='pdf') ",
    )

    parser.add_argument(
        "-afs",
        "--axis_font_size",
        nargs="?",
        default=18,
        help="Axis label font size selection (default=18)",
        type=int,
    )

    parser.add_argument(
        "-ats",
        "--axis_ticklabel_size",
        nargs="?",
        default=16,
        help="Axis tick label font size selection (default=16)",
        type=int,
    )

    parser.add_argument(
        "-x%",
        "--x_axis_percent",
        action="store_true",
        help=(
            "Used to set the x axis ticks to percentages in place of absolute \nvalues."
        ),
    )

    parser.add_argument(
        "-xm",
        "--x_axis_max",
        default="",
        help=(
            "Used to set the x value corresponding to 100 percent when \nconverting from absolute to percent values."
        ),
    )

    parser.add_argument(
        "-xr",
        "--x_axis_range",
        default="",
        help=("Used to set the range for x axis"),
    )

    parser.add_argument(
        "-y%",
        "--y_axis_percent",
        action="store_true",
        help=(
            "Used to set the y axis ticks to percentages in place of absolute \nvalues."
        ),
    )

    parser.add_argument(
        "-y2%",
        "--y_axis_2_percent",
        action="store_true",
        help=(
            "Used to set the y axis ticks to percentages in place of absolute \nvalues. For the twinned axis if present."
        ),
    )

    parser.add_argument(
        "-ym",
        "--y_axis_max",
        default="",
        help=(
            "Used to set the y value corresponding to 100 percent when \nconverting from absolute to percent values."
        ),
    )

    parser.add_argument(
        "-ym2",
        "--y_axis_2_max",
        default="",
        help=(
            "Used to set the y value corresponding to 100 percent when \nconverting from absolute to percent values."
            "For the twinned axis if present."
        ),
    )

    parser.add_argument(
        "-yr",
        "--y_axis_range",
        default="",
        help=("Used to set the range for y axis"),
    )

    parser.add_argument(
        "-yr2",
        "--y_axis_2_range",
        default="",
        help=("Used to set the range for y axis. For the twinned axis if present."),
    )

    parser.add_argument(
        "-ln",
        "--label_name",
        default="",
        help=(
            "Label names for plot legend. If multiple input files used then \n"
            "list the same number of label names eg: -nl 'leg1 leg2'\n"
            "(default = MFile file name) "
        ),
    )

    parser.add_argument(
        "-2DC",
        "--two_dimensional_contour",
        action="store_true",
        help=(
            "Option to plot 2D scans as a coloured contour plot instead of a line plot \n  "
            "Note: Non convergent points will show up with a value of zero \n "
            "Note: The scan paramters must both be in increasing orderl \n "
        ),
    )

    parser.add_argument(
        "-stc",
        "--stack_plots",
        action="store_true",
        help=(
            "Option to plot multiple 1D plots in a column of subplots \n  "
            "Variables will be plotted in order of input"
        ),
    )

    return parser.parse_args(args)


def main(args=None):
    """Main plot scans script.

    :param args: optional command-line args from test function, defaults to None
    :type args: list, optional
    """
    args = parse_args(args)

    # Parameters to be used as function input
    # ---------------------------------------
    input_files = str(args.input_files)
    output_names = str(args.y_vars)
    output_names2 = str(args.y_vars2)
    save_format = str(args.save_format)
    term_output = args.term_output
    label_name = str(args.label_name)
    x_axis_percentage = args.x_axis_percent
    x_max_input = list(filter(None, args.x_axis_max.split(" ")))
    y_axis_percentage = args.y_axis_percent
    y_max_input = list(filter(None, args.y_axis_max.split(" ")))
    y_axis_percentage2 = args.y_axis_2_percent
    y_max2_input = list(filter(None, args.y_axis_2_max.split(" ")))
    two_dimensional_contour = args.two_dimensional_contour
    stack_plots = args.stack_plots
    # ---------------------------------------

    # Input checks
    # ------------
    # Formting the inputs
    output_names = list(filter(None, output_names.split(" ")))
    output_names2 = list(filter(None, output_names2.split(" ")))
    input_files = list(filter(None, input_files.split(" ")))
    label_name = list(filter(None, label_name.split(" ")))

    # If the input file is a directory, add MFILE.DAT
    for ii in range(len(input_files)):
        if os.path.isdir(input_files[ii]):
            input_files[ii] = input_files[ii].replace("/", "")
            input_files[ii] = input_files[ii] + "/MFILE.DAT"

        # Check for the existence of the MFILE
        if not os.path.isfile(input_files[ii]):
            print(f"ERROR : The {input_files[ii]} MFILE does not exist, skipping it")
            input_files.remove(input_files[ii])

    # nsweep varible dict
    # -------------------
    # TODO WOULD BE GREAT TO HAVE IT AUTOMATICALLY GENERATED ON THE PROCESS CMAKE!
    #        THE SAME WAY THE DICTS ARE
    # This needs to be kept in sync automatically; this will break frequently
    # otherwise
    # Rem : Some variables are not in the MFILE, making the defintion rather tricky...
    nsweep_dict = {
        1: "aspect",
        2: "pflux_div_heat_load_max_mw",
        3: "p_plant_electric_net_mw",
        4: "hfact",
        5: "oacdcp",
        6: "pflux_fw_neutron_max_mw",
        7: "beamfus0",
        8: "fbig_q_plasma_min",
        9: "temp_plasma_electron_vol_avg_kev",
        10: "boundu(15)",
        11: "beta_norm_max",
        12: "f_c_plasma_bootstrap_max",
        13: "boundu(10)",
        14: "fiooic",
        15: "fjprot",
        16: "rmajor",
        17: "b_tf_inboard_peak_symmetric",  # b_tf_inboard_max the maximum T field upper limit is the scan variable
        18: "eta_cd_norm_hcd_primary_max",
        19: "boundl(16)",
        20: "cnstv.t_burn_min",
        21: "",
        22: "f_t_plant_available",
        23: "boundu(72)",
        24: "p_fusion_total_max_mw",
        25: "kappa",
        26: "triang",
        27: "tbrmin",
        28: "b_plasma_toroidal_on_axis",
        29: "radius_plasma_core_norm",
        30: "",  # OBSOLETE
        31: "f_alpha_energy_confinement_min",
        32: "epsvmc",
        33: "ttarget",
        34: "qtargettotal",
        35: "lambda_q_omp",
        36: "lambda_target",
        37: "lcon_factor",
        38: "boundu(129)",
        39: "boundu(131)",
        40: "boundu(135)",
        41: "dr_blkt_outboard",
        42: "f_nd_impurity_electrons(9)",
        43: "Obsolete",  # Removed
        44: "alstrtf",
        45: "temp_tf_superconductor_margin_min",
        46: "boundu(152)",
        47: "impurity_enrichment(9)",
        48: "n_tf_wp_pancakes",
        49: "n_tf_wp_layers",
        50: "f_nd_impurity_electrons(13)",
        51: "f_p_div_lower",
        52: "rad_fraction_sol",
        54: "b_crit_upper_nbti",
        55: "dr_shld_inboard",
        56: "p_cryo_plant_electric_max_mw",
        57: "b_plasma_toroidal_on_axis",  # Genuinly b_plasma_toroidal_on_axis lower bound
        58: "dr_fw_plasma_gap_inboard",
        59: "dr_fw_plasma_gap_outboard",
        60: "sig_tf_wp_max",
        61: "copperaoh_m2_max",
        62: "j_cs_flat_top_end",
        63: "dr_cs",
        64: "f_z_cs_tf_internal",
        65: "csfv.n_cycle_min",
        66: "pfv.f_a_cs_turn_steel",
        67: "csfv.t_crack_vertical",
        68: "inlet_temp_liq",
        69: "outlet_temp_liq",
        70: "blpressure_liq",
        71: "n_liq_recirc",
        72: "bz_channel_conduct_liq",
        73: "pnuc_fw_ratio_dcll",
        74: "f_nuc_pow_bz_struct",
        75: "dx_fw_module",
        76: "eta_turbine",
        77: "startupratio",
        78: "fkind",
        79: "eta_ecrh_injector_wall_plug",
        80: "fcoolcp",
        81: "n_tf_coil_turns",
    }
    # -------------------

    # Getting the scanned variable name
    m_file = mf.MFile(filename=input_files[-1])
    nsweep_ref = int(m_file.data["nsweep"].get_scan(-1))
    scan_var_name = nsweep_dict[nsweep_ref]
    # Get the eventual second scan variable
    nsweep_2_ref = 0
    is_2D_scan = False
    scan_2_var_name = ""
    if "nsweep_2" in m_file.data:
        is_2D_scan = True
        nsweep_2_ref = int(m_file.data["nsweep_2"].get_scan(-1))
        scan_2_var_name = nsweep_dict[nsweep_2_ref]

    # Checks
    # ------
    # Check if the nsweep dict has been updated
    if nsweep_ref > len(nsweep_dict) + 1:
        print(
            f"ERROR : nsweep = {nsweep_ref} not supported by the utility\n"
            "ERROR : Please update the 'nsweep_dict' dict"
        )
        exit()

    # Check if the scan variable is present in the
    if scan_var_name not in m_file.data:
        print(
            f"ERROR : `{scan_var_name}` does not exist in PROCESS dicts\n"
            "ERROR : The scan variable is probably an upper/lower boundary\n"
            "ERROR : Please modify 'nsweep_dict' dict with the constrained var"
        )
        exit()

    # Check if the second scan variable is present in the
    if is_2D_scan and (scan_2_var_name not in m_file.data):
        print(
            f"ERROR : `{scan_2_var_name}` does not exist in PROCESS dicts\n"
            "ERROR : The scan variable is probably an upper/lower boundary\n"
            "ERROR : Please modify 'nsweep_dict' dict with the constrained var"
        )
        exit()

    # Only one imput must be used for a 2D scan
    if is_2D_scan and len(input_files) > 1:
        print("ERROR : Only one input file can be used for 2D scans\nERROR : Exiting")
        exit()
    # ------

    # Plot settings
    # -------------
    # Plot cosmetic settings
    axis_tick_size = args.axis_ticklabel_size
    legend_size = 12
    axis_font_size = args.axis_font_size
    x_axis_range = list(filter(None, args.x_axis_range.split(" ")))
    if x_axis_range != []:
        x_axis_range = list(np.float64(x_axis_range))
    y_axis_range = list(filter(None, args.y_axis_range.split(" ")))
    if y_axis_range != []:
        y_axis_range = list(np.float64(y_axis_range))
    y_axis_range2 = list(filter(None, args.y_axis_2_range.split(" ")))
    if y_axis_range2 != []:
        y_axis_range2 = list(np.float64(y_axis_range2))

    if len(x_max_input) != len(output_names):
        x_max = []
        for i in range(len(output_names)):
            if x_max_input != []:
                j = 0
                try:
                    x_max += [float(x_max_input[i])]
                    j += 1
                except IndexError:
                    x_max += [float(x_max_input[j])]
            else:
                x_max += [None]
    else:
        x_max = np.float64(x_max_input)
    if len(y_max_input) != len(output_names):
        y_max = []
        for i in range(len(output_names)):
            if y_max_input != []:
                j = 0
                try:
                    y_max += [float(y_max_input[i])]
                    j += 1
                except IndexError:
                    y_max += [float(y_max_input[j])]
            else:
                y_max += [None]
    else:
        y_max = np.float64(y_max_input)
    if output_names2 != []:
        if len(y_max2_input) != len(output_names):
            y_max2 = []
            for i in range(len(output_names)):
                if y_max2_input != []:
                    j = 0
                    try:
                        y_max2 += [float(y_max2_input[i])]
                        j += 1
                    except IndexError:
                        y_max2 += [float(y_max2_input[j])]
                else:
                    y_max2 += [None]
        else:
            y_max2 = np.float64(y_max2_input)
    else:
        y_max2 = y_max2_input
    # -------------

    # Case of a set of 1D scans
    # ----------------------------------------------------------------------------------------------
    if not is_2D_scan:
        # Loop over the MFILEs
        output_arrays = {}
        output_arrays2 = {}
        scan_var_array = {}
        for input_file in input_files:
            # Opening the MFILE.DAT
            m_file = mf.MFile(filename=input_file)

            # Check if the the scan variable is the same for all inputs
            # ---
            # Same scan var
            nsweep = int(m_file.data["nsweep"].get_scan(-1))
            if nsweep != nsweep_ref:
                print(
                    "ERROR : You must use inputs files with the same scan variables\n"
                    "ERROR : Exiting"
                )
                exit()

            # No D scans
            if "nsweep_2" in m_file.data:
                print("ERROR : You cannot mix 1D with 2D scans\nERROR : Exiting")
                exit()
            # ---

            # Only selecting the scans that has converged
            # ---
            # Number of scan points
            n_scan = int(m_file.data["isweep"].get_scan(-1))

            # Converged indexes
            conv_i = []
            for ii in range(n_scan):
                ifail = m_file.data["ifail"].get_scan(ii + 1)
                if ifail == 1:
                    conv_i.append(ii + 1)
                else:
                    failed_value = m_file.data[scan_var_name].get_scan(ii + 1)
                    print(
                        f"Warning : Non-convergent scan point : {scan_var_name} = {failed_value}\n"
                        "Warning : This point will not be shown."
                    )

            # Updating the number of scans
            n_scan = len(conv_i)
            # ---
            # Scanned variable
            scan_var_array[input_file] = np.zeros(n_scan)

            for ii in range(n_scan):
                scan_var_array[input_file][ii] = m_file.data[scan_var_name].get_scan(
                    conv_i[ii]
                )
            # output list declaration
            output_arrays[input_file] = {}
            output_arrays2[input_file] = {}
            # First variable scan
            for output_name in output_names:
                ouput_array = np.zeros(n_scan)

                # Check if the output variable exists in the MFILE
                if output_name not in m_file.data:
                    print(
                        f"Warning : `{output_name}` does not exist in PROCESS dicts\n"
                        f"Warning : `{output_name}` will not be output"
                    )
                    continue

                for ii in range(n_scan):
                    ouput_array[ii] = m_file.data[output_name].get_scan(conv_i[ii])
                output_arrays[input_file][output_name] = ouput_array
            # Second variable scan
            for output_name2 in output_names2:
                ouput_array2 = np.zeros(n_scan)

                # Check if the output variable exists in the MFILE
                if output_name2 not in m_file.data:
                    print(
                        f"Warning : `{output_name2}` does not exist in PROCESS dicts\n"
                        f"Warning : `{output_name2}` will not be output"
                    )
                    continue

                for ii in range(n_scan):
                    ouput_array2[ii] = m_file.data[output_name2].get_scan(conv_i[ii])
                output_arrays2[input_file][output_name2] = ouput_array2
            # Terminal output
            if term_output:
                print(
                    f"\n{input_file} scan output\n\nX-axis:\n"
                    f"scan var {scan_var_name} : {scan_var_array[input_file]}\n\nY-axis:"
                )
                for output_name in output_names:
                    # Check if the output variable exists in the MFILE
                    if output_name not in m_file.data:
                        continue

                    print(f"{output_name} : {output_arrays[input_file][output_name]}")
                print()
                if output_names2 != []:
                    print(
                        f"Y2-Axis\n  {output_name2} : {output_arrays2[input_file][output_name2]}\n"
                    )
        # Plot section
        # -----------
        for index, output_name in enumerate(output_names):
            if stack_plots:
                # check stack plots will work
                if len(output_names) <= 1:
                    raise ValueError(
                        "For stack plots to be used need more than 1 output variable"
                    )
                fig, axs = plt.subplots(
                    len(output_names),
                    1,
                    figsize=(8.0, (3.5 + (1 * len(output_names)))),
                    sharex=True,
                )
                fig.subplots_adjust(hspace=0.0)
            else:
                fig, ax = plt.subplots()
                if output_names2 != []:
                    ax2 = ax.twinx()
            # reset counter for label_name
            kk = 0

            # Check if the output variable exists in the MFILE
            if output_name not in m_file.data:
                continue

            # Loop over inputs
            for input_file in input_files:
                # Legend label formating
                if label_name == []:
                    labl = input_file
                    if "/MFILE.DAT" in input_file:
                        labl = input_file[:-10]
                    elif "MFILE.DAT" in input_file:
                        labl = input_file[:-9]
                    labl = labl.replace("_", " ")
                else:
                    labl = label_name[kk]
                    kk = kk + 1

                # Plot the graph
                if output_names2 != [] and not stack_plots:
                    ax.plot(
                        scan_var_array[input_file],
                        output_arrays[input_file][output_name],
                        "--o",
                        color="blue" if len(input_files) == 1 else None,
                        label=labl,
                    )
                    if y_axis_range != []:
                        y_divisions = (y_axis_range[1] - y_axis_range[0]) / 10
                    if y_axis_percentage:
                        if y_max[index] is None:
                            y_max[index] = max(
                                np.abs(output_arrays[input_file][output_name])
                            )
                        yticks = mtick.PercentFormatter(y_max[index])
                        if y_axis_range != []:
                            y_divisions = (
                                5 * math.ceil(y_divisions / 5) * y_max[index] / 100
                            )
                            y_range = (
                                y_axis_range[0] * y_max[index] / 100,
                                y_axis_range[1] * y_max[index] / 100,
                            )
                        ax.yaxis.set_major_formatter(yticks)
                    if y_axis_range != []:
                        if y_axis_percentage is False:
                            y_range = y_axis_range
                        ax.set_ylim(y_range[0], y_range[1])
                        ax.yaxis.set_major_locator(mtick.MultipleLocator(y_divisions))
                    if x_axis_range != []:
                        x_divisions = (x_axis_range[1] - x_axis_range[0]) / 10
                    if x_axis_percentage:
                        if x_max[index] is None:
                            x_max[index] = max(np.abs(scan_var_array[input_file]))
                        xticks = mtick.PercentFormatter(x_max[index])
                        ax.xaxis.set_major_formatter(xticks)
                        if x_axis_range != []:
                            x_divisions = (
                                5 * math.ceil(x_divisions / 5) * x_max[index] / 100
                            )
                            x_range = (
                                x_axis_range[0] * x_max[index] / 100,
                                x_axis_range[1] * x_max[index] / 100,
                            )
                    plt.rc("xtick", labelsize=axis_tick_size)
                    plt.rc("ytick", labelsize=axis_tick_size)
                    if x_axis_range != []:
                        if x_axis_percentage is False:
                            x_range = x_axis_range
                        plt.xlim(x_range[0], x_range[1])
                        ax.xaxis.set_major_locator(mtick.MultipleLocator(x_divisions))
                    plt.tight_layout()
                else:
                    if stack_plots:
                        axs[output_names.index(output_name)].plot(
                            scan_var_array[input_file],
                            output_arrays[input_file][output_name],
                            "--o",
                            color="blue" if output_names2 != [] else None,
                            label=labl,
                        )
                        if y_axis_range != []:
                            y_divisions = (y_axis_range[1] - y_axis_range[0]) / 10
                        if y_axis_percentage:
                            if y_max[index] is None:
                                y_max[index] = max(
                                    np.abs(output_arrays[input_file][output_name])
                                )
                            yticks = mtick.PercentFormatter(y_max[index])
                            if y_axis_range != []:
                                y_divisions = (
                                    5 * math.ceil(y_divisions / 5) * y_max[index] / 100
                                )
                                y_range = (
                                    y_axis_range[0] * y_max[index] / 100,
                                    y_axis_range[1] * y_max[index] / 100,
                                )
                            axs[
                                output_names.index(output_name)
                            ].yaxis.set_major_formatter(yticks)
                        if y_axis_range != []:
                            if y_axis_percentage is False:
                                y_range = y_axis_range
                            axs[output_names.index(output_name)].set_ylim(
                                y_range[0], y_range[1]
                            )
                            axs[
                                output_names.index(output_name)
                            ].yaxis.set_major_locator(
                                mtick.MultipleLocator(y_divisions)
                            )
                        if x_axis_range != []:
                            x_divisions = (x_axis_range[1] - x_axis_range[0]) / 10
                        if x_axis_percentage:
                            if x_max[index] is None:
                                x_max[index] = max(np.abs(scan_var_array[input_file]))
                            xticks = mtick.PercentFormatter(x_max[index])
                            if x_axis_range != []:
                                x_divisions = (
                                    5 * math.ceil(x_divisions / 5) * x_max[index] / 100
                                )
                                x_range = (
                                    x_axis_range[0] * x_max[index] / 100,
                                    x_axis_range[1] * x_max[index] / 100,
                                )
                            axs[
                                output_names.index(output_name)
                            ].xaxis.set_major_formatter(xticks)
                        if x_axis_range != []:
                            if x_axis_percentage is False:
                                x_range = x_axis_range
                            plt.xlim(x_range[0], x_range[1])
                            axs[
                                output_names.index(output_name)
                            ].xaxis.set_major_locator(
                                mtick.MultipleLocator(x_divisions)
                            )
                        plt.rc("xtick", labelsize=axis_tick_size)
                        plt.rc("ytick", labelsize=axis_tick_size)
                        plt.tight_layout()
                    else:
                        ax.plot(
                            scan_var_array[input_file],
                            output_arrays[input_file][output_name],
                            "--o",
                            color="blue" if output_names2 != [] else None,
                            label=labl,
                        )
                        if y_axis_range != []:
                            y_divisions = (y_axis_range[1] - y_axis_range[0]) / 10
                        if y_axis_percentage:
                            if y_max[index] is None:
                                y_max[index] = max(
                                    np.abs(output_arrays[input_file][output_name])
                                )
                            yticks = mtick.PercentFormatter(y_max[index])
                            if y_axis_range != []:
                                y_divisions = (
                                    5 * math.ceil(y_divisions / 5) * y_max[index] / 100
                                )
                                y_range = (
                                    y_axis_range[0] * y_max[index] / 100,
                                    y_axis_range[1] * y_max[index] / 100,
                                )
                            ax.yaxis.set_major_formatter(yticks)
                        if y_axis_range != []:
                            if y_axis_percentage is False:
                                y_range = y_axis_range
                            ax.set_ylim(y_range[0], y_range[1])
                            ax.yaxis.set_major_locator(
                                mtick.MultipleLocator(y_divisions)
                            )
                        if x_axis_range != []:
                            x_divisions = (x_axis_range[1] - x_axis_range[0]) / 10
                        if x_axis_percentage:
                            if x_max[index] is None:
                                x_max[index] = max(np.abs(scan_var_array[input_file]))
                            xticks = mtick.PercentFormatter(x_max[index])
                            if x_axis_range != []:
                                x_divisions = (
                                    5 * math.ceil(x_divisions / 5) * x_max[index] / 100
                                )
                                x_range = (
                                    x_axis_range[0] * x_max[index] / 100,
                                    x_axis_range[1] * x_max[index] / 100,
                                )
                            ax.xaxis.set_major_formatter(xticks)
                        if x_axis_range != []:
                            if x_axis_percentage is False:
                                x_range = x_axis_range
                            plt.xlim(x_range[0], x_range[1])
                            ax.xaxis.set_major_locator(
                                mtick.MultipleLocator(x_divisions)
                            )
                        plt.rc("xtick", labelsize=axis_tick_size)
                        plt.rc("ytick", labelsize=axis_tick_size)
                        plt.tight_layout()
                if output_names2 != []:
                    ax2.plot(
                        scan_var_array[input_file],
                        output_arrays2[input_file][output_name2],
                        "--o",
                        color="red" if len(input_files) == 1 else None,
                        label=labl,
                    )
                    ax2.set_ylabel(
                        (
                            meta[output_name2].latex
                            if output_name2 in meta
                            else f"{output_name2}"
                        ),
                        fontsize=axis_font_size,
                        color="red" if len(input_files) == 1 else "black",
                    )
                    if y_axis_range2 != []:
                        y_divisions2 = (y_axis_range2[1] - y_axis_range2[0]) / 10
                    if y_axis_percentage2:
                        if y_max2[index] is None:
                            y_max2[index] = max(
                                np.abs(output_arrays2[input_file][output_name2])
                            )
                        yticks2 = mtick.PercentFormatter(y_max2[index])
                        if y_axis_range2 != []:
                            y_divisions2 = (
                                5 * math.ceil(y_divisions2 / 5) * y_max2[index] / 100
                            )
                            y_range2 = (
                                y_axis_range2[0] * y_max2[index] / 100,
                                y_axis_range2[1] * y_max2[index] / 100,
                            )
                        ax2.yaxis.set_major_formatter(yticks2)
                    if y_axis_range2 != []:
                        if y_axis_percentage2 is False:
                            y_range2 = y_axis_range2
                        ax2.set_ylim(y_range2[0], y_range2[1])
                        ax2.yaxis.set_major_locator(mtick.MultipleLocator(y_divisions2))
                    plt.rc("xtick", labelsize=axis_tick_size)
                    plt.rc("ytick", labelsize=axis_tick_size)
                    plt.tight_layout()
            if output_names2 != []:
                ax2.yaxis.grid(True)
                ax.xaxis.grid(True)
                ax.set_ylabel(
                    (
                        meta[output_name].latex
                        if output_name in meta
                        else f"{output_name}"
                    ),
                    fontsize=axis_font_size,
                    color="blue" if len(input_files) == 1 else "black",
                )
                ax.set_xlabel(
                    (
                        meta[scan_var_name].latex
                        if scan_var_name in meta
                        else f"{scan_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.rc("xtick", labelsize=axis_tick_size)
                plt.rc("ytick", labelsize=axis_tick_size)
                if len(input_files) != 1:
                    plt.legend(loc="best", fontsize=legend_size)
                plt.tight_layout()
            elif stack_plots:
                axs[output_names.index(output_name)].minorticks_on()
                axs[output_names.index(output_name)].grid(True)
                axs[output_names.index(output_name)].set_ylabel(
                    (
                        meta[output_name].latex
                        if output_name in meta
                        else f"{output_name}"
                    ),
                )

                plt.xlabel(
                    (
                        meta[scan_var_name].latex
                        if scan_var_name in meta
                        else f"{scan_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.rc("xtick", labelsize=axis_tick_size)
                plt.rc("ytick", labelsize=axis_tick_size)
                if len(input_files) > 1:
                    plt.legend(
                        loc="lower center",
                        fontsize=legend_size,
                        bbox_to_anchor=(0.5, -0.5 - (0.1 * len(output_names))),
                        # bbox_to_anchor=(0.5, -1.4),
                        fancybox=True,
                        shadow=False,
                        ncol=len(input_files),
                        columnspacing=0.8,
                    )
                plt.tight_layout()
                ymin, ymax = axs[output_names.index(output_name)].get_ylim()
                if ymin < 0 and ymax > 0:
                    axs[output_names.index(output_name)].set_ylim(
                        ymin * 1.1, ymax * 1.1
                    )
                elif ymin >= 0:
                    axs[output_names.index(output_name)].set_ylim(
                        ymin * 0.9, ymax * 1.1
                    )
                else:
                    axs[output_names.index(output_name)].set_ylim(
                        ymin * 1.1, ymax * 0.9
                    )
            else:
                plt.grid(True)
                plt.ylabel(
                    (
                        meta[output_name].latex
                        if output_name in meta
                        else f"{output_name}"
                    ),
                    fontsize=axis_font_size,
                    color="red" if output_names2 != [] else "black",
                )
                plt.xlabel(
                    (
                        meta[scan_var_name].latex
                        if scan_var_name in meta
                        else f"{scan_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.rc("xtick", labelsize=axis_tick_size)
                plt.rc("ytick", labelsize=axis_tick_size)
                plt.title(
                    f"{meta[output_name].latex if output_name in meta else {output_name}} vs "
                    f"{meta[scan_var_name].latex if scan_var_name in meta else {scan_var_name}}",
                    fontsize=axis_font_size,
                )
                plt.tight_layout()
                if len(input_files) != 1:
                    plt.legend(loc="best", fontsize=legend_size)
                    plt.tight_layout()

            # Output file naming
            if output_name == "plasma_current_MA":
                extra_str = f"plasma_current{f'_vs_{output_name2}' if output_names2 != [] else ''}"
            elif stack_plots and output_names[-1] == output_name:
                extra_str = f"{output_name}{f'_vs_{output_name2}' if output_names2 != [] else '_vs_'.join(output_names)}"
            else:
                extra_str = f"{output_name}{f'_vs_{output_name2}' if output_names2 != [] else ''}"

            plt.savefig(
                f"{args.outputdir}/scan_{scan_var_name}_vs_{extra_str}.{save_format}",
                dpi=300,
            )
            if not stack_plots:  # Display plot (used in Jupyter notebooks)
                plt.show()
                plt.clf()
        # ------------

    # In case of a 2D scan
    # ----------------------------------------------------------------------------------------------
    else:
        # Opening the MFILE.DAT
        m_file = mf.MFile(filename=input_files[0])

        # Number of scan points
        n_scan_1 = int(m_file.data["isweep"].get_scan(-1))
        n_scan_2 = int(m_file.data["isweep_2"].get_scan(-1))
        # Selecting the converged runs only
        contour_conv_ij = []  # List of non-converged scan point numbers
        conv_ij = []  # 2D array of converged scan point numbers (sweep = rows, sweep_2 = columns)
        ii_jj = 0
        for ii in range(n_scan_1):
            conv_ij.append([])
            for _jj in range(n_scan_2):
                ii_jj += 1  # Represents the scan point number in the MFILE
                ifail = m_file.data["ifail"].get_scan(ii_jj)
                if ifail == 1:
                    conv_ij[ii].append(
                        ii_jj
                    )  # Only appends scan number if scan converged
                    contour_conv_ij.append(ii_jj)
                else:
                    failed_value_1 = m_file.data[scan_var_name].get_scan(ii_jj)
                    failed_value_2 = m_file.data[scan_2_var_name].get_scan(ii_jj)
                    print(
                        f"Warning : Non-convergent scan point : ({scan_var_name},{scan_2_var_name}) "
                        f"= ({failed_value_1},{failed_value_2})\n"
                        "Warning : This point will not be shown."
                    )

        # Looping over requested outputs
        for index, output_name in enumerate(output_names):
            # Check if the output variable exists in the MFILE
            if output_name not in m_file.data:
                print(
                    f"Warning : `{output_name}` does not exist in PROCESS dicts\n"
                    f"Warning : `{output_name}` will not be output"
                )
                continue

            # Declaring the outputs
            output_arrays = []

            if two_dimensional_contour:
                output_contour_z = np.zeros((n_scan_1, n_scan_2))
                x_contour = [
                    m_file.data[scan_2_var_name].get_scan(i + 1)
                    for i in range(n_scan_2)
                ]
                y_contour = [
                    m_file.data[scan_var_name].get_scan(i + 1)
                    for i in range(1, n_scan_1 * n_scan_2, n_scan_2)
                ]
                for i in contour_conv_ij:
                    output_contour_z[((i - 1) // n_scan_2)][
                        (
                            ((i - 1) % n_scan_2)
                            if ((i - 1) // n_scan_2) % 2 == 0
                            else (-((i - 1) % n_scan_2) - 1)
                        )
                    ] = m_file.data[output_name].get_scan(i)

                flat_output_z = output_contour_z.flatten()
                flat_output_z.sort()
                fig, ax = plt.subplots()
                levels = np.linspace(
                    next(filter(lambda i: i > 0.0, flat_output_z)),
                    flat_output_z.max(),
                    50,
                )
                contour = ax.contourf(
                    x_contour,
                    y_contour,
                    output_contour_z,
                    levels=levels,
                )

                fig.colorbar(contour).set_label(
                    label=(
                        meta[output_name].latex
                        if output_name in meta
                        else f"{output_name}"
                    ),
                    size=axis_font_size,
                )
                plt.ylabel(
                    (
                        meta[scan_var_name].latex
                        if scan_var_name in meta
                        else f"{scan_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.xlabel(
                    (
                        meta[scan_2_var_name].latex
                        if scan_2_var_name in meta
                        else f"{scan_2_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                if y_axis_range != []:
                    y_divisions = (y_axis_range[1] - y_axis_range[0]) / 10
                if y_axis_percentage:
                    if y_max[index] is None:
                        y_max[index] = max(np.abs(y_contour))
                    yticks = mtick.PercentFormatter(y_max[index])
                    if y_axis_range != []:
                        y_divisions = (
                            5 * math.ceil(y_divisions / 5) * y_max[index] / 100
                        )
                        y_range = (
                            y_axis_range[0] * y_max[index] / 100,
                            y_axis_range[1] * y_max[index] / 100,
                        )
                    ax.yaxis.set_major_formatter(yticks)
                if y_axis_range != []:
                    if y_axis_percentage is False:
                        y_range = y_axis_range
                    ax.set_ylim(y_range[0], y_range[1])
                    ax.yaxis.set_major_locator(mtick.MultipleLocator(y_divisions))
                if x_axis_range != []:
                    x_divisions = (x_axis_range[1] - x_axis_range[0]) / 10
                if x_axis_percentage:
                    if x_max[index] is None:
                        x_max[index] = max(np.abs(x_contour))
                    xticks = mtick.PercentFormatter(x_max[index])
                    if x_axis_range != []:
                        x_divisions = (
                            5 * math.ceil(x_divisions / 5) * x_max[index] / 100
                        )
                        x_range = (
                            x_axis_range[0] * x_max[index] / 100,
                            x_axis_range[1] * x_max[index] / 100,
                        )
                    ax.xaxis.set_major_formatter(xticks)
                if x_axis_range != []:
                    if x_axis_percentage is False:
                        x_range = x_axis_range
                    plt.xlim(x_range[0], x_range[1])
                    ax.xaxis.set_major_locator(mtick.MultipleLocator(x_divisions))
                plt.rc("xtick", labelsize=axis_tick_size)
                plt.rc("ytick", labelsize=axis_tick_size)
                plt.tight_layout()
                plt.savefig(
                    f"{args.outputdir}/scan_{output_name}_vs_{scan_var_name}_{scan_2_var_name}.{save_format}"
                )
                plt.grid(True)
                plt.show()
                plt.clf()

            else:
                # Converged indexes, for normal 2D line plot
                fig, ax = plt.subplots()
                for conv_j in (
                    conv_ij
                ):  # conv_j is an array element containing the converged scan numbers
                    # Scanned variables
                    scan_1_var_array = np.zeros(len(conv_j))
                    scan_2_var_array = np.zeros(len(conv_j))
                    output_array = np.zeros(len(conv_j))
                    for jj in range(len(conv_j)):
                        scan_1_var_array[jj] = m_file.data[scan_var_name].get_scan(
                            conv_j[jj]
                        )
                        scan_2_var_array[jj] = m_file.data[scan_2_var_name].get_scan(
                            conv_j[jj]
                        )
                        output_array[jj] = m_file.data[output_name].get_scan(conv_j[jj])

                    # Label formating
                    labl = f"{meta[scan_var_name].latex if scan_var_name in meta else {scan_var_name}} = {scan_1_var_array[0]}"

                    # Plot the graph
                    ax.plot(scan_2_var_array, output_array, "--o", label=labl)

                plt.grid(True)
                plt.ylabel(
                    (
                        meta[output_name].latex
                        if output_name in meta
                        else f"{output_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.xlabel(
                    (
                        meta[scan_2_var_name].latex
                        if scan_2_var_name in meta
                        else f"{scan_2_var_name}"
                    ),
                    fontsize=axis_font_size,
                )
                plt.legend(loc="best", fontsize=legend_size)
                y_data = [
                    m_file.data[output_name].get_scan(i + 1) for i in range(n_scan_2)
                ]
                if y_axis_range != []:
                    y_divisions = (y_axis_range[1] - y_axis_range[0]) / 10
                if y_axis_percentage:
                    if y_max[index] is None:
                        y_max[index] = max(np.abs(y_data))
                    yticks = mtick.PercentFormatter(y_max[index])
                    if y_axis_range != []:
                        y_divisions = (
                            5 * math.ceil(y_divisions / 5) * y_max[index] / 100
                        )
                        y_range = (
                            y_axis_range[0] * y_max[index] / 100,
                            y_axis_range[1] * y_max[index] / 100,
                        )
                    ax.yaxis.set_major_formatter(yticks)
                if y_axis_range != []:
                    if y_axis_percentage is False:
                        y_range = y_axis_range
                    ax.set_ylim(y_range[0], y_range[1])
                    ax.yaxis.set_major_locator(mtick.MultipleLocator(y_divisions))
                x_data = [
                    m_file.data[scan_2_var_name].get_scan(i + 1)
                    for i in range(n_scan_2)
                ]
                if x_axis_range != []:
                    x_divisions = (x_axis_range[1] - x_axis_range[0]) / 10
                if x_axis_percentage:
                    if x_max[index] is None:
                        x_max[index] = max(np.abs(x_data))
                    xticks = mtick.PercentFormatter(x_max[index])
                    if x_axis_range != []:
                        x_divisions = (
                            5 * math.ceil(x_divisions / 5) * x_max[index] / 100
                        )
                        x_range = (
                            x_axis_range[0] * x_max[index] / 100,
                            x_axis_range[1] * x_max[index] / 100,
                        )
                    ax.xaxis.set_major_formatter(xticks)
                if x_axis_range != []:
                    if x_axis_percentage is False:
                        x_range = x_axis_range
                    plt.xlim(x_range[0], x_range[1])
                    ax.xaxis.set_major_locator(mtick.MultipleLocator(x_divisions))
                plt.rc("xtick", labelsize=8)
                plt.rc("ytick", labelsize=8)
                plt.tight_layout()
                plt.savefig(
                    f"{args.outputdir}/scan_{output_name}_vs_{scan_var_name}_{scan_2_var_name}.{save_format}"
                )

                # Display plot (used in Jupyter notebooks)
                plt.show()
                plt.clf()


if __name__ == "__main__":
    main()
