#!/usr/bin/env python

"""

  PROCESS make make_plot_dat.out out of MFILE.DAT

  James Morris
  CCFE

  Notes:
    + 20/02/2014: Initial version created
    + 12/03/2014: Updated to put into PyCharm

  Compatible with PROCESS version 274

"""

import os
import argparse
import process.io.mfile as mf
from process.io.mfile import make_plot_dat

# from process_io_lib.process_dicts import PARAMETER_DEFAULTS
from dict_tools import proc_dict


if __name__ == "__main__":

    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Process MFILE.DAT into " "make_plot_dat.out."
    )

    parser.add_argument(
        "-p",
        metavar="p",
        type=str,
        nargs="+",
        help="new variables for the output plot.dat",
    )

    parser.add_argument("-f", metavar="f", type=str, help="File to read as MFILE.DAT")

    parser.add_argument(
        "--defaults", help="run with default params", action="store_true"
    )

    parser.add_argument(
        "--reset-config", help="Reset make_plot_dat.conf", action="store_true"
    )

    parser.add_argument(
        "--columns", help="Write make_plot_dat.out in columns", action="store_true"
    )

    parser.add_argument("--all", help="Output entire MFILE", action="store_true")

    args = parser.parse_args()

    # If user has specified a file that isn't MFILE.DAT pass the filename to
    # MFILE() class.
    if args.f:
        M = mf.MFile(filename=args.f)
    else:
        M.data["rmajor"].get_scan(-1)

    # Get files in current directory to check for the config file.
    current_directory = os.listdir(".")
    if "make_plot_dat.conf" not in current_directory or args.reset_config:
        mf.write_mplot_conf()

    # Read the config file.
    INPUT_CONFIG = mf.read_mplot_conf()

    # If the user added new parameters in the command line add them to the
    # INPUT_CONFIG list to pass to make_plot_dat()
    if args.p:
        conf_file = open("make_plot_dat.conf", "a")
        for item in args.p:
            if item not in INPUT_CONFIG:
                conf_file.write(item + "\n")
        conf_file.close()
        INPUT_CONFIG = mf.read_mplot_conf()

    # If the user has requested the default parameters get default list from
    # process_io_lib
    if args.all:
        make_plot_dat(M, M.data.keys())
    else:
        if args.defaults:

            # If user has specified column format
            if args.columns:
                make_plot_dat(
                    M,
                    proc_dict.get_dicts()[proc_dict.PARAMETER_DEFAULTS],
                    file_format="column",
                )

            # If user has specified row format
            else:
                make_plot_dat(M, proc_dict.get_dicts()[proc_dict.PARAMETER_DEFAULTS])
        else:

            # If user has specified column format
            if args.columns:
                make_plot_dat(M, INPUT_CONFIG, file_format="column")

            # If user has specified row format
            else:
                make_plot_dat(M, INPUT_CONFIG)
