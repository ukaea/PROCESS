#!/usr/bin/env python
"""
Code to display the results of an evaluate_uncertainties.py run

Author: H. Lux (Hanni.Lux@ccfe.ac.uk)

Input file:
uncertainties.nc

"""


#######################
# imported libraries

from matplotlib import rc
import argparse
from matplotlib.ticker import NullFormatter
from pylab import figure, axes, show, savefig
from sys import stderr

rc("font", size=20)
rc(("xtick", "ytick"), labelsize=15)


def fig_2dscatter_and_hist(xarr, yarr, labelx, labely):

    """function to create a 2d scatter plot with histograms"""
    nullfmt = NullFormatter()

    figsize = (8, 8)
    left, width = 0.14, 0.61
    bottom, height = 0.14, 0.61
    bottom_h = left_h = left + width

    # start with a rectangular Figure
    figure(figsize=figsize)

    axscatter = axes([left, bottom, width, height])
    axhistx = axes([left, bottom_h, width, 0.2])
    axhisty = axes([left_h, bottom, 0.2, height])

    # no labels
    axhistx.xaxis.set_major_formatter(nullfmt)
    axhistx.yaxis.set_major_formatter(nullfmt)
    axhisty.xaxis.set_major_formatter(nullfmt)
    axhisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axscatter.scatter(xarr, yarr, edgecolors="None")
    axscatter.set_xlabel(labelx)
    axscatter.set_ylabel(labely)

    bins = 10
    axhistx.hist(xarr, bins=bins)
    axhisty.hist(yarr, bins=bins, orientation="horizontal")

    axhistx.set_xlim(axscatter.get_xlim())
    axhisty.set_ylim(axscatter.get_ylim())


if __name__ == "__main__":

    ###########################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to display\
     uncertainties in a given PROCESS design point."
    )

    PARSER.add_argument(
        "-e",
        "--end",
        default="pdf",
        help="file format default =\
pdf",
    )

    PARSER.add_argument(
        "variables",
        metavar="v",
        type=str,
        default="all",
        nargs="*",
        help="list of variables to be plotted; \
default = all",
    )

    PARSER.add_argument(
        "-f",
        "--filename",
        default="uncertainties.nc",
        help="uncertainties data file, default =\
uncertainties.nc",
    )
    ARGS = PARSER.parse_args()

    FILENAME = ARGS.filename

    ############################################################
    # main program

    from process.io.process_netcdf import NetCDFReader

    if ARGS.variables == "all":
        with NetCDFReader(FILENAME) as ncdf_reader:
            DICTS = ncdf_reader.get_run_dicts(0)

            if len(DICTS) > 6:
                print(
                    "Error: There are too many variables stored."
                    "Choose 2 for plotting!",
                    file=stderr,
                )
                print("Stored variables", list(DICTS.keys()), file=stderr)
                exit()

            DATA = {}
            LABELS = list(DICTS.keys())
            for i in range(len(LABELS)):
                DATA[str(i)] = []

            for datadict in ncdf_reader.run_dicts(start_run=0):
                for i in range(len(LABELS)):
                    DATA[str(i)] += [datadict[LABELS[i]]]

        for i in range(len(LABELS) - 1):
            XARR = DATA[str(i)]
            YARR = DATA[str(i + 1)]
            fig_2dscatter_and_hist(XARR, YARR, LABELS[i], LABELS[i + 1])
            savefig("Uncertainties_" + LABELS[i] + "_" + LABELS[i + 1] + "." + ARGS.end)
            print("Stored variables", list(DICTS.keys()))
            print("Number of successful runs", len(XARR))

    else:
        XARR = []
        YARR = []

        with NetCDFReader(FILENAME) as ncdf_reader:

            # Get multiple runs
            for datadict in ncdf_reader.run_dicts(start_run=0):
                try:
                    XARR += [datadict[ARGS.variables[0]]]
                except KeyError:
                    print(
                        "Error: Variable",
                        ARGS.variables[0],
                        "not in",
                        FILENAME,
                        "Choose from",
                        list(datadict.keys()),
                        file=stderr,
                    )
                    exit()
                try:
                    YARR += [datadict[ARGS.variables[1]]]
                except KeyError:
                    print(
                        "Error: Variable",
                        ARGS.variables[1],
                        "not in",
                        FILENAME,
                        "Choose from",
                        list(datadict.keys()),
                        file=stderr,
                    )
                    exit()

        fig_2dscatter_and_hist(XARR, YARR, ARGS.variables[0], ARGS.variables[1])
        savefig(
            "Uncertainties_"
            + ARGS.variables[0]
            + "_"
            + ARGS.variables[1]
            + "."
            + ARGS.end
        )
        print("Stored variables", list(datadict.keys()))
        print("Number of successful runs", len(XARR))


show()
