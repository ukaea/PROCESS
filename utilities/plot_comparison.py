#!/usr/bin/env python
"""
Code to display the evolution of two variables in a selection of MFILES

Author: H. Lux (Hanni.Lux@ukaea.uk)

Input files:
List of MFILE.DATs
"""


#######################
# imported libraries

import argparse
from matplotlib import rc
from pylab import figure, xlabel, ylabel, plot, show, savefig
import process.io.mfile as mf

rc("font", size=25)
rc("lines", lw=1.5)  # is overwritten by setting coulours/linestylies
rc(("xtick", "ytick"), labelsize=20)
rc("figure", figsize=(8, 6))
rc("figure.subplot", bottom=0.18)  # was 0.12
rc("figure.subplot", left=0.19)
rc("figure.subplot", right=0.81)


if __name__ == "__main__":

    ###########################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to display\
     the evolution of two variables in a selection of MFILEs."
    )

    PARSER.add_argument(
        "-x", "--xaxis", default="rmajor", help="x-axis, default=rmajor"
    )

    PARSER.add_argument(
        "-y", "--yaxis", default="powfmw", help="y-axis, default=powfmw"
    )

    PARSER.add_argument(
        "mfiles",
        metavar="f",
        type=str,
        default="MFILE.DAT",
        nargs="*",
        help="list of MFiles to be plotted; \
                        default = MFILE.DAT",
    )

    PARSER.add_argument(
        "-e",
        "--end",
        default="pdf",
        help="file format default =\
                        pdf",
    )

    ARGS = PARSER.parse_args()

    XARR = []
    YARR = []

    for mfile in ARGS.mfiles:

        m_file = mf.MFile(mfile)
        scan = -1

        XARR += [m_file.data[ARGS.xaxis.lower()].get_scan(scan)]
        YARR += [m_file.data[ARGS.yaxis.lower()].get_scan(scan)]

    figure()
    xlabel(ARGS.xaxis)
    ylabel(ARGS.yaxis)
    plot(XARR, YARR, "ks")
    savefig(
        "Comparison_"
        + ARGS.xaxis.replace("/", "_")
        + "_"
        + ARGS.yaxis.replace("/", "_")
        + "."
        + ARGS.end
    )

    show()
