#!/usr/bin/env python
"""
Code to display the power flow of a PROCESS run in a Sankey diagram

Author: H. Lux (Hanni.Lux@ukaea.uk)
Updated 20/08/219: A. Brown (adam.brown@ukaea.uk)

Input file:
MFILE.DAT

"""
import matplotlib
import argparse
from pylab import show, savefig
from process.io.sankey_funcs import plot_sankey

matplotlib.use("Agg")


def main(args=None):
    ###########################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to plot\
     the power flow in PROCESS using a Sankey diagram."
    )

    PARSER.add_argument("-e", "--end", default="pdf", help="file format, default = pdf")

    PARSER.add_argument(
        "-m", "--mfile", default="MFILE.DAT", help="mfile name, default = MFILE.DAT"
    )

    ARGS = PARSER.parse_args(args)

    #########################################################
    # main program

    plot_sankey(ARGS.mfile)
    savefig("SankeyPowerFlow." + ARGS.end)

    show()


if __name__ == "__main__":
    main()
