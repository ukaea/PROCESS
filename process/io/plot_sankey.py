"""
Code to display the power flow of a PROCESS run in a Sankey diagram

Author: H. Lux (Hanni.Lux@ukaea.uk)
Updated 20/08/219: A. Brown (adam.brown@ukaea.uk)

Input file:
MFILE.DAT

"""

import argparse

import matplotlib as mpl
from pylab import savefig, show

from process.io.sankey_funcs import plot_sankey

mpl.use("Agg")


def main(args=None):
    ###########################################################
    # Usage

    parser = argparse.ArgumentParser(
        description="Program to plot\
     the power flow in PROCESS using a Sankey diagram."
    )

    parser.add_argument("-e", "--end", default="pdf", help="file format, default = pdf")

    parser.add_argument(
        "-m", "--mfile", default="MFILE.DAT", help="mfile name, default = MFILE.DAT"
    )

    args = parser.parse_args(args)

    #########################################################
    # main program

    plot_sankey(args.mfile)
    savefig("SankeyPowerFlow." + args.end)

    show()


if __name__ == "__main__":
    main()
