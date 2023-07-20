#!/usr/bin/env python

"""
Code to produce several plots to aid the diagnosis of a PROCESS run.

Author: H. Lux (Hanni.Lux@ccfe.ac.uk)

Input files:
MFILE.DAT

Output files:


Notes:
19/08/2014 HL initial version of this program


Compatible with PROCESS version 319
"""

#######################
# imported libraries

import argparse
from process_io_lib.diagnose_funcs import plot_normalised_ixc, plot_normalised_icc_res
from pylab import show


if __name__ == "__main__":
    ############################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to diganose\
     a PROCESS run."
    )

    PARSER.add_argument(
        "-f", "--mfile", default="MFILE.DAT", help=" mfile, default = MFILE.DAT"
    )

    ARGS = PARSER.parse_args()

    ############################################################
    # main program

    plot_normalised_ixc(ARGS.mfile)

    plot_normalised_icc_res(ARGS.mfile)

    show()
