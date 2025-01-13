"""
script for reading and display HDF5 files as pdf scatter plots

Alexander J Pearce
05/03/22
CCFE

"""

import argparse

import pandas as pd
from pylab import figure, savefig


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Program to read and  plot PROCESS hdf5 output."
    )

    parser.add_argument(
        "-i",
        "--inputfile",
        default="uncertainties_data.h5",
        help="input hdf5 file (default=uncertainties_data.h5)",
    )

    parser.add_argument(
        "-v",
        "--vars",
        default="None",
        nargs="*",
        help="Select the output variables (default = None) \n More than one output can be plotted eg: -v 'var1 var2'\n A separate plot will be created in matrix for each inputs combination",
    )

    parser.add_argument(
        "-sf",
        "--save_format",
        nargs="?",
        default="pdf",
        help="Output format (default='pdf') ",
    )

    parser.add_argument(
        "-p",
        "--print",
        action="store_true",
        help="(print the PROCESS data to comannd line. default=False)",
    )

    return parser.parse_args(args)


def main(args=None):
    """reads inputfile and creats figure as outputfile

    :param args: None
    :return: None
    """
    args = parse_args(args)

    data_set = pd.read_hdf(args.inputfile)
    if args.print:
        print(data_set)

    output_names = args.vars

    # Select only the converged runs for creating KDF plots
    # see if we can only use if ifail = 1 for the KDF

    # TO DO make separate list of converged and failed runs
    # find way to display both to look at margins
    # we want to colour red for ifail = 2-6 and blue for ifail = 1
    ioptimz = data_set["ioptimz"][0]
    if ioptimz == -2:
        data_set_converge = data_set
    else:
        data_set_converge = data_set[data_set["ifail"] == 1.0]

    figsize = (8, 8)
    figure(figsize=figsize)
    pd.plotting.scatter_matrix(
        data_set_converge[output_names], alpha=0.2, diagonal="kde"
    )
    savefig("uncertainties." + args.save_format)


if __name__ == "__main__":
    main()
