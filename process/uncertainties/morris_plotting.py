#!/usr/bin/env python
"""
Code to create plots from the output of the Morris method
of elementary elements to investiage the sensistivity of
the input parameters in PROCESS

Author: A. Pearce (alexander.pearce@ukaea.uk)

Input files:
morris_method_output.txt (datafile output from morris_method.py,
                             in the same directory as this file)

Output files:
In the work directory specified in the config file
morris_output.pdf     -  scatter plot of mean and variance of
                        morris method output

"""

import argparse

import matplotlib.backends.backend_pdf as bpdf
import matplotlib.pyplot as plt
import numpy as np


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Program to plot the output of the\
        the sensistivity analysis by elementary element method at a given PROCESS design point."
    )

    parser.add_argument(
        "-f",
        "--datafile",
        default="morris_method_output.txt",
        type=str,
        help="datafile for plotting, default = morris_method_output.txt",
    )

    parser.add_argument(
        "-o",
        "--outputfile",
        default="morris_output.pdf",
        type=str,
        help="filename of outputed pdf file, default = morris_output.pdf",
    )

    return parser.parse_args(args)


def main(args=None):
    """Read Morris data and create pdf plot"""
    args = parse_args(args)

    # setput files
    input = args.datafile
    output = args.outputfile
    pdf = bpdf.PdfPages(output)
    page = plt.figure(figsize=(12, 9), dpi=80)

    # read in data
    n = np.loadtxt(input, dtype=str, usecols=[0], skiprows=1)
    z = np.loadtxt(input, usecols=[2], skiprows=1)
    y = np.loadtxt(input, usecols=[3], skiprows=1)

    plt.scatter(z, y)

    for i, txt in enumerate(n):
        plt.annotate(txt, (z[i], y[i]), fontsize=16)

    plt.ylabel(r"$\sigma$", fontsize=22)
    plt.xlabel(r"$\mu^{*}$", fontsize=22)
    plt.tick_params(labelsize=20)

    pdf.savefig(page)
    plt.clf()
    pdf.close()


if __name__ == "__main__":
    main()
