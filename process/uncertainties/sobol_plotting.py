"""
Code to create plots from the output of the Sobols
sensistivity analysis to investigate the input
parameters in PROCESS

Author: A. Pearce (alexander.pearce@ukaea.uk)

Input files:
sobol.txt (datafile output from sobol_method.py,
                             in the same directory as this file)

Output files:
In the work directory specified in the config file
sobol_output.pdf     -  bar chart of sobol indices

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
        the Sobols sensistivity analysis at a given PROCESS design point."
    )

    parser.add_argument(
        "-f",
        "--datafile",
        default="sobol.txt",
        type=str,
        help="datafile for plotting, default = sobol.txt",
    )

    parser.add_argument(
        "-o",
        "--outputfile",
        default="sobol_output.pdf",
        type=str,
        help="filename of outputed pdf file, default = sobol_output.pdf",
    )

    return parser.parse_args(args)


def main(args=None):
    """Read Sobol data and create pdf plot"""

    args = parse_args(args)

    # setput files
    input_raw = args.datafile
    output = args.outputfile
    pdf = bpdf.PdfPages(output)
    page = plt.figure(figsize=(12, 9), dpi=80)

    # read in data
    names = np.loadtxt(input_raw, dtype=str, usecols=[0], skiprows=1)
    s1 = np.loadtxt(input, usecols=[1], skiprows=1)
    s1_conf = np.loadtxt(input_raw, usecols=[2], skiprows=1)
    st = np.loadtxt(input, usecols=[3], skiprows=1)
    st_conf = np.loadtxt(input_raw, usecols=[4], skiprows=1)

    x = np.arange(len(names))
    width = 0.35

    plt.bar(x - width / 2, s1, width, label="$S_1$", yerr=s1_conf)
    plt.bar(x + width / 2, st, width, label="$S_T$", yerr=st_conf)
    plt.xticks(x, names)
    plt.tick_params(labelsize=16)
    plt.legend(fontsize=16)
    plt.ylabel(r"$S_{Sobol}(B_{\$})$", fontsize=22)
    # plt.xlabel('$\mu^{*}$', fontsize=22)

    pdf.savefig(page)
    plt.clf()
    pdf.close()


if __name__ == "__main__":
    main()
