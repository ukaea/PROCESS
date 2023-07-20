#!/usr/bin/env python
"""
Code to display the final uncertain input parameter distributions in comparison
to the specified input distributions

Author: H. Lux (Hanni.Lux@ccfe.ac.uk)

Input file:
uncertainties.nc
evaluate_uncertainties.json

"""


#######################
# imported libraries


from matplotlib import rc
import argparse
from scipy.stats import norm
from numpy import linspace, argwhere, logical_or, zeros, mean
from matplotlib.ticker import NullFormatter
from pylab import figure, hist, show, savefig, gca, plot, xlabel
from sys import stderr
from process_io_lib.process_config import UncertaintiesConfig
from process_io_lib.process_dicts import DICT_INPUT_BOUNDS

rc("font", size=25)
rc("lines", lw=1.5)  # is overwritten by setting coulours/linestylies
rc(("xtick", "ytick"), labelsize=20)
rc("figure", figsize=(8, 6))
rc("figure.subplot", bottom=0.18)  # was 0.12
rc("figure.subplot", left=0.19)
rc("figure.subplot", right=0.81)


def plot_distribution(xarr, labelx, unc_dict):

    """function to create a 2d scatter plot with histograms"""

    figure()

    # no labels
    gca().yaxis.set_major_formatter(NullFormatter())
    xlabel(labelx)

    (n_arr, bins, patches) = hist(xarr, bins=10, normed=True)

    if unc_dict["errortype"].lower() == "gaussian":
        u_mean = unc_dict["mean"]
        u_std = unc_dict["std"]
        xvalues = linspace(min(xarr), max(xarr), 500)
        yvalues = norm.pdf(xvalues, u_mean, u_std)
        # assures values are inside input bounds!
        if varname in DICT_INPUT_BOUNDS:
            args = argwhere(
                logical_or(
                    xvalues < DICT_INPUT_BOUNDS[varname]["lb"],
                    xvalues > DICT_INPUT_BOUNDS[varname]["ub"],
                )
            )
            yvalues[args] = zeros(args.shape)

        else:  # cutoff at 0 - typically negative values are meaningless
            args = argwhere(xvalues < 0.0)
            yvalues[args] = zeros(args.shape)

    elif unc_dict["errortype"].lower() == "uniform":
        lbound = unc_dict["lowerbound"]
        ubound = unc_dict["upperbound"]
        xvalues = linspace(lbound, ubound, 10)
        yvalues = [mean(n_arr)] * 10

    elif unc_dict["errortype"].lower() == "relative":
        err = unc_dict["percentage"] / 100.0
        lbound = unc_dict["mean"] * (1.0 - err)
        ubound = unc_dict["mean"] * (1.0 + err)
        xvalues = linspace(lbound, ubound, 10)
        yvalues = [mean(n_arr)] * 10

    elif unc_dict["errortype"].lower() == "lowerhalfgaussian":
        u_mean = unc_dict["mean"]
        u_std = unc_dict["std"]
        xvalues = linspace(min(xarr), max(xarr), 500)
        yvalues = norm.pdf(xvalues, u_mean, u_std)
        # to correct normalisation for half Gaussian
        yvalues = yvalues * 2.0
        if varname in DICT_INPUT_BOUNDS:
            args = argwhere(
                logical_or(xvalues < DICT_INPUT_BOUNDS[varname]["lb"], xvalues > u_mean)
            )
            yvalues[args] = zeros(args.shape)

        else:
            args = argwhere(logical_or(xvalues < 0.0, xvalues > u_mean))
            yvalues[args] = zeros(args.shape)

    elif unc_dict["errortype"].lower() == "upperhalfgaussian":
        u_mean = unc_dict["mean"]
        u_std = unc_dict["std"]
        xvalues = linspace(min(xarr), max(xarr), 500)
        yvalues = norm.pdf(xvalues, u_mean, u_std)
        # to correct normalisation for half Gaussian
        yvalues = yvalues * 2.0
        if varname in DICT_INPUT_BOUNDS:
            args = argwhere(
                logical_or(xvalues < u_mean, xvalues > DICT_INPUT_BOUNDS[varname]["ub"])
            )
            yvalues[args] = zeros(args.shape)
        else:
            args = argwhere(xvalues < u_mean)
            yvalues[args] = zeros(args.shape)

    plot(xvalues, yvalues, "k-")


if __name__ == "__main__":

    ###########################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to check\
     the final uncertainty distributions in the input parameters."
    )

    PARSER.add_argument(
        "-e",
        "--end",
        default="pdf",
        help="file format default =\
pdf",
    )

    PARSER.add_argument(
        "-u",
        "--uncertainties",
        default="evaluate_uncertainties.json",
        help="uncertainties config file \
default = evaluate_uncertainties.json",
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

    from process_io_lib.process_netcdf import NetCDFReader

    CONFIG = UncertaintiesConfig(ARGS.uncertainties)

    if CONFIG.uncertainties != []:

        with NetCDFReader(FILENAME) as ncdf_reader:

            for u_dict in CONFIG.uncertainties:
                varname = u_dict["varname"].lower()

                XARR = []
                for datadict in ncdf_reader.run_dicts(start_run=0):
                    try:
                        XARR += [datadict[varname]]
                    except KeyError:
                        if varname == "fimp(14)":
                            try:
                                XARR += [datadict["fimp(14"]]
                            except KeyError:
                                print(
                                    "Error: Variable",
                                    varname,
                                    "can currently not be treated!\n",
                                    "Check separately! \n",
                                    list(datadict.keys()),
                                    file=stderr,
                                )
                                break
                        else:
                            print(
                                "Error: Variable",
                                varname,
                                "can currently not be treated!\n",
                                "Check separately! \n",
                                list(datadict.keys()),
                                file=stderr,
                            )
                            break

                if XARR != []:
                    plot_distribution(XARR, varname, u_dict)
                    savefig("Uncertainties_Diagnostic_" + varname + "." + ARGS.end)
            print("Number of successful runs", len(XARR))


show()
