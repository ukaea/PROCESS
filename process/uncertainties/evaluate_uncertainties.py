#!/usr/bin/env python
"""
Code to assess uncertainties in the input parameters in PROCESS. Allows
for use with monte carlo, Morris method and Sobol techniques.

Author: H. Lux (Hanni.Lux@ukaea.uk)
        A.J. Pearce (Alexander.Pearce@ukaea.uk)

Input files:
config_evaluate_uncertainties.json (config file for uncertainties in same
                             directory as this file)
An IN.DAT file as specified in the config file

Output files:
OUT.DAT               - PROCESS output
MFILE.DAT             - PROCESS output
uncertainties_data.h5 - Dataframe of all PROCESS MFILE Outputs
morris_method.txt     - contains the dictorany of the mean and
                        variance of the elementary effects
                        (only if method = morris_method)
sobol.txt             - contains the dictionary of the Sobol
                        indices (only if method = sobol_method)

"""

import argparse

from SALib.analyze import sobol
from SALib.analyze import morris as morris_method
from SALib.sample import morris, saltelli

from pathlib import Path
import pandas as pd
import numpy as np
import process.io.mfile as mf
from process.io.in_dat import InDat
from process.io.process_config import UncertaintiesConfig
from process.io.process_funcs import (
    get_neqns_itervars,
    get_variable_range,
    check_input_error,
    process_stopped,
    no_unfeasible_mfile,
    set_variable_in_indat,
    vary_iteration_variables,
)


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Program to evaluate "
        "uncertainties in a given PROCESS design point."
    )

    parser.add_argument(
        "-f",
        "--configfile",
        default="config_evaluate_uncertainties.json",
        help="configuration file",
    )

    parser.add_argument(
        "-m",
        "--method",
        default="monte_carlo",
        help="Type of uncertainty analysis performed (default=monte_carlo)",
    )

    return parser.parse_args(args)


def run_monte_carlo(args):
    """Runs Monte Carlo uncertainty analysis

    :param args: None
    :type args: None
    :return: MFileDataSet, indexDataSet
    :rtype: list
    """

    config = UncertaintiesConfig(args.configfile)
    config.setup()

    NEQNS, itervars = get_neqns_itervars()

    config.factor = 1.0

    LBS, UBS = get_variable_range(itervars, config.factor)

    config.checks_before_run()
    config.set_sample_values()
    RUN_ID = 0
    MFileDataSet = []
    indexDataSet = []
    outputDataSet = []

    for j in range(config.no_samples):

        print("sample point", j, ":")
        config.go2newsamplepoint(j)

        for i in range(config.niter):

            print("  ", i, end=" ")

            # Define path for the input file to run
            input_path = Path(config.wdir) / "IN.DAT"
            config.run_process(input_path)

            # Check for produced MFILE.DAT in working dir
            check_input_error(wdir=config.wdir)

            if not process_stopped():

                no_unfeasible = no_unfeasible_mfile()

                if no_unfeasible <= config.no_allowed_unfeasible:
                    if no_unfeasible > 0:
                        print(
                            "WARNING: %i non feasible point(s) in sweep,"
                            % no_unfeasible,
                            "but finished anyway! Allowed  %i. "
                            % config.no_allowed_unfeasible,
                        )

                    # this might need to be depraitedd
                    # CONFIG.add_results2netcdf(RUN_ID)

                    # add run to index list
                    indexDataSet.append("run{}".format(RUN_ID))

                    # collect the process solution
                    mfilepath = Path(config.wdir) / "MFILE.DAT"
                    m_file = mf.MFile(mfilepath)
                    for item in m_file.data.keys():
                        outputDataSet.append(m_file.data[item].get_scan(-1))

                    # Append process data to list of all mfile data
                    MFileDataSet.append(outputDataSet)
                    outputDataSet = []

                    RUN_ID += 1

                else:
                    print(
                        "WARNING: %i non feasible point(s) in sweep!\
                    Rerunning!"
                        % no_unfeasible
                    )
            else:
                print("PROCESS has stopped without finishing!")

            if config.vary_iteration_variables is True:
                vary_iteration_variables(itervars, LBS, UBS)

        config.write_error_summary(j)

    # create list of variables used in MFILE
    columnDataSet = column_data_list(config.wdir)

    return MFileDataSet, indexDataSet, columnDataSet


def write_Morris_Method_Output(X, S):
    """Writes Morris Method output to .txt file
    :param args: X: morris method input bounds
    :param args: S: morris method output
    :type args: dict
    :return: None
    """

    with open("morris_method.txt", "w") as f:
        # create sensistivity indices header
        f.write("Parameter mu mu_star sigma mu_star_conf\n")
        # print the sensistivity indices
        for i in range(X["num_vars"]):
            f.write(
                "%s %f %f %f %f\n"
                % (
                    X["names"][i],
                    S["mu"][i],
                    S["mu_star"][i],
                    S["sigma"][i],
                    S["mu_star_conf"][i],
                )
            )


def write_Sobol_Output(X, S):
    """Writes Sobol Method output to .txt file
    :param args: X: sobol method input bounds
    :param args: S: sobol method output
    :type args: dict
    :return: None
    """

    with open("sobol.txt", "w") as f:
        # create first order header
        f.write("Parameter S1 S1_conf ST ST_conf\n")
        # print first order Sobol indices
        for i in range(X["num_vars"]):
            f.write(
                "%s %f %f %f %f\n"
                % (
                    X["names"][i],
                    S["S1"][i],
                    S["S1_conf"][i],
                    S["ST"][i],
                    S["ST_conf"][i],
                )
            )


def column_data_list(working_dir):
    """Reads the list outputs in MFILE.DAT file
    :param args: args
    :return: columnDataSet
    :rtype: list
    """
    # CONFIG = UncertaintiesConfig(args.configfile)
    # CONFIG.setup()

    # find path to MFILE data
    mfilepath = Path(working_dir) / "MFILE.DAT"
    m_file = mf.MFile(mfilepath)

    # read data in MFILE
    columnDataSet = []
    for item in m_file.data.keys():
        columnDataSet.append(item)

    return columnDataSet


def run_morris_method(args):
    """Runs Morris method uncertainty analysis

    :param args: None
    :return: MFileDataSet, indexDataSet
    :rtype: list
    """

    config = UncertaintiesConfig(args.configfile)
    config.setup()

    sols = np.array([])
    fail = np.array([])
    MFileDataSet = []
    indexDataSet = []
    outputDataSet = []

    # Set up the sample needed using latin hypercube scaling
    params_values = morris.sample(
        config.morris_uncertainties,
        config.no_samples,
        num_levels=config.latin_hypercube_level,
    )

    np.savetxt("param_values.txt", params_values)

    in_dat = InDat()

    run_max = int((config.morris_uncertainties["num_vars"] + 1.0) * config.no_samples)
    for run_id in range(run_max):
        print("run number =", run_id)
        i = 0
        for n in config.morris_uncertainties["names"]:
            set_variable_in_indat(in_dat, n, params_values[run_id][i])
            i = i + 1

        in_dat.write_in_dat(output_filename="IN.DAT")  # this is from in_dat

        # CONFIG.run_process()
        # Define path for the input file to run
        input_path = Path(config.wdir) / "IN.DAT"
        config.run_process(
            input_path
        )  # run_process in uncertaintyconfig or RunProcessConfig

        mfilepath = Path(config.wdir) / "MFILE.DAT"
        m_file = mf.MFile(mfilepath)
        for item in m_file.data.keys():
            outputDataSet.append(m_file.data[item].get_scan(-1))

        # Append process data to list of all mfile data
        MFileDataSet.append(outputDataSet)
        outputDataSet = []

        # add run to index list
        indexDataSet.append("run{}".format(run_id))

        # We need to find a way to catch failed runs
        process_status = m_file.data["ifail"].get_scan(-1)
        print("ifail =", process_status)

        if process_status == 1.0:
            # read the figure of merit from the MFILE
            FoM = m_file.data[config.figure_of_merit].get_scan(-1)
            sols = np.append(sols, FoM)
        else:
            # if run doesn't converge set FoM to previous value
            fail = np.append(fail, run_id)
            if run_id == 0:
                # if first run doesn't converge set FoM to mean
                sols = np.append(sols, config.output_mean)
            else:
                sols = np.append(sols, sols[np.size(sols) - 1])

    params_sol = morris_method.analyze(config.morris_uncertainties, params_values, sols)
    # np.savetxt("FoM_sol.txt", sols)
    # np.savetxt("error_log.txt", fail)

    # write output file
    write_Morris_Method_Output(config.morris_uncertainties, params_sol)

    # create list of variables used in MFILE
    columnDataSet = column_data_list(config.wdir)

    return MFileDataSet, indexDataSet, columnDataSet


def run_sobol_method(args):
    """Runs Sobol method uncertainty analysis
    :param args: None
    :return: MFileDataSet, indexDataSet
    :rtype: list
    """

    config = UncertaintiesConfig(args.configfile)
    config.setup()

    # Setup output arrays
    sols = np.array([])
    fail = np.array([])
    MFileDataSet = []
    indexDataSet = []
    outputDataSet = []

    # Generate samples
    params_values = saltelli.sample(
        config.sobol_uncertainties, config.no_samples, calc_second_order=False
    )

    np.savetxt("param_values.txt", params_values)

    in_dat = InDat()

    # find capcost_sols from PROCESS
    run_max = int((config.sobol_uncertainties["num_vars"] + 2.0) * config.no_samples)
    for run_id in range(run_max):
        print("run number =", run_id)
        i = 0
        for n in config.sobol_uncertainties["names"]:
            set_variable_in_indat(in_dat, n, params_values[run_id][i])
            i = i + 1

        in_dat.write_in_dat(output_filename="IN.DAT")

        # Define path for the input file to run
        input_path = Path(config.wdir) / "IN.DAT"
        config.run_process(input_path)

        mfilepath = Path(config.wdir) / "MFILE.DAT"
        m_file = mf.MFile(mfilepath)
        for item in m_file.data.keys():
            outputDataSet.append(m_file.data[item].get_scan(-1))

        # Append process data to list of all mfile data
        MFileDataSet.append(outputDataSet)
        outputDataSet = []

        # add run to index list
        indexDataSet.append("run{}".format(run_id))

        # We need to find a way to catch failed runs
        process_status = m_file.data["ifail"].get_scan(-1)
        print("ifail =", process_status)

        if process_status == 1.0:
            capcost = m_file.data[config.figure_of_merit].get_scan(-1)
            sols = np.append(sols, capcost)
        else:
            # if run doesn't converge use mean FoM value
            fail = np.append(fail, run_id)
            sols = np.append(sols, config.output_mean)

    params_sol = sobol.analyze(
        config.sobol_uncertainties, sols, calc_second_order=False, print_to_console=True
    )

    # write output file
    write_Sobol_Output(config.sobol_uncertainties, params_sol)

    # create list of variables used in MFILE
    columnDataSet = column_data_list(config.wdir)

    return MFileDataSet, indexDataSet, columnDataSet


def main(args=None):
    args = parse_args(args)

    if args.method == "monte_carlo":
        MFileDataSet, indexDataSet, columnDataSet = run_monte_carlo(args)
    elif args.method == "morris_method":
        MFileDataSet, indexDataSet, columnDataSet = run_morris_method(args)
    elif args.method == "sobol_method":
        MFileDataSet, indexDataSet, columnDataSet = run_sobol_method(args)
    else:
        print("Uncertainty method not recognised!")

    # create list of variables used in MFILE
    # columnDataSet = column_data_list(args)

    df = pd.DataFrame(data=MFileDataSet, columns=columnDataSet, index=indexDataSet)
    df.to_hdf("uncertainties_data.h5", key="df", mode="w")
    print("UQ finished!")


if __name__ == "__main__":
    main()
