#!/usr/bin/env python
""" Set of macro that plots information about the VMNCON optimization
    WARNING : YOU HAVE TO CREATE THE PROCESS DICT IN ORDER TO USE THIS MACRO !!


    PART 1 : The file OPT.DAT is read
    _________________________________
    This file contains :
     - The PROCESS id numers of
        1) The constraints : constraints_indexes
        2) The variables   : variables_indexes
     - The evolution with each VMCON interation of
                Description                                            | variable name in vmcon | variables name in the macro
        1) The figure of merit                                         |   abs(obj)             |   figures_of_merit     []
        2) The VMCON convergence criteria                              |   sum                  |   vmcon_convergence    []
        3) The squared sum of the constraints residuals                |   sqsumsq              |   constraints_quad_sum []
        4) The residuals of each constraints                           |   conf[]               |   constraints         [[]]
        5) The variables values normalized to first the interation one |   x[]                  |   variables           [[]]
    _________________________________


    PART 2 : Plot
    _____________
    The pyplot routines

    Plot 1 : Figure of merit plot
       VMCON index evolution of the figure of merit

    Plot 2 : Convergence plot :
       VMCON index evolution of
        - The VMCON convergence parameter
        - The quadratic sum of the constraints residuals
        - The maximum between the VMCON convergence parameter and the constraints residual quadratic sum
          This quantity is the one USED to decide that optimizer has converged or not

    Plot 3 : Dominant constraints plots :
       The last VMCON iteration is used order to rank the constraints with their residual values and to plot the
       VMCON index evolution of
        - The n_ploted_const dominant constraints values (n_ploted_const is used defined)
        - The quadratic sum of the dominant constraints
        - The total quadratic sum of the constraints (sqsumsq)
       The difference of the two quadratic sums allow to check of any other variables contribute to the constaints for any step of the optimization

    Plot 4 : Selected constraint plot :
       VMCON index evolution of
        - A constraint selected by its PROCESS number defined in vardes
        - The total quadratic sum of the constraints (sqsumsq)
       Help to focus on a given constraint ...

    Plot 5 : Major variable evolution
       The variation amplitude of the optimization varaibles (max(var) - min(var)) is used to rank the variables and to plot the
       VMCON index evolution of the n_var_plots dominant variables (n_var_plots is used defined)

    Plot 6 : Selected vaiable pair trajectory plot
       The x and y variables are selected using the PROCESS variable index defined in vardes
       The z axis also shows ln( max( sum, sqsumsq ) )

    The plots are saved in :
     - a used defined format (it has to be supported by pyplot)
     - a folder called OPT_plots folder
     - The output file naming convention still has to be discussed
    _____________
 """
import matplotlib
import os
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
from create_dicts import get_dicts

matplotlib.use("Agg")


# Load dicts from dicts JSON file
p_dicts = get_dicts()


if __name__ == "__main__":

    # PARSING USER PARAMETERS
    # please execute 'python plot_opti_process.py -h' for input information
    # Option definition
    # -----------------
    parser = argparse.ArgumentParser(
        description="Plot optimization information",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-p",
        "--plot_selec",
        nargs="?",
        default="all",
        help="Plot selection string :\n - If it containts 'FoM'      -> Figure of Merit plot \n - If it containts 'conv'     -> convergence criteria plot \n - If it containts 'domconst' -> dominant constraint plot\n - If it containts 'allconst' -> a plot for each constraint stored in All_Const/\n - If it containts 'domvar'   -> dominant variables plot \n - If it containts 'allvar'   -> a plot for each variable in All_Var/\n - If it containts 'all'      -> all the mentionned plots (default value)",
    )
    parser.add_argument(
        "-ndc",
        "--n_dom_const",
        nargs="?",
        default=3,
        help="number of plotted dominant constaints (default=3)",
        type=int,
    )
    parser.add_argument(
        "-ndv",
        "--n_dom_var",
        nargs="?",
        default=4,
        help="number of plotted dominant variables  (default=4)",
        type=int,
    )
    parser.add_argument(
        "-ic",
        "--i_const",
        nargs="?",
        default=-1,
        help="Selection of the constraint to be plotted (PROCESS number defined in vardes, default=-1)",
        type=int,
    )
    parser.add_argument(
        "-ixv",
        "--i_X_var",
        nargs="?",
        default=-1,
        help="X variable on pair plot selection (PROCESS number defined in vardes, default=-1)",
        type=int,
    )
    parser.add_argument(
        "-iyv",
        "--i_Y_var",
        nargs="?",
        default=-1,
        help="Y variable on pair plot selection (PROCESS number defined in vardes, default=-1)",
        type=int,
    )
    parser.add_argument(
        "-sf",
        "--save_format",
        nargs="?",
        default="pdf",
        help="output format (default='pdf') ",
    )
    parser.add_argument(
        "-as",
        "--axis_font_size",
        nargs="?",
        default=18,
        help="Axis label font size selection (default=18)",
        type=int,
    )

    # Option argument extraction
    # --------------------------
    args = parser.parse_args()
    plot_selection = str(args.plot_selec)
    n_ploted_const = int(args.n_dom_const)
    n_var_plots = int(args.n_dom_var)
    constraint_selection = int(args.i_const)
    x_variable_selection = int(args.i_X_var)
    y_variable_selection = int(args.i_Y_var)
    save_format = str(args.save_format)

    # Hardcoded input : Plot option selection (color and line style)
    plot_opt = [
        "r-",
        "b-",
        "g-",
        "m-",
        "r--",
        "b--",
        "g--",
        "m--",
        "r:",
        "b:",
        "g:",
        "m:",
    ]

    # Boolean swiches for plot selection
    # -----------------------------------
    plot_FoM = ("FoM" in plot_selection) or "all" == plot_selection
    plot_conv = ("conv" in plot_selection) or "all" == plot_selection
    plot_const_dom = ("domconst" in plot_selection) or "all" == plot_selection
    plot_var_dom = ("domvar" in plot_selection) or "all" == plot_selection
    plot_var_all = ("allvar" in plot_selection) or "all" == plot_selection
    plot_const_all = ("allconst" in plot_selection) or "all" == plot_selection

    plot_const_spe = int(-1) != constraint_selection
    plot_var_spe_pair = (
        int(-1) != x_variable_selection and int(-1) != y_variable_selection
    )
    axis_font_size = int(args.axis_font_size)
    #####################################################

    # Step 1 : Data extraction
    # ----------------------------------------------------------------------------------------------
    plot_indexes = list()
    vmcon_indexes = list()
    figures_of_merit = list()
    vmcon_convergence = list()
    constraints_quad_sum = list()
    convergence_parameter = list()
    constraints = list()
    variables = list()
    gradients = list()

    # Opening the pandora box
    with open("../bin/OPT.DAT", "r") as opt_data:
        opt_data_lines = opt_data.readlines()
        opt_data_lines = [
            line.strip("\n") for line in opt_data_lines
        ]  # Just removing the \n statment

        # Optimization setup information extraction
        # ------------------------------------------
        # Constrains
        n_constraints = int(opt_data_lines[1])
        constraints_indexes = [
            int(constraints_str)
            for constraints_str in opt_data_lines[4]
            .replace("  ", " ")
            .replace("  ", " ")[1:]
            .split(" ")
        ]

        # Variables
        n_variables = int(opt_data_lines[7])
        variables_indexes = [
            int(constraints_str)
            for constraints_str in opt_data_lines[10]
            .replace("  ", " ")
            .replace("  ", " ")[1:]
            .split(" ")
        ]

        # Number of VMCON iterations
        n_vmcon = int(len(opt_data_lines) - 15)

        # Plot indexes
        plot_indexes = [ii + 1 for ii in range(n_vmcon)]
        # ------------------------------------------

        # VMCON data extraction
        # ----------------------
        # Full data extraction
        data = list()
        data_str = list()
        for ii_vmcon in range(0, n_vmcon):
            data_str.append(opt_data_lines[15 + ii_vmcon].split(" "))
            ii = int(0)
            while ii != len(data_str[ii_vmcon]):
                if data_str[ii_vmcon][ii] == "":
                    data_str[ii_vmcon].pop(ii)
                else:
                    ii += 1

            data.append(list())
            for data_ii in data_str[ii_vmcon]:
                data[ii_vmcon].append(float(data_ii))

        # VMCON global indexes
        for ii_vmcon in range(0, n_vmcon):
            vmcon_indexes.append(int(data[ii_vmcon][0]))
            figures_of_merit.append(float(data[ii_vmcon][1]))
            vmcon_convergence.append(float(data[ii_vmcon][2]))
            constraints_quad_sum.append(float(data[ii_vmcon][3]))
            convergence_parameter.append(
                max(vmcon_convergence[ii_vmcon], constraints_quad_sum[ii_vmcon])
            )

        # Constrains
        for ii_constraint in range(0, n_constraints):
            constraints.append(list())
            for ii_vmcon in range(0, n_vmcon):
                constraints[ii_constraint].append(
                    abs(float(data[ii_vmcon][ii_constraint + 4]))
                )

        # Variables
        for ii_variables in range(0, n_variables):
            variables.append(list())
            gradients.append(list())
            for ii_vmcon in range(0, n_vmcon):
                variables[ii_variables].append(
                    float(data[ii_vmcon][4 + n_constraints + ii_variables])
                )
                gradients[ii_variables].append(
                    float(
                        data[ii_vmcon][4 + n_constraints + n_variables + ii_variables]
                    )
                )

        variables = variables[:-1]
        gradients = gradients[:-1]
        constraints = constraints[:-1]
        # ---------------------

    # Step 2 : Plotting
    # ----------------------------------------------------------------------------------------------
    tolerance = float(1e-7)
    if not os.path.isdir("OPT_plots"):
        os.mkdir("OPT_plots")

    # PLOT 1 : Figure of merit evolution
    # -----------------------------------
    if plot_FoM:

        # Plot
        plt.plot(plot_indexes, figures_of_merit, "g-")

        # Range
        max_FoM = max(figures_of_merit)
        min_FoM = min(figures_of_merit)
        y_range_FoM = max_FoM - min_FoM
        y_min = min_FoM - 0.1 * y_range_FoM
        y_max = max_FoM + 0.1 * y_range_FoM
        x_min = 0
        x_max = n_vmcon + 2

        # Cosmetics
        plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
        plt.ylabel("Figure of merit", fontsize=axis_font_size)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("OPT_plots/FoM_evolution." + save_format, format=save_format)
        plt.close()
    # -----------------------------------

    # Plot 2 : Convergence parameters
    # --------------------------------
    if plot_conv:

        # Ranges
        x_min = 0
        x_max = n_vmcon + 2
        y_min = 1e-8
        y_max = 1.0

        # Plot
        plt.plot(plot_indexes, vmcon_convergence, "b-", label="VMCON conv")
        plt.plot(plot_indexes, constraints_quad_sum, "c-", label="const quad sum")
        plt.plot(plot_indexes, convergence_parameter, "k:", label="final conv    ")
        plt.plot([x_min, x_max], [tolerance, tolerance], "k-")

        # Cosmetics
        plt.legend(loc="best")
        plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
        plt.ylabel("Criteria", fontsize=axis_font_size)
        plt.yscale("log")
        plt.axis([x_min, x_max, y_min, y_max])
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(
            "OPT_plots/convergence_evolution." + save_format, format=save_format
        )
        plt.close()
    # --------------------------------

    # Plots 3 : Domimant constraints parameters
    # ------------------------------------------
    # This graph will show the n_ploted_const constrains that dominates sqsumsq at the last iteration, their quadratic sum and sqsumsq
    # ---
    if n_ploted_const > n_constraints:
        print(
            "Impossible to print {} constraints as only {} were used for the run".format(
                n_ploted_const, n_constraints
            )
        )
        print(" -> Plotting {} constraints instead".format(n_constraints))
        n_ploted_const = n_constraints

    if plot_const_dom:
        # Avoiding trying to plot more constraints than the total number of constraints used for the run
        n_ploted_const = min(n_ploted_const, n_constraints)

        # Finiding n_ploted_const constraints that has the largest value at the last point
        dominant_constaints_indexes = []
        last_constraints = data[-1][4 : n_constraints + 4]
        last_constraints = [
            abs(ii_last_constraint) for ii_last_constraint in last_constraints
        ]

        for ii in range(0, n_ploted_const):
            ii_dominant_constaints = int(last_constraints.index(max(last_constraints)))
            dominant_constaints_indexes.append(ii_dominant_constaints)
            last_constraints[ii_dominant_constaints] = -1.0

        # Dominant constraints quadratic sum calculation
        dominant_quad_sum = []
        for ii_vmcon in range(0, n_vmcon):
            dominant_quad_sum.append(float(0))
            for ii in range(0, n_ploted_const):
                dominant_quad_sum[ii_vmcon] += (
                    constraints[dominant_constaints_indexes[ii]][ii_vmcon]
                    * constraints[dominant_constaints_indexes[ii]][ii_vmcon]
                )
            dominant_quad_sum[ii_vmcon] = pow(dominant_quad_sum[ii_vmcon], 0.5)

        # Plot option settings
        if n_ploted_const > len(plot_opt):
            for ii in range(0, n_ploted_const - len(plot_opt)):
                plot_opt.append("m:")

        # Ranges
        x_min = 0
        x_max = n_vmcon + 2
        y_min = 1e-8
        y_max = 1.0

        # Plot
        for ii in range(0, n_ploted_const):
            constraint_name = p_dicts["DICT_ICC_FULL"][
                str(constraints_indexes[dominant_constaints_indexes[ii]])
            ]["name"]
            if n_ploted_const > 5:
                constraint_name = constraint_name.replace(
                    "Central solenoid", "CS"
                ).replace("central solenoid", "CS")
                constraint_name = constraint_name[:17]

            plt.plot(
                plot_indexes,
                constraints[dominant_constaints_indexes[ii]],
                plot_opt[ii],
                label=str(ii + 1)
                + ": "
                + constraint_name
                + " ("
                + str(constraints_indexes[dominant_constaints_indexes[ii]])
                + ")",
            )
        plt.plot(
            plot_indexes, dominant_quad_sum, "k--", label="dominant const quad sum "
        )
        plt.plot(plot_indexes, constraints_quad_sum, "k:", label="total const quad sum")
        plt.plot([x_min, x_max], [tolerance, tolerance], "k-")

        # Actions if the legend gets too big
        if n_ploted_const > 5:
            x_max *= 2
            y_min = 1e-12

        # Cosmetics
        plt.legend(loc="best")
        plt.yscale("log")
        plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
        plt.ylabel("Constraint", fontsize=axis_font_size)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(
            "OPT_plots/constraints_evolution." + save_format, format=save_format
        )
        # plt.show()
        plt.close()
    # ------------------------------------------

    # Plot 4 : Specific constraint
    # ----------------------------
    # Check if the constraint is used
    if plot_const_spe and (constraint_selection not in constraints_indexes):
        print(
            "The constraint {} was not used for this PROCESS run".format(
                constraint_selection
            )
        )
        print(" ->  Try with one of these {}".format(constraints_indexes))

    elif plot_const_spe:
        # Constraint selection
        ii_constraint = constraints_indexes.index(constraint_selection)
        constraint_name = p_dicts["DICT_ICC_FULL"][str(constraint_selection)]["name"]

        # Ranges
        x_min = 0
        x_max = n_vmcon + 2
        y_min = 1e-12
        y_max = 1.0

        # Plot
        plt.plot(plot_indexes, constraints[ii_constraint], "r-", label=constraint_name)
        plt.plot(plot_indexes, constraints_quad_sum, "k:", label="total const quad sum")
        plt.plot([x_min, x_max], [tolerance, tolerance], "k-")

        # Cosmetics
        plt.legend(loc="best")
        plt.yscale("log")
        plt.ylabel("Constraints", fontsize=axis_font_size)
        plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(
            "OPT_plots/"
            + str(p_dicts["DICT_ICC_FULL"][str(constraint_selection)]["name"])
            + "_dominant_constraints_evolution."
            + save_format,
            format=save_format,
        )
        # plt.show()
        plt.close()
    # ----------------------------

    # Plot 5 : Major variable evolution
    # ---------------------------------
    if plot_var_dom and n_var_plots > n_variables:
        print(
            "Impossible to print {} variables as only {} were used for the run".format(
                n_var_plots, n_variables
            )
        )
        print(" -> Plotting {} variables instead".format(n_variables))
        n_var_plots = n_variables

    if plot_var_dom:
        if n_var_plots > len(plot_opt):
            for ii in range(0, n_var_plots - len(plot_opt)):
                plot_opt.append("m:")

        # Select the two variables that vary the most and plot them
        ranked_variables_index = []
        variables_ranges = [max(variable) - min(variable) for variable in variables]

        for ii in range(0, n_variables):
            ranked_variables_index.append(variables_ranges.index(max(variables_ranges)))
            variables_ranges[ranked_variables_index[ii]] = -1.0

        # Plot
        for ii in range(0, n_var_plots):
            leg_label = "{} : {} ({})".format(
                ii + 1,
                p_dicts["DICT_IXC_SIMPLE"][
                    str(variables_indexes[ranked_variables_index[ii]])
                ],
                variables_indexes[ranked_variables_index[ii]],
            )
            plt.plot(
                plot_indexes,
                variables[ranked_variables_index[ii]],
                plot_opt[ii],
                label=leg_label,
            )

        # Action if the legend gets too big ...
        if n_var_plots > 6:
            plt.xlim([0.0, int(1.8 * float(n_vmcon))])

        # Cosmetics
        plt.legend(loc="best")
        plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
        plt.ylabel("variable", fontsize=axis_font_size)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(
            "OPT_plots/dominant_variables_evolution." + save_format, format=save_format
        )
        # plt.show()
        plt.close()
    # ---------------------------------

    # Plot 6 : Specific variable pair path
    # ------------------------------------
    # Check if the variables are used in the considered PROCESS run
    if plot_var_spe_pair and (x_variable_selection not in variables_indexes):
        print(
            "The PROCESS variable {} used of the X axis is not used for the run".format(
                x_variable_selection
            )
        )
        print(" ->  Try with one of these {}".format(variables_indexes))
    elif plot_var_spe_pair and (y_variable_selection not in variables_indexes):
        print(
            "The PROCESS variable {} used of the X axis is not used for the run".format(
                y_variable_selection
            )
        )
        print(" ->  Try with one of these {}".format(variables_indexes))

    elif plot_var_spe_pair:
        # Retrieving the constraint index
        ii_x_variable = variables_indexes.index(x_variable_selection)
        ii_y_variable = variables_indexes.index(y_variable_selection)

        # Plot
        ln_convergence_parameter = [
            np.log10(conv_param) for conv_param in convergence_parameter
        ]  # the ln of the convergence parameter is taken for visibility isues
        scat = plt.scatter(
            variables[ii_x_variable],
            variables[ii_y_variable],
            c=ln_convergence_parameter,
        )
        plt.plot(variables[ii_x_variable], variables[ii_y_variable], "k-")
        plt.grid(True)
        plt.xlabel(
            p_dicts["DICT_IXC_SIMPLE"][str(x_variable_selection)],
            fontsize=axis_font_size,
        )
        plt.ylabel(
            p_dicts["DICT_IXC_SIMPLE"][str(y_variable_selection)],
            fontsize=axis_font_size,
        )
        plot_colorbar = plt.colorbar(scat)
        plot_colorbar.set_label("$ln_{10}$(final conv)", size=axis_font_size)
        plt.tight_layout()
        plt.savefig(
            "OPT_plots/var"
            + str(x_variable_selection)
            + "_vs_var"
            + str(y_variable_selection)
            + "."
            + save_format,
            format=save_format,
        )
        plt.close()
    # ------------------------------------

    # Plot 7 : plot all each inputs variables evolution
    # --------------------------------------------------
    # Creating the folder containing the variables plots
    if plot_var_all:
        if not os.path.isdir("OPT_plots/All_Var"):
            os.mkdir("OPT_plots/All_Var")

        ii = 0
        for variable in variables:

            variable_name = p_dicts["DICT_IXC_SIMPLE"][str(variables_indexes[ii])]
            variable_name = variable_name.replace("(", "").replace(")", "")

            plt.plot(plot_indexes, variable, "k-")
            plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
            plt.ylabel(variable_name, fontsize=axis_font_size)
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(
                "OPT_plots/All_Var/{}_evolution.{}".format(variable_name, save_format),
                format=save_format,
            )
            plt.close()
            ii += 1
    # --------------------------------------------------

    # Plot 8 : plot all each inputs constraint evolution
    # ---------------------------------------------------
    # Creating the folder containing the variables plots
    if plot_const_all:
        if not os.path.isdir("OPT_plots/All_Const"):
            os.mkdir("OPT_plots/All_Const")

        # Ranges
        x_min = 0
        x_max = n_vmcon + 2
        y_min = 1e-12
        y_max = 1.0

        ii = 0
        for constraint in constraints:

            constraint_name = p_dicts["DICT_ICC_FULL"][str(constraints_indexes[ii])][
                "name"
            ]
            plt.plot(plot_indexes, constraint, "r-", label=constraint_name)
            plt.plot(plot_indexes, constraints_quad_sum, "k--", label="const quad sum")

            # Cosmetics
            plt.xlabel("$VMCON$ iteration", fontsize=axis_font_size)
            plt.ylabel(constraint_name, fontsize=axis_font_size)
            plt.legend(loc="best")
            plt.yscale("log")
            plt.axis([x_min, x_max, y_min, y_max])
            plt.grid(True)
            plt.tight_layout()

            constraint_name_f = (
                constraint_name.replace(" ", "_").replace(")", "").replace("(", "")
            )
            plt.savefig(
                "OPT_plots/All_Const/{}_evolution.{}".format(
                    constraint_name_f, save_format
                ),
                format=save_format,
            )
            plt.close()
            ii += 1
    # ---------------------------------------------------
