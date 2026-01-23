"""Plot solution vectors from MFILEs to compare them.

This tool plots multiple solutions, or their differences, to allow comparisons
of the solution vectors and objective function values. It allows normalisation
to a given solution as well as plotting RMSEs relative to it.

If the MFILE is the result of a parameter scan, only the last point is plotted
currently.
"""

import logging
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from warnings import warn

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from process.data_structure import numerics
from process.io.mfile import MFile

# Variables of interest in mfiles and subsequent dataframes
# Be specific about exact names, patterns and regex
INITIAL_NORM_OPT_PARAM_VALUE_PATTERN = "xcm"
RANGE_NORM_OPT_PARAM_VALUE_PATTERN = "nitvar"
OPT_PARAM_VALUE_REGEX = r"^itvar\d{3}$"
NORM_OPT_PARAM_NAME_REGEX = r"itvar\d{3}_name"
NORM_OBJF_PATTERN = "objf"
NORM_OBJF_VALUE = "norm_objf"
NORM_OBJF_NAME = "objf_name"
TAG_REGEX = r"tag$"
TAG = "tag"
CON_VALUE_PATTERN = r"eq_con\d{3}"

# TODO Eventually this (incomplete) dict of short descriptions of constraints
# could be an attribute of the constraints themselves once in Python
# Short constraint descriptions for use on plots
# EQ: equality, LL: lower limit, UL: upper limit
CONSTRAINT_NUMS = {
    "eq_con001": "EQ: beta, temperature, density",
    "eq_con002": "EQ: global plasma power balance",
    "eq_con011": "EQ: radial build",
    "ineq_con030": "UL: injection power",
    "ineq_con015": "LL: L-H threshold power",
    "ineq_con016": "LL: net electric power",
    "ineq_con024": "UL: beta",
    "ineq_con025": "UL: peak toroidal field",
    "ineq_con026": "UL: CS current density at EOF",
    "ineq_con027": "UL: CS current density at BOP",
    "ineq_con033": "UL: SCTF coil operating J",
    "ineq_con034": "UL: CSTF coil dump voltage",
    "ineq_con035": "UL: SCTF coil winding pack J",
    "ineq_con036": "LL: SCTF coil SC temp margin",
    "ineq_con060": "LL: CS SC temp margin",
    "ineq_con062": "LL: alpha:energy confinement times",
    "ineq_con065": "UL: stress on VV during TF quench",
    "ineq_con072": "UL: CS Tresca yield stress",
    "ineq_con081": "LL: central density > pedestal density",
    "ineq_con068": "UL: divertor protection limit",
    "ineq_con031": "UL: SCTF coil case stress",
    "ineq_con032": "UL: SCTF coil conduit stress",
    "ineq_con005": "UL: density",
    "ineq_con008": "UL: neutron wall load",
    "ineq_con009": "UL: fusion power",
    "ineq_con013": "LL: burn time",
}

logger = logging.getLogger(__name__)


@dataclass
class RunMetadata:
    """Metadata for a given Process run.

    Define mfile for run and other information to identify it. Depending on
    what is being run/plotted, different fields can be undefined.
    """

    mfile_path: Path
    tag: str


def plot_mfile_solutions(
    runs_metadata: Sequence[RunMetadata],
    plot_title: str,
    normalising_tag: str | None = None,
    rmse: bool = False,
    normalisation_type: str | None = "init",
) -> tuple[mpl.figure.Figure, pd.DataFrame]:
    """Plot multiple solutions, optionally normalised by a given solution.

    :param runs_metadata: list of RunMetadata objects
    :type runs_metadata: Sequence[RunMetadata]
    :param plot_title: title of plot
    :type plot_title: str
    :param normalising_tag: tag for solution to normalise with. If provided,
    normalise, otherwise don't, defaults to None
    :type normalising_tag: str, optional
    :param rmse: plot RMS errors relative to reference solution, defaults to False
    :type rmse: bool, optional
    :param normalisation_type: opt param normalisation to use: one of ["init", "range", None], defaults to "init"
    :type normalisation_type: str, optional
    :return: figure and dataframe of solutions
    :rtype: Tuple[mpl.figure.Figure, pd.DataFrame]
    """
    if normalisation_type is not None and normalising_tag is not None:
        warn(
            "Double-normalising: using opt params normalised to each solution "
            "and normalising again to another solution. Are you sure?",
            stacklevel=1,
        )

    # Determine type of normalised opt params to plot (i.e. itvar, xcm or nitvar)
    if normalisation_type == "init":
        # Use values normalised to initial value for each solution: final/initial (xcm)
        opt_param_value_pattern = INITIAL_NORM_OPT_PARAM_VALUE_PATTERN
    elif normalisation_type == "range":
        # Use values normalised to allowed range of value (nitvar)
        opt_param_value_pattern = RANGE_NORM_OPT_PARAM_VALUE_PATTERN
    else:
        # No normalisation for each solution (itvar)
        opt_param_value_pattern = OPT_PARAM_VALUE_REGEX

    # Create dataframe from runs metadata: mfile data with a tag for each run
    results_df = _create_df_from_run_metadata(runs_metadata)

    # Filter for tag, optimisation parameters and objective function
    filtered_results_df = _filter_vars_of_interest(
        results_df,
        opt_param_value_pattern=opt_param_value_pattern,
        extra_var_names=["minmax"],
    )

    if normalising_tag is not None:
        # Calculate the normalised diffs relative to the tagged solution
        plot_results_df = _normalise_diffs(
            filtered_results_df,
            opt_param_value_pattern=opt_param_value_pattern,
            normalising_tag=normalising_tag,
        )
    else:
        # Don't perform any processing of optimisation parameters
        plot_results_df = filtered_results_df

    if rmse:
        if normalising_tag is None:
            raise ValueError(
                "RMSE plot requires normalising_tag to be specified: which "
                "solution are the errors relative to?"
            )

        # Calcualte RMS errors relative to normalising tag solution
        rmse_df = _rms_errors(
            results_df=results_df,
            opt_param_value_pattern=opt_param_value_pattern,
            normalising_tag=normalising_tag,
        )
    else:
        # Don't plot RMS errors
        rmse_df = None

    fig = _plot_solutions(
        plot_results_df,
        opt_param_value_pattern=opt_param_value_pattern,
        plot_title=plot_title,
        rmse_df=rmse_df,
        normalising_tag=normalising_tag,
        normalisation_type=normalisation_type,
    )

    # Return fig for optional further customisation
    # fig as varying numbers of axes depending on plot type
    return fig, filtered_results_df


def _extract_mfile_data(mfile_path: Path) -> dict:
    """Extract data from mfile and return as dict.

    Also include the names of the optimisation parameters.

    :param mfile_path: mfile to extract data from
    :type mfile_path: pathlib.Path
    :return: dict of all data in mfile
    :rtype: dict
    """
    mfile = MFile(str(mfile_path))
    mfile_data = {}

    for var in mfile.data:
        if var.startswith("itvar"):
            # Iteration variable: get name too
            mfile_data[f"{var}_name"] = mfile.data[var].var_description

        # Extract value
        mfile_data[var] = mfile.data[var].get_scan(-1)

    return mfile_data


def _create_df_from_run_metadata(runs_metadata: Sequence[RunMetadata]) -> pd.DataFrame:
    """Create a dataframe from multiple mfiles.

    Uses RunMetadata objects.
    :param runs_metadata: scenarios and solvers that have been run
    :type runs_metadata: Sequence[RunMetadata]
    :return: dataframe of all results
    :rtype: pandas.DataFrame
    """
    results = []
    for run_metadata in runs_metadata:
        if Path(run_metadata.mfile_path).exists():
            mfile_data = _extract_mfile_data(run_metadata.mfile_path)
        else:
            raise FileNotFoundError(f"The MFILE {run_metadata.mfile_path} doesn't exist")
        # Merge each run's metadata and results into one dict
        results.append({**asdict(run_metadata), **mfile_data})

    return pd.DataFrame(results)


def _separate_norm_solution(
    results_df: pd.DataFrame, normalising_tag: str
) -> tuple[pd.DataFrame]:
    """Separate solutions df into normalising row and rows to be normalised.

    :param results_df: multiple solutions dataframe
    :type results_df: pd.DataFrame
    :param normalising_tag: tag to identify row to normalise with
    :type normalising_tag: str
    :return: normalising row, rows to be normalised
    :rtype: Tuple[pd.DataFrame]
    """
    # Split results into normalising and non-normalising solutions
    normalising_soln = results_df[results_df[TAG] == normalising_tag]
    non_normalising_solns = results_df[results_df[TAG] != normalising_tag]

    return normalising_soln, non_normalising_solns


def _filter_opt_params(
    results: pd.DataFrame,
    opt_param_value_pattern: str,
    filter_in: bool = True,
) -> pd.DataFrame:
    """Filter optimsation parameters in or out of results.

    :param results: multiple solutions
    :type results: pd.DataFrame
    :param opt_param_value_pattern: normalisation type for opt params in mfile
    :type opt_param_value_pattern: str
    :param filter_in: include opt params, defaults to True
    :type filter_in: bool, optional
    :return: multiple solutions with opt params filtered in or out
    :rtype: pd.DataFrame
    """
    is_opt_param = results.columns.str.contains(opt_param_value_pattern)
    if filter_in:
        # TODO Check this loc can't be removed
        filtered_results = results.loc[:, is_opt_param]
    else:
        filtered_results = results.loc[:, ~is_opt_param]
    return filtered_results


def _normalise_diffs(
    results_df: pd.DataFrame,
    opt_param_value_pattern: str,
    normalising_tag: str,
) -> pd.DataFrame:
    """Normalise differences of multiple solutions with a normalising solution.

    :param results_df: dataframe of two solutions (same scenario, different solvers)
    :type results_df: pandas.DataFrame
    :param opt_param_value_pattern: normalisation type for opt params in mfile
    :type opt_param_value_pattern: str
    :param normalising_tag: tag to normalise other solutions with
    :type normalising_tag: str
    :return: normalised differences
    :rtype: pandas.DataFrame
    """
    normalising_soln, non_normalising_solns = _separate_norm_solution(
        results_df, normalising_tag
    )
    # Solution to normalise with
    normalising_soln_opt_params = _filter_opt_params(
        normalising_soln, opt_param_value_pattern=opt_param_value_pattern
    )
    # Solutions that need normalised diffs calculating
    non_normalising_solns_opt_params = _filter_opt_params(
        non_normalising_solns, opt_param_value_pattern=opt_param_value_pattern
    )

    # Calculate the normalised differences between multiple solutions
    # for optimisation parameters
    normalising_soln_opt_params_np = normalising_soln_opt_params.to_numpy()
    if (normalising_soln_opt_params_np == 0).any():
        zero_indexes = np.nonzero(normalising_soln_opt_params_np == 0)[0]
        raise ValueError(
            f"Can't normalise with 0-valued optimisation parameter at index {zero_indexes!s}."
        )

    normalised_solns_opt_params = (
        non_normalising_solns_opt_params - normalising_soln_opt_params_np
    ) / abs(normalising_soln_opt_params_np)

    # Combine dfs to get non-numerical info (metadata) alongside normalised diffs
    # e.g. tag, objf_name
    non_normalising_solns_metadata = _filter_opt_params(
        non_normalising_solns,
        opt_param_value_pattern=opt_param_value_pattern,
        filter_in=False,
    )
    return pd.concat(
        [non_normalising_solns_metadata, normalised_solns_opt_params], axis=1
    )


def _filter_vars_of_interest(
    results_df: pd.DataFrame,
    opt_param_value_pattern: str | None = None,
    constraints: bool = False,
    extra_var_names: list[str] | None = None,
) -> pd.DataFrame:
    """Filter variables of interest from full results for all solutions.

    :param results_df: full results for all solutions
    :type results_df: pandas.DataFrame
    :param opt_param_value_pattern: normalisation type for opt params in mfile,
    set if plotting opt params, defaults to None
    :type opt_param_value_pattern: str | None, optional
    :param constraints: filter in constraint values, defaults to False
    :type constraints: bool, optional
    :param extra_var_names: other variables of interest to filter, defaults to None
    :type extra_var_names: List[str], optional
    :return: variables of interest
    :rtype: pandas.DataFrame
    """
    if extra_var_names is None:
        extra_var_names = []

    if constraints:
        # Filter in tag, equality and inequality constraint values only
        return results_df.filter(regex=(TAG_REGEX + r"|" + CON_VALUE_PATTERN))

    # Filter for optimisation parameters (normalised to initial value
    # e.g. xcm001) values and names, objective function value and name, plus
    # any extra variables
    return results_df.filter(
        regex=(
            f"{TAG_REGEX}|{opt_param_value_pattern}|{NORM_OPT_PARAM_NAME_REGEX}|"
            f"{NORM_OBJF_PATTERN}|{'|'.join(extra_var_names)}".rstrip("|")
        )
    )


def _plot_solutions(
    diffs_df: pd.DataFrame,
    normalisation_type: str | None,
    opt_param_value_pattern: str,
    plot_title: str,
    normalising_tag: str | None,
    rmse_df: pd.DataFrame,
) -> mpl.figure.Figure:
    """Plot multiple solutions, optionally normalised by a given solution.

    :param diffs_df: normalised diffs for optimisation parameters and objective function
    :type diffs_df: pandas.DataFrame
    :param normalisation_type: opt param normalisation to use: ["init", "range", None]
    :type normalisation_type: str | None
    :param opt_param_value_pattern: normalisation type for opt params in mfile
    :type opt_param_value_pattern: str
    :param plot_title: title of plot
    :type plot_title: str
    :param normalising_tag: tag for normalising solution, if any
    :type normalising_tag: Union[str, None]
    :param rmse_df: RMS errors relative to reference solution
    :type rmse_df: pd.DataFrame
    :return: figure containing varying numbers of axes
    :rtype: mpl.figure.Figure
    """
    # Separate optimisation parameters and objective dfs
    opt_params_df = diffs_df.filter(
        regex=f"{opt_param_value_pattern}|{NORM_OPT_PARAM_NAME_REGEX}|{TAG}"
    )
    norm_objf_df = diffs_df.filter(regex=f"{NORM_OBJF_PATTERN}|{TAG}")

    # Further separate objective function values only, for plotting
    norm_objf_values_df = norm_objf_df.filter(regex=f"{NORM_OBJF_VALUE}|{TAG}")

    # Acquire objective function name(s), then check only one type is being plotted
    # The objective function is found by first checking the NORM_OBJF_NAME column.
    # If this entry does not exist for any of the MFiles (ie they are old),
    # then the code 'falls back' to the the minmax output. Which holds the same
    # information, but is less descriptive.
    if (
        NORM_OBJF_NAME in norm_objf_df.columns
        and not norm_objf_df[NORM_OBJF_NAME].isna().any()
    ):
        objf_list = norm_objf_df[NORM_OBJF_NAME].unique()
    else:
        numerics.init_numerics()
        objf_list = {
            numerics.lablmm[int(abs(minmax)) - 1] for minmax in diffs_df["minmax"]
        }

    if len(objf_list) != 1:
        raise ValueError("Can't plot different objective functions on the same plot")

    objf_list = list(objf_list)  # Convert set to list
    objf_name = str(objf_list[0])

    # Now separate optimisation parameter values from their names
    opt_params_values_df = opt_params_df.filter(regex=f"{opt_param_value_pattern}|{TAG}")
    opt_params_names_df = opt_params_df.filter(regex=NORM_OPT_PARAM_NAME_REGEX)

    # Replace xcm--- optimisation parameter column headers with actual var names
    # in sub-df: allows plotting showing actual var names
    # eg. column headers "tag, xcm001, xcm002, ..." become
    # "tag, b_plasma_toroidal_on_axis, rmajor, ..."
    # Normalising row may have been filtered out; reset index to ensure
    # opt param names in row 0
    opt_params_names_df_reset = opt_params_names_df.reset_index(drop=True)
    opt_params_names = opt_params_names_df_reset.loc[0].to_numpy().tolist()

    # Need to include tag column as first column header
    opt_params_names_with_tag = opt_params_names.copy()
    opt_params_names_with_tag.insert(0, TAG)
    opt_params_values_with_names_df = opt_params_values_df.set_axis(
        opt_params_names_with_tag, axis=1
    )

    # Define subfig and subplot layout
    # Matplotlib usually fits the subplots to the figure. However, as there are
    # varying numbers of opt params and various plot options, it is necessary
    # to fit the figure to the subplots, i.e. have different figure sizes
    # based on what's being plotted
    # Strategy: 2 subfigures (opt params and obj func/rmse), varying numbers of
    # subplots in each

    # Define optimisation parameter plot
    opt_param_count = len(opt_params_names)
    if normalisation_type is None and normalising_tag is None:
        # Actual values of opt params, not normalised to any solution
        # Separate axes required for each opt param: more space per param
        inches_per_opt_param = 0.7
        # Add extra subplot for legend
        nrows_opt_param_subplot = opt_param_count + 1
        inches_per_legend_entry = 0.5
        tag_count = diffs_df[TAG].count()
        external_legend_height = tag_count * inches_per_legend_entry
    else:
        # Single axis for all opt params: less space per param
        inches_per_opt_param = 0.20
        nrows_opt_param_subplot = 1
        external_legend_height = 0.0

    final_subfig_height_opt_params = (
        opt_param_count * inches_per_opt_param
    ) + external_legend_height
    figsize = [7, final_subfig_height_opt_params]

    # Add obj func/RMSE subfig
    base_subfig_height_obj_func = 0.75
    if rmse_df is not None:
        nrows_obj_func_rmse_subplot = 2
        final_subfig_height_obj_func = (
            nrows_obj_func_rmse_subplot * base_subfig_height_obj_func
        )
    else:
        nrows_obj_func_rmse_subplot = 1
        final_subfig_height_obj_func = base_subfig_height_obj_func

    figsize[1] += final_subfig_height_obj_func

    fig = plt.figure(layout="constrained", figsize=figsize)
    fig.suptitle(plot_title)
    # 2 subfig rows: 1 for opt params, 1 for objective function and optional RMSE
    # Use actual calculated height to define height ratios
    subfigs = fig.subfigures(
        nrows=2,
        ncols=1,
        height_ratios=[final_subfig_height_opt_params, final_subfig_height_obj_func],
    )
    axs_opt_params = subfigs[0].subplots(nrows=nrows_opt_param_subplot)

    # Adapt x axis label for normalisation type
    if normalisation_type == "init":
        x_axis_label = "Value normalised by initial value of each solution"
    elif normalisation_type == "range":
        x_axis_label = "Value normalised by parameter range"
    else:
        x_axis_label = "Value"

    if normalising_tag is not None:
        x_axis_label += f", normalised to {normalising_tag}"

    # Plot optimisation parameters
    # Melt df (wide to long-form) for seaborn plotting with jitter
    opt_params_values_with_names_df_melt = opt_params_values_with_names_df.melt(
        id_vars=TAG
    )

    # If normalisation_type is None (no normalisation of original opt param
    # values), plot strip plot with different axis for each opt param.
    # Otherwise plot all params on single axis
    if normalisation_type is None and normalising_tag is None:
        # axs_opt_params[0] reserved for legend
        for i, opt_param_name in enumerate(opt_params_names):
            sns.stripplot(
                data=opt_params_values_with_names_df_melt[
                    opt_params_values_with_names_df_melt["variable"] == opt_param_name
                ],
                x="value",
                y="variable",
                hue=TAG,
                jitter=True,
                ax=axs_opt_params[i + 1],
            )
            axs_opt_params[i + 1].set_xlabel("")
            axs_opt_params[i + 1].set_ylabel("")
            axs_opt_params[i + 1].get_legend().remove()
            axs_opt_params[i + 1].ticklabel_format(
                axis="x", style="scientific", scilimits=(-1, 1)
            )
        # Set legend for all subplots: plot in own subplot
        axs_opt_params[0].axis("off")
        handles, labels = axs_opt_params[1].get_legend_handles_labels()
        axs_opt_params[0].legend(handles=handles, labels=labels, ncols=3, loc="center")

        subfigs[0].supxlabel(x_axis_label)
        subfigs[0].supylabel("Optimisation parameter")
    else:
        sns.stripplot(
            data=opt_params_values_with_names_df_melt,
            x="value",
            y="variable",
            hue=TAG,
            jitter=True,
            ax=axs_opt_params,
        )

        axs_opt_params.set_xlabel(x_axis_label)
        axs_opt_params.set_ylabel("Optimisation parameter")
        axs_opt_params.legend()
        axs_opt_params.grid()

    # Plot objf change separately
    axs_opt_params = subfigs[1].subplots(nrows=nrows_obj_func_rmse_subplot)
    if nrows_obj_func_rmse_subplot > 1:
        ax_obj_func = axs_opt_params[0]
        ax_rmse = axs_opt_params[1]
    else:
        ax_obj_func = axs_opt_params

    # Melt for seaborn stripplot
    norm_objf_values_df_melt = norm_objf_values_df.melt(id_vars=TAG)
    sns.stripplot(
        data=norm_objf_values_df_melt,
        x="value",
        y="variable",
        hue=TAG,
        jitter=True,
        ax=ax_obj_func,
        formatter=lambda _label: objf_name,
    )

    ax_obj_func.get_legend().remove()
    # Objective function values are not normalised (no initial value): be explicit
    ax_obj_func.set_xlabel("Objective function value")
    ax_obj_func.set_ylabel("")

    if rmse_df is not None:
        # Ensure solution legend colours are the same: tags match between opt
        # params and RMS error plot. This can happen if not normalising opt
        # params, but RMS error plot requested
        rmse_df_matched = rmse_df[
            rmse_df["tag"].isin(opt_params_values_with_names_df["tag"])
        ]

        rmse_df_matched_melt = rmse_df_matched.melt(id_vars=TAG)
        sns.stripplot(
            data=rmse_df_matched_melt,
            x="value",
            y="variable",
            hue=TAG,
            jitter=True,
            ax=ax_rmse,
            formatter=lambda _label: f"Relative to\n{normalising_tag}",
        )
        ax_rmse.set_xlabel("RMS error")
        ax_rmse.set_ylabel("")
        ax_rmse.get_legend().remove()

    return fig


def _rms_errors(
    results_df: pd.DataFrame,
    opt_param_value_pattern: str,
    normalising_tag: str,
) -> pd.DataFrame:
    """Calculate RMS errors between different solutions.

    Summarise solution differences in single number.
    :param results_df: df containing multiple solutions
    :type results_df: pandas.DataFrame
    :param opt_param_value_pattern: normalisation type for opt params in mfile
    :type opt_param_value_pattern: str
    :param normalising_tag: tag to calculate RMSEs against
    :type normalising_tag: str
    :return: RMS errors for each solution
    :rtype: pandas.DataFrame
    """
    normalising_soln, _ = _separate_norm_solution(results_df, normalising_tag)
    # Reference solution to calculate RMS errors against
    ref_soln_opt_params = _filter_opt_params(
        normalising_soln, opt_param_value_pattern=opt_param_value_pattern
    )
    # Solutions that need RMS errors calculating
    # Include reference solution here, so that its RMS error is calculated (as 0)
    # This is so that the colours in the RMSE plot match those in the main
    # optimisation parameter plot's legend
    all_solns_opt_params = _filter_opt_params(
        results_df, opt_param_value_pattern=opt_param_value_pattern
    )

    ref_soln_opt_params_np = ref_soln_opt_params.to_numpy()
    solns_opt_params_to_rms_np = all_solns_opt_params.to_numpy()

    rmses = np.sqrt(
        np.mean((solns_opt_params_to_rms_np - ref_soln_opt_params_np) ** 2, axis=1)
    )

    # Create df of tags and rmses
    tag_series = results_df["tag"]
    rmse_series = pd.Series(rmses, name="rmse")
    return pd.concat([tag_series, rmse_series], axis=1)


def _plot_solutions_constraints(df: pd.DataFrame, title: str) -> mpl.figure.Figure:
    """Create plot from constraint values df.

    :param df: constraint values for each mfile
    :type df: pd.DataFrame
    :param title: plot title
    :type title: str
    :return: created figure
    :rtype: mpl.figure.Figure
    """
    # Replace constraint numbers with descriptions
    df = df.rename(mapper=CONSTRAINT_NUMS, axis="columns")
    df_melt = df.melt(id_vars=TAG)
    fig, ax = plt.subplots()
    sns.stripplot(
        data=df_melt,
        x="value",
        y="variable",
        hue=TAG,
        jitter=True,
        ax=ax,
    )

    ax.set_xlabel("Constraint value")
    ax.set_ylabel("Constraint")
    ax.legend()
    ax.grid()
    ax.set_title(title)
    return fig


def plot_mfile_solutions_constraints(
    runs_metadata: Sequence[RunMetadata], title: str
) -> tuple[mpl.figure.Figure, pd.DataFrame]:
    """Plot constraint values in mfiles.

    :param runs_metadata: list of RunMetadata objects
    :type runs_metadata: Sequence[RunMetadata]
    :param title: plot title
    :type title: str
    :return: figure and dataframe of solutions
    :rtype: Tuple[mpl.figure.Figure, pd.DataFrame]
    """
    # Create dataframe from runs metadata: mfile data with a tag for each run
    results_df = _create_df_from_run_metadata(runs_metadata)

    # Filter for tag and constraint values
    filtered_results_df = _filter_vars_of_interest(results_df, constraints=True)

    fig = _plot_solutions_constraints(df=filtered_results_df, title=title)
    return fig, filtered_results_df
