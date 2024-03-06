"""Plot solution vectors from MFILEs to compare them.

This tool plots multiple solutions, or their differences, to allow comparisons
of the solution vectors and objective function values.

This tool can:
1. Plot solutions from any number of arbitrary MFILEs, or their differences
(plot_mfile_solutions())
"""

from process.io.mfile import MFile
from pathlib import Path
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import logging
import seaborn as sns
from dataclasses import dataclass, asdict
from typing import Optional, Sequence, List, Dict


# Variables of interest in mfiles and subsequent dataframes
# Be specific about exact names, patterns and regex
NORM_OPT_PARAM_VALUE_PATTERN = "xcm"
NORM_OPT_PARAM_NAME_REGEX = r"itvar\d{3}_name"
NORM_OBJF_PATTERN = "objf"
NORM_OBJF_VALUE = "norm_objf"
NORM_OBJF_NAME = "objf_name"
TAG_REGEX = r"tag$"
TAG = "tag"

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
    normalise: Optional[bool] = False,
    normalising_tag: Optional[str] = None,
    rmse: Optional[bool] = False,
) -> pd.DataFrame:
    """Plot multiple solutions, optionally normalised by a given solution.

    :param runs_metadata: list of RunMetadata objects
    :type runs_metadata: Sequence[RunMetadata]
    :param plot_title: title of plot
    :type plot_title: str
    :param normalise: normalise optimisation parameters to a given solution, defaults to False
    :type normalise: Optional[bool], optional
    :param normalising_tag: tag for solution to normalise with, defaults to None
    :type normalising_tag: str, optional
    :param rmse: plot RMS errors relative to reference solution, defaults to False
    :type rmse: bool, optional
    """
    # TODO Should return ax too?
    # Create dataframe from runs metadata: mfile data with a tag for each run
    results_df = _create_df_from_run_metadata(runs_metadata)

    # Filter for tag, optimisation parameters and objective function
    filtered_results_df = _filter_vars_of_interest(results_df)

    if normalise:
        # Calculate the normalised diffs relative to the tagged solution
        plot_results_df = _normalise_diffs(
            filtered_results_df, normalising_tag=normalising_tag
        )
    else:
        # Don't perform any processing of optimisation parameters
        plot_results_df = filtered_results_df

    if rmse:
        # Calcualte RMS errors relative to normalising tag solution
        rmse_df = _rms_errors(results_df=results_df, normalising_tag=normalising_tag)
    else:
        # Don't plot RMS errors
        rmse_df = None

    _plot_solutions(
        plot_results_df,
        plot_title=plot_title,
        normalise=normalise,
        rmse_df=rmse_df,
        normalising_tag=normalising_tag,
    )

    return filtered_results_df


def _extract_mfile_data(mfile_path: Path) -> Dict:
    """Extract data from mfile and return as dict.

    Also include the names of the optimisation parameters.

    :param mfile_path: mfile to extract data from
    :type mfile_path: pathlib.Path
    :return: dict of all data in mfile
    :rtype: dict
    """
    mfile = MFile(str(mfile_path))
    mfile_data = {}

    for var in mfile.data.keys():
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
            raise FileNotFoundError(
                f"The MFILE {run_metadata.mfile_path} doesn't exist"
            )
        # Merge each run's metadata and results into one dict
        results.append({**asdict(run_metadata), **mfile_data})

    return pd.DataFrame(results)


def _separate_norm_solution(
    results_df: pd.DataFrame, normalising_tag: Optional[str] = None
) -> tuple[pd.DataFrame]:
    """Separate solutions df into normalising row and rows to be normalised.

    :param results_df: multiple solutions dataframe
    :type results_df: pd.DataFrame
    :param normalising_tag: tag to identify row to normalise with, defaults to None
    :type normalising_tag: Optional[str], optional
    :return: normalising row, rows to be normalised
    :rtype: tuple[pd.DataFrame]
    """
    # If no normalising tag specified, use first solution to normalise
    if normalising_tag is None:
        normalising_tag = results_df.loc[0][TAG]

    # Split results into normalising and non-normalising solutions
    normalising_soln = results_df[results_df[TAG] == normalising_tag]
    non_normalising_solns = results_df[results_df[TAG] != normalising_tag]

    return normalising_soln, non_normalising_solns


def _filter_opt_params(results: pd.DataFrame, filter_in: bool = True) -> pd.DataFrame:
    """Filter optimsation parameters in or out of results.

    :param results: multiple solutions
    :type results: pd.DataFrame
    :param filter_in: include opt params, defaults to True
    :type filter_in: bool, optional
    :return: multiple solutions with opt params filtered in or out
    :rtype: pd.DataFrame
    """
    is_opt_param = results.columns.str.contains(f"{NORM_OPT_PARAM_VALUE_PATTERN}")
    if filter_in:
        # TODO Check this loc can't be removed
        filtered_results = results.loc[:, is_opt_param]
    else:
        filtered_results = results.loc[:, ~is_opt_param]
    return filtered_results


def _normalise_diffs(
    results_df: pd.DataFrame, normalising_tag: Optional[str] = None
) -> pd.DataFrame:
    """Normalise differences of multiple solutions with a normalising solution.

    :param results_df: dataframe of two solutions (same scenario, different solvers)
    :type results_df: pandas.DataFrame
    :param normalising_tag: tag to normalise other solutions with, defaults to None
    :type normalising_tag: str, optional
    :return: normalised differences
    :rtype: pandas.DataFrame
    """
    normalising_soln, non_normalising_solns = _separate_norm_solution(
        results_df, normalising_tag
    )
    # Solution to normalise with
    normalising_soln_opt_params = _filter_opt_params(normalising_soln)
    # Solutions that need normalised diffs calculating
    non_normalising_solns_opt_params = _filter_opt_params(non_normalising_solns)

    # Calculate the normalised differences between multiple solutions
    # for optimisation parameters
    normalising_soln_opt_params_np = normalising_soln_opt_params.to_numpy()
    normalised_solns_opt_params = (
        non_normalising_solns_opt_params - normalising_soln_opt_params_np
    ) / abs(normalising_soln_opt_params_np)

    # Combine dfs to get non-numerical info (metadata) alongside normalised diffs
    # e.g. tag, objf_name
    non_normalising_solns_metadata = _filter_opt_params(
        non_normalising_solns, filter_in=False
    )
    normalised_solns = pd.concat(
        [non_normalising_solns_metadata, normalised_solns_opt_params], axis=1
    )

    return normalised_solns


def _filter_vars_of_interest(
    results_df: pd.DataFrame, extra_var_names: Optional[List[str]] = None
) -> pd.DataFrame:
    """Filter variables of interest from full results for all solutions.

    :param results_df: full results for all solutions
    :type results_df: pandas.DataFrame
    :param extra_var_names: other variables of interest to filter, defaults to None
    :type extra_var_names: List[str], optional
    :return: variables of interest
    :rtype: pandas.DataFrame
    """
    if extra_var_names is None:
        extra_var_names = []

    # Filter for optimisation parameters (normalised to initial value
    # e.g. xcm001) values and names, objective function value and name, plus
    # any extra variables
    filtered_results = results_df.filter(
        regex=(
            f"{TAG_REGEX}|{NORM_OPT_PARAM_VALUE_PATTERN}|{NORM_OPT_PARAM_NAME_REGEX}|"
            f"{NORM_OBJF_PATTERN}|{'|'.join(extra_var_names)}".rstrip("|")
        )
    )

    return filtered_results


def _plot_solutions(
    diffs_df: pd.DataFrame,
    plot_title: str,
    normalise: bool,
    normalising_tag: str,
    rmse_df: pd.DataFrame,
):
    """Plot multiple solutions, optionally normalised by a given solution.

    :param diffs_df: normalised diffs for optimisation parameters and objective function
    :type diffs_df: pandas.DataFrame
    :param plot_title: title of plot
    :type plot_title: str
    :param normalise: normalise optimisation parameters to a given solution
    :type normalise: bool
    :param normalising_tag: tag for normalising solution
    :type normalising_tag: str
    :param rmse_df: RMS errors relative to reference solution
    :type rmse_df: pd.DataFrame
    """
    # Separate optimisation parameters and objective dfs
    opt_params_df = diffs_df.filter(
        regex=f"{NORM_OPT_PARAM_VALUE_PATTERN}|{NORM_OPT_PARAM_NAME_REGEX}|{TAG}"
    )
    norm_objf_df = diffs_df.filter(regex=f"{NORM_OBJF_PATTERN}|{TAG}")

    # Further separate objective function values only, for plotting
    norm_objf_values_df = norm_objf_df.filter(regex=f"{NORM_OBJF_VALUE}|{TAG}")

    # Acquire objective function name(s), then check only one type is being plotted
    objf_list = norm_objf_df[NORM_OBJF_NAME].unique()

    if len(objf_list) != 1:
        raise ValueError("Can't plot different objective functions on the same plot")

    objf_name = objf_list[0]

    # Now separate optimisation parameter values from their names
    opt_params_values_df = opt_params_df.filter(
        regex=f"{NORM_OPT_PARAM_VALUE_PATTERN}|{TAG}"
    )
    opt_params_names_df = opt_params_df.filter(regex=NORM_OPT_PARAM_NAME_REGEX)

    # Replace xcm--- optimisation parameter column headers with actual var names
    # in sub-df: allows plotting showing actual var names
    # eg. column headers "tag, xcm001, xcm002, ..." become
    # "tag, bt, rmajor, ..."
    # Normalising row may have been filtered out; reset index to ensure
    # opt param names in row 0
    opt_params_names_df_reset = opt_params_names_df.reset_index(drop=True)
    opt_params_names = opt_params_names_df_reset.loc[0].values.tolist()

    # Need to include tag column as first column header
    opt_params_names.insert(0, TAG)
    opt_params_values_with_names_df = opt_params_values_df.set_axis(
        opt_params_names, axis=1
    )

    # Separate optimisation parameters and objective function subplots
    # Include space for RMS errors or not
    if rmse_df is not None:
        fig, ax = plt.subplots(ncols=3, width_ratios=[8, 1, 1])
    else:
        fig, ax = plt.subplots(ncols=2, width_ratios=[8, 1])
    fig.set_size_inches(10, 6)
    fig.suptitle(plot_title)

    # Strip plot for optimisation parameters
    # Melt df (wide to long-form) for seaborn plotting with jitter
    opt_params_values_with_names_df_melt = opt_params_values_with_names_df.melt(
        id_vars=TAG
    )
    sns.stripplot(
        data=opt_params_values_with_names_df_melt,
        x="variable",
        y="value",
        hue=TAG,
        jitter=True,
        ax=ax[0],
        formatter=lambda label: label.lstrip("xcm0"),
    )
    # TODO Still require this formatter?

    if normalise:
        ax[0].set_ylabel(f"Normalised difference to {normalising_tag}")
    else:
        ax[0].set_ylabel("Value normalised to initial point")

    ax[0].set_xlabel("Optimisation parameter")
    ax[0].legend()
    ax[0].tick_params("x", labelrotation=90)
    ax[0].grid(axis="x")

    # Plot objf change separately
    # Melt for seaborn stripplot
    norm_objf_values_df_melt = norm_objf_values_df.melt(id_vars=TAG)
    sns.stripplot(
        data=norm_objf_values_df_melt,
        x="variable",
        y="value",
        hue=TAG,
        jitter=True,
        ax=ax[1],
        formatter=lambda label: objf_name,
    )

    ax[1].get_legend().remove()
    # Objective function values are not normalised (no initial value): be explicit
    ax[1].set_ylabel("Value")

    ax[1].set_xlabel("Objective\nfunction")

    if rmse_df is not None:
        # Ensure solution legend colours are the same: tags match between opt
        # params and RMS error plot. This can happen if not normalising opt
        # params, but request RMS error plot
        rmse_df_matched = rmse_df[
            rmse_df["tag"].isin(opt_params_values_with_names_df["tag"])
        ]

        rmse_df_matched_melt = rmse_df_matched.melt(id_vars=TAG)
        sns.stripplot(
            data=rmse_df_matched_melt,
            x="variable",
            y="value",
            hue=TAG,
            jitter=True,
            ax=ax[2],
        )
        ax[2].set_ylabel("")
        ax[2].set_xlabel(f"RMS error\nrelative to\n{normalising_tag}")
        ax[2].get_legend().remove()

    # Ensure title doesn't overlap plots
    fig.tight_layout()


def _rms_errors(
    results_df: pd.DataFrame, normalising_tag: Optional[str] = None
) -> pd.DataFrame:
    """Calculate RMS errors between different solutions.

    Summarise solution differences in single number.
    :param results_df: df containing multiple solutions
    :type results_df: pandas.DataFrame
    :param normalising_tag: tag to calculate percentage changes against, defaults to None
    :type normalising_tag: str, optional
    :return: RMS errors for each solution
    :rtype: pandas.DataFrame
    """
    normalising_soln, _ = _separate_norm_solution(results_df, normalising_tag)
    # Reference solution to calculate RMS errors against
    ref_soln_opt_params = _filter_opt_params(normalising_soln)
    # Solutions that need RMS errors calculating
    # Include reference solution here, so that its RMS error is calculated (as 0)
    # This is so that the colours in the RMSE plot match those in the main
    # optimisation parameter plot's legend
    all_solns_opt_params = _filter_opt_params(results_df)

    ref_soln_opt_params_np = ref_soln_opt_params.to_numpy()
    solns_opt_params_to_rms_np = all_solns_opt_params.to_numpy()

    rmses = np.sqrt(
        np.mean((solns_opt_params_to_rms_np - ref_soln_opt_params_np) ** 2, axis=1)
    )

    # Create df of tags and rmses
    tag_series = results_df["tag"]
    rmse_series = pd.Series(rmses, name="rmse")
    rmse_df = pd.concat([tag_series, rmse_series], axis=1)
    return rmse_df
