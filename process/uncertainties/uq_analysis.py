import pandas as pd
import argparse
import json
import os
from SALib.analyze import rbd_fast, hdmr, rsa
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import figure
from shapely.geometry import LineString


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Calulate sensitivity analysis indices for Monte Carlo runs of PROCESS."
    )

    parser.add_argument(
        "-iv",
        "--input_sampled_vars",
        default="param_values.h5",
        help="input hdf5 file with sample points (default=param_values.h5)",
    )

    parser.add_argument(
        "-uf",
        "--uq_folder",
        default="uq_data_folder",
        help="Folder containing all UQ files",
    )

    parser.add_argument(
        "-f",
        "--configfile",
        default="config_evaluate_uncertainties.json",
        help="configuration file",
    )

    parser.add_argument(
        "-iu",
        "--input_uncertainty",
        default="uncertainties_data.h5",
        help="input hdf5 file with uncertainty data (default=param_values.h5)",
    )

    parser.add_argument(
        "-fom",
        "--figure_of_merit",
        default="rmajor",
        help="figure of merit, must be pulled from input uncertainty (default=param_values.h5)",
    )

    return parser.parse_args(args)


class UncertaintyData:
    """Collects and analyses the output from the evaluate_uncertainties.py tool. Supply an instance with the input
    folder for the run data. The tool looks for hdf files containing uncertainty data, merges them, and has functions
    to clean, analyse, and plot the data."""

    def __init__(self, path_to_uq_data_folder, figure_of_merit):
        self.path_in = path_to_uq_data_folder
        self.figure_of_merit = figure_of_merit

        self.sample_vars_h5_path = "param_values.h5"
        self.uncertainties_h5_path = "uncertainties_data.h5"
        self.uncertainties_df, self.sampled_vars_df = self.merge_hdf_files()
        self.uncertainties_df = self.uncertainties_df.drop(
            columns=[
                "procver",
                "date",
                "time",
                "username",
                "runtitle",
                "tagno",
                "branch_name",
                "commsg",
                "fileprefix",
                "tauelaw",
                # "radius_of_plasma-facing_side_of_inner_leg_should_be_[m]",
                "error_id",
            ]
        )
        # Remove columns which have the same value in every row.
        self.unique_array = unique_cols(self.uncertainties_df)
        self.uncertainties_df = self.uncertainties_df.loc[:, ~self.unique_array]
        self.input_names = self.sampled_vars_df.columns.tolist()
        self.itv = [
            "bt",
            "te",
            "beta",
            "dene",
            "tfcth",
            "wallmw",
            "ohcth",
            "bigq",
            "bore",
            "betalim",
            "coheof",
            "cohbop",
            "gapoh",
            "fvsbrnni",
            "itvar019",
            "itvar020",
            "jwptf",
            "vtfskv",
            "vdalw",
            "tdmptf",
            "thkcas",
            "thwcndut",
            "fcutfsu",
            "cpttf",
            "gapds",
            "plhthresh",
            "tmargtf",
            "tmargoh",
            "oh_steel_frac",
            "pdivt",
            "powfmw",
        ]
        # self.input_names.extend(itv)
        self.sampled_vars_to_plot = []
        self.plot_names = []
        self.number_sampled_vars = len(self.input_names)
        self.problem = {
            "num_vars": self.number_sampled_vars,
            "names": self.input_names,
        }
        # This section was used to combine the params_df and the uncertainties_df
        # However I found that this is problematic when a run doesn't converge
        # and the lengths of the dfs are different, because the params df is written
        # at the start of the MC loop.
        # Therefore, now we draw the params from the uncertainties_df.
        # You must check that your params are written to mfile or they won't appear.

        # self.converged_df = self.converged_df[
        #     self.converged_df["etath"].between(0.38, 0.39)
        # ]
        # self.unconverged_df = self.unconverged_df[
        #     self.unconverged_df["etath"].between(0.38, 0.39)
        # ]

        # Limit the dataset, useful for investigations.
        # self.converged_df = self.converged_df.iloc[:589]
        # self.unconverged_df = self.unconverged_df.iloc[:589]
        # self.combined_df = pd.concat(
        #     [
        #         self.sampled_vars_df,
        #         self.uncertainties_df.reset_index(drop=True),
        #     ],
        #     axis=1,
        # )
        # # Removed duplicated parameters as a result of combining the two DFs.
        # self.combined_df = self.combined_df.loc[
        #     :, ~self.combined_df.columns.duplicated()
        # ].copy()

        self.converged_df = self.uncertainties_df[self.uncertainties_df["ifail"] == 1.0]
        self.unconverged_df = self.uncertainties_df[
            self.uncertainties_df["ifail"] != 1.0
        ]
        self.converged_sampled_vars_df = self.converged_df[self.input_names]
        self.unconverged_sampled_vars_df = self.unconverged_df[self.input_names]
        self.unconverged_fom_df = self.uncertainties_df[self.figure_of_merit]
        self.plot_converged_df = pd.DataFrame()
        self.plot_unconverged_df = pd.DataFrame()
        self.sensitivity_df = pd.DataFrame()
        self.sumsq_sensitivity_df = pd.DataFrame()

        # Significance level is used for plotting, indicates level below which we consider
        # index values insignificant, 0.05 by default - subjective.
        self.significance_level = 0.05
        self.reliability_index = 0.0
        self.number_of_converged_runs = len(self.converged_df.index)
        self.number_of_unconverged_runs = len(self.unconverged_df)
        self.item_with_sensitivity = []
        self.sens_dict = {}
        self.list_new_dfs = []

        # Using a dict for converting param names for plotting
        self.name_dict = {
            "walalw": "Max neutron wall-load (MW/m$^{2}$)",
            "kappa": "Plasma separatrix elongation",
            "triang": "Plasma separatrix triangularity",
            "ralpne": "Thermal alpha/electron density (%)",
            "etaech": "ECH wall plug to injector eff. (%)",
            "pinjalw": "Max injected power (MW)",
            "alstroh": "Allowable hoop stress in CS (Pa)",
            "coreradius": "Normalised core radius (m)",
            "sig_tf_wp_max": "Max sheer stress in TF coil (Pa)",
            "psepbqarmax": "Divertor protection (MWT/m)",
            "sig_tf_case_max": "Max stress in TF coil case (Pa)",
            "powfmw": "Fusion power (MW)",
            "etath": "Thermal to electric eff. (%)",
        }

    def configure_data_for_plotting(self):
        """This function sorts the UQ data into dataframes for plotting with the hdf_to_scatter tool."""
        if len(self.sampled_vars_to_plot) == 0:
            print("Plotting all sampled parameters")
            self.plot_names.extend(self.input_names)
            self.plot_names.extend([self.figure_of_merit])

            self.plot_converged_df = self.converged_df[self.plot_names]
            self.plot_unconverged_df = self.unconverged_df[self.plot_names]
        else:
            print("Ploting user named parameters")
            self.plot_names.extend(self.sampled_vars_to_plot)
            self.plot_names.extend([self.figure_of_merit])

            self.plot_converged_df = self.converged_df[self.plot_names]
            self.plot_unconverged_df = self.unconverged_df[self.plot_names]

    def get_fom_converged_df(self, figure_of_merit):
        """Get the Figure of Merit (fom) from the dataframe containing converged runs.

        :param fom: Figure of Merit
        :type fom: str
        :return: Converged FOM DataFrame
        :rtype: pandas.Series
        """
        return self.converged_df[figure_of_merit]

    def calculate_sensitivity(self, figure_of_merit):
        """Calculate the sensitivity indices for a set of converged UQ runs.
        Uses the Salib rbd_fast analysis method.
        """
        converged_figure_of_merit_df = self.get_fom_converged_df(figure_of_merit)
        sampled_vars_df = self.converged_sampled_vars_df
        self.problem = {
            "num_vars": self.number_sampled_vars,
            "names": self.input_names,
        }
        sirbd_fast = rbd_fast.analyze(
            self.problem,
            sampled_vars_df.to_numpy(),
            converged_figure_of_merit_df.to_numpy(),
        )
        self.sensitivity_df = sirbd_fast.to_df()
        self.sensitivity_df = self.sensitivity_df.sort_values(by="S1", ascending=False)
        # print(" ")
        # print(figure_of_merit)
        # print(self.sensitivity_df)

    def find_significant_parameters(self, sensitivity_data, significance_level):
        """Find the parameters above a given significance level. This is for RBD Fast
        DataFrame with heading "S1".

        :param sensitivity_data: Dataframe with sensitivity indices
        :type sensitivity_data: pandas.DataFrame
        :param significance_level: Significance level (user defined)
        :type significance_level: float
        :return: Significant parameters
        :rtype: list
        """
        significant_df = sensitivity_data[sensitivity_data["S1"].ge(significance_level)]
        return significant_df.index.map(str)

    def find_influential_conv_parameters(self):
        """Find the input parameters with a value above the significance level.
        self.sumsq_sensitivity_df = pd.DataFrame()

        :return: _description_
        :rtype: _type_
        """
        significant_df = self.sumsq_sensitivity_df[
            self.sumsq_sensitivity_df["converged"].ge(self.significance_level)
        ]
        return significant_df.index.map(str)

    def find_most_sensitive_interaction(self, parameters_to_search):
        """Perform a sensitivity analysis on a list of variables, store if influential."""
        for item in parameters_to_search:
            self.figure_of_merit = item
            # The error here is caused by paramters which have the same value in every run.
            # For example q0 which is 1.0 for every run.
            self.calculate_sensitivity(item)
            if self.sensitivity_df.isnull().values.any() & self.sensitivity_df.empty:
                pass
            else:
                new_df = self.sensitivity_df[
                    self.sensitivity_df["S1"] > self.significance_level
                ]
                if new_df.empty:
                    pass
                else:
                    if len(new_df) == self.number_sampled_vars:
                        pass
                    else:
                        new_df.insert(2, "variable", item)
                        self.item_with_sensitivity.append(item)
                        self.sens_dict[item] = new_df.to_dict()
                        self.list_new_dfs.append(new_df)
        print(self.item_with_sensitivity)
        mc_run_df = pd.concat(self.list_new_dfs)
        print(mc_run_df)
        for name in self.input_names:
            if name in mc_run_df.index:
                name_df = mc_run_df.loc[[name]]
                print("Number influenced by ", name, len(name_df))
                print(name_df)
                name_df.to_json(
                    self.path_in
                    + name
                    + "_influence_"
                    + str(self.figure_of_merit)
                    + ".json",
                    orient="split",
                    compression="infer",
                    index="true",
                )

        # print(mc_run_df[mc_run.index.duplicated(keep=False)])

        with open("result_05si.json", "w") as fp:
            json.dump(self.sens_dict, fp)

    def calculate_reliability(self):
        """Calculate the Reliability Index, defined as number of converged runs divided by
        the number of failed runs.
        """
        self.reliability_index = round(
            (len(self.converged_df.index)) / len((self.uncertainties_df.index)), 2
        )

    def read_json(self, file):
        """Read and print a json file.

        :param file: Path to json
        :type file: str
        """
        df = pd.read_json(file, orient="split")
        print(df)  # .get("pdivfraction"))  # .values[0].keys())

    def plot_rbd_si_indices(self):
        """Calculate RBD FAST Sobol Indices and plot"""
        fig, ax = plt.subplots(1)
        # fig.set_size_inches(18, 5)

        ax.tick_params(labelsize=14)
        sensdf = self.sensitivity_df
        sensdf = sensdf.rename(self.name_dict, axis=0)
        # x-axis
        sensdf.plot(
            kind="barh", y="S1", xerr="S1_conf", ax=ax, align="center", capsize=3
        )
        # y-axis
        ax.set_xlabel("Sensitivity Indices: " + "Fusion power", fontsize=20)
        ax.set_ylabel("PROCESS parameter", fontsize=20)

        # in striped grey : not significant indices
        ax.fill_betweenx(
            y=[-0.5, len(self.sensitivity_df)],
            x1=0,
            x2=self.significance_level,
            color="grey",
            alpha=0.2,
            hatch="//",
            edgecolor="white",
        )
        plt.grid()
        plt.savefig("plots/sensitivity_fom.svg", bbox_inches="tight")
        plt.show()

    def plot_sumsq_sensitivity(self):
        """Find the input paramters influencing whether PROCESS converges."""
        fig, ax = plt.subplots(1)
        ax.tick_params(labelsize=16)
        sumsqnamedf = self.sumsq_sensitivity_df.rename(self.name_dict, axis=0)
        # x-axis
        sumsqnamedf.plot(kind="barh", y="converged", ax=ax, align="center", capsize=3)
        # y-axis
        ax.set_xlabel("Influence on convergance", fontsize=20)
        ax.set_ylabel("PROCESS parameter", fontsize=20)
        ax.fill_betweenx(
            y=[-0.5, len(self.sensitivity_df)],
            x1=0,
            x2=self.significance_level,
            color="grey",
            alpha=0.2,
            hatch="//",
            edgecolor="white",
            label="Not significant",
        )
        ax.legend(fontsize=12, borderpad=0.01, ncol=1)

        plt.grid()
        plt.savefig("plots/rds_indices.svg", bbox_inches="tight")
        plt.show()

    def convergence_study(self, n, sampled_inputs, process_output):
        """This function is used to calculate RBD FAST sensitivities indices for a subset.
        It draws a random sample, of a given size, from the total dataset. This is used to
        create a study of the convergence of the indices for the input parameters.

        :param n: Number of samples to select
        :type n: int
        :param sampled_inputs: Array of the sampled input parameters
        :type sampled_inputs: numpy.array()
        :param process_output: Array of the figure of merit
        :type process_output: numpy.array()
        :return: Sensitivity Indices dict, length of the subset
        :rtype: dict, int
        """
        subset = np.random.choice(len(process_output), size=n, replace=False)
        self.problem = {
            "num_vars": self.number_sampled_vars,
            "names": self.input_names,
        }
        rbd_results = rbd_fast.analyze(
            self.problem, X=sampled_inputs[subset], Y=process_output[subset]
        )
        return rbd_results, len(subset)

    def find_convergence_indices(self, step_size, figure_of_merit):
        """Calculate arrays of sensitivity indices for each input parameter to the figure of merit.
        Use the output array of indices to check for convergence (are they stable around a point)

        :param n: Calculate indices `n` times along the number of samples
        :type n: int
        :return: Array of arrays containing sensitivity indices for each sampled input
        :rtype: numpy.array()
        """
        converged_figure_of_merit_df = self.get_fom_converged_df(figure_of_merit)
        sampled_vars_df = self.converged_sampled_vars_df

        indices_df_list = []
        for n in np.arange(
            start=50, stop=len(converged_figure_of_merit_df) + 1, step=step_size
        ):
            rbd_results, len_subset = self.convergence_study(
                n=n,
                sampled_inputs=sampled_vars_df.to_numpy(),
                process_output=converged_figure_of_merit_df.to_numpy(),
            )
            rbd_results["samples"] = len_subset
            indices_df_list.append(rbd_results.to_df())
        indices_df = pd.concat(indices_df_list)

        return indices_df

    def plot_convergence_study(self, step_size, figure_of_merit):
        """Performs a convergence study and plots the results. Plots confidence levels only if
        final indices are greater than significance level.

        :param step_size: Number of samples to increment by when calculating sensitivity indices
        :type step_size: int
        """
        indices_df = self.find_convergence_indices(step_size, figure_of_merit)
        conv_fig, conv_ax = plt.subplots()
        conv_fig.set_size_inches(15, 8)

        conv_ax.set_title(
            "Convergence of the Sensitivity Indices for FOM: " + "Fusion power",
            fontsize=20,
        )

        conv_ax.set_ylim(ymin=-0.15, ymax=1)

        for name in self.input_names:
            name_df = indices_df.loc[[name]]
            name_df.plot(
                x="samples",
                y="S1",
                label=self.name_dict[name],
                ax=conv_ax,
                linewidth=3,
            )

            if name_df["S1"].iloc[-1] > self.significance_level:
                plt.fill_between(
                    name_df["samples"],
                    (name_df["S1"] - name_df["S1_conf"]),
                    (name_df["S1"] + name_df["S1_conf"]),
                    # color="tab:blue",
                    alpha=0.1,
                )

        conv_ax.fill_between(
            x=[indices_df["samples"].min(), indices_df["samples"].max()],
            y1=-0.15,
            y2=self.significance_level,
            color="grey",
            alpha=0.3,
            hatch="//",
            edgecolor="white",
            label="Not significant",
        )
        conv_ax.tick_params(labelsize=16)
        conv_ax.set_xlim(
            xmin=indices_df["samples"].min(), xmax=indices_df["samples"].max()
        )
        conv_ax.set_xlabel("Samples", fontdict={"fontsize": 20})
        conv_ax.set_ylabel("Sensitivity Index", fontdict={"fontsize": 20})
        conv_ax.legend(fontsize=12, borderpad=0.01, ncol=2)

        plt.grid()
        plt.savefig("plots/convergence.svg", bbox_inches="tight")
        plt.show()

    def merge_hdf_files(self):
        """Looks for uncertainty hdf files in the working folder and merges them into a
        single dataframe for analysis.

        :return: Uncertainties DataFrame, Parameters DataFrame
        :rtype: pandas.DataFrame, pandas.Dataframe
        """
        list_uncertainties_dfs = []
        list_params_dfs = []
        for root, dirs, files in os.walk(self.path_in):
            for file in files:
                pos_hdf = root + os.sep + file
                if pos_hdf.endswith(".h5") and "uncertainties_data" in pos_hdf:
                    extra_uncertainties_df = pd.read_hdf(pos_hdf)
                    list_uncertainties_dfs.append(extra_uncertainties_df)
                if pos_hdf.endswith(".h5") and "param_values" in pos_hdf:
                    extra_params_df = pd.read_hdf(pos_hdf)
                    list_params_dfs.append(extra_params_df)
        return pd.concat(list_uncertainties_dfs), pd.concat(list_params_dfs)

    def create_scatter_plot(
        self,
        axes=None,
        df=None,
        diagonal="hist",
        density_kwds=None,
        hist_kwds=None,
        marker="*",
        alpha=0.5,
        color="tab:blue",
        **kwds,
    ):
        """
        Create a scatter plot to visuals inputs against the figure of merit. Used to look for trends in the data at a glance.
        diagonal: either "hist", "kde" or "density"
        See def scatter_matrix in: https://github.com/pandas-dev/pandas/blob/526f40431a51e1b1621c30a4d74df9006e0274b8/pandas/plotting/_matplotlib/misc.py

        """
        range_padding = 0.05
        hist_kwds = hist_kwds or {}
        density_kwds = density_kwds or {}

        ## fix input data
        mask = pd.notna(df)

        boundaries_list = []
        for a in df.columns:
            values = df[a].values[mask[a].values]
            rmin_, rmax_ = np.min(values), np.max(values)
            rdelta_ext = (rmax_ - rmin_) * range_padding / 2.0
            boundaries_list.append((rmin_ - rdelta_ext, rmax_ + rdelta_ext))

        ## iterate over columns
        for i, a in enumerate(df.columns):
            for j, b in enumerate(df.columns):
                ax = axes[i, j]  ## to abbreviate the code

                if i == j:
                    values = df[a].values[mask[a].values]

                    # Deal with the diagonal by drawing a histogram there.
                    if diagonal == "hist":
                        ax.hist(values, color=color, alpha=alpha, **hist_kwds)

                    elif diagonal in ("kde", "density"):
                        from scipy.stats import gaussian_kde

                        y = values
                        gkde = gaussian_kde(y)
                        ind = np.linspace(y.min(), y.max(), 1000)
                        ax.plot(ind, gkde.evaluate(ind), color=color, **density_kwds)

                    ax.set_xlim(boundaries_list[i])
                    # Take the text off inbetween diagonal plots.
                    if i < 4:
                        ax.get_xaxis().set_visible(False)

                else:
                    common = (mask[a] & mask[b]).values
                    h = ax.hist2d(
                        x=df[b][common], y=df[a][common], bins=10, cmap="magma"
                    )
                    if i > j:
                        plt.colorbar(h[3], ax=ax)
                    plt.grid()
                    ax.set_xlim(boundaries_list[j])
                    ax.set_ylim(boundaries_list[i])
                ax.tick_params(axis="both", which="major", labelsize=12)
                ax.set_xlabel(self.name_dict[b], fontsize=12)
                ax.set_ylabel(self.name_dict[a], fontsize=12)

                if i < j:
                    axes[i, j].set_visible(False)
                if j != 0:
                    ax.yaxis.set_visible(False)

        return

    def plot_scatter_plot(self, plot_unconverged=False):
        """Configures the plot the scatter plot.
        :param plot_unconverged: Plot unconverged runs in a red histogram, defaults to False
        :type plot_unconverged: bool, optional
        """

        figsize = (17, 17)
        figure(figsize=figsize)
        axes = {}
        n = self.plot_converged_df.columns.size
        gs = mpl.gridspec.GridSpec(
            n,
            n,
            left=0.12,
            right=0.97,
            bottom=0.12,
            top=0.97,
            wspace=0.025,
            hspace=0.005,
        )
        for i, a in enumerate(self.plot_converged_df.columns):
            for j, b in enumerate(self.plot_converged_df.columns):
                axes[i, j] = plt.subplot(gs[i, j])
        self.create_scatter_plot(
            axes=axes,
            df=self.plot_converged_df,
            diagonal="hist",
            color="#0000ff",
            marker=".",
        )

        if plot_unconverged == True:
            self.create_scatter_plot(
                axes,
                diagonal="hist",
                df=self.plot_unconverged_df,
                color="#ff0000",
                marker=".",
            )
        plt.savefig("plots/2dhist.svg", bbox_inches="tight")
        plt.show()

    def regional_sensitivity_analysis(self):
        """Regional sensitivity anlysis to find the parameters which influence a converged solution.
        Uses modified RSA technique from SALib."""
        # x = inputs
        # y = figure of merit (sqsumsq)
        figure_of_merit_df = self.uncertainties_df["sqsumsq"]
        sampled_vars_df = self.uncertainties_df[self.input_names]

        self.problem = {
            "num_vars": self.number_sampled_vars,
            "names": self.input_names,
        }
        rsa_result = rsa.analyze(
            problem=self.problem,
            X=sampled_vars_df.to_numpy(),
            Y=figure_of_merit_df.to_numpy(),
            bins=2,
            target="convergence",
            print_to_console=False,
            seed=1,
        )
        self.sumsq_sensitivity_df = rsa_result.to_df()
        self.sumsq_sensitivity_df = self.sumsq_sensitivity_df.T
        self.sumsq_sensitivity_df.columns = ["converged", "unconverged"]
        self.sumsq_sensitivity_df = self.sumsq_sensitivity_df.sort_values(
            by="converged", ascending=False
        )
        # print(self.sumsq_sensitivity_df)

    def hdmr_analysis(self):
        """HDMR Analysis - not well explored."""
        fom_df = self.get_fom_converged_df("powfmw")
        hdmr_res = hdmr.analyze(
            self.problem, self.converged_sampled_vars_df.to_numpy(), fom_df.to_numpy()
        )
        print(hdmr_res.to_df().to_string())

    def ecdf_plot(self, figure_of_merit):
        """Plot Empirical Cumulative Distribution Functions for converged and unconverged samples.
        Additionally, plot the convergence rate and the design point value.

        :param figure_of_merit: Parameter to investigate
        :type figure_of_merit: str
        """
        fig, ax1 = plt.subplots(1)
        plt.rcParams["axes.xmargin"] = 0
        plt.rcParams["axes.ymargin"] = 0
        ax2 = ax1.twinx()
        plt.style.library["tableau-colorblind10"]
        # Arrange the data into df
        converged_fom_df = self.converged_df[figure_of_merit].to_numpy()
        uq_fom_df = self.uncertainties_df[figure_of_merit].to_numpy()
        # Calculate cumulative distribution functions
        ecdf_unconv = sm.distributions.ECDF(uq_fom_df)
        ecdf_conv = sm.distributions.ECDF(converged_fom_df)
        # Bin the data
        x = np.linspace(min(uq_fom_df), max(uq_fom_df))
        ri_t = []
        y_unconv = ecdf_unconv(x)
        y_conv = ecdf_conv(x)
        # Plot the ecdf functions
        ax1.step(x, y_unconv, color="tab:red", label="Unconverged samples")
        ax1.step(x, y_conv, color="tab:blue", label="Converged samples")
        ax1.set_ylabel("Percentage of Samples", fontsize=20)
        ax1.set_xlabel(self.name_dict[figure_of_merit], fontsize=20)
        # Calculate rate of convergence for bins in x
        for d in range(len(x) - 1):

            n_c = len(
                self.converged_df[
                    self.converged_df[figure_of_merit].between(x[d], x[d + 1])
                ].index
            )
            n_t = len(
                self.uncertainties_df[
                    self.uncertainties_df[figure_of_merit].between(x[d], x[d + 1])
                ].index
            )
            if n_t == 0:
                n_t = 0.0000001
            ri = n_c / n_t
            ri_t.append(ri)
        # Finds the edges of bins (must be a better way to do this section)
        h, edges = np.histogram(ri_t, bins=x)
        # Plot convergence rate
        ax1.stairs(ri_t, edges, color="tab:orange", label="Convergence rate")
        ax2.set_ylabel("Convergence Rate", fontsize=20)
        # Plot design point
        ypoint = 0.5
        if figure_of_merit == "kappa":
            ypoint = 1.0
        ### copy curve line y coords and set to a constant
        lines = y_unconv.copy()
        lines[:] = ypoint

        # get intersection
        first_line = LineString(np.column_stack((x, y_unconv)))
        second_line = LineString(np.column_stack((x, lines)))
        intersection = first_line.intersection(second_line)
        ax1.plot(*intersection.xy, "s", color="tab:green", markersize=7)

        # plot hline and vline
        ax1.hlines(
            y=ypoint,
            xmin=min(x),
            xmax=intersection.x,
            label="Design point value",
            clip_on=False,
            color="tab:green",
            linestyles="dashed",
            alpha=0.7,
        )
        ax1.vlines(
            x=intersection.x,
            ymin=0,
            ymax=intersection.y,
            clip_on=True,
            color="tab:green",
            linestyles="dashed",
            alpha=0.7,
        )

        ax1.legend()
        ax1.grid(axis="both")
        plt.savefig("plots/" + figure_of_merit + "ecdf.svg", bbox_inches="tight")


def main(args=None):

    args = parse_args(args)
    uq_data = UncertaintyData(args.uq_folder, args.figure_of_merit)
    uq_data.calculate_sensitivity()
    uq_data.calculate_reliability()
    print(uq_data.sensitivity_df)
    print(" ")
    print(
        "Number of converged runs: ",
        uq_data.number_of_converged_runs,
    )
    print(
        "Reliability index: ",
        uq_data.reliability_index,
    )
    print(" ")
    uq_data.configure_data_for_plotting()
    uq_data.plot_scatter_plot(plot_unconverged=False)


def unique_cols(df):
    """Find columns which are not all the same value.

    :param df: DataFrame in question
    :type df: pandas.DataFrame
    :return: DataFrame containing with no constant columns
    :rtype: pandas.DataFrame
    """
    a = df.to_numpy()
    return (a[0] == a).all(0)


if __name__ == "__main__":
    main()
