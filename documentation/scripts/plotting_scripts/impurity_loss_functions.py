"""Plotting script for the total impurity radiation loss function (L_z) profiles."""

from importlib import resources
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from process.core.io.plot.summary import read_imprad_data


def plot_line_brem_loss_function_profile(
    axis: plt.Axes,
    impp: str,
):
    """Function to plot Line and Bremsstrahlung loss function (L_z) profile.

    Parameters
    ----------
    axis : plt.Axes
        axis object to add plot to
    impp : str
        impurity path

    """
    # read in the impurity data
    imp_data = read_imprad_data(_skiprows=2, data_path=impp)

    impurity_labels = [
        "H",
        "He",
        "Be",
        "C",
        "N",
        "O",
        "Ne",
        "Si",
        "Ar",
        "Fe",
        "Ni",
        "Kr",
        "Xe",
        "W",
    ]

    for label, raw_species_data in zip(impurity_labels, imp_data, strict=True):
        species_data = np.asarray(raw_species_data)
        axis.plot(species_data[:, 0], species_data[:, 1], label=label)

    axis.legend(loc="best", ncol=4)
    axis.minorticks_on()

    axis.set_xlabel(r"$T_e$ [eV]")
    axis.set_ylabel(r"$L_z$ $[\mathrm{W}\mathrm{m}^3]$")
    axis.set_title("Line & Bremsstrahlung Loss Function ($L_z$) vs Temperature")
    axis.set_xlim([
        min(np.asarray(species_data)[:, 0].min() for species_data in imp_data),
        max(np.asarray(species_data)[:, 0].max() for species_data in imp_data),
    ])
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.minorticks_on()
    axis.xaxis.grid(True, which="both", alpha=0.2)
    axis.yaxis.grid(True, which="both", alpha=0.2)
    plt.savefig("adas_radiation.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    with resources.path(
        "process.data.lz_non_corona_14_elements", "Ar_lz_tau.dat"
    ) as imp_path:
        data_folder = str(imp_path.parent) + "/"

    if Path(data_folder).is_dir():
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_line_brem_loss_function_profile(ax, impp=data_folder)
    else:
        print(
            "\033[91m Warning : Impossible to recover impurity data, try running the "
            "macro in the main/utility folder"
        )
        print("          -> No impurity plot done\033[0m")
