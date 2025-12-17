"""
Code generating the TF coil inboard mid-plane stress/strain summary plots
The whole radial distribution is displayed

Author: S. Kahn (sebastien.kahn@ukaea.uk)

Input file:
SIG_TF.json
"""

import argparse
import json
import os
from argparse import RawTextHelpFormatter
from pathlib import Path

import matplotlib.pyplot as plt


def main(args=None):
    # PARSING USER PARAMETERS
    # please execute 'python plot_stress_tf.py -h' for input information
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
        help=(
            "Plot selection string :\n - If it containts 'sig'      -> Stress radial dependency \n#"
            "  - If it containts 'strain'   -> Strain \n - If it containts 'disp'     -> Displacement \n"
            " - If it containts 'all'      -> all the mentioned plots (default value)"
        ),
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
    parser.add_argument(
        "-out",
        "--term_output",
        action="store_true",
        help="Option to show stress on terminal output",
    )
    parser.add_argument(
        "-f",
        "--input_file",
        default="SIG_TF.json",
        help="specify input file path (default = SIG_TF.json)",
    )

    # Option argument extraction
    # --------------------------
    args = parser.parse_args(args)
    plot_selection = str(args.plot_selec)
    save_format = str(args.save_format)
    axis_font_size = int(args.axis_font_size)
    term_output = args.term_output

    # Boolean swiches for plot selection
    # -----------------------------------
    plot_sig = ("sig" in plot_selection) or plot_selection == "all"
    plot_disp = ("disp" in plot_selection) or plot_selection == "all"
    plot_strain = ("strain" in plot_selection) or plot_selection == "all"
    plot_sm_sig = ("sm_sig" in plot_selection) or plot_selection == "all"

    # Step 1 : Data extraction
    # ----------------------------------------------------------------------------------------------
    # Number of physical quantity value per coil layer
    n_radial_array_layer = 0

    # Physical quantities : full vectors
    radius = []
    radial_smeared_stress = []
    toroidal_smeared_stress = []
    vertical_smeared_stress = []
    tresca_smeared_stress = []
    radial_stress = []
    toroidal_stress = []
    vertical_stress = []
    vm_stress = []
    tresca_stress = []
    cea_tresca_stress = []
    radial_strain = []
    toroidal_strain = []
    vertical_strain = []
    radial_displacement = []

    # Physical quantity : WP stress
    wp_vertical_stress = []

    # Physical quantity : values at layer border
    bound_radius = []
    bound_radial_smeared_stress = []
    bound_toroidal_smeared_stress = []
    bound_vertical_smeared_stress = []
    bound_tresca_smeared_stress = []
    bound_radial_stress = []
    bound_toroidal_stress = []
    bound_vertical_stress = []
    bound_vm_stress = []
    bound_tresca_stress = []
    bound_cea_tresca_stress = []
    bound_radial_strain = []
    bound_toroidal_strain = []
    bound_vertical_strain = []
    bound_radial_displacement = []

    with open(args.input_file) as f:
        sig_file_data = json.load(f)

    # Getting the data to be plotted
    n_radial_array_layer = sig_file_data["Points per layers"]
    n_points = len(sig_file_data["Radius (m)"])
    n_layers = int(n_points / n_radial_array_layer)
    for ii in range(n_layers):
        # Full vector
        radius.append([])
        radial_stress.append([])
        toroidal_stress.append([])
        vertical_stress.append([])
        radial_smeared_stress.append([])
        toroidal_smeared_stress.append([])
        vertical_smeared_stress.append([])
        vm_stress.append([])
        tresca_stress.append([])
        cea_tresca_stress.append([])
        radial_displacement.append([])

        for jj in range(n_radial_array_layer):
            radius[ii].append(
                sig_file_data["Radius (m)"][ii * n_radial_array_layer + jj]
            )
            radial_stress[ii].append(
                sig_file_data["Radial stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            toroidal_stress[ii].append(
                sig_file_data["Toroidal stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            if len(sig_file_data["Vertical stress (MPa)"]) == 1:
                vertical_stress[ii].append(sig_file_data["Vertical stress (MPa)"][0])
            else:
                vertical_stress[ii].append(
                    sig_file_data["Vertical stress (MPa)"][
                        ii * n_radial_array_layer + jj
                    ]
                )
            radial_smeared_stress[ii].append(
                sig_file_data["Radial smear stress (MPa)"][
                    ii * n_radial_array_layer + jj
                ]
            )
            toroidal_smeared_stress[ii].append(
                sig_file_data["Toroidal smear stress (MPa)"][
                    ii * n_radial_array_layer + jj
                ]
            )
            vertical_smeared_stress[ii].append(
                sig_file_data["Vertical smear stress (MPa)"][
                    ii * n_radial_array_layer + jj
                ]
            )
            vm_stress[ii].append(
                sig_file_data["Von-Mises stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            tresca_stress[ii].append(
                sig_file_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            cea_tresca_stress[ii].append(
                sig_file_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer + jj]
            )
            radial_displacement[ii].append(
                sig_file_data["rad. displacement (mm)"][ii * n_radial_array_layer + jj]
            )

        # Layer lower boundaries values
        bound_radius.append(sig_file_data["Radius (m)"][ii * n_radial_array_layer])
        bound_radial_stress.append(
            sig_file_data["Radial stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_toroidal_stress.append(
            sig_file_data["Toroidal stress (MPa)"][ii * n_radial_array_layer]
        )
        if len(sig_file_data["Vertical stress (MPa)"]) == 1:
            bound_vertical_stress.append(sig_file_data["Vertical stress (MPa)"][0])
        else:
            bound_vertical_stress.append(
                sig_file_data["Vertical stress (MPa)"][ii * n_radial_array_layer]
            )
        bound_radial_smeared_stress.append(
            sig_file_data["Radial smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_toroidal_smeared_stress.append(
            sig_file_data["Toroidal smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_vertical_smeared_stress.append(
            sig_file_data["Vertical smear stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_vm_stress.append(
            sig_file_data["Von-Mises stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_tresca_stress.append(
            sig_file_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_cea_tresca_stress.append(
            sig_file_data["CEA Tresca stress (MPa)"][ii * n_radial_array_layer]
        )
        bound_radial_displacement.append(
            sig_file_data["rad. displacement (mm)"][ii * n_radial_array_layer]
        )

        # Layer upper boundaries values
        bound_radius.append(
            sig_file_data["Radius (m)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_radial_stress.append(
            sig_file_data["Radial stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_toroidal_stress.append(
            sig_file_data["Toroidal stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        if len(sig_file_data["Vertical stress (MPa)"]) == 1:
            bound_vertical_stress.append(sig_file_data["Vertical stress (MPa)"][0])
        else:
            bound_vertical_stress.append(
                sig_file_data["Vertical stress (MPa)"][
                    (ii + 1) * n_radial_array_layer - 1
                ]
            )
        bound_radial_smeared_stress.append(
            sig_file_data["Radial smear stress (MPa)"][
                (ii + 1) * n_radial_array_layer - 1
            ]
        )
        bound_toroidal_smeared_stress.append(
            sig_file_data["Toroidal smear stress (MPa)"][
                (ii + 1) * n_radial_array_layer - 1
            ]
        )
        bound_vertical_smeared_stress.append(
            sig_file_data["Vertical smear stress (MPa)"][
                (ii + 1) * n_radial_array_layer - 1
            ]
        )
        bound_vm_stress.append(
            sig_file_data["Von-Mises stress (MPa)"][(ii + 1) * n_radial_array_layer - 1]
        )
        bound_tresca_stress.append(
            sig_file_data["CEA Tresca stress (MPa)"][
                (ii + 1) * n_radial_array_layer - 1
            ]
        )
        bound_cea_tresca_stress.append(
            sig_file_data["CEA Tresca stress (MPa)"][
                (ii + 1) * n_radial_array_layer - 1
            ]
        )
        bound_radial_displacement.append(
            sig_file_data["rad. displacement (mm)"][(ii + 1) * n_radial_array_layer - 1]
        )

    # TRESCA smeared stress [MPa]
    for ii in range(n_layers):
        tresca_smeared_stress.append([])

        bound_tresca_smeared_stress.append(
            max(abs(radial_smeared_stress[ii][0]), abs(toroidal_smeared_stress[ii][0]))
            + vertical_smeared_stress[ii][0]
        )
        bound_tresca_smeared_stress.append(
            max(
                abs(radial_smeared_stress[ii][n_radial_array_layer - 1]),
                abs(toroidal_smeared_stress[ii][n_radial_array_layer - 1]),
            )
            + vertical_smeared_stress[ii][n_radial_array_layer - 1]
        )
        for jj in range(n_radial_array_layer):
            tresca_smeared_stress[ii].append(
                max(
                    abs(radial_smeared_stress[ii][jj]),
                    abs(toroidal_smeared_stress[ii][jj]),
                )
                + vertical_smeared_stress[ii][jj]
            )

    # Strains
    if len(sig_file_data) > 16:
        for ii in range(n_layers):
            radial_strain.append([])
            toroidal_strain.append([])
            vertical_strain.append([])

            bound_radial_strain.append(
                sig_file_data["Radial strain"][ii * n_radial_array_layer]
            )
            bound_toroidal_strain.append(
                sig_file_data["Toroidal strain"][ii * n_radial_array_layer]
            )
            bound_vertical_strain.append(
                sig_file_data["Vertical strain"][ii * n_radial_array_layer]
            )
            bound_radial_strain.append(
                sig_file_data["Radial strain"][(ii + 1) * n_radial_array_layer - 1]
            )
            bound_toroidal_strain.append(
                sig_file_data["Toroidal strain"][(ii + 1) * n_radial_array_layer - 1]
            )
            bound_vertical_strain.append(
                sig_file_data["Vertical strain"][(ii + 1) * n_radial_array_layer - 1]
            )
            for jj in range(n_radial_array_layer):
                radial_strain[ii].append(
                    sig_file_data["Radial strain"][ii * n_radial_array_layer + jj]
                )
                toroidal_strain[ii].append(
                    sig_file_data["Toroidal strain"][ii * n_radial_array_layer + jj]
                )
                vertical_strain[ii].append(
                    sig_file_data["Vertical strain"][ii * n_radial_array_layer + jj]
                )

                if "WP smeared stress (MPa)" in sig_file_data:
                    wp_vertical_stress.append(
                        sig_file_data["WP smeared stress (MPa)"][jj]
                    )

    # Terminal output
    # ---------------
    if term_output:
        ii_ins = 0
        ii_mids = int(0.5 * float(n_radial_array_layer))
        ii_outs = n_radial_array_layer - 1

        print()
        print()
        print("Layer stress details")
        print("____________________")

        for ii in range(n_layers):
            print(f"Layer {ii + 1}")
            print("------------------------------")
            print(
                "steel radial   stress in the inner/middle/out point:"
                f" {radial_stress[ii][ii_ins]}/{radial_stress[ii][ii_mids]}/{radial_stress[ii][ii_outs]} MPa"
            )
            print(
                "steel toroidal stress in the inner/middle/out point:"
                f" {toroidal_stress[ii][ii_ins]}/{toroidal_stress[ii][ii_mids]}/{toroidal_stress[ii][ii_outs]} MPa"
            )
            print(
                "steel vertical stress in the inner/middle/out point:"
                f" {vertical_stress[ii][ii_ins]}/{vertical_stress[ii][ii_mids]}/{vertical_stress[ii][ii_outs]} MPa"
            )
            print(
                "steel TRESCA   stress in the inner/middle/out point:"
                f" {tresca_stress[ii][ii_ins]}/{tresca_stress[ii][ii_mids]}/{tresca_stress[ii][ii_outs]} MPa"
            )
            print()
            print(
                "smeared radial   stress in the inner/middle/out point :"
                f" {radial_smeared_stress[ii][ii_ins]}/{radial_smeared_stress[ii][ii_mids]}/{radial_smeared_stress[ii][ii_outs]} MPa"
            )
            print(
                "smeared toroidal stress in the inner/middle/out point :"
                f" {toroidal_smeared_stress[ii][ii_ins]}/{toroidal_smeared_stress[ii][ii_mids]}/{toroidal_smeared_stress[ii][ii_outs]} MPa"
            )
            print(
                "smeared vertical stress in the inner/middle/out point :"
                f" {vertical_smeared_stress[ii][ii_ins]}/{vertical_smeared_stress[ii][ii_mids]}/{vertical_smeared_stress[ii][ii_outs]} MPa"
            )
            print(
                "smeared TRESCA   stress in the inner/middle/out point :"
                f" {tresca_smeared_stress[ii][ii_ins]}/{tresca_smeared_stress[ii][ii_mids]}/{tresca_smeared_stress[ii][ii_outs]} MPa"
            )
            print()

            if len(sig_file_data) > 16:
                print(
                    "radial   strain in the inner/middle/out point :"
                    f" {radial_strain[ii][ii_ins]}/{radial_strain[ii][ii_mids]}/{radial_strain[ii][ii_outs]}"
                )
                print(
                    "toroidal strain in the inner/middle/out point :"
                    f" {toroidal_strain[ii][ii_ins]}/{toroidal_strain[ii][ii_mids]}/{toroidal_strain[ii][ii_outs]}"
                )
                print(f"vertical strain : {vertical_strain[ii][0]}")
                print()

        if len(wp_vertical_stress) != 0:
            print(
                "smeared WP vertical stress in the inner/middle/out point :"
                f" {wp_vertical_stress[0]}/{wp_vertical_stress[ii_mids]}/{wp_vertical_stress[ii_outs]} MPa"
            )
        print()

    outdir = Path.cwd()
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    axis_tick_size = 16
    legend_size = 12
    mark_size = 13
    line_width = 3.5

    # PLOT 1 : Stress summary
    # ------------------------
    if plot_sig:
        for ii in range(n_layers):
            plt.plot(
                radius[ii],
                radial_stress[ii],
                "-",
                linewidth=line_width,
                color="lightblue",
            )
            plt.plot(
                radius[ii],
                toroidal_stress[ii],
                "-",
                linewidth=line_width,
                color="wheat",
            )
            plt.plot(
                radius[ii],
                vertical_stress[ii],
                "-",
                linewidth=line_width,
                color="lightgrey",
            )
            plt.plot(
                radius[ii], tresca_stress[ii], "-", linewidth=line_width, color="pink"
            )
            plt.plot(
                radius[ii], vm_stress[ii], "-", linewidth=line_width, color="violet"
            )
        plt.plot(
            radius[0],
            radial_stress[0],
            "--",
            color="dodgerblue",
            label=r"$\sigma_{rr}$",
        )
        plt.plot(
            radius[0],
            toroidal_stress[0],
            "--",
            color="orange",
            label=r"$\sigma_{\theta\theta}$",
        )
        plt.plot(
            radius[0],
            vertical_stress[0],
            "--",
            color="mediumseagreen",
            label=r"$\sigma_{zz}$",
        )
        plt.plot(
            radius[0],
            tresca_stress[0],
            "-",
            color="crimson",
            label=r"$\sigma_{TRESCA}$",
        )
        plt.plot(
            radius[0],
            vm_stress[0],
            "-",
            color="darkviolet",
            label=r"$\sigma_{Von\ mises}$",
        )
        for ii in range(1, n_layers):
            plt.plot(radius[ii], radial_stress[ii], "--", color="dodgerblue")
            plt.plot(radius[ii], toroidal_stress[ii], "--", color="orange")
            plt.plot(radius[ii], vertical_stress[ii], "--", color="mediumseagreen")
            plt.plot(radius[ii], tresca_stress[ii], "-", color="crimson")
            plt.plot(radius[ii], vm_stress[ii], "-", color="darkviolet")
        plt.plot(
            bound_radius,
            bound_radial_stress,
            "|",
            markersize=mark_size,
            color="dodgerblue",
        )
        plt.plot(
            bound_radius,
            bound_toroidal_stress,
            "|",
            markersize=mark_size,
            color="orange",
        )
        plt.plot(
            bound_radius,
            bound_vertical_stress,
            "|",
            markersize=mark_size,
            color="mediumseagreen",
        )
        plt.plot(
            bound_radius,
            bound_tresca_stress,
            "|",
            markersize=mark_size,
            color="crimson",
        )
        plt.plot(
            bound_radius, bound_vm_stress, "|", markersize=mark_size, color="darkviolet"
        )
        plt.grid(True)
        plt.ylabel(r"$\sigma$ [$MPa$]", fontsize=axis_font_size)
        plt.xlabel(r"$R$ [$m$]", fontsize=axis_font_size)
        plt.legend(loc="best", fontsize=legend_size)
        plt.xticks(size=axis_tick_size)
        plt.yticks(size=axis_tick_size)
        plt.tight_layout()
        plt.savefig(f"{outdir}/structure_stress.{save_format}")
        plt.clf()
        plt.cla()

    # PLOT 2 : Smeared stress summary
    # ------------------------
    if plot_sm_sig:
        for ii in range(n_layers):
            plt.plot(
                radius[ii],
                radial_smeared_stress[ii],
                "-",
                linewidth=line_width,
                color="lightblue",
            )
            plt.plot(
                radius[ii],
                toroidal_smeared_stress[ii],
                "-",
                linewidth=line_width,
                color="wheat",
            )
            plt.plot(
                radius[ii],
                vertical_smeared_stress[ii],
                "-",
                linewidth=line_width,
                color="lightgrey",
            )
            plt.plot(
                radius[ii],
                tresca_smeared_stress[ii],
                "-",
                linewidth=line_width,
                color="pink",
            )
        plt.plot(
            radius[0],
            radial_smeared_stress[0],
            "--",
            color="dodgerblue",
            label=r"$\sigma_{rr}^\mathrm{smeared}$",
        )
        plt.plot(
            radius[0],
            toroidal_smeared_stress[0],
            "--",
            color="orange",
            label=r"$\sigma_{\theta\theta}^\mathrm{smeared}$",
        )
        plt.plot(
            radius[0],
            vertical_smeared_stress[0],
            "--",
            color="mediumseagreen",
            label=r"$\sigma_{zz}^\mathrm{smeared}$",
        )
        plt.plot(
            radius[0],
            tresca_smeared_stress[0],
            "-",
            color="crimson",
            label=r"$\sigma_{TRESCA}^\mathrm{smeared}$",
        )
        for ii in range(1, n_layers):
            plt.plot(radius[ii], radial_smeared_stress[ii], "--", color="dodgerblue")
            plt.plot(radius[ii], toroidal_smeared_stress[ii], "--", color="orange")
            plt.plot(
                radius[ii], vertical_smeared_stress[ii], "--", color="mediumseagreen"
            )
            plt.plot(radius[ii], tresca_smeared_stress[ii], "-", color="crimson")
        plt.plot(
            bound_radius,
            bound_radial_smeared_stress,
            "|",
            markersize=mark_size,
            color="dodgerblue",
        )
        plt.plot(
            bound_radius,
            bound_toroidal_smeared_stress,
            "|",
            markersize=mark_size,
            color="orange",
        )
        plt.plot(
            bound_radius,
            bound_vertical_smeared_stress,
            "|",
            markersize=mark_size,
            color="mediumseagreen",
        )
        plt.plot(
            bound_radius,
            bound_tresca_smeared_stress,
            "|",
            markersize=mark_size,
            color="crimson",
        )
        plt.grid(True)
        plt.ylabel(r"$\sigma$ [$MPa$]", fontsize=axis_font_size)
        plt.xlabel(r"$R$ [$m$]", fontsize=axis_font_size)
        plt.legend(loc="best", fontsize=legend_size)
        plt.xticks(size=axis_tick_size)
        plt.yticks(size=axis_tick_size)
        plt.tight_layout()
        plt.savefig(f"{outdir}/smeared_stress.{save_format}")
        plt.clf()
        plt.cla()

    # PLOT 3 : Strain summary
    # ------------------------
    if plot_strain and len(sig_file_data) > 15:
        for ii in range(n_layers):
            plt.plot(
                radius[ii],
                radial_strain[ii],
                "-",
                linewidth=line_width,
                color="lightblue",
            )
            plt.plot(
                radius[ii],
                toroidal_strain[ii],
                "-",
                linewidth=line_width,
                color="wheat",
            )
            plt.plot(
                radius[ii],
                vertical_strain[ii],
                "-",
                linewidth=line_width,
                color="lightgrey",
            )
        plt.plot(
            radius[0],
            radial_strain[0],
            "--",
            color="dodgerblue",
            label=r"$\epsilon_{rr}$",
        )
        plt.plot(
            radius[0],
            toroidal_strain[0],
            "--",
            color="orange",
            label=r"$\epsilon_{\theta\theta}$",
        )
        plt.plot(
            radius[0],
            vertical_strain[0],
            "--",
            color="mediumseagreen",
            label=r"$\epsilon_{zz}$",
        )
        for ii in range(1, n_layers):
            plt.plot(radius[ii], radial_strain[ii], "--", color="dodgerblue")
            plt.plot(radius[ii], toroidal_strain[ii], "--", color="orange")
            plt.plot(radius[ii], vertical_strain[ii], "--", color="mediumseagreen")
        plt.plot(
            bound_radius,
            bound_radial_strain,
            "|",
            markersize=mark_size,
            color="dodgerblue",
        )
        plt.plot(
            bound_radius,
            bound_toroidal_strain,
            "|",
            markersize=mark_size,
            color="orange",
        )
        plt.plot(
            bound_radius,
            bound_vertical_strain,
            "|",
            markersize=mark_size,
            color="mediumseagreen",
        )
        plt.grid(True)
        plt.ylabel(r"$\epsilon$", fontsize=axis_font_size)
        plt.xlabel(r"$R$ [$m$]", fontsize=axis_font_size)
        plt.legend(loc="best", fontsize=legend_size)
        plt.xticks(size=axis_tick_size)
        plt.yticks(size=axis_tick_size)
        plt.tight_layout()
        plt.savefig(f"{outdir}/strains.{save_format}")
        plt.clf()
        plt.cla()

    # PLOT 4 : Displacement
    # ----------------------
    if plot_disp:
        plt.plot(radius[0], radial_displacement[0], color="dodgerblue")
        for ii in range(1, n_layers):
            plt.plot(radius[ii], radial_displacement[ii], color="dodgerblue")
        plt.grid(True)
        plt.ylabel(r"$u_{r}$ [mm]", fontsize=axis_font_size)
        plt.xlabel(r"$R$ [$m$]", fontsize=axis_font_size)
        plt.xticks(size=axis_tick_size)
        plt.yticks(size=axis_tick_size)
        plt.tight_layout()
        plt.savefig(f"{outdir}/displacement.{save_format}")
        plt.clf()
        plt.cla()


if __name__ == "__main__":
    main()
