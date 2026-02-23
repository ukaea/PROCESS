"""
Code generating the TF coil inboard mid-plane stress/strain summary plots
The whole radial distribution is displayed

Input file:
SIG_TF.json
"""

import json
import os
from dataclasses import dataclass
from operator import itemgetter
from pathlib import Path

import matplotlib.pyplot as plt


@dataclass
class StressPlotConfig:
    axis_font_size: float
    axis_tick_size: int = 16
    legend_size: int = 12
    mark_size: int = 13
    line_width: float = 3.5
    outdir: Path | None = None

    def __post_init__(self):
        if self.outdir is None:
            self.outdir = Path.cwd()
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)


def plot_stress(
    plot_selection,
    save_format,
    axis_font_size,
    input_file,
    term_output,
    plot_conf: StressPlotConfig | dict | None = None,
):
    if plot_conf is None:
        plot_conf = StressPlotConfig(axis_font_size)
    elif isinstance(plot_conf, dict):
        plot_conf = StressPlotConfig(axis_font_size, **plot_conf)

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

    with open(input_file) as f:
        sig_file_data = json.load(f)

    # Getting the data to be plotted
    n_radial_array_layer = sig_file_data["Points per layers"]
    n_points = len(sig_file_data["Radius (m)"])
    n_layers = int(n_points / n_radial_array_layer)

    # Assumes n_layers >= 1

    # Physical quantities : full vectors
    radius = [[] * n_layers]
    radial_smeared_stress = [[] * n_layers]
    toroidal_smeared_stress = [[] * n_layers]
    vertical_smeared_stress = [[] * n_layers]
    tresca_smeared_stress = [[] * n_layers]
    radial_stress = [[] * n_layers]
    toroidal_stress = [[] * n_layers]
    vertical_stress = [[] * n_layers]
    vm_stress = [[] * n_layers]
    tresca_stress = [[] * n_layers]
    cea_tresca_stress = [[] * n_layers]
    radial_strain = [[] * n_layers]
    toroidal_strain = [[] * n_layers]
    vertical_strain = [[] * n_layers]
    radial_displacement = [[] * n_layers]

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

    for ii in range(n_layers):
        # Full vector
        lb_ind = ii * n_radial_array_layer
        ub_ind = (ii + 1) * n_radial_array_layer - 1
        lb_ub = itemgetter(lb_ind, ub_ind)

        for jj in range(n_radial_array_layer):
            ij_ind = lb_ind + jj

            radius[ii].append(sig_file_data["Radius (m)"][ij_ind])
            radial_stress[ii].append(sig_file_data["Radial stress (MPa)"][ij_ind])
            toroidal_stress[ii].append(sig_file_data["Toroidal stress (MPa)"][ij_ind])
            vertical_stress[ii].append(
                sig_file_data["Vertical stress (MPa)"][
                    0 if len(sig_file_data["Vertical stress (MPa)"]) == 1 else ij_ind
                ]
            )
            radial_smeared_stress[ii].append(
                sig_file_data["Radial smear stress (MPa)"][ij_ind]
            )
            toroidal_smeared_stress[ii].append(
                sig_file_data["Toroidal smear stress (MPa)"][ij_ind]
            )
            vertical_smeared_stress[ii].append(
                sig_file_data["Vertical smear stress (MPa)"][ij_ind]
            )
            vm_stress[ii].append(sig_file_data["Von-Mises stress (MPa)"][ij_ind])
            tresca_stress[ii].append(sig_file_data["CEA Tresca stress (MPa)"][ij_ind])
            cea_tresca_stress[ii].append(
                sig_file_data["CEA Tresca stress (MPa)"][ij_ind]
            )
            radial_displacement[ii].append(
                sig_file_data["rad. displacement (mm)"][ij_ind]
            )

        # Layer lower/upper boundaries values
        bound_radius.extend(lb_ub(sig_file_data["Radius (m)"]))
        bound_radial_stress.extend(lb_ub(sig_file_data["Radial stress (MPa)"]))
        bound_toroidal_stress.extend(lb_ub(sig_file_data["Toroidal stress (MPa)"]))

        if len(sig_file_data["Vertical stress (MPa)"]) == 1:
            bvs_l = bvs_u = 0
        else:
            bvs_l = lb_ind
            bvs_u = ub_ind
        bound_vertical_stress.extend([
            sig_file_data["Vertical stress (MPa)"][bvs_l],
            sig_file_data["Vertical stress (MPa)"][bvs_u],
        ])

        bound_radial_smeared_stress.extend(
            lb_ub(sig_file_data["Radial smear stress (MPa)"])
        )
        bound_toroidal_smeared_stress.extend(
            lb_ub(sig_file_data["Toroidal smear stress (MPa)"])
        )
        bound_vertical_smeared_stress.extend(
            lb_ub(sig_file_data["Vertical smear stress (MPa)"])
        )
        bound_vm_stress.extend(lb_ub(sig_file_data["Von-Mises stress (MPa)"]))
        bound_tresca_stress.extend(lb_ub(sig_file_data["CEA Tresca stress (MPa)"]))
        bound_cea_tresca_stress.extend(lb_ub(sig_file_data["CEA Tresca stress (MPa)"]))
        bound_radial_displacement.extend(lb_ub(sig_file_data["rad. displacement (mm)"]))

        # Layer upper boundaries values

    # TRESCA smeared stress [MPa]
    for ii in range(n_layers):
        bound_tresca_smeared_stress.extend([
            max(abs(radial_smeared_stress[ii][0]), abs(toroidal_smeared_stress[ii][0]))
            + vertical_smeared_stress[ii][0],
            max(
                abs(radial_smeared_stress[ii][n_radial_array_layer - 1]),
                abs(toroidal_smeared_stress[ii][n_radial_array_layer - 1]),
            )
            + vertical_smeared_stress[ii][n_radial_array_layer - 1],
        ])
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
            bound_radial_strain.extend([lb_ub(sig_file_data["Radial strain"])])
            bound_toroidal_strain.extend([lb_ub(sig_file_data["Toroidal strain"])])
            bound_vertical_strain.extend([lb_ub(sig_file_data["Vertical strain"])])
            for jj in range(n_radial_array_layer):
                ij_ind = lb_ind + jj

                radial_strain[ii].append(sig_file_data["Radial strain"][ij_ind])
                toroidal_strain[ii].append(sig_file_data["Toroidal strain"][ij_ind])
                vertical_strain[ii].append(sig_file_data["Vertical strain"][ij_ind])

                if "WP smeared stress (MPa)" in sig_file_data:
                    wp_vertical_stress.append(
                        sig_file_data["WP smeared stress (MPa)"][jj]
                    )

    if term_output:
        terminal_output(
            n_layers,
            n_radial_array_layer,
            sig_file_data,
            radial_stress,
            toroidal_stress,
            vertical_stress,
            tresca_stress,
            wp_vertical_stress,
            radial_smeared_stress,
            toroidal_smeared_stress,
            vertical_smeared_stress,
            tresca_smeared_stress,
            radial_strain,
            toroidal_strain,
            vertical_strain,
        )

    if plot_sig:
        stress_summary(
            n_layers,
            radius,
            bound_radius,
            radial_stress,
            toroidal_stress,
            vertical_stress,
            tresca_stress,
            vm_stress,
            bound_radial_stress,
            bound_toroidal_stress,
            bound_vertical_stress,
            bound_tresca_stress,
            bound_vm_stress,
            save_format,
            plot_conf,
        )

    if plot_sm_sig:
        smeared_stress_summary(
            n_layers,
            radius,
            bound_radius,
            radial_smeared_stress,
            toroidal_smeared_stress,
            vertical_smeared_stress,
            tresca_smeared_stress,
            bound_radial_smeared_stress,
            bound_toroidal_smeared_stress,
            bound_vertical_smeared_stress,
            bound_tresca_smeared_stress,
            save_format,
            plot_conf,
        )

    if plot_strain and len(sig_file_data) > 15:
        strain_summary(
            n_layers,
            radius,
            bound_radius,
            radial_strain,
            bound_radial_strain,
            toroidal_strain,
            bound_toroidal_strain,
            vertical_strain,
            bound_vertical_strain,
            save_format,
            plot_conf,
        )

    if plot_disp:
        displacement(n_layers, radius, radial_displacement, save_format, plot_conf)


def terminal_output(
    n_layers,
    n_radial_array_layer,
    sig_file_data,
    radial_stress,
    toroidal_stress,
    vertical_stress,
    tresca_stress,
    wp_vertical_stress,
    radial_smeared_stress,
    toroidal_smeared_stress,
    vertical_smeared_stress,
    tresca_smeared_stress,
    radial_strain,
    toroidal_strain,
    vertical_strain,
):
    ii_ins = 0
    ii_mids = int(0.5 * float(n_radial_array_layer))
    ii_outs = n_radial_array_layer - 1
    dg = itemgetter(ii_ins, ii_mids, ii_outs)

    print("\n\nLayer stress details\n____________________")

    frame = """Layer {}
------------------------------
steel radial   stress in the inner/middle/out point: {} MPa
steel toroidal stress in the inner/middle/out point: {} MPa
steel vertical stress in the inner/middle/out point: {} MPa
steel TRESCA   stress in the inner/middle/out point: {} MPa
smeared radial   stress in the inner/middle/out point: {} MPa
smeared toroidal stress in the inner/middle/out point: {} MPa
smeared vertical stress in the inner/middle/out point: {} MPa
smeared TRESCA   stress in the inner/middle/out point: {} MPa

"""
    frame2 = """
radial   strain in the inner/middle/out point: {}
toroidal strain in the inner/middle/out point: {}
vertical strain: {}

"""
    frame3 = "smeared WP vertical stress in the inner/middle/out point: {} MPa"
    layer_line = "{}/{}/{}"
    for ii in range(n_layers):
        layer = ii + 1
        s_radial = layer_line.format(*dg(radial_stress[ii]))
        s_toro = layer_line.format(*dg(toroidal_stress[ii]))
        s_vert = layer_line.format(*dg(vertical_stress[ii]))
        s_tres = layer_line.format(*dg(tresca_stress[ii]))
        sm_rad = layer_line.format(*dg(radial_smeared_stress[ii]))
        sm_toro = layer_line.format(*dg(toroidal_smeared_stress[ii]))
        sm_vert = layer_line.format(*dg(vertical_smeared_stress[ii]))
        sm_tres = layer_line.format(*dg(tresca_smeared_stress[ii]))

        print(
            frame.format(
                layer,
                s_radial,
                s_toro,
                s_vert,
                s_tres,
                sm_rad,
                sm_toro,
                sm_vert,
                sm_tres,
            )
        )

        if len(sig_file_data) > 16:
            r_strain = layer_line.format(*dg(radial_strain[ii]))
            t_strain = layer_line.format(*dg(toroidal_strain[ii]))
            print(frame2.format(r_strain, t_strain, vertical_strain[ii][0]))

    if len(wp_vertical_stress) != 0:
        print(
            frame3.format(
                wp_vertical_stress[0],
                wp_vertical_stress[ii_mids],
                wp_vertical_stress[ii_outs],
            )
        )
    print()


def stress_summary(
    n_layers,
    radius,
    bound_radius,
    radial_stress,
    toroidal_stress,
    vertical_stress,
    tresca_stress,
    vm_stress,
    bound_radial_stress,
    bound_toroidal_stress,
    bound_vertical_stress,
    bound_tresca_stress,
    bound_vm_stress,
    save_format,
    plot_conf,
):
    lw = plot_conf.line_width
    ms = plot_conf.mark_size
    ats = plot_conf.axis_tick_size
    afs = plot_conf.axis_font_size
    for ii in range(n_layers):
        plt.plot(radius[ii], radial_stress[ii], "-", linewidth=lw, color="lightblue")
        plt.plot(radius[ii], toroidal_stress[ii], "-", linewidth=lw, color="wheat")
        plt.plot(radius[ii], vertical_stress[ii], "-", linewidth=lw, color="lightgrey")
        plt.plot(radius[ii], tresca_stress[ii], "-", linewidth=lw, color="pink")
        plt.plot(radius[ii], vm_stress[ii], "-", linewidth=lw, color="violet")
    plt.plot(
        radius[0], radial_stress[0], "--", color="dodgerblue", label=r"$\sigma_{rr}$"
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
        radius[0], tresca_stress[0], "-", color="crimson", label=r"$\sigma_{TRESCA}$"
    )
    plt.plot(
        radius[0], vm_stress[0], "-", color="darkviolet", label=r"$\sigma_{Von\ mises}$"
    )
    for ii in range(1, n_layers):
        plt.plot(radius[ii], radial_stress[ii], "--", color="dodgerblue")
        plt.plot(radius[ii], toroidal_stress[ii], "--", color="orange")
        plt.plot(radius[ii], vertical_stress[ii], "--", color="mediumseagreen")
        plt.plot(radius[ii], tresca_stress[ii], "-", color="crimson")
        plt.plot(radius[ii], vm_stress[ii], "-", color="darkviolet")
    plt.plot(bound_radius, bound_radial_stress, "|", markersize=ms, color="dodgerblue")
    plt.plot(bound_radius, bound_toroidal_stress, "|", markersize=ms, color="orange")
    plt.plot(
        bound_radius, bound_vertical_stress, "|", markersize=ms, color="mediumseagreen"
    )
    plt.plot(bound_radius, bound_tresca_stress, "|", markersize=ms, color="crimson")
    plt.plot(bound_radius, bound_vm_stress, "|", markersize=ms, color="darkviolet")
    plt.grid(True)
    plt.ylabel(r"$\sigma$ [$MPa$]", fontsize=afs)
    plt.xlabel(r"$R$ [$m$]", fontsize=afs)
    plt.legend(loc="best", fontsize=plot_conf.legend_size)
    plt.xticks(size=ats)
    plt.yticks(size=ats)
    plt.tight_layout()
    plt.savefig(f"{plot_conf.outdir}/structure_stress.{save_format}")
    plt.clf()
    plt.cla()


def smeared_stress_summary(
    n_layers,
    radius,
    bound_radius,
    radial_smeared_stress,
    toroidal_smeared_stress,
    vertical_smeared_stress,
    tresca_smeared_stress,
    bound_radial_smeared_stress,
    bound_toroidal_smeared_stress,
    bound_vertical_smeared_stress,
    bound_tresca_smeared_stress,
    save_format,
    plot_conf,
):
    lw = plot_conf.line_width
    ms = plot_conf.mark_size
    ats = plot_conf.axis_tick_size
    afs = plot_conf.axis_font_size
    for ii in range(n_layers):
        plt.plot(
            radius[ii], radial_smeared_stress[ii], "-", linewidth=lw, color="lightblue"
        )
        plt.plot(
            radius[ii], toroidal_smeared_stress[ii], "-", linewidth=lw, color="wheat"
        )
        plt.plot(
            radius[ii], vertical_smeared_stress[ii], "-", linewidth=lw, color="lightgrey"
        )
        plt.plot(radius[ii], tresca_smeared_stress[ii], "-", linewidth=lw, color="pink")
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
        plt.plot(radius[ii], vertical_smeared_stress[ii], "--", color="mediumseagreen")
        plt.plot(radius[ii], tresca_smeared_stress[ii], "-", color="crimson")
    plt.plot(
        bound_radius, bound_radial_smeared_stress, "|", markersize=ms, color="dodgerblue"
    )
    plt.plot(
        bound_radius, bound_toroidal_smeared_stress, "|", markersize=ms, color="orange"
    )
    plt.plot(
        bound_radius,
        bound_vertical_smeared_stress,
        "|",
        markersize=ms,
        color="mediumseagreen",
    )
    plt.plot(
        bound_radius, bound_tresca_smeared_stress, "|", markersize=ms, color="crimson"
    )
    plt.grid(True)
    plt.ylabel(r"$\sigma$ [$MPa$]", fontsize=afs)
    plt.xlabel(r"$R$ [$m$]", fontsize=afs)
    plt.legend(loc="best", fontsize=plot_conf.legend_size)
    plt.xticks(size=ats)
    plt.yticks(size=ats)
    plt.tight_layout()
    plt.savefig(f"{plot_conf.outdir}/smeared_stress.{save_format}")
    plt.clf()
    plt.cla()


def strain_summary(
    n_layers,
    radius,
    bound_radius,
    radial_strain,
    bound_radial_strain,
    toroidal_strain,
    bound_toroidal_strain,
    vertical_strain,
    bound_vertical_strain,
    save_format,
    plot_conf,
):
    lw = plot_conf.line_width
    ms = plot_conf.mark_size
    ats = plot_conf.axis_tick_size
    afs = plot_conf.axis_font_size
    for ii in range(n_layers):
        plt.plot(radius[ii], radial_strain[ii], "-", linewidth=lw, color="lightblue")
        plt.plot(radius[ii], toroidal_strain[ii], "-", linewidth=lw, color="wheat")
        plt.plot(radius[ii], vertical_strain[ii], "-", linewidth=lw, color="lightgrey")
    plt.plot(
        radius[0], radial_strain[0], "--", color="dodgerblue", label=r"$\epsilon_{rr}$"
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
    plt.plot(bound_radius, bound_radial_strain, "|", markersize=ms, color="dodgerblue")
    plt.plot(bound_radius, bound_toroidal_strain, "|", markersize=ms, color="orange")
    plt.plot(
        bound_radius,
        bound_vertical_strain,
        "|",
        markersize=ms,
        color="mediumseagreen",
    )
    plt.grid(True)
    plt.ylabel(r"$\epsilon$", fontsize=afs)
    plt.xlabel(r"$R$ [$m$]", fontsize=afs)
    plt.legend(loc="best", fontsize=plot_conf.legend_size)
    plt.xticks(size=ats)
    plt.yticks(size=ats)
    plt.tight_layout()
    plt.savefig(f"{plot_conf.outdir}/strains.{save_format}")
    plt.clf()
    plt.cla()


def displacement(n_layers, radius, radial_displacement, save_format, plot_conf):
    plt.plot(radius[0], radial_displacement[0], color="dodgerblue")
    for ii in range(1, n_layers):
        plt.plot(radius[ii], radial_displacement[ii], color="dodgerblue")
    plt.grid(True)
    plt.ylabel(r"$u_{r}$ [mm]", fontsize=plot_conf.axis_font_size)
    plt.xlabel(r"$R$ [$m$]", fontsize=plot_conf.axis_font_size)
    plt.xticks(size=plot_conf.axis_tick_size)
    plt.yticks(size=plot_conf.axis_tick_size)
    plt.tight_layout()
    plt.savefig(f"{plot_conf.outdir}/displacement.{save_format}")
    plt.clf()
    plt.cla()
