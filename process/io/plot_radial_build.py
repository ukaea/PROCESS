"""
Code to display bar type representation of the radial build of a PROCESS run.

Author: A. Pearce (alexander.pearce@ukaea.uk)
Updated 26/09/23: C. Ashe (christopher.ashe@ukaea.uk)

Input file:
MFILE.DAT

"""

import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# PROCESS libraries
import process.io.mfile as mf
from process.io.variable_metadata import var_dicts as meta


def parse_args(args):
    """Parse supplied arguments.

    :param args: arguments to parse
    :type args: list, None
    :return: parsed arguments
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(
        description="Plot optimization information",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "-f",
        "--input",
        default="MFILE.DAT",
        help=("Specify input file path (default = MFILE.DAT)"),
    )

    parser.add_argument(
        "-sf",
        "--save_format",
        nargs="?",
        default="pdf",
        help="Output format (default='pdf') ",
    )

    parser.add_argument(
        "-ib",
        "--inboard",
        action="store_true",
        default=False,
        help="Show inboard build only, default = False ",
    )

    parser.add_argument(
        "-nm",
        "--numbers",
        action="store_true",
        default=False,
        help="Show widths of components in the legend. Only use non-scan MFILE's as will only show last values",
    )

    return parser.parse_args(args)


def get_radial_build(m_file):
    isweep = int(m_file.data["isweep"].get_scan(-1))
    if isweep == 0:
        isweep = 1

    radial_labels = [
        "dr_bore",
        "dr_cs",
        "dr_cs_precomp",
        "dr_cs_tf_gap",
        "dr_tf_inboard",
        "dr_tf_shld_gap",
        "dr_shld_thermal_inboard",
        "dr_shld_vv_gap_inboard",
        "dr_vv_inboard",
        "dr_shld_inboard",
        "dr_shld_blkt_gap",
        "dr_blkt_inboard",
        "dr_fw_inboard",
        "dr_fw_plasma_gap_inboard",
        "rminor",
        "dr_fw_plasma_gap_outboard",
        "dr_fw_outboard",
        "dr_blkt_outboard",
        "dr_shld_blkt_gap",
        "dr_vv_outboard",
        "dr_shld_outboard",
        "dr_shld_vv_gap_outboard",
        "dr_shld_thermal_outboard",
        "dr_tf_shld_gap",
        "dr_tf_outboard",
    ]
    if int(m_file.data["i_tf_inside_cs"].get_scan(-1)) == 1:
        radial_labels[1] = "dr_tf_inboard"
        radial_labels[2] = "dr_cs_tf_gap"
        radial_labels[3] = "dr_cs"
        radial_labels[4] = "dr_cs_precomp"
        radial_labels[5] = "dr_tf_shld_gap"

    radial_build = [
        [m_file.data[rl].get_scan(ii + 1) for rl in radial_labels]
        for ii in range(isweep)
        if m_file.data["ifail"].get_scan(ii + 1) == 1
    ]

    radial_build = np.array(radial_build)

    # plasma is 2*rminor
    # Therefore we must count it again

    for kk in range(radial_build.shape[0]):
        radial_build[kk, 14] = 2.0 * radial_build[kk, 14]

    return radial_build.T, radial_build.shape[0]


def main(args=None):
    args = parse_args(args)

    input_file = str(args.input)
    save_format = str(args.save_format)

    # nsweep varible dict
    # -------------------
    # TODO WOULD BE GREAT TO HAVE IT AUTOMATICALLY GENERATED ON THE PROCESS CMAKE!
    #        THE SAME WAY THE DICTS ARE
    # This needs to be kept in sync automatically; this will break frequently
    # otherwise
    # Rem : Some variables are not in the MFILE, making the defintion rather tricky...

    nsweep_list = [
        "aspect",
        "pflux_div_heat_load_max_mw",
        "p_plant_electric_net_mw",
        "hfact",
        "oacdcp",
        "pflux_fw_neutron_max_mw",
        "beamfus0",
        "temp_plasma_electron_vol_avg_kev",
        "boundu(15)",
        "beta_norm_max",
        "f_c_plasma_bootstrap_max",
        "boundu(10)",
        "fiooic",
        "rmajor",
        "b_tf_inboard_peak_symmetric",  # b_tf_inboard_max the maximum T field upper limit is the scan variable
        "eta_cd_norm_hcd_primary_max",
        "boundl(16)",
        "cnstv.t_burn_min",
        "",
        "f_t_plant_available",
        "boundu(72)",
        "p_fusion_total_max_mw",
        "kappa",
        "triang",
        "tbrmin",
        "b_plasma_toroidal_on_axis",
        "radius_plasma_core_norm",
        "Obsolete",  # Removed
        "f_alpha_energy_confinement_min",
        "epsvmc",
        "ttarget",
        "qtargettotal",
        "lambda_q_omp",
        "lambda_target",
        "lcon_factor",
        "boundu(129)",
        "boundu(131)",
        "boundu(135)",
        "dr_blkt_outboard",
        "f_nd_impurity_electrons(9)",
        "Obsolete",  # Removed
        "alstrtf",
        "temp_tf_superconductor_margin_min",
        "boundu(152)",
        "impurity_enrichment(9)",
        "n_tf_wp_pancakes",
        "n_tf_wp_layers",
        "f_nd_impurity_electrons(13)",
        "f_p_div_lower",
        "rad_fraction_sol",
        "",
        "b_crit_upper_nbti",
        "dr_shld_inboard",
        "p_cryo_plant_electric_max_mw",
        "b_plasma_toroidal_on_axis",  # Genuinly b_plasma_toroidal_on_axis lower bound
        "dr_fw_plasma_gap_inboard",
        "dr_fw_plasma_gap_outboard",
        "sig_tf_wp_max",
        "copperaoh_m2_max",
        "j_cs_flat_top_end",
        "dr_cs",
        "f_z_cs_tf_internal",
        "csfv.n_cycle_min",
        "pfv.f_a_cs_turn_steel",
        "csfv.t_crack_vertical",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "v_plasma_loop_burn",
        "res_plasma",
    ]

    # "plasma_res_factor"
    # -------------------

    # Getting the scanned variable name
    m_file = mf.MFile(filename=input_file)
    nsweep_ref = int(m_file.data["nsweep"].get_scan(-1))
    scan_var_name = "Null" if nsweep_ref == 0 else nsweep_list[nsweep_ref - 1]

    radial_labels = [
        "Machine Bore",
        "Central Solenoid",
        "CS precompression",
        "CS Coil gap",
        "TF Coil Inboard Leg",
        "TF Coil gap",
        "Inboard Thermal Shield",
        "Gap",
        "Inboard VV",
        "Inboard Shield",
        "Gap",
        "Inboard Blanket",
        "Inboard First Wall",
        "Inboard SOL",
        "Plasma",
        "Outboard SOL",
        "Outboard First Wall",
        "Outboard Blanket",
        "Gap",
        "Outboard VV",
        "Outboard Shield",
        "Gap",
        "Outboard Thermal Shield",
        "Gap",
        "TF Coil Outboard Leg",
    ]
    if int(m_file.data["i_tf_inside_cs"].get_scan(-1)) == 1:
        radial_labels[1] = "TF Coil Inboard Leg"
        radial_labels[2] = "CS Coil gap"
        radial_labels[3] = "Central Solenoid"
        radial_labels[4] = "CS precompression"
        radial_labels[5] = "TF Coil gap"
    radial_color = [
        "lightgrey",
        "green",
        "yellow",
        "white",
        "blue",
        "white",
        "lime",
        "white",
        "dimgrey",
        "violet",
        "white",
        "goldenrod",
        "steelblue",
        "orange",
        "red",
        "orange",
        "steelblue",
        "goldenrod",
        "white",
        "dimgrey",
        "violet",
        "white",
        "lime",
        "white",
        "blue",
    ]
    if int(m_file.data["i_tf_inside_cs"].get_scan(-1)) == 1:
        radial_color[1] = "blue"
        radial_color[2] = "white"
        radial_color[3] = "green"
        radial_color[4] = "yellow"
        radial_color[5] = "white"
    radial_build, num_converged_sol = get_radial_build(m_file)

    # Get scan variable data
    if scan_var_name != "Null":
        nn = 0
        isweep = int(m_file.data["isweep"].get_scan(-1))
        scan_points = np.zeros(num_converged_sol)
        for ii in range(isweep):
            ifail = m_file.data["ifail"].get_scan(ii + 1)
            if ifail == 1:
                scan_points[nn] = m_file.data[scan_var_name].get_scan(ii + 1)
                nn += 1
    else:
        scan_points = 1
    index = []
    # need a set of checks - remove build parts equal to zero
    for ll in range(len(radial_build[:, 0])):
        if radial_build[ll, 0] == 0.0:
            # note index
            index = np.append(index, ll)

    # Plot settings
    # -------------
    # Plot cosmetic settings
    axis_tick_size = 12
    legend_size = 8
    axis_font_size = 16
    if scan_var_name != "Null":
        ind = [y for y, _ in enumerate(scan_points)]
    end_scan = radial_labels.index("Plasma") if args.inboard else len(radial_build)
    plt.figure(figsize=(8, 6))
    for kk in range(len(radial_build[:end_scan, 0])):
        if kk == 0:
            lower = np.zeros(len(radial_build[kk, :]))
        else:
            lower = lower + radial_build[kk - 1, :]
        plt.barh(
            ind if scan_var_name != "Null" else 0,
            radial_build[kk, :],
            left=lower,
            height=0.8,
            label=f"{radial_labels[kk]}"
            + f"\n {radial_build[kk][0]:.3f} m" * args.numbers,
            color=radial_color[kk],
            edgecolor="black",
            linewidth=0.05,
        )

    if scan_var_name != "Null":
        plt.yticks(ind, scan_points, fontsize=axis_tick_size)
        plt.ylabel(
            meta[scan_var_name].latex if scan_var_name in meta else f"{scan_var_name}",
            fontsize=axis_font_size,
        )
    else:
        plt.yticks([])

    plt.legend(
        bbox_to_anchor=(0.5, -0.15),
        loc="upper center",
        fontsize=legend_size,
        ncol=4,
    )
    plt.xlabel("Radius [m]")
    plt.tight_layout()
    output_path = (
        Path(args.input).parent / f"{Path(args.input).stem}_radial_build.{save_format}"
    )
    plt.savefig(output_path, bbox_inches="tight")

    # Display plot (used in Jupyter notebooks)
    plt.show()
    plt.clf()


if __name__ == "__main__":
    main()
