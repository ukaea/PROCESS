"""
Code to display bar type representation of the radial build of a PROCESS run.

Author: A. Pearce (alexander.pearce@ukaea.uk)
Updated 26/09/23: C. Ashe (christopher.ashe@ukaea.uk)

Input file:
MFILE.DAT

"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path

# PROCESS libraries
import process.io.mfile as mf


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
        "-op",
        "--outputdir",
        default=Path.cwd(),
        help="Output directory for plot, defaults to current working directory.",
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
    else:
        pass

    radial_labels = [
        "bore",
        "ohcth",
        "precomp",
        "gapoh",
        "tfcth",
        "tftsgap",
        "thshield_ib",
        "gapds",
        "d_vv_in",
        "shldith",
        "vvblgap",
        "blnkith",
        "fwith",
        "scrapli",
        "rminor",
        "scraplo",
        "fwoth",
        "blnkoth",
        "vvblgap",
        "d_vv_out",
        "shldoth",
        "gapsto",
        "thshield_ob",
        "tftsgap",
        "tfthko",
    ]
    if int(m_file.data["tf_in_cs"].get_scan(-1)) == 1:
        radial_labels[1] = "tfcth"
        radial_labels[2] = "gapoh"
        radial_labels[3] = "ohcth"
        radial_labels[4] = "precomp"
        radial_labels[5] = "tftsgap"

    radial_build = []

    for ii in range(isweep):
        if m_file.data["ifail"].get_scan(ii + 1) == 1:
            radial_build.append(
                [m_file.data[rl].get_scan(ii + 1) for rl in radial_labels]
            )

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

    # LaTeX labels
    # ------------
    # ToDo : WOULD BE GREAT TO HAVE IT STORED IN THE PROCESS !
    labels = dict()
    labels["shldith"] = r"$\Delta R_\mathrm{sh}$ [$m$]"
    labels["rmajor"] = r"$R_\mathrm{maj}$ [$m$]"
    labels["crypmw"] = r"$P_\mathrm{cryo}$ [$MW$]"
    labels["bt"] = r"$B_\mathrm{T}$ [$T$]"
    labels["tfcth"] = r"$\Delta R_\mathrm{TF}$ [$m$]"
    labels["powfmw"] = r"$P_\mathrm{fus}$ [$MW$]"
    labels["pinjemw"] = r"$P_\mathrm{inj}$ [$MW$]"
    labels["pnetelmw"] = r"$P_\mathrm{Net\ elec}$ [$MW$]"
    labels["taueff"] = r"$\tau_\mathrm{E}$ [s]"
    labels["ralpne"] = r"$f_\mathrm{\alpha}$"
    labels["te"] = r"$\left< T_\mathrm{e} \right>$"
    labels["taulimit"] = r"$max : \frac{\tau_\mathrm{\alpha}}{\tau_\mathrm{E}}$"
    labels["scrapli"] = r"$\Delta R_\mathrm{FW-sep}$ [$m$]"
    labels["scraplo"] = r"$\Delta R_\mathrm{FW-sep}^\mathrm{out}$ [$m$]"
    labels["vforce"] = r"$F_\mathrm{z}^\mathrm{in}$ [$N$]"
    labels["thkcas"] = r"$\Delta R_\mathrm{TF}^\mathrm{buck}$ [$m$]"
    labels["bmaxtf"] = r"$B_\mathrm{TF}^\mathrm{max}$ [$T$]"
    labels["ritfc"] = r"$I_\mathrm{TF}^\mathrm{tot}$ [$A$]"
    labels["dr_tf_wp"] = r"$\Delta R_\mathrm{TF}^\mathrm{WP}$ [$m$]"
    labels["aspect"] = r"$A$"
    labels["rminor"] = r"$a_\mathrm{min}$ [$m$]"
    labels["capcost"] = r"$C_\mathrm{cap}$ [$M\$ $]"
    labels["r_tf_outboard_mid"] = r"$\Delta R_\mathrm{TF}^\mathrm{out\ mid}$ [$m$]"
    labels["pgrossmw"] = r"$P_\mathrm{gross}^\mathrm{elec}$ [$MW$]"
    labels["htpmw"] = r"$P_\mathrm{Primary\ coolant}^\mathrm{elec}$ [$MW$]"
    labels["ppfmw"] = r"$P_\mathrm{PF}^\mathrm{elec}$ [$MW$]"
    labels["hmax"] = r"$z_\mathrm{TF}^\mathrm{pl\ side}$ [$m$]"
    labels["thicndut"] = r"\Delta l_\mathrm{steel\ jacket}^\mathrm{turn}"
    labels["cpttf"] = r"$I_\mathrm{TF}^\mathrm{turn}$ [$A$]"
    labels["boundl(2)"] = r"$B_\mathrm{T}^\mathrm{min}$ [$A$]"
    labels["pinjmw"] = r"$P_\mathrm{inj}$ [$MW$]"
    labels["hldivlim"] = r"$q_\mathrm{div}^\mathrm{max}$ [$MW.m^{-2}$]"
    labels["hfact"] = r"$f_\mathrm{H}$"
    labels["kappa"] = r"$\kappa_\mathrm{sep}$"
    labels["triang"] = r"$\delta_\mathrm{sep}$"
    labels["f_tf_steel"] = r"f_\mathrm{steel}^\mathrm{TF}"
    labels["plascur/1d6"] = r"$I_{\mathrm{p}}$[$MA$]"
    labels["n_cycle"] = r"$N_{\mathrm{CS},\mathrm{cycle}}$"
    labels["alstroh"] = r"$\sigma_{\mathrm{oh}}^{\mathrm{max}}$[$Pa$]"
    labels["ohcth"] = r"$\Delta R_{\mathrm{CS}}$[$m$]"
    labels["bore"] = r"$\Delta R_{\mathrm{bore}}$[$m$]"
    labels["dnla"] = r"$\bar{n}_{\mathrm{e}}$[$m^{-3}$]"
    labels["dnla_gw"] = r"$f_{\mathrm{GW}}$"
    labels["normalised_toroidal_beta"] = r"$\beta_{N,\mathrm{tor}}$"
    labels["copperaoh_m2"] = r"$\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$"
    labels["copperaoh_m2_max"] = r"$max\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$"
    labels["coreradius"] = r"$r_{core} [M]$"
    labels[
        "fcuohsu"
    ] = r"$f_{\mathrm{Cu}}^{\mathrm{CS}}$"  # copper fraction of strand in central solenoid
    labels["coheof"] = r"$J [A M^{-2}]$"
    labels["ohhghf"] = r"$Thickness_{\mathrm{CS}}[m]$"
    labels["pheat"] = r"$ P_{\mathrm{heat}}$ [$MW$]"
    labels["effcd"] = r"$\eta_{\mathrm{CD}}$[$A/W$]"
    labels["bigq"] = r"$Q$"
    labels["faccd"] = r"$f_{\mathrm{CD}}$"
    labels["facoh"] = r"$f_{\mathrm{CD,ind}}$"
    labels["bootipf"] = r"$f_{\mathrm{BS}}$"
    labels["itvar039"] = r"$f_{\mathrm{LH}}$"
    labels["itvar011"] = r"$f_{\mathrm{V.s}}$"
    labels["pdivt"] = r"$P_{\mathrm{sep}}$ [$MW$]"
    labels["pradmw"] = r"$P_{\mathrm{rad}}$ [$MW$]"
    labels[
        "pdivtbt/qar"
    ] = r"$\frac{P_{\mathrm{sep}}B_T}{q_{95}AR_{\mathrm{maj}}}$ [$MWTm^{-1}$]"
    labels["iooic"] = r"$I_{\mathrm{TF}}/I_{\mathrm{TF},\mathrm{crit}}$"
    labels["bktlife"] = r"$T_{\mathrm{blk}}$"
    labels["bktcycles"] = r"$N_{\mathrm{blk},\mathrm{cycle}}$"
    labels["zeff"] = r"$Z_{\mathrm{eff}}$"
    labels["tburn"] = r"$t_{\mathrm{burn}}$[$s$]"
    labels["vburn"] = r"$V_{\mathrm{loop}}$ [$V$]"
    labels["rli"] = r"$l_i$"
    labels["csfv.n_cycle_min"] = r"$N_{\mathrm{min}}^{\mathrm{CS}}$"
    labels["csfv.n_cycle"] = r"$N_{\mathrm{cycle}}^{\mathrm{CS}}$"
    labels["a_oh_turn"] = r"$Turn_{\mathrm{area}}^{\mathrm{CS}}[$m$^{2}]$"
    labels["tbrnmn"] = r"$t_{\mathrm{burn.min}}$[$s$]"
    labels["pfv.oh_steel_frac"] = r"$f_{\mathrm{Steel}}^{\mathrm{CS}}$"
    labels["csfv.t_structural_radial"] = r"$Turn_{\mathrm{radial}}^{\mathrm{CS}}[$m$]$"
    labels["csfv.t_crack_vertical"] = r"$Crack_{\mathrm{vertical}}^{\mathrm{CS}}[$m$]$"
    labels["sig_hoop"] = r"$\sigma_{\mathrm{oh},{\mathrm{hoop}}}$[$Pa$]"
    labels["rplas"] = r"$R_{\mathrm{plas}}$[$\Omega$]"
    # ------------

    # nsweep varible dict
    # -------------------
    # TODO WOULD BE GREAT TO HAVE IT AUTOMATICALLY GENERATED ON THE PROCESS CMAKE!
    #        THE SAME WAY THE DICTS ARE
    # This needs to be kept in sync automatically; this will break frequently
    # otherwise
    # Rem : Some variables are not in the MFILE, making the defintion rather tricky...

    nsweep_list = [
        "aspect",
        "hldivlim",
        "pnetelmw",
        "hfact",
        "oacdcp",
        "walalw",
        "beamfus0",
        "fqval",
        "te",
        "boundu(15)",
        "dnbeta",
        "bscfmax",
        "boundu(10)",
        "fiooic",
        "fjprot",
        "rmajor",
        "bmaxtf",  # bmxlim the maximum T field upper limit is the scan variable
        "gammax",
        "boundl(16)",
        "cnstv.tbrnmn",
        "",
        "cfactr",
        "boundu(72)",
        "powfmax",
        "kappa",
        "triang",
        "tbrmin",
        "bt",
        "coreradius",
        "Obsolete",  # Removed
        "taulimit",
        "epsvmc",
        "ttarget",
        "qtargettotal",
        "lambda_q_omp",
        "lambda_target",
        "lcon_factor",
        "boundu(129)",
        "boundu(131)",
        "boundu(135)",
        "blnkoth",
        "fimp(9)",
        "Obsolete",  # Removed
        "alstrtf",
        "tmargmin_tf",
        "boundu(152)",
        "impurity_enrichment(9)",
        "n_pancake",
        "n_layer",
        "fimp(13)",
        "ftar",
        "rad_fraction_sol",
        "",
        "b_crit_upper_nbti",
        "shldith",
        "crypmw_max",
        "bt",  # Genuinly bt lower bound
        "scrapli",
        "scraplo",
        "sig_tf_wp_max",
        "copperaoh_m2_max",
        "coheof",
        "ohcth",
        "ohhghf",
        "csfv.n_cycle_min",
        "pfv.oh_steel_frac",
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
        "fvs",  # actaully lower bound fvs
        "vburn",
        "rplas",
    ]

    # "plasma_res_factor"
    # -------------------

    # Getting the scanned variable name
    m_file = mf.MFile(filename=input_file)
    nsweep_ref = int(m_file.data["nsweep"].get_scan(-1))
    if nsweep_ref == 0:
        scan_var_name = "Null"
    else:
        scan_var_name = nsweep_list[nsweep_ref - 1]

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
    if int(m_file.data["tf_in_cs"].get_scan(-1)) == 1:
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
    if int(m_file.data["tf_in_cs"].get_scan(-1)) == 1:
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
    else:
        pass
    if args.inboard:
        end_scan = radial_labels.index("Plasma")
    else:
        end_scan = len(radial_build)
    plt.figure(figsize=(8, 6))
    for kk in range((len(radial_build[:end_scan, 0]))):
        if kk == 0:
            lower = np.zeros(len(radial_build[kk, :]))
        else:
            lower = lower + radial_build[kk - 1, :]
        plt.barh(
            ind if scan_var_name != "Null" else 0,
            radial_build[kk, :],
            left=lower,
            height=0.8,
            label=f"{radial_labels[kk]}" + f"\n {radial_build[kk][0]} m" * args.numbers,
            color=radial_color[kk],
            edgecolor="black",
            linewidth=0.05,
        )

    if scan_var_name != "Null":
        plt.yticks(ind, scan_points, fontsize=axis_tick_size)
        plt.ylabel(labels[scan_var_name], fontsize=axis_font_size)
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
    plt.savefig(
        f"{args.outputdir}/{Path(args.input).stem}_radial_build.{save_format}",
        bbox_inches="tight",
    )

    # Display plot (used in Jupyter notebooks)
    plt.show()
    plt.clf()


if __name__ == "__main__":
    main()
