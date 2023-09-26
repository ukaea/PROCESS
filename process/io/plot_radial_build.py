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
        "-o",
        "--output",
        default="radial_build_plot",
        help="name out outputted file (default=radial_build_plot)",
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
        help="Show inboard build only, default = No ",
    )

    return parser.parse_args(args)


def get_radial_build(m_file):
    isweep = int(m_file.data["isweep"].get_scan(-1))

    # Find out the number of converged solutions
    ifail = []
    for kk in range(isweep):
        ifail = np.append(ifail, m_file.data["ifail"].get_scan(kk + 1))
    num_converged_sol = np.count_nonzero(ifail == 1)
    radial_build = np.zeros((25, num_converged_sol))

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

    ll = 0
    for ii in range(isweep):
        if ifail[ii] == 1:
            for jj in range(len(radial_labels)):
                radial_build[jj, ll] = m_file.data[radial_labels[jj]].get_scan(ii + 1)
            ll += 1

    # plasma is 2*rminor
    # Therefore we must count it again
    for kk in range(num_converged_sol):
        radial_build[14, kk] = 2.0 * radial_build[14, kk]

    return radial_build, num_converged_sol


def main(args=None):
    args = parse_args(args)

    input_file = str(args.input)
    output_name = str(args.output)
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
    nsweep_dict = dict()
    nsweep_dict[1] = "aspect"
    nsweep_dict[2] = "hldivlim"
    nsweep_dict[3] = "pnetelmw"
    nsweep_dict[4] = "hfact"
    nsweep_dict[5] = "oacdcp"
    nsweep_dict[6] = "walalw"
    nsweep_dict[7] = "beamfus0"
    nsweep_dict[8] = "fqval"
    nsweep_dict[9] = "te"
    nsweep_dict[10] = "boundu(15)"
    nsweep_dict[11] = "dnbeta"
    nsweep_dict[12] = "bscfmax"
    nsweep_dict[13] = "boundu(10)"
    nsweep_dict[14] = "fiooic"
    nsweep_dict[15] = "fjprot"
    nsweep_dict[16] = "rmajor"
    nsweep_dict[
        17
    ] = "bmaxtf"  # bmxlim the maximum T field upper limit is the scan variable
    nsweep_dict[18] = "gammax"
    nsweep_dict[19] = "boundl(16)"
    nsweep_dict[20] = "cnstv.tbrnmn"
    nsweep_dict[21] = ""
    nsweep_dict[22] = "cfactr"
    nsweep_dict[23] = "boundu(72)"
    nsweep_dict[24] = "powfmax"
    nsweep_dict[25] = "kappa"
    nsweep_dict[26] = "triang"
    nsweep_dict[27] = "tbrmin"
    nsweep_dict[28] = "bt"
    nsweep_dict[29] = "coreradius"
    nsweep_dict[30] = "fimpvar"
    nsweep_dict[31] = "taulimit"
    nsweep_dict[32] = "epsvmc"
    nsweep_dict[33] = "ttarget"
    nsweep_dict[34] = "qtargettotal"
    nsweep_dict[35] = "lambda_q_omp"
    nsweep_dict[36] = "lambda_target"
    nsweep_dict[37] = "lcon_factor"
    nsweep_dict[38] = "boundu(129)"
    nsweep_dict[39] = "boundu(131)"
    nsweep_dict[40] = "boundu(135)"
    nsweep_dict[41] = "blnkoth"
    nsweep_dict[42] = "fimp(9)"
    nsweep_dict[43] = "rho_ecrh"
    nsweep_dict[44] = "alstrtf"
    nsweep_dict[45] = "tmargmin_tf"
    nsweep_dict[46] = "boundu(152)"
    nsweep_dict[47] = "impurity_enrichment(9)"
    nsweep_dict[48] = "n_pancake"
    nsweep_dict[49] = "n_layer"
    nsweep_dict[50] = "fimp(13)"
    nsweep_dict[51] = "ftar"
    nsweep_dict[52] = "rad_fraction_sol"
    nsweep_dict[54] = "b_crit_upper_nbti"
    nsweep_dict[55] = "shldith"
    nsweep_dict[56] = "crypmw_max"
    nsweep_dict[57] = "bt"  # Genuinly bt lower bound
    nsweep_dict[58] = "scrapli"
    nsweep_dict[59] = "scraplo"
    nsweep_dict[60] = "sig_tf_wp_max"
    nsweep_dict[61] = "copperaoh_m2_max"
    nsweep_dict[62] = "coheof"
    nsweep_dict[63] = "ohcth"
    nsweep_dict[64] = "ohhghf"
    nsweep_dict[65] = "csfv.n_cycle_min"
    nsweep_dict[66] = "pfv.oh_steel_frac"
    nsweep_dict[67] = "csfv.t_crack_vertical"
    nsweep_dict[77] = "fvs"  # actaully lower bound fvs
    nsweep_dict[78] = "vburn"
    nsweep_dict[79] = "rplas"  # "plasma_res_factor"
    # -------------------

    # Getting the scanned variable name
    m_file = mf.MFile(filename=input_file)
    nsweep_ref = int(m_file.data["nsweep"].get_scan(-1))
    try:
        scan_var_name = nsweep_dict[nsweep_ref]
    except:
        scan_var_name = "Null"

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
    nn = 0
    isweep = int(m_file.data["isweep"].get_scan(-1))
    scan_points = np.zeros(num_converged_sol)
    for ii in range(isweep):
        ifail = m_file.data["ifail"].get_scan(ii + 1)
        if ifail == 1:
            scan_points[nn] = m_file.data[scan_var_name].get_scan(ii + 1)
            nn += 1

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

    ind = [y for y, _ in enumerate(scan_points)]
    if args.inboard:
        end_scan = radial_labels.index("Plasma")
    else:
        end_scan = len(radial_build)

    for kk in range((len(radial_build[:end_scan, 0]))):
        if kk == 0:
            lower = np.zeros(len(radial_build[kk, :]))
        else:
            lower = lower + radial_build[kk - 1, :]
        plt.barh(
            ind,
            radial_build[kk, :],
            left=lower,
            height=0.8,
            label=radial_labels[kk],
            color=radial_color[kk],
            edgecolor="black",
            linewidth=0.05,
        )

    plt.yticks(ind, scan_points, fontsize=axis_tick_size)
    plt.ylabel(labels[scan_var_name], fontsize=axis_font_size)
    plt.legend(
        bbox_to_anchor=(0.5, -0.15), loc="upper center", fontsize=legend_size, ncol=4
    )
    plt.xlabel("Radius [m]")
    plt.tight_layout()
    plt.savefig(
        "{}.{}".format(
            output_name,
            save_format,
        ),
        bbox_inches="tight",
    )

    # Display plot (used in Jupyter notebooks)
    plt.show()
    plt.clf()


if __name__ == "__main__":
    main()
