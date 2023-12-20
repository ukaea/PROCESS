"""Python utility for plotting the output of a PROCESS scan.

Depending of the type of scans, different actions will be taken:
1D SCANS: a simple graph using the scanned variable for x axis
and the selected variable on the y axis.
- Any number of output variables can be selected, a plot will be
made for each
- Several inputs files can be used at the same time if the same variable
is scanned. The different runs results will be plotted in the same
graph.
- If several inputs are used, the folder name or the file is used as
a legend

- 2D SCANS: n_scan_1 graph will be plotted using the second scanned variable
as x axis and the selected output as y axis
- Only one 2D scan can be ploted at once.

Performed checks:
- Non converged points are not plotted
- Only outputs existing in the MFILE.DAT are plotted
- A LaTeX label dicts is integrated, there is a check if the requested
variable is set. Otherwise the sad and gloomy PROCESS name is used
- No plot is made if the MFILE does not exists
- If the file is a folder, the contained MFILE is used as an input.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
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
        "--input_files",
        default="MFILE.DAT",
        help=(
            "Specify input file(s) path(s) (default = MFILE.DAT)\n"
            "More than one input file can be used eg: -f 'A_MFILE.DAT "
            "B_MFILE.DAT'.\nYou can only specify the folder containing the "
            "MFILE.\nThe different files scan will be plotted on the same "
            "graph.\nThe scans must use the same scan variation."
        ),
    )

    # At least one output variable must be supplied in order to plot
    parser.add_argument(
        "-yv",
        "--y_vars",
        required=True,
        help=(
            "Select the output variables\nMore than one output can be plotted "
            "eg: -yv 'var1 var2'\nA separate plot will be created for each "
            "inputs"
        ),
    )

    parser.add_argument(
        "-yv2",
        "--y_vars2",
        default="",
        help=(
            "Select the 2nd axis output variable\n "
            "eg: -yv2 'var'\n 2nd variable will be plotted on shared figure "
            "inputs"
        ),
    )

    parser.add_argument(
        "-o",
        "--outputdir",
        default=Path.cwd(),
        help="Output directory for plots, defaults to current working directory.",
    )

    parser.add_argument(
        "-out",
        "--term_output",
        action="store_true",
        help="Option to show scans values on terminal",
    )

    parser.add_argument(
        "-sf",
        "--save_format",
        nargs="?",
        default="pdf",
        help="Output format (default='pdf') ",
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
        "-ln",
        "--label_name",
        default="",
        help=(
            "Label names for plot legend. If multiple input files used then \n"
            "list the same number of label names eg: -nl 'leg1 leg2'\n"
            "(default = MFile file name) "
        ),
    )

    parser.add_argument(
        "-2DC",
        "--two_dimensional_contour",
        action="store_true",
        help=(
            "Option to plot 2D scans as a coloured contour plot instead of a line plot \n  "
            "Note: Non convergent points will show up with a value of zero \n "
            "Note: The scan paramters must both be in increasing orderl \n "
        ),
    )

    parser.add_argument(
        "-stc",
        "--stack_plots",
        action="store_true",
        help=(
            "Option to plot multiple 1D plots in a column of subplots \n  "
            "Variables will be plotted in order of input"
        ),
    )

    return parser.parse_args(args)


def main(args=None):
    """Main plot scans script.

    :param args: optional command-line args from test function, defaults to None
    :type args: list, optional
    """
    args = parse_args(args)

    # Parameters to be used as function input
    # ---------------------------------------
    input_files = str(args.input_files)
    output_names = str(args.y_vars)
    output_names2 = str(args.y_vars2)
    save_format = str(args.save_format)
    term_output = args.term_output
    label_name = str(args.label_name)
    two_dimensional_contour = args.two_dimensional_contour
    stack_plots = args.stack_plots
    # ---------------------------------------

    # Input checks
    # ------------
    # Formting the inputs
    output_names = output_names.split(" ")
    while "" in output_names:
        output_names.remove("")

    output_names2 = output_names2.split(" ")
    while "" in output_names2:
        output_names2.remove("")

    input_files = input_files.split(" ")
    while "" in input_files:
        input_files.remove("")

    label_name = label_name.split(" ")
    while "" in label_name:
        label_name.remove("")

    # If the input file is a directory, add MFILE.DAT
    for ii in range(len(input_files)):
        if os.path.isdir(input_files[ii]):
            input_files[ii] = input_files[ii].replace("/", "")
            input_files[ii] = input_files[ii] + str("/MFILE.DAT")

        # Check for the existence of the MFILE
        if not os.path.isfile(input_files[ii]):
            print(
                "ERROR : The {} MFILE does not exist, skipping it".format(
                    input_files[ii]
                )
            )
            input_files.remove(input_files[ii])

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
    labels["te"] = r"$\left< T_\mathrm{e} \right> \ [keV]$"
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
    labels["sig_tf_wp_max"] = r"$\sigma_{TP,wp}^{max}$"
    labels["rli"] = r"$l_i$"
    labels["n_cycle_min"] = r"$MinCycles_{\mathrm{Stress.min}}^{\mathrm{CS}}$"
    labels["n_cycle"] = r"$Cycles_{\mathrm{Stress}}^{\mathrm{CS}}$"
    labels["a_oh_turn"] = r"$Turn_{\mathrm{area}}^{\mathrm{CS}}[$m$^{2}]$"
    labels["tbrnmn"] = r"$t_{\mathrm{burn.min}}$[$s$]"
    labels["pfv.oh_steel_frac"] = r"$f_{\mathrm{Steel}}^{\mathrm{CS}}$"
    labels["csfv.t_structural_radial"] = r"$Turn_{\mathrm{radial}}^{\mathrm{CS}}[$m$]$"
    labels["csfv.t_crack_vertical"] = r"$Crack_{\mathrm{vertical}}^{\mathrm{CS}}[$m$]$"
    labels["inlet_temp_liq"] = r"Breeder/coolant inlet T [K]"
    labels["outlet_temp_liq"] = r"Breeder/coolant outlet T [K]"
    labels["blpressure_liq"] = r"Breeder/coolant pressure [Pa]"
    labels["n_liq_recirc"] = r"Breeder/coolant recirculations"
    labels["bz_channel_conduct_liq"] = r"channel wall conductivity [A V-1 m-1]"
    labels["pnuc_fw_ratio_dcll"] = r"FW nuclear power fraction"
    labels["f_nuc_pow_bz_struct"] = r"BZ structure fraction"
    labels["pitch"] = r"FW pitch [m]"
    labels["coe"] = r"$\mathrm{LCOE}$ [$m\$/kWh$]"
    labels["beta"] = r"$\beta$"
    labels["fimp(13)"] = r"$Xe_{\mathrm{f}}$"
    labels["pdivmax/rmajor"] = r"$P_{\mathrm{div}}/R_\mathrm{maj}$ [MW/m]"
    labels["etath"] = r"Thermal to Electric efficiency"
    labels["fkind"] = r"N$^\mathrm{th}$ of a kind factor"
    labels["startupratio"] = r"Gyrotron Redundancy"
    labels["etaech"] = r"ECH wall plug to injector efficiency"
    labels["tauee"] = r"$\tau_E$"
    labels["dene"] = r"$n_e$"

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
    nsweep_dict[30] = ""  # OBSOLETE
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
    nsweep_dict[43] = "Obsolete"  # Removed
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
    nsweep_dict[68] = "inlet_temp_liq"
    nsweep_dict[69] = "outlet_temp_liq"
    nsweep_dict[70] = "blpressure_liq"
    nsweep_dict[71] = "n_liq_recirc"
    nsweep_dict[72] = "bz_channel_conduct_liq"
    nsweep_dict[73] = "pnuc_fw_ratio_dcll"
    nsweep_dict[74] = "f_nuc_pow_bz_struct"
    nsweep_dict[75] = "pitch"
    nsweep_dict[76] = "etath"
    nsweep_dict[77] = "startupratio"
    nsweep_dict[78] = "fkind"
    nsweep_dict[79] = "etaech"
    # -------------------

    # Getting the scanned variable name
    m_file = mf.MFile(filename=input_files[-1])
    nsweep_ref = int(m_file.data["nsweep"].get_scan(-1))
    scan_var_name = nsweep_dict[nsweep_ref]
    # Get the eventual second scan variable
    nsweep_2_ref = int(0)
    is_2D_scan = False
    scan_2_var_name = str()
    if "nsweep_2" in m_file.data.keys():
        is_2D_scan = True
        nsweep_2_ref = int(m_file.data["nsweep_2"].get_scan(-1))
        scan_2_var_name = nsweep_dict[nsweep_2_ref]

    # Checks
    # ------
    # Check if the nsweep dict has been updated
    if nsweep_ref > len(nsweep_dict) + 1:
        print("ERROR : nsweep = {} not supported by the utility".format(nsweep_ref))
        print("ERROR : Please update the 'nsweep_dict' dict")
        exit()

    # Check if the scan variable is present in the
    if scan_var_name not in m_file.data.keys():
        print("ERROR : `{}` does not exist in PROCESS dicts".format(scan_var_name))
        print("ERROR : The scan variable is probably an upper/lower boundary")
        print("ERROR : Please modify 'nsweep_dict' dict with the constrained var")
        exit()

    # Check if the (first) scan variable LaTeX label is set
    if scan_var_name not in labels:
        print(
            "WARNING: The {} variable LaTeX label is not defined".format(scan_var_name)
        )
        print("WARNING: Please update the 'label' dict")
        labels[scan_var_name] = scan_var_name

    if is_2D_scan:
        # Check if the second scan variable is present in the
        if scan_2_var_name not in m_file.data.keys():
            print(
                "ERROR : `{}` does not exist in PROCESS dicts".format(scan_2_var_name)
            )
            print("ERROR : The scan variable is probably an upper/lower boundary")
            print("ERROR : Please modify 'nsweep_dict' dict with the constrained var")
            exit()

        # Check if the second scan variable LaTeX label is set
        if scan_2_var_name not in labels:
            print("The {} variable LaTeX label is not defined".format(scan_2_var_name))
            print("Please update the 'label' dict")
            labels[scan_var_name] = scan_var_name

    # Only one imput must be used for a 2D scan
    if is_2D_scan and len(input_files) > 1:
        print("ERROR : Only one input file can be used for 2D scans")
        print("ERROR : Exiting")
        exit()
    # ------

    # Plot settings
    # -------------
    # Plot cosmetic settings
    axis_tick_size = 16
    legend_size = 12
    axis_font_size = 18
    # -------------

    # Case of a set of 1D scans
    # ----------------------------------------------------------------------------------------------
    if not is_2D_scan:
        # Loop over the MFILEs
        output_arrays = dict()
        output_arrays2 = dict()
        scan_var_array = dict()
        for input_file in input_files:
            # Opening the MFILE.DAT
            m_file = mf.MFile(filename=input_file)

            # Check if the the scan variable is the same for all inputs
            # ---
            # Same scan var
            nsweep = int(m_file.data["nsweep"].get_scan(-1))
            if nsweep != nsweep_ref:
                print("ERROR : You must use inputs files with the same scan variables")
                print("ERROR : Exiting")
                exit()

            # No D scans
            if "nsweep_2" in m_file.data.keys():
                print("ERROR : You cannot mix 1D with 2D scans")
                print("ERROR : Exiting")
                exit()
            # ---

            # Only selecting the scans that has converged
            # ---
            # Number of scan points
            n_scan = int(m_file.data["isweep"].get_scan(-1))

            # Converged indexes
            conv_i = list()
            for ii in range(n_scan):
                ifail = m_file.data["ifail"].get_scan(ii + 1)
                if ifail == 1:
                    conv_i.append(ii + 1)
                else:
                    failed_value = m_file.data[scan_var_name].get_scan(ii + 1)
                    print(
                        "Warning : Non-convergent scan point : {} = {}".format(
                            scan_var_name, failed_value
                        )
                    )
                    print("Warning : This point will not be shown.")

            # Updating the number of scans
            n_scan = len(conv_i)
            # ---
            # Scanned variable
            scan_var_array[input_file] = np.zeros(n_scan)

            for ii in range(n_scan):
                scan_var_array[input_file][ii] = m_file.data[scan_var_name].get_scan(
                    conv_i[ii]
                )
            # output list declaration
            output_arrays[input_file] = dict()
            output_arrays2[input_file] = dict()
            # First variable scan
            for output_name in output_names:
                ouput_array = np.zeros(n_scan)

                # Check if the output variable exists in the MFILE
                if output_name not in m_file.data.keys():
                    print(
                        "Warning : `{}` does not exist in PROCESS dicts".format(
                            output_name
                        )
                    )
                    print("Warning : `{}` will not be output".format(output_name))
                    continue

                # Check if the output LaTeX variable label exist
                if output_name not in labels:
                    print(
                        "Warning : The {} variable LaTeX label is not defined".format(
                            output_name
                        )
                    )
                    print("Warning : Please update the 'label' dict")
                    labels[output_name] = output_name

                for ii in range(n_scan):
                    ouput_array[ii] = m_file.data[output_name].get_scan(conv_i[ii])
                output_arrays[input_file][output_name] = ouput_array
            # Second variable scan
            if output_names2 != []:
                for output_name2 in output_names2:
                    ouput_array2 = np.zeros(n_scan)

                    # Check if the output variable exists in the MFILE
                    if output_name2 not in m_file.data.keys():
                        print(
                            "Warning : `{}` does not exist in PROCESS dicts".format(
                                output_name2
                            )
                        )
                        print("Warning : `{}` will not be output".format(output_name2))
                        continue

                    # Check if the output LaTeX variable label exist
                    if output_name2 not in labels:
                        print(
                            "Warning : The {} variable LaTeX label is not defined".format(
                                output_name2
                            )
                        )
                        print("Warning : Please update the 'label' dict")
                        labels[output_name2] = output_name2

                    for ii in range(n_scan):
                        ouput_array2[ii] = m_file.data[output_name2].get_scan(
                            conv_i[ii]
                        )
                    output_arrays2[input_file][output_name2] = ouput_array2
            # Terminal output
            if term_output:
                print(f"\n{input_file} scan output\n")
                print("X-axis:")
                print(f"scan var {scan_var_name} : {scan_var_array[input_file]}\n")
                print("Y-axis:")
                for output_name in output_names:
                    # Check if the output variable exists in the MFILE
                    if output_name not in m_file.data.keys():
                        continue

                    print(f"{output_name} : {output_arrays[input_file][output_name]}")
                print()
                if output_names2 != []:
                    print(
                        (
                            "Y2-Axis\n "
                            + f" {output_name2} : {output_arrays2[input_file][output_name2]}\n"
                        )
                    )
        # Plot section
        # -----------
        if stack_plots:
            fig, axs = plt.subplots(
                len(output_names),
                1,
                figsize=(8.0, (3.5 + (1 * len(output_names)))),
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.0)
        for output_name in output_names:
            # reset counter for label_name
            kk = 0

            # Check if the output variable exists in the MFILE
            if output_name not in m_file.data.keys():
                continue

            # Loop over inputs
            for input_file in input_files:
                # Legend label formating
                if label_name == []:
                    labl = input_file
                    if "/MFILE.DAT" in input_file:
                        labl = input_file[:-10]
                    elif "MFILE.DAT" in input_file:
                        labl = input_file[:-9]
                    labl = labl.replace("_", " ")
                else:
                    labl = label_name[kk]
                    kk = kk + 1

                # Plot the graph
                if output_names2 != [] and not stack_plots:
                    fig, ax = plt.subplots()
                    ax.plot(
                        scan_var_array[input_file],
                        output_arrays[input_file][output_name],
                        "--o",
                        color="blue",
                        label=labl,
                    )
                else:
                    if stack_plots:
                        axs[output_names.index(output_name)].plot(
                            scan_var_array[input_file],
                            output_arrays[input_file][output_name],
                            "--o",
                            color="blue" if output_names2 != [] else None,
                            label=labl,
                        )
                    else:
                        plt.plot(
                            scan_var_array[input_file],
                            output_arrays[input_file][output_name],
                            "--o",
                            color="blue" if output_names2 != [] else None,
                            label=labl,
                        )
                        plt.xticks(size=axis_tick_size)
                        plt.yticks(size=axis_tick_size)
                if output_names2 != []:
                    ax2 = ax.twinx()
                    ax2.plot(
                        scan_var_array[input_file],
                        output_arrays2[input_file][output_name2],
                        "--o",
                        color="red",
                        label=labl,
                    )
                    ax2.set_ylabel(
                        labels[output_name2],
                        fontsize=axis_font_size,
                        color="red",
                    )
            if output_names2 != []:
                ax2.yaxis.grid(True)
                ax.xaxis.grid(True)
                ax.set_ylabel(
                    labels[output_name], fontsize=axis_font_size, color="blue"
                )
                ax.set_xlabel(labels[scan_var_name], fontsize=axis_font_size)
            elif stack_plots:
                axs[output_names.index(output_name)].minorticks_on()
                axs[output_names.index(output_name)].grid(True)
                axs[output_names.index(output_name)].set_ylabel(
                    labels[output_name],
                )

                plt.xlabel(labels[scan_var_name], fontsize=axis_font_size)
                if len(input_files) > 1:
                    plt.legend(
                        loc="lower center",
                        fontsize=legend_size,
                        bbox_to_anchor=(0.5, -0.5 - (0.1 * len(output_names))),
                        # bbox_to_anchor=(0.5, -1.4),
                        fancybox=True,
                        shadow=False,
                        ncol=len(input_files),
                        columnspacing=0.8,
                    )
                plt.tight_layout()
                axs[output_names.index(output_name)].set_ylim(
                    axs[output_names.index(output_name)].get_ylim()[0]
                    * 0.9,  # Adds some spacing around the plots
                    axs[output_names.index(output_name)].get_ylim()[1] * 1.1,
                )
            else:
                plt.grid(True)
                plt.ylabel(
                    labels[output_name],
                    fontsize=axis_font_size,
                    color="red" if output_names2 != [] else "black",
                )
                plt.xlabel(labels[scan_var_name], fontsize=axis_font_size)
                plt.title(
                    f"{labels[output_name]} vs {labels[scan_var_name]}",
                    fontsize=axis_font_size,
                )
                plt.tight_layout()
                if len(input_files) != 1:
                    plt.legend(loc="best", fontsize=legend_size)
                    plt.tight_layout()

            # Output file naming
            if output_name == "plascur/1d6":
                plt.savefig(
                    f"{args.outputdir}/scan_{scan_var_name}_vs_plascur"
                    + f"_vs_{output_name2}"
                    if output_names2 != []
                    else f"{args.outputdir}/scan_{scan_var_name}_vs_plascur"
                    + f".{save_format}"
                )
            elif output_name == "pdivtbt/qar":
                plt.savefig(
                    f"{args.outputdir}/scan_{scan_var_name}_vs_pdivtbtqar"
                    + f"_vs_{output_name2}"
                    if output_names2 != []
                    else "" + f".{save_format}"
                )
            elif stack_plots and output_names[-1] == output_name:
                plt.savefig(
                    f"{args.outputdir}/scan_{scan_var_name}_vs_{output_name}"
                    + f"_vs_{output_name2}"
                    if output_names2 != []
                    else f"{args.outputdir}/scan_{scan_var_name}_vs_"
                    + "_vs_".join(output_names)
                    + f".{save_format}",
                    dpi=300,
                )

            else:
                plt.savefig(
                    f"{args.outputdir}/scan_{scan_var_name}_vs_{output_name}"
                    + f"_vs_{output_name2}"
                    if output_names2 != []
                    else f"{args.outputdir}/scan_{scan_var_name}_vs_{output_name}"
                    + f".{save_format}",
                    dpi=300,
                )
            if not stack_plots:  # Display plot (used in Jupyter notebooks)
                plt.show()
                plt.clf()
        # ------------

    # In case of a 2D scan
    # ----------------------------------------------------------------------------------------------
    else:
        # Opening the MFILE.DAT
        m_file = mf.MFile(filename=input_files[0])

        # Number of scan points
        n_scan_1 = int(m_file.data["isweep"].get_scan(-1))
        n_scan_2 = int(m_file.data["isweep_2"].get_scan(-1))
        # Selecting the converged runs only
        contour_conv_ij = []  # List of non-converged scan point numbers
        conv_ij = (
            list()
        )  # 2D array of converged scan point numbers (sweep = rows, sweep_2 = columns)
        ii_jj = 0
        for ii in range(n_scan_1):
            conv_ij.append(list())
            for jj in range(n_scan_2):
                ii_jj += 1  # Represents the scan point number in the MFILE
                ifail = m_file.data["ifail"].get_scan(ii_jj)
                if ifail == 1:
                    conv_ij[ii].append(
                        ii_jj
                    )  # Only appends scan number if scan converged
                    contour_conv_ij.append(ii_jj)
                else:
                    failed_value_1 = m_file.data[scan_var_name].get_scan(ii_jj)
                    failed_value_2 = m_file.data[scan_2_var_name].get_scan(ii_jj)
                    print(
                        f"Warning : Non-convergent scan point : ({scan_var_name},{scan_2_var_name}) = ({failed_value_1},{failed_value_2})"
                    )
                    print("Warning : This point will not be shown.")

        # Looping over requested outputs
        for output_name in output_names:
            # Check if the output variable exists in the MFILE
            if output_name not in m_file.data.keys():
                print(f"Warning : `{output_name}` does not exist in PROCESS dicts")
                print(f"Warning : `{output_name}` will not be output")
                continue

            # Check if the output LaTeX variable label exist
            if output_name not in labels:
                print(
                    f"Warning : The {output_name} variable LaTeX label is not defined"
                )
                print("Warning : Please update the 'label' dict")
                labels[output_name] = output_name

            # Declaring the outputs
            output_arrays = list()

            if two_dimensional_contour:
                output_contour_z = np.zeros((n_scan_1, n_scan_2))
                x_contour = []
                y_contour = []
                for i in range(n_scan_2):
                    x_contour.append(m_file.data[scan_2_var_name].get_scan(i + 1))
                for i in range(1, n_scan_1 * n_scan_2, n_scan_2):
                    y_contour.append(
                        m_file.data[scan_var_name].get_scan(i + 1)
                    )  # is the separte lists in the list
                for i in contour_conv_ij:
                    output_contour_z[((i - 1) // n_scan_2)][
                        ((i - 1) % n_scan_2)
                        if ((i - 1) // n_scan_2) % 2 == 0
                        else (-((i - 1) % n_scan_2) - 1)
                    ] = m_file.data[output_name].get_scan(i)

                flat_output_z = output_contour_z.flatten()
                flat_output_z.sort()
                plt.contourf(
                    x_contour,
                    y_contour,
                    output_contour_z,
                    levels=np.linspace(
                        list(filter(lambda i: i > 0.0, flat_output_z))[0],
                        flat_output_z.max(),
                        50,
                    ),
                )

                plt.colorbar().set_label(label=labels[output_name], size=axis_font_size)
                plt.ylabel(labels[scan_var_name], fontsize=axis_font_size)
                plt.xlabel(labels[scan_2_var_name], fontsize=axis_font_size)
                plt.tight_layout()
                plt.savefig(
                    f"{args.outputdir}/scan_{output_name}_vs_{scan_var_name}_{scan_2_var_name}.{save_format}"
                )
                plt.grid(True)
                plt.show()
                plt.clf()

            else:
                # Converged indexes, for normal 2D line plot
                for (
                    conv_j
                ) in (
                    conv_ij
                ):  # conv_j is an array element containing the converged scan numbers
                    # Scanned variables
                    scan_1_var_array = np.zeros(len(conv_j))
                    scan_2_var_array = np.zeros(len(conv_j))
                    output_array = np.zeros(len(conv_j))
                    for jj in range(len(conv_j)):
                        scan_1_var_array[jj] = m_file.data[scan_var_name].get_scan(
                            conv_j[jj]
                        )
                        scan_2_var_array[jj] = m_file.data[scan_2_var_name].get_scan(
                            conv_j[jj]
                        )
                        output_array[jj] = m_file.data[output_name].get_scan(conv_j[jj])

                    # Label formating
                    labl = f"{labels[scan_var_name]} = {scan_1_var_array[0]}"

                    # Plot the graph
                    plt.plot(scan_2_var_array, output_array, "--o", label=labl)

                plt.grid(True)
                plt.ylabel(labels[output_name], fontsize=axis_font_size)
                plt.xlabel(labels[scan_2_var_name], fontsize=axis_font_size)
                plt.legend(loc="best", fontsize=legend_size)
                plt.xticks(size=axis_tick_size)
                plt.yticks(size=axis_tick_size)
                plt.tight_layout()
                plt.savefig(
                    f"{args.outputdir}/scan_{output_name}_vs_{scan_var_name}_{scan_2_var_name}.{save_format}"
                )

                # Display plot (used in Jupyter notebooks)
                plt.show()
                plt.clf()


if __name__ == "__main__":
    main()
