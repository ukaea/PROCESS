#! /usr/bin/env python
"""
Utility to create csv-file for PROCESS runs database.
THIS SHOULD NOT BE EDITED WITHOUT PRIOR CONSULTATION OF FRANCESCO MAVIGLIA
WHO MAINTAINS THE PRCOESS RUNS DATABASE!!!!

Input: MFILE.DAT
Output: process_summary.txt
Optional arguments:
        # change the input file name
        create_csv4database.py -f MFILE.DAT
        # change the output file name
        create_csv4database.py -o out.txt

Based on plot_proc.py by James Morris
Initial version output_data.py by Richard Kemp,  08/10/2014
16/12/16 HL: separated utility to conserve format for PROCESS runs database
"""


import argparse
import process.io.mfile as mf
from dict_tools import proc_dict


def write_data(data, mfile_data, outfile, scan):
    """Function to write output to file.

    Arguments:
      data --> list of data to write
      m_file_data --> MFILE.DAT data to read
      outfile --> output file reference
      scan --> scan to read from MFILE.DAT

    """
    for i in range(len(data)):
        if isinstance(data[i][0], str):
            if mfile_data.data[data[i][0]].exists:
                dat = mfile_data.data[data[i][0]].get_scan(scan)
                if isinstance(dat, str):
                    value = dat
                else:
                    value = "%.4g" % mfile_data.data[data[i][0]].get_scan(scan)
                if "alpha" in data[i][0]:
                    value = str(float(value) + 1.0)
                if isinstance(dat, str):
                    outfile.write('"' + data[i][1] + '",   "' + value + '"\n')
                else:
                    outfile.write('"' + data[i][1] + '",   ' + value + "\n")
            else:
                mfile_data.data[data[i][0]].get_scan(-1)
                # Write the data label even if the data is not available.
                outfile.write('"' + data[i][1] + '",   ' + "value missing!" + "\n")
                # outfile.write("value missing!\n")
        else:
            dat = data[i][0]
            if isinstance(dat, str):
                value = dat
            else:
                value = "%.4g" % data[i][0]
            outfile.write('"' + data[i][1] + '",   ' + value + "\n")


def main(mfile_data, output_file, scan=-1):
    """Function to output summary data for insertion into spreadsheet.

    Arguments:
      m_file_data --> MFILE.DAT data to read
      output_file --> file to write to
      scan --> scan to read from MFILE.DAT

    """

    general = [
        ("procver", "PROCESS_Version"),
        ("date", "Date:"),
        ("time", "Time:"),
        ("username", "User:"),
    ]

    in_blanket_thk = mfile_data.data["shldith"].get_scan(scan) + mfile_data.data[
        "blnkith"
    ].get_scan(scan)
    out_blanket_thk = mfile_data.data["shldoth"].get_scan(scan) + mfile_data.data[
        "blnkoth"
    ].get_scan(scan)

    geometry = [
        ("rmajor", "R_0"),
        ("rminor", "a"),
        ("aspect", "A"),
        ("kappa95", "kappa_95"),
        ("triang95", "delta_95"),
        ("sarea", "Surface area"),
        ("vol", "Plasma volume"),
        ("n_tf", "No. of TF coils"),
        (in_blanket_thk, "i/b blkt/shld"),
        (out_blanket_thk, "o/b blkt/shld"),
        ("powfmw", "Fusion power"),
    ]

    nong = mfile_data.data["dnla"].get_scan(scan) / mfile_data.data[
        "dlimit(7)"
    ].get_scan(scan)
    dnz = mfile_data.data["dnz"].get_scan(scan) / mfile_data.data["dene"].get_scan(scan)
    tepeak = mfile_data.data["te0"].get_scan(scan) / mfile_data.data["te"].get_scan(
        scan
    )

    nepeak = mfile_data.data["ne0"].get_scan(scan) / mfile_data.data["dene"].get_scan(
        scan
    )

    physics = [
        ("plascur/1d6", "I_p"),
        ("bt", "Vacuum B_T"),
        ("q95", "q95"),
        ("normalised thermal beta", "beta_N, thermal"),
        ("normalised total beta", "beta_N, total"),
        ("thermal poloidal beta", "beta_P, thermal"),
        ("betap", "beta_P, total"),
        ("te", "<t_e>"),
        ("dene", "<n_e>"),
        (nong, "<n_e>/n_G"),
        (tepeak, "T_e0/<T_e>"),
        (nepeak, "n_e0/<n_e>"),
        ("zeff", "Z_eff"),
        ("zeffso", "Z_eff, SoL"),
        (dnz, "n_Z/n_e"),
        ("taueff", "tau_e"),
        ("hfact", "H-factor"),
        ("tauelaw", "Scaling law"),
    ]

    #    dnla = mfile_data.data["dnla"].get_scan(scan)/1.0e20
    #    bt = mfile_data.data["bt"].get_scan(scan)
    #    surf = mfile_data.data["sarea"].get_scan(scan)
    pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)
    #    err = pthresh - mfile_data.data["pthrmw(7)"].get_scan(scan)
    gross_eff = 100.0 * (
        mfile_data.data["pgrossmw"].get_scan(scan)
        / mfile_data.data["pthermmw"].get_scan(scan)
    )
    net_eff = 100.0 * (
        (
            mfile_data.data["pgrossmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
        / (
            mfile_data.data["pthermmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
    )
    plant_eff = 100.0 * (
        mfile_data.data["pnetelmw"].get_scan(scan)
        / mfile_data.data["powfmw"].get_scan(scan)
    )

    power_flows = [
        ("wallmw", "Av. neutron wall load"),
        ("pinnerzoneradmw", "inner zone radiation"),
        ("psyncpv*vol", "Synchrotron radiation"),
        ("pouterzoneradmw", "outer zone radiation"),
        ("pnucblkt", "Nuclear heating in blanket"),
        ("pnucshld", "Nuclear heating in shield"),
        ("pdivt", "Psep / Pdiv"),
        (pthresh, "H-mode threshold (M=2.5)"),
        ("fwbllife", "FW/Blanket life"),
        ("divlife", "Divertor life"),
        ("pthermmw", "Thermal Power"),
        (gross_eff, "Gross cycle efficiency"),
        (net_eff, "Net cycle efficiency"),
        ("pgrossmw", "Gross electric power"),
        ("pnetelmw", "Net electric power"),
        (plant_eff, "Plant efficiency Pe/Pfus"),
    ]

    pinjie = mfile_data.data["pinjmw"].get_scan(scan)
    pdivt = mfile_data.data["pdivt"].get_scan(scan)
    pdivr = pdivt / mfile_data.data["rmajor"].get_scan(scan)
    pdivnr = (
        10.0e20
        * mfile_data.data["pdivt"].get_scan(scan)
        / (
            mfile_data.data["rmajor"].get_scan(scan)
            * mfile_data.data["dene"].get_scan(scan)
        )
    )
    #    dnla = mfile_data.data["dnla"].get_scan(scan)/1.0e20
    #    bt = mfile_data.data["bt"].get_scan(scan)
    #    surf = mfile_data.data["sarea"].get_scan(scan)
    pthresh = mfile_data.data["pthrmw(6)"].get_scan(scan)
    flh = pdivt / pthresh
    powerht = mfile_data.data["powerht"].get_scan(scan)
    psync = mfile_data.data["psyncpv*vol"].get_scan(scan)
    pbrem = mfile_data.data["pinnerzoneradmw"].get_scan(scan)
    hfact = mfile_data.data["hfact"].get_scan(scan)
    hstar = hfact * (powerht / (powerht + psync + pbrem)) ** 0.31

    current_drive = [
        (pinjie, "SS auxiliary power"),
        ("pheat", "Power for heating only"),
        ("bootipf", "Bootstrap fraction"),
        ("faccd", "Auxiliary fraction"),
        ("facoh", "Ohmic fraction"),
        ("gamnb", "NB gamma"),
        ("enbeam", "NB energy"),
        ("powerht", "Assumed heating power"),
        (pdivr, "Pdiv/R"),
        (pdivnr, "Pdiv/(n R)"),
        (flh, "Pdiv/PLH"),
        (hstar, "H* (non-rad. corr.)"),
    ]

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile_data.data.keys():
        if "rpf(" in item:
            number_of_coils += 1
    pf_info = []
    for i in range(1, number_of_coils):
        if i % 2 != 0:
            pf_info.append(
                (
                    mfile_data.data["ric(%s)" % str(i).zfill(2)].get_scan(scan),
                    "PF %s" % str(i),
                )
            )
    tburn = mfile_data.data["tburn"].get_scan(scan) / 3600.0
    #    tftype = proc_dict.DICT_TF_TYPE[mfile_data.data["i_tf_sc_mat"].get_scan(scan)]
    vssoft = mfile_data.data["vsres"].get_scan(scan) + mfile_data.data[
        "vsind"
    ].get_scan(scan)

    coil_currents = [
        (pf_info[0][0], pf_info[0][1]),
        (pf_info[1][0], pf_info[1][1]),
        (pf_info[2][0], pf_info[2][1]),
        (vssoft, "Startup flux swing"),
        ("vstot", "Available flux swing"),
        (tburn, "Burn time"),
        ("bmaxtf", "Peak field at conductor"),
        ("iooic", r"I/I$_{\mathrm{crit}}$"),
        ("tmarg", "Temperature margin"),
        ("sig_tf_case", "TF case maximum shear stress (Tresca criterion)"),
        ("sig_tf_wp", "TF conduit maximum shear stress (Tresca criterion)"),
        (
            "sig_tf_case_max",
            "Allowable maximum shear stress in the TF case (Tresca criterion)",
        ),
        (
            "sig_tf_wp_max",
            "Allowable maximum shear stress in the TF conduit (Tresca criterion)",
        ),
    ]

    costs = [
        ("coe", "Cost of electricity"),
        ("concost", "Constructed cost"),
        ("capcost", "Total capex"),
    ]

    # open file for writing
    ofile = open(output_file, "w")

    opt = proc_dict.DICT_OPTIMISATION_VARS[
        abs(int(mfile_data.data["minmax"].get_scan(scan)))
    ]
    ofile.write('"GENERAL"\n')
    write_data(general, mfile_data, ofile, scan)
    ofile.write('"optimising:",   "' + opt + '"\n')
    ofile.write("\n")

    ofile.write('"GEOMETRY"\n')
    write_data(geometry, mfile_data, ofile, scan)
    ofile.write("\n")

    ofile.write('"PHYSICS"\n')
    write_data(physics, mfile_data, ofile, scan)
    ofile.write("\n")

    ofile.write('"POWER_FLOWS"\n')
    write_data(power_flows, mfile_data, ofile, scan)
    ofile.write("\n")

    ofile.write('"CURRENT_DRIVE"\n')
    write_data(current_drive, mfile_data, ofile, scan)
    ofile.write("\n")

    ofile.write('"COIL_CURRENTS"\n')
    write_data(coil_currents, mfile_data, ofile, scan)
    ofile.write("\n")

    ofile.write('"COSTS"\n')
    write_data(costs, mfile_data, ofile, scan)
    ofile.write("\n")

    # close file
    ofile.close()


if __name__ == "__main__":

    # Setup command line arguments
    PARSER = argparse.ArgumentParser(
        description="Produce a single-column "
        "comma-separated (.txt) summary. "
        "For info contact "
        "rich.kemp@ukaea.ac.uk or "
        "james.morris2@ukaea.ac.uk"
    )

    PARSER.add_argument(
        "-f",
        metavar="FILENAME",
        type=str,
        default="MFILE.DAT",
        help="specify input filename",
    )

    PARSER.add_argument(
        "-o",
        metavar="OUTPUT",
        type=str,
        default="process_summary.txt",
        help="specify output filename",
    )

    ARGS = PARSER.parse_args()

    # read MFILE
    MFILE = mf.MFile(ARGS.f)

    # run main
    main(MFILE, ARGS.o)
