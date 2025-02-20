"""

Python tool for comparing MFILE and outputting differences.
The tool does not work for MFiles that are not the result of
a full PROCESS run (ie if an error or exception occured).

James Morris
14/04/15

CCFE

Notes:
    + 24/11/2021: Global dictionary variables moved within the functions
                to avoid cyclic dependencies. This is because the dicts
                generation script imports, and inspects, process.

"""

import argparse
import sys

import numpy as np
from numpy import isfinite

import process.io.mfile as mf
from process.io.python_fortran_dicts import get_dicts

# Dictionary for parameter descriptions
DICT_DESCRIPTIONS = get_dicts()["DICT_DESCRIPTIONS"]

DEFAULT_COMPARE_PARAMS = [
    "rmajor",
    "rminor",
    "aspect",
    "kappa",
    "kappa95",
    "triang",
    "triang95",
    "fimp(01",
    "fimp(02",
    "fimp(03",
    "fimp(04",
    "fimp(05",
    "fimp(06",
    "fimp(07",
    "fimp(08",
    "fimp(09",
    "fimp(10",
    "fimp(11",
    "fimp(12",
    "fimp(13",
    "fimp(14",
    "a_plasma_surface",
    "vol_plasma",
    "n_tf_coils",
    "dr_shld_inboard",
    "dr_shld_outboard",
    "dr_blkt_inboard",
    "dr_blkt_outboard",
    "fusion_power",
    "plasma_current_MA",
    "bt",
    "q95",
    "beta_poloidal",
    "te",
    "dene",
    "hfact",
    "vs_cs_pf_total_pulse",
    "bt",
    "bmaxtfrp",
    "tmarg",
    "iooic",
    "sig_tf_case",
    "sig_tf_wp",
    "sig_tf_case_max",
    "sig_tf_wp_max",
    "pgrossmw",
    "htpmw",
    "pnetelmw",
    "wallmw",
    "f_nd_alpha_electron",
    "p_plasma_inner_rad_mw",
    "p_plasma_rad_mw",
    "pnucblkt",
    "pnucshld",
    "pdivt",
    "pheat",
    "bootstrap_current_fraction",
    "aux_current_fraction",
    "inductive_current_fraction",
    "gamnb",
    "beam_energy",
    "p_plasma_loss_mw",
]

BASELINE_LIST = [
    "procver",
    "time",
    "username",
    "tagno",
    "commsg",
    "ifail",
    "rmajor",
    "rminor",
    "aspect",
    "kappa",
    "kappa95",
    "triang",
    "triang95",
    "a_plasma_surface",
    "vol_plasma",
    "n_tf_coils",
    "fusion_power",
    "plasma_current_MA",
    "bt",
    "q95",
    "beta",
    "normalised_thermal_beta",
    "beta_norm_total",
    "thermal_beta",
    "thermal_poloidal_beta",
    "te",
    "te0",
    "dene",
    "ne0",
    "dnla_gw",
    "tesep",
    "nesep",
    "teped",
    "neped",
    "zeff",
    "nd_impurities",
    "t_energy_confinement",
    "hfact",
    "tauelaw",
    "f_nd_alpha_electron",
    "wallmw",
    "p_plasma_inner_rad_mw",
    "pden_plasma_sync_mw*vol_plasma",
    "p_plasma_rad_mw",
    "pnucblkt",
    "pnucshld",
    "pdivt",
    "divlife",
    "pthermmw",
    "dr_bore",
    "dr_cs",
    "dr_cs_precomp",
    "dr_cs_tf_gap",
    "dr_tf_inboard",
    "deltf",
    "dr_shld_thermal_inboard",
    "dr_shld_thermal_outboard",
    "dz_shld_thermal",
    "dr_shld_vv_gap_inboard",
    "dr_vv_inboard",
    "dr_vv_outboard",
    "dz_vv_upper",
    "d_vv_bot",
    "dr_shld_inboard",
    "dr_shld_blkt_gap",
    "dr_blkt_inboard",
    "dr_fw_inboard",
    "dr_fw_plasma_gap_inboard",
    "dr_fw_plasma_gap_outboard",
    "dr_fw_outboard",
    "dr_blkt_outboard",
    "dr_shld_outboard",
    "dr_shld_vv_gap_outboard",
    "dr_tf_shld_gap",
    "dr_tf_outboard",
    "etath",
    "pgrossmw",
    "pnetelmw",
    "pinjmw",
    "pheat",
    "bootstrap_current_fraction",
    "aux_current_fraction",
    "inductive_current_fraction",
    "gamnb",
    "beam_energy",
    "p_plasma_loss_mw",
    "pdivt",
    "vssoft",
    "vs_cs_pf_total_pulse",
    "t_burn",
    "bmaxtf",
    "iooic",
    "tmarg",
    "tftmp",
    "qtarget",
    "qtargetcomplete",
    "totalpowerlost",
]

BLANKET_COMPARE_PARAMS = [
    "dr_blkt_inboard",
    "dr_blkt_outboard",
    "fusion_power",
    "pnucblkt",
    "p_fw_nuclear_heat_total_mw",
    "ptfnuc",
    "pnucshld",
    "pnucdiv",
    "tbr",
    "li6enrich",
    "a_fw_total",
    "emult",
]

GENERIC_LIST = [
    "rmajor",
    "rminor",
    "aspect",
    "kappa",
    "kappa95",
    "triang",
    "triang95",
    "fusion_power",
    "plasma_current_MA",
    "bt",
    "q95",
    "beta",
    "te",
    "dene",
    "pinjmw",
    "pnetelmw",
    "wallmw",
    "pflux_fw_rad_mw",
    "f_nd_alpha_electron",
    "p_plasma_inner_rad_mw",
    "p_plasma_rad_mw",
    "bootstrap_current_fraction",
    "pdivmax_over_rmajor",
    "fimp(14",
    "etath",
    "capcost",
    "coe",
    "dr_bore",
    "dr_cs",
    "dr_cs_precomp",
    "dr_cs_tf_gap",
    "dr_tf_inboard",
    "dr_tf_shld_gap",
    "dr_shld_thermal_inboard",
    "dr_shld_thermal_outboard",
    "dz_shld_thermal",
    "dr_shld_vv_gap_inboard",
    "dr_vv_inboard",
    "dr_shld_inboard",
    "dr_shld_blkt_gap",
    "dr_blkt_inboard",
    "dr_fw_inboard",
    "dr_fw_plasma_gap_inboard",
    "dr_fw_plasma_gap_outboard",
    "dr_fw_outboard",
    "dr_blkt_outboard",
    "dr_shld_outboard",
    "dr_vv_outboard",
    "dr_shld_vv_gap_outboard",
    "dr_tf_shld_gap",
    "dr_tf_outboard",
    "dz_xpoint_divertor",
    "dz_divertor",
    "d_vv_bot",
    "dz_shld_lower",
    "dz_shld_vv_gap",
]


class BColors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"


def main(arg):
    """Main function for comparing MFILEs

    :param arg: List of arguments
    :return:
    """

    print_counter = 0
    n = 2
    mfile_list = []
    for item in arg.f:
        mfile = mf.MFile(filename=item)
        if mfile.data["error_status"].get_scan(-1) == 3:
            raise RuntimeError(
                f"{item} is an MFile from a PROCESS run that did not converge"
                " and instead results from an error during the run"
            )

        mfile_list.append(mfile)

    var_list = []
    missing_vars = []
    diff_list = []
    within_list = []

    key_list = mfile_list[0].data.keys()
    for var in key_list:
        store = True
        for f in mfile_list:
            if var not in f.data:
                store = False
                missing_vars.append(var)

        if store:
            var_list.append(var)

    if arg.defaults:
        var_list = DEFAULT_COMPARE_PARAMS

    if arg.blanket:
        var_list = BLANKET_COMPARE_PARAMS

    if arg.baseline:
        var_list = BASELINE_LIST

    if arg.generic:
        var_list = GENERIC_LIST

    for v in var_list:
        if "normres" in v:
            continue

        values = np.zeros(n)  # replaced scipy with numpy

        if v not in get_dicts()["DICT_VAR_TYPE"]:
            try:
                eval(mfile_list[0].data[v].get_scan(-1))
            except NameError:
                pass
            except TypeError:
                for m in range(len(mfile_list)):
                    values[m] = mfile_list[m].data[v].get_scan(-1)
            except SyntaxError:
                pass

        elif (
            get_dicts()["DICT_VAR_TYPE"][v] == "real_variable"
            or get_dicts()["DICT_VAR_TYPE"][v] == "int_variable"
        ):
            for m in range(len(mfile_list)):
                values[m] = mfile_list[m].data[v].get_scan(-1)

        norm_vals = []
        if values[0] != 0 and isfinite(values[0]):
            norm_vals = values / values[0]
        # else:
        #    print(key, values[0])

        if len(norm_vals) >= 1:
            key = v.strip(".").strip(" ")
            des = get_dicts()["DICT_DESCRIPTIONS"].get(key, "-")
            a = norm_vals >= 1.0 + arg.acc / 100.0
            b = norm_vals <= 1.0 - arg.acc / 100.0
            if a[1]:
                diff_list.append(v)
                line = (
                    BColors.ENDC
                    + v
                    + "\t"
                    + des
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + BColors.FAIL
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                wline = (
                    v
                    + "\t"
                    + des
                    + "\t"
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                print(line)
                print_counter += 1
                if arg.save:
                    with open("comp.txt", "a") as ofile:
                        ofile.write(wline + "\n")
            elif b[1]:
                diff_list.append(v)
                line = (
                    BColors.ENDC
                    + v
                    + "\t"
                    + des
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + BColors.FAIL
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                wline = (
                    v
                    + "\t"
                    + des
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                print(line)
                print_counter += 1
                if arg.save:
                    with open("comp.txt", "a") as ofile:
                        ofile.write(wline + "\n")
            else:
                within_list.append(v)
                line = (
                    BColors.ENDC
                    + v
                    + "\t"
                    + des
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                wline = (
                    v
                    + "\t"
                    + des
                    + "\t"
                    + str(values[0])
                    + "\t"
                    + str(values[1])
                    + "\t"
                    + str(round((norm_vals[1] - 1) * 100.0, 2))
                    + " %"
                )
                if arg.verbose:
                    print(line)
                    print_counter += 1
                    with open("comp.txt", "a") as ofile:
                        ofile.write(wline + "\n")

    if arg.baseline and arg.acc >= 10.0:
        if print_counter == 0:
            sys.exit(0)
        else:
            sys.exit(f"Differences in baseline output by more than {arg.acc}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Produce a comparison "
        "between two PROCESS "
        "MFILEs. User Can specify "
        "level of differences to show. "
        "For info contact "
        "james.morris2@ccfe.ac.uk"
    )

    parser.add_argument("-f", metavar="f", type=str, nargs="+", help="Files to compare")

    parser.add_argument(
        "-s", "--save", help="Save output to file called comp.txt", action="store_true"
    )

    parser.add_argument("--acc", type=float, default=5.0)

    parser.add_argument("--verbose", type=float, default=0.0)

    parser.add_argument("--defaults", action="store_true")

    parser.add_argument("--baseline", action="store_true")

    parser.add_argument("--blanket", action="store_true")

    parser.add_argument("--generic", action="store_true")

    args = parser.parse_args()

    main(args)
    print(BColors.ENDC)
