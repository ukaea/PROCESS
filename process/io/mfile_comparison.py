#!/usr/bin/env python
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

import sys
import numpy
import argparse
import process.io.mfile as mf
from numpy import isfinite
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
    "sarea",
    "vol",
    "n_tf",
    "shldith",
    "shldoth",
    "blnkith",
    "blnkoth",
    "powfmw",
    "plascur/1d6",
    "bt",
    "q95",
    "betap",
    "te",
    "dene",
    "hfact",
    "vstot",
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
    "ralpne",
    "pinnerzoneradmw",
    "pradmw",
    "pnucblkt",
    "pnucshld",
    "pdivt",
    "pheat",
    "bootipf",
    "faccd",
    "facoh",
    "gamnb",
    "enbeam",
    "powerht",
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
    "sarea",
    "vol",
    "n_tf",
    "powfmw",
    "plascur/1d6",
    "bt",
    "q95",
    "beta",
    "normalised_thermal_beta",
    "normalised_total_beta",
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
    "dnz",
    "taueff",
    "hfact",
    "tauelaw",
    "ralpne",
    "wallmw",
    "pinnerzoneradmw",
    "psyncpv*vol",
    "pradmw",
    "pnucblkt",
    "pnucshld",
    "pdivt",
    "divlife",
    "pthermmw",
    "bore",
    "ohcth",
    "precomp",
    "gapoh",
    "tfcth",
    "deltf",
    "thshield_ib",
    "thshield_ob",
    "thshield_vb",
    "gapds",
    "d_vv_in",
    "d_vv_out",
    "d_vv_top",
    "d_vv_bot",
    "shldith",
    "vvblgap",
    "blnkith",
    "fwith",
    "scrapli",
    "scraplo",
    "fwoth",
    "blnkoth",
    "shldoth",
    "gapsto",
    "tftsgap",
    "tfthko",
    "etath",
    "pgrossmw",
    "pnetelmw",
    "pinjmw",
    "pheat",
    "bootipf",
    "faccd",
    "facoh",
    "gamnb",
    "enbeam",
    "powerht",
    "pdivt",
    "vssoft",
    "vstot",
    "tburn",
    "bmaxtf",
    "iooic",
    "tmarg",
    "tftmp",
    "qtarget",
    "qtargetcomplete",
    "totalpowerlost",
]

BLANKET_COMPARE_PARAMS = [
    "blnkith",
    "blnkoth",
    "powfmw",
    "pnucblkt",
    "pnucfw",
    "ptfnuc",
    "pnucshld",
    "pnucdiv",
    "tbr",
    "li6enrich",
    "fwarea",
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
    "powfmw",
    "plascur/1d6",
    "bt",
    "q95",
    "beta",
    "te",
    "dene",
    "pinjmw",
    "pnetelmw",
    "wallmw",
    "photon_wall",
    "ralpne",
    "pinnerzoneradmw",
    "pradmw",
    "bootipf",
    "pdivmax/rmajor",
    "fimp(14",
    "etath",
    "capcost",
    "coe",
    "bore",
    "ohcth",
    "precomp",
    "gapoh",
    "tfcth",
    "tftsgap",
    "thshield_ib",
    "thshield_ob",
    "thshield_vb",
    "gapds",
    "d_vv_in",
    "shldith",
    "vvblgap",
    "blnkith",
    "fwith",
    "scrapli",
    "scraplo",
    "fwoth",
    "blnkoth",
    "shldoth",
    "d_vv_out",
    "gapsto",
    "tftsgap",
    "tfthko",
    "vgap",
    "divfix",
    "d_vv_bot",
    "shldlth",
    "vgap2",
]


class BColors(object):
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
    mfile_list = list()
    for item in arg.f:
        mfile = mf.MFile(filename=item)
        if mfile.data["error_status"].get_scan(-1) == 3:
            raise RuntimeError(
                f"{item} is an MFile from a PROCESS run that did not converge"
                " and instead results from an error during the run"
            )

        mfile_list.append(mfile)

    var_list = list()
    missing_vars = list()
    diff_list = list()
    within_list = list()

    key_list = mfile_list[0].data.keys()
    for var in key_list:
        store = True
        for f in mfile_list:
            if var not in f.data.keys():
                store = False
                missing_vars.append(var)

        if store:
            var_list.append(var)

    if arg.save:
        ofile = open("comp.txt", "w")

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

        values = numpy.zeros(n)  # replaced scipy with numpy

        if v not in get_dicts()["DICT_VAR_TYPE"].keys():
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

        norm_vals = list()
        if values[0] != 0 and isfinite(values[0]):
            norm_vals = values / values[0]
        # else:
        #    print(key, values[0])

        if len(norm_vals) >= 1:
            key = v.strip(".").strip(" ")
            if key in get_dicts()["DICT_DESCRIPTIONS"].keys():
                des = get_dicts()["DICT_DESCRIPTIONS"][key]
            else:
                des = "-"
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
                    ofile.write(wline + "\n")

    if arg.save:
        ofile.close()

    if arg.baseline:
        if arg.acc >= 10.0:
            if print_counter == 0:
                sys.exit(0)
            else:
                sys.exit(
                    "Differences in baseline output by more than {0}%".format(arg.acc)
                )


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
