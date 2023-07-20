"""
  Reference check comparison script

  Takes reference from JSON file

  Script compares output MFILE with the required parameters

  James Morris 24/11/2015
  CCFE
"""

# Third party libraries
import sys
import json
import argparse
import collections

# PROCESS libraries
from process.io.mfile import MFile

# Constants

OK = "\033[92m"
FAIL = "\033[91m"
ENDC = "\033[0m"

# *********************************** #


def param_comparison(cargs, m_file, loc):
    """Parameter comparison for CFETR small tool

    :param cargs: commands line arguments
    :param m_file: MFile object
    :param loc: output location
    """

    mfile_keys = m_file.data.keys()

    param_keys = PARAMS.keys()

    params_not_in_mfile = set(param_keys).difference(mfile_keys)

    params_in_mfile = set(param_keys).difference(params_not_in_mfile)

    same = [
        key for key in params_in_mfile if PARAMS[key] == m_file.data[key].get_scan(-1)
    ]
    diff = set(params_in_mfile).difference(same)

    print("\nParameters", file=loc)

    print(
        FAIL + "".join("{0}\n".format(k) for k in params_not_in_mfile) + ENDC, file=loc
    )

    print(
        FAIL
        + "".join(
            "{0:<15} \t {1:<5} = {2:<10.3g} \t {3:<5} = {4:.3g}\n".format(
                k, cargs.r, PARAMS[k], cargs.f, m_file.data[k].get_scan(-1)
            )
            for k in diff
        )
        + ENDC,
        file=loc,
    )

    print(
        OK
        + "".join(
            "{0:<15} \t {1:<5} = {2:<10.3g} \t {3:<5} = {4:.3g}\n".format(
                k, cargs.r, PARAMS[k], cargs.f, m_file.data[k].get_scan(-1)
            )
            for k in same
        )
        + ENDC,
        file=loc,
    )


def limit_comparison(cargs, m_file, loc):
    """Limit comparison for reference check tool

    :param cargs: commands line arguments
    :param m_file: MFile object
    :param loc: output location
    """

    LIMIT_KEYS = LIMITS.keys()

    EXCEEDED = [
        key
        for key in LIMIT_KEYS
        if LIMITS[key].type == "-"
        and LIMITS[key].limit < m_file.data[LIMITS[key].parameter].get_scan(-1)
    ]

    NOT_EXCEEDED = set([k for k in LIMIT_KEYS if LIMITS[k].type == "-"]).difference(
        EXCEEDED
    )

    NOT_MET = [
        key
        for key in LIMIT_KEYS
        if LIMITS[key].type == "+"
        and LIMITS[key].limit > m_file.data[LIMITS[key].parameter].get_scan(-1)
    ]

    MET = set([k for k in LIMIT_KEYS if LIMITS[k].type == "+"]).difference(NOT_MET)

    print("\nLimits - Upper", file=loc)

    print(
        FAIL
        + "".join(
            "{0:<10}\t{1:<20}\t{2:<5}={3:<5.3g}\t{4:<5}={5:.3g}\n".format(
                k,
                LIMITS[k].parameter,
                cargs.r,
                LIMITS[k].limit,
                cargs.f,
                m_file.data[LIMITS[k].parameter].get_scan(-1),
            )
            for k in EXCEEDED
        )
        + ENDC,
        file=loc,
    )

    print(
        OK
        + "".join(
            "{0:<10}\t{1:<20}\t{2:<5}={3:<5.3g}\t {4:<5}={5:.3g}\n".format(
                k,
                LIMITS[k].parameter,
                cargs.r,
                LIMITS[k].limit,
                cargs.f,
                m_file.data[LIMITS[k].parameter].get_scan(-1),
            )
            for k in NOT_EXCEEDED
        )
        + ENDC,
        file=loc,
    )

    print("\nLimits - Lower", file=loc)

    print(
        FAIL
        + "".join(
            "{0:<10}\t{1:<20}\t{2:<5}={3:<5.3g}\t{4:<5}={5:.3g}\n".format(
                k,
                LIMITS[k].parameter,
                cargs.r,
                LIMITS[k].limit,
                cargs.f,
                m_file.data[LIMITS[k].parameter].get_scan(-1),
            )
            for k in NOT_MET
        )
        + ENDC,
        file=loc,
    )

    print(
        OK
        + "".join(
            "{0:<10}\t{1:<20}\t{2:<5}={3:<5.3g}\t{4:<5}={5:.3g}\n".format(
                k,
                LIMITS[k].parameter,
                cargs.r,
                LIMITS[k].limit,
                cargs.f,
                m_file.data[LIMITS[k].parameter].get_scan(-1),
            )
            for k in MET
        )
        + ENDC,
        file=loc,
    )


def main(cmd_args):
    """Main function for CFETR small comparison

    :param cmd_args: command line arguments
    """

    if cmd_args.save:
        location = open("cfetr_comp.txt", "w")
    else:
        location = sys.stdout

    print("\nCFETR comparison with output file '{0}':\n".format(cmd_args.f))

    # Get input file names
    file_name = cmd_args.f

    # Read MFILE files
    mfile = MFile(filename=file_name)

    # Compare MFILE to PARAMS
    param_comparison(cmd_args, mfile, location)

    # Compare MFILE to LIMITS
    limit_comparison(cmd_args, mfile, location)


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        description="Compare a PROCESS output MFILE and a reference reference. "
        "For info contact james.morris2@ukaea.uk"
    )

    PARSER.add_argument(
        "-r",
        metavar="REFNAME",
        type=str,
        default="",
        help="specify reference JSON file",
    )

    PARSER.add_argument(
        "-f", metavar="MFILENAME", type=str, default="", help="specify PROCESS MFILE"
    )

    PARSER.add_argument(
        "-s",
        "--save",
        help="Save output to file called" "cfetr_comp.txt",
        action="store_true",
    )

    COMMAND_ARGS = PARSER.parse_args()

    # read json
    if COMMAND_ARGS.r:
        print("Reading reference file {0}".format(COMMAND_ARGS.r))
        REF = json.load(open(COMMAND_ARGS.r), object_pairs_hook=collections.OrderedDict)
        PARAMS = REF["PARAMS"]
        LIMITS_NO_TUPLES = REF["LIMITS"]

        LIMITS = {
            key: collections.namedtuple("GenericDict", LIMITS_NO_TUPLES[key].keys())(
                **LIMITS_NO_TUPLES[key]
            )
            for key in list(LIMITS_NO_TUPLES)
        }
    else:
        print("Please enter a reference JSON file with -r argument!")

    if COMMAND_ARGS.f:
        main(COMMAND_ARGS)
    else:
        print("Please enter a reference MFILE with -f argument!")

    # Make sure terminal returns to regular colours
    print(ENDC)
