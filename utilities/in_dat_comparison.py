#!/usr/bin/env python
"""
  PROCESS IN.DAT comparison tool

  Tool to compare two IN.DAT files and list which variable are in one and not
  another and also compare parameters with different values. Finally list parameters
  that are in both and the same

  James Morris 15/11/2016
  CCFE
"""

# Third party libraries
import sys
import argparse

# PROCESS libraries
from process.io.in_dat import InDat

# Constants

# *********************************** #


def compare_input_files(f1, f2):
    """Compare

    :param f1: file 1 data object
    :param f2: file 2 data object
    """

    file_1_keys = f1.data.keys()
    file_2_keys = f2.data.keys()

    both = set(file_1_keys).intersection(file_2_keys)
    both_same = [key for key in both if f1.data[key].value == f2.data[key].value]
    both_diff = set(both).difference(both_same)

    one_not_two = set(file_1_keys).difference(file_2_keys)
    two_not_one = set(file_2_keys).difference(file_1_keys)

    return both, both_same, both_diff, one_not_two, two_not_one


def output(kb, kbs, kbd, k1n2, k2n1, f1, f2, c_args):
    """Output the differences between to the two input files to either
    terminal or input_comp.txt

    :param kb: keys in both input files
    :param kbs: keys in both input files with same value
    :param kbd: keys in both input files with different value
    :param k1n2: keys in input file 1 but not file 2
    :param k2n1: keys in input file 2 but not file 1
    :param f1: input file 1
    :param f2: input file 2
    :param c_args: command arguments
    """

    if c_args.save:
        location = open("input_comp.txt", "w")
    else:
        location = sys.stdout

    print("\nPROCESS Input file comparison tool\n", file=location)

    print(
        "\nThere are {0} variables in {1} NOT IN {2}:\n".format(
            len(k1n2), c_args.f[0], c_args.f[1]
        ),
        file=location,
    )
    print("\n".join("{0} : {1}".format(*k) for k in enumerate(k1n2)), file=location)

    print(
        "\nThere are {0} variables in {1} NOT IN {2}:\n".format(
            len(k2n1), c_args.f[1], c_args.f[0]
        ),
        file=location,
    )
    print("\n".join("{0} : {1}".format(*k) for k in enumerate(k2n1)), file=location)

    print(
        "\nThere are {0} variables in BOTH {1} AND {2} with the DIFFERENT values:\n".format(
            len(kbd), c_args.f[0], c_args.f[1]
        ),
        file=location,
    )
    print(
        "\n".join(
            "{0} : {1} \n\t {2} = {3} \n\t {4} = {5}\n".format(
                k[0],
                k[1],
                c_args.f[0],
                f1.data[k[1]].value,
                c_args.f[1],
                f2.data[k[1]].value,
            )
            for k in enumerate(kbd)
        ),
        file=location,
    )

    print(
        "\nThere are {0} variables in BOTH {1} AND {2} with the SAME values:\n".format(
            len(kbs), c_args.f[0], c_args.f[1]
        ),
        file=location,
    )
    print(
        "\n".join(
            "{0} : {1} \n\t {2} = {3} \n\t {4} = {5}\n".format(
                k[0],
                k[1],
                c_args.f[0],
                f1.data[k[1]].value,
                c_args.f[1],
                f2.data[k[1]].value,
            )
            for k in enumerate(kbs)
        ),
        file=location,
    )


def main(cmd_args):
    """Main

    :param cmd_args: command line arguments
    """

    # Get input file names
    file_names = cmd_args.f

    # Read input files
    file_1 = InDat(filename=file_names[0])
    file_2 = InDat(filename=file_names[1])

    keys_both, keys_both_same, keys_both_diff, keys_1n2, keys_2n1 = compare_input_files(
        file_1, file_2
    )

    output(
        keys_both,
        keys_both_same,
        keys_both_diff,
        keys_1n2,
        keys_2n1,
        file_1,
        file_2,
        cmd_args,
    )

    return


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(
        description="Produce a comparison "
        "between two PROCESS "
        "input files."
        "For info contact "
        "james.morris2@ukaea.uk"
    )

    PARSER.add_argument("-f", metavar="f", type=str, nargs="+", help="Files to compare")

    PARSER.add_argument(
        "-s",
        "--save",
        help="Save output to file called input_comp.txt",
        action="store_true",
    )

    ARGS = PARSER.parse_args()

    main(ARGS)
