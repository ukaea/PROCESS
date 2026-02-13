"""

Modifies the PROCESS input file IN.DAT so all the iteration variables are
given their values from the output file MFILE.DAT.

James Morris 30/04/2014 based on code by Michael Kovari 9/8/13 and
J C Rivas, 16/7/2013
"""

import argparse
import re

import process.io.mfile as mf
from process.io.in_dat import InDat


def feasible_point(filename, position):
    """Function to check for feasible solution before creating new IN.DAT, or to determine the first or last feasible point in a scan

    Parameters
    ----------
    filename :
        name of MFILE.DAT to read
    position :
        e.g first or last

    Returns
    -------
    scanPoint:
        scan number to use when writing new file
    """
    mfile_data = mf.MFile(filename)
    finished = False
    scan_point = 0
    num_scans = 1

    for value in mfile_data.data:
        if "isweep" in value:
            num_scans = int(mfile_data.data["isweep"].get_scan(-1))
            break

    # Assign start point
    check_point = 1 if position == "first" else num_scans

    while finished is False:
        for value in mfile_data.data:
            # Look for feasible scan points (with ifail = 1)
            if "ifail" in value and "vmcon_error_flag_(ifail)" not in value:
                if mfile_data.data[value].get_scan(check_point) == 1:
                    finished = True
                    scan_point = check_point
                else:
                    if position == "last":
                        if check_point == 1:
                            finished = True
                        else:
                            check_point = check_point - 1
                    elif position == "first":
                        if check_point == num_scans:
                            finished = True
                        else:
                            check_point = check_point + 1
    return scan_point


def get_iteration_variables(filename, scan):
    """Function to get a list of the iteration variables and their values from
    MFILE.DAT

    Parameters
    ----------
    filename:
        name of MFILE.DAT to read
    scan:
        scan number to use

    Returns
    -------
    iteration_vars:
        dictionary of iteration variables in MFILE and their
        values.
    """
    mfile_data = mf.MFile(filename)
    iteration_vars = {}

    for value in mfile_data.data:
        if "itvar" in value and "nitvar" not in value:
            variable_name = mfile_data.data[value].var_description
            variable_value = mfile_data.data[value].get_scan(scan)
            iteration_vars[variable_name] = variable_value

    return iteration_vars


def replace_iteration_variables(iteration_vars, in_data):
    """Function to replace the iteration variables in IN.DAT if the variable
    is not defined in IN.DAT it will add the variable to the end of the file.

    Parameters
    ----------
    iteration_vars:
        dictionary of iteration variables from MFILE.DAT and
        their values
    in_data:
        IN.DAT data object.
    """

    for variable_name, variable_value in iteration_vars.items():
        if (match := re.search(r"([a-zA-Z0-9_]+)\(([0-9]+)\)", variable_name)) is None:
            in_data.add_parameter(variable_name.lower(), variable_value)
        else:
            in_data.change_array(match.group(1), int(match.group(2)) - 1, variable_value)

    return in_data


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Creates a new IN.DAT using iteration variable values from MFILE.DAT."
    )

    parser.add_argument(
        "-f",
        metavar="f",
        type=str,
        default="MFILE.DAT",
        help='File to read as MFILE.DAT (default="MFILE.DAT")',
    )

    parser.add_argument(
        "-i",
        metavar="i",
        type=str,
        default="IN.DAT",
        help='File to read as IN.DAT (default="IN.DAT")',
    )

    parser.add_argument(
        "-o",
        metavar="o",
        type=str,
        default="new_IN.DAT",
        help='File to write as new IN.DAT (default="new_IN.DAT")',
    )

    parser.add_argument(
        "-lfp",
        "--lfp",
        help="use last feasible point in a scan (default)",
        action="store_true",
    )

    parser.add_argument(
        "-ffp", "--ffp", help="use first feasible point in a scan", action="store_true"
    )

    args = parser.parse_args(args)

    if args.ffp:
        # Determine first feasible scan point
        scan = feasible_point(args.f, "first")
    else:
        # Determine last feasible scan point - default
        scan = feasible_point(args.f, "last")
    if scan == 0:
        print("No feasible points in scan")
        raise SystemExit
    print("Using scan number = ", scan)

    # Get iteration variables from MFILE.DAT
    it_vars = get_iteration_variables(args.f, scan)

    # Read IN.DAT
    in_dat_data = InDat(args.i)

    # Amend the values for the iteration variables
    in_dat_data = replace_iteration_variables(it_vars, in_dat_data)

    # Write a new IN.DAT
    in_dat_data.write_in_dat(output_filename=args.o)


if __name__ == "__main__":
    main()
