"""

Modifies the PROCESS input file IN.DAT so all the iteration variables are
given their values from the output file MFILE.DAT.

"""

import re

import process.core.io.mfile.mfile as mf
from process.core.io.in_dat.base import InDat


def feasible_point(filename, position: int):
    """Function to check for feasible solution before creating new IN.DAT, or to determine the first or last feasible point in a scan

    Parameters
    ----------
    filename :
        name of MFILE.DAT to read
    position :
        feasible position index

    Returns
    -------
    scan_point:
        scan number to use when writing new file
    """
    mfile_data = mf.MFile(filename)
    scan_point = 0
    num_scans = int(mfile_data.get("isweep", scan=-1) or 1)

    if position == -1:
        position = num_scans

    position = max(1, position)

    if position > num_scans:
        position = num_scans
        print(f"Only {num_scans} in mfile selecting last feasible_point")

    check_point = 1

    for value in mfile_data.data:
        # Look for feasible scan points (with ifail = 1)
        if "ifail" in value and "vmcon_error_flag_(ifail)" not in value:
            if mfile_data.get(value, scan=check_point) == 1:
                scan_point = check_point
            if check_point == position:
                break
            check_point += 1
    else:
        raise ValueError("No feasible point found")
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


def write_indat(mfile, indat, output, feasible_point_index):
    scan = feasible_point(mfile, feasible_point_index)
    print("Using scan number = ", scan)

    # Get iteration variables from MFILE.DAT
    it_vars = get_iteration_variables(mfile, scan)

    # Read IN.DAT
    in_dat_data = InDat(indat)

    # Amend the values for the iteration variables
    in_dat_data = replace_iteration_variables(it_vars, in_dat_data)

    # Write a new IN.DAT
    in_dat_data.write_in_dat(output_filename=output)
