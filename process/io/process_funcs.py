"""
A selection of functions for using the PROCESS code

Author: Hanni Lux (Hanni.Lux@ccfe.ac.uk)

Compatible with PROCESS version 368

24/11/2021: Global dictionary variables moved within the functions
            to avoid cyclic dependencies. This is because the dicts
            generation script imports, and inspects, process.
"""

import logging
from os.path import join as pjoin
from pathlib import Path
from sys import stderr
from time import sleep

from process.data_structure import numerics
from process.io.data_structure_dicts import get_dicts
from process.io.in_dat import InDat
from process.io.mfile import MFile

logger = logging.getLogger(__name__)


def get_neqns_itervars(wdir="."):
    """
    returns the number of equations and a list of variable
    names of all iteration variables
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()
    in_dat = InDat(pjoin(wdir, "IN.DAT"))

    ixc_list = in_dat.data["ixc"].get_value

    itervars = []
    for var in ixc_list:
        if var != "":
            itervars += [dicts["DICT_IXC_SIMPLE"][str(var)]]

    assert in_dat.number_of_itvars == len(itervars)

    return in_dat.number_of_constraints, itervars


###############################


def update_ixc_bounds(wdir="."):
    """
    updates the lower and upper bounds in DICT_IXC_BOUNDS
    from IN.DAT
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()
    in_dat = InDat(pjoin(wdir, "IN.DAT"))

    bounds = in_dat.data["bounds"].get_value

    for key, value in bounds.items():
        name = dicts["DICT_IXC_SIMPLE"][key]

        if "l" in value:
            dicts["DICT_IXC_BOUNDS"][name]["lb"] = float(value["l"])
        if "u" in value:
            dicts["DICT_IXC_BOUNDS"][name]["ub"] = float(value["u"])


###############################


def get_variable_range(itervars, factor, wdir="."):
    """
    Returns the lower and upper bounds of the variable range
    for each iteration variable.

    itervars - string list of all iteration variable names
    factor   - defines the variation range for non-f-values by
               setting them to value * factor and value / factor
               respectively while taking their process bounds
               into account.

    For f-values the allowed range is equal to their process bounds.

    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()
    in_dat = InDat(pjoin(wdir, "IN.DAT"))

    lbs = []
    ubs = []

    iteration_variables = numerics.lablxc

    for varname in itervars:
        iteration_variable_index = iteration_variables.index(varname)
        lb = numerics.boundl[iteration_variable_index]
        ub = numerics.boundu[iteration_variable_index]
        # for f-values we set the same range as in process
        if varname[0] == "f" and (varname not in dicts["NON_F_VALUES"]):
            lbs += [lb]
            ubs += [ub]

        # for non-f-values we modify the range with the factor
        else:
            value = get_from_indat_or_default(in_dat, varname)

            if value is None:
                print(f"Error: Iteration variable {varname} has None value!")
                exit()

            # to allow the factor to have some influence
            if value == 0.0:
                value = 1.0

            # assure value is within bounds!
            if value < lb:
                value = lb
            elif value > ub:
                value = ub

            if value > 0:
                lbs += [max(value / factor, lb)]
                ubs += [min(value * factor, ub)]
            else:
                lbs += [min(value / factor, lb)]
                ubs += [max(value * factor, ub)]

        if lbs[-1] > ubs[-1]:
            print(
                f"Error: Iteration variable {varname} has BOUNDL={lbs[-1]} >"
                f"BOUNDU={ubs[-1]}\n Update process_dicts or input file!",
                file=stderr,
            )

            exit()
        # assert lbs[-1] < ubs[-1]

    return lbs, ubs


###############################


def check_in_dat():
    """Tests IN.DAT during setup:
    1)Are ixc bounds outside of allowed input ranges?"""
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    in_dat = InDat()

    # Necessary for the correct use of this function as well as
    # get_variable_range
    update_ixc_bounds()

    # 1) Are ixc bounds outside of allowed input ranges?

    ixc_list = in_dat.data["ixc"].get_value

    for itervarno in ixc_list:
        itervarname = dicts["DICT_IXC_SIMPLE"][str(itervarno)]
        try:
            lowerinputbound = dicts["DICT_INPUT_BOUNDS"][itervarname]["lb"]
        except KeyError as err:
            # arrays do not have input bound checks
            if "(" in itervarname:
                continue

            print("Error in check_in_dat():")
            print("There seems to be some information missing from the dicts.")
            print("Please flag this up for a developer to investigate!")
            print(itervarname, err)
            print(dicts["DICT_INPUT_BOUNDS"][itervarname])
            exit()

        if dicts["DICT_IXC_BOUNDS"][itervarname]["lb"] < lowerinputbound:
            print(
                "Warning: boundl for ",
                itervarname,
                " lies out of allowed input range!\n Reset boundl(",
                itervarno,
                ") to ",
                lowerinputbound,
                file=stderr,
            )
            dicts["DICT_IXC_BOUNDS"][itervarname]["lb"] = lowerinputbound
            set_variable_in_indat(
                in_dat, "boundl(" + str(itervarno) + ")", lowerinputbound
            )
            sleep(1)

        upperinputbound = dicts["DICT_INPUT_BOUNDS"][itervarname]["ub"]

        if dicts["DICT_IXC_BOUNDS"][itervarname]["ub"] > upperinputbound:
            print(
                "Warning: boundu for",
                itervarname,
                f"lies out of allowed input range!\n Reset boundu({itervarno}) to",
                upperinputbound,
                file=stderr,
            )
            dicts["DICT_IXC_BOUNDS"][itervarname]["ub"] = upperinputbound
            set_variable_in_indat(
                in_dat, "boundu(" + str(itervarno) + ")", upperinputbound
            )
            sleep(1)

    in_dat.write_in_dat(output_filename="IN.DAT")


###############################


def check_logfile(logfile="process.log"):
    """
    Checks the log file of the PROCESS output.
    Stops, if an error occured that needs to be
    fixed before rerunning.
    XXX should be deprecated!! and replaced by check_input_error!
    """

    with open(logfile) as outlogfile:
        errormessage = "Please check the output file for further information."
        for line in outlogfile:
            if errormessage in line:
                print(
                    "An Error has occured. Please check the output file for more information.",
                    file=stderr,
                )
                exit()


def check_input_error(wdir="."):
    """
    Checks, if an input error has occurred.
    Stops as a consequence.
    Will also fail if the MFILE.DAT isn't found.
    """
    try:
        mfile_path = Path(wdir) / "MFILE.DAT"

        if mfile_path.exists():
            mfile_path_str = str(mfile_path)
            mfile = MFile(filename=mfile_path_str)
        else:
            raise FileNotFoundError("MFile doesn't exist")

        error_id = mfile.data["error_id"].get_scan(-1)

        if error_id == 130:
            print(
                "Error in input file. Please check OUT.DAT for more information.",
                file=stderr,
            )
            exit()
    except Exception:
        logger.exception("Check input error exception")
        raise


########################################


def process_stopped(wdir="."):
    """
    Checks the process Mfile whether it has
    prematurely stopped.
    """
    # Check for MFILE
    try:
        m_file = MFile(filename=pjoin(wdir, "MFILE.DAT"))
    except FileNotFoundError as err:
        print(f"No MFILE has been found! FYI:qn {err}", file=stderr)
        print("Code continues to run!", file=stderr)
        return True

    # Get error status; missing key indicates premature exit of Process
    # (usually a STOP 1)
    try:
        error_status = m_file.data["error_status"].get_scan(-1)
    except KeyError:
        return True

    # Process did prematurely exit
    return error_status >= 3


########################################


def process_warnings(wdir="."):
    """
    Checks the process Mfile whether any
    warnings have occurred.
    """

    m_file = MFile(filename=pjoin(wdir, "MFILE.DAT"))
    error_status = m_file.data["error_status"].get_scan(-1)

    return error_status >= 2


############################################


def mfile_exists():
    """checks whether MFILE.DAT exists"""

    try:
        with open("MFILE.DAT") as m_file:
            m_file.close()
        return True

    except FileNotFoundError:
        return False


############################################


def no_unfeasible_mfile(wdir="."):
    """
    returns the number of unfeasible points
    in a scan in MFILE.DAT
    """

    m_file = MFile(filename=pjoin(wdir, "MFILE.DAT"))

    # no scans
    if not m_file.data["isweep"].exists:
        if m_file.data["ifail"].get_scan(-1) == 1:
            return 0
        return 1

    ifail = m_file.data["ifail"].get_scans()
    try:
        return len(ifail) - ifail.count(1)
    except TypeError:
        # This seems to occur, if ifail is not in MFILE!
        # This probably means in the mfile library a KeyError
        # should be raised not only a message to stdout!
        return 100000


################################


def vary_iteration_variables(itervars, lbs, ubs, generator):
    """
    Routine to change the iteration variables in IN.DAT
    within given bounds.
    itervars  - string list of all iteration variable names
    lbs       - float list of lower bounds for variables
    ubs       - float list of upper bounds for variables
    generator - Generator numpy generator to create random numbers
    """

    in_dat = InDat()

    new_values = []

    for varname, lbnd, ubnd in zip(itervars, lbs, ubs, strict=False):
        new_value = generator.uniform(lbnd, ubnd)
        new_values += [new_value]
        in_dat.add_parameter(varname, new_value)

    in_dat.write_in_dat(output_filename="IN.DAT")

    return new_values


###################################


def get_solution_from_mfile(neqns, nvars, wdir="."):
    """
    returns
    ifail - error_value of VMCON/PROCESS
    the objective functions
    the square root of the sum of the squares of the constraints
    a list of the final iteration variable values
    a list of the final constraint residue values

    If the run was a scan, the values of the last scan point
    will be returned.
    """
    m_file = MFile(filename=pjoin(wdir, "MFILE.DAT"))

    ifail = m_file.data["ifail"].get_scan(-1)

    # figure of merit objective function
    objective_function = m_file.data["f"].get_scan(-1)

    # estimate of the constraints
    constraints = m_file.data["sqsumsq"].get_scan(-1)

    table_sol = [
        m_file.data[f"itvar{var_no + 1:03}"].get_scan(-1) for var_no in range(nvars)
    ]
    table_res = [
        m_file.data[f"normres{con_no + 1:03}"].get_scan(-1) for con_no in range(neqns)
    ]

    if ifail != 1:
        return ifail, "0", "0", ["0"] * nvars, ["0"] * neqns

    return ifail, objective_function, constraints, table_sol, table_res


############################################


def get_from_indat_or_default(in_dat, varname):
    """quick function to get variable value from IN.DAT
    or PROCESS default value"""
    dicts = get_dicts()
    if varname in in_dat.data:
        return in_dat.data[varname].get_value
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    return dicts["DICT_DEFAULT"][varname]


def set_variable_in_indat(in_dat, varname, value):
    """quick function that sets a variable value in
    IN.DAT and creates it if necessary"""

    varname = varname.lower()
    if "bound" in varname:
        number = (varname.split("("))[1].split(")")[0]
        if "boundu" in varname:
            in_dat.add_bound(number, "u", value)
        else:
            in_dat.add_bound(number, "l", value)

    elif "(" in varname:
        name = varname.split("(")[0]
        # Fortran numbering converted to Python numbering!
        identity = int((varname.split("("))[1].split(")")[0]) - 1
        in_dat.change_array(name, identity, float(value))

    else:
        in_dat.add_parameter(varname, value)
