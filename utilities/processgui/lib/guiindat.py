#!/usr/bin/env python3

"""Contains a class that is used to read and write IN.DAT files.
   files. Can be used to convert an old IN.DAT file to dictionary
   or vice versa.
   When this file is run as a program, converts IN.DAT file to new
   format.

   Tom Miller 09/14

"""

from process_io_lib.process_dicts import (
    DICT_DEFAULT,
    DICT_MODULE,
    DICT_DESCRIPTIONS,
    DICTIONARY_VERSION,
    DICT_IXC_SIMPLE,
)
import copy
import argparse

# variables in input file but not in dictionaries.
OMISSIONS = list()

# ioptimz values
ioptimz_des = {
    "-1": "for no optimisation HYBRD only",
    "0": "for HYBRD and VMCON (not recommended)",
    "1": "for optimisation VMCON only",
}


class BColours:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def is_valid(varname):
    """Tests to see if a varname is a valid variable to have in the IN.DAT
    by seeing if it appears in DICT_DEFAULT

    """
    if varname not in DICT_DEFAULT:
        error_msg = (
            BColours.FAIL + "Unrecognised input variable: {0}. "
            "Variable will be skipped".format(varname) + BColours.ENDC
        )
        print(error_msg)
        OMISSIONS.append(varname)
        # raise ValueError(error_msg)
    return True


def get_type(varname):
    """Gets the type of a varname. Undestands arrays ie. ixc(2) is int,
    dcond(3) is float

    """
    if "(" in varname:
        array_name = varname.split("(")[0]
        is_valid(array_name)
        return type(DICT_DEFAULT[array_name][0])
    else:
        is_valid(varname)
        return type(DICT_DEFAULT[varname])


def mimic_type(in_val, tar):
    """Converts in_val to have the same type as tar"""
    if isinstance(in_val, str):
        in_val = in_val.strip(", ")
    if isinstance(tar, list):
        # allow trailing commas for lists
        in_val = str(in_val)
        in_val = in_val.strip("[]")
        ret = []
        for i in in_val.split(","):
            # cast every variable in the list
            ret.append(mimic_type(i, tar[0]))
        return ret
    elif isinstance(tar, int):
        return int(in_val)
    elif isinstance(tar, float):
        return float(str(in_val).lower().replace("d", "e"))
    elif isinstance(tar, str):
        # remove all quote marks and spaces
        return str(in_val).strip("\"' ")
    else:
        raise ValueError("Unknown type for variable passed to mimic_type")


def cast(varname, val):
    """Casts a string input from user 'val' to a python variable of
    the same type as DICT_DEFAULT[varname]

    """
    try:
        if "(" in varname:
            array_name = varname.split("(")[0]
            return mimic_type(val, DICT_DEFAULT[array_name][0])
        else:
            if varname not in OMISSIONS:
                return mimic_type(val, DICT_DEFAULT[varname])
    except ValueError:
        error_msg = (
            BColours.FAIL + "Could not cast {0} to type {1} for"
            "variable {2}. Please check DICT_DEFAULT for parameter type.".format(
                str(val), str(type(DICT_DEFAULT[varname])), varname
            )
            + BColours.ENDC
        )
        print(error_msg)
        # raise ValueError("Could not cast " + str(val) + " to type " +
        #                  str(type(DICT_DEFAULT[varname])) + " for variable "+
        #                  varname)


def interpret_array_var(st):
    """Interprets a fortran style array access eg boundl(4). Returns
    the array name and python style index - eg boundl, 3
    """
    st = st.split(")")[0]
    array_name, array_st_index = st.split("(")
    is_valid(array_name)
    return array_name, int(array_st_index) - 1


def parse_int_line(line, char):
    """Helper function for get_comment. Splits the line by char,
    attempts to convert the first chunk to an int, returns the number
    and rest of the line if succesful, None otherwise
    """
    assert len(char) == 1
    li = line.split(char)
    try:
        num = int(li[0])
    except ValueError:
        return None
    return num, char.join(li[1:]).strip()


def get_line_by_int(num, desc):
    """Looks through the line in desc for a line that looks like it's
    describing the value of the variable. See if the line starts with
    = #, (#), # where # represents an integer. If this integer = num or
    -num, return the rest of the line.
    """
    for line in desc.split("\n"):
        line = line.strip()
        # look for '= #' type lines
        if line[0] == "=":
            line = line[1:].strip()
            t = parse_int_line(line, " ")
            if t and (t[0] == num or t[0] == -num):
                return t[1].strip("* ")
        # look for (#) style lines
        if line[0] == "(":
            line = line[1:].strip()
            t = parse_int_line(line, ")")
            if t and (t[0] == num or t[0] == -num):
                return t[1].strip("* ")
        # look for lines that start with an int
        if line[0].isdigit():
            t = parse_int_line(line, " ")
            if t and (t[0] == num or t[0] == -num):
                return t[1].strip("* ")
    return ""


def get_comment(val, desc):
    """Gets the description that should be printed in the IN.DAT for
    #a variable. For non-integer variables, this is the first line
    #of GUI_DESCRIPTION[var]. For integer switches, try to work out
    #which switch is selected
    """

    if desc.strip() == "":
        return ""
    # commas, full stops, colons  must not appear in comments
    desc = desc.replace(",", ";")
    desc = desc.replace(".", ";")
    desc = desc.replace(":", "*")
    firstline = desc.split("\n")[0].strip("* ")
    if not isinstance(val, int):
        # for non ints, just use the first line of the desc
        return "* " + firstline

    # for ints, get rid of parenthesis
    ret = firstline.split("(")[0].strip("* ")

    # then try and describe the value taken
    int_desc = get_line_by_int(val, desc)
    if int_desc:
        return "* " + ret + " * " + int_desc
    else:
        return "* " + ret


def make_mod_header(st):
    """Makes the header line for a module"""
    return "*" + st.center(50, "-") + "*"


def make_line(var, val, comment=True):
    """Returns a var = val * comment line
    if comment = False, don't try and find a comment

    """
    assignstr = (str(var) + " = " + str(val)).ljust(15)
    if comment and var in DICT_DESCRIPTIONS:
        if var == "ioptimz":
            desc = ioptimz_des[str(val)]
        else:
            desc = get_comment(val, DICT_DESCRIPTIONS[var])
    else:
        desc = ""

    return assignstr + "  " + desc


def strip_line(line):
    """Strips of irrelevant parts of a line such as comments, old-style
    module headings

    """
    line = line.split("*")[0]
    line = line.strip()
    line = line.rstrip("\r")
    line = line.rstrip(",")
    if line == "" or line[0] == "$":
        return ""
    return line


class GuiInDat(object):
    """Class that handles reading IN.DATs and writing them in a readable
    format
    """

    def __init__(self, filename=""):
        self.__data__ = {}
        self.Run_Description = ""
        if filename:
            self.read_file(filename)

    def __getitem__(self, key):
        key = key.lower()
        if "(" in key:
            # work out which array and which index is being accessed
            array_name, array_index = interpret_array_var(key)
            is_valid(array_name)
            if array_name not in self.__data__:
                # if it's a valid array but not in the dictionary, copy the
                # dictionary value
                self.__data__[array_name] = copy.copy(DICT_DEFAULT[array_name])
            return self.__data__[array_name][array_index]
        else:
            is_valid(key)
            if key not in self.__data__:
                self.__data__[key] = copy.copy(DICT_DEFAULT[key])
            return self.__data__[key]

    def __setitem__(self, key, val):
        key = key.lower()
        if "(" in key:
            array_name, array_index = interpret_array_var(key)
            is_valid(array_name)
            # cast the value to the correct type
            val = cast(key, val)
            if array_name not in self.__data__:
                self.__data__[array_name] = copy.copy(DICT_DEFAULT[array_name])
            self.__data__[array_name][array_index] = val
        else:
            is_valid(key)
            val = cast(key, val)
            self.__data__[key] = val

    def __delitem__(self, key):
        key = key.lower()
        if key in self.__data__:
            del self.__data__[key]

    def is_diff(self, key):
        """Returns true if a variable has a value in this IN.DAT
        different to the default

        """
        if key in self.__data__:
            if self.__data__[key] != DICT_DEFAULT[key]:
                return True
        return False

    def clip(self):
        """Chops elements beyond nvar and neqns off the end of ixc
        and icc
        """
        self["ixc"] = self["ixc"][: self["nvar"]]
        self["icc"] = self["icc"][: self["neqns"]]

    def __str__(self):
        """Convets self to text format"""
        LB = "\n"  # line break character
        ret = ""

        # print the descritption
        for line in self.Run_Description.split("\n"):
            ret += "*****" + line.rstrip() + LB
        ret += LB + make_mod_header("-") + LB + LB

        # constraints
        ret += make_mod_header("Constraint Equations") + LB
        neqns = self["neqns"]
        ret += make_line("neqns", neqns, False) + LB

        # always print in sorted order
        icc_temp = sorted(self["icc"][:neqns])
        for num in range(self["neqns"]):
            const_num = icc_temp[num]
            ret += make_line("icc(" + str(num + 1) + ")", const_num, False)
            ret += " * "
            # get the description of the constraint from lablcc
            ret += get_line_by_int(const_num, DICT_DESCRIPTIONS["lablcc"])
            ret += LB
        ret += LB

        # itervars
        ret += make_mod_header("Iteration Variables") + LB
        nvar = self["nvar"]
        ret += make_line("nvar", nvar, False) + LB
        # always print in sorted order
        ixc_temp = sorted(self["ixc"][:nvar])
        for num in range(self["nvar"]):
            iter_num = ixc_temp[num]
            ubname = "boundu(" + str(iter_num) + ")"
            lbname = "boundl(" + str(iter_num) + ")"
            ret += make_line("ixc(" + str(num + 1) + ")", iter_num, False)
            ret += " * "

            # write a description
            varname = DICT_IXC_SIMPLE[str(iter_num)]
            ret += varname + " "
            ret += get_comment(0.1, DICT_DESCRIPTIONS[varname])
            ret += LB

            # print the lower and upper bounds below
            ret += make_line(lbname, self[lbname], False) + LB
            ret += make_line(ubname, self[ubname], False) + LB
            ret += LB

        # don't print these, they've already been printed
        neverprint = {"neqns", "nvar", "ixc", "icc", "boundl", "boundu"}
        # these are important, always print them
        alwaysprint = {"ioptimz", "minmax"}

        # if isweep isn't 0, always print sweep information
        if self["isweep"] != 0:
            alwaysprint.add("nsweep")
            alwaysprint.add("sweep")

        # print the main body
        for module, arr in DICT_MODULE.items():
            ret += make_mod_header(module) + LB

            for varname in arr:
                # don't print these, they have already been printed
                if varname in neverprint:
                    continue
                # if the value is different or is an important variable
                if self.is_diff(varname) or varname in alwaysprint:
                    val = self[varname]
                    if isinstance(val, list):
                        # don't print the python style [ ] array characters
                        ret += make_line(varname, str(val)[1:-1]) + LB
                    else:
                        ret += make_line(varname, val) + LB

            ret += LB + LB

        return ret

    def add_dict(self, di):
        """Given a dictionary of values, adds them to the IN.DAT. Values can be
        given in string format and will be casted automatically

        """
        for key, value in di.items():
            self[key] = value
        self.clip()

    def readlines(self, inputlines):
        """Iterates through a list of lines and adds values to the dictionary
        or the header

        """
        lines = []
        currentline = ""
        # first we need to make a pass and merge together run-on lines
        # eg IXC = 1, 2, 3, 4
        # 5, 6, 7
        # becomes IXC = 1,2,3,4,5
        for line in inputlines:
            if line == "":
                continue
            if line[:5] == "*****":
                line = line.rstrip()
                self.Run_Description += line[5:] + "\n"
                continue
            line = strip_line(line)
            if line == "":
                continue
            # if a new variable assignment starts on this line
            if "=" in line:
                # save the accumulated line
                lines.append(currentline)
                # start a new line
                currentline = line
            else:
                # otherwise concat this line onto the last with a comma
                # seperator if needed
                if currentline:
                    currentline += "," + line
                else:
                    currentline = line

        lines.append(currentline)
        if lines[0] and "=" not in lines[0]:
            raise ValueError(
                "First line '" + lines[0] + "'" + "does not contain '=' sign"
            )

        self.Run_Description = self.Run_Description.strip()

        # now make the assignments
        for line in lines:
            if line == "":
                continue
            varname, varval = line.split("=")[:2]
            self[varname.strip()] = varval.strip()

        # cut off the end of ixc and icc if neqns and nvar are set low enough
        self.clip()

    def read_file(self, filename):
        """Opens a file and reads it in"""
        self.readlines(open(filename).readlines())

    def write_file(self, filename):
        """Writes to a given filename"""
        file_handle = open(filename, "w")
        file_handle.write(str(self))


# if called as main function, convert an IN.DAT
if __name__ == "__main__":
    PROGDESC = (
        "Reads an IN.DAT file and prints an IN.DAT file of the same "
        + "format as used by the GUI. Redirect to output file using '>'"
    )

    PARSER = argparse.ArgumentParser(description=PROGDESC)
    PARSER.add_argument("inputfile", help="IN.DAT file to read from")
    PARSER.add_argument("outputfile", help="IN.DAT file to read from")
    ARGS = PARSER.parse_args()

    print("Using dictionary for PROCESS v." + str(DICTIONARY_VERSION))
    in_dat_obj = GuiInDat(ARGS.inputfile)
    file_handle = open(ARGS.outputfile, "w")
    file_handle.write(str(in_dat_obj))
    file_handle.close()
