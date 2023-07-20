#!/usr/bin/env python3
"""This script creates several dictionaries and lists for use by the GUI.

   GUI_MODULE --> OrderedDict that sets the ordering of the main body of the UI.
                  Keys are module names, values are lists of
                  (variable_type, variable_name, shortdescription, description)
                  tuples

                  variable_type can be one of:
                  "float", "int", "array", "string"

                  shortdescription is the text displayed next to the variable

                  description is the text displayed when the shortdescription
                  if hovered over

                  description and shortdescription are obtained from
                  DICT_DESCRIPTION in process_dicts, but with removal of
                  characters such as " that could cause problems

   GUI_LABLXC -->  OrderedDict that sets iteration variable ordering. Keys
                   are ixc no., values are
                   (name, shortdesc, desc) tuples where 'name' is name
                   of variable

   GUI_LABLCC -->  OrderedDict that sets constraint equation ordering. Keys
                   are icc no., values are constraint descriptions.

   Tom Miller 09/14

"""

import re
from collections import OrderedDict
import logging
import copy
import argparse
from create_dicts import print_dict
from process_io_lib.process_dicts import (
    DICT_DEFAULT,
    DICT_IXC_FULL,
    DICT_DESCRIPTIONS,
    DICT_MODULE,
)


def print_list(li, name, comment=""):
    """Prints a list with format:
    #comment
    name = [
      li[0],
      li[1],
      ...
    ]

    """
    if len(li) == 0:
        return
    print("\n#" + comment)
    print(name + " = [")
    for value in li:
        print("\t", value, ",")
    print("]")


def parse_labl(labl):
    """Parses a 'labl' format description from DICT_DESCRIPTION
    eg ( 1) * aspect
       ( 2) * bt
    Returns a list of (number, checked, desc) tuples where 'number' is
    number in brackets, 'checked' indicates whether an asterix is present
    and 'desc' is everything that follows

    """
    ret = []
    for line in labl.split("\n"):
        regex = r"""\(([\d ]+)\)     #capture the digits and whitespace
                                    #in brackets at the beginning

                   \s*(\*?)\s*      #is there a * before the description?
                   (.*)
                """
        ma = re.match(regex, line, re.VERBOSE)
        if ma:
            num = int(ma.group(1))
            if ma.group(2) == "*":
                checked = True
            else:
                checked = False
            desc = ma.group(3).strip()
            if desc == "UNUSED":
                continue
            ret.append((num, checked, desc))
    return ret


def get_desc(name):
    """Gets an entry from DICT_DESCRIPTIONS and returns the short style
    description and the long style description with character " removed

    """

    desc = DICT_DESCRIPTIONS[name]
    # make sure no " characters are passed to html
    if '"' in desc:
        logging.warning('Removing " characters from ' + name + " description")
        desc = desc.replace('"', "")

    shortdesc = desc.split("\n")[0]
    return shortdesc, desc


def main():
    """Main program routine. Creates and prints the dictionaries"""
    temp_dict_module = copy.deepcopy(DICT_MODULE)
    dict_module = {}

    for modulename, value in temp_dict_module.items():
        if modulename.endswith(" Variables"):
            modulename = modulename[:-10]
        dict_module[modulename] = value

    # get rid of inertial confinement and reverse field pinch settings
    del dict_module["Ife"]
    del dict_module["Rfp"]

    # get rid of ixc, icc, boundl, boundu
    # these are dealt with seperately
    dict_module["Numerics"].remove("ixc")
    dict_module["Numerics"].remove("icc")
    dict_module["Numerics"].remove("neqns")
    dict_module["Numerics"].remove("nvar")
    dict_module["Numerics"].remove("boundl")
    dict_module["Numerics"].remove("boundu")

    # get rid of runtitle
    # dict_module["Global"].remove("runtitle")

    gui_module = OrderedDict()
    # make gui_module
    for module, varlist in dict_module.items():
        tuplist = []
        for varname in varlist:
            desctup = get_desc(varname)
            if isinstance(DICT_DEFAULT[varname], list):
                tuplist.append(("array", varname) + desctup)
            elif isinstance(DICT_DEFAULT[varname], str):
                tuplist.append(("string", varname) + desctup)
            elif isinstance(DICT_DEFAULT[varname], float):
                tuplist.append(("float", varname) + desctup)
            elif isinstance(DICT_DEFAULT[varname], int):
                tuplist.append(("int", varname) + desctup)
            else:
                continue

        gui_module[module] = tuplist

    # set of all variable names in main body, not printed
    gui_module_var_set = set()
    for varlist in dict_module.values():
        gui_module_var_set.update(varlist)

    gui_lablxc = OrderedDict()
    # make the LABLXC dict
    for num, dummy, dummy in parse_labl(DICT_DESCRIPTIONS["lablxc"]):
        name = DICT_IXC_FULL[str(num)]["name"]
        # check that the name is included in the description given in lablxc
        assert name in dummy

        if "obsolete" in dummy.lower():
            logging.warning(
                " Iteration variable "
                + name
                + " is marked obsolete in lablxc so is being excluded"
            )
            continue

        # tftort should be marked obsolete
        if name == "tftort":
            continue

        if name not in gui_module_var_set:
            logging.warning(
                " Iteration variable "
                + name
                + " does not appear in main body of gui_module so is being excluded"
            )
            continue

        # The description given in DICT_DESCRIPTION is better than in lablxc
        shortdesc, desc = get_desc(name)
        # descriptions like (iteration variable 2) are redundant when in
        # the GUI
        shortdesc = shortdesc.split(" (iter")[0]
        desctup = (shortdesc, desc)
        gui_lablxc[num] = (name,) + desctup

    gui_lablcc = OrderedDict()
    # Create GUI_LABLCC
    for num, dummy, desc in parse_labl(DICT_DESCRIPTIONS["lablcc"]):
        desc = desc.capitalize().replace('"', "")
        if "Unused" in desc:
            logging.warning(
                " Constraint "
                + str(num)
                + " is marked unused "
                + "so is being excluded"
            )
            continue
        gui_lablcc[num] = desc

    # Print the header
    header = """
\"\"\"
This file contains dictionaries for use by the PROCESS GUI.
GUI_MODULE      : Sets the ordering of main body variable
GUI_LABLXC      : Sets the ordering of iteration variables
GUI_LABLCC      : Sets the ordering of constraint equations

\"\"\"

from collections import OrderedDict
    """
    print(header)

    # print everything
    print("\n#List that sets order of main body variables")
    print("GUI_MODULE = OrderedDict()")
    for modname, varlist in sorted(gui_module.items()):
        print("GUI_MODULE['" + modname + "'] = [")
        for tup in varlist:
            print("\t", tup, ",")
        print("\t]")

    comment = "List that sets ixc ordering"
    print_dict(gui_lablxc, "GUI_LABLXC", comment, dict_type="OrderedDict()")

    comment = "List that sets icc ordering"
    print_dict(gui_lablcc, "GUI_LABLCC", comment, dict_type="OrderedDict()")


if __name__ == "__main__":

    desc = (
        "Produces dictionaries for use by the GUI. Prints to stdout. "
        + "Redirect to file using '>'"
    )
    PARSER = argparse.ArgumentParser(description=desc)

    PARSER.parse_args()

    main()
