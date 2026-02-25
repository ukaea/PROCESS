"""Access the dictionaries for variable information.

This ultimately provides access to variable information that is included in
the python source (e.g. docstrings) or that cannot be dynamically accessed
(e.g. variable initial values).
"""

import ast
import inspect
import logging
import re
from functools import cache
from importlib import import_module
from itertools import pairwise

import numpy as np

from process.core.optimisation.iteration_variables import ITERATION_VARIABLES
from process.init import init_all_module_vars
from process.input import INPUT_VARIABLES
from process.scan import ScanVariables

INPUT_TYPE_MAP = {int: "int", float: "real", str: "string"}

NON_F_VALUES = [
    "f_j_cs_start_pulse_end_flat_top",
    "f_c_plasma_non_inductive",
    "feffcd",
    "f_a_tf_turn_cable_copper",
]
logger = logging.getLogger(__name__)

output_dict = {}
# Dict of nested dicts e.g. output_dict['DICT_DESCRIPTIONS'] =
# {descriptions_dict}
# Dicts stored in output_dict are used to create other derivative dicts


# Classes for the various dictionary types
class Dictionary:
    # Base Dictionary class for all dicts
    def __init__(self, name):
        self.name = name  # Dict name
        self.dict = {}  # Contains the dict
        self.dict[self.name] = {}  # Structures this dict: key = dict name,
        # value = nested dict of variable info

    def make_dict(self):
        # Make the dictionary
        pass

    def post_process(self):
        # Perform any processing after making the dict
        pass

    def publish(self):
        # Add the finished dictionary to the output dict
        output_dict.update(self.dict)


class SourceDictionary(Dictionary):
    # Dictionary created from Fortran source
    def __init__(self, name, dict_creator_func):
        Dictionary.__init__(self, name)
        # Function that creates the dict
        self.dict_creator_func = dict_creator_func

    def make_dict(self):
        # Make entire nested dict from function
        self.dict[self.name] = self.dict_creator_func()


class HardcodedDictionary(Dictionary):
    # Dictionary created from a hardcoded dict in this file
    def __init__(self, name, hardcoded_dict):
        Dictionary.__init__(self, name)
        self.dict[self.name] = None
        # Hardcoded value isn't always a dict; override to None to allow the
        # value to be set to any type
        self.hardcoded_dict = hardcoded_dict

    def make_dict(self):
        # Set the nested value to a hardcoded int, list or dict
        self.dict[self.name] = self.hardcoded_dict


def to_type(string):
    """Given a string, attempts to convert the string to a numerical
    value. If the string can't be simply converted, the function
    looks to see if it begins with a integer and returns that (since
    some lines have clarification text after the number). If this
    also fails, return string.strip()

    Parameters
    ----------
    string :


    Returns
    -------
    :
        Either a float, int or string depending on input
    """

    try:
        if "." in string:
            # try a float conversion
            string_mod = string.strip().lower().replace("d", "e")
            return float(string_mod)
        # try an int conversion
        return int(string.strip())
    except ValueError:
        match = re.match(r"\s*(\d+)", string)
        if match:
            # if the string starts with an integer return that
            return int(match.group(1))

        # otherwise return the string unchanged but with whitespace removed
        return string.strip()


def grep(file, regexp, flags=re.UNICODE):
    """Implements an in-python grep. Returns the lines that match
    as a list.

    Parameters
    ----------
    file:
        Name of file to be read
    regexp:
        Regular expression to search for
    flags:
        re flags to use in search. Default is re.U which has

    Returns
    -------
    lines:
        List of matching lines
    """

    lines = []

    try:
        with open(file, encoding="utf-8") as file_open:
            lines = [line for line in file_open if re.search(regexp, line, flags)]

    except OSError:
        logger.warning("File : %s not found\n", file)
    return lines


def slice_file(file, re1, re2):
    """Returns a slice of a file that is bounded by lines containing a
    substring matching the given regular expressions. The first match
    to re1 in the file marks the start of the slice, the first match
    to re2 after that marks the end of the slice. The list of lines
    returned includes the two bounding lines.

    Parameters
    ----------
    file :
        Name of file to read through
    re1 :
        Starting regular expression
    re2 :
        Ending regular expression

    Returns
    -------
    :
        lines:List of lines from file between re1 and re2 inclusive
    """

    with open(file, encoding="utf-8") as file:
        filetext = file.readlines()
    start = None
    for i in range(len(filetext)):
        # look for first match
        if re.search(re1, filetext[i]):
            start = i
            break
    if start is None:
        logger.warning("Could not match %s in file %s\n", re1, file)
        return ""
    end = None
    for i in range(start, len(filetext)):
        # look for second match
        if re.search(re2, filetext[i]):
            end = i
            break
    if end is None:
        logger.warning("Could not match %s in file %s\n", re2, file)
        return ""
    # return slice
    return filetext[start : end + 1]


def remove_comments(line):
    """Function to remove comments from a fortran line. Works by simply
    removing everything after the first '!' in the line. This will
    cause problems in the case of '!' characters contained within strings
    so am assuming this won't happen. Need to change this.

    Parameters
    ----------
    line:
        Line to strip comments from

    Returns
    -------
    modified_line:
        Line with comments removed

    """
    if "!" in line:
        line = line[: line.find("!")]
    line = line.strip()
    if line == "":
        return line
    if line[-1] == "&":
        line = line[:-1]
    return line


def dict_ixc2nsweep():
    """Returns a dict mapping ixc_no to nsweep_no, if both exist for a
    particular variable.

    Example dictionary entry:
        DICT_IXC2NSWEEP['1'] = '1'
    """

    name_to_nsweep = {
        sv.value.variable_name: sv.value.variable_num for sv in ScanVariables
    }

    # create a dictionary that maps iteration variable names to ixc_no
    ixc_full = dict_ixc_full()
    name_to_ixc = {
        name: ixc
        for ixc, data in ixc_full.items()
        if (name := data["name"]) in name_to_nsweep
    }

    return {ixc: name_to_nsweep[name] for name, ixc in name_to_ixc.items()}


def dict_nsweep2ixc():
    """Returns a dict mapping nsweep_no to ixc_no; the inverse of
    dict_ixc2nsweep
    """

    # Use dict_ixc2nsweep from output_dict to produce dict_nsweep2ixc
    ixc2nsweep = output_dict["DICT_IXC2NSWEEP"]
    return {b: a for a, b in ixc2nsweep.items()}


def dict_var_type():
    """Function to return a dictionary mapping variable name to variable type
    eg. 'real_variable' or 'int_array'. Looks in input.f90 at the process
    functions that read in variables from IN.DAT.

    Example of line we are looking for:
        call parse_real_variable('BETA', beta, 0.0D0, 1.0D0, &

    Example dictionary entry:
        DICT_VAR_TYPE['beta'] = 'real_variable'
    """
    di = {}

    for var_name, config in INPUT_VARIABLES.items():
        var_type = (
            f"{INPUT_TYPE_MAP[config.type]}_{'array' if config.array else 'variable'}"
        )

        di[var_name] = var_type

    return di


def dict_input_bounds():
    """Returns a dictionary matching variable names to dictionary containing
    upper and lower bounds that PROCESS checks variable lies between when
    reading IN.DAT. Looks in input.f90 for parse_real_variable and
    parse_int_variable.

    Example of a line we are looking for:
         call parse_real_variable('BETA', beta, 0.0D0, 1.0D0, &

    Example dictionary entry:
         DICT_INPUT_BOUNDS['beta'] = {'lb' : 0.0, 'ub' : 1.0}
    """
    di = {}

    for var_name, config in INPUT_VARIABLES.items():
        lb, ub = None, None
        if config.range is not None:
            lb, ub = config.range

        elif config.choices is not None and config.type in [int, float]:
            lb = min(config.choices)
            ub = max(config.choices)

        if lb is not None:
            di[var_name] = {"lb": lb, "ub": ub}

    return di


def dict_nsweep2varname():
    return {sv.value.variable_num: sv.value.variable_name for sv in ScanVariables}


def dict_ixc_full():
    """Function to return a dictionary matching str(ixc_no) to a dictionary
    containing the name, lower and upper bounds of that variable.

    Example dictionary entry:
        DICT_IXC_FULL['5'] = {'name' : 'beta', 'lb' : 0.001, 'ub' : 1.0}
    """

    return {
        str(k): {"name": v.name, "lb": v.lower_bound, "ub": v.upper_bound}
        for k, v in ITERATION_VARIABLES.items()
    }


def dict_ixc_bounds():
    # Returns dictionary mapping iteration variable name to bounds
    ixc_full = output_dict["DICT_IXC_FULL"]
    ixc_bounds = {}
    for value in ixc_full.values():
        lb = value["lb"]
        ub = value["ub"]
        temp = {"lb": lb, "ub": ub}
        ixc_bounds[value["name"]] = temp

    return ixc_bounds


def dict_ixc_simple():
    # Returns dictionary mapping ixc no to iteration variable name
    ixc_simple = {}
    ixc_full = output_dict["DICT_IXC_FULL"]
    for key, value in ixc_full.items():
        ixc_simple[key] = value["name"]

    return ixc_simple


def dict_ixc_simple_rev():
    # Returns dictionary mapping iteration variable name to ixc no
    ixc_simple = output_dict["DICT_IXC_SIMPLE"]
    return {b: a for a, b in ixc_simple.items()}


# cache the output of get_dicts so that it is never re-calculated in a given
# process run.
@cache
def get_dicts():
    """Constructs the dictionaries which contain information about every PROCESS variable.

    WARNING: this function must be used carefully because it re-initialises the PROCESS state
    """
    dict_objects = []
    # Different dict objects, e.g. variable descriptions

    init_all_module_vars()
    # Make dict objects
    # Some dicts depend on other dicts already existing in output_dicts, so
    # be careful if changing the order!
    dict_objects.extend([
        HardcodedDictionary("NON_F_VALUES", NON_F_VALUES),
        HardcodedDictionary("DICT_DEFAULT", {}),
        HardcodedDictionary("DICT_MODULE", {}),
        HardcodedDictionary("DICT_DESCRIPTIONS", {}),
        SourceDictionary("DICT_INPUT_BOUNDS", dict_input_bounds),
        SourceDictionary("DICT_NSWEEP2VARNAME", dict_nsweep2varname),
        SourceDictionary("DICT_VAR_TYPE", dict_var_type),
        SourceDictionary("DICT_IXC2NSWEEP", dict_ixc2nsweep),
        SourceDictionary("DICT_NSWEEP2IXC", dict_nsweep2ixc),
        SourceDictionary("DICT_IXC_FULL", dict_ixc_full),
        SourceDictionary("DICT_IXC_BOUNDS", dict_ixc_bounds),
        SourceDictionary("DICT_IXC_SIMPLE", dict_ixc_simple),
        SourceDictionary("DICT_IXC_SIMPLE_REV", dict_ixc_simple_rev),
    ])

    # Make individual dicts within dict objects, process, then add to output_dict
    for dict_object in dict_objects:
        dict_object.make_dict()
        dict_object.post_process()
        dict_object.publish()

    for module_name in import_module("process.data_structure").__all__:
        if module_name == "__init__.py":
            continue
        module = import_module(f"process.data_structure.{module_name.split('.', 1)[0]}")

        module_tree = ast.parse(inspect.getsource(module))
        initial_values_dict = {}
        variable_names = []
        var_names_and_descriptions = {}
        dict_module_entry = {}
        variable_types = {}

        # get the variable names and initial values
        for node in ast.walk(module_tree):
            if isinstance(node, ast.AnnAssign):
                # for each variable in the file, get the initial value
                # (either is None, or value initialised in init_variables fn)
                # set default to be None if variable is not being initialised eg if you
                # just have `example_double: float` instead of `example_double: float = None`
                initial_value = getattr(module, node.target.id)
                # JSON doesn't like np arrays
                if type(initial_value) is np.ndarray:
                    initial_value = initial_value.tolist()
                initial_values_dict[node.target.id] = initial_value
                # get the variable name and add to variable_names list
                var_name = node.target.id
                variable_names.append(var_name)
                # Now want to get the types of these variables
                if isinstance(node.annotation, ast.Subscript):
                    if node.annotation.value.id == "list":
                        if node.annotation.slice.id == "str":
                            var_type = "string_array"
                        elif node.annotation.slice.id == "float":
                            var_type = "real_array"
                        elif node.annotation.slice.id == "int":
                            var_type = "int_array"
                        elif node.annotation.slice.id == "bool":
                            var_type = "bool_array"
                        else:
                            raise TypeError(
                                f"The type annotation of variable {node.target.id} is "
                                f"{node.annotation.value.id}[{node.annotation.slice.id}], and "
                                "this is not recognised. Please change your type annotation for "
                                "this variable. PROCESS recognises the following type annotations: "
                                "list[float], list[int], list[str], list[bool]."
                            )
                else:
                    if node.annotation.id == "float":
                        var_type = "real_variable"
                    elif node.annotation.id == "int":
                        var_type = "int_variable"
                    elif node.annotation.id == "str":
                        var_type = "str_variable"
                    elif node.annotation.id == "bool":
                        var_type = "bool_variable"
                    else:
                        raise TypeError(
                            f"The type annotation of variable {node.target.id} is "
                            f"{node.annotation.id}, and this is not recognised. Please change your "
                            "type annotation for this variable. PROCESS recognises the following "
                            "type annotations: float, int, str, bool."
                        )

                variable_types[node.target.id] = var_type
        # get the variable descriptions
        # need to check for pairs of ast.AnnAssign followed by an ast.Expr - this is the form of
        # a variable being declared followed by a docstring expression. can get these var descriptions
        # from here, and if there is no ast.Expr immediately after an ast.AnnAssign then this var does not
        # have a docstring and so set the description to be ""
        for node1, node2 in pairwise(module_tree.body):
            if isinstance(node1, ast.AnnAssign) and isinstance(node2, ast.Expr):
                # if docstring immediately follows the variable declaration, add docstring to descriptions dict
                var_names_and_descriptions[node1.target.id] = node2.value.value
            if isinstance(node1, ast.AnnAssign) and not isinstance(node2, ast.Expr):
                # if no docstring for variable, have a blank description
                var_names_and_descriptions[node1.target.id] = ""
        # check if last entry of ast.body is declaring a var. if it is then this var has no description and will be missing
        # from var_names_and_descriptions. need to add to var_names_and_descriptions dict
        lastvar = module_tree.body[-1]

        if (
            isinstance(lastvar, ast.AnnAssign)
            and lastvar not in var_names_and_descriptions
        ):
            var_names_and_descriptions[lastvar.target.id] = ""

        dict_module_entry[module_name] = variable_names

        output_dict["DICT_MODULE"].update(dict_module_entry)
        output_dict["DICT_DEFAULT"].update(initial_values_dict)
        output_dict["DICT_DESCRIPTIONS"].update(var_names_and_descriptions)
        output_dict["DICT_VAR_TYPE"].update(variable_types)

    return output_dict
