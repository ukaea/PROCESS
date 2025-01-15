"""Create Python-Fortran dictionaries JSON.

This module produces the python_fortran_dicts.json file used by Python in
Process. Make directs the Ford documentation program to read in the Process
source, create a project object which contains the structure of Process used for
documenting the code within Ford, and then creates a pickled file from that
project object. Make then calls this module to create dictionaries of variables
which are used by Python in Process. The dicts are created from hardcoded
values, Process source parsing and the Ford project object, and then dumped into
the JSON file for later use.

Process Python can then call process.io.python_fortran_dicts.get_dicts() to load
the dicts from the saved JSON file and use them.

This ultimately provides Process Python with the ability to access variable
information in the Process Fortran source code.
"""

import argparse
import json
import logging
import pickle
import re
from pathlib import Path

import create_dicts_config
import numpy as np
from python_dicts import get_python_variables

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


class ProjectDictionary(Dictionary):
    # Dicts that rely on the Ford project object
    def __init__(self, name, project, python_variables, value_type):
        Dictionary.__init__(self, name)
        self.project = project  # The Ford project object
        self.python_variables = (
            python_variables  # List of AnnotatedVariables from Python PhysEng models
        )
        self.value_type = value_type
        # The attribute in the project to make a dict for

    def make_dict(self):
        # Assign the variable name key to the value of an attribute of
        # the var object in the project ([var_name] = some_value_of_var)
        for module in self.project.modules:
            for var in module.variables:
                self.dict[self.name][var.name] = getattr(var, self.value_type)

        for annotated_variable in self.python_variables:
            self.dict[self.name][annotated_variable.name] = getattr(
                annotated_variable, self.value_type
            )


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


class VariableDescriptions(ProjectDictionary):
    # Dictionary of variable descriptions
    def __init__(self, project, python_variables):
        ProjectDictionary.__init__(
            self, "DICT_DESCRIPTIONS", project, python_variables, "doc"
        )

    def make_dict(self):
        # Assign the variable name key to the var description
        for module in self.project.modules:
            for var in module.variables:
                desc = getattr(var, self.value_type)

                # If var key doesn't exist, add it
                if (
                    var.name not in self.dict[self.name]
                    or not self.dict[self.name][var.name]
                    and desc
                ):
                    # Only overwrite description if it's falsey and we're
                    # overwriting with something truthy
                    # Guards against multiple declarations in different modules,
                    # when only one declaration is commented
                    self.dict[self.name][var.name] = desc

        for annotated_variable in self.python_variables:
            self.dict[self.name][annotated_variable.name] = annotated_variable.units

    def post_process(self):
        for var_name, var_desc in self.dict[self.name].items():
            # Strip the <p> and </p> tags from the description; not
            # required in output
            desc_sub = re.sub(r"</?p>", r"", var_desc)
            # Some characters are modified in Ford; convert them back here
            # Convert <code> tags back to backticks (`)
            desc_sub = re.sub(r"<\/*code>", r"`", desc_sub)
            # Convert > and < back
            desc_sub = re.sub(r"&gt;", r">", desc_sub)
            desc_sub = re.sub(r"&lt;", r"<", desc_sub)
            self.dict[self.name][var_name] = desc_sub


class DefaultValues(ProjectDictionary):
    """This is a nightmare. It takes combined declared/initialised values from
    Ford's project and combines them with Fortran source regex parsing for
    values that are declared separately and initialised only in an init
    subroutine. e.g.:

    Case 1:
    integer, parameter :: n_radial_array = 50
    Declared and initialised in same line: parameter, does not require
    re-initialisation. Use Ford project's value for this variable.

    Case 2:
    real(kind(1.0D0)) :: fcutfsu
    Declared only (requires re-initialisation)
    fcutfsu = 0.69D0
    Initialised separately in init subroutine
    The Ford project value will be None, but we want to extract the value from
    the init subroutine. This requires Fortran source regex.
    """

    # Dictionary of default values of variables
    def __init__(self, project, python_variables):
        ProjectDictionary.__init__(
            self, "DICT_DEFAULT", project, python_variables, "initial"
        )

    def make_dict(self):
        # Assign the variable name key to the initial value of the variable
        # [var_name] = initial_value
        # Use Ford's project object to attempt to fetch initial values first
        for module in self.project.modules:
            for var in module.variables:
                self.dict[self.name][var.name] = self.process_initial_value(var)

        # Now fetch values in init subroutines to overwrite Ford's initial
        # values if necessary
        self.parse_init_subroutines()

        for annotated_variable in self.python_variables:
            self.dict[self.name][annotated_variable.name] = annotated_variable.obj
            # obj is the initial value also

    def process_initial_value(self, var):
        """Process the initial value of a var from Ford.

        :param var: Ford variable object
        :type var: obj
        :return: initial value of var
        :rtype: None, int, float, str
        """
        # The initial value could be an array that hasn't been picked up by Ford
        # Ford can't handle implicit array initialisation
        # e.g. real(kind(1.0D0)), dimension(14) :: impurity_enrichment = 5.0D0
        # Ford's initial variable value is 5.0D0, but should be an array of
        # 14 5.0D0 values

        # The array dimension attribute
        size_arg = None

        # Check array declaration style, if any
        if var.dimension:
            # "impurity_enrichment(14)" style
            size_arg = var.dimension
        else:
            # "dimension(14) :: impurity_enrichment" style
            for attrib in var.attribs:
                # Only interested in attribute of the form "dimension(size)"
                if attrib.find("dimension") >= 0:
                    size_arg = attrib

        # If it's an array, parse the size
        if size_arg:
            # Check for size being a hardcoded int first
            size = self.find_int_array_dimension(size_arg)

            if size is None:
                # Check for size being a variable
                size = self.find_var_array_dimension(size_arg)

            if size:
                # The variable is an array
                var.dimension = size

                # Modify the initial value if we have one
                if var.initial:
                    if "(/" in var.initial:
                        # Initial value is already an array; convert to str to
                        # list
                        var.initial = var.initial.replace("(/", "")
                        var.initial = var.initial.replace("/)", "")
                        var.initial = var.initial.split("',")
                        # Try to get single quotes back
                        var.initial = [i.strip() + "'" for i in var.initial]
                        var.initial[-1] = var.initial[-1][:-1]
                    else:
                        # Single initial value: replace with array of identical
                        # initial values of length size
                        var.initial = [var.initial for i in range(size)]

        # If the var wasn't an array, return the original initial value
        # If it was an array and had an initial value, return an array of the
        # initial value
        # If it was an array and had no initial value, return None (but with
        # an updated dimension on the var). This will be overwritten by the
        # init subroutine parsing in overwrite_ford_init_value()
        return var.initial

    def find_int_array_dimension(self, attrib):
        # Attempt to find the size of an array given its dimension attribute
        # with a hardcoded integer argument
        size = None
        match = re.search(r"\((\d+)\)", attrib)
        # 1st capturing group matches any number of digits for the size
        if match:
            # dimension argument is a hardcoded number
            size = int(match.group(1))

        return size

    def find_var_array_dimension(self, attrib):
        # Attempt to find the size of an array given its dimension attribute
        # with a variable argument
        size = None

        # Try to match "dimension(size)" or "(size)" style
        match = re.search(r"\((\w+)\)", attrib)
        # 1st capturing group matches any number of word characters
        # for size: i.e. a variable name
        if match:
            # The array size is a variable or an expression
            size_arg = match.group(1)
            # Now look up the initial value of that variable within the current
            # module in the Ford project object
            # TODO expand to all modules
            for module in self.project.modules:
                for var in module.variables:
                    if var.name == size_arg:
                        try:
                            size = int(var.initial)
                            break
                        except ValueError:
                            pass
                        except TypeError:
                            # var.initial likely to be None: caused by arrays
                            # with calculated sizes
                            pass

            # If size_arg doesn't match a var.name or var.initial can't be
            # converted to int: they probably aren't a numerical value and
            # require further evaluation, e.g. var.initial = ngc+2
            # Too complicated and probably not worth it: ignore
        return size

    def post_process(self):
        # Most default values are numbers saved as strings, but some
        # are lists of these
        # Attempt to convert strings to floats: 1.57812D3 to 1578.12
        working_dict = self.dict[self.name]

        for key, value in working_dict.items():
            if value is not None and isinstance(value, str):
                # Guard against None
                # Is it a list?
                if type(value) is list or value[0:2] == "(/":
                    # Ford's arrays begin with "(/"
                    value = self.convert_list_to_floats(value)
                else:
                    value = self.convert_value_to_float(value)

                working_dict[key] = value
            elif value is not None and isinstance(value, np.ndarray):
                working_dict[key] = value.tolist()

    def convert_value_to_float(self, value):
        # Convert a value to float: 1D3 to 1000
        original_value = value

        # Might not convert successfully; in which case return original value
        try:
            value = value.lower()
            value = value.replace("d", "E")
            value = float(value)
            return value
        except ValueError:
            # Failed conversion; don't change anything
            return original_value

    def convert_list_to_floats(self, working_list):
        # Change the formatting for lists
        # Change Ford string "(/1D3, 4D2, 2D7/)" or regular list [1D3, 4D2, 2D7]
        # to Python list of floats [1000, 400, 20000000]
        processed_list = []

        if type(working_list) is not list:
            # Ford list: convert string to an array of strings
            working_list = working_list.replace("(/", "")
            working_list = working_list.replace("/)", "")
            working_list = working_list.split(", ")

        # Convert list values to floats
        for value in working_list:
            value = self.convert_value_to_float(value)
            processed_list.append(value)
        return processed_list

    def parse_init_subroutines(self):
        """Find contents of init subroutines in source files.

        Once found, send on for parsing individually.
        """
        # Fetch all Fortran source files
        sources = sorted(Path(SOURCEDIR).glob("*.f90"))

        for source in sources:
            with open(source) as f:
                lines = f.read()

            # Parse file for init subroutine
            # Match on a "subroutine init_" to "end subroutine" block
            init_match = re.search(
                r"(?:subroutine\sinit_)"  # non-capturing group
                r"(\w+)"  # capture "name" of "subroutine init_name"
                r"(.*?end\ssubroutine)",  # capture all subroutine contents
                lines,
                re.DOTALL,  # Dot matches newline characters (allows multiline matches)
            )

            if init_match:
                module_init_name = init_match.group(1)
                module_init_contents = init_match.group(2)

                # Now process the contents of the individual init subroutine
                self.parse_init_subroutine(module_init_name, module_init_contents)

    def parse_init_subroutine(self, init_name, init_contents):
        """Extract the initial values of the init subroutine variables.

        :param init_name: the module name
        :type init_name: str
        :param init_contents: the init subroutine contents
        :type init_contents: str
        """
        # Extract the init values of the variables from the init subroutine
        # contents
        matches = re.finditer(
            r"\s*(\w+)\s*="  # capture var name before =
            r"\s*(.+?)\s*\n"  # lazily capture value (can be multiline), then \n
            r"(?=(?:(?:\s+\w+\s*=)|"  # positive lookahead, non-capturing group
            # for next var assignment
            r"(?:\s*end)|"  # or end of subroutine
            r"(?:\s*!))|"  # or a comment
            r"(?:\s*if\s))",  # or an if statement
            init_contents,
            re.DOTALL,  # Allow multiline matches
        )

        # Now process captured value
        for match in matches:
            # Capture group 1: variable name
            # Capture group 2: variable value
            var = match.group(1)
            value = match.group(2)

            # Slice off any in-line comment
            value = value.split("!")[0]

            # Remove any \n and \& characters from the value
            # Used mainly for multiline arrays
            value = value.replace("\n", "")
            value = value.replace("&", "")

            # Remove space either side of value
            value = value.strip()

            # If the variable value is "" or '', value will be '""' or "''"
            # Adjust this to be an empty string
            if value in ['""', "''"]:
                value = ""

            # If the value is a Fortran array, convert it to a list
            match = re.search(r"\(\/(.*)\/\)", value)
            # Capture anything between (/ and /)
            # This also deals with shape and reshape() functions
            if match:
                value = match.group(1)
                value = value.replace(" ", "")

                # Some arrays of strings contain commas inside the strings
                # Try to avoid splitting these mid-string by splitting using ','
                comma_list = value.split("','")
                if type(comma_list) is list and len(comma_list) > 1:
                    # Remove stray single quotes
                    comma_list = [i.replace("'", "") for i in comma_list]
                    value = comma_list
                else:
                    # Other just split with a comma
                    # This is used for all number arrays
                    value = value.split(",")

            # Now overwrite the Ford project value, if required
            self.overwrite_ford_init_value(var, value)

    def overwrite_ford_init_value(self, var, value):
        """Overwrite Ford's init value with init subroutine value.

        If required, overwrite Ford's initial value with the value picked up in
        the init subroutine.
        :param var: name of variable
        :type var: str
        :param value: value of var in init subroutine
        :type value: int, float, str, list
        """
        if var not in self.dict[self.name]:
            # If the var is not in the dict (picked up by Ford), discard it
            # This probably means that something was picked up by the init
            # subroutine regex in error
            return
        elif self.dict[self.name][var] is None:
            # Only overwrite the value if Ford has produced a None, which is
            # stored on self.dict
            # Find the var in the Ford project again
            for module in self.project.modules:
                for mod_var in module.variables:
                    if var == mod_var.name and mod_var.dimension:
                        # The var has a dimension, so needs to be
                        # initialised as a list
                        if type(value) is list:
                            assert mod_var.dimension == len(value), (
                                "Array"
                                f" {var} has length {mod_var.dimension} "
                                f"according to Ford, but {len(value)} "
                                "in the init subroutine. Perhaps the "
                                "ford_project.pickle needs to be updated?"
                            )
                            # The value list length and Ford variable dimension
                            # match; just use the list in value
                            # Update Ford project var and self.dict value
                            mod_var.initial = value
                            self.dict[self.name][var] = value
                        else:
                            # value needs to be spread over the length of
                            # the list
                            # Set the Ford project var and the dictionary values
                            mod_var.initial = [value] * mod_var.dimension
                            self.dict[self.name][var] = mod_var.initial

            # If it's not an array, set it to the value in the init subroutine
            if self.dict[self.name][var] is None:
                self.dict[self.name][var] = value


class Modules(ProjectDictionary):
    # Dictionary mapping modules to arrays of its module-level variables
    def __init__(self, project, python_variables):
        ProjectDictionary.__init__(
            self, "DICT_MODULE", project, python_variables, "name"
        )

    def make_dict(self):
        for module in self.project.modules:
            # Create individual module dict
            self.dict[self.name][module.name] = []
            for var in module.variables:
                # Add module-level variables
                self.dict[self.name][module.name].append(var.name)

        for annotated_variable in self.python_variables:
            if annotated_variable.parent not in self.dict[self.name]:
                self.dict[self.name][annotated_variable.parent] = []
            self.dict[self.name][annotated_variable.parent].append(
                annotated_variable.name
            )


def to_type(string):
    """Given a string, attempts to convert the string to a numerical
    value. If the string can't be simply converted, the function
    looks to see if it begins with a integer and returns that (since
    some lines have clarification text after the number). If this
    also fails, return string.strip()
    Args:
         string --> String to be converted
    Returns:
         value --> Either a float, int or string depending on input
    """

    try:
        if "." in string:
            # try a float conversion
            string_mod = string.strip().lower().replace("d", "e")
            return float(string_mod)
        else:
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
    Args:
         file --> Name of file to be read
         regexp --> Regular expression to search for
         flags --> re flags to use in search. Default is re.U which has
                     no effect
    Returns:
         lines --> List of matching lines
    """

    lines = []

    try:
        with open(file, encoding="utf-8") as file_open:
            for line in file_open:
                if re.search(regexp, line, flags):
                    lines.append(line)
            file_open.close()
    except OSError:
        logging.warning("File : %s not found\n", file)
    return lines


def slice_file(file, re1, re2):
    """Returns a slice of a file that is bounded by lines containing a
    substring matching the given regular expressions. The first match
    to re1 in the file marks the start of the slice, the first match
    to re2 after that marks the end of the slice. The list of lines
    returned includes the two bounding lines.
    Args:
         file --> Name of file to read through
         re1 --> Starting regular expression
         re2 --> Ending regular expression
    Returns:
         lines --> List of lines from file between re1 and re2 inclusive
    """

    with open(file, "r", encoding="utf-8") as file:
        filetext = file.readlines()
    start = None
    for i in range(len(filetext)):
        # look for first match
        if re.search(re1, filetext[i]):
            start = i
            break
    if start is None:
        logging.warning("Could not match %s in file %s\n", re1, file)
        return ""
    end = None
    for i in range(start, len(filetext)):
        # look for second match
        if re.search(re2, filetext[i]):
            end = i
            break
    if end is None:
        logging.warning("Could not match %s in file %s\n", re2, file)
        return ""
    # return slice
    return filetext[start : end + 1]


def remove_comments(line):
    """Function to remove comments from a fortran line. Works by simply
    removing everything after the first '!' in the line. This will
    cause problems in the case of '!' characters contained within strings
    so am assuming this won't happen. Need to change this.
    Args:
         line --> Line to strip comments from
    Returns
         modified_line --> Line with comments removed
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
    particular variable. Looks in scan.f90 for mapping from nsweep_no to
    iteration variable name, and uses ixc_full to map variable name to ixc_no.

    Example of a fragment we are looking for:
        case (1)
            aspect = swp(iscn)
            vlabel = 'aspect = ' ; xlabel = 'Aspect_ratio'

    Example dictionary entry:
        DICT_IXC2NSWEEP['1'] = '1'
    """

    ixc2nsweep = {}
    file = SOURCEDIR + "/scan.f90"
    # slice the file to get the switch statement relating to nsweep
    lines = slice_file(file, r"select case \(nwp\)", r"case default")

    # remove extra lines that aren't case(#) or varname = sweep(iscan) lines
    modlines = []
    for line in lines[1:-1]:
        if "case" in line or "swp(iscn)" in line:
            line = remove_comments(line).replace(" ", "")
            modlines.append(line)

    # create a dictionary that maps iteration variable names to ixc_no
    ixc_full = dict_ixc_full()
    ixc_simple_rev = {}
    for num, value in ixc_full.items():
        ixc_simple_rev[value["name"]] = num

    for i in range(len(modlines)):
        # get the number from the case statement
        match = re.match(r"case\((\d+)\)", modlines[i])
        if match:
            num = match.group(1)
            # if the case statement matched, get the variable name
            # from the next line
            match_2 = re.match(r"(.*?)=swp\(iscn\)", modlines[i + 1])
            if not match_2:
                logging.warning("Error in dict_ixc2nsweep\n")
            else:
                name = match_2.group(1)
                if name in ixc_simple_rev:
                    ixcnum = ixc_simple_rev[name]
                    ixc2nsweep[ixcnum] = num

    return ixc2nsweep


def dict_nsweep2ixc():
    """Returns a dict mapping nsweep_no to ixc_no; the inverse of
    dict_ixc2nsweep"""

    # Use dict_ixc2nsweep from output_dict to produce dict_nsweep2ixc
    ixc2nsweep = output_dict["DICT_IXC2NSWEEP"]
    nsweep2ixc = {b: a for a, b in ixc2nsweep.items()}
    return nsweep2ixc


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
    regexp = r"call parse_(real|int|string)_(array|variable)\("
    lines = grep(SOURCEDIR + "/input.f90", regexp)
    for line in lines:
        args = line.split("(")[1]
        name = args.split(",")[1].strip()
        var_type = re.search(regexp, line).group(1)
        scalar = re.search(regexp, line).group(2)
        di[name] = var_type + "_" + scalar
    return di


def dict_icc_full():
    """Function to return a dictionary matching str(icc_no) to a dictionary
    containing the name of that constraint equation. Looks in
    numerics.f90 at !+ad_varc lines in lablcc to get icc_no and
    variable names.

    Example of a lablxc line we are looking for:
        !+ad_varc  <LI> ( 5) * beta

    Example dictionary entry:
        DICT_IXC_FULL['5'] = {'name' : 'beta'}
    """

    di = {}

    # get slice of file from ":: lablcc" to a blank line
    lcctext = slice_file(SOURCEDIR + "/numerics.f90", r"::\slablcc", r"^$")

    regexp = r"""
               !!               #var comment begins with !!

               .*?              #irrelevant stuff until open brackets

               \(\s*(\d+)\s*\)  #an integer in brackets possibly bounded by
                                #whitespace. Capture the number in group 1

               \s*\*?\s*        #whitespace and a possible asterix

               ([\w ]+)        #the name of the variable should be captured
                               #in group 2
              """
    lcc = []
    # ignore first and last lines
    for line in lcctext[1:-1]:
        match = re.search(regexp, line, re.VERBOSE)
        if match:
            num = int(match.group(1))
            name = match.group(2).strip()
            lcc.append(name)
            assert num == len(lcc)

    for i in range(len(lcc)):
        assign = {"name": lcc[i]}
        di[str(i + 1)] = assign

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
    failedlines = []
    regexp = r"call parse_(real|int)_variable\((.*)"
    lines = grep(SOURCEDIR + "/input.f90", regexp)

    for line in lines:
        match = re.search(regexp, line)
        try:
            name = match.group(2).split(",")[1].strip()
            lb = to_type(match.group(2).split(",")[2])
            ub = to_type(match.group(2).split(",")[3])
            if match.group(1) == "real":
                assert isinstance(lb, float)
            else:
                assert isinstance(lb, int)
            assert ub >= lb
            di[name] = {"lb": lb, "ub": ub}
        except (IndexError, AttributeError, AssertionError, TypeError):
            failedlines.append(line)

    if len(failedlines) != 0:
        warn_string = "dict_input_bounds failed to parse:\n"
        for line in failedlines:
            warn_string += f"{line.strip()}\n"
        logging.warning(warn_string)

    return di


def dict_nsweep2varname():
    # This function creates the nsweep2varname dictionary from the fortran code
    # It maps the sweep variable number to its variable name

    di = {}
    file = SOURCEDIR + "/scan.f90"

    # slice the file to get the switch statement relating to nsweep
    lines = slice_file(file, r"select case \(nwp\)", r"case default")

    # remove extra lines that aren't case(#) or varname = sweep(iscan) lines
    modlines = []
    for line in lines[1:-1]:
        if "case" in line or "swp(iscn)" in line:
            line = remove_comments(line).replace(" ", "")
            modlines.append(line)

    for i in range(len(modlines) // 2):
        line1 = modlines[i * 2]
        no = line1.replace("case(", "")
        no = no.replace(")", "")
        line2 = modlines[i * 2 + 1]
        varname = line2.replace("=swp(iscn)", "")
        di[no] = varname

    return di


def dict_ixc_full():
    """Function to return a dictionary matching str(ixc_no) to a dictionary
    containing the name, lower and upper bounds of that variable. Looks in
    numerics.f90 at !+ad_varc lines in lablxc to get ixc_no and
    variable names, and looks at boundu and boundl for upper and
    lower bounds.

    Example of a lablxc line we are looking for:
        lablxc(1) = 'aspect        '

    Example of a boundl line we are looking for:
        boundl(1) = 0.0D0

    Example of a boundu line we are looking for:
        boundu(1) = 1.0D0

    Example dictionary entry:
        DICT_IXC_FULL['5'] = {'name' : 'beta', 'lb' : 0.001, 'ub' : 1.0}
    """

    with open(SOURCEDIR + "/iteration_variables.f90") as myFile:
        lines = myFile.readlines()

    ixc_full = {}

    for lline in lines:
        if "subroutine init_itv_" in lline and "end" not in lline:
            itv_num = lline.split("_")[-1].strip("\n").replace(" ", "")
            ixc_full[itv_num] = {}

    for line in lines:
        if ("lablxc" in line and "=" in line) and (
            "lablxc(i)" not in line and "lablxc(ixc(i))" not in line
        ):
            labl_num = line.split("(")[1].split(")")[0]
            labl = line.split("=")[-1].strip("\n").replace(" ", "").replace("'", "")
            ixc_full[labl_num]["name"] = labl

        if ("boundl(" in line and "=" in line) and (
            "boundl(i)" not in line and "boundl(ixc(i))" not in line
        ):
            boundl_num = line.split("(")[1].split(")")[0]
            boundl_val = line.split("=")[-1].strip("\n").lower().replace("d", "e")
            ixc_full[boundl_num]["lb"] = float(boundl_val)

        if ("boundu(" in line and "=" in line) and (
            "boundu(i)" not in line and "boundu(ixc(i))" not in line
        ):
            boundu_num = line.split("(")[1].split(")")[0]
            boundu_val = line.split("=")[-1].strip("\n").lower().replace("d", "e")
            ixc_full[boundu_num]["ub"] = float(boundu_val)

    return ixc_full


def dict_ixc_bounds():
    # Returns dictionary mapping iteration variable name to bounds
    ixc_full = output_dict["DICT_IXC_FULL"]
    ixc_bounds = {}
    for key, value in ixc_full.items():
        lb = value["lb"]
        ub = value["ub"]
        temp = {"lb": lb, "ub": ub}
        ixc_bounds[value["name"]] = temp

    return ixc_bounds


def dict_ixc_default():
    # Returns dictionary mapping iteration variable name to default value
    ixc_default = {}
    default = output_dict["DICT_DEFAULT"]
    ixc_full = output_dict["DICT_IXC_FULL"]

    for key, value in ixc_full.items():
        name = value["name"]

        if name in default:
            ixc_default[name] = default[name]
        else:
            logging.warning("print_dict_ixc could not find %s in DICT_DEFAULT\n", name)

    return ixc_default


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
    ixc_simple_rev = {b: a for a, b in ixc_simple.items()}

    return ixc_simple_rev


def create_dicts(project):
    # There are 3 sources of dicts: from the Ford project object, from the
    # Fortran source and hardcoded ones from this file

    dict_objects = []
    # Different dict objects, e.g. variable descriptions

    python_variables = get_python_variables()

    # Make dict objects
    # Some dicts depend on other dicts already existing in output_dicts, so
    # be careful if changing the order!
    dict_objects.extend([
        VariableDescriptions(project, python_variables),
        DefaultValues(project, python_variables),
        Modules(project, python_variables),
        HardcodedDictionary("DICT_TF_TYPE", create_dicts_config.DICT_TF_TYPE),
        HardcodedDictionary("DICT_FIMP", create_dicts_config.DICT_FIMP),
        HardcodedDictionary(
            "DICT_OPTIMISATION_VARS", create_dicts_config.DICT_OPTIMISATION_VARS
        ),
        HardcodedDictionary("IFAIL_SUCCESS", create_dicts_config.IFAIL_SUCCESS),
        HardcodedDictionary(
            "PARAMETER_DEFAULTS", create_dicts_config.PARAMETER_DEFAULTS
        ),
        HardcodedDictionary("NON_F_VALUES", create_dicts_config.NON_F_VALUES),
        SourceDictionary("DICT_INPUT_BOUNDS", dict_input_bounds),
        SourceDictionary("DICT_NSWEEP2VARNAME", dict_nsweep2varname),
        SourceDictionary("DICT_VAR_TYPE", dict_var_type),
        SourceDictionary("DICT_ICC_FULL", dict_icc_full),
        SourceDictionary("DICT_IXC2NSWEEP", dict_ixc2nsweep),
        SourceDictionary("DICT_NSWEEP2IXC", dict_nsweep2ixc),
        SourceDictionary("DICT_IXC_FULL", dict_ixc_full),
        SourceDictionary("DICT_IXC_BOUNDS", dict_ixc_bounds),
        SourceDictionary("DICT_IXC_DEFAULT", dict_ixc_default),
        SourceDictionary("DICT_IXC_SIMPLE", dict_ixc_simple),
        SourceDictionary("DICT_IXC_SIMPLE_REV", dict_ixc_simple_rev),
    ])

    # Make individual dicts within dict objects, process, then add to output_dict
    for dict_object in dict_objects:
        dict_object.make_dict()
        dict_object.post_process()
        dict_object.publish()

    # Save output_dict as JSON, to be used by utilities scripts
    with open(DICTS_FILENAME, "w") as dicts_file:
        json.dump(output_dict, dicts_file, indent=4, sort_keys=True)


if __name__ == "__main__":
    # TODO This has been written to cause minimal disruption to the original
    # create_dicts.py. This module would benefit from more class structuring

    # Called from make; parse arguments from make
    parser = argparse.ArgumentParser(description="Create Fortran-Python dictionaries")
    parser.add_argument("fortran_source", help="Fortran source dir")
    parser.add_argument("ford_project", help="The pickled Ford project filename")
    parser.add_argument("dicts_filename", help="The output dicts filename")
    args = parser.parse_args()

    # Load the pickled Ford project
    with open(args.ford_project, "rb") as f:
        project = pickle.load(f)

    # Process source directory
    SOURCEDIR = args.fortran_source
    # The name of the output file (python_fortran_dicts.json)
    DICTS_FILENAME = args.dicts_filename

    # Create the dicts JSON output file
    create_dicts(project)
