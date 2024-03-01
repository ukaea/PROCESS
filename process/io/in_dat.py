"""

    File for reading IN.DATs
    Version 2 (mainly for use with IN.DAT created from UI)

    James Morris
    CCFE
    11/12/14

    Notes:
        + 24/11/2021: Global dictionary variables moved within the functions
                    to avoid cyclic dependencies. This is because the dicts
                    generation script imports, and inspects, process.
"""

from re import sub
import subprocess
from sys import stderr
from process.io.python_fortran_dicts import get_dicts

# ioptimz values
ioptimz_des = {
    "-2": "for no optimisation, no VMCOM or HYBRD",
    "-1": "for no optimisation HYBRD only",
    "0": "for HYBRD and VMCON (not recommended)",
    "1": "for optimisation VMCON only",
}


def fortran_python_scientific(var_value):
    """Convert FORTRAN scientific notation to Python notation

    :param var_value: variable value
    :return: value with 'd/D' notation swapped for 'e/E' notation
    """
    return var_value.replace("D", "e").replace("d", "e")


def remove_empty_lines(lines):
    """Function to remove empty lines from list.

    :param lines: list of lines (type=list)
    :return: list of lines with empty lines removed (type=list)
    """
    return [line for line in lines if line != "\n"]


def is_title(line):
    """Function to determine if line is title line

    :param line: line from IN.DAT
    :return: True/False if line is a title or header.
    """
    return line[:2] == "*-" or line[:3] == "***" or line[0] == "$"


def is_comment(line):
    """Function to determine if line is a commented line
    :param line: line from IN.DAT
    :return: True/False if line is a commented line
    """
    return line[0] == "*"


def is_iteration_variable(name):
    """Function to determine if item is an iteration variable

    :param name: Name of the variable
    :return: True/False if 'ixc' is in name
    """
    return "ixc" in name


def is_constraint_equation(name):
    """Function to determine if item is constraint equation

    :param name: Name of the variable
    :return: True/False if 'icc' is in name
    """
    return "icc" in name


def is_bound(name):
    """Function to determine if item is a bound

    :param name: Name of the variable
    :return: True/False if 'bound' is in name
    """
    return "bound" in name


def is_array(name):
    """Function to determine if item is an array

    :param name: Name of the variable
    :return: True/False if name is of an array
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()
    if "(" in name:
        name = name.split("(")[0]

    try:
        return "array" in dicts["DICT_VAR_TYPE"][name]
    except KeyError:
        # print("Warning:", name, "is not in DICT_VAR_TYPE")
        return False


def find_line_type(line):
    """Function to find line type

    :param line: Line to find type of
    :return: Return string describing the line type
    """

    # Split variable name from line
    name = line.split("=")[0].strip("")

    # If the line is just the title for a section
    if is_title(line):
        return "Title"

    # If the line is a commented line
    elif is_comment(line):
        return "Comment"

    # Else if the line contains a constraint equation
    elif is_constraint_equation(name):
        return "Constraint Equation"

    # Else if the line contains an iteration variable
    elif is_iteration_variable(name):
        return "Iteration Variable"

    # Else if the line contains a bound statement
    elif is_bound(name):
        return "Bound"

    # Else all other arrays
    elif is_array(name):
        return "Array"

    # Else the line contains an regular parameter
    else:
        return "Parameter"


def find_parameter_group(name):
    """Function to find the module which the parameter belongs to.

    :param name: Parameter name
    :return: Return the module the parameter belongs to
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Search DICT_MODULES for parameter
    for key in dicts["DICT_MODULE"].keys():
        if name in dicts["DICT_MODULE"][key]:
            return key


def write_title(title, out_file):
    """Function to write title line to file with fixed width

    :param title: The name of the section
    :param out_file: Output file for new IN.DAT
    :return: Nothing
    """

    out_file.write("\n")
    # Insert title name into line of fixed width
    formatted_title = "*" + title.center(50, "-") + "*\n"

    # Output to file
    out_file.write(formatted_title)
    out_file.write("\n")


def get_constraint_equations(data):
    """Create the constraint equation information.

    Use the constraint equation numbers from IN.DAT to find the comment
    associated with them in the source dictionary, then return both.

    :param dict data: Data dictionary for the IN.DAT information
    :return: dict of the constraint numbers and their comments
    :rtype: dict
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    constraints = {}

    # List of constraint equation numbers in IN.DAT
    constraint_numbers = data["icc"].value

    # Find associated comments and create constraint dict
    for constraint_number in constraint_numbers:
        comment = dicts["DICT_ICC_FULL"][str(constraint_number)]["name"]
        constraints[constraint_number] = comment

    return constraints


def write_constraint_equations(data, out_file):
    """Function to write constraint equation information to file

    :param data: Data dictionary for the IN.DAT information
    :param out_file: Output file for new IN.DAT
    :return: Nothing
    """

    # Header
    write_title("Constraint Equations", out_file)

    # Fetch dict of constraint equation information
    constraints = get_constraint_equations(data)

    for number, comment in constraints.items():
        constraint_line = "icc = {0} * {1}\n".format(number, comment)
        out_file.write(constraint_line)


def get_iteration_variables(data):
    """Create the iteration variable information.

    Use the iteration variable numbers from IN.DAT to find the comment
    associated with them in the source dictionary. Then check the information
    from the IN.DAT file to see if upper and/or lower bounds are present, and
    what the value is. Return all this information for each variable.

    :param dict data: Data dictionary for the IN.DAT information
    :return: variable number, comment, upper and/or lower bounds if present
    :rtype: dict
    """
    variables = {}

    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # List of variable numbers in IN.DAT
    variable_numbers = data["ixc"].value

    # Create variable dicts
    for variable_number in variable_numbers:
        variable = {}

        comment = dicts["DICT_IXC_SIMPLE"][
            str(variable_number).replace(",", ";").replace(".", ";").replace(":", ";")
        ]
        variable["comment"] = comment

        # Set bounds if there are any
        if str(variable_number) in data["bounds"].value:
            # Lower bound
            if "l" in data["bounds"].value[str(variable_number)].keys():
                variable["lower_bound"] = (
                    data["bounds"].value[str(variable_number)]["l"].replace("e", "d")
                )

            # Upper bound
            if "u" in data["bounds"].value[str(variable_number)].keys():
                variable["upper_bound"] = (
                    data["bounds"].value[str(variable_number)]["u"].replace("e", "d")
                )

        variables[variable_number] = variable

    return variables


def write_iteration_variables(data, out_file):
    """Function to write iteration variable information to file

    :param data: Data dictionary for the IN.DAT information
    :param out_file: Output file for new IN.DAT
    :return: Nothing
    """

    # Header
    write_title("Iteration Variables", out_file)

    # Fetch dict of iteration variable information
    variables = get_iteration_variables(data)

    for number, info in variables.items():
        variable_line = "ixc = {0} * {1}\n".format(number, info["comment"])
        out_file.write(variable_line)

        if "lower_bound" in info:
            lower_bound_line = "boundl({0}) = {1}\n".format(number, info["lower_bound"])
            out_file.write(lower_bound_line)

        if "upper_bound" in info:
            upper_bound_line = "boundu({0}) = {1}\n".format(number, info["upper_bound"])
            out_file.write(upper_bound_line)


def get_parameters(data, use_string_values=True):
    """Create the parameter information.

    Use the parameters from IN.DAT to produce a dict of name, value and
    comment (optional) for each parameter. This takes the form:
    parameters[module][param_name] = param_value, or {value: param_value,
    comment: param_comment}.

    If use_string_values == True, in the returned dict, ensure that
    the values of the parameters are of type string. This is used for writing
    new IN.DAT files. If False, store the values as their original type, e.g.
    int, float etc. This is used for validating the input file.

    :param dict data: Data dictionary for the IN.DAT information
    :param bool use_string_values: If True, store all parameter values as
    strings. If False, preserve parameter value type.
    :return: dict of parameters containing names, values and comments
    :rtype: dict
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    source_variables = {}
    # dict of all module-level variables in source, grouped by module
    parameters = {}
    # dict of all parameters set in input file, grouped by module
    # Include neqns to allow eq and ineq constraints to be defined in produced IN.DAT
    exclusions = ["nvar", "icc", "ixc"]
    # Parameters to exclude

    # Change module keys from DICT_MODULE: replace spaces with underscores and
    # lower the case for consistency in the formatted_input_data_dict
    # Store the key-modified dict in source_variables
    for old_module_key, variables in dicts["DICT_MODULE"].items():
        new_module_key = old_module_key.replace(" ", "_").lower()
        source_variables[new_module_key] = variables

    # Store parameters in order defined in DICT_MODULE
    # TODO: is order important? Not using ordered dicts any more
    for module, module_variables in source_variables.items():
        parameters[module] = {}

        # Loop over all module-level variables in source for given module
        for item in module_variables:
            # Store a variable in parameters dict if it's in the IN.DAT file
            # (and not in the exclusion list). Store parameter name and value
            if item not in exclusions and item in data.keys():
                if item == "fimp":
                    for k in range(len(data["fimp"].get_value)):
                        name = "fimp({0})".format(str(k + 1).zfill(1))
                        value = data["fimp"].get_value[k]
                        parameters[module][name] = value

                elif item == "ioptimz":
                    name = item
                    ioptimz = {}
                    iop_val = data["ioptimz"].get_value
                    iop_comment = ioptimz_des[str(iop_val)]
                    ioptimz["value"] = iop_val
                    ioptimz["comment"] = iop_comment
                    parameters[module][name] = ioptimz

                elif item == "zref":
                    for j in range(len(data["zref"].get_value)):
                        name = "zref({0})".format(str(j + 1).zfill(1))
                        value = data["zref"].get_value[j]
                        parameters[module][name] = value

                elif item == "impurity_enrichment":
                    for m in range(len(data["impurity_enrichment"].get_value)):
                        name = "impurity_enrichment({0})".format(str(m + 1).zfill(1))
                        value = data["impurity_enrichment"].get_value[m]
                        parameters[module][name] = value

                elif "vmec" in item:
                    name = item
                    value = data[item].value
                    parameters[module][name] = value

                else:
                    parameter = {}

                    if use_string_values:
                        # Store the parameter value as a string
                        # (data[item].value is a string)
                        line_value = data[item].value
                        line_string = ""
                        # if parameter is a list only output values comma separated
                        if isinstance(line_value, list):
                            for val in line_value:
                                line_string += str(val) + ", "
                            line_value = line_string.rstrip(", ")

                        if isinstance(line_value, str):
                            split_line = line_value.split(" ")
                        try:
                            float(split_line[0])
                            if len(split_line) > 1:
                                line_value = ", ".join([entry for entry in split_line])
                        except Exception:
                            pass

                    else:
                        # Store the parameter value preserving its type
                        # (data[item].get_value preserves the parameter's type,
                        # e.g. float)
                        line_value = data[item].get_value

                    name = item
                    parameter["value"] = line_value
                    parameter["comment"] = data[item].comment.split("\n")[0]
                    # Only use first line of comment to avoid lots of info
                    parameters[module][name] = parameter

    return parameters


def write_parameters(data, out_file):
    """Write parameters to file

    :param data: Data dictionary for the IN.DAT information
    :param out_file: Output file for new IN.DAT
    :return: Nothing
    """
    filter_list = ["fimp(", "zref(", "imp_rich", "vmec"]
    # Special parameters that require different formatting
    parameters = get_parameters(data)

    for module in parameters:
        # Write module heading: format to be more readable again
        formatted_module = module.replace("_", " ").title()
        write_title("{0}".format(formatted_module), out_file)

        # Write out parameters for this module
        for parameter, info in parameters[module].items():
            if any(var_name in parameter for var_name in filter_list):
                # No justification formatting if parameter is in filter list
                parameter_line = "{0} = {1}\n".format(parameter, info)
            else:
                # All other parameters
                # Left justification set to 8 to allow easier reading
                # info can currently be either a value or a dict
                if (
                    type(info) is dict
                    and "value" in info
                    and (type(info.get("comment")) is str)
                ):
                    parameter_line = "{0} = {1} * {2}\n".format(
                        parameter.ljust(8), info["value"], info["comment"]
                    )
                else:
                    parameter_line = "{0} = {1}\n".format(parameter.ljust(8), info)

            # Finally write the line
            out_file.write(parameter_line)


def add_iteration_variable(data, variable_number):
    """Function to add iteration variable to IN.DAT data dictionary

    :param data: Data dictionary for the IN.DAT information
    :param variable_number: Iteration variable number to add
    :return: Nothing
    """

    # Check the variable number is not already in the iteration variable list
    if variable_number not in data["ixc"].value:
        data["ixc"].value.append(variable_number)
        data["ixc"].value.sort()

    else:
        print(
            "Variable number {0} already in iteration variable list".format(
                variable_number
            )
        )


def remove_iteration_variable(data, variable_number):
    """Function to remove iteration variable from the IN.DAT data dictionary

    :param data: Data dictionary for the IN.DAT information
    :param variable_number: Iteration variable number to remove
    :return: Nothing
    """

    # Check the variable is in the iteration variable list
    if variable_number in data["ixc"].value:
        data["ixc"].value.remove(variable_number)
        data["ixc"].value.sort()
    else:
        print(
            "Variable number {0} not in iteration variable list".format(variable_number)
        )


def add_constraint_equation(data, equation_number):
    """Function to add constraint equation to the IN.DAT data dictionary

    :param data: Data dictionary for the IN.DAT information
    :param equation_number: Constraint equation number to add
    :return: Nothing
    """

    # Check the constraint is not already in the constraint equation list
    if equation_number not in data["icc"].value:
        data["icc"].value.append(equation_number)
        data["icc"].value.sort()

    else:
        print(
            "Equation number {0} already in constraint equations list".format(
                equation_number
            )
        )


def remove_constraint_equation(data, equation_number):
    """Function to remove a constraint equation from the IN.DAT data
    dictionary

    :param data: Data dictionary for the IN.DAT information
    :param equation_number: Constraint equation number to remove
    :return: Nothing
    """

    # Check the constraint is in
    if equation_number in data["icc"].value:
        data["icc"].value.remove(equation_number)
        data["icc"].value.sort()

    else:
        print(
            "Equation number {0} not in constraint equations list".format(
                equation_number
            )
        )


def add_parameter(data, parameter_name, parameter_value):
    """Function to add/change parameter to the IN.DAT data dictionary

    :param data: Data dictionary for the IN.DAT information
    :param parameter_name: Name of the parameter to add
    :param parameter_value: Value of the parameter to add
    :return: Nothing
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Check that the parameter is not already in the dictionary
    if parameter_name not in data.keys():
        parameter_group = find_parameter_group(parameter_name)
        if "fimp" in parameter_name:
            comment = dicts["DICT_DESCRIPTIONS"]["fimp"]
        else:
            try:
                comment = dicts["DICT_DESCRIPTIONS"][parameter_name]
            except KeyError:
                # The dictionary doesn't recognise the variable name
                print(
                    "Warning: Description for {0}".format(parameter_name),
                    "specified in IN.DAT not in dictionary.",
                    file=stderr,
                )
                comment = ""

        param_data = INVariable(
            parameter_name, parameter_value, "Parameter", parameter_group, comment
        )

        data[parameter_name] = param_data

    # If it is already in there change the value to the new value
    else:
        data[parameter_name].value = parameter_value

    # def __delitem__(self, key):
    #    del self.__dict__[key]

    # def keys(self):
    #    return self.__dict__.keys()


def remove_parameter(data, parameter_name):
    """Function to remove parameter from the IN.DAT data dictionary

    :param data: Data dictionary for the IN.DAT information
    :param parameter_name: Name of the parameter to remove
    :return: Nothing
    """

    # Check that the parameter exists in the data dictionary
    if parameter_name in data.keys():
        del data[parameter_name]

    # Inform the user that the parameter requested for deletion isn;t in the
    # data dictionary
    else:
        print("Parameter {0} not in IN.DAT".format(parameter_name))


def change_array(data, name, array_id, array_val):
    """Function to change value in generic array

    :param data: Data dictionary for the IN.DAT information
    :param name: generic array name
    :param array_id: generic array index
    :param array_val: generic array value
    :return:
    """

    data[name].value[array_id] = array_val
    if isinstance(data[name].value[array_id], str):
        tmp = list(data[name].value.split(","))
        tmp[array_id] = array_val
        new_val = str(tmp).strip("[").strip("]").replace("'", "")
        data[name].value = new_val
    else:
        data[name].value[array_id] = array_val


def add_bound(data, bound, bound_type, bound_value):
    """Function to add/change a bound to the bounds entry in the IN.dat data
    dictionary

    :param data: Data dictionary for the IN.DAT information
    :param bound: Bound number associated with iteration variable number to
                  change
    :param bound_type: States whether bound is upper of lower bound
    :param bound_value: New value of the bound
    :return: Nothing
    """

    # Put bound type into lower cases for consistency
    bound_type = bound_type.lower()

    # if the bound is not in the bounds dictionary initialise an empty
    # dictionary and assign new bound
    if bound not in data["bounds"].value.keys():
        data["bounds"].value[bound] = dict()
        data["bounds"].value[bound][bound_type] = str(bound_value)

    # If bound already exists change value
    elif bound in data["bounds"].value.keys():
        data["bounds"].value[bound][bound_type] = str(bound_value)

    # Bound not recognised.
    else:
        print("Bound {0} not recognised. Check type == string".format(bound))


def remove_bound(data, bound, bound_type):
    """Function to remove a bound from the bounds entry in the IN.DAT data
    dictionary

    :param data: Data dictionary for the IN.DAT information
    :param bound: Bound number associated with iteration variable number to
                  change
    :param bound_type: States whether bound is upper or lower bound
    :return: Nothing
    """

    # use local variable for cleanliness
    bounds = data["bounds"].value

    # If the bound exists (and is of the correct type) in the bounds dictionary
    if bound in bounds.keys() and bound_type in bounds[bound].keys():
        del bounds[bound][bound_type]

        # if the bound number is now an empty dictionary delete it also
        if len(bounds[bound].keys()) == 0:
            del bounds[bound]


def fortran_float_to_py(f: str) -> str:
    if not isinstance(f, str):
        return f
    p = sub(r"([0-9]+\.[0-9]+)(?:D|d)([0-9]+)", r"\1e\2", f)

    return p


def parameter_type(name, value):
    """Function to return value in correct format for altering values etc.

    :param name: Name of parameter to check type
    :param value: Value of parameter to format
    :return: Formatted value
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Find parameter type from PROCESS dictionary
    param_type = dicts["DICT_VAR_TYPE"][name]

    # Check if parameter is a list
    if isinstance(value, list):
        if value[-1] == "":
            value = value[:-1]

        # Real array parameter
        if "real_array" in param_type:
            return [
                item if item is None else float(fortran_float_to_py(item))
                for item in value
            ]
            # Convert list to floats, but not if the value is None

        # Integer array parameter
        elif "int_array" in param_type:
            return [item if item is None else int(item) for item in value]
            # Convert list to ints, but not if the value is None

    # Check if parameter is a string
    elif isinstance(value, str):
        # If a real variable just convert to float
        if "real_variable" in param_type:
            # Prepare so float conversion succeeds
            value = value.lower()
            value = value.replace("d", "e")
            return float(value)

        # If a real array split and make a float list
        elif "real_array" in param_type:
            # Prepare so float conversion succeeds
            value = value.lower()
            value = value.replace("d", "e")
            value = value.split(",")
            if value[-1] == "":
                value = value[:-1]
            return [float(item) for item in value]

        # If an integer variable convert to integer
        elif "int_variable" in param_type:
            return int(value)

        # If an integer array split and make an integer list
        elif "int_array" in param_type:
            value = value.split(",")
            if value[-1] == "":
                value = value[:-1]
            return [int(item) for item in value]

        # If type unknown return original value
        else:
            return value

    # If type is other return original value
    else:
        return value


def variable_constraint_type_check(item_number, var_type):
    """Function to put input into correct format for altering values etc.

    :param item_number: Number associated with variable or constraint
    :param var_type: States whether item is iteration variable or constraint
                     equation
    :return: Formatted item_number
    """

    # Check if item is in string format
    if isinstance(item_number, str):
        # Try evaluate and convert to an integer. Warning if number is rounded
        try:
            # eval should produce int of float otherwise raise the ValueError
            item_number = eval(item_number)

            # Integer
            if isinstance(item_number, int):
                return item_number

            # number must be float if exception not raised
            elif item_number.is_integer():
                return int(item_number)

            # rounded float number with warning
            else:
                print(
                    "Value {0} for {1} not an integer. Value rounded to {2}."
                    " Check!".format(item_number, var_type, int(item_number))
                )
                return int(item_number)

        except ValueError:
            print(
                "Value {0} for {1} not valid. Check value!".format(
                    item_number, var_type
                ),
                file=stderr,
            )

    # Check if item is in float format
    elif isinstance(item_number, float):
        # If integer convert to float and return
        if item_number.is_integer():
            return int(item_number)

        # If not an integer warn of rounding and return rounded integer
        else:
            print(
                "Value {0} for {1} not an integer. Value rounded to {2}. "
                "Check!".format(item_number, var_type, int(item_number))
            )
            return int(item_number)

    # If already an integer return unchanged
    elif isinstance(item_number, int):
        return item_number

    # Value not recognised
    else:
        print(
            "Value {0} for {1} not a recognised format. Check value!".format(
                item_number, var_type
            )
        )


def variable_bound_check(bound_number, bound_type):
    """Function to put bound_number and bound_type into correct format

    :param bound_number: Bound number to check
    :param bound_type: States whether bound is upper or lower bound
    :return: Formatted bound number and bound type
    """

    # put bound type into lower case for consistency
    bound_type = bound_type.lower()

    # check if bound is one of the allowed values if not warn user
    if bound_type not in ["l", "u", "upper", "lower"]:
        print(
            "Bound type '{0}' not recognised. Must be one of "
            "['u', 'l', 'U', 'L', 'lower', 'upper', 'LOWER', 'UPPER']".format(
                bound_type
            )
        )

    # if bound is given as full word shorten for consistency for dictionary
    # keys
    elif bound_type in ["upper", "lower"]:
        if bound_type == "upper":
            bound_type = "u"
        elif bound_type == "lower":
            bound_type = "l"

    # Format bound number value
    # If a string return unchanged
    if isinstance(bound_number, str):
        return bound_number, bound_type

    # If an int convert to string
    elif isinstance(bound_number, int):
        return str(bound_number), bound_type

    # If a float convert to str but warn of rounding when changing from float
    # to int
    elif isinstance(bound_number, float):
        if bound_number.is_integer():
            return int(bound_number), bound_type
        else:
            bound_number = int(bound_number)
            print(
                "Bound number {0} not an integer. "
                "Value rounded to {1}".format(bound_number, int(bound_number))
            )
            return bound_number, bound_type


class INVariable(object):
    def __init__(self, name, value, v_type, parameter_group, comment):
        """Class to stores the information of a single variable from the
        IN.DAT file

        :param name: Item name
        :param value: Item value
        :param v_type: Type of item
        :param parameter_group: PROCESS variable group item belongs to
        :param comment: Comment for item
        :return: Nothing
        """

        self.name = name
        self.value = value
        self.v_type = v_type
        self.parameter_group = parameter_group
        self.comment = comment

    @property
    def get_value(self):
        """Return value in correct format"""
        if self.v_type != "Bound":
            return parameter_type(self.name, self.value)
        else:
            return self.value


class InDat(object):
    """
    Class 'InDat' for handling IN.DAT data. It handles

        - Reading IN.DAT files
        - Writing IN.DAT files
        - Storing information in dictionary for use in other codes
        - Alterations to IN.DAT
    """

    def __init__(self, filename="IN.DAT", start_line=0):
        """Initialise class

        :param filename: Name of input IN.DAT
        :param start_line: Line to start reading from
        :return: Nothing
        """

        self.filename = filename
        self.start_line = start_line

        # Initialise parameters
        self.in_dat_lines = list()
        self.data = dict()
        self.unrecognised_vars = []
        self.duplicates = []  # Duplicate variables

        # read in IN.DAT
        if filename is not None:
            self.read_in_dat()

    def read_in_dat(self):
        """Function to read in 'self.filename' and put data into dictionary
        'self.data'
        """

        # Read in IN.DAT
        with open(self.filename) as indat:
            self.in_dat_lines = indat.readlines()
            self.in_dat_lines = self.in_dat_lines[self.start_line :]

        # Remove empty lines from the file
        self.in_dat_lines = remove_empty_lines(self.in_dat_lines)

        for line in self.in_dat_lines:
            # Put everything in lower case
            if "vmec" not in line.split("=")[0].lower():
                l_line = line.lower()
            else:
                l_line = line

            # find the type of the line:
            # [constraint equation, iteration variable, bound, parameter]
            line_type = find_line_type(l_line)

            # Ignore title, header and commented lines
            if line_type != "Title" and line_type != "Comment":
                try:
                    # for non-title lines process line and store data.
                    self.process_line(line_type, l_line)
                except KeyError:
                    print(
                        "Warning: Line below is causing a problem. Check "
                        "that line in IN.DAT is valid. Line skipped!\n{0}".format(line),
                        file=stderr,
                    )

                    # Store the first part of the unrecognised line (probably a
                    # variable name) as an unrecognised var
                    unrecognised_var = line.split("=")[0].strip()
                    self.unrecognised_vars.append(unrecognised_var)

    def process_line(self, line_type, line):
        """Function to process the line and return the appropriate INVariable
        object

        :param line_type: Type of information the line contains
        :param line: Line from IN.DAT to process
        :return: Nothing
        """

        # Load dicts from dicts JSON file
        dicts = get_dicts()

        # Create bound variable class using INVariable class if the bounds entry
        # doesn't exist
        if "bounds" not in self.data.keys():
            self.data["bounds"] = INVariable(
                "bounds", dict(), "Bound", "Bound", "Bounds"
            )

        # Constraint equations
        if line_type == "Constraint Equation":
            self.process_constraint_equation(line)

        # Iteration_variables
        elif line_type == "Iteration Variable":
            self.process_iteration_variables(line)

        # Bounds
        elif line_type == "Bound":
            self.process_bound(line)

        # Arrays
        elif line_type == "Array":
            # Create geneneric array variable class using INVariable class,
            # if it does not yet exist
            line_commentless = line.split("*")[0]
            array_name = line_commentless.split("(")[0]
            empty_array = dicts["DICT_DEFAULT"][array_name]
            # empty_array is what the array is initialised to when it is
            # declared in the Fortran. If the array is declared but not
            # initialised until later in a separate "init" subroutine, then
            # empty_array will be None. This is a side-effect of the need to
            # re-initialise Fortran arrays in a separate subroutine.

            if empty_array is None:
                # Array is declared but not initialised at the same time;
                # convert it to an empty list because this better reflects
                # what it is, and allows further list operations
                empty_array = []

            if array_name not in self.data.keys():
                parameter_group = find_parameter_group(array_name)

                # Get parameter comment/description from dictionary
                comment = (
                    dicts["DICT_DESCRIPTIONS"][array_name]
                    .replace(",", ";")
                    .replace(".", ";")
                    .replace(":", ";")
                )

                # Copy the default array from the dicts
                empty_array_copy = empty_array[:]
                # Copy empty_array to decouple reference of self.data to
                # DICT_DEFAULT; don't want changes to data to change the
                # defaults
                self.data[array_name] = INVariable(
                    array_name, empty_array_copy, array_name, parameter_group, comment
                )

            self.process_array(line, empty_array)

        # Parameter
        else:
            self.process_parameter(line)

    def process_parameter(self, line):
        """Function to process parameter entries in IN.DAT

        :param line: Line from IN.DAT to process
        :return: Nothing
        """
        # Load dicts from dicts JSON file
        dicts = get_dicts()

        # Remove comment from line to make things easier
        no_comment_line = line.split("*")[0].split("=")

        # Parameter name
        name = no_comment_line[0].strip()

        # Parameter value
        if len(no_comment_line[-1].split(",")) > 1:
            try:
                value = no_comment_line[1].strip()
            except IndexError:
                print(
                    "Error when reading IN.DAT file on line",
                    no_comment_line,
                    "\n Please note, that our Python Library cannot cope with",
                    " variable definitions on multiple lines.",
                    file=stderr,
                )
                exit()
        else:
            try:
                value = no_comment_line[1].strip().replace(",", "")
            except IndexError:
                print(
                    "Error when reading IN.DAT file on line",
                    no_comment_line,
                    "\n Please note, that our Python Library cannot cope with",
                    " variable definitions on multiple lines.",
                    file=stderr,
                )
                exit()

        # Find group of variables the parameter belongs to
        parameter_group = find_parameter_group(name)

        # Get parameter comment/description from dictionary
        comment = (
            dicts["DICT_DESCRIPTIONS"][name]
            .replace(",", ";")
            .replace(".", ";")
            .replace(":", ";")
        )

        # Check that the parameter isn't a duplicate; does the key already
        # exist?
        if self.data.get(name):
            self.add_duplicate_variable(name)

        # Populate the IN.DAT dictionary with the information
        self.data[name] = INVariable(name, value, "Parameter", parameter_group, comment)

    def process_constraint_equation(self, line):
        """Function to process constraint equation entry in IN.DAT

        :param line: Line from IN.DAT to process
        :return: Nothing
        """

        # Remove comment from line to make things easier
        no_comment_line = line.split("*")[0].split("=")

        # If the line contains a constraint equation in the form ICC(#)
        if "(" in no_comment_line[0] and ")" in no_comment_line[0]:
            constraints = [no_comment_line[1].strip()]

        # Else the line contains a list of constraint equations icc = #, #, #
        else:
            constraints = no_comment_line[1].strip().split(",")
            if "" in constraints:
                constraints.remove("")

        # List of new constraints read in
        value = [int(item.strip()) for item in constraints]

        # Populate data dictionary with constraint equations
        # If constraint equation list not already in data dictionary initialise
        # INVariable class
        if "icc" not in self.data.keys():
            self.data["icc"] = INVariable(
                "icc",
                value,
                "Constraint Equation",
                "Constraint Equation",
                "Constraint Equations",
            )

        else:
            # Add constraint equation numbers to list
            for item in constraints:
                if int(item) not in self.data["icc"].value:
                    self.data["icc"].value.append(int(item))
                else:
                    # Duplicate constraint equation number
                    self.add_duplicate_variable("icc = {0}".format(item))
            # Don't sort the constraints! Preserves what's eq, what's ineq;
            # first neqns are eqs, rest are ineqs
            # self.data["icc"].value.sort()

    def process_iteration_variables(self, line):
        """Function to process iteration variables entry in IN.DAT

        :param line: Line from IN.DAT to process
        :return: Nothing
        """

        # Remove comment from line to make things easier
        no_comment_line = line.split("*")[0].split("=")

        # If the line contains an iteration variable in the form IXC(#)
        if "(" in no_comment_line[0] and ")" in no_comment_line[0]:
            iteration_variables = [no_comment_line[1].strip()]

        # Else the line contains a list of iteration variables IXC = #, #, #
        else:
            iteration_variables = no_comment_line[1].strip().split(",")
            if "" in iteration_variables:
                iteration_variables.remove("")

        # List of new constraints read in
        value = [int(item.strip()) for item in iteration_variables]

        # Populate data dictionary with iteration variables
        # If iteration variables list not already in data dictionary initialise
        # INVariable class
        if "ixc" not in self.data.keys():
            self.data["ixc"] = INVariable(
                "ixc",
                value,
                "Iteration Variable",
                "Iteration Variable",
                "Iteration Variables",
            )

        else:
            # Add iteration variable to list
            for item in iteration_variables:
                if int(item) not in self.data["ixc"].value:
                    self.data["ixc"].value.append(int(item))
                else:
                    # Duplicate iteration variable
                    self.add_duplicate_variable("ixc = {0}".format(item))
            self.data["ixc"].value.sort()

    def process_bound(self, line):
        """Function to process bound entries in IN.DAT

        :param line: Line from IN.DAT to process
        :return: Nothing
        """

        # Initialise bound type
        bound_type = None

        # Remove comment from line to make things easier
        no_comment_line = line.split("*")[0].split("=")

        # If upper bound
        if "boundu" in no_comment_line[0]:
            bound_type = "u"

        # If lower bound
        elif "boundl" in no_comment_line[0]:
            bound_type = "l"

        # Get bound information
        bound = (
            no_comment_line[0].strip("boundl").replace("(", "").replace(")", "").strip()
        )
        bound_value = (
            no_comment_line[1]
            .strip()
            .replace(",", "")
            .replace("d", "e")
            .replace("D", "e")
        )

        # If bound not in the bound dictionary then add entry for bound with an
        # empty dictionary
        if bound not in self.data["bounds"].value.keys():
            self.data["bounds"].value[bound] = dict()
        elif self.data["bounds"].value[bound].get(bound_type):
            # Duplicate bound
            self.add_duplicate_variable("bound{0}({1})".format(bound_type, bound))

        # Populate self.data dictionary with bound information
        self.data["bounds"].value[bound][bound_type] = bound_value

    def process_array(self, line, empty_array):
        """Function to process generic array

        :param line: Line from IN.DAT to process
        :param empty_array: Default array for this array name
        :return: nothing
        """
        # TODO This is a mess after the Fortran module variable
        # re-initialisation work. It is hard to see how this can be improved
        # as the regex method of Fortran inspection (i.e. Python-Fortran
        # dictionaries) is increasingly untenable. An attempt is made here,
        # but with a view to the dictionaries method being dropped in future
        # in light of increasing Python f2py conversion.

        if "*" in line:
            line_commentless = line.split("*")[0]
        else:
            line_commentless = line

        name = line_commentless.split("(")[0]
        index = int(line_commentless.split("(")[1].split(")")[0]) - 1
        # The -1 assumes all Fortran arrays begin at 1: they don't in Process!
        # This bug would cause a Fortran index of 0 to be interpreted as a
        # Python index of -1 (last element). This didn't cause any exceptions,
        # but would obviously set the wrong list index. This now throws
        # exceptions when trying to access [-1] of an empty list []. Hence it
        # must be guarded against with the following:
        if index == -1:
            index = 0
        # This will cause Fortran index 0 and 1 to overwrite the same Python
        # index of 0, which is clearly awful. However, it is equally bad to the
        # previous case above, where Fortran [0] would overwrite Python [-1].
        # There isn't a way of reconciling this without knowing whether the
        # Fortran array begins at 0 or 1.
        # TODO Either enforce Fortran arrays that start at 1 throughout, or
        # find a way of determining the starting index of the array

        value = line_commentless.split("=")[-1].replace(",", "")

        # Many arrays are now declared but not initialised until the "init"
        # subroutine is run in each Fortran module, to allow re-initialisation
        # on demand for a new run etc.
        # In Fortran, we have a declared but uninitialised array of a given
        # length. This results in an empty list in the value attribute here.
        # However, we need to set the value for a given index. Therefore
        # make a list of Nones, so that we can set the value for a given
        # index. We don't know the length of the Fortran array (its value is
        # []), so we have to make the Python list as long as this Fortran index.
        # This way, at least the list is long enough for this Fortran index.
        # TODO Again, this requires a more robust solution

        list_len = len(self.data[name].value)
        # Length of Python list in self.data
        max_list_index = list_len - 1
        # The greatest index in that Python list
        array_len = index + 1
        # The Fortran array must be at least this long

        # If the default array is an empty list, make a list of Nones of
        # matching length to the Fortran array
        if len(empty_array) == 0:
            empty_array = [None] * array_len

        # If the Pyhton list is an empty list, make a list of Nones of
        # matching length to the Fortran array
        if list_len == 0:
            self.data[name].value = [None] * array_len
        # Check the Python list is long enough to store the new index
        elif max_list_index < index:
            # The Python list already has a length > 0, but it's not long enough
            # for this new index. Copy the old list and extend it with Nones so
            # that it is as long as the new index
            len_diff = index - max_list_index
            self.data[name].value.extend([None] * len_diff)
        # Array has already been set to default values (empty_array)
        # Need a way of checking for duplicate initialisations
        # If value is changing from default to custom value, then interpret as
        # first initialisation. If value is changing from one custom value to
        # another, then interpret as a duplicate initialisation
        elif self.data[name].value[index] != empty_array[index]:
            # This array index is already not its default value; any further
            # change must be a duplicate initialisation
            fortran_index = index + 1
            # Index begins at 1!
            self.add_duplicate_variable("{0}({1})".format(name, fortran_index))

        # Now we are sure that the Python list index exists, set its value
        self.data[name].value[index] = eval(fortran_python_scientific(value))

    def add_duplicate_variable(self, name):
        """Records duplicate variables in the input file.

        If a var is initialised more than once in the input file, the last
        value persists, but the overwriting is recorded here.
        :param name: The name of the variable being duplicated
        :type var: str
        """
        self.duplicates.append(name)

    def add_iteration_variable(self, variable_number):
        """Function to add iteration variable to IN.DAT data dictionary

        :param variable_number: Iteration variable number to add
        :return: Nothing
        """

        # format iteration variable number
        variable_number = variable_constraint_type_check(
            variable_number, "iteration variable"
        )
        # add iteration variable to IN.DAT data dictionary
        add_iteration_variable(self.data, variable_number)

    def remove_iteration_variable(self, variable_number):
        """Function to remove iteration variable to IN.DAT data dictionary

        :param variable_number: Iteration variable number to remove
        :return: Nothing
        """

        # format iteration variable number
        variable_number = variable_constraint_type_check(
            variable_number, "iteration variable"
        )
        # remove iteration variable from IN.DAT data dictionary
        remove_iteration_variable(self.data, variable_number)

    def add_constraint_equation(self, equation_number):
        """Function to add constraint equation to IN.DAT data dictionary

        :param equation_number: Constraint equation number to add
        :return: Nothing
        """

        # format constraint equation number
        equation_number = variable_constraint_type_check(
            equation_number, "constraint equation"
        )

        # add constraint equation to IN.DAT data dictionary
        add_constraint_equation(self.data, equation_number)

    def remove_constraint_equation(self, equation_number):
        """Function to remove a constraint equation from IN.DAT data
        dictionary

        :param equation_number: Constraint equation number to remove
        :return: Nothing
        """

        # format constraint equation number
        equation_number = variable_constraint_type_check(
            equation_number, "constraint equation"
        )

        # remove constraint equation from IN.DAT data dictionary
        remove_constraint_equation(self.data, equation_number)

    def add_parameter(self, parameter_name, parameter_value):
        """Function to add/change parameter to IN.DAT data dictionary

        :param parameter_name: Name of parameter to add/change
        :param parameter_value: Value of parameter to add/change
        :return: Nothing
        """

        # add/change parameter to/in IN.DAT data dictionary
        add_parameter(self.data, parameter_name, parameter_value)

    def change_array(self, array_name, array_index, array_value):
        """Function to change value in generic array

        :param array_name: name of generic array to change
        :param array_index: index of generic array to change
        :param array_value: value to change array entry to
        :return:
        """

        # change generic array value in IN.DAT data dictionary
        change_array(self.data, array_name, array_index, array_value)

    def remove_parameter(self, parameter_name):
        """Function to remove parameter from IN.DAT data dictionary

        :param parameter_name: Name of parameter to remove
        :return: Nothing
        """

        # remove parameter from IN.DAT data dictionary
        remove_parameter(self.data, parameter_name)

    def add_bound(self, bound, bound_type, bound_value):
        """Function to add/change a bound in IN.DAT data dictionary

        :param bound: Bound number to add/change
        :param bound_type: States whether bound is upper or lower bound
        :param bound_value: Value of bound to add/change
        :return: Nothing
        """

        # format bound number and bound type
        bound, bound_type = variable_bound_check(bound, bound_type)

        # add/change bound to/in IN.DAT data dictionary
        add_bound(self.data, bound, bound_type, bound_value)

    def remove_bound(self, bound, bound_type):
        """Function to remove a bound from IN.DAT data dictionary

        :param bound: Bound number to remove
        :param bound_type: States whether bound is upper or lower bound
        :return: Nothing
        """

        # format bound number and bound type
        bound, bound_type = variable_bound_check(bound, bound_type)

        # remove bound from IN.DAT data dictionary
        remove_bound(self.data, bound, bound_type)

    def write_in_dat(self, output_filename="new_IN.DAT"):
        """Function to write data to output file called 'output_filename'"""

        # create and open output file
        output = open(output_filename, "w")

        # Write Header
        write_title("", output)

        # Write Constraint Equations
        write_constraint_equations(self.data, output)

        # Write Iteration Variables
        write_iteration_variables(self.data, output)

        # Write parameters
        write_parameters(self.data, output)

        # close file
        output.close()

    @property
    def number_of_constraints(self):
        """
        Return number of itvars
        """
        return len(self.data["icc"].value)

    @property
    def number_of_itvars(self):
        """
        Return number of itvars
        """
        return len(self.data["ixc"].value)


def test(f):
    """Test function

    :param f: file name to test
    """
    try:
        i = InDat(filename=f)
        i.write_in_dat(output_filename="test_out_IN.DAT")
        subprocess.call(["rm", "test_out_IN.DAT"])
        return True
    except Exception:
        return False


class StructuredInputData:
    """Combines structured input file data and methods for accessing it.

    This class uses the InDat class to read in an IN.DAT input file and create
    a structured dict from it, combined with methods for accessing the data in
    that dict. This is useful for exporting the IN.DAT data to other modules,
    for example the input_validator module.

    The data dict stored within this class differs from the INDat.data dict in
    that it has a hierarchical structure that more closely resembles the IN.DAT
    input file, and hence it is more human-readable and ready to output to file.

    An example of the structure:
    self.data["parameters"]["physics_variables"]["ishape"]["value"] = 0
    """

    def __init__(self, filename="IN.DAT"):
        """Use InDat to create the data dict.

        Use InDat to read in an input file and store the data, then construct a
        structured input data dict from it.

        :param filename: input data filename, defaults to "IN.DAT"
        :type filename: str, optional
        """
        self.data = {}
        # Structured input data dict

        in_dat = InDat(filename)

        self.data["constraint_equations"] = get_constraint_equations(in_dat.data)
        self.data["iteration_variables"] = get_iteration_variables(in_dat.data)
        self.data["parameters"] = get_parameters(in_dat.data, use_string_values=False)

        self.unrecognised_vars = in_dat.unrecognised_vars
        self.duplicates = in_dat.duplicates
        # Duplicate initialisations in the input file

    def get_param(self, var_name):
        """Get a parameter's dict from the data.

        :param var_name: The name of the parameter to be returned
        :type var_name: str
        :return: Dictionary of parameter's information, or None if not found
        :rtype: dict
        """
        modules = self.data["parameters"]

        for module_dict in modules.values():
            var_dict = module_dict.get(var_name)
            if var_dict is not None:
                break

        return var_dict

    def is_param_defined(self, var_name):
        """Check if a parameter is defined or not in the input data.

        :param var_name: Name of the parameter to be checked
        :type var_name: str
        :return: True if defined, False if not
        :rtype: bool
        """
        defined = False

        if self.get_param(var_name):
            defined = True

        return defined

    def get_param_value(self, var_name):
        """Gets the value of a parameter from the input data.

        :param var_name: The name of the parameter
        :type var_name: str
        :return: Value of that parameter; can be any type. None if not found.
        :rtype: int, float
        """
        value = None
        var_dict = self.get_param(var_name)
        if var_dict:
            value = var_dict.get("value")

        return value


if __name__ == "__main__":
    # i = InDat(filename="../../modified_demo1_a31_rip06_2014_12_15.IN.DAT")
    i = InDat(filename="IN.DAT")
    # print(i.data["ixc"].value)
    # print(i.number_of_constraints)
    # print(i.number_of_itvars)
    # print(i.data["fimp"].value)
    # print(i.data["ipfloc"].value)
    # i.change_fimp(3, 0.5)
    # print(i.data["zref"].value)
    # i.change_zref(3, 0.5)
    # i.remove_constraint_equation(2.5)
    # i.add_constraint_equation("3.0")
    # i.add_constraint_equation("2")
    # i.add_iteration_variable(103)
    # i.add_iteration_variable("2")
    # i.add_iteration_variable(7.5)
    # i.add_iteration_variable("5.5")
    # i.remove_iteration_variable(2)
    # i.remove_iteration_variable("3")
    # i.remove_iteration_variable(4.5)
    # i.remove_iteration_variable("6.5")
    # # Add bound will change the bound value if it already exists
    # i.add_bound(103, "upper", 5.0)
    # i.remove_bound(2, "upper")
    # # Add parameter will change the parameter value if it already exists
    # i.add_parameter("blnktthdsd", 0.5)
    # i.add_parameter("iavail", 1)
    # i.remove_parameter("blnkithsddd")
    # i.remove_parameter("blnkith")
    # i.add_parameter("sweep", [3.0, 3.0])
    # print(i.data["bounds"].get_value)
    i.write_in_dat()
