"""

File for reading IN.DATs
Version 2 (mainly for use with IN.DAT created from UI)

Notes:
    + 24/11/2021: Global dictionary variables moved within the functions
                to avoid cyclic dependencies. This is because the dicts
                generation script imports, and inspects, process.
"""

from re import sub
from sys import stderr

from process.core.exceptions import ProcessValidationError
from process.core.solver.constraints import ConstraintManager
from process.io.data_structure_dicts import get_dicts

# ioptimz values
ioptimz_des = {
    "-2": "for no optimisation, no VMCOM or HYBRD",
    "-1": "for no optimisation HYBRD only",
    "0": "for HYBRD and VMCON (not recommended)",
    "1": "for optimisation VMCON only",
}


def fortran_python_scientific(var_value):
    """Convert FORTRAN scientific notation to Python notation

     Parameters
     ----------
     var_value:
         variable value

     Returns
     -------
    :
         value with 'd/D' notation swapped for 'e/E' notation
    """
    return var_value.replace("D", "e").replace("d", "e")


def remove_empty_lines(lines):
    """Function to remove empty lines from list.

     Parameters
     ----------
     lines:
         list of lines (type=list)

     Returns
     -------
    :
         list of lines with empty lines removed (type=list)
    """
    return [line for line in lines if line.strip(" ") != "\n"]


def is_title(line):
    """Function to determine if line is title line

     Parameters
     ----------
     line:
         line from IN.DAT

     Returns
     -------
    :
         True/False if line is a title or header.
    """
    return line[:2] == "*-" or line[:3] == "***" or line[0] == "$"


def is_comment(line):
    """Function to determine if line is a commented line

     Parameters
     ----------
     line:
         line from IN.DAT

     Returns
     -------
    :
         True/False if line is a commented line
    """
    return line[0] == "*"


def is_iteration_variable(name):
    """Function to determine if item is an iteration variable

     Parameters
     ----------
     name:
         Name of the variable

     Returns
     -------
    :
         True/False if 'ixc' is in name
    """
    return "ixc" in name


def is_constraint_equation(name):
    """Function to determine if item is constraint equation

     Parameters
     ----------
     name:
         Name of the variable

     Returns
     -------
    :
         True/False if 'icc' is in name
    """
    return "icc" in name


def is_bound(name):
    """Function to determine if item is a bound

     Parameters
     ----------
     name:
         Name of the variable

     Returns
     -------
    :
         True/False if 'bound' is in name
    """
    return "bound" in name


def is_array(name):
    """Function to determine if item is an array

     Parameters
     ----------
     name:
         Name of the variable

     Returns
     -------
    :
         True/False if name is of an array
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

     Parameters
     ----------
     line:
         Line to find type of

     Returns
     -------
    :
         Return string describing the line type
    """

    # Split variable name from line
    name = line.split("=")[0].strip("")

    # If the line is just the title for a section
    if is_title(line):
        return "Title"

    # If the line is a commented line
    if is_comment(line):
        return "Comment"

    # Else if the line contains a constraint equation
    if is_constraint_equation(name):
        return "Constraint Equation"

    # Else if the line contains an iteration variable
    if is_iteration_variable(name):
        return "Iteration Variable"

    # Else if the line contains a bound statement
    if is_bound(name):
        return "Bound"

    # Else all other arrays
    if is_array(name):
        return "Array"

    # Else the line contains an regular parameter
    return "Parameter"


def find_parameter_group(name):
    """Function to find the module which the parameter belongs to.

     Parameters
     ----------
     name:
         Parameter name

     Returns
     -------
    :
         Return the module the parameter belongs to
    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Search DICT_MODULES for parameter
    for key in dicts["DICT_MODULE"]:
        if name in dicts["DICT_MODULE"][key]:
            return key
    return None  # Explicit return


def write_title(title, out_file):
    """Function to write title line to file with fixed width

    Parameters
    ----------
    title:
        The name of the section
    out_file:
        Output file for new IN.DAT

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

    Parameters
    ----------
    dict:
        data: Data dictionary for the IN.DAT information
    data :


    Returns
    -------
    dict
        dict of the constraint numbers and their comments
    """
    constraints = {}

    # List of constraint equation numbers in IN.DAT
    constraint_numbers = [int(i) for i in data["icc"].value]

    # Find associated comments and create constraint dict
    for constraint_number in constraint_numbers:
        constraint = ConstraintManager.get_constraint(constraint_number)

        if constraint is None:
            raise ProcessValidationError(
                f"Constraint equation {constraint_number} requested but has not been registered"
            )

        # TODO: we should store a short description of each constraint that we can use here.
        # Currently, no such information exists.
        constraints[constraint_number] = ""

    return constraints


def write_constraint_equations(data, out_file):
    """Function to write constraint equation information to file

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    out_file:
        Output file for new IN.DAT

    """

    # Header
    write_title("Constraint Equations", out_file)

    # Fetch dict of constraint equation information
    constraints = get_constraint_equations(data)

    for number, comment in constraints.items():
        constraint_line = f"icc = {number} * {comment}\n"
        out_file.write(constraint_line)


def get_iteration_variables(data: dict):
    """Create the iteration variable information.

    Use the iteration variable numbers from IN.DAT to find the comment
    associated with them in the source dictionary. Then check the information
    from the IN.DAT file to see if upper and/or lower bounds are present, and
    what the value is. Return all this information for each variable.

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information

    Returns
    -------
    :
        variable number, comment, upper and/or lower bounds if present
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
            if "l" in data["bounds"].value[str(variable_number)]:
                variable["lower_bound"] = (
                    data["bounds"].value[str(variable_number)]["l"].replace("e", "d")
                )

            # Upper bound
            if "u" in data["bounds"].value[str(variable_number)]:
                variable["upper_bound"] = (
                    data["bounds"].value[str(variable_number)]["u"].replace("e", "d")
                )

        variables[variable_number] = variable

    return variables


def write_iteration_variables(data, out_file):
    """Function to write iteration variable information to file

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    out_file:
        Output file for new IN.DAT
    """

    # Header
    write_title("Iteration Variables", out_file)

    # Fetch dict of iteration variable information
    variables = get_iteration_variables(data)

    for number, info in variables.items():
        variable_line = f"ixc = {number} * {info['comment']}\n"
        out_file.write(variable_line)

        if "lower_bound" in info:
            lower_bound_line = f"boundl({number}) = {info['lower_bound']}\n"
            out_file.write(lower_bound_line)

        if "upper_bound" in info:
            upper_bound_line = f"boundu({number}) = {info['upper_bound']}\n"
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

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    use_string_values:
        If True, store all parameter values as
        strings. If False, preserve parameter value type.

    Returns
    -------
    dict
        dict of parameters containing names, values and comments
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
            if item not in exclusions and item in data:
                if item == "f_nd_impurity_electrons":
                    for k in range(len(data["f_nd_impurity_electrons"].get_value)):
                        name = f"f_nd_impurity_electrons({str(k + 1).zfill(1)})"
                        value = data["f_nd_impurity_electrons"].get_value[k]
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
                        name = f"zref({str(j + 1).zfill(1)})"
                        value = data["zref"].get_value[j]
                        parameters[module][name] = value

                elif item == "impurity_enrichment":
                    for m in range(len(data["impurity_enrichment"].get_value)):
                        name = f"impurity_enrichment({str(m + 1).zfill(1)})"
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
                        except ValueError:
                            pass
                        else:
                            if len(split_line) > 1:
                                line_value = ", ".join(list(split_line))

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

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    out_file:
        Output file for new IN.DAT

    """
    filter_list = ["f_nd_impurity_electrons(", "zref(", "imp_rich", "vmec"]
    # Special parameters that require different formatting
    parameters = get_parameters(data)

    for module in parameters:
        # Write module heading: format to be more readable again
        formatted_module = module.replace("_", " ").title()
        write_title(f"{formatted_module}", out_file)

        # Write out parameters for this module
        for parameter, info in parameters[module].items():
            if any(var_name in parameter for var_name in filter_list):
                # No justification formatting if parameter is in filter list
                parameter_line = f"{parameter} = {info}\n"
            else:
                # All other parameters
                # Left justification set to 8 to allow easier reading
                # info can currently be either a value or a dict
                if (
                    type(info) is dict
                    and "value" in info
                    and (type(info.get("comment")) is str)
                ):
                    parameter_line = (
                        f"{parameter.ljust(8)} = {info['value']} * {info['comment']}\n"
                    )
                else:
                    parameter_line = f"{parameter.ljust(8)} = {info}\n"

            # Finally write the line
            out_file.write(parameter_line)


def add_iteration_variable(data, variable_number):
    """Function to add iteration variable to IN.DAT data dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    variable_number:
        Iteration variable number to add

    """

    # Check the variable number is not already in the iteration variable list
    if variable_number not in data["ixc"].value:
        data["ixc"].value.append(variable_number)
        data["ixc"].value.sort()

    else:
        print(f"Variable number {variable_number} already in iteration variable list")


def remove_iteration_variable(data, variable_number):
    """Function to remove iteration variable from the IN.DAT data dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    variable_number:
        Iteration variable number to remove

    """

    # Check the variable is in the iteration variable list
    if variable_number in data["ixc"].value:
        data["ixc"].value.remove(variable_number)
        data["ixc"].value.sort()
    else:
        print(f"Variable number {variable_number} not in iteration variable list")


def add_constraint_equation(data, equation_number):
    """Function to add constraint equation to the IN.DAT data dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    equation_number:
        Constraint equation number to add

    """

    # Check the constraint is not already in the constraint equation list
    if equation_number not in data["icc"].value:
        data["icc"].value.append(equation_number)
        data["icc"].value.sort()

    else:
        print(f"Equation number {equation_number} already in constraint equations list")


def remove_constraint_equation(data, equation_number):
    """Function to remove a constraint equation from the IN.DAT data
    dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    equation_number:
        Constraint equation number to remove

    """

    # Check the constraint is in
    if equation_number in data["icc"].value:
        data["icc"].value.remove(equation_number)
        data["icc"].value.sort()

    else:
        print(f"Equation number {equation_number} not in constraint equations list")


def add_parameter(data, parameter_name, parameter_value):
    """Function to add/change parameter to the IN.DAT data dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    parameter_name:
        Name of the parameter to add
    parameter_value:
        Value of the parameter to add

    """
    # Load dicts from dicts JSON file
    dicts = get_dicts()

    # Check that the parameter is not already in the dictionary
    if parameter_name not in data:
        parameter_group = find_parameter_group(parameter_name)
        if "f_nd_impurity_electrons" in parameter_name:
            comment = dicts["DICT_DESCRIPTIONS"]["f_nd_impurity_electrons"]
        else:
            try:
                comment = dicts["DICT_DESCRIPTIONS"][parameter_name]
            except KeyError:
                # The dictionary doesn't recognise the variable name
                print(
                    f"Warning: Description for {parameter_name}",
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

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    parameter_name:
        Name of the parameter to remove

    """

    # Check that the parameter exists in the data dictionary
    if parameter_name in data:
        del data[parameter_name]

    # Inform the user that the parameter requested for deletion isn;t in the
    # data dictionary
    else:
        print(f"Parameter {parameter_name} not in IN.DAT")


def change_array(data, name, array_id, array_val):
    """Function to change value in generic array

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    name:
        generic array name
    array_id:
        generic array index
    array_val:
        generic array value
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

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    bound:
        Bound number associated with iteration variable number to
        change
    bound_type:
        States whether bound is upper of lower bound
    bound_value:
        New value of the bound

    """

    # Put bound type into lower cases for consistency
    bound_type = bound_type.lower()

    # if the bound is not in the bounds dictionary initialise an empty
    # dictionary and assign new bound
    if bound not in data["bounds"].value:
        data["bounds"].value[bound] = {}
        data["bounds"].value[bound][bound_type] = str(bound_value)

    # If bound already exists change value
    elif bound in data["bounds"].value:
        data["bounds"].value[bound][bound_type] = str(bound_value)

    # Bound not recognised.
    else:
        print(f"Bound {bound} not recognised. Check type == string")


def remove_bound(data, bound, bound_type):
    """Function to remove a bound from the bounds entry in the IN.DAT data
    dictionary

    Parameters
    ----------
    data:
        Data dictionary for the IN.DAT information
    bound:
        Bound number associated with iteration variable number to
        change
    bound_type:
        States whether bound is upper or lower bound

    """

    # use local variable for cleanliness
    bounds = data["bounds"].value

    # If the bound exists (and is of the correct type) in the bounds dictionary
    if bound in bounds and bound_type in bounds[bound]:
        del bounds[bound][bound_type]

        # if the bound number is now an empty dictionary delete it also
        if len(bounds[bound].keys()) == 0:
            del bounds[bound]


def fortran_float_to_py(f: str) -> str:
    if not isinstance(f, str):
        return f
    return sub(r"([0-9]+\.[0-9]+)(?:D|d)([0-9]+)", r"\1e\2", f)


def parameter_type(name, value):
    """Function to return value in correct format for altering values etc.

    Parameters
    ----------
    name:
        Name of parameter to check type
    value:
        Value of parameter to format

    Returns
    -------
    :
        Formatted value
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
        if "int_array" in param_type:
            return [item if item is None else int(item) for item in value]
            # Convert list to ints, but not if the value is None

        # otherwise, return value
        return value

    # Check if parameter is a string
    if isinstance(value, str):
        # If a real variable just convert to float
        if "real_variable" in param_type:
            # Prepare so float conversion succeeds
            value = value.lower()
            value = value.replace("d", "e")
            return float(value)

        # If a real array split and make a float list
        if "real_array" in param_type:
            # Prepare so float conversion succeeds
            value = value.lower()
            value = value.replace("d", "e")
            value = value.split(",")
            if value[-1] == "":
                value = value[:-1]
            return [float(item) for item in value]

        # If an integer variable convert to integer
        if "int_variable" in param_type:
            return int(value)

        # If an integer array split and make an integer list
        if "int_array" in param_type:
            value = value.split(",")
            if value[-1] == "":
                value = value[:-1]
            return [int(item) for item in value]

        # If type unknown return original value
        return value

    # If type is other return original value
    return value


def variable_constraint_type_check(item_number, var_type):
    """Function to put input into correct format for altering values etc.

    Parameters
    ----------
    item_number:
        Number associated with variable or constraint
    var_type:
        States whether item is iteration variable or constraint
        equation

    Returns
    -------
    :
        Formatted item_number
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
            if item_number.is_integer():
                return int(item_number)

            # rounded float number with warning
            print(
                f"Value {item_number} for {var_type} not an integer. Value rounded to {int(item_number)}."
                " Check!"
            )
            return int(item_number)

        except ValueError:
            print(
                f"Value {item_number} for {var_type} not valid. Check value!",
                file=stderr,
            )

    # Check if item is in float format
    elif isinstance(item_number, float):
        # If integer convert to float and return
        if item_number.is_integer():
            return int(item_number)

        # If not an integer warn of rounding and return rounded integer
        print(
            f"Value {item_number} for {var_type} not an integer. Value rounded to {int(item_number)}. Check!"
        )
        return int(item_number)

    # If already an integer return unchanged
    elif isinstance(item_number, int):
        return item_number

    # Value not recognised
    else:
        raise ValueError(
            f"Value {item_number} for {var_type} not a recognised format. Check value!"
        )


def variable_bound_check(bound_number, bound_type):
    """Function to put bound_number and bound_type into correct format

    Parameters
    ----------
    bound_number:
        Bound number to check
    bound_type:
        States whether bound is upper or lower bound

    Returns
    -------
    :
         Formatted bound number and bound type
    """

    # put bound type into lower case for consistency
    bound_type = bound_type.lower()

    # check if bound is one of the allowed values if not warn user
    if bound_type not in ["l", "u", "upper", "lower"]:
        print(
            f"Bound type '{bound_type}' not recognised. Must be one of "
            "['u', 'l', 'U', 'L', 'lower', 'upper', 'LOWER', 'UPPER']"
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
    if isinstance(bound_number, int):
        return str(bound_number), bound_type

    # If a float convert to str but warn of rounding when changing from float
    # to int
    if isinstance(bound_number, float):
        if bound_number.is_integer():
            return int(bound_number), bound_type
        bound_number = int(bound_number)
        print(
            f"Bound number {bound_number} not an integer. Value rounded to {int(bound_number)}"
        )
        return bound_number, bound_type

    raise TypeError(f"Unsupported type for bound_number: {type(bound_number)}")


class INVariable:
    """Class to stores the information of a single variable from the
        IN.DAT file

    Parameters
    ----------
    name:
        Item name
    value:
        Item value
    v_type:
        Type of item
    parameter_group:
        PROCESS variable group item belongs to
    comment:
        Comment for item
    """

    def __init__(self, name, value, v_type, parameter_group, comment):
        self.name = name
        self.value = value
        self.v_type = v_type
        self.parameter_group = parameter_group
        self.comment = comment

    def __eq__(self, value):
        # intentionally missing .comment, this is not necessary for the variables to be equal
        return (
            self.name == value.name
            and self.value == value.value
            and self.v_type == value.v_type
            and self.parameter_group == value.parameter_group
        )

    def __repr__(self):
        return (
            f"{type(self).__name__}(name={self.name!r}, value={self.value!r}, v_type={self.v_type!r}, "
            f"parameter_group={self.parameter_group!r}, comment={self.comment!r})"
        )

    @property
    def get_value(self):
        """ """
        if self.v_type != "Bound":
            return parameter_type(self.name, self.value)
        return self.value


class InDat:
    """Class 'InDat' for handling IN.DAT data. It handles

    - Reading IN.DAT files
    - Writing IN.DAT files
    - Storing information in dictionary for use in other codes
    - Alterations to IN.DAT
    """

    def __init__(self, filename="IN.DAT", start_line=0):
        """Initialise class

        Parameters
        ----------
        filename:
            Name of input IN.DAT
        start_line:
            Line to start reading from
        """

        self.filename = filename
        self.start_line = start_line

        # Initialise parameters
        self.in_dat_lines = []
        self.data = {}
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
            l_line = line.lower() if "vmec" not in line.split("=")[0].lower() else line

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
                        f"that line in IN.DAT is valid. Line skipped!\n{line}",
                        file=stderr,
                    )

                    # Store the first part of the unrecognised line (probably a
                    # variable name) as an unrecognised var
                    unrecognised_var = line.split("=")[0].strip()
                    self.unrecognised_vars.append(unrecognised_var)

    def process_line(self, line_type, line):
        """Function to process the line and return the appropriate INVariable
        object

        Parameters
        ----------
        line_type:
            Type of information the line contains
        line:
            Line from IN.DAT to process

        """

        # Load dicts from dicts JSON file
        dicts = get_dicts()

        # Create bound variable class using INVariable class if the bounds entry
        # doesn't exist
        if "bounds" not in self.data:
            self.data["bounds"] = INVariable("bounds", {}, "Bound", "Bound", "Bounds")

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

            if array_name not in self.data:
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

        Parameters
        ----------
        line:
            Line from IN.DAT to process

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

        Parameters
        ----------
        line:
            Line from IN.DAT to process

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
        if "icc" not in self.data:
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
                    self.add_duplicate_variable(f"icc = {item}")
            # Don't sort the constraints! Preserves what's eq, what's ineq;
            # first neqns are eqs, rest are ineqs
            # self.data["icc"].value.sort()

    def process_iteration_variables(self, line):
        """Function to process iteration variables entry in IN.DAT

        Parameters
        ----------
        line:
            Line from IN.DAT to process

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
        if "ixc" not in self.data:
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
                    self.add_duplicate_variable(f"ixc = {item}")
            self.data["ixc"].value.sort()

    def process_bound(self, line):
        """Function to process bound entries in IN.DAT

        Parameters
        ----------
        line:
            Line from IN.DAT to process

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
        if bound not in self.data["bounds"].value:
            self.data["bounds"].value[bound] = {}
        elif self.data["bounds"].value[bound].get(bound_type):
            # Duplicate bound
            self.add_duplicate_variable(f"bound{bound_type}({bound})")

        # Populate self.data dictionary with bound information
        self.data["bounds"].value[bound][bound_type] = bound_value

    def process_array(self, line, empty_array):
        """Function to process generic array

        Parameters
        ----------
        line:
            Line from IN.DAT to process
        empty_array:
            Default array for this array name

        Returns
        -------
        type
            nothing
        """
        # TODO This is a mess after the Fortran module variable
        # re-initialisation work. It is hard to see how this can be improved
        # as the regex method of Fortran inspection (i.e. Python-Fortran
        # dictionaries) is increasingly untenable. An attempt is made here,
        # but with a view to the dictionaries method being dropped in future
        # in light of increasing Python f2py conversion.

        line_commentless = line.split("*")[0] if "*" in line else line

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
            self.add_duplicate_variable(f"{name}({fortran_index})")

        # Now we are sure that the Python list index exists, set its value
        self.data[name].value[index] = eval(fortran_python_scientific(value))

    def add_duplicate_variable(self, name):
        """Records duplicate variables in the input file.

        If a var is initialised more than once in the input file, the last
        value persists, but the overwriting is recorded here.

        Parameters
        ----------
        name:
            The name of the variable being duplicated
        """
        self.duplicates.append(name)

    def add_iteration_variable(self, variable_number):
        """Function to add iteration variable to IN.DAT data dictionary

        Parameters
        ----------
        variable_number:
            Iteration variable number to add

        """

        # format iteration variable number
        variable_number = variable_constraint_type_check(
            variable_number, "iteration variable"
        )
        # add iteration variable to IN.DAT data dictionary
        add_iteration_variable(self.data, variable_number)

    def remove_iteration_variable(self, variable_number):
        """Function to remove iteration variable to IN.DAT data dictionary

        Parameters
        ----------
        variable_number:
            Iteration variable number to remove

        """

        # format iteration variable number
        variable_number = variable_constraint_type_check(
            variable_number, "iteration variable"
        )
        # remove iteration variable from IN.DAT data dictionary
        remove_iteration_variable(self.data, variable_number)

    def add_constraint_equation(self, equation_number):
        """Function to add constraint equation to IN.DAT data dictionary

        Parameters
        ----------
        equation_number:
            Constraint equation number to add

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

        Parameters
        ----------
        equation_number:
            Constraint equation number to remove

        """

        # format constraint equation number
        equation_number = variable_constraint_type_check(
            equation_number, "constraint equation"
        )

        # remove constraint equation from IN.DAT data dictionary
        remove_constraint_equation(self.data, equation_number)

    def add_parameter(self, parameter_name, parameter_value):
        """Function to add/change parameter to IN.DAT data dictionary

        Parameters
        ----------
        parameter_name:
            Name of parameter to add/change
        parameter_value:
            Value of parameter to add/change

        """

        # add/change parameter to/in IN.DAT data dictionary
        add_parameter(self.data, parameter_name, parameter_value)

    def change_array(self, array_name, array_index, array_value):
        """Function to change value in generic array

        Parameters
        ----------
        array_name:
            name of generic array to change
        array_index:
            index of generic array to change
        array_value:
            value to change array entry to
        """

        # change generic array value in IN.DAT data dictionary
        change_array(self.data, array_name, array_index, array_value)

    def remove_parameter(self, parameter_name):
        """Function to remove parameter from IN.DAT data dictionary

        Parameters
        ----------
        parameter_name:
            Name of parameter to remove

        """

        # remove parameter from IN.DAT data dictionary
        remove_parameter(self.data, parameter_name)

    def add_bound(self, bound, bound_type, bound_value):
        """Function to add/change a bound in IN.DAT data dictionary

        Parameters
        ----------
        bound:
            Bound number to add/change
        bound_type:
            States whether bound is upper or lower bound
        bound_value:
            Value of bound to add/change

        """

        # format bound number and bound type
        bound, bound_type = variable_bound_check(bound, bound_type)

        # add/change bound to/in IN.DAT data dictionary
        add_bound(self.data, bound, bound_type, bound_value)

    def remove_bound(self, bound, bound_type):
        """Function to remove a bound from IN.DAT data dictionary

        Parameters
        ----------
        bound:
            Bound number to remove
        bound_type:
            States whether bound is upper or lower bound

        """

        # format bound number and bound type
        bound, bound_type = variable_bound_check(bound, bound_type)

        # remove bound from IN.DAT data dictionary
        remove_bound(self.data, bound, bound_type)

    def write_in_dat(self, output_filename="new_IN.DAT"):
        """Function to write data to output file called 'output_filename'

        Parameters
        ----------
        output_filename:
             (Default value = "new_IN.DAT")
        """

        # create and open output file
        with open(output_filename, "w") as output:
            # Write Header
            write_title("", output)

            # Write Constraint Equations
            write_constraint_equations(self.data, output)

            # Write Iteration Variables
            write_iteration_variables(self.data, output)

            # Write parameters
            write_parameters(self.data, output)

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
    self.data["parameters"]["physics_variables"]["i_plasma_geometry"]["value"] = 0
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

        Parameters
        ----------
        var_name:
            The name of the parameter to be returned

        Returns
        -------
        dict
            Dictionary of parameter's information, or None if not found
        """
        modules = self.data["parameters"]

        for module_dict in modules.values():
            var_dict = module_dict.get(var_name)
            if var_dict is not None:
                break

        return var_dict

    def is_param_defined(self, var_name):
        """Check if a parameter is defined or not in the input data.

        Parameters
        ----------
        var_name:
            Name of the parameter to be checked

        Returns
        -------
        bool
            True if defined, False if not
        """
        defined = False

        if self.get_param(var_name):
            defined = True

        return defined

    def get_param_value(self, var_name):
        """Gets the value of a parameter from the input data.

        Parameters
        ----------
        var_name:
            The name of the parameter

        Returns
        -------
        int, float
            Value of that parameter; can be any type. None if not found.
        """
        value = None
        var_dict = self.get_param(var_name)
        if var_dict:
            value = var_dict.get("value")

        return value


if __name__ == "__main__":
    i = InDat(filename="IN.DAT")
    i.write_in_dat()
