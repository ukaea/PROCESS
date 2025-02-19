"""Handle parsing, validation, and actioning of a PROCESS input file (*IN.DAT)."""

import copy
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from warnings import warn

from process import fortran
from process.exceptions import ProcessValidationError, ProcessValueError
from process.utilities.f2py_string_patch import (
    f2py_compatible_to_string,
    string_to_f2py_compatible,
)

NumberType = int | float
ValidInputTypes = NumberType | str

DataTypes = (int, float, str)
"""Valid input variable types alias"""


@dataclass
class InputVariable:
    """A variable to be parsed from the input file."""

    module: Any
    """The Fortran module that this variable should be set on."""
    type: type
    """The expected type of the variable"""
    range: tuple[NumberType, NumberType] | None = None
    """The valid range of values (int and float only) inclusive of the endpoints."""
    choices: list[ValidInputTypes] | None = None
    """The valid values of this input."""
    array: bool = False
    """Is this input assigning values to an array?"""
    additional_validation: (
        Callable[[str, ValidInputTypes, int | None, "InputVariable"], ValidInputTypes]
        | None
    ) = None
    """A function that takes the input variable: name, value, array index, and config (this dataclass)
    as input and returns a cleaned version of the input value. May raise
    a ProcessValidationError.

    NOTE: The value passed in has already been cleaned in the default way and has
    been cast to the specified `type`.
    """
    additional_actions: (
        Callable[[str, ValidInputTypes, int | None, "InputVariable"], None] | None
    ) = None
    """A function that takes the input variable: name, value, array index, and config (this dataclass)
    as input and performs some additional action in addition to the default actions prescribed by the variables
    config. May raise a ProcessValidationError.

    NOTE: The value passed in has already been cleaned in the default
    way and has been cast to the specified `type`.

    NOTE: the input name is the cleaned name as in the input file NOT the `target_variable`.
    """
    target_name: str | None = None
    """Indicates the parsed variable name is different to its target on the `module`."""
    set_variable: bool = True
    """Do not attempt to set the variable on the module. Instead only do any `additional_actions`."""

    def __post_init__(self):
        if self.type not in DataTypes:
            error_msg = (
                f"Cannot parse variable of type {self.type}, must be one of {DataTypes}"
            )
            raise ProcessValueError(error_msg)

        if self.type not in (int, float) and self.range is not None:
            error_msg = "Can only apply a range to integer or float variables"
            raise ProcessValueError(error_msg)

        if not self.set_variable and self.additional_actions is None:
            warn(
                "Not setting variable or performing an additional action. Why are you parsing this variable?",
                stacklevel=2,
            )


INPUT_VARIABLES = {
    "runtitle": InputVariable(fortran.global_variables, str),
    "ioptimz": InputVariable(fortran.numerics, int, choices=[1, -2]),
    "epsvmc": InputVariable(fortran.numerics, float, range=(0.0, 1.0)),
    "boundl": InputVariable(fortran.numerics, float, array=True),
}


def parse_input_file():
    input_file = f2py_compatible_to_string(fortran.global_variables.fileprefix)

    input_file_path = Path("IN.DAT")
    if input_file != "":
        input_file_path = Path(input_file)

    with input_file_path.open("r") as f:
        lines = f.readlines()

    for line_no, line in enumerate(lines, start=1):
        stripped_line = line.strip()

        # don't bother trying to process blank lines
        # or comment lines
        if stripped_line == "" or stripped_line[0] == "*":
            continue

        # matches (variable name, array index, value)
        # NOTE: array index is Fortran-based hence starts at 1.
        line_match = re.match(
            r"([a-zA-Z0-9_]+)(?:\(([0-9]+)\))?[ ]*=[ ]*([ +\-a-zA-Z0-9.]+).*",
            stripped_line,
        )

        if line_match is None:
            error_msg = (
                f"Unable to parse line {line_no} of the input file ({stripped_line})"
            )
            raise ProcessValueError(error_msg)

        variable_name, array_index, variable_value = line_match.groups()

        variable_config = INPUT_VARIABLES.get(variable_name)

        if variable_config is None:
            error_msg = (
                f"Unrecognised input '{variable_name}' at line {line_no} of input file."
            )
            warn(error_msg, stacklevel=1)
            continue
            # TODO: re-enable error instead
            # raise ProcessValidationError(error_msg)

        # validate the variable and also clean it (cast to correct type)

        clean_variable_value = validate_variable(
            variable_name, variable_value.strip(), array_index, variable_config, line_no
        )

        # check if the target name in the module is different to the variable name
        # in the input file
        variable_name_in_module = variable_config.target_name or variable_name

        if variable_config.set_variable:
            if variable_config.array:
                set_array_variable(
                    variable_name_in_module,
                    clean_variable_value,
                    int(array_index),
                    variable_config,
                )
            else:
                set_scalar_variable(
                    variable_name_in_module, clean_variable_value, variable_config
                )

        if variable_config.additional_actions is not None:
            # intentionally passing the variable name as in the input file,
            # not the module target name.
            variable_config.additional_actions(
                variable_name,
                clean_variable_value,
                None if array_index is None else int(array_index),
                variable_config,
            )


def validate_variable(
    name: str,
    value: str,
    array_index: int | None,
    config: InputVariable,
    line_number: int,
) -> ValidInputTypes:
    """Validate an input.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param array_index: the array index of the variable in the input file.
    :param config: the config of the variable that describes how to validate and process it.

    :returns: the input value with the correct type.
    """
    # check its an array if and only if the config says so
    if array_index is None and config.array:
        error_msg = f"Expected '{name}' at line {line_number} to be an array."
        raise ProcessValidationError(error_msg)

    if array_index is not None and not config.array:
        error_msg = f"Not expecting '{name}' at line {line_number} to be an array."
        raise ProcessValidationError(error_msg)

    if config.type in [float, int]:
        value = value.lower().replace("d", "e")

    try:
        clean_value = config.type(value)
    except ValueError as e:
        error_msg = f"Cannot cast variable name '{name}' at line {line_number} to a {config.type} (value = {value})"
        raise ProcessValidationError(error_msg) from e

    if (
        config.range is not None
        and not config.range[0] <= clean_value <= config.range[1]
    ):
        error_msg = (
            f"Variable '{name}' at line {line_number} is not on the prescribed range"
            f" {config.range} (value = {value})"
        )
        raise ProcessValidationError(error_msg)

    if config.choices is not None and clean_value not in config.choices:
        error_msg = f"Variable '{name}' at line {line_number} is not one of {config.choices} (value = {value})"
        raise ProcessValidationError(error_msg)

    if config.additional_validation is not None:
        clean_value = config.additional_validation(
            name, clean_value, int(array_index), config
        )

    return clean_value


def set_scalar_variable(name: str, value: ValidInputTypes, config: InputVariable):
    """Set a scalar (not part of an array) variable in the `config.module`.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param config: the config of the variable that describes how to validate and process it.
    """
    current_value = getattr(config.module, name, ...)

    # use ... sentinel because None is probably a valid return from Fortran
    # and definately will be when moving to a Python data structure
    if current_value is ...:
        error_msg = (
            f"Fortran module '{config.module}' does not have a variable '{name}'."
        )
        raise ProcessValueError(error_msg)

    if config.type is str:
        setattr(config.module, name, string_to_f2py_compatible(current_value, value))
    else:
        setattr(config.module, name, value)


def set_array_variable(name: str, value: str, array_index: int, config: InputVariable):
    """Set an array variable in the `config.module`.

    The way PROCESS input files are structured, each element of the array is provided on one line
    so this function just needs to set the `value` at `array_index` (-1) because of Fortran-based indexing.

    :param name: the name of the input.
    :param value: the value of the input variable.
    :param array_index: the array index of the variable in the input file.
    :param config: the config of the variable that describes how to validate and process it.
    """
    current_array = getattr(config.module, name, ...)

    # use ... sentinel for same reason as above
    if current_array is ...:
        error_msg = f"Fortran module '{config.module}' does not have an array '{name}'."
        raise ProcessValueError(error_msg)

    if config.type is str:
        raise NotImplementedError

    new_array = copy.deepcopy(current_array)
    new_array[array_index - 1] = value
    setattr(config.module, name, new_array)
