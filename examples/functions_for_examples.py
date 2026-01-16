"""Contains functions needed for running the examples files in a temporary directory
and for getting initial values from input files."""

import re
from pathlib import Path
from shutil import copy
from tempfile import TemporaryDirectory

import process.data_structure as data_structure
from process import iteration_variables
from process.main import SingleRun


def copy_to_temp_dir(input_rel, proj_dir):
    """Copy an input file to a new temp dir and return its new path.

    The new TemporaryDirectory object is returned to avoid destruction of the
    object, which results in deletion of the directory prematurely. This way
    the cleanup() method can be used to delete the directory when required.
    :param input_rel: file path relative to project root dir
    :type input_rel: str
    :return: temporary dir and absolute path to file in temp dir
    :rtype: (TemporaryDirectory, pathlib.Path)
    """
    # Create temporary dir to contain the run's outputs
    temp_dir = TemporaryDirectory()
    temp_dir_path = Path(temp_dir.name)

    # Define absolute path for input file
    input_rel_path = Path(input_rel)
    input_abs_path = proj_dir / input_rel_path

    try:
        assert input_abs_path.exists()
    except AssertionError as err:
        raise FileNotFoundError("Input file doesn't exist.") from err

    # Copy input file to temp dir
    copy(input_abs_path, temp_dir_path)
    temp_input_path = temp_dir_path / input_abs_path.name

    return temp_dir, temp_input_path, temp_dir_path


def get_initial_values(file_path):
    """
    Initialise the input file and obtain the initial values of the iteration variables

    :param file_path: Path to the input file
    :type file_path: pathlib.Path
    """
    SingleRun(str(file_path))
    iteration_variables.load_iteration_variables()

    iteration_variable_names = []
    iteration_variable_values = []

    for i in range(data_structure.numerics.nvar):
        ivar = data_structure.numerics.ixc[i].item()

        itv = iteration_variables.ITERATION_VARIABLES[ivar]

        iteration_variable_names.append(itv.name)
        if array := re.match(r"(\w+)\(([0-9]+)\)", itv.name):
            var_name = array.group(1)
            index = array.group(2)
            iteration_variable_values.append(
                getattr(itv.module, var_name)[int(index) - 1]
            )
        else:
            iteration_variable_values.append(getattr(itv.module, itv.name))

    return iteration_variable_names, iteration_variable_values
