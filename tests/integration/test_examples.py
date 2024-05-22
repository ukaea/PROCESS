"""Integration test for example notebooks in examples/ dir"""

import os
from pathlib import Path
from shutil import copy, copytree
import pytest
import pandas
import numpy as np
from testbook import testbook


@pytest.fixture
def examples_temp_data(tmp_path):
    """Copy examples dir contents into temp dir for testing.

    Any changes are discarded on fixture teardown.
    :param tmp_path: temporary path fixture
    :type tmp_path: Path
    :return: temporary path containing examples files
    :rtype: Path
    """
    data_path = Path(__file__).parent.parent.parent / "examples"
    copytree(data_path, tmp_path / "examples")
    csv_json_path = (
        Path(__file__).parent.parent.parent / "process/io/mfile_to_csv_vars.json"
    )
    copy(csv_json_path, tmp_path)

    # Return tmp_path/examples, now containing files copied from examples dir
    return tmp_path / "examples"


def test_examples(examples_temp_data):
    """Run the examples.ipynb and check no exceptions are raised.

    examples.ipynb uses temp dirs to clean up any produced files itself.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    with testbook("examples.ipynb", execute=True, timeout=600):
        pass


@pytest.fixture
def scan_cleanup(examples_as_cwd):
    """Delete any files produced by scan.ipynb.

    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: None
    """
    yield

    # Teardown: delete produced files
    scan_files = Path.cwd().glob("*scan_*")
    for file in scan_files:
        if "IN.DAT" not in file.name:
            os.remove(file)


def test_scan(scan_cleanup):
    """Run scan.ipynb notebook check no exceptions are raised and that an MFILE is created.

    scan.ipynb intentionally produces files when running the notebook, but remove
    them when testing.
    :param scan_cleanup: fixture to delete any produced files
    :type scan_cleanup: None
    """
    with testbook("scan.ipynb", execute=True, timeout=600):
        # Run entire scan.ipynb notebook and assert an MFILE is created
        assert os.path.exists("a_scan_input_file_MFILE.DAT")


@pytest.fixture
def csv_cleanup(examples_as_cwd):
    """Delete any files produced by csv_output.ipynb.

    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: None
    """
    yield

    # Teardown: delete produced files
    csv_files = Path.cwd().glob("*csv_output_*")
    for file in csv_files:
        if "_MFILE.DAT" not in file.name:
            os.remove(file)


def test_csv(csv_cleanup):
    """Run csv_output.ipynb, check no exceptions are raised, check a csv file exists and check the csv file contains data.

    csv_output.ipynb intentionally produces files when running the notebook, but remove
    them when testing.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    with testbook("csv_output.ipynb", execute=True, timeout=600):
        # Check csv file is created
        assert os.path.exists("csv_output_large_tokamak_MFILE.csv")

        # Read in the csv file created by test and check it contains positive floats
        readcsv = pandas.read_csv("csv_output_large_tokamak_MFILE.csv")
        values = readcsv["Value"]
        value_array = np.array(values)
        check_float = False
        check_positive = False
        value_array_type = value_array.dtype
        if value_array_type.kind == "f":
            check_float = True
        assert check_float

        check_positive_count = np.sum(value_array > 0)
        if check_positive_count == len(value_array):
            check_positive = True
        assert check_positive


def test_plot_solutions(examples_temp_data):
    """Run the plot_solutions.ipynb and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    with testbook("plot_solutions.ipynb", execute=True, timeout=600):
        pass
