"""Integration test for example notebooks in examples/ dir"""

import os
from pathlib import Path
from shutil import copy, copytree

import numpy as np
import pandas as pd
import pytest
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

    # This change of directory is undone by the return_to_root fixture, hence we do not need to change back directories here
    os.chdir(tmp_path / "examples")

    # Return tmp_path/examples, now containing files copied from examples dir
    return tmp_path / "examples"


def test_examples(examples_temp_data):
    """Run the examples.ipynb and check no exceptions are raised.

    examples.ipynb uses temp dirs to clean up any produced files itself.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path

    """
    example_notebook_location = examples_temp_data / "examples.ipynb"
    with testbook(example_notebook_location, execute=True, timeout=600):
        # Check csv file is created
        assert os.path.exists(examples_temp_data / "data/large_tokamak_1_MFILE.csv")

        # Read in the csv file created by test and check it contains positive floats
        readcsv = pd.read_csv(examples_temp_data / "data/large_tokamak_1_MFILE.csv")
        value_array = np.array(readcsv["Value"])
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


def test_scan(examples_temp_data):
    """Run scan.ipynb notebook check no exceptions are raised and that an MFILE is created.

    scan.ipynb intentionally produces files when running the notebook, but remove
    them when testing.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    scan_notebook_location = examples_temp_data / "scan.ipynb"
    with testbook(scan_notebook_location, execute=True, timeout=1200):
        # Run entire scan.ipynb notebook and assert an MFILE is created
        assert os.path.exists(examples_temp_data / "data/scan_example_file_MFILE.DAT")


def test_plot_solutions(examples_temp_data):
    """Run plot_solutions.ipynb and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
     :type examples_temp_data: Path
    """
    plot_solutions_notebook_location = examples_temp_data / "plot_solutions.ipynb"
    with testbook(plot_solutions_notebook_location, execute=True, timeout=600):
        pass


def test_single_model_evaluation(examples_temp_data):
    """Run single_model_evaluation.ipynb and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path

    """
    single_model_evaluation_notebook_location = (
        examples_temp_data / "single_model_evaluation.ipynb"
    )
    with testbook(single_model_evaluation_notebook_location, execute=True, timeout=600):
        pass


def test_varyrun_example(examples_temp_data):
    """Run VaryRun-example.ipynb and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
     :type examples_temp_data: Path
    """
    varyrun_example_notebook_location = examples_temp_data / "VaryRun-example.ipynb"
    with testbook(varyrun_example_notebook_location, execute=True, timeout=600):
        pass
