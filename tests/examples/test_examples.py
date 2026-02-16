"""Integration test for example notebooks in examples/ dir"""

import os
from pathlib import Path
from shutil import copy, copytree, ignore_patterns

import jupytext
import numpy as np
import pandas as pd
import pytest
from testbook import testbook


@pytest.fixture(scope="module")
def examples_temp_data(tmp_path_factory):
    """Copy examples dir contents into temp dir for testing.

    Any changes are discarded on fixture teardown.
    :param tmp_path: temporary path fixture
    :type tmp_path: Path
    :return: temporary path containing examples files
    :rtype: Path
    """
    data_path = Path(__file__).parent.parent.parent / "examples"
    tmp_path = tmp_path_factory.mktemp("examples")
    copytree(
        data_path,
        tmp_path / "examples",
        ignore=ignore_patterns("*.md", "*log", "__pycache__", "*.ipynb*"),
    )
    csv_json_path = (
        Path(__file__).parent.parent.parent / "process/io/mfile_to_csv_vars.json"
    )
    copy(csv_json_path, tmp_path)

    # This change of directory is undone by the return_to_root fixture, hence we do not need to change back directories here
    os.chdir(tmp_path / "examples")

    # Return tmp_path/examples, now containing files copied from examples dir
    return tmp_path / "examples"


def _get_location(loc, name):
    name = Path(name).stem + "{}"
    notebook = jupytext.read(loc / name.format(".ex.py"))
    jupytext.write(notebook, loc / name.format(".ex.ipynb"), fmt="ipynb")
    return loc / name.format(".ex.ipynb")


def test_introductory_examples(examples_temp_data):
    """Run the introduction.ex.py and check no exceptions are raised.

    introudction.ex.py uses temp dirs to clean up any produced files itself.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path

    """
    example_notebook_location = _get_location(examples_temp_data, "introduction")
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
    """Run scan.ex.py notebook check no exceptions are raised and that an MFILE is created.

    scan.ex.py intentionally produces files when running the notebook, but remove
    them when testing.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    scan_notebook_location = _get_location(examples_temp_data, "scan")
    with testbook(scan_notebook_location, execute=True, timeout=1200):
        # Run entire scan.ex.py notebook and assert an MFILE is created
        assert os.path.exists(examples_temp_data / "data/scan_example_file_MFILE.DAT")


@pytest.mark.parametrize(
    "name", ("plot_solutions", "single_model_evaluation", "vary_run_example")
)
def test_no_assertion_solutions(name, examples_temp_data):
    """Run examples and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
    """
    plot_solutions_notebook_location = _get_location(examples_temp_data, name)
    with testbook(plot_solutions_notebook_location, execute=True, timeout=600):
        pass
