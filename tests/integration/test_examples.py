"""Integration test for example notebooks in examples/ dir"""

import os
from pathlib import Path
import pytest
import pandas
import numpy as np
from testbook import testbook


@pytest.fixture
def examples_as_cwd():
    """Change the cwd to the examples dir for the duration of the fixture.

    When running a Jupyter notebook, the cwd is set to the notebook's dir. The
    examples.ipynb notebook relies on relative paths to files in the repository,
    due to it being difficult to consistently get the location of a notebook
    from within the notebook itself, so absolute paths can't be used.

    When pytest is used to run the examples.py script (created directly from the
    notebook), the script uses the actual cwd instead, which is (usually) the
    project root dir. Therefore the examples.ipynb notebook and pytest-run
    examples.py script both rely on the cwd being the same, but without
    intervention it is different in each case.

    Hence the test needs to set the cwd to the notebook's dir before running so
    that the examples.py script uses the same cwd as the notebook.
    """
    # Set up by storing cwd, which is usually the project root dir, then
    # changing to examples/
    cwd = Path.cwd()
    os.chdir("examples")
    yield

    # Teardown by reverting cwd change
    os.chdir(cwd)


@pytest.fixture
def delete_plot_procs(examples_as_cwd):
    """Delete any plot_proc files produced by examples.ipynb.

    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: None
    """
    yield
    plot_proc_1 = Path("../examples/plot_proc_1")
    plot_proc_2 = Path("../examples/plot_proc_2")
    plot_proc_3 = Path("../examples/plot_proc_3")
    plot_proc_1.unlink(missing_ok=True)
    plot_proc_2.unlink(missing_ok=True)
    plot_proc_3.unlink(missing_ok=True)


def test_examples(delete_plot_procs):
    """Run examples.ipynb and check no exceptions are raised.

    examples.ipynb uses temp dirs to clean up any produced files itself.
    """
    with testbook("examples.ipynb", execute=True):
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
    with testbook("scan.ipynb", execute=True, timeout=120):
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
    :param csv_cleanup: fixture to delete any produced files
    :type csv_cleanup: None
    """
    with testbook("csv_output.ipynb", execute=True):
        # Check csv file is created
        print(os.getcwd())
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


def test_plot_solutions(examples_as_cwd):
    """Run plot_solutions.ipynb and check no exceptions are raised.

    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: NoneType
    """
    with testbook("plot_solutions.ipynb", execute=True):
        pass
