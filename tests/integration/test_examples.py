"""Integration test for examples.py.

examples.py is created by exporting examples.ipynb as a Python script.
"""
import runpy
import os
from pathlib import Path
import pytest

# TODO How to reliably run examples/examples.py script? Relies on relative path
# from project root dir and pytest having project root dir as cwd. Could this be
# improved?


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
def delete_plot_procs():
    yield
    plot_proc_1 = Path("../examples/plot_proc_1")
    plot_proc_2 = Path("../examples/plot_proc_2")
    plot_proc_1.unlink(missing_ok=True)
    plot_proc_2.unlink(missing_ok=True)


def test_examples(examples_as_cwd, delete_plot_procs):
    """Run the examples.py script and check no exceptions are raised.

    examples.py uses temp dirs to clean up any produced files itself.
    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: None
    """
    # runpy used to run entire examples.py script
    runpy.run_path("examples.py")


@pytest.fixture
def scan_cleanup(examples_as_cwd):
    """Delete any files produced by scan.py.

    :param examples_as_cwd: fixture to set examples dir as cwd
    :type examples_as_cwd: None
    """
    yield

    # Teardown: delete produced files
    scan_files = Path.cwd().glob("*_scan_*")
    for file in scan_files:
        if "IN.DAT" not in file.name:
            os.remove(file)


def test_scan(scan_cleanup):
    """Run the scan.py script and check no exceptions are raised.

    scan.py intentionally produces files when running the notebook, but remove
    them when testing.
    :param scan_cleanup: fixture to delete any produced files
    :type scan_cleanup: None
    """
    # Run entire scan.py script
    runpy.run_path("scan.py")
