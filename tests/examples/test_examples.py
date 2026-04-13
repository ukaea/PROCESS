"""Integration test for example notebooks in examples/ dir"""

import os
from pathlib import Path
from shutil import copytree, ignore_patterns

import jupytext
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

    # This change of directory is undone by the return_to_root fixture, hence we do not need to change back directories here
    os.chdir(tmp_path / "examples")

    # Return tmp_path/examples, now containing files copied from examples dir
    return tmp_path / "examples"


def _get_location(loc, name):
    name = Path(name).stem + "{}"
    notebook = jupytext.read(loc / name.format(".ex.py"))
    jupytext.write(notebook, loc / name.format(".ex.ipynb"), fmt="ipynb")
    return loc / name.format(".ex.ipynb")


def test_scan(examples_temp_data):
    """Run scan.ex.py notebook check no exceptions are raised and that an MFILE is created.

    scan.ex.py intentionally produces files when running the notebook, but remove
    them when testing.
    :param examples_temp_data: temporary dir containing examples files
    :type examples_temp_data: Path
    """
    scan_notebook_location = _get_location(examples_temp_data, "scan")

    with (
        testbook(scan_notebook_location, execute=False, timeout=1200) as tb,
        tb.patch(
            "process.core.repository._PROCESS_ROOT",
            new=scan_notebook_location.parent.resolve().as_posix(),
        ),
    ):
        tb.execute()
        # Run entire scan.ex.py notebook and assert an MFILE is created
        assert Path(examples_temp_data / "data/scan_example_file_MFILE.DAT").exists()


@pytest.mark.parametrize(
    "name",
    (
        "introduction",
        "single_model_evaluation",
        "vary_run_example",
        "optimum_solutions_comparison",
    ),
)
def test_no_assertion_solutions(name, examples_temp_data):
    """Run examples and check no exceptions are raised.

    :param examples_temp_data: temporary dir containing examples files
    """
    notebook_location = _get_location(examples_temp_data, name)
    with (
        testbook(notebook_location, execute=False, timeout=600) as tb,
        tb.patch(
            "process.core.repository._PROCESS_ROOT",
            new=notebook_location.parent.resolve().as_posix(),
        ),
    ):
        tb.execute()
