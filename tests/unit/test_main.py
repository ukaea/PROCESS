"""Unit tests for the main.py module."""

import shutil
from pathlib import Path

import pytest

from process import data_structure
from process.main import SingleRun, VaryRun


def mock_init(*args, **kwargs):
    """Used to mock out __init__ methods on classes.

    :return: Nothing
    :rtype: Nonetype
    """
    return


@pytest.fixture
def single_run(monkeypatch, input_file, tmp_path):
    """Fixture for a SingleRun object.

    :param monkeypatch: monkeypath fixture
    :type monkeypatch: object
    :return: SingleRun object
    :rtype: SingleRun
    """
    monkeypatch.setattr(SingleRun, "__init__", mock_init)
    single_run = SingleRun()

    temp_input_file = shutil.copy(input_file, tmp_path / Path(input_file).name)

    single_run.input_file = str(temp_input_file)
    single_run.models = None
    single_run.set_filenames()
    single_run.initialise()
    return single_run


def test_single_run(single_run):
    """Assert SingleRun objects can be created.

    :param single_run: single_run fixture
    :type single_run: SingleRun
    """
    assert type(single_run) is SingleRun


def test_set_input(single_run, monkeypatch, input_file):
    """Check the input file validation and setting of path.

    :param single_run: single_run fixture
    :type single_run: SingleRun
    :param monkeypatch: monkeypatch fixture
    :type monkeypatch: object
    :param input_file: fixture for input file
    :type input_file: str
    """
    expected = input_file
    # Mock the input file path to isolate this test from the other Process
    # methods (don't have to run Process.parse_args() first to set up this way)
    monkeypatch.setattr(single_run, "input_file", input_file, raising=False)

    # Mocking undo trys to set the value as none

    # Mocks set up, can now run set_input()
    single_run.set_input()
    # Check path has been set
    assert data_structure.global_variables.fileprefix == expected


def test_set_output(single_run, monkeypatch):
    """Check output filename set correctly.

    :param single_run: single_run fixture
    :type single_run: SingleRun
    :param monkeypatch: monkeypatch fixture
    :type monkeypatch: object
    """
    # Expected output prefix
    expected = "output_prefix"
    # Mock self.filename_prefix on single_run with the value of expected
    monkeypatch.setattr(single_run, "filename_prefix", expected, raising=False)

    # Mocking undo trys to set the value as none
    # monkeypatch.setattr(data_structure.global_variables, "output_prefix", None)
    # Run the method, and extract the value
    single_run.set_output()

    assert data_structure.global_variables.output_prefix == expected


def test_initialise(single_run, monkeypatch):
    """Test that the init_module runs without crashing

    :param single_run: single_run fixture
    :type single_run: SingleRun
    """
    # Mock the init subroutine with a lambda function
    # Run initialise method; this will fail on a raised exception
    single_run.initialise()


def test_set_mfile(single_run, monkeypatch):
    """Check the mfile filename is being stored correctly.

    :param single_run: single_run fixture
    :type single_run: SingleRun
    :param monkeypatch: monkeypatch fixture
    :type monkeypatch: object
    """
    prefix = "test"
    expected = Path(prefix + "MFILE.DAT")
    # Mock filename_prefix and run
    monkeypatch.setattr(single_run, "filename_prefix", prefix, raising=False)
    single_run.set_mfile()
    assert single_run.mfile_path == expected


def test_finish(single_run, monkeypatch):
    """Check that the finish subroutine is called.

    :param single_run: single_run fixture
    :type single_run: SingleRun
    :param monkeypatch: monkeypatch fixture
    :type monkeypatch: object
    """
    single_run.finish()


@pytest.fixture
def vary_run(monkeypatch):
    """Fixture to return a VaryRun object.

    :param monkeypatch: monkeypatch fixture
    :type monkeypatch: object
    :return: VaryRun object
    :rtype: VaryRun
    """
    monkeypatch.setattr(VaryRun, "__init__", mock_init)
    return VaryRun()


def test_vary_run(vary_run):
    """Assert VaryRun object can be created.

    :param vary_run: vary_run fixture
    :type vary_run: VaryRun
    """
    assert type(vary_run) is VaryRun
