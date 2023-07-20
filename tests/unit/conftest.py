"""pytest config file for unit tests.

Define fixtures that will be shared across unit test modules.
"""
import pytest
from process import fortran
from pathlib import Path


@pytest.fixture(scope="module", autouse=True)
def reinit_fix():
    """Re-initialise Fortran module variables before each test module is run.

    This is run once before each module's unit tests are run (module scope),
    ensuring that all Fortran module variables are set to initial values. The
    individual test functions in each module then use mocking to avoid changing
    the module variable values. autouse ensures that this fixture is used
    automatically by any test function in the unit directory.
    """
    fortran.init_module.init_all_module_vars()


@pytest.fixture()
def input_file():
    """Input file for testing.

    :return: path to input file
    :rtype: str
    """
    data_path = Path(__file__).parent / "data"
    input_file = data_path / "large_tokamak_IN.DAT"
    # Convert input file path to absolute and string
    input_file_path = str(Path(input_file).resolve())
    return input_file_path
