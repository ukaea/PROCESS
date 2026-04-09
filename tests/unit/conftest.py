"""pytest config file for unit tests.

Define fixtures that will be shared across unit test modules.
"""

from pathlib import Path

import pytest

from process.core.init import init_all_module_vars


@pytest.fixture(scope="function", autouse=True)
def reinit_fix():
    """Re-initialise the data structure before each test is run.

    This is run once before each unit test (function scope),
    ensuring that all of the module variables are set to their initial values.
    'autouse' ensures that this fixture is used automatically by any test
    function in the unit directory.
    """
    init_all_module_vars()


@pytest.fixture()
def input_file():
    """Input file for testing.

    :return: path to input file
    :rtype: str
    """
    data_path = Path(__file__).parent / "data"
    input_file = data_path / "large_tokamak_IN.DAT"
    # Convert input file path to absolute and string
    return str(Path(input_file).resolve())
