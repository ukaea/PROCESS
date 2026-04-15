"""pytest config file for unit tests.

Define fixtures that will be shared across unit test modules.
"""

import os
from pathlib import Path
from shutil import copy

import pytest

from process.core.init import init_all_module_vars


@pytest.fixture(autouse=True)
def reinit_fix():
    """Re-initialise the data structure before each test is run.

    This is run once before each unit test (function scope),
    ensuring that all of the module variables are set to their initial values.
    'autouse' ensures that this fixture is used automatically by any test
    function in the unit directory.
    """
    init_all_module_vars()


@pytest.fixture
def input_file():
    """Input file for testing.

    :return: path to input file
    :rtype: str
    """
    data_path = Path(__file__).parent / "data"
    input_file = data_path / "large_tokamak_IN.DAT"
    # Convert input file path to absolute and string
    return str(Path(input_file).resolve())


@pytest.fixture
def temp_data(tmp_path: Path) -> Path:
    """Copy data dir contents into temp dir for testing.

    Any changes are discarded on fixture teardown.

    Parameters
    ----------
    tmp_path:
        temporary path fixture

    Returns
    -------
    :
        temporary path containing data files
    """
    data_path = Path(__file__).parent.parent / "integration" / "data"

    for data_file in data_path.glob("*"):
        dst = tmp_path / data_file.name
        copy(data_file, dst)

    # Return tmp_path, now containing files copied from data dir
    return tmp_path


@pytest.fixture
def temp_data_cwd(temp_data: Path):
    """Change cwd to temp_data dir, then yield it.

    Used when testing command-line args that look for files in the cwd.

    Parameters
    ----------
    temp_data:
        temporary path containing data files

    Yields
    ------
    :
        temporary path containing data files
    """
    # Setup by changing cwd to temp_data and yielding it
    old_wd = Path.cwd()
    os.chdir(temp_data)
    yield temp_data

    # Teardown by changing back to previous dir
    os.chdir(old_wd)
