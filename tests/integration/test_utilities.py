"""Integration tests for utilities.

These tests check the utilities that PROCESS uses, mainly for file IO. They run
on each of the regression test scenarios.
"""

import pytest
import logging
import process.io.mfile as mf
import process.io.in_dat as indat
import process.io.plot_proc as pp

logger = logging.getLogger(__name__)

logger.info("Running utilities integration tests")

# TODO More utilities tests to be implemented
# test_convert_in_dat


@pytest.fixture
def mfile_path(temp_data, mfile_name):
    """Create a path to a scenario's MFile.

    :param temp_data: temporary path containing data files
    :type temp_data: Path
    :param mfile_name: name of the reference MFile
    :type mfile_name: str
    :return: Path to that scenario's MFile
    :rtype: Path
    """
    mfile_path = temp_data / mfile_name
    return mfile_path


@pytest.fixture
def input_file_path(temp_data):
    """Create a path to a scenario's input file.

    :param temp_data: temporary path containing data files
    :type temp_data: Path
    :return: Path to that scenario's IN.DAT
    :rtype: Path
    """
    return temp_data / "large_tokamak_IN.DAT"


def test_mfile_lib(mfile_path):
    """Test the PROCESS mfile library.

    :param mfile_path: Path to the scenario's MFile
    :type mfile_path: Path
    """
    logger.info("Testing mfile.py")

    # Test MFile for this scenario
    # This try/except is not necessary, but allows additional logging to be
    # added for clarity in addition to pytest's own logging
    try:
        assert mf.test(str(mfile_path)) is True
        # mf.test returns True on success
    except AssertionError:
        logger.exception(f"mfile test for {mfile_path.name} has failed")
        raise


def test_in_dat_lib(input_file_path):
    """Test the PROCESS in_dat library.

    :param input_file_path: Path to a scenario's input file
    :type input_file_path: Path
    """
    logger.info("Testing in_dat")

    # Test MFile for this scenario
    try:
        assert indat.test(str(input_file_path)) is True
    except AssertionError:
        logger.error(f"in_dat test for {input_file_path.name} has failed")
        raise


def test_plot_proc(mfile_path):
    """Test the PROCESS plot_proc script.

    :param mfile_path: Path to the scenario's MFile
    :type mfile_path: Path
    """
    logger.info("Testing plot_proc.py")

    # Test plot_proc on an MFile
    try:
        assert pp.test(str(mfile_path)) is True
    except AssertionError:
        logger.exception(f"plot_proc test for {mfile_path.name} has failed")
        raise
