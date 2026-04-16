"""Tests for InDat and INVariable classes"""

from pathlib import Path

import pytest

import process.core.io.in_dat as indat
from process.core.io.in_dat import InDat, INVariable


@pytest.mark.parametrize("value", ["1", "1.0", "1,2,3", "1.0, 2.0", "string"])
def test_invariable_equality(value):
    """A test to check equality between INVariable's"""
    name = "test"
    v_type = "Parameter"
    parameter_group = "test_group"

    v1 = INVariable(name, value, v_type, parameter_group, "")
    v2 = INVariable(name, value, v_type, parameter_group, "different comment")

    assert v1 == v2


def test_rewritten_indat_identical(temp_int_data_cwd):
    indat = InDat(filename=str(temp_int_data_cwd / "large_tokamak_IN.DAT"))
    indat.write_in_dat("new.IN.DAT")

    new_indat = InDat(filename=str(temp_int_data_cwd / "new.IN.DAT"))

    assert indat.data == new_indat.data


@pytest.fixture
def input_file_path(temp_int_data: Path):
    """Create a path to a scenario's input file using integration data.

    Parameters
    ----------
    temp_int_data:
        temporary path containing data files
    """
    return temp_int_data / "large_tokamak_IN.DAT"


def test_in_dat_lib(input_file_path, tmp_path):
    """Test the PROCESS in_dat library.

    :param input_file_path: Path to a scenario's input file
    :type input_file_path: Path
    """

    # Test MFile for this scenario
    outfile = Path(tmp_path, "test_out_IN.DAT")
    i = indat.InDat(filename=str(input_file_path))
    i.write_in_dat(output_filename=outfile.as_posix())

    assert outfile.is_file()
