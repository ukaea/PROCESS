import json
import shutil
import tempfile
from pathlib import Path

import pytest

from process.core.io.mfile import MFile
from process.core.io.mfile_utils import get_mfile_initial_ixc_values


def test_get_mfile_initial_ixc_values(input_file, tmp_path):
    tmp_input_file = tmp_path / "IN.DAT"
    shutil.copy(input_file, tmp_input_file)

    iteration_variable_names, iteration_variable_values = get_mfile_initial_ixc_values(
        Path(tmp_input_file)
    )

    assert iteration_variable_names[0] == "b_plasma_toroidal_on_axis"
    assert iteration_variable_values[0] == 5.7

    assert iteration_variable_names[1] == "rmajor"
    assert iteration_variable_values[1] == 8.0

    assert iteration_variable_names[-1] == "dr_tf_wp_with_insulation"
    assert iteration_variable_values[-1] == 0.5

    # A default not provided in the MFile
    assert iteration_variable_names[-4] == "f_nd_alpha_electron"
    assert iteration_variable_values[-4] == 0.1


@pytest.fixture(scope="module")
def read_mfile():
    """Read-in MFILE for testing.

    :return: parsed mfile
    :rtype: mfile2dict.MFILEParser
    """
    data_path = Path(__file__).parent / "data"

    return MFile(data_path / "large_tokamak_MFILE.DAT")


@pytest.fixture(scope="module")
def temporary_dir():
    return tempfile.mkdtemp()


def test_write_json(read_mfile, temporary_dir):
    json_f = Path(temporary_dir, "2017_baseline.json")
    read_mfile.to_json(json_f)
    assert json_f.is_file()
    with open(json_f) as file:
        assert json.load(file)
