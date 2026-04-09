import json
import os
import pickle
import tempfile
from pathlib import Path

import pytest
import yaml

from process.core.io import mfile2dict


@pytest.fixture(scope="module")
def read_mfile():
    """Read-in MFILE for testing.

    :return: parsed mfile
    :rtype: mfile2dict.MFILEParser
    """
    data_path = Path(__file__).parent / "data"
    mfile_path = data_path / "large_tokamak_MFILE.DAT"
    mfile_path_str = str(Path(mfile_path).resolve())

    return mfile2dict.MFILEParser(mfile_path_str)


@pytest.fixture(scope="module")
def temporary_dir():
    return tempfile.mkdtemp()


def test_parser_succeed(read_mfile):
    assert read_mfile._mfile_data  # noqa: SLF001


def test_value_read(read_mfile):
    assert isinstance(read_mfile.get_parameter_value("lsa"), int)


def test_write_json(read_mfile, temporary_dir):
    json_f = os.path.join(temporary_dir, "2017_baseline.json")
    read_mfile.write(json_f)
    assert os.path.exists(json_f)
    with open(json_f) as file:
        assert json.load(file)


def test_write_yaml(read_mfile, temporary_dir):
    pckl_f = os.path.join(temporary_dir, "2017_baseline.yml")
    read_mfile.write(pckl_f)
    assert os.path.exists(pckl_f)
    with open(pckl_f) as file:
        assert yaml.load(file, Loader=yaml.BaseLoader)


def test_write_pickle(read_mfile, temporary_dir):
    pckl_f = os.path.join(temporary_dir, "2017_baseline.pckl")
    read_mfile.write(pckl_f)
    assert os.path.exists(pckl_f)
    with open(pckl_f, "rb") as file:
        assert pickle.load(file)
