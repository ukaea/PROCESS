from typing import NamedTuple
import pytest

from process.fortran import read_radiation, impurity_radiation_module


@pytest.fixture(scope="function")
def init_read_radiation():
    read_radiation.init_read_radiation()
    impurity_radiation_module.init_impurity_radiation_module()


class ReadLzParam(NamedTuple):
    element: str
    te: float
    netau: float
    expected: float
    mean_z: bool = False
    mean_qz: bool = False


@pytest.mark.parametrize(
    "read_lz_param",
    (
        ReadLzParam("He", 5, 0.5, 5.958312853655463e-35),
        ReadLzParam("He", 124.45, 0.5, 1.9852736835073994, mean_z=True),
        ReadLzParam("He", 124.45, 0.5, 3.95858840328217, mean_qz=True),
    ),
)
def test_read_lz(read_lz_param: ReadLzParam):
    ret = read_radiation.read_lz(
        read_lz_param.element,
        read_lz_param.te,
        read_lz_param.netau,
        read_lz_param.mean_z,
        read_lz_param.mean_qz,
        False,
    )

    assert ret == pytest.approx(read_lz_param.expected)
