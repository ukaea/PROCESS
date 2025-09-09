import numpy as np
import pytest

import process.data_structure as data_structure
import process.init as init
import process.input as process_input
from process.exceptions import ProcessValidationError


def _create_input_file(directory, content: str):
    file_location = directory / "IN.DAT"
    with open(file_location, "w") as f:
        f.write(content)

    return str(file_location)


@pytest.mark.parametrize(
    ["epsvmc", "expected"],
    [
        (string, 1.0)
        for string in [
            "1.0",
            "1.0D0",
            "1.0d0",
            "1.0e0",
            "1.0E0",
            "1.0D+0",
            "1.0d+0",
            "1.0e+0",
            "1.0E+0",
            "1.0D-0",
            "1.0d-0",
            "1.0e-0",
            "1.0E-0",
            "0.10D1",
            "0.10d1",
            "0.10E1",
            "0.10e1",
            "0.10D+1",
            "0.10d+1",
            "0.10E+1",
            "0.10e+1",
            "10.0D-1",
            "10.0d-1",
            "10.0E-1",
            "10.0e-1",
        ]
    ]
    + [
        (string, 0.0080000000000000002)
        for string in ["0.008", "8.0E-3", "8.0D-3", "8.0d-3", "8.0e-3"]
    ]
    + [("0.546816593988753", 0.546816593988753)],
)
def test_parse_real(epsvmc, expected, tmp_path):
    """Tests the parsing of real numbers into PROCESS.

    Program to get the expected value for 0.008 provided at https://github.com/ukaea/PROCESS/pull/3067
    """
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, f"epsvmc = {epsvmc}"
    )
    init.init_process()

    assert data_structure.numerics.epsvmc == expected


@pytest.mark.parametrize(
    ["value"],
    [
        [0.546816593988753],
        [0.13134204235647895],
        [0.75],
        [0.7],
        [0.3],
        [0.1293140904093427],
    ],
)
def test_exact_parsing(value, tmp_path):
    """Tests the parsing of real numbers into PROCESS.

    These tests failed using the old input parser and serve to show that the Python parser generally
    produces more accurate floats and accumulates less error.
    """
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, f"epsvmc = {value}"
    )
    init.init_process()

    assert data_structure.numerics.epsvmc == value


def test_parse_input(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path,
        ("runtitle = my run title\nioptimz = -2\nepsvmc = 0.6\nboundl(1) = 0.5"),
    )
    init.init_process()

    assert data_structure.global_variables.runtitle == "my run title"
    assert data_structure.numerics.ioptimz == -2
    assert pytest.approx(data_structure.numerics.epsvmc) == 0.6
    assert pytest.approx(data_structure.numerics.boundl[0]) == 0.5


def test_input_choices(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, ("ioptimz = -1")
    )

    with pytest.raises(ProcessValidationError):
        init.init_process()


@pytest.mark.parametrize(
    ("input_file_value"), ((-0.01,), (1.1,)), ids=("violate lower", "violate upper")
)
def test_input_range(tmp_path, input_file_value):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, (f"epsvmc = {input_file_value}")
    )

    # check that the test data doesn't change
    assert process_input.INPUT_VARIABLES["epsvmc"].range == (0.0, 1.0)

    with pytest.raises(ProcessValidationError):
        init.init_process()


def test_input_array_when_not(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, ("epsvmc(1) = 0.5")
    )

    with pytest.raises(ProcessValidationError):
        init.init_process()


def test_input_not_array_when_is(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, ("boundl = 0.5")
    )

    with pytest.raises(ProcessValidationError):
        init.init_process()


def test_input_float_when_int(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, ("ioptimz = 0.5")
    )

    with pytest.raises(ProcessValidationError):
        init.init_process()


def test_input_array(tmp_path):
    data_structure.global_variables.fileprefix = _create_input_file(
        tmp_path, ("boundl = 0.1, 0.2, 1.0, 0.0, 1.0e2")
    )

    init.init_process()
    np.testing.assert_array_equal(
        data_structure.numerics.boundl[:6], [0.1, 0.2, 1.0, 0.0, 1.0e2, 0]
    )
