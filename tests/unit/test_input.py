import pytest
from process import fortran


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
    ],
)
def test_parse_real(epsvmc, expected, tmp_path):
    """Tests the parsing of real numbers into PROCESS.

    Program to get the expected value for 0.008 provided at https://github.com/ukaea/PROCESS/pull/3067
    """
    fortran.global_variables.fileprefix = _create_input_file(
        tmp_path, f"epsvmc = {epsvmc}"
    )
    fortran.init_module.init()

    assert fortran.numerics.epsvmc.item() == expected
