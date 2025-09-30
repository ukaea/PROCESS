"""Integration tests for the main.py module."""

import json
from shutil import copy

from process import main


def test_single_run(temp_data):
    """Test a SingleRun Process run with CLI args.

    This will just check that an exception isn't thrown.
    :param temp_data: temporary dir containing data files
    :type temp_data: Path
    """
    # Set input file path in temp_data dir
    input_path = temp_data / "large_tokamak_IN.DAT"
    input_file = str(input_path.resolve())

    # Run a SingleRun with an explicitly defined IN.DAT
    main.main(args=["-i", input_file])


def test_single_run_cwd(temp_data_cwd):
    """SingleRun without defining an input file.

    Try running without a defined input file (no args). This will look for
    an IN.DAT in the cwd.
    :param temp_data_cwd: temporary data dir, which is also the cwd
    :type temp_data_cwd: Path
    """
    # Copy input file to make a file named "IN.DAT"
    copy(temp_data_cwd / "large_tokamak_IN.DAT", temp_data_cwd / "IN.DAT")
    # Run: args must be emptylist; if None, argparse tries to use CLI args
    main.main(args=[])


def test_vary_run(temp_data):
    """Test a VaryRun with CLI args.

    :param temp_data: temporary dir containing data files
    :type temp_data: Path
    """
    # Set run_process.conf path in temp dir
    # Chosen because it's the only VaryRun in the test suite, and is fast
    conf_path = temp_data / "run_process.conf"
    conf_file = str(conf_path.resolve())

    # Run a VaryRun with an explicit conf file name
    main.main(args=["--varyiterparams", "--varyiterparamsconfig", conf_file])


def test_vary_run_cwd(temp_data_cwd):
    """Test VaryRun without explicitly defining the conf file name.

    This will look for a run_process.conf in the cwd.
    :param temp_data_cwd: temporary data dir, which is also the cwd
    :type temp_data_cwd: Path
    """
    main.main(args=["--varyiterparams"])


def test_plot_proc(temp_data, mfile_name):
    """Run plot proc via CLI.

    Currently, Process needs to run on an input file, then it can run plot_proc
    on a produced MFILE.DAT.
    :param temp_data: temporary dir containing data files
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type temp_data: str
    """
    # Specify input and mfiles
    input_file = temp_data / "large_tokamak_IN.DAT"
    input_file_str = str(input_file.resolve())

    # Run on input, then plot custom mfile name
    main.main(args=["-i", input_file_str, "--full-output"])

    # Assert a pdf has been created
    assert len(list(temp_data.glob("*.pdf")))


def test_single_run_with_mfilejson(temp_data):
    """Test a SingleRun Process run with CLI args including --mfilejson.

    This will check that the process runs without throwing an exception
    and a JSON output is produced in the working directory, then checks if
    the JSON is valid and contains expected keys.

    :param temp_data: temporary dir containing data files
    :type temp_data: Path
    """
    # Set input file path in temp_data dir.
    input_path = temp_data / "large_tokamak_eval.IN.DAT"
    mfile_path = temp_data / "large_tokamak_eval.MFILE.DAT"
    input_file = str(input_path.resolve())
    mfile = str(mfile_path.resolve())

    # Run a SingleRun with the --mfilejson flag.
    main.main(args=["-i", input_file, "--mfilejson", "-m", mfile])

    # Assert that 'large_tokamak_eval.MFILE.DAT.json' has been produced in the temp_data directory.
    expected_json = temp_data / "large_tokamak_eval.MFILE.DAT.json"
    assert expected_json.exists(), "large_tokamak_eval.MFILE.DAT.json was not found"

    # Check if the file contains valid JSON.
    try:
        with open(expected_json) as f:
            json_data = json.load(f)
    except json.JSONDecodeError as err:
        raise AssertionError("The JSON file is not valid JSON") from err

    # Check if the JSON contains expected outputs.
    expected_keys = ["rmajor", "b_plasma_toroidal_on_axis", "beta_total_vol_avg"]
    for key in expected_keys:
        assert key in json_data, f"Expected key '{key}' not found in the JSON file"
