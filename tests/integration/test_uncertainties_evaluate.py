"""Integration tests for uncertainties/evaluate.py."""

import json
from shutil import move

import pytest

from process.uncertainties import (
    evaluate_uncertainties,
    hdf_to_scatter_plot,
    morris_plotting,
    sobol_plotting,
)

# Uncertainties tests currently take too long for a 1h CI job
pytest.skip(
    "Uncertainties tests currently time out integration test jobs",
    allow_module_level=True,
)


def test_evaluate_uncertainties_monte_carlo(temp_data):
    """Run evaluate uncertainties with a config file argument
       and with the monte carlo technique flag and check for output .h5 file.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    # Rename the configfile to something custom
    old_config_path = temp_data / "config_evaluate_uncertainties.json"
    new_config_path = temp_data / "customconfig.json"
    move(old_config_path, new_config_path)

    # Set temp_data dir as working dir, and use input file in it
    # Change input file and working dir paths in config file to the temp_data
    # dir path, only known at runtime
    with open(new_config_path) as config_file:
        config_obj = json.load(config_file)
        config_obj["config"]["IN.DAT_path"] = str(
            temp_data / config_obj["config"]["IN.DAT_path"]
        )
        config_obj["config"]["working_directory"] = str(
            temp_data / config_obj["config"]["working_directory"]
        )

    with open(new_config_path, "w") as config_file:
        json.dump(config_obj, config_file)

    # Run with custom config file path (custom input file path in config file)
    config_path_str = str(new_config_path)
    evaluate_uncertainties.main(
        args=["--configfile", config_path_str, "--method", "monte_carlo"]
    )

    # assert uncertainties_data.h5 file has been created
    assert len(list(temp_data.glob("*.h5")))  # can assert the created text file?


def test_evaluate_uncertainties_monte_carlo_non_optimised(temp_data):
    """Run evaluate uncertainties with a config file argument
       and with the monte carlo technique flag and check for output .h5 file
       where input file is for a non optimised run

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    # Rename the configfile to something custom
    old_config_path = temp_data / "config_evaluate_uncertainties_nonopt.json"
    new_config_path = temp_data / "customconfig.json"
    move(old_config_path, new_config_path)

    # Set temp_data dir as working dir, and use input file in it
    # Change input file and working dir paths in config file to the temp_data
    # dir path, only known at runtime
    with open(new_config_path) as config_file:
        config_obj = json.load(config_file)
        config_obj["config"]["IN.DAT_path"] = str(
            temp_data / config_obj["config"]["IN.DAT_path"]
        )
        config_obj["config"]["working_directory"] = str(
            temp_data / config_obj["config"]["working_directory"]
        )

    with open(new_config_path, "w") as config_file:
        json.dump(config_obj, config_file)

    # Run with custom config file path (custom input file path in config file)
    config_path_str = str(new_config_path)
    evaluate_uncertainties.main(
        args=["--configfile", config_path_str, "--method", "monte_carlo"]
    )

    # assert uncertainties_data.h5 file has been created
    assert len(list(temp_data.glob("*.h5")))  # can assert the created text file?


def test_evaluate_uncertainties_morris_method(temp_data):
    """Run evaluate uncerainties with a configfile argument
       and with the morris method technique flag and check for output .txt file.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    # Rename the configfile to something custom
    old_config_path = temp_data / "config_evaluate_uncertainties.json"
    new_config_path = temp_data / "customconfig.json"
    move(old_config_path, new_config_path)

    # Set temp_data dir as working dir, and use input file in it
    # Change input file and working dir paths in config file to the temp_data
    # dir path, only known at runtime
    with open(new_config_path) as config_file:
        config_obj = json.load(config_file)
        config_obj["config"]["IN.DAT_path"] = str(
            temp_data / config_obj["config"]["IN.DAT_path"]
        )
        config_obj["config"]["working_directory"] = str(
            temp_data / config_obj["config"]["working_directory"]
        )

    with open(new_config_path, "w") as config_file:
        json.dump(config_obj, config_file)

    # Run with custom config file path (custom input file path in config file)
    config_path_str = str(new_config_path)
    evaluate_uncertainties.main(
        args=["--configfile", config_path_str, "--method", "morris_method"]
    )

    # assert morris_method.txt file has been created
    assert len(list(temp_data.glob("*.txt")))


def test_evaluate_uncertainties_sobol_method(temp_data):
    """Run evaluate uncerainties with a configfile argument
       and with the sobol method technique flag and check for output .txt file.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    # Rename the configfile to something custom
    old_config_path = temp_data / "config_evaluate_uncertainties.json"
    new_config_path = temp_data / "customconfig.json"
    move(old_config_path, new_config_path)

    # Set temp_data dir as working dir, and use input file in it
    # Change input file and working dir paths in config file to the temp_data
    # dir path, only known at runtime
    with open(new_config_path) as config_file:
        config_obj = json.load(config_file)
        config_obj["config"]["IN.DAT_path"] = str(
            temp_data / config_obj["config"]["IN.DAT_path"]
        )
        config_obj["config"]["working_directory"] = str(
            temp_data / config_obj["config"]["working_directory"]
        )

    with open(new_config_path, "w") as config_file:
        json.dump(config_obj, config_file)

    # Run with custom config file path (custom input file path in config file)
    config_path_str = str(new_config_path)
    evaluate_uncertainties.main(
        args=["--configfile", config_path_str, "--method", "sobol_method"]
    )

    # assert sobol.txt file has been created
    assert len(list(temp_data.glob("*.txt")))


def test_morris_plotting(temp_data):
    """Run morris_plotting on an output morris_method_output.txt
       and check for an output pdf.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    input_raw = temp_data / "morris_method_output.txt"
    input_str = str(input_raw)

    output = temp_data / "morris_output.pdf"
    output_str = str(output)

    morris_plotting.main(args=["-f", input_str, "-o", output_str])

    # Assert the default output pdf has been created
    assert len(list(temp_data.glob("morris_output.pdf")))


def test_sobol_plotting(temp_data):
    """Run plot_proc on an output sobol.txt file
       and check for an output pdf.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    input_raw = temp_data / "sobol.txt"
    input_str = str(input_raw)

    output = temp_data / "sobol_output.pdf"
    output_str = str(output)

    sobol_plotting.main(args=["-f", input_str, "-o", output_str])

    # Assert the default output pdf has been created
    assert len(list(temp_data.glob("sobol_output.pdf")))


def test_hdf_to_scatter_plot(temp_data):
    """Run read_hdf on an output uncertainties_data.h5 file
       and check for an output pdf.

    :param temp_data: temporary data dir
    :type temp_data: Path
    """
    output = temp_data / "uncertainties_data.h5"
    output_str = str(output)
    hdf_to_scatter_plot.main(args=["-i", output_str, "-v", "rmajor"])
