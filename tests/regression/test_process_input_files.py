"""Regression tests for various PROCESS input files.

Run the input file and compare the resulting MFILE.DAT
to a reference MFile that is tracked on a remote repository.
This will indicate any differences in the MFile contents caused
by changes made off of main.
"""

import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest
from regression_test_assets import RegressionTestAssetCollector

from process.core.io.mfile import MFile
from process.main import main

logger = logging.getLogger(__name__)

INPUT_FILES_FOLDER = Path(__file__).resolve().parent / "input_files"
EXCLUSIONS = {
    "itvar",
    "xcm",
    "convergence_parameter",
    "sqsumsq",
    "nviter",
    "commsg",
    "procver",
    r"sig_tf_r_max\(1\)",  # weird value, flips between 0 and very low?
    r"normres[0-9]+",
    r"nitvar[0-9]+",
    "process_runtime",
}


@dataclass
class MFileVariableDifference:
    name: str
    ref: float
    new: float
    percentage_change: float


class RegressionTestScenario:
    def __init__(self, input_file: Path):
        """
        Represents an input scenario (input file) to PROCESS that is to be regression tested.

        :param input_file: absolute path of the input file (`<scenario_name>.IN.DAT`)
        :type input_file: Path
        """
        self.input_file = input_file
        self.scenario_name = input_file.name.replace(".IN.DAT", "")

    def run(self, solver: str):
        """Runs the scenario input file using PROCESS"""
        logger.info(
            f"Running regression test {self.scenario_name} using input file {self.input_file}"
        )
        try:
            main(["--input", str(self.input_file), "--solver", solver])
        except Exception as e:
            raise RuntimeError(
                f"\033[1;101m An error occured while running PROCESS: {e}\033[0m"
            ) from e

    def compare(
        self, reference_mfile_location: Path, tolerance: float, opt_params_only: bool
    ):
        """Runs assertions about the MFile with respect to a reference MFile

        :param reference_mfile_location: path to the downloaded reference MFile.
        :type reference_mfile_location: Path
        :param tolerance: percentage differences under this threshold will not cause
        the test to fail.
        :type tolerance: float
        :param opt_params_only: if `True`, only compare changes in optimisation
        parameters from the reference MFile.
        :type opt_params_only: bool
        """
        mfile_location = self.input_file.parent / f"{self.scenario_name}.MFILE.DAT"

        assert mfile_location.exists(), (
            f"PROCESS has not produced an MFile at the expected location {mfile_location}. "
            "Ensure the Scenario has been run!"
        )

        with open(mfile_location) as f:
            assert len(f.readlines()) > 0, (
                "An MFile has been created, but it is empty, "
                "indicating PROCESS did not run the input file successfully!"
            )

        mfile = MFile(str(mfile_location))
        reference_mfile = MFile(str(reference_mfile_location))

        ifail = mfile.data["ifail"].get_scan(-1)

        assert ifail == 1 or mfile.data["ioptimz"].get_scan(-1) == -2, (
            f"\033[0;36m ifail of {ifail} indicates PROCESS did not solve successfully\033[0m"
        )

        mfile_keys = set(mfile.data.keys())
        reference_mfile_keys = set(reference_mfile.data.keys())
        key_mfile_not_ref = mfile_keys - reference_mfile_keys
        key_ref_not_mfile = reference_mfile_keys - mfile_keys

        key_ref_not_mfile_msg = (
            "\033[0;35m Reference MFile contains variables that are not present in "
            f"the MFILE: {key_ref_not_mfile} \033[0m"
        )
        if key_ref_not_mfile:
            logger.warning(key_ref_not_mfile_msg)

        key_mfile_not_ref_msg = (
            "\033[0;35m MFile contains variables that are not present in "
            f"the reference MFILE: {key_mfile_not_ref} \033[0m"
        )
        if key_mfile_not_ref:
            logger.warning(key_mfile_not_ref_msg)

        differences = self.mfile_value_changes(
            reference_mfile, mfile, tolerance, opt_params_only
        )
        if differences:
            differences = sorted(
                differences, key=lambda i: abs(i.percentage_change), reverse=True
            )

            logger.warning(
                f"{'Variable':20}\t{'Ref':>10}\t{'New':>10}\t{'% Change':>10}"
            )
            logger.warning("-" * 60)
            for diff in differences:
                logger.warning(
                    f"{diff.name:20}\t{diff.ref:10.3g}\t{diff.new:10.3g}\t{diff.percentage_change:10.2f}"
                )

            assert len(differences) == 0, (
                f"\033[0;32m {len(differences)} differences: the reference MFile contains different values "
                "for some of the variables. See the warnings for a breakdown of the differences.\033[0m"
            )

        assert not key_ref_not_mfile, key_ref_not_mfile_msg
        assert not key_mfile_not_ref, key_mfile_not_ref_msg

    @staticmethod
    def mfile_value_changes(
        ref: MFile, new: MFile, tolerance: float, opt_params_only: bool
    ) -> list[MFileVariableDifference]:
        """Calculates the differences between two MFiles.

        :param ref: the reference MFile
        :type ref: MFile
        :param new: the MFile generated by running the input file with the
        latest changes, to be compared with the reference MFile.
        :type new: MFile
        :param tolerance: percentage differences under this threshold will not cause
        the test to fail.
        :type tolerance: float
        :param opt_params_only: if `True`, only compare changes in optimisation
        parameters from the reference MFile.
        :type opt_params_only: bool
        """
        common_keys = set(ref.data.keys()).intersection(set(new.data.keys()))

        # exclude any variable that matches an element of EXCLUSIONS
        exclusions = r"^(?:" + r"|".join(EXCLUSIONS) + r")$"
        if opt_params_only:
            # exclude any variable thats not the optimisation parameters
            # NOTE: this excludes anything in EXCLUSIONS
            exclusions += r"|(?!norm_objf|itvar[0-9]+)"

        diffs = []
        for key in common_keys:
            # the variable name matches an exclusion? Ignore it
            if re.match(exclusions, key) is not None:
                continue

            ref_value = ref.data[key].get_scan(-1)
            new_value = new.data[key].get_scan(-1)

            try:
                ref_value = float(ref_value)
                new_value = float(new_value)
            except ValueError:
                # only compare float-able values
                continue

            # Define relative tolerance
            # Use pytest's default relative tolerance (1e-6)
            # 0 tolerance causes floating-point discrepancies
            # between local and CI runs or tolerance
            # is a percentage, rel arg takes a fraction
            rel_tolerance = None if (tolerance == 0) else (tolerance / 100)

            try:
                # Use pytest.approx for relative and absolute comparisons:
                # handles values close to 0
                assert new_value == pytest.approx(ref_value, rel=rel_tolerance)
            except AssertionError:
                # NOTE: the percentage change for a value that was originally 0
                # is 100% NOT 0% because this indicates the change better
                percentage_change = (
                    100.0 * (new_value - ref_value) / abs(ref_value)
                    if ref_value != 0
                    else 100.0
                )

                diffs.append(
                    MFileVariableDifference(
                        key,
                        ref_value,
                        new_value,
                        percentage_change,
                    )
                )

        return diffs


@pytest.fixture(scope="session")
def tracked_regression_test_assets():
    """Session fixture providing a RegressionTestAssetCollector
    for finding remote tracked MFiles.

    This fixture creates one asset collector that is shared
    between all regression tests and reduces the number of
    API calls made to the remote repository."""
    return RegressionTestAssetCollector()


@pytest.mark.parametrize(
    ["input_file"],
    [[f] for f in INPUT_FILES_FOLDER.glob("*.IN.DAT")],
    ids=lambda v: v.stem.replace(".IN", ""),
)
def test_input_file(
    input_file: Path,
    tmp_path,
    solver_name: str,
    tracked_regression_test_assets,
    reg_tolerance: float,
    opt_params_only: bool,
):
    """Tests each input file in the 'input_files' directory.

    The test will locate and download a remote reference MFile that was
    generated by running the input file on the 'main' branch.

    The input file will then be run locally and compared to the reference file.
    The test will fail if:
    * The input file fails to run successfully locally
    * Any variable has a different value between the two files
    * The reference MFile contains variables that the new file does not
    * The new MFile contains variables that the reference file does not

    :param input_file: a file in 'input_files' matching the '*.IN.DAT' pattern
    :type input_file: Path
    :param tmp_path: temporary directory for test assets to
    be created and downloaded into.
    :type tmp_path: Path
    :param solver_name: specifies the PROCESS solver to use, 'vmcon' unless
    specified on the command line.
    :type solver_name: str
    :param tracked_regression_test_assets: object providing access to tracked MFiles
    :type tracked_regression_test_assets: RegressionTestAssetCollector
    :param reg_tolerance: user specified tolerance, percentage differences below which
    are ignored.
    :type reg_tolerance: float
    :param opt_params_only: if True, user specificied that only optimisation parameters
    should be compared in the test.
    :type opt_params_only: bool
    """
    new_input_file = tmp_path / input_file.name
    shutil.copy(input_file, new_input_file)

    stellarator_config = (
        input_file / f"../{input_file.name.replace('IN.DAT', 'stella_conf.json')}"
    ).resolve()
    if stellarator_config.exists():
        shutil.copy(stellarator_config, tmp_path / stellarator_config.name)

    scenario = RegressionTestScenario(new_input_file)

    reference_mfile = tracked_regression_test_assets.get_reference_mfile(
        scenario.scenario_name
    )

    scenario.run(solver_name)

    # reference MFile cannot be found?
    # raise an error after the file is run so that any errors while running the input file
    # are raised first.
    if reference_mfile is None:
        raise RuntimeError(
            "\033[0;36m No reference input file exists (so cannot compare results). The input file ran without any exceptions.\033[0m"
        )

    scenario.compare(reference_mfile, reg_tolerance, opt_params_only)
