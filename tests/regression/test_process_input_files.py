from pathlib import Path
from dataclasses import dataclass
from typing import List
import shutil
import logging
import re

import pytest
from process.main import main
from process.io.mfile import MFile

from regression_test_assets import RegressionTestAssetCollector


logger = logging.getLogger(__name__)

INPUT_FILES_FOLDER = Path(__file__).resolve().parent / "input_files"
EXCLUSIONS = {
    "itvar",
    "xcm",
    "convergence_parameter",
    "sqsumsq",
    "nviter",
    r"sig_tf_r_max\(1\)",  # weird value, flips between 0 and very low?
    r"normres[0-9]+",
    r"nitvar[0-9]+",
}


@dataclass
class MFileVariableDifference:
    name: str
    ref: float
    new: float
    percentage_change: float


class RegressionTestScenario:
    def __init__(self, input_file: Path) -> None:
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
        main(["--input", str(self.input_file), "--solver", solver])

    def compare(
        self, reference_mfile_location: Path, tolerance: float, opt_params_only: bool
    ):
        """Runs assertions about the MFile with respect to a reference MFile"""
        mfile_location = self.input_file.parent / f"{self.scenario_name}.MFILE.DAT"

        assert mfile_location.exists(), (
            f"PROCESS has not produced an MFile at the expected location {mfile_location}. "
            "Ensure the Scenario has been run!"
        )

        with open(mfile_location, "r") as f:
            assert len(f.readlines()) > 0, (
                "An MFile has been created, but it is empty, "
                "indicating PROCESS did not run the input file successfully!"
            )

        mfile = MFile(str(mfile_location))
        reference_mfile = MFile(str(reference_mfile_location))

        assert (
            ifail := mfile.data["ifail"].get_scan(-1)
        ) == 1, f"ifail of {ifail} indicates PROCESS did not solve successfully"

        differences = self.mfile_value_changes(
            reference_mfile, mfile, tolerance, opt_params_only
        )
        if differences:
            logger.warning(
                f"{'Variable':20}\t{'Ref':>10}\t{'New':>10}\t{'% Change':>10}"
            )
            logger.warning("-" * 60)
            for diff in differences:
                logger.warning(
                    f"{diff.name:20}\t{diff.ref:10.3g}\t{diff.new:10.3g}\t{diff.percentage_change:10.2f}"
                )

            assert False, (
                "The reference MFile contains different values for some of "
                "the variables. See the warnings for a breakdown of the differences."
            )

        mfile_keys = set(mfile.data.keys())
        reference_mfile_keys = set(reference_mfile.data.keys())
        key_mfile_not_ref = mfile_keys - reference_mfile_keys
        key_ref_not_mfile = reference_mfile_keys - mfile_keys

        assert not key_ref_not_mfile, (
            "Reference MFile contains variables that are not present in "
            f"the MFILE: {key_ref_not_mfile}"
        )

        assert not key_mfile_not_ref, (
            "MFile contains variables that are not present in "
            f"the reference MFILE: {key_mfile_not_ref}"
        )

    @staticmethod
    def mfile_value_changes(
        ref: MFile, new: MFile, tolerance: float, opt_params_only: bool
    ) -> List[MFileVariableDifference]:
        common_keys = set(ref.data.keys()).intersection(set(new.data.keys()))

        if opt_params_only:
            exclusions = r"^(?:norm_objf|itvar[0-9]+)$"
        else:
            exclusions = r"^(?:" + r"|".join(EXCLUSIONS) + r")$"

        diffs = []
        for key in common_keys:
            if (re.search(exclusions, key) is not None and not opt_params_only) or (
                re.search(exclusions, key) is None and opt_params_only
            ):
                continue

            if (ref_value := ref.data[key].get_scan(-1)) == (
                new_value := new.data[key].get_scan(-1)
            ):
                continue

            try:
                ref_value = float(ref_value)
                new_value = float(new_value)
            except ValueError:
                # only compare float-able values
                continue

            # NOTE: the percentage change for a value that was
            # originally 0 is 100% NOT 0% because this allows it
            # to get picked up by various checks below
            percentage_change = (
                100 * (new_value - ref_value) / abs(ref_value)
                if ref_value != 0
                else 100
            )

            if percentage_change <= tolerance:
                continue

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
    return RegressionTestAssetCollector()


def id_input_file(val: Path):
    return val.stem.replace(".IN", "")


@pytest.mark.parametrize(
    ["input_file"],
    [[f] for f in INPUT_FILES_FOLDER.glob("*.IN.DAT")],
    ids=id_input_file,
)
def test_input_file(
    input_file,
    tmp_path,
    solver_name: str,
    tracked_regression_test_assets,
    reg_tolerance: float,
    opt_params_only: bool,
):
    new_input_file = tmp_path / input_file.name
    shutil.copy(input_file, new_input_file)

    scenario = RegressionTestScenario(new_input_file)

    reference_mfile = tracked_regression_test_assets.get_reference_mfile(
        scenario.scenario_name, tmp_path
    )

    # reference MFile cannot be found?
    # should the file be allowed to run just to test it converges (with a warning about no comparison)?
    if reference_mfile is None:
        pytest.skip(
            reason=f"A reference input file cannot be found for {scenario.scenario_name}"
        )

    scenario.run(solver_name)
    scenario.compare(reference_mfile, reg_tolerance, opt_params_only)
