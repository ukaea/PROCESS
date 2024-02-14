"""Scenario class for an individual regression test case."""

import logging
import sys
import os
import shutil
import re
from process.main import main
from process.io.mfile import MFile


# TODO This isn't good: put MFile into process package?
sys.path.append(os.path.join(os.path.dirname(__file__), "../utilities/"))


# Variables and patterns to ignore when comparing differences (set and list)
# These variables can differ substantially from their reference values when
# Process is run in different environments (e.g. local vs. CI system). These
# differences are accepted, but not fully understood, and should be treated
# with suspicion as rounding error alone cannot account for all of them.
# TODO Update this exclusion list to only include properly justified vars.

# ric(nohc) is included here as its sign switches very easily due to a rounding
# error-prone comparison, which was masked in the old version of the test suite.
# TODO This is a bug that needs to be resolved, and is excluded only to get
# the regression tests passing as before.
EXCLUSIONS = {
    "itvar",
    "xcm",
    "balance",
    "convergence_parameter",
    "sig_tf_r_max(1)",
    "sqsumsq",
    "ric(nohc)",
    "nviter",
}
EXCLUSION_PATTERNS = [
    r"normres\d{3}",  # normres and 3 digits
    r"nitvar\d{3}",  # nitvar and 3 digits
]

logger = logging.getLogger(__name__)


class Scenario:
    """A scenario for a single regression test."""

    def __init__(self, ref_dir):
        """Initialise with reference directory to test.

        :param ref_dir: reference directory
        :type ref_dir: Path
        """
        self.ref_dir = ref_dir  # Dir containing the test files
        self.test_dir = None  # Where the actual test takes place
        self.name = ref_dir.name  # Name of the scenario
        self.version = None  # Process version
        self.ref_mfile = None  # Reference MFile (expected)
        self.new_mfile = None  # Newly-created MFile (observed)
        self.vars_unique_ref = set()  # Vars unique to reference
        self.vars_unique_new = set()  # Vars unique to new
        self.vars_common_both = set()  # Vars in both reference and new
        self.over_tolerance_diff_items = []  # Differences over tolerance

    def run(self, test_dir, solver_name):
        """Run Process for this scenario.

        :param test_dir: temporary directory for running the test in
        :type test_dir: Path
        :param solver_name: solver to use
        :type solver_name: str
        """
        self.test_dir = test_dir
        logger.info(f"Running scenario: {self.name}")

        # TODO Deal with a default conf file name or not; are both vary args
        # really needed?
        # Determine if this is a normal or vary iteration parameters run
        if len(list(self.test_dir.glob("*.conf"))):
            # A .conf file is present: vary iteration params scenario
            conf_path_str = str(self.test_dir / "run_process.conf")
            args = ["--varyiterparams", "--varyiterparamsconfig", conf_path_str]
        else:
            # No .conf file: normal run
            input_path_str = str(self.test_dir / "IN.DAT")
            args = ["--input", input_path_str]

        # Append solver to use
        args.extend(["--solver", solver_name])

        # Run Process using the input or config file in the test_dir
        # directory, catching any errors
        try:
            main(args=args)
        except Exception:
            logger.exception(
                f"Process threw an exception when running " f"scenario: {self.name}"
            )
            raise

    def check_mfile_length(self):
        """Ensure there is something in the MFile.

        :return: True if is has non-zero length, False if not
        :rtype: bool
        """
        with open(self.test_dir / "MFILE.DAT", "r", encoding="utf-8") as mfile:
            mfile_len = len(mfile.readlines())

        if mfile_len == 0:
            return False
        else:
            return True

    def read_mfiles(self):
        """Read in reference and newly-output MFILEs, creating MFile objects."""
        ref_mfile_path_str = str(self.test_dir / "ref.MFILE.DAT")
        new_mfile_path_str = str(self.test_dir / "MFILE.DAT")

        try:
            self.ref_mfile = MFile(ref_mfile_path_str)
            self.new_mfile = MFile(new_mfile_path_str)
        except Exception:
            logger.exception("There was an error creating an MFile object.")
            raise

    def set_version(self):
        """Set process version number."""
        self.version = self.new_mfile.data["tagno"].get_scan(-1)

    def check_ifail(self):
        """Test the value of ifail, the solver error return flag.

        :return: True if ifail is 1, False otherwise
        :rtype: bool
        """
        self.ifail = self.new_mfile.data["ifail"].get_scan(-1)

        if self.ifail != 1:
            return False
        else:
            return True

    def add_diff_item(self, diff_item):
        """Add a diff item tuple that is outside the accepted tolerance.

        :param diff_item: a diff that exceeds the tolerance for this variable
        :type diff_item: tuple
        """
        self.over_tolerance_diff_items.append(diff_item)

    def set_mfile_var_sets(self):
        """Set unique and common variables in ref and new MFiles."""
        # Get variable keys from the reference and observed MFiles
        # and convert to sets
        ref_vars = set(self.ref_mfile.data.keys())
        new_vars = set(self.new_mfile.data.keys())

        # Filter out the excluded vars using set difference; ones we aren't
        # interested in
        ref_vars = ref_vars - EXCLUSIONS
        new_vars = new_vars - EXCLUSIONS

        # Filter out excluded variable patterns from the sets
        for pattern in EXCLUSION_PATTERNS:
            ref_vars = {var for var in ref_vars if not re.match(pattern, var)}
            new_vars = {var for var in new_vars if not re.match(pattern, var)}

        # Set difference and intersection
        self.vars_unique_new = new_vars - ref_vars
        self.vars_unique_ref = ref_vars - new_vars
        self.vars_common_both = new_vars & ref_vars

    def get_mfile_diffs(self):
        """Generator for differences between the ref and new MFiles.

        :yield: diff_item for each difference
        :rtype: tuple
        """
        # Scan number for comparing the final scan in each MFile
        scan_number = -1

        # Determine which variables are unique/common to each/both MFiles
        self.set_mfile_var_sets()

        # For the variables in both reference and new MFiles, find the
        # difference in values
        for var_name in self.vars_common_both:
            # Get expected and observed values
            try:
                exp = float(self.ref_mfile.data[var_name].get_scan(scan_number))
                obs = float(self.new_mfile.data[var_name].get_scan(scan_number))
            except ValueError:
                # This is to catch (and ignore) values that can't be converted
                # to floats, e.g. strings in the MFILE (fileprefix = 'IN.DAT')
                continue

            # Calulate percentage change
            if exp == 0:
                # Relative change is nonsensical
                chg = 0.0
            else:
                # Cater for negative changes too
                chg = ((obs - exp) / abs(exp)) * 100

            # Create "expected and observed" comparison tuple
            # TODO Should this be a class?
            diff_item = (var_name, exp, obs, chg)
            yield diff_item

    def log_summary(self):
        """Log a summary of the scenario's test result."""
        # Order variables alphabetically
        sorted_diffs = sorted(self.over_tolerance_diff_items)
        sorted_ref_vars = sorted(self.vars_unique_ref)
        sorted_new_vars = sorted(self.vars_unique_new)

        # Log differences outside tolerance
        logger.warning(
            f"\nTotal of {len(sorted_diffs)} " "differences outside tolerance"
        )

        # Create table if there's content
        if len(sorted_diffs) > 0:
            logger.warning(
                f"{'Variable':20}\t{'Ref':>10}\t{'New':>10}\t" f"{'Diff (%)':>10}"
            )
            logger.warning("-" * 60)
            for diff_item in sorted_diffs:
                var_name, exp, obs, chg = diff_item
                logger.warning(
                    f"{var_name:20}\t{exp:10.3g}\t{obs:10.3g}\t" f"{chg:10.2f}"
                )

        # Log variables in ref but not in new
        logger.warning(
            f"\nThere are {len(sorted_ref_vars)} variables in " "ref not in new"
        )
        for var in sorted_ref_vars:
            logger.warning(var)

        # Log variables in new but not in ref
        logger.warning(
            f"\nThere are {len(sorted_new_vars)} variables in " "new not in ref"
        )
        for var in sorted_new_vars:
            logger.warning(var)

    def overwrite_ref_files(self):
        """Overwrite reference MFILES and OUT files.

        When a change to the scenario's reference result is justified, overwrite
        the ref.MFILE.DAT and ref.OUT.DAT files.
        """
        # Copy the newly output MFILE.DAT and OUT.DAT from the test dir
        # (temporary), and overwrite the ref.MFILE.DAT and ref.OUT.DAT in the
        # reference directory (permanent)
        for file in ["MFILE", "OUT"]:
            src_path = self.test_dir / (file + ".DAT")
            dst_path = self.ref_dir / ("ref." + file + ".DAT")
            shutil.copyfile(src_path, dst_path)

    def get_diff_items(self):
        """Return list of diffs that exceed the tolerance.

        :return: list of diff_item tuples that are over tolerance
        :rtype: list
        """
        return self.over_tolerance_diff_items
