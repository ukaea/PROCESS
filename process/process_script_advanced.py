#!/usr/bin/env python
"""Script to build and run PROCESS simply.

This compiles, runs and displays output from PROCESS. The aim is to provide a
common entry point to the code in Python, and to automate a standard workflow.
"""
# TODO This requires converting to run the Python-wrapped version of Process
import subprocess
import argparse
from shutil import copy
from pathlib import Path
from difflib import unified_diff
from utilities.process_io_lib import input_validator

# Set required paths outside the Python package directory
ROOT_DIR = Path(__file__).parent.parent
# The root project directory (for running cmake)
PROCESS_EXE_PATH = ROOT_DIR.joinpath("bin/process.exe")
# The Process executable (for running Process)


class Process(object):
    """A Process workflow based on command line arguments."""

    def __init__(self):
        """Parse command line arguments and run accordingly."""
        # File paths
        self.run_dir = None
        self.input = None
        self.ref_input = None

        # Run actions
        self.parse_args()
        if self.args.build:
            self.build_process()
            self.create_dicts()

        # Find the input file and look for changes to it
        self.set_input()
        self.set_ref_input()
        self.input_file_diff()

        self.validate_input()
        self.run_process()
        if self.args.util:
            # Only create dicts if necessary for utilities
            self.create_dicts()
            self.run_utils()

    def parse_args(self):
        """Parse the command line arguments."""
        parser = argparse.ArgumentParser(
            description="Script to automate " "running PROCESS"
        )

        # Optional arguments
        parser.add_argument(
            "--input",
            "-i",
            metavar="input_file_path",
            help="The path to the input file that Process runs on",
        )
        parser.add_argument(
            "--ref_input",
            "-r",
            metavar="reference_input",
            help="Reference input file to record changes against",
        )
        parser.add_argument(
            "--build", "-b", action="store_true", help="Rebuilds Process if present"
        )
        parser.add_argument(
            "--util",
            "-u",
            metavar="util_name",
            help="Utility to run after Process runs",
        )

        self.args = parser.parse_args()
        # Store namespace object of the args

    def build_process(self):
        """Build Process"""
        # cmake must be run from the project root dir
        subprocess.run(["cmake", "-H.", "-Bbuild"], cwd=ROOT_DIR)
        subprocess.run(["cmake", "--build", "build"], cwd=ROOT_DIR)

    def set_input(self):
        """Try to find or validate the input file path, then store it.

        Also set the run directory according to the input file path.
        """
        if self.args.input:
            # Input file path specified as CLI argument
            input_path = Path(self.args.input)
            if input_path.exists():
                self.input = input_path
                # Set input as Path object
            else:
                raise FileNotFoundError("Input file not found on this path.")
        else:
            # No input file specified; try to find one in the current dir
            input_files = [
                path
                for path in Path.cwd().iterdir()
                if "IN.DAT" in path.name and "REF_IN.DAT" not in path.name
            ]

            if len(input_files) == 0:
                # No input file found
                raise FileNotFoundError("Input file not found in this dir.")
            elif len(input_files) == 1:
                # Only one input file in this dir; use it
                self.input = input_files[0]
            else:
                # More than one input file; ambiguous which one to use
                raise Exception(
                    "More than one IN.DAT found in this dir. "
                    'Specifiy which one to use with the "-i" option, or '
                    "remove the other ones."
                )

        # self.input is now defined
        self.run_dir = self.input.parent
        # Set run directory as dir that contains the input file

    def set_ref_input(self):
        """Find an input file to use as a reference for changes.

        Find or create a reference input file to compare with the current
        input file, to allow the changes to be tracked.
        """
        # self.args.ref_input: reference input file path; can already be set as
        # command line arg
        if self.args.ref_input:
            # Ref input specified as command line arg. Check it exists
            path = Path(self.args.ref_input)
            if path.exists():
                # Store path object
                self.ref_input = path
            else:
                raise FileNotFoundError(
                    "Reference input file not found; " "check the path."
                )

            # Input file exists
            if "REF_IN.DAT" in self.ref_input.name:
                # Ref file specified; got a ref file
                pass
            elif "IN.DAT" in self.ref_input.name:
                # Regular input file specified. Make a reference copy as
                # REF_IN.DAT and update self.ref_input
                self.create_ref_input()
            else:
                # File isn't an IN.DAT or a REF_IN.DAT
                raise ValueError(
                    "Reference input file should end in IN.DAT or " "REF_IN.DAT"
                )
        else:
            # Ref input not specified: try to find one in the input file's dir
            ref_inputs = list(self.run_dir.glob("*REF_IN.DAT"))
            if len(ref_inputs) == 0:
                # No ref input files: create one from input file
                self.create_ref_input()
            elif len(ref_inputs) == 1:
                # One pre-existing ref input file: set as ref input
                self.ref_input = ref_inputs[0]
            else:
                # More than one ref file: error
                raise Exception(
                    "More than one REF_IN.DAT found in this dir. "
                    'Specifiy which one to use with the "-r" option, or '
                    "remove the other ones."
                )

        # Now either have raised an exception or have a ref input file set,
        # with suffix "REF_IN.DAT"

    def create_ref_input(self):
        """Create a reference input file from a pre-existing input file."""
        # First, ensure that self.ref_input is set as an IN.DAT to copy
        if not self.ref_input:
            # No IN.DAT specified to use as ref: use input file instead
            self.ref_input = self.input

        # Copy the ref input and save it with the REF_IN.DAT extension
        # Then set this as the new ref input
        new_ref_input_name = self.ref_input.name.replace("IN.DAT", "REF_IN.DAT")
        new_ref_input = Path(self.run_dir).joinpath(new_ref_input_name)
        copy(self.ref_input, new_ref_input)
        self.ref_input = new_ref_input

    def input_file_diff(self):
        """Perform a diff between the reference input and input files.

        This compares the REF_IN.DAT and IN.DAT files to highlight the
        modifications.
        """
        # Run a diff on the files
        before_file = open(self.ref_input)
        after_file = open(self.input)
        before = before_file.readlines()
        after = after_file.readlines()
        before_file.close()
        after_file.close()
        diff = unified_diff(
            before, after, fromfile=self.ref_input.name, tofile=self.input.name
        )

        # Write diff to output file
        diff_path = self.run_dir.joinpath("input.diff")
        with open(diff_path, "w") as diff_file:
            diff_file.writelines(diff)

    def create_dicts(self):
        """Create Python dictionaries"""
        # cmake must be run from the project root directory
        subprocess.run(["cmake", "--build", "build", "--target", "dicts"], cwd=ROOT_DIR)

    def validate_input(self):
        """Validate the input file using the input_validator module."""
        input_validator.validate(self.input)

    def run_process(self):
        """Run Process using the input file path."""
        subprocess.run([PROCESS_EXE_PATH, self.input])

    def run_utils(self):
        """Run a utility if specified on the command line."""
        # TODO allow multiple utils to run
        # TODO allow options to be passed to utils
        subprocess.run(["python", "./utilities/" + self.args.util + ".py"])


def main():
    """Run Process."""
    print("Running process.py")
    Process()


if __name__ == "__main__":
    main()
