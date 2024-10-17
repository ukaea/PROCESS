"""Run the tracked files and move into tracking/ directory."""

import argparse
from pathlib import Path
import subprocess
import shutil


def run_and_move_tracked_files():
    # Stellarator doesn't converge so don't want to track it
    exclusions = ["stellarator"]

    path_to_tracking_files = Path("tests/regression/input_files")
    all_files = path_to_tracking_files.iterdir()

    all_files = []
    for file in list(path_to_tracking_files.iterdir()):
        if file.is_file():
            all_files.append(file.name)

    input_files = []
    for file in all_files:
        for excl in exclusions:
            if excl not in file and file.endswith("IN.DAT"):
                input_files.append(file)

    mfile_names = []
    for in_dat in input_files:
        file_prefix = in_dat.removesuffix("IN.DAT")
        mfile_name = file_prefix + "MFILE.DAT"
        mfile_names.append(mfile_name)

    new_location_mfile_names = []
    for current_name in mfile_names:
        prefix = current_name.removesuffix(".MFILE.DAT")
        new_name = prefix + "_MFILE.DAT"
        new_location_mfile_names.append(new_name)

    process_dir = Path(__file__).resolve().parent.parent
    input_file_dir = process_dir / "tests/regression/input_files"

    for count, file in enumerate(input_files):
        input_file_path = input_file_dir / file
        output_file_path = input_file_dir / mfile_names[count]
        subprocess.run(f"process -i {input_file_path}", shell=True)
        new_mfile_location = process_dir / "tracking" / new_location_mfile_names[count]
        shutil.move(output_file_path, new_mfile_location)


# def tracking():
# ^ will use the parsed args of track hash etc
# set up parse, parse args, then if sub parser is run call run_and_move_tracked_files, if other args then do
# the other fn which will be in here
if __name__ == "__main__":
    run_and_move_tracked_files()


# def main(args=None):
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-r", "--run")
#     parser.add_argument(
#         "-t",
#         "--track",
#     )
