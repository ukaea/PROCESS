"""Run the tracked files and move into tracking/ directory."""

import os
from pathlib import Path
import subprocess
import shutil

# Stellarator doesn't converge so don't want to track it
exclusions = ["stellarator"]

path_to_tracking_files = "tests/regression/input_files"
all_files = os.listdir(path_to_tracking_files)
file_names = all_files.copy()

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

print(f"new mfile names = {new_location_mfile_names}")

process_dir = Path(__file__).resolve().parent.parent
input_file_dir = process_dir / "tests/regression/input_files"

for count, file in enumerate(input_files):
    input_file_path = input_file_dir / file
    output_file_path = input_file_dir / mfile_names[count]
    subprocess.run(f"process -i {input_file_path}", shell=True)
    new_mfile_location = process_dir / "tracking" / new_location_mfile_names[count]
    shutil.move(output_file_path, new_mfile_location)
