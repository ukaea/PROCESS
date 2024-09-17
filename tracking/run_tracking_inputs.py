"""TODO"""

# Exclusion list contains stellarator
# Stellarator doesn't converge so don't want to track it
import os

from pathlib import Path
import subprocess
import shutil


# def run_and_move_tracked_files():
exclusions = ["stellarator"]

# want to find files in dir tests/regression/input_files
# then just want IN.DAT files
# then exclude any containing the word from exclusion list (stellarator)
# then want to run process on that file
# then move the created MFILE.DAT to the folder tracking/ (move from tests/regression/input_files)

path_to_tracking_files = "tests/regression/input_files"
all_files = os.listdir(path_to_tracking_files)
file_names = all_files.copy()
# print(f"all files = {all_files}")

input_files = []

for file in all_files:
    # print(file)
    for excl in exclusions:
        if excl not in file and file.endswith("IN.DAT"):
            # print(f"file {file} is an input file")
            input_files.append(file)
# print(f"input files = {input_files}")

mfile_names = []
for in_dat in input_files:
    # print(f"currently looking at file {in_dat}")
    file_prefix = in_dat.removesuffix("IN.DAT")
    # print(f"file prefix = {file_prefix}")
    mfile_name = file_prefix + "MFILE.DAT"
    # print(f"new mfile name = {mfile_name}")
    mfile_names.append(mfile_name)
# print(f"input files = {input_files}")
# print(f"new mfile names = {mfile_names}")

process_dir = Path(__file__).resolve().parent.parent
# print(f"process dir = {process_dir}")
input_file_dir = process_dir / "tests/regression/input_files"
# print(input_file_dir)
# new_location = script_dir / "tests/regresstion"
for count, file in enumerate(input_files):
    # print("in loop")
    # print(f"current input file = {file}")
    # print(f"matching mfile = {mfile_names[count]}")
    input_file_path = input_file_dir / file
    output_file_path = input_file_dir / mfile_names[count]
    # print(f"currently looking at {file}")
    # print(f"input file path = {input_file_path}")
    # print(f"output file path = {output_file_path}")
    subprocess.run(f"process -i {input_file_path}", shell=True)
    new_mfile_location = process_dir / "tracking" / mfile_names[count]
    # print(f"new mfile loc = {new_mfile_location}")
    shutil.move(output_file_path, new_mfile_location)


# if __name__ == "__main__":
#     print(f"Running script: {run_and_move_tracked_files()}")
