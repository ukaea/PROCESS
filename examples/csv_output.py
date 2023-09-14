# %% [markdown]
# # Output to csv
#
# Routine to read from a PROCESS MFILE and write specified values into a csv.
#
# Input files:
# - MFILE.DAT as output from PROCESS
# - .json variable list as defined by user (defaults to local `mfile_to_csv_vars.json`)
#
# Instructions:
# - from command line: `python mfile_to_csv.py -f </path/to/mfile.dat> -v </path/to/varfile.json>`
# - from this Jupyter notebook: run the cell below
#
# Output file:
# - .csv will be saved to the directory of the input file

# %%
from pathlib import Path
from process.io import mfile_to_csv

# Project directory for example result file and default .json list;
# not needed if you replace both target filepaths below.
proj_dir = Path.cwd().parent

# Replace this path/to/MFILE.DAT with your target file:
mfilename = proj_dir / "examples/csv_output_large_tokamak_MFILE.DAT"

# Either replace this with your own path/to/file.json target,
# or add your required variables into the identified file:
varfilename = proj_dir / "process/io/mfile_to_csv_vars.json"
# This routine attempts to find every variable in the given list and
# writes the variable name, description and value to the output csv.
# Any listed variable that isn't in that MFILE will be skipped.

# call to function:
mfile_to_csv.main(args=["-f", str(mfilename), "-v", str(varfilename)])


# %%
