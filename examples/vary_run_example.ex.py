# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,ipynb
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Demonstration of VaryRun

# %% [markdown]
# A Jupyter notebook to demonstrate usage of `VaryRun` in `PROCESS`.
#
# `VaryRun` is a tool which takes an input file that does not converge and varies the initial values of the
# iteration variables, within a tolerance, to find an initial point that converges, and creates
# a new input file using these values.
#
# `VaryRun` requires a `.conf` file which specifies certain parameters needed for `VaryRun`.
# In this file, you specify the original input file, the maximum number of iterations to be performed, and a factor within which the iteration variables are changed.
#
# If `VaryRun` is able to find a new initial point within the maximum number of iterations, it produces a new input file, called `IN.DAT`, in the same directory as your initial input file. This new file will now converge when you run `PROCESS`.
#
# ## Setup
# We run the examples in a temporary directory so all the inputs are copied there and the outputs contained there before the directory is removed when the example has finished
# running. This keeps the examples directory tidy and does not permanently modify any data files. This use of temporary directories is not needed for regular use of `PROCESS`.

# %% [markdown]
# # VaryRun setup
#
# VaryRun requires a `.conf` file in order to run. An example `.conf` file is displayed below.
#
# - `WDIR` specifies the path to the working directory where `PROCESS` is being run, here is is `.` to run in the directory of the `.conf` file, which is in `examples/data`
# - `ORIGINAL_IN_DAT` is providing the name of the original input file
# - `NITER` is specifying the maximum number of iterations to perform
# - `SEED` is specifying a random number generator. Here it is set to 5 so the results of running `VaryRun` on the original input file will always be the same
# - `FACTOR` is specifying a factor within which the iteration variables will change. Here it is set to 1.5, so this means the iteration variables can be varied within 50% of their initial values in the original input file

# %% [markdown]
# ```ini
# * CONFIG FILE FOR RUNNING PROCESS WITH MODIFIED IN.DAT
# *
#
# * Path to working directory in which PROCESS is run.
# WDIR = .
#
# * original IN.DAT name (should not be called IN.DAT!)
# ORIGINAL_IN_DAT = ORIGINAL_IN.DAT
#
# * Max no. iterations
# NITER = 30
#
# * integer seed for random number generator; use None for random seed
# SEED = 5
#
# * factor within which the iteration variables are changed
# FACTOR = 1.5
# ```

# %% [markdown]
# ## Run `VaryRun`
#
# Run `PROCESS` on an input file using the `VaryRun` class. The initial input file, `large_tokamak_varyrun_IN.DAT` does not converge. `VaryRun` will vary the initial values of the iteration variables to find an initial point that converges, and will create a new input file, `IN.DAT`, with these new values.

# %%
# %load_ext autoreload
# %autoreload 2

import os
from pathlib import Path
from shutil import copy

from functions_for_examples import copy_to_temp_dir, get_initial_values

from process.main import VaryRun

# Define project root dir; when running a notebook, the cwd is the dir the notebook is in
PROJ_DIR = Path.cwd().parent

# Path to .conf file
script_dir = Path("__file__").parent.resolve()
conf_file = script_dir / "data/run_process.conf"
temp_dir, temp_input_path, temp_dir_path = copy_to_temp_dir(conf_file, PROJ_DIR)

# .conf file relies on a separate input file too; copy this as well
# TODO This double input file requirement needs to be removed
input_file = script_dir / "data/large_tokamak_varyrun_IN.DAT"
copy(PROJ_DIR / input_file, temp_dir.name)


# VaryRun uses process_config.py, which changes the current working directory
# via os.chdir() to the temporary dir. Apart from being bad practice, once the
# temp dir is removed, this causes Path.cwd() (as used in plot_scans.py) to
# throw an exception when trying to return the (now deleted) CWD. Hence it
# needs to be set back after VaryRun()
# TODO Remove the os.chdir() from VaryRun
cwd = Path.cwd()

vary_run = VaryRun(temp_input_path.as_posix())
vary_run.run()
os.chdir(cwd)


# Get the initial values from the original input file
iteration_variable_names, original_iteration_variable_values = get_initial_values(
    input_file
)

# Get the initial values from the new input file produced by VaryRun
path_to_new_input = (temp_dir_path / "IN.DAT").as_posix()
_, updated_iteration_variable_values = get_initial_values(path_to_new_input)

# %% [markdown]
# ## Compare iteration variable values
#
# `VaryRun` has changed the values of some iteration variables.

# %%
import pandas as pd

# Use pandas to display the values of the iteration variables before and after running VaryRun
df = pd.DataFrame({
    "Iteration variable names": iteration_variable_names,
    "Original values": original_iteration_variable_values,
    "Updated values": updated_iteration_variable_values,
})
df

# %%
temp_dir.cleanup()
