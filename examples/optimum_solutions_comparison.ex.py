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
# # Optimum solutions comparison notebook
#
# %% [markdown]
# A Jupyter notebook to demonstrate how changing input parameters changes the optimum solution found by PROCESS.
#
# %% [markdown]
# <div class="alert alert-block alert-info">
# <b>NOTE</b> We run the examples in a temporary directory so all the inputs are copied there and the outputs contained there before the directory is removed when the example has finished running. This keeps the examples directory tidy and does not permanently modify any data files. The use of temporary directories is not needed for regular use of PROCESS.
# </div>

# %% [markdown]
# This notebook demonstrates how the optimum solution found by `PROCESS` changes as we vary input parameters. We will use the large tokamak example input file to do this. The figure of merit for this example is to minimise the major radius, `rmajor`.
#
# We use the functionality from `plot_solutions.py` and from `plot_proc.py` to demonstrate this.
#
# These tools plot the solution vectors (i.e. final values of optimisation parameters) for different runs of `PROCESS`. This allows visual comparisons of different solution points.
#
# It can use different intra-solution optimisation parameter normalisations (e.g. initial value, parameter range) and inter-solution normalisations (e.g. normalise to a certain solution).
#
# ### Known Limitations
#
# - The solution vectors (optimisation parameter values at the solution) currently plotted are normalised to the initial point (from the `IN.DAT`) of each solution: each element of the vector is the $x_{final}/x_{initial}$, the `xcmxxx` values in the `MFILE.DAT`. This allows all optimisation parameters to be plotted on the same axis, showing the relative changes from their initial values across multiple solutions.
# - Solutions being plotted together must also have the same optimisation parameters.
# - The solutions plotted in this example are fictitious.

# %% [markdown]
# ## Setup
#
# First we need to generate the necessary MFILES to be used in the rest of this notebook.
#

# %%
# %load_ext autoreload
# %autoreload 2

import shutil
import tempfile
from pathlib import Path

from process.core.io.plot.plot_solutions import (
    RunMetadata,
    plot_mfile_solutions,
)
from process.core.repository import get_process_root
from process.main import SingleRun

working_dir = Path.cwd()

# Define input file name relative to project dir, then copy to temp dir
script_dir = Path("__file__").parent.resolve()

data_dir = get_process_root() / "../examples/data/"
input_file = data_dir / "large_tokamak_IN.DAT"
input_file2 = data_dir / "large_tokamak_varied_min_net_electric_IN.DAT"

# %%
# Copy the file to avoid polluting the project directory with example files
temp_dir = tempfile.TemporaryDirectory()
working_dir = Path(temp_dir.name)
input_path = Path(temp_dir.name) / "large_tokamak_IN.DAT"
input_path2 = Path(temp_dir.name) / "large_tokamak_varied_min_net_electric_IN.DAT"
shutil.copy(input_file, input_path)
shutil.copy(input_file2, input_path2)

# %%
# Run process on these input files in a temporary directory
single_run = SingleRun(input_path.as_posix())
single_run.run()
single_run = SingleRun(input_path2.as_posix())
single_run.run()


# %% [markdown]
# # Plot single solution
#
# First we will look at the original large tokamak optimum solution. We will plot its solution, showing optimisation parameters normalised to their initial values.

# %%
large_tokamak_mfile = working_dir / "large_tokamak_MFILE.DAT"
# Plot the solution
runs_metadata = [
    RunMetadata(large_tokamak_mfile, "large tokamak"),
]

# Figure and dataframe returned for optional further modification
fig1, df1 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak solution",
)
df1

# %% [markdown]
# # Comparing optimum solutions
#
# Now we will see the effect that varying an input parameter in the large tokamak input file has on the optimum solution found.
#
# Here, the minimum allowable value for net electric power, `p_plant_electric_net_required_mw`, has been changed from 400MW to 200MW, and PROCESS has found a different optimum solution.
#
# We can plot the two MFILEs together, showing normalised values of the optimisation parameters at the solution points, as well as the objective function values.

# %%
large_tokamak_varied_min_net_electric_mfile = (
    working_dir / "large_tokamak_varied_min_net_electric_MFILE.DAT"
)
runs_metadata = [
    RunMetadata(large_tokamak_mfile, "original"),
    RunMetadata(
        large_tokamak_varied_min_net_electric_mfile,
        "changed min net electric",
    ),
]

fig2, df2 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="2 large tokamak solutions",
)
df2

# %% [markdown]
# We can compare the inequality constraint equations for both solutions.

# %%
import matplotlib.pyplot as plt

import process.core.io.mfile as mf
from process.core.io.plot.plot_proc import plot_inequality_constraint_equations

original_mfile = mf.MFile((large_tokamak_mfile).as_posix())
new_mfile = mf.MFile((large_tokamak_varied_min_net_electric_mfile).as_posix())

f, axs = plt.subplots(1, 2)
axs[0].set_position([0.0, 0.0, 1.2, 1.5])
axs[1].set_position([1.9, 0.0, 1.2, 1.5])
plot_inequality_constraint_equations(axis=axs[0], m_file=original_mfile, scan=-1)
plot_inequality_constraint_equations(axis=axs[1], m_file=new_mfile, scan=-1)
axs[0].set_title("Original optimum solution")
axs[1].set_title("New optimum solution when changing min net electric")
f.suptitle("Inequality Constraint Equations", y=1.6, x=1.4)

# %% [markdown]
# To have lower net electric, `PROCESS` has found a solution where:
# - the fusion power has dropped, therefore the neutron wall load has gone down
# - the toroidal field required has slightly dropped, therefore the case stress limits are lower as less current in the coils is needed

# %% [markdown]
# # Other solution comparison plots
#
# There are some other ways that you can compare solutions in `PROCESS`. We will demonstrate these for the same two MFILEs as above, but you can add in more MFILEs if you want.
#
# Here we refer to the original large tokamak file as `Large tokamak 1`, and the new solution obtained by varying `p_plant_electric_net_required_mw` as `Large tokamak 2`.

# %% [markdown]
# ## Plot one solution normalised to another
#
# Normalised differences, relative to the a given solution, can also be plotted.

# %%
runs_metadata = [
    RunMetadata(large_tokamak_mfile, "large tokamak 1"),
    RunMetadata(
        large_tokamak_varied_min_net_electric_mfile,
        "large tokamak 2",
    ),
]

fig3, df3 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak 2 solution, relative to large tokamak 1",
    normalising_tag="large tokamak 1",
)
df3

# %% [markdown]
# ## RMS Errors
#
# Plot RMS errors of multiple solutions relative to a reference solution.

# %%
fig5, df5 = plot_mfile_solutions(
    runs_metadata,
    "Large tokamak 2 solution with RMS errors normalised to large tokamak 1",
    normalising_tag="large tokamak 1",
    rmse=True,
)
df5

# %% [markdown]
# ## Solutions normalised by range
#
# Use `nitvar` values instead; the solution optimisation parameters are normalised to the range of their upper and lower bounds.

# %%
fig6, df6 = plot_mfile_solutions(
    runs_metadata,
    "Large tokamak 2 solution normalised to the range of the optimisation parameters",
    normalisation_type="range",
)
df6

# %% [markdown]
# ## Actual values

# %%
fig7, df7 = plot_mfile_solutions(
    runs_metadata,
    "Actual values of optimisation parameters for large tokamak 1 and 2 solutions",
    normalisation_type=None,
)
df7

# %%
temp_dir.cleanup()
