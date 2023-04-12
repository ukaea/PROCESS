# %% [markdown]
# # `plot_solutions` Solution Comparison Tool
#
# This tool plots the solution vectors (i.e. final values of optimisation parameters) for different runs of PROCESS. This allows visual comparisons of different solution points.

# %%
# Reload Process each time (keep editable install up-to-date)
# %load_ext autoreload
# %autoreload 2

from process.io.plot_solutions import RunMetadata, plot_mfile_solutions
from pathlib import Path

# %% [markdown]
# ## Plot single solution
#
# Plot a single solution, showing optimisation parameters normalised to their initial values.

# %%
data_dir = Path("data")
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
]

plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak solution 1",
)

# %% [markdown]
# ## Plot two solutions
#
# Plot two MFILEs together, showing normalised values of the optimisation parameters at the solution points, as well as the objective function values.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
]

plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="2 large tokamak solutions",
)

# %% [markdown]
# ## Plot one solution normalised to another
#
# Normalised differences, relative to the a given solution, can also be plotted:

# %%
plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak 2 solution, relative to large tokamak 1",
    normalise=True,
    normalising_tag="large tokamak 1",
)

# %% [markdown]
# ## Plot multiple solutions normalised by one
#
# Plot two MFILEs, normalised by a third MFILE.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
    RunMetadata(data_dir / "large_tokamak_3_MFILE.DAT", "large tokamak 3"),
]

plot_mfile_solutions(
    runs_metadata,
    "2 large tokamak solutions, relative to large tokamak 1",
    normalise=True,
    normalising_tag="large tokamak 1",
)

# %% [markdown]
# ## RMS Errors
#
# Plot RMS errors of multiple solutions relative to a reference solution.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
    RunMetadata(data_dir / "large_tokamak_3_MFILE.DAT", "large tokamak 3"),
    RunMetadata(data_dir / "large_tokamak_4_MFILE.DAT", "large tokamak 4"),
]

plot_mfile_solutions(
    runs_metadata,
    "4 large tokamak solutions with RMS errors normalised to large tokamak 2",
    normalise=False,
    normalising_tag="large tokamak 2",
    rmse=True,
)
