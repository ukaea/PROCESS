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
# # Comparing multiple PROCESS solutions (MFiles)
#
# This tool plots the solution vectors (i.e. final values of optimisation parameters) for different runs of PROCESS.
# This allows visual comparisons of different solution points.
#
# It can use different intra-solution optimisation parameter normalisations
# (e.g. initial value, parameter range) and inter-solution normalisations (e.g. normalise to a certain solution).
#
# ### Known Limitations
#
# - The solution vectors (optimisation parameter values at the solution) currently
#   plotted are normalised to the initial point (from the `IN.DAT`) of each solution:
#   each element of the vector is the $x_{final}/x_{initial}$, the `xcmxxx` values in the `MFILE.DAT`.
#   This allows all optimisation parameters to be plotted on the same axis, showing the relative changes
#   from their initial values across multiple solutions.
# - Solutions being plotted together must also have the same optimisation parameters.
# - The solutions plotted in this example are fictitious.

# %%
# Reload Process each time (keep editable install up-to-date)
# %load_ext autoreload
# %autoreload 2

from pathlib import Path

from process.io.plot_solutions import (
    RunMetadata,
    plot_mfile_solutions,
    plot_mfile_solutions_constraints,
)

# %% [markdown]
# ## Plot single solution
#
# Plot a single solution, showing optimisation parameters normalised to their initial values.

# %%
data_dir = Path("data")
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
]

# Figure and dataframe returned for optional further modification
fig1, df1 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak solution 1",
)
df1

# %% [markdown]
# ## Plot two solutions
#
# Plot two MFILEs together, showing normalised values of the optimisation parameters at the solution points,
# as well as the objective function values.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
]

fig2, df2 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="2 large tokamak solutions",
)
df2

# %% [markdown]
# ## Plot one solution normalised to another
#
# Normalised differences, relative to the a given solution, can also be plotted:

# %%
fig3, df3 = plot_mfile_solutions(
    runs_metadata=runs_metadata,
    plot_title="Large tokamak 2 solution, relative to large tokamak 1",
    normalising_tag="large tokamak 1",
)
df3

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

fig4, df4 = plot_mfile_solutions(
    runs_metadata,
    "2 large tokamak solutions, relative to large tokamak 1",
    normalising_tag="large tokamak 1",
)
df4

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

fig5, df5 = plot_mfile_solutions(
    runs_metadata,
    "3 large tokamak solutions with RMS errors normalised to large tokamak 1",
    normalising_tag="large tokamak 1",
    rmse=True,
)
df5

# %% [markdown]
# ## Solutions normalised by range
#
# Use `nitvar` values instead; the solution optimisation parameters are normalised to the range of their upper and lower bounds.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
    RunMetadata(data_dir / "large_tokamak_3_MFILE.DAT", "large tokamak 3"),
    RunMetadata(data_dir / "large_tokamak_4_MFILE.DAT", "large tokamak 4"),
]

fig6, df6 = plot_mfile_solutions(
    runs_metadata,
    "4 large tokamak solutions normalised to the range of the optimisation parameters",
    normalisation_type="range",
)
df6

# %% [markdown]
# ## Actual values

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
    RunMetadata(data_dir / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
    RunMetadata(data_dir / "large_tokamak_3_MFILE.DAT", "large tokamak 3"),
    RunMetadata(data_dir / "large_tokamak_4_MFILE.DAT", "large tokamak 4"),
]

fig7, df7 = plot_mfile_solutions(
    runs_metadata,
    "Actual values of optimisation parameters for 4 large tokamak solutions",
    normalisation_type=None,
)
df7

# %% [markdown]
# ## Plot constraints
# Plot constraint values for 2 output mfiles.

# %%
runs_metadata = [
    RunMetadata(data_dir / "large_tokamak_nof_MFILE.DAT", "large tokamak"),
    RunMetadata(data_dir / "large_tokamak_nof_2_MFILE.DAT", "large tokamak 2"),
]

fig, eval_df = plot_mfile_solutions_constraints(
    runs_metadata=runs_metadata,
    title="Constraint values at large tokamak solution",
)

# %% [markdown]
# The constraint values are the actual constraint values, where a positive value corresponds to a satisfied constraint.
