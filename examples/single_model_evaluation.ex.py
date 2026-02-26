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
# # Evaluating a single PROCESS model
# When understanding or investigating an individual model within Process, it can be useful to run the model in isolation and plot some responses. This is done here to investigate the effect of tungsten impurity concentration on radiated power and power incident on the divertor.

# %%
import matplotlib.pyplot as plt
import numpy as np

import process
from process.main import SingleRun

# %% [markdown]
# ## Set up
# First, inspect a variable to check its uninitialised value:

# %%
print(
    f"p_plasma_separatrix_mw = {process.data_structure.physics_variables.p_plasma_separatrix_mw}"
)

# %% [markdown]
# In order to initialise all variables in Process with their values at a given point (design parameter vector), run an evaluation input file (one with no optimisation) to initialise values in all models. The "large tokamak" regression test solution is used here.

# %%
from process.repository import get_process_root

data_dir = get_process_root() / "../examples/data/"
single_run = SingleRun((data_dir / "large_tokamak_eval_IN.DAT").as_posix())
single_run.run()


# %%
# Kernel can crash when running physics without correctly initialised values
# Doesn't crash after running a once-through
# Print initial values of interest
def print_values():
    print(
        f"W frac = {process.data_structure.impurity_radiation_module.f_nd_impurity_electron_array[13]:.3e}"
    )
    print(
        f"p_plasma_rad_mw = {process.data_structure.physics_variables.p_plasma_rad_mw:.3e}"
    )
    print(
        f"p_plasma_separatrix_mw = {process.data_structure.physics_variables.p_plasma_separatrix_mw:.3e}"
    )


print_values()

# %% [markdown]
# Now try increasing the tungsten impurity fraction to see if there's a change in the divertor power.

# %%
process.data_structure.impurity_radiation_module.f_nd_impurity_electron_array[13] = (
    5.0e-5
)
single_run.models.physics.physics()
print_values()

# %% [markdown]
# With a higher W impurity fraction, the radiated power has increased and the power incident on the divertor has decreased.
#
# ## Parameter study of W impurity
# Now investigate effect of varying W impurity on impurity radiation power, divertor power and constraint 15 (L-H threshold constraint).

# %%
from process.constraints import ConstraintManager


def run_impurities(w_imp_fracs):
    """Calculate responses to W impurities."""
    n = w_imp_fracs.shape[0]
    p_plasma_rad_mw = np.empty(n)
    p_plasma_separatrix_mw = np.empty(n)
    p_l_h_threshold_mw = np.empty(n)
    con15 = np.empty(n)

    # Loop over W impurity values, evaluate model and store responses at each point
    for i, imp_frac in enumerate(w_imp_fracs):
        # Set W impurity fraction, then run physics model
        process.data_structure.impurity_radiation_module.f_nd_impurity_electron_array[
            13
        ] = imp_frac
        single_run.models.physics.physics()

        # Evaluate constraint equation 15 (L-H threshold constraint)
        con15_value = ConstraintManager.evaluate_constraint(15).normalised_residual

        # Need to copy values
        p_plasma_rad_mw[i] = (
            process.data_structure.physics_variables.p_plasma_rad_mw.item()
        )
        p_plasma_separatrix_mw[i] = (
            process.data_structure.physics_variables.p_plasma_separatrix_mw.item()
        )
        p_l_h_threshold_mw[i] = (
            process.data_structure.physics_variables.p_l_h_threshold_mw.item()
        )
        # Need to flip sign of constraint so negative means violated
        con15[i] = -con15_value

    return p_plasma_rad_mw, p_plasma_separatrix_mw, p_l_h_threshold_mw, con15


# %%
# %matplotlib inline
# Run W impurity parameter study
w_imp_fracs = np.linspace(1.0e-6, 1.0e-4, 50)
p_plasma_rad_mw, p_plasma_separatrix_mw, p_l_h_threshold_mw, con15 = run_impurities(
    w_imp_fracs
)

fig, ax = plt.subplots()
ax.scatter(w_imp_fracs, p_plasma_rad_mw, label="p_plasma_rad_mw")
ax.scatter(w_imp_fracs, p_plasma_separatrix_mw, label="p_plasma_separatrix_mw")
ax.plot(w_imp_fracs, p_l_h_threshold_mw, "r", label="p_l_h_threshold_mw")
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
ax.set_title("W impurity fraction against radiated and divertor power")
ax.set_xlabel("W impurity fraction")
ax.set_ylabel("Power (MW)")
ax.legend()

# %% [markdown]
# How does the L-H threshold constraint vary?

# %%
# %matplotlib inline
fig, ax = plt.subplots()
ax.scatter(w_imp_fracs, con15, label="con15")
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
ax.set_title("W impurity fraction against L-H threshold constraint (15)")
ax.set_xlabel("W impurity fraction")
ax.set_ylabel("L-H threshold constraint value")
ax.hlines(0.0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors="r")
ax.annotate("Satisfied", (0.0, 0.1))
ax.annotate("Violated", (0.0, -0.15))

# %% [markdown]
# The constraint becomes violated for W fraction values $> 6\times10^{-5}$.
#
# This can easily be modified to investigate behaviour of any model in Process in isolation, without running other models or optimising.
