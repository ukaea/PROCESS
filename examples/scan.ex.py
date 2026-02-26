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

# %% [markdown] slideshow={"slide_type": "slide"}
# # Running and visualising a PROCESS scan
#
# Perform a parameter scan for a given input file and plot the results.
#
# ## Scan details
#
# The input file is a scan-enabled version of the large tokamak `IN.DAT`, as found in the `tests` directory. The scan-relevant values are:
# ```
# nsweep = 17 * b_tf_inboard_max, maximum peak toroidal field (T) (constraint equation 25)
# isweep = 6
# sweep = 10.5, 10.4, 10.3, 10.2, 10.1, 10.0
# ```
#
# - `nsweep`: integer denoting the variable to scan (see `scan_module` for options). Here `17` corresponds to `b_tf_inboard_max` being scanned
# - `isweep`: the number of scan points to run
# - `sweep`: array of values for the scanned variable to take; one for each run. Should be of length `isweep`

# %% slideshow={"slide_type": "subslide"}
from process.core.repository import get_process_root
from process.main import SingleRun

data_dir = get_process_root() / "../examples/data/"
input_name = data_dir / "scan_example_file_IN.DAT"
# Perform a SingleRun on a scan-enabled input file
single_run = SingleRun(str(input_name), solver="vmcon_bounded")
single_run.run()

# %% [markdown] slideshow={"slide_type": "slide"}
# ## Plot scan results
# Use `plot_scans.py` to plot the resulting `MFILE.DAT`.

# %% slideshow={"slide_type": "subslide"}
# %matplotlib inline
from process.core.io import plot_scans

# Define working directory relative to project dir and input file name
mfile_name = data_dir / "scan_example_file_MFILE.DAT"
output_dir = data_dir

plot_scans.main(
    args=[
        "-f",
        str(mfile_name),
        "-yv",
        "b_plasma_toroidal_on_axis rmajor p_plant_electric_net_mw p_fusion_total_mw capcost",
        "--outputdir",
        str(output_dir),
    ]
)
