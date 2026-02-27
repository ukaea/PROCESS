"""
Code to read from a PROCESS MFILE and write values into a csv

Input files:
mfile (default MFILE.DAT) as output from PROCESS
variable list (default mfile_to_csv_vars.json) as defined by user

Instructions:
- command line call: python mfile_to_csv.py -f </path/to/mfile.dat> -v </path/to/varfile.json>

Output file:
.csv will be saved to the directory of the input file
"""

import json
from collections.abc import Sequence
from pathlib import Path, PurePath

import numpy as np

from process.core.io.mfile.mfile import MFile

default_vars = (
    "minmax",
    "p_hcd_injected_max",
    "p_plant_electric_net_required_mw",
    "ripple_b_tf_plasma_edge_max",
    "t_burn_min",
    "alstroh",
    "sig_tf_wp_max",
    "dx_tf_turn_steel",
    "f_j_cs_start_pulse_end_flat_top",
    "alstroh",
    "rmajor",
    "dr_tf_inboard",
    "dr_cs",
    "c_tf_turn",
    "dr_tf_wp_with_insulation",
    "dr_cryostat",
    "dr_shld_outboard",
    "dz_divertor",
    "rmajor",
)


def get_vars(vfile="mfile_to_csv_vars.json"):
    """Returns variable names from identified file.

    Parameters
    ----------
    args : string
        input JSON filename
    vfile :
         (Default value = "mfile_to_csv_vars.json")

    Returns
    -------
    list
        variable names
    """
    print("Fetching list of variables from", vfile)

    return json.loads(Path(vfile).read_text())["vars"]


def read_mfile(mfilename="MFILE.DAT", variables=None):
    """Returns specified variable values from identified file.

    Parameters
    ----------
    args : string, list
        input filename, variable names
    mfilename :
         (Default value = "MFILE.DAT")
    variables :
         (Default value = None)

    Returns
    -------
    list of tuples
        variable descriptions, names, and values
    """
    if variables is None:
        variables = []
    print("Reading from MFILE:", mfilename)

    m_file = MFile(mfilename)

    output_vars = []

    # for each variable named in the input varfile, get the description and data value
    for var_name in variables:
        if var_name not in m_file.data:
            print(f"Variable '{var_name}' not in MFILE. Skipping and moving on...")
        else:
            # In case of a file containing multiple scans, (scan = -1) uses the last scan value
            output_vars.append((
                m_file.data[var_name].var_description,
                var_name,
                m_file.get(var_name, scan=-1),
            ))

    return output_vars


def get_savenamepath(mfilename="MFILE.DAT"):
    """Returns path/filename.csv for file saving.

    Parameters
    ----------
    args : string
        input filename
    mfilename :
         (Default value = "MFILE.DAT")

    Returns
    -------
    pathlib.PurePosixPath
        output filename
    """

    # Either save it locally or output the csv file to the directory of the input file
    dirname = Path.cwd() if mfilename == "MFILE.DAT" else PurePath(mfilename).parent

    csv_filename = PurePath(mfilename).stem
    return PurePath(dirname, csv_filename + ".csv")


def write_to_csv(csv_outfile, output_data=None):
    """Write to csv file.

    Parameters
    ----------
    args : string, list of tuples
        input filename, variable data
    csv_outfile :

    output_data :
         (Default value = None)
    """
    print("Writing to csv file:", csv_outfile)
    np.savetxt(
        csv_outfile,
        output_data or [],
        fmt="%.5e",
        delimiter=",",
        header="Description, Varname, Value",
        footer="",
        comments="",
    )


def to_csv(mfile, variables: Sequence[str] | str = default_vars):
    """Extract certain variables from an MFILE.DAT and output to CSV.

    Parameters
    ----------
    mfile:
        Mfile to convert
    variables:
        variable file with variables to extract
    """
    write_to_csv(
        get_savenamepath(mfile),
        read_mfile(
            mfile, get_vars(variables) if isinstance(variables, str) else variables
        ),
    )

    # write final line to screen
    print("Complete.")
