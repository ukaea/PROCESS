import re
from pathlib import Path

import process.data_structure as data_structure
import process.core.solver.iteration_variables as iteration_variables
from process.main import SingleRun


def get_mfile_initial_ixc_values(file_path: Path):
    """Initialise the input file and obtain the initial values of the iteration variables

    Parameters
    ----------
    file_path :
        The path to the MFile to get the initial iteration variable values from.

    Notes
    -----
    This method initialises a SingleRun. At present, this involves mutating the global
    data structure so it is not safe to run this method during a PROCESS run.
    """
    SingleRun(file_path.as_posix())
    iteration_variables.load_iteration_variables()

    iteration_variable_names = []
    iteration_variable_values = []

    for i in range(data_structure.numerics.nvar):
        ivar = data_structure.numerics.ixc[i].item()

        itv = iteration_variables.ITERATION_VARIABLES[ivar]

        iteration_variable_names.append(itv.name)
        if array := re.match(r"(\w+)\(([0-9]+)\)", itv.name):
            var_name = array.group(1)
            index = array.group(2)
            iteration_variable_values.append(
                getattr(itv.module, var_name)[int(index) - 1]
            )
        else:
            iteration_variable_values.append(getattr(itv.module, itv.name))

    return iteration_variable_names, iteration_variable_values
