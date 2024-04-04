from pathlib import Path
from typing import Any, List
from enum import Enum
import pickle
import dataclasses
from copy import copy
import numpy as np
import jinja2

from process import fortran


class VariableTypes(str, Enum):
    PARAMETER = "Parameter"
    INPUT = "Input"
    OUTPUT = "Output"
    VARIABLE = "Variable"

    def __str__(self) -> str:
        return self.value


@dataclasses.dataclass
class FortranVariable:
    name: str
    module: str
    description: str
    typ: VariableTypes
    datatype: str
    initial: Any
    private: bool


def get_variables_and_modules(ford_project: Path):
    with open(ford_project, "rb") as ford_pickle:
        ford_project = pickle.load(ford_pickle)

    variables: List[FortranVariable] = []
    modules: List[str] = []

    for mod in ford_project.modules:
        for var in mod.variables:
            permission = "".join(var.permission).lower()
            modules.append(mod.name.lower())
            variables.append(
                FortranVariable(
                    var.name.lower(),
                    mod.name.lower(),
                    var.doc,
                    (
                        VariableTypes.PARAMETER
                        if var.parameter
                        else VariableTypes.VARIABLE
                    ),
                    var.vartype,
                    var.initial,
                    permission == "private",
                )
            )

    return variables, modules


def get_input_output_variables(variables: List[FortranVariable]):
    fortran_initial_values = {}
    for var in variables:
        if var.typ == VariableTypes.PARAMETER:
            continue

        if var.initial is not None:
            continue

        # cant get the value of private variables
        if var.private:
            continue

        try:
            fortran_initial_values[f"{var.module}.{var.name}"] = copy(
                getattr(getattr(fortran, var.module), var.name)
            )
        except AttributeError:
            continue

    fortran.init_module.init_all_module_vars()

    for var in variables:
        current_values_entry = f"{var.module}.{var.name}"
        if current_values_entry not in fortran_initial_values:
            continue

        try:
            new_value = getattr(getattr(fortran, var.module), var.name)
        except ValueError:
            continue

        try:
            is_input = np.any(new_value != fortran_initial_values[current_values_entry])
        except (TypeError, AttributeError):
            is_input = np.any(new_value != fortran_initial_values[current_values_entry])

        if is_input:
            var.typ = VariableTypes.INPUT
            var.initial = new_value
        else:
            var.typ = VariableTypes.OUTPUT

    return variables


if __name__ == "__main__":
    variables, modules = get_variables_and_modules(
        Path(__file__).resolve().parent.parent / "build/ford_project.pickle"
    )

    variables = get_input_output_variables(variables)

    mods = {module: [] for module in modules}
    for var in variables:
        mods[var.module].append(var)

    loader = jinja2.FileSystemLoader(searchpath=Path(__file__).resolve().parent)
    env = jinja2.Environment(loader=loader)

    vardes_template = env.get_template("vardes.jinja2")

    with open(
        Path(__file__).resolve().parent / "../documentation/proc-pages/io/vardes.md",
        "w",
    ) as f:
        f.write(vardes_template.render(mods=mods))
