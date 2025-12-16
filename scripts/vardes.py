import dataclasses
import pickle
from enum import Enum
from pathlib import Path
from typing import Any

import jinja2


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

    variables: list[FortranVariable] = []
    modules: list[str] = []

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


if __name__ == "__main__":
    mods = {"vardes will be fixed in https://github.com/ukaea/PROCESS/issues/3834": []}

    loader = jinja2.FileSystemLoader(searchpath=Path(__file__).resolve().parent)
    env = jinja2.Environment(loader=loader)

    vardes_template = env.get_template("vardes.jinja2")

    (Path(__file__).resolve().parent / "../documentation/io/vardes.md").write_text(
        vardes_template.render(mods=mods)
    )
