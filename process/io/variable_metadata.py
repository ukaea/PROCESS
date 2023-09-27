from dataclasses import dataclass


@dataclass
class VariableMetadata:
    latex: str
    description: str


var_dicts = {
    "coe": VariableMetadata(latex=r"$COE$", description="Levelised COE"),
}
