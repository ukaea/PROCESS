import abc
from dataclasses import dataclass, fields

from process.data_structure.cost_2015_variables import Cost2015Data
from process.data_structure.cs_fatigue_variables import CSFatigueData
from process.data_structure.vacuum_variables import VacuumData
from process.data_structure.water_usage_variables import WaterUseData

initialise_later = object()


@dataclass(kw_only=True)
class DataStructure:
    water_use: WaterUseData = initialise_later
    costs_2015: Cost2015Data = initialise_later
    cs_fatigue: CSFatigueData = initialise_later
    vacuum: VacuumData = initialise_later

    def __post_init__(self):
        for f in fields(self):
            if getattr(self, f.name) is initialise_later:
                setattr(self, f.name, f.type())


class Model(abc.ABC):
    data: DataStructure

    @abc.abstractmethod
    def run(self) -> None:
        """Run the model.

        The run method is responsible for 'running' the model, ensuring it updates the data
        structure with variables that subsequent models will require.
        """

    @abc.abstractmethod
    def output(self) -> None:
        """Output model data.

        This method will always be called after run method and should output the model data to the
        output files.
        """
