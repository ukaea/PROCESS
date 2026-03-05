import abc
from dataclasses import dataclass, fields

from process.data_structure.water_usage_variables import WaterUseData

initialise_later = object()


@dataclass(kw_only=True)
class DataStructure:
    water_use: WaterUseData = initialise_later

    def __post_init__(self):
        for f in fields(self):
            if getattr(self, f.name) is initialise_later:
                setattr(self, f.name, f.type())


class Model(abc.ABC):
    data: DataStructure

    @abc.abstractmethod
    def run(self) -> None:
        """Run the model.

        The run method is resposible for 'running' the model, ensuring it updates the data
        structure with variables that subsequent models will require.
        """

    @abc.abstractmethod
    def output(self) -> None:
        """Output model data.

        This method will always be called after run method and should output the model data to the
        output files.
        """
