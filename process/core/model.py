import abc
from dataclasses import dataclass, fields

from process.data_structure.blanket_variables import BlanketData
from process.data_structure.build_variables import BuildData
from process.data_structure.buildings_variables import BuildingsData
from process.data_structure.ccfe_hcpb_variables import CCFEHCPBData
from process.data_structure.constraint_variables import ConstraintData
from process.data_structure.cost_2015_variables import Cost2015Data
from process.data_structure.cost_variables import CostData
from process.data_structure.cs_fatigue_variables import CSFatigueData
from process.data_structure.current_drive_variables import CurrentDriveData
from process.data_structure.dcll_variables import DCLLData
from process.data_structure.first_wall_variables import FirstWallData
from process.data_structure.fwbs_variables import FWBSData
from process.data_structure.heat_transport_variables import HeatTransportData
from process.data_structure.primary_pumping_variables import PrimaryPumpingData
from process.data_structure.pulse_variables import PulseData
from process.data_structure.reinke_variables import ReinkeData
from process.data_structure.structure_variables import StructureData
from process.data_structure.times_variables import TimesData
from process.data_structure.vacuum_variables import VacuumData
from process.data_structure.water_usage_variables import WaterUseData

initialise_later = object()


@dataclass(kw_only=True)
class DataStructure:
    water_use: WaterUseData = initialise_later
    costs_2015: Cost2015Data = initialise_later
    cs_fatigue: CSFatigueData = initialise_later
    vacuum: VacuumData = initialise_later
    costs: CostData = initialise_later
    first_wall: FirstWallData = initialise_later
    fwbs: FWBSData = initialise_later
    blanket: BlanketData = initialise_later
    structure: StructureData = initialise_later
    times: TimesData = initialise_later
    reinke: ReinkeData = initialise_later
    ccfe_hcpb: CCFEHCPBData = initialise_later
    pulse: PulseData = initialise_later
    build: BuildData = initialise_later
    primary_pumping: PrimaryPumpingData = initialise_later
    buildings: BuildingsData = initialise_later
    constraints: ConstraintData = initialise_later
    dcll: DCLLData = initialise_later
    current_drive: CurrentDriveData = initialise_later
    heat_transport: HeatTransportData = initialise_later

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
