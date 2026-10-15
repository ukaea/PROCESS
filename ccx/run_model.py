import process
from pathlib import Path
from process.core.data_structure.base import DataStructure
from process.core.log import logging_model_handler, show_errors
from process.core.process_output import OutputFileManager, oheadr
from process.core.data_structure import obsolete_vars as ov
from process.main import setup_loggers
from process.core.solver import iteration_variables
import process.core.init as init
from process.core.input import parse_input_file
from process.core.init import set_active_constraints, set_device_type, check_process, run_summary
from process.core import constants, process_output
from process.models.physics.bootstrap_current import (
    PlasmaBootstrapCurrent,
    SauterBootstrapCurrent,
)
from process.models.physics.confinement_time import PlasmaConfinementTime
from process.models.physics.current_drive import (
    CurrentDrive,
    ElectronBernstein,
    ElectronCyclotron,
    IonCyclotron,
    LowerHybrid,
    NeutralBeam,
)
from process.models.physics.density_limit import PlasmaDensityLimit
from process.models.physics.exhaust import PlasmaExhaust
from process.models.physics.impurity_radiation import (
    initialise_imprad,
)
from process.models.physics.l_h_transition import PlasmaConfinementTransition
from process.models.physics.physics import (
    Physics,
    PlasmaBeta,
    PlasmaInductance,
)
from process.models.physics.plasma_current import (
    PlasmaCurrent,
    PlasmaDiamagneticCurrent,
)
from process.models.physics.plasma_fields import PlasmaFields
from process.models.physics.plasma_geometry import PlasmaGeom
from process.models.physics.fusion_reactions import BeamReactions, PlasmaReactions
from process.models.physics.plasma_profiles import PlasmaProfile
from process.models.physics.profiles import (
    ElectronDensityProfile,
    ElectronTemperatureProfile,
    IonTemperatureProfile,
)
from process.models.physics.scrape_off_layer import ScrapeOffLayer
from process.main import SingleRun

class Models:
    def __init__(self, data: DataStructure):
        self.data = data

        self.plasma_geom = PlasmaGeom()
 
        self.ne_profile = ElectronDensityProfile()
        self.te_profile = ElectronTemperatureProfile()
        self.ti_profile = IonTemperatureProfile(
            electron_temperature_profile=self.te_profile
        )
        self.plasma_profile = PlasmaProfile(
            ne_profile=self.ne_profile,
            te_profile=self.te_profile,
            ti_profile=self.ti_profile,
        )

        self.neutral_beam = NeutralBeam(plasma_profile=self.plasma_profile)
        self.electron_cyclotron = ElectronCyclotron(plasma_profile=self.plasma_profile)
        self.lower_hybrid = LowerHybrid(plasma_profile=self.plasma_profile)
        self.current_drive = CurrentDrive(
            plasma_profile=self.plasma_profile,
            electron_cyclotron=self.electron_cyclotron,
            ion_cyclotron=IonCyclotron(plasma_profile=self.plasma_profile),
            lower_hybrid=self.lower_hybrid,
            neutral_beam=self.neutral_beam,
            electron_bernstein=ElectronBernstein(plasma_profile=self.plasma_profile),
        )
        self.plasma_beta = PlasmaBeta()
        self.plasma_inductance = PlasmaInductance()
        self.plasma_density_limit = PlasmaDensityLimit()
        self.plasma_exhaust = PlasmaExhaust()
        self.sauter_bootstrap_current = SauterBootstrapCurrent()
        self.plasma_bootstrap_current = PlasmaBootstrapCurrent(
            plasma_profile=self.plasma_profile,
            sauter_bootstrap=self.sauter_bootstrap_current,
        )
        self.plasma_confinement = PlasmaConfinementTime()
        self.plasma_transition = PlasmaConfinementTransition()
        self.plasma_current = PlasmaCurrent()
        self.plasma_fields = PlasmaFields()
        self.plasma_dia_current = PlasmaDiamagneticCurrent()
        self.scrape_off_layer = ScrapeOffLayer()
        self.plasma_reactions = PlasmaReactions(plasma_profile=self.plasma_profile)
        self.beam_reactions = BeamReactions()
        self.physics = Physics(
            plasma_profile=self.plasma_profile,
            current_drive=self.current_drive,
            plasma_beta=self.plasma_beta,
            plasma_inductance=self.plasma_inductance,
            plasma_density_limit=self.plasma_density_limit,
            plasma_exhaust=self.plasma_exhaust,
            plasma_bootstrap_current=self.plasma_bootstrap_current,
            plasma_confinement=self.plasma_confinement,
            plasma_transition=self.plasma_transition,
            plasma_current=self.plasma_current,
            plasma_fields=self.plasma_fields,
            plasma_dia_current=self.plasma_dia_current,
            plasma_geometry=self.plasma_geom,
            scrape_off_layer=self.scrape_off_layer,
            plasma_reactions=self.plasma_reactions,
            beam_reactions=self.beam_reactions,
        )

        self.setup_data_structure()

    @property
    def models(self):
        """Set up the models"""
        # At the moment, this property just returns models
        # that implement the Model interface.
        # Eventually every Model will comply and then
        # this method can be used as the caller/outputter!
        return (
            self.physics,
            self.plasma_geom,
            self.plasma_confinement,
            self.plasma_beta,
            self.current_drive,
            self.neutral_beam,
            self.plasma_density_limit,
            self.plasma_profile,
            self.plasma_dia_current,
            self.plasma_bootstrap_current,
            self.plasma_exhaust,
            self.plasma_current,
            self.scrape_off_layer,
            self.plasma_inductance,
            self.ne_profile,
            self.te_profile,
            self.ti_profile,
            self.plasma_reactions,
            self.beam_reactions,
            self.plasma_fields,
            self.sauter_bootstrap_current,
            self.plasma_transition,
            self.electron_cyclotron,
            self.lower_hybrid,
        )

    def setup_data_structure(self):
        """Set up the data structure"""
        # This Models class should be replaced with a dataclass so we can
        # iterate over the `fields`.
        # This can be a disgusting temporary measure :(
        for model in self.models:
            model.data = self.data


class RunModel(SingleRun):


    def run(self):
        """Evaluate the models once, then write the output files."""
        self.run_models()
        self.write_model_output()
        self.finish()
        self.append_input()

    def run_models(self):
        """Run the tokamak models once.

        This does not start VMCON or the equality-constraint solver.
        Geometry is evaluated first because the physics model reads its results.
        """
        self.models.plasma_geom.run()
        self.models.physics.run()

    def write_model_output(self):
        """Write the results of the models run above."""
        self.models.plasma_geom.output()
        self.models.physics.output()

