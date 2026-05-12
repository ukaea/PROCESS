from process import data_structure
from process.core.log import logging_model_handler
from process.data_structure.blanket_variables import BlktModelTypes
from process.models.tfcoil.base import TFConductorModel
from process.models.tfcoil.superconducting import (
    SuperconductingTFTurnType,
)


def write(models, data, _outfile):
    """Write the results to the main output file (OUT.DAT).

    Write the program results to a file, in a tidy format.

    Parameters
    ----------
    models : process.main.Models
        physics and engineering model objects
    _outfile : int
        Fortran output unit identifier

    """
    # ensure we are capturing warnings that occur in the 'output' stage as these are warnings
    # that occur at our solution point. So we clear existing warnings
    logging_model_handler.start_capturing()
    logging_model_handler.clear_logs()

    # Call stellarator output routine instead if relevant
    if data_structure.stellarator_variables.istell != 0:
        models.stellarator.output()
        return

    #  Call IFE output routine instead if relevant
    if data_structure.ife_variables.ife != 0:
        models.ife.output()
        return

    # Costs model
    # Cost switch values
    # No.  |  model
    # ---- | ------
    # 0    |  1990 costs model
    # 1    |  2015 Kovari model
    # 2    |  Custom model
    models.costs.output()

    # Availability model
    models.availability.output()

    # Physics model
    models.physics.output()

    # Detailed physics, currently only done at final point as values are not used
    # by any other functions
    models.physics_detailed.output()

    # TODO what is this? Not in caller.py?
    models.current_drive.output()

    # Pulsed reactor model
    models.pulse.output()

    models.divertor.output()

    # Machine Build Model
    models.build.output()

    # Cryostat build
    models.cryostat.output()

    # Toroidal field coil copper model
    if data_structure.tfcoil_variables.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
        models.copper_tf_coil.output()

    # Toroidal field coil superconductor model
    if data_structure.tfcoil_variables.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
        tf_turn_type = SuperconductingTFTurnType(
            data_structure.superconducting_tf_coil_variables.i_tf_turn_type
        )
        if tf_turn_type == SuperconductingTFTurnType.CABLE_IN_CONDUIT:
            models.cicc_sctfcoil.output()
        elif tf_turn_type == SuperconductingTFTurnType.CROSS_CONDUCTOR:
            models.croco_sctfcoil.output()
        else:
            raise ValueError(
                "Unsupported superconducting TF turn type: "
                f"{data_structure.superconducting_tf_coil_variables.i_tf_turn_type}"
            )

    # Toroidal field coil aluminium model
    if (
        data_structure.tfcoil_variables.i_tf_sup
        == TFConductorModel.HELIUM_COOLED_ALUMINIUM
    ):
        models.aluminium_tf_coil.output()

    # Tight aspect ratio machine model
    if (
        data_structure.physics_variables.itart == 1
        and data_structure.tfcoil_variables.i_tf_sup != TFConductorModel.SUPERCONDUCTING
    ):
        models.tfcoil.output()

    # Poloidal field coil model
    models.pfcoil.output()

    # Structure Model
    models.structure.output()

    # Blanket model
    # Blanket switch values
    # No.  |  model
    # ---- | ------
    # 1    |  CCFE HCPB model
    # 2    |  KIT HCPB model
    # 3    |  CCFE HCPB model with Tritium Breeding Ratio calculation
    # 4    |  KIT HCLL model
    # 5    |  DCLL model

    models.shield.output()
    models.vacuum_vessel.output()

    # First wall geometry
    models.fw.output()

    if data.fwbs.i_blanket_type == BlktModelTypes.CCFE_HCPB:
        # CCFE HCPB model
        models.ccfe_hcpb.output()

    elif data.fwbs.i_blanket_type == BlktModelTypes.DCLL:
        # DCLL model
        models.dcll.output()

    # FISPACT and LOCA model (not used)- removed

    # Power model
    models.power.output()

    # Vacuum model
    models.vacuum.output()

    # Buildings model
    models.buildings.output()

    # Water usage in secondary cooling system
    models.water_use.output()

    # stop capturing warnings so that Outfile does not end up with
    # a lot of non-model logs
    logging_model_handler.stop_capturing()
