from process import data_structure
from process.log import logging_model_handler


def write(models, _outfile):
    """Write the results to the main output file (OUT.DAT).

    Write the program results to a file, in a tidy format.

    :param models: physics and engineering model objects
    :type models: process.main.Models
    :param outfile: Fortran output unit identifier
    :type outfile: int
    """
    # ensure we are capturing warnings that occur in the 'output' stage as these are warnings
    # that occur at our solution point. So we clear existing warnings
    logging_model_handler.start_capturing()
    logging_model_handler.clear_logs()

    # Call stellarator output routine instead if relevant
    if data_structure.stellarator_variables.istell != 0:
        models.stellarator.run(output=True)
        return

    #  Call IFE output routine instead if relevant
    if data_structure.ife_variables.ife != 0:
        models.ife.run(output=True)
        return

    # Costs model
    # Cost switch values
    # No.  |  model
    # ---- | ------
    # 0    |  1990 costs model
    # 1    |  2015 Kovari model
    # 2    |  Custom model
    models.costs.run()
    models.costs.output()

    # Availability model
    models.availability.run(output=True)

    # Writing the output from physics.f90 into OUT.DAT + MFILE.DAT
    models.physics.outplas()

    # TODO what is this? Not in caller.f90?
    models.current_drive.output_current_drive()

    # Pulsed reactor model
    models.pulse.run(output=True)
    models.physics.outtim()

    models.divertor.run(output=True)

    # Machine Build Model
    # Radial build
    models.build.calculate_radial_build(output=True)

    # Vertical build
    models.build.calculate_vertical_build(output=True)

    # Cryostat build
    models.cryostat.cryostat_output()

    # Toroidal field coil copper model
    if data_structure.tfcoil_variables.i_tf_sup == 0:
        models.copper_tf_coil.run(output=True)

    # Toroidal field coil superconductor model
    if data_structure.tfcoil_variables.i_tf_sup == 1:
        models.sctfcoil.run(output=True)
        models.sctfcoil.output_tf_superconductor_info()

    # Toroidal field coil aluminium model
    if data_structure.tfcoil_variables.i_tf_sup == 2:
        models.aluminium_tf_coil.run(output=True)

    # Tight aspect ratio machine model
    if (
        data_structure.physics_variables.itart == 1
        and data_structure.tfcoil_variables.i_tf_sup != 1
    ):
        models.tfcoil.iprint = 1
        models.tfcoil.cntrpst()
        models.tfcoil.iprint = 0

    # Poloidal field coil model
    models.pfcoil.output()

    # Structure Model
    models.structure.run(output=True)

    # Poloidal field coil inductance calculation
    models.pfcoil.output_induct()

    # Blanket model
    # Blanket switch values
    # No.  |  model
    # ---- | ------
    # 1    |  CCFE HCPB model
    # 2    |  KIT HCPB model
    # 3    |  CCFE HCPB model with Tritium Breeding Ratio calculation
    # 4    |  KIT HCLL model
    # 5    |  DCLL model

    # First wall geometry
    models.fw.output_fw_geometry()

    # First wall pumping
    models.fw.output_fw_pumping()

    if data_structure.fwbs_variables.i_blanket_type == 1:
        # CCFE HCPB model
        models.ccfe_hcpb.run(output=True)
    # i_blanket_type = 2, KIT HCPB removed
    # i_blanket_type = 3, CCFE HCPB with TBR calculation removed
    # i_blanket_type = 4, KIT HCLL removed
    elif data_structure.fwbs_variables.i_blanket_type == 5:
        # DCLL model
        models.dcll.run(output=True)

    # FISPACT and LOCA model (not used)- removed

    # Toroidal field coil power model
    models.power.tfpwr(output=True)

    # Poloidal field coil power model !
    models.power.pfpwr(output=True)

    # Vacuum model
    models.vacuum.run(output=True)

    # Buildings model
    models.buildings.run(output=True)

    # Plant AC power requirements
    models.power.acpow(output=True)

    # Plant heat transport pt 2 & 3
    models.power.output_cryogenics()
    models.power.output_plant_thermal_powers()
    models.power.output_plant_electric_powers()

    # Water usage in secondary cooling system
    models.water_use.run(output=True)

    # stop capturing warnings so that Outfile does not end up with
    # a lot of non-model logs
    logging_model_handler.stop_capturing()
