import process.fortran as fortran


def init_process():
    """Routine that calls the initialisation routines
    author: P J Knight, CCFE, Culham Science Centre
    None
    This routine calls the main initialisation routines that set
    the default values for the global variables, reads in data from
    the input file, and checks the run parameters for consistency.
    """
    # Initialise error handling
    fortran.error_handling.initialise_error_list()

    # Initialise the program variables
    fortran.init_module.initial()

    # Initialise the Fortran file specifiers
    # (creating and opening the files in the process)
    fortran.init_module.open_files()

    # Input any desired new initial values
    fortran.process_input.input()

    # Initialise the Stellarator
    fortran.stellarator_module.stinit()

    # Check input data for errors/ambiguities
    fortran.init_module.check()

    fortran.main_module.run_summary()


def init_all_module_vars():
    """Initialise all module variables
    This is vital to ensure a 'clean' state of Process before a new run starts,
    otherwise components of the previous run's state can persist into the new
    run. This matters ever since Process is used as a shared library, rather
    than a 'run-once' executable.
    """
    fortran.numerics.init_numerics()
    fortran.process_input.init_input()
    fortran.buildings_variables.init_buildings_variables()
    fortran.cost_variables.init_cost_variables()
    fortran.divertor_variables.init_divertor_variables()
    fortran.error_handling.init_error_handling()
    fortran.fwbs_variables.init_fwbs_variables()
    fortran.global_variables.init_global_variables()
    fortran.ccfe_hcpb_module.init_ccfe_hcpb_module()
    fortran.heat_transport_variables.init_heat_transport_variables()
    fortran.ife_variables.init_ife_variables()
    fortran.impurity_radiation_module.init_impurity_radiation_module()
    fortran.pfcoil_module.init_pfcoil_module()
    fortran.physics_module.init_physics_module()
    fortran.physics_variables.init_physics_variables()
    fortran.scan_module.init_scan_module()
    fortran.sctfcoil_module.init_sctfcoil_module()
    fortran.stellarator_module.init_stellarator_module()
    fortran.stellarator_variables.init_stellarator_variables()
    fortran.tfcoil_variables.init_tfcoil_variables()
    fortran.times_variables.init_times_variables()
    fortran.constants.init_constants()
    fortran.current_drive_variables.init_current_drive_variables()
    fortran.primary_pumping_variables.init_primary_pumping_variables()
    fortran.pfcoil_variables.init_pfcoil_variables()
    fortran.structure_variables.init_structure_variables()
    fortran.vacuum_variables.init_vacuum_variables()
    fortran.pf_power_variables.init_pf_power_variables()
    fortran.build_variables.init_build_variables()
    fortran.constraint_variables.init_constraint_variables()
    fortran.pulse_variables.init_pulse_variables()
    fortran.rebco_variables.init_rebco_variables()
    fortran.reinke_variables.init_reinke_variables()
    fortran.define_iteration_variables.init_define_iteration_variables()
    fortran.reinke_module.init_reinke_module()
    fortran.water_usage_variables.init_watuse_variables()
    fortran.cs_fatigue_variables.init_cs_fatigue_variables()
    fortran.blanket_library.init_blanket_library()
    fortran.dcll_module.init_dcll_module()

    fortran.init_module.init_fortran_modules()
