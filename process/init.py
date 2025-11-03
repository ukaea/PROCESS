import datetime
import getpass
import socket
import subprocess
from pathlib import Path
from warnings import warn

import process
import process.fortran as fortran
import process.iteration_variables as iteration_variables
import process.process_output as process_output
from process import data_structure
from process.blanket_library import init_blanket_library
from process.build import init_build_variables
from process.buildings import init_buildings_variables
from process.constraints import ConstraintManager, init_constraint_variables
from process.current_drive import init_current_drive_variables
from process.data_structure.build_python_variables import (
    init_build_variables as init_python_build_variables,
)
from process.data_structure.cost_2015_variables import init_cost_2015_variables
from process.data_structure.cost_variables import init_cost_variables
from process.data_structure.cs_fatigue_variables import init_cs_fatigue_variables
from process.data_structure.divertor_variables import init_divertor_variables
from process.data_structure.neoclassics_variables import init_neoclassics_variables
from process.data_structure.power_variables import init_power_variables
from process.data_structure.primary_pumping_variables import (
    init_primary_pumping_variables,
)
from process.data_structure.pulse_variables import init_pulse_variables
from process.data_structure.rebco_variables import init_rebco_variables
from process.data_structure.reinke_variables import init_reinke_variables
from process.data_structure.structure_variables import init_structure_variables
from process.data_structure.times_variables import init_times_variables
from process.data_structure.vacuum_variables import init_vacuum_variables
from process.data_structure.water_usage_variables import init_watuse_variables
from process.dcll import init_dcll_module
from process.exceptions import ProcessValidationError
from process.fw import init_fwbs_variables
from process.hcpb import init_ccfe_hcpb_module
from process.ife import init_ife_variables
from process.impurity_radiation import init_impurity_radiation_module
from process.input import parse_input_file
from process.pfcoil import init_pfcoil_module, init_pfcoil_variables
from process.physics import (
    init_physics_module,
    init_physics_variables,
)
from process.power import init_heat_transport_variables, init_pf_power_variables
from process.scan import init_scan_module
from process.stellarator.initialization import init_stellarator_module, init_stellarator_variables
from process.stellarator.initialization import (
    st_init,
)
from process.superconducting_tf_coil import init_sctfcoil_module
from process.tf_coil import init_tfcoil_variables
from process.utilities.f2py_string_patch import f2py_compatible_to_string


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
    iteration_variables.initialise_iteration_variables()

    # Initialise the Fortran file specifiers
    # (creating and opening the files in the process)
    fortran.init_module.open_files()

    # Input any desired new initial values
    inputs = parse_input_file()

    # Set active constraints
    set_active_constraints()

    # set the device type (icase)
    set_device_type()

    # Initialise the Stellarator
    st_init()

    # Check input data for errors/ambiguities
    check_process(inputs)

    run_summary()


def get_git_summary() -> tuple[str, str]:
    try:
        directory = Path(process.__file__).parent

        git_branch = (
            subprocess.run(
                "git rev-parse --abbrev-ref HEAD",
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

        git_tag = (
            subprocess.run(
                "git describe --tags",
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

        return git_branch, git_tag
    except (subprocess.CalledProcessError, AttributeError):
        return "", ""


def run_summary():
    """Write a summary of the PROCESS run to the output file and MFile"""
    # Outfile and terminal #
    for outfile in [fortran.constants.nout, fortran.constants.iotty]:
        # PROCESS code header
        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)
        process_output.ocentr(outfile, "PROCESS", 110)
        process_output.ocentr(outfile, "Power Reactor Optimisation Code", 110)
        process_output.ostars(outfile, 110)
        process_output.oblnkl(outfile)

        # Run execution details
        version = process.__version__
        process_output.ocmmnt(outfile, f"Version : {version}")

        git_branch, git_tag = get_git_summary()

        process_output.ocmmnt(outfile, f"Git Tag : {git_tag}")
        process_output.ocmmnt(outfile, f"Git Branch : {git_branch}")

        date_string = datetime.datetime.now(datetime.timezone.utc).strftime(
            "%d/%m/%Y %Z"
        )
        time_string = datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M")

        process_output.ocmmnt(outfile, f"Date : {date_string}")
        process_output.ocmmnt(outfile, f"Time : {time_string}")

        user = getpass.getuser()
        machine = socket.gethostname()

        process_output.ocmmnt(outfile, f"User : {user}")
        process_output.ocmmnt(outfile, f"Computer : {machine}")
        process_output.ocmmnt(outfile, f"Directory : {Path.cwd()}")

        fileprefix = f2py_compatible_to_string(fortran.global_variables.fileprefix)
        process_output.ocmmnt(
            outfile,
            f"Input : {fileprefix}",
        )
        runtitle = f2py_compatible_to_string(fortran.global_variables.runtitle)
        process_output.ocmmnt(
            outfile,
            f"Run title : {runtitle}",
        )

        process_output.ocmmnt(
            outfile,
            f"Run type : Reactor concept design: {f2py_compatible_to_string(fortran.global_variables.icase)}, (c) UK Atomic Energy Authority",
        )

        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)
        process_output.oblnkl(outfile)

        process_output.ocmmnt(
            outfile, f"Equality constraints : {fortran.numerics.neqns.item()}"
        )
        process_output.ocmmnt(
            outfile, f"Inequality constraints : {fortran.numerics.nineqns.item()}"
        )
        process_output.ocmmnt(
            outfile,
            f"Total constraints : {fortran.numerics.nineqns.item() + fortran.numerics.neqns.item()}",
        )
        process_output.ocmmnt(
            outfile, f"Iteration variables : {fortran.numerics.nvar.item()}"
        )
        # If optimising, write objective function and convergence parameter
        if fortran.numerics.ioptimz == 1:
            process_output.ocmmnt(
                outfile, f"Max iterations : {fortran.global_variables.maxcal.item()}"
            )

            if fortran.numerics.minmax > 0:
                minmax_string = "  -- minimise "
                minmax_sign = "+"
            else:
                minmax_string = "  -- maximise "
                minmax_sign = "-"

            fom_string = f2py_compatible_to_string(
                fortran.numerics.lablmm[abs(fortran.numerics.minmax) - 1]
            )
            process_output.ocmmnt(
                outfile,
                f"Figure of merit : {minmax_sign}{abs(fortran.numerics.minmax)}{minmax_string}{fom_string}",
            )
            process_output.ocmmnt(
                outfile,
                f"Convergence parameter : {fortran.numerics.epsvmc}",
            )

        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)

    # MFile #
    mfile = fortran.constants.mfile

    process_output.ovarst(mfile, "PROCESS version", "(procver)", f'"{version}"')
    process_output.ovarst(mfile, "Date of run", "(date)", f'"{date_string}"')
    process_output.ovarst(mfile, "Time of run", "(time)", f'"{time_string}"')
    process_output.ovarst(mfile, "User", "(username)", f'"{user}"')
    process_output.ovarst(mfile, "PROCESS run title", "(runtitle)", f'"{runtitle}"')
    process_output.ovarst(mfile, "PROCESS git tag", "(tagno)", f'"{git_tag}"')
    process_output.ovarst(
        mfile, "PROCESS git branch", "(branch_name)", f'"{git_branch}"'
    )
    process_output.ovarst(mfile, "Input filename", "(fileprefix)", f'"{fileprefix}"')

    process_output.ovarin(
        mfile, "Optimisation switch", "(ioptimz)", fortran.numerics.ioptimz
    )
    # If optimising, write figure of merit switch
    if fortran.numerics.ioptimz == 1:
        process_output.ovarin(
            mfile, "Figure of merit switch", "(minmax)", fortran.numerics.minmax
        )


def init_all_module_vars():
    """Initialise all module variables
    This is vital to ensure a 'clean' state of Process before a new run starts,
    otherwise components of the previous run's state can persist into the new
    run. This matters ever since Process is used as a shared library, rather
    than a 'run-once' executable.
    """
    fortran.numerics.init_numerics()
    init_buildings_variables()
    init_cost_variables()
    init_divertor_variables()
    fortran.error_handling.init_error_handling()
    init_fwbs_variables()
    fortran.global_variables.init_global_variables()
    init_ccfe_hcpb_module()
    init_heat_transport_variables()
    init_ife_variables()
    init_impurity_radiation_module()
    init_pfcoil_module()
    init_physics_module()
    init_physics_variables()
    init_scan_module()
    init_sctfcoil_module()
    init_stellarator_module()
    init_stellarator_variables()
    init_tfcoil_variables()
    init_times_variables()
    fortran.constants.init_constants()
    init_current_drive_variables()
    init_primary_pumping_variables()
    init_pfcoil_variables()
    init_structure_variables()
    init_vacuum_variables()
    init_pf_power_variables()
    init_build_variables()
    init_constraint_variables()
    init_pulse_variables()
    init_rebco_variables()
    init_reinke_variables()
    init_watuse_variables()
    init_cs_fatigue_variables()
    init_blanket_library()
    init_dcll_module()
    init_cost_2015_variables()
    init_power_variables()
    init_python_build_variables()
    init_neoclassics_variables()

    fortran.init_module.init_fortran_modules()


def check_process(inputs):  # noqa: ARG001
    """Routine to reset specific variables if certain options are
    being used
    author: P J Knight, CCFE, Culham Science Centre
    None
    This routine performs a sanity check of the input variables
    and ensures other dependent variables are given suitable values.
    """

    # Check that there are sufficient iteration variables
    if fortran.numerics.nvar < fortran.numerics.neqns:
        raise ProcessValidationError(
            "Insufficient iteration variables to solve the problem! NVAR < NEQNS",
            nvar=fortran.numerics.nvar,
            neqns=fortran.numerics.neqns,
        )

    # Check that sufficient elements of ixc and icc have been specified
    if (fortran.numerics.ixc[: fortran.numerics.nvar] == 0).any():
        raise ProcessValidationError(
            "The number of iteration variables specified is smaller than the number stated in ixc",
            nvar=fortran.numerics.nvar,
        )

    # Check that dr_tf_wp_with_insulation (ixc = 140) and dr_tf_inboard (ixc = 13) are not being used simultaneously as iteration variables
    if (fortran.numerics.ixc[: fortran.numerics.nvar] == 13).any() and (
        fortran.numerics.ixc[: fortran.numerics.nvar] == 140
    ).any():
        raise ProcessValidationError(
            "Iteration variables 13 and 140 cannot be used simultaneously",
        )

    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 0
    ).any():
        raise ProcessValidationError(
            "The number of constraints specified is smaller than the number stated in neqns+nineqns",
            neqns=fortran.numerics.neqns,
            nineqns=fortran.numerics.nineqns,
        )

    # Deprecate constraints
    for depcrecated_constraint in [3, 4, 10, 74, 42]:
        if (
            fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns]
            == depcrecated_constraint
        ).any():
            raise ProcessValidationError(
                "Constraint equation is no longer available", icc=depcrecated_constraint
            )

    # MDK Report error if constraint 63 is used with old vacuum model
    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 63
    ).any() and data_structure.vacuum_variables.vacuum_model != "simple":
        raise ProcessValidationError(
            "Constraint 63 is requested without the correct vacuum model (simple)"
        )

    #  Fuel ion fractions must add up to 1.0
    if (
        abs(
            1.0
            - fortran.physics_variables.f_deuterium
            - fortran.physics_variables.f_tritium
            - fortran.physics_variables.f_helium3
        )
        > 1e-6
    ):
        raise ProcessValidationError(
            "Fuel ion fractions do not sum to 1.0",
            f_deuterium=fortran.physics_variables.f_deuterium,
            f_tritium=fortran.physics_variables.f_tritium,
            f_helium3=fortran.physics_variables.f_helium3,
        )

    if fortran.physics_variables.f_tritium < 1.0e-3:  # tritium fraction is negligible
        fortran.buildings_variables.triv = 0.0
        fortran.heat_transport_variables.p_tritium_plant_electric_mw = 0.0

    if fortran.impurity_radiation_module.fimp[1] != 0.1:
        raise ProcessValidationError(
            "The thermal alpha/electron density ratio should be controlled using f_nd_alpha_electron (itv 109) and not fimp(2)."
            "fimp(2) should be removed from the input file, or set to the default value 0.1D0."
        )

    # Impurity fractions
    for imp in range(fortran.impurity_radiation_module.n_impurities):
        fortran.impurity_radiation_module.impurity_arr_frac[imp] = (
            fortran.impurity_radiation_module.fimp[imp]
        )

    # Stop the run if oacdcp is used as an optimisation variable
    # As the current density is now calculated from bt without constraint 10

    if (fortran.numerics.ixc[: fortran.numerics.nvar] == 12).any():
        raise ProcessValidationError(
            "The 1/R toroidal B field dependency constraint is being depreciated"
        )

    # Plasma profile consistency checks
    if fortran.ife_variables.ife != 1 and fortran.physics_variables.ipedestal == 1:
        # Temperature checks
        if fortran.physics_variables.teped < fortran.physics_variables.tesep:
            raise ProcessValidationError(
                "Pedestal temperature is lower than separatrix temperature",
                teped=fortran.physics_variables.teped,
                tesep=fortran.physics_variables.tesep,
            )

        if (abs(fortran.physics_variables.rhopedt - 1.0) <= 1e-7) and (
            (fortran.physics_variables.teped - fortran.physics_variables.tesep) >= 1e-7
        ):
            warn(
                f"Temperature pedestal is at plasma edge, but teped "
                f"({fortran.physics_variables.teped}) differs from tesep "
                f"({fortran.physics_variables.tesep})",
                stacklevel=2,
            )

        # Core temperature should always be calculated (later) as being
        # higher than the pedestal temperature, if and only if the
        # volume-averaged temperature never drops below the pedestal
        # temperature. Prevent this by adjusting te, and its lower bound
        # (which will only have an effect if this is an optimisation run)
        if fortran.physics_variables.te <= fortran.physics_variables.teped:
            warn(
                f"Volume-averaged temperature ({fortran.physics_variables.te}) has been "
                f"forced to exceed input pedestal height ({fortran.physics_variables.teped}). "
                "Changing to te = teped*1.001",
                stacklevel=2,
            )
            fortran.physics_variables.te = fortran.physics_variables.teped * 1.001

        if (
            fortran.numerics.ioptimz >= 0
            and (fortran.numerics.ixc[: fortran.numerics.nvar] == 4).any()
            and fortran.numerics.boundl[3] < fortran.physics_variables.teped * 1.001
        ):
            warn(
                "Lower limit of volume averaged electron temperature (te) has been raised to ensure te > teped",
                stacklevel=2,
            )
            fortran.numerics.boundl[3] = fortran.physics_variables.teped * 1.001
            fortran.numerics.boundu[3] = max(
                fortran.numerics.boundu[3], fortran.numerics.boundl[3]
            )

        # Density checks
        # Case where pedestal density is set manually
        if (
            fortran.physics_variables.fgwped < 0
            or not (fortran.numerics.ixc[: fortran.numerics.nvar] == 145).any()
        ):
            # Issue #589 Pedestal density is set manually using neped but it is less than nesep.
            if fortran.physics_variables.neped < fortran.physics_variables.nesep:
                raise ProcessValidationError(
                    "Density pedestal is lower than separatrix density",
                    neped=fortran.physics_variables.neped,
                    nesep=fortran.physics_variables.nesep,
                )

            # Issue #589 Pedestal density is set manually using neped,
            # but pedestal width = 0.
            if (
                abs(fortran.physics_variables.rhopedn - 1.0) <= 1e-7
                and (fortran.physics_variables.neped - fortran.physics_variables.nesep)
                >= 1e-7
            ):
                warn(
                    "Density pedestal is at plasma edge "
                    f"({fortran.physics_variables.rhopedn = }), but neped "
                    f"({fortran.physics_variables.neped}) differs from "
                    f"nesep ({fortran.physics_variables.nesep})",
                    stacklevel=2,
                )

        # Issue #862 : Variable ne0/neped ratio without constraint eq 81 (ne0>neped)
        #  -> Potential hollowed density profile
        if (
            fortran.numerics.ioptimz >= 0
            and not (
                fortran.numerics.icc[
                    : fortran.numerics.neqns + fortran.numerics.nineqns
                ]
                == 81
            ).any()
        ):
            if (fortran.numerics.ixc[: fortran.numerics.nvar] == 145).any():
                warn(
                    "neped set with fgwped without constraint eq 81 (neped<ne0)",
                    stacklevel=2,
                )
            if (fortran.numerics.ixc[: fortran.numerics.nvar] == 6).any():
                warn(
                    "dene used as iteration variable without constraint 81 (neped<ne0)",
                    stacklevel=2,
                )

    # Cannot use Psep/R and PsepB/qAR limits at the same time
    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 68
    ).any() and (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 56
    ).any():
        raise ProcessValidationError(
            "Cannot use Psep/R and PsepB/qAR constraint equations at the same time"
        )

    # if lower bound of fgwped < fgwsep
    if (
        fortran.numerics.ixc[: fortran.numerics.nvar] == 145
    ).any() and fortran.numerics.boundl[144] < fortran.physics_variables.fgwsep:
        raise ProcessValidationError(
            "Set lower bound of iteration variable 145, fgwped, to be greater than fgwsep",
            boundl_145=fortran.numerics.boundl[144],
            fgwsep=fortran.physics_variables.fgwsep,
        )

    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 78
    ).any():
        # If Reinke criterion is used tesep is calculated and cannot be an
        # iteration variable
        if (fortran.numerics.ixc[: fortran.numerics.nvar] == 119).any():
            raise ProcessValidationError(
                "REINKE IMPURITY MODEL: tesep is calculated and cannot be an "
                "iteration variable for the Reinke model"
            )

        # If Reinke criterion is used need to enforce LH-threshold
        # using Martin scaling for consistency
        if (fortran.physics_variables.i_l_h_threshold != 6) or (
            not (
                fortran.numerics.icc[
                    : fortran.numerics.neqns + fortran.numerics.nineqns
                ]
                == 15
            ).any()
            and fortran.physics_variables.ipedestal
        ):
            warn(
                "REINKE IMPURITY MODEL: The Martin LH threshold scale is not being used and is recommned for the Reinke model",
                stacklevel=2,
            )

    if fortran.physics_variables.i_single_null == 0:
        fortran.physics_variables.n_divertors = 2
        fortran.build_variables.dz_fw_plasma_gap = (
            fortran.build_variables.dz_xpoint_divertor
        )
        fortran.build_variables.dz_shld_upper = fortran.build_variables.dz_shld_lower
        fortran.build_variables.dz_vv_upper = fortran.build_variables.dz_vv_lower
        warn("Double-null: Upper vertical build forced to match lower", stacklevel=2)
    else:  # i_single_null == 1
        fortran.physics_variables.n_divertors = 1

    #  Tight aspect ratio options (ST)
    if fortran.physics_variables.itart == 1:
        fortran.global_variables.icase = "Tight aspect ratio tokamak model"

        # Disabled Forcing that no inboard breeding blanket is used
        # Disabled i_blkt_inboard = 0

        # Check if the choice of plasma current is addapted for ST
        # 2 : Peng Ip scaling (See STAR code documentation)
        # 9 : Fiesta Ip scaling
        if (
            fortran.physics_variables.i_plasma_current != 2
            and fortran.physics_variables.i_plasma_current != 9
        ):
            warn(
                "Usual current scaling for TARTs (i_plasma_current=2 or 9) is not being used",
                stacklevel=2,
            )

        # If using Peng and Strickler (1986) model (itartpf == 0)
        # Overwrite the location of the TF coils
        # 2 : PF coil on top of TF coil
        # 3 : PF coil outside of TF coil
        if fortran.physics_variables.itartpf == 0:
            fortran.pfcoil_variables.i_pf_location[0] = 2
            fortran.pfcoil_variables.i_pf_location[1] = 3
            fortran.pfcoil_variables.i_pf_location[2] = 3

        # Water cooled copper magnets initalisation / checks
        if fortran.tfcoil_variables.i_tf_sup == 0:
            # Check if the initial centrepost coolant loop adapted to the magnet technology
            # Ice cannot flow so tcoolin > 273.15 K
            if fortran.tfcoil_variables.tcoolin < 273.15:
                raise ProcessValidationError(
                    "Coolant temperature (tcoolin) cannot be < 0 C (273.15 K) for water cooled copper magents"
                )

            # Temperature of the TF legs cannot be cooled down
            if (
                fortran.tfcoil_variables.temp_tf_legs_outboard > 0
                and fortran.tfcoil_variables.temp_tf_legs_outboard < 273.15
            ):
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) cannot be < 0 C (273.15 K) for water cooled magents"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                fortran.numerics.ixc[: fortran.numerics.nvar] == 20
            ).any() and fortran.numerics.boundu[19] < 273.15:
                raise ProcessValidationError(
                    "Too low CP conductor temperature (temp_cp_average). Lower limit for copper > 273.15 K"
                )

        # Call a lvl 3 error if superconductor magnets are used
        elif fortran.tfcoil_variables.i_tf_sup == 1:
            warn(
                "Joints res not cal. for SC (itart = 1) TF (fortran.tfcoil_variables.i_tf_sup = 1)",
                stacklevel=2,
            )

        # Aluminium magnets initalisation / checks
        # Initialize the CP conductor temperature to cryogenic temperature for cryo-al magnets (20 K)
        elif fortran.tfcoil_variables.i_tf_sup == 2:
            # Call a lvl 3 error if the inlet coolant temperature is too large
            # Motivation : ill-defined aluminium resistivity fit for T > 40-50 K
            if fortran.tfcoil_variables.tcoolin > 40.0:
                raise ProcessValidationError(
                    "Coolant temperature (tcoolin) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if the leg average temperature is low enough for the resisitivity fit
            if fortran.tfcoil_variables.temp_tf_legs_outboard > 50.0:
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                fortran.numerics.ixc[: fortran.numerics.nvar] == 20
            ).any() and fortran.numerics.boundu[19] > 50.0:
                raise ProcessValidationError(
                    "Too large CP conductor temperature (temp_cp_average). Upper limit for cryo-al < 50 K"
                )

            # Otherwise intitialise the average conductor temperature at
            fortran.tfcoil_variables.temp_cp_average = fortran.tfcoil_variables.tcoolin

        # Check if the boostrap current selection is addapted to ST
        if fortran.physics_variables.i_bootstrap_current == 1:
            raise ProcessValidationError(
                "Invalid boostrap current law for ST, do not use i_bootstrap_current = 1"
            )

        # Check if a single null divertor is used in double null machine
        if fortran.physics_variables.i_single_null == 0 and (
            fortran.physics_variables.f_p_div_lower == 1.0
            or fortran.physics_variables.f_p_div_lower == 0.0
        ):
            warn("Operating with a single null in a double null machine", stacklevel=2)

        # Set the TF coil shape to picture frame (if default value)
        if fortran.tfcoil_variables.i_tf_shape == 0:
            fortran.tfcoil_variables.i_tf_shape = 2

        # Warning stating that the CP fast neutron fluence calculation
        # is not addapted for cryoaluminium calculations yet
        if (
            fortran.tfcoil_variables.i_tf_sup == 2
            and (
                fortran.numerics.icc[
                    : fortran.numerics.neqns + fortran.numerics.nineqns
                ]
                == 85
            ).any()
            and fortran.physics_variables.itart == 1
        ):
            raise ProcessValidationError(
                "Al TF coil fluence not calculated properly for Al CP, do not use constraint 85"
            )

        # Setting the CP joints default options :
        #  0 : No joints for superconducting magents (fortran.tfcoil_variables.i_tf_sup = 1)
        #  1 : Sliding joints for resistive magnets (fortran.tfcoil_variables.i_tf_sup = 0, 2)
        if fortran.tfcoil_variables.i_cp_joints == -1:
            if fortran.tfcoil_variables.i_tf_sup == 1:
                fortran.tfcoil_variables.i_cp_joints = 0
            else:
                fortran.tfcoil_variables.i_cp_joints = 1

        # Checking the CP TF top radius
        if (
            abs(fortran.build_variables.r_cp_top) > 0
            or (fortran.numerics.ixc[: fortran.numerics.nvar] == 174).any()
        ) and fortran.build_variables.i_r_cp_top != 1:
            raise ProcessValidationError(
                "To set the TF CP top value, you must use i_r_cp_top = 1"
            )

    # Conventionnal aspect ratios specific
    else:
        if (
            fortran.physics_variables.i_plasma_current == 2
            or fortran.physics_variables.i_plasma_current == 9
        ):
            raise ProcessValidationError(
                "i_plasma_current=2,9 is not a valid option for a non-TART device"
            )

        # Set the TF coil shape to PROCESS D-shape (if default value)
        if fortran.tfcoil_variables.i_tf_shape == 0:
            fortran.tfcoil_variables.i_tf_shape = 1

        # Check PF coil configurations
        j = 0
        k = 0
        for i in range(fortran.pfcoil_variables.n_pf_coil_groups):
            if (
                fortran.pfcoil_variables.i_pf_location[i] != 2
                and fortran.pfcoil_variables.n_pf_coils_in_group[i] != 2
            ):
                raise ProcessValidationError(
                    "n_pf_coils_in_group(i) .ne. 2 is not a valid option except for (i_pf_location = 2)"
                )

            if fortran.pfcoil_variables.i_pf_location[i] == 2:
                j = j + 1
                k = k + fortran.pfcoil_variables.n_pf_coils_in_group[i]

        if k == 1:
            raise ProcessValidationError(
                "Only 1 divertor coil (i_pf_location = 2) is not a valid configuration"
            )
        if k > 2:
            raise ProcessValidationError(
                "More than 2 divertor coils (i_pf_location = 2) is not a valid configuration"
            )
        if fortran.physics_variables.i_single_null == 1 and j < 2:
            raise ProcessValidationError(
                "If i_single_null=1, use 2 individual divertor coils (i_pf_location = 2, 2; n_pf_coils_in_group = 1, 1)"
            )

        # Constraint 10 is dedicated to ST designs with demountable joints
        if (
            fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns]
            == 10
        ).any():
            raise ProcessValidationError(
                "Constraint equation 10 (CP lifetime) to used with ST desing (itart=1)"
            )

    #  Pulsed power plant model
    if data_structure.pulse_variables.i_pulsed_plant == 1:
        fortran.global_variables.icase = "Pulsed tokamak model"
    else:
        fortran.buildings_variables.esbldgm3 = 0.0

    # TF coil
    # -------
    # TF stress model not defined of r_tf_inboard = 0
    # Unless i_tf_stress_model == 2
    # -> If dr_bore + dr_cs_tf_gap + dr_cs = 0 and fixed and stress constraint is used
    #    Generate a lvl 3 error proposing not to use any stress constraints
    if (
        (
            not (
                (fortran.numerics.ixc[: fortran.numerics.nvar] == 16).any()
                or (fortran.numerics.ixc[: fortran.numerics.nvar] == 29).any()
                or (fortran.numerics.ixc[: fortran.numerics.nvar] == 42).any()
            )
        )  # No dr_bore,dr_cs_tf_gap, dr_cs iteration
        and (
            abs(
                fortran.build_variables.dr_bore
                + fortran.build_variables.dr_cs_tf_gap
                + fortran.build_variables.dr_cs
                + fortran.build_variables.dr_cs_precomp
            )
            <= 0
        )  # dr_bore + dr_cs_tf_gap + dr_cs = 0
        and (
            (
                fortran.numerics.icc[
                    : fortran.numerics.neqns + fortran.numerics.nineqns
                ]
                == 31
            ).any()
            or (
                fortran.numerics.icc[
                    : fortran.numerics.neqns + fortran.numerics.nineqns
                ]
                == 32
            ).any()
        )  # Stress constraints (31 or 32) is used
        and (
            fortran.tfcoil_variables.i_tf_stress_model != 2
        )  # TF stress model can't handle no dr_bore
    ):
        raise ProcessValidationError(
            "Invalid stress model if dr_bore + dr_cs_tf_gap + dr_cs = 0. Don't use constraint 31"
        )

    # Make sure that plane stress model is not used for resistive magnets
    if (
        fortran.tfcoil_variables.i_tf_stress_model == 1
        and fortran.tfcoil_variables.i_tf_sup != 1
    ):
        raise ProcessValidationError(
            "Use generalized plane strain for resistive magnets (i_tf_stress_model = 0 or 2)"
        )

    # bucking cylinder default option setting
    # - bucking (casing) for SC i_tf_bucking ( i_tf_bucking = 1 )
    # - No bucking for copper magnets ( i_tf_bucking = 0 )
    # - Bucking for aluminium magnets ( i_tf_bucking = 1 )
    if fortran.tfcoil_variables.i_tf_bucking == -1:
        if fortran.tfcoil_variables.i_tf_sup == 0:
            fortran.tfcoil_variables.i_tf_bucking = 0
        else:
            fortran.tfcoil_variables.i_tf_bucking = 1

    # Ensure that the TF isnt placed against the
    # CS which is now outside it
    if (
        fortran.tfcoil_variables.i_tf_bucking >= 2
        and fortran.build_variables.i_tf_inside_cs == 1
    ):
        raise ProcessValidationError(
            "Cannot have i_tf_bucking >= 2 when i_tf_inside_cs = 1"
        )

    # Ensure that no pre-compression structure
    # is used for bucked and wedged design
    if (
        fortran.tfcoil_variables.i_tf_bucking >= 2
        and fortran.build_variables.i_cs_precomp == 1
    ):
        raise ProcessValidationError(
            "No CS precompression structure for bucked and wedged, use i_cs_precomp = 0"
        )

    # Number of stress calculation layers
    # +1 to add in the inboard TF coil case on the plasma side, per Issue #1509
    fortran.tfcoil_variables.n_tf_stress_layers = (
        fortran.tfcoil_variables.i_tf_bucking
        + fortran.tfcoil_variables.n_tf_graded_layers
        + 1
    )

    # If TFC sidewall has not been set by user
    if fortran.tfcoil_variables.dx_tf_side_case_min < 0.1e-10:
        fortran.tfcoil_variables.tfc_sidewall_is_fraction = True

    # If inboard TF coil case plasma side thickness has not been set by user
    if fortran.tfcoil_variables.dr_tf_plasma_case < 0.1e-10:
        fortran.tfcoil_variables.i_f_dr_tf_plasma_case = True

    # Setting the default cryo-plants efficiencies
    if abs(fortran.tfcoil_variables.eff_tf_cryo + 1) < 1e-6:
        # The ITER cyoplant efficiency is used for SC
        if fortran.tfcoil_variables.i_tf_sup == 1:
            fortran.tfcoil_variables.eff_tf_cryo = 0.13

        # Strawbrige plot extrapolation is used for Cryo-Al
        elif fortran.tfcoil_variables.i_tf_sup == 2:
            fortran.tfcoil_variables.eff_tf_cryo = 0.40

    # Cryo-plane efficiency must be in [0-1.0]
    elif (
        fortran.tfcoil_variables.eff_tf_cryo > 1.0
        or fortran.tfcoil_variables.eff_tf_cryo < 0.0
    ):
        raise ProcessValidationError(
            "TF cryo-plant efficiency `eff_tf_cryo` must be within [0-1]"
        )

    # Integer turns option not yet available for REBCO taped turns

    if (
        fortran.tfcoil_variables.i_tf_sc_mat == 6
        and fortran.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Integer turns (i_tf_turns_integer = 1) not supported for REBCO (i_tf_sc_mat = 6)"
        )

    # Setting up insulation layer young modulae default values [Pa]

    if fortran.tfcoil_variables.eyoung_ins <= 1.0e8:
        # Copper magnets, no insulation material defined
        # But use the ITER design by default
        if (
            fortran.tfcoil_variables.i_tf_sup == 0
            or fortran.tfcoil_variables.i_tf_sup == 1
        ):
            # SC magnets
            # Value from DDD11-2 v2 2 (2009)
            fortran.tfcoil_variables.eyoung_ins = 20.0e9

        # Cryo-aluminum magnets (Kapton polymer)
        elif fortran.tfcoil_variables.i_tf_sup == 2:
            fortran.tfcoil_variables.eyoung_ins = 2.5e9

    # Setting the default WP geometry

    if fortran.tfcoil_variables.i_tf_wp_geom == -1:
        if fortran.tfcoil_variables.i_tf_turns_integer == 0:
            fortran.tfcoil_variables.i_tf_wp_geom = 1
        if fortran.tfcoil_variables.i_tf_turns_integer == 1:
            fortran.tfcoil_variables.i_tf_wp_geom = 0

    # Setting the TF coil conductor elastic properties

    if fortran.tfcoil_variables.i_tf_cond_eyoung_axial == 0:
        # Conductor stiffness is not considered
        fortran.tfcoil_variables.eyoung_cond_axial = 0
        fortran.tfcoil_variables.eyoung_cond_trans = 0
    elif fortran.tfcoil_variables.i_tf_cond_eyoung_axial == 2:
        # Select sensible defaults from the literature
        if fortran.tfcoil_variables.i_tf_sc_mat in [1, 4, 5]:
            # Nb3Sn: Nyilas, A et. al, Superconductor Science and Technology 16, no. 9 (2003): 1036-42. https://doi.org/10.1088/0953-2048/16/9/313.
            fortran.tfcoil_variables.eyoung_cond_axial = 32e9
        elif fortran.tfcoil_variables.i_tf_sc_mat == 2:
            # Bi-2212: Brown, M. et al, IOP Conference Series: Materials Science and Engineering 279 (2017): 012022. https://doi.org/10.1088/1757-899X/279/1/012022.
            fortran.tfcoil_variables.eyoung_cond_axial = 80e9
        elif fortran.tfcoil_variables.i_tf_sc_mat in [3, 7]:
            # NbTi: Vedrine, P. et. al, IEEE Transactions on Applied Superconductivity 9, no. 2 (1999): 236-39. https://doi.org/10.1109/77.783280.
            fortran.tfcoil_variables.eyoung_cond_axial = 6.8e9
        elif fortran.tfcoil_variables.i_tf_sc_mat in [6, 8, 9]:
            # REBCO: Fujishiro, H. et. al, Physica C: Superconductivity, 426-431 (2005): 699-704. https://doi.org/10.1016/j.physc.2005.01.045.
            fortran.tfcoil_variables.eyoung_cond_axial = 145e9

        if fortran.tfcoil_variables.i_tf_cond_eyoung_trans == 0:
            # Transverse stiffness is not considered
            fortran.tfcoil_variables.eyoung_cond_trans = 0
        else:
            # Transverse stiffness is significant
            fortran.tfcoil_variables.eyoung_cond_trans = (
                fortran.tfcoil_variables.eyoung_cond_axial
            )

    # Check if the WP/conductor radial thickness (dr_tf_wp_with_insulation) is large enough
    # To contains the insulation, cooling and the support structure
    # Rem : Only verified if the WP thickness is used
    if (fortran.numerics.ixc[: fortran.numerics.nvar] == 140).any():
        # Minimal WP thickness
        if fortran.tfcoil_variables.i_tf_sup == 1:
            dr_tf_wp_min = 2.0 * (
                fortran.tfcoil_variables.dx_tf_wp_insulation
                + fortran.tfcoil_variables.dx_tf_wp_insertion_gap
                + fortran.tfcoil_variables.dx_tf_turn_insulation
                + fortran.tfcoil_variables.dia_tf_turn_coolant_channel
            )

            # Steel conduit thickness (can be an iteration variable)
            if (fortran.numerics.ixc[: fortran.numerics.nvar] == 58).any():
                dr_tf_wp_min = dr_tf_wp_min + 2.0 * fortran.numerics.boundl[57]
            else:
                dr_tf_wp_min = (
                    dr_tf_wp_min + 2.0 * fortran.tfcoil_variables.dx_tf_turn_steel
                )

        # Minimal conductor layer thickness
        elif (
            fortran.tfcoil_variables.i_tf_sup == 0
            or fortran.tfcoil_variables.i_tf_sup == 2
        ):
            dr_tf_wp_min = (
                2.0
                * (
                    fortran.tfcoil_variables.dx_tf_turn_insulation
                    + fortran.tfcoil_variables.dx_tf_wp_insulation
                )
                + 4.0 * fortran.tfcoil_variables.rcool
            )

        if fortran.numerics.boundl[139] < dr_tf_wp_min:
            raise ProcessValidationError(
                "The TF coil WP thickness (dr_tf_wp_with_insulation) must be at least",
                dr_tf_wp_min=dr_tf_wp_min,
            )

    # Setting t_turn_tf_is_input to true if t_turn_tf is an input
    fortran.tfcoil_variables.t_turn_tf_is_input = (
        abs(fortran.tfcoil_variables.t_turn_tf) > 0
    )

    # Impossible to set the turn size of integer turn option
    if (
        fortran.tfcoil_variables.t_turn_tf_is_input
        and fortran.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    if (
        fortran.tfcoil_variables.i_tf_wp_geom != 0
        and fortran.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Can only have i_tf_turns_integer = 1 with i_tf_wp_geom = 0"
        )

    if (
        fortran.physics_variables.i_bootstrap_current == 5
        and fortran.physics_variables.i_diamagnetic_current != 0
    ):
        raise ProcessValidationError(
            "i_diamagnetic_current = 0 should be used with the Sakai plasma current scaling"
        )

    # Setting t_cable_tf_is_input to true if t_cable_tf is an input
    fortran.tfcoil_variables.t_cable_tf_is_input = (
        abs(fortran.tfcoil_variables.t_cable_tf) > 0
    )

    # Impossible to set the cable size of integer turn option
    if (
        fortran.tfcoil_variables.t_cable_tf_is_input
        and fortran.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    # Impossible to set both the TF coil turn and the cable dimension
    if (
        fortran.tfcoil_variables.t_turn_tf_is_input
        and fortran.tfcoil_variables.t_cable_tf_is_input
    ):
        raise ProcessValidationError(
            "Impossible to set the TF coil turn and cable size simultaneously"
        )

    # Checking the SC temperature for LTS
    if (
        fortran.tfcoil_variables.i_tf_sc_mat in [1, 3, 4, 5]
        and fortran.tfcoil_variables.tftmp > 10.0
    ):
        raise ProcessValidationError(
            "The LTS conductor temperature (tftmp) has to be lower than 10"
        )

    # PF coil resistivity is zero if superconducting
    if fortran.pfcoil_variables.i_pf_conductor == 0:
        fortran.pfcoil_variables.rho_pf_coil = 0.0

    # If there is no NBI, then hot beam density should be zero
    if fortran.current_drive_variables.i_hcd_calculations == 1:
        if (
            fortran.current_drive_variables.i_hcd_primary != 5
            and fortran.current_drive_variables.i_hcd_primary != 8
        ):
            fortran.physics_variables.f_nd_beam_electron = 0.0
    else:
        fortran.physics_variables.f_nd_beam_electron = 0.0

    # Set inboard blanket thickness to zero if no inboard blanket switch
    # used (Issue #732)
    if fortran.fwbs_variables.i_blkt_inboard == 0:
        fortran.build_variables.dr_blkt_inboard = 0.0

    # Ensure that blanket material fractions allow non-zero space for steel
    # CCFE HCPB Model

    if fortran.stellarator_variables.istell == 0 and (
        fortran.fwbs_variables.i_blanket_type == 1
    ):
        fsum = (
            fortran.fwbs_variables.breeder_multiplier
            + fortran.fwbs_variables.vfcblkt
            + fortran.fwbs_variables.vfpblkt
        )
        if fsum >= 1.0:
            raise ProcessValidationError(
                "Blanket material fractions do not sum to 1.0",
                i_blanket_type=fortran.fwbs_variables.i_blanket_type,
                breeder_multiplier=fortran.fwbs_variables.breeder_multiplier,
                vfcblkt=fortran.fwbs_variables.vfcblkt,
                vfpblkt=fortran.fwbs_variables.vfpblkt,
                fsum=fsum,
            )

    # Check that the temperature margins are not overdetermined
    if fortran.tfcoil_variables.tmargmin > 0.0001:
        # This limit has been input and will be applied to both TFC and CS
        if fortran.tfcoil_variables.tmargmin_tf > 0.0001:
            warn(
                "tmargmin_tf and tmargmin should not both be specified in IN.DAT "
                "tmargmin_tf has been ignored",
                stacklevel=2,
            )
        if fortran.tfcoil_variables.tmargmin_cs > 0.0001:
            warn(
                "tmargmin_cs and tmargmin should not both be specified in IN.DAT "
                "tmargmin_cs has been ignored",
                stacklevel=2,
            )

        fortran.tfcoil_variables.tmargmin_tf = fortran.tfcoil_variables.tmargmin
        fortran.tfcoil_variables.tmargmin_cs = fortran.tfcoil_variables.tmargmin

    if (
        fortran.physics_variables.tauee_in > 1e-10
        and fortran.physics_variables.i_confinement_time != 48
    ):
        # Report error if confinement time is in the input
        # but the scaling to use it is not selected.
        warn("tauee_in is for use with i_confinement_time=48 only", stacklevel=2)

    if (
        fortran.physics_variables.aspect > 1.7
        and fortran.physics_variables.i_confinement_time == 46
    ):
        # NSTX scaling is for A<1.7
        warn("NSTX scaling is for A<1.7", stacklevel=2)

    if (
        fortran.physics_variables.i_plasma_current == 2
        and fortran.physics_variables.i_confinement_time == 42
    ):
        raise ProcessValidationError(
            "Lang 2012 confinement scaling cannot be used for i_plasma_current=2 due to wrong q"
        )

    # Cannot use temperature margin constraint with REBCO TF coils
    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 36
    ).any() and (
        fortran.tfcoil_variables.i_tf_sc_mat == 8
        or fortran.tfcoil_variables.i_tf_sc_mat == 9
    ):
        raise ProcessValidationError(
            "turn off TF temperature margin constraint icc = 36 when using REBCO"
        )

    # Cannot use temperature margin constraint with REBCO CS coils
    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 60
    ).any() and fortran.pfcoil_variables.i_cs_superconductor == 8:
        raise ProcessValidationError(
            "turn off CS temperature margin constraint icc = 60 when using REBCO"
        )

    # Cold end of the cryocooler should be colder than the TF
    if fortran.tfcoil_variables.temp_tf_cryo > fortran.tfcoil_variables.tftmp:
        raise ProcessValidationError("temp_tf_cryo should be lower than tftmp")

    # Cannot use TF coil strain limit if i_str_wp is off:
    if (
        fortran.numerics.icc[: fortran.numerics.neqns + fortran.numerics.nineqns] == 88
    ).any() and fortran.tfcoil_variables.i_str_wp == 0:
        raise ProcessValidationError("Can't use constraint 88 if i_strain_tf == 0")


def set_active_constraints():
    """Set constraints provided in the input file as 'active'"""
    num_constraints = 0
    for i in range(ConstraintManager.num_constraints()):
        if fortran.numerics.icc[i] != 0:
            fortran.numerics.active_constraints[fortran.numerics.icc[i] - 1] = True
            num_constraints += 1

    if fortran.numerics.neqns == 0:
        # The value of neqns has not been set in the input file.  Default = 0.
        fortran.numerics.neqns = num_constraints - fortran.numerics.nineqns
    else:
        fortran.numerics.nineqns = num_constraints - fortran.numerics.neqns


def set_device_type():
    if fortran.ife_variables.ife == 1:
        fortran.global_variables.icase = "Inertial Fusion model"
    elif fortran.stellarator_variables.istell != 0:
        fortran.global_variables.icase = "Stellarator model"
