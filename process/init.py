import datetime
import getpass
import socket
import subprocess
from pathlib import Path
from warnings import warn

import process
import process.iteration_variables as iteration_variables
import process.process_output as process_output
from process import data_structure
from process.constraints import ConstraintManager
from process.core import constants
from process.data_structure.blanket_library import init_blanket_library
from process.data_structure.build_variables import init_build_variables
from process.data_structure.buildings_variables import init_buildings_variables
from process.data_structure.ccfe_hcpb_module import init_ccfe_hcpb_module
from process.data_structure.constraint_variables import init_constraint_variables
from process.data_structure.cost_2015_variables import init_cost_2015_variables
from process.data_structure.cost_variables import init_cost_variables
from process.data_structure.cs_fatigue_variables import init_cs_fatigue_variables
from process.data_structure.current_drive_variables import init_current_drive_variables
from process.data_structure.dcll_variables import init_dcll_module
from process.data_structure.divertor_variables import init_divertor_variables
from process.data_structure.fwbs_variables import init_fwbs_variables
from process.data_structure.heat_transport_variables import (
    init_heat_transport_variables,
)
from process.data_structure.ife_variables import init_ife_variables
from process.data_structure.impurity_radiation_module import (
    init_impurity_radiation_module,
)
from process.data_structure.neoclassics_variables import init_neoclassics_variables
from process.data_structure.pf_power_variables import init_pf_power_variables
from process.data_structure.pfcoil_variables import (
    init_pfcoil_module,
    init_pfcoil_variables,
)
from process.data_structure.physics_variables import (
    init_physics_module,
    init_physics_variables,
)
from process.data_structure.power_variables import init_power_variables
from process.data_structure.primary_pumping_variables import (
    init_primary_pumping_variables,
)
from process.data_structure.pulse_variables import init_pulse_variables
from process.data_structure.rebco_variables import init_rebco_variables
from process.data_structure.reinke_variables import init_reinke_variables
from process.data_structure.scan_variables import init_scan_variables
from process.data_structure.stellarator_variables import init_stellarator_variables
from process.data_structure.structure_variables import init_structure_variables
from process.data_structure.superconducting_tf_coil_variables import (
    init_superconducting_tf_coil_variables,
)
from process.data_structure.tfcoil_variables import init_tfcoil_variables
from process.data_structure.times_variables import init_times_variables
from process.data_structure.vacuum_variables import init_vacuum_variables
from process.data_structure.water_usage_variables import init_watuse_variables
from process.exceptions import ProcessValidationError
from process.input import parse_input_file
from process.log import logging_model_handler
from process.models.stellarator.initialization import st_init


def init_process():
    """Routine that calls the initialisation routines

    This routine calls the main initialisation routines that set
    the default values for the global variables, reads in data from
    the input file, and checks the run parameters for consistency.
    """
    # Initialise the program variables
    iteration_variables.initialise_iteration_variables()

    # Creating and open the files MFile and OUTFile
    process_output.OutputFileManager.open_files()

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
    for outfile in [constants.NOUT, constants.IOTTY]:
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

        fileprefix = data_structure.global_variables.fileprefix
        process_output.ocmmnt(
            outfile,
            f"Input : {fileprefix}",
        )
        runtitle = data_structure.global_variables.runtitle
        process_output.ocmmnt(
            outfile,
            f"Run title : {runtitle}",
        )

        process_output.ocmmnt(
            outfile,
            f"Run type : Reactor concept design: {data_structure.global_variables.icase}, (c) UK Atomic Energy Authority",
        )

        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)
        process_output.oblnkl(outfile)

        process_output.ocmmnt(
            outfile, f"Equality constraints : {data_structure.numerics.neqns}"
        )
        process_output.ocmmnt(
            outfile,
            f"Inequality constraints : {data_structure.numerics.nineqns}",
        )
        process_output.ocmmnt(
            outfile,
            f"Total constraints : {data_structure.numerics.nineqns + data_structure.numerics.neqns}",
        )
        process_output.ocmmnt(
            outfile, f"Iteration variables : {data_structure.numerics.nvar}"
        )
        # If optimising, write objective function and convergence parameter
        if data_structure.numerics.ioptimz == 1:
            process_output.ocmmnt(
                outfile,
                f"Max iterations : {data_structure.global_variables.maxcal}",
            )

            if data_structure.numerics.minmax > 0:
                minmax_string = "  -- minimise "
                minmax_sign = "+"
            else:
                minmax_string = "  -- maximise "
                minmax_sign = "-"

            fom_string = data_structure.numerics.lablmm[
                abs(data_structure.numerics.minmax) - 1
            ]
            process_output.ocmmnt(
                outfile,
                f"Figure of merit : {minmax_sign}{abs(data_structure.numerics.minmax)}{minmax_string}{fom_string}",
            )
            process_output.ocmmnt(
                outfile,
                f"Convergence parameter : {data_structure.numerics.epsvmc}",
            )

        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)

    # MFile #
    mfile = constants.MFILE

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
        mfile, "Optimisation switch", "(ioptimz)", data_structure.numerics.ioptimz
    )
    # If optimising, write figure of merit switch
    if data_structure.numerics.ioptimz == 1:
        process_output.ovarin(
            mfile, "Figure of merit switch", "(minmax)", data_structure.numerics.minmax
        )


def init_all_module_vars():
    """Initialise all module variables
    This is vital to ensure a 'clean' state of Process before a new run starts,
    otherwise components of the previous run's state can persist into the new
    run. This matters ever since Process is used as a shared library, rather
    than a 'run-once' executable.
    """
    logging_model_handler.clear_logs()
    data_structure.numerics.init_numerics()
    init_buildings_variables()
    init_cost_variables()
    init_divertor_variables()
    init_fwbs_variables()
    data_structure.global_variables.init_global_variables()
    init_ccfe_hcpb_module()
    init_heat_transport_variables()
    init_ife_variables()
    init_impurity_radiation_module()
    init_pfcoil_module()
    init_physics_module()
    init_physics_variables()
    init_scan_variables()
    init_superconducting_tf_coil_variables()
    init_stellarator_variables()
    init_tfcoil_variables()
    init_times_variables()
    constants.init_constants()
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
    init_neoclassics_variables()


def check_process(inputs):  # noqa: ARG001
    """Routine to reset specific variables if certain options are
    being used

    This routine performs a sanity check of the input variables
    and ensures other dependent variables are given suitable values.
    """

    # Check that there are sufficient iteration variables
    if data_structure.numerics.nvar < data_structure.numerics.neqns:
        raise ProcessValidationError(
            "Insufficient iteration variables to solve the problem! NVAR < NEQNS",
            nvar=data_structure.numerics.nvar,
            neqns=data_structure.numerics.neqns,
        )

    # Check that sufficient elements of ixc and icc have been specified
    if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 0).any():
        raise ProcessValidationError(
            "The number of iteration variables specified is smaller than the number stated in ixc",
            nvar=data_structure.numerics.nvar,
        )

    # Check that dr_tf_wp_with_insulation (ixc = 140) and dr_tf_inboard (ixc = 13) are not being used simultaneously as iteration variables
    if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 13).any() and (
        data_structure.numerics.ixc[: data_structure.numerics.nvar] == 140
    ).any():
        raise ProcessValidationError(
            "Iteration variables 13 and 140 cannot be used simultaneously",
        )

    # Can't use c_tf_turn as interation var, constraint or input if i_tf_turns_integer == 1
    if (
        data_structure.numerics.ixc[: data_structure.numerics.nvar] == 60
    ).any() and data_structure.tfcoil_variables.i_tf_turns_integer == 1:
        raise ProcessValidationError(
            "Iteration variable 60 (TF current per turn, c_tf_turn) cannot be used with the TF coil integer turn model (i_tf_turns_integer == 1) as it is a calculated output instead for this model. However, the maximum current per turn can be constrained with constraint 77."
        )

    # Can't have icc 77 and ixc 60 at the same time
    if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 60).any() and (
        data_structure.numerics.icc[: data_structure.numerics.nvar] == 77
    ).any():
        raise ProcessValidationError(
            "Cannot use iteration variable 60 (TF coil current per turn, c_tf_turn) and constraint 77 (maximum TF current per turn) simultaneously."
        )

    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 0
    ).any():
        raise ProcessValidationError(
            "The number of constraints specified is smaller than the number stated in neqns+nineqns",
            neqns=data_structure.numerics.neqns,
            nineqns=data_structure.numerics.nineqns,
        )

    # Deprecate constraints
    for depcrecated_constraint in [3, 4, 10, 74, 42]:
        if (
            data_structure.numerics.icc[
                : data_structure.numerics.neqns + data_structure.numerics.nineqns
            ]
            == depcrecated_constraint
        ).any():
            raise ProcessValidationError(
                "Constraint equation is no longer available", icc=depcrecated_constraint
            )

    # MDK Report error if constraint 63 is used with old vacuum model
    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 63
    ).any() and data_structure.vacuum_variables.i_vacuum_pumping != "simple":
        raise ProcessValidationError(
            "Constraint 63 is requested without the correct vacuum model (simple)"
        )

    #  Fuel ion fractions must add up to 1.0
    if (
        abs(
            1.0
            - data_structure.physics_variables.f_plasma_fuel_deuterium
            - data_structure.physics_variables.f_plasma_fuel_tritium
            - data_structure.physics_variables.f_plasma_fuel_helium3
        )
        > 1e-6
    ):
        raise ProcessValidationError(
            "Fuel ion fractions do not sum to 1.0",
            f_plasma_fuel_deuterium=data_structure.physics_variables.f_plasma_fuel_deuterium,
            f_plasma_fuel_tritium=data_structure.physics_variables.f_plasma_fuel_tritium,
            f_plasma_fuel_helium3=data_structure.physics_variables.f_plasma_fuel_helium3,
        )

    if (
        data_structure.physics_variables.f_plasma_fuel_tritium < 1.0e-3
    ):  # tritium fraction is negligible
        data_structure.buildings_variables.triv = 0.0
        data_structure.heat_transport_variables.p_tritium_plant_electric_mw = 0.0

    if data_structure.impurity_radiation_module.f_nd_impurity_electrons[1] != 0.1:
        raise ProcessValidationError(
            "The thermal alpha/electron density ratio should be controlled using f_nd_alpha_electron (itv 109) and not f_nd_impurity_electrons(2)."
            "f_nd_impurity_electrons(2) should be removed from the input file, or set to the default value 0.1D0."
        )

    # Impurity fractions
    for imp in range(data_structure.impurity_radiation_module.N_IMPURITIES):
        data_structure.impurity_radiation_module.f_nd_impurity_electron_array[imp] = (
            data_structure.impurity_radiation_module.f_nd_impurity_electrons[imp]
        )

    # Stop the run if oacdcp is used as an optimisation variable
    # As the current density is now calculated from b_plasma_toroidal_on_axis without constraint 10

    if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 12).any():
        raise ProcessValidationError(
            "The 1/R toroidal B field dependency constraint is being depreciated"
        )

    # Plasma profile consistency checks
    if (
        data_structure.ife_variables.ife != 1
        and data_structure.physics_variables.i_plasma_pedestal == 1
    ):
        # Temperature checks
        if (
            data_structure.physics_variables.temp_plasma_pedestal_kev
            < data_structure.physics_variables.temp_plasma_separatrix_kev
        ):
            raise ProcessValidationError(
                "Pedestal temperature is lower than separatrix temperature",
                temp_plasma_pedestal_kev=data_structure.physics_variables.temp_plasma_pedestal_kev,
                temp_plasma_separatrix_kev=data_structure.physics_variables.temp_plasma_separatrix_kev,
            )

        if (
            abs(data_structure.physics_variables.radius_plasma_pedestal_temp_norm - 1.0)
            <= 1e-7
        ) and (
            (
                data_structure.physics_variables.temp_plasma_pedestal_kev
                - data_structure.physics_variables.temp_plasma_separatrix_kev
            )
            >= 1e-7
        ):
            warn(
                f"Temperature pedestal is at plasma edge, but temp_plasma_pedestal_kev "
                f"({data_structure.physics_variables.temp_plasma_pedestal_kev}) differs from temp_plasma_separatrix_kev "
                f"({data_structure.physics_variables.temp_plasma_separatrix_kev})",
                stacklevel=2,
            )

        # Core temperature should always be calculated (later) as being
        # higher than the pedestal temperature, if and only if the
        # volume-averaged temperature never drops below the pedestal
        # temperature. Prevent this by adjusting te, and its lower bound
        # (which will only have an effect if this is an optimisation run)
        if (
            data_structure.physics_variables.temp_plasma_electron_vol_avg_kev
            <= data_structure.physics_variables.temp_plasma_pedestal_kev
        ):
            warn(
                f"Volume-averaged temperature ({data_structure.physics_variables.te}) has been "
                f"forced to exceed input pedestal height ({data_structure.physics_variables.temp_plasma_pedestal_kev}). "
                "Changing to te = temp_plasma_pedestal_kev*1.001",
                stacklevel=2,
            )
            data_structure.physics_variables.temp_plasma_electron_vol_avg_kev = (
                data_structure.physics_variables.temp_plasma_pedestal_kev * 1.001
            )

        if (
            data_structure.numerics.ioptimz >= 0
            and (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 4).any()
            and data_structure.numerics.boundl[3]
            < data_structure.physics_variables.temp_plasma_pedestal_kev * 1.001
        ):
            warn(
                "Lower limit of volume averaged electron temperature (temp_plasma_electron_vol_avg_kev) has been raised to ensure temp_plasma_electron_vol_avg_kev > temp_plasma_pedestal_kev",
                stacklevel=2,
            )
            data_structure.numerics.boundl[3] = (
                data_structure.physics_variables.temp_plasma_pedestal_kev * 1.001
            )
            data_structure.numerics.boundu[3] = max(
                data_structure.numerics.boundu[3], data_structure.numerics.boundl[3]
            )

        # Density checks
        # Case where pedestal density is set manually
        if (
            data_structure.physics_variables.f_nd_plasma_pedestal_greenwald < 0
            or not (
                data_structure.numerics.ixc[: data_structure.numerics.nvar] == 145
            ).any()
        ):
            # Issue #589 Pedestal density is set manually using nd_plasma_pedestal_electron but it is less than nd_plasma_separatrix_electron.
            if (
                data_structure.physics_variables.nd_plasma_pedestal_electron
                < data_structure.physics_variables.nd_plasma_separatrix_electron
            ):
                raise ProcessValidationError(
                    "Density pedestal is lower than separatrix density",
                    nd_plasma_pedestal_electron=data_structure.physics_variables.nd_plasma_pedestal_electron,
                    nd_plasma_separatrix_electron=data_structure.physics_variables.nd_plasma_separatrix_electron,
                )

            # Issue #589 Pedestal density is set manually using nd_plasma_pedestal_electron,
            # but pedestal width = 0.
            if (
                abs(
                    data_structure.physics_variables.radius_plasma_pedestal_density_norm
                    - 1.0
                )
                <= 1e-7
                and (
                    data_structure.physics_variables.nd_plasma_pedestal_electron
                    - data_structure.physics_variables.nd_plasma_separatrix_electron
                )
                >= 1e-7
            ):
                warn(
                    "Density pedestal is at plasma edge "
                    f"({data_structure.physics_variables.radius_plasma_pedestal_density_norm = }), but nd_plasma_pedestal_electron "
                    f"({data_structure.physics_variables.nd_plasma_pedestal_electron}) differs from "
                    f"nd_plasma_separatrix_electron ({data_structure.physics_variables.nd_plasma_separatrix_electron})",
                    stacklevel=2,
                )

        # Issue #862 : Variable nd_plasma_electron_on_axis/nd_plasma_pedestal_electron ratio without constraint eq 81 (nd_plasma_electron_on_axis>nd_plasma_pedestal_electron)
        #  -> Potential hollowed density profile
        if (
            data_structure.numerics.ioptimz >= 0
            and not (
                data_structure.numerics.icc[
                    : data_structure.numerics.neqns + data_structure.numerics.nineqns
                ]
                == 81
            ).any()
        ):
            if (
                data_structure.numerics.ixc[: data_structure.numerics.nvar] == 145
            ).any():
                warn(
                    "nd_plasma_pedestal_electron set with f_nd_plasma_pedestal_greenwald without constraint eq 81 (nd_plasma_pedestal_electron<nd_plasma_electron_on_axis)",
                    stacklevel=2,
                )
            if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 6).any():
                warn(
                    "nd_plasma_electrons_vol_avg used as iteration variable without constraint 81 (nd_plasma_pedestal_electron<nd_plasma_electron_on_axis)",
                    stacklevel=2,
                )

    # Cannot use Psep/R and PsepB/qAR limits at the same time
    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 68
    ).any() and (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 56
    ).any():
        raise ProcessValidationError(
            "Cannot use Psep/R and PsepB/qAR constraint equations at the same time"
        )

    # if lower bound of f_nd_plasma_pedestal_greenwald < f_nd_plasma_separatrix_greenwald
    if (
        data_structure.numerics.ixc[: data_structure.numerics.nvar] == 145
    ).any() and data_structure.numerics.boundl[
        144
    ] < data_structure.physics_variables.f_nd_plasma_separatrix_greenwald:
        raise ProcessValidationError(
            "Set lower bound of iteration variable 145, f_nd_plasma_pedestal_greenwald, to be greater than f_nd_plasma_separatrix_greenwald",
            boundl_145=data_structure.numerics.boundl[144],
            f_nd_plasma_separatrix_greenwald=data_structure.physics_variables.f_nd_plasma_separatrix_greenwald,
        )

    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 78
    ).any():
        # If Reinke criterion is used temp_plasma_separatrix_kev is calculated and cannot be an
        # iteration variable
        if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 119).any():
            raise ProcessValidationError(
                "REINKE IMPURITY MODEL: temp_plasma_separatrix_kev is calculated and cannot be an "
                "iteration variable for the Reinke model"
            )

        # If Reinke criterion is used need to enforce LH-threshold
        # using Martin scaling for consistency
        if (data_structure.physics_variables.i_l_h_threshold != 6) or (
            not (
                data_structure.numerics.icc[
                    : data_structure.numerics.neqns + data_structure.numerics.nineqns
                ]
                == 15
            ).any()
            and data_structure.physics_variables.i_plasma_pedestal
        ):
            warn(
                "REINKE IMPURITY MODEL: The Martin LH threshold scale is not being used and is recommended for the Reinke model",
                stacklevel=2,
            )

    if data_structure.physics_variables.i_single_null == 0:
        data_structure.divertor_variables.n_divertors = 2
        data_structure.build_variables.dz_fw_plasma_gap = (
            data_structure.build_variables.dz_xpoint_divertor
        )
        data_structure.build_variables.dz_shld_upper = (
            data_structure.build_variables.dz_shld_lower
        )
        data_structure.build_variables.dz_vv_upper = (
            data_structure.build_variables.dz_vv_lower
        )
        warn("Double-null: Upper vertical build forced to match lower", stacklevel=2)
    else:  # i_single_null == 1
        data_structure.divertor_variables.n_divertors = 1

    #  Tight aspect ratio options (ST)
    if data_structure.physics_variables.itart == 1:
        data_structure.global_variables.icase = "Tight aspect ratio tokamak model"

        # Disabled Forcing that no inboard breeding blanket is used
        # Disabled i_blkt_inboard = 0

        # Check if the choice of plasma current is addapted for ST
        # 2 : Peng Ip scaling (See STAR code documentation)
        # 9 : Fiesta Ip scaling
        if (
            data_structure.physics_variables.i_plasma_current != 2
            and data_structure.physics_variables.i_plasma_current != 9
        ):
            warn(
                "Usual current scaling for TARTs (i_plasma_current=2 or 9) is not being used",
                stacklevel=2,
            )

        # If using Peng and Strickler (1986) model (itartpf == 0)
        # Overwrite the location of the TF coils
        # 2 : PF coil on top of TF coil
        # 3 : PF coil outside of TF coil
        if data_structure.physics_variables.itartpf == 0:
            data_structure.pfcoil_variables.i_pf_location[0] = 2
            data_structure.pfcoil_variables.i_pf_location[1] = 3
            data_structure.pfcoil_variables.i_pf_location[2] = 3

        # Water cooled copper magnets initalisation / checks
        if data_structure.tfcoil_variables.i_tf_sup == 0:
            # Check if the initial centrepost coolant loop adapted to the magnet technology
            # Ice cannot flow so temp_cp_coolant_inlet > 273.15 K
            if data_structure.tfcoil_variables.temp_cp_coolant_inlet < 273.15:
                raise ProcessValidationError(
                    "Coolant temperature (temp_cp_coolant_inlet) cannot be < 0 C (273.15 K) for water cooled copper magents"
                )

            # Temperature of the TF legs cannot be cooled down
            if (
                data_structure.tfcoil_variables.temp_tf_legs_outboard > 0
                and data_structure.tfcoil_variables.temp_tf_legs_outboard < 273.15
            ):
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) cannot be < 0 C (273.15 K) for water cooled magents"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                data_structure.numerics.ixc[: data_structure.numerics.nvar] == 20
            ).any() and data_structure.numerics.boundu[19] < 273.15:
                raise ProcessValidationError(
                    "Too low CP conductor temperature (temp_cp_average). Lower limit for copper > 273.15 K"
                )

        # Call a lvl 3 error if superconductor magnets are used
        elif data_structure.tfcoil_variables.i_tf_sup == 1:
            warn(
                "Joints res not cal. for SC (itart = 1) TF (data_structure.tfcoil_variables.i_tf_sup = 1)",
                stacklevel=2,
            )

        # Aluminium magnets initalisation / checks
        # Initialize the CP conductor temperature to cryogenic temperature for cryo-al magnets (20 K)
        elif data_structure.tfcoil_variables.i_tf_sup == 2:
            # Call a lvl 3 error if the inlet coolant temperature is too large
            # Motivation : ill-defined aluminium resistivity fit for T > 40-50 K
            if data_structure.tfcoil_variables.temp_cp_coolant_inlet > 40.0:
                raise ProcessValidationError(
                    "Coolant temperature (temp_cp_coolant_inlet) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if the leg average temperature is low enough for the resisitivity fit
            if data_structure.tfcoil_variables.temp_tf_legs_outboard > 50.0:
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                data_structure.numerics.ixc[: data_structure.numerics.nvar] == 20
            ).any() and data_structure.numerics.boundu[19] > 50.0:
                raise ProcessValidationError(
                    "Too large CP conductor temperature (temp_cp_average). Upper limit for cryo-al < 50 K"
                )

            # Otherwise intitialise the average conductor temperature at
            data_structure.tfcoil_variables.temp_cp_average = (
                data_structure.tfcoil_variables.temp_cp_coolant_inlet
            )

        # Check if the boostrap current selection is addapted to ST
        if data_structure.physics_variables.i_bootstrap_current == 1:
            raise ProcessValidationError(
                "Invalid boostrap current law for ST, do not use i_bootstrap_current = 1"
            )

        # Check if a single null divertor is used in double null machine
        if data_structure.physics_variables.i_single_null == 0 and (
            data_structure.physics_variables.f_p_div_lower == 1.0
            or data_structure.physics_variables.f_p_div_lower == 0.0
        ):
            warn("Operating with a single null in a double null machine", stacklevel=2)

        # Set the TF coil shape to picture frame (if default value)
        if data_structure.tfcoil_variables.i_tf_shape == 0:
            data_structure.tfcoil_variables.i_tf_shape = 2

        # Warning stating that the CP fast neutron fluence calculation
        # is not addapted for cryoaluminium calculations yet
        if (
            data_structure.tfcoil_variables.i_tf_sup == 2
            and (
                data_structure.numerics.icc[
                    : data_structure.numerics.neqns + data_structure.numerics.nineqns
                ]
                == 85
            ).any()
            and data_structure.physics_variables.itart == 1
        ):
            raise ProcessValidationError(
                "Al TF coil fluence not calculated properly for Al CP, do not use constraint 85"
            )

        # Setting the CP joints default options :
        #  0 : No joints for superconducting magents (data_structure.tfcoil_variables.i_tf_sup = 1)
        #  1 : Sliding joints for resistive magnets (data_structure.tfcoil_variables.i_tf_sup = 0, 2)
        if data_structure.tfcoil_variables.i_cp_joints == -1:
            if data_structure.tfcoil_variables.i_tf_sup == 1:
                data_structure.tfcoil_variables.i_cp_joints = 0
            else:
                data_structure.tfcoil_variables.i_cp_joints = 1

        # Checking the CP TF top radius
        if (
            abs(data_structure.build_variables.r_cp_top) > 0
            or (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 174).any()
        ) and data_structure.build_variables.i_r_cp_top != 1:
            raise ProcessValidationError(
                "To set the TF CP top value, you must use i_r_cp_top = 1"
            )

    # Conventionnal aspect ratios specific
    else:
        if (
            data_structure.physics_variables.i_plasma_current == 2
            or data_structure.physics_variables.i_plasma_current == 9
        ):
            raise ProcessValidationError(
                "i_plasma_current=2,9 is not a valid option for a non-TART device"
            )

        # Set the TF coil shape to PROCESS D-shape (if default value)
        if data_structure.tfcoil_variables.i_tf_shape == 0:
            data_structure.tfcoil_variables.i_tf_shape = 1

        # Check PF coil configurations
        j = 0
        k = 0
        for i in range(data_structure.pfcoil_variables.n_pf_coil_groups):
            if (
                data_structure.pfcoil_variables.i_pf_location[i] != 2
                and data_structure.pfcoil_variables.n_pf_coils_in_group[i] != 2
            ):
                raise ProcessValidationError(
                    "n_pf_coils_in_group(i) .ne. 2 is not a valid option except for (i_pf_location = 2)"
                )

            if data_structure.pfcoil_variables.i_pf_location[i] == 2:
                j = j + 1
                k = k + data_structure.pfcoil_variables.n_pf_coils_in_group[i]

        if k == 1:
            raise ProcessValidationError(
                "Only 1 divertor coil (i_pf_location = 2) is not a valid configuration"
            )
        if k > 2:
            raise ProcessValidationError(
                "More than 2 divertor coils (i_pf_location = 2) is not a valid configuration"
            )
        if data_structure.physics_variables.i_single_null == 1 and j < 2:
            raise ProcessValidationError(
                "If i_single_null=1, use 2 individual divertor coils (i_pf_location = 2, 2; n_pf_coils_in_group = 1, 1)"
            )

        # Constraint 10 is dedicated to ST designs with demountable joints
        if (
            data_structure.numerics.icc[
                : data_structure.numerics.neqns + data_structure.numerics.nineqns
            ]
            == 10
        ).any():
            raise ProcessValidationError(
                "Constraint equation 10 (CP lifetime) to used with ST desing (itart=1)"
            )

    #  Pulsed power plant model
    if data_structure.pulse_variables.i_pulsed_plant == 1:
        data_structure.global_variables.icase = "Pulsed tokamak model"
    else:
        data_structure.buildings_variables.esbldgm3 = 0.0

    # TF coil
    # -------
    # TF stress model not defined of r_tf_inboard = 0
    # Unless i_tf_stress_model == 2
    # -> If dr_bore + dr_cs_tf_gap + dr_cs = 0 and fixed and stress constraint is used
    #    Generate a lvl 3 error proposing not to use any stress constraints
    if (
        (
            not (
                (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 16).any()
                or (
                    data_structure.numerics.ixc[: data_structure.numerics.nvar] == 29
                ).any()
                or (
                    data_structure.numerics.ixc[: data_structure.numerics.nvar] == 42
                ).any()
            )
        )  # No dr_bore,dr_cs_tf_gap, dr_cs iteration
        and (
            abs(
                data_structure.build_variables.dr_bore
                + data_structure.build_variables.dr_cs_tf_gap
                + data_structure.build_variables.dr_cs
                + data_structure.build_variables.dr_cs_precomp
            )
            <= 0
        )  # dr_bore + dr_cs_tf_gap + dr_cs = 0
        and (
            (
                data_structure.numerics.icc[
                    : data_structure.numerics.neqns + data_structure.numerics.nineqns
                ]
                == 31
            ).any()
            or (
                data_structure.numerics.icc[
                    : data_structure.numerics.neqns + data_structure.numerics.nineqns
                ]
                == 32
            ).any()
        )  # Stress constraints (31 or 32) is used
        and (
            data_structure.tfcoil_variables.i_tf_stress_model != 2
        )  # TF stress model can't handle no dr_bore
    ):
        raise ProcessValidationError(
            "Invalid stress model if dr_bore + dr_cs_tf_gap + dr_cs = 0. Don't use constraint 31"
        )

    # Make sure that plane stress model is not used for resistive magnets
    if (
        data_structure.tfcoil_variables.i_tf_stress_model == 1
        and data_structure.tfcoil_variables.i_tf_sup != 1
    ):
        raise ProcessValidationError(
            "Use generalized plane strain for resistive magnets (i_tf_stress_model = 0 or 2)"
        )

    # bucking cylinder default option setting
    # - bucking (casing) for SC i_tf_bucking ( i_tf_bucking = 1 )
    # - No bucking for copper magnets ( i_tf_bucking = 0 )
    # - Bucking for aluminium magnets ( i_tf_bucking = 1 )
    if data_structure.tfcoil_variables.i_tf_bucking == -1:
        if data_structure.tfcoil_variables.i_tf_sup == 0:
            data_structure.tfcoil_variables.i_tf_bucking = 0
        else:
            data_structure.tfcoil_variables.i_tf_bucking = 1

    # Ensure that the TF isnt placed against the
    # CS which is now outside it
    if (
        data_structure.tfcoil_variables.i_tf_bucking >= 2
        and data_structure.build_variables.i_tf_inside_cs == 1
    ):
        raise ProcessValidationError(
            "Cannot have i_tf_bucking >= 2 when i_tf_inside_cs = 1"
        )

    # Ensure that no pre-compression structure
    # is used for bucked and wedged design
    if (
        data_structure.tfcoil_variables.i_tf_bucking >= 2
        and data_structure.build_variables.i_cs_precomp == 1
    ):
        raise ProcessValidationError(
            "No CS precompression structure for bucked and wedged, use i_cs_precomp = 0"
        )

    # Number of stress calculation layers
    # +1 to add in the inboard TF coil case on the plasma side, per Issue #1509
    data_structure.tfcoil_variables.n_tf_stress_layers = (
        data_structure.tfcoil_variables.i_tf_bucking
        + data_structure.tfcoil_variables.n_tf_graded_layers
        + 1
    )

    # If TFC sidewall has not been set by user
    if data_structure.tfcoil_variables.dx_tf_side_case_min < 0.1e-10:
        data_structure.tfcoil_variables.tfc_sidewall_is_fraction = True

    # If inboard TF coil case plasma side thickness has not been set by user
    if data_structure.tfcoil_variables.dr_tf_plasma_case < 0.1e-10:
        data_structure.tfcoil_variables.i_f_dr_tf_plasma_case = True

    # Setting the default cryo-plants efficiencies
    if abs(data_structure.tfcoil_variables.eff_tf_cryo + 1) < 1e-6:
        # The ITER cyoplant efficiency is used for SC
        if data_structure.tfcoil_variables.i_tf_sup == 1:
            data_structure.tfcoil_variables.eff_tf_cryo = 0.13

        # Strawbrige plot extrapolation is used for Cryo-Al
        elif data_structure.tfcoil_variables.i_tf_sup == 2:
            data_structure.tfcoil_variables.eff_tf_cryo = 0.40

    # Cryo-plane efficiency must be in [0-1.0]
    elif (
        data_structure.tfcoil_variables.eff_tf_cryo > 1.0
        or data_structure.tfcoil_variables.eff_tf_cryo < 0.0
    ):
        raise ProcessValidationError(
            "TF cryo-plant efficiency `eff_tf_cryo` must be within [0-1]"
        )

    # Integer turns option not yet available for REBCO taped turns

    if (
        data_structure.tfcoil_variables.i_tf_sc_mat == 6
        and data_structure.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Integer turns (i_tf_turns_integer = 1) not supported for REBCO (i_tf_sc_mat = 6)"
        )

    # Setting up insulation layer young modulae default values [Pa]

    if data_structure.tfcoil_variables.eyoung_ins <= 1.0e8:
        # Copper magnets, no insulation material defined
        # But use the ITER design by default
        if (
            data_structure.tfcoil_variables.i_tf_sup == 0
            or data_structure.tfcoil_variables.i_tf_sup == 1
        ):
            # SC magnets
            # Value from DDD11-2 v2 2 (2009)
            data_structure.tfcoil_variables.eyoung_ins = 20.0e9

        # Cryo-aluminum magnets (Kapton polymer)
        elif data_structure.tfcoil_variables.i_tf_sup == 2:
            data_structure.tfcoil_variables.eyoung_ins = 2.5e9

    # Setting the default WP geometry

    if data_structure.tfcoil_variables.i_tf_wp_geom == -1:
        if data_structure.tfcoil_variables.i_tf_turns_integer == 0:
            data_structure.tfcoil_variables.i_tf_wp_geom = 1
        if data_structure.tfcoil_variables.i_tf_turns_integer == 1:
            data_structure.tfcoil_variables.i_tf_wp_geom = 0

    # Setting the TF coil conductor elastic properties

    if data_structure.tfcoil_variables.i_tf_cond_eyoung_axial == 0:
        # Conductor stiffness is not considered
        data_structure.tfcoil_variables.eyoung_cond_axial = 0
        data_structure.tfcoil_variables.eyoung_cond_trans = 0
    elif data_structure.tfcoil_variables.i_tf_cond_eyoung_axial == 2:
        # Select sensible defaults from the literature
        if data_structure.tfcoil_variables.i_tf_sc_mat in [1, 4, 5]:
            # Nb3Sn: Nyilas, A et. al, Superconductor Science and Technology 16, no. 9 (2003): 1036-42. https://doi.org/10.1088/0953-2048/16/9/313.
            data_structure.tfcoil_variables.eyoung_cond_axial = 32e9
        elif data_structure.tfcoil_variables.i_tf_sc_mat == 2:
            # Bi-2212: Brown, M. et al, IOP Conference Series: Materials Science and Engineering 279 (2017): 012022. https://doi.org/10.1088/1757-899X/279/1/012022.
            data_structure.tfcoil_variables.eyoung_cond_axial = 80e9
        elif data_structure.tfcoil_variables.i_tf_sc_mat in [3, 7]:
            # NbTi: Vedrine, P. et. al, IEEE Transactions on Applied Superconductivity 9, no. 2 (1999): 236-39. https://doi.org/10.1109/77.783280.
            data_structure.tfcoil_variables.eyoung_cond_axial = 6.8e9
        elif data_structure.tfcoil_variables.i_tf_sc_mat in [6, 8, 9]:
            # REBCO: Fujishiro, H. et. al, Physica C: Superconductivity, 426-431 (2005): 699-704. https://doi.org/10.1016/j.physc.2005.01.045.
            data_structure.tfcoil_variables.eyoung_cond_axial = 145e9

        if data_structure.tfcoil_variables.i_tf_cond_eyoung_trans == 0:
            # Transverse stiffness is not considered
            data_structure.tfcoil_variables.eyoung_cond_trans = 0
        else:
            # Transverse stiffness is significant
            data_structure.tfcoil_variables.eyoung_cond_trans = (
                data_structure.tfcoil_variables.eyoung_cond_axial
            )

    # Check if the WP/conductor radial thickness (dr_tf_wp_with_insulation) is large enough
    # To contains the insulation, cooling and the support structure
    # Rem : Only verified if the WP thickness is used
    if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 140).any():
        # Minimal WP thickness
        if data_structure.tfcoil_variables.i_tf_sup == 1:
            dr_tf_wp_min = 2.0 * (
                data_structure.tfcoil_variables.dx_tf_wp_insulation
                + data_structure.tfcoil_variables.dx_tf_wp_insertion_gap
                + data_structure.tfcoil_variables.dx_tf_turn_insulation
                + data_structure.tfcoil_variables.dia_tf_turn_coolant_channel
            )

            # Steel conduit thickness (can be an iteration variable)
            if (data_structure.numerics.ixc[: data_structure.numerics.nvar] == 58).any():
                dr_tf_wp_min = dr_tf_wp_min + 2.0 * data_structure.numerics.boundl[57]
            else:
                dr_tf_wp_min = (
                    dr_tf_wp_min + 2.0 * data_structure.tfcoil_variables.dx_tf_turn_steel
                )

        # Minimal conductor layer thickness
        elif (
            data_structure.tfcoil_variables.i_tf_sup == 0
            or data_structure.tfcoil_variables.i_tf_sup == 2
        ):
            dr_tf_wp_min = (
                2.0
                * (
                    data_structure.tfcoil_variables.dx_tf_turn_insulation
                    + data_structure.tfcoil_variables.dx_tf_wp_insulation
                )
                + 4.0 * data_structure.tfcoil_variables.radius_cp_coolant_channel
            )

        if data_structure.numerics.boundl[139] < dr_tf_wp_min:
            raise ProcessValidationError(
                "The TF coil WP thickness (dr_tf_wp_with_insulation) must be at least",
                dr_tf_wp_min=dr_tf_wp_min,
            )

    # Setting i_dx_tf_turn_general_input to true if dx_tf_turn_general is an input
    data_structure.tfcoil_variables.i_dx_tf_turn_general_input = (
        abs(data_structure.tfcoil_variables.dx_tf_turn_general) > 0
    )

    # Impossible to set the turn size of integer turn option
    if (
        data_structure.tfcoil_variables.i_dx_tf_turn_general_input
        and data_structure.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    if (
        data_structure.tfcoil_variables.i_tf_wp_geom != 0
        and data_structure.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Can only have i_tf_turns_integer = 1 with i_tf_wp_geom = 0"
        )

    if (
        data_structure.physics_variables.i_bootstrap_current == 5
        and data_structure.physics_variables.i_diamagnetic_current != 0
    ):
        raise ProcessValidationError(
            "i_diamagnetic_current = 0 should be used with the Sakai plasma current scaling"
        )

    # Setting i_dx_tf_turn_cable_space_general_input to true if dx_tf_turn_cable_space_general is an input
    data_structure.tfcoil_variables.i_dx_tf_turn_cable_space_general_input = (
        abs(data_structure.tfcoil_variables.dx_tf_turn_cable_space_general) > 0
    )

    # Impossible to set the cable size of integer turn option
    if (
        data_structure.tfcoil_variables.i_dx_tf_turn_cable_space_general_input
        and data_structure.tfcoil_variables.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    # Impossible to set both the TF coil turn and the cable dimension
    if (
        data_structure.tfcoil_variables.i_dx_tf_turn_general_input
        and data_structure.tfcoil_variables.i_dx_tf_turn_cable_space_general_input
    ):
        raise ProcessValidationError(
            "Impossible to set the TF coil turn and cable size simultaneously"
        )

    # Checking the SC temperature for LTS
    if (
        data_structure.tfcoil_variables.i_tf_sc_mat in [1, 3, 4, 5]
        and data_structure.tfcoil_variables.tftmp > 10.0
    ):
        raise ProcessValidationError(
            "The LTS conductor temperature (tftmp) has to be lower than 10"
        )

    # PF coil resistivity is zero if superconducting
    if data_structure.pfcoil_variables.i_pf_conductor == 0:
        data_structure.pfcoil_variables.rho_pf_coil = 0.0

    # If there is no NBI, then hot beam density should be zero
    if data_structure.current_drive_variables.i_hcd_calculations == 1:
        if (
            data_structure.current_drive_variables.i_hcd_primary != 5
            and data_structure.current_drive_variables.i_hcd_primary != 8
        ):
            data_structure.physics_variables.f_nd_beam_electron = 0.0
    else:
        data_structure.physics_variables.f_nd_beam_electron = 0.0

    # Set inboard blanket thickness to zero if no inboard blanket switch
    # used (Issue #732)
    if data_structure.fwbs_variables.i_blkt_inboard == 0:
        data_structure.build_variables.dr_blkt_inboard = 0.0

    # Ensure that blanket material fractions allow non-zero space for steel
    # CCFE HCPB Model

    if data_structure.stellarator_variables.istell == 0 and (
        data_structure.fwbs_variables.i_blanket_type == 1
    ):
        fsum = (
            data_structure.fwbs_variables.breeder_multiplier
            + data_structure.fwbs_variables.vfcblkt
            + data_structure.fwbs_variables.vfpblkt
        )
        if fsum >= 1.0:
            raise ProcessValidationError(
                "Blanket material fractions do not sum to 1.0",
                i_blanket_type=data_structure.fwbs_variables.i_blanket_type,
                breeder_multiplier=data_structure.fwbs_variables.breeder_multiplier,
                vfcblkt=data_structure.fwbs_variables.vfcblkt,
                vfpblkt=data_structure.fwbs_variables.vfpblkt,
                fsum=fsum,
            )

    # Check that the temperature margins are not overdetermined
    if data_structure.tfcoil_variables.tmargmin > 0.0001:
        # This limit has been input and will be applied to both TFC and CS
        if data_structure.tfcoil_variables.temp_tf_superconductor_margin_min > 0.0001:
            warn(
                "temp_tf_superconductor_margin_min and tmargmin should not both be specified in IN.DAT "
                "temp_tf_superconductor_margin_min has been ignored",
                stacklevel=2,
            )
        if data_structure.tfcoil_variables.temp_cs_superconductor_margin_min > 0.0001:
            warn(
                "temp_cs_superconductor_margin_min and tmargmin should not both be specified in IN.DAT "
                "temp_cs_superconductor_margin_min has been ignored",
                stacklevel=2,
            )

        data_structure.tfcoil_variables.temp_tf_superconductor_margin_min = (
            data_structure.tfcoil_variables.tmargmin
        )
        data_structure.tfcoil_variables.temp_cs_superconductor_margin_min = (
            data_structure.tfcoil_variables.tmargmin
        )

    if (
        data_structure.physics_variables.tauee_in > 1e-10
        and data_structure.physics_variables.i_confinement_time != 48
    ):
        # Report error if confinement time is in the input
        # but the scaling to use it is not selected.
        warn("tauee_in is for use with i_confinement_time=48 only", stacklevel=2)

    if (
        data_structure.physics_variables.aspect > 1.7
        and data_structure.physics_variables.i_confinement_time == 46
    ):
        # NSTX scaling is for A<1.7
        warn("NSTX scaling is for A<1.7", stacklevel=2)

    if (
        data_structure.physics_variables.i_plasma_current == 2
        and data_structure.physics_variables.i_confinement_time == 42
    ):
        raise ProcessValidationError(
            "Lang 2012 confinement scaling cannot be used for i_plasma_current=2 due to wrong q"
        )

    # Cannot use temperature margin constraint with REBCO TF coils
    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 36
    ).any() and (
        data_structure.tfcoil_variables.i_tf_sc_mat == 8
        or data_structure.tfcoil_variables.i_tf_sc_mat == 9
    ):
        raise ProcessValidationError(
            "turn off TF temperature margin constraint icc = 36 when using REBCO"
        )

    # Cannot use temperature margin constraint with REBCO CS coils
    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 60
    ).any() and data_structure.pfcoil_variables.i_cs_superconductor == 8:
        raise ProcessValidationError(
            "turn off CS temperature margin constraint icc = 60 when using REBCO"
        )

    # Cold end of the cryocooler should be colder than the TF
    if (
        data_structure.tfcoil_variables.temp_tf_cryo
        > data_structure.tfcoil_variables.tftmp
    ):
        raise ProcessValidationError("temp_tf_cryo should be lower than tftmp")

    # Cannot use TF coil strain limit if i_str_wp is off:
    if (
        data_structure.numerics.icc[
            : data_structure.numerics.neqns + data_structure.numerics.nineqns
        ]
        == 88
    ).any() and data_structure.tfcoil_variables.i_str_wp == 0:
        raise ProcessValidationError("Can't use constraint 88 if i_strain_tf == 0")


def set_active_constraints():
    """Set constraints provided in the input file as 'active'"""
    num_constraints = 0
    for i in range(ConstraintManager.num_constraints()):
        if data_structure.numerics.icc[i] != 0:
            data_structure.numerics.active_constraints[
                data_structure.numerics.icc[i] - 1
            ] = True
            num_constraints += 1

    if data_structure.numerics.neqns < 0:
        # The value of neqns has not been set in the input file.  Default = 0.
        data_structure.numerics.neqns = num_constraints - data_structure.numerics.nineqns
    else:
        data_structure.numerics.nineqns = num_constraints - data_structure.numerics.neqns


def set_device_type():
    if data_structure.ife_variables.ife == 1:
        data_structure.global_variables.icase = "Inertial Fusion model"
    elif data_structure.stellarator_variables.istell != 0:
        data_structure.global_variables.icase = "Stellarator model"
