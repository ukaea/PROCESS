from __future__ import annotations

import datetime
import getpass
import socket
import subprocess  # noqa: S404
from pathlib import Path
from typing import TYPE_CHECKING
from warnings import warn

import process
from process.core import constants, process_output
from process.core.exceptions import ProcessValidationError
from process.core.input import parse_input_file
from process.core.log import logging_model_handler
from process.core.solver import iteration_variables
from process.core.solver.constraints import ConstraintManager
from process.data_structure.blanket_variables import BlktModelTypes
from process.data_structure.build_variables import (
    InboardBlanketConfiguration,
    TFCSRadialConfiguration,
)
from process.data_structure.impurity_radiation_variables import N_IMPURITIES
from process.data_structure.numerics import FiguresOfMerit, PROCESSRunMode
from process.data_structure.pfcoil_variables import PFConductorModel
from process.data_structure.physics_variables import DivertorNumberModels
from process.models.pfcoil import PFLocationTypes
from process.models.physics.profiles import DensityProfilePedestalType
from process.models.stellarator.initialization import st_init
from process.models.superconductors import (
    SuperconductorMaterial,
    SuperconductorModel,
    SuperconductorType,
)
from process.models.tfcoil.base import TFCoilShapeModel, TFConductorModel
from process.models.tfcoil.superconducting import SuperconductingTFWPShapeType

if TYPE_CHECKING:
    from process.core.model import DataStructure


def init_process(data: DataStructure):
    """Routine that calls the initialisation routines

    This routine calls the main initialisation routines that set
    the default values for the global variables, reads in data from
    the input file, and checks the run parameters for consistency.
    """
    # Initialise the program variables
    iteration_variables.initialise_iteration_variables(data)

    # Creating and open the files MFile and OUTFile
    process_output.OutputFileManager.open_files(data.globals.output_prefix)

    # Input any desired new initial values
    inputs = parse_input_file(data)

    # Set active constraints
    set_active_constraints(data)

    # set the device type (icase)
    set_device_type(data)

    # Initialise the Stellarator
    st_init(data)

    # Check input data for errors/ambiguities
    check_process(inputs, data)

    run_summary(data)


def get_git_summary() -> tuple[str, str]:
    directory = Path(process.__file__).parent
    try:
        git_branch = (
            subprocess  # noqa: S602
            .run(
                "git rev-parse --abbrev-ref HEAD",  # noqa: S607
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

        git_tag = (
            subprocess  # noqa: S602
            .run(
                "git describe --tags",  # noqa: S607
                shell=True,
                capture_output=True,
                cwd=directory,
                check=True,
            )
            .stdout.decode()
            .strip()
        )

    except (subprocess.CalledProcessError, AttributeError):
        return "", ""
    else:
        return git_branch, git_tag


def run_summary(data: DataStructure):
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

        fileprefix = data.globals.fileprefix
        process_output.ocmmnt(
            outfile,
            f"Input : {fileprefix}",
        )
        runtitle = data.globals.runtitle
        process_output.ocmmnt(
            outfile,
            f"Run title : {runtitle}",
        )

        process_output.ocmmnt(
            outfile,
            f"Run type : Reactor concept design: {data.globals.icase}, (c) UK Atomic Energy Authority",
        )

        process_output.oblnkl(outfile)
        process_output.ostars(outfile, 110)
        process_output.oblnkl(outfile)

        process_output.ocmmnt(outfile, f"Equality constraints : {data.numerics.neqns}")
        process_output.ocmmnt(
            outfile,
            f"Inequality constraints : {data.numerics.nineqns}",
        )
        process_output.ocmmnt(
            outfile,
            f"Total constraints : {data.numerics.nineqns + data.numerics.neqns}",
        )
        process_output.ocmmnt(outfile, f"Iteration variables : {data.numerics.nvar}")
        # If optimising, write objective function and convergence parameter
        if data.numerics.ioptimz == PROCESSRunMode.OPTIMISATION:
            process_output.ocmmnt(
                outfile,
                f"Max iterations : {data.globals.maxcal}",
            )

            if data.numerics.minmax > 0:
                minmax_string = "  -- minimise "
                minmax_sign = "+"
            else:
                minmax_string = "  -- maximise "
                minmax_sign = "-"

            fom_string = FiguresOfMerit(abs(data.numerics.minmax)).description
            process_output.ocmmnt(
                outfile,
                f"Figure of merit : {minmax_sign}{abs(data.numerics.minmax)}{minmax_string}{fom_string}",
            )
            process_output.ocmmnt(
                outfile,
                f"Convergence parameter : {data.numerics.epsvmc}",
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
        mfile, "Optimisation switch", "(ioptimz)", data.numerics.ioptimz
    )
    # If optimising, write figure of merit switch
    if data.numerics.ioptimz == PROCESSRunMode.OPTIMISATION:
        process_output.ovarin(
            mfile, "Figure of merit switch", "(minmax)", data.numerics.minmax
        )


def init_all_module_vars():
    """Initialise all module variables
    This is vital to ensure a 'clean' state of Process before a new run starts,
    otherwise components of the previous run's state can persist into the new
    run. This matters ever since Process is used as a shared library, rather
    than a 'run-once' executable.
    """
    logging_model_handler.clear_logs()
    constants.init_constants()


def check_process(inputs, data):  # noqa: ARG001
    """Routine to reset specific variables if certain options are
    being used

    This routine performs a sanity check of the input variables
    and ensures other dependent variables are given suitable values.
    """
    # Check that there are sufficient iteration variables
    if data.numerics.nvar < data.numerics.neqns:
        raise ProcessValidationError(
            "Insufficient iteration variables to solve the problem! NVAR < NEQNS",
            nvar=data.numerics.nvar,
            neqns=data.numerics.neqns,
        )

    # Check that sufficient elements of ixc and icc have been specified
    if (data.numerics.ixc[: data.numerics.nvar] == 0).any():
        raise ProcessValidationError(
            "The number of iteration variables specified is smaller than the number stated in ixc",
            nvar=data.numerics.nvar,
        )

    # Check that dr_tf_wp_with_insulation (ixc = 140) and dr_tf_inboard (ixc = 13) are not being used simultaneously as iteration variables
    if (data.numerics.ixc[: data.numerics.nvar] == 13).any() and (
        data.numerics.ixc[: data.numerics.nvar] == 140
    ).any():
        raise ProcessValidationError(
            "Iteration variables 13 and 140 cannot be used simultaneously",
        )

    # Can't use c_tf_turn as interation var, constraint or input if i_tf_turns_integer == 1
    if (
        data.numerics.ixc[: data.numerics.nvar] == 60
    ).any() and data.tfcoil.i_tf_turns_integer == 1:
        raise ProcessValidationError(
            "Iteration variable 60 (TF current per turn, c_tf_turn) cannot be used with the TF coil integer turn model (i_tf_turns_integer == 1) as it is a calculated output instead for this model. However, the maximum current per turn can be constrained with constraint 77."
        )

    # Can't have icc 77 and ixc 60 at the same time
    if (data.numerics.ixc[: data.numerics.nvar] == 60).any() and (
        data.numerics.icc[: data.numerics.nvar] == 77
    ).any():
        raise ProcessValidationError(
            "Cannot use iteration variable 60 (TF coil current per turn, c_tf_turn) and constraint 77 (maximum TF current per turn) simultaneously."
        )

    if (data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 0).any():
        raise ProcessValidationError(
            "The number of constraints specified is smaller than the number stated in neqns+nineqns",
            neqns=data.numerics.neqns,
            nineqns=data.numerics.nineqns,
        )

    # Deprecate constraints
    for depcrecated_constraint in [3, 4, 10, 74, 42]:
        if (
            data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns]
            == depcrecated_constraint
        ).any():
            raise ProcessValidationError(
                "Constraint equation is no longer available", icc=depcrecated_constraint
            )

    # MDK Report error if constraint 63 is used with old vacuum model
    if (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 63
    ).any() and data.vacuum.i_vacuum_pumping != "simple":
        raise ProcessValidationError(
            "Constraint 63 is requested without the correct vacuum model (simple)"
        )

    #  Fuel ion fractions must add up to 1.0
    if (
        abs(
            1.0
            - data.physics.f_plasma_fuel_deuterium
            - data.physics.f_plasma_fuel_tritium
            - data.physics.f_plasma_fuel_helium3
        )
        > 1e-6
    ):
        raise ProcessValidationError(
            "Fuel ion fractions do not sum to 1.0",
            f_plasma_fuel_deuterium=data.physics.f_plasma_fuel_deuterium,
            f_plasma_fuel_tritium=data.physics.f_plasma_fuel_tritium,
            f_plasma_fuel_helium3=data.physics.f_plasma_fuel_helium3,
        )

    if data.physics.f_plasma_fuel_tritium < 1.0e-3:  # tritium fraction is negligible
        data.buildings.triv = 0.0
        data.heat_transport.p_tritium_plant_electric_mw = 0.0

    if data.impurity_radiation.f_nd_impurity_electrons[1] != 0.1:  # noqa: RUF069
        raise ProcessValidationError(
            "The thermal alpha/electron density ratio should be controlled using f_nd_alpha_electron (itv 109) and not f_nd_impurity_electrons(2)."
            "f_nd_impurity_electrons(2) should be removed from the input file, or set to the default value 0.1D0."
        )

    # Impurity fractions
    for imp in range(N_IMPURITIES):
        data.impurity_radiation.f_nd_impurity_electron_array[imp] = (
            data.impurity_radiation.f_nd_impurity_electrons[imp]
        )

    # Stop the run if j_tf_coil_full_area is used as an optimisation variable
    # As the current density is now calculated from b_plasma_toroidal_on_axis without constraint 10

    if (data.numerics.ixc[: data.numerics.nvar] == 12).any():
        raise ProcessValidationError(
            "The 1/R toroidal B field dependency constraint is being depreciated"
        )

    # Plasma profile consistency checks
    if data.ife.ife != 1 and data.physics.i_plasma_pedestal == 1:
        # Temperature checks
        if (
            data.physics.temp_plasma_pedestal_kev
            < data.physics.temp_plasma_separatrix_kev
        ):
            raise ProcessValidationError(
                "Pedestal temperature is lower than separatrix temperature",
                temp_plasma_pedestal_kev=data.physics.temp_plasma_pedestal_kev,
                temp_plasma_separatrix_kev=data.physics.temp_plasma_separatrix_kev,
            )

        if (abs(data.physics.radius_plasma_pedestal_temp_norm - 1.0) <= 1e-7) and (
            (
                data.physics.temp_plasma_pedestal_kev
                - data.physics.temp_plasma_separatrix_kev
            )
            >= 1e-7
        ):
            warn(
                f"Temperature pedestal is at plasma edge, but temp_plasma_pedestal_kev "
                f"({data.physics.temp_plasma_pedestal_kev}) differs from temp_plasma_separatrix_kev "
                f"({data.physics.temp_plasma_separatrix_kev})",
                stacklevel=2,
            )

        # Core temperature should always be calculated (later) as being
        # higher than the pedestal temperature, if and only if the
        # volume-averaged temperature never drops below the pedestal
        # temperature. Prevent this by adjusting te, and its lower bound
        # (which will only have an effect if this is an optimisation run)
        if (
            data.physics.temp_plasma_electron_vol_avg_kev
            <= data.physics.temp_plasma_pedestal_kev
        ):
            warn(
                f"Volume-averaged temperature ({data.physics.te}) has been "
                f"forced to exceed input pedestal height ({data.physics.temp_plasma_pedestal_kev}). "
                "Changing to te = temp_plasma_pedestal_kev*1.001",
                stacklevel=2,
            )
            data.physics.temp_plasma_electron_vol_avg_kev = (
                data.physics.temp_plasma_pedestal_kev * 1.001
            )

        if (
            data.numerics.ioptimz == PROCESSRunMode.OPTIMISATION
            and (data.numerics.ixc[: data.numerics.nvar] == 4).any()
            and data.numerics.boundl[3] < data.physics.temp_plasma_pedestal_kev * 1.001
        ):
            warn(
                "Lower limit of volume averaged electron temperature (temp_plasma_electron_vol_avg_kev) has been raised to ensure temp_plasma_electron_vol_avg_kev > temp_plasma_pedestal_kev",
                stacklevel=2,
            )
            data.numerics.boundl[3] = data.physics.temp_plasma_pedestal_kev * 1.001
            data.numerics.boundu[3] = max(
                data.numerics.boundu[3], data.numerics.boundl[3]
            )

        # Density checks
        # Issue #589: Pedestal density is lower than separatrix density
        pedestal_type = DensityProfilePedestalType(
            data.physics.i_nd_plasma_pedestal_separatrix
        )
        if (
            pedestal_type == DensityProfilePedestalType.USER_INPUT
            and data.physics.nd_plasma_pedestal_electron
            < data.physics.nd_plasma_separatrix_electron
        ) or (
            pedestal_type == DensityProfilePedestalType.GREENWALD_FRACTION
            and data.physics.f_nd_plasma_pedestal_greenwald
            < data.physics.f_nd_plasma_separatrix_greenwald
        ):
            raise ProcessValidationError(
                "Density pedestal is lower than separatrix density",
                **(
                    {
                        "nd_plasma_pedestal_electron": data.physics.nd_plasma_pedestal_electron,
                        "nd_plasma_separatrix_electron": data.physics.nd_plasma_separatrix_electron,
                    }
                    if pedestal_type == DensityProfilePedestalType.USER_INPUT
                    else {
                        "f_nd_plasma_pedestal_greenwald": data.physics.f_nd_plasma_pedestal_greenwald,
                        "f_nd_plasma_separatrix_greenwald": data.physics.f_nd_plasma_separatrix_greenwald,
                    }
                ),
            )

        if (
            abs(data.physics.radius_plasma_pedestal_density_norm - 1.0) <= 1e-7
            and (
                data.physics.nd_plasma_pedestal_electron
                - data.physics.nd_plasma_separatrix_electron
            )
            >= 1e-7
        ):
            warn(
                "Density pedestal is at plasma edge "
                f"({data.physics.radius_plasma_pedestal_density_norm = }), but nd_plasma_pedestal_electron "
                f"({data.physics.nd_plasma_pedestal_electron}) differs from "
                f"nd_plasma_separatrix_electron ({data.physics.nd_plasma_separatrix_electron})",
                stacklevel=2,
            )

        # Issue #862 : Variable nd_plasma_electron_on_axis/nd_plasma_pedestal_electron ratio without constraint eq 81 (nd_plasma_electron_on_axis>nd_plasma_pedestal_electron)
        #  -> Potential hollowed density profile
        if (
            data.numerics.ioptimz == PROCESSRunMode.OPTIMISATION
            and not (
                data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 81
            ).any()
        ):
            if (data.numerics.ixc[: data.numerics.nvar] == 145).any():
                warn(
                    "nd_plasma_pedestal_electron set with f_nd_plasma_pedestal_greenwald without constraint eq 81 (nd_plasma_pedestal_electron<nd_plasma_electron_on_axis)",
                    stacklevel=2,
                )
            if (data.numerics.ixc[: data.numerics.nvar] == 6).any():
                warn(
                    "nd_plasma_electrons_vol_avg used as iteration variable without constraint 81 (nd_plasma_pedestal_electron<nd_plasma_electron_on_axis)",
                    stacklevel=2,
                )

    # Cannot use Psep/R and PsepB/qAR limits at the same time
    if (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 68
    ).any() and (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 56
    ).any():
        raise ProcessValidationError(
            "Cannot use Psep/R and PsepB/qAR constraint equations at the same time"
        )

    # if lower bound of f_nd_plasma_pedestal_greenwald < f_nd_plasma_separatrix_greenwald
    if (data.numerics.ixc[: data.numerics.nvar] == 145).any() and data.numerics.boundl[
        144
    ] < data.physics.f_nd_plasma_separatrix_greenwald:
        raise ProcessValidationError(
            "Set lower bound of iteration variable 145, f_nd_plasma_pedestal_greenwald, to be greater than f_nd_plasma_separatrix_greenwald",
            boundl_145=data.numerics.boundl[144],
            f_nd_plasma_separatrix_greenwald=data.physics.f_nd_plasma_separatrix_greenwald,
        )

    if (data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 78).any():
        # If Reinke criterion is used temp_plasma_separatrix_kev is calculated and cannot be an
        # iteration variable
        if (data.numerics.ixc[: data.numerics.nvar] == 119).any():
            raise ProcessValidationError(
                "REINKE IMPURITY MODEL: temp_plasma_separatrix_kev is calculated and cannot be an "
                "iteration variable for the Reinke model"
            )

        # If Reinke criterion is used need to enforce LH-threshold
        # using Martin scaling for consistency
        if (data.physics.i_l_h_threshold != 6) or (
            not (
                data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 15
            ).any()
            and data.physics.i_plasma_pedestal
        ):
            warn(
                "REINKE IMPURITY MODEL: The Martin LH threshold scale is not being used and is recommended for the Reinke model",
                stacklevel=2,
            )
    i_single_null = DivertorNumberModels(data.physics.i_single_null)
    if i_single_null == DivertorNumberModels.DOUBLE_NULL:
        data.divertor.n_divertors = 2
        data.build.dz_fw_plasma_gap = data.build.dz_xpoint_divertor
        data.build.dz_shld_upper = data.build.dz_shld_lower
        data.build.dz_vv_upper = data.build.dz_vv_lower
        warn("Double-null: Upper vertical build forced to match lower", stacklevel=2)
    else:  # i_single_null == DivertorNumberModels.SINGLE_NULL
        data.divertor.n_divertors = 1

    #  Tight aspect ratio options (ST)
    if data.physics.itart == 1:
        data.globals.icase = "Tight aspect ratio tokamak model"

        # Disabled Forcing that no inboard breeding blanket is used
        # Disabled i_blkt_inboard = 0

        # Check if the choice of plasma current is addapted for ST
        # 2 : Peng Ip scaling (See STAR code documentation)
        # 9 : Fiesta Ip scaling
        if data.physics.i_plasma_current not in {2, 9}:
            warn(
                "Usual current scaling for TARTs (i_plasma_current=2 or 9) is not being used",
                stacklevel=2,
            )

        # If using Peng and Strickler (1986) model (itartpf == 0)
        # Overwrite the location of the TF coils
        # 2 : PF coil on top of TF coil
        # 3 : PF coil outside of TF coil
        if data.physics.itartpf == 0:
            data.pf_coil.i_pf_location[0] = PFLocationTypes.ABOVE_TF
            data.pf_coil.i_pf_location[1] = PFLocationTypes.OUTSIDE_TF
            data.pf_coil.i_pf_location[2] = PFLocationTypes.OUTSIDE_TF

        # Water cooled copper magnets initalisation / checks
        if data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
            # Check if the initial centrepost coolant loop adapted to the magnet technology
            # Ice cannot flow so temp_cp_coolant_inlet > 273.15 K
            if data.tfcoil.temp_cp_coolant_inlet < 273.15:
                raise ProcessValidationError(
                    "Coolant temperature (temp_cp_coolant_inlet) cannot be < 0 C (273.15 K) for water cooled copper magents"
                )

            # Temperature of the TF legs cannot be cooled down
            if (
                data.tfcoil.temp_tf_legs_outboard > 0
                and data.tfcoil.temp_tf_legs_outboard < 273.15
            ):
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) cannot be < 0 C (273.15 K) for water cooled magents"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                data.numerics.ixc[: data.numerics.nvar] == 20
            ).any() and data.numerics.boundu[19] < 273.15:
                raise ProcessValidationError(
                    "Too low CP conductor temperature (temp_cp_average). Lower limit for copper > 273.15 K"
                )

        # Call a lvl 3 error if superconductor magnets are used
        elif data.tfcoil.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
            warn(
                "Joints res not cal. for SC (itart = 1) TF (data.tfcoil.i_tf_sup = 1)",
                stacklevel=2,
            )

        # Aluminium magnets initalisation / checks
        # Initialize the CP conductor temperature to cryogenic temperature for cryo-al magnets (20 K)
        elif data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            # Call a lvl 3 error if the inlet coolant temperature is too large
            # Motivation : ill-defined aluminium resistivity fit for T > 40-50 K
            if data.tfcoil.temp_cp_coolant_inlet > 40.0:
                raise ProcessValidationError(
                    "Coolant temperature (temp_cp_coolant_inlet) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if the leg average temperature is low enough for the resisitivity fit
            if data.tfcoil.temp_tf_legs_outboard > 50.0:
                raise ProcessValidationError(
                    "TF legs conductor temperature (temp_tf_legs_outboard) should be < 40 K for the cryo-al resistivity to be defined"
                )

            # Check if conductor upper limit is properly set to 50 K or below
            if (
                data.numerics.ixc[: data.numerics.nvar] == 20
            ).any() and data.numerics.boundu[19] > 50.0:
                raise ProcessValidationError(
                    "Too large CP conductor temperature (temp_cp_average). Upper limit for cryo-al < 50 K"
                )

            # Otherwise intitialise the average conductor temperature at
            data.tfcoil.temp_cp_average = data.tfcoil.temp_cp_coolant_inlet

        # Check if the boostrap current selection is addapted to ST
        if data.physics.i_bootstrap_current == 1:
            raise ProcessValidationError(
                "Invalid boostrap current law for ST, do not use i_bootstrap_current = 1"
            )

        # Check if a single null divertor is used in double null machine
        if i_single_null == DivertorNumberModels.DOUBLE_NULL and (
            data.physics.f_p_div_lower in {1.0, 0.0}
        ):
            warn("Operating with a single null in a double null machine", stacklevel=2)

        # Set the TF coil shape to picture frame (if default value)
        if data.tfcoil.i_tf_shape == TFCoilShapeModel.DEFAULT:
            data.tfcoil.i_tf_shape = TFCoilShapeModel.PICTURE_FRAME

        # Warning stating that the CP fast neutron fluence calculation
        # is not addapted for cryoaluminium calculations yet
        if (
            data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM
            and (
                data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 85
            ).any()
            and data.physics.itart == 1
        ):
            raise ProcessValidationError(
                "Al TF coil fluence not calculated properly for Al CP, do not use constraint 85"
            )

        # Setting the CP joints default options :
        #  0 : No joints for superconducting magents (data.tfcoil.i_tf_sup = 1)
        #  1 : Sliding joints for resistive magnets (data.tfcoil.i_tf_sup = 0, 2)
        if data.tfcoil.i_cp_joints == -1:
            if data.tfcoil.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
                data.tfcoil.i_cp_joints = 0
            else:
                data.tfcoil.i_cp_joints = 1

        # Checking the CP TF top radius
        if (
            abs(data.build.r_cp_top) > 0
            or (data.numerics.ixc[: data.numerics.nvar] == 174).any()
        ) and data.build.i_r_cp_top != 1:
            raise ProcessValidationError(
                "To set the TF CP top value, you must use i_r_cp_top = 1"
            )

    # Conventionnal aspect ratios specific
    else:
        if data.physics.i_plasma_current in {2, 9}:
            raise ProcessValidationError(
                "i_plasma_current=2,9 is not a valid option for a non-TART device"
            )

        # Set the TF coil shape to PROCESS D-shape (if default value)
        if data.tfcoil.i_tf_shape == TFCoilShapeModel.DEFAULT:
            data.tfcoil.i_tf_shape = TFCoilShapeModel.D_SHAPE

        # Check PF coil configurations
        j = 0
        k = 0
        for i in range(data.pf_coil.n_pf_coil_groups):
            if (
                PFLocationTypes(data.pf_coil.i_pf_location[i])
                != PFLocationTypes.ABOVE_TF
                and data.pf_coil.n_pf_coils_in_group[i] != 2
            ):
                raise ProcessValidationError(
                    "n_pf_coils_in_group(i) .ne. 2 is not a valid option except for (i_pf_location = 2)"
                )

            if data.pf_coil.i_pf_location[i] == PFLocationTypes.ABOVE_TF:
                j += 1
                k += data.pf_coil.n_pf_coils_in_group[i]

        if k == 1:
            raise ProcessValidationError(
                "Only 1 divertor coil (i_pf_location = 2) is not a valid configuration"
            )
        if k > 2:
            raise ProcessValidationError(
                "More than 2 divertor coils (i_pf_location = 2) is not a valid configuration"
            )
        if i_single_null == DivertorNumberModels.SINGLE_NULL and j < 2:
            raise ProcessValidationError(
                "If i_single_null=1, use 2 individual divertor coils (i_pf_location = 2, 2; n_pf_coils_in_group = 1, 1)"
            )

        # Constraint 10 is dedicated to ST designs with demountable joints
        if (
            data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 10
        ).any():
            raise ProcessValidationError(
                "Constraint equation 10 (CP lifetime) to used with ST desing (itart=1)"
            )

    #  Pulsed power plant model
    if data.pulse.i_pulsed_plant == 1:
        data.globals.icase = "Pulsed tokamak model"
    else:
        data.buildings.esbldgm3 = 0.0

    # TF coil
    # -------
    # TF stress model not defined of r_tf_inboard = 0
    # Unless i_tf_stress_model == 2
    # -> If dr_bore + dr_cs_tf_gap + dr_cs = 0 and fixed and stress constraint is used
    #    Generate a lvl 3 error proposing not to use any stress constraints
    if (
        (
            not (
                (data.numerics.ixc[: data.numerics.nvar] == 16).any()
                or (data.numerics.ixc[: data.numerics.nvar] == 29).any()
                or (data.numerics.ixc[: data.numerics.nvar] == 42).any()
            )
        )  # No dr_bore,dr_cs_tf_gap, dr_cs iteration
        and (
            abs(
                data.build.dr_bore
                + data.build.dr_cs_tf_gap
                + data.build.dr_cs
                + data.build.dr_cs_precomp
            )
            <= 0
        )  # dr_bore + dr_cs_tf_gap + dr_cs = 0
        and (
            (
                data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 31
            ).any()
            or (
                data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 32
            ).any()
        )  # Stress constraints (31 or 32) is used
        and (
            data.tfcoil.i_tf_stress_model != 2
        )  # TF stress model can't handle no dr_bore
    ):
        raise ProcessValidationError(
            "Invalid stress model if dr_bore + dr_cs_tf_gap + dr_cs = 0. Don't use constraint 31"
        )

    # Make sure that plane stress model is not used for resistive magnets
    if (
        data.tfcoil.i_tf_stress_model == 1
        and data.tfcoil.i_tf_sup != TFConductorModel.SUPERCONDUCTING
    ):
        raise ProcessValidationError(
            "Use generalized plane strain for resistive magnets (i_tf_stress_model = 0 or 2)"
        )

    # bucking cylinder default option setting
    # - bucking (casing) for SC i_tf_bucking ( i_tf_bucking = 1 )
    # - No bucking for copper magnets ( i_tf_bucking = 0 )
    # - Bucking for aluminium magnets ( i_tf_bucking = 1 )
    if data.tfcoil.i_tf_bucking == -1:
        if data.tfcoil.i_tf_sup == TFConductorModel.WATER_COOLED_COPPER:
            data.tfcoil.i_tf_bucking = 0
        else:
            data.tfcoil.i_tf_bucking = 1

    # Ensure that the TF isnt placed against the
    # CS which is now outside it
    if (
        data.tfcoil.i_tf_bucking >= 2
        and data.build.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS
    ):
        raise ProcessValidationError(
            "Cannot have i_tf_bucking >= 2 when i_tf_inside_cs = 1"
        )

    # Ensure that no pre-compression structure
    # is used for bucked and wedged design
    if data.tfcoil.i_tf_bucking >= 2 and data.build.i_cs_precomp == 1:
        raise ProcessValidationError(
            "No CS precompression structure for bucked and wedged, use i_cs_precomp = 0"
        )

    # Number of stress calculation layers
    # +1 to add in the inboard TF coil case on the plasma side, per Issue #1509
    data.tfcoil.n_tf_stress_layers = (
        data.tfcoil.i_tf_bucking + data.tfcoil.n_tf_graded_layers + 1
    )

    # If TFC sidewall has not been set by user
    if data.tfcoil.dx_tf_side_case_min < 0.1e-10:
        data.tfcoil.tfc_sidewall_is_fraction = True

    # If inboard TF coil case plasma side thickness has not been set by user
    if data.tfcoil.dr_tf_plasma_case < 0.1e-10:
        data.tfcoil.i_f_dr_tf_plasma_case = True

    # Setting the default cryo-plants efficiencies
    if abs(data.tfcoil.eff_tf_cryo + 1) < 1e-6:
        # The ITER cyoplant efficiency is used for SC
        if data.tfcoil.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
            data.tfcoil.eff_tf_cryo = 0.13

        # Strawbrige plot extrapolation is used for Cryo-Al
        elif data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            data.tfcoil.eff_tf_cryo = 0.40

    # Cryo-plane efficiency must be in [0-1.0]
    elif data.tfcoil.eff_tf_cryo > 1.0 or data.tfcoil.eff_tf_cryo < 0.0:
        raise ProcessValidationError(
            "TF cryo-plant efficiency `eff_tf_cryo` must be within [0-1]"
        )

    # Integer turns option not yet available for REBCO taped turns

    if (
        data.tfcoil.i_tf_sc_mat == SuperconductorModel.CROCO_REBCO
        and data.tfcoil.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Integer turns (i_tf_turns_integer = 1) not supported for REBCO (i_tf_sc_mat = 6)"
        )

    # Setting up insulation layer young modulae default values [Pa]

    if data.tfcoil.eyoung_ins <= 1.0e8:
        # Copper magnets, no insulation material defined
        # But use the ITER design by default
        if data.tfcoil.i_tf_sup in {
            TFConductorModel.WATER_COOLED_COPPER,
            TFConductorModel.SUPERCONDUCTING,
        }:
            # SC magnets
            # Value from DDD11-2 v2 2 (2009)
            data.tfcoil.eyoung_ins = 20.0e9

        # Cryo-aluminum magnets (Kapton polymer)
        elif data.tfcoil.i_tf_sup == TFConductorModel.HELIUM_COOLED_ALUMINIUM:
            data.tfcoil.eyoung_ins = 2.5e9

    # Setting the default WP geometry
    i_tf_wp_geom = SuperconductingTFWPShapeType(data.tfcoil.i_tf_wp_geom)
    if i_tf_wp_geom == SuperconductingTFWPShapeType.UNSET:
        if data.tfcoil.i_tf_turns_integer == 0:
            data.tfcoil.i_tf_wp_geom = SuperconductingTFWPShapeType.DOUBLE_RECTANGULAR
        if data.tfcoil.i_tf_turns_integer == 1:
            data.tfcoil.i_tf_wp_geom = SuperconductingTFWPShapeType.RECTANGULAR

    # Setting the TF coil conductor elastic properties

    if data.tfcoil.i_tf_cond_eyoung_axial == 0:
        # Conductor stiffness is not considered
        data.tfcoil.eyoung_cond_axial = 0
        data.tfcoil.eyoung_cond_trans = 0
    elif data.tfcoil.i_tf_cond_eyoung_axial == 2:
        # Select sensible defaults from the literature
        if (
            SuperconductorModel(data.tfcoil.i_tf_sc_mat).material
            == SuperconductorMaterial.NB3SN
        ):
            # Nb3Sn: Nyilas, A et. al, Superconductor Science and Technology 16, no. 9 (2003): 1036-42. https://doi.org/10.1088/0953-2048/16/9/313.
            data.tfcoil.eyoung_cond_axial = 32e9
        elif (
            SuperconductorModel(data.tfcoil.i_tf_sc_mat).material
            == SuperconductorMaterial.BI2212
        ):
            # Bi-2212: Brown, M. et al, IOP Conference Series: Materials Science and Engineering 279 (2017): 012022. https://doi.org/10.1088/1757-899X/279/1/012022.
            data.tfcoil.eyoung_cond_axial = 80e9
        elif (
            SuperconductorModel(data.tfcoil.i_tf_sc_mat).material
            == SuperconductorMaterial.NBTI
        ):
            # NbTi: Vedrine, P. et. al, IEEE Transactions on Applied Superconductivity 9, no. 2 (1999): 236-39. https://doi.org/10.1109/77.783280.
            data.tfcoil.eyoung_cond_axial = 6.8e9
        elif (
            SuperconductorModel(data.tfcoil.i_tf_sc_mat).material
            == SuperconductorMaterial.REBCO
        ):
            # REBCO: Fujishiro, H. et. al, Physica C: Superconductivity, 426-431 (2005): 699-704. https://doi.org/10.1016/j.physc.2005.01.045.
            data.tfcoil.eyoung_cond_axial = 145e9

        if data.tfcoil.i_tf_cond_eyoung_trans == 0:
            # Transverse stiffness is not considered
            data.tfcoil.eyoung_cond_trans = 0
        else:
            # Transverse stiffness is significant
            data.tfcoil.eyoung_cond_trans = data.tfcoil.eyoung_cond_axial

    # Check if the WP/conductor radial thickness (dr_tf_wp_with_insulation) is large enough
    # To contains the insulation, cooling and the support structure
    # Rem : Only verified if the WP thickness is used
    if (data.numerics.ixc[: data.numerics.nvar] == 140).any():
        # Minimal WP thickness
        if data.tfcoil.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
            dr_tf_wp_min = 2.0 * (
                data.tfcoil.dx_tf_wp_insulation
                + data.tfcoil.dx_tf_wp_insertion_gap
                + data.tfcoil.dx_tf_turn_insulation
                + data.tfcoil.dia_tf_turn_coolant_channel
            )

            # Steel conduit thickness (can be an iteration variable)
            if (data.numerics.ixc[: data.numerics.nvar] == 58).any():
                dr_tf_wp_min += 2.0 * data.numerics.boundl[57]
            else:
                dr_tf_wp_min += 2.0 * data.tfcoil.dx_tf_turn_steel

        # Minimal conductor layer thickness
        elif data.tfcoil.i_tf_sup in {
            TFConductorModel.WATER_COOLED_COPPER,
            TFConductorModel.HELIUM_COOLED_ALUMINIUM,
        }:
            dr_tf_wp_min = (
                2.0
                * (data.tfcoil.dx_tf_turn_insulation + data.tfcoil.dx_tf_wp_insulation)
                + 4.0 * data.tfcoil.radius_cp_coolant_channel
            )

        if data.numerics.boundl[139] < dr_tf_wp_min:
            raise ProcessValidationError(
                "The TF coil WP thickness (dr_tf_wp_with_insulation) must be at least",
                dr_tf_wp_min=dr_tf_wp_min,
            )

    # Setting i_dx_tf_turn_general_input to true if dx_tf_turn_general is an input
    data.tfcoil.i_dx_tf_turn_general_input = abs(data.tfcoil.dx_tf_turn_general) > 0

    # Impossible to set the turn size of integer turn option
    if data.tfcoil.i_dx_tf_turn_general_input and data.tfcoil.i_tf_turns_integer == 1:
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    if (
        data.tfcoil.i_tf_wp_geom != SuperconductingTFWPShapeType.RECTANGULAR
        and data.tfcoil.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Can only have i_tf_turns_integer = 1 with i_tf_wp_geom = 0"
        )

    if data.physics.i_bootstrap_current == 5 and data.physics.i_diamagnetic_current != 0:
        raise ProcessValidationError(
            "i_diamagnetic_current = 0 should be used with the Sakai plasma current scaling"
        )

    # Setting i_dx_tf_turn_cable_space_general_input to true if dx_tf_turn_cable_space_general is an input
    data.tfcoil.i_dx_tf_turn_cable_space_general_input = (
        abs(data.tfcoil.dx_tf_turn_cable_space_general) > 0
    )

    # Impossible to set the cable size of integer turn option
    if (
        data.tfcoil.i_dx_tf_turn_cable_space_general_input
        and data.tfcoil.i_tf_turns_integer == 1
    ):
        raise ProcessValidationError(
            "Impossible to set the TF turn/cable size with the integer turn option (i_tf_turns_integer: 1)"
        )

    # Impossible to set both the TF coil turn and the cable dimension
    if (
        data.tfcoil.i_dx_tf_turn_general_input
        and data.tfcoil.i_dx_tf_turn_cable_space_general_input
    ):
        raise ProcessValidationError(
            "Impossible to set the TF coil turn and cable size simultaneously"
        )

    # Checking the SC temperature for LTS
    if (
        SuperconductorModel(data.tfcoil.i_tf_sc_mat).sc_type
        == SuperconductorType.LOW_TEMPERATURE
        and data.tfcoil.tftmp > 10.0
    ):
        raise ProcessValidationError(
            "The LTS conductor temperature (tftmp) has to be lower than 10"
        )

    # PF coil resistivity is zero if superconducting
    if data.pf_coil.i_pf_conductor == PFConductorModel.SUPERCONDUCTING:
        data.pf_coil.rho_pf_coil = 0.0

    # If there is no NBI, then hot beam density should be zero
    if data.current_drive.i_hcd_calculations == 1:
        if data.current_drive.i_hcd_primary not in {5, 8}:
            data.physics.f_nd_beam_electron = 0.0
    else:
        data.physics.f_nd_beam_electron = 0.0

    # Set inboard blanket thickness to zero if no inboard blanket switch
    # used (Issue #732)
    if data.build.i_blkt_inboard == InboardBlanketConfiguration.NO_INBOARD_BLANKET:
        data.build.dr_blkt_inboard = 0.0

    # Ensure that blanket material fractions allow non-zero space for steel
    # CCFE HCPB Model

    if data.stellarator.istell == 0 and (
        data.fwbs.i_blanket_type == BlktModelTypes.CCFE_HCPB
    ):
        fsum = data.fwbs.breeder_multiplier + data.fwbs.vfcblkt + data.fwbs.vfpblkt
        if fsum >= 1.0:
            raise ProcessValidationError(
                "Blanket material fractions do not sum to 1.0",
                i_blanket_type=data.fwbs.i_blanket_type,
                breeder_multiplier=data.fwbs.breeder_multiplier,
                vfcblkt=data.fwbs.vfcblkt,
                vfpblkt=data.fwbs.vfpblkt,
                fsum=fsum,
            )

    # Check that the temperature margins are not overdetermined
    if data.tfcoil.tmargmin > 0.0001:
        # This limit has been input and will be applied to both TFC and CS
        if data.tfcoil.temp_tf_superconductor_margin_min > 0.0001:
            warn(
                "temp_tf_superconductor_margin_min and tmargmin should not both be specified in IN.DAT "
                "temp_tf_superconductor_margin_min has been ignored",
                stacklevel=2,
            )
        if data.tfcoil.temp_cs_superconductor_margin_min > 0.0001:
            warn(
                "temp_cs_superconductor_margin_min and tmargmin should not both be specified in IN.DAT "
                "temp_cs_superconductor_margin_min has been ignored",
                stacklevel=2,
            )

        data.tfcoil.temp_tf_superconductor_margin_min = data.tfcoil.tmargmin
        data.tfcoil.temp_cs_superconductor_margin_min = data.tfcoil.tmargmin

    if data.physics.tauee_in > 1e-10 and data.physics.i_confinement_time != 48:
        # Report error if confinement time is in the input
        # but the scaling to use it is not selected.
        warn("tauee_in is for use with i_confinement_time=48 only", stacklevel=2)

    if data.physics.aspect > 1.7 and data.physics.i_confinement_time == 46:
        # NSTX scaling is for A<1.7
        warn("NSTX scaling is for A<1.7", stacklevel=2)

    if data.physics.i_plasma_current == 2 and data.physics.i_confinement_time == 42:
        raise ProcessValidationError(
            "Lang 2012 confinement scaling cannot be used for i_plasma_current=2 due to wrong q"
        )

    # Cannot use temperature margin constraint with REBCO TF coils
    if (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 36
    ).any() and (
        SuperconductorModel(data.tfcoil.i_tf_sc_mat).sc_type
        == SuperconductorMaterial.REBCO
    ):
        raise ProcessValidationError(
            "turn off TF temperature margin constraint icc = 36 when using REBCO"
        )

    # Cannot use temperature margin constraint with REBCO CS coils
    if (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 60
    ).any() and data.pf_coil.i_cs_superconductor == 8:
        raise ProcessValidationError(
            "turn off CS temperature margin constraint icc = 60 when using REBCO"
        )

    # Cold end of the cryocooler should be colder than the TF
    if data.tfcoil.temp_tf_cryo > data.tfcoil.tftmp:
        raise ProcessValidationError("temp_tf_cryo should be lower than tftmp")

    # Cannot use TF coil strain limit if i_str_wp is off:
    if (
        data.numerics.icc[: data.numerics.neqns + data.numerics.nineqns] == 88
    ).any() and data.tfcoil.i_str_wp == 0:
        raise ProcessValidationError("Can't use constraint 88 if i_strain_tf == 0")


def set_active_constraints(data: DataStructure):
    """Set constraints provided in the input file as 'active'"""
    num_constraints = 0
    for i in range(ConstraintManager.num_constraints()):
        if data.numerics.icc[i] != 0:
            data.numerics.active_constraints[data.numerics.icc[i] - 1] = True
            num_constraints += 1

    if data.numerics.neqns < 0:
        # The value of neqns has not been set in the input file.  Default = 0.
        data.numerics.neqns = num_constraints - data.numerics.nineqns
    else:
        data.numerics.nineqns = num_constraints - data.numerics.neqns


def set_device_type(data: DataStructure):
    if data.ife.ife == 1:
        data.globals.icase = "Inertial Fusion model"
    elif data.stellarator.istell != 0:
        data.globals.icase = "Stellarator model"
