import logging
import time
from dataclasses import astuple, dataclass
from enum import Enum

import numpy as np
from tabulate import tabulate

import process.core.optimisation.constraints as constraints
import process.process_output as process_output
from process.core import constants
from process.core.caller import write_output_files
from process.data_structure import (
    build_variables,
    constraint_variables,
    cost_variables,
    cs_fatigue_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    impurity_radiation_module,
    numerics,
    pfcoil_variables,
    physics_variables,
    rebco_variables,
    scan_variables,
    tfcoil_variables,
)
from process.exceptions import ProcessValueError
from process.log import logging_model_handler, show_errors
from process.solver_handler import SolverHandler

logger = logging.getLogger(__name__)


@dataclass
class ScanVariable:
    variable_name: str
    variable_description: str
    variable_num: int

    def __iter__(self):
        return iter(astuple(self)[:2])


class ScanVariables(Enum):
    @classmethod
    def _missing_(cls, var):
        if isinstance(var, int):
            for sv in cls:
                if sv.value.variable_num == var:
                    return sv
        return super()._missing_(var)

    aspect = ScanVariable("aspect", "Aspect_ratio", 1)
    pflux_div_heat_load_max_mw = ScanVariable(
        "pflux_div_heat_load_max_mw", "Div_heat_limit_(MW/m2)", 2
    )
    p_plant_electric_net_required_mw = ScanVariable(
        "p_plant_electric_net_required_mw", "Net_electric_power_(MW)", 3
    )
    hfact = ScanVariable("hfact", "Confinement_H_factor", 4)
    oacdcp = ScanVariable("oacdcp", "TF_inboard_leg_J_(MA/m2)", 5)
    pflux_fw_neutron_max_mw = ScanVariable(
        "pflux_fw_neutron_max_mw", "Allow._wall_load_(MW/m2)", 6
    )
    beamfus0 = ScanVariable("beamfus0", "Beam_bkgrd_multiplier", 7)
    temp_plasma_electron_vol_avg_kev = ScanVariable(
        "temp_plasma_electron_vol_avg_kev", "Electron_temperature_keV", 9
    )
    boundu15 = ScanVariable("boundu(15)", "Volt-second_upper_bound", 10)
    beta_norm_max = ScanVariable("beta_norm_max", "Beta_coefficient", 11)
    f_c_plasma_bootstrap_max = ScanVariable(
        "f_c_plasma_bootstrap_max", "Bootstrap_fraction", 12
    )
    boundu10 = ScanVariable("boundu(10)", "H_factor_upper_bound", 13)
    fiooic = ScanVariable("fiooic", "TFC_Iop_/_Icrit_f-value", 14)
    rmajor = ScanVariable("rmajor", "Plasma_major_radius_(m)", 16)
    b_tf_inboard_max = ScanVariable("b_tf_inboard_max", "Max_toroidal_field_(T)", 17)
    eta_cd_norm_hcd_primary_max = ScanVariable(
        "eta_cd_norm_hcd_primary_max", "Maximum_CD_gamma", 18
    )
    boundl16 = ScanVariable("boundl(16)", "CS_thickness_lower_bound", 19)
    t_burn_min = ScanVariable("t_burn_min", "Minimum_burn_time_(s)", 20)
    f_t_plant_available = ScanVariable(
        "f_t_plant_available", "Plant_availability_factor", 22
    )
    p_fusion_total_max_mw = ScanVariable(
        "p_fusion_total_max_mw", "Fusion_power_limit_(MW)", 24
    )
    kappa = ScanVariable("kappa", "Plasma_elongation", 25)
    triang = ScanVariable("triang", "Plasma_triangularity", 26)
    tbrmin = ScanVariable("tbrmin", "Min_tritium_breed._ratio", 27)
    b_plasma_toroidal_on_axis = ScanVariable(
        "b_plasma_toroidal_on_axis", "Tor._field_on_axis_(T)", 28
    )
    coreradius = ScanVariable("coreradius", "Core_radius", 29)
    f_alpha_energy_confinement_min = ScanVariable(
        "f_alpha_energy_confinement_min", "t_alpha_confinement/taueff_lower_limit", 31
    )
    epsvmc = ScanVariable("epsvmc", "VMCON error tolerance", 32)
    boundu129 = ScanVariable("boundu(129)", " Neon upper limit", 38)
    boundu131 = ScanVariable("boundu(131)", " Argon upper limit", 39)
    boundu135 = ScanVariable("boundu(135)", " Xenon upper limit", 40)
    dr_blkt_outboard = ScanVariable("dr_blkt_outboard", "Outboard blanket thick.", 41)
    f_nd_impurity_electrons9 = ScanVariable(
        "f_nd_impurity_electrons(9)", "Argon fraction", 42
    )
    sig_tf_case_max = ScanVariable(
        "sig_tf_case_max", "Allowable_stress_in_tf_coil_case_Tresca_(pa)", 44
    )
    temp_tf_superconductor_margin_min = ScanVariable(
        "temp_tf_superconductor_margin_min", "Minimum_allowable_temperature_margin", 45
    )
    boundu152 = ScanVariable(
        "boundu(152)", "Max allowable f_nd_plasma_separatrix_greenwald", 46
    )
    n_tf_wp_pancakes = ScanVariable("n_tf_wp_pancakes", "TF Coil - n_tf_wp_pancakes", 48)
    n_tf_wp_layers = ScanVariable("n_tf_wp_layers", "TF Coil - n_tf_wp_layers", 49)
    f_nd_impurity_electrons13 = ScanVariable(
        "f_nd_impurity_electrons(13)", "Xenon fraction", 50
    )
    f_p_div_lower = ScanVariable("f_p_div_lower", "lower_divertor_power_fraction", 51)
    rad_fraction_sol = ScanVariable("rad_fraction_sol", "SoL radiation fraction", 52)
    boundu157 = ScanVariable("boundu(157)", "Max allowable fvssu", 53)
    Bc2_0K = ScanVariable("Bc2(0K)", "GL_NbTi Bc2(0K)", 54)
    dr_shld_inboard = ScanVariable("dr_shld_inboard", "Inboard neutronic shield", 55)
    p_cryo_plant_electric_max_mw = ScanVariable(
        "p_cryo_plant_electric_max_mw", "max allowable p_cryo_plant_electric_mw", 56
    )
    boundl2 = ScanVariable("boundl(2)", "b_plasma_toroidal_on_axis minimum", 57)
    dr_fw_plasma_gap_inboard = ScanVariable(
        "dr_fw_plasma_gap_inboard", "Inboard FW-plasma sep gap", 58
    )
    dr_fw_plasma_gap_outboard = ScanVariable(
        "dr_fw_plasma_gap_outboard", "Outboard FW-plasma sep gap", 59
    )
    sig_tf_wp_max = ScanVariable(
        "sig_tf_wp_max", "Allowable_stress_in_tf_coil_conduit_Tresca_(pa)", 60
    )
    copperaoh_m2_max = ScanVariable(
        "copperaoh_m2_max", "Max CS coil current / copper area", 61
    )
    coheof = ScanVariable("coheof", "CS coil current density at EOF (A/m2)", 62)
    dr_cs = ScanVariable("dr_cs", "CS coil thickness (m)", 63)
    ohhghf = ScanVariable("ohhghf", "CS height (m)", 64)
    n_cycle_min = ScanVariable("n_cycle_min", "CS stress cycles min", 65)
    oh_steel_frac = ScanVariable("oh_steel_frac", "CS steel fraction", 66)
    t_crack_vertical = ScanVariable(
        "t_crack_vertical", "Initial crack vertical size (m)", 67
    )
    inlet_temp_liq = ScanVariable(
        "inlet_temp_liq", "Inlet Temperature Liquid Metal Breeder/Coolant (K)", 68
    )
    outlet_temp_liq = ScanVariable(
        "outlet_temp_liq", "Outlet Temperature Liquid Metal Breeder/Coolant (K)", 69
    )
    blpressure_liq = ScanVariable(
        "blpressure_liq", "Blanket liquid metal breeder/coolant pressure (Pa)", 70
    )
    n_liq_recirc = ScanVariable(
        "n_liq_recirc",
        "Selected number of liquid metal breeder recirculations per day",
        71,
    )
    bz_channel_conduct_liq = ScanVariable(
        "bz_channel_conduct_liq",
        "Conductance of liquid metal breeder duct walls (A V-1 m-1)",
        72,
    )
    pnuc_fw_ratio_dcll = ScanVariable(
        "pnuc_fw_ratio_dcll",
        "Ratio of FW nuclear power as fraction of total (FW+BB)",
        73,
    )
    f_nuc_pow_bz_struct = ScanVariable(
        "f_nuc_pow_bz_struct",
        "Fraction of BZ power cooled by primary coolant for dual-coolant blanket",
        74,
    )
    dx_fw_module = ScanVariable(
        "dx_fw_module", "dx_fw_module of first wall cooling channels (m)", 75
    )
    eta_turbine = ScanVariable("eta_turbine", "Thermal conversion eff.", 76)
    startupratio = ScanVariable("startupratio", "Gyrotron redundancy", 77)
    fkind = ScanVariable("fkind", "Multiplier for Nth of a kind costs", 78)
    eta_ecrh_injector_wall_plug = ScanVariable(
        "eta_ecrh_injector_wall_plug", "ECH wall plug to injector efficiency", 79
    )
    fcoolcp = ScanVariable("fcoolcp", "Coolant fraction of TF", 80)
    n_tf_coil_turns = ScanVariable("n_tf_coil_turns", "Number of turns in TF", 81)


class Scan:
    """Perform a parameter scan using the Fortran scan module."""

    def __init__(self, models, solver):
        """Immediately run the run_scan() method.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param solver: which solver to use, as specified in solver.py
        :type solver: str
        """
        self.models = models
        self.solver = solver
        self.solver_handler = SolverHandler(models, solver)
        self.run_scan()

    def run_scan(self):
        """Call a solver over a range of values of one of the variables.

        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.
        """

        if scan_variables.isweep == 0:
            # Solve single problem, rather than an array of problems (scan)
            # doopt() can also run just an evaluation
            start_time = time.time()
            ifail = self.doopt()
            write_output_files(
                models=self.models, ifail=ifail, runtime=time.time() - start_time
            )
            show_errors(constants.NOUT)
            return

        if scan_variables.isweep > scan_variables.IPNSCNS:
            raise ProcessValueError(
                "Illegal value of isweep",
                isweep=scan_variables.isweep,
                IPNSCNS=scan_variables.IPNSCNS,
            )

        if scan_variables.scan_dim == 2:
            self.scan_2d()
        else:
            self.scan_1d()

    def doopt(self):
        """Run the optimiser or solver."""
        ifail = self.solver_handler.run()
        self.post_optimise(ifail)

        return ifail

    def post_optimise(self, ifail: int):
        """Called after calling the optimising equation solver from Python.

        ifail   : input integer : error flag

        Parameters
        ----------
        ifail: int :

        """
        numerics.sqsumsq = sum(r**2 for r in numerics.rcm[: numerics.neqns]) ** 0.5

        process_output.oheadr(constants.NOUT, "Numerics")
        if self.solver == "fsolve":
            process_output.ocmmnt(
                constants.NOUT, "PROCESS has performed an fsolve (evaluation) run."
            )
        else:
            process_output.ocmmnt(
                constants.NOUT, "PROCESS has performed a VMCON (optimisation) run."
            )
        if ifail != 1:
            process_output.ovarin(constants.NOUT, "Error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.IOTTY, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            process_output.oblnkl(constants.IOTTY)

            logger.critical(f"Solver returns with ifail /= 1. {ifail=}")

            # Error code handler for VMCON
            if self.solver == "vmcon":
                self.verror(ifail)
            process_output.oblnkl(constants.NOUT)
            process_output.oblnkl(constants.IOTTY)
        else:
            # Solution found
            if self.solver != "fsolve":
                process_output.ocmmnt(
                    constants.NOUT, "and found a feasible set of parameters."
                )
                process_output.oheadr(
                    constants.IOTTY, "PROCESS found a feasible solution"
                )
            else:
                process_output.ocmmnt(
                    constants.NOUT, "and found a consistent set of parameters."
                )
                process_output.oheadr(
                    constants.IOTTY, "PROCESS found a consistent solution"
                )
            process_output.oblnkl(constants.NOUT)
            process_output.ovarin(constants.NOUT, "Error flag", "(ifail)", ifail)

            if numerics.sqsumsq >= 1.0e-2:
                process_output.oblnkl(constants.NOUT)
                process_output.ocmmnt(
                    constants.NOUT,
                    "WARNING: Constraint residues are HIGH; consider re-running",
                )
                process_output.ocmmnt(
                    constants.NOUT,
                    "   with lower values of EPSVMC to confirm convergence...",
                )
                process_output.ocmmnt(
                    constants.NOUT,
                    "   (should be able to get down to about 1.0E-8 okay)",
                )
                process_output.oblnkl(constants.NOUT)
                process_output.ocmmnt(
                    constants.IOTTY,
                    "WARNING: Constraint residues are HIGH; consider re-running",
                )
                process_output.ocmmnt(
                    constants.IOTTY,
                    "   with lower values of EPSVMC to confirm convergence...",
                )
                process_output.ocmmnt(
                    constants.IOTTY,
                    "   (should be able to get down to about 1.0E-8 okay)",
                )
                process_output.oblnkl(constants.IOTTY)

                logger.warning(f"High final constraint residues. {numerics.sqsumsq=}")

        process_output.ovarin(
            constants.NOUT, "Number of iteration variables", "(nvar)", numerics.nvar
        )
        process_output.ovarin(
            constants.NOUT,
            "Number of constraints (total)",
            "(neqns+nineqns)",
            numerics.neqns + numerics.nineqns,
        )
        process_output.ovarin(
            constants.NOUT, "Optimisation switch", "(ioptimz)", numerics.ioptimz
        )
        # Objective function output: none for fsolve
        if self.solver != "fsolve":
            process_output.ovarin(
                constants.NOUT, "Figure of merit switch", "(minmax)", numerics.minmax
            )

            objf_name = f'"{numerics.lablmm[abs(numerics.minmax) - 1]}"'

            numerics.objf_name = objf_name

            process_output.ovarst(
                constants.NOUT,
                "Objective function name",
                "(objf_name)",
                numerics.objf_name,
            )
            process_output.ovarre(
                constants.NOUT,
                "Normalised objective function",
                "(norm_objf)",
                numerics.norm_objf,
                "OP ",
            )

        process_output.ovarre(
            constants.NOUT,
            "Square root of the sum of squares of the constraint residuals",
            "(sqsumsq)",
            numerics.sqsumsq,
            "OP ",
        )
        if self.solver != "fsolve":
            process_output.ovarre(
                constants.NOUT,
                "VMCON convergence parameter",
                "(convergence_parameter)",
                global_variables.convergence_parameter,
                "OP ",
            )
            process_output.ovarin(
                constants.NOUT,
                "Number of optimising solver iterations",
                "(nviter)",
                numerics.nviter,
                "OP ",
            )
        process_output.oblnkl(constants.NOUT)

        if self.solver == "fsolve":
            if ifail == 1:
                msg = "PROCESS has solved using fsolve."
            else:
                msg = "PROCESS failed to solve using fsolve."
            process_output.write(
                constants.NOUT,
                f"{msg}\n",
            )
        else:
            if ifail == 1:
                string1 = "PROCESS has successfully optimised"
            else:
                string1 = "PROCESS has failed to optimise"

            string2 = "minimise" if numerics.minmax > 0 else "maximise"

            process_output.write(
                constants.NOUT,
                f"{string1} the optimisation parameters to {string2} the objective function: {objf_name}\n",
            )

        written_warning = False

        # Output optimisation parameters
        solution_vector_table = []
        for i in range(numerics.nvar):
            numerics.xcs[i] = numerics.xcm[i] * numerics.scafc[i]

            name = numerics.lablxc[numerics.ixc[i] - 1]
            solution_vector_table.append([name, numerics.xcs[i], numerics.xcm[i]])

            xminn = 1.01 * numerics.itv_scaled_lower_bounds[i]
            xmaxx = 0.99 * numerics.itv_scaled_upper_bounds[i]

            # Write to output file if close to optimisation parameter bounds
            if numerics.xcm[i] < xminn or numerics.xcm[i] > xmaxx:
                if not written_warning:
                    written_warning = True
                    process_output.ocmmnt(
                        constants.NOUT,
                        (
                            "Certain operating limits have been reached,"
                            "\n as shown by the following optimisation parameters that are"
                            "\n at or near to the edge of their prescribed range :\n"
                        ),
                    )

                xcval = numerics.xcm[i] * numerics.scafc[i]

                if numerics.xcm[i] < xminn:
                    location, bound = "below", "lower"
                    bounds = numerics.itv_scaled_lower_bounds
                else:
                    location, bound = "above", "upper"
                    bounds = numerics.itv_scaled_upper_bounds
                process_output.write(
                    constants.NOUT,
                    f"   {name:<30}= {xcval} is at or {location} its {bound} bound:"
                    f" {bounds[i] * numerics.scafc[i]}",
                )

            # Write optimisation parameters to mfile
            process_output.ovarre(
                constants.MFILE,
                numerics.lablxc[numerics.ixc[i] - 1],
                f"(itvar{i + 1:03d})",
                numerics.xcs[i],
            )

            if numerics.boundu[i] == numerics.boundl[i]:
                xnorm = 1.0
            else:
                xnorm = min(
                    max(
                        (numerics.xcm[i] - numerics.itv_scaled_lower_bounds[i])
                        / (
                            numerics.itv_scaled_upper_bounds[i]
                            - numerics.itv_scaled_lower_bounds[i]
                        ),
                        0.0,
                    ),
                    1.0,
                )

            process_output.ovarre(
                constants.MFILE,
                f"{name} (final value/initial value)",
                f"(xcm{i + 1:03d})",
                numerics.xcm[i],
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} (range normalised)",
                f"(nitvar{i + 1:03d})",
                xnorm,
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} (upper bound)",
                f"(boundu{i + 1:03d})",
                numerics.itv_scaled_upper_bounds[i] * numerics.scafc[i],
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} (lower bound)",
                f"(boundl{i + 1:03d})",
                numerics.itv_scaled_lower_bounds[i] * numerics.scafc[i],
            )

        # Write optimisation parameter headings to output file
        process_output.osubhd(
            constants.NOUT, "The solution vector is comprised as follows :"
        )
        process_output.write(
            constants.NOUT,
            tabulate(
                solution_vector_table,
                headers=["", "Final value", "Final / initial"],
                numalign="left",
            ),
        )

        process_output.osubhd(
            constants.NOUT,
            "The following equality constraint residues should be close to zero:",
        )

        con1, con2, err, sym, lab = constraints.constraint_eqns(
            numerics.neqns + numerics.nineqns, -1
        )

        # Write equality constraints to mfile
        equality_constraint_table = []
        for i in range(numerics.neqns):
            name = numerics.lablcc[numerics.icc[i] - 1]

            equality_constraint_table.append([
                name,
                sym[i],
                f"{con2[i]} {lab[i]}",
                f"{err[i]} {lab[i]}",
                con1[i],
            ])
            process_output.ovarre(
                constants.MFILE,
                f"{name:<33} normalised residue",
                f"(eq_con{numerics.icc[i]:03d})",
                con1[i],
            )

            process_output.ovarre(
                constants.MFILE,
                f"{name:<33} residual",
                f"(res_eq_con{numerics.icc[i]:03d})",
                err[i],
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} constraint value",
                f"(val_eq_con{numerics.icc[i]:03d})",
                con2[i],
            )

            process_output.ovarre(
                constants.MFILE,
                f"{name} units",
                f"(eq_units_con{numerics.icc[i]:03d})",
                f"'{lab[i]}'",
            )

        # Write equality constraints to output file
        process_output.write(
            constants.NOUT,
            tabulate(
                equality_constraint_table,
                headers=[
                    "",
                    "",
                    "Physical constraint",
                    "Constraint residue",
                    "Normalised residue",
                ],
                numalign="left",
            ),
        )

        # Write inequality constraints
        if numerics.nineqns > 0:
            inequality_constraint_table = []
            # Inequalities not necessarily satisfied when evaluating
            process_output.osubhd(
                constants.NOUT,
                "Negative inequality constraint (normalised) residuals indicate a constraint is satisfied.",
            )
            if self.solver == "fsolve":
                process_output.osubhd(
                    constants.NOUT,
                    "This MFile was produced via an evaluation, not an optimisation, and so the constraints "
                    "might be violated.",
                )

            for i in range(numerics.neqns, numerics.neqns + numerics.nineqns):
                name = numerics.lablcc[numerics.icc[i] - 1]
                constraint = constraints.ConstraintManager.evaluate_constraint(
                    int(numerics.icc[i])
                )

                inequality_constraint_table.append([
                    name,
                    f"{constraint.constraint_value} {constraint.units}",
                    constraint.symbol,
                    f"{constraint.constraint_bound} {constraint.units}",
                    f"{constraint.residual} {constraint.units}",
                    f"{constraint.normalised_residual}",
                ])
                process_output.ovarre(
                    constants.MFILE,
                    f"{name} normalised residue",
                    f"(ineq_con{numerics.icc[i]:03d})",
                    -constraint.normalised_residual,
                )
                process_output.ovarre(
                    constants.MFILE,
                    f"{name} physical value",
                    f"(ineq_value_con{numerics.icc[i]:03d})",
                    constraint.constraint_value,
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} symbol",
                    f"(ineq_symbol_con{numerics.icc[i]:03d})",
                    f"'{constraint.symbol}'",
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} units",
                    f"(ineq_units_con{numerics.icc[i]:03d})",
                    f"'{constraint.units}'",
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} physical bound",
                    f"(ineq_bound_con{numerics.icc[i]:03d})",
                    constraint.constraint_bound,
                )

            process_output.write(
                constants.NOUT,
                tabulate(
                    inequality_constraint_table,
                    headers=[
                        "",
                        "Physical constraint",
                        "",
                        "Physical constraint bound",
                        "Constraint residue",
                        "Normalised residue",
                    ],
                    numalign="left",
                ),
            )

    def verror(self, ifail: int):
        """Routine to print out relevant messages in the case of an
        unfeasible result from a VMCON (optimisation) run

        ifail  : input integer : error flag
        This routine prints out relevant messages in the case of
        an unfeasible result from a VMCON (optimisation) run.

        Parameters
        ----------
        ifail: int :

        """
        if ifail == -1:
            process_output.ocmmnt(constants.NOUT, "User-terminated execution of VMCON.")
            process_output.ocmmnt(constants.IOTTY, "User-terminated execution of VMCON.")
        elif ifail == 0:
            process_output.ocmmnt(
                constants.NOUT, "Improper input parameters to the VMCON routine."
            )
            process_output.ocmmnt(constants.NOUT, "PROCESS coding must be checked.")

            process_output.ocmmnt(
                constants.IOTTY, "Improper input parameters to the VMCON routine."
            )
            process_output.ocmmnt(constants.IOTTY, "PROCESS coding must be checked.")
        elif ifail == 2:
            process_output.ocmmnt(
                constants.NOUT,
                "The maximum number of calls has been reached without solution.",
            )
            process_output.ocmmnt(
                constants.NOUT,
                "The code may be stuck in a minimum in the residual space that is significantly above zero.",
            )
            process_output.oblnkl(constants.NOUT)
            process_output.ocmmnt(
                constants.NOUT, "There is either no solution possible, or the code"
            )
            process_output.ocmmnt(
                constants.NOUT, "is failing to escape from a deep local minimum."
            )
            process_output.ocmmnt(
                constants.NOUT,
                "Try changing the variables in IXC, or modify their initial values.",
            )

            process_output.ocmmnt(
                constants.IOTTY,
                "The maximum number of calls has been reached without solution.",
            )
            process_output.ocmmnt(
                constants.IOTTY,
                "The code may be stuck in a minimum in the residual space that is significantly above zero.",
            )
            process_output.oblnkl(constants.NOUT)
            process_output.oblnkl(constants.IOTTY)
            process_output.ocmmnt(
                constants.IOTTY, "There is either no solution possible, or the code"
            )
            process_output.ocmmnt(
                constants.IOTTY, "is failing to escape from a deep local minimum."
            )
            process_output.ocmmnt(
                constants.IOTTY,
                "Try changing the variables in IXC, or modify their initial values.",
            )
        elif ifail == 3:
            process_output.ocmmnt(
                constants.NOUT, "The line search required the maximum of 10 calls."
            )
            process_output.ocmmnt(
                constants.NOUT, "A feasible solution may be difficult to achieve."
            )
            process_output.ocmmnt(
                constants.NOUT, "Try changing or adding variables to IXC."
            )

            process_output.ocmmnt(
                constants.IOTTY, "The line search required the maximum of 10 calls."
            )
            process_output.ocmmnt(
                constants.IOTTY, "A feasible solution may be difficult to achieve."
            )
            process_output.ocmmnt(
                constants.IOTTY, "Try changing or adding variables to IXC."
            )
        elif ifail == 4:
            process_output.ocmmnt(
                constants.NOUT, "An uphill search direction was found."
            )
            process_output.ocmmnt(
                constants.NOUT, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.NOUT, "adding new variables to IXC.")

            process_output.ocmmnt(
                constants.IOTTY, "An uphill search direction was found."
            )
            process_output.ocmmnt(
                constants.IOTTY, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.IOTTY, "adding new variables to IXC.")
        elif ifail == 5:
            process_output.ocmmnt(
                constants.NOUT, "The quadratic programming technique was unable to"
            )
            process_output.ocmmnt(constants.NOUT, "find a feasible point.")
            process_output.oblnkl(constants.NOUT)
            process_output.ocmmnt(
                constants.NOUT, "Try changing or adding variables to IXC, or modify"
            )
            process_output.ocmmnt(
                constants.NOUT,
                "their initial values (especially if only 1 optimisation",
            )
            process_output.ocmmnt(constants.NOUT, "iteration was performed).")

            process_output.ocmmnt(
                constants.IOTTY, "The quadratic programming technique was unable to"
            )
            process_output.ocmmnt(constants.IOTTY, "find a feasible point.")
            process_output.oblnkl(constants.IOTTY)
            process_output.ocmmnt(
                constants.IOTTY, "Try changing or adding variables to IXC, or modify"
            )
            process_output.ocmmnt(
                constants.IOTTY,
                "their initial values (especially if only 1 optimisation",
            )
            process_output.ocmmnt(constants.IOTTY, "iteration was performed).")
        elif ifail == 6:
            process_output.ocmmnt(
                constants.NOUT, "The quadratic programming technique was restricted"
            )
            process_output.ocmmnt(
                constants.NOUT, "by an artificial bound, or failed due to a singular"
            )
            process_output.ocmmnt(constants.NOUT, "matrix.")
            process_output.ocmmnt(
                constants.NOUT, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.NOUT, "adding new variables to IXC.")

            process_output.ocmmnt(
                constants.IOTTY, "The quadratic programming technique was restricted"
            )
            process_output.ocmmnt(
                constants.IOTTY, "by an artificial bound, or failed due to a singular"
            )
            process_output.ocmmnt(constants.IOTTY, "matrix.")
            process_output.ocmmnt(
                constants.IOTTY, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.IOTTY, "adding new variables to IXC.")

    def scan_1d(self):
        """Run a 1-D scan."""
        # initialise dict which will contain ifail values for each scan point
        scan_1d_ifail_dict = {}

        for iscan in range(1, scan_variables.isweep + 1):
            self.scan_1d_write_point_header(iscan)
            start_time = time.time()
            ifail = self.doopt()
            scan_1d_ifail_dict[iscan] = ifail
            write_output_files(
                models=self.models, ifail=ifail, runtime=time.time() - start_time
            )

            show_errors(constants.NOUT)
            logging_model_handler.clear_logs()

        # outvar now contains results
        self.scan_1d_write_plot()
        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_values = scan_variables.sweep[: scan_variables.isweep]
        nsweep_var_name, _ = self.scan_select(
            scan_variables.nsweep, scan_variables.sweep, scan_variables.isweep
        )
        converged_count = 0
        # offsets for aligning the converged/unconverged column
        max_sweep_value_length = len(str(np.max(sweep_values)).replace(".", ""))
        offsets = [
            max_sweep_value_length - len(str(sweep_val).replace(".", ""))
            for sweep_val in sweep_values
        ]
        for iscan in range(1, scan_variables.isweep + 1):
            if scan_1d_ifail_dict[iscan] == 1:
                converged_count += 1
                print(
                    f"Scan {iscan:02d}: {nsweep_var_name} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[32mCONVERGED \u001b[0m"
                )
            else:
                print(
                    f"Scan {iscan:02d}: {nsweep_var_name} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[31mUNCONVERGED \u001b[0m"
                )
        converged_percentage = converged_count / scan_variables.isweep * 100
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d(self):
        """Run a 2-D scan."""
        # Initialise intent(out) arrays
        self.scan_2d_init()
        iscan = 1

        # initialise array which will contain ifail values for each scan point
        scan_2d_ifail_list = np.zeros(
            (scan_variables.NOUTVARS, scan_variables.IPNSCNS),
            dtype=np.float64,
            order="F",
        )
        for iscan_1 in range(1, scan_variables.isweep + 1):
            for iscan_2 in range(1, scan_variables.isweep_2 + 1):
                self.scan_2d_write_point_header(iscan, iscan_1, iscan_2)
                start_time = time.time()
                ifail = self.doopt()
                write_output_files(
                    models=self.models, ifail=ifail, runtime=time.time() - start_time
                )

                show_errors(constants.NOUT)
                logging_model_handler.clear_logs()
                scan_2d_ifail_list[iscan_1][iscan_2] = ifail
                iscan = iscan + 1

        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_1_values = scan_variables.sweep[: scan_variables.isweep]
        sweep_2_values = scan_variables.sweep_2[: scan_variables.isweep_2]
        nsweep_var_name, _ = self.scan_select(
            scan_variables.nsweep, scan_variables.sweep, scan_variables.isweep
        )
        nsweep_2_var_name, _ = self.scan_select(
            scan_variables.nsweep_2, scan_variables.sweep_2, scan_variables.isweep_2
        )
        converged_count = 0
        scan_point = 1
        # offsets for aligning the converged/unconverged column
        max_sweep1_value_length = len(str(np.max(sweep_1_values)).replace(".", ""))
        max_sweep2_value_length = len(str(np.max(sweep_2_values)).replace(".", ""))
        offsets = np.zeros(
            (scan_variables.isweep, scan_variables.isweep_2), dtype=int, order="F"
        )
        for count1, sweep1 in enumerate(sweep_1_values):
            for count2, sweep2 in enumerate(sweep_2_values):
                offsets[count1][count2] = (
                    max_sweep1_value_length
                    - len(str(sweep1).replace(".", ""))
                    + max_sweep2_value_length
                    - len(str(sweep2).replace(".", ""))
                )

        for iscan_1 in range(1, scan_variables.isweep + 1):
            for iscan_2 in range(1, scan_variables.isweep_2 + 1):
                if scan_2d_ifail_list[iscan_1][iscan_2] == 1:
                    converged_count += 1
                    print(
                        f"Scan {scan_point:02d}: ({nsweep_var_name} = {sweep_1_values[iscan_1 - 1]}, {nsweep_2_var_name} = {sweep_2_values[iscan_2 - 1]}) "
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[32mCONVERGED \u001b[0m"
                    )
                    scan_point += 1
                else:
                    print(
                        f"Scan {scan_point:02d}: ({nsweep_var_name} = {sweep_1_values[iscan_1 - 1]}, {nsweep_2_var_name} = {sweep_2_values[iscan_2 - 1]}) "
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[31mUNCONVERGED \u001b[0m"
                    )
                    scan_point += 1
        converged_percentage = (
            converged_count / (scan_variables.isweep * scan_variables.isweep_2) * 100
        )
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d_init(self):
        process_output.ovarin(
            constants.MFILE,
            "Number of first variable scan points",
            "(isweep)",
            scan_variables.isweep,
        )
        process_output.ovarin(
            constants.MFILE,
            "Number of second variable scan points",
            "(isweep_2)",
            scan_variables.isweep_2,
        )
        process_output.ovarin(
            constants.MFILE,
            "Scanning first variable number",
            "(nsweep)",
            scan_variables.nsweep,
        )
        process_output.ovarin(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_variables.nsweep_2,
        )
        process_output.ovarin(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_variables.nsweep_2,
        )
        process_output.ovarin(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_variables.nsweep_2,
        )

    def scan_1d_write_point_header(self, iscan: int):
        global_variables.iscan_global = iscan
        global_variables.vlabel, global_variables.xlabel = self.scan_select(
            scan_variables.nsweep, scan_variables.sweep, iscan
        )

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** Scan point {iscan} of {scan_variables.isweep} : {global_variables.xlabel}"
            f", {global_variables.vlabel} = {scan_variables.sweep[iscan - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarin(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan} of {scan_variables.isweep} : "
            f"{global_variables.xlabel} , {global_variables.vlabel}"
            f" = {scan_variables.sweep[iscan - 1]}"
        )

    def scan_2d_write_point_header(self, iscan, iscan_1, iscan_2):
        iscan_r = scan_variables.isweep_2 - iscan_2 + 1 if iscan_1 % 2 == 0 else iscan_2

        # Makes iscan available globally (read-only)
        global_variables.iscan_global = iscan

        global_variables.vlabel, global_variables.xlabel = self.scan_select(
            scan_variables.nsweep, scan_variables.sweep, iscan_1
        )
        global_variables.vlabel_2, global_variables.xlabel_2 = self.scan_select(
            scan_variables.nsweep_2, scan_variables.sweep_2, iscan_r
        )

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** 2D Scan point {iscan} of {scan_variables.isweep * scan_variables.isweep_2} : "
            f"{global_variables.vlabel} = {scan_variables.sweep[iscan_1 - 1]} and"
            f" {global_variables.vlabel_2} = {scan_variables.sweep_2[iscan_r - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarin(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {global_variables.xlabel}, "
            f"{global_variables.vlabel} = {scan_variables.sweep[iscan_1 - 1]}"
            f" and {global_variables.xlabel_2}, "
            f"{global_variables.vlabel_2} = {scan_variables.sweep_2[iscan_r - 1]} "
        )

        return iscan_r

    def scan_1d_write_plot(self):
        if scan_variables.first_call_1d:
            process_output.ovarin(
                constants.MFILE,
                "Number of scan points",
                "(isweep)",
                scan_variables.isweep,
            )
            process_output.ovarin(
                constants.MFILE,
                "Scanning variable number",
                "(nsweep)",
                scan_variables.nsweep,
            )

            scan_variables.first_call_1d = False

    def scan_select(self, nwp, swp, iscn):
        match nwp:
            case 1:
                physics_variables.aspect = swp[iscn - 1]
            case 2:
                divertor_variables.pflux_div_heat_load_max_mw = swp[iscn - 1]
            case 3:
                constraint_variables.p_plant_electric_net_required_mw = swp[iscn - 1]
            case 4:
                physics_variables.hfact = swp[iscn - 1]
            case 5:
                tfcoil_variables.oacdcp = swp[iscn - 1]
            case 6:
                constraint_variables.pflux_fw_neutron_max_mw = swp[iscn - 1]
            case 7:
                physics_variables.beamfus0 = swp[iscn - 1]
            case 9:
                physics_variables.temp_plasma_electron_vol_avg_kev = swp[iscn - 1]
            case 10:
                numerics.boundu[14] = swp[iscn - 1]
            case 11:
                physics_variables.beta_norm_max = swp[iscn - 1]
            case 12:
                current_drive_variables.f_c_plasma_bootstrap_max = swp[iscn - 1]
            case 13:
                numerics.boundu[9] = swp[iscn - 1]
            case 16:
                physics_variables.rmajor = swp[iscn - 1]
            case 17:
                constraint_variables.b_tf_inboard_max = swp[iscn - 1]
            case 18:
                constraint_variables.eta_cd_norm_hcd_primary_max = swp[iscn - 1]
            case 19:
                numerics.boundl[15] = swp[iscn - 1]
            case 20:
                constraint_variables.t_burn_min = swp[iscn - 1]
            case 22:
                if cost_variables.i_plant_availability == 1:
                    raise ProcessValueError(
                        "Do not scan f_t_plant_available if i_plant_availability=1"
                    )
                cost_variables.f_t_plant_available = swp[iscn - 1]
            case 24:
                constraint_variables.p_fusion_total_max_mw = swp[iscn - 1]
            case 25:
                physics_variables.kappa = swp[iscn - 1]
            case 26:
                physics_variables.triang = swp[iscn - 1]
            case 27:
                constraint_variables.tbrmin = swp[iscn - 1]
            case 28:
                physics_variables.b_plasma_toroidal_on_axis = swp[iscn - 1]
            case 29:
                impurity_radiation_module.coreradius = swp[iscn - 1]
            case 31:
                constraint_variables.f_alpha_energy_confinement_min = swp[iscn - 1]
            case 32:
                numerics.epsvmc = swp[iscn - 1]
            case 38:
                numerics.boundu[128] = swp[iscn - 1]
            case 39:
                numerics.boundu[130] = swp[iscn - 1]
            case 40:
                numerics.boundu[134] = swp[iscn - 1]
            case 41:
                build_variables.dr_blkt_outboard = swp[iscn - 1]
            case 42:
                impurity_radiation_module.f_nd_impurity_electrons[8] = swp[iscn - 1]
                impurity_radiation_module.f_nd_impurity_electron_array[8] = (
                    impurity_radiation_module.f_nd_impurity_electrons[8]
                )
            case 44:
                tfcoil_variables.sig_tf_case_max = swp[iscn - 1]
            case 45:
                tfcoil_variables.temp_tf_superconductor_margin_min = swp[iscn - 1]
            case 46:
                numerics.boundu[151] = swp[iscn - 1]
            case 48:
                tfcoil_variables.n_tf_wp_pancakes = int(swp[iscn - 1])
            case 49:
                tfcoil_variables.n_tf_wp_layers = int(swp[iscn - 1])
            case 50:
                impurity_radiation_module.f_nd_impurity_electrons[12] = swp[iscn - 1]
                impurity_radiation_module.f_nd_impurity_electron_array[12] = (
                    impurity_radiation_module.f_nd_impurity_electrons[12]
                )
            case 51:
                physics_variables.f_p_div_lower = swp[iscn - 1]
            case 52:
                physics_variables.rad_fraction_sol = swp[iscn - 1]
            case 53:
                numerics.boundu[156] = swp[iscn - 1]
            case 54:
                tfcoil_variables.b_crit_upper_nbti = swp[iscn - 1]
            case 55:
                build_variables.dr_shld_inboard = swp[iscn - 1]
            case 56:
                heat_transport_variables.p_cryo_plant_electric_max_mw = swp[iscn - 1]
            case 57:
                numerics.boundl[1] = swp[iscn - 1]
            case 58:
                build_variables.dr_fw_plasma_gap_inboard = swp[iscn - 1]
            case 59:
                build_variables.dr_fw_plasma_gap_outboard = swp[iscn - 1]
            case 60:
                tfcoil_variables.sig_tf_wp_max = swp[iscn - 1]
            case 61:
                rebco_variables.copperaoh_m2_max = swp[iscn - 1]
            case 62:
                pfcoil_variables.coheof = swp[iscn - 1]
            case 63:
                build_variables.dr_cs = swp[iscn - 1]
            case 64:
                pfcoil_variables.ohhghf = swp[iscn - 1]
            case 65:
                cs_fatigue_variables.n_cycle_min = swp[iscn - 1]
            case 66:
                pfcoil_variables.oh_steel_frac = swp[iscn - 1]
            case 67:
                cs_fatigue_variables.t_crack_vertical = swp[iscn - 1]
            case 68:
                fwbs_variables.inlet_temp_liq = swp[iscn - 1]
            case 69:
                fwbs_variables.outlet_temp_liq = swp[iscn - 1]
            case 70:
                fwbs_variables.blpressure_liq = swp[iscn - 1]
            case 71:
                fwbs_variables.n_liq_recirc = swp[iscn - 1]
            case 72:
                fwbs_variables.bz_channel_conduct_liq = swp[iscn - 1]
            case 73:
                fwbs_variables.pnuc_fw_ratio_dcll = swp[iscn - 1]
            case 74:
                fwbs_variables.f_nuc_pow_bz_struct = swp[iscn - 1]
            case 75:
                fwbs_variables.dx_fw_module = swp[iscn - 1]
            case 76:
                heat_transport_variables.eta_turbine = swp[iscn - 1]
            case 77:
                cost_variables.startupratio = swp[iscn - 1]
            case 78:
                cost_variables.fkind = swp[iscn - 1]
            case 79:
                current_drive_variables.eta_ecrh_injector_wall_plug = swp[iscn - 1]
            case 80:
                tfcoil_variables.fcoolcp = swp[iscn - 1]
            case 81:
                tfcoil_variables.n_tf_coil_turns = swp[iscn - 1]
            case _:
                raise ProcessValueError("Illegal scan variable number", nwp=nwp)

        return ScanVariables(int(nwp)).value
