from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from enum import Enum
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numpy as np
from tabulate import tabulate

from process.core import constants, process_output
from process.core.caller import write_output_files
from process.core.exceptions import ProcessValueError
from process.core.log import logging_model_handler, show_errors
from process.core.solver import constraints
from process.core.solver.solver_handler import SolverHandler
from process.data_structure.numerics import FiguresOfMerit, PROCESSRunMode
from process.data_structure.scan_variables import IPNSCNS, NOUTVARS

if TYPE_CHECKING:
    from process.core.model import DataStructure, Model

logger = logging.getLogger(__name__)


@dataclass
class ScanVariable:
    variable_area: SVE
    variable_num: int


class SVE(Enum):
    P = "physics"
    D = "divertor"
    C = "constraints"
    T = "tfcoil"
    TR = "rebco"
    CD = "current_drive"
    NUM = "numerics"
    CST = "costs"
    IR = "impurity_radiation"
    B = "build"
    HT = "heat_transport"
    PF = "pf_coil"
    CS = "cs_fatigue"
    FWBS = "fwbs"


class ScanVariables(Enum):
    @classmethod
    def _missing_(cls, value):
        if isinstance(value, int):
            for sv in cls:
                if sv.number == value:
                    return sv
        raise ProcessValueError("Illegal scan variable number", nwp=value)

    def full_name(self):
        if "__" in self.name:
            return self.name.replace("__", "(") + ")"
        return self.name

    def set(self, data, sweep):
        var_area = getattr(data, self.area.value)

        if self.value.variable_num == 22 and var_area.i_plant_availability == 1:
            raise ProcessValueError(
                "Do not scan f_t_plant_available if i_plant_availability=1"
            )

        if "__" in self.name:
            name, index = self.name.split("__")
            getattr(var_area, name)[index] = sweep
            if name == "f_nd_impurity_electrons":
                var_area.f_nd_impurity_electron_array[int(index - 1)] = sweep
        else:
            setattr(var_area, self.name, sweep)

        self._data_ = getattr(var_area, name)

    @DynamicClassAttribute
    def area(self):
        return self.value.variable_area

    @DynamicClassAttribute
    def number(self):
        return self.value.variable_num

    @DynamicClassAttribute
    def data(self):
        if hasattr(self, "_data_"):
            return self._data_
        raise ValueError("Data not available")

    aspect = ScanVariable(SVE.P, 1)
    pflux_div_heat_load_max_mw = ScanVariable(SVE.D, 2)
    p_plant_electric_net_required_mw = ScanVariable(SVE.C, 3)
    hfact = ScanVariable(SVE.P, 4)
    j_tf_coil_full_area = ScanVariable(SVE.T, 5)
    pflux_fw_neutron_max_mw = ScanVariable(SVE.C, 6)
    beamfus0 = ScanVariable(SVE.P, 7)
    temp_plasma_electron_vol_avg_kev = ScanVariable(SVE.P, 9)
    boundu__14 = ScanVariable(SVE.NUM, 10)
    beta_norm_max = ScanVariable(SVE.P, 11)
    f_c_plasma_bootstrap_max = ScanVariable(SVE.CD, 12)
    boundu__10 = ScanVariable(SVE.NUM, 13)
    rmajor = ScanVariable(SVE.P, 16)
    b_tf_inboard_max = ScanVariable(SVE.C, 17)
    eta_cd_norm_hcd_primary_max = ScanVariable(SVE.C, 18)
    boundl__16 = ScanVariable(SVE.NUM, 19)
    t_burn_min = ScanVariable(SVE.C, 20)
    f_t_plant_available = ScanVariable(SVE.CST, 22)
    p_fusion_total_max_mw = ScanVariable(SVE.C, 24)
    kappa = ScanVariable(SVE.P, 25)
    triang = ScanVariable(SVE.P, 26)
    tbrmin = ScanVariable(SVE.C, 27)
    b_plasma_toroidal_on_axis = ScanVariable(SVE.P, 28)
    coreradius = ScanVariable(SVE.IR, 29)
    f_alpha_energy_confinement_min = ScanVariable(SVE.C, 31)
    epsvmc = ScanVariable(SVE.NUM, 32)
    boundu__129 = ScanVariable(SVE.NUM, 38)
    boundu__131 = ScanVariable(SVE.NUM, 39)
    boundu__135 = ScanVariable(SVE.NUM, 40)
    dr_blkt_outboard = ScanVariable(SVE.B, 41)
    f_nd_impurity_electrons__9 = ScanVariable(SVE.IR, 42)
    sig_tf_case_max = ScanVariable(SVE.T, 44)
    temp_tf_superconductor_margin_min = ScanVariable(SVE.T, 45)
    boundu__152 = ScanVariable(SVE.NUM, 46)
    n_tf_wp_pancakes = ScanVariable(SVE.T, 48)
    n_tf_wp_layers = ScanVariable(SVE.T, 49)
    f_nd_impurity_electrons__13 = ScanVariable(SVE.IR, 50)
    f_p_div_lower = ScanVariable(SVE.P, 51)
    rad_fraction_sol = ScanVariable(SVE.P, 52)
    boundu__157 = ScanVariable(SVE.NUM, 53)
    b_crit_upper_nbti = ScanVariable(SVE.T, 54)
    dr_shld_inboard = ScanVariable(SVE.B, 55)
    p_cryo_plant_electric_max_mw = ScanVariable(SVE.HT, 56)
    boundl__2 = ScanVariable(SVE.NUM, 57)
    dr_fw_plasma_gap_inboard = ScanVariable(SVE.B, 58)
    dr_fw_plasma_gap_outboard = ScanVariable(SVE.B, 59)
    sig_tf_wp_max = ScanVariable(SVE.T, 60)
    copperaoh_m2_max = ScanVariable(SVE.TR, 61)
    coheof = ScanVariable(SVE.PF, 62)
    dr_cs = ScanVariable(SVE.B, 63)
    ohhghf = ScanVariable(SVE.PF, 64)
    n_cycle_min = ScanVariable(SVE.CS, 65)
    oh_steel_frac = ScanVariable(SVE.PF, 66)
    t_crack_vertical = ScanVariable(SVE.CS, 67)
    inlet_temp_liq = ScanVariable(SVE.FWBS, 68)
    outlet_temp_liq = ScanVariable(SVE.FWBS, 69)
    blpressure_liq = ScanVariable(SVE.FWBS, 70)
    n_liq_recirc = ScanVariable(SVE.FWBS, 71)
    bz_channel_conduct_liq = ScanVariable(SVE.FWBS, 72)
    pnuc_fw_ratio_dcll = ScanVariable(SVE.FWBS, 73)
    f_nuc_pow_bz_struct = ScanVariable(SVE.FWBS, 74)
    dx_fw_module = ScanVariable(SVE.FWBS, 75)
    eta_turbine = ScanVariable(SVE.HT, 76)
    startupratio = ScanVariable(SVE.CST, 77)
    fkind = ScanVariable(SVE.CST, 78)
    eta_ecrh_injector_wall_plug = ScanVariable(SVE.CD, 79)
    fcoolcp = ScanVariable(SVE.T, 80)
    n_tf_coil_turns = ScanVariable(SVE.T, 81)


class Scan:
    """Perform a parameter scan using the Fortran scan module."""

    def __init__(self, models: Model, solver: str, data: DataStructure):
        """Immediately run the run_scan() method.

        Parameters
        ----------
        models :
            Physics and engineering model objects
        solver :
            Which solver to use, as specified in solver.py
        data :
            Data structure object
        """
        self.models = models
        self.solver = solver
        self.data = data
        self.solver_handler = SolverHandler(models, solver, data)
        self.run_scan()

    def run_scan(self):
        """Call a solver over a range of values of one of the variables.

        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.
        """
        if self.data.scan.isweep == 0:
            # Solve single problem, rather than an array of problems (scan)
            # doopt() can also run just an evaluation
            start_time = time.time()
            ifail = self.doopt()
            write_output_files(
                models=self.models,
                data=self.data,
                ifail=ifail,
                runtime=time.time() - start_time,
            )
            show_errors(constants.NOUT)
            return

        if self.data.scan.isweep > IPNSCNS:
            raise ProcessValueError(
                "Illegal value of isweep",
                isweep=self.data.scan.isweep,
                IPNSCNS=IPNSCNS,
            )

        if self.data.scan.scan_dim == 2:
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
        self.data.numerics.sqsumsq = (
            sum(r**2 for r in self.data.numerics.rcm[: self.data.numerics.neqns]) ** 0.5
        )

        process_output.oheadr(constants.NOUT, "Numerics")
        process_output.ocmmnt(
            constants.NOUT,
            f"PROCESS has performed a {'fsolve' if self.solver == 'fsolve' else 'VMCON'} (optimisation) run.",
        )
        if ifail != 1:
            process_output.ovarin(constants.NOUT, "Error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.IOTTY, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            print()

            logger.critical("Solver returns with ifail /= 1. %s", ifail)

            # Error code handler for VMCON
            if self.solver == "vmcon":
                self.verror(ifail)
            process_output.oblnkl(constants.NOUT)
            print()
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

            if self.data.numerics.sqsumsq >= 1.0e-2:
                string = (
                    "WARNING: Constraint residues are HIGH; consider re-running\n"
                    "   with lower values of EPSVMC to confirm convergence...\n"
                    "   (should be able to get down to about 1.0E-8 okay)\n"
                )
                process_output.ocmmnt(constants.NOUT, ("\n" + string))
                print(string)

                logger.warning(
                    f"High final constraint residues. {self.data.numerics.sqsumsq=}"
                )

        process_output.ovarin(
            constants.NOUT,
            "Number of iteration variables",
            "(nvar)",
            self.data.numerics.nvar,
        )
        process_output.ovarin(
            constants.NOUT,
            "Number of constraints (total)",
            "(neqns+nineqns)",
            self.data.numerics.neqns + self.data.numerics.nineqns,
        )
        process_output.ovarin(
            constants.NOUT,
            "Optimisation switch",
            "(ioptimz)",
            self.data.numerics.ioptimz,
        )
        process_output.ocmmnt(
            constants.NOUT,
            f"     {PROCESSRunMode(self.data.numerics.ioptimz).description}",
        )
        # Objective function output: none for fsolve
        if self.solver != "fsolve":
            process_output.ovarin(
                constants.NOUT,
                "Figure of merit switch",
                "(minmax)",
                self.data.numerics.minmax,
            )

            objf_name = f'"{FiguresOfMerit(abs(self.data.numerics.minmax)).description}"'

            self.data.numerics.objf_name = objf_name

            process_output.ovarst(
                constants.NOUT,
                "Objective function name",
                "(objf_name)",
                self.data.numerics.objf_name,
            )
            process_output.ovarre(
                constants.NOUT,
                "Normalised objective function",
                "(norm_objf)",
                self.data.numerics.norm_objf,
                "OP ",
            )

        process_output.ovarre(
            constants.NOUT,
            "Square root of the sum of squares of the constraint residuals",
            "(sqsumsq)",
            self.data.numerics.sqsumsq,
            "OP ",
        )
        if self.solver != "fsolve":
            process_output.ovarre(
                constants.NOUT,
                "VMCON convergence parameter",
                "(convergence_parameter)",
                self.data.globals.convergence_parameter,
                "OP ",
            )
            process_output.ovarin(
                constants.NOUT,
                "Number of optimising solver iterations",
                "(nviter)",
                self.data.numerics.nviter,
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

            string2 = "minimise" if self.data.numerics.minmax > 0 else "maximise"

            process_output.write(
                constants.NOUT,
                f"{string1} the optimisation parameters to {string2} the objective function: {objf_name}\n",
            )

        written_warning = False

        # Output optimisation parameters
        solution_vector_table = []
        for i in range(self.data.numerics.nvar):
            self.data.numerics.xcs[i] = (
                self.data.numerics.xcm[i] * self.data.numerics.scafc[i]
            )

            name = self.data.numerics.lablxc[self.data.numerics.ixc[i] - 1]
            solution_vector_table.append([
                name,
                self.data.numerics.xcs[i],
                self.data.numerics.xcm[i],
            ])

            xminn = 1.01 * self.data.numerics.itv_scaled_lower_bounds[i]
            xmaxx = 0.99 * self.data.numerics.itv_scaled_upper_bounds[i]

            # Write to output file if close to optimisation parameter bounds
            if self.data.numerics.xcm[i] < xminn or self.data.numerics.xcm[i] > xmaxx:
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

                xcval = self.data.numerics.xcm[i] * self.data.numerics.scafc[i]

                if self.data.numerics.xcm[i] < xminn:
                    location, bound = "below", "lower"
                    bounds = self.data.numerics.itv_scaled_lower_bounds
                else:
                    location, bound = "above", "upper"
                    bounds = self.data.numerics.itv_scaled_upper_bounds
                process_output.write(
                    constants.NOUT,
                    f"   {name:<30}= {xcval} is at or {location} its {bound} bound:"
                    f" {bounds[i] * self.data.numerics.scafc[i]}",
                )

            # Write optimisation parameters to mfile
            process_output.ovarre(
                constants.MFILE,
                self.data.numerics.lablxc[self.data.numerics.ixc[i] - 1],
                f"(itvar{i + 1:03d})",
                self.data.numerics.xcs[i],
            )

            if self.data.numerics.boundu[i] == self.data.numerics.boundl[i]:
                xnorm = 1.0
            else:
                xnorm = min(
                    max(
                        (
                            self.data.numerics.xcm[i]
                            - self.data.numerics.itv_scaled_lower_bounds[i]
                        )
                        / (
                            self.data.numerics.itv_scaled_upper_bounds[i]
                            - self.data.numerics.itv_scaled_lower_bounds[i]
                        ),
                        0.0,
                    ),
                    1.0,
                )

            process_output.ovarre(
                constants.MFILE,
                f"{name} (final value/initial value)",
                f"(xcm{i + 1:03d})",
                self.data.numerics.xcm[i],
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
                self.data.numerics.itv_scaled_upper_bounds[i]
                * self.data.numerics.scafc[i],
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} (lower bound)",
                f"(boundl{i + 1:03d})",
                self.data.numerics.itv_scaled_lower_bounds[i]
                * self.data.numerics.scafc[i],
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

        con1, con2, err, _, lab = constraints.constraint_eqns(
            self.data.numerics.neqns + self.data.numerics.nineqns, -1, self.data
        )

        # Write equality constraints to mfile
        equality_constraint_table = []
        for i in range(self.data.numerics.neqns):
            name = self.data.numerics.lablcc[self.data.numerics.icc[i] - 1]

            equality_constraint_table.append([
                name,
                "=",
                f"{con2[i]} {lab[i]}",
                f"{err[i]} {lab[i]}",
                con1[i],
            ])
            process_output.ovarre(
                constants.MFILE,
                f"{name:<33} normalised residue",
                f"(eq_con{self.data.numerics.icc[i]:03d})",
                con1[i],
            )

            process_output.ovarre(
                constants.MFILE,
                f"{name:<33} residual",
                f"(res_eq_con{self.data.numerics.icc[i]:03d})",
                err[i],
            )
            process_output.ovarre(
                constants.MFILE,
                f"{name} constraint value",
                f"(val_eq_con{self.data.numerics.icc[i]:03d})",
                con2[i],
            )

            process_output.ovarre(
                constants.MFILE,
                f"{name} units",
                f"(eq_units_con{self.data.numerics.icc[i]:03d})",
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
        if self.data.numerics.nineqns > 0:
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

            for i in range(
                self.data.numerics.neqns,
                self.data.numerics.neqns + self.data.numerics.nineqns,
            ):
                name = self.data.numerics.lablcc[self.data.numerics.icc[i] - 1]
                constraint = constraints.ConstraintManager.evaluate_constraint(
                    int(self.data.numerics.icc[i]), self.data
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
                    f"(ineq_con{self.data.numerics.icc[i]:03d})",
                    -constraint.normalised_residual,
                )
                process_output.ovarre(
                    constants.MFILE,
                    f"{name} physical value",
                    f"(ineq_value_con{self.data.numerics.icc[i]:03d})",
                    constraint.constraint_value,
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} symbol",
                    f"(ineq_symbol_con{self.data.numerics.icc[i]:03d})",
                    f"'{constraint.symbol}'",
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} units",
                    f"(ineq_units_con{self.data.numerics.icc[i]:03d})",
                    f"'{constraint.units}'",
                )

                process_output.ovarre(
                    constants.MFILE,
                    f"{name} physical bound",
                    f"(ineq_bound_con{self.data.numerics.icc[i]:03d})",
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

    @staticmethod
    def verror(ifail: int):
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
            strings = ("User-terminated execution of VMCON.",)
        elif ifail == 0:
            strings = (
                "Improper input parameters to the VMCON routine.",
                "PROCESS coding must be checked.",
            )
        elif ifail == 2:
            strings = (
                "The maximum number of calls has been reached without solution.",
                "The code may be stuck in a minimum in the residual space that is significantly above zero.",
                "",
                "There is either no solution possible, or the code",
                "is failing to escape from a deep local minimum.",
                "Try changing the variables in IXC, or modify their initial values.",
            )
        elif ifail == 3:
            strings = (
                "The line search required the maximum of 10 calls.",
                "A feasible solution may be difficult to achieve.",
                "Try changing or adding variables to IXC.",
            )
        elif ifail == 4:
            strings = (
                "An uphill search direction was found.",
                "Try changing the equations in ICC, or",
                "adding new variables to IXC.",
            )
        elif ifail == 5:
            strings = (
                "The quadratic programming technique was unable to",
                "find a feasible point.",
                "",
                "Try changing or adding variables to IXC, or modify",
                "their initial values (especially if only 1 optimisation",
                "iteration was performed).",
            )

        elif ifail == 6:
            strings = (
                "The quadratic programming technique was restricted",
                "by an artificial bound, or failed due to a singular",
                "matrix.",
                "Try changing the equations in ICC, or",
                "adding new variables to IXC.",
            )

        strings = "\n".join(strings)
        process_output.ocmmnt(constants.NOUT, strings)
        print(strings)

    def scan_1d(self):
        """Run a 1-D scan."""
        # initialise dict which will contain ifail values for each scan point
        scan_1d_ifail_dict = {}

        for iscan in range(1, self.data.scan.isweep + 1):
            self.scan_1d_write_point_header(iscan)
            start_time = time.time()
            ifail = self.doopt()
            scan_1d_ifail_dict[iscan] = ifail
            write_output_files(
                models=self.models,
                data=self.data,
                ifail=ifail,
                runtime=time.time() - start_time,
            )

            show_errors(constants.NOUT)
            logging_model_handler.clear_logs()

        # outvar now contains results
        self.scan_1d_write_plot(self.data.scan)
        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_values = self.data.scan.sweep[: self.data.scan.isweep]
        nsweep_var = self.scan_select(
            self.data.scan.nsweep, self.data.scan.sweep, self.data.scan.isweep
        )
        converged_count = 0
        # offsets for aligning the converged/unconverged column
        max_sweep_value_length = len(str(np.max(sweep_values)).replace(".", ""))
        offsets = [
            max_sweep_value_length - len(str(sweep_val).replace(".", ""))
            for sweep_val in sweep_values
        ]
        for iscan in range(self.data.scan.isweep):
            pstring = (
                f"Scan {iscan:02d}: {nsweep_var.name} = {sweep_values[iscan]} "
                + " " * offsets[iscan]
                + "\u001b[32m{}CONVERGED \u001b[0m"
            )
            if scan_1d_ifail_dict[iscan + 1] == 1:
                converged_count += 1
                pstring.format("")
            else:
                pstring.format("UN")
            print(pstring)
        converged_percentage = converged_count / self.data.scan.isweep * 100
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d(self):
        """Run a 2-D scan."""
        # Initialise intent(out) arrays
        self.scan_2d_init(self.data.scan)
        iscan = 1

        # initialise array which will contain ifail values for each scan point
        scan_2d_ifail_list = np.zeros(
            (NOUTVARS, IPNSCNS),
            dtype=np.float64,
            order="F",
        )
        for iscan_1 in range(1, self.data.scan.isweep + 1):
            for iscan_2 in range(1, self.data.scan.isweep_2 + 1):
                self.scan_2d_write_point_header(iscan, iscan_1, iscan_2)
                start_time = time.time()
                ifail = self.doopt()
                write_output_files(
                    models=self.models,
                    data=self.data,
                    ifail=ifail,
                    runtime=time.time() - start_time,
                )

                show_errors(constants.NOUT)
                logging_model_handler.clear_logs()
                scan_2d_ifail_list[iscan_1][iscan_2] = ifail
                iscan += 1

        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_1_values = self.data.scan.sweep[: self.data.scan.isweep]
        sweep_2_values = self.data.scan.sweep_2[: self.data.scan.isweep_2]
        nsweep_var = self.scan_select(
            self.data.scan.nsweep, self.data.scan.sweep, self.data.scan.isweep
        )
        nsweep_2_var = self.scan_select(
            self.data.scan.nsweep_2, self.data.scan.sweep_2, self.data.scan.isweep_2
        )
        converged_count = 0
        scan_point = 1
        # offsets for aligning the converged/unconverged column
        max_sweep1_value_length = len(str(np.max(sweep_1_values)).replace(".", ""))
        max_sweep2_value_length = len(str(np.max(sweep_2_values)).replace(".", ""))
        offsets = np.zeros(
            (self.data.scan.isweep, self.data.scan.isweep_2), dtype=int, order="F"
        )
        for count1, sweep1 in enumerate(sweep_1_values):
            for count2, sweep2 in enumerate(sweep_2_values):
                offsets[count1][count2] = (
                    max_sweep1_value_length
                    - len(str(sweep1).replace(".", ""))
                    + max_sweep2_value_length
                    - len(str(sweep2).replace(".", ""))
                )

        for iscan_1 in range(1, self.data.scan.isweep + 1):
            for iscan_2 in range(1, self.data.scan.isweep_2 + 1):
                string = (
                    f"Scan {scan_point:02d}: ({nsweep_var.name} = {sweep_1_values[iscan_1 - 1]},"
                    f" {nsweep_2_var.name} = {sweep_2_values[iscan_2 - 1]}) "
                    + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                    + "\u001b[32m{}CONVERGED \u001b[0m"
                )
                if scan_2d_ifail_list[iscan_1][iscan_2] == 1:
                    converged_count += 1
                    print(string.format())
                else:
                    print(string.format("UN"))
                scan_point += 1
        converged_percentage = (
            converged_count / (self.data.scan.isweep * self.data.scan.isweep_2) * 100
        )
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d_init(self):
        sv = self.data.scan
        for d, n, v in (
            ("Number of first variable scan points", "(isweep)", sv.isweep),
            ("Number of second variable scan points", "(isweep_2)", sv.isweep_2),
            ("Scanning first variable number", "(nsweep)", sv.nsweep),
            ("Scanning second variable number", "(nsweep_2)", sv.nsweep_2),
            ("Scanning second variable number", "(nsweep_2)", sv.nsweep_2),
            ("Scanning second variable number", "(nsweep_2)", sv.nsweep_2),
        ):
            process_output.ovarin(constants.MFILE, d, n, v)

    def _set_v_x_label(self, iscan, twod=False):
        if twod:
            sv = self.scan_select(self.data.scan.nsweep_2, self.data.scan.sweep_2, iscan)
        else:
            sv = self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, iscan)
        self.data.globals.vlabel.vlabel = sv.name
        self.data.globals.xlabel = sv.data.description

    def scan_1d_write_point_header(self, iscan: int):
        self.data.globals.iscan_global = iscan
        self._set_v_x_label(iscan)

        process_output.oblnkl(constants.NOUT)
        process_output.oblnkl(constants.MFILE)

        process_output.write(
            constants.NOUT,
            f"Scan point {iscan} of {self.data.scan.isweep} : {self.data.globals.xlabel}"
            f", {self.data.globals.vlabel} = {self.data.scan.sweep[iscan - 1]} "
            "",
        )
        process_output.ovarin(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan} of {self.data.scan.isweep} : "
            f"{self.data.globals.xlabel} , {self.data.globals.vlabel}"
            f" = {self.data.scan.sweep[iscan - 1]}"
        )

    def scan_2d_write_point_header(self, iscan, iscan_1, iscan_2):
        iscan_r = self.data.scan.isweep_2 - iscan_2 + 1 if iscan_1 % 2 == 0 else iscan_2

        # Makes iscan available globally (read-only)
        self.data.globals.iscan_global = iscan

        self._set_v_x_label(iscan_1)
        self._set_v_x_label(iscan_r)

        process_output.oblnkl(constants.NOUT)
        process_output.oblnkl(constants.MFILE)

        process_output.write(
            constants.NOUT,
            f"***** 2D Scan point {iscan} of {self.data.scan.isweep * self.data.scan.isweep_2} : "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]} and"
            f" {self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
            "*****",
        )
        process_output.ovarin(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {self.data.globals.xlabel}, "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]}"
            f" and {self.data.globals.xlabel_2}, "
            f"{self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
        )

        return iscan_r

    def scan_1d_write_plot(self):
        if self.data.scan.first_call_1d:
            for d, n, v in (
                ("Number of scan points", "(isweep)", self.data.scan.isweep),
                ("Scanning variable number", "(nsweep)", self.data.scan.nsweep),
            ):
                process_output.ovarin(constants.MFILE, d, n, v)

            self.data.scan.first_call_1d = False

    def scan_select(self, nwp, swp, iscn):
        sv = ScanVariables(int(nwp))
        sweep = swp[iscn - 1]
        sv.set(self.data, sweep)
        return sv
