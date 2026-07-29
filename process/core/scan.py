"""Scanning mechanics"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from enum import Enum
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numpy as np
from tabulate import tabulate

from process.core import constants, process_output
from process.core.caller import write_output_files
from process.core.exceptions import ProcessValueError
from process.core.io.data_structure_dicts import get_dicts
from process.core.log import logging_model_handler, show_errors
from process.core.solver import constraints
from process.core.solver.solver_handler import SolverHandler
from process.data_structure.numerics import FiguresOfMerit, PROCESSRunMode
from process.data_structure.scan_variables import IPNSCNS, NOUTVARS, ScanData
from process.models.availability import AvailabilityModel

if TYPE_CHECKING:
    from process.core.model import DataStructure, Model


logger = logging.getLogger(__name__)


@dataclass
class ScanVariable:
    """Scan variable container"""

    number: int
    area: Area = field(repr=False)
    _out_name_: str | None = None


class Area(Enum):
    """Scan variable data structure area shorthand

    this mirrors data_structure.area
    """

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


class ScanVariables(ScanVariable, Enum):
    """Scan variable options"""

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, int):
            for sv in cls:
                if sv.number == value:
                    return sv
            raise ProcessValueError("Illegal scan variable number", nwp=value)

        if isinstance(value, str):
            return cls[value.replace("(", "__").replace(")", "")]

        return super()._missing_(value)

    @DynamicClassAttribute
    def fname(self):
        """Full name"""
        if "__" in self.name:
            return self.name.replace("__", "(") + ")"
        return self.name

    @DynamicClassAttribute
    def out_name(self):
        """Output name"""
        if self._out_name_ is None:
            return self.fname
        return self._out_name_

    def set(self, data: DataStructure, sweep_val: float):
        """Set value of scan variable

        Raises
        ------
        ProcessValueError
            Scanning f_t_plant_available if i_plant_availability=1"

        """
        var_area = getattr(data, self.area.value)

        if (
            self.number == 22
            and AvailabilityModel(var_area.i_plant_availability)
            != AvailabilityModel.USER_INPUT
        ):
            raise ProcessValueError(
                "Do not scan f_t_plant_available if i_plant_availability=1"
            )

        if "__" in self.name:
            name, index = self.name.split("__")
            getattr(var_area, name)[int(index) - 1] = sweep_val
            if name == "f_nd_impurity_electrons":
                var_area.f_nd_impurity_electron_array[int(index - 1)] = sweep_val
        else:
            setattr(var_area, self.name, sweep_val)
            name = self.name

        self._data_ = getattr(var_area, name)
        self._description_ = get_dicts()["DICT_DESCRIPTIONS"][name]

    @DynamicClassAttribute
    def data(self):
        """Variable value

        Raises
        ------
        ValueError
            data not set
        """
        if hasattr(self, "_data_"):
            return self._data_
        raise ValueError("Data not available")

    @DynamicClassAttribute
    def description(self):
        """Variable description

        Raises
        ------
        ValueError
            description not set
        """
        if hasattr(self, "_description_"):
            return self._description_
        raise ValueError("Description not available")

    def get_val(self, mfile, scan):
        """Get value from mfile"""
        # TODO this will fail for boundu/l we should write the scan variable to
        # the mfile and use that directly (also replacign output names)
        return mfile.get(self.out_name, scan=scan)

    aspect = (1, Area.P)
    pflux_div_heat_load_max_mw = (2, Area.D)
    p_plant_electric_net_required_mw = (3, Area.C)
    hfact = (4, Area.P)
    j_tf_coil_full_area = (5, Area.T)
    pflux_fw_neutron_max_mw = (6, Area.C)
    beamfus0 = (7, Area.P)
    temp_plasma_electron_vol_avg_kev = (9, Area.P)
    boundu__15 = (10, Area.NUM)
    beta_norm_max = (11, Area.P)
    f_c_plasma_bootstrap_max = (12, Area.CD)
    boundu__10 = (13, Area.NUM)
    f_j_tf_wp_critical_max = (14, Area.C)  # TODO is this needed
    rmajor = (16, Area.P)
    b_tf_inboard_max = (17, Area.C, "b_tf_inboard_peak_symmetric")
    eta_cd_norm_hcd_primary_max = (18, Area.C)
    boundl__16 = (19, Area.NUM)
    t_burn_min = (20, Area.C)
    f_t_plant_available = (22, Area.CST)
    p_fusion_total_max_mw = (24, Area.C)
    kappa = (25, Area.P)
    triang = (26, Area.P)
    tbrmin = (27, Area.C)
    b_plasma_toroidal_on_axis = (28, Area.P)
    coreradius = (29, Area.IR)
    f_alpha_energy_confinement_min = (31, Area.C)
    epsvmc = (32, Area.NUM)
    boundu__129 = (38, Area.NUM)
    boundu__131 = (39, Area.NUM)
    boundu__135 = (40, Area.NUM)
    dr_blkt_outboard = (41, Area.B)
    f_nd_impurity_electrons__9 = (42, Area.IR)
    sig_tf_case_max = (44, Area.T)
    temp_tf_superconductor_margin_min = (45, Area.T)
    boundu__152 = (46, Area.NUM)
    n_tf_wp_pancakes = (48, Area.T)
    n_tf_wp_layers = (49, Area.T)
    f_nd_impurity_electrons__13 = (50, Area.IR)
    f_p_div_lower = (51, Area.P)
    rad_fraction_sol = (52, Area.P)
    boundu__157 = (53, Area.NUM)
    b_crit_upper_nbti = (54, Area.T)
    dr_shld_inboard = (55, Area.B)
    p_cryo_plant_electric_max_mw = (56, Area.HT)
    boundl__2 = (57, Area.NUM)
    dr_fw_plasma_gap_inboard = (58, Area.B)
    dr_fw_plasma_gap_outboard = (59, Area.B)
    sig_tf_wp_max = (60, Area.T)
    copperaoh_m2_max = (61, Area.TR)
    coheof = (62, Area.PF)
    dr_cs = (63, Area.B)
    ohhghf = (64, Area.PF)
    n_cycle_min = (65, Area.CS)
    oh_steel_frac = (66, Area.PF)
    t_crack_vertical = (67, Area.CS)
    inlet_temp_liq = (68, Area.FWBS)
    outlet_temp_liq = (69, Area.FWBS)
    blpressure_liq = (70, Area.FWBS)
    n_liq_recirc = (71, Area.FWBS)
    bz_channel_conduct_liq = (72, Area.FWBS)
    pnuc_fw_ratio_dcll = (73, Area.FWBS)
    f_nuc_pow_bz_struct = (74, Area.FWBS)
    dx_fw_module = (75, Area.FWBS)
    eta_turbine = (76, Area.HT)
    startupratio = (77, Area.CST)
    fkind = (78, Area.CST)
    eta_ecrh_injector_wall_plug = (79, Area.CD)
    fcoolcp = (80, Area.T)
    n_tf_coil_turns = (81, Area.T)


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

        Raises
        ------
        ProcessValueError
            isweep value greater than IPNSCNS
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
        if self.solver == "fsolve":
            process_output.ocmmnt(
                constants.NOUT, "PROCESS has performed an fsolve (evaluation) run."
            )
        else:
            process_output.ocmmnt(
                constants.NOUT, "PROCESS has performed a VMCON (optimisation) run."
            )
        if ifail != 1:
            process_output.ovarre(constants.NOUT, "Error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.IOTTY, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            process_output.oblnkl(constants.IOTTY)

            logger.critical("Solver returns with ifail /= 1. %s", ifail)

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
            process_output.ovarre(constants.NOUT, "Error flag", "(ifail)", ifail)

            if self.data.numerics.sqsumsq >= 1.0e-2:
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

                logger.warning(
                    f"High final constraint residues. {self.data.numerics.sqsumsq=}"
                )

        process_output.ovarre(
            constants.NOUT,
            "Number of iteration variables",
            "(nvar)",
            self.data.numerics.nvar,
        )
        process_output.ovarre(
            constants.NOUT,
            "Number of constraints (total)",
            "(neqns+nineqns)",
            self.data.numerics.neqns + self.data.numerics.nineqns,
        )
        process_output.ovarre(
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
            process_output.ovarre(
                constants.NOUT,
                "Figure of merit switch",
                "(minmax)",
                self.data.numerics.minmax,
            )

            objf_name = f'"{FiguresOfMerit(abs(self.data.numerics.minmax)).description}"'

            self.data.numerics.objf_name = objf_name

            process_output.ovarre(
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
            process_output.ovarre(
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
                f"{string1} the optimisation parameters to {string2} "
                f"the objective function: {objf_name}\n",
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
                            "\n as shown by the following optimisation parameters"
                            " that are"
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
                "Negative inequality constraint (normalised) residuals "
                "indicate a constraint is satisfied.",
            )
            if self.solver == "fsolve":
                process_output.osubhd(
                    constants.NOUT,
                    "This MFile was produced via an evaluation, not an optimisation, "
                    "and so the constraints might be violated.",
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
                "The code may be stuck in a minimum in the residual space that is "
                "significantly above zero.",
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
                "The code may be stuck in a minimum in the residual space that is "
                "significantly above zero.",
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
        print("Scan Convergence Summary \n")
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
        for iscan in range(1, self.data.scan.isweep + 1):
            if scan_1d_ifail_dict[iscan] == 1:
                converged_count += 1
                print(
                    f"Scan {iscan:02d}: {nsweep_var.fname} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[32mCONVERGED \u001b[0m"
                )
            else:
                print(
                    f"Scan {iscan:02d}: {nsweep_var.fname} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[31mUNCONVERGED \u001b[0m"
                )
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

        print("Scan Convergence Summary\n")
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
                if scan_2d_ifail_list[iscan_1][iscan_2] == 1:
                    converged_count += 1
                    print(
                        (
                            f"Scan {scan_point:02d}: ({nsweep_var.fname} = "
                            f"{sweep_1_values[iscan_1 - 1]}, {nsweep_2_var.fname} "
                            f"= {sweep_2_values[iscan_2 - 1]}) "
                        )
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[32mCONVERGED \u001b[0m"
                    )
                    scan_point += 1
                else:
                    print(
                        (
                            f"Scan {scan_point:02d}: ({nsweep_var.fname} = "
                            f"{sweep_1_values[iscan_1 - 1]}, {nsweep_2_var.fname} = "
                            f"{sweep_2_values[iscan_2 - 1]}) "
                        )
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[31mUNCONVERGED \u001b[0m"
                    )
                    scan_point += 1
        converged_percentage = (
            converged_count / (self.data.scan.isweep * self.data.scan.isweep_2) * 100
        )
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    @staticmethod
    def scan_2d_init(scan_data: ScanData):
        """Scan 2d initialisation"""
        process_output.ovarre(
            constants.MFILE,
            "Number of first variable scan points",
            "(isweep)",
            scan_data.isweep,
        )
        process_output.ovarre(
            constants.MFILE,
            "Number of second variable scan points",
            "(isweep_2)",
            scan_data.isweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning first variable number",
            "(nsweep)",
            scan_data.nsweep,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )

    def scan_1d_write_point_header(self, iscan: int):
        """Scan 1d header"""
        self.data.globals.iscan_global = iscan
        sv = self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, iscan)

        self.data.globals.vlabel = sv.fname
        self.data.globals.xlabel = sv.description

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** Scan point {iscan} of {self.data.scan.isweep} : "
            f"{self.data.globals.xlabel}"
            f", {self.data.globals.vlabel} = {self.data.scan.sweep[iscan - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarre(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan} of {self.data.scan.isweep} : "
            f"{self.data.globals.xlabel} , {self.data.globals.vlabel}"
            f" = {self.data.scan.sweep[iscan - 1]}"
        )

    def scan_2d_write_point_header(self, iscan, iscan_1, iscan_2):
        """Scan 2d header"""
        iscan_r = self.data.scan.isweep_2 - iscan_2 + 1 if iscan_1 % 2 == 0 else iscan_2

        # Makes iscan available globally (read-only)
        self.data.globals.iscan_global = iscan
        sv_1 = self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, iscan_1)

        self.data.globals.vlabel = sv_1.fname
        self.data.globals.xlabel = sv_1.data.description

        sv_2 = self.scan_select(self.data.scan.nsweep_2, self.data.scan.sweep_2, iscan_r)

        self.data.globals.vlabel_2 = sv_2.fname
        self.data.globals.xlabel_2 = sv_2.data.description

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** 2D Scan point {iscan} of "
            f"{self.data.scan.isweep * self.data.scan.isweep_2} : "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]} and"
            f" {self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarre(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {self.data.globals.xlabel}, "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]}"
            f" and {self.data.globals.xlabel_2}, "
            f"{self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
        )

        return iscan_r

    @staticmethod
    def scan_1d_write_plot(scan_data: ScanData):
        """Scan 1d plotter"""
        if scan_data.first_call_1d:
            process_output.ovarre(
                constants.MFILE,
                "Number of scan points",
                "(isweep)",
                scan_data.isweep,
            )
            process_output.ovarre(
                constants.MFILE,
                "Scanning variable number",
                "(nsweep)",
                scan_data.nsweep,
            )

            scan_data.first_call_1d = False

    def scan_select(self, nsweep, sweep, iscan) -> ScanVariables:
        """Select a scan"""
        sv = ScanVariables(nsweep)
        sv.set(self.data, sweep[iscan - 1])
        return sv
