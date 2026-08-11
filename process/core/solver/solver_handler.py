"""Module containing solver handler routines"""

import logging

from tabulate import tabulate

from process.core import constants, process_output
from process.core.solver.evaluators import Evaluators
from process.core.solver.iteration_variables import (
    load_iteration_variables,
    load_scaled_bounds,
)
from process.core.solver.solver import get_solver
from process.data_structure.numerics import (
    FiguresOfMerit,
    PROCESSRunMode,
    SolverOutputCondition,
)

logger = logging.getLogger(__name__)


class SolverHandler:
    """Creates and runs a solver instance.

    This may be an optimiser (e.g. VMCON) or an equation solver (e.g. fsolve).

    Parameters
    ----------
    models : process.main.Models
        physics and engineering model objects
    solver_name : str
        which solver to use, as specified in solver.py
    data: DataStructure
        data structure object for providing objective/constraint data to
        the solver
    """

    def __init__(self, models, solver_name, data):
        self.models = models
        self.solver_name = solver_name
        self.data = data

    def run(self):
        """Run solver and retry if it fails in certain ways."""
        # Initialise iteration variables and bounds in Fortran
        load_iteration_variables(self.data)
        load_scaled_bounds(self.data)

        # Initialise iteration variables and bounds in Python: relies on Fortran
        # iteration variables being defined above
        # Trim maximum size arrays down to actually used size
        n = self.data.numerics.nvar
        x = self.data.numerics.xcm[:n]
        bndl = self.data.numerics.itv_scaled_lower_bounds[:n]
        bndu = self.data.numerics.itv_scaled_upper_bounds[:n]

        # Define total number of constraints and equality constraints
        m = self.data.numerics.neqns + self.data.numerics.nineqns
        meq = self.data.numerics.neqns

        # Evaluators() calculates the objective and constraint functions and
        # their gradients for a given vector x
        evaluators = Evaluators(self.models, self.data, x)

        # Configure solver for problem
        self.solver = get_solver(self.data, self.solver_name)
        self.solver.set_evaluators(evaluators)
        self.solver.set_bounds(bndl, bndu)
        self.solver.set_opt_params(x)
        self.solver.set_constraints(m, meq)
        ifail = self.solver.solve()

        # If VMCON optimisation has failed then try altering value of epsfcn
        if self.solver_name == "vmcon":
            if ifail != SolverOutputCondition.CONVERGED:
                print("Trying again with new epsfcn")
                # epsfcn is only used in evaluators.Evaluators()
                # TODO epsfcn could be set in Evaluators instance now, don't need to
                # set/unset in self.data.numerics module
                self.data.numerics.epsfcn *= 10  # try new larger value
                print("new epsfcn = ", self.data.numerics.epsfcn)

                ifail = self.solver.solve()
                # First solution attempt failed
                # (ifail != SolverOutputCondition.CONVERGED): supply ifail value
                # to next attempt
                self.data.numerics.epsfcn /= 10  # reset value

            if ifail != SolverOutputCondition.CONVERGED:
                print("Trying again with new epsfcn")
                self.data.numerics.epsfcn /= 10  # try new smaller value
                print("new epsfcn = ", self.data.numerics.epsfcn)
                ifail = self.solver.solve()
                self.data.numerics.epsfcn *= 10  # reset value

            # If VMCON has exited with error code 5
            # (ifail = SolverOutputCondition.NO_SOLUTION) try another run using a
            # multiple of the identity matrix as input for the Hessian b(n,n)
            # Only do this if VMCON has not iterated (nviter=1)
            if (
                ifail == SolverOutputCondition.NO_SOLUTION
                and self.data.numerics.nviter < 2
            ):
                print(
                    "VMCON error code = 5 (SolverOutputCondition.NO_SOLUTION). "
                    "Rerunning VMCON with a new initial estimate of the second "
                    "derivative matrix."
                )
                self.solver.set_b(2.0)
                ifail = self.solver.solve()

        self.output()
        return ifail

    def output(self):
        """Store results back in self.data.numerics module.

        Objective function value, solution vector and constraints vector.
        """
        self.data.numerics.norm_objf = self.solver.objf
        # Slicing required due to Fortran arrays being maximum possible, rather
        # than required, size
        self.data.numerics.xcm[: self.solver.x.shape[0]] = self.solver.x
        self.data.numerics.rcm[: self.solver.conf.shape[0]] = self.solver.conf

        self._numerics_output()
        self._optimisation_parameters_output()

    def _numerics_output(self):
        nums = self.data.numerics

        process_output.oheadr(constants.NOUT, "Numerics")
        s_type = (
            "fsolve (evaluation)" if self.solver == "fsolve" else "VMCON (optimisation)"
        )
        process_output.ocmmnt(
            constants.NOUT,
            f"PROCESS has performed a {s_type} run",
        )
        ifail = self.solver.info
        if ifail != SolverOutputCondition.CONVERGED:
            process_output.ovarre(constants.NOUT, "Error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.IOTTY, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            print()

            logger.critical("Solver returns with ifail /= 1. %s", ifail)

            if self.solver_name == "vmcon":
                self.solver.verror()

            process_output.oblnkl(constants.NOUT)
            print()
        else:
            # Solution found
            descr = "consistent" if self.solver == "fsolve" else "feasible"
            process_output.ocmmnt(
                constants.NOUT, f"and found a {descr} set of parameters."
            )
            process_output.oheadr(constants.IOTTY, f"PROCESS found a {descr} solution")
            process_output.oblnkl(constants.NOUT)
            process_output.ovarre(constants.NOUT, "Error flag", "(ifail)", ifail)

            if nums.sqsumsq >= 1.0e-2:
                string = (
                    "WARNING: Constraint residues are HIGH; consider re-running\n"
                    "   with lower values of EPSVMC to confirm convergence...\n"
                    "   (should be able to get down to about 1.0E-8 okay)\n"
                )
                process_output.ocmmnt(constants.NOUT, ("\n" + string))
                print(string)

                logger.warning(f"High final constraint residues. {nums.sqsumsq=}")

        for d, var, v in (
            ("Number of iteration variables", "(nvar)", nums.nvar),
            (
                "Number of constraints (total)",
                "(neqns+nineqns)",
                nums.neqns + nums.nineqns,
            ),
            ("Optimisation switch", "(ioptimz)", nums.ioptimz),
        ):
            process_output.ovarre(constants.NOUT, d, var, v)

        process_output.ocmmnt(
            constants.NOUT,
            f"     {PROCESSRunMode(nums.ioptimz).description}",
        )

        # Objective function output: none for fsolve
        if self.solver_name != "fsolve":
            process_output.ovarre(
                constants.NOUT,
                "Figure of merit switch",
                "(minmax)",
                nums.minmax,
            )

            nums.objf_name = f'"{FiguresOfMerit(abs(nums.minmax)).description}"'

            for d, var, v, o in (
                ("Objective function name", "(objf_name)", nums.objf_name, ""),
                ("Normalised objective function", "(norm_objf)", nums.norm_objf, "OP "),
                (
                    "VMCON convergence parameter",
                    "(convergence_parameter)",
                    self.data.globals.convergence_parameter,
                    "OP ",
                ),
                (
                    "Number of optimising solver iterations",
                    "(nviter)",
                    nums.nviter,
                    "OP ",
                ),
            ):
                process_output.ovarre(constants.NOUT, d, var, v, o)

        process_output.ovarre(
            constants.NOUT,
            "Square root of the sum of squares of the constraint residuals",
            "(sqsumsq)",
            nums.sqsumsq,
            "OP ",
        )

        process_output.oblnkl(constants.NOUT)

        if self.solver_name == "fsolve":
            process_output.write(
                constants.NOUT,
                "PROCESS has solved using fsolve.\n"
                if ifail == SolverOutputCondition.CONVERGED
                else "PROCESS failed to solve using fsolve.\n",
            )
        else:
            process_output.write(
                constants.NOUT,
                (
                    (
                        "PROCESS has successfully optimised"
                        if ifail == SolverOutputCondition.CONVERGED
                        else "PROCESS has failed to optimise"
                    )
                    + " the optimisation parameters to"
                    + ("minimise" if nums.minmax > 0 else "maximise")
                    + f" the objective function: {nums.objf_name}\n"
                ),
            )

    def _optimisation_parameters_output(self):
        nums = self.data.numerics

        written_warning = False

        # Output optimisation parameters
        solution_vector_table = []
        for i in range(nums.nvar):
            nums.xcs[i] = nums.xcm[i] * nums.scafc[i]

            name = nums.lablxc[nums.ixc[i] - 1]
            solution_vector_table.append([name, nums.xcs[i], nums.xcm[i]])

            xminn = 1.01 * nums.itv_scaled_lower_bounds[i]
            xmaxx = 0.99 * nums.itv_scaled_upper_bounds[i]

            # Write to output file if close to optimisation parameter bounds
            if nums.xcm[i] < xminn or nums.xcm[i] > xmaxx:
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

                xcval = nums.xcm[i] * nums.scafc[i]

                if nums.xcm[i] < xminn:
                    location, bound = "below", "lower"
                    bounds = nums.itv_scaled_lower_bounds
                else:
                    location, bound = "above", "upper"
                    bounds = nums.itv_scaled_upper_bounds
                process_output.write(
                    constants.NOUT,
                    f"   {name:<30}= {xcval} is at or {location} its {bound} bound:"
                    f" {bounds[i] * nums.scafc[i]}",
                )

            if nums.boundu[i] == nums.boundl[i]:
                xnorm = 1.0
            else:
                xnorm = min(
                    max(
                        (nums.xcm[i] - nums.itv_scaled_lower_bounds[i])
                        / (
                            nums.itv_scaled_upper_bounds[i]
                            - nums.itv_scaled_lower_bounds[i]
                        ),
                        0.0,
                    ),
                    1.0,
                )

            # Write optimisation parameters to mfile
            for d, var, v in (
                (nums.lablxc[nums.ixc[i] - 1], f"(itvar{i + 1:03d})", nums.xcs[i]),
                (
                    f"{name} (final value/initial value)",
                    f"(xcm{i + 1:03d})",
                    nums.xcm[i],
                ),
                (f"{name} (range normalised)", f"(nitvar{i + 1:03d})", xnorm),
                (
                    f"{name} (upper bound)",
                    f"(boundu{i + 1:03d})",
                    nums.itv_scaled_upper_bounds[i] * nums.scafc[i],
                ),
                (
                    f"{name} (lower bound)",
                    f"(boundl{i + 1:03d})",
                    nums.itv_scaled_lower_bounds[i] * nums.scafc[i],
                ),
            ):
                process_output.ovarre(constants.MFILE, d, var, v)

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
