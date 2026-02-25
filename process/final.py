"""Final output at the end of a scan."""

from tabulate import tabulate

import process.constraints as constraints
from process import output as op
from process import process_output as po
from process.core import constants
from process.data_structure import numerics
from process.objectives import objective_function


def finalise(models, ifail: int, non_idempotent_msg: str | None = None):
    """Routine to print out the final point in the scan.

    Writes to OUT.DAT and MFILE.DAT.

    Parameters
    ----------
    models : process.main.Models
        physics and engineering model objects
    ifail : int
        error flag
    non_idempotent_msg : None | str, optional
        warning about non-idempotent variables, defaults to None
    """
    if ifail == 1:
        po.oheadr(constants.NOUT, "Final Feasible Point")
    else:
        po.oheadr(constants.NOUT, "Final UNFEASIBLE Point")

    # Output relevant to no optimisation
    if numerics.ioptimz == -2:
        output_evaluation()

    # Print non-idempotence warning to OUT.DAT only
    if non_idempotent_msg:
        po.oheadr(constants.NOUT, "NON-IDEMPOTENT VARIABLES")
        po.ocmmnt(constants.NOUT, non_idempotent_msg)

    # Write output to OUT.DAT and MFILE.DAT
    op.write(models, constants.NOUT)


def output_evaluation():
    """Write output for an evaluation run of PROCESS"""
    po.oheadr(constants.NOUT, "Numerics")
    po.ocmmnt(constants.NOUT, "PROCESS has performed an evaluation run.")
    po.oblnkl(constants.NOUT)

    # Evaluate objective function
    norm_objf = objective_function(numerics.minmax)
    po.ovarre(constants.MFILE, "Normalised objective function", "(norm_objf)", norm_objf)

    # Print the residuals of the constraint equations

    residual_error, value, residual, symbols, units = constraints.constraint_eqns(
        numerics.neqns + numerics.nineqns, -1
    )

    labels = [
        numerics.lablcc[j]
        for j in [i - 1 for i in numerics.icc[: numerics.neqns + numerics.nineqns]]
    ]
    physical_constraint = [f"{c} {u}" for c, u in zip(value, units, strict=False)]
    physical_residual = [f"{c} {u}" for c, u in zip(residual, units, strict=False)]

    table_data = {
        "Constraint Name": labels,
        "Constraint Type": symbols,
        "Physical constraint": physical_constraint,
        "Constraint residual": physical_residual,
        "Normalised residual": residual_error,
    }

    po.write(constants.NOUT, tabulate(table_data, headers="keys"))

    for i in range(numerics.neqns):
        constraint_id = numerics.icc[i]
        po.ovarre(
            constants.MFILE,
            f"{labels[i]} normalised residue",
            f"(eq_con{constraint_id:03d})",
            residual_error[i],
        )

    for i in range(numerics.nineqns):
        constraint_id = numerics.icc[numerics.neqns + i]
        po.ovarre(
            constants.MFILE,
            f"{labels[numerics.neqns + i]}",
            f"(ineq_con{constraint_id:03d})",
            residual_error[numerics.neqns + i],
        )
