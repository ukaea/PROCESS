"""Final output at the end of a scan."""

from tabulate import tabulate

from process import output as op
from process.fortran import (
    constants,
    constraints,
    numerics,
)
from process.fortran import (
    process_output as po,
)
from process.objectives import objective_function
from process.utilities.f2py_string_patch import f2py_compatible_to_string


def finalise(models, ifail: int, non_idempotent_msg: None | str = None):
    """Routine to print out the final point in the scan.

    Writes to OUT.DAT and MFILE.DAT.

    :param models: physics and engineering model objects
    :type models: process.main.Models
    :param ifail: error flag
    :type ifail: int
    :param non_idempotent_msg: warning about non-idempotent variables, defaults to None
    :type non_idempotent_msg: None | str, optional
    """
    if ifail == 1:
        po.oheadr(constants.nout, "Final Feasible Point")
    else:
        po.oheadr(constants.nout, "Final UNFEASIBLE Point")

    # Output relevant to no optimisation
    if numerics.ioptimz == -2:
        output_once_through()

    # Print non-idempotence warning to OUT.DAT only
    if non_idempotent_msg:
        po.oheadr(constants.nout, "NON-IDEMPOTENT VARIABLES")
        po.ocmmnt(constants.nout, non_idempotent_msg)

    # Write output to OUT.DAT and MFILE.DAT
    op.write(models, constants.nout)


def output_once_through():
    """Write output for a once-through run of PROCESS"""
    po.oheadr(constants.nout, "Numerics")
    po.ocmmnt(constants.nout, "PROCESS has performed a run witout optimisation.")
    po.oblnkl(constants.nout)

    # Evaluate objective function
    norm_objf = objective_function(numerics.minmax)
    po.ovarre(
        constants.mfile, "Normalised objective function", "(norm_objf)", norm_objf
    )

    # Print the residuals of the constraint equations

    residual_error, value, residual, symbols, units = constraints.constraint_eqns(
        numerics.neqns + numerics.nineqns, -1
    )

    labels = [
        f2py_compatible_to_string(i)
        for i in numerics.lablcc[numerics.icc[: numerics.neqns + numerics.nineqns] - 1]
    ]
    units = [f2py_compatible_to_string(i) for i in units]
    physical_constraint = [
        f"{c} {u}" for c, u in zip(value.tolist(), units, strict=False)
    ]
    physical_residual = [
        f"{c} {u}" for c, u in zip(residual.tolist(), units, strict=False)
    ]

    table_data = {
        "Constraint Name": labels,
        "Constraint Type": symbols.tolist(),
        "Physical constraint": physical_constraint,
        "Constraint residual": physical_residual,
        "Normalised residual": residual_error.tolist(),
    }

    po.write(constants.nout, tabulate(table_data, headers="keys"))

    for i in range(numerics.neqns):
        constraint_id = numerics.icc[i]
        po.ovarre(
            constants.mfile,
            f"{labels[i]} normalised residue",
            f"(eq_con{constraint_id:03d})",
            residual_error[i],
        )

    for i in range(numerics.nineqns):
        constraint_id = numerics.icc[numerics.neqns + i]
        po.ovarre(
            constants.mfile,
            f"{labels[numerics.neqns + i]}",
            f"(ineq_con{constraint_id:03d})",
            residual_error[numerics.neqns + i],
        )
