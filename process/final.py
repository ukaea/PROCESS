"""Final output at the end of a scan."""

from process import fortran as ft
from process.fortran import final_module as fm
from process import output as op
from process.fortran import process_output as po


def finalise(models, ifail, non_idempotent_msg: None | str = None):
    """Routine to print out the final point in the scan.

    Writes to OUT.DAT and MFILE.DAT.

    :param models: physics and engineering model objects
    :type models: process.main.Models
    :param ifail: error flag
    :type ifail: int
    :param non_idempotent_msg: warning about non-idempotent variables, defaults to None
    :type non_idempotent_msg: None | str, optional
    """
    fm.final_header(ifail)

    # Output relevant to no optimisation
    if ft.numerics.ioptimz == -2:
        fm.no_optimisation()

    # Print non-idempotence warning to OUT.DAT only
    if non_idempotent_msg:
        po.oheadr(ft.constants.nout, "NON-IDEMPOTENT VARIABLES")
        po.ocmmnt(ft.constants.nout, non_idempotent_msg)

    # Write output to OUT.DAT and MFILE.DAT
    op.write(models, ft.constants.nout)

    fm.final_output()
