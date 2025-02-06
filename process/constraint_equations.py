import numpy as np

from process.exceptions import ProcessValueError
from process.fortran import constraints, numerics


def constraint_eqns(m, ieqn):
    """Routine that formulates the constraint equations
    author: P J Knight (UKAEA)
    author: J Morris (UKAEA)
    if `ieqn` is zero or negative, evaluate all the constraint equations, otherwise
    evaluate only the `ieqn`th equation.
    """

    cc = []
    con = []
    err = []
    symbol = []
    units = []

    # If ieqn is positive, only evaluate the 'ieqn'th constraint residue,
    # otherwise evaluate all m constraint residues
    if ieqn > 0:
        i1 = ieqn - 1
        i2 = ieqn
    else:
        i1 = 0
        i2 = m

    for i in range(i1, i2):
        constraint_equation = getattr(
            constraints, f"constraint_eqn_{numerics.icc[i]:03d}", None
        )

        if constraint_equation is None:
            raise ProcessValueError(
                "Invalid constraint equation number", icc=numerics.icc[i]
            )

        tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units = constraint_equation()

        tmp_symbol = tmp_symbol.decode().strip()
        tmp_units = tmp_units.decode().strip()

        # Issue 505 Reverse the sign so it works as an inequality constraint (cc(i) > 0)
        # This will have no effect if it is used as an equality constraint because it will be squared.
        cc.append(-tmp_cc)
        con.append(tmp_con)
        err.append(tmp_err)
        symbol.append(tmp_symbol)
        units.append(tmp_units)

        if abs(tmp_cc) > 9.99e99 or not np.isfinite(tmp_cc):
            raise ValueError(
                "Invalid number in constraint equation", icc=numerics.icc[i], cc=tmp_cc
            )

    return cc, con, err, symbol, units
