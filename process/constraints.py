from collections.abc import Callable, Hashable
from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np

import process.fortran as fortran
from process.exceptions import ProcessError, ProcessValueError

ConstraintSymbolType = Literal["=", ">=", "<="]


@dataclass
class ConstraintResult:
    """The constraint quantities given the current state of the code
    (aka given an evaluation at the point x).
    """

    cc: float
    """The normalised residual of the constraint."""
    con: float
    """The value of the constraint (in the physical units)."""
    err: float
    """The residual error of the constraint (in the physical units)."""


@dataclass
class ConstraintRegistration:
    """Contains information about a constraint.

    E.g. how to call the constraint, its units, and its symbol (=, >=, <=)
    """

    name: Hashable
    result: Callable[[], ConstraintResult]
    units: str
    symbol: ConstraintSymbolType


class ConstraintManager:
    _constraint_registry: ClassVar[dict[Hashable, ConstraintRegistration]] = {}
    """An internal registry of the PROCESS constraint equations"""

    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__} cannot be instantiated.")

    @classmethod
    def register_constraint(
        cls, name: Hashable, units: str, symbol: ConstraintSymbolType
    ) -> Callable[[], Callable[[], ConstraintResult]]:
        def wrapper(wrapped_func: Callable[[], ConstraintResult]):
            if name in cls._constraint_registry:
                raise ValueError(f"Constraint {name} already exists.")
            cls._constraint_registry[name] = ConstraintRegistration(
                name, wrapped_func, units, symbol
            )

            return wrapped_func

        return wrapper

    @classmethod
    def get_constraint(cls, name: Hashable):
        return cls._constraint_registry.get(name)


@ConstraintManager.register_constraint(1, "", "=")
def constraint_equation_1():
    """Relationship between beta, temperature (keV) and density

    author: J Morris

    beta: total plasma beta
    beta_{ft}: fast alpha beta component
    beta_{NBI}: neutral beam beta component
    n_e: electron density [/m3]
    n_i: total ion density [/m3]
    T_e: density weighted average electron temperature [keV]
    T_i: density weighted average ion temperature [keV]
    B_{tot}: total toroidal + poloidal field [T]
    """
    cc = (
        1.0
        - (
            fortran.physics_variables.beta_fast_alpha
            + fortran.physics_variables.beta_beam
            + 2.0e3
            * fortran.constants.rmu0
            * fortran.constants.electron_charge
            * (
                fortran.physics_variables.dene * fortran.physics_variables.ten
                + fortran.physics_variables.nd_ions_total
                * fortran.physics_variables.tin
            )
            / fortran.physics_variables.btot**2
        )
        / fortran.physics_variables.beta
    )
    return ConstraintResult(
        cc=cc,
        con=(fortran.physics_variables.beta * (1.0 - cc)),
        err=(fortran.physics_variables.beta * cc),
    )


@ConstraintManager.register_constraint(5, "/m3", "<")
def constraint_equation_5():
    """Equation for density upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fdene: f-value for density limit
    dene: electron density (/m3)
    dnelimt: density limit (/m3)
    dnla: line averaged electron density (m-3)

    i_density_limit:
    - 1 old ASDEX;
    - 2 Borrass model for ITER (I);
    - 3 Borrass model for ITER (II);
    - 4 JET edge radiation;
    - 5 JET simplified;
    - 6 Hugill-Murakami Mq limit;
    - 7 Greenwald limit
    """
    # Apply Greenwald limit to line-averaged density
    if fortran.physics_variables.i_density_limit == 7:
        return ConstraintResult(
            fortran.physics_variables.dnla / fortran.physics_variables.dnelimt
            - 1.0 * fortran.constraint_variables.fdene,
            fortran.constraint_variables.fdene * fortran.physics_variables.dnelimt,
            fortran.constraint_variables.fdene * fortran.physics_variables.dnelimt
            - fortran.physics_variables.dnla,
        )

    cc = (
        fortran.physics_variables.dene / fortran.physics_variables.dnelimt
        - 1.0 * fortran.constraint_variables.fdene
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.dnelimt * (1.0 - cc),
        fortran.physics_variables.dene * cc,
    )


@ConstraintManager.register_constraint(12, "V.sec", ">")
def constraint_equation_12():
    """Equation for volt-second capability lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    vs_plasma_total_required: total V-s needed (Wb)
    vs_plasma_total_required (lower limit) is positive; vs_cs_pf_total_pulse (available) is negative
    fvs: f-value for flux-swing (V-s) requirement (STEADY STATE)
    vs_cs_pf_total_pulse: total flux swing for pulse (Wb)
    """
    # vs_cs_pf_total_pulse is negative, requires sign change
    cc = (
        1.0
        - fortran.constraint_variables.fvs
        * (-fortran.pfcoil_variables.vs_cs_pf_total_pulse)
        / fortran.physics_variables.vs_plasma_total_required
    )

    return ConstraintResult(
        cc,
        fortran.pfcoil_variables.vs_plasma_total_required * (1.0 - cc),
        fortran.pfcoil_variables.vs_plasma_total_required * cc,
    )


@ConstraintManager.register_constraint(13, "sec", ">")
def constraint_equation_13():
    """Equation for burn time lower limit

    author: P B Lloyd, CCFE, Culham Science Centre

    ft_burn: f-value for minimum burn time
    t_burn: burn time (s) (calculated if i_pulsed_plant=1)
    t_burn_min:  minimum burn time (s)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.ft_burn
        * fortran.times_variables.t_burn
        / fortran.constraint_variables.t_burn_min,
        fortran.constraint_variables.t_burn_min / fortran.constraint_variables.ft_burn,
        fortran.constraint_variables.t_burn_min / fortran.constraint_variables.ft_burn
        - fortran.times_variables.t_burn,
    )


def constraint_eqns(m: int, ieqn: int):
    """Evaluates the constraints given the current state of PROCESS.

    :param m: The number of constraints to evaluate
    :param ieqn: Evaluates the 'ieqn'th constraint equation (index starts at 1)
    or all equations if <= 0
    """

    if ieqn > 0:
        i1 = ieqn - 1
        i2 = ieqn
    else:
        i1 = 0
        i2 = m

    cc, con, err, symbol, units = [], [], [], [], []

    for i in range(i1, i2):
        constraint_id = fortran.numerics.icc[i].item()
        try:
            constraint = ConstraintManager.get_constraint(constraint_id)

            if constraint is None:
                tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units = getattr(
                    fortran.constraints, f"constraint_eqn_{constraint_id:03d}"
                )()
            else:
                result = constraint.result()
                tmp_cc, tmp_con, tmp_err = result.cc, result.con, result.err
                tmp_symbol, tmp_units = constraint.symbol, constraint.units

        except AttributeError as e:
            error_msg = f"Constraint equation {i + 1} cannot be found"
            raise ProcessError(error_msg) from e

        if np.isnan(tmp_cc) or np.isinf(tmp_cc) or abs(tmp_cc) > 9.99e99:
            error_msg = (
                f"Constraint equation {constraint_id} returned an invalid residual"
            )

            raise ProcessValueError(error_msg, cc=tmp_cc)

        # Reverse the sign so it works as an inequality constraint (cc(i) > 0)
        # This will have no effect if it is used as an equality constraint because it will be squared.
        cc.append(-tmp_cc)
        con.append(tmp_con)
        err.append(tmp_err)
        symbol.append(tmp_symbol)
        units.append(tmp_units)

    return cc, con, err, symbol, units
