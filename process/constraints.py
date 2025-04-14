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


@ConstraintManager.register_constraint(2, "MW/m3", "=")
def constraint_equation_2():
    """author: J. Morris

     i_rad_loss: switch for radiation loss term usage in power balance (see User Guide):
    -  0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
    -  1 total power lost is scaling power plus core radiation only
    -  2 total power lost is scaling power only, with no additional
        allowance for radiation. This is not recommended for power plant models.

    i_plasma_ignited: switch for ignition assumption:
    -  0 do not assume plasma ignition;
    -  1 assume ignited (but include auxiliary power in costs)

     pden_electron_transport_loss_mw: electron transport power per volume (MW/m3)
     pden_ion_transport_loss_mw:  ion transport power per volume (MW/m3)
     pden_plasma_rad_mw: total radiation power per volume (MW/m3)
     pden_plasma_core_rad_mw: total core radiation power per volume (MW/m3)
     f_alpha_plasma: fraction of alpha power deposited in plasma
     alpha_power_density_total: alpha power per volume (MW/m3)
     charged_power_density: non-alpha charged particle fusion power per volume (MW/m3)
     pden_plasma_ohmic_mw: ohmic heating power per volume (MW/m3)
     p_hcd_injected_total_mw: total auxiliary injected power (MW)
     vol_plasma: plasma volume (m3)
    """
    # pscaling: total transport power per volume (MW/m3)

    pscaling = (
        fortran.physics_variables.pden_electron_transport_loss_mw
        + fortran.physics_variables.pden_ion_transport_loss_mw
    )
    # Total power lost is scaling power plus radiation:
    if fortran.physics_variables.i_rad_loss == 0:
        pnumerator = pscaling + fortran.physics_variables.pden_plasma_rad_mw
    elif fortran.physics_variables.i_rad_loss == 1:
        pnumerator = pscaling + fortran.physics_variables.pden_plasma_core_rad_mw
    else:
        pnumerator = pscaling

    # if plasma not ignited include injected power
    if fortran.physics_variables.i_plasma_ignited == 0:
        pdenom = (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_density_total
            + fortran.physics_variables.charged_power_density
            + fortran.physics_variables.pden_plasma_ohmic_mw
            + fortran.current_drive_variables.p_hcd_injected_total_mw
            / fortran.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_density_total
            + fortran.physics_variables.charged_power_density
            + fortran.physics_variables.pden_plasma_ohmic_mw
        )

    cc = 1.0 - pnumerator / pdenom

    return ConstraintResult(cc, pdenom * (1.0 - cc), pdenom * cc)


@ConstraintManager.register_constraint(3, "MW/m3", "=")
def constraint_equation_3():
    """Global power balance equation for ions
    i_plasma_ignited: switch for ignition assumption
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)

    pden_ion_transport_loss_mw: ion transport power per volume (MW/m3)
    piepv: ion/electron equilibration power per volume (MW/m3)
    f_alpha_plasma: fraction of alpha power deposited in plasma
    alpha_power_ions_density: alpha power per volume to ions (MW/m3)
    p_hcd_injected_ions_mw: auxiliary injected power to ions (MW)
    vol_plasma: plasma volume (m3)
    """
    # No assume plasma ignition:
    if fortran.physics_variables.i_plasma_ignited == 0:
        cc = 1.0 - (
            fortran.physics_variables.pden_ion_transport_loss_mw
            + fortran.physics_variables.piepv
        ) / (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_ions_density
            + fortran.current_drive_variables.p_hcd_injected_ions_mw
            / fortran.physics_variables.vol_plasma
        )
        return ConstraintResult(
            cc,
            (
                fortran.physics_variables.f_alpha_plasma
                * fortran.physics_variables.alpha_power_ions_density
                + fortran.current_drive_variables.p_hcd_injected_ions_mw
                / fortran.physics_variables.vol_plasma
            )
            * (1.0 - cc),
            (
                fortran.physics_variables.f_alpha_plasma
                * fortran.physics_variables.alpha_power_ions_density
                + fortran.current_drive_variables.p_hcd_injected_ions_mw
                / fortran.physics_variables.vol_plasma
            )
            * cc,
        )

    # Plasma ignited:
    cc = 1.0 - (
        fortran.physics_variables.pden_ion_transport_loss_mw
        + fortran.physics_variables.piepv
    ) / (
        fortran.physics_variables.f_alpha_plasma
        * fortran.physics_variables.alpha_power_ions_density
    )
    return ConstraintResult(
        cc,
        (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_ions_density
        )
        * (1.0 - cc),
        (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_ions_density
        )
        * cc,
    )


@ConstraintManager.register_constraint(4, "MW/m3", "=")
def constraint_equation_4():
    """Global power balance equation for electrons
    author: P B Lloyd, CCFE, Culham Science Centre

    i_rad_loss: switch for radiation loss term usage in power balance
    - 0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
    - 1 total power lost is scaling power plus core radiation only
    - 2 total power lost is scaling power only, with no additional
        allowance for radiation. This is not recommended for power plant models.

    i_plasma_ignited: switch for ignition assumption:
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)

    pden_electron_transport_loss_mw: electron transport power per volume (MW/m3)
    pden_plasma_rad_mw: total radiation power per volume (MW/m3)
    pden_plasma_core_rad_mw: total core radiation power per volume (MW/m3)
    f_alpha_plasma: fraction of alpha power deposited in plasma
    alpha_power_electron_density: alpha power per volume to electrons (MW/m3)
    piepv: ion/electron equilibration power per volume (MW/m3)
    p_hcd_injected_electrons_mw: auxiliary injected power to electrons (MW)
    vol_plasma: plasma volume (m3)
    """
    # pscaling: total transport power per volume (MW/m3)

    pscaling = fortran.physics_variables.pden_electron_transport_loss_mw
    # Total power lost is scaling power plus radiation:
    if fortran.physics_variables.i_rad_loss == 0:
        pnumerator = pscaling + fortran.physics_variables.pden_plasma_rad_mw
    elif fortran.physics_variables.i_rad_loss == 1:
        pnumerator = pscaling + fortran.physics_variables.pden_plasma_core_rad_mw
    else:
        pnumerator = pscaling

    # if plasma not ignited include injected power
    if fortran.physics_variables.i_plasma_ignited == 0:
        pdenom = (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_electron_density
            + fortran.physics_variables.piepv
            + fortran.current_drive_variables.p_hcd_injected_electrons_mw
            / fortran.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            fortran.physics_variables.f_alpha_plasma
            * fortran.physics_variables.alpha_power_electron_density
            + fortran.physics_variables.piepv
        )

    cc = 1.0 - pnumerator / pdenom
    return ConstraintResult(cc, pdenom * (1.0 - cc), pdenom * cc)


@ConstraintManager.register_constraint(5, "/m3", "<=")
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


@ConstraintManager.register_constraint(7, "/m3", "=")
def constraint_equation_7():
    """Equation for hot beam ion density

    i_plasma_ignited: switch for ignition assumption:
    - 0 do not assume plasma ignition
    - 1 assume ignited (but include auxiliary power in costs)
    O
    bviously, i_plasma_ignited must be zero if current drive is required.
    If i_plasma_ignited=1, any auxiliary power is assumed to be used only
    during plasma start-up, and is excluded from all steady-state
    power balance calculations.
    beam_density_out: hot beam ion density from calculation (/m3)
    nd_beam_ions: hot beam ion density, variable (/m3)
    """
    if fortran.physics_variables.i_plasma_ignited == 1:
        raise ProcessValueError(
            "Do not use constraint equation 7 if i_plasma_ignited=1"
        )

    cc = (
        1.0
        - fortran.physics_variables.beam_density_out
        / fortran.physics_variables.nd_beam_ions
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.nd_beam_ions * (1.0 - cc),
        fortran.physics_variables.nd_beam_ions * cc,
    )


@ConstraintManager.register_constraint(11, "m", "=")
def constraint_equation_11():
    """Equation for radial build
    author: P B Lloyd, CCFE, Culham Science Centre

    rbld: sum of thicknesses to the major radius (m)
    rmajor: plasma major radius (m)
    """
    cc = 1.0 - fortran.build_variables.rbld / fortran.physics_variables.rmajor
    return ConstraintResult(
        cc,
        fortran.physics_variables.rmajor * (1.0 - cc),
        fortran.physics_variables.rmajor * cc,
    )


@ConstraintManager.register_constraint(12, "V.sec", ">=")
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


@ConstraintManager.register_constraint(13, "sec", ">=")
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


@ConstraintManager.register_constraint(14, "", "=")
def constraint_equation_14():
    """Equation to fix number of NBI decay lengths to plasma centre
    author: P B Lloyd, CCFE, Culham Science Centre

    n_beam_decay_lengths_core: neutral beam e-decay lengths to plasma centre
    tbeamin: permitted neutral beam e-decay lengths to plasma centre
    """
    cc = (
        1.0
        - fortran.current_drive_variables.n_beam_decay_lengths_core
        / fortran.current_drive_variables.tbeamin
    )
    return ConstraintResult(
        cc,
        fortran.current_drive_variables.tbeamin * (1.0 - cc),
        fortran.current_drive_variables.tbeamin * cc,
    )


def constraint_equation_29():
    """Equation for inboard major radius: This is a consistency equation
    author: P B Lloyd, CCFE, Culham Science Centre

    rmajor: plasma major radius (m) (iteration variable 3)
    rminor: plasma minor radius (m)
    rinboard: plasma inboard radius (m)
    """
    cc = (
        1.0
        - (fortran.physics_variables.rmajor - fortran.physics_variables.rminor)
        / fortran.build_variables.rinboard
    )
    return ConstraintResult(
        cc,
        fortran.build_variables.rinboard * (1.0 - cc),
        fortran.build_variables.rinboard * cc,
    )


@ConstraintManager.register_constraint(43, "deg C", "=")
def constraint_equation_43():
    """Equation for average centrepost temperature: This is a consistency equation (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

    temp_cp_average: average temp of TF coil inboard leg conductor (C)e
    tcpav2: centrepost average temperature (C) (for consistency)
    itart: switch for spherical tokamak (ST) models:
    - 0 use conventional aspect ratio models;
    - 1 use spherical tokamak models
    """
    if fortran.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 43 if itart=0")

    if fortran.tfcoil_variables.i_tf_sup == 0:
        temp_cp_average = fortran.tfcoil_variables.temp_cp_average - 273.15
        tcpav2 = fortran.tfcoil_variables.tcpav2 - 273.15
    else:
        temp_cp_average = fortran.tfcoil_variables.temp_cp_average
        tcpav2 = fortran.tfcoil_variables.tcpav2

    cc = 1.0 - temp_cp_average / tcpav2

    return ConstraintResult(cc, tcpav2 * (1.0 - cc), tcpav2 * cc)


@ConstraintManager.register_constraint(51, "V.s", "=")
def constraint_equation_51():
    """Equation to enforce startup flux = available startup flux
    author: P B Lloyd, CCFE, Culham Science Centre

    vs_plasma_res_ramp: resistive losses in startup V-s (Wb)
    vs_plasma_ind_ramp:  internal and external plasma inductance V-s (Wb))
    vs_cs_pf_total_ramp:  total flux swing for startup (Wb)
    """
    cc = 1.0 - fortran.pfcoil_variables.fvs_cs_pf_total_ramp * abs(
        (
            fortran.physics_variables.vs_plasma_res_ramp
            + fortran.physics_variables.vs_plasma_ind_ramp
        )
        / fortran.pfcoil_variables.vs_cs_pf_total_ramp
    )
    return ConstraintResult(
        cc,
        fortran.pfcoil_variables.vs_cs_pf_total_ramp * (1.0 - cc),
        fortran.pfcoil_variables.vs_cs_pf_total_ramp * cc,
    )


@ConstraintManager.register_constraint(85, "years", "=")
def constraint_equation_85():
    """Equality constraint for the centerpost (CP) lifetime
    Author : S Kahn

    Depending on the chosen option i_cp_lifetime:
    - 0 : The CP full power year lifelime is set by the user (cplife_input)
    - 1 : The CP lifelime is equal to the divertor one
    - 2 : The CP lifetime is equal to the breeding blankets one
    - 3 : The CP lifetime is equal to the plant one

    cplife: calculated CP full power year lifetime (years)
    life_blkt_fpy: calculated first wall/blanket power year lifetime (years)
    divlife: calculated divertor  power year lifetime (years)
    i_cp_lifetime: switch chosing which plant element the CP
        the CP lifetime must equate
    """
    # The CP lifetime is equal to the the divertor one
    if fortran.cost_variables.i_cp_lifetime == 0:
        cc = 1.0 - fortran.cost_variables.cplife / fortran.cost_variables.cplife_input

    elif fortran.cost_variables.i_cp_lifetime == 1:
        cc = 1.0 - fortran.cost_variables.cplife / fortran.cost_variables.divlife

    # The CP lifetime is equal to the tritium breeding blankets / FW one
    elif fortran.cost_variables.i_cp_lifetime == 2:
        cc = 1.0 - fortran.cost_variables.cplife / fortran.fwbs_variables.life_blkt_fpy

    elif fortran.cost_variables.i_cp_lifetime == 3:
        cc = 1.0 - fortran.cost_variables.cplife / fortran.cost_variables.tlife

    return ConstraintResult(
        cc,
        fortran.cost_variables.divlife * (1.0 - cc),
        fortran.cost_variables.divlife * cc,
    )


@ConstraintManager.register_constraint(92, "", "=")
def constraint_equation_92():
    """Equation for checking is D/T ratio is consistent, and sums to 1.
    author: G Turkington, UKAEA

    f_deuterium: fraction of deuterium ions
    f_tritium: fraction of tritium ions
    f_helium3: fraction of helium-3 ions
    """
    f_deuterium = 1.0 - (
        fortran.physics_variables.f_tritium + fortran.physics_variables.f_helium3
    )
    cc = 1.0 - (
        f_deuterium
        + fortran.physics_variables.f_tritium
        + fortran.physics_variables.f_helium3
    )

    return ConstraintResult(cc, 1.0, cc)


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
