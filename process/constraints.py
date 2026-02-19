import logging
from collections.abc import Callable, Hashable
from dataclasses import asdict, dataclass
from typing import ClassVar, Literal

import numpy as np

import process.data_structure as data_structure
from process import constants
from process.exceptions import ProcessError, ProcessValueError

ConstraintSymbolType = Literal["=", ">=", "<="]

logger = logging.getLogger(__name__)


@dataclass
class ConstraintRegistration:
    """Contains the constraint equation and metadata about the constraint."""

    name: Hashable
    """The name (often a number) of the constraint. It can be any hashable e.g. a string."""
    constraint_equation: Callable[["ConstraintRegistration"], "ConstraintResult"]
    """The constraint equation that, when called, returns the normalised resiudal,
    constraint value, and constraint error.
    """
    units: str
    """The units of the constraint value and error."""
    symbol: ConstraintSymbolType
    """The type of constraint (<=, >=, ==). Only used for writing output diagnostics,
    this does not impact the calculations.
    """


@dataclass
class ConstraintResult(ConstraintRegistration):
    """The constraint quantities given the current state of the code
    (aka given the current iteration variables).
    """

    residual: float
    """The residual of the constraint. Ie the amount a constraint is or is not violated.
    A negative residual indicates infeasibility.
    A positive residual indicates feasibility.
    A 0 residual indicates that a constraint value is exactly equal to its bound (feasible, just).

    The residual will have the same physical units as the constraint value/bound.
    """
    normalised_residual: float
    """The normalised residual of the constraint. The sign of the normalised residual is interpreted
    identical to the residual.
    """
    constraint_value: float
    """The value of the constraint (in the physical units)."""
    constraint_bound: float
    """The bound of the constraint (in the physical units)."""


class ConstraintManager:
    """A singleton class that manages the registration of constraint equations
    and metadata.

    This class maintains an internal registry of constraints indexed by their names.
    Classmethods provide access to this registry or to directly evaluate constraints.
    """

    _constraint_registry: ClassVar[dict[Hashable, ConstraintRegistration]] = {}
    """An internal registry of the PROCESS constraint equations"""

    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__} cannot be instantiated.")

    @classmethod
    def num_constraints(cls):
        """Return the number of constraints currently in the registry"""
        return len(cls._constraint_registry)

    @classmethod
    def register_constraint(
        cls, name: Hashable, units: str, symbol: ConstraintSymbolType
    ) -> Callable[[], Callable[[], ConstraintResult]]:
        """A decorator to add a constraint equation with metadata to the registry.

        The decorator should wrap a function with no argument which returns a
        ConstraintResult.

        Parameters
        ----------
        name : Hashable
            the name of the constraint and how it can be indexed from the registry
        units : str
            the units of the constraint written to the output files
        symbol : str
            the symbol of the constraint written to the output files

        """

        def wrapper(func: Callable[[], ConstraintResult]):
            if name in cls._constraint_registry:
                raise ValueError(f"Constraint {name} already exists.")
            cls._constraint_registry[name] = ConstraintRegistration(
                name, func, units, symbol
            )

            return func

        return wrapper

    @classmethod
    def get_constraint(cls, name: Hashable):
        """Retrieves a constraint registration from the registry given its name.
        Returns None if no constraint with the name exists.

        Parameters
        ----------
        name : Hashable
            the name of the constraint

        Returns
        -------
        ConstraintRegistration | None
            the constraint registration object
        """
        return cls._constraint_registry.get(name)

    @classmethod
    def evaluate_constraint(cls, name: Hashable):
        """Evalutes a constraint with a given name.

        Parameters
        ----------
        name : Hashable
            the name of the constraint

        Returns
        -------
        ConstraintResult | None
            the result of evaluating the constraint
        """
        registration = cls.get_constraint(name)

        if registration is None:
            error_msg = f"Constraint '{name}' cannot be found."
            raise ProcessError(error_msg)

        return registration.constraint_equation(registration)


def leq(value: float, bound: float, registration: ConstraintRegistration):
    """The equation `value <= bound`."""
    residual = value - bound
    normalised_residual = (value / bound) - 1.0
    return ConstraintResult(
        **asdict(registration),
        constraint_value=value,
        constraint_bound=bound,
        normalised_residual=normalised_residual,
        residual=residual,
    )


def geq(value: float, bound: float, registration: ConstraintRegistration):
    """The equation `value >= bound`."""
    residual = bound - value
    normalised_residual = 1.0 - (value / bound)
    return ConstraintResult(
        **asdict(registration),
        constraint_value=value,
        constraint_bound=bound,
        normalised_residual=normalised_residual,
        residual=residual,
    )


def eq(value: float, bound: float, registration: ConstraintRegistration):
    """The equation `value = bound`."""
    residual = value - bound
    normalised_residual = 1.0 - (value / bound)
    return ConstraintResult(
        **asdict(registration),
        constraint_value=value,
        constraint_bound=bound,
        normalised_residual=normalised_residual,
        residual=residual,
    )


@ConstraintManager.register_constraint(1, "", "=")
def constraint_equation_1(constraint_registration):
    """Relationship between beta, temperature (keV) and density

    beta_total_vol_avg: total plasma beta
    beta_{ft}: fast alpha beta component
    beta_{NBI}: neutral beam beta component
    n_e: electron density [/m3]
    n_i: total ion density [/m3]
    T_e: density weighted average electron temperature [keV]
    T_i: density weighted average ion temperature [keV]
    B_{tot}: total toroidal + poloidal field [T]
    """
    return eq(
        # Density weighted temperature is used here as 〈nT〉 != 〈n〉_V * 〈T〉_V
        (
            data_structure.physics_variables.beta_fast_alpha
            + data_structure.physics_variables.beta_beam
            + 2.0e3
            * constants.RMU0
            * constants.ELECTRON_CHARGE
            * (
                data_structure.physics_variables.nd_plasma_electrons_vol_avg
                * data_structure.physics_variables.temp_plasma_electron_density_weighted_kev
                + data_structure.physics_variables.nd_plasma_ions_total_vol_avg
                * data_structure.physics_variables.temp_plasma_ion_density_weighted_kev
            )
            / data_structure.physics_variables.b_plasma_total**2
        ),
        data_structure.physics_variables.beta_total_vol_avg,
        constraint_registration,
    )


@ConstraintManager.register_constraint(2, "MW/m3", "=")
def constraint_equation_2(constraint_registration):
    """

     i_rad_loss: switch for radiation loss term usage in power balance (see User Guide):
    -  0 total power lost is scaling power plus radiation (needed for i_plasma_pedestal=2,3)
    -  1 total power lost is scaling power plus core radiation only
    -  2 total power lost is scaling power only, with no additional
        allowance for radiation. This is not recommended for power plant models.

    i_plasma_ignited: switch for ignition assumption:
    -  0 do not assume plasma ignition;
    -  1 assume ignited (but include auxiliary power in costs)

     pden_electron_transport_loss_mw: electron transport power per volume (MW/m3)
     pden_ion_transport_loss_mw: ion transport power per volume (MW/m3)
     pden_plasma_rad_mw: total radiation power per volume (MW/m3)
     pden_plasma_core_rad_mw: total core radiation power per volume (MW/m3)
     f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
     pden_alpha_total_mw: alpha power per volume (MW/m3)
     pden_non_alpha_charged_mw: non-alpha charged particle fusion power per volume (MW/m3)
     pden_plasma_ohmic_mw: ohmic heating power per volume (MW/m3)
     p_hcd_injected_total_mw: total auxiliary injected power (MW)
     vol_plasma: plasma volume (m3)
    """
    # pscaling: total transport power per volume (MW/m3)

    pscaling = (
        data_structure.physics_variables.pden_electron_transport_loss_mw
        + data_structure.physics_variables.pden_ion_transport_loss_mw
    )
    # Total power lost is scaling power plus radiation:
    if data_structure.physics_variables.i_rad_loss == 0:
        pnumerator = pscaling + data_structure.physics_variables.pden_plasma_rad_mw
    elif data_structure.physics_variables.i_rad_loss == 1:
        pnumerator = pscaling + data_structure.physics_variables.pden_plasma_core_rad_mw
    else:
        pnumerator = pscaling

    # if plasma not ignited include injected power
    if data_structure.physics_variables.i_plasma_ignited == 0:
        pdenom = (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.pden_alpha_total_mw
            + data_structure.physics_variables.pden_non_alpha_charged_mw
            + data_structure.physics_variables.pden_plasma_ohmic_mw
            + data_structure.current_drive_variables.p_hcd_injected_total_mw
            / data_structure.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.pden_alpha_total_mw
            + data_structure.physics_variables.pden_non_alpha_charged_mw
            + data_structure.physics_variables.pden_plasma_ohmic_mw
        )

    return eq(pnumerator, pdenom, constraint_registration)


@ConstraintManager.register_constraint(3, "MW/m3", "=")
def constraint_equation_3(constraint_registration):
    """Global power balance equation for ions
    i_plasma_ignited: switch for ignition assumption
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)

    pden_ion_transport_loss_mw: ion transport power per volume (MW/m3)
    pden_ion_electron_equilibration_mw: ion/electron equilibration power per volume (MW/m3)
    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    f_pden_alpha_ions_mw: alpha power per volume to ions (MW/m3)
    p_hcd_injected_ions_mw: auxiliary injected power to ions (MW)
    vol_plasma: plasma volume (m3)
    """
    # No assume plasma ignition:
    if data_structure.physics_variables.i_plasma_ignited == 0:
        return eq(
            (
                data_structure.physics_variables.pden_ion_transport_loss_mw
                + data_structure.physics_variables.pden_ion_electron_equilibration_mw
            ),
            (
                data_structure.physics_variables.f_p_alpha_plasma_deposited
                * data_structure.physics_variables.f_pden_alpha_ions_mw
                + data_structure.current_drive_variables.p_hcd_injected_ions_mw
                / data_structure.physics_variables.vol_plasma
            ),
            constraint_registration,
        )

    # Plasma ignited
    return eq(
        (
            data_structure.physics_variables.pden_ion_transport_loss_mw
            + data_structure.physics_variables.pden_ion_electron_equilibration_mw
        ),
        (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_ions_mw
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(4, "MW/m3", "=")
def constraint_equation_4(constraint_registration):
    """Global power balance equation for electrons

    i_rad_loss: switch for radiation loss term usage in power balance
    - 0 total power lost is scaling power plus radiation (needed for i_plasma_pedestal=2,3)
    - 1 total power lost is scaling power plus core radiation only
    - 2 total power lost is scaling power only, with no additional
        allowance for radiation. This is not recommended for power plant models.

    i_plasma_ignited: switch for ignition assumption:
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)

    pden_electron_transport_loss_mw: electron transport power per volume (MW/m3)
    pden_plasma_rad_mw: total radiation power per volume (MW/m3)
    pden_plasma_core_rad_mw: total core radiation power per volume (MW/m3)
    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    f_pden_alpha_electron_mw: alpha power per volume to electrons (MW/m3)
    pden_ion_electron_equilibration_mw: ion/electron equilibration power per volume (MW/m3)
    p_hcd_injected_electrons_mw: auxiliary injected power to electrons (MW)
    vol_plasma: plasma volume (m3)
    """
    # pscaling: total transport power per volume (MW/m3)

    pscaling = data_structure.physics_variables.pden_electron_transport_loss_mw
    # Total power lost is scaling power plus radiation:
    if data_structure.physics_variables.i_rad_loss == 0:
        pnumerator = pscaling + data_structure.physics_variables.pden_plasma_rad_mw
    elif data_structure.physics_variables.i_rad_loss == 1:
        pnumerator = pscaling + data_structure.physics_variables.pden_plasma_core_rad_mw
    else:
        pnumerator = pscaling

    # if plasma not ignited include injected power
    if data_structure.physics_variables.i_plasma_ignited == 0:
        pdenom = (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_electron_mw
            + data_structure.physics_variables.pden_ion_electron_equilibration_mw
            + data_structure.current_drive_variables.p_hcd_injected_electrons_mw
            / data_structure.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_electron_mw
            + data_structure.physics_variables.pden_ion_electron_equilibration_mw
        )

    return eq(pnumerator, pdenom, constraint_registration)


@ConstraintManager.register_constraint(5, "/m3", "<=")
def constraint_equation_5(constraint_registration):
    """Equation for density upper limit

    fdene: density limit scale
    nd_plasma_electrons_vol_avg: electron density (/m3)
    nd_plasma_electrons_max: density limit (/m3)
    nd_plasma_electron_line: line averaged electron density (m-3)

    i_density_limit:
    - 1 old ASDEX;
    - 2 Borrass model for ITER (I);
    - 3 Borrass model for ITER (II);
    - 4 JET edge radiation;
    - 5 JET simplified;
    - 6 Hugill-Murakami Mq limit;
    - 7 Greenwald limit

    fdene scales the constraint such that:
    nd_plasma_electrons_vol_avg / nd_plasma_electrons_max <= fdene.
    (Except when i_density_limit=7 when nd_plasma_electron_line is used, not nd_plasma_electrons_vol_avg)
    """
    # Apply Greenwald limit to line-averaged density
    if data_structure.physics_variables.i_density_limit == 7:
        return leq(
            data_structure.physics_variables.nd_plasma_electron_line,
            (
                data_structure.physics_variables.nd_plasma_electrons_max
                * data_structure.constraint_variables.fdene
            ),
            constraint_registration,
        )

    return leq(
        data_structure.physics_variables.nd_plasma_electrons_vol_avg,
        (
            data_structure.physics_variables.nd_plasma_electrons_max
            * data_structure.constraint_variables.fdene
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(6, "", "<=")
def constraint_equation_6(constraint_registration):
    """Equation for epsilon beta-poloidal upper limit

    beta_poloidal_eps_max: maximum (eps*beta_poloidal)
    eps: inverse aspect ratio
    beta_poloidal: poloidal beta
    """
    return leq(
        (
            data_structure.physics_variables.eps
            * data_structure.physics_variables.beta_poloidal_vol_avg
        ),
        data_structure.physics_variables.beta_poloidal_eps_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(7, "/m3", "=")
def constraint_equation_7(constraint_registration):
    """Equation for hot beam ion density

    i_plasma_ignited: switch for ignition assumption:
    - 0 do not assume plasma ignition
    - 1 assume ignited (but include auxiliary power in costs)
    Obviously, i_plasma_ignited must be zero if current drive is required.
    If i_plasma_ignited=1, any auxiliary power is assumed to be used only
    during plasma start-up, and is excluded from all steady-state
    power balance calculations.
    nd_beam_ions_out: hot beam ion density from calculation (/m3)
    nd_beam_ions: hot beam ion density, variable (/m3)
    """
    if data_structure.physics_variables.i_plasma_ignited == 1:
        raise ProcessValueError("Do not use constraint equation 7 if i_plasma_ignited=1")

    return eq(
        data_structure.physics_variables.nd_beam_ions_out,
        data_structure.physics_variables.nd_beam_ions,
        constraint_registration,
    )


@ConstraintManager.register_constraint(8, "MW/m2", "<=")
def constraint_equation_8(constraint_registration):
    """Equation for neutron wall load upper limit

    pflux_fw_neutron_max_mw: allowable wall-load (MW/m2)
    pflux_fw_neutron_mw: average neutron wall load (MW/m2)
    """
    return leq(
        data_structure.physics_variables.pflux_fw_neutron_mw,
        data_structure.constraint_variables.pflux_fw_neutron_max_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(9, "MW", "<=")
def constraint_equation_9(constraint_registration):
    """Equation for fusion power upper limit

    p_fusion_total_max_mw: maximum fusion power (MW)
    p_fusion_total_mw: fusion power (MW)
    """
    return leq(
        data_structure.physics_variables.p_fusion_total_mw,
        data_structure.constraint_variables.p_fusion_total_max_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(11, "m", "=")
def constraint_equation_11(constraint_registration):
    """Equation for radial build

    rbld: sum of thicknesses to the major radius (m)
    rmajor: plasma major radius (m)
    """
    return eq(
        data_structure.build_variables.rbld,
        data_structure.physics_variables.rmajor,
        constraint_registration,
    )


@ConstraintManager.register_constraint(12, "V.sec", ">=")
def constraint_equation_12(constraint_registration):
    """Equation for volt-second capability lower limit

    vs_plasma_total_required: total V-s needed (Wb)
    vs_plasma_total_required (lower limit) is positive; vs_cs_pf_total_pulse (available) is negative
    vs_cs_pf_total_pulse: total flux swing for pulse (Wb)
    """
    # vs_cs_pf_total_pulse is negative, requires sign change
    return geq(
        -data_structure.pfcoil_variables.vs_cs_pf_total_pulse,
        data_structure.physics_variables.vs_plasma_total_required,
        constraint_registration,
    )


@ConstraintManager.register_constraint(13, "sec", ">=")
def constraint_equation_13(constraint_registration):
    """Equation for burn time lower limit

    t_plant_pulse_burn: burn time (s) (calculated if i_pulsed_plant=1)
    t_burn_min: minimum burn time (s)
    """
    return geq(
        data_structure.times_variables.t_plant_pulse_burn,
        data_structure.constraint_variables.t_burn_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(14, "", "=")
def constraint_equation_14(constraint_registration):
    """Equation to fix number of NBI decay lengths to plasma centre

    n_beam_decay_lengths_core: neutral beam e-decay lengths to plasma centre
    n_beam_decay_lengths_core_required: permitted neutral beam e-decay lengths to plasma centre
    """
    return eq(
        data_structure.current_drive_variables.n_beam_decay_lengths_core,
        data_structure.current_drive_variables.n_beam_decay_lengths_core_required,
        constraint_registration,
    )


@ConstraintManager.register_constraint(15, "MW", ">=")
def constraint_equation_15(constraint_registration):
    """Equation for L-H power threshold limit to enforce H-mode

    f_h_mode_margin: a margin on the constraint
    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)

    Setting f_h_mode_margin != 1.0 enforces a margin on the constraint:
    I.e.  p_plasma_separatrix_mw >= f_h_mode_margin * p_l_h_threshold_mw

    For example, f_h_mode_margin = 1.2 will ensure that
    p_plasma_separatrix_mw is at least 1.2*p_l_h_threshold_mw (ie in H-mode).
    """
    return geq(
        data_structure.physics_variables.p_plasma_separatrix_mw,
        (
            data_structure.physics_variables.p_l_h_threshold_mw
            * data_structure.constraint_variables.f_h_mode_margin
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(16, "MW", ">=")
def constraint_equation_16(constraint_registration):
    """Equation for net electric power lower limit

    p_plant_electric_net_mw: net electric power (MW)
    p_plant_electric_net_required_mw: required net electric power (MW)
    """
    return geq(
        data_structure.heat_transport_variables.p_plant_electric_net_mw,
        data_structure.constraint_variables.p_plant_electric_net_required_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(17, "MW/m3", "<=")
def constraint_equation_17(constraint_registration):
    """Equation for radiation power upper limit

    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    vol_plasma: plasma volume (m3)
    pden_alpha_total_mw: alpha power per volume (MW/m3)
    pden_non_alpha_charged_mw: non-alpha charged particle fusion power per volume (MW/m3)
    pden_plasma_ohmic_mw: ohmic heating power per volume (MW/m3)
    pden_plasma_rad_mw: total radiation power per volume (MW/m3)
    fradpwr: core radiation power limit scale

    fradpwr adds a margin to the constraint constraint such that

    pden_plasma_rad_mw / pradmaxpv <= fradpwr
    """
    # Maximum possible power/vol_plasma that can be radiated (local)
    pradmaxpv = (
        data_structure.current_drive_variables.p_hcd_injected_total_mw
        / data_structure.physics_variables.vol_plasma
        + data_structure.physics_variables.pden_alpha_total_mw
        * data_structure.physics_variables.f_p_alpha_plasma_deposited
        + data_structure.physics_variables.pden_non_alpha_charged_mw
        + data_structure.physics_variables.pden_plasma_ohmic_mw
    )

    return leq(
        data_structure.physics_variables.pden_plasma_rad_mw / pradmaxpv,
        data_structure.constraint_variables.fradpwr,
        constraint_registration,
    )


@ConstraintManager.register_constraint(18, "MW/m2", "<=")
def constraint_equation_18(constraint_registration):
    """Equation for divertor heat load upper limit

    pflux_div_heat_load_max_mw: heat load limit (MW/m2)
    pflux_div_heat_load_mw: divertor heat load (MW/m2)
    """
    return leq(
        data_structure.divertor_variables.pflux_div_heat_load_mw,
        data_structure.divertor_variables.pflux_div_heat_load_max_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(19, "MVA", "<=")
def constraint_equation_19(constraint_registration):
    """Equation for MVA (power) upper limit: resistive TF coil set

    p_cp_resistive_mw: peak resistive TF coil inboard leg power (total) (MW)
    p_tf_leg_resistive_mw: TF coil outboard leg resistive power (total) (MW)
    mvalim: MVA limit for resistive TF coil set (total) (MW)
    """
    totmva = (
        data_structure.tfcoil_variables.p_cp_resistive_mw
        + data_structure.tfcoil_variables.p_tf_leg_resistive_mw
    )

    return leq(
        totmva, data_structure.constraint_variables.mvalim, constraint_registration
    )


@ConstraintManager.register_constraint(20, "m", "<=")
def constraint_equation_20(constraint_registration):
    """Equation for neutral beam tangency radius upper limit

    radius_beam_tangency_max: maximum tangency radius for centreline of beam (m)
    radius_beam_tangency: neutral beam centreline tangency radius (m)
    """
    return leq(
        data_structure.current_drive_variables.radius_beam_tangency,
        data_structure.current_drive_variables.radius_beam_tangency_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(21, "", ">=")
def constraint_equation_21(constraint_registration):
    """Equation for minor radius lower limit

    rminor: plasma minor radius (m)
    aplasmin: minimum minor radius (m)
    """
    return geq(
        data_structure.physics_variables.rminor,
        data_structure.build_variables.aplasmin,
        constraint_registration,
    )


@ConstraintManager.register_constraint(22, "MW", ">=")
def constraint_equation_22(constraint_registration):
    """Equation for L-H power threshold limit to enforce L-mode

    f_l_mode_margin: a margin on the constraint
    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)

    Setting f_l_mode_margin != 1.0 enforces a margin on the constraint:
    I.e.  p_l_h_threshold_mw >= f_l_mode_margin * p_plasma_separatrix_mw

    For example, f_l_mode_margin = 1.2 will ensure that
    p_l_h_threshold_mw is at least 1.2*p_plasma_separatrix_mw (ie in L-mode).
    """
    return geq(
        data_structure.physics_variables.p_l_h_threshold_mw,
        (
            data_structure.constraint_variables.f_l_mode_margin
            * data_structure.physics_variables.p_plasma_separatrix_mw
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(23, "m", "<=")
def constraint_equation_23(constraint_registration):
    """Equation for conducting shell radius / rminor upper limit

    rminor: plasma minor radius (m)
    dr_fw_plasma_gap_outboard: gap between plasma and first wall, outboard side (m)
    dr_fw_outboard: outboard first wall thickness, initial estimate (m)
    dr_blkt_outboard: outboard blanket thickness (m)
    f_r_conducting_wall: maximum ratio of conducting wall distance to plasma minor radius for vertical stability
    """
    # conducting shell radius (m)
    rcw = (
        data_structure.physics_variables.rminor
        + data_structure.build_variables.dr_fw_plasma_gap_outboard
        + data_structure.build_variables.dr_fw_outboard
        + data_structure.build_variables.dr_blkt_outboard
    )
    return leq(
        rcw,
        (
            data_structure.physics_variables.f_r_conducting_wall
            * data_structure.physics_variables.rminor
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(24, "", "<=")
def constraint_equation_24(constraint_registration):
    """Equation for beta upper limit

    i_beta_component: switch for beta limit scaling (constraint equation  24):
    - 0 apply limit to total beta;
    - 1 apply limit to thermal beta;
    - 2 apply limit to thermal + neutral beam beta
    - 3 apply limit to toroidal beta
    istell: switch for stellarator option:
    - 0 use tokamak model;
    - 1 use stellarator model
    beta_vol_avg_max: allowable beta
    beta_total_vol_avg: total plasma beta (calculated if i_plasma_pedestal =3)
    beta_fast_alpha: fast alpha beta component
    beta_beam: neutral beam beta component
    b_plasma_toroidal_on_axis: toroidal field
    b_plasma_total: total field
    """
    # Include all beta components: relevant for both tokamaks and stellarators
    if (
        data_structure.physics_variables.i_beta_component == 0
        or data_structure.stellarator_variables.istell != 0
    ):
        value = data_structure.physics_variables.beta_total_vol_avg
    # Here, the beta limit applies to only the thermal component, not the fast alpha or neutral beam parts
    elif data_structure.physics_variables.i_beta_component == 1:
        value = (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
            - data_structure.physics_variables.beta_beam
        )
    # Beta limit applies to thermal + neutral beam: components of the total beta, i.e. excludes alphas
    elif data_structure.physics_variables.i_beta_component == 2:
        value = (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
        )
    # Beta limit applies to toroidal beta
    elif data_structure.physics_variables.i_beta_component == 3:
        value = (
            data_structure.physics_variables.beta_total_vol_avg
            * (
                data_structure.physics_variables.b_plasma_total
                / data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )

    return leq(
        value,
        data_structure.physics_variables.beta_vol_avg_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(25, "T", "<=")
def constraint_equation_25(constraint_registration):
    """Equation for peak toroidal field upper limit

    b_tf_inboard_max: maximum peak toroidal field (T)
    b_tf_inboard_peak_symmetric: mean peak field at TF coil (T)
    """
    return leq(
        data_structure.tfcoil_variables.b_tf_inboard_peak_symmetric,
        data_structure.constraint_variables.b_tf_inboard_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(26, "A/m2", "<=")
def constraint_equation_26(constraint_registration):
    """Equation for Central Solenoid current density upper limit at EOF

    fjohc: margin for central solenoid current at end-of-flattop
    j_cs_critical_flat_top_end: allowable central solenoid current density at end of flat-top (A/m2)
    j_cs_flat_top_end: central solenoid overall current density at end of flat-top (A/m2)
    """
    return leq(
        (
            data_structure.pfcoil_variables.j_cs_flat_top_end
            / data_structure.pfcoil_variables.j_cs_critical_flat_top_end
        ),
        data_structure.constraint_variables.fjohc,
        constraint_registration,
    )


@ConstraintManager.register_constraint(27, "A/m2", "<=")
def constraint_equation_27(constraint_registration):
    """Equation for Central Solenoid current density upper limit at BOP

    fjohc0: margin for central solenoid current at beginning of pulse
    j_cs_critical_pulse_start: allowable central solenoid current density at beginning of pulse (A/m2)
    j_cs_pulse_start: central solenoid overall current density at beginning of pulse (A/m2)
    """
    return leq(
        (
            data_structure.pfcoil_variables.j_cs_pulse_start
            / data_structure.pfcoil_variables.j_cs_critical_pulse_start
        ),
        data_structure.constraint_variables.fjohc0,
        constraint_registration,
    )


@ConstraintManager.register_constraint(28, "", ">=")
def constraint_equation_28(constraint_registration):
    """Equation for fusion gain (big Q) lower limit

    big_q_plasma: Fusion gain; P_fusion / (P_injection + P_ohmic)
    big_q_plasma_min: minimum fusion gain Q
    i_plasma_ignited : input integer : switch for ignition assumption:
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)
    Obviously, ignite must be zero if current drive is required.
    If i_plasma_ignited=1, any auxiliary power is assumed to be used only
    during plasma start-up, and is excluded from all steady-state
    power balance calculations.
    """
    if data_structure.physics_variables.i_plasma_ignited != 0:
        raise ProcessValueError("Do not use constraint 28 if i_plasma_ignited=1")

    return geq(
        data_structure.current_drive_variables.big_q_plasma,
        data_structure.constraint_variables.big_q_plasma_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(29, "m", "=")
def constraint_equation_29(constraint_registration):
    """Equation for inboard major radius: This is a consistency equation

    rmajor: plasma major radius (m) (iteration variable 3)
    rminor: plasma minor radius (m)
    rinboard: plasma inboard radius (m)
    """
    return eq(
        (
            data_structure.physics_variables.rmajor
            - data_structure.physics_variables.rminor
        ),
        data_structure.build_variables.rinboard,
        constraint_registration,
    )


@ConstraintManager.register_constraint(30, "MW", "<=")
def constraint_equation_30(constraint_registration):
    """Equation for injection power upper limit

    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    p_hcd_injected_max: Maximum allowable value for injected power (MW)
    """
    return leq(
        data_structure.current_drive_variables.p_hcd_injected_total_mw,
        data_structure.current_drive_variables.p_hcd_injected_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(31, "Pa", "<=")
def constraint_equation_31(constraint_registration):
    """Equation for TF coil case stress upper limit (SCTF)

    sig_tf_case_max: Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)
    sig_tf_case: Constrained stress in TF coil case (Pa)
    """
    return leq(
        data_structure.tfcoil_variables.sig_tf_case,
        data_structure.tfcoil_variables.sig_tf_case_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(32, "Pa", "<=")
def constraint_equation_32(constraint_registration):
    """Equation for TF coil conduit stress upper limit (SCTF)

    sig_tf_wp_max: Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)
    sig_tf_wp: Constrained stress in TF conductor conduit (Pa)
    """
    return leq(
        data_structure.tfcoil_variables.sig_tf_wp,
        data_structure.tfcoil_variables.sig_tf_wp_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(33, "A/m2", "<=")
def constraint_equation_33(constraint_registration):
    """Equation for TF coil operating/critical J upper limit (SCTF)

    args : output structure : residual error; constraint value;

    fiooic: margin for TF coil operating current / critical
    j_tf_wp_critical: critical current density for winding pack (A/m2)
    j_tf_wp: winding pack current density (A/m2)

    fiooic scales the constraint such that:
    j_tf_wp / j_tf_wp_critical <= fiooic.
    """
    return leq(
        data_structure.tfcoil_variables.j_tf_wp,
        (
            data_structure.tfcoil_variables.j_tf_wp_critical
            * data_structure.constraint_variables.fiooic
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(34, "V", "<=")
def constraint_equation_34(constraint_registration):
    """Equation for TF coil dump voltage upper limit (SCTF)

    v_tf_coil_dump_quench_max_kv: max voltage across TF coil during quench (kV)
    v_tf_coil_dump_quench_kv: voltage across a TF coil during quench (kV)
    """
    return leq(
        data_structure.tfcoil_variables.v_tf_coil_dump_quench_kv,
        data_structure.tfcoil_variables.v_tf_coil_dump_quench_max_kv,
        constraint_registration,
    )


@ConstraintManager.register_constraint(35, "A/m2", "<=")
def constraint_equation_35(constraint_registration):
    """Equation for TF coil J_wp/J_prot upper limit (SCTF)

    j_tf_wp_quench_heat_max: allowable TF coil winding pack current density, for dump temperature
    rise protection (A/m2)
    j_tf_wp: winding pack current density (A/m2)
    """
    return leq(
        data_structure.tfcoil_variables.j_tf_wp,
        data_structure.tfcoil_variables.j_tf_wp_quench_heat_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(36, "K", ">=")
def constraint_equation_36(constraint_registration):
    """Equation for TF coil s/c temperature margin lower limit (SCTF)

    temp_tf_superconductor_margin: TF coil temperature margin (K)
    temp_tf_superconductor_margin_min: minimum allowable temperature margin : TF coils (K)
    """
    return geq(
        data_structure.tfcoil_variables.temp_tf_superconductor_margin,
        data_structure.tfcoil_variables.temp_tf_superconductor_margin_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(37, "1E20 A/Wm2", "<=")
def constraint_equation_37(constraint_registration):
    """Equation for current drive gamma upper limit

    eta_cd_norm_hcd_primary_max: maximum current drive gamma
    eta_cd_norm_hcd_primary: normalised current drive efficiency (1.0e20 A/W-m2)
    """
    return leq(
        data_structure.current_drive_variables.eta_cd_norm_hcd_primary,
        data_structure.constraint_variables.eta_cd_norm_hcd_primary_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(39, "K", "<=")
def constraint_equation_39(constraint_registration):
    """Equation for first wall temperature upper limit

    temp_fw_max: maximum temperature of first wall material (K) (i_thermal_electric_conversion>1)
    temp_fw_peak: peak first wall temperature (K)
    """
    if data_structure.fwbs_variables.temp_fw_peak < 1.0:
        raise ProcessValueError(
            "temp_fw_peak = 0 implies i_pulsed_plant=0; do not use constraint 39 if i_pulsed_plant=0"
        )

    return leq(
        data_structure.fwbs_variables.temp_fw_peak,
        data_structure.fwbs_variables.temp_fw_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(40, "MW", ">=")
def constraint_equation_40(constraint_registration):
    """Equation for auxiliary power lower limit

    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    p_hcd_injected_min_mw: minimum auxiliary power (MW)
    """
    return geq(
        data_structure.current_drive_variables.p_hcd_injected_total_mw,
        data_structure.constraint_variables.p_hcd_injected_min_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(41, "sec", ">=")
def constraint_equation_41(constraint_registration):
    """Equation for plasma current ramp-up time lower limit

    t_plant_pulse_plasma_current_ramp_up: plasma current ramp-up time for current initiation (s)
    t_current_ramp_up_min: minimum plasma current ramp-up time (s)
    """
    return geq(
        data_structure.times_variables.t_plant_pulse_plasma_current_ramp_up,
        data_structure.constraint_variables.t_current_ramp_up_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(42, "sec", ">=")
def constraint_equation_42(constraint_registration):
    """Equation for cycle time lower limit

    t_plant_pulse_total: full cycle time (s)
    t_cycle_min: minimum cycle time (s)
    """
    if data_structure.constraint_variables.t_cycle_min < 1.0:
        raise ProcessValueError(
            "t_cycle_min = 0 implies that i_pulsed_plant=0; do not use constraint 42 if i_pulsed_plant=0"
        )

    return geq(
        data_structure.times_variables.t_plant_pulse_total,
        data_structure.constraint_variables.t_cycle_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(43, "deg C", "=")
def constraint_equation_43(constraint_registration):
    """Equation for average centrepost temperature: This is a consistency equation (TART)

    temp_cp_average: average temp of TF coil inboard leg conductor (C)e
    tcpav2: centrepost average temperature (C) (for consistency)
    itart: switch for spherical tokamak (ST) models:
    - 0 use conventional aspect ratio models;
    - 1 use spherical tokamak models
    """
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 43 if itart=0")

    if data_structure.tfcoil_variables.i_tf_sup == 0:
        temp_cp_average = data_structure.tfcoil_variables.temp_cp_average - 273.15
        tcpav2 = data_structure.tfcoil_variables.tcpav2 - 273.15
    else:
        temp_cp_average = data_structure.tfcoil_variables.temp_cp_average
        tcpav2 = data_structure.tfcoil_variables.tcpav2

    return eq(temp_cp_average, tcpav2, constraint_registration)


@ConstraintManager.register_constraint(44, "deg C", "<=")
def constraint_equation_44(constraint_registration):
    """Equation for centrepost temperature upper limit (TART)

    temp_cp_max: maximum peak centrepost temperature (K)
    temp_cp_peak: peak centrepost temperature (K)
    itart: switch for spherical tokamak (ST) models:
    - 0: use conventional aspect ratio models;
    - 1: use spherical tokamak models
    """
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 44 if itart=0")

    if data_structure.tfcoil_variables.i_tf_sup == 0:  # ! Copper case
        temp_cp_max = data_structure.tfcoil_variables.temp_cp_max - 273.15
        temp_cp_peak = data_structure.tfcoil_variables.temp_cp_peak - 273.15
    else:
        temp_cp_max = data_structure.tfcoil_variables.temp_cp_max
        temp_cp_peak = data_structure.tfcoil_variables.temp_cp_peak

    return leq(temp_cp_peak, temp_cp_max, constraint_registration)


@ConstraintManager.register_constraint(45, "", ">=")
def constraint_manager_45(constraint_registration):
    """Equation for edge safety factor lower limit (TART)

    q95 : safety factor 'near' plasma edge
    (unless i_plasma_current = 2 (ST current scaling), in which case q = mean edge safety factor qbar)
    q95_min: lower limit for edge safety factor
    itart : input integer : switch for spherical tokamak (ST) models:
    - 0 use conventional aspect ratio models;
    - 1 use spherical tokamak models
    """
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 45 if itart=0")

    return geq(
        data_structure.physics_variables.q95,
        data_structure.physics_variables.q95_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(46, "", "<=")
def constraint_equation_46(constraint_registration):
    """Equation for Ip/Irod upper limit (TART)

    eps: inverse aspect ratio
    c_tf_total: total (summed) current in TF coils (A)
    plasma_current: plasma current (A)
    itart: switch for spherical tokamak (ST) models:
    - 0: use conventional aspect ratio models;
    - 1: use spherical tokamak models
    """
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 46 if itart=0")

    # maximum ratio of plasma current to centrepost current
    cratmx = 1.0 + 4.91 * (data_structure.physics_variables.eps - 0.62)

    return leq(
        (
            data_structure.physics_variables.plasma_current
            / data_structure.tfcoil_variables.c_tf_total
        ),
        cratmx,
        constraint_registration,
    )


@ConstraintManager.register_constraint(48, "", "<=")
def constraint_equation_48(constraint_registration):
    """Equation for poloidal beta upper limit

    beta_poloidal_max: maximum poloidal beta
    beta_poloidal: poloidal beta
    """
    return leq(
        data_structure.physics_variables.beta_poloidal_vol_avg,
        data_structure.constraint_variables.beta_poloidal_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(50, "Hz", "<=")
def constraint_equation_50(constraint_registration):
    """IFE option: Equation for repetition rate upper limit"""
    return leq(
        data_structure.ife_variables.reprat,
        data_structure.ife_variables.rrmax,
        constraint_registration,
    )


@ConstraintManager.register_constraint(51, "V.s", "=")
def constraint_equation_51(constraint_registration):
    """Equation to enforce startup flux = available startup flux

    vs_plasma_res_ramp: resistive losses in startup V-s (Wb)
    vs_plasma_ind_ramp: internal and external plasma inductance V-s (Wb))
    vs_cs_pf_total_ramp: total flux swing for startup (Wb)
    """
    return eq(
        abs(
            data_structure.physics_variables.vs_plasma_res_ramp
            + data_structure.physics_variables.vs_plasma_ind_ramp
        ),
        data_structure.pfcoil_variables.vs_cs_pf_total_ramp,
        constraint_registration,
    )


@ConstraintManager.register_constraint(52, "", ">=")
def constraint_equation_52(constraint_registration):
    """Equation for tritium breeding ratio lower limit

    The tritium breeding ratio is only calculated when using the IFE model.

    tbr: tritium breeding ratio
    tbrmin: minimum tritium breeding ratio
    """
    if data_structure.ife_variables.ife != 1:
        raise ProcessValueError(
            "Constraint 52 is only supported when running the IFE model"
        )

    return geq(
        data_structure.fwbs_variables.tbr,
        data_structure.constraint_variables.tbrmin,
        constraint_registration,
    )


@ConstraintManager.register_constraint(53, "neutron/m2", "<=")
def constraint_equation_53(constraint_registration):
    """Equation for fast neutron fluence on TF coil upper limit

    nflutfmax: max fast neutron fluence on TF coil (n/m2)
    nflutf: peak fast neutron fluence on TF coil superconductor (n/m2)
    """
    return leq(
        data_structure.fwbs_variables.nflutf,
        data_structure.constraint_variables.nflutfmax,
        constraint_registration,
    )


@ConstraintManager.register_constraint(54, "MW/m3", "<=")
def constraint_equation_54(constraint_registration):
    """Equation for peak TF coil nuclear heating upper limit

    ptfnucmax: maximum nuclear heating in TF coil (MW/m3)
    ptfnucpm3: nuclear heating in the TF coil (MW/m3) (blktmodel>0)
    """
    return leq(
        data_structure.fwbs_variables.ptfnucpm3,
        data_structure.constraint_variables.ptfnucmax,
        constraint_registration,
    )


@ConstraintManager.register_constraint(56, "MW/m", "<=")
def constraint_equation_56(constraint_registration):
    """Equation for power through separatrix / major radius upper limit

    pseprmax: maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    rmajor: plasma major radius (m)
    """
    return leq(
        (
            data_structure.physics_variables.p_plasma_separatrix_mw
            / data_structure.physics_variables.rmajor
        ),
        data_structure.constraint_variables.pseprmax,
        constraint_registration,
    )


@ConstraintManager.register_constraint(59, "", "<=")
def constraint_equation_59(constraint_registration):
    """Equation for neutral beam shine-through fraction upper limit

    f_p_beam_shine_through_max: maximum neutral beam shine-through fraction
    f_p_beam_shine_through: neutral beam shine-through fraction
    """
    return leq(
        data_structure.current_drive_variables.f_p_beam_shine_through,
        data_structure.constraint_variables.f_p_beam_shine_through_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(60, "K", ">=")
def constraint_equation_60(constraint_registration):
    """Equation for Central Solenoid s/c temperature margin lower limit

    temp_cs_superconductor_margin: Central solenoid temperature margin (K)
    temp_cs_superconductor_margin_min: Minimum allowable temperature margin : CS (K)
    """
    return geq(
        data_structure.pfcoil_variables.temp_cs_superconductor_margin,
        data_structure.tfcoil_variables.temp_cs_superconductor_margin_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(61, "", ">=")
def constraint_equation_61(constraint_registration):
    """Equation for availability lower limit

    f_t_plant_available: Total plant availability fraction
    avail_min: Minimum availability
    """
    return geq(
        data_structure.cost_variables.f_t_plant_available,
        data_structure.cost_variables.avail_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(62, "", ">=")
def constraint_equation_62(constraint_registration):
    """Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times

    t_alpha_confinement: alpha particle confinement time (s)
    t_energy_confinement: global thermal energy confinement time (sec)
    f_alpha_energy_confinement_min: Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
    """
    return geq(
        data_structure.physics_variables.t_alpha_confinement
        / data_structure.physics_variables.t_energy_confinement,
        data_structure.constraint_variables.f_alpha_energy_confinement_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(63, "", "<=")
def constraint_equation_63(constraint_registration):
    """Upper limit on n_iter_vacuum_pumps (i_vacuum_pumping = simple)

    tfno: number of TF coils (default = 50 for stellarators)
    n_iter_vacuum_pumps: number of high vacuum pumps (real number), each with the throughput
    """
    return leq(
        data_structure.vacuum_variables.n_iter_vacuum_pumps,
        data_structure.tfcoil_variables.n_tf_coils,
        constraint_registration,
    )


@ConstraintManager.register_constraint(64, "", "<=")
def constraint_equation_64(constraint_registration):
    """Upper limit on Zeff

    zeff_max: maximum value for Zeff
    n_charge_plasma_effective_vol_avg: plasma effective charge
    """
    return leq(
        data_structure.physics_variables.n_charge_plasma_effective_vol_avg,
        data_structure.constraint_variables.zeff_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(65, "Pa", "<=")
def constraint_equation_65(constraint_registration):
    """Upper limit on stress of the vacuum vessel that occurs when the TF coil quenches.

    max_vv_stress: Maximum permitted stress of the VV (Pa)
    vv_stress_quench: Stress of the VV (Pa)
    """
    return leq(
        data_structure.superconducting_tf_coil_variables.vv_stress_quench,
        data_structure.tfcoil_variables.max_vv_stress,
        constraint_registration,
    )


@ConstraintManager.register_constraint(66, "MW", "<=")
def constrain_equation_66(constraint_registration):
    """Upper limit on rate of change of energy in poloidal field

    maxpoloidalpower: Maximum permitted absolute rate of change of stored energy in poloidal field (MW)
    peakpoloidalpower: Peak absolute rate of change of stored energy in poloidal field (MW) (11/01/16)
    """
    return leq(
        data_structure.pf_power_variables.peakpoloidalpower,
        data_structure.pf_power_variables.maxpoloidalpower,
        constraint_registration,
    )


@ConstraintManager.register_constraint(67, "MW/m2", "<=")
def constraint_equation_67(constraint_registration):
    """Simple upper limit on radiation wall load

    pflux_fw_rad_max: Maximum permitted radiation wall load (MW/m^2)
    pflux_fw_rad_max_mw: Peak radiation wall load (MW/m^2)
    """
    return leq(
        data_structure.constraint_variables.pflux_fw_rad_max_mw,
        data_structure.constraint_variables.pflux_fw_rad_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(68, "MWT/m", "<=")
def constraint_equation_68(constraint_registration):
    """Upper limit on Psep scaling (PsepB/qAR)

    psepbqarmax: maximum permitted value of ratio of Psep*Bt/qAR (MWT/m)
    p_plasma_separatrix_mw: Power to conducted to the divertor region (MW)
    b_plasma_toroidal_on_axis: toroidal field on axis (T) (iteration variable 2)
    q95: safety factor q at 95% flux surface
    aspect: aspect ratio (iteration variable 1)
    rmajor: plasma major radius (m) (iteration variable 3)
    i_q95_fixed: Switch that allows for fixing q95 only in this constraint.
    q95_fixed: fixed safety factor q at 95% flux surface
    """
    if data_structure.constraint_variables.i_q95_fixed == 1:
        return leq(
            (
                (
                    data_structure.physics_variables.p_plasma_separatrix_mw
                    * data_structure.physics_variables.b_plasma_toroidal_on_axis
                )
                / (
                    data_structure.constraint_variables.q95_fixed
                    * data_structure.physics_variables.aspect
                    * data_structure.physics_variables.rmajor
                )
            ),
            data_structure.constraint_variables.psepbqarmax,
            constraint_registration,
        )

    return leq(
        (
            (
                data_structure.physics_variables.p_plasma_separatrix_mw
                * data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            / (
                data_structure.physics_variables.q95
                * data_structure.physics_variables.aspect
                * data_structure.physics_variables.rmajor
            )
        ),
        data_structure.constraint_variables.psepbqarmax,
        constraint_registration,
    )


@ConstraintManager.register_constraint(72, "Pa", "<=")
def constraint_equation_72(constraint_registration):
    """Upper limit on central Solenoid Tresca yield stress

    In the case if the bucked and wedged option ( i_tf_bucking >= 2 ) the constrained
    stress is the largest the largest stress of the
     - CS stress at maximum current (conservative as the TF inward pressure is not taken
       into account)
     - CS stress at flux swing (no current in CS) from the TF inward pressure
    This allow to cover the 2 worst stress scenario in the bucked and wedged design
    Otherwise (free standing TF), the stress limits are only set by the CS stress at max current
    Reverse the sign so it works as an inequality constraint (tmp_cc > 0)
    This will have no effect if it is used as an equality constraint because it will be squared.

    alstroh: allowable hoop stress in Central Solenoid structural material (Pa)
    s_shear_cs_peak: Maximum shear stress coils/central solenoid (Pa)
    sig_tf_cs_bucked: Maximum shear stress in CS case at flux swing (no current in CS)
                          can be significant for the bucked and weged design
    i_tf_bucking: switch for TF structure design
    """
    # bucked and wedged desing
    if (
        data_structure.tfcoil_variables.i_tf_bucking >= 2
        and data_structure.build_variables.i_tf_inside_cs == 0
    ):
        return leq(
            max(
                data_structure.pfcoil_variables.s_shear_cs_peak,
                data_structure.tfcoil_variables.sig_tf_cs_bucked,
            ),
            data_structure.pfcoil_variables.alstroh,
            constraint_registration,
        )
    # Free standing CS
    return leq(
        data_structure.pfcoil_variables.s_shear_cs_peak,
        data_structure.pfcoil_variables.alstroh,
        constraint_registration,
    )


@ConstraintManager.register_constraint(73, "MW", ">=")
def constraint_equation_73(constraint_registration):
    """Lower limit to ensure separatrix power is greater than the L-H power + auxiliary power
    Related to constraint 15

    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    p_hcd_injected_total_mw : inout real : total auxiliary injected power (MW)
    """
    return geq(
        data_structure.physics_variables.p_plasma_separatrix_mw,
        (
            data_structure.physics_variables.p_l_h_threshold_mw
            + data_structure.current_drive_variables.p_hcd_injected_total_mw
        ),
        constraint_registration,
    )


@ConstraintManager.register_constraint(74, "K", "<=")
def constraint_equation_74(constraint_registration):
    """Upper limit to ensure TF coil quench temperature < temp_croco_quench_max
    ONLY used for croco HTS coil

    temp_croco_quench: CroCo strand: Actual temp reached during a quench (K)
    temp_croco_quench_max: CroCo strand: maximum permitted temp during a quench (K)
    """
    return leq(
        data_structure.tfcoil_variables.temp_croco_quench,
        data_structure.tfcoil_variables.temp_croco_quench_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(75, "A/m2", "<=")
def constraint_equation_75(constraint_registration):
    """Upper limit to ensure that TF coil current / copper area < Maximum value
    ONLY used for croco HTS coil

    copperA_m2: TF coil current / copper area (A/m2)
    copperA_m2_max: Maximum TF coil current / copper area (A/m2)
    """
    return leq(
        data_structure.rebco_variables.coppera_m2,
        data_structure.rebco_variables.coppera_m2_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(76, "m-3", "<=")
def constraint_equation_76(constraint_registration):
    """Upper limit for Eich critical separatrix density model: Added for issue 558

    Eich critical separatrix density model
    Added for issue 558 with ref to http://iopscience.iop.org/article/10.1088/1741-4326/aaa340/pdf

    alpha_crit: critical ballooning parameter value
    nd_plasma_separatrix_electron_eich_max: critical electron density at separatrix [m-3]
    kappa: plasma separatrix elongation (calculated if i_plasma_geometry = 1-5, 7 or 9)
    triang: plasma separatrix triangularity (calculated if i_plasma_geometry = 1, 3-5 or 7)
    aspect: aspect ratio (iteration variable 1)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)
    nd_plasma_electron_max_array(7)array : density limit (/m3) as calculated using various models
    """
    # TODO: why on earth are these variables being set here!? Should they be local?
    data_structure.physics_variables.alpha_crit = (
        data_structure.physics_variables.kappa**1.2
    ) * (1.0 + 1.5 * data_structure.physics_variables.triang)
    data_structure.physics_variables.nd_plasma_separatrix_electron_eich_max = (
        5.9
        * data_structure.physics_variables.alpha_crit
        * (data_structure.physics_variables.aspect ** (-2.0 / 7.0))
        * (((1.0 + (data_structure.physics_variables.kappa**2.0)) / 2.0) ** (-6.0 / 7.0))
        * (
            (data_structure.physics_variables.p_plasma_separatrix_mw * 1.0e6)
            ** (-11.0 / 70.0)
        )
        * data_structure.physics_variables.nd_plasma_electron_max_array[6]
    )

    return leq(
        data_structure.physics_variables.nd_plasma_separatrix_electron,
        data_structure.physics_variables.nd_plasma_separatrix_electron_eich_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(77, "A/turn", "<=")
def constraint_equation_77(constraint_registration):
    """Equation for maximum TF current per turn upper limit

    c_tf_turn_max : allowable TF coil current per turn [A/turn]
    c_tf_turn : TF coil current per turn [A/turn]
    """
    return leq(
        data_structure.tfcoil_variables.c_tf_turn,
        data_structure.tfcoil_variables.c_tf_turn_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(78, "", ">=")
def constraint_equation_78(constraint_registration):
    """Equation for Reinke criterion, divertor impurity fraction lower limit

    fzmin : input : minimum impurity fraction from Reinke model
    fzactual : input : actual impurity fraction
    """
    return geq(
        data_structure.reinke_variables.fzactual,
        data_structure.reinke_variables.fzmin,
        constraint_registration,
    )


@ConstraintManager.register_constraint(79, "A/turn", "<=")
def constraint_equation_79(constraint_registration):
    """Equation for maximum CS field

    b_cs_limit_max: Central solenoid max field limit [T]
    b_cs_peak_pulse_start: maximum field in central solenoid at beginning of pulse (T)
    b_cs_peak_flat_top_end: maximum field in central solenoid at end of flat-top (EoF) (T)
    (Note: original code has "b_cs_peak_flat_top_end/b_cs_peak_pulse_start |  peak CS field [T]".)
    """
    return leq(
        max(
            data_structure.pfcoil_variables.b_cs_peak_flat_top_end,
            data_structure.pfcoil_variables.b_cs_peak_pulse_start,
        ),
        data_structure.pfcoil_variables.b_cs_limit_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(80, "MW", ">=")
def constraint_equation_80(constraint_registration):
    """Equation for p_plasma_separatrix_mw lower limit

    args : output structure : residual error; constraint value; residual error in physical units;
    output string; units string
    Lower limit p_plasma_separatrix_mw

    p_plasma_separatrix_min_mw : input : Minimum power crossing separatrix p_plasma_separatrix_mw [MW]
    p_plasma_separatrix_mw : input : Power crossing separatrix [MW]
    """
    return geq(
        data_structure.physics_variables.p_plasma_separatrix_mw,
        data_structure.constraint_variables.p_plasma_separatrix_min_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(81, "m-3", ">=")
def constraint_equation_81(constraint_registration):
    """Lower limit to ensure central density is larger that the pedestal one

    args : output structure : residual error; constraint value;
    residual error in physical units; output string; units string
    Lower limit nd_plasma_electron_on_axis > nd_plasma_pedestal_electron

    nd_plasma_electron_on_axis   : input : Central electron density [m-3]
    nd_plasma_pedestal_electron : input : Electron density at pedestal [m-3]
    """
    return geq(
        data_structure.physics_variables.nd_plasma_electron_on_axis,
        data_structure.physics_variables.nd_plasma_pedestal_electron,
        constraint_registration,
    )


@ConstraintManager.register_constraint(82, "m", ">=")
def constraint_equation_82(constraint_registration):
    """Equation for toroidal consistency of stellarator build

    toroidalgap: minimal gap between two stellarator coils
    dx_tf_inboard_out_toroidal: total toroidal width of a tf coil
    """
    return geq(
        data_structure.tfcoil_variables.toroidalgap,
        data_structure.tfcoil_variables.dx_tf_inboard_out_toroidal,
        constraint_registration,
    )


@ConstraintManager.register_constraint(83, "m", ">=")
def constraint_equation_83(constraint_registration):
    """Equation for radial consistency of stellarator build

    available_radial_space: avaible space in radial direction as given by each s.-configuration
    required_radial_space: required space in radial direction
    """
    return geq(
        data_structure.build_variables.available_radial_space,
        data_structure.build_variables.required_radial_space,
        constraint_registration,
    )


@ConstraintManager.register_constraint(84, "", ">=")
def constraint_equation_84(constraint_registration):
    """Equation for the lower limit of beta

    beta_vol_avg_min: Lower limit for beta
    beta: plasma beta
    """
    return geq(
        data_structure.physics_variables.beta_total_vol_avg,
        data_structure.physics_variables.beta_vol_avg_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(85, "years", "=")
def constraint_equation_85(constraint_registration):
    """Equality constraint for the centerpost (CP) lifetime

    Depending on the chosen option i_cp_lifetime:
    - 0 : The CP full power year lifelime is set by the user (cplife_input)
    - 1 : The CP lifelime is equal to the divertor one
    - 2 : The CP lifetime is equal to the breeding blankets one
    - 3 : The CP lifetime is equal to the plant one

    cplife: calculated CP full power year lifetime (years)
    life_blkt_fpy: calculated first wall/blanket power year lifetime (years)
    life_div_fpy: calculated divertor  power year lifetime (years)
    i_cp_lifetime: switch chosing which plant element the CP
        the CP lifetime must equate
    """
    # The CP lifetime is equal to the the divertor one
    if data_structure.cost_variables.i_cp_lifetime == 0:
        bound = data_structure.cost_variables.cplife_input

    elif data_structure.cost_variables.i_cp_lifetime == 1:
        bound = data_structure.cost_variables.life_div_fpy

    # The CP lifetime is equal to the tritium breeding blankets / FW one
    elif data_structure.cost_variables.i_cp_lifetime == 2:
        bound = data_structure.fwbs_variables.life_blkt_fpy

    elif data_structure.cost_variables.i_cp_lifetime == 3:
        bound = data_structure.cost_variables.life_plant

    return eq(data_structure.cost_variables.cplife, bound, constraint_registration)


@ConstraintManager.register_constraint(86, "m", "<=")
def constraint_equation_86(constraint_registration):
    """Upper limit on the turn edge length in the TF winding pack

    dx_tf_turn_general: TF coil turn edge length including turn insulation [m]
    t_turn_tf_max: TF turn edge length including turn insulation upper limit [m]
    """
    return leq(
        data_structure.tfcoil_variables.dx_tf_turn_general,
        data_structure.tfcoil_variables.t_turn_tf_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(87, "MW", "<=")
def constraint_equation_87(constraint_registration):
    """Equation for TF coil cryogenic power upper limit

    p_cryo_plant_electric_mw: cryogenic plant power (MW)
    p_cryo_plant_electric_max_mw: Maximum cryogenic plant power (MW)
    """
    return leq(
        data_structure.heat_transport_variables.p_cryo_plant_electric_mw,
        data_structure.heat_transport_variables.p_cryo_plant_electric_max_mw,
        constraint_registration,
    )


@ConstraintManager.register_constraint(88, "", "<=")
def constraint_equation_88(constraint_registration):
    """Equation for TF coil vertical strain upper limit (absolute value)

    str_wp_max: Allowable maximum TF coil vertical strain
    str_wp: Constrained TF coil vertical strain
    """
    return leq(
        abs(data_structure.tfcoil_variables.str_wp),
        data_structure.tfcoil_variables.str_wp_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(89, "A/m2", "<=")
def constraint_equation_89(constraint_registration):
    """Upper limit to ensure that the Central Solenoid [OH] coil current / copper area < Maximum value

    copperaoh_m2: CS coil current at EOF / copper area [A/m2]
    copperaoh_m2_max: maximum coil current / copper area [A/m2]
    """
    return leq(
        data_structure.rebco_variables.copperaoh_m2,
        data_structure.rebco_variables.copperaoh_m2_max,
        constraint_registration,
    )


@ConstraintManager.register_constraint(90, "", ">=")
def constraint_equation_90(constraint_registration):
    """Lower limit for CS coil stress load cycles

    n_cycle: Allowable number of cycles for CS
    n_cycle_min: Minimum required cycles for CS
    """
    if (
        data_structure.cost_variables.ibkt_life == 1
        and data_structure.cs_fatigue_variables.bkt_life_csf == 1
    ):
        data_structure.cs_fatigue_variables.n_cycle_min = (
            data_structure.cost_variables.bktcycles
        )

    return geq(
        data_structure.cs_fatigue_variables.n_cycle,
        data_structure.cs_fatigue_variables.n_cycle_min,
        constraint_registration,
    )


@ConstraintManager.register_constraint(91, "MW", ">=")
def constraint_equation_91(constraint_registration):
    """Lower limit to ensure ECRH te is greater than required te for ignition
    at lower values for n and B. Or if the design point is ECRH heatable (if i_plasma_ignited==0)
    stellarators only (but in principle usable also for tokamaks).

    max_gyrotron_frequency: Max. av. gyrotron frequency
    te0_ecrh_achievable: Max. achievable electron temperature at ignition point
    """
    # Achievable ECRH te needs to be larger than needed te for igntion
    if data_structure.physics_variables.i_plasma_ignited == 0:
        value = (
            data_structure.stellarator_variables.powerht_constraint
            + data_structure.current_drive_variables.p_hcd_primary_extra_heat_mw
        )
    else:
        value = data_structure.stellarator_variables.powerht_constraint

    return geq(
        value,
        data_structure.stellarator_variables.powerscaling_constraint,
        constraint_registration,
    )


@ConstraintManager.register_constraint(92, "", "=")
def constraint_equation_92(constraint_registration):
    """Equation for checking is D/T ratio is consistent, and sums to 1.

    f_plasma_fuel_deuterium: fraction of deuterium ions
    f_plasma_fuel_tritium: fraction of tritium ions
    f_plasma_fuel_helium3: fraction of helium-3 ions
    """
    return eq(
        data_structure.physics_variables.f_plasma_fuel_deuterium
        + data_structure.physics_variables.f_plasma_fuel_tritium
        + data_structure.physics_variables.f_plasma_fuel_helium3,
        1.0,
        constraint_registration,
    )


def constraint_eqns(m: int, ieqn: int):
    """Evaluates the constraints given the current state of PROCESS.

    Parameters
    ----------
    m :
        The number of constraints to evaluate
    ieqn :
        Evaluates the 'ieqn'th constraint equation (index starts at 1)
        or all equations if <= 0
    m: int :

    ieqn: int :

    """

    if ieqn > 0:
        i1 = ieqn - 1
        i2 = ieqn
    else:
        i1 = 0
        i2 = m

    cc, con, err, symbol, units = [], [], [], [], []

    for i in range(i1, i2):
        constraint_id = data_structure.numerics.icc[i]
        result = ConstraintManager.evaluate_constraint(constraint_id)

        tmp_cc, tmp_con, tmp_err = (
            result.normalised_residual,
            result.constraint_value,
            result.residual,
        )
        tmp_symbol, tmp_units = result.symbol, result.units

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

    return np.array(cc), np.array(con), np.array(err), symbol, units
