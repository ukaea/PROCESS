import logging
from collections.abc import Callable, Hashable
from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np

import process.data_structure as data_structure
from process import constants
from process.exceptions import ProcessError, ProcessValueError

ConstraintSymbolType = Literal["=", ">=", "<="]

logger = logging.getLogger(__name__)


@dataclass
class ConstraintResult:
    """The constraint quantities given the current state of the code
    (aka given an evaluation at the point x).
    """

    normalised_residual: float
    """The normalised residual of the constraint."""
    constraint_value: float
    """The value of the constraint (in the physical units)."""
    constraint_error: float
    """The residual error of the constraint (in the physical units)."""


@dataclass
class ConstraintRegistration:
    """Contains the constraint equation and metadata about the constraint."""

    name: Hashable
    """The name (often a number) of the constraint. It can be any hashable e.g. a string."""
    constraint_equation: Callable[[], ConstraintResult]
    """The constraint equation that, when called, returns the normalised resiudal,
    constraint value, and constraint error.
    """
    units: str
    """The units of the constraint value and error."""
    symbol: ConstraintSymbolType
    """The type of constraint (<=, >=, ==). Only used for writing output diagnostics,
    this does not impact the calculations.
    """


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

        :param name: the name of the constraint and how it can be indexed from the registry
        :type name: Hashable
        :param units: the units of the constraint written to the output files
        :type units: str
        :param symbol: the symbol of the constraint written to the output files
        :type symbol: str
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

        :param name: the name of the constraint
        :type name: Hashable
        :returns: the constraint registration object
        :rtype: ConstraintRegistration | None
        """
        return cls._constraint_registry.get(name)

    @classmethod
    def evaluate_constraint(cls, name: Hashable):
        """Evalutes a constraint with a given name.
        :param name: the name of the constraint
        :type name: Hashable
        :returns: the result of evaluating the constraint
        :rtype: ConstraintResult | None
        """
        registration = cls.get_constraint(name)

        if registration is not None:
            return registration.constraint_equation()

        return None


@ConstraintManager.register_constraint(1, "", "=")
def constraint_equation_1():
    """Relationship between beta, temperature (keV) and density

    author: J Morris

    beta_total_vol_avg: total plasma beta
    beta_{ft}: fast alpha beta component
    beta_{NBI}: neutral beam beta component
    n_e: electron density [/m3]
    n_i: total ion density [/m3]
    T_e: density weighted average electron temperature [keV]
    T_i: density weighted average ion temperature [keV]
    B_{tot}: total toroidal + poloidal field [T]
    """

    # Density weighted temperature is used here as 〈nT〉 != 〈n〉_V * 〈T〉_V
    cc = (
        1.0
        - (
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
        )
        / data_structure.physics_variables.beta_total_vol_avg
    )
    return ConstraintResult(
        normalised_residual=cc,
        constraint_value=(
            data_structure.physics_variables.beta_total_vol_avg * (1.0 - cc)
        ),
        constraint_error=(data_structure.physics_variables.beta_total_vol_avg * cc),
    )


@ConstraintManager.register_constraint(2, "MW/m3", "=")
def constraint_equation_2():
    """author: J. Morris

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

    cc = 1.0 - pnumerator / pdenom

    return ConstraintResult(cc, pdenom * (1.0 - cc), pdenom * cc)


@ConstraintManager.register_constraint(3, "MW/m3", "=")
def constraint_equation_3():
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
        cc = 1.0 - (
            data_structure.physics_variables.pden_ion_transport_loss_mw
            + data_structure.physics_variables.pden_ion_electron_equilibration_mw
        ) / (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_ions_mw
            + data_structure.current_drive_variables.p_hcd_injected_ions_mw
            / data_structure.physics_variables.vol_plasma
        )
        return ConstraintResult(
            cc,
            (
                data_structure.physics_variables.f_p_alpha_plasma_deposited
                * data_structure.physics_variables.f_pden_alpha_ions_mw
                + data_structure.current_drive_variables.p_hcd_injected_ions_mw
                / data_structure.physics_variables.vol_plasma
            )
            * (1.0 - cc),
            (
                data_structure.physics_variables.f_p_alpha_plasma_deposited
                * data_structure.physics_variables.f_pden_alpha_ions_mw
                + data_structure.current_drive_variables.p_hcd_injected_ions_mw
                / data_structure.physics_variables.vol_plasma
            )
            * cc,
        )

    # Plasma ignited:
    cc = 1.0 - (
        data_structure.physics_variables.pden_ion_transport_loss_mw
        + data_structure.physics_variables.pden_ion_electron_equilibration_mw
    ) / (
        data_structure.physics_variables.f_p_alpha_plasma_deposited
        * data_structure.physics_variables.f_pden_alpha_ions_mw
    )
    return ConstraintResult(
        cc,
        (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_ions_mw
        )
        * (1.0 - cc),
        (
            data_structure.physics_variables.f_p_alpha_plasma_deposited
            * data_structure.physics_variables.f_pden_alpha_ions_mw
        )
        * cc,
    )


@ConstraintManager.register_constraint(4, "MW/m3", "=")
def constraint_equation_4():
    """Global power balance equation for electrons
    author: P B Lloyd, CCFE, Culham Science Centre

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

    cc = 1.0 - pnumerator / pdenom
    return ConstraintResult(cc, pdenom * (1.0 - cc), pdenom * cc)


@ConstraintManager.register_constraint(5, "/m3", "<=")
def constraint_equation_5():
    """Equation for density upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

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
    """
    # Apply Greenwald limit to line-averaged density
    if data_structure.physics_variables.i_density_limit == 7:
        return ConstraintResult(
            data_structure.physics_variables.nd_plasma_electron_line
            / data_structure.physics_variables.nd_plasma_electrons_max
            - 1.0,
            data_structure.physics_variables.nd_plasma_electrons_max,
            data_structure.physics_variables.nd_plasma_electrons_max
            - data_structure.physics_variables.nd_plasma_electron_line,
        )

    cc = (
        data_structure.physics_variables.nd_plasma_electrons_vol_avg
        / data_structure.physics_variables.nd_plasma_electrons_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.nd_plasma_electrons_max * (1.0 - cc),
        data_structure.physics_variables.nd_plasma_electrons_vol_avg * cc,
    )


@ConstraintManager.register_constraint(6, "", "<=")
def constraint_equation_6():
    """Equation for epsilon beta-poloidal upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    beta_poloidal_eps_max: maximum (eps*beta_poloidal)
    eps: inverse aspect ratio
    beta_poloidal: poloidal beta
    """
    cc = (
        data_structure.physics_variables.eps
        * data_structure.physics_variables.beta_poloidal_vol_avg
    ) / data_structure.physics_variables.beta_poloidal_eps_max - 1.0
    return ConstraintResult(
        cc,
        data_structure.physics_variables.beta_poloidal_eps_max * (1.0 - cc),
        (
            data_structure.physics_variables.eps
            * data_structure.physics_variables.beta_poloidal_vol_avg
        )
        * cc,
    )


@ConstraintManager.register_constraint(7, "/m3", "=")
def constraint_equation_7():
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

    cc = (
        1.0
        - data_structure.physics_variables.nd_beam_ions_out
        / data_structure.physics_variables.nd_beam_ions
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.nd_beam_ions * (1.0 - cc),
        data_structure.physics_variables.nd_beam_ions * cc,
    )


@ConstraintManager.register_constraint(8, "MW/m2", "<=")
def constraint_equation_8():
    """Equation for neutron wall load upper limit

    pflux_fw_neutron_max_mw: allowable wall-load (MW/m2)
    pflux_fw_neutron_mw: average neutron wall load (MW/m2)
    """
    return ConstraintResult(
        (
            data_structure.physics_variables.pflux_fw_neutron_mw
            / data_structure.constraint_variables.pflux_fw_neutron_max_mw
            - 1.0
        ),
        data_structure.constraint_variables.pflux_fw_neutron_max_mw,
        data_structure.constraint_variables.pflux_fw_neutron_max_mw
        - data_structure.physics_variables.pflux_fw_neutron_mw,
    )


@ConstraintManager.register_constraint(9, "MW", "<=")
def constraint_equation_9():
    """Equation for fusion power upper limit

    p_fusion_total_max_mw: maximum fusion power (MW)
    p_fusion_total_mw: fusion power (MW)
    """
    cc = (
        data_structure.physics_variables.p_fusion_total_mw
        / data_structure.constraint_variables.p_fusion_total_max_mw
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.p_fusion_total_max_mw * (1.0 - cc),
        data_structure.physics_variables.p_fusion_total_mw * cc,
    )


@ConstraintManager.register_constraint(11, "m", "=")
def constraint_equation_11():
    """Equation for radial build
    author: P B Lloyd, CCFE, Culham Science Centre

    rbld: sum of thicknesses to the major radius (m)
    rmajor: plasma major radius (m)
    """
    cc = (
        1.0
        - data_structure.build_variables.rbld / data_structure.physics_variables.rmajor
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.rmajor * (1.0 - cc),
        data_structure.physics_variables.rmajor * cc,
    )


@ConstraintManager.register_constraint(12, "V.sec", ">=")
def constraint_equation_12():
    """Equation for volt-second capability lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    vs_plasma_total_required: total V-s needed (Wb)
    vs_plasma_total_required (lower limit) is positive; vs_cs_pf_total_pulse (available) is negative
    vs_cs_pf_total_pulse: total flux swing for pulse (Wb)
    """
    # vs_cs_pf_total_pulse is negative, requires sign change
    cc = (
        1.0
        - (-data_structure.pfcoil_variables.vs_cs_pf_total_pulse)
        / data_structure.physics_variables.vs_plasma_total_required
    )

    return ConstraintResult(
        cc,
        data_structure.pfcoil_variables.vs_plasma_total_required * (1.0 - cc),
        data_structure.pfcoil_variables.vs_plasma_total_required * cc,
    )


@ConstraintManager.register_constraint(13, "sec", ">=")
def constraint_equation_13():
    """Equation for burn time lower limit

    author: P B Lloyd, CCFE, Culham Science Centre

    t_plant_pulse_burn: burn time (s) (calculated if i_pulsed_plant=1)
    t_burn_min: minimum burn time (s)
    """
    return ConstraintResult(
        1.0
        - data_structure.times_variables.t_plant_pulse_burn
        / data_structure.constraint_variables.t_burn_min,
        data_structure.constraint_variables.t_burn_min,
        data_structure.constraint_variables.t_burn_min
        - data_structure.times_variables.t_plant_pulse_burn,
    )


@ConstraintManager.register_constraint(15, "MW", ">=")
def constraint_equation_15():
    """Equation for L-H power threshold limit
    author: P B Lloyd, CCFE, Culham Science Centre

    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)
    """
    return ConstraintResult(
        1.0
        - data_structure.physics_variables.p_plasma_separatrix_mw
        / data_structure.physics_variables.p_l_h_threshold_mw,
        data_structure.physics_variables.p_l_h_threshold_mw,
        data_structure.physics_variables.p_l_h_threshold_mw
        - data_structure.physics_variables.p_plasma_separatrix_mw,
    )


@ConstraintManager.register_constraint(16, "MW", ">=")
def constraint_equation_16():
    """Equation for net electric power lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    p_plant_electric_net_mw: net electric power (MW)
    p_plant_electric_net_required_mw: required net electric power (MW)
    """
    return ConstraintResult(
        1.0
        - data_structure.heat_transport_variables.p_plant_electric_net_mw
        / data_structure.constraint_variables.p_plant_electric_net_required_mw,
        data_structure.constraint_variables.p_plant_electric_net_required_mw,
        data_structure.heat_transport_variables.p_plant_electric_net_mw
        - data_structure.constraint_variables.p_plant_electric_net_required_mw,
    )


@ConstraintManager.register_constraint(14, "", "=")
def constraint_equation_14():
    """Equation to fix number of NBI decay lengths to plasma centre
    author: P B Lloyd, CCFE, Culham Science Centre

    n_beam_decay_lengths_core: neutral beam e-decay lengths to plasma centre
    n_beam_decay_lengths_core_required: permitted neutral beam e-decay lengths to plasma centre
    """
    cc = (
        1.0
        - data_structure.current_drive_variables.n_beam_decay_lengths_core
        / data_structure.current_drive_variables.n_beam_decay_lengths_core_required
    )
    return ConstraintResult(
        cc,
        data_structure.current_drive_variables.n_beam_decay_lengths_core_required
        * (1.0 - cc),
        data_structure.current_drive_variables.n_beam_decay_lengths_core_required * cc,
    )


@ConstraintManager.register_constraint(17, "MW/m3", "<=")
def constraint_equation_17():
    """Equation for radiation power upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    vol_plasma: plasma volume (m3)
    pden_alpha_total_mw: alpha power per volume (MW/m3)
    pden_non_alpha_charged_mw: non-alpha charged particle fusion power per volume (MW/m3)
    pden_plasma_ohmic_mw: ohmic heating power per volume (MW/m3)
    pden_plasma_rad_mw: total radiation power per volume (MW/m3)
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

    cc = data_structure.physics_variables.pden_plasma_rad_mw / pradmaxpv - 1.0
    return ConstraintResult(
        cc,
        pradmaxpv * (1.0 - cc),
        data_structure.physics_variables.pden_plasma_rad_mw * cc,
    )


@ConstraintManager.register_constraint(18, "MW/m2", "<=")
def constraint_equation_18():
    """Equation for divertor heat load upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    pflux_div_heat_load_max_mw: heat load limit (MW/m2)
    pflux_div_heat_load_mw: divertor heat load (MW/m2)
    """
    cc = (
        data_structure.divertor_variables.pflux_div_heat_load_mw
        / data_structure.divertor_variables.pflux_div_heat_load_max_mw
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.divertor_variables.pflux_div_heat_load_max_mw * (1.0 - cc),
        data_structure.divertor_variables.pflux_div_heat_load_mw * cc,
    )


@ConstraintManager.register_constraint(19, "MVA", "<=")
def constraint_equation_19():
    """Equation for MVA (power) upper limit: resistive TF coil set
    author: P B Lloyd, CCFE, Culham Science Centre

    p_cp_resistive_mw: peak resistive TF coil inboard leg power (total) (MW)
    p_tf_leg_resistive_mw: TF coil outboard leg resistive power (total) (MW)
    mvalim: MVA limit for resistive TF coil set (total) (MW)
    """
    totmva = (
        data_structure.tfcoil_variables.p_cp_resistive_mw
        + data_structure.tfcoil_variables.p_tf_leg_resistive_mw
    )

    cc = totmva / data_structure.constraint_variables.mvalim - 1.0
    return ConstraintResult(
        cc, data_structure.constraint_variables.mvalim * (1.0 - cc), totmva * cc
    )


@ConstraintManager.register_constraint(20, "m", "<=")
def constraint_equation_20():
    """Equation for neutral beam tangency radius upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    radius_beam_tangency_max: maximum tangency radius for centreline of beam (m)
    radius_beam_tangency: neutral beam centreline tangency radius (m)
    """
    cc = (
        data_structure.current_drive_variables.radius_beam_tangency
        / data_structure.current_drive_variables.radius_beam_tangency_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.current_drive_variables.radius_beam_tangency_max * (1.0 - cc),
        data_structure.current_drive_variables.radius_beam_tangency * cc,
    )


@ConstraintManager.register_constraint(21, "", ">=")
def constraint_equation_21():
    """Equation for minor radius lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    rminor: plasma minor radius (m)
    aplasmin: minimum minor radius (m)
    """
    cc = (
        1.0
        - data_structure.physics_variables.rminor
        / data_structure.build_variables.aplasmin
    )
    return ConstraintResult(
        cc,
        data_structure.build_variables.aplasmin * (1.0 - cc),
        data_structure.build_variables.aplasmin * cc,
    )


@ConstraintManager.register_constraint(23, "m", "<=")
def constraint_equation_23():
    """Equation for conducting shell radius / rminor upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

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

    cc = (
        rcw
        / (
            data_structure.physics_variables.f_r_conducting_wall
            * data_structure.physics_variables.rminor
        )
        - 1.0
    )
    return ConstraintManager(
        cc,
        data_structure.physics_variables.f_r_conducting_wall
        * data_structure.physics_variables.rminor
        * (1.0 - cc),
        rcw * cc,
    )


@ConstraintManager.register_constraint(24, "", "<=")
def constraint_equation_24():
    """Equation for beta upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

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
        cc = (
            data_structure.physics_variables.beta_total_vol_avg
            / data_structure.physics_variables.beta_vol_avg_max
            - 1.0
        )
        con = data_structure.physics_variables.beta_vol_avg_max
        err = (
            data_structure.physics_variables.beta_vol_avg_max
            - data_structure.physics_variables.beta_total_vol_avg
        )
    # Here, the beta limit applies to only the thermal component, not the fast alpha or neutral beam parts
    elif data_structure.physics_variables.i_beta_component == 1:
        cc = (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
            - data_structure.physics_variables.beta_beam
        ) / data_structure.physics_variables.beta_vol_avg_max - 1.0
        con = data_structure.physics_variables.beta_vol_avg_max
        err = data_structure.physics_variables.beta_vol_avg_max - (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
            - data_structure.physics_variables.beta_beam
        )
    # Beta limit applies to thermal + neutral beam: components of the total beta, i.e. excludes alphas
    elif data_structure.physics_variables.i_beta_component == 2:
        cc = (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
        ) / data_structure.physics_variables.beta_vol_avg_max - 1.0
        con = data_structure.physics_variables.beta_vol_avg_max * (1.0 - cc)
        err = (
            data_structure.physics_variables.beta_total_vol_avg
            - data_structure.physics_variables.beta_fast_alpha
        ) * cc
    # Beta limit applies to toroidal beta
    elif data_structure.physics_variables.i_beta_component == 3:
        cc = (
            data_structure.physics_variables.beta_total_vol_avg
            * (
                data_structure.physics_variables.b_plasma_total
                / data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        ) / data_structure.physics_variables.beta_vol_avg_max - 1.0
        con = data_structure.physics_variables.beta_vol_avg_max
        err = data_structure.physics_variables.beta_vol_avg_max - (
            data_structure.physics_variables.beta_total_vol_avg
            * (
                data_structure.physics_variables.b_plasma_total
                / data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )

    return ConstraintResult(cc, con, err)


@ConstraintManager.register_constraint(25, "T", "<=")
def constraint_equation_25():
    """Equation for peak toroidal field upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    b_tf_inboard_max: maximum peak toroidal field (T)
    b_tf_inboard_peak_symmetric: mean peak field at TF coil (T)
    """
    cc = (
        data_structure.tfcoil_variables.b_tf_inboard_peak_symmetric
        / data_structure.constraint_variables.b_tf_inboard_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.b_tf_inboard_max * (1.0 - cc),
        data_structure.tfcoil_variables.b_tf_inboard_peak_symmetric * cc,
    )


@ConstraintManager.register_constraint(26, "A/m2", "<=")
def constraint_equation_26():
    """Equation for Central Solenoid current density upper limit at EOF
    author: P B Lloyd, CCFE, Culham Science Centre

    j_cs_critical_flat_top_end: allowable central solenoid current density at end of flat-top (A/m2)
    j_cs_flat_top_end: central solenoid overall current density at end of flat-top (A/m2)
    """
    return ConstraintResult(
        data_structure.pfcoil_variables.j_cs_flat_top_end
        / data_structure.pfcoil_variables.j_cs_critical_flat_top_end
        - 1.0,
        data_structure.pfcoil_variables.j_cs_critical_flat_top_end,
        data_structure.pfcoil_variables.j_cs_critical_flat_top_end
        - data_structure.pfcoil_variables.j_cs_flat_top_end,
    )


@ConstraintManager.register_constraint(27, "A/m2", "<=")
def constraint_equation_27():
    """Equation for Central Solenoid current density upper limit at BOP
    author: P B Lloyd, CCFE, Culham Science Centre

    j_cs_critical_pulse_start: allowable central solenoid current density at beginning of pulse (A/m2)
    j_cs_pulse_start: central solenoid overall current density at beginning of pulse (A/m2)
    """
    return ConstraintResult(
        data_structure.pfcoil_variables.j_cs_pulse_start
        / data_structure.pfcoil_variables.j_cs_critical_pulse_start
        - 1.0,
        data_structure.pfcoil_variables.j_cs_critical_pulse_start,
        data_structure.pfcoil_variables.j_cs_critical_pulse_start
        - data_structure.pfcoil_variables.j_cs_pulse_start,
    )


@ConstraintManager.register_constraint(28, "", ">=")
def constraint_equation_28():
    """Equation for fusion gain (big Q) lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

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

    cc = (
        1.0
        - data_structure.current_drive_variables.big_q_plasma
        / data_structure.constraint_variables.big_q_plasma_min
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.big_q_plasma_min * (1.0 - cc),
        data_structure.constraint_variables.big_q_plasma_min * cc,
    )


@ConstraintManager.register_constraint(29, "m", "=")
def constraint_equation_29():
    """Equation for inboard major radius: This is a consistency equation
    author: P B Lloyd, CCFE, Culham Science Centre

    rmajor: plasma major radius (m) (iteration variable 3)
    rminor: plasma minor radius (m)
    rinboard: plasma inboard radius (m)
    """
    cc = (
        1.0
        - (
            data_structure.physics_variables.rmajor
            - data_structure.physics_variables.rminor
        )
        / data_structure.build_variables.rinboard
    )
    return ConstraintResult(
        cc,
        data_structure.build_variables.rinboard * (1.0 - cc),
        data_structure.build_variables.rinboard * cc,
    )


@ConstraintManager.register_constraint(30, "MW", "<=")
def constraint_equation_30():
    """Equation for injection power upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    p_hcd_injected_max: Maximum allowable value for injected power (MW)
    """
    return ConstraintResult(
        data_structure.current_drive_variables.p_hcd_injected_total_mw
        / data_structure.current_drive_variables.p_hcd_injected_max
        - 1.0,
        data_structure.current_drive_variables.p_hcd_injected_max,
        data_structure.current_drive_variables.p_hcd_injected_max
        - data_structure.current_drive_variables.p_hcd_injected_total_mw,
    )


@ConstraintManager.register_constraint(31, "Pa", "<=")
def constraint_equation_31():
    """Equation for TF coil case stress upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    sig_tf_case_max: Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)
    sig_tf_case: Constrained stress in TF coil case (Pa)
    """
    return ConstraintResult(
        data_structure.tfcoil_variables.sig_tf_case
        / data_structure.tfcoil_variables.sig_tf_case_max
        - 1.0,
        data_structure.tfcoil_variables.sig_tf_case_max,
        data_structure.tfcoil_variables.sig_tf_case_max
        - data_structure.tfcoil_variables.sig_tf_case,
    )


@ConstraintManager.register_constraint(32, "Pa", "<=")
def constraint_equation_32():
    """Equation for TF coil conduit stress upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    sig_tf_wp_max: Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)
    sig_tf_wp: Constrained stress in TF conductor conduit (Pa)
    """
    return ConstraintResult(
        data_structure.tfcoil_variables.sig_tf_wp
        / data_structure.tfcoil_variables.sig_tf_wp_max
        - 1.0,
        data_structure.tfcoil_variables.sig_tf_wp_max,
        data_structure.tfcoil_variables.sig_tf_wp_max
        - data_structure.tfcoil_variables.sig_tf_wp,
    )


@ConstraintManager.register_constraint(33, "A/m2", "<=")
def constraint_equation_33():
    """Equation for TF coil operating/critical J upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre
    args : output structure : residual error; constraint value;

    j_tf_wp_critical: critical current density for winding pack (A/m2)
    j_tf_wp: winding pack current density (A/m2)
    """

    cc = (
        data_structure.tfcoil_variables.j_tf_wp
        / data_structure.tfcoil_variables.j_tf_wp_critical
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.j_tf_wp_critical * (1.0 - cc),
        data_structure.tfcoil_variables.j_tf_wp * cc,
    )


@ConstraintManager.register_constraint(34, "V", "<=")
def constraint_equation_34():
    """Equation for TF coil dump voltage upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    v_tf_coil_dump_quench_max_kv: max voltage across TF coil during quench (kV)
    v_tf_coil_dump_quench_kv: voltage across a TF coil during quench (kV)
    """
    return ConstraintResult(
        data_structure.tfcoil_variables.v_tf_coil_dump_quench_kv
        / data_structure.tfcoil_variables.v_tf_coil_dump_quench_max_kv
        - 1.0,
        data_structure.tfcoil_variables.v_tf_coil_dump_quench_max_kv,
        data_structure.tfcoil_variables.v_tf_coil_dump_quench_max_kv
        - data_structure.tfcoil_variables.v_tf_coil_dump_quench_kv,
    )


@ConstraintManager.register_constraint(35, "A/m2", "<=")
def constraint_equation_35():
    """Equation for TF coil J_wp/J_prot upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    j_tf_wp_quench_heat_max: allowable TF coil winding pack current density, for dump temperature
    rise protection (A/m2)
    j_tf_wp: winding pack current density (A/m2)
    """
    return ConstraintResult(
        data_structure.tfcoil_variables.j_tf_wp
        / data_structure.tfcoil_variables.j_tf_wp_quench_heat_max
        - 1.0,
        data_structure.tfcoil_variables.j_tf_wp_quench_heat_max,
        data_structure.tfcoil_variables.j_tf_wp
        - data_structure.tfcoil_variables.j_tf_wp_quench_heat_max,
    )


@ConstraintManager.register_constraint(36, "K", ">=")
def constraint_equation_36():
    """Equation for TF coil s/c temperature margin lower limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    temp_tf_superconductor_margin: TF coil temperature margin (K)
    temp_tf_superconductor_margin_min: minimum allowable temperature margin : TF coils (K)
    """
    return ConstraintResult(
        1.0
        - data_structure.tfcoil_variables.temp_tf_superconductor_margin
        / data_structure.tfcoil_variables.temp_tf_superconductor_margin_min,
        data_structure.tfcoil_variables.temp_tf_superconductor_margin_min,
        data_structure.tfcoil_variables.temp_tf_superconductor_margin_min
        - data_structure.tfcoil_variables.temp_tf_superconductor_margin,
    )


@ConstraintManager.register_constraint(37, "1E20 A/Wm2", "<=")
def constraint_equation_37():
    """Equation for current drive gamma upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    eta_cd_norm_hcd_primary_max: maximum current drive gamma
    eta_cd_norm_hcd_primary: normalised current drive efficiency (1.0e20 A/W-m2)
    """
    cc = (
        data_structure.current_drive_variables.eta_cd_norm_hcd_primary
        / data_structure.constraint_variables.eta_cd_norm_hcd_primary_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.eta_cd_norm_hcd_primary_max * (1.0 - cc),
        data_structure.current_drive_variables.eta_cd_norm_hcd_primary * cc,
    )


@ConstraintManager.register_constraint(39, "K", "<=")
def constraint_equation_39():
    """Equation for first wall temperature upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    temp_fw_max: maximum temperature of first wall material (K) (i_thermal_electric_conversion>1)
    temp_fw_peak: peak first wall temperature (K)
    """
    if data_structure.fwbs_variables.temp_fw_peak < 1.0:
        raise ProcessValueError(
            "temp_fw_peak = 0 implies i_pulsed_plant=0; do not use constraint 39 if i_pulsed_plant=0"
        )
    cc = (
        data_structure.fwbs_variables.temp_fw_peak
        / data_structure.fwbs_variables.temp_fw_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.fwbs_variables.temp_fw_max * (1.0 - cc),
        data_structure.fwbs_variables.temp_fw_peak * cc,
    )


@ConstraintManager.register_constraint(40, "MW", ">=")
def constraint_equation_40():
    """Equation for auxiliary power lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    p_hcd_injected_min_mw: minimum auxiliary power (MW)
    """
    cc = (
        1.0
        - data_structure.current_drive_variables.p_hcd_injected_total_mw
        / data_structure.constraint_variables.p_hcd_injected_min_mw
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.p_hcd_injected_min_mw * (1.0 - cc),
        data_structure.constraint_variables.p_hcd_injected_min_mw * cc,
    )


@ConstraintManager.register_constraint(41, "sec", ">=")
def constraint_equation_41():
    """Equation for plasma current ramp-up time lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    t_plant_pulse_plasma_current_ramp_up: plasma current ramp-up time for current initiation (s)
    t_current_ramp_up_min: minimum plasma current ramp-up time (s)
    """
    cc = (
        1.0
        - data_structure.times_variables.t_plant_pulse_plasma_current_ramp_up
        / data_structure.constraint_variables.t_current_ramp_up_min
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.t_current_ramp_up_min * (1.0 - cc),
        data_structure.constraint_variables.t_current_ramp_up_min * cc,
    )


@ConstraintManager.register_constraint(42, "sec", ">=")
def constraint_equation_42():
    """Equation for cycle time lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    t_plant_pulse_total: full cycle time (s)
    t_cycle_min: minimum cycle time (s)
    """
    if data_structure.constraint_variables.t_cycle_min < 1.0:
        raise ProcessValueError(
            "t_cycle_min = 0 implies that i_pulsed_plant=0; do not use constraint 42 if i_pulsed_plant=0"
        )

    cc = (
        1.0
        - data_structure.times_variables.t_plant_pulse_total
        / data_structure.constraint_variables.t_cycle_min
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.t_cycle_min * (1.0 - cc),
        data_structure.constraint_variables.t_cycle_min * cc,
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
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 43 if itart=0")

    if data_structure.tfcoil_variables.i_tf_sup == 0:
        temp_cp_average = data_structure.tfcoil_variables.temp_cp_average - 273.15
        tcpav2 = data_structure.tfcoil_variables.tcpav2 - 273.15
    else:
        temp_cp_average = data_structure.tfcoil_variables.temp_cp_average
        tcpav2 = data_structure.tfcoil_variables.tcpav2

    cc = 1.0 - temp_cp_average / tcpav2

    return ConstraintResult(cc, tcpav2 * (1.0 - cc), tcpav2 * cc)


@ConstraintManager.register_constraint(44, "deg C", "<=")
def constraint_equation_44():
    """Equation for centrepost temperature upper limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

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

    cc = temp_cp_peak / temp_cp_max - 1.0
    return ConstraintResult(cc, temp_cp_max * (1.0 - cc), temp_cp_peak * cc)


@ConstraintManager.register_constraint(45, "", ">=")
def constraint_manager_45():
    """Equation for edge safety factor lower limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

    q95 : safety factor 'near' plasma edge
    (unless i_plasma_current = 2 (ST current scaling), in which case q = mean edge safety factor qbar)
    q95_min: lower limit for edge safety factor
    itart : input integer : switch for spherical tokamak (ST) models:
    - 0 use conventional aspect ratio models;
    - 1 use spherical tokamak models"""
    if data_structure.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 45 if itart=0")

    cc = (
        1.0
        - data_structure.physics_variables.q95 / data_structure.physics_variables.q95_min
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.q95_min * (1.0 - cc),
        data_structure.physics_variables.q95_min * cc,
    )


@ConstraintManager.register_constraint(46, "", "<=")
def constraint_equation_46():
    """Equation for Ip/Irod upper limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

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
    cc = (
        data_structure.physics_variables.plasma_current
        / data_structure.tfcoil_variables.c_tf_total
    ) / cratmx - 1.0

    return ConstraintResult(
        cc,
        cratmx * (1.0 - cc),
        data_structure.physics_variables.plasma_current
        / data_structure.tfcoil_variables.c_tf_total
        * cc,
    )


@ConstraintManager.register_constraint(48, "", "<=")
def constraint_equation_48():
    """Equation for poloidal beta upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    beta_poloidal_max: maximum poloidal beta
    beta_poloidal: poloidal beta
    """
    cc = (
        data_structure.physics_variables.beta_poloidal_vol_avg
        / data_structure.constraint_variables.beta_poloidal_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.beta_poloidal_max * (1.0 - cc),
        data_structure.physics_variables.beta_poloidal_vol_avg * cc,
    )


@ConstraintManager.register_constraint(50, "Hz", "<=")
def constraint_equation_50():
    """IFE option: Equation for repetition rate upper limit
    author: P B Lloyd, CCFE, Culham Science Centre
    author: S I Muldrew, CCFE, Culham Science Centre
    """
    cc = data_structure.ife_variables.reprat / data_structure.ife_variables.rrmax - 1.0
    return ConstraintResult(
        cc,
        data_structure.ife_variables.rrmax * (1.0 - cc),
        data_structure.ife_variables.reprat * cc,
    )


@ConstraintManager.register_constraint(51, "V.s", "=")
def constraint_equation_51():
    """Equation to enforce startup flux = available startup flux
    author: P B Lloyd, CCFE, Culham Science Centre

    vs_plasma_res_ramp: resistive losses in startup V-s (Wb)
    vs_plasma_ind_ramp: internal and external plasma inductance V-s (Wb))
    vs_cs_pf_total_ramp: total flux swing for startup (Wb)
    """
    cc = 1.0 - abs(
        (
            data_structure.physics_variables.vs_plasma_res_ramp
            + data_structure.physics_variables.vs_plasma_ind_ramp
        )
        / data_structure.pfcoil_variables.vs_cs_pf_total_ramp
    )
    return ConstraintResult(
        cc,
        data_structure.pfcoil_variables.vs_cs_pf_total_ramp * (1.0 - cc),
        data_structure.pfcoil_variables.vs_cs_pf_total_ramp * cc,
    )


@ConstraintManager.register_constraint(52, "", ">=")
def constraint_equation_52():
    """Equation for tritium breeding ratio lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    tbr: tritium breeding ratio (i_blanket_type=2,3 (KIT HCPB/HCLL))
    tbrmin: minimum tritium breeding ratio (If i_blanket_type=1, tbrmin=minimum 5-year time-averaged tritium breeding ratio)
    """
    cc = (
        1.0
        - data_structure.fwbs_variables.tbr / data_structure.constraint_variables.tbrmin
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.tbrmin * (1.0 - cc),
        data_structure.constraint_variables.tbrmin * cc,
    )


@ConstraintManager.register_constraint(53, "neutron/m2", "<=")
def constraint_equation_53():
    """Equation for fast neutron fluence on TF coil upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    nflutfmax: max fast neutron fluence on TF coil (n/m2)
    nflutf: peak fast neutron fluence on TF coil superconductor (n/m2)
    """
    cc = (
        data_structure.fwbs_variables.nflutf
        / data_structure.constraint_variables.nflutfmax
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.nflutfmax * (1.0 - cc),
        data_structure.fwbs_variables.nflutf * cc,
    )


@ConstraintManager.register_constraint(54, "MW/m3", "<=")
def constraint_equation_54():
    """Equation for peak TF coil nuclear heating upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ptfnucmax: maximum nuclear heating in TF coil (MW/m3)
    ptfnucpm3: nuclear heating in the TF coil (MW/m3) (blktmodel>0)
    """
    cc = (
        data_structure.fwbs_variables.ptfnucpm3
        / data_structure.constraint_variables.ptfnucmax
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.ptfnucmax * (1.0 - cc),
        data_structure.fwbs_variables.ptfnucpm3 * cc,
    )


@ConstraintManager.register_constraint(56, "MW/m", "<=")
def constraint_equation_56():
    """Equation for power through separatrix / major radius upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    pseprmax: maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    rmajor: plasma major radius (m)
    """
    cc = (
        data_structure.physics_variables.p_plasma_separatrix_mw
        / data_structure.physics_variables.rmajor
    ) / data_structure.constraint_variables.pseprmax - 1.0
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.pseprmax * (1.0 - cc),
        (
            data_structure.physics_variables.p_plasma_separatrix_mw
            / data_structure.physics_variables.rmajor
        )
        * cc,
    )


@ConstraintManager.register_constraint(59, "", "<=")
def constraint_equation_59():
    """Equation for neutral beam shine-through fraction upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    f_p_beam_shine_through_max: maximum neutral beam shine-through fraction
    f_p_beam_shine_through: neutral beam shine-through fraction
    """
    cc = (
        data_structure.current_drive_variables.f_p_beam_shine_through
        / data_structure.constraint_variables.f_p_beam_shine_through_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.f_p_beam_shine_through_max * (1.0 - cc),
        data_structure.current_drive_variables.f_p_beam_shine_through * cc,
    )


@ConstraintManager.register_constraint(60, "K", ">=")
def constraint_equation_60():
    """Equation for Central Solenoid s/c temperature margin lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    temp_cs_superconductor_margin: Central solenoid temperature margin (K)
    temp_cs_superconductor_margin_min: Minimum allowable temperature margin : CS (K)
    """
    return ConstraintResult(
        1.0
        - data_structure.pfcoil_variables.temp_cs_superconductor_margin
        / data_structure.tfcoil_variables.temp_cs_superconductor_margin_min,
        data_structure.tfcoil_variables.temp_cs_superconductor_margin_min,
        data_structure.tfcoil_variables.temp_cs_superconductor_margin_min
        - data_structure.pfcoil_variables.temp_cs_superconductor_margin,
    )


@ConstraintManager.register_constraint(61, "", ">=")
def constraint_equation_61():
    """Equation for availability lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    f_t_plant_available: Total plant availability fraction
    avail_min: Minimum availability
    """
    cc = (
        1.0
        - data_structure.cost_variables.f_t_plant_available
        / data_structure.cost_variables.avail_min
    )
    return ConstraintResult(
        cc,
        data_structure.cost_variables.avail_min * (1.0 - cc),
        data_structure.cost_variables.f_t_plant_available * cc,
    )


@ConstraintManager.register_constraint(62, "", ">=")
def constraint_equation_62():
    """Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
    author: P B Lloyd, CCFE, Culham Science Centre

    t_alpha_confinement: alpha particle confinement time (s)
    t_energy_confinement: global thermal energy confinement time (sec)
    f_alpha_energy_confinement_min: Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
    """
    cc = (
        1.0
        - (
            data_structure.physics_variables.t_alpha_confinement
            / data_structure.physics_variables.t_energy_confinement
        )
        / data_structure.constraint_variables.f_alpha_energy_confinement_min
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.f_alpha_energy_confinement_min,
        (
            data_structure.physics_variables.t_alpha_confinement
            / data_structure.physics_variables.t_energy_confinement
        )
        * cc,
    )


@ConstraintManager.register_constraint(63, "", "<=")
def constraint_equation_63():
    """Upper limit on n_iter_vacuum_pumps (i_vacuum_pumping = simple)
    author: P B Lloyd, CCFE, Culham Science Centre

    tfno: number of TF coils (default = 50 for stellarators)
    n_iter_vacuum_pumps: number of high vacuum pumps (real number), each with the throughput
    """
    cc = (
        data_structure.vacuum_variables.n_iter_vacuum_pumps
        / data_structure.tfcoil_variables.n_tf_coils
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.n_tf_coils,
        data_structure.tfcoil_variables.n_tf_coils * cc,
    )


@ConstraintManager.register_constraint(64, "", "<=")
def constraint_equation_64():
    """Upper limit on Zeff
    author: P B Lloyd, CCFE, Culham Science Centre

    zeff_max: maximum value for Zeff
    n_charge_plasma_effective_vol_avg: plasma effective charge
    """
    cc = (
        data_structure.physics_variables.n_charge_plasma_effective_vol_avg
        / data_structure.constraint_variables.zeff_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.zeff_max,
        data_structure.constraint_variables.zeff_max * cc,
    )


@ConstraintManager.register_constraint(65, "Pa", "<=")
def constraint_equation_65():
    """Upper limit on stress of the vacuum vessel that occurs when the TF coil quenches.
    author: Timothy Nunn, UKAEA

    max_vv_stress: Maximum permitted stress of the VV (Pa)
    vv_stress_quench: Stress of the VV (Pa)
    """
    cc = (
        data_structure.superconducting_tf_coil_variables.vv_stress_quench
        / data_structure.tfcoil_variables.max_vv_stress
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.max_vv_stress,
        data_structure.tfcoil_variables.max_vv_stress * cc,
    )


@ConstraintManager.register_constraint(66, "MW", "<=")
def constrain_equation_66():
    """Upper limit on rate of change of energy in poloidal field
    author: P B Lloyd, CCFE, Culham Science Centre

    maxpoloidalpower: Maximum permitted absolute rate of change of stored energy in poloidal field (MW)
    peakpoloidalpower: Peak absolute rate of change of stored energy in poloidal field (MW) (11/01/16)
    """
    cc = (
        data_structure.pf_power_variables.peakpoloidalpower
        / data_structure.pf_power_variables.maxpoloidalpower
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.pf_power_variables.maxpoloidalpower,
        data_structure.pf_power_variables.maxpoloidalpower * cc,
    )


@ConstraintManager.register_constraint(67, "MW/m2", "<=")
def constraint_equation_67():
    """Simple upper limit on radiation wall load
    author: P B Lloyd, CCFE, Culham Science Centre

    pflux_fw_rad_max: Maximum permitted radiation wall load (MW/m^2)
    pflux_fw_rad_max_mw: Peak radiation wall load (MW/m^2)
    """
    cc = (
        data_structure.constraint_variables.pflux_fw_rad_max_mw
        / data_structure.constraint_variables.pflux_fw_rad_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.pflux_fw_rad_max,
        data_structure.constraint_variables.pflux_fw_rad_max * cc,
    )


@ConstraintManager.register_constraint(68, "MWT/m", "<=")
def constraint_equation_68():
    """Upper limit on Psep scaling (PsepB/qAR)
    author: P B Lloyd, CCFE, Culham Science Centre

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
        cc = (
            (
                data_structure.physics_variables.p_plasma_separatrix_mw
                * data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            / (
                data_structure.constraint_variables.q95_fixed
                * data_structure.physics_variables.aspect
                * data_structure.physics_variables.rmajor
            )
        ) / data_structure.constraint_variables.psepbqarmax - 1.0
        err = (
            data_structure.physics_variables.p_plasma_separatrix_mw
            * data_structure.physics_variables.b_plasma_toroidal_on_axis
        ) / (
            data_structure.constraint_variables.q95_fixed
            * data_structure.physics_variables.aspect
            * data_structure.physics_variables.rmajor
        ) - data_structure.constraint_variables.psepbqarmax
    else:
        cc = (
            (
                data_structure.physics_variables.p_plasma_separatrix_mw
                * data_structure.physics_variables.b_plasma_toroidal_on_axis
            )
            / (
                data_structure.physics_variables.q95
                * data_structure.physics_variables.aspect
                * data_structure.physics_variables.rmajor
            )
        ) / data_structure.constraint_variables.psepbqarmax - 1.0
        err = (
            data_structure.physics_variables.p_plasma_separatrix_mw
            * data_structure.physics_variables.b_plasma_toroidal_on_axis
        ) / (
            data_structure.physics_variables.q95
            * data_structure.physics_variables.aspect
            * data_structure.physics_variables.rmajor
        ) - data_structure.constraint_variables.psepbqarmax

    return ConstraintResult(cc, data_structure.constraint_variables.psepbqarmax, err)


@ConstraintManager.register_constraint(72, "Pa", "<=")
def constraint_equation_72():
    """Upper limit on central Solenoid Tresca yield stress
    author: P B Lloyd, CCFE, Culham Science Centre

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
        cc = (
            max(
                data_structure.pfcoil_variables.s_shear_cs_peak,
                data_structure.tfcoil_variables.sig_tf_cs_bucked,
            )
            / data_structure.pfcoil_variables.alstroh
            - 1.0
        )
        err = data_structure.pfcoil_variables.alstroh - max(
            data_structure.pfcoil_variables.s_shear_cs_peak,
            data_structure.tfcoil_variables.sig_tf_cs_bucked,
        )
    # Free standing CS
    else:
        cc = (
            data_structure.pfcoil_variables.s_shear_cs_peak
            / data_structure.pfcoil_variables.alstroh
            - 1.0
        )
        err = (
            data_structure.pfcoil_variables.alstroh
            - data_structure.pfcoil_variables.s_shear_cs_peak
        )

    return ConstraintResult(cc, data_structure.pfcoil_variables.alstroh, err)


@ConstraintManager.register_constraint(73, "MW", ">=")
def constraint_equation_73():
    """Lower limit to ensure separatrix power is greater than the L-H power + auxiliary power
    Related to constraint 15
    author: P B Lloyd, CCFE, Culham Science Centre

    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    p_hcd_injected_total_mw : inout real : total auxiliary injected power (MW)
    """
    cc = 1.0 - data_structure.physics_variables.p_plasma_separatrix_mw / (
        data_structure.physics_variables.p_l_h_threshold_mw
        + data_structure.current_drive_variables.p_hcd_injected_total_mw
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.p_plasma_separatrix_mw,
        data_structure.physics_variables.p_plasma_separatrix_mw * cc,
    )


@ConstraintManager.register_constraint(74, "K", "<=")
def constraint_equation_74():
    """Upper limit to ensure TF coil quench temperature < temp_croco_quench_max
    ONLY used for croco HTS coil
    author: P B Lloyd, CCFE, Culham Science Centre

    temp_croco_quench: CroCo strand: Actual temp reached during a quench (K)
    temp_croco_quench_max: CroCo strand: maximum permitted temp during a quench (K)
    """
    cc = (
        data_structure.tfcoil_variables.temp_croco_quench
        / data_structure.tfcoil_variables.temp_croco_quench_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.temp_croco_quench,
        data_structure.tfcoil_variables.temp_croco_quench * cc,
    )


@ConstraintManager.register_constraint(75, "A/m2", "<=")
def constraint_equation_75():
    """Upper limit to ensure that TF coil current / copper area < Maximum value
    ONLY used for croco HTS coil
    author: P B Lloyd, CCFE, Culham Science Centre

    copperA_m2: TF coil current / copper area (A/m2)
    copperA_m2_max: Maximum TF coil current / copper area (A/m2)
    """
    cc = (
        data_structure.rebco_variables.coppera_m2
        / data_structure.rebco_variables.coppera_m2_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.rebco_variables.coppera_m2,
        data_structure.rebco_variables.coppera_m2 * cc,
    )


@ConstraintManager.register_constraint(76, "m-3", "<=")
def constraint_equation_76():
    """Upper limit for Eich critical separatrix density model: Added for issue 558
    author: P B Lloyd, CCFE, Culham Science Centre

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

    cc = (
        data_structure.physics_variables.nd_plasma_separatrix_electron
        / data_structure.physics_variables.nd_plasma_separatrix_electron_eich_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.nd_plasma_separatrix_electron,
        data_structure.physics_variables.nd_plasma_separatrix_electron * cc,
    )


@ConstraintManager.register_constraint(77, "A/turn", "<=")
def constraint_equation_77():
    """Equation for maximum TF current per turn upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    c_tf_turn_max : allowable TF coil current per turn [A/turn]
    c_tf_turn : TF coil current per turn [A/turn]
    """
    cc = (
        data_structure.tfcoil_variables.c_tf_turn
        / data_structure.tfcoil_variables.c_tf_turn_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.c_tf_turn_max,
        data_structure.tfcoil_variables.c_tf_turn_max * cc,
    )


@ConstraintManager.register_constraint(78, "", ">=")
def constraint_equation_78():
    """Equation for Reinke criterion, divertor impurity fraction lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fzmin : input : minimum impurity fraction from Reinke model
    fzactual : input : actual impurity fraction
    """
    cc = (
        1.0
        - data_structure.reinke_variables.fzactual
        / data_structure.reinke_variables.fzmin
    )
    return ConstraintResult(
        cc,
        data_structure.reinke_variables.fzmin * (1.0 - cc),
        data_structure.reinke_variables.fzmin * cc,
    )


@ConstraintManager.register_constraint(79, "A/turn", "<=")
def constraint_equation_79():
    """Equation for maximum CS field
    author: P B Lloyd, CCFE, Culham Science Centre

    b_cs_limit_max: Central solenoid max field limit [T]
    b_cs_peak_pulse_start: maximum field in central solenoid at beginning of pulse (T)
    b_cs_peak_flat_top_end: maximum field in central solenoid at end of flat-top (EoF) (T)
    (Note: original code has "b_cs_peak_flat_top_end/b_cs_peak_pulse_start |  peak CS field [T]".)
    """
    cc = (
        max(
            data_structure.pfcoil_variables.b_cs_peak_flat_top_end,
            data_structure.pfcoil_variables.b_cs_peak_pulse_start,
        )
        / data_structure.pfcoil_variables.b_cs_limit_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.pfcoil_variables.b_cs_limit_max,
        max(
            data_structure.pfcoil_variables.b_cs_peak_flat_top_end,
            data_structure.pfcoil_variables.b_cs_peak_pulse_start,
        )
        * cc,
    )


@ConstraintManager.register_constraint(80, "MW", ">=")
def constraint_equation_80():
    """Equation for p_plasma_separatrix_mw lower limit
    author: J Morris, Culham Science Centre
    args : output structure : residual error; constraint value; residual error in physical units;
    output string; units string
    Lower limit p_plasma_separatrix_mw

    p_plasma_separatrix_min_mw : input : Minimum power crossing separatrix p_plasma_separatrix_mw [MW]
    p_plasma_separatrix_mw : input : Power crossing separatrix [MW]
    """
    cc = (
        1.0
        - data_structure.physics_variables.p_plasma_separatrix_mw
        / data_structure.constraint_variables.p_plasma_separatrix_min_mw
    )
    return ConstraintResult(
        cc,
        data_structure.constraint_variables.p_plasma_separatrix_min_mw,
        data_structure.constraint_variables.p_plasma_separatrix_min_mw * cc,
    )


@ConstraintManager.register_constraint(81, "m-3", ">=")
def constraint_equation_81():
    """Lower limit to ensure central density is larger that the pedestal one
    author: S Kahn, Culham Science Centre
    args : output structure : residual error; constraint value;
    residual error in physical units; output string; units string
    Lower limit nd_plasma_electron_on_axis > nd_plasma_pedestal_electron

    nd_plasma_electron_on_axis   : input : Central electron density [m-3]
    nd_plasma_pedestal_electron : input : Electron density at pedestal [m-3]
    """
    cc = (
        1.0
        - data_structure.physics_variables.nd_plasma_electron_on_axis
        / data_structure.physics_variables.nd_plasma_pedestal_electron
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.nd_plasma_electron_on_axis,
        (
            data_structure.physics_variables.nd_plasma_electron_on_axis
            - data_structure.physics_variables.nd_plasma_pedestal_electron
        ),
    )


@ConstraintManager.register_constraint(82, "m", ">=")
def constraint_equation_82():
    """Equation for toroidal consistency of stellarator build
    author: J Lion, IPP Greifswald

    toroidalgap: minimal gap between two stellarator coils
    dx_tf_inboard_out_toroidal: total toroidal width of a tf coil
    """
    return ConstraintResult(
        1.0
        - data_structure.tfcoil_variables.toroidalgap
        / data_structure.tfcoil_variables.dx_tf_inboard_out_toroidal,
        data_structure.tfcoil_variables.toroidalgap,
        data_structure.tfcoil_variables.toroidalgap
        - data_structure.tfcoil_variables.dx_tf_inboard_out_toroidal,
    )


@ConstraintManager.register_constraint(83, "m", ">=")
def constraint_equation_83():
    """Equation for radial consistency of stellarator build
    author: J Lion, IPP Greifswald

    available_radial_space: avaible space in radial direction as given by each s.-configuration
    required_radial_space: required space in radial direction
    """
    cc = (
        1.0
        - data_structure.build_variables.available_radial_space
        / data_structure.build_variables.required_radial_space
    )
    return ConstraintResult(
        cc,
        data_structure.build_variables.available_radial_space * (1.0 - cc),
        data_structure.build_variables.required_radial_space * cc,
    )


@ConstraintManager.register_constraint(84, "", ">=")
def constraint_equation_84():
    """Equation for the lower limit of beta
    author: J Lion, IPP Greifswald

    beta_vol_avg_min: Lower limit for beta
    beta: plasma beta
    """
    cc = (
        1.0
        - data_structure.physics_variables.beta_total_vol_avg
        / data_structure.physics_variables.beta_vol_avg_min
    )
    return ConstraintResult(
        cc,
        data_structure.physics_variables.beta_vol_avg_min * (1.0 - cc),
        data_structure.physics_variables.beta_total_vol_avg * cc,
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
    life_div_fpy: calculated divertor  power year lifetime (years)
    i_cp_lifetime: switch chosing which plant element the CP
        the CP lifetime must equate
    """
    # The CP lifetime is equal to the the divertor one
    if data_structure.cost_variables.i_cp_lifetime == 0:
        cc = (
            1.0
            - data_structure.cost_variables.cplife
            / data_structure.cost_variables.cplife_input
        )

    elif data_structure.cost_variables.i_cp_lifetime == 1:
        cc = (
            1.0
            - data_structure.cost_variables.cplife
            / data_structure.cost_variables.life_div_fpy
        )

    # The CP lifetime is equal to the tritium breeding blankets / FW one
    elif data_structure.cost_variables.i_cp_lifetime == 2:
        cc = (
            1.0
            - data_structure.cost_variables.cplife
            / data_structure.fwbs_variables.life_blkt_fpy
        )

    elif data_structure.cost_variables.i_cp_lifetime == 3:
        cc = (
            1.0
            - data_structure.cost_variables.cplife
            / data_structure.cost_variables.life_plant
        )

    return ConstraintResult(
        cc,
        data_structure.cost_variables.life_div_fpy * (1.0 - cc),
        data_structure.cost_variables.life_div_fpy * cc,
    )


@ConstraintManager.register_constraint(86, "m", "<=")
def constraint_equation_86():
    """Upper limit on the turn edge length in the TF winding pack
    Author : S Kahn

    dx_tf_turn_general: TF coil turn edge length including turn insulation [m]
    t_turn_tf_max: TF turn edge length including turn insulation upper limit [m]
    """
    cc = (
        data_structure.tfcoil_variables.dx_tf_turn_general
        / data_structure.tfcoil_variables.t_turn_tf_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.tfcoil_variables.t_turn_tf_max * (1.0 - cc),
        data_structure.tfcoil_variables.t_turn_tf_max * cc,
    )


@ConstraintManager.register_constraint(87, "MW", "<=")
def constraint_equation_87():
    """Equation for TF coil cryogenic power upper limit
    author: S. Kahn, CCFE, Culham Science Centre

    p_cryo_plant_electric_mw: cryogenic plant power (MW)
    p_cryo_plant_electric_max_mw: Maximum cryogenic plant power (MW)
    """
    cc = (
        data_structure.heat_transport_variables.p_cryo_plant_electric_mw
        / data_structure.heat_transport_variables.p_cryo_plant_electric_max_mw
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.heat_transport_variables.p_cryo_plant_electric_max_mw
        * (1.0 - cc),
        data_structure.heat_transport_variables.p_cryo_plant_electric_mw * cc,
    )


@ConstraintManager.register_constraint(88, "", "<=")
def constraint_equation_88():
    """Equation for TF coil vertical strain upper limit (absolute value)
    author: CPS Swanson, PPPL, USA

    str_wp_max: Allowable maximum TF coil vertical strain
    str_wp: Constrained TF coil vertical strain
    """
    return ConstraintResult(
        abs(data_structure.tfcoil_variables.str_wp)
        / data_structure.tfcoil_variables.str_wp_max
        - 1.0,
        data_structure.tfcoil_variables.str_wp_max,
        data_structure.tfcoil_variables.str_wp_max
        - abs(data_structure.tfcoil_variables.str_wp),
    )


@ConstraintManager.register_constraint(89, "A/m2", "<=")
def constraint_equation_89():
    """Upper limit to ensure that the Central Solenoid [OH] coil current / copper area < Maximum value
    author: G Turkington, CCFE, Culham Science Centre

    copperaoh_m2: CS coil current at EOF / copper area [A/m2]
    copperaoh_m2_max: maximum coil current / copper area [A/m2]
    """
    cc = (
        data_structure.rebco_variables.copperaoh_m2
        / data_structure.rebco_variables.copperaoh_m2_max
        - 1.0
    )
    return ConstraintResult(
        cc,
        data_structure.rebco_variables.copperaoh_m2,
        data_structure.rebco_variables.copperaoh_m2 * cc,
    )


@ConstraintManager.register_constraint(90, "", ">=")
def constraint_equation_90():
    """Lower limit for CS coil stress load cycles
    author: A. Pearce, G Turkington CCFE, Culham Science Centre

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

    cc = (
        1.0
        - data_structure.cs_fatigue_variables.n_cycle
        / data_structure.cs_fatigue_variables.n_cycle_min
    )
    return ConstraintResult(
        cc,
        data_structure.cs_fatigue_variables.n_cycle_min * (1.0 - cc),
        data_structure.cs_fatigue_variables.n_cycle * cc,
    )


@ConstraintManager.register_constraint(91, "MW", ">=")
def constraint_equation_91():
    """Lower limit to ensure ECRH te is greater than required te for ignition
    at lower values for n and B. Or if the design point is ECRH heatable (if i_plasma_ignited==0)
    stellarators only (but in principle usable also for tokamaks).
    author: J Lion, IPP Greifswald

    max_gyrotron_frequency: Max. av. gyrotron frequency
    te0_ecrh_achievable: Max. achievable electron temperature at ignition point
    """
    # Achievable ECRH te needs to be larger than needed te for igntion
    if data_structure.physics_variables.i_plasma_ignited == 0:
        cc = (
            1.0
            - (
                data_structure.stellarator_variables.powerht_constraint
                + data_structure.current_drive_variables.p_hcd_primary_extra_heat_mw
            )
            / data_structure.stellarator_variables.powerscaling_constraint
        )
    else:
        cc = (
            1.0
            - data_structure.stellarator_variables.powerht_constraint
            / data_structure.stellarator_variables.powerscaling_constraint
        )

    return ConstraintResult(
        cc,
        data_structure.stellarator_variables.powerscaling_constraint * (1.0 - cc),
        data_structure.stellarator_variables.powerht_constraint * cc,
    )


@ConstraintManager.register_constraint(92, "", "=")
def constraint_equation_92():
    """Equation for checking is D/T ratio is consistent, and sums to 1.
    author: G Turkington, UKAEA

    f_plasma_fuel_deuterium: fraction of deuterium ions
    f_plasma_fuel_tritium: fraction of tritium ions
    f_plasma_fuel_helium3: fraction of helium-3 ions
    """
    f_plasma_fuel_deuterium = 1.0 - (
        data_structure.physics_variables.f_plasma_fuel_tritium
        + data_structure.physics_variables.f_plasma_fuel_helium3
    )
    cc = 1.0 - (
        f_plasma_fuel_deuterium
        + data_structure.physics_variables.f_plasma_fuel_tritium
        + data_structure.physics_variables.f_plasma_fuel_helium3
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
        constraint_id = data_structure.numerics.icc[i]
        constraint = ConstraintManager.get_constraint(constraint_id)

        if constraint is None:
            error_msg = f"Constraint equation {constraint_id} cannot be found"
            raise ProcessError(error_msg)

        result = constraint.constraint_equation()
        tmp_cc, tmp_con, tmp_err = (
            result.normalised_residual,
            result.constraint_value,
            result.constraint_error,
        )
        tmp_symbol, tmp_units = constraint.symbol, constraint.units

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
