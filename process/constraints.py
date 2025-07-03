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
        normalised_residual=cc,
        constraint_value=(fortran.physics_variables.beta * (1.0 - cc)),
        constraint_error=(fortran.physics_variables.beta * cc),
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
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.pden_alpha_total_mw
            + fortran.physics_variables.pden_non_alpha_charged_mw
            + fortran.physics_variables.pden_plasma_ohmic_mw
            + fortran.current_drive_variables.p_hcd_injected_total_mw
            / fortran.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.pden_alpha_total_mw
            + fortran.physics_variables.pden_non_alpha_charged_mw
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
    pden_ion_electron_equilibration_mw: ion/electron equilibration power per volume (MW/m3)
    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    f_pden_alpha_ions_mw: alpha power per volume to ions (MW/m3)
    p_hcd_injected_ions_mw: auxiliary injected power to ions (MW)
    vol_plasma: plasma volume (m3)
    """
    # No assume plasma ignition:
    if fortran.physics_variables.i_plasma_ignited == 0:
        cc = 1.0 - (
            fortran.physics_variables.pden_ion_transport_loss_mw
            + fortran.physics_variables.pden_ion_electron_equilibration_mw
        ) / (
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.f_pden_alpha_ions_mw
            + fortran.current_drive_variables.p_hcd_injected_ions_mw
            / fortran.physics_variables.vol_plasma
        )
        return ConstraintResult(
            cc,
            (
                fortran.physics_variables.f_p_alpha_plasma_deposited
                * fortran.physics_variables.f_pden_alpha_ions_mw
                + fortran.current_drive_variables.p_hcd_injected_ions_mw
                / fortran.physics_variables.vol_plasma
            )
            * (1.0 - cc),
            (
                fortran.physics_variables.f_p_alpha_plasma_deposited
                * fortran.physics_variables.f_pden_alpha_ions_mw
                + fortran.current_drive_variables.p_hcd_injected_ions_mw
                / fortran.physics_variables.vol_plasma
            )
            * cc,
        )

    # Plasma ignited:
    cc = 1.0 - (
        fortran.physics_variables.pden_ion_transport_loss_mw
        + fortran.physics_variables.pden_ion_electron_equilibration_mw
    ) / (
        fortran.physics_variables.f_p_alpha_plasma_deposited
        * fortran.physics_variables.f_pden_alpha_ions_mw
    )
    return ConstraintResult(
        cc,
        (
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.f_pden_alpha_ions_mw
        )
        * (1.0 - cc),
        (
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.f_pden_alpha_ions_mw
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
    f_p_alpha_plasma_deposited: fraction of alpha power deposited in plasma
    f_pden_alpha_electron_mw: alpha power per volume to electrons (MW/m3)
    pden_ion_electron_equilibration_mw: ion/electron equilibration power per volume (MW/m3)
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
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.f_pden_alpha_electron_mw
            + fortran.physics_variables.pden_ion_electron_equilibration_mw
            + fortran.current_drive_variables.p_hcd_injected_electrons_mw
            / fortran.physics_variables.vol_plasma
        )
    else:
        # if plasma ignited
        pdenom = (
            fortran.physics_variables.f_p_alpha_plasma_deposited
            * fortran.physics_variables.f_pden_alpha_electron_mw
            + fortran.physics_variables.pden_ion_electron_equilibration_mw
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
    nd_electron_line: line averaged electron density (m-3)

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
            fortran.physics_variables.nd_electron_line
            / fortran.physics_variables.dnelimt
            - 1.0 * fortran.constraint_variables.fdene,
            fortran.constraint_variables.fdene * fortran.physics_variables.dnelimt,
            fortran.constraint_variables.fdene * fortran.physics_variables.dnelimt
            - fortran.physics_variables.nd_electron_line,
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


@ConstraintManager.register_constraint(6, "", "<=")
def constraint_equation_6():
    """Equation for epsilon beta-poloidal upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fbeta_poloidal_eps: f-value for epsilon beta-poloidal
    beta_poloidal_eps_max: maximum (eps*beta_poloidal)
    eps: inverse aspect ratio
    beta_poloidal: poloidal beta
    """
    cc = (
        (fortran.physics_variables.eps * fortran.physics_variables.beta_poloidal)
        / fortran.physics_variables.beta_poloidal_eps_max
        - 1.0 * fortran.constraint_variables.fbeta_poloidal_eps
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.beta_poloidal_eps_max * (1.0 - cc),
        (fortran.physics_variables.eps * fortran.physics_variables.beta_poloidal) * cc,
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
    nd_beam_ions_out: hot beam ion density from calculation (/m3)
    nd_beam_ions: hot beam ion density, variable (/m3)
    """
    if fortran.physics_variables.i_plasma_ignited == 1:
        raise ProcessValueError(
            "Do not use constraint equation 7 if i_plasma_ignited=1"
        )

    cc = (
        1.0
        - fortran.physics_variables.nd_beam_ions_out
        / fortran.physics_variables.nd_beam_ions
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.nd_beam_ions * (1.0 - cc),
        fortran.physics_variables.nd_beam_ions * cc,
    )


@ConstraintManager.register_constraint(8, "MW/m2", "<=")
def constraint_equation_8():
    """Equation for neutron wall load upper limit

    fpflux_fw_neutron_max_mw: f-value for maximum wall load
    pflux_fw_neutron_max_mw: allowable wall-load (MW/m2)
    pflux_fw_neutron_mw: average neutron wall load (MW/m2)
    """
    return ConstraintResult(
        (
            fortran.physics_variables.pflux_fw_neutron_mw
            / fortran.constraint_variables.pflux_fw_neutron_max_mw
            - 1.0 * fortran.constraint_variables.fpflux_fw_neutron_max_mw
        ),
        fortran.constraint_variables.fpflux_fw_neutron_max_mw
        * fortran.constraint_variables.pflux_fw_neutron_max_mw,
        fortran.constraint_variables.fpflux_fw_neutron_max_mw
        * fortran.constraint_variables.pflux_fw_neutron_max_mw
        - fortran.physics_variables.pflux_fw_neutron_mw,
    )


@ConstraintManager.register_constraint(9, "MW", "<=")
def constraint_equation_9():
    """Equation for fusion power upper limit

    fp_fusion_total_max_mw: f-value for maximum fusion power
    p_fusion_total_max_mw: maximum fusion power (MW)
    p_fusion_total_mw: fusion power (MW)
    """
    cc = (
        fortran.physics_variables.p_fusion_total_mw
        / fortran.constraint_variables.p_fusion_total_max_mw
        - 1.0 * fortran.constraint_variables.fp_fusion_total_max_mw
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.p_fusion_total_max_mw * (1.0 - cc),
        fortran.physics_variables.p_fusion_total_mw * cc,
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
    fvs_plasma_total_required: f-value for flux-swing (V-s) requirement (STEADY STATE)
    vs_cs_pf_total_pulse: total flux swing for pulse (Wb)
    """
    # vs_cs_pf_total_pulse is negative, requires sign change
    cc = (
        1.0
        - fortran.constraint_variables.fvs_plasma_total_required
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

    ft_burn_min: f-value for minimum burn time
    t_burn: burn time (s) (calculated if i_pulsed_plant=1)
    t_burn_min: minimum burn time (s)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.ft_burn_min
        * fortran.times_variables.t_burn
        / fortran.constraint_variables.t_burn_min,
        fortran.constraint_variables.t_burn_min
        / fortran.constraint_variables.ft_burn_min,
        fortran.constraint_variables.t_burn_min
        / fortran.constraint_variables.ft_burn_min
        - fortran.times_variables.t_burn,
    )


@ConstraintManager.register_constraint(15, "MW", ">=")
def constraint_equation_15():
    """Equation for L-H power threshold limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fl_h_threshold: f-value for L-H power threshold
    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.fl_h_threshold
        * fortran.physics_variables.p_plasma_separatrix_mw
        / fortran.physics_variables.p_l_h_threshold_mw,
        fortran.physics_variables.p_l_h_threshold_mw,
        fortran.physics_variables.p_l_h_threshold_mw
        - fortran.physics_variables.p_plasma_separatrix_mw
        / fortran.constraint_variables.fl_h_threshold,
    )


@ConstraintManager.register_constraint(16, "MW", ">=")
def constraint_equation_16():
    """Equation for net electric power lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fp_plant_electric_net_required_mw: f-value for net electric power
    p_plant_electric_net_mw: net electric power (MW)
    p_plant_electric_net_required_mw: required net electric power (MW)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.fp_plant_electric_net_required_mw
        * fortran.heat_transport_variables.p_plant_electric_net_mw
        / fortran.constraint_variables.p_plant_electric_net_required_mw,
        fortran.constraint_variables.p_plant_electric_net_required_mw,
        fortran.heat_transport_variables.p_plant_electric_net_mw
        - fortran.constraint_variables.p_plant_electric_net_required_mw
        / fortran.constraint_variables.fp_plant_electric_net_required_mw,
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
        - fortran.current_drive_variables.n_beam_decay_lengths_core
        / fortran.current_drive_variables.n_beam_decay_lengths_core_required
    )
    return ConstraintResult(
        cc,
        fortran.current_drive_variables.n_beam_decay_lengths_core_required * (1.0 - cc),
        fortran.current_drive_variables.n_beam_decay_lengths_core_required * cc,
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
    fradpwr: f-value for core radiation power limit
    pden_plasma_rad_mw: total radiation power per volume (MW/m3)
    """
    # Maximum possible power/vol_plasma that can be radiated (local)
    pradmaxpv = (
        fortran.current_drive_variables.p_hcd_injected_total_mw
        / fortran.physics_variables.vol_plasma
        + fortran.physics_variables.pden_alpha_total_mw
        * fortran.physics_variables.f_p_alpha_plasma_deposited
        + fortran.physics_variables.pden_non_alpha_charged_mw
        + fortran.physics_variables.pden_plasma_ohmic_mw
    )

    cc = (
        fortran.physics_variables.pden_plasma_rad_mw / pradmaxpv
        - 1.0 * fortran.constraint_variables.fradpwr
    )
    return ConstraintResult(
        cc, pradmaxpv * (1.0 - cc), fortran.physics_variables.pden_plasma_rad_mw * cc
    )


@ConstraintManager.register_constraint(18, "MW/m2", "<=")
def constraint_equation_18():
    """Equation for divertor heat load upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fpflux_div_heat_load_mw: f-value for divertor heat load
    pflux_div_heat_load_max_mw: heat load limit (MW/m2)
    pflux_div_heat_load_mw: divertor heat load (MW/m2)
    """
    cc = (
        fortran.divertor_variables.pflux_div_heat_load_mw
        / fortran.divertor_variables.pflux_div_heat_load_max_mw
        - 1.0 * fortran.constraint_variables.fpflux_div_heat_load_mw
    )
    return ConstraintResult(
        cc,
        fortran.divertor_variables.pflux_div_heat_load_max_mw * (1.0 - cc),
        fortran.divertor_variables.pflux_div_heat_load_mw * cc,
    )


@ConstraintManager.register_constraint(19, "MVA", "<=")
def constraint_equation_19():
    """Equation for MVA (power) upper limit: resistive TF coil set
    author: P B Lloyd, CCFE, Culham Science Centre

    p_cp_resistive_mw: peak resistive TF coil inboard leg power (total) (MW)
    p_tf_leg_resistive_mw: TF coil outboard leg resistive power (total) (MW)
    fmva: f-value for maximum MVA
    mvalim: MVA limit for resistive TF coil set (total) (MW)
    """
    totmva = (
        fortran.tfcoil_variables.p_cp_resistive_mw
        + fortran.tfcoil_variables.p_tf_leg_resistive_mw
    )

    cc = (
        totmva / fortran.constraint_variables.mvalim
        - 1.0 * fortran.constraint_variables.fmva
    )
    return ConstraintManager(
        cc, fortran.constraint_variables.mvalim * (1.0 - cc), totmva * cc
    )


@ConstraintManager.register_constraint(20, "m", "<=")
def constraint_equation_20():
    """Equation for neutral beam tangency radius upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fradius_beam_tangency: f-value for neutral beam tangency radius limit
    radius_beam_tangency_max: maximum tangency radius for centreline of beam (m)
    radius_beam_tangency: neutral beam centreline tangency radius (m)
    """
    cc = (
        fortran.current_drive_variables.radius_beam_tangency
        / fortran.current_drive_variables.radius_beam_tangency_max
        - 1.0 * fortran.constraint_variables.fradius_beam_tangency
    )
    return ConstraintResult(
        cc,
        fortran.current_drive_variables.radius_beam_tangency_max * (1.0 - cc),
        fortran.current_drive_variables.radius_beam_tangency * cc,
    )


@ConstraintManager.register_constraint(21, "", ">=")
def constraint_equation_21():
    """Equation for minor radius lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    frminor: f-value for minor radius limit
    rminor: plasma minor radius (m)
    aplasmin: minimum minor radius (m)
    """
    cc = (
        1.0
        - fortran.constraint_variables.frminor
        * fortran.physics_variables.rminor
        / fortran.build_variables.aplasmin
    )
    return ConstraintResult(
        cc,
        fortran.build_variables.aplasmin * (1.0 - cc),
        fortran.build_variables.aplasmin * cc,
    )


@ConstraintManager.register_constraint(23, "m", "<=")
def constraint_equation_23():
    """Equation for conducting shell radius / rminor upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    rminor: plasma minor radius (m)
    dr_fw_plasma_gap_outboard: gap between plasma and first wall, outboard side (m)
    dr_fw_outboard: outboard first wall thickness, initial estimate (m)
    dr_blkt_outboard: outboard blanket thickness (m)
    fr_conducting_wall: f-value for conducting wall radius / rminor limit
    f_r_conducting_wall: maximum ratio of conducting wall distance to plasma minor radius for vertical stability
    """
    # conducting shell radius (m)
    rcw = (
        fortran.physics_variables.rminor
        + fortran.build_variables.dr_fw_plasma_gap_outboard
        + fortran.build_variables.dr_fw_outboard
        + fortran.build_variables.dr_blkt_outboard
    )

    cc = (
        rcw
        / (
            fortran.physics_variables.f_r_conducting_wall
            * fortran.physics_variables.rminor
        )
        - 1.0 * fortran.constraint_variables.fr_conducting_wall
    )
    return ConstraintManager(
        cc,
        fortran.physics_variables.f_r_conducting_wall
        * fortran.physics_variables.rminor
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
    istell: switch for stellarator option (set via <CODE>device.dat</CODE>):
    - 0 use tokamak model;
    - 1 use stellarator model
    fbeta_max: f-value for beta limit
    beta_max: allowable beta
    beta: total plasma beta (calculated if ipedestal =3)
    beta_fast_alpha: fast alpha beta component
    beta_beam: neutral beam beta component
    bt: toroidal field
    btot: total field
    """
    # Include all beta components: relevant for both tokamaks and stellarators
    if (
        fortran.physics_variables.i_beta_component == 0
        or fortran.stellarator_variables.istell != 0
    ):
        cc = (
            fortran.physics_variables.beta / fortran.physics_variables.beta_max
            - 1.0 * fortran.constraint_variables.fbeta_max
        )
        con = fortran.physics_variables.beta_max
        err = (
            fortran.physics_variables.beta_max
            - fortran.physics_variables.beta / fortran.constraint_variables.fbeta_max
        )
    # Here, the beta limit applies to only the thermal component, not the fast alpha or neutral beam parts
    elif fortran.physics_variables.i_beta_component == 1:
        cc = (
            (
                fortran.physics_variables.beta
                - fortran.physics_variables.beta_fast_alpha
                - fortran.physics_variables.beta_beam
            )
            / fortran.physics_variables.beta_max
            - 1.0 * fortran.constraint_variables.fbeta_max
        )
        con = fortran.physics_variables.beta_max
        err = (
            fortran.physics_variables.beta_max
            - (
                fortran.physics_variables.beta
                - fortran.physics_variables.beta_fast_alpha
                - fortran.physics_variables.beta_beam
            )
            / fortran.constraint_variables.fbeta_max
        )
    # Beta limit applies to thermal + neutral beam: components of the total beta, i.e. excludes alphas
    elif fortran.physics_variables.i_beta_component == 2:
        cc = (
            (fortran.physics_variables.beta - fortran.physics_variables.beta_fast_alpha)
            / fortran.physics_variables.beta_max
            - 1.0 * fortran.constraint_variables.fbeta_max
        )
        con = fortran.physics_variables.beta_max * (1.0 - cc)
        err = (
            fortran.physics_variables.beta - fortran.physics_variables.beta_fast_alpha
        ) * cc
    # Beta limit applies to toroidal beta
    elif fortran.physics_variables.i_beta_component == 3:
        cc = (
            (
                fortran.physics_variables.beta
                * (fortran.physics_variables.btot / fortran.physics_variables.bt) ** 2
            )
            / fortran.physics_variables.beta_max
            - 1.0 * fortran.constraint_variables.fbeta_max
        )
        con = fortran.physics_variables.beta_max
        err = (
            fortran.physics_variables.beta_max
            - (
                fortran.physics_variables.beta
                * (fortran.physics_variables.btot / fortran.physics_variables.bt) ** 2
            )
            / fortran.constraint_variables.fbeta_max
        )

    return ConstraintResult(cc, con, err)


@ConstraintManager.register_constraint(25, "T", "<=")
def constraint_equation_25():
    """Equation for peak toroidal field upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fb_tf_inboard_max: f-value for maximum toroidal field
    b_tf_inboard_max: maximum peak toroidal field (T)
    b_tf_inboard_peak: mean peak field at TF coil (T)
    """
    cc = (
        fortran.tfcoil_variables.b_tf_inboard_peak
        / fortran.constraint_variables.b_tf_inboard_max
        - 1.0 * fortran.constraint_variables.fb_tf_inboard_max
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.b_tf_inboard_max * (1.0 - cc),
        fortran.tfcoil_variables.b_tf_inboard_peak * cc,
    )


@ConstraintManager.register_constraint(26, "A/m2", "<=")
def constraint_equation_26():
    """Equation for Central Solenoid current density upper limit at EOF
    author: P B Lloyd, CCFE, Culham Science Centre

    fjohcreal: f-value for central solenoid current at end-of-flattop
    j_cs_critical_flat_top_end: allowable central solenoid current density at end of flat-top (A/m2)
    j_cs_flat_top_end: central solenoid overall current density at end of flat-top (A/m2)
    """
    return ConstraintResult(
        fortran.pfcoil_variables.j_cs_flat_top_end
        / fortran.pfcoil_variables.j_cs_critical_flat_top_end
        - 1.0 * fortran.constraint_variables.fjohc,
        fortran.pfcoil_variables.j_cs_critical_flat_top_end,
        fortran.pfcoil_variables.j_cs_critical_flat_top_end
        - fortran.pfcoil_variables.j_cs_flat_top_end
        / fortran.constraint_variables.fjohc,
    )


@ConstraintManager.register_constraint(27, "A/m2", "<=")
def constraint_equation_27():
    """Equation for Central Solenoid current density upper limit at BOP
    author: P B Lloyd, CCFE, Culham Science Centre

    fjohc0: f-value for central solenoid current at beginning of pulse
    j_cs_critical_pulse_start: allowable central solenoid current density at beginning of pulse (A/m2)
    j_cs_pulse_start: central solenoid overall current density at beginning of pulse (A/m2)
    """
    return ConstraintResult(
        fortran.pfcoil_variables.j_cs_pulse_start
        / fortran.pfcoil_variables.j_cs_critical_pulse_start
        - 1.0 * fortran.constraint_variables.fjohc0,
        fortran.pfcoil_variables.j_cs_critical_pulse_start,
        fortran.pfcoil_variables.j_cs_critical_pulse_start
        - fortran.pfcoil_variables.j_cs_pulse_start
        / fortran.constraint_variables.fjohc0,
    )


@ConstraintManager.register_constraint(28, "", ">=")
def constraint_equation_28():
    """Equation for fusion gain (big Q) lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fqval: pf-value for Q
    bigq: Fusion gain; P_fusion / (P_injection + P_ohmic)
    bigqmin: minimum fusion gain Q
    i_plasma_ignited : input integer : switch for ignition assumption:
    - 0 do not assume plasma ignition;
    - 1 assume ignited (but include auxiliary power in costs)
    Obviously, ignite must be zero if current drive is required.
    If i_plasma_ignited=1, any auxiliary power is assumed to be used only
    during plasma start-up, and is excluded from all steady-state
    power balance calculations.
    """
    if fortran.physics_variables.i_plasma_ignited != 0:
        raise ProcessValueError("Do not use constraint 28 if i_plasma_ignited=1")

    cc = (
        1.0
        - fortran.constraint_variables.fqval
        * fortran.current_drive_variables.bigq
        / fortran.constraint_variables.bigqmin
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.bigqmin * (1.0 - cc),
        fortran.constraint_variables.bigqmin * cc,
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
        - (fortran.physics_variables.rmajor - fortran.physics_variables.rminor)
        / fortran.build_variables.rinboard
    )
    return ConstraintResult(
        cc,
        fortran.build_variables.rinboard * (1.0 - cc),
        fortran.build_variables.rinboard * cc,
    )


@ConstraintManager.register_constraint(30, "MW", "<=")
def constraint_equation_30():
    """Equation for injection power upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    fp_hcd_injected_max: f-value for injection power
    p_hcd_injected_max: Maximum allowable value for injected power (MW)
    """
    return ConstraintResult(
        fortran.current_drive_variables.p_hcd_injected_total_mw
        / fortran.current_drive_variables.p_hcd_injected_max
        - 1.0 * fortran.constraint_variables.fp_hcd_injected_max,
        fortran.current_drive_variables.p_hcd_injected_max,
        fortran.current_drive_variables.p_hcd_injected_max
        - fortran.current_drive_variables.p_hcd_injected_total_mw
        / fortran.constraint_variables.fp_hcd_injected_max,
    )


@ConstraintManager.register_constraint(31, "Pa", "<=")
def constraint_equation_31():
    """Equation for TF coil case stress upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    fstrcase: f-value for TF coil case stress
    sig_tf_case_max: Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)
    sig_tf_case: Constrained stress in TF coil case (Pa)
    """
    return ConstraintResult(
        fortran.tfcoil_variables.sig_tf_case / fortran.tfcoil_variables.sig_tf_case_max
        - 1.0 * fortran.constraint_variables.fstrcase,
        fortran.tfcoil_variables.sig_tf_case_max,
        fortran.tfcoil_variables.sig_tf_case_max
        - fortran.tfcoil_variables.sig_tf_case / fortran.constraint_variables.fstrcase,
    )


@ConstraintManager.register_constraint(32, "Pa", "<=")
def constraint_equation_32():
    """Equation for TF coil conduit stress upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    fstrcond: f-value for TF coil conduit stress
    sig_tf_wp_max: Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)
    sig_tf_wp: Constrained stress in TF conductor conduit (Pa)
    """
    return ConstraintResult(
        fortran.tfcoil_variables.sig_tf_wp / fortran.tfcoil_variables.sig_tf_wp_max
        - 1.0 * fortran.constraint_variables.fstrcond,
        fortran.tfcoil_variables.sig_tf_wp_max,
        fortran.tfcoil_variables.sig_tf_wp_max
        - fortran.tfcoil_variables.sig_tf_wp / fortran.constraint_variables.fstrcond,
    )


@ConstraintManager.register_constraint(33, "A/m2", "<=")
def constraint_equation_33():
    """Equation for TF coil operating/critical J upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre
    args : output structure : residual error; constraint value;

    fiooic: f-value for TF coil operating current / critical
    j_tf_wp_critical: critical current density for winding pack (A/m2)
    j_tf_wp: winding pack current density (A/m2)
    """
    if fortran.constraint_variables.fiooic > 0.7:
        fortran.error_handling.report_error(285)

    cc = (
        fortran.tfcoil_variables.j_tf_wp / fortran.tfcoil_variables.j_tf_wp_critical
        - 1.0 * fortran.constraint_variables.fiooic
    )
    return ConstraintResult(
        cc,
        fortran.tfcoil_variables.j_tf_wp_critical * (1.0 - cc),
        fortran.tfcoil_variables.j_tf_wp * cc,
    )


@ConstraintManager.register_constraint(34, "V", "<=")
def constraint_equation_34():
    """Equation for TF coil dump voltage upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    fvdump: f-value for dump voltage
    vdalw: max voltage across TF coil during quench (kV)
    vtfskv: voltage across a TF coil during quench (kV)
    """
    return ConstraintResult(
        fortran.tfcoil_variables.vtfskv / fortran.tfcoil_variables.vdalw
        - 1.0 * fortran.constraint_variables.fvdump,
        fortran.tfcoil_variables.vdalw,
        fortran.tfcoil_variables.vdalw - fortran.tfcoil_variables.vtfskv,
    )


@ConstraintManager.register_constraint(35, "A/m2", "<=")
def constraint_equation_35():
    """Equation for TF coil J_wp/J_prot upper limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    fjprot: f-value for TF coil winding pack current density
    jwdgpro: allowable TF coil winding pack current density, for dump temperature
    rise protection (A/m2)
    j_tf_wp: winding pack current density (A/m2)
    """
    return ConstraintResult(
        fortran.tfcoil_variables.j_tf_wp / fortran.tfcoil_variables.jwdgpro
        - 1.0 * fortran.constraint_variables.fjprot,
        fortran.tfcoil_variables.jwdgpro,
        fortran.tfcoil_variables.j_tf_wp - fortran.tfcoil_variables.jwdgpro,
    )


@ConstraintManager.register_constraint(36, "K", ">=")
def constraint_equation_36():
    """Equation for TF coil s/c temperature margin lower limit (SCTF)
    author: P B Lloyd, CCFE, Culham Science Centre

    ftmargtf: f-value for TF coil temperature margin
    tmargtf: TF coil temperature margin (K)
    tmargmin_tf: minimum allowable temperature margin : TF coils (K)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.ftmargtf
        * fortran.tfcoil_variables.tmargtf
        / fortran.tfcoil_variables.tmargmin_tf,
        fortran.tfcoil_variables.tmargmin_tf,
        fortran.tfcoil_variables.tmargmin_tf - fortran.tfcoil_variables.tmargtf,
    )


@ConstraintManager.register_constraint(37, "1E20 A/Wm2", "<=")
def constraint_equation_37():
    """Equation for current drive gamma upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fgamcd: f-value for current drive gamma
    eta_cd_norm_hcd_primary_max: maximum current drive gamma
    eta_cd_norm_hcd_primary: normalised current drive efficiency (1.0e20 A/W-m2)
    """
    cc = (
        fortran.current_drive_variables.eta_cd_norm_hcd_primary
        / fortran.constraint_variables.eta_cd_norm_hcd_primary_max
        - 1.0 * fortran.constraint_variables.fgamcd
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.eta_cd_norm_hcd_primary_max * (1.0 - cc),
        fortran.current_drive_variables.eta_cd_norm_hcd_primary * cc,
    )


@ConstraintManager.register_constraint(39, "K", "<=")
def constraint_equation_39():
    """Equation for first wall temperature upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ftpeak: f-value for first wall peak temperature
    temp_fw_max: maximum temperature of first wall material (K) (i_thermal_electric_conversion>1)
    temp_fw_peak: peak first wall temperature (K)
    """
    if fortran.fwbs_variables.temp_fw_peak < 1.0:
        raise ProcessValueError(
            "temp_fw_peak = 0 implies i_pulsed_plant=0; do not use constraint 39 if i_pulsed_plant=0"
        )
    cc = (
        fortran.fwbs_variables.temp_fw_peak / fortran.fwbs_variables.temp_fw_max
        - 1.0 * fortran.constraint_variables.ftpeak
    )
    return ConstraintResult(
        cc,
        fortran.fwbs_variables.temp_fw_max * (1.0 - cc),
        fortran.fwbs_variables.temp_fw_peak * cc,
    )


@ConstraintManager.register_constraint(40, "MW", ">=")
def constraint_equation_40():
    """Equation for auxiliary power lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fauxmn: f-value for minimum auxiliary power
    p_hcd_injected_total_mw: total auxiliary injected power (MW)
    auxmin: minimum auxiliary power (MW)
    """
    cc = (
        1.0
        - fortran.constraint_variables.fauxmn
        * fortran.current_drive_variables.p_hcd_injected_total_mw
        / fortran.constraint_variables.auxmin
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.auxmin * (1.0 - cc),
        fortran.constraint_variables.auxmin * cc,
    )


@ConstraintManager.register_constraint(41, "sec", ">=")
def constraint_equation_41():
    """Equation for plasma current ramp-up time lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ft_current_ramp_up: f-value for plasma current ramp-up time
    t_current_ramp_up: plasma current ramp-up time for current initiation (s)
    t_current_ramp_up_min: minimum plasma current ramp-up time (s)
    """
    cc = (
        1.0
        - fortran.constraint_variables.ft_current_ramp_up
        * fortran.times_variables.t_current_ramp_up
        / fortran.constraint_variables.t_current_ramp_up_min
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.t_current_ramp_up_min * (1.0 - cc),
        fortran.constraint_variables.t_current_ramp_up_min * cc,
    )


@ConstraintManager.register_constraint(42, "sec", ">=")
def constraint_equation_42():
    """Equation for cycle time lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ftcycl: f-value for cycle time
    t_cycle: full cycle time (s)
    tcycmn: minimum cycle time (s)
    """
    if fortran.constraint_variables.tcycmn < 1.0:
        raise ProcessValueError(
            "tcycmn = 0 implies that i_pulsed_plant=0; do not use constraint 42 if i_pulsed_plant=0"
        )

    cc = (
        1.0
        - fortran.constraint_variables.ftcycl
        * fortran.times_variables.t_cycle
        / fortran.constraint_variables.tcycmn
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.tcycmn * (1.0 - cc),
        fortran.constraint_variables.tcycmn * cc,
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


@ConstraintManager.register_constraint(44, "deg C", "<=")
def constraint_equation_44():
    """Equation for centrepost temperature upper limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

    fptemp: f-value for peak centrepost temperature
    ptempalw: maximum peak centrepost temperature (K)
    tcpmax: peak centrepost temperature (K)
    itart: switch for spherical tokamak (ST) models:
    - 0: use conventional aspect ratio models;
    - 1: use spherical tokamak models
    """
    if fortran.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 44 if itart=0")

    if fortran.tfcoil_variables.i_tf_sup == 0:  # ! Copper case
        ptempalw = fortran.tfcoil_variables.ptempalw - 273.15
        tcpmax = fortran.tfcoil_variables.tcpmax - 273.15
    else:
        ptempalw = fortran.tfcoil_variables.ptempalw
        tcpmax = fortran.tfcoil_variables.tcpmax

    cc = tcpmax / ptempalw - 1.0 * fortran.constraint_variables.fptemp
    return ConstraintResult(cc, ptempalw * (1.0 - cc), tcpmax * cc)


@ConstraintManager.register_constraint(45, "", ">=")
def constraint_manager_45():
    """Equation for edge safety factor lower limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

    fq: f-value for edge safety factor
    q95 : safety factor 'near' plasma edge
    (unless i_plasma_current = 2 (ST current scaling), in which case q = mean edge safety factor qbar)
    q95_min: lower limit for edge safety factor
    itart : input integer : switch for spherical tokamak (ST) models:
    - 0 use conventional aspect ratio models;
    - 1 use spherical tokamak models"""
    if fortran.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 45 if itart=0")

    cc = (
        1.0
        - fortran.constraint_variables.fq
        * fortran.physics_variables.q95
        / fortran.physics_variables.q95_min
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.q95_min * (1.0 - cc),
        fortran.physics_variables.q95_min * cc,
    )


@ConstraintManager.register_constraint(46, "", "<=")
def constraint_equation_46():
    """Equation for Ip/Irod upper limit (TART)
    author: P B Lloyd, CCFE, Culham Science Centre

    eps: inverse aspect ratio
    fipir: f-value for Ip/Irod upper limit
    c_tf_total: total (summed) current in TF coils (A)
    plasma_current: plasma current (A)
    itart: switch for spherical tokamak (ST) models:
    - 0: use conventional aspect ratio models;
    - 1: use spherical tokamak models
    """
    if fortran.physics_variables.itart == 0:
        raise ProcessValueError("Do not use constraint 46 if itart=0")

    # maximum ratio of plasma current to centrepost current
    cratmx = 1.0 + 4.91 * (fortran.physics_variables.eps - 0.62)
    cc = (
        fortran.physics_variables.plasma_current / fortran.tfcoil_variables.c_tf_total
    ) / cratmx - 1.0 * fortran.constraint_variables.fipir

    return ConstraintResult(
        cc,
        cratmx * (1.0 - cc),
        fortran.physics_variables.plasma_current
        / fortran.tfcoil_variables.c_tf_total
        * cc,
    )


@ConstraintManager.register_constraint(48, "", "<=")
def constraint_equation_48():
    """Equation for poloidal beta upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fbeta_poloidal: rf-value for poloidal beta
    beta_poloidal_max: maximum poloidal beta
    beta_poloidal: poloidal beta
    """
    cc = (
        fortran.physics_variables.beta_poloidal
        / fortran.constraint_variables.beta_poloidal_max
        - 1.0 * fortran.constraint_variables.fbeta_poloidal
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.beta_poloidal_max * (1.0 - cc),
        fortran.physics_variables.beta_poloidal * cc,
    )


@ConstraintManager.register_constraint(50, "Hz", "<=")
def constraint_equation_50():
    """IFE option: Equation for repetition rate upper limit
    author: P B Lloyd, CCFE, Culham Science Centre
    author: S I Muldrew, CCFE, Culham Science Centre
    """
    cc = (
        fortran.ife_variables.reprat / fortran.ife_variables.rrmax
        - 1.0 * fortran.ife_variables.frrmax
    )
    return ConstraintResult(
        cc, fortran.ife_variables.rrmax * (1.0 - cc), fortran.ife_variables.reprat * cc
    )


@ConstraintManager.register_constraint(51, "V.s", "=")
def constraint_equation_51():
    """Equation to enforce startup flux = available startup flux
    author: P B Lloyd, CCFE, Culham Science Centre

    vs_plasma_res_ramp: resistive losses in startup V-s (Wb)
    vs_plasma_ind_ramp: internal and external plasma inductance V-s (Wb))
    vs_cs_pf_total_ramp: total flux swing for startup (Wb)
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


@ConstraintManager.register_constraint(52, "", ">=")
def constraint_equation_52():
    """Equation for tritium breeding ratio lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ftbr: f-value for minimum tritium breeding ratio
    tbr: tritium breeding ratio (i_blanket_type=2,3 (KIT HCPB/HCLL))
    tbrmin: minimum tritium breeding ratio (If i_blanket_type=1, tbrmin=minimum 5-year time-averaged tritium breeding ratio)
    """
    cc = (
        1.0
        - fortran.constraint_variables.ftbr
        * fortran.fwbs_variables.tbr
        / fortran.constraint_variables.tbrmin
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.tbrmin * (1.0 - cc),
        fortran.constraint_variables.tbrmin * cc,
    )


@ConstraintManager.register_constraint(53, "neutron/m2", "<=")
def constraint_equation_53():
    """Equation for fast neutron fluence on TF coil upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fflutf: f-value for maximum TF coil nuclear heating
    nflutfmax: max fast neutron fluence on TF coil (n/m2)
    nflutf: peak fast neutron fluence on TF coil superconductor (n/m2)
    """
    cc = (
        fortran.fwbs_variables.nflutf / fortran.constraint_variables.nflutfmax
        - 1.0 * fortran.constraint_variables.fflutf
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.nflutfmax * (1.0 - cc),
        fortran.fwbs_variables.nflutf * cc,
    )


@ConstraintManager.register_constraint(54, "MW/m3", "<=")
def constraint_equation_54():
    """Equation for peak TF coil nuclear heating upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fptfnuc: f-value for maximum TF coil nuclear heating
    ptfnucmax: maximum nuclear heating in TF coil (MW/m3)
    ptfnucpm3: nuclear heating in the TF coil (MW/m3) (blktmodel>0)
    """
    cc = (
        fortran.fwbs_variables.ptfnucpm3 / fortran.constraint_variables.ptfnucmax
        - 1.0 * fortran.constraint_variables.fptfnuc
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.ptfnucmax * (1.0 - cc),
        fortran.fwbs_variables.ptfnucpm3 * cc,
    )


@ConstraintManager.register_constraint(56, "MW/m", "<=")
def constraint_equation_56():
    """Equation for power through separatrix / major radius upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fpsepr: f-value for maximum Psep/R limit
    pseprmax: maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    rmajor: plasma major radius (m)
    """
    cc = (
        (
            fortran.physics_variables.p_plasma_separatrix_mw
            / fortran.physics_variables.rmajor
        )
        / fortran.constraint_variables.pseprmax
        - 1.0 * fortran.constraint_variables.fpsepr
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.pseprmax * (1.0 - cc),
        (
            fortran.physics_variables.p_plasma_separatrix_mw
            / fortran.physics_variables.rmajor
        )
        * cc,
    )


@ConstraintManager.register_constraint(59, "", "<=")
def constraint_equation_59():
    """Equation for neutral beam shine-through fraction upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fnbshinef: f-value for maximum neutral beam shine-through fraction
    f_p_beam_shine_through_max: maximum neutral beam shine-through fraction
    f_p_beam_shine_through: neutral beam shine-through fraction
    """
    cc = (
        fortran.current_drive_variables.f_p_beam_shine_through
        / fortran.constraint_variables.f_p_beam_shine_through_max
        - 1.0 * fortran.constraint_variables.fnbshinef
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.f_p_beam_shine_through_max * (1.0 - cc),
        fortran.current_drive_variables.f_p_beam_shine_through * cc,
    )


@ConstraintManager.register_constraint(60, "K", ">=")
def constraint_equation_60():
    """Equation for Central Solenoid s/c temperature margin lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    ftmargoh: f-value for central solenoid temperature margin
    temp_cs_margin: Central solenoid temperature margin (K)
    tmargmin_cs: Minimum allowable temperature margin : CS (K)
    """
    return ConstraintResult(
        1.0
        - fortran.constraint_variables.ftmargoh
        * fortran.pfcoil_variables.temp_cs_margin
        / fortran.tfcoil_variables.tmargmin_cs,
        fortran.tfcoil_variables.tmargmin_cs,
        fortran.tfcoil_variables.tmargmin_cs - fortran.pfcoil_variables.temp_cs_margin,
    )


@ConstraintManager.register_constraint(61, "", ">=")
def constraint_equation_61():
    """Equation for availability lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    favail: F-value for minimum availability
    cfactr: Total plant availability fraction
    avail_min: Minimum availability
    """
    cc = (
        1.0
        - fortran.cost_variables.favail
        * fortran.cost_variables.cfactr
        / fortran.cost_variables.avail_min
    )
    return ConstraintResult(
        cc,
        fortran.cost_variables.avail_min * (1.0 - cc),
        fortran.cost_variables.cfactr * cc,
    )


@ConstraintManager.register_constraint(62, "", ">=")
def constraint_equation_62():
    """Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
    author: P B Lloyd, CCFE, Culham Science Centre

    falpha_energy_confinement: f-value for lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement
    t_alpha_confinement: alpha particle confinement time (s)
    t_energy_confinement: global thermal energy confinement time (sec)
    f_alpha_energy_confinement_min: Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement times
    """
    cc = (
        1.0
        - fortran.constraint_variables.falpha_energy_confinement
        * (
            fortran.physics_variables.t_alpha_confinement
            / fortran.physics_variables.t_energy_confinement
        )
        / fortran.constraint_variables.f_alpha_energy_confinement_min
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.f_alpha_energy_confinement_min,
        (
            fortran.physics_variables.t_alpha_confinement
            / fortran.physics_variables.t_energy_confinement
        )
        * cc,
    )


@ConstraintManager.register_constraint(63, "", "<=")
def constraint_equation_63():
    """Upper limit on niterpump (vacuum_model = simple)
    author: P B Lloyd, CCFE, Culham Science Centre

    fniterpump: f-value for constraint that number of pumps < tfno
    tfno: number of TF coils (default = 50 for stellarators)
    niterpump: number of high vacuum pumps (real number), each with the throughput
    """
    cc = (
        fortran.vacuum_variables.niterpump / fortran.tfcoil_variables.n_tf_coils
        - 1.0 * fortran.constraint_variables.fniterpump
    )
    return ConstraintResult(
        cc,
        fortran.tfcoil_variables.n_tf_coils,
        fortran.tfcoil_variables.n_tf_coils * cc,
    )


@ConstraintManager.register_constraint(64, "", "<=")
def constraint_equation_64():
    """Upper limit on Zeff
    author: P B Lloyd, CCFE, Culham Science Centre

    fzeffmax: f-value for maximum zeff
    zeffmax: maximum value for Zeff
    zeff: plasma effective charge
    """
    cc = (
        fortran.physics_variables.zeff / fortran.constraint_variables.fzeffmax
        - 1.0 * fortran.constraint_variables.fzeffmax
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.fzeffmax,
        fortran.constraint_variables.fzeffmax * cc,
    )


@ConstraintManager.register_constraint(65, "Pa", "<=")
def constraint_equation_65():
    """Upper limit on stress of the vacuum vessel that occurs when the TF coil quenches.
    author: Timothy Nunn, UKAEA

    fmaxvvstress: f-value for constraint on maximum VV stress
    max_vv_stress: Maximum permitted stress of the VV (Pa)
    vv_stress_quench: Stress of the VV (Pa)
    """
    cc = (
        fortran.sctfcoil_module.vv_stress_quench
        / fortran.tfcoil_variables.max_vv_stress
        - 1.0 * fortran.constraint_variables.fmaxvvstress
    )
    return ConstraintResult(
        cc,
        fortran.tfcoil_variables.max_vv_stress,
        fortran.tfcoil_variables.max_vv_stress * cc,
    )


@ConstraintManager.register_constraint(66, "MW", "<=")
def constrain_equation_66():
    """Upper limit on rate of change of energy in poloidal field
    author: P B Lloyd, CCFE, Culham Science Centre

    fpoloidalpower: f-value for constraint on rate of change of energy in poloidal field
    maxpoloidalpower: Maximum permitted absolute rate of change of stored energy in poloidal field (MW)
    peakpoloidalpower: Peak absolute rate of change of stored energy in poloidal field (MW) (11/01/16)
    """
    cc = (
        fortran.pf_power_variables.peakpoloidalpower
        / fortran.pf_power_variables.maxpoloidalpower
        - 1.0 * fortran.constraint_variables.fpoloidalpower
    )
    return ConstraintResult(
        cc,
        fortran.pf_power_variables.maxpoloidalpower,
        fortran.pf_power_variables.maxpoloidalpower * cc,
    )


@ConstraintManager.register_constraint(67, "MW/m2", "<=")
def constraint_equation_67():
    """Simple upper limit on radiation wall load
    author: P B Lloyd, CCFE, Culham Science Centre

    fradwall: f-value for upper limit on radiation wall load
    pflux_fw_rad_max: Maximum permitted radiation wall load (MW/m^2)
    pflux_fw_rad_max_mw: Peak radiation wall load (MW/m^2)
    """
    cc = (
        fortran.constraint_variables.pflux_fw_rad_max_mw
        / fortran.constraint_variables.pflux_fw_rad_max
        - 1.0 * fortran.constraint_variables.fradwall
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.pflux_fw_rad_max,
        fortran.constraint_variables.pflux_fw_rad_max * cc,
    )


@ConstraintManager.register_constraint(68, "MWT/m", "<=")
def constraint_equation_68():
    """Upper limit on Psep scaling (PsepB/qAR)
    author: P B Lloyd, CCFE, Culham Science Centre

    fpsepbqar: f-value for upper limit on psepbqar, maximum Psep*Bt/qAR limit
    psepbqarmax: maximum permitted value of ratio of Psep*Bt/qAR (MWT/m)
    p_plasma_separatrix_mw: Power to conducted to the divertor region (MW)
    bt: toroidal field on axis (T) (iteration variable 2)
    q95: safety factor q at 95% flux surface
    aspect: aspect ratio (iteration variable 1)
    rmajor: plasma major radius (m) (iteration variable 3)
    i_q95_fixed: Switch that allows for fixing q95 only in this constraint.
    q95_fixed: fixed safety factor q at 95% flux surface
    """
    if fortran.constraint_variables.i_q95_fixed == 1:
        cc = (
            (
                (
                    fortran.physics_variables.p_plasma_separatrix_mw
                    * fortran.physics_variables.bt
                )
                / (
                    fortran.constraint_variables.q95_fixed
                    * fortran.physics_variables.aspect
                    * fortran.physics_variables.rmajor
                )
            )
            / fortran.constraint_variables.psepbqarmax
            - 1.0 * fortran.constraint_variables.fpsepbqar
        )
        err = (
            fortran.physics_variables.p_plasma_separatrix_mw
            * fortran.physics_variables.bt
        ) / (
            fortran.constraint_variables.q95_fixed
            * fortran.physics_variables.aspect
            * fortran.physics_variables.rmajor
        ) - fortran.constraint_variables.psepbqarmax
    else:
        cc = (
            (
                (
                    fortran.physics_variables.p_plasma_separatrix_mw
                    * fortran.physics_variables.bt
                )
                / (
                    fortran.physics_variables.q95
                    * fortran.physics_variables.aspect
                    * fortran.physics_variables.rmajor
                )
            )
            / fortran.constraint_variables.psepbqarmax
            - 1.0 * fortran.constraint_variables.fpsepbqar
        )
        err = (
            fortran.physics_variables.p_plasma_separatrix_mw
            * fortran.physics_variables.bt
        ) / (
            fortran.physics_variables.q95
            * fortran.physics_variables.aspect
            * fortran.physics_variables.rmajor
        ) - fortran.constraint_variables.psepbqarmax

    return ConstraintResult(cc, fortran.constraint_variables.psepbqarmax, err)


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

    foh_stress: f-value for Tresca yield criterion in Central Solenoid
    alstroh: allowable hoop stress in Central Solenoid structural material (Pa)
    s_shear_cs_peak: Maximum shear stress coils/central solenoid (Pa)
    sig_tf_cs_bucked: Maximum shear stress in CS case at flux swing (no current in CS)
                          can be significant for the bucked and weged design
    i_tf_bucking: switch for TF structure design
    """
    # bucked and wedged desing
    if (
        fortran.tfcoil_variables.i_tf_bucking >= 2
        and fortran.build_variables.i_tf_inside_cs == 0
    ):
        cc = (
            max(
                fortran.pfcoil_variables.s_shear_cs_peak,
                fortran.tfcoil_variables.sig_tf_cs_bucked,
            )
            / fortran.pfcoil_variables.alstroh
            - 1.0 * fortran.constraint_variables.foh_stress
        )
        err = fortran.pfcoil_variables.alstroh - max(
            fortran.pfcoil_variables.s_shear_cs_peak,
            fortran.tfcoil_variables.sig_tf_cs_bucked,
        )
    # Free standing CS
    else:
        cc = (
            fortran.pfcoil_variables.s_shear_cs_peak / fortran.pfcoil_variables.alstroh
            - 1.0 * fortran.constraint_variables.foh_stress
        )
        err = (
            fortran.pfcoil_variables.alstroh - fortran.pfcoil_variables.s_shear_cs_peak
        )

    return ConstraintResult(cc, fortran.pfcoil_variables.alstroh, err)


@ConstraintManager.register_constraint(73, "MW", ">=")
def constraint_equation_73():
    """Lower limit to ensure separatrix power is greater than the L-H power + auxiliary power
    Related to constraint 15
    author: P B Lloyd, CCFE, Culham Science Centre

    fplhsep: F-value for Psep >= Plh + Paux : for consistency of two values of separatrix power
    p_l_h_threshold_mw: L-H mode power threshold (MW)
    p_plasma_separatrix_mw: power to be conducted to the divertor region (MW)
    p_hcd_injected_total_mw : inout real : total auxiliary injected power (MW)
    """
    cc = (
        1.0
        - fortran.physics_variables.fplhsep
        * fortran.physics_variables.p_plasma_separatrix_mw
        / (
            fortran.physics_variables.p_l_h_threshold_mw
            + fortran.current_drive_variables.p_hcd_injected_total_mw
        )
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.p_plasma_separatrix_mw,
        fortran.physics_variables.p_plasma_separatrix_mw * cc,
    )


@ConstraintManager.register_constraint(74, "K", "<=")
def constraint_equation_74():
    """Upper limit to ensure TF coil quench temperature < tmax_croco
    ONLY used for croco HTS coil
    author: P B Lloyd, CCFE, Culham Science Centre

    fcqt: f-value: TF coil quench temparature remains below tmax_croco
    croco_quench_temperature: CroCo strand: Actual temp reached during a quench (K)
    tmax_croco: CroCo strand: maximum permitted temp during a quench (K)
    """
    cc = (
        fortran.tfcoil_variables.croco_quench_temperature
        / fortran.tfcoil_variables.tmax_croco
        - 1.0 * fortran.constraint_variables.fcqt
    )
    return ConstraintResult(
        cc,
        fortran.tfcoil_variables.croco_quench_temperature,
        fortran.tfcoil_variables.croco_quench_temperature * cc,
    )


@ConstraintManager.register_constraint(75, "A/m2", "<=")
def constraint_equation_75():
    """Upper limit to ensure that TF coil current / copper area < Maximum value
    ONLY used for croco HTS coil
    author: P B Lloyd, CCFE, Culham Science Centre

    copperA_m2: TF coil current / copper area (A/m2)
    copperA_m2_max: Maximum TF coil current / copper area (A/m2)
    f_coppera_m2: f-value for TF coil current / copper area < copperA_m2_max
    """
    cc = (
        fortran.rebco_variables.coppera_m2 / fortran.rebco_variables.coppera_m2_max
        - 1.0 * fortran.rebco_variables.f_coppera_m2
    )
    return ConstraintResult(
        cc, fortran.rebco_variables.coppera_m2, fortran.rebco_variables.coppera_m2 * cc
    )


@ConstraintManager.register_constraint(76, "m-3", "<=")
def constraint_equation_76():
    """Upper limit for Eich critical separatrix density model: Added for issue 558
    author: P B Lloyd, CCFE, Culham Science Centre

    Eich critical separatrix density model
    Added for issue 558 with ref to http://iopscience.iop.org/article/10.1088/1741-4326/aaa340/pdf

    alpha_crit: critical ballooning parameter value
    nesep_crit: critical electron density at separatrix [m-3]
    kappa: plasma separatrix elongation (calculated if i_plasma_geometry = 1-5, 7 or 9)
    triang: plasma separatrix triangularity (calculated if i_plasma_geometry = 1, 3-5 or 7)
    aspect: aspect ratio (iteration variable 1)
    p_plasma_separatrix_mw: power to conducted to the divertor region (MW)
    dlimit(7)array : density limit (/m3) as calculated using various models
    fnesep: f-value for Eich critical separatrix density
    """
    # TODO: why on earth are these variables being set here!? Should they be local?
    fortran.physics_variables.alpha_crit = (fortran.physics_variables.kappa**1.2) * (
        1.0 + 1.5 * fortran.physics_variables.triang
    )
    fortran.physics_variables.nesep_crit = (
        5.9
        * fortran.physics_variables.alpha_crit
        * (fortran.physics_variables.aspect ** (-2.0 / 7.0))
        * (((1.0 + (fortran.physics_variables.kappa**2.0)) / 2.0) ** (-6.0 / 7.0))
        * ((fortran.physics_variables.p_plasma_separatrix_mw * 1.0e6) ** (-11.0 / 70.0))
        * fortran.physics_variables.dlimit[6]
    )

    cc = (
        fortran.physics_variables.nesep / fortran.physics_variables.nesep_crit
        - 1.0 * fortran.constraint_variables.fnesep
    )
    return ConstraintResult(
        cc, fortran.physics_variables.nesep, fortran.physics_variables.nesep * cc
    )


@ConstraintManager.register_constraint(77, "A/turn", "<=")
def constraint_equation_77():
    """Equation for maximum TF current per turn upper limit
    author: P B Lloyd, CCFE, Culham Science Centre

    fcpttf: f-value for TF coil current per turn
    cpttf_max : allowable TF coil current per turn [A/turn]
    c_tf_turn : TF coil current per turn [A/turn]
    """
    cc = (
        fortran.tfcoil_variables.c_tf_turn / fortran.tfcoil_variables.cpttf_max
        - 1.0 * fortran.constraint_variables.fcpttf
    )
    return ConstraintResult(
        cc, fortran.tfcoil_variables.cpttf_max, fortran.tfcoil_variables.cpttf_max * cc
    )


@ConstraintManager.register_constraint(78, "", ">=")
def constraint_equation_78():
    """Equation for Reinke criterion, divertor impurity fraction lower limit
    author: P B Lloyd, CCFE, Culham Science Centre

    freinke : input : f-value for Reinke criterion (itv 147)
    fzmin : input : minimum impurity fraction from Reinke model
    fzactual : input : actual impurity fraction
    """
    cc = (
        1.0
        - fortran.constraint_variables.freinke
        * fortran.reinke_variables.fzactual
        / fortran.reinke_variables.fzmin
    )
    return ConstraintResult(
        cc,
        fortran.reinke_variables.fzmin * (1.0 - cc),
        fortran.reinke_variables.fzmin * cc,
    )


@ConstraintManager.register_constraint(79, "A/turn", "<=")
def constraint_equation_79():
    """Equation for maximum CS field
    author: P B Lloyd, CCFE, Culham Science Centre

    fb_cs_limit_max: F-value for CS mmax field (cons. 79, itvar 149)
    b_cs_limit_max: Central solenoid max field limit [T]
    b_cs_peak_pulse_start: maximum field in central solenoid at beginning of pulse (T)
    b_cs_peak_flat_top_end: maximum field in central solenoid at end of flat-top (EoF) (T)
    (Note: original code has "b_cs_peak_flat_top_end/b_cs_peak_pulse_start |  peak CS field [T]".)
    """
    cc = (
        max(
            fortran.pfcoil_variables.b_cs_peak_flat_top_end,
            fortran.pfcoil_variables.b_cs_peak_pulse_start,
        )
        / fortran.pfcoil_variables.b_cs_limit_max
        - 1.0 * fortran.pfcoil_variables.fb_cs_limit_max
    )
    return ConstraintResult(
        cc,
        fortran.pfcoil_variables.b_cs_limit_max,
        max(
            fortran.pfcoil_variables.b_cs_peak_flat_top_end,
            fortran.pfcoil_variables.b_cs_peak_pulse_start,
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
    #=# physics
    #=#=# fp_plasma_separatrix_min_mw, p_plasma_separatrix_mw
    Logic change during pre-factoring: err, symbol, units will be assigned only if present.
    fp_plasma_separatrix_min_mw : input : F-value for lower limit on p_plasma_separatrix_mw (cons. 80, itvar 153)
    p_plasma_separatrix_min_mw : input : Minimum power crossing separatrix p_plasma_separatrix_mw [MW]
    p_plasma_separatrix_mw : input : Power crossing separatrix [MW]
    """
    cc = (
        1.0
        - fortran.physics_variables.fp_plasma_separatrix_min_mw
        * fortran.physics_variables.p_plasma_separatrix_mw
        / fortran.constraint_variables.p_plasma_separatrix_min_mw
    )
    return ConstraintResult(
        cc,
        fortran.constraint_variables.p_plasma_separatrix_min_mw,
        fortran.constraint_variables.p_plasma_separatrix_min_mw * cc,
    )


@ConstraintManager.register_constraint(81, "m-3", ">=")
def constraint_equation_81():
    """Lower limit to ensure central density is larger that the pedestal one
    author: S Kahn, Culham Science Centre
    args : output structure : residual error; constraint value;
    residual error in physical units; output string; units string
    Lower limit ne0 > neped
    !#=# physics
    !#=#=# ne0, neped
    Logic change during pre-factoring: err, symbol, units will be
    assigned only if present.
    fne0  : input : F-value for constraint on ne0 > neped
    ne0   : input : Central electron density [m-3]
    neped : input : Electron density at pedestal [m-3]
    """
    cc = (
        1.0
        - fortran.physics_variables.fne0
        * fortran.physics_variables.ne0
        / fortran.physics_variables.neped
    )
    return ConstraintResult(
        cc, fortran.physics_variables.fne0, fortran.physics_variables.fne0 * cc
    )


@ConstraintManager.register_constraint(82, "m", ">=")
def constraint_equation_82():
    """Equation for toroidal consistency of stellarator build
    author: J Lion, IPP Greifswald

    ftoroidalgap: f-value for constraint toroidalgap > dx_tf_inboard_out_toroidal
    toroidalgap: minimal gap between two stellarator coils
    dx_tf_inboard_out_toroidal: total toroidal width of a tf coil
    """
    return ConstraintResult(
        1.0
        - fortran.tfcoil_variables.ftoroidalgap
        * fortran.tfcoil_variables.toroidalgap
        / fortran.tfcoil_variables.dx_tf_inboard_out_toroidal,
        fortran.tfcoil_variables.toroidalgap,
        fortran.tfcoil_variables.toroidalgap
        - fortran.tfcoil_variables.dx_tf_inboard_out_toroidal
        / fortran.tfcoil_variables.ftoroidalgap,
    )


@ConstraintManager.register_constraint(83, "m", ">=")
def constraint_equation_83():
    """Equation for radial consistency of stellarator build
    author: J Lion, IPP Greifswald

    f_avspace: f-value for constraint available_radial_space > required_radial_space
    available_radial_space: avaible space in radial direction as given by each s.-configuration
    required_radial_space: required space in radial direction
    """
    cc = (
        1.0
        - fortran.build_variables.f_avspace
        * fortran.build_variables.available_radial_space
        / fortran.build_variables.required_radial_space
    )
    return ConstraintResult(
        cc,
        fortran.build_variables.available_radial_space * (1.0 - cc),
        fortran.build_variables.required_radial_space * cc,
    )


@ConstraintManager.register_constraint(84, "", ">=")
def constraint_equation_84():
    """Equation for the lower limit of beta
    author: J Lion, IPP Greifswald

    fbeta_min: f-value for constraint beta-beta_fast_alpha > beta_min
    beta_min: Lower limit for beta
    beta: plasma beta
    """
    cc = (
        1.0
        - fortran.constraint_variables.fbeta_min
        * fortran.physics_variables.beta
        / fortran.physics_variables.beta_min
    )
    return ConstraintResult(
        cc,
        fortran.physics_variables.beta_min * (1.0 - cc),
        fortran.physics_variables.beta * cc,
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


@ConstraintManager.register_constraint(86, "m", "<=")
def constraint_equation_86():
    """Upper limit on the turn edge length in the TF winding pack
    Author : S Kahn

    t_turn_tf: TF coil turn edge length including turn insulation [m]
    f_t_turn_tf: f-value for TF turn edge length constraint
    t_turn_tf_max: TF turn edge length including turn insulation upper limit [m]
    """
    cc = (
        fortran.tfcoil_variables.t_turn_tf / fortran.tfcoil_variables.t_turn_tf_max
        - 1.0 * fortran.tfcoil_variables.f_t_turn_tf
    )
    return ConstraintResult(
        cc,
        fortran.tfcoil_variables.t_turn_tf_max * (1.0 - cc),
        fortran.tfcoil_variables.t_turn_tf_max * cc,
    )


@ConstraintManager.register_constraint(87, "MW", "<=")
def constraint_equation_87():
    """Equation for TF coil cryogenic power upper limit
    author: S. Kahn, CCFE, Culham Science Centre

    p_cryo_plant_electric_mw: cryogenic plant power (MW)
    f_crypmw: f-value for maximum cryogenic plant power
    p_cryo_plant_electric_max_mw: Maximum cryogenic plant power (MW)
    """
    cc = (
        fortran.heat_transport_variables.p_cryo_plant_electric_mw
        / fortran.heat_transport_variables.p_cryo_plant_electric_max_mw
        - 1.0 * fortran.heat_transport_variables.f_crypmw
    )
    return ConstraintResult(
        cc,
        fortran.heat_transport_variables.p_cryo_plant_electric_max_mw * (1.0 - cc),
        fortran.heat_transport_variables.p_cryo_plant_electric_mw * cc,
    )


@ConstraintManager.register_constraint(88, "", "<=")
def constraint_equation_88():
    """Equation for TF coil vertical strain upper limit (absolute value)
    author: CPS Swanson, PPPL, USA

    fstr_wp: f-value for TF coil strain
    str_wp_max: Allowable maximum TF coil vertical strain
    str_wp: Constrained TF coil vertical strain
    """
    return ConstraintResult(
        abs(fortran.tfcoil_variables.str_wp) / fortran.tfcoil_variables.str_wp_max
        - 1.0 * fortran.constraint_variables.fstr_wp,
        fortran.tfcoil_variables.str_wp_max,
        fortran.tfcoil_variables.str_wp_max
        - abs(fortran.tfcoil_variables.str_wp) / fortran.constraint_variables.fstr_wp,
    )


@ConstraintManager.register_constraint(89, "A/m2", "<=")
def constraint_equation_89():
    """Upper limit to ensure that the Central Solenoid [OH] coil current / copper area < Maximum value
    author: G Turkington, CCFE, Culham Science Centre

    copperaoh_m2: CS coil current at EOF / copper area [A/m2]
    copperaoh_m2_max: maximum coil current / copper area [A/m2]
    f_copperaoh_m2: f-value for CS coil current / copper area
    """
    cc = (
        fortran.rebco_variables.copperaoh_m2 / fortran.rebco_variables.copperaoh_m2_max
        - 1.0 * fortran.rebco_variables.f_copperaoh_m2
    )
    return ConstraintResult(
        cc,
        fortran.rebco_variables.copperaoh_m2,
        fortran.rebco_variables.copperaoh_m2 * cc,
    )


@ConstraintManager.register_constraint(90, "", ">=")
def constraint_equation_90():
    """Lower limit for CS coil stress load cycles
    author: A. Pearce, G Turkington CCFE, Culham Science Centre

    fncycle: f-value for constraint n_cycle > n_cycle_min
    n_cycle: Allowable number of cycles for CS
    n_cycle_min: Minimum required cycles for CS
    """
    if (
        fortran.cost_variables.ibkt_life == 1
        and fortran.cs_fatigue_variables.bkt_life_csf == 1
    ):
        fortran.cs_fatigue_variables.n_cycle_min = fortran.cost_variables.bktcycles

    cc = (
        1.0
        - fortran.constraint_variables.fncycle
        * fortran.cs_fatigue_variables.n_cycle
        / fortran.cs_fatigue_variables.n_cycle_min
    )
    return ConstraintResult(
        cc,
        fortran.cs_fatigue_variables.n_cycle_min * (1.0 - cc),
        fortran.cs_fatigue_variables.n_cycle * cc,
    )


@ConstraintManager.register_constraint(91, "MW", ">=")
def constraint_equation_91():
    """Lower limit to ensure ECRH te is greater than required te for ignition
    at lower values for n and B. Or if the design point is ECRH heatable (if i_plasma_ignited==0)
    stellarators only (but in principle usable also for tokamaks).
    author: J Lion, IPP Greifswald

    fecrh_ignition: f-value for constraint powerht_local > powerscaling
    max_gyrotron_frequency: Max. av. gyrotron frequency
    te0_ecrh_achievable: Max. achievable electron temperature at ignition point
    """
    # Achievable ECRH te needs to be larger than needed te for igntion
    if fortran.physics_variables.i_plasma_ignited == 0:
        cc = (
            1.0
            - fortran.constraint_variables.fecrh_ignition
            * (
                fortran.stellarator_variables.powerht_constraint
                + fortran.current_drive_variables.p_hcd_primary_extra_heat_mw
            )
            / fortran.stellarator_variables.powerscaling_constraint
        )
    else:
        cc = (
            1.0
            - fortran.constraint_variables.fecrh_ignition
            * fortran.stellarator_variables.powerht_constraint
            / fortran.stellarator_variables.powerscaling_constraint
        )

    return ConstraintResult(
        cc,
        fortran.stellarator_variables.powerscaling_constraint * (1.0 - cc),
        fortran.stellarator_variables.powerht_constraint * cc,
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


def init_constraint_variables():
    """Initialise the constraint variables"""
    fortran.constraint_variables.auxmin = 0.1
    fortran.constraint_variables.beta_poloidal_max = 0.19
    fortran.constraint_variables.bigqmin = 10.0
    fortran.constraint_variables.b_tf_inboard_max = 12.0
    fortran.constraint_variables.fauxmn = 1.0
    fortran.constraint_variables.fbeta_poloidal_eps = 1.0
    fortran.constraint_variables.fbeta_poloidal = 1.0
    fortran.constraint_variables.fbeta_max = 1.0
    fortran.constraint_variables.fbeta_min = 1.0
    fortran.constraint_variables.fcpttf = 1.0
    fortran.constraint_variables.fr_conducting_wall = 1.0
    fortran.constraint_variables.fdene = 1.0
    fortran.constraint_variables.fdtmp = 1.0
    fortran.constraint_variables.fflutf = 1.0
    fortran.constraint_variables.fp_fusion_total_max_mw = 1.0
    fortran.constraint_variables.fgamcd = 1.0
    fortran.constraint_variables.fpflux_div_heat_load_mw = 1.0
    fortran.constraint_variables.fiooic = 0.5
    fortran.constraint_variables.fipir = 1.0
    fortran.constraint_variables.q95_fixed = 3.0
    fortran.constraint_variables.fjohc = 1.0
    fortran.constraint_variables.fjohc0 = 1.0
    fortran.constraint_variables.fjprot = 1.0
    fortran.constraint_variables.fl_h_threshold = 1.0
    fortran.constraint_variables.fmva = 1.0
    fortran.constraint_variables.fnbshinef = 1.0
    fortran.constraint_variables.fncycle = 1.0
    fortran.constraint_variables.fnesep = 1.0
    fortran.constraint_variables.foh_stress = 1.0
    fortran.constraint_variables.fb_tf_inboard_max = 1.0
    fortran.constraint_variables.fp_hcd_injected_max = 1.0
    fortran.constraint_variables.fp_plant_electric_net_required_mw = 1.0
    fortran.constraint_variables.fradius_beam_tangency = 1.0
    fortran.constraint_variables.fpsepbqar = 1.0
    fortran.constraint_variables.fpsepr = 1.0
    fortran.constraint_variables.fptemp = 1.0
    fortran.constraint_variables.fptfnuc = 1.0
    fortran.constraint_variables.fq = 1.0
    fortran.constraint_variables.fqval = 1.0
    fortran.constraint_variables.fradpwr = 0.99
    fortran.constraint_variables.fradwall = 1.0
    fortran.constraint_variables.freinke = 1.0
    fortran.constraint_variables.frminor = 1.0
    fortran.constraint_variables.fstrcase = 1.0
    fortran.constraint_variables.fstrcond = 1.0
    fortran.constraint_variables.fstr_wp = 1.0
    fortran.constraint_variables.fmaxvvstress = 1.0
    fortran.constraint_variables.ftbr = 1.0
    fortran.constraint_variables.ft_burn_min = 1.0
    fortran.constraint_variables.ftcycl = 1.0
    fortran.constraint_variables.ftmargoh = 1.0
    fortran.constraint_variables.ftmargtf = 1.0
    fortran.constraint_variables.ft_current_ramp_up = 1.0
    fortran.constraint_variables.ftpeak = 1.0
    fortran.constraint_variables.fvdump = 1.0
    fortran.constraint_variables.fvs_plasma_total_required = 1.0
    fortran.constraint_variables.fvvhe = 1.0
    fortran.constraint_variables.fpflux_fw_neutron_max_mw = 1.0
    fortran.constraint_variables.fzeffmax = 1.0
    fortran.constraint_variables.eta_cd_norm_hcd_primary_max = 2.0
    fortran.constraint_variables.i_q95_fixed = 0
    fortran.constraint_variables.pflux_fw_rad_max = 1.0
    fortran.constraint_variables.mvalim = 40.0
    fortran.constraint_variables.f_p_beam_shine_through_max = 1e-3
    fortran.constraint_variables.nflutfmax = 1.0e23
    fortran.constraint_variables.p_plasma_separatrix_min_mw = 150.0
    fortran.constraint_variables.f_fw_rad_max = 3.33
    fortran.constraint_variables.pflux_fw_rad_max_mw = 0.0
    fortran.constraint_variables.p_plant_electric_net_required_mw = 1.0e3
    fortran.constraint_variables.p_fusion_total_max_mw = 1.5e3
    fortran.constraint_variables.psepbqarmax = 9.5
    fortran.constraint_variables.pseprmax = 25.0
    fortran.constraint_variables.ptfnucmax = 1e-3
    fortran.constraint_variables.tbrmin = 1.1
    fortran.constraint_variables.t_burn_min = 1.0
    fortran.constraint_variables.tcycmn = 0.0
    fortran.constraint_variables.t_current_ramp_up_min = 1.0
    fortran.constraint_variables.vvhealw = 1.0
    fortran.constraint_variables.pflux_fw_neutron_max_mw = 1.0
    fortran.constraint_variables.f_alpha_energy_confinement_min = 5.0
    fortran.constraint_variables.falpha_energy_confinement = 1.0
    fortran.constraint_variables.fniterpump = 1.0
    fortran.constraint_variables.zeffmax = 3.6
    fortran.constraint_variables.fpoloidalpower = 1.0
    fortran.constraint_variables.fpsep = 1.0
    fortran.constraint_variables.fcqt = 1.0
    fortran.constraint_variables.fecrh_ignition = 1.0
