from __future__ import annotations

import logging
import math
from enum import IntEnum
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numba as nb
import numpy as np

import process.models.physics.fusion_reactions as reactions
import process.models.physics.radiation_power as physics_funcs
from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure import (
    constraint_variables,
    current_drive_variables,
    divertor_variables,
    impurity_radiation_module,
    numerics,
    physics_variables,
    pulse_variables,
    reinke_variables,
    stellarator_variables,
    times_variables,
)
from process.models.physics import impurity_radiation
from process.models.physics.profiles import PlasmaProfileShapeType

if TYPE_CHECKING:
    from process.models.physics.bootstrap_current import PlasmaBootstrapCurrent
    from process.models.physics.confinement_time import (
        PlasmaConfinementTime,
    )
    from process.models.physics.density_limit import PlasmaDensityLimit
    from process.models.physics.exhaust import PlasmaExhaust
    from process.models.physics.l_h_transition import PlasmaConfinementTransition
    from process.models.physics.plasma_current import (
        PlasmaCurrent,
        PlasmaDiamagneticCurrent,
    )
    from process.models.physics.plasma_fields import PlasmaFields
    from process.models.physics.plasma_geometry import PlasmaGeom

logger = logging.getLogger(__name__)


@nb.jit(nopython=True, cache=True)
def calculate_cylindrical_safety_factor(
    rmajor: float,
    rminor: float,
    plasma_current: float,
    b_plasma_toroidal_on_axis: float,
    kappa95: float,
    triang95: float,
) -> float:
    """Calculate the cylindrical safety factor from the IPDG89 guidelines.

    Parameters
    ----------
    rmajor : float
        Major radius of the tokamak in meters.
    rminor : float
        Minor radius of the tokamak in meters.
    plasma_current : float
        Plasma current in amperes.
    b_plasma_toroidal_on_axis : float
        Toroidal magnetic field on axis in tesla.
    kappa95 : float
        Elongation at 95% of the plasma boundary.
    triang95 : float
        Triangularity at 95% of the plasma boundary.

    Returns
    -------
    float
        Cylindrical safety factor (dimensionless).

    Notes
    -----
    The cylindrical safety factor is calculated following the IPDG89 guidelines.
    The formula accounts for plasma elongation and triangularity effects on the
    safety factor through the kappa95 and triang95 parameters.

    """
    # Calculate cyclindrical safety factor from IPDG89
    return (
        ((2 * np.pi) / constants.RMU0)
        * rminor**2
        / (rmajor * plasma_current / b_plasma_toroidal_on_axis)
        * 0.5
        * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
    )


@nb.jit(nopython=True, cache=True)
def rether(
    alphan,
    alphat,
    nd_plasma_electrons_vol_avg,
    dlamie,
    te,
    temp_plasma_ion_vol_avg_kev,
    n_charge_plasma_effective_mass_weighted_vol_avg,
):
    """Routine to find the equilibration power between the
    ions and electrons
    This routine calculates the equilibration power between the
    ions and electrons.
    Unknown origin

    Parameters
    ----------
    alphan :
        density profile index
    alphat :
        temperature profile index
    nd_plasma_electrons_vol_avg :
        electron density (/m3)
    dlamie :
        ion-electron coulomb logarithm
    te :
        electron temperature (keV)
    temp_plasma_ion_vol_avg_kev :
        ion temperature (keV)
    n_charge_plasma_effective_mass_weighted_vol_avg :
        mass weighted plasma effective charge

    Returns
    -------
    pden_ion_electron_equilibration_mw  :
        ion/electron equilibration power (MW/m3)

    """
    profie = (1.0 + alphan) ** 2 / (
        (2.0 * alphan - 0.5 * alphat + 1.0) * np.sqrt(1.0 + alphat)
    )
    conie = (
        2.42165e-41
        * dlamie
        * nd_plasma_electrons_vol_avg**2
        * n_charge_plasma_effective_mass_weighted_vol_avg
        * profie
    )

    return conie * (temp_plasma_ion_vol_avg_kev - te) / (te**1.5)


# -----------------------------------------------------
# Diamagnetic and Pfirsch-Schlüter Current Calculations
# -----------------------------------------------------


@nb.jit(nopython=True, cache=True)
def ps_fraction_scene(beta: float) -> float:
    """Calculate the Pfirsch-Schlüter fraction based on the SCENE fit by Tim Hender 2019.

    Parameters
    ----------
    beta :
        the plasma beta value


    Returns
    -------
    :
        the Pfirsch-Schlüter current fraction

    """
    return -9e-2 * beta


class Physics(Model):
    def __init__(
        self,
        plasma_profile,
        current_drive,
        plasma_beta,
        plasma_inductance,
        plasma_density_limit: PlasmaDensityLimit,
        plasma_exhaust: PlasmaExhaust,
        plasma_bootstrap_current: PlasmaBootstrapCurrent,
        plasma_confinement: PlasmaConfinementTime,
        plasma_transition: PlasmaConfinementTransition,
        plasma_current: PlasmaCurrent,
        plasma_fields: PlasmaFields,
        plasma_dia_current: PlasmaDiamagneticCurrent,
        plasma_geometry: PlasmaGeom,
    ):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile
        self.current_drive = current_drive
        self.beta = plasma_beta
        self.inductance = plasma_inductance
        self.density_limit = plasma_density_limit
        self.exhaust = plasma_exhaust
        self.plasma_bootstrap_current = plasma_bootstrap_current
        self.confinement = plasma_confinement
        self.plasma_transition = plasma_transition
        self.current = plasma_current
        self.fields = plasma_fields
        self.dia_current = plasma_dia_current
        self.geometry = plasma_geometry

    def output(self):
        self.calculate_effective_charge_ionisation_profiles()
        self.outplas()
        self.outtim()

    def run(self):
        """Routine to calculate tokamak plasma physics information
        This routine calculates all the primary plasma physics parameters for a tokamak fusion reactor.

        References
        ----------
        - M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants - Part 1: Physics
          https://www.sciencedirect.com/science/article/pii/S0920379614005961
        - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
          https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
        - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
          https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
        """
        # Calculate plasma composition
        # Issue #261 Remove old radiation model (imprad_model=0)
        self.plasma_composition()

        (
            physics_variables.m_plasma_fuel_ions,
            physics_variables.m_plasma_ions_total,
            physics_variables.m_plasma_alpha,
            physics_variables.m_plasma_electron,
            physics_variables.m_plasma,
        ) = self.calculate_plasma_masses(
            physics_variables.m_fuel_amu,
            physics_variables.m_ions_total_amu,
            physics_variables.nd_plasma_ions_total_vol_avg,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.nd_plasma_alphas_vol_avg,
            physics_variables.vol_plasma,
            physics_variables.nd_plasma_electrons_vol_avg,
        )

        # Define coulomb logarithm
        # (collisions: ion-electron, electron-electron)
        physics_variables.dlamee = (
            31.0
            - (np.log(physics_variables.nd_plasma_electrons_vol_avg) / 2.0)
            + np.log(physics_variables.temp_plasma_electron_vol_avg_kev * 1000.0)
        )
        physics_variables.dlamie = (
            31.3
            - (np.log(physics_variables.nd_plasma_electrons_vol_avg) / 2.0)
            + np.log(physics_variables.temp_plasma_electron_vol_avg_kev * 1000.0)
        )

        # Calculate plasma current
        physics_variables.plasma_current = self.current.calculate_plasma_current(
            physics_variables.alphaj,
            physics_variables.alphap,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.eps,
            physics_variables.i_plasma_current,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pres_plasma_thermal_on_axis,
            physics_variables.len_plasma_poloidal,
            physics_variables.q95,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.triang,
            physics_variables.triang95,
        )

        physics_variables.qstar = calculate_cylindrical_safety_factor(
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            plasma_current=physics_variables.plasma_current,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            kappa95=physics_variables.kappa95,
            triang95=physics_variables.triang95,
        )

        # Calculate the poloidal field generated by the plasma current
        physics_variables.b_plasma_surface_poloidal_average = (
            self.fields.calculate_surface_averaged_poloidal_field(
                i_plasma_current=physics_variables.i_plasma_current,
                ip=physics_variables.plasma_current,
                q95=physics_variables.q95,
                aspect=physics_variables.aspect,
                eps=physics_variables.eps,
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                kappa=physics_variables.kappa,
                delta=physics_variables.triang,
                perim=physics_variables.len_plasma_poloidal,
            )
        )

        # -----------------------------------------------------
        # Plasma Current Profile
        # -----------------------------------------------------

        physics_variables.alphaj_wesson = self.calculate_current_profile_index_wesson(
            qstar=physics_variables.qstar, q0=physics_variables.q0
        )

        # Map calculation methods to a dictionary
        alphaj_calculations = {
            0: physics_variables.alphaj,
            1: physics_variables.alphaj_wesson,
        }

        # Calculate alphaj based on i_alphaj
        if int(physics_variables.i_alphaj) in alphaj_calculations:
            physics_variables.alphaj = alphaj_calculations[
                int(physics_variables.i_alphaj)
            ]
        else:
            raise ProcessValueError(
                "Illegal value of i_alphaj",
                i_alphaj=physics_variables.i_alphaj,
            )

        # ==================================================

        # -----------------------------------------------------
        # Plasma Normalised Internal Inductance
        # -----------------------------------------------------

        self.inductance.run()

        # ===================================================

        # Calculate density and temperature profile quantities
        # If physics_variables.i_plasma_pedestal = 1 then set pedestal density to
        #   physics_variables.f_nd_plasma_pedestal_greenwald * Greenwald density limit
        # Note: this used to be done before plasma current
        if (
            PlasmaProfileShapeType(physics_variables.i_plasma_pedestal)
            == PlasmaProfileShapeType.PEDESTAL_PROFILE
        ) and (physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0):
            physics_variables.nd_plasma_pedestal_electron = (
                physics_variables.f_nd_plasma_pedestal_greenwald
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        if (
            PlasmaProfileShapeType(physics_variables.i_plasma_pedestal)
            == PlasmaProfileShapeType.PEDESTAL_PROFILE
        ) and (physics_variables.f_nd_plasma_separatrix_greenwald >= 0e0):
            physics_variables.nd_plasma_separatrix_electron = (
                physics_variables.f_nd_plasma_separatrix_greenwald
                * 1.0e14
                * physics_variables.plasma_current
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        self.plasma_profile.run()

        # Calculate total magnetic field [T]
        physics_variables.b_plasma_total = self.fields.calculate_total_magnetic_field(
            b_plasma_toroidal=physics_variables.b_plasma_toroidal_on_axis,
            b_plasma_poloidal=physics_variables.b_plasma_surface_poloidal_average,
        )

        # Calculate the inboard and outboard toroidal field
        physics_variables.b_plasma_inboard_toroidal = (
            self.fields.calculate_plasma_inboard_toroidal_field(
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
            )
        )

        physics_variables.b_plasma_outboard_toroidal = (
            self.fields.calculate_plasma_outboard_toroidal_field(
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
            )
        )

        # Calculate the toroidal field across the plasma
        # Calculate the toroidal field profile across the plasma (1/R dependence)
        # Double element size to include both sides of the plasma
        physics_variables.b_plasma_toroidal_profile = (
            self.fields.calculate_toroidal_field_profile(
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                rmajor=physics_variables.rmajor,
                rminor=physics_variables.rminor,
                n_plasma_profile_elements=physics_variables.n_plasma_profile_elements,
            )
        )

        # ============================================

        # -----------------------------------------------------
        # Beta Components
        # -----------------------------------------------------

        self.beta.run()

        # =======================================================

        # Set PF coil ramp times
        if pulse_variables.i_pulsed_plant != 1:
            if times_variables.i_t_current_ramp_up == 0:
                times_variables.t_plant_pulse_plasma_current_ramp_up = (
                    physics_variables.plasma_current / 5.0e5
                )
                times_variables.t_plant_pulse_coil_precharge = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )
                times_variables.t_plant_pulse_plasma_current_ramp_down = (
                    times_variables.t_plant_pulse_plasma_current_ramp_up
                )

        elif times_variables.pulsetimings == 0.0e0:  # noqa: RUF069
            # times_variables.t_plant_pulse_coil_precharge is input
            times_variables.t_plant_pulse_plasma_current_ramp_up = (
                physics_variables.plasma_current / 1.0e5
            )
            times_variables.t_plant_pulse_plasma_current_ramp_down = (
                times_variables.t_plant_pulse_plasma_current_ramp_up
            )

        else:
            # times_variables.t_plant_pulse_plasma_current_ramp_up is set either in INITIAL or INPUT, or by being
            # iterated using limit equation 41.
            times_variables.t_plant_pulse_coil_precharge = max(
                times_variables.t_plant_pulse_coil_precharge,
                times_variables.t_plant_pulse_plasma_current_ramp_up,
            )
            # t_plant_pulse_plasma_current_ramp_down = max(t_plant_pulse_plasma_current_ramp_down,t_plant_pulse_plasma_current_ramp_up)
            times_variables.t_plant_pulse_plasma_current_ramp_down = (
                times_variables.t_plant_pulse_plasma_current_ramp_up
            )

        # Reset second times_variables.t_plant_pulse_burn value (times_variables.t_burn_0).
        # This is used to ensure that the burn time is used consistently;
        # see convergence loop in fcnvmc1, evaluators.f90
        times_variables.t_burn_0 = times_variables.t_plant_pulse_burn

        # Time during the pulse in which a plasma is present
        times_variables.t_plant_pulse_plasma_present = (
            times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_fusion_ramp
            + times_variables.t_plant_pulse_burn
            + times_variables.t_plant_pulse_plasma_current_ramp_down
        )
        times_variables.t_plant_pulse_no_burn = (
            times_variables.t_plant_pulse_coil_precharge
            + times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_plasma_current_ramp_down
            + times_variables.t_plant_pulse_dwell
            + times_variables.t_plant_pulse_fusion_ramp
        )

        # Total cycle time
        times_variables.t_plant_pulse_total = (
            times_variables.t_plant_pulse_coil_precharge
            + times_variables.t_plant_pulse_plasma_current_ramp_up
            + times_variables.t_plant_pulse_fusion_ramp
            + times_variables.t_plant_pulse_burn
            + times_variables.t_plant_pulse_plasma_current_ramp_down
            + times_variables.t_plant_pulse_dwell
        )

        # ***************************** #
        #      DIAMAGNETIC CURRENT      #
        # ***************************** #

        self.dia_current.run()

        # ***************************** #
        #    PFIRSCH-SCHLÜTER CURRENT   #
        # ***************************** #

        # Pfirsch-Schlüter scaling for diamagnetic current
        current_drive_variables.f_c_plasma_pfirsch_schluter_scene = ps_fraction_scene(
            physics_variables.beta_total_vol_avg
        )

        if physics_variables.i_pfirsch_schluter_current == 1:
            current_drive_variables.f_c_plasma_pfirsch_schluter = (
                current_drive_variables.f_c_plasma_pfirsch_schluter_scene
            )

        self.plasma_bootstrap_current.run()

        physics_variables.err242 = 0
        if (
            current_drive_variables.f_c_plasma_bootstrap
            > current_drive_variables.f_c_plasma_bootstrap_max
        ) and physics_variables.i_bootstrap_current != 0:
            current_drive_variables.f_c_plasma_bootstrap = min(
                current_drive_variables.f_c_plasma_bootstrap,
                current_drive_variables.f_c_plasma_bootstrap_max,
            )
            physics_variables.err242 = 1

        current_drive_variables.f_c_plasma_internal = (
            current_drive_variables.f_c_plasma_bootstrap
            + current_drive_variables.f_c_plasma_diamagnetic
            + current_drive_variables.f_c_plasma_pfirsch_schluter
        )

        # Plasma driven current fraction (Bootstrap + Diamagnetic
        # + Pfirsch-Schlüter) constrained to be less than
        # or equal to the total fraction of the plasma current
        # produced by non-inductive means (which also includes
        # the current drive proportion)
        physics_variables.err243 = 0
        if (
            current_drive_variables.f_c_plasma_internal
            > physics_variables.f_c_plasma_non_inductive
        ):
            current_drive_variables.f_c_plasma_internal = min(
                current_drive_variables.f_c_plasma_internal,
                physics_variables.f_c_plasma_non_inductive,
            )
            physics_variables.err243 = 1

        # Fraction of plasma current produced by inductive means
        physics_variables.f_c_plasma_inductive = max(
            1.0e-10, (1.0e0 - physics_variables.f_c_plasma_non_inductive)
        )
        #  Fraction of plasma current produced by auxiliary current drive
        physics_variables.f_c_plasma_auxiliary = (
            physics_variables.f_c_plasma_non_inductive
            - current_drive_variables.f_c_plasma_internal
        )

        # Auxiliary current drive power calculations

        if current_drive_variables.i_hcd_calculations != 0:
            self.current_drive.cudriv()

        # ***************************** #
        #        FUSION REACTIONS       #
        # ***************************** #

        # Calculate fusion power + components

        fusion_reactions = reactions.FusionReactionRate(self.plasma_profile)
        fusion_reactions.deuterium_branching(
            physics_variables.temp_plasma_ion_vol_avg_kev
        )
        fusion_reactions.calculate_fusion_rates()
        fusion_reactions.set_physics_variables()

        # This neglects the power from the beam
        physics_variables.p_plasma_dt_mw = (
            physics_variables.dt_power_density_plasma * physics_variables.vol_plasma
        )
        physics_variables.p_dhe3_total_mw = (
            physics_variables.dhe3_power_density * physics_variables.vol_plasma
        )
        physics_variables.p_dd_total_mw = (
            physics_variables.dd_power_density * physics_variables.vol_plasma
        )

        # Calculate neutral beam slowing down effects
        # If ignited, then ignore beam fusion effects

        if (current_drive_variables.c_beam_total != 0.0e0) and (  # noqa: RUF069
            physics_variables.i_plasma_ignited == 0
        ):
            (
                physics_variables.beta_beam,
                physics_variables.nd_beam_ions_out,
                physics_variables.p_beam_alpha_mw,
            ) = reactions.beam_fusion(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.b_plasma_surface_poloidal_average,
                physics_variables.b_plasma_toroidal_on_axis,
                current_drive_variables.c_beam_total,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_fuel_ions_vol_avg,
                physics_variables.dlamie,
                current_drive_variables.e_beam_kev,
                physics_variables.f_plasma_fuel_deuterium,
                physics_variables.f_plasma_fuel_tritium,
                current_drive_variables.f_beam_tritium,
                physics_variables.sigmav_dt_average,
                physics_variables.temp_plasma_electron_density_weighted_kev,
                physics_variables.temp_plasma_ion_density_weighted_kev,
                physics_variables.vol_plasma,
                physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            )
            physics_variables.fusden_total = (
                physics_variables.fusden_plasma
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.fusden_alpha_total = (
                physics_variables.fusden_plasma_alpha
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.p_dt_total_mw = (
                physics_variables.p_plasma_dt_mw
                + (1.0 / (1.0 - constants.DT_NEUTRON_ENERGY_FRACTION))
                * physics_variables.p_beam_alpha_mw
            )
            physics_variables.p_beam_neutron_mw = physics_variables.p_beam_alpha_mw * (
                constants.DT_NEUTRON_ENERGY_FRACTION
                / (1 - constants.DT_NEUTRON_ENERGY_FRACTION)
            )
            physics_variables.p_beam_dt_mw = physics_variables.p_beam_alpha_mw * (
                1 / (1 - constants.DT_NEUTRON_ENERGY_FRACTION)
            )
        else:
            # If no beams present then the total alpha rates and power are the same as the plasma values
            physics_variables.fusden_total = physics_variables.fusden_plasma
            physics_variables.fusden_alpha_total = physics_variables.fusden_plasma_alpha
            physics_variables.p_dt_total_mw = physics_variables.p_plasma_dt_mw

        physics_variables.fusrat_total = (
            physics_variables.fusden_total * physics_variables.vol_plasma
        )

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.pden_neutron_total_mw,
            physics_variables.p_plasma_alpha_mw,
            physics_variables.p_alpha_total_mw,
            physics_variables.p_plasma_neutron_mw,
            physics_variables.p_neutron_total_mw,
            physics_variables.p_non_alpha_charged_mw,
            physics_variables.pden_alpha_total_mw,
            physics_variables.f_pden_alpha_electron_mw,
            physics_variables.f_pden_alpha_ions_mw,
            physics_variables.p_charged_particle_mw,
            physics_variables.p_fusion_total_mw,
        ) = reactions.set_fusion_powers(
            physics_variables.f_alpha_electron,
            physics_variables.f_alpha_ion,
            physics_variables.p_beam_alpha_mw,
            physics_variables.pden_non_alpha_charged_mw,
            physics_variables.pden_plasma_neutron_mw,
            physics_variables.vol_plasma,
            physics_variables.pden_plasma_alpha_mw,
        )

        physics_variables.beta_fast_alpha = self.beta.fast_alpha_beta(
            b_plasma_poloidal_average=physics_variables.b_plasma_surface_poloidal_average,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
            nd_plasma_fuel_ions_vol_avg=physics_variables.nd_plasma_fuel_ions_vol_avg,
            nd_plasma_ions_total_vol_avg=physics_variables.nd_plasma_ions_total_vol_avg,
            temp_plasma_electron_density_weighted_kev=physics_variables.temp_plasma_electron_density_weighted_kev,
            temp_plasma_ion_density_weighted_kev=physics_variables.temp_plasma_ion_density_weighted_kev,
            pden_alpha_total_mw=physics_variables.pden_alpha_total_mw,
            pden_plasma_alpha_mw=physics_variables.pden_plasma_alpha_mw,
            i_beta_fast_alpha=physics_variables.i_beta_fast_alpha,
        )

        # Calculate ion/electron equilibration power

        physics_variables.pden_ion_electron_equilibration_mw = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.dlamie,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.temp_plasma_ion_vol_avg_kev,
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
        )

        # Calculate radiation power

        radpwrdata = physics_funcs.calculate_radiation_powers(
            self.plasma_profile,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.rminor,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.aspect,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.tbeta,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.f_sync_reflect,
            physics_variables.rmajor,
            physics_variables.kappa,
            physics_variables.vol_plasma,
        )
        physics_variables.pden_plasma_sync_mw = radpwrdata.pden_plasma_sync_mw
        physics_variables.pden_plasma_core_rad_mw = radpwrdata.pden_plasma_core_rad_mw
        physics_variables.pden_plasma_outer_rad_mw = radpwrdata.pden_plasma_outer_rad_mw
        physics_variables.pden_plasma_rad_mw = radpwrdata.pden_plasma_rad_mw

        physics_variables.p_plasma_sync_mw = (
            physics_variables.pden_plasma_sync_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_inner_rad_mw = (
            physics_variables.pden_plasma_core_rad_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_outer_rad_mw = (
            physics_variables.pden_plasma_outer_rad_mw * physics_variables.vol_plasma
        )
        physics_variables.p_plasma_rad_mw = (
            physics_variables.pden_plasma_rad_mw * physics_variables.vol_plasma
        )

        # Calculate ohmic power
        (
            physics_variables.pden_plasma_ohmic_mw,
            physics_variables.p_plasma_ohmic_mw,
            physics_variables.f_res_plasma_neo,
            physics_variables.res_plasma,
        ) = self.plasma_ohmic_heating(
            physics_variables.f_c_plasma_inductive,
            physics_variables.kappa95,
            physics_variables.plasma_current,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.vol_plasma,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        # Calculate L- to H-mode power threshold for different scalings
        self.plasma_transition.run()

        # Power transported to the divertor by charged particles,
        # i.e. excludes neutrons and radiation, and also NBI orbit loss power,
        # which is assumed to be absorbed by the first wall
        pinj = (
            current_drive_variables.p_hcd_injected_total_mw
            if physics_variables.i_plasma_ignited == 0
            else 0.0
        )

        physics_variables.p_plasma_separatrix_mw = (
            self.exhaust.calculate_separatrix_power(
                f_p_alpha_plasma_deposited=physics_variables.f_p_alpha_plasma_deposited,
                p_alpha_total_mw=physics_variables.p_alpha_total_mw,
                p_non_alpha_charged_mw=physics_variables.p_non_alpha_charged_mw,
                p_hcd_injected_total_mw=pinj,
                p_plasma_ohmic_mw=physics_variables.p_plasma_ohmic_mw,
                p_plasma_rad_mw=physics_variables.p_plasma_rad_mw,
            )
        )

        physics_variables.p_plasma_separatrix_rmajor_mw = (
            self.exhaust.calculate_psep_over_r_metric(
                p_plasma_separatrix_mw=physics_variables.p_plasma_separatrix_mw,
                rmajor=physics_variables.rmajor,
            )
        )

        physics_variables.p_div_bt_q_aspect_rmajor_mw = (
            self.exhaust.calculate_eu_demo_re_attachment_metric(
                p_plasma_separatrix_mw=physics_variables.p_plasma_separatrix_mw,
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                q95=physics_variables.q95,
                aspect=physics_variables.aspect,
                rmajor=physics_variables.rmajor,
            )
        )

        physics_variables.pflux_plasma_surface_neutron_avg_mw = (
            physics_variables.p_neutron_total_mw / physics_variables.a_plasma_surface
        )

        # KLUDGE: Ensure p_plasma_separatrix_mw is continuously positive (physical, rather than
        # negative potential power), as required by other models (e.g.
        # Physics.calculate_density_limit())
        physics_variables.p_plasma_separatrix_mw /= 1 - np.exp(
            -physics_variables.p_plasma_separatrix_mw
        )

        # if double null configuration share the power
        # over the upper and lower divertor, where physics_variables.f_p_div_lower gives
        # the factor of power conducted to the lower divertor
        if divertor_variables.n_divertors == 2:
            physics_variables.p_div_lower_separatrix_mw = (
                physics_variables.f_p_div_lower
                * physics_variables.p_plasma_separatrix_mw
            )
            physics_variables.p_div_upper_separatrix_mw = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_variables.p_plasma_separatrix_mw
            physics_variables.p_div_separatrix_max_mw = max(
                physics_variables.p_div_lower_separatrix_mw,
                physics_variables.p_div_upper_separatrix_mw,
            )

        # Resistive diffusion time = current penetration time ~ mu0.a^2/resistivity
        physics_variables.t_plasma_res_diffusion = res_diff_time(
            physics_variables.rmajor,
            physics_variables.res_plasma,
            physics_variables.kappa95,
        )

        self.density_limit.run()

        # Calculate transport losses and energy confinement time using the
        # chosen scaling law
        (
            physics_variables.pden_electron_transport_loss_mw,
            physics_variables.pden_ion_transport_loss_mw,
            physics_variables.t_electron_energy_confinement,
            physics_variables.t_energy_confinement,
            physics_variables.t_ion_energy_confinement,
            physics_variables.p_plasma_loss_mw,
        ) = self.confinement.calculate_confinement_time(
            m_fuel_amu=physics_variables.m_fuel_amu,
            p_alpha_total_mw=physics_variables.p_alpha_total_mw,
            aspect=physics_variables.aspect,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            nd_plasma_ions_total_vol_avg=physics_variables.nd_plasma_ions_total_vol_avg,
            nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
            nd_plasma_electron_line=physics_variables.nd_plasma_electron_line,
            eps=physics_variables.eps,
            hfact=physics_variables.hfact,
            i_confinement_time=physics_variables.i_confinement_time,
            i_plasma_ignited=physics_variables.i_plasma_ignited,
            kappa=physics_variables.kappa,
            kappa95=physics_variables.kappa95,
            p_non_alpha_charged_mw=physics_variables.p_non_alpha_charged_mw,
            p_hcd_injected_total_mw=current_drive_variables.p_hcd_injected_total_mw,
            plasma_current=physics_variables.plasma_current,
            pden_plasma_core_rad_mw=physics_variables.pden_plasma_core_rad_mw,
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            temp_plasma_electron_density_weighted_kev=physics_variables.temp_plasma_electron_density_weighted_kev,
            temp_plasma_ion_density_weighted_kev=physics_variables.temp_plasma_ion_density_weighted_kev,
            q95=physics_variables.q95,
            qstar=physics_variables.qstar,
            vol_plasma=physics_variables.vol_plasma,
            zeff=physics_variables.n_charge_plasma_effective_vol_avg,
        )

        # Total transport power from scaling law (MW)
        physics_variables.p_electron_transport_loss_mw = (
            physics_variables.pden_electron_transport_loss_mw
            * physics_variables.vol_plasma
        )
        physics_variables.p_ion_transport_loss_mw = (
            physics_variables.pden_ion_transport_loss_mw * physics_variables.vol_plasma
        )

        # Calculate Volt-second requirements
        (
            physics_variables.vs_plasma_internal,
            physics_variables.ind_plasma,
            physics_variables.vs_plasma_burn_required,
            physics_variables.vs_plasma_ramp_required,
            physics_variables.vs_plasma_ind_ramp,
            physics_variables.vs_plasma_res_ramp,
            physics_variables.vs_plasma_total_required,
            physics_variables.v_plasma_loop_burn,
        ) = self.inductance.calculate_volt_second_requirements(
            physics_variables.csawth,
            physics_variables.eps,
            physics_variables.f_c_plasma_inductive,
            physics_variables.ejima_coeff,
            physics_variables.kappa,
            physics_variables.rmajor,
            physics_variables.res_plasma,
            physics_variables.plasma_current,
            times_variables.t_plant_pulse_fusion_ramp,
            times_variables.t_plant_pulse_burn,
            physics_variables.ind_plasma_internal_norm,
        )

        physics_variables.e_plasma_magnetic_stored = (
            0.5e0 * physics_variables.ind_plasma * physics_variables.plasma_current**2
        )

        # Calculate auxiliary physics related information
        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.figmer,
            physics_variables.fusrat,
            physics_variables.molflow_plasma_fuelling_required,
            physics_variables.rndfuel,
            physics_variables.t_alpha_confinement,
            physics_variables.f_alpha_energy_confinement,
        ) = self.phyaux(
            physics_variables.aspect,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.fusden_total,
            physics_variables.fusden_alpha_total,
            physics_variables.plasma_current,
            sbar,
            physics_variables.nd_plasma_alphas_vol_avg,
            physics_variables.t_energy_confinement,
            physics_variables.vol_plasma,
        )

        physics_variables.ntau, physics_variables.nTtau = (
            self.confinement.calculate_double_and_triple_product(
                nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                t_energy_confinement=physics_variables.t_energy_confinement,
                temp_plasma_electrons_vol_avg_kev=physics_variables.temp_plasma_electron_vol_avg_kev,
            )
        )

        # Total transport power from scaling law (MW)
        physics_variables.pscalingmw = (
            physics_variables.p_electron_transport_loss_mw
            + physics_variables.p_ion_transport_loss_mw
        )

        # ============================================================

        # Calculate the target imbalances
        # find the total power into the targets
        physics_variables.ptarmw = physics_variables.p_plasma_separatrix_mw * (
            1.0e0 - physics_variables.rad_fraction_sol
        )
        # use physics_variables.f_p_div_lower to find deltarsep
        # Parameters taken from double null machine
        # D. Brunner et al
        physics_variables.lambdaio = 1.57e-3

        # Issue #1559 Infinities in physics_module.drsep when running single null in a double null machine
        # C W Ashe
        if physics_variables.f_p_div_lower < 4.5e-5:
            physics_variables.drsep = 1.5e-2
        elif physics_variables.f_p_div_lower > (1.0e0 - 4.5e-5):
            physics_variables.drsep = -1.5e-2
        else:
            physics_variables.drsep = (
                -2.0e0
                * 1.5e-3
                * math.atanh(2.0e0 * (physics_variables.f_p_div_lower - 0.5e0))
            )
        # Model Taken from D3-D paper for conventional divertor
        # Journal of Nuclear Materials
        # Volumes 290-293, March 2001, Pages 935-939
        # Find the innner and outer lower target imbalance

        physics_variables.fio = 0.16e0 + (0.16e0 - 0.41e0) * (
            1.0e0
            - (
                2.0e0
                / (
                    1.0e0
                    + np.exp(
                        -((physics_variables.drsep / physics_variables.lambdaio) ** 2)
                    )
                )
            )
        )
        if divertor_variables.n_divertors == 2:
            # Double Null configuration
            # Find all the power fractions accross the targets
            # Taken from D3-D conventional divertor design
            physics_variables.fli = (
                physics_variables.f_p_div_lower * physics_variables.fio
            )
            physics_variables.flo = physics_variables.f_p_div_lower * (
                1.0e0 - physics_variables.fio
            )
            physics_variables.fui = (
                1.0e0 - physics_variables.f_p_div_lower
            ) * physics_variables.fio
            physics_variables.fuo = (1.0e0 - physics_variables.f_p_div_lower) * (
                1.0e0 - physics_variables.fio
            )
            # power into each target
            physics_variables.plimw = physics_variables.fli * physics_variables.ptarmw
            physics_variables.plomw = physics_variables.flo * physics_variables.ptarmw
            physics_variables.puimw = physics_variables.fui * physics_variables.ptarmw
            physics_variables.puomw = physics_variables.fuo * physics_variables.ptarmw
        else:
            # Single null configuration
            physics_variables.fli = physics_variables.fio
            physics_variables.flo = 1.0e0 - physics_variables.fio
            # power into each target
            physics_variables.plimw = physics_variables.fli * physics_variables.ptarmw
            physics_variables.plomw = physics_variables.flo * physics_variables.ptarmw

        # Calculate some derived quantities that may not have been defined earlier
        physics_variables.p_plasma_heating_total_mw = 1e6 * (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            + current_drive_variables.p_hcd_injected_total_mw
        )
        physics_variables.f_p_plasma_separatrix_rad = (
            1.0e6
            * physics_variables.p_plasma_rad_mw
            / physics_variables.p_plasma_heating_total_mw
        )
        physics_variables.rad_fraction_total = (
            physics_variables.f_p_plasma_separatrix_rad
            + (1.0e0 - physics_variables.f_p_plasma_separatrix_rad)
            * physics_variables.rad_fraction_sol
        )
        physics_variables.pradsolmw = (
            physics_variables.rad_fraction_sol * physics_variables.p_plasma_separatrix_mw
        )

        if 78 in numerics.icc:
            po.write(
                self.outfile,
                (
                    f"reinke t and fz, physics = {physics_variables.temp_plasma_separatrix_kev} , {reinke_variables.fzmin}"
                ),
            )
            fgw = (
                physics_variables.nd_plasma_electron_max_array[6]
                / physics_variables.nd_plasma_electrons_vol_avg
            )
            # calculate separatrix temperature, if Reinke criterion is used
            physics_variables.temp_plasma_separatrix_kev = reinke_tsep(
                physics_variables.b_plasma_toroidal_on_axis,
                physics_variables.p_plasma_separatrix_mw
                / physics_variables.p_l_h_threshold_mw,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.eps,
                fgw,
                physics_variables.kappa,
                reinke_variables.lhat,
            )

            if reinke_variables.fzmin >= 1.0e0:
                logger.error(
                    "REINKE IMPURITY MODEL: fzmin is greater than or equal to 1.0, this is at least notable"
                )

            po.write(
                self.outfile,
                (
                    f" 'fzactual, frac, reinke_variables.impvardiv = {reinke_variables.fzactual},"
                    f" {impurity_radiation_module.f_nd_impurity_electron_array(reinke_variables.impvardiv)},"
                    f" {reinke_variables.impvardiv}"
                ),
            )

    @staticmethod
    def calculate_current_profile_index_wesson(qstar: float, q0: float) -> float:
        """Calculate the Wesson current profile index.

        Parameters
        ----------
        qstar : float
            Cylindrical safety factor.
        q0 : float
            Safety factor on axis.

        Returns
        -------
        float
            The Wesson current profile index.

        Notes
        -----
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              normalised internal inductance.
            - This relation is only true for the cyclindrical plasma approximation.

        References
        ----------
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.

        """
        return qstar / q0 - 1.0

    @staticmethod
    def plasma_composition():
        """Calculates various plasma component fractional makeups.

        This subroutine determines the various plasma component fractional makeups.
        It is the replacement for the original routine betcom(), and is used in conjunction
        with the new impurity radiation model.

        This function performs the following calculations:
        - Determines the alpha ash portion.
        - Calculates the proton density.
        - Calculates the beam hot ion component.
        - Sums the ion densities for all impurity ions with charge greater than helium.
        - Ensures charge neutrality by adjusting the fuel portion.
        - Calculates the total ion density.
        - Sets the relative impurity densities for radiation calculations.
        - Calculates the effective charge.
        - Defines the Coulomb logarithm for ion-electron and electron-electron collisions.
        - Calculates the fraction of alpha energy to ions and electrons.
        - Calculates the average atomic masses of injected fuel species and neutral beams.
        - Calculates the density weighted mass and mass weighted plasma effective charge.
        """
        # Alpha ash portion
        physics_variables.nd_plasma_alphas_vol_avg = (
            physics_variables.nd_plasma_electrons_vol_avg
            * physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Protons
        # This calculation will be wrong on the first call as the particle
        # production rates are evaluated later in the calling sequence
        # Issue #557 Allow f_nd_protium_electrons impurity to be specified: 'f_nd_protium_electrons'
        # This will override the calculated value which is a minimum.
        if physics_variables.fusden_alpha_total < 1.0e-6:  # not calculated yet...
            physics_variables.nd_plasma_protons_vol_avg = max(
                physics_variables.f_nd_protium_electrons
                * physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_alphas_vol_avg
                * (physics_variables.f_plasma_fuel_helium3 + 1.0e-3),
            )  # rough estimate
        else:
            physics_variables.nd_plasma_protons_vol_avg = max(
                physics_variables.f_nd_protium_electrons
                * physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_alphas_vol_avg
                * physics_variables.proton_rate_density
                / physics_variables.fusden_alpha_total,
            )

        # ======================================================================

        # Beam hot ion component
        # If ignited, prevent beam fusion effects
        if physics_variables.i_plasma_ignited == 0:
            physics_variables.nd_beam_ions = (
                physics_variables.nd_plasma_electrons_vol_avg
                * physics_variables.f_nd_beam_electron
            )
        else:
            physics_variables.nd_beam_ions = 0.0

        # ======================================================================

        # Sum of Zi.ni for all impurity ions (those with charge > helium)
        znimp = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                znimp += impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.temp_plasma_electron_vol_avg_kev])
                ).squeeze() * (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * physics_variables.nd_plasma_electrons_vol_avg
                )

        # ======================================================================

        # Fuel portion - conserve charge neutrality
        # znfuel is the sum of Zi.ni for the three fuel ions
        znfuel = (
            physics_variables.nd_plasma_electrons_vol_avg
            - 2.0 * physics_variables.nd_plasma_alphas_vol_avg
            - physics_variables.nd_plasma_protons_vol_avg
            - physics_variables.nd_beam_ions
            - znimp
        )

        # ======================================================================

        # Fuel ion density, nd_plasma_fuel_ions_vol_avg
        # For D-T-He3 mix, nd_plasma_fuel_ions_vol_avg = nD + nT + nHe3, while znfuel = nD + nT + 2*nHe3
        # So nd_plasma_fuel_ions_vol_avg = znfuel - nHe3 = znfuel - f_plasma_fuel_helium3*nd_plasma_fuel_ions_vol_avg
        physics_variables.nd_plasma_fuel_ions_vol_avg = znfuel / (
            1.0 + physics_variables.f_plasma_fuel_helium3
        )

        # ======================================================================

        # Set hydrogen and helium relative impurity densities for
        # radiation calculations
        impurity_radiation_module.f_nd_impurity_electron_array[
            impurity_radiation.element2index("H_")
        ] = (
            physics_variables.nd_plasma_protons_vol_avg
            + (
                physics_variables.f_plasma_fuel_deuterium
                + physics_variables.f_plasma_fuel_tritium
            )
            * physics_variables.nd_plasma_fuel_ions_vol_avg
            + physics_variables.nd_beam_ions
        ) / physics_variables.nd_plasma_electrons_vol_avg

        impurity_radiation_module.f_nd_impurity_electron_array[
            impurity_radiation.element2index("He")
        ] = (
            physics_variables.f_plasma_fuel_helium3
            * physics_variables.nd_plasma_fuel_ions_vol_avg
            / physics_variables.nd_plasma_electrons_vol_avg
            + physics_variables.f_nd_alpha_electron
        )

        # ======================================================================

        # Total impurity density
        physics_variables.nd_plasma_impurities_vol_avg = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.nd_plasma_impurities_vol_avg += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * physics_variables.nd_plasma_electrons_vol_avg
                )

        # ======================================================================

        # Total ion density
        physics_variables.nd_plasma_ions_total_vol_avg = (
            physics_variables.nd_plasma_fuel_ions_vol_avg
            + physics_variables.nd_plasma_alphas_vol_avg
            + physics_variables.nd_plasma_protons_vol_avg
            + physics_variables.nd_beam_ions
            + physics_variables.nd_plasma_impurities_vol_avg
        )

        # ======================================================================

        # Set some relative impurity density variables
        # for the benefit of other routines
        physics_variables.f_nd_plasma_carbon_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("C_")
            ]
        )
        physics_variables.f_nd_plasma_oxygen_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("O_")
            ]
        )
        physics_variables.f_nd_plasma_iron_argon_electron = (
            impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("Fe")
            ]
            + impurity_radiation_module.f_nd_impurity_electron_array[
                impurity_radiation.element2index("Ar")
            ]
        )

        # ======================================================================

        # Effective charge
        # Calculation should be sum(ni.Zi^2) / sum(ni.Zi),
        # but ne = sum(ni.Zi) through quasineutrality
        physics_variables.n_charge_plasma_effective_vol_avg = 0.0
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            physics_variables.n_charge_plasma_effective_vol_avg += (
                impurity_radiation_module.f_nd_impurity_electron_array[imp]
                * impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.temp_plasma_electron_vol_avg_kev])
                ).squeeze()
                ** 2
            )

        # ======================================================================

        # Fraction of alpha energy to ions and electrons
        # From Max Fenstermacher
        # (used with electron and ion power balance equations only)
        # No consideration of pden_non_alpha_charged_mw here...

        # f_temp_plasma_electron_density_vol_avg now calculated in plasma_profiles, after the very first
        # call of plasma_composition; use old parabolic profile estimate
        # in this case
        if physics_variables.first_call == 1:
            pc = (
                (1.0 + physics_variables.alphan)
                * (1.0 + physics_variables.alphat)
                / (1.0 + physics_variables.alphan + physics_variables.alphat)
            )
            physics_variables.first_call = 0
        else:
            pc = physics_variables.f_temp_plasma_electron_density_vol_avg

        physics_variables.f_alpha_electron = 0.88155 * np.exp(
            -physics_variables.temp_plasma_electron_vol_avg_kev * pc / 67.4036
        )
        physics_variables.f_alpha_ion = 1.0 - physics_variables.f_alpha_electron

        # ======================================================================

        # Average atomic masses of injected fuel species
        physics_variables.m_fuel_amu = (
            (constants.M_DEUTERON_AMU * physics_variables.f_plasma_fuel_deuterium)
            + (constants.M_TRITON_AMU * physics_variables.f_plasma_fuel_tritium)
            + (constants.M_HELION_AMU * physics_variables.f_plasma_fuel_helium3)
        )

        # ======================================================================

        # Average atomic masses of injected fuel species in the neutral beams
        # Only deuterium and tritium in the beams
        physics_variables.m_beam_amu = (
            constants.M_DEUTERON_AMU * (1.0 - current_drive_variables.f_beam_tritium)
        ) + (constants.M_TRITON_AMU * current_drive_variables.f_beam_tritium)

        # ======================================================================

        # Average mass of all ions
        physics_variables.m_ions_total_amu = (
            (
                physics_variables.m_fuel_amu
                * physics_variables.nd_plasma_fuel_ions_vol_avg
            )
            + (constants.M_ALPHA_AMU * physics_variables.nd_plasma_alphas_vol_avg)
            + (physics_variables.nd_plasma_protons_vol_avg * constants.M_PROTON_AMU)
            + (physics_variables.m_beam_amu * physics_variables.nd_beam_ions)
        )
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.m_ions_total_amu += (
                    physics_variables.nd_plasma_electrons_vol_avg
                    * impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation_module.m_impurity_amu_array[imp]
                )

        physics_variables.m_ions_total_amu /= (
            physics_variables.nd_plasma_ions_total_vol_avg
        )

        # ======================================================================

        # Mass weighted plasma effective charge
        # Sum of (Zi^2*n_i) / m_i
        physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg = (
            (
                physics_variables.f_plasma_fuel_deuterium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_DEUTERON_AMU
            )
            + (
                physics_variables.f_plasma_fuel_tritium
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_TRITON_AMU
            )
            + (
                4.0
                * physics_variables.f_plasma_fuel_helium3
                * physics_variables.nd_plasma_fuel_ions_vol_avg
                / constants.M_HELION_AMU
            )
            + (4.0 * physics_variables.nd_plasma_alphas_vol_avg / constants.M_ALPHA_AMU)
            + (physics_variables.nd_plasma_protons_vol_avg / constants.M_PROTON_AMU)
            + (
                (1.0 - current_drive_variables.f_beam_tritium)
                * physics_variables.nd_beam_ions
                / constants.M_DEUTERON_AMU
            )
            + (
                current_drive_variables.f_beam_tritium
                * physics_variables.nd_beam_ions
                / constants.M_TRITON_AMU
            )
        ) / physics_variables.nd_plasma_electrons_vol_avg
        for imp in range(impurity_radiation_module.N_IMPURITIES):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation.zav_of_te(
                        imp,
                        np.array([physics_variables.temp_plasma_electron_vol_avg_kev]),
                    ).squeeze()
                    ** 2
                    / impurity_radiation_module.m_impurity_amu_array[imp]
                )

        # ======================================================================

    @staticmethod
    @nb.njit(cache=True)
    def phyaux(
        aspect: float,
        nd_plasma_fuel_ions_vol_avg: float,
        fusden_total: float,
        fusden_alpha_total: float,
        plasma_current: float,
        sbar: float,
        nd_plasma_alphas_vol_avg: float,
        t_energy_confinement: float,
        vol_plasma: float,
    ) -> tuple[float, float, float, float, float, float, float, float]:
        """Auxiliary physics quantities

        Parameters
        ----------
        aspect : float
            Plasma aspect ratio.
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (/m3).
        fusden_total : float
            Fusion reaction rate from plasma and beams (/m3/s).
        fusden_alpha_total : float
            Alpha particle production rate (/m3/s).
        plasma_current : float
            Plasma current (A).
        sbar : float
            Exponent for aspect ratio (normally 1).
        nd_plasma_alphas_vol_avg : float
            Alpha ash density (/m3).
        t_energy_confinement : float
            Global energy confinement time (s).
        vol_plasma : float
            Plasma volume (m3).

        Returns
        -------
        tuple
            A tuple containing:
            - burnup (float): Fractional plasma burnup.
            - figmer (float): Physics figure of merit.
            - fusrat (float): Number of fusion reactions per second.
            - molflow_plasma_fuelling_required (float): Fuelling rate for D-T (nucleus-pairs/sec).
            - rndfuel (float): Fuel burnup rate (reactions/s).
            - t_alpha_confinement (float): Alpha particle confinement time (s).
            - f_alpha_energy_confinement (float): Fraction of alpha energy confinement.
            This subroutine calculates extra physics related items needed by other parts of the code.

        """
        figmer = 1e-6 * plasma_current * aspect**sbar

        # Fusion reactions per second
        fusrat = fusden_total * vol_plasma

        # Alpha particle confinement time (s)
        # Number of alphas / alpha production rate
        # only likely if DD is only active fusion reaction
        t_alpha_confinement = (
            0.0
            if fusden_alpha_total == 0.0  # noqa: RUF069
            else nd_plasma_alphas_vol_avg / fusden_alpha_total
        )

        # Fractional burnup
        # (Consider detailed model in: G. L. Jackson, V. S. Chan, R. D. Stambaugh,
        # Fusion Science and Technology, vol.64, no.1, July 2013, pp.8-12)
        # The ratio of ash to fuel particle confinement times is given by
        # tauratio
        # Possible logic...
        # burnup = fuel ion-pairs burned/m3 / initial fuel ion-pairs/m3;
        # fuel ion-pairs burned/m3 = alpha particles/m3 (for both D-T and D-He3 reactions)
        # initial fuel ion-pairs/m3 = burnt fuel ion-pairs/m3 + unburnt fuel-ion pairs/m3
        # Remember that unburnt fuel-ion pairs/m3 = 0.5 * unburnt fuel-ions/m3
        if physics_variables.burnup_in <= 1.0e-9:
            burnup = (
                nd_plasma_alphas_vol_avg
                / (nd_plasma_alphas_vol_avg + 0.5 * nd_plasma_fuel_ions_vol_avg)
                / physics_variables.tauratio
            )
        else:
            burnup = physics_variables.burnup_in

        # Fuel burnup rate (reactions/second) (previously Amps)
        rndfuel = fusrat

        # Required fuelling rate (fuel ion pairs/second) (previously Amps)
        molflow_plasma_fuelling_required = rndfuel / burnup

        f_alpha_energy_confinement = t_alpha_confinement / t_energy_confinement

        return (
            burnup,
            figmer,
            fusrat,
            molflow_plasma_fuelling_required,
            rndfuel,
            t_alpha_confinement,
            f_alpha_energy_confinement,
        )

    @staticmethod
    def plasma_ohmic_heating(
        f_c_plasma_inductive: float,
        kappa95: float,
        plasma_current: float,
        rmajor: float,
        rminor: float,
        temp_plasma_electron_density_weighted_kev: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float]:
        """Calculate the ohmic heating power and related parameters.

        Parameters
        ----------
        f_c_plasma_inductive : float
            Fraction of plasma current driven inductively.
        kappa95 : float
            Plasma elongation at 95% surface.
        plasma_current : float
            Plasma current (A).
        rmajor : float
            Major radius (m).
        rminor : float
            Minor radius (m).
        temp_plasma_electron_density_weighted_kev : float
            Density weighted average electron temperature (keV).
        vol_plasma : float
            Plasma volume (m^3).
        zeff : float
            Plasma effective charge.

        Returns
        -------
        Tuple[float, float, float, float]
            Tuple containing:
            - pden_plasma_ohmic_mw (float): Ohmic heating power per unit volume (MW/m^3).
            - p_plasma_ohmic_mw (float): Total ohmic heating power (MW).
            - f_res_plasma_neo (float): Neo-classical resistivity enhancement factor.
            - res_plasma (float): Plasma resistance (ohm).

        References
        ----------
        - ITER Physics Design Guidelines
            1989 [IPDG89], N. A. Uckan et al,

        """
        # Density weighted electron temperature in 10 keV units
        t10 = temp_plasma_electron_density_weighted_kev / 10.0

        # Plasma resistance, from loop voltage calculation in ITER Physics Design Guidelines: 1989
        res_plasma = (
            physics_variables.plasma_res_factor
            * 2.15e-9
            * zeff
            * rmajor
            / (kappa95 * rminor**2 * t10**1.5)
        )

        # Neo-classical resistivity enhancement factor
        # Taken from ITER Physics Design Guidelines: 1989
        # The expression is valid for aspect ratios in the range 2.5 to 4.0

        f_res_plasma_neo = (
            1.0 if 2.5 >= rmajor / rminor <= 4.0 else 4.3 - 0.6 * rmajor / rminor
        )

        res_plasma *= f_res_plasma_neo

        # Check to see if plasma resistance is negative
        # (possible if aspect ratio is too high)
        if res_plasma <= 0.0:
            logger.error(
                f"Negative plasma resistance res_plasma. {res_plasma=} {physics_variables.aspect=}"
            )

        # Ohmic heating power per unit volume
        # Corrected from: pden_plasma_ohmic_mw = (f_c_plasma_inductive*plasma_current)**2 * ...

        pden_plasma_ohmic_mw = (
            f_c_plasma_inductive * plasma_current**2 * res_plasma * 1.0e-6 / vol_plasma
        )

        # Total ohmic heating power
        p_plasma_ohmic_mw = pden_plasma_ohmic_mw * vol_plasma

        return pden_plasma_ohmic_mw, p_plasma_ohmic_mw, f_res_plasma_neo, res_plasma

    def outtim(self):
        po.oheadr(self.outfile, "Times")

        po.ovarrf(
            self.outfile,
            "Initial charge time for CS from zero current (s)",
            "(t_plant_pulse_coil_precharge)",
            times_variables.t_plant_pulse_coil_precharge,
        )
        po.ovarrf(
            self.outfile,
            "Plasma current ramp-up time (s)",
            "(t_plant_pulse_plasma_current_ramp_up)",
            times_variables.t_plant_pulse_plasma_current_ramp_up,
        )
        po.ovarrf(
            self.outfile,
            "Heating time (s)",
            "(t_plant_pulse_fusion_ramp)",
            times_variables.t_plant_pulse_fusion_ramp,
        )
        po.ovarre(
            self.outfile,
            "Burn time (s)",
            "(t_plant_pulse_burn)",
            times_variables.t_plant_pulse_burn,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Reset time to zero current for CS (s)",
            "(t_plant_pulse_plasma_current_ramp_down)",
            times_variables.t_plant_pulse_plasma_current_ramp_down,
        )
        po.ovarrf(
            self.outfile,
            "Time between pulses (s)",
            "(t_plant_pulse_dwell)",
            times_variables.t_plant_pulse_dwell,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plant cycle time (s)",
            "(t_plant_pulse_total)",
            times_variables.t_plant_pulse_total,
            "OP ",
        )

    def calculate_effective_charge_ionisation_profiles(self):
        """Calculate the effective charge profiles for ionisation calculations."""
        # Calculate the effective charge (zeff) profile across the plasma
        # Returns an array of zeff at each radial point
        zeff_profile = np.zeros_like(self.plasma_profile.teprofile.profile_y)
        for i in range(len(zeff_profile)):
            zeff_profile[i] = 0.0
            for imp in range(impurity_radiation_module.N_IMPURITIES):
                zeff_profile[i] += (
                    impurity_radiation_module.f_nd_impurity_electron_array[imp]
                    * impurity_radiation.zav_of_te(
                        imp, np.array([self.plasma_profile.teprofile.profile_y[i]])
                    ).squeeze()
                    ** 2
                )
        physics_variables.n_charge_plasma_effective_profile = zeff_profile

        # Assign the charge profiles of each species
        n_impurities = impurity_radiation_module.N_IMPURITIES
        te_profile = self.plasma_profile.teprofile.profile_y
        n_points = len(te_profile)
        # Create a 2D array: (n_impurities, n_points)
        charge_profiles = np.zeros((n_impurities, n_points))
        for imp in range(n_impurities):
            for i in range(n_points):
                charge_profiles[imp, i] = impurity_radiation.zav_of_te(
                    imp, np.array([te_profile[i]])
                ).squeeze()
        impurity_radiation_module.n_charge_impurity_profile = charge_profiles

    def outplas(self):
        """Subroutine to output the plasma physics information
        self.outfile : input integer : Fortran output unit identifier
        This routine writes the plasma physics information
        to a file, in a tidy format.
        """
        # Dimensionless plasma parameters. See reference below.
        physics_variables.nu_star = (
            1
            / constants.RMU0
            * (15.0e0 * constants.ELECTRON_CHARGE**4 * physics_variables.dlamie)
            / (4.0e0 * np.pi**1.5e0 * constants.EPSILON0**2)
            * physics_variables.vol_plasma**2
            * physics_variables.rmajor**2
            * physics_variables.b_plasma_toroidal_on_axis
            * np.sqrt(physics_variables.eps)
            * physics_variables.nd_plasma_electron_line**3
            * physics_variables.kappa
            / (physics_variables.e_plasma_beta**2 * physics_variables.plasma_current)
        )

        physics_variables.rho_star = np.sqrt(
            2.0e0
            * constants.PROTON_MASS
            * physics_variables.m_ions_total_amu
            * physics_variables.e_plasma_beta
            / (
                3.0e0
                * physics_variables.vol_plasma
                * physics_variables.nd_plasma_electron_line
            )
        ) / (
            constants.ELECTRON_CHARGE
            * physics_variables.b_plasma_toroidal_on_axis
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_variables.beta_mcdonald = (
            4.0e0
            / 3.0e0
            * constants.RMU0
            * physics_variables.e_plasma_beta
            / (
                physics_variables.vol_plasma
                * physics_variables.b_plasma_toroidal_on_axis**2
            )
        )

        self.geometry.output()

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "Current and Field :")

            if stellarator_variables.istell == 0:
                po.oblnkl(self.outfile)
                po.ovarin(
                    self.outfile,
                    "Plasma current scaling law used",
                    "(i_plasma_current)",
                    physics_variables.i_plasma_current,
                )

                po.ovarrf(
                    self.outfile,
                    "Plasma current (MA)",
                    "(plasma_current_MA)",
                    physics_variables.plasma_current / 1.0e6,
                    "OP ",
                )

                self.current.output_plasma_current_models()
                po.oblnkl(self.outfile)

                if physics_variables.i_alphaj == 1:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                        "OP ",
                    )
                else:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                    )
                po.ocmmnt(self.outfile, "Current profile index scalings:")
                po.oblnkl(self.outfile)

                po.ovarrf(
                    self.outfile,
                    "J. Wesson plasma current profile index",
                    "(alphaj_wesson)",
                    physics_variables.alphaj_wesson,
                    "OP ",
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "On-axis plasma current density (A/m2)",
                    "(j_plasma_on_axis)",
                    physics_variables.j_plasma_on_axis,
                    "OP ",
                )

        if stellarator_variables.istell == 0:
            po.ovarrf(
                self.outfile, "Safety factor on axis", "(q0)", physics_variables.q0
            )

            if physics_variables.i_plasma_current == 2:
                po.ovarrf(
                    self.outfile,
                    "Mean edge safety factor",
                    "(q95)",
                    physics_variables.q95,
                )

            po.ovarrf(
                self.outfile,
                "Safety factor at 95% flux surface",
                "(q95)",
                physics_variables.q95,
            )

            po.ovarrf(
                self.outfile,
                "Cylindrical safety factor (qcyl)",
                "(qstar)",
                physics_variables.qstar,
                "OP ",
            )

            if physics_variables.i_plasma_geometry == 1:
                po.ovarrf(
                    self.outfile,
                    "Lower limit for edge safety factor q95",
                    "(q95_min)",
                    physics_variables.q95_min,
                    "OP ",
                )
            po.ovarrf(
                self.outfile,
                "Plasma normalised internal inductance",
                "(ind_plasma_internal_norm)",
                physics_variables.ind_plasma_internal_norm,
                "OP ",
            )

            po.oblnkl(self.outfile)
            self.fields.output()
            po.oblnkl(self.outfile)
            po.ostars(self.outfile, 110)
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Plasma normalised internal inductance scalings:")
            po.oblnkl(self.outfile)
            po.ovarrf(
                self.outfile,
                "J. Wesson plasma normalised internal inductance",
                "(ind_plasma_internal_norm_wesson)",
                physics_variables.ind_plasma_internal_norm_wesson,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "J. Menard plasma normalised internal inductance",
                "(ind_plasma_internal_norm_menard)",
                physics_variables.ind_plasma_internal_norm_menard,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "ITER li(3) plasma normalised internal inductance",
                "(ind_plasma_internal_norm_iter_3)",
                physics_variables.ind_plasma_internal_norm_iter_3,
                "OP ",
            )

        else:
            po.ovarrf(
                self.outfile,
                "Rotational transform",
                "(iotabar)",
                stellarator_variables.iotabar,
            )
        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        # Output beta information
        self.beta.output_beta_information()

        po.osubhd(self.outfile, "Temperature and Density (volume averaged) :")
        po.ovarrf(
            self.outfile,
            "Number of radial points in plasma profiles",
            "(n_plasma_profile_elements)",
            physics_variables.n_plasma_profile_elements,
        )
        po.ovarrf(
            self.outfile,
            "Volume averaged electron temperature (keV)",
            "(temp_plasma_electron_vol_avg_kev)",
            physics_variables.temp_plasma_electron_vol_avg_kev,
        )
        po.ovarrf(
            self.outfile,
            "Ratio of ion to electron volume-averaged temperature",
            "(f_temp_plasma_ion_electron)",
            physics_variables.f_temp_plasma_ion_electron,
            "IP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temperature on axis (keV)",
            "(temp_plasma_electron_on_axis_kev)",
            physics_variables.temp_plasma_electron_on_axis_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ion temperature (keV)",
            "(temp_plasma_ion_vol_avg_kev)",
            physics_variables.temp_plasma_ion_vol_avg_kev,
        )
        po.ovarrf(
            self.outfile,
            "Ion temperature on axis (keV)",
            "(temp_plasma_ion_on_axis_kev)",
            physics_variables.temp_plasma_ion_on_axis_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temp., density weighted (keV)",
            "(temp_plasma_electron_density_weighted_kev)",
            physics_variables.temp_plasma_electron_density_weighted_kev,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ratio of electron density weighted temp. to volume averaged temp.",
            "(f_temp_plasma_electron_density_vol_avg)",
            physics_variables.f_temp_plasma_electron_density_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged electron number density (/m3)",
            "(nd_plasma_electrons_vol_avg)",
            physics_variables.nd_plasma_electrons_vol_avg,
        )
        po.ovarre(
            self.outfile,
            "Electron number density on axis (/m3)",
            "(nd_plasma_electron_on_axis)",
            physics_variables.nd_plasma_electron_on_axis,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line-averaged electron number density (/m3)",
            "(nd_plasma_electron_line)",
            physics_variables.nd_plasma_electron_line,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma pressure on axis (Pa)",
            "(pres_plasma_thermal_on_axis)",
            physics_variables.pres_plasma_thermal_on_axis,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged plasma pressure (Pa)",
            "(pres_plasma_thermal_vol_avg)",
            physics_variables.pres_plasma_thermal_vol_avg,
            "OP ",
        )
        # As array output is not currently supported, each element is output as a float instance
        # Output plasma pressure profiles to mfile
        for i in range(len(physics_variables.pres_plasma_thermal_total_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma pressure at point {i}",
                f"(pres_plasma_thermal_total_profile{i})",
                physics_variables.pres_plasma_thermal_total_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma electron pressure at point {i}",
                f"(pres_plasma_electron_profile{i})",
                physics_variables.pres_plasma_electron_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_ion_total_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma ion pressure at point {i}",
                f"(pres_plasma_ion_total_profile{i})",
                physics_variables.pres_plasma_ion_total_profile[i],
            )
        for i in range(len(physics_variables.pres_plasma_fuel_profile)):
            po.ovarre(
                self.mfile,
                f"Total plasma fuel pressure at point {i}",
                f"(pres_plasma_fuel_profile{i})",
                physics_variables.pres_plasma_fuel_profile[i],
            )

        if stellarator_variables.istell == 0:
            po.ovarre(
                self.outfile,
                "Line-averaged electron density / Greenwald density",
                "(dnla_gw)",
                physics_variables.nd_plasma_electron_line
                / physics_variables.nd_plasma_electron_max_array[6],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total Ion number density (/m3)",
            "(nd_plasma_ions_total_vol_avg)",
            physics_variables.nd_plasma_ions_total_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel ion number density (/m3)",
            "(nd_plasma_fuel_ions_vol_avg)",
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total impurity number density with Z > 2 (no He) (/m3)",
            "(nd_plasma_impurities_vol_avg)",
            physics_variables.nd_plasma_impurities_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion number density (thermalised ions only) (/m3)",
            "(nd_plasma_alphas_vol_avg)",
            physics_variables.nd_plasma_alphas_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) / electron number density",
            "(f_nd_alpha_electron)",
            physics_variables.f_nd_alpha_electron,
        )
        po.ovarre(
            self.outfile,
            "Plasma volume averaged proton number density (/m3)",
            "(nd_plasma_protons_vol_avg)",
            physics_variables.nd_plasma_protons_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Proton number density / electron number density",
            "(f_nd_protium_electrons)",
            physics_variables.f_nd_protium_electrons,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Hot beam ion number density (/m3)",
            "(nd_beam_ions)",
            physics_variables.nd_beam_ions,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hot beam ion number density / electron density",
            "(f_nd_beam_electron)",
            physics_variables.f_nd_beam_electron,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Impurities:")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Plasma ion densities / electron density:")

        for imp in range(impurity_radiation_module.N_IMPURITIES):
            # MDK Update f_nd_impurity_electrons, as this will make the ITV output work correctly.
            impurity_radiation_module.f_nd_impurity_electrons[imp] = (
                impurity_radiation_module.f_nd_impurity_electron_array[imp]
            )
            str1 = (
                impurity_radiation_module.impurity_arr_label[imp].item()
                + " concentration"
            )
            str2 = f"(f_nd_impurity_electrons({imp + 1:02}))"
            # MDK Add output flag for H which is calculated.
            if imp == 0:
                po.ovarre(
                    self.outfile,
                    str1,
                    str2,
                    impurity_radiation_module.f_nd_impurity_electrons[imp],
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    str1,
                    str2,
                    impurity_radiation_module.f_nd_impurity_electrons[imp],
                )
                if impurity_radiation_module.f_nd_impurity_electrons[imp] != 0.0:  # noqa: RUF069
                    for i in range(physics_variables.n_plasma_profile_elements):
                        po.ovarre(
                            self.mfile,
                            str1 + f" at point {i}",
                            f"(f_nd_impurity_electrons{imp}_{i})",
                            impurity_radiation_module.f_nd_impurity_electrons[imp]
                            * self.plasma_profile.neprofile.profile_y[i],
                            "OP ",
                        )

        po.ovarre(
            self.outfile,
            "Average mass of all ions (amu)",
            "(m_ions_total_amu)",
            physics_variables.m_ions_total_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all ions in plasma (kg)",
            "(m_plasma_ions_total)",
            physics_variables.m_plasma_ions_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Average mass of all fuel ions (amu)",
            "(m_fuel_amu)",
            physics_variables.m_fuel_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all fuel ions in plasma (kg)",
            "(m_plasma_fuel_ions)",
            physics_variables.m_plasma_fuel_ions,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Average mass of all beam ions (amu)",
            "(m_beam_amu)",
            physics_variables.m_beam_amu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all alpha particles in plasma (kg)",
            "(m_plasma_alpha)",
            physics_variables.m_plasma_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of all electrons in plasma (kg)",
            "(m_plasma_electron)",
            physics_variables.m_plasma_electron,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total mass of the plasma (kg)",
            "(m_plasma)",
            physics_variables.m_plasma,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Volume averaged plasma effective charge",
            "(n_charge_plasma_effective_vol_avg)",
            physics_variables.n_charge_plasma_effective_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.n_charge_plasma_effective_profile)):
            po.ovarre(
                self.mfile,
                "Effective charge at point",
                f"(n_charge_plasma_effective_profile{i})",
                physics_variables.n_charge_plasma_effective_profile[i],
                "OP ",
            )

        for imp in range(impurity_radiation_module.N_IMPURITIES):
            for i in range(physics_variables.n_plasma_profile_elements):
                po.ovarre(
                    self.mfile,
                    "Impurity charge at point",
                    f"(n_charge_plasma_profile{imp}_{i})",
                    impurity_radiation_module.n_charge_impurity_profile[imp][i],
                    "OP ",
                )

        po.ovarrf(
            self.outfile,
            "Mass-weighted Effective charge",
            "(n_charge_plasma_effective_mass_weighted_vol_avg)",
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            "OP ",
        )

        po.ovarrf(
            self.outfile, "Density profile factor", "(alphan)", physics_variables.alphan
        )
        po.ovarin(
            self.outfile,
            "Plasma profile model",
            "(i_plasma_pedestal)",
            physics_variables.i_plasma_pedestal,
        )

        if (
            PlasmaProfileShapeType(physics_variables.i_plasma_pedestal)
            == PlasmaProfileShapeType.PEDESTAL_PROFILE
        ):
            if (
                physics_variables.nd_plasma_electron_on_axis
                < physics_variables.nd_plasma_pedestal_electron
            ):
                logger.error("Central density is less than pedestal density")

            po.ocmmnt(self.outfile, "Pedestal profiles are used.")
            po.ovarrf(
                self.outfile,
                "Density pedestal r/a location",
                "(radius_plasma_pedestal_density_norm)",
                physics_variables.radius_plasma_pedestal_density_norm,
            )
            if physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(nd_plasma_pedestal_electron)",
                    physics_variables.nd_plasma_pedestal_electron,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(nd_plasma_pedestal_electron)",
                    physics_variables.nd_plasma_pedestal_electron,
                )

            # must be assigned to their exisiting values#
            fgwped_out = (
                physics_variables.nd_plasma_pedestal_electron
                / physics_variables.nd_plasma_electron_max_array[6]
            )
            fgwsep_out = (
                physics_variables.nd_plasma_separatrix_electron
                / physics_variables.nd_plasma_electron_max_array[6]
            )
            if physics_variables.f_nd_plasma_pedestal_greenwald >= 0e0:
                physics_variables.f_nd_plasma_pedestal_greenwald = (
                    physics_variables.nd_plasma_pedestal_electron
                    / physics_variables.nd_plasma_electron_max_array[6]
                )
            if physics_variables.f_nd_plasma_separatrix_greenwald >= 0e0:
                physics_variables.f_nd_plasma_separatrix_greenwald = (
                    physics_variables.nd_plasma_separatrix_electron
                    / physics_variables.nd_plasma_electron_max_array[6]
                )

            po.ovarre(
                self.outfile,
                "Electron density at pedestal / nGW",
                "(fgwped_out)",
                fgwped_out,
            )
            po.ovarrf(
                self.outfile,
                "Temperature pedestal r/a location",
                "(radius_plasma_pedestal_temp_norm)",
                physics_variables.radius_plasma_pedestal_temp_norm,
            )

            po.ovarrf(
                self.outfile,
                "Electron temp. pedestal height (keV)",
                "(temp_plasma_pedestal_kev)",
                physics_variables.temp_plasma_pedestal_kev,
            )
            if 78 in numerics.icc:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(temp_plasma_separatrix_kev)",
                    physics_variables.temp_plasma_separatrix_kev,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(temp_plasma_separatrix_kev)",
                    physics_variables.temp_plasma_separatrix_kev,
                )

            po.ovarre(
                self.outfile,
                "Electron density at separatrix (/m3)",
                "(nd_plasma_separatrix_electron)",
                physics_variables.nd_plasma_separatrix_electron,
            )
            po.ovarre(
                self.outfile,
                "Electron density at separatrix / nGW",
                "(fgwsep_out)",
                fgwsep_out,
            )

        # Issue 558 - addition of constraint 76 to limit the value of nd_plasma_separatrix_electron, in proportion with the ballooning parameter and Greenwald density
        if 76 in numerics.icc:
            po.ovarre(
                self.outfile,
                "Critical ballooning parameter value",
                "(alpha_crit)",
                physics_variables.alpha_crit,
            )
            po.ovarre(
                self.outfile,
                "Critical electron density at separatrix (/m3)",
                "(nd_plasma_separatrix_electron_eich_max)",
                physics_variables.nd_plasma_separatrix_electron_eich_max,
            )

        po.ovarrf(
            self.outfile,
            "Temperature profile index",
            "(alphat)",
            physics_variables.alphat,
        )
        po.ovarrf(
            self.outfile,
            "Temperature profile index beta",
            "(tbeta)",
            physics_variables.tbeta,
        )
        po.ovarrf(
            self.outfile,
            "Pressure profile index",
            "(alphap)",
            physics_variables.alphap,
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        if stellarator_variables.istell == 0:
            self.density_limit.output_density_limit_information()

        po.oheadr(self.outfile, "Plasma Reactions :")

        po.osubhd(self.outfile, "Fuel Constituents :")
        po.ovarrf(
            self.outfile,
            "Deuterium fuel fraction",
            "(f_plasma_fuel_deuterium)",
            physics_variables.f_plasma_fuel_deuterium,
        )
        po.ovarrf(
            self.outfile,
            "Tritium fuel fraction",
            "(f_plasma_fuel_tritium)",
            physics_variables.f_plasma_fuel_tritium,
        )
        po.ovarrf(
            self.outfile,
            "3-Helium fuel fraction",
            "(f_plasma_fuel_helium3)",
            physics_variables.f_plasma_fuel_helium3,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion rates :")
        po.ovarre(
            self.outfile,
            "Fusion rate: total (reactions/sec)",
            "(fusrat_total)",
            physics_variables.fusrat_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion rate density: total (reactions/m3/sec)",
            "(fusden_total)",
            physics_variables.fusden_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion rate density: plasma (reactions/m3/sec)",
            "(fusden_plasma)",
            physics_variables.fusden_plasma,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Fusion Powers :")
        po.ocmmnt(
            self.outfile,
            "Fusion power totals from the main plasma and beam-plasma interactions (if present)\n",
        )

        po.ovarre(
            self.outfile,
            "Total fusion power (MW)",
            "(p_fusion_total_mw)",
            physics_variables.p_fusion_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-T fusion power: total (MW)",
            "(p_dt_total_mw)",
            physics_variables.p_dt_total_mw,
            "OP ",
        )
        for i in range(len(physics_variables.fusrat_plasma_dt_profile)):
            po.ovarre(
                self.mfile,
                f"DT fusion rate at point {i}",
                f"fusrat_plasma_dt_profile{i}",
                physics_variables.fusrat_plasma_dt_profile[i],
            )

        for i in range(len(physics_variables.fusrat_plasma_dd_triton_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> T fusion rate at point {i}",
                f"fusrat_plasma_dd_triton_profile{i}",
                physics_variables.fusrat_plasma_dd_triton_profile[i],
            )
        for i in range(len(physics_variables.fusrat_plasma_dd_helion_profile)):
            po.ovarre(
                self.mfile,
                f"D-D -> 3He fusion rate at point {i}",
                f"fusrat_plasma_dd_helion_profile{i}",
                physics_variables.fusrat_plasma_dd_helion_profile[i],
            )
        for i in range(len(physics_variables.fusrat_plasma_dhe3_profile)):
            po.ovarre(
                self.mfile,
                f"D-3He fusion rate at point {i}",
                f"fusrat_plasma_dhe3_profile{i}",
                physics_variables.fusrat_plasma_dhe3_profile[i],
            )
        po.ovarre(
            self.outfile,
            "D-T fusion power: plasma (MW)",
            "(p_plasma_dt_mw)",
            physics_variables.p_plasma_dt_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-T fusion power: beam (MW)",
            "(p_beam_dt_mw)",
            physics_variables.p_beam_dt_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-D fusion power (MW)",
            "(p_dd_total_mw)",
            physics_variables.p_dd_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-D branching ratio for tritium producing reactions",
            "(f_dd_branching_trit)",
            physics_variables.f_dd_branching_trit,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "D-He3 fusion power (MW)",
            "(p_dhe3_total_mw)",
            physics_variables.p_dhe3_total_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Alpha Powers :")
        po.ovarre(
            self.outfile,
            "Alpha rate density: total (particles/m3/sec)",
            "(fusden_alpha_total)",
            physics_variables.fusden_alpha_total,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha rate density: plasma (particles/m3/sec)",
            "(fusden_plasma_alpha)",
            physics_variables.fusden_plasma_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: total (MW)",
            "(p_alpha_total_mw)",
            physics_variables.p_alpha_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power density: total (MW/m^3)",
            "(pden_alpha_total_mw)",
            physics_variables.pden_alpha_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: plasma only (MW)",
            "(p_plasma_alpha_mw)",
            physics_variables.p_plasma_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power density: plasma (MW/m^3)",
            "(pden_plasma_alpha_mw)",
            physics_variables.pden_plasma_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: beam-plasma (MW)",
            "(p_beam_alpha_mw)",
            physics_variables.p_beam_alpha_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to electrons (MW/m3)",
            "(f_pden_alpha_electron_mw)",
            physics_variables.f_pden_alpha_electron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power per unit volume transferred to ions (MW/m3)",
            "(f_pden_alpha_ions_mw)",
            physics_variables.f_pden_alpha_ions_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Neutron Powers :")
        po.ovarre(
            self.outfile,
            "Neutron power: total (MW)",
            "(p_neutron_total_mw)",
            physics_variables.p_neutron_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power density: total (MW/m^3)",
            "(pden_neutron_total_mw)",
            physics_variables.pden_neutron_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: plasma only (MW)",
            "(p_plasma_neutron_mw)",
            physics_variables.p_plasma_neutron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power density: plasma (MW/m^3)",
            "(pden_plasma_neutron_mw)",
            physics_variables.pden_plasma_neutron_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power: beam-plasma (MW)",
            "(p_beam_neutron_mw)",
            physics_variables.p_beam_neutron_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")

        po.osubhd(self.outfile, "Charged Particle Powers :")

        po.ovarre(
            self.outfile,
            "Charged particle power (p, 3He, T) (excluding alphas) (MW)",
            "(p_non_alpha_charged_mw)",
            physics_variables.p_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Charged particle power density (p, 3He, T) (excluding alphas) (MW)",
            "(pden_non_alpha_charged_mw)",
            physics_variables.pden_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total charged particle power (including alphas) (MW)",
            "(p_charged_particle_mw)",
            physics_variables.p_charged_particle_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        po.osubhd(self.outfile, "Plasma radiation powers (excluding SOL):")
        po.ovarre(
            self.outfile,
            "Plasma total synchrotron radiation power (MW)",
            "(p_plasma_sync_mw)",
            physics_variables.p_plasma_sync_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma total synchrotron radiation power density (MW/m^3)",
            "(pden_plasma_sync_mw)",
            physics_variables.pden_plasma_sync_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Synchrotron wall reflectivity factor",
            "(f_sync_reflect)",
            physics_variables.f_sync_reflect,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma normalised minor radius defining 'core' region",
            "(radius_plasma_core_norm)",
            impurity_radiation_module.radius_plasma_core_norm,
        )
        po.ovarre(
            self.outfile,
            "Fractional scaling of radiation power along core profile",
            "(f_p_plasma_core_rad_reduction)",
            impurity_radiation_module.f_p_plasma_core_rad_reduction,
        )
        po.ovarre(
            self.outfile,
            "Plasma total radiation power from core region (MW)",
            "(p_plasma_inner_rad_mw)",
            physics_variables.p_plasma_inner_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma total radiation power from edge region (MW)",
            "(p_plasma_outer_rad_mw)",
            physics_variables.p_plasma_outer_rad_mw,
            "OP ",
        )

        if stellarator_variables.istell != 0:
            po.ovarre(
                self.outfile,
                "SOL radiation power as imposed by f_rad (MW)",
                "(psolradmw)",
                physics_variables.psolradmw,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Plasma total radiation power from inside last closed flux surface (MW)",
            "(p_plasma_rad_mw)",
            physics_variables.p_plasma_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Separatrix radiation fraction (MW)",
            "(f_p_plasma_separatrix_rad)",
            physics_variables.f_p_plasma_separatrix_rad,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Average neutron flux at plasma surface (MW/m^2)",
            "(pflux_plasma_surface_neutron_avg_mw)",
            physics_variables.pflux_plasma_surface_neutron_avg_mw,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power incident on the divertor targets (MW)",
                "(ptarmw)",
                physics_variables.ptarmw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power to the lower divertor",
                "(f_p_div_lower)",
                physics_variables.f_p_div_lower,
                "IP ",
            )
            po.ovarre(
                self.outfile,
                "Outboard side heat flux decay length (m)",
                "(lambdaio)",
                physics_variables.lambdaio,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Midplane seperation of the two magnetic closed flux surfaces (m)",
                    "(drsep)",
                    physics_variables.drsep,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Fraction of power on the inner targets",
                "(fio)",
                physics_variables.fio,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower inner target",
                "(fLI)",
                physics_variables.fli,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower outer target",
                "(fLO)",
                physics_variables.flo,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper inner target",
                    "(fUI)",
                    physics_variables.fui,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper outer target",
                    "(fUO)",
                    physics_variables.fuo,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Power incident on the lower inner target (MW)",
                "(pLImw)",
                physics_variables.plimw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Power incident on the lower outer target (MW)",
                "(pLOmw)",
                physics_variables.plomw,
                "OP ",
            )
            if divertor_variables.n_divertors == 2:
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper innner target (MW)",
                    "(pUImw)",
                    physics_variables.puimw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper outer target (MW)",
                    "(pUOmw)",
                    physics_variables.puomw,
                    "OP ",
                )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Ohmic heating power (MW)",
            "(p_plasma_ohmic_mw)",
            physics_variables.p_plasma_ohmic_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power deposited in plasma",
            "(f_p_alpha_plasma_deposited)",
            physics_variables.f_p_alpha_plasma_deposited,
            "IP",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to electrons",
            "(f_alpha_electron)",
            physics_variables.f_alpha_electron,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to ions",
            "(f_alpha_ion)",
            physics_variables.f_alpha_ion,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ion transport (MW)",
            "(p_ion_transport_loss_mw)",
            physics_variables.p_ion_transport_loss_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electron transport (MW)",
            "(p_electron_transport_loss_mw)",
            physics_variables.p_electron_transport_loss_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to ions (MW)",
            "(p_hcd_injected_ions_mw)",
            current_drive_variables.p_hcd_injected_ions_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to electrons (MW)",
            "(p_hcd_injected_electrons_mw)",
            current_drive_variables.p_hcd_injected_electrons_mw,
            "OP ",
        )
        if physics_variables.i_plasma_ignited == 1:
            po.ocmmnt(self.outfile, "  (Injected power only used for start-up phase)")

        po.ovarin(
            self.outfile,
            "Ignited plasma switch (0=not ignited, 1=ignited)",
            "(i_plasma_ignited)",
            physics_variables.i_plasma_ignited,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma separatrix power (MW)",
            "(p_plasma_separatrix_mw)",
            physics_variables.p_plasma_separatrix_mw,
            "OP ",
        )

        if physics_variables.p_plasma_separatrix_mw <= 0.001e0:
            logger.error(
                "Possible problem with high radiation power, forcing p_plasma_separatrix_mw to odd values. "
                f"{physics_variables.p_plasma_separatrix_mw=}"
            )
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(self.outfile, "          Power into divertor zone is unrealistic;")
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if divertor_variables.n_divertors == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Pdivt / R ratio (MW/m) (On peak divertor)",
                "(p_plasma_separatrix_rmajor_mw)",
                physics_variables.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)",
                "(p_div_bt_q_aspect_rmajor_mw)",
                physics_variables.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Psep / R ratio (MW/m)",
                "(p_plasma_separatrix_rmajor_mw)",
                physics_variables.p_plasma_separatrix_rmajor_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Psep Bt / qAR ratio (MWT/m)",
                "(p_div_bt_q_aspect_rmajor_mw)",
                physics_variables.p_div_bt_q_aspect_rmajor_mw,
                "OP ",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)

        if stellarator_variables.istell == 0:
            self.plasma_transition.output_l_h_threshold_powers()

        self.confinement.output_confinement_time_info()

        if stellarator_variables.istell == 0:
            # Issues 363 Output dimensionless plasma parameters MDK
            po.osubhd(self.outfile, "Dimensionless plasma parameters")
            po.ocmmnt(self.outfile, "For definitions see")
            po.ocmmnt(
                self.outfile,
                "Recent progress on the development and analysis of the ITPA global H-mode confinement database",
            )
            po.ocmmnt(
                self.outfile,
                "D.C. McDonald et al, 2007 Nuclear Fusion v47, 147. (nu_star missing 1/mu0)",
            )
            po.ovarre(
                self.outfile,
                "Normalized plasma pressure beta as defined by McDonald et al",
                "(beta_mcdonald)",
                physics_variables.beta_mcdonald,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized ion Larmor radius",
                "(rho_star)",
                physics_variables.rho_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized collisionality",
                "(nu_star)",
                physics_variables.nu_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER Physics Basis definition of elongation",
                "(kappa_ipb)",
                physics_variables.kappa_ipb,
                "OP ",
            )

            po.oblnkl(self.outfile)
            po.ostars(self.outfile, 110)
            po.oblnkl(self.outfile)

            self.inductance.output_volt_second_information()
        if stellarator_variables.istell == 0:
            self.plasma_bootstrap_current.output()
            self.dia_current.output()

        po.osubhd(self.outfile, "Fuelling :")
        po.ovarre(
            self.outfile,
            "Ratio of He and pellet particle confinement times",
            "(tauratio)",
            physics_variables.tauratio,
        )
        po.ovarre(
            self.outfile,
            "Fuelling rate (nucleus-pairs/s)",
            "(molflow_plasma_fuelling_required)",
            physics_variables.molflow_plasma_fuelling_required,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel burn-up rate (reactions/s)",
            "(rndfuel)",
            physics_variables.rndfuel,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Burn-up fraction",
            "(burnup)",
            physics_variables.burnup,
            "OP ",
        )

        if 78 in numerics.icc:
            po.osubhd(self.outfile, "Reinke Criterion :")
            po.ovarin(
                self.outfile,
                "index of impurity to be iterated for divertor detachment",
                "(impvardiv)",
                reinke_variables.impvardiv,
            )
            po.ovarre(
                self.outfile,
                "Minimum Impurity fraction from Reinke",
                "(fzmin)",
                reinke_variables.fzmin,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual Impurity fraction",
                "(fzactual)",
                reinke_variables.fzactual,
            )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_masses(
        m_fuel_amu: float,
        m_ions_total_amu: float,
        nd_plasma_ions_total_vol_avg: float,
        nd_plasma_fuel_ions_vol_avg: float,
        nd_plasma_alphas_vol_avg: float,
        vol_plasma: float,
        nd_plasma_electrons_vol_avg: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma masses.

        Parameters
        ----------
        m_fuel_amu : float
            Average mass of fuel (amu).
        m_ions_total_amu : float
            Average mass of all ions (amu).
        nd_plasma_ions_total_vol_avg : float
            Total ion density (/m3).
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (/m3).
        nd_plasma_alphas_vol_avg : float
            Alpha ash density (/m3).
        vol_plasma : float
            Plasma volume (m3).
        nd_plasma_electrons_vol_avg : float
            Volume averaged electron density (/m3).

        Returns
        -------
        tuple[float, float, float, float, float]
            A tuple containing:

        """
        # Calculate mass of fuel ions
        m_plasma_fuel_ions = (m_fuel_amu * constants.ATOMIC_MASS_UNIT) * (
            nd_plasma_fuel_ions_vol_avg * vol_plasma
        )

        m_plasma_ions_total = (m_ions_total_amu * constants.ATOMIC_MASS_UNIT) * (
            nd_plasma_ions_total_vol_avg * vol_plasma
        )

        m_plasma_alpha = (nd_plasma_alphas_vol_avg * vol_plasma) * constants.ALPHA_MASS

        m_plasma_electron = constants.ELECTRON_MASS * (
            nd_plasma_electrons_vol_avg * vol_plasma
        )

        m_plasma = m_plasma_electron + m_plasma_ions_total

        return (
            m_plasma_fuel_ions,
            m_plasma_ions_total,
            m_plasma_alpha,
            m_plasma_electron,
            m_plasma,
        )


def res_diff_time(rmajor, res_plasma, kappa95):
    """Calculates resistive diffusion time

    Parameters
    ----------
    rmajor :
        plasma major radius (m)
    res_plasma :
        plasma resistivity (Ohms)
    kappa95 :
        plasma elongation at 95% flux surface

    """
    return 2 * constants.RMU0 * rmajor / (res_plasma * kappa95)


def reinke_tsep(b_plasma_toroidal_on_axis, flh, qstar, rmajor, eps, fgw, kappa, lhat):
    """Function for calculating upstream temperature(keV) in Reinke model
    This function calculates the upstream temperature in the
    divertor/SoL model used for the Reinke citerion.
    Issue #707
    M.L. Reinke 2017 Nucl. Fusion 57 034004

    Parameters
    ----------
    ---------_
    b_plasma_toroidal_on_axis :
        toroidal field on axis (T)
    flh :
        fraction of Psep/P_LH
    qstar :
        safety factor similar to q95 (see #707)
    rmajor :
        major radius (m)
    eps :
        inverse aspect ratio
    fgw :
        ratio of volume averaged density to n_GW
    kappa :
        elongation
    lhat :
        connection length factor
    """
    kappa_0 = 2.0e3  # Stangeby W/m/eV^(7/2)

    return (
        (
            b_plasma_toroidal_on_axis**0.72
            * flh**0.29
            * fgw**0.21
            * qstar**0.08
            * rmajor**0.33
        )
        * (eps**0.15 * (1.0 + kappa**2.0) ** 0.34)
        * (lhat**0.29 * kappa_0 ** (-0.29) * 0.285)
    )


class BetaNormMaxModel(IntEnum):
    """Beta norm max (β_N_max) model types"""

    USER_INPUT = (0, "User Input")
    WESSON = (1, "Wesson Scaling")
    ORIGINAL_SCALING = (2, "Original Scaling")
    MENARD = (3, "Menard Scaling")
    THLOREUS = (4, "Thloreus Scaling")
    STAMBAUGH = (5, "Stambaugh Scaling")

    def __new__(cls, value, full_name):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        return self._full_name_


class BetaComponentLimits(IntEnum):
    """Beta component to apply limit types"""

    TOTAL = (0, "Total Beta")
    THERMAL = (1, "Thermal Beta")
    THERMAL_AND_BEAM = (2, "Thermal and Beam Beta")
    TOROIDAL = (3, "Toroidal Beta")

    def __new__(cls, value, full_name):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        return self._full_name_


class PlasmaBeta:
    """Class to hold plasma beta calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def get_beta_norm_max_value(self, model: BetaNormMaxModel) -> float:
        """Get the beta norm max value (β_N_max) for the specified model.

        Parameters
        ----------
        model: BetaNormMaxModel :
        """
        model_map = {
            BetaNormMaxModel.USER_INPUT: physics_variables.beta_norm_max,
            BetaNormMaxModel.WESSON: physics_variables.beta_norm_max_wesson,
            BetaNormMaxModel.ORIGINAL_SCALING: physics_variables.beta_norm_max_original_scaling,
            BetaNormMaxModel.MENARD: physics_variables.beta_norm_max_menard,
            BetaNormMaxModel.THLOREUS: physics_variables.beta_norm_max_thloreus,
            BetaNormMaxModel.STAMBAUGH: physics_variables.beta_norm_max_stambaugh,
        }
        return model_map[model]

    def run(self):
        """Calculate plasma beta values."""
        # -----------------------------------------------------
        # Normalised Beta Limit
        # -----------------------------------------------------

        # Normalised beta from Troyon beta limit
        physics_variables.beta_norm_total = self.calculate_normalised_beta(
            beta=physics_variables.beta_total_vol_avg,
            rminor=physics_variables.rminor,
            c_plasma=physics_variables.plasma_current,
            b_field=physics_variables.b_plasma_toroidal_on_axis,
        )

        # Define beta_norm_max calculations

        physics_variables.beta_norm_max_wesson = self.calculate_beta_norm_max_wesson(
            ind_plasma_internal_norm=physics_variables.ind_plasma_internal_norm
        )

        # Original scaling law
        physics_variables.beta_norm_max_original_scaling = (
            self.calculate_beta_norm_max_original(eps=physics_variables.eps)
        )

        # J. Menard scaling law
        physics_variables.beta_norm_max_menard = self.calculate_beta_norm_max_menard(
            eps=physics_variables.eps
        )

        # E. Tholerus scaling law
        physics_variables.beta_norm_max_thloreus = self.calculate_beta_norm_max_thloreus(
            c_beta=physics_variables.c_beta,
            pres_plasma_on_axis=physics_variables.pres_plasma_thermal_on_axis,
            pres_plasma_vol_avg=physics_variables.pres_plasma_thermal_vol_avg,
        )

        # R. D. Stambaugh scaling law
        physics_variables.beta_norm_max_stambaugh = (
            self.calculate_beta_norm_max_stambaugh(
                f_c_plasma_bootstrap=current_drive_variables.f_c_plasma_bootstrap,
                kappa=physics_variables.kappa,
                aspect=physics_variables.aspect,
            )
        )

        # Calculate beta_norm_max based on i_beta_norm_max
        try:
            model = BetaNormMaxModel(int(physics_variables.i_beta_norm_max))
            physics_variables.beta_norm_max = self.get_beta_norm_max_value(model)
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_beta_norm_max",
                i_beta_norm_max=physics_variables.i_beta_norm_max,
            ) from None

        # calculate_beta_limit() returns the beta_vol_avg_max for beta
        physics_variables.beta_vol_avg_max = self.calculate_beta_limit_from_norm(
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            beta_norm_max=physics_variables.beta_norm_max,
            plasma_current=physics_variables.plasma_current,
            rminor=physics_variables.rminor,
        )

        physics_variables.beta_toroidal_vol_avg = (
            physics_variables.beta_total_vol_avg
            * physics_variables.b_plasma_total**2
            / physics_variables.b_plasma_toroidal_on_axis**2
        )

        # Calculate physics_variables.beta poloidal [-]
        physics_variables.beta_poloidal_vol_avg = self.calculate_poloidal_beta(
            b_plasma_total=physics_variables.b_plasma_total,
            b_plasma_poloidal_average=physics_variables.b_plasma_surface_poloidal_average,
            beta=physics_variables.beta_total_vol_avg,
        )

        physics_variables.beta_thermal_vol_avg = (
            physics_variables.beta_total_vol_avg
            - physics_variables.beta_fast_alpha
            - physics_variables.beta_beam
        )

        physics_variables.beta_poloidal_eps = (
            physics_variables.beta_poloidal_vol_avg * physics_variables.eps
        )

        physics_variables.beta_thermal_poloidal_vol_avg = (
            physics_variables.beta_thermal_vol_avg
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_surface_poloidal_average
            )
            ** 2
        )
        physics_variables.beta_thermal_toroidal_vol_avg = (
            physics_variables.beta_thermal_vol_avg
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )

        # =======================================================

        # Mirror the pressure profiles to match the doubled toroidal field profile
        pres_profile_total = np.concatenate([
            physics_variables.pres_plasma_thermal_total_profile[::-1],
            physics_variables.pres_plasma_thermal_total_profile,
        ])

        physics_variables.beta_thermal_toroidal_profile = np.array([
            self.calculate_plasma_beta(
                pres_plasma=pres_profile_total[i],
                b_field=physics_variables.b_plasma_toroidal_profile[i],
            )
            for i in range(len(physics_variables.b_plasma_toroidal_profile))
        ])

        # =======================================================

        physics_variables.beta_norm_thermal = self.calculate_normalised_beta(
            beta=physics_variables.beta_thermal_vol_avg,
            rminor=physics_variables.rminor,
            c_plasma=physics_variables.plasma_current,
            b_field=physics_variables.b_plasma_toroidal_on_axis,
        )

        physics_variables.beta_norm_toroidal = (
            physics_variables.beta_norm_total
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_toroidal_on_axis
            )
            ** 2
        )
        physics_variables.beta_norm_poloidal = (
            physics_variables.beta_norm_total
            * (
                physics_variables.b_plasma_total
                / physics_variables.b_plasma_surface_poloidal_average
            )
            ** 2
        )

        physics_variables.f_beta_alpha_beam_thermal = (
            physics_variables.beta_fast_alpha + physics_variables.beta_beam
        ) / physics_variables.beta_thermal_vol_avg

        # Plasma thermal energy derived from the thermal beta
        physics_variables.e_plasma_beta_thermal = self.calculate_plasma_energy_from_beta(
            beta=physics_variables.beta_thermal_vol_avg,
            b_field=physics_variables.b_plasma_total,
            vol_plasma=physics_variables.vol_plasma,
        )

        # Plasma thermal energy derived from the total beta
        physics_variables.e_plasma_beta = self.calculate_plasma_energy_from_beta(
            beta=physics_variables.beta_total_vol_avg,
            b_field=physics_variables.b_plasma_total,
            vol_plasma=physics_variables.vol_plasma,
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_beta(
        pres_plasma: float | np.ndarray, b_field: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the plasma beta (β) for a given pressure and field.

        Plasma beta is the ratio of plasma pressure to magnetic pressure.

        Parameters
        ----------
        pres_plasma : float | np.ndarray
            Plasma pressure (in Pascals).
        b_field : float | np.ndarray
            Magnetic field strength (in Tesla).

        Returns
        -------
        float | np.ndarray
            The plasma beta (dimensionless).
        """
        return 2 * constants.RMU0 * pres_plasma / (b_field**2)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_norm_max_wesson(ind_plasma_internal_norm: float) -> float:
        """Calculate the Wesson normalsied beta upper limit.

        Parameters
        ----------
        ind_plasma_internal_norm : float
            Plasma normalised internal inductance


        Returns
        -------
        float
            The Wesson normalised beta upper limit.

        Notes
        -----
             - It is recommended to use this method with the other Wesson relations for normalsied internal
             inductance and current profile index.
             - This fit is derived from the DIII-D database for β_N >= 2.5

        References
        ----------
             - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
             International Series of Monographs on Physics, Volume 149.

             - T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,”
             Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).

        """
        return 4 * ind_plasma_internal_norm

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_norm_max_original(eps: float) -> float:
        """Calculate the original scaling law normalsied beta upper limit.

        Parameters
        ----------
        eps : float
            Plasma normalised internal inductance

        Returns
        -------
        float

        References
        ----------
            The original scaling law normalised beta upper limit.

        """
        return 2.7 * (1.0 + 5.0 * eps**3.5)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_norm_max_menard(eps: float) -> float:
        """Calculate the Menard normalsied beta upper limit.

        Parameters
        ----------
        eps : float
            Plasma normalised internal inductance

        Returns
        -------
        float
            The Menard normalised beta upper limit.

        Notes
        -----
            - Found as a reasonable fit to the computed no wall limit at f_BS ≈ 50%
            - Uses maximum κ data from NSTX at A = 1.45, A = 1.75. Along with record
              β_T data from DIII-D at A = 2.9 and high κ.

        References
        ----------
            - # J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.

        """
        return 3.12 + 3.5 * eps**1.7

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_norm_max_thloreus(
        c_beta: float, pres_plasma_on_axis: float, pres_plasma_vol_avg: float
    ) -> float:
        """Calculate the E. Tholerus normalized beta upper limit.

        Parameters
        ----------
        c_beta : float
            Pressure peaking factor coefficient.
        pres_plasma_on_axis : float
            Central plasma pressure (Pa).
        pres_plasma_vol_avg : float
            Volume-averaged plasma pressure (Pa).
        c_beta: float :

        pres_plasma_on_axis: float :

        pres_plasma_vol_avg: float :


        Returns
        -------
        float
            The E. Tholerus normalized beta upper limit.

        Notes
        -----
            - This method calculates the normalized beta upper limit based on the pressure peaking factor (Fp),
              which is defined as the ratio of the peak pressure to the average pressure.
            - The formula is derived from operational space studies of flat-top plasma in the STEP power plant.

        References
        ----------
            - E. Tholerus et al., “Flat-top plasma operational space of the STEP power plant,”
              Nuclear Fusion, Aug. 2024, doi: https://doi.org/10.1088/1741-4326/ad6ea2.

        """
        return 3.7 + (
            (c_beta / (pres_plasma_on_axis / pres_plasma_vol_avg))
            * (12.5 - 3.5 * (pres_plasma_on_axis / pres_plasma_vol_avg))
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_norm_max_stambaugh(
        f_c_plasma_bootstrap: float,
        kappa: float,
        aspect: float,
    ) -> float:
        """Calculate the Stambaugh normalized beta upper limit.

        Parameters
        ----------
        f_c_plasma_bootstrap : float
            Bootstrap current fraction.
        kappa : float
            Plasma separatrix elongation.
        aspect : float
            Plasma aspect ratio.


        Returns
        -------
        float
            The Stambaugh normalized beta upper limit.

        Notes
        -----
            - This method calculates the normalized beta upper limit based on the Stambaugh scaling.
            - The formula is derived from empirical fits to high-performance, steady-state tokamak equilibria.

        References
        ----------
            - R. D. Stambaugh et al., “Fusion Nuclear Science Facility Candidates,”
              Fusion Science and Technology, vol. 59, no. 2, pp. 279-307, Feb. 2011,
              doi: https://doi.org/10.13182/fst59-279.

            - Y. R. Lin-Liu and R. D. Stambaugh, “Optimum equilibria for high performance, steady state tokamaks,”
              Nuclear Fusion, vol. 44, no. 4, pp. 548-554, Mar. 2004,
              doi: https://doi.org/10.1088/0029-5515/44/4/009.

        """
        return (
            f_c_plasma_bootstrap
            * 10
            * (-0.7748 + (1.2869 * kappa) - (0.2921 * kappa**2) + (0.0197 * kappa**3))
            / (aspect**0.5523 * np.tanh((1.8524 + (0.2319 * kappa)) / aspect**0.6163))
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_normalised_beta(
        beta: float, rminor: float, c_plasma: float, b_field: float
    ) -> float:
        """Calculate normalised beta (β_N).

        Parameters
        ----------
        beta : float
            Plasma beta (fraction).
        rminor : float
            Plasma minor radius (m).
        c_plasma : float
            Plasma current (A).
        b_field : float
            Magnetic field (T).

        Returns
        -------
        float
            Normalised beta.

        Notes
        -----
        - 1.0e8 is a conversion factor to get beta_N in standard units, as plasma current is normally in MA and
        beta is in percentage instead of fraction.

        """
        return 1.0e8 * (beta * rminor * b_field) / c_plasma

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_energy_from_beta(
        beta: float, b_field: float, vol_plasma: float
    ) -> float:
        """Calculate plasma thermal energy from beta.

        E_plasma = 1.5 * β * B² / (2 * μ_0) * V

        Parameters
        ----------
        beta : float
            Plasma beta (fraction).
        b_field : float
            Magnetic field (T).
        vol_plasma : float
            Plasma volume (m³).

        Returns
        -------
        float
            Plasma energy (J).
        """
        return (1.5e0 * beta * b_field**2) / (2.0e0 * constants.RMU0) * vol_plasma

    @staticmethod
    @nb.njit(cache=True)
    def calculate_beta_limit_from_norm(
        b_plasma_toroidal_on_axis: float,
        beta_norm_max: float,
        plasma_current: float,
        rminor: float,
    ) -> float:
        """Calculate the maximum allowed beta (β) from a given normalised (β_N).

        This subroutine calculates the beta limit using the algorithm documented in AEA FUS 172.
        The limit applies to beta defined with respect to the total B-field.
        Switch i_beta_component determines which components of beta to include.

        Parameters
        ----------
        b_plasma_toroidal_on_axis : float
            Toroidal B-field on plasma axis [T].
        beta_norm_max : float
            Troyon-like g coefficient.
        plasma_current : float
            Plasma current [A].
        rminor : float
            Plasma minor axis [m].

        Returns
        -------
        float
            Beta limit as defined below.

        Notes
        -----
            - If i_beta_component = 0, then the limit is applied to the total beta.
            - If i_beta_component = 1, then the limit is applied to the thermal beta only.
            - If i_beta_component = 2, then the limit is applied to the thermal + neutral beam beta components.
            - If i_beta_component = 3, then the limit is applied to the toroidal beta.

            - The default value for the g coefficient is beta_norm_max = 3.5.

        References
        ----------
            - F. Troyon et.al,  “Beta limit in tokamaks. Experimental and computational status,”
            Plasma Physics and Controlled Fusion, vol. 30, no. 11, pp. 1597-1609, Oct. 1988,
            doi: https://doi.org/10.1088/0741-3335/30/11/019.

            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        # Multiplied by 0.01 to convert from % to fraction
        return (
            0.01
            * beta_norm_max
            * (plasma_current / 1.0e6)
            / (rminor * b_plasma_toroidal_on_axis)
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_poloidal_beta(
        b_plasma_total: float, b_plasma_poloidal_average: float, beta: float
    ) -> float:
        """Calculates total poloidal beta (β_p)

        Parameters
        ----------
        b_plasma_poloidal_average : float
            The average poloidal magnetic field of the plasma (in Tesla).
        beta : float
            The plasma beta, a dimensionless parameter representing the ratio of plasma pressure to magnetic pressure.

        Returns
        -------
        float
            The calculated total poloidal beta.

        References
        ----------
        - J.P. Freidberg, "Plasma physics and fusion energy", Cambridge University Press (2007)
        Page 270 ISBN 0521851076

        """
        return beta * (b_plasma_total / b_plasma_poloidal_average) ** 2

    @staticmethod
    @nb.njit(cache=True)
    def fast_alpha_beta(
        b_plasma_poloidal_average: float,
        b_plasma_toroidal_on_axis: float,
        nd_plasma_electrons_vol_avg: float,
        nd_plasma_fuel_ions_vol_avg: float,
        nd_plasma_ions_total_vol_avg: float,
        temp_plasma_electron_density_weighted_kev: float,
        temp_plasma_ion_density_weighted_kev: float,
        pden_alpha_total_mw: float,
        pden_plasma_alpha_mw: float,
        i_beta_fast_alpha: int,
    ) -> float:
        """Calculate the fast alpha beta (β_fast_alpha) component.

        This function computes the fast alpha beta contribution based on the provided plasma parameters.

        Parameters
        ----------
        b_plasma_poloidal_average : float
            Poloidal field (T).
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        nd_plasma_electrons_vol_avg : float
            Electron density (m⁻³).
        nd_plasma_fuel_ions_vol_avg : float
            Fuel ion density (m⁻³).
        nd_plasma_ions_total_vol_avg : float
            Total ion density (m⁻³).
        temp_plasma_electron_density_weighted_kev : float
            Density-weighted electron temperature (keV).
        temp_plasma_ion_density_weighted_kev : float
            Density-weighted ion temperature (keV).
        pden_alpha_total_mw : float
            Alpha power per unit volume, from beams and plasma (MW/m³).
        pden_plasma_alpha_mw : float
            Alpha power per unit volume just from plasma (MW/m³).
        i_beta_fast_alpha : int
            Switch for fast alpha pressure method.

        Returns
        -------
        float
            Fast alpha beta component.

        Notes
        -----
            - For IPDG89 scaling applicability is Z_eff = 1.5, T_i/T_e = 1, 〈T〉 = 5-20 keV


        References
        ----------
            - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
            https://inis.iaea.org/collection/NCLCollectionStore/_Public/21/068/21068960.pdf

            - Uckan, N. A., Tolliver, J. S., Houlberg, W. A., and Attenberger, S. E.
            Influence of fast alpha diffusion and thermal alpha buildup on tokamak reactor performance.
            United States: N. p., 1987. Web.https://www.osti.gov/servlets/purl/5611706

        """
        # Determine average fast alpha density
        if physics_variables.f_plasma_fuel_deuterium < 1.0:
            beta_thermal = (
                2.0
                * constants.RMU0
                * constants.KILOELECTRON_VOLT
                * (
                    nd_plasma_electrons_vol_avg
                    * temp_plasma_electron_density_weighted_kev
                    + nd_plasma_ions_total_vol_avg * temp_plasma_ion_density_weighted_kev
                )
                / (b_plasma_toroidal_on_axis**2 + b_plasma_poloidal_average**2)
            )

            # jlion: This "fact" model is heavily flawed for smaller temperatures! It is unphysical for a stellarator (high n low T)
            # IPDG89 fast alpha scaling
            if i_beta_fast_alpha == 0:
                fact = min(
                    0.3,
                    0.29
                    * (nd_plasma_fuel_ions_vol_avg / nd_plasma_electrons_vol_avg) ** 2
                    * (
                        (
                            temp_plasma_electron_density_weighted_kev
                            + temp_plasma_ion_density_weighted_kev
                        )
                        / 20.0
                        - 0.37
                    ),
                )

            # Modified scaling, D J Ward
            else:
                fact = min(
                    0.30,
                    0.26
                    * (nd_plasma_fuel_ions_vol_avg / nd_plasma_electrons_vol_avg) ** 2
                    * np.sqrt(
                        max(
                            0.0,
                            (
                                (
                                    temp_plasma_electron_density_weighted_kev
                                    + temp_plasma_ion_density_weighted_kev
                                )
                                / 20.0
                                - 0.65
                            ),
                        )
                    ),
                )

            fact = max(fact, 0.0)
            fact2 = pden_alpha_total_mw / pden_plasma_alpha_mw
            beta_fast_alpha = beta_thermal * fact * fact2

        else:  # negligible alpha production, alpha_power_density = p_beam_alpha_mw = 0
            beta_fast_alpha = 0.0

        return beta_fast_alpha

    def output_beta_information(self):
        """Output beta information to file."""
        po.oheadr(self.outfile, "Plasma Beta:")

        po.ovarin(
            self.outfile,
            "Beta component for limits",
            "(i_beta_component)",
            physics_variables.i_beta_component,
        )

        if physics_variables.i_beta_component == BetaComponentLimits.TOTAL:
            po.ovarrf(
                self.outfile,
                "Upper limit on total beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )
        elif physics_variables.i_beta_component == BetaComponentLimits.THERMAL:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )
        elif physics_variables.i_beta_component == BetaComponentLimits.THERMAL_AND_BEAM:
            po.ovarrf(
                self.outfile,
                "Upper limit on thermal + NB beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )
        elif physics_variables.i_beta_component == BetaComponentLimits.TOROIDAL:
            po.ovarrf(
                self.outfile,
                "Upper limit on toroidal beta",
                "(beta_vol_avg_max)",
                physics_variables.beta_vol_avg_max,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total plasma beta",
            "(beta_total_vol_avg)",
            physics_variables.beta_total_vol_avg,
        )
        if physics_variables.i_beta_component == BetaComponentLimits.TOTAL:
            po.ovarrf(
                self.outfile,
                "Lower limit on total beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        elif physics_variables.i_beta_component == BetaComponentLimits.THERMAL:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Lower limit on thermal + NB beta",
                "(beta_vol_avg_min)",
                physics_variables.beta_vol_avg_min,
                "IP",
            )
        po.ovarre(
            self.outfile,
            "Upper limit on poloidal beta",
            "(beta_poloidal_max)",
            constraint_variables.beta_poloidal_max,
            "IP",
        )
        po.ovarre(
            self.outfile,
            "Total poloidal beta",
            "(beta_poloidal_vol_avg)",
            physics_variables.beta_poloidal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume averaged toroidal beta",
            "(beta_toroidal_vol_avg)",
            physics_variables.beta_toroidal_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.beta_thermal_toroidal_profile)):
            po.ovarre(
                self.mfile,
                f"Beta toroidal profile at point {i}",
                f"beta_thermal_toroidal_profile{i}",
                physics_variables.beta_thermal_toroidal_profile[i],
            )

        po.ovarre(
            self.outfile,
            "Fast alpha beta",
            "(beta_fast_alpha)",
            physics_variables.beta_fast_alpha,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutral Beam ion beta",
            "(beta_beam)",
            physics_variables.beta_beam,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ratio of fast alpha and beam beta to thermal beta",
            "(f_beta_alpha_beam_thermal)",
            physics_variables.f_beta_alpha_beam_thermal,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Volume averaged thermal beta",
            "(beta_thermal_vol_avg)",
            physics_variables.beta_thermal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal poloidal beta",
            "(beta_thermal_poloidal_vol_avg)",
            physics_variables.beta_thermal_poloidal_vol_avg,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal toroidal beta",
            "(beta_thermal_toroidal_vol_avg)",
            physics_variables.beta_thermal_toroidal_vol_avg,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Poloidal beta and inverse aspect ratio",
            "(beta_poloidal_eps)",
            physics_variables.beta_poloidal_eps,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Poloidal beta and inverse aspect ratio upper limit",
            "(beta_poloidal_eps_max)",
            physics_variables.beta_poloidal_eps_max,
        )
        po.osubhd(self.outfile, "Normalised Beta Information :")

        po.ovarin(
            self.outfile,
            "Maximum normalised beta model",
            "(i_beta_norm_max)",
            physics_variables.i_beta_norm_max,
        )

        if stellarator_variables.istell == 0:
            if (
                BetaNormMaxModel(physics_variables.i_beta_norm_max)
                != BetaNormMaxModel.USER_INPUT
            ):
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(beta_norm_max)",
                    physics_variables.beta_norm_max,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(beta_norm_max)",
                    physics_variables.beta_norm_max,
                )
            po.ovarrf(
                self.outfile,
                "Normalised total beta",
                "(beta_norm_total)",
                physics_variables.beta_norm_total,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Normalised thermal beta",
                "(beta_norm_thermal) ",
                physics_variables.beta_norm_thermal,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised toroidal beta",
                "(beta_norm_toroidal) ",
                physics_variables.beta_norm_toroidal,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised poloidal beta",
                "(beta_norm_poloidal) ",
                physics_variables.beta_norm_poloidal,
                "OP ",
            )

            po.osubhd(self.outfile, "Maximum normalised beta scalings :")
            po.ovarrf(
                self.outfile,
                "J. Wesson normalised beta upper limit",
                "(beta_norm_max_wesson) ",
                physics_variables.beta_norm_max_wesson,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Original normalsied beta upper limit",
                "(beta_norm_max_original_scaling) ",
                physics_variables.beta_norm_max_original_scaling,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "J. Menard normalised beta upper limit",
                "(beta_norm_max_menard) ",
                physics_variables.beta_norm_max_menard,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "E. Thloreus normalised beta upper limit",
                "(beta_norm_max_thloreus) ",
                physics_variables.beta_norm_max_thloreus,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "R. Stambaugh normalised beta upper limit",
                "(beta_norm_max_stambaugh) ",
                physics_variables.beta_norm_max_stambaugh,
                "OP ",
            )

        po.osubhd(self.outfile, "Plasma energies derived from beta :")
        po.ovarre(
            self.outfile,
            "Plasma thermal energy derived from thermal beta (J)",
            "(e_plasma_beta_thermal) ",
            physics_variables.e_plasma_beta_thermal,
            "OP",
        )

        po.ovarre(
            self.outfile,
            "Plasma thermal energy derived from the total beta (J)",
            "(e_plasma_beta)",
            physics_variables.e_plasma_beta,
            "OP",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)


class IndInternalNormModel(IntEnum):
    """Normalised internal inductance (l_i) model types"""

    USER_INPUT = 0
    WESSON = 1
    MENARD = 2


class PlasmaInductance:
    """Class to hold plasma inductance calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        physics_variables.ind_plasma_internal_norm_wesson = (
            self.calculate_internal_inductance_wesson(alphaj=physics_variables.alphaj)
        )

        # Spherical Tokamak relation for internal inductance
        # Menard et al. (2016), Nuclear Fusion, 56, 106023
        physics_variables.ind_plasma_internal_norm_menard = (
            self.calculate_internal_inductance_menard(kappa=physics_variables.kappa)
        )

        physics_variables.ind_plasma_internal_norm_iter_3 = self.calculate_normalised_internal_inductance_iter_3(
            b_plasma_poloidal_vol_avg=physics_variables.b_plasma_surface_poloidal_average,
            c_plasma=physics_variables.plasma_current,
            vol_plasma=physics_variables.vol_plasma,
            rmajor=physics_variables.rmajor,
        )

        # Calculate ind_plasma_internal_norm based on i_ind_plasma_internal_norm
        try:
            model = IndInternalNormModel(
                int(physics_variables.i_ind_plasma_internal_norm)
            )
            physics_variables.ind_plasma_internal_norm = (
                self.get_ind_internal_norm_value(model)
            )
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_ind_plasma_internal_norm",
                i_ind_plasma_internal_norm=physics_variables.i_ind_plasma_internal_norm,
            ) from None

    def get_ind_internal_norm_value(self, model: IndInternalNormModel) -> float:
        """Get the normalised internal inductance (l_i) for the specified model.

        Parameters
        ----------
        model: IndInternalNormModel :
        """
        model_map = {
            IndInternalNormModel.USER_INPUT: physics_variables.ind_plasma_internal_norm,
            IndInternalNormModel.WESSON: physics_variables.ind_plasma_internal_norm_wesson,
            IndInternalNormModel.MENARD: physics_variables.ind_plasma_internal_norm_menard,
        }
        return model_map[model]

    @staticmethod
    @nb.njit(cache=True)
    def calculate_volt_second_requirements(
        csawth: float,
        eps: float,
        f_c_plasma_inductive: float,
        ejima_coeff: float,
        kappa: float,
        rmajor: float,
        res_plasma: float,
        plasma_current: float,
        t_plant_pulse_fusion_ramp: float,
        t_plant_pulse_burn: float,
        ind_plasma_internal_norm: float,
    ) -> tuple[float, float, float, float, float, float]:
        """Calculate the volt-second requirements and related parameters for plasma physics.

        Parameters
        ----------
        csawth : float
            Coefficient for sawteeth effects
        eps : float
            Inverse aspect ratio
        f_c_plasma_inductive : float
            Fraction of plasma current produced inductively
        ejima_coeff : float
            Ejima coefficient for resistive start-up V-s component
        kappa : float
            Plasma elongation
        rmajor : float
            Plasma major radius (m)
        res_plasma : float
            Plasma resistance (ohm)
        plasma_current : float
            Plasma current (A)
        t_plant_pulse_fusion_ramp : float
            Heating time (s)
        t_plant_pulse_burn : float
            Burn time (s)
        ind_plasma_internal_norm : float
            Plasma normalized internal inductance


        Returns
        -------
        tuple[float, float, float, float, float, float]
            A tuple containing:
            - vs_plasma_internal: Internal plasma volt-seconds (Wb)
            - ind_plasma_internal: Plasma inductance (H)
            - vs_plasma_burn_required: Volt-seconds needed during flat-top (heat+burn) (Wb)
            - vs_plasma_ramp_required: Volt-seconds needed during ramp-up (Wb)
            - ind_plasma_total,: Internal and external plasma inductance V-s (Wb)
            - vs_res_ramp: Resistive losses in start-up volt-seconds (Wb)
            - vs_plasma_total_required: Total volt-seconds needed (Wb)

        References
        ----------
            - S. Ejima, R. W. Callis, J. L. Luxon, R. D. Stambaugh, T. S. Taylor, and J. C. Wesley,
            “Volt-second analysis and consumption in Doublet III plasmas,”
            Nuclear Fusion, vol. 22, no. 10, pp. 1313-1319, Oct. 1982, doi:
            https://doi.org/10.1088/0029-5515/22/10/006.

            - S. C. Jardin, C. E. Kessel, and N Pomphrey,
            “Poloidal flux linkage requirements for the International Thermonuclear Experimental Reactor,”
            Nuclear Fusion, vol. 34, no. 8, pp. 1145-1160, Aug. 1994,
            doi: https://doi.org/10.1088/0029-5515/34/8/i07.

            - S. P. Hirshman and G. H. Neilson, “External inductance of an axisymmetric plasma,”
            The Physics of Fluids, vol. 29, no. 3, pp. 790-793, Mar. 1986,
            doi: https://doi.org/10.1063/1.865934.

        """
        # Plasma internal inductance

        ind_plasma_internal = constants.RMU0 * rmajor * ind_plasma_internal_norm / 2.0

        # Internal plasma flux (V-s) component
        vs_plasma_internal = ind_plasma_internal * plasma_current

        # Start-up resistive component
        # Uses ITER formula without the 10 V-s add-on

        vs_res_ramp = ejima_coeff * constants.RMU0 * plasma_current * rmajor

        # ======================================================================

        # Hirshman and Neilson fit for external inductance

        aeps = (1.0 + 1.81 * np.sqrt(eps) + 2.05 * eps) * np.log(8.0 / eps) - (
            2.0 + 9.25 * np.sqrt(eps) - 1.21 * eps
        )
        beps = 0.73 * np.sqrt(eps) * (1.0 + 2.0 * eps**4 - 6.0 * eps**5 + 3.7 * eps**6)

        ind_plasma_external = (
            rmajor * constants.RMU0 * aeps * (1.0 - eps) / (1.0 - eps + beps * kappa)
        )

        # ======================================================================

        ind_plasma_total = ind_plasma_external + ind_plasma_internal

        # Inductive V-s component

        vs_self_ind_ramp = ind_plasma_total * plasma_current
        vs_plasma_ramp_required = vs_res_ramp + vs_self_ind_ramp

        # Plasma loop voltage during flat-top
        # Include enhancement factor in flattop V-s requirement
        # to account for MHD sawtooth effects.

        v_plasma_loop_burn = plasma_current * res_plasma * f_c_plasma_inductive

        v_burn_resistive = v_plasma_loop_burn * csawth

        # N.B. t_plant_pulse_burn on first iteration will not be correct
        # if the pulsed reactor option is used, but the value
        # will be correct on subsequent calls.

        vs_plasma_burn_required = v_burn_resistive * (
            t_plant_pulse_fusion_ramp + t_plant_pulse_burn
        )
        vs_plasma_total_required = vs_plasma_ramp_required + vs_plasma_burn_required

        return (
            vs_plasma_internal,
            ind_plasma_total,
            vs_plasma_burn_required,
            vs_plasma_ramp_required,
            vs_self_ind_ramp,
            vs_res_ramp,
            vs_plasma_total_required,
            v_plasma_loop_burn,
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_normalised_internal_inductance_iter_3(
        b_plasma_poloidal_vol_avg: float,
        c_plasma: float,
        vol_plasma: float,
        rmajor: float,
    ) -> float:
        """Calculate the normalised internal inductance using ITER-3 scaling li(3).

        Parameters
        ----------
        b_plasma_poloidal_vol_avg : float
            Volume-averaged poloidal magnetic field (T).
        c_plasma : float
            Plasma current (A).
        vol_plasma : float
            Plasma volume (m^3).
        rmajor : float
            Plasma major radius (m).

        Returns
        -------
        float
            The li(3) normalised internal inductance.

        References
        ----------
            - T. C. Luce, D. A. Humphreys, G. L. Jackson, and W. M. Solomon,
            “Inductive flux usage and its optimization in tokamak operation,”
            Nuclear Fusion, vol. 54, no. 9, p. 093005, Jul. 2014,
            doi: https://doi.org/10.1088/0029-5515/54/9/093005.

            - G. L. Jackson et al., “ITER startup studies in the DIII-D tokamak,”
            Nuclear Fusion, vol. 48, no. 12, p. 125002, Nov. 2008,
            doi: https://doi.org/10.1088/0029-5515/48/12/125002.

        """
        return (
            2
            * vol_plasma
            * b_plasma_poloidal_vol_avg**2
            / (constants.RMU0**2 * c_plasma**2 * rmajor)
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_internal_inductance_menard(kappa: float) -> float:
        """Calculate the Menard plasma normalized internal inductance.

        Parameters
        ----------
        kappa : float
            Plasma separatrix elongation.

        Returns
        -------
        float
            The Menard plasma normalised internal inductance.

        Notes
        -----
            - This relation is based off of data from NSTX for l_i in the range of 0.4-0.85
            - This model is only recommended to be used for ST's with kappa > 2.5

        References
        ----------
            - J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,”
            Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016,
            doi: https://doi.org/10.1088/0029-5515/56/10/106023.

        """
        return 3.4 - kappa

    @staticmethod
    @nb.njit(cache=True)
    def calculate_internal_inductance_wesson(alphaj: float) -> float:
        """Calculate the Wesson plasma normalized internal inductance.

        Parameters
        ----------
        alphaj : float
            Current profile index.

        Returns
        -------
        float
            The Wesson plasma normalised internal inductance.

        Notes
        -----
            - It is recommended to use this method with the other Wesson relations for normalised beta and
              current profile index.
            - This relation is only true for the cyclindrical plasma approximation with parabolic profiles.

        References
        ----------
            - Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,
            International Series of Monographs on Physics, Volume 149.

        """
        return np.log(1.65 + 0.89 * alphaj)

    def output_volt_second_information(self):
        """Output volt-second information to file."""
        po.osubhd(self.outfile, "Plasma Volt-second Requirements :")
        po.ovarre(
            self.outfile,
            "Total plasma volt-seconds required for pulse (Wb)",
            "(vs_plasma_total_required)",
            physics_variables.vs_plasma_total_required,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plasma inductive flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_ind_ramp)",
            physics_variables.vs_plasma_ind_ramp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma resistive flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_res_ramp)",
            physics_variables.vs_plasma_res_ramp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total flux consumption for plasma current ramp-up (Wb)",
            "(vs_plasma_ramp_required)",
            physics_variables.vs_plasma_ramp_required,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ejima coefficient",
            "(ejima_coeff)",
            physics_variables.ejima_coeff,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Internal plasma V-s",
            "(vs_plasma_internal)",
            physics_variables.vs_plasma_internal,
        )

        po.ovarre(
            self.outfile,
            "Plasma volt-seconds needed for flat-top (heat + burn times) (Wb)",
            "(vs_plasma_burn_required)",
            physics_variables.vs_plasma_burn_required,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Plasma loop voltage during burn (V)",
            "(v_plasma_loop_burn)",
            physics_variables.v_plasma_loop_burn,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Coefficient for sawtooth effects on burn V-s requirement",
            "(csawth)",
            physics_variables.csawth,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma resistance (ohm)",
            "(res_plasma)",
            physics_variables.res_plasma,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Plasma resistive diffusion time (s)",
            "(t_plasma_res_diffusion)",
            physics_variables.t_plasma_res_diffusion,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma inductance (H)",
            "(ind_plasma)",
            physics_variables.ind_plasma,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma magnetic energy stored (J)",
            "(e_plasma_magnetic_stored)",
            physics_variables.e_plasma_magnetic_stored,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Plasma normalised internal inductance",
            "(ind_plasma_internal_norm)",
            physics_variables.ind_plasma_internal_norm,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)


class DetailedPhysics(Model):
    """Class to hold detailed physics models for plasma processing."""

    def __init__(self, plasma_profile):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile

    def output(self):
        self.run()
        self.output_detailed_physics()

    def run(self):
        # ---------------------------
        #  Debye length calculation
        # ---------------------------

        physics_variables.len_plasma_debye_electron_vol_avg = self.calculate_debye_length(
            temp_plasma_species_kev=physics_variables.temp_plasma_electron_vol_avg_kev,
            nd_plasma_species=physics_variables.nd_plasma_electrons_vol_avg,
        )

        physics_variables.len_plasma_debye_electron_profile = (
            self.calculate_debye_length(
                temp_plasma_species_kev=self.plasma_profile.teprofile.profile_y,
                nd_plasma_species=self.plasma_profile.neprofile.profile_y,
            )
        )

        # ============================
        # Particle relativistic speeds
        # ============================

        physics_variables.vel_plasma_electron_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT,
                mass=constants.ELECTRON_MASS,
            )
        )

        physics_variables.vel_plasma_deuteron_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT
                * physics_variables.f_temp_plasma_ion_electron,
                mass=constants.DEUTERON_MASS,
            )
        )

        physics_variables.vel_plasma_triton_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT
                * physics_variables.f_temp_plasma_ion_electron,
                mass=constants.TRITON_MASS,
            )
        )

        physics_variables.vel_plasma_alpha_thermal_profile = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=self.plasma_profile.teprofile.profile_y
                * constants.KILOELECTRON_VOLT
                * physics_variables.f_temp_plasma_ion_electron,
                mass=constants.ALPHA_MASS,
            )
        )

        physics_variables.vel_plasma_alpha_birth = (
            self.calculate_relativistic_particle_speed(
                e_kinetic=constants.DT_ALPHA_ENERGY,
                mass=constants.ALPHA_MASS,
            )
        )

        # ============================
        # Plasma frequencies
        # ============================

        physics_variables.freq_plasma_electron_profile = self.calculate_plasma_frequency(
            nd_particle=self.plasma_profile.neprofile.profile_y,
            m_particle=constants.ELECTRON_MASS,
            z_particle=1.0,
        )

        # ============================
        # Larmor frequencies
        # ============================

        physics_variables.freq_plasma_larmor_toroidal_electron_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.ELECTRON_MASS,
                z_particle=1.0,
            )
        )

        physics_variables.freq_plasma_larmor_toroidal_deuteron_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.DEUTERON_MASS,
                z_particle=1.0,
            )
        )

        physics_variables.freq_plasma_larmor_toroidal_triton_profile = (
            self.calculate_larmor_frequency(
                b_field=physics_variables.b_plasma_toroidal_profile,
                m_particle=constants.TRITON_MASS,
                z_particle=1.0,
            )
        )

        # ============================
        # Upper hybrid frequencies
        # ============================

        physics_variables.freq_plasma_upper_hybrid_profile = self.calculate_upper_hybrid_frequency(
            freq_plasma=np.concatenate([
                physics_variables.freq_plasma_electron_profile[::-1],
                physics_variables.freq_plasma_electron_profile,
            ]),
            freq_larmor=physics_variables.freq_plasma_larmor_toroidal_electron_profile,
        )

        # ============================
        # Larmor radii
        # ============================

        # Since isotropic (v⟂)² = 2(v)² for a Maxwellian distribution,
        # we can use the total velocity to calculate the Larmor radius for an isotropic profile
        physics_variables.radius_plasma_deuteron_toroidal_larmor_isotropic_profile = self.calculate_larmor_radius(
            vel_perp=np.sqrt(
                2
                * np.concatenate([
                    physics_variables.vel_plasma_deuteron_profile[::-1],
                    physics_variables.vel_plasma_deuteron_profile,
                ])
                ** 2
            ),
            freq_larmor=physics_variables.freq_plasma_larmor_toroidal_deuteron_profile
            * (2 * np.pi),
        )

        # Since isotropic (v⟂)² = 2(v)² for a Maxwellian distribution,
        # we can use the total velocity to calculate the Larmor radius for an isotropic profile
        physics_variables.radius_plasma_triton_toroidal_larmor_isotropic_profile = (
            self.calculate_larmor_radius(
                vel_perp=np.sqrt(
                    2
                    * np.concatenate([
                        physics_variables.vel_plasma_triton_profile[::-1],
                        physics_variables.vel_plasma_triton_profile,
                    ])
                    ** 2
                ),
                freq_larmor=physics_variables.freq_plasma_larmor_toroidal_triton_profile
                * (2 * np.pi),
            )
        )

        # ============================
        # Coulomb logarithm
        # ============================

        physics_variables.plasma_coulomb_log_electron_electron_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_electron_deuteron_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.DEUTERON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_deuteron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.DEUTERON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_deuteron_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_electron_triton_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.ELECTRON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_electron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_deuteron_triton_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=1,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.DEUTERON_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_deuteron_profile[i],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.TRITON_MASS,
                            mass2=constants.DEUTERON_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_triton_profile[i],
                            velocity_2=physics_variables.vel_plasma_deuteron_profile[i],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        physics_variables.plasma_coulomb_log_electron_alpha_thermal_profile = np.array([
            self.calculate_coulomb_log_from_impact(
                impact_param_max=physics_variables.len_plasma_debye_electron_profile[i],
                impact_param_min=max(
                    self.calculate_classical_distance_of_closest_approach(
                        charge1=1,
                        charge2=2,
                        m_reduced=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ALPHA_MASS,
                        ),
                        vel_relative=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_alpha_thermal_profile[
                                i
                            ],
                        ),
                    ),
                    self.calculate_debroglie_wavelength(
                        mass=self.calculate_reduced_mass(
                            mass1=constants.ELECTRON_MASS,
                            mass2=constants.ALPHA_MASS,
                        ),
                        velocity=self.calculate_average_relative_velocity(
                            velocity_1=physics_variables.vel_plasma_electron_profile[i],
                            velocity_2=physics_variables.vel_plasma_alpha_thermal_profile[
                                i
                            ],
                        ),
                    ),
                ),
            )
            for i in range(len(physics_variables.len_plasma_debye_electron_profile))
        ])

        # ============================
        # Collision times
        # ============================

        physics_variables.t_plasma_electron_electron_collision_profile = self.calculate_electron_electron_collision_time(
            temp_plasma_electron_kev=self.plasma_profile.teprofile.profile_y,
            nd_plasma_electrons=self.plasma_profile.neprofile.profile_y,
            plasma_coulomb_log_electron_electron=physics_variables.plasma_coulomb_log_electron_electron_profile,
        )

        physics_variables.t_plasma_electron_deuteron_collision_profile = self.calculate_electron_ion_collision_time(
            temp_plasma_electron_kev=self.plasma_profile.teprofile.profile_y,
            nd_plasma_ions=(
                self.plasma_profile.neprofile.profile_y
                * physics_variables.f_plasma_fuel_deuterium
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            ),
            plasma_coulomb_log_electron_ion=physics_variables.plasma_coulomb_log_electron_deuteron_profile,
            n_charge_ion=1,
        )

        physics_variables.t_plasma_electron_triton_collision_profile = self.calculate_electron_ion_collision_time(
            temp_plasma_electron_kev=self.plasma_profile.teprofile.profile_y,
            nd_plasma_ions=(
                self.plasma_profile.neprofile.profile_y
                * physics_variables.f_plasma_fuel_tritium
                * (
                    physics_variables.nd_plasma_fuel_ions_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            ),
            plasma_coulomb_log_electron_ion=physics_variables.plasma_coulomb_log_electron_triton_profile,
            n_charge_ion=1,
        )

        physics_variables.t_plasma_electron_alpha_thermal_collision_profile = self.calculate_electron_ion_collision_time(
            temp_plasma_electron_kev=self.plasma_profile.teprofile.profile_y,
            nd_plasma_ions=(
                self.plasma_profile.neprofile.profile_y
                * (
                    physics_variables.nd_plasma_alphas_vol_avg
                    / physics_variables.nd_plasma_electrons_vol_avg
                )
            ),
            plasma_coulomb_log_electron_ion=physics_variables.plasma_coulomb_log_electron_alpha_thermal_profile,
            n_charge_ion=2,
        )

        # ============================
        # Collision frequencies
        # ============================

        physics_variables.freq_plasma_electron_electron_collision_profile = (
            1 / physics_variables.t_plasma_electron_electron_collision_profile
        )

        physics_variables.freq_plasma_electron_deuteron_collision_profile = (
            1 / physics_variables.t_plasma_electron_deuteron_collision_profile
        )

        physics_variables.freq_plasma_electron_triton_collision_profile = (
            1 / physics_variables.t_plasma_electron_triton_collision_profile
        )

        physics_variables.freq_plasma_electron_alpha_thermal_collision_profile = (
            1 / physics_variables.t_plasma_electron_alpha_thermal_collision_profile
        )

        # ============================
        # Mean free paths
        # ============================

        physics_variables.len_plasma_electron_electron_mean_free_path_profile = (
            physics_variables.t_plasma_electron_electron_collision_profile
            * physics_variables.vel_plasma_electron_profile
        )

        physics_variables.len_plasma_electron_deuteron_mean_free_path_profile = (
            physics_variables.t_plasma_electron_deuteron_collision_profile
            * physics_variables.vel_plasma_electron_profile
        )

        physics_variables.len_plasma_electron_triton_mean_free_path_profile = (
            physics_variables.t_plasma_electron_triton_collision_profile
            * physics_variables.vel_plasma_electron_profile
        )

        physics_variables.len_plasma_electron_alpha_thermal_mean_free_path_profile = (
            physics_variables.t_plasma_electron_alpha_thermal_collision_profile
            * physics_variables.vel_plasma_electron_profile
        )

        # ============================
        # Spitzer slow down time
        # ============================

        physics_variables.t_plasma_electron_alpha_spitzer_slow_profile = self.calculate_spitzer_ion_slowing_down_time(
            m_ion=constants.ALPHA_MASS,
            plasma_coulomb_log_electron_ion=physics_variables.plasma_coulomb_log_electron_alpha_thermal_profile,
            n_charge_ion=2,
        )

        # ============================
        # Resistivites
        # ============================

        physics_variables.res_plasma_spitzer_profile = self.calculate_spitzer_resistivity(
            n_charge=physics_variables.n_charge_plasma_effective_profile,
            electron_ion_coulomb_log=physics_variables.plasma_coulomb_log_electron_deuteron_profile,
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_debye_length(
        temp_plasma_species_kev: float | np.ndarray,
        nd_plasma_species: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the Debye length for a plasma.

        Parameters
        ----------
        temp_plasma_species_kev : float | np.ndarray
            Species temperature in keV.
        nd_plasma_species : float | np.ndarray
            Species number density (/m^3).

        Returns
        -------
        float | np.ndarray
            Debye length in meters.

        """
        return (
            (constants.EPSILON0 * temp_plasma_species_kev * constants.KILOELECTRON_VOLT)
            / (nd_plasma_species * constants.ELECTRON_CHARGE**2)
        ) ** 0.5

    @staticmethod
    @nb.njit(cache=True)
    def calculate_lorentz_factor(velocity: float | np.ndarray) -> float | np.ndarray:
        """Calculate the Lorentz factor for a given velocity.

        Parameters
        ----------
        velocity : float | np.ndarray
            Velocity in m/s.

        Returns
        -------
        float | np.ndarray
            Lorentz factor (dimensionless).

        """
        return 1 / (1 - (velocity / constants.SPEED_LIGHT) ** 2) ** 0.5

    @staticmethod
    @nb.njit(cache=True)
    def calculate_relativistic_particle_speed(
        e_kinetic: float | np.ndarray, mass: float
    ) -> float | np.ndarray:
        """Calculate the speed of a particle given its kinetic energy and mass using relativistic mechanics.

        Parameters
        ----------
        e_kinetic : float | np.ndarray
            Kinetic energy in Joules.
        mass : float
            Mass of the particle in kg.

        Returns
        -------
        float | np.ndarray
            Speed of the particle in m/s.

        """
        return (
            constants.SPEED_LIGHT
            * (1 - (1 / ((e_kinetic / (mass * constants.SPEED_LIGHT**2)) + 1) ** 2))
            ** 0.5
        )

    def calculate_coulomb_log_from_impact(
        self, impact_param_max: float, impact_param_min: float
    ) -> float:
        """Calculate the Coulomb logarithm from maximum and minimum impact parameters.

        Parameters
        ----------
        impact_param_max : float
            Maximum impact parameter in meters.
        impact_param_min : float
            Minimum impact parameter in meters.

        Returns
        -------
        float
            Coulomb logarithm (dimensionless).

        """
        return np.log(impact_param_max / impact_param_min)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_classical_distance_of_closest_approach(
        charge1: float,
        charge2: float,
        m_reduced: float,
        vel_relative: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the classical distance of closest approach for two charged particles.

        Parameters
        ----------
        charge1 :
            Charge of particle 1 in units of elementary charge.
        charge2 :
            Charge of particle 2 in units of elementary charge.
        m_reduced:
            Reduced mass of the two-particle system in kg.
        vel_relative:
            Relative velocity of the two particles in m/s.

        Returns
        -------
        float | np.ndarray
            Distance of closest approach in meters.
        """
        return (charge1 * charge2 * constants.ELECTRON_CHARGE**2) / (
            2 * np.pi * constants.EPSILON0 * m_reduced * vel_relative**2
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_debroglie_wavelength(
        mass: float, velocity: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the de Broglie wavelength of a particle.

        Parameters
        ----------
        mass : float
            Mass of the particle in kg.
        velocity : float | np.ndarray
            Velocity of the particle in m/s.

        Returns
        -------
        float | np.ndarray
            de Broglie wavelength in meters.

        :note: Reduced Planck constant (h-bar) is used in the calculation as this is for scattering.

        """
        return (constants.PLANCK_CONSTANT / (2 * np.pi)) / (mass * velocity)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_plasma_frequency(
        nd_particle: float | np.ndarray, m_particle: float, z_particle: float
    ) -> float | np.ndarray:
        """Calculate the plasma frequency for a particle species.

        Parameters
        ----------
        nd_particle : float | np.ndarray
            Number density of the particle species (/m^3).
        m_particle : float
            Mass of the particle species (kg).
        Z_particle : float
            Charge state of the particle species (dimensionless).

        Returns
        -------
        float | np.ndarray
            Plasma frequency in Hz.

        """
        return (
            (
                (nd_particle * z_particle**2 * constants.ELECTRON_CHARGE**2)
                / (m_particle * constants.EPSILON0)
            )
            ** 0.5
        ) / (2 * np.pi)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_larmor_frequency(
        b_field: float | np.ndarray, m_particle: float, z_particle: float
    ) -> float | np.ndarray:
        """Calculate the Larmor frequency for a particle species.

        Parameters
        ----------
        b_field : float | np.ndarray
            Magnetic field strength (T).
        m_particle : float
            Mass of the particle species (kg).
        Z_particle : float
            Charge state of the particle species (dimensionless).

        Returns
        -------
        float | np.ndarray
            Larmor frequency in Hz.

        """
        return (z_particle * constants.ELECTRON_CHARGE * b_field) / (
            2 * np.pi * m_particle
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_upper_hybrid_frequency(
        freq_plasma: float | np.ndarray, freq_larmor: float | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the upper hybrid frequency for a particle species.

        Parameters
        ----------
        freq_plasma : float | np.ndarray
            Plasma frequency of the particle species (Hz).
        freq_larmor : float | np.ndarray
            Larmor frequency of the particle species (Hz).

        Returns
        -------
        float | np.ndarray
            Upper hybrid frequency in Hz.
        """
        return np.sqrt(freq_plasma**2 + freq_larmor**2)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_larmor_radius(
        vel_perp: float | np.ndarray,
        freq_larmor: float,
    ) -> float | np.ndarray:
        """Calculate the Larmor radius for a particle species.

        Parameters
        ----------
        vel_perp : float | np.ndarray
            Perpendicular velocity of the particle to the magnetic field (m/s).
        freq_larmor : float
            Larmor frequency of the particle (Hz).

        Returns
        -------
        float | np.ndarray
            Larmor radius in meters.
        """
        return vel_perp / (freq_larmor)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_reduced_mass(mass1: float, mass2: float) -> float:
        """
        Calculate the reduced mass of two particles.

        Parameters
        ----------
        mass1:
            Mass of particle 1 (kg).
        mass2:
            Mass of particle 2 (kg).

        Returns
        -------
            Reduced mass (kg).
        """
        return (mass1 * mass2) / (mass1 + mass2)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_average_relative_velocity(
        velocity_1: float | np.ndarray, velocity_2: float | np.ndarray
    ) -> float | np.ndarray:
        """
        Calculate the average relative velocity between two particles.

        Parameters
        ----------
        velocity_1:
            Velocity of particle 1 (m/s).
        velocity_2:
            Velocity of particle 2 (m/s).

        Returns
        -------
            Average relative velocity (m/s).
        """
        return np.sqrt(velocity_1**2 + velocity_2**2)

    @staticmethod
    @nb.njit(cache=True)
    def calculate_electron_electron_collision_time(
        temp_plasma_electron_kev: float | np.ndarray,
        nd_plasma_electrons: float | np.ndarray,
        plasma_coulomb_log_electron_electron: float | np.ndarray,
    ) -> float | np.ndarray:
        """Calculate the electron-electron collision time for a plasma (τₑₑ).

        Parameters
        ----------
        temp_plasma_electron_kev : float | np.ndarray
            Electron temperature in keV.
        nd_plasma_electrons : float | np.ndarray
            Electron density (m^-3).
        plasma_coulomb_log_electron_electron : float | np.ndarray
            Coulomb logarithm for electron-electron collisions.

        Returns
        -------
        float | np.ndarray
            Electron-electron collision time (s).
        """
        return (
            12
            * np.sqrt(2)
            * np.pi**1.5
            * constants.EPSILON0**2
            * constants.ELECTRON_MASS**0.5
            * (temp_plasma_electron_kev * constants.KILOELECTRON_VOLT) ** 1.5
        ) / (
            plasma_coulomb_log_electron_electron
            * constants.ELECTRON_CHARGE**4
            * nd_plasma_electrons
        )

    @staticmethod
    @nb.njit(cache=True)
    def calculate_electron_ion_collision_time(
        temp_plasma_electron_kev: float | np.ndarray,
        nd_plasma_ions: float | np.ndarray,
        plasma_coulomb_log_electron_ion: float | np.ndarray,
        n_charge_ion: float = 1.0,
    ) -> float | np.ndarray:
        """Calculate the electron-ion collision time for a plasma (τₑᵢ).

        Parameters
        ----------
        temp_plasma_electron_kev : float | np.ndarray
            Electron temperature in keV.
        nd_plasma_ions : float | np.ndarray
            Ion density (m^-3).
        plasma_coulomb_log_electron_ion : float | np.ndarray
            Coulomb logarithm for electron-ion collisions.
        n_charge_ion : float
            Charge number (Z) of the ion.

        Returns
        -------
        float | np.ndarray
            Electron-ion collision time (s).
        """
        return (
            12
            * np.pi**1.5
            * constants.EPSILON0**2
            * constants.ELECTRON_MASS**0.5
            * (temp_plasma_electron_kev * constants.KILOELECTRON_VOLT) ** 1.5
        ) / (
            np.sqrt(2)
            * nd_plasma_ions
            * n_charge_ion**2
            * plasma_coulomb_log_electron_ion
            * constants.ELECTRON_CHARGE**4
        )

    def calculate_spitzer_ion_slowing_down_time(
        self,
        m_ion: float,
        plasma_coulomb_log_electron_ion: float | np.ndarray,
        n_charge_ion: float = 1.0,
    ) -> float:
        """
        Calculate the Spitzer slowing down time for ions in a plasma.

        Parameters
        ----------
        m_ion : float
            Mass of the ion (kg).
        plasma_coulomb_log_electron_ion : float | np.ndarray
            Coulomb logarithm for electron-ion collisions.
        n_charge_ion : float
            Charge number (Z) of the ion.

        Returns
        -------
        float
            Spitzer slowing down time (s).
        """

        return (
            (3 * (2 * np.pi) ** 1.5 * constants.EPSILON0**2)
            * (
                m_ion
                * (self.plasma_profile.teprofile.profile_y * constants.KILOELECTRON_VOLT)
                ** 1.5
            )
            / (
                self.plasma_profile.neprofile.profile_y
                * constants.ELECTRON_CHARGE**4
                * plasma_coulomb_log_electron_ion
                * n_charge_ion**2
                * np.sqrt(constants.ELECTRON_MASS)
            )
        )

    def calculate_spitzer_resistivity(
        self, n_charge: float | np.ndarray, electron_ion_coulomb_log: float | np.ndarray
    ) -> float | np.ndarray:
        """
        Calculate the classical Spitzer resistivity for a plasma.

        Parameters
        ----------
        n_charge : float | np.ndarray
            Charge number (Z) of the ion.
        electron_ion_coulomb_log : float | np.ndarray
            Coulomb logarithm for electron-ion collisions.

        Returns
        -------
        float | np.ndarray
            Spitzer resistivity (Ohm m).
        """

        return (
            (4 * np.sqrt(2 * np.pi) / 3)
            * (1 / (4 * np.pi * constants.EPSILON0) ** 2)
            * (
                n_charge
                * constants.ELECTRON_CHARGE**2
                * constants.ELECTRON_MASS**0.5
                * electron_ion_coulomb_log
            )
            / (self.plasma_profile.teprofile.profile_y * constants.KILOELECTRON_VOLT)
            ** 1.5
        )

    def output_detailed_physics(self):
        """Outputs detailed physics variables to file."""
        po.oheadr(self.outfile, "Detailed Plasma")

        po.osubhd(self.outfile, "Debye lengths:")

        po.ovarrf(
            self.outfile,
            "Plasma volume averaged electron Debye length (m)",
            "(len_plasma_debye_electron_vol_avg)",
            physics_variables.len_plasma_debye_electron_vol_avg,
            "OP ",
        )
        for i in range(len(physics_variables.len_plasma_debye_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron Debye length at point {i}",
                f"(len_plasma_debye_electron_profile{i})",
                physics_variables.len_plasma_debye_electron_profile[i],
            )

        po.osubhd(self.outfile, "Larmor radii:")

        for i in range(
            len(
                physics_variables.radius_plasma_deuteron_toroidal_larmor_isotropic_profile
            )
        ):
            po.ovarre(
                self.mfile,
                f"Plasma deuteron isotropic Larmor radius at point {i}",
                f"(radius_plasma_deuteron_toroidal_larmor_isotropic_profile{i})",
                physics_variables.radius_plasma_deuteron_toroidal_larmor_isotropic_profile[
                    i
                ],
            )
        for i in range(
            len(physics_variables.radius_plasma_triton_toroidal_larmor_isotropic_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma triton isotropic Larmor radius at point {i}",
                f"(radius_plasma_triton_toroidal_larmor_isotropic_profile{i})",
                physics_variables.radius_plasma_triton_toroidal_larmor_isotropic_profile[
                    i
                ],
            )

        po.osubhd(self.outfile, "Velocities:")

        for i in range(len(physics_variables.vel_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron thermal velocity at point {i}",
                f"(vel_plasma_electron_profile{i})",
                physics_variables.vel_plasma_electron_profile[i],
            )
        for i in range(len(physics_variables.vel_plasma_deuteron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma deuteron thermal velocity at point {i}",
                f"(vel_plasma_deuteron_profile{i})",
                physics_variables.vel_plasma_deuteron_profile[i],
            )
        for i in range(len(physics_variables.vel_plasma_triton_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma triton thermal velocity at point {i}",
                f"(vel_plasma_triton_profile{i})",
                physics_variables.vel_plasma_triton_profile[i],
            )
        for i in range(len(physics_variables.vel_plasma_alpha_thermal_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma alpha thermal velocity at point {i}",
                f"(vel_plasma_alpha_thermal_profile{i})",
                physics_variables.vel_plasma_alpha_thermal_profile[i],
            )
        po.ovarre(
            self.outfile,
            "Plasma alpha birth velocity (m/s)",
            "(vel_plasma_alpha_birth)",
            physics_variables.vel_plasma_alpha_birth,
        )

        po.osubhd(self.outfile, "Frequencies:")

        for i in range(len(physics_variables.freq_plasma_electron_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma electron frequency at point {i}",
                f"(freq_plasma_electron_profile{i})",
                physics_variables.freq_plasma_electron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_electron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma electron toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_electron_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_electron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_deuteron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma deuteron toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_deuteron_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_deuteron_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_larmor_toroidal_triton_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Plasma triton toroidal Larmor frequency at point {i}",
                f"(freq_plasma_larmor_toroidal_triton_profile{i})",
                physics_variables.freq_plasma_larmor_toroidal_triton_profile[i],
            )

        for i in range(len(physics_variables.freq_plasma_upper_hybrid_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma upper hybrid frequency at point {i}",
                f"(freq_plasma_upper_hybrid_profile{i})",
                physics_variables.freq_plasma_upper_hybrid_profile[i],
            )

        po.osubhd(self.outfile, "Coulomb Logarithms:")

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_electron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-electron Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_electron_profile{i})",
                physics_variables.plasma_coulomb_log_electron_electron_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_deuteron_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-deuteron Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_deuteron_profile{i})",
                physics_variables.plasma_coulomb_log_electron_deuteron_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_triton_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-triton Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_triton_profile{i})",
                physics_variables.plasma_coulomb_log_electron_triton_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_deuteron_triton_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Deuteron-triton Coulomb log at point {i}",
                f"(plasma_coulomb_log_deuteron_triton_profile{i})",
                physics_variables.plasma_coulomb_log_deuteron_triton_profile[i],
            )

        for i in range(
            len(physics_variables.plasma_coulomb_log_electron_alpha_thermal_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-alpha thermal Coulomb log at point {i}",
                f"(plasma_coulomb_log_electron_alpha_thermal_profile{i})",
                physics_variables.plasma_coulomb_log_electron_alpha_thermal_profile[i],
            )

        po.osubhd(self.outfile, "Collision Times:")

        for i in range(
            len(physics_variables.t_plasma_electron_electron_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-electron collision time at point {i}",
                f"(t_plasma_electron_electron_collision_profile{i})",
                physics_variables.t_plasma_electron_electron_collision_profile[i],
            )

        for i in range(
            len(physics_variables.t_plasma_electron_deuteron_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-deuteron collision time at point {i}",
                f"(t_plasma_electron_deuteron_collision_profile{i})",
                physics_variables.t_plasma_electron_deuteron_collision_profile[i],
            )

        for i in range(
            len(physics_variables.t_plasma_electron_triton_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-triton collision time at point {i}",
                f"(t_plasma_electron_triton_collision_profile{i})",
                physics_variables.t_plasma_electron_triton_collision_profile[i],
            )

        for i in range(
            len(physics_variables.t_plasma_electron_alpha_thermal_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-alpha thermal collision time at point {i}",
                f"(t_plasma_electron_alpha_thermal_collision_profile{i})",
                physics_variables.t_plasma_electron_alpha_thermal_collision_profile[i],
            )

        po.osubhd(self.outfile, "Collision Frequencies:")

        for i in range(
            len(physics_variables.freq_plasma_electron_electron_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-electron collision frequency at point {i}",
                f"(freq_plasma_electron_electron_collision_profile{i})",
                physics_variables.freq_plasma_electron_electron_collision_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_electron_deuteron_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-deuteron collision frequency at point {i}",
                f"(freq_plasma_electron_deuteron_collision_profile{i})",
                physics_variables.freq_plasma_electron_deuteron_collision_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_electron_triton_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-triton collision frequency at point {i}",
                f"(freq_plasma_electron_triton_collision_profile{i})",
                physics_variables.freq_plasma_electron_triton_collision_profile[i],
            )
        for i in range(
            len(physics_variables.freq_plasma_electron_alpha_thermal_collision_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-alpha thermal collision frequency at point {i}",
                f"(freq_plasma_electron_alpha_thermal_collision_profile{i})",
                physics_variables.freq_plasma_electron_alpha_thermal_collision_profile[
                    i
                ],
            )

        po.osubhd(self.outfile, "Mean Free Paths:")

        for i in range(
            len(physics_variables.len_plasma_electron_electron_mean_free_path_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-electron mean free path at point {i}",
                f"(len_plasma_electron_electron_mean_free_path_profile{i})",
                physics_variables.len_plasma_electron_electron_mean_free_path_profile[i],
            )
        for i in range(
            len(physics_variables.len_plasma_electron_deuteron_mean_free_path_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-deuteron mean free path at point {i}",
                f"(len_plasma_electron_deuteron_mean_free_path_profile{i})",
                physics_variables.len_plasma_electron_deuteron_mean_free_path_profile[i],
            )
        for i in range(
            len(physics_variables.len_plasma_electron_triton_mean_free_path_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-triton mean free path at point {i}",
                f"(len_plasma_electron_triton_mean_free_path_profile{i})",
                physics_variables.len_plasma_electron_triton_mean_free_path_profile[i],
            )
        for i in range(
            len(
                physics_variables.len_plasma_electron_alpha_thermal_mean_free_path_profile
            )
        ):
            po.ovarre(
                self.mfile,
                f"Electron-alpha thermal mean free path at point {i}",
                f"(len_plasma_electron_alpha_thermal_mean_free_path_profile{i})",
                physics_variables.len_plasma_electron_alpha_thermal_mean_free_path_profile[
                    i
                ],
            )

        po.osubhd(self.outfile, "Spitzer slowing down times:")

        for i in range(
            len(physics_variables.t_plasma_electron_alpha_spitzer_slow_profile)
        ):
            po.ovarre(
                self.mfile,
                f"Electron-alpha Spitzer slowing down time at point {i}",
                f"(t_plasma_electron_alpha_spitzer_slow_profile{i})",
                physics_variables.t_plasma_electron_alpha_spitzer_slow_profile[i],
            )

        po.osubhd(self.outfile, "Resistivities:")

        for i in range(len(physics_variables.res_plasma_spitzer_profile)):
            po.ovarre(
                self.mfile,
                f"Plasma Spitzer resistivity at point {i}",
                f"(res_plasma_spitzer_profile{i})",
                physics_variables.res_plasma_spitzer_profile[i],
            )
