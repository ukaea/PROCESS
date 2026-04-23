import logging
from enum import IntEnum

import numpy as np
from scipy.optimize import root_scalar

from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.data_structure import (
    constraint_variables,
    current_drive_variables,
    physics_variables,
    stellarator_variables,
)

logger = logging.getLogger(__name__)


class ConfinementTimeModel(IntEnum):
    """Confinement time (τ_E) model types"""

    USER_INPUT = (0, "User input electron confinement   ")
    NEO_ALCATOR = (1, "Neo-Alcator                (Ohmic)")
    MIRNOV = (2, "Mirnov                         (H)")
    MEREZHKIN_MUHKOVATOV = (3, "Merezkhin-Muhkovatov    (Ohmic)(L)")
    SHIMOMURA = (4, "Shimomura                      (H)")
    KAYE_GOLDSTON = (5, "Kaye-Goldston                  (L)")
    ITER_89P = (6, "ITER 89-P                      (L)")
    ITER_89_0 = (7, "ITER 89-O                      (L)")
    REBUT_LALLIA = (8, "Rebut-Lallia                   (L)")
    GOLDSTON = (9, "Goldston                       (L)")
    T_10 = (10, "T10                            (L)")
    JAERI = (11, "JAERI / Odajima-Shimomura      (L)")
    KAYE_BIG = (12, "Kaye-Big Complex               (L)")
    ITER_H90_P = (13, "ITER H90-P                     (H)")
    MINIMUM_OF_ITER_89P_AND_ITER_89_0 = (14, "ITER 89-P & 89-O min           (L)")
    RIEDEL_L = (15, "Riedel                         (L)")
    CHRISTIANSEN = (16, "Christiansen                   (L)")
    LACKNER_GOTTARDI = (17, "Lackner-Gottardi               (L)")
    NEO_KAYE = (18, "Neo-Kaye                       (L)")
    RIEDEL_H = (19, "Riedel                         (H)")
    ITER_H90_P_AMENDED = (20, "ITER H90-P amended             (H)")
    SUDO_ET_AL = (21, "LHD                        (Stell)")
    GYRO_REDUCED_BOHM = (22, "Gyro-reduced Bohm          (Stell)")
    LACKNER_GOTTARDI_STELLARATOR = (23, "Lackner-Gottardi           (Stell)")
    ITER_93H = (24, "ITER-93H  ELM-free             (H)")
    TITAN_REMOVED = (25, "TITAN RFP OBSOLETE                ")
    ITER_H97P = (26, "ITER H-97P ELM-free            (H)")
    ITER_H97P_ELMY = (27, "ITER H-97P ELMy                (H)")
    ITER_96P = (28, "ITER-96P (ITER-97L)            (L)")
    VALOVIC_ELMY = (29, "Valovic modified ELMy          (H)")
    KAYE = (30, "Kaye 98 modified               (L)")
    ITER_PB98P_Y = (31, "ITERH-PB98P(y)                 (H)")
    IPB98_Y = (32, "IPB98(y)                       (H)")
    ITER_IPB98Y1 = (33, "IPB98(y,1)                     (H)")
    ITER_IPB98Y2 = (34, "IPB98(y,2)                     (H)")
    ITER_IPB98Y3 = (35, "IPB98(y,3)                     (H)")
    ITER_IPB98Y4 = (36, "IPB98(y,4)                     (H)")
    ISS95_STELLARATOR = (37, "ISS95                      (Stell)")
    ISS04_STELLARATOR = (38, "ISS04                      (Stell)")
    DS03 = (39, "DS03 beta-independent          (H)")
    MURARI = (40, 'Murari "Non-power law"         (H)')
    PETTY08 = (41, "Petty 2008                 (ST)(H)")
    LANG_HIGH_DENSITY = (42, "Lang high density              (H)")
    HUBBARD_NOMINAL = (43, "Hubbard 2017 - nominal         (I)")
    HUBBARD_LOWER = (44, "Hubbard 2017 - lower           (I)")
    HUBBARD_UPPER = (45, "Hubbard 2017 - upper           (I)")
    MENARD_NSTX = (46, "Menard NSTX                (ST)(H)")
    MENARD_NSTX_PETTY08_HYBRID = (47, "Menard NSTX-Petty08 hybrid (ST)(H)")
    NSTX_GYRO_BOHM = (48, "Buxton NSTX gyro-Bohm      (ST)(H)")
    ITPA20 = (49, "ITPA20                         (H)")
    ITPA20_IL = (50, "ITPA20-IL                      (H)")

    def __new__(cls, value, full_name):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.full_name = full_name
        return obj


class ConfinementRadiationLossModel(IntEnum):
    """Confinement radiation loss model types"""

    FULL_RADIATION = (0, "All radiation included in loss power term")
    CORE_ONLY = (1, "Only core radiation included in loss power term")
    NO_RADIATION = (2, "No radiation included in loss power term")

    def __new__(cls, value, description):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.description = description
        return obj


class PlasmaConfinementTime:
    """Class to calculate plasma confinement time using various empirical scaling laws"""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def calculate_confinement_time(
        self,
        m_fuel_amu: float,
        p_alpha_total_mw: float,
        aspect: float,
        b_plasma_toroidal_on_axis: float,
        nd_plasma_ions_total_vol_avg: float,
        nd_plasma_electrons_vol_avg: float,
        nd_plasma_electron_line: float,
        eps: float,
        hfact: float,
        i_confinement_time: int,
        i_plasma_ignited: int,
        kappa: float,
        kappa95: float,
        p_non_alpha_charged_mw: float,
        p_hcd_injected_total_mw: float,
        plasma_current: float,
        pden_plasma_core_rad_mw: float,
        rmajor: float,
        rminor: float,
        temp_plasma_electron_density_weighted_kev: float,
        temp_plasma_ion_density_weighted_kev: float,
        q95: float,
        qstar: float,
        vol_plasma: float,
        zeff: float,
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate the confinement times and the transport power loss terms.

        Parameters
        ----------
        m_fuel_amu :
            Average mass of fuel (amu)
        p_alpha_total_mw :
            Alpha particle power (MW)
        aspect :
            Aspect ratio
        b_plasma_toroidal_on_axis :
            Toroidal field on axis (T)
        nd_plasma_ions_total_vol_avg :
            Total ion density (/m3)
        nd_plasma_electrons_vol_avg :
            Volume averaged electron density (/m3)
        nd_plasma_electron_line :
            Line-averaged electron density (/m3)
        eps :
            Inverse aspect ratio
        hfact :
            H factor on energy confinement scalings
        i_confinement_time :
            Switch for energy confinement scaling to use
        i_plasma_ignited :
            Switch for ignited calculation
        kappa :
            Plasma elongation
        kappa95 :
            Plasma elongation at 95% surface
        p_non_alpha_charged_mw :
            Non-alpha charged particle fusion power (MW)
        p_hcd_injected_total_mw :
            Auxiliary power to ions and electrons (MW)
        plasma_current :
            Plasma current (A)
        pden_plasma_core_rad_mw :
            Total core radiation power (MW/m3)
        q95 :
            Edge safety factor (tokamaks), or rotational transform iotabar (stellarators)
        qstar :
            Equivalent cylindrical edge safety factor
        rmajor :
            Plasma major radius (m)
        rminor :
            Plasma minor radius (m)
        temp_plasma_electron_density_weighted_kev :
            Density weighted average electron temperature (keV)
        temp_plasma_ion_density_weighted_kev :
            Density weighted average ion temperature (keV)
        vol_plasma :
            Plasma volume (m3)
        a_plasma_poloidal :
            Plasma cross-sectional area (m2)
        zeff :
            Plasma effective charge

        Returns
        -------
        type
            Tuple containing:
            - pden_electron_transport_loss_mw (float): Electron transport power (MW/m3)
            - pden_ion_transport_loss_mw (float): Ion transport power (MW/m3)
            - t_electron_energy_confinement (float): Electron energy confinement time (s)
            - t_ion_energy_confinement (float): Ion energy confinement time (s)
            - t_energy_confinement (float): Global energy confinement time (s)
            - p_plasma_loss_mw (float): Heating power (MW) assumed in calculation

        """
        # ========================================================================

        # Calculate heating power (MW)
        p_plasma_loss_mw = (
            physics_variables.f_p_alpha_plasma_deposited * p_alpha_total_mw
            + p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
        )

        # If the device is not ignited, add the injected auxiliary power
        if i_plasma_ignited == 0:
            p_plasma_loss_mw += p_hcd_injected_total_mw

        # Include the radiation as a loss term based on radiation model
        try:
            model = ConfinementRadiationLossModel(int(physics_variables.i_rad_loss))

            if model == ConfinementRadiationLossModel.FULL_RADIATION:
                p_plasma_loss_mw -= physics_variables.pden_plasma_rad_mw * vol_plasma
            elif model == ConfinementRadiationLossModel.CORE_ONLY:
                p_plasma_loss_mw -= pden_plasma_core_rad_mw * vol_plasma
            # NO_RADIATION: do not adjust p_plasma_loss_mw for radiation
        except ValueError as e:
            raise ProcessValueError(
                "Illegal value of i_rad_loss",
                i_rad_loss=physics_variables.i_rad_loss,
            ) from e

        # Ensure heating power is positive (shouldn't be necessary)
        p_plasma_loss_mw = max(p_plasma_loss_mw, 1.0e-3)

        # ========================================================================

        # Line averaged electron density in scaled units
        dnla20 = nd_plasma_electron_line * 1.0e-20
        dnla19 = nd_plasma_electron_line * 1.0e-19

        # Volume averaged electron density in units of 10**20 m**-3
        n20 = nd_plasma_electrons_vol_avg / 1.0e20

        # Plasma current in MA
        pcur = plasma_current / 1.0e6

        # Separatrix kappa defined with plasma volume for IPB scalings
        # Updated version of kappa used by the IPB98 scalings correction in:

        # None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy,
        # “Corrections to a sequence of papers in Nuclear Fusion,” Nuclear Fusion,
        # vol. 48, no. 9, pp. 099801099801, Aug. 2008,
        # doi: https://doi.org/10.1088/0029-5515/48/9/099801.

        physics_variables.kappa_ipb = (vol_plasma / (2.0 * np.pi * rmajor)) / (
            np.pi * rminor**2
        )

        # Electron energy confinement times

        try:
            model = ConfinementTimeModel(i_confinement_time)
        except ValueError as e:
            raise ProcessValueError(
                "Illegal value for i_confinement_time",
                i_confinement_time=i_confinement_time,
            ) from e

        # ========================================================================

        # User defined confinement time
        if (
            model == ConfinementTimeModel.USER_INPUT
        ):  # t_electron_energy_confinement is an input
            t_electron_confinement = physics_variables.tauee_in

        # ========================================================================

        # Nec-Alcator(NA) OH scaling
        if (
            model == ConfinementTimeModel.NEO_ALCATOR
        ):  # t_electron_energy_confinement is an input
            t_electron_confinement = self.neo_alcator_confinement_time(
                n20, rminor, rmajor, qstar
            )

        # ========================================================================

        # "Mirnov"-like scaling (H-mode)
        elif model == ConfinementTimeModel.MIRNOV:  # Mirnov scaling (H-mode)
            t_electron_confinement = self.mirnov_confinement_time(rminor, kappa95, pcur)

        # ========================================================================

        # Merezhkin-Mukhovatov (MM) OH/L-mode scaling
        elif model == ConfinementTimeModel.MEREZHKIN_MUHKOVATOV:
            t_electron_confinement = self.merezhkin_muhkovatov_confinement_time(
                rmajor,
                rminor,
                kappa95,
                qstar,
                dnla20,
                m_fuel_amu,
                temp_plasma_electron_density_weighted_kev,
            )

        # ========================================================================

        # Shimomura (S) optimized H-mode scaling
        elif model == ConfinementTimeModel.SHIMOMURA:
            t_electron_confinement = self.shimomura_confinement_time(
                rmajor, rminor, b_plasma_toroidal_on_axis, kappa95, m_fuel_amu
            )

        # ========================================================================

        # Kaye-Goldston scaling (L-mode)
        elif model == ConfinementTimeModel.KAYE_GOLDSTON:
            t_electron_confinement = self.kaye_goldston_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER Power scaling - ITER 89-P (L-mode)
        elif model == ConfinementTimeModel.ITER_89P:
            t_electron_confinement = self.iter_89p_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER Offset linear scaling - ITER 89-O (L-mode)
        elif model == ConfinementTimeModel.ITER_89_0:
            t_electron_confinement = self.iter_89_0_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )
        # ========================================================================

        # Rebut-Lallia offset linear scaling (L-mode)
        elif model == ConfinementTimeModel.REBUT_LALLIA:
            t_electron_confinement = self.rebut_lallia_confinement_time(
                rminor,
                rmajor,
                kappa,
                m_fuel_amu,
                pcur,
                zeff,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Goldston scaling (L-mode)
        elif model == ConfinementTimeModel.GOLDSTON:  # Goldston scaling (L-mode)
            t_electron_confinement = self.goldston_confinement_time(
                pcur, rmajor, rminor, kappa95, m_fuel_amu, p_plasma_loss_mw
            )

        # ========================================================================

        # T-10 scaling (L-mode)
        elif model == ConfinementTimeModel.T_10:
            t_electron_confinement = self.t10_confinement_time(
                dnla20,
                rmajor,
                qstar,
                b_plasma_toroidal_on_axis,
                rminor,
                kappa95,
                p_plasma_loss_mw,
                zeff,
                pcur,
            )

        # ========================================================================

        # JAERI / Odajima-Shimomura L-mode scaling
        elif model == ConfinementTimeModel.JAERI:  # JAERI scaling
            t_electron_confinement = self.jaeri_confinement_time(
                kappa95,
                rminor,
                m_fuel_amu,
                n20,
                pcur,
                b_plasma_toroidal_on_axis,
                rmajor,
                qstar,
                zeff,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Kaye "big"  L-mode scaling (based only on big tokamak data)
        elif model == ConfinementTimeModel.KAYE_BIG:
            t_electron_confinement = self.kaye_big_confinement_time(
                rmajor,
                rminor,
                b_plasma_toroidal_on_axis,
                kappa95,
                pcur,
                n20,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # ITER H90-P H-mode scaling
        elif model == ConfinementTimeModel.ITER_H90_P:
            t_electron_confinement = self.iter_h90_p_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Minimum of ITER 89-P and ITER 89-O
        elif model == ConfinementTimeModel.MINIMUM_OF_ITER_89P_AND_ITER_89_0:
            t_electron_confinement = min(
                self.iter_89p_confinement_time(
                    pcur,
                    rmajor,
                    rminor,
                    kappa,
                    dnla20,
                    b_plasma_toroidal_on_axis,
                    m_fuel_amu,
                    p_plasma_loss_mw,
                ),
                self.iter_89_0_confinement_time(
                    pcur,
                    rmajor,
                    rminor,
                    kappa,
                    dnla20,
                    b_plasma_toroidal_on_axis,
                    m_fuel_amu,
                    p_plasma_loss_mw,
                ),
            )

        # ========================================================================

        # Riedel scaling (L-mode)
        elif model == ConfinementTimeModel.RIEDEL_L:
            t_electron_confinement = self.riedel_l_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Christiansen et al scaling (L-mode)
        elif model == ConfinementTimeModel.CHRISTIANSEN:
            t_electron_confinement = self.christiansen_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                m_fuel_amu,
            )

        # ========================================================================

        # Lackner-Gottardi scaling (L-mode)
        elif model == ConfinementTimeModel.LACKNER_GOTTARDI:
            t_electron_confinement = self.lackner_gottardi_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Neo-Kaye scaling (L-mode)
        elif model == ConfinementTimeModel.NEO_KAYE:
            t_electron_confinement = self.neo_kaye_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ======== ================================================================

        # Riedel scaling (H-mode)
        elif model == ConfinementTimeModel.RIEDEL_H:
            t_electron_confinement = self.riedel_h_confinement_time(
                pcur,
                rmajor,
                rminor,
                kappa95,
                dnla20,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ========================================================================

        # Amended version of ITER H90-P law
        elif model == ConfinementTimeModel.ITER_H90_P_AMENDED:
            t_electron_confinement = self.iter_h90_p_amended_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                m_fuel_amu,
                rmajor,
                p_plasma_loss_mw,
                kappa,
            )

        # ==========================================================================

        # Sudo et al. scaling (stellarators/heliotron)
        elif model == ConfinementTimeModel.SUDO_ET_AL:
            t_electron_confinement = self.sudo_et_al_confinement_time(
                rmajor,
                rminor,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Gyro-reduced Bohm scaling
        elif model == ConfinementTimeModel.GYRO_REDUCED_BOHM:
            t_electron_confinement = self.gyro_reduced_bohm_confinement_time(
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
                rminor,
                rmajor,
            )

        # ==========================================================================

        # Lackner-Gottardi stellarator scaling
        elif model == ConfinementTimeModel.LACKNER_GOTTARDI_STELLARATOR:
            t_electron_confinement = self.lackner_gottardi_stellarator_confinement_time(
                rmajor,
                rminor,
                dnla20,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                q95,
            )

        # ==========================================================================

        # ITER_93 ELM-free H-mode scaling
        elif model == ConfinementTimeModel.ITER_93H:
            t_electron_confinement = self.iter_93h_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                m_fuel_amu,
                rmajor,
                dnla20,
                aspect,
                kappa,
            )

        # ==========================================================================
        # Scaling removed
        elif model == ConfinementTimeModel.TITAN_REMOVED:
            raise ProcessValueError("Scaling removed")
        # ==========================================================================

        # ELM-free: ITERH-97P
        elif model == ConfinementTimeModel.ITER_H97P:
            t_electron_confinement = self.iter_h97p_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                rmajor,
                aspect,
                kappa,
                m_fuel_amu,
            )

        # ==========================================================================

        # ELMy: ITERH-97P(y)
        elif model == ConfinementTimeModel.ITER_H97P_ELMY:
            t_electron_confinement = self.iter_h97p_elmy_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                rmajor,
                aspect,
                kappa,
                m_fuel_amu,
            )

        # ==========================================================================

        # ITER-96P (= ITER-97L) L-mode scaling
        elif model == ConfinementTimeModel.ITER_96P:
            t_electron_confinement = self.iter_96p_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                kappa95,
                rmajor,
                aspect,
                dnla19,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Valovic modified ELMy-H mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif model == ConfinementTimeModel.VALOVIC_ELMY:
            t_electron_confinement = self.valovic_elmy_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                m_fuel_amu,
                rmajor,
                rminor,
                kappa,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Kaye PPPL Workshop April 1998 L-mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif model == ConfinementTimeModel.KAYE:
            t_electron_confinement = self.kaye_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                kappa,
                rmajor,
                aspect,
                dnla19,
                m_fuel_amu,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # ITERH-PB98P(y), ELMy H-mode scaling
        # WARNING: No reference found for this scaling. This may not be its real name
        elif model == ConfinementTimeModel.ITER_PB98P_Y:
            t_electron_confinement = self.iter_pb98py_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y), ELMy H-mode scaling
        elif model == ConfinementTimeModel.IPB98_Y:
            t_electron_confinement = self.iter_ipb98y_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,1), ELMy H-mode scaling
        elif model == ConfinementTimeModel.ITER_IPB98Y1:
            t_electron_confinement = self.iter_ipb98y1_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,2), ELMy H-mode scaling
        elif model == ConfinementTimeModel.ITER_IPB98Y2:
            t_electron_confinement = self.iter_ipb98y2_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,3), ELMy H-mode scaling
        elif model == ConfinementTimeModel.ITER_IPB98Y3:
            t_electron_confinement = self.iter_ipb98y3_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # IPB98(y,4), ELMy H-mode scaling
        elif model == ConfinementTimeModel.ITER_IPB98Y4:
            t_electron_confinement = self.iter_ipb98y4_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # ISS95 stellarator scaling
        elif model == ConfinementTimeModel.ISS95_STELLARATOR:
            # dummy argument q95 is actual argument iotabar for stellarators
            iotabar = q95
            t_electron_confinement = self.iss95_stellarator_confinement_time(
                rminor,
                rmajor,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                iotabar,
            )

        # ==========================================================================

        # ISS04 stellarator scaling
        elif model == ConfinementTimeModel.ISS04_STELLARATOR:
            # dummy argument q95 is actual argument iotabar for stellarators
            iotabar = q95
            t_electron_confinement = self.iss04_stellarator_confinement_time(
                rminor,
                rmajor,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                iotabar,
            )

        # ==========================================================================

        # DS03 beta-independent H-mode scaling
        elif model == ConfinementTimeModel.DS03:
            t_electron_confinement = self.ds03_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa95,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        #  Murari "Non-power law" scaling
        elif model == ConfinementTimeModel.MURARI:
            t_electron_confinement = self.murari_confinement_time(
                pcur,
                rmajor,
                physics_variables.kappa_ipb,
                dnla19,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Petty08, beta independent dimensionless scaling
        elif model == ConfinementTimeModel.PETTY08:
            t_electron_confinement = self.petty08_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
            )

        # ==========================================================================

        # Lang high density relevant confinement scaling
        elif model == ConfinementTimeModel.LANG_HIGH_DENSITY:
            t_electron_confinement = self.lang_high_density_confinement_time(
                plasma_current,
                b_plasma_toroidal_on_axis,
                nd_plasma_electron_line,
                p_plasma_loss_mw,
                rmajor,
                rminor,
                q95,
                qstar,
                aspect,
                m_fuel_amu,
                physics_variables.kappa_ipb,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - nominal
        elif model == ConfinementTimeModel.HUBBARD_NOMINAL:
            t_electron_confinement = self.hubbard_nominal_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - lower
        elif model == ConfinementTimeModel.HUBBARD_LOWER:
            t_electron_confinement = self.hubbard_lower_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Hubbard 2017 I-mode confinement time scaling - upper
        elif model == ConfinementTimeModel.HUBBARD_UPPER:
            t_electron_confinement = self.hubbard_upper_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla20,
                p_plasma_loss_mw,
            )

        # ==========================================================================

        # Menard NSTX, ELMy H-mode scaling
        elif model == ConfinementTimeModel.MENARD_NSTX:
            t_electron_confinement = self.menard_nstx_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # Menard NSTX-Petty08 Hybrid
        elif model == ConfinementTimeModel.MENARD_NSTX_PETTY08_HYBRID:
            t_electron_confinement = self.menard_nstx_petty08_hybrid_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.kappa_ipb,
                aspect,
                m_fuel_amu,
            )

        # ==========================================================================

        # NSTX gyro-Bohm (Buxton)
        elif model == ConfinementTimeModel.NSTX_GYRO_BOHM:
            t_electron_confinement = self.nstx_gyro_bohm_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                rmajor,
                dnla20,
            )

        # ==========================================================================

        # ITPA20 H-mode scaling
        elif model == ConfinementTimeModel.ITPA20:
            t_electron_confinement = self.itpa20_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                physics_variables.triang,
                physics_variables.kappa_ipb,
                eps,
                physics_variables.m_ions_total_amu,
            )

        # ==========================================================================

        # ITPA20-IL confinement time scaling
        elif model == ConfinementTimeModel.ITPA20_IL:
            t_electron_confinement = self.itpa20_il_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                p_plasma_loss_mw,
                dnla19,
                physics_variables.m_ions_total_amu,
                rmajor,
                physics_variables.triang,
                physics_variables.kappa_ipb,
            )

        # ==========================================================================

        else:
            raise ProcessValueError(
                "Illegal value for i_confinement_time",
                i_confinement_time=i_confinement_time,
            )

        # Apply H-factor correction to chosen scaling
        t_electron_energy_confinement = hfact * t_electron_confinement

        # Ion energy confinement time
        t_ion_energy_confinement = t_electron_energy_confinement

        # Calculate H* non-radiation corrected H factor
        # Note: we will assume the IPB-98y2 scaling.
        if physics_variables.i_rad_loss == 1:
            physics_variables.hstar = (
                hfact
                * (
                    p_plasma_loss_mw
                    / (
                        p_plasma_loss_mw
                        + physics_variables.pden_plasma_sync_mw * vol_plasma
                        + physics_variables.p_plasma_inner_rad_mw
                    )
                )
                ** 0.31
            )
        elif physics_variables.i_rad_loss == 0:
            physics_variables.hstar = (
                hfact
                * (
                    p_plasma_loss_mw
                    / (
                        p_plasma_loss_mw
                        + physics_variables.pden_plasma_rad_mw * vol_plasma
                    )
                )
                ** 0.31
            )
        elif physics_variables.i_rad_loss == 2:
            physics_variables.hstar = hfact

        # Calculation of the transport power loss terms
        # Transport losses in Watts/m3 are 3/2 * n.e.T / tau , with T in eV
        # (here, temp_plasma_ion_density_weighted_kev and temp_plasma_electron_density_weighted_kev are in keV, and pden_electron_transport_loss_mw and pden_ion_transport_loss_mw are in MW/m3)

        # The transport losses is just the electron and ion thermal energies divided by the confinement time.
        pden_ion_transport_loss_mw = (
            (3 / 2)
            * (constants.ELECTRON_CHARGE / 1e3)
            * nd_plasma_ions_total_vol_avg
            * temp_plasma_ion_density_weighted_kev
            / t_ion_energy_confinement
        )
        pden_electron_transport_loss_mw = (
            (3 / 2)
            * (constants.ELECTRON_CHARGE / 1e3)
            * nd_plasma_electrons_vol_avg
            * temp_plasma_electron_density_weighted_kev
            / t_electron_energy_confinement
        )

        ratio = (nd_plasma_ions_total_vol_avg / nd_plasma_electrons_vol_avg) * (
            temp_plasma_ion_density_weighted_kev
            / temp_plasma_electron_density_weighted_kev
        )

        # Global energy confinement time

        t_energy_confinement = (ratio + 1.0e0) / (
            ratio / t_ion_energy_confinement + 1.0e0 / t_electron_energy_confinement
        )

        # For comparison directly calculate the confinement time from the stored energy calculated
        # from the total plasma beta and the loss power used above.
        physics_variables.t_energy_confinement_beta = (
            physics_variables.e_plasma_beta / 1e6
        ) / p_plasma_loss_mw

        return (
            pden_electron_transport_loss_mw,
            pden_ion_transport_loss_mw,
            t_electron_energy_confinement,
            t_ion_energy_confinement,
            t_energy_confinement,
            p_plasma_loss_mw,
        )

    @staticmethod
    def calculate_double_and_triple_product(
        nd_plasma_electrons_vol_avg: float,
        temp_plasma_electrons_vol_avg_kev: float,
        t_energy_confinement: float,
    ) -> tuple[float, float]:
        """Calculate the plasma double (nτ_E) and triple product (nTτ_E)

        Parameters
        ----------
        nd_plasma_electrons_vol_avg :
            Volume averaged electron density [m⁻³]
        temp_plasma_electrons_vol_avg_kev :
            Volume averaged electron temperature [keV]
        t_energy_confinement :
            Energy confinement time [s]

        Returns
        -------
        :
            Tuple[float, float]: (ntau, nTtau) where ntau is the plasma double product (n * τ_E) in units of m⁻³ * s,
            and nTtau is the plasma triple product (n * T * τ_E) in units of keV * s * m⁻³

        """
        ntau = t_energy_confinement * nd_plasma_electrons_vol_avg
        nTtau = ntau * temp_plasma_electrons_vol_avg_kev

        return ntau, nTtau

    def find_other_h_factors(self, i_confinement_time: int) -> float:
        """Function to find H-factor for the equivalent confinement time in other scalings.

        Parameters
        ----------
        i_confinement_time : int
            Index of the confinement time scaling to use.

        Returns
        -------
        float
            The calculated H-factor.

        """

        def fhz(hfact: float) -> float:
            """Function used to find power balance.

            Parameters
            ----------
            hfact : float
                H-factor to be used in the calculation.
            hfact: float :


            Returns
            -------
            float
                The difference between the calculated power and the required power for balance.

            """
            (
                ptrez,
                ptriz,
                _,
                _,
                _,
                _,
            ) = self.calculate_confinement_time(
                m_fuel_amu=physics_variables.m_fuel_amu,
                p_alpha_total_mw=physics_variables.p_alpha_total_mw,
                aspect=physics_variables.aspect,
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                nd_plasma_ions_total_vol_avg=physics_variables.nd_plasma_ions_total_vol_avg,
                nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                nd_plasma_electron_line=physics_variables.nd_plasma_electron_line,
                eps=physics_variables.eps,
                hfact=hfact,
                i_confinement_time=i_confinement_time,
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

            # At power balance, fhz is zero.
            fhz_value = (
                ptrez
                + ptriz
                - physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_non_alpha_charged_mw
                - physics_variables.pden_plasma_ohmic_mw
            )

            # Take into account whether injected power is included in tau_e calculation (i.e. whether device is ignited)
            if physics_variables.i_plasma_ignited == 0:
                fhz_value -= (
                    current_drive_variables.p_hcd_injected_total_mw
                    / physics_variables.vol_plasma
                )

            # Include the radiation power if requested
            if physics_variables.i_rad_loss == 0:
                fhz_value += physics_variables.pden_plasma_rad_mw
            elif physics_variables.i_rad_loss == 1:
                fhz_value += physics_variables.pden_plasma_core_rad_mw

            return fhz_value

        return root_scalar(fhz, bracket=(0.01, 150), xtol=0.001).root

    def output_confinement_time_info(self):
        po.oheadr(self.outfile, "Energy Confinement")

        if physics_variables.i_plasma_ignited == 1:
            po.ocmmnt(
                self.outfile,
                "Device is assumed to be ignited for the calculation of confinement time",
            )
            po.oblnkl(self.outfile)

        tauelaw = ConfinementTimeModel(physics_variables.i_confinement_time).full_name

        po.ocmmnt(
            self.outfile,
            f"Confinement scaling law: {tauelaw}",
        )

        po.ovarst(
            self.outfile,
            "Confinement scaling law",
            "(tauelaw)",
            f'"{tauelaw.strip().split(" ")[0]}"',
        )

        po.ovarrf(
            self.outfile, "Confinement H factor", "(hfact)", physics_variables.hfact
        )
        po.ovarrf(
            self.outfile,
            "Global thermal energy confinement time, from scaling (s)",
            "(t_energy_confinement)",
            physics_variables.t_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Directly calculated total energy confinement time (s)",
            "(t_energy_confinement_beta)",
            physics_variables.t_energy_confinement_beta,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "(Total thermal energy derived from total plasma beta / loss power)",
        )
        po.ovarrf(
            self.outfile,
            "Ion energy confinement time, from scaling (s)",
            "(t_ion_energy_confinement)",
            physics_variables.t_ion_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron energy confinement time, from scaling (s)",
            "(t_electron_energy_confinement)",
            physics_variables.t_electron_energy_confinement,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fusion double product (s/m3)",
            "(ntau)",
            physics_variables.ntau,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Lawson Triple product (keV s/m3)",
            "(nTtau)",
            physics_variables.nTtau,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Transport loss power assumed in scaling law (MW)",
            "(p_plasma_loss_mw)",
            physics_variables.p_plasma_loss_mw,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "Switch for radiation loss term usage in power balance",
            "(i_rad_loss)",
            physics_variables.i_rad_loss,
        )
        model = ConfinementRadiationLossModel(physics_variables.i_rad_loss)

        if model == ConfinementRadiationLossModel.FULL_RADIATION:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma heating power balance (MW)",
                "",
                physics_variables.p_plasma_rad_mw,
                "OP ",
            )
        elif model == ConfinementRadiationLossModel.CORE_ONLY:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma heating power balance (MW)",
                "",
                physics_variables.p_plasma_inner_rad_mw,
                "OP ",
            )
        else:  # NO_RADIATION
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma heating power balance (MW)",
                "",
                0.0e0,
            )

        po.ocmmnt(self.outfile, f"  (Radiation correction: {model.description})")
        po.ovarrf(
            self.outfile,
            "H* non-radiation corrected",
            "(hstar)",
            physics_variables.hstar,
            "OP",
        )
        po.ocmmnt(self.outfile, "  (H* assumes IPB98(y,2), ELMy H-mode scaling)")
        po.ovarrf(
            self.outfile,
            "Alpha particle confinement time (s)",
            "(t_alpha_confinement)",
            physics_variables.t_alpha_confinement,
            "OP ",
        )
        # Note alpha confinement time is no longer equal to fuel particle confinement time.
        po.ovarrf(
            self.outfile,
            "Alpha particle/energy confinement time ratio",
            "(f_alpha_energy_confinement)",
            physics_variables.f_alpha_energy_confinement,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Lower limit on f_alpha_energy_confinement",
            "(f_alpha_energy_confinement_min)",
            constraint_variables.f_alpha_energy_confinement_min,
        )

        # Plot table of al the H-factor scalings and coparison values
        self.output_confinement_comparison(istell=stellarator_variables.istell)

    def output_confinement_comparison(self, istell: int):
        """Routine to calculate ignition margin for different confinement scalings and equivalent confinement times for H=1.

        This routine calculates the ignition margin at the final point with different scalings and outputs the results to a file.

        The output includes:
        - Energy confinement times
        - Required H-factors for power balance

        The routine iterates over a range of confinement times, skipping the first user input and a specific index (25). For each confinement time, it calculates various parameters related to confinement and ignition using the `calculate_confinement_time` method. It then calculates the H-factor for when the plasma is ignited using the `find_other_h_factors` method and writes the results to the output file.

        Output format:
        - Header: "Energy confinement times, and required H-factors :"
        - Columns: "Scaling law", "Confinement time [s]", "H-factor for power balance"

        Methods used:
        - `calculate_confinement_time`: Calculates confinement-related parameters.
        - `find_other_h_factors`: Calculates the H-factor for a given confinement time.

        Parameters
        ----------
        istell :
            Indicator for stellarator (0 for tokamak, >=1 for stellarator).

        """
        po.oheadr(self.outfile, "Energy confinement times, and required H-factors :")
        po.ocmmnt(
            self.outfile,
            f"{'':>2}{'Scaling law':<27}{'Electron confinement time [s]':<32}Equivalent H-factor for",
        )
        po.ocmmnt(
            self.outfile,
            f"{'':>38}{'for H = 1':<23}same confinement time",
        )
        po.oblnkl(self.outfile)

        # List of key values for stellarator scalings
        stellarator_scalings = [21, 22, 23, 37, 38]

        # Plot all of the confinement scalings for comparison when H = 1
        # Start from range 1 as the first i_confinement_time is a user input
        # If stellarator, use the stellarator scalings
        for i_confinement_time in (
            range(1, physics_variables.N_CONFINEMENT_SCALINGS)
            if istell == 0
            else stellarator_scalings
        ):
            if i_confinement_time == 25:
                continue
            (
                _,
                _,
                taueez,
                _,
                _,
                _,
            ) = self.calculate_confinement_time(
                m_fuel_amu=physics_variables.m_fuel_amu,
                p_alpha_total_mw=physics_variables.p_alpha_total_mw,
                aspect=physics_variables.aspect,
                b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                nd_plasma_ions_total_vol_avg=physics_variables.nd_plasma_ions_total_vol_avg,
                nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                nd_plasma_electron_line=physics_variables.nd_plasma_electron_line,
                eps=physics_variables.eps,
                hfact=1.0,
                i_confinement_time=i_confinement_time,
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

            try:
                # Calculate the H-factor for the same confinement time in other scalings
                physics_variables.hfac[i_confinement_time - 1] = (
                    self.find_other_h_factors(i_confinement_time)
                )
            except ValueError:
                # This is only used for a table in the OUT.DAT so if it fails
                # just write a NaN--its not worth crashing PROCESS over.
                physics_variables.hfac[i_confinement_time - 1] = np.nan

            scaling_name = ConfinementTimeModel(i_confinement_time).full_name

            po.ocmmnt(
                self.outfile,
                f"{'':>2}{scaling_name:<38}"
                f"{taueez:<28.3f}{physics_variables.hfac[i_confinement_time - 1]:.3f}",
            )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)

    @staticmethod
    def neo_alcator_confinement_time(
        dene20: float, rminor: float, rmajor: float, qstar: float
    ) -> float:
        """Calculate the Nec-Alcator(NA) OH scaling confinement time

        Parameters
        ----------
        dene20 :
            Volume averaged electron density in units of 10**20 m**-3
        rminor :
            Plasma minor radius [m]
        rmajor :
            Plasma major radius [m]
        qstar :
            Equivalent cylindrical edge safety factor

        Returns
        -------
        :
            float: Neo-Alcator confinement time [s]


        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return 0.07e0 * dene20 * rminor * rmajor * rmajor * qstar

    @staticmethod
    def mirnov_confinement_time(rminor: float, kappa95: float, pcur: float) -> float:
        """Calculate the Mirnov scaling (H-mode) confinement time

        Parameters
        ----------
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        pcur :
            Plasma current [MA]

        Returns
        -------
        :
            float: Mirnov scaling confinement time [s]

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return 0.2e0 * rminor * np.sqrt(kappa95) * pcur

    @staticmethod
    def merezhkin_muhkovatov_confinement_time(
        rmajor: float,
        rminor: float,
        kappa95: float,
        qstar: float,
        dnla20: float,
        afuel: float,
        ten: float,
    ) -> float:
        """Calculate the Merezhkin-Mukhovatov (MM) OH/L-mode scaling confinement time

        Parameters
        ----------
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        qstar :
            Equivalent cylindrical edge safety factor
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        afuel :
            Fuel atomic mass number
        ten :
            Electron temperature [keV]

        Returns
        -------
        :
            float: Merezhkin-Mukhovatov confinement time [s]


        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            3.5e-3
            * rmajor**2.75e0
            * rminor**0.25e0
            * kappa95**0.125e0
            * qstar
            * dnla20
            * np.sqrt(afuel)
            / np.sqrt(ten / 10.0e0)
        )

    @staticmethod
    def shimomura_confinement_time(
        rmajor: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        kappa95: float,
        afuel: float,
    ) -> float:
        """Calculate the  Shimomura (S) optimized H-mode scaling confinement time

        Parameters
        ----------
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        kappa95 :
            Plasma elongation at 95% flux surface
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: Shimomura confinement time [s]

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            0.045e0
            * rmajor
            * rminor
            * b_plasma_toroidal_on_axis
            * np.sqrt(kappa95)
            * np.sqrt(afuel)
        )

    @staticmethod
    def kaye_goldston_confinement_time(
        kappa95: float,
        pcur: float,
        n20: float,
        rmajor: float,
        afuel: float,
        b_plasma_toroidal_on_axis: float,
        rminor: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Kaye-Goldston (KG) L-mode scaling confinement time

        Parameters
        ----------
        kappa95 :
            Plasma elongation at 95% flux surface
        pcur :
            Plasma current [MA]
        n20 :
            Line averaged electron density in units of 10**20 m**-3
        rmajor :
            Plasma major radius [m]
        afuel :
            Fuel atomic mass number
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        rminor :
            Plasma minor radius [m]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Kaye-Goldston confinement time [s]

        Notes
        -----
            - An isotope correction factor (M_i/1.5)^0.5 is added to the original scaling to reflect the fact
            that the empirical fits to the data were from experiments with H and D mixture, M_i = 1.5

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            0.055e0
            * kappa95**0.28e0
            * pcur**1.24e0
            * n20**0.26e0
            * rmajor**1.65e0
            * np.sqrt(afuel / 1.5e0)
            / (
                b_plasma_toroidal_on_axis**0.09e0
                * rminor**0.49e0
                * p_plasma_loss_mw**0.58e0
            )
        )

    @staticmethod
    def iter_89p_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the ITER Power scaling - ITER 89-P (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa :
            Plasma elongation
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: ITER 89-P confinement time [s]


        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            0.048e0
            * pcur**0.85e0
            * rmajor**1.2e0
            * rminor**0.3e0
            * np.sqrt(kappa)
            * dnla20**0.1e0
            * b_plasma_toroidal_on_axis**0.2e0
            * np.sqrt(afuel)
            / np.sqrt(p_plasma_loss_mw)
        )

    @staticmethod
    def iter_89_0_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the ITER Offset linear scaling - ITER 89-O (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa :
            Plasma elongation
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: ITER 89-O confinement time [s]

        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        term1 = (
            0.04e0
            * pcur**0.5e0
            * rmajor**0.3e0
            * rminor**0.8e0
            * kappa**0.6e0
            * afuel**0.5e0
        )
        term2 = (
            0.064e0
            * pcur**0.8e0
            * rmajor**1.6e0
            * rminor**0.6e0
            * kappa**0.5e0
            * dnla20**0.6e0
            * b_plasma_toroidal_on_axis**0.35e0
            * afuel**0.2e0
            / p_plasma_loss_mw
        )
        return term1 + term2

    @staticmethod
    def rebut_lallia_confinement_time(
        rminor: float,
        rmajor: float,
        kappa: float,
        afuel: float,
        pcur: float,
        zeff: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Rebut-Lallia offset linear scaling (L-mode) confinement time

        Parameters
        ----------
        rminor :
            Plasma minor radius [m]
        rmajor :
            Plasma major radius [m]
        kappa :
            Plasma elongation at 95% flux surface
        afuel :
            Fuel atomic mass number
        pcur :
            Plasma current [MA]
        zeff :
            Effective charge
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Rebut-Lallia confinement time [s]


        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        rll = (rminor**2 * rmajor * kappa) ** (1.0e0 / 3.0e0)
        term1 = 1.2e-2 * pcur * rll**1.5e0 / np.sqrt(zeff)
        term2 = (
            0.146e0
            * dnla20**0.75e0
            * np.sqrt(pcur)
            * np.sqrt(b_plasma_toroidal_on_axis)
            * rll**2.75e0
            * zeff**0.25e0
            / p_plasma_loss_mw
        )
        return 1.65e0 * np.sqrt(afuel / 2.0e0) * (term1 + term2)

    @staticmethod
    def goldston_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Goldston scaling (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Goldston confinement time [s]

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            0.037e0
            * pcur
            * rmajor**1.75e0
            * rminor ** (-0.37e0)
            * np.sqrt(kappa95)
            * np.sqrt(afuel / 1.5e0)
            / np.sqrt(p_plasma_loss_mw)
        )

    @staticmethod
    def t10_confinement_time(
        dnla20: float,
        rmajor: float,
        qstar: float,
        b_plasma_toroidal_on_axis: float,
        rminor: float,
        kappa95: float,
        p_plasma_loss_mw: float,
        zeff: float,
        pcur: float,
    ) -> float:
        """Calculate the T-10 scaling confinement time

        Parameters
        ----------
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        rmajor :
            Plasma major radius [m]
        qstar :
            Equivalent cylindrical edge safety factor
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        p_plasma_loss_mw :
            Net Heating power [MW]
        zeff :
            Effective charge
        pcur :
            Plasma current [MA]

        Returns
        -------
        :
            float: T-10 confinement time [s]

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        denfac = dnla20 * rmajor * qstar / (1.3e0 * b_plasma_toroidal_on_axis)
        denfac = min(1.0e0, denfac)
        return (
            0.095e0
            * rmajor
            * rminor
            * b_plasma_toroidal_on_axis
            * np.sqrt(kappa95)
            * denfac
            / p_plasma_loss_mw**0.4e0
            * (zeff**2 * pcur**4 / (rmajor * rminor * qstar**3 * kappa95**1.5e0))
            ** 0.08e0
        )

    @staticmethod
    def jaeri_confinement_time(
        kappa95: float,
        rminor: float,
        afuel: float,
        n20: float,
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        qstar: float,
        zeff: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the JAERI / Odajima-Shimomura L-mode scaling confinement time

        Parameters
        ----------
        kappa95 :
            Plasma elongation at 95% flux surface
        rminor :
            Plasma minor radius [m]
        afuel :
            Fuel atomic mass number
        n20 :
            Line averaged electron density in units of 10**20 m**-3
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        rmajor :
            Plasma major radius [m]
        qstar :
            Equivalent cylindrical edge safety factor
        zeff :
            Effective charge
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: JAERI confinement time [s]

        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        gjaeri = (
            zeff**0.4e0
            * ((15.0e0 - zeff) / 20.0e0) ** 0.6e0
            * (3.0e0 * qstar * (qstar + 5.0e0) / ((qstar + 2.0e0) * (qstar + 7.0e0)))
            ** 0.6e0
        )
        return (
            0.085e0 * kappa95 * rminor**2 * np.sqrt(afuel)
            + 0.069e0
            * n20**0.6e0
            * pcur
            * b_plasma_toroidal_on_axis**0.2e0
            * rminor**0.4e0
            * rmajor**1.6e0
            * np.sqrt(afuel)
            * gjaeri
            * kappa95**0.2e0
            / p_plasma_loss_mw
        )

    @staticmethod
    def kaye_big_confinement_time(
        rmajor: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        kappa95: float,
        pcur: float,
        n20: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Kaye-Big scaling confinement time

        Parameters
        ----------
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        kappa95 :
            Plasma elongation at 95% flux surface
        pcur :
            Plasma current [MA]
        n20 :
            Line averaged electron density in units of 10**20 m**-3
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Kaye-Big confinement time [s]


        References
        ----------
            - N. A. Uckan, International Atomic Energy Agency, Vienna (Austria)and ITER Physics Group,
            "ITER physics design guidelines: 1989", no. No. 10. Feb. 1990.
        """
        return (
            0.105e0
            * np.sqrt(rmajor)
            * rminor**0.8e0
            * b_plasma_toroidal_on_axis**0.3e0
            * kappa95**0.25e0
            * pcur**0.85e0
            * n20**0.1e0
            * np.sqrt(afuel)
            / np.sqrt(p_plasma_loss_mw)
        )

    @staticmethod
    def iter_h90_p_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the ITER H-mode scaling - ITER H90-P confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa :
            Plasma elongation
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: ITER H90-P confinement time [s]


        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.064e0
            * pcur**0.87e0
            * rmajor**1.82e0
            * rminor ** (-0.12e0)
            * kappa**0.35e0
            * dnla20**0.09e0
            * b_plasma_toroidal_on_axis**0.15e0
            * np.sqrt(afuel)
            / np.sqrt(p_plasma_loss_mw)
        )

    @staticmethod
    def riedel_l_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Riedel scaling (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Riedel confinement time [s]

        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.044e0
            * pcur**0.93e0
            * rmajor**1.37e0
            * rminor ** (-0.049e0)
            * kappa95**0.588e0
            * dnla20**0.078e0
            * b_plasma_toroidal_on_axis**0.152e0
            / p_plasma_loss_mw**0.537e0
        )

    @staticmethod
    def christiansen_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        afuel: float,
    ) -> float:
        """Calculate the Christiansen et al scaling (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: Christiansen confinement time [s]

        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.24e0
            * pcur**0.79e0
            * rmajor**0.56e0
            * rminor**1.46e0
            * kappa95**0.73e0
            * dnla20**0.41e0
            * b_plasma_toroidal_on_axis**0.29e0
            / (p_plasma_loss_mw**0.79e0 * afuel**0.02e0)
        )

    @staticmethod
    def lackner_gottardi_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Lackner-Gottardi scaling (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Lackner-Gottardi confinement time [s]

        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        qhat = (
            (1.0e0 + kappa95**2)
            * rminor**2
            * b_plasma_toroidal_on_axis
            / (0.4e0 * pcur * rmajor)
        )
        return (
            0.12e0
            * pcur**0.8e0
            * rmajor**1.8e0
            * rminor**0.4e0
            * kappa95
            * (1.0e0 + kappa95) ** (-0.8e0)
            * dnla20**0.6e0
            * qhat**0.4e0
            / p_plasma_loss_mw**0.6e0
        )

    @staticmethod
    def neo_kaye_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Neo-Kaye scaling (L-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Neo-Kaye confinement time [s]


        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.063e0
            * pcur**1.12e0
            * rmajor**1.3e0
            * rminor ** (-0.04e0)
            * kappa95**0.28e0
            * dnla20**0.14e0
            * b_plasma_toroidal_on_axis**0.04e0
            / p_plasma_loss_mw**0.59e0
        )

    @staticmethod
    def riedel_h_confinement_time(
        pcur: float,
        rmajor: float,
        rminor: float,
        kappa95: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Riedel scaling (H-mode) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Riedel H-mode confinement time [s]

        References
        ----------
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.1e0
            * np.sqrt(afuel)
            * pcur**0.884e0
            * rmajor**1.24e0
            * rminor ** (-0.23e0)
            * kappa95**0.317e0
            * b_plasma_toroidal_on_axis**0.207e0
            * dnla20**0.105e0
            / p_plasma_loss_mw**0.486e0
        )

    @staticmethod
    def iter_h90_p_amended_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        afuel: float,
        rmajor: float,
        p_plasma_loss_mw: float,
        kappa: float,
    ) -> float:
        """Calculate the amended ITER H90-P confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        afuel :
            Fuel atomic mass number
        rmajor :
            Plasma major radius [m]
        p_plasma_loss_mw :
            Net Heating power [MW]
        kappa :
            Plasma elongation

        Returns
        -------
        :
            float: Amended ITER H90-P confinement time [s]

        References
        ----------
            - J. P. Christiansen et al., “Global energy confinement H-mode database for ITER,”
            Nuclear Fusion, vol. 32, no. 2, pp. 291-338, Feb. 1992,
            doi: https://doi.org/10.1088/0029-5515/32/2/i11.
        """
        return (
            0.082e0
            * pcur**1.02e0
            * b_plasma_toroidal_on_axis**0.15e0
            * np.sqrt(afuel)
            * rmajor**1.60e0
            / (p_plasma_loss_mw**0.47e0 * kappa**0.19e0)
        )

    @staticmethod
    def sudo_et_al_confinement_time(
        rmajor: float,
        rminor: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Sudo et al. scaling confinement time

        Parameters
        ----------
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Sudo et al. confinement time [s]

        References
        ----------
            - S. Sudo et al., “Scalings of energy confinement and density limit in stellarator/heliotron devices,”
            Nuclear Fusion, vol. 30, no. 1, pp. 11-21, Jan. 1990,
            doi: https://doi.org/10.1088/0029-5515/30/1/002.
        """
        return (
            0.17e0
            * rmajor**0.75e0
            * rminor**2
            * dnla20**0.69e0
            * b_plasma_toroidal_on_axis**0.84e0
            * p_plasma_loss_mw ** (-0.58e0)
        )

    @staticmethod
    def gyro_reduced_bohm_confinement_time(
        b_plasma_toroidal_on_axis: float,
        dnla20: float,
        p_plasma_loss_mw: float,
        rminor: float,
        rmajor: float,
    ) -> float:
        """Calculate the Gyro-reduced Bohm scaling confinement time

        Parameters
        ----------
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rminor :
            Plasma minor radius [m]
        rmajor :
            Plasma major radius [m]

        Returns
        -------
        :
            float: Gyro-reduced Bohm confinement time [s]

        References
        ----------
            - Goldston, R. J., H. Biglari, and G. W. Hammett. "E x B/B 2 vs. μ B/B as the Cause of Transport in Tokamaks."
            Bull. Am. Phys. Soc 34 (1989): 1964.
        """
        return (
            0.25e0
            * b_plasma_toroidal_on_axis**0.8e0
            * dnla20**0.6e0
            * p_plasma_loss_mw ** (-0.6e0)
            * rminor**2.4e0
            * rmajor**0.6e0
        )

    @staticmethod
    def lackner_gottardi_stellarator_confinement_time(
        rmajor: float,
        rminor: float,
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        q: float,
    ) -> float:
        """Calculate the Lackner-Gottardi stellarator scaling confinement time

        Parameters
        ----------
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        q :
            Edge safety factor

        Returns
        -------
        :
            float: Lackner-Gottardi stellarator confinement time [s]

        References
        ----------
            - K. Lackner and N. A. O. Gottardi, “Tokamak confinement in relation to plateau scaling,”
            Nuclear Fusion, vol. 30, no. 4, pp. 767-770, Apr. 1990,
            doi: https://doi.org/10.1088/0029-5515/30/4/018.
        """
        return (
            0.17e0
            * rmajor
            * rminor**2
            * dnla20**0.6e0
            * b_plasma_toroidal_on_axis**0.8e0
            * p_plasma_loss_mw ** (-0.6e0)
            * q**0.4e0
        )

    @staticmethod
    def iter_93h_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        afuel: float,
        rmajor: float,
        dnla20: float,
        aspect: float,
        kappa: float,
    ) -> float:
        """Calculate the ITER-93H scaling ELM-free confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        afuel :
            Fuel atomic mass number
        rmajor :
            Plasma major radius [m]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        aspect :
            Aspect ratio
        kappa :
            Plasma elongation

        Returns
        -------
        :
            float: ITER-93H confinement time [s]

        References
        ----------
            - K. Thomsen et al., “ITER H mode confinement database update,”
            vol. 34, no. 1, pp. 131-167, Jan. 1994, doi: https://doi.org/10.1088/0029-5515/34/1/i10.
        """
        return (
            0.036e0
            * pcur**1.06e0
            * b_plasma_toroidal_on_axis**0.32e0
            * p_plasma_loss_mw ** (-0.67e0)
            * afuel**0.41e0
            * rmajor**1.79e0
            * dnla20**0.17e0
            * aspect**0.11e0
            * kappa**0.66e0
        )

    @staticmethod
    def iter_h97p_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        dnla19: float,
        rmajor: float,
        aspect: float,
        kappa: float,
        afuel: float,
    ) -> float:
        """Calculate the ELM-free ITER H-mode scaling - ITER H97-P confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        rmajor :
            Plasma major radius [m]
        aspect :
            Aspect ratio
        kappa :
            Plasma elongation
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: ITER H97-P confinement time [s]

        References
        ----------
            - I. C. Database and M. W. G. (presented Cordey), “Energy confinement scaling and the extrapolation to ITER,”
            Plasma Physics and Controlled Fusion, vol. 39, no. 12B, pp. B115-B127, Dec. 1997,
            doi: https://doi.org/10.1088/0741-3335/39/12b/009.
        """
        return (
            0.031e0
            * pcur**0.95e0
            * b_plasma_toroidal_on_axis**0.25e0
            * p_plasma_loss_mw ** (-0.67e0)
            * dnla19**0.35e0
            * rmajor**1.92e0
            * aspect ** (-0.08e0)
            * kappa**0.63e0
            * afuel**0.42e0
        )

    @staticmethod
    def iter_h97p_elmy_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        dnla19: float,
        rmajor: float,
        aspect: float,
        kappa: float,
        afuel: float,
    ) -> float:
        """Calculate the ELMy ITER H-mode scaling - ITER H97-P(y) confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        rmajor :
            Plasma major radius [m]
        aspect :
            Aspect ratio
        kappa :
            Plasma elongation
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: ITER H97-P(y) confinement time [s]

        References
        ----------
            - I. C. Database and M. W. G. (presented Cordey), “Energy confinement scaling and the extrapolation to ITER,”
            Plasma Physics and Controlled Fusion, vol. 39, no. 12B, pp. B115-B127, Dec. 1997,
            doi: https://doi.org/10.1088/0741-3335/39/12b/009.

            - International Atomic Energy Agency, Vienna (Austria), "Technical basis for the ITER final design report, cost review and safety analysis (FDR)",
            no.16. Dec. 1998.
        """
        return (
            0.029e0
            * pcur**0.90e0
            * b_plasma_toroidal_on_axis**0.20e0
            * p_plasma_loss_mw ** (-0.66e0)
            * dnla19**0.40e0
            * rmajor**2.03e0
            * aspect ** (-0.19e0)
            * kappa**0.92e0
            * afuel**0.2e0
        )

    @staticmethod
    def iter_96p_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        kappa95: float,
        rmajor: float,
        aspect: float,
        dnla19: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the ITER-96P (= ITER-97L) L-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        kappa95 :
            Plasma elongation at 95% flux surface
        rmajor :
            Plasma major radius [m]
        aspect :
            Aspect ratio
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: ITER-96P confinement time [s]

        Notes
        -----
            - The thermal energy confinement time is given below

        References
        ----------
            - S. B. Kaye et al., “ITER L mode confinement database,”
            Nuclear Fusion, vol. 37, no. 9, pp. 1303-1328, Sep. 1997,
            doi: https://doi.org/10.1088/0029-5515/37/9/i10.
        """
        return (
            0.023e0
            * pcur**0.96e0
            * b_plasma_toroidal_on_axis**0.03e0
            * kappa95**0.64e0
            * rmajor**1.83e0
            * aspect**0.06e0
            * dnla19**0.40e0
            * afuel**0.20e0
            * p_plasma_loss_mw ** (-0.73e0)
        )

    @staticmethod
    def valovic_elmy_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        afuel: float,
        rmajor: float,
        rminor: float,
        kappa: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Valovic modified ELMy-H mode scaling confinement time

        Parameters
        ----------
        hfact :
            H-factor
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        afuel :
            Fuel atomic mass number
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        kappa :
            Plasma elongation
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Valovic modified ELMy-H mode confinement time [s]
        """
        return (
            0.067e0
            * pcur**0.9e0
            * b_plasma_toroidal_on_axis**0.17e0
            * dnla19**0.45e0
            * afuel**0.05e0
            * rmajor**1.316e0
            * rminor**0.79e0
            * kappa**0.56e0
            * p_plasma_loss_mw ** (-0.68e0)
        )

    @staticmethod
    def kaye_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        rmajor: float,
        aspect: float,
        dnla19: float,
        afuel: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Kaye PPPL Workshop April 1998 L-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        kappa :
            Plasma elongation
        rmajor :
            Plasma major radius [m]
        aspect :
            Aspect ratio
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        afuel :
            Fuel atomic mass number
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Kaye PPPL Workshop confinement time [s]

        References
        ----------
            - Kaye PPPL Workshop April 1998
        """
        return (
            0.021e0
            * pcur**0.81e0
            * b_plasma_toroidal_on_axis**0.14e0
            * kappa**0.7e0
            * rmajor**2.01e0
            * aspect ** (-0.18e0)
            * dnla19**0.47e0
            * afuel**0.25e0
            * p_plasma_loss_mw ** (-0.73e0)
        )

    @staticmethod
    def iter_pb98py_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the ITERH-PB98P(y) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa :
            Plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: ITERH-PB98P(y) ELMy H-mode confinement time [s]
        """
        return (
            0.0615e0
            * pcur**0.9e0
            * b_plasma_toroidal_on_axis**0.1e0
            * dnla19**0.4e0
            * p_plasma_loss_mw ** (-0.66e0)
            * rmajor**2
            * kappa**0.75e0
            * aspect ** (-0.66e0)
            * afuel**0.2e0
        )

    @staticmethod
    def iter_ipb98y_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the IPB98(y) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa :
            IPB sprcific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: IPB98(y) ELMy H-mode confinement time [s]

        Notes
        -----
            - Unlike the other IPB98 scaling laws, the IPB98(y) scaling law uses the true separatrix elongation.
            - See correction paper below for more information

        References
        ----------
            - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
            Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0365e0
            * pcur**0.97e0
            * b_plasma_toroidal_on_axis**0.08e0
            * dnla19**0.41e0
            * p_plasma_loss_mw ** (-0.63e0)
            * rmajor**1.93e0
            * kappa**0.67e0
            * aspect ** (-0.23e0)
            * afuel**0.2e0
        )

    @staticmethod
    def iter_ipb98y1_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the IPB98(y,1) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: IPB98(y,1) ELMy H-mode confinement time [s]

        Notes
        -----
            - See correction paper below for more information about the re-definition of the elongation used.

        References
        ----------
            - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
            Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0503e0
            * pcur**0.91e0
            * b_plasma_toroidal_on_axis**0.15e0
            * dnla19**0.44e0
            * p_plasma_loss_mw ** (-0.65e0)
            * rmajor**2.05e0
            * kappa_ipb**0.72e0
            * aspect ** (-0.57e0)
            * afuel**0.13e0
        )

    @staticmethod
    def iter_ipb98y2_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the IPB98(y,2) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: IPB98(y,2) ELMy H-mode confinement time [s]

        Notes
        -----
            - See correction paper below for more information about the re-definition of the elongation used.

        References
        ----------
            - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
            Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0562e0
            * pcur**0.93e0
            * b_plasma_toroidal_on_axis**0.15e0
            * dnla19**0.41e0
            * p_plasma_loss_mw ** (-0.69e0)
            * rmajor**1.97e0
            * kappa_ipb**0.78e0
            * aspect ** (-0.58e0)
            * afuel**0.19e0
        )

    @staticmethod
    def iter_ipb98y3_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the IPB98(y,3) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: IPB98(y,3) ELMy H-mode confinement time [s]

        Notes
        -----
            - See correction paper below for more information about the re-definition of the elongation used.

        References
        ----------
            - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
            Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0564e0
            * pcur**0.88e0
            * b_plasma_toroidal_on_axis**0.07e0
            * dnla19**0.40e0
            * p_plasma_loss_mw ** (-0.69e0)
            * rmajor**2.15e0
            * kappa_ipb**0.78e0
            * aspect ** (-0.64e0)
            * afuel**0.20e0
        )

    @staticmethod
    def iter_ipb98y4_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the IPB98(y,4) ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: IPB98(y,4) ELMy H-mode confinement time [s]

        Notes
        -----
            - See correction paper below for more information about the re-definition of the elongation used.

        References
        ----------
            - I. P. E. G. on C. Transport, I. P. E. G. on C. Database, and I. P. B. Editors, “Chapter 2: Plasma confinement and transport,”
            Nuclear Fusion, vol. 39, no. 12, pp. 2175-2249, Dec. 1999, doi: https://doi.org/10.1088/0029-5515/39/12/302.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0587e0
            * pcur**0.85e0
            * b_plasma_toroidal_on_axis**0.29e0
            * dnla19**0.39e0
            * p_plasma_loss_mw ** (-0.70e0)
            * rmajor**2.08e0
            * kappa_ipb**0.76e0
            * aspect ** (-0.69e0)
            * afuel**0.17e0
        )

    @staticmethod
    def iss95_stellarator_confinement_time(
        rminor: float,
        rmajor: float,
        dnla19: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        iotabar: float,
    ) -> float:
        """Calculate the ISS95 stellarator scaling confinement time

        Parameters
        ----------
        rminor :
            Plasma minor radius [m]
        rmajor :
            Plasma major radius [m]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        iotabar :
            Rotational transform

        Returns
        -------
        :
            float: ISS95 stellarator confinement time [s]

        References
        ----------
            - U. Stroth et al., “Energy confinement scaling from the international stellarator database,”
            vol. 36, no. 8, pp. 1063-1077, Aug. 1996, doi: https://doi.org/10.1088/0029-5515/36/8/i11.
        """
        return (
            0.079e0
            * rminor**2.21e0
            * rmajor**0.65e0
            * dnla19**0.51e0
            * b_plasma_toroidal_on_axis**0.83e0
            * p_plasma_loss_mw ** (-0.59e0)
            * iotabar**0.4e0
        )

    @staticmethod
    def iss04_stellarator_confinement_time(
        rminor: float,
        rmajor: float,
        dnla19: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        iotabar: float,
    ) -> float:
        """Calculate the ISS04 stellarator scaling confinement time

        Parameters
        ----------
        rminor :
            Plasma minor radius [m]
        rmajor :
            Plasma major radius [m]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        iotabar :
            Rotational transform

        Returns
        -------
        :
            float: ISS04 stellarator confinement time [s]

        References
        ----------
            - H. Yamada et al., “Characterization of energy confinement in net-current free plasmas using the extended International Stellarator Database,”
            vol. 45, no. 12, pp. 1684-1693, Nov. 2005, doi: https://doi.org/10.1088/0029-5515/45/12/024.
        """
        return (
            0.134e0
            * rminor**2.28e0
            * rmajor**0.64e0
            * dnla19**0.54e0
            * b_plasma_toroidal_on_axis**0.84e0
            * p_plasma_loss_mw ** (-0.61e0)
            * iotabar**0.41e0
        )

    @staticmethod
    def ds03_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa95: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the DS03 beta-independent H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa95 :
            Plasma elongation at 95% flux surface
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: DS03 beta-independent H-mode confinement time [s]

        References
        ----------
            - T. C. Luce, C. C. Petty, and J. G. Cordey, “Application of dimensionless parameter scaling techniques to the design and interpretation of magnetic fusion experiments,”
            Plasma Physics and Controlled Fusion, vol. 50, no. 4, p. 043001, Mar. 2008,
            doi: https://doi.org/10.1088/0741-3335/50/4/043001.
        """
        return (
            0.028e0
            * pcur**0.83e0
            * b_plasma_toroidal_on_axis**0.07e0
            * dnla19**0.49e0
            * p_plasma_loss_mw ** (-0.55e0)
            * rmajor**2.11e0
            * kappa95**0.75e0
            * aspect ** (-0.3e0)
            * afuel**0.14e0
        )

    @staticmethod
    def murari_confinement_time(
        pcur: float,
        rmajor: float,
        kappa_ipb: float,
        dnla19: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Murari H-mode energy confinement scaling time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Murari confinement time [s]

        Notes
        -----
            - This scaling uses the IPB defintiion of elongation, see reference for more information.

        References
        ----------
            - A. Murari, E. Peluso, Michela Gelfusa, I. Lupelli, and P. Gaudio, “A new approach to the formulation and validation of scaling expressions for plasma confinement in tokamaks,”
            Nuclear Fusion, vol. 55, no. 7, pp. 073009-073009, Jun. 2015, doi: https://doi.org/10.1088/0029-5515/55/7/073009.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.0367
            * pcur**1.006
            * rmajor**1.731
            * kappa_ipb**1.450
            * p_plasma_loss_mw ** (-0.735)
            * (
                dnla19**0.448
                / (1.0 + np.exp(-9.403 * (dnla19 / b_plasma_toroidal_on_axis) ** -1.365))
            )
        )

    @staticmethod
    def petty08_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
    ) -> float:
        """Calculate the beta independent dimensionless Petty08 confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio

        Returns
        -------
        :
            float: Petty08 confinement time [s]

        Notes
        -----
            - This scaling uses the IPB defintiion of elongation, see reference for more information.

        References
        ----------
            - C. C. Petty, “Sizing up plasmas using dimensionless parameters,”
            Physics of Plasmas, vol. 15, no. 8, Aug. 2008, doi: https://doi.org/10.1063/1.2961043.

            - None Otto Kardaun, N. K. Thomsen, and None Alexander Chudnovskiy, “Corrections to a sequence of papers in Nuclear Fusion,”
            Nuclear Fusion, vol. 48, no. 9, pp. 099801-099801, Aug. 2008, doi: https://doi.org/10.1088/0029-5515/48/9/099801.
        """
        return (
            0.052e0
            * pcur**0.75e0
            * b_plasma_toroidal_on_axis**0.3e0
            * dnla19**0.32e0
            * p_plasma_loss_mw ** (-0.47e0)
            * rmajor**2.09e0
            * kappa_ipb**0.88e0
            * aspect ** (-0.84e0)
        )

    @staticmethod
    def lang_high_density_confinement_time(
        plasma_current: float,
        b_plasma_toroidal_on_axis: float,
        nd_plasma_electron_line: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        rminor: float,
        q: float,
        qstar: float,
        aspect: float,
        afuel: float,
        kappa_ipb: float,
    ) -> float:
        """Calculate the high density relevant confinement time

        Parameters
        ----------
        plasma_current :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        nd_plasma_electron_line :
            Line averaged electron density [m**-3]
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        rminor :
            Plasma minor radius [m]
        q :
            Safety factor
        qstar :
            Equivalent cylindrical edge safety factor
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number
        kappa_ipb :
            Plasma elongation at 95% flux surface

        Returns
        -------
        :
            float: High density relevant confinement time [s]

        References
        ----------
            - P. T. Lang, C. Angioni, R. M. M. Dermott, R. Fischer, and H. Zohm, “Pellet Induced High Density Phases during ELM Suppression in ASDEX Upgrade,”
            24th IAEA Conference Fusion Energy, 2012, Oct. 2012,
            Available: https://www.researchgate.net/publication/274456104_Pellet_Induced_High_Density_Phases_during_ELM_Suppression_in_ASDEX_Upgrade
        """
        qratio = q / qstar
        n_gw = 1.0e14 * plasma_current / (np.pi * rminor * rminor)
        nratio = nd_plasma_electron_line / n_gw
        return (
            6.94e-7
            * plasma_current**1.3678e0
            * b_plasma_toroidal_on_axis**0.12e0
            * nd_plasma_electron_line**0.032236e0
            * (p_plasma_loss_mw * 1.0e6) ** (-0.74e0)
            * rmajor**1.2345e0
            * kappa_ipb**0.37e0
            * aspect**2.48205e0
            * afuel**0.2e0
            * qratio**0.77e0
            * aspect ** (-0.9e0 * np.log(aspect))
            * nratio ** (-0.22e0 * np.log(nratio))
        )

    @staticmethod
    def hubbard_nominal_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla20: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Hubbard 2017 I-mode confinement time scaling - nominal

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Hubbard confinement time [s]

        References
        ----------
            - A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,”
            Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.
        """
        return (
            0.014e0
            * pcur**0.68e0
            * b_plasma_toroidal_on_axis**0.77e0
            * dnla20**0.02e0
            * p_plasma_loss_mw ** (-0.29e0)
        )

    @staticmethod
    def hubbard_lower_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla20: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Hubbard 2017 I-mode confinement time scaling - lower

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Hubbard confinement time [s]

        References
        ----------
            - A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,”
            Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.
        """
        return (
            0.014e0
            * pcur**0.60e0
            * b_plasma_toroidal_on_axis**0.70e0
            * dnla20 ** (-0.03e0)
            * p_plasma_loss_mw ** (-0.33e0)
        )

    @staticmethod
    def hubbard_upper_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla20: float,
        p_plasma_loss_mw: float,
    ) -> float:
        """Calculate the Hubbard 2017 I-mode confinement time scaling - upper

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]

        Returns
        -------
        :
            float: Hubbard confinement time [s]

        References
        ----------
            - A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,”
            Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.
        """
        return (
            0.014e0
            * pcur**0.76e0
            * b_plasma_toroidal_on_axis**0.84e0
            * dnla20**0.07
            * p_plasma_loss_mw ** (-0.25e0)
        )

    @staticmethod
    def menard_nstx_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the Menard NSTX ELMy H-mode scaling confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: Menard NSTX ELMy H-mode confinement time [s]

        Notes
        -----
            - "The leading NSTX conﬁnement scaling coefﬁcient is chosen such that the ITER and ST energy conﬁnement times are
            identical for a reference NSTX scenario"
            - Assumes IPB98(y,2) exponents are applicable where the ST exponents are not yet determined, i.e.
            the species mass, major radius, inverse aspect ratio and elongation. Hence here we use the IPB98(y,2) definition
            of elongation.

        References
        ----------
            - J. E. Menard, “Compact steady-state tokamak performance dependence on magnet and core physics limits,”
            Philosophical Transactions of the Royal Society A, vol. 377, no. 2141, pp. 20170440-20170440, Feb. 2019,
            doi: https://doi.org/10.1098/rsta.2017.0440.
        """
        return (
            0.095e0
            * pcur**0.57e0
            * b_plasma_toroidal_on_axis**1.08e0
            * dnla19**0.44e0
            * p_plasma_loss_mw ** (-0.73e0)
            * rmajor**1.97e0
            * kappa_ipb**0.78e0
            * aspect ** (-0.58e0)
            * afuel**0.19e0
        )

    @classmethod
    def menard_nstx_petty08_hybrid_confinement_time(
        cls,
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        kappa_ipb: float,
        aspect: float,
        afuel: float,
    ) -> float:
        """Calculate the Menard NSTX-Petty hybrid confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Line averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        kappa_ipb :
            IPB specific plasma separatrix elongation
        aspect :
            Aspect ratio
        afuel :
            Fuel atomic mass number

        Returns
        -------
        :
            float: Menard NSTX-Petty hybrid confinement time [s]

        Notes
        -----
            - Assuming a linear interpolation in (1/aspect) between the two scalings

        References
        ----------
            - J. E. Menard, “Compact steady-state tokamak performance dependence on magnet and core physics limits,”
            Philosophical Transactions of the Royal Society A, vol. 377, no. 2141, pp. 20170440-20170440, Feb. 2019,
            doi: https://doi.org/10.1098/rsta.2017.0440.
        """
        # Equivalent to A > 2.5, use Petty scaling
        if (1.0e0 / aspect) <= 0.4e0:
            return cls.petty08_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa_ipb,
                aspect,
            )

        #  Equivalent to A < 1.7, use NSTX scaling
        if (1.0e0 / aspect) >= 0.6e0:
            return cls.menard_nstx_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa_ipb,
                aspect,
                afuel,
            )
        return (((1.0e0 / aspect) - 0.4e0) / (0.6e0 - 0.4e0)) * (
            cls.menard_nstx_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa_ipb,
                aspect,
                afuel,
            )
        ) + ((0.6e0 - (1.0e0 / aspect)) / (0.6e0 - 0.4e0)) * (
            cls.petty08_confinement_time(
                pcur,
                b_plasma_toroidal_on_axis,
                dnla19,
                p_plasma_loss_mw,
                rmajor,
                kappa_ipb,
                aspect,
            )
        )

    @staticmethod
    def nstx_gyro_bohm_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        dnla20: float,
    ) -> float:
        """Calculate the NSTX gyro-Bohm confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Net Heating power [MW]
        rmajor :
            Plasma major radius [m]
        dnla20 :
            Line averaged electron density in units of 10**20 m**-3

        Returns
        -------
        :
            float: NSTX gyro-Bohm confinement time [s]

        References
        ----------
            - P. F. Buxton, L. Connor, A. E. Costley, Mikhail Gryaznevich, and S. McNamara,
            “On the energy confinement time in spherical tokamaks: implications for the design of pilot plants and fusion reactors,”
            vol. 61, no. 3, pp. 035006-035006, Jan. 2019, doi: https://doi.org/10.1088/1361-6587/aaf7e5.
        """
        return (
            0.21e0
            * pcur**0.54e0
            * b_plasma_toroidal_on_axis**0.91e0
            * p_plasma_loss_mw ** (-0.38e0)
            * rmajor**2.14e0
            * dnla20 ** (-0.05e0)
        )

    @staticmethod
    def itpa20_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        dnla19: float,
        p_plasma_loss_mw: float,
        rmajor: float,
        triang: float,
        kappa_ipb: float,
        eps: float,
        aion: float,
    ) -> float:
        """Calculate the ITPA20 Issue #3164 confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        dnla19 :
            Central line-averaged electron density in units of 10**19 m**-3
        p_plasma_loss_mw :
            Thermal power lost due to transport through the LCFS [MW]
        rmajor :
            Plasma major radius [m]
        triang :
            Triangularity
        kappa_ipb :
            IPB specific plasma separatrix elongation
        eps :
            Inverse aspect ratio
        aion :
            Average mass of all ions (amu)

        Returns
        -------
        :
            float: ITPA20 confinement time [s]

        Notes
        -----
            - Mass term is the effective mass of the plasma, so we assume the total ion mass here
            - This scaling uses the IPB defintiion of elongation, see reference for more information.

        References
        ----------
            - G. Verdoolaege et al., “The updated ITPA global H-mode confinement database: description and analysis,”
            Nuclear Fusion, vol. 61, no. 7, pp. 076006-076006, Jan. 2021, doi: https://doi.org/10.1088/1741-4326/abdb91.
        """
        return (
            0.0534
            * pcur**0.976
            * b_plasma_toroidal_on_axis**0.218
            * dnla19**0.2442
            * p_plasma_loss_mw ** (-0.6687)
            * rmajor**1.710
            * (1 + triang) ** 0.362
            * kappa_ipb**0.799
            * eps**0.354
            * aion**0.195
        )

    @staticmethod
    def itpa20_il_confinement_time(
        pcur: float,
        b_plasma_toroidal_on_axis: float,
        p_plasma_loss_mw: float,
        dnla19: float,
        aion: float,
        rmajor: float,
        triang: float,
        kappa_ipb: float,
    ) -> float:
        """Calculate the ITPA20-IL Issue #1852 confinement time

        Parameters
        ----------
        pcur :
            Plasma current [MA]
        b_plasma_toroidal_on_axis :
            Toroidal magnetic field [T]
        p_plasma_loss_mw :
            Thermal power lost due to transport through the LCFS [MW]
        dnla19 :
            Central line-averaged electron density in units of 10**19 m**-3
        aion :
            Average mass of all ions (amu)
        rmajor :
            Plasma major radius [m]
        triang :
            Triangularity
        kappa_ipb :
            IPB specific plasma separatrix elongation

        Returns
        -------
        :
            float: ITPA20-IL confinement time [s]

        Notes
        -----
            - Mass term is the effective mass of the plasma, so we assume the total ion mass here
            - This scaling uses the IPB defintiion of elongation, see reference for more information.

        References
        ----------
            - T. Luda et al., “Validation of a full-plasma integrated modeling approach on ASDEX Upgrade,”
            Nuclear Fusion, vol. 61, no. 12, pp. 126048-126048, Nov. 2021, doi: https://doi.org/10.1088/1741-4326/ac3293.
        """
        return (
            0.0670
            * pcur**1.291
            * b_plasma_toroidal_on_axis**-0.134
            * dnla19**0.1473
            * p_plasma_loss_mw ** (-0.6442)
            * rmajor**1.194
            * (1 + triang) ** 0.560
            * kappa_ipb**0.673
            * aion**0.302
        )
