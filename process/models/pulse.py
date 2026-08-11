"""Module containing the Pulse class for pulsed reactor calculations."""

import logging
from dataclasses import dataclass, fields
from typing import ClassVar

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure.pfcoil_variables import PFConductorModel

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class PulseTimings:
    """Dataclass to hold the timing parameters for a pulsed reactor."""

    t_plant_pulse_coil_precharge: float
    """Time for coil precharge (s)"""
    t_plant_pulse_plasma_current_ramp_up: float
    """Time for plasma current ramp-up (s)"""
    t_plant_pulse_fusion_ramp: float
    """Time for fusion ramp (s)"""
    t_plant_pulse_burn: float
    """Time for burn (s)"""
    t_plant_pulse_plasma_current_ramp_down: float
    """Time for plasma current ramp-down (s)"""
    t_plant_pulse_dwell: float
    """Time for dwell (s)"""

    point_abbreviations: ClassVar[tuple[str, ...]] = (
        "BOP",
        "EOR",
        "BOF",
        "EOF",
        "EOP",
        "Dwell",
    )

    point_labels: ClassVar[tuple[str, ...]] = (
        "Coil precharge",
        "$I_{\\text{p}}$ Ramp-Up",
        "Fusion ramp",
        "Burn",
        "$I_{\\text{p}}$ ramp-down",
        "Dwell",
    )

    def __post_init__(self) -> None:
        """Validate class metadata against the timing fields.

        Raises
        ------
        ValueError
            If the number of point labels or abbreviations does not match the number
            of timing fields.
        """
        n_timing_fields = len(fields(self))
        if len(self.point_labels) != n_timing_fields:
            raise ValueError(
                "PulseTimings.point_labels must contain exactly "
                f"{n_timing_fields} entries; got {len(self.point_labels)}."
            )
        if len(self.point_abbreviations) != n_timing_fields:
            raise ValueError(
                "PulseTimings.point_abbreviations must contain exactly "
                f"{n_timing_fields} entries; got {len(self.point_abbreviations)}."
            )

    @property
    def plasma_present(self) -> float:
        """Calculate the total time during which plasma is present in the reactor."""
        return (
            self.t_plant_pulse_plasma_current_ramp_up
            + self.t_plant_pulse_fusion_ramp
            + self.t_plant_pulse_burn
            + self.t_plant_pulse_plasma_current_ramp_down
        )

    @property
    def no_burn(self) -> float:
        """Calculate the total time excluding the burn phase."""
        return (
            self.t_plant_pulse_coil_precharge
            + self.t_plant_pulse_plasma_current_ramp_up
            + self.t_plant_pulse_plasma_current_ramp_down
            + self.t_plant_pulse_dwell
            + self.t_plant_pulse_fusion_ramp
        )

    @property
    def total(self) -> float:
        """Calculate the total time including the burn phase."""
        return self.no_burn + self.t_plant_pulse_burn

    @property
    def total_pulse_cumulative(
        self,
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate the cumulative timing points for all pulse phases."""
        t0 = 0.0
        t1 = t0 + self.t_plant_pulse_coil_precharge
        t2 = t1 + self.t_plant_pulse_plasma_current_ramp_up
        t3 = t2 + self.t_plant_pulse_fusion_ramp
        t4 = t3 + self.t_plant_pulse_burn
        t5 = t4 + self.t_plant_pulse_plasma_current_ramp_down
        t6 = t5 + self.t_plant_pulse_dwell
        return (t0, t1, t2, t3, t4, t5, t6)

    @property
    def n_pulse_points_total(self) -> int:
        """Calculate the total number of timing points for all pulse phases."""
        return len(self.total_pulse_cumulative)

    @property
    def pf_active_cumulative(self) -> tuple[float, float, float, float, float, float]:
        """Calculate the cumulative timing points for PF coil active phases."""
        return self.total_pulse_cumulative[:-1]  # Exclude the last point (dwell)

    @property
    def n_pf_active_points_total(self) -> int:
        """Calculate the total number of timing points for PF coil active phases."""
        return len(self.pf_active_cumulative)

    @property
    def n_pf_active_points_intervals(self) -> int:
        """Calculate the total number of timing intervals for PF coil active phases."""
        return int(self.n_pf_active_points_total - 1)


class Pulse(Model):
    """Class containing pulsed reactor calculations"""

    def __init__(self):
        self.outfile = constants.NOUT

    def output(self):
        """Write the results to the main output file (OUT.DAT)."""
        self.run(output=True)

    def run(self, output: bool = False):
        """Caller for the pulsed reactor model

        This calls the routines relevant to a pulsed reactor scenario.
        Work File Notes F/MPE/MOD/CAG/PROCESS/PULSE

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """
        if self.data.pulse.i_pulsed_plant == 1:
            self.tohswg(output=output)

            #  Burn time calculation

            self.data.times.t_plant_pulse_burn = self.calculate_burn_time(
                vs_cs_pf_total_burn=self.data.pf_coil.vs_cs_pf_total_burn,
                v_plasma_loop_burn=self.data.physics.v_plasma_loop_burn,
                t_plant_pulse_fusion_ramp=self.data.times.t_plant_pulse_fusion_ramp,
            )

    def tohswg(self, output: bool):
        """Routine to calculate the plasma current ramp-up time

        This routine calculates the plasma current ramp-up time
        for a pulsed reactor.
        Work File Note F/MPE/MOD/CAG/PROCESS/PULSE/0013
        Work File Note F/PL/PJK/PROCESS/CODE/050

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """
        if self.data.pulse.i_pulsed_plant != 1:
            return

        #  Current/turn in Central Solenoid at beginning of pulse (A/turn)

        ioht1 = self.data.pf_coil.c_pf_coil_turn[self.data.pf_coil.n_cs_pf_coils - 1, 1]

        #  Current/turn in Central Solenoid at start of flat-top (A/turn)

        ioht2 = self.data.pf_coil.c_pf_coil_turn[self.data.pf_coil.n_cs_pf_coils - 1, 2]

        #  Central Solenoid resistance (ohms)

        if self.data.pf_coil.i_pf_conductor == PFConductorModel.SUPERCONDUCTING:
            r = 0.0e0
        else:
            r = (
                self.data.pf_coil.p_cs_resistive_flat_top
                / (
                    1.0e6
                    * self.data.pf_coil.c_pf_cs_coils_peak_ma[
                        self.data.pf_coil.n_cs_pf_coils - 1
                    ]
                )
                ** 2
            )

        #  Central Solenoid bus resistance (ohms) (assumed to include power supply)
        #  Bus parameters taken from routine PFPWR.

        pfbusl = 8.0e0 * self.data.physics.rmajor + 140.0e0
        albusa = (
            abs(
                self.data.pf_coil.c_pf_coil_turn_peak_input[
                    self.data.pf_coil.n_cs_pf_coils - 1
                ]
            )
            / 100.0e0
        )

        # rho = 1.5e0 * 2.62e-4 * pfbusl / albusa
        # I have removed the fudge factor of 1.5 but included it in the value of
        #  rhopfbus
        rho = self.data.pf_coil.rhopfbus * pfbusl / (albusa / 10000)

        #  Central Solenoid power source emf (volts)

        v = self.data.pf_power.vpfskv * 1.0e3

        #  Mutual inductance between Central Solenoid and plasma (H)

        m = self.data.pf_coil.ind_pf_cs_plasma_mutual[
            self.data.pf_coil.n_cs_pf_coils - 1,
            self.data.pf_coil.n_pf_cs_plasma_circuits - 1,
        ]

        #  Self inductance of Central Solenoid (H)

        loh = self.data.pf_coil.ind_pf_cs_plasma_mutual[
            self.data.pf_coil.n_cs_pf_coils - 1, self.data.pf_coil.n_cs_pf_coils - 1
        ]

        #  Maximum rate of change of plasma current (A/s)
        #  - now a function of the plasma current itself (previously just 0.5e6)

        ipdot = 0.0455e0 * self.data.physics.plasma_current

        #  Minimum plasma current ramp-up time (s)
        #  - corrected (bus resistance is not a function of self.data.pf_coil.turns)

        self.data.constraints.t_current_ramp_up_min = (
            loh
            * (ioht2 - ioht1)
            / (
                ioht2
                * (
                    r
                    * self.data.pf_coil.n_pf_coil_turns[
                        self.data.pf_coil.n_cs_pf_coils - 1
                    ]
                    + rho
                )
                - v
                + m * ipdot
            )
        )

        #  Output section

        if output == 1 and self.data.numerics.active_constraints[40]:
            po.osubhd(self.outfile, "Central solenoid considerations:")
            po.ovarre(
                self.outfile,
                "Minimum plasma current ramp-up time (s)",
                "(t_current_ramp_up_min)",
                self.data.constraints.t_current_ramp_up_min,
            )

    @staticmethod
    def calculate_burn_time(
        vs_cs_pf_total_burn: float,
        v_plasma_loop_burn: float,
        t_plant_pulse_fusion_ramp: float,
    ) -> float:
        """Calculate the burn time for a pulsed reactor.

        This routine computes the burn time for a pulsed reactor scenario,
        based on the total Vs available in the CS and PF coils and the
        plasma loop voltage during burn. It also checks for negative burn time
        and reports an error if encountered.

        Parameters
        ----------
        vs_cs_pf_total_burn : float
            Total volt-seconds in CS and PF coils available for burn (V·s)
        v_plasma_loop_burn : float
            Plasma loop voltage during burn (V)
        t_plant_pulse_fusion_ramp : float
            Time for fusion ramp (s)

        Returns
        -------
        float
            Calculated burn time (s)
        """
        t_plant_pulse_burn = (
            abs(vs_cs_pf_total_burn) / v_plasma_loop_burn
        ) - t_plant_pulse_fusion_ramp

        if t_plant_pulse_burn < 0.0e0:
            logger.error(
                "Negative burn time available; reduce t_plant_pulse_fusion_ramp or "
                "raise PF coil V-s capability. %s %s %s %s",
                t_plant_pulse_burn,
                vs_cs_pf_total_burn,
                v_plasma_loop_burn,
                t_plant_pulse_fusion_ramp,
            )

        return t_plant_pulse_burn
