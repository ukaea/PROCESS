import logging

from process.core import constants
from process.core import process_output as po
from process.data_structure import (
    constraint_variables,
    numerics,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    pulse_variables,
    times_variables,
)

logger = logging.getLogger(__name__)


class Pulse:
    def __init__(self):
        self.outfile = constants.NOUT

    def run(self, output: bool):
        """Caller for the pulsed reactor model

        This calls the routines relevant to a pulsed reactor scenario.
        Work File Notes F/MPE/MOD/CAG/PROCESS/PULSE

        Parameters
        ----------
        output :
            indicate whether output should be written to the output file, or not
        """
        if pulse_variables.i_pulsed_plant == 1:
            self.tohswg(output=output)

            #  Burn time calculation

            times_variables.t_plant_pulse_burn = self.calculate_burn_time(
                vs_cs_pf_total_burn=pfcoil_variables.vs_cs_pf_total_burn,
                v_plasma_loop_burn=physics_variables.v_plasma_loop_burn,
                t_plant_pulse_fusion_ramp=times_variables.t_plant_pulse_fusion_ramp,
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
        if pulse_variables.i_pulsed_plant != 1:
            return

        #  Current/turn in Central Solenoid at beginning of pulse (A/turn)

        ioht1 = pfcoil_variables.c_pf_coil_turn[pfcoil_variables.n_cs_pf_coils - 1, 1]

        #  Current/turn in Central Solenoid at start of flat-top (A/turn)

        ioht2 = pfcoil_variables.c_pf_coil_turn[pfcoil_variables.n_cs_pf_coils - 1, 2]

        #  Central Solenoid resistance (ohms)

        if pfcoil_variables.i_pf_conductor == 0:
            r = 0.0e0
        else:
            r = (
                pfcoil_variables.p_cs_resistive_flat_top
                / (
                    1.0e6
                    * pfcoil_variables.c_pf_cs_coils_peak_ma[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                )
                ** 2
            )

        #  Central Solenoid bus resistance (ohms) (assumed to include power supply)
        #  Bus parameters taken from routine PFPWR.

        pfbusl = 8.0e0 * physics_variables.rmajor + 140.0e0
        albusa = (
            abs(
                pfcoil_variables.c_pf_coil_turn_peak_input[
                    pfcoil_variables.n_cs_pf_coils - 1
                ]
            )
            / 100.0e0
        )

        # rho = 1.5e0 * 2.62e-4 * pfbusl / albusa
        #  I have removed the fudge factor of 1.5 but included it in the value of rhopfbus
        rho = pfcoil_variables.rhopfbus * pfbusl / (albusa / 10000)

        #  Central Solenoid power source emf (volts)

        v = pf_power_variables.vpfskv * 1.0e3

        #  Mutual inductance between Central Solenoid and plasma (H)

        m = pfcoil_variables.ind_pf_cs_plasma_mutual[
            pfcoil_variables.n_cs_pf_coils - 1,
            pfcoil_variables.n_pf_cs_plasma_circuits - 1,
        ]

        #  Self inductance of Central Solenoid (H)

        loh = pfcoil_variables.ind_pf_cs_plasma_mutual[
            pfcoil_variables.n_cs_pf_coils - 1, pfcoil_variables.n_cs_pf_coils - 1
        ]

        #  Maximum rate of change of plasma current (A/s)
        #  - now a function of the plasma current itself (previously just 0.5e6)

        ipdot = 0.0455e0 * physics_variables.plasma_current

        #  Minimum plasma current ramp-up time (s)
        #  - corrected (bus resistance is not a function of pfcoil_variables.turns)

        constraint_variables.t_current_ramp_up_min = (
            loh
            * (ioht2 - ioht1)
            / (
                ioht2
                * (
                    r
                    * pfcoil_variables.n_pf_coil_turns[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                    + rho
                )
                - v
                + m * ipdot
            )
        )

        #  Output section

        if output == 1 and numerics.active_constraints[40]:
            po.osubhd(self.outfile, "Central solenoid considerations:")
            po.ovarre(
                self.outfile,
                "Minimum plasma current ramp-up time (s)",
                "(t_current_ramp_up_min)",
                constraint_variables.t_current_ramp_up_min,
            )

    def calculate_burn_time(
        self,
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
                "Negative burn time available; reduce t_plant_pulse_fusion_ramp or raise PF coil V-s capabilit. "
                f"{t_plant_pulse_burn=} {vs_cs_pf_total_burn=} {v_plasma_loop_burn=} {t_plant_pulse_fusion_ramp=}"
            )

        return t_plant_pulse_burn
