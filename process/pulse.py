from process import process_output as po
from process.fortran import (
    constants,
    constraint_variables,
    error_handling,
    numerics,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    pulse_variables,
    times_variables,
)


class Pulse:
    def __init__(self):
        self.outfile = constants.nout

    def run(self, output: bool) -> None:
        """Caller for the pulsed reactor model
        author: C A Gardner, AEA Fusion, Culham Laboratory
        author: P J Knight, CCFE, Culham Science Centre

        This calls the routines relevant to a pulsed reactor scenario.
        Work File Notes F/MPE/MOD/CAG/PROCESS/PULSE

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        self.tohswg(output=output)

        #  Burn time calculation

        times_variables.t_burn = self.calculate_burn_time(
            i_pulsed_plant=pulse_variables.i_pulsed_plant,
            vs_cs_pf_total_burn=pfcoil_variables.vs_cs_pf_total_burn,
            v_plasma_loop_burn=physics_variables.v_plasma_loop_burn,
            t_fusion_ramp=times_variables.t_fusion_ramp,
        )

    def tohswg(self, output: bool) -> None:
        """Routine to calculate the plasma current ramp-up time
        author: C A Gardner, AEA Fusion, Culham Laboratory
        author: P J Knight, CCFE, Culham Science Centre

        This routine calculates the plasma current ramp-up time
        for a pulsed reactor.
        Work File Note F/MPE/MOD/CAG/PROCESS/PULSE/0013
        Work File Note F/PL/PJK/PROCESS/CODE/050

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
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
        i_pulsed_plant: int,
        vs_cs_pf_total_burn: float,
        v_plasma_loop_burn: float,
        t_fusion_ramp: float,
    ) -> float:
        """
        Calculate the burn time for a pulsed reactor.

        This routine computes the burn time for a pulsed reactor scenario,
        based on the total Vs available in the CS and PF coils and the
        plasma loop voltage during burn. It also checks for negative burn time
        and reports an error if encountered.

        :param i_pulsed_plant: Indicator for pulsed plant (1 for pulsed, 0 otherwise)
        :type i_pulsed_plant: int
        :param vs_cs_pf_total_burn: Total volt-seconds in CS and PF coils available for burn (V·s)
        :type vs_cs_pf_total_burn: float
        :param v_plasma_loop_burn: Plasma loop voltage during burn (V)
        :type v_plasma_loop_burn: float
        :param t_fusion_ramp: Time for fusion ramp (s)
        :type t_fusion_ramp: float
        :return: Calculated burn time (s)
        :rtype: float

        :raises: Reports error 93 if calculated burn time is negative.

        """
        if i_pulsed_plant != 1:
            return None

        t_burn = (abs(vs_cs_pf_total_burn) / v_plasma_loop_burn) - t_fusion_ramp

        if t_burn < 0.0e0:
            error_handling.fdiags[0] = t_burn
            error_handling.fdiags[1] = vs_cs_pf_total_burn
            error_handling.fdiags[2] = v_plasma_loop_burn
            error_handling.fdiags[3] = t_fusion_ramp
            error_handling.report_error(93)

        return t_burn


def init_pulse_variables():
    """Initialise the pulse variables"""
    pulse_variables.bctmp = 320.0
    pulse_variables.dtstor = 300.0
    pulse_variables.istore = 1
    pulse_variables.itcycl = 1
    pulse_variables.i_pulsed_plant = 0
