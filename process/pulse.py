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

        self.burn(output=output)

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
                pfcoil_variables.powohres
                / (1.0e6 * pfcoil_variables.ric[pfcoil_variables.n_cs_pf_coils - 1])
                ** 2
            )

        #  Central Solenoid bus resistance (ohms) (assumed to include power supply)
        #  Bus parameters taken from routine PFPWR.

        pfbusl = 8.0e0 * physics_variables.rmajor + 140.0e0
        albusa = (
            abs(pfcoil_variables.cptdin[pfcoil_variables.n_cs_pf_coils - 1]) / 100.0e0
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

    def burn(self, output: bool):
        """Routine to calculate the burn time for a pulsed reactor
        author: C A Gardner, AEA Fusion, Culham Laboratory
        author: P J Knight, CCFE, Culham Science Centre
        author: R Kemp, CCFE, Culham Science Centre

        This routine calculates the burn time for a pulsed reactor.
        Work File Note F/MPE/MOD/CAG/PROCESS/PULSE/0012

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        if pulse_variables.i_pulsed_plant != 1:
            return

        #  Volt-seconds required to produce plasma current during start-up
        #  (i.e. up to start of flat top)

        vssoft = (
            physics_variables.vs_plasma_res_ramp + physics_variables.vs_plasma_ind_ramp
        )

        #  Total volt-seconds available during flat-top (heat + burn)
        #  (Previously calculated as (abs(pfcoil_variables.vstot) - vssoft) )

        vsmax = (
            -pfcoil_variables.vsbn
        )  # pfcoil_variables.vsbn is (or should be...) negative

        #  Loop voltage during flat-top (including MHD sawtooth enhancement)

        v_plasma_loop_burn = (
            physics_variables.plasma_current
            * physics_variables.res_plasma
            * physics_variables.inductive_current_fraction
            * physics_variables.csawth
        )

        #  Burn time (s)

        tb = vsmax / v_plasma_loop_burn - times_variables.t_fusion_ramp
        if tb < 0.0e0:
            error_handling.fdiags[0] = tb
            error_handling.fdiags[1] = vsmax
            error_handling.fdiags[2] = v_plasma_loop_burn
            error_handling.fdiags[3] = times_variables.t_fusion_ramp
            error_handling.report_error(93)

        times_variables.t_burn = max(0.0e0, tb)

        #  Output section

        if output:
            po.osubhd(self.outfile, "Volt-second considerations:")

            po.ovarre(
                self.outfile,
                "Total V-s capability of Central Solenoid/PF coils (Wb)",
                "(abs(vstot))",
                abs(pfcoil_variables.vstot),
            )
            po.ovarre(
                self.outfile,
                "Required volt-seconds during start-up (Wb)",
                "(vssoft)",
                vssoft,
            )
            po.ovarre(
                self.outfile,
                "Available volt-seconds during burn (Wb)",
                "(vsmax)",
                vsmax,
            )

            if tb <= 0.0e0:
                po.ocmmnt(
                    self.outfile,
                    "   Error... burn time is zero; insufficient volt-seconds#",
                )
