from process.fortran import physics_variables
from process.fortran import times_variables
from process.fortran import pfcoil_variables
from process.fortran import constraint_variables
from process.fortran import constants
from process.fortran import pulse_variables
from process.fortran import numerics
from process.fortran import pf_power_variables
from process.fortran import process_output as po
from process.fortran import error_handling


class Pulse:
    def __init__(self):
        self.outfile = constants.nout

    def run(self, output: bool) -> None:
        """Caller for the pulsed reactor model
        author: C A Gardner, AEA Fusion, Culham Laboratory
        author: P J Knight, CCFE, Culham Science Centre

        This calls the routines relevant to a pulsed reactor scenario.
        Work File Notes F/MPE/MOD/CAG/PROCESS/PULSE
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        if pulse_variables.lpulse != 1:
            return

        #  Current/turn in Central Solenoid at beginning of pulse (A/turn)

        ioht1 = pfcoil_variables.cpt[pfcoil_variables.nohc - 1, 1]

        #  Current/turn in Central Solenoid at start of flat-top (A/turn)

        ioht2 = pfcoil_variables.cpt[pfcoil_variables.nohc - 1, 2]

        #  Central Solenoid resistance (ohms)

        if pfcoil_variables.ipfres == 0:
            r = 0.0e0
        else:
            r = (
                pfcoil_variables.powohres
                / (1.0e6 * pfcoil_variables.ric[pfcoil_variables.nohc - 1]) ** 2
            )

        #  Central Solenoid bus resistance (ohms) (assumed to include power supply)
        #  Bus parameters taken from routine PFPWR.

        pfbusl = 8.0e0 * physics_variables.rmajor + 140.0e0
        albusa = abs(pfcoil_variables.cptdin[pfcoil_variables.nohc - 1]) / 100.0e0

        rho = 1.5e0 * 2.62e-4 * pfbusl / albusa

        #  Central Solenoid power source emf (volts)

        v = pf_power_variables.vpfskv * 1.0e3

        #  Mutual inductance between Central Solenoid and plasma (H)

        m = pfcoil_variables.sxlg[pfcoil_variables.nohc - 1, pfcoil_variables.ncirt - 1]

        #  Self inductance of Central Solenoid (H)

        loh = pfcoil_variables.sxlg[
            pfcoil_variables.nohc - 1, pfcoil_variables.nohc - 1
        ]

        #  Maximum rate of change of plasma current (A/s)
        #  - now a function of the plasma current itself (previously just 0.5e6)

        ipdot = 0.0455e0 * physics_variables.plascur

        #  Minimum plasma current ramp-up time (s)
        #  - corrected (bus resistance is not a function of pfcoil_variables.turns)

        constraint_variables.tohsmn = (
            loh
            * (ioht2 - ioht1)
            / (
                ioht2 * (r * pfcoil_variables.turns[pfcoil_variables.nohc - 1] + rho)
                - v
                + m * ipdot
            )
        )

        #  Output section

        if output == 1:
            if numerics.active_constraints[40]:
                po.osubhd(self.outfile, "Central solenoid considerations:")
                po.ovarre(
                    self.outfile,
                    "Minimum plasma current ramp-up time (s)",
                    "(tohsmn)",
                    constraint_variables.tohsmn,
                )

    def burn(self, output: bool):
        """Routine to calculate the burn time for a pulsed reactor
        author: C A Gardner, AEA Fusion, Culham Laboratory
        author: P J Knight, CCFE, Culham Science Centre
        author: R Kemp, CCFE, Culham Science Centre

        This routine calculates the burn time for a pulsed reactor.
        Work File Note F/MPE/MOD/CAG/PROCESS/PULSE/0012
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """
        if pulse_variables.lpulse != 1:
            return

        #  Volt-seconds required to produce plasma current during start-up
        #  (i.e. up to start of flat top)

        vssoft = physics_variables.vsres + physics_variables.vsind

        #  Total volt-seconds available during flat-top (heat + burn)
        #  (Previously calculated as (abs(pfcoil_variables.vstot) - vssoft) )

        vsmax = (
            -pfcoil_variables.vsbn
        )  # pfcoil_variables.vsbn is (or should be...) negative

        #  Loop voltage during flat-top (including MHD sawtooth enhancement)

        vburn = (
            physics_variables.plascur
            * physics_variables.rplas
            * physics_variables.facoh
            * physics_variables.csawth
        )

        #  Burn time (s)

        tb = vsmax / vburn - times_variables.t_fusion_ramp
        if tb < 0.0e0:
            error_handling.fdiags[0] = tb
            error_handling.fdiags[1] = vsmax
            error_handling.fdiags[2] = vburn
            error_handling.fdiags[3] = times_variables.t_fusion_ramp
            error_handling.report_error(93)

        times_variables.tburn = max(0.0e0, tb)

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
